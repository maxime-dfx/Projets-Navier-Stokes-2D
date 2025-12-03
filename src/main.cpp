#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <chrono>   
#include <iomanip>  
#include <filesystem> 
#include <vector>
#include <cstdlib> // Pour srand, rand
#include <ctime>   // Pour time

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "DataFile.h"
#include "Function.h"
#include "MACgrid.h"
#include "Laplacian.h"
#include "TimeScheme.h"
#include "ParticleSystem.h" // [NOUVEAU]

using namespace std;
using namespace Eigen;
namespace fs = std::filesystem;

int main(int argc, char** argv)
{
    // Capture des erreurs globales
    try {
        // -------------------------------------------------------------------------
        // 1. Initialisation et Vérifications
        // -------------------------------------------------------------------------
        // [NOUVEAU] Initialisation de la graine aléatoire pour le bruit
        std::srand(static_cast<unsigned int>(std::time(nullptr)));

        if (argc < 2) {
            cerr << "Erreur : Fichier d'entree manquant." << endl;
            cout << "Usage : ./ns_solver_2d input/input.toml" << endl;
            return 1;
        }

        auto start_cpu = std::chrono::high_resolution_clock::now();

        // Lecture des paramètres
        cout << "--------------------------------------------------" << endl;
        cout << "  SOLVEUR NAVIER-STOKES 2D - VON KARMAN  " << endl;
        cout << "--------------------------------------------------" << endl;
        
        DataFile* df = new DataFile(argv[1]);
        
        // Création du dossier de sortie
        string resultsPath = df->Get_results();
        if (!fs::exists(resultsPath)) {
            cout << "Creation du dossier de resultats : " << resultsPath << endl;
            fs::create_directories(resultsPath);
        }

        // -------------------------------------------------------------------------
        // 2. Initialisation des Objets Physiques
        // -------------------------------------------------------------------------
        Function* fct = new Function(df);
        MACgrid* grid = new MACgrid(fct, df);
        Laplacian* lap = new Laplacian(fct, df, grid);

        // [NOUVEAU] Système de Particules
        ParticleSystem* particles = new ParticleSystem(grid);
        particles->InitParticles(40); // On lance 40 particules

        cout << "Construction de la matrice de Poisson (Cholesky)... ";
        lap->BuildMatrix();
        cout << "[OK]" << endl;

        // -------------------------------------------------------------------------
        // 3. Choix du Schéma Temporel
        // -------------------------------------------------------------------------
        TimeScheme* time_scheme = nullptr;
        string scheme_name = df->Get_scheme();

        if (scheme_name == "ExplicitEuler" || scheme_name == "Euler") {
            cout << "Schema Temporel : Euler Explicite" << endl;
            time_scheme = new EulerScheme(df, lap, grid);
        } else {
            cout << "/!\\ Schema '" << scheme_name << "' inconnu. Fallback -> Euler Explicite." << endl;
            time_scheme = new EulerScheme(df, lap, grid);
        }

        // -------------------------------------------------------------------------
        // 4. Vérification de Stabilité (CFL)
        // -------------------------------------------------------------------------
        double h_min = std::min(df->Get_hx(), df->Get_hy());
        double nu = df->Get_nu();
        double dt = df->Get_dt();
        
        double cfl_diff = nu * dt / (h_min * h_min);
        cout << "CFL Diffusion (nu*dt/h^2) : " << cfl_diff << endl;
        if (cfl_diff > 0.25) {
            cout << "/!\\ ATTENTION : Risque d'instabilite (Viscosite)." << endl;
        }

        // -------------------------------------------------------------------------
        // 5. Boucle Temporelle
        // -------------------------------------------------------------------------
        double t = df->Get_t0();
        double t_final = df->Get_tfinal();
        int iteration = 0;

        // Sauvegarde tous les ~1% du temps total
        int total_iters = static_cast<int>((t_final - t) / dt);
        int save_freq = std::max(1, total_iters / 100); 
        
        cout << "Simulation de t=" << t << " a t=" << t_final << endl;
        cout << "Iterations prevues : " << total_iters << endl;
        cout << "Sauvegarde toutes les " << save_freq << " iterations." << endl;
        cout << "--------------------------------------------------" << endl;

        // Fichier historique
        string historyPath = resultsPath + "/history.dat";
        ofstream history_file(historyPath.c_str()); 
        history_file << "t kinetic_energy max_velocity" << endl;

        // Sauvegarde initiale (t=0)
        time_scheme->SaveSolution(0);
        particles->Save(0, resultsPath); // Save particules initiales

        while (t < t_final)
        {
            // --- 1. Avance d'un pas FLUIDE ---
            time_scheme->Advance();
            t = time_scheme->GetTime(); 
            iteration++;

            // --- 2. Avance d'un pas PARTICULES ---
            particles->Advance(dt);

            // --- 3. Calculs de contrôle (Stats) ---
            if (iteration % 10 == 0) { // On ne calcule pas les stats à chaque pas pour aller vite
                const VectorXd& U = grid->GetU();
                const VectorXd& V = grid->GetV();

                double max_u = U.lpNorm<Infinity>();
                double max_v = V.lpNorm<Infinity>();
                double max_vel = std::max(max_u, max_v);

                double cell_vol = df->Get_hx() * df->Get_hy();
                double kinetic_energy = 0.5 * (U.squaredNorm() + V.squaredNorm()) * cell_vol;

                history_file << t << " " << kinetic_energy << " " << max_vel << endl;

                // Affichage Console
                if (iteration % save_freq == 0) {
                    float progress = (float)iteration / total_iters * 100.0;
                    cout << "[ " << fixed << setprecision(1) << setw(5) << progress << "% ] "
                         << "t: " << t << " | Ek: " << scientific << setprecision(3) << kinetic_energy << endl;
                    
                    // --- 4. SAUVEGARDE ---
                    time_scheme->SaveSolution(iteration);
                    particles->Save(iteration, resultsPath);
                }
            }
        }
        
        history_file.close();

        // -------------------------------------------------------------------------
        // 6. Finalisation
        // -------------------------------------------------------------------------
        auto end_cpu = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_cpu - start_cpu;

        cout << "--------------------------------------------------" << endl;
        cout << "Simulation terminee." << endl;
        cout << "Temps de calcul : " << fixed << setprecision(2) << elapsed.count() << " s" << endl;
        cout << "Resultats dans : " << resultsPath << endl;
        cout << "--------------------------------------------------" << endl;

        // Nettoyage mémoire
        delete particles;
        delete time_scheme;
        delete lap;
        delete grid;
        delete fct;
        delete df;

    } catch (const std::exception& e) {
        cerr << "\n!!! ERREUR CRITIQUE !!!" << endl;
        cerr << e.what() << endl;
        return 1;
    }

    return 0;
}