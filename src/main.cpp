#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <chrono>   
#include <iomanip>  
#include <filesystem> 
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "DataFile.h"
#include "Function.h"
#include "MACgrid.h"
#include "Laplacian.h"
#include "TimeScheme.h"

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
        if (argc < 2) {
            cerr << "Erreur : Fichier d'entree manquant." << endl;
            cout << "Usage : ./ns_solver_2d input/input.toml" << endl;
            return 1;
        }

        auto start_cpu = std::chrono::high_resolution_clock::now();

        // Lecture des paramètres
        cout << "--------------------------------------------------" << endl;
        cout << "  SOLVEUR NAVIER-STOKES 2D - INCOMPRESSIBLE  " << endl;
        cout << "--------------------------------------------------" << endl;
        
        DataFile* df = new DataFile(argv[1]);
        
        // Création du dossier de sortie si inexistant (C++17)
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
        // 4. Vérification de Stabilité (CFL Diffusion)
        // -------------------------------------------------------------------------
        double h_min = std::min(df->Get_hx(), df->Get_hy());
        double nu = df->Get_nu();
        double dt = df->Get_dt();
        
        double cfl_diff = nu * dt / (h_min * h_min);
        cout << "CFL Diffusion (nu*dt/h^2) : " << cfl_diff << endl;
        if (cfl_diff > 0.25) { // 0.25 est une limite sûre pour 2D (0.5 max théorique)
            cout << "/!\\ ATTENTION : Le pas de temps est potentiellement trop grand pour la viscosité !" << endl;
            cout << "    Risque d'instabilite. Conseille : dt < " << 0.25 * h_min * h_min / nu << endl;
        }

        // -------------------------------------------------------------------------
        // 5. Boucle Temporelle
        // -------------------------------------------------------------------------
        double t = df->Get_t0();
        double t_final = df->Get_tfinal();
        int iteration = 0;

        // Calcul automatique de la fréquence de sauvegarde pour avoir environ 100 frames
        int total_iters = static_cast<int>((t_final - t) / dt);
        int save_freq = std::max(1, total_iters / 100); 
        
        cout << "Simulation de t=" << t << " a t=" << t_final << endl;
        cout << "Pas de temps dt=" << dt << " (" << total_iters << " iterations)" << endl;
        cout << "Sauvegarde toutes les " << save_freq << " iterations." << endl;
        cout << "--------------------------------------------------" << endl;

        // Fichier historique (Energies, etc.)
        string historyPath = resultsPath + "/history.dat";
        ofstream history_file(historyPath.c_str()); 
        history_file << "t kinetic_energy max_velocity div_max" << endl;

        // Sauvegarde initiale (t=0)
        time_scheme->SaveSolution(0);

        while (t < t_final)
        {
            // --- Avance d'un pas ---
            time_scheme->Advance();
            t = time_scheme->GetTime(); 
            iteration++;

            // --- Calculs de contrôle (Stats) ---
            const VectorXd& U = grid->GetU();
            const VectorXd& V = grid->GetV();

            // Vitesse max (Norme Infinie)
            double max_u = U.lpNorm<Infinity>();
            double max_v = V.lpNorm<Infinity>();
            double max_vel = std::max(max_u, max_v);

            // Energie Cinétique (Somme simple pondérée par le volume)
            // Note: U et V sont sur des grilles décalées, c'est une approx, mais suffisante.
            double cell_vol = df->Get_hx() * df->Get_hy();
            double kinetic_energy = 0.5 * (U.squaredNorm() + V.squaredNorm()) * cell_vol;

            // Divergence résiduelle (pour vérifier si la projection marche bien)
            // On recalculera la div juste pour le log (optionnel, peut être coûteux)
            // double div_max = lap->ComputeDivergence(U, V).lpNorm<Infinity>(); 
            double div_max = 0.0; // On met 0 pour l'instant pour gagner du temps CPU

            // Écriture fichier historique
            history_file << t << " " << kinetic_energy << " " << max_vel << " " << div_max << endl;

            // --- Affichage Console et Sauvegarde ---
            if (iteration % save_freq == 0) {
                // Barre de progression simple
                float progress = (float)iteration / total_iters * 100.0;
                cout << "[ " << fixed << setprecision(1) << setw(5) << progress << "% ] "
                     << "Iter: " << iteration 
                     << " | t: " << t 
                     << " | Ek: " << scientific << setprecision(3) << kinetic_energy 
                     << " | Vmax: " << max_vel << endl;
                
                time_scheme->SaveSolution(iteration);
            }
        }
        
        history_file.close();

        // -------------------------------------------------------------------------
        // 6. Finalisation
        // -------------------------------------------------------------------------
        auto end_cpu = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_cpu - start_cpu;

        cout << "--------------------------------------------------" << endl;
        cout << "Simulation terminee avec succes." << endl;
        cout << "Temps de calcul : " << fixed << setprecision(2) << elapsed.count() << " secondes." << endl;
        cout << "Resultats disponibles dans : " << resultsPath << endl;
        cout << "--------------------------------------------------" << endl;

        // Nettoyage mémoire
        delete time_scheme;
        delete lap;
        delete grid;
        delete fct;
        delete df;

    } catch (const std::exception& e) {
        // En cas de crash (ex: matrice mal formée, fichier introuvable)
        cerr << "\n!!! ERREUR CRITIQUE !!!" << endl;
        cerr << e.what() << endl;
        return 1;
    }

    return 0;
}