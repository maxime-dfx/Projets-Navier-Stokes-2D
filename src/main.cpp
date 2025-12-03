#include <iostream>
#include <string>
#include <cmath>
#include <fstream>  
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "DataFile.h"
#include "Function.h"
#include "MACgrid.h"
#include "Laplacian.h"
#include "TimeScheme.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
    // 1. Vérification des arguments
    if (argc < 2) {
        cout << "Usage : ./bin/ns_solver_2d input/input.toml" << endl;
        return 1;
    }

    // 2. Initialisation des objets
    // Lecture des paramètres
    DataFile* df = new DataFile(argv[1]);
    
    // Fonctions mathématiques (conditions initiales)
    Function* fct = new Function(df);
    
    // Grille MAC (s'initialise toute seule via Function désormais)
    MACgrid* grid = new MACgrid(fct, df);
    
    // Opérateur Laplacien
    Laplacian* lap = new Laplacian(fct, df, grid);
    
    // Construction de la matrice de Poisson (une seule fois au début)
    cout << "Construction de la matrice de Poisson..." << endl;
    lap->BuildMatrix();
    cout << "Matrice construite." << endl;

    // =========================================================================
    // 3. Choix du Schéma Temporel (CORRECTION CRITIQUE ICI)
    // =========================================================================
    std::string scheme_name = df->Get_scheme();
    TimeScheme* time_scheme = nullptr;

    if (scheme_name == "ExplicitEuler" || scheme_name == "Euler") {
        cout << "Schema choisi : Euler Explicite" << endl;
        time_scheme = new EulerScheme(df, lap, grid);
    } 
    else {
        // Fallback par défaut si le nom est inconnu ou mal écrit dans le .toml
        cout << "ATTENTION : Schema '" << scheme_name << "' inconnu." << endl;
        cout << "Utilisation de Euler Explicite par defaut." << endl;
        time_scheme = new EulerScheme(df, lap, grid);
    }

    // =========================================================================
    // 4. Boucle Temporelle
    // =========================================================================
    double t = df->Get_t0();
    double t_final = df->Get_tfinal();
    int iteration = 0;
    int save_freq = 10;

    // Création du fichier de log ---
    string resultsPath = df->Get_results();
    string historyPath = resultsPath + "/history.dat";
    ofstream history_file(historyPath.c_str()); 

    // En-tête du fichier pour Excel/Gnuplot
    history_file << "t kinetic_energy max_velocity" << endl;

    cout << "Lancement de la simulation..." << endl;
    time_scheme->SaveSolution(0);

    while (t < t_final)
    {
        time_scheme->Advance();
        t = time_scheme->GetTime(); 
        iteration++;

        // On le fait à chaque pas (ou tous les 10 pas) pour avoir une belle courbe
        const VectorXd& U = grid->GetU();
        const VectorXd& V = grid->GetV();
        
        // 1. Energie Cinétique Totale (Approximation discrete)
        // E_k = 0.5 * rho * dx * dy * sum(u^2 + v^2)
        // On simplifie ici par la somme des carrés (suffisant pour voir la forme de la courbe)
        double square_sum = U.squaredNorm() + V.squaredNorm();
        double kinetic_energy = 0.5 * square_sum * df->Get_hx() * df->Get_hy();
        double p_min = grid->GetP().minCoeff();
        double p_max = grid->GetP().maxCoeff();
        double p_amplitude = p_max - p_min;

        if (iteration % 100 == 0) {
        cout << "Iter: " << iteration 
         << " | P_amplitude: " << p_amplitude << endl;
}
        // 2. Vitesse max (Norme infinie)
        double max_u = U.lpNorm<Infinity>();
        double max_v = V.lpNorm<Infinity>();
        double max_vel = std::max(max_u, max_v);

        // Écriture dans le fichier
        history_file << t << " " << kinetic_energy << " " << max_vel << endl;

        // Affichage console (optionnel)
        if (iteration % 100 == 0) {
            cout << "Iter: " << iteration << " t=" << t << " Ek=" << kinetic_energy << endl;
        }

        if (iteration % save_freq == 0) {
            time_scheme->SaveSolution(iteration);
        }
    }
    
    history_file.close();
    cout << "--------------------------------------------------" << endl;
    cout << "Simulation terminee avec succes." << endl;
    cout << "Resultats disponibles dans : " << df->Get_results() << endl;
    cout << "--------------------------------------------------" << endl;

    // 5. Nettoyage de la mémoire
    delete time_scheme;
    delete lap;
    delete grid;
    delete fct;
    delete df;

    return 0;
}