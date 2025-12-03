#include "Laplacian.h"
#include "DataFile.h" 
#include "MACgrid.h"
#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;

Laplacian::Laplacian(Function* function, DataFile* data_file, MACgrid*& grid) :
_fct(function), _df(data_file), _grid(grid)
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    _H.resize(Nx * Ny, Nx * Ny);
}

void Laplacian::BuildMatrix()
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();
    double Cx = 1.0 / (hx * hx);
    double Cy = 1.0 / (hy * hy);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(5 * Nx * Ny); 

    // Point fixé pour rendre la matrice inversible (Pinning)
    int k_pinned = _grid->GetPIndex(0, 0);

    for (int i = 0; i < Ny; ++i) {       
        for (int j = 0; j < Nx; ++j) {   
            int k = _grid->GetPIndex(i, j);

            // --- 1. PINNING : On force 1.0 sur la diagonale pour le point (0,0) ---
            if (k == k_pinned) {
                triplets.push_back(Triplet<double>(k, k, 1.0));
                continue;
            }

            // --- 2. CHANGEMENT DE SIGNE ---
            // On construit la matrice opposée (-Laplacien) pour qu'elle soit Positive Définie
            // Diagonale POSITIVE
            double diag = 2.0 * (Cx + Cy);

            // Voisins NEGATIFS
            
            // Ouest
            if (j > 0) {
                int k_W = _grid->GetPIndex(i, j-1);
                if (k_W != k_pinned) triplets.push_back(Triplet<double>(k, k_W, -Cx));
            } else {
                // Mur Neumann : La dérivée nulle supprime le terme,
                // MAIS on doit réduire la diagonale pour refléter la perte du voisin
                // (Sinon on simule un Dirichlet implicite = 0)
                diag -= Cx; 
            }
            
            // Est
            if (j < Nx - 1) {
                int k_E = _grid->GetPIndex(i, j+1);
                if (k_E != k_pinned) triplets.push_back(Triplet<double>(k, k_E, -Cx));
            } else {
                diag -= Cx; // Mur Neumann
            }

            // Sud
            if (i > 0) {
                int k_S = _grid->GetPIndex(i-1, j);
                if (k_S != k_pinned) triplets.push_back(Triplet<double>(k, k_S, -Cy));
            } else {
                diag -= Cy; // Mur Neumann
            }

            // Nord
            if (i < Ny - 1) {
                int k_N = _grid->GetPIndex(i+1, j);
                if (k_N != k_pinned) triplets.push_back(Triplet<double>(k, k_N, -Cy));
            } else {
                diag -= Cy; // Mur Neumann
            }

            triplets.push_back(Triplet<double>(k, k, diag));
        }
    }
    _H.setFromTriplets(triplets.begin(), triplets.end());

    _solver.analyzePattern(_H);
    _solver.compute(_H);
    
    if(_solver.info() != Success) {
        cout << "ERREUR FATALE: Factorisation Cholesky impossible." << endl;
        exit(1);
    }
}

Eigen::VectorXd Laplacian::ComputeDivergence(const Eigen::VectorXd& U, const Eigen::VectorXd& V)
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();
    
    Eigen::VectorXd div(Nx * Ny);
    
    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k = _grid->GetPIndex(i, j);
            
            double u_right = U(_grid->GetUIndex(i, j + 1));
            double u_left  = U(_grid->GetUIndex(i, j));
            double v_top   = V(_grid->GetVIndex(i + 1, j));
            double v_bott  = V(_grid->GetVIndex(i, j));
            
            div(k) = (u_right - u_left)/hx + (v_top - v_bott)/hy;
        }
    }
    return div;
}

void Laplacian::Solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& p_sol)
{
    // IMPORTANT : Comme on a inversé le signe de la matrice (M = -Laplacien),
    // on doit inverser le signe du RHS pour résoudre M*p = -rhs
    // Sinon la pression sera inversée (Haute pression au lieu de basse).
    
    Eigen::VectorXd b = -rhs;

    // Pour le point pinné (0,0), l'équation est 1*p = 0.
    // On force donc b(0) à 0 explicitement.
    int k_pinned = _grid->GetPIndex(0, 0);
    b(k_pinned) = 0.0;

    p_sol = _solver.solve(b);
    
    // On remet la moyenne à 0 pour la propreté physique (optionnel mais recommandé)
    double p_moy = p_sol.sum() / p_sol.size();
    for(int i=0; i<p_sol.size(); ++i) p_sol(i) -= p_moy;
}

void Laplacian::ComputeGradient(const Eigen::VectorXd& p, Eigen::VectorXd& gradPx, Eigen::VectorXd& gradPy)
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();

    gradPx.resize((Nx + 1) * Ny);
    gradPy.resize(Nx * (Ny + 1));
    gradPx.setZero();
    gradPy.setZero();

    // GradPx sur les faces U (internes)
    for (int i = 0; i < Ny; ++i) {
        for (int j = 1; j < Nx; ++j) { 
            int k_u = _grid->GetUIndex(i, j);
            int k_p_right = _grid->GetPIndex(i, j);
            int k_p_left  = _grid->GetPIndex(i, j - 1);
            gradPx(k_u) = (p(k_p_right) - p(k_p_left)) / hx;
        }
    }

    // GradPy sur les faces V (internes)
    for (int i = 1; i < Ny; ++i) { 
        for (int j = 0; j < Nx; ++j) {
            int k_v = _grid->GetVIndex(i, j);
            int k_p_top  = _grid->GetPIndex(i, j);
            int k_p_bott = _grid->GetPIndex(i - 1, j);
            gradPy(k_v) = (p(k_p_top) - p(k_p_bott)) / hy;
        }
    }
}