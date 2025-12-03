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
    // On redimensionne la matrice
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

    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k = _grid->GetPIndex(i, j);

            // --- 1. SORTIE (Droite) : Pression fixée à 0 (Dirichlet) ---
            // Condition : j = Nx - 1
            if (j == Nx - 1) {
                // Équation triviale : 1 * P_k = 0
                triplets.push_back(Eigen::Triplet<double>(k, k, 1.0));
                continue; // On passe au point suivant, pas de voisins connectés
            }

            // --- 2. POINTS INTERNES & AUTRES BORDS ---
            double diag = 2.0 * (Cx + Cy);

            // OUEST (Gauche)
            if (j > 0) {
                int k_W = _grid->GetPIndex(i, j - 1);
                triplets.push_back(Eigen::Triplet<double>(k, k_W, -Cx));
            } else {
                // Mur Gauche / Entrée : Neumann (dP/dx = 0)
                // Le gradient de pression normal est nul à l'entrée car la vitesse est imposée.
                diag -= Cx;
            }

            // EST (Droite)
            // On vérifie si le voisin de droite est un point interne.
            // Si j+1 correspond à la frontière de sortie (j+1 == Nx-1), P_{est} vaut 0.
            // Le terme -Cx * P_{est} s'annule. On ne met donc PAS de coefficient dans la matrice.
            // Cela préserve la symétrie car la ligne (Nx-1) n'a pas non plus de lien vers (Nx-2).
            if (j + 1 < Nx - 1) {
                int k_E = _grid->GetPIndex(i, j + 1);
                triplets.push_back(Eigen::Triplet<double>(k, k_E, -Cx));
            }

            // SUD (Bas)
            if (i > 0) {
                int k_S = _grid->GetPIndex(i - 1, j);
                triplets.push_back(Eigen::Triplet<double>(k, k_S, -Cy));
            } else {
                // Mur Bas : Neumann (dP/dy = 0)
                diag -= Cy;
            }

            // NORD (Haut)
            if (i < Ny - 1) {
                int k_N = _grid->GetPIndex(i + 1, j);
                triplets.push_back(Eigen::Triplet<double>(k, k_N, -Cy));
            } else {
                // Mur Haut : Neumann (dP/dy = 0)
                diag -= Cy;
            }

            // Ajout de la diagonale
            triplets.push_back(Eigen::Triplet<double>(k, k, diag));
        }
    }

    // Construction de la matrice Sparse
    _H.setFromTriplets(triplets.begin(), triplets.end());

    // Factorisation (LLT pour matrice symétrique définie positive)
    _solver.compute(_H);
}

Eigen::VectorXd Laplacian::ComputeDivergence(const Eigen::VectorXd& U, const Eigen::VectorXd& V)
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();

    Eigen::VectorXd div(Nx * Ny);
    div.setZero();

    for (int i = 0; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k = _grid->GetPIndex(i, j);

            int k_u_E = _grid->GetUIndex(i, j + 1);
            int k_u_W = _grid->GetUIndex(i, j);
            int k_v_N = _grid->GetVIndex(i + 1, j);
            int k_v_S = _grid->GetVIndex(i, j);

            double du_dx = (U(k_u_E) - U(k_u_W)) / hx;
            double dv_dy = (V(k_v_N) - V(k_v_S)) / hy;

            div(k) = du_dx + dv_dy;
        }
    }
    return div;
}

void Laplacian::Solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& p_sol)
{
    // Le système est H * P = - div(u) / dt
    // Donc b = -rhs
    Eigen::VectorXd b = -rhs;

    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();

    // IMPORTANT : Imposer le second membre à 0 pour la frontière de sortie.
    // Car l'équation y est : 1 * P = 0.
    for (int i = 0; i < Ny; ++i) {
        int k_out = _grid->GetPIndex(i, Nx - 1);
        b(k_out) = 0.0;
    }

    // Résolution
    p_sol = _solver.solve(b);
}

void Laplacian::ComputeGradient(const Eigen::VectorXd& p, Eigen::VectorXd& gradPx, Eigen::VectorXd& gradPy)
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();

    // Redimensionnement
    gradPx.resize((Nx + 1) * Ny);
    gradPy.resize(Nx * (Ny + 1));
    gradPx.setZero();
    gradPy.setZero();

    // GradPx sur les faces U (internes : de j=1 à Nx-1)
    for (int i = 0; i < Ny; ++i) {
        for (int j = 1; j < Nx; ++j) {
            int k_u = _grid->GetUIndex(i, j);
            int k_p_E = _grid->GetPIndex(i, j);
            int k_p_W = _grid->GetPIndex(i, j - 1);

            gradPx(k_u) = (p(k_p_E) - p(k_p_W)) / hx;
        }
    }

    // GradPy sur les faces V (internes : de i=1 à Ny-1)
    for (int i = 1; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k_v = _grid->GetVIndex(i, j);
            int k_p_N = _grid->GetPIndex(i, j);
            int k_p_S = _grid->GetPIndex(i - 1, j);

            gradPy(k_v) = (p(k_p_N) - p(k_p_S)) / hy;
        }
    }
}