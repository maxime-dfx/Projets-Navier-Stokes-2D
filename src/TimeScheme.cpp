#include "TimeScheme.h"
#include "Function.h" 
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm> 
#include <locale>

using namespace Eigen;
using namespace std;

TimeScheme::TimeScheme(DataFile* data_file, Laplacian* lap, MACgrid* grid) :
_df(data_file), _lap(lap), _grid(grid), _t(data_file->Get_t0())
{}

TimeScheme::~TimeScheme() {}

EulerScheme::EulerScheme(DataFile* data_file, Laplacian* lap, MACgrid* grid) :
TimeScheme(data_file, lap, grid)
{}

// Application des conditions aux limites
void TimeScheme::ApplyBoundaryConditions()
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    
    VectorXd U = _grid->GetU();
    VectorXd V = _grid->GetV();
    Function* fct = _grid->GetFunction();

    // MURS VERTICAUX (Agit sur U)
    for (int i = 0; i < Ny; ++i) {
        int k_left = _grid->GetUIndex(i, 0);
        if (fct->IsDirichletLeft()) U(k_left) = fct->GetLeftU(_grid->GetUcoord(i, 0)(1)); 
        else U(k_left) = U(_grid->GetUIndex(i, 1)); 

        int k_right = _grid->GetUIndex(i, Nx);
        if (fct->IsDirichletRight()) U(k_right) = fct->GetRightU(_grid->GetUcoord(i, Nx)(1));
        else U(k_right) = U(_grid->GetUIndex(i, Nx - 1));
    }

    // MURS HORIZONTAUX (Agit sur V)
    for (int j = 0; j < Nx; ++j) {
        int k_bott = _grid->GetVIndex(0, j);
        if (fct->IsDirichletBottom()) V(k_bott) = fct->GetBottomV(_grid->GetVcoord(0, j)(0));
        else V(k_bott) = V(_grid->GetVIndex(1, j));

        int k_top = _grid->GetVIndex(Ny, j);
        if (fct->IsDirichletTop()) V(k_top) = fct->GetTopV(_grid->GetVcoord(Ny, j)(0));
        else V(k_top) = V(_grid->GetVIndex(Ny - 1, j));
    }
    _grid->SetU(U);
    _grid->SetV(V);
}

// =========================================================================
// Fonction utilitaire pour le schéma Upwind (A mettre avant Advance)
// =========================================================================
// Stabilise le calcul en regardant "d'où vient le vent".
// Si la vitesse est positive, on prend l'information de gauche (Back).
// Si la vitesse est négative, on prend l'information de droite (Front).
inline double Upwind(double vel, double val_minus, double val_center, double val_plus, double inv_h) {
    if (vel > 0) return vel * (val_center - val_minus) * inv_h;
    else return vel * (val_plus - val_center) * inv_h;
}

// =========================================================================
// EULER SCHEME : ADVANCE
// =========================================================================

void EulerScheme::Advance()
{
    // 0. Récupération des paramètres
    double dt = _df->Get_dt();
    double nu = _df->Get_nu();
    double rho = _df->Get_rho();
    
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();

    // Pré-calcul des inverses pour la vitesse
    double odx = 1.0 / hx;
    double ody = 1.0 / hy;
    double odx2 = 1.0 / (hx * hx);
    double ody2 = 1.0 / (hy * hy);

    // 1. Application des BC (Important pour avoir des valeurs fantômes correctes)
    ApplyBoundaryConditions();

    // Copie des champs actuels (u^n, v^n)
    VectorXd u_n = _grid->GetU();
    VectorXd v_n = _grid->GetV();

    // Vecteurs prédiction (u*, v*) initialisés à u^n, v^n
    VectorXd u_star = u_n;
    VectorXd v_star = v_n;
    
    Function* fct = _grid->GetFunction();

    // -----------------------------------------------------------------------
    // 2. ÉTAPE DE PRÉDICTION : u* = u^n + dt * (Diffusion - Advection)
    // -----------------------------------------------------------------------

    // --- A. Calcul pour U (Faces verticales internes : j=1 à Nx-1) ---
    for (int i = 0; i < Ny; ++i) {
        for (int j = 1; j < Nx; ++j) { 
            int k = _grid->GetUIndex(i, j);
            double u_curr = u_n(k);

            // Voisins Est/Ouest (existent toujours car j est interne)
            double u_E = u_n(_grid->GetUIndex(i, j + 1));
            double u_W = u_n(_grid->GetUIndex(i, j - 1));

            // Voisins Nord/Sud (Gestion Ghost Cells pour Dirichlet)
            double u_N, u_S;
            
            // SUD (Bas)
            if (i > 0) u_S = u_n(_grid->GetUIndex(i - 1, j));
            else { 
                // Si mur Dirichlet (vitesse=0), le point fantôme est l'opposé : -u_curr
                // (moyenne au mur = 0). Sinon Neumann (glissement) : u_curr.
                u_S = fct->IsDirichletBottom() ? -u_curr : u_curr; 
            }
            
            // NORD (Haut)
            if (i < Ny - 1) u_N = u_n(_grid->GetUIndex(i + 1, j));
            else {
                u_N = fct->IsDirichletTop() ? -u_curr : u_curr;
            }

            // --- DIFFUSION (Laplacien) ---
            double diffusion = nu * ((u_E - 2*u_curr + u_W)*odx2 + (u_N - 2*u_curr + u_S)*ody2);

            // --- ADVECTION (Upwind) ---
            
            // Terme 1 : u * du/dx
            double adv_x = Upwind(u_curr, u_W, u_curr, u_E, odx);

            // Terme 2 : v * du/dy (Moyenne de 4 V voisins pour avoir v au point u)
            double v_avg = 0.25 * (
                (i < Ny ? v_n(_grid->GetVIndex(i + 1, j)) : 0.0) +     // NE
                (i < Ny ? v_n(_grid->GetVIndex(i + 1, j - 1)) : 0.0) + // NW
                v_n(_grid->GetVIndex(i, j)) +                          // SE
                v_n(_grid->GetVIndex(i, j - 1))                        // SW
            );
            
            double adv_y = Upwind(v_avg, u_S, u_curr, u_N, ody);

            // Bilan U
            u_star(k) = u_curr + dt * (diffusion - (adv_x + adv_y));
        }
    }

    // --- B. Calcul pour V (Faces horizontales internes : i=1 à Ny-1) ---
    for (int i = 1; i < Ny; ++i) { 
        for (int j = 0; j < Nx; ++j) {
            int k = _grid->GetVIndex(i, j);
            double v_curr = v_n(k);

            // Voisins Nord/Sud (existent toujours)
            double v_N = v_n(_grid->GetVIndex(i + 1, j));
            double v_S = v_n(_grid->GetVIndex(i - 1, j));

            // Voisins Est/Ouest (Ghost Cells)
            double v_E, v_W;

            // OUEST (Gauche)
            if (j > 0) v_W = v_n(_grid->GetVIndex(i, j - 1));
            else {
                v_W = fct->IsDirichletLeft() ? -v_curr : v_curr;
            }

            // EST (Droite)
            if (j < Nx - 1) v_E = v_n(_grid->GetVIndex(i, j + 1));
            else {
                v_E = fct->IsDirichletRight() ? -v_curr : v_curr;
            }

            // --- DIFFUSION ---
            double diffusion = nu * ((v_E - 2*v_curr + v_W)*odx2 + (v_N - 2*v_curr + v_S)*ody2);

            // --- ADVECTION (Upwind) ---

            // Terme 1 : v * dv/dy
            double adv_y = Upwind(v_curr, v_S, v_curr, v_N, ody);

            // Terme 2 : u * dv/dx (Moyenne de 4 U voisins pour avoir u au point v)
            double u_avg = 0.25 * (
                (j < Nx ? u_n(_grid->GetUIndex(i, j + 1)) : 0.0) +     // NE
                u_n(_grid->GetUIndex(i, j)) +                          // NW
                (j < Nx ? u_n(_grid->GetUIndex(i - 1, j + 1)) : 0.0) + // SE
                u_n(_grid->GetUIndex(i - 1, j))                        // SW
            );

            double adv_x = Upwind(u_avg, v_W, v_curr, v_E, odx);

            // Bilan V
            v_star(k) = v_curr + dt * (diffusion - (adv_x + adv_y));
        }
    }

    // -----------------------------------------------------------------------
    // 3. PROJECTION (Pression)
    // -----------------------------------------------------------------------
    
    // Div(u*)
    VectorXd div_u_star = _lap->ComputeDivergence(u_star, v_star);
    
    // Solve H*p = (rho/dt) * div
    VectorXd rhs = (rho / dt) * div_u_star;

    // --- DEBUG: Vérification de l'activité du solveur ---
    // Si cette valeur est proche de 0 (ex: 1e-16), alors P restera 0.
    // Si cette valeur est grande, P doit apparaitre.
    // std::cout << "Time: " << _t << " | RHS Norm: " << rhs.norm() << std::endl;
    // ----------------------------------------------------

    VectorXd p_next;
    _lap->Solve(rhs, p_next);

    // -----------------------------------------------------------------------
    // 4. CORRECTION : u^{n+1} = u* - (dt/rho) * Grad(p)
    // -----------------------------------------------------------------------
    VectorXd gradPx, gradPy;
    _lap->ComputeGradient(p_next, gradPx, gradPy);

    // On modifie les vecteurs prédiction pour obtenir u_next
    // Attention : On ne touche qu'aux faces internes, les bords sont gérés par ApplyBoundaryConditions
    
    VectorXd u_next = u_star;
    VectorXd v_next = v_star;

    // Correction U
    for (int i = 0; i < Ny; ++i) {
        for (int j = 1; j < Nx; ++j) {
            int k = _grid->GetUIndex(i, j);
            u_next(k) -= (dt / rho) * gradPx(k);
        }
    }

    // Correction V
    for (int i = 1; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k = _grid->GetVIndex(i, j);
            v_next(k) -= (dt / rho) * gradPy(k);
        }
    }

    // Envoi à la grille
    _grid->SetU(u_next); 
    _grid->SetV(v_next);
    _grid->SetP(p_next); // Mise à jour de la pression pour l'affichage

    // -----------------------------------------------------------------------
    // 5. Finalisation
    // -----------------------------------------------------------------------
    // On ré-applique les conditions aux limites pour le pas suivant
    ApplyBoundaryConditions();

    _t += dt;
}

void TimeScheme::SaveSolution(int n_iteration)
{
    stringstream ss_dat;
    string resultsPath = _df->Get_results();
    ss_dat << resultsPath << "/sol_" << n_iteration << ".dat";
    ofstream dat(ss_dat.str());
    if (dat.is_open()) {
        dat.imbue(std::locale("C"));
        int Nx = _df->Get_Nx();
        int Ny = _df->Get_Ny();
        const VectorXd& P = _grid->GetP();
        const VectorXd& U = _grid->GetU();
        const VectorXd& V = _grid->GetV();

        for (int i = 0; i < Ny; ++i) {
            for (int j = 0; j < Nx; ++j) {
                // Coordonnées du centre
                double x = _df->Get_xmin() + (j + 0.5) * _df->Get_hx();
                double y = _df->Get_ymin() + (i + 0.5) * _df->Get_hy();
                
                // Valeurs
                double p_val = P(_grid->GetPIndex(i, j));
                double u_val = 0.5 * (U(_grid->GetUIndex(i, j)) + U(_grid->GetUIndex(i, j+1)));
                double v_val = 0.5 * (V(_grid->GetVIndex(i, j)) + V(_grid->GetVIndex(i+1, j)));

                // Écriture : X Y P U V
                dat << x << " " << y << " " << p_val << " " << u_val << " " << v_val << "\n";
            }
            // IMPORTANT pour Gnuplot : Saut de ligne à chaque fin de ligne Y
            // pour qu'il comprenne que c'est une grille 2D
            dat << "\n"; 
        }
    }
    dat.close();
}