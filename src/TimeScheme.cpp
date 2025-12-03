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

// =========================================================================
// APPLICATION DES CONDITIONS AUX LIMITES (Vitesse Normale / Étanchéité)
// =========================================================================
void TimeScheme::ApplyBoundaryConditions()
{
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();
    
    VectorXd U = _grid->GetU();
    VectorXd V = _grid->GetV();
    Function* fct = _grid->GetFunction();

    // 1. MURS VERTICAUX (On impose U, la vitesse normale traversante)
    for (int i = 0; i < Ny; ++i) {
        // Mur GAUCHE (Indice 0)
        int k_left = _grid->GetUIndex(i, 0);
        double y_left = _grid->GetUcoord(i, 0)(1); // Coord Y
        if (fct->IsDirichletLeft()) 
            U(k_left) = fct->GetLeftU_Normal(y_left); 
        else 
            U(k_left) = U(_grid->GetUIndex(i, 1)); 

        // Mur DROIT (Indice Nx)
        int k_right = _grid->GetUIndex(i, Nx);
        double y_right = _grid->GetUcoord(i, Nx)(1);
        if (fct->IsDirichletRight()) 
            U(k_right) = fct->GetRightU_Normal(y_right);
        else 
            U(k_right) = U(_grid->GetUIndex(i, Nx - 1));
    }

    // 2. MURS HORIZONTAUX (On impose V, la vitesse normale traversante)
    for (int j = 0; j < Nx; ++j) {
        // Mur BAS (Indice 0)
        int k_bott = _grid->GetVIndex(0, j);
        double x_bott = _grid->GetVcoord(0, j)(0); // Coord X
        if (fct->IsDirichletBottom()) 
            V(k_bott) = fct->GetBottomV_Normal(x_bott);
        else 
            V(k_bott) = V(_grid->GetVIndex(1, j));

        // Mur HAUT (Indice Ny)
        int k_top = _grid->GetVIndex(Ny, j);
        double x_top = _grid->GetVcoord(Ny, j)(0);
        if (fct->IsDirichletTop()) 
            V(k_top) = fct->GetTopV_Normal(x_top);
        else 
            V(k_top) = V(_grid->GetVIndex(Ny - 1, j));
    }
    _grid->SetU(U);
    _grid->SetV(V);
}

// Helper Upwind inchangé
inline double Upwind(double vel, double val_minus, double val_center, double val_plus, double inv_h) {
    if (vel > 0) return vel * (val_center - val_minus) * inv_h;
    else return vel * (val_plus - val_center) * inv_h;
}

// =========================================================================
// EULER SCHEME : ADVANCE
// =========================================================================

void EulerScheme::Advance()
{
    double dt = _df->Get_dt();
    double nu = _df->Get_nu();
    double rho = _df->Get_rho();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();
    int Nx = _df->Get_Nx();
    int Ny = _df->Get_Ny();

    double odx = 1.0 / hx; double ody = 1.0 / hy;
    double odx2 = 1.0 / (hx*hx); double ody2 = 1.0 / (hy*hy);

    // 1. Appliquer BC Normales
    ApplyBoundaryConditions();

    VectorXd u_n = _grid->GetU();
    VectorXd v_n = _grid->GetV();
    VectorXd u_star = u_n;
    VectorXd v_star = v_n;
    
    Function* fct = _grid->GetFunction();

    // -----------------------------------------------------------------------
    // 2. PRÉDICTION U (Interne)
    // -----------------------------------------------------------------------
    for (int i = 0; i < Ny; ++i) {
        for (int j = 1; j < Nx; ++j) { 
            int k = _grid->GetUIndex(i, j);
            double u_curr = u_n(k);
            double x_curr = _grid->GetUcoord(i,j)(0);
            
            // --- DIFFUSION (Viscosité) ---
            double u_E = u_n(_grid->GetUIndex(i, j + 1));
            double u_W = u_n(_grid->GetUIndex(i, j - 1));
            double u_N, u_S;
            
            // Gestion Ghost Cells VERTICALES pour U (Frottement sur mur haut/bas)
            // U est tangent au mur du bas et du haut.
            
            // SUD (Mur Bas)
            if (i > 0) u_S = u_n(_grid->GetUIndex(i - 1, j));
            else { 
                if (fct->IsDirichletBottom()) {
                   // Formule: (u_ghost + u_curr)/2 = u_wall  => u_ghost = 2*u_wall - u_curr
                   double u_wall = fct->GetBottomU_Tangent(x_curr);
                   u_S = 2.0 * u_wall - u_curr;
                } else {
                   u_S = u_curr; // Neumann
                }
            }
            
            // NORD (Mur Haut - Ex: Lid Driven)
            if (i < Ny - 1) u_N = u_n(_grid->GetUIndex(i + 1, j));
            else {
                if (fct->IsDirichletTop()) {
                    double u_wall = fct->GetTopU_Tangent(x_curr);
                    u_N = 2.0 * u_wall - u_curr;
                } else {
                    u_N = u_curr;
                }
            }

            double diffusion = nu * ((u_E - 2*u_curr + u_W)*odx2 + (u_N - 2*u_curr + u_S)*ody2);

            // --- ADVECTION ---
            double adv_x = Upwind(u_curr, u_W, u_curr, u_E, odx);

            // Moyenne de V aux 4 coins pour l'avoir au centre de la face U
            double v_avg = 0.25 * (
                (i < Ny ? v_n(_grid->GetVIndex(i + 1, j)) : 0.0) +     
                (i < Ny ? v_n(_grid->GetVIndex(i + 1, j - 1)) : 0.0) + 
                v_n(_grid->GetVIndex(i, j)) +                          
                v_n(_grid->GetVIndex(i, j - 1))                        
            );
            double adv_y = Upwind(v_avg, u_S, u_curr, u_N, ody);

            u_star(k) = u_curr + dt * (diffusion - (adv_x + adv_y));
        }
    }

    // -----------------------------------------------------------------------
    // 3. PRÉDICTION V (Interne)
    // -----------------------------------------------------------------------
    for (int i = 1; i < Ny; ++i) { 
        for (int j = 0; j < Nx; ++j) {
            int k = _grid->GetVIndex(i, j);
            double v_curr = v_n(k);
            double y_curr = _grid->GetVcoord(i, j)(1);

            double v_N = v_n(_grid->GetVIndex(i + 1, j));
            double v_S = v_n(_grid->GetVIndex(i - 1, j));
            double v_E, v_W;

            // Gestion Ghost Cells HORIZONTALES pour V (Frottement sur murs latéraux)
            // V est tangent aux murs gauche/droite.

            // OUEST (Mur Gauche)
            if (j > 0) v_W = v_n(_grid->GetVIndex(i, j - 1));
            else {
                if (fct->IsDirichletLeft()) {
                    double v_wall = fct->GetLeftV_Tangent(y_curr);
                    v_W = 2.0 * v_wall - v_curr;
                } else {
                    v_W = v_curr;
                }
            }

            // EST (Mur Droit)
            if (j < Nx - 1) v_E = v_n(_grid->GetVIndex(i, j + 1));
            else {
                if (fct->IsDirichletRight()) {
                    double v_wall = fct->GetRightV_Tangent(y_curr);
                    v_E = 2.0 * v_wall - v_curr;
                } else {
                    v_E = v_curr;
                }
            }

            double diffusion = nu * ((v_E - 2*v_curr + v_W)*odx2 + (v_N - 2*v_curr + v_S)*ody2);

            double adv_y = Upwind(v_curr, v_S, v_curr, v_N, ody);

            double u_avg = 0.25 * (
                (j < Nx ? u_n(_grid->GetUIndex(i, j + 1)) : 0.0) +     
                u_n(_grid->GetUIndex(i, j)) +                          
                (j < Nx ? u_n(_grid->GetUIndex(i - 1, j + 1)) : 0.0) + 
                u_n(_grid->GetUIndex(i - 1, j))                        
            );
            double adv_x = Upwind(u_avg, v_W, v_curr, v_E, odx);

            v_star(k) = v_curr + dt * (diffusion - (adv_x + adv_y));
        }
    }

    // -----------------------------------------------------------------------
    // 4. PROJECTION & CORRECTION (Standard)
    // -----------------------------------------------------------------------
    
    VectorXd div_u_star = _lap->ComputeDivergence(u_star, v_star);
    VectorXd rhs = (rho / dt) * div_u_star;

    VectorXd p_next;
    _lap->Solve(rhs, p_next);

    VectorXd gradPx, gradPy;
    _lap->ComputeGradient(p_next, gradPx, gradPy);
    
    VectorXd u_next = u_star;
    VectorXd v_next = v_star;

    for (int i = 0; i < Ny; ++i) {
        for (int j = 1; j < Nx; ++j) {
            int k = _grid->GetUIndex(i, j);
            u_next(k) -= (dt / rho) * gradPx(k);
        }
    }

    for (int i = 1; i < Ny; ++i) {
        for (int j = 0; j < Nx; ++j) {
            int k = _grid->GetVIndex(i, j);
            v_next(k) -= (dt / rho) * gradPy(k);
        }
    }

    _grid->SetU(u_next); 
    _grid->SetV(v_next);
    _grid->SetP(p_next); 

    ApplyBoundaryConditions();
    _t += dt;
}

// SaveSolution reste inchangé par rapport à ta version, il est correct.
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
        double hx = _df->Get_hx();
        double hy = _df->Get_hy();
        
        const VectorXd& P = _grid->GetP();
        const VectorXd& U = _grid->GetU();
        const VectorXd& V = _grid->GetV();

        auto get_u_center = [&](int r, int c) {
            return 0.5 * (U(_grid->GetUIndex(r, c)) + U(_grid->GetUIndex(r, c+1)));
        };
        auto get_v_center = [&](int r, int c) {
            return 0.5 * (V(_grid->GetVIndex(r, c)) + V(_grid->GetVIndex(r+1, c)));
        };

        for (int i = 0; i < Ny; ++i) {
            for (int j = 0; j < Nx; ++j) {
                double x = _df->Get_xmin() + (j + 0.5) * hx;
                double y = _df->Get_ymin() + (i + 0.5) * hy;
                
                double p_val = P(_grid->GetPIndex(i, j));
                double u_val = get_u_center(i, j);
                double v_val = get_v_center(i, j);

                // Calcul vorticité approximé (différences finies)
                double dv_dx, du_dy;

                if (j == 0) dv_dx = (get_v_center(i, j+1) - v_val) / hx;
                else if (j == Nx - 1) dv_dx = (v_val - get_v_center(i, j-1)) / hx;
                else dv_dx = (get_v_center(i, j+1) - get_v_center(i, j-1)) / (2.0 * hx);

                if (i == 0) du_dy = (get_u_center(i+1, j) - u_val) / hy;
                else if (i == Ny - 1) du_dy = (u_val - get_u_center(i-1, j)) / hy;
                else du_dy = (get_u_center(i+1, j) - get_u_center(i-1, j)) / (2.0 * hy);

                double omega = dv_dx - du_dy;

                dat << x << " " << y << " " << p_val << " " << u_val << " " << v_val << " " << omega << "\n";
            }
            dat << "\n";
        }
    }
    dat.close();
}