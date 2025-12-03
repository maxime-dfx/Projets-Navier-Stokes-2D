#include "MACgrid.h"
#include <iostream>
#include <cassert> 

using namespace std;
using namespace Eigen;

MACgrid::MACgrid(Function* function, DataFile* data_file) :
_fct(function), _df(data_file)
{
    _Nx = _df->Get_Nx();
    _Ny = _df->Get_Ny();
    
    // Allocations
    _p.resize(_Nx * _Ny);
    _U.resize((_Nx + 1) * _Ny);     
    _V.resize(_Nx * (_Ny + 1));     
    _is_solid_u.resize((_Nx + 1) * _Ny, false);
    _is_solid_v.resize(_Nx * (_Ny + 1), false);
    this->BuildObstacles(); // On construit l'obstacle dès le début
    _p.setZero();
    // _U et _V seront initialisés par les boucles ci-dessous, 
    // pas besoin de setZero() qui serait écrasé de toute façon.

    // --- Initialisation de U ---
    for (int i = 0; i < _Ny; ++i) {
        for (int j = 0; j <= _Nx; ++j) {
            int k = GetUIndex(i, j);
            VectorXd coord = GetUcoord(i, j);
            
            // Appel à la fonction corrigée (renvoie 1.0 désormais)
            _U(k) = _fct->InitialConditionU(coord(0), coord(1));
        }
    }

    // --- Initialisation de V ---
    for (int i = 0; i <= _Ny; ++i) {
        for (int j = 0; j < _Nx; ++j) {
            int k = GetVIndex(i, j);
            VectorXd coord = GetVcoord(i, j);
            
            // Appel à la fonction corrigée (renvoie 0.0)
            _V(k) = _fct->InitialConditionV(coord(0), coord(1));
        }
    }
}

int MACgrid::GetPIndex(int i, int j) const {
    return j + i * _Nx; 
}

int MACgrid::GetUIndex(int i, int j) const {
    return j + i * (_Nx + 1); 
}

int MACgrid::GetVIndex(int i, int j) const {
    return j + i * _Nx; 
}

VectorXd MACgrid::GetPcoord(int i, int j) const {
    VectorXd coord(2);
    coord(0) = _df->Get_xmin() + (j + 0.5) * _df->Get_hx(); 
    coord(1) = _df->Get_ymin() + (i + 0.5) * _df->Get_hy();
    return coord;
}

VectorXd MACgrid::GetUcoord(int i, int j) const {
    VectorXd coord(2);
    coord(0) = _df->Get_xmin() + j * _df->Get_hx(); 
    coord(1) = _df->Get_ymin() + (i + 0.5) * _df->Get_hy();
    return coord;
}

VectorXd MACgrid::GetVcoord(int i, int j) const {
    VectorXd coord(2);
    coord(0) = _df->Get_xmin() + (j + 0.5) * _df->Get_hx();
    coord(1) = _df->Get_ymin() + i * _df->Get_hy();
    return coord;
}

void MACgrid::BuildObstacles() {
    double x_min = _df->Get_xmin();
    double y_min = _df->Get_ymin();
    double hx = _df->Get_hx();
    double hy = _df->Get_hy();

    // --- DÉFINITION DE L'OBSTACLE (Cylindre) ---
    double cx = 0.5;  // Position X du centre (au début du tunnel)
    double cy = 0.5;  // Position Y du centre (milieu hauteur)
    double R  = 0.06; // Rayon du cylindre (assez gros pour perturber)
    double R2 = R * R;

    // Remplissage du masque pour U
    for (int i = 0; i < _Ny; ++i) {
        for (int j = 0; j <= _Nx; ++j) {
            // Coordonnées du point U(i,j)
            double x = x_min + j * hx;
            double y = y_min + (i + 0.5) * hy;

            // Équation du cercle : (x-cx)^2 + (y-cy)^2 < R^2
            if ( (x-cx)*(x-cx) + (y-cy)*(y-cy) <= R2 ) {
                _is_solid_u[GetUIndex(i, j)] = true;
                _U(GetUIndex(i, j)) = 0.0; // On met à 0 tout de suite
            }
        }
    }

    // Remplissage du masque pour V
    for (int i = 0; i <= _Ny; ++i) {
        for (int j = 0; j < _Nx; ++j) {
            // Coordonnées du point V(i,j)
            double x = x_min + (j + 0.5) * hx;
            double y = y_min + i * hy;

            if ( (x-cx)*(x-cx) + (y-cy)*(y-cy) <= R2 ) {
                _is_solid_v[GetVIndex(i, j)] = true;
                _V(GetVIndex(i, j)) = 0.0;
            }
        }
    }
    
    cout << "Obstacle construit : Cylindre en (" << cx << "," << cy << ") R=" << R << endl;
}