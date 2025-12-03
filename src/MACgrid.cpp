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
    
    _p.setZero();
    _U.setZero();
    _V.setZero(); // IMPORTANT: Initialise à 0 pour être sûr

    // Conditions initiales
    for (int i = 0; i < _Ny; ++i) {
        for (int j = 0; j <= _Nx; ++j) {
            int k = GetUIndex(i, j);
            VectorXd coord = GetUcoord(i, j);
            _U(k) = _fct->InitialConditionU(coord(0), coord(1));
        }
    }

    for (int i = 0; i <= _Ny; ++i) {
        for (int j = 0; j < _Nx; ++j) {
            int k = GetVIndex(i, j);
            VectorXd coord = GetVcoord(i, j);
            _V(k) = _fct->InitialConditionV(coord(0), coord(1));
        }
    }
}

int MACgrid::GetPIndex(int i, int j) const {
    return j + i * _Nx; 
}

int MACgrid::GetUIndex(int i, int j) const {
    // Largeur ligne U = Nx + 1
    return j + i * (_Nx + 1); 
}

int MACgrid::GetVIndex(int i, int j) const {
    // Largeur ligne V = Nx
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