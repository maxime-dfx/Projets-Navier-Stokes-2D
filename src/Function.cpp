#include "Function.h"
#include <cmath>

Function::Function(DataFile* data_file) : _df(data_file) {}

// --- Fonction Helper (Privée) ---
// Calcule la vitesse induite par un vortex Gaussien à la position (x,y)
// x0, y0 : Centre du vortex
// R : Rayon caractéristique
// A : Amplitude (Intensité)
// is_u : Si true renvoie U (dPsi/dy), sinon renvoie V (-dPsi/dx)
double Function::ComputeSingleVortex(double x, double y, double x0, double y0, double R, double A, bool is_u) {
    double r2 = (x - x0)*(x - x0) + (y - y0)*(y - y0);
    
    // Formule Gaussienne : Psi = A * exp(-r^2 / R^2)
    double common_factor = A * exp(-r2 / (R * R));
    
    // Dérivée du terme exponentiel : -2/R^2
    double derivative_factor = -2.0 / (R * R);

    if (is_u) {
        // u = dPsi/dy = A * exp(...) * (-2*(y-y0)/R^2)
        return common_factor * derivative_factor * (y - y0);
    } else {
        // v = -dPsi/dx = - [ A * exp(...) * (-2*(x-x0)/R^2) ]
        // Le signe moins de la formule s'annule avec le signe moins de la dérivée
        return - (common_factor * derivative_factor * (x - x0));
    }
}

// --- Conditions Initiales U ---
double Function::InitialConditionU(double x, double y) {
    double Lx = _df->Get_xmax() - _df->Get_xmin();
    double Ly = _df->Get_ymax() - _df->Get_ymin();

    // Vortex 1 : Gros, en bas à gauche
    double u1 = ComputeSingleVortex(x, y, 
                                    0.35 * Lx, 0.40 * Ly, // Position
                                    0.15,                 // Rayon
                                    0,                  // Amplitude
                                    true);                // On veut U

    // Vortex 2 : Petit, en haut à droite (Asymétrique)
    double u2 = ComputeSingleVortex(x, y, 
                                    0.65 * Lx, 0.65 * Ly, 
                                    0.10, 
                                    0, 
                                    true);

    return u1 + u2; // Superposition linéaire
}

// --- Conditions Initiales V ---
double Function::InitialConditionV(double x, double y) {
    double Lx = _df->Get_xmax() - _df->Get_xmin();
    double Ly = _df->Get_ymax() - _df->Get_ymin();

    
    // Vortex 1
    double v1 = ComputeSingleVortex(x, y, 
                                    0.35 * Lx, 0.40 * Ly, 
                                    0.15, 
                                    0, 
                                    false); // On veut V

    // Vortex 2
    double v2 = ComputeSingleVortex(x, y, 
                                    0.65 * Lx, 0.65 * Ly, 
                                    0.10, 
                                    0, 
                                    false);

    return v1 + v2;
}