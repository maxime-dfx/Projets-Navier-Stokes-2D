#include "Function.h"
#include <cmath>

Function::Function(DataFile* data_file) : _df(data_file) {}

// --- Conditions Initiales U ---
double Function::InitialConditionU(double x, double y) {
    // Initialisation uniforme à la vitesse du vent entrant (1.0)
    // Cela évite le "bang" de démarrage (choc 1.0 -> 0.0)
    if (_df->Get_BC_Left() == "Dirichlet") {
        return _df->Get_BC_Left_dir();
    }
    return 0.0;
}

// --- Conditions Initiales V ---
double Function::InitialConditionV(double x, double y) {
    // Pas de vitesse verticale initiale
    return 0.0;
}