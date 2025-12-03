#include "Function.h"
#include <cmath>

Function::Function(DataFile* data_file) : _df(data_file) {}

// Taylor-Green "Classique"
// C'est celui qui va créer de la pression car il "cogne" contre les murs au début
double Function::InitialConditionU(double x, double y) {
    double Lx = _df->Get_xmax() - _df->Get_xmin();
    double Ly = _df->Get_ymax() - _df->Get_ymin();
    double X = (x - _df->Get_xmin()) / Lx;
    double Y = (y - _df->Get_ymin()) / Ly;

    const double PI = 3.14159265358979323846;
    return sin(2.0 * PI * X) * cos(2.0 * PI * Y);
}

double Function::InitialConditionV(double x, double y) {
    double Lx = _df->Get_xmax() - _df->Get_xmin();
    double Ly = _df->Get_ymax() - _df->Get_ymin();
    double X = (x - _df->Get_xmin()) / Lx;
    double Y = (y - _df->Get_ymin()) / Ly;

    const double PI = 3.14159265358979323846;
    return -cos(2.0 * PI * X) * sin(2.0 * PI * Y);
}