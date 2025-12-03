#ifndef _TIME_SCHEME_H_
#define _TIME_SCHEME_H_

#include "DataFile.h"
#include "Laplacian.h"
#include "MACgrid.h"
#include <string>

class TimeScheme
{
protected:
    DataFile* _df;
    Laplacian* _lap;
    MACgrid* _grid;
    double _t;
    
    // Applique les conditions Dirichlet/Neumann aux frontières
    void ApplyBoundaryConditions();

public:
    TimeScheme(DataFile* data_file, Laplacian* lap, MACgrid* grid);
    virtual ~TimeScheme();
    virtual void Advance() = 0;
    
    double GetTime() const { return _t; }
    void SaveSolution(int n_iteration);
};

class EulerScheme : public TimeScheme
{
public:
    EulerScheme(DataFile* data_file, Laplacian* lap, MACgrid* grid);
    void Advance() override;
};

#endif