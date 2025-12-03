#ifndef _MACGRID_H_
#define _MACGRID_H_

#include <Eigen/Dense>
#include "Function.h"
#include "DataFile.h"

class MACgrid
{
private:
    Function* _fct; 
    DataFile* _df;
    
    int _Nx; 
    int _Ny;

    Eigen::VectorXd _p; // Pression aux centres
    Eigen::VectorXd _U; // Vitesse X aux faces verticales
    Eigen::VectorXd _V; // Vitesse Y aux faces horizontales
    std::vector<bool> _is_solid_u; // Masque pour la grille U
    std::vector<bool> _is_solid_v; // Masque pour la grille V

public:
    MACgrid(Function* function, DataFile* data_file);
    
    Function* GetFunction() const { return _fct; }
    
    const Eigen::VectorXd& GetP() const { return _p; }
    const Eigen::VectorXd& GetU() const { return _U; }
    const Eigen::VectorXd& GetV() const { return _V; }

    void SetP(const Eigen::VectorXd& p) { _p = p; }
    void SetU(const Eigen::VectorXd& u) { _U = u; }
    void SetV(const Eigen::VectorXd& v) { _V = v; }
    
    // Indexation
    int GetPIndex(int i, int j) const;
    int GetUIndex(int i, int j) const;
    int GetVIndex(int i, int j) const;

    // Coordonnées
    Eigen::VectorXd GetPcoord(int i, int j) const;
    Eigen::VectorXd GetUcoord(int i, int j) const;
    Eigen::VectorXd GetVcoord(int i, int j) const;

    void BuildObstacles(); // Fonction pour créer les formes
    bool IsSolidU(int i, int j) const { return _is_solid_u[GetUIndex(i,j)]; }
    bool IsSolidV(int i, int j) const { return _is_solid_v[GetVIndex(i,j)]; }
};

#endif