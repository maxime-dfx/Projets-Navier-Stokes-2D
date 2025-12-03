#ifndef _FUNCTION_H_
#define _FUNCTION_H_

#include "DataFile.h" 
#include <string>
#include <cmath>

class Function {
private:
    DataFile* _df;

public:
    Function(DataFile* data_file); 

    // Conditions Initiales (Vitesse à t=0)
    double InitialConditionU(double x, double y);
    double InitialConditionV(double x, double y);

    // --- LOGIQUE CONDITIONNELLE ---
    // Vérifie si le bord est un mur solide (Dirichlet) ou une sortie (Neumann)
    bool IsDirichletLeft()   const { return _df->Get_BC_Left() == "Dirichlet"; } 
    bool IsDirichletRight()  const { return _df->Get_BC_Right() == "Dirichlet"; }
    bool IsDirichletBottom() const { return _df->Get_BC_Bottom() == "Dirichlet"; }
    bool IsDirichletTop()    const { return _df->Get_BC_Top() == "Dirichlet"; } 

    // --- VALEURS IMPOSÉES AUX MURS (Dirichlet) ---
    // Pour l'instant : murs fixes (vitesse = 0)
    // Si vous voulez une "Cavité Entrainée" (Lid-Driven), mettez GetTopU = 1.0 (attention à TimeScheme)
    double GetLeftU(double y)   const { return 0.0; }
    double GetRightU(double y)  const { return 0.0; }
    double GetBottomV(double x) const { return 0.0; }
    double GetTopV(double x)    const { return 0.0; } 
};

#endif