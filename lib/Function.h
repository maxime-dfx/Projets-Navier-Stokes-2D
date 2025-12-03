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

    // Conditions Initiales (t=0)
    double InitialConditionU(double x, double y);
    double InitialConditionV(double x, double y);

    // --- LOGIQUE CONDITIONNELLE ---
    bool IsDirichletLeft()   const { return _df->Get_BC_Left() == "Dirichlet"; } 
    bool IsDirichletRight()  const { return _df->Get_BC_Right() == "Dirichlet"; }
    bool IsDirichletBottom() const { return _df->Get_BC_Bottom() == "Dirichlet"; }
    bool IsDirichletTop()    const { return _df->Get_BC_Top() == "Dirichlet"; } 

    // --- VITESSES NORMALES (Traversée du mur) ---
    // Pour des murs solides, c'est toujours 0.0.
    double GetLeftU_Normal(double y)   const { return _df->Get_BC_Left_dir(); }    
    double GetRightU_Normal(double y)  const { return 0.0; } // Non utilisé si Neumann à droite
    double GetBottomV_Normal(double x) const { return 0.0; } // Mur étanche
    double GetTopV_Normal(double x)    const { return 0.0; } // Mur étanche

    // --- VITESSES TANGENTIELLES (Glissement du mur) ---
    // C'est ici qu'on utilise les valeurs du fichier input.
    // Ex: Lid-Driven Cavity -> Le mur du haut glisse.
    
    // Mur Gauche/Droit : La vitesse tangentielle est V (Verticale)
    double GetLeftV_Tangent(double y)   const { return 0.0; }
    double GetRightV_Tangent(double y)  const { return 0.0; } 

    // Mur Bas/Haut : La vitesse tangentielle est U (Horizontale)
    double GetBottomU_Tangent(double x) const { return _df->Get_BC_Bottom_dir(); } // Correction: Bottom
    double GetTopU_Tangent(double x)    const { return _df->Get_BC_Top_dir(); }    // Correction: Top
};

#endif