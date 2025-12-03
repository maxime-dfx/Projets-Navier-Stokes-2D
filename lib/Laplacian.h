#ifndef _LAPLACIAN_H
#define _LAPLACIAN_H

#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Function.h"
#include "MACgrid.h"

class Laplacian
{
private:
	Function* _fct;
	DataFile* _df;
	MACgrid* _grid;
	
    Eigen::SparseMatrix<double> _H; 
    
    // Changement de solveur : SimplicialLLT est plus rapide et stable 
    // pour les matrices symétriques définies positives (Laplacien)
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > _solver;

public:
	Laplacian(Function* function, DataFile* data_file, MACgrid*& grid);

	void BuildMatrix();
	Eigen::VectorXd ComputeDivergence(const Eigen::VectorXd& U, const Eigen::VectorXd& V);
	void Solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& p_sol);
	void ComputeGradient(const Eigen::VectorXd& p, Eigen::VectorXd& gradPx, Eigen::VectorXd& gradPy);
    const Eigen::MatrixXd& Get_H() const { return _H; }
};

#endif