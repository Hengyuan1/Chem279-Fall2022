#if !defined EH_H
#define EH_H
#include <armadillo>
#include <cassert>
#include "AO.h"

void Generate_Hmat(arma::mat &OV_mat, std::vector<AO> &AOs, arma::mat &H_mat);

double Solve_EH(arma::mat &OV_mat, arma::mat &H_mat, arma::mat &C_mat, arma::vec &energy_vec, int num_ele);


#endif // EH_H