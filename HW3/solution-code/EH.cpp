#include "EH.h"
#include <stdio.h>
#include <math.h>
#include <cassert>
#include <unordered_map>
#include <armadillo>

using namespace std;
using namespace arma;

double K_constant = 1.75;
std::map<std::string, double> h_map { {"H1s", -13.6}, {"C2s", -21.4}, {"C2p", -11.4}, };


void Generate_Hmat(mat &OV_mat, vector<AO> &AOs, mat &H_mat){
  int dim = AOs.size();
  assert(OV_mat.n_rows == dim && OV_mat.n_cols == dim);
  assert(H_mat.n_rows == dim && H_mat.n_cols == dim);
  
  for (size_t k = 0; k <dim; k++){
    double h_k = h_map[AOs[k].get_lable()];
    for (size_t j = 0; j < k; j++){
      double h_j = h_map[AOs[j].get_lable()];
      double H_uv = K_constant * (h_j +h_k) * OV_mat(k, j) /2.0;
      H_mat(k,j) = H_uv;
      H_mat(j,k) = H_uv;
    }
      H_mat(k,k) = h_k;
  }
}

double Solve_EH(arma::mat &OV_mat, arma::mat &H_mat, arma::mat &C_mat, arma::vec &energy_vec, int num_ele){
  //Calculate X_mat
  mat U;
  vec S_eigenvalue;
  arma::eig_sym(S_eigenvalue, U, OV_mat);
  mat S_invsqrt = arma::inv( arma::diagmat( arma::sqrt(S_eigenvalue)) );
  // mat X_mat = U * S_invsqrt;
  mat X_mat = U * S_invsqrt * U.t();
  X_mat.print("X_mat");

  mat H_new = X_mat.t() * H_mat * X_mat;
  // H_new.print("H_new");
  
  arma::mat V;
  arma::eig_sym(energy_vec, V, H_new);
  C_mat = X_mat * V;

  if(num_ele % 2 == 1){
    printf("Warn:: this job is unrestricted");
  }
  double Energy = 2* arma::accu(energy_vec.subvec(0, num_ele/2 - 1));
  return Energy;
}
