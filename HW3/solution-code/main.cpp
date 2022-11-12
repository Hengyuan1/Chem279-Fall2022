#include <stdlib.h>
#include <stdexcept>
#include <stdio.h>
#include <armadillo>
#include <vector>
#include "AO.h"
#include "EH.h"

using namespace std;

int main(int argc, char* argv[])
{
  arma::mat H_STO3G, C_STO3G;
  H_STO3G.load("H_STO3G.txt");
  C_STO3G.load("C_STO3G.txt");
  // H_STO3G.print("H_STO3G");
  // C_STO3G.print("C_STO3G");

  vector<AO> MoleculeAOs;
  int num_ele;
  
  // // To debug 
  // arma::mat O_STO3G;
  // O_STO3G.load("H_STO3G.txt");
  // GenerateAOs(MoleculeAOs, "HO.txt", H_STO3G, O_STO3G);
  // AO H1s = MoleculeAOs[0];
  // AO Cpz = MoleculeAOs[4];
  // H1s.printinfo();
  // Cpz.printinfo();
  // arma::uvec lmna = H1s.get_lmn();
  // arma::uvec lmnb = Cpz.get_lmn();
  // arma::vec Ra = H1s.get_R0();
  // arma::vec Rb = Cpz.get_R0();
  // double alphaa =3.4253;
  // double alphab =5.0332;
  // double temp1 = Overlap_3d(Ra, Rb, alphaa, alphab, lmna, lmnb);
  // printf("Overlap of H1s and C2pz %1.10f\n", temp1);
  // double temp2 = Eval_Ov_AOs(H1s, Cpz);
  // printf("Overlap of H1s and C2pz %1.10f\n", temp2);

  if (argc !=2)
  {
  printf("usage hw3 filename, for example hw3 example.txt\n");
  return EXIT_FAILURE;
  }
  string fname(argv[1]);
  try
  {
    num_ele = GenerateAOs(MoleculeAOs, fname, H_STO3G, C_STO3G);
  }
  catch (invalid_argument &e)
  {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  // PrintAOs(MoleculeAOs);


  int dim = MoleculeAOs.size();
  arma::mat OV_mat(dim, dim);
  Eval_OV_mat(MoleculeAOs, OV_mat);
  OV_mat.print("OV_mat");

  arma::mat H_mat(dim, dim);
  Generate_Hmat(OV_mat, MoleculeAOs, H_mat);
  H_mat.print("H_mat");

  arma::mat C_mat(dim, dim);
  arma::vec energy_vec(dim);
  double Energy = Solve_EH(OV_mat, H_mat, C_mat, energy_vec, num_ele);
  C_mat.print("C_mat");
  // check the MOs are orthonormal
  arma::mat MO_ov = C_mat.t()*OV_mat*C_mat;
  MO_ov.print("MO overlap");
  // energy_vec.print("energy_vec");
  printf("The molecule in file %s has energy %f\n", argv[1], Energy);

  return EXIT_SUCCESS;
}
