#ifndef AO_H
#define AO_H

#include <iostream>
#include <armadillo>

double Overlap_onedim(double xa, double xb, double alphaa, double alphab, int la, int lb);


class AO
{
    private:
        arma::vec R0;
        arma::uvec lmn;
        arma::vec alpha;
        arma::vec d_coe;
        int len;
        std::string lable;
    public:
        AO(arma::vec &R0_input, arma::vec &alpha_input, arma::vec &d_input, arma::uvec &lmn_input, std::string lable_input);
        ~AO(){}
        void printinfo();
        arma::uvec get_lmn(){ return lmn;}
        arma::vec get_alpha(){ return alpha;}
        arma::vec get_d_coe(){ return d_coe;}
        arma::vec get_R0(){ return R0;}
        int get_len(){ return len;}
        std::string get_lable(){ return lable;}
};
int GenerateAOs(std::vector<AO> &AOs, std::string &fname, arma::mat &H_basis, arma::mat &C_basis);
void PrintAOs(std::vector<AO> &AOs);
double Overlap_3d(arma::vec &Ra, arma::vec &Rb, double alphaa, double alphab, arma::uvec &lmna, arma::uvec &lmnb);

double Eval_Ov_AOs(AO& sh1, AO& sh2);

void Eval_OV_mat(std::vector<AO> &MoleculeAOs, arma::mat &OV_mat);

#endif // AO_H