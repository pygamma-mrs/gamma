/* relaxNMR.i */

%{
#include "BWRRelax/relaxNMR.h"
%}

%feature("autodoc", "1" );

void RlxNMRerror(int eidx, int noret=0);
volatile void RlxNMRfatal(int eidx);

void R_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12);

double R_4(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	         int rank, int a, int b, int aa, int bb);

void R_3(super_op& LOp, double* w, int rank, gen_op* T1s, gen_op* T2s,
				          matrix& J12, double cutoff=1.e-2);

void R_2(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12);


double R_2(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	         int rank, int a, int b, int aa, int bb);

double Rodiag_2(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	                      int rank, int a, int b);

double Rdiag_2(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	                      int rank, int a, int aa);

void R_0(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, const complex& J12);

double R_0(int hs, gen_op* T1s, gen_op* T2s,
		 	         int rank, int a, int b, int aa, int bb);


void R_4s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12);

void R_3s(super_op& LOp, double* w, int rank, gen_op* T1s, gen_op* T2s, matrix& J12);

void R_2s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12);

super_op R_AC_0(spin_T &T);

void R_AC_0(spin_T &T, super_op &LOp1, gen_op &Op, double xisq = 1);

void R_AC_0(gen_op *Ts, super_op &LOp, int rank=2, double xisq=1);

void R_CC_0(spin_T &T1, spin_T &T2, super_op &LOp1,
		                     gen_op &Op, double xisq = 1);

void R_CC_0(gen_op *T1s, gen_op* T2s, super_op &LOp1,
		             int rank=2, double xisq = 1);

void R_CC_0_trans(gen_op *T1s, gen_op* T2s, super_op &LOp1,
		                   int rank=2, double xisq = 1);

void R_AC_1(spin_T &T, super_op &LOp1, gen_op &Op,
		             	double J0, double J1, double J2);

void R_AC_1(gen_op *Ts, super_op &LOp1, int rank,
		             	double J0, double J1, double J2);

void R_CC_1(spin_T &T1, spin_T &T2, super_op &LOp1,
		        gen_op &Op, double J0, double J1, double J2);

void R_CC_1(gen_op *T1s, gen_op* T2s, super_op &LOp1,
		        int rank, double J0, double J1, double J2);

void Rmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
               double w0, double w1, double w2, int l, int level=4, int autoc=0, int het=0);

void Rmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w,
                      int hs, double tau, double xi1xi2, double w0,
                            double w1, double w2, int l, int level=4, int autoc=0);

void Rmu1mu2(super_op& LOp, const spin_system& sys, gen_op& Ho, double* w,
         double* xi1s, double n1, double* xi2s, double n2, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int l, int type=0, int level=4);

void Rijkl(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

void Rijkl(super_op& LOp, const spin_system& sys, gen_op& Ho, double* w,
                  matrix& xi1s, matrix& xi2s, spin_T* T1, spin_T* T2,
                                    double tau, int type=0, int level=4);

void Rij(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

void Rij(super_op& LOp, const spin_system& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, spin_T* T1, spin_T* T2,
                       double tau, int l, int type=0, int level=4);

void Rijk(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

void Rkij(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

void Rmumuds(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
               double w0, double w1, double w2, int level=4, int autoc=0, int het=0);

void Rijklds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

void Rijds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

void Rijkds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

void Rkijds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

