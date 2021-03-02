/* relaxRF.i */

%{
#include "BWRRelax/relaxRF.h"
%}


void Rrf_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12);

complex Rrf_4(int hs, gen_op* T1s, gen_op* T2s, matrix* J12,
		 	         int rank, int a, int b, int aa, int bb);

void Rrf_3(super_op& LOp, double* w, int rank,
                   gen_op* T1s, gen_op* T2s, matrix* J12, double cutoff=1.e-2);

void Rrf_2(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12);


double Rrf_2(int hs, gen_op* T1s, gen_op* T2s, matrix* J12,
		 	         int rank, int a, int b, int aa, int bb);

void Rrf_0(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, const complex& J12);

double Rrf_0(int hs, gen_op* T1s, gen_op* T2s,
		 	         int rank, int a, int b, int aa, int bb);


void Rrf_4s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12);

void Rrf_3s(super_op& LOp, double* w, int rank, gen_op* T1s, gen_op* T2s, matrix* J12);

void Rrf_2s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12);

void Rrfmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, matrix* J12,
                 double* J, double* w, int rank=2, int level=4, int autoc=0, int het=0);

void Rrfijkl(super_op &LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

void Rrfij(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

void Rrfijk(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

void Rrfkij(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

void Rrfijklds(super_op &LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

void Rrfijds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

void Rrfijkds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level);

void Rrfkijds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level);
 
gen_op sigma_ss(spin_system& sys, super_op& L, super_op& R);

gen_op sigma_ss_it(spin_system& sys, super_op& L, super_op& Heff, super_op& R);

