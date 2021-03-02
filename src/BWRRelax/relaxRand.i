/* relaxRand.i */

%{
#include "BWRRelax/relaxRand.h"
%}


void RRRx(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                                     double tau, int type=0, int level=4);

super_op RRRx(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

void RRR(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int type=0, int level=4);

super_op RRR(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

void Rij_rdm(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

void Rmumu_rdm(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
               double w0, double w1, double w2, int level=4, int autoc=0);

row_vector R1_RR(const sys_dynamic& sys);

double R1_RR(const sys_dynamic& sys, int i);

double R1_RR_max(const sys_dynamic& sys);

row_vector T1_RR(const sys_dynamic& sys);

double T1_RR(const sys_dynamic& sys, int i);

double T1_RR_max(const sys_dynamic& sys);

row_vector R2_RR(const sys_dynamic& sys);

double R2_RR(const sys_dynamic& sys, int i);

double R2_RR_max(const sys_dynamic& sys);

row_vector T2_RR(const sys_dynamic& sys);

double T2_RR(const sys_dynamic& sys, int i);

double T2_RR_max(const sys_dynamic& sys);

row_vector LWhh_RR(const sys_dynamic& sys);

double LWhh_RR(const sys_dynamic& sys, int i);

double LWhh_RR_max(const sys_dynamic& sys, int i);
 
double LWhh_RR_max(const sys_dynamic& sys, const std::string& Iso);
 
double LWhh_RR_max(const sys_dynamic& sys);

matrix xiRDM(const sys_dynamic& dsys);

double xiRDM(const sys_dynamic& dsys, int i);

