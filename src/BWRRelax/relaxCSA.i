/* relaxCSA.i */

%{
#include "BWRRelax/relaxCSA.h"
%}


void RCC(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w, double* taus, double chi, int type=0, int level=4);

super_op RCC(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

void RCCrf(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w, double Wrflab, double* taus, double chi, int type=0, int level=4);

super_op RCCrf(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type=0, int level=4);

void RCCds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w, double* taus, double chi, int type=0, int level=4);

super_op RCCds(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

void RCCrfds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w, double Wrflab, double* taus, double chi, int type=0, int level=4);

super_op RCCrfds(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type=0, int level=4);

row_vector R1_CC(const sys_dynamic& sys);

double R1_CC(const sys_dynamic& sys, int I);

double R1_CC_max(const sys_dynamic& sys, int i);

double R1_CC_max(const sys_dynamic& sys, const std::string& Iso);

double R1_CC_max(const sys_dynamic& sys);

row_vector T1_CC(const sys_dynamic& sys);

double T1_CC(const sys_dynamic& sys, int I);

double T1_CC_max(const sys_dynamic& sys, int i);

double T1_CC_max(const sys_dynamic& sys, const std::string& Iso);

double T1_CC_max(const sys_dynamic& sys);

row_vector R2_CC(const sys_dynamic& sys);

double R2_CC(const sys_dynamic& sys, int I);

double R2_CC_max(const sys_dynamic& sys, int i);

double R2_CC_max(const sys_dynamic& sys, const std::string& Iso);

double R2_CC_max(const sys_dynamic& sys);

row_vector T2_CC(const sys_dynamic& sys);

double T2_CC(const sys_dynamic& sys, int I);

double T2_CC_max(const sys_dynamic& sys, int i);

double T2_CC_max(const sys_dynamic& sys, const std::string& Iso);

double T2_CC_max(const sys_dynamic& sys);

row_vector LWhh_CC(const sys_dynamic& sys);

double LWhh_CC(const sys_dynamic& sys, int I);

double LWhh_CC_max(const sys_dynamic& sys, int i);

double LWhh_CC_max(const sys_dynamic& sys, const std::string& Iso);

double LWhh_CC_max(const sys_dynamic& sys);

matrix xiCSA(const sys_dynamic& dsys);

double xiCSA(const sys_dynamic& dsys, int i);

matrix xiCSA(const spin_system& sys, double* CSAs);

double xiCSA(const spin_system& sys, int i, double csa);

row_vector CSA(const sys_dynamic& dsys);

double CSA(const sys_dynamic& dsys, int i);
