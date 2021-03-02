/* relax_Quad.i */

%{
#include "BWRRelax/relaxQuad.h"
%}

void RQQ(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
                     double* taus, double chi, int type=0, int level=4);

super_op RQQ(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

void RQQrf(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int type=0, int level=4);

super_op RQQrf(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type=0, int level=4);

void RQQds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int type=0, int level=4);

super_op RQQds(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

RQQrfds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int type=0, int level=4);

super_op RQQrfds(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type=0, int level=4);

row_vector R1_QQ(const sys_dynamic& sys);

double R1_QQ(const sys_dynamic& sys, int i);

row_vector T1_QQ(const sys_dynamic& sys);

double R1_QQ_max(const sys_dynamic& sys, int i);

double R1_QQ_max(const sys_dynamic& sys, const std::string& Iso);

double R1_QQ_max(const sys_dynamic& sys);

row_vector T1_QQ(const sys_dynamic& sys);

double T1_QQ(const sys_dynamic& sys, int i);

double T1_QQ_max(const sys_dynamic& sys, int i);
 
double T1_QQ_max(const sys_dynamic& sys, const std::string& Iso);

double T1_QQ_max(const sys_dynamic& sys);
 
row_vector R2_QQ(const sys_dynamic& sys);

double R2_QQ(const sys_dynamic& sys, int i);

double R2_QQ(const sys_dynamic& sys, int I);
 
double R2_QQ_max(const sys_dynamic& sys, int i);
 
double R2_QQ_max(const sys_dynamic& sys, const std::string& Iso);

row_vector T2_QQ(const sys_dynamic& sys);

double T2_QQ(const sys_dynamic& sys, int i);

double T2_QQ_max(const sys_dynamic& sys, int i);
 
double T2_QQ_max(const sys_dynamic& sys, const std::string& Iso);

double T2_QQ_max(const sys_dynamic& sys);
 
row_vector LWhh_QQ(const sys_dynamic& sys);

double LWhh_QQ(const sys_dynamic& sys, int i);

double LWhh_QQ_max(const sys_dynamic& sys, int i);

double LWhh_QQ_max(const sys_dynamic& sys, const std::string& Iso);
 
double LWhh_QQ_max(const sys_dynamic& sys);

matrix xiQ(const sys_dynamic& sys);

double xiQ(const sys_dynamic& sys, int i);
