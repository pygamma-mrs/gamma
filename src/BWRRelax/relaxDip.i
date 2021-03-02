/* relax_Dip.i */

%{
#include "BWRRelax/relaxDip.h"
%}


void RDD(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
                      double* taus, double chi, int type=0, int level=4);

super_op RDD(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

super_op RDD(const spin_system& sys, gen_op& Ho, double tau,
				 matrix& dist, int type=0, int level=4);

super_op RDD_Jgen(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

void RDDrf(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         double Wrflab, double* taus, double chi, int type=0, int level=4);

super_op RDDrf(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type=0, int level=4);

void RDDds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
                          double* taus, double chi, int type=0, int level=4);

super_op RDDds(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

void RDDrfds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
                 double Wrflab, double* taus, double chi, int type=0, int level=4);

super_op RDDrfds(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type=0, int level=4);

matrix xiD(const sys_dynamic& dsys, double cutoff=0.0);
double xiD(const sys_dynamic& dsys, int i, int j);
matrix xiD(const spin_system& sys, matrix& dist, int angs=0, double cutoff=0.0);
double xiD(double gi, double gj, double rij, int angs=0);

matrix DCC(const sys_dynamic& dsys);
double DCC(const sys_dynamic& dsys, int i, int j);
double DCC(double gi, double gj, double rij, int angs=0);

row_vector R1_DD(const sys_dynamic& sys);
double     R1_DD(const sys_dynamic& sys, int i);
double     R1_DD(const sys_dynamic& sys, int i, int j);

double R1_DD_max(const sys_dynamic& sys);
double R1_DD_max(const sys_dynamic& sys, int i);
double R1_DD_max(const sys_dynamic& sys, const std::string& Iso);

row_vector T1_DD(const sys_dynamic& sys);
double     T1_DD(const sys_dynamic& sys, int i);
double     T1_DD(const sys_dynamic& sys, int i, int j);

double T1_DD_max(const sys_dynamic& sys);
double T1_DD_max(const sys_dynamic& sys, int i);
double T1_DD_max(const sys_dynamic& sys, const std::string& Iso);

row_vector R2_DD(const sys_dynamic& sys);
double     R2_DD(const sys_dynamic& sys, int i);
double     R2_DD(const sys_dynamic& sys, int i, int j);

double R2_DD_max(const sys_dynamic& sys);
double R2_DD_max(const sys_dynamic& sys, int i); 
double R2_DD_max(const sys_dynamic& sys, const std::string& Iso); 

row_vector T2_DD(const sys_dynamic& sys);
double     T2_DD(const sys_dynamic& sys, int i);
double     T2_DD(const sys_dynamic& sys, int i, int j);

double T2_DD_max(const sys_dynamic& sys, int i); 
double T2_DD_max(const sys_dynamic& sys, const std::string& Iso); 
double T2_DD_max(const sys_dynamic& sys);

row_vector LWhh_DD(const sys_dynamic& sys);
double LWhh_DD(const sys_dynamic& sys, int i);
double LWhh_DD(const sys_dynamic& sys, int i, int j);

double LWhh_DD_max(const sys_dynamic& sys);
double LWhh_DD_max(const sys_dynamic& sys, int i);
double LWhh_DD_max(const sys_dynamic& sys, const std::string& Iso);

double NOE(const sys_dynamic& sys, int i, int j, double eta=0);

row_vector R2_DDMQT(const sys_dynamic& sys, int MQC);
double     R2_DDMQT(const sys_dynamic& sys, int MQC, int i);
double     R2_DDMQT(const sys_dynamic& sys, int MQC, int i, int j);
