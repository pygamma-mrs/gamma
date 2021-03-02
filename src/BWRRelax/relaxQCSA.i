/* relax_QCSA.i */

%{
#include "BWRRelax/relaxQCSA.h"
%}


super_op RQCX(const sys_dynamic& sys, gen_op& Ho, int level=4);

void RQCX(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                                double* taus, double chi, int level=4);

super_op RQC(const sys_dynamic& sys, gen_op& Ho, int level=4);

void RQC(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int level=4);

super_op RCQ(const sys_dynamic& sys, gen_op& Ho, int level=4);

void RCQ(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                              double* taus, double chi, int level=4);

super_op RCQrf(const sys_dynamic& sys, gen_op& Heff, double Wrflab, int level=4);

void RCQrf(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int level=4);
