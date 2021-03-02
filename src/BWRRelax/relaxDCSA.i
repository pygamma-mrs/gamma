/* relaxDCSA.i*/

%{
#include "BWRRelax/relaxDCSA.h"
%}


super_op RDCX(const sys_dynamic& sys, gen_op& Ho, int level=4);

void RDCX(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                                double* taus, double chi, int level=4);

super_op RDC(const sys_dynamic& sys, gen_op& Ho, int level=4);

void RDC(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int level=4);

super_op RCD(const sys_dynamic& sys, gen_op& Ho, int level=4);

void RCD(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                              double* taus, double chi, int level=4);

super_op RCDrf(const sys_dynamic& sys, gen_op& Heff, double Wrflab, int level=4);

void RCDrf(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int level=4);
