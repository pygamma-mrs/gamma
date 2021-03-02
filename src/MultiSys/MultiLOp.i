/* MultiLOp.i */

%{
#include "MultiSys/MultiLOp.h"
%}

%feature("autodoc", "1" );

super_op Hsuper(const multi_sys& msys, const gen_op& Heff);
super_op Lo(const multi_sys &msys);
super_op U_LS(gen_op& H);
super_op Uinv_LS(gen_op& H);
super_op Op_Ebase(super_op& L, gen_op& H);
