/* MultiHSLib.i */

%{
#include "MultiSys/MultiHSLib.h"
%}

%feature("autodoc", "1" );

gen_op sigma_eq(const multi_sys &msys);
