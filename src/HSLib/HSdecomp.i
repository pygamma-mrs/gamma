/* HSdecomp.i */

%{
#include "HSLib/HSdecomp.h"
%}

%feature("autodoc", "1" );

void Prod_base_dec(const spin_sys &sys, const gen_op &Op, double thres=BD_SMALL);
