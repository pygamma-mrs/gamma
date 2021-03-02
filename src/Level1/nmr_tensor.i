/* nmr_tensor.i */
// Swig interface file.

%{
#include "Level1/nmr_tensor.h"
%}

%feature("autodoc", "1" );

spin_T T_D(const spin_sys &sys, int spin1, int spin2);

spin_T T_D(const spin_sys &sys, spin_op &Im1, spin_op &Iz1, spin_op &Ip1,
				spin_op &Im2, spin_op &Iz2, spin_op &Ip2);

spin_op T_D(const spin_sys &sys, int spin1, int spin2, int m);

spin_T T_CSA(const spin_sys &sys, int spin);

spin_T T_CS2(const spin_sys &sys, int spin);

spin_T T_CS2(const spin_sys &sys, int spin, coord &B);

spin_op T_CS2(const spin_sys &sys, int spin, coord &B, int l, int m);

spin_T T_CS(const spin_sys &sys, int spin);

spin_op T_CS(const spin_sys &sys, int spin, int m);

spin_T T_RF(const spin_sys &sys, int spin);

spin_op T_RF(const spin_sys &sys, int spin, int l, int m);

spin_T T_Q(const spin_sys &sys, int spin);

spin_op T_Q(const spin_sys &sys, int spin, int l, int m);
