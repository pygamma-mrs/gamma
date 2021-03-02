/* MultiSOp.i */

%{
#include "MultiSys/MultiSOp.h"
%}

%feature("autodoc", "1" );

gen_op Fx(const multi_sys &msys);
gen_op Fy(const multi_sys &msys);
gen_op Fz(const multi_sys &msys);
gen_op Fe(const multi_sys &msys);
gen_op Fp(const multi_sys &msys);
gen_op Fm(const multi_sys &msys);
 
gen_op Fx(const multi_sys& msys, const std::string& iso); 
gen_op Fy(const multi_sys& msys, const std::string& iso); 
gen_op Fz(const multi_sys& msys, const std::string& iso); 
gen_op Fe(const multi_sys& msys, const std::string& iso); 
gen_op Fm(const multi_sys& msys, const std::string& iso); 
gen_op Fp(const multi_sys& msys, const std::string& iso);

gen_op Rz(const multi_sys& msys, double beta, int icomp=-1);

