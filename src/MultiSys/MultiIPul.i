/* MultiIPul.i */

%{
#include "MultiSys/MultiIPul.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

gen_op Iypuls(const multi_sys& msys, gen_op& sigma, double beta, int icomp=-1);
gen_op Iypuls(const multi_sys& msys, gen_op& sigma, 
                                         int nspin, double beta, int icomp=-1);
gen_op Iypuls(const multi_sys& msys, gen_op& sigma, 
                            const std::string& iso, double beta, int icomp=-1);

gen_op Ixpuls_U(const multi_sys& mys, int spin, double beta, int icomp=-1);

gen_op Ixpuls_U(const multi_sys& mys, double beta, int icomp=-1);
 
gen_op Iypuls_U(const multi_sys& mys, int spin, double beta, int icomp=-1);
 
gen_op Iypuls_U(const multi_sys& mys, const std::string& iso, double beta, int icomp=-1);
 
gen_op Iypuls_U(const multi_sys& mys, double beta, int icomp=-1);
 
gen_op Ixypuls_U(const multi_sys& msys, double phi, double beta, int icomp=-1);
 
