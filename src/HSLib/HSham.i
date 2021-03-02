// HSham.i
// Swig interface file.

%{
#include "HSLib/HSham.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

// defines a set of functions that convert spin_system to gen_op(s).

gen_op Hcs(const spin_system &sys);
gen_op Hcs_lab(const spin_system &sys);

gen_op HJw(const  spin_system &sys);
gen_op HJ(const   spin_system &sys);
gen_op HJwh(const spin_system &sys);
gen_op HJd(const spin_system& sys, const std::string& iso);

gen_op Ho(const     spin_system &ss);
gen_op How(const    spin_system &ss);
gen_op Ho_lab(const spin_system &ss);

gen_op Hz(const spin_system& sys);
gen_op Hz(const spin_system& sys, const std::string& I);

gen_op H1(const spin_system& sys, const std::string& iso,
                                    double gamB1=2.5e4, double phi=0.0);

gen_op Heff(spin_sys& sys, gen_op &H0, const std::string& iso,
                                  double Wrf=0, double gamB1=2.5e4, double phi=0.0);

gen_op Hg(const     spin_system& sys);
gen_op Hg_lab(const spin_system& sys);
                                                                                
gen_op HAw(const spin_system& sys);

gen_op HQsec(const spin_system& sys, double wQ, int i);

