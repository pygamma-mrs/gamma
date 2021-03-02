/* PulseI.h */

%{
#include "HSLib/PulseI.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

gen_op Ixpuls(const spin_sys& sys, const gen_op& sigma, int spin, double beta);
gen_op Ixpuls(const spin_sys& sys, const gen_op& sigma, const std::string& iso,
                                                                  double beta);

gen_op Ixpuls(const spin_sys& sys, const gen_op& sigma, double beta);
//gen_op Ixpuls(const spin_sys& sys, const gen_op& sigma, 
//                                       const flagvec& flags, double beta);

gen_op Ixpuls_sp(const spin_sys& sys, const gen_op& sigma, double beta);

gen_op Iypuls(const spin_sys& sys, const gen_op& sigma, int spin, double beta);
gen_op Iypuls(const spin_sys& sys, const gen_op& sigma, const std::string& iso,
                                                                  double beta);

gen_op Iypuls(const spin_sys& sys, const gen_op& sigma, double beta);

//gen_op Iypuls(const spin_sys& sys, const gen_op& sigma, 
//                                       const flagvec& flags, double beta);

gen_op Iypuls_sp(const spin_sys& sys, const gen_op& sigma, double beta);

gen_op Ixypuls(const spin_sys& sys, const gen_op& sigma,
                                            int spin, double phi, double beta);

gen_op Ixypuls(const spin_sys& sys,const gen_op& sigma,
                                   const std::string& iso, double phi, double beta);

gen_op Ixypuls(const spin_sys& sys, const gen_op& sigma,
                                                      double phi, double beta);

//gen_op Ixypuls(const spin_sys& sys, const gen_op& sigma, 
//                           const flagvec& flags, double phi, double beta);

gen_op Ixypuls_sp(const spin_sys& sys, const gen_op& sigma,
                                                      double phi, double beta);

gen_op Ixpuls_U(const spin_sys& sys, int spin, double beta);
gen_op Ixpuls_U(const spin_sys& sys, const std::string& iso, double beta);
gen_op Ixpuls_U(const spin_sys& sys, double beta);

//gen_op Ixpuls_U(const spin_sys& sys, const flagvec& flags, double beta);

gen_op Ixpuls_sp_U(const spin_sys& sys, double beta);
 
gen_op Iypuls_U(const spin_sys& sys, int spin, double beta);
gen_op Iypuls_U(const spin_sys& sys, const std::string& iso, double beta);
gen_op Iypuls_U(const spin_sys& sys, double beta);
//gen_op Iypuls_U(const spin_sys& sys, const flagvec& flags, double beta);
gen_op Iypuls_sp_U(const spin_sys& sys, double beta);
 
gen_op Ixypuls_U(const spin_sys& sys, int spin, double phi, double beta);
gen_op Ixypuls_U(const spin_sys& sys, const std::string& I,double phi, double beta);
gen_op Ixypuls_U(const spin_sys& sys, double phi, double beta);
//gen_op Ixypuls_U(const spin_sys& sys, const flagvec& flags,
//                                                      double phi, double beta);

gen_op Ixypuls_U_sp(const spin_sys& sys, double phi, double beta);

