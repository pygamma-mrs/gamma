/* PulseS.i */

%{
#include "HSLib/PulseS.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

void PulSerror(int eidx, int noret=0);

void volatile PulSfatality(int eidx);

gen_op Sxpuls(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
       const std::string& iso, double freq=0.0, double tp=1.e-5, double theta=90.0);

gen_op SxpulsB(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
      const std::string& iso, double freq=0.0, double tp=1.e-5, double gamB1=2.5e4);

gen_op Sypuls(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
       const std::string& iso, double freq=0.0, double tp=1.e-5, double theta=90.0);

gen_op SypulsB(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
      const std::string& iso, double freq=0.0, double tp=1.e-5, double gamB1=2.5e4);

gen_op Sxypuls(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
                const std::string& iso,double freq=0.0,
                            double tp=1.e-5,double theta=90.0, double phi=0.0);

gen_op SxypulsB(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
                const std::string& iso,double freq=0.0,
                            double tp=1.e-5,double gamB1=2.5e4,double phi=0.0);

gen_op Sxpuls_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
                          double freq=0.0, double tp=1.e-5, double theta=90.0);

gen_op SxpulsB_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
                         double freq=0.0, double tp=1.e-5, double gamB1=2.5e4);

gen_op Sypuls_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
                          double freq=0.0, double tp=1.e-5, double theta=90.0);

gen_op SypulsB_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
                         double freq=0.0, double tp=1.e-5, double gamB1=2.5e4);

gen_op Sxypuls_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
          double freq=0.0, double tp=1.e-5, double theta=90.0, double phi=0.0);

gen_op SxypulsB_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
         double freq=0.0, double tp=1.e-5, double gamB1=2.5e4, double phi=0.0);

gen_op Sxypuls(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
     const std::string& iso1, double freq1, const std::string& iso2, double freq2,
                           double tp=1.e-5, double theta=90.0, double phi=0.0);

gen_op SxypulsB(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
     const std::string& iso1, double freq1, const std::string& iso2, double freq2,
                          double tp=1.e-5, double gamB1=2.5e4, double phi=0.0);

gen_op Sxypuls_U(const spin_sys& sys, const gen_op& H,
       const std::string& iso1, double freq1, const std::string& iso2, double freq2,
                           double tp=1.e-5, double theta=90.0, double phi=0.0);

gen_op SxypulsB_U(const spin_sys& sys, const gen_op& H,
	const std::string& iso1, double freq1, const std::string& iso2, double freq2,
                          double tp=1.e-5, double gamB1=2.5e4, double phi=0.0);

gen_op Spul_axis(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
            const std::string& iso, double freq, double tp, double fact, char axis);

gen_op Spul_U_axis(const spin_sys& sys, const gen_op& H, const std::string& iso,
   			       double freq, double tp, double fact, char axis);

gen_op Spul_plane(const spin_sys& sys, const gen_op& sigma, const gen_op& H, 
           const std::string& iso, double freq, double tp, double fact, double phi);

gen_op Spul_plane(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
            const std::string& iso1, double freq1, const std::string& iso2,
                             double freq2, double tp, double fact, double phi);

gen_op Spul_U_plane(const spin_sys& sys, const gen_op& H, const std::string& iso,
   			      double freq, double tp, double fact, double phi);

gen_op Spul_U_plane(const spin_sys& sys, const gen_op& H, const std::string& iso1,
	 double freq1, const std::string& iso2, double freq2, 
                                           double tp, double fact, double phi);


