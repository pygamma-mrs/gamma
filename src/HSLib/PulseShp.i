/* PulseShp.i */

%{
#include "HSLib/PulseShp.h"
%}

%feature("autodoc", "1" );

void PulSherror(int eidx, int noret);

void volatile PulShfatality(int eidx);

gen_op Shxpuls(const spin_sys& sys, row_vector& BLK, gen_op sigma, gen_op& H,
	 const std::string& iso, double freq=0.0, double time=1.0e-5, double theta=90.0);

gen_op Shxpuls_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
		double freq=0.0, double time=1.0e-5, double theta=90.0);

gen_op Shypuls(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
	 const std::string& iso, double freq=0.0, double time=1.0e-5, double theta=90.0);

gen_op Shypuls_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
		double freq=0.0, double time=1.0e-5, double theta=90.0);

gen_op Shxypuls(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
		 const std::string& iso, double freq=0.0, double time=1.0e-5,
					 double theta=90.0, double phi=0.0);

gen_op Shxypuls(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
		 const std::string& iso1, double freq1, const std::string& iso2, double freq2,
			 double time=1.0e-5, double theta=90.0, double phi=0.0);

gen_op Shxypuls_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
	double freq=0.0, double time=1.0e-5, double theta=90.0, double phi=0.0);

gen_op Shxypuls_U(const spin_sys &sys, row_vector &BLK, gen_op &H,
		 const std::string& iso1, double freq1, const std::string& iso2, double freq2,
			 double time=1.0e-5, double theta=90.0, double phi=0.0);

gen_op ShxpulsB(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
	 const std::string& iso, double freq=0.0, double time=1.0e-5, double gamB1=2.5e4);

gen_op ShxpulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
			double freq=0.0, double time=1.0e-5, double gamB1=2.5e4);
//gen_op ShxpulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H, char *iso,
//			double freq=0.0, double time=1.0e-5, double gamB1=2.5e4);

gen_op ShypulsB(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
	 const std::string& iso, double freq=0.0, double time=1.0e-5, double gamB1=2.5e4);

gen_op ShypulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
			double freq=0.0, double time=1.0e-5, double gamB1=2.5e4);

gen_op ShxypulsB(const spin_sys &sys, row_vector &BLK, gen_op sigma,
	 	gen_op &H, const std::string& iso, double freq=0.0,
			 double time=1.0e-5, double gamB1=2.5e4, double phi=0.0);

gen_op ShxypulsB(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
	 const std::string& iso1, double freq1, const std::string& iso2, double freq2,
			 double time=1.0e-5, double gamB1=2.5e4, double phi=0.0);

gen_op ShxypulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
	double freq=0.0, double time=1.0e-5, double gamB1=2.5e4, double phi=0.0);

gen_op ShxypulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H,
	 const std::string& iso1, double freq1, const std::string& iso2, double freq2,
			 double time=1.0e-5, double gamB1=2.5e4, double phi=0.0);

gen_op Shpul_axis(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
		 const std::string& iso, double freq, double time, double fact, char axis);

gen_op Shpul_U_axis(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
			double freq, double time, double fact, char axis);

gen_op Shpul_plane(const spin_sys &sys, row_vector &BLK, gen_op sigma,
	 gen_op &H, const std::string& iso, double freq, double time, double fact, double phi);

gen_op Shpul_plane(const spin_sys &sys, row_vector &BLK, gen_op sigma,
		 gen_op &H, const std::string& iso1, double freq1, const std::string& iso2,
			 double freq2, double time, double fact, double phi);

gen_op Shpul_U_plane(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
			 double freq, double time, double fact, double phi);

//  gen_op Shpul_U_plane(const spin_sys &sys, row_vector &BLK, gen_op &H, char *iso1,
//	double freq1, char *iso2, double freq2, double time, double fact, double phi);

gen_op Shpul_U_plane(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso1,
 	double freq1, const std::string& iso2, double freq2, double time, double fact, double phi);


