/* HSauxil.h */

%{
#include "HSLib/HSauxil.h"
%}

#include "std_string.i"

%feature("autodoc", "1" );

gen_op sigma_eq(const spin_sys& sys);
gen_op sigma_eq(const spin_sys& sys, const Isotope& I);

void zero_mqc(const spin_sys &sys, gen_op &Op, int order, int type);
gen_op st_Op(gen_op &Ham, int lev1, int lev2, char axis);
void sqt_v(gen_op &Ham);

int* sort_super_op_basis (const spin_sys& sys);

int* sort_LOp_basis (const spin_sys& sys);
int* sort_Op_basis (const spin_sys& sys);

void mqt_v(const spin_sys& sys, gen_op &Ham, int qn, int type, int ncols);

void wavefunction(const spin_sys& sys, gen_op &Op, int wf, int pbf);
void wavefunctions(const spin_sys& sys, gen_op &Op, int pbf);
//void eigensystem(std::ostream& ostr, gen_op Op);

double vecmax(row_vector &vx);
complex integral(const row_vector& vx);
void lwhh(row_vector &vx, int& i1, int& i2);
int query_isotope(const spin_sys& sys, std::string& Isotype);

int query_isotope(const spin_sys& sys, std::string& Isotype, const std::string& Query);
int query_isotope(int argc, char* argv[], int argn, const spin_sys& sys, std::string& Isotype);

double query_offset(spin_system& sys, int isoset, int askit=0);
double query_offset(spin_system& sys, std::string& Isotype, int ask=0);

// Need to comment out the next two lines for the _pygamma.so to load into a python session.
//double query_Nyquist(spin_system& sys, int isoset, double lw=0, double fact=1.2);
//double query_Nyquist(const spin_system& sys, std::string& Isotype, double lw=0, double fact=1.2);

void query_file1D(std::string& filename, int& type);

//void query_FelixFile1D(std::ostream& ostr, std::string& filename);

//next line was already commented out before creating ".i" file
//void query_output1D(const spin_system& sys, std::string& filename, int& type, int FFT);

//void Felix1D_params(std::ostream& ostr, double Omega, double sw, int npts, double offset);
//void Felix2D_params(std::ostream& ostr, double Omega, double sw, int npts, double offset, int PPM);
//void Felix2D_params(std::ostream& ostr, double O2, double sw2, int npts2, double off2,
//                                   double O1, double sw1, int npts1, double off1, int PPM);

