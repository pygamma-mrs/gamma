/* Lorentzian.i */
// Swig interface file.

%{
#include "Level1/Lorentzian.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

complex Lorentzian(double Weval, double W, double R);

row_vector Lorentzian(int npts, double Wi, double Wf, double W, double R);
 
row_vector Lorentzian(int npts, double W, double R, int max1=0);

row_vector Lorentzian(double R, double W,
           double wst, double wfi, int N, complex& I, double Lcut, int Lint);

void Lorentzian(row_vector& data, double R, double W,
                  double wst, double wfi, complex& I, double Lcut=1.e-4, int Lint=5);

void Lorentzian(int* Lint, const matrix& mx, double winc, int N=5);

void Lorentz_cut(int& ilo, int& ihi, double R, double W,
                     double w0, double winc, int npts, double cutoff=1.e-4);

void Lorentz_cut(int* ilo, int* ihi, const matrix& R,
                    double w0, double winc, int npts, double cutoff=1.e-4);

void Lorentz_int(int& Lint, double R, double winc, int N=5);

void Lorentz_int(int* Lint, const matrix& R, double winc, int N=5);
 
void ask_Lorentzian(int argc, char* argv[], int& qn, int& N,
                  double& wst, double& wfi, double& W, double& R,
                                                 double& fact, int& pplw);
 
void read_Lorentzian(const std::string& filein, int& N, double& wst, double& wfi,
                 double& W, double& R, double& fact, int& pplw, int idx=-1);

complex DLorentzian(double Weval, double W, double R);

row_vector DLorentzian(int npts, double Wi, double Wf, double W, double R);

row_vector DLorentzian(int npts, double offset, double R, int max1=0);

