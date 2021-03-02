/* PulGauss.i */
// Swig interface file.

%{
#include "Pulses/PulGauss.h"
%}

%include "std_string.i"

struct Gpuldat                          // Gpulse info
{
  int N;                                // Number of steps
  double Wrf;                           // Field frequency
  std::string Iso;                      // Field selectivity
  double gamB1;                         // Field strength (Hz)
  double tau;                           // Pulse length (sec)
  double fact;                          // Cutoff values (%gamB1)
  double phi;                           // Pulse phase (degrees)
};

// :FIXME: This code gives an LNK2019 (linker error)
// on Windows Visual Studio Express, and causes problems 
// when the final pygamma library is loaded into Python
// on Linux.
// After further study the function signature included 
// in the .i and .h is not defined in the .cc. Seems
// like this was never working in the first place, so
// will leave it commented out, until further notice.
//
//gen_op Gaussian(spin_system& sys, gen_op& H, std::string& Iso,
//                        double td, double theta, double phi=0.0);

                        
void Gpulse_Hs(gen_op* Hs, gen_op& H0, gen_op& Fxy,
                          int N, double ang, double tp, double fact);                          
                          
void Gpulse_Us(gen_op* Us, gen_op& H0, gen_op& Fxy,
                          int N, double ang, double tp, double fact);

gen_op Gpulse_U(const spin_system& sys, const Gpuldat& Gdata);

gen_op Gpulse_U(const spin_sys& sys, gen_op& H0, const Gpuldat& Gdata);

gen_op Gpulse_U(gen_op& H0rot, gen_op& Fxy, const Gpuldat& Gdata);

gen_op Gpulse_U(gen_op& H0rot, gen_op& Fxy,
                          int N, double ang, double tp, double fact);
      
gen_op Gpulse_UX(gen_op& H0rot, gen_op& Fxy,
                          int N, double ang, double tp, double fact);


double Gangle(double gamB1, double tau, int N, double fact=0.025);

double GgamB1(double angle, double tau, int N, double fact=0.025);

double Gtime(double angle, double gamB1, int N, double fact=0.025);

row_vector GNvect(int N, double fact);

row_vector Gvect(double gamB1, int N, double fact=0.025);

row_vector GIntvec(double gamB1, double tau, int Npts, double fact=0.05);

row_vector Ghistogram(double gamB1, double tau, int N, double fact=0.05);

void ask_Gpulse(int argc, char* argv[], int& qn, int& N,
                 double& val1, double& val2, double& fact, const int type=0);

Gpuldat read_Gpulse(std::string& filein, spin_system& sys, int idx=-1, int pflag=1);

void read_Gpulse(std::string& filein, int& N, double& val1, double& val2,
                         double& fact, int& spin, int type=0, int idx=-1);
 
//void print_Gpulse(std::ostream& ostr, Gpuldat& Gdata, int indx=-1);
