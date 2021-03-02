/* Bloch.i */
// Swig interface file.

%{
#include "Bloch/Bloch.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

int DoubleMag(double x);

std::string SecUnits(int mag, double& sf);
std::string HzUnits(int  mag, double& sf);

matrix Mo_vector(double Mox=0, double Moy=0, double Moz=1);
///matrix Mo_vector(int argc, char* argv[], matrix& Meq, int& qn);

// Commented out as not defined in code as of 2011.05.16 (DCT)
//matrix Mss_vector(matrix& K, matrix& R, matrix& Meq);

matrix Mo_vector(int argc, char* argv[], matrix& Meq, int& qn);


void analyze(double tinc, int& ntimes,
          int& do_ss, int& qn, double T1, double gamB1, double w);

//void bloch_T1T2(const ParameterSet& pset, std::ostream& ofstr, double& T1, double& T2);

//void bloch_Mo(const ParameterSet& pset,   std::ostream& ofstr, double& Mx, double& My, double& Mz);

//void bloch_B1(const ParameterSet& pset,   std::ostream& ofstr, double& gamB1, double& phi);

//void bloch_Woff(const ParameterSet& pset, std::ostream& ofstr, double& Woff);


void TrajTiming(int argc, char* argv[],
             double& tinc, int& N, int& qn, double T1, double gamB1, double w);

