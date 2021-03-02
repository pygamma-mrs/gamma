/* PulSinc.i */
// Swig interface file.

%{
#include "Pulses/PulSinc.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

struct SincPulDat 			// Sinc Pulse Info
{
  int N;                    // Number of steps
  double Wrf;               // Field frequency
  std::string Iso;          // Field selectivity
  double gamB1;             // Field strength (Hz)
  double tau;               // Pulse length (sec)
  int node; 				// Cutoff node
  double phi;               // Pulse phase (degrees)
};

void SincPulseHs(gen_op* Hs, gen_op& H0, gen_op& Fxy,
                        int N, double ang, double tp, int node=3);

void SincPulseUs(gen_op* Us, gen_op& H0, gen_op& Fxy,
                        int N, double ang, double tp, int node=3);

gen_op SincPulseU(gen_op& H0rot, gen_op& Fxy,
                        int N, double ang, double tp, int node=3);

gen_op SincPulseU(gen_op& H0rot, gen_op& Fxy, const SincPulDat& SD);

double SincAngle(double gamB1, double tau,   int N, int node=3);
double SincGamB1(double angle, double tau,   int N, int node=3);
double SincTime(double angle,  double gamB1, int N, int node=3);

row_vector SincNVect(const SincPulDat& SD,                     int endzero=0);
row_vector SincNVect(                       int N, int node=3, int endzero=0);
row_vector SincVect(const SincPulDat& SD,                      int endzero=0);
row_vector SincVect(double gamB1,           int N, int node=3, int endzero=0);
row_vector SincIntVec(const SincPulDat& SD,                    int endzero=0);
row_vector SincIntVec(double gB1,double tp, int N, int node=3, int endzero=0);
  
int SincSteps(const ParameterSet& pset, int idx=-1, int pf=0);
int SincNode(const ParameterSet& pset,  int idx=-1, int pf=0);
int SincStrength(const ParameterSet& pset,SincPulDat& SD,int idx=-1,int pf=0);
int SincSelectivity(const ParameterSet& pset, const spin_system& sys,
                                     SincPulDat& Sdata, int idx=-1, int pf=0);

double SincPhase(const ParameterSet& pset, int idx, int pf=0);
SincPulDat ReadSinc(const std::string& filein,
                                const spin_system& sys, int idx=-1, int pf=0);

row_vector SincHistogram(double gamB1, double tau, int N, int node=3);

void SincPts(int argc,   char* argv[], int& qn, SincPulDat& SD);
void SincNode(int argc,  char* argv[], int& qn, SincPulDat& SD);
void SincTime(int argc,  char* argv[], int& qn, SincPulDat& SD);
void SincGamB1(int argc, char* argv[], int& qn, SincPulDat& SD);
void SincAngle(int argc, char* argv[], int& qn, double& pang);
void SincIso(int argc,   char* argv[], int& qn, SincPulDat& SD);
void SincWrf(int argc,   char* argv[], int& qn, SincPulDat& SD);
void SincPhi(int argc,   char* argv[], int& qn, SincPulDat& SD);
SincPulDat SincAsk(int argc, char* argv[], int& qn, int type=0);

//void SincPrint(std::ostream& ostr, const SincPulDat& Gdata, int indx=-1);

void SincZero(SincPulDat& SincData);
void SincPrep(int& N, int& node, double& den, double& Z);
