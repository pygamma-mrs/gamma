/* PulDANTE.i */
// Swig interface file.

%{
#include "Pulses/PulDANTE.h"
%}

%include "std_string.i"

%rename(__assign__) DANTE::operator=;

class PulWaveform;
class PulComposite;
class PulTrain;

class DANTE
{
  int N;                                // Number of steps
  std::string Iso;			// Applied pulse channel
  double te;				// Delay evolution time   (sec)
  double gamB1;                         // Applied pulse strength (Hz)
  double tp;				// Applied pulse length   (sec)
  double ang;				// Applied pulse angle	  (deg)
  double phi;                           // Applied pulse phase    (deg)
  double Wrf;                           // Applied pulse offset   (Hz)
  double tt;				// Full DANTE step length (sec)
  double F;				// DANTE synch frequency  (Hz)

private:

void DANTEerror(int eidx, int noret=0) const;
void volatile DANTEfatality(int error) const;
void DANTEerror(int eidx, const std::string& pname, int noret=0) const;
void SetSteps(const ParameterSet& pset, int idx=-1);
void SetPhase(const ParameterSet& pset, int idx=-1);
void SetChannel(const ParameterSet& pset, int idx=-1);

int SetAngle(const ParameterSet& pset, int idx=-1);
int SetPulLen(const ParameterSet& pset, int idx=-1);
int SetGamB1(const ParameterSet& pset, int idx=-1);
int SetEvLen(const ParameterSet& pset, int idx=-1);
int SetFreq(const ParameterSet& pset, int idx=-1);


public:

DANTE();

DANTE(const DANTE& PT1);

~DANTE();
                                                                               
DANTE& operator = (const DANTE& DANTE1);

int        steps()    const;
std::string     channel()  const;
double     dlength()  const;
double     strength() const;
double     plength()  const;
double     angle()    const;
double     phase()    const;
double     offset()   const;
double     length()   const;
 
//friend PulWaveform WF_DANTE(const DANTE& D);
PulWaveform WF( ) const;

//friend PulComposite CP_DANTE(const spin_system& sys, const DANTE& D);
PulComposite CP(const spin_system& sys) const;

//friend PulTrain PT_DANTE(const spin_system& sys, const DANTE& D);

PulTrain PT(const spin_system& sys) const;
 
void read(const std::string &filename, int idx=-1);
  
void read(const ParameterSet& pset, int idx=-1);

void ask_read(int argc, char* argv[], int argn);

//std::ostream& printBase(std::ostream &ostr) const;
//std::ostream& printInfo(std::ostream &ostr) const;                                                                         
//std::ostream& print(std::ostream &ostr, int full=0) const;
//friend std::ostream &operator << (std::ostream &ostr, const DANTE &DT);

};

PulWaveform WF_DANTE(double td, double gamB1, double tpul, double phi=0);

PulComposite CP_DANTE(const spin_system& sys, const std::string& Iso,
                           double td, double gamB1, double tpul, double phi=0);

PulComposite CP_DANTE(const spin_system& sys, const DANTE& D);

PulTrain PT_DANTE(const spin_system& sys, const std::string& Iso,
                           double td, double gamB1, double tpul, double phi=0);

// Call to this next function can be ambiguous in SWIG, so left out.

//gen_op UDANTE(const spin_system& sys, gen_op& H, const std::string& Iso,
//                                     double td, double theta, double phi=0.0);

gen_op UDANTE(const spin_system& sys, gen_op& H, const std::string& Iso,
                      double td, double gamB1, double tpul, double phi);

gen_op UDANTE(const spin_system& sys, gen_op& H, const DANTE& D);

double ask_DANTE(const spin_system& sys, const std::string& Iso, gen_op& H, double cutoff=1.e-10);

void set_DANTE(double gamB1, double& tmix, double& tpul, double tdel, int& numb, int& type);

