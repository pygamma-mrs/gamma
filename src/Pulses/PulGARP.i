/* PulGARP.i */
// Swig interface file.

%{
#include "Pulses/PulGARP.h"
%}

%include "std_string.i"

%rename(__assign__) GARP::operator=;


class GARP
{
  std::string Iso;                           // Applied pulse channel
  double gamB1;                         // Applied pulse strength (Hz)
  double phi;                           // Applied pulse phase    (deg)
  double Wrf;                           // Applied pulse offset   (Hz)

private:

void GARPerror(int eidx, int noret=0) const;
                                                               
void volatile GARPfatality(int error) const;

void GARPerror(int eidx, const std::string& pname, int noret=0) const;

void SetPhase(const ParameterSet& pset, int idx=-1);

void SetChannel(const ParameterSet& pset, int idx=-1);
 
int SetGamB1(const ParameterSet& pset, int idx=-1);


public:

GARP();
                                              
GARP(double gB1, const std::string& ch, double ph=0, double off=0);

GARP(const GARP& PT1);
         
~GARP();

GARP& operator = (const GARP& GARP1);

std::string channel()  const;
double strength() const;
void   strength(double gB1);
double phase()    const;
double offset()   const;

PulWaveform WF( ) const;
PulWaveform WF_GARP( ) const;
 
PulComposite PCmp(const spin_system& sys) const;
PulComposite PCmpGARP(const spin_system& sys) const;
 
PulComposite PCmp(const spin_system& sys, const super_op& LOp) const;

PulCycle CycGARP1(const spin_system& sys) const;

void read(const std::string &filename, int idx=-1);

void read(const ParameterSet& pset, int idx=-1);
 
void ask_read(int argc, char* argv[], int argn, int idx=-1);

//std::ostream& printBase(std::ostream &ostr) const;
//std::ostream& print(std::ostream &ostr) const;
//friend std::ostream &operator << (std::ostream &ostr, const GARP &GP);

};

