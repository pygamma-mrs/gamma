/* Pulse.i */
// Swig interface file.

%{
#include "Pulses/Pulse.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

%rename(__assign__) Pulse::operator=;



class Pulse
{

friend class WALTZ;			// WALTZ pulses are our friends
friend class MLEV;			// MLEV pulses are our friends too
std::string Iso;			// Applied pulse channel
double gamB1;               // Applied pulse strength (Hz)
double tp;				    // Applied pulse length   (sec)
double ang;				    // Applied pulse angle	  (deg)
double phi;                 // Applied pulse phase    (deg)
double Wrf;                 // Applied pulse offset   (Hz)


private:

void Pulerror(int eidx, int noret=0) const;

void volatile Pulsefatality(int error) const;

void Pulerror(int eidx, const std::string& pname, int noret=0) const;

void SetPhase(const ParameterSet& pset, int idx=-1);

void SetChannel(const ParameterSet& pset, int idx=-1);
 
int SetAngle(const ParameterSet& pset, int idx=-1);

int SetPulLen(const ParameterSet& pset, int idx=-1);
 
int SetGamB1(const ParameterSet& pset, int idx=-1);

 
public:
 
 
Pulse();

Pulse(const std::string& ch, double gB1, double len, double ph=0, double off=0);

Pulse(const Pulse& PT1);
 
virtual ~Pulse();
                                                                                
Pulse& operator = (const Pulse& Pulse1);

std::string channel()  const;
double      strength() const;
double      angle()    const;
double      phase()    const;
double      offset()   const;
double      length()   const;

void strength(double gB1);

///HSprop GetU(const spin_system& sys, gen_op& H);
 
void read(const std::string &filename, int idx=-1);

void read(const ParameterSet& pset, int idx=-1);

std::string ask_read(int argc, char* argv[], int argn);

//std::ostream& printBase(std::ostream &ostr) const;
//virtual std::ostream& print(std::ostream &ostr, int full=0) const;  
//friend std::ostream &operator << (std::ostream &ostr, const Pulse &Pul);

};

