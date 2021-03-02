/* PulSupCycle.h */
// Swig interface file.

%{
#include "Pulses/PulSupCycle.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

%rename(__assign__) PulSupCycle::operator=;


class PulSupCycle
{

friend class PulTrainSCyc;	// Allow this class full access
std::string  SCycname;	    // Pulse supercycle name
int          SCycnosteps;	// Pulse supercycle # of steps
row_vector   SCycsteps;		// Pulse supercycle step phase (deg)


private:

void SCycerror(int eidx, int noret=0) const;

void volatile SCycfatality(int error) const;


public:

PulSupCycle();
 
PulSupCycle(const row_vector& ptsteps, const std::string& ptname="");
 
PulSupCycle(const PulSupCycle& PT1);

PulSupCycle(const PulCycle& Cyc);

~PulSupCycle();

PulSupCycle& operator = (const PulSupCycle& SCyc1);

//void operator = (const PulCycle& Cyc);

int             steps()   const;
std::string     name()    const;
row_vector      values()  const;
 
complex value(int i)  const;
double  phase(int i)  const;

///void GP(int split=0, int ends=0) const;
///void GP(const PulWaveform& PW, int split, int ends) const;

//std::ostream& printBase(std::ostream &ostr, double SCyclen=0) const;
//std::ostream& printSteps(std::ostream &ostr) const;
//std::ostream& print(std::ostream &ostr, int full=0) const;
//friend std::ostream &operator << (std::ostream &ostr, const PulSupCycle &SCyc);

};

