/* PulTrainSCyc.i */
// Swig interface file.

%{
#include "Pulses/PulTrainSCyc.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

%rename(__assign__) PulTrainSCyc::operator =;


class PulTrainSCyc : public PulSupCycle
{
  double      tp;			// Pulse train cycle length (sec)
  int         Ucount;			// Pulse train propagator count
  HSprop*     Usteps;			// Pulse train cycle step propagators
  HSprop*     Usums;			// Pulse train cycle summed propagators
  LSprop*     Gsteps;			// Pulse sequence step superpropagators


private:


void PTSCerror(int eidx, int noret=0) const;

void volatile PTSCfatality(int error) const;

void deleteHams();

void deleteUprops();

void deleteGprops();

void copyHams(const PulTrainSCyc& PTSC1);

void copyUprops(const PulTrainSCyc& PTSC1);

void copyGprops(const PulTrainSCyc& PTSC1);

void SetUs(PulCycle& PTC);


public:

PulTrainSCyc();

PulTrainSCyc(PulCycle& PTC, const row_vector& S, std::string N);

PulTrainSCyc(PulCycle& PTC, const PulSupCycle& PCYC);

PulTrainSCyc(const PulTrainSCyc& PT1);

~PulTrainSCyc();

PulTrainSCyc& operator = (const PulTrainSCyc& PTSC1);

HSprop GetU(int i=-1) const;

HSprop GetUsum(int i=-1) const;

HSprop GetUmult(int N) const;

void SetGs(PulCycle& PTC);

LSprop GetG(int i=-1) const;

///int        steps()   const;
///std::string     name()    const;
double     length()  const;
///row_vector values()  const;

int        steps()   const;

///double  phase(int i)    const;

double steps(double td) const;

int fullSCYCs(double td=-1) const;
int fullsteps(double td=-1) const;

//std::ostream& printInfo(std::ostream &ostr) const;
//std::ostream& printBase(std::ostream &ostr) const;
//std::ostream& print(std::ostream &ostr, int full=0) const;
//friend std::ostream &operator << (std::ostream &ostr, const PulTrainSCyc &PTSC);

};

