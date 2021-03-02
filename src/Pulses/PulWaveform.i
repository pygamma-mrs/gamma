/* PulWaveform.i */
// Swig interface file.

%{
#include "Pulses/PulWaveform.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

%rename(__assign__) PulWaveform::operator=;


class PulWaveform
{

public:

PulWaveform();
PulWaveform(const row_vector& wfsteps, const row_vector& wftimes,
                         const std::string& wfname="", int wfrad=0);
PulWaveform(const PulWaveform& PT1);

virtual ~PulWaveform();

PulWaveform& operator = (const PulWaveform& PWF1);
virtual int         steps()   const;
virtual std::string name()    const;
virtual double      length()  const;
virtual row_vector  values()  const;
virtual row_vector  lengths() const;

        double  strength(int i) const;
virtual double  phase(int    i) const;
        double  length(int   i) const;
virtual complex value(int    i) const;

double maxlength( ) const;

double minlength(double cutoff=1.e-13) const;

double maxgamB1( ) const;

double mingamB1( ) const;

bool gamB1const() const;
bool phaseconst() const;
bool timeconst()  const;

double steps(double td) const;

int fullsteps(double td) const;

double WFs(double td) const;

int fullWFs(double td, double cut=1.e-13) const;

double sumlength(int i) const;

virtual void scalegB1(double sf);

void getIdeal(double& gB1, double& ptt, int i) const;

row_vector IvsT(int split=0, int ends=0, int N=1) const;

row_vector PvsT(int spl=0, int ends=0, int N=1, double p=0) const;

void GP(int type=0, int split=0, int ends=0, int N=1) const;

void FM(int type=0, int split=0, int ends=0, int N=1 ) const;

//std::ostream& printBase(std::ostream& ostr) const;

//std::ostream& printSteps(std::ostream& ostr, int full=0) const;

//virtual std::ostream& print(std::ostream& ostr, int full=0) const;

//friend std::ostream& operator << (std::ostream& ostr, const PulWaveform& PWF);

};


