/* PulWALTZ.i */
// Swig interface file.

%{
#include "Pulses/PulWALTZ.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

%rename(__assign__) PulWALTZ::operator=;


class WALTZ : public Pulse
{
 
private:
                                                                               
void WALTZerror(int eidx, int noret=0) const;

void volatile WALTZfatality(int error) const;
                                                                       
void WALTZerror(int eidx, const std::string& pname, int noret=0) const;

void SetPhase(const ParameterSet& pset, int idx=-1);

void SetChannel(const ParameterSet& pset, int idx=-1);
                                                     
int SetGamB1(const ParameterSet& pset, int idx=-1);


public:

WALTZ();

WALTZ(double gB1, const std::string& ch, double ph=0, double off=0);

WALTZ(const WALTZ& WP1);

~WALTZ();
                                                                              
WALTZ& operator = (const WALTZ& WALTZ1);

PulWaveform WF(int even=0) const;
PulWaveform WF_WALTZR(int even=0) const;

PulWaveform WF_WALTZK(int even=0) const;

PulWaveform WF_WALTZQ(int even=0) const;

PulComposite PCmp(const spin_system& sys, int even=0) const;
PulComposite PCmpWALTZR(const spin_system& sys, int even=0) const;
 
PulComposite PCmpWALTZK(const spin_system& sys, int even=0) const;

PulComposite PCmpWALTZQ(const spin_system& sys, int even=0) const;

PulCycle CycWALTZ4(const spin_system& sys, int even=0) const;

PulCycle CycWALTZ8(const spin_system& sys, int even=0) const;

PulCycle CycWALTZ16(const spin_system& sys, int even=0) const;

void read(const std::string &filename, int idx=-1);

void read(const ParameterSet& pset, int idx=-1);

void ask_read(int argc, char* argv[], int argn);
                                                             
//std::ostream& print(std::ostream &ostr) const;
//friend std::ostream &operator << (std::ostream &ostr, const WALTZ &GP);

};


row_vector CYC_WALTZ4(double phi=0);

row_vector CYC_WALTZ8(double phi=0);

