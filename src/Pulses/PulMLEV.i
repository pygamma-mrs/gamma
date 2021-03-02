/* PulMLEV.i */
// Swig interface file.

%{
#include "Pulses/PulMLEV.h"
%}

%include "std_string.i"

%rename(__assign__) MLEV::operator=;


class MLEV : public Pulse
{

private:
 
void MLEVerror(int eidx, int noret=0) const;
                            
void volatile MLEVfatality(int error) const;

void MLEVerror(int eidx, const std::string& pname, int noret=0) const;
 
void volatile MLEVfatality(int eidx, const std::string& pname, int noret=0) const;
                                                                                
void SetPhase(const ParameterSet& pset, int idx=-1);
                                                         
void SetChannel(const ParameterSet& pset, int idx=-1);

int SetGamB1(const ParameterSet& pset, int idx=-1);


public:                                                                        
                                                                                
MLEV();
                                              
MLEV(double gB1, const std::string& ch, double ph=0, double off=0);
                                                         
MLEV(const MLEV& PT1);

~MLEV();

MLEV& operator = (const MLEV& MLEV1);
 

///std::string channel()  const;                 INHERITED
///double strength()      const;                 INHERITED
///double phase()         const;                 INHERITED
///double offset()        const;                 INHERITED
 
PulWaveform WF( ) const;
PulWaveform WF_C180( ) const;
 
PulComposite PCmp(const spin_system& sys) const;
PulComposite PCmp_C180(const spin_system& sys) const;
 
PulCycle CycMLEV4(const spin_system& sys) const;
 
PulCycle CycMLEV8(const spin_system& sys) const;
 
PulCycle CycMLEV16(const spin_system& sys) const;
       
void read(const std::string &filename, int idx=-1);
                                                         
void read(const ParameterSet& pset, int idx=-1);
 
void ask_read(int argc, char* argv[], int argn);
         
//std::ostream& printBase(std::ostream &ostr) const;
//std::ostream& print(std::ostream &ostr) const;       
//friend std::ostream &operator << (std::ostream &ostr, const MLEV &GP);

};

 
row_vector CYC_MLEV4(double phi=0);
 
row_vector CYC_MLEV8(double phi=0);
 
row_vector CYC_MLEV16(double phi=0);

