/* PulCHIRP.i */
// Swig interface file.

%{
#include "Pulses/PulCHIRP.h"
%}

PulWaveform WF_CHIRP95(int N, double tp, double delW, double gB1, int scale=0);

PulComposite CP_CHIRP95(const spin_system& sys, const std::string& IsoC,
                       int N, double tp, double delW, double gB1, int scale=0);
  
PulCycle CYC_CHIRP95();
