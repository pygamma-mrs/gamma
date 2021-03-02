/* PulAuxil.i */
// Swig interface file.

%{
#include "Pulses/PulAuxil.h"
%}

%include "std_string.i"


row_vector pulseshift(row_vector& p, row_vector& ptime, const double& FreqOffset);
 
