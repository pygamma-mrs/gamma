/* HSdetect.i */
// Swig interface file.

%{
#include "HSLib/HSdetect.h"
%}

%feature("autodoc", "1" );

class spin_op;				// Know spin_op is a class
class spin_sys;				// Know spin_sys is a class

spin_op GenericD(const spin_sys& sys, spin_op D, double beta=0);

spin_op detector(const spin_sys &S, double B=0);
spin_op Mxy(const spin_sys &sys, double beta=0);

spin_op detector(const spin_sys &S, char *I, double B=0);
spin_op Mxy(const spin_sys &sys, char *iso, double beta=0);

spin_op detector(const spin_sys &S, int N, double B=0);
spin_op Mxy(const spin_sys &sys, int spin, double beta=0);

spin_op detector_sp(const spin_sys &S, double B=0);
spin_op Mxy_sp(const spin_sys &sys, double beta=0);
