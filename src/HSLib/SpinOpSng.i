/* SpinOpSng.i */
// Swig interface file

%{
#include "HSLib/SpinOpSng.h"
%}

%feature("autodoc", "1" );

matrix Ie(int qn);		// Identity matrix of dimension qn
matrix Ix(int qn);		// Ix matrix of dimension qn
matrix Iy(int qn);		// Iy matrix of dimension qn
matrix Iz(int qn);		// Iz matrix of dimension qn
matrix Ip(int qn);		// I+ matrix of dimension qn
matrix Im(int qn);		// I- matrix of dimension qn

matrix Raxis(int qn, double beta, char axis);
 
