/* Exponential.h */
// Swig interface file.

%{
#include "Level1/Exponential.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

class matrix;				// Knowledge of class matrix
class row_vector;			// Know about row vectors

row_vector Exponential(int npts, double W, double R);

row_vector Exponential(int npts, double time, double w, double RT, int type=0);

row_vector DExponential(int npts, double W, double R);

row_vector DExponential(int npts, double time, double w, double RT, int typ=0);
 
int Exponen_cut(int npts, double time, double w, double R, double cutoff=1.e-10);

void Exponen_cut(int* ihi,const matrix& mx,double tinc,int npts,double cutoff);

