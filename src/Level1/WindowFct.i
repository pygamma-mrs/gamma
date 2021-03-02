/* WindowFct.i */

%{
#include "Level1/WindowFct.h"
%}

%feature("autodoc", "1" );

void exponential_multiply(col_vector &db, double em=-4, int offset=0);
void exponential_multiply(row_vector &db, double em=-4, int offset=0);
row_vector exponential (int size, int offset=0, double alpha=0);
row_vector Gaussian (int size, int offset=0, double sigma=42.0);
row_vector Hamming (int size, int offset=0);
row_vector Hanning (int size, int offset=0);
row_vector hyperbol_sec (int size, int offset=0, double alpha=38.0);
row_vector Kaiser (int size, double theta=PI, int offset=0);
row_vector Lorentzian (int size, int offset=0, double alpha=1.0);
row_vector sin_square (int size, int offset=0);
row_vector sinc (int size, int offset, int inc);
row_vector square_wave (int size, int start, int finish);
row_vector Noise(int npts,         double maxN);
void       Noise(row_vector& data, double maxN);
