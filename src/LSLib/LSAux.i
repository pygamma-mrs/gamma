/* LSAux.i */

%{
#include "LSLib/LSAux.h"
%}

%feature("autodoc", "1" );

%rename(lsprint) print;

void print(const super_op& LOp, double cutoff=1.e-6, int nc=4, int ri=0);

void eigenvalues(super_op& LOp, int sort=1, int nc=4, int ri=0);

matrix UOrderMQC(const spin_sys &sys);

super_op OrderMQC(const super_op &LOp, const matrix &U);

super_op OrderMQC(super_op &LOp, spin_sys &sys);

matrix solve_it(matrix& X, matrix& Uguess, matrix& b, int lim=10);

matrix invert_it(matrix& X);

int LU_decomp(matrix& A, int* indx);

void LU_backsub(matrix &ALU, int* indx, matrix& b);

matrix LU_invert(matrix& A);
