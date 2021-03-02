/* Wigner.i */

%{
#include "Level1/Wigner.h"
%}

%feature("autodoc", "1" );


void          Wigner_error    (int error);
void volatile Wigner_fatality (int error);
double d0();
double d1half(int n, double beta );
double dm1half(int n, double beta );
double d1half(int m, int n, double beta );
double d11(int n, double beta );
double d10(int n, double beta );
double d1m1(int n, double beta);
double d1(int m, int n, double beta);
double d22(int n, double beta);
double d21(int n, double beta );
double d20(int n, double beta );
double d2m1(int n, double beta );
double d2m2(int n, double beta );
double d2(int m, int n, double beta );
double fact(int a);
double dJ_int(int J, int m, int n, double beta );
double dJ_half_int(int J, int m, int n, double beta );
double dJ(int J, int m, int n, double beta );
matrix dJ(int J, double beta);
complex DJ(int J, int m, int n, double alpha, double beta, double gamma);
double D0();
complex D1half(int m, int n, double alpha, double beta, double gamma);
complex D1(int m, int n, double alpha, double beta, double gamma);
complex D2(int m, int n, double alpha, double beta, double gamma);
matrix DJ(int J, double alpha, double beta, double gamma);
matrix DJ(const matrix& dJbeta, int J, double alpha);
