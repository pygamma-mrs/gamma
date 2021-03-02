/* relaxJ.h */

%{
#include "BWRRelax/relaxJ.h"
%}

%feature("autodoc", "1" );


void J_error (int i);

void volatile J_fatality (int error);

double J_gen(double tau, double omega, int hertz=0);

matrix J_gen(double tau, double* w, int hs, int hertz=0);

matrix J_gen_shft(double tau, double* w, double shift, int hs, int hertz=0);

void tausD(double* tau, double Dx, double Dy, double Dz);

void tausD(double* tau, const coord& Ds);

double chiD(double Dx, double Dy, double Dz);

double chiD(const coord& Ds);

void taust(double* tau, double taux, double tauy, double tauz);

void taust(double* tau, const coord& taus);

double chit(double taux, double tauy, double tauz);

double chit(const coord& taus);

void Jcoeffs(double* c, double alpha, double beta, double gamma,
						    double chi=0, double eta=0);

void Jcoeffs(double* c, const coord& EA, double chi=0, double eta=0);

double J_reduced(double *tau, double *c1, double *c2, double omega, int hertz=0);

matrix J_reduced(gen_op& Op, double *tau, double *c1, double *c2, int hertz=0);

matrix J_reduced(double* w, int size, double *tau, double *c1, double *c2, int hertz=0);

matrix J_red_shft(double* w, double shift, int size, double *tau,
						 double *c1, double *c2, int hertz=0);

matrix J_reduced(row_vector w, double *tau, double *c1, double *c2, int hertz=0);

matrix J_red_shft(row_vector w, double shift, double *tau,
					 double *c1, double *c2, int hertz=0);

double J_reduced(const coord& taus, const coord& EA1, double eta1,
			 const coord& EA2, double eta2, double omega, int hertz=0);

double J_reduced(const coord& taus, space_T& T1, space_T& T2, double omega, int hertz=0);

double Q_reduced(double *tau, double *c1, double *c2, double omega, int hertz=0);

matrix Q_reduced(gen_op& Op, double *tau, double *c1, double *c2, int hertz=0);

matrix Q_reduced(double* w, int size, double *tau, double *c1, double *c2, int hertz=0);

matrix Q_red_shft(double* w, double shift, int size, double *tau,
						 double *c1, double *c2, int hertz=0);

matrix Q_reduced(row_vector w, double *tau, double *c1, double *c2, int hertz=0);

matrix Q_red_shft(row_vector w, double shift, double *tau,
					 double *c1, double *c2, int hertz=0);

double Q_reduced(const coord& taus, const coord& EA1, double eta1, 
                   const coord& EA2, double eta2, double omega, int hertz=0);

double Q_reduced(const coord& taus, space_T& T1, space_T& T2, double omega, int hertz=0);

double J_LZ_iso(double S, double tauM, double taue, double omega);
  
double J_LZ_aniso(double S, double A, double tau1,
				double tau2, double taue, double omega);
