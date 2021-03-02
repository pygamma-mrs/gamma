/* SpinT.i */
// Swig interface file.

%{
#include "Level1/SpinT.h"
%}

%feature("autodoc", "1" );

%rename(__assign__) spin_T::operator=;

class spin_T;


spin_T T1(const spin_sys &sys, int spin);
 
spin_T T11(const spin_sys &sys, int spin);

spin_op T1(const spin_sys &sys, int spin, int l, int m);

spin_op T10(const spin_sys &sys, int spin, int m);

spin_op T10(spin_op &Ie, int m);

spin_op T11(const spin_sys &sys, int spin, int m);

spin_op T11(spin_op &Im, spin_op &Iz, spin_op &Ip, int m);

spin_T T2(const spin_sys &sys, int spin1, int spin2);

spin_T T22wh(const spin_sys &sys, int spin1, int spin2);

spin_T T22(const spin_sys &sys, int spin1, int spin2);

spin_T T22(const spin_sys &sys, spin_op &Im1, spin_op &Iz1, spin_op &Ip1,
			           spin_op &Im2, spin_op &Iz2, spin_op &Ip2);

spin_op T2(const spin_sys &sys, int spin1, int spin2, int l, int m);

spin_op T2(spin_op &Im1, spin_op &Iz1, spin_op &Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int l, int m);

spin_op T20(const spin_sys &sys, int spin1, int spin2, int m);

spin_op T20(spin_op &Im1, spin_op &Iz1, spin_op& Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int m);

spin_op T21(const spin_sys &sys, int spin1, int spin2, int m);

spin_op T21(spin_op &Im1, spin_op &Iz1, spin_op& Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int m);

spin_op T22(const spin_sys &sys, int spin1, int spin2, int m);

spin_op T22(spin_op &Im1, spin_op &Iz1, spin_op& Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int m);

spin_T T2(const spin_sys &sys, int spin, const coord &vect);
spin_T T2(const spin_sys &sys, spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect);

spin_T T2SS(const spin_sys &sys, int spin, const coord &vect, int rev=0);

spin_T T2SS(const spin_sys &sys, spin_op &Im, spin_op &Iz,
				      spin_op &Ip, const coord &vect, int rev=0);

spin_T T22(const spin_sys &sys, int spin, const coord &vect);
spin_T T22(const spin_sys &sys, spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect);

spin_T T22SSirr(const spin_sys &sys, int spin, const coord &vect, int rev=0);

spin_T T22SSirr(const spin_sys &sys, spin_op &Im, spin_op &Iz,
					 spin_op &Ip, const coord &vect, int rev=0);

spin_op T2(const spin_sys &sys, int spin, const coord &vect, int l, int m);
spin_op T2(spin_op &Im, spin_op &Iz, spin_op &Ip,
 					   const coord &vect, int l, int m);
spin_op T20(const spin_sys &sys, int spin, const coord &vect, int m);
spin_op T20(spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect, int m);
spin_op T21(const spin_sys &sys, int spin, const coord &vect, int m);
spin_op T21(spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect, int m);
spin_op T22(const spin_sys &sys, int spin, const coord &vect, int m);
spin_op T22(spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect, int m);

spin_op T2SS(const spin_sys &sys, int spin, const coord &vect, int l, int m, int rev=0);

spin_op T2SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
 				   const coord &vect, int l, int m, int rev=0);

spin_op T20SS(const spin_sys &sys, int spin, const coord &vect, int m, int rev=0);

spin_op T20SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
						 const coord &vect, int m, int rev=0);

spin_op T21SS(const spin_sys &sys, int spin, const coord &vect, int m, int rev=0);

spin_op T21SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
					 const coord &vect, int m, int rev=0);

spin_op T22SS(const spin_sys &sys, int spin, const coord &vect, int m, int rev=0);

spin_op T22SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
						const coord &vect, int m, int rev=0);

spin_op T_comp(spin_T &SphT, int L, int M);

spin_T T_mult(spin_T &SphT1, spin_T &SphT2);

spin_op T_mult(spin_T &SphT1, spin_T &SphT2, int L, int M);

spin_T T_rot(spin_T &SphT1, double alpha, double beta, double gamma);

spin_op T_rot(spin_T &SphT1, int l, int m,
			 double alpha, double beta, double gamma);

spin_op T_prod(spin_T &SphT, space_T &SphA, int l, int m);
spin_op T_prod(space_T &SphA, spin_T &SphT, int l, int m);

spin_op T_prod(spin_T &SphT, space_T &SphA, int l);
spin_op T_prod(space_T &SphA, spin_T &SphT, int l);

spin_op T_prod(spin_T &SphT, space_T &SphA);
spin_op T_prod(space_T &SphA, spin_T &SphT);

double Clebsch_Gordan(int a, int b, int alpha, int beta, int c, int gamma);
double Wigner_3j(int a, int b, int c, int alpha, int beta, int gamma);



class spin_T
{

spin_sys* sys; 		// Pointer to a spin system
int rank;			// Rank of spin tensor
spin_op*** pr;		// pointer to pointer to SOp pointers

void volatile spin_T_error(const int i);

void volatile spin_T_fatality (int error);


public:

spin_T();

spin_T(const spin_sys &sys);

spin_T(const spin_T &SphT);

spin_T(const spin_T &SphT, int l);

~spin_T();

spin_T & operator = (const spin_T & SphT);


/* Commenting out friend functions... (DCT)
friend spin_T T1(const spin_sys &sys, int spin);
friend spin_T T11(const spin_sys &sys, int spin);
friend spin_op T1(const spin_sys &sys, int spin, int l, int m);
friend spin_op T10(const spin_sys &sys, int spin, int m);
friend spin_op T10(spin_op &Ie, int m);
friend spin_op T11(const spin_sys &sys, int spin, int m);
friend spin_op T11(spin_op &Im, spin_op &Iz, spin_op &Ip, int m);
friend spin_T T2(const spin_sys &sys, int spin1, int spin2);
friend spin_T T22wh(const spin_sys &sys, int spin1, int spin2);
friend spin_T T22(const spin_sys &sys, int spin1, int spin2);
friend spin_T T22(const spin_sys &sys, spin_op &Im1, spin_op &Iz1, spin_op &Ip1,
			           spin_op &Im2, spin_op &Iz2, spin_op &Ip2);
friend spin_op T2(const spin_sys &sys, int spin1, int spin2, int l, int m);
friend spin_op T2(spin_op &Im1, spin_op &Iz1, spin_op &Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int l, int m);
friend spin_op T20(const spin_sys &sys, int spin1, int spin2, int m);
friend spin_op T20(spin_op &Im1, spin_op &Iz1, spin_op& Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int m);
friend spin_op T21(const spin_sys &sys, int spin1, int spin2, int m);
friend spin_op T21(spin_op &Im1, spin_op &Iz1, spin_op& Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int m);
friend spin_op T22(const spin_sys &sys, int spin1, int spin2, int m);
friend spin_op T22(spin_op &Im1, spin_op &Iz1, spin_op& Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int m);
friend spin_T T2(const spin_sys &sys, int spin, const coord &vect);
friend spin_T T2(const spin_sys &sys, spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect);
friend spin_T T2SS(const spin_sys &sys, int spin, const coord &vect, int rev);
friend spin_T T2SS(const spin_sys &sys, spin_op &Im, spin_op &Iz,
				      spin_op &Ip, const coord &vect, int rev);
friend spin_T T22(const spin_sys &sys, int spin, const coord &vect);
friend spin_T T22(const spin_sys &sys, spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect);
friend spin_T T22SSirr(const spin_sys &sys, int spin, const coord &vect, int rev);
friend spin_T T22SSirr(const spin_sys &sys, spin_op &Im, spin_op &Iz,
					 spin_op &Ip, const coord &vect, int rev);
friend spin_op T2(const spin_sys &sys, int spin, const coord &vect, int l, int m);
friend spin_op T2(spin_op &Im, spin_op &Iz, spin_op &Ip,
 					   const coord &vect, int l, int m);
friend spin_op T20(const spin_sys &sys, int spin, const coord &vect, int m);
friend spin_op T20(spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect, int m);
friend spin_op T21(const spin_sys &sys, int spin, const coord &vect, int m);
friend spin_op T21(spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect, int m);
friend spin_op T22(const spin_sys &sys, int spin, const coord &vect, int m);
friend spin_op T22(spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect, int m);
friend spin_op T2SS(const spin_sys &sys, int spin, const coord &vect, int l, int m, int rev);
friend spin_op T2SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
 				   const coord &vect, int l, int m, int rev);
friend spin_op T20SS(const spin_sys &sys, int spin, const coord &vect, int m, int rev);
friend spin_op T20SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
						 const coord &vect, int m, int rev);
friend spin_op T21SS(const spin_sys &sys, int spin, const coord &vect, int m, int rev);
friend spin_op T21SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
					 const coord &vect, int m, int rev);
friend spin_op T22SS(const spin_sys &sys, int spin, const coord &vect, int m, int rev);
friend spin_op T22SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
						const coord &vect, int m, int rev);
friend spin_op T_comp(spin_T &SphT, int L, int M);
*/

spin_op component(int l, int m);

//friend spin_T T_mult(spin_T &SphT1, spin_T &SphT2);
//friend spin_op T_mult(spin_T &SphT1, spin_T &SphT2, int L, int M);
//friend spin_T T_rot(spin_T &SphT1, double alpha, double beta, double gamma);
//friend spin_op T_rot(spin_T &SphT1, int l, int m,
//			 double alpha, double beta, double gamma);

spin_T rotate(double alpha, double beta, double gamma);

spin_T rotate(const coord &EA);

spin_op rotate(int l, int m, double alpha, double beta, double gamma);

spin_op rotate(int l, int m, const coord &EA);

/*
friend spin_op T_prod(spin_T &SphT, space_T &SphA, int l, int m);
friend spin_op T_prod(space_T &SphA, spin_T &SphT, int l, int m);
friend spin_op T_prod(spin_T &SphT, space_T &SphA, int l);
friend spin_op T_prod(space_T &SphA, spin_T &SphT, int l);
friend spin_op T_prod(spin_T &SphT, space_T &SphA);
friend spin_op T_prod(space_T &SphA, spin_T &SphT);
*/

int Rank();

//friend double Clebsch_Gordan(int a, int b, int alpha, int beta, int c, int gamma);
//friend double Wigner_3j(int a, int b, int c, int alpha, int beta, int gamma);

//friend std::ostream& operator<< (std::ostream& ostr, const spin_T &SphT);

};
