/* Wigner.h *****************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Wigner Rotation Matrices                    Interface		**
**                                                                      **
**      Copyright (c) 1990                                              **
**      Scott Smith                                                     **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** This module provides functions to access Wigner rotation matrices    **
** and Wigner rotation matrix elements.                                 *
**                                                                      **
*************************************************************************/

#ifndef   GWigner_h_			// Is file already included?
#  define GWigner_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/complex.h>
#include <Matrix/matrix.h>
#include <iostream>
#include <Basics/Gconstants.h>


// ____________________________________________________________________________
// i                  WIGNER FUNCTION ERROR HANDLING
// ____________________________________________________________________________


	// Input		int  : Error Flag
	// Output		none : Error Message Output
	// Output		none : Stops Execution & Error Message Output

void          Wigner_error    (int error);
void volatile Wigner_fatality (int error);

// ______________________________________________________________________
// A         REDUCED WIGNER ROTATION MATRIX ELEMENTS: RANK 0 
// ______________________________________________________________________


MSVCDLL double d0();

	// Input		none  :
	// Output		r     : rank 0 reduced Wigner rotation
	//				matrix element
 

// ______________________________________________________________________
// B        REDUCED WIGNER ROTATION MATRIX ELEMENTS: RANK 1/2 
// ______________________________________________________________________


MSVCDLL double d1half(int n, double beta );

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 1/2, m = 1/2 reduced Wigner
	//			        rotation matrix element
	// Note 		      : Euler angle input in degrees
 

MSVCDLL double dm1half(int n, double beta );

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 1/2, m = -1/2 reduced Wigner
	//				rotation matrix element
	// Note 		      : Euler angle input in degrees
 

MSVCDLL double d1half(int m, int n, double beta );

	// Input		m     : momentum index
	// 			n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 1/2 reduced Wigner
	//			        rotation matrix elements
	// Note 		      : Euler angle input in degrees
 

// ______________________________________________________________________
// C         REDUCED WIGNER ROTATION MATRIX ELEMENTS: RANK 1 
// ______________________________________________________________________


MSVCDLL double d11(int n, double beta );

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 1, m = 1 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 

MSVCDLL double d10(int n, double beta );

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 1, m = 0 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 

MSVCDLL double d1m1(int n, double beta);

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 1, m = -1 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 

MSVCDLL double d1(int m, int n, double beta);

	// Input		beta  : Euler angle
	// 			m     : momentum index
	// 			n     : momentum index
	// Output		r     : rank 1 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 

// ______________________________________________________________________
// D         REDUCED WIGNER ROTATION MATRIX ELEMENTS: RANK 2 
// ______________________________________________________________________


MSVCDLL double d22(int n, double beta);

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 2, m = 2 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 

MSVCDLL double d21(int n, double beta );

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 2, m = 1 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 

MSVCDLL double d20(int n, double beta );

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 2, m = 0 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 

MSVCDLL double d2m1(int n, double beta );

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 2, m = -1 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 

MSVCDLL double d2m2(int n, double beta );

	// Input		n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : rank 2, m = -2 reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
 

MSVCDLL double d2(int m, int n, double beta );

	// Input		beta  : Euler angle
	// 			m     : momentum index
	// 			n     : momentum index
	// Output		r     : rank 2 reduced Wigner rotation
	//				matrix element
 

// ______________________________________________________________________
// E         REDUCED WIGNER ROTATION MATRIX ELEMENTS: RANK J 
// ______________________________________________________________________


MSVCDLL double fact(int a);

	// Input		a : an integer
	// Output		r : a double value of a!
	// Note			  : return is double precision to avoid
	//			    overflow at about 13! with int


MSVCDLL double dJ_int(int J, int m, int n, double beta );

	// Input		J     : rank
	// 			m     : momentum index
	// 			n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
	// Note 		      : Only Integral (non-negative) J
 

MSVCDLL double dJ_half_int(int J, int m, int n, double beta );

	// Input		J     : rank
	// 			m     : momentum index
	// 			n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
	// Note 		      : Here J implies 1/2 units
	//				J=1->1/2; m,n=1->1/2; m,n=-1->-1/2
	//				J=2->3/2; m,n=2->3/2; m,n=1->1/2
	//				          m,n=-1->-1/2; m,n=-2->-3/2
	//				In general, J=a->J=a-1/2 and m,n then
	//				span [a-1/2, -(a-1/2)] incremented by 1
 
MSVCDLL double dJ(int J, int m, int n, double beta );

	// Input		J     : rank
	// 			m     : momentum index
	// 			n     : momentum index
	// 			beta  : Euler angle
	// Output		r     : reduced Wigner rotation
	//				matrix element
	// Note 		      : Euler angle input in degrees
	// Note 		      : Negative J implies 1/2 units
	//				J=-1->1/2; m,n=1->1/2; m,n=-1->-1/2
	//				J=-2->3/2; m,n=2->3/2; m,n=1->1/2
	//				           m,n=-1->-1/2; m,n=-2->-3/2
	//				In general, J=-a->J=a-1/2 and m,n then
	//				span [a-1/2, -(a-1/2)] incremented by 1
 

// ______________________________________________________________________
// F                   REDUCED WIGNER ROTATION MATRIX
// ______________________________________________________________________


MSVCDLL matrix dJ(int J, double beta);

	// Input		J     : rank
	// 			beta  : Euler angle
	// Output		mx    : rank J reduced Wigner rotation matrix
	// Note 		      : Euler angle beta input in degrees


// ______________________________________________________________________
// G                   WIGNER ROTATION MATRIX ELEMENTS
// ______________________________________________________________________


MSVCDLL complex DJ(int J, int m, int n, double alpha, double beta, double gamma);

	// Input		J     : rank
	// 			m     : momentum index
	// 			n     : momentum index
	// 			alpha : Euler angle
	// 			beta  : Euler angle
	// 			gamma : Euler angle
	// Output		z     : rank 2 Wigner rotation matrix element
	// Note 		      : Euler angles input in degrees
	// Note 		      : Negative J implies 1/2 units
	//				J=-1->1/2; m,n=1->1/2; m,n=-1->-1/2
	//				J=-2->3/2; m,n=2->3/2; m,n=1->1/2
	//				           m,n=-1->-1/2; m,n=-2->-3/2
	//				In general, J=-a->J=a-1/2 and m,n then
	//				span [a-1/2, -(a-1/2)] incremented by 1


MSVCDLL double D0();

	// Input		none  :
	// Output		r     : rank 0 Wigner rotation matrix element


MSVCDLL complex D1half(int m, int n, double alpha, double beta, double gamma);

	// Input		m     : momentum index
	// 			n     : momentum index
	// 			alpha : Euler angle
	// 			beta  : Euler angle
	// 			gamma : Euler angle
	// Output		z     : rank 1/2 Wigner rotation matrix element
	// Note 		      : Euler angles input in degrees


MSVCDLL complex D1(int m, int n, double alpha, double beta, double gamma);

	// Input		m     : momentum index
	// 			n     : momentum index
	// 			alpha : Euler angle
	// 			beta  : Euler angle
	// 			gamma : Euler angle
	// Output		z     : rank 1 Wigner rotation matrix element
	// Note 		      : Euler angles input in degrees


MSVCDLL complex D2(int m, int n, double alpha, double beta, double gamma);

	// Input		m     : momentum index
	// 			n     : momentum index
	// 			alpha : Euler angle
	// 			beta  : Euler angle
	// 			gamma : Euler angle
	// Output		z     : rank 2 Wigner rotation matrix element
	// Note 		      : Euler angles input in degrees


// ______________________________________________________________________
// H                   WIGNER ROTATION MATRICES
// ______________________________________________________________________

MSVCDLL matrix DJ(int J, double alpha, double beta, double gamma);

	// Input		J     : rank
	// 			alpha : Euler angle
	// 			beta  : Euler angle
	// 			gamma : Euler angle
	// Output		mx    : rank J Wigner rotation matrix
	// Note 		      : Euler angles input in degrees

MSVCDLL matrix DJ(const matrix& dJbeta, int J, double alpha);

        // Input                dj(beta) : rank J Wigner array for beta
        //                      J        : J rank
        //                      alpha    : Euler angle (degrees)
        // Output               mx    : rank J Wigner rotation matrix


#endif							// Wigner.h
