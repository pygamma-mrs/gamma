/* SphHarmic.h **************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Spherical Harmonic Functions              Interface		**
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
** This module provides spherical harmonics to GAMMA.                   **
**                                                                      **
*************************************************************************/

#ifndef   Gsphharm_h_			// Is file already included?
#  define Gsphharm_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/complex.h>
#include <Matrix/matrix.h>
#include <iostream>
#include <Basics/Gconstants.h>


// ____________________________________________________________________________
// i                  SPHERICAL HARMONICS ERROR HANDLING
// ____________________________________________________________________________


void Y_error(int eidx, int noret=0);

        // Input                eidx    : Error index
	// Output		none	: Error Message Output


void Y_error(int eidx, const std::string& pname, int noret=0);
 
        // Input                eidx    : Error index
        //                      pname   : Additional error message 
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error Message Output


void volatile Y_fatality(int error);

	// Input		none :
	// Output		none : Stops Execution & Error Message Output

// ____________________________________________________________________________
//                 SPHERICAL HARMONICS EMPLOYING RADIANS
// ____________________________________________________________________________
 

MSVCDLL double Y00rad();

	// Input		none : unnecessary
	// Output		z    : rank zero normalized spherical harmonic


MSVCDLL double Y10rad(double theta);

	// Input		theta : spherical angle
	// Output		z     : rank one normalized spherical harmonic
	// Note			      : angle theta input in radians


MSVCDLL complex Y11rad(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 1 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians


MSVCDLL complex Y1m1rad(double  theta, double  phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 1 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians


MSVCDLL double Y20rad(double theta);

	// Input		theta : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angle theta input in radians


MSVCDLL complex Y21rad(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians


MSVCDLL complex Y2m1rad(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians


MSVCDLL complex Y22rad(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians


MSVCDLL complex Y2m2rad(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians


MSVCDLL double Y30rad(double theta);

	// Input		theta : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angle theta input in radians


MSVCDLL complex Y31rad(double  theta, double  phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians
 

MSVCDLL complex Y3m1rad(double theta, double  phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians


MSVCDLL complex Y32rad(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians


MSVCDLL complex Y3m2rad(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians


MSVCDLL complex Y33rad(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians


MSVCDLL complex Y3m3rad(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in radians


MSVCDLL complex Ylmrad(int l, int m, double  theta, double  phi);

	// Input		l     : angular momentum
	// 			m     : angular momentum component
	// 			theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank l normalized spherical harmonic
	// Note			      : angles theta and phi input in radians


// ____________________________________________________________________________
//                 SPHERICAL HARMONICS EMPLOYING DEGREES
// ____________________________________________________________________________


MSVCDLL double Y00();

	// Input		none : unnecessary
	// Output		z    : rank zero normalized spherical harmonic


MSVCDLL double Y10(double theta);

	// Input		theta : spherical angle
	// Output		z     : rank one normalized spherical harmonic
	// Note			      : angle theta input in degrees


MSVCDLL complex Y11(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 1 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees


MSVCDLL complex Y1m1(double  theta, double  phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 1 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees


MSVCDLL double Y20(double  theta);

	// Input		theta : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angle theta input in degrees


MSVCDLL complex Y21(double  theta, double  phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees


MSVCDLL complex Y2m1(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees


MSVCDLL complex Y22(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees


MSVCDLL complex Y2m2(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 2 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees


MSVCDLL double Y30(double theta);

	// Input		theta : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angle theta input in radians


MSVCDLL complex Y31(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees
 

MSVCDLL complex Y3m1(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees


MSVCDLL complex Y32(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees


MSVCDLL complex Y3m2(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees


MSVCDLL complex Y33(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees


MSVCDLL complex Y3m3(double theta, double phi);

	// Input		theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank 3 normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees


MSVCDLL complex Ylm(int l, int m, double theta, double phi);

	// Input		l     : angular momentum
	// 			m     : angular momentum component
	// 			theta : spherical angle
	// 			phi   : spherical angle
	// Output		z     : rank l normalized spherical harmonic
	// Note			      : angles theta and phi input in degrees

#endif						// SphHarmic.h
