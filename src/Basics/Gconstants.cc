/* Gconstants.cc *************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      GAMMA Constants				Implementation						**
**                                                                      **
**      Copyright (c) 2002												**
**      Scott A. Smith													**		
**      National High Magnetic Field Laboratory                         ** 
**      1800 E. Paul Dirac Drive                                        ** 
**      Tallahassee Florida, 32306-4005                                 ** 
**      email: gamma@magnet.fsu.edu                                     **
**      www: http://gamma.magnet.fsu.edu                                **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  This file contains many constants that are of generaly use to the   **
**  GAMMA simulation platform.  Some of those herein are likely to be   **
**  defined in other parts of the C and/or C++ library.  For example    **
**  the C math.h file will likely have some definitions for PI.         **
**                                                                      **
*************************************************************************/

#ifndef   Gconstants_cc			// Is this file already included?
#  define Gconstants_cc 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <Basics/Gconstants.h>		// Include the header
using std::string;					// Using libstdc++ strings

// ____________________________________________________________________________ 
//                       Various Values Related To PI
// ____________________________________________________________________________

/* Some of these are included in case they are not part of the math library!
   Normally PI and PI2 are defined constants, but we'll check them anyway     */

const double PIx2 = 2.0*PI;			// 2*pi

// ____________________________________________________________________________ 
//                       Various Conversion Factors
// ____________________________________________________________________________

const double DEG2RAD = PI/180.0;			// Degrees -> radians
const double RAD2DEG = 180.0/PI;			// Radians -> degrees

const double HZ2RAD = PIx2;					// Cycles/sec -> rad/sec
const double RAD2HZ = 1.0/PIx2;				// Rad/sec -> cycles/sec

const double HZ2GAUSS  = 7.1447751747e-7;	// h/beta         (G/Hz)
const double GAUSS2HZ  = 1.0/HZ2GAUSS;		// beta/h         (Hz/G)
const double GHZ2GAUSS = 7.1447751747e2;	// (h/beta)*1.e9  (G/GHz)
const double GAUSS2GHZ = 1.0/GHZ2GAUSS;		// (beta/h)*1.e-9 (GHz/G)

// ____________________________________________________________________________ 
//                       Various Spin Isotope Constants
// ____________________________________________________________________________

const double MU_E    = -9.2847701e-24;	// Electron magnetic moment  (J/T)
const double BOHRMAG =  9.2740154e-24;	// Bohr magneton             (J/T)
const double GFREE   =  2.002319304386;	// Free electron g factor    (unitless)

const double GAMMAe  =  1.7608592e+11;	// Free e gyromagnetic ratio (1/T-sec)
const double GAMMA1H =  2.67515255e+8; 	// 1H H20 gyromagnetic ratio (1/T-sec)	
const string DEFISO = "1H";		// Default spin isotope symbol

// ____________________________________________________________________________ 
//                          A Few Other Handy Constants
// ____________________________________________________________________________

const double PLANCK  = 6.6260755e-34;	// Plancks constant (h)         (J-sec)
const double HBAR    = 1.05457266e-34;	// Plancks constant (h/2PI)     (J-sec)

#endif                  		                // Gconstants.cc
