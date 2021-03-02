/* Gconstants.h *************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      GAMMA Constants				Interface		**
**                                                                      **
**      Copyright (c) 1990, 1991, 1992, 1993, 1994                      **
**      Tilo Levante, Scott A. Smith 					**
**      Eidgenoessische Technische Hochschule                           **
**      Labor fur physikalische Chemie                                  **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**      email: gamma@magnet.fsu.edu					**
**      www: http://gamma.magnet.fsu.edu				**
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  This file contains many constants that are of generaly use to the   **
**  GAMMA simulation platform.  Many common constants are likely to be  **
**  defined in other parts of the C and/or C++ library.  For example    **
**  the C math.h file will likely have some definitions for PI.  I see  **
**  one in /lib/gcc-lib/arch/ver/include/math which reads               **
**                                                                      **
**           #define M_PI            3.14159265358979323846             **
**                                                                      **
**  A particular problem is that PI &/or PI2 may be defined elsewhere.	**
**  There is a check herein for those two cases but the real trouble    **
**  will be if PI2 is defined as PI/2 when GAMMA wants it PI*2!         **
**  That's exactly what happened in a recent compilation under Linux!   **
**  Now, PIx2 is defined herein.                                        **
**                                                                      **
*************************************************************************/

#ifndef   Gconstants_h			// Is this file already included?
#  define Gconstants_h 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <string>			// Know libstdc++ strings
#include <cmath>			// Has some constants we want
#include <GamGen.h>			// Know MSVCDLL (__declspec)

// ____________________________________________________________________________ 
//            TypeDefs To Keep Older Class Names Relatively Viable
// ____________________________________________________________________________


//extern const MSVCDLL double roundoff;	// Used for equality testing

// ____________________________________________________________________________ 
//                   Various Values Related To PI
// ____________________________________________________________________________

/* Some of these are included in case they are not part of the math library!
   Normally PI and PI2 are defined constants, but we'll check them anyway     */
   
#ifndef PI				// Define PI once
#define PI 3.14159265358979323846
#endif
 
#ifndef PI2                             // Define 2PI once so we'll
#define PI2 6.28318530717958647692      // always have it around
#endif
 
extern const MSVCDLL double PIx2;	// 2*pi

// ____________________________________________________________________________ 
//                       Various Conversion Factors
// ____________________________________________________________________________

extern const MSVCDLL double DEG2RAD;            // Degrees -> radians
extern const MSVCDLL double RAD2DEG;            // Radians -> degrees

extern const MSVCDLL double HZ2RAD;             // Cycles/sec -> rad/sec
extern const MSVCDLL double RAD2HZ;             // Rad/sec -> cycles/sec

extern const MSVCDLL double HZ2GAUSS;		// cycles/sec -> Gauss
extern const MSVCDLL double GAUSS2HZ;		// Gauss -> cycles/sec
extern const MSVCDLL double GHZ2GAUSS;		// h/beta*1.e9
extern const MSVCDLL double GAUSS2GHZ;		// Gauss -> cycles/sec * 1.e-9

// ____________________________________________________________________________ 
//                       Various Spin Isotope Constants
// ____________________________________________________________________________

extern const MSVCDLL double      MU_E; 		// Electron magnetic moment  (J/T)
extern const MSVCDLL double      BOHRMAG;	// Bohr magneton             (J/T)
extern const MSVCDLL double      GFREE;		// Free electron g factor
extern const MSVCDLL double      GAMMAe;	// Free e gyromagnetic ratio (1/T-sec)
extern const MSVCDLL double      GAMMA1H;	// Proton gyromagnetic ratio (1/T-sec)	
extern const MSVCDLL std::string DEFISO;	// Default spin isotope symbol

// ____________________________________________________________________________
//                          A Few Other Handy Constants
// ____________________________________________________________________________

extern const MSVCDLL double PLANCK;		// Plancks constant (h)         (J-sec)
extern const MSVCDLL double HBAR;   		// Plancks constant (h/2PI)     (J-sec)

#endif						// Gconstants.h
