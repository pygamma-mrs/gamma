/* PulAuxil.h ***************************************************-*-c++-*-
**			         					**
** 	                        G A M M A				**
**									**
**	Pulse Auxiliary Functions			Interface	**
**									**
**      Copyright (c) 1998                                              **
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**                                                                      **
** This file contains auxiliary functions that support GAMMA pulse 	**
** related classes.							**
**                                                                      **
*************************************************************************/

#ifndef   GPulAux_h_			// Is this file already included?
#  define GPulAux_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/row_vector.h>
#include <iostream>			// Include filestreams
#include <string>			// Include stdlibc++ strings

// ____________________________________________________________________________
// A                  PULSE AUXILIARY WAVEFORM FUNCTIONS
// ____________________________________________________________________________


MSVCDLL row_vector pulseshift(row_vector& p, row_vector& ptime, const double& FreqOffset);


// ____________________________________________________________________________
// Z                     PULSE AUXILIARY I/O FUNCTIONS
// ____________________________________________________________________________

 
MSVCDLL std::ostream& printIndx(std::ostream& ostr, int i);
 
        // Input                ostr    : Output stream 
        //                      i       : Index to print 
        // Output               ostr    : Modified to contain printed 
        //                                index with next output spot back
        //                                at index start.


MSVCDLL std::ostream& printStr(std::ostream& ostr, const std::string& Str, int plen=5);
 
        // Input                ostr    : Output stream
        //                      Str     : String to print
        //                      plen    : String print length
        // Output               ostr    : Modified to contain printed info
        //                                regarding the string
 


MSVCDLL std::ostream& printTime(std::ostream& ostr, double td);

	// Input 		ostr    : Output stream
	//			td	: A delay time
        // Output               ostr    : Modified to contain printed info
        //                                regarding the delay length


MSVCDLL std::ostream& printHz(std::ostream& ostr, double F);

	// Input 		ostr    : Output stream
	//			F       : A frequency
        // Output               ostr    : Modified to contain printed info
        //                                regarding the frequency

#endif						// PulAuxil.h
