/* PulAuxil.cc **************************************************-*-c++-*-
**			         					**
** 	                        G A M M A				**
**									**
**	Pulse Auxiliary Functions		Implementation		**
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

#ifndef _PulAux_cc_			// Is this file already included?
#define _PulAux_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation 			// This is the implementation
#endif

#include <Pulses/PulAuxil.h>		// Include our own header
#include <Basics/StringCut.h>		// Include string cutting functions
#include <Basics/Gconstants.h>
#include <stdlib.h>

// ____________________________________________________________________________
// A                  PULSE AUXILIARY WAVEFORM FUNCTIONS
// ____________________________________________________________________________


row_vector pulseshift(row_vector& p, row_vector& ptime, const double& FreqOffset)
//created by LK, Dec. 2007

// * FreqOffset is in Hz
// * ptime contains dwell time in seconds. 
//   This is not cumulative, just the dwell time for that step.

// * p is a vector of complex number that describes a pulse profile. 
//   First column is the norm or absolute value,
//   Second column is the angle or phase(in degrees). 

{
	row_vector vector_points(p.size());
	row_vector data1(p.size());
		 
	  // Loop through points
	  for(int i=0; i<p.size(); i++)	
		{
      // Store integer count (1,2,3...p.size)
			vector_points.put(i+1, i);
		}

	// Theta is the argument of the freq. offset	
	row_vector theta1=2*PI*FreqOffset*ptime(0)*vector_points; 
		  
	// Now make it complex
	// Loop through points
	for(int i=0; i<p.size(); i++)
  {
	 		double evalue1 = Re(theta1(i)) ;
			double evalue2 = Im(theta1(i)) ;
			complex a(evalue2,evalue1);

			complex rev=exp(a);
	
			// Store this value
      data1.put(complex(norm(rev),RAD2DEG*phase(rev)), i);
   }

	//data1 contains the exponential evolution factor to shift frequency
	
	//data is the final vector
	row_vector data(p.size());
	
	// Loop through points
	for(int i=0; i<p.size(); i++)
  {
		// Store this value
		data.put(complex(Re(p(i))*Re(data1(i)),Im(p(i))+Im(data1(i))), i);
  }

	return data;
}

// ____________________________________________________________________________
// Z                     PULSE AUXILIARY I/O FUNCTIONS
// ____________________________________________________________________________


        // Input                ostr    : Output stream
        //                      i	: Index to print
        // Output               ostr    : Modified to contain printed 
        //             			  index with next output spot back
	//				  at index start.

std::ostream& printIndx(std::ostream& ostr, int i)
  {
  ostr << i;
             ostr << "\b";
  if(abs(i)>9)    ostr << "\b";
  if(abs(i)>99)   ostr << "\b";
  if(abs(i)>999)  ostr << "\b";
  if(abs(i)>9999) ostr << "\b";
  if(i < 0)       ostr << "\b";
  ostr.flush();
  return ostr;
  }

	// Input 		ostr    : Output stream
	//			Str	: String to print
	//			plen    : String print length
        // Output               ostr    : Modified to contain printed info
        //                                regarding the string

std::ostream& printStr(std::ostream& ostr, const std::string& Str, int plen)
  {
  int len = Str.length();
  ostr << std::string(plen-len, ' ')  << Str;
  return ostr;
  }


	// Input 		ostr    : Output stream
	//			td	: A delay time
        // Output               ostr    : Modified to contain printed info
        //                                regarding the delay length

std::ostream& printTime(std::ostream& ostr, double td)
  {
  if(td>0.1 || !td)    ostr << Gform("%8.3f", td)       << " sec ";
  else if(td > 1.e-4)  ostr << Gform("%8.3f", td*1.e3)  << " msec";
  else if(td > 1.e-7)  ostr << Gform("%8.3f", td*1.e6)  << " usec";
  else if(td > 1.e-10) ostr << Gform("%8.3f", td*1.e9)  << " nsec";
  else if(td > 1.e-13) ostr << Gform("%8.3f", td*1.e12) << " psec";
  else if(td > 1.e-16) ostr << Gform("%8.3f", td*1.e15) << " fsec";
  else                 ostr << td << " sec ";
  return ostr;
  }


	// Input 		ostr    : Output stream
	//			F       : A frequency
        // Output               ostr    : Modified to contain printed info
        //                                regarding the frequency

std::ostream& printHz(std::ostream& ostr, double F)
  {
  if(F<1000 || !F)    ostr << Gform("%8.3f", F) << " Hz";
  else if(F < 100000) ostr << Gform("%8.3f", F*1.e-3) << " KHz";
  else                ostr << Gform("%8.3f", F*1.e-6) << " MHz";
  return ostr;
  }


#endif								// PulAux.cc
