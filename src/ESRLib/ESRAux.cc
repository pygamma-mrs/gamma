/* ESRAux.cc ****************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**   ESR Auxiliary                                  Implementation	**
**                                                                      **
**  Copyright (c) 2000							**
**  Scott A. Smith							**
**  National High Magnetic Field Laboratory                           	**
**  1800 E. Paul Dirac Drive                                          	**
**  Box 4005                                                          	**
**  Tallahassee, FL 32306                                             	**
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/
     
/*************************************************************************
**									**
**  Description								**
**									**
** ESR Auxiliary module has functions of general use in ESR/EPR/EMR     **
** computations. These are function that don't fit into any particular  **
** category.                                                            **
**                                                                      **
*************************************************************************/

#ifndef   ESRAux_cc_			// Is this file already included?
#  define ESRAux_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <ESRLib/ESRAux.h>		// Include the interface

// ____________________________________________________________________________
// A                        LS Coupling Functions
// ____________________________________________________________________________
 
gen_op gJLande(double L, double S, double J)
  { return 1.5 - (L*(L+1) - S*(S+1))/(2*J(J+1)) }


#endif							// ESRAux.cc
