/* ESRHam.h *****************************************************-*-c++-*-
**																		**
** 	                          	  G A M M A								**
**																		**
**	ESR Hamiltonians                               		Interface	    **
**								 										**
**	Copyright (c) 2000					 								**
**	Scott A. Smith				 										**
**  National High Magnetic Field Laboratory                           	**
**  1800 E. Paul Dirac Drive                                          	**
**  Box 4005                                                          	**
**  Tallahassee, FL 32306                                             	**
**						 												**
**      $Header:$
**								 										**
*************************************************************************/

/*************************************************************************
**							 											**
**  Description						 									**
**							 											**
**  Module ESR Hamiltonians contains functions that return commonly     **
**  used Hamiltonians in ESR/EPR/EMR computations.                      **
**																		**
*************************************************************************/

#ifndef _ESRHam_h_						// Is this file already included?
#define _ESRHam_h_ 1					// If no, then remember it
#  if defined(GAMPRAGMA)							// Using the GNU compiler?
#    pragma interface						// Then this is the interface
#endif

#include <HSLib/GenOp.h>				// Include general operator
#include <ESRLib/CubicSys.h>			// Include cubic systems


// ____________________________________________________________________________
// A                          Zeeman Hamiltonians
// ____________________________________________________________________________


// ____________________________________________________________________________
// B
// ____________________________________________________________________________

gen_op O40(const CubicSys& CuSys);
gen_op O44(const CubicSys& CuSys);
gen_op O60(const CubicSys& CuSys);
gen_op O64(const CubicSys& CuSys);

#endif															// ESRHam.h

