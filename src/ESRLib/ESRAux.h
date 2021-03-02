/* ESRAux.h *****************************************************-*-c++-*-
**																		**
** 	                          	  G A M M A								**
**																		**
**	ESR Auxiltonians                               		Interface	    **
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
**  Module ESR Auxiltonians contains functions that return commonly     **
**  used Auxiltonians in ESR/EPR/EMR computations.                      **
**																		**
*************************************************************************/

#ifndef _ESRAux_h_						// Is this file already included?
#define _ESRAux_h_ 1					// If no, then remember it
#  if defined(GAMPRAGMA)							// Using the GNU compiler?
#    pragma interface						// Then this is the interface
#endif

#include <ESRLib/CubicSys.h>			// Include cubic systems

// ____________________________________________________________________________
// A
// ____________________________________________________________________________

gen_op gJLande(double L, double S, double J);

#endif															// ESRAux.h

