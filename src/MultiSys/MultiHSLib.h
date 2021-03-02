/* MultiHSLib.h ***************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Multiple System Hilber Space Functions 		Interface	**
**                                                                      **
**      Copyright (c) 1995                                              **
**      Nikolai Skrynnikov                                              **
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: 
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** This file of functions support multi_sys, the GAMMA class for        **
** handling mulitple spin systems.  The routines herein generally       **
** involve such a spin system and build up common spin operators. 	**
** The functions will mirror those defined in the Hilbert space library	**
** (see module HSLib) but will return operators that exist in the 	**
** composite Hilbert spin space spanning the components in multi_sys.	**
** The Hilbert space is a direct product space of the systems involved. **
**                                                                      **
*************************************************************************/

#ifndef   MultiHSLib_h			// Is the file already included?
#  define MultiHSLib_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <MultiSys/MultiSys.h>		// Knowledge of multiple systems
#include <HSLib/GenOp.h>		// Knowledge of operators




// ____________________________________________________________________________
//         GENERAL PRODUCTION OF MULTISPIN SYSTEM DENSITY OPERATORS
// ____________________________________________________________________________
 
	       
MSVCDLL gen_op sigma_eq(const multi_sys &msys);

        // Input                msys    : A multi_spin spin system
        // Output               sigma0  : Density operator for the
	//				  system msys at equilibrium

#endif						// MultiHSLib.h
