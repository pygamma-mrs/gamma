/* MultiWBR.cc **************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Multiple System BWR Relaxation 			Implementation 	**
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
** involve such a spin system & build up BWR relaxation superoperators. **
** The functions will mirror those defined in the Hilbert space library **
** (see module HSLib) but will return operators that exist in the       **
** composite Hilbert spin space spanning the components in multi_sys.   **
** The Hilbert space is a direct product space of the systems involved. **
**                                                                      **
*************************************************************************/

#ifndef _MultiWBR_cc			// Is the file already included?
#define _MultiWBR_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// Then this is the implementation
#endif

#if defined(_MSC_VER)			// If we are using MSVC++
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif

#include <MultiSys/MultiWBR.h>		// Include our header file
#include <HSLib/SpinSys.h>		// Include base spin ssytems
#include <MultiSys/MultiLib.h>		// Include the header file
#include <MultiSys/MultiSys.h>		// Include multi_sys spin systems
#include <HSLib/GenOp.h>		// Inlcude general operators
#include <BWRRelax/relaxDip.h>               // Include dipolar relaxation
#include <BWRRelax/relaxCSA.h>               // Include CSA relaxation
#include <BWRRelax/relaxQuad.h>              // Include quadrupole relaxation
#include <BWRRelax/relaxQCSA.h>              // Include CSA-Quad X-correlation

// ____________________________________________________________________________
// G    PRODUCTION OF MULTISPIN SYSTEM PULSE RELAXATION SUPEROPERATORS
// ____________________________________________________________________________

 
/*	   Input		msys    : A multi_spin spin system
	  			H       : Isotropic Hamiltonian
	  			type    : Computation type
	  			level   : Computation level
	   Output		R	: Relaxation superoperator
	   Note				: The multize function called must
	  				  be defined for this to work

    Function					Output
    ---------	        -------------------------------------------------------
      RDD		Dipole-Dipole relaxation superoperator
      RCC		CSA-CSA relaxation superoperator
      RQQ		Quadrupole-Quadrupole relaxation superoperator
      RCQ,RQC		Quadrupole-CSA cross-correlation superoperator        */

super_op RQQ(const multi_sys& msys, gen_op& H, int type, int level)
  { return (multize(RQQ, H, type, level, msys)); }

super_op RCC(const multi_sys& msys, gen_op& H, int type, int level)
  { return (multize(RCC, H, type, level, msys)); }

super_op RDD(const multi_sys& msys, gen_op& H, int type, int level)
  { return (multize(RDD, H, type, level, msys)); }

//super_op RDQ(const multi_sys& msys, gen_op& H, int type, int level)
//  { return (multize(RDQ, H, type, level, msys)); }

//super_op RQD(const multi_sys& msys, gen_op& H, int type, int level)
//  { return (multize(RQD, H, type, level, msys)); }

//super_op RCD(const multi_sys& msys, gen_op& H, int type, int level)
//  { return (multize(RCD, H, level, msys)); }

//super_op RDC(const multi_sys& msys, gen_op& H, int type, int level)
//  { return (multize(RDC, H, level, msys)); }
 

super_op RCQ(const multi_sys& msys, gen_op& H, int level)
  { return (multize(RCQ, H, level, msys)); }
 
super_op RQC(const multi_sys& msys, gen_op& H, int level)
  { return (multize(RQC, H, level, msys)); }

#endif							// MultiWBR.cc
