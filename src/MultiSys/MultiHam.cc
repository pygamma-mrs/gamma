/* MultiHam.cc **************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Multiple System Spin Hamiltonians 	Implementation	 	**
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
** involve such a spin system and build up common Hamiltonians.		**
** The functions will mirror those defined in the Hilbert space library **
** (see module HSLib) but will return operators that exist in the       **
** composite Hilbert spin space spanning the components in multi_sys.   **
** The Hilbert space is a direct product space of the systems involved. **
**                                                                      **
*************************************************************************/

#ifndef   MultiHam_cc			// Is the file already included?
#  define MultiHam_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <MultiSys/MultiLib.h>		// Include the header file
#include <MultiSys/MultiSys.h>		// Include multi_sys spin systems
#include <HSLib/HSham.h>		// Inlcude H.S. Hamiltonians
#include <HSLib/GenOp.h>		// Inlcude general operators

// ____________________________________________________________________________
// C         GENERAL PRODUCTION OF MULTISPIN SYSTEM HAMILTONIANS
// ____________________________________________________________________________

gen_op Ho(const      multi_sys &msys) { return (multize(Ho,      msys) ); }
gen_op Hcs(const     multi_sys &msys) { return (multize(Hcs,     msys) ); }
gen_op HJ(const      multi_sys &msys) { return (multize(HJ,      msys) ); }
gen_op Hcs_lab(const multi_sys &msys) { return (multize(Hcs_lab, msys) ); }

#endif							// MultiHam.cc
