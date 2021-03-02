/* MultiHam.h ***************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Multiple System Spin Hamiltonians 	Interface		**
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

#ifndef   GMultiHam_h			// Is the file already included?
#  define GMultiHam_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <MultiSys/MultiSys.h>		// Include multi_sys spin systems
#include <HSLib/GenOp.h>		// Inlcude general operators

// ____________________________________________________________________________
// C         GENERAL PRODUCTION OF MULTISPIN SYSTEM HAMILTONIANS
// ____________________________________________________________________________


MSVCDLL gen_op Ho(const multi_sys &msys);

	// Input		msys    : A multi_spin spin system
	// Output		Ho	: Isotropic Hamiltonian


MSVCDLL gen_op Hcs(const multi_sys &msys);

	// Input		msys    : A multi_spin spin system
	// Output		Hcs	: Isotropic shift Hamiltonian


MSVCDLL gen_op HJ(const multi_sys &msys);

	// Input		msys    : A multi_spin spin system
	// Output		HJ	: Isotropic scalar coupling Hamiltonian


MSVCDLL gen_op Hcs_lab(const multi_sys &msys);

	// Input		msys    : A multi_spin spin system
	// Output		Hcs	: Isotropic shift Hamiltonian (lab)

#endif							// MultiHam.h
