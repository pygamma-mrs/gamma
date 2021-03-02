/* relaxDCSA.h **********************************-*-c++-*-
**							**
**	               G A M M A			**
**							**
**	NMR Library Dipole-CSA Relaxation Functions	**
**							**
**	Interface definition				**
**							**
**	Copyright (c) 1993				**
**	Scott A. Smith					**
**	University of California, Santa Barbara		**
**	Department of Chemistry				**
**	Santa Barbara CA. 93106 USA			**
**							**
**						 	**
**      $Header: $
**							**
*********************************************************/

/*********************************************************
**							**
** 	Description					**
**							**
** Herein are functions dealing with dipole-CSA cross	**
** phenomena.						**
**							**
*********************************************************/

///Chapter Dipole-CSA Relaxation
///Section Overview
///Body The ...
///Section Available Dipole-CSA Relaxation Functions

#ifndef   Relax_DCSA_h_		// Is this file already included?
#  define Relax_DCSA_h_ 1	// If no, then remember it
#  if defined(GAMPRAGMA)	// Using the GNU compiler?
#    pragma interface		// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <BWRRelax/relaxNMR.h>
#include <BWRRelax/relaxDip.h>
#include <BWRRelax/relaxCSA.h>
#include <BWRRelax/relaxRF.h>

// ______________________________________________________________________
// ************* Dipole & CSA Cross Relaxation Superoperators ***********
// ______________________________________________________________________


MSVCDLL super_op RDCX(const sys_dynamic& sys, gen_op& Ho, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : Dipole & CSA cross relaxation
	//			        superoperator
	// Note			      :	Computed in the eigenbasis of Ho


MSVCDLL void RDCX(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                                double* taus, double chi, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Dipole & CSA cross relaxation superop
	// Note			      :	Computed in the eigenbasis of Ho


// ______________________________________________________________________
// ************** Dipole-CSA Relaxation Superoperators ******************
// ______________________________________________________________________


MSVCDLL super_op RDC(const sys_dynamic& sys, gen_op& Ho, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//			level : Relaxation treatment level
	// Output		LOp   : Dipole-CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


MSVCDLL void RDC(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Dipole-CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


// ______________________________________________________________________
// ************** Dipole-CSA Relaxation Superoperators ******************
// ______________________________________________________________________

MSVCDLL super_op RCD(const sys_dynamic& sys, gen_op& Ho, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Dipole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


MSVCDLL void RCD(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                              double* taus, double chi, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Dipole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


// -------------- CSA-Dipole, With An Applied RF-Field ------------------

MSVCDLL super_op RCDrf(const sys_dynamic& sys, gen_op& Heff, double Wrflab, int level=4);

	// Input		sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian (Hertz)
	//			Wrflab: RF-Field frequency, lab frame (Hertz)
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Dipole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff


MSVCDLL void RCDrf(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			Wrflab: Field frequency lab frame (Hz)
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Dipole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff


#endif /* __RELAX_DCSA_H__ */
