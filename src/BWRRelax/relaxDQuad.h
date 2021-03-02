/* relaxDQuad.h *************************************************-*-c++-*-
**									**
** NMR Dipole-Quadrupole X-Correlation               Interface 		**
**                                                                      **
**      Copyright (c) 1997                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** These are functions dealing with dipole-quadrupole cross correlation **
** effects in NMR as according to WBR relaxation theory.                **
**                                                                      **
*************************************************************************/

///Chapter Dipole-Quadrupole Cross Correlation
///Section Overview
///Body The ...
///Section Available Dipole-Quadrupole Cross Correlation Functions

#ifndef _relax_DQuad_h_		// Is this file already included?
#define _relax_DQuad_h_ 1	// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface		// This is the interface
#endif

#include <BWRRelax/relaxNMR.h>
#include <BWRRelax/relaxDip.h>	// Include dipolar relaxation
#include <BWRRelax/relaxQuad.h>	// Include quadrupolar relaxation 
#include <BWRRelax/relaxRF.h>	// Include RF field effects

// ______________________________________________________________________
// ******** Dipole & Quadrupolar Cross Correlation Superoperators *******
// ______________________________________________________________________


super_op RDQX(sys_dynamic& sys, gen_op& Ho, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : Dipole & Quad cross relaxation
	//			        superoperator
	// Note			      :	Computed in the eigenbasis of Ho


void RDQX(super_op& LOp, sys_dynamic& sys, gen_op& Ho, double*w,
                                double* taus, double chi, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Dipole & Quad cross relaxation superop
	// Note			      :	Computed in the eigenbasis of Ho


// ______________________________________________________________________
// ************** Dipole-Quad Relaxation Superoperators ******************
// ______________________________________________________________________


super_op RDQ(sys_dynamic& sys, gen_op& Ho, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//			level : Relaxation treatment level
	// Output		LOp   : Dipole-Quad relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


void RDQ(super_op& LOp, sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Dipole-Quad relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


// ______________________________________________________________________
// ************** Dipole-Quad Relaxation Superoperators ******************
// ______________________________________________________________________

super_op RQD(sys_dynamic& sys, gen_op& Ho, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : Quad-Dipole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


void RQD(super_op& LOp, sys_dynamic& sys, gen_op& Ho, double*w,
                              double* taus, double chi, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Quad-Dipole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


// -------------- Quad-Dipole, With An Applied RF-Field ------------------

super_op RQDrf(sys_dynamic& sys, gen_op& Heff, double Wrflab, int level=4);

	// Input		sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian (Hertz)
	//			Wrflab: RF-Field frequency, lab frame (Hertz)
	//			level : Relaxation treatment level
	// Output		LOp   : Quad-Dipole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff


void RQDrf(super_op& LOp, sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			Wrflab: Field frequency lab frame (Hz)
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Quad-Dipole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff


#endif 						// relaxDipQ.h
