/* relax_QCSA.h *************************************************-*-c++-*-
**									**
**	                            G A M M A				**
**									**
**	NMR Quad-CSA X-Correlation 	        Interface definition	**
**									**
**	Copyright (c) 1997						**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**								 	**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** Description								**
**									**
** Herein are functions dealing with quadrupole-CSA cross correlation	**
** effects in NMR as according to WBR relaxation theory.		**
**									**
*************************************************************************/

///Chapter Quadrupole-CSA Relaxation
///Section Overview
///Body The ...
///Section Available Quadrupole-CSA Relaxation Functions

#ifndef   Relax_QCSA_h_			// Is this file already included?
#  define Relax_QCSA_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <BWRRelax/relaxNMR.h>		// Include base WBR functions
#include <BWRRelax/relaxQuad.h>		// Include quadrupolar relaxation
#include <BWRRelax/relaxCSA.h>		// Include CSA relaxation
#include <BWRRelax/relaxRF.h>		// Include rf-field functions

// ______________________________________________________________________
// *********** Quadrupole & CSA Cross Relaxation Superoperators *********
// ______________________________________________________________________


MSVCDLL super_op RQCX(const sys_dynamic& sys, gen_op& Ho, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : Quadrupole & CSA cross relaxation
	//			        superoperator
	// Note			      :	Computed in the eigenbasis of Ho


MSVCDLL void RQCX(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                                double* taus, double chi, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Quad. & CSA X-relaxation superop
	// Note			      :	Computed in the eigenbasis of Ho


// ______________________________________________________________________
// ************ Quadrupole-CSA Relaxation Superoperators ****************
// ______________________________________________________________________


MSVCDLL super_op RQC(const sys_dynamic& sys, gen_op& Ho, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//			level : Relaxation treatment level
	// Output		LOp   : Quad. & CSA X-relaxation superop
	// Note			      :	Computed in the eigenbasis of Ho


MSVCDLL void RQC(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Quad. & CSA X-relaxation superop
	// Note			      :	Computed in the eigenbasis of Ho


// ______________________________________________________________________
// ************ Quadrupole-CSA Relaxation Superoperators ****************
// ______________________________________________________________________

MSVCDLL super_op RCQ(const sys_dynamic& sys, gen_op& Ho, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Quad X-relaxation superop
	// Note			      :	Computed in the eigenbasis of Ho


MSVCDLL void RCQ(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                              double* taus, double chi, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Quad X-relaxation superop
	// Note			      :	Computed in the eigenbasis of Ho


// ------------ CSA-Quadrupole, With An Applied RF-Field ----------------

MSVCDLL super_op RCQrf(const sys_dynamic& sys, gen_op& Heff, double Wrflab, int level=4);

	// Input		sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian (Hertz)
	//			Wrflab: RF-Field frequency, lab frame (Hertz)
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Quad X-relaxation superop
	// Note			      :	Computed in the eigenbasis of Heff


MSVCDLL void RCQrf(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			Wrflab: Field frequency lab frame (Hz)
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Quad X-relaxation superop
	// Note			      :	Computed in the eigenbasis of Heff


#endif 						// __RELAX_QCSA_H__ 
