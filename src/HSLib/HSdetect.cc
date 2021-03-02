/* HSdetect.cc **************************************************-*-c++-*-
**									**
**                              G A M M A				**
**								 	**
**	Detection, NMR Library                    Implementation	**
**								 	**
**	Copyright (c) 1990			 			**
**	Scott Smith				 			**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fuer physikalische Chemie		 			**
**	8092 Zurich / Switzerland		 			**
**						 			**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**							 		**
** The GAMMA Platform Provides Functions for Simulation of Magnetic     **
** Resonance Experiments and Other Associated Mathematical 		**
** Capabilities.  The Set of Functions Herein Allows for the Generation	**
** of Detection Operators for Obtaining Transverse Magnetization in 	**
** Spin Hilbert Space. These Functions Provide the Operators for	**
** Detection of xy-plane Magnetization.					**
**								 	**
*************************************************************************/

#ifndef   HSdetect_cc_			// Is file already included?
#  define HSdetect_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// this is the implementation
#  endif

#include <HSLib/HSdetect.h>		// Include header info
#include <HSLib/SpinOpCmp.h>		// Include composite spin operators
#include <HSLib/SpinOpRot.h>		// Include spin rotation operators

// ____________________________________________________________________________
// A			     DETECTION OPERATORS
// ____________________________________________________________________________

	// Input		sys : A spin system
	//			beta: A phase angle (degrees)
	// Output		D   : A Hilbert space detection operator
	// Note			    : Selectivity is determined by
	//			      additional function arguments
	// Note			    : Mxy implies an operator for bulk
	//			      magnetization in the xy-plane

spin_op GenericD(const spin_sys& sys, spin_op D, double beta)

  {
  spin_op RZ = Rz(sys,beta);
  return  RZ*D*RZ.adjoint();
  }


spin_op detector(const spin_sys &S, double B) { return Mxy(S,B); }
spin_op Mxy(const spin_sys &sys, double beta)

	// Output		D   : spin operator I- affecting all
	//			      spins in the spin system and
	//			      rotated to a phase angle beta

  {
  if(!beta) return Fm(sys);
  return GenericD(sys,Fm(sys),beta);
  }


spin_op detector(const spin_sys &S, char *I, double B) { return Mxy(S,I,B); }
spin_op Mxy(const spin_sys &sys, char *iso, double beta)

	// Input  		iso : isotope type
	// Output		D   : spin operator I+ affecting only
	//			      spins of the isotopes specified
	//			      and rotated to a phase angle beta

  {
  if(!beta) return Fm(sys,iso);
  return GenericD(sys,Fm(sys,iso),beta);
  }


spin_op detector(const spin_sys &S, int N, double B) { return Mxy(S,N,B); }
spin_op Mxy(const spin_sys &sys, int spin, double beta)

	//			spin: spin number
	// Output		D   : spin operator I- affecting the
	//			      specified spin in the spin system
	//			      & rotated to a phase angle beta

  {
  if(!beta) return Fm(sys,spin);
  return GenericD(sys,Fm(sys,spin),beta);
  }


spin_op detector_sp(const spin_sys &S, double B) { return Mxy_sp(S,B); }
spin_op Mxy_sp(const spin_sys &sys, double beta)

	// Output		D   : spin operator I- affecting the
	//			      spins in the spin system with
	//			      spin flags set to TRUE at function
	//			      call & rotated to a phase angle beta

  {
  if(!beta) return Fm_sp(sys);
  return GenericD(sys,Fm_sp(sys,beta));
  }

#endif 						// HSdetect.cc
