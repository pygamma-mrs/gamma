/* HSdetect.h **************************************************-*-c++-*-
**									**
**                              G A M M A				**
**								 	**
**	Detection, NMR Library                    Interface		**
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
**                                                                      **
** Description                                                          **
**                                                                      **
** The GAMMA Platform Provides Functions for Simulation of Magnetic     ** 
** Resonance Experiments and Other Associated Mathematical              **
** Capabilities.  The Set of Functions Herein Allows for the Generation **
** of Detection Operators for Obtaining Transverse Magnetization in     **
** Spin Hilbert Space. These Functions Provide the Operators for        **
** Detection of xy-plane Magnetization.                                 **
**                                                                      **
*************************************************************************/

#ifndef   HSdetect_h_			// Is file already included?
#  define HSdetect_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
class spin_op;				// Know spin_op is a class
class spin_sys;				// Know spin_sys is a class

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

// sosi - want to pass D as a constant reference below.....
MSVCDLL spin_op GenericD(const spin_sys& sys, spin_op D, double beta=0);


MSVCDLL spin_op detector(const spin_sys &S, double B=0);
MSVCDLL spin_op Mxy(const spin_sys &sys, double beta=0);

	// Output		D   : spin operator I- affecting all
	//			      spins in the spin system and
	//			      rotated to a phase angle beta


MSVCDLL spin_op detector(const spin_sys &S, char *I, double B=0);
MSVCDLL spin_op Mxy(const spin_sys &sys, char *iso, double beta=0);

	// Input  		iso : isotope type
	// Output		D   : spin operator I+ affecting only
	//			      spins of the isotopes specified
	//			      and rotated to a phase angle beta


MSVCDLL spin_op detector(const spin_sys &S, int N, double B=0);
MSVCDLL spin_op Mxy(const spin_sys &sys, int spin, double beta=0);

	//			spin: spin number
	// Output		D   : spin operator I- affecting the
	//			      specified spin in the spin system
	//			      & rotated to a phase angle beta


MSVCDLL spin_op detector_sp(const spin_sys &S, double B=0);
MSVCDLL spin_op Mxy_sp(const spin_sys &sys, double beta=0);

	// Output		D   : spin operator I- affecting the
	//			      spins in the spin system with
	//			      spin flags set to TRUE at function
	//			      call & rotated to a phase angle beta


#endif 							// HSdetect.h
