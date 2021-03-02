/* HSdecomp.h ***************************************************-*-c++-*-
**                                  					**
**                                   G A M M A				** 
**									**
**      Base Decomposition			Interface		**
**									**
**      Copyright (c) 1991						**
**      Beat H. Meier, Scott A. Smith					**
**      Eidgenoessische Technische Hochschule				**
**      Labor fuer physikalische Chemie					**
**      8092 Zurich / Switzerland					**
**									**
**      $Header: $
**                		              				**
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The GAMMA Platform Provides Functions for Simulation of Magnetic     ** 
** Resonance Experiments and Other Associated Mathematical              **
** Capabilities.  The Set of Functions Herein Provides For The		**
** Decomposition of a General Operator Into its Components Proportional	**
** to an Orthogonal Basis Set.   		                        **
**									**
*************************************************************************/

#ifndef   HSdecomp_h_ 			// Is file already included?
#  define HSdecomp_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)               	// Using the GNU compiler
#    pragma interface               	// this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/SpinSys.h>		// Include base spin systems
#include <HSLib/GenOp.h>		// Include general operators
#define  BD_SMALL 1.0e-15		// Floats < SMALL are zero

MSVCDLL void Prod_base_dec(const spin_sys &sys, const gen_op &Op,double thres=BD_SMALL);

	// Input	sys	: Spin_system
	//		Op	: Operator to be decomposed
        //		thres	: Coefficients below thres (in absolute value) will be set to 0 
	// Output	Decomposition is written to stdout
 
#endif							// HSdecomp.h
