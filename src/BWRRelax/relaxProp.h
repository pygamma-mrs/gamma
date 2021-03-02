/* relaxProp.h **********************************-*-c++-*-
**							**
**	               G A M M A			**
**							**
**	Relaxation Propagator Functions			**
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
** Herein are functions dealing with some of the more	**
** common superoperator propagators.			**
**							**
*********************************************************/

///Chapter Superoperator Propagators
///Section Overview
///Body The ...
///Section Available Superoperator Propagator Functions

#ifndef   Relax_prop_h_		// Is this file already included?
#  define Relax_prop_h_ 1	// If no, then remember it
#  if defined(GAMPRAGMA)	// Using the GNU compiler?
#    pragma interface		// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <BWRRelax/relaxNMR.h>

// sosi: this is now in LSLib module in R_prop
//void set_trace(gen_op& sigma, double tr);

// ______________________________________________________________________
// *************** Specific Density Matrix Superoperators ***************
// ______________________________________________________________________


MSVCDLL super_op LOp_sigma(gen_op& sigma);

// Input		sigma   : Density matrix
// Output		LOp     : Superoperator which fullfills
//	                              sigma = LOp * sigma0

//				  for any valid sigma0
//
// Note				: LOp constructed in default Liouville
//				  space and Hilbert space of sigma

//  <a,aa|LOp|b,bb> = del    * <a,aa|sigma|1> = del    * <a|sigma|aa>
//			 b,bb                      b,bb



//super_op R_prop(super_op& eLt, gen_op sigmaeq);

// Input		eLt	: Exponential Liouvillian to relax
//				  the density matrix
//			sigmaeq : Density matrix at equilibrium
//				  (or steady state density matrix)
//			Ho      : Operator in Hilbert space
// Output		LOp     : Superoperator which fullfills
//
//	   eLt * {sigma-sigmaeq} + sigmaeq  = LOp * sigma
//
// Note				: LOp is constructed in the current basis
//				  of sigmaeq


//super_op R_prop(super_op& L, gen_op& sigmaeq, double t);

// Input		eLt	: Liouvillian for evolution
//			sigmaeq : Density matrix at equilibrium
//				  (or steady state density matrix)
//			t       : Evolution time
// Output		LOp     : Superoperator which fullfills
//
//	   eLt * {sigma-sigmaeq} + sigmaeq  = LOp * sigma


//super_op gammaAx(super_op& eLt, gen_op& sigmaeq);

// Input		eLt	: Exponential Liouvillian
//			sigmaeq : Density matrix at equilibrium
//			Ho      : Operator in Hilbert space
// Output		LOp      : Superoperator fullfilling the
//				  following relationship for any
//				  valid density matrix sigma 
//
//	   eLt * sigmaeq = A * sigma   for any sigma
//
// Note				: A is constructed in the default
//				  Liouville space and in the Hilbert
//				  space of the input operator Ho


//super_op gammaBx(gen_op& sigmaeq);

// Input		sigmaeq : Density matrix at equilibrium
// Output		B       : Superoperator fullfilling the
//				  following relationship for any
//				  valid density matrix sigma 
//
//	     sigmaeq = B * sigma   for any sigma
//
// Note				: B is constructed in the default
//				  Liouville space and in the Hilbert
//				  space of Ho


#endif                                                 // relaxProp.h
