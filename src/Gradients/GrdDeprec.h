/* GrdDeprec.h *************************************************-*-c++-***
**									**
**	                         G A M M A 				**
**						 			**
**	Field Gradients Deprecated Functions               Interface	**
**						 			**
**      Copyright (c) 1996                                              **
**      Dr. Scott A. Smith             					**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**							 		**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  This GAMMA module provides functions for simulation of field        **
**  gradients in NMR spectroscopy.                                      **
**                                                                      **
*************************************************************************/

#ifndef   GrdDep_h_			// Is file already included?
#  define GrdDep_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/GenOp.h>		// Know general operators

// ____________________________________________________________________________
// A                        Time Evolution Functions
// ____________________________________________________________________________

MSVCDLL void evolve(int N, gen_op& s0, gen_op* Hs, double t, gen_op* sigs);// DEPRECATED!
MSVCDLL void evolve(int N, gen_op& sigma0, gen_op* Us, gen_op* sigmas);    // DEPRECATED!

MSVCDLL void evolve(int N, gen_op* s0, gen_op& H, double t, gen_op *sigs); // DEPRECATED!
MSVCDLL void evolve(int N, gen_op* s0, gen_op& U, gen_op* sigs);           // DEPRECATED!

MSVCDLL void evolve(int N, gen_op* s0, gen_op* Hs, double t, gen_op* sigs);// DEPRECATED!
MSVCDLL void evolve(int N, gen_op* s0, gen_op* U, gen_op* sigs);           // DEPRECATED!

#endif                                                          // GradDeprec.h
