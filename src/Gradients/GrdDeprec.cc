/* GrdDeprec.cc ************************************************-*-c++-***
**									**
**	                         G A M M A 				**
**						 			**
**	Field Gradients Deprecated Functions 	Implementation		**
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
** This GAMMA module provides functions for simulation of field		**
** gradients in NMR spectroscopy.					**
**                                                                      **
*************************************************************************/

#ifndef   GrdDep_cc_			// Is file already included?
#  define GrdDep_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Gradients/GrdDeprec.h>	// Include the interface
#include <Gradients/sys_gradz.h>	// Include z-gradient spin systems
#include <HSLib/HSLibIF.h>		// Include Hilbert space stuff
#include <vector>			// Include libstdc++ STL vectors
#include <Level2/RelaxBas.h>		// Include phenomenological relaxation

// ____________________________________________________________________________
//                         Time Evolution Functions
// ____________________________________________________________________________

void evolve(int N, gen_op& s0, gen_op* Hs, double t, gen_op* sigs) // DEPRECATED!
  { for(int i=0; i<N; i++) sigs[i] = evolve(s0, Hs[i], t); }
void evolve(int N, gen_op& sigma0, gen_op* Us, gen_op* sigmas)	   // DEPRECATED!
  { for(int i=0; i<N; i++) sigmas[i] = evolve(sigma0, Us[i]); }

void evolve(int N, gen_op* s0, gen_op& H, double t, gen_op *sigs)  // DEPRECATED!
  { gen_op U = prop(H,t); evolve(N, s0, U, sigs); }

void evolve(int N, gen_op* s0, gen_op& U, gen_op* sigs)            // DEPRECATED!
  { for(int i=0; i<N; i++) sigs[i] = evolve(s0[i], U); }

void evolve(int N, gen_op* s0, gen_op* Hs, double t, gen_op* sigs) // DEPRECATED!
  { for(int i=0; i<N; i++) sigs[i] = evolve(s0[i], Hs[i], t); }
void evolve(int N, gen_op* s0, gen_op* U, gen_op* sigs)		   // DEPRECATED!
  { for(int i=0; i<N; i++) sigs[i] = evolve(s0[i], U[i]); }

#endif 								// GradDeprec.cc

