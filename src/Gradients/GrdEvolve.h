/* GrdEvolve.h *************************************************-*-c++-***
**									**
**	                         G A M M A 				**
**						 			**
**	Field Gradients Time Evolution 	         	Interface	**
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

#ifndef   GrdEvolve_h_			// Is file already included?
#  define GrdEvolve_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Gradients/sys_gradz.h>	// Know z-gradient systems
#include <HSLib/GenOp.h>		// Know general operators
#include <Level2/RelaxBas.h>		// Know phenomenological relaxation
#include <vector>			// Know libstdc++ STL vectors

// ____________________________________________________________________________
// A            Evolve A Single State Into An Array of States
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//                     Hilbert Space Evolution Functions
//-----------------------------------------------------------------------------

/* These functions assume that we have an single density operator that is
   appropriate for all sub-systems in the gradient. However, the each will
   evolve in time under a different Hamiltonian (or equivalently are affected
   by a different propagator.) So the functions return an array of evolved
   density operators, one for each of the subsystems.  The arguments demand
   an initial density operator along with as set of evolution Hamiltonians
   with an evolution time or (equivalently) a set of evolution propagators.
   They then produce a set of evolved density operators.

                                  t
      sigma (t) = U (t)*sigma(0)*U (t) = exp(-i*H *t) * sigma * exp(i*H *t)
           i       i              i              i                     i     */

MSVCDLL std::vector<gen_op> evolve(gen_op& sigma0, std::vector<gen_op>& Hs, double t);
MSVCDLL std::vector<gen_op> evolve(gen_op& sigma0, std::vector<gen_op>& Us);

        // Input                sigma0  : Intitial density operator
        //                      Hs      : Array of Hamiltonians
        //                      t       : Evolution time
        //                      Us      : Array of propagators
        // Output               sigmas  : Array of density operators
        //                                evolved from sigma0 for time t
	//				  under the array of Hamiltonians Hs
	//				  or evolved under the propagators Us

        // Input    DEPRECATED  NSS     : Number of sub-systems


//-----------------------------------------------------------------------------
//                Evolution Under Phenomenological Relaxation
//-----------------------------------------------------------------------------

MSVCDLL std::vector<gen_op> evolve(gen_op& sigma0, std::vector<gen_op>& Hs, RBasic& R, double t);

// -----------------------------------------------------------------------
//  Evolve Array of States Into Array of States Under A Single Propagator
// -----------------------------------------------------------------------

	// Input		NSS	: Number of sub-systems
	//			sigmas0	: Intitial density operators
	//			H	: A Hamiltonian
	//			t	: Evolution time
        //                      sigmas  : Array of density operators
	// Output		void	: The density operators sigmas0 are
	//				  evolved under the same Hamiltonian
	//				  H for the time t into a new array
	//				  of density operators

MSVCDLL std::vector<gen_op> evolve(std::vector<gen_op>& sigmas0, gen_op& H, double t);
MSVCDLL std::vector<gen_op> evolve(std::vector<gen_op>& sigmas0, gen_op& U);


// -----------------------------------------------------------------------
//            Evolve States Array Under Array Of Propagators
// -----------------------------------------------------------------------


	// Input		NSS	: Number of sub-systems
	//			sigmas0	: Intitial density operators
	//			Hs	: An array of  Hamiltonians
	//			t	: Evolution time
        //                      sigmas  : Array of density operators
	// Output		void	: The density operators sigmas0 are
	//				  evolved under their associated 
	//				  Hamiltonians Hs for the time t


	// Input		NSS	: Number of sub-systems
	//			sigmas0	: Intitial density operators
	//			U	: A propagators
        //                      sigmas  : Array of density operators
	// Output		void	: The density operators sigmas0 are
	//				  evolved under the same propagator
	//				  U into a new array of density
	//				  operators

MSVCDLL std::vector<gen_op> evolve(std::vector<gen_op>& sigma0, std::vector<gen_op>& H, double t);
MSVCDLL std::vector<gen_op> evolve(std::vector<gen_op>& sigma0, std::vector<gen_op>& U);

MSVCDLL std::vector<gen_op> evolve(std::vector<gen_op>& sigma0, std::vector<gen_op>& Hs, RBasic& R, double t);

#endif 							// GrdEvolve.h

