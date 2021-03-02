/* Gradients2.h ************************************************-*-c++-***
**									**
**	                         G A M M A 				**
**						 			**
**	Field Gradients Module 2 	         	Interface	**
**						 			**
**      Copyright (c) 1996                                              **
**      Dr. Scott A. Smith             
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

#ifndef   Gradients2_h_			// Is file already included?
#  define Gradients2_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Gradients/sys_gradz.h>	// Know z-gradient systems
#include <HSLib/GenOp.h>		// Know general operators
#include <vector>			// Know libstdc++ STL vectors

//************************************************************************
//                         Hamiltonian Functions
//************************************************************************

// These routines assume that a gradient is applied along a particular
// axis for a set time.  The spin system both before and after the
// applied gradient can be represented by a single density operator        

MSVCDLL void Hzgrad(const sys_gradz& sys, gen_op& H0, gen_op* H);
MSVCDLL std::vector<gen_op> Hzgrad(const sys_gradz& sys, gen_op& H0);

	// Input		sys	: System in a Bo z-Gradient
	//			Ho	: Active Hamiltonian, no Gradient
	//			H	: Array of Hamiltonians
	// Output		void	: The array of Hamiltonains H
	//				  will be altered according to
	//				   H(i) = H0 + HZ(gradz,i) 


//************************************************************************
//                         Propagator Functions
//************************************************************************

MSVCDLL void Props(int NSS, gen_op* Hs, double t, gen_op* Us);
MSVCDLL std::vector<gen_op> Props(std::vector<gen_op>& Hs, double t);
 
        // Input                NSS     : Number of sub-systems
        //                      Hs      : Array of Hamiltonians
        //                      t       : Evolution time
	//			Us	: Array of propagators
        // Output               void	: The array of propagators that
        //                                correspond to evolution under
        //                                the assoicated Hamiltonian
        //                                for time t 
 
#endif 							// Gradients2.h

