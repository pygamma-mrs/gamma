/* MutExch.h  ****************************************************
**								**
**                           G A M M A 				**
**                                				**
**      Mutual Exchange                         Interface	**
**								**
**      Copyright (c) 2001 					**
**      Scott A. Smith						**
**      Nikolai Skrynnikov					**
**                                				**
**      National High Magnetic Field Laboratory			**
**      1800 E. Paul Dirac Drive				**
**      Tallahassee, Florida, 32310				**
**                                				**
**      $Header: $
**                                				**
*****************************************************************/

/*****************************************************************
**								**
** Description							**
**								**
** The following functions herein provide easy access to mutual **
** exchange matrices in spin Hilbert space. These functions are **
** based on defined mutual exchange processes (class ExchProcM)	**
** and any GAMMA spin system.					**
**								**
*****************************************************************/


#ifndef   MutExch_h_			// Is this file already included?
#  define MutExch_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/matrix.h>		// Know GAMMA matrices
#include <vector>			// Know libstdc++ STL vectors
#include <LSLib/SuperOp.h>		// Know superoperators
#include <LSLib/sys_dynamic.h>		// Know anisotropic systems
#include <Level1/ExProcessM.h>		// Know mutual exchange processes

//____________________________________________________________________________
// A                   Exchange Permutation Matrices
//____________________________________________________________________________

/* This function sets up a "pre-conditioned" exchange matrix based on an input
   spin system (to define a spin Hilbert space) and a mutual exchange process
   (to indicate which spins are exchanging, and how). There is no checking to
   insure the spins indicated in the exchange are valid with respect to the
   spin system, nor that the spins in exchange are even of the same isotope
   type. Essentially, the returned array, K, contains the value of 1/2 at each
   <a|K|b> and <b|K|a> where the states a & b are in exchange.

           Input        sys     : A spin system
                        XP      : A mutual exchange process
           Output       Kmx     : Mutual exchange matrix:w
           Note                 : There is not checking done
                                  of spin indices or even
                                  spin conservation.                        */

MSVCDLL matrix Kex(const spin_sys& sys, const ExchProcM& XP);


//____________________________________________________________________________
// B                         Exchange Superoperators
//____________________________________________________________________________

/* These functions return exchange superoperators for a specified spin system
   spin system (to define a spin Liouville space) & mutual exchange process
   or array of such exchange processes (to indicate the spin in exchange &
   how they are exchanging). There is no checking to insure any of the spins
   indicated in the exchange(s) are valid with respect to the spin system,
   nor that the spins in exchange are even of the same isotope type.

   The routines will generate the exchange superoperator for a particular
   exchange process according to

                     ^                              ^   ^
                     ^            ^   ^   ^   ^     ^   ^
                     X   = K  * [ K x K - E x E ] = K - E
                      ex

    and the total exchange superoperator will sum over all single
    exchange superoperators.

           Input        sys     : A spin system
                        XPs     : An array of mutual exchange process
                        Bs      : Hilbert space basis
           Output       Kmx     : Mutual exchange superoperator             */

MSVCDLL super_op Kex(const spin_sys& sys,
                          const std::vector<ExchProcM>& XPs, const basis& Bs);

MSVCDLL super_op Kex(const sys_dynamic& sys, const basis& Bs);

#endif						// MutExch.cc
