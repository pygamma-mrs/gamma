/* MutExch.cc ****************************************************
**                                                              **
**                           G A M M A                          **
**                                                              **
**      Mutual Exchange                         Interface       **
**                                                              **
**      Copyright (c) 2001                                      **
**      Scott A. Smith                                          **
**      Nikolai Skrynnikov                                      **
**                                                              **
**      National High Magnetic Field Laboratory                 **
**      1800 E. Paul Dirac Drive                                **
**      Tallahassee, Florida, 32310                             **
**                                                              **
**      $Header: $
**                                                              **
*****************************************************************/

/*****************************************************************
**                                                              **
** Description                                                  **
**                                                              **
** The following functions herein provide easy access to mutual **
** exchange matrices in spin Hilbert space. These functions are **
** based on defined mutual exchange processes (class ExchProcM) **
** and any GAMMA spin system.                                   **
**                                                              **
*****************************************************************/

#ifndef   MutExch_cc_			// This file already included?
#  define MutExch_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Level2/MutExch.h>		// Include the header
#include <Level1/ExProcessM.h>		// Include mutual exchange
#include <HSLib/SpinSys.h>		// Include base spin systems
#include <LSLib/SuperOp.h>		// Include superoperators
#include <LSLib/sys_dynamic.h>		// Include anisotropic systems

using std::vector;			// Using libstdc++ STL vectors

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

	   Input	sys	: A spin system
	  		XP	: A mutual exchange process
	   Output	Kmx	: Mutual exchange matrix:w
	   Note			: There is not checking done
	  			  of spin indices or even
	  			  spin conservation.                        */

matrix Kex(const spin_sys& sys, const ExchProcM& XP)
  {
  int    ns = sys.spins();			// # of system spins
  matrix P(ns,ns, i_matrix_type);  		// Permutation array
  int    i=0, j=0;				// Spin indices
  for(int l=0; l<XP.NSpins()-1; l++)		// Loop involved spin pairs
    {
    i = XP.Comp(l);				//   Get 1st spin in pair
    j = XP.Comp(l+1);				//   Get 2nd spin in pair
    P.put(complex1,j,i);			//   Set <j|P|i> = 1
    P.put(complex0,i,i);			//   Set <i|P|i> = 0
    }
  i = j;					// Close cycle by setting
  j = XP.Comp(0);				// last 1st, 1st last
  P.put(complex1,j,i);				// Set <j|P|i> = 1
  P.put(complex0,i,i);				// Set <i|P|i> = 0

/* At this point, P is a "permutation" array for all spins exchanging in this
   particular (kth) mutual exchange process, and K is the rate at which they
   are exchanging (1/sec). For A 2-spin exchange, P is a true permutation
   array. For a multi-spin exchange process, P may have multiple elements on
   a single row and a single column, each row and column being normalized.
   This implies that any one spin is actually going to multiple sites.      */

  matrix Mz  = sys.qStates();			// <bf|Mz|spin> array
  matrix PMz = Mz*P;				// Permuted Mz (basis fncts.)
  int hs     = sys.HS();			// Get spin Hilbert space
  matrix K(hs, hs, complex0); 			// Exchange mx (Hilb. space)
  int bf, pbf;					// Basis function indices
  bool match;					// Flag for basis match
  for(bf=0; bf<hs; bf++)			// Loop basis functions
    {						// (these are pbf indices)
    for(pbf=0; pbf<hs; pbf++)			//   Loop basis functions
      {						//   (permuted pbf indices)
      match = true;				//   Assume bf & pbf don't match
      for(i=0; i<ns && match; i++)		// Loop through spins and see if
        if(Mz.get(bf,i) != PMz.get(pbf,i))	// they do indeed match.  All spin
	  match=false;				// states must be the same
      if(match)					// For a match, add the appropriate
        { 					// terms to exchange mx (Hilbert)
        K.put_h(K.get(pbf,bf)+0.5,pbf,bf);	// This should keep K as Hermitian 
        if(pbf == bf)
          K.put_h(K.get(pbf,bf)+0.5,pbf,bf);
        }
      }						// Next permuted basis fnct.
    }						// Try next basis function
  return K;					// Return exchange matrix
  }

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

	   Input	sys	: A spin system
	  		XPs	: An array of mutual exchange process
			Bs	: Hilbert space basis
	   Output	Kmx	: Mutual exchange superoperator             */


super_op Kex(const spin_sys& sys,const vector<ExchProcM>& XPs,const basis& Bs)
  {
/*                                                           ^
                                         ^                   ^
                First construct operator E and superoperator E              */

  int hs = sys.HS();				// Get Hilbert space
  gen_op E(matrix(hs,hs,i_matrix_type));	// Get E (identity) operator
  E.Op_base(Bs);				// Set E basis that of Bs
  super_op UtransE = U_transform(E);		// Get ExE superoperator

/*                                                     ^
                                                       ^
                Now Loop Over Exchange Processes & Sum X
                                                        ex                  */

  gen_op   KOp;					// Temporary operator
  super_op Xp, X;				// Working superoperators
  double K;					// Exchange rate (1/sec)
  for(unsigned k=0; k<XPs.size(); k++)		// Loop over exchange rates
    {
    K = XPs[k].Kex();				//   Exchange rate process k
    KOp = gen_op(Kex(sys, XPs[k]));		//   Pre-X for process k
    KOp.Op_base(Bs);				//   Put Pre-X in basis of Op
    Xp = -K*(U_transform(KOp) - UtransE); 	//   Set kth exchange superop
    X += Xp;					//   Add to full exchange LOp
    }
  return X;					// Return exchange superop
  }

super_op Kex(const sys_dynamic& sys, const basis& Bs)
  { return Kex(sys, sys.MExProcs(), Bs); }

#endif						// MutExch.cc
