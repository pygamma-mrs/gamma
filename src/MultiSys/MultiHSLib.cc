/* MultiHSLib.cc ************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Multiple System Hilbert Space Routines 	    Implementation	**
**                                                                      **
**      Copyright (c) 1995                                              **
**      Nikolai Skrynnikov                                              **
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header:
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** This file of functions support multi_sys, the GAMMA class for        **
** handling mulitple spin systems.  The routines herein generally       **
** involve such a spin system and build up common Hilbert Space library	**
** funcitons. These functions mirror those defined in the HS library 	**
** (see module HSLib) but will return operators that exist in the       **
** composite Hilbert spin space spanning the components in multi_sys.   **
** The Hilbert space is a direct product space of the systems involved. **
**                                                                      **
*************************************************************************/

#ifndef   MultiHSLib_cc			// Is the file already included?
#  define MultiHSLib_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <MultiSys/MultiHSLib.h>	// Include the header file
#include <MultiSys/MultiSys.h>		// Include multi_sys spin systems
#include <HSLib/HSauxil.h>		// Include H.S. Lib sigma_eq
#include <HSLib/GenOp.h>		// Include general operators
#include <Matrix/matrix.h>		// Include GAMMA matrices
#include <LSLib/sys_dynamic.h>		// Include dynamic systems

using std::vector;			// Using libstdc++ STL vectors

// ____________________________________________________________________________
// A                    EQUILIBRIUM DENSITY OPERATORS
// ____________________________________________________________________________
 
/* The density operator associated with a multipe spin system has a slightly
   different normalization that that of a single spin system. The single system
   density operator is approximated by (high-field limit)

                            ---  gamma
                         ~  \         i     homonuclear
                 sigmaeq =  /    ------ Iz  -----------> Fz
                            ---  gamma    i
                             i        0

    Scaling by the relative gyromagnetic ratios produces the proper relative
    polarizations of the spins involved. When multiple systems are blended 
    together we need to insure that three additional aspects are satisfied.

    1.) System density ops must all be relative to the same type (same gamma )
    2.) System density ops must be normalized to their own spin Hilbert spaces.  
    3.) System density ops must be scaled by their respective populations.

    Thus, for a multiple spin system containing components indexed with c,
    the densit operator for each component is given by

                        pop   ---  gamma                       pop
                      ~    c  \         i,c       homonuclear     c 
             sigmaeq  = ----  /    -------- Iz    -----------> ---- Fz
                    c   HS    ---  gamma      i,c              HS
                          c    i        0                        c

    And the system wide density operator places these about the diagonal of the
    operator (in the full composite space) in the order in which the components
    are specified in the system.

	   Input	msys	: A multi_sys spin system
	   Output	MOp     : A density operator for the system at
	  			  equilibrium                                */

gen_op sigma_eq(const multi_sys& msys)
  {
  int nc = msys.NComps();			// Number of components
  vector<matrix> mxc;				// Array of matrices
  vector<matrix> bsc;				// Array of bases
  gen_op Op;					// Scratch operator
  Isotope I   = (msys.Comp(0)).isotope(0);	// 1st comp, 1st spin type
  double  HS0 = (msys.Comp(0)).HS();		// 1st comp. Hilbert space
  double  sf;					// Scaling factor
  for(int i=0; i<nc; i++)			// Loop over spin components
    {
    sf = msys.pop(i)*HS0/(msys.Comp(i)).HS();	//   Scaling component i
    Op = sigma_eq(msys.Comp(i), I);		//   Requested Op, comp i
    mxc.push_back(sf*Op.get_mx());		//   Store its matrix 
    bsc.push_back((Op.get_basis()).U());	//   Store its basis
    }
  return gen_op(mxc, bsc);
  }

#endif							// MultiHSLib.cc
