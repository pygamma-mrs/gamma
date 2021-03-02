/* HDipolar.h ***************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Oriented Dipolar Hamiltonians 		     Implementation	**
**                                                                      **
**                                                                      **
**      Scott Smith                                                     **
**      Copyright (c) 2001                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
** This module of the GAMMA MR platform provides Hamiltonians for       **
** dipolar interations.  The functions are typically designed to take   **
** spin system containing dipolar interactions as an input argument.    **
** Differing functions allow users to obtain the Hamiltonian under      **
** different levels of approximation, at any orientation, an for either **
** all system dipolar interations or just that for a spin pair. Note    **
** that the returned Hamiltonians will reside in the Hilbert space of   **
** the spin system (i.e. the composite Hilbert space.)                  **
**                                                                      **
*************************************************************************/

#ifndef   HDipolar_h_			// Is file already included?
#  define HDipolar_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <IntRank2/SolidSys.h>		// Knowledge of spin systems
#include <HSLib/GenOp.h>		// Knowledge of operators

// ____________________________________________________________________________
// A                1st Order Dipolar Interaction Hamiltonians
// ____________________________________________________________________________

/*       These are SECULAR (Rotationally Invariant About Bo Field Axis)
   Applicable When The Dipolar Interaction Is A Small Perturbation To Zeeman

   The returned Hamiltonian exists in a (possibly multiple) rotating frame, as
   set up by the input spin system.  However, dipolar interactions are not
   rotationally invariant about the externally applied field axis (+z in GAMMA)
   and are actually time dependent when represented in the rotating frame(s).  
   Here, we assume that all dipolar interactions are weak relative to Zeeman
   interactions - i.e. direct spin interactions with the externally applied
   magnetic field. Thus, we can ignore any components of the dipolar 
   Hamiltonian that become time dependent in the rotating frame(s). Exactly
   what terms these are depend upon whether the spin pairs involved are homo-
   or hetero-nuclear.

   For example, consider a heteronuclear spin pair in the case that the 
   returned dipolar Hamiltonian is to be in the rotating frame of both I and
   S (or just that of I or S).  The dipolar Hamiltonian is time dependent 
   due to the I+S- and the I-S+ terms of T20 as well as from any terms in
   T2{1,-1,2,-2}.  So, we simply assume there are insignificant relative to
   the Zeeman Hamiltonian and ignore them, leaving only the IzSz term of T20.
   A homonuclear spin pair in the same rotating frame(s) has time dependent
   dipolar Hamiltonian terms only in the T2{1,-1,2,-1} terms so in that case
   we keep everything from T20 but nothing else.  In short, we only keep
   terms that commute with the z-axis rotation operator.

   In more explicit mathematical terms, the Hamiltonian returned in this
   section is given by

       spins spins          spins spins                   spins spins
  [1]   ---   ---  [1]       ---   ---    D   D   'D       ---   ---   (0)
 H   =  \     \   H  (i,j) = \     \   [Xi * A  * T  ]   = \     \    H   (i,j)
  D     /     /    D         /     /          20   20 ij   /     /     D       
        ---   ---            ---   ---                     ---   --- 
         i    j>i             i    j>i                      i    j>i


   where
                  [1]         (0)              D       'D
                 H   (i,j) = H   (i,j) = Xi   A  (i,j) T   (i,j)
                  D           D            ij  2,0      2,0

   and                    1/2                                   1/2
           ' D         [1]                          hetero   [4]
           T   (i,j) = |-| * [3* I  * I   - I . I ] -------> |-| * I  * I
            2,0        [6]        zi   zj    i   j  nuclear  [6]    zi   zj


           Input                sys	: Spin system
           Output               HD0	: The secular part of the dipolar
                                          Hamiltonian (default basis, Hz)
                                          dipolar interaction
           Note				: Also called the 1st order dipolar
                                          interaction (perturbation theory)
           Note                      	: Rotationally invariant about z
	   Note				: This will return in the composite
					  spin Hilbert space of the system   */
 
MSVCDLL gen_op HD0(const solid_sys& sys);
MSVCDLL gen_op HD0(const solid_sys& sys, int i, int j);

// ____________________________________________________________________________
// B                2nd Order Dipolar Interaction Hamiltonians
// ____________________________________________________________________________
 
/*       These are SECULAR (Rotationally Invariant About Bo Field Axis)
      Applicable When The Dipolar Interaction Is A Perturbation To Zeeman
            This Term Is Dependent Upon The External Field Strength

  As is dicussed in section A, the dipolar Hamiltonian is time-dependent in the
  rotating (or multiply rotating) frame. For the 1st order Hamiltonian we 
  simply throw out any time dependencies, assuming they will be insignificant
  relative to the Zeeman interaction the spins have with an externally applied
  magnetic field. This function assumes that, although weak, the time-dependent
  terms will make a significant contribution to the Hamiltonian and should not
  be ignored. The additonal term is given by
 
    [2]              -1   [ 1    ]2 [                          2     2      
   H   (theta,phi) = -- * | - w  |  | 2*A A  (theta,phi)*Iz*(4*I - 8Iz - 1)
    D                Om   [ 3  D ]  [    1 -1
  
                                                                2     2     ]
                                    + 2*A A  (theta,phi)*Iz*(2*I - 2Iz - 1) |
                                         2 -2                               ]
   
  Note that this must be added to the 1st order dipolar Hamiltonian to produce
  a proper Hamiltonian for use in computations. Also note that the above is
  field dependent as a result of its competition with the Zeeman interaction.

           Input                Om      : Field Strength (Larmor in Hz)
           Output               HD1     : The 2nd order secular part of the
                                          Dipolar Hamiltonian (default
                                          basis, Hz) for the interaction
           Note                         : Also called the 1st order Dipolar
                                          interaction (perturbation theory)  */


//matrix IntDip::HD1(const solid_sys& sys) const;
//matrix IntDip::HD1(const solid_sys& sys, double theta, double phi=0) const;
 
// ____________________________________________________________________________
// D                  Full Dipolar Interaction Hamiltonians
// ____________________________________________________________________________

/*        These Have No Approximations & Are Independent Of Field Strength
        (They Are Correct in The Laboratory Frame, NOT in a Rotating Frame)

  In this case we simply sum over 5 irreducible rank 2 spatial an spin tensor
  components to produce a dipolar Hamiltonian for a paticular spin pair in the
  system.  The full dipolar Hamiltonian is the sum over all spin pair dipolar
  Hamiltonians.

      spins spins         spins spins -2,2
       ---   ---           ---   ---  ---   D       m    D           D
 H  =  \     \   H (i,j) = \     \    \   Xi  * (-1)  * A   (i,j) * T    (i,j)
  D    /     /    D        /     /    /     ij           2,m         2,-m
       ---   ---           ---   ---  ---
        i    j>i            i    j>i   m

           Input                D       : Dipolar interaction
           Output                       : Full system dipolar Hamiltonian
           Note                         : NOT rotationally invariant about z
           Note                         : This will return in the composite
                                          spin Hilbert space of the system   */

MSVCDLL gen_op HD(const solid_sys& sys);
MSVCDLL gen_op HD(const solid_sys& sys, int i, int j);


#endif 							//  HDipolar.h

