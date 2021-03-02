/* HElectronG.h *************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Electron G Hamiltonians 	 	    Interface		**
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
** electron g interations.  These functions are typically designed to	**
** take spin system containing electron g interactions as an input	**
** argument. Functions allow users to obtain the Hamiltonian under      **
** different levels of approximation, at any orientation, an for either **
** all system electron g interations or just that for a particular	**
** spin. Note that the returned Hamiltonians will reside in the Hilbert **
** space of the spin system (i.e. the composite spin Hilbert space.)    **
**                                                                      **
*************************************************************************/

#ifndef   GHG_h_				// Is file already included?
#  define GHG_h_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the interface
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <IntRank2/SolidSys.h>			// Knowledge of spin systems
#include <HSLib/GenOp.h>			// Knowledge of operators

// ____________________________________________________________________________
// A                 First Order Electron G Hamiltonians
// ____________________________________________________________________________

/*
         These are SECULAR (Rotationally Invariant About Bo Field Axis)
    These Are For When The G factor Interaction Is A Perturbation To Zeeman

  The secular part of the G Hamiltonian, returned by functions in this section, 
  contain only those Hamiltonian components which explicitly commute with z-axis
  rotations (i.e. time-independent when expressed in a rotating frame). Note 
  that this Hamiltonian is field depenent since it is the electron response to
  the externally applied field.

  In more explicit mathematical terms, the Hamiltonian returned in this
  section is given by

              spins          spins                    spins
         [1]   ---   [1]      ---     G   G    G       ---   (0)
        H   =  \    H   (i) = \    [Xi * A  * T  ]   = \    H   (i)
         G     /     G        /           20   20 i    /     G
               ---            ---                      ---
                i              i                        i

   where
                    [1]       (0)           Q       Q
                   H   (i) = H   (i) = Xi  A   (i) T   (i)
                    Q         Q          i  2,0     2,0

 
   Note that this is the anisotropic G Hamiltonian in a high field limit!
   Also note that this must be ADDED to the isotropic G Hamiltonian. 

           Input                sys     : Spin system
           Output               HG0     : The secular part of the electron G
           Output               H0	: The secular part of the electron G
                                          Hamiltonian (default basis, Hz)
                                          electron G interaction
                                          for orientation {phi, theta}
                                          from the tensor PAS if set
           Note				: Also called the 1st order electron G
                                          interaction (perturbation theory)
           Note                      	: Rotationally invariant about z
	   Note				: Return is in the spin Hilbert
	  				  space of dimension 2I+1            */
 
MSVCDLL gen_op HG0(const solid_sys& sys);
MSVCDLL gen_op HG0(const solid_sys& sys, int i);

// ____________________________________________________________________________
// B                   2nd Order Electron G Hamiltonians
// ____________________________________________________________________________
 
/*       These are SECULAR (Rotationally Invariant About Bo Field Axis)
    Applicable When The Electron G Interaction Is A Perturbation To Zeeman

   As mentioned in section A, the electron G Hamiltonian is not rotationally
   invariant about the z-axis. In that section we simply ignored all terms
   which are time dependent in the rotating frame. In thise section we build
   Hamiltonians that compensate for the neglect of such terms while still 
   allowing one to express the Hamiltonian in a rotating frame.

   In more explicit mathematical terms, the Hamiltonian returned in this
   section is given by
 
    [2]              -1   [ Xi ]2 [                           2    2   
   H   (theta,phi) = -- * | -- |  | 2*A A  (theta,phi)*I *(4*I - 8I  - 1)
    G                Om   [ 2  ]  [    1 -1             z          z
                           
                                                                2    2      ]
                                    + 2*A A  (theta,phi)*I *(2*I - 2I  - 1) |
                                         2 -2             z          z      ]

   Note that this term will need to be added to the 1st order electron G
   Hamiltonian (as well as the Zeeman Hamiltonian) in order to produce something
   that will be useful in affecting a system.

           Input                sys     : Spin system
           Output               H1	: The 2nd order secular part of the
                                          electron G Hamiltonian (default
                                          basis, Hz) for the interaction
           Note				: Also called the 1st order electron G
                                          interaction (perturbation theory)
	   Note 			: This will be zero in PAS if no eta
	  				  and small if Om >> GCC, i.e. when
					  Zeeman >> Electron G!             */

/*
gen_op HG1(const solid_sys& sys)
gen_op HG1(const solid_sys& sys, int i)
*/

 
// ____________________________________________________________________________
// C         Summed First & Second Order Electron G Hamiltonians
// ____________________________________________________________________________
 
/*       These are SECULAR (Rotationally Invariant About Bo Field Axis)
     Applicable When The Electron G Interaction Is A Perturbation To Zeeman
 
   This is just a blend of the two Hamiltonians returned in sections A and B.
   It should be added to the Zeeman interaction Hamiltonian (in the rotating
   or multiply rotating frame).

                                  (0)    (1)      [1]    [2]
                         H    =  H    + H     =  H    + H
                          G       G      G        G      G

        // Input                sys  : Spin system
	//			i    : Spin index
        // Output               HGw  : The secular part of the electron G
        //                             Hamiltonian (default basis, Hz)
	//			       for the spin i
	// Note			     : This is the sum of the 1st & 2nd order
	//			       electron G interactions (pert. theory)
	// Note			     : This is rotationally invariant about z
	// Note			     : This will return in the spin Hilbert
	//			       space of the spin system              */

/*
gen_op HGw(const solid_sys& sys)
gen_op HGw(const solid_sys& sys, int i)
*/
 
// ____________________________________________________________________________
// D                      Full Electron G Hamiltonians
// ____________________________________________________________________________
 
/*                 These Hamiltonians Assume No Approximations 
       They Are Independent Of The Strength Of Any Externally Applied Field
         They Are Correct In The Laboratory Frame, NOT In A Rotating Frame

                                 spins
                                  --- 
                            H   = \   H (theta , phi )
                             G    /    G      i     i
                                  --- 
                                   i

                              G       m    G                   G
          H (theta ,phi ) = Xi  * (-1)  * A   (theta ,phi ) * T    (i)
           G      i    i      i            2,m      i    i     2,-m   
 
           Input                sys  : Spin system
	   Note				: This will return in the spin
	  				  Hilbert space of spin system       */

MSVCDLL gen_op HG(const solid_sys& sys);
MSVCDLL gen_op HG(const solid_sys& sys, int i);
 
#endif 							//  HElectronG.h

