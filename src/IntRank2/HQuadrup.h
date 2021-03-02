/* HQuadrup.h ***************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Quadrupolar Hamiltonians 	 	    Implementation	**
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
** quadrupolar interations.  These functions are typically designed to  **
** take spin system containing quadrupolar interactions as an input     **
** argument. Functions allow users to obtain the Hamiltonian under      **
** different levels of approximation, at any orientation, an for either **
** all system quadrupolar interations or just that for a particular     **
** spin. Note that the returned Hamiltonians will reside in the Hilbert **
** space of the spin system (i.e. the composite spin Hilbert space.)    **
**                                                                      **
*************************************************************************/

#ifndef   HQuad_h_				// Is file already included?
#  define HQuad_h_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the interface
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <IntRank2/SolidSys.h>			// Knowledge of spin systems
#include <HSLib/GenOp.h>			// Knowledge of operators

// ____________________________________________________________________________
// A                   1st Order Quadrupolar Hamiltonians
// ____________________________________________________________________________

/*        These are SECULAR (Rotationally Invariant About Bo Field Axis)
    Applicable When Quadrupolar Interaction Is A Small Perturbation To Zeeman

   The returned Hamiltonian exists in a single ) rotating frame, as set up by
   the input spin system.  However, quadrupolar interactions are not
   rotationally invariant about the externally applied field axis (+z in GAMMA)
   and are actually time dependent when represented in the rotating frame(s).
   Here, we assume that all quadrupolar interactions are weak relative to
   Zeeman interactions - i.e. direct spin interactions with the externally
   applied magnetic field. Thus, we can ignore any components of the
   quadrupolar Hamiltonian that become time dependent in the rotating frame.
   These will be the Hamiltonian components from the T2{1,-1,2,-2} terms,
   so we keep everything from T20 but nothing else.  In short, we only keep
   terms that commute with the z-axis rotation operator.

   In more explicit mathematical terms, the Hamiltonian returned in this
   section is given by

                [1]               Q   Q               Q      (0)
               H  (theta,phi) = Xi * A (theta,phi) * T    = H
                Q                     2,0             2,0    Q

           Input                sys     : Spin system
           Output               HQ0     : The secular part of the quadrupolar
                                          Hamiltonian (default basis, Hz)
                                          quadrupolar interaction
                                          for orientation {phi, theta}
                                          from the tensor PAS if set
           Note                         : Also called the 1st order quadrupolar
                                          interaction (perturbation theory)
           Note                         : Rotationally invariant about z
           Note                         : Return is in the spin Hilbert
                                          space of dimension 2I+1            */


MSVCDLL gen_op HQ0(const solid_sys& sys);
MSVCDLL gen_op HQ0(const solid_sys& sys, int i);

// ____________________________________________________________________________
// B                   2nd Order Quadrupolar Hamiltonians
// ____________________________________________________________________________

/*       These are SECULAR (Rotationally Invariant About Bo Field Axis)
    Applicable When The Quadrupolar Interaction Is A Perturbation To Zeeman

   As mentioned in section A, the quadrupolar Hamiltonian is not rotationally
   invariant about the z-axis. In that section we simply ignored all terms
   which are time dependent in the rotating frame. In thise section we build
   Hamiltonians that compensate for the neglect of such terms while still
   allowing one to express the Hamiltonian in a rotating frame.

   In more explicit mathematical terms, the Hamiltonian returned in this
   section is given by

    [2]              -1   [ Xi ]2 [                           2    2
   H   (theta,phi) = -- * | -- |  | 2*A A  (theta,phi)*I *(4*I - 8I  - 1)
    Q                Om   [ 2  ]  [    1 -1             z          z

                                                                2    2      ]
                                    + 2*A A  (theta,phi)*I *(2*I - 2I  - 1) |
                                         2 -2             z          z      ]


   Note that this term will need to be added to the 1st order quadrupolar
   Hamiltonian (as well as the Zeeman Hamiltonian) in order to produce something
   that will be useful in affecting a system.

           Input                sys  : Spin system
           Output               H1      : The 2nd order secular part of the
                                          quadrupolar Hamiltonian (default
                                          basis, Hz) for the interaction
           Note                         : Also called the 1st order quadrupolar
                                          interaction (perturbation theory)
           Note                         : This will be zero in PAS if no eta
                                          and small if Om >> QCC, i.e. when
                                          Zeeman >> Quadrupolar!             */

MSVCDLL gen_op HQ1(const solid_sys& sys);
MSVCDLL gen_op HQ0(const solid_sys& sys, int i);


// ____________________________________________________________________________
// C         Summed First & Second Order Quadrupolar Hamiltonians
// ____________________________________________________________________________

/*       These are SECULAR (Rotationally Invariant About Bo Field Axis)
     Applicable When The Quadrupolar Interaction Is A Perturbation To Zeeman

   This is just a blend of the two Hamiltonians returned in sections A and B.
   It should be added to the Zeeman interaction Hamiltonian (in the rotating
   or multiply rotating frame).

                                  (0)    (1)      [1]    [2]
                         H    =  H    + H     =  H    + H
                          Q       Q      Q        Q      Q

        // Input                sys  : Spin system
        //                      i    : Spin index
        // Output               HQw  : The secular part of the quadrupolar
        //                             Hamiltonian (default basis, Hz)
        //                             for the spin i
        // Note                      : This is the sum of the 1st & 2nd order
        //                             quadrupolar interactions (pert. theory)
        // Note                      : This is rotationally invariant about z
        // Note                      : This will return in the spin Hilbert
        //                             space of the spin system              */

MSVCDLL gen_op HQw(const solid_sys& sys);
MSVCDLL gen_op HQw(const solid_sys& sys, int i);


// ____________________________________________________________________________
// D                      Full Quadrupolar Hamiltonians
// ____________________________________________________________________________

/*                 These Hamiltonians Assume No Approximations
       They Are Independent Of The Strength Of Any Externally Applied Field
         They Are Correct In The Laboratory Frame, NOT In A Rotating Frame

                                 spins
                                  ---
                            H   = \   H (theta , phi )
                             Q    /    Q      i     i
                                  ---
                                   i

                              Q       m    Q                   Q
          H (theta ,phi ) = Xi  * (-1)  * A   (theta ,phi ) * T    (i)
           Q      i    i      i            2,m      i    i     2,-m

           Input                sys  : Spin system
           Note                         : This will return in the spin
                                          Hilbert space of spin system       */

MSVCDLL gen_op HQ(const solid_sys& sys);
MSVCDLL gen_op HQ(const solid_sys& sys, int i);

#endif							//  HQuadrup.h

