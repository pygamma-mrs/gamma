/* HDipolar.cc *************************************************-*-c++-*-*
**									**
**	                        G A M M A	 			**
**								 	**
**	Dipolar Hamiltonians                        Implementation 	**
**								 	**
**      Scott Smith                                                     **
**      Copyright (c) 2001                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**								 	**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
**  Description							 	**
**								 	**
** This module of the GAMMA MR platform provides Hamiltonians for 	**
** dipolar interations.  The functions are typically designed to take   **
** spin system containing dipolar interactions as an input argument.    **
** Differing functions allow users to obtain the Hamiltonian under 	**
** different levels of approximation, at any orientation, an for either **
** all system dipolar interations or just that for a spin pair. Note	**
** that the returned Hamiltonians will reside in the Hilbert space of	**
** the spin system (i.e. the composite Hilbert space.)			**
**						 			**
*************************************************************************/

#ifndef _HDipolar_cc_			// Is file already included?
#define _HDipolar_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// this is the implementation
#endif

#if defined(_MSC_VER)			// If using MSVC++ then we
 #pragma warning (disable : 4786)	// Kill STL namelength warnings too
#endif

#include <IntRank2/HDipolar.h>		// Include the file headers
#include <IntRank2/SolidSys.h>		// Knowledge of solid state systems
#include <string>			// Knowledge of libstdc++ strings 
#include <HSLib/SpinSystem.h>		// Knowldege of spin systems
#include <HSLib/GenOp.h>		// Knowldege of general operators
#include <HSLib/SpinOp.h>		// Knowledge of spin operators
#include <HSLib/SpinOpCmp.h>		// Knowledge of composite spin ops.
#include <HSLib/SpinOpRot.h>		// Knowledge of rotated spin ops.

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
 
gen_op HD0(const solid_sys& sys)
  {
  int i,j,k;					// Spin and dipole indices
  bool het;					// Heteronuclear flag
  int ns = sys.spins();				// Number of spins
  std::vector<int> HSs = sys.HSvect();		// Array of spin Hilbert spaces
  matrix H;					// Hamiltonian matrix
  IntDip D;					// Dipolar interaction
  for(i=0,k=0; i<ns-1; i++)			// Loop over all spin pairs
    for(j=i+1; j<ns; j++,k++)
      {
      D = sys.getDipInt(i, j);			//    Get dipolar interaction
      het = (sys.isotope(i) != sys.isotope(j));	//    Check if heteronuclear
      H += D.H0(HSs, i, j, het);		//    Add spin pair H component
      }
  return gen_op(H);
  }

gen_op HD0(const solid_sys& sys, int i, int j)
  {
  std::vector<int> HSs = sys.HSvect();		// Array of spin Hilbert spaces
  bool het = (sys.isotope(i) != sys.isotope(j));// Flag if heteronuclear
  IntDip D = sys.getDipInt(i, j);		// Dipolar interaction
  return gen_op(D.H0(HSs, i, j, het));		// Spin pair dip. Hamiltonian
  }

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
 
/*
gen_op HD1(const solid_sys& sys)
  {
  int i,j,k;					// Spin and dipole indices
  bool het;					// Heteronuclear flag
  int ns = sys.spins();				// Number of spins
  vector<int> HSs = sys.HSvect();		// Array of spin Hilbert spaces
  matrix H;					// Hamiltonian matrix
  IntDip D;					// Dipolar interaction
  for(i=0,k=0; i<ns-1; i++)			// Loop over all spin pairs
    for(j=i+1; j<ns; j++,k++)
      {
      D = sys.getDipInt(i, j);			//    Get dipolar interaction
      het = (sys.isotope(i) != sys.isotope(j));	//    Check if heteronuclear
      H += D.H0(HSs, i, j, het);		//    Add spin pair H component
      }
  return gen_op(H);
  }


gen_op HD1(const solid_sys& sys) const;
gen_op HD1(const solid_sys& sys, double theta, double phi=0) const;
*/
 
// ____________________________________________________________________________
// C       Summed First & Second Order Dipolar Interaction Hamiltonians
// ____________________________________________________________________________

/*       These are SECULAR (Rotationally Invariant About Bo Field Axis)
    These Are For When The Dipolar Interaction Is A Perturbation To Zeeman
*/ 

//matrix IntDip::Hw(double Om) const

        // Input                sys  : Spin system
        //                      wD   : Dipolar frequency (Hz)
	//			i    : Spin index
        // Output               HDw  : The secular part of the dipolar
        //                             Hamiltonian (default basis, Hz)
	//			       for the spin i
	// Note			     : No asymmetry is considered here
	// Note			     : This is the sum of the 1st & 2nd order
	//			       dipolar interactions (pert. theory)
	// Note			     : This is rotationally invariant about z

//  This function returns only the secular part of the second order dipolar
//  Hamiltonian (from perturbation theory).  Note that this still assumes that
//  the dipolar interaction is a perturbation to to the Zeeman interaction.

//                                (0)    (1)      [1]    [2]
//                       H    =  H    + H     =  H    + H
//                        D       D      D        D      D

//  { return H0() + H1(Om); }

 
//matrix IntDip::Hw(double Om, double theta, double phi) const
 
        // Input                Om	: Field Strength (Larmor in Hz)
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HD	: The 2nd order secular part of the
        //                                dipolar Hamiltonian (default
        //                                basis, Hz) for the interaction
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note				: Also called the 1st order dipolar
        //                                interaction (perturbation theory)
        // Note                      	: Rotationally invariant about z
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1

//  { return H0(theta, phi) + H1(Om, theta, phi); }


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

	   Input		D	: Dipolar interaction
	   Output			: Full system dipolar Hamiltonian
	   Note				: NOT rotationally invariant about z
	   Note				: This will return in the composite
					  spin Hilbert space of the system   */
 
gen_op HD(const solid_sys& sys)
  {
  int i,j;					// Spin and dipole indices
  int ns = sys.spins();				// Number of spins
  std::vector<int> HSs = sys.HSvect();		// Array of spin Hilbert spaces
  matrix H;					// Hamiltonian matrix
  IntDip D;					// Dipolar interaction
  for(i=0; i<ns-1; i++)				// Loop over all spin pairs
    for(j=i+1; j<ns; j++)
      {
      D = sys.getDipInt(i, j);			//    Get dipolar interaction
      H += D.H(HSs, i, j);			//    Add spin pair H component
      }
  return gen_op(H);
  }

gen_op HD(const solid_sys& sys, int i, int j)
  {
  std::vector<int> HSs = sys.HSvect();		// Array of spin Hilbert spaces
  IntDip D = sys.getDipInt(i, j);		// Dipolar interaction
  return gen_op(D.H(HSs, i, j));		// Spin pair dip. Hamiltonian
  }

// ____________________________________________________________________________
// E             DIPOLE DIPOLE HAMILTONIAN GENERIC FUNCTIONS
// ____________________________________________________________________________

/* These Functions Don't Rely On GAMMA's Rank 2 Interaction Classes At All.
   The Functions Just Mimic The Other Hamiltonians Generateed Herein For A
   Single Spin Pair (As May Be Done Also In Class IntRank 2). The Difference 
   Is That The User May Obtain Them Quickly At Specific Orientations Without 
   Having To Create An Interaction Or Spin System.                           */

/*
matrix HD0(double qn, double wDo, double eta, double theta, double phi)
 
	// Input		qn	: Duantum number (1, 1.5, 2.5,...)
	//			wDo     : PAS Dipolar frequency
	//			eta     : Dipolar asymmetry [0,1]
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               H0	: The secular part of the dipolar
        //                                Hamiltonian (default basis, Hz)
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note				: Also called the 1st order dipolar
        //                                interaction (perturbation theory)
        // Note                      	: Rotationally invariant about z
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1
 
//  The secular part of the dipolar Hamiltonian is that which returns
//  only those components which commute with z axis rotations.  Here we have
 
//            [1]             1                   2             (0)
//           H  (theta,phi) = - w (the,phi) * [3*I - I(I+1)] = H
//            D               6  D                z             D
 
// where 

//                      [ 1      2               1        2                   ]
// w (theta,phi) = W    | - [3cos (theta) - 1] + - eta sin (theta)*cos(2*phi) |
//  D               D,o [ 2                      2                            ]

// and                                                                          
//                                    3*DCC
//                            w    = --------
//                             D,o   2I(2I-1)
 
  {
  int Ival = int(2.*qn + 1);			// For 1 spin SOp functions
  matrix IE = Ie(Ival);				// The identity operator
  matrix IZ = Iz(Ival);				// The Fz operator
  double Ctheta = cos(theta*DEG2RAD);		// Need cos(theta)
  double Stheta = sin(theta*DEG2RAD);		// Need sin(theta)
  double C2phi = cos(2.0*phi*DEG2RAD);		// Need cos(2phi)
  double wD = wDo * ((3.0*Ctheta*Ctheta-1.0)	// "Oriented" wD
            + eta*Stheta*Stheta*C2phi);
  return (wD/12.0)*(3*(IZ*IZ) - (qn*(qn+1))*IE);// 1st Order Hamiltonian
  }  


matrix HD1(double Om, double qn, double wDo, double eta,
                                                      double theta, double phi)
 
        // Input                Om	: Field Strength (Larmor in Hz)
	// 			qn	: Duantum number (1, 1.5, 2.5,...)
	//			wDo     : PAS Dipolar frequency
	//			eta     : Dipolar asymmetry [0,1]
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HD1	: The 2nd order secular part of the
        //                                dipolar Hamiltonian (default
        //                                basis, Hz) for the interaction
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note				: Also called the 1st order dipolar
        //                                interaction (perturbation theory)
        // Note                      	: Rotationally invariant about z
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1
	// Note				: This will be zero in PAS if no eta
	//				  and small if Om >> DCC

//  [2]              -1   [ 1    ]2 [                          2     2   
// H   (theta,phi) = -- * | - w  |  | 2*V V  (theta,phi)*I *(4*I - 8I  - 1)
//  D                Om   [ 6  D ]  [    1 -1             z          z
//
//                                                              2    2      ]
//                                  + 2*V V  (theta,phi)*I *(2*I - 2I  - 1) |
//                                       2 -2             z          z      ]
// where
//
//                   1 [ [    2    2                               ]    4
// 2V V  (the,phi) = - | | eta *cos (2*phi) - 6*eta*cos(2*phi) + 9 | cos (the)
//   1 -1            2 [ [                                         ]
//
//                  [      2    2                                2   ]   2
//                - | 2*eta *cos (2*phi) - 6*eta*cos(2*phi) - eta + 9|cos (the)
//                  [                                                ]
//
//                  [    2      2             ] ]
//                + | eta * [cos (2*phi) - 1] | |
//                  [                         ] ]
//
// and
//                   3 [[  1    2    2          1                  3]    4
//  V V  (the,phi) = - || --*eta *cos (2*phi) - -*eta*cos(2*phi) + -| cos (the)
//   2 -2            4 [[ 12                    2                  4]
//
//                   [ -1    2    2          1    2   3 ]    2
//                 + | --*eta *cos (2*phi) + -*eta  - - | cos (theta)
//                   [  6                    3        2 ]
//
//                   [  1    2    2          1                  3 ] ]
//                 + | --*eta *cos (2*phi) + -*eta*cos(2*phi) + - | |
//                   [ 12                    2                  4 ] ]
 
//  in accordance with the article by P.P. Man "Dipolar Interactions" in
//  the Encyclopedia of Magnetic Resonance by Grant and Harris, Vol 6, Ped-Rel,
//  page 3840, Eq. (19), coupled with page 3841, Eq. (32).

  {
  double C2phi = cos(2.0*phi*DEG2RAD);		// cos(2.0*phi)
  double C2phisq = C2phi*C2phi;			// cos(2*phi)*cos(2*phi)
  double Ctheta = cos(theta*DEG2RAD);		// cos(theta)
  double Cthetasq = Ctheta*Ctheta;		// cos(theta)*cos(theta)
  double Ctheta4 = Cthetasq*Cthetasq;		// (cos(theta))^4
  double etasq = eta*eta;			// eta*eta
  double etasqC2phisq = etasq*C2phisq;		// (eta^2)*(cos(2*phi)^2)
  double etaC2phi = eta*C2phi;			// eta*cos(2*phi)

  double V1Vm11 = (etasqC2phisq - 6.0*etaC2phi + 9.0)*Ctheta4;
  double V1Vm12 = (2.0*etasqC2phisq - 6.0*etaC2phi + - etasq + 9.0)*Cthetasq;
  double V1Vm13 = etasq * (C2phisq - 1.0);
  double twoV1Vm1 = 0.5*(V1Vm11 - V1Vm12 + V1Vm13);

  double V2Vm21 = (etasqC2phisq/12.0 - 0.5*etaC2phi + 0.75)*Ctheta4;
  double V2Vm22 = (-etasqC2phisq/6.0 + etasq/3.0 - 1.5)*Cthetasq;
  double V2Vm23 = (etasqC2phisq/12.0 + 0.5*etaC2phi + 0.75);
  double twoV2Vm2 = 1.5*(V2Vm21 + V2Vm22 + V2Vm23);

  int Ival = int(2.*qn + 1);			// For 1 spin SOp functions
  matrix IZ = Iz(Ival);				// The operator Iz
  matrix IE = Ie(Ival);				// The operator E
  matrix Hmx =  twoV1Vm1*IZ*(4.0*qn*(qn+1)*IE - 8.0*IZ*IZ - IE);
  Hmx       += twoV2Vm2*IZ*(2.0*qn*(qn+1)*IE - 2.0*IZ*IZ - IE);
  Hmx       *= (-wDo*wDo/(36.0*Om));
  return Hmx;
  }
*/ 

#endif 								// HDipolar.cc
