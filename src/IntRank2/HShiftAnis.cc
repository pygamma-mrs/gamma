/* HShiftAnis.cc ***********************************************-*-c++-*-*
**									**
**	                        G A M M A	 			**
**								 	**
**	Shift Anisotropy Hamiltonians 			Implementation 	**
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

#ifndef _HShiftAnis_cc_			// Is file already included?
#define _HShiftAnis_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// this is the implementation
#endif

#include <IntRank2/HShiftAnis.h>	// Include the file headers
#include <string>			// Knowledge of libstdc++ strings 
#include <HSLib/SpinSystem.h>		// Knowldege of spin systems
#include <HSLib/GenOp.h>		// Knowldege of general operators
#include <HSLib/SpinOp.h>		// Knowledge of spin operators
#include <HSLib/SpinOpCmp.h>		// Knowledge of composite spin ops.
#include <HSLib/SpinOpRot.h>		// Knowledge of rotated spin ops.

// ____________________________________________________________________________
// B           1st Order Shift Anisotropy Interaction Hamiltonians
// ____________________________________________________________________________

/*       These are SECULAR (Rotationally Invariant About Bo Field Axis)
    These Are For When The CSA Interaction Is A Perturbation To Zeeman

  The returned Hamiltonian exists in a rotating frame, as set up by the input
  spin system.  However, shift anisotropy interactions are not rotationally
  invariant about the externally applied field axis (+z in GAMMA) and are 
  actually time dependent when represented in the rotating frame(s).
  Here, we assume that all shift anisotropy interactions are weak relative to 
  Zeeman interactions - i.e. direct spin interactions with the externally 
  applied magnetic field. Thus, we can ignore any components of the shift 
  anisotropy Hamiltonian that become time dependent in the rotating frame(s). 
  What remains is called the secular part of the shift anisotropy Hamiltonian.
  This funciton returns the secular part of the Hamiltonian, that which will
  commute with z axis rotations.

  In more explicit mathematical terms, the Hamiltonian returned in this
  section is given by

              spins         spins                     spins 
         [1]   ---   [1]     ---     SA   SA    SA     ---   (0)
        H   =  \    H  (i) = \    [Xi  * A   * T  ]  = \    H   (i)
         SA    /     SA      /            20    20 i   /     SA
               ---           ---                       --- 
                i             i                         i 


	   Input		SA	: CSA interaction
                                theta	: Orientation angle (degrees)
                                phi     : Orientation angle (degrees)
           Output               H0	: Secular part of shift anisotropy
                                          Hamiltonian (default basis, Hz)
           Note				: Also called the 1st order SA
                                          interaction (perturbation theory)
           Note                      	: Rotationally invariant about z
	   Note				: Return is in the spin Hilbert
	  				  space of the full spin system      */
 

gen_op HSA0(const solid_sys& sys)
  {
  int ns = sys.spins();                         // Number of spins
  std::vector<int> HSs = sys.HSvect();               // Array of spin Hilbert spaces
  double Om;					// Larmor frequency
  matrix H;                                     // Hamiltonian matrix
  IntCSA SA;                                    // Shift anisotropy interaction
  for(int i=0; i<ns; i++)                       // Loop over all spins
    {
    SA = sys.getCSAInt(i);			//    Get SA interact.
    Om = sys.Omega(i);				//    Larmor in MHz
    H += SA.H0(HSs, i, Om);			//    Add spin pair H component
    }
  return gen_op(H);
  }

gen_op HSA0(const solid_sys& sys, int i)
  {
  std::vector<int> HSs = sys.HSvect();		// Array of spin Hilbert spaces
  IntCSA SA = sys.getCSAInt(i);			// Shift anisotropy interaction
  double Om = sys.Omega(i);			//    Larmor in MHz
  return gen_op(SA.H0(HSs, i, Om));		// Spin SA Hamiltonian
  }

 
// ----------------------------------------------------------------------------
//              Second Order CSA Interaction Hamiltonians
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
//  These Are For When The CSA Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------

//matrix IntCSA::H1(double Om) const
 
        // Input                Om	: Field Strength (Larmor in Hz)
        // Output               H1	: The 2nd order secular part of the
        //                                SA Hamiltonian (default
        //                                basis, Hz) for the interaction
        // Note				: Also called the 1st order SA
        //                                interaction (perturbation theory)
	// Note				: This will be zero in PAS if no eta
	//				  and small if Om >> CSA

//  [2]              -1   [ Xi ]2 [                           2    2   
// H   (theta,phi) = -- * | -- |  | 2*A A  (theta,phi)*I *(4*I - 8I  - 1)
//  SA               Om   [ 2  ]  [    1 -1             z          z
//                         
//                                                              2    2      ]
//                                  + 2*A A  (theta,phi)*I *(2*I - 2I  - 1) |
//                                       2 -2             z          z      ]

//  {
//  int Ival = int(2.*I + 1);			// For 1 spin SOp functions
//  matrix IZ = Iz(Ival);				// The operator Iz
//  matrix IE = Ie(Ival);				// The operator E
//  complex twoA1Am1 = 2.0*Asph[1]*Asph[2];
//  complex twoA2Am2 = 2.0*Asph[3]*Asph[4];
//  double Xi = xi();
//  matrix Hmx =  twoA1Am1*IZ*(4.0*I*(I+1)*IE - 8.0*IZ*IZ - IE);
//  Hmx += twoA2Am2*IZ*(2.0*I*(I+1)*IE - 2.0*IZ*IZ - IE);
//  Hmx *= (-Xi*Xi/(4.*Om));
//  return Hmx;
//  }

 
//matrix IntCSA::H1(double Om, double theta, double phi) const
 
        // Input                Om	: Field Strength (Larmor in Hz)
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HC	: The 2nd order secular part of the
        //                                shift anisotropy Hamiltonian (default
        //                                basis, Hz) for the interaction
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note				: Also called the 1st order shift anisotropy
        //                                interaction (perturbation theory)
        // Note                      	: Rotationally invariant about z
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1
	// Note				: This will be zero in PAS if no eta
	//				  and small if Om >> QCC

//  [2]              -1   [ Xi ]2 [                           2    2   
// H   (theta,phi) = -- * | -- |  | 2*A A  (theta,phi)*I *(4*I - 8I  - 1)
//  Q                Om   [ 2  ]  [    1 -1             z          z
//                         
//                                                              2    2      ]
//                                  + 2*A A  (theta,phi)*I *(2*I - 2I  - 1) |
//                                       2 -2             z          z      ]

//  {
//  int Ival = int(2.*I + 1);			// For 1 spin SOp functions
//  matrix IZ = Iz(Ival);				// The operator Iz
//  matrix IE = Ie(Ival);				// The operator E
//  complex Asph1 = A1(theta, phi);
//  complex Asph2 = A2(theta, phi);
//  complex twoA1Am1 = -2.0*Asph1*conj(Asph1);
//  complex twoA2Am2 = 2.0*Asph2*conj(Asph2);
//  double Xi = xi();
//  matrix Hmx =  twoA1Am1*IZ*(4.0*I*(I+1)*IE - 8.0*IZ*IZ - IE);
//  Hmx += twoA2Am2*IZ*(2.0*I*(I+1)*IE - 2.0*IZ*IZ - IE);
//  Hmx *= (-Xi*Xi/(4.0*Om));
//  return Hmx;
//  }

 
// ----------------------------------------------------------------------------
//      Summed First & Second Order CSA Interaction Hamiltonians
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
//  These Are For When The CSA Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------
 

//matrix IntCSA::Hw(double Om) const

        // Input                sys  : Spin system
        //                      wQ   : CSA frequency (Hz)
	//			i    : Spin index
        // Output               HCw  : The secular part of the shift anisotropy
        //                             Hamiltonian (default basis, Hz)
	//			       for the spin i
	// Note			     : No asymmetry is considered here
	// Note			     : This is the sum of the 1st & 2nd order
	//			       shift anisotropy interactions (pert. theory)
	// Note			     : This is rotationally invariant about z

//  This function returns only the secular part of the second order shift anisotropy
//  Hamiltonian (from perturbation theory).  Note that this still assumes that
//  the shift anisotropy interaction is a perturbation to to the Zeeman interaction.

//                                (0)    (1)      [1]    [2]
//                       H    =  H    + H     =  H    + H
//                        Q       Q      Q        Q      Q

//  { return H0() + H1(Om); }

 
//matrix IntCSA::Hw(double Om, double theta, double phi) const
 
        // Input                Om	: Field Strength (Larmor in Hz)
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HC	: The 2nd order secular part of the
        //                                shift anisotropy Hamiltonian (default
        //                                basis, Hz) for the interaction
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note				: Also called the 1st order shift anisotropy
        //                                interaction (perturbation theory)
        // Note                      	: Rotationally invariant about z
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1

//  { return H0(theta, phi) + H1(Om, theta, phi); }


// ----------------------------------------------------------------------------
//                      Full CSA Interaction Hamiltonians
//       These Have No Approximations & Are Independent of Field Strength
//     (They Are Correct in The Laboratory Frame, NOT in a Rotating Frame) 
// ----------------------------------------------------------------------------

/* Here we make no approximations. As such, these Hamiltonians are not those
   used in the rotating frame. Here we simply sum over the 5 spatial and spin
   tensor components to produce the Hamiltonian. Recall that the 5 spatial
   tensor components, indexed by the spin angular momentum component
   m = {-2,-1,0,1,2}, map into the Asph array by the following:

                        {  0->0, 1->1, -1->2, 2->2, -2->4 }

                       -2,2
                        ---   SA       m    SA                  SA
     H  (theta, phi) =  \   Xi   * (-1)  * A   (theta, phi) * T
      SA                /                   2,m                2,-m
                        ---
                         m
   Note: There are no m = +/-2 spin components in rank 2 CSA interactions

	   Input		Q	: CSA interaction
	   Note				: Return is in the spin Hilbert
	  				  space of dimension 2I+1            */

/*
matrix IntCSA::H( ) const
  {
  matrix Hmx = Acomp(0)*T0;			// First set the m=0 terms
  if(norm(Acomp(1)))				// Then add in m=+/-1 terms
    Hmx -= (Acomp(1)*T1+Acomp(2)*Tm1);		// if they are present
  return xi()*Hmx;
  }

matrix IntCSA::H(double theta, double phi) const
  {
  matrix Hmx  = A0( theta,phi)*T0;
         Hmx -= Am1(theta,phi)*T1 + A1(theta,phi)*Tm1;
  return xi()*Hmx;
  }

matrix IntCSA::H(const vector<int>& HSs, int i) const
  {
  matrix Hmx  = Acomp(0)*T20(HSs,i);
         Hmx -= Acomp(1)*T21(HSs,i) + Acomp(2)*T2m1(HSs,i);
  return xi()*Hmx;
  }

matrix IntCSA::H(const vector<int>& HSs, int i, double T, double P) const
  {
  matrix Hmx  = A0( T,P)*T20(HSs,i);
         Hmx -= Am1(T,P)*T21(HSs,i) + A1(T,P)*T2m1(HSs,i);
  return xi()*Hmx;
  }
*/

#endif 							// HShiftAnis.cc
