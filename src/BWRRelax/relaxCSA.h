/* relaxCSA.h ***************************************************-*-c++-*-
**									**
**	                           G A M M A				**
**									**
**	CSA Relaxation Functions                       Interface	**
**									**
**	Copyright (c) 1993						**
**	Scott A. Smith							**
**	University of California, Santa Barbara				**
**	Department of Chemistry						**
**	Santa Barbara CA. 93106 USA					**
**								 	**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** Description								**
**									**
** This module deals with chemical shift anisotropy relaxation, mostly	**
** with respect to the treatment by Wangness, Bloch, and Redfield.	**
**									**
*************************************************************************/

///Chapter CSA-CSA Relaxation
///Section Overview
///Body The ...
///Section Available CSA Relaxation Functions

#ifndef   Relax_CSA_h_			// Is this file already included?
#  define Relax_CSA_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <BWRRelax/relaxNMR.h>		// Include WBR relaxation
#include <BWRRelax/relaxRF.h>		// Include WBR relaxation with rf
#include <HSLib/HSham.h>		// Include common Hamiltonians

// ____________________________________________________________________________
// A			CSA-CSA Relaxation Superoperators
// ____________________________________________________________________________


MSVCDLL void RCC(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                            double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


MSVCDLL super_op RCC(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


// ------------- Dipole-Dipole, With An Applied RF-Field ----------------

MSVCDLL void RCCrf(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			Wrflab: Field freqency in the lab frame (Hz)
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


MSVCDLL super_op RCCrf(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type=0, int level=4);

	// Input                sys   : Dynamic spin system
	//                      Heff  : Effective Hamiltonian (Hertz)
	//                      Wrf   : RF-Field frequecy, lab frame (Hertz)
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff


// ----- CSA-CSA, No Applied RF-Field, Dynamic Frequency Shift  ------


MSVCDLL void RCCds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho



MSVCDLL super_op RCCds(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


// ---- CSA-CSA, With An Applied RF-Field, Dynamic Frequency Shift  -----

MSVCDLL void RCCrfds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			Wrflab: Field freqency in the lab frame (Hz)
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


MSVCDLL super_op RCCrfds(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type=0, int level=4);

	// Input                sys   : Dynamic spin system
	//                      Heff  : Effective Hamiltonian (Hertz)
	//                      Wrf   : RF-Field frequecy, lab frame (Hertz)
	//                      phi   : Field phase angle (Degrees)
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff

// ______________________________________________________________________
// ********************* CLASSICAL CSA RELAXATION ***********************
// ______________________________________________________________________

// ------------------ Longitudinal Relaxation, T1 -----------------------

// Farrar & Becker, "Pulse and Fourier Transform NMR", Academic Press, New York, 1971 
//                         Eq. 4.28 page 59 and Eq. 4.9 page 51.
//
//
//  I    1            2    2         2  tau [       2         ]
// R  = --- = (gamma ) (B ) (s  - s )   ___ | _______________ |
//  1    T          I    0    ||   |        |         2     2 |
//        1                   ||  ---   15  | 1 + (w ) (tau)  |
//             				    [	    I         ]
//
//
//                    2     I  2   tau  [       3         ]
//          = (Omega )  (del  )   _____ | _______________ |
//                  I       zz          |       2       2 |
//                                 10   | 1 + (w ) (tau)  |
//             				[       I         ]
//             
//                                3
// where  (s  - s ) =  /\ sigma = _ del  , the first two are what
//          ||   |    /__\        2    zz
//          ||  ---
//
// is commonly called the CSA, the latter what is stored in GAMMA!


MSVCDLL row_vector R1_CC(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		R1 : R1 values for all system spins relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons


MSVCDLL double R1_CC(const sys_dynamic& sys, int I);

	// Input 	       sys : Spin system
	// 			 I : Spin I
	// Output		R1 : R1 value for spin I relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons


MSVCDLL double R1_CC_max(const sys_dynamic& sys, int i);

        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            R1max : Maximum R1 value due to CSA relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.


MSVCDLL double R1_CC_max(const sys_dynamic& sys, const std::string& Iso);

        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            R1max : Maximum R1 value due to CSA relaxation
        //                           over all spins in the system of type Iso
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.


MSVCDLL double R1_CC_max(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		R1 : Maximum of all R1 values of the system spins
	//			     relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


MSVCDLL row_vector T1_CC(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output	       T1s : T1 values for all system spins relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons


MSVCDLL double T1_CC(const sys_dynamic& sys, int I);

	// Input 	       sys : Spin system
	// 			 I : Spin I
	// Output		T1 : T1 value for spin I relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons


MSVCDLL double T1_CC_max(const sys_dynamic& sys, int i);

        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            T1max : Maximum T1 value due to CSA relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.


MSVCDLL double T1_CC_max(const sys_dynamic& sys, const std::string& Iso);

        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            T1max : Maximum T1 value due to CSA relaxation
        //                           over all spins in the system of type Iso
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.


MSVCDLL double T1_CC_max(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		T1 : Maximum of all T1 values of the system spins
	//			     relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


// ------------------- Transverse Relaxation, T2 ------------------------


// Farrar & Becker, "Pulse and Fourier Transform NMR", Academic Press, New York, 1971 
//                               Eq. 4.29 page 59.
//
//
//  I    1            2    2         2  tau [       6              ]
// R  = --- = (gamma ) (H ) (s  - s )   ___ | _______________  + 8 |
//  2    T          I    0    ||   |        |         2     2      |
//        2                   ||  ---   90  | 1 + (w ) (tau)       |
//             				    [	    I              ]
//
//
//                    2     I  2   tau  [       3              ]
//          = (Omega )  (del  )   _____ | _______________  + 4 |
//                  I       zz          |       2       2      |
//                                 20   | 1 + (w ) (tau)       |
//             				[       I              ]
//             
//

MSVCDLL row_vector R2_CC(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		R2 : R2 values for system spins relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons


MSVCDLL double R2_CC(const sys_dynamic& sys, int I);

	// Input 	       sys : Spin system
	// 			 I : Spin I
	// Output		R2 : R2 value for spin I relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons


MSVCDLL double R2_CC_max(const sys_dynamic& sys, int i);

        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            R2max : Maximum R2 value due to CSA relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.


MSVCDLL double R2_CC_max(const sys_dynamic& sys, const std::string& Iso);

        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            R2max : Maximum R2 value due to CSA relaxation
        //                           over all spins in the system of type Iso
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.


MSVCDLL double R2_CC_max(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		R2 : Maximum of all R2 values of the system spins
	//			     relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


MSVCDLL row_vector T2_CC(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		T2 : T2 values for system spins relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons


MSVCDLL double T2_CC(const sys_dynamic& sys, int I);

	// Input 	       sys : Spin system
	// 			 I : Spin I
	// Output		T2 : T2 value for spin I relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons


MSVCDLL double T2_CC_max(const sys_dynamic& sys, int i);

        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            T2max : Maximum T2 value due to CSA relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.


MSVCDLL double T2_CC_max(const sys_dynamic& sys, const std::string& Iso);

        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            T2max : Maximum T2 value due to CSA relaxation
        //                           over all spins in the system of type Iso
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.


MSVCDLL double T2_CC_max(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		T2 : Maximum of all T2 values of the system spins
	//			     relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


// ------------------------- CSA Linewidths -----------------------------
//
//                                 CC
//                                R
//                           CC    2      1.0
//                         LW   = --- = -------
//                           hh    pi    CC 
//                                     	T  * pi
//                                       2
//

MSVCDLL row_vector LWhh_CC(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output	       LWs : Expected linewidths at half-height
	//			     for system spins relaxed by CSA


MSVCDLL double LWhh_CC(const sys_dynamic& sys, int I);

	// Input 	       sys : Spin system
	// 			 I : Spin I
	// Output		LW : Expected linewidth at half-height
	//			     for spin I relaxed by CSA


MSVCDLL double LWhh_CC_max(const sys_dynamic& sys, int i);

        // Input               sys : Spin system
        //                       i : Spin index
        // Output               LW : Maximum expected linewidth at
        //                           half-height of all spins of the
        //                           same type i, relaxed by CSA


MSVCDLL double LWhh_CC_max(const sys_dynamic& sys, const std::string& Iso);

        // Input               sys : Spin system
        //                     Iso : Spin isotope type
        // Output               LW : Maximum expected linewidth at
        //                           half-height of all spins of isotope
        //                           type Iso relaxed by CSA


MSVCDLL double LWhh_CC_max(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		LW : Maximum expected linewidth at half-height
	//			     of all spins relaxed by CSA in the system


// ______________________________________________________________________
// *************** CSA-CSA Relaxation Auxiliary Functions ***************
// ______________________________________________________________________


MSVCDLL matrix xiCSA(const sys_dynamic& dsys);

	// Input		dsys  : A dynamic system
	// Return		xis   : A matrix of CSA interaction
	//				constants (rad/sec)
	// Note			      : The PPM units of delzz are
	//				counteracted by the MHz units
	//				on Omega
	//
	//		 1/2
	//   CSA   [6*pi]
	// xi    = |----| * gamma * B * del  (i) = K * Omega  * del  (i)
	//   i     [ 5  ]        i   0     zz               i      zz


MSVCDLL double xiCSA(const sys_dynamic& dsys, int i);

	// Input		dsys  : A dynamic system
	//			i     : Spin index
	// Return		xii   : The gamma CSA interaction
	//				constant (rad/sec)
	// Note			      : The PPM units of delzz are
	//				counteracted by the MHz units
	//				on Omega
	//
	//		 1/2
	//   CSA   [6*pi]
	// xi    = |----| * gamma * B * del  (i) = K * Omega  * del  (i)
	//   i     [ 5  ]        i   0     zz               i      zz


MSVCDLL matrix xiCSA(const spin_system& sys, double* CSAs);

	// Input		dsys  : A dynamic system
	//			CSAs  : Vector of CSA values (PPM)
	// Return		xis   : A matrix of CSA interaction
	//				constants (rad/sec)
	// Note			      : The PPM units of CSA's are
	//				counteracted by the MHz units
	//				on Omega
	//
	//		 1/2
	//   CSA   [6*pi]               2  /\          2 
	// xi    = |----| * gamma * B * - /__\sigma  = _ * K * Omega  * CSA
	//   i     [ 5  ]        i   0  3          i   3            i      i


MSVCDLL double xiCSA(const spin_system& sys, int i, double csa);

	// Input		dsys  : A dynamic system
	//			i     : Spin index
	//			csa   : CSA value of spin i (PPM)
	// Return		xii   : The gamma CSA interaction
	//				constant (rad/sec)
	// Note			      : The PPM units of CSA is
	//				counteracted by the MHz units
	//				on Omega
	//
	//		 1/2
	//   CSA   [6*pi]               2  /\          2 
	// xi    = |----| * gamma * B * - /__\sigma  = _ * K * Omega  * CSA
	//   i     [ 5  ]        i   0  3          i   3            i      i


MSVCDLL row_vector CSA(const sys_dynamic& dsys);

	// Input		dsys  : A dynamic system
	// Return		csas  : A vector of CSA values in PPM
	//
	//   /\                2
	//  /__\sigma = CSA  = - del  (i)
	//           i     i   3    zz 


MSVCDLL double CSA(const sys_dynamic& dsys, int i);

	// Input		dsys  : A dynamic system
	// Return		csa   : CSA value of spin i in PPM
	//
	//   /\                2
	//  /__\sigma = CSA  = - del  (i)
	//           i     i   3    zz 


#endif							// relaxCSA.h
