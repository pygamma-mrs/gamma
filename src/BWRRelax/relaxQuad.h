/* relax_Quad.h *********************************-*-c++-*-
**							**
**	               G A M M A			**
**							**
**	NMR Library Quadrupolar Relaxation Functions	**
**							**
**	Interface definition				**
**							**
**	Copyright (c) 1991, 1992, 1993			**
**	Scott Smith					**
**	University of California, Santa Barbara		**
**	Department of Chemistry				**
**	Santa Barbara CA. 93106 USA			**
**							**
**						 	**
**      $Header:
**							**
**      Modifications:					**
**							**
*********************************************************/

/*********************************************************
**							**
** 	Description					**
**							**
** Herein are functions dealing with some of the more	**
** common quadrupolar relaxation phenomena.		**
**							**
*********************************************************/

///Chapter Quadrupolar-Quadrupolar Relaxation
///Section Overview
///Body The ...
///Section Available Quadrupolar Relaxation Functions

#ifndef   RelaxQuad_h_			// Is this file already included?
#  define RelaxQuad_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <BWRRelax/relaxNMR.h>		// Include WBR base module
#include <BWRRelax/relaxRF.h>		// Include relaxation with RF
#include <HSLib/HSham.h>		// Include common MR Hamiltonians

// ______________________________________________________________________
// ********* Quadrupole-Quadrupole Relaxation Superoperators ************
// ______________________________________________________________________

// ----------------- Quad-Quad, No Applied RF-Field ------------------

MSVCDLL void RQQ(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
                     double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Q-Q AC & CC
	//				   + = Q-Q AC
	//				   - = Q-Q CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quadrupolar relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


MSVCDLL super_op RQQ(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = Q-Q AC & CC
	//				   + = Q-Q AC
	//				   - = Q-Q CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quadrupolar relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


// --------------- Quad-Quad, With An Applied RF-Field ------------------

MSVCDLL void RQQrf(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			Wrflab: Field freqency in the lab frame (Hz)
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Quad-Quad AC & CC
	//				   + = Quad-Quad AC
	//				   - = Quad-Quad CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quad relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


MSVCDLL super_op RQQrf(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type=0, int level=4);

	// Input                sys   : Dynamic spin system
	//                      Heff  : Effective Hamiltonian (Hertz)
	//                      Wrf   : RF-Field frequecy, lab frame (Hertz)
	//                      phi   : Field phase angle (Degrees)
	// 			type  : Relaxation type to compute
	//				   0 = Quad-Quad AC & CC
	//				   + = Quad-Quad AC
	//				   - = Quad-Quad CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quad relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff


// ---- Quad-Quad, No Applied RF-Field, Dynamic Frequency Shift  -----


MSVCDLL void RQQds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Q-Q AC & CC
	//				   + = Q-Q AC
	//				   - = Q-Q CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quadrupolar relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho



MSVCDLL super_op RQQds(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = Q-Q AC & CC
	//				   + = Q-Q AC
	//				   - = Q-Q CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quadrupolar relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


// -- Quad-Quad, With An Applied RF-Field, Dynamic Frequency Shift  --


MSVCDLL void RQQrfds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			Wrflab: Field freqency in the lab frame (Hz)
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Quad-Quad AC & CC
	//				   + = Quad-Quad AC
	//				   - = Quad-Quad CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quad relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


MSVCDLL super_op RQQrfds(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type=0, int level=4);

	// Input                sys   : Dynamic spin system
	//                      Heff  : Effective Hamiltonian (Hertz)
	//                      Wrf   : RF-Field frequecy, lab frame (Hertz)
	//                      phi   : Field phase angle (Degrees)
	// 			type  : Relaxation type to compute
	//				   0 = Quad-Quad AC & CC
	//				   + = Quad-Quad AC
	//				   - = Quad-Quad CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quad relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff


// ______________________________________________________________________
// **************** CLASSICAL QUADRUPOLAR RELAXATION ********************
// ______________________________________________________________________

// 		       Quadrupolar Longitudinal Relaxation
// Sudmeier, Anderson, and Frye, "Calculation of Nuclear Spin Relaxation Times",
//        Conc. Magn. Reson., 1990, 2, 197-212: page 201 equation [25]
//
//                                     2 
//       1      3(2I+3)      2 [    eta  ][   2*tau         8*tau     ]
// R  = --- = ----------- QCC  |1 + ---  || ---------- + ------------ |
//  1    T        2            [     3   ]|      2   2         2   2  |
//        1   400I (2I-1)                 [ 1 + w tau    1 + 4w tau   ]
//
// where QCC is the quadrupolar coupling constant in angular frequency units


MSVCDLL row_vector R1_QQ(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		R1 : R1 values for all spins under quad. relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


MSVCDLL double R1_QQ(const sys_dynamic& sys, int i);

	// Input	       sys : Spin system
	// 	 	         i : Spin index
	// Output		R1 : R1 value for spin due to quadrupolar relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


MSVCDLL row_vector T1_QQ(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		T1 : T1 values for all spins under quad. relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


MSVCDLL double R1_QQ_max(const sys_dynamic& sys, int i);

        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            R1max : Maximum R1 value due to quad relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.


MSVCDLL double R1_QQ_max(const sys_dynamic& sys, const std::string& Iso);

        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            R1max : Maximum R1 value due to quad relaxation
        //                           over all spins in the system of type Iso
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.


MSVCDLL double R1_QQ_max(const sys_dynamic& sys);

        // Input               sys : Spin system
        // Output               R1 : Maximum of all R1 values of the system spins
        //                           relaxed by a quadrupolar interaction
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.


MSVCDLL row_vector T1_QQ(const sys_dynamic& sys);

        // Input               sys : Spin system
        // Output              T1s : T1 values for all system spins relaxed by quad
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
        // Note                    : This routine does not use the quadrupolar interaction
        //                           constant as it is used for comparisons


MSVCDLL double T1_QQ(const sys_dynamic& sys, int i);

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		T1 : T1 value for spin i under quad. relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


MSVCDLL double T1_QQ_max(const sys_dynamic& sys, int i);
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            T1max : Maximum T1 value due to quad relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.

 
MSVCDLL double T1_QQ_max(const sys_dynamic& sys, const std::string& Iso);

        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            T1max : Maximum T1 value due to quad relaxation
        //                           over all spins in the system of type Iso
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
 
MSVCDLL double T1_QQ_max(const sys_dynamic& sys);
 
        // Input               sys : Spin system
        // Output               T1 : Maximum of all T1 values of the system spins
        //                           relaxed by quadrupolar interacions
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.

 
// ------------------- Transverse Relaxation, T2 ------------------------
 
// 		       Quadrupolar Tranverse  Relaxation
// Sudmeier, Anderson, and Frye, "Calculation of Nuclear Spin Relaxation Times",
//        Conc. Magn. Reson., 1990, 2, 197-212: page 201 equation [26]
//
//                                     2 
//       1    3(2I+3)*tau    2 [    eta  ][          5            2       ]
// R  = --- = ----------- QCC  |1 + ---  || 3 + ---------- + ------------ |
//  2    T        2            [     3   ]|          2   2         2   2  |
//        2   400I (2I-1)                 [     1 + w tau    1 + 4w tau   ]
//
// where QCC is the quadrupolar coupling constant in angular frequency units


MSVCDLL row_vector R2_QQ(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		R2 : R2 values for spins via quadrupolar relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


MSVCDLL double R2_QQ(const sys_dynamic& sys, int i);

	// Input	       sys : Spin system
	// 	 	         i : Spin index
	// Output		R2 : R2 value for spin due to quadrupolar relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


MSVCDLL double R2_QQ(const sys_dynamic& sys, int I);
 
        // Input               sys : Spin system
        //                       I : Spin I
        // Output               R2 : R2 value for spin I relaxed by quad
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
        // Note                    : This routine does not use the quad interaction
        //                           constant as it is used for comparisons
 
 
MSVCDLL double R2_QQ_max(const sys_dynamic& sys, int i);
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            R2max : Maximum R2 value due to quad relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 

MSVCDLL double R2_QQ_max(const sys_dynamic& sys, const std::string& Iso);

        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            R2max : Maximum R2 value due to quad relaxation
        //                           over all spins in the system of type Iso
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and


MSVCDLL row_vector T2_QQ(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		T2 : T2 values for spins via quadrupolar relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


MSVCDLL double T2_QQ(const sys_dynamic& sys, int i);

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		T2 : T2 value for spini via quadrupolar relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

 
MSVCDLL double T2_QQ_max(const sys_dynamic& sys, int i);
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            T2max : Maximum T2 value due to quad relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 

MSVCDLL double T2_QQ_max(const sys_dynamic& sys, const std::string& Iso);

        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            T2max : Maximum T2 value due to quad relaxation
        //                           over all spins in the system of type Iso
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
 
MSVCDLL double T2_QQ_max(const sys_dynamic& sys);
 
        // Input               sys : Spin system
        // Output               T2 : Maximum of all T2 values of the system spins
        //                           relaxed by quad
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
 
// --------------------- Quadrupolar Linewidths -------------------------
//
// 		     LineWidth via Quadrupolar T2 Relaxation

//                                 QQ
//                                R
//                           QQ    2      1.0
//                         LW   = --- = -------
//                           hh    pi    QQ 
//                                     	T  * pi
//                                       2
//

MSVCDLL row_vector LWhh_QQ(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output	       LWs : Expected linewidths at half-height
	//			     for system spins under quad relaxtion


MSVCDLL double LWhh_QQ(const sys_dynamic& sys, int i);

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		LW : Expected linewidth at half-height
	//			     for spin i under quad relaxtion

 
MSVCDLL double LWhh_QQ_max(const sys_dynamic& sys, int i);

        // Input               sys : Spin system
        //                       i : Spin index
        // Output               LW : Maximum expected linewidth at
        //                           half-height of all spins of the
        //                           same type i, relaxed by quad
 
 
MSVCDLL double LWhh_QQ_max(const sys_dynamic& sys, const std::string& Iso);
 
        // Input               sys : Spin system
        //                     Iso : Spin isotope type
        // Output               LW : Maximum expected linewidth at
        //                           half-height of all spins of isotope
        //                           type Iso relaxed by quad
 

MSVCDLL double LWhh_QQ_max(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		LW : Expected linewidth at half-height
	//			     for spin i under quad relaxtion

// ______________________________________________________________________
// ********** Quadrupole-Quadrupole Relaxation Auxiliary Functions ******
// ______________________________________________________________________


MSVCDLL matrix xiQ(const sys_dynamic& sys);

	// Input		dsys  : A dynamic system (this)
	// Return		dximx : A matrix of quadrupolar interaction
	//				constants (rad/sec)
	//
	//		 1/2   QCC   	        1/2  del  (i)
	//   Q     [6*pi]         i       [6*pi]        zz
	// xi    = |----| * ----------  = |----|  * ----------
	//   i     [ 5  ]   2I (2I -1)    [ 5  ]    2I (2I -1)
	//                    i   i                   i   i


MSVCDLL double xiQ(const sys_dynamic& sys, int i);

	// Input		dsys  : A dynamic system (this)
	//			i     : Spin index
	// Return		xii   : Quadrupolar interaction constant (rad/sec)
	//
	//		 1/2   QCC   	        1/2  del  (i)
	//   Q     [6*pi]         i       [6*pi]        zz
	// xi    = |----| * ----------  = |----|  * ----------
	//   i     [ 5  ]   2I (2I -1)    [ 5  ]    2I (2I -1)
	//                    i   i                   i   i


#endif						// relaxQuad.h

