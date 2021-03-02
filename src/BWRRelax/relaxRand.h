/* relaxRand.h **********************************-*-c++-*-
**							**
**	               G A M M A			**
**							**
**	NMR Library Random Field Relaxation Functions	**
**							**
**	Interface definition				**
**							**
**	Copyright (c) 1993				**
**	Scott A. Smith					**
**	University of California, Santa Barbara		**
**	Department of Chemistry				**
**	Santa Barbara CA. 93106 USA			**
**							**
**						 	**
**      $Header: $
**							**
*********************************************************/

/*********************************************************
**							**
** 	Description					**
**							**
** Herein are functions dealing with some of the more	**
** common Random Field relaxation phenomena.		**
**							**
*********************************************************/

///Chapter Random Field Relaxation
///Section Overview
///Body The ...
///Section Available Random Field Relaxation Functions

#ifndef   Relax_Random_h_		// Is this file already included?
#  define Relax_Random_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <BWRRelax/relaxNMR.h>
#include <HSLib/HSham.h>


// ______________________________________________________________________
// ****************** RDM-RDM Relaxation Superoperators *****************
// ______________________________________________________________________


MSVCDLL void RRRx(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                                     double tau, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			tau   : System random field correlation times
	// 			type  : Relaxation type to compute
	//				   0 = RDM-RDM AC & CC
	//				   + = RDM-RDM AC
	//				   - = RDM-RDM CC
	//			level : Relaxation treatment level
	// Output		LOp   : Random Field relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho
	// Note			      :	Uses a spherical top spectral density


MSVCDLL super_op RRRx(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = RDM-RDM AC & CC
	//				   + = RDM-RDM AC
	//				   - = RDM-RDM CC
	//			level : Relaxation treatment level
	// Output		LOp   : Random Field relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho
	// Note			      :	Uses a spherical top spectral density

MSVCDLL void RRR(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = RDM-RDM AC & CC
	//				   + = RDM-RDM AC
	//				   - = RDM-RDM CC
	//			level : Relaxation treatment level
	// Output		LOp   : RDM relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho



MSVCDLL super_op RRR(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = RDM-RDM AC & CC
	//				   + = RDM-RDM AC
	//				   - = RDM-RDM CC
	//			level : Relaxation treatment level
	// Output		LOp   : RDM relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


MSVCDLL void Rij_rdm(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			xi1s  : Interaction constant, mu1
	//			xi2s  : Interaction constant, mu2
	//			A1    : Spatial tensor, mu1
	//			A2    : Spatial tensor, mu2
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Auto & Cross Correlation
	//				   + = Auto Correlation Only
	//				   - = Cross Correlation Only
	//			level : Relaxation treatment level

/*	              Two Single Spin Mechanisms

                      --- ---             --- ---
                      \   \               \   \
               LOp += /   /    R        = /   /   R
  	              --- ---   mu1,mu2   --- ---  i,j
                      mu1 mu2              i   j                       */

MSVCDLL void Rmumu_rdm(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
               double w0, double w1, double w2, int level=4, int autoc=0);

// ______________________________________________________________________
// ***************** CLASSICAL RANDOM FIELD RELAXATION ******************
// ______________________________________________________________________

// ------------------ Longitudinal Relaxation, T1 -----------------------


MSVCDLL row_vector R1_RR(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		R1 : R1 values for all system spins relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons


MSVCDLL double R1_RR(const sys_dynamic& sys, int i);

	// Input 	       sys : Spin system
	// 			 i : Spin i
	// Output		R1 : R1 value for spin i relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons


MSVCDLL double R1_RR_max(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		R1 : Maximum of all R1 values of the system spins
	//			     relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


MSVCDLL row_vector T1_RR(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output	       T1s : T1 values for all system spins relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons


MSVCDLL double T1_RR(const sys_dynamic& sys, int i);

	// Input 	       sys : Spin system
	// 			 i : Spin i
	// Output		T1 : T1 value for spin i relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons


MSVCDLL double T1_RR_max(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		T1 : Maximum of all T1 values of the system spins
	//			     relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


// ------------------- Transverse Relaxation, T2 ------------------------


MSVCDLL row_vector R2_RR(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		R2 : R2 values for system spins relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons


MSVCDLL double R2_RR(const sys_dynamic& sys, int i);

	// Input 	       sys : Spin system
	// 			 i : Spin i
	// Output		R2 : R2 value for spin i relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons


MSVCDLL double R2_RR_max(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		R2 : Maximum of all R2 values of the system spins
	//			     relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


MSVCDLL row_vector T2_RR(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		T2 : T2 values for system spins relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons


MSVCDLL double T2_RR(const sys_dynamic& sys, int i);

	// Input 	       sys : Spin system
	// 			 i : Spin i
	// Output		T2 : T2 value for spin i relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons


MSVCDLL double T2_RR_max(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		T2 : Maximum of all T2 values of the system spins
	//			     relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


// ------------------------- R.F. Linewidths -----------------------------
//
//                                 RR
//                                R
//                           RR    2      1.0
//                         LW   = --- = -------
//                           hh    pi    RR 
//                                     	T  * pi
//                                       2
//

MSVCDLL row_vector LWhh_RR(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output	       LWs : Expected linewidths at half-height
	//			     for system spins relaxed by R.F.


MSVCDLL double LWhh_RR(const sys_dynamic& sys, int i);

	// Input 	       sys : Spin system
	// 			 i : Spin i
	// Output		LW : Expected linewidth at half-height
	//			     for spin i relaxed by R.F.


MSVCDLL double LWhh_RR_max(const sys_dynamic& sys, int i);
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output          LWhhmax : Max. LWhh value from Rand. Field relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
 

MSVCDLL double LWhh_RR_max(const sys_dynamic& sys, const std::string& Iso);
 
        // Input               sys : Spin system
        //                     Iso : A string for an Isotope type
        // Output          LWhhmax : Max. LWhh value due to Rand.Field relaxation
        //                           over all spins in the system of isotope
        //                           type as Iso
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 

MSVCDLL double LWhh_RR_max(const sys_dynamic& sys);

	// Input 	       sys : Spin system
	// Output		LW : Maximum expected linewidth at half-height
	//			     of all spins relaxed by Rand.Field in the system


// ______________________________________________________________________
// *************** RDM-RDM Relaxation Auxiliary Functions ***************
// ______________________________________________________________________


MSVCDLL matrix xiRDM(const sys_dynamic& dsys);

	// Input		dsys  : A dynamic system
	// Return		xis   : A matrix of RDM interaction
	//				constants (rad/sec)
	//


MSVCDLL double xiRDM(const sys_dynamic& dsys, int i);

	// Input		dsys  : A dynamic system
	// Return		xi    : RDM interaction constant (rad/sec)
	//


#endif /* __RELAX_Random_H__ */
