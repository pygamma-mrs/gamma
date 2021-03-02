/* relax_Dip.h **********************************-*-c++-*-
**							**
**	               G A M M A			**
**							**
**	NMR Library Dipolar Relaxation Functions	**
**							**
**	Interface definition				**
**							**
**	Copyright (c) 1991, 1992, 1993			**
**	Scott Smith					**
**	Eidgenoessische Technische Hochschule		**
**	Labor fuer physikalische Chemie			**
**	8092 Zuerich / Switzerland			**
**							**
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
** common dipolar relaxation phenomena.			**
**							**
*********************************************************/

///Chapter Dipole-Dipole Relaxation
///Section Overview
///Body The ...
///Section Available Dipolar Relaxation Functions

#ifndef   Relax_Dip_h_			// Is this file already included?
#  define Relax_Dip_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

//#include <BWRRelax/relaxNMR.h>		// Include Base WBR module
//#include <BWRRelax/relaxRF.h>		// Include Relaxation with RF
//#include <HSLib/HSham.h>		// Include common Hamiltonians
#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <LSLib/SuperOp.h>		 // Include superoperators
#include <LSLib/sys_dynamic.h>		// Include dynamic spin systems

// ______________________________________________________________________
// A        Dipole-Dipole Relaxation Superoperators, No RF Field
// ______________________________________________________________________

MSVCDLL void RDD(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
                      double* taus, double chi, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note				Computed in the eigenbasis of Ho


MSVCDLL super_op RDD(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note				Computed in the eigenbasis of Ho


MSVCDLL super_op RDD(const spin_system& sys, gen_op& Ho, double tau,
				 matrix& dist, int type=0, int level=4);

	// Input		sys   : Spin system
	//			Ho    : General operator
	// 			tau   : Correlation time (in seconds)
	// 			dist  : Distance matrix (in meters)
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	//				computed in the basis of Ho


MSVCDLL super_op RDD_Jgen(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho
	// Note			      : Uses generic spectral density functions


// ______________________________________________________________________
// B        Dipole-Dipole Relaxation Superoperators, With RF Field
// ______________________________________________________________________

MSVCDLL void RDDrf(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         double Wrflab, double* taus, double chi, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note				Computed in the eigenbasis of Ho


MSVCDLL super_op RDDrf(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian (Hertz)
	//			Wrf   : RF-Field frequecy,lab frame (Hertz)
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff


// ______________________________________________________________________
// C Dipole-Dipole Dynamic Frequency Shift Superoperators, No RF Field
// ______________________________________________________________________


MSVCDLL void RDDds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
                          double* taus, double chi, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note				Computed in the eigenbasis of Ho


MSVCDLL super_op RDDds(const sys_dynamic& sys, gen_op& Ho, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note				Computed in the eigenbasis of Ho


// ______________________________________________________________________
// D Dipole-Dipole Dynamic Frequency Shift Superoperators, With RF Field
// ______________________________________________________________________


MSVCDLL void RDDrfds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
                 double Wrflab, double* taus, double chi, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note				Computed in the eigenbasis of Heff


MSVCDLL super_op RDDrfds(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type=0, int level=4);

	// Input		sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian (Hertz)
	//			Wrf   : RF-Field frequecy, lab frame (Hertz)
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff

// ______________________________________________________________________
// E              Dipole-Dipole Interaction & Coupling Constants
// ______________________________________________________________________

/* In GAMMA the dipolar interaction constant is defined to be

	  		          1/2
	  	            [6*pi]      mu.
	  	       -2 * |----|    * --- * hbar * gamma  * gamma
	  	  D         [ 5  ]      4pi               i        j
	  	xi   = _____________________________________________
	  	  ij		         3
	  			        r
	  			         ij                            */

	// Input		cvec  : A coordinate vector(this)
	// 			cutoff: Cutoff value to zero xi's
	// Input		sys   : Spin system
	//       		dsys  : A dynamic spin system
	// 			i,j   : Spin indices
	// 		        dist  : Distances between spins
	// Input		gi    : Gamma of 1st spin
	// Input		gj    : Gamma of 2nd spin
	// Input		rij   : Distance between spins
	// 			angs  : Flag for angstrom distances
	// Return		xi    : Dipolar interaction constant
	//				for spin pair i,j (rad/sec)
	//					  -1      -2
	// Note			      : 1T = 1 J-C  -sec-m	

MSVCDLL matrix xiD(const sys_dynamic& dsys, double cutoff=0.0);
MSVCDLL double xiD(const sys_dynamic& dsys, int i, int j);
MSVCDLL matrix xiD(const spin_system& sys, matrix& dist, int angs=0, double cutoff=0.0);
MSVCDLL double xiD(double gi, double gj, double rij, int angs=0);

/* In GAMMA the dipolar coupling constant is defined to be
 
	  	       mu.
	  	       --- * hbar * gamma  * gamma
	  	       4pi               i        j
	       DCC   = ____________________________
	  	  ij		  3
	  			 r
	  			  ij                                    */

	// Input		cvec  : A coordinate vector(this)
	// 			angs  : Flag for angstrom distances
	// Input		dsys  : A dynamic spin system
	// 			i,j   : Spin indices
	// Input		gi    : Gamma of 1st spin
	// Input		gj    : Gamma of 2nd spin
	// Input		rij   : Distance between spins
	// 			angs  : Flag for angstrom distances
	// Return		DCCmx : A matrix of dipolar coupling
	//				constants (unit rad/sec)
	// Return		DCCij : Dipolar coupling constant
	//				for spin pair i,j (rad/sec)
	// Return		DCCij : Dipolar coupling constant
	//				(unit rad/sec)
	//					  -1      -2
	// Note			      : 1T = 1 J-C  -sec-m

MSVCDLL matrix DCC(const sys_dynamic& dsys);
MSVCDLL double DCC(const sys_dynamic& dsys, int i, int j);
MSVCDLL double DCC(double gi, double gj, double rij, int angs=0);


// ____________________________________________________________________________
// ************************** CLASSICAL RELAXATION ****************************
// ____________________________________________________________________________

// ____________________________________________________________________________
// F                 Longitudinal Relaxation, T1
// ____________________________________________________________________________

//  Farrar & Becker, "Pulse & F.T. NMR", Academic Press, New York, 1971 
//
//  For Spin I Relaxed by Unlike Spin S:
//
//  DDU   1     2   2      2        [ 1               3          3            ]
// R   = --- = g * g * hbar * S(S+1)| __ J (w - w ) + _ J (w ) + _ J (w + w ) |
//  1     T     I   S               [ 12  0  I   S    2  1  I    4  2  I   S  ]
//	   1
//
//              2   2      2          tau  [         2               6                  12      ]
//           = g * g * hbar  * S(S+1) ____ | _________________ + __________ + _________________ |
//              I   S                    6 |             2   2        2   2               2   2 |
//                                    15r  | 1 + (w - w ) tau    1 + w tau    1 + (w + w ) tau  |
//             		         	   [	   I   S              I             I   S       ]
//
//
//  For Spin I Relaxed by Like Spin:
//
//
//  DDL   1    3  4      2         [                ]
// R   = --- = _ g * hbar * I(I+1) | J (w) + J (2w) |
//  1     T    2                   [  1       2     ]
//         1
//
//                  4      2         tau [     1              4      ]
//           = 2 * g * hbar * I(I+1) ___ | __________ +  ___________ |
//                                     6 |      2   2          2   2 |
//                                   5r  [ 1 + w tau     1 + 4w tau  ]
//             
//
// FOR SI UNITS (WHICH GAMMA USES) WE NEED 2 mu /(4*pi) FACTORS ALSO IN FRONT!
//                                             0
//                                                 -1      -2
//                          RECALL THAT 1 T = 1 J C   sec m
//

	// Input	       sys : Spin system
	//		         i : Spin index
	//			 j : Spin index
	// Output		R1 : R1 values for spins in the system
	//			     via dipolar relaxation.
	// Output		R1 : R1 value for spin i dipolar relaxed by all
	//			     other spins in the system
	// Output		R1 : R1 value for spin i dipolar relaxed
	//			     by spin j
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

MSVCDLL row_vector R1_DD(const sys_dynamic& sys);
MSVCDLL double     R1_DD(const sys_dynamic& sys, int i);
MSVCDLL double     R1_DD(const sys_dynamic& sys, int i, int j);

        // Input               sys : Spin system
        //                       i : A spin index
        //                     Iso : A string for spin type 
        // Output            R1max : Maximum R1 value of all spins in
        //                           the system of isotope the same as the
        //                           input spin i
        // Output            R1max : Maximum R1 value of all spins in
        //                           the system of isotope the same as Iso
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.

MSVCDLL double R1_DD_max(const sys_dynamic& sys);
MSVCDLL double R1_DD_max(const sys_dynamic& sys, int i);
MSVCDLL double R1_DD_max(const sys_dynamic& sys, const std::string& Iso);

	// Input	       sys : Spin system
	//		         i : Spin index
	//			 j : Spin index
	// Output		T1 : T1 values for spins in the system
	//			     via dipolar relaxation.
	// Output		T1 : T1 value for spin i dipolar relaxed by all
	//			     other spins in the system
	// Output		T1 : T1 value for spin i dipolar relaxed
	//			     by spin j
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

MSVCDLL row_vector T1_DD(const sys_dynamic& sys);
MSVCDLL double     T1_DD(const sys_dynamic& sys, int i);
MSVCDLL double     T1_DD(const sys_dynamic& sys, int i, int j);

        // Input               sys : Spin system
        //                       i : A spin index
        //                     Iso : A string for spin type
	// Output	     T1max : Maximum T1 value of all spin spins in
	//			     the system under dipolar relaxation    
        // Output            T1max : Maximum T1 value of all spins in
        //                           the system of isotope the same as the
        //                           input spin i
        // Output            T1max : Maximum T1 value of all spins in
        //                           the system of isotope the same as Iso
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.

MSVCDLL double T1_DD_max(const sys_dynamic& sys);
MSVCDLL double T1_DD_max(const sys_dynamic& sys, int i);
MSVCDLL double T1_DD_max(const sys_dynamic& sys, const std::string& Iso);

// ____________________________________________________________________________
// G                 Transverse Relaxation, T2
// ____________________________________________________________________________

//  Farrar & Becker, "Pulse & F.T. NMR", Academic Press, New York, 1971 
//
//  For Spin I Relaxed by Unlike, Uncoupled Spin S
//        Eq. 4.19 page 55 & Eq. 4.9 page 51
//
//  DDU   1     2   2      2         [ 1          1           3          3          3            ]
// R   = --- = g * g * hbar * S(S+1) | _ J (0) + __ J (w  ) + _ J (w ) + _ J (w ) + _ J (w + w ) |
//  2     T     I   S                [ 6  0      24  0  IS    4  1  I    2  1  S    8  2  I   S  ]
//	   2
//
//
//              2   2      2         tau  [          1            3            6               6          ]
//           = g * g * hbar * S(S+1) ____ | 4 + ___________ + __________ + __________ + _________________ |
//              I   S                   6 |          2    2        2   2        2   2               2   2 |
//                                   15r  |     1 + w  tau    1 + w tau    1 + w tau    1 + (w + w ) tau  |
//             			          [	     IS            I            S             I   S       ]
//
//
//  For Spin I Relaxed by Unlike, Coupled Spin S
//
//  DDJ   1     2   2      2         tau  [          1            3              3              6         ]
// R   = --- = g * g * hbar * S(S+1) ____ | 4 + ___________ + __________ + __________ + _________________ |
//  2     T     I   S                   6 |          2    2        2   2        2   2               2   2 |
//         2                         15r  |     1 + w  tau    1 + w tau    1 + w tau    1 + (w + w ) tau  |
//             				  [	     IS            I            S             I           ]
//             
//
//  For Spin I Relaxed by Like Spin
//  Eq. 4.16 page 55,  Eq. 4.9 page 51, and Eq. 4.21 page 56
//
//  DDL   1     4      2         [ 3         15          3         ]
// R   = --- = g * hbar * I(I+1) | _ J (0) + __ J (w ) + _ J (2w ) |
//  2     T                      [ 8  0       4  1  I    8  2   I  ]
//	   2
//
//              4      2         tau [           5            2     ]
//           = g * hbar * I(I+1) ___ | 3 + __________ + ___________ |
//                                 6 |          2   2         2   2 |
//                               5r  [     1 + w tau    1 + 4w tau  ]
//             
//
// FOR SI UNITS (WHICH GAMMA USES) WE NEED 2 mu /(4*pi) FACTORS ALSO IN FRONT!
//
//

	// Input	       sys : Spin system
	//		         i : Spin index
	//			 j : Spin index
	// Output (sys)		R2 : R2 values for spins in the system
	//			     via dipolar relaxation.
	// Output (sys, i)	R2 : R2 value for spin i dipolar relaxed by all
	//			     other spins in the system
	// Output (sys, i, j)	R2 : R2 value for spin i dipolar relaxed
	//			     by spin j
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

MSVCDLL row_vector R2_DD(const sys_dynamic& sys);
MSVCDLL double     R2_DD(const sys_dynamic& sys, int i);
MSVCDLL double     R2_DD(const sys_dynamic& sys, int i, int j);

	// Input	       sys : Spin system
        //                       i : A spin index
        //                     Iso : A string for the spin isotope type
	// Output (sys)	     R2max : Maximum R2 value of all spin spins in
	//			     the system under dipolar relaxation
        // Output (sys,i)    R2max : Maximum R2 (under dipolar relaxation) value
        //                           of all spins in the system of isotope type 
        //                           the same as the input spin i
        // Output (sys,Iso)  R2max : Maximum R2 (under dipolar relaxation) value
        //                           of all spins in the system of isotope type Iso
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

MSVCDLL double R2_DD_max(const sys_dynamic& sys);
MSVCDLL double R2_DD_max(const sys_dynamic& sys, int i); 
MSVCDLL double R2_DD_max(const sys_dynamic& sys, const std::string& Iso); 


	// Input	       sys : Spin system
	// Output		T2 : T2 values for spins in the system
	//			     via dipolar relaxation.
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

	// Input	       sys : Spin system
	//		         i : Spin index
	// Output		T2 : T2 value for spin i dipolar relaxed by all
	//			     other spins in the system
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


MSVCDLL row_vector T2_DD(const sys_dynamic& sys);
MSVCDLL double     T2_DD(const sys_dynamic& sys, int i);
MSVCDLL double     T2_DD(const sys_dynamic& sys, int i, int j);

	// Input	       sys : Spin system
	//		         i : Spin index
	//			 j : Spin index
	// Output		T2 : T2 value for spin i dipolar relaxed
	//			     by spin j
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            T2max : Maximum T2 (under dipolar relaxation) value
        //                           of all spins in the system of isotope type 
        //                           the same as the input spin i
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
        //                           the system moving as an isotropic top.
        // Input               sys : Spin system
        //                     Iso : A string for the spin isotope type
        // Output            T2max : Maximum T2 (under dipolar relaxation) value
        //                           of all spins in the system of isotope type Iso
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
	// Input	       sys : Spin system
	// Output	     T2max : Maximum T2 value of all spin spins in
	//			     the system under dipolar relaxation
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

MSVCDLL double T2_DD_max(const sys_dynamic& sys, int i); 
MSVCDLL double T2_DD_max(const sys_dynamic& sys, const std::string& Iso); 
MSVCDLL double T2_DD_max(const sys_dynamic& sys);

// ____________________________________________________________________________
// H                       Dipolar Linewidths
// ____________________________________________________________________________

//                                 DD
//                                R
//                           DD    2      1.0
//                         LW   = --- = -------
//                           hh    pi    DD 
//                                     	T  * pi
//                                       2
//
	// Input	       sys : Spin system
	//		         i : Spin index
	//			 j : Spin index
	// Output		LW : Expected linewidths at half-height
 	//			     for spins in the system via dipolar relaxation.
	// Output		LW : Expected linewidth at half-height
	//			     for spin i dipolar relaxed by the system
	// Output		LW : Expected linewidth at half-height
	//			     for spin i dipolar relaxed by spin j
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

MSVCDLL row_vector LWhh_DD(const sys_dynamic& sys);
MSVCDLL     double LWhh_DD(const sys_dynamic& sys, int i);
MSVCDLL     double LWhh_DD(const sys_dynamic& sys, int i, int j);

        // Input               sys : Spin system
        //                       i : A spin index       
        // Output               LW : Maximum linewidth at half-height
        //                           under dipolar relaxation over all
        //                           spins in the system of isotope the
        //                           same as the input spin i
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and

        // Input               sys : Spin system
        //                     Iso : A string for isotope type
        // Output               LW : Maximum linewidth at half-height
        //                           under dipolar relaxation over all
        //                           spins in the system of isotope of type Iso
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and

	// Input 	       sys : Spin system
	// Output		LW : Maximum linewidth at half-height
	//			     over all spins under dipolar relaxation

MSVCDLL double LWhh_DD_max(const sys_dynamic& sys);
MSVCDLL double LWhh_DD_max(const sys_dynamic& sys, int i);
MSVCDLL double LWhh_DD_max(const sys_dynamic& sys, const std::string& Iso);

// ____________________________________________________________________________
// I                    Nuclear Overhauser Enhancements
// ____________________________________________________________________________

//  NOE Enhancement For Spin Relaxed by Spin S

//                [    -1                 6         ]
//                | ___________ + _________________ |
//             g  |      2    2               2   2 |
//              I | 1 + w  tau    1 + (w + w ) tau  |
//                [	   IS           S   I       ]		   
// NOE = ___________________________________________________  = rho
//          [     1             3                6         ]       NOE
//          | ___________ + __________ + _________________ |
//       g  |      2    2        2   2               2   2 |
//        S | 1 + w  tau    1 + w tau    1 + (w + w ) tau  |
//          [	   IS            I             I   S       ]
//             
//             
//                       eta    = 1 + rho
//                          NOE          NOE


MSVCDLL double NOE(const sys_dynamic& sys, int i, int j, double eta=0);

	// Input	       sys : Spin system
	//		         i : Spin index
	//			 j : Spin index
	//		       eta : flag for % enhancement
	// Output	       NOE : Expected NOE for for spin i due to
	//			     dipolar relaxation by spin j
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


// ____________________________________________________________________________
// I               Multiple Quantum Transition Relaxation Rates
// ____________________________________________________________________________

//  ZQT - For Spin I Relaxed by Unlike, Coupled Spin S
//
//  DDJ   1     2   2      2         tau  [     2             3            3      ]
// R   = --- = g * g * hbar * S(S+1) ____ | ___________ + __________ + __________ |
//  2     T     I   S                   6 |      2    2        2   2        2   2 |
//         2                         15r  | 1 + w  tau    1 + w tau    1 + w tau  |
//             				  [	 IS            I            S     ]
//             
//
//  SQT - For Spin I Relaxed by Unlike, Coupled Spin S
//
//  DDJ   1     2   2      2         tau  [          1            3              3              6         ]
// R   = --- = g * g * hbar * S(S+1) ____ | 4 + ___________ + __________ + __________ + _________________ |
//  2     T     I   S                   6 |          2    2        2   2        2   2               2   2 |
//         2                         15r  |     1 + w  tau    1 + w tau    1 + w tau    1 + (w + w ) tau  |
//             				  [	     IS            I            S             I           ]
//             
//
//  DQT - For Spin I Relaxed by Unlike, Coupled Spin S
//
//  DDJ   1     2   2      2         tau  [     3              3              12        ]
// R   = --- = g * g * hbar * S(S+1) ____ | __________ + __________ + _________________ |
//  2     T     I   S                   6 |      2   2        2   2               2   2 |
//         2                         15r  | 1 + w tau    1 + w tau    1 + (w + w ) tau  |
//             				  [	 I            S             I   S	]
//             
//
// FOR SI UNITS (WHICH GAMMA USES) WE NEED 2 mu /(4*pi) FACTORS ALSO IN FRONT!
//
//


	// Input	       sys : Spin system
	//		       MQC : Multiple quantum coherence
	//		         i : Spin index
	//			 j : Spin index
	// Output (sys,MQC)	R2 : R2 values for spins in the system
	//			     via dipolar relaxation.
	// Output (sys,MQC,i)	R2 : R2 value for spin i dipolar relaxed by all
	//			     other spins in the system
	// Output (sys,MQC,i,j)	R2 : R2 value for MQC involving spin i dipolar relaxed
	//			     by spin j
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

MSVCDLL row_vector R2_DDMQT(const sys_dynamic& sys, int MQC);
MSVCDLL double     R2_DDMQT(const sys_dynamic& sys, int MQC, int i);
MSVCDLL double     R2_DDMQT(const sys_dynamic& sys, int MQC, int i, int j);


#endif						// relaxDip.h

