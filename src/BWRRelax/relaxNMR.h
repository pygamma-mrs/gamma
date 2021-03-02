/* relaxNMR.h ***************************************************-*-c++-*-
**									**
**	                         G A M M A				**
**									**
**	NMR Relaxation Functions 		     Interface		**
**									**
**	Copyright (c) 1991, 1992, 1993					**
**	Scott Smith							**
**	Eidgenoessische Technische Hochschule				**
**	Labor fuer physikalische Chemie					**
**	8092 Zuerich / Switzerland					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** Description								**
**									**
** Herein are functions for producing some of the more common		**
** relaxation superoperators.						**
**									**
*************************************************************************/


///Chapter Relaxation Superoperators
///Section Overview
///Body The ...
///Section Available Relaxation Superoperators

#ifndef   RelaxMR_h_			// Is file already included?
#  define RelaxMR_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <BWRRelax/relaxJ.h>		// Include rigid diffusion J
#include <HSLib/SpinSystem.h>		// Include base spin systems
#include <Level1/SpaceT.h>		// Include spatial tensors
#include <Level1/SpinT.h>		// Include spin tensors
#include <LSLib/SuperOp.h>		// Include superoperators
#include <Matrix/matrix.h>		// Include GAMMA matrices
#include <HSLib/GenOp.h>		// Include general operators
#include <LSLib/sys_dynamic.h>		// Include anisotropic systems

// ____________________________________________________________________________
// i               RELAXATION SUPEROPERATOR ERROR HANDLING
// ____________________________________________________________________________

	// Input		eidx	: Error index
	//			noret	: Flag for return (0=linefeed)
	// Output		none 	: Error Message Output
	// Output		none 	: Stops execution (fatal)

         void RlxNMRerror(int eidx, int noret=0);
volatile void RlxNMRfatal(int eidx);

// ____________________________________________________________________________
// ********* RELAXATION SUPEROPERATORS VIA ELEMENT CALCULATIONS *********
// ____________________________________________________________________________

// ---------------- Level 4 via Element Calculation ---------------------

MSVCDLL void R_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12);

	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	// Return		void  : Relaxation superoperator for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R
	//				          out	   in    12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases


//                  <a,a'| LOp |b,b'> += <a,a'| R   |b,b'>
//		                                 1,2

// where a, a', b, b' are basis funciton indices.


MSVCDLL double R_4(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	         int rank, int a, int b, int aa, int bb);

	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//			rank  : Rank of the two interactions
	//			a, b  : 1st transisiton indices
	//			aa, bb: 2nd transisiton indices
	// Output		Rel   : Relaxation superoperator element
	//
	//				    <a,aa|R  |b,bb>
	//					   12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases

/*                    --- [ ---
                      \   | \                 m       m 
   <a,a'|R   |b,b'> = /   | /   delta     <a|T |g><b|T |g> J  (w  )
  	  1,2	      --- | ---      a',b'    1       2     12  gb
  		       m  [  g
  
                                 m        m                       m        m
                           - <a|T |b><a'|T |b'> J  (w    ) - <b'|T |a'><b|T |a> J  (w  )
  	     	                 1        2      12  b'a'         1        2     12  ab
  
                                                     ---                                     ]
                                                     \               m        m              |
                                                  +  /   delta   <g|T |a'><g|T |b'> J  (b'g) |
  	     	                                     ---      a,b    1        2      12      |
                                                      g                                      ]
*/


// ---------------- Level 3 via Element Calculation ---------------------
//                  (Applies Secular Approximation)

MSVCDLL void R_3(super_op& LOp, double* w, int rank, gen_op* T1s, gen_op* T2s,
				          matrix& J12, double cutoff=1.e-2);

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector of system energies
	//				in rad/sec in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//			cutoff: Secular approximation cutoff (Hz)
	// Return		void  : Relaxation superoperator for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R
	//				          out	   in    12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases


//                  <a,a'| LOp |b,b'> += <a,a'| R   |b,b'>
//		                                 1,2

// where a, a', b, b' are basis funciton indices.


// ---------------- Level 2 via Element Calculation ---------------------
//                    (No Degenerate Transitions)

MSVCDLL void R_2(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12);


	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	// Return		void  : Relaxation superoperator for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R
	//				          out	   in    12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases


//                  <a,a'| LOp |b,b'> += <a,a'| R   |b,b'>
//		                                 1,2

// where a, a', b, b' are basis funciton indices.


MSVCDLL double R_2(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	         int rank, int a, int b, int aa, int bb);

	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//			rank  : Rank of the 2 interactions
	//			a, b  : 1st transition indices
	//			aa, bb: 2nd transition indices
	// Output		Rel   : Level 2 relaxation superoperator
	//				element for interactions 1 & 2
	//
	//				    <a,b|R2  |aa,bb>
	//					   12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases

/*			      Level 2 Diagonal Elements
  
                       --- [ ---
                       \   | \       m       m 
   <a,a'|R2   |a,a'> = /   | /   <a|T |g><a|T |g>*J(w  )
  	   1,2	       --- | ---     1       2       ga
  		        m  [  g
  
                                                   m        m
                                       - 2.0 * <a|T |a><a'|T |a'> * J(w  )
  	     	                                   1        2          aa
  
                                                     ---                            ]
                                                     \       m        m             |
                                                  +  /   <g|T |a'><g|T |a'>*J(w   ) |
  	     	                                     ---     1        2        ga'  |
  										    ]
  
  			Level 2 Non-Zero Off Diagonal Elements
  
                     ---
                     \              m       m 
   <a,a|R2   |b,b> = /   -2.0 * <a|T |b><a|T |b> * J  (w  )
  	  1,2	     ---            1       2       12  ba

   where m sums over angular momentum components, g over basis functions,
   and J values include xi1 and xi2
*/


MSVCDLL double Rodiag_2(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	                      int rank, int a, int b);

	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components of
	//				first interaction
	// 			T2s   : Spin tensor components of
	//				second interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//			rank  : Rank of the 2 interactions
	//			a, b  : Transisiton indices
	// Output		Rel   : Relaxation superoperator element
	//				for interactions 1 & 2 for component m
	//
	//				    <a,a|R  |b,b>
	//					  12
	//
	// Note				All operators contained in both
	//				T1s and T2s are assumed to be in
	//				the static Hamiltonian eigenbasis
	// Note				Would be nice if we could use double
	//				rather than complex here!! 


// sosi swithc this to complex
MSVCDLL double Rdiag_2(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	                      int rank, int a, int aa);

	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components of
	//				first interaction
	// 			T2s   : Spin tensor components of
	//				second interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//			rank  : Rank of the 2 interactions
	//			a, aa : Transisiton indices
	// Output		Rel   : Relaxation superoperator element
	//				for interactions 1 & 2 for component m
	//
	//				    <a,a|R  |b,b>
	//					  12
	//
	// Note				All operators contained in both
	//				T1s and T2s are assumed to be in
	//				the static Hamiltonian eigenbasis
	// Note				Would be nice if we could use double
	//				rather than complex here!! 


// ---------------- Level 0 via Element Calculation ---------------------
//			  (Extreme Narrowing)

MSVCDLL void R_0(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, const complex& J12);

	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Reduced spectral density at 0 Hz 
	//				(extreme narrowing) for 2 interactions
	//				times the two interaction constants
	// Return		void  : Relaxation superoperator for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + LOp
	//				          out	   in      12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases


//                  <a,a'| LOp |b,b'> += <a,a'| R   |b,b'>
//		                                 1,2

// where a, a', b, b' are basis funciton indices.  This assumes extreme narrowing.


// sosi switch this to complex
MSVCDLL double R_0(int hs, gen_op* T1s, gen_op* T2s,
		 	         int rank, int a, int b, int aa, int bb);


	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			rank  : Rank of the two interactions
	//			a, b  : 1st transisiton indices
	//			aa, bb: 2nd transisiton indices
	// Output		Rel   : Relaxation superoperator element
	//
	//				    <a,aa|R  |b,bb>
	//					   12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases
	// Note				LOp element return is unitless!


/*                    --- [ ---
                      \   | \                 m       m 
   <a,a'|R   |b,b'> = /   | /   delta     <a|T |g><b|T |g>
  	  1,2	      --- | ---      a',b'    1       2
  		       m  [  g
  
                                 m        m                     m 
                           - <a|T |b><a'|T |b'> - <b'|T |a'><b|T |a>
  	     	                 1        2            1        2
  
                                                     ---                            ]
                                                     \               m        m     |
                                                  +  /   delta   <g|T |a'><g|T |b'> |
  	     	                                     ---      a,b    1        2     |
                                                      g                             ]
*/


// ______________________________________________________________________
// ************ RELAXATION SUPEROPERATORS VIA DOUBLE COMMUTATORS ********
// ______________________________________________________________________


MSVCDLL void R_4s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12);

	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the interactions involved
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	// Return		void  : Relaxation superoperator for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + LOp
	//				          out	   in      12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases

/*                    --- ---
                      \   \       m    m     -m
  	     R(1,2) = /   /   (-1)  [ T , [ T   ,  ] ] J  (w )
  	    	      --- ---          1     2,p        12  p
  		       p   m
*/

// where p sums over transitions, m sums over angular momentum components,
// and J values include xi1 and xi2 (J12


// ------------------ Level 3 via Double Commutators --------------------

//                    (Applies Secular Approximation)

MSVCDLL void R_3s(super_op& LOp, double* w, int rank, gen_op* T1s, gen_op* T2s, matrix& J12);

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector of system energies
	//				in rad/sec in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	// Return		void  : Relaxation superoperator for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + LOp
	//				          out	   in      12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases

/*             --- --- ---
               \   \   \       m                 m       -m
      R(1,2) = /   /   /   (-1)  delta        [ T   , [ T    ,  ] ] J  (w  )
  	       --- --- ---            w , -w     1,p     2,p'        12  p'
  		p   p'  m              p    p'
*/

// where p and p' sum over transitions, m sums over angular momentum components,
// and J values include xi1 and xi2.  The secular approximation restricts the
// sum to contain only terms where the transition frequencies for p & p' cancel.


// ------------------ Level 2 via Double Commutators --------------------

//                      (No Degenerate Transitions)

MSVCDLL void R_2s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12);

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector of system energies
	//				in rad/sec in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	// Return		void  : Relaxation superoperator for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + LOp
	//				          out	   in      12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases

/*                          --- ---
                            \   \       m    m       -m
                   R(1,2) = /   /   (-1)  [ T   , [ T   ,  ] ] J  (w )
  	                    --- ---          1,p     2,p        12  p
  		             p   m
*/

// where p sums over transitions, m sums over angular momentum components,
// and J values include xi1 and xi2.  The Level 2 treatment assumes that there
// are no degenerate transitions in the system.


// ------------------ Level 1 via Double Commutators --------------------

//              (Spin Tensor Component Constant Oscillation)

//void R_1s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12);

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector of system energies
	//				in rad/sec in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	// Return		void  : Relaxation superoperator for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + LOp
	//				          out	   in      12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases

/*                          ---
                            \       m    m     -m
                   R(1,2) = /   (-1)  [ T , [ T  ,  ] ] J  (mw )
  	                    ---          1     2         12   0
  		             m
*/

// m sums over angular momentum components, and J values include xi1 and xi2.
// The Level 1 treatment assumes that there are no degenerate transitions in the system.


// ------------------ Level 0 via Double Commutators --------------------

//                         (Extreme Narrowing)

//void R_0s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, double J12);

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector of system energies
	//				in rad/sec in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Reduced spectral density function value
	//				at 0 Hz times the 2 interaction constants
	// Return		void  : Relaxation superoperator for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + LOp
	//				          out	   in      12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases

/*                                ---
                                  \       m    m     -m
                   R(1,2) = J (0) /   (-1)  [ T , [ T  ,  ] ]
  	                    12    ---          1     2
  		                   m
*/

// m sums over angular momentum components, and J value include xi1 and xi2.
// The Level 1 treatment assumes extreme narrowing.


// ______________________________________________________________________
// ***************** USEFUL DOUBLE COMMUTATOR SUPEROPERATORS ************
// ______________________________________________________________________

// ----------------------- Level 0 Functions ----------------------------

MSVCDLL super_op R_AC_0(spin_T &T);

	// Input		T     : Spin tensor
	// 			LOp1  : Super operator
	//			Op    : General operator
	// 			xisq  : Scaling factor
	// Output		LOp   : Relaxation superoperator

/*                   --- [  ]m [   m  [   -m   ] ]
               LOp = \   |-1|  | T1 , | T1  ,  | |
                     /   [  ]  [   l  [   l    ] ]
  	             ---
                      m
*/


MSVCDLL void R_AC_0(spin_T &T, super_op &LOp1, gen_op &Op, double xisq = 1);

	// Input		T     : Spin tensor
	// 			LOp1  : Super operator
	//			Op    : General operator
	// 			xisq  : Scaling factor
	// Output		None  : LOp has relaxation superoperator
	//				added to it. Computed in the basis
	//				of Op.

/*                    --- [  ]m [   m  [   -m   ] ]
               LOp += \   |-1|  | T1 , | T1  ,  | |
                      /   [  ]  [   l  [   l    ] ]
  	              ---
                       m
*/


MSVCDLL void R_AC_0(gen_op *Ts, super_op &LOp, int rank=2, double xisq=1);

	// Input		Ts    : Spin tensor 1 components
	// 			LOp   : Super operator
	// 			rank  : Interaction rank
	// 			xisq  : Scaling factor
	// Output		None  : LOp has relaxation superoperator
	//				added to it.
	// Note				All spin tensor components in T1s
	//				are assumed in computation basis

/*                    --- [  ]m [   m  [   -m   ] ]
               LOp += \   |-1|  | T1 , | T1  ,  | |
                      /   [  ]  [   l  [   l    ] ]
  	              ---
                       m
*/


MSVCDLL void R_CC_0(spin_T &T1, spin_T &T2, super_op &LOp1,
		                     gen_op &Op, double xisq = 1);

	// Input		T1    : Spin tensor 1
	// 			T2    : Spin tensor 2
	// 			LOp1  : Super operator
	//			Op    : General operator
	// 			xisq  : Scaling factor
	// Output		None  : LOp has relaxation superoperator
	//				added to it. Computed in the basis
	//				of Op.

/*                    --- [  ]m [   m  [   -m   ] ]
               LOp += \   |-1|  | T1 , | T2  ,  | |
                      /   [  ]  [   l  [   l    ] ]
  	              ---
                       m
*/


MSVCDLL void R_CC_0(gen_op *T1s, gen_op* T2s, super_op &LOp1,
		             int rank=2, double xisq = 1);

	// Input		T1s   : Spin tensor 1 components
	// 			T2    : Spin tensor 2 components
	// 			LOp1  : Super operator
	// 			rank  : Interaction rank
	// 			xisq  : Scaling factor
	// Output		None  : LOp has relaxation superoperator
	//				added to it. Computed in the basis
	//				of Op.

//                    --- [  ]m [   m  [   -m   ] ]
//             LOp += \   |-1|  | T1 , | T2  ,  | |
//                    /   [  ]  [   l  [   l    ] ]
//	              ---
//                     m


MSVCDLL void R_CC_0_trans(gen_op *T1s, gen_op* T2s, super_op &LOp1,
		                   int rank=2, double xisq = 1);

	// Input		T1s   : Spin tensor 1 components
	// 			T2    : Spin tensor 2 components
	// 			LOp1  : Super operator
	// 			rank  : Interaction rank
	// 			xisq  : Scaling factor
	// Output		None  : LOp has relaxation superoperator
	//				added to it. Computed in the basis
	//				of Op.
	// Note			None  : It is assumed the vector T2s
	//				contains the transposed operators

//                    --- [   m  [ {  m}t   ] ]
//             LOp += \   | T1 , | |T2 | ,  | |
//                    /   [   l  [ {  l}    ] ]
//	              ---
//                     m


// ----------------------- Level 1 Functions ----------------------------


MSVCDLL void R_AC_1(spin_T &T, super_op &LOp1, gen_op &Op,
		             	double J0, double J1, double J2);

	// Input		T     : Spin tensor
	// 			LOp1  : Super operator
	//			Op    : General operator
	// 			J0    : Spectral density at 0 Hz
	// 			J1    : Spectral density at Larmor frequency
	// 			J2    : Spectral density at 2*Larmor
	// Output		None  : LOp has relaxation superoperator
	//				added to it. Computed in the basis
	//				of Op.

//                    ---        [  ]m [   m  [   -m   ] ]
//             LOp += \   J(mw ) |-1|  | T1 , | T1  ,  | |
//                    /       0  [  ]  [   l  [   l    ] ]
//	              ---
//                     m



MSVCDLL void R_AC_1(gen_op *Ts, super_op &LOp1, int rank,
		             	double J0, double J1, double J2);

	// Input		Ts   : Spin tensor 1 components
	// 			LOp1  : Super operator
	// 			rank  : Interaction rank
	// 			J0    : Spectral density at 0 Hz
	// 			J1    : Spectral density at Larmor frequency
	// 			J2    : Spectral density at 2*Larmor
	// Output		None  : LOp has relaxation superoperator
	//				added to it.
	// Note				All operators (T1s & T2s) should be
	//				in the perfered computation basis

//                    ---        [  ]m [   m  [   -m   ] ]
//             LOp += \   J(mw ) |-1|  | T1 , | T1  ,  | |
//                    /       0  [  ]  [   l  [   l    ] ]
//	              ---
//                     m


MSVCDLL void R_CC_1(spin_T &T1, spin_T &T2, super_op &LOp1,
		        gen_op &Op, double J0, double J1, double J2);

	// Input		T1    : Spin tensor 1
	// 			T2    : Spin tensor 2
	// 			LOp1  : Super operator
	//			Op    : General operator
	// 			J0    : Spectral density at 0 Hz
	// 			J1    : Spectral density at Larmor frequency
	// 			J2    : Spectral density at 2*Larmor
	// Output		None  : LOp has relaxation superoperator
	//				added to it. Computed in the basis
	//				of Op.

//                    ---        [  ]m [   m  [   -m   ] ]
//             LOp += \   J(mw ) |-1|  | T1 , | T2  ,  | |
//                    /       0  [  ]  [   l  [   l    ] ]
//	              ---
//                     m



MSVCDLL void R_CC_1(gen_op *T1s, gen_op* T2s, super_op &LOp1,
		        int rank, double J0, double J1, double J2);

	// Input		T1s   : Spin tensor 1 components
	// 			T2    : Spin tensor 2 components
	// 			LOp1  : Super operator
	// 			rank  : Interaction rank
	// 			J0    : Spectral density at 0 Hz
	// 			J1    : Spectral density at Larmor frequency
	// 			J2    : Spectral density at 2*Larmor
	// Output		None  : LOp has relaxation superoperator
	//				added to it.
	// Note				All operators (T1s & T2s) should be
	//				in the perfered computation basis

//                    ---        [  ]m [   m  [   -m   ] ]
//             LOp += \   J(mw ) |-1|  | T1 , | T2  ,  | |
//                    /       0  [  ]  [   l  [   l    ] ]
//	              ---
//                     m


// ______________________________________________________________________
// ************* RELAXATION SUPEROPERATOR GENERAL FUNCTIONS *************
// ______________________________________________________________________

// ---------------- Routing to Specific Level Computation ---------------

 
MSVCDLL void Rmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
               double w0, double w1, double w2, int l, int level=4, int autoc=0, int het=0);

	// Input                LOp   : Superoperator
	//                      T1s   : Spin tensor components of mu1
	// 			T2s   : Spin tensor components of mu2
	//			w     : Transisiton frequencies
	//			hs    : Hilbert space size
	//			taus  : 5 effective correlation times
	//			c1s   : J coefficients of mu1
	//			c2s   : J coefficients of mu2
	//			xi1xi2: Interaction constant product 
	//			w0    : Zero quantum transition frequency
	//			w1    : Single quantum transition frequency
	//			w2    : Double quantum transition frequency
	//			l     : Rank of the relaxation mechanisms
	//			level : Relaxation treatment level
	//			autoc : Flag for auto correlation vs cross correlation
	// Output		LOp   : Relaxation superoperator
	// Note			      :	Computed in the basis T1s and T2s


 
MSVCDLL void Rmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w,
                      int hs, double tau, double xi1xi2, double w0,
                            double w1, double w2, int l, int level=4, int autoc=0);

	// Input                LOp   : Superoperator
	//                      T1s   : Spin tensor components of mu1
	// 			T2s   : Spin tensor components of mu2
	//			w     : Transisiton frequencies
	//			hs    : Hilbert space size
	//			tau   : Correlation time
	//			xi1xi2: Interaction constant product 
	//			w0    : Zero quantum transition frequency
	//			w1    : Single quantum transition frequency
	//			w2    : Double quantum transition frequency
	//			l     : Rank of the relaxation mechanisms
	//			level : Relaxation treatment level
	//			autoc : Flag for auto correlation vs cross correlation
	// Output		LOp   : Relaxation superoperator
	// Note			      :	Computed in the basis T1s and T2s
	// Note			      :	Uses spherical top spectral density functions!


MSVCDLL void Rmu1mu2(super_op& LOp, const spin_system& sys, gen_op& Ho, double* w,
         double* xi1s, double n1, double* xi2s, double n2, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int l, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			xi1s  : Interaction constants, mu1
	//			n1    : Number of mu1 spins or spin pairs
	//			xi2s  : Interaction constants, mu2
	//			n1    : Number of mu2 spins or spin pairs
	//			A1    : Spatial tensor, mu1
	//			A2    : Spatial tensor, mu2
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			l     : Rank of the relaxation mechanisms
	// 			type  : Relaxation type to compute
	//				   0 = Auto & Cross Correlation
	//				   + = Auto Correlation Only
	//				   - = Cross Correlation Only
	//			level : Relaxation treatment level


// --------- Spin-Pair with Spin-Pair Interaction Functions -------------

MSVCDLL void Rijkl(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
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



MSVCDLL void Rijkl(super_op& LOp, const spin_system& sys, gen_op& Ho, double* w,
                  matrix& xi1s, matrix& xi2s, spin_T* T1, spin_T* T2,
                                    double tau, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			xi1s  : Interaction constant, mu1
	//			xi2s  : Interaction constant, mu2
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			tau   : System correlation time
	// 			type  : Relaxation type to compute
	//				   0 = Auto & Cross Correlation
	//				   + = Auto Correlation Only
	//				   - = Cross Correlation Only
	//			level : Relaxation treatment level


//	      Start the Relaxation Superoperator Computation
//	             Two Spin-Pair Rank 2 Mechanisms


// -------------------- Spin with Spin Functions -------------------

MSVCDLL void Rij(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
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

/*                    --- ---             --- ---
                      \   \               \   \
               LOp += /   /    R        = /   /   R
  	              --- ---   mu1,mu2   --- ---  i,j
                      mu1 mu2              i   j
*/


MSVCDLL void Rij(super_op& LOp, const spin_system& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, spin_T* T1, spin_T* T2,
                       double tau, int l, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			xi1s  : Interaction constant, mu1
	//			xi2s  : Interaction constant, mu2
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			tau   : System correlation time
	//			l     : Rank of the relaxation mechanism(s)
	// 			type  : Relaxation type to compute
	//				   0 = Auto & Cross Correlation
	//				   + = Auto Correlation Only
	//				   - = Cross Correlation Only
	//			level : Relaxation treatment level
	// Note			      : Uses spherical top spectral densities


/*                    --- ---             --- ---
                      \   \               \   \
               LOp += /   /    R        = /   /   R
  	              --- ---   mu1,mu2   --- ---  i,j
                      mu1 mu2              i   j
*/


// --------------- Spin-Pair with Spin Functions -------------------

MSVCDLL void Rijk(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
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


MSVCDLL void Rkij(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
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

// ______________________________________________________________________
// *RELAXATION SUPEROPERATOR GENERAL FUNCTIONS   DYNAMIC FREQUENCY SHIFT*
// ______________________________________________________________________

// ---------------- Routing to Specific Level Computation ---------------

MSVCDLL void Rmumuds(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
               double w0, double w1, double w2, int level=4, int autoc=0, int het=0);

	// Input                LOp   : Superoperator
	//                      T1s   : Spin tensor components of mu1
	// 			T2s   : Spin tensor components of mu2
	//			w     : Transisiton frequencies
	//			hs    : Hilbert space size
	//			taus  : 5 effective correlation times
	//			c1s   : J coefficients of mu1
	//			c2s   : J coefficients of mu2
	//			xi1xi2: Interaction constant product 
	//			w0    : Zero quantum transition frequency
	//			w1    : Single quantum transition frequency
	//			w2    : Double quantum transition frequency
	//			level : Relaxation treatment level
	//			autoc : Flag for auto correlation vs cross correlation
	// Output		LOp   : Relaxation superoperator
	// Note			      :	Computed in the basis T1s and T2s

//                            LOp += R   (Level)
//                                    1,2  


// --------- Spin-Pair with Spin-Pair Interaction Functions -------------

MSVCDLL void Rijklds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
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

//	             Two Spin-Pair Rank 2 Mechanisms

/*                    --- ---             --- --- --- ---
                      \   \               \   \   \   \
               LOp += /   /    R        = /   /   /   /   R
  	              --- ---   mu1,mu2   --- --- --- ---  ij,kl
                      mu1 mu2              i   j   k   l
*/


MSVCDLL void Rijds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
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

//	          Two Single Spin Rank 2 Mechanisms

/*                    --- ---             --- ---
                      \   \               \   \
               LOp += /   /    R        = /   /   R
  	              --- ---   mu1,mu2   --- ---  i,j
                      mu1 mu2              i   j
*/


// --------------- Spin-Pair with Spin Functions -------------------


MSVCDLL void Rijkds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
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

//	       Mechanism 1 = Spin-Pair; Mechanism 2 = Spin

/*                    --- ---             --- --- ---
                      \   \               \   \   \
               LOp += /   /    R        = /   /   /   R
  	              --- ---   mu1,mu2   --- --- ---  ij,k
                      mu1 mu2              i   j   k
*/


MSVCDLL void Rkijds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
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

//	       Mechanism 1 = Spin; Mechanism 2 = Spin-Pair

/*                    --- ---             --- --- ---
                      \   \               \   \   \
               LOp += /   /    R        = /   /   /   R
  	              --- ---   mu1,mu2   --- --- ---  k,ij
                      mu1 mu2              k   i   j 
*/

#endif 
