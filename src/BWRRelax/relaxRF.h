/* relaxRF.h ************************************-*-c++-*-
**							**
**	               G A M M A			**
**							**
**	NMR Library RF-Field Relaxation Functions	**
**							**
**	Interface definition				**
**							**
**	Copyright (c) 1993				**
**	Scott A. Smith					**
**	University of California, Santa Barbara		**
**	Department of Chemistry				**
**	Santa Barbara CA. 93106 USA			**
**							**
**      $
**							**
*********************************************************/

/*********************************************************
**							**
** 	Description					**
**							**
** Herein are functions dealing with some of the more	**
** common relaxation with an external rf-field problems	**
**							**
*********************************************************/

///Chapter RF-Field Relaxation
///Section Overview
///Body The ...
///Section Available RF-Field Relaxation Functions

#ifndef   Relax_rf_h_			// Is this file already included?
#  define Relax_rf_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <LSLib/SuperOp.h>		// Include superoperators
#include <LSLib/sys_dynamic.h>		// Include dynamic systems
#include <Level1/SpinT.h>		// Include spin tensors

// ____________________________________________________________________________
// ******* RF-FIELD RELAXATION SUPEROPERATORS VIA ELEMENT CALCULATIONS ********
// ____________________________________________________________________________

// -------------------- Level 4 via Element Calculation -----------------------


	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st mu
	// 			T2s   : Spin tensor components, 2nd mu
	//			J12   : Matrices of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//				evaluated at m based frequencies
	//			rank  : Rank of the 2 interactions
	// Return		void  : Relaxation superoperator, interactions
	//				1 & 2 added to input superoperator LOp
	//
	//				LOp    = LOp   + LOp
	//				   out	    in      12
	//
	// Note				T1s, T2s, & LOp assumed in proper bases


	//                  <a,a'| LOp |b,b'> += <a,a'| R   |b,b'>
	//                                               1,2
   
	// where a, a', b, b' are basis funciton indices.


MSVCDLL void Rrf_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12);
MSVCDLL complex Rrf_4(int hs, gen_op* T1s, gen_op* T2s, matrix* J12,
		 	         int rank, int a, int b, int aa, int bb);

	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrices of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//				evaluated at m based frequencies
	//				times the two interaction constants
	//			rank  : Rank of the 2 interactions
	//			a, b  : 1st transisiton indices
	//			aa, bb: 2nd transisiton indices
	// Output		Rel   : Relaxation superoperator element
	//				for interactions 1 & 2 for component m
	//
	//				    <a,aa|R  |b,bb>
	//					   12
	//
	// Note				T1s, T2s assumed in proper bases

          
/*                    --- [ ---
                      \   | \                 m       m  
   <a,a'|R   |b,b'> = /   | /   delta     <a|T |g><b|T |g> J  (w  - mW  )
          1,2         --- | ---      a',b'    1       2     12  gb    rf
                       m  [  g
  
                                 m        m
                           - <a|T |b><a'|T |b'> J  (w    - mW  )
                                 1        2      12  b'a'    rf

                                           m        m           
                                    - <b'|T |a'><b|T |a> J  (w  - mW  )
                                           1        2     12  ab    rf

                              ---                                            ]
                              \               m        m                     |   
                           +  /   delta   <g|T |a'><g|T |b'> J  (w   - mW  ) |
                              ---      a,b    1        2      12  b'g    rf  |
                               g                                             ]

   where a, a', b, b' and g are basis function indices.  The index m sums over
   components of angular momentum and the spectral densities are scaled by the
   two interaction constants.						     */


// -------------------- Level 3 via Element Calculation -----------------------
//                      (Applies Secular Approximation)

MSVCDLL void Rrf_3(super_op& LOp, double* w, int rank,
                   gen_op* T1s, gen_op* T2s, matrix* J12, double cutoff=1.e-2);

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector system enery levels of Heff
	//				(rad/sec) in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st mu
	// 			T2s   : Spin tensor components, 2nd mu
	//			J12   : Matrices of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//				evaluated at m based frequencies
        //                      cutoff: Secular approximation cutoff value (Hz)
	// Return		void  : Relaxation superoperator for interactions
	//				1 & 2 added to input superoperator LOp
	//
	//				LOp    = LOp   + LOp
	//				   out	    in      12
	//
	// Note				T1s, T2s, & LOp assumed in proper basis (Heff)


//                  <a,a'| LOp |b,b'> += <a,a'| R   |b,b'>
//                                               1,2
   
// where a, a', b, b' are basis funciton indices.


 
// -------------------- Level 2 via Element Calculation -----------------------
//                        (No Degenerate Transitions)

MSVCDLL void Rrf_2(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12);


	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrices of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//				evaluated at m based frequencies
	// Return		void  : Level 2 relaxation superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R2
	//				          out	   in     12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases


//                  <a,a'| LOp |b,b'> += <a,a'| R2   |b,b'>
//		                                  1,2

// where a, a', b, b' are basis function indices.


MSVCDLL double Rrf_2(int hs, gen_op* T1s, gen_op* T2s, matrix* J12,
		 	         int rank, int a, int b, int aa, int bb);

	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrices of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//				evaluated at m based frequencies
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
   <a,a'|R2   |a,a'> = /   | /   <a|T |g><a|T |g>*J(w  - m*W  )
  	   1,2	       --- | ---     1       2       ga     rf
  		        m  [  g
  
                                                   m        m
                                       - 2.0 * <a|T |a><a'|T |a'> * J(w  - m*W  )
  	     	                                   1        2          aa     rf
  
                                              ---                                   ]
                                              \       m        m                    |
                                           +  /   <g|T |a'><g|T |a'>*J(w   - m*W  ) |
  	     	                              ---     1        2        ga'     rf  |
  								                    ]
  
  			Level 2 Non-Zero Off Diagonal Elements
  
                     ---
                     \              m       m 
   <a,a|R2   |b,b> = /   -2.0 * <a|T |b><a|T |b> * J  (w   - m*W  )
  	  1,2	     ---            1       2       12  ba      rf
                      m
*/

// where m sums over angular momentum components, g over basis functions,
// and J values include xi1 and xi2

// ---------------- Level 1 via Element Calculation ---------------------


// ---------------- Level 0 via Element Calculation ---------------------
//			  (Extreme Narrowing)

MSVCDLL void Rrf_0(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, const complex& J12);

	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Reduced spectral density at 0 Hz 
	//				(extreme narrowing) for 2 interactions
	//				times the two interaction constants
	// Return		void  : Level 0 relaxation superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R0
	//				          out	   in     12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases


//                  <a,a'| LOp |b,b'> += <a,a'| R0   |b,b'>
//		                                  1,2

// where a, a', b, b' are basis function indices.  This assumes extreme narrowing.

MSVCDLL double Rrf_0(int hs, gen_op* T1s, gen_op* T2s,
		 	         int rank, int a, int b, int aa, int bb);


	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			rank  : Rank of the two interactions
	//			a, b  : 1st transition indices
	//			aa, bb: 2nd transition indices
	// Output		Rel   : Level 0 relaxation superop element
	//
	//				    <a,aa|R0  |b,bb>
	//					    12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases
	// Note				LOp element return is unitless!


/*                     --- [ ---
                       \   | \                 m       m 
   <a,a'|R0   |b,b'> = /   | /   delta     <a|T |g><b|T |g>
  	   1,2	       --- | ---      a',b'    1       2
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

// ----------------- Level 4 by Double Commutators ----------------------

MSVCDLL void Rrf_4s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12);

	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components 1st interaction
	// 			T2s   : Spin tensor components 2nd interaction
	//			J12   : Matrices of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//				evaluated at m based frequencies
	// Void			LOp   : Relaxation superoperator elements
	//				for interactions 1 & 2 added to input
	//				superoperator elements
	//
	//				LOp    = LOp   + LOp
	//				   out	    in      12
	//
	// Note				All operators contained in Ho,
	//				T1s, and T2s are assumed to be in
	//				the static Hamiltonian eigenbasis

/*                    --- ---
                      \   \       m    m     -m
  	   R  (1,2) = /   /   (-1)  [ T , [ T   ,  ] ] J  (w - mW  )
  	    rf	      --- ---          1     2,p        12  p    rf
  		       p   m
*/


// ----------------- Level 3 by Double Commutators ----------------------
//                  (Applies Secular Approximation)

MSVCDLL void Rrf_3s(super_op& LOp, double* w, int rank, gen_op* T1s, gen_op* T2s, matrix* J12);

	// Input		LOp   : Current relaxation superoperator
	//			w     : Double vector of enery levels
	//				(rad/sec) for static Hamiltonian
	//				in the LAB frame (range in 10**8 Hz)
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components of
	//				first interaction
	// 			T2s   : Spin tensor components of
	//				second interaction
	//			J12   : Matrices of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//				evaluated at m based frequencies
	// Void			LOp   : Relaxation superoperator elements
	//				for interactions 1 & 2 added to input
	//				superoperator elements
	//
	//				LOp    = LOp   + LOp
	//				   out	    in      12
	//
	// Note				All operators contained in Ho,
	//				T1s, and T2s are assumed to be in
	//				the static Hamiltonian eigenbasis
	// Note				Would be nice if we could use double
	//				rather than complex here!! 


// ------------------ Level 2 via Double Commutators --------------------
//                      (No Degenerate Transitions)

MSVCDLL void Rrf_2s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12);

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector of system energies
	//				in rad/sec in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	// Return		void  : Level 2 relaxation superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R2
	//				          out	   in     12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases

/*                   --- ---
                     \   \       m    m       -m
             R2    = /   /   (-1)  [ T   , [ T   ,  ] ] J  (w - m*W  )
  	       1,2   --- ---          1,p     2,p        12  p     rf
  		      p   m
*/

// where p sums over transitions, m sums over angular momentum components,
// and J values include xi1 and xi2.  The Level 2 treatment assumes that there
// are no degenerate transitions in the system.


// ------------------ Level 1 via Double Commutators --------------------
//              (Spin Tensor Component Constant Oscillation)

//void Rrf_1s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12)

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector of system energies
	//				in rad/sec in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	// Return		void  : Level 1 relaxation superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R1
	//				          out	   in     12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases

/*                         ---
                           \       m    m     -m
                   R1    = /   (-1)  [ T , [ T  ,  ] ] J  (mw )
  	             1,2   ---          1     2         12   0
  		            m
*/

// m sums over angular momentum components, and J values include xi1 and xi2.
// The Level 1 treatment assumes that there are no degenerate transitions in the system.


// ------------------ Level 0 via Double Commutators --------------------
//                         (Extreme Narrowing)

//void Rrf_0s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, double J12)

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector of system energies
	//				in rad/sec in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Reduced spectral density function value
	//				at 0 Hz times the 2 interaction constants
	// Return		void  : Level 0 relaxation superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R0
	//				          out	   in     12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases

/*                               ---
                                 \       m    m     -m
                   R0 ,  = J (0) /   (-1)  [ T , [ T  ,  ] ]
  	             1,2    12   ---          1     2
  		                  m
*/

// m sums over angular momentum components, and J value include xi1 and xi2.
// The Level 1 treatment assumes extreme narrowing.


// ______________________________________________________________________
// ************* RELAXATION SUPEROPERATOR GENERAL FUNCTIONS *************
// ______________________________________________________________________

// ---------------- Routing to Specific Level Computation ---------------

MSVCDLL void Rrfmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, matrix* J12,
                 double* J, double* w, int rank=2, int level=4, int autoc=0, int het=0);

	// Input                LOp   : Superoperator
	//                      T1s   : Spin tensor components of mu1
	// 			T2s   : Spin tensor components of mu2
	//			J12   : Array of reduced spectral density
	//				matricies scaled by xi1*xi2
	//			J     : Vector containing J0, J1, & J2
	//				scaled by xi1*xi2
	//			rank  : Rank of the two interactions
	//			level : Relaxation treatment level
	//			autoc : Flag auto- vs cross- correlation
	// Output		void  : Relaxation superoperator LOp altered
	// Note			      :	Computed in the basis T1s and T2s

//                            LOp += R   (Level)
//                                    1,2  


// --------- Spin-Pair with Spin-Pair Interaction Functions -------------

MSVCDLL void Rrfijkl(super_op &LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
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


// -------------------- Spin with Spin Functions -------------------


MSVCDLL void Rrfij(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
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

MSVCDLL void Rrfijk(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			w     : Transition frequencies
	//			Wrflab: Lab frame rf-field frequency
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



MSVCDLL void Rrfkij(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			w     : Transition frequencies
	//			Wrflab: Lab frame rf-field frequency
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


// ______________________________________________________________________
// *RELAXATION SUPEROPERATOR GENERAL FUNCTIONS   DYNAMIC FREQUENCY SHIFT*
// ______________________________________________________________________


// --------- Spin-Pair with Spin-Pair Interaction Functions -------------

MSVCDLL void Rrfijklds(super_op &LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
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


// -------------------- Spin with Spin Functions -------------------

MSVCDLL void Rrfijds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type=0, int level=4);

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
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

MSVCDLL void Rrfijkds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level);

        // Input                LOp   : Superoperator
        //                      sys   : Dynamic spin system
        //                      Heff  : Effective Hamiltonian
        //                      w     : Transition frequencies
        //                      Wrflab: Lab frame rf-field frequency
        //                      xi1s  : Interaction constant, mu1
        //                      xi2s  : Interaction constant, mu2
        //                      A1    : Spatial tensor, mu1
        //                      A2    : Spatial tensor, mu2
        //                      T1    : Spin tensor, mu1
        //                      T2    : Spin tensor, mu2
        //                      taus  : Sytem 5 effective correlation times
        //                      chi   : Spin system chi value
        //                      type  : Relaxation type to compute
        //                                 0 = Auto & Cross Correlation
        //                                 + = Auto Correlation Only
        //                                 - = Cross Correlation Only
        //                      level : Relaxation treatment level

//             Mechanism 1 = Spin-Pair; Mechanism 2 = Spin

/*                    --- ---             --- --- ---
                      \   \               \   \   \
               LOp += /   /    R        = /   /   /   R
                      --- ---   mu1,mu2   --- --- ---  ij,k
                      mu1 mu2              i   j   k
*/

 
MSVCDLL void Rrfkijds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level);
 
        // Input                LOp   : Superoperator
        //                      sys   : Dynamic spin system
        //                      Heff  : Effective Hamiltonian
        //                      w     : Transition frequencies
        //                      Wrflab: Lab frame rf-field frequency
        //                      xi1s  : Interaction constant, mu1
        //                      xi2s  : Interaction constant, mu2
        //                      A1    : Spatial tensor, mu1
        //                      A2    : Spatial tensor, mu2
        //                      T1    : Spin tensor, mu1
        //                      T2    : Spin tensor, mu2
        //                      taus  : Sytem 5 effective correlation times
        //                      chi   : Spin system chi value
        //                      type  : Relaxation type to compute
        //                                 0 = Auto & Cross Correlation
        //                                 + = Auto Correlation Only
        //                                 - = Cross Correlation Only
        //                      level : Relaxation treatment level
 
//             Mechanism 1 = Spin; Mechanism 2 = Spin-Pair
 
/*                    --- ---             --- --- ---
                      \   \               \   \   \
               LOp += /   /    R        = /   /   /   R
                      --- ---   mu1,mu2   --- --- ---  k,ij
                      mu1 mu2              k   i   j
*/
 

// ______________________________________________________________________
// ***************** STEADY STATE DENSITY MATRIX FUNCTIONS **************
// ______________________________________________________________________

// ----------------- Level 4 by Double Commutators ----------------------

MSVCDLL gen_op sigma_ss(spin_system& sys, super_op& L, super_op& R);

	// Input		sys   : Spin system
	// 			L     : Full Liouvillian
	//			R     : Relaxation Superoperator
        // Return		sss   : Steady state density matrix
	// Note			      : Both R and L must be in the same
	//				units (rad/sec).  Be careful because
	//		 		usually Ho is generated in Hz!
	// Note			      : Because L is singular, the steady
	//				state density matrix must be
	//				determined via a reduced Liouville
	//				space.  In turn, this function is
	//				forced to use matrix routines rather
	//				than the (code simple) superoperator
	//				functions!
	// Note			      : Because this algorithm assumes the
	//				trace of all density matrices is 1
	//				we must first set the equilibrium matrix
	//				to satisfy this condition them reset the
	//				final answer back to trace = 0 as is used
	//				the rest of the GAMMA routines

//	L|s  > = R|s  >    -------> [L - L|S><E|]|s  > = R|s  > - |L1>
//         ss       eq		                   ss       eq


MSVCDLL gen_op sigma_ss_it(spin_system& sys, super_op& L, super_op& Heff, super_op& R);

	// Input		sys   : Spin system
	// 			L     : Full Liouvillian
	//			Heff  : Hamiltonian Superoperator
	//			R     : Relaxation Superoperator
        // Return		sss   : Steady state density matrix
	// Note			      : Because L is singular, the steady
	//				state density matrix must be
	//				determined via a reduced Liouville
	//				space.  In turn, this function is
	//				forced to use matrix routines rather
	//				than the (code simple) superoperator
	//				functions!


#endif /* __RELAX_RF_H__ */
