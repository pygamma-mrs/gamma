/* relaxRF.cc ********************************************
**							**
** 	            G A M M A				**
**							**
**	Relaxation with an External RF-Field	 	**
**						 	**
**	Implementation   				**
**						 	**
**	Copyright (c) 1991, 1992, 1993		 	**
**	Scott A. Smith				 	**
**							**
**	University of California, Santa Barbara		**
**	Department of Chemistry				**
**	Santa Barbara CA. 93106 USA			**
**						 	**
**      $Header: $
**						 	**
** Modifications:					**
**						 	**
*********************************************************/

/*********************************************************
**						 	**
** 	Description				 	**
**						 	**
**	The following functions are used to generate	**
**	relaxation superoperators and related species	**
**	when there is an externally applied rf-field.	**
**						 	**
*********************************************************/

#ifndef _relax_rf_cc_		// Is this file already included?
#define _relax_rf_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation          // This is the implementation
#endif

#include <BWRRelax/relaxRF.h>		// Include the header info
#include <BWRRelax/relaxNMR.h>		// Include relaxation controls
#include <LSLib/DensOp.h>		// Include density operators
#include <LSLib/LSAux.h>		// Include inversion routines
#include <HSLib/HSauxil.h>		// Inlcude sigma_eq functio
#include <stdlib.h>

// ____________________________________________________________________________ 
// ******* RF-FIELD RELAXATION SUPEROPERATORS VIA ELEMENT CALCULATIONS ********
// ____________________________________________________________________________ 

// -------------------- Level 4 via Element Calculation -----------------------

void Rrf_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12)

	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st mu
	// 			T2s   : Spin tensor components, 2nd mu
	//			J12   : Matrices of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//				evaluated at m based frequencies
	//			rank  : Rank of the 2 interactions
	// Return		void  : Level 4 relax. superop for interactions
	//				1 & 2 added to input superoperator LOp
	//
	//				LOp    = LOp   + R4
	//				   out	    in     12
	//
	// Note				T1s, T2s, & LOp assumed in proper bases


	//                  <a,a'| LOp |b,b'> += <a,a'|R4   |b,b'>
	//                                               1,2
   
	// where a, a', b, b' are basis function indices.

   {
   int hs = T1s[0].dim();			// Get Hilbert space size
   int aaa=0, bbb;
   complex Rel;
   for(int a=0; a<hs; a++)			// Sum over transition a-aa
     for(int aa=0; aa<hs; aa++)
       {
       bbb = 0;
       for(int b=0; b<hs; b++)			// Sum over transition b-bb
         for(int bb=0; bb<hs; bb++)
           {
           Rel = LOp.get(aaa,bbb);		// Get the current LOp element
           Rel += Rrf_4(hs,T1s,T2s,J12,		// Add the R component
			       rank,a,b,aa,bb);
           LOp.put(aaa, bbb, Rel);		// Put in the new element
           bbb++;
           }
       aaa++;
       }
   return;
   }


complex Rrf_4(int hs, gen_op* T1s, gen_op* T2s, matrix* J12,
		 	         int rank, int a, int b, int aa, int bb)
//return Rel(0,0);

	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrices of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//				evaluated at m based frequencies
	//				times the two interaction constants
	//			rank  : Rank of the 2 interactions
	//			a, b  : 1st transition indices
	//			aa, bb: 2nd transition indices
	// Output		Rel   : Level 4 relaxation superop element
	//				for interactions 1 & 2 for component m
	//
	//				    <a,aa|R4  |b,bb>
	//					    12
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
   
*/
// where a, a', b, b' and g are basis function indices.  The index m sums over
// components of angular momentum and the spectral densities are scaled by the
// two interaction constants.


   {
   complex Rel(0,0);
   int k=0, g=0;
   for(int m=-rank; m<=rank; m++)		// Sum over tensor components
     {
     Rel -= J12[k].get(bb,aa)*			// Add in terms RII
             T1s[k].get(a,b)*
              conj(T2s[k].get(aa,bb));
     Rel -= J12[k].get(a,b)* 			// Add in terms RIII
             T1s[k].get(bb,aa)*
              conj(T2s[k].get(b,a));
     for(g=0; g<hs; g++)
       {
       if(aa == bb)				// Add terms RI over gamma sum
         Rel += J12[k].get(g,b)*
		T1s[k].get(a,g)*
                conj(T2s[k].get(b,g));
       if(a == b) 				// Add terms RIV over gamma sum
         Rel += J12[k].get(bb,g)*
                T1s[k].get(g,aa)*
                conj(T2s[k].get(g,bb));
       }
     k++;
     }
   return Rel;
   }


 
// -------------------- Level 3 via Element Calculation -----------------------
//                      (Applies Secular Approximation)
 
void Rrf_3(super_op& LOp, double* w, int rank, gen_op* T1s,
                                         gen_op* T2s, matrix* J12, double cutoff)

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
	// Return		void  : Level 3 relax. superop for interactions
	//				1 & 2 added to input superoperator LOp
	//
	//				LOp    = LOp   + R3
	//				   out	    in     12
	//
	// Note				T1s, T2s, & LOp assumed in
	//				proper bases (eigenbasis of Heff)


// <a,a'| LOp |b,b'> += <a,a'| R3   |b,b'> == del           * <a,a'|R4   |b,b'>
//                               1,2             w    , w             1,2
//                                                a,a'   b,b'
   
// where a, a', b, b' are basis function indices.


   {
   int hs = T1s[0].dim();			// Get Hilbert space size
   double delwaa, delwbb;
   int aaa=0, bbb=0;
   complex Rel;
   for(int a=0; a<hs; a++)			// Sum over transition a-aa
     for(int aa=0; aa<hs; aa++)
       {
       delwaa = w[a] - w[aa];			// Transition a-aa frequency
       bbb = 0;
       for(int b=0; b<hs; b++)			// Sum over transition b-bb
         for(int bb=0; bb<hs; bb++)
           {
           delwbb = w[b] - w[bb];		// Transition b-bb frequency
	   if(fabs(delwaa-delwbb) < cutoff)	// Apply secular approximation
             {
             Rel = LOp.get(aaa,bbb);
             Rel += Rrf_4(hs,T1s,T2s,J12,
			       rank,a,b,aa,bb);
             LOp.put(aaa, bbb, Rel);		// Add to the relax. mx. element
             }
           bbb++;
           }
       aaa++;
       }
   return;
   }

  
// -------------------- Level 2 via Element Calculation -----------------------
//                        (No Degenerate Transitions)

void Rrf_2(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12)


	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrices of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//				evaluated at m based frequencies
	// Return		void  : Level 2 relax. superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R2
	//				          out	   in     12
	//
	// Note				T1s, T2s, & LOp assumed in proper bases


//                  <a,a'| LOp |b,b'> += <a,a'| R2   |b,b'>
//		                                  1,2

// where a, a', b, b' are basis function indices.

   {
   int hs = T1s[0].dim();			// Get Hilbert space size
   int aaa=0, bbb=0;
   complex Rel;
   for (int a=0; a<hs; a++)			// Sum over transitions a-aa
     for (int aa=0; aa<hs; aa++)
       {
       bbb = 0;
       for (int b=0; b<hs; b++)			// Sum over transitions b-bb
         for (int bb=0; bb<hs; bb++)
           {
           Rel = LOp.get(aaa,bbb);
           Rel += Rrf_2(hs,T1s,T2s,
			   J12,rank,a,b,aa,bb);
           LOp.put(aaa, bbb, Rel);
           bbb++;
           }
       aaa++;
       }
   return;
   }


//sosi switch this to complex
double Rrf_2(int hs, gen_op* T1s, gen_op* T2s, matrix* J12,
		 	         int rank, int a, int b, int aa, int bb)

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
	// Note				T1s, T2s, & LOp assumed in proper bases

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

   where m sums over angular momentum components, g over basis functions,
   and J values include xi1 and xi2                                         */

   {
   complex Rel = 0;
   int k=0, g=0, m=0;
   matrix J;
   complex J0;
   if((a==b) && (aa==bb))			// Algorithm for diagonals
     {
     for(m=-rank, k=0; m<=rank; m++, k++)	// Sum over tensor components
       {
       J = J12[k];
       J0 = 2.0*J.get(a,a);
       Rel -= J0*T1s[k].get(a,a)*T2s[k].get(aa,aa);
       for(g=0; g<hs; g++)
         {
         Rel += 2.0*J.get(g,a)*
		    T1s[k].get(a,g)*T2s[k].get(a,g);
         Rel += 2.0*J.get(aa,g)*
                  T1s[k].get(g,aa)*T2s[k].get(g,aa);
         }
       }
     }
   else if((a==aa) && (b==bb) && (a!=b))	// Algorithm non-zero off-diags
     {
     for(m=-rank, k=0; m<=rank; m++, k++)	// Sum over tensor components
       {
       J = J12[k];
       Rel -= 2.0*J.get(b,a)*
                T1s[k].get(a,b)*T2s[k].get(a,b);
       }
     }
   return Re(Rel);
   }


// ---------------- Level 1 via Element Calculation ---------------------


// ---------------- Level 0 via Element Calculation ---------------------
//                        (Extreme Narrowing)

void Rrf_0(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, const complex& J12)

	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Reduced spectral density at 0 Hz 
	//				(extreme narrowing) for 2 interactions
	//				times the two interaction constants
	// Return		void  : Level 0 relax superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R0
	//				          out	   in     12
	//
	// Note				T1s, T2s, & LOp assumed in proper bases


//                  <a,a'| LOp |b,b'> += <a,a'| R0   |b,b'>
//		                                  1,2

// where a, a', b, b' are basis function indices.  Assumes extreme narrowing.

   { R_0(LOp, rank, T1s, T2s, J12); } 		// Same as no rf-field equation

// sosi switch this to complex
double Rrf_0(int hs, gen_op* T1s, gen_op* T2s,
		 	         int rank, int a, int b, int aa, int bb)


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
	// Note				T1s, T2s, & LOp assumed in proper bases
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

   { return R_0(hs,T1s,T2s,rank,a,b,aa,bb); } 	// Same as no rf-field equation

// ______________________________________________________________________
// ************ RELAXATION SUPEROPERATORS VIA DOUBLE COMMUTATORS ********
// ______________________________________________________________________

// ----------------- Level 4 by Double Commutators ----------------------

void Rrf_4s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12)

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
	// Note				T1s, T2s, & LOp assumed in proper bases

/*                    --- ---
                      \   \       m    m     -m
  	   R  (1,2) = /   /   (-1)  [ T , [ T   ,  ] ] J  (w - mW  )
  	    rf	      --- ---          1     2,p        12  p    rf
  		       p   m                                                 */

   {
   int hs = T1s[0].dim();			// Get Hilbert space size
   matrix mx(hs, hs, 0.0);			// Null matrix
   matrix mxtmp;
   basis bs = T1s[0].get_basis();		// Current Basis
   gen_op nullOp;
   int k=0;
   double Jwp = 0;
   complex Taaa;
   gen_op *T2ps;				// Storage for 5 T    cmpnts.
   T2ps = new gen_op[2*rank+1];			//                2,p

   for(int a=0; a<hs; a++)			// Sum over transition a-aa (p)
     for(int aa=0; aa<hs; aa++)			// and generate the 5 T   cmps.
       {					//		       2,p
       k=0;
       for(int m=-rank; m<=rank; m++)		// Sum all tensor components
         {
         T2ps[k] = nullOp;			// Spin tensor operator component
         Taaa = T2s[k].get(a,aa);		// T2ps[k] is either null or has
         if(Re(Taaa) || Im(Taaa))		// 1 non-zero element!  To avoid
           {					// some of the needless matrix
           mx.put(Taaa, a, aa);			// copying mx stores then rezeroes
           T2ps[k] = gen_op(mx,bs);		// the element
           mx.put(0.0, a, aa);
           }
//         Jwp = Re(J12[k].get(a,aa));		// This line doesn't work unfortunately
         mxtmp = J12[k];
         Jwp = Re(mxtmp.get(a,aa));
// Should this be restricted to sum only over the current m of T1s?
//         R_CC_0(T1s, T2ps, LOp, rank, Jwp);	// This will sum over m's again
         k++;
         }
       R_CC_0(T1s, T2ps, LOp, rank, Jwp);	// This will sum over m's again
       }
   return;
   }

// ----------------- Level 3 by Double Commutators ----------------------
//                  (Applies Secular Approximation)

void Rrf_3s(super_op& LOp, double* w, int rank, gen_op* T1s, gen_op* T2s, matrix* J12)

	// Input		LOp   : Current relaxation superoperator
	//			w     : Double vector of enery levels
	//				(rad/sec) for static Hamiltonian
	//				in the LAB frame (range in 10**8 Hz)
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
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
	// Note				T1s, T2s, & LOp assumed in proper bases

   {
// sosi UNDER CONSTRUCTION: THIS  DOESN'T WORK YET
   int hs = T1s[0].dim();			// Get Hilbert space size
   matrix mx1(hs, hs, 0.0);			// Two zero matrices of the
   matrix mx2(hs, hs, 0.0);			// Hilbert space size
   basis bs = T1s[0].get_basis();		// Current Basis
   gen_op nullOp;				// A NULL operator
   double delwaa, delwbb;
   int k=0, kk=0;
//   double Jwp;
   complex Ta_aa, Tb_bb;
   gen_op *T1ps;
   T1ps = new gen_op[2*rank+1];
   gen_op *T2ps;
   T2ps = new gen_op[2*rank+1];
double Js[5];
   for(int a=0; a<hs; a++)			// Sum over transitions a-aa: p
     for(int aa=0; aa<hs; aa++)
       {
       delwaa = w[a] - w[aa];			// Transition a-aa: p frequency
//       Jwp = Re(J12.get(a,aa));
       k=0;
       for(int m=-rank; m<=rank; m++)		// Sum over all tensor components
         {					// and generate the T1ps eigenoperators
         T1ps[k] = nullOp;
         Ta_aa = T1s[k].get(a,aa);
         if(Re(Ta_aa) || Im(Ta_aa))
           {
           mx1.put(Ta_aa, a, aa);
           T1ps[k] = gen_op(mx1,bs);
           mx1.put(0.0, a, aa);
           }
Js[m+2] = Re(J12[m+2].get(a,aa));
         k++;
         }
       for (int b=0; b<hs; b++)			// Sum over transitions b-bb: p'
         for (int bb=0; bb<hs; bb++)
           {
           delwbb = w[b] - w[bb];		// Transition b-bb: p' frequency
	   if(fabs(delwaa+delwbb) < 1.e-9)	// Secular approximation: p = -p'
             {
             kk=0;
             for (int mm=-rank; mm<=rank; mm++)	// Sum over all tensor components
               {				// and generate the T2ps eigenoperators
               T2ps[kk] = nullOp;
               Tb_bb = T2s[kk].get(b,bb);
               if(Re(Tb_bb) || Im(Tb_bb))
                 {
                 mx2.put(Tb_bb, b, bb);
                 T2ps[kk] = gen_op(mx2,bs);
                 mx2.put(0.0, b, bb);
                 }
               kk++;
               }
// R_CC_0(T1ps, T2ps, LOp, rank, Js);
//               R_CC_0(T1ps, T2ps, LOp, rank, Jwp);
// the above step should be reviewed.
             }
           }
       }
   return;
   if(LOp.dim() && J12==NULL) hs=0;		// Compiler likes these used 
   }

// ------------------ Level 2 via Double Commutators --------------------
//                      (No Degenerate Transitions)

void Rrf_2s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12)

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector of system energies
	//				in rad/sec in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	// Return		void  : Level 2 relax superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R2
	//				          out	   in     12
	//
	// Note				T1s, T2s, & LOp assumed in proper bases

/*                   --- ---
                     \   \       m    m       -m
             R2    = /   /   (-1)  [ T   , [ T   ,  ] ] J  (w - m*W  )
  	       1,2   --- ---          1,p     2,p        12  p     rf
  		      p   m                                                  */

// where p sums over transitions, m sums over angular momentum components,
// and J values include xi1 and xi2.  The Level 2 treatment assumes that there
// are no degenerate transitions in the system.

// sosi UNDER CONSTRUCTION: THIS  DOESN'T WORK YET
   {
   int hs = T1s[0].dim();			// Get Hilbert space size
   matrix mx1(hs, hs, 0.0);			// Two Zeroed Hilbert space mxs
   matrix mx2(hs, hs, 0.0);
   basis bs = T1s[0].get_basis();		// Current working basis
   gen_op nullOp;				// A NULL operator
//   double Jwp;				// Spectral density, transition p
   complex Ta_aa;
   gen_op *T1ps, *T2ps;				// Components of the two spin tensors
   T1ps = new gen_op[2*rank+1]; 		// (2*l+1 of these for each p)
   T2ps = new gen_op[2*rank+1];			// (2*l+1 of these for each p)
   int m=0, k=0;				// Angular momentum indicies
   for(int a=0; a<hs; a++)			// Sum over transitions a-aa: p
     for(int aa=0; aa<hs; aa++)
       {
       for(m=-rank, k=0; m<=rank; m++, k++)	// Perform pre-summation over m
         {					// & set up spin tensor operator
         T1ps[k] = nullOp;			// components this transition
         Ta_aa = T1s[k].get(a,aa);		// 
         if(Re(Ta_aa) || Im(Ta_aa))		//                       m
           { 					//           T1ps[k] == T
           mx1.put(Ta_aa, a, aa); 		//                       1,p
           T1ps[k] = gen_op(mx1,bs);
           mx1.put(0.0, a, aa);
           }
         T2ps[k] = nullOp;
         Ta_aa = T2s[k].get(a,aa);		//                       m
         if(Re(Ta_aa) || Im(Ta_aa))		//           T2ps[k] == T
           {					//                       2,p
           mx2.put(Ta_aa, aa, a);		//
           T2ps[k] = gen_op(mx2,bs);		// These are sparse or null
           mx2.put(0.0, aa, a);			// so insure null components are
           }					// recognized (& not used).
         }
// Need a vector of Jwp[m] here instead
//Jwp = Re(J12.get(a,aa));			// Get spectral density for this transition
//R_CC_0_trans(T1ps, T2ps, LOp, rank, Jwp);// Now, perform sum over m, and add to LOp
// sosi - the above is still goofy 
       }
   return;
   if(LOp.dim() && J12==NULL) hs=0;		// Compiler likes these used 
   }

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
	// Note				T1s, T2s, & LOp assumed in proper bases

/*                         ---
                           \       m    m     -m
                   R1    = /   (-1)  [ T , [ T  ,  ] ] J  (mw )
  	             1,2   ---          1     2         12   0
  		            m                                                */

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
	// Note				T1s, T2s, & LOp assumed in proper bases

/*                               ---
                                 \       m    m     -m
                   R0 ,  = J (0) /   (-1)  [ T , [ T  ,  ] ]
  	             1,2    12   ---          1     2
  		                  m                                          */

// m sums over angular momentum components, and J value include xi1 and xi2.
// The Level 1 treatment assumes extreme narrowing.


// ______________________________________________________________________
// ************* RELAXATION SUPEROPERATOR GENERAL FUNCTIONS *************
// ______________________________________________________________________

// ---------------- Routing to Specific Level Computation ---------------

   void Rrfmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, matrix* J12,
                  double* J, double* w, int rank, int level, int autoc, int het)

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
   {
het=0;
   double cutoff = 1.e-6;
   switch(level)
     {
     case 4:					// Level 4 mu1-mu2: element by element
//       Rrf_4(LOp, rank, T1s, T2s, J12);
       Rrf_3(LOp, w, rank, T1s, T2s, J12, 1.e6);
       break;
     case -4:					// Level 4 mu1-mu2: double commutator
       Rrf_4s(LOp, rank, T1s, T2s, J12);
       break;
     case 3:					// Level 3 mu1-mu2: element by element
       Rrf_3(LOp, w, rank, T1s, T2s, J12);
       break;
     case -3:					// Level 3 mu1-mu2: double commutator
//       Rrf_3s(LOp, w, rank, T1s, T2s, J12);
       break;
     case 2:					// Level 2 mu1-mu2: element by element
//       Rrf_2(LOp, rank, T1s, T2s, J12);
       break;
     case -2:					// Level 2 mu1-mu2: double commutator
//       Rrf_2s(LOp, rank, T1s, T2s, J12);
       break;
     case 1: 					// Level 1 mu1-mu2: double commutator
       if(autoc)
         R_AC_1(T1s,LOp,rank,J[0],J[1],J[2]);
       else
         R_CC_1(T1s,T2s,LOp,rank,J[0],J[1],J[2]);
       break;
     case 0:					// Level 0 mu1-mu2: element by element
         Rrf_0(LOp,rank,T1s,T2s,complex(J[0]));
       break;
     default:					// Level 0 mu1-mu2: double commutator
       if(fabs(J[0]) > cutoff)
       { if(autoc)
           R_AC_0(T1s, LOp, rank, J[0]);
         else
           R_CC_0(T1s,T2s,LOp,rank,J[0]);
       }
       break;
     }
   return;
   }


// --------- Spin-Pair with Spin-Pair Interaction Functions -------------

   void Rrfijkl(super_op &LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

   {
// sosi added 3/11/96
int het = sys.heteronuclear();                  // Flag if system is heteronuclear
int rank=2;
// sosi: this is still dipole specific!!
matrix theta = sys.thetas();			// Get dipole theta values (radians) 
matrix phi = sys.phis(); 			// Get dipole phi values (radians)
double alphaij, betaij;			// Two Euler angles for ij dipole
double alphakl, betakl;			// Two Euler angles for kl dipole
   double xi1, xi2, xi1xi2=0;			// Specific interaction constants
   coord EA1, EA2;				// Specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
//   int ls = hs*hs;				// Total system Liouville space
   int m;					// z-component of spin ang. momentum
   int ij = 0;                                  // Sum over all dipolar pairs
   int kl = 0;                                  // and build relaxation matrix

//
//		    Prepare the Interaction Constants
//		     and Spectral Density Components

  double mWrf=0;
  double J[3];
  matrix* J12 = NULL;                          // Matrices of spectral densities
  if(abs(level) > 1)                           // Needed for higher level computations
    {
    J12 = new matrix[5];
    Heff.eigvals(w);                           // Frequencies w in field rotating frame
    }

//	      Start the Relaxation Superoperator Computation

   for(int i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)
       {
       xi1 = Re(xi1s.get(i,j));			// 1st interaction constant, i&j
alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
/*      Jcoeffs(c1s, alphaij, betaij, 0.0, chi);	// Set the 5 coefficients for ij*/
EA1.xyz(alphaij, betaij, 0.0);
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for ij
       for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
         {					// vector of operators in basis of Fz
         T1s[m+2] = gen_op((T1[ij].component(2,m)));
//         T1s[m+2] = gen_op((T1[ij].component(2,m)).matrix());
//         T1s[m+2] = gen_op(T1[ij].component(2,m),
//                          sys.get_basis());
         T1s[m+2].Op_base(Heff);
         }
       kl = 0;
       for(int k=0; k<ns-1; k++)		// Sum spin pairs k&l: dipole kl
         for(int l=k+1; l<ns; l++)
           {
           if((ij == kl) && (type >= 0))	// Auto-correlation term
             {
             xi1xi2 = xi1*xi1;			// xi(ij)*xi(ij)
             if(abs(level) > 1)
               for(m=-rank; m<=rank; m++)	// Set up all the reduced spectral
                 {				// density functions J(w-m*Wrf) 
                 mWrf = 1.0*double(m)*Wrflab;
//                 mWrf = -1.0*double(m)*Wrflab;
                 J12[m+2]=J_red_shft(w,mWrf,
 		             hs,taus,c1s,c1s,1);
                 J12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
                 }
              else
                 {				// density functions J(w-m*Wrf) 
                 w0 = 0.0;			// For Level < 2, need only zero,
                 w1 = sys.gamma(i)/GAMMA1H; 	// single, and double quantum
                 w1 *= sys.Omega()*1.0e6; 	// transition frequencies
                 J[0] = xi1xi2*J_reduced(taus,	// J(0)
                                 c1s, c1s,w0,1);
                 J[1] = xi1xi2*J_reduced(taus,	// J(w0)
 				 c1s, c1s,w1,1);
                 J[2] = xi1xi2*J_reduced(taus,	// J(2w0)
         		      c1s,c1s,2.0*w2,1);
                 }
             Rrfmumu(LOp,T1s,T1s,J12,J,w,rank,level,1,het);
             }
           else if((ij != kl) && (type <= 0))	// Cross-correlation term
             {
             xi2 = Re(xi2s.get(k,l));		// Dipolar interaction constant k&l
             xi1xi2 = xi1*xi2;			// xi(ij)*xi(kl)
alphakl = Re(phi.get(k,l));	// Get the alpha Euler angle kl 
betakl = Re(theta.get(k,l));	// Get the beta Euler angle kl
/*            Jcoeffs(c2s, alphakl, betakl, 	// Set the 5 coefficients for kl*/
/*			   0.0,  chi);  // where gammakl = 0 dipolar*/
             EA2.xyz(alphakl, betakl, 0.0);
             Jcoeffs(c2s, EA2, chi); 		// Set the 5 coefficients for kl
             for(m=-2; m<3; m++)		// Put spin tensor for k,l into a
               {				// vector of operators in basis of Fz
               T2s[m+2] = gen_op((T2[kl].component(2,m)));
//               T2s[m+2] = gen_op((T2[kl].component(2,m)).matrix());
//               T2s[m+2] = gen_op(T2[kl].
//                            component(2,m),
//                              sys.get_basis());
               T2s[m+2].Op_base(Heff);
               if(abs(level) > 1)
                 {
                 mWrf = 1.0*double(m)*Wrflab;
//                 mWrf = -1.0*double(m)*Wrflab;
                 J12[m+2]=J_red_shft(w, mWrf,	// Get all reduced spectral densities
 	         hs, taus, c1s, c2s, 1);
                 J12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
                 }
               }
             if(abs(level) < 2)
               {
               wi = sys.gamma(i)/GAMMA1H;	// Need only zero, single and double
               wi *= sys.Omega()*1.0e6; 	// quantum transition frequencies
               wj = sys.gamma(j)/GAMMA1H;
               wj *= sys.Omega()*1.0e6; 	// !! Don't know which to use here!!
               w0 = wi-wj;
               w1 *= wi;
               w2 = wi+wj;
               J[0] = xi1xi2*J_reduced(taus,	// J(0)
                                 c1s, c2s,w0,1);
               J[1] = xi1xi2*J_reduced(taus,	// J(w0)
 				 c1s, c2s,w1,1);
               J[2] = xi1xi2*J_reduced(taus,	// J(2w0)
         		      c1s,c2s,2.0*w2,1);
               }
             Rrfmumu(LOp,T1s,T2s,J12,J,w,rank,level,0,het);
             }
           kl++;				// Increment second dipole
           }
       ij++;					// increment first dipole
       }

     gen_op Op;					// Psuedo destruction of operators
     if(T1s)
       {
       for(int ii=0; ii<5; ii++)
         T1s[ii] = Op;
       T1s = NULL;
       }
     if(T2s)
       {
       for(int jj=0; jj<5; jj++)
         T2s[jj] = Op;
        T2s = NULL;
       }
   return;
   if(A1==NULL && A2==NULL) type=0;		// Compiler likes these used
   }

// -------------------- Spin with Spin Functions -------------------

void Rrfij(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

   {
/* sosi added 3/11/96*/
int het = sys.heteronuclear();                  // Flag if system is heteronuclear
   double xi1, xi2, xi1xi2;			// For specific interaction constants
   coord EA1, EA2;				// For specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
//   int ls = hs*hs;				// Total system Liouville space
   int m;					// z-component of ang. momentum
   double cutoff = 1e-12;
   int rank = 2;
   double mWrf=0;
   double J[3];
   matrix* J12 = NULL;
   if(abs(level) > 1)                           // Needed for higher level computations
     {
     J12 = new matrix[5];
     Heff.eigvals(w);                           // Frequencies w in field rotating frame
     }

   for(int i=0; i<ns; i++)			// Sum over spins i (mu1)
     {
     xi1 = Re(xi1s.get(i,i));			// Get spin i (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       EA1 = (A1[i]).PASys_EA();	 		// Get spin i (mu1) space tensor Euler angles
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for i
       for(m=-2; m<3; m++)			// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of Heff
         T1s[m+2] = gen_op(T1[i].component(2,m));
//         T1s[m+2] = gen_op((T1[i].component(2,m)).matrix());
//         T1s[m+2] = gen_op(T1[i].component(2,m),
//                               sys.get_basis());
         T1s[m+2].Op_base(Heff);
         }
       for(int j=0; j<ns; j++)			// Sum over spins j
         {
         if((i==j) && (type >= 0))		// Auto-correlation term
           {
           xi1xi2 = xi1*xi1;			// xi(mu1,i)*xi(mu1,i)
           if(abs(level) > 1)
             for(m=-rank; m<=rank; m++)		// Set up all the reduced spectral
               {				// density functions J(w-m*Wrf)
               mWrf = -1.0*double(m)*Wrflab;
               J12[m+2]=J_red_shft(w,mWrf,
               hs,taus,c1s,c1s,1);
               J12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
               }
           if(abs(level) <2)
             {
             w0 = 0.0;				// Need only zero, single, & double
             w1 = sys.gamma(i)/GAMMA1H;	// quantum transition frequencies
             w1 *= sys.Omega()*1.0e6;
             J[0] = xi1xi2*J_reduced(taus,	// J(0)
                                 c1s, c1s,w0,1);
             J[1] = xi1xi2*J_reduced(taus,	// J(w0)
                                 c1s, c1s,w1,1);
             J[2] = xi1xi2*J_reduced(taus,	// J(2w0)
                                 c1s,c1s,2.0*w2,1);
             }
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             Rrfmumu(LOp,T1s,T1s,J12,J,w,rank,level,1,het);
           }
         else if((i!=j) && (type<=0))		// Cross-correlation term
           {
           xi2 = Re(xi2s.get(j,j));		// Get spin j (mu2) interaction constant
           xi1xi2 = xi1*xi2;			// xi(mu1,i)*xi(mu2,j)
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             {
             EA2 = (A2[j]).PASys_EA();		// Get spin j (mu2) space tensor Euler angles
             Jcoeffs(c2s, EA2, chi);		// Set the 5 coefficients for j
             for(m=-2; m<3; m++)		// Put spin tensor for j into a
               {				// vector of operators in basis of Heff
               T2s[m+2] = gen_op(T2[j].component(2,m));
//               T2s[m+2] = gen_op((T2[j].component(2,m)).matrix());
//               T2s[m+2] = gen_op(T2[j].component(2,m),
//                                   sys.get_basis());
               T2s[m+2].Op_base(Heff);
               if(abs(level) > 1)
                 {
                 mWrf = -1.0*double(m)*Wrflab;
                 J12[m+2]=J_red_shft(w, mWrf,   // Get all reduced spectral densities
                         hs, taus, c1s, c2s, 1);
                 J12[m+2] *= complex(xi1xi2);   // Scale by the interaction constants
                 }
               }
             if(abs(level) <2)
               {
               wi = sys.gamma(i)/GAMMA1H; 	// Need only zero, single, & double
               wi *= sys.Omega()*1.0e6; 	// quantum transition frequencies
               wj = sys.gamma(j)/GAMMA1H;
               wj *= sys.Omega()*1.0e6;
               w0 = wi-wj;
               w1 *= wi;			// !! Don't know which to use here!!
               w2 = wi+wj;
               J[0] = xi1xi2*J_reduced(taus,	// J(0)
                                 c1s, c2s,w0,1);
               J[1] = xi1xi2*J_reduced(taus,	// J(w0)
                                 c1s, c2s,w1,1);
               J[2] = xi1xi2*J_reduced(taus,	// J(2w0)
                                 c1s,c2s,2.0*w2,1);
               }
// sosi - added this to insure that heteronuclear terms are not included
//if(sys.element(i) == sys.element(j))
             Rrfmumu(LOp,T1s,T2s,J12,J,w,rank,level,0,het);
             }
           }
         }					// Increment second spin
       }
     }						// Increment first spin
   return;
   }


// --------------- Spin-Pair with Spin Functions -------------------

void Rrfijk(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

   {
// sosi added 3/11/96
int het = sys.heteronuclear();                  // Flag if system is heteronuclear
   int rank=2;
// sosi: these three variables are still dipole specific!!
matrix theta = sys.thetas();			// Get dipole theta values (radians) 
matrix phi = sys.phis(); 			// Get dipole phi values (radians)
double alphaij, betaij;				// Two Euler angles for ij dipole
   double xi1, xi2, xi1xi2=0;			// For specific interaction constants
   coord EA1, EA2;				// For specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
//   int ls = hs*hs;				// Total system Liouville space
   int m;					// z-component of ang. momentum
   double cutoff = 1e-12;
   double mWrf=0;
   double J[3];
   matrix* J12 = NULL;
   if(abs(level) > 1)                           // Needed for higher level computations
     {
     J12 = new matrix[5];
     Heff.eigvals(w);                           // Frequencies w in field rotating frame
     }
   int ij = 0;					// Sum over all dipolar pairs
   for(int i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)
       {
       xi1 = Re(xi1s.get(i,j));			// Dipolar interaction constant i&j
alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
//       Jcoeffs(c1s, alphaij, betaij, 0.0, chi);	// Set the 5 coefficients for ij
EA1.xyz(alphaij, betaij, 0.0);
       Jcoeffs(c1s, EA1, chi);
       for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
         {					// vector of operators in basis of Heff
         T1s[m+2] = gen_op(T1[ij].component(2,m));
//         T1s[m+2] = gen_op((T1[ij].component(2,m)).matrix());
//         T1s[m+2] = gen_op(T1[ij].component(2,m),
//                          sys.get_basis());
         T1s[m+2].Op_base(Heff);
         }
       for(int k=0; k<ns; k++)			// Sum over spins k
         {
         xi2 = Re(xi2s.get(k,k));		// Get spin k (mu2) interaction constant
         xi1xi2 = xi1*xi2;			// xi(mu1,i)*xi(mu2,j)
         if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
           {
           EA2 = (A2[k]).PASys_EA();		// Get spin k (mu2) space tensor Euler angles
           Jcoeffs(c2s, EA2, chi);		// Set the 5 coefficients for k
           for(m=-2; m<3; m++)			// Put spin tensor for k into a
             {					// vector of operators in basis of Heff
             T2s[m+2] = gen_op(T2[k].component(2,m));
//             T2s[m+2] = gen_op((T2[k].component(2,m)).matrix());
//             T2s[m+2] = gen_op(T2[k].component(2,m),
//                                 sys.get_basis());
             T2s[m+2].Op_base(Heff);
             }
           if(abs(level) > 1)
             for(m=-rank; m<=rank; m++)		// Set up all the reduced spectral
               {				// density functions J(w-m*Wrf) 
               mWrf = -1.0*double(m)*Wrflab;
               J12[m+2]=J_red_shft(w,mWrf,
 		             hs,taus,c1s,c1s,1);
               J12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
               }
           else
             {
             wi = sys.gamma(i)/GAMMA1H; 	// Need only zero, single, & double
             wi *= sys.Omega()*1.0e6;	 	// quantum transition frequencies
             wj = sys.gamma(j)/GAMMA1H;
             wj *= sys.Omega()*1.0e6;
             w0 = wi-wj;
             w1 *= wi;				// !! Don't know which to use here!!
             w2 = wi+wj;
             }
           Rrfmumu(LOp,T1s,T2s,J12,J,w,rank,level,0,het);
           }
         } 					// Increment spin k (mu2)
       }					// Increment pair ij (mu1)
   return;
   if(A1==NULL) type=0;				// Compiler likes these used
   }


void Rrfkij(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

   {
// sosi added 3-11-96
int het = sys.heteronuclear();                  // Flag if system is heteronuclear
// sosi: these three variables are still dipole specific!!
   int rank=2;
matrix theta = sys.thetas();			// Get dipole theta values (radians) 
matrix phi = sys.phis(); 			// Get dipole phi values (radians)
double alphaij, betaij;			// Two Euler angles for ij dipole

   double xi1, xi2, xi1xi2;			// For specific interaction constants
   coord EA1, EA2;				// For specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
//   int ls = hs*hs;				// Total system Liouville space
   int m;					// z-component of ang. momentum
   int ij=0;					// z-component of ang. momentum
   double cutoff = 0;
//   double cutoff = 1e-12;
   double mWrf=0;
   double J[3];
   matrix* J12 = NULL;
   if(abs(level) > 1)                           // Needed for higher level computations
     {
     J12 = new matrix[5];
     Heff.eigvals(w);                           // Frequencies w in field rotating frame
     }
   for(int k=0; k<ns; k++)			// Sum over spins k (mu1)
     {
     xi1 = Re(xi1s.get(k,k));			// Get spin k (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       EA1 = (A1[k]).PASys_EA(); 		// Get spin k (mu1) space tensor Euler angles
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for k
       for(m=-2; m<3; m++)			// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of effH
         T1s[m+2] = gen_op(T1[k].component(2,m));
//         T1s[m+2] = gen_op((T1[k].component(2,m)).matrix());
//         T1s[m+2] = gen_op(T1[k].component(2,m),
//                               sys.get_basis());
         T1s[m+2].Op_base(Heff);
         }
       ij = 0;					// Set dipole count to zero
       for(int i=0; i<ns-1; i++)		// Sum spin pairs i&j: dipole ij (mu2)
         for(int j=i+1; j<ns; j++)
           {
           xi2 = Re(xi2s.get(i,j));		// Dipolar interaction constant k&l
           xi1xi2 = xi1*xi2;			// xi(k)*xi(ij)
alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
//Jcoeffs(c2s,alphaij,betaij,0.0,chi);	// Set the 5 coefficients for ij
EA2.xyz(alphaij, betaij, 0.0);
           Jcoeffs(c2s,EA2,chi);	// Set the 5 coefficients for ij
           for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
             {					// vector of operators in basis of Heff
             T2s[m+2] = gen_op(T2[ij].component(2,m));
//             T2s[m+2] = gen_op((T2[ij].component(2,m)).matrix());
//             T2s[m+2] = gen_op(T2[ij].
//                            component(2,m),
//                              sys.get_basis());
               T2s[m+2].Op_base(Heff);
             }
           if(abs(level) > 1)
             for(m=-rank; m<=rank; m++)		// Set up all the reduced spectral
               {				// density functions J(w-m*Wrf) 
               mWrf = -1.0*double(m)*Wrflab;
               J12[m+2]=J_red_shft(w,mWrf,
 		             hs,taus,c1s,c1s,1);
               J12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
               }
           else
             {
             wi = sys.gamma(i)/GAMMA1H;	// Need only zero, single and double
             wi *= sys.Omega()*1.0e6;	 	// quantum transition frequencies
             wj = sys.gamma(j)/GAMMA1H;
             wj *= sys.Omega()*1.0e6;	 	// !! Don't know which to use here!!
             w0 = wi-wj;
             w1 *= wi;
             w2 = wi+wj;
             }
           Rrfmumu(LOp,T1s,T2s,J12,J,w,rank,level,0,het);
           ij++;
           }					// Increment second dipole (mu2)
       }
     }						// Increment spin (mu1)
   gen_op Op;					// Psuedo destruction of operators
   if(T1s)
     {
     for(int ii=0; ii<5; ii++)
       T1s[ii] = Op;
     T1s = NULL;
     }
   if(T2s)
     {
     for(int jj=0; jj<5; jj++)
       T2s[jj] = Op;
     T2s = NULL;
     }
   return;
   if(A2==NULL) type=0;				// Compiler likes these used
   }


// ______________________________________________________________________
// *RELAXATION SUPEROPERATOR GENERAL FUNCTIONS   DYNAMIC FREQUENCY SHIFT*
// ______________________________________________________________________


// --------- Spin-Pair with Spin-Pair Interaction Functions -------------

   void Rrfijklds(super_op &LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

   {
// sosi: this is still dipole specific!!
matrix theta = sys.thetas();			// Dipole theta values (radians) 
matrix phi = sys.phis(); 			// Dipole phi values (radians)
double alphaij, betaij;			// Two Euler angles for ij dipole
double alphakl, betakl;			// Two Euler angles for kl dipole
   double xi1, xi2, xi1xi2=0;			// Specific interaction consts.
   coord EA1, EA2;				// Specific Euler angles
   double c1s[5];				// Set up 5 coefficients int. 1
   double c2s[5];				// Set up 5 coefficients int. 2
   gen_op *T1s;					// Compiler hates gen_op T1s[5]
   gen_op *T2s;					// These for spin tensor comps.
   T1s = new gen_op[5];
   T2s = new gen_op[5];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for Q's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
//   int ls = hs*hs;				// Total system Liouville space
   int m;					// z-comp. spin ang. momentum
   int ij=0, kl=0;				// Spin pair indices, 2 mechanisms

//
//		    Prepare the Interaction Constants
//		     and Spectral Density Components

  double mWrf=0;
  int rank = 2;
  double Q[3];
  matrix* Q12 = NULL;                          // Matrices of spectral densities
  if(abs(level) > 1)                           // Needed for higher level computations
    {
    Q12 = new matrix[5];
    Heff.eigvals(w);                           // Frequencies w in field rotating frame
    }

//	      Start the Relaxation Superoperator Computation

   for(int i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)
       {
       xi1 = Re(xi1s.get(i,j));			// 1st interaction constant, i&j
alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
/*      Jcoeffs(c1s, alphaij, betaij, 0.0, chi);	// Set the 5 coefficients for ij*/
EA1.xyz(alphaij, betaij, 0.0);
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for ij
       for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
         {					// vector of operators in basis of Fz
         T1s[m+2] = gen_op(T1[ij].component(2,m));
//         T1s[m+2] = gen_op((T1[ij].component(2,m)).matrix());
//         T1s[m+2] = gen_op(T1[ij].component(2,m),
//                          sys.get_basis());
         T1s[m+2].Op_base(Heff);
         }
       kl = 0;
       for(int k=0; k<ns-1; k++)		// Sum spin pairs k&l: dipole kl
         for(int l=k+1; l<ns; l++)
           {
           if((ij == kl) && (type >= 0))	// Auto-correlation term
             {
             xi1xi2 = xi1*xi1;			// xi(ij)*xi(ij)
             if(abs(level) > 1)
               for(m=-rank; m<=rank; m++)	// Set up all the reduced spectral
                 {				// density functions Q(w-m*Wrf) 
                 mWrf = -1.0*double(m)*Wrflab;
                 Q12[m+2]=Q_red_shft(w,mWrf,
 		             hs,taus,c1s,c1s,1);
                 Q12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
                 }
              else
                 {				// density functions Q(w-m*Wrf) 
                 w0 = 0.0;			// For Level < 2, need only zero,
                 w1 = sys.gamma(i)/GAMMA1H; 	// single, and double quantum
                 w1 *= sys.Omega()*1.0e6; 	// transition frequencies
                 Q[0] = xi1xi2*Q_reduced(taus,	// Q(0)
                                 c1s, c1s,w0,1);
                 Q[1] = xi1xi2*Q_reduced(taus,	// Q(w0)
 				 c1s, c1s,w1,1);
                 Q[2] = xi1xi2*Q_reduced(taus,	// Q(2w0)
         		      c1s,c1s,2.0*w2,1);
                 }
             Rrfmumu(LOp,T1s,T1s,Q12,Q,w,rank,level,1);
             }
           else if((ij != kl) && (type <= 0))	// Cross-correlation term
             {
             xi2 = Re(xi2s.get(k,l));		// Dipolar interaction constant k&l
             xi1xi2 = xi1*xi2;			// xi(ij)*xi(kl)
alphakl = Re(phi.get(k,l));	// Get the alpha Euler angle kl 
betakl = Re(theta.get(k,l));	// Get the beta Euler angle kl
/*            Jcoeffs(c2s, alphakl, betakl, 	// Set the 5 coefficients for kl
//			   0.0,  chi);  // where gammakl = 0 dipolar*/
             EA2.xyz(alphakl, betakl, 0.0);
             Jcoeffs(c2s, EA2, chi); 		// Set the 5 coefficients for kl
             for(m=-2; m<3; m++)		// Put spin tensor for k,l into a
               {				// vector of operators in basis of Fz
               T2s[m+2] = gen_op(T2[kl].component(2,m));
//               T2s[m+2] = gen_op((T2[kl].component(2,m)).matrix());
//               T2s[m+2] = gen_op(T2[kl].
//                            component(2,m),
//                              sys.get_basis());
               T2s[m+2].Op_base(Heff);
               if(abs(level) > 1)
                 {
                 mWrf = -1.0*double(m)*Wrflab;
                 Q12[m+2]=Q_red_shft(w, mWrf,	// Get all reduced spectral densities
 	         hs, taus, c1s, c2s, 1);
                 Q12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
                 }
               }
             if(abs(level) < 2)
               {
               wi = sys.gamma(i)/GAMMA1H;	// Need only zero, single and double
               wi *= sys.Omega()*1.0e6; 	// quantum transition frequencies
               wj = sys.gamma(j)/GAMMA1H;
               wj *= sys.Omega()*1.0e6; 	// !! Don't know which to use here!!
               w0 = wi-wj;
               w1 *= wi;
               w2 = wi+wj;
               Q[0] = xi1xi2*Q_reduced(taus,	// Q(0)
                                 c1s, c2s,w0,1);
               Q[1] = xi1xi2*Q_reduced(taus,	// Q(w0)
 				 c1s, c2s,w1,1);
               Q[2] = xi1xi2*Q_reduced(taus,	// Q(2w0)
         		      c1s,c2s,2.0*w2,1);
               }
             Rrfmumu(LOp,T1s,T2s,Q12,Q,w,rank,level,0);
             }
           kl++;				// Increment second dipole
           }
       ij++;					// increment first dipole
       }

     gen_op Op;					// Psuedo destruction of operators
     if(T1s)
       {
       for(int ii=0; ii<5; ii++)
         T1s[ii] = Op;
       T1s = NULL;
       }
     if(T2s)
       {
       for(int jj=0; jj<5; jj++)
         T2s[jj] = Op;
        T2s = NULL;
       }
   return;
   if(A2==NULL && A1==NULL) type=0;		// Compiler likes these used
   }

// -------------------- Spin with Spin Functions -------------------

void Rrfijds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

   {
   double xi1, xi2, xi1xi2;			// For specific interaction constants
   coord EA1, EA2;				// For specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for Q's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
/*   int ls = hs*hs;				// Total system Liouville space*/
   int m;					// z-component of ang. momentum
   double cutoff = 1e-12;
   int rank = 2;
   double mWrf=0;
   double Q[3];
   matrix* Q12 = NULL;
   if(abs(level) > 1)                           // Needed for higher level computations
     {
     Q12 = new matrix[5];
     Heff.eigvals(w);                           // Frequencies w in field rotating frame
     }

   for(int i=0; i<ns; i++)			// Sum over spins i (mu1)
     {
     xi1 = Re(xi1s.get(i,i));			// Get spin i (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       EA1 = (A1[i]).PASys_EA();	 		// Get spin i (mu1) space tensor Euler angles
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for i
       for(m=-2; m<3; m++)			// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of Heff
         T1s[m+2] = gen_op(T1[i].component(2,m));
//         T1s[m+2] = gen_op((T1[i].component(2,m)).matrix());
//         T1s[m+2] = gen_op(T1[i].component(2,m),
//                               sys.get_basis());
         T1s[m+2].Op_base(Heff);
         }
       for(int j=0; j<ns; j++)			// Sum over spins j
         {
         if((i==j) && (type >= 0))		// Auto-correlation term
           {
           xi1xi2 = xi1*xi1;			// xi(mu1,i)*xi(mu1,i)
           if(abs(level) > 1)
             for(m=-rank; m<=rank; m++)		// Set up all the reduced spectral
               {				// density functions Q(w-m*Wrf)
               mWrf = -1.0*double(m)*Wrflab;
               Q12[m+2]=Q_red_shft(w,mWrf,
               hs,taus,c1s,c1s,1);
               Q12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
               }
           if(abs(level) <2)
             {
             w0 = 0.0;				// Need only zero, single, & double
             w1 = sys.gamma(i)/GAMMA1H;	// quantum transition frequencies
             w1 *= sys.Omega()*1.0e6;
             Q[0] = xi1xi2*Q_reduced(taus,	// Q(0)
                                 c1s, c1s,w0,1);
             Q[1] = xi1xi2*Q_reduced(taus,	// Q(w0)
                                 c1s, c1s,w1,1);
             Q[2] = xi1xi2*Q_reduced(taus,	// Q(2w0)
                                 c1s,c1s,2.0*w2,1);
             }
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             Rrfmumu(LOp,T1s,T1s,Q12,Q,w,rank,level,1);
           }
         else if((i!=j) && (type<=0))		// Cross-correlation term
           {
           xi2 = Re(xi2s.get(j,j));		// Get spin j (mu2) interaction constant
           xi1xi2 = xi1*xi2;			// xi(mu1,i)*xi(mu2,j)
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             {
             EA2 = (A2[j]).PASys_EA();		// Get spin j (mu2) space tensor Euler angles
             Jcoeffs(c2s, EA2, chi);		// Set the 5 coefficients for j
             for(m=-2; m<3; m++)		// Put spin tensor for j into a
               {				// vector of operators in basis of Heff
               T2s[m+2] = gen_op(T2[j].component(2,m));
//               T2s[m+2] = gen_op((T2[j].component(2,m)).matrix());
//               T2s[m+2] = gen_op(T2[j].component(2,m),
//                                   sys.get_basis());
               T2s[m+2].Op_base(Heff);
               if(abs(level) > 1)
                 {
                 mWrf = -1.0*double(m)*Wrflab;
                 Q12[m+2]=Q_red_shft(w, mWrf,   // Get all reduced spectral densities
                         hs, taus, c1s, c2s, 1);
                 Q12[m+2] *= complex(xi1xi2);   // Scale by the interaction constants
                 }
               }
             if(abs(level) <2)
               {
               wi = sys.gamma(i)/GAMMA1H; 	// Need only zero, single, & double
               wi *= sys.Omega()*1.0e6; 	// quantum transition frequencies
               wj = sys.gamma(j)/GAMMA1H;
               wj *= sys.Omega()*1.0e6;
               w0 = wi-wj;
               w1 *= wi;			// !! Don't know which to use here!!
               w2 = wi+wj;
               Q[0] = xi1xi2*Q_reduced(taus,	// Q(0)
                                 c1s, c2s,w0,1);
               Q[1] = xi1xi2*Q_reduced(taus,	// Q(w0)
                                 c1s, c2s,w1,1);
               Q[2] = xi1xi2*Q_reduced(taus,	// Q(2w0)
                                 c1s,c2s,2.0*w2,1);
               }
             Rrfmumu(LOp,T1s,T2s,Q12,Q,w,rank,level,0);
             }
           }
         }					// Increment second spin
       }
     }						// Increment first spin
   return;
   }


// --------------- Spin-Pair with Spin Functions -------------------

void Rrfijkds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

   {
   int rank=2;
// sosi: these three variables are still dipole specific!!
matrix theta = sys.thetas();			// Get dipole theta values (radians) 
matrix phi = sys.phis(); 			// Get dipole phi values (radians)
double alphaij, betaij;				// Two Euler angles for ij dipole
   double xi1, xi2, xi1xi2=0;			// For specific interaction constants
   coord EA1, EA2;				// For specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
/*   int ls = hs*hs;				// Total system Liouville space*/
   int m;					// z-component of ang. momentum
   double cutoff = 1e-12;
   double mWrf=0;
   double J[3];
   matrix* J12 = NULL;
   if(abs(level) > 1)                           // Needed for higher level computations
     {
     J12 = new matrix[5];
     Heff.eigvals(w);                           // Frequencies w in field rotating frame
     }
   int ij = 0;					// Sum over all dipolar pairs
   for(int i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)
       {
       xi1 = Re(xi1s.get(i,j));			// Dipolar interaction constant i&j
alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
/*       Jcoeffs(c1s, alphaij, betaij, 0.0, chi);	// Set the 5 coefficients for ij*/
EA1.xyz(alphaij, betaij, 0.0);
       Jcoeffs(c1s, EA1, chi);
       for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
         {					// vector of operators in basis of Heff
         T1s[m+2] = gen_op(T1[ij].component(2,m));
//         T1s[m+2] = gen_op((T1[ij].component(2,m)).matrix());
//         T1s[m+2] = gen_op(T1[ij].component(2,m),
//                          sys.get_basis());
         T1s[m+2].Op_base(Heff);
         }
       for(int k=0; k<ns; k++)			// Sum over spins k
         {
         xi2 = Re(xi2s.get(k,k));		// Get spin k (mu2) interaction constant
         xi1xi2 = xi1*xi2;			// xi(mu1,i)*xi(mu2,j)
         if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
           {
           EA2 = (A2[k]).PASys_EA();		// Get spin k (mu2) space tensor Euler angles
           Jcoeffs(c2s, EA2, chi);		// Set the 5 coefficients for k
           for(m=-2; m<3; m++)			// Put spin tensor for k into a
             {					// vector of operators in basis of Heff
             T2s[m+2] = gen_op(T2[k].component(2,m));
//             T2s[m+2] = gen_op((T2[k].component(2,m)).matrix());
//             T2s[m+2] = gen_op(T2[k].component(2,m),
//                                 sys.get_basis());
             T2s[m+2].Op_base(Heff);
             }
           if(abs(level) > 1)
             for(m=-rank; m<=rank; m++)		// Set up all the reduced spectral
               {				// density functions J(w-m*Wrf) 
               mWrf = -1.0*double(m)*Wrflab;
// sosi 3-10-96 why isnt this Q and c1s with c2s??
               J12[m+2]=J_red_shft(w,mWrf,
 		             hs,taus,c1s,c1s,1);
               J12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
               }
           else
             {
             wi = sys.gamma(i)/GAMMA1H; 	// Need only zero, single, & double
             wi *= sys.Omega()*1.0e6;	 	// quantum transition frequencies
             wj = sys.gamma(j)/GAMMA1H;
             wj *= sys.Omega()*1.0e6;
             w0 = wi-wj;
             w1 *= wi;				// !! Don't know which to use here!!
             w2 = wi+wj;
             }
           Rrfmumu(LOp,T1s,T2s,J12,J,w,rank,level,0);
           }
         } 					// Increment spin k (mu2)
       }					// Increment pair ij (mu1)
   return;
   if(A1==NULL) type=0;				// Compiler likes these used
   }


void Rrfkijds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

   {
// sosi: these three variables are still dipole specific!!
   int rank=2;
matrix theta = sys.thetas();			// Get dipole theta values (radians) 
matrix phi = sys.phis(); 			// Get dipole phi values (radians)
double alphaij, betaij;			// Two Euler angles for ij dipole

   double xi1, xi2, xi1xi2;			// For specific interaction constants
   coord EA1, EA2;				// For specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
/*   int ls = hs*hs;				// Total system Liouville space*/
   int m;					// z-component of ang. momentum
   int ij=0;					// z-component of ang. momentum
   double cutoff = 0;
//   double cutoff = 1e-12;
   double mWrf=0;
   double J[3];
   matrix* J12 = NULL;
   if(abs(level) > 1)                           // Needed for higher level computations
     {
     J12 = new matrix[5];
     Heff.eigvals(w);                           // Frequencies w in field rotating frame
     }
   for(int k=0; k<ns; k++)			// Sum over spins k (mu1)
     {
     xi1 = Re(xi1s.get(k,k));			// Get spin k (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       EA1 = (A1[k]).PASys_EA(); 		// Get spin k (mu1) space tensor Euler angles
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for k
       for(m=-2; m<3; m++)			// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of effH
         T1s[m+2] = gen_op(T1[k].component(2,m));
//         T1s[m+2] = gen_op((T1[k].component(2,m)).matrix());
//         T1s[m+2] = gen_op(T1[k].component(2,m),
//                               sys.get_basis());
         T1s[m+2].Op_base(Heff);
         }
       ij = 0;					// Set dipole count to zero
       for(int i=0; i<ns-1; i++)		// Sum spin pairs i&j: dipole ij (mu2)
         for(int j=i+1; j<ns; j++)
           {
           xi2 = Re(xi2s.get(i,j));		// Dipolar interaction constant k&l
           xi1xi2 = xi1*xi2;			// xi(k)*xi(ij)
alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
//Jcoeffs(c2s,alphaij,betaij,0.0,chi);	// Set the 5 coefficients for ij
EA2.xyz(alphaij, betaij, 0.0);
           Jcoeffs(c2s,EA2,chi);	// Set the 5 coefficients for ij
           for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
             {					// vector of operators in basis of Heff
             T2s[m+2] = gen_op(T2[ij].component(2,m));
//             T2s[m+2] = gen_op((T2[ij].component(2,m)).matrix());
//             T2s[m+2] = gen_op(T2[ij].
//                            component(2,m),
//                              sys.get_basis());
               T2s[m+2].Op_base(Heff);
             }
           if(abs(level) > 1)
             for(m=-rank; m<=rank; m++)		// Set up all the reduced spectral
               {				// density functions J(w-m*Wrf) 
               mWrf = -1.0*double(m)*Wrflab;
// sosi 3-10-96 why isnt this Q and c1s with c2s??
               J12[m+2]=J_red_shft(w,mWrf,
 		             hs,taus,c1s,c1s,1);
               J12[m+2] *= complex(xi1xi2);	// Scale by the interaction constants
               }
           else
             {
             wi = sys.gamma(i)/GAMMA1H;	// Need only zero, single and double
             wi *= sys.Omega()*1.0e6;	 	// quantum transition frequencies
             wj = sys.gamma(j)/GAMMA1H;
             wj *= sys.Omega()*1.0e6;	 	// !! Don't know which to use here!!
             w0 = wi-wj;
             w1 *= wi;
             w2 = wi+wj;
             }
           Rrfmumu(LOp,T1s,T2s,J12,J,w,rank,level,0);
           ij++;
           }					// Increment second dipole (mu2)
       }
     }						// Increment spin (mu1)
   gen_op Op;					// Psuedo destruction of operators
   if(T1s)
     {
     for(int ii=0; ii<5; ii++)
       T1s[ii] = Op;
     T1s = NULL;
     }
   if(T2s)
     {
     for(int jj=0; jj<5; jj++)
       T2s[jj] = Op;
     T2s = NULL;
     }
   return;
   if(A2==NULL) type=0;				// Compiler likes these used
   }


// ______________________________________________________________________
// ***************** STEADY STATE DENSITY MATRIX FUNCTIONS **************
// ______________________________________________________________________

// ----------------- Level 4 by Double Commutators ----------------------

gen_op sigma_ss(spin_system& sys, super_op& L, super_op& R)
// sosi - THIS ROUTINE NEEDS TO BE PATCHED!!! IT ASSUMES R & L ARE IN
//        THE SAME LIOUVILLE BASIS CURRENTLY!!!

	// Input		sys   : Spin system
	// 			L     : Full Liouvillian (rad/sec)
	//			R     : Relaxation Superoperator (rad/sec)
        // Return		sss   : Steady state density matrix
	// Note			      : Careful R & L units, usually Ho is Hz!
	// Note			      : L is singular, thus steady state must
	//				be determined in a reduced Liouville
	//				space.  Best done with matrix routines
	//				not (code simple) superop. functions!
	// Note			      : Algorithm needs density matrix trace = 1
	//				Function insures this then returns trace
	//				back to the original value.

//	L|s  > = R|s  >    -------> [L - L|S><E|]|s  > = R|s  > - |L1>
//         ss       eq		                   ss       eq

  {
  gen_op seq;
  double d = 1.e-9; 
  R.LOp_base(L);			// Insure R in Liouville base of L
  if(R.below(d))			// If R==0, |s  > = 0
    {					//            ss
    std::cout << "\n\n\tWarning relax_rf: Steady State Does Not Exist";
    return seq;
    }
  else if(L == R)			// If R==L, |s  > = s  >
    return sigma_eq(sys);		//	      ss     eq	
  else
    {
    basis Lbs = L.get_basis();		// Store the Hilbert basis of L
    seq = sigma_eq(sys);		// Equilibrium density matrix
    int hs = seq.dim();			// Get the Hilbert space size
    int ls = hs*hs;			// Compute the Liouville space
    matrix Emx(hs,hs,i_matrix_type);	// Temporary Identity matrix
    matrix trace1 = Emx/complex(hs);	// Scale by dimension
    seq += trace1;			// Set trace to 1 for this
    seq.Op_base(Lbs);			// Put operator into L basis

//	Calculate Equation Right Hand Side of Steady State Equation
//
//			  '
//			|s  > = R|s  > - |L1>
//			  eq       eq

    matrix L1, ssp, sub_ssp;
    L1=(L.get_mx()).get_block(0,0,ls,1);// L1 matrix = 1st column of L
    ssp=((R*seq).get_mx()).resize(ls,1);// Modified seq, basis of R
    ssp -= L1;				// |ssp> = R|seq> - |L1>
    sub_ssp = ssp.get_block(1,0,ls-1,1);// Reduced space seq (as matrix)

//	Calculate Equation Left Hand Side of Steady State Equation
//
//			X = L - L|S><E|
//

    matrix Smx(hs, hs, 0.0);		// Temproray Null matrix
    Smx.put(1.0, 0, 0);			// 1st element of Smx set to 1 
    gen_op S(Smx, Lbs);			// Create the S operator
    gen_op E(Emx, Lbs);			// Create the S operator
    super_op SE(S,E);			// SE = |S><E|
    super_op X = L;
    X -= L*SE; 				// X = L - L|S><E|
    matrix sub_X((X.get_mx()).		// Reduced space X
	    get_block(1,1,ls-1,ls-1));

// Form the LU Decomposition of the Reduced & Modified Liouville Matrix
  
/*                      '                    -1  '
  	     X|s  > = |s  >  ---->  |s  > = X  |s  >
             - -ss     -eq           -ss    -   -eq
  
                       '                       -1  '
           LU|s  > = |s  >  ---->  |s  > = (LU)  |s  >
           -- -ss     -eq           -ss     --    -eq
*/

    matrix ALU = sub_X;
    int *indx;
    indx = new int[ls-1];
    LU_decomp(ALU,indx);

/*		Calculate the Steady State Density Matrix
  
                         -1  '               trace
  	     |s  > = (LU)  |s  >  ;   |s  > ------> |s  >
              -ss     --   -eq         -ss            ss
*/


    matrix sub_sss = sub_ssp;		// Steady state in subspace
    LU_backsub(ALU,indx,sub_sss);	// Back solve for steady state
    matrix sssmx(hs, hs, 0.0);		// Reform the full steady state
    int k=0;				// density matrix in Hilbert
    for(int i=0; i<hs; i++)		// space.
      for(int j=0; j<hs; j++)
        if(i!=0 || j!=0)		// All but first element
          {
          sssmx.put(sub_sss.get(k,0),i,j);
          k++;
          }
    sssmx.put(1.0-trace(sssmx),0,0);	// First element from trace
    sssmx.resize(hs,hs);		// Back to a square array
    sssmx -= trace1;			// Set trace back to 0
    gen_op sss(sssmx, seq.get_basis()); // Put this into a gen. op.
    delete [] indx;
    return sss;
    }
  }


gen_op sigma_ss_it(spin_system& sys, super_op& L, super_op& Heff, super_op& R)

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

  {
  gen_op seq = sigma_eq(sys);		// Equilibrium density matrix
  int hs = seq.dim();			// Get the Hilbert space size
  int ls = hs*hs;			// Compute the Liouville space

//	Calculate Equation Right Hand Side of Steady State Equation
//
//			  '
//			|s  > = R|s  > - |L1>
//			  eq       eq

  matrix L1, ssp, sub_ssp, sub_seq;
  L1=(L.get_mx()).get_block(0,0,ls,1);	// L1 matrix is first column of L
  ssp=((R*seq).get_mx()).resize(ls,1);	// Modified seq, basis of R
  ssp -= L1;				// |ssp> = R|seq> - |L1>
  sub_ssp = ssp.get_block(1,0,ls-1,1);	// Reduced space seq (as matrix)
  ssp=(seq.get_mx()).resize(ls,1);	// Only done for next line
  sub_seq=ssp.get_block(1,0,ls-1,1);	// sub_seq used later in algorithm

//	Calculate Equation Left Hand Side of Steady State Equation
  
//  			X = L - L|S><E|
  

  basis Lbs = L.get_basis();		// Store the Hilbert basis of L
  matrix Smx(hs, hs, 0.0);		// Temproray Null matrix
  Smx.put(1.0, 0, 0);			// 1st element of Smx set to 1 
  gen_op S(Smx, Lbs);			// Create the S operator
  matrix Emx(hs,hs,i_matrix_type);	// Temproray Null matrix
  gen_op E(Emx, Lbs);			// Create the S operator
  super_op SE(S,E);			// SE = |S><E|
  super_op X = L;
  X -= L*SE;		 		// X = L - L|S><E|
  matrix sub_X((X.get_mx()).		// Reduced space X
	    get_block(1,1,ls-1,ls-1));

/*	Compute MX Subspace Array for Iterative Solution
  
                             -1        
  	    	      MX = [M  ][M - X]
                      --    -    -   -
*/

  matrix sub_M = (R.get_mx()).		// Attempt to use R+I for M
              get_block(1,1,ls-1,ls-1);
  matrix sub_I(ls-1,ls-1,i_matrix_type);
std::cout << "\nThis the the sub matrix?" << sub_M;
std::cout << "\nIs this the identity matrix?" << sub_I;
//  sub_M += sub_I;			// !!!! This doesn't work!!!  GAMMA error!!
for(int jj=0; jj<ls-1; jj++)
  sub_M.put(sub_M(jj,jj) + 1.0, jj,jj);
std::cout << "\nDoes this have the identity matrix added?" << sub_M;
  matrix sub_Minv, sub_MX;
  sub_Minv = invert_it(sub_M);
  sub_MX = sub_Minv*(sub_M-sub_X);

/*	Iterate for the Subspace Steady State Density Matrix
  
              -1                    -1  '                   -1  '
   |s(N)> = [M  ][M - X]|s(N-1)> + M  |s  > = MX|s(N-1)> + M  |s  >
    -ss      -    -   -  -ss       -   -eq    -- -ss       -   -eq
*/

  matrix b, sub_sss;
  b = sub_Minv * sub_ssp;
  sub_sss=solve_it(sub_MX, sub_seq, b,4);// Iterate to sub steady state
  matrix sssmx(hs, hs, 0.0);		// Reform the steady state
  int k=0;				// density matrix in Hilbert
  for(int i=0; i<hs; i++)		// space.
    for(int j=0; j<hs; j++)
      if(i!=0 || j!=0)			// All but first element
        {
        sssmx.put(sub_sss.get(k,0),i,j);
        k++;
        }
  sssmx.put(1.0-trace(sssmx),0,0);	// First element from trace
  gen_op sss(sssmx.resize(hs,hs),	// Put this into a gen. op.
  		      seq.get_basis());
  return(sss);
  if(Heff.dim()) hs=0;			// Compiler likes Heff used
  }


#endif /* __RELAX_RF_CC__ */

