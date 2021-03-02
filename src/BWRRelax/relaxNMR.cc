/* relaxNMR.cc ***********************************************************
**									**
** 	                           G A M M A				**
**									**
**	NMR Relaxation Related Functions            Implementation 	**
**						 			**
**	Copyright (c) 1991, 1992, 1993		 			**
**	Scott A. Smith				 			**
**									**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fuer physikalische Chemie		 			**
**	8092 Zurich / Switzerland		 			**
**									**
**	University of California, Santa Barbara				**
**	Department of Chemistry						**
**	Santa Barbara CA. 93106 USA					**
**						 			**
**      $Header: $
**						 			**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**						 			**
** The following functions provide easy access to several generalized	**
** relaxation superoperator matrices and other related quantities.	**
**						 			**
*************************************************************************/

#ifndef RelaxMR_cc_			// Is this file already included?
#  define RelaxMR_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <BWRRelax/relaxNMR.h>		// Include the header file
#include <Basics/Gconstants.h>		// Include GAMMA1H definitio
#include <stdlib.h>
#include <string>			// Include libstdc++ strings

using std::string;			// Using libstdc++ strings

// ____________________________________________________________________________
// i                RELAXATION SUPEROPERATOR ERROR HANDLING
// ____________________________________________________________________________

	// Input		eidx	: Error index
	//			noret   : Flag for return (0=linefeed)
	// Output		none 	: Error Message Output

void RlxNMRerror(int eidx, int noret)
  {
  string hdr("RelaxNMR");
  switch(eidx)
    {
    case 0:
    default:GAMMAerror(hdr, eidx, noret); break;
    }
  }

void volatile RlxNMRfatal(int eidx)
  {
  RlxNMRerror(eidx,1);
  if(eidx) RlxNMRerror(0);
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// A           RELAXATION SUPEROPERATORS VIA ELEMENT CALCULATIONS
// ____________________________________________________________________________

// ------------------- Level 4 via Element Calculation ------------------------

void R_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12)

	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	// Return		void  : Level 4 relax. superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R4 
	//				          out	   in     12
	//
	// Note				T1s, T2s, & LOp assumed in proper bases


/*                  <a,a'| LOp |b,b'> += <a,a'| R4   |b,b'>
  		                                  1,2

   where a, a', b, b' are basis function indices.                            */

  {
  int hs = T1s[0].dim();			// Get Hilbert space size
  int aaa=0, bbb=0;
  complex Rel;
  for(int a=0; a<hs; a++)				// Sum over trans. a-aa
    for(int aa=0; aa<hs; aa++)				// with LS index aaa
      {
      bbb = 0;						// Transition LS index
      for(int b=0; b<hs; b++)				// Sum over trans. b-bb
        for(int bb=0; bb<hs; bb++)
          {
          Rel = LOp.get(aaa,bbb);			// Get LOp element
          Rel += R_4(hs,T1s,T2s,J12,rank,a,b,aa,bb); 	// Add the R component
          LOp.put(aaa, bbb, Rel);			// Put in new element
          bbb++;					// Next element index
          }
      aaa++;						// Next element row
      }
  }


double R_4(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
                                        int rank, int a, int b, int aa, int bb)

// sosi - why I switched indices is a mystery.  This appears correct but
//        it would probably be better to use a,aa,b,bb as the input order
//        because those are the transitions cross relaxing.  The a,b and
//        aa,bb correspond to tensor element indices but that is due to
//        the correlation of the transitions

	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//			rank  : Rank of the two interactions
	//			a, b  : 1st transition indices
	//			aa, bb: 2nd transition indices
	// Output		Rel   : Level 4 Relaxation superoperator element
	//
	//				    <a,aa|R4  |b,bb>
	//					    12
	//
	// Note				T1s, T2s, & LOp assumed in proper bases

/*                     --- [ ---
                       \   | \                 m       m 
   <a,a'|R4   |b,b'> = /   | /   delta     <a|T |g><b|T |g> J  (w  )
  	   1,2	       --- | ---      a',b'    1       2     12  gb
  		        m  [  g
  
                       m        m                       m        m
                 - <a|T |b><a'|T |b'> J  (w    ) - <b'|T |a'><b|T |a> J  (w  )
  	               1        2      12  b'a'         1        2     12  ab

                                      ---                                     ]
                                      \               m        m              |
                                   +  /   delta   <g|T |a'><g|T |b'> J  (b'g) |
                                      ---      a,b    1        2      12      |
                                       g                                      ]

 where m sums over angular momentum components, g over basis functions, and
 J values include xi1 and xi2                                                */

  {
  complex Rel = 0;
  complex J2 = J12.get(bb,aa);
  complex J3 = J12.get(a,b);
  int k=0, g=0;
  for(int m = -rank; m<=rank; m++)		// Sum over tensor components
     {
     Rel-=J2*T1s[k].get(a,b)*T2s[k].get(aa,bb);	// Add in terms RII
     Rel-=J3*T1s[k].get(bb,aa)*T2s[k].get(b,a);	// Add in terms RIII
     for(g=0; g<hs; g++)
       {
       if(aa == bb)				// Add terms RI over gamma sum
         Rel += J12.get(g,b)*
		T1s[k].get(a,g)*T2s[k].get(b,g);
       if(a == b) 				// Add terms RIV over gamma sum
         Rel += J12.get(bb,g)*
              T1s[k].get(g,aa)*T2s[k].get(g,bb);
       }
     k++;
     }
   return Re(Rel);
   }

// ------------------- Level 3 via Element Calculation ------------------------
//                     (Applies Secular Approximation)

void R_3(super_op& LOp, double* w, int rank, gen_op* T1s, gen_op* T2s,
	                                            matrix& J12, double cutoff)

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector of system energies
	//				in rad/sec in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//			cutoff: Secular approximation cutoff value (Hz)
	// Return		void  : Level 3 relax. superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R3
	//				          out	   in     12
	//
	// Note				T1s, T2s, & LOp assumed in proper bases


/*  <a,a'| LOp |b,b'> += <a,a'| R3   |b,b'> == del         * <a,a'| R4   |b,b'>
                                  1,2             w   ,w              1,2
                                                    aa'  bb'

   where a, a', b, b' are basis function indices.                            */

  {
  int hs = T1s[0].dim();			// Get Hilbert space size
  double delwaa, delwbb;			// Transition frequencies (lab)
  int aaa=0, bbb=0;				// Element indices in LS
  complex Rel;					// Temporary element storage
  for(int a=0; a<hs; a++)			// Sum over transition a-aa
    for(int aa=0; aa<hs; aa++)
      {
      delwaa = w[a] - w[aa];			// Transition a-aa frequency
      bbb = 0;					// Set column index to 0
      for (int b=0; b<hs; b++)			// Sum over transition b-bb
        for (int bb=0; bb<hs; bb++)		// (or R matrix columns)
          {
          delwbb = w[b] - w[bb];		// Transition b-bb frequency
	  if(fabs(delwaa-delwbb) < cutoff)	// Apply secular approximation
            {
            Rel = LOp.get(aaa,bbb);		// Get current R matrix element
            Rel += R_4(hs,T1s,T2s,J12,		// Add in the contribution here
	                        rank,a,b,aa,bb);
            LOp.put(aaa, bbb, Rel);		// Add to the R matrix element
            }
          bbb++;				// Next element (next column)
          }
      aaa++;
      }
  }

// ---------------- Level 2 via Element Calculation ---------------------
//                    (No Degenerate Transitions)

void R_2(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12)


	// Input		LOp   : Current relaxation superoperator
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
           if((a==b) && (aa==bb))		// Algorithm for diagonals
             {					//     <a,aa|R|a,aa>
             Rel = LOp.get(aaa,bbb);
             Rel += Rdiag_2(hs, T1s, T2s,
			      J12, rank, a, aa);
             LOp.put(aaa, bbb, Rel);
             }
           else if((a==aa) && (b==bb) && (a!=b))// Algorithm for off-diagonals
             {					//      <a,a|R|b,b>
             Rel = LOp.get(aaa,bbb);
             Rel += Rodiag_2(hs, T1s, T2s,
			      J12, rank, a, b);
             LOp.put(aaa, bbb, Rel);
             }
           bbb++;
           }
       aaa++;
       }
   return;
  }


//sosi switch this to complex
double R_2(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	         int rank, int a, int b, int aa, int bb)

// sosi - why I switched indices is a mystery.  This appears correct but
//        it would probably be better to use a,aa,b,bb as the input order
//        because those are the transitions cross relaxing.  The a,b and
//        aa,bb correspond to tensor element indices but that is due to
//        the correlation of the transitions (or the initial and final
//        eigenstates of the transitions)

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

//			      Level 2 Diagonal Elements
//
//                     --- [ ---
//                     \   | \       m       m 
// <a,a'|R2   |a,a'> = /   | /   <a|T |g><a|T |g>*J(w  )
//	   1,2	       --- | ---     1       2       ga
//		        m  [  g
//
//                                                 m        m
//                                     - 2.0 * <a|T |a><a'|T |a'> * J(w  )
//	     	                                   1        2          aa
//
//                                                   ---                            ]
//                                                   \       m        m             |
//                                                +  /   <g|T |a'><g|T |a'>*J(w   ) |
//	     	                                     ---     1        2        ga'  |
//										    ]
//
//			Level 2 Non-Zero Off Diagonal Elements
//
//                   ---
//                   \              m       m 
// <a,a|R2   |b,b> = /   -2.0 * <a|T |b><a|T |b> * J  (w  )
//	  1,2	     ---            1       2       12  ba

// where m sums over angular momentum components, g over basis functions,
// and J values include xi1 and xi2

   {
   complex Rel = 0;
   int k=0, g=0, m=0;
   complex J0 = 2.0;
   complex Jba = 2.0;
   if((a==b) && (aa==bb))			// Algorithm for diagonals
     {
     J0 *= J12.get(a,a);
     for(m=-rank, k=0; m<=rank; m++, k++)	// Sum over all tensor components
       {
       Rel -= J0*T1s[k].get(a,a)*T2s[k].get(aa,aa);
       for(g=0; g<hs; g++)
         {
         Rel += J12.get(g,a)*
		    T1s[k].get(a,g)*T2s[k].get(a,g);
         Rel += J12.get(aa,g)*
                  T1s[k].get(g,aa)*T2s[k].get(g,aa);
         }
       }
     }
   else if((a==aa) && (b==bb) && (a!=b))	// Algorithm for non-zero off-diagonals
     {
     Jba *= J12.get(b,a);
     for(m=-rank, k=0; m<=rank; m++, k++)	// Sum over all tensor components
       Rel -= Jba*T1s[k].get(a,b)*T2s[k].get(a,b);
     }
   return Re(Rel);
   }


//sosi switch this to complex
double Rodiag_2(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	                      int rank, int a, int b)

	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//			rank  : Rank of the 2 interactions
	//			a, b  : Transisiton indices
	// Output		Rel   : Level 2 off-diagonal relaxation superop
	//				element for interactions 1 & 2
	//
	//				    <a,a|R2  |b,b>
	//					   12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases

//			Level 2 Non-Zero Off Diagonal Elements
//
//                             ---
//                             \              m       m 
//           <a,a|R2   |b,b> = /   -2.0 * <a|T |b><a|T |b> * J  (w  )
//	            1,2	       ---            1       2       12  ba
//		                m

   {
   complex Rel = 0;
   complex Jba = 2.0*J12.get(b,a);
   int k=0;
   for(int m=-rank; m<=rank; m++, k++)		// Sum over all tensor components
     Rel -= Jba*T1s[k].get(a,b)*T2s[k].get(a,b);
   return Re(Rel);
   hs = 0;					// Compiler likes this used
   }


// sosi switch this to complex
double Rdiag_2(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	                      int rank, int a, int aa)

	// Input		hs    : Spin system Hilbert space
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	//			rank  : Rank of the 2 interactions
	//			a, aa : Transisiton indices
	// Output		Rel   : Level 2 diagonal relaxation superop element
	//				for interactions 1 & 2
	//
	//				    <a,aa|R2  |a,aa>
	//					    12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases

//			      Level 2 Diagonal Elements
//
//                     --- [ ---
//                     \   | \       m       m 
// <a,a'|R2   |a,a'> = /   | /   <a|T |g><a|T |g>*J(w  )
//	   1,2	       --- | ---     1       2       ga
//		        m  [  g
//
//                                                 m        m
//                                     - 2.0 * <a|T |a><a'|T |a'> * J(w  )
//	     	                                   1        2          aa
//
//                                                   ---                            ]
//                                                   \       m        m             |
//                                                +  /   <g|T |a'><g|T |a'>*J(w   ) |
//	     	                                     ---     1        2        ga'  |
//                                                    g                             ]

// where m sums over angular momentum components, g over basis functions,
// and J values include xi1 and xi2

   {
   complex Rel = 0.0;
   complex J0 = 2.0*J12.get(a,a);
   int k=0, g=0;
   for (int m = -rank; m<=rank; m++)		// Sum over all tensor components
     {
     Rel -= J0*T1s[k].get(a,a)*T2s[k].get(aa,aa);
     for(g=0; g<hs; g++)
       {
       Rel += J12.get(g,a)*
		T1s[k].get(a,g)*T2s[k].get(a,g);
       Rel += J12.get(aa,g)*
              T1s[k].get(g,aa)*T2s[k].get(g,aa);
       }
     k++;
     }
   return Re(Rel);
   }


// ---------------- Level 0 via Element Calculation ---------------------
//			  (Extreme Narrowing)

void R_0(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, const complex& J12)

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

   {
   int hs = T1s[0].dim();			// Get Hilbert space size
   int aaa=0, bbb=0;
   complex Rel;
   for(int a=0; a<hs; a++)			// Sum over transition a-aa
     for(int aa=0; aa<hs; aa++)
       {
       bbb = 0;
       for(int b=0; b<hs; b++)			// Sum over transition b-bb
         for(int bb=0; bb<hs; bb++)
           {
           Rel = LOp.get(aaa,bbb);
           Rel += J12*R_0(hs, T1s, T2s,
	                   rank, a, b, aa, bb);
           LOp.put(aaa, bbb, Rel);		// Add to the relaxation matrix element
           bbb++;
           }
       aaa++;
       }
   return;
   }

// sosi switch this to complex
double R_0(int hs, gen_op* T1s, gen_op* T2s,
		 	         int rank, int a, int b, int aa, int bb)

// sosi - why I switched indices is a mystery.  This appears correct but
//        it would probably be better to use a,aa,b,bb as the input order
//        because those are the transitions cross relaxing.  The a,b and
//        aa,bb correspond to tensor element indices but that is due to
//        the correlation of the transitions (or the initial and final
//        eigenstates of the transitions)

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


//                     --- [ ---
//                     \   | \                 m       m 
// <a,a'|R0   |b,b'> = /   | /   delta     <a|T |g><b|T |g>
//	   1,2	       --- | ---      a',b'    1       2
//		        m  [  g
//
//                               m        m                     m 
//                         - <a|T |b><a'|T |b'> - <b'|T |a'><b|T |a>
//	     	                 1        2            1        2
//
//                                                   ---                            ]
//                                                   \               m        m     |
//                                                +  /   delta   <g|T |a'><g|T |b'> |
//	     	                                     ---      a,b    1        2     |
//                                                    g                             ]

   {
   complex Rel = 0;
   int k=0, g=0;
   for (int m = -rank; m<=rank; m++)		// Sum over all tensor components
     {
     Rel -= T1s[k].get(a,b)*T2s[k].get(aa,bb);	// Add in terms RII
     Rel -= T1s[k].get(bb,aa)*T2s[k].get(b,a);	// Add in terms RIII
     for(g=0; g<hs; g++)
       {
       if(aa == bb)				// Add in terms RI over gamma sum
         Rel += T1s[k].get(a,g)*T2s[k].get(b,g);
       if(a == b) 				// Add in terms RIV over gamma sum
         Rel += T1s[k].get(g,aa)*T2s[k].get(g,bb);
       }
     k++;
     }
   return Re(Rel);
   }


// ____________________________________________________________________________
// B            RELAXATION SUPEROPERATORS VIA DOUBLE COMMUTATORS
// ____________________________________________________________________________


void R_4s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12)

	// Input		LOp   : Current relaxation superoperator
	//			rank  : Rank of the interactions involved
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	// Return		void  : Level 4 relax. superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R4
	//				          out	   in     12
	//
	// Note				T1s, T2s, & LOp assumed in proper bases

//                   --- ---
//                   \   \       m    m     -m
//	     R4    = /   /   (-1)  [ T , [ T   ,  ] ] J  (w )
//	       1,2   --- ---          1     2,p        12  p
//		       p   m

// where p sums over transitions, m sums over angular momentum components,
// and J values include xi1 and xi2

   {
   int hs = T1s[0].dim();			// Get Hilbert space size
   matrix mx(hs, hs, 0.0);			// Null Hilbert space matrix
   basis bs = T1s[0].get_basis();		// Current basis of spin tensors
   gen_op nullOp;				// Set up a null operator
   int k=0, m=0;				// Angular momentum indices
   double Jwp;
   complex Taaa;
   gen_op *T2ps;				// Components for specific transitions
   T2ps = new gen_op[2*rank+1];			//     (2*l+1 of these for each p)
   for(int a=0; a<hs; a++)			// Sum over transitions a-aa: p
     for(int aa=0; aa<hs; aa++)
       {
       for(m=-rank, k=0; m<=rank; m++, k++)	// Perform a pre-summation over m
         {					// and set up the spin tensor operator
         T2ps[k] = nullOp;			// components for this transition
         Taaa = T2s[k].get(a,aa);		//                       m
         if(Re(Taaa) || Im(Taaa))		//           T2ps[k] == T
           {					//                       2,p
           mx.put(Taaa, a, aa);			// 
           T2ps[k] = gen_op(mx,bs);		// These are either sparse or null arrays
           mx.put(0.0, a, aa);			// so we insure that null components are
           }					// recognized (and not used in computation).
         }
       Jwp = Re(J12.get(a,aa));			// Get spectral density for this transition
       R_CC_0(T1s, T2ps, LOp, rank, Jwp);	// Now, perform sum over m, and add to LOp
       }
   return;
   }


// ------------------ Level 3 via Double Commutators --------------------
//                    (Applies Secular Approximation)

void R_3s(super_op& LOp, double* w, int rank, gen_op* T1s, gen_op* T2s, matrix& J12)

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector of system energies
	//				in rad/sec in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	// Return		void  : Level 3 relaxation superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R3
	//				          out	   in     12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases

//            --- --- ---
//            \   \   \       m                 m       -m
//    R3    = /   /   /   (-1)  delta        [ T   , [ T    ,  ] ] J  (w  )
//	1,2   --- --- ---            w , -w     1,p     2,p'        12  p'
//	       p   p'  m              p    p'

// where p and p' sum over transitions, m sums over angular momentum components,
// and J values include xi1 and xi2.  The secular approximation restricts the
// sum to contain only terms where the transition frequencies for p & p' cancel.

   {
   double cutoff = 1.e-9;			// Secular approximation cutoff
   int hs = T1s[0].dim();			// Get Hilbert space size
   matrix mx1(hs, hs, 0.0);			// Two Zeroed Hilbert space matrices
   matrix mx2(hs, hs, 0.0);
   basis bs = T1s[0].get_basis();		// Current working basis
   gen_op nullOp;				// A NULL operator
   double wp, wpp;				// Transition frequencies for p & p'
   double Jwp;					// Spectral density, transition p (& p')
   complex Ta_aa, Tb_bb;
   gen_op *T1ps, *T2ps;				// Components of the two spin tensors
   T1ps = new gen_op[2*rank+1]; 		// (2*l+1 of these for each p)
   T2ps = new gen_op[2*rank+1];			// (2*l+1 of these for each p')
   int m=0, k=0;				// Angular momentum indicies
   for(int a=0; a<hs; a++)			// Sum over transitions a-aa: p
     for(int aa=0; aa<hs; aa++)
       {
       wp = w[a] - w[aa];			// Transition frequency for a-aa: p
       Jwp = Re(J12.get(a,aa));			// Get spectral density for transition p
       for(m=-rank, k=0; m<=rank; m++, k++)	// Perform a pre-summation over m
         {					// and set up the spin tensor operator
         T1ps[k] = nullOp;			// components for this transition
         Ta_aa = T1s[k].get(a,aa);		//                       m
         if(Re(Ta_aa) || Im(Ta_aa))		//           T1ps[k] == T
           {					//                       1,p
           mx1.put(Ta_aa, a, aa);		// 
           T1ps[k] = gen_op(mx1,bs);		// These are either sparse or null arrays
           mx1.put(0.0, a, aa);			// so we insure that null components are
           }					// recognized (and not used in computation).
         }
       for(int b=0; b<hs; b++)			// Sum over transitions b-bb: p'
         for(int bb=0; bb<hs; bb++)
           {
           wpp = w[b] - w[bb];			// Transition b-bb: p' frequency
	   if(fabs(wp+wpp) < cutoff)		// Apply the secular approximation: p = -p'
             {
             for(m=-rank,k=0; m<=rank; m++,k++)	// Perform a pre-summation over m
               {				// and set up the spin tensor operator
               T2ps[k] = nullOp;		// components for this transition
               Tb_bb = T2s[k].get(b,bb);	//                       m
               if(Re(Tb_bb) || Im(Tb_bb))	//           T2ps[k] == T
                 {				//                       2,p'
                 mx2.put(Tb_bb, b, bb);		//
                 T2ps[k] = gen_op(mx2,bs);	// These are either sparse or null arrays
                 mx2.put(0.0, b, bb);		// so we insure that the null components are
                 }				// recognized (and not used in computaiton).
               }
             R_CC_0(T1ps, T2ps, LOp, rank, Jwp);// Perform sum over m, add LOp    to LOp
             }					//			      pp'
           }					// Next transition p'
       }					// Next transition p
   return;
   }

// --------------------- Level 2 via Double Commutators -----------------------
//                         (No Degenerate Transitions)

void R_2s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12)

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector of system energies
	//				in rad/sec in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Matrix of reduced spectral density 
	//				function values for 2 interactions
	//				times the two interaction constants
	// Return		void  : Level 2 relax. superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R2
	//				          out	   in     12
	//
	// Note				T1s, T2s, & LOp assumed in proper bases

//                         --- ---
//                         \   \       m    m       -m
//                 R2    = /   /   (-1)  [ T   , [ T   ,  ] ] J  (w )
//	             1,2   --- ---          1,p     2,p        12  p
//		            p   m

// where p sums over transitions, m sums over angular momentum components,
// and J values include xi1 and xi2.  The Level 2 treatment assumes that there
// are no degenerate transitions in the system.

   {
   int hs = T1s[0].dim();			// Get Hilbert space size
   matrix mx1(hs, hs, 0.0);			// Two Zeroed Hilbert space matrices
   matrix mx2(hs, hs, 0.0);
   basis bs = T1s[0].get_basis();		// Current working basis
   gen_op nullOp;				// A NULL operator
   double Jwp;					// Spectral density, transition p
   complex Ta_aa;
   gen_op *T1ps, *T2ps;				// Components of the two spin tensors
   T1ps = new gen_op[2*rank+1]; 		// (2*l+1 of these for each p)
   T2ps = new gen_op[2*rank+1];			// (2*l+1 of these for each p)
   int m=0, k=0;				// Angular momentum indicies
   for(int a=0; a<hs; a++)			// Sum over transitions a-aa: p
     for(int aa=0; aa<hs; aa++)
       {
       for(m=-rank, k=0; m<=rank; m++, k++)	// Perform a pre-summation over m
         {					// and set up the spin tensor operator
         T1ps[k] = nullOp;			// components for this transition
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
           T2ps[k] = gen_op(mx2,bs);		// These are either sparse or null arrays
           mx2.put(0.0, aa, a);			// so we insure that the null components are
           }					// recognized (and not used in computaiton).
         }
       Jwp = Re(J12.get(a,aa));			// Get spectral density for this transition
       R_CC_0_trans(T1ps, T2ps, LOp, rank, Jwp);// Now, perform sum over m, and add to LOp
       }
   return;
   }

// --------------------- Level 1 via Double Commutators -----------------------
//              (Spin Tensor Component Constant Oscillation)

//void R_1s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12)

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

//                         ---
//                         \       m    m     -m
//                 R1    = /   (-1)  [ T , [ T  ,  ] ] J  (mw )
//	             1,2   ---          1     2         12   0
//		            m

// m sums over angular momentum components, and J values include xi1 and xi2.
// The Level 1 treatment assumes that there are no degenerate transitions in the system.


// ------------------ Level 0 via Double Commutators --------------------
//                         (Extreme Narrowing)

//void R_0s(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, double J12)

	// Input		LOp   : Current relaxation superoperator
	//			w     : Vector of system energies
	//				in rad/sec in the LAB frame
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			J12   : Reduced spectral density function value
	//				at 0 Hz times the 2 interact. constants
	// Return		void  : Level 0 relax. superop for interactions
	//			        1 & 2 added to input superoperator LOp
	//
	//				       LOp    = LOp   + R0
	//				          out	   in     12
	//
	// Note				T1s, T2s, and LOp assumed in proper bases

//                               ---
//                               \       m    m     -m
//                 R01,  = J (0) /   (-1)  [ T , [ T  ,  ] ]
//	             1,2    12   ---          1     2
//		                  m

// m sums over angular momentum components, and J value include xi1 and xi2.
// The Level 1 treatment assumes extreme narrowing.


// ____________________________________________________________________________
//                      USEFUL DOUBLE COMMUTATOR SUPEROPERATORS
// ____________________________________________________________________________

// -------------------------- Level 0 Functions -------------------------------

super_op R_AC_0(spin_T &T)

	// Input		T     : Spin tensor
	// Output		LOp   : Prototype autocorrelation relaxation
	//				superoperator

//                   --- [  ]m [   m  [   -m   ] ]
//             LOp = \   |-1|  | T1 , | T1  ,  | |
//                   /   [  ]  [   l  [   l    ] ]
//	             ---
//                    m

   {
   super_op LOp;
   gen_op Op1, Op2;
   int l = T.Rank();
   Op1 = T.component(l,0);		// m = 0 component
   LOp = d_commutator(Op1);
   Op1 = T.component(l,1);
   Op2 = T.component(l,-1);
   LOp -= d_commutator(Op1, Op2);	// m=1 with m=-1
   LOp -= d_commutator(Op2, Op1);	// m=-1 with m=1
   if(l > 1)
     {
     Op1 = T.component(l,2);
     Op2 = T.component(l,-2);
     LOp += d_commutator(Op1, Op2);	// m=2 with m=-2
     LOp += d_commutator(Op2, Op1);	// m=-2 with m=2
     }
   return LOp;
   }




	// Input		T     : Spin tensor
	// 			LOp1  : Super operator
	//			Op    : General operator
	// 			xisq  : Scaling factor
	// Output		None  : LOp has the prototype autocorrelation
	//				relaxation superoperator added to it.
	// Note			      : Computed in the basis of Op.

//                    --- [  ]m [   m  [   -m   ] ]
//             LOp += \   |-1|  | T1 , | T1  ,  | |
//                    /   [  ]  [   l  [   l    ] ]
//	              ---
//                     m

void R_AC_0(spin_T &T, super_op &LOp1, gen_op &Op, double xisq)
   {
   super_op LOp;
   gen_op Op1, Op2;
   int l = T.Rank();
   Op1 = T.component(l,0);		// m = 0 component
   Op1.Op_base(Op);
   LOp += d_commutator(Op1);
   Op1 = T.component(l,1);
   Op1.Op_base(Op);
   Op2 = T.component(l,-1);
   Op2.Op_base(Op);
   LOp -= d_commutator(Op1, Op2);	// m=1 with m=-1
   LOp -= d_commutator(Op2, Op1);	// m=-1 with m=1
   if(l > 1)
     {
     Op1 = T.component(l,2);
     Op1.Op_base(Op);
     Op2 = T.component(l,-2);
     Op2.Op_base(Op);
     LOp += d_commutator(Op1, Op2);	// m=2 with m=-2
     LOp += d_commutator(Op2, Op1);	// m=-2 with m=2
     }
   if(xisq != 1)
     LOp *= xisq;
   LOp1 += LOp;
   return;
   }


void R_AC_0(gen_op *Ts, super_op &LOp1, int rank, double xisq)

	// Input		Ts    : Spin tensor components
	// 			LOp1  : Super operator
	// 			rank  : Interaction rank
	// 			xisq  : Scaling factor
	// Output		None  : LOp has prototype autocorrelation
	//				relaxation superoperator added to it.
	// Note			      : Computed in the basis of Op.

//                    --- [  ]m [   m  [   -m   ] ]
//             LOp += \   |-1|  | T1 , | T1  ,  | |
//                    /   [  ]  [   l  [   l    ] ]
//	              ---
//                     m

   {
   super_op LOp;
   if(rank == 1)
     {
     LOp = d_commutator(Ts[1]);			// m=0 with m=0
     LOp -= d_commutator(Ts[2], Ts[0]);		// m=1 with m=-1
     LOp -= d_commutator(Ts[0], Ts[2]);		// m=-1 with m=1
     }
   else if(rank == 2)
     {
     LOp = d_commutator(Ts[2]);			// m=0 with m=0
     LOp -= d_commutator(Ts[3], Ts[1]);		// m=1 with m=-1
     LOp -= d_commutator(Ts[1], Ts[3]);		// m=-1 with m=1
     LOp += d_commutator(Ts[4], Ts[0]);		// m=2 with m=-2
     LOp += d_commutator(Ts[0], Ts[4]);		// m=-2 with m=2
     }
   if(xisq != 1)
     LOp *= xisq;
   LOp1 += LOp;
   return;
  }


void R_CC_0(spin_T &T1, spin_T &T2, super_op &LOp1,
		                     gen_op &Op, double xisq)

	// Input		T1    : Spin tensor 1
	// 			T2    : Spin tensor 2
	// 			LOp1  : Super operator
	//			Op    : General operator
	// 			xisq  : Scaling factor
	// Output		None  : LOp has prototype cross-correlation
	//				relaxation superoperator added to it.
	// Note			      : Computed in the basis of Op.

//                    --- [  ]m [   m  [   -m   ] ]
//             LOp += \   |-1|  | T1 , | T2  ,  | |
//                    /   [  ]  [   l  [   l    ] ]
//	              ---
//                     m

   {
   gen_op Op1, Op2;
   super_op LOp;
   int l = T1.Rank();
   Op1 = T1.component(l,0);
   Op1.Op_base(Op);
   Op2 = T2.component(l,0);
   Op2.Op_base(Op);
   LOp = d_commutator(Op1, Op2); 		// m=0 term
   Op1 = T1.component(l,1);
   Op1.Op_base(Op);
   Op2 = T2.component(l,-1);
   Op2.Op_base(Op);
   LOp = d_commutator(Op1, Op2);		// m = (+/-)1 terms
   LOp -= d_commutator(Op1, Op2);		// m = (+/-)1 terms
   LOp -= d_commutator(Op2, Op1);
   if(l > 1)
     {
     Op1 = T1.component(l,2);
     Op1.Op_base(Op);
     Op2 = T2.component(l,-2);
     Op2.Op_base(Op);
     LOp += d_commutator(Op1, Op2);		// m = (+/-)2 terms
     LOp += d_commutator(Op2, Op1);
     }
   if(xisq != 1)
     LOp *= xisq;
   LOp1 += LOp;
   return;
  }


void R_CC_0(gen_op *T1s, gen_op* T2s, super_op &LOp1, int rank, double xisq)

	// Input		T1s   : Spin tensor 1 components
	// 			T2    : Spin tensor 2 components
	// 			LOp1  : Super operator
	// 			rank  : Interaction rank
	// 			xisq  : Scaling factor
	// Output		None  : LOp has prototype cross-corrleation
	//			        relaxation superoperator added to it.
	// Note		              : Computed in the basis of Op.

//                    --- [  ]m [   m  [   -m   ] ]
//             LOp += \   |-1|  | T1 , | T2  ,  | |
//                    /   [  ]  [   l  [   l    ] ]
//	              ---
//                     m

   {
   super_op LOp;
   if(rank == 1)
     {
     LOp = d_commutator(T1s[1], T2s[1]);	// m=0 with m=0
     LOp -= d_commutator(T1s[2], T2s[0]);	// m=1 with m=-1
     LOp -= d_commutator(T1s[0], T2s[2]);	// m=-1 with m=1
     }
   else if(rank == 2)
     {
     LOp = d_commutator(T1s[2], T2s[2]);	// m=0 with m=0
     LOp -= d_commutator(T1s[3], T2s[1]);	// m=1 with m=-1
     LOp -= d_commutator(T1s[1], T2s[3]);	// m=-1 with m=1
     LOp += d_commutator(T1s[4], T2s[0]);	// m=2 with m=-2
     LOp += d_commutator(T1s[0], T2s[4]);	// m=-2 with m=2
     }
   if(xisq != 1)
     LOp *= xisq;
   LOp1 += LOp;
   return;
  }


void R_CC_0_trans(gen_op *T1s, gen_op* T2s, super_op &LOp1, int rank, double xisq)

	// Input		T1s   : Spin tensor 1 components
	// 			T2    : Spin tensor 2 components
	// 			LOp1  : Super operator
	// 			rank  : Interaction rank
	// 			xisq  : Scaling factor
	// Output		None  : LOp has prototype cross-correlation
	//				relaxation superoperator added to it.
	// Note			      : Computed in the basis of Op.
	// Note			None  : It is assumed the vector T2s
	//				contains the transposed operators

//                    --- [   m  [ {  m}t   ] ]
//             LOp += \   | T1 , | |T2 | ,  | |
//                    /   [   l  [ {  l}    ] ]
//	              ---
//                     m

   {
   super_op LOp;
   if(rank == 1)
     {
     LOp = d_commutator(T1s[0], T2s[0]);	// m=-1 with m=-1
     LOp += d_commutator(T1s[1], T2s[1]);	// m=0 with m=0
     LOp += d_commutator(T1s[2], T2s[2]);	// m=1 with m=1
     }
   else if(rank == 2)
     {
     LOp = d_commutator(T1s[0], T2s[0]);	// m=-2 with m=-2
     LOp += d_commutator(T1s[1], T2s[1]);	// m=-1 with m=-1
     LOp += d_commutator(T1s[2], T2s[2]);	// m=0 with m=0
     LOp += d_commutator(T1s[3], T2s[3]);	// m=1 with m=1
     LOp += d_commutator(T1s[4], T2s[4]);	// m=2 with m=2
     }
   if(xisq != 1)
     LOp *= xisq;
   LOp1 += LOp;
   return;
  }

// ----------------------- Level 1 Functions ----------------------------


void R_AC_1(spin_T &T, super_op &LOp1, gen_op &Op,
		             	double J0, double J1, double J2)

	// Input		T     : Spin tensor
	// 			LOp1  : Super operator
	//			Op    : General operator
	// 			J0    : Spectral density at 0 Hz
	// 			J1    : Spectral density at Larmor frequency
	// 			J2    : Spectral density at 2*Larmor
	// Output		None  : LOp has prototype relaxation superop
	//				added to it.
	// Note			      : Computed in the basis of Op.


//                    ---        [  ]m [   m  [   -m   ] ]
//             LOp += \   J(mw ) |-1|  | T1 , | T1  ,  | |
//                    /       0  [  ]  [   l  [   l    ] ]
//	              ---
//                     m

   {
   super_op LOp;
   gen_op Op1, Op2;
   int l = T.Rank();
   Op1 = T.component(l,0);		// m = 0 component
   Op1.Op_base(Op);
   LOp = J0*d_commutator(Op1);
   Op1 = T.component(l,1);
   Op1.Op_base(Op);
   Op2 = T.component(l,-1);
   Op2.Op_base(Op);
   LOp -= d_commutator(Op1, Op2, J1);	// m=1 with m=-1
   LOp -= d_commutator(Op2, Op1, J1);	// m=-1 with m=1
   if(l > 1)
     {
     Op1 = T.component(l,2);
     Op1.Op_base(Op);
     Op2 = T.component(l,-2);
     Op2.Op_base(Op);
     LOp += d_commutator(Op1, Op2, J2);	// m=2 with m=-2
     LOp += d_commutator(Op2, Op1, J2);	// m=-2 with m=2
     }
   LOp1 += LOp;
   return;
   }


void R_AC_1(gen_op *Ts, super_op &LOp1, int rank,
		             	double J0, double J1, double J2)

	// Input		Ts   : Spin tensor 1 components
	// 			LOp1  : Super operator
	// 			rank  : Interaction rank
	// 			J0    : Spectral density at 0 Hz
	// 			J1    : Spectral density at Larmor frequency
	// 			J2    : Spectral density at 2*Larmor
	// Output		None  : LOp has prototype auto-correlation
	//				relaxation superoperator added to it.
	// Note				All operators (T1s & T2s) should be
	//				in the perfered computation basis

//                    ---        [  ]m [   m  [   -m   ] ]
//             LOp += \   J(mw ) |-1|  | T1 , | T1  ,  | |
//                    /       0  [  ]  [   l  [   l    ] ]
//	              ---
//                     m

   {
   super_op LOp;
   if(rank == 1)
     {
     LOp = d_commutator(Ts[1], J0);		// m=0 with m=0
     LOp -= d_commutator(Ts[2], Ts[0], J1);	// m=1 with m=-1
     LOp -= d_commutator(Ts[0], Ts[2], J1);	// m=-1 with m=1
     }
   else if(rank == 2)
     {
     LOp = d_commutator(Ts[2], J0);		// m=0 with m=0
     LOp -= d_commutator(Ts[3], Ts[1], J1);	// m=1 with m=-1
     LOp -= d_commutator(Ts[1], Ts[3], J1);	// m=-1 with m=1
     LOp += d_commutator(Ts[4], Ts[0], J2);	// m=2 with m=-2
     LOp += d_commutator(Ts[0], Ts[4], J2);	// m=-2 with m=2
     }
   LOp1 += LOp;
   return;
  }


void R_CC_1(spin_T &T1, spin_T &T2, super_op &LOp1,
		        gen_op &Op, double J0, double J1, double J2)

	// Input		T1    : Spin tensor 1
	// 			T2    : Spin tensor 2
	// 			LOp1  : Super operator
	//			Op    : General operator
	// 			J0    : Spectral density at 0 Hz
	// 			J1    : Spectral density at Larmor frequency
	// 			J2    : Spectral density at 2*Larmor
	// Output		None  : LOp has cross-correlation relaxation
	//			        superoperator added to it.
	// Note			      : Computed in the basis of Op.

//                    ---        [  ]m [   m  [   -m   ] ]
//             LOp += \   J(mw ) |-1|  | T1 , | T2  ,  | |
//                    /       0  [  ]  [   l  [   l    ] ]
//	              ---
//                     m

   {
   gen_op Op1, Op2;
   super_op LOp;
   int l = T1.Rank();
   Op1 = T1.component(l,0);
   Op1.Op_base(Op);
   Op2 = T2.component(l,0);
   Op2.Op_base(Op);
   LOp = d_commutator(Op1, Op2, J0);
   Op1 = T1.component(l,1);
   Op1.Op_base(Op);
   Op2 = T2.component(l,-1);
   Op2.Op_base(Op);
   LOp -= d_commutator(Op1, Op2, J1);
   LOp -= d_commutator(Op2, Op1, J1);
   if(l > 1)
     {
     Op1 = T1.component(l,2);
     Op1.Op_base(Op);
     Op2 = T2.component(l,-2);
     Op2.Op_base(Op);
     LOp += d_commutator(Op1, Op2, J2);
     LOp += d_commutator(Op2, Op1, J2);
     }
   LOp1 += LOp;
   return;
   }


void R_CC_1(gen_op *T1s, gen_op* T2s, super_op &LOp1,
		        int rank, double J0, double J1, double J2)

	// Input		T1s   : Spin tensor 1 components
	// 			T2    : Spin tensor 2 components
	// 			LOp1  : Super operator
	// 			rank  : Interaction rank
	// 			J0    : Spectral density at 0 Hz
	// 			J1    : Spectral density at Larmor frequency
	// 			J2    : Spectral density at 2*Larmor
	// Output		None  : LOp has prototype cross-correlation
	//				relaxation superoperator added to it.
	// Note				All operators (T1s & T2s) should be
	//				in the perfered computation basis

//                    ---        [  ]m [   m  [   -m   ] ]
//             LOp += \   J(mw ) |-1|  | T1 , | T2  ,  | |
//                    /       0  [  ]  [   l  [   l    ] ]
//	              ---
//                     m

   {
   super_op LOp;
   if(rank == 1)
     {
     LOp = d_commutator(T1s[1], T2s[1], J0);	// m=0 with m=0
     LOp -= d_commutator(T1s[2], T2s[0], J1);	// m=1 with m=-1
     LOp -= d_commutator(T1s[0], T2s[2], J2);	// m=-1 with m=1
// sosi - are the above J's  & T's the proper combination?
// here, m = i-rank, so 0 --> 1-1; 1-->2-1; -1-->0-1
// I don't think J2 is correct, should be J1?
     }
   else if(rank == 2)
     {
     LOp = d_commutator(T1s[2], T2s[2], J0);	// m=0 with m=0
     LOp -= d_commutator(T1s[3], T2s[1], J1);	// m=1 with m=-1
     LOp -= d_commutator(T1s[1], T2s[3], J1);	// m=-1 with m=1
     LOp += d_commutator(T1s[4], T2s[0], J2);	// m=2 with m=-2
     LOp += d_commutator(T1s[0], T2s[4], J2);	// m=-2 with m=2
     }
   LOp1 += LOp;
   return;
  }

// ____________________________________________________________________________
//                RELAXATION SUPEROPERATOR GENERAL FUNCTIONS
// ____________________________________________________________________________

// ------------------- Routing to Specific Level Computation ------------------
 
void Rmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
            double* taus, double* c1s, double* c2s, double xi1xi2,
         double w0, double w1, double w2, int l, int level, int autoc, int het)

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
	//			autoc : Flag for auto-corr vs cross-correlation
	// Output		LOp   : Relaxation superoperator
	// Note			      :	Computed in the basis T1s and T2s

/*                            LOp += R   (Level)
                                      1,2                                    */

  {
  matrix J12;					// Used for J's levels > 1
  double J0, J1, J2;				// Used for J's levels 0 & 1
  switch(level)
    {
    case 4:					// Lev 4 mu1-mu2: elem. by elem.
       J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spect. dens.
       J12 *= complex(xi1xi2);			// Scale Js by interact. const.
       if(!het) R_4(LOp, l, T1s, T2s, J12);	// No secular approx if homonuc
       else R_3(LOp,w,l,T1s,T2s, J12, 1.e6);	// Use secular approx if hetero 
       break;
     case -4:					// Lev4 mu1-mu2: double comm.
       J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spect. dens.
       J12 *= complex(xi1xi2);			// Scale Js by interact. const.
       R_4s(LOp, l, T1s, T2s, J12);		// Generate Rmumu contribution
       break;
     case 3:					// Lev 3 mu1-mu2: elem. by elem.
       J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spect.dens.
       J12 *= complex(xi1xi2);			// Scale Js by interact. const.
       R_3(LOp, w, l, T1s, T2s, J12);		// Generate Rmumu elementwise
       break;
     case -3:					// Lev 3 mu1-mu2: double comm.
       J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spect. dens.
       J12 *= complex(xi1xi2);			// Scale Js by interact. const.
       R_3s(LOp, w, l, T1s, T2s, J12);		// Generate Rmumu dbl. commut.
       break;
     case 2:					// Lev 2 mu1-mu2: elem. by elem.
       J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spect. dens.
       J12 *= complex(xi1xi2);			// Scale Js by interact. const.
       R_2(LOp, l, T1s, T2s, J12);		// Generate Rmumu elementwise
       break;
     case -2:					// Lev 2 mu1-mu2: double comm.
       J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spect. dens.
       J12 *= complex(xi1xi2);			// Scale Js by interact. const.
       R_2s(LOp, l, T1s, T2s, J12);		// Generate Rmumu dbl. commut.
       break;
     case 1: 					// Lev 1 mu1-mu2: double comm.
       J0 = J_reduced(taus, c1s, c1s, w0, 1);	// J(0)  - spect. density @ 0
       J1 = J_reduced(taus, c1s, c1s, w1, 1);	// J(w0) - spect. density @ W
       J2 = J_reduced(taus, c1s, c1s, w2, 1);	// J(2w0)- sepct. density @ 2W
       if(autoc)
         R_AC_1(T1s,LOp,l,J0*xi1xi2,J1*xi1xi2,J2*xi1xi2);
       else
         R_CC_1(T1s, T2s,LOp,l,J0*xi1xi2,J1*xi1xi2,J2*xi1xi2);
       break;
     case 0:					// Lev 0 mu1-mu2: elem. by elem.
       J0 = J_reduced(taus,c1s,c2s,0.0,1);	// Need J(0) only, ext. narrow.
       if(fabs(xi1xi2*J0) > 1.e-6)
         R_0(LOp,l,T1s,T2s,complex(xi1xi2*J0));
       break;
     default:					// Lev 0 mu1-mu2: double comm.
       J0 = J_reduced(taus,c1s,c2s,0.0,1);	// Need J(0) only, ext. narrow.
       if(fabs(xi1xi2*J0) > 1.e-6)		// Only calc. if it contributes
       { if(autoc) R_AC_0(T1s,    LOp,l,xi1xi2*J0);
         else      R_CC_0(T1s,T2s,LOp,l,xi1xi2*J0);
       }
       break;
     }
   return;
   }

 
void Rmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w,
                      int hs, double tau, double xi1xi2, double w0,
                             double w1, double w2, int l, int level, int autoc)

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
	//			autoc : Flag for auto-corr. vs cross-corr
	// Output		LOp   : Relaxation superoperator
	// Note			      :	Computed in the basis T1s and T2s
	// Note			      :	Uses spherical top J(w)'s

//                            LOp += R   (Level)
//                                    1,2  

// sosi: the rank is set to 2 which is inappropriate (in general)!
  {
  matrix J12;					// Used for J's levels > 1
  double J0, J1, J2;				// Used for J's levels 0 & 1
  const double pi4 = 12.566370614;		// Used to scale generic J to reduced J

//	      Start the Relaxation Superoperator Computation

  switch(level)
    {
     case 4:					// Lev 4 mu1-mu2: element by element
       J12 = J_gen(tau, w, hs, 1);	 	// Get all reduced spectral densities
       J12 *= complex(xi1xi2)/pi4;		// Scale by the interaction constants
       R_4(LOp, l, T1s, T2s, J12);
       break;
     case -4:					// Lev 4 mu1-mu2: double commutator
       J12 = J_gen(tau, w, hs, 1);	 	// Get all reduced spectral densities
       J12 *= complex(xi1xi2)/pi4;		// Scale by the interaction constants
       R_4s(LOp, l, T1s, T2s, J12);
       break;
     case 3:					// Lev 3 mu1-mu2: element by element
       J12 = J_gen(tau, w, hs, 1);	 	// Get all reduced spectral densities
       J12 *= complex(xi1xi2)/pi4;		// Scale by the interaction constants
       R_3(LOp, w, l, T1s, T2s, J12);
       break;
     case -3:					// Lev 3 mu1-mu2: double commutator
       J12 = J_gen(tau, w, hs, 1);	 	// Get all reduced spectral densities
       J12 *= complex(xi1xi2)/pi4;		// Scale by the interaction constants
       R_3s(LOp, w, l, T1s, T2s, J12);
       break;
     case 2:					// Lev 2 mu1-mu2: element by element
       J12 = J_gen(tau, w, hs, 1);	 	// Get all reduced spectral densities
       J12 *= complex(xi1xi2)/pi4;		// Scale by the interaction constants
       R_2(LOp, l, T1s, T2s, J12);
       break;
     case -2:					// Lev 2 mu1-mu2: double commutator
       J12 = J_gen(tau, w, hs, 1);	 	// Get all reduced spectral densities
       J12 *= complex(xi1xi2)/pi4;		// Scale by the interaction constants
       R_2s(LOp, l, T1s, T2s, J12);
       break;
     case 1: 					// Lev 1 mu1-mu2: double commutator
       J0 = J_gen(tau, w0, 1)/pi4;		// J(0)
       J1 = J_gen(tau, w1, 1)/pi4;		// J(w0)
       J2 = J_gen(tau, w2, 1)/pi4;		// J(2w0)
       if(autoc)
         R_AC_1(T1s, LOp, l,
                     J0*xi1xi2,J1*xi1xi2,J2*xi1xi2);
       else
         R_CC_1(T1s, T2s, LOp, l,
                   J0*xi1xi2,J1*xi1xi2,J2*xi1xi2);
       break;
     case 0:					// Lev 0 mu1-mu2: element by element
       J0 = J_gen(tau, 0.0, 1)/pi4;		// Need J(0) only, extreme narrowing
       if(fabs(xi1xi2*J0) > 1.e-6)
         R_0(LOp,l,T1s,T2s,complex(xi1xi2*J0));
       break;
     default:					// Lev 0 mu1-mu2: double commutator
       J0 = J_gen(tau, 0.0, 1)/pi4;		// Need J(0) only, extreme narrowing
       if(fabs(xi1xi2*J0) > 1.e-6)
       { if(autoc)
           R_AC_0(T1s, LOp, l, xi1xi2*J0);
         else
           R_CC_0(T1s,T2s,LOp,l,xi1xi2*J0);
       }
       break;
     }
   return;
   }


void Rmu1mu2(super_op& LOp, const spin_system& sys, gen_op& Ho, double* w,
         double* xi1s, double n1, double* xi2s, double n2, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int l, int type, int level)

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
	//			l     : Rank of the mechanisms
	// 			type  : Relaxation type to compute
	//				   0 = Auto & Cross Correlation
	//				   + = Auto Correlation Only
	//				   - = Cross Correlation Only
	//			level : Relaxation treatment level


//                            LOp += R   (Level)
//                                    1,2  

   {
   double cutoff = 1e-12;
   double xi1, xi2, xi1xi2;			// For specific interaction constants
   coord EA1, EA2;				// For specific Euler angles
   double c1s[5];				// Set up 5 coefficients interaction 1
   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[2*l+1]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[(2*l)+1];
   T2s = new gen_op[(2*l)+1];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int hs = sys.HS();				// Total system Hilbert space
   int m=0, k=0;				// Angular momentum component indices
   for(int i=0; i<n1; i++)			// Sum over spins or dipoles i (mu1)
     {
     xi1 = xi1s[i];				// Get i (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       EA1 = (A1[i]).PASys_EA();	 	// Get i (mu1) space tensor Euler angles
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for i
       for(m=-l,k=0; m<=l; m++,k++)		// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of Ho
         T1s[k] = gen_op((T1[i]).component(2,m));
//         T1s[k] = gen_op(((T1[i]).component(2,m)).matrix());
//         T1s[k] = gen_op(T1[i].component(2,m),
//                               sys.get_basis());
         T1s[k].Op_base(Ho);
         }
       for(int j=0; j<n2; j++)			// Sum over spins or dipoles j (mu2)
         {
         if((i==j) && (type >= 0))		// Auto-correlation term
           {
           xi1xi2 = xi1*xi1;			// xi(mu1,i)*xi(mu1,i)
           if(abs(level) <2)
             {
             w0 = 0.0;				// Need only zero, single, & double
             w1 = sys.gamma(i)/GAMMA1H;	 	// quantum transition frequencies
             w1 *= sys.Omega()*1.0e6;
             w2 = 2.0*w1;
             }
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             Rmumu(LOp, T1s, T1s, w, hs, taus,
                c1s,c1s,xi1xi2,w0,w1,w2,l,level,1);
           }
         else if((i!=j) && (type<=0))		// Cross-correlation term
           {
           xi2 = xi2s[j];;			// Get j (mu2) interaction constant
           xi1xi2 = xi1*xi2;			// xi(mu1,i)*xi(mu2,j)
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             {
             EA2 = (A2[j]).PASys_EA();		// Get spin j (mu2) space tensor Euler angles
             Jcoeffs(c2s, EA2, chi);		// Set the 5 coefficients for j
             for(m=-l,k=0; m<=l; m++,k++)	// Put spin tensor for j into a
               {				// vector of operators in basis of Ho
               T2s[k] = gen_op(T2[j].component(2,m));
//               T2s[k] = gen_op(T2[j].component(2,m),
//                                   sys.get_basis());
               T2s[k].Op_base(Ho);
               }
             if(abs(level) <2)
               {
// sosi Now trouble with these frequencies
               wi = sys.gamma(i)/GAMMA1H; 	// Need only zero, single, & double
               wi *= sys.Omega()*1.0e6; 	// quantum transition frequencies
               wj = sys.gamma(j)/GAMMA1H;
               wj *= sys.Omega()*1.0e6;
               w0 = wi-wj;
               w1 *= wi;			// !! Don't know which to use here!!
               w2 = wi+wj;
               }
             Rmumu(LOp, T1s, T2s, w, hs, taus,
               c1s,c2s,xi1xi2,w0,w1,w2,l,level,0);
             }
           }
         }					// Increment second spin
       }
     }						// Increment first spin
   return;
   }

// --------- Spin-Pair with Spin-Pair Interaction Functions -------------

void Rijkl(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

/*	             Two Spin-Pair Rank 2 Mechanisms

                                          >=i     >=k
                      --- ---             --- --- --- ---
                      \   \               \   \   \   \
               LOp += /   /    R        = /   /   /   /   R
  	              --- ---   mu1,mu2   --- --- --- ---  ij,kl
                      mu1 mu2              i   j   k   l                   */

   {
// sosi added 3/11/96
int het = sys.heteronuclear();			// Flag if system is heteronuclear
int rank=2;
// sosi: this is still dipole specific!!
matrix theta = sys.thetas();			// Get dipole theta values (radians) 
matrix phi = sys.phis(); 			// Get dipole phi values (radians)
double alphaij, betaij;			// Two Euler angles for ij dipole
double alphakl, betakl;			// Two Euler angles for kl dipole
   double xi1, xi2, xi1xi2;			// Specific interact. constants
   coord EA1, EA2;				// For specific Euler angles
   double c1s[5];				// The 5 coeffs. interaction 1
   double c2s[5];				// The 5 coeffs. interaction 2
   gen_op *T1s, *T2s;				// The compiler didn't like
   T1s = new gen_op[5]; 			// gen_op T1s[5] so we do these
   T2s = new gen_op[5];				// steps, spin tensor components
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
   int m;					// z-component of ang. momentum
   int ij = 0;					// Sum over all dipolar pairs
   int kl = 0;					// and build relaxation matrix
   for(int i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)
       {
       xi1 = Re(xi1s.get(i,j));			// Dipolar interact. const. i&j
alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
//       Jcoeffs(c1s, alphaij, betaij, 0.0, chi);	// Set the 5 coefficients for ij
EA1.xyz(alphaij, betaij, 0.0);
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for ij
       for(m=-2; m<3; m++)			// Put i,j spin tensor into a
         {					// vector of Ops in basis of Ho
         T1s[m+2] = gen_op(T1[ij].component(2,m));
//         T1s[m+2] = gen_op(T1[ij].component(2,m),
//                          sys.get_basis());
         T1s[m+2].Op_base(Ho);
         }
       kl = 0;					// Set dipole kl count to zero
       for(int k=0; k<ns-1; k++)		// Sum spin pairs k&l: dipole kl
         for(int l=k+1; l<ns; l++)
           {
           if((ij == kl) && (type >= 0))	// Auto-correlation term
             {
             xi1xi2 = xi1*xi1;			// xi(ij)*xi(ij)
             if(abs(level) <2)
               {
               w0 = 0.0;			// Need only 0, 1, & 2
               w1 = sys.gamma(i)/GAMMA1H; 	// quantum transition freqs.
               w1 *= sys.Omega()*1.0e6;		// (in the laboratory frame)
               w2 = 2.0*w1;
               }
             Rmumu(LOp,T1s,T1s,w,hs,taus,
                c1s,c1s,xi1xi2,w0,w1,w2,rank,level,1,het);
             }
           else if((ij != kl) && (type <= 0))	// Cross-correlation term
             {
//           xi2 = Re(xi1s.get(k,l));		// Dip. interact. constant k&l
             xi2 = Re(xi2s.get(k,l));		// Dip. interact. constant k&l
             xi1xi2 = xi1*xi2;			// xi(ij)*xi(kl)
alphakl = Re(phi.get(k,l));	// Get the alpha Euler angle kl 
betakl = Re(theta.get(k,l));	// Get the beta Euler angle kl
//             Jcoeffs(c2s, alphakl, betakl, 	// Set the 5 coefficients for kl
//				   0.0,  chi);  // where gammakl = 0 dipolar
EA2.xyz(alphakl, betakl, 0.0);
             Jcoeffs(c2s, EA2, chi);		// Set the 5 coefficients for kl
             for(m=-2; m<3; m++)		// Put k,l spin tensor into a
               {				// vector of Ops in basis of Ho
               T2s[m+2] = gen_op(T2[kl].component(2,m));
//               T2s[m+2] = gen_op(T2[kl].
//                            component(2,m),
//                              sys.get_basis());
               T2s[m+2].Op_base(Ho);
               }
             if(abs(level) <2)
               {
               wi = sys.gamma(i)/GAMMA1H;	// Need only 0, 1, and 2
               wi *= sys.Omega()*1.0e6; 	// quantum transition freqs
               wj = sys.gamma(j)/GAMMA1H;	// (in the laboratory frame)
               wj *= sys.Omega()*1.0e6; 	// !! Which to use here??!!
               w0 = wi-wj;
               w1 *= wi;
               w2 = wi+wj;
               }
             Rmumu(LOp,T1s,T2s,w,hs,taus,
               c1s,c2s,xi1xi2,w0,w1,w2,rank,level,0, het);
             }
           kl++;				// Increment second dipole
           }
       ij++;					// Increment first dipole
       }
   gen_op Op;					// Psuedo destruction of Ops
   if(T1s)					// Delete the first 
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
   if(A1 != NULL) type = 0;			// Compiler likes this used
   if(A2 != NULL) type = 0;			// Compiler likes this used
   }


void Rijkl(super_op& LOp, const spin_system& sys, gen_op& Ho, double* w,
                  matrix& xi1s, matrix& xi2s, spin_T* T1, spin_T* T2,
                                    double tau, int type, int level)

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
	// Note			      : Uses spherical top density functions

/*	             Two Spin-Pair Rank 2 Mechanisms

                      --- ---             --- --- --- ---
                      \   \               \   \   \   \
               LOp += /   /    R        = /   /   /   /   R
  	              --- ---   mu1,mu2   --- --- --- ---  ij,kl
                      mu1 mu2              i   j   k   l                  */

   {
   int rank=2;
   double xi1, xi2, xi1xi2=0;			// For specific interaction constants
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[5];
   T2s = new gen_op[5];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
//   int ls = hs*hs;				// Total system Liouville space
   int m;					// z-component of ang. momentum
   int ij = 0;					// Sum over all dipolar pairs
   int kl = 0;					// and build relaxation matrix
   for(int i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)
       {
       xi1 = Re(xi1s.get(i,j));			// Get the interaction constant for i&j
       for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
         {					// vector of operators in basis of Ho
         T1s[m+2] = gen_op(T1[ij].component(2,m));
//         T1s[m+2] = gen_op(T1[ij].component(2,m),
//                          sys.get_basis());
         T1s[m+2].Op_base(Ho);
         }
       kl = 0;					// Set the spin pair kl count to zero
       for (int k=0; k<ns-1; k++)		// Sum spin pairs k&l: pair kl
         for (int l=k+1; l<ns; l++)
           {
           if((ij == kl) && (type >= 0))	// Auto-correlation term
             {
             xi1xi2 = xi1*xi1;			// xi(ij)*xi(ij)
             if(abs(level) <2)
               {
               w1 = sys.gamma(i)/GAMMA1H; 	// Larmor frequency
               w1 *= sys.Omega()*1.0e6;
               }
             Rmumu(LOp,T1s,T1s,w,hs,tau,
                   xi1xi2,0.0,w1,2*w1,rank,level,1);
             }
           else if((ij != kl) && (type <= 0))	// Cross-correlation term
             {
             xi2 = Re(xi2s.get(k,l));		// Dipolar interaction constant k&l
//             xi2 = Re(xi1s.get(k,l));		// Dipolar interaction constant k&l
             xi1xi2 = xi1*xi2;			// xi(ij)*xi(kl)
             for(m=-2; m<3; m++)		// Put spin tensor for k,l into a
               {				// vector of operators in basis of Ho
               T2s[m+2] = gen_op(T2[kl].
                            component(2,m));
//               T2s[m+2] = gen_op(T2[kl].
//                            component(2,m),
//                              sys.get_basis());
               T2s[m+2].Op_base(Ho);
               }
             if(abs(level) <2)
               {
               wi = sys.gamma(i)/GAMMA1H;	// Need only zero, single and double
               wi *= sys.Omega()*1.0e6; 	// quantum transition frequencies
               wj = sys.gamma(j)/GAMMA1H;
               wj *= sys.Omega()*1.0e6; 	// !! Don't know which to use here!!
               w0 = wi-wj;
               w1 *= wi;
               w2 = wi+wj;
               }
             Rmumu(LOp,T1s,T2s,w,hs,tau,xi1xi2,
                              w0,w1,w2,rank,level,0);
             }
           kl++;				// Increment second dipole
           }
       ij++;					// Increment first dipole
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
   }

// -------------------- Spin with Spin Functions -------------------

void Rij(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

/*	          Two Single Spin Rank 2 Mechanisms

                      --- ---             --- ---
                      \   \               \   \
               LOp += /   /    R        = /   /   R
  	              --- ---   mu1,mu2   --- ---  i,j
                      mu1 mu2              i   j                       */

   {
// sosi added 3/11/96
int het = sys.heteronuclear();			// Flag if system is heteronuclear
   int rank=2;
   double cutoff = 1e-12;
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
   for(int i=0; i<ns; i++)			// Sum over spins i (mu1)
     {
     xi1 = Re(xi1s.get(i,i));			// Get spin i (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       EA1 = (A1[i]).PASys_EA();	 		// Get spin i (mu1) space tensor Euler angles
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for i
       for(m=-2; m<3; m++)			// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of Ho
         T1s[m+2] = gen_op(T1[i].component(2,m));
//         T1s[m+2] = gen_op(T1[i].component(2,m),
//                               sys.get_basis());
         T1s[m+2].Op_base(Ho);
         }
// sosi this should be j=i, don't have ij & ji in sum!
       for(int j=0; j<ns; j++)			// Sum over spins j
         {
         if((i==j) && (type >= 0))		// Auto-correlation term
           {
           xi1xi2 = xi1*xi1;			// xi(mu1,i)*xi(mu1,i)
           if(abs(level) <2)
             {
             w0 = 0.0;				// Need only zero, single, & double
             w1 = sys.gamma(i)/GAMMA1H;	// quantum transition frequencies
             w1 *= sys.Omega()*1.0e6;
             w2 = 2.0*w1;
             }
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             Rmumu(LOp, T1s, T1s, w, hs, taus,
                c1s,c1s,xi1xi2,w0,w1,w2,rank,level,1,het);
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
               {				// vector of operators in basis of Ho
               T2s[m+2] = gen_op(T2[j].component(2,m));
//               T2s[m+2] = gen_op(T2[j].component(2,m),
//                                   sys.get_basis());
               T2s[m+2].Op_base(Ho);
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
               }
             Rmumu(LOp, T1s, T2s, w, hs, taus,
               c1s,c2s,xi1xi2,w0,w1,w2,rank,level,0,het);
             }
           }
         }					// Increment second spin
       }
     }						// Increment first spin
   return;
   }


void Rij(super_op& LOp, const spin_system& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, spin_T* T1, spin_T* T2,
                             double tau, int l, int type, int level)

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

/*	          Two Single Spin Rank 2 Mechanisms

                      --- ---             --- ---
                      \   \               \   \
               LOp += /   /    R        = /   /   R
  	              --- ---   mu1,mu2   --- ---  i,j
                      mu1 mu2              i   j                     */

   {
   double cutoff = 1e-12;			// Strength cutoff value
   double xi1, xi2, xi1xi2;			// For specific interaction constants
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[2*l+1];
   T2s = new gen_op[2*l+1];
   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
//   int ls = hs*hs;				// Total system Liouville space
   int m;					// z-component of ang. momentum
   for(int i=0; i<ns; i++)			// Sum over spins i (mu1)
     {
     xi1 = Re(xi1s.get(i,i));			// Get spin i (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       for(m=-l; m<=l; m++)			// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of Ho
         T1s[m+l] = gen_op(T1[i].component(l,m));
//         T1s[m+l] = gen_op(T1[i].component(l,m),
//                               sys.get_basis());
         T1s[m+l].Op_base(Ho);
         }
       for(int j=0; j<ns; j++)			// Sum over spins j
         {
         if((i==j) && (type >= 0))		// Auto-correlation term
           {
           xi1xi2 = xi1*xi1;			// xi(mu1,i)*xi(mu1,i)
           if(abs(level) <2)
             {
             w0 = 0.0;				// Need only zero, single, & double
             w1 = sys.gamma(i)/GAMMA1H;	 	// quantum transition frequencies
             w1 *= sys.Omega()*1.0e6;
             w2 = 2.0*w1;
             }
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             Rmumu(LOp, T1s, T1s, w, hs, tau,
                  xi1xi2,0.0,w1,2*w1,l,level,1);
           }
         else if((i!=j) && (type<=0))		// Cross-correlation term
           {
           xi2 = Re(xi2s.get(j,j));		// Get spin j (mu2) interaction constant
           xi1xi2 = xi1*xi2;			// xi(mu1,i)*xi(mu2,j)
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             {
             for(m=-l; m<=l; m++)		// Put spin tensor for j into a
               {				// vector of operators in basis of Ho
               T2s[m+l] = gen_op(T2[j].component(l,m));
//               T2s[m+l] = gen_op(T2[j].component(l,m),
//                                   sys.get_basis());
               T2s[m+l].Op_base(Ho);
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
               }
             Rmumu(LOp, T1s, T2s, w, hs, tau,
                     xi1xi2,w0,w1,w2,l,level,0);
             }
           }
         }					// Increment second spin
       }
     }						// Increment first spin
   return;
   }


// --------------- Spin-Pair with Spin Functions -------------------

void Rijk(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

// sosi 11-7-93: BUG FOUND - the index ij was not being incremented!!!

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

/*	       Mechanism 1 = Spin-Pair; Mechanism 2 = Spin

                      --- ---             --- --- ---
                      \   \               \   \   \
               LOp += /   /    R        = /   /   /   R
  	              --- ---   mu1,mu2   --- --- ---  ij,k
                      mu1 mu2              i   j   k              */

   {
// sosi added 3/11/96
int het = sys.heteronuclear();			// Flag if system is heteronuclear
   int rank=2;
// sosi: these three variables are still dipole specific!!
   matrix theta = sys.thetas();			// Get dipole theta values (radians) 
   matrix phi = sys.phis(); 			// Get dipole phi values (radians)
   double alphaij, betaij;			// Two Euler angles for ij dipole

   double xi1, xi2, xi1xi2;			// For specific interaction constants
//   coord EA1, EA2;				// For specific Euler angles
   coord EA2;					// For specific Euler angles
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
   double cutoff = 0;
//   double cutoff = 1e-12;
   int ij = 0;					// Sum over all dipolar pairs
   for(int i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
// sosi - this is the 11/7/93 bug: for(int j=i+1; j<ns; j++)
     for(int j=i+1; j<ns; j++, ij++)
       {
       xi1 = Re(xi1s.get(i,j));			// Dipolar interaction constant i&j
       alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
       betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
       Jcoeffs(c1s, alphaij, betaij, 0.0, chi);	// Set the 5 coefficients for ij
       for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
         {					// vector of operators in basis of Ho
         T1s[m+2] = gen_op(T1[ij].component(2,m));
//         T1s[m+2] = gen_op(T1[ij].component(2,m),
//                          sys.get_basis());
         T1s[m+2].Op_base(Ho);
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
             {					// vector of operators in basis of Ho
             T2s[m+2] = gen_op(T2[k].component(2,m));
//             T2s[m+2] = gen_op(T2[k].component(2,m),
//                                 sys.get_basis());
             T2s[m+2].Op_base(Ho);
             }
           if(abs(level) <2)
             {
             wi = sys.gamma(i)/GAMMA1H; 	// Need only zero, single, & double
             wi *= sys.Omega()*1.0e6;	 	// quantum transition frequencies
             wj = sys.gamma(j)/GAMMA1H;
             wj *= sys.Omega()*1.0e6;
             w0 = wi-wj;
             w1 *= wi;				// !! Don't know which to use here!!
             w2 = wi+wj;
             }
           Rmumu(LOp, T1s, T2s, w, hs, taus,
               c1s,c2s,xi1xi2,w0,w1,w2,rank,level,0,het);
           }
         } 					// Increment spin k (mu2)
       }					// Increment pair ij (mu1)
   return;
   if(A1 != NULL) type = 0;			// Compiler likes these used
   }


void Rkij(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

/*	       Mechanism 1 = Spin; Mechanism 2 = Spin-Pair

                      --- ---             --- --- ---
                      \   \               \   \   \
               LOp += /   /    R        = /   /   /   R
  	              --- ---   mu1,mu2   --- --- ---  k,ij
                      mu1 mu2              k   i   j                   */

   {
// sosi added 3/11/96
int het = sys.heteronuclear();			// Flag if system is heteronuclear
// sosi: these three variables are still dipole specific!!
   int rank=2;
   matrix theta = sys.thetas();			// Get dipole theta values (radians) 
   matrix phi = sys.phis(); 			// Get dipole phi values (radians)
   double alphaij, betaij;			// Two Euler angles for ij dipole

   double xi1, xi2, xi1xi2;			// For specific interaction constants
   coord EA1 ;					// For specific Euler angles
//   coord EA1, EA2;				// For specific Euler angles
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
   for(int k=0; k<ns; k++)			// Sum over spins k (mu1)
     {
     xi1 = Re(xi1s.get(k,k));			// Get spin k (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       EA1 = (A1[k]).PASys_EA(); 		// Get spin k (mu1) space tensor Euler angles
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for k
       for(m=-2; m<3; m++)			// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of Ho
         T1s[m+2] = gen_op(T1[k].component(2,m));
//         T1s[m+2] = gen_op(T1[k].component(2,m),
//                               sys.get_basis());
         T1s[m+2].Op_base(Ho);
         }
       ij = 0;					// Set dipole count to zero
       for(int i=0; i<ns-1; i++)		// Sum spin pairs i&j: dipole ij (mu2)
         for(int j=i+1; j<ns; j++)
           {
           xi2 = Re(xi2s.get(i,j));		// Dipolar interaction constant k&l
           xi1xi2 = xi1*xi2;			// xi(k)*xi(ij)
           alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
           betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
           Jcoeffs(c2s,alphaij,betaij,0.0,chi);	// Set the 5 coefficients for ij
           for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
             {					// vector of operators in basis of Ho
             T2s[m+2] = gen_op(T2[ij].
                            component(2,m));
//             T2s[m+2] = gen_op(T2[ij].
//                            component(2,m),
//                              sys.get_basis());
               T2s[m+2].Op_base(Ho);
             }
           if(abs(level) <2)
             {
             wi = sys.gamma(i)/GAMMA1H;	// Need only zero, single and double
             wi *= sys.Omega()*1.0e6;	 	// quantum transition frequencies
             wj = sys.gamma(j)/GAMMA1H;
             wj *= sys.Omega()*1.0e6;	 	// !! Don't know which to use here!!
             w0 = wi-wj;
             w1 *= wi;
             w2 = wi+wj;
             }
           Rmumu(LOp, T1s, T2s, w, hs, taus,
               c1s,c2s,xi1xi2,w0,w1,w2,rank,level,0,het);
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
   if(A2 != NULL) type = 0;			// Compiler likes these used
   }

// ______________________________________________________________________
// *RELAXATION SUPEROPERATOR GENERAL FUNCTIONS   DYNAMIC FREQUENCY SHIFT*
// ______________________________________________________________________

// ---------------- Routing to Specific Level Computation ---------------

void Rmumuds(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
               double w0, double w1, double w2, int level, int autoc, int het)
// sosi: the rank is set to 2 which is inappropriate!

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

   {
   double secapp = 1.e6;			// Heteronuclear secular approx. (Hz)
   int rank = 2;
   matrix J12;					// Used for J's levels > 1
   double J0, J1, J2;				// Used for J's levels 0 & 1
   switch(level)
     {
     case 4:					// Level 4 mu1-mu2: element by element
       J12 = Q_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spectral densities
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       if(!het) R_4(LOp, rank, T1s, T2s, J12);	// This if homonuclear spin system
       else R_3(LOp,w,rank,T1s,T2s,J12,secapp);	// This if heteronuclear spin system
       break;
     case -4:					// Level 4 mu1-mu2: double commutator
       J12 = Q_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spectral densities
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       R_4s(LOp, rank, T1s, T2s, J12);
       break;
     case 3:					// Level 3 mu1-mu2: element by element
       J12 = Q_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spectral densities
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       R_3(LOp, w, rank, T1s, T2s, J12);
       break;
     case -3:					// Level 3 mu1-mu2: double commutator
       J12 = Q_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spectral densities
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       R_3s(LOp, w, rank, T1s, T2s, J12);
       break;
     case 2:					// Level 2 mu1-mu2: element by element
       J12 = Q_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spectral densities
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       R_2(LOp, rank, T1s, T2s, J12);
       break;
     case -2:					// Level 2 mu1-mu2: double commutator
       J12 = Q_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spectral densities
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       R_2s(LOp, rank, T1s, T2s, J12);
       break;
     case 1: 					// Level 1 mu1-mu2: double commutator
       J0 = Q_reduced(taus, c1s, c1s, w0, 1);	// J(0)
       J1 = Q_reduced(taus, c1s, c1s, w1, 1);	// J(w0)
       J2 = Q_reduced(taus, c1s, c1s, w2, 1);	// J(2w0)
       if(autoc)
         R_AC_1(T1s, LOp, rank,
                     J0*xi1xi2,J1*xi1xi2,J2*xi1xi2);
       else
         R_CC_1(T1s, T2s, LOp, rank,
                   J0*xi1xi2,J1*xi1xi2,J2*xi1xi2);
       break;
     case 0:					// Level 0 mu1-mu2: element by element
       J0 = Q_reduced(taus,c1s,c2s,0.0,1);	// Need J(0) only, extreme narrowing
       if(fabs(xi1xi2*J0) > 1.e-6)
         R_0(LOp,rank,T1s,T2s,complex(xi1xi2*J0));
       break;
     default:					// Level 0 mu1-mu2: double commutator
       J0 = Q_reduced(taus,c1s,c2s,0.0,1);	// Need J(0) only, extreme narrowing
       if(fabs(xi1xi2*J0) > 1.e-6)
       { if(autoc)
           R_AC_0(T1s, LOp, rank, xi1xi2*J0);
         else
           R_CC_0(T1s,T2s,LOp,rank,xi1xi2*J0);
       }
       break;
     }
   return;
   }

// --------- Spin-Pair with Spin-Pair Interaction Functions -------------

void Rijklds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

/*	             Two Spin-Pair Rank 2 Mechanisms

                      --- ---             --- --- --- ---
                      \   \               \   \   \   \
               LOp += /   /    R        = /   /   /   /   R
  	              --- ---   mu1,mu2   --- --- --- ---  ij,kl
                      mu1 mu2              i   j   k   l              */

   {
// sosi added 3/11/96
int het = sys.heteronuclear();                  // Flag if system is heteronuclear
//int rank=2;
// sosi: this is still dipole specific!!
matrix theta = sys.thetas();			// Get dipole theta values (radians) 
matrix phi = sys.phis(); 			// Get dipole phi values (radians)
double alphaij, betaij;			// Two Euler angles for ij dipole
double alphakl, betakl;			// Two Euler angles for kl dipole
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
   int ij = 0;					// Sum over all dipolar pairs
   int kl = 0;					// and build relaxation matrix
   int i = 0;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)
       {
       xi1 = Re(xi1s.get(i,j));			// Dipolar interaction constant i&j
alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
//       Jcoeffs(c1s, alphaij, betaij, 0.0, chi);	// Set the 5 coefficients for ij
EA1.xyz(alphaij, betaij, 0.0);
         Jcoeffs(c1s, EA1, chi);		// Set the 5 coefficients for ij
       for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
         {					// vector of operators in basis of Ho
         T1s[m+2] = gen_op(T1[ij].component(2,m));
//         T1s[m+2] = gen_op(T1[ij].component(2,m),
//                          sys.get_basis());
         T1s[m+2].Op_base(Ho);
         }
       kl = 0;					// dipole kl count to zero
       for (int k=0; k<ns-1; k++)		// Sum spin pairs k&l: dipole kl
         for (int l=k+1; l<ns; l++)
           {
           if((ij == kl) && (type >= 0))	// Auto-correlation term
             {
             xi1xi2 = xi1*xi1;			// xi(ij)*xi(ij)
             if(abs(level) <2)
               {
               w0 = 0.0;			// Need only zero, single, & double
               w1 = sys.gamma(i)/GAMMA1H; 	// quantum transition frequencies
               w1 *= sys.Omega()*1.0e6;
               w2 = 2.0*w1;
                }
             Rmumuds(LOp,T1s,T1s,w,hs,taus,
                c1s,c1s,xi1xi2,w0,w1,w2,level,1,het);
             }
           else if((ij != kl) && (type <= 0))	// Cross-correlation term
             {
             xi2 = Re(xi2s.get(k,l));		// Dipolar interaction constant k&l
// sosi - switched above to xi2s from xi1s 1/31/96
             xi1xi2 = xi1*xi2;			// xi(ij)*xi(kl)
alphakl = Re(phi.get(k,l));	// Get the alpha Euler angle kl 
betakl = Re(theta.get(k,l));	// Get the beta Euler angle kl
//             Jcoeffs(c2s, alphakl, betakl, 	// Set the 5 coefficients for kl
//				   0.0,  chi);  // where gammakl = 0 dipolar
EA2.xyz(alphakl, betakl, 0.0);
             Jcoeffs(c2s, EA2, chi);		// Set the 5 coefficients for kl
             for(m=-2; m<3; m++)		// Put spin tensor for k,l into a
               {				// vector of operators in basis of Ho
               T2s[m+2] = gen_op(T2[kl].
                            component(2,m));
//               T2s[m+2] = gen_op(T2[kl].
//                            component(2,m),
//                              sys.get_basis());
               T2s[m+2].Op_base(Ho);
               }
             if(abs(level) <2)
               {
               wi = sys.gamma(i)/GAMMA1H;	// Need only zero, single and double
               wi *= sys.Omega()*1.0e6; 	// quantum transition frequencies
               wj = sys.gamma(j)/GAMMA1H;
               wj *= sys.Omega()*1.0e6; 	// !! Don't know which to use here!!
               w0 = wi-wj;
               w1 *= wi;
               w2 = wi+wj;
               }
             Rmumuds(LOp,T1s,T2s,w,hs,taus,
               c1s,c2s,xi1xi2,w0,w1,w2,level,0,het);
             }
           kl++;				// Increment second dipole
           }
       ij++;					// Increment first dipole
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
   if(A1 != NULL) i = 0;			// Compiler likes this used
   if(A2 != NULL) i = 0;			// Compiler likes this used
   }


// -------------------- Spin with Spin Functions -------------------


void Rijds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

/*	          Two Single Spin Rank 2 Mechanisms

                      --- ---             --- ---
                      \   \               \   \
               LOp += /   /    R        = /   /   R
  	              --- ---   mu1,mu2   --- ---  i,j
                      mu1 mu2              i   j                   */

   {
// sosi added 3/11/96
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
   for(int i=0; i<ns; i++)			// Sum over spins i (mu1)
     {
     xi1 = Re(xi1s.get(i,i));			// Get spin i (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       EA1 = (A1[i]).PASys_EA();		// Get spin i (mu1) space tensor Euler angles
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for i
       for(m=-2; m<3; m++)			// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of Ho
         T1s[m+2] = gen_op(T1[i].component(2,m));
//         T1s[m+2] = gen_op(T1[i].component(2,m),
//                               sys.get_basis());
         T1s[m+2].Op_base(Ho);
         }
       for(int j=0; j<ns; j++)			// Sum over spins j
         {
         if((i==j) && (type >= 0))		// Auto-correlation term
           {
           xi1xi2 = xi1*xi1;			// xi(mu1,i)*xi(mu1,i)
           if(abs(level) <2)
             {
             w0 = 0.0;				// Need only zero, single, & double
             w1 = sys.gamma(i)/GAMMA1H;	 	// quantum transition frequencies
             w1 *= sys.Omega()*1.0e6;
             w2 = 2.0*w1;
             }
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             Rmumuds(LOp, T1s, T1s, w, hs, taus,
                c1s,c1s,xi1xi2,w0,w1,w2,level,1,het);
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
               {				// vector of operators in basis of Ho
               T2s[m+2] = gen_op(T2[j].component(2,m));
//               T2s[m+2] = gen_op(T2[j].component(2,m),
//                                   sys.get_basis());
               T2s[m+2].Op_base(Ho);
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
               }
             Rmumuds(LOp, T1s, T2s, w, hs, taus,
               c1s,c2s,xi1xi2,w0,w1,w2,level,0,het);
             }
           }
         }					// Increment second spin
       }
     }						// Increment first spin
   return;
   }


// --------------- Spin-Pair with Spin Functions -------------------


void Rijkds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

/*	       Mechanism 1 = Spin-Pair; Mechanism 2 = Spin

                      --- ---             --- --- ---
                      \   \               \   \   \
               LOp += /   /    R        = /   /   /   R
  	              --- ---   mu1,mu2   --- --- ---  ij,k
                      mu1 mu2              i   j   k                  */

   {
// sosi added 3/11/96
int het = sys.heteronuclear();                  // Flag if system is heteronuclear
//   int rank=2;
// sosi: these three variables are still dipole specific!!
   matrix theta = sys.thetas();			// Get dipole theta values (radians) 
   matrix phi = sys.phis(); 			// Get dipole phi values (radians)
   double alphaij, betaij;			// Two Euler angles for ij dipole

   double xi1, xi2, xi1xi2;			// For specific interaction constants
//   coord EA1, EA2;				// For specific Euler angles
   coord EA2;					// For specific Euler angles
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
   double cutoff = 0;
//   double cutoff = 1e-12;
   int ij = 0;					// Sum over all dipolar pairs
   for(int i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)
       {
       xi1 = Re(xi1s.get(i,j));			// Dipolar interaction constant i&j
       alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
       betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
       Jcoeffs(c1s, alphaij, betaij, 0.0, chi);	// Set the 5 coefficients for ij
       for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
         {					// vector of operators in basis of Ho
         T1s[m+2] = gen_op(T1[ij].component(2,m));
//       T1s[m+2] = gen_op(T1[ij].component(2,m),
//                          sys.get_basis());
         T1s[m+2].Op_base(Ho);
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
             {					// vector of operators in basis of Ho
             T2s[m+2] = gen_op(T2[k].component(2,m));
//             T2s[m+2] = gen_op(T2[k].component(2,m),
//                                 sys.get_basis());
             T2s[m+2].Op_base(Ho);
             }
           if(abs(level) <2)
             {
             wi = sys.gamma(i)/GAMMA1H; 	// Need only zero, single, & double
             wi *= sys.Omega()*1.0e6;	 	// quantum transition frequencies
             wj = sys.gamma(j)/GAMMA1H;
             wj *= sys.Omega()*1.0e6;
             w0 = wi-wj;
             w1 *= wi;				// !! Don't know which to use here!!
             w2 = wi+wj;
             }
           Rmumuds(LOp, T1s, T2s, w, hs, taus,
               c1s,c2s,xi1xi2,w0,w1,w2,level,0,het);
           }
         } 					// Increment spin k (mu2)
       }					// Increment pair ij (mu1)
   return;
   if(A1 != NULL) type = 0;			// Compiler likes these used
   }


void Rkijds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

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

/*	       Mechanism 1 = Spin; Mechanism 2 = Spin-Pair

                      --- ---             --- --- ---
                      \   \               \   \   \
               LOp += /   /    R        = /   /   /   R
  	              --- ---   mu1,mu2   --- --- ---  k,ij
                      mu1 mu2              k   i   j                 */

   {
// sosi added 3/11/96
int het = sys.heteronuclear();                  // Flag if system is heteronuclear
// sosi: these three variables are still dipole specific!!
//   int rank=2;
   matrix theta = sys.thetas();			// Get dipole theta values (radians) 
   matrix phi = sys.phis(); 			// Get dipole phi values (radians)
   double alphaij, betaij;			// Two Euler angles for ij dipole

   double xi1, xi2, xi1xi2;			// For specific interaction constants
   coord EA1 ;					// For specific Euler angles
//   coord EA1, EA2;				// For specific Euler angles
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
   for(int k=0; k<ns; k++)			// Sum over spins k (mu1)
     {
     xi1 = Re(xi1s.get(k,k));			// Get spin k (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       EA1 = (A1[k]).PASys_EA(); 		// Get spin k (mu1) space tensor Euler angles
       Jcoeffs(c1s, EA1, chi);			// Set the 5 coefficients for k
       for(m=-2; m<3; m++)			// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of Ho
         T1s[m+2] = gen_op(T1[k].component(2,m));
//         T1s[m+2] = gen_op(T1[k].component(2,m),
//                               sys.get_basis());
         T1s[m+2].Op_base(Ho);
         }
       ij = 0;					// Set dipole count to zero
       for(int i=0; i<ns-1; i++)		// Sum spin pairs i&j: dipole ij (mu2)
         for(int j=i+1; j<ns; j++)
           {
           xi2 = Re(xi2s.get(i,j));		// Dipolar interaction constant k&l
           xi1xi2 = xi1*xi2;			// xi(k)*xi(ij)
           alphaij = Re(phi.get(i,j));		// Alpha Euler angle: dipole ij
           betaij = Re(theta.get(i,j));		// Beta Euler angle: dipole kl
           Jcoeffs(c2s,alphaij,betaij,0.0,chi);	// Set the 5 coefficients for ij
           for(m=-2; m<3; m++)			// Put spin tensor for i,j into a
             {					// vector of operators in basis of Ho
             T2s[m+2] = gen_op(T2[ij].
                            component(2,m));
//             T2s[m+2] = gen_op(T2[ij].
//                            component(2,m),
//                              sys.get_basis());
               T2s[m+2].Op_base(Ho);
             }
           if(abs(level) <2)
             {
             wi = sys.gamma(i)/GAMMA1H;	// Need only zero, single and double
             wi *= sys.Omega()*1.0e6;	 	// quantum transition frequencies
             wj = sys.gamma(j)/GAMMA1H;
             wj *= sys.Omega()*1.0e6;	 	// !! Don't know which to use here!!
             w0 = wi-wj;
             w1 *= wi;
             w2 = wi+wj;
             }
           Rmumuds(LOp, T1s, T2s, w, hs, taus,
               c1s,c2s,xi1xi2,w0,w1,w2,level,0,het);
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
   if(A2 != NULL) type = 0;			// Compiler likes these used
   }

#endif						// relaxNMR.cc
