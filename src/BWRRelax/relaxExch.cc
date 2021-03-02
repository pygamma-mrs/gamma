/* relaxExch.cc **************************************************
**								**
**                           G A M M A 				**
**                                				**
**      NMR Exchange                    Implementation          **
**								**
**      Copyright (c) 1995					**
**      Scott A. Smith						**
**      Nikolai Skrynnikov					**
**                                				**
**      National High Magnetic Field Laboratory			**
**      1800 E. Paul Dirac Drive				**
**      Tallahassee, Florida, 32310				**
**                                				**
**      $Header: $
**                                				**
**                                				**
*****************************************************************/

/*****************************************************************
**								**
** Description							**
**								**
** The following functions herein provide easy access to many	**
** common exchange superoperators.  Typically, mutual exchange	**
** functions will demand a sys_dynamic variable as a argument.	**
** For non-mutual exchange functions a multi_sys variable is	**
** necessary.  							**
**								**
*****************************************************************/

#ifndef _relax_EXCH_cc_			// Is this file already included?
#define _relax_EXCH_cc_ 1		// If no, then remember that it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#endif

# include <BWRRelax/relaxExch.h>	// Include the header
# include <LSLib/SuperOp.h>		// Include superoperators
# include <LSLib/sys_dynamic.h>		// Include dynamic systems

super_op Rex(const sys_dynamic& sys)

	// Input		sys	: Dynamic spin system
	// Output		Rex     : Exchange superoperator
	// Note			        : The exchange information
	//				  matrix provided by sys
	//				  has an unique structure
	// Note				: The superoperator returned
	//				  from this function is for
	//				  intramolecular mutual exchange.

//           ^	                                  ^      ^
//           ^             ^     ^      ^   ^     ^      ^
//           X    = K  * [ K   x K    - E x E ] = K    - E
//            i,j    ex     i,j   i,j              i,j

   {
//                                                      ^
//                                  ^                   ^
//               Construct operator E and superoperator E

   int hs = sys.HS();				// Get Hilbert space
   gen_op E(matrix(hs,hs,i_matrix_type));	// Get E (identity) operator
   return Rex(sys, E);				// Use function overload
   }


super_op Rex(const sys_dynamic& sys, gen_op Op)

	// Input		sys	: Dynamic spin system
	// 			Op	: Operator providing basis
	// Output		Rex     : Exchange superoperator
	// Note			        : The exchange information
	//				  matrix provided by sys
	//				  has an unique structure
	// Note				: The superoperator returned
	//				  from this function is for
	//				  intramolecular mutual exchange.

//           ^	                                  ^      ^
//           ^             ^     ^      ^   ^     ^      ^
//           X    = K  * [ K   x K    - E x E ] = K    - E
//            i,j    ex     i,j   i,j              i,j

   {
//                                                      ^
//                                  ^                   ^
//               Construct operator E and superoperator E

   int hs = sys.HS();				// Get Hilbert space
   int ns = sys.size();				// Number of spins
   gen_op E(matrix(hs,hs,i_matrix_type));	// Get E (identity) operator
   E.Op_base(Op);				// Set E basis that of Op
   super_op UtransE = U_transform(E);		// Get ExE superoperator

//    Construct Permutation Matrices K    Based on Kex Info Matrix
//                                    i,j

   matrix Kmx = sys.Kex();			// Get Exchange info
   int numbK = Kmx.rows();			// Number of K's defined
   int ninex;					// Number spins in exchange
   super_op X, Xa;				// Tot., Kex exchange superops
   double Kex;					// Rate for exchange process
   int i=0,j=0,l=0;					// Spin indices
   matrix Mz = sys.qStates();			// <bf|Mz|spin> array
   matrix Pid = matrix (ns,ns, i_matrix_type);  // Identity matrix
   matrix P, PMz; 				// Permutation mx, permuted Mz
   matrix K, K0=matrix(hs,hs,0.0,d_matrix_type);// Exch. mx in Hilbert space
   gen_op KOp;

   for(int k=0; k<numbK; k++)			// Loop over exchange rates
     {
     Kex = Kmx.getRe(k,ns);			// Get exchange process rate
     ninex = int(Kmx.getRe(k,ns+1));		// Get # of spins in exchange
     P = Pid;					// Start with P=I
     for(l=0; l<ninex-1; l++)			// Loop involved spin pairs
       {
       i = int(Kmx.getRe(k,l));			// Get 1st spin in pair
       j = int(Kmx.getRe(k,l+1));		// Get 2nd spin in pair
       P.put(complex1,j,i);			// Set the P elements
       P.put(complex0,i,i);			// on upper half (i->j)
       }
     i = j;					// Close (1,2,..n) cycle by
     j = int(Kmx.getRe(k,0)); 			// setting last 1st, 1st last
     P.put(complex1,j,i);			// Set the P elements for this
     P.put(complex0,i,i);			// Last exchange

//   P Is A "Permutation" Array For All Spins Exchanging Under Process With Rate Kex
//   For A 2-Spin Exchange, P Is A True Permutation Array.  For A Multi-Spin Exchange
//   Process P May Have Multiple Elements On a Single Row And Column, Each Row and
//   Column Being Normalized.  This Implies One Spin Going To Multiple Sites.

     PMz = Mz*P;				// Get permuted Mz (basis functions)
     K = K0;					// Zero exchange mx (Hilbert space)
     int match, isp, bf, pbf;
     for(bf=0; bf<hs; bf++)			// Loop over the basis functions
       for(pbf=0; pbf<hs; pbf++)		// Loop permuted basis functions
         {
         match=1;				// Assume functions bf & pbf match
         for(isp=0; isp<ns && match; isp++)	// Loop through spins and see if
           if(Mz.get(bf,isp)!=PMz.get(pbf,isp))	// they do indeed match.  All spin
	     match=0;				// states must be the same
         if(match)				// For a match, add the appropriate
           { 					// terms to exchange mx (Hilbert)
           K.put_h(K.get(pbf,bf)+0.5,pbf,bf);	// This should keep K as Hermitian 
           if(pbf == bf)
             K.put_h(K.get(pbf,bf)+0.5,pbf,bf);
           }
         }

     KOp = gen_op(K);				// Switch K to gen_op
     KOp.Op_base(Op);				// Set its basis to that of Op
     Xa = U_transform(KOp);			// Form exchange superoperator
     Xa -= UtransE;				// Form exchange superoperator
     Xa *= -Kex; 				// Multiply with exchange rate
     X += Xa;					// Add to total exchange superop
     }
   return X;
   }


#endif						// relaxExch.cc

