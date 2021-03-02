/* relaxProp.cc ******************************************
**							**
**			G A M M A			**
**						 	**
**	Relaxation Propagators      Implementation	**
**						 	**
**	Copyright (c) 1993			 	**
**	Scott A. Smith				 	**
**							**
**	University of California, Santa Barbara		**
**	Department of Chemistry				**
**	Santa Barbara CA. 93106 USA			**
**						 	**
**      $Header: #
**						 	**
*********************************************************/

/*********************************************************
**						 	**
** 	Description				 	**
**						 	**
**  The functions herein provide relaxation propagators	**
**  and functions which utilize them.			**
**						 	**
*********************************************************/

#ifndef _relax_prop_cc_				// Is this file already included?
#define _relax_prop_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)	// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <BWRRelax/relaxProp.h>		// Include the file header
#include <LSLib/LSprop.h>		// Include Liouville space props

// ______________________________________________________________________
// *************** Specific Density Matrix Superoperators ***************
// ______________________________________________________________________


super_op LOp_sigma(gen_op& sigma)

// Input		sigma   : Density matrix
// Output		LOp     : Superoperator which fullfills
//	                              sigma = LOp * sigma0

//				  for any valid sigma0
//
// Note				: LOp constructed in default Liouville
//				  space and Hilbert space of sigma

//  <a,aa|LOp|b,bb> = del    * <a,aa|sigma|1> = del    * <a|sigma|aa>
//			 b,bb                      b,bb

  {
  int hs = sigma.dim();	                    // Get the Hilbert space
  int ls = hs*hs;                           // Set the Liouville space
  matrix mx(ls, ls, 0.0);                   // Set null Liouville matrix
  super_op LOp(mx, sigma.get_basis());		// Construct empty LOp
  int a,aa,b,bb,aaa=0,bbb=0;
  if(sigma.in_EBR())
    {
    for(a=0; a<hs; a++)				// Calculation for diagonal sigma
      for(aa=0; aa<hs; aa++)
        {
        if(a==aa)
          {
          bbb = 0;
          for(b=0; b<hs; b++)
            for(bb=0; bb<hs; bb++)
              {
              if(b==bb)
                LOp.put(aaa,bbb,sigma(a,a));
              bbb++;
              }
          }
      aaa++;
      }
    return LOp;
    }
  else						// Calculation for general sigma
    {
    for(a=0; a<hs; a++)
      for(aa=0; aa<hs; aa++)
        {
        bbb = 0;
        for(b=0; b<hs; b++)
          for(bb=0; bb<hs; bb++)
            {
            if(b==bb)
              LOp.put(aaa, bbb, sigma(a,aa));
            bbb++;
            }
          }
        aaa++;
        }
    return LOp;
    }


/*
This function is now in R_prop in LSLib module

void set_trace(gen_op& sigma, double tr)

  {
  complex z;
  double acttr = Re(trace(sigma));
  double difftr = acttr - tr;
  double fact;
  int i;
  int hs = sigma.dim();
  if(difftr)
    {
    fact = difftr/double(hs);
    for(i=0; i<hs; i++)
      {
      z = sigma.get(i,i) - fact;
      sigma.put(z, i,i); 
      }
    }
  acttr = Re(trace(sigma));
  }
*/
 
/* This is also now defined in LSLib/LSprop
super_op R_prop(super_op& eLt, gen_op sigmaeq)

// Input		eLt	: Exponential Liouvillian to relax
//				  the density matrix
//			sigmaeq : Density matrix at equilibrium
//				  (or steady state density matrix)
//			Ho      : Operator in Hilbert space
// Output		LOp     : Superoperator which fullfills
//
//	   eLt * {sigma-sigmaeq} + sigmaeq  = LOp * sigma
//
// Note				: LOp constructed in the current basis of
//				  sigmaeq & in the default Liouville basis
// Note				: sigmaeq is copied so trace is unchanged

  {
  eLt.set_HBR();				// Put eLt in default Liouville basis
  super_op LOp(eLt);				// Copy the original superoperator
  eLt.LOp_base(sigmaeq);			// Put sigma into eLt Hilbert basis
  set_trace(sigmaeq, (double)1.0);			// Set the trace to 1
  int hs = sigmaeq.dim();			// Get the Hilbert space size
  int a,aa,b,bb,g,gg,aaa=0,bbb=0,ggg=0;
  complex Rel;
  int nd = sigmaeq.in_EBR();			// Test if sigmaeq is diagonal
// *should add a check for i-matrix so prop is i-matrix returned
//  if(mx.stored_type() != i_matrix_type)
// *actually need to test (1-eLt)*sigmaeq and see if it is diagonal, not sigmaeq!
// should set nd=0;

  for(b=0; b<hs; b++)				// For nd=1, diagonal 
    for(bb=0; bb<hs; bb++)
      {
      if(b==bb) 				// <aaa|LOp|bbb> = <aaa|eLt|bbb> if |bb> != |b>
        {					// so do not make any modifications
        aaa = 0;		
        for(a=0; a<hs; a++)
          for(aa=0; aa<hs; aa++)
            {
            Rel = LOp.get(aaa,bbb);
            if(!nd || a==aa)			// If non-diagonal sigmaeq or |a> = |aa> must
              Rel += sigmaeq.get(a,aa);		// add <aaa|sigmaeq|1> to each <aaa|LOp|bbb>
            ggg=0;
            for(g=0; g<hs; g++)
              for(gg=0; gg<hs; gg++)
                {
                Rel -= eLt.get(aaa,ggg)*sigmaeq.get(g,gg);
                ggg++;
                }
            LOp.put(aaa,bbb,Rel);
            aaa++;
            }
        }
      bbb++;
      }
  return LOp;
  }
*/
/* This is also now defined in LSLib/LSprop
super_op R_prop(super_op& L, gen_op& sigmaeq, double t)

// Input		eLt	: Liouvillian for evolution
//			sigmaeq : Density matrix at equilibrium
//				  (or steady state density matrix)
//			t       : Evolution time
// Output		LOp     : Superoperator which fullfills
//
//	   eLt * {sigma-sigmaeq} + sigmaeq  = LOp * sigma

  {
  super_op LOp = exp(L, -t);		// Exponential Liouvillian
  return R_prop(LOp, sigmaeq);		// Use function overload
  }
*/
//THIS WORKS
//void FID0(gen_op& sigma, gen_op& det, super_op &Lt, block_1D &fid, int np=0)
void FID0(gen_op& sigma, gen_op& det, super_op &Lt, row_vector &fid, int np=0)

	// Input		sigma : Operator to be propagated (initial density mx)
	// 			det   : Detection operator (in trace computation)
	//			Lt     : Relaxation propagator superoperator
	//			fid   :	Data block_1D containing the result
	//	 		np    : Number of points to generate
	// Output		      : None, fid data block_1D filled
	// Note			      : The block_1D fid size is assumed >= np
	// Note			      : Lt should be an exponential Liouvillian
	//				relaxation propagator for the dwell time
	// Note			      : The trace of sigma should be one here

  // fid(k) = Trace { det *  sigma(tk) },  where tk = td*k
  //        = Trace { det *  [exp(-L*tk) * sigma] }
  //        = Trace { det * [ S * exp(-D*tk) * S-1 * sigma] }
  //        = Trace { det * [ S * exp(-D*t)^k * S-1 * sigma] }

  {
  if(np <= 0) np = fid.size();		// Set up the block size if not specified
  Lt.set_EBR();				// Put relaxation propagator into its eigenbase
  Lt.LOp_base(sigma);			// Put sigma into LOp Hilbert basis
  Lt.LOp_base(det);			// Put det into LOp Hilbert basis
  int hs = sigma.dim();			// Get the Hilbert space dimension
  int ls = Lt.dim();			// Get the Liouville space dimension
  basis S = Lt.get_Lbasis();		// Get superoperator basis
  matrix Sm1 = inv(S.U());	// Get superoperator basis inverse
//  matrix S = Lt.get_Lbasis();		// Get superoperator basis
//  matrix Sm1 = inv(S);			// Get superoperator basis inverse
  matrix mx = sigma.get_mx();		// Get the initial density matrix
  matrix Sm1sig = Sm1*mx.resize(ls, 1);	// Sm1|sigma> is ls x 1 vector
  mx = adjoint(det.get_mx());		// Adjoint, since <det| not |det>
  matrix emDtk(ls, ls, i_matrix_type);	// Initial exp(-D*t)^k = I for k=0 
  gen_op Optmp;
  for(int k=0; k<np; k++)		// Begin computing the FID points 
    {
//    mx = (S * emDtk) * Sm1sig;		// mx = S * exp(-D*t)^k * S-1 * sigma
    mx=((S.U())*emDtk)*Sm1sig;// mx = S * exp(-D*t)^k * S-1 * sigma
    mx = mx.resize(hs, hs);
    Optmp = gen_op(mx);
    fid.put(trace(det, Optmp),k);	// Compute FID point k
    emDtk *= Lt.get_mx(); 		// Determine exp(-D*t)^k
//cout << "\n" << k << ". Trace = " << fid(k);
    }
  }


// THIS WORKS
// to rowvector
void FID1(gen_op& sigma, gen_op& det, super_op &L, row_vector &fid, int np=0)

	// Input		sigma : Operator to be propagated (initial density mx)
	// 			det   : Detection operator (in trace computation)
	//			L     : Relaxation propagator superoperator
	//			fid   :	Data block_1D containing the result
	//	 		np    : Number of points to generate
	// Output		      : None, fid data block_1D filled
	// Note			      : The block_1D fid size is assumed >= np
	// Note			      : L should be an exponential Liouvillian
	//				relaxation propagator for the dwell time
	// Note			      : The trace of sigma should be one here

  // fid(k) = Trace { det *  sigma(tk) },  where tk = td*k
  //        = Trace { det *  [exp(-L*tk) * sigma] }
  //        = Trace { det * [ S * exp(-D*tk) * S-1 * sigma] }
  //        = Trace { det * [ S * exp(-D*t)^k * S-1 * sigma] }

  {
  if(np <= 0) np = fid.size();		// Set up the block size if not specified
  L.set_EBR();				// Put relaxation propagator into its eigenbase
  L.LOp_base(sigma);			// Put sigma into LOp Hilbert basis
  L.LOp_base(det);			// Put det into LOp Hilbert basis
//  int hs = sigma.dim();		// Get the Hilbert space dimension
  int ls = L.dim();			// Get the Liouville space dimension
  basis S = L.get_Lbasis();		// Get superoperator basis
  matrix Sm1 = inv(S.U());	// Get superoperator basis inverse
//  matrix S = L.get_Lbasis();		// Get superoperator basis
//  matrix Sm1 = inv(S);			// Get superoperator basis inverse

  matrix mx = sigma.get_mx();
  matrix Sm1sig = Sm1*mx.resize(ls, 1);	// Sm1|sigma> is ls x 1 vector

  matrix emDtk(ls, ls, i_matrix_type);
  complex tr;
  matrix dmx = adjoint(det.get_mx());	// Adjoint, since <det| not |det>
  dmx = dmx.resize(ls,1);
  for(int k=0; k<np; k++)
    {
    mx=(S.U())*emDtk*Sm1sig;	// |mx> = S*emDtk*Sm1*sigma
//    mx = S * emDtk * Sm1sig;		// |mx> = S*emDtk*Sm1*sigma
    tr = 0.0;
    for(int l=0; l<ls; l++)
      tr += dmx.get(l,1)*mx.get(l,1);
    fid.put(tr,k);
//cout << "\n" << k << ". Trace = " << tr;
    emDtk *= L.get_mx(); 
    }
  }


//THIS WORKS
// to rowvec
/*This function is now in LSLib/LSacquire
void FID(gen_op& sigma, gen_op& det, super_op &L, row_vector &fid, int np=0)

	// Input		sigma : Operator to be propagated (initial density mx)
	// 			det   : Detection operator (in trace computation)
	//			L     : Relaxation propagator superoperator
	//			fid   :	Data block_1D containing the result
	//	 		np    : Number of points to generate
	// Output		      : None, fid data block_1D filled
	// Note			      : The block_1D fid size is assumed >= np
	// Note			      : L should be an exponential Liouvillian
	//				relaxation propagator for the dwell time
	// Note			      : For sigma, Tr{sigma} = 1 in this routine
	// Note			      : In superoperator formalism, the trace is
	//				defined as the following whre t=adjoint and
	//				(|A>,|B>) is the scalar product of |A> and |B>
	//				            t         t
	//				Tr(A*B) = <A |B> = (|A >, |B>)

  // fid(k) = Trace { det * sigma(tk) },  where tk = td*k
  //              t                   t
  //        = <det | sigma(tk)> = <det | exp(-L*tk) * sigma>
  //              t
  //        = <det | S * exp(-D*tk) * S-1 |sigma>
  //              t
  //        = <det * S| * exp(-D*t)^k |S-1 * sigma> 

  {
  if(np <= 0) np = fid.size();		// Set up the block size if not specified
  L.set_EBR();				// Put relaxation propagator into its eigenbase
  L.LOp_base(sigma);			// Put sigma into LOp Hilbert basis
  L.LOp_base(det);			// Put det into LOp Hilbert basis
//  int hs = sigma.dim();		// Get the Hilbert space dimension
  int ls = L.dim();			// Get the Liouville space dimension
//  matrix S = L.get_Lbasis();		// Get superoperator basis
//  matrix Sm1 = inv(S);			// Get superoperator basis inverse
  basis S = L.get_Lbasis();		// Get superoperator basis
  matrix Sm1 = inv(S.U());	// Get superoperator basis inverse

  matrix mx = sigma.get_mx();
  matrix Sm1sig = Sm1*mx.resize(ls, 1);	// Sm1|sigma> is ls x 1 vector

  mx = adjoint(det.get_mx());		// Start the formation of the detection ket
  mx = mx.resize(ls,1);			//          t                        t 
  mx = transpose(mx);			// Form <det | from transpose of |det >
					//          t            t
//  matrix detS = mx*S;			// Form <det  * S| = <det |S
  matrix detS = mx*S.U();	// Form <det  * S| = <det |S
  matrix emDtk(ls, ls, i_matrix_type);
  for(int k=0; k<np; k++)
    {      				//               -1
    mx = detS * emDtk * Sm1sig;		// |mx> = emDtk|S  * sigma>
    fid(k) = mx.get(0,0);
//cout << "\n" << k << ". Trace = " << fid(k);
    emDtk *= L.get_mx(); 
    }
  }
*/

#endif /* __RELAX_PROP_CC__ */

