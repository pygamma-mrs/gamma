/* acquire.cc *******************************************-*-c++-*-
**								**
**      	        G A M M A				**
**								**
**	  Class Acquire	               Implementation		**
**						 		**
**	  Copyright (c) 1991, 1992, 1993		 	**
**	  Tilo Levante, Scott Smith		 		**
**	  Eidgenoessische Technische Hochschule		 	**
**	  Labor fur physikalische Chemie		 	**
**	  8092 Zurich / Switzerland			 	**
**								**
**      $Header: 
**								**
*****************************************************************/

/*****************************************************************
**								**
**  This file contains the implementation of class acquire, a	**
**  data type which facilitates the repetitive computation of	**
**  expectation values.  Use of this class can drastically	**
**  reduce computation time when multiple expectation values	**
**  of the type <Op(t)> = Trace<Op*sigma(t)> are desired for	**
**  over evenly spaced time increments.				**
**								**
**  Class acquire is designed around the following reasoning:	**
**								**
**  If one calculates the expectation value			**
**						     -1		**
**  <Op(t )> = Tr<Op*sigma(t )> = Tr<Op*U*sigma(t )*U		**
**       1                  1                    o		**
**								**
**  where t = t  + t , then the values at t  = t + nt  can be	**
**	   1    o    d	                   n    o    d		**
**								**
**  more readily computed than the first.			**
**					       n	    -n	**
**  <Op(t )> = Tr<Op*U*sigma(t   )*U> = Tr<Op*U *sigma(t )*U  >	**
**       n	              n-1                       o	**
**		   -n	  n					**
**  <Op(t )> = Tr<U  *Op*U *sigma(t )>				**
**       n	                   o				**
**				     n	    -n			**
**  Working in the eigenbasis of U, U  and U   are easily	**
**  determined.							**
**	       ---       -n	           n			**
**  <Op(t )> = \   Tr<a|U  |a><a|Op|b>*<b|U |b>*<b|sigma(t )|a>	**
**       n     /					  o	**
**             ---						**
**             a,b						**
**								**
**	       ---       n					**
**  <Op(t )> = \   A  * B   where A and B are now the vectors	**
**       n     /    ab   ab					**
**             ---						**
**             ab						**
**								**
**	     -1							**
**  B  = <a|U  |a><b|U|b> and A  = <a|Op|b><b|sigma(t )|a>	**
**   ab                        ab                    o		**
**								**
**  Note that the vector B accounts entirely for the prop-	**
**  agator U. It need be determined only once. Similary, the	**
**  vector is computed only once as well.  The expectation	**
**  at differing times is easily determined by taking powers	**
**  of the elements of B in a scalar product with A.		**
**								**
**  The working form of for class acquires comes from the	**
**  recognition that vector B will contain many zeros.  Using	**
**  the index p and ingnoring all elements of B which are 0,	**
**								**
**	                       ---      n			**
**                  <Op(t )> = \   A * B			**
**                       n     /    p   p			**
**                             ---				**
**                              p				**
**								**
**  Thus, class require contains the vector B the matrix Op	**
**  which will combine with a sigma to form the vector A, and	**
**  the indexing p which will ignore contributions from zeros.	**
**								**
**  An example of class acquires use would be in the simulation	**
**  of a multi-dimensional NMR spectrum.  In such a case, there	**
**  are many FIDs calculated using equally spaced increments	**
**  between FID points (the dwell time).  Furthermore, each	**
**  FID uses the same detection operator and propagator U for	**
**  the delay between FID points.  Only the initial density	**
**  matrix varies. Thus, class acquire can be used to generate	**
**  the vector B only once.  Then, at the start of each FID the	**
**  vector A is generated from an input density matrix.  The	**
**  desired FID is then calculated with a minimum of effort.	**
**								**
*****************************************************************/

#ifndef   acquire_cc_			// Is this file already included?
#  define acquire_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Basics/Gconstants.h>			// Include PI
#include <Deprecated/acquire.h>
#include <HSLib/GenOp.h>
#include <LSLib/SuperOp.h>
#include <Matrix/matrix.h>
#include <Matrix/row_vector.h>
#include <Basics/StringCut.h>		  // Inlude Gform

  const double cutoffx = 1.0e-24;		// Element cutoff magnitude
//  const double cutoffx = 1.0e-12;		// Element cutoff magnitude

// ----------------------------------------------------------------------
// ------------------------ PRIVATE FUNCTIONS ---------------------------
// ----------------------------------------------------------------------

// ______________________________________________________________________
//                    CLASS ACQUIRE ERROR HANDLING
// ______________________________________________________________________


void acquire::create(gen_op &det, gen_op &U, int debug)

	// Input	det   : Detection operator (in trace)
	//	 	U     : Propagation operator (in trace)
	// Output	this  : An acquire, parameters to facilitate
	//			repetetive acquisition calculations
	// Note		      : Hilbert acquire - ls=0, J!=0, Sm1=0
	// Note		      : Setting ls=0 flags its Hilbert space

  {
debug=0;
//  if(debug)
//    {
//    cout << "\nDetection: " << det;
//    cout << "\nPropagator: " << U;
//    }
  LS = 0;				// Set the Liouville space to zero
  trinf=0;				// Set infinite trace to zero
  U.set_EBR();				// Insure U in its eigenbase
  det.Op_base(U);			// Insure det in eigenbase of U
  int hs = det.size();			// Get Hilbert space size
  int *I2Dtmp;			// Temporary index storage
  I2Dtmp = new int[hs*hs];
  int k=0;
  int i,j;
			// Determine Minimum Size
  
  pos=0;
  for(i=0; i<hs; i++)			// Loop over all det elements
    for(j=0; j<hs; j++)			// Only take non-zeros
      {
      if(square_norm(det(i,j))>cutoffx)
        {
	I2Dtmp[k] = 1;			// Temporary flag good element 
	pos++;				// Counter for good elements
        }
      else
	I2Dtmp[k] = 0;
      k++;
      }
			// Allocate Storage

  I = new int[pos];			// Minimum i index array
  J = new int[pos];			// Minimum j index array
  A = new complex[pos];			// Minimum A array (detector)
  B = new complex[pos];			// Minimum B array (propagator)
  bs = U.get_basis();			// Propagator basis

			// Compute A,B,I,J
  k = 0;
  int p = 0;
  for(i=0; i<hs; i++)
    for(j=0; j<hs; j++)
      {
      if(I2Dtmp[k])			// Only use good elements
	{
	I[p] = i;			// Store i index
	J[p] = j;			// Store j index
	A[p] = det(i,j);		// Store A & B for (i,j)
	B[p] = conj(U(i,i),U(j,j));
	p++;
	}
      k++;
      }
  delete [] I2Dtmp;
  }


void acquire::create(gen_op& det, super_op& eLt, double cutoff)

	// Input	det   : Detection operator (in trace)
	//	 	eLt   : Propagation superoperator (in trace)
	// Output	this  : An acquire, parameters to facilitate
	//			repetetive acquisition calculations
	// Note		      : eLt must be a superoperator which
	//			operates directly on the density matrix
	//			for the incrementation time
	// Note		      : Liouville acquire - ls!=0, J=NULL, Sm1!=0

  {
  trinf=0;				// Set infinite trace to zero
  LS = eLt.size();			// Set the Liouville space dimension
  eLt.set_EBR();			// Relaxation propagator into eigenbasis
  eLt.LOp_base(det);			// Put det into LOp Hilbert basis
//  basis S = eLt.get_Lbasis();		// Get superoperator basis
//  Sm1 = inv(S);				// Set superoperator basis inverse
  basis S = eLt.get_Lbasis();		// Get superoperator basis
  Sm1 = inv(S.U());				// Set superoperator basis inverse

  matrix mx = adjoint(det.get_mx());	// Start the formation of the detection ket
  mx = mx.resize(LS,1);			//          t                        t 
  mx = transpose(mx);			// Form <det | from transpose of |det >
					//          t            t
  matrix detS = mx*S.U();			// Form <det  * S| = <det |S
//  matrix detS = mx*S;			// Form <det  * S| = <det |S

  int *I2Dtmp;			// Temporary index storage
  I2Dtmp = new int[LS];
//  int k=0;
			// Determine Minimum Size

  pos=0;			
  int i; 				//                     t	
  for(i=0; i<LS; i++)			// Loop over all <1|det * S|i> elements
    {
    if(square_norm(detS.get(0,i))>cutoff)
      {
      I2Dtmp[i] = 1;			// Temporary flag good element 
      pos++;				// Counter for good elements
      }
    else
      I2Dtmp[i] = 0;
    }
			// Allocate Storage

  I = new int[pos];			// Minimum i index array
  J = new int[pos];			// Unused, consistent w/ Hilbert treatment
  A = new complex[pos];			// Minimum A array (detector)
  B = new complex[pos];			// Minimum B array (propagator)
  bs = eLt.get_basis();			// Propagator Hilbert space basis

			// Compute A,B,I
//  k = 0;
  int p = 0;
  for(i=0; i<LS; i++)
    if(I2Dtmp[i])			// Only use good elements
      {
      I[p] = i;				// Store "good" i index
      A[p] = detS.get(0,i);		// Store A vector
      B[p] = eLt.get(i,i);		// Store B vector
      p++;
      }
  delete [] I2Dtmp;
  }


void acquire::create(gen_op& det, super_op& eLt, gen_op sigma, int debug)

	// Input	det   : Detection operator (in trace)
	//	 	eLt   : Propagation superoperator (in trace)
	//		sigma : Infinite time density operator
	// Output	this  : An acquire, parameters to facilitate
	//			repetetive acquisition calculations
	// Note		      : eLt must be a superoperator which
	//			operates directly on the density matrix
	//			for the incrementation time
	// Note		      : Liouville acquire - ls!=0, J=NULL, Sm1!=0

  {
  LS = eLt.size();			// Set the Liouville space dimension
  eLt.set_EBR();			// Relaxation propagator into eigenbasis
  eLt.LOp_base(det);			// Put det into LOp Hilbert basis
  eLt.LOp_base(sigma);			// Put sigma into LOp Hilbert basis
  basis S = eLt.get_Lbasis();		// Get superoperator basis
  Sm1 = inv(S.U());		// Set superoperator basis inverse
  siginf = sigma.get_mx();		// Set infinite time matrix
  trinf = trace(det,sigma);		// Set infinite time trace

  trinf=0;				// Set infinite trace to zero

  matrix mx = adjoint(det.get_mx());	// Start the formation of the detection ket
  mx = mx.resize(LS,1);			//          t                        t 
  mx = transpose(mx);			// Form <det | from transpose of |det >
					//          t            t
  matrix detS = mx*S.U();			// Form <det  * S| = <det |S
  if(debug)
    std::cout << "\n\n\tdetS from acquire" << detS;

  int *I2Dtmp;			// Temporary index storage
  I2Dtmp = new int[LS];
  int k=0;
			// Determine Minimum Size

  pos=0;
  int i;				//                     t	
  for(i=0; i<LS; i++)			// Loop over all <1|det * S|i> elements
    {
    if(square_norm(detS.get(0,i))>cutoffx)
      {
      I2Dtmp[i] = 1;			// Temporary flag good element 
      pos++;				// Counter for good elements
      }
    else
      I2Dtmp[i] = 0;
    }
			// Allocate Storage

  I = new int[pos];			// Minimum i index array
  J = new int[pos];			// Unused, consistent w/ Hilbert treatment
  A = new complex[pos];			// Minimum A array (detector)
  B = new complex[pos];			// Minimum B array (propagator)
  bs = eLt.get_basis();			// Propagator Hilbert space basis

			// Compute A,B,I
  k = 0;
  int p = 0;
  for(i=0; i<LS; i++)
    if(I2Dtmp[i])			// Only use good elements
      {
      I[p] = i;				// Store "good" i index
      A[p] = detS(0,i);			// Store A vector
      B[p] = eLt.get(i,i);		// Store B vector
      p++;
      }
  delete [] I2Dtmp;
  }


// ----------------------------------------------------------------------
// ------------------------- PUBLIC FUNCTIONS ---------------------------
// ----------------------------------------------------------------------

// ______________________________________________________________________
//                CLASS ACQUIRE CONSTRUCTION, DESTRUCTION
// ______________________________________________________________________

// ------------------------- Null Constructor ---------------------------

acquire::acquire() : A(NULL), B(NULL), I(NULL)

	// Input	none  :
	// Output	ACQ   : A NULL acquire (this)

  { LS = 0; }				// Insure Liouville space is 0


// --------------- Hilbert Space Treatment Constructors -----------------


acquire::acquire(gen_op &det, gen_op &U)

	// Input	det   : Detection operator
	//		U     : A propagator for time increment
	// Output	ACQ   : Acquire (this) is constructed

  { create(det, U); }			// Fill ACQ 


acquire::acquire(matrix &det, gen_op &U)

	// Input	det   : Detection matrix
	//		U     : A propagator for time increment
	// Output	ACQ   : Acquire (this) is constructed

  {
  gen_op Op(det);			// Set det to an operator 
  create(Op, U);			// Fill ACQ 
  }


acquire::acquire(gen_op &det, gen_op &H, double dw)

	// Input	det   : Detection operator (trace computation)
	//	 	H     : Evolution Hamiltonian (trace computation)
	//		dw    : Dwell time between acquisition points
	// None		ACQ   : Acquire (this) is constructed

  {
  H.set_EBR();				// Hamiltonian into its eigenbase
  complex z(0,-2.0*PI*dw);		// Exponential factor (Hz to rad/sec)
  gen_op U = H.exp(z);			// Compute evolution propagator
//  gen_op U = exp(z*H);			// Compute evolution propagator
  create(det, U);			// Fill up the ACQ
  }


acquire::acquire(matrix &det, gen_op &H, double dw)

	// Input	det   : Detection matrix (trace computation)
	//	 	H     : Evolution Hamiltonian (trace computation)
	//		dw    : Dwell time between acquisition points
	// None		ACQ   : Acquire (this) is constructed
  {
  H.set_EBR();				// Hamiltonian into its eigenbase
  complex z = complex(0,-2.0*PI*dw);
  gen_op U = H.exp(z);			// Compute evolution propagator
//  gen_op U =				// Compute evolution propagator
//         exp( complex(0,-2.0*PI*dw)*H );
  gen_op Op(det);			// Set det matrix to an operator
  create(Op, U);			// Fill ACQ
  }


// -------------- Liouville Space Treatment Constructors ----------------


acquire::acquire(gen_op &det, super_op &L)

	// Input	det   : Detection operator
	//		L     : A superoperator for time increment
	// Output	ACQ   : Acquire (this) is constructed

// sosi - this is a VERY temporary patch
  { create(det, L, 1.e-6);  }		// Fill ACQ 


acquire::acquire(matrix &det, super_op &L)

	// Input	det   : Detection matrix
	//		L     : A superoperator for time increment
	// Output	ACQ   : Acquire (this) is constructed

  {
  gen_op Op(det);			// Set det to an operator 
  create(Op, L);			// Fill ACQ 
  }


acquire::acquire(gen_op &det, super_op &R, double dw)

	// Input	det   : Detection operator (trace computation)
	//	 	R     : Evolution Superoperator (trace computation)
	//		dw    : Dwell time between acquisition points
	// None		ACQ   : Acquire (this) is constructed

  {
  R.set_EBR();				// Hamiltonian into its eigenbase
  super_op L = exp(R, -dw); 		// Compute evolution propagator
  create(det, L);			// Fill ACQ
  }


acquire::acquire(matrix &det, super_op &R, double dw)

	// Input	det   : Detection matrix (trace computation)
	//	 	R     : Evolution Superoperator (trace computation)
	//		dw    : Dwell time between acquisition points
	// None		ACQ   : Acquire (this) is constructed
  {
  R.set_EBR();				// Hamiltonian into its eigenbase
  super_op L = exp(R, -dw); 		// Compute evolution propagator
  gen_op Op(det);			// Set det to an operator 
  create(Op, L);			// Fill ACQ
  }


// ------------------ Redfield Treatment Constructors -------------------


acquire::acquire(gen_op &det, super_op &L, gen_op& sigma)

	// Input	det   : Detection operator
	//		L     : A superoperator for time increment
	//		sigma : Infinite time density matrix
	// Output	ACQ   : Acquire (this) is constructed

  { create(det, L, sigma); }		// Fill ACQ 


acquire::acquire(matrix &det, super_op &L, gen_op& sigma)

	// Input	det   : Detection matrix
	//		sigma : Infinite time density matrix
	//		L     : A superoperator for time increment
	// Output	ACQ   : Acquire (this) is constructed

  {
  gen_op Op(det);			// Set det to an operator 
  create(Op, L, sigma);			// Fill ACQ 
  }


acquire::acquire(gen_op &det, super_op &R, gen_op& sigma, double dw)

	// Input	det   : Detection operator (trace computation)
	//		sigma : Infinite time density matrix
	//	 	R     : Evolution Superoperator (trace computation)
	//		dw    : Dwell time between acquisition points
	// None		ACQ   : Acquire (this) is constructed
  {
  R.set_EBR();				// Hamiltonian into its eigenbase
  super_op L = exp(R, -dw); 		// Compute evolution propagator
  create(det, L, sigma);		// Fill ACQ
  }


acquire::acquire(matrix &det, super_op &R, gen_op& sigma, double dw)

	// Input	det   : Detection matrix (trace computation)
	//		sigma : Infinite time density matrix
	//	 	R     : Evolution Superoperator (trace computation)
	//		dw    : Dwell time between acquisition points
	// None		ACQ   : Acquire (this) is constructed
  {
  R.set_EBR();				// Hamiltonian into its eigenbase
  super_op L = exp(R, -dw); 		// Compute evolution propagator
  gen_op Op(det);			// Set det to an operator 
  create(Op, L, sigma);			// Fill ACQ
  }


// ------------------------ Self Construction ---------------------------


// sosi: Memeory allocation and copying avoidable with class integer matrix!
acquire::acquire (const acquire& ACQ1)

	// Input	ACQ1  : Acquire
	// None		ACQ   : Acquire (this) constructed from ACQ1

  {
  LS = ACQ1.LS;			// Copy the Liouville space size
  pos = ACQ1.pos;		// Copy number of necessary points
  I = new int[pos];		// Minimum i index array
  J = new int[pos];		// Minimum j index array
  A = new complex[pos];		// Minimum A array (detector)
  B = new complex[pos];		// Minimum B array (propagator)
  for(int i=0; i<pos; i++)
    {
    A[i] = ACQ1.A[i];		// Copy A array (detector equivalent)
    B[i] = ACQ1.B[i];		// Copy B array (time prop. equivalent)
    I[i] = ACQ1.I[i]; 		// Copy array for index i
    J[i] = ACQ1.J[i];		// Copy array for index j
    }
  bs = ACQ1.bs;			// Copy basis (eigenbasis of prop)
  Sm1 = ACQ1.Sm1;		// Copy the Liouville inverse eigenvectors
  siginf = ACQ1.siginf;		// Copy the infinite time density matrix
  trinf = ACQ1.trinf;		// Copy the infinite time trace
  }

// --------------------------- Destruction ------------------------------

acquire::~acquire()

	// Input	ACQ   : An acquire (this)
	// Output	none  : ACQ is destructed

  { 
  if (A) delete [] A;
  if (B) delete [] B;
  if (I) delete [] I;
  if (J) delete [] J;
  }

// ---------------------------- Assignment ------------------------------


void acquire::operator = (const acquire& ACQ1)

	// Input	ACQ1  : Acquire
	// None		ACQ   : Acquire (this) copied from ACQ1
  {
  if(pos)			// Check for filled ACQ, destruct old
    {
    if (A) delete [] A;
    if (B) delete [] B;
    if (I) delete [] I;
    if (J) delete [] J;
    }
  LS = ACQ1.LS;			// Copy the Liouville dimension
  pos = ACQ1.pos;		// Copy number of necessary points
// sosi: Memeory allocation and copying avoidable with class integer matrix!
  I = new int[pos];		// Minimum i index array
  J = new int[pos];		// Minimum j index array
  A = new complex[pos];		// Minimum A array (detector)
  B = new complex[pos];		// Minimum B array (propagator)
  for(int i=0; i<pos; i++)
    {
    A[i] = ACQ1.A[i];		// Copy A array (detector equivalent)
    B[i] = ACQ1.B[i];		// Copy B array (time prop. equivalent)
    I[i] = ACQ1.I[i]; 		// Copy array for index i
    J[i] = ACQ1.J[i];		// Copy array for index j
    }
  bs = ACQ1.bs;			// Copy Hilbert basis
  Sm1 = ACQ1.Sm1;		// Copy the Liouville inverse eigenvectors
  siginf = ACQ1.siginf;		// Copy the infinite time density matrix
  trinf = ACQ1.trinf;		// Copy the infinite time trace
  }

// ______________________________________________________________________
//                    EXPECTATION VALUE CALCULATIONS
// ______________________________________________________________________


void acquire::operator () (gen_op &sigma, row_vector& fid, int np)

	// Input	ACQ   : An acquire (this)
	// 		sigma : Density matrix (operator propagated)
	//		fid   :	Data row_vector containing the result
	//		np    : Number of points to compute
	// Output	None  : Data row_vector "fid" filled
 
// This routine fills the data block fid according to the formula
 
//               t                    -1              t                  n  -1
//  fid(n) = <det *S|exp(L    *n*td)|S  *simga> = <det *S|[exp(L    *td)] |S  *sigma>
//                        diag                                  diag
 
// where det is a detection operator, sigma the state of the system at the
// start of some series of evolution steps by time increment td, and the
// exponential Liouvillian L is
 
//                                                    -1  
//                         exp(L*t) = S*exp(L    *t)*S
//                                           diag
 
// In terms of class acquire, this is then
 
//                         n                                th              -1
//  fid(n) = Sum { A[i] * B [i] * C[i] } where C[i] is the i   element of |S  * sigma>
//
 
// Variation 1: Hilbert space treatment - L is already diagonal, S = inv(S) = I
// Variation 2: Redfield treatment      - It is not sigma which evolved but sigma-sigmainf
//					  thus component from this must be added and subtracted.
//					  (Strictly true only if <det|exp(L*t)|sigmainf> !=0)
 
  {
  if(np<1) np = fid.size(); 		// Get number of points if <=0
  sigma.Op_base(bs);			// Put sigma in proper basis
  matrix sig = sigma.get_mx();
  matrix Sm1sig;	
  int k, pt = 0;
  complex *A1 = new complex[pos];
  if(LS==0)
    {					// Hilbert space acquisition
    for(k=0; k<pos; k++)		// Generate full Hilbert A array
      A1[k] = A[k]*sigma(J[k],I[k]);
    }
  else if(!siginf.rows())		// Liouville space acquisition
    {					//	 -1
    Sm1sig = Sm1*sig.resize(LS,1);	// Set |S  * sigma>
    for(k=0; k<pos; k++)		// Generate full Liouville A array
      A1[k] = A[k]*Sm1sig(I[k],0);
    }
  else					// Redfield acquisition
    {
    sig -= siginf;			//	 -1
    Sm1sig = Sm1*sig.resize(LS,1);	// Set |S  * {sigma-sigmainf}>
    for(k=0; k<pos; k++)		// Generate full Liouville A array
      A1[k] = A[k]*Sm1sig(I[k],0);
    }
  for(k=0; k<np; k++)			// Compute acquisition from reduced
    {					// arrays A1 and B.  Since B contains
    complex z(0);			// evolution over the dwell time, by
    for(pt=0; pt<pos; pt++)		// re-applying B (by looping over k)
      {					// we evolve to each successive point.
      z += A1[pt];
      A1[pt] *= B[pt];
      }
    fid.put(z, k);
    } 
  if(LS && siginf.rows())		// For Redfield acquisiton add in the
    if(square_norm(trinf)>cutoffx) 	// infinite trace component if it is
      for(k=0; k<np; k++)		// significant
        fid.put(fid.get(k)+trinf, k);
  delete [] A1;
  }
 

void acquire::operator () (gen_op &sigma, matrix& mx, double cutoff)

	// Input	ACQ   : An acquire (this)
	// 		sigma : Density matrix (operator propagated)
        //              mx    : A matrix for output data
	//		cutoff: An intensity cutoff, contributions with
	//			norms below this value are not output
        // Output       void  : The matrix mx is output as a posx2 array.
	//			Each row corresponds to one of the pos components
	//			detected according to det, evolving according to
	//			exp(Lt), and beginning in a state described by sigma.
	//			Component i is a complex exponential oscillating
	//			at frequency Im(mx(i,0)) & decaying at rate Re(mx(i,0))
	//			with complex intensity mx(i,1)
        // Note               : No sorting is performed!  No components are combined!
        // Note               : If ACQ was formed with exp(L*td) with td != 1, then
	//			          THIS ROUTINE WILL BE INAPPROPRIATE!
        // Note               : An ACQ formed with (detect,H) then this routine
	//			will not do what it is meant to.  In Hilbert space,
	//			use constructor (detect,H,1) and then the eigenvalues
	//			in B will be the exponential transition frequencies
 
// This routine fills the matrix mx with components of the time evolved signal
// fid(t) as defined by the formula
 
//                                t                 -1
//                   fid(t) = <det *S|exp(L    *t)|S  *simga>
//                                         diag
 
// where det is a detection operator, sigma the state of the system at the start
// of some evolution, and the exponential Liouvillian L which describes the evolution is
 
//                                                    -1  
//                        exp(L*t) = S*exp(L    *t)*S
//                                          diag
 
// In terms of class acquire, this is then
 
//                   fid(t) = Sum { A[i] * B[i]*t * C[i] }
//                             i
//                    th              -1
// where C[i] is the i   element of |S  * sigma>.  The matrix mx contains a row
// for each signal component, i.e. for each i.  Each row, or each component, has
// two parts: column 0 = B[i], column 1 = A[i]*C[i].
 
//               mx(i,0) = B[i]           mx(i,1) = A[i] * C[i]

//  Since B[i]*t are complex exponentials, column 0 contains exponential precession
// (imaginaries) and decay (real) rates.  Column 1 contains the complex exponential
// intensities (magnitudes and phases).  

// Note that an "analog" frequency domain signal can be obtained by Fourier transformation
// of the above equations.  The result would be that the values in the matrix then
// correspond to complex Lorentzians with frequencies and witdhs in column 0 and the
// intensities in column 1.
 
// Variation 1: Hilbert space treatment - L is already diagonal, S = inv(S) = I
// Variation 2: Redfield treatment      - It is not sigma which evolved but sigma-sigmainf
//					  thus component from this must be added and subtracted.
//					  (Strictly true only if <det|exp(L*t)|sigmainf> !=0)
 
  {
  sigma.Op_base(bs);                    // Put sigma in proper basis
  matrix sig, Sm1sig;
  sig = sigma.get_mx();			// Get the sigma matrix (for resize)
  complex *A1 = new complex[pos];	// Array to combined A and C 
  int* I1 = new int[pos];		// Array to flag sizeable components
  int pos1 = 0;				// Number of sizeable components
  int i=0, j=0;
  if(LS==0)
    {					// Hilbert space acquisition
    for(i=0; i<pos; i++)		// Loop over all non-zero components
      {
      A1[i] = A[i]*sigma(J[i],I[i]);	// Store the complex magnitudes
      if(norm(A1[i]) >= cutoff)		// Test for magnitude vs. cutoff
        {				// Only pos1 of pos will be the
        I1[i] = 1;			// contributions large enough
        pos1++;
        }
      else
        I1[i] = 0;
      }
    }
  else 
    {
    Sm1sig = Sm1*sig.resize(LS,1);	// Set |S  * sigma>
    for(i=0; i<pos; i++)		// Loop over all non-zero components
      {
      A1[i] = A[i]*Sm1sig(I[i],0);	// Store the complex magnitudes
      if(norm(A1[i]) >= cutoff)		// Test for magnitude vs. cutoff
        {				// Only pos1 of pos will be the
        I1[i] = 1;			// contributions large enough
        pos1++;
        }
      else
        I1[i] = 0;
      }
    }
  mx = matrix(pos1,2);			// Allocate array for these components
  for(i=0, j=0; i<pos; i++)		// Loop over all non-zero components
    {					// Index j will track good components
    if(I1[i])				// This flags what is good vs. bad
      {
      mx.put(B[i],j,0);			// Store precession & decay rates
      mx.put(A1[i], j,1);			// Store the complex magnitudes
      j++;				// Increment to good component count
      }
    }
  delete [] A1;
  delete [] I1;
  return;
  }

  
// ______________________________________________________________________
//                  CLASS ACQUIRE AUXILIARY FUNCTIONS
// ______________________________________________________________________


int acquire::ls() const { return LS; }

	// Input	ACQ   : An acquire (this)
	// Output	size  : Size of default Liouville space


int acquire::size() { return pos; }

	// Input	ACQ   : An acquire (this)
	// Output	size  : Size of the arrays A & B


int acquire::full_size()

	// Input	ACQ   : An acquire (this)
	// Output	size  : Full size of operators
	//			the Hilbert space squared

  {
  int HS = bs.size();
  return HS*HS;
  }


// ____________________________________________________________________________
//                        CLASS ACQUIRE I/O FUNCTIONS
// ____________________________________________________________________________

// ----------------- ASCII Output Functions (For Debugging) -------------------


std::ostream& acquire::print(std::ostream &ostr)

	// Input		ACQ   : Acquire
        //                      ostr  : Output stream
        // Output               none  : Acquire vector is sent
	//				to the output stream

  {
  ostr << "\n" << pos << " Non-zero"			// Output Points
       << " Points out of " << full_size() << " Possible.";
  ostr << "\nA,B Pairs: ";				// Output Vectors
  int i=0;
  if(LS == 0)
    {
    for(i=0; i<pos; i++)
      ostr << "\n" << i << ". A = <"
           << I[i] << "|Op|" << J[i]
           << "> : " << A[i] << "; B = <"
           << I[i] << "|U|" << I[i]
           << "><" << J[i] << "|U*|"
           << J[i] << "> : " << B[i];
    }
  else
    {
    for(i=0; i<pos; i++)
      ostr << "\n" << i << ". A = <1|Opt*S|"
           << I[i] << "> : " << A[i]
           << "; B = <" << I[i] << "|D|"
           << I[i] << "> : " << B[i];
    }
  ostr << "\nHilbert Space Basis: " << bs;	// Output Basis
  if((siginf).rows())
    {
    ostr << "\nInfinite Time Density Matrix:" << siginf;
    ostr << "\nInfinite Time Matrix Trace:" << trinf;
    }
  return ostr;
  }
  

std::ostream& acquire::print(std::ostream &ostr, gen_op& sigma)

	// Input		ACQ   : Acquire
        //                      ostr  : Output stream
        // Output               none  : Acquire vector is sent
	//				to the output stream

  {
  sigma.Op_base(bs);				// Put sigma in proper basis
  matrix sig = sigma.get_mx();			// Get the sigma matrix (for resize)
  matrix Sm1sig; 				// Set |inv(S) * sigma>
  if(LS) Sm1sig = Sm1*sig.resize(LS,1);		// if Liouville space exists

  ostr << "\n" << pos << " Non-zero"		// Output Points
       << " Points out of " << full_size()
       << " Possible.";
  ostr << "\nA,B Pairs: ";			// Output Vectors
  int i=0;
  double nsum1=0, nsum2=0, nsum3=0;
  if(LS == 0)
    {
    for(i=0; i<pos; i++)
      {
      ostr << "\n" << i << ". A = <"
           << I[i] << "|Op|" << J[i]
           << "> : " << A[i] << "; B = <"
           << I[i] << "|U|" << I[i]
           << "><" << J[i] << "|U*|"
           << J[i] << "> : " << B[i]
           << ";C = <"
           << J[i] << "|sig|" << I[i]
           << ">: " << sigma(J[i],I[i])
           << "; AC = " << A[i]*sigma(J[i],I[i]);
      nsum1 += norm(A[i]);
      nsum2 += norm(sigma(J[i],I[i]));
      nsum3 += norm(A[i]*sigma(J[i],I[i]));
      }
    }
  else
    {
    for(i=0; i<pos; i++)
      {
      ostr << "\n" << i << ". A = <1|Opt*S|"
           << I[i] << "> : " << A[i]
           << "; B = <" << I[i] << "|D|"
           << I[i] << "> : " << B[i]
           << ";C = <" << ">: " << Sm1sig(I[i],0)
           << "; AC = " << A[i]*Sm1sig(I[i],0);
      nsum1 += norm(A[i]);
      nsum2 += norm(Sm1sig.get(I[i],0));
      nsum3 += norm(A[i]*Sm1sig.get(I[i],0));
      }
    }
  ostr << "\n" << "Total norm(A): "
       << nsum1 << "; Total norm(C): "
       << nsum2 << "; Total norm(AC): "
       << Gform("%16.10f", nsum3) << "\n";

  ostr << "\nHilbert Space Basis: " << bs;	// Output Basis
  if((siginf).rows())
    {
    ostr << "\nInfinite Time Density Matrix:" << siginf;
    ostr << "\nInfinite Time Matrix Trace:" << trinf;
    }
  return ostr;
  }


std::ostream &operator << (std::ostream &ostr, acquire &ACQ)

	// Input		ACQ   : Acquire
        //                      ostr  : Output stream
        // Output               none  : Acquire vector is sent
	//				to the output stream

  {
  ostr << "\n" << ACQ.pos << " Non-zero"	// Output Points
       << " Points out of " << ACQ.full_size()
       << " Possible.";
  ostr << "\nA,B Pairs: ";			// Output Vectors
  int i=0;
  if(ACQ.LS == 0)
    for(i=0; i<ACQ.pos; i++)
      ostr << "\n" << i << ". A = <"
           << ACQ.I[i] << "|Op|" << ACQ.J[i]
           << "> : " << ACQ.A[i] << "; B = <"
           << ACQ.I[i] << "|U|" << ACQ.I[i]
           << "><" << ACQ.J[i] << "|U*|"
           << ACQ.J[i] << "> : " << ACQ.B[i];
  else
    for(i=0; i<ACQ.pos; i++)
      ostr << "\n" << i << ". A = <1|Opt*S|"
           << ACQ.I[i] << "> : " << ACQ.A[i]
           << "; B = <" << ACQ.I[i] << "|D|"
           << ACQ.I[i] << "> : " << ACQ.B[i];
  ostr << "\nHilbert Space Basis: " << ACQ.bs;	// Output Basis
  if((ACQ.siginf).rows())
    {
    ostr << "\nInfinite Time Density Matrix:" << ACQ.siginf;
    ostr << "\nInfinite Time Matrix Trace:" << ACQ.trinf;
    }
  return ostr;
  }


std::ostream& acquire::printT(std::ostream &ostr, gen_op& sigma, int np)

	// Input		ACQ   : Acquire
        //                      ostr  : Output stream
	// 			sigma : Density matrix (operator propagated)
	//			np    : Number of points to compute
        // Output               ostr  : The output stream modified by
        //                                the acquisition time domain
        //                                calculation 
// sosi - I've left out the steady state magnetization for convenience herein

  {
  print(ostr, sigma);			// 1st print computation core
  sigma.Op_base(bs);			// Put sigma in proper basis
  matrix sig = sigma.get_mx();		// Need sigma in matrix form
  matrix Sm1sig;			// For final vector (Liouville)	
  int k, pt = 0;
  complex *A1 = new complex[pos];	// This will be the FID state 
  complex *exps;			// for exponential factors
  complex *A10;
  exps = new complex[pos];
  A10 = new complex[pos];
  if(LS==0)
    {					// Hilbert space acquisition
    for(k=0; k<pos; k++)		// Generate full Hilbert A array
      A1[k] = A[k]*sigma(J[k],I[k]);
    }
  else if(!siginf.rows())		// Liouville space acquisition
    {					//	 -1
    Sm1sig = Sm1*sig.resize(LS,1);	// Set |S  * sigma>
    for(k=0; k<pos; k++)		// Generate full Liouville A array
      A1[k] = A[k]*Sm1sig(I[k],0);
    }
  else					// Redfield acquisition
    {
    sig -= siginf;			//	 -1
    Sm1sig = Sm1*sig.resize(LS,1);	// Set |S  * {sigma-sigmainf}>
    for(k=0; k<pos; k++)		// Generate full Liouville A array
      A1[k] = A[k]*Sm1sig(I[k],0);
    }
  ostr << "\n\n\t\tPoint Calculations";	// Output point calculation title
  for(k=0; k<pos; k++) A10[k] = A1[k]; 	// Get base intensities
  double tpt=0; 			// Time at particular FID point
  complex z;				// Point & exponential factor
  for(pt=0; pt<pos; pt++)		//   Loop the transitions
    exps[pt] = complex1;
  for(k=0; k<np; k++)			// Loop the points & compute them
    { 					// from reduced A1 and B
    z = complex0;			//   Initialize point to zer0
    ostr << "\n\npt " << k		//   Output point index
         << ". time = " << tpt;

// Since B contains the evolution over the dwell time, by re-applying B
// (by looping over k) we evolve every transition to each successive point.
//  Each FID point is obtained by summing up transition contributions

    for(pt=0; pt<pos; pt++)		//   Loop the transitions
      {
      z += A1[pt];			//     Add trans. contribution
      ostr << "\n\ttr " << pt << ". I="	//     Output intensity, frequency,
           << A10[pt] << ", W=" << B[pt]//     & the point time from start
           << ", E=" << B[pt];
      ostr << ", e=" << exps[pt];	//     Output exponential factor
      ostr << ", c=" << A1[pt]		//     Output contribution & point
           << ", pt=" << z;
      A1[pt] *= B[pt];			//     Evolve transition to next time
      exps[pt] *= B[pt];		//     Evolve exponential to next time
      }
//  tpt += tdw;				//  Increment point time (not avail.)
    ostr << "\nPoint " << k 		//  Output point value
         << " Final Value (x10^10): "
         << z*1.e10;
    } 
  delete [] A1;
  delete [] A10;
  delete [] exps;
  return ostr;
  }
  
// ______________________________________________________________________
//                  CLASS ACQUIRE OUTDATED FUNCTIONS
// ______________________________________________________________________


void acquire::operator () (gen_op &sigma, gen_op &sigma0, super_op &R,
						 row_vector& fid, int np)

	// Input	ACQ   : An acquire (this)
	// 		sigma : Density matrix (operator propagated)
	// 		sigma0: Equilibrium matrix
	//		R     : Relaxation superoperator
	//		fid   :	Data row_vector containing the result
	//		np    : Number of points to compute
	// Output	None  : Acquisition "fid" data row_vector filled
	// Note		      : This routine should utilze a relaxation
	//			superoperator R which utilizes the
	//			secular approximation!
  
{
// sosi: must check that R shares a basis with the acquire
// sosi: must check that the dwell time is not zero (if acq formed from prop)
  if(!np)				// Get number of points if 0
    np = fid.size();
  sigma.Op_base(bs);			// Put sigma in proper acquisiton basis
  sigma0.Op_base(bs);			// Put sigma0 in proper acquisiton basis
  complex *A1 = new complex[pos];	// Set up array space for A1, A2
  complex* A2 = new complex[pos];
  int k;
  for(k=0; k<pos; k++)			// Copy array A
    A2[k] = A[k];
  gen_op sigmadel = sigma-sigma0;	// Difference density matrix
  gen_op sigmat;			// For decaying sigma
  double t2;				// For total t2 time
//sosi- this now sucks too
double tdw=1.0;
  super_op eRt;				// For exponential superoperator
  eRt = exp(R, -tdw);			// Compute the exponential
  for (int i=0; i<np; i++)		// Compute acquisition point i
    {
//sosi: this now sucks
double tdw = 1.0;
    t2 = -double(i)*tdw; 		// Compute t2 time
// sosi: Done faster if only the appropriate elements relax!!!
// sosi: Next 3 lines replaced by faster method (perhaps)
//    eRt = exp(R, t2);			// Compute the exponential
//    sigmat = eRt*sigmadel;		// Relax the difference density matrix
//    sigmat += sigma0;			// Add back in equilibruim
    sigmadel = eRt*sigmadel;		// Relax the difference density matrix
    sigmat = sigmadel + sigma0;		// Add back in equilibruim
    for (k=0; k<pos; k++)		// Generate an A1 array
      A1[k] = A2[k]*sigmat(J[k],I[k]);
    complex z(0);
    for (k=0; k<pos; k++)
      {
      z += A1[k];			// Generate point
      A2[k] *= B[k];			// Update A2 elements for tine
      }
    fid.put(z,i);			// Copy point into acquisition
    }
  delete [] A1;
  delete [] A2;
  }


void acquire::operator () (gen_op &sigma, gen_op &sigmass, super_op &R,
				row_vector& fid, double offset, int np)

	// Input	ACQ    : An acquire (this)
	// 		sigma  : Density matrix (operator propagated)
	// 		sigmass: Equilibrium matrix
	//		R      : Relaxation superoperator
	//		fid    : Data row_vector containing the result
	//		offset : RF field frequency offset
	//		np     : Number of points to compute
	// Output	None   : Acquisition "fid" data row_vector filled
	// Note		       : This routine should utilze a relaxation
	//		 	 superoperator R which utilizes the
	//		 	 secular approximation!
	// Note		       : This routine assumes relaxation is
	//			 performed in a rotating frame which
	//			 not that desired for acquisition.
  
{
// sosi: must check that R shares a basis with the acquire
// sosi: must check that the dwell time is not zero (if acq formed from prop)
  if(!np)				// Get number of points if 0
    np = fid.size();
  sigma.Op_base(bs);			// Put sigma in proper basis
  sigmass.Op_base(bs);			// Put sigmass in proper basis
  complex *A1 = new complex[pos];	// Set up array space for A1, A2
  complex* A2 = new complex[pos];
  int k;
  for(k=0; k<pos; k++)			// Copy array A
    A2[k] = A[k];
  gen_op sigmadel = sigma-sigmass;	// Difference density matrix
  gen_op sigmat;			// For decaying sigma
  double t2;				// For total t2 time
//sosi: this now sucks
double tdw=1.0;
super_op Rtmp;
  for (int i=0; i<np; i++)		// Compute acquisition point i
    {
    t2 = -double(i)*tdw; 		// Compute t2 time
// sosi: Done faster if only the appropriate elements relax!!!
// sosi: Isn't this way of getting the exponential exceedingly slow you idiot!!!
    Rtmp = R.exp(t2);
    sigmat = Rtmp*sigmadel;		// Relax the difference density matrix
//    sigmat = exp(R*t2)*sigmadel;	// Relax the difference density matrix
    sigmat += sigmass;			// Add back in equilibruim
// sosi: kick sigmat out of the interaction rep frame here
// sosi: kick sigmat out of the rotating frame here
// sosi: kick sigmat out of the odd basis here
    for (k=0; k<pos; k++)		// Generate an A1 array
      A1[k] = A2[k]*sigmat(J[k],I[k]);
    complex z(0);
    for (k=0; k<pos; k++)
      {
      z += A1[k];			// Generate point
      A2[k] *= B[k];			// Update A2 elements for tine
      }
    fid.put(z, i);			// Copy point into acquisition
    }
  delete [] A1;
  delete [] A2;
offset = 0;				// compiler likes this used
}


#endif		// acquire.cc
