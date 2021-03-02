/* acquire1D.cc *************************************************-*-c++-*-
**									**
**      	                G A M M A				**
**									**
**	  Class Acquire1D	                   Implementation	**
**						 			**
**	  Copyright (c) 1994				 		**
**	  Dr. Scott A. Smith			 			**
**        National High Magnetic Field Laboratory 			**
**        1800 E. Paul Dirac Drive                                	**
**        Tallahassee Florida, 32306-4005                         	**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** This file contains the implementation of class acquire1D which	**
** facilitates the repetitive computation of expectation values.	**
**									**
*************************************************************************/

#ifndef   acquire1D_cc_			// Is this file already included?
#  define acquire1D_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Basics/Gutils.h>		// Know GAMMA error messages
#include <Level2/acquire1D.h>		// Include the header
#include <Level2/TrnsTable1D.h>
#include <iostream>			// Include filestreams (read/write)
#include <Level1/Lorentzian.h>		// Include Lorentzians
#include <Level1/Exponential.h>		// Include Exponentials
#include <Matrix/matrix.h>		// Need to know arrays
#include <Matrix/row_vector.h>		// Need to know vectors 
#include <HSLib/GenOp.h>		// Need to know operators
#include <LSLib/SuperOp.h>		// Need to know superopertors
#include <Basics/StringCut.h>
#include <Basics/Gconstants.h>

using std::string;			// Using libstdc++ strings
using std::cout;			// Using libstdc++ standard output
using std::vector;			// Using libstdc++ STL vectors
using std::ostream;			// Using libstdc++ output streams
using std::ofstream;			// Using libstdc++ output file streams
using std::ifstream;			// Using libstdc++ input file streams
using std::ios;

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                 CLASS ACQUIRE1D ERROR HANDLING
// ____________________________________________________________________________

	// Input	ACQ1D : An acquire1D (this)
	//		eidx  : An error index
	//		noret : A flag for a carriage return
	// Output	void  : A error message is sent to std output

void acquire1D::ACQerror(int eidx, int noret) const
  {
  string hdr("1D Acquisition");
  string msg;
  switch (eidx)
    {
    case 1:  msg = string("Prepared Density Operator of Wrong Dimension");
             GAMMAerror(hdr,msg,noret); break; 				// (1)
    case 2:  msg = string("Acquisiton Attempt With NO Prepared DensOp Set");
             GAMMAerror(hdr,msg,noret); break; 				// (2)
    default: GAMMAerror(hdr,eidx,noret); break;
    }
  }

volatile void acquire1D::ACQfatal(int eidx) const
  {
  ACQerror(eidx, 1);
  if(eidx) ACQerror(0);
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                CLASS ACQUIRE1D PRIVATE CONSTRUCTORS
// ____________________________________________________________________________
 
/* These are functions that contain code commonly used by many of the routines
   in this class.  Since they deal with the base structure of the class they
   MUST be kept private and used in a fashion that will not corrupt the
   integrity of any 1D acquisition.  There are two types of these auxiliary
   construction functions. The first handles those which utilize either a 
   static Hamiltonian or static Liouvillian to evolve a system for an 
   unspecified time. The second handles those which utilize either a Hilbert
   space or Liouville space propagator defined for evolution over a specific
   time. In both cases the functions set up the "pre-acqusition", i.e. the base
   construction that maintains how the system will evolve should an aquisition
   take place.

	// Input	ACQ1D : An acquire1D (this)
	// Output	void  : Major components of ACQ1D are set
	//			  ls    the Liouville space
	//			  Sm1   The inverse eigenvector of L
	// Note		      : The following must be preset in the
	//			constructors:
	//			  L      - The system Liouvillian
	//			  det    - The detection superoperator
	//			  siginf - The infinite time density matrix
	//			           (it can be left NULL)
	//			  trinf  - The infinte time trace
	//			  cutoff - An intensity cutoff level         */

// ----------------- This When LOp Is The System Liouvillian ------------------

void acquire1D::create()
  {
  _LS = LOp.size();			// Set the Liouville space dimension
  LOp.set_EBR();			// Liouvillian into eigenbasis
  LOp.SetHSBaseOf(det);			// Put det into LOp Hilbert basis
  matrix S = (LOp.LBs()).U();		// Get superoperator basis
  Sm1 = inv(S);				// Set superoperator basis inverse
  col_vector mk = det.superket();
  mk = adjoint(S)*mk; 	                //               t         t
  row_vector detS = adjoint(mk);	// <det| * S = {S  * |det>} 

//        Determine Minimum Size Based on Detection Operator Only

  I.clear();				// Zero array of viable indices
  int i, k;				// Temp variables
  for(i=0; i<_LS; i++)			// Loop over <det * S|i> elements
    if(square_norm(detS(i)) > DCUTOFF)	// & use if they are large enough
      I.push_back(i);			//   Store the index i in I

//       Allocate Storage Based on the Reduced Size, Fill A & B 

  pos = I.size();			// Reduced detector dimension 
  A   = row_vector(pos);		// Minimum A array (detector)
  B = row_vector(pos);			// Minimum B array (propagator)
  for(i=0; i<pos; i++)			// Loop over Liouville space
    {
    k = I[i];				//   Viable index in full space
    A.put(detS(k),  i);			//   Store <A|k>   (partial intensity)
    B.put(LOp(k,k), i);			//   Store <k|B|k> (frequency, rate)
    }
  }

// ------------------ This When LOp Is A Superoperator Propagator -------------

void acquire1D::createU()
  {
  double MINR = 1.e-8;			// Minimum rate allowed
  _LS = LOp.size();			// Set the Liouville space dimension
  LOp.set_EBR();			// Liouvillian into eigenbasis
  LOp.SetHSBaseOf(det);			// Put det into LOp Hilbert basis
  matrix S = (LOp.LBs()).U();		// Get superoperator basis
  Sm1 = inv(S);				// Set superoperator basis inverse
  col_vector mk = det.superket();
  mk = adjoint(S)*mk; 	                //               t         t
  row_vector detS = adjoint(mk);	// <det| * S = {S  * |det>} 

//        Determine Minimum Size Based on Detection Operator Only

  I.clear();				// Zero array of viable indices
  int i, k;				// Temp variables
  for(i=0; i<_LS; i++)			// Loop over <det * S|i> elements
    if(square_norm(detS(i)) > DCUTOFF)	// & use if they are large enough
      I.push_back(i);			//   Store the index i in I

//       Allocate Storage Based on the Reduced Size, Fill A & B 

  pos = I.size();			// Reduced detector dimension 
  A   = row_vector(pos);		// Minimum A array (detector)
  B = row_vector(pos);			// Minimum B array (propagator)
  complex sf = 1.0/delt;		// To undo exponential values
  complex z;				// Value in exponent
  for(i=0; i<pos; i++)			// Loop over Liouville space
    {
    k = I[i];				//   Viable index in full space
    A.put(detS(k),  i);			//   Store <A|k>   (partial intensity)
    z = log(LOp(k,k))*sf;		//   Store <k|B|k> (frequency, rate)
    if(fabs(Re(z)) < MINR)
        z = complex(0,Im(z));
    B.put(z, i);			//   Store <k|B|k> (frequency, rate)
    }
  }


// ____________________________________________________________________________
// iii            CLASS ACQUIRE1D TRANSITION TABLE INTERFACE
// ____________________________________________________________________________
 
/* These are functions that contain code commonly used by many of the routines
   in this class.  Since they deal with the base structure of the class they
   MUST be kept private and used in a fashion that will not corrupt the
   integrity of any 1D acquisition.  There are two types of functions. First
   are those which set up the "pre-acqusition", i.e. the base construction that
   maintains how the system will evolve should an aquisiton take place.   The
   second type are "acqusition" functions that set up a transitions table when
   the user has specified the state of the system (prepared density operator)
   at the acqusition start.                                                  */




	// Input	ACQ1D	: An acquire1D (this)
	//		Sp	: Density matrix (operator propagated)
	// Void		TTab   	: The internal prepared density operator is
	//			  set to Sp and the transitions matrix TTab
	//			  is reconstructed if necessary

void acquire1D::make_table(const gen_op& Sp)
  {
  if(sigmap==Sp && TTab.cols())			// If nothing new & we know the
    return;					// trans. table then just return
  sigmap = Sp;					// Reset the prepared system
  int ls = sigmap.LS();				// Liouville space of sigmap
  if(ls !=  LOp.size()) ACQfatal(1);		// Insure sigmap proper size
  basis bs = LOp.Bs();				// Acquisition basis
  sigmap.Op_base(bs);				// Put sigmap in proper basis
  col_vector delsig = sigmap.superket(); 	// Get matrix for delsig
  if(siginf.dim())				// Subtract off siginf only if
    {						// there is anything there
    siginf.Op_base(bs);
    col_vector mk = siginf.superket();
    delsig -= mk;
    }
  matrix Sm1sig = Sm1*delsig;           	// Set |inv(S)* delsigma>
  complex *A1 = new complex[pos];		// Intensity vector replacement
  int k;					// Dummy index
  for(k=0; k<pos; k++)				// Generate new intensity vector
    A1[k] = A(k)*Sm1sig.get(I[k],0);
  int itr, ntr=0;				// # of non-zero transitions
  for(k=0; k<pos; k++)				// Find # of non-zero trans.
    if(square_norm(A1[k])>DCUTOFF) ntr++;
  if(square_norm(trinf)) ntr++;			// Add 0 transition if trinf !0
  matrix mx(ntr,2);				// Matrix for table
  for(k=0,itr=0; k<pos; k++)			// Fill table with the spectrum
    if(square_norm(A1[k])>DCUTOFF)		// Only if intensity > DCUTOFF
      {
      mx.put(B.get(k),itr,0);			//	w & R in column 0
      mx.put(A1[k],itr,1);			//	Intensity in column 1
      itr++;					//	Increase transit. index
      }
  if(square_norm(trinf))			// Add 0 transition if trinf !0
    {
    mx.put(complex0,itr,0);			//	w & R in column 0
    mx.put(trinf,itr, 1);			//	Intensity in column 1
    }
  TTab = TTable1D(mx);				// Store the new spectrum map
  delete [] A1;					// Destruct temp. vector herein
  return;
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                  CLASS ACQUIRE1D CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
// --------------------------- Simple Constructors ----------------------------
// ----------------------------------------------------------------------------

acquire1D::acquire1D()
  {
  _LS       = 0; 			// Insure Liouville space is 0
  pos       = 0;			// No acquire 1D dimension either
  trinf     = complex0;			// Set <adjoint(det)|siginf>=0 
  DCUTOFF   = 1.e-12;			// Set the detection cutoff level
  TTab.ICUT = 1.e-12;			// Set intensity cutoff (in T.Table)
  delt      = 0;			// No set time increment
  }

acquire1D::acquire1D(const acquire1D& ACQ1)
  {
  _LS     = ACQ1._LS;		// Copy the Liouville space size
  pos     = ACQ1.pos;		// Copy number of necessary points
  I       = ACQ1.I;		// Copy minimum I index array
  A       = ACQ1.A;		// Copy minimum A array (detector)
  B       = ACQ1.B;		// Copy minimum B array (propagator)
  DCUTOFF = ACQ1.DCUTOFF;	// Copy the detection cutoff value
  LOp     = ACQ1.LOp;		// Copy the system Liouvillian
  Sm1     = ACQ1.Sm1;		// Copy the Liouville inverse eigenvectors
  det     = ACQ1.det;		// Copy the system detection operator
  siginf  = ACQ1.siginf;	// Copy the infinite time density matrix
  trinf   = ACQ1.trinf;		// Copy the infinite time trace
  delt    = ACQ1.delt;		// Copy the time increment
  }

// ----------------------------------------------------------------------------
// ------------------ Hilbert Space Treatment Constructors --------------------
// ----------------------------------------------------------------------------

/* These constructors set up the 1D acquisition EXCLUDING the effects of
   relaxation.  That is, the system is to be evolved either under a static 
   Hamiltonian H or under a specified propagator in Hilbert space. In the
   latter case, U must be for a specific time increment.

	   Input	detect: Detection operator
	  		H     : The system Hamiltonian
	                U     : The evolution propagator (for td)
	  		cutoff: Detection cutoff level
	   Output	ACQ1D : Acquire1D (this) is constructed
	   Note		      : When U is used in construction we store
				a superoperator propagator as LOp. If H is
				used in construction we store the Liouvillian
				as LOp.
	   Note		      : In this case, siginf is NULL                */

acquire1D::acquire1D(gen_op& detect, gen_op& H)
  {
  H.set_EBR();					// Put H into its eigenbasis
  det     = detect;				// Set the detection operator
  LOp     = (-complexi)*Hsuper(H);		// Set Liouvillian (diagonal)
  trinf   = complex0;				// Set <adjoint(det)|siginf>=0 
  DCUTOFF = 1.e-12;				// Set detection cutoff level
  delt    = 0;					// Set no specific time inc.
  create();					// Set up the rest of ACQ1D
  }

acquire1D::acquire1D(gen_op& detect, gen_op& H, double cutoff)
  {
  H.set_EBR();					// Put H into its eigenbasis
  det     = detect;				// Set the detection operator
  LOp     = (-complexi)*Hsuper(H);		// Set Liouvillian (diagonal)
  trinf   = complex0;				// Set <adjoint(det)|siginf>=0 
  DCUTOFF = cutoff;				// Set detection cutoff level
  delt    = 0;					// Set no specific time inc.
  create();					// Set up the rest of ACQ1D
  }

acquire1D::acquire1D(gen_op& detect, HSprop& U)
  {
  U.SetEBR();					// Put U into its eigenbasis
  det     = detect;                             // Set the detection operator
// sosi - this is not correct because LOp is the Liouvillian whereas
//        this is a propagator. Must figure out adjusted class structure
//        and fix the below or dealings with LOp will be messed up!
//        yet this still seems to work, why??? Because createU fixes the
//        the differences but NOT when it comes to storage of LOp.
LOp = U_transform(U.Op());
  trinf   = complex0;                           // Set <adjoint(det)|siginf>=0
  DCUTOFF = 1.e-12;				// Set detection cutoff level
  delt = U.length();				// Set specific time increment
  createU(); 					// Set up the rest of ACQ1D
  }

acquire1D::acquire1D(gen_op& detect, HSprop& U, double cutoff)
  {
  U.SetEBR();					// Put U into its eigenbasis
  det     = detect;                             // Set the detection operator
// sosi - this is not correct because LOp is the Liouvillian whereas
//        this is a propagator. Must figure out adjusted class structure
//        and fix the below or dealings with LOp will be messed up!
//        yet this still seems to work, why??? Because createU fixes the
//        the differences but NOT when it comes to storage of LOp.
LOp = U_transform(U.Op());
  trinf   = complex0;                           // Set <adjoint(det)|siginf>=0
  DCUTOFF = cutoff;                             // Set detection cutoff level
  delt = U.length();				// Set specific time increment
  createU(); 					// Set up the rest of ACQ1D
  }

// ----------------------------------------------------------------------------
// ----------------- Liouville Space Treatment Constructors -------------------
// ----------------------------------------------------------------------------

/* These constructors set up the 1D acquisition with the possibility of
   INCLUDING the effects of relaxation.  That is, the system is evolved under
   a static Liouvillian superoperator in Liouville space.  Furthermore the
   acquistion tracks the "infinite time" density opertator which the system
   will ultimately evolve into (if that is of any consequence).

           Input        det   : Detection operator
                        L     : The system Liouvillian
                        siginf: Infinite time density matrix
	  		cut   : Intensity cutoff level
           Output       ACQ1D : Acquire1D (this) is constructed
	   Note		      : Although the trace calculation
	  			is basis independent, both det
	  			& siginf must share the same one
	  			for <adjoint(det)|siginf>. So, here
	  			it is explicitly set to that of LOp
	  			since it will be done later anyway.

    Here L should be the system Liouvillian. This is something akin to
    the following:

        1. L = -i [H, ] where H is the commutation supoeroperator
                        from say Hsuper(H) or -i*2*PI*commutator(H)
        2. L = -i [Heff, ] + R + X; 

    it can NOT be the unitary transformation superoperator formed from
    the Hamiltonian (L = U_transform(H)) because that is a propagator not
    a Liouvillian, so the math will get messed up.                           */

acquire1D::acquire1D(gen_op& detect, super_op& L, double cut)
  {
  det     = detect;				// Set the detection operator
  LOp     = L;					// Set the system Liouvillian
  trinf   = complex0;				// Set <adjoint(det)|siginf>=0 
  DCUTOFF = cut;				// Set detection cutoff level
  delt    = 0;					// Set no specific time inc.
  create();					// Setup rest of ACQ1D
  }

acquire1D::acquire1D(gen_op& detect, super_op& L, gen_op& sigmaf, double cut)
  { 
  L.SetHSBaseOf(sigmaf);			// Put sigmaf 2 LOp Hilb. bs
  L.SetHSBaseOf(detect);			// Put detect 2 LOp Hilb. bs
  det = detect;					// Set the detection operator
  LOp  = L;					// Set the system Liouvillian
  siginf = sigmaf;				// Store inf. time dens. op.
  trinf = trace(det,siginf);			// Set <adjoint(det)|siginf>=0 
  DCUTOFF = cut;				// Set the cutoff level
  delt    = 0;					// Set no specific time inc.
  create();					// Setup rest of ACQ1D
  }

acquire1D::acquire1D(gen_op& detect, LSprop& G, double cutoff)
  {
  G.SetEBR();					// Put U into its eigenbasis
  det     = detect;                             // Set the detection operator
LOp = G.LOp();
  trinf   = complex0;                           // Set <adjoint(det)|siginf>=0
  DCUTOFF = cutoff;                             // Set detection cutoff level
  delt = G.length();				// Set specific time increment
  createU(); 					// Set up the rest of ACQ1D
  }


// ----------------------------------------------------------------------------
// ------------------------- Destruction, Assignment --------------------------
// ----------------------------------------------------------------------------

acquire1D::~acquire1D() { }

acquire1D& acquire1D::operator= (const acquire1D& ACQ1)
  {
  if(this == &ACQ1) return *this;	// Self assignment does nothing
  _LS     = ACQ1._LS;			// Copy the Liouville dimension
  pos     = ACQ1.pos;			// Copy number of necessary points
  I       = ACQ1.I;			// Copy minimum i index array
  A       = ACQ1.A;			// Copy minimum A array (detector)
  B       = ACQ1.B;			// Copy minimum B array (propagator)
  DCUTOFF = ACQ1.DCUTOFF;		// Copy the detection cutoff value
  LOp     = ACQ1.LOp;			// Copy the system Liouvillian
//  bs      = ACQ1.bs;			// Copy Hilbert basis
  Sm1     = ACQ1.Sm1;			// Copy Liouville inverse eigenvectors
  det     = ACQ1.det;			// Copy detection operator
  siginf  = ACQ1.siginf;		// Copy infinite time density matrix
  trinf   = ACQ1.trinf;			// Copy infinite time trace
  delt    = ACQ1.delt;			// Copy specific time increment
  return *this;
  }


// ____________________________________________________________________________
// B                      ALTERATION AND ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions allow users to access the internals of the 1D acquisition.

   Function   Return                         Value Returned
   --------  --------   -------------------------------------------------------
      L      super_op   Liouvillian evolution superoperator
      D      gen_op     Detection operator
      S      matrix     Eigenvector array in Liouville space
      Sinv   matrix     Eigenvector inverse array in Liouville space
      TTable TTable1D   Transitions table                                    */

// sosi - add LSprop & G to obtain it?
// sosi - SetHSBaseOf function should be adjusted to use LOp or G...
const super_op& acquire1D::L()      const { return LOp; } 
const gen_op&   acquire1D::D()      const { return det; }
const matrix    acquire1D::S()      const { return (LOp.LBs()).U(); }
const matrix&   acquire1D::Sinv()   const { return Sm1; }
const TTable1D& acquire1D::TTable() const { return TTab; }

/* These functions allow users to reset some of the internals of a 1D 
   acquisition.  In doing so much of the stored acquisition paramters must
   be recomputed.

   Function   Return                         Value Returned
   --------  --------   -------------------------------------------------------
   Detector    void     New detection operator
      D      gen_op     Detection operator
      S      matrix     Eigenvector array in Liouville space
      Sinv   matrix     Eigenvector inverse array in Liouville space         */


// sosi - must adjust when using propagators....
//       probably should fix this ot use a general function as code repeats...
void acquire1D::Detector(const gen_op& detect)
  {
  det = detect;				// Set detector to detect
  LOp.SetHSBaseOf(det);			// Put det into LOp Hilbert basis
  matrix S = (LOp.LBs()).U();		// Get superoperator basis

  col_vector mk = det.superket();       // Start forming detection superbra
                                        //                    t        t
  matrix detS = adjoint(S)*mk;          // Form <det| * S = {S * |det>} 
  detS = adjoint(detS);

//             Determine Minimum Size Based on Detector Only

  int *I2Dtmp;			// Temporary index storage
  I2Dtmp = new int[_LS];
  pos=0;				//               t	
  int i;
  for(i=0; i<_LS; i++)		// Loop over <det * S|i> elements
    {
    if(square_norm(detS.get(0,i))>DCUTOFF)
      {
      I2Dtmp[i] = 1;			// Temporary flag good element 
      pos++;				// Counter for good elements
      }
    else I2Dtmp[i] = 0;
    }

//               Allocate Storage Based on the Size pos

  I = vector<int> (pos);		// Minimum i index array
  A = row_vector(pos);			// Minimum A array (detector)
  B = row_vector(pos);			// Minimum B array (propagator)

//                Compute Needed Vector Arrays A,B,I

  int p = 0;
  for(i=0; i<_LS; i++)			// Loop over Liouville space
    if(I2Dtmp[i])			// Only use "good" elements
      {
      I[p] = i;				// Store "good" i index
      A.put(detS.get(0,i),p);		// Store A vector (partial intensity)
      B.put(LOp.get(i,i),   p);		// Store B vector (frequency, rate)
      p++;				// Update dimension index counter
      }

  if(siginf.size())			// If inf. time dens. op.
    trinf = trace(det,siginf);		// Reset <adjoint(det)|siginf>=0
  delete [] I2Dtmp;
  return;
  }

// ____________________________________________________________________________
// C                    Time Domain Spectra Generation
// ____________________________________________________________________________

/* These functions take the 1D transitions table and generate a time domain
   spectrum.  The user just inputs the prepared system (density operator) that
   will be evolved during the spectrum acquisition, the number of points
   desired and the time increment to use between them.  There are a few other
   factors which will influence the spectra generated by these functions.

	ICUT	This is an intensity cutoff value.  Any transition whose
                intensity (norm) is below this value will not contribute to
                the output spectrum.
        PERCUT  This is a percentage cutoff value.  For each transition,
                points at times past where the transition intensity has
                decayed below PERCUT*Io will not be included.  The value of
		PERCUT is decimal %, i.e. PERCUT = [0, 1].                   */

	// Input	ACQ1D	: An acquire1D (this)
	// 		sigmap	: Density matrix (operator propagated)
	//		npts   	: Number of points desired
        // or
	//		data  : 1D row vector for time interferrogram
	//		tinc  : Time increment between points (sec)
	// Output	data  : 1D data vector is filled with a 
	//			time domain spectrum.
	// Note		      : The input time increment is adjusted to
	//			account for "Hz" output.
	// Note		      : Value PERCUT is a % in decimal form
	//			of the maximum exponential intensity
	//			considered worth adding to the spectrum.
	// Note		      : Each exponential has intensity above
	//			ICUT between the points [0-ifi]. Only
	//			those points used in generating the spectrum

// sosi - if tinc is set & using propagators should warn if that time is
//       not equal to the propagator time. Also should allow for these
//      function to not need a time if a propagator is used

row_vector acquire1D::T(const gen_op& sigmap, double tinc, int npts)
  { return T(sigmap, npts, tinc); }

row_vector acquire1D::T(const gen_op& sigmap, int npts,    double tinc)
  {
  make_table(sigmap);				// Construct transitions table
TTab.setFreqRev();				// Reverse the frequencies
  return TTab.T(npts,tinc);			// Use table to make spectrum
  }

void acquire1D::T(const gen_op& sigmap, row_vector& data, double tinc)
  {
  make_table(sigmap);				// Construct transitions table
TTab.setFreqRev();				// Reverse the frequencies
  TTab.T(data,tinc);				// Use table to make spectrum
  }

// ____________________________________________________________________________ 
// D                        Frequency Domain Spectra
// ____________________________________________________________________________

/* These functions take the 1D transitions table and generate a frequency
   domain spectrum.  The user just inputs the prepared system (density
   operator) that will be evolved during the spectrum acquisition, the 
   number of points desired and the spectral range. There are a few other
   factors which will influence the spectra generated by these functions.

	ICUT	This is an intensity cutoff value.  Any transition whose
                intensity (norm) is below this value will not contribute to
                the output spectrum.
        PERCUT  This is a percentage cutoff value.  For each transition,
                points at frequencies outside where the transition intensity
                has fallen below PERCUT*Io will not be included.  The value
		of PERCUT is decimal %, i.e. PERCUT = [0, 1]. 

	   Input	ACQ1D	: An acquire1D (this)
	   		sigmap	: Density matrix (operator propagated)
	  		npts   	: Number of points desired
	  		data	: 1D row vector for frequency spectrum
	  		Fst	: Starting frequency (Hz)
	  		Ffi	: Final frequency (Hz)
	   Output	data  	: Row vector filled with a frequency spectrum.
	   Note		      	: Value PERCUT is a % in decimal form
	  			  of the maximum Lorentzian intensity
	  			  considered worth adding to the spectrum.
	   Note		      	: Each Lorentzian (transition) will have an
	  			  intensity above ICUT*Imax between the points
	  			  bound by the indices ist and ifi.  Only those
	  			  points are used in generating the spectrum
	   Note		      	: To avoid poor resolution, if the frequency
	  			  between points, winc, does not fullfill the
	  			  relationship 5*winc < lwhh, then the
	  			  integrated Lorentzian intensities are used
	  		          rather than the point value.                */

row_vector acquire1D::F(const gen_op& sigmap,int npts,double Fst,double Ffi)
  {
  make_table(sigmap);				// Construct transitions table
  return TTab.F(npts,Fst,Ffi);			// Use table to make spectrum
  }

void acquire1D::F(const gen_op& sigmap,row_vector& data,double Fst,double Ffi)
  {
  make_table(sigmap);				// Construct transitions table
  TTab.F(data,Fst,Ffi);				// Use table to make spectrum
  }


// ____________________________________________________________________________ 
// E             Frequency Domain Spectra In Derivative Mode
// ____________________________________________________________________________

/* These functions take the 1D transitions table and generate a frequency
   domain spectrum in derivative mode (ESR/EPR).  The user just inputs the
   prepared system (density operator) that will be evolved during the spectrum
   acquisition, the number of points desired and the spectral range. There are
   a few other factors which will influence the spectra generated by these 
   functions.

	ICUT	This is an intensity cutoff value.  Any transition whose
                intensity (norm) is below this value will not contribute to
                the output spectrum.
        PERCUT  This is a percentage cutoff value.  For each transition,
                points at frequencies outside where the transition intensity
                has fallen below PERCUT*Io will not be included.  The value of
		PERCUT is decimal %, i.e. PERCUT = [0, 1].                   */

	// Input	ACQ1D	: An acquire1D (this)
	// 		sigmap	: Density matrix (operator propagated)
	//		npts   	: Number of points desired
	//		data	: 1D row vector for frequency spectrum
	//		Fst	: Starting frequency (Hz)
	//		Ffi	: Final frequency (Hz)
	// Output	data  	: Row vector filled with a frequency spectrum.
	// Note		      	: Value PERCUT is a % in decimal form
	//			  of the maximum Lorentzian intensity
	//			  considered worth adding to the spectrum.
	// Note		      	: Each Lorentzian (transition) will have an
	//			  intensity above ICUT*Imax between the points
 
row_vector acquire1D::FD(const gen_op& sigmap,int npts,double Fst,double Ffi)
  {
  make_table(sigmap);				// Construct transitions table
  return TTab.FD(npts,Fst,Ffi);			// Use table to make spectrum
  }

void acquire1D::FD(const gen_op& sigmap,row_vector& data,double Fst,double Ffi)

  {
  make_table(sigmap);				// Construct transitions table
  TTab.FD(data,Fst,Ffi);			// Use table to make spectrum
  }

// ____________________________________________________________________________ 
// F                     DIRECT TRANSITION TABLE GENERATION
// ____________________________________________________________________________

	// Input	ACQ1D : An acquire1D (this)
	// 		sigmap: Density matrix (operator propagated)
	// Output	TTab  : Transitions table
	//			Re(<i|TTab|0>) = transition i decay rate
	//			Im(<i|TTab|0>) = transition i frequency
	//			   <i|TTab|1>  = transition i intensity
	// Note		      : Frequencies and Rates are in 1/sec
	// Note		      : This sets TTab and sigmap internally

// sosi - trigger the table to know if propagator is used?
const TTable1D& acquire1D::table(const gen_op& sigmap)
  { make_table(sigmap); return TTab; }

const TTable1D& acquire1D::table() const { return TTab; }


const TTable1D acquire1D::table_snapshot(const gen_op& sigmap)
{
    make_table(sigmap); 
    return TTab; 
}

const TTable1D acquire1D::table_snapshot() const
{ 
   return TTab; 
}

// ____________________________________________________________________________ 
// G                     TRANSITION TABLE MANIPULATIONS
// ____________________________________________________________________________

/* These functions allow for the direct manipulation of the transitions table.
   By and large these are handled by the class TTable1D, so that the functions
   here simply relay to the class TTable1D functions using the internal table
   TTab.

   Function   Return                         Application
   --------  --------   -------------------------------------------------------
   offset      void     Transition(s) frequency offset by specified amount
   FRscale     void	Transition(s) frequency scaled by specified amount
   Iscale      void     Transition(s) intensity scaled by specified amount
   broaden     void     Transition(s) linewidth altered by specified amount
   resolution  void     Transitions blended if specified as unresolved 
   pcorrect   void/z    Transitions are phase corrected as specified         */

void acquire1D::offset(double F, int inHz)         { TTab.offset(F,inHz); }
void acquire1D::offset(double F, int tr, int inHz) { TTab.offset(F,tr,inHz); }
void acquire1D::FRscale(double Fscf)               { TTab.FRscale(Fscf); }
void acquire1D::FRscale(double Fscf, int tr)       { TTab.FRscale(Fscf, tr); }
void acquire1D::Iscale(double Iscf)                { TTab.Iscale(Iscf); }
void acquire1D::Iscale(double Iscf, int tr)        { TTab.Iscale(Iscf, tr); }
void acquire1D::broaden(double LWR, int inHz)      { TTab.broaden(LWR, inHz); }
void acquire1D::broaden(double LWR,int tr,int Hz)  { TTab.broaden(LWR,tr,Hz); }
void acquire1D::resolution(double res)             { TTab.resolution(res); }
void acquire1D::pcorrect(double Wpivot,complex& P) { TTab.pcorrect(Wpivot,P); }
complex acquire1D::pcorrect(double& w0,double w1,int ord)
                                           { return TTab.pcorrect(w0,w1,ord); }

/* These functions either return values regarding the transitions table or
   allow the user to affect the transitions table output format.  Again, these
   are handled by the class TTable1D, so that the functions here simply relay
   to the class TTable1D functions using the internal table TTab.

   Function   Return                         Application
   --------  --------   -------------------------------------------------------
   Wmax       double    Maximum transition frequency in the table
   LWmax      double    Maximum transition linewidth in the table
   setSort    void      Sets a flag for output sorting by frequency
   setConv    void      Sets a conversion factor for out frequencies         */

double acquire1D::Wmax()  const   { return TTab.FRmax(); }
double acquire1D::LWmax() const   { return TTab.LWmax(); }
void   acquire1D::setSort(int sf) { TTab.setSort(sf);    }
void   acquire1D::setConv(int cf) { TTab.setConv(cf);    }

/* The above function setConv sets the internal value of FRQCONV in class
   TrnsTable1D.  It is only active when printing the transitions table and then
   only when the frequencies are ouput in either PPM (FRQTYPE=1) or in Gauss
   (FRQTYPE=2) For the former users should set FRQCONV to be the Larmor 
   frequency in MHz as this value will be used for Hz -> PPM. For the latter
   set FRQCONV to be the electron g value                                    */

// ____________________________________________________________________________ 
// H                TRANSITION TABLE FORMATTED OUTPUT
// ____________________________________________________________________________

/* These functions allow users to output the transitions table to any output
   stream with nice formatting.  This ability is provided for directly in
   class TTable1D, so these functions here simply relay to the class TTable1D
   functions with the additional possibility of setting the initial spin
   system in the function as well.                                           */


// GAMMA's Lorentzian relationship between the linewidth at half-height (lwhh)
// and its R value (a decay rate of its related exponential) according to

//                    R     1
//          lwhh   = -- = -----  ---> lwhh  = 2R   &  lwhh     = 2R
//              Hz   PI   PI*T            Hz    Hz        1/sec    1/sec
//                            2
//                     -1
// where R & T  are sec  , lwhh in Hz. Switching both lwhh & R into the same
//            2                        units implies the relationship lw = 2R.

	// Input	ACQ1D	: An acquire1D (this)
	// 		sigmap	: Density matrix (operator propagated)
	// 		ostr	: An output stream
	// Output	TTab  : Matrix representation of the 1D spectrum
	//			Re(<i|TTab|0>) = transition i decay rate
	//			Im(<i|TTab|0>) = transition i frequency
	//			   <i|TTab|1>  = transition i intensity
	//			is placed into the output stream
	// Note		      : This sets TTab, sigmap internally

	// Input	ostr  : An output stream
	// Output       void  : A table of transition infomation
	//			is placed into the output stream
	//			Re(<i|TTab|0>) = transition i decay rate
	//			Im(<i|TTab|0>) = transition i frequency
	//			   <i|TTab|1>  = transition i intensity
	//		type  : Flag for printout type
	//				0 = Hertz (default)
	//				1 = PPM
	//				2 = Gauss
	//		cf    : Conversion factor if needed
	//			    type=1 ==> cf = Om in MHz
	//			    type=2 ==> cf = g (unitless)
        // Note               : Matrix freqs & rates are in 1/sec
	// Note		      : The array TTab should have been filled
	//			prior to this function call

//                                   R     1
//                         lwhh   = -- = -----
//                             Hz   PI   PI*T
//                                           2

void acquire1D::table(const gen_op& sigmap, ostream& ostr)
  {
  make_table(sigmap);
  TTab.print(ostr);
  }
 
void acquire1D::table(ostream& ostr) const
  {
  if(!TTab.rows())				// Check transitions exist
    { ACQerror(2); return; } 			// If not, warning & return
  TTab.print(ostr);
  }


// ____________________________________________________________________________
// F                  CLASS ACQUIRE1D AUXILIARY FUNCTIONS
// ____________________________________________________________________________


int acquire1D::ls()          const { return _LS; }
int acquire1D::size()        const { return pos; }
int acquire1D::full_size()   const { return LOp.size(); }
int acquire1D::transitions() const { return TTab.size(); }

// sosi - what does this now do?
// error as LW is redefined below
void acquire1D::parameters(matrix& mx, double& SW,
                    double& LW, double& dt, int N, int pf) const 

        // Input        ACQ1D : An acquire1D (this)
        //              mx    : Array of acquire parameters
        //              SW    : Suggested spectral width (Hz)
        //              LW    : Suggested linewidth (Hz)
        //              dt    : Suggested dwell time (sec)
        //              N     : Number of FID points
        //              pf    : Print flag
	// Output	void  : Full size of operators
	//			the Hilbert space squared
	// Note		X     : The final FID intensity used
	//			to determine required apodiczation is
	//			arbitrarily set to 10-3.  

  {
  double X = 1.e-5;				// Final FID intensity
//  double R=0, LW=0, W=0;
  double R=0, W=0;
  double LWmax = 0;
  double Wmax = 0;
  for(int i=0; i<mx.rows(); i++)		// Loop transitions
    {
    R = mx.getRe(i,0);				// Transition decay rate	
    LW = R/PI;					// Transition linewidth
    W = mx.getIm(i,0)/(2.0*PI);			// Transition freq (Hz)
    if(LW > LWmax) LWmax = LW;			// Maximum linewidth
    if(fabs(W) > Wmax) Wmax = fabs(W);		// Maximum frequency
    }
  SW = 2.0*Wmax;				// Suggested spectral width
  dt = 1.0/SW;					// Suggested dwell time
  LW = LWmax;					// Suggested linewidth
  double Rsugg = -log(X)/(dt*double(N-1));	// Required rate for X
  double LWsugg = Rsugg/PI;			// Required linewidth for X
  if(LWsugg > LW) LW = LWsugg; 			// Set suggested linewidth
  if(pf)
    {
    cout << "\n\tSuggested Spectral Width is " << SW << " Hz";
    cout << "\n\tSuggested Dwell Time is " << dt << " Seconds";
    if(Rsugg < R)
      cout << "\n\tApodization is Not Required";
    else
      cout << "\n\tSuggested Apodization is " << LW << " Hz";
    cout << "\n\t(Apodization Based on " << N 
         << " Point FID, Decay to " << X << " Intensity)";
    }
  }


// ____________________________________________________________________________ 
// G                 CLASS ACQUIRE1D ASCII BASED I/O FUNCTIONS
// ____________________________________________________________________________ 

// ----------------- ASCII Output Functions (For Debugging) -------------------

/*              Input           ACQ   : Acquire
                ostr                  : An output stream for info
                                sigmap: State of system at acquisition
                                tinc  : Dwell time in seconds
                                npts  : Number of points
                                         0 = Evolve From t=0 To Point
                                        !0 = Evolve From Point To Point
                Output          ostr  : The output stream modified by
                                        acquisition parameters/computation   
 
        Function                                Output
  --------------------------------------------------------------------------
    print(ostr)         Core elements of the acquisition computation
    print(ostr,sigmap)  Core elements of the acquisiton from state sigmap
    printT(.....)       Info on acquisiton time tomain calculation
    <<                  Same as print(ostr) function   
                                                                             */

// sosi - these must be adjusted for use of propagators...
ostream& acquire1D::print(ostream& ostr) const
  {
  string hdr;
  string spc;
  if(!_LS)
    {
    hdr = string("Empty 1D Acquisition");
    spc = string(40-hdr.length()/2, ' ');
    ostr << "\n" << spc << hdr << "\n";
    return ostr;
    }

  hdr = string("1D Acquisition");
  spc = string(40-hdr.length()/2, ' ');
  ostr << "\n" << spc << hdr;
  ostr << "\n" << spc << string(hdr.length(), '=');
  hdr = Gdec(pos) + string(" Transition");
  if(pos > 1) hdr += string("s");
  hdr += string(" From Liouville Dimension Of ");
  hdr += Gdec(full_size());
  spc = string(40-hdr.length()/2, ' ');
  ostr << "\n\n" << spc << hdr;

  if(!_LS)
    {
    ostr << "\n\nThe " << pos << " A,B Pairs";	// Output Vectors
    ostr << " (Hilbert Space Viewpoint)";
    ostr << ":\n";					
    }
  else 
    {
    ostr << "\n\n\t           t                -1                        -1";
    ostr << "\n\tTr(i) = <Op *S|i><i|D|i><i|S  *sigma> = A(i)*B(i)*<i|S  *sigma>";
    ostr << "\n";
    }
  string tmp = Gdec(pos);
  int ilen = tmp.length();
  tmp = Gdec(full_size());
  int lslen = tmp.length();
  int i=0;
  bool Areal = A.is_real();
  bool Bimag = B.is_imaginary();
  if(_LS == 0)					// This for a Hilbert
    { 						// space acquisition
    for(i=0; i<pos; i++)
      ostr << "\n\t" << Gdec(i,ilen) << ". A = <"
           << "> : " << A.get(i) << "; B = <"
           << I[i] << "|U|" << I[i]
           << "> : " << B.get(i);
    }
  else						// This for a Liouville
    { 						// space acquisition
    for(i=0; i<pos; i++)
      {
      ostr << "\n\t" << Gdec(i,ilen) << ". A = <Opt*S|"
           << Gdec(I[i], lslen) << "> : ";
           if(Areal) ostr << A.getRe(i);
           else      ostr << A.get(i);
      ostr << "    B = <" << Gdec(I[i], lslen) << "|D|"
           << Gdec(I[i], lslen) << "> : ";
           if(Bimag) ostr << B.getIm(i) << "i";
           else      ostr << B.get(i);
      }

    ostr << "\n\nHilbert Space Basis: "		// Output Basis
         << LOp.Bs();
    if((siginf).size())
      {
      ostr << "\nInfinite Time Density Operator:\n" << siginf;
      cout << "\nInfinite Time Matrix Trace:" << trinf;
      }
    }
  return ostr;
  }


ostream& acquire1D::print(ostream& ostr, gen_op& sigp)
  {
  if(!_LS)					// This if empty
    {
    string hdr("Empty 1D Acquisition\n");
    string spacer = string(40-hdr.length()/2, ' ');
    ostr << "\n\n" << spacer << hdr << "\n";
    return ostr;
    }

  string hdr = Gdec(pos) + " Contributing Transitions";
  string spacer = string(40-hdr.length()/2, ' ');
  ostr << "\n\n" << spacer << hdr;
  hdr = "(Liouville Space Is " + Gdec(full_size()) + ")";
  spacer = string(40-hdr.length()/2, ' ');
  ostr << "\n" << spacer << hdr << "\n";
// sosi - this does not possibly work! _LS=0 --> Liouville!
  if(!_LS) hdr = "Hilbert Space Viewpoint";
  else     hdr = "Liouville Space Viewpoint";
  spacer = string(40-hdr.length()/2, ' ');
  ostr << "\n\n" << spacer << hdr;
  ostr << "\n"   << spacer << string(hdr.length(), '=') << "\n";
  
  basis bs = LOp.Bs();				// Acquisition basis
  sigp.Op_base(bs);				// Put sigma in proper basis
  col_vector delsig=sigp.superket();		// Get matrix for delsig
  if(siginf.dim())				// Subtract off siginf only if
    {						// the infinite time density
    siginf.Op_base(bs);				// operator has been defined
    col_vector mk = siginf.superket();
    delsig -= mk;
    }
  matrix Sm1sig = Sm1*delsig;			// Set |inv(S)* delsigma>

  double nsum1=0, nsum2=0, nsum3=0;
  ostr << ":";					
  int i=0;
  if(_LS == 0)					// This for a Hilbert
    {
    for(i=0; i<pos; i++)			// space acquisition
      ostr << "\n\t" << i << ". A = <"
           << "> : " << A.get(i) << "; B = <"
           << I[i] << "|U|" << I[i]
           << "> : " << B.get(i) << "; C = <"
           << I[i] ;
    }
  else						// This for a Liouville
    { 						// space acquisition
    string kfmt, ifmt;				// First we set the i&k
    if(pos<10)        kfmt = "%1i";		// indices output format
    else if(pos<100)  kfmt = "%2i";		// so that we get nice
    else if(pos<1000) kfmt = "%3i";		// columns
    else              kfmt = "%4i";
    if(_LS<10)        ifmt = "%1i";
    else if(_LS<100)  ifmt = "%2i";
    else if(_LS<1000) ifmt = "%3i";
    else              ifmt = "%4i";
    for(i=0; i<pos; i++)
      {
      ostr << "\n" << Gform(kfmt.c_str(), i) << ". A=<Opt*S|"
           << Gform(ifmt.c_str(), I[i]) << ">: " << A.get(i)
           << "; B=<" << Gform(ifmt.c_str(), I[i]) << "|D|"
           << Gform(ifmt.c_str(), I[i]) << ">: " << B.get(i)
           << "; C=<" << Gform(ifmt.c_str(), I[i]) << "|Sm1*delsig>: "
           << Sm1sig.get(I[i],0)
           << "; AC=" << A.get(i)*Sm1sig.get(I[i],0);

      nsum1 += norm(A.get(i));
      nsum2 += norm(Sm1sig.get(I[i],0));
      nsum3 += norm(A.get(i)*Sm1sig.get(I[i],0));
      }
    ostr << "\n" << "Total norm(A): "
         << nsum1 << "; Total norm(C): "
         << nsum2 << "; Total norm(AC): "
         << Gform("%16.10f", nsum3);

    ostr << "\n\nHilbert Space Basis: "		// Output Basis
         << LOp.Bs();
    if((siginf).size())
      {
      ostr << "\nInfinite Time Density Operator:\n" << siginf;
      cout << "\nInfinite Time Matrix Trace:"       << trinf;
      }
    }
  return ostr;
  }

	// Input		ACQ1D	: Acquire1D (this)
	//			sigp	: State of system at detection
	//			tinc	: Dwell time in seconds
	//			npts	: Number of points
	//			P2P     : Computation control
	//				    0 = Evolve From t=0 To Point 
	//				   !0 = Evolve From Point To Point 
        // Output               ostr	: The output stream modified by
        //                                the acquisition time domain
	//				  calculation

ostream& acquire1D::printT(ostream& ostr, gen_op& sigp, double tinc, int npts, int P2P)
  {
  print(ostr, sigp);				// First use the overload
  int ntr = TTab.rows();			// Get number of transitions
  if(!ntr && !TTab.cols())			// Insure there are transitions
    {						// present or there can't be
    ACQerror(2); 				// any acquisition performed
    return ostr;
    }
  return TTab.printT(ostr, tinc, npts, P2P);
  }


ostream& operator << (ostream& ostr, acquire1D& ACQ)
  { ACQ.print(ostr); return ostr; }

// ____________________________________________________________________________
// H                 CLASS ACQUIRE1D BINARY BASED I/O FUNCTIONS
// ____________________________________________________________________________

// sosi - all binary I/O must be adjusted if using propagators

// ---------------------------- Binary Output ---------------------------------

	// Input		ACQ1D : Acquire1D (this)
        //                      fn    : Filename
        // Return               void  : Acquisition  ACQ1D is written
        //                              to a file called fn
        // Note                       : The file format is BINARY
 
	// Input		ACQ1D : Acquire1D (this)
        //                      fp    : File stream (pointing at ACQ1D)
        // Return               void  : Acquisition  ACQ1D is written
        //                              into file fp at current location
        // Note                       : The file format is BINARY and the file
	//				point fp is advanced
 

void acquire1D::write(const string& fn) const
  {
  ofstream fp;					// Construct a file
  fp.open(fn.c_str(),ios::out|ios::binary);	// Open file
  write(fp);                                    // Write Op
  fp.close();                                   // Close file
  }

ofstream& acquire1D::write(ofstream& fp) const
  {
  fp.write((char*)&_LS,sizeof(int));		// Write Liouville space dim.
  fp.write((char*)&pos,sizeof(int));		// Write ACQ1D dim. (pos <= ls)
  A.write(fp);					// Write ACQ1D complex array A
  B.write(fp);					// Write ACQ1D complex array B
  for(int i=0; i<pos; i++) {
    int val = I[i];
    fp.write((char*)&val,sizeof(int));		// Write ACQ1D int array I
  }
  fp.write((char*)&DCUTOFF,sizeof(double));	// Write ACQ1D cutoff
  LOp.write(fp);					// Write system Liouvillian
  Sm1.write(fp);			// Write inverse eigenverctor array
  det.write(fp);			// Write detection operator
  int inff;
  if(siginf.size())		
    {
    inff = 1;
    fp.write((char*)&inff,sizeof(int));
    siginf.write(fp);			// Write infinite time dens. op.
    trinf.write(fp); 			// Write infinite time Tr(det*siginf)
    }
  else
    {
    inff = 0;
    fp.write((char*)&inff,sizeof(int));
    }
  return fp;
  }
 
// ---------------------------- Binary Input ----------------------------------
 
	// Input		ACQ1D : Acquire1D (this)
        //                      fn    : Filename
        // Return               void  : ACQ1D is read from file fn
        // Note                       : No basis sharing is assumed!
 
	// Input		ACQ1D : Acquire1D (this)
        //                      fp    : File stream (pointing at ACQ1D spot)
        // Return               void  : ACQ1D is read from fp at current
	//			        location. Pointer fp is advanced.
        // Note                       : No basis sharing is assumed!
 
 
void acquire1D::read(const string& fn)
  {
  ifstream fp;  				// Construct a file
  fp.open(fn.c_str(),ios::in|ios::binary);	// Open file
  read(fp);					// Read ACQ1D, use overload
  fp.close();					// Close file
  }
 
ifstream& acquire1D::read(ifstream& fp)
  {
  fp.read((char*)&_LS,sizeof(int));	// Read Liouville space dim.
  fp.read((char*)&pos,sizeof(int));	// Read ACQ1D dim. (pos <= ls)
  I = vector<int> (pos);		// Minimum i index array
  A = row_vector(pos);			// Minimum A array (detector)
  B = row_vector(pos);			// Minimum B array (propagator)
  A.read(fp);				// Read ACQ1D complex array A
  B.read(fp);				// Read ACQ1D complex array B
  int itmp;
  for(int i=0; i<pos; i++)
    {
    fp.read((char*)&itmp,sizeof(int));		// Read ACQ1D int array I
    I[i] = itmp;
    }
  fp.read((char*)&DCUTOFF,sizeof(double));	// Read ACQ1D detection cutoff
  LOp.read(fp);				// Read system Liouvillian
basis bs = LOp.Bs();			// Propagator Hilbert space basis
  Sm1.read(fp);				// Read inverse eigenvector array
  det.read(fp);				// Read detection operator
det.put_basis(bs);			// Set det basis same as L
  int inff;				// Flag for siginf existance
  fp.read((char*)&inff,sizeof(int));
  if(inff)		
    {
    siginf.read(fp);			// Read infinite time dens. op.
siginf.put_basis(bs);			// Set det basis same as L
    trinf.read(fp);			// Read infinite time Tr(det*siginf)
    }
  else trinf = 0;
  return fp;
  }

#endif						// acquire1D.cc
