/* BlochAcq.cc *************************************************-*-c++-*-
**									**
**      	                G A M M A				**
**									**
**	  Bloch Acquisition Class                  Implementation	**
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
**                                                                      **
** This file contains the implementation of class BlochAcq. The purpose **
** of this class is to simplify the repeated computation of magneti-    **
** zation values over evenly spaced time increments.                    **
**                                                                      **
*************************************************************************/

#ifndef   BlochAcq_cc_			// Is this file already included?
#  define BlochAcq_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Bloch/BlochAcq.h>		// Include the header
#include <Basics/Gutils.h>		// Know GAMMA error messages
#include <Basics/StringCut.h>		// Know GAMMA Gdec function

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                 CLASS BLOCH ACQUISITION ERROR HANDLING
// ____________________________________________________________________________

	// Input	BAcq  : A Bloch acquire (this)
	//		eidx  : An error index
	//		noret : A flag for a carriage return
	// Output	void  : A error message is sent to std output

void BlochAcq::ACQerror(int eidx, int noret) const
  {
  std::string hdr("Bloch Acquisition");
  std::string msg;
  switch (eidx)
    {
    case 1:  msg = std::string("Prepared Density Operator of Wrong Dimension");
             GAMMAerror(hdr,msg,noret); break; 				// (1)
    case 2:  msg = std::string("Acquisiton Attempt With NO Prepared DensOp Set");
             GAMMAerror(hdr,msg,noret); break; 				// (2)
    default: GAMMAerror(hdr,eidx,noret); break;
    }
  }

volatile void BlochAcq::ACQfatality(int eidx) const
  {
  ACQerror(eidx, 1);
  if(eidx) ACQerror(0);
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                CLASS BLOCH ACQUISITION PRIVATE CONSTRUCTORS
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


	// Input	BAcq : A Bloch acquire (this)
	// Output	void  : Major components of BAcq are set
	//			  ls    the Bloch dimension
	//			  Sm1   The inverse eigenvector of L
	// Note		      : The following must be preset in the
	//			constructors:
	//			  L      - The system Liouvillian
	//			  det    - The detection superoperator
	//			  Minf - Infinite time magnetization vector
	//			           (it can be left NULL)
	//			  trinf  - The infinte time trace
	//			  cutoff - An intensity cutoff level

void BlochAcq::create()
  {
  BS = L.rows();			// Set the Bloch dimension
  matrix LD, S;				// For eigensystem of L
  L.Diagonalize(LD, S);			// Generate eigensystem
  L               = LD;			// Only store the eigenvalues
  Sm1             = inv(S);		// Set L basis inverse
  row_vector detS = det*S;		// <det| * S = <detS| 

//             Determine Minimum Size Based on Prepared System

  int *I2Dtmp;			// Temporary index storage
  I2Dtmp = new int[BS];
  pos=0;				//               t	
  int i;
  for(i=0; i<BS; i++)			// Loop over <det * S|i> elements
    {
    if(square_norm(detS.get(i))>DCUTOFF)
      {
      I2Dtmp[i] = 1;			// Temporary flag good element 
      pos++;				// Counter for good elements
      }
    else I2Dtmp[i] = 0;
    }

//               Allocate Storage Based on the Size pos

  I = std::vector<int> (pos);		// Minimum i index array
  A = row_vector(pos);			// Minimum A array (detector)
  B = row_vector(pos);			// Minimum B array (propagator)

//                Compute Needed Vector Arrays A,B,I

  int p = 0;
  for(i=0; i<BS; i++)			// Loop over Bloch dimension
    if(I2Dtmp[i])			// Only use "good" elements
      {
      I[p] = i;				// Store "good" i index
      A.put(detS.get(i), p);		// Store A vector (partial intensity)
      B.put(L.get(i,i),    p);		// Store B vector (frequency, rate)
      p++;				// Update dimension index counter
      }
  delete [] I2Dtmp;
  }

	// Input	BAcq	: A Bloch acquire (this)
	//		Sp	: Density matrix (operator propagated)
	// Void		TTab   	: The internal prepared density operator is
	//			  set to Sp and the transitions matrix TTab
	//			  is reconstructed if necessary

void BlochAcq::make_table(const col_vector& Sp)
  {
  if(sigmap==Sp && TTab.size())			// If nothing new & we know the
    return;					// trans. table then just return
  sigmap = Sp;					// Reset the prepared system
  int ls = sigmap.size();			// Bloch dimension of sigmap
  if(ls !=  L.rows()) ACQfatality(1);		// Insure sigmap proper size
  col_vector delsig = sigmap; 			// Get matrix for delsig
  if(Minf.size()) delsig -= Minf;		// Subtract Minf if present
  col_vector Sm1sig = Sm1*delsig;           	// Set |inv(S)* delsigma>
  complex *A1 = new complex[pos];		// Intensity vector replacement
  int k;					// Dummy index
  for(k=0; k<pos; k++)				// Generate new intensity vector
    A1[k] = A(k)*Sm1sig.get(I[k]);
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
// A                  CLASS BLOCH ACQUISITION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ---------------------------- Null Constructor ------------------------------

	// Input		none  :
	// Output		BAcq : A NULL BlochAcq (this)

BlochAcq::BlochAcq()
  {
  BS        = 0; 			// Insure Bloch dimension is 0
  pos       = 0;			// No dimension either
  trinf     = complex0;			// Set <adjoint(det)|Minf>=0 
  DCUTOFF   = 1.e-12;			// Set the detection cutoff level
//  TTab.ICUT = 1.e-12;			// Set intensity cutoff (in T.Table)
  }

BlochAcq::BlochAcq(const row_vector& Dvec, const matrix& Lmx, double cutoff)
  {
  det     = Dvec;				// Set the detection vector
  L       = Lmx;				// Set the system Liouvillian
  trinf   = complex0;				// Set <det|Minf>=0 
  DCUTOFF = cutoff;				// Set detection cutoff level
  create();					// Setup rest of BAcq
  }

BlochAcq::BlochAcq(const row_vector& Dvec, const matrix& Lmx,
                                    const col_vector& MVinf,  double cutoff)
  { 
  det     = Dvec;				// Set the detection vector
  L       = Lmx;				// Set the system Liouvillian
  Minf    = MVinf;				// Store inf. time dens. op.
  trinf   = trace(det,Minf);			// Set <det|Minf> = trinf 
  DCUTOFF = cutoff;				// Set the cutoff level
  create();					// Setup rest of BAcq
  }

// ----------------------------------------------------------------------------
// ---------------- Self Construction, Destruction, Assignment ----------------
// ----------------------------------------------------------------------------

BlochAcq::BlochAcq(const BlochAcq& BAcq1)
  {
  BS      = BAcq1.BS;		// Copy the Bloch dimension size
  pos     = BAcq1.pos;		// Copy number of necessary points
  I       = BAcq1.I;		// Copy minimum I index array
  A       = BAcq1.A;		// Copy minimum A array (detector)
  B       = BAcq1.B;		// Copy minimum B array (propagator)
  DCUTOFF = BAcq1.DCUTOFF;	// Copy the detection cutoff value
  L       = BAcq1.L;		// Copy the system Liouvillian
  Sm1     = BAcq1.Sm1;		// Copy the Liouville inverse eigenvectors
  det     = BAcq1.det;		// Copy the system detection operator
  Minf    = BAcq1.Minf;		// Copy the infinite time density matrix
  trinf   = BAcq1.trinf;	// Copy the infinite time trace
  }

BlochAcq::~BlochAcq() { }

BlochAcq& BlochAcq::operator = (const BlochAcq& BAcq1)
{
  BS      = BAcq1.BS;		// Copy the Bloch dimension
  pos     = BAcq1.pos;		// Copy number of necessary points
  I       = BAcq1.I;		// Copy minimum i index array
  A       = BAcq1.A;		// Copy minimum A array (detector)
  B       = BAcq1.B;		// Copy minimum B array (propagator)
  DCUTOFF = BAcq1.DCUTOFF;	// Copy the detection cutoff value
  L       = BAcq1.L;		// Copy the system Liouvillian
  Sm1     = BAcq1.Sm1;		// Copy the Liouville inverse eigenvectors
  det     = BAcq1.det;		// Copy the system detection operator
  Minf    = BAcq1.Minf;		// Copy the infinite time density matrix
  trinf   = BAcq1.trinf;	// Copy the infinite time trace

  return (*this);
}

// ____________________________________________________________________________
// B                      ALTERATION AND ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions allow users to access the internals of the 1D acquisition.

   Function   Return                         Value Returned
   --------  --------   -------------------------------------------------------
      L2     super_op   Liouvillian evolution superoperator
      D      gen_op     Detection operator
      S      matrix     Eigenvector array in Bloch dimension
      Sinv   matrix     Eigenvector inverse array in Bloch dimension
      TTable TTable1D   Transitions table                                    */

/*
const super_op& BlochAcq::L2()     const { return L; } 
const gen_op&   BlochAcq::D()      const { return det; }
const matrix    BlochAcq::S( )     const { return (L.get_Lbasis()).U(); }
const matrix&   BlochAcq::Sinv()   const { return Sm1; }
const TTable1D& BlochAcq::TTable() const { return TTab; }
*/

/* These functions allow users to reset some of the internals of a 1D 
   acquisition.  In doing so much of the stored acquisition paramters must
   be recomputed.

   Function   Return                         Value Returned
   --------  --------   -------------------------------------------------------
   Detector    void     New detection operator
      D      gen_op     Detection operator
      S      matrix     Eigenvector array in Bloch dimension
      Sinv   matrix     Eigenvector inverse array in Bloch dimension         */


/*
void BlochAcq::Detector(const gen_op& detect)
  {
  det = detect;				// Set detector to detect
  L.LOp_base(det);			// Put det into LOp Hilbert basis
  matrix S = (L.get_Lbasis()).U();	// Get superoperator basis

  col_vector mk = det.superket();       // Start forming detection superbra
                                        //                    t        t
  matrix detS = adjoint(S)*mk;          // Form <det| * S = {S * |det>} 
  detS = adjoint(detS);

//             Determine Minimum Size Based on Detector Only

  int *I2Dtmp;			// Temporary index storage
  I2Dtmp = new int[BS];
  pos=0;				//               t	
  int i;
  for(i=0; i<BS; i++)		// Loop over <det * S|i> elements
    {
    if(square_norm(detS.get(0,i))>DCUTOFF)
      {
      I2Dtmp[i] = 1;			// Temporary flag good element 
      pos++;				// Counter for good elements
      }
    else I2Dtmp[i] = 0;
    }

//               Allocate Storage Based on the Size pos

  I = std::vector<int> (pos);		// Minimum i index array
  A = row_vector(pos);			// Minimum A array (detector)
  B = row_vector(pos);			// Minimum B array (propagator)

//                Compute Needed Vector Arrays A,B,I

  int p = 0;
  for(i=0; i<BS; i++)			// Loop over Bloch dimension
    if(I2Dtmp[i])			// Only use "good" elements
      {
      I[p] = i;				// Store "good" i index
      A.put(detS.get(0,i),p);		// Store A vector (partial intensity)
      B.put(L.get(i,i),   p);		// Store B vector (frequency, rate)
      p++;				// Update dimension index counter
      }

  if(Minf.size())			// If inf. time dens. op.
    trinf = trace(det,Minf);		// Reset <adjoint(det)|Minf>=0
  delete [] I2Dtmp;
  return;
  }
*/

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

	// Input	BAcq	: A Bloch acquire (this)
	// 		M0	: Initial magnetization vector
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

row_vector BlochAcq::T(const col_vector& M0, int npts, double tinc)
  {
  make_table(M0);				// Construct transitions table
  return TTab.T(npts,tinc);			// Use table to make spectrum
  }

void BlochAcq::T(const col_vector& M0, row_vector& data, double tinc)
  {
  make_table(M0);				// Construct transitions table
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
                has fallen below PERCUT*Io will not be included.  The value of
		PERCUT is decimal %, i.e. PERCUT = [0, 1].                   */

	// Input	BAcq	: A Bloch acquire (this)
	// 		M0	: Initial magnetization vector
	//		npts   	: Number of points desired
        // or
	//		data    : 1D row vector for time interferrogram
	//		Fst     : Plot start frequency (Hz)
	//		Ffi     : Plot end frequency (Hz)
	// Output	data  	: Row vector filled with a frequency spectrum.
	// Note		      	: Value PERCUT is a % in decimal form
	//			  of the maximum Lorentzian intensity
	//			  considered worth adding to the spectrum.
	// Note		      	: Each Lorentzian (transition) will have an
	//			  intensity above ICUT*Imax between the points
	//			  bound by the indices ist and ifi.  Only those
	//			  points are used in generating the spectrum
	// Note		      	: To avoid poor resolution, if the frequency
	//			  between points, winc, does not fullfill the
	//			  relationship 5*winc < lwhh, then the
	//			  integrated Lorentzian intensities are used
	//		          rather than the point value.

row_vector BlochAcq::F(const MagVec& M0, int npts, double Fst, double Ffi)
  {
  make_table(M0);				// Construct transitions table
  return TTab.F(npts,Fst,Ffi);			// Use table to make spectrum
  }

void BlochAcq::F(const MagVec& M0, row_vector& data ,double Fst,double Ffi)
  {
  make_table(M0);				// Construct transitions table
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

	// Input	BAcq	: A Bloch acquire (this)
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
 
/*
row_vector BlochAcq::FD(const gen_op& sigmap,int npts,double Fst,double Ffi)
  {
  make_table(sigmap);				// Construct transitions table
  return TTab.FD(npts,Fst,Ffi);			// Use table to make spectrum
  }

void BlochAcq::FD(const gen_op& sigmap,row_vector& data,double Fst,double Ffi)

  {
  make_table(sigmap);				// Construct transitions table
  TTab.FD(data,Fst,Ffi);			// Use table to make spectrum
  }
*/


// ____________________________________________________________________________ 
// F                     DIRECT TRANSITION TABLE GENERATION
// ____________________________________________________________________________

	// Input	BAcq : A Bloch acquire (this)
	//              M0   : Starting magnetization vector
	// Output	TTab  : Transitions table
	//			Re(<i|TTab|0>) = transition i decay rate
	//			Im(<i|TTab|0>) = transition i frequency
	//			   <i|TTab|1>  = transition i intensity
	// Note		      : Frequencies and Rates are in 1/sec
	// Note		      : This sets TTab and sigmap internally

const TTable1D& BlochAcq::table(const col_vector& M0)
  { make_table(M0); return TTab; }

const TTable1D& BlochAcq::table() const { return TTab; }

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

/*
void BlochAcq::offset(double F, int inHz)         { TTab.offset(F,inHz); }
void BlochAcq::offset(double F, int tr, int inHz) { TTab.offset(F,tr,inHz); }
void BlochAcq::FRscale(double Fscf)               { TTab.FRscale(Fscf); }
void BlochAcq::FRscale(double Fscf, int tr)       { TTab.FRscale(Fscf, tr); }
void BlochAcq::Iscale(double Iscf)                { TTab.Iscale(Iscf); }
void BlochAcq::Iscale(double Iscf, int tr)        { TTab.Iscale(Iscf, tr); }
void BlochAcq::broaden(double LWR, int inHz)      { TTab.broaden(LWR, inHz); }
void BlochAcq::broaden(double LWR,int tr,int Hz)  { TTab.broaden(LWR,tr,Hz); }
void BlochAcq::resolution(double res)             { TTab.resolution(res); }
void BlochAcq::pcorrect(double Wpivot,complex& P) { TTab.pcorrect(Wpivot,P); }
complex BlochAcq::pcorrect(double& w0,double w1,int ord)
                                           { return TTab.pcorrect(w0,w1,ord); }
*/

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

/*
double BlochAcq::Wmax()  const   { return TTab.FRmax(); }
double BlochAcq::LWmax() const   { return TTab.LWmax(); }
void   BlochAcq::setSort(int sf) { TTab.setSort(sf); }
void   BlochAcq::setConv(int cf) { TTab.setConv(cf); }
*/

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



	// Input	BAcq	: A Bloch acquire (this)
	// 		sigmap	: Density matrix (operator propagated)
	// 		ostr	: An output stream
	// Output	TTab  : Matrix representation of the 1D spectrum
	//			Re(<i|TTab|0>) = transition i decay rate
	//			Im(<i|TTab|0>) = transition i frequency
	//			   <i|TTab|1>  = transition i intensity
	//			is placed into the output stream
	// Note		      : This sets TTab, sigmap internally

/*
void BlochAcq::table(const gen_op& sigmap, std::ostream& ostr)
  {
  make_table(sigmap);
  TTab.print(ostr);
  }
*/



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
 
/*
void BlochAcq::table(std::ostream& ostr) const
  {
  if(!TTab.rows())				// Check transitions exist
    {						// If not, print a warning
    ACQerror(2); 				// then return and continue
    return;
    }
  TTab.print(ostr);
  }
*/


// ____________________________________________________________________________
// F                  CLASS BLOCH ACQUISITION AUXILIARY FUNCTIONS
// ____________________________________________________________________________


/*
int BlochAcq::ls() const { return BS; }

	// Input	BAcq : A Bloch acquire (this)
	// Output	BS    : Bloch dimension dimension


int BlochAcq::size() const { return pos; }

	// Input	BAcq : A Bloch acquire (this)
	// Output	size  : Size of the arrays A & B
*/


int BlochAcq::full_size() const { return L.rows(); }

	// Input	BAcq : A Bloch acquire (this)
	// Output	size  : Full size of operators
	//			the Hilbert space squared

/*

int BlochAcq::transitions() const { return TTab.size(); }

	// Input	BAcq : A Bloch acquire (this)
	// Output	ntr   : Number of transitions in 1D spectrum

// sosi - what does this now do?
// error as LW is redefined below
void BlochAcq::parameters(matrix& mx, double& SW,
                    double& LW, double& dt, int N, int pf) const 

        // Input        BAcq : A Bloch acquire (this)
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

*/

// ____________________________________________________________________________ 
// G                 CLASS BLOCH ACQUISITION ASCII BASED I/O FUNCTIONS
// ____________________________________________________________________________ 

// ----------------- ASCII Output Functions (For Debugging) -------------------

/*              Input           BAcq  : Bloch acquire
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
    <<                  Same as print(ostr) function   
                                                                             */

std::ostream& BlochAcq::print(std::ostream& ostr, bool pf) const
  {
  std::string hdr("Bloch Acqusition");			// Basic header
  std::string spc;						// String to center
  if(!BS)						// If empty acquisition
    {							// do this
    hdr = std::string("Empty ") + hdr;			//   Appropriate header
    spc = std::string(40-hdr.length()/2, ' ');		//   Space to center
    ostr << "\n" << spc << hdr << "\n";			//   Output centered
    return ostr;					//   Exit
    }

  spc = std::string(40-hdr.length()/2, ' ');			// String to center
  ostr << "\n" << spc << hdr;				// Output header center
  ostr << "\n" << spc << std::string(hdr.length(), '=');	// Ouput underline cent
  hdr = Gdec(pos) + std::string(" Component");		// Line for # comps
  if(pos > 1) hdr += std::string("s");
  hdr += std::string(" From Bloch Vector Dimension Of ");
  hdr += Gdec(full_size());
  spc = std::string(40-hdr.length()/2, ' ');			// Length to center
  ostr << "\n\n" << spc << hdr;				// Ouput centerd comps.
  spc = std::string("\n        ");				// Space center lines

  ostr << "\n";
  ostr << spc << "                           -1                        -1";
  ostr << spc << "Tr(i) = <D*S|i><i|LD|i><i|S  *sigma> = A(i)*B(i)*<i|S  *sigma>";
  ostr << "\n";

  std::string tmp = Gdec(pos);
  int ilen = tmp.length();
  tmp = Gdec(full_size());
  int lslen = tmp.length();
  int i=0;
  bool Areal = A.is_real();				// Flag if A real
  bool Bimag = B.is_imaginary();			// Flag if B imaginary
  for(i=0; i<pos; i++)
    {
    ostr << "\n\t" << Gdec(i,ilen) << ". A = <D*S|"
         << Gdec(I[i], lslen) << "> : ";
         if(Areal) ostr << A.getRe(i);
         else      ostr << A.get(i);
    ostr << "    B = <" << Gdec(I[i], lslen) << "|D|"
         << Gdec(I[i], lslen) << "> : ";
         if(Bimag) ostr << B.getIm(i) << "i";
         else      ostr << B.get(i);
    }

  if(Minf.size())
    {
    hdr = "Infinite Time Contribution";
    if(Minf.norm() < 1.e-10)
      {
      hdr = "No " + hdr;
      spc = std::string(40-hdr.length()/2, ' ');
      ostr << "\n\n" << spc << hdr << "\n";
      }
    else 
      ostr << "\n\n\t" << hdr << ": " << trinf;
    if(pf)
      ostr << "\nInfinite Time Density Operator:\n" << Minf;
    }
  return ostr;
  }


std::ostream& BlochAcq::print(std::ostream& ostr, const col_vector& sigp, bool pf) const
  {
  std::string hdr;
  std::string spc;
  if(!BS)
    {
    hdr = std::string("Empty Bloch Acquisition");
    spc = std::string(40-hdr.length()/2, ' ');
    ostr << "\n" << spc << hdr << "\n";
    return ostr;
    }
  col_vector delsig = sigp;			// Get matrix for delsig
  if(Minf.size()) delsig -= Minf;		// Subtract Minf if present
  matrix Sm1sig = Sm1*delsig;			// Set |inv(S)* delsigma>

  double nsum1=0, nsum2=0, nsum3=0;
  ostr << "\n\t\t" << pos << " Non-zero"	// Output total transitions
       << " Components From A Bloch "		// used in the calculation
       << " Dimension Of " << full_size()
       << "\n";
  int i=0;
  std::string kfmt, ifmt;				// First we set the i&k
  if(pos<10)        kfmt = "%1i";		// indices output format
  else if(pos<100)  kfmt = "%2i";		// so that we get nice
  else if(pos<1000) kfmt = "%3i";		// columns
  else              kfmt = "%4i";
  if(BS<10)         ifmt = "%1i";
  else if(BS<100)   ifmt = "%2i";
  else if(BS<1000)  ifmt = "%3i";
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

  if(pf)
    {
    ostr << "\n\n\t" << "Total norm(A): "
         << nsum1 << "; Total norm(C): "
         << nsum2 << "; Total norm(AC): "
         << Gform("%16.10f", nsum3);
    }

  if(Minf.size())
    {
    hdr = "Infinite Time Contribution";
    if(Minf.norm() < 1.e-10)
      {
      hdr = "No " + hdr;
      spc = std::string(40-hdr.length()/2, ' ');
      ostr << "\n\n" << spc << hdr << "\n";
      }
    else 
      ostr << "\n\n\t" << hdr << ": " << trinf;
    if(pf)
      ostr << "\nInfinite Time Density Operator:\n" << Minf;
    }
  return ostr;
  }

std::ostream &operator << (std::ostream &ostr, const BlochAcq &ACQ)
  { ACQ.print(ostr); return ostr; }


// ____________________________________________________________________________
// H                 CLASS BLOCH ACQUISITION BINARY BASED I/O FUNCTIONS
// ____________________________________________________________________________

// ---------------------------- Binary Output ---------------------------------


	// Input		BAcq : Acquire1D (this)
        //                      fn    : Filename
        // Return               void  : Acquisition  BAcq is written
        //                              to a file called fn
        // Note                       : The file format is BINARY

/*
void BlochAcq::write(const std::string& fn) const
  {
  ofstream fp;					// Construct a file
  fp.open(fn.c_str(),ios::out|ios::binary);	// Open file
  write(fp);                                    // Write Op
  fp.close();                                   // Close file
  }
*/


 
	// Input		BAcq : Acquire1D (this)
        //                      fp    : File stream (pointing at BAcq)
        // Return               void  : Acquisition  BAcq is written
        //                              into file fp at current location
        // Note                       : The file format is BINARY and the file
	//				point fp is advanced
 
/*
ofstream& BlochAcq::write(ofstream& fp) const
  {
  fp.write((char*)&BS,sizeof(int));		// Write Bloch dimension dim.
  fp.write((char*)&pos,sizeof(int));		// Write BAcq dim. (pos <= ls)
  A.write(fp);					// Write BAcq complex array A
  B.write(fp);					// Write BAcq complex array B
  for(int i=0; i<pos; i++)
    fp.write((char*)I[i],sizeof(int));		// Write BAcq int array I
  fp.write((char*)&DCUTOFF,sizeof(double));	// Write BAcq cutoff
  L.write(fp);				// Write system Liouvillian
  Sm1.write(fp);			// Write inverse eigenverctor array
  det.write(fp);			// Write detection operator
  int inff;
  if(Minf.size())		
    {
    inff = 1;
    fp.write((char*)&inff,sizeof(int));
    Minf.write(fp);			// Write infinite time dens. op.
    trinf.write(fp); 			// Write infinite time Tr(det*Minf)
    }
  else
    {
    inff = 0;
    fp.write((char*)&inff,sizeof(int));
    }
  return fp;
  }
*/
 
// ---------------------------- Binary Input ----------------------------------
 
 
	// Input		BAcq : Acquire1D (this)
        //                      fn    : Filename
        // Return               void  : BAcq is read from file fn
        // Note                       : No basis sharing is assumed!
 
/*
void BlochAcq::read(const std::string& fn)
  {
  ifstream fp;  				// Construct a file
  fp.open(fn.c_str(),ios::in|ios::binary);	// Open file
  read(fp);					// Read BAcq, use overload
  fp.close();					// Close file
  }
*/
 
 
 
	// Input		BAcq : Acquire1D (this)
        //                      fp    : File stream (pointing at BAcq spot)
        // Return               void  : BAcq is read from fp at current
	//			        location. Pointer fp is advanced.
        // Note                       : No basis sharing is assumed!
 
/*
ifstream& BlochAcq::read(ifstream &fp)
  {
asdfasdf: not found
  fp.read((char*)&pos,sizeof(int));		// Read BAcq dim. (pos <= ls)
  I = std::vector<int> (pos);			// Minimum i index array
  A = row_vector(pos);			// Minimum A array (detector)
  B = row_vector(pos);			// Minimum B array (propagator)
  A.read(fp);				// Read BAcq complex array A
  B.read(fp);				// Read BAcq complex array B
  int itmp;
  for(int i=0; i<pos; i++)
    {
    fp.read((char*)&itmp,sizeof(int));		// Read BAcq int array I
    I[i] = itmp;
    }
  fp.read((char*)&DCUTOFF,sizeof(double));	// Read BAcq detection cutoff
  L.read(fp);				// Read system Liouvillian
basis bs = L.get_basis();		// Propagator Hilbert space basis
  Sm1.read(fp);				// Read inverse eigenvector array
  det.read(fp);				// Read detection operator
det.put_basis(bs);			// Set det basis same as L
  int inff;				// Flag for Minf existance
  fp.read((char*)&inff,sizeof(int));
  if(inff)		
    {
    Minf.read(fp);			// Read infinite time dens. op.
Minf.put_basis(bs);			// Set det basis same as L
    trinf.read(fp);			// Read infinite time Tr(det*Minf)
    }
  else trinf = 0;
  return fp;
  }
*/
 
#endif							// BlochAcq.cc
