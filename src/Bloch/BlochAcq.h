/* BlochAcq.h ***************************************************-*-c++-*-
**									**
** 	                       G A M M A				**
**									**
**	Bloch Acquisiiton Class                          Interface	**
**									**
**	Copyright (c) 2001 						**
**	Dr. Scott A. Smith				 		**
**      National High Magnetic Field Laboratory 	                **
**      1800 E. Paul Dirac Drive                        	        **
**      Tallahassee Florida, 32306-4005                         	**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** This file contains the implementation of class BlochAcq. The purpose	**
** of this class is to simplify the repeated computation of magneti-	**
** zation values over evenly spaced time increments.			**
**									**
*************************************************************************/

#ifndef   BlochAcq_h_			// Is this file already included?
#  define BlochAcq_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Bloch/MagVec.h>		// Know Bloch magnetization vectors
#include <Matrix/complex.h>		// Knowledge of complex numbers
#include <Matrix/matrix.h>		// Know GAMMA matrices
#include <Matrix/row_vector.h>		// Know GAMMA row vectors
#include <Matrix/col_vector.h>		// Know GAMMA column vectors
#include <Level2/TrnsTable1D.h>		// Know GAMMA transition tables
#include <vector>			// Know libstdc++ STL vectors

class BlochAcq 
  {
  int BS;			// Bloch dimension

  matrix L;			// System "Liouvillian" matrix
  matrix Sm1;			// L matrix inverse eigenvectors

  row_vector det;		// The detection operator
  col_vector Minf;		// Infinite time magnetization vector
  complex trinf;		// Trace at infinite time Tr(det*Minf)

  int pos;			// Acquire1D dimensions (pos <= ls)
  row_vector A;			// A array (detector equivalent)
  row_vector B;			// B array (time propagator equivalent)
  std::vector<int> I;		// Array for indexing
  double DCUTOFF;		// Detection cut value

  TTable1D TTab;		// Transitions table for 1D spectrum
  col_vector sigmap;		// Prepared density operator
  double Icut;			// Transitiion intensity cut

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS BLOCH ACQUISITION ERROR HANDLING
// ____________________________________________________________________________

void ACQerror(int eidx, int noret=0) const;

        // Input        ACQ1D : An BlochAcq (this)
	//		eidx  : An error index
	//		noret : Flag for carriage return
        // Output       void  : Error Message Output

volatile void ACQfatality(int eidx=0) const;

        // Input        ACQ1D : An BlochAcq (this)
	//		eidx  : An error index
        // Output       void  : Error Message Output, Execution Stopped

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

void create( );

	// Input	ACQ1D : An BlochAcq (this)
	// Output	void  : Major components of ACQ1D are set
	//			  ls    the Liouville space
	//			  Sm1   The inverse eigenvector of L
	// Note		      : The following MUST be preset in the
	//			constructors:
	//			  L      - The system Liouvillian
	//			  det    - The detection superoperator
	//			  Minf - The infinite time density matrix
	//			           (it can be left NULL)
	//			  trinf  - The infinte time trace
	//			  ACUTOFF - An intensity cut level


void make_table(const col_vector& Sp);

        // Input	ACQ1D	: An BlochAcq (this)
        //		Sp	: Density matrix (operator propagated)
        // Output	void	: The internal prepared density operator is
	//			  set to Sp and the transitions matrix TTab
	//			  is reconstructed if necessary


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 public:
// friend acquire2D;			// Let acquire2D have everything

// ____________________________________________________________________________
// A               CLASS BLOCH ACQUISITION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
// --------------------------- Simple Constructors ----------------------------
// ----------------------------------------------------------------------------

MSVCDLC BlochAcq();
MSVCDLC BlochAcq(const row_vector& det, const matrix& L, 
                                   const col_vector& Minf, double cut=1.e-12);
MSVCDLC BlochAcq(const row_vector& det, const matrix& L, double cut=1.e-12);

// ----------------------------------------------------------------------------
// ---------------- Self Construction, Destruction, Assignment ----------------
// ----------------------------------------------------------------------------

MSVCDLC      BlochAcq(const BlochAcq& ACQ1);
MSVCDLC      ~BlochAcq();
MSVCDLL BlochAcq& operator = (const BlochAcq& ACQ1);

// ____________________________________________________________________________
// B                    ALTERATION AND ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions allow users to access the internals of the 1D acquisition.

   Function   Return                         Value Returned
   --------  --------   -------------------------------------------------------
      L2     super_op   Liouvillian evolution superoperator
      D      gen_op     Detection operator
      S      matrix     Eigenvector array in Liouville space
      Sinv   matrix     Eigenvector inverse array in Liouville space
      TTable TTable1D   Transitions table                                    */

/*
const super_op& L2()     const;
const gen_op&   D()      const;
const matrix    S()      const;
const matrix&   Sinv()   const;
const TTable1D& TTable() const;

void Detector(const gen_op& detect);
*/

        // Input        ACQ1D : An BlochAcq (this)
        // Input        det   : Detection operator
        // Output       void  : Acquire1D (this) is modified for
        //                      use of new detection operator

// ____________________________________________________________________________
// C                    Time Domain Spectra Generation
// ____________________________________________________________________________

/* These functions take the 1D transitions table and generate a time domain
   spectrum.  The user just inputs the prepared system (magnetization vector)
   that will be evolved during the spectrum acquisition, the number of points
   desired and the time increment to use between them.  There are a few other
   factors which will influence the spectra generated by these functions.

	ICUT	This is an intensity cut value.  Any transition whose
                intensity (norm) is below this value will not contribute to
                the output spectrum.
        PERCUT  This is a percentage cut value.  For each transition,
                points at times past where the transition intensity has
                decayed below PERCUT*Io will not be included.  The value of
		PERCUT is decimal %, i.e. PERCUT = [0, 1].                   */

	// Input	ACQ1D	: An BlochAcq (this)
	// 		M0	: Magnetization vector evolved
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

MSVCDLL row_vector T(const col_vector& M0, int npts,         double tinc);
MSVCDLL void       T(const col_vector& M0, row_vector& data, double tinc);

// ____________________________________________________________________________ 
// D                        Frequency Domain Spectra
// ____________________________________________________________________________

/* These functions take the 1D transitions table and generate a frequency
   domain spectrum.  The user just inputs the prepared system (magnetization
   vector) that will be evolved during the spectrum acquisition, the
   number of points desired and the spectral range. There are a few other
   factors which will influence the spectra generated by these functions.

        ICUT    This is an intensity cutoff value.  Any transition whose
                intensity (norm) is below this value will not contribute to
                the output spectrum.
        PERCUT  This is a percentage cutoff value.  For each transition,
                points at frequencies outside where the transition intensity
                has fallen below PERCUT*Io will not be included.  The value of
                PERCUT is decimal %, i.e. PERCUT = [0, 1].                   */

        // Input        BAcq    : A Bloch acquire (this)
	// 		M0	: Magnetization vector evolved
        //              npts    : Number of points desired
        // Input        BAcq    : A Bloch acquire (this)
        //              M0      : Initial magnetization vector
        //              npts    : Number of points desired
        // or
        //              data    : 1D row vector for time interferrogram
        //              Fst     : Plot start frequency (Hz)
        //              Ffi     : Plot end frequency (Hz)
        // Output       data    : Row vector filled with a frequency spectrum.
        // Note                 : Value PERCUT is a % in decimal form
        //                        of the maximum Lorentzian intensity
        //                        considered worth adding to the spectrum.
        // Note                 : Each Lorentzian (transition) will have an
        //                        intensity above ICUT*Imax between the points
        //                        bound by the indices ist and ifi.  Only those
        //                        points are used in generating the spectrum
        // Note                 : To avoid poor resolution, if the frequency
        //                        between points, winc, does not fullfill the
        //                        relationship 5*winc < lwhh, then the
        //                        integrated Lorentzian intensities are used
        //                        rather than the point value.

row_vector F(const MagVec& M0, int npts,         double Fst, double Ffi);
void       F(const MagVec& M0, row_vector& data, double Fst, double Ffi);

// ____________________________________________________________________________ 
// E             Frequency Domain Spectra In Derivative Mode
// ____________________________________________________________________________

/* These functions take the 1D transitions table and generate a frequency
   domain spectrum in derivative mode (ESR/EPR).  The user just inputs the
   prepared system (density operator) that will be evolved during the spectrum
   acquisition, the number of points desired and the spectral range. There are
   a few other factors which will influence the spectra generated by these 
   functions.

	ICUT	This is an intensity cut value.  Any transition whose
                intensity (norm) is below this value will not contribute to
                the output spectrum.
        PERCUT  This is a percentage cut value.  For each transition,
                points at frequencies outside where the transition intensity
                has fallen below PERCUT*Io will not be included.  The value of
		PERCUT is decimal %, i.e. PERCUT = [0, 1].                   */

	// Input	ACQ1D	: An BlochAcq (this)
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
row_vector FD(const gen_op& sigma,int npts,double Fst,double Ffi);
void FD(const gen_op& sigma,row_vector& data,double Fst,double Ffi);
*/


// ____________________________________________________________________________
// F                   DIRECT TRANSITION TABLE GENERATION
// ____________________________________________________________________________

/* These functions return a 1D transitions table given a prepared system
   (density operator) that will be evolved during the an acquisition.  There
   are a few other factors which will influence the table generated by these 
   functions.

	ICUT	This is an intensity cut value.  Any transition whose
                intensity (norm) is below this value will be removed from the
  		returned table.
*/

MSVCDLL const TTable1D& table(const col_vector& M0);
MSVCDLL const TTable1D& table() const;

	// Input	ACQ1D : An BlochAcq (this)
	// 		sigma : Density matrix (operator propagated)
	// Output	mx    : Matrix representation of the 1D spectrum
	//			Re(<i|mx|0>) = transition i decay rate
	//			Im(<i|mx|0>) = transition i frequency
	//			   <i|mx|1>  = transition i intensity
        // Note               : Frequencies and Rates are in 1/sec

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
void offset(double     F,            int inHz=1);
void offset(double     F,    int tr, int inHz);
void FRscale(double    Fscf);
void FRscale(double    Fscf, int tr);
void Iscale(double     Iscf);
void Iscale(double     Iscf, int tr);
void broaden(double    LWR,          int inHz=1);
void broaden(double    LWR,  int tr, int inHz);
void resolution(double res);
void pcorrect(double   Wpivot, complex& P);
complex pcorrect(double& w0, double w1, int order=5);
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
double Wmax()  const;
double LWmax() const;
void   setSort(int sf);
void   setConv(int cf);
*/

/* The above function setConv sets the internal value of FRQCONV in class
   TrnsTable1D.  It is only active when printing the transitions table and then
   only when the frequencies are ouput in either PPM (FRQTYPE=1) or in Gauss
   (FRQTYPE=2) For the former users should set FRQCONV to be the Larmor 
   frequency in MHz as this value will be used for Hz -> PPM. For the latter
   set FRQCONV to be the electron g value                                    */

//void table(const gen_op& sigma, ostream& ostr);

        // Input        ACQ1D : An BlochAcq (this)
        //              sigma : Density matrix (operator propagated)
        //              ostr  : An output stream
        //              Icut  : Intensity cut
        //              type  : Flag for printout type 
        //                              0 = Hertz (default) 
        //                              1 = PPM 
        //                              2 = Gauss 
        //              Om    : Omega needed for conversion to PPM
        // Output       dtab  : Matrix representation of the 1D spectrum
        //                      Re(<i|dtab|0>) = transition i decay rate
        //                      Im(<i|dtab|0>) = transition i frequency
        //                         <i|dtab|1>  = transition i intensity
        //                      is placed into the output stream
        // Note               : This sets dtab, sigmap, and Icut
        //                      internal to the class


//void table(ostream& ostr) const;

        // Input        ostr  : An output stream
        //              type  : Flag for printout type 
        //                              0 = Hertz (default) 
        //                              1 = PPM 
        //                              2 = Gauss 
        //              Om    : Omega needed for conversion to PPM
        // Output       void  : A table of transition infomation
        //                      is placed into the output stream
        //                      Re(<i|dtab|0>) = transition i decay rate
        //                      Im(<i|dtab|0>) = transition i frequency
        //                         <i|dtab|1>  = transition i intensity
        //              cf    : Conversion factor if needed
        //                          type=1 ==> cf = Om in MHz
        //                          type=2 ==> cf = g (unitless)
        // Note               : Frequencies and Rates are in 1/sec
        // Note               : The array dtab should have been filled
        //                      prior to this function call

//                                   R     1
//                         lwhh   = -- = -----
//                             Hz   PI   PI*T
//                                           2
 
// ____________________________________________________________________________
// F                   CLASS BLOCH ACQUISITION AUXILIARY FUNCTIONS
// ____________________________________________________________________________

///F_list size        - Liouville space size
///F_list size        - Reduced space size
///F_list full_size   - Full space size


//int ls() const;

	// Input	ACQ   : An BlochAcq (this)
        // Output       ls    : Liouville space size


//int size() const;

	// Input	ACQ   : An BlochAcq (this)
	// Output	size  : Size of the arrays A & B


MSVCDLL int full_size() const;

	// Input	ACQ   : An BlochAcq (this)
	// Output	size  : Full size of operators
	//			the Hilbert space squared

//int transitions() const;

	// Input	ACQ   : An BlochAcq (this)
	// Output	ntr   : Number of transtions in 1D spectrum


//void parameters(matrix& mx, double& SW, double& LW,
//                                     double& dt, int N, int pf=0) const;

        // Input        ACQ1D : An BlochAcq (this)
        //              mx    : Array of acquire parameters
	//		SW    : Suggested spectral width (Hz)
        //              LW    : Suggested linewidth (Hz)
        //              dt    : Suggested dwell time (sec)
	//		N     : Number of FID points
        //              pf    : Print flag
        // Output       size  : Full size of operators
        //                      the Hilbert space squared

  
// ____________________________________________________________________________
// G                 CLASS BLOCH ACQUISITION ASCII BASED I/O FUNCTIONS
// ____________________________________________________________________________

///F_list print                  - Send acquisition to output stream.

/* 		Input		ACQ   : Acquire
        	ostr		      : An output stream for info
        			sigmap: State of system at acquisition
        			tinc  : Dwell time in seconds
        			npts  : Number of points
        				 0 = Evolve From t=0 To Point
        				!0 = Evolve From Point To Point
		Output		ostr  : The output stream modified by
				        acquisition parameters/computation 

	Function				Output
  --------------------------------------------------------------------------
    print(ostr)		Core elements of the acquisition computation
    print(ostr,sigmap)  Core elements of the acquisiton from state sigmap
    printT(.....)       Info on acquisiton time tomain calculation
    <<			Same as print(ostr) function
                                                                             */

MSVCDLL        std::ostream& print(std::ostream& ostr,                         bool pf=false) const;
MSVCDLL        std::ostream& print(std::ostream& ostr, const col_vector& sigp, bool pf=false) const;
MSVCDLL friend std::ostream& operator << (std::ostream&    ostr, const BlochAcq &ACQ);
};

#endif								// BlochAcq.h
