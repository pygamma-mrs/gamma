/* acquire1D.h **************************************************-*-c++-*-
**									**
** 	                       G A M M A				**
**									**
**	Class Acquire1D	                           Interface		**
**									**
**	Copyright (c) 1994				 		**
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
** This file contains the implementation of class acquire1D.  The 	**
** purpose of this class is to simplify the computation	of multiple	**
** expectation values in a system which evolves	through a single delay.	**
**									**
*************************************************************************/

#ifndef   acquire1D_h_			// Is this file already included?
#  define acquire1D_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/complex.h>		// Knowledge of complex numbers
#include <HSLib/GenOp.h>		// Know Hilbert space operators
#include <HSLib/HSprop.h>		// Know Hilbert space propagators
#include <LSLib/SuperOp.h>		// Knowledge of superoperators
#include <LSLib/LSprop.h>		// Know Liouville space props
#include <Level2/TrnsTable1D.h>		// Know 1D transitions tables
#include <string>			// Knowledge of libstdc++ strings
#include <vector>			// Know libstdc++ STL vectors
//#include <basis.h>			// Knowledge of bases

//class acquire2D;			// Knowledge of 2D acquisitions

class acquire1D 
  {
  int _LS;			// Liouville space dimension
  int pos;			// Acquire1D dimension (pos <= ls)
  row_vector A;			// A array (detector equivalent)
  row_vector B;			// B array (time propagator equivalent)
  std::vector<int> I;		// Indexing array (constributors to pos)
  double DCUTOFF;		// Detection cutoff value

  super_op LOp;			// System Liouvillian or evolution propagator
  matrix Sm1;			// Liouville space inverse eigenvectors

  gen_op det;			// The detection operator
  gen_op siginf;		// Infinite time density matrix
  complex trinf;		// Trace at infinite time Tr(det*siginf)

  TTable1D TTab;		// Transitions table for 1D spectrum
  gen_op sigmap;		// Prepared density operator
  double Icutoff;		// Transitiion intensity cutoff

  double delt;			// Time increment




// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS ACQUIRE ERROR HANDLING
// ____________________________________________________________________________


        // Input        ACQ1D : An acquire1D (this)
	//		eidx  : An error index
	//		noret : Flag for carriage return
        // Output       void  : Error Message Output

         void ACQerror(int eidx, int noret=0) const;
volatile void ACQfatal(int eidx=0)            const;

        // Input        ACQ1D : An acquire1D (this)
	//		eidx  : An error index
        // Output       void  : Error Message Output, Execution Stopped

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

        // Input        ACQ1D : An acquire1D (this)
        // Output       void  : Major components of ACQ1D are set
        //                        ls    the Liouville space
        //                        Sm1   The inverse eigenvector of L
        // Note               : The following must be preset in the
        //                      constructors:
        //                        L      - The system Liouvillian
        //                        det    - The detection superoperator
        //                        siginf - The infinite time density matrix
        //                                 (it can be left NULL)
        //                        trinf  - The infinte time trace
        //                        cutoff - An intensity cutoff level         */

void create();
void createU();

// ____________________________________________________________________________
// iii                CLASS ACQUIRE1D PRIVATE CONSTRUCTORS
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


	// Input	ACQ1D : An acquire1D (this)
	// Output	void  : Major components of ACQ1D are set
	//			  ls    the Liouville space
	//			  Sm1   The inverse eigenvector of L
	// Note		      : The following MUST be preset in the
	//			constructors:
	//			  L      - The system Liouvillian
	//			  det    - The detection superoperator
	//			  siginf - The infinite time density matrix
	//			           (it can be left NULL)
	//			  trinf  - The infinte time trace
	//			  ACUTOFF - An intensity cutoff level


void make_table(const gen_op& Sp);

        // Input	ACQ1D	: An acquire1D (this)
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
// A                 CLASS ACQUIRE CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
// --------------------------- Simple Constructors ----------------------------
// ----------------------------------------------------------------------------

MSVCDLC acquire1D();				// Null constructor
MSVCDLC acquire1D(const acquire1D& ACQ1);	// Self constructor

// ----------------------------------------------------------------------------
// ------------------ Hilbert Space Treatment Constructors --------------------
// ----------------------------------------------------------------------------

/* These constructors set up the 1D acquisition EXCLUDING the effects of
   relaxation.  That is, the system is evolved either under a static 
   Hamiltonian (no specific time increment needed) or under a static
   propagator (keyed to a specific time increment). Both are performed in
   Hilbert space.

           Input        det   : Detection operator
                        H     : The system Hamiltonian
	  		U     : The system Propagator
                        cutoff: Intensity cutoff level
           Output       ACQ1D : Acquire1D (this) is constructed              */

MSVCDLC acquire1D(gen_op& det, gen_op& H);
MSVCDLC acquire1D(gen_op& det, gen_op& H, double cutoff);
MSVCDLC acquire1D(gen_op& det, HSprop& U);
MSVCDLC acquire1D(gen_op& det, HSprop& U, double cutoff);

// ----------------------------------------------------------------------------
// ----------------- Liouville Space Treatment Constructors -------------------
// ----------------------------------------------------------------------------

/* These constructors set up the 1D acquisition with the possibility of
   INCLUDING the effects of relaxation.  That is, the system is evolved under
   a static Liouvillian superoperator in Liouville space.  Furthermore the
   acquistion tracks the "infinite time" density opertator which the system
   will ultimately evolve into.

        // Input        D     : Detection operator
        //              L     : The system Liouvillian
        //              sigi  : Infinite time density matrix
        //              cut   : Intensity cutoff level
        // Output       ACQ1D : Acquire1D (this) is constructed              */


MSVCDLC acquire1D(gen_op& D, super_op& L, gen_op& sigi, double cut=1.e-12);
MSVCDLC acquire1D(matrix& D, super_op& L, gen_op& sigi, double cut=1.e-12);
MSVCDLC acquire1D(gen_op& D, super_op& L,               double cut=1.e-12);
MSVCDLC acquire1D(matrix& D, super_op& L,               double cut=1.e-12);
MSVCDLC acquire1D(gen_op& D, LSprop&   G,               double cut=1.e-12);

// ----------------------------------------------------------------------------
// ------------------------- Destruction, Assignment --------------------------
// ----------------------------------------------------------------------------

MSVCDLC            ~acquire1D();
MSVCDLL acquire1D& operator= (const acquire1D& ACQ1);

// ____________________________________________________________________________
// B                    ALTERATION AND ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions allow users to access the internals of the 1D acquisition.

   Function   Return                         Value Returned
   --------  --------   -------------------------------------------------------
      L      super_op   Liouvillian evolution superoperator
      D      gen_op     Detection operator
      S      matrix     Eigenvector array in Liouville space
      Sinv   matrix     Eigenvector inverse array in Liouville space
      TTable TTable1D   Transitions table                                    */

MSVCDLL const super_op& L()	    const;
MSVCDLL const gen_op&   D()      const;
MSVCDLL const matrix    S()      const;
MSVCDLL const matrix&   Sinv()   const;
MSVCDLL const TTable1D& TTable() const;
MSVCDLL       void      Detector(const gen_op& detect);

        // Input        ACQ1D : An acquire1D (this)
        // Input        det   : Detection operator
        // Output       void  : Acquire1D (this) is modified for
        //                      use of new detection operator

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

MSVCDLL row_vector T(const gen_op& sigmap, double tinc,      int npts);
MSVCDLL row_vector T(const gen_op& sigmap, int npts,         double tinc);
MSVCDLL void       T(const gen_op& sigmap, row_vector& data, double tinc);

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
	//			  bound by the indices ist and ifi.  Only those
	//			  points are used in generating the spectrum
	// Note		      	: To avoid poor resolution, if the frequency
	//			  between points, winc, does not fullfill the
	//			  relationship 5*winc < lwhh, then the
	//			  integrated Lorentzian intensities are used
	//		          rather than the point value.

MSVCDLL row_vector F(const gen_op& sigmap,int npts,double Fst,double Ffi);
MSVCDLL void F(const gen_op& sigmap,row_vector& data,double Fst,double Ffi);

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
 
MSVCDLL row_vector FD(const gen_op& sigma,int npts,double Fst,double Ffi);
MSVCDLL void       FD(const gen_op& sigma,row_vector& data,double Fst,double Ffi);


// ____________________________________________________________________________
// F                   DIRECT TRANSITION TABLE GENERATION
// ____________________________________________________________________________

/* These functions return a 1D transitions table given a prepared system
   (density operator) that will be evolved during the an acquisition.  There
   are a few other factors which will influence the table generated by these 
   functions.

	ICUT	This is an intensity cutoff value.  Any transition whose
                intensity (norm) is below this value will be removed from the
  		returned table.
*/

MSVCDLL const TTable1D& table(const gen_op& sigmap);
MSVCDLL const TTable1D& table() const;

MSVCDLL const TTable1D table_snapshot(const gen_op& sigmap);
MSVCDLL const TTable1D table_snapshot() const;

	// Input	ACQ1D : An acquire1D (this)
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

MSVCDLL void offset(double     F,            int inHz=1);
MSVCDLL void offset(double     F,    int tr, int inHz);
MSVCDLL void FRscale(double    Fscf);
MSVCDLL void FRscale(double    Fscf, int tr);
MSVCDLL void Iscale(double     Iscf);
MSVCDLL void Iscale(double     Iscf, int tr);
MSVCDLL void broaden(double    LWR,          int inHz=1);
MSVCDLL void broaden(double    LWR,  int tr, int inHz);
MSVCDLL void resolution(double res);
MSVCDLL void pcorrect(double   Wpivot, complex& P);
MSVCDLL complex pcorrect(double& w0, double w1, int order=5);

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

MSVCDLL double Wmax()  const;
MSVCDLL double LWmax() const;
MSVCDLL void   setSort(int sf);
MSVCDLL void   setConv(int cf);

/* The above function setConv sets the internal value of FRQCONV in class
   TrnsTable1D.  It is only active when printing the transitions table and then
   only when the frequencies are ouput in either PPM (FRQTYPE=1) or in Gauss
   (FRQTYPE=2) For the former users should set FRQCONV to be the Larmor 
   frequency in MHz as this value will be used for Hz -> PPM. For the latter
   set FRQCONV to be the electron g value                                    */

MSVCDLL void table(const gen_op& sigma, std::ostream& ostr);

        // Input        ACQ1D : An acquire1D (this)
        //              sigma : Density matrix (operator propagated)
        //              ostr  : An output stream
        //              Icut  : Intensity cutoff
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
        // Note               : This sets dtab, sigmap, and Icutoff
        //                      internal to the class


MSVCDLL void table(std::ostream& ostr) const;

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
// F                   CLASS ACQUIRE AUXILIARY FUNCTIONS
// ____________________________________________________________________________

///F_list size        - Liouville space size
///F_list size        - Reduced space size
///F_list full_size   - Full space size

MSVCDLL int ls()          const;		// Liouville space size
MSVCDLL int size()        const;		// Acquisition size
MSVCDLL int full_size()   const;		// Full acquisition size
MSVCDLL int transitions() const;		// Number of transitions

MSVCDLL void parameters(matrix& mx, double& SW, double& LW,
                                     double& dt, int N, int pf=0) const;

        // Input        ACQ1D : An acquire1D (this)
        //              mx    : Array of acquire parameters
	//		SW    : Suggested spectral width (Hz)
        //              LW    : Suggested linewidth (Hz)
        //              dt    : Suggested dwell time (sec)
	//		N     : Number of FID points
        //              pf    : Print flag
        // Output       size  : Full size of operators
        //                      the Hilbert space squared

  
// ____________________________________________________________________________
// G                 CLASS ACQUIRE1D ASCII BASED I/O FUNCTIONS
// ____________________________________________________________________________

///F_list print                  - Send acquisition to output stream.

// ----------------- ASCII Output Functions (For Debugging) -------------------
 
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
    printT(.....)       Info on acquisiton time domain calculation
    <<			Same as print(ostr) function
                                                                             */

MSVCDLL std::ostream& print(std::ostream&  ostr) const;
MSVCDLL std::ostream& print(std::ostream&  ostr, gen_op& sigmap);
MSVCDLL std::ostream& printT(std::ostream& ostr, gen_op& sigmap,
                                         double tinc, int npts=10, int P2P=0);
MSVCDLL friend std::ostream& operator << (std::ostream& ostr, acquire1D& ACQ);

 
// ____________________________________________________________________________ 
// H                 CLASS ACQUIRE1D BINARY BASED I/O FUNCTIONS 
// ____________________________________________________________________________

// ------------------ Binary Output Functions (For Storage) -------------------


	// Input		ACQ1D : Acquire1D (this)
        //                      fn    : Filename
        //                      fp    : File stream (pointing at ACQ1D)
        // Return               void  : The acquisition  ACQ1D is written
        //                              to a file called filename or into the
	//				ouptut file stream fp at its current
	//				location.
        // Note                       : The file format is BINARY

MSVCDLL void           write(const std::string& fn) const;
MSVCDLL std::ofstream& write(std::ofstream&     fp) const;
 

// ------------------- Binary Read Functions (For Storage) --------------------

 
//  void read(string& fn, gen_op& Op1);
 
	// Input		ACQ1D : Acquire1D (this)
        //                      fn   : Filename
        //                      Op1  : A 2nd operator
        // Return               void : Operator Op is read from
        //                             a file called filename. Its
        //                             basis will be set to that of
        //                             Op1 if they are found to agree
        // Note                       : The file format is BINARY
 

//  void read(string& fn, basis& bs);
 
	// Input		ACQ1D : Acquire1D (this)
        //                      fn   : Filename
        // Return               void : Operator Op is read from
        //                             a file called filename. Its
        //                             basis will be set to bs if
        //                             they are found to agree
        // Note                       : The file format is BINARY
 
 

 
	// Input		ACQ1D : Acquire1D (this)
        //                      fn    : Filename
        //                      fp    : File stream (pointing at ACQ1D spot)
        // Return               void  : Operator Op is read from
        //                              a file called filename or from
	//				filestream fp at current location 
        // Note                       : Only WBR is read
        // Note                       : No basis sharing is assumed!
 
MSVCDLL void           read(const std::string& fn);
MSVCDLL std::ifstream& read(std::ifstream&     fp);

};

/*************************************************************************
**									**
** This file contains the implementation of class acquire1D.  The 	**
** purpose of this class is to simplify the computation	of multiple	**
** expectation values in a system which evolves	through a single delay.	**
** It is assumed that system evolution	proceeds according to		**
**									**
**             |sigma(t)> = exp(-Lt)|delsig> + |sigmainf>		**
**									**
** where L is the active Liouvillian.  An expection value for operator	**
** M would then be given by						**
**									**
**           t               t                      t			**
**         <M |sigma(t)> = <M |exp(-Lt)|delsig> + <M |sigmainf>		**
**									**
** where L is the active Liouvillian.  Working in the eigenbasis of L,	**
** where								**
**                                          -1				**
**	                       L = S * D * S				**
**									**
** our working equation becomes						**
**									**
**        t               t            -1             t			**
**      <M |sigma(t)> = <M S|exp(-Dt)|S  *delsig> + <M |sigmainf>	**
**									**
** Noting that the trace over the infinite time density matrix in a	**
** constant, K, and that the superoperator D is diagonal, it is clear 	**
** that the expection value will oscillate in time and be a sum of 	**
** decaying exponentials.						**
**									**
**                              ls					**
**                              ---					**
**                t             \ 			            	**
**              <M |sigma(t)> = /   I  * exp(-D*t)   + K		**
**                              ---	 i         i			**
**                               i 					**
**	th	   							**
** The i   exponential has magnitude given by				**
**									**
**                                    -1             			**
**                    I = <M*S|i>*<i|S  *delsig>			**
**                     i						**
**									**
** and oscillates with the frequency and magnitude taken the eigen-	**
** values contained in D.						**
**									**
**            w  = Im(D )                      R = Re(D )		**
**             i       i                        i      i		**
**									**
** Most of the intensities are zero because the operator M acts	as a	**
** filter to select out specific properties of the system & the system	**
** itself is often prepared in such a way as to zero particular comp- 	**
** onents. Therefore the sum will likely not span the full Liouville	**
** space, and we can restrict it to sum only over elements which have	**
** appreciable intensities, i.e. elements with norm(I ) > cutoff.	**
**                                                   p			**
**									**
**                                 <ls					**
**                                 ---					**
**                   t             \ 			            	**
**                 <M |sigma(t)> = /   I  * exp(-D*t)   + K		**
**                                 ---	 p         p			**
**                                  p 					**
**									**
** Assuming expectation values correspond to experimentally measured 	**
** data, the collective values form a time domain spectrum, and this 	**
** relates to a frequency domain spectrum through a Fourier transform.	**
**									**
**                     t						**
**            S(t) = <M |sigma(t)>         S(w) = FT{S(t)}		**
**									**
** It is assumed that the Fourier integral must be performed on values	**
** ranging from time t=0 to infinity, and that our data decays to zero	**
** by the time we are through measuring	except for constant components.	**
** Since the Fourier transform	is defined over (-inf,inf), the actual 	**
** transform taken is that of a step function multiplied into the 	**
** values we "could" have taken at t<0 but are not using.  Rather than	**
** use a discrete Fourier transform, we can theoretically use a formal	**
** one and the transformation outcome well known.			**
**									**
**               inf             inf [ <ls		      ]		**
**              /                /   | ---                  |		**
**              |       -iwt	   |   | \                    | -iwt	**
**       S(w) = | S(t)e     dt = |   | /   I exp(-D *t) + K |e    dt	**
**              /                /   | ---  p      p	      |		**
**              0                0   [  p 		      ]		**
**									**
**									**
**              inf    [ R  + i(w-w ) ]					**
**              ---    |  p	     p	|     [          i ]		**
**              \      | ------------ | + K | del(w) - - |		**
**       S(w) = /   I  |  2         2 |     [          w ]		**
**              ---  p | R  + (w-w )  |					**
**               p 	 [  p	    p   ]				**
**									**
*************************************************************************/

/*************************************************************************
**									**
** Additional Comments:							**
**									**
** Under certain approximations (e.g. liquid NMR ignoring relaxation)	**
** density operator evolution is vastly simplified, and the Liouville	**
** equation takes the form						**
**									**
**      	      |sigma(t)> = exp(-Lt)|sigma>			**
**									**
** L is still the active Liouvillian but the transformation is unitary.	**
** The superoperator exp(-Lt) where the Liouvillian is the commutation 	**
** superoperator of H, L = H (x) E - E (x) H, is then equivalent to a 	**
** unitary transformation superoperator	of exp(-iHt), which is given by	**
** exp(-iHt) (X) exp(+iHt).  System evolution in this case could be 	**
** written in Hilbert space without superoperators.  With this app-	**
** roximation								**
**                   t               t					**
**                 <M |sigma(t)> = <M |exp(-Lt)|sigma>			**
**									**
** and in L will be diagonal if constructed in the eigenbasis of H	**
**                                   -1					**
**	                L = S * D * S   = D				**
**									**
** The equation becomes							**
**									**
**                   t               t					**
**                 <M |sigma(t)> = <M |exp(-Dt)|sigma>			**
**									**
*************************************************************************/

#endif
