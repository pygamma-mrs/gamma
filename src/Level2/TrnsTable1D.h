/* TrnsTable1D.h ************************************************-*-c++-*-
**									**
** 	                       G A M M A				**
**									**
**	Class Transitions Table 1D		Interface		**
**									**
**	Copyright (c) 1999				 		**
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
** This file contains the interface of class TrnsTable1D.  The 		**
** purpose of this class is to embody a digitial representation of a	**
** 1D nmr/epr spectrum.  In this instance the spectrum is stored as a	**
** ntr x 2 complex array where ntr is the number of transtions.  Each	**
** row corresponds to a transition.  The first column contains both the	**
** transtition frequency (imaginaries, in 1/sec) and the transition	**
** linewidths (reals, 1/sec). The second column contains the transition	**
** intensities. The intensity is stored complex (a+ib --> Rexp[iphi])	**
** in order to accomodate to handle phased lines.			**
**									**
** This class may be used independently, but it also serves as an 	**
** auxiliary class to the GAMMA acquisition classes. Often GAMMA 	**
** programs deal with these table indirectly through other classes.	**
**									**
*************************************************************************/

#ifndef   TrnsTab1D_h_			// Is this file already included?
#  define TrnsTab1D_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Basics/ParamSet.h>		// Knowledge of GAMMA parameter sets
#include <Matrix/matrix.h>		// Knowledge of GAMMA matrices
#include <Matrix/complex.h>		// Knowledge of GAMMA complex numbers
#include <Matrix/row_vector.h>		// Knowledge of GAMMA row vectors
#include <string>			// Knowledge of stdlibc++ strings
#include <vector>			// Knowledge of stdlibc++ STL vectors
#include <list>

//forward definitions
class TTable1D;
MSVCDLL TTable1D sum(const TTable1D& TT1, const TTable1D& TT2, double res=1.e-6);
void offset(matrix& mx, double F, double LWR, int inHz=0);

class TTable1D: private matrix 
  {
         double ICUT;			// Spectrum intensity cutoff
         double INORM;			// Spectrum intensity normalization
         double PERCUT;			// Spectrum % contribution cutoff
         int    FRQTYPE;		// Frequency (output) type flag
         int    FRQSORT;		// Frequency (output) sort flag
         double FRQCONV;		// Frequency (output) conversion factor
         double SN;			// Signal to noise ratio
         int    HP, RP, LWP, T2P, PHP; 	// Various print flags
  static bool   FRQREV;			// Flag to reverse frequencies

  friend class acquire1D;		// Allow 1D acquisitions full access

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                CLASS TRANSITIONS TABLE 1D ERROR HANDLING
// ____________________________________________________________________________


        // Input        TTab1D: A transitions table (this)
	//		eidx  : An error index
	//		noret : Flag for carriage return
        // Output       void  : Error Message Output

 
        // Input        TTab1D: A transitions table (this)
        //              eidx    : Error index
        //              pname   : Additional error message
        //              noret   : Flag for linefeed (0=linefeed)
        // Output       void    : Error message output
     
/* The following error messages use the defaults set in the Gutils package
 
                Case                          Error Message
 
                (1)                     Problems With File PNAME
                (2)                     Cannot Read Parameter PNAME
                (3)                     Invalid Use Of Function PNAME
                (4)                     Use Of Deprecated Function PNAME
                (5)                     Please Use Class Member PNAME
                default                 Unknown Error - PNAME                */  


        // Input        TTab1D: A transitions table (this)
	//		eidx  : An error index
        // Output       void  : Error Message Output, Execution Stopped
 
         void TTaberror(int    eidx,                           int noret=0) const;
         void TTaberror(int    eidx, const std::string& pname, int noret=0) const;
volatile void TTabfatality(int eidx=0)                                      const;
volatile void TTabfatality(int eidx, const std::string& pname)              const;
 
// ____________________________________________________________________________
// ii              CLASS TRANSITIONS TABLE 1D CHECKING FUNCTIONS
// ____________________________________________________________________________

int CheckType(int typ, int warn=0) const;

	// Input	TTab1D: A transitions table (this)
	//		typ   : A frequency output type
	//				0 = Hertz (default)
	//				1 = PPM
	//				2 = Gauss
	//		warn  : A warning flag
	// Output	TF    : TRUE if value of typ is known

// ____________________________________________________________________________
// iii              CLASS TRANSITIONS TABLE COMMON CODE
// ____________________________________________________________________________

void setDefaults();

	// Input	TTab1D: A transitions table (this)
	// Output	void  : This sets the default internals


void copySettings(const TTable1D& TTab1);

	// Input	TTab1D: A transitions table (this)
	//		TTab1 : A second transitions table
	// Output	void  : Copies internal settings from TT1


// ____________________________________________________________________________
// iv                  TRANSITIONS TABLE 1D SETUP FUNCTIONS
// ____________________________________________________________________________
 
/* These functions set up specific aspects of a transitions table.  Since
   they make assumptions about the order in which the table is set up, the
   functions MUST be private, their misuse may make an inconsistent table.   */

int SetNTrans(const ParameterSet& pset, int warn=2);
 
        // Input                TTab1D  : A transitions table (this)
        //                      pset    : A parameter set
        //                                      0 = no warning
        //                                      1 = warning, non-fatal
        //                                     >1 = fatal
        // Output               none    : Table number of transitions is
        //                                set from parameters in pset


int SetConv(const ParameterSet& pset, int warn=0);

        // Input                TTab1D  : A transitions table (this)
        //                      pset    : A parameter set
        //                                      0 = no warning
        //                                      1 = warning, non-fatal
        //                                     >1 = fatal
        // Output               none    : Table conversion factor used in
        //                                PPM or Gauss units is set from
        //                                parameters in pset


double GetFreq(const ParameterSet& pset, int& idx, int warn=1);

        // Input                TTab1D  : A transitions table (this)
        //                      pset    : A parameter set
        //                      idx     : Transition index
        //                      warn    : Warning level
        //                                      0 = no warning
        //                                      1 = warning, non-fatal
        //                                     >1 = fatal
        // Output               W       : Table transition idx frequency
	// Note				: Return idx is -1 if failure


double GetLWR2T2(const ParameterSet& pset, int idx, int warn=0);

        // Input                TTab1D  : A transitions table (this)
        //                      pset    : A parameter set
        //                      idx     : Transition index
        //                      warn    : Warning level
        //                                      0 = no warning
        //                                      1 = warning, non-fatal
        //                                     >1 = fatal
        // Output               R       : Table transition idx rate


double GetPhase(const ParameterSet& pset, int idx, int warn=0);

        // Input                TTab1D  : A transitions table (this)
        //                      pset    : A parameter set
        //                      idx     : Transition index
        //                      warn    : Warning level
        //                                      0 = no warning
        //                                      1 = warning, non-fatal
        //                                     >1 = fatal
        // Output               Phi     : Table transition idx phase


double GetIntensity(const ParameterSet& pset, int idx, int warn=1);

        // Input                TTab1D  : A transitions table (this)
        //                      pset    : A parameter set
        //                      idx     : Transition index
        //                      warn    : Warning level
        //                                      0 = no warning
        //                                      1 = warning, non-fatal
        //                                     >1 = fatal
        // Output               Inorm   : Table transition idx intensity


int SetTrans(const ParameterSet& pset, int warn=2);

        // Input                TTab1D  : A transitions table (this)
        //                      pset    : A parameter set
        //                                      0 = no warning
        //                                      1 = warning, non-fatal
        //                                     >1 = fatal
        // Output               none    : Table of transitions is filled from
        //                                transition parameters in pset

// ____________________________________________________________________________
// v                TRANSITIONS TABLE 1D CUTOFF FUNCTIONS
// ____________________________________________________________________________

/* Consider the request for a time domain signal associated with a transitions
   table. If those transitions decay, their time domain signals will be zero
   after a certain amount of time, after which they need not be computed since
   they are effectively zero. The function ExponentialCutoffs looks at each
   transition and determines the last point to calculate, i.e. the point
   after which the transition will no longer generate any intensity to
   within a cutoff value. These indices are returned in a vector, on index
   per transition.

   If k is the last point in the exponential that should be calculated, and
   cutoff is a value (<0) afterwhich the intensity is essentially zero, then
   we know that
                          k >= -ln(cutoff)/[R*tinc]

   where tinc is the time between points and R is the transition decay rate.

                                                           */

std::vector<int> ExpCutoffs(double tinc, int N, double C) const;

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 public:

// ____________________________________________________________________________
// A                 CLASS TRANSITIONS TABLE 1D CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///Center Transitions table Algebraic
///F_list TTable1D     - Constructor
///F_list =           - Assignment

 
/* These constructors set up a new transitions table.  There are several ways
   to make an transitions tabel as listed below:

             Input Arguments                Resulting Operator
             ---------------       -----------------------------------------
                   -               Empty Table, No Transitions
                  mx               Table directly from nx2 complex array
                 TTab              Table replicating the input table         */

MSVCDLC      TTable1D();
MSVCDLC      TTable1D(const matrix& mx);
MSVCDLC      TTable1D(const matrix& mx, int warn);
MSVCDLC      TTable1D(const TTable1D& TTab1);
MSVCDLC      ~TTable1D();

MSVCDLL TTable1D& operator = (const TTable1D& TTab1);


// ____________________________________________________________________________ 
// B                     TRANSITION TABLE MANIPULATIONS
// ____________________________________________________________________________

// __________________________ Frequency Alterations ___________________________

/*             Function              Return Type & Value
               --------   -----------------------------------------------------
                center	   double:  The center frequency of the transitions
                offset     void:    Shift one or all transition frequencies
                FRscale    void:    Scale one or all transition frequencies
                  BC	   void:    Baseline correction by removal of 
				    frequencies near zero.

	   Input		TTab1D	: Transitions table (this)
                                wa      : Flag for weighted average (def=true)
	  			F   	: Frequency offset (Hz=default, 1/sec)
	  			inHz	: Flag whether F is in Hz or 1/sec
	  			tr	: Transition index
          			Fscf	: Frequency scaling factor
	  			res	: Frequency resolution
	   Output center        double  : The center frequency is returned
	   Output offset	void	: All frequency values or that of the
	  				  the specified transitions tr are
	  				  adjusted as Im(<i|mx|0>) += F 
	   Output FRscale	void	: All frequency values or that of the
	  				  the specified transitions tr are
	  				  adjusted as Im(<i|mx|0>) *= Fscf 
	   Output BC		void	: All transitions havign frequencies
					  within 0+/-res will be removed     */
 
MSVCDLL double center(bool wa=true);
MSVCDLL void   offset(double F, int inHz=1);
MSVCDLL void   offset(double F, int tr, int inHz);
MSVCDLL void   FRscale(double Fscf);
MSVCDLL void   FRscale(double Fscf, int tr);
MSVCDLL void   BC(double res=-1.0);

// __________________________ Intensity Alterations ___________________________

	// Input		TTab1D	: Transitions table (this)
        //			Iscf	: Intensity scaling factor
	// Output		void	: The transition intensities are all
	//				  adjusted as <i|mx|1> *= Iscf
 
MSVCDLL void Iscale(double         Iscf);
MSVCDLL void Iscale(double         Iscf, int tr);
MSVCDLL void Iscale(const complex& Iscf);
MSVCDLL void Iscale(const complex& Iscf, int tr);
MSVCDLL void Iremove(double        dcut=-1.0);

// __________________________ Linewidth Alterations ___________________________

	// Input		TTab1D	: Transitions table (this)
	//			LWR	: Linewidth/Rate adjustment (Hz)
	//			inHz	: Flag whether LWR is a linewidth (Hz)
	//				  or a rate (1/sec) <-- default
	// Output		void	: The linewidth values are all adjusted
	//				  as Re(<i|mx|0>) += F
	// Note				: Negative linewidths are not allowed,


MSVCDLL void broaden(double LWR, int inHz=1);
MSVCDLL void broaden(double LWR, int tr, int inHz);

// __________________________ Resolution Alterations __________________________

MSVCDLL void resolution(double res);

// ____________________________________________________________________________
// C                    Time Domain Spectra Generation
// ____________________________________________________________________________

/* These functions take the 1D transitions table and generate a time domain
   spectrum.  The user just inputs the number of points desired and the time
   increment to use between them.  There are a few other factors which will
   influence the spectra generated by these functions.

	ICUT	This is an intensity cutoff value.  Any transition whose
                intensity (norm) is below this value will not contribute to
                the output spectrum.
        PERCUT  This is a percentage cutoff value.  For each transition,
                points at times past where the transition intensity has
                decayed below PERCUT*Io will not be included.  The value of
		PERCUT is decimal %, i.e. PERCUT = [0, 1].                   */

	// Input	TTab1D: A transitions table (this)
	//		npts     : Number of points desired
        // or
	//		data  : 1D row vector for time interferrogram

	//		tinc  : Time increment between points (sec)
	// Output	data  : 1D data vector is filled with a 
	//			time domain spectrum.
	// Note		      : The input time increment is adjusted to
	//			account for "Hz" output.

MSVCDLL row_vector              T(int npts,         double tinc) const;
MSVCDLL void                    T(row_vector& data, double tinc) const;
MSVCDLL std::vector<row_vector> Ts(int npts,        double tinc) const;
MSVCDLL std::vector<int>        TCutoffs(int npts,  double tinc) const;

// ____________________________________________________________________________ 
// D                        Frequency Domain Spectra
// ____________________________________________________________________________

/* These functions take the 1D transitions table and generate a frequency
   domain spectrum.  The user just inputs the number of points desired and the
   spectral range. There are a few other factors which will influence the
   spectra generated by these functions.

	ICUT	This is an intensity cutoff value.  Any transition whose
                intensity (norm) is below this value will not contribute to
                the output spectrum.
        PERCUT  This is a percentage cutoff value.  For each transition,
                points at frequencies outside where the transition intensity
                has fallen below PERCUT*Io will not be included.  The value of
		PERCUT is decimal %, i.e. PERCUT = [0, 1].                   */

	// Input	TTab1D	: A transitions table (this)
	//		N     	: Number of points desired
	//		data	: 1D row vector for frequency spectrum
	//		Fst	: Starting frequency (Hz)
	//		Ffi	: Final frequency (Hz)
	//		Icut	: Intensity cutoff
	//		Per	: % below max. to ignore [0,1]
	// Output	data  	: Row vector filled with a frequency spectrum.
	// Note		      	: Value Per is a % in decimal form
	//			  of the maximum Lorentzian intensity
	//			  considered worth adding to the spectrum.
	// Note		      	: Each Lorentzian has intensity above
	//			  Icut*Imax between the points bound
	//			  by the indices ist and ifi.  Only those
	//			  points are used in generating the spectrum
	// Note		      	: To avoid poor resolution, if the frequency
	//			  between points, winc, does not fullfill the
	//			  relationship 5*winc < lwhh, then the
	//			  integrated Lorentzian intensities are used
	//		          rather than the point value.
	// Note			: Users must insure data is first zeroed!


MSVCDLL row_vector              F(int npts,        double Fst,double Ffi) const;
MSVCDLL void                    F(row_vector& data,double Fst,double Ffi) const;
MSVCDLL std::vector<row_vector> Fs(int npts,       double Fst,double Ffi) const;

// ____________________________________________________________________________ 
// E             Frequency Domain Spectra In Derivative Mode
// ____________________________________________________________________________

	// Input	TTab1D: A transitions table (this)
	// 		sigma : Density matrix (operator propagated)

	//		N     : Number of points desired
        // or
	//		data  : 1D row vector for frequency spectrum
	// or
	// 		mx    : Matrix of transition data

	//		fstart: Starting frequency (Hz)
	//		fend  : Final frequency (Hz)
	//		cutoff: Intensity cutoff
	//		cutper: % below max. to ignore [0,1]
	// Output	data  : Row vector filled with a
	//			frequency spectrum.
	// Note		      : Value cutper is a % in decimal form
	//			of the maximum Lorentzian intensity
	//			considered worth adding to the spectrum.
	// Note		      : Each Lorentzian has intensity above
	//			cutoff*Imax between the points bound
	//			by the indices ist and ifi.  Only those
	//			points are used in generating the spectrum
	// Note		      : To avoid poor resolution, if the frequency
	//			between points, winc, does not fullfill the
	//			relationship 5*winc < lwhh, then the
	//			integrated Lorentzian intensities are used
	//		        rather than the point value.
	// Note		      : Users must insure data is first zeroed!


MSVCDLL row_vector FD(int N,            double fstart, double fend) const;
MSVCDLL void       FD(row_vector& data, double fstart, double fend) const;

// ____________________________________________________________________________
// F                            PHASE CORRECTION
// ____________________________________________________________________________

/* These routines apply phase corrections to the transition table.           */


        // Input        TTab1D: A transitions table (this)
	//		w0    : Zero point pivot frequency (Hz)
	//		w1    : 1st order pivot frequency (Hz)
	//		order : Maximum allowed correction order 
	//           OR
	//		wpivot: 1st order pivot frequency (Hz)
	//		P     : Phase correction factors
	//			  Re(P) = 0th order p.c. in radians
	//			  Im(P) = 1st order time const. (sec)
	// Output             : The transition phases in table are
	//			all phased according to the correction
	//			specified in Wpivot and P
	// 		z     : A complex number whose real part
	//			is the zeroth order phase correction
	//			and whose imaginary is the 1st order
	//			phase correction
	//		w0    : Set to 1st order pivot (rad/sec)


/* The first routine phase corrects all the transitions in the table using a
   0th order and 1st order phase correction.  The 0th order finds the
   transition closest to w0, takes the phase of that transition, and then
   removes that phase component from all transitions.  The 1st order finds
   the transition closest to w1 and looks at its phase.  It then assumes that
   there is some constant tp such that

                             (w1-w0)*tp = phase(w1) 

   and calculates tp, noting that and addition of 360 degrees to the phase will
   suffice, i.e. that there is another factor, ntp, where

                               (w1-w0)*ntp = 360
   such that
                      (w1-w0)*(tp+n*ntp) = phase(w1) + n*360 

   where n is an integer.  Any value of n will produce a proper phase value
   of (w1-w0)*(tp+n*ntp) that would undo the phase at w1.  The only requirement
   is to select which value of n to use.  To to that, this routine picks the
   biggest transition, say Wtr, and sees which value of n minimizes

                      (Wtr-w0)*(tp+n*ntp) - phase(Wtr) -----> 0

   where n is kept between the range [0, order].  Once found, all transitions
   are phase corrected by adjusting their complex intensities according to

  		    Inew(W) = Iold(W) * exp[-i*(W-w0)*(tp+n*ntp)]

   --------------------------------------------------------------------------

   The second function phase corrects all the transitions in the table using a
   specified 0th order and 1st order phase correction as well as a pivot where
   the 1st order correctionis zero.  Using phi0 as the zeroth order phase
   [phi0 = Re(P)] and t1 as the 1st order phase correction time constant
   [t1 = Im(P)], the correction to a transition at frequency W is given by

  				    -i*phi0    -(W-Wpivot)*t1
  		I   (W) = I  (W) * e        * e
  		 out       in
   
   where Wpivot is the pivot frequency where no 1st order phase correction
   is applied.                                                               */

MSVCDLL complex pcorrect(double& w0,    double w1, int order);
MSVCDLL void    pcorrect(double Wpivot, complex& P);

// ____________________________________________________________________________ 
// G            CLASS TRANSITIONS TABLE 1D AUXILIARY FUNCTIONS
// ____________________________________________________________________________ 

/* These functions allow users to access individual elements of the transitions
   table.  They each take a transition index and will return 0 if the index
   provided is out of range.     

               Function              Return Type & Value
               --------   ------------------------------------------
                  R2       double:  The transition relaxation rate
                  FR       double:  The transition frequency
                  I        complex: The transition intensity 
                  Tr	   vector:  The transiton (in array format)          */

MSVCDLL double     R2(int tr) const;
MSVCDLL double     Fr(int tr) const;
MSVCDLL complex    I(int  tr) const;
MSVCDLL row_vector Tr(int tr) const;

/* These functions allow users to access values that are global over all
   transitons in the table.  No transition index is required.

        Function                   Return Type & Value
       -----------   ------------------------------------------------------
       LineWidths    int: TRUE if any transitions have non-zero linewidths
       Intensities   int: TRUE if any transitions have intensity >= ICUT
       Phases	     int: TRUE if any transtions have phases > .0001 rad
       size          int: Total number of transitions in the table
       FRmax         double: Largest transition frequency in the table
       FRmin         double: Smallest transition frequency in the table
       Tdmin         double: Smallest dwell time needed for time domain spec	
       LWmax         double: Largest transition linewidth in the table
       Imax          double: Largest transition intensity in the table
       Noisemax      double: Largest noise signal present in table spectrum
       mx            matrix: Table as a matrix                               */

MSVCDLL bool   LineWidths()  const;
MSVCDLL bool   Intensities() const;
MSVCDLL bool   Phases()      const;
MSVCDLL int    size()        const;
MSVCDLL double FRmax()       const;
MSVCDLL double FRmin()       const;
MSVCDLL double Tdmin()       const;
MSVCDLL double LWmax()       const;
MSVCDLL double LWmin()       const;
MSVCDLL double Imax()        const;
MSVCDLL double Noisemax()    const;
MSVCDLL matrix mx()	       const;

MSVCDLL friend TTable1D sum(const TTable1D& TT1, const TTable1D& TT2, double res);

	// Input	TT1	: A transitions table
	// 		TT2	: Another transitions table
        //              res	: Resolution (in radians/sec)
        // Output	TT	: Transitions table which is a blend of
	//			  the transitions from TT1 & TT2.  Transitions
	//			  within res frequency from one another are
	//			  blended into a single transition with a 
	//			  summed intensity and weighted frequency
 

MSVCDLL std::vector<int> Sort(int k, int type, int colf) const;

        // Input	TTab1D	: A transitions table (this)
 	//		k     	: A row/column index to sort
	//		type	: Basis for the sorting
	//			    0 = sort real values (default)
	//			   >0 = sort norms
	//		           <0 = sort imaginaries
	//		colf	: Row or column sort
	//			    0 = sort rows (default)
	//			   !0 = sort columns
	// Output	vector  : A vector of integers will be return that
	//			  contains indices sorted as specified
	// Note		      	: The matrix is left unaltered

// sosi - General for any matrix.  Put into class matrix when time permits.

  
// ____________________________________________________________________________
// H       CLASS TRANSITIONS TABLE 1D WITH PARAMETERS & PARAMETER SETS
// ____________________________________________________________________________

/* These functions allow transition tables to be built from GAMMA parameter
   sets. The overloaded form allows for direct reading of an external ASCII
   file or gleaning the table from a pre-existing parameter set.           
 
           Input                TTab1D  : A transitions table (this)
                                filein  : An input (ASCII) file
                                pset    : A parameter set
                                indx    : Table index
                                warn    : Warning output level
                                            0 = no warnings
                                            1 = non-fatal warnings
                                            2 = fatal warnings
           Output               none    : Transitions table filled with
                                          transitions specified in filein    */
 
MSVCDLL bool readPSet(const std::string& filein, int indx=-1, int warn=2);
MSVCDLL bool readPSet(const ParameterSet& pset,  int indx,    int warn=2);
  
// ____________________________________________________________________________
// H             CLASS TRANSITIONS TABLE 1D I/O FUNCTIONS
// ____________________________________________________________________________

/* Class TTable1D maintains both the transition frequencies and decay rates
   in units of 1/sec.  For output, these values are converted to appropriate
   units for display.  In particular, the half-height linewidths and T2
   relaxation times are related to the transition decay rate according to

                                     R     1
                           lwhh   = -- = -----
                               Hz   PI   PI*T
                                             2

   Several other internal flags will affect the appearance of the output 
   transitions table.  These may be set prior to output using the provided
   access functions.  The flags and access functions are listed below.

    Flag     Access Function               Values & Meaning
   -------   ---------------    ------------------------------------------
   FRQTYPE      setType         0:Hz, 1:PPM, 2:Gauss (output freq. units)
   FRQSORT      setSort         0:No sort, !0:Sort output by frequency
   FRQCONV      setConv         FRQTYPE=1: Om(MHz), FRQTYPE=2: g(unitless)
   ICUT         setIcut         #>0: Intensity cutoff value for output
   SN           setSN           #>=0: Signal/Noise ratio, 0==no noise
   HP           setHprint       !0=print header(default), 0=no header
   RP           setRprint       1: Flag print R values \ 
   LWP          setLWprint      1: Flag print linewidths\  <0 = never
   T2P          setT2print      1: Flag print T2 values /   0 = as needed
   PHP          setPHprint      1: Flag print phases   /    1 = always       */

MSVCDLL void setType(int     typ);
MSVCDLL void setSort(int     sf);
MSVCDLL void setConv(double  cf);
MSVCDLL void setIcut(double  ct);
MSVCDLL void setInorm(double in=0.0);
MSVCDLL void setSN(double    S2N);
MSVCDLL void setHprint(int   hp);
MSVCDLL void setRprint(int   rp);
MSVCDLL void setLWprint(int  lwp);
MSVCDLL void setT2print(int  t2p);
MSVCDLL void setPHprint(int  php);
MSVCDLL void setFreqRev();
 
MSVCDLL int    getType()    const;
MSVCDLL int    getSort()    const;
MSVCDLL double getConv()    const;
MSVCDLL double getIcut()    const;
MSVCDLL double getInorm()   const;
MSVCDLL double getSN()      const;
MSVCDLL int    getHprint()  const;
MSVCDLL int    getRprint()  const;
MSVCDLL int    getLWprint() const;
MSVCDLL int    getT2print() const;
MSVCDLL int    getPHprint() const;
MSVCDLL bool   getFreqRev() const;

// ---------------------- ASCII Output Functions ------------------------------

MSVCDLL std::vector<std::string> printStrings()           const;
MSVCDLL std::ostream&            print(std::ostream& out) const;

        // Input                out      : output stream;
        // Output               none     : modifies output stream
        ///F_list print                  - Send spin system to output stream.

MSVCDLL friend std::ostream& operator << (std::ostream &ostr, const TTable1D& TTab);

	// Input		TTab   : Transitions table
        //                      ostr  : Output stream
        // Output               none  : Transitions table vector is sent
	//				to the output stream

MSVCDLL std::ostream& status(std::ostream& ostr) const;

        // Input		TTab1D	: A transitions table (this)
        //                      ostr    : An output stream
        // Output               ostr    : The output stream modified by
        //                                a view of the table status


// ------------------ Binary Output Functions (For Storage) -------------------

	// Input		TTab1D	: Transitions table1D (this)
        //                      fn	: Filename
        //                      fp	: File stream (pointing at TTab1D)
        // Return               void	: The  table TTab1D is written
        //				  to a file called fn
	//				  or to an output file stream fp
	//				  starting at it current position
        // Note				: The file format is BINARY

MSVCDLL void           write(const std::string& fn) const;
MSVCDLL std::ofstream& write(std::ofstream&     fp) const;

MSVCDLL void       dbwrite_old( const std::string& fileName, 
                                const std::string& compname,  // metabolite name			 
                                const double& lowppm, 
                                const double& highppm, 
                                const double& specfreq, 
                                const double& reffreq,
                                const int& loop,
                                const std::vector<std::string> & header) const;

MSVCDLL void        dbwrite(    const std::string& fileName, 
                                const std::string& compname, // metabolite name.
                                const double& specfreq,
                                const int& numberspins,
                                const int& loop,
                                const std::vector<std::string> & header) const;


MSVCDLL unsigned int calc_spectra(  std::vector<double> & freq,
                                    std::vector<double> & ampl,
                                    std::vector<double> & phase,
                                    double specfreq,
                                    int numberspins,
                                    double freqtol = 0.001, // was 0.1/specfreq
                                    double phasetol = 45.0, // was 50.0
                                    double lowppm = -1000.0,
                                    double highppm = 1000.0) const;

MSVCDLL unsigned int calc_spectra2( std::vector<double> & freq,
                                    std::vector<double> & ampl,
                                    std::vector<double> & phase,
                                    double specfreq,
                                    int numberspins,
                                    double freqtol = 0.001, // was 0.1/specfreq
                                    double phasetol = 45.0, // was 50.0
                                    double lowppm = -1000.0,
                                    double highppm = 1000.0,
                                    double normal = 0.0,
                                    bool verbose = false) const;


// ------------------- Binary Read Functions (For Storage) --------------------
 
    // Input            TTab1D  : Transitions table1D (this)
    //                      fn  : Filename
    //                      fp  : File stream (pointing at TTab1D spot)
    // Return             void  : Table is read from file filename or
    //                            from filestream fp taken from its
    //                            current location
 
MSVCDLL void           read(const std::string& fn);
MSVCDLL std::ifstream& read(std::ifstream&     fp);

// ____________________________________________________________________________ 
// I         CLASS TRANSITIONS TABLE 1D DEBUGGING HELPER FUNCTIONS
// ____________________________________________________________________________ 

 
MSVCDLL std::ostream& printT(std::ostream& ostr, double tinc, int npts, int P2P=0);

    // Input        TTab1D  : Transitions table (this)
    //              tinc	: Dwell time in seconds
    //              npts	: Number of points
    //              P2P     : Computation control
    //                     0 = Evolve From t=0 To Point 
    //                    !0 = Evolve From Point To Point 
    // Output       ostr	: The output stream modified by
    //                                the acquisition time domain
    //                                calculation with debugging info
 
MSVCDLL std::ostream& printF(std::ostream& ostr, int npts, double Fst, double Ffi);

	// Input		TTab1D : Transitions table (this)
	//			npts	: Number of points
	//			Fst	: Initial frequency
	//			Ffi	: Final frequency
        // Output               ostr	: The output stream modified by
        //                                the acquisition frequency domain
	//				  calculation with debugging info

// ____________________________________________________________________________
// J         CLASS TRANSITIONS TABLE OLD GAMMA SUPPORT FUNCTIONS
// ____________________________________________________________________________

// sosi  this is deprecated.... ever used?

friend void offset(matrix& mx, double F, double LWR, int inHz);

  };

#endif							// TrnsTable1D.h
