/* PulWaveform.h ************************************************-*-c++-*-
**			         					**
** 	                        G A M M A				**
**									**
**	Class Pulse Waveform                          Interface		**
**									**
**	Copyright (c) 1998			 			**
**	Dr. Scott A. Smith			 			**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**                                                                      **
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**                                                                      **
** This class sets up a pulse waveform for implementation in GAMMA.     **
** The waveform in designated in steps, each step having a specified    **
** strength, rf-phase, and applied time.                                **
**                                                                      **
*************************************************************************/

#ifndef   GPulWaveform_h_		// Is this file already included?
#  define GPulWaveform_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/row_vector.h>		// Know GAMMA row vectors
#include <Matrix/complex.h>		// Know GAMMA complex numbers
#include <fstream>			// Know libstdc++ file streams
#include <string>			// Know libstdc++ strings

class PulWaveform
  {
  friend class PulComposite;		// Allow this class full access
  friend class PulCycle;		// Allow this class full access 
  friend class PulTrain;		// Allow this class full access
  std::string  WFname;			// Pulse waveform name
  int          WFsteps;			// Pulse waveform # of steps
  row_vector   WFvals;			// Pulse waveform step {gamB1s,phases}
  row_vector   WFtimes;			// Pulse waveform step lengths (sec)
  double       WFtp;			// Pulse waveform total length (sec)
  int          rad;			// Flag for step phases in radians

private:

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                CLASS PULSE WAVEFORM ERROR HANDLING
// ____________________________________________________________________________

        // Input		PT	: Pulse Waveform (this)
	//			eidx 	: Error flag
	//			noret	: Return flag
        // Output               none	: Error Message Output
	//				: Stops execution if fatal

         void PWFerror(int eidx, int noret=0) const;
volatile void PWFfatality(int eidx=0)         const;

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 public:

// ____________________________________________________________________________
// A              CLASS PULSE WAVEFORM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///Center Pulse Waveform Algebraic
///F_list 	PulWaveform	- Constructor
///F_list	=		- Assignment

        // Input        wfsteps : Pulse waveform step points (gB1,phase)
        //              wftimes : Pulse waveform step lengths
        //              wfname  : Pulse waveform name    
        //              wfiso   : Pulse waveform isotope channel    
        //              wfrad   : Pulse waveform radians flag
        // Output       PWF	: A new Pulse waveform (this)

MSVCDLC PulWaveform();
MSVCDLC PulWaveform(const row_vector& wfsteps, const row_vector& wftimes,
                         const std::string& wfname="", int wfrad=0);
MSVCDLC PulWaveform(const PulWaveform& PT1);


// ------------------------ Destruction & Assignment --------------------------

MSVCDLC virtual ~PulWaveform();
MSVCDLL PulWaveform&  operator = (const PulWaveform& PWF1);

// ____________________________________________________________________________
// B               CLASS PULSE WAVEFORM ACCESS FUNCTIONS
// ____________________________________________________________________________
 
// --------------------- Functions Over Full Waveform -------------------------

MSVCDLL virtual int         steps()   const;
MSVCDLL virtual std::string name()    const;
MSVCDLL virtual double      length()  const;
MSVCDLL virtual row_vector  values()  const;
MSVCDLL virtual row_vector  lengths() const;

	// Input	PWF	: A pulse waveform (this)
        // Output       steps   : PWF steps
        //              name    : PWF name
        //              length  : PWF length (sec)
	//		values  : Array of { gamB1, phi } values
	//		lengths : Array of times

// ----------------- Functions For Specific Waveform Step ---------------------

MSVCDLL         double  strength(int i) const;
MSVCDLL virtual double  phase(int    i) const;
MSVCDLL         double  length(int   i) const;
MSVCDLL virtual complex value(int    i) const;
 
        // Input        PWF     : A pulse waveform (this)
        // Output       strength: Step rf-field strength (Hz)
        //              phase   : Step phase value (degrees)
        //              value   : Step { gamB1, phi } values
        //              length  : Step length (sec)

 
// ____________________________________________________________________________
// C                CLASS PULSE WAVEFORM AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 
// ----------------------------- Step Lengths ---------------------------------
 
MSVCDLL double maxlength( ) const;
 
        // Input        PWF     : A pulse waveform (this)
        // Output       tsum    : Length from longest waveform step
 
 
MSVCDLL double minlength(double cutoff=1.e-13) const;
 
        // Input        PWF     : A pulse waveform (this)
        //              cutoff  : Minimum length to consider
        // Output       tsum    : Length from shortest waveform step


// ---------------------------- Step Strengths --------------------------------
 
 
MSVCDLL double maxgamB1( ) const;
 
        // Input        PWF     : A pulse waveform (this)
        // Output       gB1     : Strength of strongest waveform step
 
 
MSVCDLL double mingamB1( ) const;
 
        // Input        PWF     : A pulse waveform (this)
        // Output       gB1     : Strength of weakest waveform step
        // Output       tsum    : Length from shortest waveform step

// --------------------------- Variation Checks -------------------------------

        // Input        PWF    : A pulse waveform (this)
        // Output       T/F     : Returns true if the RF field strength
        //                        is constant through out the waveform
        // Output       T/F     : Returns true if the RF field phase
        //                        is constant through out the waveform
        // Output       T/F     : Returns true if the step time
        //                        is constant through out the waveform

MSVCDLL bool gamB1const() const;
MSVCDLL bool phaseconst() const;
MSVCDLL bool timeconst()  const;

// ---------------------------- Step Counting ---------------------------------

 
MSVCDLL double steps(double td) const;
  
        // Input        PWF     : A pulse waveform (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of waveform steps needed
        //                        to evolve for time td

 
MSVCDLL int fullsteps(double td) const;
 
        // Input        PWF     : A pulse waveform (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of full waveform steps
        //                        that can occur in the time td
 
// -------------------------- Waveform Counting -------------------------------

 
MSVCDLL double WFs(double td) const;
  
        // Input        PWF     : A pulse waveform (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of waveforms needed
        //                        to evolve for time td


MSVCDLL int fullWFs(double td, double cut=1.e-13) const;

        // Input        PWF     : Pulse Waveform
        //              td      : An evolution time (sec)
	//		cut     : Time zero cutoff (sec)
        // Output       steps   : Number of full waveforms that
        //                        are within time td


MSVCDLL double sumlength(int i) const;
  
        // Input        PWF     : A pulse waveform (this)
        //              i       : Step in waveform 
        // Output       tsum    : Length of first i waveform steps (sec)


// --------------------------- Waveform Scaling -------------------------------

MSVCDLL virtual void scalegB1(double sf);
 
        // Input        PWF     : A pulse waveform (this)
        //              sf      : A scaling factor
        // Output       void    : All step field strengths in PWF
        //                        are multiplied by sf.  The exception are
        //                        steps of zero length (ideal pulses)
 
  
// ____________________________________________________________________________
// D                 CLASS PULSE WAVEFORM PLOTTING FUNCTIONS
// ____________________________________________________________________________

//------------------------- Generic Plotting Functions ------------------------
 

MSVCDLL void getIdeal(double& gB1, double& ptt, int i) const;

        // Input                PWF     : Pulse Waveform
        //                      gB1     : Step intensity
        //                      ptt     : Step length
        //                      i       : Step index
        // Output               none    : gB1 and ptt are modified to
        //                                reflect how an ideal pulse should
        //                                be plotted for step i.
 

MSVCDLL row_vector IvsT(int split=0, int ends=0, int N=1) const;

        // Input                PWF     : Pulse Waveform
        //                      split   : Flag to split steps
        //                                 0: Don't split apart
        //                                 #: Split by #*.1*1st pulse length
        //                      ends    : Flag to add ends
        //                                 0: Don't put on ends
        //                                 #: Add ends length #*1st pulse
	//			N       : Number of waveforms
        // Output               none    : Pulse Waveform plot vector
        //                                is made interactively of the
        //                                RF-intensity vs time.


MSVCDLL row_vector PvsT(int spl=0, int ends=0, int N=1, double p=0) const;

        // Input                PWF     : Pulse Waveform
        //                      spl	: Flag to split steps
        //                                 0: Don't split apart
        //                                 #: Split by #*.1*1st pulse length
        //                      ends    : Flag to add ends
        //                                 0: Don't put on ends
        //                                 #: Add ends length #*1st pulse
	//			N       : Number of waveforms
        //                      p       : Added phase 
        // Output               none    : Pulse Waveform plot vector
        //                                is made interactively of the
        //                                RF-intensity vs time.


//------------------------- Gnuplot Plotting Functions ------------------------
 
 
MSVCDLL void GP(int type=0, int split=0, int ends=0, int N=1) const;

        // Input                PWF	: Pulse Waveform
        //                      type    : Type of plot to create
        //                                 0: Waveform time vs phase
        //                                 1: Waveform time vs gamB1
        //                      split   : Flag to split steps
        //                                 0: Don't split apart
        //                                 #: Split by #*.1*1st pulse length
        //                      ends    : Flag to add ends
        //                                 0: Don't put on ends
        //                                 #: Add ends length #*1st pulse
	//			N       : Number of waveforms
        // Output               none    : Pulse Waveform plot
        //                                is made interactively using
        //                                Gnuplot.


//---------------------- FrameMaker Plotting Functions ------------------------


MSVCDLL void FM(int type=0, int split=0, int ends=0, int N=1 ) const;

        // Input                PWF     : Pulse Waveform
        //                      type    : Type of plot to create
        //                                 0: Waveform time vs phase
        //                                 1: Waveform time vs gamB1
        //                      split   : Flag to split steps
        //                                 0: Don't split apart
        //                                 #: Split by #*.1*1st pulse length
        //                      ends    : Flag to add ends
        //                                 0: Don't put on ends
        //                                 #: Add ends length #*1st pulse
	//			N       : Number of waveforms
        // Output               none    : Pulse Waveform plot
        //                                is output in Framemaker MIF format

 
// ____________________________________________________________________________
// Z                 CLASS PULSE WAVEFORM I/O FUNCTIONS
// ____________________________________________________________________________


MSVCDLL std::ostream& printBase(std::ostream& ostr) const;
 
        // Input                PWF     : Pulse Waveform
        //                      ostr    : Output stream    
        // Output               none    : Pulse Waveform base parameters are
        //                                sent into the output stream
 

MSVCDLL std::ostream& printSteps(std::ostream& ostr, int full=0) const;
 
        // Input                PWF     : Pulse Waveform
        //                      ostr    : Output stream
        // Output               none    : Pulse Waveform steps are sent
        //                                to the output stream


MSVCDLL virtual std::ostream& print(std::ostream& ostr, int full=0) const;

        // Input                PWF	: Pulse Waveform
        //                      ostr	: Output stream
        //                      full	: Flag for output amount
        // Output               none	: Pulse Waveform info is sent
        //                                to the output stream


MSVCDLL friend std::ostream& operator << (std::ostream& ostr, const PulWaveform& PWF);

	// Input		PWF	: Pulse Waveform
        //                      ostr	: Output stream
        // Output               none	: Pulse waveform waveform is sent
	//				  to the output stream
};

#endif							// PulWaveform.h
