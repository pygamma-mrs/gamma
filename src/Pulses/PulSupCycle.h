/* PulSupCycle.h ************************************************-*-c++-*-
**			         					**
** 	                        G A M M A				**
**									**
**	Class Pulse Supercycle			   Interface		**
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
**									**
** This class sets up a pulse supercycle for GAMMA. A pulse supercycle	**
** consists of a defined pulse cycle (a cycled pulse waveform) which is	**
** repeated (i.e. cycled again) over a specified number of steps with	**
** defined phase chages.						**
**									**
*************************************************************************/

#ifndef   GPulSupCycle_h_		// Is this file already included?
#  define GPulSupCycle_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/row_vector.h>		// Include row vectors
#include <Pulses/PulWaveform.h>		// Include pulse waveforms
#include <Pulses/PulCycle.h>		// Include pulse cycles
#include <string>			// Include libstdc++ strings

class PulSupCycle
  {
  friend class PulTrainSCyc;		// Allow this class full access
  std::string       SCycname;		// Pulse supercycle name
  int          SCycnosteps;		// Pulse supercycle # of steps
  row_vector   SCycsteps;		// Pulse supercycle step phase (deg)

private:

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                CLASS PULSE SUPERSCYCLE ERROR HANDLING
// ____________________________________________________________________________


void SCycerror(int eidx, int noret=0) const;

        // Input		SCyc	: Pulse Supercycle (this)
	//			eidx 	: Error flag
	//			noret	: Return flag
        // Output               none	: Error Message Output


void volatile SCycfatality(int error) const;

        // Input		SCyc	: Pulse Supercycle (this)
	//			eidx 	: Error flag
        // Output               none	: Stops execution & error Message

 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 public:

// ____________________________________________________________________________
// A                 CLASS PULSE SUPERSCYCLE CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///Center Pulse Supercycle Algebraic


MSVCDLC PulSupCycle();

	// Input	none		:
	// Output	SCyc		: A NULL pulse supercycle (this)
	///F_list 	PulSupCycle	- Constructor

 
MSVCDLC PulSupCycle(const row_vector& ptsteps, const std::string& ptname="");

        // Input        ptsteps : Pulse Supercycle step points (phase)
        //              ptname  : Pulse Supercycle name
        // Output       SCyc     : A new pulse cycle (this)
        // Note                 : Pulse Supercycle will NOT contain
        //                        a working pulse waveform 

 
MSVCDLC PulSupCycle(const PulSupCycle& PT1);

	// Input	SCyc1	: Pulse Supercycle
	// None		SCyc	: Pulse Supercycle (this) made from SCyc1


MSVCDLC PulSupCycle(const PulCycle& Cyc);

	// Input	SCyc	: Pulse supercycle (this)
	//		Cyc	: Pulse cycle
	// None		SCyc	: Pulse Supercycle (this) made from SCyc1


// ------------------------------ Destruction ---------------------------------


MSVCDLC ~PulSupCycle();

	// Input	SCyc	: A pulse cycle (this)
	// Output	none	: SCyc is destructed


// ------------------------------- Assignment ---------------------------------
 

MSVCDLL PulSupCycle& operator = (const PulSupCycle& SCyc1);

        // Input        SCyc1	: Pulse Supercycle
        // None         SCyc	: Pulse Supercycle (this) copied from SCyc1
        ///F_list = 		- Assignment


MSVCDLL void operator = (const PulCycle& Cyc);

        // Input        SCyc    : Pulse Supercycle (this)
        //              Cyc     : Pulse Cycle
        // None         SCyc    : Pulse Supercycle copied from cycle

// ____________________________________________________________________________
// E               CLASS PULSE SUPERSCYCLE ACCESS FUNCTIONS
// ____________________________________________________________________________
 
// -------------------- Functions Over Full Supercycle ------------------------

MSVCDLL int        steps()   const;
MSVCDLL std::string     name()    const;
MSVCDLL row_vector values()  const;

	// Input	SCyc	: A pulse supercycle (this)
        // Output       steps   : SCyc steps
        //              name    : SCyc name
        //              length  : SCyc length (sec)
	//		values  : Array of phi values

// ----------------- Functions For Specific Supercycle Step -------------------
 
MSVCDLL complex value(int i)  const;
MSVCDLL double  phase(int i)  const;


// ____________________________________________________________________________
// F                CLASS PULSE SUPERSCYCLE AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 

// ____________________________________________________________________________
// D                 CLASS PULSE WAVEFORM PLOTTING FUNCTIONS
// ____________________________________________________________________________


//void GP(int split=0, int ends=0) const;

        // Input                SCyc    : Pulse Supercycle
        //                      split   : Flag to split steps
        //                                 0: Don't split apart
        //                                 #: Split by #*.1*1st pulse length
        //                      ends    : Flag to add ends
        //                                 0: Don't put on ends
        //                                 #: Add ends length #*1st pulse
        // Output               none    : Pulse Supercycle plot
        //                                is made interactively using
        //                                Gnuplot.


//void GP(const PulWaveform& PW, int split, int ends) const;

        // Input                SCyc     : Pulse Supercycle
        //                      PWF     : Pulse Waveform 
        //                      split   : Flag to split steps
        //                                 0: Don't split apart
        //                                 #: Split by #*step length/10
        //                      ends    : Flag to add ends
        //                                 0: Don't put on ends
        //                                 #: Ends, length #*cyclelength/10
        // Output               none    : Pulse Supercycle plot
        //                                is made interactively using
        //                                Gnuplot. Plot is phase vs time

// ____________________________________________________________________________
// Z                    CLASS PULSE SUPERSCYCLE I/O FUNCTIONS
// ____________________________________________________________________________


MSVCDLL std::ostream& printBase(std::ostream &ostr, double SCyclen=0) const;

        // Input                SCyc     : Pulse Supercycle
        //                      ostr    : Output stream
        //                      SCyclen  : Length of 1 cycle (sec)
        //                      full    : Flag for output amount
        // Output               none    : Pulse Supercycle base info is sent
        //                                to the output stream


MSVCDLL std::ostream& printSteps(std::ostream &ostr) const;
 
        // Input                SCyc     : Pulse Supercycle
        //                      ostr    : Output stream
        // Output               none    : Pulse Supercycle steps are sent
        //                                to the output stream


MSVCDLL std::ostream& print(std::ostream &ostr, int full=0) const;

        // Input                SCyc    : Pulse Supercycle
        //                      ostr	: Output stream
        //                      full	: Flag for output amount
        // Output               none	: Pulse Supercycle info is sent
        //                                to the output stream


MSVCDLL friend std::ostream &operator << (std::ostream &ostr, const PulSupCycle &SCyc);

	// Input		SCyc	: Pulse Supercycle
        //                      ostr	: Output stream
        // Output               none	: Pulse Supercycle waveform is sent
	//				  to the output stream
};

#endif
