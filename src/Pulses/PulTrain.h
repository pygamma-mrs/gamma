/* PulTrain.h ***************************************************-*-c++-*-
**			         					**
** 	                        G A M M A				**
**									**
**	Class Pulse Train                        Interface		**
**									**
**	Copyright (c) 1995			 			**
**	Dr. Scott A. Smith			 			**
**      Philippe Pelupessy						**
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
** This file contains a definition of Pulse Trains in GAMMA.  These 	**
** serve to facilitate application of pulse trains in the simulation 	**
** of NMR experiments.							**
**									**
** Each pulse train consists of a defined Waveform, Cycle, and Super-	**
** Cycle.  The Waveform is repeated with phase changes as specified by	**
** the Cycle and the Cycle is repeated with phase changes as specified	**
** by the SuperCycle.  Neither the Cycle or SuperCycle is mandatory.	**
** The Waveform (function and propagators) is handled by Composite	**
** Pulse class (PulComposite) which serves as the base class herein.	**
**									**
*************************************************************************/

#ifndef   GPulTrain_h_			// Is this file already included?
#  define GPulTrain_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Pulses/PulComposite.h>	// Include base class knowledge
#include <Pulses/PulTrain.h>		// Know about pulse cycles 
#include <Pulses/PulTrainSCyc.h>	// Know about pulse train supercycles 
#include <string>			// Know stdlibc++ std::strings

class PulTrain : public PulComposite
  {
  std::string Name;				// Pulse train name
  PulCycle PTCyc;			// Pulse train cycle
  PulTrainSCyc PTSCyc;			// Pulse train supercycle
  int Cycles;				// Flag if cycle defined
  int SCycles;				// Flag if supercycle defined

private:

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                      CLASS PULSE TRAIN ERROR HANDLING
// ____________________________________________________________________________


void PTerror(int eidx, int noret=0) const;

        // Input		PT	: Pulse Train (this)
	//			eidx 	: Error flag
	//			noret	: Return flag
        // Output               none	: Error Message Output


void volatile PTfatality(int error);

        // Input                error: Error flag
        // Output               none : Stops execution & error Message

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 public:

// ____________________________________________________________________________
//                   CLASS PULSE TRAIN CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///Center Pulse Train Algebraic


MSVCDLC PulTrain();

	// Input	none  	:
	// Output	PT	: A NULL PulTrain (this)
	///F_list 	PulTrain- Constructor

 
MSVCDLC PulTrain(const PulComposite& CPul, std::string N="");
 
        // Input        CPul    : Composite pulse
	//		N	: Composite pulse name
        // Output       PT      : Pulse train (this) constructed
        // Note                 : There will be no defined cycles  
        //                        or supercycles in the train 
 

MSVCDLC PulTrain(const PulComposite& CPul, const PulCycle& Cyc, std::string N="");

        // Input        CPul    : Composite pulse
	//		Cyc	: Pulse cycle
	//		N	: Composite pulse name
        // Output       PT	: Pulse train (this) constructed


MSVCDLC PulTrain(const PulComposite& CPul, const PulCycle& Cyc,
                                   const PulSupCycle& SCyc, std::string N="");

        // Input        CPul    : Composite pulse
        //              Cyc     : Pulse cycle
        //              SCyc	: Pulse supercycle
        //              N       : Pulse train name
        // Output       PT      : Pulse train (this) constructed


MSVCDLC PulTrain(const PulTrain& PT1);

	// Input	PT1  : Pulse train
	// None		PT   : Pulse train (this) constructed from PT1


MSVCDLC ~PulTrain();

	// Input	PT   : A PulTrain (this)
	// Output	none : PT is destructed


MSVCDLL PulTrain& operator = (const PulTrain& PT1);

        // Input        PT1  : Pulse train
        // None         PT   : Pulse train (this) copied from PT1
        ///F_list =           - Assignment



// ____________________________________________________________________________
// C                  CLASS PULSE TRAIN PROPAGATOR FUNCTIONS
// ____________________________________________________________________________

/* These functions allow access and production of propagators while the pulse
   train is active.  Remember, the waveform sets individual steps of the
   composite pulse.  It is cycled through with the phase changes specified
   in the pulse cycle.  The full cycle is cycled through with phase changes
   specified in the pulse supercycle.                                        */

MSVCDLL HSprop GetU(double td);

        // Input                PT      : A pulse train (this)
        //                      td      : A delay time(sec)
        // Output               U       : The propagator for evolving a system
        //                                under the pulse train for time td


// ____________________________________________________________________________
//                   CLASS PULSE TRAIN ACQUISITION FUNCTIONS
// ____________________________________________________________________________


MSVCDLL std::ostream& info(std::ostream& ostr, double td) const;

        // Input        PT    : An pultrain (this)
        //              td    : A delay time(sec)
        // Output       ostr  : Information concerning
        //                      the pulse train evolution are
        //                      written into the output stream



MSVCDLL std::ostream& info(std::ostream& ostr, double td, int npts) const;

        // Input        PT    : An pultrain (this)
        //              td    : A delay time(sec)
	//		npts  : Number of successive evolutions
        // Output       ostr  : Information concerning
        //                      the pulse train evolution are
        //                      written into the output stream
	// Note		      : Assumes one wishes npts
	//			successive evolutions over length td

// ____________________________________________________________________________
//                   CLASS PULSE TRAIN ACQUISITION FUNCTIONS
// ____________________________________________________________________________



MSVCDLL row_vector FID(int npts, double td, gen_op &D, gen_op& sigmap);
 
        // Input        PT    : An pultrain (this)
        //              npts  : Number of FID points
        //              td    : Dwell time between FID points
        //              D     : Detection operator
        //              sigmap: Prepared density operator
        // Output       data  : A row vector contiaing an npts point FID
        //                      that was generated by detection with D
        //                      of the evolution of sigmap under the pulse
        //                      train PT.
 
 
MSVCDLL row_vector FIDR(int npts, double td, gen_op &D, gen_op& sigmap);
 
        // Input        PT    : An pultrain (this)
        //              npts  : Number of FID points
        //              td    : Dwell time between FID points
        //              D     : Detection operator
        //              sigmap: Prepared density operator
        // Output       data  : A row vector contiaing an npts point FID
        //                      that was generated by detection with D
        //                      of the evolution of sigmap under the pulse
        //                      train PT.
        // Note               : Assumes Relaxation Active!
 

// ____________________________________________________________________________
//                      CLASS PULSE TRAIN I/O FUNCTIONS
// ____________________________________________________________________________


MSVCDLL std::ostream& printEvolve(std::ostream &ostr, double td, int full=0) const;

        // Input                PT      : Pulse Train
        //                      td      : Evolution time
        //                      ostr    : Output stream
        //                      full    : Flag for output amount
        // Output               none    : Pul. Train cycle evolution info
        //                                is sent to the output stream


MSVCDLL std::ostream& printCycle(std::ostream &ostr, int full=0) const;

        // Input                PT    : Pulse Train
        //                      ostr  : Output stream
        //                      full  : Flag for output amount
        // Output               none  : Pulse Train cycle info is sent
        //                              to the output stream

 
MSVCDLL std::ostream& printSCycle(std::ostream &ostr, int full=0) const;

        // Input                PT    : Pulse Train
        //                      ostr  : Output stream
        //                      full  : Flag for output amount
        // Output               none  : Pulse Train supercycle info is
        //                              sent to the output stream


MSVCDLL std::ostream& print(std::ostream &ostr, int full=0) const;

        // Input                PT    : Pulse Train
        //                      ostr  : Output stream
        //                      full  : Flag for output amount
        // Output               none  : Pulse Train info is sent
        //                              to the output stream


MSVCDLL friend std::ostream &operator << (std::ostream &ostr, const PulTrain &PT);

	// Input		PT   : Pulse train
        //                      ostr  : Output stream
        // Output               none  : Pulse train vector is sent
	//				to the output stream
 
};

#endif
