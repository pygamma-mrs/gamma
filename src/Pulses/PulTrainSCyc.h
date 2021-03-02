/* PulTrainSCyc.h ***********************************************-*-c++-*-
**			         					**
** 	                        G A M M A				**
**									**
**	Class Pulse Train Supercycle		   Interface		**
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
** This class sets up a supercycle for GAMMA's pulse train class.  A    **
** pulse cycle consists of a defined waveform which is repeated (cycled)**
** over a specified number of steps with a defined phase.  In turn, a   **
** pulse supercycle is an additonal loop which repeats the pulse cycle  **
** over a specified number of steps with a defined phase.  Information  **
** regarding the supercycle steps and phases is engulfed in class       **
** PulSupCycle, the class which serves as a base to this class. However **
** PulSupCycle has only minimal functionality and knows nothing about   **
** the applied pulse cycle nor the cycle step waveform.  This class     **
** also knows the pulse supercycle, but in addition handles the prop-   **
** agators for spin system evolution as required in NMR simulations.    **
** Note that it still has little to no knowledge of the applied cycle,  **
** the cycle step waveform, and is merely meant as a auxiliary class to **
** GAMMA's Pulse Train class (PulTrain) which will know specifics of    **
** the waveform, cycle, and supercycle including applicable Hamitonians **
** & propagators.  By and large it is class PulTrain that is used in    **
** GAMMA programs, not this class.                                      **
**                                                                      **
*************************************************************************/

#ifndef   GPulTrainSCyc_h_		// Is this file already included?
#  define GPulTrainSCyc_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Pulses/PulSupCycle.h>		// Know about the base clase
#include <Pulses/PulCycle.h>		// Know about pulse cycles
#include <HSLib/HSprop.h>		// Know about Hilbert propagators
#include <LSLib/LSprop.h>		// Know about Liouville propagators

class PulTrainSCyc : public PulSupCycle
  {
  double      tp;			// Pulse train cycle length (sec)
  int         Ucount;			// Pulse train propagator count
  HSprop*     Usteps;			// Pulse train cycle step propagators
  HSprop*     Usums;			// Pulse train cycle summed propagators
  LSprop*     Gsteps;			// Pulse sequence step superpropagators

private:

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                CLASS PULSE TRAIN SUPERCYCLE ERROR HANDLING
// ____________________________________________________________________________


void PTSCerror(int eidx, int noret=0) const;

        // Input		PTSC	: Pulse Train Cycle (this)
	//			eidx 	: Error flag
	//			noret	: Return flag
        // Output               none	: Error Message Output


void volatile PTSCfatality(int error) const;

        // Input		PTSC	: Pulse Train Cycle (this)
	//			eidx 	: Error flag
        // Output               none	: Stops execution & error Message

// ____________________________________________________________________________
// ii            CLASS PULSE TRAIN SUPERCYCLE DESTRUCTION FACILITATORS 
// ____________________________________________________________________________


void deleteHams();

        // Input        PTSC	: A pulse train cycle (this)
        // Output       none    : PTSC step Hilbert operators
        //                        are deleted if they exist


void deleteUprops();
 
        // Input        PTSC	: A pulse train cycle (this)
        // Output       none    : PTSC step Hilbert propagators
        //                        are deleted if they exist
 

void deleteGprops();
 
        // Input        PTSC     : A pulse train cycle (this)
        // Output       none    : PTSC step superoperator propagators
        //                        are deleted if they exist


// ____________________________________________________________________________
// iii             CLASS PULSE TRAIN SUPERCYCLE COPY FACILITATORS 
// ____________________________________________________________________________


void copyHams(const PulTrainSCyc& PTSC1);

        // Input        PTSC	: A pulse train cycle (this)
        //              PTSC1	: A second pulse train cycle
        // Output       none    : PTSC step Hamiltonians & their indices
        //                        are copied from PTSC1
        // Note                 : Assumed arrays Hsteps and Hindex
        //                        are currently empty


void copyUprops(const PulTrainSCyc& PTSC1);
 
        // Input        PTSC	: A pulse train cycle (this)
        //              PTSC1	: A second pulse train cycle
        // Output       none    : PTSC step Hilbert propagators
        //                        are copied from PTSC1
        // Note                 : Assumed array Usteps is empty
 
 
void copyGprops(const PulTrainSCyc& PTSC1);
 
        // Input        PTSC	: A pulse train cycle (this)
        //              PTSC1	: A second pulse train cycle
        // Output       none    : PTSC step Liouville propagators
        //                        are copied from PTSC1
        // Note                 : Assumed array Gsteps is empty
 

// ____________________________________________________________________________
// iv     CLASS PULSE TRAIN SUPERCYCLE HILBERT SPACE PROPAGATOR GENERATORS
// ____________________________________________________________________________

/* These functions allow for the filling up two arrays of propagators for each
   step in the waveform (Usteps and Usums).  These will usually be generated
   during construction of the pulse train cycle, they are system dependent.  */

void SetUs(PulCycle& PTC);

        // Input                PTSC	: A pulse train cycle (this)
        //                      PTC	: A pulse train cycle
        // Output               void    : Propagators for each step of the
        //                                pulse train cycle are calculated
        // Note                         : This assumes that the propagators
        //                                in PTC have already been calculated!
 
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 public:

// ____________________________________________________________________________
// A              CLASS PULSE TRAIN SUPERCYCLE CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///Center Pulse Train Cycle Algebraic


MSVCDLC PulTrainSCyc();

	// Input	none		:
	// Output	PTSC		: A NULL pulse train supercycle (this)
	///F_list 	PulTrainSCyc	- Constructor

 
MSVCDLC PulTrainSCyc(PulCycle& PTC, const row_vector& S, std::string N);
 
        // Input        PTC	: Pulse train cycle
        //              S       : Pulse train supercycle steps (phase)
        //              N       : Pulse train supercycle name
        // Output       PTSC    : A new pulse train supercycle (this) 

 
MSVCDLC PulTrainSCyc(PulCycle& PTC, const PulSupCycle& PCYC);
 
        // Input        PTC	: Pulse train cycle
        //              PCYC	: Pulse supercycle 
        // Output       PTSC	: A new pulse train supercycle (this)
 


MSVCDLC PulTrainSCyc(const PulTrainSCyc& PT1);

	// Input	PTSC1	: Pulse train supercycle
	// None		PTSC	: Pulse train supercycle (this) made from PTSC1


// ------------------------------ Destruction ---------------------------------


MSVCDLC ~PulTrainSCyc();

	// Input	PTSC	: A pulse train supercycle (this)
	// Output	none	: PTSC is destructed


// ------------------------------- Assignment ---------------------------------
 

MSVCDLL PulTrainSCyc& operator = (const PulTrainSCyc& PTSC1);

        // Input        PTSC1	: Pulse train supercycle
        // None         PTSC	: Pulse train supercycle (this) copied from PTSC1
        ///F_list = 		- Assignment


// ____________________________________________________________________________
// B               CLASS PULSE TRAIN SUPERCYCLE HAMILTONIAN FUNCTIONS
// ____________________________________________________________________________

 
/* These functions allow accessing average Hamiltonians over the pulse train
   cycle steps.  They are generated from an input pulse waveform.  It is the
   waveform which actually contains individual step Hamiltonians.           */
 
//gen_op GetH(PulCycle& PTC, int i=-1) const
 
        // Input                PTSC    : A pulse train supercycle (this)
        //                      i     : Step in pulse train supercycle
        // Output               H     : The effective Hamiltonian
        //                              for pulse train supercycle step i

// ____________________________________________________________________________
// C          CLASS PULSE TRAIN SUPERCYCLE HILBERT SPACE PROPAGATOR FUNCTIONS
// ____________________________________________________________________________

/* These functions allow access and production of propagators over the pulse
   train supercycle steps.  They are generated from an input pulse train cycle.
   It is the cycle which actually contains individual step propagators.      */ 

MSVCDLL HSprop GetU(int i=-1) const;

	// Input                PTSC	: A pulse train supercycle (this)
        //                      i       : Step in pulse train supercycle
        // Output               U       : The propagator for pulse train supercycle
        //                                that evolves from start of step i
        //                                             to the end of step i
        // Note                         : An exception is made for default
        //                                then the return will be the
        //                                propagator for the entire pulse
        //                                waveform
        // Note                         : If requested and not existing,
        //                                this will trigger ALL propators
        //                                to be generated if possible

 
MSVCDLL HSprop GetUsum(int i=-1) const;
 
        // Input                PTSC    : A pulse train supercycle (this)
        //                      i       : Steps in pulse train supercycle
        //                                where each step is a cycle
        // Output               U       : Propagator that evolves through
        //                                i steps of the pulse train supercycle
        // Note                         : If requested and not existing,
        //                                this will trigger ALL propators
        //                                to be generated if possible
        // Note                         : Default propagator is for full
        //                                supercycle (last in Usums array)

   
MSVCDLL HSprop GetUmult(int N) const;
 
        // Input                PTSC    : A pulse train supercycle (this)
        //                      N       : Number of pulse train supercycles
        // Output               U       : The propagator for N pulse train
        //                                supercycles applied in succession
 

 
// ____________________________________________________________________________
// D        CLASS PULSE TRAIN SUPERCYCLE SUPEROPERATOR PROPAGATOR FUNCTIONS
// ____________________________________________________________________________


MSVCDLL void SetGs(PulCycle& PTC);

        // Input                PTSC	: A pulse train supercycle (this)
        //                      PTC	: A pulse waveform
        // Output               void    : Superpropagators for each step of
        //                                the pulse train supercycle are calculated
        // Note                         : This demands that the Liouville
        //                                propagators have already been made
        //                                in the waveform.


MSVCDLL LSprop GetG(int i=-1) const;
 
        // Input                PTSC   : A pulse train supercycle (this)
        //                      i     : Step in pulse train supercycle
        //                              Default is last superpropagator -
        //                              i.e. for the entire pulse train supercycle
        // Output               G     : The superpropagator for pulse train supercycle
        //                              that evolves from the start of
        //                              the train to the end of step i

// ____________________________________________________________________________
// E               CLASS PULSE TRAIN SUPERCYCLE ACCESS FUNCTIONS
// ____________________________________________________________________________
 
// ------------------------ Functions Over Full Cycle -------------------------

//int        steps()   const;
//std::string     name()    const;
MSVCDLL double     length()  const;
//row_vector values()  const;

MSVCDLL int        steps()   const;

	// Input	PTSC	: A pulse train supercycle (this)
        // Output       steps   : PTSC steps
        //              name    : PTSC name
        //              length  : PTSC length (sec)
	//		values  : Array of phi values

// ------------------ Functions For Specific Cycle Step -----------------------

//double  phase(int i)    const;
 
        // Input        PTSC     : A pulse train supercycle (this)
        // Output       phase   : Step phase value (degrees)

 
MSVCDLL double steps(double td) const;
 
        // Input        PTSC    : A pulse train supercycle (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of supercycles steps needed
        //                        to evolve for time td
 


MSVCDLL int fullSCYCs(double td=-1) const;
MSVCDLL int fullsteps(double td=-1) const;

        // Input        PTSC    : A pulse train supercycle (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of supercycles steps needed
        //                        to evolve for time td

 

 
// ____________________________________________________________________________
// F                CLASS PULSE TRAIN SUPERCYCLE AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 

// ____________________________________________________________________________
// G          CLASS PULSE TRAIN SUPERCYCLE HAMILTONIAN & PROPAGATOR FUNCTIONS
// ____________________________________________________________________________


// ____________________________________________________________________________
// Z                 CLASS PULSE TRAIN SUPERCYCLE I/O FUNCTIONS
// ____________________________________________________________________________


MSVCDLL std::ostream& printInfo(std::ostream &ostr) const;

        // Input                PTSC     : Pulse Train Cycle
        //                      ostr    : Output stream
        //                      full    : Flag for output amount
        // Output               none    : PTSC storage info is sent
        //                                to the output stream



MSVCDLL std::ostream& printBase(std::ostream &ostr) const;

        // Input                PTSC     : Pulse Train Cycle
        //                      ostr    : Output stream
        // Output               none    : Pulse Train base info is
        //                                sent to the output stream

 
MSVCDLL std::ostream& print(std::ostream &ostr, int full=0) const;

        // Input                PTSC    : Pulse Train Cycle
        //                      ostr	: Output stream
        //                      full	: Flag for output amount
        // Output               none	: Pulse Train Cycle info is sent
        //                                to the output stream


MSVCDLL friend std::ostream &operator << (std::ostream &ostr, const PulTrainSCyc &PTSC);

	// Input		PTSC	: Pulse Train Cycle
        //                      ostr	: Output stream
        // Output               none	: Pulse train supercycle waveform is sent
	//				  to the output stream
};

#endif
