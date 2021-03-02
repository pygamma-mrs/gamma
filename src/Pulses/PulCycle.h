/* PulCycle.h ***************************************************-*-c++-*-
**			         					**
** 	                        G A M M A				**
**									**
**	Class Pulse Cycle			   Interface		**
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
** This class sets up a pulse cycle for GAMMA.  A pulse	cycle consists	**
** of a defined composite pulse (waveform) which is repeated (cycled) 	**
** over	a specified number of steps with a defined phase.  The class	**
** PulComposite handles aspects of the waveform and its propgators. 	**
** This contains a composite pulse, and tracks the pulse cycle and	**
** handles the propagators for spin system evolution as required in NMR	**
** simulations.								**
**									**
** Note: Because this module contains functions that plot the pulse	**
**       cycle in Gnuplot and FrameMaker, this will depend upon the	**
**       GamIO library.							**
**									**
*************************************************************************/

#ifndef   GPulCycle_h_			// Is this file already included?
#  define GPulCycle_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
# include <Pulses/PulComposite.h>	// Know about composite pulses
# include <HSLib/HSprop.h>		// Know about Hilbert propagators
# include <LSLib/LSprop.h>		// Know about Liouville propagators
# include <Matrix/row_vector.h>		// Know about row vectors
# include <string>			// Know about libstdc++ strings

class PulCycle : public PulComposite
  {
/*
  std::string      WFname;   INHERITED       // Pulse waveform name
  int         WFsteps;  INHERITED       // Pulse waveform # of steps
  row_vector  WFvals;   INHERITED       // Pulse waveform step {gamB1s,phases}
  row_vector  WFtimes;  INHERITED       // Pulse waveform step lengths (sec)
  double      WFtp;     INHERITED       // Pulse waveform total length (sec)
  int         rad;      INHERITED       // Flag for step phases in radians 
  std::string      Iso;	INHERITED	// Composite pulse channel
  gen_op*     Hsteps;	INHERITED	// Effective Hamiltonians
  int*        Hindex;	INHERITED	// Hamiltonian indexing
  int         Hcount;	INHERITED	// Hamiltonian count
  HSprop*     Usteps;	INHERITED	// Composite pulse step propagators
  HSprop*     Usums;	INHERITED	// Composite pulse summed propagators
  int         Ucount;	INHERITED	// Step propagator count
  LSprop*     Gsteps;	INHERITED	// Pulse sequence step superpropagators
  gen_op      Fzed;	INHERITED	// Fz operator, for z-axis rotations
  super_op    R;	INHERITED	// Relaxation/exchange superoperator
  densop      sigmaeq;	INHERITED	// Equilibrium density operator      */

  std::string      CYCname;			// Pulse cycle name
  int         CYCsteps;			// Pulse cycle # of steps
  row_vector  CYCvals;			// Pulse cycle step phase (deg)
  double      CYCtp;			// Pulse cycle length (sec)
  int*        Pindex;			// Pulse cycle propagator indexing
  int         Pcount;			// Pulse cycle propagator count
  int         CUcount;			// Pulse cycle propagator count
  HSprop*     CUsteps;			// Pulse cycle step propagators
  HSprop*     CUsums;			// Pulse cycle summed propagators
  int         CGcount;			// Pulse cycle super propagator count
  LSprop*     CGsteps;			// Pulse cycle step superprops
  LSprop*     CGsums;			// Pulse cycle summed superprops

private:

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                CLASS PULSE CYCLE ERROR HANDLING
// ____________________________________________________________________________


void CYCerror(int eidx, int noret=0) const;

        // Input		CYC	: Pulse Cycle (this)
	//			eidx 	: Error flag
	//			noret	: Return flag
        // Output               none	: Error Message Output


void volatile CYCfatality(int error) const;

        // Input		CYC	: Pulse Cycle (this)
	//			eidx 	: Error flag
        // Output               none	: Stops execution & error Message

// ____________________________________________________________________________
// ii            CLASS PULSE CYCLE DESTRUCTION FACILITATORS 
// ____________________________________________________________________________


void deleteCIndxs();

        // Input        CYC     : A pulse cycle (this)
        // Output       none    : CYC propagator indices
        //                        are deleted if they exist


void deleteCUprops();
 
        // Input        CYC     : A pulse cycle (this)
        // Output       none    : CYC step Hilbert propagators
        //                        are deleted if they exist
 

void deleteCGprops();
 
        // Input        CYC     : A pulse cycle (this)
        // Output       none    : CYC step superoperator propagators
        //                        are deleted if they exist


// ____________________________________________________________________________
// iii             CLASS PULSE CYCLE COPY FACILITATORS 
// ____________________________________________________________________________


void copyCIndxs(const PulCycle& CYC1);
 
        // Input        CYC     : A pulse cycle (this)
        //              CYC1    : A second pulse cycle
        // Output       none    : CYC1 propagator indices are
        //                        are copied from CYC1
        // Note                 : Assumed array Pindex is currently empty


void copyCUprops(const PulCycle& CYC1);
 
        // Input        CYC     : A pulse cycle (this)
        //              CYC1    : A second pulse cycle
        // Output       none    : CYC step Hilbert propagators
        //                        are copied from CYC1
        // Note                 : Assumed array Usteps is empty
 
 
void copyCGprops(const PulCycle& CYC1);
 
        // Input        CYC     : A pulse cycle (this)
        //              CYC1    : A second pulse cycle
        // Output       none    : CYC step Liouville propagators
        //                        are copied from CYC1
        // Note                 : Assumed array Gsteps is empty
 
 
// ____________________________________________________________________________
// iv        CLASS PULSE CYCLE PROPAGATOR INDEX GENERATORS
// ____________________________________________________________________________


void SetCIndxs( );

        // Input        CYC     : A pulse cycle (this)
        // Output       none    : CYC1 propagator indices are
        //                        are set up
        // Note                 : Assumed array Pindex is currently empty



// ____________________________________________________________________________
// v     CLASS PULSE CYCLE HILBERT SPACE PROPAGATOR GENERATORS
// ____________________________________________________________________________

/* These functions allow for the filling up two arrays of propagators for each
   step in the waveform (Usteps and Usums).  These will usually be generated
   during construction of the pulse cycle, they are system dependent.  */


void SetCUs();

        // Input                CYC     : A pulse cycle (this)
        // Output               void    : Propagators for each step of the
        //                                pulse cycle are calculated


// ____________________________________________________________________________
// vi    CLASS PULSE CYCLE LIOUVILLE SPACE PROPAGATOR GENERATORS
// ____________________________________________________________________________

/* These functions allow for the filling up two arrays of propagators for each
   step in the cycle (Gsteps and Gsums).  These are typically NOT generated
   during construction of the pulse cycle but may be generated at any
   time using the combination of the composite pulse & pulse cycle.  Note that
   these are indeed system dependent.                                        */

void SetCGs( );
 
        // Input                CYC     : A pulse cycle (this)
        //                      PW      : A pulse waveform
        // Output               void    : Superpropagators for each step of
        //                                the pulse cycle are calculated
        // Note                         : This demands that the Liouville
        //                                propagators have already been made
        //                                in the waveform.

     
 
// ____________________________________________________________________________
// vii                  CLASS PULSE CYCLE BASIS HANDLING
// ____________________________________________________________________________

/* These functions are only active on propagators.  They will have no effect 
   if there have been no propagators generated.                              */
 
void SetBasis(gen_op& Op);
     
        // Input                CYC     : A pulse cycle (this)
        // Output               void    : Hilbert space propagators for all
        //                                steps of the pulse cycle and its
        //                                composite pulse are put
        //                                into the current working basis of Op
 

 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 public:

// ____________________________________________________________________________
// A              CLASS PULSE CYCLE CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///Center Pulse Cycle Algebraic


MSVCDLC PulCycle();

	// Input	none		:
	// Output	CYC		: A NULL pulse cycle (this)
	///F_list 	PulCycle	- Constructor

 
MSVCDLC PulCycle(const PulComposite& P, const row_vector& S, std::string N);
 
        // Input        P       : Pulse cycle step composite pulse
        //              S       : Pulse cycle cycle steps (phase)
        //              N       : Pulse cycle cycle name
        // Output       CYC     : A new pulse cycle (this) 

 
MSVCDLC PulCycle(const PulCycle& PT1);

	// Input	CYC1	: Pulse cycle cycle
	// None		CYC	: Pulse cycle cycle (this) made from CYC1


// ------------------------------ Destruction ---------------------------------


MSVCDLC ~PulCycle();

	// Input	CYC	: A pulse cycle (this)
	// Output	none	: CYC is destructed


// ------------------------------- Assignment ---------------------------------
 

MSVCDLL PulCycle& operator = (const PulCycle& CYC1);

        // Input        CYC1	: Pulse cycle cycle
        // None         CYC	: Pulse cycle cycle (this) copied from CYC1
        ///F_list = 		- Assignment


// ____________________________________________________________________________
// B               CLASS PULSE CYCLE HAMILTONIAN FUNCTIONS
// ____________________________________________________________________________

 
/* These functions allow accessing average Hamiltonians over the 
   cycle steps.  They are generated from an input pulse waveform.  It is the
   waveform which actually contains individual step Hamiltonians.           */
 
//gen_op GetH(PulComposite& PC, int i=-1) const
 
        // Input                CYC    : A pulse cycle (this)
        //                      i     : Step in pulse cycle
        // Output               H     : The effective Hamiltonian
        //                              for pulse cycle step i

// ____________________________________________________________________________
// C          CLASS PULSE CYCLE HILBERT SPACE PROPAGATOR FUNCTIONS
// ____________________________________________________________________________

/* These functions allow access and production of propagators over the pulse
   train cycle steps.  They are generated from an input pulse waveform.  It
   is the waveform which actually contains individual step propagators.      */

MSVCDLL virtual HSprop GetCU(int i=-1);

	// Input                CYC     : A pulse cycle (this)
        //                      i       : Step in pulse cycle
        // Output               U       : The propagator for pulse cycle
        //                                that evolves from start of step i
        //                                             to the end of step i
        // Note                         : An exception is made for default
        //                                then the return will be the
        //                                propagator for the entire pulse
        //                                waveform
        // Note                         : If requested and not existing,
        //                                this will trigger ALL propators
        //                                to be generated if possible

 
MSVCDLL LSprop GetCG(int i, int j, double td);
 
        // Input                CYC     : A pulse cycle (this)
        //                      i       : Pulse cycle step index
        //                                where each step is a waveform
        //                      j       : Waveform step index
        //                      td      : Waveform step evolution time
        // Output               G       : Superpropagator that evolves through
        //                                WF step j of pulse cycle step i
        //                                for a time ot td 
        // Note                         : Indices are:  j=[0,WFsteps-1] 
        //                                              i=[0,CYCsteps-1] 

 
MSVCDLL virtual HSprop GetCUsum(int i=-1);
 
        // Input                CYC     : A pulse cycle (this)
        //                      i       : Step in pulse cycle
        //                                Default is last propagator -
        //                                i.e. for the entire pulse cycle
        // Output               U       : The propagator for pulse cycle
        //                                that evolves from the start of
        //                                waveform to the end of step i
        // Note                         : If requested and not existing,
        //                                this will trigger ALL propators
        //                                to be generated if possible


MSVCDLL HSprop GetCUsum(int i, int j);

        // Input                CYC     : A pulse cycle (this)
        //                      i       : Pulse cycle step index
        //                                where each step is a waveform
        //                      j       : Waveform step index
        // Output               U       : Propagator that evolves through
        //                                j WF steps of pulse cycle step i



MSVCDLL virtual HSprop GetCUmult(int N);

        // Input                CYC     : A pulse cycle (this)
        //                      N       : Number of pulse cycles
        // Output               U       : The propagator for N  full
        //                                cycles applied in succession


//virtual HSprop GetCU(double td);

        // Input                CYC     : A pulse cycle (this)
        //                      td      : Evolution time
        // Output               U       : The propagator that evolves a system
        //                                for time td under the pulse cycle


 
// ____________________________________________________________________________
// D        CLASS PULSE CYCLE SUPEROPERATOR PROPAGATOR FUNCTIONS
// ____________________________________________________________________________


MSVCDLL LSprop GetCG(int i=-1);
 
        // Input                CYC   : A pulse cycle (this)
        //                      i     : Step in pulse cycle
        //                              Default is last superpropagator -
        //                              i.e. for the entire pulse cycle
        // Output               G     : Superpropagator for pulse cycle
        //                              that evolves from the start of
        //                              the train to the end of step i


MSVCDLL LSprop GetCGsum(int i=-1);

        // Input                CYC     : A pulse cycle (this)
        //                      i       : Steps in pulse cycle
        //                                where each step is a waveform
        // Output               G       : Superpropagator that evolves through
        //                                i steps of the pulse cycle
        // Note                         : If requested and not existing,
        //                                this will trigger ALL propators
        //                                to be generated if possible
        // Note                         : Default propagator is for full
        //                                cycle (last in Gsums array)


MSVCDLL LSprop GetCGsum(int i, int j);

        // Input                CYC     : A pulse cycle (this)
        //                      i       : Pulse cycle step index
        //                                where each step is a waveform
        //                      j       : Waveform step index
        // Output               G       : Superpropagator that evolves through
        //                                j WF steps of pulse cycle step i
 
 
MSVCDLL LSprop GetCGmult(int N);
 
        // Input                CYC     : A pulse cycle (this)
        //                      N       : Number of pulse cycles
        // Output               G       : The superpropagator for N
        //                                cycles applied in succession
 

// ____________________________________________________________________________
// E               CLASS PULSE CYCLE ACCESS FUNCTIONS
// ____________________________________________________________________________
 
// ------------------------ Functions Over Full Cycle -------------------------

MSVCDLL int          steps()		const;
MSVCDLL int          WF_steps()		const;
MSVCDLL std::string  name()		const;
MSVCDLL std::string  WF_name()		const;
MSVCDLL row_vector   values()		const;
MSVCDLL row_vector   WF_values()	const;
MSVCDLL double       length()		const;
MSVCDLL double       WF_length()	const;

	// Input	CYC	: A pulse cycle (this)
        // Output       steps   : CYC steps
        //              name    : CYC name
        //              length  : CYC length (sec)
	//		values  : Array of phi values

// ------------------ Functions For Specific Cycle Step -----------------------

MSVCDLL complex value(int i)  const;
MSVCDLL double  phase(int i)  const;
 
        // Input        CYC     : A pulse cycle (this)
        // Output       phase   : Step phase value (degrees)

// ____________________________________________________________________________
// E               CLASS PULSE CYCLE AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 

MSVCDLL double steps(double td) const;
 
        // Input        CYC     : A pulse cycle (this) 
        //              td      : An evolution time (sec) 
        // Output       steps   : Number of cycle steps needed
        //                        to evolve for time td 


//int fullsteps(double td) const;

        // Input        CYC     : A pulse cycle (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of full cycle steps
        //                        that can occur in the time td
        // Note                 : For negative time we return
        //                        the total number of steps

 
MSVCDLL double cycles(double td) const;
 
        // Input        CYC     : A pulse cycle (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of cycles needed
        //                        to evolve for time td

 
MSVCDLL int fullcycles(double td=-1) const;
 
        // Input        CYC     : A pulse cycle (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of full cycles that can
        //                        occur in the time td
	// Note			: default is total # of steps 



MSVCDLL void scalegB1(double sf);

        // Input        CYC     : A pulse cycle (this)
        //              sf      : A scaling factor
        // Output       void    : All step field strengths in the waveform
        //                        are multiplied by sf.  The exception are
        //                        steps of zero length (ideal pulses)


// ____________________________________________________________________________
// F               CLASS PULSE CYCLE PLOTTING FUNCTIONS
// ____________________________________________________________________________

 
MSVCDLL row_vector IvsT(int split, int ends, int N=1) const;
 
        // Input                CYC     : A pulse cycle (this)
        //                      split   : Flag to split steps
        //                                 0: Don't split apart
        //                                 #: Split by #*10%*biggest step
        //                      ends    : Flag to add ends
        //                                 0: Don't put on ends
        //                                 #: Add ends #*.01*cycle length
        //                      N       : Number of cycles
        // Output               none    : Pulse Train Waveform plot vector
        //                                is made interactively of the
        //                                RF-intensity vs time.
        // Note                         : Ideal pulses are set by a step
        //                                having zero length but having
        //                                an non-zero "gamB1" setting.  In
        //                                such cases "gamB1" is the pulse angle
        // Note                         : For plotting purposes, the length of
        //                                an ideal 90 pulse will be taken to
        //                                be the same as the shortest non-zero
        //                                waveform step.  Others scale with
        //                                ideal pulse angle.
 
 
MSVCDLL row_vector PvsT(int split, int ends, int N=1, double ph=0) const;
 
        // Input                PWF     : Pulse Cycle Waveform
        //                      split   : Flag to split steps
        //                                 0: Don't split apart
        //                                 #: Split by #*.1*1st pulse length
        //                      ends    : Flag to add ends
        //                                 0: Don't put on ends
        //                                 #: Add ends length #*1st pulse
        //                      N       : Number of waveforms
        //                      ph      : An additional phase factor
        // Output               none    : Pulse Cycle Waveform plot vector
        //                                is made interactively of the
        //                                RF-phase vs time.



MSVCDLL void GP(int ty=1, int spl=0, int ed=0, int N=1, double p=0) const;

        // Input                CYC	: A pulse cycle (this)
	//                      ty	: Type of plot to create
        //                                 0: Waveform time vs phase
        //                                 1: Waveform time vs gamB1
        //                      spl	: Flag to split steps
        //                                 0: Don't split apart
        //                                 #: Split by #*.1*1st pulse length
        //                      ed	: Flag to add ends
        //                                 0: Don't put on ends
        //                                 #: Add ends length #*1st pulse
	//			N	: Number of cycles to plot
        // Output               none    : Pulse cycle cycle plot
	//			p	: Added phase (degrees)
        //                                is made interactively using
        //                                Gnuplot.

 
MSVCDLL void FM(int ty=1, int spl=0, int ed=0, int N=1, double p=0) const;
 
        // Input                CYC     : A pulse cycle (this)
	//                      ty	: Type of plot to create
        //                                 0: Waveform time vs phase
        //                                 1: Waveform time vs gamB1
        //                      spl	: Flag to split steps
        //                                 0: Don't split apart
        //                                 #: Split by #*.1*1st pulse length
        //                      ed	: Flag to add ends
        //                                 0: Don't put on ends
        //                                 #: Add ends length #*1st pulse
        //                      N       : Number of cycles to plot
	//			p	: Added phase (degrees)
        // Output               none    : Pulse cycle cycle plot
        //                                is made in FrameMaker MIF format
 
 
// ____________________________________________________________________________
// G               CLASS COMPOSITE PULSE EVOLUTION FUNCTIONS
// ____________________________________________________________________________
 
// ------------------------ Acquisiton Helper Functions 0----------------------
 
 
MSVCDLL double FIDsync(double& SW) const;
 
        // Input        CYC     : A pulse cycle (this)
        //              SW      : Desired spectral width
        // Output       SWsync  : Spectral width synchronized with the
        //                        pulse cycle length (at best) or
        //                        with the pulse step length (2nd best)
        // Note                 : If no synchronization possible the input
        //                        value is just returned.

 
MSVCDLL int FIDtest(double td, int& nCYs, int& nWFs, int& nSTPs, double& tr) const;
 
        // Input        CPul    : A composite pulse (this)
        //              td      : Desired dwell time (sec)
	//		nCYs	: Number of cycles
        //              nWFs    : Number of full waveforms
        //              nSTPs   : Number of full step
        //              tr      : Remainder time (sec)
        // Output       Fsync   : Synchronization flag
        //                             >1: Synchronized to Cycle
        //                              1: Synchronized to Waveform
        //                              0: Not synchronized
        //                              -: Synchronized to Step
        // Note                 : Values of nCYs, nWFs, nSTPs, and tr altered

 
// ------------------------ Without Relaxaton & Exchange ----------------------
 

MSVCDLL row_vector FIDsynchCYC(int npts, int nCYs,
                                       gen_op &D, gen_op& sigmap, int track=0);

        // Input        CYC     : A pulse cycle (this)
        //              npts    : Number of FID points
        //              nCYs    : Number of cycles between points
        //              D       : Detection operator
        //              sigmap  : Prepared density operator
        //              track   : Flag for tracking the computation
        // Output       data    : A row vector contiaing an npts point FID
        //                        that was generated by detection with D of
        //                        the evolution of sigmap under the composite
        //                        pulse CPul.


MSVCDLL row_vector FIDWFsynch(int npts, int nWFs,
                                       gen_op &D, gen_op& sigmap, int track=0);

        // Input        CPul    : A composite pulse (this)
        //              npts    : Number of FID points
        //              nWFs    : Number of waveforms between points
        //              D       : Detection operator
        //              sigmap  : Prepared density operator
        //              track   : Flag for tracking the computation
        // Output       data    : A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under composite
        //                        composite pulse CPul.

 
MSVCDLL row_vector FIDSTsynch(int npts, int nSTs,
                                       gen_op &D, gen_op& sigmap, int track=0);
 
        // Input        CPul    : A composite pulse (this)
        //              npts    : Number of FID points
        //              nSTs    : Number of steps between points
        //              D       : Detection operator
        //              sigmap  : Prepared density operator
        //              track   : Flag for tracking the computation
        // Output       data    : A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
        //                        composite pulse CPul.

 
MSVCDLL row_vector FID(int N, double td, gen_op &D, gen_op& sp, int F=0);
 
        // Input        CYC     : A pulse cycle (this)
        //              N       : Number of FID points
        //              td      : Dwell time between FID points
        //              D       : Detection operator
        //              sp      : Prepared density operator
        // Output       data    : A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
        //                        composite pulse CPul.
 
 
// ------------------------- With Relaxaton & Exchange ------------------------
 
MSVCDLL virtual row_vector FIDRsynchCYC(int npts, int nCYs,
                                       gen_op &D, gen_op& sigmap, int track=0);
 
        // Input        CYC     : A pulse cycle (this)
        //              npts    : Number of FID points
        //              nCYs    : Number of cycles between points
        //              D       : Detection operator
        //              sigmap  : Prepared density operator
        //              track   : Flag for tracking the computation
        // Output       data    : A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
        //                        composite pulse CPul.
        // Note                 : Assumes Relaxation Active!
 

MSVCDLL virtual row_vector FIDRWFsynch(int npts, int nWFs,
                                       gen_op &D, gen_op& sigmap, int track=0);
 
        // Input        CPul    : A composite pulse (this)
        //              npts    : Number of FID points
        //              nWFs    : Number of waveforms between points
        //              D       : Detection operator
        //              sigmap  : Prepared density operator
        //              track   : Flag for tracking the computation
        // Output       data    : A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
        //                        composite pulse CPul.
        // Note                 : Assumes Relaxation Active!


MSVCDLL virtual row_vector FIDRSTsynch(int npts, int nSTs,
                                       gen_op &D, gen_op& sigmap, int track=0);
 
        // Input        CPul    : A composite pulse (this)
        //              npts    : Number of FID points
        //              nSTs    : Number of steps between points
        //              D       : Detection operator
        //              sigmap  : Prepared density operator
        //              track   : Flag for tracking the computation
        // Output       data    : A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
        //                        composite pulse CPul.
        // Note                 : Assumes Relaxation Active!

 
MSVCDLL row_vector FIDR(int N, double td, gen_op &D, gen_op& sp, int F=0);

        // Input        CYC     : A pulse cycle (this)                    
        //              N       : Number of FID points
        //              td      : Dwell time between FID points
        //              D       : Detection operator
        //              sp      : Prepared density operator
        //              F       : Print evolution flag
        // Output       data    : A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
        //                        composite pulse CPul.
        // Note                 : Assumes Relaxation Active!


// ____________________________________________________________________________
// G                  CLASS PULSE CYCLE I/O FUNCTIONS
// ____________________________________________________________________________


MSVCDLL std::ostream& printEvolve(std::ostream &ostr, double td) const;
 
        // Input                CYC     : A pulse cycle (this)           
        //                      td      : Evolution time
        //                      ostr    : Output stream
        //                      full    : Flag for output amount
        // Output               none    : Pul. Train cycle evolution info
        //                                is sent to the output stream


MSVCDLL virtual std::ostream& printFID(std::ostream &ostr, double td, int npts) const;

        // Input        CYC     : A pulse cycle (this)
        //              ostr    : An output stream
        //              td      : Dwell time between FID points
        //              npts    : Number of FID points
        // Output               : Information regarding the FID generation
        //                        is set to the output stream ostr


MSVCDLL std::ostream& printSteps(std::ostream &ostr) const;

        // Input                CYC     : Pulse Cycle
        //                      ostr    : Output stream
        // Output               none    : Pulse Cycle steps are sent
        //                                to the output stream


MSVCDLL std::ostream& printInfo(std::ostream &ostr) const;

        // Input                CYC     : Pulse Cycle
        //                      ostr    : Output stream
        //                      full    : Flag for output amount
        // Output               none    : CYC storage info is sent
        //                                to the output stream


MSVCDLL std::ostream& printBase(std::ostream &ostr) const;

        // Input                CYC     : Pulse Cycle
        //                      ostr    : Output stream
        // Output               none    : Pulse Train base info is
        //                                sent to the output stream

 
MSVCDLL std::ostream& print(std::ostream &ostr, int full=0) const;

        // Input                CYC    : Pulse Cycle
        //                      ostr	: Output stream
        //                      full	: Flag for output amount
        // Output               none	: Pulse Cycle info is sent
        //                                to the output stream


MSVCDLL friend std::ostream &operator << (std::ostream &ostr, const PulCycle &CYC);

	// Input		CYC	: Pulse Cycle
        //                      ostr	: Output stream
        // Output               none	: Pulse cycle cycle waveform is sent
	//				  to the output stream
};

#endif
