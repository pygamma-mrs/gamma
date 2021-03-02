/* PulComposite.h ************************************************-*-c++-*-
**			         					**
** 	                        G A M M A				**
**									**
**	Class Composite Pulse 			     Interface		**
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
** This class handles composite pulses for use in GAMMA.   A composite	**
** pulse will contain a composite pulse (composite pulse steps) and can	**
** generate bot step and pulse Hamilonians & propagators.		**
**									**
*************************************************************************/

#ifndef   GPulComposite_h_		// Is this file already included?
#  define GPulComposite_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Pulses/PulWaveform.h>		// Include the header file
#include <HSLib/SpinSystem.h>		// Include isotropic systems
#include <HSLib/GenOp.h>		// Know about operators
#include <HSLib/HSprop.h>		// Know about Hilbert propagators
#include <LSLib/LSprop.h>		// Know about Liouville props
#include <LSLib/SuperOp.h>		// Know about superoperators
#include <LSLib/DensOp.h>		// Know about density operators
#include <string>			// Know about libstdc++ strings

class PulComposite : public PulWaveform
  {
  friend class PulCycle;		// Allow this class full access
  friend class PulTrain;		// Allow this class full access
/*
  std::string      WFname;	INHERITED	// Pulse waveform name
  int         WFsteps;	INHERITED	// Pulse waveform # of steps
  row_vector  WFvals;	INHERITED	// Pulse waveform step {gamB1s,phases}
  row_vector  WFtimes;	INHERITED	// Pulse waveform step lengths (sec)
  double      WFtp;	INHERITED	// Pulse waveform total length (sec)
  int         rad;	INHERITED	// Flag for step phases in radians   */
  std::string      Iso;			// Composite pulse channel
  gen_op*     Hsteps;			// Effective Hamiltonians
  int*        Hindex; 			// Hamiltonian indexing
  int         Hcount;			// Hamiltonian count
  HSprop*     Usteps;			// Composite pulse step propagators
  HSprop*     Usums;			// Composite pulse summed propagators
  int*        Uindex; 			// Propagator indexing
  int         Ucount;			// Step propagator count
  super_op*   Lsteps;			// Effective Hamiltonian superoperators
  LSprop*     Gsteps;			// Pulse sequence step superpropagators
  LSprop*     Gsums;			// Pulse sequence summed superprops
  densop*     SigSSs;			// Steady-state density operators
  gen_op      Fzed;			// Fz operator, for z-axis rotations
  super_op    R;			// Relaxation/exchange superoperator
  densop      sigmaeq;			// Equilibrium density operator
  double      cutzero;			// Zero length cutoff (precision)

private:

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                CLASS COMPOSITE PULSE ERROR HANDLING
// ____________________________________________________________________________


void CPulerror(int eidx, int noret=0) const;

        // Input		CPul	: Composite Pulse (this)
	//			eidx 	: Error flag
	//			noret	: Return flag
        // Output               none	: Error Message Output


void volatile CPulfatality(int error) const;

        // Input		CPul	: Composite Pulse (this)
	//			eidx 	: Error flag
        // Output               none	: Stops execution & error Message

// ____________________________________________________________________________
// ii            CLASS COMPOSITE PULSE DESTRUCTION FACILITATORS 
// ____________________________________________________________________________


void deleteHams();

        // Input        CPul	: A composite pulse (this)
        // Output       none	: CPul step Hilbert operators
        //                        are deleted if they exist


void deleteHIndxs();

        // Input        CPul    : A composite pulse (this)
        // Output       none    : CPul Hamiltonian indices
        //                        are deleted if they exist

 
void deleteUIndxs();
 
        // Input        CPul    : A composite pulse (this)
        // Output       none    : CPul propagator indices 
        //                        are deleted if they exist 
 

virtual void deleteUprops();
 
        // Input        CPul	: A composite pulse (this)
        // Output       none    : CPul step Hilbert propagators
        //                        are deleted if they exist
 

void deleteLOps();
 
        // Input        CPul    : A composite pulse (this)
        // Output       none    : CPul step Liouville space operators
        //                        are deleted if they exist
 

virtual void deleteGprops();
 
        // Input        CPul	: A composite pulse (this)
        // Output       none    : CPul step superoperator propagators
        //                        are deleted if they exist
 
 
void deleteSSs();
 
        // Input        CPul    : A composite pulse (this)
        // Output       none    : CPul step steady-state operators
        //                        are deleted if they exist
 

// ____________________________________________________________________________
// iii             CLASS COMPOSITE PULSE COPY FACILITATORS 
// ____________________________________________________________________________


void copyHams(const PulComposite& CPul1);

        // Input        CPul     : A composite pulse (this)
        //              CPul1    : A second composite pulse
        // Output       none    : CPul step Hamiltonians & their indices
        //                        are copied from CPul1
        // Note                 : Assumed arrays Hsteps and Hindex
        //                        are currently empty


void copyHIndxs(const PulComposite& CPul1);
 
        // Input        CPul    : A composite pulse (this)
        //              CPul1   : A second composite pulse
        // Output       none    : CPul1 Hamiltonian indices are
        //                        are copied from CPul1
        // Note                 : Assumed array Hindex is empty
 
 
void copyUIndxs(const PulComposite& CPul1);
 
        // Input        CPul    : A composite pulse (this)
        //              CPul1   : A second composite pulse
        // Output       none    : CPul1 propagator indices are
        //                        are copied from CPul1
        // Note                 : Assumed array Pindex is currently empty
 

virtual void copyUprops(const PulComposite& CPul1);
 
        // Input        CPul     : A composite pulse (this)
        //              CPul1    : A second composite pulse
        // Output       none    : CPul step Hilbert propagators
        //                        are copied from CPul1
        // Note                 : Assumed array Usteps is empty

     
void copyLOps(const PulComposite& CPul1);
 
        // Input        CPul    : A composite pulse (this)
        //              CPul1    : A second PulComposite
        // Output       none    : CPul step Hamiltonians superoperators
        //                        are copied from CPul1
        // Note                 : Assumed arrays Lsteps is currently
        //                        currently empty but Hindex is filled
        // Note                 : Recall that some Lsteps[n] may be empty
        //                        as the actual operator for n will reside
        //                        in Lsteps[Hindex[n]]!

 
virtual void copyGprops(const PulComposite& CPul1);
 
        // Input        CPul     : A composite pulse (this)
        //              CPul1    : A second composite pulse
        // Output       none    : CPul step Liouville propagators
        //                        are copied from CPul1
        // Note                 : Assumed array Gsteps is empty


void copySSs(const PulComposite& CPul1);
 
        // Input        CPul    : A composite pulse (this)
        //              CPul1   : A second PulComposite
        // Output       none    : CPul step steady state operators 
        //                        are copied from CPul1
        // Note                 : Assumed array SigSSs is currently empty
        // Note                 : Recall that some SigSSs[n] will be empty
        //                        as the actual operator for n will reside
        //                        in SigSSs[Hindex[n]]!


// ____________________________________________________________________________
// iv       CLASS COMPOSITE PULSE HAMILTONIAN & PROPAGATOR INDICES
// ____________________________________________________________________________


void SetHIndxs( );

        // Input        CPul    : A pulse cycle (this)
        // Output       none    : CPul Hamiltonian indices are
        //                        are set up
        // Note                 : Assumed array Hindex is currently empty



void SetUIndxs( );

        // Input        CPul    : A pulse cycle (this)
        // Output       none    : CPul propagator indices are
        //                        are set up
        // Note                 : Assumed array Uindex is currently empty
        // Note                 : Assumed Hindex is current


// ____________________________________________________________________________
// v            CLASS COMPOSITE PULSE HAMILTONIAN GENERATORS
// ____________________________________________________________________________
 
/* These functions allow for the filling up an array of Hamiltonians (Hsteps)
   for each step in the waveform.  These are typically generated during
   construction of the composite pulse as they are spin system dependent.
   Equivalent steps in a waveform will share the same Hamiltonians as indexed
   in the integer array Hindex.                                              */


void SetHs(const spin_system& sys);

        // Input                CPul    : A composite pulse (this)
        //                      sys     : A spin system
        // Output               void    : Effective Hamiltonians for all of the
        //                                composite pulse steps are generated
        // Note                         : Will delete any stored propagators
        //                                propagators for the waveform
        // Note                         : Assumes H is Ho and that it needs
        //                                F(X,Y,Z) for offsets and phasing.
        // Note                         : ASSUMES PULES CHANNEL SET


void SetHs(gen_op& H, gen_op& FX, gen_op& FY, gen_op& FZ);

        // Input                CPul    : A composite pulse (this)
        //                      H       : Active isotropic Hamiltonian
        //                      FX      : Active Fx operator 
        //                      FY      : Active Fy operator 
        //                      FZ      : Active Fz operator 
        // Output               void    : Effective Hamiltonians for all of the
        //                                composite pulse steps are generated
        // Note                         : This will delete any stored
        //                                propagators for the waveform
        // Note                         : Uses F(X,Y,Z) for phasing & rf
        // Note                         : ASSUMES PULES CHANNEL SET


void SetLs();

        // Input                CPul    : A composite pulse (this)
        // Output               void    : Effective Hamiltonian superops for
        //                                all composite pulse steps generated
        // Note                         : Will delete any stored superprops
        // Note                         : Assumes H is Ho and that it needs
        //                                AND that Hsteps already present



// ____________________________________________________________________________
// vi     CLASS COMPOSITE PULSE HILBERT SPACE PROPAGATOR GENERATORS
// ____________________________________________________________________________

/* These functions allow for the filling up two arrays of propagators for each
   step in the waveform (Usteps and Usums).  These well NOT be generated during
   construction of the composite pulse although as they are system dependent.
   Rather, they will be generated when requested if the step Hamiltonians are
   present.                                                                  */

virtual void SetUs();

        // Input                CPul    : A composite pulse (this)
        // Output               void  : Propagators for each step
        //                              of the composite pulse are calculated
        // Note                       : This assumes that the Hamiltonians
        //                              have already been calculated!
        // Note                       : U[i] is the propagator that evolves
        //                              the system from the start of the
        //                              composite pulse to the end of step i
        // Note                       : This does not destroy superpropagators
        //                              Only Hamiltonian changes do that.

   
// ____________________________________________________________________________
// vii             CLASS COMPOSITE PULSE STEADY-STATE HANDLING
// ____________________________________________________________________________
 
 
void SetSSs();
 
        // Input                CPul    : A composite pulse (this)
        // Output               void    : Steady-state operators for each
        //                                step of the composite pulse are
        //                                calculated

   
// ____________________________________________________________________________
// viii    CLASS COMPOSITE PULSE LIOUVILLE SPACE PROPAGATOR GENERATORS
// ____________________________________________________________________________
 
/* These functions allow for access of individual waveform step Liouville space
   propagators.  These will be propagators for the evolution from pulse start
   to step end.                                                              */

 
virtual void SetGs( );
 
        // Input                CPul    : A composite pulse (this)
        // Output               void    : Superpropagators for each step
        //                                of the composite pulse are calculated
        // Note                         : This assumes that the Hamiltonians
        //                                have already been calculated!
        // Note                         : Assumes the relaxation superoperator
        //                                has already been set
        // Note                         : G[i] is superpropagator that evolves
        //                                the system from the start of the
        //                                composite pulse to the end of step i
        // Note                         : Doesn't destroy any HS propagators


// ____________________________________________________________________________
// ix               CLASS COMPOSITE PULSE BASIS HANDLING
// ____________________________________________________________________________
 
/* These functions are only active on propagators.  They will have no effect
   if there have been no propagators generated.                              */
 
void SetUBasis(gen_op& Op);
 
        // Input                CPul    : A composite pulse (this)
        // Output               void    : Hilbert space propagators for all
        //                                steps of the composite pulse are put
        //                                into the current working basis of Op 
        // Note                         : This assumes that the Hamiltonians
        //                                have already been calculated!
 
void SetGBasis(super_op& LOp);
 
        // Input                CPul    : A composite pulse (this)
        // Output               void    : Liouville space propagators for all
        //                                steps of the composite pulse are put
        //                                into the current working basis of LOp 
        // Note                         : This assumes that the Hamiltonians
        //                                have already been calculated!


// ____________________________________________________________________________
// x                  CLASS COMPOSITE CHECKING FUNCTIONS
// ____________________________________________________________________________


int CheckH(int eflag, int eflag2=0) const;

        // Input                CPul    : A composite pulse (this)
        //                      eflag   : Error index
        //                      eflag2  : Second error index
        // Output               T/F     : Check to insure the composite pulse
        //                                step Hamiltonians are present


int CheckCH(const spin_sys& sys, const std::string& ch, int ef) const;

        // Input                CPul    : A composite pulse (this)
        //                      sys     : A spin system
        //                      ch      : An isotope channel
        //                      ef      : Error index
        // Output               T/F     : Check to insure channel O.K.


int CheckStep(int stp, int eflag, int eflag2=0) const;

        // Input                CPul    : A composite pulse (this)
        //                      stp     : Composite pulse step
        //                      eflag   : Error index
        //                      eflag2  : Second error index
        // Output               T/F     : Check to insure the composite pulse
        //                                step exists


void SetNULL();

        // Input                CPul    : A composite pulse (this)
        // Output               void    : Make everything in CPul NULL
        // Note                         : Call only during construction
        //                                (or after full destruction!)


LSprop LSIprop() const;
 
        // Input                CPul    : A composite pulse (this)
        // Output               GI      : An identity superpropagator
 

 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 public:

// ____________________________________________________________________________
// A                 CLASS COMPOSITE PULSE CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///Center Composite Pulse Algebraic


MSVCDLC PulComposite();

	// Input	none		:
	// Output	CPul		: A NULL composite pulse (this)
	///F_list 	PulComposite	- Constructor

 
        // Input        pulwf   : Pulse waveform
        //              sys     : A spin system
        //              isoch   : An isotope channel 
        // Output       CPul    : A new composite pulse (this)

 
MSVCDLC PulComposite(const PulWaveform& pulwf,
                                  const spin_system& sys, const std::string& isoch);
 
        // Input        pulwf   : Pulse waveform
        //              sys     : A spin system
        //              isoch   : An isotope channel 
        // Output       CPul    : A new composite pulse (this)

 
MSVCDLC PulComposite(const PulWaveform& pulwf,
           gen_op& H, gen_op& FX, gen_op& FY, gen_op& FZ, const std::string& isoch);
 
        // Input        pulwf   : Pulse waveform
        //		H       : Active isotropic Hamiltonian
        //		FX      : Active Fx operator 
        //		FY      : Active Fy operator 
        //		FZ      : Active Fz operator 
        //              isoch   : An isotope channel 
        // Output       CPul    : A new composite pulse (this)
	// Note			: Assumes H,FX,FY,FZ set for isoch



MSVCDLC PulComposite(const PulWaveform& pulwf,
             const spin_system& sys, const super_op& LOp, const std::string& isoch);

        // Input        pulwf   : Pulse waveform
        //              sys     : A spin system
        //              LOp     : Relaxation/exchange superoperator
        //              isoch   : An isotope channel
        // Output       CPul    : A new composite pulse (this)


MSVCDLC PulComposite(const PulComposite& PT1);

	// Input	CPul1	: Composite pulse
	// None		CPul	: Composite pulse (this) constructed from CPul1


// ------------------------------ Destruction ---------------------------------


MSVCDLC virtual ~PulComposite();

	// Input	CPul	: A composite pulse (this)
	// Output	none	: CPul is destructed


// ------------------------------- Assignment ---------------------------------
 

MSVCDLL PulComposite & operator = (const PulComposite& CPul1);

        // Input        CPul1	: Composite pulse
        // None         CPul	: Composite pulse (this) copied from CPul1
        ///F_list =           - Assignment


// ____________________________________________________________________________
// B        CLASS COMPOSITE PULSE HAMILTONIAN/LIOUVILLIAN FUNCTIONS
// ____________________________________________________________________________

/* These functions allow for access of individual waveform step Hamiltonians
   and Liouvillians.                                                         */

MSVCDLL gen_op GetH(int i) const;

        // Input                PT    : A composite pulse (this)
        //                      i     : Step in composite pulse
        // Output               H     : The effective Hamiltonian
        //                              for composite pulse step i


MSVCDLL super_op L0(int i);
MSVCDLL super_op GetL0(int i);

        // Input                CPul    : A composite pulse (this)
        //                      i       : Step in composite pulse waveform
        // Output               L       : The effective Hamiltonian superop
        //                                for composite pulse step i

MSVCDLL super_op Leff(int i);
MSVCDLL super_op GetLeff(int i);

       // Input                CPul    : A composite pulse (this)
        //                      i       : Step in composite pulse waveform
        // Output               L       : The effective Liouvillian superop
        //                                for composite pulse step i
 
/*                      L = -i*[Heff, ] + R    (rad/sec)                     */   


MSVCDLL densop SigSS(int i);

        // Input                CPul    : A composite pulse (this)
        //                      i       : Step in composite pulse waveform
        // Output               SigmaSS : The steady-state density operator
        //                                for step i
 


// ____________________________________________________________________________
// C          CLASS COMPOSITE PULSE HILBERT SPACE PROPAGATOR FUNCTIONS
// ____________________________________________________________________________

/* These functions allow for access of individual waveform step Hilbert space
   propagators.  These can either be for a single step evolution or for the
   evolution from pulse start to step end.                                   */


MSVCDLL virtual HSprop GetU(int i=-1);

	// Input                CPul	: A composite pulse (this)
        //                      i       : Step in composite pulse
        // Output               U       : The propagator for composite pulse
        //                                that evolves from start of step i
        //                                             to the end of step i
        // Note                         : An exception is made for default (-1)
        //                                which returns U for full waveform
        // Note                         : If requested and not existing, this
        //                                triggers generation of ALL propators
	// Note				: Indices here span [0, Wfsteps-1]!


MSVCDLL virtual HSprop GetU(int i, double td);
 
        // Input                CPul    : A composite pulse (this)
	//			i	: Pulse step index
        //                      td      : Evolution time
        // Output               U       : The propagator that evolves a system
        //                                under composite pulse step i for
	//				  time td.
	// Note				: Indices here span [0, Wfsteps-1]!
 

MSVCDLL virtual HSprop GetU(double td);
 
        // Input                CPul    : A composite pulse (this)
        //                      td      : Evolution time
        // Output               U       : The propagator that evolves a system
        //                                for time td under the composite pulse
 

 
MSVCDLL virtual HSprop GetUsum(int i=-1);
 
        // Input                CPul     : A composite pulse (this)
        //                      i       : Step in composite pulse
        //                                Default is last propagator -
        //                                i.e. for the entire composite pulse
        // Output               U       : The propagator for composite pulse
        //                                that evolves from the start of
        //                                waveform to the end of step i
        // Note                         : If requested and not existing,
        //                                this will trigger ALL propators
        //                                to be generated if possible

 
MSVCDLL virtual HSprop GetUmult(int N);

        // Input                CPul    : A composite pulse (this)
        //                      N       : Number of composite pulses
        // Output               U       : The propagator for N composite
        //                                pulses applied in succession

 
// ____________________________________________________________________________
// D        CLASS COMPOSITE PULSE SUPEROPERATOR PROPAGATOR FUNCTIONS
// ____________________________________________________________________________


/* These functions allow for access of individual waveform step Liouville space
   propagators.  These will be for the evolution pulse start to step end.    */

MSVCDLL virtual LSprop GetG(int i=-1);

        // Input                CPul    : A composite pulse (this)
        //                      i       : Step in composite pulse
        //                                Default is last superpropagator -
        //                                i.e. for the entire composite pulse
        // Output               G       : Superpropagator for composite pulse
        //                                that evolves from the start of
        //                                the waveform to the end of step i
        // Note                         : Assumes relaxation operator set


MSVCDLL LSprop GetG(int i, double td);

        // Input                CPul    : A composite pulse (this)
        //                      i       : Pulse step index
        //                      td      : Evolution time
        // Output               G       : The propagator that evolves a system
        //                                under composite pulse step i for
        //                      

MSVCDLL LSprop GetG(double td);

        // Input                CPul    : A composite pulse (this)
        //                      td      : Evolution time
        // Output               G       : Superpropagator that evolves a system
        //                                for time td under the composite pulse
 
// ____________________________________________________________________________
// E                CLASS COMPOSITE PULSE ACCESS FUNCTIONS
// ____________________________________________________________________________
 
 
MSVCDLL virtual LSprop GetGsum(int i=-1);
  
        // Input                CPul    : A composite pulse (this)
        //                      i       : Step in composite pulse 
        //                                Default is last propagator -
        //                                i.e. for the entire composite pulse 
        // Output               G       : Superpropagator for composite pulse 
        //                                that evolves from the start of 
        //                                waveform to the end of step i 
        // Note                         : If requested and not existing, 
        //                                this will trigger ALL propators
        //                                to be generated if possible 
 
 
MSVCDLL LSprop GetGmult(int N);
 
        // Input                CPul    : A composite pulse (this)
        //                      N       : Number of composite pulses
        // Output               G       : Superpropagator for N composite
        //                                pulses applied in succession

// --------------------- Functions Over Full Composite ------------------------

//int        steps()   const;
//std::string     name()    const;
//double     length()  const;
MSVCDLL std::string     channel() const;
//row_vector values()  const;

	// Input	CPul	: A composite pulse (this)
        // Output       steps   : CPul steps
        //              name    : CPul name
        //              length  : CPul length (sec)
        // 		channel	: CPul channel
	//		values  : Array of { gamB1, phi } values

// ----------------- Functions For Specific Composite Step --------------------

//double  strength(int i) const;
//double  phase(int i)    const;
//double  length(int i)   const;
//complex value(int i)    const;
 
        // Input        CPul     : A composite pulse (this)
        // Output       strength: Step rf-field strength (Hz)
        //              phase   : Step phase value (degrees)
        //              value   : Step { gamB1, phi } values
        //              length  : Step length (sec)
 
// --------------------- Other Pulse Composite Access -------------------------
 
MSVCDLL gen_op       FZ()        const;
MSVCDLL super_op     ROp()       const;
MSVCDLL densop       SigEq()     const;
MSVCDLL double       Precision() const;
 
        // Input        CPul    : A composite pulse (this)
        // Output       FZ      : Active Fz operator
        //              waveform: Pulse waveform   
        // Output       R       : Active Relaxation/exchange superoperator
        //              sigmaeq : Equilibrium density operator 
	//		cutzero : Precision with which time is known
 
//double steps(double td) const;			INHERITED
//int fullsteps(double td) const;			INHERITED
//double sumlength(int i) const;			INHERITED
 
  

// ____________________________________________________________________________
// F                CLASS COMPOSITE PULSE AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 

//int gamB1const() const;
//int phaseconst() const;
//int timeconst() const;

        // Input        CPul    : A composite pulse (this)
        // Output       gamB1   : Returns true if the RF field strength
        //                        is constant through out the waveform
        // Output       phase   : Returns true if the RF field phase
        //                        is constant through out the waveform
        // Output       time    : Returns true if the step time
        //                        is constant through out the waveform


MSVCDLL virtual void scalegB1(double sf); 
  
        // Input        CPul    : A composite pulse (this)
        //              sf      : A scaling factor
        // Output       void    : All step field strengths in PWF
        //                        are multiplied by sf.  The exception are 
        //                        steps of zero length (ideal pulses) 


MSVCDLL void setRelax(const spin_system& sys, const super_op& LOp);
  
        // Input        CPul    : A composite pulse (this)
        //              sys     : A spin system
        //              LOp     : Relaxation/exchange superoperator
        // Output       void    : Relaxation specifications are put int
        //                        CPul


// ____________________________________________________________________________
// G               CLASS COMPOSITE PULSE EVOLUTION FUNCTIONS
// ____________________________________________________________________________
 

// ----------------------- Acquisiton Tracking Functions ----------------------
 
MSVCDLL void FIDheader(int typ, int rlx=0) const;
  
        // Input        CPul    : A composite pulse (this)
        //              typ     : FID header type 
	//		rlx     : Relaxation flag
        // Output       void    : Outputs header for FID tracking
 

MSVCDLL void FIDpoint(int typ, int pt, int iWFs, int iSTs) const;

        // Input        CPul    : A composite pulse (this)
        //              typ     : FID header type
        //              pt      : Point index
        //              iWFs    : Number of waveforms
        //              iSTs    : Number of steps
        // Output       void    : Outputs header for FID tracking


MSVCDLL void FIDvalue(int typ, double td, const complex& z) const;

        // Input        CPul    : A composite pulse (this)
        //              typ     : FID header type
        //              td      : Delay length
        //              z       : Data point
        // Output       void    : Outputs point info for FID tracking


MSVCDLL void FIDtell(double SW) const;

        // Input        CPul    : A composite pulse (this)
        //              SW      : Desired spectral width
 

MSVCDLL virtual double FIDsync(double& SW, int warn=0) const;

        // Input        CPul    : A composite pulse (this)
        //              SW      : Desired spectral width
	//		warn	: Flag for warning output
        // Output       SWsync  : Spectral width synchronized with the
        //                        composite pulse length (at best) or
        //                        with the pulse step length (2nd best)
        // Note                 : If no synchronization possible the input
        //                        value is just returned.

 
MSVCDLL virtual int FIDtest(double td, int& nWFs, int& nSTPs, double& tr) const;
 
        // Input        CPul    : A composite pulse (this)
        //              td      : Desired dwell time (sec)
        //              nWFs    : Number of full waveforms
        //              nSTPs   : Number of full step
        //              tr      : Remainder time (sec)
        // Output       Fsync   : Synchronization flag
        //                              +: Synchronized to Waveform
        //                              0: Not synchronized
        //                              -: Synchronized to Step
        // Note                 : Values of nWFs, nSTPs, and tr are altered


// ----------------------- Without Relaxaton & Exchange -----------------------


MSVCDLL row_vector FIDsynchWF(int npts, int nWFs,
                                       gen_op &D, gen_op& sigmap, int track=0);

        // Input        CPul    : A composite pulse (this)
        //              npts    : Number of FID points
        //              nWFs    : Waveforms between FID points
        //              D       : Detection operator
        //              sigmap  : Prepared density operator
        //              track   : Print output to track calculation
        // Output       data    : A row vector containing npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
        //                        composite pulse CPul.
 

MSVCDLL row_vector FIDsynchST(int npts, int nSTs,
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

 
MSVCDLL row_vector FIDsynchFR(int npts, int nFRs,
                                       gen_op &D, gen_op& sigmap, int track=0);

        // Input        CPul    : A composite pulse (this)
        //              npts    : Number of FID points
        //              nFRs    : Number of WF fractions between points
        //              D       : Detection operator
        //              sigmap  : Prepared density operator
        //              track   : Flag for tracking the computation
        // Output       data    : A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
        //                        composite pulse CPul.


MSVCDLL virtual row_vector FID(int N, double td, gen_op &D, gen_op& sp, int track=0);
 
        // Input        CPul    : A composite pulse (this)
        //              N	: Number of FID points
        //              td	: Dwell time between FID points
        //              D	: Detection operator
        //              sp	: Prepared density operator
	//		track   : Print output to track calculation
        // Output       data 	: A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
        //                        composite pulse CPul.


// ------------------------- With Relaxaton & Exchange ------------------------

MSVCDLL virtual row_vector FIDRsynchWF(int npts, int nWFs,
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
 
 
MSVCDLL virtual row_vector FIDRsynchST(int npts, int nSTs,
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

 
MSVCDLL virtual row_vector FIDRsynchFR(int npts, int nFRs,
                                       gen_op &D, gen_op& sigmap, int track=0);
 
        // Input        CPul    : A composite pulse (this)        
        //              npts    : Number of FID points 
        //              nFRs    : Number of WF fractions between points
        //              D       : Detection operator 
        //              sigmap  : Prepared density operator 
        //              track   : Flag for tracking the computation 
        // Output       data    : A row vector contiaing an npts point FID 
        //                        that was generated by detection with D 
        //                        of the evolution of sigmap under the 
        //                        composite pulse CPul. 

 
MSVCDLL virtual row_vector FIDR(int N, double td, gen_op &D, gen_op& sp, int track=0);

        // Input        CPul    : A composite pulse (this)
        //              N	: Number of FID points
        //              td	: Dwell time between FID points
        //              D	: Detection operator
        //              sp	: Prepared density operator
	//		track   : Print output to track calculation
        // Output       data	: A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
        //                        composite pulse CPul.
        // Note			: Assumes Relaxation Active!


// ____________________________________________________________________________
// H                 CLASS COMPOSITE PULSE PLOTTING FUNCTIONS
// ____________________________________________________________________________


//void GP(int type, int split, int ends) const;         INHERITED

// ____________________________________________________________________________
// I                 CLASS COMPOSITE PULSE I/O FUNCTIONS
// ____________________________________________________________________________
 
 
MSVCDLL virtual std::ostream& printEvolve(std::ostream &ostr, double td) const;
 
        // Input                PWF     : Pulse Waveform 
        //                      td      : Evolution time 
        //                      ostr    : Output stream 
        // Output               none    : Pulse Waveform evolution info 
        //                                is sent to the output stream 
 
 
MSVCDLL virtual std::ostream& printFID(std::ostream &ostr, double td, int n) const;
 
        // Input        PT    : A pulse train (this) 
        //              ostr  : An output stream 
        //              td    : Dwell time between FID points 
        //              n     : Number of FID points
        // Output             : Information regarding the FID generation 
        //                      is set to the output stream ostr
 

MSVCDLL virtual std::ostream& printInfo(std::ostream &ostr) const;

        // Input                CPul     : Composite Pulse
        //                      ostr    : Output stream
        //                      full    : Flag for output amount
        // Output               none    : CPul storage info is sent
        //                                to the output stream

 
MSVCDLL virtual std::ostream& print(std::ostream &ostr, int full=0) const;

        // Input                CPul    : Composite Pulse
        //                      ostr	: Output stream
        //                      full	: Flag for output amount
        // Output               none	: Composite Pulse info is sent
        //                                to the output stream


MSVCDLL friend std::ostream &operator << (std::ostream &ostr, const PulComposite &CPul);

	// Input		CPul	: Composite Pulse
        //                      ostr	: Output stream
        // Output               none	: Composite composite pulse is sent
	//				  to the output stream
};

#endif						// PulComposite.h
