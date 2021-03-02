/* PulTrainSCyc.cc **********************************************-*-c++-*-
**									**
**      	                G A M M A				**
**									**
**	Class Pulse Train Supercycle		  Implementation	**
**                                                                      **
**      Copyright (c) 1998                                              **
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**                                                                      **
** This class sets up a supercycle for GAMMA's pulse train class.  A	**
** pulse cycle consists of a defined waveform which is repeated (cycled)**
** over a specified number of steps with a defined phase.  In turn, a	**
** pulse supercycle is an additonal loop which repeats the pulse cycle	**
** over a specified number of steps with a defined phase.  Information	**
** regarding the supercycle steps and phases is engulfed in class 	**
** PulSupCycle, the class which serves as a base to this class. However	**
** PulSupCycle has only minimal functionality and knows nothing about 	**
** the applied pulse cycle nor the cycle step waveform.  This class	**
** also knows the pulse supercycle, but in addition handles the prop- 	**
** agators for spin system evolution as required in NMR simulations. 	**
** Note that it still has little to no knowledge of the applied cycle,	**
** the cycle step waveform, and is merely meant as a auxiliary class to	**
** GAMMA's Pulse Train class (PulTrain) which will know specifics of	**
** the waveform, cycle, and supercycle including applicable Hamitonians	**
** & propagators.  By and large it is class PulTrain that is used in 	**
** GAMMA programs, not this class.					**
**                                                                      **
*************************************************************************/

#ifndef _PulTrainSCyc_cc_		// Is this file already included?
#define _PulTrainSCyc_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Pulses/PulTrainSCyc.h>	// Include the header
#include <Pulses/PulSupCycle.h>		// Know about pulse supercycles
#include <Pulses/PulCycle.h>		// Know about pulse cycles
#include <HSLib/SpinSystem.h>		// Must know about spin systems 
#include <HSLib/HSprop.h>		// Know Hilbert space props
#include <LSLib/LSprop.h>		// Know Liouville space props
#include <string>			// Must know about strings


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                 CLASS PULSE TRAIN SUPERCYCLE ERROR HANDLING
// ____________________________________________________________________________

 
void PulTrainSCyc::PTSCerror(int eidx, int noret) const
 
        // Input                PTSC	: Pulse Train (this)
        //                      eidx    : Error flag
        //                      noret   : Return flag 
        // Output               none    : Error Message Output
 

  {
  std::cout << "\nClass Pulse Train Supercycle: ";
  switch(eidx)
    {
    case 0:								// (0)
      std::cout << "Program Aborting....";
      break;
    case 1:								// (1)
      std::cout << "Error During Construction";
      break;
    case 6:								// (6)
      std::cout << "Build Step Propagators Before Requesting Propagators";
      break;
    case 30:								// (30)
      std::cout << "Step Hamiltonian Access, Hamiltonian Does Not Exist";
      break;
    case 31:								// (31)
      std::cout << "Build Step Hamiltonians Before Their Access";
      break;
    case 32:								// (32)
      std::cout << "Step Propagator Access, Propagator Does Not Exist";
      break;
    case 52:								// (52)
      std::cout << "Step Propagator Access, Superop. Propagator Does Not Exist";
      break;
    default:
      std::cout << "Unknown Error (Number " << eidx << ")";
    }
  if(!noret) std::cout << ".\n";
  else       std::cout << ".";
  }
 
 
void volatile PulTrainSCyc::PTSCfatality(int eidx) const
 
        // Input                PTSC	: Pulse Train (this)
        //                      eidx    : Error flag
        // Output               none	: Stops execution & error Message
 
  {
  PTSCerror(eidx, 1);
  if(eidx) PTSCerror(0);
  GAMMAfatal();					// Clean exit from program
  }


// ____________________________________________________________________________
// ii            CLASS PULSE TRAIN SUPERCYCLE DESTRUCTION FACILITATORS 
// ____________________________________________________________________________


void PulTrainSCyc::deleteUprops()

	// Input	PTSC	: A pulse train supercycle (this)
	// Output	none	: PTSC step Hilbert propagators
	//			  are deleted if they exist

  {
  if(Usteps) delete [] Usteps;		// Delete current propagators
  Usteps = NULL;			// Set them to NULL
  if(Usums) delete [] Usums;		// Delete current sum propagators
  Usums = NULL;				// Set them to NULL
  Ucount = 0;
  }


void PulTrainSCyc::deleteGprops()

	// Input	PTSC	: A pulse train supercycle (this)
	// Output	none	: PTSC step superoperator propagators
	//			  are deleted if they exist

  {
  if(Gsteps) delete [] Gsteps;		// Delete current superpropagators
  Gsteps = NULL;			// Set them to NULL
  }

// ____________________________________________________________________________
// iii              CLASS PULSE TRAIN SUPERCYCLE COPY FACILITATORS 
// ____________________________________________________________________________


void PulTrainSCyc::copyUprops(const PulTrainSCyc& PTSC1)

	// Input	PTSC	: A pulse train supercycle (this)
	//		PTSC1    : A second PulTrainSCyc
	// Output	none	: PTSC step Hilbert propagators
	//			  are copied from PTSC1
	// Note			: Assumed array Usteps is empty

  {
  Usteps = NULL; 
  Usums = NULL;
  Ucount = 0;
  if(PTSC1.Usteps)			// If the propagators exist
    {
    Usteps = new HSprop[SCycnosteps];	//	Allocate space for them
    for(int i=0; i<SCycnosteps; i++)	//	Copy all of the propagators
      Usteps[i] = PTSC1.Usteps[i];
    }
  if(PTSC1.Usums)			// If sum propagators exist
    {
    Usums = new HSprop[SCycnosteps];	//	Allocate space for them
    for(int i=0; i<SCycnosteps; i++)	//	Copy all of the propagators
      Usums[i] = PTSC1.Usums[i];
    }
  }


void PulTrainSCyc::copyGprops(const PulTrainSCyc& PTSC1)

	// Input	PTSC	: A pulse train supercycle (this)
	//		PTSC1    : A second PulTrainSCyc
	// Output	none	: PTSC step Liouville propagators
	//			  are copied from PTSC1
	// Note			: Assumed array Gsteps is empty

  {
  Gsteps = NULL; 
  if(PTSC1.Gsteps)			// If the super propagators exist
    {
    Gsteps = new LSprop[SCycnosteps];	//	Allocate space for them
    for(int i=0; i<SCycnosteps; i++)	//	Copy all of the superprops
      Gsteps[i] = PTSC1.Gsteps[i];
    }
  }

    
// ____________________________________________________________________________
// iv     CLASS PULSE TRAIN SUPERCYCLE HILBERT SPACE PROPAGATOR GENERATORS
// ____________________________________________________________________________
 
/* These functions allow for the filling up two arrays of propagators for each
   step in the waveform (Usteps and Usums).  These will usually be generated
   during construction of the pulse train supercycle, they are system dependent.  */


void PulTrainSCyc::SetUs(PulCycle& PTC)

        // Input                PTSC	: A pulse train supercycle (this)
        //                      PTC     : A pulse train cycle
        // Output               void	: Propagators for each step of the
        //                                pulse train supercycle are calculated
        // Note				: This assumes that the propagators
        //                                in PTC have already been calculated!

  {
  if(!SCycnosteps) return;		// Do nothing if no cycle steps defined
  if(!Usteps)				// If no cycle step propagators
      Usteps = new HSprop[SCycnosteps];	// allocate space for them
  if(!Usums)				// If no cycle step summed propagators
      Usums = new HSprop[SCycnosteps];	// allocate space for them
  HSprop Ubase = PTC.GetU();		// Begin with pulse waveform propagator
  gen_op RZ;				// Rotation operator for phasing
  double phi;                           // Cycle step phase
  for(int i=0; i<SCycnosteps; i++)
    {
    phi = DEG2RAD*SCycsteps.getRe(i);	// Cycle step i phase (rad)
    RZ = exp((-complexi*phi)*PTC.FZ());// Rotation operator for phi
    Usteps[i] = Ubase.sim_trans(RZ);	// Phase the base propagator
    if(i) Usteps[i] *= Usteps[i-1];	// Add in previous cycle steps
    }
  }
 
// sosix


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A            CLASS PULSE TRAIN SUPERCYCLE CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________


PulTrainSCyc::PulTrainSCyc() : PulSupCycle()

	// Input	none	:
	// Output	PTSC	: A NULL pulse train supercycle (this)

  {
  tp        = 0;			// No pulse train supercycle length
  Usteps    = NULL;			// No pulse train supercycle propagators
  Usums     = NULL;			// No pulse train supercycle sum propagators
  Gsteps    = NULL;			// No pulse train supercycle superprops
  Ucount    = 0;			// No stored propagators
  }


PulTrainSCyc::PulTrainSCyc(PulCycle& P, const row_vector& S, std::string N)
             :PulSupCycle(S, N)

	// Input	P	: Composite pulse
	// 		S	: Pulse train cycle steps (phase)
	//		N	: Pulse train cycle name
	// Output	PTSC	: A new pulse train supercycle (this)

  {
  tp         = P.length()*SCycnosteps;	// No pulse train supercycle length
  Ucount     = 0;			// No stored propagators
  Usteps     = NULL;			// No pulse train supercycle propagators
  Usums      = NULL;			// No pulse train supercycle sum propagators
  Gsteps     = NULL;			// No pulse train supercycle superpropagators
  SetUs(P);
  }


PulTrainSCyc::PulTrainSCyc(PulCycle& P, const PulSupCycle& PCY)
             :PulSupCycle(PCY)

	// Input	P	: Composite pulse
	// 		PCY	: Pulse cycle 
	// Output	PTSC	: A new pulse train supercycle (this)

  {
  tp         = P.length()*SCycnosteps;	// No pulse train supercycle length
  Ucount     = 0;			// No stored propagators
  Usteps     = NULL;			// No pulse train supercycle propagators
  Usums      = NULL;			// No pulse train supercycle sum propagators
  Gsteps     = NULL;			// No pulse train supercycle superpropagators
  SetUs(P);
  }



PulTrainSCyc::PulTrainSCyc(const PulTrainSCyc& PTSC1) : PulSupCycle(PTSC1)

	// Input	PTSC1	: Pulse Train Supercycle
	// None		PTSC	: Pulse Train Supercycle (this), identical
	//			  copy of PTSC1

  {
  tp = PTSC1.tp;			// Copy pulse train supercycle length
  copyUprops(PTSC1);			// Copy any Hilbert Propagators
  copyGprops(PTSC1);			// Copy any Liouville Propagators
  }


// ------------------------------ Destruction ---------------------------------

PulTrainSCyc::~PulTrainSCyc()

	// Input	PTSC	: A pulse train supercycle (this)
	// Output	none	: PTSC is destructed

  { 
  if(Usteps) delete [] Usteps;		// Delete any current propagators
  if(Usums)  delete [] Usums;		// Delete any current sum propagators
  if(Gsteps) delete [] Gsteps;		// Delete any current superpropagators
  }

// ------------------------------- Assignment ---------------------------------


PulTrainSCyc& PulTrainSCyc::operator = (const PulTrainSCyc& PTSC1)

	// Input	PTSC1	: Pulse Train
	// None		PTSC	: Pulse Train (this) copied from PTSC1

{
  PulSupCycle::operator=(PTSC1);		// Copy pulse train definition
  deleteUprops();			// Delete current propagators
  deleteGprops();			// Delete current superpropagators
  tp = PTSC1.tp;				// Copy pulse train supercycle length
  copyUprops(PTSC1);			// Copy any Hilbert Propagators
  copyGprops(PTSC1);			// Copy any Liouville Propagators

  return (*this);
}

// ____________________________________________________________________________
// B                CLASS PULSE TRAIN SUPERCYCLE HAMILTONIAN FUNCTIONS
// ____________________________________________________________________________

/* These functions allow accessing average Hamiltonians over the pulse train
   cycle steps.  They are generated from an input pulse waveform.  It is the
   waveform which actually contains individual step Hamiltonians.           */

// sosi
//gen_op PulTrainSCyc::GetH(PulCycle& PW, int i) const

        // Input                PTSC    : A pulse train supercycle (this)
        //                      i     : Step in pulse train supercycle
        // Output               H     : The effective Hamiltonian
	//				for pulse train supercycle step i

//  {
//  if(!Hsteps)				// If Hamiltonians dont exist
//    {					// then this is a fata error
//    PTSCerror(30);
//    PTSCfatality(31);
//    }
//  return Hsteps[(Hindex[i])];
//  }


// ____________________________________________________________________________
// C                CLASS PULSE TRAIN SUPERCYCLE PROPAGATOR FUNCTIONS
// ____________________________________________________________________________

/* These functions allow access and production of propagators over the pulse
   train cycle steps.  They are generated from an input pulse waveform.  It
   is the waveform which actually contains individual step propagators.      */ 
   
HSprop PulTrainSCyc::GetU(int i) const
 
        // Input                PTSC	: A pulse train supercycle (this)
        //                      i       : Step in pulse train supercycle
        // Output               U       : The propagator for pulse train supercycle
        //                                that evolves from start of step i
        //                                             to the end of step i
	// Note				: Steps are indexed [0, nsteps-1)
        // Note                         : An exception is made for default
        //                                then the return will be the
        //                                propagator for the entire pulse
        //                                train cycle

  {
  if(i<0) return GetUsum(i);
  if(!Usteps)                           // If propagators dont exist
    {                                   // then end of story
    PTSCerror(32, 1);
    PTSCerror(30, 1);
    PTSCfatality(31);
    }
  return Usteps[i];                     // Just return the requested U
  }  
 
 
HSprop PulTrainSCyc::GetUsum(int i) const

        // Input                PTSC	: A pulse train supercycle (this)
        //                      i	: Steps in pulse train supercycle
	//				  where each step is a cycle
        // Output               U	: Propagator that evolves through
	//				  i steps of the pulse train supercycle
	// Note				: If requested and not existing,
	//				  this will trigger ALL propators 
	//				  to be generated if possible
        // Note				: Default propagator is for full
	//				  supercycle (last in Usums array)

  {
  if(!Usums && i)			// If propagators dont exist
    {					// end of story, we'll quit
    PTSCerror(32, 1);
    PTSCerror(30, 1);
    PTSCfatality(31);
    }
  else if(!i) 				// Identity if no steps
    if(!Usums) return HSprop(1);
    else       return HSprop(Usteps[0].dim());
  else if(i<0) i=SCycnosteps;		// Set for total cycle prop
  return Usums[i-1];			// Just return the requested U
  }

   
HSprop PulTrainSCyc::GetUmult(int N) const
 
        // Input                PTSC	: A pulse train supercycle (this)
        //                      N       : Number of pulse train supercycles
        // Output               U       : The propagator for N pulse train
	//				  supercycles applied in succession

  {
  if(!Usteps)					// If propagators dont exist
    {						// we can't do anything
    PTSCerror(32, 1);
    PTSCerror(30, 1);
    PTSCfatality(31);
    }
  if(N<=0) return HSprop(Usteps[0].dim());	// Identity if no steps
  HSprop USC = Usums[SCycnosteps-1];		// Prop for 1 supercycle
  HSprop UNSC = USC;				// Set for 1 supercycle
  for(int i=1; i<N; i++) UNSC *= USC;		// Add in next suupercycles
  return UNSC;
  }


// ____________________________________________________________________________
// D       CLASS PULSE TRAIN SUPERCYCLE SUPEROPERATOR PROPAGATOR FUNCTIONS
// ____________________________________________________________________________

/* These functions allow access and production of superoperator propagators
   over the pulse train supercycle steps.  They are generated from an input pulse
   waveform.  Its the waveform which contains individual step propagators.   */ 


void PulTrainSCyc::SetGs(PulCycle& PW)

        // Input                PTSC	: A pulse train supercycle (this)
	//			PW	: A pulse waveform
        // Output               void	: Superpropagators for each step of
	//				  the pulse train supercycle are calculated
        // Note				: This demands that the Liouville
        //                                propagators have already been made
	//				  in the waveform.

  {
  if(!SCycnosteps) return;		// Do nothing if no cycle steps defined
  if(!Gsteps)				// If no cycle step propagators
    Gsteps = new LSprop[SCycnosteps];	// allocate space for them

// sosix
LSprop Gbase;
//  LSprop Gbase = PW.GetG();		// Begin with waveform propagator
  gen_op RZ;                            // Rotation operator
  super_op GZ;                          // Rotation superoperator
  double phi;                           // Cycle step phase
  for(int i=0; i<SCycnosteps; i++)
    {
    phi = DEG2RAD*SCycsteps.getRe(i);	// Cycle step i phase (rad)
    RZ = exp((-complexi*phi)*PW.FZ());	// Rotation operator for phi
    GZ = U_transform(RZ);               // Rotation superopertor for phi
    Gsteps[i] = GZ*Gbase;		// Phase the base propagator
    if(i) Gsteps[i] *= Gsteps[i-1];	// Add in previous cycle steps
    }
  }


LSprop PulTrainSCyc::GetG(int i) const

        // Input                PTSC    : A pulse train supercycle (this)
        //                      i     : Step in pulse train supercycle
        //                              Default is last superpropagator -
        //                              i.e. for the entire pulse train supercycle
        // Output               G     : The superpropagator for pulse train supercycle
        //                              that evolves from the start of
        //                              the waveform to the end of step i
	// Note			      : Assumes relaxation operator set

  {
  if(!Gsteps)				// If propagators dont exist
    {					// check if the Hamiltonians exist
    PTSCerror(52);
    PTSCerror(30);
    PTSCfatality(31);
    }
  if(i<0) i=SCycnosteps-1;		// Set for total cycle prop
  return Gsteps[i];			// Just return the requested G
  }


// ____________________________________________________________________________
// D               CLASS PULSE TRAIN SUPERCYCLE ACCESS FUNCTIONS
// ____________________________________________________________________________

// ---------------------- Functions Over Full Cycle ---------------------------

//int        PulTrainSCyc::steps()   const { return SCycnosteps; }    INHERITED
//std::string     PulTrainSCyc::name()    const { return SCycname; }       INHERITED
//row_vector PulTrainSCyc::values()  const { return SCycsteps; }      INHERITED

int        PulTrainSCyc::steps()   const { return SCycnosteps; }

double     PulTrainSCyc::length()  const { return tp; }

	// Input	PTSC	: A pulse train supercycle (this)
	// Output	steps	: PTSC steps
	// 		name	: PTSC name
	//		length  : PTSC length (sec)
        //              values  : Array of phi values


// ------------------ Functions For Specific Cycle Step ----------------------

//double  PulTrainSCyc::phase(int i)    const { return SCycsteps.getRe(i); }

	// Input	PTSC	: A pulse train supercycle (this)
	// Output	phase	: Step phase value (degrees)

 
double PulTrainSCyc::steps(double td) const

	// Input	PTSC	: A pulse train supercycle (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of supercycles steps needed
        //                        to evolve for time td

  { 
  if(td <= 0)       return 0;
  else if(tp == 0)  return 0;
  else              return td/tp;
  }


int PulTrainSCyc::fullSCYCs(double td) const { return fullsteps(td); }
int PulTrainSCyc::fullsteps(double td) const

	// Input	PTSC	: A pulse train supercycle (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of full supercycles that can
        //                        occur in the time td

  { 
  if(td < 0)       return PulSupCycle::steps();
  else if(tp == 0)  return 0;
  int ns=-1;
  while(td>=0)
    {
    td -= tp;
    ns++;
    }
  return ns;
  }


// ____________________________________________________________________________
// F               CLASS PULSE TRAIN SUPERCYCLE AUXILIARY FUNCTIONS
// ____________________________________________________________________________





// ____________________________________________________________________________
// Z                 CLASS PULSE TRAIN SUPERCYCLE I/O FUNCTIONS
// ____________________________________________________________________________


std::ostream& PulTrainSCyc::printInfo(std::ostream &ostr) const

	// Input		PTSC	: Pulse Train Supercycle
        //                      ostr	: Output stream
	//			full	: Flag for output amount
        // Output               none	: PTSC storage info is sent
	//				  to the output stream

  {
  ostr << "\n\tCycle Propagators:               ";
  if(Usteps || Usums)
    {
    ostr << "Present (";
    if(Usteps)
      {
      ostr << "steps";
      if(Usums) ostr << ", ";
      }
    if(Usums) ostr << "sums";
    ostr << ")";
    }
  else                ostr << "Absent";
  ostr << "\n\tCycle SuperPropagators:          ";
  if(Gsteps) ostr << "Present";
  else       ostr << "Absent";
  return ostr;
  }


std::ostream& PulTrainSCyc::printBase(std::ostream &ostr) const

	// Input		PTSC	: Pulse Train Supercycle
        //                      ostr	: Output stream
        // Output               none	: Pulse Train base info is
	//				  sent to the output stream

  {
  ostr << "\n\tCycle Steps:                     " << SCycnosteps;
  if(tp)
    {
    ostr << "\n\tCycle Length:                    ";
    if(tp > 0.1)         ostr << tp      << " sec";
    else if(tp > 0.0001) ostr << tp*1.e3 << " msec";
    else                 ostr << tp*1.e6 << " nsec";
    ostr << "\n\tCycle Spectral Width:            ";
    double SW = 1.0/tp;
    if(SW < 1000)        ostr << SW << " Hz";
    else if(SW < 100000) ostr << SW*1.e-3 << " KHz";  
    else                 ostr << SW*1.e-6 << " MHz";

    ostr << "\n\tCycle Step Length:               ";
    double cslen = tp/SCycnosteps;
    if(cslen > 0.1)         ostr << cslen      << " sec";
    else if(cslen > 0.0001) ostr << cslen*1.e3 << " msec";
    else                    ostr << cslen*1.e6 << " nsec";
    ostr << "\n\tCycle Step Spectral Width:       ";
    double cssw = 1.0/cslen;
    if(cssw < 1000)        ostr << cssw << " Hz";
    else if(cssw < 100000) ostr << cssw*1.e-3 << " KHz";
    else                   ostr << cssw*1.e-6 << " MHz";
    }
  return ostr;
  }


std::ostream& PulTrainSCyc::print(std::ostream &ostr, int full) const

	// Input		PTSC	: Pulse Train Supercycle
        //                      ostr	: Output stream
	//			full	: Flag for output amount
        // Output               none	: Pulse Train info is sent
	//				  to the output stream

  {
  if(!SCycnosteps)
    {
    ostr << "\n\n\t\t\tEmpty Pulse Train Supercycle\n\n";
    return ostr;
    }
  ostr << "\n\n\t\t\t";
  if(!full) ostr << "\t";
  ostr << "Pulse Train Supercycle " << SCycname << "\n";
  printBase(ostr);
  if(full)
    {
    printInfo(ostr);
    if(full) printSteps(ostr);		// Print cycle steps
    }
  ostr << "\n\n"; 
  return ostr;
  }

            
std::ostream &operator << (std::ostream &ostr, const PulTrainSCyc &PTSC)

	// Input		PTSC	: Pulse Train Supercycle
        //                      ostr	: Output stream
        // Output               none	: Pulse Train info is sent
	//				  to the output stream

  {
  PTSC.print(ostr);
  return ostr;
  }
 
#endif						// PulTrainSCyc.cc
