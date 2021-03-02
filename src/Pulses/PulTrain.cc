/* PulTrain.cc **************************************************-*-c++-*-
**									**
**      	                G A M M A				**
**									**
**	Class Pulse Train                 	   Implementation	**
**                                                                      **
**      Copyright (c) 1995                                              **
**      Dr. Scott A. Smith                                              **
**      Philippe Pelupessy                                              **
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
** This file contains a definition of Pulse Trains in GAMMA.  These     **
** serve to facilitate application of pulse trains in the simulation    **
** of NMR experiments.                                                  **
**                                                                      **
** Each pulse train consists of a defined Waveform, Cycle, and Super-   **
** Cycle.  The Waveform is repeated with phase changes as specified by  **
** the Cycle and the Cycle is repeated with phase changes as specified  **
** by the SuperCycle.  Neither the Cycle or SuperCycle is mandatory.    **
** The Waveform (function and propagators) is handled by Composite      **
** Pulse class (PulComposite) which serves as the base class herein.    **
**                                                                      **
*************************************************************************/
 
#ifndef _PulTrain_cc_			// Is this file already included?
#define _PulTrain_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Pulses/PulTrain.h>		// Include the header
#include <Pulses/PulAuxil.h>		// Include auxiliary functions 
#include <Pulses/PulCycle.h>		// Know about pulse cycles
#include <Pulses/PulSupCycle.h>		// Know about pulse supercycles
#include <Pulses/PulTrainSCyc.h>	// Know pulse train supercycles
#include <HSLib/SpinSystem.h>		// Must know about spin systems 
#include <HSLib/HSprop.h>		// Know about Hilbert space propagators
#include <HSLib/HSauxil.h>		// Know about prop function
#include <string>			// Must know about strings
#include <Basics/StringCut.h>


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS PULSE TRAIN ERROR HANDLING
// ____________________________________________________________________________


 
void PulTrain::PTerror(int eidx, int noret) const
 
        // Input                PT      : Pulse Train (this)
        //                      eidx    : Error flag
        //                      noret   : Return flag 
        // Output               none    : Error Message Output
 

  {
  std::cout << "\nClass Pulse Train: ";
  switch(eidx)
    {
    case 0:                                                            // (0)
      std::cout << "Program Aborting....";
      break;
    case 1:                                                            // (1)
      std::cout << "Error During Construction";
      break;
    case 3:                                                            // (3)
      std::cout << "Construction From Improper GAMMA Parameter";
      break;
    case 5:                                                            // (5)
      std::cout << "Step Propagators Requested Without Step Hamiltonians";
      break;
    case 6:                                                            // (6)
      std::cout << "Build Step Hamiltonians Before Requesting Propagators";
      break;
//    case 10:                                                           // (10)
//      std::cout << "Unknown Cycle Requested - " << _cname;
//      break;
//    case 20:                                                           // (20)
//      std::cout << "Unknown SuperCycle Requested - " << _scname;
//      break;
    case 30:								// (30)
      std::cout << "Step Hamiltonian Access, Hamiltonian Does Not Exist";
      break;
    case 31:								// (31)
      std::cout << "Build Step Hamiltonians Before Their Access";
      break;
    case 32:								// (32)
      std::cout << "Step Propagator Access, Propagator Does Not Exist";
      break;
    case 33:								// (33)
      std::cout << "Cycle Propagator Access, Propagator Does Not Exist";
      break;
    case 34:								// (34)
      std::cout << "Supercycle Propagator Access, Propagator Does Not Exist";
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
 
 
void volatile PulTrain::PTfatality(int eidx)
 
        // Input                eidx: Error flag
        // Output               none : Stops execution & error Message
 
  {
  PTerror(eidx);
  if(eidx) PTerror(0);
  GAMMAfatal();					// Clean exit from program
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A               CLASS PULSE TRAIN CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________


PulTrain::PulTrain() : PulComposite() { Name = ""; Cycles=0; SCycles=0; }

	// Input	none	:
	// Output	PT	: A NULL pulse train (this)
	// Note			: Pulse train waveform automatically 
	//			  initialized to null


PulTrain::PulTrain(const PulComposite& CPul, std::string N) : PulComposite(CPul)

        // Input        CPul	: Composite pulse
        //              N	: Pulse train name
        // Output       PT	: Pulse train (this) constructed
	// Note			: There will be no defined cycles
	//			  or supercycles in the train 

  { Name = N; Cycles=0; SCycles=0; }


PulTrain::PulTrain(const PulComposite& CPul, const PulCycle& Cyc, std::string N)
         : PulComposite(CPul)

        // Input        CPul	: Composite pulse
        //              Cyc	: Pulse cycle
        //              N	: Pulse train name
        // Output       PT	: Pulse train (this) constructed

  {
  Name = N;				// Set pulse train cycle name
  PTCyc  = Cyc;				// Set pulse train cycle
  Cycles = 1;				// Flag cycles as defined
  SCycles = 0;				// Flag no supercycles defined
  }


PulTrain::PulTrain(const PulComposite& CPul, const PulCycle& Cyc,
                                     const PulSupCycle& SCyc, std::string N)
         : PulComposite(CPul)

        // Input        CPul	: Composite pulse
        //              Cyc	: Pulse cycle
        //              N	: Pulse train name
        // Output       PT	: Pulse train (this) constructed

  {
  Name    = N;				// Set pulse train cycle name
  PTCyc   = Cyc;			// Set pulse train cycle
  PTSCyc  = PulTrainSCyc(PTCyc, SCyc);	// Set pulse train supercycle
  Cycles  = 1;				// Flag cycles as defined
  SCycles = 1;				// Flag supercycles defined
  }



PulTrain::PulTrain(const PulTrain& PT1) : PulComposite(PT1)

	// Input	PT1  : Pulse Train
	// None		PT    : Pulse Train (this), an identical
	//			copy of PT1

  {
  Name    = PT1.Name;			// Copy pulse train name
  PTCyc   = PT1.PTCyc;			// Copy pulse train cycle
  PTSCyc  = PT1.PTSCyc;			// Copy pulse train supercycle
  Cycles  = PT1.Cycles;			// Copy pulse train cycle flag 
  SCycles = PT1.SCycles;		// Copy pulse train supercycle flag 
  }


// ------------------------------ Destruction ---------------------------------

PulTrain::~PulTrain() { }

	// Input	PT    : An PulTrain (this)
	// Output	none  : PT is destructed



// ------------------------------- Assignment ---------------------------------


PulTrain& PulTrain::operator = (const PulTrain& PT1)

	// Input	PT1   : Pulse Train
	// None		PT    : Pulse Train (this) copied from PT1

{
  Name    = PT1.Name;			// Copy pulse train name
  PulComposite::operator=(PT1); 	// Copy the composite pulse
  PTCyc   = PT1.PTCyc;			// Copy the pulse train cycle
  PTSCyc  = PT1.PTSCyc;			// Copy the pulse train supercycle
  Cycles  = PT1.Cycles;			// Copy pulse train cycle flag 
  SCycles = PT1.SCycles;		// Copy pulse train supercycle flag 

  return (*this);
}
 

// ____________________________________________________________________________
// C                  CLASS PULSE TRAIN PROPAGATOR FUNCTIONS
// ____________________________________________________________________________

/* These functions allow access and production of propagators while the pulse
   train is active.  Remember, the waveform sets individual steps of the
   composite pulse.  It is cycled through with the phase changes specified
   in the pulse cycle.  The full cycle is cycled through with phase changes
   specified in the pulse supercycle.                                        */  

HSprop PulTrain::GetU(double td)

	// Input		PT 	: A pulse train (this)
	//			td	: A delay time(sec)
        // Output               U	: The propagator for evolving a system
	//				  under the pulse train for time td
 
  {
//                      Determine Evolution Steps
 
  double te = td;                       // Evolution time
  int nFSC  = PTSCyc.fullSCYCs(te);	// Number of full supercycles
  te -= nFSC * PTSCyc.length();		// Time to evolve remaining
  int nFC   = PTCyc.fullcycles(te);	// Number of full cycles
  te -= nFC * PTCyc.length();		// Time to evolve remaining
  int nFWF = fullWFs(te);               // Number of full waveforms
  te -= nFWF*WFtp;			// Time to evolve reamaining
//  int nFWFS = fullsteps(te);            // Number of full waveform steps
//  te -= length(nFWFS);                  // Time to evolve reamaining
 
//                      Build HS Propagators For Steps
//	 	(We Must Also Account For Cycle Phase Changes)
 
  HSprop USCs;				// For full supercycles
  HSprop UCs;				// For full cycles
  HSprop UWFs;				// For full waveforms
  HSprop UWFsteps;			// For waveform steps
  double SCphi;				// for supercycle phase
  double Cphi;				// For cycle phase
  gen_op RZ;

  if(SCycles)				// If supercycles defined
    {
    USCs     = PTSCyc.GetUmult(nFSC);	// Evolve nFSC full supercycles
    UCs      = PTSCyc.GetUsum(nFC);	// Evolve nFC full cycles
    UWFs     = PTCyc.GetUsum(nFWF);	// Evolve nFWF full waveforms
    UWFsteps = PulComposite::GetU(te);	// Evolve waveform steps
    SCphi    = PTSCyc.phase(nFC);	// S.C. phase, last part cycle
    RZ=exp((-complexi*SCphi*DEG2RAD)*FZ());// Rotation op. for this phi
    }
  else if(Cycles)			// If cycles defined
    {
    UCs      = PTSCyc.GetUmult(nFC);	// Evolve nFC full cycles
    UWFs     = PTCyc.GetUsum(nFWF);	// Evolve nFWF full waveforms
    UWFsteps = PulComposite::GetU(te);	// Evolve waveform steps
    Cphi    = PTSCyc.phase(nFC);	// S.C. phase, last part cycle
    RZ = exp((-complexi*Cphi*DEG2RAD)*FZ());// Rotation op. for this phi
    }

  double phi = PTSCyc.phase(nFC);		// S.C. phase, last part cycle
  gen_op Rz = exp((-complexi*phi*DEG2RAD)*FZ());// Rotation op. for this phi

// sosix need to phase UWFs & UWFsteps

//             Build Full HS Propagator For Steps, REVERSE Order
//	 	(These Already Account For Cycle Phase Changes)
 
  HSprop   U  = UWFsteps;		// Evolve over waveform steps
std::cout << "\nSteps & Partial Steps: " << U.length();
  if(nFWF) U *= UWFs;		// Evolve over waveforms
std::cout << "\nFull Waveforms: " << U.length();
  if(nFC)  U *= UCs;		// Evolve over cycles
std::cout << "\nFull Cycles: " << U.length();
  if(nFSC) U *= USCs;		// Evolve over supercycles
std::cout << "\nFull Supercycles: " << U.length();
  return U;           
  }

// ____________________________________________________________________________
//                   CLASS PULSE TRAIN AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 

std::ostream& PulTrain::info(std::ostream& ostr, double td) const

	// Input	PT    : A pulse train (this)
	//		td    : A delay time(sec)
	// Output	ostr  : Information concerning
	//			the pulse train evolution over the
	//			time td are placed into the output stream

  { 
  if(!steps())
    {
    ostr << "\n\n\t\tEmpty Pulse Train, No Evolution Possible\n\n";
    return ostr;
    }
  ostr << "\n\n\t\t";
  if(Name.length()) ostr << Name;
  else              ostr << "\t"; 
  ostr << " Pulse Train Evolution Info\n"; 		// Output Header
  ostr << "\n\tSpecified Evolution Time:          ";
  printTime(ostr, td);
  ostr << "\n\tEvolution Spectral Width:          ";
  double SW = 1.0/td;
  printHz(ostr, SW);

//			     Output Supercycle Info

  if(PTSCyc.length())
    {
    ostr << "\n\tPulse Train Supercycle:            " << PTSCyc.name();
    double tscyc = PTSCyc.length();
    ostr << "\n\tSupercycle Length:                 ";
    printTime(ostr, tscyc);
    double tSW = 1/tscyc;
    ostr << "\n\tSupercycle Spectral Width:         ";
    printHz(ostr, tSW);
    double cpe = td/tscyc;
    ostr << "\n\tSupercycles Spanning Evolution:    " << cpe;
    }

//                              Output Cycle Info

  if(PTCyc.length())
    {
    ostr << "\n\tPulse Train Cycle:                   " << PTCyc.name();
    double tcyc = PTCyc.length();
    ostr << "\n\tCycle Length:                      ";
    printTime(ostr, tcyc);
    double tSW = 1/tcyc;
    ostr << "\n\tCycle Spectral Width:              ";
    printHz(ostr, tSW);
    double cpe = td/tcyc;
    ostr << "\n\tCycles Spanning Evolution:         "
         << Gform("%8.3f", cpe);
    }

//                             Output Waveform Info

  ostr << "\n\tPulse Train Waveform:                " << name();
  ostr << "\n\tWaveform Strength:                 ";
  if(!gamB1const()) ostr << "Variable";
  else              printHz(ostr, strength(0));
  ostr << "\n\tWaveform Length:                   ";
  double twf = length();
  printTime(ostr, twf);
  ostr << "\n\tWaveform Steps:                      " << steps();
  ostr << "\n\tWaveform Steps Spanning Evolution:   " << steps(td);

//                             Output Evoluiton Info

  ostr << "\n\tGAMMA Evolution Steps:\n";
  printEvolve(ostr, td);
  ostr << "\n";
  return ostr;
  }


std::ostream& PulTrain::info(std::ostream& ostr, double td, int npts) const

	// Input	PT    : A pulse train (this)
	//		ostr  : An output stream
	//		td    : Dwell time between FID points
        //              npts  : Number of FID points
        // Output	      : Information regarding the FID generation
	//			is set to the output stream ostr

  {
  if(!steps())
    {
    ostr << "\n\n\t\tEmpty Pulse Train, No Acquisition Possible\n\n";
    return ostr;
    }
  ostr << "\n\n\t\t";
  if(Name.length()) ostr << Name;
  else              ostr << "\t"; 
  ostr << " Pulse Train Acquisition Info\n"; 		// Output Header
  ostr << "\n\tSpecified Evolution Time:";
  ostr << std::string(16, ' ');
  printTime(ostr, td);
  ostr << "\n\tEvolution Spectral Width:";
  ostr << std::string(16, ' ');
  double SW = 1.0/td;
  printHz(ostr, SW);
  ostr << "\n\n"; 
 
  double tacq=0;
  for(int i=0; i<npts; i++)
    {
    printEvolve(ostr, tacq);
    ostr << "\n\tSample"
         << Gdec("%d4", i) << " Acquisition Point Time";
    ostr << std::string(14, ' ');
    printTime(ostr, tacq);
    tacq += td;
    ostr << "\n";
    }
  return ostr;
  }  


// ____________________________________________________________________________
//                   CLASS PULSE TRAIN ACQUISITION FUNCTIONS
// ____________________________________________________________________________

 
row_vector PulTrain::FID(int npts, double td, gen_op &D, gen_op& sigmap)

        // Input        PT    : A pulse train (this)
        //              npts  : Number of FID points
        //              td    : Dwell time between FID points
        //              D     : Detection operator
        //              sigmap: Prepared density operator
        // Output       data  : A row vector contiaing an npts point FID
        //                      that was generated by detection with D
        //                      of the evolution of sigmap under the pulse
        //                      train PT.

  {
  row_vector data(npts);                        // Vector for FID points
  double twf   = length();			// Length of waveform
  double tcyc  = PTCyc.length();		// Length of 1 cycle in wform
  double tscyc = PTSCyc.length();		// Length of 1 supcyc in wform

  int ncycl  = PTCyc.steps();			// Number of cycle steps
  int nscycl = PTSCyc.steps();			// Number of supercycle steps

  int i;
  int nsscycl = int(npts*td/(twf*ncycl*nscycl));// Supercycles in the FID
  HSprop *USSCycl;
  USSCycl = new HSprop[nsscycl];
//  HSprop USSCycl[nsscycl];                      // Connected supercycle props
  USSCycl[0] = PTSCyc.GetU();                   // First for 1 full supercycle
  for(i=1 ; i<nsscycl; i++)
    USSCycl[i] = USSCycl[i-1]*USSCycl[0];
// tREM, NPT
	                             
  int NSC, NC, NWF, NS;				// Various counters
  double phi;					// Phase and time tracking
  HSprop Ut, Urest;				// Needed propagators
  gen_op sigma; 	                     // Needed operators
  gen_op Uz2, Uz3;                              // Needed rotations
 
// sosix
  double te = 0;				// FID evolution time
  for(i=0; i<npts; i++)
    {
//		First Determine Which Propagation Steps Are Required
//		   (For Evolution From FID Start to Current Point)

    NSC = PTSCyc.fullsteps(te);			// Supercycles within te
    te -= NSC * tscyc;				// Adjust evolution time
    NC  = PTCyc.fullsteps(te);			// Cycles within te
    te -= NC * tcyc;				// Adjust evolution time
    NWF = fullWFs(te);				// Waveforms within te
    te -= NWF * twf;				// Adjust evolution time
    NS = fullsteps(te); 			// Waveform steps within te
    te -= sumlength(NS);			// Adjust evolution time

//		Next We Must Build A Propagator In REVERSE Order
//		(And We Must Track Phase Changes Along The Way)

    phi = PTSCyc.phase(NC%nscycl);		// PT phase at FID point (deg)
    Uz2 = exp((-complexi*phi*DEG2RAD)*FZ());    // Rotation op. for this phi
    phi += PTCyc.phase(NWF%ncycl);		// PT phase at FID point (deg)
    Uz3 = exp((-complexi*phi*DEG2RAD)*FZ());    // Rotation op. for this phi
    Urest = HSprop(GetH(NS), te);		// Propagator for partial step
/*
    Ut = Uz3*Urest*adjoint(Uz3);                // First evolve partial step
    if(NS>0)  Ut*= Uz3*GetU(NS-1)*adjoint(Uz3); // Next evolve all full steps
    if(NWF>0)                                   // Next evolve all full trains
      Ut*= Uz2*PTCyc.GetU(NWF-1)*adjoint(Uz2);
    if(NC>0)  Ut*= PTSCyc.GetU(NC-1);           // Next evolve all full cycles
    if(NSC>0) Ut*= USSCycl[NSC-1];              // Next evolve all full supcycs

//            Finally We Evolve Input Sigma To FID Point And Sample

    sigma = evolve(sigmap, Ut);                 // This does the evolution
*/
    data.put(trace(D,sigma),i);                 // Store FID point
    te += td;					// Update FID evolution time
    }
  delete [] USSCycl;
  return data;
  }  

 

row_vector PulTrain::FIDR(int npts, double td, gen_op &D, gen_op& sigmap)

	// Input	PT    : A pulse train (this)
        //              npts  : Number of FID points
	//		td    : Dwell time between FID points
	//		D     : Detection operator
	//		sigmap: Prepared density operator
        // Output	data  : A row vector contiaing an npts point FID
	//			that was generated by detection with D
	//			of the evolution of sigmap under the pulse
	//			train PT.
	// Note		      : Assumes Relaxation Active!

  {
  row_vector data(npts);			// Vector for FID points
//  double tSTP = steplength();			// Length of 1 step in wform
//  double tp = length();				// Length of pulse train
//  double tC = cyclelength();			// Length of 1 cycle in wform
//  double tSC = scyclelength();			// Length of 1 supcyc in wform

//  int ncycl = cycles();				// Number of cycle steps
//  int nscycl = supercycles();			// Number of supercycle steps
//  int nsscycl = int(npts*td/(tp*ncycl*nscycl));	// Supercycles in the FID

//  super_op GSSCycl[nsscycl];                   	// Connected supercycle props
//  GSSCycl[0] = Gsupercycle();			// First for 1 full supercycle
//  for(int i=1 ; i<nsscycl; i++)
//    GSSCycl[i] = GSSCycl[i-1]*GSSCycl[0];

//  int NSC,NC,NPT,NS;                            // Various counters
//  double phi, tREM;                             // Phase and time tracking
//  super_op Gt, Grest; 				// Needed superoperators
//  super_op Lrest, eLrest;			// More needed superops
//  gen_op sigmass;				// For steady state stuff
//  gen_op Uz2, Uz3;                              // Needed rotations
//  super_op Gz2, Gz3;				// Needed superop. rotations
//  gen_op sigma;					// Working density operator
//  double deg2rad = PI/180.0;                    // Conversion factor
 
// sosi: comparisons with previous points could greatly
//       speed this up!  The next point may not need any
//       new supercycle, cycle, train, rotation, ...
//       operators that were not calculated on the previous
//       step.
 
//  for(i=0; i<npts; i++)
//    {
//    NSC = int(i*td/tSC);			// # supercycles to FID point
//    NC = int((i*td-tSC*NSC)/tC);		// # cycles to FID point
//    NPT = int((i*td-tSC*NSC-tC*NC)/tp);		// # pulse trains to FID point
//    NS = int((i*td-tSC*NSC-tC*NC-tp*NPT)/tSTP);	// # individual steps to point
//    tREM = i*td-tSC*NSC-tC*NC-tp*NPT-tSTP*NS;	// Remaining time to point
   
//    phi = scphase(NC%nscycl);			// PT phase at FID point (deg)
//    Uz2 = exp((-complexi*phi*deg2rad)*FZ());	// Rotation op. for this phi
//    Gz2 = U_transform(Uz2);			// Rotation superop. for phi
//    phi += cphase(NPT%ncycl);			// PT phase at FID point (deg)
//    Uz3 = exp((-complexi*phi*deg2rad)*FZ());	// Rotation op. for this phi
//    Gz3 = U_transform(Uz3);			// Rotation superop. for phi
 
//    Lrest = complexi*Hsuper(getH(NS));		// L = -i*[Heff, ] (rad/sec)
//    Lrest += R;					// Full Liouvillian part. step
//    sigmass = sigma_ss(Lrest,R,sigmaeq);	// Steady state dens.op.
//    eLrest = exp(Lrest, -tREM);			// Exp. Liouvillian, part. step
//    Grest = R_prop(eLrest, sigmass);		// Relaxation propagator, jth step
//    Gt = Gz3*Grest;				// First evolve partial step
//    if(NS>0)  Gt*= Gz3*G(NS-1);			// Next evolve all full steps
//    if(NPT>0) Gt*= Gz2*Gcycle(NPT-1);		// Next evolve all full trains
//    if(NC>0)  Gt*= Gsupercycle(NC-1); 		// Next evolve all full cycles
//    if(NSC>0) Gt*= GSSCycl[NSC-1];              // Next evolve all full supcycs
//    sigma = evolve(sigmap, Gt);                 // This does the evolution
//    data.put(trace(D,sigma),i);			// Store FID point
//    }
  return data;
  }  
 


// ____________________________________________________________________________
// Z                    CLASS PULSE TRAIN I/O FUNCTIONS
// ____________________________________________________________________________


std::ostream& PulTrain::printEvolve(std::ostream &ostr, double td, int full) const

	// Input		PT	: Pulse Train
	//			td	: Evolution time
        //                      ostr	: Output stream
	//			full	: Flag for output amount
        // Output               none	: Pul. Train cycle evolution info
	//				  is sent to the output stream

  {
  std::string lstart = "\n\t"; 
  int lgap = 30; 
  int stp = 1;                                  // Evolution Step 
  double tev = 0;                               // Evolution time 
  int plen;                                     // Print length 
 
//		     Output Evolution Under Supercycle(s)

  double tscycle = PTSCyc.length();		// Supercycle length (s)
  int numbsc = 0;				// # of supercycles
  if(tscycle)					// Get # of supercycles
    {
    numbsc = PTSCyc.fullsteps(td);		//	this many scycles
    td -= PTSCyc.length()*numbsc;		//	adjust ev. time
    }
  if(numbsc)					// Output any evolution
    {						// under supercycles
    ostr << lstart << "Step " << stp << "."	//	Evolve step
         << Gdec("%d3", numbsc) << " ";
    ostr << PTSCyc.name() << " Supercycle";   	//	Supercycle name
    plen = (PTSCyc.name()).length() + 11;	//	Add a spacer
    if(numbsc > 1) { ostr << "s"; plen++; }
    ostr << std::string(lgap-plen, ' ');
    printTime(ostr, numbsc*tscycle);		//	Evolve time
    stp++;					//	Adjust step
    tev += numbsc*tscycle;			//	Time evolved
    }

//		        Output Evolution Under Cycle(s)

  double tcycle  = PTCyc.length();		// Cycle length (sec)
  int numbc = 0;				// # of cycles
  if(tcycle)					// Get # of cycle steps
    {
    numbc = PTCyc.fullsteps(td);		//	this many
    td -= PTCyc.length()*numbc;			//	adjust ev. time 
    }
  if(numbc)					// Output any evolution
    {						// under cycles
    ostr << "\n\tStep " << stp << "."
         << Gdec("%d3", numbc) <<  " ";
    ostr << PTCyc.name() << " Cycle";		//	Cycle name
    plen = (PTCyc.name()).length() + 6;		//	Add a spacer
    if(numbc > 1) { ostr << "s"; plen++; }
    ostr << std::string(lgap-plen, ' ');
    printTime(ostr, numbc*tcycle);		//	Evolve time
    stp++;					//	Adjust step 	
    tev += numbc*tcycle;			//	Time evolved
    }

//		        Output Evolution Under Waveform

  int numbw = 0;				// # of waveforms
  while(td >= length())				// Get number of WFs
    {
    td -= length();				//	Evolve this
    numbw++;					//	Adjust WF count
    }
  if(numbw)					// Output any evolution
    {						// Under the waveform steps
    ostr << "\n\tStep " << stp << "."
         << Gdec("%d3", numbw) << " ";
    ostr << name() << " Waveform";		//	Waveform name
    plen = (name()).length() + 8;		//	Add a spacer
    if(numbw > 1) { ostr << "s"; plen++; }
    ostr << std::string(lgap-plen, ' ');
    printTime(ostr, numbw*length());
    stp++;					//	Adjust step
    tev += numbw*length();			//	Time evolved
    }

//		      Output Evolution Under Waveform Steps

  int numb = 0;					// # of waveform steps
  while(td >= length(numb))			// Get number of WF steps
    {
    td -= length(numb);				//	Evolvee this
    numb++;					//	Adjust step count
    }
  if(numb)					// Output any evolution
    {						// Under the waveform steps
    ostr << "\n\tStep " << stp << "."
         << Gdec("%d3", numb) << " ";
    ostr << name() << " Waveform Step";		//	Waveform name
    plen = (name()).length() + 14;		//	Add a spacer
    if(numb > 1) { ostr << "s"; plen++; }
    ostr << std::string(lgap-plen, ' ');
    printTime(ostr, length(numb));
    stp++;					//	Adjust step
    tev += length(numb);			//	Time evolved
    }

//	      Output Evolution Under Partial Waveform Steps

  if(td)			// Get number of WF steps
    {						// Under the waveform steps
    ostr << "\n\tStep " << stp << "."
         << Gdec("%d3", 1) << " ";
    ostr << name() << " Partial Waveform Step ";	//	Waveform name
    plen = (name()).length() + 23;		//	Add a spacer
    ostr << std::string(lgap-plen, ' ');
    printTime(ostr, td);			//	Time evolved
    }


  return ostr;
  }



std::ostream& PulTrain::printCycle(std::ostream &ostr, int full) const

	// Input		PT    : Pulse Train
        //                      ostr  : Output stream
	//			full  : Flag for output amount
        // Output               none  : Pulse Train cycle info is sent
	//				to the output stream

  {
//  if(!PTCyc.steps())
//    {
//    ostr << "\n\tPulse Train Cycle:                None";
//    return ostr;
//    }
  ostr << "\n\tPulse Train Cycle:                " << PTCyc.name();
  PTCyc.printBase(ostr);
  if(full)
    {
    PTCyc.printInfo(ostr);
    PTCyc.printSteps(ostr);
    ostr << "\n";
    }
  return ostr;
  }


std::ostream& PulTrain::printSCycle(std::ostream &ostr, int full) const

	// Input		PT    : Pulse Train
        //                      ostr  : Output stream
	//			full  : Flag for output amount
        // Output               none  : Pulse Train supercycle info is
	//				sent to the output stream

  {
//  if(!PTSCyc.steps())
//    {
    ostr << "\n\tPulse Train SuperCycle:           None";
//    return ostr;
//    }
  return ostr;
  }


std::ostream& PulTrain::print(std::ostream &ostr, int full) const

	// Input		PT    : Pulse Train
        //                      ostr  : Output stream
	//			full  : Flag for output amount
        // Output               none  : Pulse Train info is sent
	//				to the output stream

  {
  if(!steps())
    {
    ostr << "\n\n\t\t\tEmpty Pulse Train\n\n";
    return ostr;
    }
  ostr << "\n\n\t\t\t   Pulse Train "; 		// Output Header
  if(Name.length()) ostr << Name; 
  else              ostr << "Unnamed";
  ostr << "\n\n\tPulse Train Channel:              " << Iso;

  ostr << "\n\tComposite Pulse Name:             " << PulComposite::name();
  PulWaveform::printBase(ostr);	// Print waveform base info
  if(full-1)			// If desired, also print out
    {				// individual waveform steps
    PulWaveform::printSteps(ostr);// and Hamiltonian/propagator
    ostr << "\n";		// status for each step in the
    printInfo(ostr);		// composite pulse
    }

  printCycle(ostr, full-1);	// Print pulse train cycle
  printSCycle(ostr, full-1);	// Print pulse train supercycle 

  ostr << "\n\n";
  return ostr;
  }

            
std::ostream &operator << (std::ostream &ostr, const PulTrain &PT)

	// Input		PT    : Pulse Train
        //                      ostr  : Output stream
        // Output               none  : Pulse Train info is sent
	//				to the output stream

  {
  PT.print(ostr);
  return ostr;
  }
 
#endif						// PulTrain.cc
