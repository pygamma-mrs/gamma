/* PulCycle.cc **************************************************-*-c++-*-
**									**
**      	                G A M M A				**
**									**
**	Class Pulse Cycle			  Implementation	**
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
** This class sets up a pulse cycle for GAMMA.  A pulse cycle consists  **
** of a defined composite pulse (waveform) which is repeated (cycled)   **
** over a specified number of steps with a defined phase.  The class    **
** PulComposite handles aspects of the waveform and its propgators.     **
** This contains a composite pulse, and tracks the pulse cycle and      **
** handles the propagators for spin system evolution as required in NMR **
** simulations.                                                         **
**                                                                      **
** Note: Because this module contains functions that plot the pulse	**
**       cycle in Gnuplot and FrameMaker, this will depend upon the	**
**       GamIO library.							**
**									**
*************************************************************************/

#ifndef PulCycle_cc_			// Is this file already included?
#  define PulCycle_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

# include <Pulses/PulCycle.h>		// Include the header
# include <Pulses/PulAuxil.h>		// Include auxiliary functions 
# include <Pulses/PulComposite.h>	// Know about composite pulses
# include <HSLib/SpinSystem.h>		// Must know about spin systems 
# include <HSLib/HSprop.h>		// Know Hilbert space props
# include <LSLib/LSprop.h>		// Know Liouville space props
# include <GamIO/Ggnuplot.h>		// Include gnuplot output
# include <GamIO/FrameMaker.h>		// Include FrameMaker output
#include <stdlib.h>
# include <string>			// Must know stdlibc++ strings
# include <Basics/StringCut.h>


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                 CLASS PULSE CYCLE ERROR HANDLING
// ____________________________________________________________________________

 
void PulCycle::CYCerror(int eidx, int noret) const
 
        // Input                CYC	: Pulse Cycle (this)
        //                      eidx    : Error flag
        //                      noret   : Return flag 
        // Output               none    : Error Message Output
 

  {
  std::cout << "\nClass Pulse Cycle: ";
  switch(eidx)
    {
    case 0:								// (0)
      std::cout << "Program Aborting....";
      break;
    case 1:								// (1)
      std::cout << "Error During Construction";
      break;
    case 30:								// (30)
      std::cout << "Step Hamiltonian Access, Hamiltonian Does Not Exist";
      break;
    case 31:								// (31)
      std::cout << "Build Step Hamiltonians Before Their Access";
      break;
    case 50:								// (50)
      std::cout << "Evolution For Negative Time Requested";
      break;
    case 51:								// (51)
      std::cout << "Problems In FID Step Timing!  Report Bug Please!";
      break;
    case 60:                                                            // (60)
      std::cout << "Step Synchronized Acquisition With Non-Constant Step Lengths";
      break;
    case 61:                                                            // (61)
      std::cout << "Acquistion Step Synchronization Not Possible!";
      break;
    default:
      std::cout << "Unknown Error (Number " << eidx << ")";
    }
  if(!noret) std::cout << ".\n";
  else       std::cout << ".";
  }
 
 
void volatile PulCycle::CYCfatality(int eidx) const
 
        // Input                CYC	: Pulse Cycle (this)
        //                      eidx    : Error flag
        // Output               none	: Stops execution & error Message
 
  {
  CYCerror(eidx);
  if(eidx) CYCerror(0);
  GAMMAfatal();					// Clean exit from program
  }


// ____________________________________________________________________________
// ii            CLASS PULSE CYCLE DESTRUCTION FACILITATORS 
// ____________________________________________________________________________


void PulCycle::deleteCIndxs()

	// Input	CYC	: A pulse cycle (this)
        // Output       none    : CYC propagator indices
        //                        are deleted if they exist

  {
  if(Pindex) delete [] Pindex;          // Delete current P indices
  Pindex = NULL;			// Set prop indices to NULL
  Pcount = 0;				// Set prop count to NULL
  }


	// Input	CYC	: A pulse cycle (this)
	// Output	none	: CYC step Hilbert propagators
	//			  are deleted if they exist

void PulCycle::deleteCUprops()
  {
  if(CUsteps) delete [] CUsteps;		// Delete current propagators
  CUsteps = NULL;				// Set them to NULL
  if(CUsums) delete [] CUsums;			// Delete current sum propagators
  CUsums = NULL;				// Set them to NULL
  CUcount = 0;					// Zero propagator count
  if(!CGsteps && !CGsums) deleteCIndxs();	// Delete propagator indices
  }

	// Input	CYC	: A pulse cycle (this)
	// Output	none	: CYC step superoperator propagators
	//			  are deleted if they exist

void PulCycle::deleteCGprops()
  {
  if(CGsteps) delete [] CGsteps;		// Delete current superprops
  CGsteps = NULL;				// Set them to NULL
  if(CGsums) delete [] CGsums;			// Delete any sum superprops
  CGsums = NULL;				// Set them to NULL
  CGcount = 0;					// Set prop count to zero
  if(!CUsteps && !CUsums) deleteCIndxs();	// Delete propagator indices
  }


// ____________________________________________________________________________
// iii              CLASS PULSE CYCLE COPY FACILITATORS 
// ____________________________________________________________________________

 
void PulCycle::copyCIndxs(const PulCycle& CYC1)
 
	// Input	CYC	: A pulse cycle (this)
        //              CYC1    : A second pulse cycle
        // Output       none    : CYC1 propagator indices are
        //                        are copied from CYC1
        // Note                 : Assumed array Pindex is currently empty
 
  {
  Pindex = NULL;
  Pcount = CYC1.Pcount;			// Copy the propagator count
  if(CYC1.Pindex)			// If the indices exist
    {
    Pindex = new int[CYCsteps];		//      Allocate space em
    for(int i=0; i<CYCsteps; i++)	//      Copy all indices
      Pindex[i] = CYC1.Pindex[i];
    }  
  }
     


void PulCycle::copyCUprops(const PulCycle& CYC1)

	// Input	CYC	: A pulse cycle (this)
	//		CYC1    : A second PulCycle
	// Output	none	: CYC step Hilbert propagators
	//			  are copied from CYC1
	// Note			: Assumed array CUsteps is empty

  {
  CUsteps = NULL; 
  CUsums = NULL;
  CUcount = 0;
  if(CYC1.CUsteps)			// If the propagators exist
    {
    CUsteps = new HSprop[CYCsteps];	//	Allocate space for them
    for(int i=0; i<CYCsteps; i++)	//	Copy all of the propagators
      CUsteps[i] = CYC1.CUsteps[i];
    }
  if(CYC1.CUsums)			// If sum propagators exist
    {
    CUsums = new HSprop[CYCsteps];	//	Allocate space for them
    for(int i=0; i<CYCsteps; i++)	//	Copy all of the propagators
      CUsums[i] = CYC1.CUsums[i];
    }
  }


void PulCycle::copyCGprops(const PulCycle& CYC1)

	// Input	CYC	: A pulse cycle (this)
	//		CYC1    : A second PulCycle
	// Output	none	: CYC step Liouville propagators
	//			  are copied from CYC1
	// Note			: Assumed array CGsteps is empty

  {
  CGsteps = NULL; 
  CGsums = NULL;
  CGcount = 0;
  if(CYC1.CGsteps)			// If the super propagators exist
    {
    CGsteps = new LSprop[CYCsteps];	//	Allocate space for them
    for(int i=0; i<CYCsteps; i++)	//	Copy all of the superprops
      CGsteps[i] = CYC1.CGsteps[i];
    }
  if(CYC1.CGsums)			// If sum super propagators exist
    {
    CGsums = new LSprop[CYCsteps];	//	Allocate space for them
    for(int i=0; i<CYCsteps; i++)	//	Copy all of the superprops
      CGsums[i] = CYC1.CGsums[i];
    }
  }

    
// ____________________________________________________________________________
// iv        CLASS PULSE CYCLE PROPAGATOR INDEX GENERATORS
// ____________________________________________________________________________
 
 
void PulCycle::SetCIndxs( )
 
	// Input	CYC	: A pulse cycle (this)
        // Output       none    : CYC1 propagator indices are
        //                        are set up
        // Note                 : Assumed array Pindex is currently empty

  {
  if(!CYCsteps) return;			// Do nothing if no cycle steps defined
  if(!Pindex)				// If no cycle step propagators
    Pindex = new int[CYCsteps]; 	// allocate space for them
  int i, j, found;			// Flag if found phase repeat
  double ph = 0;				// Current phase
  Pcount = 0;				// Zero the propagator count
  for(i=0; i<CYCsteps; i++)		// Loop steps & look for repeat phases
    {
    found = 0;                                  // Assume phase not found
    ph = CYCvals.getRe(i);			// Step i phase
    for(j=0; j<i && !found; j++)                // Loop previous steps
      {
      if(ph == CYCvals.getRe(j))		//      If same value then
        {                                       //      flag its been found
        found = 1;                              //      then just store index
        Pindex[i] = j;                          //      where it is kept
        }
      }
    if(!found)                                  // If not found then calculate
      {
      Pindex[i] = i;                            //      Store its index
      Pcount++;					//	Update the counter
      }
    }
  }

    
// ____________________________________________________________________________
// v         CLASS PULSE CYCLE HILBERT SPACE PROPAGATOR GENERATORS
// ____________________________________________________________________________
 
/* These functions allow for the filling up two arrays of propagators for each
   step in the cycle (CUsteps and CUsums).  These are typically NOT generated
   during construction of the pulse cycle but may be generated at any
   time using the combination of the composite pulse & pulse cycle.  Note that
   these are indeed system dependent.                                        */


void PulCycle::SetCUs()

        // Input                CYC	: A pulse cycle (this)
        // Output               void	: Propagators for each step of the
        //                                pulse cycle are calculated
        // Note				: This assumes that the Hamiltonians
        //                                in PC have already been calculated!

  {
  if(!CYCsteps) return;			// Do nothing if no cycle steps defined
  if(!CUsteps)				// If no cycle step propagators
      CUsteps = new HSprop[CYCsteps];	// allocate space for them
  if(!CUsums)				// If no cycle step summed propagators
      CUsums = new HSprop[CYCsteps];	// allocate space for them
  if(!Pindex) SetCIndxs( );		// Set propagator indices if needed
  HSprop Ubase = GetU();		// Begin with pulse waveform propagator
  gen_op RZ;				// Rotation operator for phasing
  double phi=0; 			// Cycle step phase
  int i, idx;
  for(i=0; i<CYCsteps; i++)
    {
    idx = Pindex[i];
    if(i == idx)
      {
      phi = DEG2RAD*CYCvals.getRe(i);	// Cycle step i phase (rad)
      RZ = exp((-complexi*phi)*Fzed);	// Rotation operator for phi
      CUsteps[i] = Ubase.sim_trans(RZ);	// Phase the base propagator
      }
    CUsums[i] = CUsteps[idx];		// Start prop step i full evolve
    if(i) CUsums[i] *= CUsums[i-1];	// Add in previous cycle sums
    }
  }

    
// ____________________________________________________________________________
// vi    CLASS PULSE CYCLE LIOUVILLE SPACE PROPAGATOR GENERATORS
// ____________________________________________________________________________
 
/* These functions allow for the filling up two arrays of propagators for each
   step in the cycle (CGsteps and CGsums).  These are typically NOT generated
   during construction of the pulse cycle but may be generated at any
   time using the combination of the composite pulse & pulse cycle.  Note that
   these are indeed system dependent.                                        */

        // Input                CYC	: A pulse cycle (this)
	//			PW	: A pulse waveform
        // Output               void	: Superpropagators for each step of
	//				  the pulse cycle are calculated
        // Note				: This demands that the Liouville
        //                                propagators have already been made
	//				  in the waveform.

void PulCycle::SetCGs( )
  {
  if(!CYCsteps) return;			// Do nothing if no cycle steps defined
  if(!CGsteps)				// If no cycle step propagators
    CGsteps = new LSprop[CYCsteps];	// allocate space for them
  if(!CGsums)				// If no cycle step summed propagators
    CGsums = new LSprop[CYCsteps];	// allocate space for them
  if(!Pindex) SetCIndxs( );		// Set propagator indices if needed
  LSprop Gbase = GetG();		// Begin with waveform superpropagator
  gen_op RZ;                            // Rotation operator
  super_op GZ;                          // Rotation superoperator
  double phi;                           // Cycle step phase
  int i, idx;
  for(i=0; i<CYCsteps; i++)
    {
    idx = Pindex[i];
    if(i == idx)
      {
      phi = DEG2RAD*CYCvals.getRe(i);	// Cycle step i phase (rad)
      RZ = exp((-complexi*phi)*Fzed);	// Rotation operator for phi
      GZ = U_transform(RZ);		// Rotation superopertor for phi
      CGsteps[i] = GZ*Gbase;		// Phase the base propagator
      }
    CGsums[i] = CGsteps[idx];		// Start prop step i full evolve
    if(i) CGsums[i] *= CGsums[i-1];	// Add in previous cycle sums
    }
  }


// ____________________________________________________________________________
// vii             CLASS COMPOSITE PULSE BASIS HANDLING
// ____________________________________________________________________________
 
/* These functions are only active on propagators.  They will have no effect
   if there have been no propagators generated.                              */
 
void PulCycle::SetBasis(gen_op& Op)
 
        // Input                CYC	: A pulse cycle (this)
        // Output               void    : Hilbert space propagators for all
        //                                steps of the pulse cycle and its
	//				  composite pulse are put
        //                                into the current working basis of Op
 
  {
  SetUBasis(Op);			// Set all of the props
  int i=0;
  if(CUsteps)                           // If step propagators exist
    {                                   // set the basis for each one
    for(i=0; i<CYCsteps; i++)
      CUsteps[Pindex[i]].SetBasis(Op);
    }
  if(CUsums)                            // If sum propagators exist
    {                                   // set the basis for each one
    for(i=0; i<CYCsteps; i++)
      CUsums[i].SetBasis(Op);
    }
  }                 


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A            CLASS PULSE CYCLE CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________


PulCycle::PulCycle() : PulComposite()

	// Input	none	:
	// Output	CYC	: A NULL pulse cycle (this)

  {
  CYCname   = "";			// No pulse cycle name
  CYCsteps  = 0;			// No pulse cycle steps
  CYCtp     = 0;			// No pulse cycle length
  Pindex    = NULL;			// No pulse cycle prop. indices
  Pcount    = 0;			// No pulse cycle prop count
  CUsteps    = NULL;			// No pulse cycle propagators
  CUsums     = NULL;			// No pulse cycle sum propagators
  CGsteps    = NULL;			// No pulse cycle superprops
  CGsums     = NULL;			// No pulse cycle sum superprops
  CUcount    = 0;			// No stored propagators
  CGcount    = 0;			// No stored super propagators
  }


PulCycle::PulCycle(const PulComposite& P, const row_vector& S, std::string N)
         :PulComposite(P)

	// Input	P	: Composite pulse
	// 		S	: Pulse train cycle steps (phase)
	//		N	: Pulse train cycle name
	// Output	CYC	: A new pulse cycle (this)

  {
  CYCsteps   = S.size();		// Set pulse cycle # steps
  CYCvals    = S;			// Set pulse cycle steps
  CYCname    = N;			// Set pulse cycle name
  CYCtp      = WFtp*CYCsteps;		// Set pulse cycle length
  Pindex     = NULL;			// No pulse cycle prop. indices
  Pcount     = 0;			// No pulse cycle prop count
  CUcount    = 0;			// Zero stored propagators
  CUsteps    = NULL;			// Zero pulse cycle propagators
  CUsums     = NULL;			// Zero pulse cycle sum props
  CGcount    = 0;			// No stored super propagators
  CGsteps    = NULL;			// Zero pulse cycle superprops
  CGsums     = NULL;			// Zero pulse train cyc. sum superprops
  }


PulCycle::PulCycle(const PulCycle& CYC1) : PulComposite(CYC1)

	// Input	CYC1	: Pulse Cycle
	// None		CYC	: Pulse Cycle (this), identical
	//			  copy of CYC1

  {
  CYCsteps  = CYC1.CYCsteps;            // Copy the number of steps
  CYCname   = CYC1.CYCname;             // Copy pulse cycle name
  CYCvals   = CYC1.CYCvals;		// Set pulse cycle steps
  CYCtp     = CYC1.CYCtp;		// Copy pulse cycle length
  copyCIndxs(CYC1);			// Copy and propagator indices
  copyCUprops(CYC1);			// Copy any Hilbert Propagators
  copyCGprops(CYC1);			// Copy any Liouville Propagators
  }


// ------------------------------ Destruction ---------------------------------

PulCycle::~PulCycle()

	// Input	CYC	: A pulse cycle (this)
	// Output	none	: CYC is destructed

  { 
  if(Pindex)  delete [] Pindex;		// Delete any current prop indices
  if(CUsteps) delete [] CUsteps;	// Delete any current propagators
  if(CUsums)  delete [] CUsums;		// Delete any current sum propagators
  if(CGsteps) delete [] CGsteps;	// Delete any current superpropagators
  if(CGsums)  delete [] CGsums;		// Delete any current sum superprops
  Pindex     = NULL;			// No pulse cycle prop. indices
  Pcount     = 0;			// No pulse cycle prop count
  CUcount    = 0;			// Zero stored propagators
  CUsteps    = NULL;			// Zero pulse cycle propagators
  CUsums     = NULL;			// Zero pulse cycle sum props
  CGcount    = 0;			// No stored super propagators
  CGsteps    = NULL;			// Zero pulse cycle superprops
  CGsums     = NULL;			// Zero pulse cycle sum superprops
  }

// ------------------------------- Assignment ---------------------------------


PulCycle& PulCycle::operator = (const PulCycle& CYC1)

	// Input	CYC1	: Pulse Cycle
	// None		CYC	: Pulse Cycle (this) copied from CYC1

{
  PulComposite::operator = (CYC1);	// Copy the composite pulse
  CYCsteps  = CYC1.CYCsteps;            // Copy the number of steps
  CYCname   = CYC1.CYCname;             // Copy pulse cycle name
  CYCvals   = CYC1.CYCvals;		// Copy pulse cycle steps
  CYCtp     = CYC1.CYCtp;		// Copy pulse cycle length
  deleteCUprops();			// Delete current propagators
  deleteCGprops();			// Delete current superpropagators
  copyCIndxs(CYC1);			// Copy and propagator indices
  copyCUprops(CYC1);			// Copy any Hilbert Propagators
  copyCGprops(CYC1);			// Copy any Liouville Propagators

  return *this;
}

// ____________________________________________________________________________
// B                CLASS PULSE CYCLE HAMILTONIAN FUNCTIONS
// ____________________________________________________________________________

/* These functions allow accessing average Hamiltonians over the pulse train
   cycle steps.  They are generated from an input pulse waveform.  It is the
   waveform which actually contains individual step Hamiltonians.           */

// sosi
//gen_op PulCycle::GetH(PulComposite& PW, int i) const

        // Input                CYC    : A pulse cycle (this)
        //                      i     : Step in pulse cycle
        // Output               H     : The effective Hamiltonian
	//				for pulse cycle step i

//  {
//  if(!Hsteps)				// If Hamiltonians dont exist
//    {					// then this is a fata error
//    CYCerror(30);
//    CYCfatality(31);
//    }
//  return Hsteps[(Hindex[i])];
//  }


// ____________________________________________________________________________
// C                CLASS PULSE CYCLE PROPAGATOR FUNCTIONS
// ____________________________________________________________________________

/* These functions allow access and production of propagators over the pulse
   train cycle steps.  They are generated from the composite pulse definition
   in conjunction with the pulse cycle phases.  The composite pulse handles
   individual waveform step propagators.                                     */

 
        // Input                CYC	: A pulse cycle (this)
        //                      i       : Step in pulse cycle
        // Output               U       : The propagator for pulse cycle
        //                                that evolves from start of step i
        //                                             to the end of step i
        // Note                         : Steps are indexed [0, nsteps-1)
        // Note                         : An exception is made for default
        //                                then the return will be the
        //                                propagator for the entire pulse
        //                                train cycle

HSprop PulCycle::GetCU(int i)
  {
  if(i<0) return GetCUsum(i);
  if(!CUsteps) SetCUs();		// If props dont exist, build em
  return CUsteps[Pindex[i]];		// Just return the requested U
  }  

/*
 
HSprop PulCycle::GetCU(int i, int j)
 
        // Input                CYC	: A pulse cycle (this)
        //                      i       : Step in pulse cycle
        //                      j       : Step in composite pulse
        // Output               U       : Propagator for pulse cycle
        //                                that evolves from start cycle
	//				  step i, waveform step j to the
        //                                end of that waveform step
        // Note                         : An exception is made for default
        //                                then the return will be the
        //                                propagator for the entire pulse
        //                                waveform
        // Note                         : If requested and not existing,
        //                                this will trigger ALL propators
        //                                to be generated if possible
 
  {
  if(i<0) return GetUsum(i);
  if(!Hsteps)                           // Insure Hamiltonians exist
    {                                   // if not its a fatal error
    CPulerror(10,1);
    CPulerror(17,1);
    CPulfatality(16);
    }
  if(i>=WFsteps)
    {
    CPulerror(40,1);
    CPulfatality(41);
    }
  if(!Usteps) SetUs();                  // If props dont exist, build em
  return Usteps[Uindex[i]];             // Just return the requested U
  }
*/

 
 
HSprop PulCycle::GetCUsum(int i)

        // Input                CYC	: A pulse cycle (this)
        //                      i       : Steps in pulse cycle
        //                                where each step is a waveform
        // Output               U       : Propagator that evolves through
        //                                i steps of the pulse cycle
        // Note                         : If requested and not existing,
        //                                this will trigger ALL propators
        //                                to be generated if possible
        // Note                         : Default propagator is for full
        //                                cycle (last in CUsums array)

  {
  if(i==0)				// Identity if no steps
    return HSprop(Usteps[0].dim());
  if(!CUsums) SetCUs();			// If props dont exist, build em
  if(i<=0) i=CYCsteps;			// Set for total cycle prop
  return CUsums[i-1];			// Just return the requested U
  }
 
 
HSprop PulCycle::GetCUsum(int i, int j)

        // Input                CYC	: A pulse cycle (this)
        //                      i       : Pulse cycle step index
        //                                where each step is a waveform
	//			j	: Waveform step index
        // Output               U       : Propagator that evolves through
        //                                j WF steps of pulse cycle step i

  {
  if(!CUsteps) SetCUs();		// If props dont exist, build em
  HSprop Ubase = GetUsum(j);		// Base propagator over j WF steps
  double phi = DEG2RAD*CYCvals.getRe(i);// Cycle step i phase (rad)
  gen_op RZ = exp((-complexi*phi)*Fzed);// Rotation operator for phi
  return Ubase.sim_trans(RZ);		// Phase the base propagator
  }


HSprop PulCycle::GetCUmult(int N)

        // Input                CYC	: A pulse cycle (this)
        //                      N       : Number of pulse cycles
        // Output               U       : The propagator for N pulse train
        //                                cycles applied in succession

  {
  if(N<=0) return HSprop(CUsteps[0].dim());	// Identity if no steps
  if(!CUsums) SetCUs();				// If no props, build em
  return (CUsums[CYCsteps-1]).Pow(N);
  }

 
/*
 
HSprop PulCycle::GetCU(double td)
 
        // Input                CYC	: A pulse cycle (this)
        //                      td      : Evolution time
        // Output               U       : The propagator that evolves a system
        //                                for time td under the pulse cycle
 
  {
 
//  Determine Props For Cycles(s) Waveform(s), Steps(s), Last Partial Step
 
  int nCY = fullcycles(td);		// Number of complete cycles
  HSprop UCYs = GetCUmult(nCY);		// Propagator over these waveforms
  td -= nCY*CYCtp;			// Time we have left to evolve
  int nWF = fullWFs(td);                // Number of complete waveforms
  HSprop UWFs = GetCUsum(nWF);		// Propagator over these waveforms
  td -= nWF*WFtp;                       // Time we have left to evolve
  int nS = fullsteps(td);               // Number of complete WF steps
  HSprop UStps = GetUsum(nS);		// Propagator over these steps
  td -=  sumlength(nS);                 // Time we have left to evolve
  int j = Hindex[nS];                   // Last (partial) step index

//               Build HS Propagator For Steps, REVERSE Order
 
  HSprop U = HSprop(Hsteps[j],td);      // Prop of step i's evolution
  U *= UStps;                           // Evolve full waveform steps
// sosix - need phasing here
  U *= UWFs;                            // Evolve full cycle steps (waveforms)
  U *= UCYs;                            // Evolve full cycles 
  return U;
  }

*/


// ____________________________________________________________________________
// D       CLASS PULSE CYCLE SUPEROPERATOR PROPAGATOR FUNCTIONS
// ____________________________________________________________________________

/* These functions allow access and production of superoperator propagators
   over the pulse cycle steps.  They are generated from an input pulse
   waveform.  Its the waveform which contains individual step propagators.   */

LSprop PulCycle::GetCG(int i)

        // Input                CYC	: A pulse cycle (this)
        //                      i	: Step in pulse cycle
        //                                Default is last superpropagator -
        //                                i.e. for the entire pulse cycle
        // Output               G	: Superpropagator for pulse cycle
        //                                that evolves from the start of
        //                                the waveform to the end of step i
	// Note				: Assumes relaxation operator set

  {
  if(i<0) return GetCGsum(i);
  if(!CGsteps) SetCGs();		// If props dont exist, build em
  return CGsteps[i];			// Just return the requested G
  }
 
LSprop PulCycle::GetCG(int i, int j, double td)

        // Input                CYC	: A pulse cycle (this)
        //                      i       : Pulse cycle step index
        //                                where each step is a waveform
	//			j	: Waveform step index
	//			td	: Waveform step evolution time
	// Output               G       : Superpropagator that evolves through
        //                                WF step j of pulse cycle step i
	//				  for a time ot td
	// Note				: Indices are:  j=[0,WFsteps-1]
	//						i=[0,CYCsteps-1]

  {
  LSprop Gbase = GetG(j,td);		// Base prop., WF step j for td sec
  double phi = DEG2RAD*CYCvals.getRe(i);// Cycle step i phase (rad)
  gen_op RZ = exp((-complexi*phi)*Fzed);// Rotation operator for phi
  super_op GZ = U_transform(RZ);	// Rotation superopertor for phi
  return GZ*Gbase;			// Phase the base propagator
  }
 
 
LSprop PulCycle::GetCGsum(int i)

        // Input                CYC	: A pulse cycle (this)
        //                      i       : Steps in pulse cycle
        //                                where each step is a waveform
        // Output               G       : Superpropagator that evolves through
        //                                i steps of the pulse cycle
        // Note                         : If requested and not existing,
        //                                this will trigger ALL propators
        //                                to be generated if possible
        // Note                         : Default propagator is for full
        //                                cycle (last in CGsums array)

  {
  if(i==0)				// Identity if no steps
    return LSprop(CGsteps[0].dim());
  if(!CGsteps) SetCGs();		// If props dont exist, build em
  if(i<=0) i=CYCsteps;			// Set for total cycle superprop
  return CGsums[i-1];			// Just return the requested G
  }
 
 
LSprop PulCycle::GetCGsum(int i, int j)

        // Input                CYC	: A pulse cycle (this)
        //                      i       : Pulse cycle step index
        //                                where each step is a waveform
	//			j	: Waveform step index
        // Output               G       : Superpropagator that evolves through
        //                                j WF steps of pulse cycle step i

  {
  if(!CGsteps) SetCGs();		// If props dont exist, build em
  LSprop Gbase = GetGsum(j);		// Base propagator over j WF steps
  double phi = DEG2RAD*CYCvals.getRe(i);// Cycle step i phase (rad)
  gen_op RZ = exp((-complexi*phi)*Fzed);// Rotation operator for phi
  super_op GZ = U_transform(RZ);	// Rotation superopertor for phi
  return GZ*Gbase;			// Phase the base propagator
  }


LSprop PulCycle::GetCGmult(int N)

        // Input                CYC	: A pulse cycle (this)
        //                      N       : Number of pulse cycles
        // Output               G       : The superpropagator for N
        //                                cycles applied in succession

  {
  if(N<=0) return LSprop(CGsteps[0].dim());      // Identity if no steps
  if(!CGsteps) SetCGs();			// If no props, build em
  LSprop GSC = CGsums[CYCsteps-1];		// Prop for 1 cycle
  LSprop GNSC = GSC;                            // Set for 1 cycle
  for(int i=1; i<N; i++) GNSC *= GSC;           // Add in next cycles
  return GNSC;
  }


// ____________________________________________________________________________
// D                 CLASS PULSE CYCLE ACCESS FUNCTIONS
// ____________________________________________________________________________

// ---------------------- Functions Over Full Cycle ---------------------------

int          PulCycle::steps()     const { return CYCsteps; }
int          PulCycle::WF_steps()  const { return WFsteps; }
std::string       PulCycle::name()      const { return CYCname; }
std::string       PulCycle::WF_name()   const { return WFname; }
row_vector   PulCycle::values()    const { return CYCvals; }
row_vector   PulCycle::WF_values() const { return WFvals; }
double       PulCycle::length()    const { return CYCtp; }
double       PulCycle::WF_length() const { return WFtp; }

	// Input	CYC	: A pulse cycle (this)
	// Output	steps	: CYC steps
	// 		name	: CYC name
        //              values  : Array of phi values
	//		length  : CYC length, summed step lengths (sec)
	//		Fz	: The z-rotation operator
	//		PC      : The composite pulse


// ------------------ Functions For Specific Cycle Step ----------------------


complex PulCycle::value(int i)  const { return CYCvals.get(i); }
double  PulCycle::phase(int i)  const { return CYCvals.getRe(i); }

        // Input        CYC     : A pulse cycle (this)
        // Output       steps   : CYC steps
        //              name    : CYC name
        //              values  : Array of phi values

// ____________________________________________________________________________
// E               CLASS PULSE CYCLE AUXILIARY FUNCTIONS
// ____________________________________________________________________________


double PulCycle::steps(double td) const

        // Input        CYC	: A pulse cycle (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of cycle steps needed
        //                        to evolve for time td

  {
  if(td <= 0)         return 0;
  else if(CYCtp == 0) return 0;
  else                return double(steps())*td/CYCtp;
  }

/*

int PulCycle::fullsteps(double td) const

        // Input        CYC	: A pulse cycle (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of full cycle steps 
	//			  that can occur in the time td
	// Note			: For negative time we return
	//			  the total number of steps

  {
  if(td < 0)          steps();			// No. of steps
  else if(CYCtp == 0) return 0;			// No steps
  int ns=-1;					// Assume no steps
  double tstep = CYCtp/CYCsteps;		// Step length
  while(td>=0)
    {
    td -= tstep;
    ns++;
    }    
  return ns;
  }
*/


double PulCycle::cycles(double td) const

        // Input        CYC	: A pulse cycle (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of cycles needed
        //                        to evolve for time td

  {
  if(td <= 0)         return 0;
  else if(CYCtp == 0) return 0;
  else                return td/CYCtp;
  }


int PulCycle::fullcycles(double td) const

        // Input        CYC	: A pulse cycle (this)
        //              td      : An evolution time (sec)
        // Output       steps   : Number of full cycles that can
        //                        occur in the time td

  {
  if(td < 0)          steps();			// Return cycle steps if - time
  else if(CYCtp == 0) return 0;			// Return no cycles if no time
  int ns=-1;					// Set for no cycles
  while(td>=0)					// Loop over cycle lengths
    {
    td -= CYCtp;
    ns++;
    }    
  if(fabs(td) < cutzero) ns++;			// Add one more if remaining
  return ns;					// time is essentially zero
  }


 
void PulCycle::scalegB1(double sf)
 
        // Input        CYC     : A pulse cycle (this)
        //              sf      : A scaling factor 
        // Output       void    : All step field strengths in the waveform
        //                        are multiplied by sf.  The exception are 
        //                        steps of zero length (ideal pulses) 
 
  {
  deleteCUprops();			// Delete current propagators
  deleteCGprops();			// Delete current superpropagators
  PulComposite::scalegB1(sf);		// Scale composite pulse
  }


// ____________________________________________________________________________
// E                 CLASS PULSE CYCLE PLOTTING FUNCTIONS
// ____________________________________________________________________________




        // Input		CYC	: A pulse cycle (this)
	//			split   : Flag to split steps
	//				   0: Don't split apart
	//				   #: Split by #*10%*biggest step
	//			ends    : Flag to add ends
	//				   0: Don't put on ends
	//				   #: Add ends #*.01*cycle length
	//			N       : Number of cycles
        // Output               none	: A pulse cycle (this) plot vector
	//				  is made interactively of the
	//				  RF-intensity vs time.
	// Note				: Ideal pulses are set by a step
	//				  having zero length but having
	//				  an non-zero "gamB1" setting.  In
	//				  such cases "gamB1" is the pulse angle
	// Note				: For plotting purposes, the length of
	//				  an ideal 90 pulse will be taken to 
	//				  be the same as the shortest non-zero
	//				  waveform step.  Others scale with 
	//				  ideal pulse angle.

row_vector PulCycle::IvsT(int split, int ends, int N) const
  {
  if(!CYCtp) return row_vector(); 		// Return empty if no cycle set
  row_vector WF = PulComposite::IvsT(split,0,CYCsteps);	// Single cycle plot 
  int spf = 0;                                  // Assume NOT splitting steps
  split = abs(split);                           // Insure split length >= 0
  if(split) spf++;                              // If splitting set split flag
  int endf = 0;                                 // Assume NO ends plotted
  if(ends) endf++;                              // If ends wanted, set end flag
  double delt, endt;                            // Lengths of splits and ends
  double tstep1 = maxlength();			// Biggest WF step length
  delt = double(split)*0.1*tstep1;              // This will be split length
  endt = double(ends)*0.01*CYCtp;		// This will be end length
  int WFpts = WF.size();			// Plotted pts per cycle step
  int      npts = N*WFpts;			// Number of points to return
  if(endf) npts += 2;                           // More points if ends added
  row_vector data(npts, complex0);		// Return data
  double pttime = 0;				// Start at time 0
  int j=0;                                      // Global point index
  if(ends)                                      // Begin output
    {                                           // If ends to be added
    data.put(complex(pttime,0),j++);            //      Add point 0,0
    pttime += endt;                             //      Next pt   endt,0
    }
  double pt=0, ptgB1=0;

//   Add In 1 Cycle At At Time, Including Breaks Between Cycles If Desired  
//	    Time On The Horizontal Axis, gamB1 On the Vertical Axis
//                   (1*endf + N*cyclepts + 1*endf)

  for(int k=0; k<N; k++)			// Loop over desired cycles
    {
    for(int i=0; i<WFpts; i++)			// Loop cycle step points
      {
      pt    = WF.getRe(i);			// Step i length, cycle step k
      ptgB1 = WF.getIm(i);			// Step i gamB1, cycle step k
      data.put(complex(pttime+pt,ptgB1), j++);	// Add step i start point
      }
    pttime += pt;				// Advance to last point added
    if(split && k!=N-1)				// If spacers desired then add
      pttime += delt;				// some space between cycles
    }

  if(ends)					// If ends to be added
    {						// draw horizontal line
    pttime += endt;				//	Next pt   endt,0
    data.put(complex(pttime,0),j++);		//	Draw to   endt,0
    }
  return data;
  }



	// Input		CYC	: A pulse cycle (this)
	//			split   : Flag to split steps
	//				   0: Don't split apart
	//				   #: Split by #*.1*1st pulse length
	//			ends    : Flag to add ends
	//				   0: Don't put on ends
	//				   #: Add ends length #*1st pulse
	//			N       : Number of waveforms
	//			ph	: An additional phase factor
        // Output               none	: A pulse cycle (this) plot vector
	//				  is made interactively of the
	//				  RF-phase vs time.

row_vector PulCycle::PvsT(int split, int ends, int N, double ph) const
  {
  if(!CYCtp) return row_vector();		// Return empty if no cycle set
  int spf = 0;                                  // Assume NOT splitting steps
  split = abs(split);                           // Insure split length >= 0
  if(split) spf++;                              // If splitting set split flag
  int endf = 0;                                 // Assume NO ends plotted
  if(ends) endf++;                              // If ends wanted, set end flag
  double delt, endt;                            // Lengths of splits and ends
  double tstep1 = maxlength();			// Biggest WF step length
  delt = double(split)*0.1*tstep1;              // This will be split length
  endt = double(ends)*0.01*CYCtp;		// This will be end length

  double phi = ph + CYCvals.getRe(0);		// First cycle step phase
  row_vector WF=PulComposite::PvsT(split,0,1,phi);// First cycle step plot 
  int WFpts = WF.size();			// Plotted pts per cycle step
  int      npts = CYCsteps*WFpts;		// Number of points per cycle
  row_vector cycle(npts, complex0);		// Array for cycle steps

  double pttime = 0;				// Start at time 0
  int i, k, j=0;				// Global point index
  double pt=0, pp=0;				// Point time and phase

//   Add In 1 Waveform At At Time, Including Breaks Between Steps If Desired  
//	    Time On The Horizontal Axis, Phase On the Vertical Axis
//                          (cyclesteps * cyclepts)

  for(k=0; k<CYCsteps; k++)			// Loop over cycle steps
    {
    phi = ph + CYCvals.getRe(k);		// Phase for this cycle step
    WF = PulComposite::PvsT(split,0,1,phi);	// Cycle step k plot points 
    for(i=0; i<WFpts; i++)			// Loop cycle step points
      {
      pt    = WF.getRe(i);			// Step i length, cycle step k
      pp    = WF.getIm(i);			// Step i phase, cycle step k
      cycle.put(complex(pttime+pt,pp), j++);	// Add step i start point
      }
    pttime += pt;				// Advance to last point added
    if(split && k!=CYCsteps-1)			// If spacers desired then add
      pttime += delt;				// some space between cycles
    }


//    Add In 1 Cycle At At Time, Including Breaks Between Cycles If Desired  
//	    Time On The Horizontal Axis, Phase On the Vertical Axis
//                          (cyclesteps * cyclepts)

  int cycpts = cycle.size();			// Points in 1 cycle plot
  npts = N*cycpts;				// Number of points to return
  if(endf) npts += 2;                           // More points if ends added
  row_vector data (npts, complex0);		// Array for cycle steps
  pttime = 0;					// Start at time 0
  j=0;						// Reset point counter
  if(ends)                                      // Begin output
    {                                           // If ends to be added
    data.put(complex(pttime,0),j++);            //      Add point 0,0
    pttime += endt;                             //      Next pt   endt,0
    }
  for(k=0; k<N; k++)				// Loop over cycles
    {
    for(i=0; i<cycpts; i++)			// Loop cycle points
      {
      pt    = cycle.getRe(i);			// Step i length, cycle k
      pp    = cycle.getIm(i);			// Step i phase, cycle k
      pp    = 360.0*fmod(pp, 360.0);		// Keep phase in [0, 360)
      if(fabs(pp-360.0) < 1.e-10) pp=0.0;
      data.put(complex(pttime+pt,pp), j++);	// Add step i start point
      }
    pttime += pt;				// Advance to last point added
    if(split && k!=N-1)				// If spacers desired then add
      pttime += delt;				// some space between cycles
    }
  if(ends)					// If ends to be added
    {						// draw horizontal line
    pttime += endt;				//	Next pt   endt,0
    data.put(complex(pttime,0),j++);		//	Draw to   endt,0
    }
  return data;
  }


//------------------------- Gnuplot Plotting Functions ------------------------


	// Input		CYC	: A pulse cycle (this)
        //                      type    : Type of plot to create
	//				   0: Waveform time vs phase
	//				   1: Waveform time vs gamB1
	//			split   : Flag to split steps
	//				   0: Don't split apart
	//				   #: Split by #*.1*1st pulse length
	//			ends    : Flag to add ends
	//				   0: Don't put on ends
	//				   #: Add ends length #*1st pulse
	//			N       : Number of waveforms
	//			phi	: Added phase (degrees)
        // Output               none	: A pulse cycle (this) plot
	//				  is made interactively using
	//				  Gnuplot.

void PulCycle::GP(int type, int split, int ends, int N, double phi) const
  {
  switch(type)
    {
    case 0:					// Time vs. WF Phase
      {
      std::string Aname = CYCname + "_TvsP.asc";
      std::string Gname = CYCname + "_TvsP.gnu";
      row_vector data = PvsT(split, ends, N, phi);
      GP_xy(Aname, data);
      GP_xyplot(Gname, Aname, -1);
      }
      break;
    case 1:					// Time vs. WF gamB1
      {
      std::string Aname = CYCname + "_TvsB1.asc";
      std::string Gname = CYCname + "_TvsB1.gnu";
      row_vector data = IvsT(split, ends, N);
      GP_xy(Aname, data);
      GP_xyplot(Gname, Aname, -1);
      }
    break;
    default:
      {
      break;
      }
    }
  }


//----------------------- FrameMaker Plotting Functions -----------------------



	// Input		CYC	: A pulse cycle (this)
        //                      type    : Type of plot to create
	//				   0: Waveform time vs phase
	//				   1: Waveform time vs gamB1
	//			split   : Flag to split steps
	//				   0: Don't split apart
	//				   #: Split by #*.1*1st pulse length
	//			ends    : Flag to add ends
	//				   0: Don't put on ends
	//				   #: Add ends length #*1st pulse
	//			N       : Number of waveforms
	//			phi	: Added phase (degrees)
        // Output               none	: A pulse cycle (this) plot
	//				  is output in Framemaker MIF format

void PulCycle::FM(int type, int split, int ends, int N, double phi) const
  {
  switch(type)
    {
    case 0:					// Time vs. WF Phase
      {
      std::string Aname = CYCname + "_TvsP.mif";
      row_vector data = PvsT(split, ends, N, phi);
      FM_xyPlot(Aname, data);
      }
      break;
    case 1:					// Time vs. WF gamB1
      {
      std::string Aname = CYCname + "_TvsB1.mif";
      row_vector data = IvsT(split, ends, N);
      FM_xyPlot(Aname, data);
      }
    break;
    default:
      {
      break;
      }
    }
  }

 
// ____________________________________________________________________________
// G                 CLASS PULSE CYCLE EVOLUTION FUNCTIONS
// ____________________________________________________________________________
 
 

 
        // Input        CYC     : A pulse cycle (this)
        //              SW	: Desired spectral width
        // Output       SWsync  : Spectral width synchronized with the
	//			  pulse cycle length (at best) or
	//			  with the pulse step length (2nd best)
	// Note			: If no synchronization possible the input
	//			  value is just returned.

double PulCycle::FIDsync(double& SW) const 
  {
  double td = 1/SW;				// Dwell time for SW
  int nCYCs = fullcycles(td);			// Full cycles within td
  if(nCYCs > 0) return 1.0/(nCYCs * CYCtp);	// Synch to full cycle
  else
    {
    int nWFs = fullWFs(td);			// Full waveforms with td
    if(nWFs > 0) return 1.0/(nWFs * WFtp);	// Synch to full waveform
    }
  return SW;
  }

int PulCycle::FIDtest(double td, int& nCYs, int& nWFs,
                                                   int& nSTs, double& tr) const
  
        // Input        CPul    : A composite pulse (this) 
        //              td      : Desired dwell time (sec) 
        //              nCYs    : Number of cycles
        //              nWFs    : Number of full waveforms 
        //              nSTs	: Number of full step 
        //              tr      : Remainder time (sec) 
        // Output       Fsync   : Synchronization flag
        //                             >1: Synchronized to Cycle
        //                              1: Synchronized to Waveform 
        //                              0: Not synchronized
        //                              -: Synchronized to Step
        // Note                 : Values of nCYs, nWFs, nSTPs, and tr altered 
 

  {
  tr = td;                                      // Evolve time desired
  if(fabs(tr) < cutzero) tr = 0;                // Zero time if below precision
  nCYs = fullcycles(tr);			// Full cycles within tr
  tr  -= nCYs*CYCtp;				// Time left to evolve
  if(fabs(tr) < cutzero) tr = 0;                // Zero time if below precision
  nWFs = fullWFs(tr);                           // Full waveforms within tr
  tr  -= nWFs*WFtp;                             // Time left to evolve
  if(fabs(tr) < cutzero) tr = 0;                // Zero time if below precision
  nSTs = fullsteps(tr);                         // Additonal full WF steps
  tr  -= sumlength(nSTs);                       // Time left to evolve
  if(fabs(tr) < cutzero) tr = 0;                // Zero time if below precision
  if(!tr)					// If no time remaining it is
    {						// synchronized somehow
    if(!nSTs)					//	Here if WF or CYC synch
      {
      if(!nWFs) return 2;			//      Here if cycle synch.
      else      return 1;			//      Here if waveform synch.
      }
    else if(timeconst()) return -1;		//      Here if step synchron.
    }
  return 0;					//      Here if unsynchronized
  }  


row_vector PulCycle::FIDsynchCYC(int npts, int nCYs,
                                          gen_op &D, gen_op& sigmap, int track)

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

  {
  row_vector data(npts, complex0);              // Vector for FID points
  int HS = sigmap.dim();                        // Working Hilbert space
  HSprop Ut(HS);				// Needed superprop.
  gen_op sigma;                                 // Needed operators
  D.Op_base(sigmap);                            // Set D basis to sigmap's
  if(track)                                     // Output header if tracking
    {
    std::string spacer = std::string(15, ' ');
    std::cout << "\n\n\t" << spacer << "Cycle Synchronized Acquistion Tracking\n";
    std::cout << "\n\t" << "  FID       Cycle      Evolution               Point Values";
    std::cout << "\n\t" << " Point      Count         Time           Real    Imaginary   Norm";
    std::cout << "\n\t" << "-------   ---------   ------------   -------------------------------";
    }
  HSprop UCYs = GetCUsum();			// Prop. for 1 cycle
  for(int j=1; j<nCYs; j++) UCYs *= GetCUsum(); // Set prop for nCYs cycles
  double te = 0.0;                              // FID evolution time
  int iCYs = 0;                                 // Number of cycles
  for(int i=0; i<npts; i++)                     // Loop over FID points
    {                                           // Output point info
    sigma = Ut.evolve(sigmap);                  // This does the evolution
    data.put(trace(D,sigma),i);                 // Store FID point
    if(track)                                   // Output FID point if tracking
      {
      std::cout << "\n\t" << Gdec(i+1,5) << ".   ";
      std::cout << Gdec(iCYs,6) << "      ";
      printTime(std::cout, Ut.length());
      std::cout << "    " << Gform("%8.3f", data.getRe(i))
           << "  "   << Gform("%8.3f", data.getIm(i))
           << "  "   << Gform("%8.3f", norm(data.get(i)));
      std::cout.flush();
      }
    if(i != npts-1)
      {
      iCYs += nCYs;                             // Update cycle count
      te += CYCtp;				// Update FID evolution time
      Ut *= UCYs;                               // Update evolution prop
      }
    }
  return data;
  }


row_vector PulCycle::FIDWFsynch(int npts, int nWFs,
                                          gen_op &D, gen_op& sigmap, int track)

        // Input        CYC     : A pulse cycle (this)
        //              npts    : Number of FID points
        //              nWFs    : Number of waveforms between points
        //              D       : Detection operator
        //              sigmap  : Prepared density operator
        //              track   : Flag for tracking the computation
        // Output       data    : A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under composite
        //                        composite pulse CPul.

  {
  row_vector data(npts, complex0);              // Vector for FID points
  int HS = sigmap.dim();                        // Working Hilbert space
  HSprop Ut(HS);				// Needed propagator
  gen_op sigma;                                 // Needed operators
  D.Op_base(sigmap);                            // Set D basis to sigmap's
  int nCYs = fullcycles(double(nWFs)*WFtp);     // Cycles in 1 increment
  nWFs -= nCYs*CYCsteps;                        // Waveforms left in 1 inc.
  if(!nWFs)                                     // If no WFs, use cycle synch
    return FIDsynchCYC(npts,nCYs,D,sigmap,track);
 
  if(track)                                     // Output header if tracking
    {
    std::string spacer = std::string(13, ' ');
    std::cout << "\n\n\t" << spacer << "Waveform Synchronized Acquistion Tracking\n";
    std::cout << "\n\t" << "  FID       Cycle    Waveform    Evolution               Point Values";
    std::cout << "\n\t" << " Point      Count     Count        Time           Real    Imaginary   Norm";
    std::cout << "\n\t" << "-------   ---------  --------  ------------   -------------------------------";
    }    
 
  HSprop UCYs = Ut;				// Prop. for 0 cycles
  HSprop UCYtot = Ut;				// Prop. for total cycles
  for(int j=0; j<nCYs; j++) UCYs *= GetCUsum(); // Set prop for nCYs cycles
  int iCYs = 0;                                 // Number of cycles
  int iWFs = 0;                                 // Number of waveforms
  for(int i=0; i<npts; i++)                     // Loop over FID points
    {                                           // Output point info
    sigma = Ut.evolve(sigmap);                  // This does the evolution
    data.put(trace(D,sigma),i);                 // Store FID point
    if(track)                                   // Output FID point if tracking
      {
      std::cout << "\n\t" << Gdec(i+1,5) << ".    ";
      std::cout << Gdec(iCYs,6) << "    ";
      std::cout << Gdec(iWFs,6) << "    ";
      printTime(std::cout, Ut.length());
      std::cout << "    " << Gform("%8.3f", data.getRe(i))
           << "  "   << Gform("%8.3f", data.getIm(i))
           << "  "   << Gform("%8.3f", norm(data.get(i)));
      std::cout.flush();
      }
    if(i != npts-1)
      {
      iCYs   += nCYs;				// Update cycle count
      UCYtot *= UCYs;				// Update cycle prop
      iWFs += nWFs;                             // Update the waveform count
      if(iWFs >= CYCsteps)			// If more step than full cycle
        {                                       // then just use full cycle
        iCYs++;                                 //      Add 1 cycle to count
        UCYtot *= GetCUsum();			//      Add 1 cycle to prop
        iWFs -= CYCsteps;			//      Sub 1 cycle WF count
        }
      Ut = UCYtot;
      if(iWFs) Ut *= GetCUsum(iWFs);		// If WFs, update prop
      }
    }
  return data;
  }


row_vector PulCycle::FIDSTsynch(int npts, int nSTs,
                                          gen_op &D, gen_op& sigmap, int track)
 
        // Input        CYC     : A pulse cycle (this)
        //              npts    : Number of FID points
        //              nSTs    : Number of steps between points
        //              D       : Detection operator
        //              sigmap  : Prepared density operator
        //              track   : Flag for tracking the computation
        // Output       data    : A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
        //                        composite pulse CPul.

  {
  if(!timeconst())                              // There can be no step
    {                                           // synchronized acquisition
    CYCerror(60,1);				// if step lengths vary!
    CYCfatality(61);
    }
  row_vector data(npts, complex0);              // Vector for FID points
  int HS = sigmap.dim();                        // Working Hilbert space
  HSprop Ut(HS*HS);                             // Needed propagator
  gen_op sigma;                                 // Needed operators
  D.Op_base(sigmap);                            // Set D basis to sigmap's
  double tstep = PulComposite::length(0);	// Single WF step length

  int nCYs = fullcycles(double(nSTs)*tstep);	// Cycles in 1 increment
  nSTs -= nCYs*CYCsteps;			// Steps left in 1 increment
  int nWFs = fullWFs(double(nSTs)*tstep);	// Waveforms in 1 increment
  nSTs -= nWFs*WFsteps;                         // Steps left in 1 increment
  if(!nSTs)                                     // If no steps, then use either
    {						// cycle synch or WF synch
    if(!nWFs)
      return FIDsynchCYC(npts,nCYs,D,sigmap,track);
    else
      {
      nWFs += CYCsteps*nCYs;
      return FIDWFsynch(npts,nWFs,D,sigmap,track);
      }
    }
 
  if(track)                                     // Output header if tracking
    {
    std::string spacer = std::string(18, ' ');
    std::cout << "\n\n\t" << spacer << "Step Synchronized Acquistion Tracking\n";
    std::cout << "\n\t" << "  FID       Cycle    Waveform     Step     Evolution               Point Values";
    std::cout << "\n\t" << " Point      Count      Count     Count        Time           Real    Imaginary   Norm";
    std::cout << "\n\t" << "-------   --------   ---------   ------   ------------   -------------------------------";
    }
 
  HSprop UCYtot = Ut;                           // Prop. for total cycles
  HSprop UCYs = Ut;				// Prop. for cycle(s) per pt.
  for(int j=0; j<nCYs; j++) UCYs *= GetCUsum();	// Set prop for nCYs cycles
  int iCYs = 0;                                 // Number of cycles
  int iWFs = 0;                                 // Number of waveforms
  int iSTs = 0;                                 // Number of steps
  for(int i=0; i<npts; i++)                     // Loop over FID points
    {                                           // Output point info
    sigma = Ut.evolve(sigmap);                  // This does the evolution
    data.put(trace(D,sigma),i);                 // Store FID point
    if(track)                                   // Output FID point if tracking
      {
      std::cout << "\n\t" << Gdec(i+1,5) << ".   ";
      std::cout << Gdec(iCYs,6) << "      ";
      std::cout << Gdec(iWFs,6) << "      ";
      std::cout << Gdec(iSTs,6) << "      ";
      printTime(std::cout, Ut.length());
      std::cout << "    " << Gform("%8.3f", data.getRe(i))
           << "  "   << Gform("%8.3f", data.getIm(i))
           << "  "   << Gform("%8.3f", norm(data.get(i)));
      std::cout.flush();
      }

   if(i != npts-1)
      {
      iCYs   += nCYs;				// Update cycle count
      iWFs   += nWFs;				// Update waveform count
      iSTs   += nSTs;				// Update step count
      UCYtot *= UCYs;				// Update cycle superprop
      if(iSTs >= WFsteps)			// See if steps than full WF
        {					// then just use full WF
        iWFs++;                                 //      Add 1 WF to count
        iSTs -= WFsteps;                        //      Sub 1 WF from step cnt
        }
      if(iWFs >= CYCsteps)                      // If more WFs than full cycle
        {                                       // then just use full cycle
        iCYs++;                                 //      Add 1 cycle to count
        UCYtot *= GetCUsum();                   //      Add 1 cycle to prop
        iWFs -= CYCsteps;                       //      Sub 1 cycle WF count
        }
      Ut = UCYtot;				// Prop covers full cycles
      if(iWFs) Ut *= GetCUsum(iWFs);            // Add waveforms if neeeded
      if(iSTs) Ut *= GetCUsum(iWFs, iSTs);	// Add steps if needed
      }
    }
  return data;
  }

row_vector PulCycle::FID(int npts, double td,
                                          gen_op &D, gen_op& sigmap, int track)
 
        // Input        CYC     : A pulse cycle (this)
        //              npts	: Number of FID points
        //              td      : Dwell time between FID points
        //              D       : Detection operator
        //              sigmap	: Prepared density operator
	//		track	: Print evolution flag
        // Output       data    : A row vector containing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
        //                        composite pulse CPul.
 
  {
  row_vector data(npts);                        // Vector for FID points
  SetBasis(sigmap);                             // Put props in sigmap basis
  HSprop Ut, Utmp, UI(sigmap.dim());            // Needed propagators
  gen_op sigma;                                 // Needed operators
  double tr, te = 0;                            // FID evolution time
  int nCYs, nWFs, nSTs;
  int synch = FIDtest(td,nCYs, nWFs,nSTs,tr);	// Check FID synchronicity
  if(synch > 1)                                 // Here if Waveform synch.
    return FIDsynchCYC(npts,nWFs,D,sigmap,track);
  else if(synch > 0)                            // Here if Waveform synch.
    {
    int nWFtot = nCYs*CYCsteps + nWFs;
    return FIDWFsynch(npts,nWFtot,D,sigmap,track);
    }
  else if(synch < 0)                            // Here if step synch.
    {
    int nSTtot = nCYs*CYCsteps + nWFs*WFsteps + nSTs;
    return FIDSTsynch(npts,nSTtot,D,sigmap,track);
    }
  if(track)
    {
    std::cout << "\n\n\t                                     Asynchronous Acquistion Tracking\n";
    std::cout << "\n\t" << "  FID    Full     Full     Full   Partial     Partial Step    Propagation             Point Values";
    std::cout << "\n\t" << " Point  Cycles  Waveforms  Steps   Step          Length           Time          Real    Imaginary   Norm";
    std::cout << "\n\t" << "------- ------  ---------  -----  -------     ------------    -----------    -------------------------------";
    }
  double phi=0;
  gen_op RZ;
  for(int i=0; i<npts; i++)
    {
//              First Determine Which Propagation Steps Are Required
//                 (For Evolution From FID Start to Current Point)

    tr   = te;                                  // Evolve time this point
    if(fabs(tr) < cutzero) tr = 0;              // Zero time if below precision
    nCYs = fullcycles(tr);			// Full cycles within tr
    tr  -= nCYs*CYCtp;				// Time left to evolve
    if(fabs(tr) < cutzero) tr = 0;              // Zero time if below precision
    nWFs = fullWFs(tr);                         // Full waveforms within tr
    tr  -= nWFs*WFtp;                           // Time left to evolve
    if(fabs(tr) < cutzero) tr = 0;              // Zero time if below precision
    nSTs = fullsteps(tr);                       // Additonal full WF steps
    tr  -= sumlength(nSTs);                     // Time left to evolve
    if(fabs(tr) < cutzero) tr = 0;              // Zero time if below precision
//  SetBasis(sigmap);                             // Put props in sigmap basis

    if(track)
      {
      std::cout << "\n\t" << Gdec(i+1,5) << ".   ";
      std::cout << Gdec(nCYs,5) << "   ";
      std::cout << Gdec(nWFs,6) << "   ";
      std::cout << Gdec(nSTs,6) << "   ";
      if(tr)
        {
        std::cout << Gdec(nSTs,6) << "    ";
        printTime(std::cout, tr);
        std::cout << "  ";
        }
      else
        std::cout << "                         ";
      std::cout.flush();
      }  
    if(tr < 0)					// Insure evolution time
      {						// is non-negative
      CYCerror(50, 1);
      CYCfatality(51);
      }
    
//               Next We Must Build A Propagator In REVERSE Order
 
             Ut = UI;                           // Begin with identity prop
//std::cout << "\n\t\t\tU: in-" << Ut.length();
    if(tr)   Ut=HSprop(Hsteps[Hindex[nSTs]],tr);// Evolve for partial NSTs step
//std::cout << ", ps-" << Ut.length();
    if(nSTs) Ut*=GetUsum(nSTs);			// Next evolve nSTs steps
//std::cout << ", fs-" << Ut.length();
    phi = DEG2RAD*CYCvals.getRe(nWFs);		// Cycle step nWFs phase (rad)
    RZ = exp((-complexi*phi)*Fzed);		// Rotation operator for phi
    Ut = Ut.sim_trans(RZ); 			// Put proper phase on prop.
    if(nWFs) Ut*=GetCUsum(nWFs);		// Next evolve nWFs waveforms
//std::cout << ", wf-" << Ut.length();
    if(nCYs) Ut*=GetCUmult(nCYs);		// Next evolve nCYs cycles
//std::cout << ", cy-" << Ut.length();
    if(track)
      {
      printTime(std::cout, Ut.length());
      std::cout.flush();
      }
 
//            Finally We Evolve Input Sigma To FID Point And Sample
 
    sigma = Ut.evolve(sigmap);                  // This does the evolution
    data.put(trace(D,sigma),i);                 // Store FID point
    te += td;                                   // Update FID evolution time
    if(track)
      std::cout << "    " << Gform("%8.3f", data.getRe(i))
           << "  "   << Gform("%8.3f", data.getIm(i))
           << "  "   << Gform("%8.3f", norm(data.get(i)));
    }
  return data;
  }  

 
// ------------------------- With Relaxaton & Exchange ------------------------

row_vector PulCycle::FIDRsynchCYC(int npts, int nCYs,
                                          gen_op &D, gen_op& sigmap, int track)

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
        // Note                 : Assumes Relaxation Active!

  {
  double FIDcut = 1.e-6;                        // FID cutoff intensity
  int FIDptcnt = 0;                             // FID zero point count
  int FIDnilcnt = 4;                            // Allowed FID zero points
  row_vector data(npts, complex0);              // Vector for FID points
  int HS = sigmap.dim();                        // Working Hilbert space
  LSprop Gt(HS*HS);                             // Needed superprop.
  gen_op sigma;                                 // Needed operators
  D.Op_base(sigmap);                            // Set D basis to sigmap's
  if(track)                                     // Output header if tracking
    {
    std::string spacer = std::string(1, ' ');
    std::cout << "\n\n\t" << spacer << "Cycle Synchronized Acquistion Tracking, Relaxation/Exchange Active\n";
    std::cout << "\n\t" << "  FID       Cycle      Evolution               Point Values";
    std::cout << "\n\t" << " Point      Count         Time           Real    Imaginary   Norm";
    std::cout << "\n\t" << "-------   ---------   ------------   -------------------------------";
    }
  LSprop GCYs = GetCGsum();			// Prop. for 1 cycle
  for(int j=1; j<nCYs; j++) GCYs *= GetCGsum(); // Set prop for nCYs cycles
  int iCYs = 0;                                 // Number of cycles
  for(int i=0; i<npts; i++)                     // Loop over FID points
    {                                           // Output point info
    sigma = Gt.evolve(sigmap);                  // This does the evolution
    data.put(trace(D,sigma),i);                 // Store FID point
    if(track)                                   // Output FID point if tracking
      {
      std::cout << "\n\t" << Gdec(i+1,5) << ".   ";
      std::cout << Gdec(iCYs,6) << "      ";
      printTime(std::cout, Gt.length());
      std::cout << "    " << Gform("%8.3f", data.getRe(i))
           << "  "   << Gform("%8.3f", data.getIm(i))
           << "  "   << Gform("%8.3f", norm(data.get(i)));
      std::cout.flush();
      }
    if(norm(data.get(i)) < FIDcut) FIDptcnt++;
    else                           FIDptcnt=0;
    if(FIDptcnt > FIDnilcnt)
      {
      if(track)
        std::cout << "\n\n\tFID Computation Into Baseline Noise Level ......";
      return data;
      }
    if(i != npts-1)
      {
      iCYs += nCYs;                             // Update cycle count
      Gt *= GCYs;                               // Update evolution operator
      }
    }
  return data;
  }  
 

row_vector PulCycle::FIDRWFsynch(int npts, int nWFs, 
                                          gen_op &D, gen_op& sigmap, int track)
 
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
        // Note                 : Assumes Relaxation Active!

  {
  double FIDcut = 1.e-6;                        // FID cutoff intensity
  int FIDptcnt = 0;                             // FID zero point count
  int FIDnilcnt = 4;                            // Allowed FID zero points
  row_vector data(npts, complex0);              // Vector for FID points
  int HS = sigmap.dim();                        // Working Hilbert space
  LSprop Gt(HS*HS);				// Needed superprop.
  gen_op sigma;                                 // Needed operators
  D.Op_base(sigmap);                            // Set D basis to sigmap's
  int nCYs = fullcycles(double(nWFs)*WFtp);	// Cycles in 1 increment
  nWFs -= nCYs*CYCsteps;			// Waveforms left in 1 inc.
  if(!nWFs)                                     // If no WFs, use cycle synch
    return FIDRsynchCYC(npts,nCYs,D,sigmap,track);
 
  if(track)                                     // Output header if tracking
    {
    std::string spacer = std::string(4, ' ');
    std::cout << "\n\n\t" << spacer << "Waveform Synchronized Acquistion Tracking, Relaxation/Exchange Active\n";
    std::cout << "\n\t" << "  FID       Cycle    Waveform    Evolution               Point Values";
    std::cout << "\n\t" << " Point      Count     Count        Time           Real    Imaginary   Norm";
    std::cout << "\n\t" << "-------   ---------  --------  ------------   -------------------------------";
    }    
 
  LSprop GCYs = Gt;				// Prop. for 0 cycles
  LSprop GCYtot = Gt;				// Prop. for total cycles
  for(int j=0; j<nCYs; j++) GCYs *= GetCGsum();	// Set prop for nCYs cycles
  int iCYs = 0;                                 // Number of cycles
  int iWFs = 0;                                 // Number of waveforms
  for(int i=0; i<npts; i++)                     // Loop over FID points
    {                                           // Output point info
    sigma = Gt.evolve(sigmap);                  // This does the evolution
    data.put(trace(D,sigma),i);                 // Store FID point
    if(track)                                   // Output FID point if tracking
      {
      std::cout << "\n\t" << Gdec(i+1,5) << ".    ";
      std::cout << Gdec(iCYs,6) << "    ";
      std::cout << Gdec(iWFs,6) << "    ";
      printTime(std::cout, Gt.length());
      std::cout << "    " << Gform("%8.3f", data.getRe(i))
           << "  "   << Gform("%8.3f", data.getIm(i))
           << "  "   << Gform("%8.3f", norm(data.get(i)));
      std::cout.flush();
      }
    if(norm(data.get(i)) < FIDcut) FIDptcnt++;
    else                           FIDptcnt=0;
    if(FIDptcnt > FIDnilcnt)
      {
      if(track)
        std::cout << "\n\n\tFID Computation Into Baseline Noise Level ......";
      return data;
      }
    if(i != npts-1)
      {
      iCYs   += nCYs;				// Update cycle count
      GCYtot *= GCYs;				// Update cycle prop
      iWFs += nWFs;                             // Update the waveform count
      if(iWFs >= CYCsteps)			// If more step than full cycle
        {                                       // then just use full cycle
        iCYs++;                                 //      Add 1 cycle to count
        GCYtot *= GetCGsum();			//      Add 1 cycle to prop
        iWFs -= CYCsteps;			//      Sub 1 cycle WF count
        }
      Gt = GCYtot;
      if(iWFs) Gt *= GetCGsum(iWFs);		// If WFs, update prop
      }
    }    
  return data;
  }

 
row_vector PulCycle::FIDRSTsynch(int npts, int nSTs,
                                          gen_op &D, gen_op& sigmap, int track)
 
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

  {
  if(!timeconst())                              // There can be no step
    {                                           // synchronized acquisition
    CYCerror(60,1);				// if step lengths vary!
    CYCfatality(61);
    }
  double FIDcut = 1.e-6;                        // FID cutoff intensity
  int FIDptcnt = 0;                             // FID zero point count
  int FIDnilcnt = 4;                            // Allowed FID zero points
  row_vector data(npts, complex0);              // Vector for FID points
  int HS = sigmap.dim();                        // Working Hilbert space
  LSprop Gt(HS*HS);                             // Needed superprop.
  gen_op sigma;                                 // Needed operators
  D.Op_base(sigmap);                            // Set D basis to sigmap's
  double tstep = PulComposite::length(0);

  int nCYs = fullcycles(double(nSTs)*tstep);	// Cycles in 1 increment
  nSTs -= nCYs*CYCsteps;			// Steps left in 1 increment
  int nWFs = fullWFs(double(nSTs)*tstep);	// Waveforms in 1 increment
  nSTs -= nWFs*WFsteps;                         // Steps left in 1 increment
  if(!nSTs)                                     // If no steps, then use either
    {						// cycle synch or WF synch
    if(!nWFs)
      return FIDRsynchCYC(npts,nCYs,D,sigmap,track);
    else
      {
      nWFs += CYCsteps*nCYs;
      return FIDRWFsynch(npts,nWFs,D,sigmap,track);
      }
    }
 
  if(track)                                     // Output header if tracking
    {
    std::string spacer = std::string(4, ' ');
    std::cout << "\n\n\t" << spacer << "Step Synchronized Acquistion Tracking, Relaxation/Exchange Active\n";
    std::cout << "\n\t" << "  FID       Cycle    Waveform     Step     Evolution               Point Values";
    std::cout << "\n\t" << " Point      Count      Count     Count        Time           Real    Imaginary   Norm";
    std::cout << "\n\t" << "-------   --------   ---------   ------   ------------   -------------------------------";
    }
 
  LSprop GCYtot = Gt;                           // Prop. for total cycles
  LSprop GCYs = Gt;				// Prop. for cycle(s) per pt.
  for(int j=0; j<nCYs; j++) GCYs *= GetCGsum();	// Set prop for nCYs cycles
  int iCYs = 0;                                 // Number of cycles
  int iWFs = 0;                                 // Number of waveforms
  int iSTs = 0;                                 // Number of steps
  for(int i=0; i<npts; i++)                     // Loop over FID points
    {                                           // Output point info
    sigma = Gt.evolve(sigmap);                  // This does the evolution
    data.put(trace(D,sigma),i);                 // Store FID point
    if(track)                                   // Output FID point if tracking
      {
      std::cout << "\n\t" << Gdec(i+1,5) << ".   ";
      std::cout << Gdec(iCYs,6) << "      ";
      std::cout << Gdec(iWFs,6) << "      ";
      std::cout << Gdec(iSTs,6) << "      ";
      printTime(std::cout, Gt.length());
      std::cout << "    " << Gform("%8.3f", data.getRe(i))
           << "  "   << Gform("%8.3f", data.getIm(i))
           << "  "   << Gform("%8.3f", norm(data.get(i)));
      std::cout.flush();
      }
    if(norm(data.get(i)) < FIDcut) FIDptcnt++;
    else                           FIDptcnt=0;
    if(FIDptcnt > FIDnilcnt)
      {
      if(track)
        std::cout << "\n\n\tFID Computation Into Baseline Noise Level ......";
      return data;
      }

   if(i != npts-1)
      {
      iCYs   += nCYs;				// Update cycle count
      iWFs   += nWFs;				// Update waveform count
      iSTs   += nSTs;				// Update step count
      GCYtot *= GCYs;				// Update cycle superprop
      if(iSTs >= WFsteps)			// See if steps than full WF
        {					// then just use full WF
        iWFs++;                                 //      Add 1 WF to count
        iSTs -= WFsteps;                        //      Sub 1 WF from step cnt
        }
      if(iWFs >= CYCsteps)                      // If more WFs than full cycle
        {                                       // then just use full cycle
        iCYs++;                                 //      Add 1 cycle to count
        GCYtot *= GetCGsum();                   //      Add 1 cycle to prop
        iWFs -= CYCsteps;                       //      Sub 1 cycle WF count
        }
      Gt = GCYtot;				// Prop covers full cycles
      if(iWFs) Gt *= GetCGsum(iWFs);            // Add waveforms if neeeded
      if(iSTs) Gt *= GetCGsum(iWFs, iSTs);	// Add steps if needed
      }
    }
  return data;
  }  


row_vector PulCycle::FIDR(int npts, double td,
                                          gen_op &D, gen_op& sigmap, int track)
 
        // Input        CYC     : A pulse cycle (this)
        //              npts	: Number of FID points
        //              td      : Dwell time between FID points
        //              D       : Detection operator
        //              sigmap	: Prepared density operator
	//		track   : Print evolution flag
        // Output       data	: A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigp under the
        //                        composite pulse CPul.
        // Note			: Assumes Relaxation Active!
 
  {
  if(!sigmaeq.dim())                            // If no relaxation is set
    return FID(npts, td, D, sigmap, track);	// just calculate w/o it
  double tpart=0;				// FID partial step evolution time
  int nCYs, nWFs, nSTs;				// Cycles, WFs, Steps per pt
  int synch = FIDtest(td,nCYs,nWFs,nSTs,tpart);	// Check FID synchronicity
  if(synch > 1)                                 // Here if Waveform synch.
    return FIDRsynchCYC(npts,nWFs,D,sigmap,track);
  else if(synch > 0)                            // Here if Waveform synch.
    {
    int nWFtot = nCYs*CYCsteps + nWFs;
    return FIDRWFsynch(npts,nWFtot,D,sigmap,track);
    }
  else if(synch < 0)                            // Here if step synch.
    {
    int nSTtot = nCYs*CYCsteps + nWFs*WFsteps + nSTs;
    return FIDRSTsynch(npts,nSTtot,D,sigmap,track);
    }

  if(track)
    {
    std::cout << "\n\n                     Asynchronous Acquistion Tracking, Relaxation/Exchange Active\n";
    std::cout << "\n  " << "  FID       Cycle    Waveform     Step       Partial      Evolution               Point Values";
    std::cout << "\n  " << " Point      Count      Count     Count     Step Time         Time           Real    Imaginary   Norm";
    std::cout << "\n  " << "-------   --------   ---------   ------   ------------   ------------   -------------------------------";
    }

  double FIDcut = 1.e-6;                        // FID cutoff intensity
  int FIDptcnt = 0;                             // FID zero point count
  int FIDnilcnt = 4;                            // Allowed FID zero points
  row_vector data(npts, complex0);              // Vector for FID points
  SetBasis(sigmap);				// Put props in sigmap basis
  int HS = sigmap.dim();			// Spin Hilbert space 
  LSprop Gt(HS*HS);				// Needed propagator
  gen_op RZ, sigma;				// Needed operators
  double te=0, tfid=0;				// For FID evolve time

  LSprop GCYtot = Gt;                           // Prop. for total cycles
  LSprop GCYs = Gt;				// Prop. for cycle(s) per pt.
  for(int j=0; j<nCYs; j++) GCYs *= GetCGsum();	// Set prop for nCYs cycles
  int iCYs = 0;                                 // Number of cycles
  int iWFs = 0;                                 // Number of waveforms
  int iSTs = 0;                                 // Number of steps
  for(int i=0; i<npts; i++)
    {
    sigma = Gt.evolve(sigmap);                  // This does the evolution
    data.put(trace(D,sigma),i);                 // Store FID point

    if(track)
      {
      std::cout << "\n" << Gdec(i+1,5) << ".   ";
      std::cout << Gdec(iCYs,6) << "      ";
      std::cout << Gdec(iWFs,6) << "    ";
      std::cout << Gdec(iSTs,6) << "      ";
      printTime(std::cout, te);
      std::cout << " ";
      printTime(std::cout, Gt.length());
      std::cout << "    " << Gform("%8.3f", data.getRe(i))
           << "  "   << Gform("%8.3f", data.getIm(i))
           << "  "   << Gform("%8.3f", norm(data.get(i)));
      std::cout.flush();
      }
    if(te < 0)					// Insure evolution time
      {						// is non-negative
      CYCerror(50, 1);
      CYCfatality(51);
      }

    if(norm(data.get(i)) < FIDcut) FIDptcnt++;
    else                           FIDptcnt=0;
    if(FIDptcnt > FIDnilcnt)
      {
      if(track)
        std::cout << "\n\n\tFID Computation Into Baseline Noise Level ......";
      return data;
      }

//	             This Tracks Evolution To the Next FID Point

   if(i != npts-1)
      {
      tfid += td;
      te    = tfid;				// Total FID evolution time
      iCYs += nCYs;				// Update cycle count
      iWFs += nWFs;				// Update waveform count
      te   -= double(iCYs*CYCsteps+iWFs)*WFtp;	// FID time after CYCs & Wfs
      if(fullWFs(te)) 				// See if another WF needed
        {					// If so add it to counter
        iWFs++;					// and decrement time left
        te -= WFtp;				// to next acquisiton
        }
      iSTs  = fullsteps(te);			// Update step count
      te   -= sumlength(iSTs);			// Update FiD time left

      GCYtot *= GCYs;				// Update cycle superprop
      if(iWFs >= CYCsteps)                      // If more WFs than full cycle
        {                                       // then just use full cycle
        iCYs++;                                 //      Add 1 cycle to count
        GCYtot *= GetCGsum();                   //      Add 1 cycle to prop
        iWFs -= CYCsteps;                       //      Sub 1 cycle WF count
        }
      Gt = GCYtot;				// Prop covers all full cycles
/*      if(iSTs >= WFsteps)			// See if steps than full WF
        {					// then just use full WF
        iWFs++;                                 //      Add 1 WF to count
        iSTs -= WFsteps;                        //      Sub 1 WF from step cnt
        }
      if(PulComposite::length(iSTs)-te		// If partial step longer
                                     <= cutzero)// than step itself, then use
        {					// a full step
        te -= PulComposite::length(iSTs);	// New partial step length
        iSTs++; 				// Adjust full step counter
        if(iSTs >= WFsteps)			// See if now steps than full
          {					// WF, then just use full WF
          iWFs++;				//      Add 1 WF to count
          iSTs -= WFsteps;			//      Sub 1 WF from step cnt
          }
        }
*/
      if(fabs(te) < cutzero) te = 0;		// Zero time if below clock ability
      if(iWFs) Gt *= GetCGsum(iWFs);            // Add waveforms if neeeded
      if(iSTs) Gt *= GetCGsum(iWFs, iSTs);	// Add steps if needed
      if(te)   Gt *= GetCG(iWFs,iSTs,te);	// Add partial step if needed
      }
    }
  return data;
  }  

 
// ____________________________________________________________________________
// Z                  CLASS PULSE CYCLE I/O FUNCTIONS
// ____________________________________________________________________________

 
std::ostream& PulCycle::printEvolve(std::ostream &ostr, double td) const
 
        // Input                CYC	: A pulse cycle (this)
        //                      td      : Evolution time
        //                      ostr    : Output stream
        // Output               none    : Pul. Train cycle evolution info
        //                                is sent to the output stream
 
  {
  std::string lstart = "\n\t";			// Line beginning
  int lgap = 30;				// Mid-spacer
  int stp = 1;                                  // Evolution Step
  double tev = 0;                               // Evolution time
  int plen = 0;					// Print length
 
  ostr << "\n\n\t\t";
  if(CYCname.length()) ostr << CYCname;
  else              ostr << "\t";
  ostr << " Cycle Evolution Info\n";            // Output Header
  ostr << "\n\tSpecified Evolution Time:";
  ostr << std::string(16, ' ');
  printTime(ostr, td);
  ostr << "\n\tEvolution Spectral Width:";
  ostr << std::string(16, ' ');
  double SW = 1.0/td;
  printHz(ostr, SW);

//                      Output Evolution Under Cycle(s)
 

  int nCYCs = fullcycles(td);			// # of cycles
  if(nCYCs)                                     // Output any evolution
    {                                           // under full cycles
    ostr << lstart << "Step " << stp << "."
         << Gdec(nCYCs, 3) <<  " ";
    ostr << CYCname << " Cycle";		//      Cycle name
    plen = CYCname.length() + 6;		//      Add a spacer
    if(nCYCs > 1) { ostr << "s"; plen++; }
    ostr << std::string(lgap-plen, ' ');
    printTime(ostr, nCYCs*CYCtp);		//      Evolve time
    stp++;					//      Adjust step
    tev += nCYCs*CYCtp;				//      Time evolved
    }
 
//                   Output Evolution Under Full Waveform(s)
 
  int nWFs = fullWFs(td-tev);			// # of full waveforms
  if(nWFs)                                      // Output any evolution
    {
    ostr << lstart << "Step " << stp << "."     //      Evolve step
         << Gdec(nWFs, 3) << " ";
    ostr << WFname << " Waveform";              //      Supercycle name
    plen = WFname.length() + 9;                 //      Add a spacer
    if(nWFs > 1) { ostr << "s"; plen++; }
    ostr << std::string(lgap-plen, ' ');
    printTime(ostr, nWFs*WFtp);                 //      Evolve time
    stp++;                                      //      Adjust step
    tev += nWFs*WFtp;                           //      Time evolved
    }

//                    Output Evolution Under Waveform Steps
//                       (For Less Than 1 Full Waveform)

  int nS = fullsteps(td-tev);                   // # of full waveform steps
  if(nS)                                        // Output any evolution
    {                                           // Under the waveform steps
    ostr << "\n\tStep " << stp << "."
         << Gdec(nS, 3) << " ";
    ostr << WFname << " Waveform Step";         //      Waveform name
    plen = WFname.length() + 14;                //      Add a spacer
    if(nS > 1) { ostr << "s"; plen++; }
    ostr << std::string(lgap-plen, ' ');
    printTime(ostr, sumlength(nS));
    stp++;                                      //      Adjust step
    tev += sumlength(nS);                       //      Time evolved
    }

//            Output Evolution Under Partial Waveform Steps

  double tadd;                                              
  if(fabs(td-tev) > cutzero)                    // It there is a partial step
    {                                           // output any evolution
    ostr << "\n\tStep " << stp << "."
         << Gdec(1, 3) << " ";
    ostr << WFname << " Partial Step ";//      Waveform name
    plen = WFname.length() + 14;                //      Add a spacer
    ostr << std::string(lgap-plen, ' ');
    tadd = td - tev;
    printTime(ostr, tadd);                      //      Time evolved
    tev += tadd;
    }
  ostr << "\n\tSummed Evolution Time:";
  ostr << std::string(19, ' ');
  printTime(ostr, tev);
  ostr << "\n\tRemainder Evolution Time:";
  ostr << std::string(16, ' ');
  printTime(ostr, td-tev);
  ostr << "\n";
  return ostr;
  }
 

std::ostream& PulCycle::printFID(std::ostream &ostr, double td, int npts) const

        // Input        CYC	: A pulse cycle (this)
        //              ostr	: An output stream
        //              td 	: Dwell time between FID points
        //              npts	: Number of FID points
        // Output		: Information regarding the FID generation
        //                        is set to the output stream ostr
	// Note			


 
  {
  if(!WFtp)
    {
    ostr << "\n\n\t\tEmpty Cycle, Acquisition Unaffected By Cycle\n\n";
    return ostr;
    }
  ostr << "\n\n\t\t";
  if(CYCname.length()) ostr << CYCname;
  else              ostr << "\t";
  ostr << " Pulse Cycle Acquisition Info\n";            // Output Header
  ostr << "\n\tSpecified Dwell Time:";
  ostr << std::string(20, ' ');
  printTime(ostr, td);
  ostr << "\n\tSpecified Spectral Width:";
  ostr << std::string(16, ' ');
  double SW = 1.0/td;
  printHz(ostr, SW);
  ostr << "\n";
 
  std::string lstart = "\n\t";
  int lgap = 30;
  int stp = 1;                                  // Evolution Step
  double tev = 0;                               // Evolution time
  int plen;                                     // Print length
 
  double tacq=0;
  int nCYs=0, nWFs=0, nS=0;
  for(int i=0; i<npts; i++)
    {
 
//                   Output Evolution Under Full Cycles(s)
    tev = 0;                                    // Zero evolution time
    stp = 1;					// Back to first step
    ostr << lstart << "Point Evolution Time:";
    ostr << std::string(20, ' ');
    printTime(ostr, tacq);
 
    nCYs = fullcycles(tacq);			// # of full cycles
    if(nCYs)                                    // Output any evolution
      {
      ostr << lstart << "Step " << stp << "."   //      Evolve step
           << Gdec(nCYs, 3) << " ";
      ostr << CYCname << " Cycle";		//      Cycle name
      plen = CYCname.length() + 6;		//      Add a spacer
      if(nCYs > 1) { ostr << "s"; plen++; }
      ostr << std::string(lgap-plen, ' ');
      printTime(ostr, nCYs*CYCtp);		//      Evolve time
      stp++;                                    //      Adjust step
      tev += nCYs*CYCtp;			//      Time evolved
      }
 
//                   Output Evolution Under Full Waveform(s)
 
    nWFs = fullWFs(tacq-tev);			// # of full waveforms
    if(nWFs)                                    // Output any evolution
      {
      ostr << lstart << "Step " << stp << "."   //      Evolve step
           << Gdec(nWFs, 3) << " ";
      ostr << WFname << " Waveform";            //      Supercycle name
      plen = WFname.length() + 9;               //      Add a spacer
      if(nWFs > 1) { ostr << "s"; plen++; }
      ostr << std::string(lgap-plen, ' ');
      printTime(ostr, nWFs*WFtp);               //      Evolve time
      stp++;                                    //      Adjust step
      tev += nWFs*WFtp;                         //      Time evolved
      }

//                    Output Evolution Under Waveform Steps
//                       (For Less Than 1 Full Waveform)
 
    nS = fullsteps(tacq-tev);                   // # of full waveform steps
    if(nS)                                      // Output any evolution
      {                                         // Under the waveform steps
      ostr << "\n\tStep " << stp << "."
           << Gdec(nS, 3) << " ";
      ostr << WFname << " Waveform Step";       //      Waveform name
      plen = WFname.length() + 14;              //      Add a spacer
      if(nS > 1) { ostr << "s"; plen++; }
        ostr << std::string(lgap-plen, ' ');
      printTime(ostr, sumlength(nS));
      stp++;                                    //      Adjust step
      tev += sumlength(nS);                     //      Time evolved
      }
 
//                 Output Evolution Under Partial Waveform Steps

     if(fabs(tacq-tev) > cutzero)               // It there is a partial step
       {                                        // output any evolution
       ostr << "\n\tStep " << stp << "."
            << Gdec(1, 3) << " ";
       ostr << WFname << " Partial Step ";      //      Waveform name
       plen = WFname.length() + 14;             //      Add a spacer
       ostr << std::string(lgap-plen, ' ');
       printTime(ostr, tacq-tev);               //      Time evolved
       } 

    ostr << "\n\tSample"
         << Gdec(i,4) << " Acquisition Point Time";
    ostr << std::string(8, ' ');
    printTime(ostr, tacq);
    tacq += td;
    ostr << "\n";
    }
  ostr << "\n";
  return ostr;
  }


std::ostream& PulCycle::printSteps(std::ostream &ostr) const

        // Input                CYC     : Pulse Cycle
        //                      ostr    : Output stream
        //                      full    : Flag for output amount
        // Output               none    : CYC step info is sent
        //                                to the output stream

  {
  std::string startln = "\n\t";
  std::string spacer = "  ";                 // Column separator
  std::string dblsp = "    ";
  ostr << "\n\tCycle Step Phases:\n" << startln;
  int k, st=0, ncol=4;                  // Number of columns
  for(k=0; k<ncol; k++)
    ostr << "Step" << spacer << " Phase " << dblsp;
  ostr << startln;
  for(k=0; k<ncol; k++)
    ostr << "----" << spacer << "-------" << dblsp;
  ostr << startln;
  double sphase;
  for(k=0, st=0; k<CYCsteps; k++,st++)
    {
    if(st == ncol)
      {
      ostr << startln;
      st = 0;
      }
    sphase = CYCvals.getRe(k);			// Step phase
    ostr << Gform("%3i.", k) << spacer;          // Output step number
    ostr << Gform("%7.2f", sphase);              // Output step phase
    ostr << dblsp;                              // Write a column spacer
    }
  return ostr;
  }


std::ostream& PulCycle::printInfo(std::ostream &ostr) const

	// Input		CYC	: Pulse Cycle
        //                      ostr	: Output stream
	//			full	: Flag for output amount
        // Output               none	: CYC storage info is sent
	//				  to the output stream


  {
  ostr << "\n\tCycle Propagators:                ";
  if(CUsteps || CUsums)
    {
    ostr << "Present (";
    if(CUsteps)
      {
      ostr << "steps";
      if(CUsums) ostr << ", ";
      }
    if(CUsums) ostr << "sums";
    ostr << ")";
    }
  else                ostr << "Absent";

  ostr << "\n\tCycle SuperPropagators:           ";
  if(CGsteps || CGsums)
    {
    ostr << "Present (";
    if(CGsteps)
      {
      ostr << "steps";
      if(CGsums) ostr << ", ";
      }
    if(CGsums) ostr << "sums";
    ostr << ")";
    }
  else       ostr << "Absent";
  if(Pindex)			// Delete current P indices
    {
    ostr << "\n\tCycle Propagators Conserved:      ";
    int con = CYCsteps - Pcount;
    ostr << Gdec(con,4);
    }

  return ostr;
  }


std::ostream& PulCycle::printBase(std::ostream &ostr) const

	// Input		CYC	: Pulse Cycle
        //                      ostr	: Output stream
        // Output               none	: Pulse Cycle base info is
	//				  sent to the output stream

  {
  ostr << "\n\tCycle Steps:                      "
       << Gdec(CYCsteps,4);
  if(CYCtp)
    {
    ostr << "\n\tCycle Length:                     ";
    printTime(ostr, CYCtp);
    ostr << "\n\tCycle Spectral Width:             ";
    double SW = 1.0/CYCtp;
    printHz(ostr, SW);
    ostr << "\n\tCycle Step Length:                ";
    double cslen = CYCtp/CYCsteps;
    printTime(ostr, cslen);
    ostr << "\n\tCycle Step Spectral Width:        ";
    double cssw = 1.0/cslen;
    printHz(ostr, cssw);
    }
  return ostr;
  }


std::ostream& PulCycle::print(std::ostream &ostr, int full) const

	// Input		CYC	: Pulse Cycle
        //                      ostr	: Output stream
	//			full	: Flag for output amount
        // Output               none	: Pulse Cycle info is sent
	//				  to the output stream

  {
  if(!CYCsteps)
    {
    ostr << "\n\n\t\t\tEmpty Pulse Cycle\n\n";
    return ostr;
    }
  ostr << "\n\n\t\t\t";
  if(!full) ostr << "\t";
  ostr << "Pulse Cycle " << CYCname << "\n";
  ostr << "\n\tPulse Channel:";
  ostr << std::string(24-(Iso.length()), ' ');
  ostr << Iso;
  printBase(ostr);			// Print base info
  ostr << "\n\tComposite Pulse Waveform Name:    " << PulComposite::name();
  if(full)
    {
    printInfo(ostr);			// Print prop info
    printSteps(ostr);			// Print cycle steps
    ostr << "\n";
    PulComposite::printBase(ostr);	// Print composite pulse
    if(full > 1)
      {
      PulComposite::printSteps(ostr, full);
      ostr << "\n";
      PulComposite::printInfo(ostr);
      }
    }
  ostr << "\n\n"; 
  return ostr;
  }

            
std::ostream &operator << (std::ostream &ostr, const PulCycle &CYC)

	// Input		CYC	: Pulse Cycle
        //                      ostr	: Output stream
        // Output               none	: Pulse Cycle info is sent
	//				  to the output stream

  {
  CYC.print(ostr);
  return ostr;
  }
 
#endif						// PulCycle.cc
