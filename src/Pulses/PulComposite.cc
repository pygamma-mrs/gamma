/* PulComposite.cc ************************************************-*-c++-*-
**									**
**      	                G A M M A				**
**									**
**	Class Composite Pulse			   Implementation	**
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
** This class handles composite pulses for use in GAMMA.   A composite  **
** pulse will contain a composite pulse (composite pulse steps) and can	**
** generate bot step and pulse Hamilonians & propagators.               **
**                                                                      **
*************************************************************************/

#ifndef _PulComposite_cc_		// Is this file already included?
#define _PulComposite_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Pulses/PulComposite.h>	// Include the header
#include <Pulses/PulAuxil.h>		// Include auxiliary functions
#include <Pulses/PulWaveform.h>		// Know about pulse waveforms
#include <HSLib/SpinSystem.h>		// Must know about spin systems 
#include <HSLib/HSham.h>		// Include Ho function
#include <HSLib/SpinOpCmp.h>		// Include FX,FY,FZ funcitons
#include <HSLib/GenOp.h>		// Include operators
#include <LSLib/DensOp.h>		// Include density operators
#include <Basics/StringCut.h>
#include <string>			// Must know libstdc++ strings
#include <stdlib.h>

double CP_TCUT = 1.e-11;		// Time precision

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS COMPOSITE PULSE ERROR HANDLING
// ____________________________________________________________________________

 
void PulComposite::CPulerror(int eidx, int noret) const
 
        // Input                CPul	: Composite Pulse (this)
        //                      eidx    : Error flag
        //                      noret   : Return flag 
        // Output               none    : Error Message Output

  {
  std::cout << "\nClass Composite Pulse: ";
  switch(eidx)
    {
    case 0:								// (0)
      std::cout << "Program Aborting....";
      break;
    case 1:								// (1)
      std::cout << "Error During Construction";
      break;
    case 2:								// (2)
      std::cout << "No Pulse Waveform Defined";
      break;
    case 3:								// (3)
      std::cout << "Cannot Generate Step Hamiltonians";
      break;
    case 4:								// (4)
      std::cout << "No Pulse Channel Defined";
      break;
    case 10:								// (10)
      std::cout << "Waveform Step Hamiltonians NOT Present";
      break;
    case 11:								// (11)
      std::cout << "Must Build Step Hamiltonians Before Access Allowed";
      break;
    case 15:								// (15)
      std::cout << "Step HS Propagator Build Tried Without Available Hamiltonians";
      break;
    case 16:								// (16)
      std::cout << "Build Step Hamiltonians Before Requesting HS Propagators";
      break;
    case 17:								// (17)
      std::cout << "HS Step Propagator Access Without Available Hamiltonians";
      break;
    case 18:								// (18)
      std::cout << "Summed HS Step Propagator Get Without Available Hamiltonians";
      break;
    case 20:								// (20)
      std::cout << "Step LS Propagator Build Tried Without Available Hamiltonians";
      break;
    case 21:								// (21)
      std::cout << "Build Step Hamiltonians Before Requesting LS Propagators";
      break;
    case 22:								// (22)
      std::cout << "Summed LS Step Propagator Get Without Available Hamiltonians";
      break;
    case 30:								// (30)
      std::cout << "System Contains No Isotopes Affected On Chosen Pulse Channel";
      break;
    case 40:								// (40)
      std::cout << "Requested Propagator For Non-Existant Waveform Step";
      break;
    case 41:								// (41)
      std::cout << "Waveform Steps Span [0, " << WFsteps << ")";
      break;
    case 50:                                                            // (50)
      std::cout << "Evolution For Negative Time Requested";
      break;
    case 51:                                                            // (51)
      std::cout << "Problems In FID Step Timing!  Report Bug Please!";
      break;
    case 60:                                                            // (60)
      std::cout << "Step Synchronized Acquisition With Non-Constant Step Lengths";
      break;
    case 61:                                                            // (61)
      std::cout << "Acquistion Step Synchronization Not Possible!";
      break;
    case 62:                                                            // (62)
      std::cout << "Full Waveform Synchronization Enabled";
      break;
    case 63:                                                            // (63)
      std::cout << "Full Waveform Step Synchronization Enabled";
      break;
    case 65:                                                            // (65)
      std::cout << "Fraction Synchronized Acquisition With < 1 Fraction";
      break;
    default:
      std::cout << "Unknown Error (Number " << eidx << ")";
    }
  if(!noret) std::cout << ".\n";
  else       std::cout << ".";
  }
 
 
void volatile PulComposite::CPulfatality(int eidx) const
 
        // Input                CPul	: Composite Pulse (this)
        //                      eidx    : Error flag
        // Output               none	: Stops execution & error Message
 
  {
  CPulerror(eidx);
  if(eidx) CPulerror(0);
  std::cout << "\n";
  exit(-1);
  }


// ____________________________________________________________________________
// ii            CLASS COMPOSITE PULSE DESTRUCTION FACILITATORS 
// ____________________________________________________________________________


void PulComposite::deleteHams()

	// Input	CPul	: A composite pulse (this)
	// Output	none	: CPul step Hilbert operators
	//			  are deleted if they exist
	// Note			: Also deletes Hindex!

  {
  if(Hsteps) delete [] Hsteps;		// Delete current Hamiltonians
  Hsteps = NULL;			// Set Hamiltonians to NULL
  deleteHIndxs();
  }


void PulComposite::deleteHIndxs()

	// Input	CPul	: A composite pulse (this)
        // Output       none    : CPul Hamiltonian indices
        //                        are deleted if they exist

  {
  if(Hindex) delete [] Hindex;          // Delete current H indices
  Hindex = NULL;                        // Set H indices to NULL
  Hcount = 0;                           // Set H count to NULL
  }


void PulComposite::deleteUIndxs()

	// Input	CPul	: A composite pulse (this)
        // Output       none    : CPul propagator indices
        //                        are deleted if they exist

  {
  if(Uindex) delete [] Uindex;          // Delete current U indices
  Uindex = NULL;                        // Set prop indices to NULL
  Ucount = 0;                           // Set prop count to NULL
  }


void PulComposite::deleteUprops()

	// Input	CPul	: A composite pulse (this)
	// Output	none	: CPul step Hilbert propagators
	//			  are deleted if they exist
	// Note			: Also deletes Uindex (if no Gs)!

  {
  if(Usteps) delete [] Usteps;		// Delete current propagators
  if(Usums)  delete [] Usums;		// Delete current sum propagators
  Usteps = NULL;			// Set step propagators to NULL
  Usums = NULL;				// Set sum propagators to NULL
  if(!Gsteps) deleteUIndxs();		// Delete prop indices if needed
  }


void PulComposite::deleteLOps()

	// Input	CPul	: A composite pulse (this)
	// Output	none	: CPul step Liouville space operators
	//			  are deleted if they exist

  {
  if(Lsteps) delete [] Lsteps;		// Delete Hamiltonian superops
  Lsteps = NULL;			// Set Hamiltonians to NULL
  }


void PulComposite::deleteGprops()

	// Input	CPul	: A composite pulse (this)
	// Output	none	: CPul step superoperator propagators
	//			  are deleted if they exist
	// Note			: Also deletes Uindex (if no Us)!

  {
  if(Gsteps) delete [] Gsteps;		// Delete current superpropagators
  if(Gsums)  delete [] Gsums;		// Delete current sum propagators
  Gsteps = NULL;			// Set step superprops to NULL
  Gsums  = NULL;			// Set sum superprops to NULL
  if(!Usteps) deleteUIndxs();		// Delete prop indices if needed
  }


void PulComposite::deleteSSs()

	// Input	CPul	: A composite pulse (this)
	// Output	none	: CPul step steady-state operators
	//			  are deleted if they exist
	// Note			: Also deletes Sindex!

  {
  if(SigSSs) delete [] SigSSs;		// Delete steady state props
  SigSSs=NULL;
  }


// ____________________________________________________________________________
// iii              CLASS COMPOSITE PULSE COPY FACILITATORS 
// ____________________________________________________________________________


void PulComposite::copyHams(const PulComposite& CPul1)

	// Input	CPul	: A composite pulse (this)
	//		CPul1    : A second PulComposite
	// Output	none	: CPul step Hamiltonians & their indices
	//			  are copied from CPul1
	// Note			: Assumed arrays Hsteps and Hindex
	//			  are currently empty
	// Note			: Recall that some Hsteps[n] will be empty
	//			  as the actual operator for n will reside
	//			  in Hsteps[Hindex[n]]!

  {
  Hsteps = NULL; 
  if(CPul1.Hsteps)			// If the Hamiltonians exist
    {
    Hsteps = new gen_op[WFsteps];	//	Allocate space for Hs
    for(int i=0; i<WFsteps; i++)	//	Copy all of the Hamiltonians
      Hsteps[i] = CPul1.Hsteps[i];	//	(Some may be NULL)
    Hcount = CPul1.Hcount;		//	Copy Hamiltonian count 
    }
  if(!Hindex) copyHIndxs(CPul1);	// Copy H indices if needed
  }


void PulComposite::copyHIndxs(const PulComposite& CPul1)

        // Input        CPul	: A composite pulse (this)
        //              CPul1	: A second composite pulse
        // Output       none    : CPul1 Hamiltonian indices are
        //                        are copied from CPul1
        // Note                 : Assumed array Hindex is empty

  {
  Hindex = NULL;
  Hcount = CPul1.Hcount;		// Copy the Hamiltonian count
  if(CPul1.Hindex)			// If the indices exist
    {
    Hindex = new int[WFsteps];		//      Allocate space em
    for(int i=0; i<WFsteps; i++)	//      Copy all indices
      Hindex[i] = CPul1.Hindex[i];
    }
  }


void PulComposite::copyUIndxs(const PulComposite& CPul1)

        // Input        CPul	: A composite pulse (this)
        //              CPul1	: A second composite pulse
        // Output       none    : CPul1 propagator indices are
        //                        are copied from CPul1
        // Note                 : Assumed array Uindex is empty

  {
  Uindex = NULL;
  Ucount = CPul1.Ucount;		// Copy the propagator count
  if(CPul1.Uindex)			// If the indices exist
    {
    Uindex = new int[WFsteps];		//      Allocate space em
    for(int i=0; i<WFsteps; i++)	//      Copy all indices
      Uindex[i] = CPul1.Uindex[i];
    }
  }


void PulComposite::copyUprops(const PulComposite& CPul1)

	// Input	CPul	: A composite pulse (this)
	//		CPul1	: A second PulComposite
	// Output	none	: CPul step Hilbert propagators
	//			  are copied from CPul1
	// Note			: Assumed array Usteps is empty
	//			  as should be Usums
	// Note			: Will copy propagator indices
	//			  if currently not existing

  {
  Usteps = NULL; 
  Usums  = NULL;
  if(CPul1.Usteps)			// If the propagators exist
    {
    Usteps = new HSprop[WFsteps];	//	Allocate space for props
    for(int i=0; i<WFsteps; i++)	//	Copy all of the propagators
      Usteps[i] = CPul1.Usteps[i];
    }
  if(CPul1.Usums)			// If sum propagators exist
    {
    Usums = new HSprop[WFsteps];	//	Allocate space for them
    for(int i=0; i<WFsteps; i++)	//	Copy all of the propagators
      Usums[i] = CPul1.Usums[i];
    }
  if(!Uindex) copyUIndxs(CPul1);	// Copy prop indices if needed
  }


void PulComposite::copyLOps(const PulComposite& CPul1)

	// Input	CPul	: A composite pulse (this)
	//		CPul1    : A second PulComposite
	// Output	none	: CPul step Hamiltonians superoperators
	//			  are copied from CPul1
	// Note			: Assumed arrays Lsteps is currently
	//			  currently empty but Hindex is filled
	// Note			: Recall that some Lsteps[n] may be empty
	//			  as the actual operator for n will reside
	//			  in Lsteps[Hindex[n]]!

  {
  Lsteps = NULL; 
  if(CPul1.Lsteps)			// If the Hamiltonians exist
    {
    Lsteps = new super_op[WFsteps];	//	Allocate space for Hs
    for(int i=0; i<WFsteps; i++)	//	Copy Hamiltonian superops
      Lsteps[i] = CPul1.Lsteps[i];	//	(Some may be NULL)
    }
  }


void PulComposite::copyGprops(const PulComposite& CPul1)

	// Input	CPul	: A composite pulse (this)
	//		CPul1    : A second PulComposite
	// Output	none	: CPul step Liouville propagators
	//			  are copied from CPul1
	// Note			: Assumed array Gsteps is empty
	//			  as should be Gsums
	// Note			: Will copy propagator indices
	//			  if currently not existing

  {
  Gsteps = NULL; 
  Gsums = NULL;
  if(CPul1.Gsteps)			// If the super propagators exist
    {
    Gsteps = new LSprop[WFsteps];	//	Allocate space for them
    for(int i=0; i<WFsteps; i++)	//	Copy all of the superprops
      Gsteps[i] = CPul1.Gsteps[i];
    }
  if(CPul1.Gsums)			// If sum superpropagators exist
    {
    Gsums = new LSprop[WFsteps];	//	Allocate space for them
    for(int i=0; i<WFsteps; i++)	//	Copy all of the propagators
      Usums[i] = CPul1.Usums[i];
    }
  if(!Uindex) copyUIndxs(CPul1);	// Copy prop indices if needed
  }


void PulComposite::copySSs(const PulComposite& CPul1)

	// Input	CPul	: A composite pulse (this)
	//		CPul1	: A second PulComposite
	// Output	none	: CPul step steady state operators 
	//			  are copied from CPul1
	// Note			: Assumed array SigSSs is currently empty
	// Note			: Recall that some SigSSs[n] will be empty
	//			  as the actual operator for n will reside
	//			  in SigSSs[Hindex[n]]!

  {
  SigSSs = NULL; 
  if(CPul1.SigSSs)			// If the steady state ops exist
    {
    SigSSs = new densop[WFsteps];	//	Allocate space for SSs
    for(int i=0; i<WFsteps; i++)	//	Copy all of the Hamiltonians
      SigSSs[i] = CPul1.SigSSs[i];	//	(Some may be NULL)
    }
  }


// ____________________________________________________________________________
// iv       CLASS COMPOSITE PULSE HAMILTONIAN & PROPAGATOR INDICES
// ____________________________________________________________________________
 

void PulComposite::SetHIndxs( )

        // Input        CPul	: A pulse cycle (this)
        // Output       none    : CPul Hamiltonian indices are
        //                        are set up
        // Note                 : Assumed array Hindex is currently empty
 
  {
  if(!WFsteps) return;			// Do nothing if no steps defined
  if(!Hindex) Hindex = new int[WFsteps];// If no prop. indices, allocate space
  int i, j, found;			// Flag if found {gamB1, phase} repeat
  Hcount = 0;                           // Zero the propagator count
  for(i=0; i<WFsteps; i++)		// Loop steps & look for repeat Hs
    {
    found = 0;				// 	Assume Hamiltonian step unique
    for(j=0; j<i && !found; j++)	//	Loop all previous steps
      {
      if(WFvals.get(i) == WFvals.get(j))//      If same {gamB1, phase} then
        {				//      repeat H! Flag its been found 
        found = 1; 			//	and store index where it is kept
        Hindex[i] = j;
        }
      else if(!WFvals.getRe(i) && 	//	If same gamB1=0, then things
                       !WFvals.getRe(i))//	match no matter what phase &
        {				//      repeat H! Flag its been found 
        found = 1; 			//	and store index where it is kept
        Hindex[i] = j;
        }
      }
    if(!found)				// If not found then, H is unique
      {
      Hindex[i] = i;			//      Store its index
      Hcount++;				//      Update the H counter
      }
    }
  }
 

void PulComposite::SetUIndxs( )

        // Input        CPul	: A pulse cycle (this)
        // Output       none    : CPul propagator indices are
        //                        are set up
        // Note                 : Assumed array Uindex is currently empty
	// Note			: Assumed Hindex is current
 
  {
  if(!WFsteps) return;			// Do nothing if no steps defined
  if(!Uindex) Uindex=new int[WFsteps];	// If no prop. indices, allocate space
  int i, j, k, found;			// Flag if found phase repeat
  double te = 0;			// Current evolve time
  Ucount = 0;                           // Zero the propagator count
  for(i=0; i<WFsteps; i++)		// Loop steps & look for repeat Hs
    {
    found = 0;				// 	Assume propagator step unique
    j = Hindex[i];			// 	Hamiltonian index, this step
    if(i != j)				//	If Hamiltonian is repeated
      {
      te = WFtimes.getRe(i);		//	Current step i length
      for(k=0; k<i && !found; k++)	//	Loop all previous steps
        {				//	sharing this Hamiltonian
        if(j==Hindex[k]
             && te == WFtimes.getRe(k))	//      If same Hamitonian and the 
          {				//      same evolution time then
          found = 1;			//	repeate U! Flag its been found 
          Uindex[i] = k;		//      and store index where it is kept
          }
        }
      }
    if(!found)				// If not found then, U is unique
      {
      Uindex[i] = i;			//      Store its index
      Ucount++;				//      Update the prop counter
      }
    }
  }    
 

// ____________________________________________________________________________
// v              CLASS COMPOSITE PULSE HAMILTONIAN GENERATORS
// ____________________________________________________________________________

/* These functions allow for the filling up an array of Hamiltonians (Hsteps)
   for each step in the waveform.  These are typically generated during
   construction of the composite pulse as they are spin system dependent.
   Equivalent steps in a waveform will share the same Hamiltonians as indexed
   in the integer array Hindex.                                              */
 
void PulComposite::SetHs(const spin_system& sys)

        // Input                CPul    : A composite pulse (this)
        //                      sys     : A spin system
        // Output               void    : Effective Hamiltonians for all of the
        //                                composite pulse steps are generated
        // Note                         : Will delete any stored propagators
        //                                propagators for the waveform
        // Note                         : Assumes H is Ho and that it needs
        //                                F(X,Y,Z) for offsets and phasing.
        // Note                         : ASSUMES PULSE CHANNEL SET

  {
  gen_op H = Ho(sys);				// The isotropic Hamiltonian
  gen_op FX = Fx(sys,Iso);			// Spin operators for the
  gen_op FY = Fy(sys,Iso);			// rf-field on the axes
  gen_op FZ = Fz(sys,Iso);			// Spin operator for phasing
  SetHs(H, FX, FY, FZ);				// Use overload
  }


void PulComposite::SetHs(gen_op& H, gen_op& FX, gen_op& FY, gen_op& FZ)

        // Input                CPul	: A composite pulse (this)
        //                      H	: Active isotropic Hamiltonian
        //                      FX	: Active Fx operator 
        //                      FY	: Active Fy operator 
        //                      FZ	: Active Fz operator 
        // Output               void	: Effective Hamiltonians for all of the
        //                                composite pulse steps are generated
	// Note				: This will delete any stored
	//				  propagators for the waveform
	// Note				: Uses F(X,Y,Z) for phasing & rf
        // Note                         : ASSUMES PULSE CHANNEL SET

  {
  deleteHams();					// First delete existing Hs
  deleteUprops();				// Delete any step propagators
  deleteGprops();				// Delete any step superprops
  if(!WFsteps)					// If no pulse waveform defined
    {						// we can't set up Hamiltonians
    CPulerror(2, 1);
    CPulfatality(3);
    }
  if(!Iso.length())				// If no pulse channel quit
    {
    CPulerror(4, 1);
    CPulfatality(3);
    }
  Fzed = FZ; 					// Set the Fz operator
  double phi=0, gamB1=0;			// For rf phase and strength 
  complex zx, zy;				// Scratch values
  SetHIndxs( );					// Set Hamiltonian indices	
  Hsteps = new gen_op[WFsteps];			// and their indices
  for(int i=0; i<WFsteps; i++)			// Construct only unique ones
    {						// which are flagged by their
    if(i == Hindex[i])				// self pointing index
      {
      phi       = phase(i)*DEG2RAD; 		//	Step phase (radians)
      gamB1     = WFvals.getRe(i); 		//	Step strength (Hz)
      zx        = gamB1*cos(phi);		//	RF-field x-component
      zy        = gamB1*sin(phi);		//	RF-field y-component
      Hsteps[i] = H + zx*FX + zy*FY; 		// 	Ham. during step i (Hz)
      }
    }
  }

 
void PulComposite::SetLs()

        // Input                CPul    : A composite pulse (this)
        // Output               void    : Effective Hamiltonian superops for
	//				  all composite pulse steps generated
        // Note                         : Will delete any stored superprops
        // Note                         : Assumes H is Ho and that it needs
	//				  AND that Hsteps already present

  {
  Lsteps = new super_op[WFsteps];
  for(int i=0; i<WFsteps; i++)			// Construct only unique ones
    if(i == Hindex[i])				// by insuring H is unique
      Lsteps[i] = Hsuper(Hsteps[i]);		// L = [Heff, ] (1/sec)
  }


// ____________________________________________________________________________
// vi      CLASS COMPOSITE PULSE HILBERT SPACE PROPAGATOR GENERATORS
// ____________________________________________________________________________

/* These functions allow for the filling up two arrays of propagators for each
   step in the waveform (Usteps and Usums).  These well NOT be generated during
   construction of the composite pulse although as they are system dependent.
   Rather, they will be generated when requested if the step Hamiltonians are
   present.								     */
 
void PulComposite::SetUs()
 
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
 
  {
  if(!WFsteps) return;			// Exit if NO steps
  if(!CheckH(10,15)) CPulfatality(16);	// Insure Hamiltonians exist
  if(!Usteps)                           // If propagators don't yet exist
    Usteps = new HSprop[WFsteps];	// allocate some space for them
  if(!Usums)                            // If sum propagators don't yet exist
    Usums = new HSprop[WFsteps];	// allocate some space for them
  if(!Uindex) SetUIndxs( );		// Set propagator indices	
  double tinc;                          // Step time increment
  for(int i=0; i<WFsteps; i++)
    {
    if(i == Uindex[i])			// If unique step propagator
      {					// then we have to calculate it
      tinc = WFtimes.getRe(i);		// Step i length
      Usteps[i]=HSprop(Hsteps[Hindex[i]],tinc);
      }
    Usums[i] = Usteps[Uindex[i]];	// Start prop step i full evolve
    if(i) Usums[i] *= Usums[i-1];       // Add all previous steps evolutions
    }
  }
     

// ____________________________________________________________________________
// vii             CLASS COMPOSITE PULSE STEADY-STATE HANDLING
// ____________________________________________________________________________


void PulComposite::SetSSs()

        // Input                CPul	: A composite pulse (this)
        // Output               void	: Steady-state operators for each
        //				  step of the composite pulse are
	//			   	  calculated

  {
  if(!WFsteps) return;			// Exit if NO steps
  if(!CheckH(10,15)) CPulfatality(16);	// Insure Hamiltonians exist
  if(!R.size()) return;			// No can do if no relaxation set
  if(!SigSSs)                           // If propagators don't yet exist
    SigSSs = new densop[WFsteps];	// allocate some space for them
  super_op Leff;
  for(int i=0; i<WFsteps; i++)
    {
    if(i == Hindex[i])			// If unique Hamiltonian 
      {					// then calculate unique ss
      Leff = complexi*Lsteps[i] + R;	// Full Liouvillian, ith step
      SigSSs[i]=densop(Leff,R,sigmaeq);	// Steady state dens.op. ith step
      }
    }
  }
     

// ____________________________________________________________________________
// viii   CLASS COMPOSITE PULSE LIOUVILLE SPACE PROPAGATOR GENERATORS
// ____________________________________________________________________________

/* These functions allow for the filling up two arrays of Liouville space
   propagators for each step in the waveform (Gsteps and Gsums).  These will
   NOT be generated during construction of the composite pulse although they
   are system dependent.  Rather, they will be generated when requested if the
   step Hamiltonians are present.					     */

void PulComposite::SetGs( )

        // Input                CPul    : A composite pulse (this)
        // Output               void	: Superpropagators for each step
	//				  of the composite pulse are calculated
        // Note				: This assumes that the Hamiltonians
        //                                have already been calculated!
	// Note				: Assumes the relaxation superoperator
	//				  has already been set
        // Note				: G[i] is superpropagator that evolves
        //                                the system from the start of the
	//				  composite pulse to the end of step i
	// Note			        : Doesn't destroy any HS propagators

  {
  if(!WFsteps) return;			// Nothing if no steps specified
  if(!CheckH(10,20)) CPulfatality(21);	// Insure Hamiltonians exist
  if(!Lsteps) SetLs();			// Construct H superops if needed
  if(!Gsteps)Gsteps=new LSprop[WFsteps];// Allocate superprop space if needed
  if(!Gsums) Gsums=new LSprop[WFsteps];	// Allocate sum superprop space
  if(!Uindex) SetUIndxs( );		// Set propagator indices	
  if(!SigSSs) SetSSs();			// Set up steady state density ops
  double tinc;				// Time increment per step 
  int i=0, idx=0;
  densop sigss;
  for(i=0; i<WFsteps; i++)
    {
    idx = Uindex[i];			// Index of propagator
    if(i == idx)			// If unique step propagator
      {					// then we have to calculate it
      tinc = WFtimes.getRe(i);		// Step i length
      Gsteps[i] =			// Prop of step i's evolution
         LSprop(Leff(i),SigSS(i),tinc);
      }
    Gsums[i] = Gsteps[idx];		// Start prop step i full evolve
    if(i) Gsums[i] *= Gsums[i-1];       // Add all previous steps evolutions
    }
  }

// ____________________________________________________________________________
// ix              CLASS COMPOSITE PULSE BASIS HANDLING
// ____________________________________________________________________________

/* These functions are only active on propagators.  They will have no effect
   if there have been no propagators generated.                              */

void PulComposite::SetUBasis(gen_op& Op)

        // Input                CPul    : A composite pulse (this)
        // Output               void	: Hilbert space propagators for all
	//				  steps of the composite pulse are put
	//				  into the current working basis of Op
        // Note				: This assumes that the Hamiltonians
        //                                have already been calculated!

  {
  int i=0;
  if(Usteps)				// If step propagators exist
    { 					// set the basis for each one
    for(i=0; i<WFsteps; i++)
      Usteps[Uindex[i]].SetBasis(Op);
    }
  if(Usums)                            // If sum propagators exist
    {					// set the basis for each one
    for(i=0; i<WFsteps; i++)
      Usums[i].SetBasis(Op);
    }
  }


void PulComposite::SetGBasis(super_op& LOp)

        // Input                CPul    : A composite pulse (this)
        // Output               void	: Liouville space propagators for all
	//				  steps of the composite pulse are put
	//				  into the current working basis of LOp
        // Note				: This assumes that the Hamiltonians
        //                                have already been calculated!

  {
  int i=0;
  if(Gsteps)				// If step propagators exist
    { 					// set the basis for each one
    for(i=0; i<WFsteps; i++)
      Gsteps[Uindex[i]].SetBasis(LOp);
    }
  if(Gsums)                            // If sum propagators exist
    {					// set the basis for each one
    for(i=0; i<WFsteps; i++)
      Gsums[i].SetBasis(LOp);
    }
  }

// ____________________________________________________________________________
// x                  CLASS COMPOSITE CHECKING FUNCTIONS
// ____________________________________________________________________________


int PulComposite::CheckH(int eflag, int eflag2) const

        // Input                CPul    : A composite pulse (this)
	//			eflag	: Error index
        // Output               T/F	: Check to insure the composite pulse
	//				  step Hamiltonians are present

  {
  if(!Hsteps) 				// If Hamiltonians dont exist
    {					// then this is bad.  Issue error
    CPulerror(eflag, 1);
    if(eflag2) CPulerror(eflag2,1);
    return 0;
    }
  return 1;				// Hamiltonians are O.K.
  }


int PulComposite::CheckCH(const spin_sys& sys, const std::string& I, int ef) const

        // Input                CPul    : A composite pulse (this)
	//			sys	: A spin system
	//			I	: An isotope channel
	//			ef	: Error index
        // Output               T/F	: Check to insure channel O.K.

  {
  if(!sys.isotopes(I))		// If system has no isotopes of
    {					// type ch then this is bad
    CPulerror(30, 1);
    return 0;
    }
  return 1;				// Channel is O.K.
  }


int PulComposite::CheckStep(int stp, int eflag, int eflag2) const

        // Input                CPul    : A composite pulse (this)
	//			stp	: Composite pulse step
	//			eflag	: Error index
	//			eflag2	: Second error index
        // Output               T/F	: Check to insure the composite pulse
	//				  step exists

  {
  if(stp>=WFsteps || stp<0)		// Insure step within WF
    {
    CPulerror(eflag, 1);
    if(eflag2) CPulerror(eflag2,1);
    return 0;
    }
  return 1;				// Step index O.K.
  }


void PulComposite::SetNULL()

        // Input                CPul    : A composite pulse (this)
        // Output               void	: Make everything in CPul NULL
	// Note				: Call only during construction
	//				  (or after full destruction!)

  {
  Iso 	    = "";			// No composite pulse channel
  Hsteps    = NULL;			// No composite pulse Hamiltonians
  Hindex    = NULL;			// No composite pulse Hamilton. indices
  Usteps    = NULL;			// No composite pulse propagators
  Usums     = NULL;			// No composite pulse sum propagators
  Uindex    = NULL;			// No (super) propagator indices
  Lsteps    = NULL;			// No Hamiltonian superoperators
  Gsteps    = NULL;			// No composite pulse superpropagators
  Gsums     = NULL;			// No composite pulse sum superprops
  SigSSs    = NULL;			// No steady state opertors
  Hcount    = 0;			// No stored Hamiltonians
  Ucount    = 0;			// No stored propagators
  cutzero   = CP_TCUT;			// Set time precision
  }


LSprop PulComposite::LSIprop() const

        // Input                CPul    : A composite pulse (this)
        // Output               GI	: An identity superpropagator

  {
  int HS = Hsteps[0].dim();		// Hilbert space dimenstion
  return LSprop(HS*HS);			// Return I superprop
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                  CLASS COMPOSITE PULSE CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________


PulComposite::PulComposite() : PulWaveform() { SetNULL(); }

	// Input	none	:
	// Output	CPul	: A NULL PulComposite (this)


PulComposite::PulComposite(const PulWaveform& pulwf,
                                  const spin_system& sys, const std::string& I)
             :PulWaveform(pulwf)

	// Input	pulwf	: Pulse waveform
	//		sys	: A spin system
	//		I	: An isotope channel
	// Output	CPul	: A new composite pulse (this)

  {
  SetNULL();				// Zero everything
  if(!CheckCH(sys,I,30))CPulfatality(1);// Check channel
  Iso = I;				// Set composite pulse channel
  SetHs(sys);		 		// Build step Hamiltonians
  }

 
PulComposite::PulComposite(const PulWaveform& pulwf, gen_op& H,
                           gen_op& FX, gen_op& FY, gen_op& FZ, const std::string& I)
             :PulWaveform(pulwf)

        // Input        pulwf   : Pulse waveform
        //              H       : Active isotropic Hamiltonian
        //              FX      : Active Fx operator
        //              FY      : Active Fy operator
        //              FZ      : Active Fz operator
        //              I	: An isotope channel
        // Output       CPul    : A new composite pulse (this)
        // Note                 : Assumes H,FX,FY,FZ set for isoch

  {
  SetNULL();				// Zero everything
  Iso = I;				// Set composite pulse channel
  SetHs(H, FX, FY, FZ); 		// Build step Hamiltonians
  }


PulComposite::PulComposite(const PulWaveform& pulwf,
                  const spin_system& sys, const super_op& LOp, const std::string& I)
             :PulWaveform(pulwf)

	// Input	pulwf	: Pulse waveform
	//		sys	: A spin system
	//		LOp	: Relaxation/exchange superoperator
	//		isoch	: An isotope channel
	// Output	CPul	: A new composite pulse (this)

  {
  SetNULL();				// Zero everything
  if(!CheckCH(sys,I,30))CPulfatality(1);// Check channel
  Iso       = I;			// Set composite pulse channel
  Fzed      = Fz(sys,I);		// Set Fz operator for channel
  R         = LOp;			// Set relax/exchange superop
  sigmaeq   = densop(sys);		// Get equilibrium density op
  SetHs(sys);		 		// Build step Hamiltonians
  }


PulComposite::PulComposite(const PulComposite& CPul1) : PulWaveform(CPul1)

	// Input	CPul1	: Composite Pulse
	// None		CPul	: Composite Pulse (this), an identical
	//			  copy of CPul1

  {
  SetNULL();				// Zero everything
  Iso       = CPul1.Iso;		// Copy composite pulse channel
  Fzed      = CPul1.Fzed;		// Copy Fz operator for channel
  R         = CPul1.R;			// Copy relax/exchange superoperator 
  sigmaeq   = CPul1.sigmaeq;		// Copy equilib. density operator
  cutzero   = CPul1.cutzero;		// Copy time precision
  copyHams(CPul1);			// Copy any Hamiltonians & indices
  copyUprops(CPul1);			// Copy any Hilbert Propagators
  copyLOps(CPul1);			// Copy any Hamiltonian superops
  copyGprops(CPul1);			// Copy any Liouville Propagators
  copySSs(CPul1);			// Copy any steady-state operators
  }


// ------------------------------ Destruction ---------------------------------

PulComposite::~PulComposite()

	// Input	CPul	: A composite pulse (this)
	// Output	none	: CPul is destructed

  { 
  if(Hsteps) delete [] Hsteps;		// Delete any current Hamiltonians
  if(Hindex) delete [] Hindex;		// Delete any current Ham. indices
  if(Usteps) delete [] Usteps;		// Delete any current propagators
  if(Usums)  delete [] Usums;		// Delete any current sum propagators
  if(Uindex) delete [] Uindex;		// Delete any current prop indices
  if(Lsteps) delete [] Lsteps;		// Delete any Hamiltonian superops
  if(Gsteps) delete [] Gsteps;		// Delete any current superpropagators
  if(Gsums)  delete [] Gsums;		// Delete any current sum superprops
  if(SigSSs) delete [] SigSSs;		// Delete any steady state opertors
  }

// ------------------------------- Assignment ---------------------------------


PulComposite& PulComposite::operator = (const PulComposite& CPul1)

	// Input	CPul1	: Composite Pulse
	// None		CPul	: Composite Pulse (this) copied from CPul1

  {
  PulWaveform::operator=(CPul1);	// Copy pulse waveform
  deleteHams();				// Delete current Hamiltonians
  deleteUprops();			// Delete current propagators
  deleteGprops();			// Delete current superpropagators
  Iso      = CPul1.Iso;			// Copy composite pulse channel
  Fzed     = CPul1.Fzed;		// Copy Fz operator for channel
  R        = CPul1.R;			// Copy relax/exchange superoperator 
  sigmaeq  = CPul1.sigmaeq;		// Copy equilib. density operator
  copyHams(CPul1);			// Copy any Hamiltonians & indices
  copyUprops(CPul1);			// Copy any Hilbert Propagators
  copyGprops(CPul1);			// Copy any Liouville Propagators
  copyLOps(CPul1);			// Copy any Hamiltonian superops
  copySSs(CPul1);			// Copy any steady-state operators
  cutzero  = CPul1.cutzero;		// Copy time precision

  return *this;
  }

// ____________________________________________________________________________
// B  CLASS COMPOSITE PULSE HAMILTONIAN/LIOUVILLIAN/STEADY-STATE FUNCTIONS
// ____________________________________________________________________________

/* These functions allow for access of individual waveform step Hamiltonians,
   Liouvillians and Steady-state density operators.                          */
 

gen_op PulComposite::GetH(int i) const

        // Input                CPul    : A composite pulse (this)
        //                      i	: Step in composite pulse waveform
        // Output               H	: The effective Hamiltonian
	//				  for composite pulse step i

  {
  if(!CheckH(10)) CPulfatality(11);	// Insure Hamiltonians exist
  return Hsteps[(Hindex[i])];		// Return requested Hamiltonian
  }


super_op PulComposite::L0(int i) { return GetL0(i); }
super_op PulComposite::GetL0(int i)

        // Input                CPul    : A composite pulse (this)
        //                      i	: Step in composite pulse waveform
        // Output               L	: The effective Hamiltonian superop
	//				  for composite pulse step i

/* 			L = -i*[Heff, ]     (rad/sec)                        */


  {
  if(!CheckH(10)) CPulfatality(11);	// Insure Hamiltonians exist
  if(!Lsteps) SetLs();			// Construct H superop if needed
  return Lsteps[(Hindex[i])];		// Return requested H superop
  }


super_op PulComposite::Leff(int i) { return complexi*GetL0(i)+R; }
super_op PulComposite::GetLeff(int i) { return complexi*GetL0(i)+R; }

        // Input                CPul    : A composite pulse (this)
        //                      i	: Step in composite pulse waveform
        // Output               L	: The effective Liouvillian superop
	//				  for composite pulse step i

/* 			L = -i*[Heff, ] + R (rad/sec)                        */


densop PulComposite::SigSS(int i)

        // Input                CPul    : A composite pulse (this)
        //                      i	: Step in composite pulse waveform
        // Output               SigmaSS : The steady-state density operator
	//				  for step i

  {
  if(!CheckH(10)) CPulfatality(11);	// Insure Hamiltonians exist
  if(!SigSSs) SetSSs();			// Construct steady state if needed
  return SigSSs[(Hindex[i])];		// Return requested density operator
  }


// ____________________________________________________________________________
// C          CLASS COMPOSITE PULSE HILBERT SPACE PROPAGATOR FUNCTIONS
// ____________________________________________________________________________
 
/* These functions allow for access of individual waveform step Hilbert space
   propagators.  These can either be for a single step evolution or for the
   evolution from pulse start to step end.                                   */
 

HSprop PulComposite::GetU(int i)

        // Input                CPul	: A composite pulse (this)
        //                      i	: Step in composite pulse
        // Output               U	: The propagator for composite pulse
        //				  that evolves from start of step i
        //                                             to the end of step i
        // Note				: An exception is made for default
	//				  then the return will be the
	//				  propagator for the entire pulse
	//				  waveform
	// Note				: If requested and not existing,
	//				  this will trigger ALL propators 
	//				  to be generated if possible

  {
  if(i<0) return GetUsum(i);		// For negative i, return G full WF
  if(!CheckH(10, 17)) CPulfatality(16); // Insure Hs exist, fatal if not
  if(!CheckStep(i,40)) CPulfatality(41);// Insure step within WF
  if(!Usteps) SetUs(); 			// If props dont exist, build em
  return Usteps[Uindex[i]];		// Just return the requested U
  }


HSprop PulComposite::GetU(int i, double td)
 
        // Input                CPul    : A composite pulse (this)
        //                      i       : Pulse step index 
        //                      td      : Evolution time
        // Output               U       : The propagator that evolves a system
        //                                under composite pulse step i for
        //                                time td.

  {
  if(!CheckH(10, 17)) CPulfatality(16); // Insure Hs exist, fatal if not
  if(!CheckStep(i,40)) CPulfatality(41);// Insure step within WF
  return HSprop(Hsteps[Hindex[i]],td); 	// Just return the requested U
  }
 

HSprop PulComposite::GetUsum(int i)

        // Input                CPul	: A composite pulse (this)
        //                      i	: Step in composite pulse
        //                                Default is last propagator -
        //                                i.e. for the entire composite pulse
        // Output               U	: The propagator for composite pulse
        //				  that evolves through i steps of
        //                                the waveform.
	// Note				: If requested and not existing,
	//				  this will trigger ALL propators 
	//				  to be generated if possible
	// Note				: Negative i --> Full Waveform

  {
  if(!CheckH(10, 18)) CPulfatality(16); // Insure Hs exist, fatal if not
  if(!i) return HSprop(Hsteps[0].dim());// No steps return I propagator
  if(!Usums) SetUs(); 			// If props dont exist, build em
  if(i<0) i = WFsteps;			// Set for total cycle prop
  return Usums[i-1];			// Just return the requested U
  }


HSprop PulComposite::GetUmult(int N)
 
        // Input                CPul	: A composite pulse (this)
        //                      N       : Number of composite pulses
        // Output               U       : The propagator for N composite
        //                                pulses applied in succession
	// Note				: U is returned in its eigenbasis
 
  {
  if(!CheckH(10, 18)) CPulfatality(16); 	// Insure Hs exist, else fatal
  if(N<=0) return HSprop(Hsteps[0].dim());	// Identity if no steps
  if(!Usums) SetUs();				// If no props, build em
  return Usums[WFsteps-1].Pow(N);
  }


HSprop PulComposite::GetU(double td)

        // Input                CPul	: A composite pulse (this)
        //                      td	: Evolution time
        // Output               U	: The propagator that evolves a system
        //				  for time td under the composite pulse

  {

//  Determine Propagators For Waveform(s), Steps(s), Last Partial Step

  int nWF = fullWFs(td);		// Number of complete waveforms 
  HSprop UWFs = GetUmult(nWF);		// Propagator over these waveforms
  td -= UWFs.length();;			// Time we have left to evolve
  int nS = fullsteps(td);		// Number of complete steps
  HSprop UStps = GetUsum(nS);		// Propagator over these steps
  td -=  UStps.length();		// Time we have left to evolve
  int j = Hindex[nS];			// Last (partial) step index

//		 Build HS Propagator For Steps, REVERSE Order

  HSprop U = HSprop(Hsteps[j],td);	// Prop of step i's evolution
  U *= UStps; 				// Evolve over full waveform steps
  U *= UWFs;				// Evolve over full waveforms
  return U;
  }
    

// ____________________________________________________________________________
// D       CLASS COMPOSITE PULSE SUPEROPERATOR PROPAGATOR FUNCTIONS
// ____________________________________________________________________________

/* These functions allow for access of individual waveform step Liouville space
   propagators.  These will be for the evolution pulse start to step end.    */

LSprop PulComposite::GetG(int i)

        // Input                CPul	: A composite pulse (this)
        //                      i	: Step in composite pulse
        //                                Default is last superpropagator -
        //                                i.e. for the entire composite pulse
        // Output               G	: Superpropagator for composite pulse
        //                                that evolves from the start of
        //                                the waveform to the end of step i
	// Note				: Assumes relaxation operator set

  {
  if(i<0) return GetGsum(i);		// Return full WF prop if i negative
  if(!CheckH(10, 22))  CPulfatality(21);// Insure Hs exist, fatal if not
  if(!CheckStep(i,40)) CPulfatality(41);// Insure step within WF
  if(!Gsteps) SetGs();			// Build propagators if needed
  return Gsteps[Uindex[i]];		// Just return the requested G
  }

 
LSprop PulComposite::GetG(int i, double td)
 
        // Input                CPul    : A composite pulse (this)
        //                      i       : Pulse step index
        //                      td      : Evolution time
        // Output               G       : The propagator that evolves a system
        //                                under composite pulse step i for
        //                                time td.
 
  {  
  if(!CheckH(10, 22))  CPulfatality(21);// Insure Hs exist, fatal if not
  if(!CheckStep(i,40)) CPulfatality(41);// Insure step within WF
  if(!Lsteps) SetLs();			// Construct Hamiltonian superops
  if(!SigSSs) SetSSs();			// Construct steady state ops
  return LSprop(Leff(i),SigSS(i),td);	// Just return requested G
  }


LSprop PulComposite::GetGsum(int i)
 
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
 
  {
  if(!CheckH(10, 18)) CPulfatality(16); // Insure Hs exist, fatal if not
  if(!i)return LSIprop();		// No step, return I superprop
  if(!Gsums) SetGs();			// If props dont exist, build em
  if(i<0) i = WFsteps;			// Set for total cycle prop
  return Gsums[i-1];			// Just return the requested G
  }


LSprop PulComposite::GetGmult(int N)
 
        // Input                CPul	: A composite pulse (this)
        //                      N       : Number of composite pulses
        // Output               G       : Superpropagator for N composite
        //                                pulses applied in succession
	// Note				: Non-positive N returns an
	//				  identity superpropagator 	
 
  {
  if(!CheckH(10, 18)) CPulfatality(16); // Insure Hs exist, fatal if not
  if(N<=0) return LSIprop();		// No steps, return I superprop
  if(!Gsums) SetGs();			// If no props, build em
  LSprop GWF = Gsums[WFsteps-1];	// Prop for 1 waveform
  LSprop GNWF = GWF;			// Set for 1 waveform
  for(int i=1; i<N; i++) GNWF *= GWF;	// Add in next waveforms
  return GNWF;
  }


LSprop PulComposite::GetG(double td)
 
        // Input                CPul    : A composite pulse (this)
        //                      td      : Evolution time
        // Output               G       : Superpropagator that evolves a system
        //                                for time td under the composite pulse

  {

//  Determine Propagators For Waveform(s), Steps(s), Last Partial Step

  int nWF = fullWFs(td);                // Number of complete waveforms
  LSprop GWFs = GetGmult(nWF);          // Propagator over these waveforms
  td -= nWF*WFtp;                       // Time we have left to evolve
  int i = fullsteps(td);		// Number of complete steps
  LSprop GStps = GetGsum(i);		// Propagator over these steps
  td -=  sumlength(i);			// Time we have left to evolve
  if(fabs(td) < cutzero) td = 0;	// Zero time if below cutoff
 
//               Build LS Propagator For Steps, REVERSE Order
 
  LSprop G=LSprop(Leff(i),SigSS(i),td);	// Final partial step superprop
  G *= GStps;                           // Evolve over full waveform steps
  G *= GWFs;                            // Evolve over full waveforms
  return G;
  }



// ____________________________________________________________________________
// E               CLASS COMPOSITE PULSE ACCESS FUNCTIONS
// ____________________________________________________________________________

// --------------------- Functions Over Full Composite ------------------------

/*
int          PulComposite::steps()    const { return WFsteps; }  INHERITED
std::string       PulComposite::name()     const { return name(); }   INHERITED
double       PulComposite::length()   const { return length(); } INHERITED
row_vector   PulComposite::values()   const { return values(); } INHERITED   */

std::string       PulComposite::channel()  const { return Iso; }

	// Input	CPul	: A composite pulse (this)
	// Output	steps	: CPul steps
	// 		name	: CPul name
	//		length  : CPul length (sec)
	// 		channel	: CPul channel
        //              values  : Array of { gamB1, phi } values


// ----------------- Functions For Specific Composite Step --------------------

/*
double  PulComposite::strength(int i) const { return strength(i); } INHERITED
double  PulComposite::phase(int i)    const { return phase(i); }    INHERITED
complex PulComposite::value(int i)    const { return value(i); }    INHERITED
double  PulComposite::length(int i)   const { return length(i); }   INHERITED */

	// Input	CPul	: A composite pulse (this)
	// Output	strength: Step rf-field strength (Hz)
	//		phase	: Step phase value (degrees)
        //              value   : Step { gamB1, phi } values
	//              length  : Step length (sec)


// --------------------- Other Pulse Composite Access -------------------------

gen_op       PulComposite::FZ()        const { return Fzed; }
super_op     PulComposite::ROp()       const { return R; }
densop	     PulComposite::SigEq()     const { return sigmaeq; }
double       PulComposite::Precision() const { return cutzero; } 

        // Input        CPul    : A composite pulse (this)
        // Output       FZ	: Active Fz operator
	//		waveform: Pulse waveform
        // Output       R	: Active Relaxation/exchange superoperator
	//		sigmaeq : Equilibrium density operator
	//		cutzero : Precision with which time is known


// ____________________________________________________________________________
// F               CLASS COMPOSITE PULSE AUXILIARY FUNCTIONS
// ____________________________________________________________________________

/*
double PulWaveform::maxlength( )          const;		INHERITED
double PulWaveform::minlength(double)     const;		INHERITED
int    PulWaveform::gamB1const()          const;		INHERITED
int    PulWaveform::phaseconst()          const;		INHERITED
int    PulWaveform::timeconst()	          const;		INHERITED
double PulWaveform::steps(double td)      const; 		INHERITED
int    PulWaveform::fullsteps(double td)  const;		INHERITED
double PulWaveform::WFs(double td)        const;		INHERITED
int    PulWaveform::fullWFs(double td)    const;		INHERITED
double PulWaveform::sumlength(int i)      const;		INHERITED    */


void PulComposite::scalegB1(double sf)

        // Input        CPul    : A composite pulse (this)
        //              sf      : A scaling factor
        // Output       void    : All step field strengths in PWF
        //                        are multiplied by sf.  The exception are
        //                        steps of zero length (ideal pulses)

  {
  deleteUprops();			// These won't be applicable
  deleteGprops();			// These won't be applicable
  deleteHams();				// These won't apply either
  PulWaveform::scalegB1(sf);		// Adjust waveform field strengths
// sosix - this would work if Ho, FX, FY, & FZ is O..K.
  }

 
 
void PulComposite::setRelax(const spin_system& sys, const super_op& LOp)
  
        // Input        CPul    : A composite pulse (this)
        //              sys     : A spin system
        //              LOp     : Relaxation/exchange superoperator
        // Output       void    : Relaxation specifications are put int
        //                        CPul

  {
  R         = LOp;			// Copy relax/exchange superoperator 
  sigmaeq   = densop(sys);		// Set equilibrium density op
  }


// ____________________________________________________________________________
// G               CLASS COMPOSITE PULSE EVOLUTION FUNCTIONS
// ____________________________________________________________________________

// ----------------------- Acquisiton Tracking Functions ----------------------

  
        // Input        CPul    : A composite pulse (this)
	//		typ     : FID header type
	//		rlx     : Relaxation flag (0=no, !0=yes)
        // Output       void    : Outputs header for FID tracking

void PulComposite::FIDheader(int typ, int rlx) const
  {
  std::string WF(" Waveform");
  std::string SYN(" Synchronized ");
  std::string ASYN(" Asynchronous ");
  std::string AT("Acquistion Tracking");
  std::string RA(", Relaxation Active");
  std::string FPT[3] = { "  FID  ", " Point ", "-------" };
//  std::string FCY[3] = { "  Full   ", " Cycles  ", "---------" };
  std::string FWF[3] = { "  Full   ", "Waveforms", "---------" };
  std::string FST[3] = { "Full ", "Steps", "-----" };
  std::string PST[3] = { "   Partial  ", "    Step    ", "------------" };
  std::string ETS[3] = { " Evolution  ", "    Time    ", "------------" };
  std::string PVS[3] = { "          Point Values",
                    "    Real    Imaginary   Norm",
                    "-------------------------------" };

//			Set Spacing At Start Of Each Line

  std::string LS;					// Start of each line
  switch(typ)
    {
    case 0:
    case 1:
    case 3:
    case 4:
    case 5:
    default:
      LS = "\n\t"; break;
    case 2:
    case 6:
    case 7:
      LS = "\n  "; break;
    }

  int Hsp=0;				// Header spacing length 
  int asyn=0;				// Asynchronous acquire flag
  switch(typ)
    {
    case 0: Hsp=13;         break;
    case 1: Hsp=18;         break;
    case 2: Hsp=18;         break;
    case 3: Hsp=28; asyn=1; break;
    case 4: Hsp=4;  asyn=1; break;
    case 5: Hsp=9;          break;
    case 6: Hsp=6;          break;
    case 7: Hsp=20; asyn=1; break;
    }
  std::string spacer = std::string(Hsp-WFname.length()/2, ' ');
  std::string Header = "\n" + LS + spacer + WFname;
  switch(typ)
    {
    case 0: Header += WF; break;
    case 1: Header += WF + " Step"; break;
    case 2: Header += WF + " Fraction"; break;
    case 3: Header += WFname; break;
    case 4: Header += WFname; break;
    case 5: Header += WF; break;
    case 6: Header += WF + " Fraction" + SYN; break;
    }
  if(!asyn) Header += SYN;		// Either add in "Synchronized"
  else      Header += ASYN;		// Or     add in "Asynchronous"
  Header += AT;				// Add in "Acquisiton Tracking"
  if(rlx) Header += RA;			// Add in "Relaxation Active"
  Header += "\n";			// Add in a line return

//			Set Spacing Between Each Column

  std::string SP = std::string(3, ' ');

//			Form 1st 2 And Last 2 Columns

  std::string Colsi2[3], Colsf2[3];
  int i, j, nr=3;
  int iFST=0, iPST=0;
  for(i=0; i<3; i++)
    {
    Colsi2[i] = LS + FPT[i] + SP + FWF[i] + SP;
    Colsf2[i] = ETS[i] + SP + PVS[i];
    }

  std::cout << Header;			// Output the top Header
  for(j=0; j<nr; j++)
    {
    std::cout << Colsi2[j];			// Print space & 2 columns
    if(iFST) std::cout << FST[i] << SP;
    if(iPST) std::cout << PST[i] << SP;
    std::cout << Colsf2[i];
    }

  switch(typ)
    {
    case 0:
    default:
      std::cout << Header;
      iFST = 0;
      for(j=0; j<nr; j++)
        {
        std::cout << Colsi2[j];			// Print space & 2 columns
        std::cout << Colsf2[0];
        }
      break;
    case 1:
      std::cout << Header;
      std::cout << Colsi2[0] << FST[0] << SP << Colsf2[0];
      std::cout << Colsi2[1] << FST[1] << SP << Colsf2[1];
      std::cout << Colsi2[2] << FST[2] << SP << Colsf2[2];
      break;
    case 2:
      std::cout << Header;
      std::cout << Colsi2[0] << FST[0] << SP << PST[0] << SP << Colsf2[0];
      std::cout << Colsi2[1] << FST[1] << SP << PST[1] << SP << Colsf2[1];
      std::cout << Colsi2[2] << FST[2] << SP << PST[2] << SP << Colsf2[2];
      break;
    case 3:
      std::cout << Header;
      std::cout << Colsi2[0] << FST[0] << SP << "Partial" << SP << PST[0] << SP << Colsf2[0];
      std::cout << Colsi2[1] << FST[1] << SP << " Step  " << SP << PST[1] << SP << Colsf2[1];
      std::cout << Colsi2[2] << FST[2] << SP << "-------" << SP << PST[2] << SP << Colsf2[2];
      break;
    case 4:
      std::cout << Colsi2[0] << Colsf2[0];
      std::cout << Colsi2[1] << Colsf2[1];
      std::cout << LS << Colsi2[2] << Colsf2[2];
      break;
    case 5:
      std::cout << Colsi2[0] << FST[0] << SP << Colsf2[0];
      std::cout << Colsi2[1] << FST[0] << SP << Colsf2[1];
      std::cout << Colsi2[2] << FST[0] << SP << Colsf2[2];
      break;
    case 6:
      std::cout << Colsi2[0] << FST[0] << SP << PST[0] << Colsf2[0];
      std::cout << Colsi2[1] << FST[1] << SP << PST[1] << Colsf2[1];
      std::cout << Colsi2[2] << FST[2] << SP << PST[2] << Colsf2[2];
      break;
    case 7:
      std::cout << Colsi2[0] << FST[0] << SP << "Partial" << SP << PST[0] << SP << Colsf2[0];
      std::cout << Colsi2[1] << FST[1] << SP << " Step  " << SP << PST[1] << SP << Colsf2[1];
      std::cout << Colsi2[2] << FST[2] << SP << "-------" << SP << PST[2] << SP << Colsf2[2];
    }
  std::cout.flush();
  }


void PulComposite::FIDpoint(int typ, int pt, int iWFs, int iSTs) const
  
        // Input        CPul    : A composite pulse (this)
	//		typ     : FID header type
	//		pt      : Point index
	//		iWFs    : Number of waveforms
	//		iSTs	: Number of steps
        // Output       void    : Outputs header for FID tracking

  {
  std::string spacer;
  switch(typ)
    {
    case 0:
    default:
      std::cout << "\n\t" << Gdec(pt+1,5) << ".    ";
      std::cout << Gdec(iWFs,6) << "    ";
      break;
    case 1:
      std::cout << "\n\t" << Gdec(pt+1,5) << ".    ";
      std::cout << Gdec(iWFs,6) << "    ";
      std::cout << Gdec(iSTs,6) << "    ";
      break;
    case 2:
      std::cout << "\n  " << Gdec(pt+1,5) << ".    ";
      std::cout << Gdec(iWFs,6) << "    ";
      std::cout << Gdec(iSTs,6) << "    ";
      break;
    case 3:
      std::cout << "\n  " << Gdec(pt+1,5) << ".   ";
      std::cout << Gdec(iWFs,6) << "   ";
      std::cout << Gdec(iSTs,6) << " ";
      break;
    case 4:
      std::cout << "\n\t" << Gdec(pt+1,5) << ".   ";
      std::cout << Gdec(iWFs,6) << "      ";
      break;
    case 5:
      std::cout << "\n\t" << Gdec(pt+1,5) << ".   ";
      std::cout << Gdec(iWFs,6) << "    ";
      std::cout << Gdec(iSTs,6) << "    ";
      break;
    case 6:
      std::cout << "\n  " << Gdec(pt+1,5) << ".   ";
      std::cout << Gdec(iWFs,6) << "    ";
      std::cout << Gdec(iSTs,6) << "    ";
      break;
    case 7:
      std::cout << "\n  " << Gdec(pt+1,5) << ".   ";
      std::cout << Gdec(iWFs,6) << "   ";
      std::cout << Gdec(iSTs,6) << " ";
      break;
    }
  std::cout.flush();
  }


void PulComposite::FIDvalue(int typ, double td, const complex& z) const
  
        // Input        CPul    : A composite pulse (this)
	//		typ     : FID header type
	//		td	: Delay length
	//		z	: Data point
        // Output       void    : Outputs point info for FID tracking

  {
  std::string spacer;
  switch(typ)
    {
    case 0:
    case 1:
    default:
      printTime(std::cout, td);
      std::cout << "    " << Gform("%8.3f", Re(z))
           << "  "   << Gform("%8.3f", Im(z))
           << "  "   << Gform("%8.3f", norm(z));
    break;
    }
  std::cout.flush();
  }

// sosix

// ------------------------ Acquisiton Helper Functions -----------------------

 
void PulComposite::FIDtell(double SW) const
 
        // Input        CPul    : A composite pulse (this)
        //              SW      : Desired spectral width
 
  {
  if(!maxgamB1( ))
    {
    std::cout << "\n\tThe Decoupler RF-Stength Is Zero!"
         << "\n\tDecoupling Will Play No Part In FID Synchronization...";
    return;
    }
  double td = 1/SW;			// Dwell time for desired SW
  double SWWF=0, SWWFh=0, SWWFl=0;

//		   Attempt To Get A SW For Waveform Synchronization
//		      (Not Possible If Dwell Time << 1 Waveform)

  int nWFs = fullWFs(td);			// Full waveforms within td
  if(nWFs == 1)					// If multiple waveforms are
    {
    SWWF = 1.0/(nWFs * WFtp);			// Synch to full waveform
    std::cout << "\n\tWaveform Sync @ " << nWFs << " Waveforms: " << SWWF << " Hz";
    }
  else if(nWFs>1)
    {
    SWWF = 1.0/(nWFs * WFtp);			// Synch to full waveform
    SWWFl = 1.0/((nWFs+1) * WFtp);		// Synch to full waveform
    SWWFh = 1.0/((nWFs-1) * WFtp);		// Synch to full waveform
    std::cout << "\n\tWaveform Sync @ " << nWFs << " Waveforms: " << SWWF << " Hz";
    std::cout << "\n\tWaveform Sync @ " << nWFs+1 << " Waveforms: " << SWWFl << " Hz";
    std::cout << "\n\tWaveform Sync @ " << nWFs-1 << " Waveforms: " << SWWFh << " Hz";
    }
  else
    std::cout << "\n\tWaveform Synch is NOT Possible";

//		      Attempt To Get A SW For Step Synchronization

  double SWST=0, SWSTh=0, SWSTl=0;
  if(timeconst())				// If constant step length we
    { 						// maybe can synch on steps
    int nSTPs = fullsteps(td);			// Full steps within td
    if(nSTPs==1)
      {
      SWST = 1.0/(nSTPs * length(0));		// Synch to full steps
      std::cout << "\n\tWaveform Step Sync @ " << nSTPs << " Steps: " << SWST << " Hz";
      }
    else if(nSTPs>1)
      {
      SWST = 1.0/(nSTPs * length(0));		// Synch to full waveform
      SWSTl = 1.0/((nSTPs+1) * length(0));	// Synch to full waveform
      SWSTh = 1.0/((nSTPs-1) * length(0));	// Synch to full waveform
      std::cout << "\n\tWaveform Step Sync @ " << nSTPs << " Steps: " << SWWF << " Hz";
      std::cout << "\n\tWaveform Step Sync @ " << nSTPs+1 << " Steps: " << SWWFl << " Hz";
      std::cout << "\n\tWaveform Step Sync @ " << nSTPs-1 << " Steps: " << SWWFh << " Hz";
      }
    else
      std::cout << "\n\tWaveform Step Synch is NOT Possible";
    }
  else
    std::cout << "\n\tWaveform Step Synch is NOT Possible";
  }

 
double PulComposite::FIDsync(double& SW, int warn) const
 
        // Input        CPul    : A composite pulse (this)
        //              SW      : Desired spectral width
	//		warn	: Flag for warning output
        // Output       SWsync  : Spectral width synchronized with the
        //                        composite pulse length (at best) or
        //                        with the pulse step length (2nd best)
        // Note                 : If no synchronization possible the input
        //                        value is just returned.
 
  {
  double td = 1/SW;                             // Dwell time for desired SW

//		   Attempt To Get A SW For Waveform Synchronization
//		      (Not Possible If Dwell Time << 1 Waveform)

  int nWFs = fullWFs(td);			// Full waveforms within td
  if(nWFs > 0)					// If multiple waveforms are
    {						// spanned by td, shorten td
    if(warn) CPulerror(62);			// a bit (increase SW) to synch
    return 1.0/(nWFs * WFtp);			// Synch to full waveform
    }

//		      Attempt To Get A SW For Step Synchronization

  if(timeconst())				// If constant step length we
    { 						// maybe can synch on steps
    if(warn) CPulerror(63);
    int nSTPs = fullsteps(td);			// Full steps within td
    if(nSTPs> 0) return 1.0/(nSTPs * length(0));// Synch to full steps
    }
  if(warn) CPulerror(61);
  return SW;
  }

 
int PulComposite::FIDtest(double td, int& nWFs, int& nSTs, double& tr) const
 
        // Input        CPul    : A composite pulse (this)
        //              td      : Desired dwell time (sec)
	//		nWFs	: Number of full waveforms
	//		nSTs	: Number of full step
	//		tr	: Remainder time (sec)
        // Output       Fsync	: Synchronization flag
	//			        +: Synchronized to Waveform
	//				0: Not synchronized
	//				-: Synchronized to Step
        // Note                 : Values of nWFs, nSTs, and tr are altered
 
  {
  nWFs = 0;					// Set for no waveforms
  nSTs = 0;					// Set for no steps
  tr = 0;					// Set for no residual time
  if(fabs(td) < cutzero) return 0;		// If no time, return nothing
  double te = td; 				// Evolve time desired
  int iWFs = fullWFs(te);				// Full waveforms within tr
  te  -= iWFs*WFtp;				// Time left to evolve
  if(fabs(te) < cutzero) te = 0;		// Zero time if below precision
  int iSTs = fullsteps(te);				// Additonal full steps
  te  -= sumlength(iSTs);			// Time left to evolve
  if(fabs(te) < cutzero) te = 0;		// Zero time if below precision
  tr = te;
  nWFs = iWFs;
  nSTs = iSTs;
  if(!te)
    {
    if(!iSTs) return 1;				//	Here if waveform synch.
    else      return -1;			//	Here if step synchronized
    }
  else return 0;				//	Here if unsynchronized
  }

// ----------------- Acquisition Without Relaxaton & Exchange -----------------

row_vector PulComposite::FIDsynchWF(int npts, int nWFs,
                                          gen_op &D, gen_op& sigmap, int track)

	// Input	CPul	: A composite pulse (this)
        //              npts	: Number of FID points
        //              nWFs	: Waveforms between FID points
        //              D 	: Detection operator
        //              sigmap	: Prepared density operator
	//		track	: Print output to track calculation
        // Output       data	: A row vector containing npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
	//			  composite pulse CPul.

  {
  row_vector data(npts, complex0);              // Vector for FID points
  int HS = sigmap.dim();                        // Working Hilbert space
  HSprop Ut(HS);                                // Needed propagator
  gen_op sigma;                                 // Needed operators
  D.Op_base(sigmap);                            // Set D basis to sigmap's
  if(track) FIDheader(0, 0); 			// Output header if tracking
 
  HSprop UWFs = Ut;				// Prop. for 0 waveforms
  for(int j=0; j<nWFs; j++) UWFs *= GetUsum();	// Set prop for nWFs waveforms
  int iWFs = 0;                                 // Number of waveforms
  for(int i=0; i<npts; i++)                     // Loop over FID points
    {                                           // Output point info
    sigma = Ut.evolve(sigmap);                  // This does the evolution
    data.put(trace(D,sigma),i);                 // Store FID point
    if(track)                                   // Output FID point if tracking
      {
      FIDpoint(0, i, iWFs, 0);
      FIDvalue(0, Ut.length(), data.get(i));
      }
    if(i != npts-1)
      {
      iWFs += nWFs;                             // Update the waveform count
      Ut *= UWFs;				// Update propagator
      }
    }
  return data;
  }


row_vector PulComposite::FIDsynchST(int npts, int nSTs,
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
 
  {
  if(!timeconst())                              // There can be no step
    {                                           // synchronized acquisition
    CPulerror(60,1);				// if step lengths vary!
    CPulfatality(61);
    }
  row_vector data(npts, complex0);              // Vector for FID points
  int HS = sigmap.dim();                        // Working Hilbert space
  HSprop Ut(HS); 				// Needed propagator
  gen_op sigma;                                 // Needed operators
  D.Op_base(sigmap);                            // Set D basis to sigmap's
  double tstep = length(0);			// Single WF step length
 
  int nWFs = fullWFs(double(nSTs)*tstep);       // Waveforms in 1 increment
  nSTs -= nWFs*WFsteps;                         // Steps left in 1 increment
  if(!nSTs)                                     // If no steps, then use the
    return FIDsynchWF(npts,nWFs,D,sigmap,track);// WF synch FID function
  if(track) FIDheader(1, 0); 			// Output header if tracking
 
  HSprop UWFtot = Ut;                           // Prop. for total waveforms
  HSprop UWFs = Ut;                             // Prop. for WF(s) per pt.
  for(int j=0; j<nWFs; j++) UWFs *= GetUsum();	// Set prop for nWFs wavefomrs
  int iWFs = 0;                                 // Number of waveforms
  int iSTs = 0;                                 // Number of steps
  for(int i=0; i<npts; i++)                     // Loop over FID points
    {                                           // Output point info
    sigma = Ut.evolve(sigmap);                  // This does the evolution
    data.put(trace(D,sigma),i);                 // Store FID point
    if(track)                                   // Output FID point if tracking
      {
      FIDpoint(1, i, iWFs, iSTs);
      FIDvalue(1, Ut.length(), data.get(i));
      }
   if(i != npts-1)
      {
      iWFs   += nWFs;                           // Update waveform count
      iSTs   += nSTs;                           // Update step count
      UWFtot *= UWFs;                           // Update waveform prop
      if(iSTs >= WFsteps)                       // See if steps than full WF
        {                                       // then just use full WF
        iWFs++;                                 //      Add 1 WF to count
        UWFtot *= GetUsum();			//	Add 1 WF to prop 
        iSTs -= WFsteps;                        //      Sub 1 WF from step cnt
        }
      if(iSTs)
        {
        Ut = GetUsum(iSTs);			// Set prop for steps
        Ut *= UWFtot;				// Add all WFs to prop
        }
      else Ut = UWFtot;				// Set prop for all WFs
      }
    }
  return data;
  }


row_vector PulComposite::FIDsynchFR(int npts, int nFRs,
                                          gen_op &D, gen_op& sigmap, int track)
 
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

  {
  row_vector data(npts, complex0);              // Vector for FID points
  int HS = sigmap.dim();                        // Working Hilbert space
  HSprop Ut(HS);				// Needed propagator
  gen_op sigma;                                 // Needed operators
  D.Op_base(sigmap);                            // Set D basis to sigmap's
  if(nFRs < 1)					// Can't allow less than 1
    {                                           // waveform fraction
    CPulerror(65,1);
    CPulfatality(61);
    }
  if(nFRs == 1)					// If only 1 fraction, use
    return FIDsynchWF(npts,1,D,sigmap,track);	// waveform synch acquire
  double td = WFtp/nFRs;			// Dwell time to use
  int nSTs;					// Number of WF steps
  if(timeconst())				// Check if step synch is
    {						// possible too
    nSTs = fullsteps(td);			//   These steps in td 
    if(fabs(nSTs*length(0)-td) < cutzero) 	//   If steps span td
    return FIDsynchST(npts,nSTs,D,sigmap,track);//   use step synch acquire
    }
  if(track) FIDheader(2, 0); 			// Output header if tracking
 
  HSprop *FrUs;
  FrUs = new HSprop[nFRs];				// Space for needed props
  HSprop Usummed = Ut;
  double tr=0, te=td;				// Evolution time
  for(int m=0; m<nFRs; m++)			// Fill up props with evolve
    {						// over waveform fractions
    FrUs[m] = GetU(te);
    te += td;
    }

  int iWFs = 0;                                 // Number of waveforms
  int iSTs = 0;                                 // Number of steps
  te = 0;					// Zero FID evolve time
  int i,j;
  for(i=0; i<npts; i++)				// Loop over FID points
    {
    j = i%nFRs;					// Fraction index
    sigma = Ut.evolve(sigmap);			// This does the evolution
    data.put(trace(D,sigma),i);			// Store FID point
    if(track)					// Output FID point if tracking
      {
      FIDpoint(2, i, iWFs, iSTs);
      if(fabs(tr) > cutzero) printTime(std::cout, tr);
      else                   printTime(std::cout, 0);
      std::cout << "  ";
      FIDvalue(2, Ut.length(), data.get(i));
      std::cout.flush();
      }
    if(i != npts-1)				// Note iWFs, iSTs
      {
      if(!j && i) Usummed *= FrUs[nFRs-1];	// Go to next waveform
      Ut = FrUs[j];				// Update evolution operator
      Ut *= Usummed;				// Add previous WFs
      if(track)
        {
        te += td;				// Update FID evolution time
        iWFs = fullWFs(te);			// Update waveform count
        tr = te-iWFs*WFtp;			// Length remaining
        iSTs = fullsteps(te-iWFs*WFtp);		// Update step count
        tr -= sumlength(iSTs);			// Length remaining
        }
      }
    }
  delete [] FrUs;
  return data;
  }


row_vector PulComposite::FID(int npts, double td,
                                    gen_op &D, gen_op& sigmap, int track)

	// Input	CPul	: A composite pulse (this)
        //              npts	: Number of FID points
        //              td	: Dwell time between FID points
        //              D 	: Detection operator
        //              sigmap	: Prepared density operator
	//		track	: Print output to track calculation
        // Output       data	: A row vector containing npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
	//			  composite pulse CPul.
	// Note			: If the composite pulse length is zero
	//			  then the FID will be a normal acquisition


  {
  row_vector data(npts);                        // Vector for FID points
  SetUBasis(sigmap);				// Put props in sigmap basis
  HSprop Ut, UI(sigmap.dim());			// Needed propagators
  gen_op sigma;                                 // Needed operators
  double tr, te = 0;                            // FID evolution time
  int nWFs, nSTs;				// Waveforms & steps used
  int synch = FIDtest(td,nWFs,nSTs,tr);		// Check FID synchronicity
  if(synch > 0)					// Here if Waveform synch.
    return FIDsynchWF(npts,nWFs,D,sigmap,track);
  else if(synch < 0)				// Here if step synch.
    {
    int nSTtot = nWFs*WFsteps + nSTs;
    return FIDsynchST(npts,nSTtot,D,sigmap,track);
    }
  if(track) FIDheader(3, 0); 			// Output header if tracking

  for(int i=0; i<npts; i++)
    {
//		First Determine Which Propagation Steps Are Required
//		   (For Evolution From FID Start to Current Point)

    tr   = te;					// Evolve time this point
    if(fabs(tr) < cutzero) tr = 0;		// Zero time if below precision
    nWFs = fullWFs(tr);				// Full waveforms within tr
    tr  -= nWFs*WFtp;				// Time left to evolve
    if(fabs(tr) < cutzero) tr = 0;		// Zero time if below precision
    nSTs = fullsteps(tr);			// Additonal full steps
    tr  -= sumlength(nSTs);			// Time left to evolve
    if(fabs(tr) < cutzero) tr = 0;		// Zero time if below precision
    SetUBasis(sigmap);				// Put props in sigmap basis
 
    if(track)
      {
      FIDpoint(3, i, nWFs, nSTs);
      if(tr)
        {
        std::cout << Gdec(nSTs,6) << "     ";
        printTime(std::cout, tr);
        std::cout << "  ";
        }
      else
        std::cout << "                          ";
      std::cout.flush();
      }
    if(tr < 0)                                  // Insure evolution time
      {                                         // is non-negative
      CPulerror(50,1);
      CPulfatality(51);
      }

//		 Next We Must Build A Propagator In REVERSE Order

             Ut = UI;				// Begin with identity prop
             Ut.SetBasis(sigmap);		// Put prop in sigmap basis
    if(tr)   Ut=HSprop(Hsteps[Hindex[nSTs]],tr);// Evolve for partial NSTs step
             Ut.SetBasis(sigmap);		// Put prop in sigmap basis
    if(nSTs) Ut*=GetUsum(nSTs);			// Next evolve nSTs steps
    if(nWFs) Ut*=GetUmult(nWFs);		// Next evolve nWFs waveforms

//            Finally We Evolve Input Sigma To FID Point And Sample

    sigma = Ut.evolve(sigmap);			// This does the evolution
    data.put(trace(D,sigma),i);                 // Store FID point
    if(track) FIDvalue(3,Ut.length(),data.get(i));
    te += td;					// Update FID evolution time
    }
  return data;
  }  


// ------------------- Acquisition With Relaxaton & Exchange ------------------
 

row_vector PulComposite::FIDRsynchWF(int npts, int nWFs,
                                          gen_op &D, gen_op& sigmap, int track)

	// Input	CPul	: A composite pulse (this)
        //              npts	: Number of FID points
	//		nWFs 	: Number of waveforms between points
        //              D 	: Detection operator
        //              sigmap	: Prepared density operator
	//		track	: Flag for tracking the computation
        // Output       data	: A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the
	//		  	  composite pulse CPul.
	// Note			: Assumes Relaxation Active!

  {
  double FIDcut = 1.e-6;			// FID cutoff intensity
  int FIDptcnt = 0;				// FID zero point count
  int FIDnilcnt = 4;				// Allowed FID zero points
  row_vector data(npts, complex0);		// Vector for FID points
  int HS = sigmap.dim();			// Working Hilbert space
  LSprop Gt(HS*HS);				// Needed superprop.
  gen_op sigma;					// Needed operators
  D.Op_base(sigmap);				// Set D basis to sigmap's
  if(track) FIDheader(4, 1); 			// Output header if tracking

  LSprop GWFs = GetGsum();			// Prop. for 1 waveform
  for(int j=1; j<nWFs; j++)			// Set prop for nWFs waveforms
    GWFs *= GetGsum();
  double te = 0.0;				// FID evolution time
  int iWFs = 0;					// Number of waveforms
  for(int i=0; i<npts; i++)			// Loop over FID points
    {						// Output point info 
    sigma = Gt.evolve(sigmap);			// This does the evolution
    data.put(trace(D,sigma),i);                 // Store FID point
    if(track)					// Output FID point if tracking
      {
      FIDpoint(4, i, iWFs, 0);
      FIDvalue(4, Gt.length(), data.get(i));
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
      iWFs += nWFs;				// Update waveform count
      te += WFtp;				// Update FID evolution time
      Gt *= GWFs;				// Update evolution operator
      }
    }
  return data;
  }


row_vector PulComposite::FIDRsynchST(int npts, int nSTs, 
                                          gen_op &D, gen_op& sigmap, int track)

	// Input	CPul	: A composite pulse (this)
        //              npts	: Number of FID points
	//		nSTs 	: Number of steps between points
        //              D 	: Detection operator
        //              sigmap	: Prepared density operator
	//		track	: Flag for tracking the computation
        // Output       data	: A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the composite
	//		  	  pulse CPul.
	// Note			: Assumes Relaxation Active!

  {
  if(!timeconst())				// There can be no step
    {						// synchronized acquisition
    CPulerror(60,1);				// if step lengths vary!
    CPulfatality(61);
    }
  double FIDcut = 1.e-6;			// FID cutoff intensity
  int FIDptcnt = 0;				// FID zero point count
  int FIDnilcnt = 4;				// Allowed FID zero points
  row_vector data(npts, complex0);		// Vector for FID points
  LSprop Gt = LSIprop();			// Needed Identity superprop.
  gen_op sigma;					// Needed operators
  D.Op_base(sigmap);				// Set D basis to sigmap's
  int nWFs = fullWFs(double(nSTs)*length(0));	// Waveforms in 1 increment
  nSTs -= nWFs*WFsteps;				// Steps left in 1 increment
  if(!nSTs)					// If no steps, use WF synch 
    return FIDRsynchWF(npts,nWFs,D,sigmap,track);
  if(track) FIDheader(5, 1); 			// Output header if tracking

  LSprop GWFs = Gt;				// Prop. for waveforms
  LSprop GWFtot = Gt;				// Prop. for total waveforms
  for(int j=0; j<nWFs; j++) GWFs *= GetGsum();	// Set prop for nWFs waveforms
  int iWFs = 0;					// Number of waveforms
  int iSTs = 0;					// Number of steps
  for(int i=0; i<npts; i++)			// Loop over FID points
    {						// Output point info 
    sigma = Gt.evolve(sigmap);			// This does the evolution
    data.put(trace(D,sigma),i);                 // Store FID point
    if(track)					// Output FID point if tracking
      {
      FIDpoint(5, i, iWFs, iSTs);
      FIDvalue(5, Gt.length(), data.get(i));
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
      iWFs += nWFs;				// Update waveform count
      GWFtot *= GWFs;				// Update WF evolution operator
      iSTs += nSTs;				// Update the step count
      if(iSTs >= WFsteps)			// If more step than full WF
        {					// then just use full WF
        iWFs++;					//	Add 1 WF to count
        GWFtot *= GetGsum();			//	Add 1 WF to prop
	iSTs -= WFsteps;			//	Sub 1 WF from step count
        }
      Gt = GWFtot;				// Set evolution under WFs
      if(iSTs) Gt *= GetGsum(iSTs);		// Add evolution under steps
      }
    }
  return data;
  }


row_vector PulComposite::FIDRsynchFR(int npts, int nFRs,
                                          gen_op &D, gen_op& sigmap, int track)
 
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

  {
  row_vector data(npts, complex0);              // Vector for FID points
  int HS = sigmap.dim();                        // Working Hilbert space
  LSprop Gt(HS*HS);				// Needed superpropagator
  gen_op sigma;                                 // Needed operators
  D.Op_base(sigmap);                            // Set D basis to sigmap's
  if(nFRs < 1)					// Can't allow less than 1
    {                                           // waveform fraction
    CPulerror(65,1);
    CPulfatality(61);
    }
  if(nFRs == 1)					// If only 1 fraction, use
    return FIDRsynchWF(npts,1,D,sigmap,track);	// waveform synch acquire
  double td = WFtp/nFRs;			// Dwell time to use
  int nSTs;					// Number of WF steps
  if(timeconst())				// Check if step synch is
    {						// possible too
    nSTs = fullsteps(td);			//   These steps in td 
    if(fabs(nSTs*length(0)-td) < cutzero) 	//   If steps span td
    return FIDRsynchST(npts,nSTs,D,sigmap,track);//   use step synch acquire
    }
  if(track) FIDheader(6, 1); 			// Output header if tracking
 
  LSprop *FrGs;
  FrGs = new LSprop[nFRs];				// Space for needed props
  LSprop Gsummed = Gt;
  double tr=0, te=td;				// Evolution time
  for(int m=0; m<nFRs; m++)			// Fill up props with evolve
    {						// over waveform fractions
    FrGs[m] = GetG(te);
    te += td;
    }

  int iWFs = 0;                                 // Number of waveforms
  int iSTs = 0;                                 // Number of steps
  te = 0;					// Zero FID evolve time
  int i,j;
  for(i=0; i<npts; i++)				// Loop over FID points
    {
    j = i%nFRs;					// Fraction index
    sigma = Gt.evolve(sigmap);			// This does the evolution
    data.put(trace(D,sigma),i);			// Store FID point
    if(track)					// Output FID point if tracking
      {
      FIDpoint(6, i, iWFs, iSTs);
      if(fabs(tr) > cutzero) printTime(std::cout, tr);
      else                   printTime(std::cout, 0);
      std::cout << "  ";
      FIDvalue(6, Gt.length(), data.get(i));
      }
    if(i != npts-1)				// Note iWFs, iSTs
      {
      if(!j && i) Gsummed *= FrGs[nFRs-1];	// Go to next waveform
      Gt = FrGs[j];				// Update evolution operator
      Gt *= Gsummed;				// Add previous WFs
      if(track)
        {
        te += td;				// Update FID evolution time
        iWFs = fullWFs(te);			// Update waveform count
        tr = te-iWFs*WFtp;			// Length remaining
        iSTs = fullsteps(te-iWFs*WFtp);		// Update step count
        tr -= sumlength(iSTs);			// Length remaining
        }
      }
    }
  delete [] FrGs;
  return data;
  }


row_vector PulComposite::FIDR(int npts, double td,
                                            gen_op &D, gen_op& sigmap, int trk)

	// Input	CPul	: A composite pulse (this)
        //              npts	: Number of FID points
        //              td	: Dwell time between FID points
        //              D 	: Detection operator
        //              sigmap	: Prepared density operator
	//		trk	: Flag for tracking the computation
        // Output       data	: A row vector contiaing an npts point FID
        //                        that was generated by detection with D
        //                        of the evolution of sigmap under the composite
	//		  	  pulse CPul.
	// Note			: Assumes Relaxation Active!
	// Note			: Calculation will be done in the basis of
	//			  the input density operator sigmap 

  {
  if(!sigmaeq.dim())				// If no relaxation is set
    return FID(npts, td, D, sigmap, trk);	// just calculate w/o it
  row_vector data(npts);			// Vector for FID points
  int HS = sigmap.dim();			// Working Hilbert space
  LSprop Gt, GI(HS*HS);				// Needed superpropagators
  gen_op sigma;					// Needed operators
  D.Op_base(sigmap);				// Set D basis to sigmap's
  double tr, te = 0;				// FID evolution time
  int nWFs, nSTs;				// Waveforms & steps used
  int synch = FIDtest(td,nWFs,nSTs,tr);		// Check FID synchronicity
  if(synch > 0)					// Here if Waveform synch.
    return FIDRsynchWF(npts,nWFs,D,sigmap,trk);
  else if(synch < 0)				// Here if step synch.
    {
    int nSTtot = nWFs*WFsteps + nSTs;
    return FIDRsynchST(npts,nSTtot,D,sigmap,trk);
    }
  if(trk) FIDheader(7, 1); 			// Output header if tracking

  double FIDcut = 1.e-6;                        // FID cutoff intensity
  int FIDptcnt = 0;                             // FID zero point count
  int FIDnilcnt = 4;                            // Allowed FID zero points

  for(int i=0; i<npts; i++)
    {
//		First Determine Which Propagation Steps Are Required
//		   (For Evolution From FID Start to Current Point)

    tr   = te;					// Evolve time this point
    if(fabs(tr) < cutzero) tr = 0;		// Zero time if below precision
    nWFs = fullWFs(tr);				// Full waveforms within tr
    tr  -= nWFs*WFtp;				// Time left to evolve
    if(fabs(tr) < cutzero) tr = 0;		// Zero time if below precision
    nSTs = fullsteps(tr);			// Additonal full steps
    tr  -= sumlength(nSTs);			// Time left to evolve
    if(fabs(tr) < cutzero) tr = 0;		// Zero time if below precision

    if(trk)
      {  
      FIDpoint(7, i, nWFs, nSTs);
      if(tr)
        {
        std::cout << Gdec(nSTs,6) << "     ";
        printTime(std::cout, tr);
        }
      else std::cout << "                        ";
      }
    if(tr < 0)                                  // Insure evolution time
      {                                         // is non-negative
      CPulerror(50,1);
      CPulfatality(51);
      }

//		Next We Must Build A Propagator In REVERSE Order
//		(And We Must Track Phase Changes Along The Way)

    Gt = GI;					// Start with I prop.
    if(tr) 
     Gt = LSprop(Leff(nSTs), SigSS(nSTs), tr);	// Prop of step i's evolution
    if(nSTs) Gt *= GetGsum(nSTs);		// Next evolve nSTs steps
    if(nWFs) Gt *= GetGmult(nWFs);		// Next evolve nWFs waveforms

//            Finally We Evolve Input Sigma To FID Point And Sample

    sigma = Gt.evolve(sigmap);			// This does the evolution
    data.put(trace(D,sigma),i);                 // Store FID point
    if(norm(data.get(i)) < FIDcut) FIDptcnt++;
    else                           FIDptcnt=0;
    if(FIDptcnt > FIDnilcnt)
      {  
      if(trk)
        std::cout << "\n\n\tFID Computation Into Baseline Noise Level ......";
      return data;
      }  
    te += td;					// Update FID evolution time
    if(trk) std::cout << "   ";
    if(trk) FIDvalue(7,Gt.length(),data.get(i));// Output data if tracking
    }
  return data;
  }  

 
// ____________________________________________________________________________
// H                 CLASS COMPOSITE PULSE PLOTTING FUNCTIONS
// ____________________________________________________________________________

                                                                                
/*
void       PulWaveform::getIdeal(double, ptt, int)  const;	INHERITED
row_vector PulWaveform::IvsT(int, int, int)         const;	INHERITED
row_vector PulWaveform::PvsT(int, int, int, double) const;	INHERITED
void       PulWaveform::GP(int, int, int)           const       INHERITED
void       PulWaveform::FM(int, int, int, int)      const;	INHERITED    */
 

// ____________________________________________________________________________
// I                 CLASS COMPOSITE PULSE I/O FUNCTIONS
// ____________________________________________________________________________

 
/*
std::ostream& PulWaveform::printBase(std::ostream &ostr)  const;		INHERITED
std::ostream& PulWaveform::printSteps(std::ostream &ostr) const;		INHERITED    */



std::ostream& PulComposite::printEvolve(std::ostream& ostr, double td) const
 
        // Input                PWF     : Pulse Waveform
        //                      td      : Evolution time (sec)
        //                      ostr    : Output stream
        // Output               none    : Pulse Waveform evolution info
        //                                is sent to the output stream

  {
  std::string lstart = "\n\t";
  int lgap = 30;   
  int stp = 1;                                  // Evolution Step
  double tev = 0;                               // Evolution time
  int plen;                                     // Print length

  ostr << "\n\n\t\t";
  if(WFname.length()) ostr << WFname;
  else              ostr << "\t";
  ostr << " Waveform Evolution Info\n";            // Output Header
  ostr << "\n\tSpecified Evolution Time:";
  ostr << std::string(16, ' ');
  printTime(ostr, td);
  ostr << "\n\tEvolution Spectral Width:";
  ostr << std::string(16, ' ');
  double SW = 1.0/td;
  printHz(ostr, SW);

//                   Output Evolution Under Full Waveform(s)

  int nWFs = fullWFs(td);			// # of full waveforms
  if(nWFs)					// Output any evolution
    {
    ostr << lstart << "Step " << stp << "."     //      Evolve step
         << Gdec(nWFs, 3) << " ";
    ostr << WFname << " Waveform";		//      Supercycle name
    plen = WFname.length() + 9;			//      Add a spacer
    if(nWFs > 1) { ostr << "s"; plen++; }
    ostr << std::string(lgap-plen, ' ');
    printTime(ostr, nWFs*WFtp);			//      Evolve time
    stp++;                                      //      Adjust step
    tev += nWFs*WFtp;				//      Time evolved
    }
 
//                    Output Evolution Under Waveform Steps
//			 (For Less Than 1 Full Waveform)
 
  int nS = fullsteps(td-tev); 			// # of full waveform steps
  if(nS)					// Output any evolution
    {                                           // Under the waveform steps
    ostr << "\n\tStep " << stp << "."
         << Gdec(nS, 3) << " ";
    ostr << WFname << " Waveform Step";         //      Waveform name
    plen = WFname.length() + 14;		//      Add a spacer
    if(nS > 1) { ostr << "s"; plen++; }
    ostr << std::string(lgap-plen, ' ');
    printTime(ostr, sumlength(nS));
    stp++;                                      //      Adjust step
    tev += sumlength(nS);			//      Time evolved
    }

//            Output Evolution Under Partial Waveform Steps
 
  double tadd;
  if(fabs(td-tev) > cutzero)			// It there is a partial step
    {                                           // output any evolution
    ostr << "\n\tStep " << stp << "."
         << Gdec(1, 3) << " ";
    ostr << WFname << " Partial Step ";//      Waveform name
    plen = WFname.length() + 14;		//      Add a spacer
    ostr << std::string(lgap-plen, ' ');
    tadd = td - tev; 
    printTime(ostr, tadd);			//      Time evolved
    tev += tadd;
    }
  ostr << "\n\tSummed Evolution Time:";
  ostr << std::string(19, ' ');
  printTime(ostr, tev);
  ostr << "\n\tRemainder Evolution Time:";
  ostr << std::string(16, ' ');
  printTime(ostr, td-tev);
  return ostr;
  }

 
std::ostream& PulComposite::printFID(std::ostream &ostr,
                                        double td, int npts) const
 
        // Input        PT    : A pulse train (this)
        //              ostr  : An output stream
        //              td    : Dwell time between FID points
        //              npts  : Number of FID points
        // Output             : Information regarding the FID generation
        //                      is set to the output stream ostr
 
  {
  if(!WFtp)
    {
    ostr << "\n\n\t\tEmpty Waveform, No Acquisition Possible\n\n";
    return ostr;
    }
  ostr << "\n\n\t\t";
  if(WFname.length()) ostr << WFname;
  else              ostr << "\t";
  ostr << " Waveform Acquisition Info\n";            // Output Header
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
  int nWFs, nS;
  for(int i=0; i<npts; i++)
    {

//                   Output Evolution Under Full Waveform(s)

    tev = 0;					// Zero evolution time
    ostr << lstart << "Point Evolution Time:";
    ostr << std::string(20, ' ');
    printTime(ostr, tacq);
    nWFs = fullWFs(tacq);			// # of full waveforms
    if(nWFs)					// Output any evolution
      {
      ostr << lstart << "Step " << stp << "."	//      Evolve step
           << Gdec(nWFs, 3) << " ";
      ostr << WFname << " Waveform";		//      Supercycle name
      plen = WFname.length() + 9;		//      Add a spacer
      if(nWFs > 1) { ostr << "s"; plen++; }
      ostr << std::string(lgap-plen, ' ');
      printTime(ostr, nWFs*WFtp);		//      Evolve time
      stp++;                                    //      Adjust step
      tev += nWFs*WFtp;				//      Time evolved
      }
 
//                    Output Evolution Under Waveform Steps
//			 (For Less Than 1 Full Waveform)
 
    nS = fullsteps(tacq-tev); 			// # of full waveform steps
    if(nS)					// Output any evolution
      {                                         // Under the waveform steps
      ostr << "\n\tStep " << stp << "."
           << Gdec(nS, 3) << " ";
      ostr << WFname << " Waveform Step";       //      Waveform name
      plen = WFname.length() + 14;		//      Add a spacer
      if(nS > 1) { ostr << "s"; plen++; }
        ostr << std::string(lgap-plen, ' ');
      printTime(ostr, sumlength(nS));
      stp++;                                    //      Adjust step
      tev += sumlength(nS);			//      Time evolved
      }

//            Output Evolution Under Partial Waveform Steps
 
     if(fabs(tacq-tev) > cutzero)		// It there is a partial step
       {                            		// output any evolution
       ostr << "\n\tStep " << stp << "."
            << Gdec(1, 3) << " ";
       ostr << WFname << " Partial Step ";	//      Waveform name
       plen = WFname.length() + 14;		//      Add a spacer
       ostr << std::string(lgap-plen, ' ');
       printTime(ostr, tacq-tev);		//      Time evolved
       }
 
    ostr << "\n\tSample"
         << Gdec(i,4) << " Acquisition Point Time";
    ostr << std::string(8, ' ');
    printTime(ostr, tacq);
    tacq += td;
    ostr << "\n";
    }
  return ostr;
  }

 
std::ostream& PulComposite::printInfo(std::ostream &ostr) const

	// Input		CPul	: Composite Pulse Composite
        //                      ostr	: Output stream
	//			full	: Flag for output amount
        // Output               none	: CPul storage info is sent
	//				  to the output stream

  {
  ostr << "\n\tComposite Hamiltonians:           ";
  if(Hsteps)
    {
    ostr << "Present (" << Hcount;
    int con = WFsteps - Hcount;
    if(con) ostr << ", " << con << " conserved";
    ostr << ")";
    }
  else       ostr << "Absent";
  ostr << "\n\tComposite Propagators:            ";
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
  ostr << "\n\tComposite SuperPropagators:       ";
  if(Gsteps) ostr << "Present";
  else       ostr << "Absent";
  return ostr;
  }


std::ostream& PulComposite::print(std::ostream &ostr, int full) const

	// Input		CPul	: Composite Pulse Composite
        //                      ostr	: Output stream
	//			full	: Flag for output amount
        // Output               none	: Composite Pulse info is sent
	//				  to the output stream

  {
  if(!WFsteps)
    {
    ostr << "\n\n\t\t\tEmpty Composite Pulse\n\n";
    return ostr;
    }
  ostr << "\n\n\t\t\t  Composite Pulse " << name() << "\n";
  ostr << "\n\tPulse Channel:                    " << Iso;
  printBase(ostr);
  if(full)
    {
    printSteps(ostr);
    ostr << "\n";
    printInfo(ostr);
    }
  ostr << "\n\n"; 
  return ostr;
  }

            
std::ostream &operator << (std::ostream &ostr, const PulComposite &CPul)

	// Input		CPul	: Composite Pulse Composite
        //                      ostr	: Output stream
        // Output               none	: Composite Pulse info is sent
	//				  to the output stream

  {
  CPul.print(ostr);
  return ostr;
  }


#endif						// PulComposite.cc
