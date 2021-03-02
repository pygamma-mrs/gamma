/* PulSinc.cc **************************************************-*-c++-*-*
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Sinc Shaped Pulses                        Implementation	**
**                                                                      **
**      Copyright (c) 1995                                              **
**      Dr. Scott A. Smith                                              **
**      1800 E. Paul Dirac Drive                                        **
**      National High Magnetic Field Laboratory                         **
**      Tallahassee FL 32306 USA                                        **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
** This GAMMA module deals with sinc shaped pulses. Herein are          **
** functions which compute propagators for the facile use sinc pulses   **
** as single steps in GAMMA programs. Included are functions to analyze **
** sinc pulse details, return the sinc waveform as a vector, and        **
** directly apply such a pulse to a spin system.                        **
**                                                                      **
*************************************************************************/

#ifndef _PulSinc_cc_			// Is file already included?
#define _PulSinc_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Pulses/PulSinc.h>		// Include the interface
#include <Basics/ParamSet.h>		// Include parameter sets
#include <HSLib/SpinSystem.h>		// Include spin systems
#include <HSLib/HSprop.h>		// Include propagators
#include <HSLib/HSauxil.h>		// Include propagators
#include <HSLib/GenOp.h>		// Include operators
#include <HSLib/HSprop.h>		// Include operators
#include <LSLib/SuperOp.h>		// Include superoperators
#include <HSLib/HSham.h>		// Knowledge of isotropic Hamiltonians
#include <Level2/acquire1D.h>		// Include knowledge of "acquisitions"
#include <HSLib/SpinOpCmp.h>		// Include composite spin ops
#include <HSLib/SpinOpRot.h>		// Include spin rotation operators
#include <string>			// Include libstdc++ strings
#include <Basics/StringCut.h>		// Include Gform and Gdec
#include <stdlib.h>

using std::cout;
using std::cin;

// ____________________________________________________________________________
//                         Sinc Pulse Hamiltonians
// ____________________________________________________________________________

/* These functions are meant to generate an array of Hamiltonians for a Sinc
   pulse waveform.  Each step of the Sinc waveform will have a Hamiltonian,
   For all of these function the Sinc waveform must be specified, that means
   that we need two of { pang, gamB1, tp } as well as the end node and number
   of pulse steps.  In addition, the function will require a static Hamiltonian
   (e.g. the isotropic Hamiltonian) and the operator for the spin component
   of the pulse rf-field (e.g. Fx)                                           */

void SincPulseHs(gen_op* Hs, gen_op& H0, gen_op& Fxy,
                                        int N, double ang, double tp, int node)

        // Input        Hs      : Array of general operators
        //              H0      : Static Hamiltonian *
        //              Fxy     : RF-Field Hamiltonian *
        //              N       : Number of Sinc steps
        //              ang     : Pulse rotation angle (degrees)
        //              tp      : Sinc pulse length (sec)
        //              node 	: Sinc function end node
        // Output       void    : Hamiltonians for a Sinc pulse
        //                        are put into array Hs

  {
  double gamB1 = SincGamB1(ang, tp, N, node);	// Get field strength
  row_vector Sinc = SincVect(gamB1, N, node);	// RF intensities (waveform)
  Fxy.Op_base(H0);				// Set Fxy basis for easy add
  for(int i=0; i<N; i++)			// Loop over Sinc steps
    {						// For the latter times we use
    if(N-1-i < i) Hs[i] = Hs[N-1-i];		// symmetry to reuse the Hs
    else          Hs[i] = H0 - Sinc.get(i)*Fxy;	// These are the first half
    }   					// of the waveform where we
  }						// must must calculate the Hs
      

// ____________________________________________________________________________
// B                        Sinc Pulse Propagators
// ____________________________________________________________________________

/* These functions are meant to generate an array of propagators for a Sinc     
   pulse waveform.  Each step of the Sinc waveform will have a propagator Ui.
   For all of these function the Sinc waveform must be specified, that means
   that we need two of { pang, gamB1, tp } as well as the end node and number
   of pulse steps.  In addition, the function will require a static Hamiltonian
   (e.g. the isotropic Hamiltonian) and the operator for the spin component  */
 
void SincPulseUs(gen_op* Us, gen_op& H0, gen_op& Fxy,
                                     int N, double ang, double tp, int node)
 
        // Input        Us      : Array of general operators
        //              H0      : Static Hamiltonian *
        //              Fxy     : RF-Field Hamiltonian *
        //              N       : Number of Sinc steps
        //              ang     : Pulse rotation angle (degrees)
        //              tp      : Sinc pulse length (sec)
        //              node 	: Sinc function end node
        // Output       void    : Propagators for a Sinc pulse
        //                        are put into array Us
 
  {
  double tdiv = tp/double(N);			// Incremental prop. time
  double gamB1 = SincGamB1(ang, tp, N, node);	// Get field strength
  row_vector Sinc = SincVect(gamB1, N, node);	// RF intensities (waveform)
  Fxy.Op_base(H0);				// Set Fxy basis for easy add
  gen_op H;					// Scratch operator
  for(int i=0; i<N; i++)			// Loop over Sinc steps
    {						// For latter times we use
    if(N-1-i < i) Us[i] = Us[N-1-i];		// symmetry to reuse the Hs
    else					// These are the first half
      { 					// of the waveform. We must
      H = H0 - Sinc.get(i)*Fxy;			// calculate the H & U for
      Us[i] = prop(H, tdiv);			// each of these pulse steps
      }
    }
  }    

// ____________________________________________________________________________
// C                        Sinc Pulse Full Propagator
// ____________________________________________________________________________

/*
  Note that these use the symmetry and only calculate step Hamiltonians and
  propagators over half of the waveform.  That makes it faster than using 
  the functions which return Hamiltonian or Propagator arrays for all steps
  then summing/multiplying over them.                                        */

gen_op SincPulseU(gen_op& H0rot, gen_op& Fxy,
                                        int N, double ang, double tp, int node)
 
        // Input        H0rot   : Static Hamiltonian *
        //              Fxy     : RF-Field Hamiltonian *
        //              N       : Number of Sinc steps
        //              ang     : Pulse rotation angle (degrees)
        //              tp      : Sinc pulse length (sec)
        //              node 	: Sinc function end node
        // Output       U       : Propagator for a Sinc pulse
        // Note                 : HOrot is the static Hamiltonian
        //                        in the rotating frame of the pulse.
        // Note                 : Fxy implicitly contains rf-pulse phase
 
/*                      U     = U       * U     * U
                         Sinc    t>tp/2    tp/2    t<tp/2                    */

  {
  double dN = double(N);			// Pulse steps
  double tdiv = tp/dN;				// Inc. prop. time
  double gamB1 = SincGamB1(ang, tp, N, node);	// Get field strength
  row_vector Sinc = SincVect(gamB1, N, node);	// Sinc RF intensities
  Fxy.Op_base(H0rot);				// Set Fxy basis for easy add
  double gB1 = Sinc.getRe(0);			// Field strength of 1st step
  gen_op H = H0rot + gB1*Fxy; 			// Total Ham., 1st step
  gen_op Utlow = prop(H, tdiv);			// Propagator for 1st step
  gen_op Uthigh = Utlow;			// Propagator for Nth step
  gen_op U, Ustep;				// Working propagators
  int i;
  for(i=1; i<N-1-i; i++)			// Loop 1st half of Sinc
    {						// steps, symmetry for 2nd half
    gB1 = Sinc.getRe(i);			// RF-field strength this step
    H = H0rot + gB1*Fxy;			// Total Hamiltonian, this step
    Ustep = prop(H, tdiv);			// Prop. for this pulse step
    Utlow &= Ustep;				// Utlow = Ustep*Utlow
    Uthigh *= Ustep;				// Uthigh = Uthigh*Ustep
    }
  if(N == 1) U = Utlow;				// For one step, U is 1st step
  else if(i == N-1-i)				// For odd number steps (odd N)
    {						// still need to do middle step
    gB1 = Sinc.getRe(i);			// RF-field of middle step
    H = H0rot + gB1*Fxy;			// Total Hamiltonian, this step
    U = prop(H,tdiv);				// Prop. for middle step
    U *= Utlow;					// Add in after earlier steps
    U &= Uthigh;				// Add later steps after middle
    }
  else U = Uthigh*Utlow;			// Even # steps(even N), this U
  return U;
  }


gen_op SincPulseU(gen_op& H0rot, gen_op& Fxy, const SincPulDat& SD)

 
        // Input        H0rot   : Static Hamiltonian
        //              Fxy     : RF-Field Hamiltonian
        //              SD	: Sinc pulse structure
        // Output       U       : Propagator for a Sinc pulse
        // Note                 : HOrot **MUST** be in the rotating frame
        //                        of the applied pulse
        // Note                 : Fxy **MUST** contain the phase of
        //                        the applied pulse
        // Note                 : In line with the previous two notes,
        //                        SD's phase and offset are not used
        //                        herein (they're assumed in H0 & Fxy)
 
   {
   double tp = SD.tau; 				// Pulse length (sec)
   int N = SD.N;				// Pulse steps
   double gamB1 = SD.gamB1;			// Pulse strength (Hz)
   int node = SD.node;			// Pulse end node
   double ang = SincAngle(gamB1, tp, N, node);	// Pulse angle (degrees)
   return SincPulseU(H0rot,Fxy,N,ang,tp,node);	// Use function overload
   }
 

// ____________________________________________________________________________
// D                      Sinc Pulse Access Functions
// ____________________________________________________________________________

/* A Sinc pulse angle is defined as the total rotation angle experienced by
   a single spin on resonance.  That angle is determined by the applied
   RF-field strength (gamB1), the length of the pulse, and (in this discrete
   treatment) the number of pulse steps.  As such, if we know any of the two
   quantities { angle, gamB1, tp } and the number of steps we have defined
   a Sinc pulse. The following functions allow users to access these values.

   Function                                Purpose
   ---------         -------------------------------------------------------
   SincAngle	     Given {gamB1, tp,    N } returns effective pulse angle
   SincGamB1	     Given {pang,  tp,    N } returns Sinc rf-field strength
   SincTime	     Given {pang,  gamB1, N } returns Sinc pulse length   */


double SincAngle(double gamB1, double tau, int N, int node)
   {
   row_vector Sinc = SincNVect(N, node);	// Normalized Sinc function
   double angle = Re(Sinc.sum())*gamB1*tau;	// On resonance pulse angle
   return angle*360.0/double(N);		// Return angle in degrees
   }

double SincGamB1(double angle, double tau, int N, int node)
   {
   row_vector Sinc = SincNVect(N, node);	// Normalized Sinc function
   double normang=Re(Sinc.sum())*tau/double(N);	// Normalized rot. angle
   return (angle/360.0)/normang;		// Return gamB1 in Hz
   }

double SincTime(double angle, double gamB1, int N, int node)
   {
   row_vector Sinc = SincNVect(N, node);	// Normalized Sinc function
   double normang = gamB1*Re(Sinc.sum());	// Get Sinc integral (Hz)
   return double(N)*angle/(360.0*normang);	// Scale to seconds
   }

// ____________________________________________________________________________
// E                     Sinc Pulse WaveForm Functions
// ____________________________________________________________________________

/* These returns a vector filled with a Sinc function that is centered about
   the vector mid-point and spanning "node" nodes on each side.  It does
   using the discrete sinc function given by

                                sin[PI*node*(2i-N+1)/N-1]
                       S  = X * -------------------------
                        i        [PI*node*(2i-N+1)/N-1]

   For a normalized Sinc function X will be 1, and for a Sinc associated with a
   pulse X will be the applied RF field strength, gamB1.

        	Input	gamB1   : The rf-field strength (Hz)
			tau	: Pulse length (sec)
        		N       : Number of Sinc steps
			node    : Cutoff node [1,2,3,...] (default 3)
			endzero : Flag whether to have ends at zero

   Function                               Purpose
   ---------         --------------------------------------------------------
   SincNVect 	     Given a # of points & a node, returns sinc with max 1
   SincVect	     Given a # of pts, a node, & a gamB1, sinc with max gamB1
   SincIntVec        Given a gamB1, tau, # pts & node, pulse ang.vs time (deg)

  An additonal aspect worthy of mention is the flag endzero.  These functions
  will end Sinc at a node, i.e. at zero intensity.  However, for sinc pulses
  it doesn't make any sense to have a zero intensity pulse step.  So, if the
  flag endzero is NOT set it will assume points -1 and N are they that reside
  at the nodes.  Then Sinc pulses will not (very often) have no zero steps.   */

row_vector SincNVect(const SincPulDat& SD, int endzero)
  { return SincNVect(SD.N, SD.node, endzero); }

row_vector SincNVect(int N, int node, int endzero)
   {
   if(!endzero) N+=2;			// Add 2 points if not node end
   row_vector Sinc(N);			// Vector of Sinc points
   int Nm1 = N-1;			// Need N-1 for each point
   double dNm1 = double(Nm1);		// Need double N-1 also
   double pin = PI*node;		// Need PI*node for each point
   double arg;				// For sin(arg)/arg
   for(int i=0; i<N; i++)		// Loop over Sinc steps
     {
     if(Nm1-i < i)			// Use symmetry on 2nd half avoid
       Sinc.put(Sinc.get(Nm1-i), i);	// recalculating same pts in 1st
     else if(double(i) == dNm1/2)	// Whereas on the mid-point set
       Sinc.put(1,i);			// the intensity to one
     else				// For the 1st half we calculate
       {				// the sinc function intensities
       arg=pin*(2*double(i)-dNm1)/dNm1;	// Indexed so Sinc is at mid-pulse
       Sinc.put(sin(arg)/arg, i);	// Store the Sinc intensity
       } 
     }
   if(!endzero)				// Clip last two points (nodes)
     Sinc = Sinc.get_block(0,1,1,N-2);	// if we don't want node at ends
   return Sinc;
   }

row_vector SincVect(const SincPulDat& SD, int endzero)
  { return SincVect(SD.gamB1, SD.N, SD.node, endzero); }

row_vector SincVect(double gamB1, int N, int node, int endzero)
   { return gamB1*SincNVect(N, node, endzero); }

row_vector SincIntVec(const SincPulDat& SD, int endzero)
  { return SincIntVec(SD.gamB1, SD.tau, SD.N, SD.node, endzero); }

row_vector SincIntVec(double gamB1, double tau, int N, int node, int endz)
  {
  row_vector Sinc = SincVect(gamB1,N,node,endz);// Vector of Sinc points
  row_vector SIvect(N+1);			// Vector for integral points
  SIvect.put(complex0, 0);			// First point of int. is 0
  double tdiv = tau/double(N);			// Single step time
  double time = tdiv;				// Time at point i=1
  double Int = 0.0;				// Integral at point i
  double Z = tdiv*360.0;			// Reuse Z for conversion
  for(int i=0; i<N; i++)			// Loop over Sinc steps
    {
    time += tdiv;				// Time at point i
    Int += Sinc.getRe(i)*Z;			// Integral (in degrees)
    SIvect.put(complex(time,Int), i+1);		// Store integral up to here
    }
  return SIvect;
  }

// ____________________________________________________________________________
// F                    Sinc Pulse Parameter Set Functions
// ____________________________________________________________________________
 
/* These functions read in Sinc pulse parameters from an external ASCII file 
   (a GAMMA parameter set file).   
 
     Function                 Purpose                     Parameters Read 
  ---------------  -------------------------------  ---------------------------
  SincSteps        Reads in # of Sinc pulse steps          SincSteps
  SincNode         Reads in Sinc pulse end node            SincNode
  SincStrength     Reads in Sinc pulse strength     SincAng, SincLen, SincGamB1
  SincSelectivity  Reads in Sinc pulse selectivity   SincSpin, SincW, SincIso
  SincPhase        Reads in Sinc pulse phase               SincPhi
  ReadSinc         Reads in Sinc pulse                                       */


int SincSteps(const ParameterSet& pset, int idx, int pf)
  {
//                 Read The Number of Sinc Pulse Steps
 
  ParameterSet::const_iterator item;         // A pix into parameter list
  std::string pname, ssfile, pstate;                 // Items in each pset entry
  std::string SI = std::string("(") + Gdec(idx)            // Name adjustment if indexed
            + std::string(")");
  int npt;
  pname = std::string("SincSteps");			// Number of Sinc steps
  if(idx >= 0) pname += SI;                     // Adjust name if indexed
  item = pset.seek(pname);                      // Pix in parameter list
  if(item != pset.end())
    (*item).parse(pname,npt,pstate);            // Get Sinc steps
  else
    {
    if(pf)
      cout << "\n\tCant Find " << pname << " in "
           << "Parameter Set.  Setting 501 Sinc Steps";
    npt = 500;
    }
  return npt;
  }


int SincNode(const ParameterSet& pset, int idx, int pf)
  {
//                 Read The Number of Sinc Pulse Steps
 
  ParameterSet::const_iterator item;         // A pix into parameter list
  std::string pname, ssfile, pstate;                 // Items in each pset entry
  std::string SI = std::string("(") + Gdec(idx)            // Name adjustment if indexed
            + std::string(")");
  int nd=3;
  pname = std::string("SincNode");			// Sinc pulse node
  if(idx >= 0) pname += SI;                     // Adjust name if indexed
  item = pset.seek(pname);                      // Pix in parameter list
  if(item != pset.end())			// Get node parameter
    (*item).parse(pname,nd,pstate);
  else
    {
    if(pf)
      cout << "\n\tCant Find " << pname << " in "
           << "Parameter Set.   Setting Sinc Node to 3";
    nd = 3;
    }
  return nd;
  }


int SincStrength(const ParameterSet& pset, SincPulDat& SD, int idx, int pf)
  {
//    Try & Read Two of Three: {Pulse Angle, RF Strength, Pulse Length} 

  ParameterSet::const_iterator item;         // A pix into parameter list
  std::string pname, ssfile, pstate;                 // Items in each pset entry
  std::string SI;					// Name adjustment if indexed
  if(idx >= 0) SI = std::string("(") + Gdec(idx)	
                  + std::string(")");
  double gamB1, time, pang;			// Sync parameters
  int G=0, T=0, A=0;				// Flags if found
  pname = std::string("SincAng") + SI;		// Sync pulse pang
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// Get pang parameter
    {
    (*item).parse(pname,pang,pstate);		//   Read in pulse pang
    A = 1;					//   Flag we know the pang
    }
  pname = std::string("SincLen") + SI;		// Sync pulse length
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// Get pang parameter
    {
    (*item).parse(pname,time,pstate);		//   Read in pulse length
    T = 1;					//   Flag we know the length
    }
  if(A+T != 2)
    {
    pname = std::string("SincGamB1") + SI;		// Sync pulse pang
    item = pset.seek(pname);			// Pix in parameter list
    if(item != pset.end())			// Get gamB1 parameter
      {
      (*item).parse(pname,gamB1,pstate);	//   Read in the pulse strength
      G = 1;					//   Flag we know the strength
      }
    }

//     We Should Now Have 2/3: {Pulse Angle, RF Strength, Pulse Length} 

  if(A+T+G == 0) return 0;			// Quit if we know none of em
  if(A+T+G >= 2)				// If we know 2 of 3 then set
    {						// {time,gamB1} appropriately
    int node = SD.node;				// Get Sinc node used
    int npt  = SD.N;				// Get Sinc steps used
    if(A && T)
      {
      SD.tau = time;				// Set the time
      SD.gamB1=SincGamB1(pang,time,npt,node);	// Set field strength
      }
    else if(G && T)
      {
      SD.tau = time;				// Set the time
      SD.gamB1 = gamB1;				// Set the field strength
      }
    else if(A && G)
      {
      SD.tau = SincTime(pang,gamB1,npt,node);// Set the time
      SD.gamB1 = gamB1;				// Set the field strength
      }
    return 1;					// All is just dandy
    }

//     Perhaps We Only Have 1/3: {Pulse Angle, RF Strength, Pulse Length} 

  if(A+T+G == 1)
   {
   int node = SD.node;				// Get Sinc node used
   int npt  = SD.N;				// Get Sinc steps used
   if(A)					// We only know pulse angle!
     {						// Issue a warning then use
     pname = std::string("SincLen") + SI;		// a default pulse length &
     if(pf)					// set {time,gamB1} 
       cout << "\n\t\tCant Find " << pname
            << " in Parameter Set. "
            << "Setting Sync Pulse Length"
            << " to 30 msec";
     SD.tau = 0.03;				// Set the time
     SD.gamB1=SincGamB1(pang,time,npt,node);	// Set field strength
     }
   if(T)
     {
     pname = std::string("SincAng");
     if(idx >= 0) pname += SI;
     if(pf)
       cout << "\n\t\tCant Find " << pname
            << " in Parameter Set. "
            << "Setting Sync Pulse Angle"
            << " to 90 degrees";
     SD.tau = time;				// Set the time
     SD.gamB1=SincGamB1(90.0,time,npt,node);	// Set field strength
     }
   if(G)
     {
     pname = std::string("SincLen");			// Sync pulse length
     if(idx >= 0) pname += SI;			// Adjust name if indexed			
     if(pf)
       cout << "\n\t\tCant Find " << pname
            << " in Parameter Set. "
            << "Setting Sync Pulse Length"
            << " to 30 msec";
     SD.tau = 0.03;				// Set the time
     SD.gamB1 = gamB1;				// Set the field strength
     }
   return 1;
   }
 cout << "\n\tProblems Reading Sync Pulse parameters\n";
 return 0;
 }


int SincSelectivity(const ParameterSet& pset, const spin_system& sys,
                                            SincPulDat& Sdata, int idx, int pf)

/*                        Read Sinc Selectivity

   If "SincSpin" is declared, apply pulse at shift(SincSpin).  If "SincSpin"
   is not present, we need both the offset frequency and isotope selectivity.
   The first can be input in either Hz or PPM, PPM being preferable with the
   parameters SincWPPM or SincW.  The second is input using the parameter
   SincIso.  Note that the isotope selectivity need not be set in a 
   homonuclear spin system.                                                  */

  {
  ParameterSet::const_iterator item;         // A pix into parameter list
  std::string pname, ssfile, pstate;                 // Items in each pset entry
  std::string SI;					// Name adjustment if indexed
  if(idx >= 0) SI = std::string("(") + Gdec(idx)	
                  + std::string(")");

//		First Look for Selectivity Via Spin Index

  int spin;					// For Sinc spin freq.
  pname = std::string("SincSpin");			// Sinc spin selectivity
  item = pset.seek(pname+SI);			// Pix in parameter list
  if(item != pset.end())                        // If we have found the spin
    {						// then set selectivity
    (*item).parse(pname,spin,pstate);		//    Get spin
    Sdata.Wrf = sys.shift(spin);                //    Sinc pulse offset
    Sdata.Iso = sys.symbol(spin);               //    Sinc isotope channel
    return 1;
    }

//	      Next Look for Selectivity Via { W , Iso } Pair
//		  First We Go For The Isotope Channel

  if(sys.homonuclear())				// On homoncuclear systems
    Sdata.Iso = sys.symbol(0);			// the selectivity is clear
  else						// We must ask for it in
    {						// heteronuclear systems
    std::string Siso;					// For Sinc isotope type
    pname = std::string("SincIso");			// Now look for selectivity
    item = pset.seek(pname+SI);			// Pix in parameter list
    if(item == pset.end())			// We found isotope channel
      {
      if(pf)
        {
        cout << "\n\tCant Find " << pname << " in Parameter Set";
        cout << "\n\tThis Must Be Specified in a Heteronuclear System!";
        }
      return 0;
      }
    (*item).parse(pname,Siso,pstate);		//   Get isotope channel
    Sdata.Iso = Siso;				//   Set the Sync channel
    if(sys.symbol(0)!=Siso)			//   Insure channel is in
      {						//   the spin system!
      cout << "\n\tYou've Specified A Sync"
           << " Pulse On The " << Siso
           << " Isotope Channel.\n\tYour Spin"
           << " System Contains No Spins Of That"
           << " Type!\n\tDon't Expect Much Effect"
           << " From The Pulse...";
      }
    }

//      At This Point We Have  Iso Of { W , Iso } Pair, Look For W 

  double W;					// For Sinc offset value
  pname = std::string("SincWPPM");			// Sync frequency in PPM
  item = pset.seek(pname+SI);			// Pix in parameter list
  if(item != pset.end()) 			// Do we have carrier freq?
    {
    (*item).parse(pname,W,pstate);		// Get the carrier frequency
    Sdata.Wrf = W*sys.Omega(Sdata.Iso);		// Set the Sync frequency
    }
  else						// If we can't get the offset
    {						// in PPM look for it in Hz
    pname = std::string("SincW");			// Sync frequency in Hz
    item = pset.seek(pname+SI);			// Pix in parameter list
    if(item != pset.end()) 			// Do we have carrier freq?
      {
      (*item).parse(pname,W,pstate);		// Get the carrier frequency
      Sdata.Wrf = W;				// Set the Sync frequency
      }
    else					// We cannot find the offset 
      {						// frequency, so we will just
      if(pf)					// set it for no offset
        cout << "\n\tCant Find " << pname 
             << " in ParameterSet."
             << " Setting Sync Frequency to 0";
      Sdata.Wrf = 0.0;
      }
    return 1;
    }
  return 0;					// Should never reach here
  }


double SincPhase(const ParameterSet& pset, int idx, int pf)
 
  {
  ParameterSet::const_iterator item;         // A pix into parameter list
  std::string pname, ssfile, pstate;                 // Items in each pset entry
  std::string SI;					// Name adjustment if indexed
  if(idx >= 0) SI = std::string("(") + Gdec(idx)	
                  + std::string(")");
  double phi;
  pname = std::string("SincPhi") + SI;		// Sync phase angle
  item = pset.seek(pname);                      // Pix in parameter list
  if(item != pset.end())
    (*item).parse(pname,phi,pstate); 		// Get Sync phase
  else
    {
    if(pf)
      cout << "\n\tCant Find " << pname 
           << " in Parameter Set."
           << " Setting Sync Pulse Phase to 0";
    phi = 0.0;
    }
  return phi;
  }


SincPulDat ReadSinc(const std::string& Fin, const spin_system& sys, int idx, int pf)
 
        // Input        Fin	: A filename
        //              sys     : Spin system
        //              idx     : Sinc pulse index (default none)
        //              pf	: Print flag (default 1=print)
        // Output       SincData: Structure defining the Sinc pulse

  {
  ParameterSet pset;                            // A parameter set
  pset.read(Fin);				// Read pset in
  SincPulDat SincData;				// Sinc parameters
  SincData.N    = SincSteps(pset, idx, pf);	// Get number of steps
  SincData.node = SincNode(pset, idx, pf);	// Get the sinc node at end
  int TF=SincStrength(pset,SincData,idx,pf);	// Get 2 of {gamB1,tp,pang}
  if(!TF)					// Fail if can't setup Sinc 
    {
    std::string pname;
    std::string SI;
    if(idx>=0) SI += std::string("(") + Gdec(idx)	// Name adjustment if indexed
                  +  std::string(")");
    cout << "\n\n\tTo Define a Sync Pulse You "
         << "MUST Set 2 of the Following 3 "
         << "Parameters in " << Fin << ":";
    pname = std::string("SincAng");
    cout << "\n\n\t\t" << pname+SI << " Sync"
         << " Pulse Angle (Degrees)";
    pname = std::string("SincLen");	
    cout << "\n\t\t" << pname+SI << " Sync"
         << " Pulse Length (Seconds)";
    pname = std::string("SincGamB1");
    cout << "\n\t\t" << pname+SI << " Sync"
         << " Pulse Strength (Hz)";
    cout << "\n\nSorry, Cannot Continue.\n\n";
    exit(-1);
    }
  SincSelectivity(pset,sys,SincData,idx,pf);	// Get sinc selectivity
  SincData.phi = SincPhase(pset, idx, pf);	// Get the sinc phase
  return SincData;
  }


// ____________________________________________________________________________
// F                       Sinc Pulse ASCII Input Functions
// ____________________________________________________________________________


row_vector SincHistogram(double gamB1, double tau, int N, int node)

        // Input        gamB1   : The rf-field strength (Hz)
        //              tau     : Sinc pulse length (sec)
        //              N       : Number of Sinc steps
	//		node    : Cutoff node [1,2,3,...]
	//			  defaults to 3
        // Output       Sshape  : Vector of point for plot of Sinc
	//			  as a histogram.

   {
double den, Z;				// Denominator & numerator
SincPrep(N, node, den, Z);			// Prep. Sinc values
   double tdiv = tau/double(N);                 // Incremental time
   int M = N+1;
   row_vector SincShape(2*N+M);			// Vector of Sinc points
   double lastv, lastt;
   int I = 0;
   row_vector Sinc = SincVect(gamB1,N,node);	// Sinc intensities
   double time = 0.0;				// Time in pulse
   complex z;
   for(int i=0; i<N; i++)			// Loop over Sinc steps
     {
     if(i)
       {
       z = complex(time,Sinc.getRe(i-1));	// Sinc function point (x,y)
       SincShape.put(z, I++);			// For horizontals in hist.
       if(M) SincShape.put(complex(time),I++);	// For verticals in hist.
       }
     else if(M)
       SincShape.put(complex(time,0),I++);	// First vertical in hist.
     z = complex(time,Sinc.getRe(i));		// Sinc function point (x,y)
     SincShape.put(z,I++);			// Sinc intensity
     lastv = Sinc.getRe(i);			// Store the previous intensity
     lastt = double(i);				// Store the previous point
     if(i==N-1)					// For the last point
       {
       time += tdiv;				// Increase time increment
       z = complex(time,Sinc.getRe(i));		// Sinc function point (x,y)
       SincShape.put(z, I++);			// For last horizontal in hist.
       if(M) SincShape.put(time, I++);		// For last vertical in hist.
       }
     time += tdiv;				// Increment to next time
     }						// (for  horizontal axis)
   return SincShape;
   }

// ____________________________________________________________________________
// H                 Sinc Pulse Interactive Setup Functions
// ____________________________________________________________________________

void SincPts(int argc, char* argv[], int& qn, SincPulDat& SD)
  { query_parameter(argc,argv,qn++,"\n\tNumber of Points in Sinc? ", SD.N);   }

void SincNode(int argc, char* argv[], int& qn, SincPulDat& SD)
  { query_parameter(argc,argv,qn++,"\n\tEnd Node (1, 2, 3, ...)? ", SD.node); }

void SincTime(int argc, char* argv[], int& qn, SincPulDat& SD)
  { query_parameter(argc,argv,qn++,"\n\tSinc Pulse Length (msec)? ", SD.tau);
    SD.tau *= 1.e-3;                                                          }

void SincGamB1(int argc, char* argv[], int& qn, SincPulDat& SD)
  { query_parameter(argc,argv,qn++,"\n\tRF-Field Stength (Hz)? ", SD.gamB1);  }

void SincAngle(int argc, char* argv[], int& qn, double& pang)
  { query_parameter(argc,argv,qn++,"\n\tSinc Pulse Angle (degrees)? ", pang); }

void SincIso(int argc, char* argv[], int& qn, SincPulDat& SD)
  { query_parameter(argc,argv,qn++,"\n\tPulse Channel (1H,13C..)? ", SD.Iso); }

void SincWrf(int argc, char* argv[], int& qn, SincPulDat& SD)
  { query_parameter(argc,argv,qn++,"\n\tSinc Pulse Offset (Hz)?  ", SD.Wrf);  }

void SincPhi(int argc, char* argv[], int& qn, SincPulDat& SD)
  { query_parameter(argc,argv,qn++,"\n\tSinc Pulse Phase (degrees)? ",SD.phi);}


SincPulDat SincAsk(int argc, char* argv[], int& qn, int type)

        // Input        argc    : Number of argments
        //              argv    : Array of arguments
        //              qn      : Argument index
        //              type    : How to characterize the Sinc
        //                         0 = angle and length (default)
        //                         1 = strength and length
        //                         2 = strength and time 
        // Output       SincData: The sinc pulse structer values set
        //                        interactively or by values in argv.

  {
  SincPulDat SD;				// Sinc parameters
  SincZero(SD);
  SincPts(argc, argv, qn, SD);			// 1.) Get number of steps
  SincNode(argc,argv, qn, SD);			// 2.) Get end node
  double pang;					// For pulse angle
  switch(type)					// 3.) Get strength and
    {						// 4.) and length 
    case 0:					// Sinc by angle & time
    default:
      SincTime(argc, argv, qn, SD);		//   Ask for/Set time
      SincAngle(argc, argv, qn, pang);		//   Ask for/get angle
      SD.gamB1 =				//   Set pulse rf
           SincGamB1(pang,SD.tau,SD.N,SD.node);
      break;
    case 1:					// Sinc by gamB1 & time
      SincGamB1(argc, argv, qn, SD);		//   Ask for/Set gamB1
      SincTime(argc, argv, qn, SD);		//   Ask for/Set time
      break;
    case 2:					// Sinc by angle & gamB1
      SincGamB1(argc, argv, qn, SD);		//   Ask for/Set gamB1
      SincAngle(argc, argv, qn, pang);		//   Ask for/get angle
      SD.tau =					//   Set pulse length
           SincTime(pang,SD.gamB1,SD.N,SD.node);
      break;
    }
  SincIso(argc,argv, qn, SD);			// 5.) Get pulse channel
  SincWrf(argc,argv, qn, SD);			// 6.) Get pulse offset
  SincPhi(argc,argv, qn, SD);			// 7.) Get pulse phase
  return SD;
  } 

// sosi - this function makes us depende upon the Level2 module
 
double ask_Sinc(spin_system& sys, std::string& Iso, gen_op& H, double cutoff)

	// Input		sys	: A spin system
	//			Iso     : String designating an isotope
	//			H	: Isotropic Hamiltonian (Hz)
	//			cutoff	: An intensity cutoff
	// Output		v       : A transition of H (Hz)
	// Note				: This routine looks over all single

  {
  gen_op detect = Fm(sys, Iso);		// Set detection operator to F-
  gen_op sigma = Fx(sys);		// Set system in transverse plane
  super_op L = complexi*Hsuper(H);	// L = -i*[Ho, ] (rad/sec)
  acquire1D ACQ(detect, L);		// Prepare for "Dirac" acquisitions
  TTable1D trans = ACQ.table(sigma);	// Determine "Dirac" Acquisition
  int ntr = trans.size();		// Get the number of transitions
  int itr = 0;
  double v;
  cout << "\n\n\tThere are " << ntr << " " << Iso << " Transitions\n";
  for(int i=0; i<ntr; i++)
    cout << "\n\t\t" << i << ". " << trans.Fr(i)/(PIx2);
  cout << "\n\n\tChoose [0-" << ntr-1 << "] To Set Sinc Repetition,"
       << " Any Other Integer To Specify A Different Value: ";
  cin >> itr;
  if(itr >= ntr || itr < 0)
    {
    cout << "\n\t\t" << "Please Input A Value: ";
    cin >> v;
    }
  else
    v = trans.Fr(itr)/(PIx2); 
  cout << "\n";
  return v;
  }


void set_Sinc(double gamB1, double& tmix, double& tpul, double tdel, int& numb, int& type)

	// Input		gamB1: RF-Field strength (Hz)
	// 			tmix : Total Sinc mixing time
	// 			tpul : Individual pulse length
	// 			tdel : Individual delay length
	// 			numb : Times Sinc needs repeating for tmix
	// 			type : Flag for Sinc type
	//				0 = no pulse, delay with no relaxation
	//				1 = ideal pulse, delay with no relaxation
	//				2 = real pulse, delay with no relaxation
	//				3 = relaxation pulse, delay with no relaxation
	//				11 = ideal pulse, delay with relaxation
	//				12 = real pulse, delay with relaxation
	//				13 = relaxation pulse, delay with relaxation
	// Output		     : Void, argument parameters are set
	// Note			     : gamB1 should be set to zero on input if
	// 			       no external value is present

    {
    cout << "\n\n\t\tSinc MIXING SEQUENCE SETUP\n";

    cout << "\n\tTotal Mixing Time Desired? ";	// Get the total mixing time
    cin >> tmix;

    double tangle;
    cout << "\n\tTotal Pulse Rotation Angle? ";	// Get the total pulse angle
    cin >> tangle;

    cout << "\n\tNumber of Pulse-Delay Steps?";	// Get the number of steps
    cin >> numb;

    double angle = tangle/numb;
    tpul = 0.0;
    int repeat = 1;
    std::string typ;
    type = 0;
    if(angle)
      {
      cout << "\n\tThe Rotation Angle of Each Pulse is " << angle << " Degrees"; 
      if(gamB1)
        {
        while(repeat)
          {
          cout << "\n\tPulse Type Desired [i=ideal pulses, p=pulse, r=relaxation]? ";
          cin >> typ;
          if(typ == "i")
            {
            cout << "\n\tPulses Set Ideal";
            repeat = 0;
            type = 1;
            }
          if(typ == "p")
            {
            tpul = angle/(360.0*gamB1);
            cout << "\n\tEach Pulse Lasts " << tpul << " Seconds at This Field Strength";
            repeat = 0;
            type = 2;
            }
          if(typ == "r")
            {
            tpul = angle/(360.0*gamB1);
            cout << "\n\tEach Pulse Lasts " << tpul << " Seconds at This Field Strength";
            repeat = 0;
            type = 3;
            }
          }
        }
      else
        {
        cout << "\n\tPulses Set Ideal";
        type = 1;
        }
      }
  
    double ttot = tmix - double(numb)*tpul;
    tdel = ttot/double(numb);
    if(tdel)
      {
      cout << "\n\tEach Delay Step Lasts " << tdel << " Seconds";
      repeat = 1;				// Get the type of Sinc
      while (repeat)
        {
        cout << "\n\n\tInclude Relaxation During Delays "
             << "[y = yes, relaxation," " n=no, real pulses]? ";
        cin >> typ;
        if(typ == "y")
          {
          repeat = 0;
          type += 10;				// Relaxation Desired
          }
        else if(typ == "n")			// No Relaxation Desired
          {
          repeat = 0;
          }
        }
      }
    return;
    }

// ____________________________________________________________________________
//                       Sinc Pulse Output Functions
// ____________________________________________________________________________

 
void SincPrint(std::ostream& ostr, const SincPulDat& SincData, int indx)

        // Input        ostr		: An output stream
	//		SincData	: Sinc pulse structure
        // Output       void		: Sinc pulse parameters
	//			  	  are put into ostr

    {
    int il=0;
    std::string ival, blanks = "         ";
    ostr << "\n\n\t\tSinc Pulse ";
    if(indx>=0) ostr << indx << " ";
    ostr << "Parameters\n";
    ostr << "\n\tSinc Pulse Frequency (Hz): " << Gform("%11.3f", SincData.Wrf);
    ostr << "\n\tSinc Pulse Length (sec):   " << Gform("%11.3f", SincData.tau);
    ostr << "\n\tSinc Pulse Angle(degrees): "
         << Gform("%11.3f", SincAngle(SincData.gamB1, SincData.tau, SincData.N, SincData.node));
    ostr << "\n\tSinc Pulse Phase (degrees):" << Gform("%11.3f", SincData.phi);
    ival = Gdec(SincData.N);
    il = 6-ival.length();
    ostr << "\n\tNumber of Sinc Pulse Steps: " << std::string(il, ' ') << SincData.N;
    il = 6-(SincData.Iso).length();
    ostr << "\n\tSinc Pulse Selectivity:     " << std::string(il, ' ') << SincData.Iso;
    ostr << "\n\tSinc Pulse Strength (Hz):  " << Gform("%11.3f", SincData.gamB1);
    ostr << "\n\tSinc Pulse End Node:       " << Gform("%7i", SincData.node);
    ostr << "\n";
    }


void SincZero(SincPulDat& SD)
  {
  SD.N     = 0;					// Set pulse steps
  SD.Wrf   = 0.0;                                // Set pulse offset
  SD.Iso   = std::string("1H");                       // Set pulse selectivity
  SD.tau   = 0;					// Set the pulse length
  SD.node  = 0;                                  // Set pulse node
  SD.phi   = 0.0;                                // Set pulse phase
  SD.gamB1 = 0.0;				// Set pulse strength
  }


void SincPrep(int& N, int& node, double& den, double& Z)

   {
   if(N<1)    N = 1;				// Insure at least 1 point
   if(node<1) node = 3;				// Make sure node is between
   else if(node > 30) node = 3;			// [1,30].  More is silly.

   den = double((N-1)*(N-1));			// Calculate the denominator 
   Z = 0.0;					// Exponential nodeor, set to
   if(N > 1) Z = log(double(node))/den;		// zero if N=1;
   return;
   }



#endif                                             		// PulSinc.cc

