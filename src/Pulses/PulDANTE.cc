/* PulDANTE.cc *************************************************-*-c++-*-*
**                                                                      **
**                                   G A M M A                          **
**                                                                      **
**      DANTE Pulse Functions                        Implementation	**
**                                                                      **
**      Copyright (c) 1998                                              **
**      Scott Smith                                                     **
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** This GAMMA module facilitates the use of DANTE pulse trains in NMR   **
** simulations.                                                         **
**                                                                      **
*************************************************************************/

#ifndef _DANTE_cc_			// Is this file already included?
#define _DANTE_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Pulses/PulDANTE.h>		// Include the interface
#include <Pulses/PulAuxil.h>		// Include auxiliary functions
#include <Pulses/PulWaveform.h>		// Include pulse waveforms
#include <Pulses/PulComposite.h>	// Include composite pulses
#include <Pulses/PulTrain.h>		// Include pulse trains
#include <HSLib/SpinOpCmp.h>		// Include composite spin operators
#include <HSLib/SpinOpRot.h>		// Include spin rotation operators
#include <HSLib/SpinSystem.h>		// Include istropic spin systems
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets
#include <Level2/acquire1D.h>		// Include 1D acquisitions
#include <Basics/StringCut.h>		// Include Gdec and Gform
#include <stdlib.h>
#include <string>                       // Include libstdc++ strings
#include <iostream>                     // Include input output streams (cout)

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams
using std::cout;			// Using libstdc++ standard output
using std::cin;				// Using libstdc++ standard input

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS DANTE ERROR HANDLING
// ____________________________________________________________________________


void DANTE::DANTEerror(int eidx, int noret) const

        // Input                DANTE	: DANTE parameters
        //                      eidx    : Error flag
        //                      noret   : Return flag
        // Output               none    : Error Message Output

  {
  cout << "\n\tDANTE Parameters: ";
  switch(eidx)
    {
    case 0:                                                             // (0)
      cout << "Program Aborting....";
      break;
    case 1:                                                             // (1)
      cout << "Error During Construction";
      break;
    case 2:								// (2)
      cout << "Cannot Find An Adequate Set of Pulse Parameters";
      break;
    case 3:								// (3)
      cout << "Cannot Have No Delay With No/Ideal Pulse";
      break;
    default:
      cout << "Unknown Error (Number " << eidx << ")";
    }
  if(!noret) cout << ".\n";
  else       cout << ".";
  }

    
void volatile DANTE::DANTEfatality(int eidx) const

        // Input                DANTE	: DANTE parameters
        //                      eidx    : Error flag
        // Output               none    : Stops execution & error Message

  {                                                                       
  DANTEerror(eidx,1);
  if(eidx) DANTEerror(0);
  cout << "\n";
  exit(-1);
  }


 
void DANTE::DANTEerror(int eidx, const string& pname, int noret) const
 
        // Input                DANTE	: DANTE parameters
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
 
  {
  cout << "\nDANTE Parameters: ";
  switch(eidx)
    {
    case 40:                                                    // (40)
      cout << "Problems with File " << pname;
      break;
    case 100:                                                   // (100)
      cout << "Can't Read Parameter " << pname;
      break;
    case 101:                                                   // (101)
      cout << "Can't Find DANTE Parameters For " << pname;
      break;
    case 130:                                                   // (130)
      cout << "Parameter " << pname << " Is The Culprit!\n";
      break;
    default:
      cout << "Unknown error";
      break;
    }
  if(!noret) cout << ".\n";
  }  

// ____________________________________________________________________________
// ii                CLASS DANTE PARAMETER SET FUNCTIONS
// ____________________________________________________________________________


void DANTE::SetSteps(const ParameterSet& pset, int idx)
 
        // Intput		DT	: DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               none    : DANTE number of steps read in
        //                                from parameter in pset
        //                                with index idx

  {
  int npt;					// Parameter variable
  string pstate;				// Dummy string variable
  string pname = string("DANTEstps");		// Number of DANTE steps
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  std::list<SinglePar>::const_iterator item;	// A pix into parameter list
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,npt,pstate);	 	// 	Get DANTE steps
    N = npt;					// 	Set DANTE steps
    }
  else N = 0;					// If not found set to 0
  }


void DANTE::SetPhase(const ParameterSet& pset, int idx)
 
        // Intput		DT	: DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               none    : DANTE pulse phase read in
        //                                from parameter in pset
        //                                with index idx

  {
  double phiin;
  string pstate;				// Dummy string variable
  string pname = string("DANTEphi");		// DANTE phase angle
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,phiin,pstate); 	// 	Get DANTE phase
    phi= phiin;					// 	Set DANTE phase
    }
  else phi = 0;					// If not found set to 0
  }


void DANTE::SetChannel(const ParameterSet& pset, int idx)
 
        // Intput		DT	: DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               none    : DANTE pulse channel read in
        //                                from parameter in pset
        //                                with index idx

  {
  string Diso;					// Variable for reading
  string pstate;				// Dummy string variable
  string pname = string("DANTEiso");		// Channel selectivity
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,Diso,pstate);	// 	Read selectivity
    Iso = Diso;					// 	Set DANTE channel
    }
  else Iso = "";				// If not found, don't set
  }


int DANTE::SetAngle(const ParameterSet& pset, int idx)
 
        // Intput		DT	: DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               TF	: DANTE pulse angle read in
        //                                from parameter in pset
        //                                with index idx

  {
  double pang;					// Variable for reading
  string pstate;				// Dummy string variable
  string pname = string("DANTEang");		// DANTE pulse angle
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,pang,pstate); 	// 	Get DANTE pulse angle
    ang = pang;					// 	Set DANTE pulse angle
    return 1;					//	Return TRUE
    }
  else ang = 0;					// If not found set to 0
  return 0;					// Return FALSE
  }


int DANTE::SetPulLen(const ParameterSet& pset, int idx)
 
        // Intput		DT	: DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               TF	: DANTE pulse length read in
        //                                from parameter in pset
        //                                with index idx

  {
  double plen;					// Variable for reading
  string pstate;				// Dummy string variable
  string pname = string("DANTEtp");		// DANTE pulse length
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,plen,pstate); 	// 	Get DANTE pulse length
    tp = plen;					// 	Set DANTE pulse length
    return 1;					//	Return TRUE
    }
  else tp = 0;					// If not found set to 0
  return 0;					// Return FALSE
  }


int DANTE::SetGamB1(const ParameterSet& pset, int idx)
 
        // Intput		DT	: DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               TF	: DANTE pulse strength read in
        //                                from parameter in pset
        //                                with index idx

  {
  double gB1;					// Variable for reading
  string pstate;				// Dummy string variable
  string pname = string("DANTEgamB1");		// DANTE pulse strength
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,gB1,pstate); 	// 	Get DANTE pul strength
    gamB1 = gB1;				// 	Set DANTE pul strength
    return 1;					//	Return TRUE
    }
  else gB1 = 0;					// If not found set to 0
  return 0;					// Return FALSE
  }


int DANTE::SetEvLen(const ParameterSet& pset, int idx)
 
        // Intput		DT	: DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               TF	: DANTE delay length read in
        //                                from parameter in pset
        //                                with index idx

  {
  double dlen;					// Variable for reading
  string pstate;				// Dummy string variable
  string pname = string("DANTEtd");		// DANTE pulse length
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,dlen,pstate); 	// 	Get DANTE delay length
    te = dlen;					// 	Set DANTE delay length
    return 1;					//	Return TRUE
    }
  else te = 0;					// If not found set to 0
  return 0;					// Return FALSE
  }


int DANTE::SetFreq(const ParameterSet& pset, int idx)
 
        // Intput		DT	: DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               TF	: DANTE frequency read in
        //                                from parameter in pset
        //                                with index idx

  {
  double Fr;					// Variable for reading
  string pstate;				// Dummy string variable
  string pname = string("DANTEF");		// DANTE frequency
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,Fr,pstate); 	// 	Get DANTE frequency
    F = Fr;					//	Set DANTE frequency
    double cycs = 1/F;				//	Full cycles of F
    te = cycs - tp;				//	Delay time for 1 cycle
    while(te < 0) te += cycs;			//	Search for evolve time 
    tt = tp + te;				//	Set DANTE step time
    return 1;					//	Return TRUE
    }
  te = 0;					// If not found set to 0
  F = 1/tt;					// This is the frequency
  return 0;					// Return FALSE
  }

 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A                  CLASS DANTE CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 
 
DANTE::DANTE()
 
        // Input        none    :
        // Output	DANTE	: DANTE parameters (this)
 
  {
  N         = 0;			// No pulse waveform steps
  te	    = 0;			// No evolution time
  tp	    = 0;			// No pulse time
  gamB1     = 0;			// No rf-field strength	
  phi       = 0;			// No rf-field phase
  Wrf       = 0;			// No rf-field offset
  tt        = 0;			// No step time
  ang       = 0;			// No pulse angle
  F         = 0;			// No synch frequency
  }
 
                                                                                
DANTE::DANTE(const DANTE& DANTE1)
 
        // Intput       DANTE1  : DANTE parameters
        // Output       DANTE   : DANTE parameters(this) from DANTE1

  {
  N     = DANTE1.N;			// Copy pulse waveform steps
  te	= DANTE1.te;			// Copy evolution time
  tp	= DANTE1.tp;			// Copy pulse time
  gamB1 = DANTE1.gamB1;			// Copy rf-field strength	
  phi   = DANTE1.phi;			// Copy rf-field phase
  Wrf   = DANTE1.Wrf;			// Copy rf-field offset
  tt	= DANTE1.tt;			// Copy step length
  ang	= DANTE1.ang;			// Copy pulse angle
  F	= DANTE1.F;			// Copy synch frequency
  }
 
                                                                              
// ------------------------------ Destruction ---------------------------------
 
                                                                                
DANTE::~DANTE() {}
 
        // Intput       DANTE1  : DANTE parameters (this)
        // Output       none    : DANTE is destructed
 
                                                    
// ------------------------------- Assignment ---------------------------------
 
                                         
DANTE& DANTE::operator = (const DANTE& DANTE1)

        // Intput       DANTE1  : DANTE parameters
        // Output       DANTE   : DANTE parameters(this) from DANTE1
        ///F_list =           - Assignment

  {
  N     = DANTE1.N;			// Copy pulse waveform steps
  te	= DANTE1.te;			// Copy evolution time
  tp	= DANTE1.tp;			// Copy pulse time
  gamB1 = DANTE1.gamB1;			// Copy rf-field strength	
  phi   = DANTE1.phi;			// Copy rf-field phase
  Wrf   = DANTE1.Wrf;			// Copy rf-field offset
  tt	= DANTE1.tt;			// Copy step length
  ang	= DANTE1.ang;			// Copy pulse angle
  F	= DANTE1.F;			// Copy synch frequency

  return (*this);
  }
 
 
// ____________________________________________________________________________
// B                     CLASS DANTE ACCESS FUNCTIONS
// ____________________________________________________________________________
 

int        DANTE::steps()    const { return N; }
string     DANTE::channel()  const { return Iso; }
double     DANTE::dlength()  const { return te; }
double     DANTE::strength() const { return gamB1; }
double     DANTE::plength()  const { return tp; }
double     DANTE::angle()    const { return ang; }
double     DANTE::phase()    const { return phi; }
double     DANTE::offset()   const { return Wrf; }
double     DANTE::length()   const { return te+tp; }

        // Intput       DANTE1  : DANTE parameters
        // Output       steps   : DANTE steps
	//		channel : DANTE isotope channel
	//		dlength	: DANTE delay length (sec)
	//		strength: DANTE pulse strength (Hz)
	//		dlength	: DANTE pulse length (sec)
	//		angle   : DANTE pulse angle (sec)
	//		phase   : DANTE pulse phase (sec)
	//		offset  : DANTE pulse offset (Hz)

 
// ____________________________________________________________________________
// C                  CLASS DANTE HILBERT SPACE PROPAGATORS
// ____________________________________________________________________________
 

/*
HSprop DANTE::GetU(const spin_system& sys, gen_op& H)

	// Input	     sys   : Spin system
	// 		     H     : Static Hamlitonian without the field
        //		     D	   : DANTE specifications
	// Output	     U     : Propagator for a DANTE step
	// Note			   : The propagator must be built in REVERSE
	//			     order for the pulse preceeds the delay

  {
  HSprop Udante(H, te);			// Delay propagator
  if(tp) Udante *= Pulse::GetU(sys,H); 	// Pulse propagator
  return Udante;
  }
*/

 
// ____________________________________________________________________________
// D                CLASS DANTE LIOUVILLE SPACE PROPAGATORS
// ____________________________________________________________________________
 

// ____________________________________________________________________________
// E                       DANTE WAVEFORM FUNCTIONS
// ____________________________________________________________________________

 
PulWaveform WF_DANTE(const DANTE& D) { return D.WF(); }
PulWaveform DANTE::WF( ) const
 
        // Input        D       : DANTE specifications
        // Output       WF      : A pulse waveform for DANTE pulse-delay 
	// Note			: Ideal pulses are indicated by placing
	//			  the pulse angle in WFsteps (rather than
	//			  the pulse strength) and setting the step
	//			  length to zero.

  {
  row_vector WFsteps(2);                        // Waveform vector{ gamB1,phi }
  if(!tp) WFsteps.put(ang,   0);		// Ideal pulse, put in angle! 	
  else    WFsteps.put(gamB1, 0);		// Set pulse strength
  WFsteps.put(0, 1);				// Set delay step
  row_vector WFtimes(2);                        // Vector for step times
  WFtimes.put(tp, 0);
  WFtimes.put(te, 1);
  return PulWaveform(WFsteps, WFtimes, "DANTE Pulse-Delay");
  }


// ____________________________________________________________________________
// F                    DANTE COMPOSITE PULSE FUNCTIONS
// ____________________________________________________________________________


PulComposite CP_DANTE(const spin_system& sys, const DANTE& D)

        // Input        sys     : Spin system                                  
        // 		D       : DANTE specifications
        // Output       CP      : A DANTE composite pulse

 { return D.CP(sys); }


PulComposite DANTE::CP(const spin_system& sys) const

        // Input        D       : DANTE specifications(this)
        // Input        sys     : Spin system                                  
        // Output       CP      : A DANTE composite pulse

  { return PulComposite(WF(), sys, Iso); }


// ____________________________________________________________________________
// G                       DANTE PULSE TRAIN FUNCTIONS
// ____________________________________________________________________________
 

PulTrain PT_DANTE(const spin_system& sys, const DANTE& D)

        // Input        sys     : Spin system                                  
        // 		D       : DANTE specifications
        // Output       CP      : A DANTE composite pulse

  { return D.PT(sys); }


PulTrain DANTE::PT(const spin_system& sys) const

        // Input        sys     : Spin system                                  
        // 		D       : DANTE specifications
        // Output       CP      : A DANTE composite pulse

  { return PulTrain(CP(sys), "DANTE"); }




// ____________________________________________________________________________
// Y                      CLASS DANTE INPUT FUNCTIONS
// ____________________________________________________________________________


void DANTE::read(const string &filename, int idx)

        // Intput		DT	: DANTE parameters
        //                      idx     : DANTE index
        // Output               none    : DANTE parameters are read in
        //                                from parameters in file filename
        //                                with index idx
              
  {
  ParameterSet pset;                // Declare a parameter set
  if(!pset.read(filename, 1))          // Read in pset from file
    {   
    DANTEerror(40, filename);		// Filename problems
    DANTEfatality(21);			// Fatal error
    }   
  read(pset, idx);			// User overloaded function
  return;
  }
 

void DANTE::read(const ParameterSet& pset, int idx)
 
        // Intput		DT	: DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               none    : DANTE parameters are read in
        //                                from parameters in pset
        //                                with index idx
 
  {
//	        First We'll Set Up The DANTE Pulse Parameters
// For The Pulse We Need Two of Three: {Pulse Angle, RF Strength, Pulse Length} 

// Note: Reading order is angle --> length --> strength for DANTE.  If only
//       the angle is specified the length will be set to zero (ideal pulse)

  int G=0, T=0, A=0;				// Flags if pulse params found
  A = SetAngle(pset, idx);			// Try to set pulse angle
  if(A)						// If pulse angle specified
    {
    T = SetPulLen(pset, idx);			// Try for pulse length
    if(T) gamB1 = ang/(360.0*tp);		// If found set pulse strength
    else					// If not try for strength
      {
      G = SetGamB1(pset, idx);			// If strength set, set pulse
      if(G) tp = ang/(360.0*gamB1);		// length from {ang, gamB1} 
      }
    }
  else						// If no angle has been set
    {						// we must demand {tp, gamB1}
    T = SetPulLen(pset, idx);			// Set pulse length
    G = SetGamB1(pset, idx);			// Set pulse strength
    if(!T || !G) DANTEfatality(2); 		// Quit if either isn't set
    }
  SetChannel(pset,idx);				// Set DANTE pulse channel
  SetPhase(pset, idx);				// Set DANTE pulse phase

//	        Next We Set Up The DANTE Delay Parameters

  if(!SetFreq(pset, idx))			// Set delay frequency based
    if(!SetEvLen(pset, idx))			// Else just set delay
      if(!tp) DANTEfatality(3);			// No/ideal pulse & no delay!

//	        Last We Set Up Any Other DANTE Parameters

  SetSteps(pset, idx);				// Set # of DANTE steps
  }


void DANTE::ask_read(int argc, char* argv[], int argn)

        // Intput		DT	: DANTE parameters
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the DANTE parameters
	//				  are read
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : DANTE parameters are modifed (filled)

  {
  string filename;                            // Name of spin system file
  query_parameter(argc, argv, argn,           // Get filename from command
    "\n\tDANTE parameters filename? ", filename); // Or ask for it
  read(filename);                             // Read system from filename
  }



// ____________________________________________________________________________
// Z                      CLASS DANTE I/O FUNCTIONS
// ____________________________________________________________________________
 
 
ostream& DANTE::printBase(ostream &ostr) const
 
        // Intput		DT	: DANTE parameters
        //                      ostr    : Output stream
        // Output               none    : DANTE basic parameters are sent
        //                                to the output stream
 
  {
  ostr << "\n\tIsotope channel:           " << Iso; 
  ostr << "\n\tDelay time:              ";
  printTime(ostr, te);
  ostr << "\n\tPulse length:            ";
  printTime(ostr, tp);
  ostr << "\n\tPulse strength:          ";
  printHz(ostr, gamB1);
  ostr << "\n\tPulse angle:             ";
  ostr << Gform("%8.3f", ang) << " degrees"; 
  ostr << "\n\tPulse phase:             ";
  ostr << Gform("%8.3f", phi) << " degrees"; 
  ostr << "\n\tPulse offset:            ";
  printHz(ostr, Wrf);
  return ostr;
  }
 
 
ostream& DANTE::printInfo(ostream &ostr) const
 
        // Intput		DT	: DANTE parameters
        //                      ostr    : Output stream
        // Output               none    : DANTE additional information sent
        //                                to the output stream
 
  {
  ostr << "\n\tDANTE step length:       ";
  printTime(ostr, tt); 
  ostr << "\n\tDANTE frequency:         ";
  printHz(ostr, F); 
  return ostr;
  }

 
ostream& DANTE::print(ostream &ostr, int full) const
 
        // Intput		DT	: DANTE parameters
        //                      ostr    : Output stream
        //                      full    : Flag for output amount
        // Output               none    : DANTE parameters are sent
        //                                to the output stream
 
  {
  ostr << "\n\n\t\t\tDANTE Parameters\n";
  printBase(ostr);
  if(full) printInfo(ostr);
  ostr << "\n\tDANTE steps:         " << Gdec(N,8); 
  ostr <<"\n\n";
  return ostr;
  }
 

ostream &operator << (ostream &ostr, const DANTE &DT)
 
        // Intput		DT	: DANTE parameters
        //                      ostr    : Output stream
        // Output               none    : DANTE parameters are sent
        //                                to the output stream
 
  {
  DT.print(ostr);
  return ostr;
  }


// ____________________________________________________________________________
// ____________________________________________________________________________
//                       END OF DANTE CLASS FUNCTIONS
// ____________________________________________________________________________
// ____________________________________________________________________________

// ____________________________________________________________________________
// A                       DANTE WAVEFORM FUNCTIONS
// ____________________________________________________________________________


PulWaveform WF_DANTE(double td, double gamB1, double tpul, double phi)

	// Input	td    : Delay following pulse application
	// 		gamB1 : RF-Field strength (Hz)
	// 		tpul  : Pulse length
	// 		phi   : Pulse phase (degrees)
        // Output	WF      : A pulse waveform for DANTE pulse-delay

  {
  row_vector WFsteps(2);                        // Waveform vector{ gamB1,phi }
  WFsteps.put(complex(gamB1,0.0),    0);        // Set waveform values
  WFsteps.put(complex(0.0, 0.0),     1);
  row_vector WFtimes(2);                        // Vector for step times
  WFtimes.put(tpul,  0);
  WFtimes.put(td,    1);
  return PulWaveform(WFsteps, WFtimes, "DANTE Pulse-Delay");
  }
 
// ____________________________________________________________________________
// B                     DANTE COMPOSITE PULSE FUNCTIONS
// ____________________________________________________________________________

 
PulComposite CP_DANTE(const spin_system& sys, const  string& Iso,
                             double td, double gamB1, double tpul, double phi)
 
        // Input        sys     : Spin system
        //              Iso     : Isotope channel the field is on
        //              td      : Delay following pulse application 
        //              gamB1   : RF-Field strength (Hz) 
        //              tpul    : Pulse length 
        //              phi     : Pulse phase (degrees)
        // Output       CP      : A DANTE composite pulse
 

  { return PulComposite(WF_DANTE(td,gamB1,tpul,phi), sys, Iso); }
 

// ____________________________________________________________________________
// C                       DANTE PULSE TRAIN FUNCTIONS
// ____________________________________________________________________________
 
 
PulTrain PT_DANTE(const spin_system& sys, const string& Iso,
                              double td, double gamB1, double tpul, double phi)
 
        // Input        sys     : Spin system
        //              Iso     : Isotope channel the field is on
        //              td      : Delay following pulse application 
        //              gamB1   : RF-Field strength (Hz) 
        //              tpul    : Pulse length 
        //              phi     : Pulse phase (degrees)
        // Output       PT      : A DANTE pulse train
	// Note			: Since DANTE is applied without any cycle
	//			  or supercycle, the DANTE composite pulse
	//			  should be functionally equivalent to this

  { return PulTrain(CP_DANTE(sys, Iso, td, gamB1, tpul, phi), "DANTE"); }

 

//PulComposite PT_DANTE(const spin_system& sys, const DANTE& D)

        // Input        sys     : Spin system                                  
        // 		D       : DANTE specifications
        // Output       CP      : A DANTE composite pulse

//  { return PulTrain(CP_DANTE(sys, D), "DANTE"); }


// ____________________________________________________________________________
// D                        DANTE PULSE PROPAGATORS
// ____________________________________________________________________________

// -------------------- Individual DANTE Pulse-Delay Steps --------------------

gen_op UDANTE(const spin_system& sys, gen_op& H, const string& Iso,
                                           double td, double theta, double phi)

	// Input	     sys   : Spin system
	// 		     H     : Static Hamlitonian without the field
	// 		     Iso   : Isotope channel the field is on
	// 		     td    : Delay following pulse application
	// 		     theta : Pulse angle (degrees)
	// 		     phi   : Pulse phase (degrees)
	// Output	     U     : Propagator for a DANTE sequence
	//			     about the axis specified by phi
	// Note			   : Ideal Pulses, No Relaxation

  {
  gen_op Udante = prop(H, td);			// Delay propagator (no relaxation)
  Udante *= Ixypuls_U(sys, Iso, phi, theta);	// Multiply into pulse propagator (ideal)
  return Udante;
  }


gen_op UDANTE(const spin_system& sys, gen_op& H, const string& Iso,
                              double td, double gamB1, double tpul, double phi)

	// Input	     sys   : Spin system
	// 		     H     : Static Hamlitonian without the field
	// 		     Iso   : Isotope channel the field is on
	// 		     td    : Delay following pulse application
	// 		     gamB1 : RF-Field strength (Hz)
	// 		     tpul  : Pulse length
	// 		     phi   : Pulse phase (degrees)
	// Output	     U     : Propagator for an DANTE sequence
	//			     about the axis specified by phi
	// Note			   : Real Pulses, No Relaxation

    {
    gen_op Udante = prop(H, td);		// Delay propagator (no relaxation)
    Udante *= 					// Multiply into pulse propagator (real)
       SxypulsB_U(sys,H,Iso,0.0,tpul,gamB1,phi);  
    return Udante;
    }


gen_op UDANTE(const spin_system& sys, gen_op& H, const DANTE& D)

	// Input	     sys   : Spin system
	// 		     H     : Static Hamlitonian without the field
        //		     D	   : DANTE specifications
	// Output	     U     : Propagator for a DANTE step

  {
  gen_op Udante = prop(H, D.dlength());		// Delay propagator (no relaxation)
  if(D.plength())
    Udante *= 					// Multiply into pulse propagator (real)
       Sxypuls_U(sys,H,D.channel(),D.offset(),D.plength(),D.angle(),D.phase());  
  else
    Udante *= 					// Multiply into pulse propagator (real)
       Ixypuls_U(sys,D.channel(),D.phase(),D.angle());  
  return Udante;
  }

/*

super_op GDANTE(spin_system& sys, gen_op& H, string& Iso,
                                   super_op& L, gen_op& sigma0,
                                         double td, double theta, double phi)

	// Input	     sys   : Spin system
	// 		     H     : Static Hamlitonian without the field
	// 		     Iso   : Isotope channel the field is on
	//		     L     : Relaxation Liouvillian
	//		     sigma0: Infinite time density matrix
	// 		     td    : Delay following pulse application
	// 		     gamB1 : RF-Field strength (Hz)
	// 		     tpul  : Pulse length
	// 		     phi   : Pulse phase (degrees)
	// Output	     U     : Propagator for an DANTE sequence
	//			     about the axis specified by phi
	// Note			   : Ideal Pulses, Relaxation in Delay

    {
    super_op eLt = exp(L, -td);			// Exponential for delay with relaxation
    super_op Dante = R_prop(eLt, sigma0);	// Propagator for delay with relaxation
    gen_op Up = Ixypuls_U(sys,Iso,phi,theta);	// Propagator for ideal pulse
    Dante.LOp_base(Up);				// Insure Up in right basis
    Dante *= U_transform(Up);			// Multiply into pulse propagator (real)
    return Dante;
    }


super_op GDANTE(spin_system& sys, gen_op& H, string& Iso,
                                   super_op& L, gen_op& sigma0,
                          double td, double gamB1, double tpul, double phi)

	// Input	     sys   : Spin system
	// 		     H     : Static Hamlitonian without the field
	// 		     Iso   : Isotope channel the field is on
	//		     L     : Relaxation Liouvillian
	//		     sigma0: Infinite time density matrix
	// 		     td    : Delay following pulse application
	// 		     gamB1 : RF-Field strength (Hz)
	// 		     tpul  : Pulse length
	// 		     phi   : Pulse phase (degrees)
	// Output	     U     : Propagator for an DANTE sequence
	//			     about the axis specified by phi
	// Note			   : Real Pulses, Relaxation in Delay

    {
    super_op eLt = exp(L, -td);			// Exponential for delay with relaxation
    super_op Dante = R_prop(eLt, sigma0);	// Propagator for delay with relaxation
    gen_op Up =                           	// Propagator for real pulse
       SxypulsB_U(sys,H,Iso,0.0,tpul,gamB1,phi);  
    Dante.LOp_base(Up);				// Insure Up in right basis
    Dante *= U_transform(Up);			// Multiply into pulse propagator (real)
    return Dante;
    }


super_op GDANTE(spin_system& sys, super_op& L, gen_op& sigma0, double td,
                               super_op& Lrf, gen_op& sigmass, double tp, double phi)

	// Input	   sys     : Spin system
	//		   L       : Relaxation Liouvillian, no field
	//		   sigma0  : Infinite time density matrix, no field
	// 		   td      : Delay following pulse application
	//		   Lrf     : Relaxation Liouvillian, with field
	//		   sigmass : Infinite time density matrix, with field
	// 		   tp      : Pulse length
	// 		   phi     : Pulse phase (degrees)
	// Output	     U     : Propagator for an DANTE sequence
	//			     about the axis specified by phi
	// Note			   : Relaxation During Both Pulse & Delay

  {
  super_op eLt = exp(L, -td);			// Exponential for delay with relaxation
  super_op Dante = R_prop(eLt, sigma0);	// Propagator for delay with relaxation
  eLt = exp(Lrf, -tp);			// Exponential for field with relaxation
  Dante *= R_prop(eLt, sigmass);		// Multiply into delay with field & relaxation
  return Dante;
  }


// ---------------------- Total DANTE Pulse Propagators -----------------------


gen_op DANTE(spin_system& sys, gen_op& H, string& Iso,
                      double v, int n, double theta, double phi, int rad)

	// Input	     sys   : Spin system
	// 		     H     : Static Hamlitonian without the field
	// 		     Iso   : Isotope channel the field is on
	// 		     v     : Dante "frequency" center (Hz)
	// 		     n     : Number of pulse-delay steps
	// 		     theta : Total DANTE Pulse angle (degrees)
	// 		     phi   : Pulse phase (degrees)
	//		     rad   : Flag if frequency input in radians/sec
	// Output	     U     : Propagator for an DANTE sequence
	//			     about the axis specified by phi
	// Note			   : Ideal Pulses, No Relaxation

    {
    double the = theta/double(n);		// Pulse angle for each step
    double td = 1/v;				// Time for 1 cycle at frequency v (Hz)
    if(rad) td *= 2.0*PI;			// Adjust if frequency given in radians
    gen_op Udante = prop(H, td);		// Delay propagator (no relaxation)
    Udante *= Ixypuls_U(sys, Iso, phi, the);	// Multiply into pulse propagator (ideal)
    return pow(Udante,n); 			// Apply propagator n times for DANTE
    }
*/

// ____________________________________________________________________________
//                         DANTE INTERACTIVE FUNCTIONS
// ____________________________________________________________________________

// sosi - this function makes us depende upon the Level2 module

double ask_DANTE(const spin_system& sys, const string& Iso, gen_op& H, double cutoff)

	// Input		sys	: A spin system
	//			Iso     : String designating an isotope
	//			H	: Isotropic Hamiltonian (Hz)
	//			cutoff	: An intensity cutoff
	// Output		v       : A transition of H (Hz)
	// Note				: This routine looks over all single

  {
  gen_op detect = Fm(sys, Iso);			// Set detection operator to F-
  gen_op sigma = Fx(sys);			// Set system in transverse plane
  super_op L = complexi*Hsuper(H);		// L = -i*[Ho, ] (rad/sec)
  acquire1D ACQ(detect, L);			// Prepare for "Dirac" acquisitions
  TTable1D trans = ACQ.table(sigma);		// Array of transitions
  int ntr = trans.size();			// Get the number of transitions
  int itr = 0;
  double v;
  cout << "\n\n\tThere are " << ntr << " " << Iso << " Transitions\n";
  for(int i=0; i<ntr; i++)
    cout << "\n\t\t" << i << ". " << trans.Fr(i)/PIx2;
  cout << "\n\n\tChoose [0-" << ntr-1 << "] To Set DANTE Repetition,"
       << " Any Other Integer To Specify A Different Value: ";
  cin >> itr;
  if(itr >= ntr || itr < 0)
    {
    cout << "\n\t\t" << "Please Input A Value: ";
    cin >> v;
    }
  else
    v = trans.Fr(itr)/PIx2; 
  cout << "\n";
  return v;
  }


void set_DANTE(double gamB1, double& tmix, double& tpul, double tdel, int& numb, int& type)

	// Input		gamB1: RF-Field strength (Hz)
	// 			tmix : Total DANTE mixing time
	// 			tpul : Individual pulse length
	// 			tdel : Individual delay length
	// 			numb : Times DANTE needs repeating for tmix
	// 			type : Flag for DANTE type
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
    cout << "\n\n\t\tDANTE MIXING SEQUENCE SETUP\n";

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
    string typ;
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
      repeat = 1;				// Get the type of DANTE
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

#endif 						// PulDANTE.cc
