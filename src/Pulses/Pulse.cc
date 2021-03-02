/* Pulse.cc ****************************************************-*-c++-*-*
**                                                                      **
**                                   G A M M A                          **
**                                                                      **
**      Pulse Functions                            Implementation	**
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
** This GAMMA module facilitates the use of generic pulses in NMR       **
** simulations.                                                         **
**                                                                      **
*************************************************************************/

#ifndef _Pulse_cc_			// Is this file already included?
#define _Pulse_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Pulses/Pulse.h>		// Include the interface
#include <Pulses/PulAuxil.h>		// Include auxiliary functions
#include <Pulses/PulWaveform.h>		// Include pulse waveforms
#include <Pulses/PulComposite.h>	// Include composite pulses
#include <Pulses/PulTrain.h>		// Include pulse trains
#include <HSLib/SpinOpCmp.h>		// Include composite spin operators
#include <HSLib/SpinOpRot.h>		// Include spin rotation operators
#include <HSLib/SpinSystem.h>
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets
#include <Basics/StringCut.h>


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                       CLASS PULSE ERROR HANDLING
// ____________________________________________________________________________


void Pulse::Pulerror(int eidx, int noret) const

        // Input                Pulse	: Pulse parameters
        //                      eidx    : Error flag
        //                      noret   : Return flag
        // Output               none    : Error Message Output

  {
  std::cout << "\n\tPulse Parameters: ";
  switch(eidx)
    {
    case 0:                                                             // (0)
      std::cout << "Program Aborting....";
      break;
    case 1:                                                             // (1)
      std::cout << "Error During Construction";
      break;
    case 2:								// (2)
      std::cout << "Cannot Find An Adequate Set of Pulse Parameters";
      break;
    case 3:								// (3)
      std::cout << "Cannot Have No Delay With No/Ideal Pulse";
      break;
    default:
      std::cout << "Unknown Error (Number " << eidx << ")";
    }
  if(!noret) std::cout << ".\n";
  else       std::cout << ".";
  }

    
void volatile Pulse::Pulsefatality(int eidx) const

        // Input                Pulse	: Pulse parameters
        //                      eidx    : Error flag
        // Output               none    : Stops execution & error Message

  {                                                                       
  Pulerror(eidx,1);
  if(eidx) Pulerror(0);
  GAMMAfatal();					// Clean exit from program
  }


 
void Pulse::Pulerror(int eidx, const std::string& pname, int noret) const
 
        // Input                Pulse	: Pulse parameters
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
 
  {
  std::cout << "\nPulse Parameters: ";
  switch(eidx)
    {
    case 40:                                                    // (40)
      std::cout << "Problems with File " << pname;
      break;
    case 100:                                                   // (100)
      std::cout << "Can't Read Parameter " << pname;
      break;
    case 101:                                                   // (101)
      std::cout << "Can't Find Pulse Parameters For " << pname;
      break;
    case 130:                                                   // (130)
      std::cout << "Parameter " << pname << " Is The Culprit!\n";
      break;
    default:
      std::cout << "Unknown error";
      break;
    }
  if(!noret) std::cout << ".\n";
  }  

// ____________________________________________________________________________
// ii                CLASS PULSE PARAMETER SET FUNCTIONS
// ____________________________________________________________________________


void Pulse::SetPhase(const ParameterSet& pset, int idx)
 
        // Input		Pul	: Pulse parameters
        //                      pset    : Parameter set
        //                      idx     : Pulse index
        // Output               none    : Pulse phase read in
        //                                from parameter in pset
        //                                with index idx

  {
  double phiin;
  std::string pstate;				// Dummy string variable
  std::string pname = std::string("Pphi");		// Pulse phase angle
  std::string SI = std::string("(") + Gdec(idx)            // Name adjustment if indexed
            + std::string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  ParameterSet::const_iterator item;    // A pix into parameter std::list
  item = pset.seek(pname);			// Parameter for Pulsephi
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,phiin,pstate); 		// 	Get Pulse phase
    phi= phiin;					// 	Set Pulse phase
    }
  else phi = 0;					// If not found set to 0
  }


void Pulse::SetChannel(const ParameterSet& pset, int idx)
 
        // Input		Pul	: Pulse parameters
        //                      pset    : Parameter set
        //                      idx     : Pulse index
        // Output               none    : Pulse channel read in
        //                                from parameter in pset
        //                                with index idx

  {
  std::string Diso;					// Variable for reading
  std::string pstate;				// Dummy string variable
  std::string pname = std::string("Piso");		// Channel selectivity
  std::string SI = std::string("(") + Gdec(idx)            // Name adjustment if indexed
            + std::string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed
  ParameterSet::const_iterator item;         // A pix into parameter std::list
  item = pset.seek(pname);			// Parameter for Pulseiso
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,Diso,pstate);		// 	Read selectivity
    Iso = Diso;					// 	Set Pulse channel
    }
  else Iso = "";				// If not found, don't set
  }


int Pulse::SetAngle(const ParameterSet& pset, int idx)
 
        // Input		Pul	: Pulse parameters
        //                      pset    : Parameter set
        //                      idx     : Pulse index
        // Output               TF	: Pulse angle read in
        //                                from parameter in pset
        //                                with index idx

  {
  double pang;					// Variable for reading
  std::string pstate;				// Dummy string variable
  std::string pname = std::string("Pang");	// Pulse angle
  std::string SI = std::string("(") + Gdec(idx)	// Name adjustment if indexed
            + std::string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  ParameterSet::const_iterator item;	// A pix into parameter std::list
  item = pset.seek(pname);			// Parameter for Pulseang
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,pang,pstate); 		// 	Get Pulse angle
    ang = pang;					// 	Set Pulse angle
    return 1;					//	Return TRUE
    }
  else ang = 0;					// If not found set to 0
  return 0;					// Return FALSE
  }


int Pulse::SetPulLen(const ParameterSet& pset, int idx)
 
        // Input		Pul	: Pulse parameters
        //                      pset    : Parameter set
        //                      idx     : Pulse index
        // Output               TF	: Pulse length read in
        //                                from parameter in pset
        //                                with index idx

  {
  double plen;					// Variable for reading
  std::string pstate;				// Dummy string variable
  std::string pname = std::string("Plen");		// Pulse length
  std::string SI = std::string("(") + Gdec(idx)            // Name adjustment if indexed
            + std::string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  ParameterSet::const_iterator item;         // A pix into parameter std::list
  item = pset.seek(pname);			// Parameter for Pulsetp
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,plen,pstate); 		// 	Get Pulse length
    tp = plen;					// 	Set Pulse length
    return 1;					//	Return TRUE
    }
  else tp = 0;					// If not found set to 0
  return 0;					// Return FALSE
  }


int Pulse::SetGamB1(const ParameterSet& pset, int idx)
 
        // Input		Pul	: Pulse parameters
        //                      pset    : Parameter set
        //                      idx     : Pulse index
        // Output               TF	: Pulse strength read in
        //                                from parameter in pset
        //                                with index idx

  {
  double gB1;					// Variable for reading
  std::string pstate;				// Dummy string variable
  std::string pname = std::string("PgamB1");		// Pulse strength
  std::string SI = std::string("(") + Gdec(idx)            // Name adjustment if indexed
            + std::string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  ParameterSet::const_iterator item;         // A pix into parameter std::list
  item = pset.seek(pname);			// Parameter for PulsegamB1
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,gB1,pstate); 		// 	Get Pulse pul strength
    gamB1 = gB1;				// 	Set Pulse pul strength
    return 1;					//	Return TRUE
    }
  else gB1 = 0;					// If not found set to 0
  return 0;					// Return FALSE
  }

 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A                  CLASS PULSE CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 
 
Pulse::Pulse()
 
        // Input        none    :
        // Output	Pulse	: Pulse parameters (this)
 
  {
  Iso       = "";			// No pulse channel
  tp	    = 0;			// No pulse time
  gamB1     = 0;			// No rf-field strength	
  phi       = 0;			// No rf-field phase
  Wrf       = 0;			// No rf-field offset
  ang       = 0; 			// No pulse angle
  }
 
         
Pulse::Pulse(const std::string& ch, double gB1, double len, double ph, double off)

        // Input        ch      : RF-channel
        //              gB1     : RF-field strength (Hz)
	//		len	: RF-pulse length (sec)
        //              ph      : RF-phase (degrees)
        //              off     : RF-offset (Hz)
        // Output       Pulse   : Pulse parameters(this) 

  {
  Iso   = ch;				// Set pulse channel
  gamB1 = gB1;				// Set rf-field strength	
  tp	= len;				// Set pulse length
  phi   = ph;				// Set rf-field phase
  Wrf   = off;				// Set rf-field offset
  ang   = gB1*len*360.0;		// Set pulse angle
  }

                                                                                
Pulse::Pulse(const Pulse& Pulse1)
 
        // Input       Pulse1  : Pulse parameters
        // Output       Pulse   : Pulse parameters(this) from Pulse1

  {
  Iso   = Pulse1.Iso;			// Copy pulse channel
  tp	= Pulse1.tp;			// Copy pulse time
  gamB1 = Pulse1.gamB1;			// Copy rf-field strength	
  phi   = Pulse1.phi;			// Copy rf-field phase
  Wrf   = Pulse1.Wrf;			// Copy rf-field offset
  ang   = Pulse1.ang; 			// Copy pulse angle
  }
 
                                                                              
// ------------------------------ Destruction ---------------------------------
 
                                                                                
Pulse::~Pulse() {}
 
        // Input       Pulse1  : Pulse parameters (this)
        // Output       none    : Pulse is destructed
 
                                                    
// ------------------------------- Assignment ---------------------------------
 
                                         
Pulse& Pulse::operator = (const Pulse& Pulse1)

        // Input       Pulse1  : Pulse parameters
        // Output       Pulse   : Pulse parameters(this) from Pulse1
        ///F_std::list =           - Assignment

{
  Iso   = Pulse1.Iso;			// Copy pulse channel
  tp	= Pulse1.tp;			// Copy pulse time
  gamB1 = Pulse1.gamB1;			// Copy rf-field strength	
  phi   = Pulse1.phi;			// Copy rf-field phase
  Wrf   = Pulse1.Wrf;			// Copy rf-field offset
  ang	= Pulse1.ang;			// Copy pulse angle

  return (*this);
}
 
 
// ____________________________________________________________________________
// B                     CLASS PULSE ACCESS FUNCTIONS
// ____________________________________________________________________________
 

std::string     Pulse::channel()  const { return Iso; }
double     Pulse::strength() const { return gamB1; }
double     Pulse::angle()    const { return ang; }
double     Pulse::phase()    const { return phi; }
double     Pulse::offset()   const { return Wrf; }
double     Pulse::length()   const { return tp; }

        // Input       Pulse1  : Pulse parameters
        // Output       channel : Pulse isotope channel
	//		strength: Pulse strength (Hz)
	//		length	: Pulse length (sec)
	//		angle   : Pulse angle (sec)
	//		phase   : Pulse phase (sec)
	//		offset  : Pulse offset (Hz)

void Pulse::strength(double gB1)

        // Input	Pul	: Pulse parameters
	//		gB1	: Pulse strength (Hz)
	// Output	void	: Pulse field strength is set
// sosix - still have to think about ideal pulses!

 { gamB1 = gB1; }

 
// ____________________________________________________________________________
// C                  CLASS PULSE HILBERT SPACE PROPAGATORS
// ____________________________________________________________________________
 

/*
HSprop Pulse::GetU(const spin_system& sys, gen_op& H)

	// Input	     sys   : Spin system
	// 		     H     : Static Hamlitonian without the field
        //		     D	   : Pulse specifications
	// Output	     U     : Propagator for a Pulse step
	// Note			   : The propagator must be built in REVERSE
	//			     order for the pulse preceeds the delay

  {
  HSprop Udante(H, te);			// Delay propagator
  if(tp) Udante *= Pulse::GetU(sys,H); 	// Pulse propagator
  return Udante;
  }
*/

 

// ____________________________________________________________________________
// Y                      CLASS PULSE INPUT FUNCTIONS
// ____________________________________________________________________________


void Pulse::read(const std::string &filename, int idx)

        // Input		Pul	: Pulse parameters
        //                      idx     : Pulse index
        // Output               none    : Pulse parameters are read in
        //                                from parameters in file filename
        //                                with index idx
              
  {
  ParameterSet pset; 	               // Declare a parameter set
  if(!pset.read(filename, 1))          // Read in pset from file
    {   
    Pulerror(40, filename);		// Filename problems
    Pulsefatality(21);			// Fatal error
    }   
  read(pset, idx);			// User overloaded function
  return;
  }
 

void Pulse::read(const ParameterSet& pset, int idx)
 
        // Input		Pul	: Pulse parameters
        //                      pset    : Parameter set
        //                      idx     : Pulse index
        // Output               none    : Pulse parameters are read in
        //                                from parameters in pset
        //                                with index idx
 
  {
//	        First We'll Set Up The Pulse Pulse Parameters
// For The Pulse We Need Two of Three: {Pulse Angle, RF Strength, Pulse Length} 

// Note: Reading order is angle --> length --> strength for Pulse.  If only
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
    if(!T || !G) Pulsefatality(2); 		// Quit if either isn't set
    }
  SetChannel(pset,idx);				// Set Pulse channel
  SetPhase(pset, idx);				// Set Pulse phase
  }


std::string Pulse::ask_read(int argc, char* argv[], int argn)

        // Input		Pul	: Pulse parameters
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the Pulse parameters
	//				  are read
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : Pulse parameters are modifed (filled)

  {
  std::string filename;                            // Name of spin system file
  query_parameter(argc, argv, argn,           // Get filename from command
    "\n\tPulse parameters filename? ", filename); // Or ask for it
  read(filename);                             // Read system from filename
  return filename;
  }



// ____________________________________________________________________________
// Z                      CLASS PULSE I/O FUNCTIONS
// ____________________________________________________________________________
 
 
std::ostream& Pulse::printBase(std::ostream &ostr) const
 
        // Input		Pul	: Pulse parameters
        //                      ostr    : Output stream
        // Output               none    : Pulse basic parameters are sent
        //                                to the output stream
 
  {
  ostr << "\n\tIsotope channel:           " << Iso; 
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
 
 
std::ostream& Pulse::print(std::ostream &ostr, int full) const
 
        // Input		Pul	: Pulse parameters
        //                      ostr    : Output stream
        //                      full    : Flag for output amount
        // Output               none    : Pulse parameters are sent
        //                                to the output stream
 
  {
  ostr << "\n\n\t\t\tPulse Parameters\n";
  printBase(ostr);
  ostr <<"\n\n";
  return ostr;
  }
 

std::ostream &operator << (std::ostream &ostr, const Pulse &Pul)
 
        // Input		Pul	: Pulse parameters
        //                      ostr    : Output stream
        // Output               none    : Pulse parameters are sent
        //                                to the output stream
 
  {
  Pul.print(ostr);
  return ostr;
  }

#endif 						// Pulse.cc
