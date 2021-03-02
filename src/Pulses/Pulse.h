/* Pulse.h *****************************************************-*-c++-*-*
**									**
**	                          G A M M A 				**
**								 	**
**	Pulse Functions					Interface 	**
**								 	**
**	Copyright (c) 1998					 	**
**	Scott Smith						 	**
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**								 	**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** This GAMMA module facilitates the use of generic pulses in NMR	**
** simulations.								**
**								 	**
*************************************************************************/

#ifndef GPulse_h_			// Is file already included?
#define GPulse_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/SpinSystem.h>		// Include knowledge of spin systems
#include <HSLib/PulseI.h>		// Inlcude knowledge of ideal pulses
#include <HSLib/PulseS.h>		// Inlcude knowledge of pulses
#include <HSLib/HSham.h>		// Knowledge of isotropic Hamiltonians
#include <HSLib/GenOp.h>		// Inlucde knowledge of operators 
//#include <BWRRelax/relaxProp.h>		// Include knowledge of propagators
#include <string>			// Include knowledge of Strings

class Pulse
  {
  friend class WALTZ;			// WALTZ pulses are our friends
  friend class MLEV;			// MLEV pulses are our friends too
  std::string Iso;			// Applied pulse channel
  double gamB1;                         // Applied pulse strength (Hz)
  double tp;				// Applied pulse length   (sec)
  double ang;				// Applied pulse angle	  (deg)
  double phi;                           // Applied pulse phase    (deg)
  double Wrf;                           // Applied pulse offset   (Hz)

private:

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                      CLASS PULSE ERROR HANDLING
// ____________________________________________________________________________


void Pulerror(int eidx, int noret=0) const;
 
        // Input                Pulse   : Pulse parameters(this)
        //                      eidx    : Error flag
        //                      noret   : Return flag
        // Output               none    : Error Message Output
 
 
void volatile Pulsefatality(int error) const;
 
        // Input                Pulse   : Pulse parameters(this)
        //                      eidx    : Error flag
        // Output               none    : Stops execution & error Message


void Pulerror(int eidx, const std::string& pname, int noret=0) const;

        // Input                Pulse   : Pulse parameters(this)
        //                      eidx    : Error index
        //                      pname   : String included in message
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message


// ____________________________________________________________________________
// ii               CLASS Pulse PARAMETER SET FUNCTIONS
// ____________________________________________________________________________


void SetPhase(const ParameterSet& pset, int idx=-1);

        // Intput               Pul      : Pulse parameters
        //                      pset    : Parameter set
        //                      idx     : Pulse index
        // Output               none    : Pulse phase read in
        //                                from parameter in pset
        //                                with index idx


void SetChannel(const ParameterSet& pset, int idx=-1);

        // Intput               Pul      : Pulse parameters
        //                      pset    : Parameter set
        //                      idx     : Pulse index
        // Output               none    : Pulse channel read in
        //                                from parameter in pset
        //                                with index idx
 
 
int SetAngle(const ParameterSet& pset, int idx=-1);

        // Intput               Pul      : Pulse parameters
        //                      pset    : Parameter set
        //                      idx     : Pulse index
        // Output               TF      : Pulse angle read in
        //                                from parameter in pset
        //                                with index idx


int SetPulLen(const ParameterSet& pset, int idx=-1);

        // Intput               Pul      : Pulse parameters
        //                      pset    : Parameter set
        //                      idx     : Pulse index
        // Output               TF      : Pulse length read in
        //                                from parameter in pset
        //                                with index idx
 
 
int SetGamB1(const ParameterSet& pset, int idx=-1);
 
        // Intput               Pul      : Pulse parameters
        //                      pset    : Parameter set
        //                      idx     : Pulse index
        // Output               TF      : Pulse strength read in
        //                                from parameter in pset
        //                                with index idx


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
 public:
 
// ____________________________________________________________________________
// A                   CLASS Pulse CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 
 
MSVCDLC Pulse();
 
        // Input        none	:
        // Output	Pulse   : Pulse parameters(this)
        ///F_list       Pulse	- Constructor


MSVCDLC Pulse(const std::string& ch, double gB1, double len, double ph=0, double off=0);

        // Input        ch      : RF-channel
        //              gB1     : RF-field strength (Hz)
        //              len     : RF-pulse length (sec)
        //              ph      : RF-phase (degrees)
        //              off     : RF-offset (Hz)
        // Output       Pulse   : Pulse parameters(this)
 

MSVCDLC Pulse(const Pulse& PT1);
 
        // Intput	Pulse1	: Pulse parameters
        // Output	Pulse   : Pulse parameters(this) from Pulse1
 
 
// ------------------------------ Destruction ---------------------------------
 
 
MSVCDLC virtual ~Pulse();
 
        // Intput	Pulse1	: Pulse parameters (this)
        // Output       none    : Pulse is destructed
 
 
// ------------------------------- Assignment ---------------------------------

                                                                                
MSVCDLL Pulse& operator = (const Pulse& Pulse1);
 
        // Intput	Pulse1	: Pulse parameters
        // Output	Pulse   : Pulse parameters(this) from Pulse1
        ///F_list 	=	- Assignment
 

// ____________________________________________________________________________
// B                     CLASS Pulse ACCESS FUNCTIONS
// ____________________________________________________________________________

MSVCDLL std::string channel()  const;
MSVCDLL double      strength() const;
MSVCDLL double      angle()    const;
MSVCDLL double      phase()    const;
MSVCDLL double      offset()   const;
MSVCDLL double      length()   const;
 
        // Intput       Pulse1  : Pulse parameters
        // Output 	channel : Pulse isotope channel  
        //              strength: Pulse strength (Hz)
        //              length	: Pulse length (sec)
        //              angle   : Pulse angle (sec)
        //              phase   : Pulse phase (sec)
        //              offset  : Pulse offset (Hz)

        //              offset  : Pulse offset (Hz)
 
MSVCDLL void strength(double gB1);
 
        // Input        Pul     : Pulse parameters
        //              gB1     : Pulse strength (Hz)
        // Output       void    : Pulse field strength is set


// ____________________________________________________________________________
// C                  CLASS Pulse HILBERT SPACE PROPAGATORS
// ____________________________________________________________________________


/*
MSVCDLL HSprop GetU(const spin_system& sys, gen_op& H);

        // Input             sys   : Spin system
        //                   H     : Static Hamlitonian without the field
        //                   D     : Pulse specifications
        // Output            U     : Propagator for a Pulse step
        // Note                    : The propagator must be built in REVERSE
        //                           order for the pulse preceeds the delay
*/

// ____________________________________________________________________________
// D                CLASS Pulse LIOUVILLE SPACE PROPAGATORS
// ____________________________________________________________________________
 
 
// ____________________________________________________________________________
// E                       Pulse WAVEFORM FUNCTIONS
// ____________________________________________________________________________

 
// ____________________________________________________________________________
// F                    Pulse COMPOSITE PULSE FUNCTIONS
// ____________________________________________________________________________


// ____________________________________________________________________________
// G                       Pulse PULSE TRAIN FUNCTIONS
// ____________________________________________________________________________

// ____________________________________________________________________________
// Y                      CLASS Pulse INPUT FUNCTIONS
// ____________________________________________________________________________

 
MSVCDLL void read(const std::string &filename, int idx=-1);

        // Intput               Pul	: Pulse parameters
        //                      idx     : Pulse index
        // Output               none    : Pulse parameters are read in
        //                                from parameters in file filename
        //                                with index idx

    
MSVCDLL void read(const ParameterSet& pset, int idx=-1);

        // Intput               Pul	: Pulse parameters
        //                      pset    : Parameter set
        //                      idx     : Pulse index
        // Output               none    : Pulse parameters are read in
        //                                from parameters in pset
        //                                with index idx
 
//                  Read Pulse Delay & Pulse Parameters
// For The Pulse We Need Two of Three: {Pulse Angle, RF Strength, Pulse Length}
 
// Note: Reading order is angle --> length --> strength for Pulse.  If only
//       the angle is specified the length will be set to zero (ideal pulse)
 
 
MSVCDLL std::string ask_read(int argc, char* argv[], int argn);
 
        // Intput               Pul	: Pulse parameters
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the Pulse parameters
        //                                are read
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : Pulse parameters are modifed (filled)


// ____________________________________________________________________________
// Z                      CLASS Pulse I/O FUNCTIONS
// ____________________________________________________________________________


MSVCDLL std::ostream& printBase(std::ostream &ostr) const;
 
        // Intput               Pul	: Pulse parameters
        //                      ostr    : Output stream
        // Output               none    : Pulse basic parameters are sent
        //                                to the output stream


MSVCDLL virtual std::ostream& print(std::ostream &ostr, int full=0) const;

        // Intput               Pul	: Pulse parameters
        //                      ostr    : Output stream
        //                      full    : Flag for output amount
        // Output               none    : Pulse parameters are sent
        //                                to the output stream
 
    
MSVCDLL friend std::ostream &operator << (std::ostream &ostr, const Pulse &Pul);
 
        // Intput               Pul	: Pulse parameters
        //                      ostr    : Output stream
        // Output               none    : Pulse parameters are sent
        //                                to the output stream



  };


#endif 						// Pulse.h

