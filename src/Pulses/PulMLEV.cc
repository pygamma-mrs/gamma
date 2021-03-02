/* PulMLEV.cc ***************************************************-*-c++-*-
**									**
**      	                G A M M A				**
**									**
**	MLEV Pulse Functions      	   Implementation		**
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
** This file contains the a variety of functions supporting the use     **
** of MLEV pulses and pulse trains in the GAMMA simulation platform.    **
** The basic MLEV waveforms are composite 180 pulses and these are      **
** cycled to build up longer MLEV based pulse trains.  MLEV sequences   **
** are used in both mixing and broad-band decoupling.  For details see  **
**                                                                      **
*************************************************************************/

#ifndef _MLEV_cc_			// Is this file already included?
#define _MLEV_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif


# include <Pulses/PulMLEV.h>		// Include header info
# include <Pulses/PulAuxil.h>		// Include pulse auxilary functions
# include <Pulses/PulWaveform.h>	// Know pulse waveforms
# include <Pulses/PulComposite.h>	// Know composite pulses
# include <Pulses/PulCycle.h>		// Know pulse cycles
//# include <PulSupCycle.h>		// Know pulse supercycles
# include <Pulses/PulTrain.h>		// Know pulse trains
#include <Matrix/row_vector.h>		// Know about row vectors
#include <Basics/StringCut.h>		// Include Gform and Gdec functions
#include <string>                       // Include libstdc++ strings
#include <list>				// Include libstdc++ STL lists

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ lists
using std::cout;			// Using libstdc++ standard output

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS MLEV ERROR HANDLING
// ____________________________________________________________________________


void MLEV::MLEVerror(int eidx, int noret) const

        // Input                MLEV   : MLEV parameters
        //                      eidx    : Error flag
        //                      noret   : Return flag
        // Output               none    : Error Message Output

  {
  cout << "\n\tMLEV Parameters: ";
  switch(eidx)
    {
    case 0:                                                             // (0)
      cout << "Program Aborting....";
      break;
    case 1:                                                             // (1)
      cout << "Error During Construction";
      break;
    case 2:                                                             // (2)
      cout << "Cannot Find An Adequate Set of Pulse Parameters";
      break;
    case 3:                                                             // (3)
      cout << "Cannot Have No Delay With No/Ideal Pulse";
      break;
    default:
      cout << "Unknown Error (Number " << eidx << ")";
    }
  if(!noret) cout << ".\n";
  else       cout << ".";
  }


void volatile MLEV::MLEVfatality(int eidx) const
 
        // Input                MLEV   : MLEV parameters
        //                      eidx    : Error flag
        // Output               none    : Stops execution & error Message

  {                                                                       
  MLEVerror(eidx,1);
  if(eidx) MLEVerror(0);
  GAMMAfatal();					// Clean exit from program
  }


void MLEV::MLEVerror(int eidx, const string& pname, int noret) const
 
        // Input                MLEV   : MLEV parameters
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
 
  {
  cout << "\n\tMLEV Parameters: ";
  switch(eidx)
    {
    case 21:                                                    // (21)
      cout << "Cannot Construct MLEV From Parameters in File " << pname;
      break;
    case 40:                                                    // (40)
      cout << "Problems with File " << pname;
      break;
    case 100:                                                   // (100)
      cout << "Can't Read Parameter " << pname;
      break;
    case 101:                                                   // (101)
      cout << "Can't Find MLEV Parameters For " << pname;
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


void volatile MLEV::MLEVfatality(int eidx, const string& pname, int noret) const
 
        // Input                MLEV	: MLEV parameters
        //                      eidx    : Error flag
	//			pname   : Part of error message
        //                      noret   : Flag for return (0=return)
        // Output               none    : Stops execution & error Message

  {                                                                       
  MLEVerror(eidx, pname, 1);
  if(eidx) MLEVerror(0);
  GAMMAfatal();					// Clean exit from program
  }


 
 
// ____________________________________________________________________________
// ii                CLASS MLEV PARAMETER SET FUNCTIONS
// ____________________________________________________________________________


void MLEV::SetPhase(const ParameterSet& pset, int idx)
 
        // Intput               MP      : MLEV parameters
        //                      pset    : Parameter set
        //                      idx     : MLEV index
        // Output               none    : MLEV pulse phase read in
        //                                from parameter in pset
        //                                with index idx
 
  {
  double phiin;
  string pstate;                                // Dummy string variable
  string pname = string("MLEVphi");            // MLEV phase angle
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;                     // Adjust name if indexed
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);                         // Parameter for MLEVphi
  if(item != pset.end())                                      // If parameter found
    {
    (*item).parse(pname,phiin,pstate);     //      Get MLEV phase
    phi= phiin;                                 //      Set MLEV phase
    }
  else phi = 0;                                 // If not found set to 0
  }  


void MLEV::SetChannel(const ParameterSet& pset, int idx)
 
        // Intput               MP      : MLEV parameters
        //                      pset    : Parameter set
        //                      idx     : MLEV index
        // Output               none    : MLEV pulse channel read in
        //                                from parameter in pset
        //                                with index idx
 
  {
  string Diso;                                  // Variable for reading
  string pstate;                                // Dummy string variable
  string pname = string("MLEViso");             // Channel selectivity
  string SI = string("(") + Gdec(idx)           // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;                     // Adjust name if indexed
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);                      // Parameter for MLEViso
  if(item != pset.end())                        // If parameter found
    {
    (*item).parse(pname,Diso,pstate);           //      Read selectivity
    Iso = Diso;                                 //      Set MLEV channel
    }
  else Iso = string("");                   // If not found, don't set
  }  

 
int MLEV::SetGamB1(const ParameterSet& pset, int idx)
 
        // Intput               DT      : MLEV parameters
        //                      pset    : Parameter set
        //                      idx     : MLEV index
        // Output               TF      : MLEV pulse strength read in
        //                                from parameter in pset
        //                                with index idx
 
  {
  double gB1;                                   // Variable for reading
  string pstate;                                // Dummy string variable
  string pname = string("MLEVgamB1");          // MLEV pulse strength
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;                     // Adjust name if indexed
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);                         // Parameter for MLEVgamB1
  if(item != pset.end())                                      // If parameter found
    {
    (*item).parse(pname,gB1,pstate);       //      Get MLEV pul strength
    gamB1 = gB1;                                //      Set MLEV pul strength
    return 1;                                   //      Return TRUE
    }
  else gB1 = 0;                                 // If not found set to 0
  return 0;                                     // Return FALSE
  }
 
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A                  CLASS MLEV CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 

MLEV::MLEV() : Pulse() { }
 
        // Input        none    :
        // Output       MLEV   : MLEV parameters (this)
 
 
MLEV::MLEV(double gB1, const string& ch, double ph, double off)
      :Pulse(ch, gB1, 0.25/gB1, ph, off) { }

        // Input        ch      : RF-channel 
        //              gB1     : RF-field strength (Hz)
        //              ph      : RF-phase (degrees)
        //              off     : RF-offset (Hz)
        //              num     : Number of waveforms
        // Note                 : The base pulse is set to a 90 pulse
        //                        of strength gB1
 

MLEV::MLEV(const MLEV& MLEV1) : Pulse(MLEV1) { }
 
        // Intput       MLEV1  : MLEV parameters
        // Output       MLEV   : MLEV parameters(this) from MLEV1

                                                                     
// ------------------------------ Destruction ---------------------------------
 

MLEV::~MLEV() {}
 
        // Intput       MLEV1  : MLEV parameters (this)
        // Output       none    : MLEV is destructed
 
// ------------------------------- Assignment ---------------------------------
 
 
MLEV& MLEV::operator = (const MLEV& MLEV1) { Pulse::operator=(MLEV1); return (*this); }

        // Intput       MLEV1  : MLEV parameters                         
        // Output       MLEV   : MLEV parameters(this) from MLEV1
        ///F_list =             - Assignment

                                             
// ____________________________________________________________________________
// B                     CLASS MLEV ACCESS FUNCTIONS
// ____________________________________________________________________________

/*                                                                              
string     MLEV::channel()  const { return Iso; }      INHERITED
double     MLEV::strength() const { return gamB1; }    INHERITED
double     MLEV::phase()    const { return phi; }      INHERITED
double     MLEV::offset()   const { return Wrf; }      INHERITED            */  

        // Intput       MLEV1  : MLEV parameters                              
        // Output       channel : MLEV isotope channel
        //              strength: MLEV pulse strength (Hz)
        //              phase   : MLEV pulse phase (sec)
        //              offset  : MLEV pulse offset (Hz)

                                                          
// ____________________________________________________________________________
// C                  CLASS MLEV HILBERT SPACE PROPAGATORS
// ____________________________________________________________________________
 

 
// ____________________________________________________________________________
// D                CLASS MLEV LIOUVILLE SPACE PROPAGATORS
// ____________________________________________________________________________

                                                                                
// ____________________________________________________________________________
// E                        MLEV WAVEFORM FUNCTIONS
// ____________________________________________________________________________

/******************************************************************************

     There are 3 steps in a MLEV composite 180 sequence, as shown below.

                        U = U (90)*U (180)*U (90)
                             x      y       x

     Each U represents a pulse of the specified angle about the axis
     indicates by the subscript.

******************************************************************************/

PulWaveform MLEV::WF( ) const { return WF_C180(); }
PulWaveform MLEV::WF_C180( ) const

        // Input        MP      : MLEV parameters
        // Output       PWF     : Pulse waveform for MLEV composite 180
        // Note                 : Constant rf-amplitude each step
        // Note                 : Ideal pulses are allowed

  {                                                                             
  row_vector WFsteps(3);                        // Waveform vector{ gamB1,phi }
  WFsteps.put(complex(gamB1,0.0),    0);        // Set waveform values
  WFsteps.put(complex(gamB1,90.0),   1);
  WFsteps.put(complex(gamB1,0.0),    2);
  row_vector WFtimes(3);                        // Vector for step times
  double tdegree = 1/(gamB1*360);               // Length per pulse degree
  WFtimes.put(90.0 *tdegree,  0);
  WFtimes.put(180.0*tdegree,  1);
  WFtimes.put(90.0 *tdegree,  2);
  return PulWaveform(WFsteps, WFtimes, "Composite 180");
  }
 

// ____________________________________________________________________________
// F                    MLEV COMPOSITE PULSE FUNCTIONS
// ____________________________________________________________________________


/******************************************************************************

     There are 3 steps in a MLEV composite 180 sequence, as shown below.

                        U = U (90)*U (180)*U (90)
                             x      y       x

     Each U represents a pulse of the specified angle about the axis
     indicates by the subscript.

******************************************************************************/

PulComposite MLEV::PCmp(const spin_system& sys) const {return PCmp_C180(sys);}
PulComposite MLEV::PCmp_C180(const spin_system& sys) const
 
        // Input        MP      : MLEV parameters
        //              sys     : A spin system
        // Output       PCmp	: Composite pulse for MLEV composite 180
        //                        applicable to sys

  { return PulComposite(WF(), sys, Iso); }


// ____________________________________________________________________________
// G                      MLEV PULSE CYCLE FUNCTIONS
// ____________________________________________________________________________

/****************************************************************************** 
 
     There are 4 cycles associated with MLEV-4, as shown below
                             _ _
                     U = R R R R        R = 90 180  90
                                              x   y   x

     Each R represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A MLEV-4 pulse train will employ R = MLEV
     which is just a composite 180 pulse.  The MLEV-4 cycle is simply the
     phi phi phi+180 phi+180 phase sequence.
 
******************************************************************************/
 

PulCycle MLEV::CycMLEV4(const spin_system& sys) const

        // Input        MP      : MLEV parameters       
        //              sys     : A spin system
        // Output       PCyc    : MLEV-4 pulse cycle for sys

  { return PulCycle(PCmp(sys), CYC_MLEV4(), "MLEV-4"); }

 
/*I****************************************************************************
 
     There are 8 cycles associated with MLEV-8, as shown below
                       _ _ _     _
               U = R R R R R R R R           R = 90 180 90
                                                   x   y  x
 
     Each R represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A MLEV-8 pulse train will employ R = MLEV
     which is just a composite 180 pulse.  The MLEV-8 cycle is simply the
     phi phi phi+180 phi+180 phi+180 phi phi phi+180 phase sequence.
 
******************************************************************************/

 
PulCycle MLEV::CycMLEV8(const spin_system& sys) const

        // Input        MP      : MLEV parameters
        //              sys     : A spin system
        // Output       PCyc    : MLEV-8 pulse cycle for sys

  { return PulCycle(PCmp(sys), CYC_MLEV8(), "MLEV-8"); }

 

/* ****************************************************************************
 
     There are 16 cycles associated with MLEV-16, as shown below
                    _ _ _     _ _ _       _ _
            U = R R R R R R R R R R R R R R R R          R = 90 180 90
                                                               x   y  x

     Each R represents a pulse waveform and the bar indicates a         
     180 degree phase shift.  A MLEV-16 pulse train will employ R = MLEV
     which is just a composite 180 pulse.  The MLEV-16 cycle is simply the
     phi phi phi+180 phi+180 phi+180 phi phi phi+180 ..... phase sequence.
 
******************************************************************************/

PulCycle MLEV::CycMLEV16(const spin_system& sys) const

        // Input        MP      : MLEV parameters
        //              sys     : A spin system
        // Output       PCyc    : MLEV-16 pulse cycle for sys

  { return PulCycle(PCmp(sys), CYC_MLEV16(), "MLEV-16"); }
 
 
// ____________________________________________________________________________
// H                    MLEV PULSE SUPERCYCLE FUNCTIONS
// ____________________________________________________________________________
 
// ____________________________________________________________________________
// I                        MLEV PULSE TRAIN FUNCTIONS
// ____________________________________________________________________________

// ____________________________________________________________________________ 
// K                      CLASS MLEV INPUT FUNCTIONS
// ____________________________________________________________________________
 

void MLEV::read(const string &filename, int idx)

        // Intput               ML      : MLEV parameters
        //                      idx     : MLEV index
        // Output               none    : MLEV parameters are read in
        //                                from parameters in file filename
        //                                with index idx

  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename, 1))		// Read in pset from file
    {
    MLEVerror(40, filename, 1);		// Filename problems
    MLEVfatality(21, filename);		// Fatal error
    }
  read(pset, idx);                      // User overloaded function
  return;
  }


void MLEV::read(const ParameterSet& pset, int idx)

        // Intput               ML      : MLEV parameters
        //                      pset    : Parameter set
        //                      idx     : MLEV index
        // Output               none    : MLEV parameters are read in
        //                                from parameters in pset
        //                                with index idx
 
  {
//              First We'll Set Up The MLEV Pulse Parameters
// For The Pulse We Need One of Two: {RF Strength, Composite Pulse Length}
 
// Note: Reading order is strength --> length for MLEV.
 
  int G = SetGamB1(pset, idx);                  // If strength set, set pulse
  if(!G) MLEVfatality(2);                      // Quit if no strength set
  SetChannel(pset,idx);                         // Set MLEV pulse channel
  SetPhase(pset, idx);                          // Set MLEV pulse phase
  }


    
void MLEV::ask_read(int argc, char* argv[], int argn)
 
        // Intput               ML      : MLEV parameters
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the MLEV parameters
        //                                are read
        //                                If the argument argn is not in argv,
        //                                the user is asked for a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : MLEV parameters are modifed (filled)

  {                                                                             
  string filename;                            // Name of spin system file
  query_parameter(argc, argv, argn,           // Get filename from command
    "\n\tMLEV parameters filename? ", filename); // Or ask for it
  read(filename);                             // Read system from filename
  }


// ____________________________________________________________________________
// L                      CLASS MLEV I/O FUNCTIONS
// ____________________________________________________________________________
 
 
std::ostream& MLEV::printBase(std::ostream &ostr) const
 
        // Intput               MP      : MLEV parameters
        //                      ostr    : Output stream
        // Output               none    : MLEV basic parameters are sent
        //                                to the output stream
 
  {
  ostr << "\n\tIsotope channel:";
  ostr << string(13-Iso.length(), ' ');
  ostr << Iso;
  ostr << "\n\tRF strength:             ";
  printHz(ostr, gamB1);
  ostr << "\n\tRF phase:                ";
  ostr << Gform("%8.3f", phi) << " degrees";
  ostr << "\n\tRF offset:               ";
  printHz(ostr, Wrf);
  return ostr;
  }
         
 
std::ostream& MLEV::print(std::ostream &ostr) const
 
        // Intput               MP      : MLEV parameters
        //                      ostr    : Output stream
        // Output               none    : MLEV parameters are sent
        //                                to the output stream
 
  {
  ostr << "\n\n\t\t\tMLEV Parameters\n";
  printBase(ostr);
  ostr <<"\n\n";
  return ostr;
  }



std::ostream &operator << (std::ostream &ostr, const MLEV &MP)
 
        // Intput               MP      : MLEV parameters
        //                      ostr    : Output stream
        // Output               none    : MLEV parameters are sent
        //                                to the output stream
 
  {
  MP.print(ostr);
  return ostr;
  }


// ____________________________________________________________________________
// AA               ADDITIONAL MLEV PHASE CYCLE FUNCTIONS
// ____________________________________________________________________________


/* **************************************************************************** 
 
     There are 4 cycles associated with MLEV-4, as shown below
                             _ _
                     U = R R R R        R = 90 180  90
                                              x   y   x

     Each R represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A MLEV-4 pulse train will employ R = MLEV
     which is just a composite 180 pulse.  The MLEV-4 cycle is simply the
     phi phi phi+180 phi+180 phase sequence.

******************************************************************************/

row_vector CYC_MLEV4(double phi)

        // Input        void    : None
        //              phi     : Phase angle (degrees)
        // Output       PCyc    : MLEV-4 pulse cycle

  {
  double phibar = phi+180.0;                    // Phase of "barred" steps 
  row_vector PTCsteps(4);                       // Cycle vector{ phi } 
  PTCsteps.put(phi,    0);                      // Set cycle phase values
  PTCsteps.put(phi,    1);
  PTCsteps.put(phibar, 2); 
  PTCsteps.put(phibar, 3); 
  return PTCsteps;
  }

 
/* ****************************************************************************
 
     There are 8 cycles associated with MLEV-8, as shown below
                       _ _ _     _
               U = R R R R R R R R           R = 90 180 90
                                                   x   y  x
 
     Each R represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A MLEV-8 pulse train will employ R = MLEV
     which is just a composite 180 pulse.  The MLEV-8 cycle is simply the
     phi phi phi+180 phi+180 phi+180 phi phi phi+180 phase sequence.
 
******************************************************************************/

row_vector CYC_MLEV8(double phi)

        // Input        void    : None
        //              phi     : Phase angle (degrees)
        // Output       PCyc    : MLEV-8 pulse cycle
 

  {
  double phibar = phi+180.0;                    // Phase of "barred" steps 
  row_vector PTCsteps(8);                       // Cycle vector{ phi } 
  PTCsteps.put(phi,    0);                      // Set cycle phase values
  PTCsteps.put(phi,    1);
  PTCsteps.put(phibar, 2); 
  PTCsteps.put(phibar, 3); 
  PTCsteps.put(phibar, 4); 
  PTCsteps.put(phi,    5);
  PTCsteps.put(phi,    6);
  PTCsteps.put(phibar, 7); 
  return PTCsteps;
  }


/* ****************************************************************************
 
     There are 16 cycles associated with MLEV-16, as shown below
                    _ _ _     _ _ _       _ _
            U = R R R R R R R R R R R R R R R R          R = 90 180 90
                                                               x   y  x

     Each R represents a pulse waveform and the bar indicates a         
     180 degree phase shift.  A MLEV-16 pulse train will employ R = MLEV
     which is just a composite 180 pulse.  The MLEV-16 cycle is simply the
     phi phi phi+180 phi+180 phi+180 phi phi phi+180 ..... phase sequence.
 
******************************************************************************/

row_vector CYC_MLEV16(double phi)

        // Input        void    : None
        //              phi     : Phase angle (degrees)
        // Output       PCyc    : MLEV-16 pulse cycle
 
  {
  double phibar = phi+180.0;                    // Phase of "barred" steps 
  row_vector PTCsteps(16);			// Cycle vector{ phi } 
  PTCsteps.put(phi,    0);                      // Set cycle phase values
  PTCsteps.put(phi,    1);
  PTCsteps.put(phibar, 2); 
  PTCsteps.put(phibar, 3); 
  PTCsteps.put(phibar, 4); 
  PTCsteps.put(phi,    5);
  PTCsteps.put(phi,    6);
  PTCsteps.put(phibar, 7); 
  PTCsteps.put(phibar, 8); 
  PTCsteps.put(phibar, 9); 
  PTCsteps.put(phi,   10);
  PTCsteps.put(phi,   11);
  PTCsteps.put(phi,   12);
  PTCsteps.put(phibar,13); 
  PTCsteps.put(phibar,14); 
  PTCsteps.put(phi,   15);
  return PTCsteps;
  }

#endif						// PulMLEV.cc
