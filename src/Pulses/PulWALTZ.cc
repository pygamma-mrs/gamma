/* PulWALTZ.cc **************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      WALTZ Pulse Functions                      Implementation	**
**                                                                      **
**      Copyright (c) 1998                                              **
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
**  This file contains the a variety of functions supporting the use	**
**  of WALTZ pulse trains in the GAMMA simulation platform.  The basic	**
**  WALTZ waveforms are composite pulses and these are cycled to build	**
**  up longer WALTZ based pulse trains.  WALTZ is primarily used as a	**
**  broad-band decoupling sequence.  For details see			**
**                                                                      **
** See J. Magn. Reson., 53, 313-340 (1983), A.J. Shaka, James Keeler,	**
** and Ray Freeman "Evaluation of a New Broadband Decoupling Sequence:	**
** WALTZ-16"								**
**                                                                      **
** WALTZ waveforms (base composite pulses): WALTZ-R, WALTZ-Q, WALTZ-K	**
** WALTZ pulse trains (cycled waveforms):   WALTZ-4, WALTZ-8, WALTZ-16	**
**                                                                      **
*************************************************************************/

#ifndef _WALTZ_cc_				// Is this file already included?
#define _WALTZ_cc_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)	// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Pulses/PulWALTZ.h>			// Include header info
#include <Pulses/PulAuxil.h>			// Include pulse auxiliary
#include <Matrix/row_vector.h>			// Know about row vectors
#include <Pulses/PulWaveform.h>			// Know pulse waveforms
#include <Pulses/PulComposite.h>		// Know composite pulses
#include <Pulses/PulCycle.h>			// Know pulse cycles
#include <Basics/ParamSet.h>	
#include <Basics/StringCut.h>			// Include Gdec and Gform
#include <string>				// Include libstdc++ strings
#include <list>					// Include libstdc++ STL lists
#include <iostream>				// Include input output streams

using std::list;				// Using libstdc++ lists
using std::string;				// Using libstdc++ strings
using std::cout;
using std::ostream;				// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS WALTZ ERROR HANDLING
// ____________________________________________________________________________


void WALTZ::WALTZerror(int eidx, int noret) const

        // Input                WALTZ	: WALTZ parameters
        //                      eidx    : Error flag
        //                      noret   : Return flag
        // Output               none    : Error Message Output

  {
  cout << "\n\tWALTZ Parameters: ";
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
 
 
void volatile WALTZ::WALTZfatality(int eidx) const
 
        // Input                WALTZ	: WALTZ parameters
        //                      eidx    : Error flag
        // Output               none    : Stops execution & error Message
 
  {
  WALTZerror(eidx,1);
  if(eidx) WALTZerror(0);
  GAMMAfatal();					// Clean exit from program
  }
 

void WALTZ::WALTZerror(int eidx, const string& pname, int noret) const

        // Input                WALTZ	: WALTZ parameters
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message

  {                                                     
  cout << "\nWALTZ Parameters: ";
  switch(eidx)
    {
    case 40:                                                    // (40)
      cout << "Problems with File " << pname;
      break;
    case 100:                                                   // (100)
      cout << "Can't Read Parameter " << pname;
      break;
    case 101:                                                   // (101)
      cout << "Can't Find WALTZ Parameters For " << pname;
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
// ii                CLASS WALTZ PARAMETER SET FUNCTIONS
// ____________________________________________________________________________


void WALTZ::SetPhase(const ParameterSet& pset, int idx)
 
        // Intput		DT	: WALTZ parameters
        //                      pset    : Parameter set
        //                      idx     : WALTZ index
        // Output               none    : WALTZ pulse phase read in
        //                                from parameter in pset
        //                                with index idx

  {
  double phiin;
  string pstate;				// Dummy string variable
  string pname = string("WALTZphi");		// WALTZ phase angle
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Parameter for WALTZphi
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,phiin,pstate); 		// 	Get WALTZ phase
    phi= phiin;					// 	Set WALTZ phase
    }
  else phi = 0;					// If not found set to 0
  }


void WALTZ::SetChannel(const ParameterSet& pset, int idx)
 
        // Intput		DT	: WALTZ parameters
        //                      pset    : Parameter set
        //                      idx     : WALTZ index
        // Output               none    : WALTZ pulse channel read in
        //                                from parameter in pset
        //                                with index idx

  {
  string Diso;				// Variable for reading
  string pstate;				// Dummy string variable
  string pname = string("WALTZiso");	// Channel selectivity
  string SI = string("(") + Gdec(idx)	// Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Parameter for WALTZiso
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,Diso,pstate);		// 	Read selectivity
    Iso = Diso;					// 	Set WALTZ channel
    }
  else Iso = string("");			// If not found, don't set
  }


int WALTZ::SetGamB1(const ParameterSet& pset, int idx)
 
        // Intput		DT	: WALTZ parameters
        //                      pset    : Parameter set
        //                      idx     : WALTZ index
        // Output               TF	: WALTZ pulse strength read in
        //                                from parameter in pset
        //                                with index idx

  {
  double gB1;					// Variable for reading
  string pstate;				// Dummy string variable
  string pname = string("WALTZgamB1");		// WALTZ pulse strength
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Parameter for WALTZgamB1
  if(item != pset.end())			// If parameter found
    {
    (*item).parse(pname,gB1,pstate); 		// 	Get WALTZ pul strength
    gamB1 = gB1;				// 	Set WALTZ pul strength
    return 1;					//	Return TRUE
    }
  else gB1 = 0;					// If not found set to 0
  return 0;					// Return FALSE
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                  CLASS WALTZ CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________


WALTZ::WALTZ() : Pulse() { }

        // Input        none    :
        // Output	WALTZ	: WALTZ parameters (this)


WALTZ::WALTZ(double gB1, const string& ch, double ph, double off)
      :Pulse(ch, gB1, 0.25/gB1, ph, off) { }
 
        // Input	ch      : RF-channel
	//	        gB1     : RF-field strength (Hz)
        //              ph      : RF-phase (degrees) 
	//		off     : RF-offset (Hz)
	//		num	: Number of waveforms
	// Note			: The base pulse is set to a 90 pulse 
	//			  of strength gB1

    
WALTZ::WALTZ(const WALTZ& WALTZ1) : Pulse(WALTZ1) { }

        // Intput       WALTZ1  : WALTZ parameters
        // Output       WALTZ   : WALTZ parameters(this) from WALTZ1
 
    
// ------------------------------ Destruction ---------------------------------

     
WALTZ::~WALTZ() {}

        // Intput       WALTZ1	: WALTZ parameters (this)
        // Output       none    : WALTZ is destructed

                                                      
// ------------------------------- Assignment ---------------------------------


WALTZ& WALTZ::operator = (const WALTZ& WALTZ1) 
{ 
  Pulse::operator=(WALTZ1); 
  return (*this);
}
 
        // Intput       WALTZ1  : WALTZ parameters
        // Output       WALTZ   : WALTZ parameters(this) from WALTZ1
        ///F_list =		- Assignment
 
 
// ____________________________________________________________________________
// B                     CLASS WALTZ ACCESS FUNCTIONS
// ____________________________________________________________________________
 
/*
string     WALTZ::channel()  const { return Iso; }	INHERITED
double     WALTZ::strength() const { return gamB1; }	INHERITED
double     WALTZ::phase()    const { return phi; }	INHERITED
double     WALTZ::offset()   const { return Wrf; }	INHERITED            */
 
        // Intput       WALTZ1	: WALTZ parameters
        // Output 	channel : WALTZ isotope channel
        //              strength: WALTZ pulse strength (Hz)
        //              phase   : WALTZ pulse phase (sec)
        //              offset  : WALTZ pulse offset (Hz)
 

// ____________________________________________________________________________
// C                  CLASS WALTZ HILBERT SPACE PROPAGATORS
// ____________________________________________________________________________

 

// ____________________________________________________________________________
// D                CLASS WALTZ LIOUVILLE SPACE PROPAGATORS
// ____________________________________________________________________________
 

// ____________________________________________________________________________
// E                        WALTZ WAVEFORM FUNCTIONS
// ____________________________________________________________________________

/* ****************************************************************************

     There are 3 steps in a WALTZ-R sequence, as shown below.
                                  _                _
                        U = U(90)*U(180)*U(270) = 123

     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift.

******************************************************************************/
 

PulWaveform WALTZ::WF(int even) const { return WF_WALTZR(even); }
PulWaveform WALTZ::WF_WALTZR(int even) const

        // Input        WP	: WALTZ parameters
	//		even	: Flag if even step sizes
        // Output       PWF	: Pulse waveform for WALTZ-1
        // Note                 : Constant rf-amplitude each step

  {
  double phibar = phi+180.0;			// Phase of "barred" steps
  int Nsteps = 3;				// Number of steps
  if(even) Nsteps = 6;				// Modify if even increments
  row_vector WFsteps(Nsteps);			// Waveform vector{ gamB1,phi }
  row_vector WFtimes(Nsteps);			// Vector for step times
  double tdegree = 0;                           // Increment time per pulse
  if(gamB1>0) tdegree = 1/(gamB1*360);          // degree
  string Name("WALTZ-R");
  if(!even)
    {
    WFsteps.put(complex(gamB1,phi),    0);	// Set waveform values
    WFsteps.put(complex(gamB1,phibar), 1);
    WFsteps.put(complex(gamB1,phi),    2);
    WFtimes.put(90.0 *tdegree,  0);		// Set step lengths
    WFtimes.put(180.0*tdegree,  1);
    WFtimes.put(270.0*tdegree, 2);
    }
  else
    {
    WFsteps.put(complex(gamB1,phi),    0);
    WFsteps.put(complex(gamB1,phibar), 1);
    WFsteps.put(complex(gamB1,phibar), 2);
    WFsteps.put(complex(gamB1,phi),    3);
    WFsteps.put(complex(gamB1,phi),    4);
    WFsteps.put(complex(gamB1,phi),    5);
    WFtimes.put(90.0*tdegree, 0);
    WFtimes.put(90.0*tdegree, 1);
    WFtimes.put(90.0*tdegree, 2);
    WFtimes.put(90.0*tdegree, 3);
    WFtimes.put(90.0*tdegree, 4);
    WFtimes.put(90.0*tdegree, 5);
    Name += string("*");
    }
  return PulWaveform(WFsteps, WFtimes, Name);
  }


/* ****************************************************************************

     There are 5 steps in a WALTZ-K sequence, as shown below.
                  _             _             _       _ _ _
              U = U(180)*U(360)*U(180)*U(270)*U(90) = 24231

     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift.

******************************************************************************/
 
PulWaveform WALTZ::WF_WALTZK(int even) const

        // Input        gamB1	: RF field strength (Hz)
	//		even	: Flag if even step sizes
	//		phi	: RF field phase
        // Output       PWF     : Composite pulse waveform
        // Note                 : Uses constant rf-amplitude.

  {
  double phibar = phi+180.0;			// Phase of "barred" steps
  int Nsteps=5;					// Set number of steps
  if(even) Nsteps = 12;				// More if even length steps 
  row_vector WFsteps(Nsteps);			// Waveform vector{ gamB1,phi }
  row_vector WFtimes(Nsteps);			// Vector for step times
  double tdegree = 0;                           // Increment time per pulse
  if(gamB1>0) tdegree = 1/(gamB1*360);          // degree
  if(!even)
    {
    WFsteps.put(complex(gamB1,phibar), 0);	// Set waveform values
    WFsteps.put(complex(gamB1,phi),    1);
    WFsteps.put(complex(gamB1,phibar), 2);
    WFsteps.put(complex(gamB1,phi),    3);
    WFsteps.put(complex(gamB1,phibar), 4);
    WFtimes.put(180.0*tdegree,  0);
    WFtimes.put(360.0*tdegree,  1);
    WFtimes.put(180.0*tdegree,  2);
    WFtimes.put(270.0*tdegree,  3);
    WFtimes.put( 90.0*tdegree,  4);
    }
  else
    {
    WFsteps.put(complex(gamB1,phibar), 0);
    WFsteps.put(complex(gamB1,phibar), 1);
    WFsteps.put(complex(gamB1,phi),    2);
    WFsteps.put(complex(gamB1,phi),    3);
    WFsteps.put(complex(gamB1,phi),    4);
    WFsteps.put(complex(gamB1,phi),    5);
    WFsteps.put(complex(gamB1,phibar), 6);
    WFsteps.put(complex(gamB1,phibar), 7);
    WFsteps.put(complex(gamB1,phi),    8);
    WFsteps.put(complex(gamB1,phi),    9);
    WFsteps.put(complex(gamB1,phi),   10);
    WFsteps.put(complex(gamB1,phibar),11);
    WFtimes.put(90.0*tdegree,  0);
    WFtimes.put(90.0*tdegree,  1);
    WFtimes.put(90.0*tdegree,  2);
    WFtimes.put(90.0*tdegree,  3);
    WFtimes.put(90.0*tdegree,  4);
    WFtimes.put(90.0*tdegree,  5);
    WFtimes.put(90.0*tdegree,  6);
    WFtimes.put(90.0*tdegree,  7);
    WFtimes.put(90.0*tdegree,  8);
    WFtimes.put(90.0*tdegree,  9);
    WFtimes.put(90.0*tdegree, 10);
    WFtimes.put(90.0*tdegree, 11);
    }
  return PulWaveform(WFsteps, WFtimes, "WALTZ-K");
  }
 

/* ****************************************************************************

     There are 9 steps in a WALTZ-Q sequence, as shown below.
      _             _             _            _             _        _ _ _ _ _
  U = U(270)*U(360)*U(180)*U(270)*U(90)*U(180)*U(360)*U(180)*U(270) = 342312423

     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift.

******************************************************************************/


PulWaveform WALTZ::WF_WALTZQ(int even) const

        // Input        tp      : Single Pi Pulse length (sec)
	//		even	: Flag if even step sizes
        // Output       Cvect   : Vector of composite pulse waveform
        // Note                 : Uses constant rf-amplitude and increment
        //                        time, thus is constructed as 4 steps

  {
  double phibar = phi+180.0;			// Phase of "barred" steps
  int Nsteps = 9;				// Number of steps
  if(even) Nsteps = 24;				// More if even length steps 
  row_vector WFsteps(Nsteps);			// Waveform vector{ gamB1,phi }
  row_vector WFtimes(Nsteps);			// Vector for step times
  double tdegree = 0;                           // Increment time per pulse
  if(gamB1>0) tdegree = 1/(gamB1*360);          // degree
  if(!even)
    {
    WFsteps.put(complex(gamB1,phibar), 0);        // Set waveform values
    WFsteps.put(complex(gamB1,phi),    1);
    WFsteps.put(complex(gamB1,phibar), 2);
    WFsteps.put(complex(gamB1,phi),    3);
    WFsteps.put(complex(gamB1,phibar), 4);
    WFsteps.put(complex(gamB1,phi),    5);
    WFsteps.put(complex(gamB1,phibar), 6);
    WFsteps.put(complex(gamB1,phi),    7);
    WFsteps.put(complex(gamB1,phibar), 8);
    WFtimes.put(270.0*tdegree,  0);
    WFtimes.put(360.0*tdegree,  1);
    WFtimes.put(180.0*tdegree,  2);
    WFtimes.put(270.0*tdegree,  3);
    WFtimes.put( 90.0*tdegree,  4);
    WFtimes.put(180.0*tdegree,  5);
    WFtimes.put(360.0*tdegree,  6);
    WFtimes.put(190.0*tdegree,  7);
    WFtimes.put(270.0*tdegree,  8);
    }
  else
    {
    WFsteps.put(complex(gamB1,phibar), 0);	// Set waveform values
    WFsteps.put(complex(gamB1,phibar), 1);	// Set waveform values
    WFsteps.put(complex(gamB1,phibar), 2);	// Set waveform values
    WFsteps.put(complex(gamB1,phi),    3);
    WFsteps.put(complex(gamB1,phi),    4);
    WFsteps.put(complex(gamB1,phi),    5);
    WFsteps.put(complex(gamB1,phi),    6);
    WFsteps.put(complex(gamB1,phibar), 7);
    WFsteps.put(complex(gamB1,phibar), 8);
    WFsteps.put(complex(gamB1,phi),    9);
    WFsteps.put(complex(gamB1,phi),   10);
    WFsteps.put(complex(gamB1,phi),   11);
    WFsteps.put(complex(gamB1,phibar),12);
    WFsteps.put(complex(gamB1,phi),   13);
    WFsteps.put(complex(gamB1,phi),   14);
    WFsteps.put(complex(gamB1,phibar),15);
    WFsteps.put(complex(gamB1,phibar),16);
    WFsteps.put(complex(gamB1,phibar),17);
    WFsteps.put(complex(gamB1,phibar),18);
    WFsteps.put(complex(gamB1,phi),   19);
    WFsteps.put(complex(gamB1,phi),   20);
    WFsteps.put(complex(gamB1,phibar),21);
    WFsteps.put(complex(gamB1,phibar),22);
    WFsteps.put(complex(gamB1,phibar),23);
    for(int i=0; i<Nsteps; i++)
      WFtimes.put(90.0*tdegree,  i);
    }
  return PulWaveform(WFsteps, WFtimes, "WALTZ-Q");
  }


// ____________________________________________________________________________
// F                     WALTZ COMPOSITE PULSE FUNCTIONS
// ____________________________________________________________________________

/* ****************************************************************************

     There are 3 steps in a WALTZ-R sequence, as shown below.
                                  _                _
                        U = U(90)*U(180)*U(270) = 123

     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift.

******************************************************************************/

PulComposite WALTZ::PCmp(const spin_system& S, int even) const
                                               { return PCmpWALTZR(S, even); }
PulComposite WALTZ::PCmpWALTZR(const spin_system& sys, int even) const

        // Input        WP	: WALTZ parameters
	//		sys	: A spin system
	//		even	: Flag if even step sizes
        // Output       CP	: Composite pulse for WALTZ-R
	//			  applicable to sys

  { return PulComposite(WF(even), sys, Iso); }


/* ****************************************************************************

     There are 5 steps in a WALTZ-K sequence, as shown below.
                  _             _             _       _ _ _
              U = U(180)*U(360)*U(180)*U(270)*U(90) = 24231

     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift.

******************************************************************************/


PulComposite WALTZ::PCmpWALTZK(const spin_system& sys, int even) const
 
        // Input        WP	: WALTZ parameters
	//		sys	: A spin system
	//		even	: Flag if even step sizes
        // None         CP	: A composite pulse is returned set
        //                        to WALTZ-K for the system sys
        //                        of strength gamB1 and phase phi
 
  { return PulComposite(WF_WALTZK(even), sys, Iso); }
 

/* ****************************************************************************

     There are 9 steps in a WALTZ-Q sequence, as shown below.
      _             _             _            _             _        _ _ _ _ _
  U = U(270)*U(360)*U(180)*U(270)*U(90)*U(180)*U(360)*U(180)*U(270) = 342312423

     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift.

******************************************************************************/
 
PulComposite WALTZ::PCmpWALTZQ(const spin_system& sys, int even) const
 
        // Input        WP	: WALTZ parameters
	//		sys	: A spin system
	//		even	: Flag if even step sizes
        // None         CP	: A composite pulse is returned set
        //                        to WALTZ-Q for the system sys
        //                        of strength gamB1 and phase phi
 
  { return PulComposite(WF_WALTZQ(even), sys, Iso); }


// ____________________________________________________________________________
// G                      WALTZ PULSE CYCLE FUNCTIONS
// ____________________________________________________________________________


/* ****************************************************************************

     There are 4 cycles associated with WALTZ-4, as shown below
                     _ _    _   _  _ _ _ _
             U = R R R R = 123 123 123 123       R = 90 180  270
                                                       x   -x   x
 
     Each R represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A WALTZ4 pulse train will employ R = WALTZ-R.
     The WALTZ-4 cycle is simply the phi phi phi+180 phi+180 sequence 

******************************************************************************/


PulCycle WALTZ::CycWALTZ4(const spin_system& sys, int even) const

        // Input        void    : None
	//		phi	: Phase angle (degrees)
	//		even	: Flag if even step sizes
        // Output       PCyc	: WALTZ-4 pulse cycle

  { return PulCycle(PCmp(sys, even), CYC_WALTZ4(), "WALTZ-4"); }


/* ****************************************************************************

     There are 4 cycles associated with WALTZ-8, as shown below
         _ _     _ _ _    _ _     _  _   _ _ _
   U = K K K K = 24231 * 24231 * 24231 * 24231     K = 180  360 180  270  90
                                                          -x   x   -x   x   -x
 
     Each K represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A WALTZ8 pulse train will employ K = WALTZ-K.
     The WALTZ-8 cycle is simply the phi phi+180 phi+180 phi sequence.

******************************************************************************/


PulCycle WALTZ::CycWALTZ8(const spin_system& sys, int even) const

        // Input        void    : None
	//		phi	: Phase angle (degrees)
	//		even	: Flag if even step sizes
        // Output       PCyc	: WALTZ-8 pulse cycle

  { return PulCycle(PCmpWALTZK(sys, even), CYC_WALTZ8(), "WALTZ-8"); }

/******************************************************************************

     There are 4 cycles associated with WALTZ-16, as shown below
               _ _     _ _ _ _ _    _ _ _ _     _ _ _ _    _ _ _ _ _
         U = Q Q Q Q = 342312423 * 342312423 * 342312423 * 342312423

                  Q = 270  360 180  270 90  180 360  180 270
                         -x   x   -x   x  -x   x   -x   x   -x
 
     Each Q represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A WALTZ16 pulse train will employ Q = WALTZ-Q

******************************************************************************/

 
PulCycle WALTZ::CycWALTZ16(const spin_system& sys, int even) const

        // Input        void    : None
	//		phi	: Phase angle (degrees)
	//		even	: Flag if even step sizes
        // Output       PCyc	: WALTZ-16 pulse cycle
 
  { return PulCycle(PCmpWALTZK(sys, even), CYC_WALTZ8(), "WALTZ-16"); }


// ____________________________________________________________________________
// H                    WALTZ PULSE SUPERCYCLE FUNCTIONS
// ____________________________________________________________________________


// ____________________________________________________________________________
// I                        WALTZ PULSE TRAIN FUNCTIONS
// ____________________________________________________________________________


/* According to the original WALTZ-1 article
                           _ _   _ __ _ __  _    _      
                   U = R R R R = PQPQ PQPQ PQPQ PQPQ

where each R is a single WALTZ-1 composite pulse, and the - indicates a
180 phase shift.  Indeed, this is just the cycle used in WALTZ-4             */


//PulTrain WALTZ::PT(const spin_system& sys) const { return PT_WALTZ1(sys); }
//PulTrain WALTZ::PT_WALTZ1(const spin_system& sys) const

        // Input        WP	: WALTZ parameters
	//		sys	: A spin system
        // Output       PT	: Pulse train for WALTZ-1
	//			  applicable to sys

//  { return PulTrain(CP(sys), CYC(), "WALTZ-1"); }


// ____________________________________________________________________________
// J                    CLASS WALTZ INTERACITVE FUNCTIONS
// ____________________________________________________________________________


/*
void WALTZ1_ask(int argc, char* argv[], int& qn, double& gamB1, string& IsoG)

	// Input	argc	: No. arguments
        //              argv    : Argument strings
        //              qn      : Query number
        //              gamB1   : RF-Field strength (kHz)
	//		IsoG	: Decoupling channel
	// None		void    : The rf-field strength for a WALTZ sequence
	//			  is interactively requested unless supplied
	//			  in the array argv.

  {
  query_parameter(argc, argv, qn++, 
    "\n\tWALTZ Decoupling Field Strength (KHz)? ",gamB1);
  gamB1 *= 1.e3;
  query_parameter(argc, argv, qn++, 
    "\n\tWALTZ Decoupling Channel (1H, 13C, ...)? ", IsoG);
  return;
  }


void read_WALTZ1(string& filein, spin_sys& sys, 
                                         double& gamB1, string& IsoG, int idx)

	// Input	filein  : Input parameter file name
        //              sys     : Active spin system
        //              gamB1   : WALTZ RF-Field strength (Hz)
	//		IsoG	: WALTZ Isotope channel
	//		idx     : Parameter name qualifier
	// None		void    : The 2 WALTZ decoupling parameters are set
	//			  from the parameter values specified in 
	//			  the input parameter file.

  {
  ParameterSet pset;				// A parameter set
  pset.read(filein);                            // Read pset in
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  string pname, ssfile, pstate;                 // Items in each pset entry
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");

//               Determine the WALTZ Decoupling Isotope Channel

  if(!sys.heteronuclear()) IsoG = sys.symbol(0);
  else
    {
    pname = string("IsoG");			// Chirp Isotope Parameter
    if(idx >= 0) pname += SI;			// Adjust name if indexed
    item = pset.seek(pname);			// Pix in parameter list
    if(item != pset.end())
      (*item).parse(pname,IsoG,pstate);// Set Detected Isotope
    else
      {
      cout << "\n\tCan't Read WALTZ Isotope Channel Parameter "
          << pname << "!\n";
      exit(-1);
      }
    }

//             Determine the Field Strength in the WALTZ Sequence
 
  pname = string("gB1WALTZ1");			// RF Field strength of WALTZ
  if(idx >= 0) pname += SI;                     // Adjust name if indexed
  item = pset.seek(pname);                        // Pix in parameter list
  if(item != pset.end())
    (*item).parse(pname,gamB1,pstate);	// Set the WALTZ field strength
  else
    {
    cout << "\n\tCan't Read WALTZ Field Strength Parameter "
          << pname << "!\n";
    exit(-1);
    }
  }
*/

// ____________________________________________________________________________
// K                      CLASS WALTZ INPUT FUNCTIONS
// ____________________________________________________________________________


void WALTZ::read(const string &filename, int idx)

        // Intput               WP      : WALTZ parameters
        //                      idx     : WALTZ index
        // Output               none    : WALTZ parameters are read in
        //                                from parameters in file filename
        //                                with index idx

  {
  ParameterSet pset;                // Declare a parameter set
  if(!pset.read(filename, 1))          // Read in pset from file
    {
    WALTZerror(40, filename);           // Filename problems
    WALTZfatality(21);                  // Fatal error
    }
  read(pset, idx);                      // User overloaded function
  return;
  } 


void WALTZ::read(const ParameterSet& pset, int idx)

        // Intput               WP      : WALTZ parameters
        //                      pset    : Parameter set
        //                      idx     : WALTZ index
        // Output               none    : WALTZ parameters are read in
        //                                from parameters in pset
        //                                with index idx

  {
//              First We'll Set Up The WALTZ Pulse Parameters
// For The Pulse We Need One of Two: {RF Strength, Composite Pulse Length}

// Note: Reading order is strength --> length for WALTZ.

  int G = SetGamB1(pset, idx);                  // If strength set, set pulse
  if(!G) WALTZfatality(2); 			// Quit if no strength set
  SetChannel(pset,idx);                         // Set WALTZ pulse channel
  SetPhase(pset, idx);                          // Set WALTZ pulse phase
  }
    
 
void WALTZ::ask_read(int argc, char* argv[], int argn)
 
        // Intput               WP      : WALTZ parameters
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the WALTZ parameters
        //                                are read
        //                                If the argument argn is not in argv,
        //                                the user is asked for a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : WALTZ parameters are modifed (filled)
 
  {
  string filename;                            // Name of spin system file
  query_parameter(argc, argv, argn,           // Get filename from command
    "\n\tWALTZ parameters filename? ", filename); // Or ask for it
  read(filename);                             // Read system from filename
  }
 

 
// ____________________________________________________________________________
// L                      CLASS WALTZ I/O FUNCTIONS
// ____________________________________________________________________________


/*
ostream& WALTZ::printBase(ostream &ostr) const

        // Intput               WP      : WALTZ parameters
        //                      ostr    : Output stream
        // Output               none    : WALTZ basic parameters are sent
        //                                to the output stream

  {                                                            
  ostr << "\n\tIsotope channel:";
  ostr << replicate(" ", 13-Iso.length());
  ostr << Iso;
  ostr << "\n\tRF strength:             ";
  printHz(ostr, gamB1);
  ostr << "\n\tRF phase:                ";
  ostr << Gform("%8.3f", phi) << " degrees";
  ostr << "\n\tRF offset:               ";
  printHz(ostr, Wrf);
  return ostr;
  }
*/


ostream& WALTZ::print(ostream &ostr) const

        // Intput               WP      : WALTZ parameters
        //                      ostr    : Output stream
        // Output               none    : WALTZ parameters are sent
        //                                to the output stream

  {                                                            
  ostr << "\n\n\t\t\tWALTZ Parameters\n";
  printBase(ostr);
  ostr <<"\n\n";
  return ostr;
  }
 
 
 
ostream &operator << (ostream &ostr, const WALTZ &WP)

        // Intput               WP      : WALTZ parameters
        //                      ostr    : Output stream
        // Output               none    : WALTZ parameters are sent
        //                                to the output stream

  {                                                            
  WP.print(ostr);
  return ostr;
  }
 
// ____________________________________________________________________________
// AA               ADDITIONAL WALTZ PHASE CYCLE FUNCTIONS
// ____________________________________________________________________________

/******************************************************************************

     A WALTZ-4 pulse cycle is given by
                     _ _    _   _  _ _ _ _
             U = R R R R = 123 123 123 123       R = 90 180  270
                                                       x   -x   x

     Each R represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A WALTZ4 pulse train will employ R = WALTZ-R
 
******************************************************************************/

row_vector CYC_WALTZ4(double phi)

        // Input        void    : None
        //              phi     : Phase angle (degrees)
        // Output       PCyc    : WALTZ-4 pulse cycle

  {
  double phibar = phi+180.0;			// Phase of "barred" steps
  row_vector PTCsteps(4);			// Cycle vector{ phi }
  PTCsteps.put(phi,    0);			// Set cycle phase values
  PTCsteps.put(phi,    1);
  PTCsteps.put(phibar, 2);
  PTCsteps.put(phibar, 3);
  return PTCsteps;
  }


/******************************************************************************

     A WALTZ-8 pulse cycle is given by
         _ _     _ _ _    _ _     _  _   _ _ _
   U = K K K K = 24231 * 24231 * 24231 * 24231     K = 180  360 180  270  90
                                                          -x   x   -x   x   -x

     Each K represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A WALTZ8 pulse train will employ K = WALTZ-K

******************************************************************************/


row_vector CYC_WALTZ8(double phi)

        // Input             sys   : Spin system
        //                   H     : Static Hamlitonian without the field
        //                   Iso   : Isotope channel the field is on
        //                   gamB1 : The rf-field strength (Hz)
        // Output            U     : Propagator for an WALTZ-8 sequence
        //                           about the axis specified by phi
        // Note                    : This propagator assumes one is in
        //                           the rotating frame of the rf-field

  {
  double phibar = phi+180.0;			// Phase of "barred" steps
  row_vector PTCsteps(4);			// Cycle vector{ phi }
  PTCsteps.put(phi,    0);			// Set cycle phase values
  PTCsteps.put(phibar, 1);
  PTCsteps.put(phibar, 2);
  PTCsteps.put(phi,    3);
  return PTCsteps;
  }






 
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________

// ____________________________________________________________________________
// E                      WALTZ PULSE TRAIN FUNCTIONS
// ____________________________________________________________________________


/*   There are 4 cycles associated with WALTZ-4, as shown below
                     _ _    _   _  _ _ _ _
             U = R R R R = 123 123 123 123       R = 90 180  270
                                                       x   -x   x
 
     Each R represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A WALTZ4 pulse train will employ R = WALTZ-R   */

 
//PulTrain PT_WALTZ4(const spin_system& sys, const string& IsoC,
//                                                     double gamB1, double phi)
 
        // Input        sys     : Active spin system
        //              IsoC    : Pulse channel
        //              gamB1   : RF field strength (Hz)
        //              phi     : RF field phase (degrees)
        // None         PT 	: A pulse train is returned set
        //                        to WALTZ-4 for the system sys
        //                        of strength gamB1 and phase phi

/*
  {
  PulComposite PCmp = PCmpWALTZR(sys, IsoC, gamB1, phi);
  PulCycle PCy = CYC_WALTZ4();
  return PulTrain(PCmp, PCy, "WALTZ-4");
  }
*/


/*   There are 4 cycles associated with WALTZ-8, as shown below
         _ _     _ _ _    _ _     _  _   _ _ _
   U = K K K K = 24231 * 24231 * 24231 * 24231     K = 180  360 180  270  90
                                                          -x   x   -x   x   -x
 
     Each K represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A WALTZ8 pulse train will employ K = WALTZ-K.
     The WALTZ-8 cycle is simply the phi phi+180 phi+180 phi sequence        */

 
//PulTrain PT_WALTZ8(const spin_system& sys, const string& IsoC,
//                                                     double gamB1, double phi)
 
        // Input        sys     : Active spin system
        //              IsoC    : Pulse channel
        //              gamB1   : RF field strength (Hz)
        //              phi     : RF field phase (degrees)
        // None         PT 	: A pulse train is returned set
        //                        to WALTZ-8 for the system sys
        //                        of strength gamB1 and phase phi

/*

  {
  PulComposite PCmp = PCmpWALTZK(sys, IsoC, gamB1, phi);
  PulCycle PCy = CYC_WALTZ8();
  return PulTrain(PCmp, PCy, "WALTZ-8");
  }
*/




/*   There are 4 cycles associated with WALTZ-16, as shown below
               _ _     _ _ _ _ _    _ _ _ _     _ _ _ _    _ _ _ _ _
         U = Q Q Q Q = 342312423 * 342312423 * 342312423 * 342312423

                  Q = 270  360 180  270 90  180 360  180 270
                         -x   x   -x   x  -x   x   -x   x   -x
 
     Each Q represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A WALTZ16 pulse train will employ Q = WALTZ-Q  */

 
//PulTrain PT_WALTZ16(const spin_system& sys, const string& IsoC,
//                                                     double gamB1, double phi)
 
        // Input        sys     : Active spin system
        //              IsoC    : Pulse channel
        //              gamB1   : RF field strength (Hz)
        //              phi     : RF field phase (degrees)
        // None         PT 	: A pulse train is returned set
        //                        to WALTZ-16 for the system sys
        //                        of strength gamB1 and phase phi

/*

  {
  PulComposite PCmp = PCmpWALTZQ(sys, IsoC, gamB1, phi);
  PulCycle PCy = CYC_WALTZ8();
  return PulTrain(PCmp, PCy, "WALTZ-16");
  }
*/


// ____________________________________________________________________________
// F                      WALTZ PROPAGATOR FUNCTIONS
// ____________________________________________________________________________



//HSprop WALTZ_R(spin_system& sys, gen_op& H, string& Iso,
//                                                    double gamB1, double phi)
 
        // Input             sys   : Spin system
        //                   H     : Static Hamlitonian without the field
        //                   Iso   : Isotope channel the field is on
        //                   gamB1 : The rf-field strength (Hz)
        // Output            U     : Propagator for a single WALTZ step
        //                           about the axis specified by phi
        // Note                    : This propagator assumes one is in
        //                           the rotating frame of the rf-field
//                                                         _
//    U = U (90, phi) * U (180, phi+180) * U (270, phi) = 123
//         p             p                  p
 
//  {
//  if(gamB1 <= 0.0) return prop(H.dim());	// If no applied field, exit
//  double Wrf = 0.0;				// Apply field on resonance
//  double t90 = 0.25/gamB1;			// Time for a 90 pulse
//  double t180 = 2.0*t90;			// Time for a 180 pulse
//  double t270 = 3.0*t90;			// Time for a 270 pulse
//  HSprop U1    = Pulse(sys, H, Iso, Wrf, t90, gamB1, phi);
//  HSprop U2bar = Pulse(sys, H, Iso, Wrf, t180, gamB1, phi+180);
//  HSprop U3    = Pulse(sys, H, Iso, Wrf, t270, gamB1, phi);
//  return U1*U2bar*U3;
// 
//  gen_op U;
//  double fact;
//  if(tp == 0.0)
//    U = Rxy(sys, iso, theta, phi);      // Rotation operator if tp=0
//  else
//    { 
//    fact = theta/(360.0*tp);
//    U = Spul_U_plane(sys, H, iso, freq, tp, fact, phi);
//    }  
//  return U;
//}    

#endif						// PulW_WALTZ.cc
