/* PulGARP.cc ***************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      GARP Pulse Waveform Functions               Implementation	**
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
**  This file contains the implementation of GARP pulse trains in the   **
**  GAMMA simulation platform.  GARP is used as a broad-band decoupling **
**  sequence.  For details see                                          **
**                                                                      **
** J. Magn. Reson., 64, 547-552 (1985), A.J. Shaka, P.B. Barker, &  	**
** Ray Freeman "Computer Optimized Decoupling Scheme for Wideband       **
** Applications and Low-Level Operation"                                **
**                                                                      **
** The 25 Steps (Pulses) In The GARP-1 Sequence Are As Follows:         **
**                                                                      **
**           ____       _____          ____       ____       _____      **
**     30.5  55.2 257.8 268.3  69.3    62.2  85.0 91.8 134.5 256.1      **
**           ____        ____         _____       ____        ____      **
**     66.4  45.9  25.5  72.7 119.5   138.2 258.4 64.9  70.9  77.2      **
**          _____        ____                                           **
**     98.2 133.6 255.9  65.6  53.4                                     **
**                                                                      **
** The above waveform is the repeated in the WALTZ-4 Cycle 	        **
**                                    _ _                               **
**                                R R R R                               **
**                                                                      **
** The total of rotation over the 25 steps will be 2857 degrees.  This	**
** implies that each GARP-1 waveform will take 15.87222 * t		**
**                                                         180		**
**                                                                      **
*************************************************************************/

#ifndef _pwgarp_cc_			// Is this file already included?
#define _pwgarp_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Pulses/PulGARP.h>		// Inlcude the header
#include <Pulses/PulAuxil.h>		// Know auxiliary functions
#include <Pulses/PulWaveform.h>		// Know pulse waveforms
#include <Pulses/PulComposite.h>	// Know composite pulses
#include <Pulses/PulCycle.h>		// Know pulse cycles
#include <Pulses/PulWALTZ.h>		// Know WALTZ pulse module (WALTZ4 cyc)
#include <Matrix/row_vector.h>		// Know about row vectors
#include <Basics/StringCut.h>		// Include Gform and Gdec functions
#include <string>                       // Include libstdc++ strings
#include <iostream>                     // Include input output streams (cout)
#include <list>				// Include libstdc++ STL lists

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ lists
using std::ostream;			// Using libstdc++ output streams
using std::cout;			// Using libstdc++ standard output


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS GARP ERROR HANDLING
// ____________________________________________________________________________


void GARP::GARPerror(int eidx, int noret) const

        // Input                GARP	: GARP parameters
        //                      eidx    : Error flag
        //                      noret   : Return flag
        // Output               none    : Error Message Output

  {
  cout << "\n\tGARP Parameters: ";
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
 
 
void volatile GARP::GARPfatality(int eidx) const
 
        // Input                GARP	: GARP parameters
        //                      eidx    : Error flag
        // Output               none    : Stops execution & error Message
 
  {
  GARPerror(eidx,1);
  if(eidx) GARPerror(0);
  GAMMAfatal();					// Clean exit from program
  }
 

void GARP::GARPerror(int eidx, const string& pname, int noret) const

        // Input                GARP	: GARP parameters
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message

  {                                                     
  cout << "\nGARP Parameters: ";
  switch(eidx)
    {
    case 40:                                                    // (40)
      cout << "Problems with File " << pname;
      break;
    case 100:                                                   // (100)
      cout << "Can't Read Parameter " << pname;
      break;
    case 101:                                                   // (101)
      cout << "Can't Find GARP Parameters For " << pname;
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
// ii                CLASS GARP PARAMETER SET FUNCTIONS
// ____________________________________________________________________________


void GARP::SetPhase(const ParameterSet& pset, int idx)
 
        // Intput		DT	: GARP parameters
        //                      pset    : Parameter set
        //                      idx     : GARP index
        // Output               none    : GARP pulse phase read in
        //                                from parameter in pset
        //                                with index idx

  {
  double phiin;
  string pstate;				// Dummy string variable
  string pname = string("GARPphi");		// GARP phase angle
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);				// Parameter for GARPphi
  if(item != pset.end())					// If parameter found
    {
    (*item).parse(pname,phiin,pstate); 	// 	Get GARP phase
    phi= phiin;					// 	Set GARP phase
    }
  else phi = 0;					// If not found set to 0
  }


void GARP::SetChannel(const ParameterSet& pset, int idx)
 
        // Intput		DT	: GARP parameters
        //                      pset    : Parameter set
        //                      idx     : GARP index
        // Output               none    : GARP pulse channel read in
        //                                from parameter in pset
        //                                with index idx

  {
  string Diso;					// Variable for reading
  string pstate;				// Dummy string variable
  string pname = string("GARPiso");		// Channel selectivity
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);				// Parameter for GARPiso
  if(item != pset.end())					// If parameter found
    {
    (*item).parse(pname,Diso,pstate);	// 	Read selectivity
    Iso = Diso;					// 	Set GARP channel
    }
  else Iso = "";				// If not found, don't set
  }


int GARP::SetGamB1(const ParameterSet& pset, int idx)
 
        // Intput		DT	: GARP parameters
        //                      pset    : Parameter set
        //                      idx     : GARP index
        // Output               TF	: GARP pulse strength read in
        //                                from parameter in pset
        //                                with index idx

  {
  double gB1;					// Variable for reading
  string pstate;				// Dummy string variable
  string pname = string("GARPgamB1");		// GARP pulse strength
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");
  if(idx >= 0) pname += SI;			// Adjust name if indexed
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);				// Parameter for GARPgamB1
  if(item != pset.end())					// If parameter found
    {
    (*item).parse(pname,gB1,pstate); 	// 	Get GARP pul strength
    gamB1 = gB1;				// 	Set GARP pul strength
    return 1;					//	Return TRUE
    }
  else gB1 = 0;					// If not found set to 0
  return 0;					// Return FALSE
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                  CLASS GARP CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________


GARP::GARP()

        // Input        none    :
        // Output	GARP	: GARP parameters (this)

  {                                                       
  Iso       = "";			// No rf-field channel
  gamB1     = 0;                        // No rf-field strength
  phi       = 0;                        // No rf-field phase
  Wrf       = 0;                        // No rf-field offset
  }


GARP::GARP(double gB1, const string& ch, double ph, double off)
 
        // Input	ch      : RF-channel
	//	        gB1     : RF-field strength (Hz)
        //              ph      : RF-phase (degrees) 
	//		off     : RF-offset (Hz)
	//		num	: Number of waveforms
        // Output       GARP    : GARP parameters(this)

  {
  Iso   = ch;				// Set rf-field channel
  gamB1 = gB1;				// Set rf-field strength
  phi   = ph;				// Set rf-field phase
  Wrf   = off;				// Set rf-field offset
  }

    
GARP::GARP(const GARP& GARP1)

        // Intput       GARP1  : GARP parameters
        // Output       GARP   : GARP parameters(this) from GARP1
 
  {
  Iso   = GARP1.Iso;			// Copy rf-field channel
  gamB1 = GARP1.gamB1;			// Copy rf-field strength
  phi   = GARP1.phi;			// Copy rf-field phase
  Wrf   = GARP1.Wrf;			// Copy rf-field offset
  }

    
// ------------------------------ Destruction ---------------------------------

     
GARP::~GARP() {}

        // Intput       GARP1	: GARP parameters (this)
        // Output       none    : GARP is destructed

                                                      
// ------------------------------- Assignment ---------------------------------


GARP& GARP::operator = (const GARP& GARP1)
 
        // Intput       GARP1  : GARP parameters
        // Output       GARP   : GARP parameters(this) from GARP1
        ///F_list =		- Assignment
 
  {
  Iso   = GARP1.Iso;			// Copy rf-field channel
  gamB1 = GARP1.gamB1;			// Copy rf-field strength
  phi   = GARP1.phi;			// Copy rf-field phase
  Wrf   = GARP1.Wrf;			// Copy rf-field offset

  return (*this);
  }

 
// ____________________________________________________________________________
// B                     CLASS GARP ACCESS FUNCTIONS
// ____________________________________________________________________________
 
 
string     GARP::channel()  const     { return Iso; }
double     GARP::strength() const     { return gamB1; }
void	   GARP::strength(double gB1) { gamB1 = gB1; }
double     GARP::phase()    const     { return phi; }
double     GARP::offset()   const     { return Wrf; }
 
        // Intput       GARP1	: GARP parameters
        // Output 	channel : GARP isotope channel
        //              strength: GARP pulse strength (Hz)
        //              phase   : GARP pulse phase (sec)
        //              offset  : GARP pulse offset (Hz)
 
 
// ____________________________________________________________________________
// C                  CLASS GARP HILBERT SPACE PROPAGATORS

// ____________________________________________________________________________
// C                  CLASS GARP HILBERT SPACE PROPAGATORS
// ____________________________________________________________________________

 
//PulTrain GARP1(string& filein, spin_sys& sys, int idx)

	// Input	filein  : Input parameter file name
        //              sys     : Active spin system
	//		idx     : Parameter name qualifier
	// None		PT      : A pulse train is returned set
	//			  to GARP-1 for the system sys
	//			  based on the parameters specified
	//			  in the input parameter file.
	//			: This pulse train has WALTZ-4 as the
	//			  supercycle, no cycle, and the sequence
	//			  set to a composite 25 step GARP-1

// According to the original GARP-1 article
//                         _ _   _ __ _ __  _    _      
//                 U = R R R R = PQPQ PQPQ PQPQ PQPQ

// where each R is a single GARP-1 composite pulse, and the - indicates a
// 180 phase shift.  However, the article leaves things a bit mysterious
// because their listed R (of 25 steps) doesn't coincide with their listing
// of P and Q (even when you combine the first PQ to make 12 steps).  So,
// for the time being I will just use the full 25 steps as the GARP-1 cycle.
// Note that this is also the one in Varian's waveform generator......

//  {
//  pultrain PT;				// Empty pulse train
//  double tp;					// Composite Pi pulse length
//  double phi;					// Composite Pi pulse phase
//  string IsoP;				// Composite Pi isotope channel
//  read_C180(filein,sys,tp,phi,IsoP);		// Read Composite Pi parameters
//  string WF = "GARP-1";			// Set waveform label
//  string CY = "None";				// Set cycle label
//  string SC = "WALTZ-4";			// Set super-cycle label
// sosi needs a serious revamp due to different pulse lengths......
//  row_vector C = C180(tp,phi);		// Get composite pulse waveform
//  PT.pultrain(C, 2*tp, IsoP, WF);		// Set composite pulse waveform
//  Set_Cycle(PT, CY);				// Set for No Cycle
//  Set_Supercycle(PT, SC);			// Set WALTZ-4 Supercycle
//  return PT;	
//  }
                                                                     
// ____________________________________________________________________________
// D                CLASS GARP LIOUVILLE SPACE PROPAGATORS
// ____________________________________________________________________________
 

// ____________________________________________________________________________
// E                        GARP WAVEFORM FUNCTIONS
// ____________________________________________________________________________
 

/* According to the original GARP-1 article
                           _ _   _ __ _ __  _    _      
                   U = R R R R = PQPQ PQPQ PQPQ PQPQ

where each R is a single GARP-1 composite pulse, and the - indicates a
180 phase shift.  However, the article leaves things a bit mysterious
because their listed R (of 25 steps) doesn't coincide with their listing
of P and Q (even when you combine the first PQ to make 12 steps).  So,
for the time being I will just use the full 25 steps as the GARP-1 cycle.
Note that this is also the one in Varian's waveform generator......

The 25 Steps (Pulses) In The GARP-1 Sequence Are As Follows: 

          ____       _____          ____       ____       _____
    30.5  55.2 257.8 268.3  69.3    62.2  85.0 91.8 134.5 256.1 
          ____        ____         _____       ____        ____
    66.4  45.9  25.5  72.7 119.5   138.2 258.4 64.9  70.9  77.2 
         _____        ____
    98.2 133.6 255.9  65.6  53.4                                             */


PulWaveform GARP::WF( ) const { return WF_GARP(); }
PulWaveform GARP::WF_GARP() const

        // Input        GP	: GARP parameters
        // Output       PWF	: Pulse waveform for GARP-1
        // Note                 : Constant rf-amplitude each step
	// Note			: Zero gamB1 produces no decoupling

  {
  double phibar = phi+180.0;			// Phase of "barred" steps
  row_vector WFsteps(25);			// Vector for waveform
  WFsteps.put(complex(gamB1,phi),    0);	// Set waveform values
  WFsteps.put(complex(gamB1,phibar), 1); 	//  { gamB1, phi }
  WFsteps.put(complex(gamB1,phi),    2);
  WFsteps.put(complex(gamB1,phibar), 3);
  WFsteps.put(complex(gamB1,phi),    4);
  WFsteps.put(complex(gamB1,phibar), 5);
  WFsteps.put(complex(gamB1,phi),    6);
  WFsteps.put(complex(gamB1,phibar), 7);
  WFsteps.put(complex(gamB1,phi),    8);
  WFsteps.put(complex(gamB1,phibar), 9);
  WFsteps.put(complex(gamB1,phi),   10);
  WFsteps.put(complex(gamB1,phibar),11);
  WFsteps.put(complex(gamB1,phi),   12);
  WFsteps.put(complex(gamB1,phibar),13);
  WFsteps.put(complex(gamB1,phi),   14);
  WFsteps.put(complex(gamB1,phibar),15);
  WFsteps.put(complex(gamB1,phi),   16);
  WFsteps.put(complex(gamB1,phibar),17);
  WFsteps.put(complex(gamB1,phi),   18);
  WFsteps.put(complex(gamB1,phibar),19);
  WFsteps.put(complex(gamB1,phi),   20);
  WFsteps.put(complex(gamB1,phibar),21);
  WFsteps.put(complex(gamB1,phi),   22);
  WFsteps.put(complex(gamB1,phibar),23);
  WFsteps.put(complex(gamB1,phi),   24);

  row_vector WFtimes(25);			// Vector for step times
  double tdegree = 0; 				// Increment time per pulse
  if(gamB1>0) tdegree = 1/(gamB1*360); 		// degree
  WFtimes.put(30.5*tdegree,  0);
  WFtimes.put(55.2*tdegree,  1);
  WFtimes.put(257.8*tdegree, 2);
  WFtimes.put(268.3*tdegree, 3);
  WFtimes.put(69.3*tdegree,  4);
  WFtimes.put(62.2*tdegree,  5);
  WFtimes.put(85.0*tdegree,  6);
  WFtimes.put(91.8*tdegree,  7);
  WFtimes.put(134.5*tdegree, 8);
  WFtimes.put(256.1*tdegree, 9);
  WFtimes.put(66.4*tdegree,  10);
  WFtimes.put(45.9*tdegree,  11);
  WFtimes.put(25.5*tdegree,  12);
  WFtimes.put(72.7*tdegree,  13);
  WFtimes.put(119.5*tdegree, 14);
  WFtimes.put(138.2*tdegree, 15);
  WFtimes.put(258.4*tdegree, 16);
  WFtimes.put(64.9*tdegree,  17);
  WFtimes.put(70.9*tdegree,  18);
  WFtimes.put(77.2*tdegree,  19);
  WFtimes.put(98.2*tdegree,  20);
  WFtimes.put(133.6*tdegree, 21);
  WFtimes.put(255.9*tdegree, 22);
  WFtimes.put(65.6*tdegree,  23);
  WFtimes.put(53.4*tdegree,  24);
  return PulWaveform(WFsteps, WFtimes, "GARP");
  }

// ____________________________________________________________________________
// F                     GARP COMPOSITE PULSE FUNCTIONS
// ____________________________________________________________________________
 
/* According to the original GARP-1 article
                           _ _   _ __ _ __  _    _      
                   U = R R R R = PQPQ PQPQ PQPQ PQPQ

where each R is a single GARP-1 composite pulse, and the - indicates a
180 phase shift.                                                             */


PulComposite GARP::PCmp(const spin_system& sys) const {return PCmpGARP(sys);}
PulComposite GARP::PCmpGARP(const spin_system& sys) const

        // Input        GP	: GARP parameters
	//		sys	: A spin system
        // Output       CP	: Composite pulse for GARP-1
	//			  applicable to sys

  { return PulComposite(WF(), sys, Iso); }

 
PulComposite GARP::PCmp(const spin_system& sys, const super_op& LOp) const
 
        // Input        GP	: GARP parameters
        //              sys     : A spin system
        //              LOp     : Relaxation/exchange superoperator
        // Output       CP	: Composite pulse for GARP-1
	//			  applicable to sys, relaxation active
 
  { return PulComposite(WF(), sys, LOp, Iso); }
 

// ____________________________________________________________________________
// G                      GARP PULSE CYCLE FUNCTIONS
// ____________________________________________________________________________


/******************************************************************************

   According to the original GARP-1 article
                           _ _   _ __ _ __  _    _      
                   U = R R R R = PQPQ PQPQ PQPQ PQPQ

where each R is a single GARP-1 composite pulse, and the - indicates a
180 phase shift.  Indeed, this is just the cycle used in WALTZ-4

******************************************************************************/

PulCycle GARP::CycGARP1(const spin_system& sys) const

        // Input        GP	: GARP parameters
	//		sys	: A spin system
        // Output       CP	: Pulse cycle for GARP-1 is returned
	// Note			: Uses definition defined in WALTZ

  { return PulCycle(PCmp(sys), CYC_WALTZ4(), "GARP-1"); }


// ____________________________________________________________________________
// H                    GARP PULSE SUPERCYCLE FUNCTIONS
// ____________________________________________________________________________


// ____________________________________________________________________________
// I                        GARP PULSE TRAIN FUNCTIONS
// ____________________________________________________________________________


/* According to the original GARP-1 article
                           _ _   _ __ _ __  _    _      
                   U = R R R R = PQPQ PQPQ PQPQ PQPQ

where each R is a single GARP-1 composite pulse, and the - indicates a
180 phase shift.  Indeed, this is just the cycle used in WALTZ-4             */


/*
PulTrain GARP::PT(const spin_system& sys) const { return PT_GARP1(sys); }
PulTrain GARP::PT_GARP1(const spin_system& sys) const

        // Input        GP	: GARP parameters
	//		sys	: A spin system
        // Output       PT	: Pulse train for GARP-1
	//			  applicable to sys

  { return PulTrain(CP(sys), CYC(), "GARP-1"); }

*/

// ____________________________________________________________________________
// J                    CLASS GARP INTERACITVE FUNCTIONS
// ____________________________________________________________________________


/*
void GARP1_ask(int argc, char* argv[], int& qn, double& gamB1, string& IsoG)

	// Input	argc	: No. arguments
        //              argv    : Argument strings
        //              qn      : Query number
        //              gamB1   : RF-Field strength (kHz)
	//		IsoG	: Decoupling channel
	// None		void    : The rf-field strength for a GARP sequence
	//			  is interactively requested unless supplied
	//			  in the array argv.

  {
  query_parameter(argc, argv, qn++, 
    "\n\tGARP Decoupling Field Strength (KHz)? ",gamB1);
  gamB1 *= 1.e3;
  query_parameter(argc, argv, qn++, 
    "\n\tGARP Decoupling Channel (1H, 13C, ...)? ", IsoG);
  return;
  }


void read_GARP1(string& filein, spin_sys& sys, 
                                         double& gamB1, string& IsoG, int idx)

	// Input	filein  : Input parameter file name
        //              sys     : Active spin system
        //              gamB1   : GARP RF-Field strength (Hz)
	//		IsoG	: GARP Isotope channel
	//		idx     : Parameter name qualifier
	// None		void    : The 2 GARP decoupling parameters are set
	//			  from the parameter values specified in 
	//			  the input parameter file.

  {
  ParameterSet pset;                         // A parameter set
  pset.read(filein);                            // Read pset in
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  string pname, ssfile, pstate;                 // Items in each pset entry
  string SI = string("(") + Gdec(idx)            // Name adjustment if indexed
            + string(")");

//               Determine the GARP Decoupling Isotope Channel

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
      cout << "\n\tCan't Read GARP Isotope Channel Parameter "
          << pname << "!\n";
      exit(-1);
      }
    }

//             Determine the Field Strength in the GARP Sequence
 
  pname = string("gB1GARP1");			// RF Field strength of GARP
  if(idx >= 0) pname += SI;                     // Adjust name if indexed
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())
    (*item).parse(pname,gamB1,pstate);	// Set the GARP field strength
  else
    {
    cout << "\n\tCan't Read GARP Field Strength Parameter "
          << pname << "!\n";
    exit(-1);
    }
  }
*/

// ____________________________________________________________________________
// K                      CLASS GARP INPUT FUNCTIONS
// ____________________________________________________________________________


void GARP::read(const string &filename, int idx)

        // Intput               GP      : GARP parameters
        //                      idx     : GARP index
        // Output               none    : GARP parameters are read in
        //                                from parameters in file filename
        //                                with index idx

  {
  ParameterSet pset;                // Declare a parameter set
  if(!pset.read(filename, 1))          // Read in pset from file
    {
    GARPerror(40, filename);           // Filename problems
    GARPfatality(21);                  // Fatal error
    }
  read(pset, idx);                      // User overloaded function
  return;
  } 


void GARP::read(const ParameterSet& pset, int idx)

        // Intput               GP      : GARP parameters
        //                      pset    : Parameter set
        //                      idx     : GARP index
        // Output               none    : GARP parameters are read in
        //                                from parameters in pset
        //                                with index idx

  {
//              First We'll Set Up The GARP Pulse Parameters
// For The Pulse We Need One of Two: {RF Strength, Composite Pulse Length}

// Note: Reading order is strength --> length for GARP.

  int G = SetGamB1(pset, idx);                  // If strength set, set pulse
  if(!G) GARPfatality(2); 			// Quit if no strength set
  SetChannel(pset,idx);                         // Set GARP pulse channel
  SetPhase(pset, idx);                          // Set GARP pulse phase
  }
    
 
void GARP::ask_read(int argc, char* argv[], int argn, int idx)
 
        // Intput               GP      : GARP parameters
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
        //                      idx     : GARP index
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the GARP parameters
        //                                are read
        //                                If the argument argn is not in argv,
        //                                the user is asked for a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : GARP parameters are modifed (filled)
 
  {
  string filename;                            // Name of spin system file
  query_parameter(argc, argv, argn,           // Get filename from command
    "\n\tGARP parameters filename? ", filename); // Or ask for it
  read(filename, idx);				// Read system from filename
  }
 

 
// ____________________________________________________________________________
// L                      CLASS GARP I/O FUNCTIONS
// ____________________________________________________________________________


ostream& GARP::printBase(ostream &ostr) const

        // Intput               GP      : GARP parameters
        //                      ostr    : Output stream
        // Output               none    : GARP basic parameters are sent
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


ostream& GARP::print(ostream &ostr) const

        // Intput               GP      : GARP parameters
        //                      ostr    : Output stream
        // Output               none    : GARP parameters are sent
        //                                to the output stream

  {                                                            
  ostr << "\n\n\t\t\tGARP Parameters\n";
  printBase(ostr);
  ostr <<"\n\n";
  return ostr;
  }
 
 
 
ostream &operator << (ostream &ostr, const GARP &GP)

        // Intput               GP      : GARP parameters
        //                      ostr    : Output stream
        // Output               none    : GARP parameters are sent
        //                                to the output stream

  {                                                            
  GP.print(ostr);
  return ostr;
  }
 
 
#endif						// PulW_GARP.cc
