/* PulCHIRP.cc **************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      CHIRP Pulse Functions                     Implementation	**
**                                                                      **
**      Copyright (c) 1998                                              **
**      Dr. Scott A. Smith                                              **
**      Philippe Pelupessy                                              **
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
** This file contains the implementation of CHIRP pulses and pulse      **
** trains in the GAMMA simulation platform.  CHIRP is used as a         **
** broad-band decoupling sequence.  For details see                     **
**									**
*************************************************************************/

#ifndef _PulCHIRP_cc_			// Is this file already included?
#define _PulCHIRP_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif
 
#include <Pulses/PulCHIRP.h>		// Include header files
#include <Pulses/PulWaveform.h>		// Know pulse waveforms
#include <Pulses/PulComposite.h>	// Know composite pulses
#include <Pulses/PulCycle.h>		// Know pulse cycles
#include <Matrix/row_vector.h>		// Know about row vectors
#include <Basics/StringCut.h>		// Know about the Gdec function
#include <string>                       // Include libstdc++ strings

using std::string;			// Using libstdc++ strings
// ____________________________________________________________________________
// A                         CHIRP-95 WAVEFORM FUNCTIONS
// ____________________________________________________________________________


PulWaveform WF_CHIRP95(int N, double tp, double delW, double gamB1, int scale)

        // Input        N       : Number of pulse steps
        //              tp      : Pulse length (sec)
        //              delW    : Sweep width (Hz)
        //              gamB1   : Field strength scaling (default 1)
	//		scale   : Flag whether to scale the phases
	//			      0 = No scaling (default)
	//			     !0 = phase is [-pi, pi]
        // Output       Cvect   : Vector of Chirp waveform
        // Note                 : This uses constant rf-amplitude
	//		          and increment time

  {
  row_vector WFsteps(N);			// Waveform vector{ gamB1,phi }
  double tinc = tp/double(N);			// Step increment time
  row_vector WFtimes(5, tinc);			// Vector for step times
  double phase, cutoff=-1.e-6;			// For phase and accuracy 
  double t=0;					// Start at time 0
  for(int i=0; i<N; i++)			// Loop over waveform points
    {
    phase = -0.5*t + t*t/(2.0*tp);
    phase *= PI2*delW;
    if(scale)
      {
      while(phase < -PI2)   phase += PI2;
      if(phase+PI < cutoff) phase += PI2;
      }
    WFsteps.put(complex(gamB1, phase),i);	// Set the value in waveform
    t += tinc;					// Increment step time 
    }
  return PulWaveform(WFsteps, WFtimes, "CHIRP");
  }

// ____________________________________________________________________________
// B                      CHIRP COMPOSITE PULSE FUNCTIONS
// ____________________________________________________________________________


PulComposite CP_CHIRP95(const spin_system& sys, const string& IsoC,
                        int N, double tp, double delW, double gamB1, int scale)
 
        // Input        sys     : Active spin system
        //              IsoC    : Pulse channel
        // 		N       : Number of pulse steps
        //              tp      : Pulse length (sec)
        //              delW    : Sweep width (Hz)
        //              gamB1   : Field strength scaling (default 1)
	//		scale   : Flag whether to scale the phases
	//			      0 = No scaling (default)
	//			     !0 = phase is [-pi, pi]
        // None         PCom    : A composite pulse is returned set
        //                        to CHIRP-95 for the system sys
 
  { return PulComposite(WF_CHIRP95(N, tp, delW, gamB1, scale), sys, IsoC); }
 
 
// ____________________________________________________________________________
// C                      CHIRP PULSE CYCLE FUNCTIONS
// ____________________________________________________________________________
 
 
PulCycle CYC_CHIRP95()
 
        // Input        void    : Not a thing
        // None         PCyc    : Pulse cycle for CHIRP-95 is returned
        // Note                 : Uses 5 steps cycle as suggested by
        // 			  Tyko and Pines in Chem. Phys. Lett.
	//                        111, (1984) 462: {0,150,60,160,0}
 
  {
  row_vector PTCsteps(5);                       // Cycle vector{ phi }
  PTCsteps.put(0.0,   0);			// Set cycle phase values
  PTCsteps.put(150.0, 1);
  PTCsteps.put(60.0,  2);
  PTCsteps.put(150.0, 3);
  PTCsteps.put(0.0,   4);
  return PulCycle();
//  return PulCycle(PTCsteps, "Tyko-Pines");
  }


// ____________________________________________________________________________
//                           CHIRP-95 Setup Functions
// ____________________________________________________________________________

/*
void CHIRP_ask(int argc, char* argv[], int& qn,
                         int& nsteps, double& tp, double& delW, double& gamB1)

	// Input	argc	: No. arguments
        //              argv    : Argument strings
        //              qn      : Query number
	//		nsteps  : Number of CHIRP steps
	//		tp      : CHIRP length (msec)
        //              delW    : CHIRP sweep width (kHz)
        //              gamB1   : CHIRP RF-Field strength (kHz)
	// None		void    : The 4 CHIRP parameters are interactively
	//			  requested unless supplied in argv. All
	//			  four values are herein set and returned
	//			  in "SI" units (i.e. sec & Hz).

  {
  query_parameter(argc, argv, qn++,             // Number chirp steps 
         "\n\tNumber of CHIRP Steps? ",nsteps);
  query_parameter(argc, argv, qn++,
    "\n\tTotal CHIRP Pulse Length (msec)? ",tp);
  tp *= 1.e-3;
  query_parameter(argc, argv, qn++,
    "\n\tTotal CHIRP Sweep Width (KHz)? ", delW);
  delW *= 1.e3;
  query_parameter(argc, argv, qn++,
    "\n\tCHIRP RF-Field Strength (KHz)? ", gamB1);
  gamB1 *= 1.e3;
  return;
  }


void read_CHIRP(string& filein, spin_sys& sys, int& nsteps,
      double& tp, double& delW, double& gamB1, string& IsoCHIRP, int idx)

	// Input	filein  : Input parameter file name
        //              sys     : Active spin system
	//		nsteps  : Number of CHIRP steps
	//		tp      : CHIRP length (sec)
        //              delW    : CHIRP sweep width (Hz)
        //              gamB1   : CHIRP RF-Field strength (Hz)
	//		IsoCHIRP: CHIRP Isotope channel
	//		idx     : Parameter name qualifier
	// None		void    : The 5 CHIRP parameters are set from the
	//			  parameter values specified in the input
	//			  parameter file.

  {
  ParameterAVLSet pset;                         // A parameter set
  pset.read(filein);                            // Read pset in
  Pix item;                                     // A pix into the pset
  SinglePar par;
  string pname, ssfile, pstate;                 // Items in each pset entry
  string SI = string("(") + Gdec(idx)		// Name adjustment if indexed
            + string(")");

//               Determine the CHIRP Decoupling Isotope Channel

  if(!sys.heteronuclear()) IsoCHIRP = sys.symbol(0);
  else
    {
    pname = string("IsoCHIRP");			// Chirp Isotope Parameter
    if(idx >= 0) pname += SI;			// Adjust name if indexed
    par = SinglePar(pname);
    item = pset.seek(par);			// Pix in parameter list
    if(item)
      (pset(item)).parse(pname,IsoCHIRP,pstate);// Set Detected Isotope
    else
      {
      cout << "\n\tCan't Read CHIRP Isotope Channel Parameter "
          << pname << "!\n";
      exit(-1);
      }
    }

//             Determine the Number of CHIRP Steps
 
  pname = string("CHIRPstps");			// Number of CHIRP steps
  if(idx >= 0) pname += SI;                     // Adjust name if indexed
  par = SinglePar(pname);
  item = pset.seek(par);                        // Pix in parameter list
  if(item)
    (pset(item)).parse(pname,nsteps,pstate);	// Set number of steps
  else
    {
    cout << "\n\tCan't Read Number of CHIRP Steps Parameter "
          << pname << "!\n";
    exit(-1);
    }

//             Determine the Length of the CHIRP Sequence
 
  pname = string("tpCHIRP");			// Length of CHIRP
  if(idx >= 0) pname += SI;                     // Adjust name if indexed
  par = SinglePar(pname);
  item = pset.seek(par);                        // Pix in parameter list
  if(item)
    (pset(item)).parse(pname,tp,pstate);	// Set CHIRP pulse length
  else
    {
    cout << "\n\tCan't Read CHIRP Length Parameter "
          << pname << "!\n";
    exit(-1);
    }

//             Determine the Sweep Width of the CHIRP Sequence
 
  pname = string("swCHIRP");			// Sweep width of CHIRP
  if(idx >= 0) pname += SI;                     // Adjust name if indexed
  par = SinglePar(pname);
  item = pset.seek(par);                        // Pix in parameter list
  if(item)
    (pset(item)).parse(pname,delW,pstate);	// Set the CHIRP sweep width
  else
    {
    cout << "\n\tCan't Read CHIRP Sweep Width Parameter "
          << pname << "!\n";
    exit(-1);
    }

//             Determine the Field Strength in the CHIRP Sequence
 
  pname = string("gB1CHIRP");			// RF Field strength of CHIRP
  if(idx >= 0) pname += SI;                     // Adjust name if indexed
  par = SinglePar(pname);
  item = pset.seek(par);                        // Pix in parameter list
  if(item)
    (pset(item)).parse(pname,gamB1,pstate);	// Set the CHIRP field strength
  else
    {
    cout << "\n\tCan't Read CHIRP Field Strength Parameter "
          << pname << "!\n";
    exit(-1);
    }
  }


pultrain CHIRP95(string& filein, spin_sys& sys, int idx)

	// Input	filein  : Input parameter file name
        //              sys     : Active spin system
	//		idx     : Parameter name qualifier
	// None		PT      : A pulse sequence is returned set
	//			  to CHIRP-95 for the system sys
	//			  based on the parameters specified
	//			  in the input parameter file.

  {
  pultrain PT;						// Empty pulse train
  int ns;						// Chirp pulse steps
  double tp;                                            // Chirp pulse length
  double delW;                                          // Chirp freq. width
  double gamB1;                                         // Chirp field strength
  string IsoP;                                          // Chirp isotope channel
  read_CHIRP(filein,sys,ns,tp,delW,gamB1,IsoP);		// Read Chirp parameters
  string WF = "Chirp-95";				// Set waveform label
  string CY = "Tyko-Pines";				// Set cycle label
  string SC = "MLEV-16";				// Set super-cycle label
  row_vector C = Chirp(ns,tp,delW,gamB1,1);		// Construct Chirp waveform
  PT.pultrain(C, tp, IsoP, WF, 1);			// Set Chirp Waveform
  Set_Cycle(PT, CY);					// Set Tyko-Pines Cycle
  Set_Supercycle(PT, SC);				// Set MLEV-16 Supercycle
  return PT;	
  }


void CHIRP95(pultrain& PT, int ns, double tpul, double delW,
                                  double gamB1, string& IsoCHIRP, int scl)

	// Input	PT	: Pulse Train
        //              ns      : Number of pulse steps
        //              tpul    : Pulse length (sec)
        //              delW    : Sweep width (Hz)
        //              gamB1   : Field strength scaling (default 1)
	//		IsoCHIRP: CHIRP Isotope channel
	//		scl	: Flag whether to scale the phases
	//			      0 = No scaling (default)
	//			     !0 = phase is [-pi, pi]
	// None		void  : Pulse Train cycle set to CHIRP-95

  {
  string WF = "Chirp-95";			// Set waveform label
  string CY = "Tyko-Pines";			// Set cycle label
  string SC = "MLEV-16";			// Set super-cycle label
  row_vector C = Chirp(ns,tpul,delW,gamB1,scl);	// Construct Chirp waveform
  PT.pultrain(C, tpul, IsoCHIRP, WF, 1);	// Set Chirp-95 Waveform
  Set_Cycle(PT, CY);				// Set Tyko-Pines Cycle
  Set_Supercycle(PT, SC);			// Set MLEV-16 Supercycle
  }
*/

#endif						// PulCHIRP.cc
