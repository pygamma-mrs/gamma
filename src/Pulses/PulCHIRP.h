/* PulCHIRP.h ***************************************************-*-c++-*-
**			         					**
** 	                        G A M M A				**
**									**
**	CHIRP Pulse Functions 		               Interface	**
**									**
**	Copyright (c) 1998			 			**
**	Dr. Scott A. Smith			 			**
**      Philippe Pelupessy						**
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
** This file contains the implementation of CHIRP pulses and pulse	**
** trains in the GAMMA simulation platform.  CHIRP is used as a 	**
** broad-band decoupling sequence.  For details see 			**
**                                                                      **
** The above waveform is the repeated in the Tyko-Pines Cycle		**
**                                    _ _                               **
**                                R R R R                               **
**                                                                      **
*************************************************************************/

#ifndef   GPulCHIRP_h_			// Is this file already included?
#  define GPulCHIRP_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Pulses/PulWaveform.h>		// Know about pulse waveforms
#include <Pulses/PulComposite.h>	// Know about composite pulses
#include <Pulses/PulCycle.h>		// Know about pulse cycles
 
// ____________________________________________________________________________
// A                         CHIRP-95 WAVEFORM FUNCTIONS
// ____________________________________________________________________________


MSVCDLL PulWaveform WF_CHIRP95(int N, double tp, double delW, double gB1, int scale=0);

        // Input        N       : Number of pulse steps
        //              tp      : Pulse length (sec)
        //              delW    : Sweep width (Hz)
        //              gB1	: Field strength scaling (default 1)
        //              scale   : Flag whether to scale the phases
        //                            0 = No scaling (default)
        //                           !0 = phase is [-pi, pi]
        // Output       PWF	: Chirp waveform
        // Note                 : This uses constant increment time

 
// ____________________________________________________________________________
// B                      CHIRP COMPOSITE PULSE FUNCTIONS
// ____________________________________________________________________________
 
 
MSVCDLL PulComposite CP_CHIRP95(const spin_system& sys, const std::string& IsoC,
                       int N, double tp, double delW, double gB1, int scale=0);
 
        // Input        sys     : Active spin system
        //              IsoC    : Pulse channel
        // 		N       : Number of pulse steps
        //              tp      : Pulse length (sec)
        //              delW    : Sweep width (Hz)
        //              gB1	: Field strength scaling (default 1)
        //              scale   : Flag whether to scale the phases
        //                            0 = No scaling (default)
        // None         PCom    : A composite pulse is returned set
        //                        to CHIRP95 for the system sys

// ____________________________________________________________________________
// C                      CHIRP PULSE CYCLE FUNCTIONS
// ____________________________________________________________________________
  
  
MSVCDLL PulCycle CYC_CHIRP95();
  
        // Input        void    : Not a thing
        // None         PCyc    : Pulse cycle for CHIRP-95 is returned 
        // Note                 : Uses 5 steps cycle as suggested by
        //                        Tyko and Pines in Chem. Phys. Lett.
        //                        111, (1984) 462: {0,150,60,160,0}
  


// ____________________________________________________________________________
//                         CHIRP-95 SETUP FUNCTIONS
// ____________________________________________________________________________

/*

MSVCDLL void CHIRP_ask(int argc, char* argv[], int& qn,
                  int& nsteps, double& tp, double& delW, double& gamB1);

        // Input        argc    : No. arguments
        //              argv    : Argument strings
        //              qn      : Query number
        //              nsteps  : Number of CHIRP steps
        //              tp      : CHIRP length (msec)
        //              delW    : CHIRP sweep width (kHz)
        //              gamB1   : CHIRP RF-Field strength (kHz)
        // None         void    : The 4 CHIRP parameters are interactively
        //                        requested unless supplied in argv. All
        //                        four values are herein set and returned
        //                        in "SI" units (i.e. sec & Hz).


MSVCDLL void read_CHIRP(std::string& filein, spin_sys& sys, int& nsteps,
      double& tp, double& delW, double& gamB1, std::string& IsoCHIRP, int idx=-1);

        // Input        filein  : Input parameter file name
        //              sys     : Active spin system
        //              nsteps  : Number of CHIRP steps
        //              tp      : CHIRP length (sec)
        //              delW    : CHIRP sweep width (Hz)
        //              gamB1   : CHIRP RF-Field strength (Hz)
        //              IsoCHIRP: CHIRP Isotope channel
        //              idx     : Parameter name qualifier
        // None         void    : The 5 CHIRP parameters are set from the
        //                        parameter values specified in the input
        //                        parameter file.


MSVCDLL pultrain CHIRP95(std::string& filein, spin_sys& sys, int idx=-1);

        // Input        filein  : Input parameter file name
        //              sys     : Active spin system
        //              idx     : Parameter name qualifier
        // None         PT      : A pulse sequence is returned set
        //                        to CHIRP-95 for the system sys
        //                        based on the parameters specified
        //                        in the input parameter file.


MSVCDLL void CHIRP95(pultrain& PT, int ns, double tpul, double delW,
                                   double gamB1, std::string& IsoCHIRP, int scl);

	// Input	PT	: Pulse Train
        //              ns      : Number of pulse steps
        //              tpul    : Pulse length (sec)
        //              delW    : Sweep width (Hz)
        //              gamB1   : Field strength scaling (default 1)
	//		IsoCHIRP: Pulse isotope channel
	//		scl	: Flag whether to scale the phases
	//			      0 = No scaling (default)
	//			     !0 = phase is [-pi, pi]
	// None		void  : Pulse Train cycle set to CHIRP-95


// ______________________________________________________________________
//                 PULTRAIN AUXILIARY WAVEFORM FUNCTIONS
// ______________________________________________________________________


MSVCDLL row_vector Chirp(int N, double tp, double delW, double gamB1, int scale);

        // Input        N       : Number of pulse steps
        //              tp      : Pulse length (sec)
        //              delW    : Sweep width (Hz)
        //              gamB1   : Field strength scaling (default 1)
	//		scale   : Flag whether to scale the phases
	//			      0 = No scaling (default)
	//			     !0 = phase is [-pi, pi]
        // Output       Cvect   : Vector of Chirp waveform
        // Note                 : This assumes constant rf-amplitude
	//		        : and increment time

*/

#endif						// pultrain_aux.cc
