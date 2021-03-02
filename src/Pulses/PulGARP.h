/* PulGARP.h **************************************************-*-c++-*-
**			         					**
** 	                        G A M M A				**
**									**
**	GARP Pulse Waveform Functions			Interface	**
**									**
**	Copyright (c) 1998			 			**
**	Dr. Scott A. Smith			 			**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
**  This file contains the implementation of GARP pulse trains in the	**
**  GAMMA simulation platform.  GARP is used as a broad-band decoupling	**
**  sequence.  For details see						**
**									**
** J. Magn. Reson., 64, 547-552 (1985), A.J. Shaka, P.B. Barker, &	**
** Ray Freeman "Computer Optimized Decoupling Scheme for Wideband	**
** Applications and Low-Level Operation"				**
**									**
** The 25 Steps (Pulses) In The GARP-1 Sequence Are As Follows:		**
**									**
**           ____       _____          ____       ____       _____	**
**     30.5  55.2 257.8 268.3  69.3    62.2  85.0 91.8 134.5 256.1 	**
**           ____        ____         _____       ____        ____	**
**     66.4  45.9  25.5  72.7 119.5   138.2 258.4 64.9  70.9  77.2 	**
**          _____        ____						**
**     98.2 133.6 255.9  65.6  53.4					**
**									**
** The above waveform is the repeated in the WALTZ-4 SuperCycle		**
**                                    _ _				**
**                                R R R R				**
**									**
*************************************************************************/

#ifndef   GGARP_h_			// Is this file already included?
#  define GGARP_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Pulses/PulWaveform.h>		// Know about pulse waveforms
#include <Pulses/PulComposite.h>	// Know about composite pulses
#include <Pulses/PulCycle.h>		// Know about pulse cycles
#include <Basics/ParamSet.h>		// Know about parameter sets 
#include <string>			// Know about stdlibc++ strings

class GARP
  {
  std::string Iso;                           // Applied pulse channel
  double gamB1;                         // Applied pulse strength (Hz)
  double phi;                           // Applied pulse phase    (deg)
  double Wrf;                           // Applied pulse offset   (Hz)

private:
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                       CLASS GARP ERROR HANDLING
// ____________________________________________________________________________
 
 
void GARPerror(int eidx, int noret=0) const;

        // Input                GARP	: GARP parameters(this)
        //                      eidx    : Error flag
        //                      noret   : Return flag
        // Output               none    : Error Message Output

                                                               
void volatile GARPfatality(int error) const;

        // Input                GARP	: GARP parameters(this)
        //                      eidx    : Error flag
        // Output               none    : Stops execution & error Message
 

void GARPerror(int eidx, const std::string& pname, int noret=0) const;
 
        // Input                GARP	: GARP parameters(this)
        //                      eidx    : Error index
        //                      pname   : String included in message
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
                                                         
// ____________________________________________________________________________
// ii               CLASS GARP PARAMETER SET FUNCTIONS
// ____________________________________________________________________________


void SetPhase(const ParameterSet& pset, int idx=-1);

        // Intput               DT      : GARP parameters
        //                      pset    : Parameter set
        //                      idx     : GARP index
        // Output               none    : GARP pulse phase read in
        //                                from parameter in pset
        //                                with index idx


void SetChannel(const ParameterSet& pset, int idx=-1);

        // Intput               DT      : GARP parameters
        //                      pset    : Parameter set
        //                      idx     : GARP index
        // Output               none    : GARP pulse channel read in
        //                                from parameter in pset
        //                                with index idx
 
 
int SetGamB1(const ParameterSet& pset, int idx=-1);
 
        // Intput               DT      : GARP parameters
        //                      pset    : Parameter set
        //                      idx     : GARP index
        // Output               TF      : GARP pulse strength read in
        //                                from parameter in pset
        //                                with index idx


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 public:

// ____________________________________________________________________________
// A                   CLASS GARP CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________


MSVCDLC GARP();

        // Input        none    :
        // Output       GARP	: GARP parameters(this)
        ///F_list       GARP	- Constructor

                                              
MSVCDLC GARP(double gB1, const std::string& ch, double ph=0, double off=0);

        // Input        gB1     : RF-field strength (Hz)
	//		ch	: RF-channel
        //              ph      : RF-phase (degrees)
        //              off     : RF-offset (Hz) 
        // Output       GARP    : GARP parameters(this)


MSVCDLC GARP(const GARP& PT1);

        // Intput       GARP1	: GARP parameters
        // Output       GARP	: GARP parameters(this) from GARP1

                                                                     
// ------------------------------ Destruction ---------------------------------
         
         
MSVCDLC ~GARP();

        // Intput       GARP1	: GARP parameters (this)
        // Output       none    : GARP is destructed

                                                      
// ------------------------------- Assignment ---------------------------------
 

MSVCDLL GARP& operator = (const GARP& GARP1);

        // Intput       GARP1	: GARP parameters
        // Output       GARP	: GARP parameters(this) from GARP1
        ///F_list       =       - Assignment


// ____________________________________________________________________________
// B                     CLASS GARP ACCESS FUNCTIONS
// ____________________________________________________________________________

MSVCDLL std::string channel()  const;
MSVCDLL double strength() const;
MSVCDLL void   strength(double gB1);
MSVCDLL double phase()    const;
MSVCDLL double offset()   const;

        // Intput       GARP1  : GARP parameters
        // Output 	channel : GARP isotope channel
        //              strength: GARP pulse strength (Hz)
        //              phase   : GARP pulse phase (sec)
        //              offset  : GARP pulse offset (Hz)


// ____________________________________________________________________________
// C                  CLASS GARP HILBERT SPACE PROPAGATORS
// ____________________________________________________________________________


// ____________________________________________________________________________
// D                CLASS GARP LIOUVILLE SPACE PROPAGATORS
// ____________________________________________________________________________


// ____________________________________________________________________________
// E                       GARP WAVEFORM FUNCTIONS
// ____________________________________________________________________________
 

MSVCDLL PulWaveform WF( ) const;
MSVCDLL PulWaveform WF_GARP( ) const;

        // Input        G       : GARP specifications
        // Output       WF      : A pulse waveform for GARP pulse-delay

                                                                         
// ____________________________________________________________________________
// F                    GARP COMPOSITE PULSE FUNCTIONS
// ____________________________________________________________________________
 
 
MSVCDLL PulComposite PCmp(const spin_system& sys) const;
MSVCDLL PulComposite PCmpGARP(const spin_system& sys) const;

        // Input        GP      : GARP parameters
        //              sys     : A spin system
        // Output       CP      : Composite pulse for GARP-1
        //                        applicable to sys
 
 
MSVCDLL PulComposite PCmp(const spin_system& sys, const super_op& LOp) const;
  
        // Input        GP      : GARP parameters
        //              sys     : A spin system 
        //              LOp     : Relaxation/exchange superoperator 
        // Output       CP      : Composite pulse for GARP-1
        //                        applicable to sys, relaxation active

// ____________________________________________________________________________
// G                      GARP PULSE CYCLE FUNCTIONS
// ____________________________________________________________________________

/* According to the original GARP-1 article
                           _ _   _ __ _ __  _    _      
                   U = R R R R = PQPQ PQPQ PQPQ PQPQ

where each R is a single GARP-1 composite pulse, and the - indicates a
180 phase shift.  Indeed, this is just the cycle used in WALTZ-4             */

MSVCDLL PulCycle CycGARP1(const spin_system& sys) const;

        // Input        GP      : GARP parameters
        //              sys     : A spin system
        // Output       CYC     : Pulse cycle for GARP-1 is returned
        // Note                 : Uses definition defined in WALTZ


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
 
 
//PulTrain PT(const spin_system& sys) const;
//PulTrain PT_GARP1(const spin_system& sys) const;
 
        // Input        GP      : GARP parameters
        //              sys     : A spin system
        // Output       PT      : Pulse train for GARP-1
        //                        applicable to sys

// ____________________________________________________________________________
// J                     GARP INTERACTIVE SETUP FUNCTIONS
// ____________________________________________________________________________


//void GARP1_ask(int argc, char* argv[], int& qn, double& gamB1, std::string& IsoG);

	// Input	argc	: No. arguments
        //              argv    : Argument strings
        //              qn      : Query number
        //              gamB1   : RF-Field strength (kHz)
	//		IsoG	: Decoupling channel
	// None		void    : The rf-field strength for a GARP sequence
	//			  is interactively requested unless supplied
	//			  in the array argv.


//void read_GARP1(string& filein, spin_sys& sys, 
//                                    double& gamB1, std::string& IsoG, int idx=-1);

	// Input	filein  : Input parameter file name
        //              sys     : Active spin system
        //              gamB1   : GARP RF-Field strength (Hz)
	//		IsoG	: GARP Isotope channel
	//		idx     : Parameter name qualifier
	// None		void    : The 2 GARP decoupling parameters are set
	//			  from the parameter values specified in 
	//			  the input parameter file.


//pultrain GARP1(std::string& filein, spin_sys& sys, int idx=-1);

        // Input        filein  : Input parameter file name
        //              sys     : Active spin system
        //              idx     : Parameter name qualifier
        // None         PT      : A pulse train is returned set
        //                        to GARP-1 for the system sys
        //                        based on the parameters specified
        //                        in the input parameter file.
        //                      : This pulse train has WALTZ-4 as the
        //                        supercycle, no cycle, and the sequence
        //                        set to a composite 25 step GARP-1
 
// According to the original GARP-1 article
//                         _ _   _ __ _ __  _    _
//                 U = R R R R = PQPQ PQPQ PQPQ PQPQ
 
// where each R is a single GARP-1 composite pulse, and the - indicates a
// 180 phase shift.  However, the article leaves things a bit mysterious
// because their listed R (of 25 steps) doesn't coincide with their listing
// of P and Q (even when you combine the first PQ to make 12 steps).  So,
// for the time being I will just use the full 25 steps as the GARP-1 cycle.
// Note that this is also the one in Varian's waveform generator.....
 
// ____________________________________________________________________________
// K                      CLASS GARP INPUT FUNCTIONS
// ____________________________________________________________________________


MSVCDLL void read(const std::string &filename, int idx=-1);

        // Intput               GP      : GARP parameters
        //                      idx     : GARP index
        // Output               none    : GARP parameters are read in
        //                                from parameters in file filename
        //                                with index idx


MSVCDLL void read(const ParameterSet& pset, int idx=-1);

        // Intput               GP      : GARP parameters
        //                      pset    : Parameter set
        //                      idx     : GARP index
        // Output               none    : GARP parameters are read in
        //                                from parameters in pset
        //                                with index idx

//                  Read GARP Delay & Pulse Parameters
// For The Pulse We Need Two of Three: {Pulse Angle, RF Strength, Pulse Length}

// Note: Reading order is angle --> length --> strength for GARP.  If only
//       the angle is specified the length will be set to zero (ideal pulse)

 
MSVCDLL void ask_read(int argc, char* argv[], int argn, int idx=-1);

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

 
// ____________________________________________________________________________
// L                      CLASS GARP I/O FUNCTIONS
// ____________________________________________________________________________
 

MSVCDLL std::ostream& printBase(std::ostream &ostr) const;
 
        // Intput               GP      : GARP parameters
        //                      ostr    : Output stream
        // Output               none    : GARP basic parameters are sent
        //                                to the output stream
 
 
MSVCDLL std::ostream& print(std::ostream &ostr) const;
 
        // Intput               GP      : GARP parameters
        //                      ostr    : Output stream
        // Output               none    : GARP parameters are sent
        //                                to the output stream

                                                               
MSVCDLL friend std::ostream &operator << (std::ostream &ostr, const GARP &GP);

        // Intput               GP      : GARP parameters    
        //                      ostr    : Output stream
        // Output               none    : GARP parameters are sent
        //                                to the output stream
  };

#endif						// PulW_GARP.h
