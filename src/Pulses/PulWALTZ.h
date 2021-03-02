/* PulWALTZ.h ***************************************************-*-c++-*-
**			         					**
** 	                        G A M M A				**
**									**
**	WALTZ Pulse Functions				Interface	**
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
**                                                                      **
**  This file contains the a variety of functions supporting the use    **
**  of WALTZ pulse trains in the GAMMA simulation platform.  The basic  **
**  WALTZ waveforms are composite pulses and these are cycled to build  **
**  up longer WALTZ based pulse trains.  WALTZ is primarily used as a   **
**  broad-band decoupling sequence.  For details see                    **
**                                                                      **
** See J. Magn. Reson., 53, 313-340 (1983), A.J. Shaka, James Keeler,   **
** and Ray Freeman "Evaluation of a New Broadband Decoupling Sequence:  **
** WALTZ-16"                                                            **
**                                                                      **
** WALTZ waveforms (base composite pulses): WALTZ-R, WALTZ-Q, WALTZ-K   **
** WALTZ pulse trains (cycled waveforms):   WALTZ-4, WALTZ-8, WALTZ-16  **
**                                                                      **
*************************************************************************/

#ifndef   GWALTZ_h_			// Is this file already included?
#  define GWALTZ_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Pulses/Pulse.h>		// Know about pulses
#include <Pulses/PulWaveform.h>		// Know pulse waveforms
#include <Pulses/PulComposite.h>	// Know composite pulses 
#include <Pulses/PulCycle.h>		// Know pulse cycles
#include <Matrix/row_vector.h>		// Know about row vectors

class WALTZ : public Pulse
  {
 
private:

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i                       CLASS WALTZ ERROR HANDLING
// ____________________________________________________________________________

                                                                                
void WALTZerror(int eidx, int noret=0) const;
 
        // Input                WALTZ    : WALTZ parameters(this)
        //                      eidx    : Error flag
        //                      noret   : Return flag
        // Output               none    : Error Message Output
 

void volatile WALTZfatality(int error) const;
 
        // Input                WALTZ    : WALTZ parameters(this)
        //                      eidx    : Error flag
        // Output               none    : Stops execution & error Message

                                                                          
void WALTZerror(int eidx, const std::string& pname, int noret=0) const;

        // Input                WALTZ    : WALTZ parameters(this)   
        //                      eidx    : Error index
        //                      pname   : String included in message
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message

// ____________________________________________________________________________
// ii               CLASS WALTZ PARAMETER SET FUNCTIONS
// ____________________________________________________________________________
 
 
void SetPhase(const ParameterSet& pset, int idx=-1);
 
        // Intput               DT      : WALTZ parameters
        //                      pset    : Parameter set
        //                      idx     : WALTZ index
        // Output               none    : WALTZ pulse phase read in
        //                                from parameter in pset
        //                                with index idx
 
 
void SetChannel(const ParameterSet& pset, int idx=-1);
         
        // Intput               DT      : WALTZ parameters
        //                      pset    : Parameter set
        //                      idx     : WALTZ index
        // Output               none    : WALTZ pulse channel read in
        //                                from parameter in pset
        //                                with index idx

                                                         
int SetGamB1(const ParameterSet& pset, int idx=-1);

        // Intput               DT      : WALTZ parameters
        //                      pset    : Parameter set
        //                      idx     : WALTZ index
        // Output               TF      : WALTZ pulse strength read in
        //                                from parameter in pset
        //                                with index idx

         
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
 public:
 
// ____________________________________________________________________________
// A                   CLASS WALTZ CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 
 
MSVCDLC WALTZ();
 
        // Input        none	:
        // Output       WALTZ	: WALTZ parameters(this)
        ///F_list       WALTZ	- Constructor
 
 
MSVCDLC WALTZ(double gB1, const std::string& ch, double ph=0, double off=0);
 
        // Input        gB1     : RF-field strength (Hz)
        //              ch      : RF-channel
        //              ph      : RF-phase (degrees)
        //              off     : RF-offset (Hz)
        // Output       WALTZ	: WALTZ parameters(this)
 
 
MSVCDLC WALTZ(const WALTZ& WP1);
 
        // Intput       WALTZ1   : WALTZ parameters
        // Output       WALTZ    : WALTZ parameters(this) from WALTZ1

         
// ------------------------------ Destruction ---------------------------------
         
 
MSVCDLC ~WALTZ();
 
        // Intput       WALTZ1   : WALTZ parameters (this)
        // Output       none    : WALTZ is destructed
 

// ------------------------------- Assignment ---------------------------------

                                                                                
MSVCDLL WALTZ& operator = (const WALTZ& WALTZ1);
 
        // Intput       WALTZ1   : WALTZ parameters
        // Output       WALTZ    : WALTZ parameters(this) from WALTZ1
        ///F_list       =       - Assignment
 
 
// ____________________________________________________________________________
// B                     CLASS WALTZ ACCESS FUNCTIONS
// ____________________________________________________________________________
 
/*
std::string channel()  const;			INHERITED
double strength() const;			INHERITED
double phase()    const;			INHERITED
double offset()   const;			INHERITED                    */
 
        // Intput       WALTZ1  : WALTZ parameters
        // Output       channel : WALTZ isotope channel
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

/******************************************************************************

     There are 3 steps in a WALTZ-R sequence, as shown below.
                                  _                _
                        U = U(90)*U(180)*U(270) = 123

     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift.

******************************************************************************/

MSVCDLL PulWaveform WF(int even=0) const;
MSVCDLL PulWaveform WF_WALTZR(int even=0) const;
 
        // Input        WP      : WALTZ parameters
	//		even	: Flag if even step sizes
        // Output       PWF     : Pulse waveform for WALTZ-1
        // Note                 : Constant rf-amplitude each step
	// Note			: Ideal pulses are allowed
 

/******************************************************************************

     There are 5 steps in a WALTZ-K sequence, as shown below.
                  _             _             _       _ _ _
              U = U(180)*U(360)*U(180)*U(270)*U(90) = 24231

     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift.
 
******************************************************************************/

MSVCDLL PulWaveform WF_WALTZK(int even=0) const;

        // Input        gamB1   : RF field strength (Hz)
	//		even	: Flag if even step sizes
        // Output       PWF     : Composite pulse waveform
        // Note                 : Uses constant rf-amplitude.


/* ****************************************************************************

     There are 9 steps in a WALTZ-Q sequence, as shown below.
      _             _             _            _             _        _ _ _ _ _
  U = U(270)*U(360)*U(180)*U(270)*U(90)*U(180)*U(360)*U(180)*U(270) = 342312423

     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift.

******************************************************************************/
 
MSVCDLL PulWaveform WF_WALTZQ(int even=0) const;

        // Input        tp      : Single Pi Pulse length (sec)
	//		even	: Flag if even step sizes
        // Output       Cvect   : Vector of composite pulse waveform
        // Note                 : Uses constant rf-amplitude and increment
        //                        time, thus is constructed as 4 steps


// ____________________________________________________________________________
// F                    WALTZ COMPOSITE PULSE FUNCTIONS
// ____________________________________________________________________________


/* **************************************************************************** 
 
     There are 3 steps in a WALTZ-R sequence, as shown below. 
                                  _                _
                        U = U(90)*U(180)*U(270) = 123 
 
     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift. 
 
******************************************************************************/
 
 
MSVCDLL PulComposite PCmp(const spin_system& sys, int even=0) const;
MSVCDLL PulComposite PCmpWALTZR(const spin_system& sys, int even=0) const;
 
        // Input        WP      : WALTZ parameters
        //              sys     : A spin system
	//		even	: Flag if even step sizes
        // Output       CP      : WALTZ-R composite pulse for sys


/* ****************************************************************************
 
     There are 5 steps in a WALTZ-K sequence, as shown below.
                  _             _             _       _ _ _
              U = U(180)*U(360)*U(180)*U(270)*U(90) = 24231
 
     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift.
 
******************************************************************************/
 

MSVCDLL PulComposite PCmpWALTZK(const spin_system& sys, int even=0) const;
 
        // Input        WP      : WALTZ parameters
        //              sys     : A spin system
	//		even	: Flag if even step sizes
        // None         CP      : A composite pulse is returned set
        //                        to WALTZ-K for the system sys
        //                        of strength gamB1 and phase phi
 
                                                  
/* ****************************************************************************
 
     There are 9 steps in a WALTZ-Q sequence, as shown below.
      _             _             _            _             _        _ _ _ _ _
  U = U(270)*U(360)*U(180)*U(270)*U(90)*U(180)*U(360)*U(180)*U(270) = 342312423
 
     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift.
 
******************************************************************************/
 
MSVCDLL PulComposite PCmpWALTZQ(const spin_system& sys, int even=0) const;
 
        // Input        WP      : WALTZ parameters
        //              sys     : A spin system
	//		even	: Flag if even step sizes
        // None         CP      : A composite pulse is returned set
        //                        to WALTZ-Q for the system sys
        //                        of strength gamB1 and phase phi

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


MSVCDLL PulCycle CycWALTZ4(const spin_system& sys, int even=0) const;

        // Input        void    : None
        //              phi     : Phase angle (degrees)
	//		even	: Flag if even step sizes
        // Output       PCyc    : WALTZ-4 pulse cycle

/* ****************************************************************************

     There are 4 cycles associated with WALTZ-8, as shown below
         _ _     _ _ _    _ _     _  _   _ _ _
   U = K K K K = 24231 * 24231 * 24231 * 24231     K = 180  360 180  270  90
                                                          -x   x   -x   x   -x

     Each K represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A WALTZ8 pulse train will employ K = WALTZ-K.
     The WALTZ-8 cycle is simply the phi phi+180 phi+180 phi sequence.

******************************************************************************/


MSVCDLL PulCycle CycWALTZ8(const spin_system& sys, int even=0) const;

        // Input        void    : None
        //              phi     : Phase angle (degrees)
	//		even	: Flag if even step sizes
        // Output       PCyc    : WALTZ-8 pulse cycle
 
/******************************************************************************

     There are 4 cycles associated with WALTZ-16, as shown below
               _ _     _ _ _ _ _    _ _ _ _     _ _ _ _    _ _ _ _ _
         U = Q Q Q Q = 342312423 * 342312423 * 342312423 * 342312423

                  Q = 270  360 180  270 90  180 360  180 270
                         -x   x   -x   x  -x   x   -x   x   -x

     Each Q represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A WALTZ16 pulse train will employ Q = WALTZ-Q

******************************************************************************/


MSVCDLL PulCycle CycWALTZ16(const spin_system& sys, int even=0) const;

        // Input        void    : None
        //              phi     : Phase angle (degrees)
	//		even	: Flag if even step sizes
        // Output       PCyc    : WALTZ-16 pulse cycle


// ____________________________________________________________________________
// H                    WALTZ PULSE SUPERCYCLE FUNCTIONS
// ____________________________________________________________________________
 
// ____________________________________________________________________________
// I                        WALTZ PULSE TRAIN FUNCTIONS
// ____________________________________________________________________________

// ____________________________________________________________________________
// K                      CLASS WALTZ INPUT FUNCTIONS
// ____________________________________________________________________________
 
 
MSVCDLL void read(const std::string &filename, int idx=-1);
 
        // Intput               GP      : WALTZ parameters
        //                      idx     : WALTZ index
        // Output               none    : WALTZ parameters are read in
        //                                from parameters in file filename
        //                                with index idx
 
 
MSVCDLL void read(const ParameterSet& pset, int idx=-1);
 
        // Intput               GP      : WALTZ parameters
        //                      pset    : Parameter set
        //                      idx     : WALTZ index
        // Output               none    : WALTZ parameters are read in
        //                                from parameters in pset
        //                                with index idx
 
//                  Read WALTZ Delay & Pulse Parameters
// For The Pulse We Need Two of Three: {Pulse Angle, RF Strength, Pulse Length}
 
// Note: Reading order is angle --> length --> strength for WALTZ.  If only
//       the angle is specified the length will be set to zero (ideal pulse)
 
 

MSVCDLL void ask_read(int argc, char* argv[], int argn);
 
        // Intput               GP      : WALTZ parameters
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
 

// ____________________________________________________________________________
// L                      CLASS WALTZ I/O FUNCTIONS
// ____________________________________________________________________________

                                                                                
//std::ostream& printBase(std::ostream &ostr) const;

        // Intput               GP      : WALTZ parameters
        //                      ostr    : Output stream
        // Output               none    : WALTZ basic parameters are sent
        //                                to the output stream

                                                               
MSVCDLL std::ostream& print(std::ostream &ostr) const;

        // Intput               GP      : WALTZ parameters
        //                      ostr    : Output stream
        // Output               none    : WALTZ parameters are sent
        //                                to the output stream
 

MSVCDLL friend std::ostream &operator << (std::ostream &ostr, const WALTZ &GP);
 
        // Intput               GP      : WALTZ parameters
        //                      ostr    : Output stream
        // Output               none    : WALTZ parameters are sent
        //                                to the output stream
};

// ____________________________________________________________________________
// AA            ADDITIONAL WALTZ PHASE CYCLE FUNCTIONS
// ____________________________________________________________________________

/******************************************************************************

     A WALTZ-4 pulse cycle is given by
                     _ _    _   _  _ _ _ _
             U = R R R R = 123 123 123 123       R = 90 180  270
                                                       x   -x   x

     Each R represents a pulse waveform and the bar indicates a   
     180 degree phase shift.  A WALTZ4 pulse train will employ R = WALTZ-R 
 
******************************************************************************/

MSVCDLL row_vector CYC_WALTZ4(double phi=0);

        // Input        void    : None
	//		phi	: Phase angle (degrees)
        // Output       PCyc	: WALTZ-4 pulse cycle
 


/******************************************************************************

     A WALTZ-8 pulse cycle is given by
         _ _     _ _ _    _ _     _  _   _ _ _
   U = K K K K = 24231 * 24231 * 24231 * 24231     K = 180  360 180  270  90
                                                          -x   x   -x   x   -x

     Each K represents a pulse waveform and the bar indicates a                
     180 degree phase shift.  A WALTZ8 pulse train will employ K = WALTZ-K
 
******************************************************************************/
 
 
MSVCDLL row_vector CYC_WALTZ8(double phi=0);
 
        // Input             sys   : Spin system
        //                   H     : Static Hamlitonian without the field
        //                   Iso   : Isotope channel the field is on
        //                   gamB1 : The rf-field strength (Hz)
        // Output            U     : Propagator for an WALTZ-8 sequence
        //                           about the axis specified by phi
        // Note                    : This propagator assumes one is in
        //                           the rotating frame of the rf-field
 



// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________




// AA                       WALTZ WAVEFORM FUNCTIONS
//PulWaveform WF_WALTZR(double gamB1, double phi=0);

        // Input        tp      : Single Pi Pulse length (sec)
        // Output       Cvect   : Vector of composite pulse waveform
        // Note                 : Uses constant rf-amplitude and increment
        //                        time, thus is constructed as 4 steps
/*
     There are 3 steps in a WALTZ-R sequence, as shown below.
                                  _                _
                        U = U(90)*U(180)*U(270) = 123
 
     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift.                                     */


//PulWaveform WF_WALTZK(double gamB1, double phi=0);

        // Input        tp      : Single Pi Pulse length (sec)
        // Output       Cvect   : Vector of composite pulse waveform
        // Note                 : Uses constant rf-amplitude and increment
        //                        time, thus is constructed as 4 steps
/*
     There are 5 steps in a WALTZ-K sequence, as shown below.
                 _             _             _       _ _ _
             U = U(180)*U(360)*U(180)*U(270)*U(90) = 24231

     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift.                                     */
 

//PulWaveform WF_WALTZQ(double gamB1, double phi=0);

        // Input        tp      : Single Pi Pulse length (sec)
        // Output       Cvect   : Vector of composite pulse waveform
        // Note                 : Uses constant rf-amplitude and increment
        //                        time, thus is constructed as 4 steps
/*
     There are 9 steps in a WALTZ-Q sequence, as shown below.
    _               _             _            _             _        _ _ _ _ _
U = U(270)*U(360,p)*U(180)*U(270)*U(90)*U(180)*U(360)*U(180)*U(270) = 342312423

     Each U represents a pulse of the specified angle and the bar
     indicates a 180 degree phase shift.                                     */

// ____________________________________________________________________________
// B                     WALTZ COMPOSITE PULSE FUNCTIONS
// ____________________________________________________________________________

/*

PulComposite PCmpWALTZR(const spin_system& sys, const std::string& IsoC,
                                                   double gamB1, double phi=0);
 
        // Input        sys     : Active spin system
        //              IsoC    : Pulse channel
        //              gamB1   : RF field strength (Hz)
        //              phi     : RF field phase (degrees)
        // None         PCom    : A composite pulse is returned set
        //                        to WALTZ-R for the system sys
        //                        of strength gamB1 and phase phi

PulComposite PCmpWALTZK(const spin_system& sys, const std::string& IsoC,
                                                   double gamB1, double phi=0);

        // Input        sys     : Active spin system
        //              IsoC    : Pulse channel
        //              gamB1   : RF field strength (Hz)
        //              phi     : RF field phase (degrees)
        // None         PCom    : A composite pulse is returned set
        //                        to WALTZ-K for the system sys
        //                        of strength gamB1 and phase phi

                                                             
PulComposite PCmpWALTZQ(const spin_system& sys, const std::string& IsoC,
                                                   double gamB1, double phi=0);
*/

        // Input        sys     : Active spin system                            
        //              IsoC    : Pulse channel
        //              gamB1   : RF field strength (Hz)
        //              phi     : RF field phase (degrees)
        // None         PCom    : A composite pulse is returned set
        //                        to WALTZ-Q for the system sys
        //                        of strength gamB1 and phase phi


// ____________________________________________________________________________
// C                       WALTZ PULSE CYCLE FUNCTIONS
// ____________________________________________________________________________

// ____________________________________________________________________________
// D                   WALTZ PULSE SUPERCYCLE FUNCTIONS
// ____________________________________________________________________________


// ____________________________________________________________________________
// E                      WALTZ PULSE TRAIN FUNCTIONS
// ____________________________________________________________________________


//PulTrain PT_WALTZ4(const spin_system& sys, const std::string& IsoC,
//                                                   double gamB1, double phi=0);
 
        // Input        sys     : Active spin system
        //              IsoC    : Pulse channel
        //              gamB1   : RF field strength (Hz)
        //              phi     : RF field phase (degrees)
        // None         PCom    : A composite pulse is returned set
        //                        to WALTZ-R for the system sys
        //                        of strength gamB1 and phase phi
 
/*   There are 4 cycles associated with WALTZ-16, as shown below
               _ _     _ _ _ _ _    _ _ _ _     _ _ _ _    _ _ _ _ _
         U = Q Q Q Q = 342312423 * 342312423 * 342312423 * 342312423

                  Q = 270  360 180  270 90  180 360  180 270
                         -x   x   -x   x  -x   x   -x   x   -x

     Each Q represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A WALTZ16 pulse train will employ Q = WALTZ-Q  */


/*   There are 4 cycles associated with WALTZ-8, as shown below
         _ _     _ _ _    _ _     _  _   _ _ _
   U = K K K K = 24231 * 24231 * 24231 * 24231     K = 180  360 180  270  90
                                                          -x   x   -x   x   -x

     Each K represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A WALTZ8 pulse train will employ K = WALTZ-K.
     The WALTZ-8 cycle is simply the phi phi+180 phi+180 phi sequence        */  
 

//PulTrain PT_WALTZ8(const spin_system& sys, const std::string& IsoC,
//                                                   double gamB1, double phi=0);

        // Input        sys     : Active spin system
        //              IsoC    : Pulse channel
        //              gamB1   : RF field strength (Hz)
        //              phi     : RF field phase (degrees)
        // None         PT      : A pulse train is returned set
        //                        to WALTZ-8 for the system sys
        //                        of strength gamB1 and phase ph


/*   There are 4 cycles associated with WALTZ-16, as shown below
               _ _     _ _ _ _ _    _ _ _ _     _ _ _ _    _ _ _ _ _
         U = Q Q Q Q = 342312423 * 342312423 * 342312423 * 342312423

                  Q = 270  360 180  270 90  180 360  180 270
                         -x   x   -x   x  -x   x   -x   x   -x
 
     Each Q represents a pulse waveform and the bar indicates a
     180 degree phase shift.  A WALTZ16 pulse train will employ Q = WALTZ-Q  */
 
 
//PulTrain PT_WALTZ16(const spin_system& sys, const std::string& IsoC,
//                                                   double gamB1, double phi=0);
 
        // Input        sys     : Active spin system
        //              IsoC    : Pulse channel
        //              gamB1   : RF field strength (Hz)
        //              phi     : RF field phase (degrees)
        // None         PT      : A pulse train is returned set
        //                        to WALTZ-16 for the system sys
        //                        of strength gamB1 and phase phi


#endif						// PulW_WALTZ.h
