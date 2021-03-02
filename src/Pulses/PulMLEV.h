/* PulMLEV.h ****************************************************-*-c++-*-
**			         					**
** 	                        G A M M A				**
**									**
**	MLEV Pulse Functions				Interface	**
**									**
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
** This file contains the a variety of functions supporting the use	**
** of MLEV pulses and pulse trains in the GAMMA simulation platform.	**
** The basic MLEV waveforms are composite 180 pulses and these are	**
** cycled to build up longer MLEV based pulse trains.  MLEV sequences	**
** are used in both mixing and broad-band decoupling.  For details see	**
**                                                                      **
**                                                                      **
*************************************************************************/

#ifndef   GMLEV_h_			// Is this file already included?
#  define GMLEV_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the implementation
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
# include <Pulses/Pulse.h> 		// Know base pulse parameters
# include <Pulses/PulWaveform.h>	// Know pulse waveforms
# include <Pulses/PulComposite.h>	// Know composite pulses
# include <Pulses/PulCycle.h>		// Know pulse cycles
# include <Matrix/row_vector.h>
# include <Basics/ParamSet.h>


class MLEV : public Pulse
  {

private:
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                       CLASS MLEV ERROR HANDLING
// ____________________________________________________________________________
 
 
void MLEVerror(int eidx, int noret=0) const;

        // Input                MLEV    : MLEV parameters(this)
        //                      eidx    : Error flag
        //                      noret   : Return flag
        // Output               none    : Error Message Output

                                                               
void volatile MLEVfatality(int error) const;

        // Input                MLEV    : MLEV parameters(this)
        //                      eidx    : Error flag
        // Output               none    : Stops execution & error Message
 

void MLEVerror(int eidx, const std::string& pname, int noret=0) const;
 
        // Input                MLEV    : MLEV parameters(this)
        //                      eidx    : Error index
        //                      pname   : String included in message
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
 
 
void volatile MLEVfatality(int eidx, const std::string& pname, int noret=0) const;

        // Input                MLEV    : MLEV parameters                        
        //                      eidx    : Error flag
        //                      pname   : Part of error message
        //                      noret   : Flag for return (0=return)
        // Output               none    : Stops execution & error Message



// ____________________________________________________________________________
// ii                  CLASS MLEV PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

                                                                                
void SetPhase(const ParameterSet& pset, int idx=-1);

        // Intput               DT      : MLEV parameters
        //                      pset    : Parameter set
        //                      idx     : MLEV index
        // Output               none    : MLEV pulse phase read in
        //                                from parameter in pset
        //                                with index idx

                                                         
void SetChannel(const ParameterSet& pset, int idx=-1);

        // Intput               DT      : MLEV parameters
        //                      pset    : Parameter set
        //                      idx     : MLEV index
        // Output               none    : MLEV pulse channel read in
        //                                from parameter in pset
        //                                with index idx
 

int SetGamB1(const ParameterSet& pset, int idx=-1);
 
        // Intput               DT      : MLEV parameters
        //                      pset    : Parameter set
        //                      idx     : MLEV index
        // Output               TF      : MLEV pulse strength read in
        //                                from parameter in pset
        //                                with index idx
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 public:                                                                        

// ____________________________________________________________________________
// A                   CLASS MLEV CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

                                                                                
MSVCDLC MLEV();

        // Input        none	:
        // Output       MLEV	: MLEV parameters(this)
        ///F_list       MLEV	- Constructor

                                              
MSVCDLC MLEV(double gB1, const std::string& ch, double ph=0, double off=0);
         
        // Input        gB1     : RF-field strength (Hz)
        //              ch      : RF-channel
        //              ph      : RF-phase (degrees)
        //              off     : RF-offset (Hz)
        // Output       MLEV	: MLEV parameters(this)

                                                         
MSVCDLC MLEV(const MLEV& PT1);

        // Intput       MLEV1   : MLEV parameters
        // Output       MLEV    : MLEV parameters(this) from MLEV1

 
// ------------------------------ Destruction ---------------------------------


MSVCDLC ~MLEV();
 
        // Intput       MLEV1   : MLEV parameters (this)
        // Output       none    : MLEV is destructed

                                                      
// ------------------------------- Assignment ---------------------------------
 

MSVCDLL MLEV& operator = (const MLEV& MLEV1);

        // Intput       MLEV1   : MLEV parameters
        // Output       MLEV    : MLEV parameters(this) from MLEV1
        ///F_list       =       - Assignment

                                             
// ____________________________________________________________________________
// B                     CLASS MLEV ACCESS FUNCTIONS
// ____________________________________________________________________________
 
/*
std::string channel()  const;                 INHERITED
double strength() const;                 INHERITED
double phase()    const;                 INHERITED
double offset()   const;                 INHERITED                     */  

        // Intput       MLEV1	: MLEV parameters                              
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
 
MSVCDLL PulWaveform WF( ) const;
MSVCDLL PulWaveform WF_C180( ) const;

        // Input        MP      : MLEV parameters
        // Output       PWF     : Pulse waveform for MLEV composite 180
        // Note                 : Constant rf-amplitude each step
        // Note                 : Ideal pulses are allowed


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
 
MSVCDLL PulComposite PCmp(const spin_system& sys) const;
MSVCDLL PulComposite PCmp_C180(const spin_system& sys) const;

        // Input        MP      : MLEV parameters           
        //              sys     : A spin system
        // Output       PCmp	: Composite pulse for MLEV composite 180
        //                        applicable to sys

 
// ____________________________________________________________________________
// G                      MLEV PULSE CYCLE FUNCTIONS
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

 
MSVCDLL PulCycle CycMLEV4(const spin_system& sys) const;
 
        // Input        MP      : MLEV parameters           
        //              sys     : A spin system
        // Output       PCyc    : MLEV-4 pulse cycle for sys

 
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
 
 
MSVCDLL PulCycle CycMLEV8(const spin_system& sys) const;
 
        // Input        MP      : MLEV parameters           
        //              sys     : A spin system
        // Output       PCyc    : MLEV-8 pulse cycle for sys

 
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
 
MSVCDLL PulCycle CycMLEV16(const spin_system& sys) const;
 
        // Input        MP      : MLEV parameters           
        //              sys     : A spin system
        // Output       PCyc    : MLEV-16 pulse cycle for sys


// ____________________________________________________________________________
// H                    MLEV PULSE SUPERCYCLE FUNCTIONS
// ____________________________________________________________________________

// ____________________________________________________________________________ 
// I                        MLEV PULSE TRAIN FUNCTIONS
// ____________________________________________________________________________
 
// ____________________________________________________________________________
// K                      CLASS MLEV INPUT FUNCTIONS
// ____________________________________________________________________________

       
MSVCDLL void read(const std::string &filename, int idx=-1);

        // Intput               MP      : MLEV parameters
        //                      idx     : MLEV index
        // Output               none    : MLEV parameters are read in
        //                                from parameters in file filename
        //                                with index idx

                                                         
MSVCDLL void read(const ParameterSet& pset, int idx=-1);

        // Intput               MP      : MLEV parameters
        //                      pset    : Parameter set
        //                      idx     : MLEV index
        // Output               none    : MLEV parameters are read in
        //                                from parameters in pset
        //                                with index idx
 
//                  Read MLEV Delay & Pulse Parameters
// For The Pulse We Need Two of Three: {Pulse Angle, RF Strength, Pulse Length}

// Note: Reading order is angle --> length --> strength for MLEV.  If only     
//       the angle is specified the length will be set to zero (ideal pulse)

 
MSVCDLL void ask_read(int argc, char* argv[], int argn);

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

                                                                                
// ____________________________________________________________________________
// L                      CLASS MLEV I/O FUNCTIONS
// ____________________________________________________________________________
 
         
MSVCDLL std::ostream& printBase(std::ostream &ostr) const;
 
        // Intput               ML      : MLEV parameters
        //                      ostr    : Output stream
        // Output               none    : MLEV basic parameters are sent
        //                                to the output stream
 
         
MSVCDLL std::ostream& print(std::ostream &ostr) const;
 
        // Intput               GP      : MLEV parameters
        //                      ostr    : Output stream
        // Output               none    : MLEV parameters are sent
        //                                to the output stream

         
MSVCDLL friend std::ostream &operator << (std::ostream &ostr, const MLEV &GP);

        // Intput               GP      : MLEV parameters    
        //                      ostr    : Output stream
        // Output               none    : MLEV parameters are sent
        //                                to the output stream
};
 

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
 
MSVCDLL row_vector CYC_MLEV4(double phi=0);
 
        // Input        void    : None
        //              phi     : Phase angle (degrees)
        // Output       PCyc    : MLEV-4 pulse cycle
 

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
 
MSVCDLL row_vector CYC_MLEV8(double phi=0);

        // Input        void    : None
        //              phi     : Phase angle (degrees)
        // Output       PCyc    : MLEV-4 pulse cycle
 


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
 
MSVCDLL row_vector CYC_MLEV16(double phi=0);
 
        // Input        void    : None
        //              phi     : Phase angle (degrees)
        // Output       PCyc    : MLEV-16 pulse cycle





#endif						// PulMLEV.h
