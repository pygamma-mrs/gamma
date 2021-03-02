/* PulDANTE.h **************************************************-*-c++-*-*
**									**
**	                             G A M M A 				**
**								 	**
**	DANTE Pulse Functions				Interface 	**
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
** This GAMMA module facilitates the use of DANTE pulse trains in NMR	**
** simulations.								**
**								 	**
*************************************************************************/

#ifndef   GDANTE_h_			// Is file already included?
#  define GDANTE_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/SpinSystem.h>		// Include knowledge of spin systems
#include <HSLib/PulseI.h>		// Inlcude knowledge of ideal pulses
#include <HSLib/PulseS.h>		// Inlcude knowledge of pulses
#include <HSLib/HSham.h>		// Knowledge of isotropic Hamiltonians
#include <HSLib/GenOp.h>		// Include knowledge of operators
#include <LSLib/LSprop.h>		// Include knowledge of superop props 
//#include <BWRRelax/relaxProp.h>		// Include knowledge of propagators
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets
#include <string>			// Include stdlibc++ strings

class PulWaveform;
class PulComposite;
class PulTrain;

class DANTE
  {
  int N;                                // Number of steps
  std::string Iso;			// Applied pulse channel
  double te;				// Delay evolution time   (sec)
  double gamB1;                         // Applied pulse strength (Hz)
  double tp;				// Applied pulse length   (sec)
  double ang;				// Applied pulse angle	  (deg)
  double phi;                           // Applied pulse phase    (deg)
  double Wrf;                           // Applied pulse offset   (Hz)
  double tt;				// Full DANTE step length (sec)
  double F;				// DANTE synch frequency  (Hz)

private:

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                CLASS PULSE WAVEFORM ERROR HANDLING
// ____________________________________________________________________________


void DANTEerror(int eidx, int noret=0) const;
 
        // Input                DANTE   : DANTE parameters(this)
        //                      eidx    : Error flag
        //                      noret   : Return flag
        // Output               none    : Error Message Output
 
 
void volatile DANTEfatality(int error) const;
 
        // Input                DANTE   : DANTE parameters(this)
        //                      eidx    : Error flag
        // Output               none    : Stops execution & error Message


void DANTEerror(int eidx, const std::string& pname, int noret=0) const;

        // Input                DANTE   : DANTE parameters(this)
        //                      eidx    : Error index
        //                      pname   : String included in message
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message


// ____________________________________________________________________________
// ii               CLASS DANTE PARAMETER SET FUNCTIONS
// ____________________________________________________________________________


void SetSteps(const ParameterSet& pset, int idx=-1);

        // Intput               DT      : DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               none    : DANTE number of steps read in
        //                                from parameter in pset
        //                                with index idx

 
void SetPhase(const ParameterSet& pset, int idx=-1);

        // Intput               DT      : DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               none    : DANTE pulse phase read in
        //                                from parameter in pset
        //                                with index idx


void SetChannel(const ParameterSet& pset, int idx=-1);

        // Intput               DT      : DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               none    : DANTE pulse channel read in
        //                                from parameter in pset
        //                                with index idx
 
 
int SetAngle(const ParameterSet& pset, int idx=-1);

        // Intput               DT      : DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               TF      : DANTE pulse angle read in
        //                                from parameter in pset
        //                                with index idx


int SetPulLen(const ParameterSet& pset, int idx=-1);

        // Intput               DT      : DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               TF      : DANTE pulse length read in
        //                                from parameter in pset
        //                                with index idx
 
 
int SetGamB1(const ParameterSet& pset, int idx=-1);
 
        // Intput               DT      : DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               TF      : DANTE pulse strength read in
        //                                from parameter in pset
        //                                with index idx


int SetEvLen(const ParameterSet& pset, int idx=-1);
 
        // Intput               DT      : DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               TF      : DANTE delay length read in
        //                                from parameter in pset
        //                                with index idx


int SetFreq(const ParameterSet& pset, int idx=-1);

        // Intput               DT      : DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               TF      : DANTE frequency read in
        //                                from parameter in pset
        //                                with index idx  
 

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
 public:
 
// ____________________________________________________________________________
// A                   CLASS DANTE CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 
 
MSVCDLC DANTE();
 
        // Input        none	:
        // Output	DANTE   : DANTE parameters(this)
        ///F_list       DANTE	- Constructor
 

MSVCDLC DANTE(const DANTE& PT1);
 
        // Intput	DANTE1	: DANTE parameters
        // Output	DANTE   : DANTE parameters(this) from DANTE1
 
 
// ------------------------------ Destruction ---------------------------------
 
 
MSVCDLC ~DANTE();
 
        // Intput	DANTE1	: DANTE parameters (this)
        // Output       none    : DANTE is destructed
 
 
// ------------------------------- Assignment ---------------------------------

                                                                                
MSVCDLL DANTE& operator = (const DANTE& DANTE1);
 
        // Intput	DANTE1	: DANTE parameters
        // Output	DANTE   : DANTE parameters(this) from DANTE1
        ///F_list 	=	- Assignment
 

// ____________________________________________________________________________
// B                     CLASS DANTE ACCESS FUNCTIONS
// ____________________________________________________________________________

MSVCDLL int        steps()    const;
MSVCDLL std::string     channel()  const;
MSVCDLL double     dlength()  const;
MSVCDLL double     strength() const;
MSVCDLL double     plength()  const;
MSVCDLL double     angle()    const;
MSVCDLL double     phase()    const;
MSVCDLL double     offset()   const;
MSVCDLL double     length()   const;
 
        // Intput       DANTE1  : DANTE parameters
        // Output       steps   : DANTE steps
        //              channel : DANTE isotope channel  
        //              dlength : DANTE delay length (sec)  
        //              strength: DANTE pulse strength (Hz)
        //              dlength : DANTE pulse length (sec)
        //              angle   : DANTE pulse angle (sec)
        //              phase   : DANTE pulse phase (sec)
        //              offset  : DANTE pulse offset (Hz)


// ____________________________________________________________________________
// C                  CLASS DANTE HILBERT SPACE PROPAGATORS
// ____________________________________________________________________________


/*
HSprop GetU(const spin_system& sys, gen_op& H);

        // Input             sys   : Spin system
        //                   H     : Static Hamlitonian without the field
        //                   D     : DANTE specifications
        // Output            U     : Propagator for a DANTE step
        // Note                    : The propagator must be built in REVERSE
        //                           order for the pulse preceeds the delay
*/

// ____________________________________________________________________________
// D                CLASS DANTE LIOUVILLE SPACE PROPAGATORS
// ____________________________________________________________________________
 
 
// ____________________________________________________________________________
// E                       DANTE WAVEFORM FUNCTIONS
// ____________________________________________________________________________

 
MSVCDLL friend PulWaveform WF_DANTE(const DANTE& D);
MSVCDLL PulWaveform WF( ) const;
 
        // Input        D       : DANTE specifications
        // Output       WF      : A pulse waveform for DANTE pulse-delay
 

// ____________________________________________________________________________
// F                    DANTE COMPOSITE PULSE FUNCTIONS
// ____________________________________________________________________________


MSVCDLL friend PulComposite CP_DANTE(const spin_system& sys, const DANTE& D);
MSVCDLL PulComposite CP(const spin_system& sys) const;

        // Input        sys     : Spin system
        //              D       : DANTE specifications
        // Output       CP      : A DANTE composite pulse


// ____________________________________________________________________________
// G                       DANTE PULSE TRAIN FUNCTIONS
// ____________________________________________________________________________


MSVCDLL friend PulTrain PT_DANTE(const spin_system& sys, const DANTE& D);
 
        // Input        sys     : Spin system
        //              D       : DANTE specifications
        // Output       CP      : A DANTE composite pulse
 


MSVCDLL PulTrain PT(const spin_system& sys) const;

        // Input        sys     : Spin system
        //              D       : DANTE specifications
        // Output       CP      : A DANTE composite pulse

 
// ____________________________________________________________________________
// Y                      CLASS DANTE INPUT FUNCTIONS
// ____________________________________________________________________________

 
MSVCDLL void read(const std::string &filename, int idx=-1);

        // Intput               DT      : DANTE parameters
        //                      idx     : DANTE index
        // Output               none    : DANTE parameters are read in
        //                                from parameters in file filename
        //                                with index idx

    
MSVCDLL void read(const ParameterSet& pset, int idx=-1);

        // Intput               DT      : DANTE parameters
        //                      pset    : Parameter set
        //                      idx     : DANTE index
        // Output               none    : DANTE parameters are read in
        //                                from parameters in pset
        //                                with index idx
 
//                  Read DANTE Delay & Pulse Parameters
// For The Pulse We Need Two of Three: {Pulse Angle, RF Strength, Pulse Length}
 
// Note: Reading order is angle --> length --> strength for DANTE.  If only
//       the angle is specified the length will be set to zero (ideal pulse)
 
 
MSVCDLL void ask_read(int argc, char* argv[], int argn);
 
        // Intput               DT      : DANTE parameters
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the DANTE parameters
        //                                are read
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : DANTE parameters are modifed (filled)


// ____________________________________________________________________________
// Z                      CLASS DANTE I/O FUNCTIONS
// ____________________________________________________________________________


MSVCDLL std::ostream& printBase(std::ostream &ostr) const;
 
        // Intput               DT      : DANTE parameters
        //                      ostr    : Output stream
        // Output               none    : DANTE basic parameters are sent
        //                                to the output stream


MSVCDLL std::ostream& printInfo(std::ostream &ostr) const;

        // Intput               DT      : DANTE parameters
        //                      ostr    : Output stream
        // Output               none    : DANTE additional information sent
        //                                to the output stream
 

                                                                                
MSVCDLL std::ostream& print(std::ostream &ostr, int full=0) const;

        // Intput               DT      : DANTE parameters
        //                      ostr    : Output stream
        //                      full    : Flag for output amount
        // Output               none    : DANTE parameters are sent
        //                                to the output stream
 
    
MSVCDLL friend std::ostream &operator << (std::ostream &ostr, const DANTE &DT);
 
        // Intput               DT      : DANTE parameters
        //                      ostr    : Output stream
        // Output               none    : DANTE parameters are sent
        //                                to the output stream



  };


// ____________________________________________________________________________
// ____________________________________________________________________________
//                       END OF DANTE CLASS DEFINITION
// ____________________________________________________________________________
// ____________________________________________________________________________
 

// ____________________________________________________________________________
// A                       DANTE WAVEFORM FUNCTIONS
// ____________________________________________________________________________


MSVCDLL PulWaveform WF_DANTE(double td, double gamB1, double tpul, double phi=0);

        // Input        td	: Delay following pulse application
        //              gamB1	: RF-Field strength (Hz)
        //              tpul	: Pulse length
        //              phi	: Pulse phase (degrees)
        // Output       WF      : A pulse waveform for DANTE pulse-delay


//PulWaveform WF_DANTE(const DANTE& D);

        // Input        D	: DANTE specifications
        // Output       WF      : A pulse waveform for DANTE pulse-delay


// ____________________________________________________________________________
// B                     DANTE COMPOSITE PULSE FUNCTIONS
// ____________________________________________________________________________


MSVCDLL PulComposite CP_DANTE(const spin_system& sys, const std::string& Iso,
                           double td, double gamB1, double tpul, double phi=0);

	// Input	sys	: Spin system
	// 		Iso	: Isotope channel the field is on
        //		td	: Delay following pulse application
        //              gamB1	: RF-Field strength (Hz)
        //              tpul	: Pulse length
        //              phi	: Pulse phase (degrees)
        // Output       CP	: A DANTE composite pulse
                                                                   

MSVCDLL PulComposite CP_DANTE(const spin_system& sys, const DANTE& D);
 
        // Input        sys     : Spin system
        //              D       : DANTE specifications
        // Output       CP      : A DANTE composite pulse
 
  
// ____________________________________________________________________________
// C                       DANTE PULSE TRAIN FUNCTIONS 
// ____________________________________________________________________________
  

MSVCDLL PulTrain PT_DANTE(const spin_system& sys, const std::string& Iso,
                           double td, double gamB1, double tpul, double phi=0);
 
        // Input        sys     : Spin system
        //              Iso     : Isotope channel the field is on
        //              td      : Delay following pulse application
        //              gamB1   : RF-Field strength (Hz)
        //              tpul    : Pulse length
        //              phi     : Pulse phase (degrees)
        // Output       PT      : A DANTE pulse train
        // Note                 : Since DANTE is applied without any cycle 
        //                        or supercycle, the DANTE composite pulse
        //                        should be functionally equivalent to this


// ____________________________________________________________________________
// D                         DANTE PULSE PROPAGATORS
// ____________________________________________________________________________

// ----------------------- Individual Pulse-Delay Steps -----------------------

MSVCDLL gen_op UDANTE(const spin_system& sys, gen_op& H, const std::string& Iso,
                                     double td, double theta, double phi=0.0);

	// Input	     sys   : Spin system
	// 		     H     : Static Hamlitonian without the field
	// 		     Iso   : Isotope channel the field is on
	// 		     td    : Delay following pulse application
	// 		     theta : Pulse angle (degrees)
	// 		     phi   : Pulse phase (degrees)
	// Output	     U     : Propagator for an DANTE sequence
	//			     about the axis specified by phi


MSVCDLL gen_op UDANTE(const spin_system& sys, gen_op& H, const std::string& Iso,
                      double td, double gamB1, double tpul, double phi=0.0);

	// Input	     sys   : Spin system
	// 		     H     : Static Hamlitonian without the field
	// 		     Iso   : Isotope channel the field is on
	// 		     td    : Delay following pulse application
	// 		     gamB1 : RF-Field strength (Hz)
	// 		     tpul  : Pulse length
	// 		     phi   : Pulse phase (degrees)
	// Output	     U     : Propagator for an DANTE sequence
	//			     about the axis specified by phi


MSVCDLL gen_op UDANTE(const spin_system& sys, gen_op& H, const DANTE& D);

        // Input             sys   : Spin system
        //                   H     : Static Hamlitonian without the field
        //                   D     : DANTE specifications
        // Output            U     : Propagator for a DANTE step



/*

MSVCDLL super_op GDANTE(spin_system& sys, gen_op& H, std::string& Iso,
                                   super_op& L, gen_op& sigma0,
                                         double td, double theta, double phi);

	// Input	     sys   : Spin system
	// 		     H     : Static Hamlitonian without the field
	// 		     Iso   : Isotope channel the field is on
	//		     L     : Relaxation Liouvillian
	//		     sigma0: Infinite time density matrix
	// 		     td    : Delay following pulse application
	// 		     gamB1 : RF-Field strength (Hz)
	// 		     tpul  : Pulse length
	// 		     phi   : Pulse phase (degrees)
	// Output	     U     : Propagator for an DANTE sequence
	//			     about the axis specified by phi
	// Note			   : Ideal Pulses, Relaxation in Delay


MSVCDLL super_op GDANTE(spin_system& sys, gen_op& H, std::string& Iso,
                                   super_op& L, gen_op& sigma0,
                          double td, double gamB1, double tpul, double phi);

	// Input	     sys   : Spin system
	// 		     H     : Static Hamlitonian without the field
	// 		     Iso   : Isotope channel the field is on
	//		     L     : Relaxation Liouvillian
	//		     sigma0: Infinite time density matrix
	// 		     td    : Delay following pulse application
	// 		     gamB1 : RF-Field strength (Hz)
	// 		     tpul  : Pulse length
	// 		     phi   : Pulse phase (degrees)
	// Output	     U     : Propagator for an DANTE sequence
	//			     about the axis specified by phi


// ---------------------- Total DANTE Pulse Propagators -----------------------


MSVCDLL gen_op UDANTE(spin_system& sys, gen_op& H, std::string& Iso,
                      double v, int n, double theta, double phi=0.0, int rad=0);

	// Input	     sys   : Spin system
	// 		     H     : Static Hamlitonian without the field
	// 		     Iso   : Isotope channel the field is on
	// 		     v     : Dante "frequency" center (Hz)
	// 		     n     : Number of pulse-delay steps
	// 		     theta : Total DANTE Pulse angle (degrees)
	// 		     phi   : Pulse phase (degrees)
	//		     rad   : Flag if frequency input in radians/sec
	// Output	     U     : Propagator for an DANTE sequence
	//			     about the axis specified by phi
	// Note			   : Ideal Pulses, No Relaxation

*/

// ____________________________________________________________________________
//                   PULSE SEQUENCE AUXILIARY FUNCTIONS
// ____________________________________________________________________________


MSVCDLL double ask_DANTE(const spin_system& sys, const std::string& Iso, gen_op& H, double cutoff=1.e-10);

	// Input		sys	: A spin system
	//			Iso     : String designating an isotope
	//			H	: Isotropic Hamiltonian
	//			cutoff	: An intensity cutoff
	// Output		v       : A transition of H
	// Note				: This routine looks over all single


MSVCDLL void set_DANTE(double gamB1, double& tmix, double& tpul, double tdel, int& numb, int& type);

	// Input		gamB1: RF-Field strength (Hz)
	// 			tmix : Total DANTE mixing time
	// 			tpul : Individual pulse length
	// 			tdel : Individual delay length
	// 			numb : Times DANTE needs repeating for tmix
	// 			type : Flag for DANTE type
	// Output		     : Void, argument parameters are set


#endif 						// PulDANTE.h

