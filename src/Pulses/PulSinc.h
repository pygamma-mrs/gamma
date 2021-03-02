/* PulSinc.h ***************************************************-*-c++-*-*
**									**
**	                        G A M M A 				**
**									**
**	Sinc Shaped Pulses				Interface 	**
**									**
**	Copyright (c) 1995				 		**
**	Dr. Scott A. Smith				 		**
**	1800 E. Paul Dirac Drive					**
**	National High Magnetic Field Laboratory				**
**	Tallahassee FL 32306 USA					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
**  Description							 	**
**									**
** This GAMMA module deals with sync shaped pulses. Herein are		**
** functions which compute propagators for the facile use sync pulses	**
** as single steps in GAMMA programs. Included are functions to analyze	**
** sync pulse details, return the sync waveform as a vector, and 	**
** directly apply such a pulse to a spin system.		        **
**								 	**
*************************************************************************/

#ifndef   GPulSinc_h_			// Is file already included?
#  define GPulSinc_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/GenOp.h>		// Include general operators 
#include <HSLib/SpinSystem.h>		// Include isotropic systems
#include <string>			// Include stdlibc++ strings

struct SincPulDat 			// Sinc Pulse Info
  {
  int N;                                // Number of steps
  double Wrf;                           // Field frequency
  std::string Iso;                      // Field selectivity
  double gamB1;                         // Field strength (Hz)
  double tau;                           // Pulse length (sec)
  int node; 				// Cutoff node
  double phi;                           // Pulse phase (degrees)
  };

// ____________________________________________________________________________
// A                        Sinc Pulse Hamiltonians
// ____________________________________________________________________________
 
/* These functions are meant to generate an array of Hamiltonians for a Sinc
   pulse waveform.  Each step of the Sinc waveform will have a Hamiltonian,
   For all of these function the Sinc waveform must be specified, that means
   that we need two of { pang, gamB1, tp } as well as the end node and number
   of pulse steps.  In addition, the function will require a static Hamiltonian
   (e.g. the isotropic Hamiltonian) and the operator for the spin component
   of the pulse rf-field (e.g. Fx)                                           */
 
MSVCDLL void SincPulseHs(gen_op* Hs, gen_op& H0, gen_op& Fxy,
                                     int N, double ang, double tp, int node=3);
 
        // Input        Hs      : Array of general operators
        //              H0      : Static Hamiltonian *
        //              Fxy     : RF-Field Hamiltonian *
        //              N       : Number of Sinc steps
        //              ang     : Pulse rotation angle (degrees)
        //              tp      : Sinc pulse length (sec)
        //              fact    : Cutoff factor [0.0,1.0]
        // Output       void    : Hamiltonians for a Sinc pulse
        //                        are put into array Hs

// ____________________________________________________________________________
// B                      Sinc Pulse Propagators
// ____________________________________________________________________________
 
/* These functions are meant to generate an array of propagators for a Sinc
   pulse waveform.  Each step of the Sinc waveform will have a propagator Ui.
   For all of these function the Sinc waveform must be specified, that means
   that we need two of { pang, gamB1, tp } as well as the end node and number
   of pulse steps.  In addition, the function will require a static Hamiltonian
   (e.g. the isotropic Hamiltonian) and the operator for the spin component  */

MSVCDLL void SincPulseUs(gen_op* Us, gen_op& H0, gen_op& Fxy,
                                     int N, double ang, double tp, int node=3);

        // Input        Us      : Array of general operators
        //              H0      : Static Hamiltonian *
        //              Fxy     : RF-Field Hamiltonian *
        //              N       : Number of Sinc steps
        //              ang     : Pulse rotation angle (degrees)
        //              tp      : Sinc pulse length (sec)
        //              fact    : Cutoff factor [0.0,1.0]
        // Output       void    : Propagators for a Sinc pulse
        //                        are put into array Us

// ____________________________________________________________________________
// C                        Sinc Pulse Full Propagator
// ____________________________________________________________________________

/*
  Note that these use the symmetry and only calculate step Hamiltonians and
  propagators over half of the waveform.  That makes it faster than using
  the functions which return Hamiltonian or Propagator arrays for all steps
  then summing/multiplying over them.                                        */


MSVCDLL gen_op SincPulseU(gen_op& H0rot, gen_op& Fxy,
                                    int N, double ang, double tp, int node=3);

        // Input        H0rot   : Static Hamiltonian *
        //              Fxy     : RF-Field Hamiltonian *
        //              N       : Number of Sinc steps
        //              ang     : Pulse rotation angle (degrees)
        //              tp      : Sinc pulse length (sec)
        //              node    : Sinc function end node
        // Output       U       : Propagator for a Sinc pulse
        // Note                 : HOrot is the static Hamiltonian
        //                        in the rotating frame of the pulse.
        // Note                 : Fxy implicitly contains rf-pulse phase

/*                      U     = U       * U     * U                      
                         Sinc    t>tp/2    tp/2    t<tp/2                    */ 



MSVCDLL gen_op SincPulseU(gen_op& H0rot, gen_op& Fxy, const SincPulDat& SD);

        // Input        H0rot   : Static Hamiltonian
        //              Fxy     : RF-Field Hamiltonian
        //              SD      : Sinc pulse structure
        // Output       U       : Propagator for a Sinc pulse
        // Note                 : HOrot **MUST** be in the rotating frame
        //                        of the applied pulse
        // Note                 : Fxy **MUST** contain the phase of
        //                        the applied pulse
        // Note                 : In line with the previous two notes,
        //                        SD's phase and offset are not used
        //                        herein (they're assumed in H0 & Fxy)


// ____________________________________________________________________________
// D                      Sinc Pulse Access Functions
// ____________________________________________________________________________
 
/* A Sinc pulse angle is defined as the total rotation angle experienced by
   a single spin on resonance.  That angle is determined by the applied
   RF-field strength (gamB1), the length of the pulse, and (in this discrete
   treatment) the number of pulse steps.  As such, if we know any of the two
   quantities { angle, gamB1, tp } and the number of steps we have defined
   a Sinc pulse. The following functions allow users to access these values.
 
   Function
   ---------         -------------------------------------------------------
   SincAngle         Given {gamB1, tp,    N } returns effective pulse angle
   SincGamB1         Given {pang,  tp,    N } returns Sinc rf-field strength
   SincTime          Given {pang,  gamB1, N } returns Sinc pulse length   */

MSVCDLL double SincAngle(double gamB1, double tau,   int N, int node=3);
MSVCDLL double SincGamB1(double angle, double tau,   int N, int node=3);
MSVCDLL double SincTime(double angle,  double gamB1, int N, int node=3);

// ____________________________________________________________________________
// E                     Sinc Pulse WaveForm Functions
// ____________________________________________________________________________

/* These returns a vector filled with a Sinc function that is centered about
   the vector mid-point and spanning "node" nodes on each side.  It does
   using the discrete sinc function given by

                                sin[PI*node*(2i-N+1)/N-1]
                       S  = X * -------------------------
                        i        [PI*node*(2i-N+1)/N-1]

   For a normalized Sinc function X will be 1, and for a Sinc associated with a
   pulse X will be the applied RF field strength, gamB1.

                Input   gamB1   : The rf-field strength (Hz)
                        N       : Number of Sinc steps
                        node    : Cutoff node [1,2,3,...] (default 3)
                Output  Sincvect: Vector of Sinc pulse RF intensities 

   Function
   ----------        ----------------------------------------------------------
   SincNVect         Returns a vector with NORMALIZED sinc intensities
   SincVect          Returns a vector with sinc maximum at gamB1
   SincIntVec        Returns a vector of sinc integral values

  An additonal aspect worthy of mention is the flag endzero.  These functions
  will end Sinc at a node, i.e. at zero intensity.  However, for sinc pulses
  it doesn't make any sense to have a zero intensity pulse step.  So, if the
  flag endzero is NOT set it will assume points -1 and N are they that reside
  at the nodes.  Then Sinc pulses will not (very often) have no zero steps.   */

MSVCDLL row_vector SincNVect(const SincPulDat& SD,                     int endzero=0);
MSVCDLL row_vector SincNVect(                       int N, int node=3, int endzero=0);
MSVCDLL row_vector SincVect(const SincPulDat& SD,                      int endzero=0);
MSVCDLL row_vector SincVect(double gamB1,           int N, int node=3, int endzero=0);
MSVCDLL row_vector SincIntVec(const SincPulDat& SD,                    int endzero=0);
MSVCDLL row_vector SincIntVec(double gB1,double tp, int N, int node=3, int endzero=0);

// ____________________________________________________________________________
// F                    Sinc Pulse Parameter Set Functions
// ____________________________________________________________________________

/* These functions read in Sinc pulse parameters from an external ASCII file
   (a GAMMA parameter set file). 

     Function                 Purpose                     Parameters Read
  ---------------  -------------------------------  ---------------------------
  SincSteps        Reads in # of Sinc pulse steps          SincSteps
  SincNode         Reads in Sinc pulse end node            SincNode
  SincStrength     Reads in Sinc pulse strength     SincAng, SincLen, SincGamB1
  SincSelectivity  Reads in Sinc pulse selectivity   SincSpin, SincW, SincIso
  SincPhase        Reads in Sinc pulse phase               SincPhi
  ReadSinc         Reads in Sinc pulse                                       */
  

MSVCDLL int SincSteps(const ParameterSet& pset, int idx=-1, int pf=0);
MSVCDLL int SincNode(const ParameterSet& pset,  int idx=-1, int pf=0);
MSVCDLL int SincStrength(const ParameterSet& pset,SincPulDat& SD,int idx=-1,int pf=0);
MSVCDLL int SincSelectivity(const ParameterSet& pset, const spin_system& sys,
                                     SincPulDat& Sdata, int idx=-1, int pf=0);
MSVCDLL double SincPhase(const ParameterSet& pset, int idx, int pf=0);
MSVCDLL SincPulDat ReadSinc(const std::string& filein,
                                const spin_system& sys, int idx=-1, int pf=0);

// ____________________________________________________________________________
// G                    Sinc Pulse Plotting Functions
// ____________________________________________________________________________


MSVCDLL row_vector SincHistogram(double gamB1, double tau, int N, int node=3);

        // Input        gamB1   : The rf-field strength (Hz)
        //              tau     : Sinc pulse length (sec)
        //              N       : Number of Sinc steps
	//		fact    : Cutoff factor
        // Output       Gshape  : Vector of point for plot of Sinc
	//			  as a histogram.


// ____________________________________________________________________________
// H                 Sinc Pulse Interactive Setup Functions
// ____________________________________________________________________________

MSVCDLL void SincPts(int argc,   char* argv[], int& qn, SincPulDat& SD);
MSVCDLL void SincNode(int argc,  char* argv[], int& qn, SincPulDat& SD);
MSVCDLL void SincTime(int argc,  char* argv[], int& qn, SincPulDat& SD);
MSVCDLL void SincGamB1(int argc, char* argv[], int& qn, SincPulDat& SD);
MSVCDLL void SincAngle(int argc, char* argv[], int& qn, double& pang);
MSVCDLL void SincIso(int argc,   char* argv[], int& qn, SincPulDat& SD);
MSVCDLL void SincWrf(int argc,   char* argv[], int& qn, SincPulDat& SD);
MSVCDLL void SincPhi(int argc,   char* argv[], int& qn, SincPulDat& SD);
MSVCDLL SincPulDat SincAsk(int argc, char* argv[], int& qn, int type=0);
 
        // Input        argc    : Number of argments
        //              argv    : Array of arguments
        //              qn      : Argument index
        //              type    : How to characterize the Sinc
        //                         0 = angle and length (default)
        //                         1 = strength and length
        //                         2 = strength and time
        // Output       SincData: The sinc pulse structer values set
        //                        interactively or by values in argv.
 
 
// _________________________________________________________________________
// I                       Sinc Pulse Output Functions
// _________________________________________________________________________


MSVCDLL void SincPrint(std::ostream& ostr, const SincPulDat& Gdata, int indx=-1);

        // Input        ostr    : An output stream
        //              Gdata   : Sinc pulse structure
        // Output       void    : Sinc pulse parameters
        //                        are put into ostr


MSVCDLL void SincZero(SincPulDat& SincData);
MSVCDLL void SincPrep(int& N, int& node, double& den, double& Z);


#endif								// PulSinc.h
