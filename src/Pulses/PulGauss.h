/* PulGauss.h ******************************************-*-c++-*-*
**								**
**	                        G A M M A 			**
**								**
**	Gaussian Shaped Pulses                Interface 	**
**								**
**	Copyright (c) 1995				 	**
**	Dr. Scott A. Smith				 	**
**	1800 E. Paul Dirac Drive				**
**	National High Magnetic Field Laboratory			**
**	Tallahassee FL 32306 USA				**
**								**
**      $Header: $
**								**
*****************************************************************/

/*****************************************************************
**								**
**  Description						 	**
**								**
**  This GAMMA module deals with Gaussian shaped pulses. Herein	**
**  are functions which compute propagators for the facile use	**
**  Gaussian pulses as single steps in GAMMA programs. Included **
**  are functions to analyze shaped pulse details, return the	**
**  shaped waveform as a vector, and directly apply such a	**
**  pulse to a spin system.				        **
**							 	**
*****************************************************************/

#ifndef   GPulGauss_h_			// Is file already included?
#  define GPulGauss_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/GenOp.h>		// Include general operators 
#include <HSLib/SpinSystem.h>		// Include isotropic systems
#include <string>			// Include stdlibc++ strings

struct Gpuldat                          // Gpulse info
  {
  int N;                                // Number of steps
  double Wrf;                           // Field frequency
  std::string Iso;                      // Field selectivity
  double gamB1;                         // Field strength (Hz)
  double tau;                           // Pulse length (sec)
  double fact;                          // Cutoff values (%gamB1)
  double phi;                           // Pulse phase (degrees)
  };

// ______________________________________________________________________
//                       Gaussian Pulse Propagators
// ______________________________________________________________________

// -------------------- Individual Pulse-Delay Steps --------------------

// :FIXME: The function declared below does not match any defined in the .cc
// file, so it will be commented out until further notice.
//
//MSVCDLL gen_op Gaussian(spin_system& sys, gen_op& H, std::string& Iso,
//                                 double td, double theta, double phi=0.0);

	// Input	 sys   : Spin system
	// 		     H     : Static Hamiltonian without the field
	// 		     Iso   : Isotope channel the field is on
	// 		     td    : Delay following pulse application
	// 		     theta : Pulse angle (degrees)
	// 		     phi   : Pulse phase (degrees)
	// Output	 U     : Propagator for an Gaussian sequence
	//			         about the axis specified by phi


MSVCDLL void Gpulse_Hs(gen_op* Hs, gen_op& H0, gen_op& Fxy,
                          int N, double ang, double tp, double fact);

        // Input        Hs      : Array of general operators
        //              H0      : Static Hamiltonian *
        //              Fxy     : RF-Field Hamiltonian *
        //              N       : Number of Gaussian steps
        //              ang     : Pulse rotation angle (degrees)
        //              tp      : Gaussian pulse length (sec)
        //              fact    : Cutoff factor [0.0,1.0]
        // Output       void    : Hamiltonians for a Gaussian pulse
        //                        are put into array Hs

// _________________________________________________________________________
//                         Gaussian Pulse Propagators
// _________________________________________________________________________

// ------------- Gaussian Pulse Propagators Without Relaxation -------------

MSVCDLL void Gpulse_Us(gen_op* Us, gen_op& H0, gen_op& Fxy,
                          int N, double ang, double tp, double fact);
 
        // Input        Us      : Array of general operators
        //              H0      : Static Hamiltonian *
        //              Fxy     : RF-Field Hamiltonian *
        //              N       : Number of Gaussian steps
        //              ang     : Pulse rotation angle (degrees)
        //              tp      : Gaussian pulse length (sec)
        //              fact    : Cutoff factor [0.0,1.0]
        // Output       void    : Propagators for a Gaussian pulse
        //                        are put into array Us


MSVCDLL gen_op Gpulse_U(const spin_system& sys, const Gpuldat& Gdata);

        // Input        sys     : Spin system
        //              Gdata   : Gaussian pulse structure
        // Output       U       : Propagator for a Gaussian pulse


MSVCDLL gen_op Gpulse_U(const spin_sys& sys, gen_op& H0, const Gpuldat& Gdata);

        // Input        sys     : Spin system
        //              H0      : Static Hamiltonian
        //              Gdata   : Gaussian pulse structure
        // Output       U       : Propagator for a Gaussian pulse
	// Note                 : The pulse offset in Gdata is taken
        //                        as relative to the rotating frame
        //                        of the input Hamiltonian H0


MSVCDLL gen_op Gpulse_U(gen_op& H0rot, gen_op& Fxy, const Gpuldat& Gdata);

        // Input        H0rot	: Static Hamiltonian
        //              Fxy     : RF-Field Hamiltonian
        //              Gdata   : Gaussian pulse structure
        // Output       U       : Propagators for a Gaussian pulse
        // Note                 : HOrot **MUST** be in the rotating frame
        //                        of the applied pulse
        // Note                 : Fxy **MUST** contain the phase of
        //                        the applied pulse
        // Note                 : In line with the previous two notes,
        //                        Gdata's phase and offset are not used
        //                        herein (they're assumed in H0 & Fxy)

 
MSVCDLL gen_op Gpulse_U(gen_op& H0rot, gen_op& Fxy,
                          int N, double ang, double tp, double fact);

        // Input        H0rot   : Static Hamiltonian
        //              Fxy     : RF-Field Hamiltonian
        //              N       : Number of Gaussian steps
        //              ang     : Pulse rotation angle (degrees)
        //              tp      : Gaussian pulse length (sec)
        //              fact    : Cutoff factor [0.0,1.0]
        // Output       U       : Propagators for a Gaussian pulse
        // Note                 : HOrot is the static Hamiltonian
        //                        in the rotating frame of the pulse.
        // Note                 : Fxy implicitly contains rf-pulse phase
 
 
//              U         = U       * U     * U
//               Gaussian    t>tp/2    tp/2    t<tp/2
 
 
// This is a replica of the previous function that is used for debugging
// It contains additional print statements which indicate how the propagator
// is built step by step.  It isn't documented and may one day disappear.
      
MSVCDLL gen_op Gpulse_UX(gen_op& H0rot, gen_op& Fxy,
                          int N, double ang, double tp, double fact);
 

// _________________________________________________________________________
//                    Gaussian Pulse Auxiliary Functions
// _________________________________________________________________________

// --------------- Gaussian Pulse On Resonance Rotation Angle --------------

MSVCDLL double Gangle(double gamB1, double tau, int N, double fact=0.025);

        // Input        gamB1   : The rf-field strength (Hz)
        //              tau     : Gaussian pulse length (sec)
        //              N       : Number of Gaussian steps
	//		fact    : Cutoff factor [0.0,1.0]
	//			  defaults to 2.5 %
        // Output       angle   : Gaussian pulse rotation angle on
	//			  on resonance

// ------ Gaussian Pulse Strength To Obtain On Resonance Rotation Angle -------

MSVCDLL double GgamB1(double angle, double tau, int N, double fact=0.025);

        // Input        angle   : The pulse rotation angle (degrees)
        //              tau     : Gaussian pulse length (sec)
        //              N       : Number of Gaussian steps
        //              fact    : Cutoff factor
        // Output       angle   : Gaussian pulse rotation angle on
        //                        on resonance

// ------- Gaussian Pulse Length To Obtain On Resonance Rotation Angle --------

MSVCDLL double Gtime(double angle, double gamB1, int N, double fact=0.025);

        // Input        angle   : The pulse rotation angle (degrees)
        //              gamB1   : Gaussian pulse rf-field strength (Hz)
        //              N       : Number of Gaussian steps
        //              fact    : Cutoff factor [0.0,1.0]
        //                        defaults to 2.5 %
        // Output       time    : Gaussian pulse length to obtain
        //                        specified rotation angle on resonance
        //                        for indicated pulse strength & cutoff


// ----------------- Gaussian Pulse RF-Field Intensity Vector ------------------

MSVCDLL row_vector GNvect(int N, double fact);

        // Input        N       : Number of Gaussian steps
        //              fact    : Cutoff factor
        // Output       Gvect   : Vector of Gaussian points
        //                        with a maximum of 1


MSVCDLL row_vector Gvect(double gamB1, int N, double fact=0.025);

        // Input        gamB1   : The rf-field strength (Hz)
        //              N       : Number of Gaussian steps
	//		fact    : Cutoff factor
        // Output       Gvect   : Vector of Gaussian pulse
	//			  rf-field intensities


// ------------------ Gaussian Pulse RF-Field Integral Vector ------------------


MSVCDLL row_vector GIntvec(double gamB1, double tau, int Npts, double fact=0.05);

        // Input        gamB1   : The rf-field strength (Hz)
        //              tau     : Gaussian pulse length (sec)
        //              Npts    : Number of Gaussian steps
        //              fact    : Cutoff factor
        // Output       Gvect   : Vector of Gaussian integral


// _________________________________________________________________________
//                   Gaussian Pulse Plotting Functions
// _________________________________________________________________________


MSVCDLL row_vector Ghistogram(double gamB1, double tau, int N, double fact=0.05);

        // Input        gamB1   : The rf-field strength (Hz)
        //              tau     : Gaussian pulse length (sec)
        //              N       : Number of Gaussian steps
	//		fact    : Cutoff factor
        // Output       Gshape  : Vector of point for plot of Gaussian
	//			  as a histogram.


// _________________________________________________________________________
//                 Gaussian Pulse Interactive Setup Functions
// _________________________________________________________________________

 
MSVCDLL void ask_Gpulse(int argc, char* argv[], int& qn, int& N,
                 double& val1, double& val2, double& fact, const int type=0);

        // Input        argc    : Number of argments
        //              argv    : Array of arguments
        //              qn      : Argument index
        //              N       : Number of points
        //              val1    : Either pulse angle or strength
        //              val2    : Either pulse length or strength
        //              fact    : Cutoff percent value
        //              type    : How to characterize the Gaussian
        //                         0 = angle and length (default)
        //                         1 = strength and length
        //                         2 = strength and time 
        // Output       void    : The values N, val1, val2, and fact
	//			  are set either interactively or by
        //                        the values in argv.

 
MSVCDLL Gpuldat read_Gpulse(std::string& filein, spin_system& sys, int idx=-1, int pflag=1);

        // Input        filein  : A filename
        //              sys     : Spin system
        //              idx     : Gaussian pulse index (default none)
        //              pflag   : Print flag (default 1=print)
        // Output       Gdata   : Structure defining the Gaussian pulse


 
MSVCDLL void read_Gpulse(std::string& filein, int& N, double& val1, double& val2,
                         double& fact, int& spin, int type=0, int idx=-1);

        // Input        filein  : A filename
        //              qn      : Argument index
        //              N       : Number of points
        //              val1    : Either pulse angle or strength
        //              val2    : Either pulse length or strength
        //              fact    : Cutoff percent value
	//		spin    : Spin affected by pulse
        //              type    : How to characterize the Gaussian
        //                         0 = angle and length (default)
        //                         1 = strength and length
        //                         2 = strength and time 
        //              idx     : Gaussian pulse index (default none)
        // Output       void    : The values N, val1, val2, and fact
        //                        are set from values in file filein



//MSVCDLL double ask_Gaussian(spin_system& sys, std::string& Iso, gen_op& H, double cutoff=1.e-10);

	// Input		sys	: A spin system
	//			Iso     : String designating an isotope
	//			H	: Isotropic Hamiltonian
	//			cutoff	: An intensity cutoff
	// Output		v       : A transition of H
	// Note				: This routine looks over all single


//MSVCDLL void set_Gaussian(double gamB1, double& tmix, double& tpul, double tdel, int& numb, int& type);

	// Input		gamB1: RF-Field strength (Hz)
	// 			tmix : Total Gaussian mixing time
	// 			tpul : Individual pulse length
	// 			tdel : Individual delay length
	// 			numb : Times Gaussian needs repeating for tmix
	// 			type : Flag for Gaussian type
	// Output		     : Void, argument parameters are set

// _________________________________________________________________________
//                    Gaussian Pulse Output Functions
// _________________________________________________________________________

 
MSVCDLL void print_Gpulse(std::ostream& ostr, Gpuldat& Gdata, int indx=-1);

        // Input        ostr    : An output stream
        //              Gdata   : Gaussian pulse structure
        // Output       void    : Gaussian pulse parameters
        //                        are put into ostr


/*************************************************************************
**									**
**              GAMMA's Gaussian Pulse Implementation			**
**									**
**  A gaussian function centered about time t  is given formally by	**
**                                           o				**
**									**
**                           [        2 /        2  ]			**
**                G(t) = exp | -(t-t ) / (2*sigma ) |			**
**                           [      0 /             ]			**
**									**
**  The discrete function is similar except that for a Gaussian pulse	**
**  we would like:							**
**									**
**	1. The Gaussian maximum centered in the middle of the Gaussian	**
**         points, i.e. t  --> (N-1)/2.					**
**                       0						**
**	2. The Gaussian maximum should be gamma*B1 (in Hz) rather than	**
**	   1 as in the above analog function.  So we scale by gamma*B1.	**
**									**
**      3. For each Gaussian point we apply an rf-field at that amp-	**
**	   litude for a duration delt, where the total pulse length is	**
**         t  = N * delt for an N point Gaussian.			**
**          p								**
**      4. We would like to define sigma in terms of a cutoff, i.e.	**
**         some % of the Gaussian maximum. That is to say, we should	**
**         like to set the Gaussian linewidth such that the first and	**
**         last points are at a set percentage of gamma*B1.		**
**									**
**  For a cutoff of fact (for a cutoff of 2% of max. fact=0.02) we	**
**   need to satisfy the following conditions.				**
**									**
**                                 2	     2				**
**                    fact = exp(-N / 8*sigma )				**
**									**
**  or									**
**                  sigma = N / sqrt[-8*ln(fact)]			**
**									**
**  where N is the number of Gaussian steps taken and N/2 is the peak	**
**  maximum.  Setting the peak maximum to be related to an rf-field	**
**  strength, this leaves us with the formula				**
**									**
**                       [           2               2 ]		**
**            G(i) = exp | [2i-(N-1)] ln(fact)/ (N-1)  |		**
**                       [                             ]		**
**                                                                      **
*************************************************************************/

#endif                                                /* __PULGAUSS_H__ */ 
