/* PulGauss.cc *****************************************-*-c++-*-*
**								**
**	                     G A M M A 				**
**								**
**	Gaussian Shaped Pulses              Implementation	**
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

#ifndef _PulGauss_cc_			// Is file already included?
#define _PulGauss_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Pulses/PulGauss.h>		// Include the interface
#include <HSLib/SpinSystem.h>		// Include spin systems
#include <HSLib/HSprop.h>		// Include propagators
#include <HSLib/HSauxil.h>		// Include propagators
#include <HSLib/GenOp.h>		// Include operators
#include <HSLib/HSprop.h>		// Include operators
#include <LSLib/SuperOp.h>		// Include superoperators
#include <HSLib/HSham.h>		// Knowledge of isotropic Hamiltonians
#include <Level2/acquire1D.h>		// Include knowledge of "acquisitions"
#include <HSLib/SpinOpCmp.h>		// Include composite spin ops
#include <HSLib/SpinOpRot.h>		// Include spin rotation operators
#include <stdlib.h>
#include <string>			// Include libstd++ strings
#include <Basics/StringCut.h>		// Include Gdec and Gform functions
#include <iostream>                     // Include input output streams
#include <list>				// Include libstdc++ STL lists

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ lists
using std::ostream;			// Using libstdc++ output streams
using std::cout;			// Using libstdc++ standard output
using std::cin;				// Using libstdc++ standard input

// ______________________________________________________________________
//                         Gaussian Pulse Propagators 
// ______________________________________________________________________

// ---------------------- Real Pulses, No Relaxation --------------------

// This function does not seem to be used, and is not in the .h file
// so will comment it out.
/*
gen_op Gaussian(spin_system& sys, string& Iso, double Wrf,
                         double tp, int N, double theta, double phi=0.0)

  // Input             sys   : Spin system
  //                   Iso   : Isotope channel of pulse
  //		     Wrf   : Freqency of the applied field
  //                   tp    : Gaussian pulse time
  //		     N     : Number of steps
  //                   theta : Gaussian pulse angle (degrees)
  //                   phi   : Gaussian pulse phase (degrees)
  // Output            U     : Propagator for a Gaussian pulse
  //			     of length tp applied at frequency
  //			     Wrf on channel Iso and phase phi.
  //			     The pulse rotates magnetization on
  //			     resonance by angle theta
  // Note                    : Real Pulses, No Relaxation

  // sosi - this isn't working yet
   {
   if(Iso.length()) Wrf=0; // for compiler
   theta = phi; // for compiler
   gen_op H0 = Ho(sys);				// Isotropic Hamiltonian
   gen_op HRF, H;				// RF, Total Hamiltonians

   double gamB1=0, fact=0;
   double tdiv = tp/double(N);			// Single step time
   row_vector Gs = Gvect(gamB1,N,fact);		// Gaussian intensities
   gen_op U, Ustep, Fxy;
   gen_op Hstep;
   for(int i=0; i<N; i++)			// Loop Gaussian steps
     {
     // should build two in reverse order then multiply em!
     Hstep = H0;
     Hstep -= Gs.get(i)*Fxy;			// Total Ham. this step
     Ustep = prop(H, tdiv);			// Prop for this step
     U *= Ustep;
     }
   return U;					// Return Gaussian prop.
   }
*/

// This next function is also not declared in the .h file
// so will be commented out for now.
/*
gen_op Gaussian(spin_system& sys, string& Iso, double Wrf,
                                double tp, double theta, double phi=0.0)

    // Input             sys   : Spin system
    //                   Iso   : Isotope channel of pulse
	//		     Wrf   : Freqency of the applied field
    //                   tp    : Gaussian pulse time
    //                   theta : Gaussian pulse angle (degrees)
    //                   phi   : Gaussian pulse phase (degrees)
    // Output            U     : Propagator for a Gaussian pulse
	//			     of length tp applied at frequency
	//			     Wrf on channel Iso and phase phi.
	//			     The pulse rotates magnetization on
	//			     resonance by angle theta
    // Note                    : Real Pulses, No Relaxation

// sosi - this isn't working yet
{
   if(sys.spins()) Wrf=0; // for compiler
   if(Iso.length()) tp = 0; // for compiler
   theta = phi; // for compiler
   gen_op U;
   return U;					// Return Gaussian prop.
}
*/

// ______________________________________________________________________
//                       Gaussian Pulse Hamiltonians
// ______________________________________________________________________

 void Gpulse_Hs(gen_op* Hs, gen_op& H0, gen_op& Fxy,
                              int N, double ang, double tp, double fact)

        // Input        Hs	: Array of general operators
	//		H0      : Static Hamiltonian *
	//		Fxy	: RF-Field Hamiltonian *
        //              N       : Number of Gaussian steps
	//		ang     : Pulse rotation angle (degrees)
        //              tp      : Gaussian pulse length (sec)
	//		fact    : Cutoff factor [0.0,1.0]
        // Output       void    : Hamiltonians for a Gaussian pulse
	//			  are put into array Hs

   {
   double gamB1 = GgamB1(ang, tp, N, fact);	// Get field strength
   row_vector Gs = Gvect(gamB1, N, fact);	// RF intensities
   Fxy.Op_base(H0);				// Set Fxy basis for easy add
   for(int i=0; i<N; i++)			// Loop over Gaussian steps
     {
     if(N-1-i < i)                              //      Use symmetry to avoid
       Hs[i] = Hs[N-1-i];			//      recalculating same Hs
     else
       {
       Hs[i] = H0;
       Hs[i] -= Gs.get(i)*Fxy;			// Total Ham., this step
       }
     }
   }

// _________________________________________________________________________
//                         Gaussian Pulse Propagators
// _________________________________________________________________________

// ------------- Gaussian Pulse Propagators Without Relaxation -------------

void Gpulse_Us(gen_op* Us, gen_op& H0, gen_op& Fxy,
                              int N, double ang, double tp, double fact)

        // Input        Us	: Array of general operators
	//		H0      : Static Hamiltonian *
	//		Fxy	: RF-Field Hamiltonian *
        //              N       : Number of Gaussian steps
	//		ang     : Pulse rotation angle (degrees)
        //              tp      : Gaussian pulse length (sec)
	//		fact    : Cutoff factor [0.0,1.0]
        // Output       void    : Propagators for a Gaussian pulse
	//			  are put into array Us

   {
   double tdiv = tp/double(N);			// Incremental prop. time
   double gamB1 = GgamB1(ang, tp, N, fact);	// Get field strength
   row_vector Gs = Gvect(gamB1, N, fact);	// RF intensities
   Fxy.Op_base(H0);				// Set Fxy basis for easy add
   gen_op H;					// Scratch operator
   for(int i=0; i<N; i++)			// Loop over Gaussian steps
     {
     if(N-1-i < i)                              //      Use symmetry to avoid
       Us[i] = Us[N-1-i];			//      recalculating same Us
     else
       {
       H = H0;
       H -= Gs.get(i)*Fxy;			// Total Ham., this step
       Us[i] = prop(H, tdiv);			// Prop. this pulse step
       }
     }
   }


 gen_op Gpulse_U(const spin_system& sys, const Gpuldat& Gdata)

        // Input 	sys     : Spin system
	//		Gdata   : Gaussian pulse structure
        // Output       U       : Propagator for a Gaussian pulse

   {
   gen_op H = Ho(sys);				// Isotropic Hamiltonian
   return Gpulse_U(sys, H, Gdata);		// Use function overload
   }

 gen_op Gpulse_U(const spin_sys& sys, gen_op& H0, const Gpuldat& Gdata)

        // Input 	sys     : Spin system
	//		H0      : Static Hamiltonian
	//		Gdata   : Gaussian pulse structure
        // Output       U       : Propagator for a Gaussian pulse
	// Note 		: The pulse offset in Gdata is taken
	//			  as relative to the rotating frame
	//			  of the input Hamiltonian H0

   {
   double tp = Gdata.tau;			// Pulse length (sec)
   int N = Gdata.N;				// Pulse steps
   double gamB1 = Gdata.gamB1;			// Pulse strength (Hz)
   double fact = Gdata.fact;			// Pulse end cutoff (%)
   double W = Gdata.Wrf;			// Pulse offset (Hz)
   string I = Gdata.Iso;			// Pulse selectivity
   double phi = Gdata.phi;			// Pulse phase
   double ang = Gangle(gamB1, tp, N, fact);	// Pulse angle
   gen_op Hrot = H0;				// H in input rotating frame
   if(W) Hrot += W*Fz(sys, I);			// H in pulse rotating frame
   gen_op FXY =  Fxy(sys, phi);			// RF field component (phased)
   return Gpulse_U(Hrot, FXY, N, ang, tp, fact);// Use function overload
   }


gen_op Gpulse_U(gen_op& H0rot, gen_op& Fxy, const Gpuldat& Gdata)

        // Input 	H0rot   : Static Hamiltonian
	//		Fxy	: RF-Field Hamiltonian
	//		Gdata   : Gaussian pulse structure
        // Output       U       : Propagator for a Gaussian pulse
	// Note			: HOrot **MUST** be in the rotating frame
	//			  of the applied pulse
	// Note			: Fxy **MUST** contain the phase of
	//			  the applied pulse
	// Note			: In line with the previous two notes,
	//			  Gdata's phase and offset are not used
	//			  herein (they're assumed in H0 & Fxy)


   {
   double tp = Gdata.tau;			// Pulse length (sec)
   int N = Gdata.N;				// Pulse steps
   double gamB1 = Gdata.gamB1;			// Pulse strength (Hz)
   double fact = Gdata.fact;			// Pulse cutoff (%)
   double ang = Gangle(gamB1, tp, N, fact);	// Pulse angle (degrees)
   return Gpulse_U(H0rot,Fxy,N,ang,tp,fact);	// Use function overload
   }


gen_op Gpulse_U(gen_op& H0rot, gen_op& Fxy,
                          int N, double ang, double tp, double fact)

        // Input 	H0rot   : Static Hamiltonian *
	//		Fxy	: RF-Field Hamiltonian *
        //              N       : Number of Gaussian steps
	//		ang     : Pulse rotation angle (degrees)
        //              tp      : Gaussian pulse length (sec)
	//		fact    : Cutoff factor [0.0,1.0]
        // Output       U       : Propagator for a Gaussian pulse
	// Note                 : HOrot is the static Hamiltonian
	//			  in the rotating frame of the pulse.
	// Note                 : Fxy implicitly contains rf-pulse phase

//		U         = U       * U     * U
//		 Gaussian    t>tp/2    tp/2    t<tp/2

   {
   double dN = double(N);			// Pulse steps
   double tdiv = tp/dN;				// Inc. prop. time
   double gamB1 = GgamB1(ang, tp, N, fact);	// Get field strength
   row_vector Gs = Gvect(gamB1, N, fact);	// RF intensities
   Fxy.Op_base(H0rot);				// Set Fxy basis for easy add
   double gB1 = Gs.getRe(0);			// Field strength first step
   gen_op H = H0rot;
   H += gB1*Fxy;				// Total Ham., 1st step
   gen_op Utlow = prop(H, tdiv);		// Propagator for 1st step
   gen_op Uthigh = Utlow;			// Propagator for Nth step
   gen_op U, Ustep;				// Working propagators
   int i;
   for(i=1; i<N-1-i; i++)			// Loop 1st half of Gaussian
     { 						// steps, symmetry for 2nd half
     gB1 = Gs.getRe(i);				// RF-field strength this step
     H = H0rot;					// Static Hamiltonian
     H += gB1*Fxy;				// Total Hamiltonian, this step
     Ustep = prop(H, tdiv);			// Prop. for this pulse step
     Utlow &= Ustep;				// Utlow = Ustep*Utlow
     Uthigh *= Ustep;				// Uthigh = Uthigh*Ustep
     }
   if(N == 1) U = Utlow;			// For one step, U is 1st step
   else if(i == N-1-i)				// For odd number steps (odd N)
     {						// still need to do middle step
     gB1 = Gs.getRe(i);				// RF-field of middle step
     H = H0rot;					// Static Hamiltonian
     H += gB1*Fxy;				// Total Ham. , middle step
     U = prop(H,tdiv);				// Prop. for middle step
     U *= Utlow;				// Add in after earlier steps
     U &= Uthigh;				// Add later steps after middle
     }
   else U = Uthigh*Utlow;			// Even # steps (even N), this U
   return U;
   }

// This is a replica of the previous function that is used for debugging
// It contains additional print statements which indicate how the propagator
// is built step by step.  It isn't documented and may one day disappear.

 gen_op Gpulse_UX(gen_op& H0rot, gen_op& Fxy,
                          int N, double ang, double tp, double fact)

        // Input 	H0rot   : Static Hamiltonian *
	//		Fxy	: RF-Field Hamiltonian *
        //              N       : Number of Gaussian steps
	//		ang     : Pulse rotation angle (degrees)
        //              tp      : Gaussian pulse length (sec)
	//		fact    : Cutoff factor [0.0,1.0]
        // Output       U       : Propagator for a Gaussian pulse
        // Note                 : HOrot is the static Hamiltonian
        //                        in the rotating frame of the pulse.
        // Note                 : Fxy implicitly contains rf-pulse phase


//		U         = U       * U     * U
//		 Gaussian    t>tp/2    tp/2    t<tp/2

   {
   int debug=1;
   double dN = double(N);
   double tdiv = tp/dN;				// Inc. prop. time
   double gamB1 = GgamB1(ang, tp, N, fact);	// Get field strength
   row_vector Gs = Gvect(gamB1, N, fact);	// RF intensities
   Fxy.Op_base(H0rot);				// Set Fxy basis for easy add
   double thang, tang = 0;
   double gB1 = Gs.getRe(0);			// Field strength first step
   gen_op H = H0rot;
   H += gB1*Fxy;				// Total Ham., 1st step
   if(debug)
     {
     cout << "\n\t\tPulse Length = " << tp;
     cout << "\n\t\tStep Length = " << tdiv;
     cout << "\n\t\tSteps = " << N;
     cout << "\n\t\tStep 1 at " << gB1 << " Hz";
     cout << ", Angle " << gB1*tdiv*360. << " Degrees";
     tang += gB1*tdiv*360.;
     cout << ", Total Angle " << tang << " Degrees";
     }
   gen_op Utlow = prop(H, tdiv);		// Propagator for 1st step
   gen_op Uthigh = Utlow;			// Propagator for Nth step
   gen_op U, Ustep;				// Working propagators
   double tpx = tdiv;				// Used for debugging
   int nstps = 1;				// Used for debugging
   complex gb = gB1/dN;				// Used for debugging
   int i;
   for(i=1; i<N-1-i; i++)			// Loop 1st half of Gaussian
     { 						// steps, symmetry for 2nd half
     gB1 = Gs.getRe(i);				// RF-field strength this step
     H = H0rot;					// Static Hamiltonian
     H += gB1*Fxy;				// Total Hamiltonian, this step
     Ustep = prop(H, tdiv);			// Prop. for this pulse step
     Utlow &= Ustep;				// Utlow = Ustep*Utlow
     Uthigh *= Ustep;				// Uthigh = Uthigh*Ustep
     if(debug)
       {
       tpx += tdiv;				// Used for debugging
       nstps++;					// Used for debugging
       gb += Gs.get(i)/dN;			// Used for debugging
       cout << "\n\t\tStep " << i+1 << " at " << gB1 << " Hz";
       cout << ", Angle " << gB1*tdiv*360. << " Degrees";
       tang += gB1*tdiv*360.;
       cout << ", Total Angle " << tang << " Degrees";
       }
     }
   thang = tang;				// Used for debugging
   if(N == 1) U = Utlow;			// For one step, U is 1st step
   else if(i == N-1-i)				// For odd number steps (odd N)
     {						// need to still do middle step
     gB1 = Gs.getRe(i);
     H = H0rot;
     H += gB1*Fxy;			// Total Ham. , middle step
     U = prop(H,tdiv);
     U *= Utlow;
     U &= Uthigh;
     if(debug)
       {
       cout << "\n\t\tStep " << i+1 << " at " << gB1 << " Hz";
       cout << ", Angle " << gB1*tdiv*360. << " Degrees";
       tang += gB1*tdiv*360.;
       cout << ", Total Angle " << tang << " Degrees (MIDDLE)";
       cout << "\n\t\tEnd Steps Add an Additional Angle of " << thang;
       cout << " Degrees, Total Angle " << thang+tang << " Degrees";
       }
     }
   else
     {
     U = Uthigh*Utlow;			// Even # steps (even N), this U
     if(debug)
       {
       cout << "\n\t\tEnd Steps Add an Additional Angle of " << thang;
       cout << " Degrees, Total Angle " << thang+tang << " Degrees";
       }
     }
  if(debug)
    {
    if(N==1)
      {
      tpx = tdiv; 
      nstps = 1;
      gb =  Gs.get(0)/dN;
      }
    else if(i == N-1-i)
      {
      tpx = 2*tpx+tdiv; 
      nstps = 2*nstps+1;
      gb =  2*gb + Gs.get(i)/dN;
      }
    else
      {
      tpx *=2;
      nstps *= 2;
      gb *= 2.0;
      }
    cout << "\n\n\tSUPPOSED PULSE LENGTH = " << tp;
    cout << "\n\tSUMMED PULSE LENGTH = " << tpx;
    cout << "\n\tREQUIRED GAMMA B1 = " << gamB1;
    cout << "\n\tINTEGRATED GAMMA B1 = " << gb;
    cout << "\n\tSUPPOSED NUMBER OF STEPS = " << N;
    cout << "\n\tACTUAL NUMBER OF STEPS = " << nstps;
    cout << "\n\tSUPPOSED PULSE ANGLE = " << ang;
    cout << "\n\tACTUAL PULSE ANGLE = " << tpx*gb*360.0 << "\n";
    }
   return U;
   }

// _________________________________________________________________________
//                    Gaussian Pulse Auxiliary Functions
// _________________________________________________________________________

// --------------- Gaussian Pulse On Resonance Rotation Angle --------------

 void Gprep(int& N, double& fact, double& den, double& Z)

   {
   if(N<1) N = 1;				// Insure at least 1 point

   if(fact>1.0) fact = 1.0;			// Make sure fact is between
   else if(fact<0.0)				// [0,1].  It can't be zero
     {                                          // if 2 or 1 point requested!
     if(N > 2) fact = 0.0;                      // In those cases just set it
     else fact = 1.0;                           // to 1!
     }   
   else if(fact < 1.e-8) fact = 1.e-8;

   den = double((N-1)*(N-1));			// Calculate the denominator 
   Z = 0.0;					// Exponential factor, set to
   if(N > 1) Z = log(fact)/den;			// zero if N=1;
   return;
   }


 double Gangle(double gamB1, double tau, int N, double fact)

        // Input        gamB1   : The rf-field strength (Hz)
        //              tau     : Gaussian pulse length (sec)
        //              N       : Number of Gaussian steps
	//		fact    : Cutoff factor [0.0,1.0]
	//			  defaults to 2.5 %
        // Output       angle   : Gaussian pulse rotation angle on
	//			  on resonance

   {
   double den, Z;				// Denominator & numerator
   Gprep(N, fact, den, Z);			// Prep. Gaussian values
   row_vector GNs = GNvect(N, fact);		// Normalized Gaussian
   double angle = Re(GNs.sum())*gamB1*tau;	// On resonance pulse angle
   return angle*360.0/double(N);		// Return angle in degrees
   }


// ------ Gaussian Pulse Strength To Obtain On Resonance Rotation Angle -------

double GgamB1(double angle, double tau, int N, double fact)
 
        // Input        angle   : The pulse rotation angle (degrees)
        //              tau     : Gaussian pulse length (sec)
        //              N       : Number of Gaussian steps
	//		fact    : Cutoff factor [0.0,1.0]
	//			  defaults to 2.5 %
        //              fact    : Cutoff factor
        // Output       gamB1   : Gaussian pulse stength to obtain
	//			  specified rotation angle on resonance
	//			  for indicated pulse length and cutoff

   {
   double den, Z;				// Denominator & numerator
   Gprep(N, fact, den, Z);			// Prep. Gaussian values
   row_vector GNs = GNvect(N, fact);		// Normalized Gaussian
   double normang = Re(GNs.sum())*tau/double(N);// Normalized rot. angle
   return (angle/360.0)/normang;		// Return gamB1 in Hz
   }


// ------- Gaussian Pulse Length To Obtain On Resonance Rotation Angle --------

   double Gtime(double angle, double gamB1, int N, double fact)
 
        // Input        angle   : The pulse rotation angle (degrees)
        //              gamB1   : Gaussian pulse rf-field strength (Hz)
        //              N       : Number of Gaussian steps
	//		fact    : Cutoff factor [0.0,1.0]
	//			  defaults to 2.5 %
        // Output       time	: Gaussian pulse length to obtain
	//			  specified rotation angle on resonance
	//			  for indicated pulse strength & cutoff

   {
   double den, Z;				// Denominator & numerator
   Gprep(N, fact, den, Z);			// Prep. Gaussian values
   row_vector GNs = GNvect(N, fact);		// Normalized Gaussian
   double normang = gamB1*Re(GNs.sum());	// Get Gaussian integral (Hz)
   return double(N)*angle/(360.0*normang);	// Scale to seconds
   }


// ----------------- Gaussian Pulse RF-Field Intensity Vector ------------------

   row_vector GNvect(int N, double fact)

        // Input        N       : Number of Gaussian steps
	//		fact    : Cutoff factor
        // Output       Gvect   : Vector of Gaussian points
	//			  (in reals) with a maximum of 1

   {
   double num, den, Z;				// Denominator & numerator
   Gprep(N, fact, den, Z);			// Prep. Gaussian values
   row_vector Gs(N);
   for(int i=0; i<N; i++)			// Loop over Gaussian steps
     {
     if(N-1-i < i)                              //      Use symmetry to avoid
       Gs.put(Gs.get(N-1-i), i);		//      recalculating same pts
     else
       {
       num = double(2*i)-double(N-1);		// Index Gaussian to middle
       Gs.put(exp(Z*num*num), i);		// Gaussian intensity
       }
     }
   return Gs;
   }


   row_vector Gvect(double gamB1, int N, double fact)

        // Input        gamB1   : The rf-field strength (Hz)
        //              N       : Number of Gaussian steps
	//		fact    : Cutoff factor
        // Output       Gvect   : Vector of Gaussian pulse
	//			  rf-field intensities

   {
   double den, Z;				// Denominator & numerator
   Gprep(N, fact, den, Z);			// Prep. Gaussian values
   double Gnorm;				// Normalized Gaussian intensity
   double num;
   row_vector Gs(N);				// Vector of Gaussian points
   for(int i=0; i<N; i++)			// Loop over Gaussian steps
     {
     num = double(2*i)-double(N-1);		// Index so Gaussian mid-pulse
     Gnorm = exp(Z*num*num);			// Normalize Gauss. amplitude
     Gs.put(complex(gamB1*Gnorm), i);		// Gaussian intensity
     }
   return Gs;
   }


// ------------------ Gaussian Pulse RF-Field Integral Vector ------------------


 row_vector GIntvec(double gamB1, double tau, int N, double fact)

        // Input        gamB1   : The rf-field strength (Hz)
        //              tau     : Gaussian pulse length (sec)
        //              N       : Number of Gaussian steps
	//		fact    : Cutoff factor
        // Output       Gvect	: Vector of Gaussian integral
	//			  in degrees!

   {
   double den, Z;				// Denominator & numerator
   Gprep(N, fact, den, Z);			// Prep. Gaussian values
   row_vector Gs = Gvect(gamB1, N, fact);	// Vector of Gaussian points
   row_vector Gvect(N+1);			// Vector for integral points
   Gvect.put(complex0, 0);			// First point of int. is 0
   double tdiv = tau/double(N);			// Single step time
   double time = tdiv;				// Time at point i=1
   double Int = 0.0;				// Integral at point i
   Z = tdiv*360.0;				// Reuse Z for conversion
   for(int i=1; i<N; i++)			// Loop over Gaussian steps
     {
     time += tdiv;				// Time at point i
     Int += Gs.getRe(i)*Z;			// Integral (in degrees)
     Gvect.put(complex(time,Int), i);		// Store integral up to here
     }
   return Gvect;
   }


// _________________________________________________________________________
//                   Gaussian Pulse Plotting Functions
// _________________________________________________________________________


 row_vector Ghistogram(double gamB1, double tau, int N, double fact)

        // Input        gamB1   : The rf-field strength (Hz)
        //              tau     : Gaussian pulse length (sec)
        //              N       : Number of Gaussian steps
	//		fact    : Cutoff factor
        // Output       Gshape  : Vector of point for plot of Gaussian
	//			  as a histogram.

   {
   double den, Z;				// Denominator & numerator
   Gprep(N, fact, den, Z);			// Prep. Gaussian values
   double tdiv = tau/double(N);                 // Incremental time
   int M = N+1;
   row_vector Gshape(2*N+M);			// Vector of Gaussian points
   double lastv, lastt;
   int I = 0;

   row_vector Gs = Gvect(gamB1,N,fact);		// Gaussian intensities

   double time = 0.0;				// Time in pulse
   for(int i=0; i<N; i++)			// Loop over Gaussian steps
     {
     if(i)
       {
       Gshape.put(complex(time,Gs.getRe(i-1)),I++);// For horizontals in hist.
       if(M) Gshape.put(complex(time),I++);	// For verticals in hist.
       }
     else if(M)
       Gshape.put(complex(time,0),I++);		// First vertical in hist.
     Gshape.put(complex(time,Gs.getRe(i)),I++);	// Gaussian intensity
     lastv = Gs.getRe(i);			// Store the previous intensity
     lastt = double(i);				// Store the previous point
     if(i==N-1)					// For the last point
       {
       time += tdiv;
       Gshape.put(complex(time,Gs.getRe(i)),I++);// For last horizontal in hist.
       if(M) Gshape.put(complex(time),I++);	// For last vertical in hist.
       }
     time += tdiv;
     }						// (evolve/acq step goes here)
   return Gshape;
   }


// _________________________________________________________________________
//                 Gaussian Pulse Interactive Setup Functions
// _________________________________________________________________________

 
  void ask_Gpulse(int argc, char* argv[], int& qn, int& N,
                  double& val1, double& val2, double& fact, const int type)

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
        //                        are set either interactively or by
        //                        the values in argv.

  {
  query_parameter(argc, argv, qn++,		// Get number of steps(pts)
       "\n\tNumber of Points in Gaussian? ", N);
  double gamB1, time, angle;
  switch(type)
    {
    case 0:					// Gaussian by angle & time
    default:
      query_parameter(argc, argv, qn++,
      "\n\tGaussian Pulse Angle (degrees)? ", angle);
      query_parameter(argc, argv, qn++,
       "\n\tGaussian Pulse Length (sec)? ", time);
      val1 = angle;
      val2 = time;
      break;
    case 1:					// Gaussian by gamB1 & time
      query_parameter(argc, argv, qn++,	
         "\n\tRF-Field Stength (Hz)? ", gamB1);
      query_parameter(argc, argv, qn++,
       "\n\tGaussian Pulse Length (sec)? ", time);
      val1 = gamB1;
      val2 = time;
      break;
    case 2:					// Gaussian by angle & gamB1
      query_parameter(argc, argv, qn++,
      "\n\tGaussian Pulse Angle (degrees)? ", angle);
      query_parameter(argc, argv, qn++,	
         "\n\tRF-Field Stength (Hz)? ", gamB1);
      val1 = angle;
      val2 = gamB1;
      break;
    }
  query_parameter(argc, argv, qn++,             // Get cutoff percent
   "\n\tPercent Intensity at Ends (0.0, 1.0]? ", fact);
  return;
  } 

 
Gpuldat read_Gpulse(string& filein, spin_system& sys, int idx, int pflag)

        // Input        filein  : A filename
	//		sys     : Spin system
	//		idx	: Gaussian pulse index (default none)
	//		pflag   : Print flag (default 1=print)
        // Output       Gdata   : Structure defining the Gaussian pulse

  {
  Gpuldat Gdata;
  ParameterSet pset;				// A parameter set
  pset.read(filein);				// Read pset in
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  string pname, ssfile, pstate;			// Items in each pset entry
  string SI = string("(") + Gdec(idx)		// Name adjustment if indexed
            + string(")");

//		   Read The Number of Gaussian Steps

  int npt;
  pname = string("Gstps");			// Number of Gaussian steps
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())
    {
    (*item).parse(pname,npt,pstate); 		// Get Gaussian steps
    Gdata.N = npt;				// Set Gaussian steps
    }
  else
    {
    if(pflag)
      cout << "\n\tCant Find " << pname << " in "
           << filein
           << ".  Setting 500 Gaussian Steps";
    Gdata.N = 500;
    }

//                   Read RF Cutoff Intensity as %

  double cutoff;
  pname = string("Gcut");			// Gaussian cutoff
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end()) 				// Get gamB1 parameter
    {
    (*item).parse(pname,cutoff,pstate);
    Gdata.fact = cutoff;
    }
  else
    {
    if(pflag)
      cout << "\n\tCant Find " << pname << " in "
           << filein
           << ".   Setting Gaussian Cutoff to 2%";
    Gdata.fact = 0.02;
    }

//     Read Two of Three: {Pulse Angle, RF Strength, Pulse Length} 

  double gamB1, time, angle;			// Gaussian parameters
  int G=0, T=0, A=0;				// Flags if found
  pname = string("Gang");			// Gaussian pulse angle
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// Get angle parameter
    {
    (*item).parse(pname,angle,pstate);
    A = 1;
    }
  pname = string("Glen");			// Gaussian pulse length
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// Get angle parameter
    {
    (*item).parse(pname,time,pstate);
    T = 1;
    }
  if(A+T != 2)
    {
    pname = string("GgamB1");			// Gaussian pulse angle
    if(idx >= 0) pname += SI;			// Adjust name if indexed			
    item = pset.seek(pname);			// Pix in parameter list
    if(item != pset.end())			// Get gamB1 parameter
      {
      (*item).parse(pname,gamB1,pstate);
      G = 1;
      }
    }
  if(A+T+G == 2)
    {
    if(A && T)
      {
      Gdata.tau = time;				// Set the time
      Gdata.gamB1=GgamB1(angle,time,npt,cutoff);// Set the field strength
      }
    else if(G && T)
      {
      Gdata.tau = time;				// Set the time
      Gdata.gamB1 = gamB1;			// Set the field strength
      }
    else if(A && G)
      {
      Gdata.tau = Gtime(angle,gamB1,npt,cutoff);// Set the time
      Gdata.gamB1 = gamB1;			// Set the field strength
      }
    }
  else if(A+T+G == 0)
    {
    cout << "\n\n\tTo Define a Gaussian Pulse You"
         << " Must Set 2 of the Following 3"
         << " Parameters in " << filein << ":";
    pname = string("Gang");
    if(idx >= 0) pname += SI;
    cout << "\n\n\t\t" << pname << " Gaussian"
         << " Pulse Angle (Degrees)";
    pname = string("Glen");			// Gaussian pulse length
    if(idx >= 0) pname += SI;			// Adjust name if indexed			
    cout << "\n\t\t" << pname << " Gaussian"
         << " Pulse Lenght (Seconds)";
    pname = string("GgamB1");			// Gaussian pulse angle
    if(idx >= 0) pname += SI;			// Adjust name if indexed			
    cout << "\n\t\t" << pname << " Gaussian"
         << " Pulse Strength (Hz)";
    cout << "\n\nSorry, Cannot Continue.\n\n";
    exit(-1);
    }
  else if(A+T+G == 1)
   {
   if(A)
     {
     pname = string("Glen");			// Gaussian pulse length
     if(idx >= 0) pname += SI;			// Adjust name if indexed			
     if(pflag)
       cout << "\n\t\tCant Find " << pname
            << " in " << filein
            << ".   Setting Gaussian Pulse Length"
            << " to 30 msec";
     Gdata.tau = 0.03;				// Set the time
     Gdata.gamB1=GgamB1(angle,0.03,npt,cutoff);	// Set the field strength
     }
   if(T)
     {
     pname = string("Gang");
     if(idx >= 0) pname += SI;
     if(pflag)
       cout << "\n\t\tCant Find " << pname
            << " in " << filein
            << ".   Setting Gaussian Pulse Angle"
            << " to 90 degrees";
     Gdata.tau = time;				// Set the time
     Gdata.gamB1=GgamB1(90.0,time,npt,cutoff);	// Set the field strength
     }
   if(G)
     {
     pname = string("Glen");			// Gaussian pulse length
     if(idx >= 0) pname += SI;			// Adjust name if indexed			
     if(pflag)
       cout << "\n\t\tCant Find " << pname
            << " in " << filein
            << ".   Setting Gaussian Pulse Length"
            << " to 30 msec";
     Gdata.tau = 0.03;				// Set the time
     Gdata.gamB1 = gamB1;			// Set the field strength
     }
   }
 else 
   {
   cout << "\n\tProblems Reading Gaussian Pulse parameters\n";
   exit(-1);
   }

//                    Read Gaussian Selectivity

// If "Gspin" is declared, apply pulse at shift(Gspin).  If "Gspin"
// is not present, need both the offset frequency and isotope
// selectivity.  The first can be input in either Hz or PPM, PPM
// being preferable with the parameters GWPPM or GW.  The second is
// input using the parameter Giso.  Note that the isotope selectivity
// need not be set in a homonuclear spin system.

  double W;
  string Giso;
  int spin;
  pname = string("Gspin");			// Gaussian cutoff
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// Get angle parameter
    {
    (*item).parse(pname,spin,pstate);	// Get spin
    Gdata.Wrf = sys.shift(spin);		// Gaussian pulse on it
    Gdata.Iso = sys.symbol(spin);		// Gaussian selectivity
    }
  else
    {
    pname = string("Giso");			// Now look for selectivity
    if(idx >= 0) pname += SI;			// based on isotope type
    item = pset.seek(pname);			// Pix in parameter list
    if(item != pset.end())					// Get isotope selectivity
      {
      (*item).parse(pname,Giso,pstate);	// Read the selectivity
      Gdata.Iso = Giso;				// Set the Gaussian selectivity
      if(sys.homonuclear() &&
                       sys.symbol(0) != Giso)
        {
        cout << "\n\tYou've Specified A Gaussian "
             << "Pulse On The " << Giso
             << " Isotope Channel.\n\tYour Spin "
             << "System Contains No Spins Of That "
             << "Type!\n\tDon't Expect Much Effect"
             << " From The Pulse...";
        }
      }
    else if(sys.homonuclear())			// For homonuclear systems
      Gdata.Iso = sys.symbol(0);		// the selectivity is clear
    else
      {
      cout << "\n\tCant Find " << pname << " in "
           << filein;
      cout << "\n\tThis Must Be Specified in a "
           << "Heteronuclear System!";
      exit(-1);
      }
    pname = string("GWPPM");			// Gaussian frequency in PPM
    if(idx >= 0) pname += SI;			// Adjust name if indexed			
    item = pset.seek(pname);			// Pix in parameter list
    if(item != pset.end()) 					// Get pulse carrier freq.
      {
      (*item).parse(pname,W,pstate);	// Get the carrier frequency
      Gdata.Wrf = W*sys.Omega(Gdata.Iso);	// Set the Gaussian frequency
      }
    else
      {
      pname = string("GW");			// Gaussian frequency in Hz
      if(idx >= 0) pname += SI;			// Adjust name if indexed			
      item = pset.seek(pname);			// Pix in parameter list
      if(item != pset.end()) 					// Get pulse carrier freq.
        {
        (*item).parse(pname,W,pstate);	// Get the carrier frequency
        Gdata.Wrf = W;				// Set the Gaussian frequency
        }
      else
        {
        if(pflag)
          cout << "\n\tCant Find " << pname << " in "
               << filein
               << ".     Setting Gaussian Frequency to 0";
        Gdata.Wrf = 0.0;
        }
      }
    }

//		   Read The Gaussian Pulse Phase

  double phi;
  pname = string("Gphi");			// Gaussian phase angle
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())
    {
    (*item).parse(pname,phi,pstate); 	// Get Gaussian phase
    Gdata.phi= phi;				// Set Gaussian phase
    }
  else
    {
    if(pflag)
      cout << "\n\tCant Find " << pname << " in "
           << filein
           << ".   Setting Gaussian Pulse Phase to 0";
    Gdata.phi = 0.0;
    }

  return Gdata;
  }


  void read_Gpulse(string& filein, int& N, double& val1, double& val2,
                              double& fact, int& spin, int type, int idx)

        // Input        filein  : A filename
        //              N       : Number of points
        //              val1    : Either pulse angle or strength
        //              val2    : Either pulse length or strength
        //              fact    : Cutoff percent value
        //              spin    : Spin affected by Gaussian
        //              type    : How to characterize the Gaussian
        //                         0 = angle and length (default)
        //                         1 = strength and length
        //                         2 = strength and time 
	//		idx	: Gaussian pulse index (default none)
        // Output       void    : The values N, val1, val2, and fact
        //                        are set from values in file filein
// sosi - this doesn't handle the frequency and selectivity well yet

  {
  ParameterSet pset;				// A parameter set
  pset.read(filein);				// Read pset in
  SinglePar par;
  string pname, ssfile, pstate;			// Items in each pset entry
  string SI = string("(") + Gdec(idx)		// Name adjustment if indexed
            + string(")");

//		   Read The Number of Gaussian Steps

  pname = string("Gstps");			// Number of Gaussian steps
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  std::list<SinglePar>::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())
     (*item).parse(pname,N,pstate); 	// Get gamB1 parameter
  else
    {
    cout << "\n\tCant Find " << pname << " in "
         << filein;
    cout << ".  Setting 500 Gaussian Steps";
    N = 500;
    }

//     Read Two of Three: {Pulse Angle, RF Strength, Pulse Length} 

  double gamB1, time, angle;
  switch(type)
    {
    case 0:					// Gaussian by angle & time
    case 2:
    default:
      pname = string("Gang");			// Gaussian pulse angle
      if(idx >= 0) pname += SI;			// Adjust name if indexed			
      item = pset.seek(pname);			// Pix in parameter list
      if(item != pset.end())					// Get angle parameter
        (*item).parse(pname,angle,pstate);
      else
        {
        cout << "\n\tCant Find " << pname << " in "
             << filein;
        cout << ".   Setting Gaussian Angle to 90 Degrees";
        angle = 90.0;
        }
      val1 = angle;
      break;
    case 1:					// Gaussian by strength & time
      pname = string("GgamB1");			// Gaussian pulse angle
      if(idx >= 0) pname += SI;			// Adjust name if indexed			
      item = pset.seek(pname);			// Pix in parameter list
      if(item != pset.end())					// Get gamB1 parameter
        (*item).parse(pname,gamB1,pstate);
      else
        {
        cout << "\n\tCant Find " << pname << " in "
             << filein;
        cout << ".  Setting Gaussian Field Strength to 50 Hz";
        gamB1 = 50.0;
        }
      val1 = gamB1;
      break;
    }

  switch(type)
    {
    case 0:					// Gaussian by angle & time
    case 1:					// Gaussian by strength & time
    default:
      pname = string("Glen");			// Gaussian pulse length
      if(idx >= 0) pname += SI;			// Adjust name if indexed			
      item = pset.seek(pname);			// Pix in parameter list
      if(item != pset.end())					// Get angle parameter
        (*item).parse(pname,time,pstate);
      else
        {
        cout << "\n\tCant Find " << pname << " in "
             << filein;
        cout << ".  Setting Gaussian Pulse Length to 30 msec";
        time = 0.03;
        }
      val2 = time;
      break;
    case 2:					// Gaussian by angle & strength
      pname = string("GgamB1");			// Gaussian pulse angle
      if(idx >= 0) pname += SI;			// Adjust name if indexed			
      item = pset.seek(pname);			// Pix in parameter list
      if(item != pset.end())					// Get gamB1 parameter
        (*item).parse(pname,gamB1,pstate);
      else
        {
        cout << "\n\tCant Find " << pname << " in "
             << filein;
        cout << ".  Setting Gaussian Field Strength to 50 Hz";
        gamB1 = 50.0;
        }
      val2 = gamB1;
      break;
    }

//                   Read RF Cutoff Intensity as %

  pname = string("Gcut");			// Gaussian cutoff
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end()) 					// Get gamB1 parameter
         (*item).parse(pname,fact,pstate);
  else
    {
    cout << "\n\tCant Find " << pname << " in "
         << filein;
    cout << ".  Setting Gaussian Cutoff to 2%";
    fact = 0.02;
    }

//                    Read Gaussian Selectivity

// If "Gspin" is declared, apply pulse at shift(Gspin).  If "Gspin"
// is not present, need both GW(frequency) and Giso(selectivity)  

double Wrf;
string Giso;
int OK = 0;
  pname = string("Gspin");			// Gaussian cutoff
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  item = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end()) 					// Get gamB1 parameter
         (*item).parse(pname,spin,pstate);
  else
    {
    pname = string("GW");
    if(idx >= 0) pname += SI;			// Adjust name if indexed			
    item = pset.seek(pname);			// Pix in parameter list
    if(item != pset.end()) 					// Get pulse carrier freq.
      {
      (*item).parse(pname,Wrf,pstate);	// Get the carrier frequency
      pname = string("Giso");			// Now look for selectivity
      if(idx >= 0) pname += SI;			// based on isotope type
      item = pset.seek(pname);			// Pix in parameter list
      if(item != pset.end())					// Get isotope selectivity
        {
        (*item).parse(pname,Giso,pstate);
        OK = 1;
        }
      }
    else if(!OK)
      {
      cout << "\n\tCant Find " << pname << " in "
           << filein;
      cout << ".  Setting Spin Index to 0";
      spin = 0;
      }
    }
  }

// sosi - this function makes us depende upon the Level2 module

double ask_Gaussian(spin_system& sys, string& Iso, gen_op& H, double cutoff)

	// Input		sys	: A spin system
	//			Iso     : String designating an isotope
	//			H	: Isotropic Hamiltonian (Hz)
	//			cutoff	: An intensity cutoff
	// Output		v       : A transition of H (Hz)
	// Note				: This routine looks over all single

  {
  gen_op detect = Fm(sys, Iso);		// Set detection operator to F-
  gen_op sigma = Fx(sys);		// Set system in transverse plane
  super_op L = complexi*Hsuper(H);	// L = -i*[Ho, ] (rad/sec)
  acquire1D ACQ(detect, L);		// Prepare for "Dirac" acquisitions
  TTable1D trans = ACQ.table(sigma);	// Determine "Dirac" Acquisition
  int ntr = trans.size();		// Get the number of transitions
  int itr = 0;
  double v;
  cout << "\n\n\tThere are " << ntr << " " << Iso << " Transitions\n";
  for(int i=0; i<ntr; i++)
    cout << "\n\t\t" << i << ". " << trans.Fr(i)/PIx2;
  cout << "\n\n\tChoose [0-" << ntr-1 << "] To Set Gaussian Repetition,"
       << " Any Other Integer To Specify A Different Value: ";
  cin >> itr;
  if(itr >= ntr || itr < 0)
    {
    cout << "\n\t\t" << "Please Input A Value: ";
    cin >> v;
    }
  else v = trans.Fr(itr)/PIx2;
  cout << "\n";
  return v;
  }


  void set_Gaussian(double gamB1, double& tmix, double& tpul, double tdel, int& numb, int& type)

	// Input		gamB1: RF-Field strength (Hz)
	// 			tmix : Total Gaussian mixing time
	// 			tpul : Individual pulse length
	// 			tdel : Individual delay length
	// 			numb : Times Gaussian needs repeating for tmix
	// 			type : Flag for Gaussian type
	//				0 = no pulse, delay with no relaxation
	//				1 = ideal pulse, delay with no relaxation
	//				2 = real pulse, delay with no relaxation
	//				3 = relaxation pulse, delay with no relaxation
	//				11 = ideal pulse, delay with relaxation
	//				12 = real pulse, delay with relaxation
	//				13 = relaxation pulse, delay with relaxation
	// Output		     : Void, argument parameters are set
	// Note			     : gamB1 should be set to zero on input if
	// 			       no external value is present

    {
    cout << "\n\n\t\tGaussian MIXING SEQUENCE SETUP\n";

    cout << "\n\tTotal Mixing Time Desired? ";	// Get the total mixing time
    cin >> tmix;

    double tangle;
    cout << "\n\tTotal Pulse Rotation Angle? ";	// Get the total pulse angle
    cin >> tangle;

    cout << "\n\tNumber of Pulse-Delay Steps?";	// Get the number of steps
    cin >> numb;

    double angle = tangle/numb;
    tpul = 0.0;
    int repeat = 1;
    string typ;
    type = 0;
    if(angle)
      {
      cout << "\n\tThe Rotation Angle of Each Pulse is " << angle << " Degrees"; 
      if(gamB1)
        {
        while(repeat)
          {
          cout << "\n\tPulse Type Desired [i=ideal pulses, p=pulse, r=relaxation]? ";
          cin >> typ;
          if(typ == "i")
            {
            cout << "\n\tPulses Set Ideal";
            repeat = 0;
            type = 1;
            }
          if(typ == "p")
            {
            tpul = angle/(360.0*gamB1);
            cout << "\n\tEach Pulse Lasts " << tpul << " Seconds at This Field Strength";
            repeat = 0;
            type = 2;
            }
          if(typ == "r")
            {
            tpul = angle/(360.0*gamB1);
            cout << "\n\tEach Pulse Lasts " << tpul << " Seconds at This Field Strength";
            repeat = 0;
            type = 3;
            }
          }
        }
      else
        {
        cout << "\n\tPulses Set Ideal";
        type = 1;
        }
      }
  
    double ttot = tmix - double(numb)*tpul;
    tdel = ttot/double(numb);
    if(tdel)
      {
      cout << "\n\tEach Delay Step Lasts " << tdel << " Seconds";
      repeat = 1;				// Get the type of Gaussian
      while (repeat)
        {
        cout << "\n\n\tInclude Relaxation During Delays "
             << "[y = yes, relaxation," " n=no, real pulses]? ";
        cin >> typ;
        if(typ == "y")
          {
          repeat = 0;
          type += 10;				// Relaxation Desired
          }
        else if(typ == "n")			// No Relaxation Desired
          {
          repeat = 0;
          }
        }
      }
    return;
    }

// _________________________________________________________________________
//                    Gaussian Pulse Output Functions
// _________________________________________________________________________

 
  void print_Gpulse(ostream& ostr, Gpuldat& Gdata, int indx)

        // Input        ostr    : An output stream
	//		Gdata   : Gaussian pulse structure
        // Output       void    : Gaussian pulse parameters
	//			  are put into ostr

    {
    int il=0;
    string ival, blanks = "         ";
    ostr << "\n\n\t\tGaussian Pulse ";
    if(indx>=0) ostr << indx << " ";
    ostr << "Parameters\n";
    ostr << "\n\tGaussian Pulse Frequency (Hz): "
         << Gform("%11.3f", Gdata.Wrf);
    ostr << "\n\tGaussian Pulse Length (sec):   "
         << Gform("%11.3f", Gdata.tau);
    ostr << "\n\tGaussian Pulse Angle(degrees): "
         << Gform("%11.3f", Gangle(Gdata.gamB1, Gdata.tau, Gdata.N, Gdata.fact));
    ostr << "\n\tGaussian Pulse Phase (degrees):"
         << Gform("%11.3f", Gdata.phi);
    ival = Gdec(Gdata.N);
    il = 6-ival.length();
    ostr << "\n\tNumber of Gaussian Pulse Steps: "
         << string(il, ' ') << Gdata.N;
    il = 6-(Gdata.Iso).length();
    ostr << "\n\tGaussian Pulse Selectivity:     "
         << string(il, ' ') << Gdata.Iso;
    ostr << "\n\tGaussian Pulse Strength (Hz):  "
         << Gform("%11.3f", Gdata.gamB1);
    ostr << "\n\tGaussian Pulse Cutoff (%gamB1):" 
         << Gform("%11.3f", Gdata.fact);
    ostr << "\n";
    }


#endif                                             		// PulGauss.cc

