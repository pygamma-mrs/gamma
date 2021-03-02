/* PulseShp.cc **************************************************-*-c++-*-
**									**
**	                        G A M M A				**
**									**
**	Shaped Pulses		                 Implementation 	**
**									**
**	Copyright (c) 1990					 	**
**	Scott Smith						 	**
**	Eidgenoessische Technische Hochschule			 	**
**	Labor fuer physikalische Chemie			 		**
**	8092 Zurich / Switzerland				 	**
**							 		**
**      $Header:
**							 		**
*************************************************************************/

/*****************************************************************
**							 	**
**  Description						 	**
**							 	**
**  The folowing functions implement shaped pulses in GAMMA.    **
**  The pulses often take a vector containing the desired	**
**  waveform and approximate the pulse by "integrating" over	**
**  the waveform points.					**
**								**
**  So far, relaxation is not included in these functions.	**
**  This can be done using relaxation propagators though.	**
**							 	**
*****************************************************************/

#ifndef _nmr_shpuls_cc_ 		// Is file already included?
#define _nmr_shpuls_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#endif

#include <HSLib/PulseShp.h>		// Include the header file
#include <HSLib/SpinOp.h>		// Knowledge of spin operators
#include <HSLib/SpinOpRot.h>		// Knowledge of spin rotation ops
#include <HSLib/PulseI.h>		// Knowledge of ideal pulses
#include <HSLib/PulseS.h>		// Knowledfe of soft pulses
#include <HSLib/HSprop.h>		// Knowledge of evolution steps
#include <Basics/Gutils.h>		// Know GAMMA error handling
#include <vector>			// Know stdlibc++ stl vectors


// ____________________________________________________________________________
// i                       Shaped Pulse Error Handling
// ____________________________________________________________________________

 
void PulSherror(int eidx, int noret)
 
        // Input                eidx    : Error index
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error Message Output
 
  {
  std::string hdr("Soft Pulse");
  std::string msg;
  switch(eidx)
    {
    case 1: msg = "Error During Pulse Computation";
            GAMMAerror(hdr, msg, noret); break; //                         (1)
    case 2: msg = "Negative Pulse Lengths Are Disallowed";
            GAMMAerror(hdr, msg, noret); break; //                         (2)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }


void volatile PulShfatality(int eidx)

        // Input                none :
        // Input                eidx    : Error index

  {                                                   
  PulSherror(eidx, 1);                           // Output error message
  if(eidx) PulSerror(0);                        // Now output it fatal
  GAMMAfatal();					// Clean exit from program
  }


// ______________________________________________________________________
// ************** Functions with the Pulse Angle Specified **************
// ______________________________________________________________________

// ______________________ Pulse Along the x_axis ________________________

// _______________ Shaped x-Pulse Density Matrix Evolution ________________

gen_op Shxpuls(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
	 const std::string& iso, double freq, double time, double theta)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			theta : pulse angle (for a spin at resonance)
	// Output		sigma : density matrix following application
	//				of a shaped pulse

  {
  double fact;
  if(time == 0.0)				// If pulse length 0, then 
    return Ixpuls(sys, sigma, iso, theta);	// simply an ideal pulse
  else
    {
    fact = theta/(360*time);			// Scaling for "axis" function
    return Shpul_axis(sys,BLK,sigma,		// Use generic "axis" function
		 H,iso,freq,time,fact,'x');
    }
  }

//_______________________ Shaped x-Pulse Propagator ______________________

gen_op Shxpuls_U(const spin_sys &sys, row_vector& BLK,
           gen_op &H, const std::string& iso, double freq, double time, double theta)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			theta : pulse angle (for a spin at resonance)
	// Output		U     : pulse propagator shaped pulse

  {
  gen_op U;
  double fact;
  if(time == 0.0)
    U = Rx(sys, iso, theta);			// Simply an ideal pulse if time=0
  else
    {
    fact = theta/(360*time);
    U = Shpul_U_axis(sys, BLK, H, iso, freq, time, fact, 'x');
    }
  return U;					// Return a gen_op not a spin_op !
  }

// ______________________ Pulse Along the y_axis ________________________

// _______________ Shaped y-Pulse Density Matrix Evolution ________________

gen_op Shypuls(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
	 const std::string& iso, double freq, double time, double theta)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			theta : pulse angle
	// Output		sigma : density matrix following application
	//				of a shaped pulse.

  {
  double fact;
  if(time == 0.0)
    return Iypuls(sys, sigma, iso, theta);	// Simply an ideal pulse if time=0
  else
    {
    fact = theta/(360*time);
    return Shpul_axis(sys, BLK, sigma, H, iso, freq, time, fact, 'y');
    }
  }

//_______________________ Shaped y-Pulse Propagator ______________________

gen_op Shypuls_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
		double freq, double time, double theta)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	//			H     : static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			theta : pulse angle (for a spin at resonance)
	// Output		U     : pulse propagator for a shaped pulse

  {
  gen_op U;
  double fact;
  if(time == 0.0)
    U = Ry(sys, iso, theta);		// Simply an ideal pulse if time=0
  else
    {
    fact = theta/(360*time);
    U = Shpul_U_axis(sys, BLK, H, iso, freq, time, fact, 'y');
    }
  return U;				// Return a gen_op not a spin_op !
  }

// ______________________ Pulse in the xy_plane ________________________

// _______________ Shaped xy-Pulse Density Matrix Evolution ______________

gen_op Shxypuls(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H, const std::string& iso,
	double freq, double time, double theta, double phi)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			theta : pulse angle
	// 			phi   : pulse phase angle
	// Output		sigma : density matrix following application
	//				of a shaped pulse.

  {
  double fact;
  if(time == 0.0)
    sigma = Ixypuls(sys, sigma, iso, theta, phi);	// Ideal pulse if time=0
  else
    {
    fact = theta/(360*time);
    sigma = Shpul_plane(sys, BLK, sigma, H, iso, freq, time, fact, phi);
    }
  return sigma;
  }


gen_op Shxypuls(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
	 const std::string& iso1, double freq1, const std::string& iso2, double freq2,
		  double time, double theta, double phi)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso1  : isotope type
	// 			freq1 : pulse frequency (offset relative to isotope 1)
	// 			iso2  : isotope type
	// 			freq2 : pulse frequency (offset relative to isotope 2)
	// 			time  : shaped pulse length
	// 			theta : pulse angle
	// 			phi   : pulse phase angle
	// Output		sigma : density matrix following application
	//				of a shaped pulse.

  {
  double fact;
  if(time == 0.0)				// Just do an ideal pulse
    {						// if the time is set zero
    flagvec flags=sys.GetFlags(iso1,1,0);	// Spin flags with iso1 TRUE
    for(int i=0; i<sys.spins(); i++)		// Also set flags for spins
      if(sys.symbol(i) == iso2) flags[i]=1;	// of type iso2 to TRUE
    sigma=Ixypuls(sys,sigma,flags,theta,phi);	// Ideal pulse if time=0
    }
  else
    {
    fact = theta/(360*time);
    sigma = Shpul_plane(sys, BLK, sigma, H, iso1, freq1, iso2, freq2, time, fact, phi);
    }
  return sigma;
  }

//______________________ Shaped xy-Pulse Propagators _____________________

gen_op Shxypuls_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
	double freq, double time, double theta, double phi)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			theta : pulse angle
	// 			phi   : pulse phase angle
	// Output		U     : pulse propagator for a shaped pulse.

  {
  gen_op U;
  double fact;
  if(time == 0.0)
    U = Rxy(sys, iso, theta, phi);	// Rotation operator if time=0
  else
    {
    fact = theta/(360*time);
    U = Shpul_U_plane(sys, BLK, H, iso, freq, time, fact, phi);
    }
  return U;
  }


gen_op Shxypuls_U(const spin_sys &sys, row_vector &BLK, gen_op &H,
	 const std::string &iso1, double freq1, const std::string& iso2, double freq2,
		  double time, double theta, double phi)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	//			H     : acting static Hamiltonian
	// 			iso1  : First isotope type
	// 			freq1 : Pulse frequency (offset relative to isotope 1)
	// 			iso2  : Second isotope type
	// 			freq2 : Pulse frequency (offset relative to isotope 2)
	// 			time  : Shaped pulse length
	// 			theta : pulse angle
	// 			phi   : pulse phase angle
	// Output		U     : propagator for a shaped pulse.

  {
  gen_op U;
  double fact;
  if(time == 0.0)
    {
    flagvec flags=sys.GetFlags(iso1,1,0);	// Spin flags with iso1 TRUE
    for(int i=0; i<sys.spins(); i++)		// Also set flags for spins
      if(sys.symbol(i) == iso2) flags[i]=1;	// of type iso2 to TRUE
    U = Rxy(sys, flags, theta, phi);		// Rotation operator if time=0
    }
  else
    {
    fact = theta/(360*time);
    U = Shpul_U_plane(sys, BLK, H, iso1, freq1, iso2, freq2, time, fact, phi);
    }
  return U;
  }

// ----------------------------------------------------------------------
// ************* Functions with the Field Strength Specified ************
// ----------------------------------------------------------------------

// ______________________ Pulse Along the x_axis ________________________

// _______________ Shaped x-Pulse Density Matrix Evolution ________________

gen_op ShxpulsB(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
	 const std::string& iso, double freq, double time, double gamB1)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// Output		sigma : density matrix following application
	//				of a shaped pulse.

{ return Shpul_axis(sys, BLK, sigma, H, iso, freq, time, gamB1, 'x'); }

//_______________________ Shaped x-Pulse Propagator ______________________


gen_op ShxpulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
		double freq, double time, double gamB1)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// Output		U     : pulse propagator for a shaped pulse

{ return Shpul_U_axis(sys, BLK, H, iso, freq, time, gamB1, 'x'); }

// ______________________ Pulse Along the y_axis ________________________

// _______________ Shaped y-Pulse Density Matrix Evolution ________________

gen_op ShypulsB(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
	 const std::string& iso, double freq, double time, double gamB1)
//	 char *iso, double freq, double time, double gamB1)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// Output		sigma : density matrix following application
	//				of a shaped pulse.

{ return Shpul_axis(sys, BLK, sigma, H, iso, freq, time, gamB1, 'y'); }

//_______________________ Shaped y-Pulse Propagator ______________________

gen_op ShypulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
		double freq, double time, double gamB1)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// Output		U     : pulse propagator for a shaped pulse

{ return Shpul_U_axis(sys, BLK, H, iso, freq, time, gamB1, 'y'); }

// ______________________ Pulse in the xy_plane ________________________

// _______________Shaped xy-Pulse Density Matrix Evolution _______________

gen_op ShxypulsB(const spin_sys &sys, row_vector &BLK,
                 gen_op sigma, gen_op &H, const std::string& iso,
	               double freq, double time, double gamB1, double phi)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// 			phi   : pulse phase angle
	// Output		sigma : density matrix following application
	//				of a shaped pulse.

{ return Shpul_plane(sys, BLK, sigma, H, iso, freq, time, gamB1, phi); }


gen_op ShxypulsB(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
	const std::string& iso1, double freq1, const std::string& iso2, double freq2,
		 double time, double gamB1, double phi)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso1  : isotope type
	// 			freq1 : pulse frequency (offset relative to isotope 1)
	// 			iso2  : isotope type
	// 			freq2 : pulse frequency (offset relative to isotope 2)
	// 			time  : shaped pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// 			phi   : pulse phase angle
	// Output		sigma : density matrix following application
	//				of a shaped pulse.

{ return Shpul_plane(sys, BLK, sigma, H, iso1, freq1, iso2, freq2, time, gamB1, phi); }

//______________________ Shaped xy-Pulse Propagators ______________________

gen_op ShxypulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
	double freq, double time, double gamB1, double phi)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// 			phi   : pulse phase angle
	// Output		U     : propagator for a shaped pulse

{ return Shpul_U_plane(sys, BLK, H, iso, freq, time, gamB1, phi); }


gen_op ShxypulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H,
	const std::string& iso1, double freq1, const std::string& iso2, double freq2,
		double time, double gamB1, double phi)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	//			H     : acting static Hamiltonian
	// 			iso1  : isotope type
	// 			freq1 : pulse frequency (offset relative to isotope 1)
	// 			iso2  : isotope type
	// 			freq2 : pulse frequency (offset relative to isotope 2)
	// 			time  : shaped pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// 			phi   : pulse phase angle
	// Output		U     : propagator for a shaped pulse

{ return Shpul_U_plane(sys, BLK, H, iso1, freq1, iso2, freq2, time, gamB1, phi); }

// ----------------------------------------------------------------------
// ********************* Generic Shaped Pulse Functions *******************
// ----------------------------------------------------------------------

// _____________________ Pulse Along a Specific Axis ____________________

// _______ Shaped Pulse on Coordinate Axis Density Matrix Evolution _______

gen_op Shpul_axis(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
		 const std::string& iso, double freq, double time, double fact, char axis)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// Output		sigma : density matrix following application
	//				of a shaped pulse.

{
  gen_op U;
  U = Shpul_U_axis(sys, BLK, H, iso, freq, time, fact, axis);
  sigma = evolve(sigma, U);
  return sigma;
}

//_______________ Shaped Pulse on Coordinate Axis Propagator ______________

gen_op Shpul_U_axis(const spin_sys &sys, row_vector& BLK, gen_op& H,
					const std::string& iso, double freq, double time, double fact, char axis)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	//			H     : active static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	//			axis  : axis along which the pulse is applied
	// Output		U     : pulse propagator for the application
	//				of a shaped pulse.

{
  gen_op U, Ui;
  double timei, facti, lastfacti; 
  int steps;					// Number of steps in shape
  if(time < 0.0)				// Evolve for a negative time?
    {
    PulSherror(1, 1);                           // Error during pulse calculation
    PulShfatality(2); 				// Negative pulse length
    }
  else if(time == 0.0)				// Evolve for zero time?
    U = I_gen_op(H.get_basis());		// Simply the unity operator
  else
    {						// Adjust input Hamiltonian to frame
    steps = BLK.size();				// Number of steps1
    timei = time/double(steps-1);		// Compute time increment per step
    fact = fact/Re((BLK.sum() - BLK(0)));	// Normalize the factor to BLK
    U = I_gen_op(H.get_basis());		// Unity operator initially
    lastfacti = 0;
    for(int i=1; i<steps; i++)
      {
      facti = Re(BLK(i))*fact*(steps-1);	// Scaling factor
      if(facti != lastfacti)					// Pulse proagator this waveform section
        Ui = Spul_U_axis(sys, H, iso, freq, timei, facti, axis);
      U = Ui * U;
      lastfacti = facti;
      }
    }
  return U;
  }

// ____________________ Pulse Applied in the xy-Plane ___________________

// _______ Shaped Pulse in the xy-Plane Density Matrix Evolution ________

gen_op Shpul_plane(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
	 	const std::string& iso, double freq, double time, double fact, double phi)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// 			phi   : pulse phase angle
	// Output		sigma : density matrix following application
	//				of a shaped pulse.

  {
  gen_op U;
  U = Shpul_U_plane(sys, BLK, H, iso, freq, time, fact, phi);
  sigma = evolve(sigma, U);
  return sigma;
  }


gen_op Shpul_plane(const spin_sys &sys, row_vector &BLK, gen_op sigma,
        gen_op &H, const std::string& iso1, double freq1, const std::string& iso2,
                               double freq2, double time, double fact, double phi)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso1  : isotope type
	// 			freq1 : pulse frequency (offset from isotope 1)
	// 			iso2  : isotope type
	// 			freq2 : pulse frequency (offset from isotope 2)
	// 			time  : shaped pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// 			phi   : pulse phase angle
	// Output		sigma : density matrix following application
	//				of a shaped pulse.

{
  gen_op U;
  U = Shpul_U_plane(sys, BLK, H, iso1, freq1, iso2, freq2, time, fact, phi);
  sigma = evolve(sigma, U);
  return sigma;
}

//__________________ Shaped Pulse in xy-Plane Propagators _________________

gen_op Shpul_U_plane(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
			 double freq, double time, double fact, double phi)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : shaped pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// 			phi   : pulse phase angle
	// Output		U     : propagator for a shaped pulse

{
  gen_op U, Ui;
  double timei, facti, lastfacti; 
  int steps;
  if(time < 0.0)
    {
    PulSherror(1, 1);                           // Error during pulse calculation
    PulShfatality(2); 				// Negative pulse length
    }
  else if(time == 0.0)
    U = I_gen_op(H.get_basis());		// Unity Operator if time = 0
  else
    {
    steps = BLK.size();				// Number of steps1
    timei = time/double(steps-1);		// Compute time increment per step
    fact = fact/Re((BLK.sum() - BLK(0)));	// Normalize the factor to BLK
    U = I_gen_op(H.get_basis());		// Unity Operator initially
    lastfacti = 0;
    for(int i=1; i<steps; i++)
      {

// Note: should test here for 0 or very small facti and skip them. 
//       do this by testing the factor relative to maximum where
//       the maximum would be computed with lines such as
//	 double maxfacti;
//	 maxfacti = Re(BLK(i))*fact*(steps-1);
//	 if(facti < maxfacti*1.e-6) .... 

      facti = Re(BLK(i))*fact*(steps-1);	// Scaling factor
      if(facti != lastfacti)
        Ui = Spul_U_plane(sys, H, iso,		// Pulse proagator this waveform section 
		 freq, timei, facti, phi);
      U = Ui * U;
      lastfacti = facti;
      }
    }
  return U;
}


gen_op Shpul_U_plane(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso1,
					 double freq1, const std::string& iso2, double freq2, double time,
					 double fact, double phi)

	// Input		sys   : spin system
	//			BLK   : data block containing pulse function
	//			H     : acting static Hamiltonian
	// 			iso1  : isotope type
	// 			freq1 : pulse frequency (offset from isotope 1)
	// 			iso2  : isotope type
	// 			freq2 : pulse frequency (offset from isotope 2)
	// 			time  : shaped pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// 			phi   : pulse phase angle
	// 			sigma : density matrix

  {
  gen_op U, Ui;
  double timei, facti, lastfacti; 
  int steps;
  if(time < 0.0)
    {
    PulSherror(1, 1);                           // Error during pulse calculation
    PulShfatality(2); 				// Negative pulse length
    }
  else if(time == 0.0)
    U = I_gen_op(H.get_basis());		// Unity Operator if time = 0
  else
    {
    steps = BLK.size();				// Number of steps
    timei = time/double(steps-1);		// Compute time increment per step
    fact = fact/Re((BLK.sum() - BLK(0)));	// Normalize the factor to BLK
    U = I_gen_op(H.get_basis());		// Unity Operator initially
    lastfacti = 0;
    for(int i=1; i<steps; i++)
      {
      facti = Re(BLK(i))*fact*(steps-1);	// Scaling factor
      if(facti != lastfacti)					// Pulse proagator this waveform section
        Ui = Spul_U_plane(sys, H, iso1, freq1, iso2, freq2, timei, facti, phi);
      U = Ui * U;
      lastfacti = facti;
      }
    }
  return U;
  }

#endif 						// PulseShp.cc

