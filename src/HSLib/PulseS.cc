/* PulseS.cc*****************************************************-*-c++-*-
**									**
**	                       G A M M A				**
**								 	**
**	NMR Library Soft Pulses                   Implementation	**
**						        	 	**
**      Copyright (c) 1990, 1991, 1992, 1993                            **
**      Scott Smith                                                     **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                   			**
** Description                                                  	**
**                                                              	**
** This module supplies functions that allow for facile use of  	**
** rectangular pulses in GAMMA.  There are three categories of 		**
** functions: Those that evolve the density operator, those that 	**
** return propagators, and those that help in pulse I/O.  Users should  **
** bear in mind the mathematical limitations of these pulses.  Simply   **
** put, pulses on a particular isotope channel ONLY influence spins of  **
** that isotope type.  These functions fail treating spins of differing	**
** type that have similar resonance frequencies.  That would occur when	**
** treating systems in very weak Bo fields and/or when two isotope 	**
** types happen to have similar gyromagntic ratios. Most of the time	**
** this isn't something to bother over! Working with typical spin types **
** (e.g. 1H, 13C, 19F, 14N, 2H, ....) at modern spectrometer fields	**
** (e.g. 1H resonance in MHz region) will not normally be any trouble   **
** whatsoever so use these function freely (& don't read on)            ** 
**                                                              	**
** You want to treat exactly this situation?  Then I hope you grasp the	**
** whole rotating frame reasoning in the first place......i.e. how the	**
** Liouville equation is solved in the rotation frame under the static	**
** and rf-field Hamiltonians combined.  If so, perhaps this will help.	**
** You just need to work in a single rotating frame, at least a single  **
** rotating frame must be used for the overlapping spin types. First	**
** set up your active Hamiltonian so that it properly reflects relative **
** frequencies of the overlapping spins. The apply the rf-field by	**
** mimicking one of the funcitons below, but make sure the "Fx" spin    **
** operator is active for both spin types (see the two channel pulse	**
** functions) That should do it.  In extreme cases just work in the	**
** lab frame and hope roundoff doesn't become a bother.			**
**									**
*************************************************************************/

#ifndef _PulseS_cc_			// Is file already included?
#define _PulseS_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <HSLib/PulseS.h>		// Include the header files
#include <HSLib/SpinOp.h>		// Knowledge of spin operators
#include <HSLib/SpinOpCmp.h>		// Knowledge of composite spin ops.
#include <HSLib/SpinOpRot.h>		// Knowledge of spin rotation ops.
#include <HSLib/PulseI.h>		// Knowledge of ideal pulses
#include <HSLib/HSprop.h>		// Knowledge of propagators
#include <Basics/Gutils.h>		// Include GAMMA error handling


// ____________________________________________________________________________
// i                       Soft Pulse Error Handling
// ____________________________________________________________________________


void PulSerror(int eidx, int noret)

        // Input                eidx    : Error index
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error Message Output

  {
  std::string hdr("Soft Pulse");
  std::string msg;
  switch(eidx)
    {
    case 1: msg = "Error During Soft Pulse Computation";
            GAMMAerror(hdr, msg, noret); break; //                         (1)
    case 2: msg = "Negative Pulse Lengths Are Disallowed";
            GAMMAerror(hdr, msg, noret); break; //                         (2)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }
     
void volatile PulSfatality(int eidx)
 
        // Input                none :
        // Input                eidx    : Error index
 
  {
  PulSerror(eidx, 1);				// Output error message
  if(eidx) PulSerror(0);			// Now output it fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// A        Rectangular Pulses Acting Directly on Density Operators
// ____________________________________________________________________________

/* The functions below apply a rectangular pulse of specified rotation 
   to an input density operator. The pulse angle of rotation beta, (for a spin
   exactly on resonance with the pulse rf frequency) is given by

                 beta = gamB1(Hz)*tp(sec)*360(deg/cycle)

   where gamB1 relates to the rf-field strength (gamma*B1) and tp is the length
   of time the pulse is applied.  Any two of the parameters {beta,gamB1,tp}
   suffice to set up the pulse strength/length and function overloads exist
   for specification of pulses with either {tp,beta} or {tp,gamB1}.  Note that
   since these pulses are of finite length they MUST be applied in the rotating
   frame of the pulse rf-field.  The offset of the pulse from the current
   rotating frame of the input density operator is given by the input argument
   freq. The functions themselves will adjust for the rf-field offset.  The
   density operator will advance by the time the pulse is applied.  Unlike
   ideal pulses, there are no function overloads which allow for spin selective
   or spin set selective pulses herein.  Pulse selectivity is entirely dictated
   by the channel the pulse is applied on, the frequency it has, and the length
   and strength of the pulse. 

	   Input		sys   : A (base) spin system
	   			sigma : Current density operator
				H     : Aciive isotropic Hamiltonian (Ho)
				iso   : Isotope channel (1H, 13C,....)
				freq  : Pulse offset (from sigma rot. frame)
				tp    : Pulse length (sec)
	   			beta  : Pulse rotation angle (degrees)
                             or gamB1 : Pulse strength (Hz)
				phi   : Pulse phase angle (degrees)
	   Output		sigma : Density operator after application
	  				of rectangular pulse on channel iso  */

// ************************* Pulse Along The X Axis ***************************

gen_op Sxpuls(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
                       const std::string& iso, double freq, double tp, double theta)
  {
  if(!tp) return Ixpuls(sys,sigma,iso,theta);	// An Ideal pulse if tp=0
  double fact = theta/(360*tp);
  return Spul_axis(sys, sigma, H, iso, freq, tp, fact, 'x');
  }

gen_op SxpulsB(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
                       const std::string& iso, double freq, double tp, double gamB1)
  { return Spul_axis(sys, sigma, H, iso, freq, tp, gamB1, 'x'); }


// ************************* Pulse Along the Y Axis ***************************

gen_op Sypuls(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
                       const std::string& iso, double freq, double tp, double theta)
  {
  if(!tp) return Iypuls(sys,sigma,iso,theta);	// An ideal pulse if tp=0
  double fact = theta/(360*tp);
  return Spul_axis(sys, sigma, H, iso, freq, tp, fact, 'y');
  }

gen_op SypulsB(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
                       const std::string& iso, double freq, double tp, double gamB1)
  { return Spul_axis(sys, sigma, H, iso, freq, tp, gamB1, 'y'); }

// ************* Pulse In the XY Plane On Axis Phi Degrees From X *************

gen_op Sxypuls(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
           const std::string& iso, double freq, double tp, double theta, double phi)
  {
  if(!tp) return Ixypuls(sys,sigma,iso,theta,phi);	// Ideal pulse if tp=0
  double fact = theta/(360.0*tp);
  return Spul_plane(sys, sigma, H, iso, freq, tp, fact, phi);
  }

gen_op SxypulsB(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
           const std::string& iso, double freq, double tp, double gamB1, double phi)
{ return Spul_plane(sys, sigma, H, iso, freq, tp, gamB1, phi); }

// ____________________________________________________________________________
// B      Rectangular Pulses Propagators Which Act on Density Operators
// ____________________________________________________________________________

/* The functions below return rectangular pulse propagators of specified
   rotation and offset.  The pulse angle of rotation beta, (for a spin
   exactly on resonance with the pulse rf frequency) is given by

                 beta = gamB1(Hz)*tp(sec)*360(deg/cycle)

   where gamB1 relates to the rf-field strength (gamma*B1) and tp is the length
   of time the pulse is applied.  Any two of the parameters {beta,gamB1,tp}
   suffice to set up the pulse strength/length and function overloads exist
   for specification of pulses with either {tp,beta} or {tp,gamB1}.  Note that
   since these pulses are of finite length they MUST be applied in the rotating
   frame of the pulse rf-field.  The offset of the pulse from the rotating
   frame of density operator to whcih the propagator applies is given by the
   input argument freq. The functions themselves will adjust for the rf-field 
   offset.  The density operator acted on by the propagator will advance by the
   time the pulse is applied.  Unlike ideal pulses, there are no function 
   overloads which allow for spin selective or spin set selective pulses here.
   Pulse selectivity is entirely dictated by the channel the pulse is applied
   on, the frequency it has, and the length and strength of the pulse. 

	   Input		sys   : A (base) spin system
				H     : Active static Hamiltonian (Ho)
				iso   : Isotope channel (1H, 13C,....)
				freq  : Pulse offset (from sigma rot. frame)
				tp    : Pulse length (sec)
	   			beta  : Pulse rotation angle (degrees)
                             or gamB1 : Pulse strength (Hz)
				phi   : Pulse phase angle (degrees)
	   Output		U     : Propagator for application
	  				of rectangular pulse on channel iso
*/

// ----------------- Propagator For Pulse Along The X Axis --------------------

gen_op Sxpuls_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
                                          double freq, double tp, double theta)
  {
  if(!tp) return Rx(sys,iso, theta);
  double fact = theta/(360*tp);
  return Spul_U_axis(sys, H, iso, freq, tp, fact, 'x');
  }

gen_op SxpulsB_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
                                          double freq, double tp, double gamB1)
  { return Spul_U_axis(sys, H, iso, freq, tp, gamB1, 'x'); }

// ----------------- Propagator For Pulse Along The Y Axis --------------------

gen_op Sypuls_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
                                          double freq, double tp, double theta)
  {
  if(!tp) return Ry(sys,iso, theta);
  double fact = theta/(360*tp);
  return Spul_U_axis(sys, H, iso, freq, tp, fact, 'y');
  }

gen_op SypulsB_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
                                          double freq, double tp, double gamB1)
  { return Spul_U_axis(sys, H, iso, freq, tp, gamB1, 'y'); }


// ----------------- Propagator For Pulse In The XY Plane ---------------------

gen_op Sxypuls_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
                              double freq, double tp, double theta, double phi)
  {
  if(!tp) return Rxy(sys,iso,theta,phi);
  double fact = theta/(360.0*tp);
  return  Spul_U_plane(sys, H, iso, freq, tp, fact, phi);
  }

gen_op SxypulsB_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
                              double freq, double tp, double gamB1, double phi)
  { return Spul_U_plane(sys, H, iso, freq, tp, gamB1, phi); }

// ____________________________________________________________________________
// C Two Channel Rectangular Pulses (Propagators & Density Operator Evolution)
// ____________________________________________________________________________

/* The functions below deal with rectangular pulses appplied to two channels
   simultaneously.  This it NOT STRICTLY CORRECT mathematically.  A tacit
   assumption is made that the Larmor frequencies of the two spin types hit
   by the pulses are far enough apart that there will be no overlap in the
   power spectrum of the two pulses.  Although often the case, it will fail
   when two isotope types have similar gyromagnetic ratios and/or the field
   strength of the spectrometer is small.  It is up to the user to keep this
   in mind when using these functions!  The functions don't bother checking
   if the pulse is valid for your system as that depends on the two isotope
   types, the Bo field strength, the pulse length, and rf-field strengths
   applied.

   The the pulse length is identical on the channels, so the pulse rotation
   (on resonance) for channel i is given by

                 beta = gamB1(Hz)*tp(sec)*360(deg/cycle)
                     i       i

   where gamB1 relates to the rf-field strength on channel i [gamma(i)*B1]
   and tp is the length of time the pulse is applied. 

	   Input	sys      : A (base) spin system
			H        : Active static Hamiltonian (Ho)
			iso(1,2) : Isotope channel (1H, 13C,....)
			freq(1,2): Pulse offset (from sigma rot. frame)
			tp       : Pulse length (sec)
	   		beta     : Pulse rotation angle (degrees)
			phi      : Pulse phase angle (degrees)
	   Output	U	 : Propagator for application
*/

gen_op Sxypuls(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
     const std::string& iso1, double freq1, const std::string& iso2, double freq2,
	   	                           double tp, double theta, double phi)
  {
  if(!tp)					// If no pulse length then
    {						// use an ideal pulse!
    flagvec flags=sys.GetFlags(iso1,1,0);	// Spin flags with iso1 TRUE
    for(int i=0; i<sys.spins(); i++)		// Also set flags for spins
      if(sys.symbol(i) == iso2) flags[i]=1;	// of type iso2 to TRUE
    return Ixypuls(sys,sigma,flags,theta,phi);	// Ideal pulse if tp=0
    }
  double fact= theta/(360.0*tp);
  return Spul_plane(sys,sigma,H,iso1,freq1,iso2,freq2,tp,fact,phi);
  }


gen_op SxypulsB(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
     const std::string& iso1, double freq1, const std::string& iso2, double freq2,
                                           double tp, double gamB1, double phi)
  { return Spul_plane(sys,sigma,H,iso1,freq1,iso2,freq2,tp,gamB1,phi); }

//______________________ Soft xy-Pulse Propagators _____________________


gen_op Sxypuls_U(const spin_sys& sys, const gen_op& H,
       const std::string& iso1, double freq1, const std::string& iso2, double freq2,
		                           double tp, double theta, double phi)
  {
  if(!tp)
    {
    flagvec flags=sys.GetFlags(iso1,1,0);	// Spin flags with iso1 TRUE
    for(int i=0; i<sys.spins(); i++)		// Also set flags for spins
      if(sys.symbol(i) == iso2) flags[i]=1;	// of type iso2 to TRUE
    return Rxy(sys, flags, theta, phi);		// Rotation operator if tp=0
    }
  double fact = theta/(360.0*tp);
  return Spul_U_plane(sys,H,iso1,freq1,iso2,freq2,tp,fact,phi);
  }

gen_op SxypulsB_U(const spin_sys& sys, const gen_op& H,
	const std::string& iso1, double freq1, const std::string& iso2, double freq2,
                                           double tp, double gamB1, double phi)
  { return Spul_U_plane(sys,H,iso1,freq1,iso2,freq2,tp,gamB1,phi); }

// ----------------------------------------------------------------------
// D                Generic Rectangular Pulse Functions 
// ----------------------------------------------------------------------

// _____________________ Pulse Along a Specific Axis ____________________

// _______ Soft Pulse on Coordinate Axis Density Matrix Evolution _______

gen_op Spul_axis(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
             const std::string& iso, double freq, double tp, double fact, char axis)

	// Input		sys   : spin system
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			tp    : soft pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// Output		sigma : density matrix following application
	//				of a soft pulse.

  {
  if(!tp) return sigma;				// No effect from 0 length pulse
  if(tp < 0.0)
    {
    PulSerror(1, 1);				// Error during soft pulse
    PulSfatality(110);				// Negative pulse length
    }
  gen_op Heff = H - freq*Fz(sys, iso);		// rotating @ rf-freq "freq" 
  switch(axis)					// Add rf-field to rot. Ham.
    {
    case 'x':
    default:  Heff += fact*Fx(sys,iso); break;// Field along the x-axis
    case 'y': Heff += fact*Fy(sys,iso); break;// Field along the y-axis
    }
  gen_op sigmap = evolve(sigma, Heff, tp);	// Evolve through soft pulse
  Heff  = freq*Fz(sys, iso);			// To back out of rotframe
  return evolve(sigmap, Heff, tp);		// Evolve out of rotfame
  }

//_______________ Soft Pulse on Coordinate Axis Propagator ______________

gen_op Spul_U_axis(const spin_sys& sys, const gen_op& H, const std::string& iso,
			double freq, double tp, double fact, char axis)

	// Input		sys   : spin system
	//			H     : active static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			tp    : soft pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	//			axis  : axis along which the pulse is applied
	// Output		U     : pulse propagator for the application
	//				of a soft pulse.

  {
  gen_op U1, Heff;
  if(tp < 0.0)
    {
    PulSerror(1, 1);				// Error during soft pulse
    PulSfatality(110);				// Negative pulse length
    }
  if(!tp)
    {
    matrix mx = matrix(H.dim(),H.dim(),i_matrix_type);
    return gen_op(mx);				// Simply the Unity Operator
//    U1 = gen_op(mx, sys.get_basis());		// Simply the Unity Operator
    }
  Heff = H - freq*Fz(sys, iso);			// H rotating @ rf freq. "freq" 
  switch(axis)					// Add rf-field to rot. Ham.
    {
    case 'x':
    default:  Heff += fact*Fx(sys,iso); break;	// Field along the x-axis
    case 'y': Heff += fact*Fy(sys,iso); break;	// Field along the y-axis
    }
  U1 = prop(Heff, tp);				// Soft pulse prop. in rotframe
  Heff = freq*Fz(sys, iso);			// Ham. to back out of rotframe
  U1 &= prop(Heff, tp);				// Total soft pulse propagator
  return U1;
  }

// ____________________ Pulse Applied in the xy-Plane ___________________

// ________ Soft Pulse in the xy-Plane Density Matrix Evolution _________

gen_op Spul_plane(const spin_sys& sys, const gen_op& sigma, const gen_op& H, const std::string& iso,
			 double freq, double tp, double fact, double phi)

	// Input		sys   : spin system
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			tp    : soft pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// 			phi   : pulse phase angle
	// Output		sigma : density matrix following application
	//				of a soft pulse.

  {
  if(!tp) return sigma;			// Nothing for 0 length pulse
  if(tp < 0.0)
    {
    PulSerror(1, 1);			// Error during soft pulse
    PulSfatality(110);			// Negative pulse length
    }
  gen_op Heff = H - freq*Fz(sys, iso);	// Rotating @ rf-field frequency "freq" 
  Heff += fact*Fxy(sys, iso, phi);	// Add rf-field to rotating Hamiltonian
  gen_op sigmap = evolve(sigma,Heff,tp);// Evolve sigma through soft pulse
  Heff = freq*Fz(sys, iso);		// Back out of frame rotating at freq
  return evolve(sigmap, Heff, tp); 	// Propagate the density operator
  }


gen_op Spul_plane(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
            const std::string& iso1, double freq1, const std::string& iso2,
                              double freq2, double tp, double fact, double phi)

	// Input		sys   : spin system
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso1  : isotope type
	// 			freq1 : pulse frequency (offset from isotope 1)
	// 			iso2  : isotope type
	// 			freq2 : pulse frequency (offset from isotope 2)
	// 			tp    : soft pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// 			phi   : pulse phase angle
	// Output		sigma : density matrix following application
	//				of a soft pulse.

  {
  if(!tp) return sigma;				// Nothing if 0 length pulse
  if(tp < 0.0)					// Problems in pulse length
    {						// is of negative time.
    PulSerror(1,1);
    PulSfatality(110);
    }
// sosi
  gen_op Hz = freq1*Fz(sys,iso1)		// For switch into two rf-field
            + freq2*Fz(sys,iso2);		// rotating frames
  gen_op Heff = H - Hz;				// H in the rotating frames
  flagvec flags=sys.GetFlags(iso1,1,0);	// Spin flags with iso1 TRUE
  for(int i=0; i<sys.spins(); i++)		// Also set flags for spins
  if(sys.symbol(i) == iso2) flags[i]=1;		// of type iso2 to TRUE
  Heff += fact*Fxy(sys,flags,phi);		// Add rf-field to rot. Ham.
// sosi OK to there
std::cout << "\n\nHeff in Function\n\n" << Heff << "\n\n";
  gen_op sigmap = evolve(sigma, Heff, tp);	// Evolve through soft pulse
  Heff  = freq1*Fz(sys, iso1);			// To back out of rotframe
  Heff -= freq2*Fz(sys, iso2);			// of both isotope channels
  return evolve(sigmap, -Hz, tp);		// Back out of rotating frames
  }

//__________________ Soft Pulse in xy-Plane Propagators _________________

gen_op Spul_U_plane(const spin_sys& sys, const gen_op& H, const std::string& iso,
			 double freq, double tp, double fact, double phi)

	// Input		sys   : spin system
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			tp    : soft pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// 			phi   : pulse phase angle
	// Output		U     : propagator for a soft pulse

{
gen_op U, Heff;
matrix mx;
if(tp<0.0)				// Can't have a pulse with 
  {					// a negative length
  PulSerror(1, 1);			// Error during soft pulse
  PulSfatality(110);			// Negative pulse length
  }
else if(tp == 0.0)			// If the pulse is of 0 length
  {					// then it won't affect anything
  int hs = H.dim();			// This is the Hilbert space size
  mx = matrix(hs,hs,i_matrix_type);	// An identity matrix
  U = gen_op(mx);			// Return the unity operator
//    U = gen_op(mx, sys.get_basis());		// Simply the Unity Operator
  }
else
  {
  Heff = H;					// Start with input Hamiltonian
  if(freq)					// Adjust Hamiltonian to frame
    Heff -= freq*Fz(sys, iso);			// rotating at rf-field frequency "freq" 
  Heff += fact*Fxy(sys, iso, phi);		// Add rf-field to rotating Hamiltonian
  U = prop(Heff, tp);				// Propagator through soft pulse
  if(freq)
    {
    Heff = freq*Fz(sys, iso);			// Back out of frame rotating at freq
    U &= prop(Heff, tp);			// Total soft pulse propagator
    }
  }
return U;
}


gen_op Spul_U_plane(const spin_sys& sys, const gen_op& H, const std::string& iso1,
	 double freq1, const std::string& iso2, double freq2, 
                                            double tp, double fact, double phi)

	// Input		sys   : spin system
	//			H     : acting static Hamiltonian
	// 			iso1  : isotope type
	// 			freq1 : pulse frequency (offset from isotope 1)
	// 			iso2  : isotope type
	// 			freq2 : pulse frequency (offset from isotope 2)
	// 			tp    : soft pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// 			phi   : pulse phase angle
	// 			sigma : density matrix

  {
  if(tp < 0.0)
    {
    PulSerror(1,1);				// Error during soft pulse
    PulSfatality(110);				// Negative pulse length
    }
  if(!tp)
    {
    int dim = H.dim();
    return gen_op(matrix(dim,dim,i_matrix_type));// Return unity operator
//    rturn gen_op(mx, sys.get_basis());	// Simply the Unity Operator
    }
  gen_op Heff = H-freq1*Fz(sys, iso1);		// Rotating @ rf-freq. "freq1" 
  Heff        -=  freq2*Fz(sys, iso2);		// & @ rf-freq. "freq2" 
  flagvec flags=sys.GetFlags(iso1,1,0);	// Spin flags with iso1 TRUE
  for(int i=0; i<sys.spins(); i++)		// Also set flags for spins
    if(sys.symbol(i) == iso2) flags[i]=1;	// of type iso2 to TRUE
  Heff += fact*Fxy(sys, flags, phi);		// Add rf-field to rot. Ham.
  gen_op U = prop(Heff, tp);			// Propagator through soft pulse
  Heff = freq1*Fz(sys, iso1);			// Back out of rot. frame @
  Heff -= freq2*Fz(sys, iso2);			// freq1 & freq2 
  U &= prop(Heff, tp);				// Total soft pulse propagator
  return U;
  }


#endif 							// PulseS.cc

