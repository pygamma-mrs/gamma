/* PulseS.h *****************************************************-*-c++-*-
**									**
**	                       G A M M A				**
**								 	**
**	Rectangular Pulses 			   Interface 		**
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

#ifndef   PulseS_h_			// Is file already included?
#  define PulseS_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/SpinSys.h>		// Include base spin system
#include <HSLib/GenOp.h>		// Include general operator

// ____________________________________________________________________________
// i                       Soft Pulse Error Handling
// ____________________________________________________________________________


void PulSerror(int eidx, int noret=0);

        // Input                eidx    : Error index
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error Message Output

     
void volatile PulSfatality(int eidx);
 
        // Input                none :
        // Input                eidx    : Error index
 

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

MSVCDLL gen_op Sxpuls(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
       const std::string& iso, double freq=0.0, double tp=1.e-5, double theta=90.0);

MSVCDLL gen_op SxpulsB(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
      const std::string& iso, double freq=0.0, double tp=1.e-5, double gamB1=2.5e4);

// ************************* Pulse Along the Y Axis ***************************

MSVCDLL gen_op Sypuls(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
       const std::string& iso, double freq=0.0, double tp=1.e-5, double theta=90.0);

MSVCDLL gen_op SypulsB(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
      const std::string& iso, double freq=0.0, double tp=1.e-5, double gamB1=2.5e4);

// ************* Pulse In the XY Plane On Axis Phi Degrees From X *************

MSVCDLL gen_op Sxypuls(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
                const std::string& iso,double freq=0.0,
                            double tp=1.e-5,double theta=90.0, double phi=0.0);

MSVCDLL gen_op SxypulsB(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
                const std::string& iso,double freq=0.0,
                            double tp=1.e-5,double gamB1=2.5e4,double phi=0.0);

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

MSVCDLL gen_op Sxpuls_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
                          double freq=0.0, double tp=1.e-5, double theta=90.0);

MSVCDLL gen_op SxpulsB_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
                         double freq=0.0, double tp=1.e-5, double gamB1=2.5e4);

// ----------------- Propagator For Pulse Along The Y Axis --------------------

MSVCDLL gen_op Sypuls_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
                          double freq=0.0, double tp=1.e-5, double theta=90.0);

MSVCDLL gen_op SypulsB_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
                         double freq=0.0, double tp=1.e-5, double gamB1=2.5e4);

// ----------------- Propagator For Pulse In The XY Plane ---------------------

MSVCDLL gen_op Sxypuls_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
          double freq=0.0, double tp=1.e-5, double theta=90.0, double phi=0.0);

MSVCDLL gen_op SxypulsB_U(const spin_sys& sys, const gen_op& H, const std::string& iso,
         double freq=0.0, double tp=1.e-5, double gamB1=2.5e4, double phi=0.0);

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
	   Output		U     : Propagator for application
*/

MSVCDLL gen_op Sxypuls(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
     const std::string& iso1, double freq1, const std::string& iso2, double freq2,
                           double tp=1.e-5, double theta=90.0, double phi=0.0);

MSVCDLL gen_op SxypulsB(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
     const std::string& iso1, double freq1, const std::string& iso2, double freq2,
                          double tp=1.e-5, double gamB1=2.5e4, double phi=0.0);

//______________________ Soft xy-Pulse Propagators _____________________


MSVCDLL gen_op Sxypuls_U(const spin_sys& sys, const gen_op& H,
       const std::string& iso1, double freq1, const std::string& iso2, double freq2,
                           double tp=1.e-5, double theta=90.0, double phi=0.0);

MSVCDLL gen_op SxypulsB_U(const spin_sys& sys, const gen_op& H,
	const std::string& iso1, double freq1, const std::string& iso2, double freq2,
                          double tp=1.e-5, double gamB1=2.5e4, double phi=0.0);

// ____________________________________________________________________________
// D                Generic Rectangular Pulse Functions 
// ____________________________________________________________________________

// _____________________ Pulse Along a Specific Axis ____________________

// _______ Soft Pulse on Coordinate Axis Density Matrix Evolution _______

MSVCDLL gen_op Spul_axis(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
            const std::string& iso, double freq, double tp, double fact, char axis);

MSVCDLL gen_op Spul_U_axis(const spin_sys& sys, const gen_op& H, const std::string& iso,
   			       double freq, double tp, double fact, char axis);

MSVCDLL gen_op Spul_plane(const spin_sys& sys, const gen_op& sigma, const gen_op& H, 
           const std::string& iso, double freq, double tp, double fact, double phi);

MSVCDLL gen_op Spul_plane(const spin_sys& sys, const gen_op& sigma, const gen_op& H,
            const std::string& iso1, double freq1, const std::string& iso2,
                             double freq2, double tp, double fact, double phi);

//__________________ Soft Pulse in xy-Plane Propagators _________________

MSVCDLL gen_op Spul_U_plane(const spin_sys& sys, const gen_op& H, const std::string& iso,
   			      double freq, double tp, double fact, double phi);

MSVCDLL gen_op Spul_U_plane(const spin_sys& sys, const gen_op& H, const std::string& iso1,
	 double freq1, const std::string& iso2, double freq2, 
                                           double tp, double fact, double phi);

 
#endif 								// PulseS.h

