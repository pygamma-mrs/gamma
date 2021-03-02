/* PulseShp.h ***********************************-*-c++-*-
**							**
**	               G A M M A			**
**						 	**
**	NMR Library                    Definitions      **
**						 	**
**	Copyright (c) 1990			 	**
**	Scott Smith				 	**
**	Eidgenoessische Technische Hochschule	 	**
**	Labor fuer physikalische Chemie		 	**
**	8092 Zurich / Switzerland		 	**
**						 	**
**      $Header:
**						 	**
** Modifications:					**
**						 	**
*********************************************************/

/*********************************************************
**						 	**
** 	Description				 	**
**						 	**
**	The folowing functions implement shaped pulses  **
**						 	**
*********************************************************/

#ifndef   Nmr_shpuls_h_			// Is file already included?
#  define Nmr_shpuls_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/GenOp.h>
#include <Matrix/row_vector.h>
#include <HSLib/SpinSys.h>


// ____________________________________________________________________________
// i                       Shaped Pulse Error Handling
// ____________________________________________________________________________
 

void PulSherror(int eidx, int noret);

        // Input                eidx    : Error index
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error Message Output
     
 
void volatile PulShfatality(int eidx);
 
        // Input                none :
        // Input                eidx    : Error index

// ______________________________________________________________________
// ************** Functions with the Pulse Angle Specified **************
// ______________________________________________________________________

// ______________________ Pulse Along the x_axis ________________________

// _______________ Shaped x-Pulse Density Matrix Evolution ________________

MSVCDLL gen_op Shxpuls(const spin_sys& sys, row_vector& BLK, gen_op sigma, gen_op& H,
	 const std::string& iso, double freq=0.0, double time=1.0e-5, double theta=90.0);

	// Input		sys   : spin system
	// 			BLK   : data block
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			theta : pulse angle (for a spin at resonance)
	// Output		sigma : density matrix following application
	//				of the following soft pulse -
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length = time
	//				pulse offset = freq (from isotope carrier)
	//				pulse axis   = x
	//				pulse angle  = theta (for all pulse lengths)

	// Note			      : only a spin on resonance will experience
	//				an absolute pulse angle of theta degrees

//_______________________ Shaped x-Pulse Propagator ______________________

//gen_op Shxpuls_U(const spin_sys &sys, row_vector &BLK, gen_op &H, char *iso,
MSVCDLL gen_op Shxpuls_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
		double freq=0.0, double time=1.0e-5, double theta=90.0);

	// Input		sys   : spin system
	// 			BLK   : data block
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			theta : pulse angle (for a spin at resonance)
	// Output		U     : propagator for the following soft pulse
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length = time
	//				pulse offset = freq (from isotope carrier)
	//				pulse axis   = x
	//				pulse angle  = theta (for all pulse lengths)

	// Note			      : only a spin on resonance will experience
	//				an absolute pulse angle of theta degrees


// ______________________ Pulse Along the y_axis ________________________

// _______________ Shaped y-Pulse Density Matrix Evolution ________________

MSVCDLL gen_op Shypuls(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
//	 char *iso, double freq=0.0, double time=1.0e-5, double theta=90.0);
	 const std::string& iso, double freq=0.0, double time=1.0e-5, double theta=90.0);

	// Input		sys   : spin system
	// 			BLK   : data block
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			theta : pulse angle
	// Output		sigma : density matrix following application
	//				of the following soft pulse -
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length = time
	//				pulse offset = freq (from isotope carrier)
	//				pulse axis   = y
	//				pulse angle  = theta (for all pulse lengths)

	// Note			      : on a spin on resonance will experience
	//				an absolute pulse angle of theta degrees

//_______________________ Shaped y-Pulse Propagator ______________________

//gen_op Shypuls_U(const spin_sys &sys, row_vector &BLK, gen_op &H, char *iso,
MSVCDLL gen_op Shypuls_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
		double freq=0.0, double time=1.0e-5, double theta=90.0);

	// Input		sys   : spin system
	// 			BLK   : data block
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			theta : pulse angle (for a spin at resonance)
	// Output		U     : propagator for the following soft pulse
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length = time
	//				pulse offset = freq (from isotope carrier)
	//				pulse axis   = y
	//				pulse angle  = theta (for all pulse lengths)

	// Note			      : only a spin on resonance will experience
	//				an absolute pulse angle of theta degrees


// ______________________ Pulse in the xy_plane ________________________

// _______________ Shaped xy-Pulse Density Matrix Evolution _______________

MSVCDLL gen_op Shxypuls(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
		 const std::string& iso, double freq=0.0, double time=1.0e-5,
					 double theta=90.0, double phi=0.0);
//		 	char *iso, double freq=0.0, double time=1.0e-5,
//					 double theta=90.0, double phi=0.0);

	// Input		sys   : spin system
	// 			BLK   : data block
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			theta : pulse angle
	// 			phi   : pulse phase angle
	// Output		sigma : density matrix following application
	//				of the following soft pulse -
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length = time
	//				pulse offset = freq (from isotope carrier)
	//				pulse axis   = y
	//				pulse angle  = theta (for all pulse lengths)

	// Note			      : on a spin on resonance will experience
	//				an absolute pulse angle of theta degrees


MSVCDLL gen_op Shxypuls(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
		 const std::string& iso1, double freq1, const std::string& iso2, double freq2,
			 double time=1.0e-5, double theta=90.0, double phi=0.0);
//		 char *iso1, double freq1, char*iso2, double freq2,
//			 double time=1.0e-5, double theta=90.0, double phi=0.0);

	// Input		sys   : spin system
	// 			BLK   : data block
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso1  : first isotope type
	// 			freq1 : pulse frequency (offset from isotope 1)
	// 			iso2  : second isotope type
	// 			freq2 : pulse frequency (offset from isotope 2)
	// 			time  : soft pulse length
	// 			theta : pulse angle
	// 			phi   : pulse phase angle
	// Output		sigma : density matrix following application
	//				of the following soft pulse -
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length = time
	//				pulse offset = freq1 (isotope 1 carrier)
	//				               freq2 (isotope 2 carrier)
	//				pulse axis   = phi degrees from x
	//				pulse angle  = theta (for all pulse lengths)

	// Note			      : on a spin on resonance will experience
	//				an absolute pulse angle of theta degrees

//______________________ Shaped xy-Pulse Propagators _____________________

//gen_op Shxypuls_U(const spin_sys &sys, row_vector &BLK, gen_op &H, char *iso,
MSVCDLL gen_op Shxypuls_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
	double freq=0.0, double time=1.0e-5, double theta=90.0, double phi=0.0);

	// Input		sys   : spin system
	// 			BLK   : data block
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			theta : pulse angle
	// 			phi   : pulse phase angle
	// Output		U     : propagator for the following soft pulse -
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length = time
	//				pulse offset = freq (from isotope carrier)
	//				pulse axis   = y
	//				pulse angle  = theta (for all pulse lengths)

	// Note			      : on a spin on resonance will experience
	//				an absolute pulse angle of theta degrees


MSVCDLL gen_op Shxypuls_U(const spin_sys &sys, row_vector &BLK, gen_op &H,
		 const std::string& iso1, double freq1, const std::string& iso2, double freq2,
			 double time=1.0e-5, double theta=90.0, double phi=0.0);
//		 char *iso1, double freq1, char*iso2, double freq2,
//			 double time=1.0e-5, double theta=90.0, double phi=0.0);

	// Input		sys   : spin system
	// 			BLK   : data block
	//			H     : acting static Hamiltonian
	// 			iso1  : first isotope type
	// 			freq1 : pulse frequency (offset from isotope 1)
	// 			iso2  : second isotope type
	// 			freq2 : pulse frequency (offset from isotope 2)
	// 			time  : soft pulse length
	// 			theta : pulse angle
	// 			phi   : pulse phase angle
	// Output		U     : propagator for the following soft pulse -
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length = time
	//				pulse offset = freq1 (isotope 1 carrier)
	//				               freq2 (isotope 2 carrier)
	//				pulse axis   = phi degrees from x
	//				pulse angle  = theta (for all pulse lengths)

	// Note			      : on a spin on resonance will experience
	//				an absolute pulse angle of theta degrees

// ______________________________________________________________________
// ************* Functions with the Field Strength Specified ************
// ______________________________________________________________________

// ______________________ Pulse Along the x_axis ________________________

// _______________ Shaped x-Pulse Density Matrix Evolution ________________

MSVCDLL gen_op ShxpulsB(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
	 const std::string& iso, double freq=0.0, double time=1.0e-5, double gamB1=2.5e4);
//	 char *iso, double freq=0.0, double time=1.0e-5, double gamB1=2.5e4);

	// Input		sys   : spin system
	// 			BLK   : data block
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// Output		sigma : density matrix following application
	//				of the following soft pulse -
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length   = time
	//				pulse offset   = freq (from isotope carrier)
	//				pulse axis     = x
	//				pulse angle    = gamB1*time*360 (in degrees)
	//				pulse strength = gamB1/gamma (isotope type iso)
	//
	// Note			      : only a spin on resonance will experience
	//				an absolute pulse angle of theta degrees

//_______________________ Shaped x-Pulse Propagator ______________________

MSVCDLL gen_op ShxpulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
			double freq=0.0, double time=1.0e-5, double gamB1=2.5e4);
//MSVCDLL gen_op ShxpulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H, char *iso,
//			double freq=0.0, double time=1.0e-5, double gamB1=2.5e4);

	// Input		sys   : spin system
	// 			BLK   : data block
	//			H     : normal Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// Output		U     : propagator for the following soft pulse
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length   = time
	//				pulse offset   = freq (from isotope carrier)
	//				pulse axis     = x
	//				pulse angle    = gamB1*time*360 (in degrees)
	//				pulse strength = gamB1/gamma (isotope type iso)
	//
	// Note			      : only a spin on resonance will experience
	//				an absolute pulse angle of theta degrees


// ______________________ Pulse Along the y_axis ________________________

// _______________ Shaped y-Pulse Density Matrix Evolution ________________

MSVCDLL gen_op ShypulsB(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
	 const std::string& iso, double freq=0.0, double time=1.0e-5, double gamB1=2.5e4);

	// Input		sys   : spin system
	// 			BLK   : data block
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// Output		sigma : density matrix following application
	//				of the following soft pulse -
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length   = time
	//				pulse offset   = freq (from isotope carrier)
	//				pulse axis     = y
	//				pulse angle    = gamB1*time*360 (in degrees)
	//				pulse strength = gamB1/gamma (isotope type iso)
	//
	// Note			      : only a spin on resonance will experience
	//				an absolute pulse angle of theta degrees

//_______________________ Shaped y-Pulse Propagator ______________________

//gen_op ShypulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H, char *iso,
MSVCDLL gen_op ShypulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
			double freq=0.0, double time=1.0e-5, double gamB1=2.5e4);

	// Input		sys   : spin system
	// 			BLK   : data block
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// Output		U     : propagator for the following soft pulse
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length   = time
	//				pulse offset   = freq (from isotope carrier)
	//				pulse axis     = y
	//				pulse angle    = gamB1*time*360 (in degrees)
	//				pulse strength = gamB1/gamma (isotope type iso)
	//
	// Note			      : only a spin on resonance will experience
	//				an absolute pulse angle of theta degrees


// ______________________ Pulse in the xy_plane ________________________

// _______________Shaped xy-Pulse Density Matrix Evolution _______________

MSVCDLL gen_op ShxypulsB(const spin_sys &sys, row_vector &BLK, gen_op sigma,
	 	gen_op &H, const std::string& iso, double freq=0.0,
			 double time=1.0e-5, double gamB1=2.5e4, double phi=0.0);
//	 	gen_op &H, char* iso, double freq=0.0,
//			 double time=1.0e-5, double gamB1=2.5e4, double phi=0.0);

	// Input		sys   : spin system
	// 			BLK   : data block
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// 			phi   : pulse phase angle
	// Output		sigma : density matrix following application
	//				of the following soft pulse -
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length   = time
	//				pulse offset   = freq (from isotope carrier)
	//				pulse axis     = x
	//				pulse angle    = gamB1*time*360 (in degrees)
	//				pulse phase    = phi (in degrees)
	//				pulse strength = gamB1/gamma (isotope type iso)
	//
	// Note			      : only a spin on resonance will experience
	//				an absolute pulse angle of theta degrees


MSVCDLL gen_op ShxypulsB(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
	 const std::string& iso1, double freq1, const std::string& iso2, double freq2,
			 double time=1.0e-5, double gamB1=2.5e4, double phi=0.0);
//	 char *iso1, double freq1, char* iso2, double freq2,
//			 double time=1.0e-5, double gamB1=2.5e4, double phi=0.0);

	// Input		sys   : spin system
	// 			BLK   : data block
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso1  : first isotope type
	// 			freq1 : pulse frequency (offset from isotope 1)
	// 			iso2  : second isotope type
	// 			freq2 : pulse frequency (offset from isotope 2)
	// 			time  : soft pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// 			phi   : pulse phase angle
	// Output		sigma : density matrix following application
	//				of the following soft pulse -
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length   = time
	//				pulse offset   = freq1 (from isotope 1)
	//						 freq2 (from isotope 2)
	//				pulse axis     = x
	//				pulse angle    = gamB1*time*360 (in degrees)
	//				pulse phase    = phi (in degrees)
	//				pulse strength = gamB1/gamma (isotope type iso)
	//
	// Note			      : only a spin on resonance will experience
	//				an absolute pulse angle of theta degrees

//______________________ Shaped xy-Pulse Propagators ______________________

//gen_op ShxypulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H, char *iso,
MSVCDLL gen_op ShxypulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
	double freq=0.0, double time=1.0e-5, double gamB1=2.5e4, double phi=0.0);

	// Input		sys   : spin system
	// 			BLK   : data block
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// 			phi   : pulse phase angle
	// Output		U     : propagator for the following soft pulse -
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length   = time
	//				pulse offset   = freq (from isotope carrier)
	//				pulse axis     = x
	//				pulse angle    = gamB1*time*360 (in degrees)
	//				pulse phase    = phi (in degrees)
	//				pulse strength = gamB1/gamma (isotope type iso)
	//
	// Note			      : only a spin on resonance will experience
	//				an absolute pulse angle of theta degrees


MSVCDLL gen_op ShxypulsB_U(const spin_sys &sys, row_vector &BLK, gen_op &H,
	 const std::string& iso1, double freq1, const std::string& iso2, double freq2,
			 double time=1.0e-5, double gamB1=2.5e4, double phi=0.0);
//	 char *iso1, double freq1, char* iso2, double freq2,
//			 double time=1.0e-5, double gamB1=2.5e4, double phi=0.0);

	// Input		sys   : spin system
	// 			BLK   : data block
	//			H     : acting static Hamiltonian
	// 			iso1  : first isotope type
	// 			freq1 : pulse frequency (offset from isotope 1)
	// 			iso2  : second isotope type
	// 			freq2 : pulse frequency (offset from isotope 2)
	// 			time  : soft pulse length
	// 			gamB1 : field strength (gamma*B1) in Hertz
	// 			phi   : pulse phase angle
	// Output		U     : propagator for the following soft pulse -
	//
	//				    Shaped Pulse Parameters
	//                              ----------------------------------------
	//				pulse length   = time
	//				pulse offset   = freq1 (from isotope 1)
	//						 freq2 (from isotope 2)
	//				pulse axis     = x
	//				pulse angle    = gamB1*time*360 (in degrees)
	//				pulse phase    = phi (in degrees)
	//				pulse strength = gamB1/gamma (isotope type iso)
	//
	// Note			      : only a spin on resonance will experience
	//				an absolute pulse angle of theta degrees

// ----------------------------------------------------------------------
// ********************* Generic Shaped Pulse Functions *******************
// ----------------------------------------------------------------------

// _____________________ Pulse Along a Specific Axis ____________________

// _______ Shaped Pulse on Coordinate Axis Density Matrix Evolution _______

MSVCDLL gen_op Shpul_axis(const spin_sys &sys, row_vector &BLK, gen_op sigma, gen_op &H,
		 const std::string& iso, double freq, double time, double fact, char axis);
//		 char *iso, double freq, double time, double fact, char axis);

	// Input		sys   : spin system
	// 			BLK   : data block
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// Output		sigma : density matrix following application
	//				of a soft pulse.
	// Note			      : Accessed by Shqpuls functions where
	//				q = axis, e.g. x or y

//_______________ Shaped Pulse on Coordinate Axis Propagator ______________

//gen_op Shpul_U_axis(const spin_sys &sys, row_vector &BLK, gen_op &H, char *iso,
MSVCDLL gen_op Shpul_U_axis(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
			double freq, double time, double fact, char axis);

	// Input		sys   : spin system
	// 			BLK   : data block
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	//			axis  : axis along which the pulse is applied
	// Output		U     : pulse propagator for a soft pulse
	// Note			      : Accessed by Shqpuls functions where
	//				q = axis, e.g. x


// ____________________ Pulse Applied in the xy-Plane ___________________

// ________ Shaped Pulse in the xy-Plane Density Matrix Evolution _________

MSVCDLL gen_op Shpul_plane(const spin_sys &sys, row_vector &BLK, gen_op sigma,
	 gen_op &H, const std::string& iso, double freq, double time, double fact, double phi);
//	 gen_op &H, char *iso, double freq, double time, double fact, double phi);

	// Input		sys   : spin system
	// 			BLK   : data block
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// 			phi   : pulse phase angle
	// Output		sigma : density matrix following application
	//				of a soft pulse.
	// Note			      : Accessed by Shxypuls functions


MSVCDLL gen_op Shpul_plane(const spin_sys &sys, row_vector &BLK, gen_op sigma,
		 gen_op &H, const std::string& iso1, double freq1, const std::string& iso2,
			 double freq2, double time, double fact, double phi);
//		 gen_op &H, char *iso1, double freq1, char *iso2,
//			 double freq2, double time, double fact, double phi);

	// Input		sys   : spin system
	// 			BLK   : data block
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso1  : isotope type
	// 			freq1 : pulse frequency (offset for isotope 1)
	// 			iso2  : isotope type
	// 			freq2 : pulse frequency (offset for isotope 2)
	// 			time  : soft pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// 			phi   : pulse phase angle
	// Output		sigma : density matrix following application
	//				of a soft pulse.
	// Note			      : Accessed by Shxypuls functions

//__________________ Shaped Pulse in xy-Plane Propagators _________________

//gen_op Shpul_U_plane(const spin_sys &sys, row_vector &BLK, gen_op &H, char *iso,
MSVCDLL gen_op Shpul_U_plane(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso,
			 double freq, double time, double fact, double phi);

	// Input		sys   : spin system
	// 			BLK   : data block
	// 			sigma : density matrix
	//			H     : acting static Hamiltonian
	// 			iso   : isotope type
	// 			freq  : pulse frequency (offset)
	// 			time  : soft pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// 			phi   : pulse phase angle
	// Output		U     : soft pulse propagator
	// Note			      : Accessed by Shxypuls_U functions


//gen_op Shpul_U_plane(const spin_sys &sys, row_vector &BLK, gen_op &H, char *iso1,
//	double freq1, char *iso2, double freq2, double time, double fact, double phi);
MSVCDLL gen_op Shpul_U_plane(const spin_sys &sys, row_vector &BLK, gen_op &H, const std::string& iso1,
 	double freq1, const std::string& iso2, double freq2, double time, double fact, double phi);

	// Input		sys   : spin system
	// 			BLK   : data block
	//			H     : acting static Hamiltonian
	// 			iso1  : isotope type
	// 			freq1 : pulse frequency (offset for isotope 1)
	// 			iso2  : isotope type
	// 			freq2 : pulse frequency (offset for isotope 2)
	// 			time  : soft pulse length
	// 			fact  : either gamma*B1 or pulse angle scaled
	// 			phi   : pulse phase angle
	// Output		U     : soft pulse propagator
	// Note			      : Accessed by Shxypuls_U functions
 
#endif 						// PulseShp.h

