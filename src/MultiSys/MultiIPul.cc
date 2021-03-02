/* MultiIPul.cc *************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Multiple System Ideal Pulses 			Implementation 	**
**                                                                      **
**      Copyright (c) 1995                                              **
**      Nikolai Skrynnikov                                              **
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header:
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** This file of functions support multi_sys, the GAMMA class for        **
** handling mulitple spin systems. The routines herein generally use	**
** such a spin system & build up common identity pulse functions.	**
** The functions will mirror those defined in the Hilbert space library **
** (see module HSLib) but will return operators that exist in the       **
** composite Hilbert spin space spanning the components in multi_sys.   **
** The Hilbert space is a direct product space of the systems involved. **
**                                                                      **
*************************************************************************/

#ifndef   MultiIPul_cc			// Is the file already included?
#  define MultiIPul_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <MultiSys/MultiIPul.h>		// Include out header file
#include <MultiSys/MultiLib.h>		// Include multize functions
#include <MultiSys/MultiSys.h>		// Include multi_sys spin systems
#include <HSLib/GenOp.h>		// Inlcude H.S. general operators
#include <HSLib/PulseI.h>		// Inlcude H.S. ideal pulses
#include <HSLib/SpinSys.h>		// Include base spin ssytems

using std::string;			// Using libstdc++ strings

// ____________________________________________________________________________
// A              Pulses Acting Directly on Density Operators
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                           Pulse Along The Y Axis
// ----------------------------------------------------------------------------
/*
    Multize Arguments          Function Form                  Example
------------------------- ------------------------- ---------------------------
OpFct,Op,dble,msys,int    OpFct(sys,Op,double)        Iypuls(sys,sigma,90.0)
OpFct,Op,int,dbl,msys,int OpFct(sys,Op,int,double)    Iypuls(sys,sigma,1,90.0)
OpFct,Op,str,dbl,msys,int OpFct(sys,Op,string,double) Iypuls(sys,sigma,"1H",90.)

           Input                msys    : A multi_spin spin system
                                sigma   : A density operator for msys
                                beta    : A pulse angle (degrees)
                                nspin   : A spin index
                                iso     : An isotope channel
                                icomp   : Spin system affected (-1 = all)
           Output               sigmap  : A density operator which is
                                          sigma following the application of
                                          an ideal pulse on the y-axis of
                                          angle beta on in component(s) icomp
                                          of specified selectivity
           Note                         : The multize function used must
                                          be defined for this to work        */

gen_op Iypuls(const multi_sys& msys, gen_op& sigma,     double beta, int icomp)
  { return (multize(Iypuls, sigma, beta, msys, icomp)); }

gen_op Iypuls(const multi_sys& msys, gen_op& sigma,
                                             int nspin, double beta, int icomp)
  { return (multize(Iypuls, sigma, nspin, beta, msys, icomp)); }

gen_op Iypuls(const multi_sys& msys, gen_op& sigma, 
                                     const string& iso, double beta, int icomp)
  { return (multize(Iypuls, sigma, iso, beta, msys, icomp)); }







// ____________________________________________________________________________
// A              Pulses Acting Directly on Density Operators
// ____________________________________________________________________________

/* These use the "multize" function that are defined in section A of this file:

     gen_op multize(gen_op Op(args.....), const multi_sys& msys, int comp)

where arg[0] is the name of the (Pulse)Function and the reset of the arguments
are those sent to the function being called.

 Op(args....) = PulseFunction(const spin_sys&,SELECT,double angle,double phase)

           Input                msys	: A multi_spin spin system
                                sigma	: Density operator for the system
                                SEL	: Pulse selectivity
                                beta	: Pulse rotation angle (degrees)
                                phi	: Pulse phase angle (degrees)
				comp	: Spin system affected (-1 = all)
           Output               sigma	: Density operator after application
                                          of ideal pulse on selected spin(s)

          OverLoads of Function where Type of Second Argument Changes

           OL Function   SEL    Label   Selectivity Result
           -- --------  -----   -----   ------------------------------
           A. Iupuls     int    spin    Spin index, only spin affected
           B. Iupuls    string  iso     Isotope label, type iso spins affected
           C. Iupuls     ---    ---     All spins in sys affected
           D. Iupuls     int*   flags   Spins with flags TRUE affected
           E. Iupuls_sp  ---            Spins with sys flags TRUE affected   */

// ----------------------------------------------------------------------------
//                            Pulse Along The X Axis 
// ----------------------------------------------------------------------------

gen_op Ixpuls(const multi_sys& msys, gen_op& sigma, const string& iso, double beta, int icomp)
  { return (multize(Ixpuls, sigma, iso, beta, msys, icomp)); }


// ************* Pulse In the XY Plane On Axis Phi Degrees From X *************

// ----------------------------------------------------------------------------
//                           Ideal Pulse Propagators
// ----------------------------------------------------------------------------

gen_op Ixpuls_U(const multi_sys& msys, int spin, double beta, int icomp)

	// Input		msys    : A multi_spin spin system
        //                      spin	: A spin label
        //                      beta	: A pulse angle (degrees)
	//			comp    : Spin system affected (-1 = all)
        // Output               U	: A propagator for an ideal pulse on
	//				  on the x-axis of angle beta to 
	//				  spin spin in component(s) icomp
	// Note				: The multize function below must
	//				  be defined for this to work

  { return (multize(Ixpuls_U, spin, beta, msys, icomp)); }


gen_op Ixpuls_U(const multi_sys& msys, double beta, int icomp)

	// Input		msys    : A multi_spin spin system
        //                      beta	: A pulse angle (degrees)
	//			comp    : Spin system affected (-1 = all)
        // Output               U	: A propagator for an ideal pulse on
	//				  on the x-axis of angle beta to 
	//				  all spins in component(s) icomp
	// Note				: The multize function below must
	//				  be defined for this to work

  { return (multize(Ixpuls_U, beta, msys, icomp)); }


gen_op Iypuls_U(const multi_sys& msys, int spin, double beta, int icomp)

	// Input		msys    : A multi_spin spin system
        //                      spin	: A spin label
        //                      beta	: A pulse angle (degrees)
	//			comp    : Spin system affected (-1 = all)
        // Output               U	: A propagator for an ideal pulse on
	//				  on the y-axis of angle beta to 
	//				  spin spin in component(s) icomp
	// Note				: The multize function below must
	//				  be defined for this to work

  { return (multize(Iypuls_U, spin, beta, msys, icomp)); }


gen_op Iypuls_U(const multi_sys& msys, const string& iso, double beta, int icomp)

	// Input		msys    : A multi_spin spin system
	//			iso	: An isotope channel
        //                      beta	: A pulse angle (degrees)
	//			comp    : Spin system affected (-1 = all)
        // Output               U	: A propagator for an ideal pulse on
	//				  on the y-axis of angle beta to 
	//				  spins of type iso in component(s)
 	//				  icomp
	// Note				: The multize function below must
	//				  be defined for this to work

  { return (multize(Iypuls_U, iso, beta, msys, icomp)); }


gen_op Iypuls_U(const multi_sys& msys, double beta, int icomp)

	// Input		msys    : A multi_spin spin system
        //                      beta	: A pulse angle (degrees)
	//			comp    : Spin system affected (-1 = all)
        // Output               U	: A propagator for an ideal pulse on
	//				  on the y-axis of angle beta to 
	//				  all spins in component(s) icomp
	// Note				: The multize function below must
	//				  be defined for this to work

  { return (multize(Iypuls_U, beta, msys, icomp)); }


gen_op Ixypuls_U(const multi_sys& msys, double phi, double beta, int icomp)

	// Input		msys    : A multi_spin spin system
        //                      phi	: A pulse phase (degrees)
        //                      beta	: A pulse angle (degrees)
	//			comp    : Spin system affected (-1 = all)
        // Output               U	: A propagator for an ideal pulse on
	//				  on the y-axis of angle beta to 
	//				  all spins in component(s) icomp

  { return (multize(Ixypuls_U, phi, beta, msys, icomp)); }

#endif							// MultiIPul.cc
