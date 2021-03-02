/* MultiIPul.h **************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Multiple System Ideal Pulses 			Interface 	**
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

#ifndef   GMultiIPul_h			// Is the file already included?
#  define GMultiIPul_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <MultiSys/MultiSys.h>		// Include multi_sys spin systems
#include <HSLib/GenOp.h>		// Include HS general operators
#include <string>			// Inlcude libstdc++ strings


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

MSVCDLL gen_op Iypuls(const multi_sys& msys, gen_op& sigma, double beta, int icomp=-1);
MSVCDLL gen_op Iypuls(const multi_sys& msys, gen_op& sigma, 
                                         int nspin, double beta, int icomp=-1);
MSVCDLL gen_op Iypuls(const multi_sys& msys, gen_op& sigma, 
                            const std::string& iso, double beta, int icomp=-1);


// ____________________________________________________________________________
//              PRODUCTION OF MULTISPIN SYSTEM PULSE PROPAGATORS
// ____________________________________________________________________________
 

MSVCDLL gen_op Ixpuls_U(const multi_sys& mys, int spin, double beta, int icomp=-1);
 
        // Input                msys    : A multi_spin spin system
        //                      spin    : A spin label
        //                      beta    : A pulse angle (degrees)
        //                      comp    : Spin system affected (-1 = all)
        // Output               U       : A propagator for an ideal pulse on
        //                                on the x-axis of angle beta to
        //                                spin spin on component(s) icomp
        // Note                         : The multize function used must
        //                                be defined for this to work
	// Uses multize(gen_op fct, int , double, multi_sys, int)
 

MSVCDLL gen_op Ixpuls_U(const multi_sys& mys, double beta, int icomp=-1);
 
        // Input                msys    : A multi_spin spin system
        //                      beta    : A pulse angle (degrees)
        //                      comp    : Spin system affected (-1 = all)
        // Output               U       : A propagator for an ideal pulse on
        //                                on the x-axis of angle beta to
        //                                all spins on component(s) icomp
        // Note                         : The multize function used must
        //                                be defined for this to work
	// Uses multize(gen_op fct, double, multi_sys, int)
 

MSVCDLL gen_op Iypuls_U(const multi_sys& mys, int spin, double beta, int icomp=-1);
 
        // Input                msys    : A multi_spin spin system
        //                      spin    : A spin label
        //                      beta    : A pulse angle (degrees)
        //                      comp    : Spin system affected (-1 = all)
        // Output               U       : A propagator for an ideal pulse on
        //                                on the y-axis of angle beta to
        //                                spin spin on component(s) icomp
        // Note                         : The multize function used must
        //                                be defined for this to work
	// Uses multize(gen_op fct, int , double, multi_sys, int)
 

MSVCDLL gen_op Iypuls_U(const multi_sys& mys, const std::string& iso, double beta, int icomp=-1);
 
        // Input                msys    : A multi_spin spin system
        //                      iso     : An isotope channel
        //                      beta    : A pulse angle (degrees)
        //                      comp    : Spin system affected (-1 = all)
        // Output               U       : A propagator for an ideal pulse on
        //                                on the y-axis of angle beta to
        //                                spins of type iso in the
	//				  component(s) icomp
        // Note                         : The multize function used must
        //                                be defined for this to work
	// Uses multize(gen_op fct, string, double, multi_sys, int)
 

MSVCDLL gen_op Iypuls_U(const multi_sys& mys, double beta, int icomp=-1);
 
        // Input                msys    : A multi_spin spin system
        //                      beta    : A pulse angle (degrees)
        //                      comp    : Spin system affected (-1 = all)
        // Output               U       : A propagator for an ideal pulse on
        //                                on the y-axis of angle beta to
        //                                spins in the component(s) icomp
        // Note                         : The multize function used must
        //                                be defined for this to work
	// Uses multize(gen_op fct, double, multi_sys, int)

 
MSVCDLL gen_op Ixypuls_U(const multi_sys& msys, double phi, double beta, int icomp=-1);
 
        // Input                msys    : A multi_spin spin system
        //                      phi     : A pulse phase (degrees)
        //                      beta    : A pulse angle (degrees)
        //                      comp    : Spin system affected (-1 = all)
        // Output               U       : A propagator for an ideal pulse of
        //                                phase phi and angle beta to
        //                                all spins in component(s) icomp
 
#endif							// MultiIPul.h
