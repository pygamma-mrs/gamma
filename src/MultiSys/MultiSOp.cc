/* MultiSOp.cc **************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Multiple System Spin Operators             Implementation 	**
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
** handling mulitple spin systems.  The routines herein generally       **
** involve such a spin system and build up common spin operators.       **
** The functions will mirror those defined in the Hilbert space library **
** (see module HSLib) but will return operators that exist in the       **
** composite Hilbert spin space spanning the components in multi_sys.   **
** The Hilbert space is a direct product space of the systems involved. **
**                                                                      **
*************************************************************************/

#ifndef _MultiSOp_cc			// Is the file already included?
#define _MultiSOp_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// Then this is the implementation
#endif

#include <HSLib/SpinSys.h>		// Include base spin ssytems
#include <MultiSys/MultiLib.h>		// Include the header file
#include <MultiSys/MultiSys.h>		// Include multi_sys spin systems
#include <HSLib/SpinOp.h>		// Include simple spin operators
#include <HSLib/SpinOpCmp.h>		// Include composite spin operators
#include <HSLib/SpinOpRot.h>		// Include rotation spin operators
#include <HSLib/GenOp.h>		// Inlcude general operators

using std::string;			// Using libstdc++ strings

// ____________________________________________________________________________
// D             PRODUCTION OF MULTISPIN SYSTEM SPIN OPERATORS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                       Total System Spin Operators
// ----------------------------------------------------------------------------

/* 	   Input                msys    : A multi_spin spin system
           Output               Fu      : Fu for the entire system
                                          where for u we have
                                             x -> x   y -> y  z -> z
                                             e -> e   m -> -  p -> + 
	   Note				: Uses multize function

     gen_op multize(spin_op op(const spin_sys&), const multi_sys &msys)      */

gen_op Fx(const multi_sys& msys) { return multize(Fx, msys); }
gen_op Fy(const multi_sys& msys) { return multize(Fy, msys); }
gen_op Fz(const multi_sys& msys) { return multize(Fz, msys); }
gen_op Fe(const multi_sys& msys) { return multize(Fe, msys); }
gen_op Fp(const multi_sys& msys) { return multize(Fp, msys); }
gen_op Fm(const multi_sys& msys) { return multize(Fm, msys); }


// ----------------------------------------------------------------------------
//                     Isotope Selective Spin Operators
// ----------------------------------------------------------------------------

/*         Input                msys	: A multi_spin spin system
                                I	: Isotope type
           Output               Fu	: F* selective for isotope
                                          where for u we have
                                             x -> x   y -> y  z -> z
                                             e -> e   m -> -  p -> + 
	   Note				: Uses multize function

 gen_op multize(spin_op op(const spin_sys&, const string&),const multi_sys&) */

gen_op Fx(const multi_sys &msys,const string& I) { return multize(Fx,I,msys); }
gen_op Fy(const multi_sys &msys,const string& I) { return multize(Fy,I,msys); }
gen_op Fz(const multi_sys &msys,const string& I) { return multize(Fz,I,msys); }
gen_op Fe(const multi_sys &msys,const string& I) { return multize(Fe,I,msys); }
gen_op Fm(const multi_sys &msys,const string& I) { return multize(Fm,I,msys); }
gen_op Fp(const multi_sys &msys,const string& I) { return multize(Fp,I,msys); }

// ____________________________________________________________________________
// F          PRODUCTION OF MULTISPIN SYSTEM ROTATION OPERATORS
// ____________________________________________________________________________

gen_op Rz(const multi_sys& msys, double beta, int icomp)

        // Input                msys    : A multi_spin spin system
        //                      beta    : A rotation angle (degrees)
        //                      comp    : Spin system affected (-1 = all)
        // Output               RZ      : A spin rotation operator (prop)
        //                                which rotates spin angular
        //                                momentum about the z-axis by
        //                                the angle beta for all spins
        //                                in component(s) icomp
        // Note                         : The multize function below must
        //                                be defined for this to work

  { return (multize(Rz, beta, msys, icomp)); }


#endif							// MultiSOp.cc
