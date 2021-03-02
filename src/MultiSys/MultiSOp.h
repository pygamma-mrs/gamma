/* MultiSOp.h ***************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Multiple System Spin Operators			Interface	**
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
** involve such a spin system and build up common spin operators. 	**
** The functions will mirror those defined in the Hilbert space library	**
** (see module HSLib) but will return operators that exist in the 	**
** composite Hilbert spin space spanning the components in multi_sys.	**
** The Hilbert space is a direct product space of the systems involved. **
**                                                                      **
*************************************************************************/

#ifndef   MultiSOp_h			// Is the file already included?
#  define MultiSOp_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <MultiSys/MultiSys.h>		// Knowledge of multiple systems
#include <HSLib/GenOp.h>		// Knowledge of operators
//#include <HSLib/SpinSystem.h>
//#include <LSLib/sys_dynamic.h>

 
// ______________________________________________________________________
// A                 Total Spin System Spin Operatators
// ______________________________________________________________________
 
MSVCDLL gen_op Fx(const multi_sys &msys);
MSVCDLL gen_op Fy(const multi_sys &msys);
MSVCDLL gen_op Fz(const multi_sys &msys);
MSVCDLL gen_op Fe(const multi_sys &msys);
MSVCDLL gen_op Fp(const multi_sys &msys);
MSVCDLL gen_op Fm(const multi_sys &msys);
 
        // Input                msys    : A multi_spin spin system
        // Output               F*      : F* for the entire system
	//			          where for * we have
	//				     x -> x   y -> y  z -> z
	//				     e -> e   m -> -  p -> +
 
 
// ----------------------------------------------------------------------
//                   Isotope Selective Spin Operators
// ----------------------------------------------------------------------


MSVCDLL gen_op Fx(const multi_sys& msys, const std::string& iso); 
MSVCDLL gen_op Fy(const multi_sys& msys, const std::string& iso); 
MSVCDLL gen_op Fz(const multi_sys& msys, const std::string& iso); 
MSVCDLL gen_op Fe(const multi_sys& msys, const std::string& iso); 
MSVCDLL gen_op Fm(const multi_sys& msys, const std::string& iso); 
MSVCDLL gen_op Fp(const multi_sys& msys, const std::string& iso);

        // Input                msys    : A multi_spin spin system
        //                      iso     : Isotope type
        // Output               F*      : F* selective for isotope
	//			          where for * we have
	//				     x -> x   y -> y  z -> z
	//				     e -> e   m -> -  p -> +

// ____________________________________________________________________________
//           PRODUCTION OF MULTISPIN SYSTEM ROTATION OPERATORS
// ____________________________________________________________________________

MSVCDLL gen_op Rz(const multi_sys& msys, double beta, int icomp=-1);

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


#endif						// MultiSOp.h
