/* MultiLOp.h ***************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Multiple System Super Operators			Interface	**
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
** involve such a spin system and build up common superoperators. 	**
** The functions mirror those defined in the Liouville space library	**
** (see module LSLib) but will return super operators that exist in the	**
** composite Liouville spin space spanning the components in multi_sys.	**
** The composite Liouville space in this case is a direct product of	**
** the Liouville spaces of the systems involved. 			**
**                                                                      **
*************************************************************************/

#ifndef   MultiLOp_h			// Is the file already included?
#  define MultiLOp_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <LSLib/SuperOp.h>		// Knowledge of superoperators
#include <MultiSys/MultiSys.h>		// Knowledge of multiple systems
#include <HSLib/GenOp.h>		// Knowledge of operators

// ____________________________________________________________________________
//                PRODUCTION OF MULTISPIN SYSTEM SUPEROPERATORS
// ____________________________________________________________________________
 

        // Input                msys    : A multi_spin spin system
        //                      Heff    : Effective Hamiltonian (Hz)
        // Output               LOp     : Hamiltonian commutation superop
        // Note                         : This superoperator is returned
        //                                in the default basis
        // Note                         : LOp returned in rad/sec
        // Uses multize(super_op SOpFct(gen_op&), gen_op& , multi_sys&)
        // Uses super_op SOpFct(gen_op&) ---> super_op Hsuper(Heff)


MSVCDLL super_op Hsuper(const multi_sys& msys, const gen_op& Heff);
MSVCDLL super_op Lo(const multi_sys &msys);

MSVCDLL super_op U_LS(gen_op& H);
MSVCDLL super_op Uinv_LS(gen_op& H);
MSVCDLL super_op Op_Ebase(super_op& L, gen_op& H);

#endif								// MultiLOp.h
