/* MultiWBR.h ***************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Multiple System Spin BWR Relaxation 		Interface	**
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
** involve such a spin system & build up BWR relaxation superoperators.	**
** The functions will mirror those defined in the Hilbert space library	**
** (see module HSLib) but will return operators that exist in the 	**
** composite Hilbert spin space spanning the components in multi_sys.	**
** The Hilbert space is a direct product space of the systems involved. **
**                                                                      **
*************************************************************************/

#ifndef   MultiWBR_h			// Is the file already included?
#  define MultiWBR_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <MultiSys/MultiSys.h>		// Knowledge of multiple systems
#include <LSLib/SuperOp.h>		// Knowledge of superoperators
#include <HSLib/GenOp.h>		// Knowledge of operators

// ____________________________________________________________________________
//          PRODUCTION OF MULTISPIN SYSTEM RELAXATION SUPEROPERATORS
// ____________________________________________________________________________

 
        // Input                msys    : A multi_spin spin system
        //                      H       : Isotropic Hamiltonian
        //                      type    : Computation type
        //                      level   : Computation level
        // Output               R       : Relaxation superoperator
        //                                for quadrupolar interactions
        // Note                         : The multize function below must
        //                                be defined for this to work
        // Output               R       : Relaxation superoperator
        //                                for dipole-dipolar interactions
        // Note                         : The multize function below must
        //                                be defined for this to work
        // Output               R       : Relaxation superoperator
        //                                for shielding anisotropy
        // Note                         : The multize function below must
        //                                be defined for this to work

MSVCDLL super_op RQQ(const multi_sys& msys, gen_op& H, int type=0, int level=4);
MSVCDLL super_op RCC(const multi_sys& msys, gen_op& H, int type=0, int level=4);
MSVCDLL super_op RDD(const multi_sys& msys, gen_op& H, int type=0, int level=4);
//super_op RDQ(const multi_sys& msys, gen_op& H, int type=0, int level=4);
//super_op RQD(const multi_sys& msys, gen_op& H, int type=0, int level=4);
//super_op RCD(const multi_sys& msys, gen_op& H, int type=0, int level=4);
//super_op RDC(const multi_sys& msys, gen_op& H, int type=0, int level=4);
MSVCDLL super_op RCQ(const multi_sys& msys, gen_op& H, int level=4);
MSVCDLL super_op RQC(const multi_sys& msys, gen_op& H, int level=4);

#endif						// MultiWBR.h
