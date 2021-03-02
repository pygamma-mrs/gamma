/* MultiExch.h ***************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Multiple System Library		               Interface	**
**                                                                      **
**      Copyright (c) 1995                                              **
**      Nikolai Skrynnikov                                              **
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** This module of function supports the multi_sys, the GAMMA class      **
** handling mulitple spin systems.  The routines herein generally       **
** involve such a spin system and build up common operators, in this    **
** case in a direct product space of the systems involved.              **
**                                                                      **
*************************************************************************/

#ifndef   MultiExch_h			// Is the file already included?
#  define MultiExch_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/GenOp.h>		// Knowledge of operators
#include <LSLib/SuperOp.h>		// Knowledge of superoperators
#include <Matrix/row_vector.h>		// Knowledge of row vectors
#include <MultiSys/MultiSys.h>
#include <HSLib/SpinSystem.h>
#include <LSLib/sys_dynamic.h>

MSVCDLL super_op Xnm( const    multi_sys& msys);
MSVCDLL matrix   Xnmp(const    multi_sys& msys, int p);
MSVCDLL void     Xnmpblk(const multi_sys& msys, const ExchProc& Pro, matrix& Xp,
            double K, int cmpI, int cmpJ, int Io, int Iend, int Jo, int Jend);
MSVCDLL bool     Xnmpelem(const multi_sys& msys, const ExchProc& Pro,
                const row_vector& braI, const row_vector& ketI,
                const row_vector& braJ, const row_vector& ketJ, 
                int cmpI, int cmpJ, int& hsnorm);

MSVCDLL super_op Xm(const multi_sys& msys);

// ----------------------------------------------------------------------------
//  This Function Mimics Xnm Above But Just Outputs Details Of The Calculation
// ----------------------------------------------------------------------------

MSVCDLL void Xnm(std::ostream&     ostr, const multi_sys& sys);
MSVCDLL void Xnmp(std::ostream&    ostr, const multi_sys& msys, int p);
MSVCDLL void Xnmpdblk(std::ostream& ostr, const multi_sys& msys, double K, int Io, int Iend);
MSVCDLL void Xnmpblk(std::ostream& ostr,  const multi_sys& msys, const ExchProc& Pro,
             double K, int cmpI, int cmpJ, int Io, int Iend, int Jo, int Jend);

MSVCDLL super_op XXnm( const    multi_sys& msys);


#endif								// MultiExch.h
