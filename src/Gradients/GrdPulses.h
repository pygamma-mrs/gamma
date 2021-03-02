/* GrdPulses.h *************************************************-*-c++-***
**									**
**	                         G A M M A 				**
**						 			**
**	Field Gradients Pulse Module  	         	Interface	**
**						 			**
**      Copyright (c) 2001                                              **
**      Dr. Scott A. Smith             					**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**							 		**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  This GAMMA module provides functions for simulation of field        **
**  gradients in NMR spectroscopy. In particular, these functions 	**
**  apply pulses to a system in a field gradiend or return pulse	**
**  propagators that are useful for applying pulses to systems in	**
**  field gradients.							**
**                                                                      **
*************************************************************************/

#ifndef   GrdPul_h_			// Is file already included?
#  define GrdPul_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Gradients/sys_gradz.h>	// Know z-gradient systems
#include <HSLib/GenOp.h>		// Know general operators
#include <vector>			// Know libstdc++ STL vectors


// ____________________________________________________________________________
// A                       Ideal Pulse Functions
// ____________________________________________________________________________

/* Note that ideal pulses, being infinitely narrow in time, are perfect pulses.
   That is, they produce the same flip angle and phase over the entire system,
   regardless of "pulse offset" (of which there is none). In gradient work this
   means that all systems in a gradient see the same pulse and hence may use
   the same pulse propagator. Some functions below do return an array of ideal
   pulse propagators for convenience, but they are all identical. Not to worry,
   GAMMA will know they are identical and not use any extra memory to store
   them.                                                                     */

MSVCDLL std::vector<gen_op> Ixpulse_Us(const sys_gradz& sys,
                                                const std::string& I, double angle);
MSVCDLL std::vector<gen_op> Iypulse_Us(const sys_gradz& sys, 
                                                const std::string& I, double angle);

// ____________________________________________________________________________
// B                       Rectangular Pulse Functions
// ____________________________________________________________________________

MSVCDLL std::vector<gen_op> Sxpuls_U(const sys_gradz& sys, std::vector<gen_op>& Hs,
                     const std::string& I, double offset, double tp, double angle);
MSVCDLL std::vector<gen_op> Sypuls_U(const sys_gradz& sys, std::vector<gen_op>& Hs,
                     const std::string& I, double offset, double tp, double angle);

// ____________________________________________________________________________
// C                         Gaussian Pulse Functions
// ____________________________________________________________________________

MSVCDLL std::vector<gen_op> Gxpulse_U(const sys_gradz& sys, std::vector<gen_op>& Hs,
               const std::string& I, double tp, double angle, int N, double cutoff);

MSVCDLL std::vector<gen_op> Gypulse_U(const sys_gradz& sys, std::vector<gen_op>& Hs,
               const std::string& I, double tp, double angle, int N, double cutoff);

// ____________________________________________________________________________
// D                         Sinc Pulse Functions
// ____________________________________________________________________________

MSVCDLL std::vector<gen_op> SincPulseXUs(const sys_gradz& sys, std::vector<gen_op>& Hs,
                      const std::string& I, double tp, double angle, int N, int NN);

MSVCDLL std::vector<gen_op> SincPulseYUs(const sys_gradz& sys, std::vector<gen_op>& Hs,
                      const std::string& I, double tp, double angle, int N, int NN);

#endif 							// GrdPulses.h

