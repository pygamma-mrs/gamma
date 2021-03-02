/* Bloch.h ******************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**      Bloch Equation Entities                   Interface		**
**									**
**      Copyright (c) 1995						**
**      S.A. Smith 							**
**      National High Magnetic Field Laboratory				**
**      1800 E. Paul Dirac Drive					**
**      Tallahassee Florida, 32306-4005					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
**  Description 							**
**									**
**  This module contains functions which aid in simulations involving 	**
**  the phenomenological Bloch equations. The routines determine common	**
**  entities involved in the Bloch approach: Magnetizaiton Vectors,	**
**  R Matrices, K Matrices,... 						**
**									**
*************************************************************************/

#ifndef   Gbloch_h_ 			// Is file already included?
#  define Gbloch_h_ 1      		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface   
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/matrix.h>		// Know GAMMA matrices
#include <Basics/ParamSet.h>		// Know GAMMA parameter sets

// ____________________________________________________________________________
// A                  Interactive I/O Auxiliary Functions
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                      Output Magnitude Functions
// ----------------------------------------------------------------------------

MSVCDLL int DoubleMag(double x);

// ----------------------------------------------------------------------------
//                        Output Unit Functions
// ----------------------------------------------------------------------------

MSVCDLL std::string SecUnits(int mag, double& sf);
MSVCDLL std::string HzUnits(int  mag, double& sf);

// ----------------------------------------------------------------------------
//                 Initial Magnetization Vector Functions
// ----------------------------------------------------------------------------

MSVCDLL matrix Mo_vector(double Mox=0, double Moy=0, double Moz=1);
//matrix Mo_vector(int argc, char* argv[], matrix& Meq, int& qn);

// ----------------------------------------------------------------------------
//           Steady-State Magnetization Vector Functions
// ----------------------------------------------------------------------------

// This next function is not currently implemented, 2011.05.16 (DCT)
// MSVCDLL matrix Mss_vector(matrix& K, matrix& R, matrix& Meq);

// This next function was implemented, but until now not
// declared 2011.05.16 (DCT)
MSVCDLL matrix Mo_vector(int argc, char* argv[], matrix& Meq, int& qn);

// ----------------------------------------------------------------------------
//                      Analysis Functions
// ----------------------------------------------------------------------------

MSVCDLL void analyze(double tinc, int& ntimes,
          int& do_ss, int& qn, double T1, double gamB1, double w);


// ----------------------------------------------------------------------------
//                          Bloch Related Parameters
// ----------------------------------------------------------------------------

MSVCDLL void bloch_T1T2(const ParameterSet& pset, std::ostream& ofstr, double& T1, double& T2);

MSVCDLL void bloch_Mo(const ParameterSet& pset,   std::ostream& ofstr, double& Mx, double& My, double& Mz);

MSVCDLL void bloch_B1(const ParameterSet& pset,   std::ostream& ofstr, double& gamB1, double& phi);

MSVCDLL void bloch_Woff(const ParameterSet& pset, std::ostream& ofstr, double& Woff);


// ----------------------------------------------------------------------------
//                        Trajectory Timing Functions
// ----------------------------------------------------------------------------


MSVCDLL void TrajTiming(int argc, char* argv[],
             double& tinc, int& N, int& qn, double T1, double gamB1, double w);


#endif								// Bloch.h
