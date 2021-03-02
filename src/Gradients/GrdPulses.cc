/* GrdPulses.cc ************************************************-*-c++-***
**									**
**	                         G A M M A 				**
**						 			**
**	Field Gradients Pulse Module 	           Implementation	**
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
**  gradients in NMR spectroscopy.                                      **
**                                                                      **
*************************************************************************/

#ifndef   GrdPulse_cc_			// Is file already included?
#  define GrdPulse_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Gradients/GrdPulses.h>	// Include the interface
#include <Gradients/sys_gradz.h>	// Include z-gradient spin systems
#include <HSLib/HSLibIF.h>		// Include Hilbert space stuff
#include <Pulses/PulGauss.h>		// Include Gaussian pulses
#include <Pulses/PulSinc.h>		// Include Sinc pulses
#include <vector>			// Include libstdc++ STL vectors

using std::vector;			// Using libstdc++ STL vectors
using std::string;			// Using libstdc++ strings

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

vector<gen_op> Ixpulse_Us(const sys_gradz& sys, const string& I, double angle)
  {
  unsigned N = sys.NSS();			// # of sub-systems
  vector<gen_op> Us;				// Array of propagators
  gen_op UIx = Ixpuls_U(sys, I, angle);		// Ideal pulse propagator
  for(unsigned i=0; i<N; i++)			// Loop over sub-systesm
    Us.push_back(UIx);				// Set same prop. each system
  return Us;					// Return the propagators
  }

vector<gen_op> Iypulse_Us(const sys_gradz& sys, const string& I, double angle)
  {
  unsigned N = sys.NSS();			// # of sub-systems
  vector<gen_op> Us;				// Array of propagators
  gen_op UIy = Iypuls_U(sys, I, angle);		// Ideal pulse propagator
  for(unsigned i=0; i<N; i++)			// Loop over sub-systesm
    Us.push_back(UIy);				// Set same prop. each system
  return Us;					// Return the propagators
  }
 
// ____________________________________________________________________________
// B                    Rectangular Pulse Functions
// ____________________________________________________________________________

vector<gen_op> Sxpuls_U(const sys_gradz& sys, vector<gen_op>& Hs,
                     const string& I, double offset, double tp, double angle)
  {
  unsigned N = Hs.size();
  vector<gen_op> Us;
  for(unsigned i=0; i<N; i++)
    Us.push_back(Sxpuls_U(sys,Hs[i],I,offset,tp,angle));
  return Us;
  }

vector<gen_op> Sypuls_U(const sys_gradz& sys, vector<gen_op>& Hs,
                     const string& I, double offset, double tp, double angle)
  {
  unsigned N = Hs.size();
  vector<gen_op> Us;
  for(unsigned i=0; i<N; i++)
    Us.push_back(Sypuls_U(sys,Hs[i],I,offset,tp,angle));
  return Us;
  }
 
// ____________________________________________________________________________
// C                    Gaussian Pulse Functions
// ____________________________________________________________________________

vector<gen_op> Gxpulse_U(const sys_gradz& sys, vector<gen_op>& Hs,
                const string& I, double tp, double angle, int N, double cutoff)
  {
  unsigned nss = Hs.size();			// Number of subsystems
  vector<gen_op> Us;				// Array of pulse props
  gen_op FX = Fx(sys, I);			// Get Fx for spins of type I
  gen_op U;					// Single propagator
  for(unsigned i=0; i<nss; i++)			// Loop over subsystems
    {
    U = Gpulse_U(Hs[i],FX,N,angle,tp,cutoff);	//  Prop for subsystem i
    Us.push_back(U);				//  Store in the array
    }
  return Us;
  }

vector<gen_op> Gypulse_U(const sys_gradz& sys, vector<gen_op>& Hs,
                const string& I, double tp, double angle, int N, double cutoff)
  {
  unsigned nss = Hs.size();			// Number of subsystems
  vector<gen_op> Us;				// Array of pulse props
  gen_op FY = Fy(sys, I);			// Get Fy for spins of type I
  gen_op U;					// Single propagator
  for(unsigned i=0; i<nss; i++)			// Loop over subsystems
    {
    U = Gpulse_U(Hs[i],FY,N,angle,tp,cutoff);	//  Prop for subsystem i
    Us.push_back(U);				//  Store in the array
    }
  return Us;
  }

// ____________________________________________________________________________
// D                         Sinc Pulse Functions
// ____________________________________________________________________________

vector<gen_op> SincPulseXUs(const sys_gradz& sys, vector<gen_op>& Hs,
                const string& I, double tp, double angle, int N, int NN)
  {
  unsigned nss = Hs.size();			// Number of subsystems
  vector<gen_op> Us;				// Array of pulse props
  gen_op FX = Fx(sys, I);			// Get Fx for spins of type I
  gen_op U;					// Single propagator
  for(unsigned i=0; i<nss; i++)			// Loop over subsystems
    {
    U = SincPulseU(Hs[i],FX,N,angle,tp,NN);
    Us.push_back(U);				//  Store in the array
    }
  return Us;
  }

vector<gen_op> SincPulseYUs(const sys_gradz& sys, vector<gen_op>& Hs,
                const string& I, double tp, double angle, int N, int NN)
  {
  unsigned nss = Hs.size();			// Number of subsystems
  vector<gen_op> Us;				// Array of pulse props
  gen_op FY = Fy(sys, I);			// Get Fy for spins of type I
  gen_op U;					// Single propagator
  for(unsigned i=0; i<nss; i++)			// Loop over subsystems
    {
    U = SincPulseU(Hs[i],FY,N,angle,tp,NN);
    Us.push_back(U);				//  Store in the array
    }
  return Us;
  }


#endif 							// GrdPulses.cc

