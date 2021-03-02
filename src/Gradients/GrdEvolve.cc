/* GrdEvolve.cc ************************************************-*-c++-***
**									**
**	                         G A M M A 				**
**						 			**
**	Field Gradients Time Evolution 	           Implementation	**
**						 			**
**      Copyright (c) 1996                                              **
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
** This GAMMA module provides functions for simulation of field		**
** gradients in NMR spectroscopy.					**
**                                                                      **
*************************************************************************/

#ifndef   GrdEvolve_cc_			// Is file already included?
#  define GrdEvolve_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Gradients/GrdEvolve.h>	// Include the interface
#include <Gradients/sys_gradz.h>	// Include z-gradient spin systems
#include <HSLib/HSLibIF.h>		// Include Hilbert space stuff
#include <vector>			// Include libstdc++ STL vectors
#include <Level2/RelaxBas.h>		// Include phenomenological relaxation

using std::vector;			// Using libstdc++ STL vectors


// ____________________________________________________________________________
// A            Evolve A Single State Into An Array of States
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//                     Hilbert Space Evolution Functions
//-----------------------------------------------------------------------------

/* These functions assume that we have an single density operator that is
   appropriate for all sub-systems in the gradient. However, the each will
   evolve in time under a different Hamiltonian (or equivalently are affected
   by a different propagator.) So the functions return an array of evolved
   density operators, one for each of the subsystems.  The arguments demand
   an initial density operator along with as set of evolution Hamiltonians
   with an evolution time or (equivalently) a set of evolution propagators. 
   They then produce a set of evolved density operators.

                                  t
      sigma (t) = U (t)*sigma(0)*U (t) = exp(-i*H *t) * sigma * exp(i*H *t)
           i       i              i              i                     i     */

vector<gen_op> evolve(gen_op& sigma0, vector<gen_op>& Hs, double t)
  {
  unsigned nss = Hs.size();
  vector<gen_op> sigmas;
  for(unsigned i=0; i<nss; i++)
    sigmas.push_back(evolve(sigma0, Hs[i], t));
  return sigmas;
  }

vector<gen_op> evolve(gen_op& sigma0, vector<gen_op>& Us)
  {
  unsigned nss = Us.size();
  vector<gen_op> sigmas;
  for(unsigned i=0; i<nss; i++)
   sigmas.push_back(evolve(sigma0, Us[i]));
  return sigmas;
  }

//-----------------------------------------------------------------------------
//                Evolution Under Phenomenological Relaxation
//-----------------------------------------------------------------------------

vector<gen_op> evolve(gen_op& sigma0, vector<gen_op>& Hs, RBasic& R, double t)
  {
  unsigned nss = Hs.size();			// Number of Hamiltonians
  vector<gen_op> sigmas;			// Array of evolved dens. ops.
  for(unsigned i=0; i<nss; i++)			// Loop over each Hamiltonian
    {
    R.SetH0(Hs[i]);				//   Set active Hamiltonian
    sigmas.push_back(R.Evolve(sigma0, t));	//   Evolve under R & Hs[i]
    }						//   Store result as simgas[i]
  return sigmas;
  }

//-----------------------------------------------------------------------------
/* Evolve An Array Of States Under A Single Propagator To An Array of States
//-----------------------------------------------------------------------------

   These functions assume that we have a set of density operators that 
   describe a spin system in a field gradient, each operator in the set
   corresponds to a single sub-system.  Each is assumed to evolve under the
   same Hamiltonian for the same amount of time (or equivalently are affected
   by the same propagator.) The functions then returns an array of evolved
   density operators, one for each of the subsystems.  The arguments demand
   an initial set of density operators along with an evolution Hamiltonian
   with an evolution time or (equivalently) an evolution propagator. 
   They then produce a set of evolved density operators.                     */

vector<gen_op> evolve(vector<gen_op>& sigmas0, gen_op& H, double t)
  {
  gen_op U = prop(H,t);
  return evolve(sigmas0, U);
  }

vector<gen_op> evolve(vector<gen_op>& sigmas0, gen_op& U)
  {
  unsigned nss = sigmas0.size();
  vector<gen_op> sigmas;
  for(unsigned i=0; i<nss; i++)
    sigmas.push_back(evolve(sigmas0[i], U));
  return sigmas;
  }

/* Evolve An Array Of States Under Array Of Propagators To New Array of States

   These functions assume that we have a set of density operators that 
   describe a spin system in a field gradient, each operator in the set
   corresponds to a single sub-system.  Each is assumed to evolve under its
   own particular Hamiltonian for the same amount of time (or equivalently is
   affected by its particular propagator.) The functions then returns an array
   of evolved density operators, one for each of the subsystems.  The input
   arguments demand an initial set of density operators along with a set of
   evolution Hamiltonians and an evolution time or (equivalently) an set of
   evolution propagators. They produce a set of evolved density operators.   */

// ----------------------------------------------------------------------------
//                               t
//  sigma (t) = U (t)*sigma (0)*U (t) = exp(-i*H *t) * sigma (0) * exp(i*H *t)
//       i       i         i     i              i           i             i
// ----------------------------------------------------------------------------

vector<gen_op> evolve(vector<gen_op>& sigma0, vector<gen_op>& H, double t)
  {
  unsigned N = sigma0.size();
  vector<gen_op> sigma;
  for(unsigned i=0; i<N; i++)
    sigma.push_back(evolve(sigma0[i], H[i], t));
  return sigma;
  }

vector<gen_op> evolve(vector<gen_op>& sigma0, vector<gen_op>& U)
  {
  unsigned N = sigma0.size();
  vector<gen_op> sigma;
  for(unsigned i=0; i<N; i++)
    sigma.push_back(evolve(sigma0[i], U[i]));
  return sigma;
  }

vector<gen_op> evolve(vector<gen_op>& sigs, vector<gen_op>& Hs, RBasic& R, double t)
  {
  unsigned nss = Hs.size();			// Get number of operators
  vector<gen_op> sigmas;			// For evovled density operators
  for(unsigned i=0; i<nss; i++)			// Loop over all operators
    {						// and evolve each independently
    R.SetH0(Hs[i]);				//   Set the Hamiltonian of i
    sigmas.push_back(R.Evolve(sigs[i], t));	//   Evolve density operator i
    }
  return sigmas;
  }

#endif 							// GradEvolve.cc

