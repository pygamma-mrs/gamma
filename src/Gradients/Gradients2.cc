/* Gradients2.cc ***********************************************-*-c++-***
**									**
**	                         G A M M A 				**
**						 			**
**	Field Gradients Module 2 	           Implementation	**
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
**  This GAMMA module provides functions for simulation of field        **
**  gradients in NMR spectroscopy.                                      **
**                                                                      **
*************************************************************************/

#ifndef   Gradients2_cc_		// Is file already included?
#  define Gradients2_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Gradients/Gradients2.h>	// Include the interface
#include <Gradients/sys_gradz.h>	// Include z-gradient spin systems
#include <HSLib/HSLibIF.h>		// Include Hilbert space stuff
#include <vector>			// Include libstdc++ STL vectors
 
// ____________________________________________________________________________
// A                           Hamiltonian Functions
// ____________________________________________________________________________

/* These routines assumes that a system in a z-gradient has an associated 
   Hamiltonian H0 without the gradient applied.  When the gradient is applied
   this Hamiltonian is modified by an Fz term depending upon the gradient
   strength at a particular position.  The returned array contains modified
   Hamiltonians which are H0's adjusted for each sub-system in the gradients
   by the Bo field gradient.                                                 */

	// Input		sys	: System in a Bo z-Gradient
	//			Ho	: Active Hamiltonian, no Gradient
	//			H	: Array of Hamiltonians
	// Output		void	: The array of Hamiltonains H
	//				  will be altered according to
	//				   H(i) = H0 + HZ(gradz,i) 

void Hzgrad(const sys_gradz& sys, gen_op& H0, gen_op* H)
  {
  int nss = sys.NSS();			// Number of subsystems
  int ns = sys.spins();			// Number of spins
  int i;
  gen_op Hcs0; 				// Raw shift Hamiltonian
  for(i=0; i<ns; i++)			// Fill raw shift Hamiltonian
    Hcs0 -= sys.gamma(i)*Iz(sys,i);	// (units are in rad/T)
  H0.set_EBR();
  Hcs0.Op_base(H0);
  double dBodm = sys.BoGrad()*RAD2HZ;	// Get applied gradient (T-Hz/rad-m)
  double Bi;				// Shift scaling from gradient
  for(i=0; i<nss; i++)			// Loop the subsystems and fill
    {					// H vector with H0 + Hi
    Bi = dBodm*sys.SysDist(i);		// 	Grad. field subsystem i
    H[i] = H0 + Bi*Hcs0;		//	Add in additional shift
    }
  return;
  }

std::vector<gen_op> Hzgrad(const sys_gradz& sys, gen_op& H0)
  {
  std::vector<gen_op> Hs;		// Array of Hamiltonians
  int i, nss = sys.NSS();		// Number of subsystems
  int ns = sys.spins();			// Number of spins
  gen_op Hcs0; 				// Raw shift Hamiltonian
  for(i=0; i<ns; i++)			// Fill raw shift Hamiltonian
    Hcs0 -= sys.gamma(i)*Iz(sys,i);	// (units are in rad/T)
  H0.set_DBR();				// Put EBR in default basis
  Hcs0.Op_base(H0);			// Set shift in same basis
  double dBodm = sys.BoGrad()*RAD2HZ;	// Get applied gradient (T-Hz/rad-m)
  double Bi;				// Shift scaling from gradient
  for(i=0; i<nss; i++)			// Loop the subsystems and fill
    {					// H vector with H0 + Hi
    Bi = dBodm*sys.SysDist(i);		// 	Grad. field subsystem i
    Hs.push_back(H0 + Bi*Hcs0);		//	Add in additional shift
    }
  return Hs;
  }

// ____________________________________________________________________________
// B                           Propagator Functions
// ____________________________________________________________________________

/* Given an array of Hamiltonians and evolution times, this function will
   produce an array of propagators.                                          */

void Props(int NSS, gen_op* Hs, double t, gen_op* Us)
  { for(int i=0; i<NSS; i++) Us[i] = prop(Hs[i], t); }

std::vector<gen_op> Props(std::vector<gen_op>& Hs, double t)
  {
  unsigned ns = Hs.size();
  std::vector<gen_op> Us;
  for(unsigned i=0; i<ns; i++)
    Us.push_back(prop(Hs[i], t));
  return Us;
  }
    
#endif 							// Gradients2.cc

