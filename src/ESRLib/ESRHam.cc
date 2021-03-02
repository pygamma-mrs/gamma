/* ESRHam.cc ****************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**   ESR hamiltonian                                  Implementation	**
**                                                                      **
**  Copyright (c) 2000							**
**  Scott A. Smith							**
**  National High Magnetic Field Laboratory                           	**
**  1800 E. Paul Dirac Drive                                          	**
**  Box 4005                                                          	**
**  Tallahassee, FL 32306                                             	**
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/
     
/*************************************************************************
**									**
**  Description								**
**									**
**  Module ESR Hamiltonians contains functions that return commonly     **
**  used Hamiltonians in ESR/EPR/EMR computations.                      **
**																		**
*************************************************************************/

#ifndef   ESRHam_cc_			// Is this file already included?
#  define ESRHam_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <ESRLib/ESRHam.h>		// Include the interface
#include <Matrix/matrix.h>		// Include matrices
#include <HSLib/SpinOpSng.h>		// Inlcude spin spin operators
 
// ____________________________________________________________________________
// F              ISOTROPIC ELECTRON G FACTOR HAMILTONIANS
// ____________________________________________________________________________

/* This Hamiltonian has a negative sign for consistency with the nuclear spin
   chemical shift Hamiltonian, likely this is NOT the sign convention ESR 
   people like.  Can't please everybody though. All electrons are "shielded"
   from the free electron (negative gyromagnetic ration, positive G).

                                H  = -v * S
                                 g         z                                */

gen_op Hg(const CubicSys& CuSys, double v)
  {
  int hs    = CuSys.HS();			// Get Hilbert space
  matrix H = -v*Iz(hs);				// Empty spin operator
  return gen_op(H);
  }

// ____________________________________________________________________________
// A               Magnetic Ion Constructors, Destructor, Assignment
// ____________________________________________________________________________
 
gen_op O40(const CubicSys& CuSys)
  {
  double J = CuSys.qn();
  int hs = CuSys.HS();
  matrix Jz = Iz(hs);
  matrix JzSq = Jz*Jz;
  matrix Jz4th = JzSq*JzSq;
  matrix I = Ie(hs);
  matrix Hmx = 35*Jz4th
             - (30*J*(J+1)-25)*JzSq
             - 6*J*(J+1)*I
             + 3*J*J*(J+1)*(J+1)*I;
  return gen_op(Hmx);
  }

gen_op O44(const CubicSys& CuSys)
  {
  int hs = CuSys.HS();
  matrix Jp = Ip(hs);
  matrix Jm = Im(hs);
  Jp *= Jp;  Jp *= Jp;
  Jm *= Jm;  Jm *= Jm;
  return gen_op(0.5*(Jp+Jm));
  }

gen_op O60(const CubicSys& CuSys)
  {
  double J = CuSys.qn();
  int hs = CuSys.HS();
  matrix Jz = Iz(hs);
  matrix I = Ie(hs);
  matrix JzSq = Jz*Jz;
  matrix Jz4th = JzSq*JzSq;
  matrix Jz6th = Jz4th*JzSq;
  matrix Hmx = 231*Jz6th
             - (105*(3*J*(J+1)-7))*Jz4th
             + (105*J*J*(J+1)*(J+1)-525*J*(J+1)+294)*JzSq
             - 5*J*J*J*(J+1)*(J+1)*(J+1)*I
             + 40*J*J*(J+1)*(J+1)*I
             - 60*J*(J+1)*I;
  return gen_op(Hmx);
  }
  
gen_op O64(const CubicSys& CuSys)
  {
  double J = CuSys.qn();
  int hs = CuSys.HS();
  matrix Jz = Iz(hs);
  matrix Jp = Ip(hs);
  matrix Jm = Im(hs);
  matrix JzSq = Jz*Jz;
  matrix Jp4th = Jp*Jp; Jp4th *= Jp4th;
  matrix Jm4th = Jm*Jm; Jm4th *= Jm4th;
  matrix I = Ie(hs);
  matrix Hmx1 = 11*JzSq-(J*(J+1)+38)*I;
  matrix Hmx2 = Jp4th + Jm4th;
  return gen_op(0.25*(Hmx1*Hmx2 + Hmx2*Hmx1));
  }

#endif							// ESRHam.cc
