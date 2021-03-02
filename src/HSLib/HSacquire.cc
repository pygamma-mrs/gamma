/* HSacquire.cc *************************************************-*-c++-*-
**									**
**	                      G A M M A 				**
**								 	**
**	Hilbert Space Acquisitions		   Implementation	**
**								 	**
**	Copyright (c) 1990, 1991, 1992			 		**
**	Tilo Levante, Scott Smith			 		**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fuer physikalische Chemie			 		**
**	8092 Zurich / Switzerland		 			**
**						 			**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** The GAMMA Platform Provides Functions for Simulation	of Magnetic	**
** Resonance Experiments and Other Associated 	Mathematical		**
** Capabilities.  The Set of Functions Herein Allow for the Simulation	**
** of Free Induction Decays(FID's) in Spin Hilbert Space & Adds 	**
** Abilities Related To Simulated Data Acquistion.			**
**								 	**
*************************************************************************/

#ifndef   HSacquire_cc_			// Is file already included?
#  define HSacquire_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <HSLib/HSacquire.h>		// Include the file header
#include <Matrix/row_vector.h>		// Knowledge of row vectors
#include <HSLib/GenOp.h>		// Knowledge of operators
//#include <HSLib/HSprop.h>		// Knowledge of nmr evolutions
#include <Basics/Gconstants.h>		// Knowledge of PI


// ____________________________________________________________________________
// A		         Generic Free Induction Decays
// ____________________________________________________________________________

/* These functions fill up a row vector with a Free Induction Decay (FID) based
   upon an initial density operator representing the initial state of the spin
   system evolving under a static Hamiltonian. The user specifies the property
   detected (usually transverse magnetization) as well as the time between FID
   points (dwell time). Note that the acquire1D class will also generate FID's
   (see the Level 1 Module in GAMMA) but is much more flexible. Also note that
   if one intends to apply a Fourier transform to the generated FID it should
   have a point size that is base 2 (e.g. 2048, 4096, .....)

	Input		sig0  : Operator propagated (initial dens. op.)
			D     : Detection operator in trace computation
	 		H     : Hamiltonian for propagation (in Hertz)
			U     : Propagator for one increment
	 		td    : Dwell time (seconds, per point)
	 		N     : Number of FID points generated
			fid   :	Data vector containing the result
	Output		      : None, FID data vector filled
	Note		      : Assumed fid is at least as large as N
	Note		      : Propagator is unitless, so 2*PI factor
				used to get radians in exponent 
        Note		      : If no block size is specified it will
                                be taken as the vector size
        Note		      : If no vector is specified on will be
                                generated and returned                      */  


void acquire(gen_op& S, gen_op& D, gen_op& H, double T, int N, row_vector& F, double CO)
  { FID(S,D,H,T,N,F,CO); }

void acquire(gen_op& sig0 , gen_op& D, gen_op& U, int N, row_vector& fid, double CO)
  { FID(sig0, D, U, 0, N, fid, CO); }

// This function was declared, but not defined, before 2011.05.16 (DCT)
void acquire(gen_op& sig0, gen_op& D, HSprop& Hp, int N, row_vector& fid, double CO)
  { FID(sig0, D, Hp, N, fid, CO); }


void FID(gen_op& sig0, gen_op& D, gen_op& H, double td, int N, row_vector& fid, double CO)
  {
  if(fid.size() < N) fid=row_vector(N);	// Insure block is large enough
  H.set_EBR();				// First put H into its eigenbase
  sig0.Op_base(H);			// Next put sig0 into eigenbase of H
  D.Op_base(H);				// Put detection op. to eigenbase of H
  complex z(0,-PIx2*td);		// Exponential factor for H -> U
  gen_op U;				// Propagator for evolution
  if(td) U = (z*H).exp();		// Set the td evolution propagator
  else   U = H;				// If td = 0, assume U was input
  int hs = dim(H);			// Hilbert space dimension
  int ls = hs*hs;			// Liouville space dimension
  complex *A = new complex[ls];		// Array for A(p)
  complex *B = new complex[ls];		// Array for B(p)

  int i,j,pos;
  for(pos=0,i=0; i<hs; i++)		// Generate the A & B arrays
    for(j=0; j<hs; j++)			// We'll only count non-zero
      {					// contributors to acquisition
      A[pos] = D(i,j)*sig0(j,i);
      B[pos] = conj(U(i,i), U(j,j));
      if(square_norm(A[pos])>CO) 	// Use elements if trace contribution
	pos ++;
      }

  for(int k=0; k<N; k++)		// Loop over desired points (& times)
    {
    z = 0;				//   Begin with zero point
    for(int p=0; p<pos; p++)		//   Loop over point contributors
      {	
      z += A[p];			//      Add contribution from p
      A[p] *= B[p];			//	                          k+1
      }					//      Adjust A so its really A*B
    fid.put(z,k);			//   Set point, move to next k (time) 
    } 
  delete [] A;				// Delete complex array A
  delete [] B;				// Delete complex array B
  }

void FID(gen_op& sig0 , gen_op& D, gen_op& U, int N, row_vector& fid, double CO)
  { FID(sig0, D, U, 0, N, fid, CO); }

void FID(gen_op& sig0 , gen_op& D, HSprop& U, int N, row_vector& fid, double CO)
  { gen_op UOp = U.Op(); FID(sig0, D, UOp, 0, N, fid, CO); }


row_vector acquire(gen_op& sig0 , gen_op& D, gen_op& H, double td, int N, double CO)
  { row_vector data(N); FID(sig0,D,H,td,N,data,CO); return data; }

row_vector acquire(gen_op& sig0 , gen_op& D, gen_op& U, int N, double CO)
  { row_vector data(N); FID(sig0,D,U,0,N,data,CO); return data; }

// This function was declared, but not defined, before 2011.05.16 (DCT)
row_vector acquire(gen_op& sig0, gen_op& D, HSprop& Hp, int N, double CO)
  { row_vector data(N); FID(sig0,D,Hp,N,data,CO); return data; }


row_vector FID(gen_op& sig0 , gen_op& D, gen_op& H, double td, int N, double CO)
  { row_vector data(N); FID(sig0,D,H,td,N,data,CO); return data; }

row_vector FID(gen_op& sig0 , gen_op& D, gen_op& U, int N, double CO)
  { row_vector data(N); FID(sig0,D,U,0,N,data,CO); return data; }

row_vector FID(gen_op& sig0 , gen_op& D, HSprop& U, int N, double CO)
  { row_vector data(N); gen_op UOp=U.Op(); FID(sig0,D,UOp,0,N,data,CO); return data; }

// ____________________________________________________________________________
// B		   Acquisiton Generation Output Functions
// ____________________________________________________________________________

/* This functions exist primarly to check the acquisition functions found
   earlier in this file. They will simply print output aspects of the
   calculation used for the specified input parameters.                      */
 
// ____________________________________________________________________________
// Z		      Free Induction Decay Applied Theory
// ____________________________________________________________________________

/* We will work in the eigenbasis of the propergator U, where U is generated
   from the static Hamiltonian H as

                                      -2*pi*i*H*td
                                 U = e

   Both U and H are represented by diagonal matrices in this eigenbasis and
   since H is Hermtitian, U will be unitary (adjoint == inverse).  An 
   acquisition point k, which corresponds to time t = k*td, is given by

   fid(k) = Trace { D * sigma(t) }

                         k                        k
          = Trace { D * U (td)  * sig0 * adjoint[U (td)] }

                                 k                            k
          = Sum i,j { <i|D|j><j|U (td)|j>*<j|sig0|i>*conj(<i|U (td)|i>) }

                                                                 k
          = Sum i,j { <i|D|j><j|sig0|i> * [<j|U|j>*conj(<i|U|i>)]  }

                                     k                        k
          = Sum i,j { A(i,j) * B(i,j) }    =  Sum p { A(p) * B (p) }

   The FID functions all implement the last form for generating a series
   of expectation values evenly incremented in time.  The routines will
   immediately generate the arrays A and B then quicky do the sum. The
   points are repeatly generated by taking powers of B elements.             */


#endif 							// HSacquire.cc
