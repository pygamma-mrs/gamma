/* GrdAcquire.cc ***********************************************-*-c++-***
**									**
**	                         G A M M A 				**
**						 			**
**	Field Gradients Acquisition Module            Implementation	**
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
**  gradients in NMR spectroscopy. In particular, these functions       **
**  perform acquisitions of evenly spaced expectation values on a       **
**  a system in a field gradient.                                       **
**                                                                      **
*************************************************************************/

#ifndef   GrdAcq_cc_			// Is file already included?
#  define GrdAcq_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Gradients/GrdAcquire.h>	// Include the interface
#include <HSLib/HSacquire.h>		// Include Hilbert space stuff
#include <Level2/acquire1D.h>		// Include class acquire1D
#include <vector>			// Include libstdc++ STL vectors

using std::vector;			// Using libstdc++ STL vectors

// ____________________________________________________________________________
// A                            FID Functions
// ____________________________________________________________________________

/* These functions generate time domain acquisitions based on the GAMMA FID
   functions. That is, they require user specification of the Hamiltonians
   that evolve the system (all subsystems), the initial state of the system
   (all subsystems), the detection operator, and the (dwell) time between
   points in the acquisition.  The functions will automatically sum over all
   sub-systems in the field gradient. Note that the format and names of these
   functions mirror those found in the GAMMA Hilbert Space Library Module
   (src/HSLib/HSacquire.*)  Note thate here each subsystem is assumed ot have
   its own Hamiltonian which may alter the system transitions by more than
   simply off-setting their frequencies.                                     */

// ----------------------------------------------------------------------------
//              Functions Which Take A row_vector As An Argument
//                               (acquire == FID)
// ----------------------------------------------------------------------------

void acquire(vector<gen_op>& sigmas, gen_op& D, vector<gen_op>& Hs,
                                         double td, int t2pts, row_vector& fid)
  { FID(sigmas, D, Hs, td, t2pts, fid); }

void acquire(vector<gen_op>& sigmas, gen_op& D, vector<gen_op>& Us,
                                                    int t2pts, row_vector& fid)
  { FID(sigmas, D, Us, t2pts, fid); }

void FID(vector<gen_op>& sigmas, gen_op& D, vector<gen_op>& Hs,
                                         double td, int t2pts, row_vector& fid)
  {
  unsigned nss = sigmas.size();			// # of sub-systems
  row_vector data0(fid.size(),complex0);	// Empty data block
  fid = data0;					// Zero input data block
  row_vector data, datasum = data0; 		// Working data vectors
  for(unsigned i=0; i<nss; i++)			// Loop over spin systems
    {
    data = data0;				//   Zero working block
    FID(sigmas[i],D,Hs[i],td,t2pts,data);	//   Perform acquisition 
    datasum += data;				//   Sum with previous
    }
  if(data.size() >= fid.size()) fid = data;	// Copy to fid if block bigger
  else						// Else copy only generated
    {						// acquisition points
    for(int j=0; j<t2pts; j++)
      fid.put(data.get(j), j);
    }
  }

void FID(vector<gen_op>& sigmas, gen_op& D, vector<gen_op>& Us,
                                                    int t2pts, row_vector& fid)
  {
  unsigned nss = sigmas.size();			// # of sub-systems
  row_vector data0(fid.size(),complex0);	// Empty data block
  fid = data0;					// Zero input data block
  row_vector data, datasum = data0; 		// Working data vectors
  for(unsigned i=0; i<nss; i++)			// Loop over spin systems
    {
    data = data0;				//   Zero working block
    FID(sigmas[i],D,Us[i],t2pts,data);		//   Perform acquisition 
    datasum += data;				//   Sum with previous
    }
  if(data.size() >= fid.size()) fid = data;	// Copy to fid if block bigger
  else						// Else copy only generated
    {						// acquisition points
    for(int j=0; j<t2pts; j++)
      fid.put(data.get(j), j);
    }
  }

// ----------------------------------------------------------------------------
//                      Functions Which Return A row_vector
//                                (acquire == FID)
// ----------------------------------------------------------------------------

row_vector acquire(vector<gen_op>& sigmas, gen_op& D, vector<gen_op>& Hs,
                                                         double td, int t2pts)
  { return FID(sigmas, D, Hs, td, t2pts); }

row_vector acquire(vector<gen_op>& sigmas, gen_op& D, vector<gen_op>& Us,
                                                                    int t2pts)
  { return FID(sigmas, D, Us, t2pts); }

row_vector FID(vector<gen_op>& sigmas, gen_op& D, vector<gen_op>& Hs,
                                                         double td, int t2pts)
  {
  unsigned nss = sigmas.size();			// # of sub-systems
  row_vector data(t2pts,complex0);		// Data block for acquire
  for(unsigned i=0; i<nss; i++)			// Loop over spin systems
    data += FID(sigmas[i],D,Hs[i],td,t2pts);	// and sum each acquisiiton 
  return data;					// Return the propagators
  }

row_vector FID(vector<gen_op>& sigmas, gen_op& D, vector<gen_op>& Us,
                                                                    int t2pts)
  {
  unsigned nss = sigmas.size();			// # of sub-systems
  row_vector data(t2pts,complex0);		// Data block for acquire
  for(unsigned i=0; i<nss; i++)			// Loop over spin systems
    data += FID(sigmas[i],D,Us[i],t2pts);	// and sum each acquisiiton 
  return data;					// Return the propagators
  }

// ----------------------------------------------------------------------------
//               Functions Which Include Phenomenological Relaxation
//                                (acquire == FID)
// ----------------------------------------------------------------------------

row_vector acquire(RBasic& RB, vector<gen_op>& sigmas,
                           gen_op& D, vector<gen_op>& Hs, double td, int t2pts)
  { return FID(RB,sigmas,D,Hs,td,t2pts); }

row_vector FID(RBasic& RB, vector<gen_op>& sigmas,
                           gen_op& D, vector<gen_op>& Hs, double td, int t2pts)
  {
  RB.SetDet(D);					// Set detection operator
  int nss = sigmas.size();			// Number of subsystems
  row_vector data(t2pts, complex0);		// Data block for FID
  for(int i=0; i<nss; i++)			// Loop over all subsystems
    {
    RB.SetH0(Hs[i]);                            //   Set Hamiltonian in FID
    data += RB.FID(sigmas[i],td,t2pts); 	//   Acquisition
    }
  return data;
  }

// ____________________________________________________________________________
// B                       Class Acquire1D Functions
// ____________________________________________________________________________

/* The funcitons in the previous section essentially duplicated the GAMMA FID
   functionality when there is no gradient. The functions in this section 
   cannot duplicate those provided by class acquire1D because we cannot be
   member functions of that class. Rather we just make use of the class as
   best we can to provide users with the main abilities if acquire1D.        */

row_vector acquire1DT(gen_op& D, vector<gen_op>& Hs,
                      vector<gen_op>& sigmas, int t2pts, double td, bool norm)
  {
  unsigned nss = sigmas.size();			// # of sub-systems
  row_vector data;				// Data block for acquire
  for(unsigned i=0; i<nss; i++)			// Loop over spin systems
    {
    acquire1D ACQ(D, Hs[i]);			//   Prepare for acquisition
    data += ACQ.T(sigmas[i],t2pts,-td);		//   Sum up sub-system FID's
    }
  if(norm) data /= nss;				// Normalized if desired
  return data;					// Return the propagators
  }

row_vector acquire1DT(gen_op& D, vector<gen_op>& Hs,
                               gen_op& sigmas, int t2pts, double td, bool norm) 
  {
  unsigned nss = Hs.size();			// # of sub-systems
  row_vector data;				// Data block for acquire
  for(unsigned i=0; i<nss; i++)			// Loop over spin systems
    {
    acquire1D ACQ(D, Hs[i]);			//   Prepare for acquisition
    data += ACQ.T(sigmas,t2pts,-td);		//   Sum up sub-system FID's
    }
  if(norm) data /= nss;				// Normalized if desired
  return data;					// Return the propagators
  }

// ____________________________________________________________________________
// X                    Single Point Detection Functions
// ____________________________________________________________________________

complex detect(gen_op& D, vector<gen_op>& sigmas)
  {
  unsigned nss = sigmas.size();			// # of sub-systems
  complex z(0);					// Point that is detected
  for(unsigned i=0; i<nss; i++)			// Loop over spin systems
    z += trace(D, sigmas[i]);
  return z;
  }

#endif 								// GrdAcquire.cc

