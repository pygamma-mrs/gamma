/* BlochAcquire.h **********************************************-*-c++-***
**									**
**	                         G A M M A 				**
**						 			**
**	Bloch Acquisition Module                  	Interface	**
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
**  This GAMMA module provides functions for MR simulations using the	**
**  Bloch equations. In particular, these functions are set to perform 	**
**  acquisitions of evenly spaced magnetization values on a system of	**
**  magnetization vectors.						**
**                                                                      **
*************************************************************************/

#ifndef _BlochAcq_h_			// Is file already included?
#define _BlochAcq_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma interface			// This is the interface
#endif

#include <Bloch/BlochSys.h>		// Know Bloch systems
#include <HSLib/GenOp.h>		// Know general operators
#include <Level2/RelaxBas.h>		// Know simple ad-hoc relaxation
#include <vector>			// Know libstdc++ STL vectors


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
   simply off-setting their frequencies.				     */

// ----------------------------------------------------------------------------
//              Functions Which Take A row_vector As An Argument
//                               (acquire == FID)
// ----------------------------------------------------------------------------

/*
void acquire(vector<gen_op>& sigmas, gen_op& D, vector<gen_op>& Hs,
                                       double td, int t2pts, row_vector& fid);
void acquire(vector<gen_op>& sigmas, gen_op& D, vector<gen_op>& Us,
                                                  int t2pts, row_vector& fid);
*/

void FID(vector<gen_op>&     sigmas, gen_op& D, vector<gen_op>& Hs,
                                       double td, int t2pts, row_vector& fid);
void FID(vector<gen_op>&     sigmas, gen_op& D, vector<gen_op>& Us,
                                                  int t2pts, row_vector& fid);

// ----------------------------------------------------------------------------
//                      Functions Which Return A row_vector
//                                (acquire == FID)
// ----------------------------------------------------------------------------

/*
row_vector acquire(vector<gen_op>& sigmas, gen_op& D, vector<gen_op>& Hs,
                                                         double td, int t2pts);
row_vector acquire(vector<gen_op>& sigmas, gen_op& D, vector<gen_op>& Us,
                                                                    int t2pts);
row_vector FID(vector<gen_op>&     sigmas, gen_op& D, vector<gen_op>& Hs,
                                                         double td, int t2pts);
row_vector FID(vector<gen_op>&     sigmas, gen_op& D, vector<gen_op>& Us,
                                                                    int t2pts);
*/

// ----------------------------------------------------------------------------
//             Functions Which Include Phenomenological Relaxation
//                                (acquire == FID)
// ----------------------------------------------------------------------------

/*

row_vector acquire(RBasic& RB, vector<gen_op>& sigmas,
                          gen_op& D, vector<gen_op>& Hs, double td, int t2pts);
row_vector FID(RBasic&     RB, vector<gen_op>& sigmas,
                          gen_op& D, vector<gen_op>& Hs, double td, int t2pts);
*/

// ____________________________________________________________________________
// B                       Class Acquire1D Functions
// ____________________________________________________________________________

/* The funcitons in the previous section essentially duplicated the GAMMA FID
   functionality when there is no gradient. The functions in this section
   cannot duplicate those provided by class acquire1D because we cannot be
   member functions of that class. Rather we just make use of the class as
   best we can to provide users with the main abilities if acquire1D. 
   As in the last section, the functions here assume that each subsystem has
   its own Hamiltonian which may alter the system transitions by more than 
   simple frequency offsets.						     */

/*
row_vector acquire1DT(gen_op& D, vector<gen_op>& Hs,
                    vector<gen_op>& sigmas, int t2pts, double td, bool norm=1);

row_vector acquire1DT(gen_op& D, vector<gen_op>& Hs,
                            gen_op& sigmas, int t2pts, double td, bool norm=1);

*/
// ____________________________________________________________________________
// D                           Simple FID Functions
// ____________________________________________________________________________

/* These functions generate time domain acquisitions based on the GAMMA FID
   functions. That is, they require user specification of the Hamiltonians
   that evolve the system (all subsystems), the initial state of the system
   (all subsystems), the detection operator, and the (dwell) time between 
   points in the acquisition.  The functions will automatically sum over all
   sub-systems in the field gradient.

   Unlike the functions in section A of this module, these functions do NOT
   assume that, although each sub-system is affected by its own Hamiltonian,
   all systems produce essentially the same transitions but shifted by some
   constant amount due location in the Gradient. This implies that we are in
   a high-field limit so that the Zeeman Hamiltonian contribution dominates
   the system Hamiltonians.                                                  */

// ____________________________________________________________________________
// X                    Single Point Detection Functions
// ____________________________________________________________________________

//complex detect(gen_op& D, vector<gen_op>& sigmas);

#endif 								// BlochAcquire.h

