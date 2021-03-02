/* BlochR.h *****************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**      Bloch Equation  Relaxation Matrix             Interface		**
**									**
**      Copyright (c) 2002						**
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
** This file contains functions which generate relaxation matrices for  **
** use in magnetization vector evolution under the phenomenological     **
** Bloch equations. In the simplest case the returned array will be     **
** a 3x3 matrix which appears as                                        **
**                                                                      **
**                               [ R   0   0  ]                         **
**                               |  2         |                         **
**                               |            |                         **
**                           R = | 0   R   0  |                         **
**                               |      2     |                         **
**                               |            |                         **
**                               | 0   0   R  |                         **
**                               [          1 ]                         **
**                                                  t                   **
** and this will act on the magnetization vector |M> = [Mx My Mz].      **
** In a more general context, the above array will be a single block    **
** on the diagonal of a larger matrix of dimension 3N x 3N where N is   **
** the number of sub-vectors in a general magnetization vector . In     **
** that case, the array appears as                                      **
**                                                                      **
**                        [ [R ]  0    0   . . .   0    ]               **
**                        | [ 0]                        |               **
**                        |      [R ]  0   . . .   0    |               **
**                        |  0   [ 1]                   |               **
**                        |           [R ] . . .   0    |               **
**                    R = |  0    0   [ 2]              |               **
**                        |                . . .   0    |               **
**                        |  .    .    .   . . .        |               **
**                        |  .    .    .   . . .   0    |               **
**                        |  .    .    .   . . .        |               **
**                        |                . . . [R   ] |               **
**                        [  0    0    0   . . . [ N-1] ]               **
**                                                                      **
**                                             t                        **
** and will act on the magnetization vector |M> = ||M >|M >....|M   >.  **
**                                                   0   1       N-1    **
**                                                                      **
*************************************************************************/

#ifndef   BlochR_h_ 			// Is file already included?
#  define BlochR_h_ 1      		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface  
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Bloch/BlochMx.h>		// Know about Bloch matrices
#include <Matrix/matrix.h>		// Know about GAMMA matrices
#include <Bloch/BlochSys.h>		// Know about GAMMA Bloch systems
#include <string>			// Know about libstdc++ strings

// ____________________________________________________________________________
// A                Single Component Relaxation Matrix Functions
// ____________________________________________________________________________

/* These are just direct functions which return 3x3 matrices applicable to
   a "single" magnetization vector with only components {Mx, My, Mz}. In such
   a case the user can just input the transverse and longitudinal relaxation
   times. The returned array will be diagonal.

           Input             T1 : Longitudinal relaxation time/rate
                             T2 : Transverse relaxation time/rate
			    rate: Flag if values are rates or times
           Output             R : 3x3 Bloch relaxation matrix
           Note                 : T1 & T2 values in seconds
           Note                 : R is output in 1/sec
           Note                 : Zero or negative relaxation
                                  rates result in a zero R                   */

MSVCDLL BlochMx AFBeforeR(double, double, bool rate=false);
MSVCDLL BlochMx BlochR(double T1, double T2, bool rate=false);

// ____________________________________________________________________________
// B                 Two Component Relaxation Matrix Functions
// ____________________________________________________________________________

/* These are just direct functions which return 6x6 matrices applicable to
   a magnetization vector with two sub-vectors: |M> = ||M1>|M2>> having the
   components {M1x, M1y, M1z, M2x, M2y, M2z}. In such a case the user can just
   respective field transverse and longitudinal rates.

           Input            T11 : 1st Longitudinal relaxation rate
                            T21 : 1st Transverse relaxation rate
                            T12 : 2nd Longitudinal relaxation rate
                            T22 : 2nd Transverse relaxation rate
                            rate: Flag if values are rates or times
                            argc: Number of arguments
                            argv: Array of arguments
                            qn  : Query index
           Output             R : 6x6 Bloch relaxation matrix
           Note                 : T1 & T2 values in seconds
           Note                 : R is output in 1/sec
           Note                 : Zero or negative relaxation
                                  rates result in a zero R                   */

MSVCDLL BlochMx BlochR(double T11, double T12, double T21, double T22, bool rate=false);

// ____________________________________________________________________________
// C                 Multi-Component Relaxation Matrix Functions
// ____________________________________________________________________________

/* These allow users to set up an R matrix for any number of components.

	   Input	    T1s : Longitudinal relaxation rate/time
	        	    T2s : Transverse relaxation rate/time
                            rate: Flag if values are rates or times
	   Output	      R	: 3Nx2N Bloch relaxation matrix where N
				  is the size of vectors T1s and T2s
	   Note			: T1 & T2 values in seconds
	   Note			: R is output in 1/sec
	   Note			: Zero or negative relaxation
	  			  rates result in a zero R                   */

MSVCDLL BlochMx BlochR(std::vector<double> T1s, std::vector<double> T2s, bool rate=false);

// ____________________________________________________________________________
// D            Relaxation Matrix Functions Defined Over A Bloch System
// ____________________________________________________________________________

/* These are functions that will return 3Nx3N matrices where N is the number
   of sub-vectors contained in the magnetization vector. In this case the
   magnetization vector contains N sub-vectors each of which has components
   {Mxi, Myi, Mzi}. In Bloch system keeps track of how many sub-vectors and
   their relaxation rates. The magnetization vector which this is intended
   to evolve appears as
                                t
                             |M>  = ||M >|M >....|M   >
                                       0   1       N-1

           Input            sys : Bloch system
           Output             R : 3Nx3N Bloch relaxation matrix
           Note                 : R is output in 1/sec & diagonal
           Note                 : Zero or negative relaxation
                                  rates result in a zero R                   */

MSVCDLL BlochMx Bloch_R(const BlochSys& sys);

// ____________________________________________________________________________
// E                        Interactive Functions
// ____________________________________________________________________________

MSVCDLL BlochMx BlochR(int argc, char* argv[], double& T1, double& T2, int& qn, bool rate=false);

// ____________________________________________________________________________
// F                         Deprecated Functions
// ____________________________________________________________________________

/* These are marked for removal in a future version. They will work but will
   issue warning messages.                                                   */

MSVCDLL matrix R_matrix(double T1, double T2);
MSVCDLL matrix R_matrix(int argc, char* argv[], double& T1, double& T2, int& qn);

#endif						// BlochR.h
