/* BlochK.h *****************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**      Bloch Equation Exchange Matrix             Interface		**
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
** This file contains functions which generate exchange matrices for	**
** use in magnetization vector evolution under the phenomenological     **
** Bloch equations. In the simplest case the returned array will be     **
** a 6x6 matrix which appears as                                        **
**                                                                      **
**                          [  k   0   0  -k   0   0   ]		**
**                          |   12          12         |		**
**                          |                          |		**
**                          |  0   k   0   0  -k   0   |		**
**                          |       12          12     |		**
**                          |                          |		**
**                          |  0   0   k   0   0  -k   |		**
**                          |           12          12 |		**
**                    K   = |                          |		**
**                     12   | -k   0   0   k   0   0   |		**
**                          |   12          12         |		**
**                          |                          |		**
**                          |  0  -k   0   0   k   0   |		**
**                          |       12          12     |		**
**                          |                          |		**
**                          |  0   0  -k   0   0   k   |		**
**                          |           12          12 |		**
**                                                                      **
**                     t						**
** This will act on |M  > = ||M >|M >> = [ M1x M1y M1z M2x M2y M2z].	**
**                    12       1   2					**
**                                                                      **
** In a more general context, elements of the above array will be	**
** coorespond to exchange between two sub-vectors & where the elements 	**
** occur depends upon placement of the two sub-vectors in the larger	**
** magnetization vector. The returned array will have dimension 3N x 3N	**
** where N is the number of sub-vectors in a general magnetization	**
** vector . This array is best viewed as a sum over arrays of the same	**
** dimension for each sub-vector pair.					**
**                                                                      **
**                                     N				**
**                                    ---				**
**                               K =  \    K				**
**                                    /     ij				**
**                                    ---      				**
**                                   i, j>i				**
**                                                                      **
** The diagonal elements of the matrix for a sub-vector pair occur	**
** in two sets of triples and are given by 				**
**                                                                      **
**    <3*i|K  |3*i> = <3*i+1|K  |3*i+1> = <3*i+2|K  |3*i+2> = k		**
**          ij                ij                  ij           ij	**
**                                                                      **
**    <3*j|K  |3*j> = <3*j+1|K  |3*j+1> = <3*j+2|K  |3*j+2> = k		**
**          ij                ij                  ij           ij	**
**                                                                      **
** The off diagonals of the sub-vector pair arrays occur where there is	**
** a match up of Mix<->Mjx, Miy<->Mjy, Miz<->Mjz and there will be six	**
** of these as given by							**
**                                                                      **
**    <3*i|K  |3*j> = <3*i+1|K  |3*j+1> = <3*i+2|K  |3*j+2> = -k	**
**          ij                ij                  ij            ij	**
**                                                                      **
**    <3*j|K  |3*i> = <3*j+1|K  |3*i+1> = <3*j+2|K  |3*i+2> = -k	**
**          ij                ij                  ij            ij	**
**                                                                      **
** All other elemnts will be zero.					**
**                                                                      **
*************************************************************************/

#ifndef   BlochK_h_ 			// Is file already included?
#  define BlochK_h_ 1      		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface  
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/matrix.h>		// Know about GAMMA matrices
#include <Bloch/BlochSys.h>		// Know about GAMMA Bloch systems
#include <string>			// Know about libstdc++ strings


// ----------------------------------------------------------------------------
//                      Simple Exchange Matrix Functions
// ----------------------------------------------------------------------------

/* These are just direct functions which return 6x6 matrices applicable to
   a magnetization vector with two sub-vectors: |M> = ||M1>|M2>> having the
   components {M1x, M1y, M1z, M2x, M2y, M2z}. In such a case the user can just
   input the exchange rate. The returned array will be real symmetric.

           Input              k : Exchange rate (1/sec)
           Output             K : 6x6 Bloch exchange matrix
           Note                 : K is output in 1/sec                       */

MSVCDLL matrix BlochK(double k);

// ----------------------------------------------------------------------------
//             Exchange Matrix Functions Defined Over A Bloch System
// ----------------------------------------------------------------------------

/* These are functions that will return 3Nx3N matrices where N is the number
   of sub-vectors contained in the magnetization vector. In this case the
   magnetization vector contains N sub-vectors each of which has components
   {Mxi, Myi, Mzi}. In Bloch system keeps track of how many sub-vectors and
   their eachange rates. The magnetization vector which this is intended
   to evolve appears as
                                t
                             |M>  = ||M >|M >....|M   >
                                       0   1       N-1

           Input            sys : Bloch system
           Output             K : 3Nx3N Bloch excahnge matrix
           Note                 : K is output in 1/sec & real symmetric     */

MSVCDLL matrix BlochK(const BlochSys& sys);

#endif						// BlochK.h
