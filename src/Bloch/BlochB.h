/* BlochB.h *****************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**      Bloch Equation B1 and Offset Matrix             Interface	**
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
** This file contains functions which generate basic evolution matrices	**
** for use in magnetization vector evolution under the phenomenological	**
** Bloch equations. In the simplest case the returned array will be	**
** a 3x3 matrix which appears as					**
**									**
**          [        0          -w + w       g*B  * sin(phi) ]		**
**          |                     0   rf        1            |		**
**          |                                                |		**
**      B = |     w - w             0       -g*B  * cos(phi) |		**
**          |      0   rf                       1            |		**
**          |                                                |		**
**          [ -gB * sin(phi)   gB * cos(phi)       0         |		**
**          [    1               1                           ]		**
**                                                  t			**
** and this will act on the magnetization vector |M> = [Mx My Mz].	**
** In a more general context, the above array will be a single block	**
** on the diagonal of a larger matrix of dimension 3N x 3N where N is	**
** the number of sub-vectors in a general magnetization vector . In	**
** that case, the array appears as					**
**									**
**                        [ [B ]  0    0   . . .   0    ]		**
**                        | [ 0]                        |		**
**                        |      [B ]  0   . . .   0    |		**
**                        |  0   [ 1]                   |		**
**                        |           [B ] . . .   0    |		**
**                    B = |  0    0   [ 2]              |		**
**                        |                . . .   0    |		**
**                        |  .    .    .   . . .        |		**
**                        |  .    .    .   . . .   0    |		**
**                        |  .    .    .   . . .        |		**
**                        |                . . . [B   ] |		**
**                        [  0    0    0   . . . [ N-1] ]		**
**									**
**                                             t			**
** and will act on the magnetization vector |M> = ||M >|M >....|M   >.	**
**                                                   0	 1       N-1	**
**									**
*************************************************************************/

#ifndef   BlochB_h_ 			// Is file already included?
#  define BlochB_h_ 1      		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface   
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Bloch/BlochMx.h>		// Know about Bloch matrices
#include <Matrix/matrix.h>		// Know about GAMMA matrices
#include <Bloch/BlochSys.h>		// Know about GAMMA Bloch systems
#include <string>			// Know about libstdc++ strings
#include <vector>			// Know about libstdc++ STL vectors

// ____________________________________________________________________________
// A                Single Component RF-Field Matrix Functions
// ____________________________________________________________________________

/* These are just direct functions which return 3x3 matrices applicable to
   a "single" magnetization vector with only components {Mx, My, Mz}. In such
   a case the user can just input the field strength, offset, and phase.

           Input          gamB1 : Applied RF field strength (Hz)
                            Wrf : Applied RF offset         (Hz)
                            phi : Applied RF phase          (degrees)
                            argc: Number of arguments
                            argv: Array of arguments
                            qn  : Query index
           Output             R : 3x3 Bloch offset & field matrix
           Note                 : B is output in 1/sec                       */

MSVCDLL BlochMx BlochB(double gamB1, double w, double phi=0.0);

// ____________________________________________________________________________
// B                 Two Component RF-Field Matrix Functions
// ____________________________________________________________________________

/* These are just direct functions which return 6x6 matrices applicable to
   a magnetization vector with two sub-vectors: |M> = ||M1>|M2>> having the
   components {M1x, M1y, M1z, M2x, M2y, M2z}. In such a case the user can just
   respective field strengths, offsets, and field phases.

           Input          gamB1 : Applied RF field strength (Hz)
                             W1 : Offset of 1st sub-vector  (Hz)
                             W2 : Offset of 2nd sub-vector  (Hz)
                            phi : Applied RF phase          (degrees)
                            argc: Number of arguments
                            argv: Array of arguments
                            qn  : Query index
           Output             R : 3x3 Bloch offset & field matrix
           Note                 : B is output in 1/sec                       */

MSVCDLL BlochMx BlochB(double gamB1,  double w1, double w2, double phi);
MSVCDLL BlochMx BlochB(double gamB11, double w1,            double phi1,
               double gamB12, double w2,            double phi2=0.0);

// ____________________________________________________________________________
// C                 Multi-Component RF-Field Matrix Functions
// ____________________________________________________________________________

MSVCDLL BlochMx BlochB(std::vector<double> gamB1s, std::vector<double> Ws);
MSVCDLL BlochMx BlochB(std::vector<double> gamB1s, std::vector<double> Ws,
                                                     std::vector<double> phis);

// ____________________________________________________________________________
// D             RF-Field Functions Defined Over A Bloch System
// ____________________________________________________________________________

/* These are functions that will return 3Nx3N matrices where N is the number
   of sub-vectors contained in the magnetization vector. In this case the
   magnetization vector contains N sub-vectors each of which has components
   {Mxi, Myi, Mzi}. In Bloch system keeps track of how many sub-vectors and
   their respective offsets, applied fields, & field phases. The magnetization
   vector which this is intended to evolve appears as

                                t
                             |M>  = ||M >|M >....|M   >
                                       0   1       N-1

           Input            sys : Bloch system
           Output             B : 3Nx3N Bloch field & offset matrix
           Note                 : B is output in 1/sec & blocked             */

MSVCDLL BlochMx BlochB(const BlochSys& sys);

// ____________________________________________________________________________
// E                        Interactive Functions
// ____________________________________________________________________________


MSVCDLL BlochMx BlochB(int argc, char* argv[], double& gamB1, double& w, double& phi, int& qn);

// ____________________________________________________________________________
// F                         Deprecated Functions
// ____________________________________________________________________________

/* These are marked for removal in a future version. They will work but will
   issue warning messages.                                                   */

MSVCDLL matrix K_matrix(matrix& R, double gamB1, double w, double phi);
MSVCDLL matrix K_matrix(int argc, char* argv[], matrix& R,
                               double& gamB1, double& w, double& phi, int& qn);

#endif								// BlochB.h
