/* BlochB.cc ****************************************************-*-c++-*-
**									**
**                             G A M M A				**
**									**
**      Bloch Equation B1 and Offset Matrix        Implementation	**
**									**
**      Copyright (c) 1995						**
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

#ifndef   BlochB_cc_ 			// Is file already included?
#  define BlochB_cc_ 1      		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation   
#  endif

#include <Bloch/BlochB.h>		// Include our interface
#include <Bloch/Bloch.h>		// Include Bloch auxiliary
#include <Basics/Gutils.h>		// Include paramter query
#include <Basics/StringCut.h>		// Include Gform function

// ____________________________________________________________________________
// A                Single Component RF-Field Matrix Functions
// ____________________________________________________________________________

/* These are just direct functions which return 3x3 matrices applicable to
   a "single" magnetization vector with only components {Mx, My, Mz}. In such
   a case the user can just input the field strength, offset, and phase.

	   Input	  gamB1 : Applied RF field strength (Hz)
	        	    Wrf : Applied RF offset         (Hz)
                            phi : Applied RF phase          (degrees)
                            argc: Number of arguments
                            argv: Array of arguments
                            qn  : Query index
	   Output	      R	: 3x3 Bloch offset & field matrix
	   Note			: B is output in 1/sec                       */

BlochMx BlochB(double gamB1, double w, double phi)
  {
  matrix B(3, 3, complex0);
  gamB1 *= HZ2RAD; 			// Switch rf strength to rad/sec
  phi   *= DEG2RAD;			// Switch phase to radians
  w     *= HZ2RAD;			// Switch offset to rad/sec
  double GBSphi = gamB1*sin(phi);
  double GBCphi = gamB1*cos(phi);
  B.put(-w,0,1);			// Set <1|K|2> to -w
  B.put(GBSphi,0,2);			// <1|K|3> to gamB1*sin(phi)
  B.put(w,1,0);				// Set <2|K|1> to w
  B.put(-GBCphi,1,2);			// <2|K|3> to -gamB1*cos(phi)
  B.put(-GBSphi,2,0);			// <3|K|1> to -gamB1*sin(phi)
  B.put( GBCphi,2,1);			// <3|K|2> to  gamB1*cos(phi)
  return BlochMx(B);
  }

// ____________________________________________________________________________
// B                 Two Component RF-Field Matrix Functions
// ____________________________________________________________________________

/* These are just direct functions which return 6x6 matrices applicable to
   a magnetization vector with two sub-vectors: |M> = ||M1>|M2>> having the
   components {M1x, M1y, M1z, M2x, M2y, M2z}. In such a case the user can just
   respective field strengths, offsets, and field phases.

	   Input	  gamB1 : Applied RF field strength (Hz)
	        	     W1 : Offset of 1st sub-vector  (Hz)
	        	     W2 : Offset of 2nd sub-vector  (Hz)
                            phi : Applied RF phase          (degrees)
                            argc: Number of arguments
                            argv: Array of arguments
                            qn  : Query index
	   Output	      R	: 3x3 Bloch offset & field matrix
	   Note			: B is output in 1/sec                       */

BlochMx BlochB(double gamB1, double w1, double w2, double phi)
  {
  matrix B(6, 6, complex0);
  gamB1 *= HZ2RAD; 			// Switch rf strength to rad/sec
  phi   *= DEG2RAD;			// Switch phase to radians
  w1    *= HZ2RAD;			// Switch 1st offset to rad/sec
  w2    *= HZ2RAD;			// Switch 2nd offset to rad/sec

  B.put(-w1,0,1);			// Set <0|K|1> to -w1
  B.put(gamB1*PIx2*sin(phi),0,2);	// <0|K|2> to gamB1*sin(phi)
  B.put(w1,1,0);			// Set <1|K|0> to w1
  B.put(-gamB1*PIx2*cos(phi),1,2);	// <1|K|2> to -gamB1*cos(phi)
  B.put(-gamB1*PIx2*sin(phi),2,0);	// <2|K|0> to -gamB1*sin(phi)
  B.put( gamB1*PIx2*cos(phi),2,1);	// <2|K|1> to  gamB1*cos(phi)

  B.put(-w2,3,4);			// Set <3|K|4> to -w2
  B.put(gamB1*PIx2*sin(phi),3,5);	// <3|K|5> to gamB1*sin(phi)
  B.put(w2,4,3);			// Set <4|K|3> to w2
  B.put(-gamB1*PIx2*cos(phi),4,5);	// <4|K|5> to -gamB1*cos(phi)
  B.put(-gamB1*PIx2*sin(phi),5,3);	// <5|K|3> to -gamB1*sin(phi)
  B.put( gamB1*PIx2*cos(phi),5,4);	// <5|K|4> to  gamB1*cos(phi)
  return BlochMx(B);
  }

BlochMx BlochB(double gamB11, double w1, double phi1,
               double gamB12, double w2, double phi2)
  {
  matrix B(6, 6, complex0);
  gamB11 *= HZ2RAD;			// Switch 1st rf strength to rad/sec
  gamB12 *= HZ2RAD;			// Switch 2nd rf strength to rad/sec
  phi1   *= DEG2RAD;			// Switch 1st phase to radians
  phi2   *= DEG2RAD;			// Switch 2nd phase to radians
  w1     *= HZ2RAD;			// Switch 1st offset to rad/sec
  w2     *= HZ2RAD; 			// Switch 2nd offset to rad/sec

  B.put(-w1,0,1);                       // Set <0|K|1> to -w1
  B.put(gamB11*PIx2*sin(phi1),0,2);	// <0|K|2> to gamB11*sin(phi1)
  B.put(w1,1,0);                        // Set <1|K|0> to w1
  B.put(-gamB11*PIx2*cos(phi1),1,2);	// <1|K|2> to -gamB11*cos(phi1)
  B.put(-gamB11*PIx2*sin(phi1),2,0);	// <2|K|0> to -gamB11*sin(phi1)
  B.put(gamB11*PIx2*cos(phi1),2,1);	// <2|K|1> to  gamB11*cos(phi1)

  B.put(-w2,3,4);                       // Set <3|K|4> to -w2
  B.put(gamB12*PIx2*sin(phi2),3,5);	// <3|K|5> to gamB12*sin(phi2)
  B.put(w2,4,3);                        // Set <4|K|3> to w2
  B.put(-gamB12*PIx2*cos(phi2),4,5);	// <4|K|5> to -gamB12*cos(phi2)
  B.put(-gamB12*PIx2*sin(phi2),5,3);	// <5|K|3> to -gamB12*sin(phi2)
  B.put(gamB12*PIx2*cos(phi2),5,4);	// <5|K|4> to  gamB12*cos(phi2)
  return BlochMx(B);
  }

// ____________________________________________________________________________
// C                 Multi-Component RF-Field Matrix Functions
// ____________________________________________________________________________

BlochMx BlochB(std::vector<double> gamB1s, std::vector<double> Ws)
  {
  int N = gamB1s.size();
  std::vector<double> phis(N,0.0);
  return BlochB(gamB1s, Ws, phis);
  }

BlochMx BlochB(std::vector<double> gamB1s, std::vector<double> Ws,
                                                     std::vector<double> phis)
  {
  int N = gamB1s.size();
  matrix B(3*N, 3*N, complex0);
  double w, gB1, phi;
  for(int i=0; i<N; i++)
    {
    w   = Ws[i]     * HZ2RAD;
    gB1 = gamB1s[i] * HZ2RAD;
    phi = phis[i]   * DEG2RAD;
    B.put(-w,                 3*i,   3*i+1);
    B.put(gB1*PIx2*sin(phi),  3*i,   3*i+2);
    B.put(w,                  3*i+1, 3*i);
    B.put(-gB1*PIx2*cos(phi), 3*i+1, 3*1+2);
    B.put(-gB1*PIx2*sin(phi), 3*i+2, 3*i);
    B.put(gB1*PIx2*cos(phi),  3*i+2, 3*i+1);
    }
  return BlochMx(B);
  }

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

BlochMx BlochB(const BlochSys& sys) { return sys.B(); }

// ____________________________________________________________________________
// E                        Interactive Functions
// ____________________________________________________________________________

BlochMx BlochB(int argc, char* argv[], double& gamB1, double& w, double& phi, int& qn)
  {
  int magB = DoubleMag(gamB1);
  int magW = DoubleMag(w);
  double sfB, sfW;
  std::string  uB,  uW;
  uB = HzUnits(magB, sfB);
  uW = HzUnits(magW, sfW);
  std::string msgB = "\n\tRF-Field Strength (gamma*B1) in " + uB
              + " [" + Gform("%8.4f", gamB1*sfB) + "]? ";
  std::string msgW = "\n\tVector Offset Frequency in " + uW
              + " [" + Gform("%8.4f", w*sfW) + "]? ";
  ask_set(argc, argv, qn++, msgB, gamB1);
  if(gamB1)
    {
    std::string msgP = "\n\tRF-Field Phase Angle (0=x) in Degrees ["
       + Gform("%8.4f", phi) + "]? ";
    ask_set(argc, argv, qn++, msgP, phi);
    }
  ask_set(argc, argv, qn, msgW, w);
  return BlochB(gamB1, w, phi);
  }

// ____________________________________________________________________________
// F                         Deprecated Functions
// ____________________________________________________________________________

/* These are marked for removal in a future version. They will work but will
   issue warning messages.                                                   */

matrix K_matrix(matrix& R, double gamB1, double w, double phi)
  {
  std::string hdr("Bloch Module");
  GAMMAerror(hdr, 4, "K_matrix", 1);
  GAMMAerror(hdr, 6, "BlochB",   1);
  GAMMAerror(hdr, 6, "BlochR");
  return R + BlochB(gamB1, w, phi);
  }

matrix K_matrix(int argc, char* argv[], matrix& R,
                     double& gamB1, double& w, double& phi, int& qn)
  {
  std::string hdr("Bloch Module");
  GAMMAerror(hdr, 4, "K_matrix", 1);
  GAMMAerror(hdr, 6, "BlochB",   1);
  GAMMAerror(hdr, 6, "BlochR");
  return R + BlochB(argc, argv, gamB1, w, phi, qn);
  }

#endif								// BlochB.cc

