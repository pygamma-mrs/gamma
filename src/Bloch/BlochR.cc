/* BlochR.cc ****************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**      Bloch Equation Relaxation Matrix             Implementation	**
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
** This file contains functions which generate relaxation matrices for	**
** use in magnetization vector evolution under the phenomenological 	**
** Bloch equations. In the simplest case the returned array will be	**
** a 3x3 matrix which appears as					**
**									**
**                               [ R   0   0  ]                         **
**                               |  2         |                         **
**                               |            |                         **
**                           R = | 0   R   0  |                         **
**                               |      2     |                         **
**                               |            |                         **
**                               | 0   0   R  |                         **
**                               [          1 ]                         **
**                                                  t			**
** and this will act on the magnetization vector |M> = [Mx My Mz].	**
** In a more general context, the above array will be a single block	**
** on the diagonal of a larger matrix of dimension 3N x 3N where N is	**
** the number of sub-vectors in a general magnetization vector . In	**
** that case, the array appears as					**
**									**
**                        [ [R ]  0    0   . . .   0    ]		**
**                        | [ 0]                        |		**
**                        |      [R ]  0   . . .   0    |		**
**                        |  0   [ 1]                   |		**
**                        |           [R ] . . .   0    |		**
**                    R = |  0    0   [ 2]              |		**
**                        |                . . .   0    |		**
**                        |  .    .    .   . . .        |		**
**                        |  .    .    .   . . .   0    |		**
**                        |  .    .    .   . . .        |		**
**                        |                . . . [R   ] |		**
**                        [  0    0    0   . . . [ N-1] ]		**
**									**
**                                             t			**
** and will act on the magnetization vector |M> = ||M >|M >....|M   >.	**
**                                                   0	 1       N-1	**
**									**
*************************************************************************/

#ifndef   BlochR_cc_ 			// Is file already included?
#  define BlochR_cc_ 1      		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation   
#  endif

#include <Bloch/BlochR.h>		// Include our interface
#include <Bloch/Bloch.h>		// Include Bloch auxiliary stuff
#include <Basics/Gutils.h>		// Include paramter query
#include <Basics/StringCut.h>		// Include Gform function

// ____________________________________________________________________________
// A                Single Component Relaxation Matrix Functions
// ____________________________________________________________________________

/* These are just direct functions which return 3x3 matrices applicable to
   a "single" magnetization vector with only components {Mx, My, Mz}. In such
   a case the user can just input the transverse and longitudinal relaxation
   times. The returned array will be diagonal.

	   Input	     T1 : Longitudinal relaxation rate
	        	     T2 : Transverse relaxation rate
                            rate: Flag if values are rates or times
                            argc: Number of arguments
                            argv: Array of arguments
                            qn  : Query index
	   Output	      R	: 3x3 Bloch relaxation matrix
	   Note			: T1 & T2 values in seconds
	   Note			: R is output in 1/sec
	   Note			: Zero or negative relaxation
	  			  rates result in a zero R                   */

BlochMx AFBeforeR(double T1, double T2, bool rate)
  {
  matrix R(3, 3, 0.0, d_matrix_type);
  if(rate) T1 = T2;
  return R;
  }

BlochMx BlochR(double T1, double T2, bool rate)
  {
  matrix R(3, 3, 0.0, d_matrix_type);
  double R1, R2;
  if(rate)
    {
    R1 = fabs(T1);
    R2 = fabs(T2);
    }
  else
    {
    if(fabs(T1) < 1.e-13) R1 = 0; 
    else                  R1 = 1.0/fabs(T1); 
    if(fabs(T2) < 1.e-13) R2 = 0; 
    else                  R2 = 1.0/fabs(T2); 
    }
  R.put(R2, 0, 0);				// Set <1|R|1> to 1/T2
  R.put(R2, 1, 1);				// Set <2|R|2> to 1/T2
  R.put(R1, 2, 2);				// Set <3|R|3> to 1/T1
  return BlochMx(R);
  }

// ____________________________________________________________________________
// B                 Two Component Relaxation Matrix Functions
// ____________________________________________________________________________

/* These are just direct functions which return 6x6 matrices applicable to
   a magnetization vector with two sub-vectors: |M> = ||M1>|M2>> having the
   components {M1x, M1y, M1z, M2x, M2y, M2z}. In such a case the user can just
   respective field transverse and longitudinal rates.

	   Input	    T11 : 1st Longitudinal relaxation rate
	        	    T21 : 1st Transverse relaxation rate
	   		    T12 : 2nd Longitudinal relaxation rate
	        	    T22 : 2nd Transverse relaxation rate
                            rate: Flag if values are rates or times
                            argc: Number of arguments
                            argv: Array of arguments
                            qn  : Query index
	   Output	      R	: 6x6 Bloch relaxation matrix
	   Note			: T1 & T2 values in seconds
	   Note			: R is output in 1/sec
	   Note			: Zero or negative relaxation
	  			  rates result in a zero R                   */

BlochMx BlochR(double T11, double T12, double T21, double T22, bool rate)
  {
  matrix R(6, 6, 0.0, d_matrix_type);
  double R11, R12, R21, R22;
  if(rate)
    {
    R11 = fabs(T11);
    R12 = fabs(T12);
    R21 = fabs(T21);
    R22 = fabs(T22);
    }
  else
    {
    if(fabs(T11) < 1.e-13) R11 = 0; 
    else                   R11 = 1.0/fabs(T11); 
    if(fabs(T12) < 1.e-13) R12 = 0; 
    else                   R12 = 1.0/fabs(T12); 
    if(fabs(T21) < 1.e-13) R21 = 0; 
    else                   R21 = 1.0/fabs(T21); 
    if(fabs(T22) < 1.e-13) R22 = 0; 
    else                   R22 = 1.0/fabs(T22); 
    }
  R.put(R12, 0, 0);				// Set <1|R|1> to 1/T12
  R.put(R12, 1, 1);				// Set <2|R|2> to 1/T12
  R.put(R11, 2, 2);				// Set <3|R|3> to 1/T11
  R.put(R22, 3, 3);				// Set <4|R|4> to 1/T22
  R.put(R22, 4, 4);				// Set <5|R|5> to 1/T22
  R.put(R21, 5, 5);				// Set <6|R|6> to 1/T21

  return BlochMx(R);
  }

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

BlochMx BlochR(std::vector<double> T1s, std::vector<double> T2s, bool rate)
  {
  int nv = T1s.size();
  matrix R(nv, nv, 0.0, d_matrix_type);
  if(rate)
    {
    for(int i=0; i<nv; i++)
      {
      R.put(fabs(T2s[i]), 3*i,   3*i);		// Set <1|Ri|1> to R2i
      R.put(fabs(T2s[i]), 3*i+1, 3*i+1);	// Set <2|Ri|2> to R2i
      R.put(fabs(T1s[i]), 3*i+2, 3*i+2);	// Set <3|Ri|3> to R1i
      }
    }
  else
    {
    double r2, r1, t2, t1;
    for(int i=0; i<nv; i++)
      {
      t2 = fabs(T2s[i]);
      t1 = fabs(T1s[i]);
      if(t1 < 1.e-13) r1 = 0; 
      else            r1 = 1.0/t1; 
      if(t2 < 1.e-13) r2 = 0; 
      else            r2 = 1.0/t2; 
      R.put(r2, 3*i,   3*i);		// Set <1|Ri|1> to 1/T2i
      R.put(r2, 3*i+1, 3*i+1);		// Set <2|Ri|2> to 1/T2i
      R.put(r1, 3*i+2, 3*i+2);		// Set <3|Ri|3> to 1/T1i
      }
    }
  return BlochMx(R);
  }

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

	   Input	    sys : Bloch system
	   Output	      R	: 3Nx3N Bloch relaxation matrix
	   Note			: R is output in 1/sec & diagonal
	   Note			: Zero or negative relaxation
	  			  rates result in a zero R                   */

BlochMx BlochR(const BlochSys& sys) { return sys.R(); }

// ____________________________________________________________________________
// E                        Interactive Functions
// ____________________________________________________________________________

BlochMx BlochR(int argc, char* argv[], double& T1, double& T2, int& qn, bool rate)
  {
  matrix R(3, 3, 0.0, d_matrix_type);		// R matrix
  std::string msg1, msg2;			// For query text
  int mag1 = DoubleMag(T1);			// Size of T1 default value
  int mag2 = DoubleMag(T2);			// Size of T2 default value
  double sf1, sf2;				// Scaling factors (I/O)
  std::string un1, un2;				// Output units
  double R1, R2;				// Relaxatio rates (1/sec)
  if(rate)
    {
    msg1 = "\n\tLongitudinal Relaxation Rate (R1) in 1/sec [";
    msg2 = "\n\tTransverse   Relaxation Rate (R2) in 1/sec [";
    ask_set(argc, argv, qn++, msg1, R1);
    ask_set(argc, argv, qn,   msg2, R2);
    R1 = fabs(R1);
    R2 = fabs(R2);
    T1 = (R1 > 1.e-13)?1.0/R1:0.0; 
    T2 = (R2 > 1.e-13)?1.0/R2:0.0; 
    }
  else
    {
    un1 = SecUnits(mag1, sf1);
    un2 = SecUnits(mag2, sf2);
    msg1 = "\n\tLongitudinal Relaxation Time (T1) in " + un1
         + "[" + Gform("%8.3f", T1*sf1) + "]? ";
    msg2 = "\n\tTransverse   Relaxation Time (T2) in " + un2
         + "[" + Gform("%8.3f", T2*sf2) + "]? ";
    ask_set(argc, argv, qn, msg1, T1);
    if(fabs(T1) > 1.e-13) 
      {
      qn++;
      T1 = fabs(T1);
      R1 = 1.0/T1;
      ask_set(argc, argv, qn, msg2, T2);
      T2 = fabs(T2);
      if(T2 > 1.e-13) { R2 = 1.0/T2;        }
      else            { R2 = 0.0; T2 = 0.0; }
      }
    else
      {
      T1 = 0.0;
      T2 = 0.0;
      R1 = 0.0;
      R2 = 0.0;
      }
    }
  R.put(R2, 0, 0);				// Set <1|R|1> to 1/T2
  R.put(R2, 1, 1);				// Set <2|R|2> to 1/T2
  R.put(R1, 2, 2);				// Set <3|R|3> to 1/T1
  return BlochMx(R);
  }


// ____________________________________________________________________________
// F                         Deprecated Functions
// ____________________________________________________________________________

/* These are marked for removal in a future version. They will work but will
   issue warning messages.                                                   */

matrix R_matrix(double T1, double T2)
  {
  std::string hdr("Bloch Module");
  GAMMAerror(hdr, 4, "R_matrix", 1);
  GAMMAerror(hdr, 6, "BlochR");
  return BlochR(T1, T2);
  }

matrix R_matrix(int argc, char* argv[], double& T1, double& T2, int& qn)
  {
  std::string hdr("Bloch Module");
  GAMMAerror(hdr, 4, "R_matrix", 1);
  GAMMAerror(hdr, 6, "BlochR");
  return BlochR(argc, argv, T1, T2, qn);
  }

#endif							// BlochR.cc

