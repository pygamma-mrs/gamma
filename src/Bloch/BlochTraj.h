/* BlochTraj.h **************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**      Bloch Magnetization Trajectories               Interface	**
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
** This functions provided herein generate magnetization trajectories	**
** over evenly spaced time increments. The users specifies the initial	**
** magnetization vector, the evolution matrix, and the infinite time	**
** magnetization vector. The functions will generate with a row_vector	**
** of magnetization values or else a coord_vec of magnetization 	**
** components.								**
**                                                                      **
** The generalized solution to the Bloch equation(s) that's implemented **
** herein is								**
**                      -Gt						**
**          |M(t)> = exp   [ |M(0)> - |M(inf)> ] + |M(inf)>		**
**                     							**
** where G is the Bloch evolution matrix, |M(0)> the initial vector and	**
** |M(inf)> the vector given infinite evolution time. The functions 	**
** take a particular detection vector <D| and generate a trajectory for	**
** <D| according to							**
**                            -n*t					**
**  D(nt ) = <D|M(t)> = <D|exp    d [ |M(0)> - |M(inf)> ] + <D|M(inf)>	**
**      d    								**
**                                                                      **
** The supplied trajectory functions typically demand that the user	**
** supply the intial and infinite time magnetization vectors, the 	**
** evolution matrix G, the number of trajectory points n, and t		**
** the time increment between points. The detection vector D   d	**
** is usually implicit in the function call.				**
**                                                                      **
*************************************************************************/

#ifndef   BlochTraj_h_ 			// Is file already included?
#  define BlochTraj_h_ 1      		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface  
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/matrix.h>		// Know about GAMMA matrices
#include <Matrix/row_vector.h>		// Know about GAMMA row vectors
#include <Bloch/MagVec.h>		// Know about magnetization vectors
#include <Bloch/DetVec.h>		// Know about detection vectors
#include <Level1/coord_vec.h>		// Know about coordiante vectors
#include <string>			// Know about libstdc++ strings

// ____________________________________________________________________________
// A                       Bloch Detection Vectors
// ____________________________________________________________________________

MSVCDLL DetVec MxDetection(const MagVec& M);		// Detect Mx
MSVCDLL DetVec MyDetection(const MagVec& M);		// Detect My
MSVCDLL DetVec MzDetection(const MagVec& M);		// Detect Mz
MSVCDLL DetVec MpDetection(const MagVec& M);		// Detect Mx+iMy
MSVCDLL DetVec MmDetection(const MagVec& M);		// Detect Mx-iMy

// ____________________________________________________________________________
// B                   Magnetization Trajectories
// ____________________________________________________________________________

MSVCDLL row_vector MTraj(const DetVec& D, const MagVec& Mo,
                                            const matrix& G, int N, double td);
MSVCDLL row_vector MTraj(const DetVec& D, const MagVec& Mo,
                        const matrix& G, const MagVec& Minf, int N, double td);

MSVCDLL row_vector MxTraj(const MagVec& Mo, const matrix& G,                     int N, double td);
MSVCDLL row_vector MxTraj(const MagVec& Mo, const matrix& G, const MagVec& Minf, int N, double td);

MSVCDLL row_vector MyTraj(const MagVec& Mo, const matrix& G,                     int N, double td);
MSVCDLL row_vector MyTraj(const MagVec& Mo, const matrix& G, const MagVec& Minf, int N, double td);

MSVCDLL row_vector MzTraj(const MagVec& Mo, const matrix& G,                     int N, double td);
MSVCDLL row_vector MzTraj(const MagVec& Mo, const matrix& G, const MagVec& Minf, int N, double td);

MSVCDLL row_vector MpTraj(const MagVec& Mo, const matrix& G,                     int N, double td);
MSVCDLL row_vector MpTraj(const MagVec& Mo, const matrix& G, const MagVec& Minf, int N, double td);

MSVCDLL row_vector MmTraj(const MagVec& Mo, const matrix& G,                     int N, double td);
MSVCDLL row_vector MmTraj(const MagVec& Mo, const matrix& G, const MagVec& Minf, int N, double td);

// ____________________________________________________________________________
// C                     M Magnetization Trajectories
// ____________________________________________________________________________

MSVCDLL coord_vec MTraj(const MagVec& Mo, const matrix& G, const MagVec& Minf,
                                                             int N, double td);

MSVCDLL row_vector MNormTraj(const MagVec& Mo, const matrix& G, const MagVec& Minf,
                                                             int N, double td);

#endif						// BlochTraj.h
