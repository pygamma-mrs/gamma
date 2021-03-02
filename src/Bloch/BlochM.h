/* BlochM.h *****************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**      Bloch Equation Specific Magnetization Vectors     Interface	**
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
** This file contains functions which generate magnetizaiton vectors 	**
** for use in the phenomenological Bloch equations. In the simplest	**
** case the returned vector will be a 3x1 column vector which appears	**
** as 									**
**                                                                      **
**                                  [ M  ]				**
**                                  |  x |                              **
**                                  |    |                              **
**                            |M> = | M  |                              **
**                                  |  y |                              **
**                                  |    |                              **
**                                  | M  |                              **
**                                  [  z ]                              **
**                                                  			**
** In a more general context, the above vector will be a single block	**
** within a larger vector of dimension 3N where N is the number of 	**
** sub-vectors in a general magnetization vector . In that case, the	**
** vector appears as							**
**                                                                      **
**                                  [  [M ]  ]                          **
**                                  |  [ 0]  |                          **
**                                  |        |                          **
**                                  |  [M ]  |                          **
**                                  |  [ 1]  |                          **
**                                  |    .   |                          **
**                            |M> = |    .   |                          **
**                                  |    .   |                          **
**                                  |    .   |                          **
**                                  |    .   |                          **
**                                  | [M   ] |                          **
**                                  [ [ N-1] ]                          **
**                                                  			**
** Since the GAMMA Bloch module contains the class MagVec to handle 	**
** most aspects of such vectors, this file contains functions that	**
** produce some of the more commonly used/encountered ones.		**
**                                                                      **
*************************************************************************/

#ifndef   BlochM_h_ 			// Is file already included?
#  define BlochM_h_ 1      		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface   
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/matrix.h>		// Know about GAMMA matrices
#include <Bloch/BlochSys.h>		// Know about GAMMA Bloch systems
#include <Bloch/MagVec.h>		// Know about magnetization vectors
#include <string>			// Know about libstdc++ strings

// ----------------------------------------------------------------------------
//                           Deprecated Functions
// ----------------------------------------------------------------------------

/* These are marked for removal in a future version. They will work but will
   issue warning messages.                                                   */

MSVCDLL matrix Meq_vector();
MSVCDLL matrix Mss_vector(const matrix& G, const matrix&R, matrix& Meq);

// ----------------------------------------------------------------------------
//                      Equilibrium Magnetization Vectors
// ----------------------------------------------------------------------------

/* These are just direct functions which return 3Nx1 column vectors (as MagVec)
   Them may contain multiple sub-vectors but they will always have all
   z-components set to 1 and all x&y components set to zero. That is, the
   equilibrium state has all vectors oriented along the +z axis.             */

MSVCDLL MagVec BlochMeq();
MSVCDLL MagVec BlochMeq(int NC);
MSVCDLL MagVec BlochMeq(const MagVec& M);
MSVCDLL MagVec BlochMeq(const BlochSys& sys);

// ----------------------------------------------------------------------------
//                      Mx, My, Mz Magnetization Vectors
// ----------------------------------------------------------------------------

/* These are just direct functions which return 3Nx1 column vectors (as MagVec)
   Them may contain multiple sub-vectors but they will always have all
   x-components set to 1 and all z&y components set to zero. That is, the
   state has all vectors oriented along the +x axis.                        */

MSVCDLL MagVec BlochMX();
MSVCDLL MagVec BlochMX(int NC);
MSVCDLL MagVec BlochMX(const MagVec& M);
MSVCDLL MagVec BlochMX(const BlochSys& sys);

MSVCDLL MagVec BlochMY();
MSVCDLL MagVec BlochMY(int NC);
MSVCDLL MagVec BlochMY(const MagVec& M);
MSVCDLL MagVec BlochMY(const BlochSys& sys);

MSVCDLL MagVec BlochMZ();
MSVCDLL MagVec BlochMZ(int NC);
MSVCDLL MagVec BlochMZ(const MagVec& M);
MSVCDLL MagVec BlochMZ(const BlochSys& sys);


// ----------------------------------------------------------------------------
//                    Steady-State Magnetization Vectors
// ----------------------------------------------------------------------------

/* These are just direct functions which return 3Nx1 column vectors (as MagVec)
   Them may contain multiple sub-vectors. They are formally generated by the
   formula
                                         R
                               |M  > =   - |M  >
                                 ss      G   eq

   but of course this will have troubles whenever the matrix R is zero or the
   matrix G cannot be easily inverted. When R is zero (no-relaxation) then
   there simply is no steady-state defined. Alternativly, if R=G (i.e. there
   are no rf-fields appled) the the steady state is one and the same as the
   equilibrium state.                                                        */

MSVCDLL MagVec BlochMss(const matrix& G, const matrix& R, const MagVec& Meq);


// ----------------------------------------------------------------------------
//              Interactive Specification Of Magnetization Vectors
// ----------------------------------------------------------------------------

MSVCDLL MagVec BlochMo(int argc, char* argv[], int& qn, double Mx=0, double My=0, double Mz=1);

#endif								// BlochM.h
