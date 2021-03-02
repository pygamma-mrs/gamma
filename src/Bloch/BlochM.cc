/* BlochM.cc ****************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**      Bloch Equation Specific Mangetization Vectors Implementation	**
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
**                                                                      **
** This file contains functions which generate magnetizaiton vectors    **
** for use in the phenomenological Bloch equations. In the simplest     **
** case the returned vector will be a 3x1 column vector which appears   **
** as                                                                   **
**                                                                      **
**                                  [ M  ]                              **
**                                  |  x |                              **
**                                  |    |                              **
**                            |M> = | M  |                              **
**                                  |  y |                              **
**                                  |    |                              **
**                                  | M  |                              **
**                                  [  z ]                              **
**									**
** In a more general context, the above vector will be a single block   **
** within a larger vector of dimension 3N where N is the number of      **
** sub-vectors in a general magnetization vector . In that case, the    **
** vector appears as                                                    **
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
**                                                                      **
** Since the GAMMA Bloch module contains the class MagVec to handle     **
** most aspects of such vectors, this file contains functions that      **
** produce some of the more commonly used/encountered ones.             **
**                                                                      **
*************************************************************************/

#ifndef   BlochM_cc_ 			// Is file already included?
#  define BlochM_cc_ 1      		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation   
#  endif

#include <Bloch/BlochM.h>		// Include our interface
#include <Bloch/BlochSys.h>		// Include Bloch systems
#include <Basics/Gutils.h>		// Include paramter query
#include <Basics/StringCut.h>		// Include function Gform

// ----------------------------------------------------------------------------
//                           Deprecated Functions
// ----------------------------------------------------------------------------

/* These are marked for removal in a future version. They will work but will
   issue warning messages.                                                   */

matrix Meq_vector()
  {
  std::string hdr("Bloch Module");
  GAMMAerror(hdr, 4, "R_matrix", 1);
  GAMMAerror(hdr, 6, "BlochMeq");
  matrix Meq(3, 1, complex0);
  Meq.put(1.0,2,0); 
  return Meq; 
  }

// sosi - this assumes a 3x3 array....
matrix Mss_vector(const matrix& G, const matrix&R, matrix& Meq)
  {
  std::string hdr("Bloch Module");
  GAMMAerror(hdr, 4, "Mss_vector", 1);
  GAMMAerror(hdr, 6, "BlochMss");
  matrix Mss(3, 1, complex0);
  if(R.getRe(0,0) || R.getRe(1,1)  		 // If gamB1=0 && R=0, Mss=0
  || (G.getRe(0,2) && G.getRe(1,2)))
    Mss = (inv(G)*R)*Meq;                       // Mss = Kinv*R*Meq
  return Mss;
  }

// ----------------------------------------------------------------------------
//                      Equilibrium Magnetization Vectors
// ----------------------------------------------------------------------------

/* These are just direct functions which return 3Nx1 column vectors (as MagVec)
   Them may contain multiple sub-vectors but they will always have all
   z-components set to 1 and all x&y components set to zero. That is, the
   equilibrium state has all vectors oriented along the +z axis.             */


MagVec BlochMeq()                    { return MagVec(0,0,1); }
MagVec BlochMeq(int NC)              { return MagVec(NC);    }
MagVec BlochMeq(const MagVec& M)     { return MagVec(M.NComps()); }
MagVec BlochMeq(const BlochSys& sys) { return sys.Mo(); }


// ----------------------------------------------------------------------------
//                      Mx, My, Mz Magnetization Vectors
// ----------------------------------------------------------------------------

/* These are just direct functions which return 3Nx1 column vectors (as MagVec)
   Them may contain multiple sub-vectors but they will always have all
   x-components set to 1 and all z&y components set to zero. That is, the
   state has all vectors oriented along the +x axis.                        */

MagVec BlochMX()                      { return MagVec(1,0,0);     }
MagVec BlochMX(int NC) 		      { return MagVec::MxVec(NC); }
MagVec BlochMX(const MagVec& M)       { return M.Mx();            }
MagVec BlochMX(const BlochSys& sys)   { return sys.Mx();          }

MagVec BlochMY()                      { return MagVec(0,1,0);     }
MagVec BlochMY(int NC) 		      { return MagVec::MyVec(NC); }
MagVec BlochMY(const MagVec& M)       { return M.My();            }
MagVec BlochMY(const BlochSys& sys)   { return sys.My();          }

MagVec BlochMZ()                      { return MagVec(0,0,1);     }
MagVec BlochMZ(int NC) 		      { return MagVec::MzVec(NC); }
MagVec BlochMZ(const MagVec& M)       { return M.Mz();            }
MagVec BlochMZ(const BlochSys& sys)   { return sys.Mz();          }

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

MagVec BlochMss(const matrix& G, const matrix& R, const MagVec& Meq)
  {
  if(R.is_zero()) return Meq;
  return (inv(G)*R)*Meq;
  }

// ----------------------------------------------------------------------------
//              Interactive Specification Of Magnetization Vectors
// ----------------------------------------------------------------------------

MagVec BlochMo(int argc, char* argv[], int& qn, double Mx, double My, double Mz)
  {
  std::string msg("\n\tMagnetization Vector ");
  std::string msgf = msg + "X Component [" + Gform("%8.2f", Mx) + "]? ";
  ask_set(argc, argv, qn++, msgf, Mx); 
  msgf = msg + "Y Component [" + Gform("%8.2f", My) + "]? ";
  ask_set(argc, argv, qn++, msgf, My); 
  msgf = msg + "Z Component [" + Gform("%8.2f", Mz) + "]? ";
  ask_set(argc, argv, qn, msgf, Mz); 
  return MagVec(Mx,My,Mz);
  }




#endif							// BlochM.cc

