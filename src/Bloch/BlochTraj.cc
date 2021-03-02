/* BlochTraj.cc *************************************************-*-c++-*-
**									**
**                             G A M M A				**
**									**
**      Bloch Magnetization Trajectories	   Implementation	**
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
**                                                                      **
** This functions provided herein generate magnetization trajectories   **
** over evenly spaced time increments. The users specifies the initial  **
** magnetization vector, the evolution matrix, and the infinite time    **
** magnetization vector. The functions will generate with a row_vector  **
** of magnetization values or else a coord_vec of magnetization         **
** components.                                                          **
**                                                                      **
** The generalized solution to the Bloch equation(s) that's implemented **
** herein is                                                            **
**                      -Gt                                             **
**          |M(t)> = exp   [ |M(0)> - |M(inf)> ] + |M(inf)>             **
**                                                                      **
** where G is the Bloch evolution matrix, |M(0)> the initial vector and **
** |M(inf)> the vector given infinite evolution time. The functions     **
** take a particular detection vector <D| and generate a trajectory for **
** <D| according to                                                     **
**                            -n*t                                      **
**  D(nt ) = <D|M(t)> = <D|exp    d [ |M(0)> - |M(inf)> ] + <D|M(inf)>  **
**      d                                                               **
**                                                                      **
** The supplied trajectory functions typically demand that the user     **
** supply the intial and infinite time magnetization vectors, the       **
** evolution matrix G, the number of trajectory points n, and t         **
** the time increment between points. The detection vector D   d        **
** is usually implicit in the function call.                            **
**									**
*************************************************************************/

#ifndef   BlochTraj_cc_ 		// Is file already included?
#  define BlochTraj_cc_ 1      		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation   
#  endif

#if defined(_MSC_VER)			// If using MSVC++ then we
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif 

#include <Bloch/BlochTraj.h>		// Include our interface
#include <Basics/Gutils.h>		// Include paramter query
#include <Level1/coord.h>       // Include GAMMA coordinates

// ____________________________________________________________________________
// A                       Bloch Detection Vectors
// ____________________________________________________________________________

DetVec MxDetection(const MagVec& M) { return DetVec(M.NComps(), 'x'); }
DetVec MyDetection(const MagVec& M) { return DetVec(M.NComps(), 'y'); }
DetVec MzDetection(const MagVec& M) { return DetVec(M.NComps(), 'z'); }
DetVec MpDetection(const MagVec& M) { return DetVec(M.NComps(), '+'); }
DetVec MmDetection(const MagVec& M) { return DetVec(M.NComps(), '-'); }

// ____________________________________________________________________________
//                       U-Magnetization Trajectories
// ____________________________________________________________________________

/* These functions will evolve an input magnetization vector |M> under the
   Bloch matrix G for N time increments of td. At each point, magnetization is
   stored either in a row vector or coordinate vector. The magnetization 
   detected may be user selected with the vector D or function specific for
   common usage. If the infinite time matrix Minf is supplied then it is the
   difference magnetization vector that is evolved.

     Function Name                                  Purpose
     -------------  -----------------------------------------------------------
         MTraj      Row vector of magnetization using detection vector D
         MxTraj     Row vector of Mx magnetization
         MyTraj     Row vector of My magnetization
         MzTraj     Row vector of Mz magnetization

   For Mx the we track <3*N|M(t)>, for My <3*N+1|M(t)>, and for Mz <3*N+2|M(t)>
   These correspond to use of detection vectors <Dx|, <Dy|, & <Dz| with elements

            <D |i> = d          <D |i> = d             <D |i> = d
              x       i,3*N       y       i,3*N+1        z       i,3*N+2

   where d is a Kronecker delta. Two other common detection vectors generated
   complex magnetization data. These are M+ (Mx+iMy) & M- (Mx-iMy) with elements

            <D |i> = d      + i*d          <D |i> = d      - i*d
              +       i,3*N      i,3*N+1     -       i,3*N      i,3*N+1       */


// ----------------------------------------------------------------------------
//                        General M Component Trajectories
// ----------------------------------------------------------------------------

row_vector MTraj(const DetVec& D, const MagVec& Mo,
                                             const matrix& G, int N, double td)
  {
  matrix GD, S;					// G eigensystem
  G.Diagonalize(GD, S);				// Diagonalize G
  matrix     Sinv = inv(S);			// Get inverse of S
  DetVec Det  = row_vector(D)*S;		// <D| --> <D*S| 
// The following line makes no sense. Changed M on RHS to Mo.
// I am not sure anybody is using the Bloch stuff at all.
// MAER 11/2014
  col_vector M    = Sinv*Mo;			// |M> --> |inv(S)*M>
  int         bd   = Mo.size();			// Bloch dimension
  std::vector<int> OK;				// Array of indices 
  int i,j;					// Temp indeices
  for(i=0; i<bd; i++)				// Loop Bloch dimension
    {						// and see what contributes
    if(norm(Det.get(i)*M.get(i) > 1.e-9))	// to the detected values
      OK.push_back(i);
    }
  int rd = OK.size();				// Reduced dimension
  if(!rd) return row_vector(N,complex0);	// Zero if no contributions
  std::vector<complex> I(rd);			// Array of intensities
  std::vector<complex> B(rd);			// Array of w's & R's
  std::vector<complex> delB(rd);			// Time incrementation
  for(i=0; i<rd; i++)				// Loop over contributors
    {						// & store in reduced form
    j       = OK[i];
    I[i]    = Det.get(j)*M.get(j);
    B[i]    = complex1;
    delB[i] = exp(-td*GD.get(j,j));
    }
  row_vector traj(N);				// Trajectory of M
  complex z;
  for(i=0; i<N; i++)				// Loop desired points
    {
    z = complex0;
    for(j=0; j<rd; j++)
      { z += I[j]*B[j]; B[j]*=delB[j]; }
    traj.put(z, i);
    }
  return traj;
  }

row_vector MTraj(const DetVec& D, const MagVec& Mo,
                         const matrix& G, const MagVec& Minf, int N, double td)
  {
  matrix GD, S;					// G eigensystem
  G.Diagonalize(GD, S);				// Diagonalize G
  matrix      Sinv = inv(S);			// Get inverse of S
  DetVec Det  = row_vector(D)*S;		// <D| --> <D*S| 
  MagVec delM = Sinv*(Mo-Minf);			// |delM> --> |inv(S)*delM>
  int         bd   = Mo.size();			// Bloch dimension
  std::vector<int> OK;				// Array of indices 
  int i,j;					// Temp indeices
  for(i=0; i<bd; i++)				// Loop Bloch dimension
    {						// and see what contributes
    if(norm(Det.get(i)*delM.get(i) > 1.e-9))	// to the detected values
      OK.push_back(i);
    }
  int rd = OK.size();				// Reduced dimension
  complex zinf = row_vector(D)*Minf;		// Infinite contribution 
  if(!rd) return row_vector(N,zinf);		// Zero if no contributions
  std::vector<complex> I(rd);			// Array of intensities
  std::vector<complex> B(rd);			// Array of w's & R's
  std::vector<complex> delB(rd);			// Time incrementation
  for(i=0; i<rd; i++)				// Loop over contributors
    {						// & store in reduced form
    j       = OK[i];
    I[i]    = Det.get(j)*delM.get(j);
    B[i]    = complex1;
    delB[i] = exp(-td*GD.get(j,j));
    }
  row_vector traj(N);				// Trajectory of Mx
  complex z;
  for(i=0; i<N; i++)				// Loop desired points
    {
    z = complex0;
    for(j=0; j<rd; j++)
      { z += I[j]*B[j]; B[j]*=delB[j]; }
    traj.put(z+zinf, i);
    }
  return traj;
  }

// ----------------------------------------------------------------------------
//                           Mx Component Trajectories
// ----------------------------------------------------------------------------

row_vector MxTraj(const MagVec& Mo, const matrix& G, int N, double td)
  {
  DetVec  D = MxDetection(Mo);			// Detection vector
  return MTraj(D, Mo, G, N, td);			// Use overload
  }

row_vector MxTraj(const MagVec& Mo, const matrix& G, const MagVec& Minf,
                                                             int N, double td)
  {
  DetVec  D = MxDetection(Mo);			// Detection vector
  return MTraj(D, Mo, G, Minf, N, td);			// Use overload
  }

// ----------------------------------------------------------------------------
//                         My Component Trajectories
// ----------------------------------------------------------------------------

row_vector MyTraj(const MagVec& Mo, const matrix& G,          int N, double td)
  {
  DetVec  D = MyDetection(Mo);			// Detection vector
  return MTraj(D, Mo, G, N, td);			// Use overload
  }

row_vector MyTraj(const MagVec& Mo, const matrix& G, const MagVec& Minf,
                                                              int N, double td)
  {
  DetVec  D = MyDetection(Mo);			// Detection vector
  return MTraj(D, Mo, G, Minf, N, td);			// Use overload
  }

// ----------------------------------------------------------------------------
//                         Mz Component Trajectories
// ----------------------------------------------------------------------------

row_vector MzTraj(const MagVec& Mo, const matrix& G,          int N, double td)
  {
  DetVec  D = MzDetection(Mo);			// Detection vector
  return MTraj(D, Mo, G, N, td);			// Use overload
  }
  
row_vector MzTraj(const MagVec& Mo, const matrix& G, const MagVec& Minf,
                                                              int N, double td)
  {
  DetVec  D = MzDetection(Mo);			// Detection vector
  return MTraj(D, Mo, G, Minf, N, td);			// Use overload
  }

// ----------------------------------------------------------------------------
//                         Mx + iMy Component Trajectories
// ----------------------------------------------------------------------------

row_vector MpTraj(const MagVec& Mo, const matrix& G,          int N, double td)
  {
  DetVec  D = MpDetection(Mo);			// Detection vector
  return MTraj(D, Mo, G, N, td);			// Use overload
  }
  
row_vector MpTraj(const MagVec& Mo, const matrix& G, const MagVec& Minf,
                                                              int N, double td)
  {
  DetVec  D = MpDetection(Mo);			// Detection vector
  return MTraj(D, Mo, G, Minf, N, td);			// Use overload
  }

// ----------------------------------------------------------------------------
//                         Mx - iMy Component Trajectories
// ----------------------------------------------------------------------------

row_vector MmTraj(const MagVec& Mo, const matrix& G,          int N, double td)
  {
  DetVec  D = MmDetection(Mo);			// Detection vector
  return MTraj(D, Mo, G, N, td);			// Use overload
  }
  
row_vector MmTraj(const MagVec& Mo, const matrix& G, const MagVec& Minf,
                                                              int N, double td)
  {
  DetVec  D = MmDetection(Mo);			// Detection vector
  return MTraj(D, Mo, G, Minf, N, td);			// Use overload
  }

// ____________________________________________________________________________
// C                     M Magnetization Trajectories
// ____________________________________________________________________________

coord_vec MTraj(const MagVec& Mo, const matrix& G, const MagVec& Minf,
                                                             int N, double td)
  {
  DetVec  Dx   = MxDetection(Mo);		// Mx Detection vector
  DetVec  Dy   = MyDetection(Mo);		// My Detection vector
  DetVec  Dz   = MzDetection(Mo);		// Mz Detection vector
  MagVec     delM = MagVec(Mo - Minf);		// Difference vector
  double xinf = Re(row_vector(Dx)*Minf);	// Infinite contribution 
  double yinf = Re(row_vector(Dy)*Minf);	// Infinite contribution 
  double zinf = Re(row_vector(Dz)*Minf);	// Infinite contribution 
  coord Mxyzinf(xinf,yinf,zinf);
  matrix GD, S;					// G eigensystem					
  G.Diagonalize(GD, S);				// Diagonalize G
  matrix Sinv = inv(S);				// Get inverse of S
  Dx          = row_vector(Dx)*S;		// <Dx| --> <Dx*S| 
  Dy          = row_vector(Dy)*S;		// <Dy| --> <Dy*S| 
  Dz          = row_vector(Dz)*S;		// <Dz| --> <Dz*S| 
  delM        = Sinv*delM;			// |delM> --> |inv(S)*delM>
  coord_vec traj(N);				// Trajectory of M
  coord pt;
  complex x,y,z;
// sosi this should be sped up
//int        NC   = Mo.NComps();		// Number of components
matrix tmp1;
  for(int i=0; i<N; i++)
    {
    tmp1 = -td*GD*i;
    x = row_vector(Dx)*tmp1.exp()*delM;
    y = row_vector(Dy)*tmp1.exp()*delM;
    z = row_vector(Dz)*tmp1.exp()*delM;
    pt = coord(Re(x), Re(y), Re(z));
    traj.put(pt+Mxyzinf, i);
    }
  return traj;
  }


row_vector MNormTraj(const MagVec& Mo, const matrix& G, const MagVec& Minf,
                                                             int N, double td)
  {
  DetVec  Dx   = MxDetection(Mo);		// Mx Detection vector
  DetVec  Dy   = MyDetection(Mo);		// My Detection vector
  DetVec  Dz   = MzDetection(Mo);		// Mz Detection vector
  MagVec     delM = MagVec(Mo - Minf);		// Difference vector
  double xinf = Re(row_vector(Dx)*Minf);	// Infinite contribution 
  double yinf = Re(row_vector(Dy)*Minf);	// Infinite contribution 
  double zinf = Re(row_vector(Dz)*Minf);	// Infinite contribution 
  coord Mxyzinf(xinf,yinf,zinf);
  matrix GD, S;					// G eigensystem					
  G.Diagonalize(GD, S);				// Diagonalize G
  matrix Sinv = inv(S);				// Get inverse of S
  Dx          = row_vector(Dx)*S;				// <Dx| --> <Dx*S| 
  Dy          = row_vector(Dy)*S;				// <Dy| --> <Dy*S| 
  Dz          = row_vector(Dz)*S;				// <Dz| --> <Dz*S| 
  delM        = Sinv*delM;			// |delM> --> |inv(S)*delM>
  row_vector NormVec(N);				// Trajectory of M
  double pt;
  complex x,y,z;
matrix tmp1;
  for(int i=0; i<N; i++)
    {
    tmp1 = -td*GD*i;
    x = row_vector(Dx)*tmp1.exp()*delM + xinf;
    y = row_vector(Dy)*tmp1.exp()*delM + yinf;
    z = row_vector(Dz)*tmp1.exp()*delM + zinf;
    pt = sqrt(Re(x*x) + Re(y*y) + Re(z*z));
    NormVec.put(pt, i);
    }
  return NormVec;
  }

#endif						// BlochTraj.cc

