/* SpinOpRot.cc *************************************************-*-c++-*-
**									**
**                               G A M M A				**
**									**
**      Spin Rotation Operators                    Implementation	**
**									**
**      Copyright (c) 1994						**
**									**
**      S.A. Smith							**
**      National High Magnetic Field Laboratory				**
**      1800 E. Paul Dirac Drive					**
**      Tallahassee Florida, 32306-4005					**
**									**
**      Z.L. Madi, T. Levante						**
**      Eidgenoessische Technische Hochschule				**
**      Labor fur physikalische Chemie					**
**      8092 Zurich / Switzerland					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
**  Description								**
**									**
** This module contains functions which produce spin rotation operators	**
** and rotated spin operators for a composite spin system.  These will	**
** span a spin Hilbert space, defined defined by a list of spins 	**
** specified by some spin system (class spin_sys).			**
**									**
**  The composite space matrix representation is in the natural basis	**
**  constructed from the direct product of single spin bases.		**
**									**
**  See class spin_sys for spin system function and structure.  	**
**  See module single_spin_op for single spin operators.       		**
**  See module single_spin_op for single spin rotation ops.    		**
**  See class spin_op (near death) for spin operator workings  		**
**  See class gen_op (future) if spin_op has been removed!     		**
**  See module comp_spin_op for functions of unrotated spin ops.	**
**									**
*************************************************************************/

#ifndef   SpinOpRot_cc_			// Is file already included?
#  define SpinOpRot_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <HSLib/SpinOpRot.h>		// Include the header file
#include <Basics/Gconstants.h>		// Include definition of PI
#include <Basics/Gutils.h>		// Include GAMMA error stuff
#include <HSLib/SpinOp.h>		// Knowledge of spin operators
#include <HSLib/SpinOpSng.h>		// Knowledge of single spin operators
#include <HSLib/SpinOpCmp.h>		// Knowledge of composite spin operators
#include <HSLib/SpinSys.h>		// Knowledge of base spin systems
#include <stdlib.h>
#include <string>			// Include libstdc++ strings

// ____________________________________________________________________________
// A             SPIN ROTATION OPERATORS ABOUT THE X,Y, OR Z AXIS
// ____________________________________________________________________________

/* The rotation operators returned here reside in the composite Hilbert space
   of the spin system given in the function call.  Overloads operate on a 
   single spin, all spins of a specified isotope type, all spins in the system,
   or any combination of spins in the system as specified by the spin flag 
   settings in the system.

  	   Input		sys  : Spin system (base)
	   			***  : Designation of affected spins
	   			beta : Rotation angle (degrees)
	   Return		SOp  : Spin rotation operator Ru for
 	  			       rotations about the u-axis
	  			       of angle beta.

	   OverLoads of Function where Type of Second Argument Changes

	   OL Fcn       Type    Label  Result
	   -- ---       ----    ----   ------------------------------
	   A. Ru  	int	spin : Spin index, only spin affected
	   B. Ru 	string 	iso  : Isotope label, type iso spins affected
	   C. Ru 	---	     : All spins in sys affected
	   D. Ru 	int*	flags: Spins with flags TRUE affected
	   E. Ru_sp	---	     : Spins with isys flags TRUE affected   */

// ---------------------- Rotations About The X-Axis --------------------------

spin_op Rx(const spin_sys& sys, int spin, double beta)
  { return Raxis(sys, spin, beta, 'x'); }

spin_op Rx(const spin_sys& sys, const std::string& iso, double beta)
  { return Raxis(sys, iso, beta, 'x'); }

spin_op Rx(const spin_sys& sys, double beta)
  { return Raxis(sys, beta, 'x'); }

spin_op Rx(const spin_sys& sys, const flagvec& flags, double beta)
  { return Raxis(sys, flags, beta, 'x'); }

spin_op Rx_sp(const spin_sys& sys, double beta)
  { return Raxis_sp(sys, beta, 'x'); }

// ---------------------- Rotations About The Y-Axis --------------------------

spin_op Ry(const spin_sys& sys, int spin, double beta)
  { return Raxis(sys, spin, beta, 'y'); }

spin_op Ry(const spin_sys& sys, const std::string& iso, double beta)
  { return Raxis(sys, iso, beta, 'y'); }

spin_op Ry(const spin_sys& sys, double beta)
  { return Raxis(sys, beta, 'y'); }

spin_op Ry(const spin_sys& sys, const flagvec& flags, double beta)
  { return Raxis(sys, flags, beta, 'y'); }

spin_op Ry_sp(const spin_sys& sys, double beta)
  { return Raxis_sp(sys, beta, 'y'); }

// ---------------------- Rotations About The Z-Axis --------------------------

spin_op Rz(const spin_sys& sys, int spin, double beta)
  { return Raxis(sys, spin, beta, 'z'); }

spin_op Rz(const spin_sys& sys, const std::string& iso, double beta)
  { return Raxis(sys, iso, beta, 'z'); }

spin_op Rz(const spin_sys& sys, double beta)
  { return Raxis(sys, beta, 'z'); }

spin_op Rz(const spin_sys& sys, const flagvec& flags, double beta)
  { return Raxis(sys, flags, beta, 'z'); }

spin_op Rz_sp(const spin_sys& sys, double beta)
  { return Raxis_sp(sys, beta, 'z'); }


// ____________________________________________________________________________
// B       GENERIC SPIN ROTATION OPERATORS ABOUT THE X,Y, OR Z AXIS
// ____________________________________________________________________________

/* The rotation operator functions in this section are generic, they handle
   construction of the rotation operators reqested in previously defined
   functions (about x,y, & z axes).  They are all similar and conserve code.

	   Input		sys  : Spin system (base)
	   			***  : Designation of affected spins
	   			beta : Rotation angle (degrees)
	   			axis : Rotation axis (x, y, or z)
	   Return		SOp  : Spin rotation operator Rx for
 	  			       rotations about the specified
	  			       Cartesian axis of angle beta.

	   OverLoads of Function where Type of Second Argument Changes

	   OL Fcn       Type    Label  Result
	   -- ---       ----    ----   ------------------------------
	   A. Raxis  	int	spin : Spin index, only spin affected
	   B. Raxis 	string	iso  : Isotope label, iso type spins affected
	   C. Raxis 	---	     : All spins in sys affected
	   D. Raxis 	int*	flags: Spins with flags TRUE affected
	   E. Raxis_sp	---	     : Spins with sys flags TRUE affected    */


spin_op Raxis(const spin_sys & sys, int spin,           double beta, char axis)
  { return Raxis(sys, sys.GetFlags(spin,1,0), beta, axis); }

spin_op Raxis(const spin_sys& sys, const std::string& iso,   double beta, char axis)
  { return Raxis(sys, sys.GetFlags(iso,1,0), beta, axis); }

spin_op Raxis(const spin_sys & sys,                     double beta, char axis)
  { return Raxis(sys, sys.GetFlags(1), beta, axis); }

spin_op Raxis(const spin_sys& sys,const flagvec& F,double beta, char axis)
  {
  int ns = sys.spins();				// Number of spins
  matrix* spr = new matrix[ns];			// Array of subspace matrices
  int hs;					// Single spin Hilbert space
  for(int spin=0; spin<ns; spin++)		// Loop spin in system
    {
    hs = sys.HS(spin);				//   This spin's Hilbert space
    if(F[spin]) spr[spin] = Raxis(hs,beta,axis);//   If active, add R from spin
    else        spr[spin] = Ie(hs);		//   If inactive, just I matrix
    } 
  spin_op SOp(ns, spr);				// Set SOp for system
  if(spr) delete [] spr;			// Delete our new arrays
  return SOp;
  }

spin_op Raxis_sp(const spin_sys& sys,                   double beta, char axis)
  { return Raxis(sys, sys.GetFlags(), beta, axis); }

// ____________________________________________________________________________
// C           ROTATION FUNCTIONS ABOUT AN AXIS IN THE XY-PLANE
// ____________________________________________________________________________

/* These rotation operators are not associated with one of the Cartesian axes.
   Rather, they are associated with one of the Cartesian planes in that they
   perform a rotation about any axis in a specific plane. Overloads operate on
   a single spin, all spins of a specified isotope type, all spins in the 
   system, or any combination of spins in the system as specified either by
   the spin flag settings in the system (see class spin_sys) or an integer
   array of spin flags in the function argument.

  	   Input		sys  : Spin system (base)
	   			***  : Designation of affected spins
	   			phi  : Phase angle (degrees)
	   			beta : Rotation angle (degrees)
	   Return		SOp  : Spin rotation operator Ruv for
 	  			       rotations about an axis in the
	  			       uv-plane phi degrees from u-axis
	  			       of angle beta (right hand rule rot).

	   OverLoads of Function where Type of Second Argument Changes

	   OL Fcn       Type    Label  Result
	   -- ---       ----    ----   ------------------------------
	   A. Ruv 	int	spin : Spin index, only spin affected
	   B. Ruv	string	iso  : Isotope label, type iso spins affected
	   C. Ruv	---	     : All spins in sys affected
	   D. Ruv	int*	flags: Spins with flags TRUE affected
	   E. Ruv_sp	---	     : Spins with sys flags TRUE affected    */

// ----------------- Rotations About An Axis In The XY-Plane ------------------

spin_op Rxy(const spin_sys& sys, int spin, double phi, double beta)
  { return Rplane(sys, spin, phi, beta, 'x'); }

spin_op Rxy(const spin_sys& sys, const std::string& iso, double phi, double beta)
  { return Rplane(sys, iso, phi, beta, 'x'); }

spin_op Rxy(const spin_sys& sys, double phi, double beta)
  { return Rplane(sys, phi, beta, 'x'); }

spin_op Rxy(const spin_sys& sys,const flagvec& flags,double phi,double B)
  { return Rplane(sys, flags, phi, B, 'x'); }

spin_op Rxy_sp(const spin_sys& sys, double phi, double beta)
{ return Rplane_sp(sys, phi, beta, 'x'); }


// ----------------- Rotations About An Axis In The YZ-Plane ------------------

spin_op Ryz(const spin_sys& sys, int spin, double phi, double beta)
  { return Rplane(sys, spin, phi, beta, 'y'); }

spin_op Ryz(const spin_sys& sys, const std::string& iso, double phi, double beta)
  { return Rplane(sys, iso, phi, beta, 'y'); }

spin_op Ryz(const spin_sys& sys, double phi, double beta)
  { return Rplane(sys, phi, beta, 'y'); }

spin_op Ryz(const spin_sys& sys,const flagvec& flags,double phi,double B)
  { return Rplane(sys, flags, phi, B, 'y'); }

spin_op Ryz_sp(const spin_sys& sys, double phi, double beta)
{ return Rplane_sp(sys, phi, beta, 'y'); }


// ----------------- Rotations About An Axis In The ZX-Plane ------------------

spin_op Rzx(const spin_sys& sys, int spin, double phi, double beta)
  { return Rplane(sys, spin, phi, beta, 'z'); }

spin_op Rzx(const spin_sys& sys, const std::string& iso, double phi, double beta)
  { return Rplane(sys, iso, phi, beta, 'z'); }

spin_op Rzx(const spin_sys& sys, double phi, double beta)
  { return Rplane(sys, phi, beta, 'z'); }

spin_op Rzx(const spin_sys& sys,const flagvec& flags,double phi,double B)
  { return Rplane(sys, flags, phi, B, 'z'); }

spin_op Rzx_sp(const spin_sys& sys, double phi, double beta)
  { return Rplane_sp(sys, phi, beta, 'z'); }


// -------------- Rotations About An Axis In A Specified Plane ----------------

/*         Input                sys  : Spin system (base)
                                ***  : Designation of affected spins
                                phi  : Phase angle (degrees)
                                beta : Rotation angle (degrees)
	   			plane: Rotation plane (xy, yz, or zx)
	   Return		SOp  : Spin rotation operator R for
 	  			       rotations of angle beta about an
	  			       axis in the plane specified (xy,yz,zx),
	  			       phi degrees from the main axis (x,y,z)

           OverLoads of Function where Type of Second Argument Changes

           OL Fcn       Type    Label  Result
           -- ---       ----    ----   ------------------------------
           A. Rplane    int     spin : Spin index, only spin affected
           B. Rplane    string	iso  : Isotope label, type iso spins affected
           C. Rplane    ---          : All spins in sys affected
           D. Rplane    vector  flags: Spins with flags TRUE affected
           E. Rplane_sp ---          : Spins with sys flags TRUE affected    */

spin_op Rplane(const spin_sys& sys, int spin, double phi, double beta, char p)
  { return Rplane(sys,sys.GetFlags(spin,1,0),phi,beta,p); }

spin_op Rplane(const spin_sys& sys, const std::string& iso,
                                               double phi, double beta, char p)
  { return Rplane(sys,sys.GetFlags(iso,1,0),phi,beta,p); }

spin_op Rplane(const spin_sys& sys, double phi, double beta, char plane)
  { return Rplane(sys, sys.GetFlags(1), phi, beta, plane); }


spin_op Rplane(const spin_sys & sys, const flagvec& flags, 
                                          double phi, double beta, char plane)

	// Input		sys  : Spin system (base)
	// 			phi  : Phase angle (degrees)
	// 			beta : Rotation angle (degrees)
	// 			plane: Rotation plane (xy, yz, or zx)
	// Return		SOp  : Spin rotation operator R for
 	//			       rotations of angle beta about an
	//			       axis in the plane specified phi
	//			       degrees from the main axis
	// Note			     : All spins whose spin flag is TRUE
	//			       are affected by R
	// sosi			     - In need of better sim.trans. routine

  {
  int dim;
  int ns = sys.spins();			// Number of spins in system
  matrix* spr = new matrix[ns];		// Array of subspace matrices
  double cb,sb,cp,sp;
  cb=cos(beta*PI/360.0);		// angle is beta/2 in radians
  sb=sin(beta*PI/360.0);		// angle is beta/2 in radians
  cp = cos(phi*PI/180.0);		// angle is phi in radians
  sp = sin(phi*PI/180.0);		// angle is phi in radians
  for(int spin=0; spin<ns; spin++)
    {
    dim = sys.HS(spin);
    if(flags[spin])
      {
      spr[spin]=matrix(dim,dim,complex0);
      matrix &p = spr[spin];
      switch (plane)
	{
	case 'x':			// axis in the xy plane
	  if(dim == 2)
	    {
	    p.put(cb,0,0);
	    p.put(sb*complex(-sp,-cp), 0, 1);
	    p.put(sb*complex(sp,-cp), 1, 0);
	    p.put(cb,1,1);
	    }
	  else
	    {
	    p = times_adjoint(Raxis(dim,beta,'x'),Raxis(dim,phi,'z'));
	    p = Raxis(dim,phi,'z')*p;
	    }
	  break;
	case 'y':			// axis in the yz plane
	  if(dim == 2)
	    {
	    p.put(complex(cb, -sb*cp), 0, 0);
	    p.put(-sb*sp, 0, 1);
	    p.put(sb*sp, 1, 0);
	    p.put(complex(cb, sb*cp), 1, 1);
	    }
	  else
	    {
	    p = Raxis(dim,beta,'z')*Raxis(dim,phi,'x');
	    p = adjoint_times(Raxis(dim,phi,'x'),p);
	    }
	  break;
	case 'z':			// axis in the zx plane
	default:
	  if (dim == 2)
	    {
	    p.put(complex(cb, -sb*cp), 0, 0);
	    p.put(complex(0.0, -sb*sp), 0, 1);
	    p.put(complex(0.0, -sb*sp), 1, 0);
	    p.put(complex(cb, sb*cp), 1, 1);
	    }
	  else
	    {
	    p = times_adjoint(Raxis(dim,beta,'z'),Raxis(dim,phi,'y'));
	    p = Raxis(dim,phi,'y')*p;
	    }
	  break;
	}
      }
    else spr[spin]=Ie(dim);
    }
  spin_op SOp(ns, spr);			// Set SOp for system
  if(spr) delete [] spr;
  return SOp;
  }

spin_op Rplane_sp(const spin_sys& sys, double phi, double beta, char plane)
  { return Rplane(sys, sys.GetFlags(), phi, beta, plane); }


// ____________________________________________________________________________
// D              ROTATION FUNCTIONS ABOUT AN ARBITRARY AXIS
// ____________________________________________________________________________

/* These rotation operators are ment to perform rotations about any specified
   axis. In this instance the rotation axis is set from the two spherical
   angles. Overloads operate on a single spin, all spins of a specified isotope
   type, all spins in the system, or any combination of spins in the system as
   specified by the spin flag settings in the system (see class spin_sys).   */

// --------------- Rotations About An Axis In Spherical Space -----------------

/*         Input                sys  : Spin system (base)
                                ***  : Designation of affected spins
                                theta: Rotation axis polar angle (degrees)
                                phi  : Rotation axis polar angle (degrees)
                                beta : Rotation angle (degrees)
           Return               SOp  : Spin rotation operator Rxyz for
                                       rotations of angle beta about an axis
                                       specified by theta and phi, affecting
                                       the spin(s) specified

           OverLoads of Function where Type of Second Argument Changes

      OL Fcn       Type    Label  Result
      -- ---       ----    ----   ------------------------------
      A. Rxyz	int     spin : Spin index, only spin affected
      B. Rxyz	string	iso  : Isotope label, spins of type iso affected
      C. Rxyz	---          : All spins in sys affected
      D. Rxyz	vector  flags: Spins with flags TRUE affected
      E. Rxyz_sp	---  : Spins with internal sys flags TRUE affected   */


spin_op Rxyz(const spin_sys& sys,int spin,double theta,double phi,double beta)
  { return Rspace(sys,sys.GetFlags(spin,1,0),theta,phi,beta); }

spin_op Rxyz(const spin_sys& sys, const std::string& iso, 
	                                 double theta, double phi, double beta)
  { return Rspace(sys,sys.GetFlags(iso,1,0),theta,phi,beta); }

spin_op Rxyz(const spin_sys& sys, double theta, double phi, double beta)
  { return Rspace(sys,sys.GetFlags(1),theta,phi,beta); }

spin_op Rxyz(const spin_sys& sys,const flagvec& FGs,
                                          double theta, double phi,double beta)
  { return Rspace(sys,FGs,theta,phi,beta); }

spin_op Rxyz_sp(const spin_sys& sys, double theta, double phi, double beta)
  { return Rspace(sys,sys.GetFlags(),theta,phi,beta); }


// ----------- Generic Rotations About An Axis In Spherical Space -------------

spin_op Rspace(const spin_sys & sys, const flagvec& flags, 
                                        double theta, double phi, double beta)

        // Input                sys  : Spin system (base)
	//			flags: Array of spin flags
        //                      theta: Rotation axis polar angle
        //                      phi  : Rotation axis polar angle
        //                      beta : Rotation angle
        // Return               SOp  : Spin rotation operator Rxyz for
        //                             rotations of angle beta about an axis
        //                             specified by theta and phi
        // Note                      : Spins whose spin flags are TRUE
        //                             in flags are affected by R
	// sosi			     - In need of better sim.trans. routine

  {
  int dim;
  int ns = sys.spins();
  double cb,sb,cp,sp,ct,st;
  matrix* spr = new matrix [ns];
  cb = cos(beta*PI/360.0);		// angle is beta/2 in radians
  sb = sin(beta*PI/360.0);		// angle is beta/2 in radians
  cp = cos(phi*PI/180.0);		// angle is phi in radians
  sp = sin(phi*PI/180.0);		// angle is phi in radians
  ct = cos(theta*PI/180.0);		// angle is theta in radians
  st = sin(theta*PI/180.0);		// angle is theta in radians
  for(int spin=0; spin<ns; spin++)	// Loop over all spins in sys
    {
    dim = sys.HS(spin);			// Current spin Hilbert space
    if(flags[spin])			// See if rotation affects spin
      {
      spr[spin]=matrix(dim,dim,complex0);
      matrix &p = spr[spin];
      if(dim == 2)			// For I=1/2, do this (easy)
	{
	 p.put(complex(cb, -ct*sb), 0,0);	// cos(B) - i*cos(T)*sin(B)
	 p.put(sb*complex(-st*sp,-st*cp), 0,1);
	 p.put(sb*complex(st*sp,-st*cp), 1,0);
	 p.put(complex(cb, ct*sb), 1,1);
	 }
       else				// For I!=1/2, do this (harder)
	 {
	 p = times_adjoint(Raxis(dim,beta,'z'),Raxis(dim,theta,'y'));
	 p = Raxis(dim,theta,'y')*p;	// p = Ryz(theta,beta) here
	 p = times_adjoint(p,Raxis(dim,phi,'z'));
	 p = Raxis(dim,phi,'z')*p;	// p = Rxyz (theta,phi,beta)
	 }
      }
    else spr[spin] = Ie(dim);		// Unaffected spin contribution
    }
  spin_op SOp(ns, spr);			// Set SOp for system
  if(spr) delete [] spr;
  return SOp;
  }

// ____________________________________________________________________________
// E                      EULER ANGLE ROTATION FUNCTIONS
// ____________________________________________________________________________

/* These rotation operators are meant to perform rotations about any specified
   axis. In this instance the rotation axis is set from the three Euler angles.
   Overloads operate on a single spin, all spins of a specified isotope type,
   all spins in the system, or any combination of spins in the system as
   specified by the spin flag settings in the system (see class spin_sys).   */

// ---------- Euler Rotations (Three Rotations, Three Euler Angles) -----------

/*         Input		sys  : Spin system (base)
	   			alpha: Euler angle alpha (or a)
	   			beta : Euler angle beta  (or b)
	   			gamma: Euler angle gamma (or g)
	   Return		SOp  : Spin rotation operator for Euler
 	  			       rotation of angles alpha, beta,
	  			       and gamma with specified selectivity

           OverLoads of Function where Type of Second Argument Changes

  OL Function   Type    Label  Result
  -- ---------- ------- -----  ------------------------------
  A. R_Euler    int     spin : Spin index, only spin affected
  B. R_Euler    string& iso  : Isotope label, spins of type iso affected
  C. R_Euler    ---          : All spins in sys affected
  D. R_Euler_sp ---          : Spins w/ internal sys flags TRUE affected     */

spin_op R_Euler(const spin_sys& sys, int spin, double a, double b, double g)
  { return R_Euler_plane(sys, sys.GetFlags(spin,1,0), a, b, g); }

spin_op R_Euler(const spin_sys& sys,const std::string& I,double a,double b,double g)
  { return R_Euler_plane(sys, sys.GetFlags(I,1,0), a, b, g); }

spin_op R_Euler(const spin_sys& sys, double a, double b, double g)
  { return R_Euler_plane(sys, sys.GetFlags(1), a, b, g); }

spin_op R_Euler_sp(const spin_sys& sys, double a, double b, double g)
  { return R_Euler_plane(sys, sys.GetFlags(), a, b, g); }

// --------------------- Generic Euler Rotation Function ----------------------

// sosi - This function remains untested !!

spin_op R_Euler_plane(const spin_sys& sys, const flagvec& flags, 
                                                  double a, double b, double g)
  {
  int dim;
  int ns = sys.spins();
  double q,cb,sb,cs,ss,cd,sd;
  matrix* spr = new matrix[ns];
  cb = cos(b/2.0);
  sb = sin(b/2.0);
  cs = cos((a+g)/2.0);
  ss = sin((a+g)/2.0);
  cd = cos((a-g)/2.0);
  sd = sin((a-g)/2.0);
  for(int j=0; j<ns; j++)
    {
    q = sys.qn(j);
    if(q != 0.5)
      {
      GAMMAerror("Spin Operator", "Euler Rotation On Spin I>1/2", 1);
      GAMMAerror("Spin Operator", "Problems in R_Euler_plane");
      exit(-1);
      }
    dim = sys.HS(j);	// get size of spin Hilbert space
    if(flags[j] == 1)
      {
      spr[j]=matrix(dim,dim,complex(0.0));
      matrix &p = spr[j];
      p.put(cb*complex(cs, ss), 0, 0);
      p.put(-sb*complex(cd, sd), 0, 1);
      p.put(sb*complex(cd,-sd), 1, 0);
      p.put(cb*complex(cs,-ss), 1, 1);
      }
    else
      spr[j]=Ie(dim);
    }
  spin_op SOp(ns, spr);				// Set SOp for system
  if(spr) delete [] spr;
  return SOp;
  }


// ____________________________________________________________________________
// F                           ROTATED SPIN OPERATORS
// ____________________________________________________________________________

/* These functions produce common spin operators which have been rotated. Most
   of the un-rotated spin operators are found in the GAMMA module SpinOpCmp.
   Overloads operate on a single spin, all spins of a specified isotope type,
   all spins in the system, or any combination of spins in the system as 
   specified by the spin flag settings in the system (see class spin_sys).   */

// sosi - most of these functions can use sim_trans when spin operators
//        start using class matrix more soundly & saving CPU time
// sosi - MAKE SURE Rz returns identity when theta is zero!!
// sosi - these suck because they copy the spin flags 3x each!!!
//        need to make spin flags then call the more specific functions

// ------------ Spin Operators Rotated About the Cartesian Z-Axis -------------

/*	   Input		sys  : Spin system (base)
	   			***  : Designation of affected spins
	   			theta: Rotation angle (degrees)
	   Return		SOp  : Spin operator Fxy which is Fx
 	  			       rotated about the z-axis of angle theta

	   OverLoads of Function where Type of Second Argument Changes

	OL Fcn     Type   Label  Result
	-- ---     ----   -----   ------------------------------
	A. Ixy     int    spin : Spin index, only spin affected
	A. Fxy 	   int	  spin : Spin index, only spin affected
	B. Fxy 	   string iso  : Isotope label, spins of type iso affected
	C. Fxy 	   ---	       : All spins in sys affected
	D. Fxy 	   vector flags: Spins with flags TRUE affected
	E. Fxy_sp  ---         : Spins with internal sys flags TRUE affected */


spin_op Ixy(const spin_sys& sys, int S, double theta)
  { return Rz(sys,S,theta) * Faxis(sys,S,'x') * adjoint(Rz(sys,S,theta)); }

spin_op Fxy(const spin_sys& sys, int S, double theta)
  { return Rz(sys,S,theta) * Faxis(sys,S,'x') * adjoint(Rz(sys,S,theta)); }

spin_op Fxy(const spin_sys& sys, const std::string& iso, double theta)
  { return Rz(sys,iso,theta)*Faxis(sys,iso,'x')*adjoint(Rz(sys,iso,theta)); }

spin_op Fxy(const spin_sys& sys, double theta)
  { return Rz(sys,theta) * Faxis(sys,'x') * adjoint(Rz(sys,theta)); }

spin_op Fxy(const spin_sys& sys, const flagvec& FGs, double theta)
  { return Rz(sys,FGs,theta)*Faxis(sys,FGs,'x')*adjoint(Rz(sys,FGs,theta)); }

spin_op Fxy_sp(const spin_sys& sys, double theta)
  { return Rz_sp(sys,theta) * Faxis_sp(sys,'x') * adjoint(Rz_sp(sys,theta)); }


// ----------------------------------------------------------------------------
//                      Spin Operators Associated With F+
// ----------------------------------------------------------------------------

/*	   OverLoads of Function where Type of Second Argument Changes

	   OL Fcn       Type    Label  Result
	   -- ---       ----    ----   ------------------------------
	   A. Ip  	int	spin : Spin index, only spin affected
	   B. Fp 	int	spin : Spin index, only spin affected
	   C. Fp 	string	iso  : Isotope label, type iso spins affected
	   D. Fp 	---	     : All spins in sys affected
	   E. Fp 	int*	flags: Spins with flags TRUE affected        */

spin_op Ip(const spin_sys& sys, int S, double theta)
  { return Rz(sys,S,theta) * Faxis(sys,S,'p') * Rz(sys,S,theta).adjoint(); }

spin_op Fp(const spin_sys& sys, int S, double theta)
  { return Rz(sys,S,theta) * Faxis(sys,S,'p') * Rz(sys,S,theta).adjoint(); }

spin_op Fp(const spin_sys& sys, const std::string& iso, double theta)
  { return Rz(sys,iso,theta)*Faxis(sys,iso,'p') * Rz(sys,iso,theta).adjoint(); }

spin_op Fp(const spin_sys& sys, double theta)
  { return Rz(sys,theta) * Faxis(sys,'p') * Rz(sys,theta).adjoint(); }

spin_op Fp(const spin_sys& sys, const flagvec& FGs, double theta)
  { return Rz(sys,FGs,theta)*Faxis(sys,FGs,'p') * Rz(sys,FGs,theta).adjoint(); }

spin_op Fp_sp(const spin_sys& sys, double theta)
  { return Rz_sp(sys,theta) * Faxis_sp(sys,'p') * Rz_sp(sys,theta).adjoint(); }


// ----------------------------------------------------------------------------
//                      Spin Operators Associated With F-
// ----------------------------------------------------------------------------

/*	   OverLoads of Function where Type of Second Argument Changes

	   OL Fcn       Type    Label  Result
	   -- ---       ----    ----   ------------------------------
	   A. Im  	int	spin : Spin index, only spin affected
	   B. Fm 	int	spin : Spin index, only spin affected
	   C. Fm 	string	iso  : Isotope label, type iso spins affected
	   D. Fm 	---	     : All spins in sys affected
	   E. Fm 	int*	flags: Spins with flags TRUE affected        */

spin_op Im(const spin_sys& sys, int S, double theta)
  { return RotSpinOp(Rz(sys,S,theta),Faxis(sys,S,'m')); }

spin_op Fm(const spin_sys& sys, int S, double theta)
  { return RotSpinOp(Rz(sys,S,theta),Faxis(sys,S,'m')); }

spin_op Fm(const spin_sys& sys, const std::string& iso, double theta)
  { return RotSpinOp(Rz(sys,iso,theta),Faxis(sys,iso,'m')); }

spin_op Fm(const spin_sys& sys, double theta)
  { return RotSpinOp(Rz(sys,theta),Faxis(sys,'m')); }

spin_op Fm(const spin_sys& sys, const flagvec& FGs, double theta)
  { return RotSpinOp(Rz(sys,FGs,theta),Faxis(sys,FGs,'m')); }

spin_op Fm_sp(const spin_sys & sys, double theta)
  { return RotSpinOp(Rz_sp(sys,theta),Faxis_sp(sys,'m')); }


// ____________________________________________________________________________
// G                    GENERIC ROTATED SPIN OPERATOR FUNCTION
// ____________________________________________________________________________


spin_op Fplane(const spin_sys& sys, double theta, char OPtype)

	// Input		sys	: Spin system (base)
	// 			OPtype	: Operator type
	// Return		SOp	: Spin operator rotated about
	//			          the z-axis by the angle theta
	// Note				: Spins included if flags = TRUE
	// sosi				: Should be able to do this better

  { 
  spin_op R = Rz_sp(sys,theta);
  spin_op F = Faxis(sys, OPtype);
  return RotSpinOp(R, F);
  }


spin_op RotSpinOp(const spin_op& R, const spin_op& F) {return R*F*R.adjoint();}

	// Input		R	: A spin rotation operator
	// 			F	: A spin operator to be rotated
	// Return		SOp	: F rotated by R
	// sosi				: Should be able to do this better

#endif 						// SpinOpRot.cc

