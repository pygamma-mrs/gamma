/* SpinOpRot.h **************************************************-*-c++-*-
**									**
**                                  G A M M A				**
**									**
**      Rotation Spin Operators                     Interface		**
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
**                                                                      **
**  Description                                                         **
**                                                                      **
** This module contains functions which produce spin rotation operators **
** and rotated spin operators for a composite spin system.  These will  **
** span a spin Hilbert space, defined defined by a list of spins        **
** specified by some spin system (class spin_sys).                      **
**                                                                      **
**  The composite space matrix representation is in the natural basis   **
**  constructed from the direct product of single spin bases.           **
**                                                                      **
**  See class spin_sys for spin system function and structure.          **
**  See module single_spin_op for single spin operators.                **
**  See module single_spin_op for single spin rotation ops.             **
**  See class spin_op (near death) for spin operator workings           **
**  See class gen_op (future) if spin_op has been removed!              **
**  See module comp_spin_op for functions of unrotated spin ops.        **
**                                                                      **
*************************************************************************/

#ifndef   SpinOpRot_h_			// Is file already included?
#  define SpinOpRot_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/SpinOp.h>		// Include spin operators
#include <HSLib/SpinSys.h>		// Include base spin systems
#include <string>			// Include libstdc++ strings
#include <vector>			// Include libstdc++ stl vectors

/* The functions defined in this section return common spin rotation
   operators as well as some common "rotated" spin operators.  The code
   which produces them all is found in the GAMMA module SpinOpRot.cc.
   This header file, SpinOpRot.h, is currently not doing anything at all.
   Current plans are to remove class spin_op from GAMMA (although it will
   still work) - letting class matrix handle all of its internal attributes.
   When that occurs, this header file will be activated and the section in
   SpinOp.h that defines these functions now will be deleted.                */


// ____________________________________________________________________________
// A             SPIN ROTATION OPERATORS ABOUT THE X,Y, OR Z AXIS
// ____________________________________________________________________________
 
/* The rotation operators returned here reside in the composite Hilbert space
   of the spin system given in the function call.  Overloads operate on a
   single spin, all spins of a specified isotope type, all spins in the system,
   or any combination of spins in the system as specified by the spin flag
   settings in the system                                                    */
 
/*         Input                sys  : Spin system (base)
                                ***  : Designation of affected spins
                                beta : Rotation angle (degrees)
           Return               SOp  : Spin rotation operator Ru for
                                       rotations about the u-axis
                                       of angle beta where u={x,y,z}.
 
           OverLoads of Function where Type of Second Argument Changes
 
           OL Fcn       Type    Label  Result
           -- ---       ----    ----   ------------------------------
           A. Ru        int     spin : Spin index, only spin affected
           B. Ru        string  iso  : Isotope label, type iso spins affected
           C. Ru        ---          : All spins in sys affected
           D. Ru        vector  flags: Spins with flags TRUE affected
           E. Ru_sp     ---          : Spins with isys flags TRUE affected   */

// ---------------------- Rotations About The X-Axis --------------------------

MSVCDLL spin_op Rx(const spin_sys& sys,    int spin,                double beta);
MSVCDLL spin_op Rx(const spin_sys& sys,    const std::string& iso,  double beta);
MSVCDLL spin_op Rx(const spin_sys& sys,                             double beta);
MSVCDLL spin_op Rx(const spin_sys& sys,    const flagvec& flags,    double beta);
MSVCDLL spin_op Rx_sp(const spin_sys& sys,                          double beta);

// ---------------------- Rotations About The Y-Axis --------------------------

MSVCDLL spin_op Ry(const spin_sys& sys,    int spin,               double beta);
MSVCDLL spin_op Ry(const spin_sys& sys,    const std::string& iso, double beta);
MSVCDLL spin_op Ry(const spin_sys& sys,                            double beta);
MSVCDLL spin_op Ry(const spin_sys& sys,    const flagvec& flags,   double beta);
MSVCDLL spin_op Ry_sp(const spin_sys& sys,                         double beta);

// ---------------------- Rotations About The Z-Axis --------------------------

MSVCDLL spin_op Rz(const spin_sys& sys,    int spin,                  double beta);
MSVCDLL spin_op Rz(const spin_sys& sys,    const std::string& iso,         double beta);
MSVCDLL spin_op Rz(const spin_sys& sys,                               double beta);
MSVCDLL spin_op Rz(const spin_sys& sys,    const flagvec& flags, double beta);
MSVCDLL spin_op Rz_sp(const spin_sys& sys,                            double beta);
 
// ____________________________________________________________________________
// B       GENERIC SPIN ROTATION OPERATORS ABOUT THE X,Y, OR Z AXIS
// ____________________________________________________________________________
 
/* The rotation operator functions in this section are generic, they handle
   construction of the rotation operators reqested in previously defined
   functions (about x,y, & z axes).  They are all similar and conserve code.
 
           Input                sys  : Spin system (base)
                                ***  : Designation of affected spins
                                beta : Rotation angle (degrees)
                                axis : Rotation axis (x, y, or z)
           Return               SOp  : Spin rotation operator Rx for
                                       rotations about the specified
                                       Cartesian axis of angle beta.
 
           OverLoads of Function where Type of Second Argument Changes
 
           OL Fcn       Type    Label  Result
           -- ---       ----    ----   ------------------------------
           A. Raxis     int     spin : Spin index, only spin affected
           B. Raxis     string  iso  : Isotope label, iso type spins affected
           C. Raxis     ---          : All spins in sys affected
           D. Raxis     vector  flags: Spins with flags TRUE affected
           E. Raxis_sp  ---          : Spins with sys flags TRUE affected    */

MSVCDLL spin_op Raxis(const spin_sys& sys,int I,             double beta, char axis);
MSVCDLL spin_op Raxis(const spin_sys& sys,const std::string& iso, double beta, char axis);
MSVCDLL spin_op Raxis(const spin_sys& sys,                   double beta, char axis);
MSVCDLL spin_op Raxis(const spin_sys& sys,const flagvec& flags,double B,char A);
MSVCDLL spin_op Raxis_sp(const spin_sys& sys,                double beta, char axis);


// ____________________________________________________________________________
// C           ROTATION FUNCTIONS ABOUT AN AXIS IN THE XY-PLANE
// ____________________________________________________________________________

// Note - These rotation operators are not associated with one of the
//	  Cartesian axes.  Rather, they are associated with one of the
//	  Cartesian planes in that they perform a rotation about any
//	  axis in a specific plane. Overloads operate on a single spin,
//        all spins of a specified isotope type, all spins in the system,
//        or any combination of spins in the system as specified by the
//	  spin flag settings in the system (see class spin_sys).

// -------------- Rotations About An Axis In The XY-Plane ---------------

        // Input                sys  : Spin system (base)
        //                      ***  : Designation of affected spins
        //                      phi  : Phase angle (degrees)
        //                      beta : Rotation angle (degrees)
        // Return               SOp  : Spin rotation operator Rxy for
        //                             rotations about an axis in the
        //                             xy-plane phi degrees from x-axis
        //                             of angle beta.

        // OverLoads of Function where Type of Second Argument Changes

        // OL Fcn       Type    Label  Result
        // -- ---       ----    ----   ------------------------------
        // A. Rxy       int     spin : Spin index, only spin affected
        // B. Rxy       string	iso  : Isotope label, spins of type iso affected
        // C. Rxy       ---          : All spins in sys affected
        // D. Rxy       vector  flags: Spins with flags TRUE affected
        // E. Rxy_sp    ---          : Spins w/ internal sys flags TRUE affected

MSVCDLL spin_op Rxy(const spin_sys& sys,int spin,              double phi,double beta);
MSVCDLL spin_op Rxy(const spin_sys& sys,const std::string& iso,     double phi,double beta);
MSVCDLL spin_op Rxy(const spin_sys& sys,                       double phi,double beta);
MSVCDLL spin_op Rxy(const spin_sys& sys,const flagvec& flags,double phi,double B);
MSVCDLL spin_op Rxy_sp(const spin_sys& sys,                    double phi,double beta);


// ----------------- Rotations About An Axis In The YZ-Plane ------------------

        // Input                sys  : Spin system (base)
        //                      ***  : Designation of affected spins
        //                      phi  : Phase angle (degrees)
        //                      beta : Rotation angle (degrees)
	// Return		SOp  : Spin rotation operator Ryz for
 	//			       rotations about an axis in the
	//			       yz-plane theta degrees over from z-axis
        //                             of angle beta.

        // OverLoads of Function where Type of Second Argument Changes

        // OL Fcn       Type    Label  Result
        // -- ---       ----    ----   ------------------------------
        // A. Ryz       int     spin : Spin index, only spin affected
        // B. Ryz       string	iso  : Isotope label, type iso spins affected
        // C. Ryz       ---          : All spins in sys affected
        // D. Ryz       vector  flags: Spins with flags TRUE affected
        // E. Ryz_sp    ---          : Spins w/ sys flags TRUE affected

MSVCDLL spin_op Ryz(const spin_sys& sys, int spin,         double theta,double beta);
MSVCDLL spin_op Ryz(const spin_sys& sys, const std::string& iso,double theta,double beta);
MSVCDLL spin_op Ryz(const spin_sys& sys,                   double theta,double beta);
MSVCDLL spin_op Ryz(const spin_sys& sys, const flagvec& flags, double T,double B);
MSVCDLL spin_op Ryz_sp(const spin_sys& sys,                double theta,double beta);


// -------------- Rotations About An Axis In The ZX-Plane ---------------

        // Input                sys  : Spin system (base)
        //                      ***  : Designation of affected spins
        //                      phi  : Phase angle (degrees)
        //                      beta : Rotation angle (degrees)
	// Return		SOp  : Spin rotation operator Rzx for
 	//			       rotations about an axis in the
	//			       zx-plane theta degrees down from z-axis
        //                             of angle beta.

        // OverLoads of Function where Type of Second Argument Changes

        // OL Fcn       Type    Label  Result
        // -- ---       ----    ----   ------------------------------
        // A. Rzx       int     spin : Spin index, only spin affected
        // B. Rzx       string	iso  : Isotope label, type iso spins affected
        // C. Rzx       ---          : All spins in sys affected
        // D. Rzx       vector  flags: Spins with flags TRUE affected
        // E. Rzx_sp    ---          : Spins with sys flags TRUE affected

MSVCDLL spin_op Rzx(const spin_sys& sys, int spin,          double theta, double beta);
MSVCDLL spin_op Rzx(const spin_sys& sys, const std::string& iso, double theta, double beta);
MSVCDLL spin_op Rzx(const spin_sys& sys,                    double theta, double beta);
MSVCDLL spin_op Rzx(const spin_sys& sys, const flagvec& flags, double P, double B);
MSVCDLL spin_op Rzx_sp(const spin_sys& sys,                 double phi, double beta);


// ----------- Rotations About An Axis In A Specified Plane -------------

/*         Input                sys  : Spin system (base)
                                ***  : Designation of affected spins
                                phi  : Phase angle (degrees)
                                beta : Rotation angle (degrees)
                                plane: Rotation plane (xy, yz, or zx)
           Return               SOp  : Spin rotation operator R for
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

MSVCDLL spin_op Rplane(const spin_sys& sys,int spin, double phi, double beta, char p);
MSVCDLL spin_op Rplane(const spin_sys& sys,const std::string& iso,
                                             double phi, double beta, char p);
MSVCDLL spin_op Rplane(const spin_sys& sys,          double phi, double beta, char p);
MSVCDLL spin_op Rplane(const spin_sys& S,const flagvec& F,
                                             double phi, double beta, char p);
MSVCDLL spin_op Rplane_sp(const spin_sys &sys,       double phi, double beta, char p);

// ____________________________________________________________________________
// D           ROTATION FUNCTIONS ABOUT AN ARBITRARY AXIS
// ____________________________________________________________________________

/* These rotation operators are ment to perform rotations about any specified
   axis. In this instance the rotation axis is set from the two spherical
   angles. Overloads operate on a single spin, all spins of a specified isotope
   type, all spins in the system, or any combination of spins in the system
   as specified by the spin flag settings in the system (see class spin_sys) */

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
      A. Rxyz      int     spin : Spin index, only spin affected
      B. Rxyz      string& iso  : Isotope label, spins of type iso affected
      C. Rxyz      ---          : All spins in sys affected
      D. Rxyz      vector  flags: Spins with flags TRUE affected
      E. Rxyz_sp   ---          : Spins w/ internal sys flags TRUE affected  */

MSVCDLL spin_op Rxyz(const spin_sys& sys,int spin,double theta,double phi,double beta);
MSVCDLL spin_op Rxyz(const spin_sys& sys,const std::string& iso,
                                                  double theta,double phi,double beta);
MSVCDLL spin_op Rxyz(const spin_sys& sys,         double theta,double phi,double beta);
MSVCDLL spin_op Rxyz(const spin_sys& sys,const flagvec& flags,
                                                  double theta,double phi,double beta);
MSVCDLL spin_op Rxyz_sp(const spin_sys& sys,      double theta,double phi,double beta);

// ----------- Generic Rotations About An Axis In Spherical Space -------------

MSVCDLL spin_op Rspace(const spin_sys & sys, const flagvec& flags,
                                        double theta, double phi, double beta);

        // Input                sys  : Spin system (base)
        //                      flags: Array of spin flags
        //                      theta: Rotation axis polar angle
        //                      phi  : Rotation axis polar angle
        //                      beta : Rotation angle
        // Return               SOp  : Spin rotation operator Rxyz for
        //                             rotations of angle beta about an axis
        //                             specified by theta and phi
        // Note                      : Spins whose spin flags are TRUE
        //                             in flags are affected by R

// ____________________________________________________________________________
// E                     EULER ANGLE ROTATION FUNCTIONS
// ____________________________________________________________________________

/* These rotation operators are ment to perform rotations about any specified
   axis. In this instance the rotation axis is set from the three Euler angles.
   Overloads operate on a single spin, all spins of a specified isotope type,
   all spins in the system, or any combination of spins in the system as
   specified by the spin flag settings in the system (see class spin_sys).   */

// ---------- Euler Rotations (Three Rotations, Three Euler Angles) -----------

/*         Input		sys  : Spin system (base)
	   			alpha: Euler angle alpha
	   			beta : Euler angle beta
	   			gamma: Euler angle gamma
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

MSVCDLL spin_op R_Euler(const spin_sys &sys, int spin, double a, double b, double g);
MSVCDLL spin_op R_Euler(const spin_sys &sys,const std::string& iso,
                                               double alpha,double beta,double gamma);
MSVCDLL spin_op R_Euler(const spin_sys &sys,   double alpha,double beta,double gamma);
MSVCDLL spin_op R_Euler_sp(const spin_sys &sys,double alpha,double beta,double gamma);

// --------------------- Generic Euler Rotation Function ----------------------

MSVCDLL spin_op R_Euler_plane(const spin_sys &sys, const flagvec& flags,
                                     double alpha, double beta, double gamma); 

// ____________________________________________________________________________
// F                      ROTATED SPIN OPERATORS
// ____________________________________________________________________________

/* Here are functions to produce common spin operators which have been rotated.
   Most of the un-rotated spin operators are found in the GAMMA module 
   SpinOpCmp. Overloads operate on a single spin, all spins of a specified 
   isotope type, all spins in the system, or any combination of spins in the 
   system as specified by the spin flag settings in the system (see class
   spin_sys).                                                                */

// ------------ Spin Operators Rotated About the Cartesian Z-Axis -------------

/*         Input                sys  : Spin system (base)
                                ***  : Designation of affected spins
                                theta: Rotation angle (degrees)
           Return               SOp  : Spin operator Fxy which is Fx
                                       rotated about the z-axis
                                       of angle theta

           OverLoads of Function where Type of Second Argument Changes

     OL Fcn       Type    Label  Result
     -- ---       ----    ----   ------------------------------
     A. Ixy       int     spin : Spin index, only spin affected
     A. Fxy       int     spin : Spin index, only spin affected
     B. Fxy       strin   iso  : Isotope label, spins of type iso affected
     C. Fxy       ---          : All spins in sys affected
     D. Fxy       vector  flags: Spins with flags TRUE affected
     E. Fxy_sp    ---          : Spins with internal sys flags TRUE affected */

MSVCDLL spin_op Ixy(const spin_sys &sys, int spin,                  double theta);
MSVCDLL spin_op Fxy(const spin_sys &sys, int spin,                  double theta);
MSVCDLL spin_op Fxy(const spin_sys &sys, const std::string& iso,    double theta);
MSVCDLL spin_op Fxy(const spin_sys &sys,                            double theta);
MSVCDLL spin_op Fxy(const spin_sys &sys, const flagvec& flags,      double theta);
MSVCDLL spin_op Fxy_sp(const spin_sys &sys,                         double theta);

// ----------------------------------------------------------------------------
//                      Spin Operators Associated With F+
// ----------------------------------------------------------------------------

//   The Function Declarations Are in the Header of Class spin_op (SpinOp.h)

/*         OverLoads of Function where Type of Second Argument Changes
 
           OL Fcn       Type    Label  Result
           -- ---       ----    ----   ------------------------------
           A. Ip        int     spin : Spin index, only spin affected
           B. Fp        int     spin : Spin index, only spin affected
           C. Fp        string  iso  : Isotope label, type iso spins affected
           D. Fp        ---          : All spins in sys affected
           E. Fp        vector  flags: Spins with flags TRUE affected        */

MSVCDLL spin_op Ip(const spin_sys &sys, int spin,                  double theta);
MSVCDLL spin_op Fp(const spin_sys &sys, int spin,                  double theta);
MSVCDLL spin_op Fp(const spin_sys &sys, const std::string& iso,         double theta);
MSVCDLL spin_op Fp(const spin_sys &sys, const flagvec& flags, double theta);
MSVCDLL spin_op Fp(const spin_sys &sys,                            double theta);
MSVCDLL spin_op Fp_sp(const spin_sys &sys,                         double theta);

// ----------------------------------------------------------------------------
//                      Spin Operators Associated With F-
// ----------------------------------------------------------------------------

// The Function Declarations Are in the Header of Class spin_op (SpinOp.h)

/*         OverLoads of Function where Type of Second Argument Changes

           OL Fcn       Type    Label  Result
           -- ---       ----    ----   ------------------------------
           A. Im        int     spin : Spin index, only spin affected
           B. Fm        int     spin : Spin index, only spin affected
           C. Fm        string  iso  : Isotope label, type iso spins affected
           D. Fm        ---          : All spins in sys affected
           E. Fm        vector  flags: Spins with flags TRUE affected        */

MSVCDLL spin_op Im(const spin_sys &sys, int spin,                  double theta);
MSVCDLL spin_op Fm(const spin_sys &sys, int spin,                  double theta);
MSVCDLL spin_op Fm(const spin_sys &sys, const std::string& iso,    double theta);
MSVCDLL spin_op Fp(const spin_sys &sys, const flagvec& flags,      double theta);
MSVCDLL spin_op Fm(const spin_sys &sys,                            double theta);
MSVCDLL spin_op Fm_sp(const spin_sys &sys,                         double theta);

// ____________________________________________________________________________
// G                 GENERIC ROTATED SPIN OPERATOR FUNCTION
// ____________________________________________________________________________


MSVCDLL spin_op Fplane(const spin_sys&sys, double theta, char OPtype);

	// Input		sys    : Spin system (base)
	// 			OPtype : Operator type
	// Return		SOp    : Spin operator rotated about
	//			         the z-axis by the angle theta
	// Note		               : Spins included if flags = TRUE

MSVCDLL spin_op RotSpinOp(const spin_op& R, const spin_op& F);

	// Input		R	: A spin rotation operator
	// 			F       : An operator to be rotated
	// Return		SOp	: F rotateb by R
 
#endif					         // SpinOpRot.h
