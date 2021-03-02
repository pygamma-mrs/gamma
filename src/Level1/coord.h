/* coord.h ******************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Coordinate                                Interface		**
**                                                                      **
**      Copyright (c) 1990, 1991, 1992                                  **
**      Scott Smith                                                     **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
** Class coordinate provides the  means with which a coordinate in      **
** 3-dimensional space can be easily manipulated.  Each coordinate      **
** consists of three real numbers which may be translated, rotated,     **
** etc. Functions are provided for component access, distance and       **
** angle computations, and various other entities.                      **
**                                                                      **
*************************************************************************/

///Chapter Class Coord
///Section Overview
///Body    The class Coord is
///Section Available Coordinate Functions

#ifndef   Coord_h_			// Is file already included?
#  define Coord_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <Matrix/matrix.h>		// Provides knowledge of matrices
#include <Basics/SinglePar.h>		// Provides GAMMA parameters
#include <Basics/ParamSet.h>		// Provides GAMMA parameter sets
#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <string>			// Know libstdc++ strings

MSVCDLL matrix Rmx(double alpha, double beta, double gamma);

class coord
{
    
// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------

  double cx;				// First ordinate
  double cy;				// Second ordinate
  double cz;				// Third ordinate
					// Output format flags
  static int         PTotype;  		//   Output type
  static int         PTscience; 	//   Scientific notation
  static int         PTdigits;		//   Digits per point
  static int         PTprecise;		//   Digits after decimal
  static std::string PTform;		//   Output format
  static coord       DefCoord;		// Default coordinate
  static double      OrdCutoff;		// Zero cutoff for ordinate

/*                      Individual Coordinate Output Format
                        -----------------------------------

        ptotpe  : 0    == output as x,y,z       Cartesian (default)
	          1    == output as r,theta,phi Spherical 
	          2    == output as r,theta,z   Cylindrical
        science: !0    == output as 5.0e5
                  0    == output as 500000.00
        digits :          number of total digits output
        precise:          digits after decimal point                         */
  
private:

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS COORDINATE ERROR HANDLING
// ____________________________________________________________________________

/* These functions handle all error messages for class coord.  They tie into
   the generic GAMMA error messaging format.

              Input             pt      : A coordinate
                                eidx    : Error index
                                pn	: string in message
                                nr	: Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */ 

         void PTerror(int eidx,                        int nr=0) const;
         void PTerror(int eidx, const std::string& pn, int nr=0) const;
volatile void PTfatal(int eidx)                                  const;
volatile void PTfatal(int eidx, const std::string& pn)           const;

// ____________________________________________________________________________
// ii                    CLASS COORDINATE SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions can translate parameters in a given GAMMA parameter set into
   a coordinate. These are priviate functions because often misuse of these
   types of functions can mess up the class structure. I do not see how that
   is possible for this class, but it is better to be safe. They are not
   needed by the user anyway.                                                */


        // Input                 pt     : A coordinate point (this)
        //                       pset   : A parameter set
        //                       indx   : Point index
        //                       warn   : Warning output level
        // Output (*Cartesian)   TF     : Coordinate point filled with
        //                                Cartesian coordinate in pset
        // Output (*Spherical)   TF     : Coordinate point filled with
        //                                Spherical coordinate in pset
        // Output (*Cylindrical) TF     : Coordinate point filled with
        //                                Cylindrical coordinate in pset
        // Note                         : Regardless of the input type, the
        //                                points will always be stored in
        //                                Cartesian space. Conversion is done
        //                                as needed.

bool SetPtCartesian(const   ParameterSet& pset, int indx, int warn=1);
bool SetPtSpherical(const   ParameterSet& pset, int indx, int warn=1);
bool SetPtCylindrical(const ParameterSet& pset, int indx, int warn=1);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A              CLASS GENERAL COORDINATE CONSTRUCTORS/DESTRUCTOR
// ____________________________________________________________________________

///Center Constructors, Destructors, and Asssignment
///F_list coord		     - Coordinate point constructor
///F_list =                  - Coordinate point assignment
	
/* These constructors set up a new coordinate.  There are several ways to make
   a coordinate as listed below:
 
          Input Arguments                       Resulting Coordinate
   --------------------------        -----------------------------------------
                -                    Zero coordinate, {0,0,0}
   double x,double y,double z        Coordinate {x,y,z}, y & z 0 by default
             coord pt                Coordinate { pt.x, pt.y, pt.z }
     ParameterSet, idx, warn         Coordinate with index i read from pset 
             SinglePar               Coordinate set from parameter

   Setting a coordinate from a parameter set may fail if the parameters that
   define the coordinate don't exist. The flag warn dictates what happens
   upon a failure: 0=no warnings, 1=non-fatal warnings, 2=fatal warnings
   Setting a coordinate from a single parameter assumes a parameter of type
   3 with value stored as a string "( #, #, # )". All I/O between GAMMA
   coordinates and single parameters must follow this format!              */

MSVCDLC        coord( );
MSVCDLC        coord(double xx, double yy=0, double zz=0);
MSVCDLC        coord(const coord& pt1);
MSVCDLC        coord(const ParameterSet& pset, int idx=-1, int warn=2);
MSVCDLC        coord(const SinglePar& par);
MSVCDLC        ~coord( );
MSVCDLL coord& operator= (const coord& pt);

// ____________________________________________________________________________
// B                CARTESIAN BASED COORDINATE ACCESS FUNCTIONS
// ____________________________________________________________________________

///Center Coordinate & Ordinate Access Functions
///F_list get		     - Ordinate access
///F_list x		     - First ordinate access
///F_list y		     - Second ordinate access
///F_list z		     - Third ordinate access
///F_list xyz		     - Set all three ordinates

/*      Function                                     Result
    -----------------             ---------------------------------------------
       get(int)                   Return ordinate specified by int
       x,x(double)                Get/Set first ordinate
       y,y(double)                Get/Set second ordinate
       z,z(double)                Get/Set third ordinate
       xyz(...)                   Set all ordinates                          */

MSVCDLL double get(int i)    const;
MSVCDLL double x()           const;
MSVCDLL void   x(double xx);
MSVCDLL double y()           const;
MSVCDLL void   y(double yy);
MSVCDLL double z()           const;
MSVCDLL void   z(double zz);
MSVCDLL void   xyz(double xx, double yy, double zz);
MSVCDLL void   xyz(const coord& pt1);
MSVCDLL double norm() const;

// ____________________________________________________________________________
// C                SPHERICAL BASED COORDINATE ACCESS FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                    Rad is The Distance From The Origin
// ----------------------------------------------------------------------------

///F_list phi		     - Access spherical angle phi
///F_list Rad		     - Access radius
///F_list theta              - Access spherical angle theta


	// Input		pt   : Coordinate (this)
	// Return		Rad  : Radius or length from (0,0,0)
 
        // Input                pt      : Coordinate (this)
        //                      pt2     : A second coordinate
        // Return               Rad     : Radius or length vector from
        //                                between the two points
 
	// Input		x,y,z: Coordinate
	// Return		Rad  : Radius or length from (0,0,0)
 
MSVCDLL double Rad() const;
MSVCDLL double Rad(const coord& pt2) const;
MSVCDLL friend double Rad(double x, double y, double z);

// ----------------------------------------------------------------------------
//               Theta is The Angle Down From The Positive Z Axis
// ----------------------------------------------------------------------------
 
        // Input                pt    : Coordinate (this)
        // Return               theta : Spherical angle, down from +z axis
        // Note                       : Two cases exist -
        //                              1.) R=0  --> theta=0
        //                              2.) R!=0 --> 0<=theta<=180
        //                              The latter case range coincides with
        //                              the principal values returned by the
        //                              arcosine function

MSVCDLL double theta() const;
MSVCDLL double theta(const coord& pt2) const;
 
        // Input                pt    : Coordinate (this)
        //                      pt2   : A second coordinate
        // Return               theta : Spherical angle, down from +z axis
        //                              of the vector connecting pt to pt2
        // Note                       : Two cases exist -
        //                              1.) R=0  --> theta=0
        //                              2.) R!=0 --> 0<=theta<=180
        //                              The latter case range coincides with
        //                              the principal values returned by the

MSVCDLL friend double theta(double x, double y, double z);

	// Input		x,y,z: Coordinate
	// Return		theta: Spherical angle, down from the
	//			       z axis
 
// ----------------------------------------------------------------------------
//               Phi is The Angle Over From The Positive X Axis
// ----------------------------------------------------------------------------

MSVCDLL double phi()                 const;
MSVCDLL double phi(const coord& pt2) const;
 
	// Input		pt   : Coordinate (this)
	// Return		phi  : Spherical longitudinal angle
        // Input                pt      : Coordinate (this)
        //                      pt2     : A second coordinate
        // Return               phi     : Spherical angle, over from the
        //                                x axis for the vector connecting
        //                                pt & pt2

MSVCDLL friend double phi(double x, double y, double z);

	// Input		x,y,z: Coordinate
	// Return		phi  : Spherical angle, over from the
	//			       x axis
	// Note			     : Several cases exist -
	//					   x-axis
	//				1.) y=0, x>=0  --> phi=0
	//				2.) y=0, x<0   --> phi=180
	//					left hemisphere
	//				3.) y<0, x=0   --> phi=270
	//				4.) y<0, x<0   --> 180<=phi<=270
	//				5.) y<0, x>0   --> 270<=phi<=360
	//					right hemisphere
	//				6.) y>0, x=0   --> phi=90
	//				7.) y>0, x<0   --> 90<=phi<=180
	//				8.) y>0, x>0   --> 0<=phi<=90

	//			       Cases 4,5, & 7 lie outside the range
	//			       covered by the arctangent function.
	//			       The principal value range of atan is
	//			       [-90, 90]


MSVCDLL void invert();

	// Input		pt   : Coordinate (this)
	// Return		none : All three coordinate points are
	//			       inverted
        ///F_list phi		     - Access spherical angle phi



// ____________________________________________________________________________
// D                     COORDINATE ROTATION FUNCTIONS
// ____________________________________________________________________________

///Center Coordinate Rotation Functions

// ----------------------------------------------------------------------------
//            Typical Rotations Active On Cartesian Coordinates
// ----------------------------------------------------------------------------

MSVCDLL static matrix Rz(double phi, int rad=0);
 
        // Input                pt   : Coordinate (this)
        //                      phi  : Spherical angle (degrees over from +x)
	//			rad  : Flag if angle in radians or degrees(0)
        // Return               mx   : Rotation matrix to rotate point
        //                             by angle theta about the z-axis
	// Note			     : This rotates the coordinate via
	//			       the right hand rule about +z
	// Note			     : The reference axes are static!


MSVCDLL static matrix Rx(double theta, int rad=0);

        // Input                pt   : Coordinate (this)
        //                      theta: Spherical angle (degrees down from +z)
        //                      rad  : Flag for angle in degrees or radians
        // Return               mx   : Rotation matrix to rotate point
        //                             by angle phi about the x-axis
	// Note			     : This rotates the coordinate via
	//			       the right hand rule about +x
	// Note			     : The reference axes are static!

 
MSVCDLL static matrix Ry(double theta, int rad=0);
 
        // Input                pt   : Coordinate (this)
        //                      theta: Spherical angle (degrees down from +z)
        //                      rad  : Flag for angle in degrees or radians
        // Return               mx   : Rotation matrix to rotate point
        //                             by angle theta about the y-axis
 
 
MSVCDLL coord xrotate(double theta, int rad=0) const;

        // Input                pt   : Coordinate (this)
        //                      theta: Spherical angle (degrees down from +z)
        //                      rad  : Flag for angle in degrees or radians
        // Return               rotpt: The pt rotated about the x-axis by
        //                             the input angle
        // Note                      : The reference axes are static!  Thus
        //                             multiple rotations should be cumulative
 
 
MSVCDLL coord yrotate(double theta, int rad=0) const;

        // Input                pt   : Coordinate (this)
        //                      theta: Spherical angle (degrees down from +z)
        //                      rad  : Flag for angle in degrees or radians
        // Return               rotpt: The pt rotated about the y-axis by
        //                             the input angle
        // Note                      : The reference axes are static!  Thus
        //                             multiple rotations should be cumulative
 
 
MSVCDLL coord zrotate(double phi, int rad=0) const;

        // Input                pt   : Coordinate (this)
        //                      phi  : Spherical angle (over from +x)
        //                      rad  : Flag for angle in degrees or radians
        // Return               rotpt: The pt rotated about the z-axis by
        //                             the input angle
        // Note                      : The reference axes are static!  Thus
        //                             multiple rotations should be cumulative


// ----------------------------------------------------------------------------
//           Euler Rotation Matrices Active On Cartesian Coordinates
// ----------------------------------------------------------------------------
 
/* For these functions we assume that the AXES ROTATE COUNTER-CLOCKWISE
   & that the COORDINATES ARE STATIC.  This is normal for Euler rotations.   */

        // Input                pt   : Coordinate (this)
        //                      alpha: Euler angle in radians
        //                      beta : Euler angle in radians
        //                      gamma: Euler angle in radians
        //                      rad  : Flag for angle in degrees or radians
        // Return  (alpha)      mx   : Euler rotation matrix to axes by angle
        //                             alpha counter-clockwise about z
        // Return  (beta)       mx   : Euler rotation matrix to axes by angle
        //                             beta counter-clockwise about y
        // Return  (gamma)      mx   : Euler rotation matrix to axes by angle
        //                             gamma counter-clockwise about z

MSVCDLL matrix Ralpha(double alpha, int rad) const; 
MSVCDLL matrix Rbeta(double  beta,  int rad) const;
MSVCDLL matrix Rgamma(double gamma, int rad) const;

MSVCDLL matrix REuler(double alpha, double beta, double gamma, int rad) const;
 
        // Input                pt   : Coordinate (this)
        //                      alpha: Euler angle in radians
        //                      beta : Euler angle in radians
        //                      gamma: Euler angle in radians
        //                      rad  : Flag for angle in degrees or radians
        // Return               mx   : Euler rotation matrix on axes by angles
        //                             {alpha, beta, gamma} in successive
        //                             manner about z-y-z axes respectrively
        //                             Rotations are counter-clockwise about
        //                             the axes.
 


MSVCDLL friend matrix Rmx1(double alpha);

	// Input		pt   : Coordinate (this)
	// 			alpha: Euler angle in radians
	// Return		mx   : Rotation matrix
        ///F_list Rmx1		     - Rotation matrix for alpha rotation


MSVCDLL friend matrix Rmx2(double beta);

	// Input		pt   : Coordinate (this)
	// 			beta : Euler angle in radians
	// Return		mx   : Rotation matrix
        ///F_list Rmx2		     - Rotation matrix for beta rotation


MSVCDLL friend matrix Rmx3(double gamma);

	// Input		pt   : Coordinate (this)
	// 			gamma: Euler angle in radians
	// Return		mx   : Rotation matrix
        ///F_list Rmx3		     - Rotation matrix for gamma rotation


MSVCDLL friend matrix Rmx(double alpha, double beta, double gamma);

	// Input		pt   : Coordinate (this)
	// 			alpha: Euler angle in radians
	// 			beta : Euler angle in radians
	// 			gamma: Euler angle in radians
	// Return		mx   : Rotation matrix
	// Note			     : Returned matrix should be
	//			       equivalent to the product
	//			       Rmx3(gamma)*Rmx2(beta)*Rmx1(alpha)
        ///F_list Rmx		     - Rotation matrix for general rotation

MSVCDLL friend matrix Rmx(coord EA);

	// Input		EA   : Euler angles (radians)
	// Return		mx   : Rotation matrix
	// Note			     : Returned matrix should be
	//			       equivalent to the product
	//			       Rmx3(gamma)*Rmx2(beta)*Rmx1(alpha)
        ///F_list Rmx		     - Rotation matrix for general rotation


// ********************** Rotation of a Coordinate ****************************


MSVCDLL coord rotate(double alpha, double beta, double gamma);

	// Input		pt   : Coordinate (this)
	// 			alpha: Euler angles (radians)
	// 			beta : Euler angles (radians)
	// 			gamma: Euler angles (radians)
	// Return		rotpt: pt rotated by input Euler angles
        ///F_list rotate	     - Rotate a coordinate point


MSVCDLL coord rotate(coord& EA);

	// Input		pt   : Coordinate (this)
	// 			EA   : Euler angles (radians)
	// Return		rotpt: pt rotated by input Euler angles


// ____________________________________________________________________________
// E               CARTESIAN COORDINATE TRANSLATION FUNCTIONS
// ____________________________________________________________________________

///Center Coordinate Translation Functions
///F_list trans_x	        - Translate along x axis
///F_list trans_x_ip	        - Translate in place along x axis
///F_list trans_y	        - Translate along y axis
///F_list trans_y_ip	        - Translate in place along y axis
///F_list trans_z	        - Translate along z axis
///F_list trans_z_ip	        - Translate in place along z axis
///F_list translate	        - Translate coordinate point
///F_list translate_ip	        - Translate coordinate point in place

	// Input		pt      : Coordinate (this)
	// 			delx    : Change in x coordinate
	// 			dely    : Change in y coordinate
	// 			delz    : Change in z coordinate
	// Return (non-ip)	transpt : A coordinate point which is pt
	//				  translated by delx, dely, delz
	//                                or
	// Return (ip)          void	: Point (this) is translated by
	//				  delx, dely, delz


MSVCDLL coord trans_x(double    delx);
MSVCDLL void  trans_x_ip(double delx);
MSVCDLL coord trans_y(double    dely);
MSVCDLL void  trans_y_ip(double dely);
MSVCDLL coord trans_z(double    delz);
MSVCDLL void  trans_z_ip(double delz);

MSVCDLL coord translate(double    delx, double dely=0, double delz=0);
MSVCDLL void  translate_ip(double delx, double dely=0, double delz=0);

MSVCDLL coord translate(const    coord& del) const;
MSVCDLL void  translate_ip(const coord& del);

// ----------------------------------------------------------------------------
//             These Member Functions Also Perform A Translation
// ----------------------------------------------------------------------------

MSVCDLL coord operator +  (const coord& del) const;
MSVCDLL coord operator -  (const coord& del) const;
MSVCDLL coord&  operator += (const coord& del);
MSVCDLL coord&  operator -= (const coord& del);

// ____________________________________________________________________________
// F                  COORDINATE WITH COORDINATE
// ____________________________________________________________________________

///Center Coordinate with Coordinate Functions
 
// ----------------------------------------------------------------------------
//                      Distance Between Two Coordinates
// ----------------------------------------------------------------------------

MSVCDLL friend double Rad(const coord& pt1, const coord& pt2);

	// Input		pt1     : Coordinate point
	// 			pt2     : Coordinate point
	// Return		dist    : Distance between pt1 & pt2
        ///F_list Rad	        	- Distance between two points


MSVCDLL friend double theta(const coord& pt1, const coord& pt2);

	// Input		pt1     : Coordinate point
	// 			pt2     : Coordinate point
	// Return		theta   : Spherical angle, down from the
	//				  z axis for the vector connecting
	// 				  pt1 & pt2
        ///F_list theta	        	- Spherical angle theta between two points
	// Note			        : Two cases exist -

	//				  1.) R=0  --> theta=0
	//				  2.) R!=0 --> 0<=theta<=180

	//				  The latter case range coincides with
	//			  	  the principal values returned by the
	//				  arcosine function  


MSVCDLL friend double phi(const coord& pt1, const coord& pt2);

	// Input		pt1     : Coordinate point
	// 			pt2     : Coordinate point
	// Return		phi     : Spherical angle, over from the
	//			          x axis for the vector connecting
	// 			          pt1 & pt2
        ///F_list phi			- Spherical angle phi between two points
	// Note			        : Several cases exist -
	//					   x-axis
	//				1.) y=0, x>=0  --> phi=0
	//				2.) y=0, x<0   --> phi=180
	//					left hemisphere
	//				3.) y<0, x=0   --> phi=270
	//				4.) y<0, x<0   --> 180<=phi<=270
	//				5.) y<0, x>0   --> 270<=phi<=360
	//					right hemisphere
	//				6.) y>0, x=0   --> phi=90
	//				7.) y>0, x<0   --> 90<=phi<=180
	//				8.) y>0, x>0   --> 0<=phi<=90

	//			       Cases 4,5, & 7 lie outside the range
	//			       covered by the arctangent function.
	//			       The principal value range of atan is
	//			       [-90, 90]


MSVCDLL friend coord cdvect(const coord& pt1, const coord& pt2);

	// Input		pt1     : Coordinate point (Cartesian)
	// 			pt2     : Coordinate point (Cartesian)
	// Return		pt      : Coordinate point which corresponds
	//				  to the vector from pt1 to pt2
	//				  (Spherical)
        ///F_list vector		- Coordinate for vector between points

// ____________________________________________________________________________
// G                         COORDINATE WITH SCALAR
// ____________________________________________________________________________

///Center Coordinate with Scalar Functions

	// Input		pt1 : A coordinate point
	//			r   : A real number
	// Return		pt  : A new coordinate point which has
	//			      all components of pt1 scaled by scalar r
	// Return		pt  : A new coordinate point which has
	//			      all components of pt1 scaled by scalar r
	// Return		pt  : Coordinate point has had all
	//			      components of scaled by scalar r
	// Return		pt  : Coordinate point has had all
	//			      components of scaled by scalar 1/r
        ///F_list *		    - Coordinate multiplication by scalar
        ///F_list *=		    - Coordinate unary multiplication by scalar
        ///F_list /=		    - Coordinate division by scalar
        ///F_list /=		    - Coordinate unary division by scalar
	
MSVCDLL friend coord  operator *  (double r, const coord& pt1);
MSVCDLL coord  operator *  (double r) const;
MSVCDLL coord&   operator *= (double r);
MSVCDLL coord  operator /  (double r) const;
MSVCDLL coord&   operator /= (double r);

// ____________________________________________________________________________
// H                       COORDINATE WITH MATRIX
// ____________________________________________________________________________

///Center Coordinate with Matrix Functions
 
MSVCDLL friend coord operator * (const matrix& mx, const coord& pt);

	// Input		pt   : A coordinate point (this)
	// 			mx   : A matrix (must be 3x3)
	// Return		pt1  : A coordinate point which is
	//			         mx *  pt  =  pt1
	//			       (3x3)*(3x1) = (3x1)
        // Note                      : Does not check for proper matrix
        //                             dimensioning (speed) nor if the
        //                             transformation makes any sense   
        ///F_list *		     - Multiplication by (rotation) matrix


// ____________________________________________________________________________
// I                 CLASS COORDINATE WITH PARAMETERS
// ____________________________________________________________________________

///Center Coordinate with GAMMA Parameters & Parameter Sets

//-----------------------------------------------------------------------------
//                       Single Parameter Functions
//-----------------------------------------------------------------------------

/* These functions allow for generation of a GAMMA parameter equivvalent to 
   the coordinate. That is useful in writing coordinates to an output ASCII
   file that is readable by GAMMA. Thus, coordinates may be written and read
   to/from an ASCII file.

        // Input               pt    : A coordinate point (this)
        //                     pname : A parameter name
        // Return              par   : A GAMMA parameter of type coordinate
        //                             with the name pname and statment ps
	//			       if one is supplied.                   */
 
MSVCDLL SinglePar param(const std::string& pname) const;
MSVCDLL SinglePar param(const std::string& pname, const std::string& ps) const;

//-----------------------------------------------------------------------------
//                          Parameter Set Functions
//-----------------------------------------------------------------------------
 
        // Input                pt      : A coordinate point (this)
        //                      filein  : An input (aSCII) file
        //                      indx    : Point index
	//			warn	: Warning level
        //                      warn    : Warning level
        //                                  0 = no warnings
        //                                  1 = non-fatal warnings
        //                                  2 = fatal warnings
        // Output               none    : Coordinate point filled with
        //                                coordinate specified in filein


MSVCDLL int read(const std::string& filein, int indx, int warn=1);
MSVCDLL int read(const ParameterSet&  pset, int indx, int warn=1);

        // Input                pt      : A coordinate point (this)
        //                      pset    : A parameter set
        //                      indx    : Point index
        //                      warn    : Warning level
        //                                  0 = no warnings
        //                                  1 = non-fatal warnings
        //                                  2 = fatal warnings
        // Output               none    : Coordinate point filled with
        //                                coordinate in pset
        // Note                         : This function is used by the read
        //                                functions as well!

// ____________________________________________________________________________
// J                      CLASS COORDINATE I/O FUNCTIONS
// ____________________________________________________________________________

///Center Coordinate I/O Functions

//-----------------------------------------------------------------------------
//           Functions To Get & Set The Coordinate Output Format
//-----------------------------------------------------------------------------

MSVCDLL static int length();

MSVCDLL friend void coord_setf (int otype, int science, int digits, int precise);

        // Input                otype  : 0    == output as x,y,z
	//			         1    == output as r,theta,phi
	//				 2    == output as r,theta,z
        //                      science:TRUE  == output as 5.0e5
        //                              FALSE == output as 500000
        //                      digits :         number of digits total
        //                      precise:         digits after decimal point 
        //                                      (negative = any number possible)
        ///F_list coord_setf	       - Set output format


MSVCDLL friend void coord_getf (int& otype, int& science, int& digits, int& precise);

        // Output               otype  : 0    == output as x,y,z
	//			         1    == output as r,theta,phi
	//				 2    == output as r,theta,z
        //                      science:TRUE  == output as 5.0e5
        //                              FALSE == output as 500000
        //                      digits :         number of digits total
        //                      precise:         digits after decimal point 
        //                                      (negative = any number possible)
        ///F_list coord_getf	       - Get output format

//-----------------------------------------------------------------------------
//                   Functions To Output The Coordinate
//-----------------------------------------------------------------------------

        // Input                pt      : A coordinate point (this)
        //                      ostr	: Output stream
	// Output 		ostr	: Returns the modified output stream
        ///F_list <<			- Standard output

MSVCDLL        std::ostream& print(std::ostream& ostr) const;
MSVCDLL friend std::ostream& operator << (std::ostream& ostr, const coord& pt);

MSVCDLL friend std::istream& operator >> (std::istream& istr, coord& pt);

	// Input 		istr : Input stream
	//			pt   : Coordinate
	// Output 		istr : Coordinate set from input
	// Bug			     : No error handling
        ///F_list >>		     - Standard input
 
// ____________________________________________________________________________
// K                  CLASS COORDINATE CONVERSION FUNCTIONS
// ____________________________________________________________________________
 
        // Input                pt      : A coordinate point (this)
        //                      rad     : Radians vs. Deg angle output
        // Output (*2Sph)	pt1     : Equivalent coordinate in Spherical
        //                                space { R, theta, phi }
        // Output (*2Cart)      pt1     : Equivalent coordinate in Cartesian
        //                                space { x, y, z }
        // Output (*2Cyl)	pt1     : Equivalent coordinate in Cylindrical
        //                                space { R, theta, z }
	// Note (Cart2*)		: Input coordinate pt assumed
        //                                to be Cartesian {x, y, z }
        // Note (Sph2*)			: Input coordinate pt assumed
        //                                to be spherical { R, theta,phi }
        // Note (Cyl2*)			: Input coordinate pt assumed
        //                                to be cylindrical { R, theta, z }

MSVCDLL coord Cart2Sph(int rad=1) const;
MSVCDLL coord Sph2Cart(int rad=1) const;
MSVCDLL coord Cart2Cyl(int rad=1) const;
MSVCDLL coord Cyl2Cart(int rad=1) const;
MSVCDLL coord Sph2Cyl(int  rad=1) const;
MSVCDLL coord Cyl2Sph(int  rad=1) const;

// ____________________________________________________________________________
// L                  CLASS COORDINATE AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 
        // Input                pt      : A coordinate point (this)
	//			dpt	: An additional coordinate
        // Output               pt1     : The default coordinate
        //                 OR   void    : Default coordinate set to dpt
	// Note				: The default coordinate is a static
	//				  coordiante shared by all coordinates.
	//				  It is part of the class structure and
	//				  initialized in the implementation
 
MSVCDLL static coord getDefCoord();
MSVCDLL static void  setDefCoord(const coord& dpt);
 
// ____________________________________________________________________________
// M                  CLASS COORDINATE COMPARISON FUNCTIONS
// ____________________________________________________________________________
 
        // Input                pt      : A coordinate point (this)
        //                      pt2     : Another coordinate
        // Output ==            T/F     : TRUE if pt2 and this are equal
        // Output !=            T/F     : FALSE if pt2 and this are equal 
        ///F_list ==                    - Equality 
        ///F_list !=                    - Inequality 
 
MSVCDLL static void SetCutoff(double co=-1);
MSVCDLL        bool operator==(const coord& pt) const;
MSVCDLL        bool operator!=(const coord& pt) const;
MSVCDLL        bool operator>(const  coord& pt) const;
MSVCDLL        bool operator<(const  coord& pt) const;

};


/*****************************************************************************/
/*****************************************************************************/
/*                           CLASS COORD CONSTANTS                           */
/*****************************************************************************/
/*****************************************************************************/

extern const MSVCDLL coord UnitX;		// Unitx  = 1i + 0j + 0k
extern const MSVCDLL coord UnitY;		// Unity  = 0i + 1j + 0k
extern const MSVCDLL coord UnitZ;		// Unitz  = 0i + 0j + 1k
extern const MSVCDLL coord coord0;		// coord0 = 0i + 0j + 0k

#endif 								// coord.h

