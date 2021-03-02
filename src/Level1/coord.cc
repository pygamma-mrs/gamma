/* coord.cc *****************************************************-*-c++-*-
**									**
**                                G A M M A 				**
**									**
**      Coordinate                                Implementation	**
**									**
**	Copyright (c) 1990, 1991, 1992					**
**	Scott Smith						 	**
**	Eidgenoessische Technische Hochschule			 	**
**	Labor fuer physikalische Chemie					**
**	8092 Zurich / Switzerland				 	**
**								 	**
**      $Header: $
**								 	**
** Thanks be to Jacco van Beek (ETHZ) for finding the rotatons that     **
** make the Euler related functions consistent with the rest of GAMMA.  **
** These modifications are included as of version 4.0.4 in functions    **
** Rmx2, Rmx, and coord::rotate. They still need test programs in the   **
** GAMMA test suite.                                                    **
**								 	**
*************************************************************************/

/*************************************************************************
**									**
**  Description								**
**									**
** Class coordinate provides the  means with which a coordinate in 	**
** 3-dimensional space can be easily manipulated.  Each coordinate 	**
** consists of three real numbers which may be translated, rotated,	**
** etc. Functions are provided for component access, distance and 	**
** angle computations, and various other entities. 			**
** 									**
*************************************************************************/

#ifndef Coord_cc_			// Is file already included?
#define Coord_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Level1/coord.h>		// Include the header file
#include <Basics/Gconstants.h>		// Include constant PI
#include <Basics/StringCut.h>		// Additional string manipulations
#include <Basics/SinglePar.h>		// Include GAMMA parameters
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <Matrix/complex.h>		// Include GAMMA complex numbers
#include <Matrix/matrix.h>		// Include GAMMA matrices
#include <stdlib.h>
#include <math.h>
#include <list>				// Include libstdc++ STL lists
#include <string>			// Include libstdc++ strings
#include <iostream>			// Include form output function

using std::list;			// Using libstdc++ lists
using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams
using std::istream;			// Using libstdc++ input streams

// ----------------------------------------------------------------------------
//                           CLASS COORD CONSTANTS
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// a                   Intenal Constants (Not Exported)
// ____________________________________________________________________________

coord  coord::DefCoord = coord(0,0,0);	// Default coordinate
double coord::OrdCutoff = 1.e-13;	// Zero cutoff for ordinate

/*                      Individual Coordinate Output Format
                        -----------------------------------
 
        PTotype  : 0    == output as x,y,z       Cartesian (default)
                  1    == output as r,theta,phi Spherical
                  2    == output as r,theta,z   Cylindrical
        science: !0    == output as 5.0e5
                  0    == output as 500000.00
        digits :          number of total digits output
        precise:          digits after decimal point                         */

int    coord::PTotype=0;		// Static variables utilized in
int    coord::PTscience=0;		// specifying the output format
int    coord::PTdigits=6;
int    coord::PTprecise=2;
string coord::PTform = string("%6.2f");	

// ____________________________________________________________________________
// b                  External Constants (Not Exported)
// ____________________________________________________________________________

/* These constants were declared as extern in the class header. The must be
   initialized herein.                                                       */

const coord UnitX(1,0,0);		// unitx  = 1i + 0j + 0k
const coord UnitY(0,1,0);		// unity  = 0i + 1j + 0k
const coord UnitZ(0,0,1);		// unitz  = 0i + 0j + 1k
const coord coord0(0,0,0);		// coord0 = 0i + 0j + 0k

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS COORDINATE ERROR HANDLING
// ____________________________________________________________________________

/* These functions handle all error messages for class coord.  They tie into
   the generic GAMMA error messaging format.

		Input		pt      : A coordinate
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
		Output		none    : Error message output
                                          Execution stopped (if fatal)       */  
 
void coord::PTerror(int eidx, int noret) const
  {
  string hdr("Coordinate");
  string msg;
  switch(eidx)
    {
    case 3:  msg = string("Construction From Bad Parameter");      // (3)
             GAMMAerror(hdr,msg,noret); break;
    case 4:  msg = string("Bad Ordinate Request! Range [0,2]");    // (4)
             GAMMAerror(hdr,msg,noret); break;
    case 5:  msg = string("Can't Construct From Parameter Set");   // (5)
             GAMMAerror(hdr,msg,noret); break;
    default: GAMMAerror(hdr, eidx, noret); break;	// Unknown Error 
    }
  }

void coord::PTerror(int eidx, const string& pname, int noret) const
  {
  string hdr("Coordinate");
  string msg;
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr, 1, pname, noret);  break; // File Problems  (1)
    case 2:  GAMMAerror(hdr, 2, pname, noret);  break; // !Read Par      (2)


    case 17:							// (17)
      msg = string("Cannot Find Any Coordinate ") + pname
          + string(" Parameters!?");	
      GAMMAerror(hdr, msg, noret); break;
    case 18:							// (18)
      msg = string("Setting Coordinate ") + pname
          + string(" To Default Coordinate?");
      GAMMAerror(hdr, msg, noret); break;
    case 101:							// (101)
      msg = string("Can't Use Parameter File ") + pname;
      GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, -1, pname, noret); break; // Unknown Error  (-1)
    }
  }  

volatile void coord::PTfatal(int eidx) const
  {
  PTerror(eidx, eidx);                          // Output error message
  if(eidx) PTerror(0);                          // State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

volatile void coord::PTfatal(int eidx, const string& pname) const
  {
  PTerror(eidx, pname, eidx);			// Output error message
  if(eidx) PTerror(0);                          // State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                    CLASS COORDINATE SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions can translate parameters in a given GAMMA parameter set into
   a coordinate. These are priviate functions because often misuse of these
   types of functions can mess up the class structure. I do not see how that
   is possible for this class, but it is better to be safe. They are not
   needed by the user anyway.                                                */
 
        // Input		 pt	: A coordinate point (this) 
        //                       pset	: A parameter set
	//			 indx	: Point index
	//			 warn	: Warning output level
        // Output (*Cartesian)   TF  	: Coordinate point filled with
        //                                Cartesian coordinate in pset
        // Output (*Spherical)	 TF  	: Coordinate point filled with
        //                                Spherical coordinate in pset
        // Output (*Cylindrical) TF  	: Coordinate point filled with
        //                                Cylindrical coordinate in pset
	// Note				: Regardless of the input type, the
	//				  points will always be stored in
	//				  Cartesian space. Conversion is done
	//				  as needed.

/*
              Parameter Name      Space              Units
              --------------   -----------   ----------------------- 
                 Coord          Cartesian    {x,y,z} in Angstroms
                 Coordn         Cartesian    {x,y,z} in Nanometers
                 Coordm         Cartesian    {x,y,z} in Meters
                 CoordSph       Spherical    R in A,  {q,f} in deg.
                 CoordSphn      Spherical    R in nm, {q,f} in deg.
                 CoordSphm      Spherical    R in m,  {q,f} in deg.
                 CoordSphR      Spherical    R in A,  {q,f} in rad.
                 CoordSphnR     Spherical    R in nm, {q,f} in rad.
                 CoordSphmR     Spherical    R in m,  {q,f} in rad.
                 CoordCyl      Cylindrical   {R,z} in A,  q in deg.
                 CoordCyln     Cylindrical   {R,z} in nm, q in deg.
                 CoordCylm     Cylindrical   {R,z} in m,  q in deg.
                 CoordCylR     Cylindrical   {R,z} in A,  q in rad.
                 CoordCylnR    Cylindrical   {R,z} in nm, q in rad.
                 CoordCylmR    Cylindrical   {R,z} in m,  q in rad.
*/

bool coord::SetPtCartesian(const ParameterSet& pset, int indx, int warn)
  {
  int nn=3;					// # allowed parameter names
  string pname, pnames[3] =			// Possible parameter names
               { "Coord", "Coordn", "Coordm" };
  string app = string("(") + Gdec(indx) + ")";	// Parameter index append
  double sf[3] = { 1.e-10, 1.e-9, 1.0 };	// Scaling factors
  ParameterSet::const_iterator item;		// Pix in parameter list
  for(int i=0; i<nn; i++)			// Loop possible parameters
    {
    pname = pnames[i]; 				// Consruct parameter name
    if(indx != -1) pname += app;		// Use suffix if not -1 index
    item  = pset.seek(pname);			// Pix in parameter list
    if(item != pset.end())			// If the parameter is found
      {						// then set the coordinate
      *this = coord(*item)*sf[i];		// in meters and exit
      return true;
      }
    }
  return false;
  }

bool coord::SetPtSpherical(const ParameterSet& pset, int indx, int warn)
  {
  int nn=6;					// # allowed parameter names
  string pname, pnames[6] =			// Possible parameter names
   { "CoordSph",  "CoordSphn",  "CoordSphm",
     "CoordSphR", "CoordSphRn", "CoordSphRm" };
  string app = string("(") + Gdec(indx) + ")";	// Parameter index append
  double sf[3] = { 1.e-10, 1.e-9, 1.0 };	// Scaling factors
  double R, theta, phi;				// Spherical values
  ParameterSet::const_iterator item;		// Pix in parameter list
  coord pt;					// Temporary point
  for(int i=0; i<nn; i++)			// Loop possible parameters
    {
    pname = pnames[i] + app;			// Consruct parameter name
    item  = pset.seek(pname);			// Pix in parameter list
    if(item != pset.end())			// If the parameter is found
      {						// then store the coordinate
      pt = coord(*item);			// for the moment.
      R     = pt.cx*sf[i%3];			// Radius in meters
      theta = pt.cy;				// Azimuthal (rad or deg)
      phi   = pt.cz;				// From +x (rad or deg)
      switch(i)					// Insure angles in radians
        {					// before conversion into
        case 0:					// Cartesian units
        case 1:
        case 2:
          theta *= DEG2RAD;
          phi   *= DEG2RAD;
          break;
        default:
          break;
        }
      pt.cx = R*sin(theta)*cos(phi);		// Convert to Cartesian
      pt.cy = R*sin(theta)*sin(phi);		// from spherical
      pt.cz = R*cos(theta);
      *this = pt;
      return true;
      }
    }
  return false;
  }

bool coord::SetPtCylindrical(const ParameterSet& pset, int indx, int warn)
  {
  int nn=6;					// # allowed parameter names
  string pname, pnames[6] =			// Possible parameter names
   { "CoordCyl",  "CoordCyln",  "CoordCylm",
     "CoordCylR", "CoordCylRn", "CoordCylRm" };
  string app = string("(") + Gdec(indx) + ")";	// Parameter index append
  double sf[3] = { 1.e-10, 1.e-9, 1.0 };	// Scaling factors
  double R, theta, Z;				// Cylindrical values
  SinglePar par;				// Single parameter
  ParameterSet::const_iterator item;		// Pix in parameter list
  coord pt;					// Temporary point
  for(int i=0; i<nn; i++)			// Loop possible parameters
    {
    pname = pnames[i] + app;			// Consruct parameter name
    par   = SinglePar(pname);			// Construct parameter
    item  = pset.seek(pname);			// Pix in parameter list
    if(item != pset.end())			// If the parameter is found
      {						// then store the coordinate
      pt = coord(*item);			// for the moment.
      R     = pt.cx*sf[i%3];			// Radius (in plane) in meters
      theta = pt.cy;				// Azimuthal (rad or deg)
      Z     = pt.cz*sf[i%3];			// Height in meters
      switch(i)					// Insure angles in radians
        {					// before conversion into
        case 0:					// Cartesian units
        case 1:
        case 2:
          theta *= DEG2RAD;
          break;
        default:
          break;
        }
      pt.cx = R*cos(theta);			// Convert to Cartesian
      pt.cy = R*sin(theta);			// from cylindrical
      pt.cz = Z;
      *this = pt;
      return true;
      }
    }
  return false;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A          CLASS COORDINATE CONSTRUCTORS/DETRUCTOR AND ASSIGNMENT
// ____________________________________________________________________________

/* These constructors set up a new coordinate.  There are several ways to make
   a coordinate as listed below:
 
          Input Arguments                       Resulting Coordinate
   --------------------------        -----------------------------------------
                -                    Zero coordinate, {0,0,0}
   double x,double y,double z        Coordinate {x,y,z}, y & z 0 by default
             coord pt                Coordinate { pt.x, pt.y, pt.z }
     ParameterSet, idx, warn	     Coordinate with index i read from pset 
             SinglePar               Coordinate set from parameter

   Setting a coordinate from a parameter set may fail if the parameters that
   define the coordinate don't exist. The flag warn dictates what happens
   upon a failure: 0=no warnings, 1=non-fatal warnings, 2=fatal warnings
   Setting a coordinate from a single parameter assumes a parameter of type
   3 with value stored as a string "( #, #, # )". All I/O between GAMMA
   coordinates and single parameters must follow this format!	           */
	
coord::coord( )                               { cx = 0;  cy = 0;  cz = 0;  }
coord::coord(double xx, double yy, double zz) { cx = xx; cy = yy; cz = zz; }
coord::coord(const coord& P)                  { cx=P.cx; cy=P.cy; cz=P.cz; }

coord::coord(const ParameterSet& pset, int idx, int warn)
  {
  if(SetPtCartesian(pset,   idx, 0)) return;	// First try Cartesian point
  if(SetPtSpherical(pset,   idx, 0)) return;	// Next try Spherical point
  if(SetPtCylindrical(pset, idx, 0)) return;	// Next try Cylindrical point
  if(warn)					// Looks like we can't read
    {						// the coordinate point in
    string sl("");
    if(idx != -1) sl = string(Gdec(idx));
    PTerror(17, sl, 1);				// Can't find coord info
    if(warn > 1)  PTfatal(5);			// Can't construct from pset
    else          PTerror(18, sl, 1);		// Setting default coord
    }
  *this = DefCoord;				// Use default coord
  }

coord::coord(const SinglePar& par)
  {
  if(par.type() != 3) PTfatal(3);	// Exit if wrong parameter type
  string val1 = par.data();		// Get the parameter data
  cutBlksXBlks(val1,"(");		// Remove blanks-(-blanks at start
  string val2 = cutDouble(val1,0);	// Remove the first double
  cx = atof(val2.c_str());		// Set the first ordinate
  cutBlksXBlks(val1,",");		// Remove blanks-,-blanks at start
  val2 = cutDouble(val1,0);		// Remove the second double
  cy = atof(val2.c_str());		// Set the second ordinate
  cutBlksXBlks(val1,",");		// Remove blanks-,-blanks at start
  val2 = cutDouble(val1,0);		// Remove the third double
  cz = atof(val2.c_str());		// Set the third ordinate
  }

coord::~coord () { }			// No need to destruct anything

coord& coord::operator= (const coord& pt) 
  { 
  if(this == &pt) return *this;		// Do not bother to self assign
  cx=pt.cx;				// Copy the x (R) ordinate
  cy=pt.cy;				// Copy the y (theta) ordinate
  cz=pt.cz;				// Copy the z (phi) ordinate
  return *this;
  }

// ____________________________________________________________________________
// B                  CARTESIAN COORDINATE ACCESS FUNCTIONS
// ____________________________________________________________________________
 
/*      Function                                     Result
    -----------------             ---------------------------------------------
       get(int)                   Return ordinate specified by int
       x,x(double)                Get/Set first ordinate
       y,y(double)                Get/Set second ordinate
       z,z(double)                Get/Set third ordinate
       xyz(...)                   Set all ordinates                          */

double coord::get(int i) const
  {
  if(i<0 || i>2) PTfatal(4);			// Bad ordinate request
  if(!i)        return cx;
  else if(i==1) return cy;
  return cz;
  }

double coord::x()          const { return cx; }
void   coord::x(double xx)       { cx = xx;   }
double coord::y()          const { return cy; }
void   coord::y(double yy)       { cy = yy;   }
double coord::z()          const { return cz; }
void   coord::z(double zz)       { cz = zz;   }
void   coord::xyz(double xx, double yy, double zz)     { cx=xx; cy=yy; cz=zz; }
void   coord::xyz(const coord& pt1)        { cx=pt1.cx; cy=pt1.cy; cz=pt1.cz; }
double coord::norm() const { return sqrt(cx*cx + cy*cy + cz*cz); }

// ____________________________________________________________________________
// C                SPHERICAL BASED COORDINATE ACCESS FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                    Rad is The Distance From The Origin
// ----------------------------------------------------------------------------

	// Input		pt	: Coordinate (this)
	//			pt2	: A second coordinate
	//			{x,y,z} : Cartesian ordinates
	// Return Rad()		Rad	: Radius or length from (0,0,0)
	// Return Rad(pt2)	Rad	: Radius or length vector from
	//				  between the two points
	// Return Rad(x,y,z)	Rad	: Radius or length from (0,0,0)

double coord::Rad() const { return sqrt(cx*cx + cy*cy + cz*cz); }
double coord::Rad(const coord& pt2) const
  { 
  double delx = pt2.cx-cx;
  double dely = pt2.cy-cy;
  double delz = pt2.cz-cz;
  return sqrt(delx*delx + dely*dely + delz*delz);
  }

double Rad(double x, double y, double z) { return sqrt(x*x + y*y + z*z); }

// ----------------------------------------------------------------------------
//	         Theta is The Angle Down From The Positive Z Axis
// ----------------------------------------------------------------------------


	// Input		pt    : Coordinate (this)
	// Return		theta : Spherical angle, down from +z axis
	// Note			      : Two cases exist -
	//				1.) R=0  --> theta=0
	//				2.) R!=0 --> 0<=theta<=180
	//				The latter case range coincides with
	//				the principal values returned by the
	//				arcosine function  

double coord::theta() const
  {
  double rad = Rad();
  if(rad) return acos(cz/rad);
  else    return 0.0;
  }



	// Input		pt    : Coordinate (this)
	//			pt2   : A second coordinate
	// Return		theta : Spherical angle, down from +z axis
	//				of the vector connecting pt to pt2
	// Note			      : Two cases exist -
	//				1.) R=0  --> theta=0
	//				2.) R!=0 --> 0<=theta<=180
	//				The latter case range coincides with
	//				the principal values returned by the
	//				arcosine function  

double coord::theta(const coord& pt2) const
  {
  double rad = Rad(pt2);
  if(rad) return acos((pt2.z()-cz)/rad);
  else    return 0.0;
  }



	// Input		pt    : Coordinate (this)
	// Return		theta : Spherical angle, down from the
	//				z axis
	// Note			      : Two cases exist -
	//				1.) R=0  --> theta=0
	//				2.) R!=0 --> 0<=theta<=180
	//				The latter case range coincides with
	//				the principal values returned by the
	//				arcosine function  

double theta(double x, double y, double z)
  {
  double rad = Rad(x,y,z);
  if(rad) return acos(z/rad);
  else    return 0.0;
  }

// ----------------------------------------------------------------------------
//	         Phi is The Angle Over From The Positive X Axis
// ----------------------------------------------------------------------------


double coord::phi() const

	// Input		pt    : Coordinate (this)
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

  {
  if(cy == 0)					// On x axis
    {
    if(cx >= 0) return 0.0; 			// 1.) On positive x-axis
    else return PI; 				// 2.) On negative x-axis
    }
  else if(cy < 0)				// In left hemisphere
    {
    if(cx == 0)     return 1.5*PI;		// 3.) On negative y-axis
    else if(cx < 0) return (atan(cy/cx)+PI);	// 4.) 180 < phi < 270
    else            return (atan(cy/cx)+2.*PI);	// 5.) 270 < phi < 360
    }
  else						// In right hemisphere
    {
    if(cx == 0) return 0.5*PI;			// 6.) On positive y-axis
    else if(cx < 0) return (atan(cy/cx)+PI);	// 7.) 90 < phi < 180
    else return atan(cy/cx);			// 8.) 0 < phi < 90
    }
  }



	// Input		pt	: Coordinate (this)
	//			pt2	: A second coordinate
	// Return		phi	: Spherical angle, over from the
	//			          x axis for the vector connecting
	// 			          pt & pt2

double coord::phi(const coord& pt2) const
  {
  double delx = pt2.x()-cx;
  double dely = pt2.y()-cy;
  coord pt(delx, dely, 0.0);
  return pt.phi();
  }



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

double phi(double x, double y, double z)
  {
  if(y == 0)					// On x axis
    {
    if(x >= 0) return 0.0; 			// 1.) On positive x-axis
    else       return PI; 			// 2.) On negative x-axis
    }
  else if(y < 0)				// In left hemisphere
    {
    if(x == 0) return 1.5*PI;			// 3.) On negative y-axis
    else if(x < 0) return (atan(y/x)+PI);	// 4.) 180 < phi < 270
    else return (atan(y/x)+2.*PI);		// 5.) 270 < phi < 360
    }
  else						// In right hemisphere
    {
    if(x == 0) return 0.5*PI;			// 6.) On positive y-axis
    else if(x < 0) return (atan(y/x)+PI);	// 7.) 90 < phi < 180
    else return atan(y/x);			// 8.) 0 < phi < 90
    }
  z=0;						// Compiler likes z to be used
  }

void coord::invert() { cx = 1/cx; cy = 1/cy; cz = 1/cz; }

	// Input		pt    : Coordinate (this)
	// Return		none  : All three coordinate values are
	//				inverted

// ____________________________________________________________________________
// D                    COORDINATE ROTATION FUNCTIONS
// ____________________________________________________________________________

/* Rotations can be a wee bit troublesome due to problems of defining what is
   being rotated and in which direction.  Rotations can either be viewed as
   directly rotating the point leaving the axes alone or as rotating the axis
   system then outputting the point relative to the new axes.  For a single
   rotation these are essentially identical except they use opposite angles.
   It is successive rotations about multiple axes where troubles can ensue.
   We'll try and be very careful in stating how each rotation function is
   applied.                                                                  */

// ----------------------------------------------------------------------------
//            Typical Rotations Active On Cartesian Coordinates
// ----------------------------------------------------------------------------

/* For these functions users can assume that the AXES REMAIN STATIC and that
   the COORDINATES ROTATE.  The rotation will occur in a COUNTER-CLOCKWISE
   fashion about the rotation axis. This is the opposite of an Euler rotation.

   If one insists on viewing the axes as rotating, the axes will rotate in a
   clockwise fashion, the opposite of Euler type of rotations.  Changing the
   sign on the rotation angle switches the rotation direction.  Thus negative
   angle input coincides with an Euler type of rotation (i.e. axes move in
   counter-clockwise fashion or points move in clockwise fashion.            */


	// Input		pt   : Coordinate (this)
	// 			phi  : Spherical angle (degrees over from +x)
	//			rad  : Flag for angle in degrees or radians
	// Return		mx   : Rotation matrix to rotate point by
	//			       angle phi c-clockwise about the x-axis

matrix coord::Rz(double phi, int rad)
  {
  matrix R(3,3,complex0);
  if(!rad) phi *= DEG2RAD;
  double cp = cos(phi);		//           [ cos(phi) -sin(phi)  0 ]
  double sp = sin(phi);		//  	     [                       ]
  R.put(cp, 0, 0);		// R (phi) = [ sin(phi)  cos(phi)  0 ]
  R.put(-sp, 0, 1);		//  z	     [                       ]
  R.put(sp, 1, 0);		// 	     [      0       0      1 ] 
  R.put(cp, 1, 1);
  R.put(1, 2, 2);
  return R;
  }

	// Input		pt   : Coordinate (this)
	// 			theta: Spherical angle (degrees down from +z)
	//			rad  : Flag for angle in degrees or radians
	// Return		mx   : Rotation matrix to rotate point by
	//			       angle theta c-clockwise about the x-axis

matrix coord::Rx(double theta, int rad)
  {
  matrix R(3,3,complex0);
  if(!rad) theta *= DEG2RAD;
  double ct = cos(theta);	//             [ 1      0          0       ]
  double st = sin(theta);	//             [                           ]
  R.put(1,  0,0);		// R (theta) = [ 0  cos(theta) -sin(theta) ]
  R.put(ct, 1,1);		//  x	       [                           ]
  R.put(-st, 1,2);		//	       [ 0  sin(theta)  cos(theta) ]
  R.put(st,2,1);
  R.put(ct, 2,2);
  return R;
  }



	// Input		pt   : Coordinate (this)
	// 			theta: Spherical angle (degrees down from +z)
	//			rad  : Flag for angle in degrees or radians
	// Return		mx   : Rotation matrix to rotate point by
	//			       angle theta c-clockwise about the y-axis

matrix coord::Ry(double theta, int rad)
  {
  matrix R(3,3,complex0);
  if(!rad) theta *= DEG2RAD;
  double ct = cos(theta);	//             [  cos(theta) 0  sin(theta) ]
  double st = sin(theta);	//             [                           ]
  R.put(ct, 0,0);		// R (theta) = [      0      1      0      ]
  R.put(st, 0,2);		//  y	       [                           ]
  R.put(1,  1,1);		//	       [ -sin(theta) 0  cos(theta) ]
  R.put(-st,2,0);
  R.put(ct, 2,2);
  return R;
  }


coord coord::xrotate(double theta, int rad) const
 
        // Input                pt   : Coordinate (this)
        //                      theta: Spherical angle (degrees down from +z)
        //                      rad  : Flag for angle in degrees or radians
        // Return               rotpt: The pt rotated about the x-axis by
	//			       the input angle
        // Note                      : The reference axes are static!  Thus
        //                             multiple rotations should be cumulative

  { return Rx(theta,rad)*(*this); }
 

coord coord::yrotate(double theta, int rad) const

        // Input                pt   : Coordinate (this)
        //                      theta: Spherical angle (degrees down from +z)
        //                      rad  : Flag for angle in degrees or radians
        // Return               rotpt: The pt rotated about the y-axis by
        //                             the input angle
        // Note                      : The reference axes are static!  Thus
        //                             multiple rotations should be cumulative

  { return Ry(theta,rad)*(*this); }
 
 
coord coord::zrotate(double phi, int rad) const { return Rz(phi,rad)*(*this); }
 
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
   & that the COORDINATES ARE STATIC.  This is normal for Euler rotations.

           Input                pt   : Coordinate (this)
                                alpha: Euler angle in radians
                                beta : Euler angle in radians
	   			gamma: Euler angle in radians
                                rad  : Flag for angle in degrees or radians
           Return               mx   : Euler rotation matrix on set axis by
                                       angle counter-clockwise about axis
                                       or about Euler angles about set axes

     Function                         Result
 -----------------   -------------------------------------------------------
      Ralpha         Rotation of angle alpha counter-clockwise about z axis
      Rbeta          Rotation of angle beta  counter-clockwise about y axis
      Rgamma         Rotation of angle gamma counter-clockwise about z axis 
      REuler         Successive rotations of angles alpha,beta,gamma above   */
                    
matrix coord::Ralpha(double alpha, int rad) const
                                              { return UnitZ.Rz(-alpha, rad); } 
matrix coord::Rbeta(double beta,   int rad) const
                                              { return UnitZ.Ry(-beta, rad);  }
matrix coord::Rgamma(double gamma, int rad) const
                                              { return UnitZ.Rz(-gamma, rad); } 

matrix coord::REuler(double alpha, double beta, double gamma, int rad) const
  { return Rgamma(gamma,rad)*Rbeta(beta,rad)*Ralpha(alpha,rad); }


matrix Rmx1(double alpha) { return UnitZ.Rz(-alpha, 1); } 

	// Input		pt   : Coordinate (this)
	// 			alpha: Euler angle in radians
	// Return		mx   : Rotation matrix
/*
  {
  matrix R(3,3,complex0);
  double ca = cos(alpha);	//            [ cos(alpha)  sin(alpha)  0 ]
  double sa = sin(alpha);	// 	      [                           ]
  R.put(ca, 0, 0);		// R(alpha) = [ -sin(alpha) cos(alpha)  0 ]
  R.put(sa, 0, 1);		// 	      [                           ]
  R.put(-sa, 1, 0);		// 	      [      0          0       1 ] 
  R.put(ca, 1, 1);
  R.put(1, 2, 2);
  return R;
  }
*/


matrix Rmx2(double beta) { return UnitX.Rx(-beta, 1); } 

	// Input		pt   : Coordinate (this)
	// 			beta : Euler angle in radians
	// Return		mx   : Rotation matrix

/*
  {
  matrix R(3,3,complex0);
  double cb = cos(beta);        //           [ cos(beta)  0  -sin(beta)]
  double sb = sin(beta);        //           [                         ]
  R.put( cb,0,0);               // R(beta) = [ 0          1      0     ]
  R.put(-sb,0,2);               //           [                         ]
  R.put(  1,1,1);               //           [ sin(beta)  0  cos(beta) ]
  R.put( sb,2,0);
  R.put( cb,2,2);
  return R;
  }
*/


matrix Rmx3(double gamma) { return UnitZ.Rz(-gamma, 1); } 

	// Input		pt   : Coordinate (this)
	// 			gamma: Euler angle in radians
	// Return		mx   : Rotation matrix

/*
  {
  matrix R(3,3,complex0);
  double cg = cos(gamma);	//            [ cos(gamma)  sin(gamma) 0 ]
  double sg = sin(gamma);	//            [                          ]
  R.put(cg,0,0);		// R(gamma) = [ -sin(gamma) cos(gamma) 0 ]
  R.put(sg,0,1);		//            [                          ]
  R.put(-sg,1,0);		//            [     0           0      1 ]                
  R.put(cg,1,1);
  R.put(1,2,2);
  return R;
  }
*/

	// Input		pt   : Coordinate (this)
	// 			alpha: Euler angle in radians
	// 			beta : Euler angle in radians
	// 			gamma: Euler angle in radians
	// Return		mx   : Rotation matrix
	// Note			     : Returned matrix should be
	//			       equivalent to the product
	//			       Rmx3(gamma)*Rmx2(beta)*Rmx1(alpha)

	// Input		pt   : Coordinate (this)
	// 			EA   : Coordinate containing Euler angle
	// Return		mx   : Rotation matrix
	// Note			     : Euler angles should be radians

	// Input		pt   : Coordinate (this)
	// 			alpha: Euler angle (radians)
	// 			beta : Euler angle (radians)
	// 			gamma: Euler angle (radians)
	// Return		rotpt: pt rotated by input Euler angles

matrix Rmx(double alpha, double beta, double gamma)
  {
  matrix mx(3,3);
  double ca = cos(alpha);       //    [                                    ]
  double sa = sin(alpha);       //    |  CaCbCg-SaSg   SaCbCg+CaSg   -SbCg |
  double cb = cos(beta);        //    |                                    |
  double sb = sin(beta);        //    |                                    |
  double cg = cos(gamma);       // R =| -SaCg-CaCbSg  -SaCbSg+CaCg    SbSg |
  double sg = sin(gamma);       //    |                                    |
  mx.put(ca*cb*cg - sa*sg,0,0); //    |                                    |
  mx.put(sa*cb*cg + ca*sg,0,1); //    |      CaSb         SaSb         Cb  |
  mx.put(-sb*cg,0,2);           //    [                                    ]
  mx.put(-(sa*cg+ca*sg*cb),1,0);
  mx.put(-sa*cb*sg+ca*cg,1,1);
  mx.put(sb*sg,1,2);
  mx.put(ca*sb,2,0);
  mx.put(sa*sb,2,1);
  mx.put(cb,2,2);
  return mx;
  }

matrix Rmx(coord EA) { return Rmx(EA.x(), EA.y(), EA.z()); }

coord coord::rotate(double alpha, double beta, double gamma)
  {
  coord rotpt;
  double ca = cos(alpha);
  double sa = sin(alpha);
  double cb = cos(beta);
  double sb = sin(beta);
  double cg = cos(gamma);
  double sg = sin(gamma);
  rotpt.cx  = (ca*cb*cg - sa*sg)*cx;		// <0|R|0>*Xo
  rotpt.cx += (sa*cb*cg + ca*sg)*cy;		// <0|R|1>*Yo
  rotpt.cx -= sb*cg*cz;				// <0|R|2>*Zo
  rotpt.cy  = -(sa*cg+ca*cb*sg)*cx;		// <1|R|0>*Xo
  rotpt.cy += (-sa*cb*sg+ca*cg)*cy;		// <1|R|1>*Yo
  rotpt.cy += sb*sg*cz;				// <1|R|2>*Zo
  rotpt.cz  = ca*sb*cx;				// <2|R|0>*Xo
  rotpt.cz += sa*sb*cy;				// <2|R|1>*Yo
  rotpt.cz += cb*cz;				// <2|R|2>*Zo
  return rotpt;
  }

	// Input		pt   : Coordinate (this)
	// 			EA   : Euler angles (radians)
	// Return		rotpt: pt rotated by input Euler angles

coord coord::rotate(coord& EA)
  { return (*this).rotate(EA.x(), EA.y(), EA.z()); }


// ____________________________________________________________________________
// E                     COORDINATE TRANSLATION FUNCTIONS
// ____________________________________________________________________________

/* These functions are used to translate a coordinate along a chosen axis.

           Input             pt      : Coordinate (this)
	   		     delu    : Change in u coordinate, u=x,y,z
           Return            transpt : Coordinate after translation
                        OR   void    : The "ip" funcitons act directly on pt

     Function                         Result
 -----------------   -------------------------------------------------------
     trans_u         Return new point translated along axis u, u={x,y,z}
    trans_u_ip       Translate point along axis u, u={x,y,z}
    translate        Return new point translated on all axes
TranslateSuccessive rotations of angles alpha,beta,gamma above   */
                    

coord coord::trans_x(double delx)
  { coord transpt(*this); transpt.cx += delx; return transpt; }
coord coord::trans_y(double dely)
  { coord transpt(*this); transpt.cy += dely; return transpt; }
coord coord::trans_z(double delz)
  { coord transpt(*this); transpt.cz += delz; return transpt; }

void coord::trans_x_ip(double delx) { cx += delx; }
void coord::trans_y_ip(double dely) { cy += dely; }
void coord::trans_z_ip(double delz) { cz += delz; }

coord coord::translate(double delx, double dely, double delz)
  { return coord(cx+delx, cy+dely, cz+delz); }

coord coord::translate(const coord& del) const
  { return coord(cx+del.x(), cy+del.y(), cz+del.z()); }

void coord::translate_ip(double delx, double dely, double delz)
  { cx += delx; cy += dely; cz += delz; }

void coord::translate_ip(const coord& del)
  { cx += del.x(); cy += del.y(); cz += del.z(); }

// ----------------------------------------------------------------------------
//             These Member Functions Also Perform A Translation
// ----------------------------------------------------------------------------

coord coord::operator + (const coord& del) const
  { return coord(cx+del.cx, cy+del.cy, cz+del.cz); }

coord coord::operator - (const coord& del) const
  { return coord(cx-del.cx, cy-del.cy, cz-del.cz); }

coord& coord::operator +=(const coord& del)
  { cx += del.x(); cy += del.y(); cz += del.z(); return (*this); }

coord& coord::operator -=(const coord& del)
  { cx -= del.x(); cy -= del.y(); cz -= del.z(); return (*this); }

// ____________________________________________________________________________
// F                        COORDINATE WITH COORDINATE
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//	                Distance Between Two Coordinates
// ----------------------------------------------------------------------------

double Rad(const coord& pt1, const coord& pt2)

	// Input		pt1     : Coordinate point
	// 			pt2     : Coordinate point
	// Return		dist    : Distance between pt1 & pt2

  {
  double delx = pt2.cx - pt1.cx;
  double dely = pt2.cy - pt1.cy;
  double delz = pt2.cz - pt1.cz;
  return Rad(delx, dely, delz);
  }

double theta(const coord& pt1, const coord& pt2)

	// Input		pt1     : Coordinate point
	// 			pt2     : Coordinate point
	// Return		theta   : Spherical angle, down from the
	//				  z axis for the vector connecting
	// 				  pt1 & pt2
	// Note			        : Two cases exist -

	//				  1.) R=0  --> theta=0
	//				  2.) R!=0 --> 0<=theta<=180

	//				  The latter case range coincides with
	//			  	  the principal values returned by the
	//				  arcosine function  

  {
  double rad = Rad(pt1, pt2);
  if(rad) return acos((pt2.z()-pt1.z())/rad);
  else    return 0.0;
  }


double phi(const coord& pt1, const coord& pt2)

	// Input		pt1  : Coordinate point
	// 			pt2  : Coordinate point
	// Return		phi  : Spherical angle, over from the
	//			       x axis for the vector connecting
	// 			       pt1 & pt2
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

  {
  double x = pt2.x()-pt1.x();			// Calculate delta x
  double y = pt2.y()-pt1.y();			// Calculate delta y
  return phi(x,y,0.0);
  }


coord cdvect(const coord& pt1, const coord& pt2)

	// Input		pt1     : Coordinate point (Cartesian)
	// 			pt2     : Coordinate point (Cartesian)
	// Return		pt      : Coordinate point which corresponds
	//				  to the vector from pt1 to pt2
	//				  (Spherical)

  {
  double rad = Rad(pt1, pt2);
  double the = theta(pt1, pt2);
  double ph  = phi(pt1, pt2);
  coord pt(rad,the,ph);
  return pt;
  }


// ____________________________________________________________________________
// G                        COORDINATE WITH SCALAR
// ____________________________________________________________________________

	// Input		pt1 : A coordinate point
	//			r   : A real number
	// Return		pt  : Coordinate point which has
	//			      components scaled by scalar r (or 1/r)

coord operator * (double r, const coord& pt1)
                              { return coord(pt1.cx*r, pt1.cy*r, pt1.cz*r); }
coord coord::operator *  (double r) const { return coord(cx*r, cy*r, cz*r); }
coord&  coord::operator *= (double r)       { cx *= r; cy *= r; cz *= r; return (*this); }
coord coord::operator /  (double r) const { return coord(cx/r, cy/r, cz/r); }
coord&  coord::operator /= (double r)       { cx /= r; cy /= r; cz /= r; return (*this); }

// ____________________________________________________________________________
// G                       COORDINATE WITH MATRIX
// ____________________________________________________________________________
 

	// Input		pt   : A coordinate point (this)
	// 			mx   : A matrix (must be 3x3)
	// Return		pt1  : A coordinate point which is
	//			         mx *  pt  =  pt1
	//			       (3x3)*(3x1) = (3x1)
	// Note			     : Does not check for proper matrix
	// 			       dimensioning (speed) nor if the
	//			       transformation makes any sense	

coord operator * (const matrix& mx, const coord& pt)
  {
  coord pt1;
  pt1.cx =  mx.getRe(0,0)*pt.cx; 
  pt1.cx += mx.getRe(0,1)*pt.cy;
  pt1.cx += mx.getRe(0,2)*pt.cz;
  pt1.cy =  mx.getRe(1,0)*pt.cx; 
  pt1.cy += mx.getRe(1,1)*pt.cy;
  pt1.cy += mx.getRe(1,2)*pt.cz;
  pt1.cz =  mx.getRe(2,0)*pt.cx; 
  pt1.cz += mx.getRe(2,1)*pt.cy;
  pt1.cz += mx.getRe(2,2)*pt.cz;
  return pt1;
  }


// ____________________________________________________________________________
// H        CLASS COORDINATE WITH PARAMETERS & PARAMETER SETS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//                       Single Parameter Functions
//-----------------------------------------------------------------------------
 
 
        // Input               pt    : A coordinate point (this)
        //                     pname : A parameter name
        // Return              par   : A GAMMA parameter of type coordinate
        //                             with the name pname

SinglePar coord::param(const string& pname) const
  {
  string pstate = "Coordinate Point";	// Default statement
  return param(pname, pstate);		// Use overload function
  }
 
 
 
        // Input               pt    : A coordinate point (this) 
        //                     pname : A parameter name
        //                     pstate: A parameter statement
        // Return              par   : A GAMMA parameter of type coordinate
        //                             with the name pname and statement pstate
 
SinglePar coord::param(const string& pname, const string& pstate) const
  {
  string pdata = "(";			// Forming coordinate data
  pdata += Gform("%g", cx);			// string for parameter output
  pdata += string(", ");			// in ASCII.  This will look
  pdata += Gform("%g", cy);			// like
  pdata += string(", ");
  pdata += Gform("%g", cz);			// ( x.xxx, y.yyy, z.zzz)
  pdata += string(") ");
  SinglePar par(pname,3,pdata,pstate);		// Now construct parameter
  return par;					// Return the parameter
  }

//-----------------------------------------------------------------------------
//                          Parameter Set Functions
//-----------------------------------------------------------------------------
 
int coord::read(const string& filein, int indx, int warn)
 
        // Input		pt	: A coordinate point (this) 
        //                      filein	: An input (aSCII) file
	//			indx    : Point index
	//			warn    : Warning output level
        // Output               none	: Coordinate point filled with
        //                                coordinate specified in filein

  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filein, warn?1:0))	// Read in pset from file
    {   
    if(warn)
      {
      PTerror(1, filein, 1);		// Problems with file filein
      if(warn>1) PTfatal(101,filein);	// Can't use parameter file
      else       PTerror(101,1);	// Cant use parameter file
      }
    return 0;  
    }   
  return read(pset, indx, warn);
  }

        // Input		pt	: A coordinate point (this) 
        //                      pset	: A parameter set
	//			indx    : Point index
        //                      warn    : Warning level
        //                                  0 = no warnings
        //                                  1 = non-fatal warnings
        //                                  2 = fatal warnings
        // Output               none	: Coordinate point filled with
        //                                coordinate in pset
        // Note				: This function is used by the read
        //				  functions as well!

int coord::read(const ParameterSet& pset, int indx, int warn)
  {
  if(SetPtCartesian(pset, indx, 0)) return 1;	// First try Cartesian point
  if(SetPtSpherical(pset, indx, 0)) return 1;	// Next try Spherical point
  if(SetPtCylindrical(pset,indx,0)) return 1; 	// Next try Cylindrical point
  if(warn)					// Looks like we can't read
    {						// the coordinate point in
    string sl("");
    if(indx != -1) sl = string(Gdec(indx));
    PTerror(17, sl, 1);				// Can't find coord info
    if(warn > 1)  PTfatal(5);			// Can't construct from pset
    else          PTerror(18, sl);		// Setting default coordinate
    }
  *this = DefCoord;				// Use default coord
  return 0;
  }

// ____________________________________________________________________________
// I                     CLASS COORDINATE I/O FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//           Functions To Get & Set The Coordinate Output Format 
//-----------------------------------------------------------------------------

int coord::length() { return 6 + 3*PTdigits; }

        // Input                otype  : 0    == output as x,y,z
        //                               1    == output as r,theta,phi
        //                               2    == output as r,theta,z
        //                      science:TRUE  == output as 5.0e5
        //                              FALSE == output as 500000
        //                      digits :         number of digits total
        //                      precise:         digits after decimal point
        //                                      (negative = any number possible)

void coord_setf(int otype, int science, int digits, int precise)
  {
  coord::PTotype=otype;
  coord::PTscience=science;
  coord::PTdigits=digits;
  coord::PTprecise=precise;
  char ef = (science)?'e':'f';
  coord::PTform = string("%");
  coord::PTform += Gdec(digits);
  if(precise>=0) 
    coord::PTform += string(".") + Gdec(precise);
  coord::PTform += ef;
  }

        // Output               otype  : 0    == output as x,y,z
        //                               1    == output as r,theta,phi
        //                               2    == output as r,theta,z
        //                      science:TRUE  == output as 5.0e5
        //                              FALSE == output as 500000
        //                      digits :         number of digits total
        //                      precise:         digits after decimal point
        //                                      (negative = any number possible)

void coord_getf(int& otype, int& science, int& digits, int& precise)
  {
  otype=coord::PTotype;
  science=coord::PTscience;
  digits=coord::PTdigits;
  precise=coord::PTprecise;
  }

//-----------------------------------------------------------------------------
//                   Functions To Output The Coordinate
//-----------------------------------------------------------------------------

        // Input		pt	: A coordinate point (this) 
	//				  assumed to be Cartesian {x, y, z }
        //                      ostr    : Output stream
        // Output               none    : Coordinate info placed into the
	//				  output stream ostr

ostream& coord::print(ostream& ostr) const
  {
  double x=0,y=0,z=0;
  switch(coord::PTotype)
    {
    case 0:				// Output Cartesian coordinates
    default:
      x = cx;
      y = cy;
      z = cz;
      break;
    case 1:				// Output Spherical coordinates
      x = Rad();
      y = theta();
      z = phi();
      break;
    case 2:				// Output Cylindrical coordinates
      x = sqrt(cx*cx + cy*cy);
      y = theta();
      z = cz;
      break;
    }
  ostr << '(' << Gform(coord::PTform, x) << ", "
              << Gform(coord::PTform, y) << ", "
              << Gform(coord::PTform, z) << ')';
  return ostr;
  }

ostream& operator<< (ostream& out, const coord& pt) { return pt.print(out); }

        // Input                istr : Input stream
        //                      pt   : Coordinate
        // Output               istr : Coordinate set from input
        // Bug                       : No error handling

istream& operator >> (istream& istr, coord& pt)
  { return istr >> pt.cx >> pt.cy >> pt.cz; }


// ____________________________________________________________________________
// J                  CLASS COORDINATE CONVERSION FUNCTIONS
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

//-----------------------------------------------------------------------------
//               Cartesian <---------> Spherical Conversion
//-----------------------------------------------------------------------------

coord coord::Cart2Sph(int rad) const
  {
  if(rad) return coord(Rad(), theta(), phi());
  return coord(Rad(), theta()*RAD2DEG, phi()*RAD2DEG);
  }

coord coord::Sph2Cart(int rad) const
  {
  double R     = cx;
  double Theta = cy;
  double Phi   = cz;
  if(!rad)
    {
    Theta *= DEG2RAD;
    Phi   *= DEG2RAD;
    }
  double x = R*sin(Theta)*cos(Phi);            // Convert to Cartesian
  double y = R*sin(Theta)*sin(Phi);            // from spherical
  double z = R*cos(Theta);
  return coord(x,y,z);
  }

//-----------------------------------------------------------------------------
//               Cartesian <---------> Cylindrical Conversion
//-----------------------------------------------------------------------------

coord coord::Cart2Cyl(int rad) const
  {
  double R = sqrt(cx*cx + cy*cy);
  double Theta;
  if(cx) Theta = atan(cy/cx);
  else
    {
    if(cy>0)      Theta = PI/2.0;
    else if(cy<0) Theta = 3.0*PI/2.0;
    else          Theta = 0.0;
    }
  if(rad) return coord(R, Theta, cz);
  return coord(R, Theta*RAD2DEG, cz);
  }

coord coord::Cyl2Cart(int rad) const
  {
  double R     = cx;
  double Theta = cy;
  double Z     = cz;
  if(!rad) Theta *= DEG2RAD;
  double X = R*cos(Theta);			// Convert to Cartesian
  double Y = R*sin(Theta);			// from cylindrical
  return coord(X,Y,Z);
  }

//-----------------------------------------------------------------------------
//               Spherical <---------> Cylindrical Conversion
//-----------------------------------------------------------------------------

coord coord::Sph2Cyl(int rad) const
  {
  coord pt = Sph2Cart(rad);
  return pt.Cart2Cyl(rad);
  }

coord coord::Cyl2Sph(int rad) const
  {
  coord pt = Cyl2Cart(rad);
  return pt.Cart2Sph(rad);
  }

// ____________________________________________________________________________
// L                  CLASS COORDINATE AUXILIARY FUNCTIONS
// ____________________________________________________________________________

coord coord::getDefCoord() { return DefCoord; }
 
        // Input		pt	: A coordinate point (this) 
        // Output               pt1     : The default coordinate

void coord::setDefCoord(const coord& dpt) { DefCoord = dpt; }
 
        // Input                pt      : A coordinate point (this)
        //                      dpt     : An additional coordinate
        // Output               void    : Default coordinate set to dpt
 
// ____________________________________________________________________________
// M                  CLASS COORDINATE COMPARISON FUNCTIONS
// ____________________________________________________________________________

        // Input                pt      : A coordinate point (this)
        //                      pt2	: Another coordinate
        // Output ==            T/F     : TRUE if pt2 and this are equal
        // Output !=            T/F     : FALSE if pt2 and this are equal

void coord::SetCutoff(double co)
  {  OrdCutoff = (co == -1.0)?1.e-13:fabs(co); }

bool coord::operator==(const coord& pt2) const
  {
  if(fabs(cx - pt2.cx) > OrdCutoff) return false;
  if(fabs(cy - pt2.cy) > OrdCutoff) return false;
  if(fabs(cz - pt2.cz) > OrdCutoff) return false;
  return true;
  }

bool coord::operator!=(const coord& pt2) const
  {
  if(fabs(cx - pt2.cx) > OrdCutoff) return true;
  if(fabs(cy - pt2.cy) > OrdCutoff) return true;
  if(fabs(cz - pt2.cz) > OrdCutoff) return true;
  return false;
  }

bool coord::operator>(const coord& pt) const
  { return (Rad() > pt.Rad()); }

bool coord::operator<(const coord& pt) const
  { return (Rad() < pt.Rad()); }

#endif						// coord.cc

