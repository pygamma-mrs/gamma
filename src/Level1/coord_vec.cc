/* coord_vec.cc *************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Coordinate Vector                      Implementation           **
**                                                                      **
**      Copyright (c) 1990, 1991, 1992                                  **
**      Scott Smith                                                     **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zuerich / Switzerland                                      **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  Class coordinate vector provides a facile means to manipulate a     **
**  general list of coordinates.  Eachvariable of type coord_vec        **
**  consists of a vector of coordinates in three dimensional space      **
**  and the number of points the vector contains.  The entire vector    **
**  can be manipulated as a whole(rotated, tranlated, etc.) in the same **
**  manner that a single  coordinate point (class coord) may be         **
**  manipulated.                                                        **
**                                                                      **
*************************************************************************/

#ifndef   Gcoord_vec_cc_		// Is this file already included?
#  define Gcoord_vec_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Level1/coord_vec.h>		// Include the interface
#include <Level1/coord.h>		// Include single coordinates
#include <Basics/Gconstants.h>		// This tells us what HUGE_VAL is
#include <Basics/Gutils.h>		// So we can use query_parameter
#include <Basics/ParamSet.h>		// We know parameter sets
#include <Basics/StringCut.h>		// Know string manipulations
#include <string>			// Know aobut libstdc++ strings
#include <list>				// Know about libstdc++ STL lists
#include <cmath>			// Include HUGE_VAL_VAL

using std::list;			// Using libstdc++ lists
using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                CLASS COORDINATE VECTOR ERROR HANDLING
// ____________________________________________________________________________

/* These functions handle all error messages for class coord_vec.  They tie
   into the generic GAMMA error messaging format.

              Input             cvec    : A coordinate vector
                                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  

void coord_vec::CVerror(int eidx, int noret) const
  {
  string hdr("Coordinate Vector");
  string msg;
  switch(eidx)
    {
    case 0: GAMMAerror(hdr, 0, noret); break;   // Program Aborting        (0)
    case 3: GAMMAerror(hdr, 3, noret); break;   // !Construct From Pset    (3)
    case 5: GAMMAerror(hdr, 5, noret); break;   // !Write To Pset          (5)
    case 6: msg = "Error Accessing Coordinate in Coordinate Vector";
            GAMMAerror(hdr, msg, noret); break; //                         (6)
    case 7: msg = "Cannot Read Coordinates from Parameter File";
            GAMMAerror(hdr, msg, noret); break; //                         (7)
    case 9: GAMMAerror(hdr, 9, noret); break;   // Construction Problems   (9)
    case 10: msg = "Negative Number of Points in Coordinate Vector";
             GAMMAerror(hdr, msg, noret); break; //                        (10)
    case 20: msg = "Set Parameter NCoords To Indicate Total # Coordinates";
             GAMMAerror(hdr, msg, noret); break; //                        (20)
    case 22: msg = "Problems Writing Coordinate Vector to Output FileStream";
             GAMMAerror(hdr, msg, noret); break; //                        (22)
    case 23: msg = "Cannot Output Coordinate Vector Parameters";
             GAMMAerror(hdr, msg, noret); break; //                        (23)
    case 26: msg = "Cannot Determine The Number of Coordinates";
             GAMMAerror(hdr, msg, noret); break; //                        (26)
    case 27: msg = "Please Set Either NCoords or NSpins Parameter";
             GAMMAerror(hdr, msg, noret); break; //                        (27)
    case 68: msg = "Coordinate Index Out of Coordinate Vector Range";
             GAMMAerror(hdr, msg, noret); break; //                        (68)
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

           
void coord_vec::CVerror(int eidx, const string& pname, int noret) const
  {      
  string hdr("Coordinate Vector");
  switch(eidx)
    {
    case 1:  GAMMAerror(hdr, 1, pname, noret);  break; // File Problems  (1)
    case 2:  GAMMAerror(hdr, 2, pname, noret);  break; // !Read Par      (2)
    default: GAMMAerror(hdr, -1, pname, noret); break; // Unknown Error  (-1)
    }
  }  

volatile void coord_vec::CVfatality(int eidx) const
  {
  CVerror(eidx, 1);				// Output error message
  if(eidx) CVerror(0);                          // State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                  COORDINATE VECTOR SETUP FUNCTIONS
// ____________________________________________________________________________
 
/* These functions set up specific aspects of a coordinate vector.  Since
   they make assumptions about the order in which the vector is set up, the
   functions MUST be private, their misuse may make an inconsistent vector.  */  


int coord_vec::SetNPoints(const ParameterSet& pset, int warn)

	// Input		cvec	: A Coordinate vector(this)
	// 			pset	: A parameter set
        //			warn	: Warning level
        //                                      0 = no warning
        //                                      1 = warning, non-fatal
        //                                     >1 = fatal
	// Output		none 	: Coord. vector # of points is set
	//				  from parameters in pset
	// Note				: If # points read, returns TRUE
	//				  whether or not points actually exist

/* ----------------------------------------------------------------------
-           Try & Read The Number Of Points In The Coordinate Vector 	-
-                                                                       -
-  Only two parameter names are currently acceptable:			-
-                                                                       -
-                         NCoords and NSpins				-
-                                                                       -
- The former will be read preferentially.  Both are read as integers.	-
-                                                                       -
-           Ncoords    (0) : 4   - Number of coordinates		-
-                                                                       -
- It neither NCoords or NSpins are found, it will attempt to read 	-
- individual coordinates (e.g. Coord(0)) beginning at index 0. If	-
- found it will set the number at the first one missing.		-
-                                                                       -
-----------------------------------------------------------------------*/

  {
  string pstate;
  int npts;
  ParameterSet::const_iterator item;		// Pix in parameter list
  string pname = string("NCoords");		// # coordinates par name
  item = pset.seek(pname);			// Try & find NCoords par
  if(item != pset.end())			// Retrieve the number of coordinates
    {
    (*item).parse(pname,npts,pstate);
    *this = coord_vec(npts);		
    return 1;
    }

//    We Can't Find NCoords, Look For NSpins As An Alternative Parameter Name

  pname = string("NSpins"); 
  item = pset.seek(pname);		// Pix in parameter list for NSpins
  if(item != pset.end())		// Retrieve the number of spins
    {
    (*item).parse(pname,npts,pstate);
    *this = coord_vec(npts);		
    return 1;
    }

// Can't Find NCoords Or NSpins. This Will Evoke Non-Fatal Warnings If Desired.
//       Then We'll Try & Just Count The Coordinates In The Parameter

  if(warn)
    {
    CVerror(26,0);			// Can't read NCoords
    CVerror(27);			// Please set NCoords!
    }
  coord pt;				// A dummy coordinate
  int nc = 0;				// A dummy index, start at 0
  while(pt.read(pset,nc,0))		// Try and read Coord(nc)
  nc++;					// until cant read one
  if(!nc)				// If cant find any coordinates
    {					// then we must issue warnings
    if(warn) 				// & perhaps die as desired!
      {
      CVerror(7, 1);			// Can't read any coordinates
      if(warn > 1) CVfatality(3);	// Can't make any vector, quit
      else         CVerror(3);		// Can't make cvec from parameters
      }
    return 0;
    }
  if(warn) CVerror(20);			// We've counted the coordinates
  *this = coord_vec(nc);		// Please set NCoords Or NSpins
  return 1;
  }


int coord_vec::SetCoords(const ParameterSet& pset, int warn)

	// Input		cvec	: A Coordinate vector(this)
	// 			pset	: A parameter set
        //                      warn    : Warning level, missing coord
	//				       -1 = warning, def 0 (DEF)
	//					0 = no warning, def 0
	//					1 = warning, def big
	//				       >2 = fatal

        //                                      0 = no warning
        //                                      1 = warning, non-fatal
        //                                     >1 = fatal

	// Output		none	: Coordinate vector points are
	//				  set from parameters in pset
	// Note				: Individual coordinates are 
	//				  set by class coord.
	// Note				: Pre-Set npts, # points read
	//				  and assumes array dimensioned

/* ----------------------------------------------------------------------
-           Try & Read The Individual Coordinates Of The  Vector 	-
-                                                                       -
- These will all be of the format, and individual ordinates are read	-
-                                                                       -
-           Coord(#)  (3) : (x,y,z)   - Coordinate x, y, z values	-
-                                                                       -
- as well as other formats supported by class coord.			-
-                                                                       -
-----------------------------------------------------------------------*/

  {
  coord pt;				// Working single coordinate
  int wl=0;				// Class coord warning level
  if(warn) wl=1;
  string pname = string("Coord(");	// For warnings only
  int count = 0;			// Counter of actual points read
  for(int i=0; i<Npts; i++)		// Retrieve the coordinates
    {					// All coordinates must be present
    if(pt.read(pset,i,wl))		// Note that this uses class coord!
      {
      put(pt,i);			// Store the point if found
      count++;				// Increment point counter
      }
    else
      {					// Can't find the coordinate
      if(warn)				// Output warnings if needed
        {
        pname += Gdec(i) + string(")");	//	Parameter name
        CVerror(2, pname, 1);		// 	Can't read parameter
        if(warn>1) CVfatality(3);	// 	Die if fatal error
        else std::cout << "\n";
        }
      put(pt,i);			// Store (default) point?
      }
    }
  if(!count)				// If we haven't read any points
    {
    if(Npts) delete [] Pts;		//   Delete all of the coordinates
    Pts = NULL;				//   Set the pointer to NULL
    Npts = 0;				//   Set the # of points to 0
    }
  return count;
  }

// ____________________________________________________________________________
// iii                 COORDINATE VECTOR CHECKING FUNCTIONS
// ____________________________________________________________________________

void coord_vec::check(int index) const
  {
  if(index<0 || index>=Npts)	// Check index range
    {
    CVerror(6, 1);	// Error accessing coordinate
    CVfatality(68);	// Coordinate access beyond vector range
    }
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A             CLASS COORDINATE VECTOR CONSTRUCTORS/DESTRUCTOR
// ____________________________________________________________________________

/* These constructors set up a new coordinate vector.  There are several ways
   to make a coordinate vector as listed below:

          Input Arguments                       Resulting Coordinate
   --------------------------        -----------------------------------------
                -                    Empty coordinate vector (no points)
             int pts                 Coordinate vector with pts points
          coord_vec cvec             Coordinate vector identical to cvec
     ParameterSet, idx, warn         Coord. vector w/ prefix i read from pset
             SinglePar               Coordinate set from parameter
 
   Setting a coordinate vector from a parameter set may fail if the parameters 
   that define coordinates don't exist. The flag warn dictates what happens
   upon a failure: 0=no warnings, 1=non-fatal warnings, 2=fatal warnings
   Individual coordinate (see class coord) parameters assume a parameter of
   type 3 with value stored as a string "( #, #, # )". All I/O between GAMMA
   coordinates and single parameters must follow this format!              */

coord_vec::coord_vec ( ) { Npts = 0; Pts = NULL; }
coord_vec::coord_vec(int pts)
  {
  if(pts<0)				// Check number of points
    {
    CVerror(10,1);			// Number of points is negative
    CVfatality(9);			// Error during construction
    }
  Npts = pts;				// Set the number of points
  if(Npts) Pts = new coord[pts]; 	// Set up array if >0 points
  else     Pts = NULL; 			// else point vector is NULL
  }

coord_vec::coord_vec(const coord_vec &cvec1)
  {
  Npts = cvec1.Npts;		// Copy number of points
  if(Npts)
    {				// Copy all points or vector
    Pts = new coord[Npts];	// array pt to NULL if no points
    for(int i=0; i<Npts; i++)
      Pts[i] = cvec1.Pts[i];
    }
  else Pts = NULL;
  }

coord_vec::coord_vec(const row_vector& X, const row_vector& Y, const row_vector& Z)
  {
  Npts = X.size();				// Length of X (Y & Z)
  if(!Npts) { Pts = NULL; return; }
  Pts  = new coord[Npts];
  for(int i=0; i<Npts; i++)
    Pts[i] = coord(X.getRe(i), Y.getRe(i), Z.getRe(i));
  }

coord_vec::coord_vec(const ParameterSet& pset, int idx, int warn)
  {
  Npts = 0;					// Start with no points
  Pts = NULL;					// and a NULL vector
  ParameterSet subpset;				// Copy parameter set
  if(idx != -1) subpset = pset.strip(idx);	// to glean out those for
  else          subpset = pset;			// the specified index

  if(!SetNPoints(subpset))			// Attempt to deterine how
    {						// many points are defined
    if(warn)					// If not able to, output
      { 					// some warnings if desired
      CVerror(3,1);				// Can't set cvec from pset
      if(warn>1) CVfatality(9);			// Error during construction
      else       CVerror(9);
      }
    return;
    }

  if(!SetCoords(subpset, warn))			// set Npts and allocate vector
    {						// Now we try & get points
    if(warn)					// If not able to, output
      { 					// some warnings if desired
      CVerror(7,1);				// Can't set cvec from pset
      if(warn>1) CVfatality(9);			// Error during construction
      else       CVerror(9);
      }
    }
  return;
  }


coord_vec::~coord_vec()
  {
  if(Npts) delete [] Pts;	// Delete all of the coordinates
  Pts = NULL;			// Set the pointer to NULL
  }


coord_vec& coord_vec::operator= (const coord_vec &cvec1)

	// Input		cvec1 : Coordinate vector(this)
	// 			cvec  : Coordinate vector
	// Output		none  : cvec1 vector identical to cvec

{
  if(Npts)			// Delete previous points if existing
    {
    delete [] Pts;
    Pts = NULL;
    Npts = 0;
    }
  if(cvec1.Npts)		// Copy number of points and points
    {
    Npts = cvec1.Npts;
    Pts = new coord[Npts];
    for(int i=0; i<Npts; i++)
      Pts[i] = cvec1.Pts[i];
    }

  return (*this);
}


// ____________________________________________________________________________
// B                  COORDINATE VECTOR ROTATION FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//            Typical Rotations Active On Cartesian Coordinates
// ----------------------------------------------------------------------------

	// Input		cvec : Coordinate vector(this)
        //                      theta: Spherical angle (degrees down from +z)
        //                      phi  : Spherical angle (degrees over from +x)
        //                      rad  : Flag for angle in degrees or radians
        // Return               rotcv: The vector rotated about the x/y-axis 
	//			       or the z-axis by the input angle
        // Note                      : The reference axes are static!  Thus
        //                             multiple rotations should be cumulative

coord_vec coord_vec::xrotate(double theta, int rad) const
  { return rotate(UnitX.Rx(theta,rad)); }

coord_vec coord_vec::yrotate(double theta, int rad) const
  { return rotate(UnitY.Ry(theta,rad)); }

coord_vec coord_vec::zrotate(double phi, int rad) const
  { return rotate(UnitZ.Rz(phi,rad)); }

coord_vec coord_vec::rotate(const matrix& aRmx) const

	// Input		cvec : Coordinate vector(this)
        //                    aRmx  : 3x3 Rotation matrix
        //                      rad  : Flag for angle in degrees or radians
        // Return               rotcv: The vector multiplied by aRmx

  {
  coord_vec rotcv(Npts);
  for(int i=0; i<Npts; i++)			// Rotate all the coords.
    rotcv.Pts[i] = aRmx * Pts[i];
  return rotcv;					// Return rotcv
  }


// ----------------------------------------------------------------------------
//           Euler Rotation Matrices Active On Cartesian Coordinates
// ----------------------------------------------------------------------------

coord_vec coord_vec::rotate(double alpha, double beta, double gamma) const

	// Input		cvec : Coordinate vector(this)
	// 			alpha: Euler angle (radians)
	// 			beta : Euler angle (radians)
	// 			gamma: Euler angle (radians)
	// Return		rotpt: Coordinate vector rotated
	//			       by input Euler angles

  {
  coord_vec rotcv(Npts);
  matrix rotmx = Rmx(alpha, beta, gamma);	// Set rotation matrix
  for(int i=0; i<Npts; i++)			// Rotate all the coords.
    rotcv.Pts[i] = rotmx * Pts[i];
  return rotcv;					// Return rotcv
  }


coord_vec coord_vec::rotate(const coord &EA) const

	// Input		cvec : Coordinate vector(this)
	// 			EA   : Coordinate containing Euler angles
	// Return		rotcv: Coordinate vector rotated
	//			       by input Euler angles
	// Note			     : Angles input in radians

  {
  coord_vec rotcv(Npts);
  matrix rotmx = Rmx(EA.x(), EA.y(), EA.z());	// Set rotation matrix
  for(int i=0; i<Npts; i++)			// Rotate all the coords.
    rotcv.Pts[i] = rotmx * Pts[i];
  return rotcv;					// Return rotcv
  }


void coord_vec::rotate_ip(double alpha, double beta, double gamma)

	// Input		cvec : Coordinate vector(this)
	// 			alpha: Euler angle (radians)
	// 			beta : Euler angle (radians)
	// 			gamma: Euler angle (radians)
	// Return		cvec : Coordinate vector rotated
	//			       by input Euler angles

  {
  matrix rotmx = Rmx(alpha, beta, gamma);	// Set rotation matrix
  for(int i=0; i<Npts; i++)			// Rotate all the coords.
    Pts[i] = rotmx * Pts[i];
  return;
  }


void coord_vec::rotate_ip(const coord &EA)

	// Input		cvec : Coordinate vector(this)
	// 			EA   : Coordinate containing Euler angles
	// Return		cvec : Coordinate vector rotated
	//			       by input Euler angles
	// Note			     : Angles input in radians

  {
  matrix rotmx = Rmx(EA.x(), EA.y(), EA.z());	// Set rotation matrix
  for(int i=0; i<Npts; i++)			// Rotate all the coords.
    Pts[i] = rotmx * Pts[i];
  return;
  }


// ____________________________________________________________________________
// C                  COORDINATE TRANSLATION FUNCTIONS
// ____________________________________________________________________________


coord_vec coord_vec::translate(double delx, double dely, double delz) const

	// Input		cvec    : Coordinate vector (this)
	// 			delx    : Change in x coordinate
	// 			dely    : Change in y coordinate
	// 			delz    : Change in z coordinate
	// Return		tcvec   : A coordinate vector which is cvec
	//				  translated by delx,dely, & delz
	// Note				: dely and delz default to zero
  {
  coord_vec tcvec(*this);
  for(int i=0; i<Npts; i++)
    (tcvec.Pts[i]).translate_ip(delx,dely,delz);
  return tcvec;
  }


coord_vec coord_vec::translate(const coord &del) const

	// Input		cvec    : Coordinate vector (this)
	// 			del     : Change in x,y,&z coordinates
	// Return		tcvec   : A coordinate vector which is cvec
	//				  translated by del
  {
  coord_vec tcvec(*this);
  for(int i=0; i<Npts; i++)
    (tcvec.Pts[i]).translate_ip(del);
  return tcvec;
  }


void coord_vec::translate_ip(double delx, double dely, double delz)

	// Input		cvec    : Coordinate vector (this)
	// 			delx    : Change in x coordinate
	// 			dely    : Change in y coordinate
	// 			delz    : Change in z coordinate
	// Return		cvec    : Coordinate vector is
	//				  translated by delx,dely, & delz
	// Note				: dely and delz default to zero

  { for(int i=0; i<Npts; i++) (Pts[i]).translate_ip(delx,dely,delz); }


void coord_vec::translate_ip(const coord &del)

	// Input		cvec    : Coordinate vector (this)
	// 			del     : Change in x,y,&z coordinates
	// Return		cvec    : Coordinate vector is
	//				  translated by del

  { for(int i=0; i<Npts; i++) (Pts[i]).translate_ip(del); }


coord_vec coord_vec::trans_x(double delx) const

	// Input		cvec    : Coordinate vector (this)
	// 			delx    : Change in x coordinate
	// Return		tcvec   : A coordinate vector which is cvec
	//				  translated by delx

  {
  coord_vec tcvec(*this);
  for(int i=0; i<Npts; i++)
   (tcvec.Pts[i]).trans_x_ip(delx);
  return tcvec;
  }


void coord_vec::trans_x_ip(double delx)

	// Input		cvec    : Coordinate vector (this)
	// 			delx    : Change in x coordinate
	// Return		cvec	: Coordinate vector is
	//				  translated along x axis by delx

  { for(int i=0; i<Npts; i++) (Pts[i]).trans_x_ip(delx); }


coord_vec coord_vec::trans_y(double dely) const

	// Input		cvec    : Coordinate vector (this)
	// 			dely    : Change in y coordinate
	// Return		tcvec   : A coordinate vector which is cvec
	//				  translated by dely

  { 
  coord_vec tcvec(*this);
  for(int i=0; i<Npts; i++)
   (tcvec.Pts[i]).trans_y_ip(dely);
  return tcvec;
  }


void coord_vec::trans_y_ip(double dely)

	// Input		cvec    : Coordinate vector (this)
	// 			dely    : Change in y coordinate
	// Return		cvec	: Coordinate vector is
	//				  translated along y axis by dely

  { for(int i=0; i<Npts; i++) (Pts[i]).trans_y_ip(dely); }


coord_vec coord_vec::trans_z(double delz) const

	// Input		cvec    : Coordinate vector (this)
	// 			delz    : Change in z coordinate
	// Return		tcvec   : A coordinate vector which is cvec
	//				  translated by delz

  {
  coord_vec tcvec(*this);
  for(int i=0; i<Npts; i++)
   (tcvec.Pts[i]).trans_z_ip(delz);
  return tcvec;
  }


void coord_vec::trans_z_ip(double delz)

	// Input		cvec    : Coordinate vector (this)
	// 			delz    : Change in z coordinate
	// Return		cvec	: Coordinate vector is
	//				  translated along z axis by delz

  { for(int i=0; i<Npts; i++) (Pts[i]).trans_z_ip(delz); }

// ____________________________________________________________________________
// D                  COORDINATE PROJECTION FUNCTIONS
// ____________________________________________________________________________


row_vector coord_vec::project(int projx, int projy) const

	// Input		cvec    : Coordinate vector (this)
	// 			projx   : Axis flag to project onto x
	// 			projy   : Axis to project onto y
	//				  1 = x axis; -1 = -x axis
	//				  2 = y axis; -2 = -y axis
	//				  3 = x axis; -3 = -z axis
	// Return		vx      : Row vector having cvec (3D)
	//		  		  projected into it (2D) via axes
	//				  mapping specified by projx & projy

  {
  row_vector vx((*this).Npts);
  double xx,yy;
  for(int i=0; i<Npts; i++)
    {
    switch(projx)
      {
      case 1:			// 3D x-axis into 2D x-axis
      default:
        xx = Pts[i].x();
        break;
      case -1:			// 3D x-axis into 2D negative x-axis
        xx = -Pts[i].x();
        break;
      case 2:			// 3D y-axis into 2D x-axis
        xx = Pts[i].y();
        break;
      case -2:			// 3D y-axis into 2D negative x-axis
        xx = -Pts[i].y();
        break;
      case 3:			// 3D z-axis into 2D x-axis
        xx = Pts[i].z();
        break;
      case -3:			// 3D z-axis into 2D negative x-axis
        xx = -Pts[i].z();
        break;
      }

    switch(projy)
      {
      case 1:			// 3D x-axis into 2D y-axis
        yy = Pts[i].x();
        break;
      case -1:			// 3D x-axis into 2D negative y-axis
        yy = -Pts[i].x();
        break;
      case 2:			// 3D y-axis into 2D y-axis
      default:
        yy = Pts[i].y();
        break;
      case -2:			// 3D y-axis into 2D negative y-axis
        yy = -Pts[i].y();
        break;
      case 3:			// 3D z-axis into 2D y-axis
        yy = Pts[i].z();
        break;
      case -3:			// 3D z-axis into 2D negative y-axis
        yy = -Pts[i].z();
        break;
      }
    vx.put(complex(xx,yy), i);
    }
  return vx;
  }


// ____________________________________________________________________________
// E                    COORDINATE VECTOR OPERATORS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                      COORDINATE VECTOR WITH SCALAR
// ----------------------------------------------------------------------------

/*  Function/Operator                         Result
    -----------------   -------------------------------------------------------
         cv*r           Multiplies all coordinates in vector by r
         r*cv           Multiplies all coordinates in vector by r
        cv*=r           Multiplies all coordinates in vector by r
         cv/r           Multiplies all coordinates in vector by 1/r
        cv/=r           Multiplies all coordinates in vector by 1/r          */
	
coord_vec operator * (double r, const coord_vec& cvec1)
  { coord_vec cvec(cvec1); cvec *= r; return cvec; }

coord_vec coord_vec::operator * (double r) const
  { coord_vec cvec(*this); cvec*=r; return cvec; }

coord_vec& coord_vec::operator *= (double r)
  { for(int i=0; i<Npts; i++) Pts[i] *= r; return (*this); }

coord_vec coord_vec::operator / (double r) const
  { coord_vec cvec(*this); cvec/=r; return cvec; }

coord_vec& coord_vec::operator /= (double r)
  { for(int i=0; i<Npts; i++) Pts[i] /= r; return (*this); }

// ----------------------------------------------------------------------------
//                     COORDINATE VECTOR WITH COORDINATE
// ----------------------------------------------------------------------------

/*  Function/Operator                         Result
    -----------------   -------------------------------------------------------
         cv+pt          Adds      point p to all coordinates in vector
         cv+=pt         Adds      point p to all coordinates in vector
         cv-pt          Subtracts point p from all coordinates in vector 
         cv-=pt         Subtracts point p from all coordinates in vector     */

coord_vec coord_vec::operator + (const coord& pt) const
  { coord_vec cvec(*this); cvec+=pt; return cvec; }

void coord_vec::operator += (const coord& pt)
  { for(int i=0; i<Npts; i++) Pts[i] += pt; }

coord_vec coord_vec::operator - (const coord& pt) const
  { coord_vec cvec(*this); cvec-=pt; return cvec; }

void coord_vec::operator -= (const coord& pt)
  { for(int i=0; i<Npts; i++) Pts[i] -= pt; }

// ----------------------------------------------------------------------------
//                  COORDINATE VECTOR WITH COORDINATE VECTOR
// ----------------------------------------------------------------------------

// sosi should add a bound check here

coord_vec coord_vec::operator + (const coord_vec& cv) const
  { coord_vec cvec(*this); cvec+=cv; return cvec; }

coord_vec& coord_vec::operator += (const coord_vec& cv)
  { for(int i=0; i<Npts; i++) Pts[i] += cv.Pts[i]; return (*this); }

coord_vec coord_vec::operator - (const coord_vec& cv) const
  { coord_vec cvec(*this); cvec-=cv; return cvec; }

coord_vec& coord_vec::operator -= (const coord_vec& cv)
  { for(int i=0; i<Npts; i++) Pts[i] -= cv.Pts[i]; return (*this); }

// ____________________________________________________________________________
// F                  COORDINATE VECTOR GENERAL FUNCTIONS
// ____________________________________________________________________________

/*  Function/Operator                         Result
    -----------------   -------------------------------------------------------
         size           Returns number of coordinates in the vector
         max_x          Returns x-ordinate with maximum magnitude
         max_y          Returns y-ordinate with maximum magnitude
         max_z          Returns z-ordinate with maximum magnitude
         maxima         Returns max_x, max_y, and max_z in a coordinate
         maxima         Returns max_x, max_y, and max_z individually
         max_R		Returns distance of point farthest from origin
         max_R		Returns distance & index of point farthest from origin
         vectors        Returns coordinates of vectors connecting all points
         cv-=pt         Subtracts point p from all coordinates in vector     */

int coord_vec::size() const { return Npts; }

double coord_vec::max_x() const
  {
  double xx, maxx;
  maxx = -HUGE_VAL;
  for(int i=0; i<Npts; i++)
    {
    xx = Pts[i].x();
    if(fabs(xx) > maxx)
      maxx = xx;
    }
  return maxx;
  }

double coord_vec::max_y() const
  {
  double yy, maxy;
  maxy = -HUGE_VAL;
  for(int i=0; i<Npts; i++)
    {
    yy = Pts[i].y();
    if(fabs(yy) > maxy)
      maxy = yy;
    }
  return maxy;
  }

double coord_vec::max_z() const
  {
  double zz, maxz;
  maxz = -HUGE_VAL;
  for(int i=0; i<Npts; i++)
    {
    zz = Pts[i].z();
    if(fabs(zz) > maxz)
      maxz = zz;
    }
  return maxz;
  }


coord coord_vec::maxima() const
  { return coord(max_x(), max_y(), max_z()); }

void coord_vec::maxima(double &x, double &y, double &z) const
  { x = max_x(); y = max_y(); z = max_z(); }

double coord_vec::max_R() const
  {
  double RR, maxR;
  maxR = -HUGE_VAL;
  for(int i=0; i<Npts; i++)
    {
    RR = Pts[i].Rad();
    if(fabs(RR) > maxR)
      maxR = RR;
    }
  return maxR;
  }

void coord_vec::max_R(int &maxi, double &maxR) const
  {
  double RR;
  int i;
  maxR = -HUGE_VAL;
  for(i=0; i<Npts; i++)
    {
    RR = Pts[i].Rad();
    if(fabs(RR) > maxR)
      {
      maxR = RR;
      maxi = i;
      }
    }
  return;
  }

coord_vec coord_vec::vectors() const
  {
  coord_vec cvect((Npts*(Npts-1))/2);
  int k=0;
  coord pt;
  for(int i=0; i<Npts-1; i++)
    for(int j=i+1; j<Npts; j++)
    {
    pt = cdvect((*this)(i), (*this)(j));
    cvect.put(pt,k);
    k++;
    }
  return cvect;
  }


coord_vec coord_vec::vectors_f() const

	// Input		cvec  : A coordinate vector(this)
	// Return		cvect : Coordinates of the vectors
	//				connecting pairs of points

  {
  coord_vec cvect(Npts*Npts);
  int k=0;
  coord pt;
  for(int i=0; i<Npts; i++)
    for(int j=0; j<Npts; j++)
    {
    pt = cdvect((*this)(i), (*this)(j));
    cvect.put(pt,k);
    k++;
    }
  return cvect;
  }


matrix coord_vec::distances(int Angs) const

	// Input		cvec  : A coordinate vector(this)
	//			Angs  : Flag for output in Angstroms
	// Return		dmx   : A distance matrix between
	//				the points of cvec
	// Note			      : Flag Angs defaults to 0

  {
  matrix dmx(Npts, Npts, complex0);
  double dist;
  double Af =1.0;
  if(Angs) Af = 1.0e10;
  int i=0, j=0;
  for(i=0; i<Npts-1; i++)
    for(j=i+1; j<Npts; j++)
      {
      dist = Af*Rad(get(i), get(j));
      dmx.put(dist, i, j);
      dmx.put(dist, j, i);
      }
  return dmx;
  }


double coord_vec::distance(int pt1, int pt2, int Angs) const

	// Input		cvec  : A coordinate vector(this)
	//			pt1   : First coordinate point
	//			pt2   : Second coordinate point
	//			Angs  : Flag for output in Angstroms
	// Return		dist  : Distance between points pt1 & pt2
	// Note			      : Flag Angs defaults to 0

  {
  double   Af = 1.0;
  if(Angs) Af = 1.0e10;
  return Af*Rad(get(pt1), get(pt2));
  }


matrix coord_vec::thetas(int deg) const

	// Input		cvec  : A coordinate vector(this)
	//			deg   : Flag for output in degrees
	// Return		the_mx: Matrix of theta angles between
	//				the points of cvec
	// Note			      : Flag def defaults to 0

  {
  matrix the_mx(Npts, Npts, complex0);
  double Theta;
  double df =1.0;
  if(deg) df = 360./(2*3.14159);
  int i=0, j=0;
  for(i=0; i<Npts; i++)
    for(j=0; j<Npts; j++)
      {
      if(i != j)
        {
        Theta = theta(get(i), get(j));
        the_mx.put(Theta*df, i,j);
        }
      }
  return the_mx;
  }


matrix coord_vec::phis(int deg) const

	// Input		cvec  : A coordinate vector(this)
	//			deg   : Flag for output in degrees
	// Return		phi_mx: Matrix of phi angles between
	//				the points of cvec
	// Note			      : Flag deg defaults to 0

  {
  matrix phi_mx(Npts, Npts, complex0);
  double Phi;
  double df =1.0;
  if(deg) df = 360./(2*3.14159);
  for(int i=0; i<Npts; i++)
    for(int j=0; j<Npts; j++)
      if(i != j)
        {
        Phi = phi(get(i), get(j));
        phi_mx.put(Phi*df, i,j);
        }
  return phi_mx;
  }


// ____________________________________________________________________________
// G                      INDIVIDUAL COORDINATE ACCESS
// ____________________________________________________________________________

///Center Coordinate Access Functions

/*  Function/Operator                         Result
    -----------------   -------------------------------------------------------
          =             Current array set equal to input mx
        (int)           Access to coordinate point at position specified
     put(pt,int)        Assigns coordinate point at position specified
    put(x,y,z,int)      Assigns coordinate point at position specified
       get(int)         Copy of coordinate at position specified
    get_block(i,np)     Returns coordinate vector of size np starting at i
    put_block(i,cv)     Sets coordinates in vector starting at i from cv
        x(int)          Get first ordinate of coordinate at index 
        y(int)          Get seconds ordinate of coordinate at index 
        z(int)          Get third ordinate of coordinate at index            */

coord coord_vec::operator() (int index) const
  { check(index); return Pts[index]; }

void coord_vec::put(const coord &pt1, int index)
  { check(index); Pts[index] = pt1; }

void coord_vec::put(double x, double y, double z, int index)
  { check(index); Pts[index].xyz(x,y,z); }

coord  coord_vec::get(int index) const { check(index); return Pts[index]; }
double coord_vec::x(int index)   const { check(index); return Pts[index].x(); }
double coord_vec::y(int index)   const { check(index); return Pts[index].y(); }
double coord_vec::z(int index)   const { check(index); return Pts[index].z(); }

coord_vec coord_vec::get_block(int index, int npts) const
  {
  check(index); 				// Insure index in range
  check(index+npts-1); 				// Insure span in range
  coord_vec cv(npts);				// New coordinate vector
  for(int i=0; i<npts; i++)			// Copy requested coordinates
    cv.Pts[i] = Pts[index+i];
  return cv;
  }
  
void coord_vec::put_block(int index, const coord_vec& cv) const
  {
  check(index); 				// Insure index in range
  check(index+cv.Npts-1);			// Insure span in range
  for(int i=0; i<cv.Npts; i++)			// Copy requested coordinates
    Pts[index+i] = cv.Pts[i];
  }

// ____________________________________________________________________________
// H                      PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//      Functions To Make A Parameter Set From A Coordinate Vector
// ----------------------------------------------------------------------------

	// Input		cvec  : Coordinate vector
	//  			pset  : Parameter set
	//			idx   : Coordinate parameter name prefix #
	// Output ()		pset  : Parameter set with
	//			        only coordinates
	// Output +=            void  : Parameter set has coordinate added
        // Output Add           void    : Coordinate vector parameters are
        //                                are added ot the parameter set
        //                                with interaction index idx
        // Note                         : The parameter names & types
        //                                output here MUST match those used
        //                                in setting the coordinate vector
        //                                from parameters sets

coord_vec::operator ParameterSet( ) const
 { ParameterSet pset; pset += *this; return pset; }

void operator+= (ParameterSet& pset, const coord_vec &cvec)
  { cvec.PSetAdd(pset); }

void coord_vec::PSetAdd(ParameterSet& pset, int idx) const
  {
  string prefx;                                 // Parameter prefix
  if(idx != -1)                                 // Only use prefix if idx
    prefx = string("[")+Gdec(idx)+string("]");	// is NOT -1
  string pname;
  string pdata;
  string pstate;
  SinglePar par;
 
  pname = prefx + string("NCoords");		// Add the number of points
  pstate = string("Number of Coordinates");
  pdata = Gdec(Npts);
  par = SinglePar(pname,1,pdata,pstate);
  pset.push_back(par);

  pstate = string("Coordinate Point");		// Add Coordinates
  for(int i=0; i<Npts; i++)
    {
    pname = prefx + string("Coord(");
    pname += Gdec(i);
    pname += string(")");
    pdata = string("( ");
    par = (Pts[i]).param(pname,pstate);
    pset.push_back(par);
    }
  return;
  } 

// ----------------------------------------------------------------------------
//      Functions To Make A Coordinate Vector From A Parameter Set
// ----------------------------------------------------------------------------

	// Input		cvec     : A Coordinate vector(this)
	// 			pset     : A parameter set
	// Output		none	 : Coordinate vector filled with
	//				   coordinates in pset
	// Note				 : Two things are gleaned from the
	//				   parameter set from a coordinate
	//				   vector
	//				   1.) The number of points
	//				   2.) Coordinates for each point
	// Note				 : Functions which place a coordinate
	//				   vector into a parameter set must
	//				   contain the information read here
	// Note				 : This function is used by the read
	//				   functions as well!

int coord_vec::operator= (const ParameterSet& pset)
   {
   if(!SetNPoints(pset)) return 0;
   if(!SetCoords(pset))  return 0;
   return 1;
   }

// ----------------------------------------------------------------------------
//    Functions To Output Coordinate Vector To ASCII From A Parameter Set 
// ---------------------------------------------------------------------------- 


	// Input		cvec   	 : Coordinate vector
	//			filename : Output file name
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
	// Output		none 	 : Coordinate vector is written as a 
	//				   parameter set to file filename

int coord_vec::write(const string &filename, int idx, int warn) const
  {
  if(!Npts) return 1;			// Nothing if no coordinates
  std::ofstream ofstr(filename.c_str());// Open filename for output
  int w2 = 0;				// Added warning level
  if(warn) w2 = 1;
  if(!write(ofstr, idx, w2))		// If file bad then exit
    {
    if(warn)
      {
      CVerror(1, filename);		// Problems with file
      if(warn>1) CVfatality(5);		// Cant write to pset file
      }  
    return 0;
    }
  ofstr.close();                        // Close it now
  return 1;
  }



	// Input		cvec   	 : Coordinate vector
        //                      ofstr   : Output file stream
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Coordinate vector is written as a
        //                                parameter set to output filestream
        // Note                         : This depends on function PSetAdd!

int coord_vec::write(std::ofstream& ofstr, int idx, int warn) const
  {
  if(!Npts) return 1;			// Nothing if no coordinates
  if(!ofstr.good())                    // If file bad then exit
    {
    if(warn) CVerror(22, 1);		// Problems with file
    if(warn>1) CVfatality(23);		// Fatal
    }
  ParameterSet pset;                 // Declare a parameter set
  PSetAdd(pset, idx);                   // Add system to parameter set
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!pset.write(ofstr, w2))            // Use parameter set to write
    {
    if(warn)
      {
      CVerror(22, 1);			// Problems writing to filestream
      if(warn>1) CVfatality(23);	// Fatal error
      }
    return 0;
    }  
  return 1;
  }  

// ____________________________________________________________________________
// I                COORDINATE VECTOR AUXILIARY FUNCTIONS
// ____________________________________________________________________________

	// Input		cvec  : Coordinate vector (this)
        // Output               TF    : True if all coordinates in
	//				cvec are (0,0,0)

bool coord_vec::is_zero( ) const
   {
   for(int i=0; i<Npts; i++)		// Loop over all coords
    {					// and see if any have non-zero
         if(Pts[i].x()) return false;	// components
    else if(Pts[i].y()) return false;
    else if(Pts[i].z()) return false;
    }
   return true;
   }


// ____________________________________________________________________________
// J                    COORDINATE VECTOR INPUT FUNCTIONS
// ____________________________________________________________________________

 
string coord_vec::ask_read(int argc, char* argv[], int argn, int idx, int warn)
 
        // Input                cvec    : Dipolar interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
	//			idx	: Vector index
	//			warn	: Warning output level
        // Output               string  : Parameter argn of array argc is used
        //                                to supply a filename from which the
        //                                coordinate vector coords. are read
        //                                If argument argn is not in argv, the
        //                                user is asked to supply a filename
        //                                The set filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing known coord_vec params.
        // Note                         : The coordinate vector is modifed
 
  {
  string filename;				// Name of parameter file
  query_parameter(argc, argv, argn,             // Get filename from command
   "\n\tCoordinate Vector filename? ",filename);// Or ask for it
  read(filename, idx, warn);			// Read system from filename
  return filename;
  }



	// Input		cvec	: Coordinate vector
	// 			filein	: Input filename
	//			idx	: Vector index
	//			warn	: Warning output level
	//				      0: no warnings
	//				      1: non-fatal warnings
	//				     >2: fatal warnings
	// Output		none	: Spin system filled with
	//				  parameters read from file

	// Input		cvec    : Coordinate vector
        //                      pset    : A parameter set
	//			idx	: Vector index
	//			warn	: Warning output level
	//				      0: no warnings
	//				      1: non-fatal warnings
	//				     >1: fatal warnings
        // Output               none    : Coordinate vector filled with
        //                                coordinates in pset
        // Note                         : This function is used by the read
        //                                functions as well!

int coord_vec::read(const string &filein, int idx, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filein, 1))             // Read in pset from file
    {					// If we can't do that:
    if(warn>0)				//   Issue warnings if
      {					//   desired
      CVerror(1, filein, 1);		//   Problems with file filein
      if(warn>1) CVfatality(7);		//   Can't read from file filein
      }	
    return 0;				//   Failed to even read pset
    }
  return read(pset, idx, warn);		// Know pset, get coord_vec
  }


int coord_vec::read(const ParameterSet& pset, int idx, int warn)
  {  
  if(Npts)					// Delete existing vector
    {
    delete [] Pts;
    Pts  = NULL;
    Npts = 0;
    }
  ParameterSet  subpset;			// Copy parameter set
  if(idx != -1) subpset = pset.strip(idx);	// to glean out those for
  else          subpset = pset;			// the specified index
  if(!SetNPoints(subpset, warn)) return 0;	// Get # of points in vector
  if(!SetCoords(subpset, warn))  return 0;	// Fill in the points
  return 1;
  }  


// ____________________________________________________________________________
// K                  COORDINATE VECTOR OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//                   Formatted Output To Screen Or File
//-----------------------------------------------------------------------------

	// Input		cvec  : Coordinate vector (this)
        //                      out   : Output stream
 	//			units : Flag for input units
	//				0 = no units, print as is 
	//				1 = Angs., scale by 1.e10
	//				2 = nm,    scale by 1.e19
	//				3 = m,     scale by 1
        // Output               none  : Coordinate vector is sent
	//				to the output stream
	// Note			      : Flag units defaults to 0

ostream& coord_vec::printf(ostream& out, int units) const
  {
  double XYZ=0.0, fact=1.0;
  out << "\nCoordinates";
  if(units)
    {
    if(units==1)
      {
      fact = 1.e10;
      out << "\n(Angstroms)";
      }
    else if(units==2)
      {
      fact = 1.e9;
      out << "\n (nmeters)";
      }
    else if(units==3)
      out << "\n (meters)";
    }
  for(int i=0; i<3; i++) 
    {
    if(i==0) out << "\nX        :";
    else if(i==1) out << "\nY        :";
    else out << "\nZ        :";
//    out << setw(10) << setprecision(2);
    for (int j=0; j<Npts; j++)
      {
      switch(i)
        {
        case 0:				// First print the x coordinates
          XYZ = fact*x(j);
          break;
        case 1:				// Now print all y coordinates
          XYZ = fact*y(j);
          break;
        case 2:				// Last print all z coordinates
          XYZ = fact*z(j);
          break;
        }
      out << Gform("%10.2f", XYZ);
      }
    }
  return out;
  }


ostream& coord_vec::printCylindrical(ostream& ostr, double sf) const

	// Input		cvec  : Coordinate vector (this)
	//				assumed Cartesian {x,y,z}
        //                      ostr  : Output stream

  {
  for(int i=0; i<Npts; i++)
    ostr << "\n" << (sf*Pts[i]).Cart2Cyl(0);
  return ostr;
  }


ostream& coord_vec::printSpherical(ostream& ostr, double sf) const

	// Input		cvec  : Coordinate vector (this)
	//				assumed Cartesian {x,y,z}
        //                      ostr  : Output stream

  {
  for(int i=0; i<Npts; i++)
    ostr << "\n" << (sf*Pts[i]).Cart2Sph(0);
  return ostr;
  }


ostream& coord_vec::printCartesian(ostream& ostr, double sf) const

	// Input		cvec  : Coordinate vector (this)
	//				assumed Cartesian {x,y,z}
        //                      ostr  : Output stream

  {
  for(int i=0; i<Npts; i++) ostr << "\n" << Pts[i]*sf;
  return ostr;
  }


	// Input		cvec  : Coordinate vector (this)
        //                      ostr  : Output stream
        // Output               none  : Coordinate vector is sent
	//				to the output stream
	// Note			      : Output format of individual
	//				coordinates is set in class
	//				coord
    
ostream& coord_vec::print(ostream& ostr, int NC, int N) const
  {
  if(N < 1)  N  = Npts;				// Insure # points OK
  if(NC < 1) NC = 1;				// Insure # columns OK
  int cp = 0;
  string lst("\n");				// Line start string
  string cs("   ");				// Column spacing
  int ilen = (Gdec(N-1)).length();		// Width of index string
  int llen = coord::length()*NC			// Width of print line
           + cs.length()*(NC-1);
  if(llen/2 < 40)				// Add spacing to center
     lst += string(40-llen/2, ' ');		// if needed
  
  for(int i=0; i<N; i++)
    {
    if(!cp) ostr << lst;			// Either write line start
    else    ostr << cs;				// or space to new column
    ostr << Gdec(i, ilen) << ". ";		// Write the point index
    ostr << Pts[i];				// Write the point itself
    cp++;					// Update column count
    if(cp >= NC) cp = 0;			// Adjust if N columns output
    }
  return ostr;
  }

ostream &operator << (ostream &ostr, const coord_vec &cvec)
  { return cvec.print(ostr); }


//-----------------------------------------------------------------------------
//                           Binary Output To File
//-----------------------------------------------------------------------------


//void coord_vec::write(const string &filename) const

	// Input		cvec   	 : Coordinate vector
	//			filename : Output file name
	// Output		none 	 : Coordinate vector is written as a 
	//				   parameter set to file filename


// ____________________________________________________________________________
// L             CLASS COORDINATE VECTOR CONVERSION FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//               Cartesian <---------> Spherical Conversion
//-----------------------------------------------------------------------------

coord_vec coord_vec::Cart2Sph(int rad) const

	// Input		cvec    : Coordinate vector
        //                                assumed to be Cartesian {x, y, z }
        //                      rad     : Radians vs. Deg angle output
        // Output               cvec1   : Equivalent vector in Spherical
        //                                coordinates { R, theta, phi }

  {  
  coord_vec cvec1(Npts);
  for(int i=0; i<Npts; i++)
    cvec1.put(Pts[i].Cart2Sph(rad),i);
  return cvec1;
  }  


coord_vec coord_vec::Sph2Cart(int rad) const

        // Input                cvec    : Coordinate vector
        //                      rad     : Radians vs. Deg angle input
        // Output               pt1     : Equivalent vector in Spherical
        //                                coordinates { R, theta, phi }

  {  
  coord_vec cvec1(Npts);
  for(int i=0; i<Npts; i++)
    cvec1.put(Pts[i].Sph2Cart(rad),i);
  return cvec1;
  }  


//-----------------------------------------------------------------------------
//               Cartesian <---------> Cylindrical Conversion
//-----------------------------------------------------------------------------

coord_vec coord_vec::Cart2Cyl(int rad) const

        // Input                cvec    : Coordinate vector
        //                                assumed to be Cartesian {x, y, z }
        //                      rad     : Radians vs. Deg angle output
        // Output               pt1     : Equivalent vector in Cylindrical
        //                                coordinates { R, theta, z }

  {  
  coord_vec cvec1(Npts);
  for(int i=0; i<Npts; i++)
    cvec1.put(Pts[i].Cart2Cyl(rad),i);
  return cvec1;
  }  


coord_vec coord_vec::Cyl2Cart(int rad) const

        // Input                cvec    : Coordinate vector
        //                                assumed cylindrical { R, theta, z }
        //                      rad     : Radians vs. Deg angle input
        // Output               pt1     : Equivalent vector in Cartesian
        //                                coordinates { x, y, z }

  {  
  coord_vec cvec1(Npts);
  for(int i=0; i<Npts; i++)
    cvec1.put(Pts[i].Cyl2Cart(rad),i);
  return cvec1;
  }  

//-----------------------------------------------------------------------------
//               Spherical <---------> Cylindrical Conversion
//-----------------------------------------------------------------------------

coord_vec coord_vec::Sph2Cyl(int rad) const

        // Input                cvec    : Coordinate vector
        //                                assumed spherical { R, theta,phi }
        //                      rad     : Radians vs. Deg angle output
        // Output               pt1     : Equivalent vector in Cylindrical
        //                                coordinates { R, theta, z }

  {  
  coord_vec cvec1(Npts);
  for(int i=0; i<Npts; i++)
    cvec1.put(Pts[i].Sph2Cyl(rad),i);
  return cvec1;
  }  


coord_vec coord_vec::Cyl2Sph(int rad) const

        // Input                cvec    : Coordinate vector
        //                                assumed cylindrical { R, theta, z }
        //                      rad     : Radians vs. Deg angle input
        // Output               pt1     : Equivalent vector in Spherical
        //                                coordinates { R, theta, phi }

  {  
  coord_vec cvec1(Npts);
  for(int i=0; i<Npts; i++)
    cvec1.put(Pts[i].Cyl2Sph(rad),i);
  return cvec1;
  }

#endif						// coord_vec.cc
