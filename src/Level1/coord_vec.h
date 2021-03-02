/* coord_vec.h **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**	Coordinate Vector                          Interface		**
**								    	**
**	Copyright (c) 1990, 1991, 1992				    	**
**	Scott Smith				  		  	**
**	Eidgenoessische Technische Hochschule	    			**
**	Labor fuer physikalische Chemie		    			**
**	8092 Zuerich / Switzerland		    			**
**						   	 		**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  Class coordinate vector provides a facile means to manipulate a     **
**  general list of coordinates.  Each variable of type coord_vec	**
**  consists of a vector of coordinates in three dimensional space      **
**  and the number of points the vector contains.  The entire vector    **
**  can be manipulated as a whole(rotated, tranlated, etc.) in the same **
**  manner that a single  coordinate point (class coord) may be         **
**  manipulated.                                                        **
**                                                                      **
*************************************************************************/

///Chapter Class Coordinate Vector
///Section Overview
///Body	   None.
///Section Available Coordinate Vector Functions

#ifndef   Coord_vec_h_			// Is file already included?
#  define Coord_vec_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Level1/coord.h>		// Provides knowledge of coordinates
#include <Matrix/row_vector.h>		// Provides knowledge of row vectors
#include <Matrix/matrix.h>		// Provides knowledge of matrices
#include <Basics/ParamSet.h>		// Provides knowledge of parameter sets
#include <string>			// Know about libstdc++ strings

class coord_vec
{
    
// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------

///Center Coordinate Vector Algebraic

  coord *Pts;			// Pointer to Array of Coordinates
  int    Npts;			// Number of Coordinates
  
private:


 
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

void CVerror(int eidx, int noret=0) const;
void CVerror(int eidx, const std::string& pname, int noret=0) const;
volatile void CVfatality(int eidx) const;

// ____________________________________________________________________________
// ii                  COORDINATE VECTOR SETUP FUNCTIONS
// ____________________________________________________________________________
 
/* These functions set up specific aspects of a coordinate vector.  Since
   they make assumptions about the order in which the vector is set up, the
   functions MUST be private, their misuse may make an inconsistent vector.  */  

        // Input                cvec	: A Coordinate vector(this)
        //                      pset	: A parameter set
        //                      warn	: Warning level
        //                      SetNPoints      0 = no warning
        //                                      1 = warning, non-fatal
        //                                     >1 = fatal
        //                      SetCoords      -1 = warning, def 0 (DEF)
        //                                      0 = no warning, def 0
        //                                      1 = warning, def big
        //                                     >2 = fatal
        // Output SetNPoints	none	: Coordinate vector number of
        //                                points is set from parameters
        //                                in pset
        // Output SetCoords	none	: Coordinate vector coordinates
        //                                set from parameters in pset

int SetNPoints(const ParameterSet& pset, int warn=2);
int SetCoords(const  ParameterSet& pset, int warn=-1);

// ____________________________________________________________________________
// iii                 COORDINATE VECTOR CHECKING FUNCTIONS
// ____________________________________________________________________________
 
void check(int index) const;

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A             CLASS COORDINATE VECTOR CONSTRUCTORS/DESTRUCTOR
// ____________________________________________________________________________

///Center Constructors, Destructors, and Asssignment
///F_list coord_vec	     - Constructor
///F_list =		      - Assignment

/* These constructors set up a new coordinate vector.  There are several ways
   to make a coordinate vector as listed below:

          Input Arguments                       Resulting Coordinate
   --------------------------        -----------------------------------------
                -                    Empty coordinate vector (no points)
             int pts		     Coordinate vector with pts points
          coord_vec cvec             Coordinate vector identical to cvec
     ParameterSet, idx, warn         Coord. vector w/ prefix i read from pset
 
   Setting a coordinate vector from a parameter set may fail if the parameters 
   that define coordinates don't exist. The flag warn dictates what happens
   upon a failure: 0=no warnings, 1=non-fatal warnings, 2=fatal warnings
   Individual coordinate (see class coord) parameters assume a parameter of
   type 3 with value stored as a string "( #, #, # )". All I/O between GAMMA
   coordinates and single parameters must follow this format!              */

MSVCDLC         coord_vec( );
MSVCDLC         coord_vec(int pts);
MSVCDLC         coord_vec (const coord_vec& cvec1);
MSVCDLC         coord_vec(const ParameterSet& pset, int idx=-1, int warn=2);
MSVCDLC         coord_vec(const row_vector& X, const row_vector& Y, const row_vector& Z);
MSVCDLC virtual ~coord_vec ();
MSVCDLL coord_vec& operator= (const coord_vec& cvec1);

// ____________________________________________________________________________
// B                  COORDINATE VECTOR ROTATION FUNCTIONS
// ____________________________________________________________________________

///Center Rotation Functions

// ----------------------------------------------------------------------------
//            Typical Rotations Active On Cartesian Coordinates
// ----------------------------------------------------------------------------

///F_list rotate	     - Rotate a coordinate vector
///F_list rotate_ip	     - Rotate a coordinate vector in place
 
        // Input                cvec : Coordinate vector(this)
        //                      theta: Spherical angle (degrees down from +z)
        //                      phi  : Spherical angle (degrees over from +x)
        //                      rad  : Flag for angle in degrees or radians
        // Return               rotcv: The vector rotated about the x/y-axis
        //                             or the z-axis by the input angle
        // Note                      : The reference axes are static!  Thus
        //                             multiple rotations should be cumulative
 
MSVCDLL coord_vec xrotate(double theta, int rad=0) const;
MSVCDLL coord_vec yrotate(double theta, int rad=0) const;
MSVCDLL coord_vec zrotate(double phi,   int rad=0) const;

        // Input                cvec : Coordinate vector(this)
        //                      Rmx  : 3x3 Rotation matrix
	// 			alpha: Euler angle (radians)
	// 			beta : Euler angle (radians)
	// 			gamma: Euler angle (radians)
	// 			EA   : Coordinate containing Euler angles
	// Return		rcvec: Coordinate vector rotated
	//			       either by input Euler angles
	//			       or by input rotation matrix Rmx
	// Note			     : Angles input in radians

MSVCDLL coord_vec rotate(const matrix& Rmx) const;
MSVCDLL coord_vec rotate(double alpha, double beta, double gamma) const;
MSVCDLL coord_vec rotate(const coord& EA) const;

	// Input		cvec : Coordinate vector(this)
	// 			alpha: Euler angle (radians)
	// 			beta : Euler angle (radians)
	// 			gamma: Euler angle (radians)
	// 			EA   : Coordinate containing Euler angles
	// Return		cvec : Coordinate vector rotated
	//			       by input Euler angles
	//			       or by input Euler angles
	// Note			     : Angles input in radians

MSVCDLL void rotate_ip(double alpha, double beta, double gamma);
MSVCDLL void rotate_ip(const coord &EA);

// ____________________________________________________________________________
// C                  COORDINATE VECTOR TRANSLATION FUNCTIONS
// ____________________________________________________________________________

///Center Translation Functions
///F_list trans_x		- Translate along the x-axis
///F_list trans_y		- Translate along the y-axis
///F_list trans_z		- Translate along the z-axis
///F_list trans_x_ip		- Translate along the x-axis in place
///F_list trans_y_ip		- Translate along the y-axis in place
///F_list trans_z_ip		- Translate along the z-axis in place
///F_list tranlate		- Translate a coordinate vector
///F_list tranlate_ip		- Translate a coordinate vector in place

	// Input		cvec    : Coordinate vector (this)
	// 			delx    : Change in x coordinate
	// 			dely    : Change in y coordinate
	// 			delz    : Change in z coordinate
	// 			delpt   : Change in x,y,& z coordinates
	// Return (non-ip)	tcvec   : A coordinate vector which is cvec
	//				  translated by delx,dely, & delz
	//				  or by delpt
	// Return (ip)		void    : Coordinate vector cvec is
	//				  translated by del or by
	//				  delx,dely, & delz

MSVCDLL coord_vec translate(double    delx, double dely=0, double delz=0) const;
MSVCDLL void      translate_ip(double delx, double dely=0, double delz=0);
MSVCDLL coord_vec translate(const     coord& delpt) const;
MSVCDLL void      translate_ip(const  coord& delpt);

	// Input		cvec    : Coordinate vector (this)
	// 			delx    : Change in x coordinate
	// 			dely    : Change in y coordinate
	// 			delz    : Change in z coordinate
	// Return (non ip)	tcvec   : A coordinate vector which is cvec
	//				  translated along u axis by delu
	//				  where u = {x,y,z}
	// Return (ip)		void	: Coordinate vector is
	//				  translated along u axis by delu
	//				  where u = {x,y,z} directly

MSVCDLL coord_vec trans_x(double delx) const;
MSVCDLL coord_vec trans_y(double dely) const;
MSVCDLL coord_vec trans_z(double delz) const;

MSVCDLL void trans_x_ip(double delx);
MSVCDLL void trans_y_ip(double dely);
MSVCDLL void trans_z_ip(double delz);


// ____________________________________________________________________________
// D                  COORDINATE VECTOR PROJECTION FUNCTIONS
// ____________________________________________________________________________

///Center Projection Functions

MSVCDLL row_vector project(int projx, int projy) const;

	// Input		cvec    : Coordinate vector (this)
	// 			projx   : Axis flag to project onto x
	// 			projy   : Axis to project onto y
	//				  1 = x axis; -1 = -x axis
	//				  2 = y axis; -2 = -y axis
	//				  3 = x axis; -3 = -z axis
	// Return		vx      : Row vector having cvec (3D)
	//		  		  projected into it (2D) via axes
	//				  mapping specified by projx & projy
	///F_list project		- Coordinate vector projection


// ____________________________________________________________________________
// E                    COORDINATE VECTOR OPERATORS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                      COORDINATE VECTOR WITH SCALAR
// ----------------------------------------------------------------------------

///Center Coordinate Vector with Scalar Functions
///F_list *		      - Multiplicaton by scalar
///F_list *=		      - Unary multiplicaton by scalar
///F_list /		      - Division by scalar
///F_list /=		      - Unary division by scalar

/*  Function/Operator                         Result
    -----------------   -------------------------------------------------------
         cv*r           Multiplies all coordinates in vector by r
         r*cv           Multiplies all coordinates in vector by r
        cv*=r           Multiplies all coordinates in vector by r
         cv/r           Multiplies all coordinates in vector by 1/r
        cv/=r           Multiplies all coordinates in vector by 1/r          */
      
MSVCDLL friend    coord_vec  operator *  (double r, const coord_vec& cvec1);

MSVCDLL coord_vec  operator *  (double r) const;
MSVCDLL coord_vec& operator *= (double r);
MSVCDLL coord_vec  operator /  (double r) const;
MSVCDLL coord_vec& operator /= (double r);

// ----------------------------------------------------------------------------
//                      COORDINATE VECTOR WITH COORDINATE
// ----------------------------------------------------------------------------

/*  Function/Operator                         Result
    -----------------   -------------------------------------------------------
         cv+pt          Adds      point p to all coordinates in vector
         cv+=pt         Adds      point p to all coordinates in vector
         cv-pt          Subtracts point p from all coordinates in vector
         cv-=pt         Subtracts point p from all coordinates in vector     */

MSVCDLL coord_vec operator +  (const coord& pt) const;
MSVCDLL void      operator += (const coord& pt);
MSVCDLL coord_vec operator -  (const coord& pt) const;
MSVCDLL void      operator -= (const coord& pt);

// ----------------------------------------------------------------------------
//                  COORDINATE VECTOR WITH COORDINATE VECTOR
// ----------------------------------------------------------------------------

MSVCDLL coord_vec  operator +  (const coord_vec& cv) const;
MSVCDLL coord_vec& operator += (const coord_vec& cv);
MSVCDLL coord_vec  operator -  (const coord_vec& cv) const;
MSVCDLL coord_vec& operator -= (const coord_vec& cv);

// ____________________________________________________________________________
// F                COORDINATE VECTOR GENERAL FUNCTIONS
// ____________________________________________________________________________

///Center Coordinate Vector General Functions
///F_list size			- Number of coordinates in vector
///F_list max_x			- Maximum x value in coordinates vector
///F_list max_y			- Maximum y value in coordinates vector
///F_list max_z			- Maximum z value in coordinates vector
///F_list distances		- Matrix of distances between points
///F_list thetas		- Matrix of thetas between points
///F_list phis			- Matrix of phis between points

	// Input		cvec  : A coordinate vector(this)
	// Return		Npts  : Number of points in cvec

MSVCDLL int size() const;

	// Input		cvec  : A coordinate vector(this)
	// Return		maxx  : Maximum x value
	// Return		mayy  : Maximum y value
	// Return		mazz  : Maximum z value

MSVCDLL double max_x() const;
MSVCDLL double max_y() const;
MSVCDLL double max_z() const;


MSVCDLL coord maxima() const;

	// Input		cvec  : A coordinate vector(this)
	// Return		pt    : Coordinate point whose x, y, & z
	//				contain maximum values over the vector
	///F_list maxima	      - Maximum x,y, & z values in vector



MSVCDLL void maxima(double &x, double &y, double &z) const;

	// Input		cvec  : A coordinate vector(this)
	// 			x     : x coordinate
	// 			y     : y coordinate
	// 			z     : z coordinate
	// Return		none  : x, y, & z filled with the
	//			        maximum values over the vector


MSVCDLL double max_R() const;

	// Input		cvec  : A coordinate vector(this)
	// Return		maxR  : Maximum radius value
	///F_list max_R	      	      - Maximum radius of coodinate in vector


MSVCDLL void max_R(int &maxi, double &maxR) const;

	// Input		cvec  : A coordinate vector(this)
	// Return		void  : Maximum radius value
	//				is copied into maxR and
	//				the index of this coordinate
	//				is copied into maxi


MSVCDLL coord_vec vectors() const;

	// Input		cvec  : A coordinate vector(this)
	// Return		cvect : Coordinates of the vectors
	//				connecting pairs of points
	///F_list vectors      	      - Vectors connecting vector points


MSVCDLL coord_vec vectors_f() const;

	// Input		cvec  : A coordinate vector(this)
	// Return		cvect : Coordinates of the vectors
	//				connecting pairs of points




MSVCDLL double distance(int pt1, int pt2, int Angs=0) const;

	// Input		cvec  : A coordinate vector(this)
	//			pt1   : First coordinate point
	//			pt2   : Second coordinate point
	//			Angs  : Flag for output in Angstroms
	// Return		dist  : Distance between points pt1 & pt2


	// Input		cvec  : A coordinate vector(this)
	//			Angs  : Flag for output in Angstroms
	//			deg   : Flag for output in degrees
	// Return		the_mx: Matrix of theta angles between
	//				the points of cvec
	// Return		phi_mx: Matrix of phi angles between
	//				the points of cvec
	// Return		dmx   : A distance matrix between
	//				the points of cvec

MSVCDLL matrix distances(int Angs=0) const;
MSVCDLL matrix thetas(int    deg=0)  const;
MSVCDLL matrix phis(int      deg=0)  const;

// ____________________________________________________________________________
// G                      INDIVIDUAL COORDINATE ACCESS
// ____________________________________________________________________________

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

///Center Coordinate Access Functions
///F_list ()      	      - Coordinate access
///F_list put      	      - Set individual coordinate
///F_list get      	      - Get individual coordinate
///F_list x      	      - X-component of individual coordinate
///F_list y      	      - y-component of individual coordinate
///F_list z      	      - z-component of individual coordinate

MSVCDLL coord                operator() (int index)                    const;
MSVCDLL void      put(const coord &pt1, int index);
MSVCDLL void      put(double x, double y, double z, int index);
MSVCDLL coord     get(int index)                            const;
MSVCDLL double    x(int index)                              const;
MSVCDLL double    y(int index)                              const;
MSVCDLL double    z(int index)                              const;
MSVCDLL coord_vec get_block(int index, int npts)            const;
MSVCDLL void      put_block(int index, const coord_vec& cv) const;

// ____________________________________________________________________________
// H                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

///Center Parameter Set Functions


// ----------------------------------------------------------------------------
//      Functions To Make A Parameter Set From A Coordinate Vector
// ----------------------------------------------------------------------------


	// Input		cvec  : Coordinate vector
	//  			pset  : Parameter set
	// Output		pset  : Parameter set with
	//			        only coordinates
	///F_list =      	      - Conversion to parameter set

MSVCDLL             operator ParameterSet( ) const;
MSVCDLL friend void operator+= (ParameterSet& pset, const coord_vec& cvec);

	// Input		cvec  : Coordinate vector
	//  			pset  : Parameter set
	// Output		pset  : Parameter set with
	//			        only spin system parameters
	///F_list +=      	      - Unary addition to parameter set


MSVCDLL void PSetAdd(ParameterSet& pset, int idx=-1) const;

        // Input                cvec    : Coordinate vector
        //                      pset    : Parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        // Output               void    : Coordinate vector parameters are
        //                                are added ot the parameter set
        //                                with interaction index idx
        // Note                         : The parameter names & types
        //                                output here MUST match those used
        //                                in setting the coordinate vector
        //                                from parameters sets


MSVCDLL int operator= (const ParameterSet& pset);

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
	///F_list =      	         - Assignment from parameter set


// ----------------------------------------------------------------------------
//    Functions To Output Coordinate Vector To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

        // Input                cvec     : Coordinate vector
        //                      filename : Output file name
        //                      ofstr	 : Output file stream
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none     : Coordinate vector is written as a
        //                                 parameter set to file filename or to
	//				   output filestream
        // Note                         : This depends on function PSetAdd!


MSVCDLL int write(const std::string &filename, int idx=-1, int warn=2) const;
MSVCDLL int write(std::ofstream& ofstr,        int idx=-1, int warn=2) const;
 
// ____________________________________________________________________________
// I                COORDINATE VECTOR AUXILIARY FUNCTIONS
// ____________________________________________________________________________

MSVCDLL bool is_zero( ) const;
 
        // Input                cvec  : Coordinate vector (this)
        // Output               TF    : False if any coordinate in
        //                              cvec is non-zero

// ____________________________________________________________________________
// J                 CLASS COORDINATE VECTOR INPUT FUNCTIONS
// ____________________________________________________________________________

///Center Coordinate Vector I/O Functions
 
MSVCDLL virtual std::string ask_read(int argc, char* argv[], int argn,
                                                       int idx=-1, int warn=2);
 
        // Input                cvec    : Dipolar interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
	//			idx	: Vector index
        //                      warn    : Warning output level
        // Output               string  : Parameter argn of array argc is used
        //                                to supply a filename from which the
        //                                coordinate vector coords. are read
        //                                If argument argn is not in argv, the
        //                                user is asked to supply a filename
        //                                The set filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing known coord_vec params.
        // Note                         : The coordinate vector is modifed



	// Input		cvec    : Coordinate vector
	// 			filename: Input filename
        //                      pset    : A parameter set
	//			idx	: Vector index
        //                      warn    : Warning output level
        // Output               T/F	: Coordinate vector filled with
        //                                coordinates in file or in pset.
	//				  Return is True if read properly.
	///F_list read      	        - Input from disk file

MSVCDLL virtual int read(const std::string& filename, int idx=-1, int warn=2);
MSVCDLL virtual int read(const ParameterSet& pset,    int idx=-1, int warn=2);

// ____________________________________________________________________________
// K               CLASS COORDINATE VECTOR OUTPUT FUNCTIONS
// ____________________________________________________________________________


	// Input		cvec  : Coordinate vector (this)
        //                      out   : Output stream
	//			units : Flag for output type
        // Output               none  : Coordinate vector is sent
	//				to the output stream
	///F_list print               - Formatted output to ostream

MSVCDLL std::ostream& printf(std::ostream& out, int units = 0) const;
MSVCDLL std::ostream& printCylindrical(std::ostream& ostr, double sf=1) const;
 
        // Input                cvec  : Coordinate vector (this)
        //                              assumed Cartesian {x,y,z}
        //                      ostr  : Output stream
 
 
MSVCDLL std::ostream& printSpherical(std::ostream& ostr, double sf=1) const;
 
        // Input                cvec  : Coordinate vector (this)
        //                              assumed Cartesian {x,y,z}
        //                      ostr  : Output stream

 
MSVCDLL std::ostream& printCartesian(std::ostream& ostr, double sf=1) const;
 
        // Input                cvec  : Coordinate vector (this)
        //                              assumed Cartesian {x,y,z} 
        //                      ostr  : Output stream

	// Input		cvec  : Coordinate vector (this)
        //                      ostr  : Output stream
	// Output 		ostr  : Returns the modified output stream
        // Output               none  : constaining cvec values
	///F_list print               - Output to ostream
	///F_list <<                  - Output to ostream

MSVCDLL        std::ostream& print(std::ostream& ostr, int ncols=2, int N=-1) const;
MSVCDLL friend std::ostream& operator << (std::ostream& ostr, const coord_vec &cvec);

// ____________________________________________________________________________
// L             CLASS COORDINATE VECTOR CONVERSION FUNCTIONS
// ____________________________________________________________________________
 
//-----------------------------------------------------------------------------
//               Cartesian <---------> Spherical Conversion
//-----------------------------------------------------------------------------

        // Input                cvec    : Coordinate vector
        //                                assumed to be Cartesian {x, y, z }
        //                      rad     : Radians vs. Deg angle output/input
	// Output               cvec1   : Equivalent vector in Spherical
        //                                coordinates { R, theta, phi }
 
MSVCDLL coord_vec Cart2Sph(int rad=1) const;
MSVCDLL coord_vec Sph2Cart(int rad=1) const;

//-----------------------------------------------------------------------------
//               Cartesian <---------> Cylindrical Conversion
//-----------------------------------------------------------------------------

        // Input                cvec    : Coordinate vector
        //                                assumed to be Cartesian {x, y, z }
        //                      rad     : Radians vs. Deg angle output/input
        // Output               pt1     : Equivalent vector in Cylindrical
        //                                coordinates { R, theta, z }
 
MSVCDLL coord_vec Cart2Cyl(int rad=1) const;
MSVCDLL coord_vec Cyl2Cart(int rad=1) const;

//-----------------------------------------------------------------------------
//               Spherical <---------> Cylindrical Conversion
//-----------------------------------------------------------------------------

        // Input                cvec    : Coordinate vector
        //                                assumed spherical { R, theta,phi }
        //                      rad     : Radians vs. Deg angle output
        // Output               pt1     : Equivalent vector in Cylindrical
        //                                coordinates { R, theta, z }

MSVCDLL coord_vec Sph2Cyl(int rad=1) const;
MSVCDLL coord_vec Cyl2Sph(int rad=1) const;


};

#endif 						// coord_vec.h
