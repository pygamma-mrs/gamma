/* FrameMakerSph.h **********************************************-*-c++-*-
**									**
**	                          G A M M A 				**
**								 	**
**	FrameMaker 3D Spherical Plots                     Interface	**
**						 			**
**	Copyright (c) 1991, 1992, 2001				 	**
**	Dr. Scott A Smith						**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fuer physikalische Chemie				 	**
**	8092 Zurich / Switzerland				 	**
**								 	**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  The FrameMaker classes and Routines provide GAMMA users with the    **
**  ability to output various data sets to Adobe FrameMaker in MIF      **
**  format.  The functions herein handle generation of 3D spherical     **
**  plots. That is, a set of 3D coordinates (coord_vec) may be plotted  **
**  with respect to a choosen set of Cartesian axes (x,y,z). The plot   **
**  will be restrained to reside within a sphere of radius R, scaling   **
**  will be done to all points to insure this is true.                  **
**                                                                      **
**  Users have the option to alter the following:                       **
**                                                                      **
**  1.) The radius of the sphere containing the plots.                  **
**  2.) How the points are plotted.                                     **
**      a.) As a continuous PolyLIne   (e.g. a trajectory)              **
**      b.) As a discrete points       (e.g. probablity density)        **
**      c.) As vectors from the origin (e.g. magnetization vectors)     **
**  3.) Include a drawing of the three axes                             **
**  4.) Include a drawing of the three plane/sphere interfaces.         **
**  5.) Specify the orientation of the coordiante axes.                 **
**                                                                      **
**  Since the output is editable in FrameMaker, much of the cosmentic   **
**  work to make these plots presentable can be done after their        **
**  production.                                                         **
**                                                                      **
*************************************************************************/

#ifndef   GFrameMakerSph_h_		// Is this file already included?
#  define GFrameMakerSph_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <GamIO/FrameMakerP.h>		// Know about FM parameters
#include <Matrix/row_vector.h>		// Know about row vectors
#include <Level1/coord_vec.h>		// Know about coordinate vectors 
#include <Level2/EAngles.h>		// Know about Euler angles
#include <string>			// Know about libstdc++ string
#include <vector>			// Know about libstdc++ STL vectors

extern const coord DefOrient;		// Default 3D Plot Orientation

class FMSph
  {
         double    R;			// Sphere radius
         EAngles   EAs;			// Orientation Euler Angles
         int       NPlane;		// No. Points To Draw Plane
         int       PType;		// Type of plot to make
         int       LType;		// Type of line to make
         double    Marg;		// Margin (on all 4 sides)
  static int       AxisID;		// Coordinate axis ID base number
  static int       PlaneID;		// Coordinate plane ID base number
  std::vector<int> AxisIDs;		// Coordinate axis ID number
  std::vector<int> PlaneIDs;		// Coordinate plane ID number

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                FrameMaker Spherical Plot Error Handling
// ____________________________________________________________________________

        // Input                sys     : Spin system (this)
        //                      ei      : Error index
        //                      nr      : Flag for linefeed (0=linefeed)
        //                      pn      : string in message

/*
         void FMSerror(int ei,                        int nr=0) const;
         void FMSerror(int ei, const std::string& pn, int nr=0) const;
volatile void FMSfatal(int ei) const;
volatile void FMSfatal(int ei, const std::string& pn) const;
*/

// ____________________________________________________________________________
// ii                FrameMaker Spherical Plot Start & End
// ____________________________________________________________________________

/* These are private because they assume that 1.) the output stream is
   viable and 2.) the radius has been set to something reasonable.           */

void start(std::ostream& out) const;
void end(std::ostream&   out) const;

// ____________________________________________________________________________
// iii                FrameMaker Spherical Coordinate Axes
// ____________________________________________________________________________

/* These are private because they assume that 1.) the output stream is
   viable and 2.) both the radius and orientation has been set to something
   reasonable.                                                               */

void axes(std::ostream& out);

// ____________________________________________________________________________
// iii                FrameMaker Spherical Coordinate Planes
// ____________________________________________________________________________

/* These are private because they assume that 1.) the output stream is
   viable and 2.) both the radius and orientation has been set to something
   reasonable.                                                               */

void plane(std::ostream&  out, char plane, int ID) const;
void planes(std::ostream& out);

// ____________________________________________________________________________
// iv                  FrameMaker Spherical Points
// ____________________________________________________________________________

row_vector project(const coord_vec& data) const;
void       draw(std::ostream& out, const row_vector& vx, int FMID=305) const;

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________
// A            FrameMaker Spherical Plot Constructors, Detructor
// ____________________________________________________________________________


MSVCDLC FMSph();
MSVCDLC FMSph(double R, const EAngles& EA);
MSVCDLC FMSph(const FMSph& FMS);

// ____________________________________________________________________________
// B            FrameMaker Spherical Plot Access Functions
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                      Sphere Radius & Plot Margin Access
// ----------------------------------------------------------------------------

MSVCDLL double Radius() const;
MSVCDLL void   Radius(double Rad);
MSVCDLL double Margin() const;
MSVCDLL void   Margin(double Mar);

// ----------------------------------------------------------------------------
//                      Sphere Orientation Access
// ----------------------------------------------------------------------------

/* Externally all single angles are expressed in degrees. Internally they are
   stored as Euler angles which are in radians. The conversion is performed
   here when needed.                                                         */

MSVCDLL const EAngles& Orient() const;
MSVCDLL void  Orient(const EAngles& EA);
MSVCDLL void  Orient(double A, double B, double G);

MSVCDLL double Alpha() const;
MSVCDLL double Beta()  const;
MSVCDLL double Gamma() const;

MSVCDLL void   Alpha(double A);
MSVCDLL void   Beta(double  B);
MSVCDLL void   Gamma(double G);

// ----------------------------------------------------------------------------
//                         Spherical Plot Type Access
// ----------------------------------------------------------------------------

MSVCDLL int  PlotType() const;
MSVCDLL void PlotType(int PT);


// ----------------------------------------------------------------------------
//                      Spherical Plot Line Type Access
// ----------------------------------------------------------------------------

MSVCDLL int  LineType() const;
MSVCDLL void LineType(int LT);

// ----------------------------------------------------------------------------
//                      Spherical Plane Points Access
// ----------------------------------------------------------------------------

MSVCDLL int  PlanePts() const;
MSVCDLL void PlanePts(int PP);

// ____________________________________________________________________________
// C            FrameMaker Spherical Plot Plotting Functions
// ____________________________________________________________________________

MSVCDLL void plot(const std::string& filename);
MSVCDLL void plot(const std::string& filename, const coord_vec& cvec);
MSVCDLL void plot(const std::string& filename, const std::vector<coord_vec>& cvs);


// ____________________________________________________________________________
// D             FrameMaker Spherical Plot ASCII Ouptut Functions
// ____________________________________________________________________________

MSVCDLL std::ostream& print(std::ostream& ostr) const;

};						// End Of Class FrameMakerSph

// ____________________________________________________________________________
// E             NonMember FrameMaker Spherical Plot Functions
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//              Simple NonMember FrameMaker Spherical Plot Functions
// ----------------------------------------------------------------------------

/* These are intended to be very simple function calls that need few very
   arguements. Thus GAMMA users can easily generate 3D plots without any
   fretting over added details. 

	   Input	filename  : Output filename
			data      : Set of {x,y,z} points
                        type      : Plot type
	   		alpha     : Euler angle (degrees)
	   		beta      : Euler angle (degrees)
	   		gamma     : Euler angle (degrees)
                    or  EA        : Euler angles (radians)
	  		radius    : Plot (x & y) dimension in cm
			points    : Number of points in plane plots
	   Return		  : Void. xy-plot is output to file          */

MSVCDLL void FrameMakerSphere(const std::string& name,
      int type=2, const coord& EAs=DefOrient, double radius=5, int points=100);

MSVCDLL void FrameMakerSphere(const std::string& name, const coord_vec& data,
      int type=2, const coord& EAs=DefOrient, double radius=5, int points=100);

MSVCDLL void FrameMakerSphere(const std::string& name, const std::vector<coord_vec>& cvs,
      int type=2, const coord& EAs=DefOrient, double radius=5, int points=100);

// ----------------------------------------------------------------------------
//              Other NonMember FrameMaker Spherical Plot Functions
// ----------------------------------------------------------------------------

MSVCDLL void FM_sphere(const std::string& filename, int type=2, double alpha=0,
           double beta=15, double gamma=-15, double radius=5, int points=100);

MSVCDLL void FM_sphere(const std::string& filename, const coord_vec& data, int type=0,
	 double alpha=0, double beta=15, double gamma=-15,
                                             double radius=5, int points=100);

MSVCDLL void FM_sphere(const std::string& filename, coord_vec& data1,
     coord_vec& data2, int type=0, double alpha=0, double beta=15,
		           double gamma=-15, double radius=5, int points=100);

#endif 							// FrameMakerSph.h
