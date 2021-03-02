/* FrameMakerSph.cc *********************************************-*-c++-*-
**									**
**	                        G A M M A 				**
**								 	**
**	FrameMaker 3D Spherical Plots	             Implementation   	**
**						 			**
**	Copyright (c) 1991, 1992, 2001		 			**
**	Scott A. Smith							**
**	Eidgenoessische Technische Hochschule			 	**
**	Labor fuer physikalische Chemie				 	**
**	8092 Zurich / Switzerland				 	**
**								 	**
**      $Header: $
**	                         				 	**
*************************************************************************/

/*************************************************************************
**								 	**
**  Description							 	**
**						 			**
**  The FrameMaker classes and Routines provide GAMMA users with the	**
**  ability to output various data sets to Adobe FrameMaker in MIF	**
**  format.  The functions herein handle generation of 3D spherical	**
**  plots. That is, a set of 3D coordinates (coord_vec) may be plotted	**
**  with respect to a choosen set of Cartesian axes (x,y,z). The plot	**
**  will be restrained to reside within a sphere of radius R, scaling	**
**  will be done to all points to insure this is true.  		**
**								 	**
**  Users have the option to alter the following:			**
**								 	**
**  1.) The radius of the sphere containing the plots.			**
**  2.) How the points are plotted.					**
**      a.) As a continuous PolyLIne   (e.g. a trajectory)		**
**      b.) As a discrete points       (e.g. probablity density)	**
**      c.) As vectors from the origin (e.g. magnetization vectors)	**
**  3.) Include a drawing of the three axes				**
**  4.) Include a drawing of the three plane/sphere interfaces.		**
**  5.) Specify the orientation of the coordiante axes.			**
**	 							 	**
**  Since the output is editable in FrameMaker, much of the cosmentic	**
**  work to make these plots presentable can be done after their	**
**  production.							 	**
**	 							 	**
**  If there is no orientation, the coordiantes are projected onto the	**
**  plane of view as a right-hand coordinate system with the z-axis	**
**  vertical, the y-axis pointing to the viewer, and the x-axis hor-	**
**  izontal to the left.  This perspective will change if an Euler	**
**  rotation is specified { alpha, beta, gamma }. The Euler rotations	**
**  are implemented as a rotation about occur as a rotation about y by	**
**	 							 	**
**            Z
**            ^
**            |
**            |
**      X <---+--- 
**            |
**            |
**            |
**	 							 	**
**         (0,0,0}
**	 							 	**
*************************************************************************/

#ifndef   FrameMakerSph_cc_			// Is file already included?
#  define FrameMakerSph_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#  endif

#include <GamIO/FrameMakerSph.h>		// Include the header
#include <GamIO/FrameMakerM.h>			// Include FM MIF functions
#include <GamIO/FrameMakerP.h>			// Include FM parameters
#include <Basics/Gconstants.h>			// Include PI
#include <Level2/EAngles.h>			// Inlcude Euler angles
#include <stdlib.h>
#include <string>				// Inlcude stdlibc++ strings
#include <iostream>				// Include iostreams (cout)
#include <fstream>				// Include filestreams
#include <vector>				// Include libstdc++ STL vectors

using std::string;				// Using libstdc++ strings
using std::ostream;				// Using libstdc++ output streams
using std::ofstream;				// Using libstdc++ output file streams
using std::vector;				// Using libstdc++ vectors
using std::cout;				// Using libstdc++ standard output

const coord DefOrient(0.0,15.0,5.0);		// Default 3D Plot Orientation
      int FMSph::AxisID  = 250;			// Initialize axis  ID number
      int FMSph::PlaneID = 260;			// Initialize plane ID number
 

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                FrameMaker Spherical Plot Error Handling
// ____________________________________________________________________________

// ____________________________________________________________________________
// ii                FrameMaker Spherical Plot Start & End
// ____________________________________________________________________________

/* These are private because they assume that 1.) the output stream is
   viable and 2.) the radius has been set to something reasonable.           */

void FMSph::start(ostream& out) const
  {
  FM_Begin(out);				// FrameMaker begin comment
  FM_AFrames_Begin(out);			// Anchored Frame begin
  FM_AFrame_Set(out, 2*R+2,2*R+2);		// Header of Frame (1cm border)
  }

void FMSph::end(ostream& out) const
  {
  FM_AFrame_End(out);				// End Anchored Frame
  FM_AFrames_End(out);				// End Anchored Frames
  FM_ParaText_End(out);				// End Frame TextFlow
  FM_End(out);					// End Frame Document
  }

// ____________________________________________________________________________
// iii                FrameMaker Spherical Coordinate Axes
// ____________________________________________________________________________

/* These are private because they assume that 1.) the output stream is
   viable and 2.) both the radius and orientation has been set to something 
   reasonable. The axes are drawn in three segments, the last being any arrow-
   head. The section that spans the sphere is cut in half, the back half drawn
   as a dashed line and the front half drawn as a solid line.                */

void FMSph::axes(ostream &out)
  {
//                 Set Up Coordinate Axes & Rotate Them
//          (Positive End Sticks Out Beyond Sphere Of Radius R)

  coord x1(-R    ,   0.0,   0.0);		// Unrotated -x
  coord x2( R    ,   0.0,   0.0);		// Unrotated  x
  coord x3( R+0.9,   0.0,   0.0);		// Unrotated  x end
  coord y1(   0.0,-R    ,   0.0);		// Unrotated -y
  coord y2(   0.0, R    ,   0.0);		// Unrotated  y
  coord y3(   0.0, R+0.9,   0.0);		// Unrotated  y end
  coord z1(   0.0,   0.0,-R    );		// Unrotated -z
  coord z2(   0.0,   0.0, R    );		// Unrotated  z
  coord z3(   0.0,   0.0, R+0.9);		// Unrotated  z end

  double alpha = EAs.alpha();
  double beta  = EAs.beta();
  double gamma = EAs.gamma();
  x1=x1.rotate(alpha,beta,gamma);		// Rotated -x
  x2=x2.rotate(alpha,beta,gamma);		// Rotated  x
  x3=x3.rotate(alpha,beta,gamma);		// Rotated  x end
  y1=y1.rotate(alpha,beta,gamma);		// Rotated -y
  y2=y2.rotate(alpha,beta,gamma);		// Rotated  y
  y3=y3.rotate(alpha,beta,gamma);		// Rotated  y end
  z1=z1.rotate(alpha,beta,gamma);		// Rotated -z
  z2=z2.rotate(alpha,beta,gamma);		// Rotated  z
  z3=z3.rotate(alpha,beta,gamma);		// Rotated  z end

//                 Translate Coodinate Axes (Centering)

  x1=x1.translate(R+Marg,0,R+Marg);		// Centered -x
  x2=x2.translate(R+Marg,0,R+Marg);		// Centered  x
  x3=x3.translate(R+Marg,0,R+Marg);		// Centered  x end
  y1=y1.translate(R+Marg,0,R+Marg);		// Centered -y
  y2=y2.translate(R+Marg,0,R+Marg);		// Centered  y
  y3=y3.translate(R+Marg,0,R+Marg);		// Centered  y end
  z1=z1.translate(R+Marg,0,R+Marg);		// Centered -z
  z2=z2.translate(R+Marg,0,R+Marg);		// Centered  z
  z3=z3.translate(R+Marg,0,R+Marg);		// Centered  z end

//                     Set Up For Output In MIF

  double ALW = 1.0;				// Axis line width
  coord zero(R+Marg,0,R+Marg);			// Origin
  AxisIDs.clear();				// No Axis ID numbers

//                  Draw In The Three Cartesian Axes
//           Now Oriented And Projected, 3D {x,z} -> 2D {-x,-y}

  double D=2.0*(R+Marg);
  double delx, delz;
  delx = x1.x() - x2.x();
  delz = x1.z() - x2.z();
  int pen;
  if(sqrt(delx*delx + delz*delz) > R/1.e4)
    {
    (x1.y() >= 0)?pen=0:pen=10;
    FM_Line(out,AxisID,0,ALW,D-x1.x(),D-x1.z(),zero.x(), zero.z(),pen);// X-axis
    (x2.y() >= 0)?pen=0:pen=10;
    FM_Line(out,AxisID,0,ALW,zero.x(),zero.z(),D-x2.x(),D-x2.z(), pen);// X-axis
    FM_Line(out,AxisID,1,ALW,D-x2.x(),D-x2.z(),D-x3.x(),D-x3.z());	// X-axis end
    FM_TextLine(out,AxisID,D-x3.x(),D-x3.z(),'X',10,0);		// X-axis label
    AxisIDs.push_back(AxisID);						// Store axis ID
    AxisID++;								// Increment ID
    }

  delx = y1.x() - y2.x();
  delz = y1.z() - y2.z();
  if(sqrt(delx*delx + delz*delz) > R/1.e4)
    {
    (y1.y() >= 0)?pen=0:pen=10;
    FM_Line(out,AxisID,0,ALW,D-y1.x(),D-y1.z(),zero.x(),zero.z(),pen);	// Y-axis
    (y2.y() >= 0)?pen=0:pen=10;
    FM_Line(out,AxisID,0,ALW,zero.x(),zero.z(),D-y2.x(),D-y2.z(),pen);	// Y-axis
    FM_Line(out,AxisID,1,ALW,D-y2.x(),D-y2.z(),D-y3.x(),D-y3.z());	// Y-axis end
    FM_TextLine(out,AxisID,D-y3.x(),D-y3.z(),'Y',10);			// Y-axis label
    AxisIDs.push_back(AxisID);						// Store axis ID
    AxisID++;								// Increment ID
    }

  delx = z1.x() - z2.x();
  delz = z1.z() - z2.z();
  if(sqrt(delx*delx + delz*delz) > R/1.e4)
    {
    (z1.y() >= 0)?pen=0:pen=10;
    FM_Line(out,AxisID,0,ALW,D-z1.x(),D-z1.z(),zero.x(),zero.z(),pen);	// Z-axis
    (z2.y() >= 0)?pen=0:pen=10;
    FM_Line(out,AxisID,0,ALW,zero.x(),zero.z(),D-z2.x(),D-z2.z(),pen);	// Z-axis
    FM_Line(out,AxisID,1,ALW,D-z2.x(),D-z2.z(),D-z3.x(),D-z3.z());	// Z-axis end
    FM_TextLine(out,AxisID,D-z3.x(),D-z3.z(),'Z',10);			// Z-axis label
    AxisIDs.push_back(AxisID);						// Store axis ID
    AxisID++;								// Increment ID
    }
  }

// ____________________________________________________________________________
// iii                FrameMaker Spherical Coordinate Planes
// ____________________________________________________________________________

/* These are private because they assume that 1.) the output stream is
   viable and 2.) both the radius and orientation has been set to something 
   reasonable.                                                               */

void FMSph::plane(ostream &out, char plane, int ID) const
  {
  vector<coord> C(NPlane);			// Coordinates for points
  double xx,yy;					// For temp x,y ordinates
  double theta;					// For temp angle
  int k;					// For point indexing
  for(k=0; k<NPlane; k++)			// Generate circle 3D data npts
    {
    theta = k*PIx2/(NPlane-1);			// Angle to point on circle
    xx = R*cos(theta);				// X ordinate of this point
    yy = R*sin(theta);				// Y ordinate of this point
    switch(plane)				// 
      {
      case 'z':					// Circle in xy-plane
        C[k].x(xx);				// Map circle {x,y} -> {x,y,0}
        C[k].y(yy);
        C[k].z(0);
        break;
      case 'y':					// Circle in xz-plane
      default:					// Map circle {x,y} -> {x,0,y}
        C[k].x(xx);
        C[k].y(0);
        C[k].z(yy);
        break;
      case 'x':					// Circle in yz-plane
        C[k].x(0);				// Map circle {x,y} -> {0,x,y}
        C[k].y(xx);
        C[k].z(yy);
        break;
      }
    }
//             Having Generated Coordinates For Circle of Plane
//        We Need To Rotate To Proper Orientation & Translate To Center
//               (Assumes 1cm Border Around Sphere of Radius R)

  double alpha = EAs.alpha();
  double beta  = EAs.beta();
  double gamma = EAs.gamma();
  for(k=0; k<NPlane; k++)			// Loop points in circle
    {		
    C[k] = C[k].rotate(alpha,beta,gamma);	//   Rotate by Euler angles
    C[k] = C[k].translate(R+1,0,R+1);		//   Translate to center
    }

//        We Now Have Proper 3D Coordinates For Desired Circle (Plane)
//            Now It Must Be Projected Onto The 2D Plane Of View
//     This Is Done Using FM PolyLines. In Addition, Points In Front Of
//      POV Should Be A Solid Line While Those Behind Whould Be Dashed

  row_vector vx(NPlane);			// Vector for 2D points
  double D = 2.0*(R+Marg);
  xx = D-C[0].x();				// Get initial point of circle
  yy = D-C[0].z();
  vx.put(complex(xx,yy), 0);			// Store 1st Point
  int np = 1;
  int pen = 0;					// Set pen solid
  if(C[0].y() > 0) pen = 10;			// Line should be dashed
  for(k=1; k<NPlane; k++)
    {
    xx = D-C[k].x();
    yy = D-C[k].z();
    vx.put(complex(xx,yy), np);
    np++;
    if(C[k].y() > 0)				// Point should be dashed
      {
      if(pen==0)				//  If pen was solid then
        {					//  output PolyLine,
        FM_PolyLine(out, vx, ID, 15, np, 0);	//  switch pens, and begin
        pen = 10;				//  new dashed PolyLine
        vx.put(vx(np-1), 0);
        np = 1;
        }
      }
    else					// Point should be solid
      {	
      if(pen)					// If pen was dashed then
        {					// output PolyLine,
        FM_PolyLine (out, vx, ID, 15, np, 10);	// switch pens, and begin
        pen = 0;				// new solid PolyLine
        vx.put(vx(np-1), 0);
        np = 1;
        }
      }
    }
  if(np>1) FM_PolyLine(out, vx, ID,15,np,pen);	// Output any remaining
  }						// points

void FMSph::planes(ostream &out)
  {
cout << "\n\tDrawing Planes";
//  FM_Polygon(out, 298, R+1, R+1, R, 0, 15);
  plane(out, 'z', PlaneID); 			// Draw circle in xy-plane
  PlaneIDs.push_back(PlaneID);
  PlaneID++;
  plane(out, 'y', PlaneID); 			// Draw circle in xz-plane
  PlaneIDs.push_back(PlaneID);
  PlaneID++;
  plane(out, 'x', PlaneID); 			// Draw circle in yz-plane
  PlaneIDs.push_back(PlaneID);
  PlaneID++;
  }

// ____________________________________________________________________________
// iv                  FrameMaker Spherical Points
// ____________________________________________________________________________

row_vector FMSph::project(const coord_vec& data) const
  {
  double alpha = EAs.alpha();
  double beta  = EAs.beta();
  double gamma = EAs.gamma();
  double Rmax  = data.max_R();			// Largest radius in data set
  double scale = R/Rmax;			// Factor so Rmax = radius
  coord_vec datproj(data);			// Copy the data points
  datproj *= scale;				// Scale so maximum radius OK
  datproj.rotate_ip(alpha,beta,gamma);		// Rotate by Euler angles
  datproj.translate_ip(-R-1,0,-R-1);		// Translate to center points
  return datproj.project(-1,-3);		// Project 3D sphere : 2D plane
  } 						//   3D {x,z} ==> {-x,-y} 2D

void FMSph::draw(ostream& out, const row_vector& vx, int FMID) const
  {
cout << "\n\tDrawing Points";
  int npts = vx.elements();			// Points in the vector
  int i;					// Temp index
  double x1,y1,x2,y2;				// Temp coordinates
  double FMLW = 0.5;				// MIF line width
  int    FMLT = 0;				// MIF line type
  int    FMFP = 15;				// MIF fill pattern
  switch(LType)	
    {
    case 0:					// Plot points as one continuous
    default: 					// line. Max length set in PL 
      FM_PolyLine(out,vx,FMID,FMFP,npts,0);	// Output solid black line
      break;
    case 1:					// Plot points discretely
      FMLT = 1;					// Set discrete line type
      for(i=0; i<npts; i++)			// although they are grouped
        {					// (Actually small line pairs)
        double xy=0.003;			// Point variance for pair
        x1 = vx.getRe(i);			// Point x coordinate
        y1 = vx.getIm(i);			// Point y coordinate
        x2 = x1+xy;				// Offset x coordinate
        y2 = y1+xy;				// Offset y coordinate
        FM_Line(out,FMID,FMLT,FMLW,x1,y1,x2,y2);// Plot point (as small line)
        }
      break;
    case 2:					// Vectors from origin to points
      x1 = R + 1.0;				// Sphere origin x ordinate
      y1 = R + 1.0;				// Sphere origin y ordinate
      FMLW = 1.0;				// Set vector line width
      FMLT = 1;					// Set vector line type
      for(i=0; i<npts; i++)			// Loop over all points
        {
        x2 = vx.getRe(i);			// Point x ordinate 
        y2 = vx.getIm(i);			// Point y ordinate
        FM_Line(out,FMID,FMLT,FMLW,x1,y1,x2,y2);// Vector from origin to point
        }
      break;
    case 3:					// Horizontal vectors
      x1 = R+1;				// Cylinder axis x ordinate
      FMLT = 1;					// Set vector line type
      for(i=0; i<npts; i++)			// Loop over all points		
        {
        x2 = vx.getRe(i);			// Point x ordinate 
        y2 = vx.getIm(i);			// Point y ordinate
        FM_Line(out,FMID,FMLT,FMLW,x1,y2,x2,y2);// Vector z-axis to point
        }
      break;
    case 4:					// Vertical vectors 
      y1 = R+1;				// XY-plain y z-ordinate
      FMLT = 1;					// Set vector line type
      for(i=0; i<npts; i++)			// Loop over points
        {
        x2 = vx.getRe(i);			// Point x ordinate 
        y2 = vx.getIm(i);			// Point y ordinate
        FM_Line(out,FMID,FMLT,FMLW,x2,y1,x2,y2);// Vector xy-plane to point
        }
      break;
    }
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A            FrameMaker Spherical Plot Constructors, Detructor
// ____________________________________________________________________________

FMSph::FMSph()
  {
  R      = 5;					// Set sphere radius 5
  EAs    = EAngles(DefOrient);			// Set plot orientation
  NPlane = 100;					// Set 100 points in plane
  PType  = 2;					// Set plot type 2
  Marg   = 1.0;					// Set plot margins
  }

FMSph::FMSph(double Rad, const EAngles& EA) 
  {
  R = Rad;					// Set plot radius
  if(R<1 || R>10) R=5; 				// Insure plot size sensible
  EAs = EA;					// Set Euler angles
  NPlane = 100;					// Set 100 points in plane
  PType  = 2;					// Set plot type 2
  Marg   = 1.0;					// Set plot margins
  }

FMSph::FMSph(const FMSph& FMS) 
  {
  R        = FMS.R;				// Copy sphere radius 5
  EAs      = FMS.EAs;				// Copy plot orientation
  NPlane   = FMS.NPlane;			// Copy 100 points in plane
  PType    = FMS.PType;				// Copy plot type 2
  Marg     = FMS.Marg;				// Copy plot margins
  AxisIDs  = FMS.AxisIDs;			// Copy any axis ID numbers
  PlaneIDs = FMS.PlaneIDs;			// Copy any plane ID numbers
  }

// ____________________________________________________________________________
// B              FrameMaker Spherical Plot Access Functions
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                           Sphere Radius Access
// ----------------------------------------------------------------------------

double FMSph::Radius() const     { return R; }
void   FMSph::Radius(double Rad) { if(Rad>=1 && Rad<10) R=Rad; }

double FMSph::Margin() const     { return Marg; }
void   FMSph::Margin(double Mar) { if(Mar>=0.25 && Mar<3) Marg=Mar; }

// ----------------------------------------------------------------------------
//                         Sphere Orientation Access
// ----------------------------------------------------------------------------

/* Externally all single angles are expressed in degrees. Internally they are
   stored as Euler angles which are in radians. The conversion is performed
   here when needed.                                                         */

const EAngles& FMSph::Orient() const              { return EAs; }
void  FMSph::Orient(const EAngles& EA)            { EAs = EA; }
void  FMSph::Orient(double A, double B, double G) { EAs = EAngles(A,B,G,true); }

double FMSph::Alpha() const { return EAs.alpha()*RAD2DEG; }
double FMSph::Beta()  const { return EAs.beta() *RAD2DEG; }
double FMSph::Gamma() const { return EAs.gamma()*RAD2DEG; }

void   FMSph::Alpha(double A) { EAs.alpha(A*DEG2RAD); }
void   FMSph::Beta(double  B) { EAs.beta(B *DEG2RAD); }
void   FMSph::Gamma(double G) { EAs.gamma(G*DEG2RAD); }

// ----------------------------------------------------------------------------
//                         Spherical Plot Type Access
// ----------------------------------------------------------------------------

int  FMSph::PlotType() const { return PType; }
void FMSph::PlotType(int PT) { PType = PT;   }

// ----------------------------------------------------------------------------
//                      Spherical Plot Line Type Access
// ----------------------------------------------------------------------------

int  FMSph::LineType() const { return LType; }
void FMSph::LineType(int LT) { LType = LT;   }

// ----------------------------------------------------------------------------
//                      Spherical Plane Points Access
// ----------------------------------------------------------------------------

int  FMSph::PlanePts() const { return NPlane; }
void FMSph::PlanePts(int PP) { if(PP>20 && PP< 1.e4) NPlane = PP; }

// ____________________________________________________________________________
// C            FrameMaker Spherical Plot Plotting Functions
// ____________________________________________________________________________

void FMSph::plot(const string& filename)
  {
  ofstream out(filename.c_str());		// Open new file for plotting
  start(out);					// Initialize FM file
  if(abs(PType) != 1) axes(out);		// Draw coordinate axes
  if(PType >  0) planes(out);			// Draw sphere & planes
  int na = AxisIDs.size();			// Number of axes output
  int i=0;
  for(; i<na; i++)			// Loop output axes, group each
    FM_Group(out, AxisIDs[i], 280);		// as one object (segmented)
  if(na) FM_Group(out, 280, 300);		// Group all axes together
  int np = PlaneIDs.size();			// Number of planes output
  for(i=0; i<np; i++)			// Loop output planes, group
    FM_Group(out, PlaneIDs[i], 290);		// each as one object
  if(np) FM_Group(out, 290, 300);		// Group all planes together
  if(na && np)					// If both axes and planes,
   FM_Group(out, 300);				// group everything together
  end(out);					// Finalize FM file
  }

void FMSph::plot(const string& filename, const coord_vec& cvec)
  {
int PID = 250;
  ofstream out(filename.c_str());		// Open new file for plotting
  start(out);					// Initialize FM file
  if(abs(PType) != 1) axes(out);		// Draw coordinate axes
  if(PType >  0) planes(out);			// Draw sphere & planes
  row_vector vx = project(cvec); 		// Get points in plane
  draw(out, vx, PID);				// Draw points on sphere
  int na = AxisIDs.size();			// Number of axes output
  int i=0;
  for(; i<na; i++)				// Loop output axes, group each
    FM_Group(out, AxisIDs[i], 280);		// as one object (segmented)
  if(na) FM_Group(out, 280, 300);		// Group all axes together
  int np = PlaneIDs.size();			// Number of planes output
  for(i=0; i<np; i++)				// Loop output planes, group
    FM_Group(out, PlaneIDs[i], 290);		// each as one object
  if(np) FM_Group(out, 290, 300);		// Group all planes together
  FM_Group(out, 300);				// group everything together
  end(out);					// Finalize FM file
  }

void FMSph::plot(const string& filename, const vector<coord_vec>& cvs)
  {
int PID = 250;
  ofstream out(filename.c_str());		// Open new file for plotting
  start(out);					// Initialize FM file
  if(abs(PType) != 1) axes(out);		// Draw coordinate axes
  if(PType >  0) planes(out);			// Draw sphere & planes
  int nv = cvs.size();				// Get number of plots
  row_vector vx;				// Temp vectors
  int i=0;
  for(; i<nv; i++)			// Loop data sets
    {
    vx = project(cvs[i]); 			// Get points in plane
    draw(out, vx, PID++);			// Draw points on sphere
    }
  int na = AxisIDs.size();			// Number of axes output
  for(i=0; i<na; i++)			// Loop output axes, group each
    FM_Group(out, AxisIDs[i], 280);		// as one object (segmented)
  if(na) FM_Group(out, 280, 300);		// Group all axes together
  int np = PlaneIDs.size();			// Number of planes output
  for(i=0; i<np; i++)			// Loop output axes
    FM_Group(out, PlaneIDs[i], 290);		// Group all planes together
  if(np) FM_Group(out, 290, 300);		// Group all planes together
//  FM_Group(out, 298);				// Group everything together
  end(out);					// Finalize FM file
  }

// ____________________________________________________________________________
// D             FrameMaker Spherical Plot ASCII Ouptut Functions
// ____________________________________________________________________________

ostream& FMSph::print(ostream& ostr) const
  {
  string hdr("FrameMaker Spherical Plot");
  int hl = hdr.length();
  ostr << "\n" << string(40-hl/2, ' ') << hdr << "\n";
  ostr << "\n\t\tRadius:          " << R << " cm";
  ostr << "\n\t\tOrientation:     " << EAs.alpha()*RAD2DEG 
       << ", "                    << EAs.beta()*RAD2DEG
       << ", "                    << EAs.gamma()*RAD2DEG;
  ostr << "\n\t\tPoints In Plane: " << NPlane;
  ostr << "\n\t\tPlot Type:       ";
  switch(PType)
    {
    case 0:  ostr << "Cartesian Axes Only";     break;
    case 1:  ostr << "Sphere And Planes Only";  break;
    default: ostr << "Sphere, Axes And Planes"; break;
    }
  ostr << "\n\t\tPoint Display:   ";
  switch(LType)
    {
    case 0: ostr << "Discrete Points"; break;
    default:
    case 1: ostr << "Joined Points"; break;
    case 2: ostr << "Vectors From Origin"; break;
    }
  ostr << "\n";
  return ostr;
  }


// ____________________________________________________________________________
// E             NonMember FrameMaker Spherical Plot Functions
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//              Simple NonMember FrameMaker Spherical Plot Functions
// ----------------------------------------------------------------------------

/* These are intended to be very simple function calls that need few very
   arguements. Thus GAMMA users can easily generate 3D plots without any
   fretting over added details.

           Input        filename  : Output filename
                        data      : Set of {x,y,z} points
	  		Ltype     : Type of plotting to use
	  				0 = discrete points
	  				1 = joined points
	  				2 = vectors from origin
                        alpha     : Euler angle (degrees)
                        beta      : Euler angle (degrees)
                        gamma     : Euler angle (degrees)
                    or  EA        : Euler angles (radians)
                        radius    : Plot (x & y) dimension in cm
	  		Ptype      : Type of plot to produce
	  			      <=0 = Plot axes only
	  				1 = Plot sphere & planes only
	  			       >1 = Plot sphere,axes, & planes
                        points    : Number of points in plane plots
           Return                 : Void. xy-plot is output to file          */

void FrameMakerSphere(const string& name,
                         int Ptype, const coord& EA, double Rad, int N)
  {
  string fname = name + string(".mif");		// Insure file is .mif
  FMSph FMP(Rad, EA);				// Start FrameMaker 3D
  FMP.PlotType(Ptype);				// Set plot type
  FMP.PlanePts(N);				// Set plane point size
  FMP.plot(fname);				// Make MIF plot file
  }

void FrameMakerSphere(const string& name, const coord_vec& data,
                         int Ptype, const coord& EA, double Rad, int N)
  {
  string fname = name + string(".mif");		// Insure file is .mif
  FMSph FMP(Rad, EA);				// Start FrameMaker 3D
  FMP.PlotType(Ptype);				// Set plot type
  FMP.PlanePts(N);				// Set plane point size
  FMP.plot(fname, data);			// Make MIF plot file
  }

void FrameMakerSphere(const string& name, const vector<coord_vec>& data,
                         int Ptype, const coord& EA, double Rad, int N)
  {
  string fname = name + string(".mif");		// Insure file is .mif
  FMSph FMP(Rad, EA);				// Start FrameMaker 3D
  FMP.PlotType(Ptype);				// Set plot type
  FMP.PlanePts(N);				// Set plane point size
  FMP.plot(fname, data);			// Make MIF plot file
  }

// ----------------------------------------------------------------------------
//              Other NonMember FrameMaker Spherical Plot Functions
// ----------------------------------------------------------------------------

void FM_sphere(const string& filename, int Ptype,
                                 double A, double B, double G, double Rad, int N)
  {
  A *= DEG2RAD;					// Switch angles to radians
  B *= DEG2RAD;
  G *= DEG2RAD;
  EAngles EA(A,B,G);
  FMSph FMP(Rad, EA);				// Start FrameMaker 3D
  FMP.PlotType(Ptype);
  FMP.PlanePts(N);
  FMP.plot(filename);
  }

void FM_sphere(const string& filename, const coord_vec& data,
                        int Ltype, double A, double B, double G, double Rad, int N)
  {
  A *= DEG2RAD;					// Switch angles to radians
  B *= DEG2RAD;
  G *= DEG2RAD;
  EAngles EA(A,B,G);
  FMSph FMP(Rad, EA);				// Start FrameMaker 3D
  FMP.LineType(Ltype);
  FMP.PlanePts(N);
  FMP.plot(filename, data);
  }


void FM_sphere(const string& filename, coord_vec& data1, coord_vec& data2,
	 int Ltype, double A, double B, double G, double Rad, int N)
  {
  A *= DEG2RAD;					// Switch angles to radians
  B *= DEG2RAD;
  G *= DEG2RAD;
  EAngles EA(A,B,G);
  FMSph FMP(Rad, EA);				// Start FrameMaker 3D
  FMP.LineType(Ltype);
  FMP.PlanePts(N);
  vector<coord_vec> cvs;
  cvs.push_back(data1);
  cvs.push_back(data2);
  FMP.plot(filename, cvs);
  }

#endif 						// FrameMakerSphere.cc
