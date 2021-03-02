/* GgnuplotStk.cc ***********************************************-*-c++-*-
**									**
**	                           G A M M A 				**
**								 	**
**	 Gnuplot 2D/3D Stack Plots		Implementation   	**
**				        	 		 	**
**      Copyright (c) 2002						**
**      National High Magnetic Field Laboratory 	                **
**      1800 E. Paul Dirac Drive                        	        **
**      Tallahassee Florida, 32306-4005                       		**
**									**
**      $Header: $
**	                         		         	 	**
*************************************************************************/

/*************************************************************************
**								 	**
**  Description						 		**
**								 	**
**  This file contains the functions that generate 2D/3D stack plots    **
**  in Gnuplot from GAMMA vectors and matrices.                         **
**							 		**
*************************************************************************/
 
#ifndef _gnuplotStk_cc_			// Is file already included?
#  define _gnuplotStk_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamIO/GgnuplotStk.h>		// Include gnuplot header
#include <Basics/StringCut.h>		// Include GAMMA Gdec function
#include <Basics/Gconstants.h>		// Include definition of RAD2DEG
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <iostream>			// Include libstdc++ file streams
#include <string>			// Include libstdc++ ANSI strings
#include <GamIO/Ggnuplot.h>
#include <cmath>			// Inlcude HUGE_VAL_VAL

#ifdef _MSC_VER				// If we areusing MSVC++
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif

void GPStack::GPSerror(int eidx, int noret)
  {
  string hdr("Gnuplot Stack Plot Controls");
  switch(eidx)
    {
    case 10: GAMMAerror(hdr,"Cannot Find Gnuplot Executable",  noret);  //(10)
      break;
    case 50: GAMMAerror(hdr,"Cannot Plot Without A Load File",  noret); //(50)
      break;
    case 51: GAMMAerror(hdr,"Cannot Plot Without A Data File",  noret); //(51)
      break;
    default: GAMMAerror(hdr,eidx, noret); break;// Usually Unknown Error  (-1)
      break;
    }
  }

void GPStack::GPSerror(int eidx, const string& pname, int noret)
  {
  string hdr("Gnuplot Stack Plot Controls");
  string msg;
  switch(eidx)
    {
    case 10: msg = string("Executable Set To " + pname);                // (10)
             GAMMAerror(hdr, msg, noret); break;
      break;
    case 50: msg = string("Cannot Open New Load File " + pname);        // (50)
             GAMMAerror(hdr, msg, noret); break;
      break;
    default: GAMMAerror(hdr,  -1, pname, noret); break; // Unknown Error   (-1)
    }
  }

volatile void GPStack::GPSfatal(int eidx)
  { GPSerror(eidx, 1); if(eidx) GPSerror(0); cout << "\n"; exit(-1); }

volatile void GPStack::GPSfatal(int eidx, const string& pname)
  { GPSerror(eidx, pname, 1); if(eidx) GPSerror(0); cout << "\n"; exit(-1); }

// ____________________________________________________________________________
// 
// ____________________________________________________________________________

void GPStack::defaults()
  {
  GPControls::defaults();			// Set generic controls
  plottitle   = string("Stack 2D/3D Plot");	// Our default title
  xlabel      = string("X Axis");
  ylabel      = string("Y Axis");
  zlabel      = string("Z Axis");
  POVtheta    = 50.0 ;				// Our default POV theta
  POVphi      = 130.0;				// Our default POV phi
  basesphere  = true;				// Plot base sphere
  sphereaxes  = false;				// No axes output by default
  spherescale = 1.0;				// Overall scaling
  zaxisscale  = 1.0;				// Z-Axis scaling
  degrees     = false;				// Output in degrees
  spherical   = true;				// Mapping is spherical     
  hidden      = false;				// Hidden is off
  parametric  = true;				// Set parametric on
  dataout     = false;
  normalize   = true;
  }

void GPStack::copy(const GPStack& GPS)
  {
  GPControls::copy(GPS);				// Set generic controls
  plottitle   = GPS.plottitle;				// Our default title
  POVtheta    = GPS.POVtheta;				// Our default POV theta
  POVphi      = GPS.POVphi;				// Our default POV phi
  basesphere  = GPS.basesphere;				// Base sphere plot flag
  sphereaxes  = GPS.sphereaxes;				// No axes output by default
  spherescale = GPS.spherescale;			// Overall scaling
  zaxisscale  = GPS.zaxisscale;				// Z-Axis scaling
  degrees     = GPS.degrees;				// Output in degrees
  spherical   = GPS.spherical;				// Mapping is spherical     
  hidden      = GPS.hidden;				// Hidden is off
  parametric  = GPS.parametric;				// Set parametric on
  dataout     = GPS.dataout;				// Data output flag
  normalize   = GPS.normalize;
  }

// ____________________________________________________________________________
// iii             Gnuplot Stack 2D/3D Plot Loadfile Output
// ____________________________________________________________________________

void GPStack::WriteTitle(ofstream& ofstr)
  {
  if(plottitle.length())
    ofstr << "set title \"" << plottitle << "\"\n"; 
  }

void GPStack::WriteAngles(ofstream& ofstr)
  {
  if(degrees)
    ofstr << "set angles degrees\n";
  }

void GPStack::WriteMapping(ofstream& ofstr)
  {
  if(spherical)
    ofstr << "set mapping spherical\n";
  else
    ofstr << "set mapping cylindrical\n";
  }

void GPStack::WriteView(ofstream& ofstr)
  {
  ofstr << "set view ";
  ofstr << POVtheta    << ", ";
  ofstr << POVphi      << ", ";
  ofstr << spherescale << ", ";
  ofstr << zaxisscale  << "\n";
  }

void GPStack::WriteSplot(ofstream& ofstr)
  {
  ofstr << "splot ";
  bool newline = false;
  if(basesphere)
    {
    ofstr << "cos(u)*cos(v),"
          << "cos(u)*sin(v),"
          << "sin(u)";
    newline = true;
    }
  if(dataout)
    {
    if(newline) ofstr << ",\\" << endl;
    if(datafile.length())
      {
      ofstr << "\'" << datafile << "\'";
      newline = true;
      }
    int N = datafiles.size();
    if(N && newline) { ofstr << ",\\" << endl; newline=false; }
    for(int i=0; i<N; i++)
      {
      if(newline) ofstr << ",\\";
      ofstr << ",\\\n\'" << datafiles[i] << "\'";
      newline = true;
      }
    }
  if(sphereaxes)
    {
    if(newline) ofstr << ",\\" << endl;
    ofstr << "\'xaxis.asc\', \'yaxis.asc\', \'zaxis.asc\'";
    }
  ofstr << endl;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A      Gnuplot Stack 2D/3D Controls Contruction, Assigment, Destruction
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

GPStack::GPStack() { defaults(); }
GPStack::GPStack(const GPStack& GPC) { copy(GPC);  }

GPStack::~GPStack() {}
GPStack& GPStack::operator= (const GPStack& GPC)
  {
  if(this == &GPC) return *this;
  copy(GPC);
  return *this;
  }

// ____________________________________________________________________________
// B              Gnuplot Stack 2D/3D Control Access
// ____________________________________________________________________________

void GPStack::SetDegrees(bool    dg)    { degrees    = dg; }
void GPStack::SetStackAxes(bool ax)    { sphereaxes = ax; }
void GPStack::SetBaseStack(bool bs)    { basesphere = bs; }
void GPStack::SetNormalization(bool nm) { normalize  = nm; }

// ____________________________________________________________________________
// C              Gnuplot Stack 2D/3D Loadfile Generation
// ____________________________________________________________________________

void GPStack::LoadFile()
  {
  if(!NewLoadFile()) GPSfatal(50); 		// Begin new load file
  WriteView(Lfp);				// Set POV (set view)
  WriteTitle(Lfp);				// Set plot title
  WriteXLabel(Lfp);				// Set x-axis label
  WriteYLabel(Lfp);				// Set y-axis label
  WriteZLabel(Lfp);				// Set z-axis label
  WriteAngles(Lfp);				// Set angle units
  Lfp  << "set data style lines\n";	
  WriteParametric(Lfp);				// Set parametric flag
  WriteKey(Lfp);				// Set legend
  WriteMapping(Lfp);				// Set mapping
//  WriteRange(Lfp);				// Write the 
  Lfp << "set samples 32\n";			// This may long. lat. lines
  Lfp << "set isosamples 9, 20\n";		// This many lat/log lines
  Lfp << "set urange [-pi/2:pi/2]\n";		// Plotting range
  Lfp << "set vrange [0:2*pi]\n";
  WriteSplot(Lfp);				// Issue plot command
  WritePause(Lfp);
  CloseLoadFile();				// Close load file
  }

// ____________________________________________________________________________
// D                  Gnuplot Stack 2D/3D Data File Production
// ____________________________________________________________________________

/* The ASCII data is set into three columns { theta, phi, R }. However, Gnuplot
   takes theta as the angle over from the +x axis - that which GAMMA calls the
   phase angle phi. Gnuplot takes the angle up/down from the XY plane as phi,
   which also contrasts with GAMMA in that we take theta as the angle down from
   the +z axis. So, we need to output {PHI, 90-THETA, R} to keep Gnuplot happy.
   Lastly, we can output out data in either degrees or radians, the later 
   being the default.                                                        */

bool GPStack::DataAxisFile(char axis, bool warn)
  {
  string fileout;
  switch(axis)
    {
    default:
    case 'x': fileout = "xaxis.asc"; break;
    case 'y': fileout = "yaxis.asc"; break;
    case 'z': fileout = "zaxis.asc"; break;
    }
  ofstream ofstr;
  ofstr.open(fileout.c_str(), ios::out);	// Open loadfile file
  if(!ofstr.is_open())				// Insure file is OK
    {
    if(warn) 
      {
      GPCerror(50, loadfile, 1);		//   Problems with file
      GPSfatal(51);				// Begin new data file	
      }
    return false;
    }
  double x1=0, y1=0, z1=1;
  double fact = degrees?1.0:DEG2RAD;
  switch(axis)
    {
    default:
    case 'x': break;
    case 'y': x1=fact*90.0; break;
    case 'z': y1=fact*90.0; break;
    }
  ofstr << 0.0 << delim << 0.0 << delim << 0.0 << endl;
  ofstr << x1  << delim << y1  << delim << z1  << endl;
  ofstr.close();
  return true;
  }

void GPStack::DataFile(const row_vector& data, int idx)
  {
  if(!NewDataFile(idx)) GPSfatal(51);		// Begin new data file
  int N         = data.size();			// Length of coord vector
  double fact   = (degrees)?RAD2DEG:1.0;	// Angle conversion factor
  double ninety = PI/2.0;			// 90 degrees in radians
  double COL1, COL2, COL3;			// For output columns
  double Rnorm = 1.0;
  if(normalize) Rnorm = data.max_R();
  coord pt;					// Temporary point
  for(int l=0; l<N; l++)			// Loop points to plot
    {
    pt   = data.get(l);				//   Current point
    COL1 = fact*pt.phi();			//   Point phi (rad/deg)
    COL2 = fact*(ninety-pt.theta()); 		//   Point theta (rad/deg)
    COL3 = pt.Rad()/Rnorm;			//   Point radius
    Dfp << COL1 << delim			//   Write line output stream
        << COL2 << delim			//   in ASCII format, three
        << COL3 << endl;			//   columns dileneated
    }
  CloseDataFile();				// Close data file
  }

void GPStack::DataFiles(const vector<row_vector>& data)
  {
  int N = data.size();				// Number of vectors
  for(int i=0; i<N; i++)			// Loop over each vector
    DataFile(data[i], i);			// Write into same data file
  }

// ____________________________________________________________________________
// E             Gnuplot Stack 2D/3D Plot Output Generation 
// ____________________________________________________________________________

/* These functions will generate the Gnuplot load file and one or more ASCII
   data files that will be plotted. Subsequently, the command is issued to
   run gnuplot using the produced load file which reads the data files.      */

void GPStack::Plot()
  {
  if(sphereaxes)				// If coordinate axes desired
    {
    DataAxisFile('x');
    DataAxisFile('y');
    DataAxisFile('z');
    }
  dataout = false;				// Set for not data out
  LoadFile();					// Create load file
  RunLoadFile();				// Run Gnuplot w/ load file
  }

void GPStack::Plot(const coord_vec& data)
  {
  DataFile(data);				// Create data file
  dataout = true;				// Set output true
  LoadFile();					// Create load file
  RunLoadFile();				// Run Gnuplot w/ load file
  }

void GPStack::Plot(const vector<coord_vec>& data)
  {
  DataFiles(data);				// Create data file
  dataout = true;				// Set output true
  LoadFile();					// Create load file
  RunLoadFile();				// Run Gnuplot w/ load file
  }


// ____________________________________________________________________________
// F         Gnuplot Stack 2D/3D Controls Standard Output Functions
// ____________________________________________________________________________


ostream& GPStack::print(ostream& ostr)
  {
  string hdr("Gnuplot Stack 2D/3D Plotting Controls");	// Output header
  int hl = hdr.length();                                // Header length
  ostr << "\n" << string(40-hl/2, ' ') << hdr;
  string On("On");
  string Off("Off");
  ostr << "\n";
  ostr << "\n\tData File Name:            " << datafile;
  ostr << "\n\tLoad File Name:            " << loadfile;
  ostr << "\n\tPlot Title:                " << plottitle;
  ostr << "\n\tView POV Angle From +X:    " << POVtheta << " Degrees";
  ostr << "\n\tView POV Angle Up From XY: " << POVphi   << " Degrees";
  ostr << "\n\tOverall Plot Scaling:      " << spherescale;
  ostr << "\n\tZ-Axis Scaling:            " << zaxisscale;
  ostr << "\n\tDraw Key Base Stack:      "; if(basesphere) ostr << On;
                                             else           ostr << Off;
  ostr << "\n\tDraw Stack Axes:          "; if(sphereaxes) ostr << On;
                                             else           ostr << Off;

  ostr << "\n\tDraw Key:                  "; if(nokey)     ostr << Off;
                                             else          ostr << On;
  ostr << "\n\tAngles Output In:          "; if(degrees)   ostr << "Degrees";
                                             else          ostr << "Radians";
  ostr << "\n\tMapping Is:                "; if(spherical) ostr << "Spherical";
                                             else          ostr << "Cylindrical";
  ostr << "\n\tHidden Flag Is:            "; if(hidden)    ostr << "On";
                                             else          ostr << "Off";
  ostr << "\n\tNormalization Is:          "; if(normalize) ostr << "On";
                                             else          ostr << "Off";
  return ostr;
  }

ostream& operator<< (ostream& ostr, GPStack& GPS)
  { return GPS.print(ostr); }

// ____________________________________________________________________________ 
// G                          Non-Member Functions
// ____________________________________________________________________________
 
/* These are functions built prior to construction of this class. Some of
   them are deprecated and will issue warnings if utilized.                  */
 

/*************************************************************************
**			       gnuplot Stack Plots			**
**									**
** Author       : S.A. Smith						**
** Date		: March 5, 1994						**
** Description  : These routines are intended to output data files in	**
**		  ASCII format which are readable by the public domain	**
**		  program gnuplot.  Their purpose is to provide an easy **
**		  means to view stack plots generated in GAMMA.		**
** Notes        : gnuplot takes ASCII files filled with {x,y,z}	points	**
**		  as input for 3D plots.  Only one point per line is	**
**		  allowed, the three ordinates separated by spaces. The **
**		  x & y ordinates are optional but are output in these	**
**		  routines. For stack plots, the output ASCII file must	**
**		  have a blank line following each row. The plots	**
**		  themselves will be displayed as follows:		**
**									**
**			[-x, +x] from left to right			**
**			[-y, +y] from front to back			**
**			[-z, +z] from top to bottom			**
**									**
**		  Useful gnuplot commands for stack plots are		**
**									**
**			set parametric       (independent rows)		**
**			set data style lines (connected points)		**
**			splot "filename"     (draw stack plot)		**
**									**
**		  Because gnuplot ASCII files of this type can become	**
**		  very large, a data reduction scheme is utilized	**
**		  based on the change in the zcoordinate within each	**
**		  row. If there is no change relative to a cutoff value	**
**		  the point is skipped in the output file.		**
**									**
*************************************************************************/


	// Input	ofstr     : Output filestream
	// 		vx        : Data vector
	//		row       : Row value or y-axis value
        //		ri	  : Flag for real versus imaginary plot
	//			    ri = 0  : reals only
	//			    ri = <0 : imaginaries only
	//			    ri = >0 : both 
        //              xmin      : Plot horizontal (x=0) point value
        //              xmax      : Plot horizontal (x=end) point value
	//		cutoff    : Data compression cutoff
	// Return		  : Void, the information in vector vx is
	//			    output to ofstr in ASCII format suitable
	//			    for use in gnuplot.  Each data point is
	//			    written as {x,y,z} where x is the point
	//			    number, y the row number, and z the vector
	//			    point value 
	// Note			  : This is essentially the same function as
	//			    GP_1D except output is {x,y,z} not {x,y}

void GP_stack(ofstream& ofstr, const row_vector &vx, double row,
                              int ri, double xmin, double xmax, double cutoff)
  {
  int np = vx.size();				// Number of points in vector
  int i = 0;					// Vector point count
  int ptaxis = 0;				// Assume x ordinates not output
  double delx = 0;				// If they are to be output,
  if(xmin != xmax)				// then set the ptaxis flag
    {						// & figure out what the x-axis
    ptaxis = 1;					// increment is
    delx = (xmax-xmin)/double(np-1);
    }
  double pt, lpt=HUGE_VAL;				// Set up for data reduction
  double ymax=HUGE_VAL,ymin=HUGE_VAL;			// (done if ordinates output)
  if(ptaxis && !cutoff)				// If ordinates output, find max
    {						// & min if no compress. cutoff
    if(ri >= 0)					// Search the real points
      {
      pt = Re(vx.get(i));
      if(ymax < pt) ymax = pt;
      if(ymin > pt) ymin = pt;
      }
    if(ri)					// Search the imaginary points
      {
      pt = Im(vx.get(i));
      if(ymax < pt) ymax = pt;
      if(ymin > pt) ymin = pt;
      }
    }
  if(!cutoff) cutoff = 1.e-3*(ymax-ymin);	// Cutoff/ data redux (resolut.)
  int skip = 0;
  if(ri >= 0)                                   // Write data for real plot
    for(i=0; i<np; i++)				//	Loop over all the points
      {
      if(ptaxis)				//	This if an x-axis needed
        {					//	We can use data reduct.
        pt = Re(vx.get(i));			//	The point (intensity)
        if(fabs(pt-lpt) > cutoff || i==np-1)	//	Plot if y differs/last
          {
          if(skip)				// 	If the previous point was not
            {					//	plotted, we must plot it now
            ofstr<<xmin+double(i-1)*delx << "  "//	or else plot won't be smooth
                  << row << " " 
                  << lpt << "\n";
            }
          ofstr<<xmin+double(i)*delx << "  "	//   First output x axis value
                << row << " " 			//   Now output y axis value
                << pt << "\n";			//   Now output intensity value
          lpt = pt;				//   Store as last plotted point
          skip = 0;				//   Flag last pt wasn't skipped
          }
        else					//   If a pt was skipped, then
          skip = 1;				//   set the skip flag
        }
      else					//  If just y axis (default)
        {					//  Must plot all pts this way
        ofstr << Gdec(i) << " " 		//  cause x vals assumed evenly 
              << row << " " 			//  incremented
              << vx.getRe(i) << "\n";
        }
      }
  if(ri > 0)                                    // Blank line to separate plots
      ofstr << "\n";				// if both real & imags output
  lpt = HUGE_VAL;
  skip = 0;
  if(ri)                                        // Write data for imaginary plot
    for(i=0; i<np; i++)
      {
      if(ptaxis)
        {
        pt = Im(vx.get(i));
        if(fabs(pt-lpt) > cutoff || i==np-1)
          {
          if(skip)
            {
            ofstr << xmin + double(i-1)*delx << "  "
                  << row << " " 
                  << lpt << "\n";
            }
          ofstr << xmin + double(i)*delx << "  "
                << row << " " 
                << pt << "\n";
          lpt = pt;
          skip = 0;
          }
        else
          skip = 1;
        }
      else
        {
        ofstr << Gdec(i) << " " 
              << row << " " 
              << Im(vx.get(i)) << "\n";
        }
      }
  ofstr << "\n";			// Another empty line, assume
  }					// more rows will follow



	// Input	filename  : Output file name
	// 		mx        : Data matrix
        //		ri	  : Flag for real versus imaginary plot
	//			    ri = 0  : reals only
	//			    ri = <0 : imaginaries only
	//			    ri = >0 : both 
        //              ymin      : Plot "vertical" (y=0) point value
        //              ymax      : Plot "vertical" (y=end) point value
        //              xmin      : Plot horizontal (x=0) point value
        //              xmax      : Plot horizontal (x=end) point value
	// Return		  : Void, the information in vector mx is
	//			    output to file filename in an ASCII
	//			    format suitable for use in gnuplot.
	//			    Each data point is written as {x,y,z}
	//			    where x is the point number, y the row 
	//			    number, and z the vector point value 

void GP_stack(const string& filename, matrix& mx, int ri,
                            double ymin, double ymax, double xmin, double xmax)
  {
  ofstream ofstr(filename.c_str());		// Begins an output file stream
  GP_stack(ofstr,mx,ri,ymin,ymax,xmin,xmax);	// Use GP_stack overload
  ofstr.close();				// Close the file stream
  return;
  }
  


	// Input	ofstr     : An output file stream
	// 		mx        : Data matrix
        //		ri	  : Flag for real versus imaginary plot
	//			    ri = 0  : reals only
	//			    ri = <0 : imaginaries only
	//			    ri = >0 : both 
        //              ymin      : Plot "vertical" (y=0) point value
        //              ymax      : Plot "vertical" (y=end) point value
        //              xmin      : Plot horizontal (x=0) point value
        //              xmax      : Plot horizontal (x=end) point value
	// Return		  : Void, the information in vector mx is
	//			    output to file filename in an ASCII
	//			    format suitable for use in gnuplot.
	//			    Each data point is written as {x,y,z}

void GP_stack(ofstream& ofstr, matrix& mx, int ri,
                            double ymin, double ymax, double xmin, double xmax)
  {
  int nc = mx.cols();
  int nr = mx.rows();
  int i, j;
  double zmin=HUGE_VAL, zmax=-HUGE_VAL;
  double pt;
  for(i=0; i<nr; i++) 				// Determine a suitable scaling
    for(j=0; j<nc; j++) 			// for data reduction by first 
      { 					// finding global maxima
      if(ri >= 0)				// Search the real points
        {
        pt = Re(mx.get(i,j));
        if(zmax < pt) zmax = pt;
        if(zmin > pt) zmin = pt;
        }
      if(ri)					// Search the imaginary points
        {
        pt = Im(mx.get(i,j));
        if(zmax < pt) zmax = pt;
        if(zmin > pt) zmin = pt;
        }
      }
  double cutoff = 1.e-3*(zmax-zmin);		// For data reduction (resolut.)
  row_vector vx;
  double y = 0;					// Set up y-axis scaling
  double dely = 0;				// if any has been specified
  int ptaxis = 0;				// Assume none has been set
  if(ymin != ymax)				// If so, then flag that it is
    {						// & figure out what the y-axis
    ptaxis = 1;					// increment is
    dely = (ymax-ymin)/double(nr-1);
    }
  for(i=0; i<nr; i++)				// Loop through the matrix rows
    {
    vx = mx.get_block(i,0,1,nc);		// Put this row as a row_vector		
    if(ptaxis)					// If non-default y-axis values
      {
      y = ymin + double(i)*dely;		// Get the y-axis value
      GP_stack(ofstr,vx,y,ri,xmin,xmax,cutoff);	// Use GP_stack overload
      }
    else					// If default y-axis values
      GP_stack(ofstr,vx,	 		// Use GP_stack overload, type
               (double)i,ri,xmin,xmax,cutoff);  // i so right fcn is used
    }
  }
 

//------------------------------ Common Stack Plots ---------------------------


void GP_stackplot(const string& gnumacro, const string& file2D)
  {
  int join=1;
  ofstream gnuload(gnumacro.c_str());		// File of gnuplot commands
  SetLineType(gnuload, join);			// Set plot line type
  gnuload << "set parametric\n";		// Plot x vs y as parametric
  gnuload << "splot \"" << file2D << "\"\n";	// Command for gnuplot plot
  CloseMacro(gnuload);				// Add pause+exit, close macro
  RunGnuplot(gnumacro);				// Run gnuplot with gnumacro
  }




#endif 							// GnuplotStk.cc
