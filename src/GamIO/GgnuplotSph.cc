/* GgnuplotSph.cc ***********************************************-*-c++-*-
**									**
**	                           G A M M A 				**
**								 	**
**	 Gnuplot 3D Spherical Plots		Implementation   	**
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
**  This file contains the functions that generate 3D spherical plots   **
**  in Gnuplot from GAMMA coordinate vectors. Normally the coordinate   **
**  vector will contain Cartesian coordinates. These are mapped into    **
**  three dimension upon output into a Gnuplot compatible ASCII file.   **
**							 		**
*************************************************************************/
 
#ifndef   GgnuplotSph_cc_		// Is file already included?
#  define GgnuplotSph_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamIO/GgnuplotSph.h>		// Include gnuplot header
#include <Basics/StringCut.h>		// Include GAMMA Gdec function
#include <Basics/Gconstants.h>		// Include definition of RAD2DEG
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <iostream>			// Include libstdc++ file streams
#include <string>			// Include libstdc++ ANSI strings
#include <vector>			// Include libstdc++ STL vectors
#include <GamIO/Ggnuplot.h>

using std::string;			// Using libstdc++ strings
using std::vector;			// Using libstdc++ vectors
using std::ostream;			// Using libstdc++ output streams
using std::ofstream;			// Using libstdc++ output file streams
using std::ios;				// Using libstdc++ file type settings
using std::endl;			// Using libstdc++ line end

void GPSphere::GPSerror(int eidx, int noret)
  {
  string hdr("Gnuplot Spherical Plot Controls");
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

void GPSphere::GPSerror(int eidx, const string& pname, int noret)
  {
  string hdr("Gnuplot Spherical Plot Controls");
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

volatile void GPSphere::GPSfatal(int eidx)
  { GPSerror(eidx, 1); if(eidx) GPSerror(0); GAMMAfatal(); }

volatile void GPSphere::GPSfatal(int eidx, const string& pname)
  { GPSerror(eidx, pname, 1); if(eidx) GPSerror(0); GAMMAfatal(); }

// ____________________________________________________________________________
// 
// ____________________________________________________________________________

void GPSphere::defaults()
  {
  GPControls::defaults();			// Set generic controls
  plottitle   = string("Spherical 3D Plot");	// Our default title
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

void GPSphere::copy(const GPSphere& GPS)
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
// iii             Gnuplot Spherical 3D Plot Loadfile Output
// ____________________________________________________________________________

void GPSphere::WriteTitle(ofstream& ofstr)
  {
  if(plottitle.length())
    ofstr << "set title \"" << plottitle << "\"\n"; 
  }

void GPSphere::WriteAngles(ofstream& ofstr)
  {
  if(degrees)
    ofstr << "set angles degrees\n";
  }

void GPSphere::WriteMapping(ofstream& ofstr)
  {
  if(spherical)
    ofstr << "set mapping spherical\n";
  else
    ofstr << "set mapping cylindrical\n";
  }

void GPSphere::WriteView(ofstream& ofstr)
  {
  ofstr << "set view ";
  ofstr << POVtheta    << ", ";
  ofstr << POVphi      << ", ";
  ofstr << spherescale << ", ";
  ofstr << zaxisscale  << "\n";
  }

void GPSphere::WriteSplot(ofstream& ofstr)
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
// A      Gnuplot Spherical 3D Controls Contruction, Assigment, Destruction
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

GPSphere::GPSphere() { defaults(); }
GPSphere::GPSphere(const GPSphere& GPC) { copy(GPC);  }

GPSphere::~GPSphere() {}
GPSphere& GPSphere::operator= (const GPSphere& GPC)
  {
  if(this == &GPC) return *this;
  copy(GPC);
  return *this;
  }

// ____________________________________________________________________________
// B              Gnuplot Spherical 3D Control Access
// ____________________________________________________________________________

void GPSphere::SetDegrees(bool    dg)    { degrees    = dg; }
void GPSphere::SetSphereAxes(bool ax)    { sphereaxes = ax; }
void GPSphere::SetBaseSphere(bool bs)    { basesphere = bs; }
void GPSphere::SetNormalization(bool nm) { normalize  = nm; }

// ____________________________________________________________________________
// C              Gnuplot Spherical 3D Loadfile Generation
// ____________________________________________________________________________

void GPSphere::LoadFile()
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
// D             Gnuplot Spherical 3D Data File Production
// ____________________________________________________________________________

/* The ASCII data is set into three columns { theta, phi, R }. However, Gnuplot
   takes theta as the angle over from the +x axis - that which GAMMA calls the
   phase angle phi. Gnuplot takes the angle up/down from the XY plane as phi,
   which also contrasts with GAMMA in that we take theta as the angle down from
   the +z axis. So, we need to output {PHI, 90-THETA, R} to keep Gnuplot happy.
   Lastly, we can output out data in either degrees or radians, the later 
   being the default.                                                        */

bool GPSphere::DataAxisFile(char axis, bool warn)
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

void GPSphere::DataFile(const coord_vec& data, int idx)
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

void GPSphere::DataFiles(const vector<coord_vec>& data)
  {
  int N = data.size();
  for(int i=0; i<N; i++)
    {
    DataFile(data[i], i);
    }
  }

// ____________________________________________________________________________
// E             Gnuplot Spherical 3D Plot Output Generation 
// ____________________________________________________________________________

/* These functions will generate the Gnuplot load file and one or more ASCII
   data files that will be plotted. Subsequently, the command is issued to
   run gnuplot using the produced load file which reads the data files.      */

void GPSphere::Plot()
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

void GPSphere::Plot(const coord_vec& data)
  {
  DataFile(data);				// Create data file
  dataout = true;				// Set output true
  LoadFile();					// Create load file
  RunLoadFile();				// Run Gnuplot w/ load file
  }

void GPSphere::Plot(const vector<coord_vec>& data)
  {
  DataFiles(data);				// Create data file
  dataout = true;				// Set output true
  LoadFile();					// Create load file
  RunLoadFile();				// Run Gnuplot w/ load file
  }


// ____________________________________________________________________________
// F         Gnuplot Spherical 3D Controls Standard Output Functions
// ____________________________________________________________________________


ostream& GPSphere::print(ostream& ostr)
  {
  string hdr("Gnuplot Spherical 3D Plotting Controls");	// Output header
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
  ostr << "\n\tDraw Key Base Sphere:      "; if(basesphere) ostr << On;
                                             else           ostr << Off;
  ostr << "\n\tDraw Sphere Axes:          "; if(sphereaxes) ostr << On;
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

ostream& operator<< (ostream& ostr, GPSphere& GPS)
  { return GPS.print(ostr); }

// ____________________________________________________________________________ 
// G                          Non-Member Functions
// ____________________________________________________________________________
 
/* These are functions built prior to construction of this class. Some of
   them are deprecated and will issue warnings if utilized.                  */
 
void GP_sphere(const string& name, const coord_vec& data, bool old)
  {
  GPSphere GPS;
  GPS.SetDataFile(name);
  if(old) GPS.SetDegrees(true);			// Set output in degrees
  GPS.DataFile(data);				// Generate data file
  }

void GP_sphereplot(const string& gnumacro, const string& Aname, bool old)
  {
  GPSphere GPS;
  GPS.SetDataFile(Aname);
  GPS.SetLoadFile(gnumacro);
  if(old) GPS.SetDegrees(true);			// Set output in degrees
  GPS.LoadFile();				// Create load file
  GPS.RunLoadFile();				// Run Gnuplot w/ load file
  }

void GP_sphereplot(const string& gnumacro, int N, string* files)
  {
  GPSphere GPS;
  GPS.SetLoadFile(gnumacro);
  vector<string> dfs;
  for(int i=0; i<N; i++)
    dfs.push_back(files[i]);
  GPS.SetDataFiles(dfs);
  GPS.LoadFile();					// Create load file
  GPS.RunLoadFile();				// Run Gnuplot w/ load file
  }

void GnuplotSphere(const string& name, const coord_vec& data)
  {
  GPSphere GPS;
  GPS.SetDataFile(name+".asc");
  GPS.SetLoadFile(name+".gnu");
  GPS.Plot(data);
  }

#endif 						// GnuplotSph.cc
