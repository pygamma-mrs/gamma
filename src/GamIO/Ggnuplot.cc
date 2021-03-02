/* Ggnuplot.cc **************************************************-*-c++-*-
**									**
**	                           G A M M A 				**
**								 	**
**	Gnuplot Support 	         	   Implementation   	**
**				        	 		 	**
**      Copyright (c) 1994, 1998					**
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
**  This is a still rather crude attempt to create an interface between	**
**  GAMMA and the public domain package called gnuplot.  Even though 	**
**  the functions herein work, they are far from what one would call	**
**  optimal!  As things progress they should work better, and perhaps	**
**  this will be made into a class which supports gnuplot command files	**
**  some day.  For now, it is stupid and just puts out data files.	**
**							 		**
*************************************************************************/
 
#ifndef   Ggnuplot_cc_			// Is file already included?
#  define Ggnuplot_cc_ 	1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamGen.h>			// Include OS specific stuff
#include <GamIO/Ggnuplot.h>		// Include gnuplot header
#include <Basics/StringCut.h>		// Include GAMMA Gdec function
#include <Basics/Gconstants.h>		// Include definition of RAD2DEG
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <stdlib.h>
#include <iostream>			// Include libstdc++ file streams
#include <fstream>
#include <string>			// Include libstdc++ ANSI strings
#include <vector>			// Include libstdc++ STL vectors
#include <cstdio>			// Include FILE
#include <cmath>			// Inlcude HUGE_VAL_VAL

using std::string;			// Using libstdc++ strings
using std::vector;			// Using libstdc++ vectors
using std::cout;			// Using libstdc++ standard output
using std::ifstream;

// ____________________________________________________________________________
// i                      GAMMA Gnuplot Error Handling
// ____________________________________________________________________________

/*         Input                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  

void GPerror(int eidx, int noret)
  {
  string hdr("GAMMA Gnuplot Routine");
  switch(eidx)
    {
    case 10: GAMMAerror(hdr,"Cannot Find Gnuplot Executable",  noret);  //(10)
      break;
    case 11: GAMMAerror(hdr,"Stopping Program, No Executable", noret);  //(11)
      break;
    case 12: GAMMAerror(hdr,"Execution Shell Unavailable?",    noret);  //(12)
      break;
    case 13: GAMMAerror(hdr,"Cannot Run Gnuplot Interactively",noret);  //(13)
      break;
    case 14: GAMMAerror(hdr,"Put Command \"gnuplot\" In PATH", noret);  //(14)
      break;
    default: GAMMAerror(hdr,eidx, noret); break;// Usually Unknown Error  (-1)
      break;
    }
  }

void GPerror(int eidx, const string& pname, int noret)
  {
  string hdr("GAMMA Gnuplot Routine");
  string msg;
  switch(eidx)
    {
    case 10: msg = string("Executable Set To " + pname);                // (10)
             GAMMAerror(hdr, msg, noret); break;
    case 11: msg = string("System Call To " + pname + " Returns 0");    // (11)
             GAMMAerror(hdr, msg, noret); break;
    case 12: msg = string("System Call " + pname + " Fails");           // (12)
             GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr,  -1, pname, noret); break; // Unknown Error   (-1)
    }
  }

volatile void GPfatality(int eidx)
  { GPerror(eidx, 1); if(eidx) GPerror(0); GAMMAfatal(); }

volatile void GPfatality(int eidx, const string& pname)              
  { GPerror(eidx, pname, 1); if(eidx) GPerror(0); GAMMAfatal(); }
 

//_____________________________________________________________________________
// A                     GAMMA Gnuplot System Interface
//_____________________________________________________________________________

        // Input        none    : Output filename
	//		warn    : Warning flag
        // Return       GPE     : A string for the Gnuplot
        //                        executable command (hopefully)
	// Note         	: Add more Gnuplot paths to GPNames in
	//			  GPSeekStrings as needed.

vector<string> GPSeekStrings()
  {
  vector<string> GPNames;
  string GP("gnuplot");
  string WGP("wgnuplot");
  string UB("/usr/bin/");
  string ULB("/usr/local/bin/");
  string UFW("/usr/freeware/bin/");
  string CW("/cygwin/bin/");
  string CWU("/cygwin/usr/bin");
  string CWUL("/cygwin/usr/local/bin");
  string MW("/mingw");
  string MWB("/mingw/bin");
  GPNames.push_back(GP);			// If command local
  GPNames.push_back(ULB + GP);			// GNUish install spot
  GPNames.push_back(UB + GP);			// Linux RPM install spot
  GPNames.push_back(UFW + GP);			// SGI tardist install spot
  GPNames.push_back(WGP);			// If command local Windoze
  GPNames.push_back(ULB + WGP);			// Gnuish for cygwin
  GPNames.push_back(UB + WGP);			// Why not for cygwin
  #if defined(_MSC_VER) || defined(__MINGW32__)
    string DC("C:");
    string DD("D:");
    GPNames.push_back(DC + ULB  + GP);
    GPNames.push_back(DC + CW   + GP);
    GPNames.push_back(DC + CWU  + GP);
    GPNames.push_back(DC + CWUL + GP);
    GPNames.push_back(DC + ULB  + WGP);
    GPNames.push_back(DC + CW   + WGP);
    GPNames.push_back(DC + CWU  + WGP);
    GPNames.push_back(DC + CWUL + WGP);
    GPNames.push_back(DD + ULB  + GP);
    GPNames.push_back(DD + CW   + GP);
    GPNames.push_back(DD + CWU  + GP);
    GPNames.push_back(DD + CWUL + GP);
    GPNames.push_back(DD + ULB  + WGP);
    GPNames.push_back(DD + CW   + WGP);
    GPNames.push_back(DD + CWU  + WGP);
    GPNames.push_back(DD + CWUL + WGP);
  #endif

  return GPNames;
  }

string GPFind(bool vocal)
  {
  vector<string> GPNames = GPSeekStrings();	// Gnuplot path names
  int N = GPNames.size();			// Number of path names
  string GPE, GPName;				// Gnuplot execution command
  bool found = false;				// Flag if found executable
  ifstream gexec;           // Gnuplot executable?
  int i=0;
  for(; i<N && !found; i++)
    {
    GPE = GPNames[i]; 				// Try this name first
    if(vocal)
      cout << "\n\t Seeking Gnuplot Executable " << GPE;

    //gexec = fopen(GPE.c_str(), "rb");		// Does this file exist?
    //if(gexec != NULL) found=true;			// If so, use this command
	  gexec.open(GPE.c_str(), ifstream::binary);
    if(gexec.good()== true)
			found = true;			

    if(vocal)
    { if(found) 
				cout << " - Success!"; 
      else      
				cout << " - Not Found"; 
    }
    }

  if(found) 
	  { 
	  gexec.close(); 
	  return GPE; 
	  } 

  for(i=0; i<N && !found; i++)
    {
    GPE = GPNames[i] + string(".exe"); 				// Then try this name
    if(vocal)
      cout << "\n\t Seeking Gnuplot Executable " << GPE;

    gexec.open(GPE.c_str(), ifstream::binary);
    if(gexec.good()== true)
	found = true;	
    
    if(vocal)
      {
      if(found) 
				cout << " - Success!"; 
      else      
				cout << " - Not Found"; 
      }
    }

  if(found) 
	  { 
		gexec.close(); 
		return GPE; 
	  } 

  for(i=0; i<N && !found; i++)
    {
    GPE = GPNames[i] + string(".out"); 				// Try this name first
    if(vocal)
      cout << "\n\t Seeking Gnuplot Executable " << GPE;

	  gexec.open(GPE.c_str(), ifstream::binary);
    if(gexec.good()== true)
			found = true;	

    if(vocal)
    { if(found) 
				cout << " - Success!"; 
      else      
				cout << " - Not Found"; 
    }
    }
  if(found) 
	  { 
	  gexec.close(); 
		return GPE; 
	  } 

  return string("");
  }

/* saved version *
string GPFind(bool vocal)
  {
  vector<string> GPNames = GPSeekStrings();	// Gnuplot path names
  int N = GPNames.size();			// Number of path names
  string GPE, GPName;				// Gnuplot execution command
  bool found = false;				// Flag if found executable
  FILE* gexec = NULL;				// Gnuplot executable?
  int i=0;
  for(; i<N && !found; i++)
    {
    GPE = GPNames[i]; 				// Try this name first
    if(vocal)
      cout << "\n\t Seeking Gnuplot Executable "
           << GPE;
    gexec = fopen(GPE.c_str(), "rb");		// Does this file exist?
    if(gexec != NULL) found=true;			// If so, use this command
    if(vocal)
      if(found) cout << " - Success!"; 
      else      cout << " - Not Found"; 
    }
  if(found) { fclose(gexec); return GPE; } 
  for(i=0; i<N && !found; i++)
    {
    GPE = GPNames[i] + string(".exe"); 				// Try this name first
    if(vocal)
      cout << "\n\t Seeking Gnuplot Executable "
           << GPE;
    gexec = fopen(GPE.c_str(), "rb");		// Does this file exist?
    if(gexec != NULL) found=true;			// If so, use this command
   if(vocal)
     if(found) cout << " - Success!"; 
     else      cout << " - Not Found"; 
    }
  if(found) { fclose(gexec); return GPE; } 
  for(i=0; i<N && !found; i++)
    {
    GPE = GPNames[i] + string(".out"); 				// Try this name first
    if(vocal)
      cout << "\n\t Seeking Gnuplot Executable "
           << GPE;
    gexec = fopen(GPE.c_str(), "rb");		// Does this file exist?
    if(gexec != NULL) found=true;			// If so, use this command
    if(vocal)
      if(found) cout << " - Success!"; 
      else      cout << " - Not Found"; 
    }
  if(found) { fclose(gexec); return GPE; } 
  return string("");
  }
*/

string GPExec(int warn)
  {
  string GPE = GPFind();			// Look for Gnuplot exec.
  if(!GPE.length())				// See if command found
    {						// 
    if(warn)                                    // If not, issue warnings as
      {						// desired
      GPerror(10, 1);                           //   Can't find gnuplot executable
      if(warn <= 1) GPerror(10, GPE, 1);        //   Setting executable to GPE
      else          GPfatality(11);             //   Stopping program, no executable
      }
    vector<string> GPNames = GPSeekStrings();	// Gnuplot path names
    GPE = GPNames[0];                           // Set GPE to default value
    }
  return GPE;
  }


void SetLineType(std::ostream& gnuload, int ltype)
  {
  string sds("set data style ");
  gnuload << sds;
  switch(ltype)
    {
    case -1: gnuload << string("points\n");      break;
    case 0:  gnuload << string("linespoints\n"); break;
    case 1:
    default: gnuload << string("lines\n");       break;
    case 2:  gnuload << string("boxes\n");       break;
    }
  }


void CloseMacro(std::ostream& gnuload)
  {
  gnuload << "pause -1 \'<Return>\' To Exit \n";// Command to pause before quit
  gnuload << "exit\n";				// Command to exit gnuplot
  gnuload.flush();
//  gnuload.close();				// Close gnuplot macro file
  }


void RunGnuplot(const string& gnumacro)
  {
  string pltcmd = GPExec() + string(" \"")	// Command to execute macro
                + gnumacro + string("\"\n");
  cout << "\n" << std::endl;				// Add lines/flush: nice screen
  if(system(pltcmd.c_str()))			// Plot to screen now
    {
    GPerror(13, 1);				// Cannot run gnuplot?
    string tmp = GPExec() + string(" \"")	// Command to execute macro
                + gnumacro + string("\"");
    GPerror(12, tmp, 1);			// System call failure
    GPfatality(14);				// Exit, need system "gnuplot"
    }
  }

//_____________________________________________________________________________
// B                      GAMMA Gnuplot 1D Output Files
//_____________________________________________________________________________

// ____________________________________________________________________________ 
//                             Gnuplot Parametric Plots
// ____________________________________________________________________________ 

/* These functions create Gnuplot compatible ASCII files in parametric form
   from input row_vectors. They are parametric in that the resulting plot will
   track points p[i] = (Re{vx[i]}, Im{vx[i]}) and these may be plotted from
   point to point.

           Input        filename  : Output filename    (single file output)
                   OR   ofstream  : Output data stream (multi-plot output)
                        vx        : Data vector
           Return                 : Void, file or ofstream is modified       */

 
void GP_xy(const string& filename, const row_vector &vx)
  {
  std::ofstream ofstr(filename.c_str());		// Open file for plotting
  GP_xy(ofstr, vx);				// Use GP_xy overload
  ofstr.close();				// Close the file stream
  }
 
void GP_xy(std::ofstream& ofstr, const row_vector& vx)
  {
  for(int i=0; i<vx.size(); i++)		// Loop over all the points
    ofstr << vx.getRe(i) << "  "                // Write them to output stream
          << vx.getIm(i) << "\n";               // (in ASCII format)
  }

 
// ____________________________________________________________________________ 
//                              Gnuplot Stack Plots 
// ____________________________________________________________________________ 

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

void GP_stack(std::ofstream& ofstr, const row_vector &vx, double row,
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
  double pt, lpt=HUGE_VAL;			// Set up for data reduction
  double ymax=HUGE_VAL,ymin=HUGE_VAL;		// (done if ordinates output)
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


// sosi - So far these matrix argument cannot be constant as
//        the matrix get_block is not constant.  When it is
//         this call should be changed too.


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
  std::ofstream ofstr(filename.c_str());		// Begins an output file stream
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

void GP_stack(std::ofstream& ofstr, matrix& mx, int ri,
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
 

void GnuplotStack(const string& filename, const vector<row_vector>& vxs,
                    int ri, double ymin, double ymax, double xmin, double xmax)
  {
  int nc = vxs[0].size();			// Number of columns
  int nr = vxs.size();				// Number of rows
  string Afile = filename + ".asc";
  std::ofstream ofstr(Afile.c_str());		// Begins an output file stream
  int i, j;
  double zmin=HUGE_VAL, zmax=-HUGE_VAL;
  double pt;
  for(i=0; i<nr; i++) 				// Determine a suitable scaling
    for(j=0; j<nc; j++) 			// for data reduction by first 
      { 					// finding global maxima
      if(ri >= 0)				// Search the real points
        {
        pt = vxs[i].getRe(j);
        if(zmax < pt) zmax = pt;
        if(zmin > pt) zmin = pt;
        }
      if(ri)					// Search the imaginary points
        {
        pt = vxs[i].getIm(j);
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
  for(i=0; i<nr; i++)				// Loop through plotted rows
    {
    if(ptaxis)					// If non-default y-axis values
      {
      y = ymin + double(i)*dely;		// Get the y-axis value
      GP_stack(ofstr,vxs[i],y,ri,xmin,xmax,cutoff);	// Use GP_stack overload
      }
    else					// If default y-axis values
      GP_stack(ofstr,vxs[i],	 		// Use GP_stack overload, type
               (double)i,ri,xmin,xmax,cutoff);  // i so right fcn is used
    }
  ofstr.close();

  int join=1;
  string Gfile = filename + ".gnu";
  std::ofstream gnuload(Gfile.c_str());           // File of gnuplot commands
  SetLineType(gnuload, join);                   // Set plot line type
  gnuload << "set parametric\n";                // Plot x vs y as parametric
  gnuload << "splot \"" << Afile << "\"\n";    // Command for gnuplot plot
  CloseMacro(gnuload);                          // Add pause+exit, close macro
  RunGnuplot(Gfile);                         // Run gnuplot with gnumacro
  }

// ____________________________________________________________________________ 
//                             Gnuplot Contour Plots 
// ____________________________________________________________________________ 

/* These routines are intended to output data files in ASCII format which are
   readable by the public domain program Gnuplot.  Their purpose is to provide
   a simple means to view GAMMA generated contour plots. Note that Gnuplot
   takes ASCII files filled with {x,y,z} points as input for 3D plots.  Only 
   one point per ASCII file line is allowed, the three ordinates separated by 	
   spaces. The x & y ordinates are optional but are output in these routines.
   For contour plots, the output ASCII file must have a blank line following 	
   each row.  Also, the points must be "gridded", meaning that all points of
   the 2D data set must be present (no skipping allowed as in the stack plots).	
   The plots themselves will be displayed as		
  									
  			[-x, +x] from left to right			
  			[-y, +y] from front to back			
  			[-z, +z] from top to bottom			
  									
   Useful gnuplot commands for stack plots are		
  									
  			set parametric       (independent rows)		
  			set data style lines (connected points)		
  			splot "filename"     (draw stack plot)               */
  									

void GP_contour(std::ofstream& ofstr, const row_vector &vx, double row,
                              int ri, double xmin, double xmax, double cutoff)

	// Input	ofstr     : Output filestream
	// 		vx        : Data vector
	//		row       : Row value or y-axis value
        //		ri	  : Flag for real versus imaginary plot
	//			    ri = 0  : reals only
	//			    ri = <0 : imaginaries only
	//			    ri = >0 : both 
        //              xmin      : Plot horizontal (x=0) point value
        //              xmax      : Plot horizontal (x=end) point value
	//		cutoff    : Data cutoff (where to set zero)
	// Return		  : Void, the information in vector vx is
	//			    output to ofstr in ASCII format suitable
	//			    for use in gnuplot.  Each data point is
	//			    written as {x,y,z} where x is the point
	//			    number, y the row number, and z the vector
	//			    point value 
	// Note			  : This is essentially the same function as
	//			    GP_stack except that all points are output

//    (xmin,row,vx(0))  (xmin+delx,row,vx(1))  ......... (xmax,row,vx(np-1))

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
  double pt;					// Set up for data reduction
  double ymax=HUGE_VAL,ymin=HUGE_VAL;			// (done if ordinates output)
  if(ptaxis && !cutoff)				// If ordinates output, find max & min
    {						// if no cutoff is set for compression
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
  if(!cutoff) cutoff = 1.e-3*(ymax-ymin);	// Cutoff for data reduction (resolution)
  if(ri >= 0)                                   // Write data for real plot
    for(i=0; i<np; i++)				//	Loop over all the points
      {
      if(ptaxis)				//	This if an x-axis is important
        {
        pt = Re(vx.get(i));			//	This is the point (intensity)
        ofstr<<xmin+double(i)*delx << "  "	//	First output the x axis value
              << row << " " 			//	Now output the y axis value
              << pt << "\n";			//	Now output the intensity value
        }
      else					// 	This if just y axis used (default)
        {					//	Must plot all points this way
        ofstr << Gdec(i) << " " 			//	because x values assumed evenly 
              << row << " " 			//	incremented
              << vx.getRe(i) << "\n";
        }
      }
  if(ri > 0)                                    // Blank line to separate plots
      ofstr << "\n";				// if both real and imaginaries output
  if(ri)                                        // Write data for imaginary plot
    for(i=0; i<np; i++)
      {
      if(ptaxis)
        {
        pt = vx.getIm(i);
        ofstr << xmin + double(i)*delx << "  "
              << row << " " 
              << pt << "\n";
        }
      else
        {
        ofstr << Gdec(i) << " " 
              << row << " " 
              << vx.getIm(i) << "\n";
        }
      }
  ofstr << "\n";			// Another empty line, assume
  }					// more rows will follow


// sosi - So far these matrix arguments cannot be constant as
//        the matrix get_block is not constant.  When it is
//         this call should be changed too.


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
	//			    where x is the point number, y the row number,
	//			    and z the vector point value 

void GP_contour(const string& filename, matrix& mx, int ri,
                 double ymin, double ymax, double xmin, double xmax)
  {
  std::ofstream ofstr(filename.c_str());		// Begins an output file stream
  GP_contour(ofstr,mx,ri,ymin,ymax,xmin,xmax);	// Use GP_contour overload
  ofstr.close();				// Close the file stream
  return;
  }
  

void GP_contour(std::ofstream& ofstr, matrix& mx, int ri,
                 double ymin, double ymax, double xmin, double xmax)

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


/*		    Here's How The Data Will Be Written
               (No Parenthesis, Line Feeds Between Points)

(xmn,ymn,   <0|mx|0>)  (xmn+delx,ymn,<0|mx|1>)   .....   (xmx,ymn,<0|mx|np-1>)
(xmn,ymn+dy,<1|mx|0>) (xmn+delx,ymn+dy,<1|mx|1>) ..... (xmx,ymin+dy,<1|mx|np-1>)
         |                       |                                |
         |                       |                                |
(xmn,ymx,<np-1|mx|0>) (xmn+delx,ymx,<np-1|mx|1>) .....  (xmx,ymx,<np-1|mx|np-1>)
 
  		      Here's How The Data Will Be Plotted

        Matrix                     Plot                  Plot w/ Axes
        ______                     ____                  ____________
  
 <0|0> --------- <0|N>     <N|0> --------- <N|N>   xmin,ymax ----- xmax,ymax
   |               |         |               |         |               |
   |               |         |               |         |               |
   |               |         |               |         |               |
   |               |         |               |         |               |
   |               |         |               |         |               |
 <N|0> --------- <N|N>     <0|0> --------- <0|N>   xmin,ymin ----- xmax,ymin */

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
  double cutoff = 1.e-3*(zmax-zmin);		// Cutoff for data reduction (resolution)
  row_vector vx;
  double y = 0;					// Set up y-axis scaling
  double dely = 0;				// if any has been specified
  int ptaxis = 0;				// Assume none has been set
  if(ymin != ymax)				// If it has then flag that it is
    {						// and figure out what the y-axis
    ptaxis = 1;					// increment is
    dely = (ymax-ymin)/double(nr-1);
    }
  for(i=0; i<nr; i++)				// Loop through the matrix rows
    {
    vx = mx.get_block(i,0,1,nc);		// Put this row as a row_vector		
    if(ptaxis)					// If non-default y-axis values
      {
      y = ymin + double(i)*dely;		// Get the y-axis value
      GP_contour(ofstr,vx,y,ri,xmin,xmax,cutoff);// Use GP_contour overload
      }
    else					// If default y-axis values
      GP_contour(ofstr,vx,	 		// Use GP_contour overload, type
               (double)i,ri,xmin,xmax,cutoff);  // i so right fcn is used
    }
  }
 
 
// ____________________________________________________________________________ 
//                            Gnuplot Spherical Plots
// ____________________________________________________________________________
 
 
/*
void GP_sphere(const string& name, const coord_vec& data, bool old)
  {
  int N = data.size();				// Length of coord vector
  row_vector proj(N);				// For projected data
  double TH, PH;
  double fact1 = 1.0;				// In radians
  double fact2 = PI/2.0;			// In radians
  if(old) { fact1 = 180.0/PI; fact2 = 90.0; }	// Older gnuplot used degrees
  coord pt;
  for(int l=0; l<N; l++)
    {
    pt = data.get(l);				// Current point
    TH = fact1*pt.theta();			// Point theta (deg)
    PH = fact1*pt.phi();			// Point phi (deg)
    proj.put(complex(PH,fact2-TH),l);		// Store point
    }
  GP_xy(name, proj);				// Create 2D ASCII file
  }


// sosi - probably good to have the projection angles as input here?
void GnuplotSphere(const string& name, const coord_vec& data)
  {
  string Afile = name + ".asc";			// Set data file name
  string Gfile = name + ".gnu";			// Set gnuplot load filename
  GP_sphere(Afile, data);			// Make data file (ASCII)
  GP_sphereplot(Gfile, Afile);			// Make gnuplot load file
  }						// and call gnuplot of plot
*/

//____________________________________________________________________________
//                        Gnuplot Array Visualization
//____________________________________________________________________________

void GP_array(const string& basename, matrix& mx,
                                           int type, bool plot, double cutoff)
  {
  vector<string> Afiles;		// Array of output ASCII files
  string fname, fnamea, fnameb;		// String for said files
  int i,j;
  switch(type)
    {
    case 0:                                     // Here for non-zeros shown
    default:					// This is also the default mode
      {
      fname = basename + string(".asc");	//   Output ASCII file name
      Afiles.push_back(fname);			//   Store the name in vector
      std::ofstream ofstr(fname.c_str());		//   Create ASCII output file
      for(i=0; i<mx.rows(); i++)		//   Loop array elements and
        for(j=0; j<mx.cols(); j++)		//   output coordinates of all
          if(norm(mx.get(i,j)) > cutoff)	//   non-zero matrix elements
            ofstr << j << "  " << i << "\n";
      ofstr.close();				//   Close the output ASCII file
      }
      break;
   case 1:					// Here for +/- split into two
     {						// ASCII files for display
     fnamea = basename + string("+.asc");	//   Output ASCII >0 elements
     fnameb = basename + string("-.asc");	//   Output ASCII <0 elements
     Afiles.push_back(fnamea);			//   Store 1st ASCII file name
     Afiles.push_back(fnameb);			//   Store 2nd ASCII file name
     std::ofstream ofstra(fnamea.c_str());		//   Open file for >0 elements
     std::ofstream ofstrb(fnameb.c_str());		//   Open file for <0 elements
     for(i=0; i<mx.rows(); i++)			//   Loop array elements and
       for(j=0; j<mx.cols(); j++)		//   output coordinates of the
         {					//   >0 elements to 1st file 
         if(mx.getRe(i,j) > cutoff)		//   while <0 element coords
           ofstra << j << "  " << i << "\n";	//   go into the 2nd file
         else if(mx.getRe(i,j) < -cutoff)
           ofstrb << j << "  " << i << "\n";
         }
     ofstra.close();				//   Close the 1st file
     ofstrb.close();				//   Close the 2nd file
     }
     break;
   }

  fname = basename + ".gnu";			// Name of gnuplot load file
  std::ofstream gnuplot(fname.c_str());              // File of gnuplot commands
  gnuplot << "set parametric\n";                // Plots pts parametrically
  gnuplot << "set xrange [0:" << mx.rows()-1    // Reverse x-axis
          << "] reverse\n";
  gnuplot << "set yrange [0:" << mx.cols()-1    // Reverse x-axis
          << "] reverse\n";
  gnuplot << "set nokey\n";                     // Don't put in any legend
  gnuplot << "plot ";
  for(unsigned k=0; k<Afiles.size(); k++)
    {
    gnuplot << "\'" << Afiles[k] << "\'";
    if(k+1 < Afiles.size()) gnuplot << ", ";
    }
  gnuplot << "\n";
  gnuplot << "pause -1 \'<Return> To Exit'\n";  // Command to pause
  gnuplot << "exit\n";                          // Exit gnuplot
  gnuplot.close();                              // Close the gnuplot load file
  string pltcmd = GPExec() + string(" \""); 			// Command to execute macro
  pltcmd += fname + string("\"\n");		// with file name fname
  system(pltcmd.c_str());                       // Plot to screen now
  }

// ____________________________________________________________________________ 
//                          Gnuplot Interactive Plots 
// ____________________________________________________________________________ 

//--------------------------- Common Parametric Plots -------------------------


	// Input	gnumacro  : A filename for gnu macro file
	//		file1D    : A filename of 1D gnuplot data
	//		join      : Flags whether to join points
	// Return	void	  : Information concerning how
	//			    to plot in gnuplot is output

void GP_xyplot(const string& gnumacro, const string& file1D, int join)
  {
  std::ofstream gnuload(gnumacro.c_str());		// File of gnuplot commands
  SetLineType(gnuload, join);			// Set plot line type
  gnuload << "set parametric\n";		// Plot x vs y as parametric
  gnuload << "plot \"" << file1D << "\"\n";	// Command for gnuplot plot
  CloseMacro(gnuload);				// Add pause+exit, close macro
  RunGnuplot(gnumacro);				// Run gnuplot with gnumacro
  }


//------------------------------ Common Stack Plots ---------------------------


void GP_stackplot(const string& gnumacro, const string& file2D)

	// Input	gnumacro  : A filename for gnu macro file
	//		file2D    : A filename of 2D gnuplot data
	// Return	void	  : Information concerning how
	//			    to plot in gnuplot is output

  {
  int join=1;
  std::ofstream gnuload(gnumacro.c_str());		// File of gnuplot commands
  SetLineType(gnuload, join);			// Set plot line type
  gnuload << "set parametric\n";		// Plot x vs y as parametric
  gnuload << "splot \"" << file2D << "\"\n";	// Command for gnuplot plot
  CloseMacro(gnuload);				// Add pause+exit, close macro
  RunGnuplot(gnumacro);				// Run gnuplot with gnumacro
  }


//--------------------------- Common Contour Plots ---------------------------

void GP_contplot(const string& gnumacro, const string& file2D, int pf)

	// Input	gnumacro  : A filename for gnu macro file
	//		file2D    : A filename of 2D gnuplot data
	//		pf	  : Plot flag (default is 0)
	//			 	0 - Only make contour plot
	//			       !0 - Plot contours & surface
	// Return	void	  : Information concerning how
	//			    to plot in gnuplot is output

  {
  int join=1;
  std::ofstream gnuload(gnumacro.c_str());		// File of gnuplot commands
  SetLineType(gnuload, join);			// Set plot line type
  gnuload << "set parametric\n";		// Plot x vs y as parametric
  gnuload << "set contour\n";			// Contours in plane
  if(!pf)
    {
    gnuload << "set nosurface\n";		// Don't include the surface
    gnuload << "set view 0,0,1\n";		// Look down into xy-plane 
    }
    
  gnuload << "splot \"" << file2D << "\"\n";	// Command for gnuplot plot
  CloseMacro(gnuload);				// Add pause+exit, close macro
  RunGnuplot(gnumacro);				// Run gnuplot with gnumacro
  }


//-------------------------- Common Spherical Plots --------------------------
   
 
/*
void GP_sphereplot(const string& gnumacro, const string& Aname, bool old)
 
        // Input        gnumacro  : A filename for gnu macro file
        //              file2D    : A filename of 2D gnuplot data
        //              pf        : Plot flag (default is 0)
        //                              0 - Only make contour plot
        //                             !0 - Plot contours & surface
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output
 
  {
  std::ofstream gnuload(gnumacro.c_str());		// File of gnuplot commands
  SetLineType(gnuload, 1);			// Set plot line type
  gnuload << "set parametric\n";                // Plots pts parametrically
  if(old)
    gnuload << "set angles degrees\n";		// Set for spherical plot
  gnuload << "set title \"3D Trajectory\"\n";   // Use this title
  gnuload << "set nokey\n";                     // Don't put in any legend
  gnuload << "set view 80,140,.6, 2.5\n";       // Set the view perspective
  gnuload << "set mapping spherical\n";         // Map cartesian->speherical
  gnuload << "set samples 32\n";                // This may long. lat. lines
  gnuload << "set isosamples 9\n";              // This many lat lines
  gnuload << "set urange [-pi/2:pi/2]\n";       // Plotting range
  gnuload << "set vrange [0:2*pi]\n";
  gnuload << "splot cos(u)*cos(v),cos(u)*sin(v)"// Now make 3D plot
          << ",sin(u) with lines 3 4,\\\n"
          << "\'" + Aname + "\' with lines 1 2\n";
  CloseMacro(gnuload);				// Add pause+exit, close macro
  RunGnuplot(gnumacro);				// Run gnuplot with gnumacro
  }

 
void GP_sphereplot(const string& gnumacro, int N, string* files)
 
        // Input        gnumacro  : A filename for gnu macro file
        //              N         : Number of plot (ASCII) files
        //              files     : Names of plot (ASCII) files
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output
 
  {
  std::ofstream gnuload(gnumacro.c_str());		// File of gnuplot commands
  SetLineType(gnuload, 1);			// Set plot line type
  gnuload << "set parametric\n";                // Plots pts parametrically
  gnuload << "set angles degrees\n";            // Set for spherical plot
  gnuload << "set title \"3D Trajectory\"\n";   // Use this title
  gnuload << "set nokey\n";                     // Don't put in any legend
  gnuload << "set view 80,140,.6, 2.5\n";       // Set the view perspective
  gnuload << "set mapping spherical\n";         // Map cartesian->speherical
  gnuload << "set samples 32\n";                // This may long. lat. lines
  gnuload << "set isosamples 9\n";              // This many lat lines
  gnuload << "set urange [-pi/2:pi/2]\n";       // Plotting range
  gnuload << "set vrange [0:2*pi]\n";
  gnuload << "splot cos(u)*cos(v),cos(u)*sin(v)"// Now make 3D plot
          << ",sin(u) with lines 3 4";
  gnuload << ",\\\n\'" << files[0] 
          << "\' with lines 1 2";
  for(int i=1; i<N; i++)
    gnuload << ",\\\n\'" << files[i]
            << "\' with lines " << 3+2*i
            << " " << 4+2*i;
  gnuload << "\n";
  CloseMacro(gnuload);				// Add pause+exit, close macro
  RunGnuplot(gnumacro);				// Run gnuplot with gnumacro
  cout << "\n\n";
  }
*/

// ____________________________________________________________________________ 
//                           Gnuplot Structure Functions
// ____________________________________________________________________________ 
  

void defaults(GPdat& G)

	// Input	G         : Plotting controls
	// Return	void	  : The plotting controls are
	//			    set for default values

  {
  string none = "";
  G.term = none;			// Terminal is default
  G.basedir = none;			// Base directory
  G.basename = none;			// Base filename
  G.cmd = none;				// System command
  G.ylabel = none;			// Label for yaxis
  G.xlabel = none;			// Label for xaxis
  G.title = none;			// Label for plot

  G.nokey = 0;				// Flag for nokey (yes/no)
  G.border = 1;				// Border flag (yes/no)
  G.join = 1;				// Line draw flag (yes/no)
  G.xtics = 1;				// X-axis tic flag (yes/no)
  G.ytics = 1;				// Y-axis tic flag (yes/no)
  G.xaxis = 1;				// Flag for x-axis draw (xes/no)
  G.yaxis = 1;				// Flag for y-axis draw (yes/no)
 
  G.xsize = 0;				// X-axis scaling factor
  G.ysize = 0;				// Y-axis scaling factor
  }

// ____________________________________________________________________________ 
//                           Gnuplot Output Blurbs
// ____________________________________________________________________________ 
  
void GP_stackblurb(std::ofstream& ostr, const string& plotname)

	// Input	ostr      : An output stream
	// 		plotname  : A gnuplot stack plotname
	// Return	void	  : Information concerning how
	//			    to plot in gnuplot is output
	//			    into the ofstream

  {
  ostr << "\n\n\t\tGNUPlot 2D ASCII Data File "
       << plotname << " Complete";
  ostr << "\n\t\tCommands To Plot: 1.) set data style lines"
       << "\n\t\t                  2.) set parametric"
       << "\n\t\t                  3.) splot \"" << plotname << "\"";
  ostr << "\n";
  }


void GP_contblurb(std::ofstream& ostr, const string& plotname)

	// Input	ostr      : An output stream
	// 		plotname  : A gnuplot contour plotname
	// Return	void	  : Information concerning how
	//			    to plot in gnuplot is output
	//			    into the ofstream

  {
  ostr << "\n\n\t\tGNUPlot 2D ASCII Data File "
       << plotname << " Complete";
  ostr << "\n\t\tCommands To Plot: 1.) set data style lines"
       << "\n\t\t                  2.) set parametric"
       << "\n\t\t                  2.) set contour both"
       << "\n\t\t                  3.) splot \"" << plotname << "\"";
  ostr << "\n";
  }

#endif 						// Gnuplot.cc
