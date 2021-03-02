/* FrameMaker.cc ************************************************-*-c++-*-
**									**
**	                        G A M M A 				**
**								 	**
**	FrameMaker 		                  Implementation   	**
**						 			**
**	Copyright (c) 1991, 1992		 			**
**	Tilo Levante, Scott Smith, Beat Meier				**
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
**  The FrameMaker Library Provides Functions to output different data	**
**  types to FrameMaker. 						**
**								 	**
**  The functions provided herein assign group ID numbers to the	**
**  graphic objects they create.  The following is a listing of the	**
**  currently used ID numbers, the function they are in, and what they	**
**  are used for.  Programmers can use this list to avoid mistakenly	**
**  using the same ID number for different objects.			**
**	 					 			**
**	              ID	Function	Use		 	**
**	              --	--------	---------------------	**
**			 					 	**
**	               8	FM_Axis		For x-axis		**
**                     9	FM_Axis		For y-axis	 	**
**	              11	FM_AFrame_Set	For default A_Frame	**
**	 							 	**
*************************************************************************/

#ifndef   GFrameMaker_cc_			// Is file already included?
#  define GFrameMaker_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#  endif

#include <GamIO/FrameMaker.h>			// Include the header
#include <GamIO/FrameMakerM.h>			// Include FM MIF functions
#include <GamIO/FrameMakerP.h>			// Include FM parameters
#include <Basics/Gconstants.h>			// Include PI
#include <Basics/Gutils.h>			// Include GAMMA errors
#include <stdlib.h>
#include <string>				// Inlcude stdlibc++ strings
#include <iostream>				// Include ostreams (cout)
#include <fstream>				// Include filestreams
#include <cstdio>				// Include FILE
#include <cmath>				// Inlcude for HUGE_VAL

using namespace std;

#ifndef _MSC_VER				// If we are not using MSVC++
  std::string LBrak("\[");			// set string for bracket here
#else						// else we must define one
  std::string LBrak("[");			// like this
 #pragma warning (disable : 4786)		// Kill STL namelength warnings
#endif
	
#include <Basics/StringCut.h>			// Include form & dec functions

// ____________________________________________________________________________
// i                        FRAMEMAKER ERROR HANDLING
// ____________________________________________________________________________


/*         Input                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */  

void FM_error(int error, int noret=1)
  {
  std::string hdr("GAMMA FrameMaker");
  switch(error)
    {
    case 1: GAMMAerror(hdr,"Stack Plot |yinc| too Large! Adjusting", noret);		// (1)
      break;
    case 2: GAMMAerror(hdr,"Stack Plot |xinc| too Large! Adjusting", noret);		// (2)
      break;
    case 3: GAMMAerror(hdr,"Difficult Stack Plot Vertical Scaling - Iterating", noret);	// (3)
      break;
    case 4: GAMMAerror(hdr,"1-D Plot Data is Constant, Horizontal Plot", noret);	// (4)
      break;
    case 5: GAMMAerror(hdr,"Data Vector Found to Contain 0 or 1 point", noret);		// (5)
      break;
    case 7: GAMMAerror(hdr,"Data Matrix Found Constant. Cannot Contour!", noret);	// (7)
      break;
    case 8: GAMMAerror(hdr,"Threshold Above Data. Adjusting to Global Min.!", noret);	// (8)
      break;
    case 9: GAMMAerror(hdr,"Threshold Below Data. Adjusting to Global Max.!", noret);	// (9)
      break;
    case 11: GAMMAerror(hdr,"Data Vector x Span is Zero. Histogram Impossible.", noret);// (11)
      break;
    case 12: GAMMAerror(hdr,"Too Few Columns in Data Matrix for Contouring.", noret);	// (12)
      break;
    case 13: GAMMAerror(hdr,"Too Few Rows in Data Matrix for Contouring.", noret);	// (13)
      break;
    case 15: GAMMAerror(hdr,"Initial Contour Threshold Bad, Adjusting.", noret);	// (15)
      break;
    case 16: GAMMAerror(hdr,"Points in Contour Exceeds Array Dimension!", noret);	// (16)
      break;
    case 17: GAMMAerror(hdr,"Fix Contouring Routine Array Dimensions!", noret);		// (17)
      break;
    case 50: GAMMAerror(hdr,"Plotted Array Contains < 2 Rows!", noret);			// (50)
      break;
    case 51: GAMMAerror(hdr,"Cannot Produce A Stack Plot.", noret);			// (51)
      break;
    case 100: GAMMAerror(hdr,"Cannot Proceed with 1D Plot", noret);			// (100)
      break;
    default: GAMMAerror(hdr,"Unkown Error", noret);
    }
  }

volatile void FM_fatality(int error)
  {
  FM_error(error);
  if(error) FM_error(0);
  GAMMAfatal();					// Clean exit from program
  }

// ______________________________________________________________________
//                 FRAMEMAKER MIF AUXILIARY GRAPIC FUNCTIONS
// ______________________________________________________________________

     	// Input	out	  : std::ostream for FrameMaker File
	//		ID	  : Tic mark FrameMaker ID
     	//         	direction : 'x' or 'y' flags which axis
     	//         	offset    : Tic mark position
	//		pmin      : Minimum value of axis
	//		pmax      : Maximum value of axis
	//		min       : Minimum for axis labeling
	//		max       : Maximum of axis labeling
     	// Output		  : None, void function. Draws an axis tics
	//			    to the output file in FrameMaker Format.
	// Note			  : Assumes file are currently in a Frame;
     	// Note			  : For x-axis tics, offset = y position &
     	//			    for y-axis tics, offset = x position
     	// Note			  : FM_Axis sets x-axis tics ID=8,
	//			    y-axis tics ID=9

void FM_Axis_tics (std::ostream &out, int ID, char direction, double offset,
		 	double pmin, double pmax, double min, double max)
  {
  double sf = 1; 
  if(max)
    sf = 1;
  if(fabs(max-min) < 1.e-9*sf)		// Warning if Plot is Horizontal
    FM_error(4);
  double delta = exp(log(10.0)*		// Estimate the distance between two tics
       (floor(log(max-min)/log(10.0))-2));
  double z[3];				// Search for a delta to yeild
  z[0] = delta;				// approx. 6 tick marks. Begin with
  z[1] = 2.0*delta;			// three guesses and keep increasing until
  z[2] = 5.0*delta;			// one is suitable
  while(((max-min)/delta)>7 )
    {
    int i=0;
    while ((((max-min)/delta)>7)&&(i<3))
      {
      delta = z[i];
      z[i]=delta*10;
      i++;
      }
    }
						// Begin tic marks
  double pdelta = delta*(pmax-pmin)/(max-min);	// Plot tic increment 
  double x0 = delta*ceil(min/delta);
  double p0 = (x0-min)*(pmax-pmin)/(max-min)	// Initiial tic mark
					 + pmin;
  for (; x0<max; x0+=delta, p0+=pdelta)
    {						// Put in tics & labels
    if (direction == 'x')
      {
      FM_Line(out, ID, 0, 0.5, p0,
                     offset, p0, offset-0.2);
      FM_TextLine(out,ID,0,p0,offset+0.5,x0);	// Text 0.5cm below x axis
      }
    else
      {
      FM_Line (out, ID, 0, 0.5,
                 offset, p0, offset+0.2, p0);
      FM_TextLine(out,ID,1,offset-0.1, p0, x0); // Text 0.1cm left y axis
      }
    }
  }



     	// Input	out	  : std::ostream for FrameMaker File
     	//         	direction : 'x' or 'y' flags which axis
     	//         	offset    : Axis position
	//		pmin      : Minimum value of axis
	//		pmax      : Maximum value of axis
	//		min       : Minimum for axis labeling
	//		max       : Maximum of axis labeling
     	// Output		  : None, void function. Draws an axis
	//			    to the output file in FrameMaker Format.
	// Note			  : Assumes file "out" currently in a Frame
     	// Note			  : For x-axis, offset = y position and
     	//			    for y-axis, offset = x position
     	// Note			  : Sets x-axis ID=8, y-axis ID=9

void FM_Axis(std::ostream &out, char direction, double offset,
		 	double pmin, double pmax, double min, double max)
  {
  int ID;
  if (direction=='x')				// Write a horizontal axis
    {
    ID = 8;					// Set the x-axis FM ID
    FM_Line(out,ID,0,1.0,			// Output the x-axis line
                    pmin,offset,pmax,offset);
    FM_Axis_tics(out, ID, direction, offset,	// Write in the tic marks
		 	pmin, pmax, min, max);
    }
  else						// Write a vertical axis
    {		
    ID = 9;					// Set the y-axis FM ID
    FM_Line(out,ID,0,1.0,			// Output the y-axis line
                     offset,pmin,offset,pmax);
    FM_Axis_tics(out, ID, direction, offset,	// Write in the tic marks
		 	pmin, pmax, min, max);
    }
  FM_Group (out, ID);				// Group line & tics together
  return;
  }


void FM_maxima(const row_vector &vx, double &min, double &max, int ri)

	// Input	vx       : Data vector
	//		min      : Minimum of vx
	//		max      : Maximum of vx
        //		ri	 : Flag for real versus imaginary
	//			   0=reals, non-zero=imaginaries
	// Return		 : Void, min and max are altered

  {
  int np = vx.elements();
  max = -HUGE_VAL;			// Compute max & min of vx
  min = HUGE_VAL;
  double t;
  if(ri)
    {
    for (int i=0; i<np; i++)	// Search all imaginary points
      {				// for maximum & minimum
      t = vx.getIm(i);
      if (t<min) min = t;
      if (t>max) max = t;
      }
    }
  else
    {
    for (int j=0; j<np; j++)	// Search all real points
      {				// for maximum & minimum
      t = vx.getRe(j);
      if (t<min) min=t;
      if (t>max) max=t;
      }
    }
  }

void FM_borders(double xsize, double ysize,
	      double &top, double &bottom, double &right, double &left)

	// Input	xsize     : Plot horizontal (x) dimension in cm
	//		ysize     : Plot vertical (y) dimension in cm
	//		top       : Actual plot top in cm
	//		bottom    : Actual plot bottom in cm
        //		right	  : Actual plot right in cm
        //		left	  : Actual plot left in cm
	// Return       void	  : Values of top, bottom, right, & left
	//			    are set.

  {
  top = 0.5;		// Top & bottom borders (relative to Frame!)
  bottom = ysize-1;	// 0.5 & 1 cm margins at top & bottom
  right = xsize-0.5;	// Left & right borders (relative to Frame!)
  left = 2;		// 2 & 0.5 cm margins at left and right
  }


void FM_Box(std::ostream &out, double xmin, double ymin, double xmax, double ymax)

     	// Input	out	  : An std::ostream for FrameMaker file
	//		xmin      : Minimum x-value of box (cm)
	//		ymin      : Minimum y-value of box (cm)
	//		xmax      : Maximum x-value of box (cm)
	//		ymax      : Maximum y-value of box (cm)
     	// Output       void	  : Draws a box to the output
	//			    file in FrameMaker format.
	// Note			  : Assumes file "out" currently in a Frame
     	// Note			  : Sets box ID to 11

    {
    int ID = 11;					// Set the box FM ID
    FM_Line(out,ID,0,1.0, xmin,ymin,xmax,ymin); 	// Output lower horizontal
    FM_Line(out,ID,0,1.0, xmax,ymin,xmax,ymax); 	// Output right vertical
    FM_Line(out,ID,0,1.0, xmax,ymax,xmin,ymax); 	// Output upper horizontal
    FM_Line(out,ID,0,1.0, xmin,ymax,xmin,ymin); 	// Output left vertical
    FM_Group (out, ID);					// Group line & tics together
    return;
    }

// ____________________________________________________________________________
//                      FRAMEMAKER EXECUTION FUNCTIONS
// ____________________________________________________________________________

        // Input        none    : Output filename
	//		warn    : Warning flag
        // Return       FME     : A string for the FrameMaker
        //                        executable command (hopefully)
	// Note         	: Add more FrameMaker paths to FMNames in
	//			  FMSeekStrings as needed.

std::vector<std::string> FMSeekStrings()
  {
  std::vector<std::string> FMNames;
  std::string FM("FrameMaker");
  std::string WFM("wframe");
  std::string UB("/usr/bin/");
  std::string ULB("/usr/local/bin/");
  std::string UFW("/usr/freeware/bin/");
  std::string CW("/cygwin/bin/");
  std::string CWU("/cygwin/usr/bin");
  std::string CWUL("/cygwin/usr/local/bin");
  std::string MW("/mingw");
  std::string MWB("/mingw/bin");
  std::string AD("/Adobe");
  std::string PF("/Program Files");
  std::string FM55("/FrameMaker5.5");
  std::string FM6("/FrameMaker6");
  FMNames.push_back(FM);			// If command local
  FMNames.push_back(ULB + FM);			// GNUish install spot
  FMNames.push_back(UB + FM);			// Linux RPM install spot
  FMNames.push_back(UFW + FM);			// SGI tardist install spot
  FMNames.push_back(WFM);			// If command local Windoze
  FMNames.push_back(ULB + WFM);			// Gnuish for cygwin
  FMNames.push_back(UB + WFM);			// Why not for cygwin
  #if defined(_MSC_VER) || defined(__MINGW32__)
    std::string DC("C:");
    std::string DD("D:");
    FMNames.push_back(DC + PF  + AD + FM55 + "/" + FM + ".exe");
    FMNames.push_back(DC + CW   + FM);
    FMNames.push_back(DC + CWU  + FM);
    FMNames.push_back(DC + CWUL + FM);
    FMNames.push_back(DC + ULB  + WFM);
    FMNames.push_back(DC + CW   + WFM);
    FMNames.push_back(DC + CWU  + WFM);
    FMNames.push_back(DC + CWUL + WFM);
    FMNames.push_back(DD + ULB  + FM);
    FMNames.push_back(DD + CW   + FM);
    FMNames.push_back(DD + CWU  + FM);
    FMNames.push_back(DD + CWUL + FM);
    FMNames.push_back(DD + ULB  + WFM);
    FMNames.push_back(DD + CW   + WFM);
    FMNames.push_back(DD + CWU  + WFM);
    FMNames.push_back(DD + CWUL + WFM);
  #endif

  return FMNames;
  }

std::string FMFind(bool vocal)
  {
  std::vector<std::string> FMNames = FMSeekStrings();	// FrameMaker path names
  int N = FMNames.size();			// Number of path names
  std::string FME, FMName;				// FrameMaker execution command
  bool found = false;				// Flag if found executable
	ifstream gexec;				    // FrameMaker executable?
  int i=0;
 
	for(; i<N && !found; i++)
    {
    FME = FMNames[i]; 				// Try this name first
    if(vocal)
      std::cout << "\n\t Seeking FrameMaker Executable " << FME;
		gexec.open(FME.c_str(), ios_base::binary);
    if(gexec.good() == true) 
			found=true;			// If so, use this command
    if(vocal)
    { if(found) 
				std::cout << " - Success!"; 
      else      
				std::cout << " - Not Found"; 
    }
    }
  if(found) 
		{ 
		gexec.close(); 
		return FME; 
		}

  for(i=0; i<N && !found; i++)
    {
    FME = FMNames[i] + std::string(".exe"); 				// Try this name first
    if(vocal)
      std::cout << "\n\t Seeking FrameMaker Executable " << FME;
		gexec.open(FME.c_str(), ios_base::binary);
    if( gexec.good() == true )
			found=true;			// If so, use this command
	if(vocal)
	{  if(found) 
		std::cout << " - Success!"; 
	  else      
		std::cout << " - Not Found"; 
        }
    }
  if(found) 
		{ 
		gexec.close(); 
		return FME; 
		} 

  for(i=0; i<N && !found; i++)
    {
    FME = FMNames[i] + std::string(".out"); 				// Try this name first
    if(vocal)
      std::cout << "\n\t Seeking FrameMaker Executable " << FME;
		gexec.open(FME.c_str(), ios_base::binary);
		if( gexec.good() == true ) 
			found=true;			// If so, use this command
    if(vocal)
    { if(found) 
	std::cout << " - Success!"; 
      else      
	std::cout << " - Not Found"; 
    }
    }
  if(found) 
	{ 
		gexec.close(); 
		return FME; 
	} 

  return std::string("");
  }


std::string FMExec(int warn)
  {
  std::string FME = FMFind();			// Look for FrameMaker exec.
  if(!FME.length())				// See if command found
    {						// 
    if(warn)                                    // If not, issue warnings as
      {						// desired
FM_error(10);
//      FM_error(10, 1);                          //   Can't find gnuplot executable
//      if(warn <= 1) FM_error(10, FME, 1);       //   Setting executable to FME
//      else          FM_fatality(11);            //   Stopping program, no executable
      }
    std::vector<std::string> FMNames = FMSeekStrings();	// FrameMaker path names
    FME = FMNames[0];                           // Set FME to default value
    }
  return FME;
  }


// ____________________________________________________________________________
//                      FRAMEMAKER 1D PLOTTING FUNCTIONS
// ____________________________________________________________________________

/* These function take a single row_vector and plot parameters then produce
   a FrameMaker MIF file containing the 1D plot, Re(i)/Im(i) vs i 

	Input 		out	: Output file stream (ASCII Frame MIF file)
			vx	: Input data vector to be plotted
			xsize	: Plot horizontal dimension in cm
			ysize	: Plot vertical dimension in cm
	  		pmax	: Plot vertical maximum (axis label)
	  		pmin	: Plot vertical minimum (axis label)
          		ri	: Flag for real versus imaginary plot
	  			    0=reals, non-zero=imaginaries

  Function                                Purpose
  --------	---------------------------------------------------------------
  FM_1D_ri      This is a auxiliary function to FM_1D that uses a PolyLine size 
                defaulting to 128 points (PLmax = 128)
*/

void FM_1D_ri(std::ostream &out, const row_vector &vx, double xsize,
                               double ysize, double pmin, double pmax, int ri)
  {
  int np = vx.elements();		// Total points
  if(np <= 1)				// Must have more than 1 point
    {
    FM_error(5);			// Error from vector of 0 or 1 pt. 
    FM_fatality(100);			// Fatal error flagged in 1D plot
    }
  double ymin, ymax, ydel;		// Begin scaling
  FM_maxima(vx, ymin, ymax, ri);	// Get vertical maxima
  ydel = ymax-ymin;
  if(ydel == 0) ydel = ymax;
  if(ydel == 0) ydel = 1;
  double top, bottom, right, left;	// Declare and set margins
  FM_borders(xsize, ysize,
             top, bottom, right, left);
  double yscale =			// Scaling in x and y directions
             (bottom-top)/(ymax-ymin);
  double xscale = (right-left)/(np-1);

  int FrameID = 11;			// Set Frame and Polyline ID numbers
  int PolyLineID = 1;
  if(ri)
    {
    FrameID++;
    PolyLineID = 2;
    }
  FM_AFrame_Set(out, xsize, ysize ,FrameID);	// Write Frame header
						// Begin plotting data
  int PLmax = 128;				// Limit size of Polyline(s)
  row_vector PL(PLmax+1);			// Vector for Polyline 128+1 pt(max.)
  double xcm, ycm, yin;
  xcm = left;				// Initial point, 1st Polyline
  yin = ymax*yscale + top;
  if(ri)
    ycm = yin - vx.getIm(0)*yscale;
  else
    ycm = yin - vx.getRe(0)*yscale;
  int npl = 0;				// Polyline count
  int pts = 0;				// Polyline point count (-1)
  int skip = 0;
  PL.put(complex(xcm, ycm), pts);
  for(int j=1; j<np; j++)		// Go through each point, get
    {					// coordinates in cm
    xcm += xscale;			// x(j) = j*xscale+left
    if(ri)				// Treat Imaginaries
      ycm = yin - vx.getIm(j)*yscale;
    else				// Treat Reals
      ycm = yin - vx.getRe(j)*yscale;
    if(pts == 0)			// Add pt if 2nd in Polyline
      {
      pts++;				
      PL.put(complex(xcm, ycm), pts);
      skip = 0;
      }
    if(j == np-1)			// Add pt if last in Polyline
      {
      pts++;				// or last in data vector
      if(skip)
        {
        PL.put(complex(xcm-xscale, Im(PL(pts-1))),pts);
        pts++;
        }
      PL.put(complex(xcm, ycm), pts);
      skip = 0;
      }
    else if((fabs(Im(PL(pts)) - ycm))	// Add pt if height changes
                               >0.0001) // (sets vertical resolution!)
      {
      pts++;
      if(skip)
        {
        PL.put(complex(xcm-xscale, Im(PL(pts-1))), pts);
        pts++;
        }
      PL.put(complex(xcm, ycm), pts);
      skip = 0;
      }
    else
      skip = 1;

    if(pts >= PLmax-1)			// Limit Polyline length to PLmax
      {
      FM_PolyLine(out, PL, PolyLineID, 15, pts+1, 0);
      npl++;
      PL.put(PL(pts), 0);
      pts = 0;
      skip = 0;
      }
    }
  if(pts >= 1)				// Output Polyline if end
    {
    FM_PolyLine(out,PL,PolyLineID,15,pts+1,0);
    npl++;
    }
  if (npl>1)				// Group all Polylines
    FM_Group (out, PolyLineID);
  FM_Axis (out, 'x', bottom, left, right, pmin, pmax);
  FM_Axis (out, 'y', left, bottom, top, ymin, ymax);
  FM_AFrame_End (out);			// End of Frame
  }


void FM_1D_ri(std::ostream &out, const row_vector &vx, FMPar& FMP)

	// Input	out       : Output stream
	// 		vx        : Data vector
	//		FMP       : FrameMaker plot parameters
	// Return		  : Void, file out is modified
	// Note			  : AUXILIARY FUNCTION FOR FM_1D
 
  {
  int np = vx.elements();			// Total points
  if(np <= 1)					// Must have > 1 point
    {
    FM_error(5);				// Error: vector <2 pts. 
    FM_fatality(100);				// Fatal: 1D plot problem
    }

//	              Set Plot Data Scaling Values

  double ymin, ymax;				// For vertical scaling
  FMP.maxima(vx, ymin, ymax);			// Get vector vertical maxima
  double ydel = ymax - ymin;			// Total vector vertical span
  if(ydel == 0) ydel = ymax;
  if(ydel == 0) ydel = 1;
  double top, bottom, right, left;		// FrameMaker boundaries
  FMP.borders(top, bottom, right, left);	// Set the FM boundaries
  double yscale = (bottom-top)/(ymax-ymin); 	// Vertical scaling factor
  double xscale = (right-left)/(np-1);		// Horizontal scaling factor

//	      Start An Anchored Frame To Contain The 1D Plot

  double hsz = FMP.HPlotSize();			// Horizontal plot size (cm) 
  double vsz = FMP.VPlotSize();			// Vertical plot size (cm) 
  int du = FMP.DataUse();			// Type of data to use 
  int FrameID = 11 + du;			// Set An Anchored Frame ID
  int PolyLineID = 1+ du;			// Set a PolyLine ID
  FM_AFrame_Set(out,hsz,vsz,FrameID);	 	// Write Frame Header

//		           Set Scaling For 1D Plot

  int plsz = FMP.PLsize(); 			// Allowed PolyLine size
  row_vector PL(plsz + 1);			// Vector for PolyLines
  double xcm, ycm, yin;				// FM plot coordinates in cm
  xcm = left;					// Initial point, 1st PolyLine
  yin = ymax*yscale + top;			// Base vertical coordinate
  switch(du)
    {
    case 0:
    default:					// Plotting reals
      ycm = yin - vx.getRe(0)*yscale; break;	// set 1st PL point
    case 1:					// Plotting imaginaries 
      ycm = yin - vx.getIm(0)*yscale; break;	// set 1st PL point
    case 3:					// Plotting  norms
      ycm = yin - norm(vx.get(0))*yscale; break;// set 1st PL point
    }

//		           Begin Plotting 1D Data

  int npl = 0;					// Polyline count
  int pts = 0;					// Polyline point count (-1)
  int skip = 0;				// Flag for point skip
  double xm1, ym1;				// Temp values
  PL.put(complex(xcm, ycm), pts);		// Set 1st point, 1st PL
  for(int j=1; j<np; j++)			// Go through each point, get
    {						// coordinates in cm
    xcm += xscale;				// x(j) = j*xscale+left
    switch(du)
      {
      case 0:
      default:					// Plotting reals
        ycm = yin - vx.getRe(0)*yscale; break;	// set 1st PL point
      case 1:					// Plotting imaginaries 
        ycm = yin - vx.getIm(0)*yscale; break;	// set 1st PL point
      case 3:					// Plotting  norms
        ycm = yin - norm(vx.get(0))*yscale;	// set 1st PL point
      }
    if(pts == 0)				// We MUST add this point
      { 					// if 2nd point in PolyLine
      pts++;					//	Increment point count
      PL.put(complex(xcm, ycm), pts);		// 	Set the point
      skip = 0;				//	Flag not skipped
      }
    if(j == np-1)				// If point is last in data
      { 					// we MUST add the point too!
      pts++;					//	Increment point count
      if(skip)					// 	If last point was
        {					//	skipped must add it
        xm1 = xcm-xscale;			//	Previous x value
	ym1 = PL.getIm(pts-1);			//	Previous y value
        PL.put(complex(xm1, ym1), pts);		//	Add previous point
        pts++;					//	Go to next point
        }
      PL.put(complex(xcm, ycm), pts);		// 	Set the point
      skip = 0;					//	Flag its not skipped
      }
    else if((fabs(Im(PL(pts)) - ycm)) >0.0001) 	// If height changes MUST add
      { 					// (sets vertical resolution!)
      pts++;					//	Increment point count
      if(skip)					// 	If last point was
        {					//	skipped must add it
        xm1 = xcm-xscale;			//	Previous x value
	ym1 = PL.getIm(pts-1);			//	Previous y value
        PL.put(complex(xm1, ym1), pts);		//	Add previous point
        pts++;					//	Go to next point
        }
      PL.put(complex(xcm, ycm), pts);		// 	Set the point
      skip = 0;					//	Flag its not skipped
      }
    else skip = 1;				// Point line previous, SKIP

    if(pts >= plsz-1)				// Check PolyLine length
      { 					// Write it if length>=PLmax
      FM_PolyLine(out,PL,PolyLineID,15,pts+1,0);//	Output PolyLine
      npl++;					//	Increment PL count
      PL.put(PL(pts), 0);			//	Reset the 1st point
      pts = 0;					//	Reset the point count
      skip = 0;				//	Flag point not skipped
      }
    }
  if(pts >= 1)					// Output last Polyline
    { 						// more than 1 point long
    FM_PolyLine(out,PL,PolyLineID,15,pts+1,0);	// 	Output PolyLine
    npl++;					//	Increment PL count
    }

//		End The Anchored Frame Containing The 1D Plot

  double pmin = FMP.HAxMin();			// H. Axis left value
  double pmax = FMP.HAxMax();			// H. Axis right value
  if(npl > 1) FM_Group(out,PolyLineID);		// Group all Polylines
  FM_Axis(out,'x',bottom,left,right,pmin,pmax);	// Draw horizontal axis
  FM_Axis(out,'y',left,bottom,top,ymin,ymax);	// Draw vertical axis
  FM_AFrame_End(out);				// End of Frame
  }



	// Input	filename  : Output filename
	// 		vx        : Data vector
	//		xsize     : Plot horizontal (x) dimension in cm
	//		ysize     : Plot vertical (y) dimension in cm
	//		pmax      : Plot horizontal maximum (axis label)
	//		pmin      : Plot horizontal minimum (axis label)
        //		ri	  : Flag for real versus imaginary plot
	//			    ri = 0  : reals only
	//			    ri = <0 : imaginaries only
	//			    ri = >0 : both 
	// Return		  : Void, file out is modified

void FM_1D(const std::string& filename, const row_vector &vx, double xsize,
                                double ysize, double pmin, double pmax, int ri)
  {
  if((xsize < 3) || (xsize > 20)) xsize = 14;	// Insure plot sizes sensible
  if((ysize < 3) || (ysize > 27)) ysize = 14;
  std::ofstream out(filename.c_str());		// Open file for plotting
  FM_Begin(out);				// FrameMaker begin comment
  FM_AFrames_Begin(out);			// Anchored Frame begin
  if(ri >= 0)
    FM_1D_ri(out,vx,xsize,ysize,pmin,pmax,0);	// Plot reals
  if(ri)
    FM_1D_ri(out,vx,xsize,ysize,pmin,pmax,1);	// Plot imaginaries
  FM_AFrames_End(out);				// End of AFrame
  FM_TextFlow_Set(out);				// Set TextFlow
  if(ri >= 0)
    {
    FM_Paragraph_Set(out);			// Set real Paragraph
    out << "   <String `real    '> \n";
    FM_Paragraph_End(out, 11);			// End real Paragraph
    }
  if(ri)
    {
    FM_Paragraph_Set(out);			// Set imag Paragraph
    out << "   <String `imaginary '> \n";
    FM_Paragraph_End(out, 12);			// End imag Paragraph
    }
  FM_TextFlow_End(out);
  FM_End(out);					// FrameMaker end comment
  }


void FM_1D(const std::string& filename, const row_vector &vx, FMPar& FMP)

	// Input	filename  : Output filename
	// 		vx        : Data vector
	//		FMP       : FrameMaker plot parameters
	// Return		  : A FrameMaker MIF file is
	//			    made called filename to contain
	//			    1D plot(s) of data in vector vx
	//			    according to parameters in FMP

  {
  FMP.plotsize(); 				// Insure plot sizes sensible
  std::ofstream out(filename.c_str());		// Open file for plotting
  FM_Begin(out);				// FrameMaker begin comment
  FM_AFrames_Begin(out);			// Anchored Frame begin
  datause du = FMP.DataUse();			// Type of data to use
  switch(du)
    {
    case plotreals:
    default:
      FM_1D_ri(out, vx, FMP);			// Make the 1D plot
      FM_AFrames_End(out);			// End of Anchored Frame
      FM_TextFlow_Set(out);			// Set TextFlow
      FM_Paragraph_Set(out);			// Set Paragraph
      out << "   <String `Reals    '> \n";	// Write its reals
      FM_Paragraph_End(out, 11);		// End Paragraph
      FM_TextFlow_End(out);
      break;
    case plotimags:
      FM_1D_ri(out, vx, FMP);			// Make the 1D plot
      FM_AFrames_End(out);			// End of Anchored Frame
      FM_TextFlow_Set(out);			// Set TextFlow
      FM_Paragraph_Set(out);			// Set Paragraph
      out << "   <String `Imaginaries   '> \n";	// Write its imaginaires
      FM_Paragraph_End(out, 11);		// End Paragraph
      FM_TextFlow_End(out);
      break;
    case plotnorms:
      FM_1D_ri(out, vx, FMP);			// Make the 1D plot
      FM_AFrames_End(out);			// End of Anchored Frame
      FM_TextFlow_Set(out);			// Set TextFlow
      FM_Paragraph_Set(out);			// Set Paragraph
      out << "   <String `Norms         '> \n";	// Write its norms
      FM_Paragraph_End(out, 11);		// End Paragraph
      FM_TextFlow_End(out);
      break;
    case plotboth:
      FMP.DataUse(plotreals);
      FM_1D_ri(out, vx, FMP);			// Make the 1D plot
      FM_AFrames_End(out);			// End of Anchored Frame
      FM_TextFlow_Set(out);			// Set TextFlow
      FM_Paragraph_Set(out);			// Set Paragraph
      out << "   <String `Reals         '> \n";	// Write its reals
      FM_Paragraph_End(out, 11);		// End Paragraph
      FM_TextFlow_End(out);

      FMP.DataUse(plotimags);
      FM_1D_ri(out, vx, FMP);			// Make the 1D plot
      FM_AFrames_End(out);			// End of Anchored Frame
      FM_TextFlow_Set(out);			// Set TextFlow
      FM_Paragraph_Set(out);			// Set Paragraph
      out << "   <String `Imaginaries   '> \n";	// Write its imaginaries
      FM_Paragraph_End(out, 11);		// End Paragraph

      FMP.DataUse(du);				// Reset data use value
      break;
    case plotall3:
      FMP.DataUse(plotreals);
      FM_1D_ri(out, vx, FMP);			// Make the 1D plot
      FM_AFrames_End(out);			// End of Anchored Frame
      FM_TextFlow_Set(out);			// Set TextFlow
      FM_Paragraph_Set(out);			// Set Paragraph
      out << "   <String `Reals         '> \n";	// Write its reals
      FM_Paragraph_End(out, 11);		// End Paragraph
      FM_TextFlow_End(out);

      FMP.DataUse(plotimags);
      FM_1D_ri(out, vx, FMP);			// Make the 1D plot
      FM_AFrames_End(out);			// End of Anchored Frame
      FM_TextFlow_Set(out);			// Set TextFlow
      FM_1D_ri(out, vx, FMP);			// Make the 1D plot
      FM_AFrames_End(out);			// End of Anchored Frame
      FM_TextFlow_Set(out);			// Set TextFlow
      FM_Paragraph_Set(out);			// Set Paragraph
      out << "   <String `Imaginaries '> \n";	// Write its imaginaries
      FM_Paragraph_End(out, 11);		// End Paragraph
      FM_TextFlow_End(out);

      FMP.DataUse(plotnorms);
      FM_1D_ri(out, vx, FMP);			// Make the 1D plot
      FM_AFrames_End(out);			// End of Anchored Frame
      FM_TextFlow_Set(out);			// Set TextFlow
      FM_Paragraph_Set(out);			// Set Paragraph
      out << "   <String `Norms         '> \n";	// Write its norms
      FM_Paragraph_End(out, 11);		// End Paragraph
      FM_TextFlow_End(out);

      FMP.DataUse(du);				// Reset data use value
      break;
    }
  FM_End(out);					// FrameMaker end comment
  }

void FrameMaker1D(const std::string& name, const row_vector& vx, 
                                              int ri, double xmin, double xmax)
  {
  std::string fname = name + std::string(".mif");	// Insure file is .mif
  double xsize = 14;					// Set width  of 14 cm
  double ysize = 14;					// Set height of 14 cm
  FM_1D(fname,vx,xsize,ysize,xmin,xmax,ri);

  std::string pltcmd = std::string("\"")
                     + FMExec()
                     + std::string("\" ")
                      + fname	// Command to execute macro
                 + std::string("\n");
  system(pltcmd.c_str());				// Plot to screen now

  }

void FrameMaker1D(const std::string& name, const std::vector<row_vector>& vxs, 
                                              int ri, double xmin, double xmax)
  {
  std::string fname = name + std::string(".mif");	// Insure file is .mif
  double xsize = 14;					// Set width  of 14 cm
  double ysize = 14;					// Set height of 14 cm
  FM_1Dm(fname,vxs,xsize,ysize,xmin,xmax,ri);
  }

// ____________________________________________________________________________
//                     FrameMaker Multiple 1D Plot Functions
// ____________________________________________________________________________


/*************************************************************************
**									**
** 		         FrameMaker 1D-MultiPlots			**
**									**
** Author	: SOSI							**
** Description  : Takes a pointer to a set of row_vectors as input	**
**		  along with plot parameters to produce	a FrameMaker	**
**		  MIF file containing a 1D plot	with all vectors 	**
**		  graphed on the same axes				**
**									**
*************************************************************************/

// Functions FM_1D_rim are derived from the function FM_1D_ri and are newer.
// Some day FM_1D_ri should be replaced by these two funcitons which 
// are more general.


	// Input	out       : Output stream
	// 		vx        : Data vector
	//		left      : Plot x start in cm
	//		top       : Plot y start in cm
	//		xscale    : Scaling factor (pts to cm)
	//		yscale    : Scaling factor (vx units to cm)
	//		ymax      : Largest y (vx units)
	//		PL	  : Vector for Polyline plotting
	//		PLmax	  : Maximum size of PL 
	//		PL_ID	  : Polyline ID number
        //		ri	  : Flag for real versus imaginary plot
	//			    0=reals, non-zero=imaginaries
	// Return		  : Void, file out is modified
	// Note			  : AUXILIARY FUNCTION FOR FM_1Dm
	// Note			  : Maximum Polyline size is currently
	//			    set to 128 points by integer PLmax
 
void FM_1D_rim(std::ostream &out, const row_vector& vx, double left, double top,
                double xscale, double yscale, double ymax,
			          row_vector &PL, int PLmax, int PL_ID, int ri)
  {
  int np = vx.elements();		// Vector size

  double xcm, ycm, yin;
  xcm = left;				// Initial point, 1st Polyline
  yin = ymax*yscale + top;
  if(ri)
    ycm = yin - vx.getIm(0)*yscale;
  else
    ycm = yin - vx.getRe(0)*yscale;

  int npl = 0;				// Polyline # count
  int pts = 0;				// Polyline pt count
  int skip = 0;
  PL.put(complex(xcm, ycm), pts);

  for(int j=1; j<np; j++)		// Go through each point, get
    {					// coordinates in cm
    xcm += xscale;			// x(j) = j*xscale+left
    if(ri)				// Treat Imaginaries
      ycm = yin - vx.getIm(j)*yscale;
    else				// Treat Reals
      ycm = yin - vx.getRe(j)*yscale;
    if(pts == 0)			// Add pt if 2nd in Polyline
      {
      pts++;				
      PL.put(complex(xcm, ycm), pts);
      skip = 0;
      }
    if(j == np-1)			// Add pt if last in Polyline
      {
      pts++;				// or last in data vector
      if(skip)
        {
        PL.put(complex(xcm-xscale, Im(PL(pts-1))), pts);
        pts++;
        }
      PL.put(complex(xcm, ycm), pts);
      skip = 0;
      }
    else if((fabs(Im(PL(pts)) - ycm))	// Add pt if height changes
                               >0.0001) // (sets vertical resolution!)
      {
      pts++;
      if(skip)
        {
        PL.put(complex(xcm-xscale, Im(PL(pts-1))), pts);
        pts++;
        }
      PL.put(complex(xcm, ycm), pts);
      skip = 0;
      }
    else
      skip = 1;

    if(pts >= PLmax-1)			// Limit Polylines to PLmax pts
      {
      FM_PolyLine(out, PL, PL_ID, 15, pts+1, 0);
      npl++;
      PL.put(PL(pts), 0);
      pts = 0;
      skip = 0;
      }
    }
  if(pts >= 1)				// Output Polyline if end
    {
    FM_PolyLine(out,PL,PL_ID,15,pts+1,0);
    npl++;
    }
  if(npl>1) FM_Group (out, PL_ID);	// Group Polylines this vector
  }



	// Input	out       : Output stream
	// 		nvec      : Number of data vectors
	// 		vx        : Data vector(s)
	//		xsize     : Plot horizontal (x) dimension in cm
	//		ysize     : Plot vertical (y) dimension in cm
	//		pmax      : Plot vertical maximum (axis label)
	//		pmin      : Plot vertical minimum (axis label)
        //		ri	  : Flag for real versus imaginary plot
	//			    0=reals, non-zero=imaginaries
	// Return		  : Void, file out is modified
	// Note			  : AUXILIARY FUNCTION FOR FM_1Dm
	// Note			  : Maximum Polyline size is currently
	//			    set to 128 points by integer PLmax
 
//void FM_1D_rim(std::ostream &out, int nvec, row_vector *vx, double xsize,
//	                   	double ysize, double pmin, double pmax, int ri)
void FM_1D_rim(std::ostream &out, const std::vector<row_vector>& vxs, double xsize,
	                   	double ysize, double pmin, double pmax, int ri)
  {
  double ymin=0, ymax=0;
  double ymini, ymaxi;
  ymini = HUGE_VAL;
  ymaxi = -HUGE_VAL;
  int xmax, xmaxi;
  xmax = 1;
  int i;
  int nvec = vxs.size();
  for(i=0; i<nvec; i++)
//    if(vx[i])
      {
      FM_maxima(vxs[i], ymini, ymaxi, ri);	// Get vertical maxima
      if(ymaxi > ymax) ymax=ymaxi;
      if(ymini < ymin) ymin=ymini;
      xmaxi = (vxs[i]).elements();		// Total points this vector
      if(xmaxi > xmax) xmax = xmaxi;
      }
  if(xmax <= 1)				// Must have more than 1 point
    {
    FM_error(5);			// Error - vector of 0 or 1 pt. 
    FM_fatality(100);			// Fatal error in 1D plot
    }

  double ydel = ymax - ymin;		// Begin scaling
  if(ydel == 0) ydel = ymax;
  int np = xmax;
  double top, bottom, right, left;	// Declare and set margins
  FM_borders(xsize, ysize,
             top, bottom, right, left);
  double yscale =			// Scaling in x and y directions
             (bottom-top)/(ymax-ymin);
  double xscale = (right-left)/(np-1);

  int FrameID = 11;			// Set Frame and Polyline ID's 
  int PL_ID = 1;
  if(ri)
    {
    FrameID++;
    PL_ID = 2;
    }
  FM_AFrame_Set(out, xsize, ysize	// Header of Frame
			     ,FrameID);
					// Begin plotting data
  int PLmax = 128;			// Limit size of Polyline(s)
  row_vector PL(PLmax+1);		// Polyline Vec. 128+1 pt(max.)

//                Finally, We Plot Each Vector Successively

  for(i=0; i<nvec; i++)					// Loop vectors
    {
    FM_1D_rim(out, vxs[i], left, top,
          xscale, yscale, ymax, PL, PLmax, PL_ID, ri);
    PL_ID++;						// New ID number
    }
  FM_Axis (out, 'x', bottom, left, right, pmin, pmax);	// Horizontal Axis
  FM_Axis (out, 'y', left, bottom, top, ymin, ymax);	// Vertical Axis
  FM_AFrame_End (out);					// End of Frame
  }

// ____________________________________________________________________________
// sosik




	// Input	filename  : Output filename
	// 		nvec      : Number of data vectors
	// 		vx        : Data vector
	//		xsize     : Plot horizontal (x) dimension in cm
	//		ysize     : Plot vertical (y) dimension in cm
	//		pmax      : Plot vertical maximum (axis label)
	//		pmin      : Plot vertical minimum (axis label)
        //		ri	  : Flag for real versus imaginary plot
	//			    ri = 0  : reals only
	//			    ri = <0 : imaginaries only
	//			    ri = >0 : both 
	// Return		  : Void, file out is modified
	// Note			  : This function cannot take an array of
	//			    block_1D's!!

void FM_1Dm(const std::string& filename, int nvec, row_vector *vx, double xsize,
	                        double ysize, double pmin, double pmax, int ri)
  {
  std::vector<row_vector> vxs;			// Array of vectors
  for(int i=0; i<nvec; i++)			// Loop over input vectors
    vxs.push_back(vx[i]);			// store in array
  FM_1Dm(filename,vxs,xsize,ysize,pmin,pmax,ri);// Use overload function
  }

void FM_1Dm(const std::string& filename, const std::vector<row_vector>& vxs, double xsize,
	                         double ysize, double pmin, double pmax, int ri)
  {
  if((xsize < 3) || (xsize > 20)) xsize = 14;	// Insure plot sizes sensible
  if((ysize < 3) || (ysize > 27)) ysize = 14;
  std::ofstream out(filename.c_str());		// Open file for plotting
  FM_Begin(out);				// FrameMaker begin comment
  FM_AFrames_Begin(out);			// Anchored Frame begin
  if(ri >= 0)
    FM_1D_rim(out,vxs,xsize,ysize,pmin,pmax,0); // Plot the reals
  if(ri)
    FM_1D_rim(out,vxs,xsize,ysize,pmin,pmax,1); // Plot the imaginaries
  FM_AFrames_End(out);				// End of AFrame
  FM_TextFlow_Set(out);				// Set TextFlow
  if(ri >= 0)
    {
    FM_Paragraph_Set(out);			// Set real Paragraph
    out << "   <String `Reals       '> \n";
    FM_Paragraph_End(out, 11);			// End real Paragraph
    }
  if(ri)
    {
    FM_Paragraph_Set(out);			// Set imag Paragraph
    out << "   <String `Imaginaries '> \n";
    FM_Paragraph_End(out, 12);			// End imag Paragraph
    }
  FM_TextFlow_End(out);
  FM_End(out);					// FrameMaker end comment
  }


// ____________________________________________________________________________
//                     FrameMaker Parametric Plot Functions
// ____________________________________________________________________________


/*************************************************************************
**									**
** 		              FrameMaker xy Plots			**
**									**
** Author	: TILO, SOSI, BEME 					**
** Description  : Takes a row_vector from input along with plot		**
**                parameters and produces a FrameMaker file containing	**
**		  an xy plot of the vector				**
**									**
*************************************************************************/

void FM_xyPlot(const std::string& filename, row_vector &vx, double xsize,
             double ysize, double xmin, double xmax, double ymin, double ymax)

	// Input	filename  : Output filename
	// 		vx        : Data vector
	//		xsize     : Plot horizontal (x) dimension in cm
	//		ysize     : Plot vertical (y) dimension in cm
        //              xmin      : Minimum value for x-scale
        //              xmax      : Maximum value for x-scale
        //              ymin      : Minimum value for y-scale
        //              ymax      : Maximum value for y-scale
	// Return		  : Void. xy-plot is output to file

  {
  int np = vx.elements();		// Total number of points
  if((xsize < 5) || (xsize > 20))
			    xsize = 14;	// Insure plot sizes sensible
  if((ysize < 5) || (ysize > 27))
			    ysize = 14;
					// Begin scaling
  if ( ymin == 0 && ymax == 0 )
	FM_maxima(vx, ymin, ymax, 1);	// Get vertical maxima if not specified
  if ( xmin == 0 && xmax == 0 )
  	FM_maxima(vx, xmin, xmax, 0);	// Get horizontal maxima if not specified
  double top, bottom, right, left;	// Declare and set margins
  FM_borders(xsize, ysize,
             top, bottom, right, left);
         				// Scaling in x and y direction
  double yscale = (top-bottom)/(ymax-ymin);
  double xscale = (right-left)/(xmax-xmin);
					// Begin FrameMaker Plot
  std::ofstream out(filename.c_str());	// Open new file for plotting
  FM_Begin(out);			// FrameMaker begin comment
  FM_AFrames_Begin(out);		// Anchored Frame begin
  FM_AFrame_Set(out,xsize,ysize);	// Header of Frame

  int j;				// Looping Index
  int npl = 128;			// Maximum points per PolyLine
  int anpl;				// Actual points per Polyline
  for (int i=0; i<np-1; i+=npl-1)	// Loop over all points plotting
    {					// in chunks of anpl points
    anpl = (npl<np-i)?npl:np-i;
    out << "  <PolyLine \n";
    out << "    <Fill 15> \n";	// Set Polyline transparent
    if (np>npl)
	out << "    <GroupID 1> \n";	// Set an ID if multiple Polylines
    out << Gform( "    <NumPoints %d>\n",anpl);
    for (j=i; j<i+anpl; j++)
      {
      out << Gform("    <Point %3.3f cm ", (vx.getRe(j)-xmin)*xscale+left);
      out << Gform("%3.3f cm>\n", (vx.getIm(j)-ymin)*yscale+bottom);
      }
    out << "    > # end of PolyLine \n";
    }
  if (np>npl)				// Group all Polylines
    FM_Group (out, 1);
  FM_Axis (out, 'x', bottom,
       left, right, xmin, xmax);	// Output x axis
  FM_Axis (out, 'y', left,
       bottom, top, ymin, ymax);	// Output y axis
  FM_AFrame_End (out);			// End of Frame
  FM_AFrames_End(out);			// Anchored Frame end
  FM_ParaText_End (out);		// Create a TextFlow with Frame
  FM_End(out);				// FrameMaker end comment
  }


void FM_xyPlot(const std::string& filename, row_vector &vx, FMxy& FMXY)

	// Input	filename  : Output filename
	// 		vx        : Data vector
	//		FMXY	  : FrameMaker xy-plot controls
	// Return		  : Void. xy-plot is output to file

  {
//	     Check Plot Parameters For Reasonable Settings

  int np = vx.elements();			// Total number of points
//  if((FMXY.hsize<5) || (FMXY.hsize>FMPAGEWIDTH))// Insure plot sizes sensible
    FMXY.hsize = FMPAGEWIDTH;
//  if((FMXY.vsize<5) || (FMXY.vsize>FMPAGEHEIGHT))
    FMXY.vsize = FMPAGEHEIGHT;
  if(FMXY.PLmax <= 3) FMXY.PLmax = FMPLmax;	// Insure reasonable PL size

//	              Set Plot Data Scaling Values

  if(FMXY.vmin==0 && FMXY.vmax==0)		// If no vertical scaling
    FM_maxima(vx,FMXY.vmin,FMXY.vmax,1);	// specified, get vx maxima
  if(FMXY.hmin==0 && FMXY.hmax==0 )		// If no horizontal scaling
    FM_maxima(vx,FMXY.hmin,FMXY.hmax,0);	// specified, get vx maxima
  double top, bot, right, left;			// Declare and set margins
  FM_borders(FMXY.hsize, FMXY.vsize,
                         top, bot, right, left);
  double vdel = FMXY.vmax - FMXY.vmin;		// Data vertical span
  double hdel = FMXY.hmax - FMXY.hmin;		// Data horizontal span
  double vscf=(top-bot)/vdel; 			// Vertical scaling factor
  double hscf=(top-bot)/hdel; 			// Horizontal scaling factor

//		Begin Output File With FrameMaker Plot

  std::ofstream out(filename.c_str());		// Open new file for plotting
  FM_Begin(out);				// FrameMaker begin comment
  FM_AFrames_Begin(out);			// Anchored Frame begin
  FM_AFrame_Set(out,FMXY.hsize,FMXY.vsize);	// Header of Frame

//	      Plot All Vector Points Into FM Output File

  int j;					// Looping Index
  int npl = FMXY.PLmax;				// Max. points per PolyLine
  int anpl;					// Actual points per Polyline
  for(int i=0; i<np-1; i+=npl-1)		// Loop over all data points
    {						// in chunks of the PL size
    anpl = (npl<np-i)?npl:np-i;			//	Length of output PL
    out << "  <PolyLine \n";			//	Start a new PL
    out << "    <Fill 15> \n";			//      Set PL transparent
    if(np>npl) out << "    <GroupID 1> \n";	// 	Set PL ID (many PLs)
    out << Gform( "    <NumPoints %d>\n",anpl);	//	Set PL size
    for(j=i; j<i+anpl; j++)			//	Loop PL points
      {
      out << Gform("    <Point %3.3f cm",	(vx.getRe(j)-FMXY.hmin)*hscf+left);
      out << Gform("%3.3f cm>\n",          (vx.getIm(j)-FMXY.vmin)*vscf+bot);
      }
    out << "    > # end of PolyLine \n";	// 	End the PL
    }

//		End Output File With FrameMaker Plot

  if(np>npl) FM_Group(out, 1); 			// Group all Polylines
  FM_Axis(out, 'x', bot, left, right, 		// Output x axis
                          FMXY.hmin, FMXY.hmax);
  FM_Axis(out, 'y', left, bot, top,		// Output y axis
                          FMXY.vmin, FMXY.vmax);
  FM_AFrame_End (out);				// End of Frame
  FM_AFrames_End(out);				// Anchored Frame end
  FM_ParaText_End(out);				// Create a TextFlow with Frame
  FM_End(out);					// FrameMaker end comment
  }

// ____________________________________________________________________________
//                     FrameMaker Scatter Plot Functions
// ____________________________________________________________________________


/*************************************************************************
**									**
** 		           FrameMaker Scatter Plots			**
**									**
** Author	: SOSI 							**
** Description  : Takes a vector from input along with plot parameters	**
**		  and produces a FrameMaker file containing a scatter 	**
**		  plot of the vector.					**
**									**
*************************************************************************/


void FM_scatter(const std::string& fout,row_vector &vx,char z,double xsz,double ysz)

	// Input	fout	: FrameMaker mif file name
	// 		vx      : Data vector containing points
	//		z	: Character to use to mark points
	// 		xsz 	: Horizontal plot size in cm.
	// 		ysz	: Vertical plot size in cm.
	// Return		: Void, FrameMaker plot file output

  {
  if((xsz < 5)||(xsz > 20)) xsz = 14;		// Insure plot sizes sensible
  if((ysz < 5)||(ysz > 27)) ysz = 14;
  double ymin, ymax, xmax, xmin;		// Calculate maxima and minima
  FM_maxima(vx, ymin, ymax, 1);			// Get vertical maxima
  FM_maxima(vx, xmin, xmax, 0);			// Get horizontal maxima
  double top, bottom, right, left;		// Declare and set margins
  FM_borders(xsz,ysz,top,bottom,right,left);	// Set plot borders
  double yscale = (top-bottom)/(ymax-ymin); 	// Compute scaling factors
  double xscale = (right-left)/(xmax-xmin); 	// Begin FrameMaker Plot
  std::ofstream out(fout.c_str());			// Open file
  FM_Begin(out);				// FrameMaker begin comment
  FM_AFrames_Begin(out);			// Anchored Frame begin
  FM_AFrame_Set (out, xsz, ysz);		// Header of Frame

  int np = vx.elements(); 			// Get total points plotted
  double x, y;					// For point x,y coordinates
  for(int i=0; i<np; i++)			// Loop over all points
    {
    x = (vx.getRe(i)-xmin)*xscale+left;		//   Point x coordinate
    y = (vx.getIm(i)-ymin)*yscale+bottom;	//   Point y coordinate
    FM_TextLine(out, 72, x, y, z);		//   Draw symbol z at point
    }
  FM_Axis (out,'x',bottom,left,right,xmin,xmax);// Output x axis
  FM_Axis (out,'y',left,bottom,top,ymin,ymax);	// Output y axis
  FM_Group(out, 72, 73);			// Group all points together
  FM_AFrame_End(out);				// Tail of Frame
  FM_AFrames_End(out);				// Anchored Frame end
  FM_ParaText_End(out);				// End Frame TextFlow
  FM_End(out);					// FrameMaker End
  }


void FM_scatter(const std::string& filename, row_vector &vx, char z, FMxy& FMXY)

	// Input	filename: Output filename
	// 		vx      : Data vector
	//		z	: Character to use to mark points
	//		FMXY	: FrameMaker xy-plot controls
	// Return		: Void, FrameMaker plot file output

  {
//	     Check Plot Parameters For Reasonable Settings

  int np = vx.elements();			// Total number of points
//  if((FMXY.hsize<5) || (FMXY.hsize>FMPAGEWIDTH))// Insure plot sizes sensible
    FMXY.hsize = FMPAGEWIDTH;
//  if((FMXY.vsize<5) || (FMXY.vsize>FMPAGEHEIGHT))
    FMXY.vsize = FMPAGEHEIGHT;
  if(FMXY.debug)
    {
    std::cout << "\n\n\tAccessing FrameMaker FM_scatter Function";
    std::cout << "\n\t\tInput Vector of "<< np << " Points";
    std::cout << "\n\t\tPlot Horizontal Width of "<< FMXY.hsize << " Centimeters";
    std::cout << "\n\t\tPlot Vertical Width of "<< FMXY.vsize << " Centimeters";
    }

//	              Set Plot Data Scaling Values

  double vxvmin, vxvmax;			// Data vertical maxima
  FM_maxima(vx, vxvmin, vxvmax, 1);		// specified, get vx maxima
  double vxhmin, vxhmax;			// Data vertical maxima
  FM_maxima(vx, vxhmin, vxhmax, 0);		// specified, get vx maxima
  double top, bot, right, left;			// Declare and set margins
  FM_borders(FMXY.hsize, FMXY.vsize,
                         top, bot, right, left);
  double vdel = vxvmax - vxvmin;		// Data vertical span
  double hdel = vxhmax - vxhmin;		// Data horizontal span
  double vscf=(top-bot)/vdel; 			// Vertical scaling factor
  double hscf=(right-left)/hdel; 		// Horizontal scaling factor

  if(FMXY.debug)
    {
    std::cout << "\n\t\tPlot Data Vertical Maximum of "<< vxvmax;
    std::cout << "\n\t\tPlot Data Vertical Minimum of "<< vxvmin;
    std::cout << "\n\t\tPlot Data Horizontal Maximum of "<< vxhmax;
    std::cout << "\n\t\tPlot Data Horizontal Minimum of "<< vxhmin;
    std::cout << "\n\t\tPlot Vertical Scaling Factor "<< vscf << " Centimeters";
    std::cout << "\n\t\tPlot Horizontal Scaling Factor "<< hscf << " Centimeters";
    }

//		Begin Output File With FrameMaker Plot

  std::ofstream out(filename.c_str());		// Open new file for plotting
  FM_Begin(out);				// FrameMaker begin comment
  FM_AFrames_Begin(out);			// Anchored Frame begin
  FM_AFrame_Set(out,FMXY.hsize,FMXY.vsize);	// Header of Frame

//	      Plot All Vector Points Into FM Output File

  double x, y;					// Use for x,y values
  for(int i=0; i<np; i++)			// Loop vector points
    {
    x = (vx.getRe(i)-vxhmin)*hscf + left;	// 	FM x coordinate
    y = (vx.getIm(i)-vxvmin)*vscf + bot;	//	FM y coordinate
    FM_TextLine(out, 72, x, y, z);		//	Output pt tagged z
    }
  FM_Axis(out, 'x', bot, left, right, 		// Output x axis
                          FMXY.hmin, FMXY.hmax);
  FM_Axis(out, 'y', left, bot, top,		// Output y axis
                          FMXY.vmin, FMXY.vmax);
  FM_Group(out, 72, 73);			// Group points together

//	             Finish Up Output MIF File & Close

  FM_AFrame_End(out);				// Tail of Frame
  FM_AFrames_End(out);				// Anchored Frame end
  FM_ParaText_End(out);				// Create a TextFlow with Frame
  FM_End(out);					// FrameMaker end comment
  }



void FM_scatterm(const std::string& filename, int nvec, row_vector *vx, char z, FMxy& FMXY)

	// Input	filename: Output filename
	//		nvec	: Number of input vectors
	// 		vx      : Array of data vectors
	//		z	: Character to use to mark points
	//		FMXY	: FrameMaker xy-plot controls
	// Return		: Void, FrameMaker plot file output

  {
//	     Check Plot Parameters For Reasonable Settings

//  if((FMXY.hsize<5) || (FMXY.hsize>FMPAGEWIDTH))// Insure plot sizes sensible
    FMXY.hsize = FMPAGEWIDTH;
//  if((FMXY.vsize<5) || (FMXY.vsize>FMPAGEHEIGHT))
    FMXY.vsize = FMPAGEHEIGHT;
  if(FMXY.debug)
    {
    std::cout << "\n\n\tAccessing FrameMaker FM_scatterm Function";
    std::cout << "\n\t\t" << nvec << " Data Vectors Input";
    std::cout << "\n\t\tPlot Horizontal Width of "<< FMXY.hsize << " Centimeters";
    std::cout << "\n\t\tPlot Vertical Width of "<< FMXY.vsize << " Centimeters";
    }

//	                  Set Plot Data Scaling Values

  int i;
  double vxvmin=0, vxvmax=0;
  double vxvmini=HUGE_VAL, vxvmaxi=-HUGE_VAL;
  for(i=0; i<nvec; i++)				// Loop over all input vectors
    {
    FM_maxima(vx[i], vxvmini, vxvmaxi, 1);	// Vertical maxima vector i
    if(vxvmaxi > vxvmax) vxvmax=vxvmaxi;	// Set vertical maxima for
    if(vxvmini < vxvmin) vxvmin=vxvmini;	// all vectors
    }
  double vxhmin=0, vxhmax=0;
  double vxhmini=HUGE_VAL, vxhmaxi=-HUGE_VAL;
  for(i=0; i<nvec; i++)				// Loop over all input vectors
    {
    FM_maxima(vx[i], vxhmini, vxhmaxi, 0);	// Horizontal maxima vector i
    if(vxhmaxi > vxhmax) vxhmax=vxhmaxi;	// Set horizontal maxima for
    if(vxhmini < vxhmin) vxhmin=vxhmini;	// all vectors
    }
  double top, bot, right, left;			// Declare and set margins
  FM_borders(FMXY.hsize, FMXY.vsize,
                         top, bot, right, left);
  double vdel = vxvmax - vxvmin;		// Data vertical span
  double hdel = vxhmax - vxhmin;		// Data horizontal span
  double vscf=(top-bot)/vdel; 			// Vertical scaling factor
  double hscf=(right-left)/hdel; 		// Horizontal scaling factor

  if(FMXY.debug)
    {
    std::cout << "\n\t\tPlot Data Vertical Maximum of "<< vxvmax;
    std::cout << "\n\t\tPlot Data Vertical Minimum of "<< vxvmin;
    std::cout << "\n\t\tPlot Data Horizontal Maximum of "<< vxhmax;
    std::cout << "\n\t\tPlot Data Horizontal Minimum of "<< vxhmin;
    std::cout << "\n\t\tPlot Vertical Scaling Factor "<< vscf << " Centimeters";
    std::cout << "\n\t\tPlot Horizontal Scaling Factor "<< hscf << " Centimeters";
    }

//		Begin Output File With FrameMaker Plot

  std::ofstream out(filename.c_str());		// Open new file for plotting
  FM_Begin(out);				// FrameMaker begin comment
  FM_AFrames_Begin(out);			// Anchored Frame begin
  FM_AFrame_Set(out,FMXY.hsize,FMXY.vsize);	// Header of Frame

//	      Plot All Vectors Points Into FM Output File

  double x, y;					// Use for x,y values
  int ID = 72;					// Base group ID
  for(int plt=0; plt<nvec; plt++, ID++, z++)
    {
    for(i=0; i<vx[plt].size(); i++)		// Loop vector points
      {
      x = (vx[plt].getRe(i)-vxhmin)*hscf + left;// 	FM x coordinate
      y = (vx[plt].getIm(i)-vxvmin)*vscf + bot;	//	FM y coordinate
      FM_TextLine(out, ID, x, y, z);		//	Output pt tagged z
      }
    FM_Group(out, ID);				// Group vector points
    }
  FM_Axis(out, 'x', bot, left, right, 		// Output x axis
                          FMXY.hmin, FMXY.hmax);
  FM_Axis(out, 'y', left, bot, top,		// Output y axis
                          FMXY.vmin, FMXY.vmax);

//	             Finish Up Output MIF File & Close

  FM_AFrame_End(out);				// Tail of Frame
  FM_AFrames_End(out);				// Anchored Frame end
  FM_ParaText_End(out);				// Create a TextFlow with Frame
  FM_End(out);					// FrameMaker end comment
  }


void FM_scatter(const std::string& filename, row_vector &vx, int sides,
                                    double PGsize, double xsize, double ysize)

	// Input	filename  : FrameMaker mif file name
	// 		vx        : Data vector containing points
	//		sides     : Number of sides in Polygon
	//		size	  : Size of Polygon in cm.
	// 		xsize     : Horizontal plot size in cm.
	// 		ysize     : Vertical plot size in cm.
	// Return		  : Void, FrameMaker plot file output

  {
  int np = vx.elements();		// Total points
  if((xsize < 5)||(xsize > 20))
		    xsize = 14;		// Insure plot sizes sensible
  if((ysize < 5)||(ysize > 27))
		    ysize = 14;
  if((PGsize > xsize) ||		// Insure Polygon size sensible
      (PGsize > ysize) || (PGsize <=.001)) 
    PGsize = xsize/100;	
  double ymin, ymax, xmax, xmin;	// Calculate maxima and minima
  FM_maxima(vx, ymin, ymax, 1);		// Get vertical maxima
  FM_maxima(vx, xmin, xmax, 0);		// Get horizontal maxima
  double top, bottom, right, left;	// Declare and set margins
  FM_borders(xsize, ysize,
             top, bottom, right, left);
  double yscale =			// Declare & compute scaling factors
	      (top-bottom)/(ymax-ymin);
  double xscale =
	      (right-left)/(xmax-xmin);
					// Begin FrameMaker Plot
  std::ofstream out(filename.c_str());	// Open file
  FM_Begin(out);			// FrameMaker begin comment
  FM_AFrames_Begin(out);		// Anchored Frame begin
  FM_AFrame_Set (out, xsize, ysize);	// Header of Frame

  double x, y;
  for (int i=0; i<np; i++)
    {
    x = (Re(vx(i))-xmin)*xscale+left;
    y = (Im(vx(i))-ymin)*yscale+bottom;
    FM_Polygon(out, 72, x, y, PGsize, sides);
    }
  FM_Axis (out, 'x', bottom,
       left, right, xmin, xmax);	// Output x axis
  FM_Axis (out, 'y', left,
       bottom, top, ymin, ymax);	// Output y axis
  FM_Group(out, 72, 73);
  FM_AFrame_End(out);			// Tail of Frame
  FM_AFrames_End(out);			// Anchored Frame end
  FM_ParaText_End(out);			// Create a TextFlow with Frame
  FM_End(out);				// FrameMaker end comment
  }


// ____________________________________________________________________________
//                     FrameMaker Histogram Plot Functions
// ____________________________________________________________________________


/*************************************************************************
**									**
** 		         FrameMaker Histogram Plots			**
**									**
** Author	: SOSI 							**
** Description  : Takes a vector from input along with plot parameters	**
**		  and produces a FrameMaker file containing a histogram **
**		  plot of the vector.					**
**									**
*************************************************************************/

void FM_histogram(const std::string& filename, row_vector &vx, int bins,
                                                   double xsize, double ysize)

	// Input	filename  : FrameMaker mif file name
	// 		vx        : Data vector containing points
	//		bin       : Number of bins in histogram
	// 		xsize     : Horizontal plot size in cm.
	// 		ysize     : Vertical plot size in cm.
	// Return		  : Void, FrameMaker plot file output

  {
  int pts = vx.elements();
  if((xsize < 5)||(xsize > 20))
		    xsize = 14;	// Insure plot sizes sensible
  if((ysize < 5)||(ysize > 27))
		    ysize = 14;
  double xmin, xmax;		// Insure data has horizontal spread
  FM_maxima(vx, xmin, xmax, 0);	// (units here still that of vx)
  if(xmax == xmin)
    FM_fatality(11);
  if(bins<=0 || bins>pts)	// Insure number of bins sensible
    bins = pts;
  double width=(xmax-xmin)/bins;// Compute bin width (vx units)

  row_vector hist(bins, xmin);	// Decare histogram array & fill
  int bin;
  for(bin=0;bin<bins;bin++)	// in bin horizontal coordinates
    hist.put(hist(bin)+(bin+0.5)*width, bin);

  double x, y, xst, xfi;	// Loop over all points of vx
  int search;
  for(int j=0; j<pts; j++) 	// and find bin(s) where it belongs
    {
    xst = xmin;
    xfi = width;
    y = Im(vx(j));
    x = Re(vx(j));
    search = 1;
    bin=0;
    while((search) && (bin<bins))
      {
      if((x > xst) && (x < xfi))
        {
        hist.put(hist(bin)+complex(0, y), bin);
        search = 0;
        }
      else if(x == xmin)
        {
        hist.put(hist(bin)+complex(0, y), bin);
        search = 0;
        }
      else if(x == xmax)
        {
        hist.put(hist(bins-1) + complex(0, y), bins-1);
        search = 0;
        }
      else if(x == xst)
        {
        hist.put(hist(bin)+0.5*complex(0, y), bin);
        hist.put(hist(bin+1)+0.5*complex(0, y), bin+1);
        search = 0;
        }
      else if(x == xfi)
        {
        hist.put(hist(bin) + 0.5*complex(0, y), bin);
        hist.put(hist(bin+1) + 0.5*complex(0, y), bin+1);
        search = 0;
        }
      xst += width;
      xfi += width;
      bin++;
      }
    }

  double ymin, ymax;			// Calculate vertical maxima
  FM_maxima(hist, ymin, ymax, 1);
  double top, bottom, right, left;	// Declare and set margins
  FM_borders(xsize, ysize,		// these are in centimeters
             top, bottom, right, left);
  double yscale =			// Declare & compute scaling factors
	      (top-bottom)/(ymax-ymin);	// xscale, yscale in cm/vx(scales)
  double xscale =
	      (right-left)/(xmax-xmin);

  std::ofstream out(filename.c_str());	// Open file, begin plot
  FM_Begin(out);			// FrameMaker begin comment
  FM_AFrames_Begin(out);		// Anchored Frame begin
  FM_AFrame_Set (out, xsize, ysize);	// Header of Frame

  double height, xinit, yinit;		// Switch from original vx units
  width *= xscale;			// to centimeters and output bins
  xinit = left-width/2;
  yinit = bottom - (ymin*yscale);
  for (int i=0; i<bins; i++)
    {
    x = xinit + Re(hist(i))*xscale;
    y = yinit + Im(hist(i))*yscale;
    height = (ymin-Im(hist(i)))*yscale;
    FM_Rectangle(out, 20, x , y, width, height);
    }

  FM_Axis(out, 'x', bottom,
       left, right, xmin, xmax);	// Output x axis
//  FM_Line (out,8,0,1,left,		// Output x axis line
// 		   bottom,right,bottom);
//  FM_Axis_tics(out, 8, 'x', bottom,	// Output x axis tics
//	      left, right, xmin, xmax);
//  double delx = Re(hist(1))-Re(hist(0));
//  FM_Axis_tics(out, 8, 'x', bottom,	// Output x axis tics
//     left+0.5*width, right-0.5*width,
//	 xmin+0.5*delx, xmax-0.5*delx);
//  FM_Group (out, 8);			// Group x axis line, tics
  FM_Axis (out, 'y', left,
       bottom, top, ymin, ymax);	// Output y axis
  FM_Group (out, 20, 21);		// Group Histogram
  FM_AFrame_End(out);			// Tail of Frame
  FM_AFrames_End(out);			// Anchored Frame end
  FM_ParaText_End(out);			// Create a TextFlow with Frame
  FM_End(out);				// FrameMaker end comment
  }

// ____________________________________________________________________________
//		           FrameMaker Contour Plots
// ____________________________________________________________________________

/*************************************************************************
**									**
**			FrameMaker Contour Plots			**
**									**
** Author	: Serge Boentges, Scott Smith				**
** Modifications: Switched to C++ and adapted for GAMMA	by SOSI		**
** Description  : Takes a matrix from input along with plot parameters	**
**		  and produces a FrameMaker file containing a contour	**
**		  plot of the matrix.					**
**									**
**		  The original code by Serge was constructed for	**
**		  C-Pheonix which does rapid contour plots inter-	**
**		  actively.  This ended up in function "contour_pts"	**
**		  included below.  The GAMMA FrameMaker contouring	**
**		  "FM_contour" functions, have sacrificed a lot of the	**
**		  speed in Serge's original routine in order to reduce	**
**		  the size of the output FrameMaker mif file while	**
**		  maintaining code clarity and the ability to manip-	**
**		  ulate the contour plot graphically within FM.		**
**									**
** Notes	: There are several routines written to	simplify the	**
**		  code to FM_contour (& also FM_fcontour).  These are	**
**		  listed below with their intended purpose.		**
**									**
**		1. contour_setup - insures that the input parameters	**
**		    are sensible.  If not it determines sensible ones.	**
**		    Also sets up some basic plot parameters.		**
**		2. contour_extr - insures that the chosen contours are	**
**		    reasonable.						**
**		3. thresh_inc - determines the amount to increment in	**
**		    order to reach the next contour level. Used in	**
**		    cases when the levels are internally determined.	**
**		4. contour_levels - determines the contour levels to	**
**		    take.  Done either by an interal calculation based	**
**		    on CLI & CLM values or by an input array of desired	**
**		    contour levels.					**
**		5. FM_same - determines whether two points in an	**
**		    xy-plane should be considered equivalent.		**
**		6. contour_lines - Links together the contour points	**
**		    from "contour_pts" into crude contour lines.	**
**		7. contour_refine - Takes the crude contour lines from	**
**		    "contour_lines" and bridges them into complete	**
**		    contour lines.					**
**		8. contour_refine - Writes all contour lines in a 	**
**		    particular contour level out to the FM file as	**
**		    grouped PolyLines.					**
**		9. group_contours - Groups the contour levels together	**
**		    to form FrameMaker graphics objects			**
**									**
** Note: These routine were redone in mid-June, 1995, by Scott as they	**
**       were pretty messy and bug ridden.  I hope I've improved it.	**
**       Things seem to work fine on the Bruker UXNMR data sets its	**
**       been tested on.  It runs slow though, a 1K x 2K NOESY data set **
**       takes about 10 min on a SPARC 10 to complete, 9 levels with	**
**	 with about 40,000 crude points in the lower contours!		**
**									**
*************************************************************************/

void contour_setup(const matrix& mx, FMcont& FMCP)

	// Input		mx	: Data matrix
	//      		FMCP	: Contour plot parameters
	// Return		  	: Void, argument FMCP may be modified
	// Note			  	: AUXILIARY FUNCTION FOR FM_contour
	// Note				: This relates the following entities:
	//				     horizontal = x = mx.cols() = j
	//				     vertical   = y = mx.rows() = i
 
 { 
 int MAXCONTS = 20;				// Max. allowed contour levels
 if(FMCP.vsize<5 || FMCP.vsize> 27) 		// Insure plot sizes sensible
   FMCP.vsize = 15;				// Minimum allowed is 5x5 cm
 if(FMCP.hsize<5 || FMCP.vsize>20)		// Maximum allowed is 20x27 cm
   FMCP.hsize = 15;				// the default is 15x15 cm
 int v_dim = mx.rows();				// Get array dimensions
 int h_dim = mx.cols(); 			// & check mx large enough
 if(h_dim < 5) FM_fatality(12);			// Minimum 5x5 points allowed
 if(v_dim < 5) FM_fatality(13);
 
 FMCP.vscale = FMCP.vsize/v_dim;		// Set vertical scaling (cm/pt)
 FMCP.hscale = FMCP.hsize/h_dim; 		// Set horizon. scaling (cm/pt)

 double r;					// Compute global max & min,
 double max = -HUGE_VAL;				// check data vertical spread
 double min = HUGE_VAL;
 for(int i=0; i<v_dim; i++)			// Loop over entire matrix
   for(int j=0; j<h_dim; j++)
     {
     r = mx.getRe(i, j);
     if(r > max) max=r;
     if(r < min) min=r;
     }
 if(max == min) FM_fatality(7);			// Insure matrix not level!
 FMCP.dmax = max;				// Store the array maximum
 FMCP.dmin = min;				// Store the array minimum

 if(FMCP.steps<=0 || FMCP.steps>MAXCONTS) 	// Insure [1,20] contours
   FMCP.steps = 3;
    
 if(FMCP.CPN < -1) FMCP.CPN = -1;		// Insure CPN flag is
 else if(FMCP.CPN > 1) FMCP.CPN = 1;		// either -1, 0, or 1

 switch(FMCP.CPN)				// Check threshold, i.e. the
   {						// initial contour level
   case -1:					// Only NEGETIVE contours
     if(FMCP.thresh < FMCP.dmin)		// Insure catches some data
       {
       FM_error(9);				// Threshold below all data
       FMCP.thresh = FMCP.dmax;			// switch to appropriate level
       if(FMCP.thresh > 0) FMCP.thresh = 0;
       }
     break;
   case 0:					// Only POSITIVE contours
     if(FMCP.thresh > FMCP.dmax)		// Insure catches some data
       {
       FM_error(8);				// Threshold above all data
       FMCP.thresh = FMCP.dmin;			// switch to appropriate level
       if(FMCP.thresh < 0) FMCP.thresh = 0;
       }
     break;
   case 1:					// Both POS & NEG contours
     if(FMCP.thresh < 0)			// Threshold must start off
       FMCP.thresh = -FMCP.thresh;		// > 0.
     else if(FMCP.thresh == 0 ||
                      FMCP.thresh > FMCP.dmax)
       {
       if(FMCP.dmin <= 0) FMCP.thresh = FMCP.dmax/20;
       else FMCP.thresh = FMCP.dmin;
       FM_error(15);
       }
     break;
   }

 if((FMCP.CLM) && (!FMCP.CLI))			// Check CLM vs. CLI input
   FMCP.CLM = 1;				// Cannot set CLI=0 & CLM!=0

 if(FMCP.CLM < 1) FMCP.CLM = 1;			// Insure >= 1 CLM
 if(FMCP.CLI < 0) FMCP.CLI = -FMCP.CLI; 	// Contour inc. must be > 0
 else if (FMCP.CLI == 0)			// This flags internal setting
   {						// of CLI and forces CLM = 1
   if(FMCP.CPN < 0) FMCP.CLI = 			// Compute CLI for the 
         fabs(FMCP.dmin-FMCP.thresh)/FMCP.steps;// steps of linear increments
   else FMCP.CLI =
          (FMCP.dmax - FMCP.thresh)/FMCP.steps;	
   }
 if(FMCP.debug)
   {
   std::cout << "\n\n\tFrom contour_setup:";
   std::cout << "\n\t   Framemaker Width is " << FMCP.vsize << " cm";
   std::cout << "\n\t   Framemaker Height is " << FMCP.hsize << " cm";
   std::cout << "\n\t   Framemaker Width Scale is " << FMCP.hscale << " cm/point";
   std::cout << "\n\t   Framemaker Height Scale is " << FMCP.vscale << " cm/point";
   std::cout << "\n\t   Input Array Maximum is " << FMCP.dmax;
   std::cout << "\n\t   Input Array Minimum is " << FMCP.dmin;
   std::cout << "\n\t   Requested Number of Contour Steps is " << FMCP.steps;
   std::cout << "\n\t   Lowest Threshold is " << FMCP.thresh;
   std::cout << "\n\t   Contour Level Increment is " << FMCP.CLI;
   std::cout << "\n\t   Contour Level Modifier is " << FMCP.CLM;
   std::cout << "\n\t   Contour Positive/Negative Flag is " << FMCP.CPN;
   }
 }


void contour_setup(const matrix& mx, FMcont& FMCP, FMclev& FMCL)


	// Input		mx	: Data matrix
	//      		FMCP	: Contour plot parameters
	// 		        FMCL	: Contour level data
	// Return		  	: Void, arguments FMCP & FMCL
	//				  may be modified
	// Note			  	: AUXILIARY FUNCTION FOR FM_contour
	// Note				: This relates the following entities:
	//				     horizontal = x = mx.cols() = j
	//				     vertical   = y = mx.rows() = i
 
 { 
 contour_setup(mx, FMCP);			// Setup/check plot parameters 
 FMCL.vdim = mx.rows()*mx.cols();		// Vector (contour points) dimenson
 FMCL.vxi = row_vector(FMCL.vdim);		// Contour level coordinates
 FMCL.vxf = row_vector(FMCL.vdim);		// Contour level coordinates
 FMCL.FMID = 0;					// Zero contour ID number
 FMCL.CLmaxnum = 10000;				// Maximum allowed contour lines
 FMCL.penneg = 1;				// Number of negative contours
 FMCL.penpos = 2;				// Number of positive contours
 FMCL.conneg = 0;				// Number of negative contours
 FMCL.conpos = 0;				// Number of positive contours
 FMCL.posID = 101;				// Base ID for positive contours
 FMCL.negID = 99;				// Base ID for negative contours
 if(FMCP.debug)
   {
   std::cout << "\n\t   Framemaker Vector Dimension is " << FMCL.vdim;
   std::cout << "\n\t   Framemaker Base ID for Positive Contours is " << FMCL.posID;
   std::cout << "\n\t   Framemaker Base ID for Negative Contours is " << FMCL.negID;
   std::cout << "\n\t   Framemaker Pen for Positive Contours is " << FMCL.penpos;
   std::cout << "\n\t   Framemaker Pen for Negative Contours is " << FMCL.penneg;
   std::cout << "\n\t   Maximum Allowed Lines Per Contour is " << FMCL.CLmaxnum;
   }
 }


int contour_extr(FMcont& FMCP, int posneg, double& threshold, double& extremum)

	// Input	FMCP      : Contour plot parameters
	// 				posneg 	  : Positive or negative contouring
	//				treshold  : Initial contour level
	//				extremum  : Largest level in contouring
	// Return       stop      : Flag for stopping contouring
	// Note			  : If requested contours are all above
	//			    (or below) the array maximum, stop is
	//			    returned true
	// Note			  : The following values are set in
	//			    this routine -

	//             FMCP.CLI   -  Contour Level Increment (+ or -)
	//	       extremum   -  Max. Value Attained in Contouring
	//	       threshold  -  Level of 1st Contour 

  {
  int stop = 0;
  switch(posneg)
    {
    case 0:					// Positive (or increasing) contours
      if(FMCP.CPN < 0) 
				stop = 1;		// Stop if only negative contours
      else
         {
         FMCP.CLI = fabs(FMCP.CLI);		// Insure CLI increases contour level 
         extremum = threshold;			// Start with extremum at threshold

         for(int i=1; i<FMCP.steps; i++)	// Calculate extremum				  
           extremum += pow(FMCP.CLM, i-1)*FMCP.CLI;
					
         if(extremum > FMCP.dmax) 		// Insure it isn't greater than max
                         extremum = FMCP.dmax;
         }
      break;
    case 1:					// Negative (or decreasing) contours
      if(FMCP.CPN == 0) stop = 1;			// Stop if only positive contours
      else
        {
        if(FMCP.CPN > 0) 
					threshold = -FMCP.thresh;
        FMCP.CLI = -fabs(FMCP.CLI);		// Insure CLI decreases contour level
        extremum = threshold;			// Start with extremum at threshold

        for(int i=1; i<FMCP.steps; i++)		// Calculate extremum
          extremum += pow(FMCP.CLM, i-1) * FMCP.CLI;

        if(extremum<FMCP.dmin) 			// Insure it isn't less than min
                         extremum = FMCP.dmin;	// Insure it isn't less than min
         }
      break;
    }
  if(fabs(threshold) > fabs(extremum)) stop=1; // Check limits
  if(FMCP.debug)
    {
    std::cout << "\n\n\tFrom contour_extr:";
    if(posneg == 0) std::cout << "\n\t   Treating Positive Contours";
    else if(posneg==1) std::cout << "\n\t   Treating Negative Contours";
    else std::cout << "\n\t   Cannot Distinguish Positive From Negative Contours!!";
    std::cout << "\n\t   Contour Level Increment Adjusted To " << FMCP.CLI;
    std::cout << "\n\t   Threshold Set To " << threshold;
    std::cout << "\n\t   Extremum Set To " << extremum;
    std::cout << "\n\t   Stop Status is " << stop;
    }
  return stop;
  }
 
 

double thresh_inc(FMcont& FMCP, FMclev& FMCL, double fact=1.0)

	// Input	FMCP      : Contour plot parameters
	// 		FMCL	  : Contour level data
	// 		Ts 	  : Array of contour levels
	//		IDs	  : Array of contour FrameMaker IDs

  {
  double CLM = FMCP.CLM;		// Contour level modifier
  double CLI = FMCP.CLI;		// Contour level increment
  int level = FMCL.level;		// Contour level
  return pow(CLM,level)*CLI*fact;
  }
 

void contour_levels(FMcont& FMCP, FMclev& FMCL, double* Ts, int* IDs)

	// Input	FMCP      : Contour plot parameters
	// 		FMCL	  : Contour level data
	// 		Ts 	  : Array of contour levels
	//		IDs	  : Array of contour FrameMaker IDs
	// Note			  : These are generated levels
	//			    from CLI and CLM values

 {
 long lev = 0;					// Contour level
 int NCLs = 0;					// Number of contour lines
 int conpos=0;					// Number of positive contours
 int conneg=0;					// Number of negative contours
 double thmax = 0;				// Maximum threshold (determined)
 double thresh = FMCP.thresh;			// Initial threshold (specified)
 int stop = 0;
 double tval=0;					// Adjustment for 0 threshold
 for(int pn=0; pn<2; pn ++)			// Loop over +/- contours
   {
   lev = 0;					// Contour level number
   stop = contour_extr(FMCP,pn,thresh,thmax);	// Set FMCP.CLI,thresh,thmax
   while((lev<FMCP.steps) && !stop)		// Loop over contours
     {
     tval = thresh;				// Exact 0 contour disallowed
     if(tval == 0)				// so may slightly adjust value
       tval += thresh_inc(FMCP, FMCL, 0.001);	// towards next threshold (1/1000th)
     if(tval > 0)				// Positive contour
       {					// Increment counter
       conpos++;				// and set ID
       IDs[NCLs] = FMCL.posID + conpos;
       Ts[NCLs] = tval;				// Set actual threshold
       NCLs++;
       }
     else if(tval < 0)				// Negative contour
       {					// Increment counter
       conneg++;				// and set ID
       IDs[NCLs] = FMCL.negID - conneg;
       Ts[NCLs] = tval;				// Set actual threshold
       NCLs++;
       }
     thresh += thresh_inc(FMCP, FMCL);		// Increment contour height
     lev = lev + 1;				// Increment contour count
     if(fabs(thresh) > fabs(thmax)) stop=1;	// Check limits
     }						// Contour level looping
   }
 FMCL.conpos = conpos;				// Set no. of + contours
 FMCL.conneg = conneg;				// Set no. of - contours
 FMCP.steps = NCLs;				// Reset no. of - contours
 if(FMCP.debug)
   {
   std::cout << "\n\n\tFrom contour_levels:";
   std::cout << "\n\t   Total Contours to Calculate is " << FMCP.steps;
   std::cout << "\n\t   Total Positve Contours is " << FMCL.conpos;
   int k;
   for(k=0; k<FMCL.conpos; k++)
     std::cout << "\n\t     Positve Contour " << k
          << " at " << Ts[k] << ", ID of " << IDs[k];
   std::cout << "\n\t   Total Negative Contours is " << FMCL.conneg;
   for(k=0; k<FMCL.conneg; k++)
     std::cout << "\n\t     Negative Contour " << k
          << " at " << Ts[k+FMCL.conpos]
          << ", ID of " << IDs[k+FMCL.conpos];
    }
  return;
  }
 

void contour_levels(FMcont& FMCP, FMclev& FMCL, double* levels, double* Ts, int* IDs)

	// Input	FMCP      : Contour plot parameters
	// 		FMCL	  : Contour level data
	//		levels    : Desired contour levels
	// 		Ts 	  : Array of used contour levels
	//		IDs	  : Array of contour FrameMaker IDs
	// Note			  : These levels are user specified
	//			    directly in levels

 {
 double tval=0;					// Just a temp value
 int cl, cl2;					// Contour indices
 for(cl=0; cl<FMCP.steps; cl++)			// Loop over input contours
   Ts[cl] = levels[cl];				// and copy them to Ts array
 for(cl=0; cl<FMCP.steps-1; cl++)		// Loop over input contours
   for(cl2=cl+1; cl2<FMCP.steps; cl2++)		// and reorder from highest to
     {						// lowest
     if(Ts[cl2] > Ts[cl])
       {
       tval = Ts[cl];
       Ts[cl] = Ts[cl2];
       Ts[cl2] = tval;
       }
     }
// Now remove ones past max & min
// and repeat input values
// for(int

 int conpos=0;					// Number of positive contours
 int conneg=0;					// Number of negative contours
 for(cl=0; cl<FMCP.steps; cl++)			// Loop over input contours
   {
   if(Ts[cl] == 0)				// so may slightly adjust value
     {
     if(cl) Ts[cl] = Ts[0]/1000.0;		// towards next threshold (1/1000th)
     else   Ts[cl] = Ts[FMCP.steps]/1000.0;
     }
   if(Ts[cl] > 0)				// Positive contour
     {						// Increment counter
     conpos++;					// and set ID
     IDs[cl] = FMCL.posID + conpos;
     }
   else						// Negative contour
     {						// Increment counter
     conneg++;					// and set ID
     IDs[cl] = FMCL.negID - conneg;
     }
   }
 FMCL.conpos = conpos;				// Set no. of + contours
 FMCL.conneg = conneg;				// Set no. of - contours
 if(FMCP.debug)
   {
   std::cout << "\n\n\tFrom contour_levels:";
   std::cout << "\n\t   Total Contours to Calculate is " << FMCP.steps;
   std::cout << "\n\t   Total Positve Contours is " << FMCL.conpos;
   int k;
   for(k=0; k<FMCL.conpos; k++)
     std::cout << "\n\t     Positve Contour " << k
          << " at " << Ts[k] << ", ID of " << IDs[k];
   std::cout << "\n\t   Total Negative Contours is " << FMCL.conneg;
   for(k=0; k<FMCL.conneg; k++)
     std::cout << "\n\t     Negative Contour " << k
          << " at " << Ts[k+FMCL.conpos]
          << ", ID of " << IDs[k+FMCL.conpos];
   }
 return;
 }


void contour_pts(const matrix &mx, FMcont& FMCP, FMclev& FMCL)

	// Input	mx	  : Data matrix
	//      	FMCP      : Contour plot parameters
	// 		FMCL	  : Contour level data
	// Note			  : AUXILIARY FUNCTION FOR FM_contour

  {
  int ip, ib;					// For triangle pt1 coords.
  int jp, jb;					// For triangle pt2 coords.
  double p,					// First point   p(i  ,j  )
  	 a,					// Second point  a(i  ,j+1) or a(i+1, j)
  	 b;					// Third point   b(i+1,j+1)
  double x, y;					// Intersection point
  complex zi, zf;				// Points of contour
  int situation;				// Flag for contour section type
  int h_dim = mx.cols();			// First get dimensions
  int v_dim = mx.rows();
  double hscale = FMCP.hscale;
  double vscale = FMCP.vscale;
  double threshold = FMCL.threshold;
  int ptcnt=0;					// Zero the point count
//  int maxv = h_dim*v_dim/3; 			// Assume < 1/3 pts in contour
//  row_vector vxi(maxv);
//  row_vector vxf(maxv);
  for(int i=0; i<v_dim-1; i++)			// Loop over entire matrix
    {						// and over both upper and
    ip = i;					// Triangle pt1 x coord.
    ib = ip + 1;				// Triangle pt2 x coord.
    for(int j=0; j<h_dim-1; j++)
      {
      jp = j;					// Triangle pt1 y coord.
      jb = jp + 1;				// Triangle pt2 y coord.
      p = mx.getRe(ip, jp);			// Triangle first point
      b = mx.getRe(ib, jb);			// Triangle second point
      for(int trian=0; trian<2; trian ++)	// Upper & lower triangles
        {
        if(trian) a = mx.getRe(ib, jp);		// Triangle third point
        else a = mx.getRe(ip, jb);
        situation = 0;				// Determine contour type, if
        if(p < threshold) situation+=4;		// any, in current triangle
        if(b < threshold) situation+=2;
        if(a < threshold) situation+=1;

        switch(situation)			// Generate two points of contour
	  {
     	  case 2:				// p&a>t, b<t or p&a<t, b>t
      	  case 5:
	     x = ((threshold - p)*1.41421356)/(b - p);
	     x *= 0.70710678;
             zi = complex((j+x)*hscale, (i+x)*vscale);
      	     y = (threshold - a)/(b - a);
      	     if(trian) zf = complex((j+y)*hscale, (i+1)*vscale);
             else zf = complex((j+1)*hscale, (i+y)*vscale);
      	     break;

       	  case 3:				// p>t, a&b<t or p<t, a&b>t
      	  case 4:
	     x = (threshold - p)/(a - p);
      	     if(trian) zi = complex(j*hscale, (i+x)*vscale);
             else zi = complex((j+x)*hscale, i*vscale);
      	     y = ((threshold - p)*1.41421356)/(b - p);
	     y *= 0.70710678;
             zf = complex((j+y)*hscale, (i+y)*vscale);
      	     break;

	  case 1:			 	// p&b>t,a<t or p&b<t,a>t
      	  case 6:
	     x = (threshold - p)/(a - p);
      	     y = (threshold - a)/(b - a);
      	     if(trian)
               {
               zi = complex(j*hscale, (i+x)*vscale);
               zf = complex((j+y)*hscale, (i+1)*vscale);
               }
             else
               {
               zi = complex((j+x)*hscale, i*vscale);
               zf = complex((j+1)*hscale, (i+y)*vscale);
               }
             break;

          default: 				// All "points" in triangle
             situation = 0; 			// either above (0) or below (8)
             break;                             // threshold, so no contour
	  }
        if(situation)				// If contour point found, add 
          { 					// it to the list
//          vxi.put(zi,ptcnt);			// These are the initial points
//          vxf.put(zf,ptcnt);			// These are the final points
          (FMCL.vxi).put(zi,ptcnt);		// These are the initial points
          (FMCL.vxf).put(zf,ptcnt);		// These are the final points
          ptcnt++;
          }
        }					// Upper/Lower triangle loop
      }						// Data matrix element loop
//    if(ptcnt > maxv)				// Insure the number of points
    if(ptcnt > h_dim*v_dim)				// Insure the number of points
      {						// in the contour doesn't exceed
      FM_error(16);				// what we have allocated.
      FM_fatality(17);
      } 
    }
  if(FMCP.debug)				// If debugging, output this
    { 
    if(FMCL.level == 0) std::cout << "\n\n\tFrom FM_pts:";
    std::cout << "\n\t   Level is " << FMCL.level;
    std::cout << "\n\t     Points in Contour Level " << 2.0*ptcnt;
    std::cout << "\n\t     Contour Level is " << FMCL.threshold;
    std::cout << "\n\t     Contour FM ID is " << FMCL.FMID;
    } 
//  row_vector vxa(ptcnt);			// Now reduce the memory used
//  row_vector vxb(ptcnt);			// for these vectors 
//  for(i=0; i<ptcnt; i++)
//    {
//    vxa.put((FMCL.vxi).get(i),i);
//    vxb.put((FMCL.vxf).get(i),i);
//    }
  FMCL.npts = ptcnt;				// Set the point count
//  FMCL.vxi = vxa;				// Set the intiial points
//  FMCL.vxf = vxb;				// Set the final points
  }


int FM_same(const complex& z1, const complex& z2, double cutoff=0, int inches=0)

	// Input	z1	  : Point 1
	//		z2        : Point 2
	// 		cutoff    : Cutoff distance
	// 		inches    : Flag for cm vs. inches
	// Return	same	  : True if z1 & z2 are separated
	//			    by <= cutoff

// Tells whether the two dimensional points z1 & z2 are the same
// relative to a specified cutoff distance between them

  {
  int same = 0;
  double x = Re(z1) - Re(z2);		// Calculate z1, z2 delta x
  double y = Im(z1) - Im(z2);		// Calculate z1, z2 delta y
  double dist = sqrt(x*x+y*y);		// Calculate z1, z2 separation
  if(inches)
    {
    dist *= 2.54;			// Convert from inches to cm 
    cutoff *= 2.54;
    }
  if(!cutoff) cutoff = 0.001;		// 1 tenth of a mm
  if(dist <= cutoff) same = 1;
  return same;
  }


void contour_lines(FMcont& FMCP, FMclev& FMCL, row_vector* RCLs,
                          row_vector* LCLs, int* NRpts, int* NLpts, int& NCLs)

	// Input 	FMCP      : Contour plot parameters
	// 		FMCL	  : Contour level data
	//		RCLs      : Array of right contour lines
	//		LCLs      : Array of left contour lines
	//		NRpts     : Points in right contour lines
	//		NLpts     : Points in left contour lines
	//		NCLs	  : Number of contour lines
	// Return		  : Contour points in FMCL are
	//			    joined together in vectors.
	// Note			  : AUXILIARY FUNCTION FOR FM_contour

  {
  complex z, zi, zf;
  int connect = 0;				// Flag for point connection
  int j,k=0;
  double scale = 0;
  int start = 25;
  int expand = 25;
  if(FMCL.npts > 1000) expand=100;
  if(FMCL.npts > 10000) expand=500;
  for(int i=0; i<FMCL.npts; i++)		// Loop over entire matrix
    {
    zi = (FMCL.vxi).get(i);			// An intial contour point
    zf = (FMCL.vxf).get(i);			// A final contour point
    connect = 0;				// Flag that not connected
    for(j=0; j<NCLs && !connect; j++)		// Perhaps {zi,zf} connects with
      { 					// a previous contour line
      z = (RCLs[j]).get(NRpts[j]-1);		// Last point this contour line
      if(FM_same(z, zi, scale))			// If its the same, add it
        {
        RCLs[j].put(zf,NRpts[j]);		// Add zf to right contour line
        NRpts[j] = NRpts[j] + 1;		// Increment the point count
        connect = 1;				// Flag that it connected
        k = j;					// Add to this contour line
        }
      else if(FM_same(z, zf, scale))
        {
        RCLs[j].put(zi,NRpts[j]);		// Add zi to right contour line
        NRpts[j] = NRpts[j] + 1;		// Increment the point count
        connect = 1;				// Flag that it connected
        k = j;					// Add to this contour line
        }
      else					// If not already connected
        {					// try the left contour lines
        z = (LCLs[j]).get(NLpts[j]-1);		// Last point this contour line
        if(FM_same(z, zi, scale))
          {
          LCLs[j].put(zf,NLpts[j]);		// Add zf to left contour line
          NLpts[j] = NLpts[j] + 1;		// Increment the point count
          connect = 1;			// Flag that it connected
          k = j;				// Add to this contour line
          }
        else if(FM_same(z, zf, scale))
          {
          LCLs[j].put(zi,NLpts[j]);		// Add zi to left contour line
          NLpts[j] = NLpts[j] + 1;		// Increment the point count
          connect = 1;			// Flag that it connected
          k = j;				// Add to this contour line
          }
        }
      }
    if(!connect)				// Existing contour lines do not connect
      {						// to {zi,zf} so start a new contour line
      RCLs[NCLs]= row_vector(start);		// Vector for right contour line
      LCLs[NCLs]= row_vector(start);		// Vector for left contour line
//      RCLs[NCLs]= row_vector(FMCL.CLmaxsize-i);	// Vector for right contour line
//      LCLs[NCLs]= row_vector(FMCL.CLmaxsize-i);	// Vector for left contour line
      (RCLs[NCLs]).put(zf,0);			// Right contour line 1st pt is zf
      (LCLs[NCLs]).put(zi,0);			// Left contour line 1st pt is zi
      NRpts[NCLs] = 1;				// Right contour line size is 1
      NLpts[NCLs] = 1;				// Left contour line size is 1
      NCLs++;					// Increment the number of lines
      }
    else					// Expand the vector size
      {						// if it has gotten too large
      if(NRpts[k] == RCLs[k].size()-1)
        {
        row_vector vx(RCLs[k].size() + expand);
        for(int l=0; l<NRpts[k]; l++)
          vx.put(RCLs[k].get(l),l);
        RCLs[k] = vx;
        }
      if(NLpts[k] == LCLs[k].size()-1)
        {
        row_vector vx(LCLs[k].size() + expand);
        for(int l=0; l<NLpts[k]; l++)
          vx.put(LCLs[k].get(l),l);
        LCLs[k] = vx;
        }
      }
    if(NCLs >= FMCL.CLmaxnum)
      {
      std::cout << "\n\n\tToo Many Contour Lines!!"; 
      exit(-1);
      }
    }
  if(FMCP.debug)                                // If debugging, output this
    { 
    std::cout << "\n\t     Number of Raw Contour Lines is " << NCLs;
    int nr=0;
    int nl=0;
    int k=0;
    for(k=0; k<NCLs; k++)
      {
      nr += NRpts[k];
      nl += NLpts[k];
      }
    std::cout << "\n\t     Number of Raw Contour Points is "
         << nr+nl << " (" << nr << " Right, " << nl
         << " Left)";
    if(FMCP.debug > 1)
      {
      for(k=0; k<NCLs; k++)
        {
        std::cout << "\n\t        Line " << k << " : ";
        std::cout << NRpts[k] << " Right, " << NLpts[k] << " Left";
        if(FMCP.debug > 2)
          std::cout << "- " << (RCLs[k]).get(NRpts[k]-1)
               << ", " << (LCLs[k]).get(NLpts[k]-1);
        } 
      } 
    } 
  return;
  }


void contour_refine(FMcont& FMCP, FMclev& FMCL, row_vector* RCLs,
                   row_vector* LCLs, int* NRpts, int* NLpts, int& NCLs)

	// Input 	FMCP      : Contour plot parameters
	// 		FMCL	  : Contour level data
	//		RCLs      : Array of right contour lines
	//		LCLs      : Array of left contour lines
	//		NRpts     : Points in right contour lines
	//		NLpts     : Points in left contour lines
	//		NCLs	  : Number of contour lines
	// Return		  : Contour points in FMCL are
	//			    joined together in vectors.

  {
  complex zri, zli;
  complex zrj, zlj;
  complex zrix, zlix;
  double scale = 0;
  int dri=0, dli=0;
  int pri=0,pli=0,prj=0,plj=0;
  int i,j,k;
  row_vector vxnull;
  int m=0, n=0;
  int reduce = 0;
  for(i=0; i<NCLs-1; i++)		// Loop over all crude lines
    {
    pri = NRpts[i];			// Points in right vector
    pli = NLpts[i];			// Points in left vector
    zri = RCLs[i].get(pri-1);		// Get end of right line
    zli = LCLs[i].get(pli-1);		// Get end of left line
    dri = RCLs[i].size();		// Right vector dimension
    dli = LCLs[i].size();		// Left vector dimension
    for(j=i+1; j<NCLs; j++)
      {
      reduce = 0;
      prj = NRpts[j];			// Points in right vector
      plj = NLpts[j];			// Points in left vector
      zrj = RCLs[j].get(prj-1);		// Get end of right ljne
      zlj = LCLs[j].get(plj-1);		// Get end of left line
      if(FM_same(zri, zrj, scale))	// Add lines if bridge exists
        {
        if(dri < pri+prj+plj-1)		//	Expand line i if too
          {				//	small
          row_vector vx(pri+prj+plj-1);
          for(m=0; m<pri; m++)
            vx.put((RCLs[i]).get(m),m);
          RCLs[i] = vx; 
          }
        for(n=prj-2,m=pri;n>=0;n--,m++)	//	Add j right line to i
          (RCLs[i]).put((RCLs[j]).get(n),m);	
        for(n=0; n<plj; n++,m++)	//	Add j left line to i
          (RCLs[i]).put((LCLs[j]).get(n),m);
        zrix = zri;			//	For debugging
        zlix = zli;
        NRpts[i] = pri+prj+plj-1; 	// 	Reset pts in right vector
        pri = NRpts[i];			//	Working pts in right vector
        zri = RCLs[i].get(pri-1);	// 	Reset end of right line
        dri = RCLs[i].size();		// 	Reset right vector dimension
        reduce = 1;			// 	Flag for data reduction
        }
      else if(FM_same(zri, zlj, scale))	// Add lines if bridge exists
        {
        if(dri < pri+prj+plj-1)		//	Expand line i if too
          {				//	small
          row_vector vx(pri+prj+plj-1);
          for(m=0; m<pri; m++)
            vx.put((RCLs[i]).get(m),m);
          RCLs[i] = vx; 
          }
        for(n=plj-2,m=pri;n>=0;n--,m++)	//	Add j left line to i
          (RCLs[i]).put((LCLs[j]).get(n),m);
        for(n=0; n<prj; n++,m++)	//	Add j right line to i
          (RCLs[i]).put((RCLs[j]).get(n),m);
        zrix = zri;			//	For debugging
        zlix = zli;
        NRpts[i] = pri+prj+plj-1; 	// 	Reset pts in right vector
        pri = NRpts[i];			//	Working pts in right vector
        zri = RCLs[i].get(pri-1);	// 	Reset end of right line
        dri = RCLs[i].size();		// 	Reset right vector dimension
        reduce = 2;			// 	Flag for data reduction
        }
      else if(FM_same(zli, zrj, scale))
        {
        if(dli < pli+prj+plj-1)
          {
          row_vector vx(pli+prj+plj-1);
          for(m=0; m<pli; m++)
            vx.put((LCLs[i]).get(m),m);
          LCLs[i] = vx; 
          }
        for(n=prj-2,m=pli;n>=0;n--,m++)	//	Add j right line to i
          (LCLs[i]).put((RCLs[j]).get(n),m);
        for(n=0; n<plj; n++,m++)	//	Add j right line to i
          (LCLs[i]).put((LCLs[j]).get(n),m);
        zrix = zri;			//	For debugging
        zlix = zli;
        NLpts[i] = pli+prj+plj-1; 	// 	Reset pts in left vector
        pli = NLpts[i];			//	Working pts in left vector
        zli = LCLs[i].get(pli-1);	// 	Reset end of left line
        dli = LCLs[i].size();		// 	Reset right vector dimension
        reduce = 3;			// 	Flag for data reduction
        }
      else if(FM_same(zli, zlj, scale))
        {
        if(dli < pli+prj+plj-1)
          {
          row_vector vx(pli+prj+plj-1);
          for(m=0; m<pli; m++)
            vx.put((LCLs[i]).get(m),m);
          LCLs[i] = vx; 
          }
        for(n=plj-2,m=pli;n>=0;n--,m++)	//	Add j left line to i
          (LCLs[i]).put((LCLs[j]).get(n),m);
        for(n=0; n<prj; n++,m++)	//	Add j right line to i
          (LCLs[i]).put((RCLs[j]).get(n),m);
        NLpts[i] = pli+prj+plj-1; 	// 	Reset pts in left vector
        zrix = zri;			//	For debugging
        zlix = zli;
        pli = NLpts[i];
        zli = LCLs[i].get(pli-1);	// 	Reset end of left line
        dli = LCLs[i].size();		// 	Reset right vector dimension
        reduce = 4;			// 	Flag for data reduction
        }
      if(reduce)			// If contour lines joined, the
        {				// reduce the number we have
        if(FMCP.debug > 1)
          {
          std::cout << "\n\tBridged crude contour lines "
               << i << " and " << j << ", Type " << reduce;
//          cout << "\n\t\tCrude line " << i << ": "
//               << pri << " Right points, "
//               << pli << " Left points";
          if(FMCP.debug > 2)
            std::cout << " - " << zrix << ", " << zlix;
//          cout << "\n\t\tCrude line " << j << ": "
//               << prj << " Right points, "
//               << plj << " Left points";
          if(FMCP.debug > 2)
            std::cout << " - " << zrj << ", " << zlj;
//          cout << "\n\t\tNew line " << i << ": "
//               << NRpts[i] << " Right points, "
//               << NLpts[i] << " Left points";
          if(FMCP.debug > 2)
            std::cout << " - " << zri << ", " << zli;
          }
        for(k=j+1; k<NCLs; k++)		//	Move all other lines up
          {				//	by one.  This involves
          RCLs[k-1] = RCLs[k];		//	Both right and left lines
          LCLs[k-1] = LCLs[k];		//	as well as the number of
          NRpts[k-1] = NRpts[k];	//	points each contains
          NLpts[k-1] = NLpts[k];
          }
        RCLs[NCLs-1] = vxnull;		//	Empty last right line
        LCLs[NCLs-1] = vxnull;		//	Empty last left line
        NRpts[NCLs-1] = 0; 		// 	No points in last line
        NLpts[NCLs-1] = 0; 		// 	No points in last line
        NCLs--;				//	Decrease line count
        j--;				//	Decrease index for refinement
        }
      }
    }
  if(FMCP.debug)                                // If debugging, output this
    { 
    std::cout << "\n\t     Number of Refined Contour Lines is "
         << NCLs;
    int nr = 0;
    int nl = 0;
    for(k=0; k<NCLs; k++)
      {
      nr += NRpts[k];
      nl += NLpts[k];
      }
    std::cout << "\n\t     Number of Refined Contour Points is "
         << nr+nl << " (" << nr << " Right, " << nl
         << " Left)";
    if(FMCP.debug > 1)
      {
      for(k=0; k<NCLs; k++)
        {
        std::cout << "\n\t        Line " << k << " : ";
        std::cout << NRpts[k] << " Right, " << NLpts[k] << " Left";
        } 
      } 
    } 
  return;
  if(FMCL.CLmaxnum) k=0;		// Compiler likes this used
  }


void contour_output(std::ostream& out, FMcont& FMCP, FMclev& FMCL, row_vector* RCLs,
                   row_vector* LCLs, int* NRpts, int* NLpts, int& NCLs)

	// Input	out	  : Output file std::ostream
	//      	FMCP      : Contour plot parameters
	// 		FMCL	  : Contour level data
	//		PLmxr	  : Matrix of contour in right PolyLines
	//		PLmxl	  : Matrix of contour in left PolyLines
	//		PLptr	  : Points in right PolyLines
	//		PLptl	  : Points in left PolyLines
	//		PL	  : PolyLine being tested for closure
	//		CLnum	  : Number of contour lines
	// Return		  : Void, file out may be modified
	// 			  : Value of CLnum may be modified
	// Note			  : AUXILIARY FUNCTION FOR FM_contour

// All current contour lines are to be output.  A final attempt is made
// to connect those in close proximity with the hope of forming closed
// contours.  The value scale is used to decide whether the lines are
// connected within the original input data resolution.
  
  {
  double scale = 0;
  complex z;
  int cons;
  int pri, pli,pts;				// Pts in right,left,total
  int i,j,k; 					// Use for indices
  int comp=0, incomp=0;
  int pen = 0;					// Default pen is black
  if(FMCL.FMID >= FMCL.posID) pen = FMCL.penpos;// Pen for positive contours
  else pen = FMCL.penneg;			// Pen for negative contours
  for(i=0; i<NCLs; i++)				// Loop over all contour lines
    {
    pri = NRpts[i];				// Points in right vector
    pli = NLpts[i];				// Points in left vector
    pts = pri+pli;
    row_vector vx(pts)	;			// Vector for contour line output
    for(j=pli-1,k=0; j>=0; j--, k++)		// Fill vx from end of left 
      vx.put(LCLs[i].get(j),k);			// contour part to "middle"
    for(j=0; j<pri; j++, k++)			// Continue filling vx from start
      vx.put(RCLs[i].get(j),k);			// ("middle") of right contour part
    z = vx.get(0);				// Check that contour line is in
    cons = 1;				// fact not a constant point
    for(k=1; k<pts && cons; k++) 
    if(!FM_same(vx.get(k), z, scale))
      cons = 0;
    if(pts>1 && !cons)
      {
      FM_PolyLine(out,vx,FMCL.FMID,15,pts,pen);	// Output PolyLine , ID, clear fill
      if(FM_same(vx.get(0),vx.get(pts-1),scale))
        comp++;
      else incomp++;
      }
    }
  if(FMCP.debug)                                // If debugging, output this
    { 
    if(comp)
      std::cout << "\n\t     " << comp
           << " Complete Contour Lines Output";
    if(incomp)
      std::cout << "\n\t     " << incomp
           << " Unjoined Contour Lines Output";
    }
  return;
  }


void group_contours(std::ostream &out, int conpos, int conneg)

	// Input	out 	  : Output file std::ostream
	// 		conpos	  : Number of positive contours
	// 		conneg	  : Number of negative contours

  {
  int FMIDpos = 101;
  int FMIDneg = 99;
  int l;
  for(l=1; l<=conpos; l++)			// Each positive contour
    FM_Group(out, FMIDpos+l, FMIDpos);		// grouped with ID = 101
  for(l=1; l<=conneg; l++)			// Each negative contour
    FM_Group(out, FMIDneg-l, FMIDneg);		// grouped with ID = 99
  if((conpos) && (conneg))
    {						// If both +/- contours group all
    if(conpos>1) FM_Group(out, FMIDpos); 	// + together and all - together
    if(conneg>1) FM_Group(out, FMIDneg);
    }      
  }


void FM_contour(const std::string& filename, const matrix &mx,
	    double low_t, int steps, double CLI, double CLM,
			           int CPN, double hsize, double vsize)

	// Input	filename  : Output filename
	// 		mx        : Data matrix
	//		low_t     : Lowest threshold
	//		steps     : Number of contour steps
	//		CLI	  : Contour level increment
	//		CLM 	  : Contour level modifier
	//		CPN 	  : Positive/Negetive contour flag
	//		hsize     : Plot horizontal (x or j) dimension in cm
	//		vsize     : Plot vertical (y or i) dimension in cm
	// Return		  : Void, file out is modified
	// Note			  : AUXILIARY FUNCTION FOR FM_contour
 
 { 
 FMcont FMCP;					// Put input info in a structure
 FMCP.thresh = low_t; 				// Set the threshold
 FMCP.steps = steps;				// Set the number of contours
 FMCP.CLI = CLI;				// Set contour level increment
 FMCP.CLM = CLM;				// Set contour level modifier
 FMCP.CPN = CPN;				// Set positive/negative contour flag
 FMCP.hsize = hsize;				// Plot horizontal size (cm)
 FMCP.vsize = vsize;				// Plot vertical size (cm)
 FMCP.vlow = 0;					// Lower vertical axis value
 FMCP.vhigh = 1;				// Upper vertical axis value
 FMCP.hlow = 0;					// Lower horizontal axis value
 FMCP.hhigh = 1;				// Upper horizontal axis value
 FMCP.debug=0;					// Set for no debugging
 FM_contour(filename, mx, FMCP);		// Use function overload
 }


void FM_contourx(const std::string& filename, matrix &mx, FMcont& FMCP)

	// Input	filename  : Output filename
	// 		mx        : Data matrix
	//		FMCP      : Contour plot parameters
	// Return		  : Void, file out is modified
	// Note			  : AUXILIARY FUNCTION FOR FM_contour
 
 
 { 
// ----------------- Open File, Set Up FrameMaker Document ---------------- //

 contour_setup(mx, FMCP);			// Setup/check parameters 
 std::ofstream out(filename.c_str());			// Define output stream
 FM_Begin(out);					// FrameMaker begin comment
 FM_AFrames_Begin(out);				// Anchored Frame begin
 FM_AFrame_Set(out, FMCP.hsize, FMCP.vsize);	// Set Anchored Frame

// --------------------- Initialize Contour Parameters -------------------- //

 FMclev FMCL;					// Values for a contour level
 FMCL.vdim = mx.rows()*mx.cols();		// Vector (contour points) dimenson
 FMCL.vxi = row_vector(FMCL.vdim);		// Contour level coordinates
 FMCL.vxf = row_vector(FMCL.vdim);		// Contour level coordinates
 FMCL.FMID = 0;					// Zero contour ID number
 FMCL.penneg = 1;				// Number of negative contours
 FMCL.penpos = 2;				// Number of positive contours
 FMCL.conneg = 0;				// Number of negative contours
 FMCL.conpos = 0;				// Number of positive contours
 FMCL.posID = 101;				// Base ID for positive contours
 FMCL.negID = 99;				// Base ID for negative contours

// -------------------- Search Begins Here For Contours ------------------- //

 FMCL.CLmaxnum = 10000;				// Maximum allowed contour lines
 int NCLs = 0;					// Number of contour lines
 row_vector vxnull;				// Empry row vector
  double *Ts;
  int *IDs;
  Ts = new double[2*FMCP.steps];			// Array of contour thresholds
  IDs = new int[2*FMCP.steps];				// Array of contour IDs
 contour_levels(FMCP, FMCL, Ts, IDs);		// Set contour levels
 FMCL.level = 0;				// Contour level number
 for(int cl=0; cl<FMCP.steps; cl++)		// Loop over all contours
   {
   if(cl == FMCL.conpos) FMCL.level=0;		// Used in debugging
   FMCL.threshold = Ts[cl];			// Set the threshold
   FMCL.FMID = IDs[cl];				// Set the contour ID
   contour_pts(mx, FMCP, FMCL);			// Find contour pts, this level
   FMCL.CLmaxsize = FMCL.npts+3;		// Set max. size of a contour line
   row_vector *RCLs, *LCLs;
   RCLs = new row_vector[FMCL.CLmaxnum];		// Array of right contour lines
   LCLs = new row_vector[FMCL.CLmaxnum];		// Array of left contour lines
   for(int k=0; k<FMCL.CLmaxnum; k++)		// Set all vectors to NULL
     {
     RCLs[k] = vxnull;
     LCLs[k] = vxnull;
     }
   int *NRpts, *NLpts;
   NRpts = new int[FMCL.CLmaxnum];			// Number pts in right contour lines
   NLpts = new int[FMCL.CLmaxnum];			// Number pts in left contour lines
   NCLs = 0;					// Number of contour lines
   contour_lines(FMCP,FMCL,RCLs,LCLs,		// Make crude contour lines 
                              NRpts,NLpts,NCLs);
   contour_refine(FMCP,FMCL,RCLs,LCLs,		// Join crude contour lines 
                              NRpts,NLpts,NCLs);
   contour_output(out, FMCP,FMCL,RCLs,LCLs,	// Output full contour lines 
                              NRpts,NLpts,NCLs);
   delete [] LCLs;
   delete [] RCLs;
   delete [] NRpts;
   delete [] NLpts;
     FMCL.level = FMCL.level + 1;		// Increment contour count
   }						// +/- contour looping
 group_contours(out, FMCL.conpos, FMCL.conneg);	// Group all contours


// ---------------------------- Draw in the Axes -------------------------- //

 double hs = FMCP.hsize;			// Horizontal plot size (cm)
 double vs = FMCP.vsize;			// Vertical plot size (cm)
 double hl = FMCP.hlow;				// Horizontal low point (cm)
 double hh = FMCP.hhigh;			// Horizontal high point (cm)
 double vl = FMCP.vlow;				// Vertical low point (cm)
 double vh = FMCP.vhigh;			// Vertical high point (cm)
 FM_Box(out, 0, 0, hs, vs);			// Box around contour plot
 FM_Axis(out, 'x', vs, 0, hs, hl, hh);		// Output horizontal axis (bottom)
 FM_Axis(out, 'y', 0, hs, 0, vl, vh);		// Output vertical axis (left)

// ------------------------ End FrameMaker Document ----------------------- //

 FM_AFrame_End(out);				// End of frame
 FM_AFrames_End(out);				// Anchored frame end
 FM_ParaText_End(out);				// Text flow
 FM_End(out);					// FrameMaker end comment
  delete Ts;
  delete [] IDs;
 }

void FM_contour(const std::string& filename, const matrix &mx, 
                                                  FMcont& FMCP, double* levels)

	// Input	filename  : Output filename
	// 		mx        : Data matrix
	//		FMCP      : Contour plot parameters
	//		levels    : Desired contour levels
	// Return		  : Void, file out is modified
	// Note			  : Uses GAMMA calculated levels
 
 
 { 
 FMclev FMCL;					// Values for a contour level
 contour_setup(mx, FMCP, FMCL); 		// Setup/check parameters 
 double* Ts;					// For contour thresholds
 int* IDs;					// For contour FM ID #s
 Ts = new double[2*FMCP.steps];			// Allocate array of thresholds
 IDs = new int[2*FMCP.steps];			// Allocate array of contour IDs
 contour_levels(FMCP, FMCL, levels, Ts, IDs);	// Set contour levels
 FM_contour(filename, mx, FMCP, FMCL, Ts, IDs); // Perform contouring
 delete [] Ts;					// Delete any thresholds array
 delete [] IDs;					// Deleta any FM ID array
 }

void FM_contour(const std::string& filename, const matrix &mx, FMcont& FMCP)

	// Input	filename  : Output filename
	// 		mx        : Data matrix
	//		FMCP      : Contour plot parameters
	// Return		  : Void, file out is modified
	// Note			  : Uses GAMMA calculated levels
 
 
 { 
 FMclev FMCL;					// Contour level parameters
 contour_setup(mx, FMCP, FMCL); 		// Setup/check parameters 
 double* Ts;					// For contour thresholds
 int* IDs;					// For contour FM ID #s
 Ts = new double[2*FMCP.steps];			// Allocate array of thresholds
 IDs = new int[2*FMCP.steps];			// Allocate array of contour IDs
 contour_levels(FMCP, FMCL, Ts, IDs);		// Set contour levels
 FM_contour(filename, mx, FMCP, FMCL, Ts, IDs); // Perform contouring
 delete [] Ts;					// Delete any thresholds array
 delete [] IDs;					// Deleta any FM ID array
 }


void FM_contour(const std::string& filename, const matrix &mx,
                           FMcont& FMCP, FMclev& FMCL, double* Ts, int* IDs)
	// Input	filename  : Output filename
	// 		mx        : Data matrix
	//		FMCP      : Contour plot parameters
	// 		FMCL	  : Contour level data
	//		Ts        : Contour levels
	//		IDs	  : Contour level IDs (FM)
	// Return		  : Void, file out is modified
	// Note			  : AUXILIARY FUNCTION FOR FM_contour
 
 
 { 
// ----------------- Open File, Set Up FrameMaker Document ---------------- //

 std::ofstream out(filename.c_str());			// Define output stream
 FM_Begin(out);					// FrameMaker begin comment
 FM_AFrames_Begin(out);				// Anchored Frame begin
 FM_AFrame_Set(out, FMCP.hsize, FMCP.vsize);	// Set Anchored Frame

// -------------------- Search Begins Here For Contours ------------------- //

 int NCLs = 0;					// Number of contour lines
 row_vector vxnull;				// Empry row vector
 FMCL.level = 0;				// Contour level number
 for(int cl=0; cl<FMCP.steps; cl++)		// Loop over all contours
   {
   if(cl == FMCL.conpos) FMCL.level=0;		// Used in debugging
   FMCL.threshold = Ts[cl];			// Set the threshold
   FMCL.FMID = IDs[cl];				// Set the contour ID
   contour_pts(mx, FMCP, FMCL);			// Find contour pts, this level
   FMCL.CLmaxsize = FMCL.npts+3;		// Set max. size of contour line
   row_vector *RCLs, *LCLs;
   RCLs = new row_vector[FMCL.CLmaxnum];		// Array of right contour lines
   LCLs = new row_vector[FMCL.CLmaxnum];		// Array of left contour lines
   for(int k=0; k<FMCL.CLmaxnum; k++)		// Set all vectors to NULL
     {
     RCLs[k] = vxnull;
     LCLs[k] = vxnull;
     }
   int *NRpts, *NLpts;
   NRpts = new int[FMCL.CLmaxnum];			// # pts in right contour lines
   NLpts = new int[FMCL.CLmaxnum];			// # pts in left contour lines
   NCLs = 0;					// Number of contour lines
   contour_lines(FMCP,FMCL,RCLs,LCLs,		// Make crude contour lines 
                              NRpts,NLpts,NCLs);
   contour_refine(FMCP,FMCL,RCLs,LCLs,		// Join crude contour lines 
                              NRpts,NLpts,NCLs);
   contour_output(out, FMCP,FMCL,RCLs,LCLs,	// Output full contour lines 
                              NRpts,NLpts,NCLs);
   delete [] RCLs;
   delete [] LCLs;
   delete [] NRpts;
   delete [] NLpts;
   FMCL.level = FMCL.level + 1;			// Increment contour count
   }						// +/- contour looping
 group_contours(out, FMCL.conpos, FMCL.conneg);	// Group all contours


// ---------------------------- Draw in the Axes -------------------------- //

 double hs = FMCP.hsize;			// Horizontal plot size (cm)
 double vs = FMCP.vsize;			// Vertical plot size (cm)
 double hl = FMCP.hlow;				// Horizontal low point (cm)
 double hh = FMCP.hhigh;			// Horizontal high point (cm)
 double vl = FMCP.vlow;				// Vertical low point (cm)
 double vh = FMCP.vhigh;			// Vertical high point (cm)
 FM_Box(out, 0, 0, hs, vs);			// Box around contour plot
 FM_Axis(out, 'x', vs, 0, hs, hl, hh);		// Output horiz. axis (bottom)
 FM_Axis(out, 'y', 0, hs, 0, vl, vh);		// Output vertical axis (left)

// ------------------------ End FrameMaker Document ----------------------- //

 FM_AFrame_End(out);				// End of frame
 FM_AFrames_End(out);				// Anchored frame end
 FM_ParaText_End(out);				// Text flow
 FM_End(out);					// FrameMaker end comment
 }

// ____________________________________________________________________________
//                           FRAMEMAKER MMF FUNCTIONS
// ____________________________________________________________________________

	// Input	out       : output stream
	// Output	none      : Function is void.  Puts out the
	//			    header for a mathematical expressions

void MMF_Head (std::ostream &out)
  {
  out << "  <Math                                              \n";
  out << "  <BRect  2.607 \" 0.168 \" 1.119 \" 0.662 \">       \n";
  out << "  <MathFullForm "                                     ;
  }


void MMF_Tail (std::ostream &out)

	// Input	out       : output stream
	// Output	none      : Function is void.  Puts out the
	//			    tail for a mathematical expressions

{
   out << "  > # end of MathFullForm                            \n";
   out << "  <MathOrigin  3.167 \" 0.535 \">                    \n";
   out << "  <MathAlignment Center >                            \n";
   out << "  <MathSize MathMedium >                             \n";
   out << "  <Angle 0>                                          \n";
   out << "  > # end of Math                                    \n";
}


void FM_PtFrame_Set (std::ostream &out,
                 double xsize, double ysize, int ID=11)
    // Set Anchored Frame of a MIF File
{
 out << "  <Frame                                       \n";
 out << Gform("  <ID %d>\n", ID);
 out << "  <Pen 15>                                     \n";
 out << "  <Fill 7>                                     \n";
 out << "  <PenWidth  0.5 pt>                           \n";
 out << "  <Separation 0>                               \n";
  out << Gform("  <BRect 0.0 pt 0.0 pt %3.0f pt ",xsize);
  out << Gform("%3.0f pt>\n", ysize);
 out << "  <FrameType Below >                           \n";
 out << "  <Float Yes>                                  \n";
 out << "  <NSOffset  0.0 \">                           \n";
 out << "  <AnchorAlign Center >                        \n";
 out << "  <Cropped No >                                \n";
 out << "  <Pen 0>                                      \n";
}

void FM_Matrix_Real_Number(std::ostream &out, double num,
		 	int prec = FM_PP, double threshold = FM_PT)

	// Input	out       : output stream
	//		num	  : real number
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out a
	//			    real number in a mathematical expression
	// Note			  : Sign is included in output

{
	 // Changed format[5] to format[51].
   char format[51];
   if (fabs(int(num) - num) > threshold) 
		{
#ifdef _MSC_VER
      sprintf_s(format, 30, "%%.%df", prec); 
#else
			sprintf(format,"%%.%df",prec);
#endif
		}
   else
		{
#ifdef _MSC_VER
      sprintf_s(format, 30, "%%.0f"); 
#else    
      sprintf(format,"%%.0f"); 
#endif
		}
   out << "num" << LBrak << Gform(format,num) << ",\"" << Gform(format,num) << "\"]";
}


void FM_Matrix_Imag_Number(std::ostream &out, double num,
			 int prec = FM_PP, double threshold = FM_PT)

	// Input	out       : output stream
	//		num	  : imaginary (real) number
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out an
	//			    imaginary number in a mathematical expression
	// Note			  : Sign is included in output

{
   if (fabs(num) > threshold) {
     if (num>0) {
         out << "cdot" << LBrak;
         FM_Matrix_Real_Number ( out, num, prec ,threshold);
         out << ",char" << LBrak << "i]]";
     }
     else {
         out << "minus" << LBrak << "cdot" << LBrak;
         FM_Matrix_Real_Number( out, -num, prec ,threshold);
         out << ",char" << LBrak << "i]]]";
     }
   }
   else 
     FM_Matrix_Imag_Number ( out, num, prec ,threshold);
}


void FM_Matrix_Complex_Number (std::ostream &out, complex num,
			 int prec = FM_PP, double threshold = FM_PT)

	// Input	out       : output stream
	//		num	  : complex number
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out a
	//			    complex number in a mathematical expression
	// Note			  : Sign is included in output

{
   double rep, imp;
   rep = Re(num);
   imp = Im(num);

   // Within threshold equal to zero. Put out a simple '0'
   if( fabs(rep) < threshold && fabs(imp) < threshold) 
      FM_Matrix_Real_Number ( out, fabs(rep), prec ,threshold);

   // Within threshold a purely real number 
   if( fabs(rep) > threshold && fabs(imp) < threshold)
      FM_Matrix_Real_Number ( out, rep, prec ,threshold);

   // Within threshold a purely immaginary number
   if( fabs(rep) < threshold && fabs(imp) > threshold)
      FM_Matrix_Imag_Number ( out, imp, prec ,threshold);
   
   // A full complex number
   if( fabs(rep) > threshold && fabs(imp) > threshold) {
         out << "id" << LBrak << "plus" << LBrak;
         FM_Matrix_Real_Number ( out, rep, prec ,threshold);
         out << "," ;
         FM_Matrix_Imag_Number ( out, imp, prec ,threshold);
         out << "]]";
    }
      
}

/* Note that the FM_Matrix & FM_Mat_Plot functions which take GAMMA operators
   and superoperators have been taken out of this module.  Use of these same
   functions involving matrices can be used by simply asking gen_op & super_op
   for their matrix then calling the functions herein.  I don't know if I'll
   put the ones with Op & LOp back into GAMMA in a higher level module.  The
   reason they were removed is to "decouple" this module for the Hilbert 
   space and Liouville space modules.                                        */


/*
void FM_Matrix(const std::string& filename, const gen_op& Op, int prec,  double threshold)

	// Input	out       : output stream
	//		Op	  : general operator
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    of an operator to a mathematical expression

  { FM_Matrix(filename, Op.get_mx(), prec, threshold); }


void FM_Matrix (const std::string& filename, const super_op& LOp, int prec,  double threshold)

	// Input	out       : output stream
	//		LOp	  : superoperator
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    of a superoperator to a mathematical expression

  { FM_Matrix (filename, LOp.get_mx(), prec,threshold); }
*/


void FM_Matrix(const std::string& filename, const matrix& mx, int prec, double threshold)

	// Input	out       : output stream
	//		mx	  : matrix
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    to a mathematical expression

  {
  std::ofstream out(filename.c_str());	//Output file
  FM_Begin(out);			//Initialize MIF
  MMF_Head(out);	  		//Initialize MMF

  out << "`matrix" << LBrak;                  // Start a FM matrix

  int nr = mx.rows();			// Number of matrix rows 
  int nc = mx.cols();			// Number of matrix columns
  out << nr << "," << nc ;
  for(int i=0; i<nr; i++ ) 		// Loop over matrix elements
    {
    for(int j=0; j<nc; j++ )
      {
      out << ","  << "\n";
      FM_Matrix_Complex_Number(out, mx.get(i,j), prec, threshold);
      }
    }
  out << "]'\n";
  MMF_Tail(out);			// Terminate MMF
  FM_End(out);				// Terminate MIF
  }


/*
void FM_Mat_Plot(const std::string& filename, const gen_op& Op, double threshold)

	// Input	out       : output stream
	//		Op	  : general operator
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    of an operator as a schematic
	//			    graphics object

  { FM_Mat_Plot(filename, Op.get_mx(), threshold); }


void FM_Mat_Plot(const std::string& filename, const super_op& LOp, double threshold)

	// Input	out       : output stream
	//		LOp	  : superoperator
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    of a superoperator as a schematic
	//			    graphics object

  { FM_Mat_Plot (filename, LOp.get_mx(), threshold); }


void FM_Mat_Plot(const std::string& filename, const gen_op& Op,
		 		      const gen_op& Op_ref, double threshold)

	// Input	out       : output stream
	//		Op	  : general operator
	//		Op_ref	  : reference general operator
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    representing the differences between
	//			     Op and Op_ref.  The output is as a schematic
	//			    graphics object

  { FM_Mat_Plot (filename, Op.get_mx(), Op_ref.get_mx(), threshold); }


void FM_Mat_Plot(const std::string& filename, const super_op& LOp,
		 		const super_op& LOp_ref, double threshold)

	// Input	out       : output stream
	//		LOp	  : superoperator
	//		LOp_ref	  : reference superoperator
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    representing the differences between
	//			    LOp and LOp_ref.  The output is as a schematic
	//			    graphics object

  { FM_Mat_Plot (filename, LOp.get_mx(), LOp_ref.get_mx(), threshold); }
*/


void FM_Mat_Plot(const std::string& filename, const matrix &mx, double threshold)

	// Input	out       : output stream
	//		mx	  : matrix
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    as a schematic graphics object

  {
   std::ofstream out(filename.c_str());	//Output file
   const float Rand = 5;
   const float Startx = Rand;	//Positioning in FrameMaker (upper left x coord)
   const float Starty = Rand;	//Positioning in FrameMaker (upper left y coord)
   const float Sizex = 512;	//x-Size of Matrix Picture
   const float Sizey = 512;	//y-Size of Matrix Picture

   FM_Begin(out);		//Title Line for MIF File
   FM_AFrames_Begin(out);                // Anchored Frame begin
   FM_PtFrame_Set(out,Sizex+2*Rand,Sizey+2*Rand);       // Header of Frame
					
// Number of rows and columns
   int nr,nc;
   float incx,incy;
   nr=mx.rows();                         
   nc=mx.cols();
   incx=Sizex/nc;
   incy=Sizey/nr;
   std::string Fill;
   // Make the empty grid
   int i;
   for(i=0; i<nr+1; i++ )
     {
     out << "<PolyLine\n"; 
     out << "<GroupID 2831>\n";	    // Group grid-lines together
     out << "<NumPoints 2>\n";
     out << "<Point " <<  Startx+i*incx << "pt " << Starty <<  "pt>\n";
     out << "<Point " <<  Startx+i*incx << "pt " << Starty+Sizey <<  "pt>\n";     
     out << "> # end of Polyline\n";    
     }
   for ( i = 0; i < nc+1; i++ ){
         out << "<PolyLine\n"; 
         out << "<NumPoints 2>\n";
         out << "<GroupID 2831>\n";	    // Group grid-lines together
         out << "<Point " <<  Startx << "pt " << Starty+i*incy <<  "pt>\n";
         out << "<Point " <<  Startx+Sizex << "pt " << Starty+i*incy <<  "pt>\n";     
         out << "> # end of Polyline\n";    
   }
  out << "<Group\n";
  out << "<ID 2831>\n";
  out << "<GroupID 2835>\n";     //Group the whole object together
  out << "> # end of Group\n";

   // and fill in the necessary squares
   for ( i = 0; i < nr; i++ ){
       for ( int j = 0; j < nc; j++ ) {
	// if the element is nonzero, make a black rectangle        Fill = "<Fill 7>";
        if(norm(mx.get(i,j)) > threshold) {
           Fill = "<Fill 0>";
           out << "<Rectangle\n" ;             // Start a rectangle for framing"
           out << "<GroupID 2833>\n";	    // Group rectangles together
           out << " <Pen 0>\n";	            // Black frame
	   out << Fill << "\n";
	   out << "<BRect " << Startx+j*incx<< " pt " << Starty+i*incy << " pt ";
           out <<             incx  << " pt " << incy << " pt >\n";
           out << "> # end of Rectangle\n";
        }
   }
  }
     out << "<Group\n";
     out << "<ID 2833>\n";
     out << "<GroupID 2835>\n";
     out << "> # end of Group\n";

     out << "<Group\n";
     out << "<ID 2835>\n";
     out << "> # end of Group\n";


    FM_AFrame_End (out);                  // End of Frame
    FM_AFrames_End(out);                  // Anchored Frame end
    FM_ParaText_End (out);                // Create a TextFlow with Frame
   FM_End(out);			// Terminate MIF
}



void FM_Mat_Plot (const std::string& filename, const matrix &mx,
				 const matrix &mx_ref, double threshold)

	// Input	out       : output stream
	//		mx	  : matrix
	//		mx_ref	  : reference general operator
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    representing the differences between Op
	//			    and Op_ref.  The output is as a schematic
	//			    graphics object

  {
  std::ofstream out(filename.c_str());	//Output file
  const float Rand = 5;
  const float Startx = Rand;	//Positioning in FM (upper left x coord)
  const float Starty = Rand;	//Positioning in FM (upper left y coord)
  const float Sizex = 512;	//x-Size of Matrix Picture
  const float Sizey = 512;	//y-Size of Matrix Picture

  FM_Begin(out);		//Title Line for MIF File
  FM_AFrames_Begin(out);                // Anchored Frame begin
  FM_PtFrame_Set(out,Sizex+2*Rand,Sizey+2*Rand);       // Header of Frame
					
// Number of rows and columns
   int nr,nc;
   float incx,incy;
   nr=mx.rows();                         
   nc=mx.cols();
   incx=Sizex/nc;
   incy=Sizey/nr;
   std::string Fill;

  int i;
  for(i=0; i<nr+1; i++)			// First make an empty grid
    {					// Start by drawing horizontals
    out << "<PolyLine\n"; 		// Each a polyline
    out << "<NumPoints 2>\n";		// containing two poings
    out << "<GroupID 2831>\n";		// Group ID for horizontals
    out << "<Point " <<  Startx+i*incx << "pt " << Starty <<  "pt>\n";
    out << "<Point " <<  Startx+i*incx << "pt " << Starty+Sizey <<  "pt>\n";     
    out << "> # end of Polyline\n";    
    }
  for(i=0; i<nc+1; i++)			// Finish grid by drawing veritcals
    {
    out << "<PolyLine\n"; 		// Each vertical a polyline
    out << "<NumPoints 2>\n";		// containing two poings
    out << "<GroupID 2831>\n";		// Group ID for verticals
    out << "<Point " <<  Startx << "pt " << Starty+i*incy <<  "pt>\n";
    out << "<Point " <<  Startx+Sizex << "pt " << Starty+i*incy <<  "pt>\n";     
    out << "> # end of Polyline\n";    
    }
  out << "<Group\n";
  out << "<ID 2831>\n";
  out << "<GroupID 2835>\n";			// Group entire grid together
  out << "> # end of Group\n";

  for(i=0; i<nr; i++)				// Loop matrix elements &
    { 						// fill in proper squares
    for(int j=0; j<nc; j++)			// Elements different from mx_ref
      {						// by "threshold" are black (<Fill 7>)
      if(norm(mx.get(i,j)-mx_ref.get(i,j))
              > threshold*norm(mx_ref.get(i,j))) 
        {
        Fill = "<Fill 0>";			// Set for dark fill
        out << "<Rectangle\n";			// Start a rectangle
        out << "<GroupID 2833>\n";		// Group id for rectangles
        out << " <Pen 0>\n";			// Black rectangle frame
        out << Fill << "\n";			// Output fill command
        out << "<BRect " << Startx+j*incx<< " pt " << Starty+i*incy << " pt ";
        out <<             incx  << " pt " << incy << " pt >\n";
        out << "> # end of Rectangle\n";
        }
      }
    }
  out << "<Group\n";
  out << "<ID 2833>\n";
  out << "<GroupID 2835>\n";
  out << "> # end of Group\n";

  out << "<Group\n";
  out << "<ID 2835>\n";
  out << "> # end of Group\n";

  FM_AFrame_End (out);                  // End of Frame
  FM_AFrames_End(out);                  // Anchored Frame end
  FM_ParaText_End (out);                // Create a TextFlow with Frame
  FM_End(out);				// Terminate MIF
  }


/*****************************************************************
**								**
** 		     FrameMaker Molecule Plots			**
**								**
** Author	: SOSI 						**
** Description  : Takes a molecule from input along with plot	**
**		  parameters and produces a FrameMaker file	**
**		  containing a plot of the molecule		**
**								**
*****************************************************************/


/* sosi
void FM_molecule (const std::string& filename, molecule mol, double radius,
	 double alpha, double beta, double gamma)
//  int type=2, int plpts=100

	// Input	filename  : Output filename
	//              mol       : Molecule
	//		radius    : Overall plot dimension in cm
	// 		alpha     : Euler angle (radians)
	// 		beta      : Euler angle (radians)
	// 		gamma     : Euler angle (radians)
	// Return		  : Void. FrameMaker MIF molecule plot
	//			    is output to the specified file

	//		type      : Type of plotting to use
	//				0 = discrete points
	//				1 = joined points
	//				2 = vectors from origin
	//		plpts     : Number of points to use when
	//			    drawing plane-sphere intersection

{
  if((radius < 1) || (radius > 10))
			   radius = 5;	// Insure plot size sensible
					// Begin FrameMaker Plot
  std::ofstream out(filename); 		// Open new file for plotting
  FM_3D_start (out, radius);		// Initialize FM file
  FM_sphere_axes(out, alpha, beta,	// Draw in three axes
	                 gamma, radius);
  row_vector vx = FM_sphere_proj(mol,	// Get atomic coordinates in plane
	      alpha, beta, gamma, radius);
  int molID = 72;
  double x, y;
  int atoms = ((coord_vec)mol).size();

  int color[atoms];			// Assign colors to atoms
  int acol = 7;
  color[0] = acol;
  int k = 0;
  int hit = 0;
  for (int j=1; j<atoms; j++)
    {
    color[j] = -1;
    k = 0;
    hit = 0;
    while((k < j) && (hit == 0))
      {
      if(mol.isotope(k) == mol.isotope(j))
        {
        color[j] = color[k];
        hit = 1;
        }
      k++;
      }
    if(color[j] == -1)			 // New atom type
      {
      acol--;
      if(acol < 0) acol = 7;
      color[j] = acol;
      }
    }
  for (int i=0; i<atoms; i++)
    {
    x = Re(vx(i));
    y = Im(vx(i));
//    x = (Re(vx(i))-xmin)*xscale+left;
//    y = (Im(vx(i))-ymin)*yscale+bottom;
    FM_Ellipse (out, x, y, radius/10.0, radius/10.0, molID, 7, color[i]);
    FM_TextLine(out, molID, x, y, mol.symbol(i));
    }
  FM_Group(out, molID);			// Group the molecule together
  FM_3D_end (out);			// Finalize FM file
}
sosi */

/*****************************************************************
**								**
** 		 FrameMaker Table Generator 			**
**								**
** Author	: SOSI 						**
** Description  : Takes a matrix from input and attempts to	**
**		  produce a FrameMaker file a table of the	**
**		  matrix elements.				**
**								**
*****************************************************************/

void FM_Mat_Tbl(const std::string& filename, const matrix& mx, double threshold)

	// Input	out       : output stream
	//		mx	  : matrix
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    as a schematic graphics object

  {
  int mxr = mx.rows(); 				// Get matrix rows
  int mxc = mx.cols(); 				// Get matrix columns
  std::ofstream out(filename.c_str());		// Begin output file
  FM_Begin(out);				// Title Line for MIF File
  FM_Tbl_Begin(out);				// Begin Table, # 1
  out << " <TblTag `Format A'>\n";		// Specify Table Format A
  out << Gform(" <TblNumColumns %d>\n", mxc);	// Set Table Columns
  FM_Tbl_Title(out);				// Write Table Title
  FM_TblBody_Begin(out);			// Begin Table Body
  double x;
//int j=0;
  for(int i=0; i<mxr; i++)
    {
    out << "  <Row \n";				// Start Table Row
    for(int j=0; j<mxc; j++)
      {
      x = mx.getRe(i,j);			// Get matrix value
      out << "   <Cell \n";			// Start Table Cell
      out << "    <CellContent \n";		// Start Cell Contents
      out << "     <Para \n";			// Cell Paragraph
      out << "      <PgfTag `CellHeading'> \n";	// Cell Paragraph Type
      out << "      <ParaLine \n";
      out << "      <String `" << x << "'> \n";	// Output value 
      out << "      > \n";			//
      out << "     > # end of Para\n";		// End fo Cell Paragraph
      out << "    > # end of CellContent\n";	// End Cell Contents
      out << "   > # end of Cell \n";		// End of Table Cell
      }
    out << "  > # end of Row\n";		// End Table Row
    }
  FM_TblBody_End(out);				// End Table Body
  FM_Tbl_End(out);				// End Table, # 1
  FM_TextFlow_Set(out);				// Begin Text Flow
  FM_Paragraph_Set(out);			// Begin Paragraph
  out << "  <ATbl 1>\n";			// Reference Table #1
  FM_Paragraph_End(out);			// End Paragraph
  FM_TextFlow_End(out);				// End of Text Flow
  FM_End(out);					// Terminate MIF
 threshold = 0;					// Compiler likes this used
 }


#endif 						// FrameMaker.cc
