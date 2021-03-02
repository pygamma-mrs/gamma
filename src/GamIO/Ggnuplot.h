/* Ggnuplot.h ***************************************************-*-c++-*-
**                                                                      **
**                                 G A M M A                            **
**                                                                      **
**      Gnuplot Support                            Interface		**
**                                                                      **
**      Scott A. Smith							**
**      Copyright (c) 1994, 1998                                        **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  This is a still rather crude attempt to create an interface between **
**  GAMMA and the public domain package called gnuplot.  Even though    **
**  the functions herein work, they are far from what one would call    **
**  optimal!  As things progress they should work better, and perhaps   **
**  this will be made into a class which supports gnuplot command files **
**  some day.  For now, it is stupid and just puts out data files.      **
**                                                                      **
*************************************************************************/

#ifndef   Ggnuplot_h_			// Is this file already included?
#  define Ggnuplot_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/row_vector.h>		// Know GAMMA row vectors
#include <Matrix/matrix.h>		// Know GAMMA matrices
#include <Level1/coord_vec.h>		// Know GAMMA coordinate vectors
#include <string>                       // Know libstdc++ strings
#include <vector>                       // Know libstdc++ STL vectors
#include <list>                         // Know libstdc++ STL lists

struct GPdat				// Gnuplot info
  {
  std::string term;			// Terminal type
  std::string basedir;			// Base directory
  std::string macro;			// Macro local filename
  std::string data;			// ASCII data local filename
  std::string outfile;			// Output filename
  std::string basename;			// Base filename
  std::string cmd;			// System command 
  std::string ylabel;			// Label for yaxis
  std::string xlabel;			// Label for xaxis
  std::string title;			// Label for plot
  std::string tics;			// Tics in vs out
 
  int nokey;				// Flag for key (yes/no)
  int border;				// Border flag (yes/no)
  int join;				// Line draw flag (yes/no)
  int xtics;				// X-axis tic flag (yes/no)
  int ytics;				// Y-axis tic flag (yes/no)
  int xaxis;				// Flag for x-axis draw (xes/no)
  int yaxis;				// Flag for y-axis draw (yes/no)

  int    compress;			// XY plot data compression

  double xmin;				// XY plot x minimum value
  double xmax;				// XY plot x maximum value
  double cutoff;			// XY plot intensity resolution
  double riflag;			// XY plot {real,imag,complex} flag
  double xsize;				// X-axis scaling factor
  double ysize;				// Y-axis scaling factor

//  bool   pause;				// Flag for pause after plot

  GPdat()
    {
    nokey    = 0;			// Keep key plotted
    border   = 1;			// Plot border
    join     = 1;			// Join plotted point
    compress = 0;			// No data compression
    xaxis    = 1;			// Plot x axis
    yaxis    = 1;			// Plot y axis

    xmin     = 0;			// No minimum on x axis
    xmax     = 0;			// No maximum on x axis
    cutoff   = 0;			// No intensity difference cutoff
    riflag   = 0;			// Plot reals only
    xsize    = 0;			// No plot x size set
    ysize    = 0;			// No plot y size set
    } 
  };


// ____________________________________________________________________________
// i                      GAMMA Gnuplot Error Handling
// ____________________________________________________________________________
 
/*         Input                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */

         void GPerror(int    eidx,                              int noret=0);
         void GPerror(int    eidx,    const std::string& pname, int noret=0);
volatile void GPfatality(int eidx);
volatile void GPfatality(int eidx,    const std::string& pname);

//_____________________________________________________________________________
// A                     GAMMA Gnuplot System Interface
//_____________________________________________________________________________

        // Input        none    : Output filename
        //              warn    : Warning flag
        // Return       GPE     : A string for the Gnuplot
        //                        executable command (hopefully)

MSVCDLL std::vector<std::string> GPSeekStrings();	// Strings for gnuplot path
MSVCDLL std::string GPFind(bool  vocal=false);		// Look for gnuplot exec.
MSVCDLL std::string GPExec(int   warn=0);		// Gnuplot executable name
MSVCDLL void RunGnuplot(const std::string& gnumacro);	// Run Gnuplot

MSVCDLL void SetLineType(std::ostream& gnuload, int ltype);
MSVCDLL void CloseMacro(std::ostream&  gnuload);

//_____________________________________________________________________________
// B                      GAMMA Gnuplot 1D Output Files
//_____________________________________________________________________________

//_____________________________________________________________________________
//                          gnuplot Parametric Plots
//_____________________________________________________________________________

 
MSVCDLL void GP_xy(const std::string& filename, const row_vector &vx);
 
        // Input        filename  : Output filename
        //              vx        : Data vector
        // Return                 : Void, file out is modified

 
MSVCDLL void GP_xy(std::ofstream& ofstr, const row_vector& vx);
 
        // Input        ofstr     : An output filestream
        //              vx        : Data vector
        // Return                 : Void, output filestream is modified


//____________________________________________________________________________
//                           Gnuplot Stack Plots 
//____________________________________________________________________________


MSVCDLL void GP_stack(std::ofstream& ofstr, const row_vector &vx, double row,
              int ri=0, double xmin=0, double xmax=0, double cutoff=0);

        // Input        ofstr     : Output filestream
        //              vx        : Data vector
        //              row       : Row value (y value)
        //              ri        : Flag for real versus imaginary plot
        //                          ri = 0  : reals only
        //                          ri = <0 : imaginaries only
        //                          ri = >0 : both 
        //              xmin      : Plot horizontal (x=0) point value
        //              xmax      : Plot horizontal (x=end) point value
        //              cutoff    : Data compression cutoff
        // Return                 : Void, the information in vector vx is
        //                          output to ofstr in ASCII format suitable
        //                          for use in gnuplot.  Each data point is
        //                          written as {x,y,z} where x is the point
        //                          number, y the row number, and z the vector
        //                          point value 
        // Note                   : This is essentially the same function as
        //                          GP_1D except that output is {x,y,z} not {x,y}
 

MSVCDLL void GP_stack(const std::string& filename, matrix& mx, int ri=0,
                 double ymin=0, double ymax=0, double xmin=0, double xmax=0);

        // Input        filename  : Output file name
        //              mx        : Data matrix
        //              ri        : Flag for real versus imaginary plot
        //                          ri = 0  : reals only
        //                          ri = <0 : imaginaries only
        //                          ri = >0 : both 
        //              ymin      : Plot "vertical" (y=0) point value
        //              ymax      : Plot "vertical" (y=end) point value
        //              xmin      : Plot horizontal (x=0) point value
        //              xmax      : Plot horizontal (x=end) point value
        // Return                 : Void, the information in vector mx is
        //                          output to file filename in an ASCII
        //                          format suitable for use in gnuplot.
        //                          Each data point is written as {x,y,z}
        //                          where x is the point number, y the row number,
        //                          and z the vector point value 
 
 
MSVCDLL void GP_stack(std::ofstream& ofstr, matrix& mx, int ri=0,
                 double ymin=0, double ymax=0, double xmin=0, double xmax=0);
 
        // Input        ofstr     : An output file stream
        //              mx        : Data matrix
        //              ri        : Flag for real versus imaginary plot
        //                          ri = 0  : reals only
        //                          ri = <0 : imaginaries only
        //                          ri = >0 : both 
        //              ymin      : Plot "vertical" (y=0) point value
        //              ymax      : Plot "vertical" (y=end) point value
        //              xmin      : Plot horizontal (x=0) point value
        //              xmax      : Plot horizontal (x=end) point value
        // Return                 : Void, the information in vector mx is
        //                          output to file filename in an ASCII
        //                          format suitable for use in gnuplot.
        //                          Each data point is written as {x,y,z}



MSVCDLL void GnuplotStack(const std::string& filename, const std::vector<row_vector>& vxs,
                 int ri=0, double ymin=0, double ymax=0, double xmin=0, double xmax=0);
 
//____________________________________________________________________________
//                            gnuplot Contour Plots 
//____________________________________________________________________________
 
/*****************************************************************
**                      gnuplot Contour Plots                   **
**                                                              **
** Author       : S.A. Smith                                    **
** Date         : October 4, 1994                               **
** Description  : These routines are intended to output data    **
**                files in ASCII format which are readable by   **
**                the public domain program gnuplot.  Thier     **
**                purpose is to provide an easy and accessible  **
**                means to view GAMMA generated contour plots.  **
** Notes        : gnuplot takes ASCII files filled with {x,y,z} **
**                points as input for 3D plots.  Only one point **
**                per line is allowed, the three ordinates      **
**                separated by spaces. The x & y ordinates are  **
**                optional but are output in these routines.    **
**                For contour plots, the output ASCII file must **
**                have a blank line following each row.  Also,  **
**                the points must be "gridded", meaning that    **
**                all points of the 2D data set must be present **
**                (no skipping allowed as in the stack plots).  **
**                The plots themselves will be displayed as     **
**                                                              **
**                      [-x, +x] from left to right             **
**                      [-y, +y] from front to back             **
**                      [-z, +z] from top to bottom             **
**                                                              **
**                Useful gnuplot commands for stack plots are   **
**                                                              **
**                      set parametric       (independent rows) **
**                      set data style lines (connected points) **
**                      splot "filename"     (draw stack plot)  **
**                                                              **
*****************************************************************/
 
MSVCDLL void GP_contour(std::ofstream& ofstr, const row_vector &vx, double row,
              int ri=0, double xmin=0, double xmax=0, double cutoff=0);
 
        // Input        ofstr     : Output filestream
        //              vx        : Data vector
        //              row       : Row value or y-axis value
        //              ri        : Flag for real versus imaginary plot
        //                          ri = 0  : reals only
        //                          ri = <0 : imaginaries only
        //                          ri = >0 : both
        //              xmin      : Plot horizontal (x=0) point value
        //              xmax      : Plot horizontal (x=end) point value
        //              cutoff    : Data cutoff (where to set zero)
        // Return                 : Void, the information in vector vx is
        //                          output to ofstr in ASCII format suitable
        //                          for use in gnuplot.  Each data point is
        //                          written as {x,y,z} where x is the point
        //                          number, y the row number, and z the vector
        //                          point value
        // Note                   : This is essentially the same function as
        //                          GP_stack except that all points are output
 

MSVCDLL void GP_contour(const std::string& filename, matrix& mx, int ri=0,
                 double ymin=0, double ymax=0, double xmin=0, double xmax=0);

        // Input        filename  : Output file name
        //              mx        : Data matrix
        //              ri        : Flag for real versus imaginary plot
        //                          ri = 0  : reals only
        //                          ri = <0 : imaginaries only
        //                          ri = >0 : both 
        //              ymin      : Plot "vertical" (y=0) point value
        //              ymax      : Plot "vertical" (y=end) point value
        //              xmin      : Plot horizontal (x=0) point value
        //              xmax      : Plot horizontal (x=end) point value
        // Return                 : Void, the information in vector mx is
        //                          output to file filename in an ASCII
        //                          format suitable for use in gnuplot.
        //                          Each data point is written as {x,y,z}
        //                          where x is the point #, y the row #,
        //                          and z the vector point value 
 
 
MSVCDLL void GP_contour(std::ofstream& ofstr, matrix& mx, int ri=0,
                 double ymin=0, double ymax=0, double xmin=0, double xmax=0);
 
 
        // Input        ofstr     : An output file stream
        //              mx        : Data matrix
        //              ri        : Flag for real versus imaginary plot
        //                          ri = 0  : reals only
        //                          ri = <0 : imaginaries only
        //                          ri = >0 : both 
        //              ymin      : Plot "vertical" (y=0) point value
        //              ymax      : Plot "vertical" (y=end) point value
        //              xmin      : Plot horizontal (x=0) point value
        //              xmax      : Plot horizontal (x=end) point value
        // Return                 : Void, the information in vector mx is
        //                          output to file filename in an ASCII
        //                          format suitable for use in gnuplot.
        //                          Each data point is written as {x,y,z}

 
//____________________________________________________________________________
//                        Gnuplot Spherical Plots
//____________________________________________________________________________
 
 
//void GP_sphere(const std::string& name, const coord_vec& data, bool old=false);
//void GnuplotSphere(const std::string& name, const coord_vec& data);

//____________________________________________________________________________
//                        Gnuplot Array Visualization
//____________________________________________________________________________

MSVCDLL void GP_array(const std::string& name, matrix& mx,
                             int type=1, bool plot=true, double cutoff=1.e-9);
 
//____________________________________________________________________________
//                        Gnuplot Interactive Plots
//____________________________________________________________________________



        // Input        gnumacro  : A filename for gnu macro file
        //              file1D    : A filename of 1D gnuplot data
	//		join      : Flag whether to join points
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output

        // Input        G         : Plotting controls
        // Return       void      : 1D plot performed with gnuplot
        //                          according to information in
        //                          structure G

        // Input        basedir   : Base directory name
        //              basename  : Base file name
        //              term      : Gnuplot terminal type
        //              join      : Flags whether to join points
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output
 


        // Input        gnumacro  : A filename for gnu macro file
        //              file1D    : A filename of 1D gnuplot data
        //              term      : Gnuplot terminal type
        //              ofile     : Gnuplot output file
        //              join      : Flags whether to join points
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output


 
        // Input        gnumacro  : A filename for gnu macro file
        //              N         : Number of plot (ASCII) files
        //              files     : Name of plot (ASCII) files
        //              join      : Flags whether to join points
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output

 
        // Input        gnumacro  : A filename for gnu macro file
        //              N         : Number of plot (ASCII) files
        //              files     : Name of plot (ASCII) files
        //              join      : Flags whether to join points
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output


//--------------------------- Common Parametric Plots -------------------------

MSVCDLL void GP_xyplot(const std::string& gnumacro, const std::string& file1D, int join=1);

        // Input        gnumacro  : A filename for gnu macro file
        //              file1D    : A filename of 1D gnuplot data
        //              join      : Flags whether to join points
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output
 

//------------------------------ Common Stack Plots ---------------------------

 
MSVCDLL void GP_stackplot(const std::string& gnumacro, const std::string& file2D);
 
        // Input        gnumacro  : A filename for gnu macro file
        //              file2D    : A filename of 2D gnuplot data
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output
 

//---------------------------- Common Contour Plots ---------------------------

MSVCDLL void GP_contplot(const std::string& gnumacro, const std::string& file2D, int pf=0);

        // Input        gnumacro  : A filename for gnu macro file
        //              file2D    : A filename of 2D gnuplot data
        //              pf        : Plot flag (default is 0)
        //                              0 - Only make contour plot
        //                             !0 - Plot contours & surface
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output
 

//--------------------------- Common Spherical Plots --------------------------
 
 
//void GP_sphereplot(const std::string& Gname, const std::string& Aname, bool old=false);

        // Input        gnumacro  : A filename for gnu macro file
        //              file2D    : A filename of 2D gnuplot data
        //              pf        : Plot flag (default is 0)
        //                              0 - Only make contour plot
        //                             !0 - Plot contours & surface
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output
 
 
//void GP_sphereplot(const std::string& gnumacro, int N, std::string* files);
 
        // Input        gnumacro  : A filename for gnu macro file 
        //              N         : Number of plot (ASCII) files
        //              files     : Names of plot (ASCII) files
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output
 

//____________________________________________________________________________
//                           Gnuplot Structure Functions
//____________________________________________________________________________


MSVCDLL void defaults(GPdat& G);

        // Input        G         : Plotting controls
        // Return       void      : The plotting controls are
        //                          set for default values


//____________________________________________________________________________
//                           Gnuplot Output Blurbs
//____________________________________________________________________________


MSVCDLL void GP_stackblurb(std::ofstream& ostr, const std::string& plotname);

        // Input        ostr      : An output stream
        //              plotname  : A gnuplot stack plotname
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output
        //                          into the ofstream


MSVCDLL void GP_contblurb(std::ofstream& ostr, const std::string& plotname);

        // Input        ostr      : An output stream
        //              plotname  : A gnuplot contour plotname
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output
        //                          into the ofstream

#endif 						// Ggnuplot.h
