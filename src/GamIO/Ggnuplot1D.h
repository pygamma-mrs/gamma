/* Ggnuplot1D.h *************************************************-*-c++-*-
**                                                                      **
**                                 G A M M A                            **
**                                                                      **
**      Gnuplot Support X Vs. Y Plots                  Interface	**
**                                                                      **
**      Scott A. Smith							**
**      Copyright (c) 2002						**
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  This is the portion of the Gnuplot related functionality which will	**
**  handle production of X versus Y type of plots in GAMMA. There are	**
**  three basic types of functions: 1.) Those which produce ASCII files	**
**  of data for use in Gnuplot plots, 2.) Those which produce Gnuplot	**
**  load files (also ASCII) and then plot ASCII data files using it,	**
**  and 3.) Those which do both, producing a screen plot interactively.	**
**                                                                      **
*************************************************************************/

#ifndef   Ggnuplot1D_h_			// Is this file already included?
#  define Ggnuplot1D_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamIO/Ggnuplot.h>		// Know GPdat structure
#include <Matrix/row_vector.h>		// Know GAMMA row vectors
#include <Matrix/matrix.h>		// Know GAMMA matrices
#include <string>                       // Know libstdc++ strings
#include <vector>                       // Know libstdc++ STL vectors

//_____________________________________________________________________________
// A                            ASCII Data Output
//_____________________________________________________________________________

/* These functions are used to place the data points into the Gnuplot ASCII
   data file that will be used for plotting. It uses the parameters found in
   the GPdat structure as an indication of how the users wishes the output data
   to be formatted.

   Here is what exists in G that affects the output data:

	{ xmax, xmin }	: Initial & final values on the x-axis. These are set
                          to zero by default and hence signal that there should
                          be no axis data values output. If the values differ
                          then the axis values will be output.
        compress        : A boolean flag which indicates if data compression
                          should be attempted. If there are no axis values
                          output then this is inactive. If active and set true
                          then any place three or more adjacent data points
                          have y-values that do not significantly differ, the
                          middle points (not 1st & last in series) are skipped.
        riflag		: Flag whether real, imaginary, or complex data is
                          plotted. riflag>=0: plot reals, riflag!=0 plot imags.
			  Thus riflag=0 plots reals only, riflag<0 plots imags
                          only and riflag >0 plots both reals and imaginaries.

		// Input		ofstr	: Output ASCII file stream
		//			vx	: Vector of data points
		//			G	: Gnuplot data controls
		// Output		void	: The file ofstr will be
		//				  modified by the addition of
		//				  data points found in vx in
		//				  accordance with flags
		//				  found in G.                */

MSVCDLL std::vector<double> GP1DXAxis(const   GPdat& G, int np);
MSVCDLL std::vector<int>    GP1DPtFlags(const GPdat& G, int np,  const row_vector& vx, bool real);
MSVCDLL void                GPXYData(std::ofstream& ofstr, const row_vector& vx, const GPdat& G);

//_____________________________________________________________________________
// B                        ASCII Load File Output
//_____________________________________________________________________________

/* These functions are used to place the Gnuplot commands into the Gnuplot
   ASCII load file that will be called on to invoke plotting. It uses the
   parameters found in the GPdat structure as an indication of how the user
   wishes the output plot to appear.                                         */

MSVCDLL void GPXYLoad(std::ofstream& ofstr, const GPdat& G);

//_____________________________________________________________________________
// C                       Single Vector XY Plots
//_____________________________________________________________________________

/* These GP_1D functions take a single vector of data desired to be plotted and
   writes it to either a new file or to an open filestream.  The data may be
   output either as only reals, only imaginaries, or real/imag back to back. 
   These functions only deal with ASCII output and they do NOT interactively
   produce any plots.
 
           Input        filename  : Output filename
        		ofstr     : An output filestream
                        vx        : Data vector
                        ri        : Flag for real versus imaginary plot
                                    ri = 0  : reals only (default)
                                    ri = <0 : imaginaries only
                                    ri = >0 : both 
                        xmin      : Plot horizontal (x=0) point value
                        xmax      : Plot horizontal (x=end) point value
                        cutoff    : Data compression cutoff
           Return                 : Void, file out is created/modified
        			    or output filestream is modified         */

MSVCDLL void GP_1D(std::ofstream&     ofstr,    const row_vector& vx, const GPdat& G);
MSVCDLL void GP_1D(const std::string& filename, const row_vector& vx, const GPdat& G);
MSVCDLL void GP_1D(const std::string& filename, const row_vector &vx,
                      int ri=0, double xmin=0, double xmax=0, double cutoff=0);
MSVCDLL void GP_1D(std::ofstream& ofstr,       const row_vector& vx,
                      int ri=0, double xmin=0, double xmax=0, double cutoff=0);

//_____________________________________________________________________________
// D                       Multi-Vector XY Plots
//_____________________________________________________________________________

/* The GP_1Dm functions take multiple vectors of data desired to be plotted
   together and writes them to either a new file or to an open filestream.
   The data may be output either as only reals, only imaginaries, or real/imag
   back to back. These functions only deal with ASCII output and they do NOT
   interactively produce any plots.
 
        // Input        filename  : Output filename
        //              vx        : Array of data vectors
        //              N         : Number of vectors
        //              ri        : Flag for real versus imaginary plot
        //                          ri = 0  : reals only (default)
        //                          ri = <0 : imaginaries only
        //                          ri = >0 : both 
        //              xmin      : Plot horizontal (x=0) point value
        //              xmax      : Plot horizontal (x=end) point value
        //              cutoff    : Data compression cutoff
        // Return                 : Void, file out is modified              */

MSVCDLL void GP_1Dm(const std::string& filename, row_vector* vxs, int N,
                     int ri=0, double xmin=0, double xmax=0, double cutoff=0);
MSVCDLL void GP_1Dm(const std::string& filename, const std::vector<row_vector>& vxs,
                     int ri=0, double xmin=0, double xmax=0, double cutoff=0);
MSVCDLL void GP_1Dm(const std::vector<std::string>& fnames, const std::vector<row_vector>& vxs,
                     int ri=0, double xmin=0, double xmax=0, double cutoff=0);

//_____________________________________________________________________________
// E                      Gnuplot Interactive Plots
//_____________________________________________________________________________


/* The function GPXYLoad (see above) is typically used to add Gnuplot commands
   into an existing ASCII output file stream that are suitable for X vs. Y type
   of plots. That function will NOT open the file stream, nor will it close the
   file stream, nor will it run any system command so that Gnuplot is invoked
   and a plot shown on screen. These functions use GPXYLoad to fill up an ASCII
   load file for Gnuplot X vs. Y plots, but these also handle opening/closing
   that file as well as running Gnuplot itself so that a plot is displayed.   */

	// Input	G	  : Gnuplot data controls
	//      	gnumacro  : A filename for gnu macro file
	//		file1D    : A filename of 1D gnuplot data
	//		join      : Flags whether to join points
	//      	basedir   : Base directory name
	//		basename  : Base file name
	//		term      : Gnuplot terminal type
	//		N	  : Number of data files           (ASCII)
	//		files	  : Pointer to array of file names (ASCII)
	//         OR   files	  : Array of data file names       (ASCII)
	// Return	void	  : Information concerning how
	//			    to plot in gnuplot is output
	//			    Gnuplot then run with load file

MSVCDLL void GP_1Dplot(const GPdat& G);
MSVCDLL void GP_1Dplot(const std::string& gnumacro, const std::string& file1D, int join=1);
MSVCDLL void GP_1Dplot(const std::string& gnumacro, const std::string& file1D,
                    const std::string& term, const std::string& ofile, int join=1);
MSVCDLL void GP_1Dplot(const std::string& basedir, const std::string& basename,
                                            const std::string& term, int join=1);

/* The next two functions will generate a single Gnuplot load file but they
   will plot multiple ASCII data files, in particular those data files whose
   names are found in the array files.                                       */

MSVCDLL void GP_1Dplot(const std::string& gnumacro, int N,             std::string* files, int join=1);
MSVCDLL void GP_1Dplot(const std::string& gnumacro, const std::vector<std::string>& files, int join=1);
 
//_____________________________________________________________________________
// F                    Interactive Output Functions
//_____________________________________________________________________________

/* These are the simplest functions to use when using Gnuplot from within GAMMA
   to produces interactive plots on the screen. These will generate an ASCII
   data file, generate a Gnuplot ASCII load file, then call Gnuplot to use the
   load file and display the plot to screen. Since many different parameters
   can be used to change both the way the data is output and the way the data
   is displayed, the most convenient function calls will make use of a GPdat
   structure with the proper setttings previously specified.                  */


MSVCDLL void Gnuplot1D(const std::string& name, const row_vector &vx,
          int ri=0, double xmin=0, double xmax=0, int join=1, double cutoff=0);

MSVCDLL void Gnuplot1D(const std::string& name, const std::vector<row_vector>& vxs,
                     int ri=0, double xmin=0, double xmax=0,
                             bool samefile=false, int join=1, double cutoff=0);

MSVCDLL void Gnuplot1D(const std::vector<std::string>& names, const std::vector<row_vector>& vxs,
                    int ri=0, double xmin=0, double xmax=0,
                            bool samefile=false, int join=1, double cutoff=0);

//____________________________________________________________________________
// G                         Gnuplot Output Blurbs
//____________________________________________________________________________
                                                                   
MSVCDLL void GP_1Dblurb(std::ofstream& ostr, const std::string& plotname);
 
        // Input        ostr      : An output stream
        //              plotname  : A gnuplot 1D plotname
        // Return       void      : Information concerning how
        //                          to plot in gnuplot is output
        //                          into the ofstream

#endif 						// Ggnuplot1D.h
