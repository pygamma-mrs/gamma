/* Ggnuplot1D.cc ************************************************-*-c++-*-
**									**
**	                           G A M M A 				**
**								 	**
**	Gnuplot Support X Vs. Y Plots	         Implementation   	**
**				        	 		 	**
**      Scott A. Smith							**
**      Copyright (c) 2002			 	                **
**									**
**      $Header: $
**	                         		         	 	**
*************************************************************************/

/**********************************************************tation   	**
**				        	 		 	**
**      Scott A. Smith							**
**      Copyright (c) 2002			 	                **
**									**
**      $Header: $
**	                         		         ***************
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
 
#ifndef   Ggnuplot1D_cc_		// Is file already included?
#  define Ggnuplot1D_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamIO/Ggnuplot1D.h>		// Include gnuplot header
#include <GamIO/Ggnuplot.h>		// Include for SetLineType function
#include <Basics/StringCut.h>		// Include GAMMA Gdec function
#include <Basics/Gconstants.h>		// Include definition of RAD2DEG
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <iostream>			// Include libstdc++ file streams
#include <string>			// Include libstdc++ ANSI strings

//_____________________________________________________________________________
// A                        ASCII Data File Output
//_____________________________________________________________________________

/* These functions are used to place the data points into the Gnuplot ASCII
   data file that will be used for plotting. It uses the parameters found in
   the GPdat structure as an indication of how the user wishes the output data
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

std::vector<double> GP1DXAxis(const GPdat& G, int np)
  {
  std::vector<double> XVals;			// X axis values for Y vs. X plot
  if(G.xmin == G.xmax) return XVals;		// Nothing if no axis specified
  XVals = std::vector<double>(np);		// Set our array size
  double delx = (G.xmax-G.xmin)/double(np-1);	// The X-axis increment
  for(int i=0; i<np; i++)			//   Loop all points and set
    XVals[i] = G.xmin + double(i)*delx;		//   x-axis values
  return XVals;					// Return our X axis values
  }

std::vector<int> GP1DPtFlags(const GPdat& G, int np, const row_vector& vx, bool real)
  {
  std::vector<int> PtFlags(np, 1);		// Point plot flags
  double Icut = G.cutoff;
  if(!Icut)					// If there is no intensity
    {						// resolution, we set it now
    double pc = 1.e-3;
    if(real) Icut = pc*(vx.maxRe()-vx.minRe());
    else     Icut = pc*(vx.maxIm()-vx.minIm());
    }

  if(real)					// Set array if real data
    {
    double lpt = vx.getRe(0);			//     Previous point intensity
    for(int i=1; i<np-1; i++)			//     Examine successive pts
      {						//     for intensity changes
      if(fabs(vx.getRe(i)-lpt) > Icut)		//	 If intensity differs
        {					//       from last plotted
        PtFlags[i-1] = 1;			//       point flag that it &
        PtFlags[i]   = 1;			//       previous pts plotted
        lpt = vx.getRe(i);			//       Set this as last pt
        }
      else PtFlags[i] = 0;			//    Don't plot if no change
      }
    }
  else
    {
    double lpt = vx.getIm(0);			//     Previous point intensity
    for(int i=1; i<np-1; i++)			//     Examine successive pts
      {						//     for intensity changes
      if(fabs(vx.getIm(i)-lpt) > Icut)		//	 If intensity differs
        {					//       from last plotted
        PtFlags[i-1] = 1;			//       point flag that it &
        PtFlags[i]   = 1;			//       previous pts plotted
        lpt = vx.getIm(i);			//       Set this as last pt
        }
      else PtFlags[i] = 0;			//    Don't plot if no change
      }
    }
  return PtFlags;
  }

void GPXYData(std::ofstream& ofstr, const row_vector& vx, const GPdat& G)
  {
  int np = vx.size();				// Number of points in vector
  std::vector<double> XVals;			// For x-axis values (if used)
  XVals = GP1DXAxis(G, np);			// Fill array with any values
  bool XAxis = XVals.size()?true:false;		// Flag if axis values output
  std::vector<int> RePtFlags;			// Real pt. compression flags
  std::vector<int> ImPtFlags;			// Imag pt. compression flags
  if(XAxis && G.compress)			// Compression only considered
    {						// if X axis values & compress
    if(G.riflag >= 0)				//   If plotting reals, get
      RePtFlags = GP1DPtFlags(G, np, vx, true);	//   flags of which pts plotted
    if(G.riflag)				//   If plotting imags, get
      ImPtFlags = GP1DPtFlags(G, np, vx, false);//   flags of which pts plotted
    }
  if(G.riflag >= 0)				// Write data for real plot
    {
    if(!XAxis)					// Just output reals if no
      {			 			// x-ordinates are needed
      for(int i=0; i<np; i++)
        ofstr << vx.getRe(i) << "\n";
      }
    else if(!RePtFlags.size())			// Output {x,y} pairs if
      for(int i=0; i<np; i++) 			// x-ordinates but no
        { 					// compression is used
        ofstr << XVals[i] << "  ";
        ofstr << vx.getRe(i) << "\n";
        }
    else					// Output {x,y} pairs with
      for(int i=0; i<np; i++) 			// x-ordinates and data
        { 					// compression
        if(RePtFlags[i])
          {
          ofstr << XVals[i] << "  ";
          ofstr << vx.getRe(i) << "\n";
          }
        }
    }
  if(G.riflag > 0) ofstr << std::endl; 		// Blank line to separate plots
  if(G.riflag)					// Write data for imaginary plot
    {
    if(!XAxis)					// Just output imags if no
      for(int i=0; i<np; i++) 			// x-ordinates are needed
        ofstr << vx.getIm(i) << "\n";
    else if(!ImPtFlags.size())			// Output {x,y} pairs if
      for(int i=0; i<np; i++) 			// x-ordinates but no
        { 					// compression is used
        ofstr << XVals[i]    << "  ";
        ofstr << vx.getIm(i) << "\n";
        }
    else					// Output {x,y} pairs with
      for(int i=0; i<np; i++) 			// x-ordinates and data
        { 					// compression
        if(ImPtFlags[i])
          {
          ofstr << XVals[i]    << "  ";
          ofstr << vx.getIm(i) << "\n";
          }
        }
    }
  }

//_____________________________________________________________________________
// B                        ASCII Load File Output
//_____________________________________________________________________________

/* These functions are used to place the Gnuplot commands into the Gnuplot
   ASCII load file that will be called on to invoke plotting. It uses the
   parameters found in the GPdat structure as an indication of how the user
   wishes the output plot to appear.                                         */


void GPXYLoad(std::ofstream& ofstr, const GPdat& G)
  {
  if(G.basedir  != "")
    ofstr << "cd \"" << G.basedir << "\"\n";		// Switch to base dir
  SetLineType(ofstr, G.join);				// Set plot line type
  ofstr << "set noparametric\n";			// Plot x vs y simple 1D
  if(G.term != "")
    ofstr << "set terminal " << G.term    << "\n";	// Set the terminal type
  if(G.outfile != "")
    ofstr << "set output \"" << G.outfile << "\"\n";	// Set the output name
  if(G.nokey)        ofstr << "set nokey\n";		// Set nokey (no auto filename)
  if(G.title!= "")   ofstr << "set title \"" 		// Label with title
                           << G.title<< "\"\n";
  if(G.xlabel != "") ofstr << "set xlabel \""		// Label x-axis
                             << G.xlabel << "\"\n";
  if(G.ylabel != "") ofstr << "set ylabel \""		// Label y-axis 
                             << G.ylabel << "\"\n";
  if(G.tics=="out")  ofstr << "set tics out\n";	// Set axis tic marks out

  if(!G.border) ofstr << "set noborder\n";		// Removes gnuplot border
  if(!G.xtics)  ofstr << "set noxtics\n";		// Stop x-axis tic marks
  if(!G.ytics)  ofstr << "set noytics\n";		// Stop y-axis tic marks
  if(!G.xaxis)  ofstr << "set noxzeroaxis\n";		// Stop x-axis draw
  if(!G.yaxis)  ofstr << "set noyzeroaxis\n";		// Stop y-axis draw
  if(G.xsize && G.ysize) 				// Scale the plot
    ofstr << "set size " << G.xsize << ", "
          << G.ysize     << "\n";
  ofstr << "plot \"" << G.data << "\"\n";		// Command for plot
  }

//_____________________________________________________________________________
// C                       Single Vector XY Plots
//_____________________________________________________________________________

/* These functions either generate an external ASCII file containing data
   points suitable for use in Gnuplot or add points to such a file which is
   previously existing. The output file can be named anything, however it is
   recommended it be called "foo.asc" where the .asc suffix indicates its
   nature. There are several flags which affect the data output to the ASCII
   file. These are as follows:

           Input        filename  : Output filename (ASCII)
        		ofstr     : An output filestream
                        vx        : Data vector
			G         : Gnuplot control parameteres
                        ri        : Flag for real versus imaginary plot
                                    ri = 0  : reals only (default)
                                    ri = <0 : imaginaries only
                                    ri = >0 : both 
                        xmin      : Plot horizontal (x=0) point value
                        xmax      : Plot horizontal (x=end) point value
	  		cutoff    : Data compression cutoff

           Return                 : Void. Output file filename is generated 
                                    or overwritten, or the ASCII data is
				    added to the output file stream.        */

void GP_1D(std::ofstream& ofstr, const row_vector& vx, const GPdat& G)
  { GPXYData(ofstr, vx, G); }

void GP_1D(const std::string& filename, const row_vector& vx, const GPdat& G)
  {
  std::ofstream ofstr(filename.c_str());	// Open file for plotting
  GPXYData(ofstr,vx,G);				// Use GP_1D overload
  ofstr.close();				// Close the file stream
  }

void GP_1D(const std::string& filename, const row_vector& vx,
                             int ri, double xmin, double xmax, double cutoff)
  {
  std::ofstream ofstr(filename.c_str());	// Open file for plotting
  GP_1D(ofstr,vx,ri,xmin,xmax,cutoff);		// Use GP_1D overload
  ofstr.close();				// Close the file stream
  }

void GP_1D(std::ofstream& ofstr, const row_vector& vx,
                             int ri, double xmin, double xmax, double cutoff)
  {
  GPdat G;					// Gnuplot infrastructure
  G.xmin   = xmin;				// Set x-axis minimum
  G.xmax   = xmax;				// Set x-axis maximum
  G.cutoff = cutoff;				// Set intensity resolution
  G.riflag = ri;				// Set {real,imag,complex}
  GPXYData(ofstr, vx, G);			// Use function overload
  }

//_____________________________________________________________________________
// D                         Multi-Vector XY Plots
//_____________________________________________________________________________

/* These functions will plot multiple vectors within a single graph. Each plot
   within the graph will be associated with single row_vector. The vector ASCII
   data can be set to reside in either a single file or in individual files. */

        // Input        filename  : Output filename
        //              vx        : Array of data vectors
	//		N	  : Number of vectors
        //              ri        : Flag for real versus imaginary plot
        //                          ri = 0  : reals only (default)
        //                          ri = <0 : imaginaries only
        //                          ri = >0 : both 
        //              xmin      : Plot horizontal (x=0) point value
        //              xmax      : Plot horizontal (x=end) point value
	//		cutoff    : Data compression cutoff
        // Return                 : Void, file out is modified
	// Note			  : All data is written into a single
	//			    ASCII file, each set separated by
	//			    a blank line. All plots will be the
	//			    same color.

void GP_1Dm(const std::string& fname, row_vector* vxs, int N,
                              int ri, double xmin, double xmax, double cutoff)
  {
  std::vector<row_vector> vs;			// Array of row_vectors
  for(int i=0; i<N; i++)			// Loop over row_vectors
    vs.push_back(vxs[i]);			// and put in STL vector
  GP_1Dm(fname,vs,ri,xmin,xmax,cutoff);		// Use overload function
  }

void GP_1Dm(const std::string& fname, const std::vector<row_vector>& vxs,
                              int ri, double xmin, double xmax, double cutoff)
  {
  std::ofstream ofstr(fname.c_str());		// Open file for plotting
  for(unsigned i=0; i<vxs.size(); i++)
    {
    GP_1D(ofstr,vxs[i],ri,xmin,xmax,cutoff);	// Use GP_1D overload
    ofstr << "\n";
    }
  ofstr.close();				// Close the file stream
  }

void GP_1Dm(const std::vector<std::string>& fnames, const std::vector<row_vector>& vxs,
                              int ri, double xmin, double xmax, double cutoff)
  {
  for(int i=0; i<int(vxs.size()); i++)
    {
    std::ofstream ofstr((fnames[i]).c_str());	// Open file for plotting
    GP_1D(ofstr,vxs[i],ri,xmin,xmax,cutoff);	// Use GP_1D overload
    ofstr.close();				// Close the file stream
    }
  return;
  }

// ____________________________________________________________________________ 
// E                        Gnuplot Interactive Plots 
// ____________________________________________________________________________

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

void GP_1Dplot(const GPdat& G)
  {
  std::string gmacro = G.basedir + G.macro; 	// Gnuplot macro file
  std::ofstream gnuload(gmacro.c_str());	// Gnuplot commands file
  GPXYLoad(gnuload, G);				// Generic plot commands
  CloseMacro(gnuload);				// Add pause+exit, close macro
  RunGnuplot(gmacro);				// Run gnuplot with gnumacro
  }

void GP_1Dplot(const std::string& gnumacro, const std::string& file1D, int join)
  {
  GPdat G;					// Gnuplot infrastructure
  G.macro  = gnumacro;				// Set load file name in G
  G.join   = join;				// Set line joining flag in G
  G.data   = file1D;				// Set data file name
  GP_1Dplot(G);
  }

void GP_1Dplot(const std::string& gnumacro, const std::string& file1D,
                             const std::string& term, const std::string& ofile, int join)
  {
  GPdat G;					// Gnuplot infrastructure
  G.macro   = gnumacro;				// Set load file name in G
  G.data    = file1D;				// Set ASCII data file
  G.term    = term;				// Set terminal type
  G.outfile = ofile;				// Set output plot file
  G.join    = join;				// Set line joining flag in G
  GP_1Dplot(G);
  }

void GP_1Dplot(const std::string& basedir, const std::string& basename,
                                                 const std::string& term, int join)
  {
  GPdat G;					// Gnuplot infrastructure
  G.join    = join;				// Set line joining flag in G
  G.term    = term;				// Set terminal type
  G.basedir = basedir + "/";
  G.macro   = basename + ".asc";
  G.data    = basename + ".gnu";
  G.outfile = basename + "." + term;
  GP_1Dplot(G);
  }

/* The next two functions will generate a single Gnuplot load file but they
   will plot multiple ASCII data files, in particular those data files whose
   names are found in the array files.                                       */

void GP_1Dplot(const std::string& gnumacro, int N, std::string* files, int join)
  {
  std::vector<std::string> fs;			// Filenames array
  for(int i=0; i<N; i++) fs.push_back(files[i]);// Fill filename array
  GP_1Dplot(gnumacro, fs, join);		// Use overload
  }

void GP_1Dplot(const std::string& gnumacro, const std::vector<std::string>& files, int join)
  {
  std::ofstream gnuload(gnumacro.c_str());		// File of gnuplot commands
  SetLineType(gnuload, join);			// Set plot line type
  gnuload << "set noparametric\n";		// Plot x vs y, as simple 1D
  gnuload << "plot "; 				// Command for gnuplot plots
  int N = int(files.size());			// Number of data files
  for(int i=0; i<N-1; i++)			// Loop all but last data file
    gnuload << "\"" << files[i] << "\", ";	// and set to be plotted
  gnuload << "\"" << files[N-1] << "\"\n";	// Set last data file plotted
  CloseMacro(gnuload);				// Add pause+exit, close macro
  RunGnuplot(gnumacro);				// Run gnuplot with gnumacro
  }

//_____________________________________________________________________________
// F                     Interactive Output Functions
//_____________________________________________________________________________

/* These are the simplest functions to use when using Gnuplot from within GAMMA
   to produces interactive plots on the screen. These will generate an ASCII
   data file, generate a Gnuplot ASCII load file, then call Gnuplot to use the
   load file and display the plot to screen. Since many different parameters
   can be used to change both the way the data is output and the way the data
   is displayed, the most convenient function calls will make use of a GPdat
   structure with the proper setttings previously specified.                  */


void Gnuplot1D(const std::string& name, const row_vector &vx,
                    int ri, double xmin, double xmax, int join, double cutoff)
  {
  std::string Afile = name + ".asc";		// Set data file name
  std::string Gfile = name + ".gnu";		// Set gnuplot load filename
  GP_1D(Afile, vx, ri, xmin, xmax, cutoff);	// Make data file (ASCII)
  GP_1Dplot(Gfile, Afile, join);		// Make gnuplot load file
  }						// and call gnuplot of plot

void Gnuplot1D(const std::string& name, const std::vector<row_vector>& vxs,
     int ri, double xmin, double xmax, bool samefile, int join, double cutoff)
  {
  if(samefile)					// If only 1 data file
    {
    std::string Afile = name + ".asc";		//   Set data file name
    std::string Gfile = name + ".gnu";		//   Set gnuplot load filename
    GP_1Dm(Afile, vxs, ri, xmin, xmax, cutoff);	//   Output single data file
    GP_1Dplot(Gfile, Afile, join);		//   Plot to screen
    }
  else						// If multiple data files
    {
    int nf = vxs.size();			//   Get number of data sets
    std::vector<std::string> Afs;		//   Array of data file names
    std::string Afile;				//   For data file name
    std::string Gfile = name + ".gnu";		//   Set gnuplot load filename
    for(int i=0; i<nf; i++)			//   Loop over data sets
      { 					//   and output a file for each
      Afile = name + Gdec(i) + ".asc";		//     Use this filename
      GP_1D(Afile,vxs[i],ri,xmin,xmax,cutoff);	//     Output single data file
      Afs.push_back(Afile);			//     Store file name in array
      }
    GP_1Dplot(Gfile, Afs, join);		//   Plot to screen
    }
  }

void Gnuplot1D(const std::vector<std::string>& names, const std::vector<row_vector>& vxs,
     int ri, double xmin, double xmax, int join, double cutoff)
  {
  int nf = vxs.size();				//   Get number of data sets
  std::vector<std::string> Afs;				//   Array of data file names
  std::string Afile;					//   For data file name
  std::string Gfile = names[0] + ".gnu";		//   Set gnuplot load filename
  for(int i=0; i<nf; i++)			//   Loop over data sets
    { 						//   and output a file for each
    Afile = names[i] + Gdec(i) + ".asc";	//     Use this filename
    GP_1D(Afile,vxs[i],ri,xmin,xmax,cutoff);	//     Output single data file
    Afs.push_back(Afile);			//     Store file name in array
    }
  GP_1Dplot(Gfile, Afs, join);			//   Plot to screen
  }

// ____________________________________________________________________________ 
// G                         Gnuplot Output Blurbs
// ____________________________________________________________________________ 

	// Input	ostr      : An output stream
	// 		plotname  : A gnuplot 1D plotname
	// Return	void	  : Information concerning how
	//			    to plot in gnuplot is output
	//			    into the ofstream

void GP_1Dblurb(std::ofstream& ostr, const std::string& plotname)
  {
  ostr << "\n\n\t\tGNUPlot 1D ASCII Data File "
       << plotname << " Complete";
  ostr << "\n\t\tCommands To Plot: 1.) set data style lines"
       << "\n\t\t                  2.) set noparametric"
       << "\n\t\t                  3.) plot \"" << plotname << "\"";
  ostr << "\n";
  }

#endif 						// Gnuplot1D.cc
