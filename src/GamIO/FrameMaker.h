/* FrameMaker.h *************************************************-*-c++-*-
**									**
**	                          G A M M A 				**
**								 	**
**	FrameMaker                                    Interface		**
**						 			**
**	Copyright (c) 1991, 1992				 	**
**	Tilo Levante, Scott Smith, Beat Meier				**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fuer physikalische Chemie				 	**
**	8092 Zurich / Switzerland				 	**
**								 	**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** The FrameMaker Library provides functions to create different 	**
** structures readable in the desktop publishing program FrameMaker.	**
**								 	**
*************************************************************************/

#ifndef   GFrameMaker_h_		// Is this file already included?
#  define GFrameMaker_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <GamIO/FrameMakerP.h>		// Know about FM parameters
#include <Matrix/row_vector.h>		// Know about row vectors
#include <Matrix/matrix.h>		// Know about matrices
#include <Level1/coord_vec.h>		// Know about coordinate vectors 
#include <HSLib/GenOp.h>		// Know about operators
#include <LSLib/SuperOp.h>		// Know about superoperators
#include <string>			// Know about libstdc++ string
//#include <Deprecated/molecule.h>	// Know about molecules

#ifndef PI
#define PI 3.141592653589
#endif

#define FMPAGEWIDTH  19.0;		// 19.05 cm =  7.5 in
#define FMPAGEHEIGHT 25.0;		// 25.40 cm = 10.0 in
#define FMPLmax      128		// Default polyline size
#define FM_PP        2			// Digits past decimal point (FMMatrix)
#define FM_PT        0.001		// Default threshold         (FMMatrix)
extern const coord DefOrient;		// Default 3D Plot Orientation


// sosi
// Next Add Class MIFFILE To Here......
// Perhaps Also Class PolyLine
 
// ----------------------------------------------------------------------------
// ------------------------- END OF FMCONTROLS CLASS --------------------------
// ----------------------------------------------------------------------------


struct FMcont				// FM contour parameters
  {
  double thresh;			// Threshold value
  int steps;				// Number of contours
  double CLI;				// Contour level increment
  double CLM;				// Contour level modifier
  int CPN;				// Postive/negative contour flag
  double hsize;				// FM horizontal plot size (cm)
  double vsize;				// FM vertical plot size (cm)
  double hscale;			// Horizontal scale (cm/pt)
  double vscale;			// Vertical scale (cm/pt)
  double vlow;				// Lower vertical axis label
  double vhigh;				// Upper vertical axis label
  double hlow;				// Lower horizontal axis label
  double hhigh;				// Upper horizontal axis label
  double dmax;				// Data array maximum
  double dmin;				// Data array minimum
  int debug;				// Debugging level
  };


struct FMclev				// FM contour level
  {
  int vdim;				// Vector dimension
  row_vector vxi;			// Vector of initial points
  row_vector vxf;			// Vector of final points
  int npts;				// Number of contour points
  double threshold;			// Contour threshold
  int level;				// Contour level
  int posID;				// Positive contours base FM ID
  int negID;				// Negative coutours base FM ID
  int FMID;				// Contour level FM ID
  int nconts;				// Number of countour lines
  int CLmaxnum;				// Max contour lines per level
  int CLmaxsize;			// Max size of a contour line
  int conpos;				// Number of positive contours
  int conneg;				// Number of negative contours
  int penpos;				// Positive contours pen index
  int penneg;				// Negative contours pen index
  };


struct FMxy				// FM xy parameters
  {
  double hsize;				// FM horizontal plot size (cm)
  double vsize;				// FM vertical plot size (cm)
  double hmin;				// Horizontal axis label start
  double hmax;				// Horizontal axis label finish
  double vmin;				// Vertical plot label start
  double vmax;				// Vertical plot label finish
  int PLmax;				// Maximum polyline size
  int debug;				// Debugging level
  };

//_____________________________________________________________________________
// A                     GAMMA FrameMaker System Interface
//_____________________________________________________________________________

        // Input        none    : Output filename
        //              warn    : Warning flag
        // Return       FME     : A string for the FrameMaker
        //                        executable command (hopefully)

MSVCDLL std::vector<std::string> FMSeekStrings();	// Strings for gnuplot path
MSVCDLL std::string FMFind(bool vocal=false);		// Look for gnuplot exec.
MSVCDLL std::string FMExec(int warn=0);			// FrameMaker executable name
//void RunGnuplot(const std::string& gnumacro);	// Run Gnuplot

// ____________________________________________________________________________
//                        FrameMaker 1D Plot Functions
// ____________________________________________________________________________


MSVCDLL void FM_1D(const std::string& filename, const row_vector &vx, double xsize=14,
                      double ysize=14, double pmin=0, double pmax=1, int ri=0);

	// Input	filename  : Output filename
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


MSVCDLL void FM_1D(const std::string& filename, const row_vector& vx, FMPar& FMP);

        // Input        filename  : Output filename
        //              vx        : Data vector
        //              FMP       : FrameMaker plot parameters
        // Return                 : A FrameMaker MIF file is
        //                          made called filename to contain
        //                          1D plot(s) of data in vector vx
        //                          according to parameters in FMP



MSVCDLL void FM_1Dm(const std::string& filename, int nvec, row_vector *vx, double xsize=14,
	       double ysize=14, double pmin=0, double pmax=1, int ri=0);

MSVCDLL void FM_1Dm(const std::string& filename, const std::vector<row_vector>& vxs,
   double xsize=14, double ysize=14, double pmin=0, double pmax=1, int ri=0);

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

MSVCDLL void FrameMaker1D(const std::string& name, const row_vector& vx,
                                       int ri=0, double xmin=0, double xmax=1);

MSVCDLL void FrameMaker1D(const std::string& name, const std::vector<row_vector>& vxs,
                                       int ri=0, double xmin=0, double xmax=1);

// ____________________________________________________________________________
// C                    FrameMaker Parametric Plot Functions
// ____________________________________________________________________________


MSVCDLL void FM_xyPlot (const std::string& filename, row_vector &vx, double xsize=14,
  double ysize=14, double xmin=0, double xmax=0, double ymin=0, double ymax=0);

	// Input	filename  : Output filename
	// 		vx        : Data vector
	//		xsize     : Plot horizontal (x) dimension in cm
	//		ysize     : Plot vertical (y) dimension in cm
	//		xmin	  : Minimum value for x-scale
	//		xmax	  : Maximum value for x-scale
	//		ymin	  : Minimum value for y-scale
	//		ymax	  : Maximum value for y-scale
	// Return		  : Void, file out is modified


MSVCDLL void FM_xyPlot(const std::string& filename, row_vector &vx, FMxy& FMXY);

        // Input        filename  : Output filename
        //              vx        : Data vector
        //              FMXY      : FrameMaker xy-plot controls
        // Return                 : Void. xy-plot is output to file

 
// ____________________________________________________________________________
// D                    FrameMaker Scatter Plot Functions
// ____________________________________________________________________________
 

MSVCDLL void FM_scatter(const std::string& filename, row_vector &vx, char z,
	                                   double xsize=14, double ysize=14);

	// Input	filename  : FrameMaker mif file name
	// 		vx        : Data vector containing points
	//		z	  : Character to mark points
	// 		xsize     : Horizontal plot size in cm.
	// 		ysize     : Vertical plot size in cm.
	// Return		  : Void, FrameMaker plot file output


MSVCDLL void FM_scatter(const std::string& filename, row_vector &vx, char z, FMxy& FMXY);

        // Input        filename: Output filename
        //              vx      : Data vector
        //              z       : Character to use to mark points
        //              FMXY    : FrameMaker xy-plot controls
        // Return               : Void, FrameMaker plot file output


MSVCDLL void FM_scatterm(const std::string& filename, int nvec, row_vector *vx, char z, FMxy& FMXY);

        // Input        filename: Output filename
        //              nvec    : Number of input vectors
        //              vx      : Array of data vectors
        //              z       : Character to use to mark points
        //              FMXY    : FrameMaker xy-plot controls
        // Return               : Void, FrameMaker plot file output



MSVCDLL void FM_scatter(const std::string& filename, row_vector &vx, int sides=0,
		        double PGsize=0, double xsize=14, double ysize=14);

	// Input	filename  : FrameMaker mif file name
	// 		vx        : Data vector containing points
	//		sides	  : Number of sides in polygon
	// 		PGsize    : Polygon plot size in cm.
	// 		xsize     : Horizontal plot size in cm.
	// 		ysize     : Vertical plot size in cm.
	// Return		  : Void, FrameMaker plot file output

 
// ____________________________________________________________________________
// E                    FrameMaker Histogram Plot Functions
// ____________________________________________________________________________

MSVCDLL void FM_histogram(const std::string& filename, row_vector &vx,
		                int bins=0, double xsize=14, double ysize=14);

	// Input	filename  : FrameMaker mif file name
	// 		vx        : Data vector containing points
	//		bin       : Number of bins in histogram
	// 		xsize     : Horizontal plot size in cm.
	// 		ysize     : Vertical plot size in cm.
	// Return		  : Void, FrameMaker plot file output

 
// ____________________________________________________________________________
// F                     FrameMaker Contour Plot Functions
// ____________________________________________________________________________


MSVCDLL void FM_contour(const std::string& filename, const matrix &mx,
		 double low_t, int steps=0, double CLI=0,
		 double CLM=1, int CPN=1, double xsize=15, double ysize=15);

	// Input	filename  : Output filename
	// 		mx        : Data matrix
	//		low_t     : Lowest threshold
	//		steps     : Number of contours to take
	//		CLI       : Contour level increment
	//		CLM 	  : Contour level modifier
	//		CPN 	  : Positive/Negetive contour flag
	//		xsize     : Plot horizontal (x) dimension in cm
	//		ysize     : Plot vertical (y) dimension in cm
	// Return		  : Void, file out is modified

 
MSVCDLL void FM_contour(const std::string& filename, const matrix &mx, FMcont& FMCP);
 
        // Input        filename  : Output filename
        //              mx        : Data matrix
        //              FMCP      : Contour plot parameters
        // Return                 : Void, file out is modified
        // Note                   : AUXILIARY FUNCTION FOR FM_contour


MSVCDLL void FM_contour(const std::string& filename, const matrix &mx, FMcont& FMCP, double* levels);

        // Input        filename  : Output filename
        //              mx        : Data matrix
        //              FMCP      : Contour plot parameters
        //              levels    : Desired contour levels 
        // Return                 : Void, file out is modified
        // Note                   : Uses GAMMA calculated levels
 

MSVCDLL void FM_contour(const std::string& filename, const matrix &mx, FMcont& FMCP, FMclev& FMCL,
                                                         double* Ts, int* IDs);
        // Input        filename  : Output filename
        //              mx        : Data matrix
        //              FMCP      : Contour plot parameters
        //              FMCL      : Contour level data
        //              Ts        : Contour levels
        //              IDs       : Contour level IDs (FM)
        // Return                 : Void, file out is modified
 
// ____________________________________________________________________________
// G                        FRAMEMAKER MMF FUNCTIONS
// ____________________________________________________________________________


MSVCDLL void FM_Matrix(const std::string& filename, const gen_op& Op,
	                         int prec = FM_PP,  double threshold = FM_PT);

	// Input	out       : output stream
	//		Op	  : general operator
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    of an operator to a mathematical expression


MSVCDLL void FM_Matrix(const std::string& filename, const super_op& LOp,
		 	int prec = FM_PP,  double threshold = FM_PT);

	// Input	out       : output stream
	//		LOp	  : superoperator
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    of a superoperator to a mathematical expression


MSVCDLL void FM_Matrix(const std::string& filename, const matrix& mx,
			 int prec = FM_PP, double threshold = FM_PT);

	// Input	out       : output stream
	//		mx	  : matrix
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    to a mathematical expression


MSVCDLL void FM_Mat_Plot(const std::string& filename, const gen_op& Op, double threshold = FM_PT);

	// Input	out       : output stream
	//		Op	  : general operator
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    of an operator as a schematic
	//			    graphics object


MSVCDLL void FM_Mat_Plot(const std::string& filename, const super_op& LOp, double threshold = FM_PT);

	// Input	out       : output stream
	//		LOp	  : superoperator
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    of a superoperator as a schematic
	//			    graphics object


MSVCDLL void FM_Mat_Plot(const std::string& filename, const gen_op& Op,
		 		const gen_op& Op_ref, double threshold = FM_PT);

	// Input	out       : output stream
	//		Op	  : general operator
	//		Op_ref	  : reference general operator
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    representing the differences between
	//			     Op and Op_ref.  The output is as a schematic
	//			    graphics object


MSVCDLL void FM_Mat_Plot(const std::string& filename, const super_op& LOp,
		 		const super_op& LOp_ref, double threshold = FM_PT);

	// Input	out       : output stream
	//		LOp	  : general superoperator
	//		LOp_ref	  : reference superoperator
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    representing the differences between
	//			    LOp and LOp_ref.  The output is as a schematic
	//			    graphics object


MSVCDLL void FM_Mat_Plot (const std::string& filename, const matrix &mx, double threshold = FM_PT);

	// Input	out       : output stream
	//		mx	  : matrix
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    as a schematic graphics object


MSVCDLL void FM_Mat_Plot (const std::string& filename, const matrix &mx,
				 const matrix &mx_ref, double threshold = FM_PT);

	// Input	out       : output stream
	//		mx	  : matrix
	//		mx_ref	  : reference general operator
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    representing the differences between
	//			     Op and Op_ref.  The output is as a schematic
	//			    graphics object



//void FM_molecule (const std::string& filename, molecule mol, double radius=5,
//	                    double alpha=0, double beta=15, double gamma=-15);


MSVCDLL void FM_Mat_Tbl(const std::string& filename, const matrix& mx, double threshold = FM_PT);

	// Input	out       : output stream
	//		mx	  : matrix
	//		prec	  : precision, default (FM_PP)
	//		threshold : magnitude threshold, default (FM_PT)
	// Output	none      : Function is void.  Puts out the matrix
	//			    as a schematic graphics object



#endif 							// FrameMaker.h
