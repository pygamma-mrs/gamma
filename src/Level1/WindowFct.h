/* WindowFct.h **************************************************-*-c++-*-
**									**
**	                         G A M M A				**
**								 	**
**	Windowing & Other Spatial Functions            Interface	**
**						 			**
**	Copyright (c) 1990, 1999				 	**
**      Scott Smith, Tilo Levante					**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fuer physikalische Chemie				 	**
**	8092 Zurich / Switzerland		 			**
**							 		**
**      $Header: $
**								 	**
**								 	**
*************************************************************************/

/*************************************************************************
**						 			**
** Description							 	**
**								 	**
** Often various "windowing" functions are used in the processing of	**
** NMR data.  This modules provides some of these functions in order	**
** to accomplish some of these same processing step on GAMMA simulated	**
** (or imported) spectra.						**
**								 	**
*************************************************************************/

#ifndef   GWindowFct_h_			// Is file already included?
#  define GWindowFct_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Basics/Gconstants.h>		// Need default PI value
#include <Matrix/row_vector.h>		// Need to know row vectors
#include <Matrix/col_vector.h>		// Need to know column vectors

MSVCDLL void exponential_multiply(col_vector &db, double em=-4, int offset=0);
MSVCDLL void exponential_multiply(row_vector &db, double em=-4, int offset=0);

      // Input      db : row_vector
      //            em : see below
      //         offset: shift for the zero freqency
      // Output          db will be modified. db contains the
      //                 the original data multiplied with an
      //                 exponential function so that the last point
      //                 is multiplied by exp(em) and the first by 1.


// ______________________________________________________________________
// _________________________ WINDOW FUNCTIONS ___________________________
// ______________________________________________________________________


MSVCDLL row_vector exponential (int size, int offset=0, double alpha=0);

	// Input	size   : data block size
	//		offset : function offset
	//		alpha  : line broadening parameter
	// Output	BLK    : data block containing an exponential 
	//			 function having maximum = 1 at the
	//			 offset point and a half-height
	//			 linewidth of 2ln2*alpha [1.3863]


MSVCDLL row_vector Gaussian (int size, int offset=0, double sigma=42.0);

	// Input	size   : data block size
	//		offset : function offset
	//		sigma  : function width factor
	// Output	BLK    : data block containing a Gaussian function
	//			 having maximum = 1 at the offset point and
	//			 a half-height linewidth of 2.25*sigma points


MSVCDLL row_vector Hamming (int size, int offset=0);

	// Input	size   : data block size
	//		offset : function offset
	// Output	BLK    : data block containing a Hamming
	//			 function having maximum = 1 at the
	//			 offset point


MSVCDLL row_vector Hanning (int size, int offset=0);

	// Input	size   : data block size
	//		offset : function offset
	// Output	BLK    : data block containing a Hanning function
	//			 function having maximum = 1 at the
	//			 offset point


MSVCDLL row_vector hyperbol_sec (int size, int offset=0, double alpha=38.0);

	// Input	size   : data block size
	//		offset : function offset
	//		alpha  : function width factor
	// Output	BLK    : data block containing a hyperbolic secant
	//			 function having maximum = 1 at the offset point
	//			 and a half-height linewidth of 2.64*alpha points


MSVCDLL row_vector Kaiser (int size, double theta=PI, int offset=0);

	// Input	size   : data block size
	//		theta  : angle
	//		offset : function offset
	// Output	BLK    : data block containing a Kaiser function


MSVCDLL row_vector Lorentzian (int size, int offset=0, double alpha=1.0);

	// Input	size   : data block size
	//		offset : function offset
	//		alpha  : function width factor
	// Output	BLK    : data block containing a Lorentzian function
	//			 having maximum = 1 at the offset point and
	//			 a half-height linewidth of 2*alpha


MSVCDLL row_vector sin_square (int size, int offset=0);

	// Input	size  : data block size
	//		offset: function offset		   2
	// Output	BLK   : data block containing a sin (x) function
	//			having zero at the offset point and pi
	//			at the final block point


MSVCDLL row_vector sinc (int size, int offset, int inc);

	// Input	size   : data block size
	//		offset : function offset
	//		inc    : point increment to first node
	// Output	BLK    : data block containing a sinc(x) function
	//			 having maximum at the offset point and its
	//			 first node inc points away


MSVCDLL row_vector square_wave (int size, int start, int finish);

	// Input	size  : data block size
	//		start : starting point
	//		finish: finishing point
	// Output	BLK   : data block containing a square wave function
	//			having a value 1 between points start to
	//			finish and zero elsewhere

// ____________________________________________________________________________
//                            Random Noise Functions
// ____________________________________________________________________________

MSVCDLL row_vector Noise(int npts,         double maxN);
MSVCDLL void       Noise(row_vector& data, double maxN);

#endif 						// WindowFct.h
