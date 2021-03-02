/* WindowFct.cc *************************************************-*-c++-*-
**									**
**	                         G A M M A				**
**								 	**
**	Windowing & Other Spatial Functions            Implementation	**
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

#ifndef   WindowFct_cc_			// Is file already included?
#  define WindowFct_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// this is the implementation
#  endif

#include <Level1/WindowFct.h>		// Include the header file
#include <Matrix/row_vector.h>		// Include row vectors
#include <Basics/Gconstants.h>		// Include constant PI
#include <cmath>			// Include Bessel functions (j0)
#include <stdlib.h>


using std::cout;			// Using libstdc++ standard output

// ____________________________________________________________________________
// ____________________________ TRANSFORMS ______________________________
// ____________________________________________________________________________

void exponential_multiply (col_vector &db, double em, int offset)
  {
  col_vector vx = db;
  int s = db.elements();
  double s1 = s-1;
  for (int i=0; i<s; i++)
    vx.put(vx(i)*exp(em*abs(i-offset)/s1), i);
  db = vx;
  }

void exponential_multiply (row_vector &db, double em, int offset)
  {
  row_vector vx = db;
  int s = db.elements();
  double s1 = s-1;
  for (int i=0; i<s; i++)
    vx.put(vx(i)*exp(em*abs(i-offset)/s1), i);
  db = vx;
  }

// ______________________________________________________________________
// _________________________ WINDOW FUNCTIONS ___________________________
// ______________________________________________________________________


row_vector exponential(int size, int offset, double alpha)

	// Input	size   : data block size
	//		offset : function offset
	//		alpha  : line broadening parameter
	// Output	BLK    : data block containing an exponential 
	//			 function having maximum = 1 at the
	//			 offset point.  The decay rate is
	//			 specified by alpha as follows -
	//
	//			 alpha            application
	//			 _____   ______________________________
	//		      1.  >0     linewidth = 2ln2*alpha [1.3863]
	//		      2.   0     last point intensity exp(-4)
	//		      3. [-1,0)  last point intensity |alpha|
	//		      4.  <1	 same as >0 using |alpha|

  {
  row_vector BLK(size);			// Assume BLK is zeroed
  double fact = 0;
  if(alpha == 0)			// 2. set alpha so last pt. is exp(-4)
    alpha = double(size-1)/4.0;
  else if((alpha<0) && (alpha >=-1))	// 3. set alpha so last pt. has
    alpha = double(size-1) 		// intensity of |alpha|
		/log(fabs(alpha));
  fact = -1.0/fabs(alpha);		// This will cover 1. & 4. also
  for (int i=0; i<size; i++)		// Compute exponential
    BLK.put(exp(fact*abs(i-offset)),i);
  return BLK;
  }

	// Input	size	: Number of points
	// 		offset	: Zero offset (in points)
	// 		sigma	: Width parameter
	// Return	BLK	: 1D-data block containing a
	//			  discrete Gaussian function

row_vector Gaussian (int size, int offset, double sigma)
  {
  row_vector BLK(size);				// Assume BLK is zeroed
  double x = 0.;
  double denom = 2*sigma*sigma;
  for (int i=0; i<size; i++)
    {
    x = double(i-offset);			// Maximum=1 at point offset
    BLK.put(exp(-x*x/denom),i);		// 2.35*sigma is half-height linewidth
    }
  return BLK;
  } 

row_vector Hanning (int size, int offset)
  {
  row_vector BLK(size);
  double x = 0.;
  double fact = PI/double(size-1);
  for (int i=0; i<size; i++)
    {
    x = double(i-offset)*fact;	
    BLK.put(0.5 + 0.5*cos(x), i);		// Hanning function
    }
  return BLK;
  } 

row_vector Hamming (int size, int offset)
{
  row_vector BLK(size);				// Assume BLK is zeroed
  double x = 0.;
  double fact = PI/double(size-1);
  for (int i=0; i<size; i++)
    {
    x = double(i-offset)*fact;
    BLK.put(0.54 + 0.46*cos(x), i);		// Hamming function
    }
  return BLK;
} 

row_vector hyperbol_sec(int size, int offset, double alpha)
  {
  row_vector BLK(size);				// Assume BLK is zeroed
  double x = 0.;
  for (int i=0; i<size; i++)
    {
    x = double(i-offset)/alpha;
    BLK.put(1.0/cosh(x), i);			// Hyperbolic Secant function
    }
  return BLK;
  } 

row_vector Kaiser (int size, double theta, int offset)
  {
  row_vector BLK(size);				// Assume BLK is zeroed
  double x = 0.;
  double fact = 1.0/double(size-1);
  for (int i=0; i<size; i++)
    {
    x = double(i-offset)*fact;
    x = 1 - x*x;
    x = theta*sqrt(x);
#ifdef _MSC_VER
    BLK.put(_j0(x)/_j0(theta), i);		// Kaiser function
#else
    BLK.put(j0(x)/j0(theta), i);		  // Kaiser function
#endif
    }						// Uses Math Lib.
  return BLK;					// zero order Bessel
  } 

row_vector Lorentzian (int size, int offset, double alpha)
  {
  row_vector BLK(size);				// Assume BLK is zeroed
  double x = 0.;
  double alpsq = alpha*alpha;
  for (int i=0; i<size; i++)
    {
    x = double(i-offset);			// Maximum=1 at point offset
    BLK.put(alpsq/(alpsq+x*x),i);		// 2*alpha is half-height linewidth
    }
  return BLK;
  } 

row_vector sinc (int size, int offset, int inc)
 {
  row_vector BLK(size);				// Assume BLK is zeroed
  double x = 0.;
  int i;
  double fact = PI/double(inc);
  for (i=0; i<offset; i++)
    {
    x = double(i-offset)*fact;			 // Node inc points from offset
    BLK.put(sin(x)/x, i);
    }
  if ((offset>=0)&&(offset<size))                // takes care of sinc(0)
    BLK.put(1, offset);
  for (i=offset+1; i<size; i++)
    {
    x = double(i-offset)*fact;			 // Node inc points from offset
    BLK.put(sin(x)/x, i);
    }
  return BLK;
  } 

row_vector sin_square (int size, int offset)
  {
  row_vector BLK(size);
  double x = 0.;
  double fact = PI/double(size-offset-1);
  for (int i=0; i<size; i++)
    {
    x = double(i-offset)*fact;
    BLK.put(sin(x)*sin(x), i);
    }
  return BLK;
  }

row_vector square_wave(int size, int start, int finish)
  {
  row_vector BLK(size);				// Assume BLK is zeroed
  if(start < 0 || start >= finish)		// Check indices
    cout << "Invalid Initial Index for Square Wave Function";
  if(finish >= size)
    cout << "Invalid Final Index for Square Wave Function";
  for (int i=start; i<=finish; i++)
    BLK.put(1.0, i);
  return BLK;
  }

// ____________________________________________________________________________
//                            Random Noise Functions
// ____________________________________________________________________________

/* These functions deal with random noise.  They either return row vectors
   filled with random noise or add random noise to an existing row vector.
   The user must provide the MAXIMUM noise value, maxN (as a + number). The
   return will have all random noise between the values [-maxN, maxN]        */

row_vector Noise(int npts, double maxN)
  {
  row_vector data(npts, complex0);		// Construct a data vector
  Noise(data, maxN);				// Fill with random noise
  return data;
  }

void Noise(row_vector& data, double maxN)
  {
  double n1, n2;				// Real & Imag noise values
  for(int i=0; i<data.size(); i++)		// Loop over the vector points
    {	
    n1=(2.*(double(rand())/RAND_MAX)-1.0)*maxN;	// Get a noise value	
    n2=(2.*(double(rand())/RAND_MAX)-1.0)*maxN;	// Another noise value
    data.put(data.get(i)+complex(n1,n2), i);	// Add in (complex) noise
    }
  }

#endif 						// WindowFct.cc
