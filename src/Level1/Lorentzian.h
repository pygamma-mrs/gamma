/* Lorentzian.h *************************************************-*-c++-*-
**                                                                      **
**                                 G A M M A                            **
**                                                                      **
**        Lorentzian                               Interface		**
**                                                                      **
**        Copyright (c) 1994                                            **
**        Dr. Scott A. Smith                                            **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** This file contains the implementation of GAMMA functions treating    **
** Lorentzians.  Function herein include standard Lorentzians as well   **
** as their derivatives.                                                **
**                                                                      **
*************************************************************************/

#ifndef GLorentzian_h_			// Is this file already included?
#  define GLorentzian_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <string>			// Include libstdc++ strings
#include <Matrix/complex.h>		// Include GAMMA complex numbers
#include <Matrix/row_vector.h>		// Include row vectors
#include <Matrix/matrix.h>		// Include matrices

// ____________________________________________________________________________
// A                    Straight Lorentzian Functions
// ____________________________________________________________________________ 
 
/* The analog complex Lorentzian used in GAMMA is given by (lwhh = 2*R)

                             R - i*(w-w )
                                       o         1
                      L(w) = ------------ = -----------
                              2         2   R + i(w-w )
                             R  + (w-w )             o   
                                      o

   which has the discrete form, using w  = w + k*w
                                       k    0     inc

                    R - i*(w -w )
                            k  o           1                 1   
            <L|k> = -------------- = ------------- -->  -------------
                     2           2   R + i(w - w )      Rpt + i(k-k )
                    R  + (w - w )           k   o                  o
                           k   o
 
   The last form (on the right) has both R and W expressed in points.        */

MSVCDLL complex Lorentzian(double Weval, double W, double R);

        // Input        Weval   : Frequency at evaluation
        //              W       : Lorentzian offset frequency
        //              R       : Lorentzian linewidth factor
        // Output       L(weval): Lorentzian evaluated at Weval


 
MSVCDLL row_vector Lorentzian(int npts, double Wi, double Wf, double W, double R);
 
        // Input        npts    : Block size in pts
        //              Wi      : Initial frequency
        //              Wf      : Final frequency
        //              W       : Lorentzian offset frequency
        //              R       : Lorentzian linewidth factor
        // Output       data    : Lorentzian centered about offset with a
        //                        half-height linewidth of 2*R
 
 
MSVCDLL row_vector Lorentzian(int npts, double W, double R, int max1=0);
 
        // Input        npts    : Block size in pts
        //              W	: Lorentzian offset in pts
        //              R	: Linewidth factor in pts
	//		max1    : Flag for normalization
        // Output       data    : Lorentzian centered about offset with a
        //                        half-height linewidth of 2*R
	// Note			: An extra R is put into the numerator 
	//			  for Lmax=1 if max1 != 0
 

// sosi - this is currently in HSproc.cc but maybe it should be move to
//        this file instead?
//row_vector Lorentzian(int size, int offset=0, double alpha=1.0);

	// Input	size  : Number of points
	//	 	offset: Offset point
	//		alpha : Width Factor
	// Output	BLK   : A 1D data block containing a
	//			complex Lorentzian function

//                  alpha + i*(k-offset)
//           BLK  = --------------------
//              k        2             2
//                  alpha  + (k-offset)


// ____________________________________________________________________________
// B                       Fancy Lorentzian Functions
// ____________________________________________________________________________ 


MSVCDLL row_vector Lorentzian(double R, double W,
           double wst, double wfi, int N, complex& I, double Lcut, int Lint);

        // Input        R     : Lorentzian rate (rad/sec)
        //              W     : Lorentzian frequency (rad/sec)
        //              wst   : Initial frequency (rad/sec)
        //              wfi   : Final frequency (rad/sec)
        //              N     : Number of points
        //              I     : General intensity scaling
        //              Lcut  : Decimal for maximum value cutoff
        //              Lint  : Flag for integrated intensities
        //                       0 - Don't use integrated values
        //                      !0 - Use integrated value, this
        //                           value is # points required to  
        //                           span the half-height linewidth
        // Output       Lvect : A row_vector containing a complex
        //                      discrete Lorentzian function


MSVCDLL void Lorentzian(row_vector& data, double R, double W,
                  double wst, double wfi, complex& I, double Lcut=1.e-4, int Lint=5);

	// Input	data  : A row vector
	//	 	R     : Lorentzian rate
	//		W     : Lorentzian frequency
	//		wst   : Initial frequency
	//		wfi   : Final frequency
	//		I     : General intensity scaling
	//		Lcut  : Decimal for maximum value cutoff
	//		Lint  : Flag for integrated intensities
	//			 0 - Don't use integrated values
	//			!0 - Use integrated value, this
	//			     value is # points required to  
	//			     span the half-height linewidth
	// Output	void  : A row_vector containing a complex
	//			discrete Lorentzian function
	// Note		      : Switching {wst,wfi} order switches
	//			the sign on the imaginary component
	// Note		      : Units on R, W, wst, and wfi should
	//			all match up (normally radians/sec)

//                                     R - i*(w - W)
//         R - i*(w-W)                         i                     wfi-wst
//  L(w) = -----------  ---> data(i) = -------------  ; w  = wst + i -------
//          2        2                  2          2     i           npts-1
//         R  + (w-W)                  R  + (w - W)
//                                             i

MSVCDLL void Lorentzian(int* Lint, const matrix& mx, double winc, int N=5);

	// Input	Lint  : Array of flags (signaling poor resolution)
	//      	mx    : Array of transitions.  For transition i
	//		        R = Re(<i|mx|0>) and w = Im(<i|mx|0>
	//			with both quantities in (1/sec)
	//		winc  : Frequency between points (1/sec)
	//		N     : Number of points acceptable per lwhh
	// Output	void  : The array Lint if filled with 0 or 1
	//			based on how many points will be used
	//			to characterize a Lorentzian function.
	//			A non-zero value implies that an integrated
	//			or average value should be used rather than
	//			the function at that particular point because
	//			the function is being poorly sampled
	// Note		      : The test is based on how many points
	//			span the Lorentzian linewidth at half-height.
	//			For a well sampled function, it is assumed
	//			           N*winc < lwhh
	//			i.e. that there will be at least N points
	//			sampled before a single linewidth is spanned.
	// Note			The default value of N is 5. A lower value
	//			may produce poor looking discrete functions but
	//			in a faster computation.  A higher value produces
	//			a more well defined discrete function but can
	//			be more CPU intensive.

// The magnitude of a complex Lorentzian is given by

//                                      R - i*(w - w0)
//         R - i*(w-w0)                         i
//  L(w) = ------------  ---> data(i) = --------------  ; w  = wst + i winc
//          2         2                  2           2     i
//         R  + (w-w0)                  R  + (w - w0)
//                                             i

// This Lorentzian form relates to it's use in NMR according to the following.

//                               R     1
//                       lwhh = -- = -----
//                              PI   PI*T
//                                       2
//
// where R is in radians/sec and T in sec and lwhh in Hz. Switching lwhh to
// seconds implies that lw = 2R

// Herein, the inequality  

//                       lwhh = 2R > N*winc

// for all transitions in matrix mx.  If it fails, then the corresponding
// integration flag is set to non-zero, meaning the resolution will be poor.


// ____________________________________________________________________________
// C                       Lorentzian Cutoff Functions
// ____________________________________________________________________________

// These functions test Lorentzian baselines, i.e. where the Lorentzian points
// are indistinguishable from baseline.  If the Lorentzian linewidth is narrow
// relative to the point spread most points will be essentially zero and need
// not be computed.  

	// Input	ilo   : Low pt or low pt array of cutoffs
	//		ihi   : High pt of high pt array of cutoffs
	//      	R     : A Lorentzian rate or matrix whose 1st column
	//			contains Lorentzian rates (1/sec)
	//	 	W     : A frequency (1/sec).  In the case of array
	//			input the 1st column imaginaries of R will
	//			be taken as frequencies.
	//		w0    : Frequency at 1st point (1/sec)
	//		winc  : Frequency increment (1/sec)
	//		npts  : Number of points allowed
	//		npts  : Number of points allowed
	//		cutoff: Cutoff factor, range is [0, 1]

// Note that these function DO NOT return Lorentzians.  They only return
// point indices that mark where the Lorentizan values fall below a specified
// cutoff value.  These flags may be subsequently used in other Lorentzian
// functions (section A) in determining the fastest way to make discrete
// Lorentzians.  Lorentzian points below ilo(s) and above ihi(s) need not
// be computed as the Lorentzian magnitude will be below the value

//			      cutoff * (1/R)

// Switching the order of wst & wfi switches the sign on the imaginary 
// component! Units on R, wo, wst, and wfi must all match up, normally 1/sec.
// The mathematical details of this cutoff are give at the end of this file. 

MSVCDLL void Lorentz_cut(int& ilo, int& ihi, double R, double W,
                     double w0, double winc, int npts, double cutoff=1.e-4);

MSVCDLL void Lorentz_cut(int* ilo, int* ihi, const matrix& R,
                    double w0, double winc, int npts, double cutoff=1.e-4);

// ____________________________________________________________________________
// D                     Lorentzian Integration Functions
// ____________________________________________________________________________

// These functions test Lorentzian "resolution", i.e. how well a Lorentzian
// is characterized by the discrete points.  If few points are taken or the
// Lorentzian is very broad it may be that the discrete function poorly
// represents the function.

	// Input	Lint  : A flag or array of flags which signal
	//			Lorentzian  poor digital resolution
	//      	R     : A Lorentzian rate or matrix whose 1st column
	//			contains Lorentzian rates (1/sec)
	//		winc  : Frequency between points (1/sec)
	//		N     : Number of points acceptable per lwhh
	// Output	void  : Lint value(s) set to 0 or 1 based
	//			on whether the discrete Lorentzian function
	//			is sufficiently digitized.
	//			 0=Well digitized; 1=Poorly digitized

// Note that these function DO NOT return Lorentzians.  They only return
// flags regarding Lorentzian digitization.  These flags may be subsequently
// used in other Lorentzian functions (section A) in determining the best way 
// to make discrete Lorentzians (section A).  The mathematical description is
// given at the end of this file.

// The test is based on how many points span the Lorentzian linewidth at
// half-height.  For a well sampled function its assumed that the condition 

//		           N*winc < lwhh = 2*R

// is satisfied, i.e. that there will be at least N points sampled before a
// single linewidth is spanned.  The default value of N is 5 implying that at
// least 5 should characterize the peak.  A lower value will allow for poorer
// digitized functions to be deemed acceptable, but this may result in faster
// Lorentzian computations if the returned flags are used to decide how GAMMA
// Lorentzians are generated.  In contrast, a higher N value will make the
// flags require well digitized functions (vs. linewidth) but may imply more
// CPU intensive calculations when used to decide computation routes.

MSVCDLL void Lorentz_int(int& Lint, double R, double winc, int N=5);

MSVCDLL void Lorentz_int(int* Lint, const matrix& R, double winc, int N=5);

// ____________________________________________________________________________
// E               Lorentzian Interactive Setup Functions
// ____________________________________________________________________________

 
MSVCDLL void ask_Lorentzian(int argc, char* argv[], int& qn, int& N,
                  double& wst, double& wfi, double& W, double& R,
                                                 double& fact, int& pplw);

        // Input        argc    : Number of argments
        //              argv    : Array of arguments
        //              qn      : Argument index
        //              N       : Number of points
        //              wst	: Initial Lorentzian frequency
        //              wfi	: Final Lorentzian frequency
	//		W       : Lorentzian offset
	//		R       : Lorentzian linewidth
        //              fact    : Cutoff percent value
        //              pplw    : Minimum allowed points per linewidth
        // Output       void    : The values N, wst, wfi, W, R,
	//			  fact and pplw are set either interactively
	//			  or by the values in argv.

 
MSVCDLL void read_Lorentzian(const std::string& filein, int& N, double& wst, double& wfi,
                 double& W, double& R, double& fact, int& pplw, int idx=-1);

        // Input        filein  : A filename
        //              N       : Number of points
        //              wst	: Initial Lorentzian frequency
        //              wfi	: Final Lorentzian frequency
	//		W       : Lorentzian offset
	//		R       : Lorentzian linewidth
        //              fact    : Cutoff percent value
        //              pplw    : Minimum allowed points per linewidth
	//		idx	: Lorentzian pulse index (default none)
        // Output       void    : The values N, val1, val2, and fact
        //                        are set from values in file filein


// ____________________________________________________________________________
// F                    Lorentzian Derivative Functions
// ____________________________________________________________________________

/* The analog form of the Differential Lorentzian used in GAMMA are given by
 
                                                  -i
                          DL(w) = L'(w) =  ----------------
                                                         2
                                           [R + i(w - w ) ]
                                                       o
 
   which has the discrete form, using w  = w + k*w
                                       k    0     inc
 
                                     -i                    -i
                      <DL|k> = ---------------- ---> ---------------
                                              2                    2
                               [R  + (w - w )]       [R  + i(k-k )]
                                       k   o           pt       o
 
   The last form (on the right) has both R and W expressed in points.        */


MSVCDLL complex DLorentzian(double Weval, double W, double R);
 
        // Input        Weval   : Frequency at evaluation
        //              W       : Lorentzian offset frequency
        //              R       : Lorentzian linewidth factor
        // Output       L(weval): Lorentzian evaluated at Weval


MSVCDLL row_vector DLorentzian(int npts, double Wi, double Wf, double W, double R);

        // Input        npts    : Block size in pts
        //              Wi      : Initial frequency
        //              Wf      : Final frequency
        //              W       : Lorentzian offset frequency
        //              R       : Lorentzian linewidth factor
        // Output       data    : Differential of a complex Lorentzian
        //                        spanning frequencies [Wi,Wf]
        //                        centered about W with lwhh of 2R

MSVCDLL row_vector DLorentzian(int npts, double offset, double R, int max1=0);
 
        // Input        npts    : Block size in pts
        //              W       : Lorentzian offset in pts
        //              R	: Linewidth factor in pts
	//		max1    : Flag for normalization
        // Output       data    : A differential Lorentzian taken from a
        //                        Lorentzian centered about W with a
        //                        half-height linewidth of 2*R
	// Note			: An extra R*R is put into the numerator 
	//			  for L'max=1 if max1 != 0
 
/*************************************************************************
**									**
**                    Additional Mathematical Details			**
**									**
*************************************************************************/

// ----------------------------------------------------------------------------
// ----------------------- Lorentzian Cutoff Functions ------------------------
// ----------------------------------------------------------------------------

/* The magnitude of a complex Lorentzian is given by

                   |   2       2  |1/2   |             |1/2
                   |  R + (w-W)   |      |      1      |
          |L(w)| = | -----------  |    = | ----------- |
                   |[ 2        2]2|      |  2        2 |
                   |[R  + (w-W) ] |      | R  + (w-W)  |

   where W is the frequency of maximum intensity. We seek a frequency
   w where the magnitude reaches the value cutoff/R. Since the largest
   Lorentzian magnitude is 1/R, the input for cutoff should [0,1], if
   if cutoff >1 no points will fit the request.

             2                                |     2        |
       cutoff        1                        |    R       2 |1/2
       ------ = -----------   OR  (w-W) = +/- | ------- - R  |
          2      2        2                   |       2      |
         R      R  + (w-W)                    | cutoff       |

   Since cutoff <= 1, the argument in the square root is always positive.
   Evidently, whenever

                            |    1        |1/2
                w = W +/- R | ------- - 1 |
                            |       2     |
                            | cutoff      |

   the Lorentzian value is at cutoff/R, or at (100xcutoff)% of the peak
   maximum.  Note that if cutoff = 1 then only at w=W will we have a
   a solution and for cutoff = 0 all points are solutions.

   A final adjustment is necessary because the concern herein is for a
   discrete, not continuous, function. Thus, the frequency at point k is
    given by (where k spans [0, npts])

                            w  = w  + k*w                
                             k    0      inc

   Thus the discrete equation will be
     
                                        |    1        |1/2
                  k*w    = (W-w ) +/- R | ------- - 1 |
                     inc       0        |       2     |
                                        | cutoff      |

   So, the discrete Lorentzian magnitude is roughly at value cutoff/R when

                    1  [                |    1        |1/2 ]
              k = ---- | (W - w ) +/- R | ------- - 1 |    |
                  w    |       0        |       2     |    |
                   inc [                | cutoff      |    ]

   which makes a bit more sense if we use K as the "index" at frequency W,
   and delk as the adjustment from there

                              k = K +/- delk

   It is only rough because k is an integer and the r.h.s. of this equation
   will not be.  Since the Lorentzian is symmetric, it is easier to see the
   range over which the Lorentzian is above the set cutoff by noting that
   the Lorentzian magnitude maximizes at
 
                              W - W
                                   0
                    k    =    ------
                     max       w 
                                inc

   and hence the cutoff level is reached when
 
                           |    1        |1/2
                       R   |             |
          k = k   +/- ___  | ------- - 1 | = k    +/- delk
              max     w    |       2     |    max
                       inc | cutoff      |

   As two k values result and both are integer values in the computer, the
   lower of the two values will need to be decremented and the upper of the
   two has its value incremented.  (The reverse is true if winc is < 0).  
   Again, if cutoff is 1 then only 1 point will suffice, the integer k is 
   int([W-w0]/winc).  Note that there is nothing to hinder a result which has 
   i be negative or larger than is desired, so this is compensated for in
   order to keep the range between [0, npts].				     */

// ----------------------------------------------------------------------------
// -------------------- Lorentzian Integration Functions ----------------------
// ----------------------------------------------------------------------------

/* The magnitude of a complex Lorentzian is given by

                                      R - i*(w - w0)
         R - i*(w-w0)                         i
  L(w) = ------------  ---> data(i) = --------------  ; w  = wst + i winc
          2         2                  2           2     i
         R  + (w-w0)                  R  + (w - w0)
                                             i

 This Lorentzian form relates to it's use in NMR according to the following.

                               R     1
                       lwhh = -- = -----
                              PI   PI*T
                                       2
                     -1                                   -1
 where R & T  are sec  , lwhh in Hz. Switching lwhh to sec  implies lw = 2R
            2

 Herein, the inequality  

                       lwhh = 2R > N*winc

 for all transitions in matrix mx.  If it fails, then the corresponding
 integration flag is set to non-zero, meaning the resolution will be poor.    */

#endif							// Lorentizian.h

