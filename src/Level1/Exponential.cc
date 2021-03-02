/* Exponential.cc ***********************************************-*-c++-*-
**                                                                      **
**                                 G A M M A                            **
**                                                                      **
**     Exponential Related Functions               Implementation	**
**                                                                      **
**      Scott A. Smith                                                  **
**      Copyright (c) 1994, 1997                                        **
**      National High Magnetic Field Laboratory                         ** 
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 ** 
**      $Header: $
**                                                                      **
**									**
*************************************************************************/

/*************************************************************************
**                                                                      **
** This file contains the implementation of GAMMA functions treating    **
** exponentials.  Function herein include standard exponentials as well	**
** as their derivatives (which are used in ESR).			**
**                                                                      **
*************************************************************************/

#ifndef   Gexponential_cc_		// Is this file already included?
#  define Gexponential_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Level1/Exponential.h>		// Include the header
#include <Basics/Gconstants.h>		// Get some conversion factors
#include <Matrix/matrix.h>		// Include matrices
#include <Matrix/row_vector.h>		// Include row vectors

// ____________________________________________________________________________
// A                    Straight Exponential Functions
// ____________________________________________________________________________

/* The analog exponential function in GAMMA is

                                         i*(w+iR)t
   [1]                           E(t) = e
  
   where w is the oscillation frequency and R is the decay rate, both in 1/sec.
   This exponential is a Fourier transform pair with a complex Lorentzian that
   is centered at freqency w and has a half-height linewidth of 2*R.  Note that
   unless R is non-negative the exponential will blow up at large t.

   Since we can't output anything but a discrete representation of this function 
   we must know the number of points and the total time the exponential lasts
   for proper scaling.  For point k, k=[0,npts), we have 

                      i*(w+iR)t
   [2]           E = e         k     where      t = k*(time/npts-1) = k*tinc
                  k                              k

   Alternatively, we could input our w and R values in 1/pt.  The equation then
   is
                                   i*(w+iR)k
   [3]                        E = e         
                               k
  
   The problem with this form is that we have to suffer in relating it to 
   the Lorenzian produced from its Fourier transform.  The transform likes 
   things in radians as well and we'll need conversion factors if we wish to
   use [3] with respect to it's Lorentzian (Fourier Transform Pair).  So, if
   users need a quick picture of a decaying exponential use the simple function
   that deals in points.  For most purposes, you'll be better off using the
   overloaded form in which you explicitly state your input units.           */
 
row_vector Exponential(int npts, double W, double R)
 
        // Input        npts    : Block size in pts
        //              W       : Exponential frequency (1/pt)
        //              R	: Decay rate (1/pt)
        // Output       data    : Exponential with oscillation frequency w
	//			  and decay rate R
	// Note			: We have no relationship between time and
	//			  point increment, so all input in pt & 1/pt
	//			  i.e. implements [3] above.
 
  {
  row_vector data(npts);                        // Data block for return
  complex expon(-R,W);				// Exponent term i(w+iR)
  complex efact = exp(expon);			// Exponential increment
  complex evalue = complex1; 			// Start with value of 1
  data.put(evalue, 0);				// Initial point is 1
  for(int i=1; i<npts; i++)			// Loop through points
    {
    evalue *= efact;				// Increment in complex plane
    data.put(evalue, i);			// Store this value
    }
  return data;
  }  

row_vector Exponential(int npts, double time, double w, double RT, int type)
 
        // Input        npts    : Block size in pts
	//		time    : Length of the exponential (sec)
        //              w       : Exponential frequency 
        //              RT      : Decay rate, decay time
        //              units   : Flag for input units & type
	//			      0 : w=1/sec R=1/sec
	//			      1 : w=Hz, T=sec
	//			      2 : w=cycles/point R=1/point
        // Output       data    : Exponential with oscillation frequency w
	//			  and decay rate R
	// Note			: Since we know time, we use equation [2]
	//			  above and adjust {w,RT} to be rad/sec
 
  {
  row_vector data(npts);                        // Data block for return
  double W, R;					// Freqeuncy, decay rate, 
  switch(type)
    {
    case 0:
    default:					// Input: rad/sec & 1/sec
      W = w;					// Already 1/sec
      R = RT;					// Already 1/sec
      break;
    case 1:					// Input: Hz & sec
      W = w*HZ2RAD;				// W now in radians/sec
      R = 1/RT;					// R now in 1/sec
      break;
    case 2:					// Input: cycles/pt & decay/pt
      time = double(npts-1);			// Set time so axis is points
      W = w*HZ2RAD;				// W now in radians/point
      R = RT;					// R now in decay/point
      break;
    }
  double tinc = time/double(npts-1);		// Time increment per point
  complex expon = complex(-R,W);		// Exponent term i(w+iR)
  complex efact = exp(expon*tinc);		// Exponential increment
  complex evalue = complex1; 			// Start with value of 1
  data.put(evalue, 0);				// Initial point is 1
  for(int i=1; i<npts; i++)			// Loop through points
    {
    evalue *= efact;				// Increment in complex plane
    data.put(evalue, i);			// Store this value
    }
  return data;
  }  


// ____________________________________________________________________________
// B                    Exponential Derivative Functions
// ____________________________________________________________________________

/* The analog derivative exponential function in GAMMA is

                          d  [  i(W+iR)t]        i(W+iR)t
   [1]            E'(t) = -- [ e        ] = -it*e
                          dW
  
   where the derivative is taken with respect to the frequency W.  Both the
   oscillation frequency W and the decay rate R will typically have units 1/sec
   This derivative exponential is a Fourier transform pair with a complex
   derivative Lorentzian. Whereas the derivative of the normal exponential will
   produce a Lorentzian centered at frequency W with a half-height linewidth
   of 2*R, the differential exponential's transform will produce a differential
   Lorentzian that crosses the axis at W and whose peak-to-peak linewidth is
   given by 2*R/sqrt(3).  Note that unless R is non-negative the derivative
   exponential will blow up at large t.

   Since we can't output anything but a discrete representation of this function
   we must know the number of points and the total time the differential 
   exponential lasts for proper scaling.  For point k, k=[0,npts), we have
 
                          i(W+iR)t
   [2]        E' = -it * e        k  where      t = k*(time/npts-1) = k*tinc
               k      k                          k

   Alternatively, we could input our W and R values in 1/pt.  The equation then
   is
                                         i*(w+iR)k
   [3]                        E = -ik * e         
                               k      
  
   The problem with this form is that we have to suffer in relating it to 
   the Lorenzian produced from its Fourier transform.  The transform likes 
   things in radians as well and we'll need conversion factors if we wish to
   use [3] with respect to it's Lorentzian (Fourier Transform Pair).  So, if
   users need a quick picture of a decaying differential exponential use the
   simple function that deals in points.  For most purposes, you'll be better
   off using the overloaded form in which you explicitly state input units.  */
 
row_vector DExponential(int npts, double W, double R)
 
        // Input        npts    : Block size in pts
        //              W       : Exponential frequency (1/pt)
        //              R	: Decay rate (1/pt)
        // Output       data    : Derivative exponential with oscillation
	//			  frequency W and decay rate R
	// Note			: We have no relationship between time and
	//			  point increment, so all input in pt & 1/pt
	//			  i.e. implements [3] above
 
  {
  row_vector data(npts);                        // Data block for return
  complex expon(-R,W);				// Exponent term i(w+iR)
  complex efact = exp(expon);			// Exponential increment
  complex evalue = complex1; 			// Start with value of 0
  data.put(complex0, 0);			// Initial point is 1
  for(int i=1; i<npts; i++)			// Loop through points
    {
    evalue *= efact;				// Increment in complex plane
    data.put(complex(0,-double(i))*evalue, i);	// Store this value
    }
  return data;
  }  

 
row_vector DExponential(int npts, double time, double w, double RT, int type)
 
        // Input        npts    : Block size in pts
        //              time    : Length of the exponential (sec) 
        //              w       : Exponential frequency
        //              RT      : Decay rate, decay time
        //              typ     : Flag for input units & type
        //                            0 : w=1/sec R=1/sec
        //                            1 : w=Hz, T=sec
        //                            2 : w=cycles/point R=1/point
        // Output       data    : Derivative exponential with oscillation
        //                        frequency w and decay rate R
        // Note                 : Since we know time, we use equation [2]
        //                        above and adjust {w,RT} to be rad/sec

  {
  row_vector data(npts);                        // Data block for return
  double W, R;					// Freqeuncy, decay rate, 
  switch(type)
    {
    case 0:
    default:					// Input: rad/sec & 1/sec
      W = w;					// Already 1/sec
      R = RT;					// Already 1/sec
      break;
    case 1:					// Input: Hz & sec
      W = w*HZ2RAD;				// W now in radians/sec
      R = 1/RT;					// R now in 1/sec
      break;
    case 2:					// Input: cycles/pt & decay/pt
      time = double(npts-1);			// Set time so axis is points
      W = w*HZ2RAD;				// W now in radians/point
      R = RT;					// R now in decay/point
      break;
    }
  double tinc = time/double(npts-1);		// Time increment per point
  complex expon = complex(-R,W);		// Exponent term i(w+iR)
  complex efact = exp(expon*tinc);		// Exponential increment
  complex evalue = complex1; 			// Start @ exponent val of 1
  data.put(0, 0);					// Initial point is 0
  double t = 0;					// Initial time is zero
  for(int i=1; i<npts; i++)			// Loop through points
    {
    t -= tinc;					// Time * -1 at this point
    evalue *= efact;				// Increment in complex plane
    data.put(complex(0,t)*evalue, i);		// Store this value
    }
  return data;
  }  


// ____________________________________________________________________________
// C                      Exponential Cutoff Functions
// ____________________________________________________________________________

//                       k       >= -ln(cutoff)/[R*tinc]
//                        cutoff 


        // Input        npts    : Block size in pts
	//		time    : Length of the exponential (sec)
        //              w       : Exponential frequency 
        //              R       : Decay rate
	//		cutoff	: Cutoff factor, range is [0, 1]
	// Output	k	: The point at which the exponential defined
	//		 cut	  by {npts,time,w,R} falls below cutoff

int Exponen_cut(int npts, double time, double w, double R, double cutoff)
  {
  if(cutoff>1 || cutoff<0) cutoff = 1.e-4; 	// Insure cutoff reasonable
  else if(!cutoff)         return npts;		// If no cutoff, last point
  else if(cutoff == 1)     return 0;		// If cutoff 1, 1st point
  double tinc = time/double(npts-1);		// Time increment per point
  double tcutoff = -log(cutoff)/R;		// Cutoff time desired
  return int(ceil(tcutoff/tinc));		// Cutoff point
  }


void Exponen_cut(int* ihi,const matrix& mx,double tinc,int npts,double cutoff)

	// Input	ihi   : Array of upper indices
	//      	mx    : Array of transitions.  For transition i
	//		        R = Re(<i|mx|0>) and w = Im(<i|mx|0>
	//			with both quantities in (1/sec)
	//		tinc  : Time between points (sec)
	//		npts  : Number of points allowed
	//		cutoff: Cutoff factor, range is [0, 1]
	// Output	void  : The points in ihi are set.  At points >ihi
	//			the exponential magnitude value will be below
	//			cutoff. 

 {
 if(cutoff>1 || cutoff<1.e-10) cutoff = 1.e-4;	// Insure cutoff = (0,1]
 double mod = -log(cutoff)/tinc;
 double R;
 int ntr = mx.rows();				// Number of transitions
 for(int i=0, tr=0; tr<ntr; tr++)		// Loop transiitons
   {
   R = mx.getRe(tr,0);				// Get the decay rate
   if(R <=0) ihi[tr] = npts;			// No decay, all points O.K.
   else						// If point does decay, find
     {						// where intensity < cutoff
     i = int(mod/R);				// Point where cutoff is
     if(i>npts-1) ihi[tr] = npts;		// Insure point within range
     else         ihi[tr] = i+1;		// or adjust to next point 
     }
   }
 return;
 }

/*************************************************************************
**			Mathematical Details				**
*************************************************************************/

// A complex exponential is given by

//                                i(w+iR)t
//                        E(t) = e

// where w is the oscillation frequency (in 1/sec) and R is the decay
// rate (in 1/sec).  R will be a non-negative number or else the function
// will blow up for large t.

// The magnitude of a decaying complex exponential is given by

//                                 -Rt
//                       |E(t)| = e

// We seek a time where the magnitude reaches the value cutoff where
// cutoff spans [0.0,1.0].

//                  -Rt
//        cutoff = e     --> t = -ln(cutoff)/R

// For the discrete function, the time is that at which the value of
// cutoff is reached is given by (at for i=0, t=0)

//                      t = i * t
//                               inc

// Thus the discrete equation will be roughly

//                  i = -ln(cutoff)/Rt

// Since i is integer, it is best to avoid problems from double
// roundoff and at one to the computed value

#endif

