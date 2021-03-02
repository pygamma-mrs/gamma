/* Lorentzian.cc ************************************************-*-c++-*-
**									**
**                                 G A M M A                            **
**									**
**      Lorentzian Related Functions 	           Implementation       **
**									**
**      Scott A. Smith							**
**      Copyright (c) 1994, 1997					**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** This file contains the implementation of GAMMA functions treating	**
** Lorentzians.  Function herein include standard Lorentzians as well	**
** as their derivatives (which are used in ESR).			**
**									**
*************************************************************************/

#ifndef GLorentzian_cc_			// Is this file already included?
#define GLorentzian_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Level1/Lorentzian.h>		// Include the header
#include <Matrix/matrix.h>		// Include matrices
#include <Matrix/row_vector.h>		// Include row vectors
#include <Basics/Gutils.h>		// Include query parameter
#include <Basics/Gconstants.h>		// Include constant PI
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/StringCut.h>		// Include GAMMA string manipulation
#include <string>			// Include libstdc++ strings
#include <list>				// Include libstdc++ lists

using std::string;			// Using libstdc++ strings
using std::cout;			// Using libstdc++ cout function
using std::list;			// Using libstdc++ STL lists

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


        // Input        Weval	: Frequency at evaluation
        //              W	: Lorentzian offset frequency
        //              R	: Lorentzian linewidth factor
        // Output       L(weval): Lorentzian evaluated at Weval

complex Lorentzian(double Weval, double W, double R)
  { return 1.0/complex(R, Weval-W); }


        // Input        npts    : Block size in pts
        //              Wi	: Initial frequency
        //              Wf	: Final frequency
        //              W	: Lorentzian offset frequency
        //              R	: Lorentzian linewidth factor
        // Output       data    : Lorentzian spanning frequencies [Wi,Wf]
	//			  centered about W with lwhh of 2R

row_vector Lorentzian(int npts, double Wi, double Wf, double W, double R)
  {
  row_vector data(npts);                        // Data block for return
  double Winc = (Wf-Wi)/double(npts-1);		// Frequency increment per pt
  double Wk = Wi;
  for(int k=0; k<npts; k++)
    {
    data.put(1.0/complex(R, Wk-W), k);
    Wk += Winc;
    }
  return data;
  }

        // Input        npts    : Block size in pts
        //              W	: Lorentzian offset in pts
        //              R	: Linewidth factor in pts
	//		max1	: Flag whether to normalize
        // Output       data    : Lorentzian centered about offset with a
        //                        half-height linewidth of 2*R
	// Note			: An extra R is put into the numerator 
	//			  for Lmax=1 if requested

row_vector Lorentzian(int npts, double W, double R, int max1)
  {
  row_vector data(npts);                        // Data block for return
  double num = 1;				// Set the numerator
  if(max1) num = R;				// (its R if max1 is set)
  for(double i=0; i<npts; i++)
    data.put(num/complex(R, i-W), int(i));
  return data;
  }
 

// sosi - this is currently in nmr_proc.cc but maybe it should be move to
//        this file instead?
//row_vector Lorentzian(int size, int offset=0, double alpha=1.0)

	// Input	size  : Number of points
	//	 	offset: Offset point
	//		alpha : Width Factor
	// Output	BLK   : A 1D data block containing a
	//			complex Lorentzian function

//                  alpha + i*(k-offset)
//           BLK  = --------------------
//              k        2             2
//                  alpha  + (k-offset)

//{
//  row_vector BLK(size);				// Assume BLK is zeroed
//  double x = 0.;
//  double alpsq = alpha*alpha;
//  for (int i=0; i<size; i++)
//    {
//    x = double(i-offset);			// Maximum=1 at point offset
//    BLK(i) = alpsq/(alpsq+x*x);			// 2*alpha is half-height linewidth
//    }
//  return BLK;
//} 


// ____________________________________________________________________________
// B                       Fancy Lorentzian Functions
// ____________________________________________________________________________

row_vector Lorentzian(double R, double W,
           double wst, double wfi, int N, complex& I, double Lcut, int Lint)

	// Input	R     : Lorentzian rate (rad/sec)
	//		W     : Lorentzian frequency (rad/sec)
	//		wst   : Initial frequency (rad/sec)
	//		wfi   : Final frequency (rad/sec)
	//		N     : Number of points
	//		I     : General intensity scaling
	//		Lcut  : Decimal for maximum value cutoff
	//		Lint  : Flag for integrated intensities
	//			 0 - Don't use integrated values
	//			!0 - Use integrated value, this
	//			     value is # points required to  
	//			     span the half-height linewidth
	// Output	Lvect : A row_vector containing a complex
	//			discrete Lorentzian function

  {
  row_vector Lvect(N);				// Vector for Lorentzian pts
  Lorentzian(Lvect,R,W,wst,wfi,I,Lcut,Lint); 	// Use function overload
  return Lvect;
  }


void Lorentzian(row_vector& data, double R, double W,
           double wst, double wfi, complex& I, double Lcut, int Lint)

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

/*                                     R - i*(w - W)
           R - i*(w-W)                         i                     wfi-wst
    L(w) = -----------  ---> data(i) = -------------  ; w  = wst + i -------
            2        2                  2          2     i           npts-1
           R  + (w-W)                  R  + (w - W)
                                               i                             */
  {
  int npts = data.size();			// Number of points in data
  double winc = (wfi-wst)/(npts-1);		// Freq. between pts (1/sec)
  int ist = 0;					// This is the first point
  int ifi = npts;				// This is (beyond) the last pt
//sosi
//  int ifi = npts;				// This is (beyond) the last pt
  if(Lcut) 					// If cutoff indicated, then
    Lorentz_cut(ist,ifi,R,W,wst,winc,npts,Lcut);// find computation pt range
  if(Lint) 					// If integrated values allowed
    {						// then find if integral needed
    int M = Lint;
    Lorentz_int(Lint, R, winc, M);
    }
  W = (wst * winc*double(ist)) - W;		// Starting frequency
  double Wa, Wb;				// Freqs for integrated values
  double re,im;					// Values for integration
  int i;					// Used for point index
  double R2 = R*R;				// Need the rate squared
  double W2 = W*W;				// Need frequency squared
  if(!Lint)
    for(i=ist;i<ifi; i++,W+=winc)		// Loop over contributing pts
      { 					// using Lorentzian values
      W2 = W*W;					// Use freq. squared in denom.
      data.put(I*complex(R,-W)/(R2+W2),i);
      }
  else
    for(i=ist;i<ifi; i++,W+=winc)		// Loop over contributing pts
      {						// using integrated values
      Wa = W - winc/2.0;			//	Pt 1/2 increment back
      Wb = W + winc/2.0;			//	Pt 1/2 increment ahead
      re = atan(Wb/R)-atan(Wa/R);		//	Integrated real part
      im = -0.5*log((R2+Wb*Wb)/(R2+Wa*Wa));	//	Integrated imag part
      data.put(I*complex(re/winc,im/winc),i);	//	Use average integrated
      } 					//	value 
  return;
  }


void Lorentzian(int* Lint, const matrix& mx, double winc, int N)

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
   integration flag is set to non-zero, meaning the resolution will be poor. */

  {
  int ntr = mx.rows();				// Number of frequencies
  double nw = double(N)*winc;			// Minimal frequency span
  for(int tr=0; tr<ntr; tr++)			// Loop over all frequencies
    { 						// If linewidth is less than
    if(2.0*mx.getRe(tr,0) < nw) 		// requested frequency span
      Lint[tr] = 1;				// then flag for integration
    else					// If not, then flag for
      Lint[tr] = 0;				// use of direct value
//cout << "\nTransition " << tr
//     << " Linewidth " << 2.0*mx.getRe(tr,0)/(2.0*PI) << ": "
//     << " Hz: Integral Flag = " << Lint[tr];
    }
  return;
  }

// ____________________________________________________________________________
// C                       Lorentzian Cutoff Functions
// ____________________________________________________________________________
 
// These functions test Lorentzian baselines, i.e. where the Lorentzian points
// are indistinguishable from baseline.  If the Lorentzian linewidth is narrow
// relative to the point spread most points will be essentially zero and need
// not be computed.
 
        // Input        ilo   : Low pt or low pt array of cutoffs
        //              ihi   : High pt of high pt array of cutoffs
        //              R     : A Lorentzian rate or matrix whose 1st column
        //                      contains Lorentzian rates (1/sec)
        //              W     : A frequency (1/sec).  In the case of array
        //                      input the 1st column imaginaries of R will
        //                      be taken as frequencies.
        //              w0    : Frequency at 1st point (1/sec)
        //              winc  : Frequency increment (1/sec)
        //              npts  : Number of points allowed
        //              npts  : Number of points allowed
        //              cutoff: Cutoff factor, range is [0, 1]
 
// Note that these function DO NOT return Lorentzians.  They only return
// point indices that mark where the Lorentizan values fall below a specified
// cutoff value.  These flags may be subsequently used in other Lorentzian
// functions (section A) in determining the fastest way to make discrete
// Lorentzians.  Lorentzian points below ilo(s) and above ihi(s) need not
// be computed as the Lorentzian magnitude will be below the value
 
//                            cutoff * (1/R)

// Switching the order of wst & wfi switches the sign on the imaginary
// component! Units on R, wo, wst, and wfi must all match up, normally 1/sec.
// Mathematical details of this are give at the end of the header file.


void Lorentz_cut(int& ilo, int& ihi, double R, double W,
                        double w0, double winc, int npts, double cutoff)

 {
 if(cutoff>1.0 || cutoff<1.e-10)		// Set a decent cutoff
   cutoff = 1.e-4;				// if none is specified
 double mod = R*sqrt(1.0/(cutoff*cutoff)-1.0);	// Modification from middle
 double delw = W-w0;				// Between 1st point & middle  
 double w1 = (1/winc)*(delw + mod);		// One extrema, (double)pts
 double w2 = (1/winc)*(delw - mod);		// Other extrema, (double)pts
//cout << "\nCutoff Level = " << cutoff;
//cout << "\nFrequency of First Point = " << w0/(2.0*PI);
//cout << "\nFrequency of Lorentzian Offset = " << W/(2.0*PI);
//cout << "\nFrequency Difference = " << delw/(2.0*PI);
//cout << "\nThe Linewidth is = " << 2*R/(2.0*PI);
//cout << "\nThe Point Modifier = " << mod << "(" << int(mod) << ")";
//cout << "\nw1 (Points as double) = " << w1;
//cout << "\nw2 (Points as double) = " << w2;
 if(w1 > w2)					// Here if w2 is low frequency end
   {
   if(w2 <= 0) w2 = 0;				// 	Insure low end !< 0
   else if(w2>npts) w2=npts+1;			//	Insure low end !> npts
   else w2 -= 1;				//	Drop a point for safety

   if(w1 <= 0) w1 = 0;				// 	Insure high end !<0
   else if(w1>npts) w1=npts+1;			//	Insure hign end !> npts
   else w1 += 1;				//	Add a point for safety

   ilo = int(w2);				// 	Set low point as integer
   ihi = int(w1);				//	Set high point as integer
   }
 else
   {
   if(w1 <= 0) w1 = 0;				// 	Insure low end !<0
   else if(w1>npts) w1=npts+1;			//	Insure low end !> npts
   else w1 -= 1;				//	Drop a point for safety

   if(w2 <= 0) w2 = 0;				// 	Insure high end !<0
   else if(w2>npts) w2=npts+1;			//	Insure high end !> npts
   else w2 += 1;				//	Add a point for safety

   ilo = int(w1);				//	Set low point as integer
   ihi = int(w2);				//	Set high point a integer
   }
 return;
 }


void Lorentz_cut(int* ilo, int* ihi, const matrix& mx,
                   double w0, double winc, int npts, double cutoff)

 {
 if(cutoff>1.0 || cutoff<1.e-10) cutoff = 1.e-4;// Insure a decent cutoff
 double sqrtmod = sqrt(1.0/(cutoff*cutoff)-1.0);// Modification from middle
 double mod, delw;
 int kmax, delk, k1, k2;

 int ntr = mx.rows();				// Number of components
 for(int tr=0; tr<ntr; tr++)			// Loop over all transitions
   {
   mod = mx.getRe(tr,0)*sqrtmod;		// Pt modification from middle
   delw = mx.getIm(tr,0) - w0;			// Freq. difference from start
   kmax = int(delw/winc);			// Point at Lorentzian maximum
   delk = int(mod/winc);			// Point adjustment to cutoff
   if(delk>0) delk++;
   else delk--;
   k1 = kmax - delk;				// First point extrema
   k2 = kmax + delk;				// Second point extrema
   if(k1 <= 0) k1 = 0;				// Insure k1 !< 0
   else if(k1>=npts-1) k1=npts-1;		// Insure k1 !> npts
   if(k2 <= 0) k2 = 0;				// Insure k2 !< 0
   else if(k2>=npts-1) k2=npts-1;		// Insure k1 !> npts
   if(k1 > k2)					// If k2 is low frequency
     {
     ilo[tr] = k2;				// Set k2 as low point
     ihi[tr] = k1;				// Set k1 as high point
     }
   else						// If k1 is low frequency
     {
     ilo[tr] = k1;				// Set k1 as low point
     ihi[tr] = k2;				// Set k2 as high point
     }
   if(ilo[tr] >= 1) ilo[tr] = ilo[tr]-1;	// Expand extrema to account
   if(ihi[tr] <= npts-1) ihi[tr] = ihi[tr]+1;	// for any double->int roundoff
//cout << "\nTransition " << tr
//     << " at " << mx.getIm(tr,0)/(2.0*PI) << ": "
//     << ilo[tr] << ", " << ihi[tr];
   }
//cout << "\n\n";
 return;
 }


// ____________________________________________________________________________
// D                     Lorentzian Integration Functions
// ____________________________________________________________________________

/* These functions test Lorentzian "resolution", i.e. how well a Lorentzian
   is characterized by the discrete points.  If few points are taken or the
   Lorentzian is very broad it may be that the discrete function poorly
   represents the function.

           Input        Lint  : A flag or array of flags which signal
                                Lorentzian  poor digital resolution
                        R     : A Lorentzian rate or matrix whose 1st column
                                contains Lorentzian rates (1/sec)
                        winc  : Frequency between points (1/sec)
                        N     : Number of points acceptable per lwhh
           Output       void  : Lint value(s) set to 0 or 1 based
                                on whether the discrete Lorentzian function
                                is sufficiently digitized.
                                 0=Well digitized; 1=Poorly digitized
 
   Note that these function DO NOT return Lorentzians.  They only return
   flags regarding Lorentzian digitization.  These flags may be subsequently
   used in other Lorentzian functions (section A) in determining the best way
   to make discrete Lorentzians (section A).  The mathematical description is
   given at the end of the header file.
 
   The test is based on how many points span the Lorentzian linewidth at
   half-height.  For a well sampled function its assumed that the condition
 
                           N*winc < lwhh = 2*R
 
   is satisfied, i.e. that there will be at least N points sampled before a
   single linewidth is spanned.  The default value of N is 5 implying that at
   least 5 should characterize the peak.  A lower value will allow for poorer
   digitized functions to be deemed acceptable, but this may result in faster
   Lorentzian computations if the returned flags are used to decide how GAMMA
   Lorentzians are generated.  In contrast, a higher N value will make the
   flags require well digitized functions (vs. linewidth) but may imply more
   CPU intensive calculations when used to decide computation routes.        */

void Lorentz_int(int& Lint, double R, double winc, int N)
  {
  double nw = double(N)*winc;		// Compute the N*winc value
  if(2.0*R < nw) Lint = 1;		// If 2R<N*winc, should integrate
  else Lint = 0;			// If not, integration not needed
  }


void Lorentz_int(int* Lint, const matrix& mx, double winc, int N)
  {
  int ntr = mx.rows();			// Number of Lorentzian functions
  double nw = double(N)*winc;		// Compute the N*winc value
  for(int tr=0; tr<ntr; tr++)		// Loop over all Lorentzians
    {
    if(2.0*Re(mx.get(tr,0)) < nw)	//	Test if 2*R < N*winc
      Lint[tr] = 1;			//	If so, should integrate
    else Lint[tr] = 0;			//	If not, not needed
    }
  }

// ____________________________________________________________________________ 
// E                Lorentzian Interactive Setup Functions
// ____________________________________________________________________________ 

 
  void ask_Lorentzian(int argc, char* argv[], int& qn, int& N,
                  double& wst, double& wfi, double& W, double& R,
                                                 double& fact, int& pplw)

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

  {
  query_parameter(argc, argv, qn++,		// Get number of steps(pts)
       "\n\tNumber of Points in Lorentzian? ", N);
  query_parameter(argc, argv, qn++,             // Plot initial frequency
       "\n\n\tInitial Frequency(Hz)? ", wst);
  query_parameter(argc, argv, qn++,             // Plot final frequency
       "\n\n\tFinal Frequency(Hz)? ", wfi);
  double lw;
  query_parameter(argc, argv, qn++,             // Lorentzian frequency
      "\n\n\tLorentzian Frequency Offset(Hz)? ", W);
  query_parameter(argc, argv, qn++,             // Lorentzian linewidth
      "\n\n\tLorentzian Linewidth (Hz)? ", lw);
  R = lw/2.0;
  query_parameter(argc, argv, qn++,             // Minimum Points per width
       "\n\n\tMinimum Allowed Points per Linewidth? ", pplw);
  query_parameter(argc, argv, qn++,             // Intensity cutoff
       "\n\n\tPercent of Maximum Cutoff Value [0.0, 1.0]? ", fact);
  return;
  } 

 
void read_Lorentzian(const string& filein, int& N, double& wst, double& wfi,
                       double& W, double& R, double& fact, int& pplw, int idx)

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


  {
  ParameterSet pset;				// A parameter set
  pset.read(filein);				// Read pset in
  ParameterSet::const_iterator item;		// A pix into the pset
  SinglePar par;				// A single parameter
  string pname, ssfile, pstate;			// Items in each pset entry
  string SI = string("(") + Gdec(idx)		// Name adjustment if indexed
            + string(")");

//		   Read The Number of Lorentzian Steps

  pname = string("LzPts");			// Number of Lorentzian steps
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  par = SinglePar(pname);
  item = pset.seek(par);			// Pix in parameter list
  if(item != pset.end()) 			// Get number of points
    (*item).parse(pname,N,pstate);	 	// if the parameter is found
  else
    {
    cout << "\n\tCant Find " << pname << " in "
         << filein;
    cout << "\n\tSetting 1000 Lorentzian Points";
    N = 1000;
    }

//		   Read The Initial Frequency (Hz)

  pname = string("LzWst");			// Initial frequency
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  par = SinglePar(pname);
  item = pset.seek(par);			// Pix in parameter list
  if(item != pset.end()) 			// Get initial frequency
    (*item).parse(pname,wst,pstate);	 	// if the parameter is found
  else
    {
    cout << "\n\tCant Find " << pname << " in "
         << filein;
    cout << "\n\tSetting Initial Frequency to 0";
    wst = 0.0;
    }

//		   Read The Final Frequency (Hz)

  pname = string("LzWfi");			// Final frequency
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  par = SinglePar(pname);
  item = pset.seek(par);			// Pix in parameter list
  if(item != pset.end()) 			// Get final frequency
    (*item).parse(pname,wfi,pstate);	 	// if the parameter is found
  else
    {
    cout << "\n\tCant Find " << pname << " in "
         << filein;
    cout << "\n\tSetting Final Frequency to 1000";
    wst = 1000.0;
    }

//		   Read The Frequency Offset (Hz)

  pname = string("LzWo");			// Frequency offset
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  par = SinglePar(pname);
  item = pset.seek(par);			// Pix in parameter list
  if(item != pset.end()) 			// Get frequency offset
    (*item).parse(pname,W,pstate);	 	// if the parameter is found
  else
    {
    cout << "\n\tCant Find " << pname << " in "
         << filein;
    cout << "\n\tSetting Frequency to 500";
    W = 500.0;
    }

//              Read One of Two: {Line Width, Rate}

// GAMMA's Lorentzian relationship between the linewidth at half-height (lwhh)
// and its R value (a decay rate of its related exponential) according to

//                    R     1
//          lwhh   = -- = -----  ---> lwhh  = 2R   &  lwhh     = 2R
//              Hz   PI   PI*T            Hz    Hz        1/sec    1/sec
//                            2
//                     -1
// where R & T  are sec  , lwhh in Hz. Switching both lwhh & R into the same
//            2                        units implies the relationship lw = 2R.

  double lwhh;					// Linewidth value (Hz)
  pname = string("Lzlw");			// Line width parameter
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  par = SinglePar(pname);
  item = pset.seek(par);			// Pix in parameter list
  if(item != pset.end()) 			// Get half-height linewidth
    { 						// if the parameter is found
    (*item).parse(pname,lwhh,pstate);		// This is set in Hz
    R = PI*lwhh;				// Set R value in 1/sec
    }
  else
    {
    pname = string("LzR");			// Rate parameter (1/sec)
    if(idx >= 0) pname += SI;			// Adjust name if indexed			
    par = SinglePar(pname);
    item = pset.seek(par);			// Pix in parameter list
    if(item != pset.end()) 			// Get associated decay rate
      { 					// if the parameter is found
      (*item).parse(pname,R,pstate);		// This is set in (1/sec)
      }
    else
      {
      cout << "\n\tCant Find " << pname << " in "
         << filein;
      R = 25.0;
      cout << "\n\tSetting Linewidth to "
           << R/PI << " Hz";
      }
    }

//                     Read Cutoff Intensity as %

  pname = string("Lzcut");			// Lorentzian cutoff
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  par = SinglePar(pname);
  item = pset.seek(par);			// Pix in parameter list
  if(item != pset.end()) 			// Get intensity cutoff parameter
    (*item).parse(pname,fact,pstate);		// Set as a % [0-100] of maximum
  else
    {
    cout << "\n\tCant Find " << pname << " in "
         << filein;
    cout << "\n\tSetting Lorentzian Cutoff to 2%";
    fact = 0.02;
    }

//                    Read Lorentzian Resolution

  pname = string("Lzpplw");			// Lorentzian resolution
  if(idx >= 0) pname += SI;			// Adjust name if indexed			
  par = SinglePar(pname);
  item = pset.seek(par);			// Pix in parameter list
  if(item != pset.end()) 			// Get Lorentzian resolution
    (*item).parse(pname,pplw,pstate);
  else
    {
    cout << "\n\tCant Find " << pname << " in "
         << filein;
    cout << "\n\tSetting Lorentzian Resolution for"
         << " a Minimum of 5 Points Per Linewidth";
    pplw = 5;
    }
  }


// ____________________________________________________________________________ 
// F                    Lorentzian Derivative Functions
// ____________________________________________________________________________

/* The analog form of the Differential Lorentzian used in GAMMA are given by

                                                  -i
                          DL(w) = L'(w) =  ----------------
                                                          2
                                           [R + i(w - w )]
                                                       o

   which has the discrete form, using w  = w + k*w
                                       k    0     inc

                                     -i                    -i
                      <DL|k> = ---------------- ---> ---------------
                                              2                    2
                               [R  + (w - w )]       [R  + i(k-k )]
                                       k   o           pt       o

   The last form (on the right) has both R and W expressed in points.        */


complex DLorentzian(double Weval, double W, double R)

        // Input        Weval	: Frequency at evaluation
        //              W	: Lorentzian offset frequency
        //              R	: Lorentzian linewidth factor
        // Output       L(weval): Lorentzian evaluated at Weval

  { complex z(R,Weval-W); return -complexi/(z*z); }


row_vector DLorentzian(int npts, double Wi, double Wf, double W, double R)
 
        // Input        npts    : Block size in pts
        //              Wi      : Initial frequency
        //              Wf      : Final frequency
        //              W       : Lorentzian offset frequency
        //              R       : Lorentzian linewidth factor
        // Output       data    : Differential of a complex Lorentzian
        // 			  spanning frequencies [Wi,Wf]
	//			  centered about W with lwhh of 2R
 
  {
  row_vector data(npts);                        // Data block for return
  double Winc = (Wf-Wi)/double(npts-1);         // Frequency increment per pt
  double delW = Wi-W;
  complex numer(0,-1);
  complex denom;
  for(int k=0; k<npts; k++)
    {
    denom = complex(R, delW);
    data.put(numer/(denom*denom), k);
    delW += Winc;
    }
  return data;
  }  
 

row_vector DLorentzian(int npts, double W, double R, int max1)

        // Input        npts    : Block size in pts
        //              W	: Lorentzian offset in pts
        //              R 	: Linewidth factor in pts
	//		max1    : Flag whether to normalize
        // Output       data	: A differential Lorentzian taken from a
	//			  Lorentzian centered about W with a
	//			  half-height linewidth of 2*R
	// Note			: An extra R*R in the numerator is used 
	//			  to set L'max = -i  
 
  {
  row_vector data(npts);			// Data block for return
  double x;
  complex denom;
  complex mR2i = -complexi;
  if(max1) mR2i *= (R*R);
  for(double i=0; i<npts; i++)
    {
    x = i-W;					// Maximum=1 at point offset
    denom = complex(R, x);
    data.put(mR2i/(denom*denom), int(i));
    }
  return data;
  }    

/*************************************************************************
**			Mathematical Details				**
*************************************************************************/

// ----------------------------------------------------------------------
// -------------------- Lorentzian Cutoff Functions ---------------------
// ----------------------------------------------------------------------

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

                             w = w  + k*w                
                                  0      inc

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

   It is only rough because k is an integer and the r.h.s. of this equation
   will not be.  As two k values result, the lower of the two has its value
   decremented and the upper of the two has its value incremented.  (The
   reverse is true if winc is < 0).  Again, if cutoff is 1 then only 1 point
   will suffice, the integer closest to (W-w0)/winc.  Note that there is no-
   thing to hinder a result which has i be negative or larger than is desired,
   so this is compensated for in order to keep the range between [0, npts]. */


// ----------------------------------------------------------------------
// ----------------- Lorentzian Integration Functions -------------------
// ----------------------------------------------------------------------

/* The magnitude of a complex Lorentzian is given by

                                      R - i*(w - w0)
         R - i*(w-w0)                         i
  L(w) = ------------  ---> data(i) = --------------  ; w  = wst + i winc
          2         2                  2           2     i
         R  + (w-w0)                  R  + (w - w0)
                                               i

  This Lorentzian form relates to it's use in NMR via to the following:

                                   R     1
                           lwhh = -- = -----
                                  PI   PI*T
                                           2
                       -1                                  -1
  where R & T  are sec  , lwhh in Hz. Switching lwhh to sec   ==> lw = 2R
             2

  Herein, the inequality  
                            lwhh = 2R > N*winc

  for all transitions in matrix mx. If it fails, then the corresponding
  integration flag is set !0, meaning the resolution will be poor.      */

#endif							// Lorentzian.cc
       
