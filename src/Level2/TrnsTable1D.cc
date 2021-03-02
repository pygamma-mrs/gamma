/* TrnsTable.cc *************************************************-*-c++-*-
**									**
** 	                       G A M M A				**
**									**
**	Class Transitions Table 1D		Implementation		**
**									**
**	Copyright (c) 1999				 		**
**	Dr. Scott A. Smith				 		**
**      National High Magnetic Field Laboratory 	                **
**      1800 E. Paul Dirac Drive                        	        **
**      Tallahassee Florida, 32306-4005                         	**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** This file contains the implementation of class TrnsTable1D.  The 	**
** purpose of this class is to embody a digitial representation of a	**
** 1D nmr/epr spectrum.  In this instance the spectrum is stored as a	**
** ntr x 2 complex array where ntr is the number of transtions.  Each	**
** row corresponds to a transition.  The first column contains both the	**
** transtition frequency (imaginaries, in 1/sec) and the transition	**
** linewidths (reals, 1/sec). The second column contains the transition	**
** intensities. The intensity is stored complex (a+ib --> Rexp[iphi])	**
** in order to accomodate to handle phased lines.			**
**									**
** Note that although one can use this class independently, it is by &	**
** large just an auxiliary class to the GAMMA acquisition classes. That	**
** implies that most users won't ever deal with these transition arrays	**
** directly in a GAMMA program.						**
**									**
*************************************************************************/

#ifndef _TrnsTab1D_cc_			// Is this file already included?
#  define _TrnsTab1D_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Level2/TrnsTable1D.h>		// Include the header
#include <Basics/Gutils.h>		// Know GAMMA error messages
#include <Basics/Gconstants.h>		// Need to know HZ2GAUSS & PI
#include <iostream>			// Include filestreams (read/write)
#include <Level1/Lorentzian.h>		// Include Lorentzians
#include <Level1/Exponential.h>		// Include Exponentials
#include <Level1/WindowFct.h>		// Include Random Noise
#include <Matrix/matrix.h>		// Need to know arrays
#include <Matrix/row_vector.h>		// Need to know vectors 
#include <vector>			// Include libstdc++ STL vectors
#include <Basics/StringCut.h>		// Include GAMMA string manipulations
#include <stdlib.h>
#include <cmath>			// Include HUGE_VAL_VAL

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ STL lists
using std::cout;			// Using libstdc++ standard output
using std::vector;			// Using libstdc++ STL vectors
using std::ostream;			// Using libstdc++ output streams
using std::ofstream;			// Using libstdc++ output file streams
using std::ifstream;			// Using libstdc++ input file streams
using std::ios;

bool TTable1D::FRQREV = false;		// Default Reverse Frequencies Flag

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                 CLASS TRANSITIONS TABLE 1D ERROR HANDLING
// ____________________________________________________________________________
 
/*      Input	                TTab1D  : A transitions table (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void TTable1D::TTaberror(int eidx, int noret) const
  {
  string hdr("1D Transition Table");
  string msg;
  switch(eidx)
    {
    case 49: msg=string("Frequency Output Type Unknown");		//(49)
             GAMMAerror(hdr,msg,noret); break;
    case 50: msg=string("Wrong Number of Columns in Table");		//(50)
             GAMMAerror(hdr,msg,noret); break;
    case 51: msg=string("Transition Linewidths Cannot Be Negative");	//(51)
             GAMMAerror(hdr,msg,noret); break;
    case 52: msg=string("Setting Transition Linewidth To Zero");	//(52)
             GAMMAerror(hdr,msg,noret); break;
    case 60: msg=string("Can't Determine # of Transitions");		//(60)
             GAMMAerror(hdr,msg,noret); break;
    case 61: msg=string("Should Set TrN in ASCII Parameter File");	//(61)
             GAMMAerror(hdr,msg,noret); break;
    default: GAMMAerror(hdr, eidx, noret); break;
    }
  }

void TTable1D::TTaberror(int eidx, const string& pname, int noret) const
  {
  string hdr("1D Transition Table");
  string msg;
  switch(eidx)
    {
    case 4:
    case 5: GAMMAerror(hdr, msg+pname, noret); break;
    case 49: msg=string("Frequency Output Type Unknown");		//(49)
             GAMMAerror(hdr,msg,noret); break;
    case 50: msg=string("Cannot Read Table From ");			//(50)
             GAMMAerror(hdr, msg+pname, noret); break;
    case 53: msg=string("Linewidth Value Of ")+pname+string("/sec");	//(53)
             GAMMAerror(hdr, msg, noret); break;
    case 60: msg=string("Cannot Determine Transition ");		//(60)
             GAMMAerror(hdr, msg+pname, noret); break;
    default: GAMMAerror(hdr,eidx,pname,noret); break;
    }  
  }  

volatile void TTable1D::TTabfatality(int eidx) const
  {
  TTaberror(eidx, 1);
  if(eidx) TTaberror(0);
  GAMMAfatal();					// Clean exit from program
  }
    
volatile void TTable1D::TTabfatality(int eidx, const string& pname) const
  {
  TTaberror(eidx, pname, 1);
  if(eidx) TTaberror(0);
  GAMMAfatal();					// Clean exit from program
  }


// ____________________________________________________________________________
// ii              CLASS TRANSITIONS TABLE 1D CHECKING FUNCTIONS
// ____________________________________________________________________________

	// Input	TTab1D: A transitions table (this)
	//		typ   : A frequency output type
	//				0 = Hertz (default)
	//				1 = PPM
	//				2 = Gauss
	//		warn  : A warning flag
	// Output	TF    : TRUE if value of typ is known

int TTable1D::CheckType(int typ, int warn) const
  {
  if(typ<0 || typ>2)
    {
    if(warn)
      {
      if(warn>1) TTabfatality(49);		// Unknown freq. output type
      else       TTaberror(49);
      }
    return 0;
    }
  return 1;
  }

// ____________________________________________________________________________
// iii              CLASS TRANSITIONS TABLE COMMON CODE
// ____________________________________________________________________________

	// Input		TTab1D: A transitions table (this)
	//			TTab1 : A second transitions table
	// Output (Defaults)	void  : This sets the default internals
	// Output (Copy)	void  : Copies internal settings from TT1

void TTable1D::setDefaults()
  {
  ICUT    = 1.e-12;			// Set transition intensity cutoff
  INORM   = 1;				// Set transition normalization
  PERCUT  = 0.001;			// Set transition intensity % cutoff
  FRQTYPE = 0;				// Set for output in Hz
  FRQSORT = 1;				// Set for frequency sorted output
  FRQCONV = 1.0;			// Set frequency output conversion fact 
  SN  =  0;				// Set the signal to noise (0=no noise)
  PHP =  0;				// Set phases to not print
  LWP =  0;				// Set linewidths to print as needed
  RP  = -1;				// Set R values to print as needed
  T2P = -1;				// Set T2 values to not print 
  HP = 1;				// Set Header to print 
  }

void TTable1D::copySettings(const TTable1D& TTab1)
  {
  ICUT    = TTab1.ICUT;			// Copy transition intensity cutoff
  INORM   = TTab1.INORM;		// Copy transition normaliation
  PERCUT  = TTab1.PERCUT;		// Copy transition intensity % cutoff
  FRQTYPE = TTab1.FRQTYPE;		// Copy frequency output type flag
  FRQSORT = TTab1.FRQSORT;		// Copy frequency output sort flag
  FRQCONV = TTab1.FRQCONV;		// Copy frequency output conv. factor
  SN      = TTab1.SN;			// Copy signal to noise ratio
  RP      = TTab1.RP;			// Copy R values print flag
  HP      = TTab1.HP;			// Copy Header print flag
  LWP     = TTab1.LWP;			// Copy linewidths print flag
  T2P     = TTab1.T2P;			// Copy T2 values print flag
  PHP     = TTab1.PHP;			// Copy phases print flag
  }

// ____________________________________________________________________________
// iv                  TRANSITIONS TABLE 1D SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions set up specific aspects of a transitions table.  Since       
   they make assumptions about the order in which the table is set up, the
   functions MUST be private, their misuse may make an inconsistent table.

   Parameter Name(s)    Optional            Function and Priority
 =====================  ========  ============================================
          TrN             No      The number of transitions in the table
   TrOmega, Trgfact       Yes     The conversion factor for PPM or Gauss input
  TrLW,TrT2,TrR,TrPP      No      The transition decay rate
          TrI             No      The transition intensity
          TrPh            Yes     The transition phase
 TrHz, TrW, TrPPM, TrB    No      The transition frequency

   Note that specification of the conversion factor, either TrOmega or Trgfact,
   will limit the valid frequency parameter names. If TrOmega has been set then
   we can take PPM values for input frequencies. If Trgfact has been specified
   then we can take Gauss values for input frequencies.  If neither of these
   has been set then we can ONLY take frequencies in Hz or 1/sec (TrHz, TrW)
   This is because we always store values in 1/sec and must be able to convert.

   0. 1/SEC ===> W : no conversion needed                 
   1. Hz ======> W : cycles/sec *     1      * PIx2 radians/cycle
   2. PPM =====> W : ppm        * Om Hz/ppm  * PIx2 radians/cycle
   3. GAUSS ===> W : Gauss      * GAUSS2HZ   * PIx2 radians/cycle

   So it is only for PPM and GAUSS frequency input that we must know some
   conversion factor. For the former this will be a spectrometer frequency in
   MHz (e.g 500 MHz), and that must of course be based upon the nuclear spin
   type whose transitions we are representing.  For the latter we need to know
   how to convert from GAUSS to Hz.  The formula for this is g*beta*H/hbar = w,
   so one needs an appropriate g-factor since H is specified and both beta and 
   hbar are constants.

   Function   Parameters                     What It Does
 ------------ ---------- ------------------------------------------------------
   SetNTrans    TrN      Sets the number of transitions. Must be present.
    SetConv    TrOmega   Sets conversion factor, FRQCONV, for PPM/Gauss input.
               Trgfact   This will also set the frequency type flag, FRQTYPE.
   GetLWR2T   TrLW, TrT2 Set transition decay rate. May be rate R, the lwhh,
              TrR, TrPP  the relaxation time T2, or peak-2-peak width, lwpp
   GetPhase     TrPh     Sets the transition phase (default is zero)
 GetIntensity   TrI      Sets the transition intensity, this must be present.
 GetFrequency TrHz, TrW  Sets the transition frequency, this must be present.
              TrPPM, TrB Can be set in Hz, 1/sec, PPM or Gauss.
   SetTrans     ---      Sets all of the transitions, one by one, using above

	   Input		TTab1D	: A transitions table (this)
                                pset	: A parameter set
				warn    : Warning level flag
                                                0 = no warning
                                                1 = warning, non-fatal
                                               >1 = fatal
           Output               TF	: Boolean whether we worked or not   */

int TTable1D::SetNTrans(const ParameterSet& pset, int warn)
  {
  string pstate;
  int npts;
  ParameterSet::const_iterator item;		// Pix in parameter list
  string pname = string("TrN");			// # transitions par name
  item = pset.seek(pname);                      // Try & find TrN par
  if(item != pset.end())                        // Retrieve the number of
    { 						// transitions and set the
    (*item).parse(pname,npts,pstate);		// table size for this
    matrix mx(npts, 2, complex0);
    *this = TTable1D(mx);
    return 1;
    }
  if(warn)					// Cant read TrN
    {
    TTaberror(60, 1);				// Can't read # transitions
    TTaberror(61, 1);				// Set Trn in parameter set
    if(warn > 1) TTabfatality(9);		// Problems during construction
    else         TTaberror(9);  		
    }
  return 0;
  }  
 
int TTable1D::SetConv(const ParameterSet& pset, int warn)
  {
  string pstate;
  ParameterSet::const_iterator item;		// Pix in parameter list

  string pname = string("TrOmega");             // Transitions base freq.
  item = pset.seek(pname);                      // Try & find TrOmega par
  if(item != pset.end())                        // Retrieve the base freq.
    {                                           // of transitions and set
    (*item).parse(pname,FRQCONV,pstate);        // it as FRQCONV
    FRQTYPE = 1;				// Frequency type is PPM
    return 1;
    }

  pname = string("Trgfac");			// Transitions g factor
  item = pset.seek(pname);			// Try & find Trgfac par
  if(item != pset.end()) 			// Retrieve the g factor
    {						// of transitions and set
    (*item).parse(pname,FRQCONV,pstate);        // it as FRQCONV
    FRQTYPE = 2;				// Frequency type is Gauss
    return 1;
    }
  return 0;
  }
 
double TTable1D::GetLWR2T2(const ParameterSet& pset, int idx, int warn)
  {
  string pstate;
  double pval = 0;
  ParameterSet::const_iterator item;		// Pix in parameter list
  string pidx=string("(")+Gdec(idx)+string(")");// Transition index
  string pname = string("TrLW") + pidx;		// Transition linewidth
  item = pset.seek(pname);                      // Try & find TrLW par
  if(item != pset.end())                        // Retrieve the Linewidth
    {						// which will be input in
    (*item).parse(pname,pval,pstate);		// Hz, so we also convert
    pval *= PI;					// it to rad/sec
    return pval;
    }
  pname = string("TrT2") + pidx; 		// Transition T2 value
  item = pset.seek(pname);                      // Try & find TrLW par
  if(item != pset.end())                        // Retrieve the T2 value
    {						// which will be input in
    (*item).parse(pname,pval,pstate);		// sec, so we also convert
    pval = 1/pval;				// it to rad/sec
    return pval;
    }
  pname = string("TrR") + pidx; 		// Transition R value
  item = pset.seek(pname);                      // Try & find TrR par
  if(item != pset.end())                        // Retrieve the R value
    {
    (*item).parse(pname,pval,pstate); 		// which will in rad/sec
    return pval;				// (or whatever units)
    }
  pname = string("TrPP") + pidx; 		// Transition peak2peak lw
  item = pset.seek(pname);                      // Try & find TrPP par
  if(item != pset.end())                        // Retrieve the lwpp value
    {
    (*item).parse(pname,pval,pstate); 		// which will input in Gauss,
    pval *= sqrt(0.75);				// so we convert to rad/sec
    pval *= PIx2/HZ2GAUSS; 			// R = sqrt(3/4)*lwpp
    }
  return pval;	
  }
 
double TTable1D::GetPhase(const ParameterSet& pset, int idx, int warn)
  {
  string pstate;
  double pval = 0;
  ParameterSet::const_iterator item;         // Pix in parameter list
  string pidx=string("(")+Gdec(idx)+string(")");// Transition index
  string pname = string("TrPh") + pidx;		// Transition phase 
  item = pset.seek(pname);                      // Try & find TrPh par
  if(item != pset.end())                        // Retrieve the phase
    (*item).parse(pname,pval,pstate);
  return pval;
  }

double TTable1D::GetIntensity(const ParameterSet& pset, int idx, int warn)
  {
  string pstate;
  double pval = 0;
  ParameterSet::const_iterator item;         // Pix in parameter list
  string pidx=string("(")+Gdec(idx)+string(")");// Transition index
  string pname = string("TrI") + pidx;		// Transition intensity
  item = pset.seek(pname);                      // Try & find TrI par
  if(item != pset.end())                        // Retrieve intenstity
    (*item).parse(pname,pval,pstate);
  else
    {
    if(warn)					// Cant read TrI
      {
      TTaberror(2, pname, 1);			// Can't read parameter pname
      if(warn > 1) TTabfatality(9);		// Problems during construction
      }
    }
  return pval;
  }

double TTable1D::GetFreq(const ParameterSet& pset, int& idx, int warn)
  {
  string pstate;				// Dummy for statement
  string pname;					// For the par name
  double pval = 0;				// For the value
  ParameterSet::const_iterator item;		// Pix in parameter list
  string pidx=string("(")+Gdec(idx)+string(")");// Transition index
  switch(FRQTYPE)
    {
    case 0:					// Here if in Hz or 1/sec
    default:					// or any other units
      {
      string pname = string("TrHz") + pidx;	// Transition freq name
      item = pset.seek(pname);			// Try & find TrHz par
      if(item != pset.end())			// Retrieve the freq
        { 					// in Hz if possible
        (*item).parse(pname,pval,pstate);
        pval *= HZ2RAD;	
        return pval;
        }
      pname = string("TrW") + pidx;		// Transition freq name
      item = pset.seek(pname);			// Try & find TrW par
      if(item != pset.end())			// Retrieve the freq
        { 					// in rad/sec if possible
        (*item).parse(pname,pval,pstate);
        return pval;
        }
      idx = -1;					// Flag we have failed
      if(warn)					// Cant read TrHz
        {
        pname = string("TrHz") + pidx;		// Transition freq name
        TTaberror(2, pname, 1);			// Can't read parameter pname
        if(warn > 1) TTabfatality(9);		// Problems during construction
        }
      }
      break;
    case 1:					// Here if in PPM
      {
      string pname = string("TrPPM") + pidx;	// Transition phase 
      item = pset.seek(pname);			// Try & find TrPPM par
      if(item != pset.end())			// Retrieve the freq
        { 					// in PPM if possible
        (*item).parse(pname,pval,pstate);	// (but return in rad/sec)
        return pval*FRQCONV;
        }
      pname = string("TrHz") + pidx;		// Transition freq name
      item = pset.seek(pname);			// Try & find TrHz par
      if(item != pset.end())			// Retrieve the freq
        { 					// in Hz if possible
        (*item).parse(pname,pval,pstate);
        pval *= HZ2RAD;	
        return pval;
        }
      pname = string("TrW") + pidx;		// Transition freq name
      item = pset.seek(pname);			// Try & find TrW par
      if(item != pset.end())			// Retrieve the freq
        { 					// in rad/sec if possible
        (*item).parse(pname,pval,pstate);
        return pval;
        }
      idx = -1;					// Flag we have failed
      if(warn)					// Cant read TrHz
        {
        pname = string("TrPPM") + pidx;		// Transition freq name
        TTaberror(2, pname, 1);			// Can't read parameter pname
        if(warn > 1) TTabfatality(9);		// Problems during construction
        }
      }
      break;
    case 2:					// Here if in Gauss
      {
      string pname = string("TrB") + pidx;	// Transition field
      item = pset.seek(pname);			// Try & find TrB par
      if(item != pset.end())			// Retrieve the field
        { 					// in Gauss if possible
        (*item).parse(pname,pval,pstate);	// but we'll convert the
        return (pval/HZ2GAUSS)*PIx2;		// result to rad/sec
        }
      idx = -1;					// Flag we have failed
      if(warn)					// Cant read TrB
        {
        TTaberror(2, pname, 1);			// Can't read parameter pname
        if(warn > 1) TTabfatality(9);		// Problems during construction
        }
      }
      break;
    }
  return pval;
  }

int TTable1D::SetTrans(const ParameterSet& pset, int warn)
  {
  int ntr = size();			// Number of transitions
  double W, R, Inorm, Iphase;		// Transition values
  int TF = 1;
  int tmp;
  complex ztmp;
  for(int tr=0; tr<ntr; tr++)		// Loop over transitions
    {
    tmp = tr;
    W =      GetFreq(pset, tmp, 1);	// Get the transition frequency
    if(tmp < 0 && warn)			// Failed if tmp < 0
      {
      TF = 0;				// Flag we have failed
      string pname = Gdec(tr);		// For error message
      TTaberror(60, pname, 1);		// Can't read transition
      if(warn > 1) TTabfatality(9);	// Problems during construction
      }
    R =      GetLWR2T2(pset, tr, 0);	// Get the transition rate
    Iphase = GetPhase(pset, tr, 0);	// Get the transition phase
    Inorm =  GetIntensity(pset, tr, 1);	// Get the transition intensity
    if(!Inorm && warn)			// Failed if no intensity
      {
      TF = 0;				// Flag we have failed
      string pname = Gdec(tr);		// For error message
      TTaberror(60, pname, 1);		// Can't read transition
      if(warn > 1) TTabfatality(9);	// Problems during construction
      }
    put(complex(R,W), tr, 0);		// Set {R,W} this transition 
    Iphase *= DEG2RAD;			// Switch to radians
    ztmp = (complexi*Iphase).Zexp();
    put(Inorm*ztmp,tr,1);		// Set intensity this transition 
    }
  return TF;
  } 

// ____________________________________________________________________________
// v                TRANSITIONS TABLE 1D CUTOFF FUNCTIONS
// ____________________________________________________________________________

/* Consider the request for a time domain signal associated with a transitions
   table. If those transitions decay, their time domain signals will be zero
   after a certain amount of time, after which they need not be computed since
   they are effectively zero. The function ExpCutoffs looks at each
   transition and determines the last point to calculate, i.e. the point
   after which the transition will no longer generate any intensity to
   within a cutoff value. These indices are returned in a vector, on index 
   per transition.
 
   If k is the last point in the exponential that should be calculated, and
   cutoff is a value (<0) afterwhich the intensity is essentially zero, then
   we know that
                          k >= -ln(cutoff)/[R*tinc]

   where tinc is the time between points & R the transition decay rate.      */

vector<int> TTable1D::ExpCutoffs(double tinc, int npts, double CO) const
  {
  if(CO>1 || CO<1.e-10) CO = 1.e-4;		// Insure reasonable cutoff
  double mod = -log(CO)/tinc;			// Base in cutoff formula
  double Rtr;					// Transition decay rate
  int i, ntr = rows(); 				// Number of transitions
  vector<int> ifis(ntr);			// Array of final points 
  for(int tr=0; tr<ntr; tr++)			// Loop transiitons
    {
    Rtr = getRe(tr,0);				//   Transition decay rate
    if(Rtr <= 0) ifis[tr] = npts;		//   No decay, all points O.K.
    else					//   If point does decay, find
      {						//   where intensity < CO 
      i = int(mod/Rtr);				//   Point where cutoff is
      if(i>npts-1) ifis[tr] = npts;		//   Insure point within range
      else         ifis[tr] = i+1;		//   or adjust to next point
      }
    }
  return ifis;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A             CLASS TRANSITIONS TABLE 1D CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ---------------------------- Null Constructor ------------------------------

TTable1D::TTable1D() : matrix() { setDefaults(); }

	// Input	none  :
	// Output	TTab1D : A NULL TTable1D (this)


// ----------------------- Construction From Matrix ---------------------------


	// Input	none  :
	// Output	TTab1D : A NULL TTable1D (this)

TTable1D::TTable1D(const matrix& mx, int warn)
  {
  if(!mx.rows()){setDefaults(); return;}// If the matrix is empty
  if(mx.cols() != 2)			// The table may only have two columns
    {
    TTaberror(50, 1);			// Wrong columns in transition table
    TTabfatality(9);			// Problems during construction	
    }
  matrix::operator=(mx);		// Copy the table array
  double fr, lw;			// For frequency & linewidth
  double lwcut = 1.e-12;		// Linewidth cutoff
  string msg;				// For error messages
  for(int i=0; i<mx.rows(); i++)	// Check all relaxation rates >= 0
    {					// (or all transition linewidth)
    lw = getRe(i,0);			//   Linewidth transition i
    if(lw < 0)				//   If negative linewidth, act
      {					//   to correct this problem
      if(fabs(lw) > lwcut) 		//   If above linewidth cutoff
        {
        msg = Gform("%10.3e",lw);	//     For error message
        TTaberror(51, 1);		//     Negative trans. linewidth
        TTaberror(53,msg,1);		//     Linewidth value is msg
        if(warn > 1) TTabfatality(9);	//     Fatal, bad construction	
        TTaberror(52, 1);		//     We set linewith to zero
        }
      fr = getIm(i,0);			//     This is transition frequency
      put(complex(0,fr),i,0);		//     Set linewidth to zero
      }
    }
  setDefaults();			// Set default internal values
  }

TTable1D::TTable1D(const matrix& mx) { *this = TTable1D(mx, 1); }

// --------------------------- Self Construction ------------------------------

TTable1D::TTable1D(const TTable1D& TTab1) : matrix(TTab1) {copySettings(TTab1);}

	// Input	TTab1  : Transitions table
	// None		TTab1D : Transitions table (this), an identical
	//			 copy of TTab1

// ------------------------------ Destruction ---------------------------------

TTable1D::~TTable1D() { }

	// Input	TTab1D: A transitions table (this)
	// Output	none  : TTab is destructed (done in matrix class)

// ------------------------------- Assignment ---------------------------------

	// Input	TTab1  : Transitions table
	// None		TTab1D : Transitions table (this) copied from TTab1

TTable1D& TTable1D::operator = (const TTable1D& TTab1)
{
  matrix::operator=(TTab1);		// Copy the table array
  copySettings(TTab1);			// Copy the internal settings

  return (*this);
}

// ____________________________________________________________________________ 
// B                     TRANSITION TABLE MANIPULATIONS
// ____________________________________________________________________________

// __________________________ Frequency Alterations ___________________________

/*             Function              Return Type & Value
               --------   -----------------------------------------------------
                center     double:  The center frequency of the transitions
                offset     void:    Shift one or all transition frequencies
                FRscale    void:    Scale one or all transition frequencies
                  BC       void:    Baseline correction by removal of
                                    frequencies near zero.

           Input                TTab1D  : Transitions table (this)
                                wa      : Flag for weighted average (def=true)
                                F       : Frequency offset (Hz=default, 1/sec)
                                inHz    : Flag whether F is in Hz or 1/sec
                                tr      : Transition index
                                Fscf    : Frequency scaling factor
                                res     : Frequency resolution
           Output center        double  : The center frequency is returned
           Output offset        void    : All frequency values or that of the
                                          the specified transitions tr are
                                          adjusted as Im(<i|mx|0>) += F
           Output FRscale       void    : All frequency values or that of the
                                          the specified transitions tr are
                                          adjusted as Im(<i|mx|0>) *= Fscf
           Output BC            void    : All transitions havign frequencies
                                          within 0+/-res will be removed     */

double TTable1D::center(bool wa)
  {
  int ntr = rows();			// Number of transitions
  if(!ntr) return 0;			// If no transitions, exit
  double Fcenter;			// Center transition frequency
  if(!wa)				// If not weighted average
    {
    double Fmax = FRmax(); 		// Largest transition frequency
    double Fmin = FRmin(); 		// Smallest transition frequency
    Fcenter = (Fmin-Fmax)/2.0;		// Frequency in spectrum center
    }
  else
    {
    double Fri;				// For transition i's frequency
    double Ii;				// For transition i's intensity
    double Isum = 0.0;			// For summed intensities
    double WFsum = 0.0;			// For summed weighted frequencies
    for(int i=0; i<ntr; i++)		// Loop over all transitons
      {
      Ii  = norm(I(i));			//   Transition i intensity
      Fri = getIm(i,0);			//   Transition i frequency
      WFsum += Ii*Fri; 			//   Weighted frequency sum
      Isum  += Ii;			//   Intensity sum
      }
    Fcenter = WFsum/Isum;		// Weighed freqeuncy center
    }
  offset(Fcenter, 0);			// Offset frequencies, center now ~0
  return Fcenter;			// Return this center
  }

void TTable1D::offset(double F, int inHz)
  {
  double W = F;				// Offset frequency
  if(inHz) W *= PIx2;			// Convert to 1/sec if needed
  complex zadd(0,W);			// Value to add to frequencies
  for(int i=0; i<rows(); i++)		// Add value to all transitions
    put(get(i,0)+zadd,i,0);
  }

void TTable1D::offset(double F, int tr, int inHz)
  {
  if(tr<0 || tr>=rows()) return;	// Do nothing if no transition tr
  double W = F;				// Offset frequency
  complex z(0,W);			// Value to add to frquency
  if(inHz) W *= PIx2;			// Convert to 1/sec if needed
  put(get(tr,0)+z,tr,0);		// Add value to transition
  }

void TTable1D::FRscale(double Fscf)
  {
  complex zF;					// Adjusted frquency
  for(int i=0; i<rows(); i++)			// Loop over transitions
    {   
    zF = complex(getRe(i,0), getIm(i,0)*Fscf);	//   Adjusted frequency
    put(zF, i, 0);				//   Set new frequency
    }   
  } 

void TTable1D::FRscale(double Fscf, int tr)
  {
  if(tr<0 || tr>=rows()) return;	// Do nothing if no transition tr
  complex zF = complex(getRe(tr,0), getIm(tr,0)*Fscf);
  put(zF, tr, 0);
  } 

void TTable1D::BC(double res)
  {
  if(res < 0) res = 0;				// Minimum range
  for(int i=0; i<rows(); i++)			// Loop over transitions
    {
    if(fabs(getIm(i,0)) <= res)			// If frequency 0 within res
      put(complex0, i, 1);			// zero the transition
    }
  }

// __________________________ Intensity Alterations ___________________________

	// Input		TTab1D	: Transitions table (this)
        //			Iscf	: Intensity scaling factor
	//			tr      : Transition index
	// Output		void	: Transition intensities are
	//				  adjusted as <i|mx|1> *= Iscf
	//				  where i spans all transitions or just
	//				  a single one specified by tr
 
void TTable1D::Iscale(double Iscf)
  { for(int i=0; i<rows(); i++) put(Iscf*get(i,1), i, 1); } 
 
void TTable1D::Iscale(double Iscf, int tr)
  {
  if(tr<0 || tr>=rows()) return;	// Do nothing if no transition tr
  put(Iscf*get(tr,1), tr, 1);		// Adjust intensity of transition tr
  } 
 
void TTable1D::Iscale(const complex& Iscf)
  { for(int i=0; i<rows(); i++) put(Iscf*get(i,1), i, 1); } 
 
void TTable1D::Iscale(const complex& Iscf, int tr)
  {
  if(tr<0 || tr>=rows()) return;	// Do nothing if no transition tr
  put(Iscf*get(tr,1), tr, 1);		// Adjust intensity of transition tr
  } 

void TTable1D::Iremove(double dcut)
  {
  int nr = rows();			// Input transitions
  if(!nr) return;			// Quit if no transitions
  if(dcut <= 0.0) dcut = ICUT;		// Use ICUT if dcut not specified
  vector<int> trs;			// Transition indices
  int i,j;				// Transitions index
  double Ii;				// Transition intensity (norm)
  for(i=0; i<nr; i++)			// Loop the transitions
    {
    Ii = norm(get(i,1));		// Transition intensity
    if(Ii >= dcut) trs.push_back(i);	// Store if intensity OK
    }
  nr = trs.size();			// New transition count
  if(!nr) { *this=TTable1D(); return; }	// If not transitions, table empty
  matrix mx(nr, 2);			// New array for transitions
  for(i=0; i<nr; i++)			// Loop rows of array to add
    {					// transitions with intensity > dcut
    j = trs[i];
    mx.put(get(j,0), i,0);
    mx.put(get(j,1), i,1);
    }
  *this = TTable1D(mx);			// We are now a new table
  }

// __________________________ Linewidth Alterations ___________________________

	// Input		TTab1D	: Transitions table (this)
	//			LWR	: Linewidth/Rate adjustment (Hz)
	//			inHz	: Flag whether LWR is a linewidth (Hz)
	//				  or a rate (1/sec) <-- default
	// Output		void	: The linewidth values are all adjusted
	//				  as Re(<i|mx|0>) += F
	// Note				: Negative linewidths are not allowed,


void TTable1D::broaden(double LWR, int inHz)
  {
  double R = LWR; 				// Rate in 1/sec
  if(inHz>0)					// If this >0 (default) then 
    {						// take input value as lwhh
    R *= PI;					// in Hz. Convert to rad/sec
    }						// by relationship R = PI*lwhh
  else if(inHz<0)				// If this is <0 then take the
    {						// input value as lwpp (Gauss)
    R *= sqrt(0.75);				// Convert to rad/sec using
    R *= HZ2RAD/HZ2GAUSS;			// R = sqrt(3/4)*lwpp
    }
  for(int i=0; i<rows(); i++)			// We will now add in this 
    {						// amount of linewidth, making 
    (*this).put((*this)(i,0) + R, i,0);		// sure we never go below 0
    if(getRe(i,0)<0) (*this).put(0,i,0);
    }
  }


void TTable1D::broaden(double LWR, int tr, int inHz)
  {
  if(tr<0 || tr>=rows()) return;		// Nothing if no transition tr
  double R = LWR; 				// Rate in 1/sec
  if(inHz) R *= PI;				// Adjust if lwhh(Hz) input
    (*this).put((*this)(tr,0) + R, tr, 0);
  if(getRe(tr,0)<0) (*this).put(0,tr,0);
  }


// __________________________ Resolution Alterations __________________________


	// Input		TTab1D	: Transitions table (this)
        //			res     : A resolution in Hz
	// Output		void	: The transitions matrix is
	//			  	  modified by the blending of
	//			  	  transitions that are within res
	//			  	  from each another in frequency

void TTable1D::resolution(double res)
  {
  int nr = rows();			// Input transitions
  if(!nr) return;			// Quit if no transitions
  int i, I, tr=1, resolved=0;		// Count transitions, resolved
  double vi, vI;
  res *= 2.0*PI;			// Convert res to rad/sec
  for(i=1; i<nr; i++)			// Loop the transitions
    {
    vi = getIm(i,0);			// Trans. frequency (rad/sec)
    resolved = 1;                       // Flag its resolved
    for(I=0; I<tr && resolved; I++)     // Loop resolved transitions
      {
      vI = getIm(I,0);			// Trans. I freqeuncy (rad/sec)
      if(fabs(vi-vI) <= res)            // If i isn't resolved from I
        {
        put(get(I,1)+get(i,1),I,1);
        resolved = 0;
        }
      }
    if(resolved)                        // If resolved transition then
      {                                 // put it at the next spot in
      put(get(i,0),tr,0);		// the transitions array
      put(get(i,1),tr,1);
      tr++;
      } 
    }   
  *this = get_block(0,0,tr,2);
  }

// ____________________________________________________________________________
// C                    Time Domain Spectra Generation
// ____________________________________________________________________________

// sosi - The cutoff points, ifi array below, seem to be bad sometimes
//        by being large negative numbers.... why so? My guess is that the
//        linewidths are goofy (i.e. negative too).


	// Input	TTab1D: A transitions table (this)
	// 		sigma : Density matrix (operator propagated)

	//		N     : Number of point desired
        // or
	//		data  : 1D row vector for time interferrogram

	//		tinc  : Time increment between points (sec)
	//		ICUT  : Intensity cutoff
	//		cutper: % below max. to ignore [0,1]
	// Output	data  : 1D data vector is filled with a 
	//			time domain spectrum.
	// Note		      : Value cutper is a % in decimal form
	//			of the maximum exponential intensity
	//			considered worth adding to the spectrum.
	// Note		      : Each exponential has intensity above
	//			cutoff between the points [0-ifi]. Only
	//			those points used in generating the spectrum
	// Note		      : Users must insure data is first zeroed!
	// Note		      : The input time increment is adjusted to
	//			account for "Hz" output.

row_vector TTable1D::T(int npts, double tinc) const
  {
  row_vector data(npts, complex0);		// Array for the spectrum
  T(data, tinc);				// Use function overload
  return data;
  }

void TTable1D::T(row_vector& data, double tinc) const
  {
  int ntr = rows();				// Get number of transitions
  if(!ntr) return;				// Exit if no transitions
  int npts = data.size();			// Get the number of points
  if(!npts) return;				// Exit if no points
  complex* newvec = NULL;			// New time domain spectrum
  newvec = new complex[npts];			// Make this the same size
  vector<int> ifi=ExpCutoffs(tinc,npts,PERCUT);	// Set point cutoffs
  complex Itr, z, expatit;			// Transition tr intensity
  int i, lst, tr;
  double mtinc = -tinc;				// Need negative in exponent
  
  complex *nvo = newvec;			// Start of new vector
  complex *nvf = nvo + npts;			// End of new vector
  complex *nvi;					// Counter in vector
  for(nvi=nvo; nvi<nvf; nvi++)			// Zero all points in vector
    *nvi = complex0;

  for(tr=0; tr<ntr; tr++)			// Loop over all transitions
    {
    Itr = get(tr,1);				// Intensity of transition tr
    if(norm(Itr) > ICUT)			// Use trans. if intensity O.K.
      {
      z = FRQREV?get(tr,0).conj():get(tr,0);
      expatit = (z*mtinc).Zexp();		// Evolution for time inc.
      lst = ifi[tr];				// Last point tr contributes to
      if(lst < 0) lst = npts;			// If negative, use all points
      nvf = nvo + lst;				// Last point in new vector
      z = complex1;				// Start at time = 0
      for(nvi=nvo; nvi<nvf; nvi++)		// Loop over contributing pts
        {					// from this transition
        *nvi += Itr*z;				//   Add intensity from trans.
        z *= expatit;				//   Evolve trans. to next time
        }					// Next spectrum point, this tr
      }						// Finished spec from this point
    }						// Next transition contribution

  for(i=0,nvi=nvo; i<npts; i++, nvi++)		// Add new vector (time domain
    data.put(data.get(i)+(*nvi),i);		// spectrum) to input vector
  if(SN) Noise(data, Noisemax());		// Add random noise if desired
  delete [] newvec;				// This was only temporary
  }

vector<row_vector> TTable1D::Ts(int npts, double tinc) const
  {
  vector<row_vector> vects;			// Vector of vectors
  int ntr = rows();				// Get number of transitions
  if(!ntr) return vects;			// Exit if no transitions
  row_vector data, data0(npts, complex0);	// Array for single spectra
  vector<int> ifi=ExpCutoffs(tinc,npts,PERCUT);	// Set point cutoffs
  complex Itr, z, expatit;			// Transition tr intensity
  int i, tr;					// Point, transition indices
  if(!FRQREV)
    {
    for(tr=0; tr<ntr; tr++)			// Loop over all transitions
      {
      Itr = get(tr,1);				// Intensity of transition tr
      if(norm(Itr) > ICUT)			// Use trans. if intensity O.K.
        {
        data = data0;				// Begin with zero spectrum
        expatit = (-get(tr,0)*tinc).Zexp();	// Evolution for time inc.
        z = complex1;				// Start at time = 0
        for(i=0; i<ifi[tr]; i++)		// Loop over contributing pts
          {					// from this transition
          data.put(Itr*z,i);			//   Set intensity from trans.
          z *= expatit;				//   Evolve trans. to next time
          }
        if(SN) Noise(data, Noisemax());		// Add random noise if desired
        vects.push_back(data);			// Store this spectrum
        }
      }
    }
  else
    {
    for(tr=0; tr<ntr; tr++)			// Loop over all transitions
      {
      Itr = get(tr,1);				// Intensity of transition tr
      if(norm(Itr) > ICUT)			// Use trans. if intensity O.K.
        {
        data = data0;				// Begin with zero spectrum
        z = -get(tr,0);				// We'll reverse the frequencies
        expatit = (z.conj()*tinc).Zexp();	// Evolution for time inc.
        z = complex1;				// Start at time = 0
        for(i=0; i<ifi[tr]; i++)		// Loop over contributing pts
          {					// from this transition
          data.put(Itr*z,i);			//   Set intensity from trans.
          z *= expatit;				//   Evolve trans. to next time
          }
        if(SN) Noise(data, Noisemax());		// Add random noise if desired
        vects.push_back(data);			// Store this spectrum
        }
      }
    }
  return vects;
  }

vector<int> TTable1D::TCutoffs(int npts, double tinc) const
  { return ExpCutoffs(tinc,npts,PERCUT); }

// ____________________________________________________________________________ 
// D                        Frequency Domain Spectra
// ____________________________________________________________________________

	// Input	TTab1D	: A transitions table (this)
	//		N     	: Number of points desired
	//		data	: 1D row vector for frequency spectrum
	//		Fst	: Starting frequency (Hz)
	//		Ffi	: Final frequency (Hz)
	//		ICUT	: Intensity cutoff
	//		PERCUT  : % below max. to ignore [0,1]
	// Output	data  	: Row vector filled with a frequency spectrum.
	// Note		      	: Value PerCUT is a % in decimal form
	//			  of the maximum Lorentzian intensity
	//			  considered worth adding to the spectrum.
	// Note		      	: Each Lorentzian has intensity above
	//			  ICUT*Imax between the points bound
	//			  by the indices ist and ifi.  Only those
	//			  points are used in generating the spectrum
	// Note		      	: To avoid poor resolution, if the frequency
	//			  between points, winc, does not fullfill the
	//			  relationship 5*winc < lwhh, then the
	//			  integrated Lorentzian intensities are used
	//		          rather than the point value.
	// Note			: Users must insure data is first zeroed!


row_vector TTable1D::F(int npts, double Fst, double Ffi) const
  {
  row_vector data(npts, complex0);		// Array for the spectrum
  F(data, Fst, Ffi);				// Use function overload
  return data;
  }

void TTable1D::F(row_vector& data, double Fst, double Ffi) const
  {
  int ntr = rows();				// Get number of transitions
  int npts = data.size();			// Get number of points
  double wo =   PIx2*Fst;			// Init. frequency (1/seconds)
  double winc = PIx2*(Ffi-Fst)/(npts-1);	// Frequency between pts(1/sec)
  int *ist, *ifi;				// Init, final indices, each tr
  ist = new int[ntr];
  ifi = new int[ntr];
  Lorentz_cut(ist,ifi,*this,wo,winc,npts,PERCUT);// Fill index arrays ist & ifi
  int *Lint;					// Flags if integration needed
  Lint = new int[ntr];
  Lorentzian(Lint,*this,winc);			// Fill Lint integration flags

  double w, dw, dw2;				// Frequency, difference freq.
  double Rtr, Rtr2, Wtr;			// Transition tr rate & freq.
  complex Itr;					// Transition tr intensity
  double dwa,dwb,re,im;
  int i;					// Point indices
  complex z;
  for(int tr=0; tr<ntr; tr++)			// Loop over all transitions
    {
    Itr = get(tr,1);				// Intensity of transition tr
    if(norm(Itr) > ICUT)			// Use trans. if intensity O.K.
      {
      Rtr = getRe(tr,0);			// Rate for transition tr
      if(Rtr < 1.e-4) Rtr = 1.e-4;		// (Minimimum line width)
      Rtr2 = Rtr*Rtr;
      Wtr = getIm(tr,0);			// Frequency of transition tr
      w = wo + winc*ist[tr];			// Set w 1st contributing pt
      dw = w - Wtr;				// Set dw of first contributor
      dw2 = dw*dw;
      if(!Lint[tr])				// This section if NO integral
        for(i=ist[tr];i<ifi[tr]; i++,dw+=winc)	// Loop over contributing pts
          {					// No integration (averaging)
          z = data.get(i);			// required in this case
          dw2 = dw*dw;				// Square the frequency
          z += Itr*complex(Rtr,dw)/(Rtr2+dw2);	// Use Lorentzian value at dw
          data.put(z,i);			// Set the data point
          }
      else					//This section if integrating
        for(i=ist[tr];i<ifi[tr]; i++,dw+=winc)	// Loop over contributing pts
          {					// using ave. integrated value
          z = data.get(i);			//	Current value
          dwa = dw - winc/2.0;			//	Lower frequency
          dwb = dw + winc/2.0;			//	Upper frequency
// sosi - need to check for when Rtr == 0!
if(!Rtr) Rtr = 1.e-4;
          re = atan(dwb/Rtr)-atan(dwa/Rtr);	//	Real integral
          im = -0.5*log((Rtr2+dwb*dwb)/(Rtr2+dwa*dwa));//
          z += Itr*complex(re/winc,im/winc);	//	Add in integral value
          data.put(z,i);			//	Store new value
          }
      }
    }
  if(SN) Noise(data, Noisemax());		// Add random noise if desired
  delete [] ifi;
  delete [] ist;
  delete [] Lint;
  return;
  }

vector<row_vector> TTable1D::Fs(int npts, double Fst, double Ffi) const
  {
  vector<row_vector> vects;			// Vector of vectors
  row_vector data, data0(npts, complex0);	// Array for single spectra
  int ntr = rows();				// Get number of transitions
  double wo =   PIx2*Fst;			// Init. frequency (1/seconds)
  double winc = PIx2*(Ffi-Fst)/(npts-1);	// Frequency between pts(1/sec)
  int *ist, *ifi;				// Init, final indices, each tr
  ist = new int[ntr];
  ifi = new int[ntr];
  Lorentz_cut(ist,ifi,*this,wo,winc,npts,PERCUT);// Fill index arrays ist & ifi
  int *Lint;					// Flags if integration needed
  Lint = new int[ntr];
  Lorentzian(Lint,*this,winc);			// Fill Lint integration flags

  double w, dw, dw2;				// Frequency, difference freq.
  double Rtr, Rtr2, Wtr;			// Transition tr rate & freq.
  complex Itr;					// Transition tr intensity
  double dwa,dwb,re,im;
  int i;					// Point indices
  complex z;
  for(int tr=0; tr<ntr; tr++)			// Loop over all transitions
    {
    Itr = get(tr,1);				// Intensity of transition tr
    if(norm(Itr) > ICUT)			// Use trans. if intensity O.K.
      {
      data = data0;				// Zero spectrum this transition
      Rtr = getRe(tr,0);			// Rate for transition tr
if(Rtr < 1.e-4) Rtr = 1.e-4;
      Rtr2 = Rtr*Rtr;
      Wtr = getIm(tr,0);			// Frequency of transition tr
      w = wo + winc*ist[tr];			// Set w 1st contributing pt
      dw = w - Wtr;				// Set dw of first contributor
      dw2 = dw*dw;
      if(!Lint[tr])				// This section if NO integral
        for(i=ist[tr];i<ifi[tr]; i++,dw+=winc)	// Loop over contributing pts
          {					// No integration (averaging)
          dw2 = dw*dw;				// Square the frequency
          z = Itr*complex(Rtr,dw)/(Rtr2+dw2);	// Use Lorentzian value at dw
          data.put(z,i);			// Set the data point
          }
      else					//This section if integrating
        for(i=ist[tr];i<ifi[tr]; i++,dw+=winc)	// Loop over contributing pts
          {					// using ave. integrated value
          dwa = dw - winc/2.0;			//	Lower frequency
          dwb = dw + winc/2.0;			//	Upper frequency
// sosi - need to check for when Rtr == 0!
if(!Rtr) Rtr = 1.e-4;
          re = atan(dwb/Rtr)-atan(dwa/Rtr);	//	Real integral
          im = -0.5*log((Rtr2+dwb*dwb)/(Rtr2+dwa*dwa));//
          z = Itr*complex(re/winc,im/winc);	//	Add in integral value
          data.put(z,i);			//	Store new value
          }
      if(SN) Noise(data, Noisemax());		// Add random noise if desired
      vects.push_back(data);			// Store this spectrum
      }
    }
  delete [] ifi;
  delete [] ist;
  delete [] Lint;
  return vects;
  }

// ____________________________________________________________________________ 
// E             Frequency Domain Spectra In Derivative Mode
// ____________________________________________________________________________

	// Input	TTab1D: A transitions table (this)
	// 		sigma : Density matrix (operator propagated)

	//		N     : Number of points desired
        // or
	//		data  : 1D row vector for frequency spectrum
	// or
	// 		mx    : Matrix of transition data

	//		fstart: Starting frequency (Hz)
	//		fend  : Final frequency (Hz)
	//		cutoff: Intensity cutoff
	//		cutper: % below max. to ignore [0,1]
	// Output	data  : Row vector filled with a
	//			frequency spectrum.
	// Note		      : Value cutper is a % in decimal form
	//			of the maximum Lorentzian intensity
	//			considered worth adding to the spectrum.
	// Note		      : Each Lorentzian has intensity above
	//			cutoff*Imax between the points bound
	//			by the indices ist and ifi.  Only those
	//			points are used in generating the spectrum
	// Note		      : To avoid poor resolution, if the frequency
	//			between points, winc, does not fullfill the
	//			relationship 5*winc < lwhh, then the
	//			integrated Lorentzian intensities are used
	//		        rather than the point value.
	// Note		      : Users must insure data is first zeroed!


row_vector TTable1D::FD(int N, double fstart, double fend) const
  {
  row_vector data(N, complex0);			// Array for the spectrum
  FD(data, fstart, fend);			// Use function overload
  return data;
  }

void TTable1D::FD(row_vector& data, double fstart, double fend) const

  {
  int ntr     = rows();				// Get number of transitions
  int npts    = data.size();			// Get number of points
  double wo   = fstart*HZ2RAD;			// Init. frequency (1/seconds)
  double winc = HZ2RAD*(fend-fstart)/(npts-1);	// Frequency between pts(1/sec)
// sosi - Still crude for esr, PERCUT here someday too
  double Rtr, Wtr;				// Transition tr rate & freq.
  complex Itr;					// Transition tr intensity
  for(int tr=0; tr<ntr; tr++)			// Loop over all transitions
    {
    Itr = get(tr,1);				// Intensity of transition tr
    if(norm(Itr) > ICUT)			// Use trans. if intensity O.K.
      {
      Rtr = getRe(tr,0);			// Rate for transition tr
      if(Rtr < 1.e-4) Rtr = 1.e-4;
      Wtr = getIm(tr,0);			// Frequency of transition tr
Rtr /= winc;					// Convert rate to points
Wtr /= winc;					// Convert frequency to points
Wtr -= wo/winc;					// Offset to proper point
      data += Itr*DLorentzian(npts, int(Wtr), Rtr);
      }
    }
  if(SN) Noise(data, Noisemax());		// Add random noise if desired
  return;
  }


// ____________________________________________________________________________
// F                            PHASE CORRECTION
// ____________________________________________________________________________

/* These routines apply phase corrections to the transition table.           */


        // Input        TTab1D: A transitions table (this)
	//		w0    : Zero point pivot frequency (Hz)
	//		w1    : 1st order pivot frequency (Hz)
	//		order : Maximum allowed correction order 
	// Output       z     : A complex number whose real part
	//			is the zeroth order phase correction
	//			and whose imaginary is the 1st order
	//			phase correction
	// Note		      : Transition table phase values adjusted 
	//		w0    : Set to 1st order pivot (rad/sec)

/* This routine phase corrects all the transitions in matrix mx using a
   0th order and 1st order phase correction.  The 0th order finds the
   transition closest to w0, takes the phase of that transition, and then
   removes that phase component from all transitions.  The 1st order finds
   the transition closest to w1 and looks at its phase.  It then assumes that
   there is some constant tp such that

                             (w1-w0)*tp = phase(w1) 

   and calculates tp, noting that and addition of 360 degrees to the phase will
   suffice, i.e. that there is another factor, ntp, where

                               (w1-w0)*ntp = 360

   such that

                      (w1-w0)*(tp+n*ntp) = phase(w1) + n*360 

   where n is an integer.  Any value of n will produce a proper phase value
   of (w1-w0)*(tp+n*ntp) that would undo the phase at w1.  The only requirement
   is to select which value of n to use.  To to that, this routine picks the
   biggest transition, say Wtr, and sees which value of n minimizes

                      (Wtr-w0)*(tp+n*ntp) - phase(Wtr) -----> 0

   where n is kept between the range [0, order].  Once found, all transitions
   are phase corrected by adjusting their complex intensities according to

  		    Inew(W) = Iold(W) * exp[-i*(W-w0)*(tp+n*ntp)]            */

complex TTable1D::pcorrect(double& w0, double w1, int order)
  {
//		     First Nearest Transitions To w0 and w1

  double pi2 = 2.0*PI;				// Required constant
  double w0rad = w0*pi2;			// Zero pivot in radians
  double w1rad = w1*pi2;			// 1st pivot in radians
  double Wtr;					// Transition frequency
  double Wdel0, Wdelmax0 = HUGE_VAL;		// Use to look for w0
  double Wdel1, Wdelmax1 = HUGE_VAL;		// Use to look for w1
  int i, i0max=0, i1max=0;
  for(i=0; i<rows(); i++)			// Loop rows
    {
    Wtr = getIm(i,0);				// 	Transition i freq.
    Wdel0 = w0rad - Wtr;			//	Freq. Difference
    if(fabs(Wdel0) < Wdelmax0)			//	Look for a 0th order
       {					//	pivot frequency close
       Wdelmax0 = fabs(Wdel0);			//	by w0rad	
       i0max = i;
       }
    Wdel1 = w1rad - Wtr; 			//	Freq. Difference
    if(fabs(Wdel1) < Wdelmax1)			//	Look for 1st order
       {					//	pivot frequency close
       Wdelmax1 = fabs(Wdel1);			//	by w1rad
       i1max = i;
       }
    }
  w0rad = getIm(i0max, 0);			// Set zero pivot transition
  if(i0max == i1max)				// Insure both pivots aren't
    {						// the same point
    i1max = i0max-1;
    if(i1max < 0) i1max = i0max+1;
    }
  w1rad = getIm(i1max, 0);			// Set 1st pivot transition

//		   Now Apply The Zeroth Order Phase Correction
//	         (Transition w0 Will Have Zero Phase Hereafter)

  complex z = get(i0max,1);			// Zero order intensity
  double phase0 = phase(z);			// Zero order phase
  complex pc = norm(z)/z;			// Zero order phase correction
  for(i=0; i<rows(); i++)			// Loop rows
    {
    z = get(i,1);				// 	Transition intensity
    z *= pc;					// 	Adjust transition phase
    put(z,i,1);					// 	Reset the intensity
    }

//		    Set Up For First Order Phase Correction

  z = get(i1max,1);				// 1st order intensity
  double phase1 = phase(z);			// 1st order phase
  if(phase1 < 0) phase1 += pi2;			// Keep phase positive
  double tp = phase1/(w1rad-w0rad);		// Base time constant for phase
  double ntp = pi2/(w1rad-w0rad);		// Adjustable time constant

//        Find Largest Transition Which is Not One Of The Two Pivots
//	      (We'll Deal With Case Of <= Two Transitions Later)

  int maxi=0;					// Index of large transition
  double Imax = -HUGE_VAL;				// Largest transition intensity
  for(i=0; i<rows(); i++)			// Loop transitions
    {
    if(i!=i0max && i!=i1max)			// 	Don't use a pivot
      {
      if(norm(get(i,1)) > Imax)			//	Check intensity
        {					//	and store if large
        Imax = norm(get(i,1));
        maxi = i;
        }
      }
    }
  double wrad = getIm(maxi,0);			// Frequency of large transition
  double prad = phase(get(maxi,1));		// Phase of large transition
  if(prad < 0) prad += pi2;			// Keep phase positive

//                 Find Best Value of "n" For 1st Order Phasing

  double delp, pp, pmin=HUGE_VAL;
  int nmin = 0;
  for(i=0; i<order; i++)			// Loop through allowed orders
    {
    pp = (wrad-w0rad)*(tp+i*ntp);		//	Phase made at wrad, n=i
    delp = fabs(prad-pp); 			//	See if good adjustment
    if(delp < pmin)				//	If its good, store this 
      {						//	value of i & use it if
      pmin = delp;				//	nothing better is found
      nmin = i;
      }
    }

//		     Apply The First Order Phase Correction

  double tconst = tp+(double(nmin)*ntp);
  double delW;
  for(i=0; i<rows(); i++)			// Loop rows
    {
    Wtr = getIm(i,0);				// 	Transition i freq.
    delW = Wtr - w0rad;				//	Freq. Difference
    pp = delW*(tp+nmin*ntp);			// 	Needed phase correction
    z = complex(0, pp);				//	-i*phi
    z = z.Zexp();				//	z = exp(-i*phi)
    pc = norm(z)/z;				//	Phase correction factor
    z = get(i,1);				// 	Transition intensity
    z *= pc;					// 	Adjust transition phase
    put(z,i,1);					// 	Reset (phase) intensity
    }
  w0 = w0rad;					// Set w0 to be pivot point
  return complex(phase0, tconst);		// Return phase factors
  }



        // Input        TTab1D: A transitions table (this)
	//		wpivot: 1st order pivot frequency (Hz)
	//		P     : Phase correction factors
	//			  Re(P) = 0th order p.c. in radians
	//			  Im(P) = 1st order time const. (sec)
	// Output       void  : The transition phases in table are
	//			all phased according to the correction
	//			specified in Wpivot and P

/* This routine phase corrects all the transitions in matrix mx using a
   0th order and 1st order phase correction.  Using phi0 as the zeroth
   order phase [phi0 = Re(P)] and t1 as the 1st order correction time
   constant [t1 = Im(P)], the correction to a transition at frequency W
   is given by

  				    -i*phi0    -(W-Wpivot)*t1
  		I   (W) = I  (W) * e        * e
  		 out       in                                                */

void TTable1D::pcorrect(double Wpivot, complex& P)
  {
//	   First Apply The Zeroth Order Phase Correction

  double phase0 = zRe(P);			// Zero order phase (radians)
  complex pc0 = (-complexi*phase0).Zexp();	// Zero order phase correction
  int i=0;
  complex z;
  for(i=0; i<rows(); i++)			// Loop rows
    {
    z = get(i,1);				// 	Transition intensity
    z *= pc0;					// 	Adjust transition phase
    put(z,i,1);				// 	Reset the intensity
    }

//		    Now Apply The First Order Phase Correction

  double t1 = zIm(P);				// 1st order time const. (sec)
  double delW = 0;				// Diff. frequency (rad/sec)
  complex pc1;					// 1st order phase correction
  for(i=0; i<rows(); i++)                    // Loop rows
    {
    delW = getIm(i,0) - Wpivot;		//      Freq. Difference
    pc1 = (-complexi*delW*t1).Zexp();		//      Phase correction
    z = get(i,1);                            //      Transition intensity
    z *= pc1;					//      Adjust transition phase
    put(z,i,1);                              //      Reset (phase) intensity
    }
  }


// ____________________________________________________________________________ 
// G            CLASS TRANSITIONS TABLE 1D AUXILIARY FUNCTIONS
// ____________________________________________________________________________ 

/* These functions allow users to access individual elements of the transitions
   table.  They each take a transition index and will return 0 if the index
   provided is out of range.     

               Function              Return Type & Value
               --------   ------------------------------------------
                  R2       double:  The transition relaxation rate
                  Fr       double:  The transition frequency
                  I        complex: The transition intensity
                  Tr       vector:  The transition itself                    */

double TTable1D::R2(int tr) const
  { if(tr<0 || tr>= rows()) return 0; return getRe(tr,0); }

double TTable1D::Fr(int tr) const
  { if(tr<0 || tr>= rows()) return 0; return getIm(tr,0); }

complex TTable1D::I(int tr) const
  { if(tr<0 || tr>= rows()) return complex0; return get(tr,1); }

row_vector TTable1D::Tr(int tr) const
  { 
  if(tr<0 || tr>= rows()) 
    return row_vector(2,complex0);
  matrix tmx = mx();
  return row_vector(tmx.get_block(tr,0,1,2));
  }

/* These functions allow users to access values that are global over all
   transitons in the table.  No transition index is required.

        Function                   Return Type & Value
       -----------   ------------------------------------------------------
       LineWidths    int: TRUE if any transitions have non-zero linewidths
       Intensities   int: TRUE if any transitions have intensity >= ICUT
       Phases	     int: TRUE if any transtions have phases > .0001 rad
       size          int: Total number of transitions in the table
       FRmax         double: Largest transition frequency in the table
       FRmin         double: Smallest transition frequency in the table
       Tdmin         double: Smallest dwell time needed for time domain spec	
       LWmax         double: Largest transition linewidth in the table
       Imax          double: Largest transition intensity in the table
       Noisemax      double: Largest noise signal present in table spectrum
       mx            matrix: Table as a matrix                               */


bool TTable1D::LineWidths() const
  {
  bool TF = false;				// Assume none have width
  for(int i=0; i<rows() && !TF; i++)		// Determine if any
    if(getRe(i,0) > 1.e-6) TF=true; 		// lines have finite width
  return TF;
  }

bool TTable1D::Intensities() const
  {
  bool TF = false;				// Assume none have width
  for(int tr=0; tr<rows() && !TF; tr++)		// Determine if any
    if(norm(get(tr,1)) > ICUT) TF=true; 	// lines have finite intensity
  return TF;
  }

bool TTable1D::Phases() const
  {
  bool TF = false;				// Assume none have phase
  for(int i=0; i<rows() && !TF; i++)		// Determine if any
    if(phase(getRe(i,1)) > 1.e-4) TF=true; 	// lines have finite phase
  return TF;
  }

int TTable1D::size() const { return rows(); }

double TTable1D::FRmax() const
  {
  int ntr = rows();			// Number of transitions
  if(!ntr) return 0;			// Exit if no transitions
  double W=0;				// For transition frequency
  double I=0;				// For transition intensity
  int    i=0;				// For transition index
  double Wmax = 0;			// For max. transition freq.
  while(norm(get(i,1)<ICUT) && i<ntr)  // Look for 1st non-zero
    i++;
  Wmax = i<ntr?getIm(i,0):0;		// Set max to 1st nonzero freq.
  for(; i<ntr; i++)			// Loop rest of transitions
    { 
    I = norm(get(i,1));			//   Transition i intensity
    W = getIm(i,0);			//   Transition i frequency
    if(I> ICUT && W>Wmax) Wmax = W;	//   Adjust maximum frequency
    }
  return Wmax/PIx2;			// Return maximum frequency (Hz)
  }

double TTable1D::FRmin() const
  {
  int ntr = rows();			// Number of transitions
  if(!ntr) return 0;			// Exit if no transitions
  double W=0;				// For transition frequency
  double I=0;				// For transition intensity
  int    i=0;				// For transition index
  double Wmin = 0;			// For min. transition freq.
  while(norm(get(i,1)<ICUT) && i<ntr)  // Look for 1st non-zero
    i++;
  Wmin = i<ntr?getIm(i,0):0;		// Set min to 1st nonzero freq.
  for(; i<ntr; i++)			// Loop rest of transitions
    { 
    I = norm(get(i,1));			//   Transition i intensity
    W = getIm(i,0);			//   Transition i frequency
    if(I> ICUT && W<Wmin) Wmin = W;	//   Adjust minimum frequency
    }
  return Wmin/PIx2;			// Return minimum frequency (Hz)
  }

double TTable1D::Tdmin() const
  {
  int ntr = rows();			// Number of transitions
  if(!ntr) return 0;			// Exit if no transitions
  double W=0, Wmax = fabs(getIm(0,0));	// Set max to 1st transition freq.
  for(int i=1; i<ntr; i++)		// Loop rest of transitions
    { 
    W = fabs(getIm(i,0));		// Transition freq
    if(W > Wmax) Wmax = W;		// Maximum frequency
    }
  return 1/(2.0*Wmax);			// Return minimum dwell time (sec)
  }

double TTable1D::LWmax() const
  {
  int ntr = rows();			// Number of transitions
  if(!ntr) return 0;			// Exit if no transitions
  double LW, LWmax=getRe(0,0);
  for(int i=1; i<rows(); i++)		// Loop transitions
    {
    LW = getRe(i,0);			// Transition linewidth
    if(LW > LWmax) LWmax = LW;		// Maximum linewidth
    }
  return LWmax/PI;			// Return max. linewidth (Hz)
  }

double TTable1D::LWmin() const
  {
  int ntr = rows();			// Number of transitions
  if(!ntr) return 0;			// Exit if no transitions
  double LW, LWmin=getRe(0,0);
  for(int i=1; i<rows(); i++)		// Loop transitions
    {
    LW = getRe(i,0);			// Transition linewidth
    if(LW < LWmin) LWmin = LW;		// Maximum linewidth
    }
  return LWmin/PI;			// Return min. linewidth (Hz)
  }

double TTable1D::Imax() const
  {
  int ntr = rows();			// Number of transitions
  if(!ntr) return 0;			// Exit if no transitions
  double I, Imax = norm(get(0,1));	// Set I max to 1st transition
  for(int k=1; k<rows(); k++)		// Loop transitions
    {
    I = norm(get(k,1));			// Transition intensity
    if(I > Imax) Imax = I;		// Maximum linewidth
    }
  return Imax;				// Return maximum intensity
  }

double TTable1D::Noisemax() const
  {
  if(!SN) return 0;			// No maximum noise if no S/N value set
  double maxI = Imax();			// Get the maximum intensity
  return maxI/SN;			// Maximum noise possible
  }

matrix TTable1D::mx() const { return matrix(*this); }

	// Input	TT1	: A transitions table
	// 		TT2	: Another transitions table
        //              res	: Resolution (in radians/sec)
        // Output	TT	: Transitions table which is a blend of
	//			  the transitions from TT1 & TT2.  Transitions
	//			  within res frequency from one another are
	//			  blended into a single transition with a 
	//			  summed intensity and weighted frequency
 
TTable1D sum(const TTable1D& TT1, const TTable1D& TT2, double res)
  {
  int nr1 = TT1.rows();                         // Rows of original array
  if(!nr1) return TT2;                          // 1st array empty, sum is TT2
  int nr2 = TT2.rows();                         // Rows of array to add in
  if(!nr2) return TT1;                          // 2nd array empty, sum is TT1
  vector<int> aflag(nr2);			// Add flag
  int added=0, newtr=0;                         // Flag if to be added
  int i1, i2;                                   // For array row indexing
  for(i2=0; i2<nr2; i2++)                       // Loop rows of array to add
    {
    aflag[i2] = -1;
    for(i1=0, added=0; i1<nr1 && !added; i1++)  // Loop rows of array to add to
      if(fabs(TT1.getIm(i1,0)-TT2.getIm(i2,0)) <= res)
        {
        added++;
        aflag[i2] = i1;
        }
    if(!added) newtr++;
    }
  matrix mx(nr1+newtr, 2);
  for(i1=0; i1<nr1; i1++)                       // Loop rows of array to add
    {
    mx.put(TT1.get(i1,0), i1,0);
    mx.put(TT1.get(i1,1), i1,1);
    }
  int newrow = nr1;
  int oldrow;
  for(i2=0; i2<nr2; i2++)
    {
    oldrow = aflag[i2];
    if(oldrow >= 0)
      mx.put(mx.get(oldrow,1) + TT2.get(i2,1), oldrow,1);
    else
      {
      mx.put(TT2.get(i2,0), newrow,0);
      mx.put(TT2.get(i2,1), newrow,1);
      newrow++;
      }
    }
  TTable1D TT(mx);
  return TT;
// sosi need to set the other parameters to defaults
  }  



        // Input	TTab1D	: A transitions table (this)
 	//		k     	: A row/column index to sort
	//		type	: Basis for the sorting
	//			    0 = sort real values (default)
	//			   >0 = sort norms
	//		           <0 = sort imaginaries
	//		colf	: Row or column sort
	//			    0 = sort rows (default)
	//			   !0 = sort columns
	// Output	vector  : A vector of integers will be return that
	//			  contains indices sorted as specified
	// Note		      	: The matrix is left unaltered

// sosi - General for any matrix.  Put into class matrix when time permits.

vector<int> TTable1D::Sort(int k, int type, int colf) const
  {
//			Initial Setup For Sorting

  int      nvals = rows();		// Number of values to be sorted
  if(colf) nvals = cols();		// Sort the columns based on values
  double *vals;				// Array of values to be sorted
  vals = new double[nvals];		// Allocate the array of values
  vector<int> indx;			// Allocate vector to return
  if(type > 0)  type = 1;		// Insure type = { 0, 1, 2 }
  else if(type) type = 2;

//		Copy Column Or Row That Will Be Sorting Key

  int i=0,j=0;				// Array indices we'll use
  for(i=0; i<nvals; i++)		// First fill array vals with values 
    { 					// we want to sort & initialize the
    indx.push_back(i); 			// vector of sorted indices we return
    switch(type)
      {
      case 0:					// Sorting real values
      default:					// so copy either the specified
        if(colf) vals[i] = getRe(k,i);		// row or column real values
        else     vals[i] = getRe(i,k);
        break;
      case 1:					// Sorting norms
        if(colf) vals[i] = norm(get(k,i));	// so copy either the specified
        else     vals[i] = norm(get(i,k));	// row or column norm values
        break;
      case 2:					// Sorting imaginary values
        if(colf) vals[i] = getIm(k,i);		// so copy either the specified
        else     vals[i] = getIm(i,k);		// row or column imaginary vals
        break;
      }
    }

//	   Sort Column Or Row And Put Sorted Indices Into Return Vector

  double maxval;			// Value for maximum
  int maxind;				// Value of index at maximum
  for(i=0; i<nvals; i++)		// Now begin the sorting
    {					// by looping over the values
    maxval = vals[i];			// Assume this is maximum
    maxind = i;				// This is it's index
    for(j=i+1; j<nvals; j++)		// Now compare with all others
      {
      if(vals[j] > maxval)		// 	If we find a larger one
        {				//	then we set it as the maximum
        maxval = vals[j];		//	--> Store the bigger value
        maxind = j;			//	--> Store the big value index
        }
      }
    vals[maxind] = vals[i];		// Move value we started with to max
    j = indx[maxind];			// Here is the index of largest value
    indx[maxind] = indx[i];		// Copy current index to largest spot
    indx[i] = j;			// Copy index of largest to current spot
    }

//       Sorting Is Complete, Clean Up & Return Vector Of Sorted Indices

  delete [] vals;			// Delete the values
  return indx;				// Return vector of indices
  }

// ____________________________________________________________________________
// H       CLASS TRANSITIONS TABLE 1D WITH PARAMETERS & PARAMETER SETS
// ____________________________________________________________________________

/* These functions allow users to fill up transition tables from external
   paramters via a GAMMA parameter set.  The parameters can either be in an
   external ASCII file or already in a ParameterSet.  If a prefix has been set
   then the parameter names are taken to have a [#] before the names.

           Input		TTab1D	: A transitions table (this)
                                filein  : An input (ASCII) file
                                pset    : An input parameter files
                                indx    : Table index
                                warn    : Warning output level
                                            0 = no warnings
                                            1 = non-fatal warnings
                                            2 = fatal warnings
           Output               none    : Transitions table filled with
                                          transitions specified in filein

    Note that the work, i.e. the parameter parsing and table filling, is
    done in these 3 private functions: SetNTrans, SetConv, and SetTrans      */

bool TTable1D::readPSet(const string& filein, int indx, int warn)                       
  {
  ParameterSet pset;				// Declare a parameter set
  if(!pset.read(filein, warn?1:0))		// Read in pset from file
    {
    if(warn)
      {
      TTaberror(1, filein, 1);			// Problems with file filein
      if(warn>1) TTabfatality(50,filein);	// Can't read table
      else       TTaberror(50,filein,1);
      }
    return false;
    }  
  return readPSet(pset, indx, warn);
  }

bool TTable1D::readPSet(const ParameterSet& pset, int idx, int warn)
  {
  ParameterSet subpset;                         // Copy parameter set
  if(idx != -1) subpset = pset.strip(idx);      // to glean out those for
  else          subpset = pset;                 // the specified index [idx]
  if(!SetNTrans(subpset, warn)) return false;	// Get # of transitions
  TTable1D::SetConv(subpset, warn);		// Get type & conversion fac.
  if(!SetTrans(subpset, warn)) return false;	// Fill in the transitions
  return true;
  }

// ____________________________________________________________________________ 
// H                CLASS TRANSITIONS TABLE 1D I/O FUNCTIONS
// ____________________________________________________________________________ 

// ------------------------- ASCII Output Functions ---------------------------

/* Class TTable1D maintains both the transition frequencies and decay rates
   in units of 1/sec.  For output, these values are converted to appropriate
   units for display.  In particular, the half-height linewidths and T2
   relaxation times are related to the transition decay rate according to

                                     R     1
                           lwhh   = -- = -----
                               Hz   PI   PI*T
                                             2

   Several other internal flags will affect the appearance of the output 
   transitions table.  These may be set prior to output using the provided
   access functions.  The flags and access functions are listed below.

    Flag     Access Function               Values & Meaning
   -------   ---------------    ------------------------------------------
   FRQTYPE      setType         0:Hz, 1:PPM, 2:Gauss (output freq. units)
   FRQSORT      setSort         0:No sort, !0:Sort output by frequency
   FRQCONV      setConv         FRQTYPE=1: Om(MHz), FRQTYPE=2: g(unitless)
   FRQREV       set             T/F: Switch sign on frequency output
   ICUT         setIcut         #>0: Intensity cutoff value for output
   INORM        setInorm        #>0: Intensity noramlization for output
   SN           setSN           #>=0: Signal/Noise ratio, 0==no noise 
   HP           setHprint       !0=print header(default), 0=no header
   RP           setRprint       1: Flag print R values \ 
   LWP          setLWprint      1: Flag print linewidths\  <0 = never
   T2P          setT2print      1: Flag print T2 values /   0 = as needed
   PHP          setPHprint      1: Flag print phases   /    1 = always       */


void TTable1D::setType(int typ)
  {
  if(!CheckType(typ, 1))		// Checking the freq. output type
    {
    TTaberror(52, 1);			// Setting freq. output to Hz
    typ = 0;				// Set flag for freq. output in Hz
    }
  FRQTYPE = typ;			// Set the internal freq. type flag
  }

void TTable1D::setSort(int sf)     { sf?FRQSORT=1:FRQSORT=0;          }
void TTable1D::setConv(double cf)  { FRQCONV=cf;                      }
void TTable1D::setIcut(double ct)  { if(ct>0) ICUT = ct;              }
void TTable1D::setInorm(double in) { in?INORM=fabs(in)   :INORM=1;    }
void TTable1D::setSN(double S2N)   { SN = fabs(S2N);                  }
void TTable1D::setHprint(int hp)   { hp?HP=1             :HP=0;       }
void TTable1D::setRprint(int rp)   { rp?RP=rp/abs(rp)    :RP=0;       }
void TTable1D::setLWprint(int lwp) { lwp?LWP=lwp/abs(lwp):LWP=0;      }
void TTable1D::setT2print(int t2p) { t2p?T2P=t2p/abs(t2p):T2P=0;      }
void TTable1D::setPHprint(int php) { php?PHP=php/abs(php):PHP=0;      }
void TTable1D::setFreqRev()        { FRQREV?FRQREV=false:FRQREV=true; }

int    TTable1D::getType()    const { return FRQTYPE; }
int    TTable1D::getSort()    const { return FRQSORT; }
double TTable1D::getConv()    const { return FRQCONV; }
double TTable1D::getIcut()    const { return ICUT;    }
double TTable1D::getInorm()   const { return INORM;   }
double TTable1D::getSN()      const { return SN;      }
int    TTable1D::getHprint()  const { return HP;      }
int    TTable1D::getRprint()  const { return RP;      }
int    TTable1D::getLWprint() const { return LWP;     }
int    TTable1D::getT2print() const { return T2P;     }
int    TTable1D::getPHprint() const { return PHP;     }
bool   TTable1D::getFreqRev() const { return FRQREV;  }

// -------------------------- ASCII Output Functions --------------------------

        // Input		TTab1D	: A transitions table (this)
        //                      ostr    : An output stream
        // Output               ostr    : The output stream modified by
        //                                a listing of transitions
	// Note				: The value of FRQTYPE dictates how	
	//				  the transition frequencies are output
	//					0 = Hertz (default)
	//					1 = PPM
	//					2 = Gauss

std::vector<std::string> TTable1D::printStrings() const
  {
  std::vector<std::string> PStrings;
//		     First Decide Which Items Need Printing

  int r2p = 0;					// Flag to print R values
  int lwp = 0;					// Flag to print linewidths
  int t2p = 0;					// Flag to print T2 values
  int php = 0;					// Flag to print phases
  switch(RP)
    {
    case -1: r2p = 0;            break;		// Don't print relax. rates
    case  0: r2p = LineWidths(); break;		// Print relax. rates as needed
    case  1: r2p = 1;            break;		// Print relax. rates 
    }
  switch(LWP)
    {
    case -1: lwp = 0;            break;		// Don't print Linewidths
    case  0: lwp = LineWidths(); break;		// Print Linewidths as needed
    case  1: lwp = 1;            break;		// Print Linewidths 
    }
  switch(T2P)
    {
    case -1: t2p = 0;            break;		// Don't print T2 times
    case  0: t2p = LineWidths(); break;		// Print T2 times as needed
    case  1: t2p = 1;            break;		// Print T2 times 
    }
  switch(PHP)
    {
    case -1: php = 0;        break;		// Don't print phases
    case  0: php = Phases(); break;		// Print phases as needed
    case  1: php = 1;        break;		// Print phases 
    }
  
//		         Now We'll Print A Table Header

  string hdrS[6] = { " Frequency ", "Intensity",// Labels in the header
                     " Phase",      "Linewidth",
                     "T2(sec)",     "R2(1/sec)" };
  if(FRQTYPE==2 && FRQCONV)			// Use field lable if EPR
    hdrS[0] = string("   Field   ");		// and we know the conversion
  int i = 0;					// We'll resuse this index
  int lablens[6];				// Label string lengths
  for(; i<6; i++) lablens[i] = hdrS[i].length();// Store the string lengths
  int sprs[5] = { 6, 6, 7, 5, 5 };		// Spacer after header string 

  string  hdr =  hdrS[0];			// Add frequency lab. to header
          hdr += string(sprs[0],' ') + hdrS[1];	// Add intensity lab. to header
  if(php) hdr += string(sprs[1],' ') + hdrS[2];	// Add phase label to header
  if(lwp) hdr += string(sprs[2],' ') + hdrS[3];	// Add linewidth label to header
  if(t2p) hdr += string(sprs[3],' ') + hdrS[4];	// Add T2 label to header
  if(r2p) hdr += string(sprs[4],' ') + hdrS[5];	// Add R2 label to header
  int nblks = 40-hdr.length()/2;		// Number of blanks to center
  if(nblks<=0) nblks = 1;			// Insure this is positive
  string sline;
  sline = string(nblks, ' ') + hdr;		// Output 1st line of header
  PStrings.push_back(sline);

  string  hdr2 =  string(lablens[0], '-');	// Add frequency lab. underline
          hdr2 += string(sprs[0], ' ')		// Add intensity lab. underline
               +  string(lablens[1], '-');
  if(php) hdr2 += string(sprs[1], ' ')		// Add phase lab. underline
               +  string(lablens[2], '-');
  if(lwp) hdr2 += string(sprs[2], ' ')		// Add linewidth lab. underline
               +  string(lablens[3], '-');
  if(t2p) hdr2 += string(sprs[3], ' ')		// Add T2 lab. underline
               +  string(lablens[4], '-');
  if(r2p) hdr2 += string(sprs[4], ' ')		// Add R2 lab. underline
               +  string(lablens[5], '-');
  sline = string(nblks, ' ') + hdr2;		// Output 2nd line of header
  PStrings.push_back(sline);
//		   Finally We Print The Table Transitions

  vector<int> indx;				// Vector for sorted indices
  if(FRQSORT) indx = Sort(0, -1, 0);		// Get freq. sored indices
  double R, LW, T2, W, I, P;			// Variables for values printed
  int nrows = rows();				// Number of transitions
  int j=0;					// Transition index
  for(j=0; j<nrows; j++)			// Loop the transitions
    {
    if(!FRQSORT) i = j;				// Unsorted output
    else         i = indx[j];			// Sorted index
    I = norm(get(i,1));				// Get the intensity
    if(I >= ICUT)
      {

/* Begin line (1st column) by writing the frequency/field value.  The header
   label is 11 chars wide & we'll print the value using 10 characters followed
   by a space then a 3 character wide unit (which aligns with 3 chars of
   the 6 char wide space after the frequency header label.                   */

      W = getIm(i,0)/(PIx2);			// Frequency (Hz)
      sline =  string(nblks, ' ');		// Start a new line
      if(FRQTYPE==1 && FRQCONV)			// An NMR spectrum
        {					// with output in PPM
        W /= FRQCONV;				// Switch to PPM
        sline += Gform("%10.3f ppm", W);	// Output the frequency
        }
      else if(FRQTYPE==2 && FRQCONV)		// An EPR spectrum
        {					// with field units output
        W *= HZ2GAUSS/FRQCONV; 			// Switch to Gauss
        if(fabs(W) > 1.e11)     sline += Gform("%10.3f KT ", W*1.e-7);
        else if(fabs(W) > 1.e8) sline += Gform("%10.3f T  ", W*1.e-4);
        else if(fabs(W) < 0.1)  sline += Gform("%10.3f mG ", W*1.e3);
        else                    sline += Gform("%10.3f G  ", W);
        }
      else					// An NMR or EPR spectrum
        {					// with frequency units output
        if(fabs(W) > 1.e11)     sline += Gform("%10.3f GHz", W*1.e-9);
        else if(fabs(W) > 1.e8) sline += Gform("%10.3f MHz", W*1.e-6);
        else if(fabs(W) > 1.e5) sline += Gform("%10.3f KHz", W*1.e-3);
        else                    sline += Gform("%10.3f Hz ", W);
        }

/*  Now output the intensity, both value & header lable are 9 chars wide.
    We need a 3 char wide space though fill the columns after the.           */

      sline += string(3, ' ') + Gform("%9.3f", I/INORM); // Output intensity

/* Now output the phase if needed. The phase header is six characters wide
   and six unused colums (spaces) before the header.  We'll print 10 chars
   for the value so we begin with 2 blanks.                                  */

      if(php)
        {
        P = phase(get(i,1))*RAD2DEG;			// Get the phase (deg)
        if(fabs(P) < 1.e-10) P = 0.0;			// Set low ones to 0
        sline += string(2, ' ') + Gform("%10.3f", P);	// Print the phase
        }

/* Next we output the linewidths as desired.  The label is 9 characters wide
   and there are seven unused columns (spaces) before the header.  We'll
   use 10 chars to print the value so we need six preceding blanks.          */

      R = getRe(i,0);				// Here is the decay rate
      string SP(6, ' ');			// Spacer
      if(lwp)					// Output linewidths if needed
        {
        if(FRQTYPE==2)				// Assume its a EPR spectrum
          {					// and output lwpp
          LW = sqrt(4.0/3.0)*R;			//   Peak-peak linewidth (1/s)
          LW *= RAD2HZ;				//   Peak-peak linewidth (Hz)
          LW *= HZ2GAUSS;			//   Peak-peak linewidth (G)
          sline += SP + Gform("%10.3f G", LW); 	//   Ouput linewidth (G)
          }
        else					// Assume its an NMR spectrum
          { 					// and output lwhh
          LW = R/PI;				//   Half-height linewidth (Hz)
          sline += SP + Gform("%10.3f Hz", LW); //   Ouput linewidth (Hz)
          }
        }

/* Next we output the T2 relaxation times.  The label is 7 characters wide
   and there are five unused columns (spaces) before the header.  We'll
   use 10 chars to print the value so we need two preceding blanks.          */

      if(t2p)					// Output T2 times if needed
        {
        T2 = 1.0/R;				//	This is the T2 Time
        if(R) sline += "  " + Gform("%10.3f",T2);//	Output T2 time
        else  sline += "       inf  ";		// 	Output T2 time units
        }

/* Lastly we output the R2 relaxation rates.  The label is 9 characters wide
   and there are five unused columns (spaces) before the header.  We'll
   use 10 chars to print the value so we need four preceding blanks.         */

      if(r2p)						// Output R2 times if
        sline += string(4, ' ') + Gform("%10.3f", R);	// needed
      PStrings.push_back(sline);
      }
    }
  return PStrings;
  }

ostream& TTable1D::print(ostream& ostr) const
  {
  if(!rows())
    {
    if(HP)
      ostr << "\n\t1D Transitions Table: No Transitions Found\n";
    return ostr;
    }
  if(!Intensities())
    {
    if(HP)
      ostr << "\n\t1D Transitions Table: "
           << "No Transitions Found With Intensity Above " << ICUT << "\n";
    return ostr;
    }

//		     First Decide Which Items Need Printing

  int r2p = 0;					// Flag to print R values
  int lwp = 0;					// Flag to print linewidths
  int t2p = 0;					// Flag to print T2 values
  int php = 0;					// Flag to print phases
  switch(RP)
    {
    case -1: r2p = 0;            break;		// Don't print relax. rates
    case  0: r2p = LineWidths(); break;		// Print relax. rates as needed
    case  1: r2p = 1;            break;		// Print relax. rates 
    }
  switch(LWP)
    {
    case -1: lwp = 0;            break;		// Don't print Linewidths
    case  0: lwp = LineWidths(); break;		// Print Linewidths as needed
    case  1: lwp = 1;            break;		// Print Linewidths 
    }
  switch(T2P)
    {
    case -1: t2p = 0;            break;		// Don't print T2 times
    case  0: t2p = LineWidths(); break;		// Print T2 times as needed
    case  1: t2p = 1;            break;		// Print T2 times 
    }
  switch(PHP)
    {
    case -1: php = 0;        break;		// Don't print phases
    case  0: php = Phases(); break;		// Print phases as needed
    case  1: php = 1;        break;		// Print phases 
    }

//		         Now We'll Print A Table Header

  string hdrS[6] = { " Frequency ", "Intensity",// Labels in the header
                     " Phase",      "Linewidth",
                     "T2(sec)",     "R2(1/sec)" };

//setFreqRev

  if(FRQTYPE==2 && FRQCONV)			// Use field lable if EPR
    hdrS[0] = string("   Field   ");		// and we know the conversion
  int i = 0;					// We'll resuse this index
  int lablens[6];				// Label string lengths
  for(; i<6; i++) lablens[i] = hdrS[i].length();// Store the string lengths
  int sprs[5] = { 6, 6, 7, 5, 5 };		// Spacer after header string 

  string  hdr =  hdrS[0];			// Add frequency lab. to header
          hdr += string(sprs[0],' ') + hdrS[1];	// Add intensity lab. to header
  if(php) hdr += string(sprs[1],' ') + hdrS[2];	// Add phase label to header
  if(lwp) hdr += string(sprs[2],' ') + hdrS[3];	// Add linewidth label to header
  if(t2p) hdr += string(sprs[3],' ') + hdrS[4];	// Add T2 label to header
  if(r2p) hdr += string(sprs[4],' ') + hdrS[5];	// Add R2 label to header
  int nblks = 40-hdr.length()/2;		// Number of blanks to center
  if(nblks<=0) nblks = 1;			// Insure this is positive
  ostr << string(nblks, ' ') << hdr;		// Output 1st line of header

  string  hdr2 =  string(lablens[0], '-');	// Add frequency lab. underline
          hdr2 += string(sprs[0], ' ')		// Add intensity lab. underline
               +  string(lablens[1], '-');
  if(php) hdr2 += string(sprs[1], ' ')		// Add phase lab. underline
               +  string(lablens[2], '-');
  if(lwp) hdr2 += string(sprs[2], ' ')		// Add linewidth lab. underline
               +  string(lablens[3], '-');
  if(t2p) hdr2 += string(sprs[3], ' ')		// Add T2 lab. underline
               +  string(lablens[4], '-');
  if(r2p) hdr2 += string(sprs[4], ' ')		// Add R2 lab. underline
               +  string(lablens[5], '-');
  ostr << "\n" << string(nblks, ' ') << hdr2;	// Output 2nd line of header

//		   Finally We Print The Table Transitions

  vector<int> indx;				// Vector for sorted indices
  if(FRQSORT) indx = Sort(0, -1, 0);		// Get freq. sored indices
  double R, LW, T2, W, I, P;			// Variables for values printed
  int nrows = rows();				// Number of transitions
  int j=0;					// Transition index
  for(j=0; j<nrows; j++)			// Loop the transitions
    {
    if(!FRQSORT) i = j;				// Unsorted output
    else         i = indx[j];			// Sorted index
    I = norm(get(i,1));				// Get the intensity
    if(I >= ICUT)
      {

/* Begin line (1st column) by writing the frequency/field value.  The header
   label is 11 chars wide & we'll print the value using 10 characters followed
   by a space then a 3 character wide unit (which aligns with 3 chars of
   the 6 char wide space after the frequency header label.                   */

      W = getIm(i,0)/(PIx2);			// Frequency (Hz)
      if(FRQREV) W *= -1.0;
      ostr << "\n" << string(nblks, ' ');	// Start a new line
      if(FRQTYPE==1 && FRQCONV)			// An NMR spectrum
        {					// with output in PPM
        W /= FRQCONV;				// Switch to PPM
        ostr << Gform("%10.3f ppm", W);		// Output the frequency
        }
      else if(FRQTYPE==2 && FRQCONV)		// An EPR spectrum
        {					// with field units output
        W *= HZ2GAUSS/FRQCONV; 			// Switch to Gauss
        if(fabs(W) > 1.e11)     ostr << Gform("%10.3f KT ", W*1.e-7);
        else if(fabs(W) > 1.e8) ostr << Gform("%10.3f T  ", W*1.e-4);
        else if(fabs(W) < 0.1)  ostr << Gform("%10.3f mG ", W*1.e3);
        else                    ostr << Gform("%10.3f G  ", W);
        }
      else					// An NMR or EPR spectrum
        {					// with frequency units output
        if(fabs(W) > 1.e11)     ostr << Gform("%10.3f GHz", W*1.e-9);
        else if(fabs(W) > 1.e8) ostr << Gform("%10.3f MHz", W*1.e-6);
        else if(fabs(W) > 1.e5) ostr << Gform("%10.3f KHz", W*1.e-3);
        else                    ostr << Gform("%10.3f Hz ", W);
        }

/*  Now output the intensity, both value & header lable are 9 chars wide.
    We need a 3 char wide space though fill the columns after the.           */

      ostr << string(3, ' ') << Gform("%9.3f", I/INORM); // Output intensity

/* Now output the phase if needed. The phase header is six characters wide
   and six unused colums (spaces) before the header.  We'll print 10 chars
   for the value so we begin with 2 blanks.                                  */

      if(php)
        {
        P = phase(get(i,1))*RAD2DEG;			// Get the phase (deg)
        if(fabs(P) < 1.e-10) P = 0.0;			// Set low ones to 0
        ostr << string(2, ' ') << Gform("%10.3f", P);	// Print the phase
        }

/* Next we output the linewidths as desired.  The label is 9 characters wide
   and there are seven unused columns (spaces) before the header.  We'll
   use 10 chars to print the value so we need six preceding blanks.          */

      R = getRe(i,0);				// Here is the decay rate
      string SP(6, ' ');			// Spacer
      if(lwp)					// Output linewidths if needed
        {
        if(FRQTYPE==2)				// Assume its a EPR spectrum
          {					// and output lwpp
          LW = sqrt(4.0/3.0)*R;			//   Peak-peak linewidth (1/s)
          LW *= RAD2HZ;				//   Peak-peak linewidth (Hz)
          LW *= HZ2GAUSS;			//   Peak-peak linewidth (G)
          ostr << SP << Gform("%10.3f G", LW); 	//   Ouput linewidth (G)
          }
        else					// Assume its an NMR spectrum
          { 					// and output lwhh
          LW = R/PI;				//   Half-height linewidth (Hz)
          ostr << SP << Gform("%10.3f Hz", LW); //   Ouput linewidth (Hz)
          }
        }

/* Next we output the T2 relaxation times.  The label is 7 characters wide
   and there are five unused columns (spaces) before the header.  We'll
   use 10 chars to print the value so we need two preceding blanks.          */

      if(t2p)					// Output T2 times if needed
        {
        T2 = 1.0/R;				//	This is the T2 Time
        if(R) ostr << "  " << Gform("%10.3f",T2);//	Output T2 time
        else  ostr << "       inf  ";		// 	Output T2 time units
        }

/* Lastly we output the R2 relaxation rates.  The label is 9 characters wide
   and there are five unused columns (spaces) before the header.  We'll
   use 10 chars to print the value so we need four preceding blanks.         */

      if(r2p)						// Output R2 times if
        ostr << string(4, ' ') << Gform("%10.3f", R);	// needed
      }
    }
  ostr << "\n";
  return ostr;
  }

	// Input		TTab1D	: Transitions table
        //                      ostr 	: Output stream
        // Output               none	: Transitions table is sent
	//				  to the output stream

ostream& operator << (ostream& ostr, const TTable1D& TTab)
  { TTab.print(ostr); return ostr; }



        // Input		TTab1D	: A transitions table (this)
        //                      ostr    : An output stream
        // Output               ostr    : The output stream modified by
        //                                a view of the table status

ostream& TTable1D::status(ostream& ostr) const
  {
  ostr << "\n\t\t\t1D Transitions Table\n";
  ostr << "\n\t\tNumber of Transitions:   " << rows();
  ostr << "\n\t\tIntensity Cutoff Value:  " << ICUT;
  ostr << "\n\t\tIntensity Normalization: " << INORM;
  ostr << "\n\t\tIntensity % Cutoff Val:  " << PERCUT;
  ostr << "\n\t\tSignal To Noise Value:   " << SN;
  ostr << "\n\t\tFrequency Output Label:  "; 
  switch(FRQTYPE)
    {
    case 0:
    default: cout << "Hz";    break;
    case 1:  cout << "PPM";   break;
    case 2:  cout << "Gauss"; break;
    }
  ostr << "\n\t\tOutput Frequency Sort:   ";
  if(FRQSORT) cout << "Yes";
  else        cout << "No";
  if(FRQTYPE == 1 || FRQTYPE == 2) 
    {
    ostr << "\n\t\tFrequency Conversion:    ";
    ostr << FRQCONV;
    if(FRQTYPE == 1) ostr << " Hz/ppm";
    if(FRQTYPE == 2) ostr << " electron g factor";
    }
  ostr << "\n\t\tReverse Frequency Sign:  ";
  if(FRQREV) cout << "Yes";
  else       cout << "No";
  string ptypes[3] = { "(Do Not Print)", "(Print As Needed)", "(Print)" };
  ostr << "\n\t\tPhase Print Flag:        " << PHP << " " << ptypes[PHP+1];
  ostr << "\n\t\tLinewidth Print Flag:    " << LWP << " " << ptypes[LWP+1];
  ostr << "\n\t\tT2 Time Print Flag:      " << T2P << " " << ptypes[T2P+1];
  ostr << "\n\t\tRelax. Rate Print Flag:  " << RP  << " " << ptypes[RP+1];
  return ostr;
  }

// ------------------------- Binary Output Functions --------------------------

	// Input		TTab1D	: Transitions table (this)
        //                      fn	: Filename
        //                      fp	: File stream (pointing at TTab1D)
        // Return               void	: The transitions table is written
        //                                to a file called filename or into
        //				  file fp at current location
        // Note                       	: The file format is BINARY
 
void TTable1D::write(const string& fn) const
  {
  ofstream fp;					// Construct a file
  fp.open(fn.c_str(),ios::out|ios::binary);	// Open file
  write(fp);                                    // Write Op
  fp.close();                                   // Close file
  }
 
ofstream& TTable1D::write(ofstream& fp) const
  {
  matrix::write(fp);
  return fp;
  }
 

/* Added 09/03/09 by DCT */
void TTable1D::dbwrite_old( const string& fileName, 
		                    const string& compname,  // metabolite name			 
		                    const double& lowppm, 
					        const double& highppm, 
					        const double& specfreq, 
					        const double& reffreq,
					        const int& loop,
							const vector<string> & header) const
{
    ofstream ofs;					        // Construct a file
    ofs.open(fileName.c_str(), ios::out);	// Open file

	// write out header lines

	ofs << "; " << "\n";
	vector<string>::const_iterator vs;
	for(vs = header.begin(); vs != header.end(); vs++)
	{
		ofs << "; " << (*vs).c_str() << "\n";
	}
	ofs << "; " << "\n";


    // write the file output in the standard format for continued processing
	std::vector<double>	freqs;
	std::vector<int>	mx_index;
	std::vector<double>	freqout;
	std::vector<double>	ampout;
	std::vector<double>	phaseout;

	double			freq;
	double			normal = 1.0;
	int			    normindex = 0;
	double			amptemp, ampsum, phasetemp;
	unsigned long	bincount = 0;
	const double	freqtol = 0.1/specfreq;		// Use something like half the minimum
	const double	phasetol = 50.0;				// coupling const. divided by field
	unsigned long	i,k;						// strength for freqtol
	int				foundone = 0;
	double			refhigh, reflow;

	// this is 0 to .2 ppm for protons

	refhigh = reffreq + 0.1;
	reflow  = reffreq - 0.1;

	// cout << "reflow " << reflow << endl;
	// cout << "refhigh " << refhigh << endl;

	/* Get index array in PPM order -------------------------------------*/

	mx_index = this->Sort(0,-1,0);     

	
	/* Convert frequencies to PPM ---------------------------------------*/
	
	for(i=0; i<static_cast<unsigned long>(this->size()); i++)      
	{
		freqs.push_back(-this->Fr(mx_index[i])/(2.0*PI*specfreq));
	}

	/* Sum over the normalizer peak -------------------------------------*/

	std::vector<double>::iterator itf;

	for(i=0, itf=freqs.begin(); i<static_cast<unsigned long>(this->size()); i++, itf++)
	{
		freq = *itf;
		// cout << freq << endl;
		if (freq > reflow && freq < refhigh) 
		{
			if (normindex == 0) 
			{
				normal =  norm(this->I(mx_index[i]));
				++normindex;
			}
			else 
			{
				normal += norm(this->I(mx_index[i])); 
			}
		}
	}

	std::cout << "Normalization in dbwrite: " << normal << " with " << normindex << " line indices\n";
	/* Simple peak blending based on Freqtol and Phasetol ---------------*/
	
	for(i=0, itf=freqs.begin(); i<static_cast<unsigned long>(this->size()); i++, itf++)
	{
		freq = *itf;
		
		if ((freq > lowppm && freq < reflow) || (freq > refhigh && freq < highppm)) 
		{
			amptemp   =  norm(this->I(mx_index[i]))/normal;
			phasetemp = -RAD2DEG*phase(this->I(mx_index[i]));
			
			if (bincount == 0) 
			{
				freqout.push_back(freq);
				ampout.push_back(amptemp);
				phaseout.push_back(phasetemp);
				++bincount;
			}
			else 
			{  
				for(k=0; k < bincount && foundone == 0; k++) 
				{
					if(freq >= freqout[k]-freqtol && freq <= freqout[k]+freqtol)
					{
						if (phasetemp >= phaseout[k]-phasetol && phasetemp <= phaseout[k]+phasetol) 
						{
							ampsum      =  ampout[k]+amptemp;
							freqout[k]  = (ampout[k]*freqout[k]  + amptemp*freq)/ampsum;
							phaseout[k] = (ampout[k]*phaseout[k] + amptemp*phasetemp)/ampsum;
							ampout[k]  +=  amptemp;
							foundone = 1;
						}
					}
				}
				if (foundone == 0) 
				{
					freqout.push_back(freq);
					ampout.push_back(amptemp);
					phaseout.push_back(phasetemp);
					++bincount;
				}
				foundone = 0;
			}
		}
	}

	std::vector<double>::iterator Fr, Am, Ph;

	for(i=0, Fr=freqout.begin(), Am=ampout.begin(), Ph=phaseout.begin(); i<bincount; i++, Fr++, Am++, Ph++)
	  {
		ofs << compname;					// metabolite name
		ofs << "\t" <<  loop;				// loop number
		ofs << "\t" <<  0;					// group number in metabolite
		ofs << "\t" <<  i;					// line index in metabolite
		ofs << "\t" <<  *Fr;				// frequency in ppm
		ofs << "\t" <<  *Am;				// intensity
		ofs << "\t" <<  ((*Ph) * (-1));		// phase
		ofs << "\n";		
	  }

    ofs.close();     // Close file
}


// dbwrite.cpp ----------------------------------------------------------- //



// #include <iostream>
// using namespace std;
// modified by LK, Dec. 2007
// imported into TTable1D with slight modifications on interface
// by DCT 10/07/09

void TTable1D::dbwrite( const string& fileName, 
                        const std::string& compname, 
                        const double& specfreq,
                        const int& numberspins,
                        const int& loop,
                        const vector<string> & header) const
// Input
// Note       : Frequencies and Rates are in 1/sec

{

    // Added to LK's version so that only the filename needs 
    // to be passed in (no iostream type). This will work better
    // in the python implementation.  DCT (10/06/09).
    ofstream ostr;					        // Construct a file
    ostr.open(fileName.c_str(), ios::out);	// Open file

    // write out header lines

    ostr << "; " << "\n";
    vector<string>::const_iterator vs;
    for(vs = header.begin(); vs != header.end(); vs++)
    {
        ostr << "; " << (*vs).c_str() << "\n";
    }
    ostr << "; " << "\n";


    std::vector<double>    freqs;
    std::vector<int>       mx_index; 
    std::vector<double>    freqout;
    std::vector<double>    ampout;
    std::vector<double>    phaseout;

    double           freq;
    double           normal = 1.0;
    double           amptemp, ampsum, phasetemp;
    unsigned long    bincount = 0;
    const double     freqtol = 0.1/specfreq;    // Use something like half the minimum
    const double     phasetol = 50.0;           // coupling const. divided by field
    unsigned long    k;                         // strength for freqtol
    int              foundone = 0;

    int ns=numberspins;


    /* Get index array in PPM order -------------------------------------*/

    mx_index = this->Sort(0,-1,0);     
       

    /* Convert frequencies to PPM ---------------------------------------*/

    for(long ii=0; ii<this->size(); ii++)
    {
        freqs.push_back(-this->Fr(mx_index[ii])/(2.0*PI*specfreq));
    }

    // there is no more reference spins, normalization is taken
    // place by the formula below, where ns is the number of spins
    normal=pow(2.0,(ns-1));
    normal=normal/2;

    //std::cout << "normal is " << normal<< std::endl;
    //std::cout << "numberofspins is " << ns<< std::endl;

    // DCT - added on 11/18/09
    std::cout << "Normalization in dbwrite: " << normal << " \n";

    /* Simple peak blending based on Freqtol and Phasetol ---------------*/

    {
        std::vector<double>::iterator itf;
        long i;

        for(i=0, itf=freqs.begin(); i<this->size(); i++, itf++)
        {
		    freq = *itf;
		   
			  //LK no need for low and high    
			  //if ((freq > lowppm && freq < reflow) || (freq > refhigh && freq < highppm))

		    {
		        amptemp   =  norm(this->I(mx_index[i]))/normal;
		        phasetemp = -RAD2DEG*phase(this->I(mx_index[i]));

	   			 // std::cout << "phase is " << phase(this->I(mx_index[i]))<< std::endl;
		       		       
		        if (bincount == 0) //this is to write out the very first set of numbers
		        {
		            freqout.push_back(freq);
		            ampout.push_back(amptemp);
		            phaseout.push_back(phasetemp);
		            ++bincount;
		        }
		        else
		        { 
		            for(k=0; k < bincount && foundone == 0; k++) //this is equivalent to while loop in matlab
		            {
		                if(freq >= freqout[k]-freqtol && freq <= freqout[k]+freqtol)
		                {
		                    if (phasetemp >= phaseout[k]-phasetol && phasetemp <= phaseout[k]+phasetol)
		                    {
		                        ampsum      =  ampout[k]+amptemp;
		                        freqout[k]  = (ampout[k]*freqout[k]  + amptemp*freq)/ampsum;
		                        phaseout[k] = (ampout[k]*phaseout[k] + amptemp*phasetemp)/ampsum;
		                       
		                        ampout[k]  +=  amptemp;
		                        foundone = 1;                           
		                    }
		                }
		            }
		            if (foundone == 0)
		            {
		                freqout.push_back(freq);
		                ampout.push_back(amptemp);
		                phaseout.push_back(phasetemp);
		                ++bincount;
		            }
		            foundone = 0;
		        }//end of else   
		  }
		}
	}


    /* Output values to text file --------------------------------------------*/


	{
		std::vector<double>::iterator Fr, Am, Ph;
		unsigned long i;
		char tempFr[100];
		char tempAm[100];
		char tempPh[100];

		for(i=0, Fr=freqout.begin(), Am=ampout.begin(), Ph=phaseout.begin(); i<bincount; i++, Fr++, Am++, Ph++)
		{
		    ostr << compname;                    // metabolite name
		    ostr << "\t" <<  loop;                // loop number
		    ostr << "\t" <<  0;                    // group number in metabolite
		    ostr << "\t" <<  i;                    // line index in metabolite
#ifdef _MSC_VER
			sprintf_s(tempFr, "%.8f", *Fr);
			sprintf_s(tempAm, "%.8f", *Am);
		    sprintf_s(tempPh, "%.8f", ((*Ph) * (-1.0)) );
#else
			sprintf(tempFr, "%.8f", *Fr);
			sprintf(tempAm, "%.8f", *Am);
			sprintf(tempPh, "%.8f", ((*Ph) * (-1.0)) );
#endif
		    ostr << "\t" << tempFr;     // frequency in ppm
			ostr << "\t" << tempAm;     // intensity
		    ostr << "\t" << tempPh ;    // phase
		    ostr << "\n";       
		}
        //ostr << "bincount = " << bincount << "\n";
	}

    ostr.close();     // Close file
}


unsigned int    TTable1D::calc_spectra( vector<double> & frq,
                                        vector<double> & amp,
                                        vector<double> & phz,
                                        double specfreq,
                                        int numberspins,
                                        double freqtol,
                                        double phasetol,
                                        double lowppm,
                                        double highppm) const

// Note       : Frequencies and Rates are in 1/sec
{
    vector<double>    freqs;
    vector<int>       mx_index;

    vector<double> freqout;
    vector<double> ampout;
    vector<double> phaseout;

    double           freq;
    double           normal = 1.0;
    double           amptemp, ampsum, phasetemp;
    unsigned long    bincount = 0;
    unsigned long    k;                        // strength for freqtol
    int              foundone = 0;

    int ns=numberspins;


    /* Get index array in PPM order -------------------------------------*/
    mx_index = this->Sort(0,-1,0);     
       
    /* Convert frequencies to PPM ---------------------------------------*/
    for(long ii=0; ii<this->size(); ii++)
    {
        freqs.push_back(-this->Fr(mx_index[ii])/(2.0*PI*specfreq));
    }

    // there is no more reference spins, normalization is taken
    // place by the formula below, where ns is the number of spins
    normal = pow(2.0,(ns-1));
    normal = normal/2;

    /* ----- Simple peak blending based on Freqtol and Phasetol -----*/
    {
        std::vector<double>::iterator itf;
        long i;

        for(i=0, itf=freqs.begin(); i<this->size(); i++, itf++)
        {
            freq = *itf;  
            
            if ((freq > lowppm ) && ( freq < highppm))
            {
                amptemp   =  norm(this->I(mx_index[i]))/normal;
                phasetemp = -RAD2DEG*phase(this->I(mx_index[i]));		       
               
                if (bincount == 0) //this is to calculate the very first set of numbers
                {
                    freqout.push_back(freq);
                    ampout.push_back(amptemp);
                    phaseout.push_back(phasetemp);

                    bincount += 1;
                }
                else
                { 
                    for(k=0; k < bincount && foundone == 0; k++) //this is equivalent to while loop in matlab
                    {
                        if(freq >= freqout[k]-freqtol && freq <= freqout[k]+freqtol)
                        {
                            if (phasetemp >= phaseout[k]-phasetol && phasetemp <= phaseout[k]+phasetol)
                            {
                                ampsum      =  ampout[k]+amptemp;
                                freqout[k]  = (ampout[k]*freqout[k]  + amptemp*freq)/ampsum;
                                phaseout[k] = (ampout[k]*phaseout[k] + amptemp*phasetemp)/ampsum;
                                ampout[k]  +=  amptemp;

                                foundone = 1;                           
                            }
                        }
                    }
                    if (foundone == 0)
                    {
                        freqout.push_back(freq);
                        ampout.push_back(amptemp);
                        phaseout.push_back(phasetemp);

                        bincount += 1;
                    }

                    foundone = 0;
                } //end of else   
            }
        }
    }
    

    frq.resize(bincount);
    amp.resize(bincount);
    phz.resize(bincount);

    for(unsigned int i=0; i<bincount; i++)
    {
        frq[i] = freqout[i];
        amp[i] = ampout[i];
        phz[i] = -1.0*phaseout[i];
    }

    return bincount;
}


unsigned int    TTable1D::calc_spectra2(vector<double> & frq,
                                        vector<double> & amp,
                                        vector<double> & phz,
                                        double specfreq,
                                        int numberspins,
                                        double freqtol,
                                        double phasetol,
                                        double lowppm,
                                        double highppm,
                                        double normal,
                                        bool verbose) const

// Note       : Frequencies and Rates are in 1/sec
{
    vector<double>    freqs;
    vector<int>       mx_index;

    vector<double> freqout;
    vector<double> ampout;
    vector<double> phaseout;

    double           freq;
    double           normal_v1 = 1.0;
    double           amptemp, ampsum, phasetemp;
    unsigned long    bincount = 0;
    unsigned long    k;                        // strength for freqtol
    int              foundone = 0;

    int ns=numberspins;


    /* Get index array in PPM order -------------------------------------*/
    mx_index = this->Sort(0,-1,0);     
       
    /* Convert frequencies to PPM ---------------------------------------*/
    for(long ii=0; ii<this->size(); ii++)
    {
        freqs.push_back(-this->Fr(mx_index[ii])/(2.0*PI*specfreq));
    }

    // there is no more reference spins, normalization is taken
    // place by the formula below, where ns is the number of spins
    normal_v1 = pow(2.0,(ns-1));
    normal_v1 = normal_v1/2;

    if (normal == 0.0) normal = normal_v1;      // assume external norm NOT calculated

    if (verbose)
    {
        std::cout << "normal v1: " << normal_v1 << " \n";
        std::cout << "normal   : " << normal    << " \n";
    }

    /* ----- Simple peak blending based on Freqtol and Phasetol -----*/
    {
        std::vector<double>::iterator itf;
        long i;

        for(i=0, itf=freqs.begin(); i<this->size(); i++, itf++)
        {
            freq = *itf;  
            
            if ((freq > lowppm ) && ( freq < highppm))
            {
                amptemp   =  norm(this->I(mx_index[i]))/normal;
                phasetemp = -RAD2DEG*phase(this->I(mx_index[i]));		       
               
                if (verbose)
                {
                    std::cout << "fre, amp, pha: " << freq << ", " << amptemp << ", " << phasetemp << " \n";
                }               
               
                if (bincount == 0) //this is to calculate the very first set of numbers
                {
                    freqout.push_back(freq);
                    ampout.push_back(amptemp);
                    phaseout.push_back(phasetemp);

                    bincount += 1;
                }
                else
                { 
                    for(k=0; k < bincount && foundone == 0; k++) //this is equivalent to while loop in matlab
                    {
                        if(freq >= freqout[k]-freqtol && freq <= freqout[k]+freqtol)
                        {
                            if (phasetemp >= phaseout[k]-phasetol && phasetemp <= phaseout[k]+phasetol)
                            {
                                ampsum      =  ampout[k]+amptemp;
                                freqout[k]  = (ampout[k]*freqout[k]  + amptemp*freq)/ampsum;
                                phaseout[k] = (ampout[k]*phaseout[k] + amptemp*phasetemp)/ampsum;
                                ampout[k]  +=  amptemp;

                                foundone = 1;                           
                            }
                        }
                    }
                    if (foundone == 0)
                    {
                        freqout.push_back(freq);
                        ampout.push_back(amptemp);
                        phaseout.push_back(phasetemp);

                        bincount += 1;
                    }

                    foundone = 0;
                } //end of else   
            }
        }
    }
    

    frq.resize(bincount);
    amp.resize(bincount);
    phz.resize(bincount);

    for(unsigned int i=0; i<bincount; i++)
    {
        frq[i] = freqout[i];
        amp[i] = ampout[i];
        phz[i] = -1.0*phaseout[i];
    }

    return bincount;
}


// -------------------------- Binary Input Functions --------------------------
 
	// Input		TTab1D	: Transitions table (this)
        //                      fn	: Filename
        //                      fp	: File stream
        // Return               void	: Table TTab1D is read from
        //                                a binary file called filename or from
	//				  a binary filestream called fp
	// Note				: Read is performed in filestream
        //                                fp at current location or from start
	//				  of file fn

void TTable1D::read(const string& fn)
  {
  ifstream fp;  				// Construct a file
  fp.open(fn.c_str(),ios::in|ios::binary);	// Open file
  read(fp);					// Read TTab1D, use overload
  fp.close();					// Close file
  }
 
ifstream& TTable1D::read(ifstream &fp)
  { matrix::read(fp); return fp; }

// ____________________________________________________________________________ 
// I         CLASS TRANSITIONS TABLE 1D DEBUGGING HELPER FUNCTIONS
// ____________________________________________________________________________ 

	// Input		TTab1D : Transitions table (this)
	//			tinc	: Dwell time in seconds
	//			npts	: Number of points
	//			P2P     : Computation control
	//				    0 = Evolve From t=0 To Point 
	//				   !0 = Evolve From Point To Point 
        // Output               ostr	: The output stream modified by
        //                                the acquisition time domain
	//				  calculation with debugging info
 
std::ostream& TTable1D::printT(std::ostream& ostr,double tinc,int npts,int P2P)
  {
  std::string hdr("Time Domain Spectrum Calculation");
  int hl = hdr.length();
  int cl = 40-hl/2;
  ostr << string(cl,' ') << hdr << "\n";
  ostr << string(cl,' ') << string(hl,'=') 
       << "\n";
  int ntr = rows();				// Get number of transitions
  if(!ntr)					// If no transitions we cannot
    {						// calculate any spectrum
    hdr = "No Transitions Present, No Spectrum";
    cl = 40 - hdr.length()/2;
    ostr << "\n" << string(cl, ' ') << hdr;
    return ostr;
    }

//	 First Check Which Transitions Are Used & Their Point Range

  std::vector<std::string> hdrs;
  hdrs.push_back(std::string("Index"));
  hdrs.push_back(std::string("        Base Intensity      "));
  hdrs.push_back(std::string("  Frequency "));
  hdrs.push_back(std::string("    Rate    "));
  hdrs.push_back(std::string("Range"));
  hdrs.push_back(std::string("Utilized")); 

  std::string csp("  ");
  int i,j,k;
  hdr = std::string("");
  for(i=0; i<int(hdrs.size())-1; i++)
    hdr += hdrs[i] + csp;
  hdr += hdrs[i];
  cl = 40 - hdr.length()/2;
  ostr << "\n" << string(cl, ' ') << hdr;
  hdr = std::string("");
  for(i=0; i<int(hdrs.size())-1; i++)
    hdr += std::string(hdrs[i].length(), '-') + csp;
  hdr += std::string(hdrs[i].length(), '-');
  cl = 40 - hdr.length()/2;
  ostr << "\n" << string(cl, ' ') << hdr;
  std::string st(cl, ' ');
  double rval, ival;
  

  int *ifi;					// Final indices, each tr
  ifi = new int[ntr];
  Exponen_cut(ifi,*this,tinc,npts,PERCUT);	// Fill index array ifi

  int trlen = Gdec(ntr).length();
  for(k=0; k<ntr; k++)				// Now loop the transitons
    {						// and output how they would
    ostr << "\n" << st;
    ostr << Gdec(k,trlen) << "." << csp;	// Output the transition index

    rval = getRe(k,1);
    ival = getIm(k,1);
    ostr << Gform("%12.5f", rval);
    if(ival<0) ostr << " - ";
    else       ostr << " + ";
    ostr << Gform("%12.5f", fabs(ival));
    ostr << "i" << csp;

    rval = getRe(k,0);				// Output the frequency
    ival = getIm(k,0);				// and linewidth
    ostr << Gform("%12.5f", ival) << csp;
    ostr << Gform("%12.5f", rval) << csp;

    if(ifi[k]<0) ostr << " All " << csp;	// Output last point the
    else         ostr << Gdec(ifi[k], 5) << csp;// transition contributes to

    if(norm(get(k,1)) > ICUT)			// Yes/No if transition has
      ostr << "   Yes  ";			// intensity above cutoff
    else
      ostr << "   No   ";
    }



  if(P2P) hdr = "Point To Point Computation [t-1 => t]";
  else    hdr = "Start To Point Computation [t0 => t]";
  cl = 40 - hdr.length()/2;
  ostr << "\n\n" << string(cl, ' ') << hdr;
  ostr << "\n"   << string(cl, ' ') << string(hdr.length(), '-');

  complex *expwt=NULL, *dxi=NULL, *exptit=NULL;
  if(P2P)
    {
    dxi = new complex[ntr];
    expwt = new complex[ntr];
    exptit = new complex[ntr];
    for(int tr=0; tr<ntr; tr++)
      {
      expwt[tr] = (get(tr,0)*tinc).Zexp();	// Incremental evolutions
      dxi[tr] = get(tr,1);			// Transition intensities
      exptit[tr] = complex1;			// Track exponentials
      }
    }

//	      This Calculation Is Performed From Start To Point

  int pl = (Gdec(npts)).length();

  std::string RS("\n\n\n");
  std::string CAP("Calculation Of Point ");

//            Header For Individual Transition Constributions

  hdrs.clear();
  hdrs.push_back(std::string("   Exp[(-iW-R)t]    "));
  hdrs.push_back(std::string("        Contribution        "));
  hdrs.push_back(std::string("      Summed Intensity      "));

  std::string AT(" At Time ");
  std::string TRN(trlen+2, ' ');
  std::string TRU(trlen+2, '-');

  complex Itr, z, expatit, ptval, Icont;	// Transition tr intensity
  double time = 0;				// Track point time
  for(i=0; i<npts; i++)				// Loop over all desired pts
    {
    ostr << RS    << CAP << Gdec(i,pl) << AT << time << "\n";
    ostr << "\n" << TRN << " ";
    for(j=0; j<int(hdrs.size()); j++) ostr << hdrs[j] << csp;
    ostr << "\n" << TRN << " ";
    for(j=0; j<int(hdrs.size()); j++) ostr << string(hdrs[j].length(), '-') << csp;
    ptval = complex0;				//   Set this point to zero
    if(!P2P)
      {
      for(int tr=0; tr<ntr; tr++)		// Loop over all transitions
        {
        z = complex1;				//   Exponential factor is 1

        ostr << "\n" << Gdec(tr,trlen) << "." << csp;

        expatit = (get(tr,0)*time).Zexp();	//    This is the precessed
        rval = zRe(expatit);			//    exponential of transition
        ival = zIm(expatit);			//    frequency/linewith for
        ostr << Gform("%8.5f", rval);
        if(ival<0) ostr << " - ";
        else       ostr << " + ";
        ostr << Gform("%8.5f", fabs(ival));
        ostr << "i" << csp;

        if(norm(get(tr,0)) > ICUT)		//    Use the transition only 
          { 					//    if the intensity is O.K.
          Icont = get(tr,1)*expatit;		//    This is the contribution
          rval = zRe(Icont);			//    to the point intensity from
          ival = zIm(Icont);			//    this transition. Output its
          ostr << Gform("%12.5f", rval);
          if(ival<0) ostr << " - ";
          else       ostr << " + ";
          ostr << Gform("%12.5f", fabs(ival));
          ostr << "i" << csp;
       
          ptval += Icont;			//    Add precessed transition
          rval = zRe(ptval);			//    to the total intensity of
          ival = zIm(ptval);			//    the point and output its
          ostr << Gform("%12.5f", rval);
          if(ival<0) ostr << " - ";
          else       ostr << " + ";
          ostr << Gform("%12.5f", fabs(ival));
          ostr << "i";
          }
        }
      }
    else
      {
      for(int tr=0; tr<ntr; tr++)		// Loop over all transitions
        {
        Itr = get(tr,1);			//   Intensity of transition tr
        ostr << "\n\ttr " << tr << ". I="	//   Output transition index,
             << Itr << ", W=" << get(tr,0)	//   intensity, frequency,
             << ", E=" << expwt[tr];		//   and exponential increment
        if(norm(Itr) > ICUT)			//   Use trans. if intensity OK
          {
          ptval += dxi[tr];			//     Add trans. contribution
          ostr << ", e=" << exptit[tr];		//     Output exponent. advanc2
          ostr << ", c=" << Itr*exptit[tr] 	//     Output contribution and
               << ", pt=" << ptval;		//     current point state
          }
        dxi[tr] *= expwt[tr];			//   Advance point to next time
        exptit[tr] *= expwt[tr];		//   Track the exponentials
        }
      }
      time += tinc;				// Increment the point time
      ostr << "\n\nFinal Intensity Value (x10^10): "	// Output point value
           << ptval*1.e10; 

    }

  if(P2P)
    {
    delete [] dxi;
    delete [] expwt;
    delete [] exptit;
    }
 
  delete [] ifi;
  return ostr;
  }

 
ostream& TTable1D::printF(ostream& ostr, int npts, double Fst, double Ffi)

	// Input		TTab1D : Transitions table (this)
	//			npts	: Number of points
	//			Fst	: Initial frequency
	//			Ffi	: Final frequency
        // Output               ostr	: The output stream modified by
        //                                the acquisition frequency domain
	//				  calculation with debugging info

  {
  int ntr = rows();				// Get number of transitions
  if(!ntr) return ostr;				// Can't print, no transitions

//	 First Check Which Transitions Are Used & The Point Range

  double wo =   PIx2*Fst;			// Init. frequency (1/seconds)
  double winc = PIx2*(Ffi-Fst)/(npts-1);	// Frequency between pts(1/sec)
  int *ist, *ifi;			// Init, final indices, each tr
  ist = new int[ntr];
  ifi = new int[ntr];
  Lorentz_cut(ist,ifi,*this,wo,winc,npts,PERCUT);// Fill index arrays ist & ifi
  int *Lint;				// Flags if integration needed
  Lint = new int[ntr];
  Lorentzian(Lint,*this,winc);			// Fill index arrays ist & ifi

  for(int k=0; k<ntr; k++)			// Now loop the transitons
    {						// and output how they would
    ostr << "\n Trans. " << k << ". "; 		// be used in makeing the
    if(norm(get(k,1)) > ICUT)			// time domain spectrum
      {
      ostr << "USE     ";
      ostr << ", Range: " << ist[k]
           << "," << ifi[k];
      if(Lint[k]) ostr << " INTEGRATE";
      else        ostr << " NO INTEGRAL";
      }
    else ostr << "NEGLECT ";
    }


/*
  double w, dw, dw2;				// Frequency, difference freq.
  double Rtr, Rtr2, Wtr;			// Transition tr rate & freq.
  complex Itr;					// Transition tr intensity
  double dwa,dwb,re,im;
  int i;					// Point indices
*/
  delete [] ist;
  delete [] ifi;
  delete [] Lint;
  return ostr;
  }


// ____________________________________________________________________________ 
// J         CLASS TRANSITIONS TABLE OLD GAMMA SUPPORT FUNCTIONS
// ____________________________________________________________________________ 

// sosi
 
void offset(matrix& mx, double F, double LWR, int inHz)
 
        // Input        mx    : Table of transition information
        //                      Re(<i|mx|0>) = transition i decay rate
        //                      Im(<i|mx|0>) = transition i frequency
        //                         <i|mx|1>  = transition i intensity
        //              F     : Frequency offset (Hz)
        //              LWR   : Linewidth/Rate adjustment (Hz)
        //              inHz  : Flag whether LWR is a linewidth (Hz)
        //                      or a rate (1/sec) <-- default
        // Note         mx    : Matrix values adjusted
 
// GAMMA's Lorentzian relationship between the linewidth at half-height (lwhh)
// and its R value (a decay rate of its related exponential) according to
 
//                    R     1
//          lwhh   = -- = -----  ---> lwhh  = 2R   &  lwhh     = 2R
//              Hz   PI   PI*T            Hz    Hz        1/sec    1/sec
//                            2  
//                     -1
// where R & T  are sec  , lwhh in Hz. Switching both lwhh & R into the same
//            2                        units implies the relationship lw = 2R.
 
  {
  TTable1D tmp;
  string fct("offset");
  tmp.TTaberror(4, fct, 1);
  tmp.TTaberror(5, fct);
  double W = F*2.0*PI;                          // Offset in 1/sec
  double R = LWR;                               // Rate in 1/sec
  if(inHz) R *= PI;                             // Adjust if lwhh(Hz) input
  complex zadd(R,W);
  complex z;
  for(int i=0; i<mx.rows(); i++)
    {
    z = mx.get(i,0) + zadd;
    mx.put(z, i, 0);
    }
  }

#endif							// TrnsTable1D.cc
