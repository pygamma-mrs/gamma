/* FrameMakerP.cc ***********************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**	FrameMaker Parameters 	                  Implementation   	**
**                                                                      **
**      Copyright (c) 1997                                              **
**      Scott A. Smith                                                  **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The Class FMPar is a simple data type to facility the passing of     **
** parameters to GAMMA FrameMaker functions.                            **
**                                                                      **
*************************************************************************/

#ifndef   FrameMakerP_cc_			// Is file already included?
#  define FrameMakerP_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#  endif

#include <GamIO/FrameMakerP.h>			// Include the header
#include <Basics/Gconstants.h>			// Include PI
#include <Basics/Gutils.h>			// Include GAMMA errors
#include <iostream>				// Include input output streams
#include <string>				// Include libstdc++ STL strings
#include <cmath>				// Include HUGE_VAL

using std::string;				// Using libstdc++ strings
using std::ostream;				// Using libstdc++ output streams
using std::cout;				// Using libstdc++ standard output

const double FMPAGEWIDTH  = 19.0;		// Page width 19.05cm = 7.5in
const double FMPAGEHEIGHT = 22.8;		// 22.86 cm = 9in
const int    FMPLmax = 128;			// Default polyline size

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
//                       CLASS FMPARAMETERS ERROR HANDLING
// ____________________________________________________________________________

                                                                                
void FMPar::FMPerror(int eidx, int noret) const
 
        // Input                FMP	: FrameMaker Parameters (this)
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
 
  {
  cout << "\nFrameMaker Parameters: ";
  switch(eidx)
    {
    case 0:                                                             // (0)
      cout << "Program Aborting.....";
      break;
    case 50:                                                            // (50)
      cout << "Plotted Array Contains < 2 Rows!";
      break;
    case 51:                                                            // (51)
      cout << "Cannot Produce A Stack Plot.";
      break;
    }
  if(!noret) cout << ".\n";
  }

 
volatile void FMPar::FMPfatality(int eidx) const

        // Input                FMP	: FrameMaker Parameters (this)
        //                      eidx    : Flag for error type
        // Output               none    : Error message output
        //                                Program execution stopped

  {                                                                 
  FMPerror(eidx);				// Output error message
  if(eidx) FMPerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }
 



// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                   FMPARAMETERS CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________


FMPar::FMPar ( )
 
        // Input                none    : No arguments required
        // Output               none    : FrameMaker parameters
 
  {
  hsize = FMPAGEWIDTH;			// FM horizontal plot size (cm)
  vsize = FMPAGEHEIGHT;			// FM vertical plot size   (cm)
  haxmin = 0.0;				// Horizontal axis label start
  haxmax = 1.0;				// Horizontal axis label finish
  vaxmin = 0.0;				// Vertical axis label start
  vaxmax = 1.0;				// Vertical axis label finish
  duse = plotreals;			// Assume plot reals
  PLmax = FMPLmax;			// Default polyline size
  VRES = 0.0001;			// Vertical resolution (1 um)
  HRES = 0.0001;			// Horizontal resolution (1 um)
  FMdebug = 0;				// No debugging
  }


FMPar::~FMPar ( ) {}

        // Input                LOp  : A super operator (this)
        // Output               none : Destroys super operator LOp


void FMPar::operator=(const FMPar& FMP1)

        // Input                FMP     : FrameMaker parameters (this)
        //                      FMP1    : FrameMaker parameters
        // Output               void    : FrameMaker parameters FMP1
        //                                copied into FMP

  {
  hsize   = FMP1.hsize;			// Copy horizontal plot size (cm)
  vsize   = FMP1.vsize;			// Copy vertical plot size   (cm)
  haxmin  = FMP1.haxmin;		// Copy horizontal axis label start 
  haxmax  = FMP1.haxmax;		// Copy horizontal axis label finish
  vaxmin  = FMP1.vaxmin;		// Copy vertical axis label start
  vaxmax  = FMP1.vaxmax;		// Copy vertical axis label finish
  VRES    = FMP1.VRES;			// Copy vertical plot resolution
  HRES    = FMP1.HRES;			// Copy horizontal plot resolution
  vaxmax  = FMP1.vaxmax;		// Copy vertical axis label finish
  duse    = FMP1.duse;			// Copy type of data to plot
  PLmax   = FMP1.PLmax;			// Copy polyline size
  FMdebug = FMP1.FMdebug;		// Copy debugging flag
  }

 
// ____________________________________________________________________________
//                        PARAMETER ACCESS FUNCTIONS
// ____________________________________________________________________________
 
// ------------------------- Overall Plot Dimensions --------------------------

double FMPar::HPlotSize() const { return hsize; }

        // Input                FMP     : FrameMaker parameters (this)
        // Output               hsize   : Horizontal plot size (cm)


void FMPar::HPlotSize(double hps) { hsize = hps; }

        // Input                FMP     : FrameMaker parameters (this)
	//			hps	: Plotting dimension (cm)
        // Output               void	: Horizontal plot size set to hps



double FMPar::VPlotSize() const { return vsize; }

        // Input                FMP     : FrameMaker parameters (this)
        // Output               vsize   : Vertical plot size (cm)


void FMPar::VPlotSize(double vps) { vsize = vps; }

        // Input                FMP     : FrameMaker parameters (this)
	//			vps	: Plotting dimension (cm)
        // Output               void	: Vertical plot size set to vps


// ------------------------------ Axis Labeling -------------------------------


double FMPar::HAxMin() const { return haxmin; }

        // Input                FMP     : FrameMaker parameters (this)
        // Output               haxmin	: Horizontal axis left label


 
void FMPar::HAxMin(double hmval) { haxmin = hmval; }
 
        // Input                FMP     : FrameMaker parameters (this)
        //                      hmval   : Horizontal axis left value
        // Output               void    : Horizontal axis left value set

 

double FMPar::HAxMax() const { return haxmax; }

        // Input                FMP     : FrameMaker parameters (this)
        // Output               haxmax  : Horizontal axis right label



void FMPar::HAxMax(double hmval) { haxmax = hmval; }
 
        // Input                FMP     : FrameMaker parameters (this)
        //                      hmval   : Horizontal axis right value
        // Output               void    : Horizontal axis right value set


// ------------------------------ PolyLine Size -------------------------------


int FMPar::PLsize() const

        // Input                FMP     : FrameMaker parameters (this)
        // Output               PLmax   : Maximum allowed PolyLine size

  { return PLmax; }


void FMPar::PLsize(int PLS)

        // Input                FMP     : FrameMaker parameters (this)
	//			PLS	: Polyline size
        // Output               void    : Maximum allowed PolyLine size
	//				  is set in FMP

  { PLmax = PLS; }


// ----------------------------- Plot Resolution ------------------------------


double FMPar::PVRes() const

        // Input                FMP     : FrameMaker parameters (this)
        // Output               VRES 	: Vertical plot resolution (cm)

  { return VRES; }


void FMPar::PVRes(double vcm)

        // Input                FMP     : FrameMaker parameters (this)
	//			vcm	: Vertical plot resolution (cm)
        // Output               void    : Vertical plot resolution is
	//				  is set to vcm in FMP

  { VRES = vcm; }


double FMPar::PHRes() const

        // Input                FMP     : FrameMaker parameters (this)
        // Output               HRES 	: Horizontal plot resolution (cm)

  { return HRES; }


void FMPar::PHRes(double hcm)

        // Input                FMP     : FrameMaker parameters (this)
	//			hcm	: Horizontal plot resolution (cm)
        // Output               void    : Horizontal plot resolution is
	//				  is set to hcm in FMP

  { HRES = hcm; }




// ------------------------------- Data Usage ---------------------------------


datause FMPar::DataUse() const { return duse; }

        // Input                FMP     : FrameMaker parameters (this)
        // Output               duse	: The type of data used is returned


void FMPar::DataUse(datause du) { duse = du; }

        // Input                FMP     : FrameMaker parameters (this)
	//			du	: Type of data useage
        // Output               void	: Type of data to use is set

 

// -------------------------------- Debugging ---------------------------------


int FMPar::DebugLev() const

        // Input                FMP     : FrameMaker parameters (this)
        // Output               FMdebug	: The debugging level is returned

  { return FMdebug; }


void FMPar::DebugLev(int dl)

        // Input                FMP     : FrameMaker parameters (this)
	//			dl	: Debugging level
        // Output               void	: The debugging level is set to dl

  { FMdebug = dl; }
 

// ____________________________________________________________________________
//                      PARAMETER AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 

void FMPar::maxima(const row_vector& vx, double &min, double& max)

        // Input                FMP     : FrameMaker parameters (this)
	// 			vx	: Data vector
	//			min	: Minimum of vx
	//			max	: Maximum of vx
	// Return		void	: Values of min & max are altered

  {
  int np = vx.elements();		// Get the number of points
  max = -HUGE_VAL;			// Set values of max & min so
  min = HUGE_VAL;			// that we can find the maxima
  double tmp;
  int i;
  switch(duse)
    {
    case plotreals:
    default:
      for(i=0; i<np; i++)		// Search all real points
        {				// for maximum & minimum
        tmp = vx.getRe(i);
        if (tmp<min) min = tmp;
        if (tmp>max) max = tmp;
        }
      break;
    case plotimags:
      for(i=0; i<np; i++)		// Search all imaginary points
        {				// for maximum & minimum
        tmp = vx.getIm(i);
        if (tmp<min) min = tmp;
        if (tmp>max) max = tmp;
        }
      break;
    case plotnorms:
      for(i=0; i<np; i++)		// Search all point norms
        {				// for maximum & minimum
        tmp = norm(vx.get(i));
        if (tmp<min) min = tmp;
        if (tmp>max) max = tmp;
        }
      break;
    }
  }


void FMPar::borders(double &top, double &bottom, double &right, double &left)

        // Input	FMP     : FrameMaker parameters (this)
	//		top	: Actual plot top in cm
	//		bottom	: Actual plot bottom in cm
        //		right	: Actual plot right in cm
        //		left	: Actual plot left in cm
	// Return       void	: Values of top, bottom, right, & left set

  {
  top = 0.5;		// Top & bottom borders (relative to Frame!)
  bottom = vsize - 1.0;	// 0.5 & 1 cm margins at top & bottom
  right = hsize - 0.5;	// Left & right borders (relative to Frame!)
  left = 2.0;		// 2 & 0.5 cm margins at left and right
  }


// ____________________________________________________________________________
//                         PARAMETER CHECKING FUNCTIONS
// ____________________________________________________________________________


void FMPar::plotsize()

        // Input                FMP     : FrameMaker parameters (this)
        // Output               void    : The overall plotsize is checked
        // Note				: If unreasonable plot sizes have
	//				  been set they are adjusted

  {
  if((vsize<3) || (vsize>FMPAGEHEIGHT))
     vsize = FMPAGEHEIGHT;
  if((hsize<3) || (hsize>FMPAGEWIDTH))
    hsize = FMPAGEWIDTH;
  }


// ____________________________________________________________________________
//                           STANDARD I/O FUNCTIONS
// ____________________________________________________________________________


ostream& FMPar::print(ostream& ostr) const

        // Input                FMP     : FrameMaker parameters (this)
        //                      ostr	: An output stream
        // Output               ostr	: The output stream modified by
        //                                the FrameMaker parameters

  {
  ostr << "\n\tFrameMaker Plot Output Controls\n";
  ostr << "\n\t\tHorizontal Plot Size:       " << hsize << " cm";
  ostr << "\n\t\tVertical Plot Size:         " << vsize << " cm";
  ostr << "\n\t\tHorizontal Axis Start:      " << haxmin;
  ostr << "\n\t\tHorizontal Axis Finish:     " << haxmax;
  ostr << "\n\t\tVertical Axis Start:        " << vaxmin;
  ostr << "\n\t\tVertical Axis Finish:       " << vaxmax;
  ostr << "\n\t\tVertical Plot Resolution:   " << VRES << " cm";
  ostr << "\n\t\tHorizontal Plot Resolution: " << HRES << " cm";
  ostr << "\n\t\tData Useage In Plot:        " << int(duse);
  ostr << "\n\t\tMaximum PolyLine Size:      " << PLmax << " Points";
  ostr << "\n\t\tDebugging Level:            " << FMdebug;
  ostr << "\n";
  return ostr;
  }

 
ostream& operator<< (ostream& ostr, const FMPar& FMP)
 
        // Input                ostr	: An output stream
        // 			FMP	: Spin system (this)
        // Output               ostr	: The output stream modified by
        //                                the FrameMaker parameters
     
  {
  FMP.print(ostr);				// Call the overload
  return ostr;					// Return the stream
  }

#endif 					// FrameMakerP.cc
