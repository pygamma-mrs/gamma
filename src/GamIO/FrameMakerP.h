/* FrameMakerP.h ************************************************-*-c++-*-
**									**
**	                          G A M M A 				**
**								 	**
**	FrameMaker Parameters 			     Interface		**
**						 			**
**	Copyright (c) 1997						**
**	Scott A. Smith							**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**								 	**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** The Class FMPar is a simple data type to facility the passing of	**
** parameters to GAMMA FrameMaker functions.				**
**								 	**
*************************************************************************/

#ifndef   GFrameMakerP_h_		// Is this file already included?
#  define GFrameMakerP_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/row_vector.h>		// Know about row_vectors
#include <Matrix/matrix.h>		// Know about matrices

extern const double FMPAGEWIDTH;	// Default FM page width cm
extern const double FMPAGEHEIGHT;	// Default FM page height cm
extern const int    FMPLmax;		// Default FM PolyLine size

#define FM_PP        2			// Digits past decimal point (FMMatrix)
#define FM_PT        0.001		// Default threshold         (FMMatrix)


enum datause {
               plotreals,		// Plot real data points
               plotimags,		// Plot imaginary data points
               plotboth,		// Plot both real and imaginary points
               plotnorms,		// Plot point norms
               plotall3			// Plot reals, imaginaries, & norms
             };
   

class FMPar
  {
  double hsize;				// FM horizontal plot size (cm)
  double vsize;				// FM vertical plot size   (cm)
  double haxmin;			// Horizontal axis label start
  double haxmax;			// Horizontal axis label finish
  double vaxmin;			// Vertical axis label start
  double vaxmax;			// Vertical axis label finish
  double VRES;				// Vertical plot resolution (cm)
  double HRES;				// Horizontal plot resolution (cm)
  datause duse;				// Type of data to plot
  int PLmax;				// Maximum polyline size
  int FMdebug;				// Debugging level

friend class FMStack;			// Allow FMStack full access

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
//                       CLASS FMPARAMETERS ERROR HANDLING
// ____________________________________________________________________________

         
void FMPerror(int eidx, int noret=0) const;
 
        // Input                FMP     : FrameMaker Parameters (this)
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message


volatile void FMPfatality(int eidx) const;
 
        // Input                FMP     : FrameMaker Parameters (this)
        //                      eidx    : Flag for error type
        // Output               none    : Error message output
        //                                Program execution stopped

 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:                                                                         
 
// ____________________________________________________________________________
// A                   FMPARAMETERS CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

MSVCDLC FMPar();
 
        // Input		none	: No arguments required
        // Output               none	: FrameMaker parameters
        ///F_list FMP			- Constructor
 
 
MSVCDLC ~FMPar();
 
        // Input                FMP	: FrameMaker parameters (this)
        // Output               none	: FMP is destructed
 
 
MSVCDLL void operator=(const FMPar &FMP1);
 
        // Input                FMP	: FrameMaker parameters (this)
        // 			FMP1	: FrameMaker parameters
        // Output               void 	: FrameMaker parameters FMP1
        //                                copied into FMP
        ///F_list =			- Assignment
 
 
// ____________________________________________________________________________
//                        PARAMETER ACCESS FUNCTIONS
// ____________________________________________________________________________


// ------------------------- Overall Plot Dimensions --------------------------


MSVCDLL double HPlotSize() const;
 
        // Input                FMP     : FrameMaker parameters (this)
        // Output               hsize   : Horizontal plot size (cm)
 
 
MSVCDLL void HPlotSize(double hps);

        // Input                FMP     : FrameMaker parameters (this)
        //                      hps     : Plotting dimension (cm)
        // Output               void    : Horizontal plot size set to hps
 

MSVCDLL double VPlotSize() const;

        // Input                FMP     : FrameMaker parameters (this)
        // Output               vsize   : Vertical plot size (cm)

 
MSVCDLL void VPlotSize(double vps);

        // Input                FMP     : FrameMaker parameters (this)
        //                      vps     : Plotting dimension (cm)
        // Output               void    : Vertical plot size set to vps
 

// ------------------------------ Axis Labeling -------------------------------


MSVCDLL double HAxMin() const;

        // Input                FMP     : FrameMaker parameters (this)
        // Output               haxmin  : Horizontal axis left label
 

MSVCDLL void HAxMin(double hmval);

        // Input                FMP     : FrameMaker parameters (this)
        // 			hmval   : Horizontal axis left value
        // Output               void 	: Horizontal axis left value set
 
 
MSVCDLL double HAxMax() const;
 
        // Input                FMP     : FrameMaker parameters (this)
        // Output               haxmax  : Horizontal axis right label
 

MSVCDLL void HAxMax(double hmval);

        // Input                FMP     : FrameMaker parameters (this)
        // 			hmval   : Horizontal axis right value
        // Output               void 	: Horizontal axis right value set
 
 
// ------------------------------ PolyLine Size -------------------------------
 
 
MSVCDLL int PLsize() const;
 
        // Input                FMP     : FrameMaker parameters (this)
        // Output               PLmax   : Maximum allowed PolyLine size
 
 
MSVCDLL void PLsize(int PLS);

        // Input                FMP     : FrameMaker parameters (this)
        //                      PLS     : Polyline size
        // Output               void    : Maximum allowed PolyLine size
        //                                is set in FMP 


// ----------------------------- Plot Resolution ------------------------------
 

MSVCDLL double PVRes() const;
 
        // Input                FMP     : FrameMaker parameters (this)
        // Output               VRES    : Vertical plot resolution (cm)
 
 
MSVCDLL void PVRes(double vcm);
 
        // Input                FMP     : FrameMaker parameters (this)
        //                      vcm     : Vertical plot resolution (cm)
        // Output               void    : Vertical plot resolution is
        //                                is set to vcm in FMP
 
 
MSVCDLL double PHRes() const;
 
        // Input                FMP     : FrameMaker parameters (this)
        // Output               HRES    : Horizontal plot resolution (cm)
 
 
MSVCDLL void PHRes(double hcm);
 
        // Input                FMP     : FrameMaker parameters (this)
        //                      hcm     : Horizontal plot resolution (cm)
        // Output               void    : Horizontal plot resolution is
        //                                is set to hcm in FMP

                   
// ------------------------------- Data Usage ---------------------------------
 
 
MSVCDLL datause DataUse() const;

        // Input                FMP     : FrameMaker parameters (this)
        // Output               duse    : The type of data used is returned
 
 
 
MSVCDLL void DataUse(datause du);
 
        // Input                FMP     : FrameMaker parameters (this)
        //                      du      : Type of data useage
        // Output               void    : Type of data to use is set
 
 
// -------------------------------- Debugging ---------------------------------

 
MSVCDLL int DebugLev() const;
 
        // Input                FMP     : FrameMaker parameters (this)
        // Output               FMdebug : The debugging level is returned
 
 
MSVCDLL void DebugLev(int dl);
 
        // Input                FMP     : FrameMaker parameters (this)
        //                      dl      : Debugging level
        // Output               void    : The debugging level is set to dl


// ____________________________________________________________________________
//                      PARAMETER AUXILIARY FUNCTIONS
// ____________________________________________________________________________

                                                                                
MSVCDLL void maxima(const row_vector& vx, double &min, double& max);
   
        // Input                FMP     : FrameMaker parameters (this)
        //                      vx      : Data vector
        //                      min     : Minimum of vx
        //                      max     : Maximum of vx
        // Return               void    : Values of min & max are altered
 

MSVCDLL void borders(double &top, double &bottom, double &right, double &left);

        // Input        FMP     : FrameMaker parameters (this)
        //              top     : Actual plot top in cm
        //              bottom  : Actual plot bottom in cm
        //              right   : Actual plot right in cm
        //              left    : Actual plot left in cm
        // Return       void    : Values of top, bottom, right, & left set

// ____________________________________________________________________________
//                      PARAMETER CHECKING FUNCTIONS
// ____________________________________________________________________________
 
 
MSVCDLL void plotsize();
 
        // Input                FMP     : FrameMaker parameters (this)
        // Output               void    : The overall plotsize is checked
        // Note                         : If unreasonable plot sizes have
        //                                been set they are adjusted
 

// ____________________________________________________________________________
// I                     FMPARAMETERS STANDARD I/O FUNCTIONS
// ____________________________________________________________________________


MSVCDLL std::ostream& print(std::ostream& ostr) const;

        // Input                FMP	: FrameMaker parameters (this)
        // 			ostr	: An output stream
        // Output               none	: Modifies output stream
        ///F_list print			- FrameMaker parameters sent to ostr


MSVCDLL friend std::ostream& operator<<(std::ostream& ostr, const FMPar& FMP);

        // Input                ostr	: An output stream
	//			FMP	: FrameMaker parameters
        // Output               none	: Modifies output stream
        ///F_list << 			- FrameMaker parameters sent to ostr
 
  };


#endif 							// FrameMakerP.h
