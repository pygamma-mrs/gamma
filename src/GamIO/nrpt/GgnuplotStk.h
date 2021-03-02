/* GgnuplotStk.h ************************************************-*-c++-*-
**                                                                      **
**                                 G A M M A                            **
**                                                                      **
**      Gnuplot 2D/3D Stack Plots                   Interface		**
**                                                                      **
**      Scott A. Smith							**
**      Copyright (c) 2002						**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  This file contains the functions that generate 2D/3D stack plots	**
**  in Gnuplot from GAMMA vectors and matrices.				**
**                                                                      **
*************************************************************************/

#ifndef _gnuplotStk_h_			// Is this file already included?
#define _gnuplotStk_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma interface			// This is the interface
#endif

#include <GamIO/GgnuplotC.h>		// Know GAMMA Gnuplot controls
#include <Level1/coord_vec.h>		// Know GAMMA coordinate vectors
#include <string>                       // Know libstdc++ strings
#include <vector>                       // Know libstdc++ STL vectors

class GPStack: public GPControls	// Gnuplot info
  {
  bool   degrees;			// Flag angles in degrees/radians
  double POVtheta;			// Point of view theta angle
  double POVphi;			// Point of view phi angle
  bool   basestack;
  bool   stackaxes;
  double stackscale;
  double zaxisscale;
  bool   spherical;
  bool   hidden;
  bool   dataout;
  bool   normalize;

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i           GAMMA Gnuplot Stackical 2D/3D Plot Error Handling
// ____________________________________________________________________________

/*         Input                eidx    : Error index
                                pname   : string in message
                                noret   : Flag for return (0=linefeed)
           Output               none    : Error message output
                                          Execution stopped (if fatal)       */

         void GPStack::GPSerror(int eidx,                            int noret=0);
         void GPStack::GPSerror(int eidx,  const std::string& pname, int noret=0);
volatile void GPStack::GPSfatal(int eidx);
volatile void GPStack::GPSfatal(int eidx,  const std::string& pname);

//____________________________________________________________________________
// ii                        Gnuplot Structure Functions
//____________________________________________________________________________

void GPStack::defaults();
void GPStack::copy(const GPStack& GPS);

// ____________________________________________________________________________
// iii             Gnuplot Stack 2D/3D Plot Loadfile Output
// ____________________________________________________________________________

void GPStack::WriteTitle(std::ofstream&   ofstr);
void GPStack::WriteView(std::ofstream&    ofstr);
void GPStack::WriteAngles(std::ofstream&  ofstr);
void GPStack::WriteMapping(std::ofstream& ofstr);
void GPStack::WriteSplot(std::ofstream&   ofstr);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A     Gnuplot Stack 2D/3D Controls Contruction, Assigment, Destruction
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

GPStack::GPStack();				// Default controls
GPStack::GPStack(const GPStack& GPS);	// Duplicate controls

virtual   GPStack::~GPStack();
GPStack& GPStack::operator= (const GPStack& GPS);

// ____________________________________________________________________________
// B              Gnuplot Stack 2D/3D Control Access
// ____________________________________________________________________________

void GPStack::SetDegrees(bool       dg);
void GPStack::SetStackAxes(bool    ax);
void GPStack::SetBaseStack(bool    bs);
void GPStack::SetNormalization(bool nm);

// ____________________________________________________________________________
// C              Gnuplot Stack 2D/3D Loadfile Generation
// ____________________________________________________________________________

void GPStack::LoadFile();

// ____________________________________________________________________________
// D              Gnuplot Stack 2D/3D Datafile Generation
// ____________________________________________________________________________

/* The ASCII data is set into three columns { theta, phi, R }. However, Gnuplot
   takes theta as the angle over from the +x axis - that which GAMMA calls the
   phase angle phi. Gnuplot takes the angle up/down from the XY plane as phi,
   which also contrasts with GAMMA in that we take theta as the angle down from
   the +z axis. So, we need to output {PHI, 90-THETA, R} to keep Gnuplot happy.
   Lastly, we can output out data in either degrees or radians, the later
   being the default.                                                        */

bool GPStack::DataAxisFile(char axis, bool warn=true);
void GPStack::DataFile(const coord_vec& data, int idx = -1);
void GPStack::DataFiles(const std::vector<coord_vec>& data);

// ____________________________________________________________________________
// E             Gnuplot Stack 2D/3D Plot Output Generation
// ____________________________________________________________________________

/* These functions will generate the Gnuplot load file and one or more ASCII
   data files that will be plotted. Subsequently, the command is issued to
   run gnuplot using the produced load file which reads the data files.      */

void GPStack::Plot();
void GPStack::Plot(const coord_vec& data);
void GPStack::Plot(const std::vector<coord_vec>& data);

// ____________________________________________________________________________
// Z         Gnuplot Stack 2D/3D Controls Standard Output Functions
// ____________________________________________________________________________

virtual std::ostream& GPStack::print(std::ostream& ostr);
friend  std::ostream& operator<< (std::ostream&     ostr, GPStack& GPS);

// ____________________________________________________________________________
// G                          Non-Member Functions
// ____________________________________________________________________________

/* These are functions built prior to construction of this class. Some of
   them are deprecated and will issue warnings if utilized.                  */

friend void GP_stack(std::ofstream& ofstr, const row_vector &vx, double row,
              int ri=0, double xmin=0, double xmax=0, double cutoff=0);

friend void GP_stack(const std::string& filename, matrix& mx, int ri=0,
                 double ymin=0, double ymax=0, double xmin=0, double xmax=0);
 
friend void GP_stack(std::ofstream& ofstr, matrix& mx, int ri=0,
                 double ymin=0, double ymax=0, double xmin=0, double xmax=0);
 
//------------------------------ Common Stack Plots ---------------------------

 
friend void GP_stackplot(const std::string& gnumacro, const std::string& file2D);
friend void GP_stackblurb(std::ofstream& ostr, const std::string& plotname);

};


#endif 								// GgnuplotStk.h
