/* CoordFrame.h *************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Coordinate Frame 				Interface	**
**                                                                      **
**      Scott Smith                                                     **
**      Copyright (c) 2002                                              **
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
** A variable of type CoordFrame defines a transformation between two	**
** sets	of coordinate axes over a series of components. Often, such a	**
** transformation is defined by a single set of Euler angles, but that	**
** assumes that all components share initial and final axes. This	**
** class makes no such assumption, allowing the rotation which takes	**
** one frame into another to vary. 					**
**                                                                      **
*************************************************************************/

#ifndef _CoordFrame_h_			// Is file already included?
#define _CoordFrame_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma interface 			// This is the interface
#endif

#include <Basics/ParamSet.h>		// Include parameter sets
#include <Level1/coord.h>		// Include coordinates
#include <Level2/EAngles.h>		// Include Euler angles
#include <Level2/Quaternion.h>		// Include quaternions
#include <string>			// Include stdlibc++ strings

class CoordFrame
  {
  std::string Axesi;			// Initial axes 
  std::string Axesf;			// Final axes 
  std::vector<EAngles> EAs;		// Vector of Euler angles

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                  COORDINATE FRAME ERROR HANDLING
// ____________________________________________________________________________

/*       Input 		      CFrm    : Coordinate Frame (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

         void CFRerror(int eidx,                           int noret=0) const;
         void CFRerror(int eidx, const std::string& pname, int noret=0) const;
volatile void CFRfatal(int eidx)                                        const;
volatile void CFRfatal(int eidx, const std::string& pname)              const;
 
// ____________________________________________________________________________
// ii                COORDINATE FRAME CHECKING FUNCTIONS
// ____________________________________________________________________________

bool ChkIdx(int i, bool warn=true) const;	// Check index i OK

// ____________________________________________________________________________
// iii         COORDINATE FRAME PARAMETER SET SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions try and set a coordinate frame from parameters found in
   a particular parameter set.  The initial and final axes names are input as

                 CoordFrmi  (2) : AxesName - Initial Coordinate Axes
                 CoordFrmf  (2) : AxesName - Final Coordinate Axes

   and individual obect rotations are defined by the Euler angles that take the
   final frame (associated with this coordinate frame) into the initial frame.

              FrmRot(#)  (3) : ( alpha, beta, gamma) - Axesf Into Axesi

   where the # is used to indicate the object that uses the rotation. 

           Input        CFrm    : Coordinate Frame (this)
                        pset    : A parameter set
                        pfx     : Prefix on rotation parameters
                        warn    : Warning level
                                    0 - no warnings
                                    1 - warnings
                                   >1 - fatal warnings
       Output           TF      : Coordinate Frame is set
                                  from parameters in pset		     */

bool SetCoordFrm(const    ParameterSet& pset,int pfx=-1,  int warn=2);
bool SetRotation(const    ParameterSet& pset,             int idx=-1);
bool GetEulerAngles(const ParameterSet& pset,EAngles& EA, int idx=-1);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A                 COORDINATE FRAME CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/* Rotation construction. Since each rotation may contain a series of single
   rotations we allow use arrays for their construction. Thus the rotation may
   be constructed empty (no rotations), with a single rotation (either using
   Euler angles or a Quaternion), or with a vector of rotations (again
   either using Euler angles or Quaternions).                                */

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

CoordFrame();

// ----------------------------------------------------------------------------
//                   Constructors Using Single Rotations
// ----------------------------------------------------------------------------

CoordFrame(double alpha, double beta=0, double gamma=0);
CoordFrame(const  EAngles& EA);
CoordFrame(const  quatern& Q);
 
// ----------------------------------------------------------------------------
//                  Constructors Using Multiple Rotations
// ----------------------------------------------------------------------------

CoordFrame(const std::vector<EAngles>& EAvec);
CoordFrame(const std::vector<quatern>& Qvec);

// ----------------------------------------------------------------------------
//                    Construction Using Parameter Sets
// ----------------------------------------------------------------------------

CoordFrame(const ParameterSet& pset, int idx=-1, int warn=2);

// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

CoordFrame& operator= (const CoordFrame& CFrm);
     ~CoordFrame();

// ____________________________________________________________________________
// B                    COORDINATE FRAME ACCESS FUNCTIONS
// ____________________________________________________________________________

EAngles EA(int i) const;
void    EA(const EAngles& ea, int i);

std::string  InitialAxes() const;
std::string  FinalAxes()   const;

std::string  InitialAxes(const std::string& Ai);
std::string  FinalAxes(const   std::string& Af);

int      size() const;

// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//          Functions To Make A Parameter Set From A Coordinate Frame
// ----------------------------------------------------------------------------

/* These functions will construct a parameter set from a coordinate frame.
   Individual rotations will be placed into the parameter set as Euler
   angle in degrees.                                                         */

            operator ParameterSet( ) const;
friend void operator+= (ParameterSet& pset, const CoordFrame& CFrm);
       void PSetAdd(ParameterSet& pset,   int pfx=-1) const;

// ----------------------------------------------------------------------------
//      Functions To Make A ASCII Parameter File From A Coordinate Frame
// ----------------------------------------------------------------------------

/* These functions will construct an ASCII parameter file from a coordinate
   frame.  Individual rotations will be placed into the file parameters as
   Euler angle in degrees.                                                   */

void write(const std::string &filename, int pfx=-1) const;

// ____________________________________________________________________________
// D                      COORDINATE FRAME INPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions are used to specify the coordinate frame from either a
   working GAMMA parameter set or an external ASCII file (in parameter set
   format).  The functions do NOT allow for the suffix (#), since it is
   used/reserved for specific object rotation indices. The functions do allow
   for the prefix [#] so that multiple coordinate frames may be defined in the
   same file by switching the prefix indices.

           Input                CFrm    : Coordinate Frame
                                file    : Input file name
                                pset    : Parameter set
                                argv    : Vector of argc arguments
                                argn    : Argument index
                                pfx     : Coordinate frame index (default -1->none)
           Output               none    : Coordinate frame is read in
                                          from parameters in file filename
                                file    : Name of file used to set frame    */

bool read(const std::string& file,  int idx=-1, int warn=2);
bool read(const ParameterSet& pset, int idx=-1, int warn=2);
std::string ask_read(int argc, char* argv[], int argn,
                                                                   int idx=-1);
std::string ask_read(int argc, char* argv[], int argn,
                                                const std::string& def, int idx=-1);
 
// ____________________________________________________________________________
// E                       COORDINATE FRAME OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions will output information concerning the coordiante frame 
   to any specified output stream.

           Input                CFrm	: Coordinate Frame (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Sum rotation only
                                           !0 - Full composite rotation
           Output               none    : Coordinate Frame information
                                          placed into the output stream      */

       std::ostream& print(std::ostream& out, int fflag=-1) const;
friend std::ostream& operator<<  (std::ostream& out, const CoordFrame& CFrm);

  };

#endif							// CoordFrame.h
