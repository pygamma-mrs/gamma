/* RefFrames.h **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Reference Frames 				Interface	**
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
** A variable of type RefFrame defines a transformation between a 	**
** vector of coordinate axes in one frame into a vector of coordinate 	**
** axes in another frame. A variable of type RefFrames is a vector of	**
** multiple RefFrame objects.						**
**                                                                      **
*************************************************************************/

#ifndef _RefFrames_h_			// Is file already included?
#define _RefFrames_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma interface 			// This is the interface
#endif

#include <Basics/ParamSet.h>		// Include parameter sets
#include <Level1/coord.h>		// Include coordinates
#include <Level2/EAngles.h>		// Include Euler angles
#include <IntRank2/RefFrame>
#include <string>			// Include stdlibc++ strings

class RefFrames
  {
  std::vector<RefFrame> RFs;		// Vector of refrence frames

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                  REFERENCE FRAMES ERROR HANDLING
// ____________________________________________________________________________

/*       Input 		      RFrms    : Reference Frame (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

         void RFSerror(int eidx,                           int noret=0) const;
         void RFSerror(int eidx, const std::string& pname, int noret=0) const;
volatile void RFSfatal(int eidx)                                        const;
volatile void RFSfatal(int eidx, const std::string& pname)              const;
 
// ____________________________________________________________________________
// ii                REFERENCE FRAMES CHECKING FUNCTIONS
// ____________________________________________________________________________

bool ChkIdx(int i, bool warn=true) const;	// Check index i OK

// ____________________________________________________________________________
// iii         REFERENCE FRAMES PARAMETER SET SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions try and set a coordinate frame from parameters found in
   a particular parameter set.  The initial and final axes names are input as

                 RefFramesi  (2) : AxesName - Initial Reference Axes
                 RefFramesf  (2) : AxesName - Final Reference Axes

   and individual obect rotations are defined by the Euler angles that take the
   final frame (associated with this coordinate frame) into the initial frame.

              FrmRot(#)  (3) : ( alpha, beta, gamma) - Framef Into Framei

   where the # is used to indicate the object that uses the rotation. 

           Input        RFrms    : Reference Frame (this)
                        pset    : A parameter set
                        pfx     : Prefix on rotation parameters
                        warn    : Warning level
                                    0 - no warnings
                                    1 - warnings
                                   >1 - fatal warnings
       Output           TF      : Reference Frame is set
                                  from parameters in pset		     */

bool SetRefFrm(const  ParameterSet& pset,int pfx=-1,  int warn=2);
bool SetNames(const   ParameterSet& pset,         bool warn=true);
bool SetNAxes(const   ParameterSet& pset,         bool warn=true);
bool SetEAngles(const ParameterSet& pset,         bool warn=true);
bool SetEAngVec(const ParameterSet& pset,         bool warn=true);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A                 REFERENCE FRAMES CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/* Rotation construction. Since each rotation may contain a series of single
   rotations we allow use arrays for their construction. Thus the rotation may
   be constructed empty (no rotations), with a single rotation (either using
   Euler angles or a Quaternion), or with a vector of rotations (again
   either using Euler angles or Quaternions).                                */

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

RefFrames();

// ----------------------------------------------------------------------------
//                    Constructors Using The Same Rotation
// ----------------------------------------------------------------------------

RefFrames(const std::string& Fi, const std::string& Ff,
                                                     const EAngles& EA, int N);

// ----------------------------------------------------------------------------
//                    Constructors Using Multiple Rotations
// ----------------------------------------------------------------------------

RefFrames(const std::string& Fi, const std::string& Ff,
                                            const std::vector<EAngles>& EAvec);

// ----------------------------------------------------------------------------
//                    Construction Using Parameter Sets
// ----------------------------------------------------------------------------

RefFrames(const ParameterSet& pset, int idx=-1, int warn=2);

// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

RefFrames& operator= (const RefFrames& RFrms);
     ~RefFrames();

// ____________________________________________________________________________
// B                    REFERENCE FRAMES ACCESS FUNCTIONS
// ____________________________________________________________________________

std::string  InitialAxes() const;
std::string  FinalAxes()   const;

void         InitialAxes(const std::string& Ai);
void         FinalAxes(const   std::string& Af);

EAngles      EA(int i=-1) const;
void         EA(const EAngles& ea, int i);

int           size() const;

// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//          Functions To Make A Parameter Set From A Reference Frame
// ----------------------------------------------------------------------------

/* These functions will construct a parameter set from a coordinate frame.
   Individual rotations will be placed into the parameter set as Euler
   angle in degrees.                                                         */

            operator ParameterSet( ) const;
friend void operator+= (ParameterSet& pset, const RefFrames& RFrms);
       void PSetAdd(ParameterSet& pset,   int pfx=-1) const;

// ----------------------------------------------------------------------------
//      Functions To Make A ASCII Parameter File From A Reference Frame
// ----------------------------------------------------------------------------

/* These functions will construct an ASCII parameter file from a coordinate
   frame. Individual rotations will be placed into the file parameters as
   Euler angle in degrees.                                                   */

void write(const std::string &filename, int pfx=-1) const;

// ____________________________________________________________________________
// D                      REFERENCE FRAMES INPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions are used to specify the coordinate frame from either a
   working GAMMA parameter set or an external ASCII file (in parameter set
   format).  The functions do NOT allow for the suffix (#), since it is
   used/reserved for specific object rotation indices. The functions do allow
   for the prefix [#] so that multiple coordinate frames may be defined in the
   same file by switching the prefix indices.

           Input                RFrms    : Reference Frame
                                file    : Input file name
                                pset    : Parameter set
                                argv    : Vector of argc arguments
                                argn    : Argument index
                                pfx     : Reference frame index (default -1->none)
           Output               none    : Reference frame is read in
                                          from parameters in file filename
                                file    : Name of file used to set frame    */

bool read(const std::string& file,  int idx=-1, int warn=2);
bool read(const ParameterSet& pset, int idx=-1, int warn=2);
std::string ask_read(int argc, char* argv[], int argn,
                                                                   int idx=-1);
std::string ask_read(int argc, char* argv[], int argn,
                                                const string& def, int idx=-1);
 
// ____________________________________________________________________________
// E                       REFERENCE FRAMES OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions will output information concerning the coordiante frame 
   to any specified output stream.

           Input                RFrms	: Reference Frame (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Sum rotation only
                                           !0 - Full composite rotation
           Output               none    : Reference Frame information
                                          placed into the output stream      */

       std::ostream& print(std::ostream& out, int fflag=-1) const;
friend std::ostream& operator<<      (std::ostream& out, const RefFrames& RFrms);

  };

#endif							// RefFrames.h
