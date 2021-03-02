/* FrameMap.h ***************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Reference Frame Mapping 			Interface	**
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
** A variable of type FrameMap defines a transformation between a       **
** vector of coordinate axes in one frame into a vector of coordinate   **
** axes in another frame. Each coordinate axes {Xi,Yi,Zi} is assigned   **
** a set of Euler angles {alphai, betai, gammai} that will rotate them  **
** into another set of axes {xi,yi,zi}. The initial axes are designated **
** as {xi,yi,zi} and are in the initial frame. The final axes are       **
** {Xi, Yi, Zi) and in the final frame. There can be any number of axes **
** sets mapped in the class and each may have a different set of Euler  **
** angles that rotate them between reference frames.                    **
**                                                                      **
*************************************************************************/

#ifndef   FrameMap_h_			// Is file already included?
#  define FrameMap_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Level1/coord.h>		// Include coordinates
#include <Level2/EAngles.h>		// Include Euler angles
#include <string>			// Include stdlibc++ strings

class FrameMap
  {
  std::string          Framei;		// Initial frame 
  std::string          Framef;		// Final frame 
  std::vector<EAngles> _EAs;		// Vector of Euler angles
  EAngles              _EA;		// Single set of Euler angles
  int                  NAx;		// Number of axes sets

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                  REFERENCE FRAME MAP ERROR HANDLING
// ____________________________________________________________________________

/*       Input 		      FrmMap  : Reference Frame mapping (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

         void FMerror(int eidx,                           int noret=0) const;
         void FMerror(int eidx, const std::string& pname, int noret=0) const;
volatile void FMfatal(int eidx)                                        const;
volatile void FMfatal(int eidx, const std::string& pname)              const;
 
// ____________________________________________________________________________
// ii                REFERENCE FRAME MAP CHECKING FUNCTIONS
// ____________________________________________________________________________

bool ChkIdx(int i, bool warn=true) const;	// Check index i OK

// ____________________________________________________________________________
// iii         REFERENCE FRAME MAP PARAMETER SET SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions try and set a coordinate frame from parameters found in
   a particular parameter set.  The initial and final axes names are input as

                 FrameMapi  (2) : AxesName - Initial Reference Axes
                 FrameMapf  (2) : AxesName - Final Reference Axes

   and individual obect rotations are defined by the Euler angles that take the
   final frame (associated with this coordinate frame) into the initial frame.

              FrmRot(#)  (3) : ( alpha, beta, gamma) - Framef Into Framei

   where the # is used to indicate the object that uses the rotation. 

           Input        FrmMap    : Reference Frame (this)
                        pset    : A parameter set
                        pfx     : Prefix on rotation parameters
                        warn    : Warning level
                                    0 - no warnings
                                    1 - warnings
                                   >1 - fatal warnings
       Output           TF      : Reference Frame is set
                                  from parameters in pset		     */

bool SetFrmMap(const  ParameterSet& p,int i, int f,
                                                   int pfx=-1, int  warn=2);
bool SetNames(const   ParameterSet& p, int i, int f, bool warn=true);
bool SetNAxes(const   ParameterSet& p, int i, int f, bool warn=true);
bool SetEAngles(const ParameterSet& p, int i, int f, bool warn=true);
bool SetEAngVec(const ParameterSet& p, int i, int f, bool warn=true);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A                 REFERENCE FRAME MAP CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/* Rotation construction. Since each rotation may contain a series of single
   rotations we allow use arrays for their construction. Thus the rotation may
   be constructed empty (no rotations), with a single rotation (either using
   Euler angles or a Quaternion), or with a vector of rotations (again
   either using Euler angles or Quaternions).                                */

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

MSVCDLC FrameMap();

// ----------------------------------------------------------------------------
//                    Constructors Using The Same Rotation
// ----------------------------------------------------------------------------

MSVCDLC FrameMap(const std::string& Fi, const std::string& Ff,
                                                     const EAngles& EA, int N);

// ----------------------------------------------------------------------------
//                    Constructors Using Multiple Rotations
// ----------------------------------------------------------------------------

MSVCDLC FrameMap(const std::string& Fi, const std::string& Ff,
                                            const std::vector<EAngles>& EAvec);

// ----------------------------------------------------------------------------
//                    Construction Using Parameter Sets
// ----------------------------------------------------------------------------

MSVCDLC FrameMap(const ParameterSet& pset,int i,int f,int pfx=-1,int warn=2);

// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

MSVCDLL FrameMap& operator= (const FrameMap& FrmMap);
MSVCDLC      ~FrameMap();

// ____________________________________________________________________________
// B                    REFERENCE FRAME MAP ACCESS FUNCTIONS
// ____________________________________________________________________________

MSVCDLL std::string  InitialFrame() const;
MSVCDLL std::string  FinalFrame()   const;

MSVCDLL void         InitialFrame(const std::string& Ai);
MSVCDLL void         FinalFrame(const   std::string& Af);

MSVCDLL EAngles      EA(int i=-1) const;
MSVCDLL void         EA(const EAngles& ea, int i);

MSVCDLL int           size() const;

// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//          Functions To Make A Parameter Set From A Reference Frame
// ----------------------------------------------------------------------------

/* These functions will construct a parameter set from a coordinate frame.
   Individual rotations will be placed into the parameter set as Euler
   angle in degrees.                                                         */

MSVCDLL             operator ParameterSet( ) const;
MSVCDLL friend void operator+= (ParameterSet& pset, const FrameMap& FrmMap);
MSVCDLL        void PSetAdd(ParameterSet& pset,int i,int f,int pfx=-1) const;

// ----------------------------------------------------------------------------
//      Functions To Make A ASCII Parameter File From A Reference Frame
// ----------------------------------------------------------------------------

/* These functions will construct an ASCII parameter file from a coordinate
   frame.  Individual rotations will be placed into the file parameters as
   Euler angle in degrees.                                                   */

MSVCDLL void write(const std::string &filename, int i, int f) const;

// ____________________________________________________________________________
// D                      REFERENCE FRAME MAP INPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions are used to specify the coordinate frame from either a
   working GAMMA parameter set or an external ASCII file (in parameter set
   format).  The functions do NOT allow for the suffix (#), since it is
   used/reserved for specific object rotation indices. The functions do allow
   for the prefix [#] so that multiple coordinate frames may be defined in the
   same file by switching the prefix indices.

           Input                FrmMap  : Reference Frame
                                file    : Input file name
                                pset    : Parameter set
                                argv    : Vector of argc arguments
                                argn    : Argument index
                                i       : Initial frame index (0)
                                f       : Final frame index (1)
                                pfx     : Reference frame prefix (default -1->none)
           Output               none    : Reference frame is read in
                                          from parameters in file filename
                                file    : Name of file used to set frame    */

MSVCDLL bool read(const std::string&  file,int i=0,int f=1,int pfx=-1,int w=2);
MSVCDLL bool read(const ParameterSet& pset,int i=0,int f=1,int pfx=-1,int w=2);
MSVCDLL std::string ask_read(int argc, char* argv[], int argn,
                                                 int i=0, int f=1, int idx=-1);
MSVCDLL std::string ask_read(int argc, char* argv[], int argn,
                         const std::string& def, int i=0, int f=1, int idx=-1);
 
// ____________________________________________________________________________
// E                       REFERENCE FRAME MAP OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions will output information concerning the coordiante frame 
   to any specified output stream.

           Input                FrmMap	: Reference Frame (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Sum rotation only
                                           !0 - Full composite rotation
           Output               none    : Reference Frame information
                                          placed into the output stream      */

MSVCDLL        std::ostream& print(std::ostream& out, int fflag=-1) const;
MSVCDLL friend std::ostream& operator<<  (std::ostream& out, const FrameMap& FrmMap);

  };

#endif							// FrameMap.h
