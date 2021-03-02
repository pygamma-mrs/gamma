/* CoordFrame.cc ************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Coordinate Frame  			Implementation		**
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
** A variable of type CoordFrame defines a transformation between two   **
** sets of coordinate axes over a series of components. Often, such a   **
** transformation is defined by a single set of Euler angles, but that  **
** assumes that all components share initial and final axes. This       **
** class makes no such assumption, allowing the rotation which takes    **
** one frame into another to vary.                                      **
**                                                                      **
*************************************************************************/

#ifndef   CoordFrame_cc_		// Is file already included?
#  define CoordFrame_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <IntRank2/CoordFrame.h>	// Include interface definition
#include <Basics/Gutils.h>		// Include GAMMA std errors
#include <Basics/Gconstants.h>		// Include GAMMA constants
#include <Level2/EAngles.h>		// Include Euler angles
#include <Basics/StringCut.h>		// Include form & dec functions
#ifdef _MSC_VER				// If we are using MSVC++
 #pragma warning (disable : 4786)       //   Kill STL namelength warnings
#endif


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                CLASS COORDINATE FRAME ERROR HANDLING
// ____________________________________________________________________________

/*       Input                CFrm    : Coordinate Frame (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void CoordFrame::CFRerror(int eidx, int noret) const
  {
  std::string hdr("Coordinate Frame");
  switch(eidx)
    {
    case  0: GAMMAerror(hdr,"Program Aborting.....",       noret);break; // (0)
    case  2: GAMMAerror(hdr,"Problems During Construction",noret);break; // (2)
    case 11: GAMMAerror(hdr,"Cannot Set Any Rotations",    noret);break; //(11)
    case 12: GAMMAerror(hdr,"Insufficient Parameters",     noret);break; //(12)
    case 20: GAMMAerror(hdr,"Accessed Object Is Absent",   noret);break; //(20)
    case 33: GAMMAerror(hdr,"Cannot Read Parameters",      noret);break; //(33)
    case 34: GAMMAerror(hdr,"Cannot Set From Parameters",  noret);break; //(34)
    case 40: GAMMAerror(hdr,"Cannot Write To File",        noret);break; //(40)
    case 41: GAMMAerror(hdr,"Cannot Read From File",       noret);break; //(41)
    case 42: GAMMAerror(hdr,"Cannot Read File Parameters", noret);break; //(42)
    }
  }

volatile void CoordFrame::CFRfatal(int eidx) const
  {
  CFRerror(eidx, 1);				// Output error message
  if(eidx) CFRerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }


/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Funciton pname
      5     Please Use Class Member Funciton pname
   default  Unknown Error - pname                                           */

void CoordFrame::CFRerror(int eidx, const std::string& pname, int noret) const 
  {
  std::string hdr("Coordinate Frame");
  std::string msg;
  switch(eidx)
    {
    case  21: msg = std::string("Can't Access Object ")     + pname; 	// (21)
              GAMMAerror(hdr, msg, noret); break;
    case 101: msg = std::string("Can't Find Parameters For ") + pname;	//(101)
              GAMMAerror(hdr, msg, noret); break;
    default:  GAMMAerror(hdr, eidx, pname, noret);         break;
    }
  }

volatile void CoordFrame::CFRfatal(int eidx, const std::string& pname) const
  {
  CFRerror(eidx, 1);				// Output error message
  if(eidx) CFRfatal(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                COORDINATE FRAME CHECKING FUNCTIONS
// ____________________________________________________________________________

bool CoordFrame::ChkIdx(int i, bool warn) const
  {
  if(i<0 || i>=int(size()))		// Check that object exists
    {					// If not we will fail
    if(warn)				//   If warnings desired we output
      {					//   message(s) as desired
      CFRerror(20, 1);			//     Rotation not present
      CFRerror(21,Gdec(i),1);		//     Cannot access object i
      }
    return false;			//   Rotation i is not present
    }
  return true;				// AOK, rotation is present
  }

// ____________________________________________________________________________
// iii           COORDINATE FRAME PARAMETER SET SETUP FUNCTIONS 
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
                                  from parameters in pset                    */

bool CoordFrame::SetCoordFrm(const ParameterSet& pset, int pfx, int warn)
  {
  ParameterSet  subpset;                        // Copy parameter set
  if(pfx != -1) subpset = pset.strip(pfx);      // to glean out those for
  else          subpset = pset;                 // the specified prefix only
  EAs.clear();					// Remove all Euler angles
  int i=0;					// Temporary rotation index
  while(SetRotation(subpset, i)) { i++; }	// Read successive rotations
  if(!size())					// If we found some, set the
    {						// as desired, perhaps stop
    CFRerror(11, 1);				//   Cant set any rotations
    if(warn>1) CFRfatal(12);			//   Insufficient parameters
    else       CFRerror(12, 1);
    }
  return false;
  }


bool CoordFrame::SetAxes(const ParameterSet& pset, bool warn)
  {
  std::string pname("CoordFrmi");			// Initial axes name
  std::string pname("CoordFrmi");			// Initial axes name
  }

bool CoordFrame::SetRotation(const ParameterSet& pset, int idx)
  {
  EAngles EA;
  if(GetEulerAngles(pset, EA, idx))	// Try & read Euler angle set
    {
    EAs.push_back(EA);			//   Store Euler angles
    Qs.push_back(quatern(EA));		//   Store Quaternion
    return true;			//   Flag we were succesful
    }
  quatern Q;
  if(GetQuaternion(pset, Q, idx))	// Try & read Quaternion
    {
    EAs.push_back(Q.EA());		//   Store 1st Euler angles
    Qs.push_back(Q);			//   Store 1st Quaternion
    return true;			//   Flag we were succesful
    }
  return false;
  }
  
bool CoordFrame::GetEulerAngles(const ParameterSet& pset, EAngles& EA, int idx)
  { return EA.read(pset, idx); }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                   COORDINATE FRAME CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/* Each coordinate frame contains the name of the initial frame, the name of
   the final frame, and an array of Euler angles defined over (unspecified)
   objects that rotate from the final frame into the initial frame.          */

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

CoordFrame::CoordFrame() { Axesi = "Initial"; Axesf = "Final"; }

// ----------------------------------------------------------------------------
//                    Constructors Using The Same Rotation
// ----------------------------------------------------------------------------

CoordFrame::CoordFrame(const std::string& Ai, const std::string& Af, 
                                                      const EAngles& EA, int N)
  {
  Axesi = Ai;				// Set initial axes name
  Axesf = Af;				// Set final axes name
  EAs = vector<EAngles> (N, EA);	// Set rotations same
  }

// ----------------------------------------------------------------------------
//                    Constructors Using Multiple Rotations
// ----------------------------------------------------------------------------

CoordFrame::CoordFrame(const std::string& Ai, const std::string& Af, 
                                                  const vector<EAngles>& EAvec)
  {
  Axesi = Ai;				// Set initial axes name
  Axesf = Af;				// Set final axes name
  EAs = EAvec;				// Set rotations from Af into Ai
  }

// ----------------------------------------------------------------------------
//                    Construction Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions will construct the coordinate frame from parameters found
   in a specified GAMMA parameter set. This allows the coordinate frame to be
   generated from a series of individual rotations as specified in an 
   external ASCII file.                                                      */

CoordFrame::CoordFrame(const ParameterSet& pset, int pfx, int warn)
  { 
  bool TF = SetCoordFrm(pset, pfx, warn?1:0);	// Try to construct
  if(!TF && warn)				// If construction fails &
    { 						// warnings desired 
    CFRerror(2, 1);				//   Problems in construction
    if(warn>1) CFRfatal(22,Gdec(pfx));		//   Cant construct rotation
    else       CFRerror(22,Gdec(pfx),1);
    }
  }

// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

CoordFrame& CoordFrame::operator=(const CoordFrame& CFrm) 
  {
  if(this == &CFrm) return *this;	// Nothing if self-assignment
  Axesi = CFrm.Axesi;			// Copy name of initial frame
  Axesf = CFrm.Axesf;			// Copy name of final frame
  EAs   = CFrm.EAs;			// Copy vector of Euler angles
  return *this;				// We are not same as CFrm
  }

CoordFrame::~CoordFrame () { }


// ____________________________________________________________________________
// B                     OORDINATE FRAME ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions allow users to either get or set individual rotations in
   the vector of rotations. They also allow for obtaining or setting of the
   coordinate axes names.                                                    */
 
EAngles CoordFrame::EA(int i) const               { return EAs[i]; }
void    CoordFrame::EA(const EAngles& ea, int i)  { EAs[i] = ea;   }

std::string  CoordFrame::InitialAxes() const           { return Axesi; }
std::string  CoordFrame::FinalAxes()   const           { return Axesf; }

std::string  CoordFrame::InitialAxes(const std::string& Ai) { Axesi = Ai; }
std::string  CoordFrame::FinalAxes(const   std::string& Af) { Axesf = Af; }

int     CoordFrame::size()                        { return EAs.size(); }

// ____________________________________________________________________________
// C                        PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//          Functions To Make A Parameter Set From A Coordinate Frame
// ----------------------------------------------------------------------------

/* These functions will construct a parameter set from a coordinate frame.
   Individual rotations will be placed into the parameter set as Euler 
   angle in degrees.                                                         */

CoordFrame::operator ParameterSet( ) const
   { ParameterSet pset; pset += *this; return pset; }	

void operator+= (ParameterSet& pset, const CoordFrame& CFrm)
  { CFrm.PSetAdd(pset); }

void CoordFrame::PSetAdd(ParameterSet& pset, int pfx) const
  {
  int nr = size();				// Number of objects
  if(!nr) return;				// Nothing if no objects
  std::string prefx;                                 // Parameter prefix
  if(pfx != -1)                                 // Only use prefix if pfx
    prefx = std::string("[")+Gdec(pfx)+std::string("]");	// is NOT -1
  std::string    pname = prefx + "CoordFrmi";	// Parameter name
  std::string    pvalue = Axesi;			// Parameter value
  std::string    pstate("Inital Coordinate Axes");	// Parameter statement
  SinglePar par(pname,pvalue,pstate);		// Parameter for initial axes
  pset.push_back(par);				// Add parameter to set
  pname  = prefx + "CoordFrmf";			// Parameter name
  pvalue = Axesf;				// Parameter value
  pstate("Final Coordinate Axes");		// Parameter statement
  par = SinglePar(pname,pvalue,pstate);		// Parameter for final state
  pset.push_back(par);				// Add parameter to set
  
  pstate = Axesf+" Into "+Axesi+" (Degrees)";	// Paraemter statement
  std::string pbi("FrmRot(");			// Parameter base name
  coord ABG;					// Workign coordinate

  for(int i=0; i<nr; i++)			// Loop over individual
    {						// objects
    pname = pst + pb + Gdec(i) + ")";		//   Parameter name ith object
    ABG.x(EAs[i].alpha()*RAD2DEG);		//   We use a coordinate for
    ABG.y(EAs[i].beta()*RAD2DEG);		//   output and set the angles
    ABG.z(EAs[i].gamma()*RAD2DEG);		//   to be in degrees
    par = ABG.param(pname, pstate);		//   Parameter for Euler angles
    pset.push_back(par);			//   Add parameter to set
    } 
  }

// ----------------------------------------------------------------------------
//      Functions To Make A ASCII Parameter File From A Coordinate Frame
// ----------------------------------------------------------------------------

/* These functions will construct an ASCII parameter file from a coordinate 
   frame.  Individual rotations will be placed into the file parameters as
   Euler angle in degrees.                                                   */

void CoordFrame::write(const std::string& filename, int pfx) const
   {
   ofstream ofstr(filename.c_str());	// Open filename for input
   if(!ofstr.good())                    // If file bad then exit
     {
     CFRerror(1, filename);		// Filename problems
     CFRfatal(40);			// Cannot write to file
     }
   ofstr.close();			// Close it now
   ParameterSet pset;			// Declare a parameter set
   PSetAdd(pset,pfx);			// Add ourself into parameter set
   pset.write(filename);		// Write parameter set to filename
   return;
   }

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

bool CoordFrame::read(const std::string& filename, int pfx, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename, warn?1:0))	// Read in pset from file
    {
    CFRerror(1, filename,1);		// Filename problems
    if(warn > 1) CFRfatal(41);		// Cant read from file, fatal 
    else         CFRerror(41);		// or non-fatal warning
    return false;
    }
  if(!read(pset,pfx,warn?1:0))		// User overloaded function
    {
    CFRerror(1, filename,1);		// Filename problems
    if(warn > 1) CFRfatal(42);		// Cannot read file parameters
    else         CFRerror(42);		// or non-fatal one
    return false;
    }
  return true;
  }

bool CoordFrame::read(const ParameterSet& pset, int pfx, int warn)
  { 
  bool TF = SetCoordFrm(pset,pfx,warn?1:0);	// Try & set ourself up
  if(!TF)                                       // If we didn't handle
    {                                           // setting ourself properly
    if(warn)                                    // then we'll issue some
      {                                         // warnings if desired
                   CFRerror(33, 1);		//   Cant read parameters
      if(warn > 1) CFRfatal(34);		//   Can't set from parameters
      else         CFRerror(34,1);		//   or a warning issued
      }
    return false; 				// Return that we failed
    }
  return TF;
  }

std::string CoordFrame::ask_read(int argc, char* argv[], int argn, int idx)
  {
  std::string filename;				// Name of parameter file
  query_parameter(argc, argv, argn,		// Get filename from command
  "\n\tCoordinate Frame Filename? ", filename); // Or ask for it
  read(filename, idx); 				// Read rotation from filename
  return filename;
  }

std::string CoordFrame::ask_read(int argc, char* argv[], int argn, const std::string& def, int idx)
  {
  std::string msg = "\n\tCoordinate Frame Filename ["// Query we will ask if
             + def + "]? ";                     // it is needed
  std::string filename = def;                        // Name of spin system file
  ask_set(argc,argv,argn,msg,filename);         // Or ask for it
  read(filename, idx);                          // Read system from filename
  return filename;                              // Return filename
  }

// ____________________________________________________________________________
// E                  COORDINATE FRAME FORMATTED OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions will output information concerning the coordiante frame
   to any specified output stream.

           Input                CFrm    : Coordinate Frame (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Sum rotation only
                                           !0 - Full composite rotation
           Output               none    : Coordinate Frame information
                                          placed into the output stream      */


std::ostream& CoordFrame::print(std::ostream& ostr, int fflag) const
  {
  int nrot = size();
  if(!nrot)
    {
    std::string hdr("Empty Coordinate Frame\n");
    std::string spacer = std::string(40-hdr.length()/2, ' ');
    ostr << "\n\n" << spacer << hdr << "\n";
    return ostr;
    }

/*          Output Some Printed Lines Which Will Look Like The Following

                         Coordinate Frame: Axesf --> Axesi

                                          [D  , D  , D  ]
   Spin Quantum Numbers:      x.xxxxxx    [ xx   xy   xz]   [ x.x, x.x, x.x]
   DCC (kHz):                 x.xxxxxx    [D  , D  , D  ] = [ x.x, x.x, x.x]
   Xi Value (rad/sec):        x.xxxxxx    [ yx   yy   yz]   [ x.x, x.x, x.x]
                                          [D  , D  , D  ]
                                          [ zx   zy   zz]                    */

  std::string hdr = "Coordinate Frame: " + Axesi + " ===> " + Axesf;
  std::string spacer = std::string((40-hdr.length()/2), ' ');
  ostr << "\n" << spacer << hdr << "\n" ;
  ostr << "\nObject            Euler Angles             Object           Euler Angles"
  ostr << "\n Index      alpha      beta     gamma      Index     alpha      beta     gamma";
  ostr << "\n"
         << "--------   --------  --------  --------   --------  --------  --------  --------";
  for(int i=0; i<nrot; i++)
    {
    ostr << "\n   " << Gdec(i) << "       ";
    ostr << Gform("%8.4f", EAs[i].alpha()) << "  "
    << Gform("%8.4f", EAs[i].beta())  << "  "
    << Gform("%8.4f", EAs[i].gamma()) << "   ";
    }
  ostr << "\n\n";
  return ostr;
  }

std::ostream& operator<< (std::ostream& out, const CoordFrame& CFrm) { return CFrm.print(out); }

#endif							// CoordFrame.cc
