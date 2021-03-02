/* FrameMap.cc **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Reference Frame Mapping  		Implementation		**
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
** A variable of type FrameMap defines a transformation between a	**
** vector of coordinate axes in one frame into a vector of coordinate   **
** axes in another frame. Each coordinate axes {Xi,Yi,Zi} is assigned	**
** a set of Euler angles {alphai, betai, gammai} that will rotate them	**
** into another set of axes {xi,yi,zi}. The initial axes are designated **
** as {xi,yi,zi} and are in the initial frame. The final axes are 	**
** {Xi, Yi, Zi) and in the final frame. There can be any number of axes **
** sets mapped in the class and each may have a different set of Euler	**
** angles that rotate them between reference frames.			**
**                                                                      **
*************************************************************************/

#ifndef   FrameMap_cc_			// Is file already included?
#  define FrameMap_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <IntRank2/FrameMap.h>		// Include interface definition
#include <Basics/Gutils.h>		// Include GAMMA std errors
#include <Basics/Gconstants.h>		// Include GAMMA constants
#include <Level2/EAngles.h>		// Include Euler angles
#include <Basics/StringCut.h>		// Include form & dec functions
#if defined(_MSC_VER) || defined(__SUNPRO_CC)				// If we are using MSVC++
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ STL lists
using std::vector;			// Using libstdc++ STL vectors
using std::ofstream;			// Using libstdc++ output file streams
using std::ostream;			// Using libstdc++ output streams


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                CLASS REFERENCE FRAME MAP ERROR HANDLING
// ____________________________________________________________________________

/*       Input                FrmMap  : Frame Mapping (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void FrameMap::FMerror(int eidx, int noret) const
  {
  string hdr("Frame Mapping");
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
    case 50: GAMMAerror(hdr,"Cannot Set Frame Names",      noret);break; //(50)
    case 51: GAMMAerror(hdr,"Cannot Set Number Axes Sets", noret);break; //(51)
    case 52: GAMMAerror(hdr,"Cannot Set All Euler Angles", noret);break; //(52)
    }
  }

volatile void FrameMap::FMfatal(int eidx) const
  {
  FMerror(eidx, 1);				// Output error message
  if(eidx) FMerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }


/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Funciton pname
      5     Please Use Class Member Funciton pname
   default  Unknown Error - pname                                           */

void FrameMap::FMerror(int eidx, const string& pname, int noret) const 
  {
  string hdr("Frame Mapping");
  string msg;
  switch(eidx)
    {
    case  21: msg = string("Can't Access Object ")     + pname; 	// (21)
              GAMMAerror(hdr, msg, noret); break;
    case 101: msg = string("Can't Find Parameters For ") + pname;	//(101)
              GAMMAerror(hdr, msg, noret); break;
    case 102: msg = string("Can't Find Parameter ") + pname;		//(102)
              GAMMAerror(hdr, msg, noret); break;
    default:  GAMMAerror(hdr, eidx, pname, noret);         break;
    }
  }

volatile void FrameMap::FMfatal(int eidx, const string& pname) const
  {
  FMerror(eidx, 1);				// Output error message
  if(eidx) FMfatal(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                REFERENCE FRAME MAP CHECKING FUNCTIONS
// ____________________________________________________________________________

bool FrameMap::ChkIdx(int i, bool warn) const
  {
  if(i<0 || i>=int(size()))		// Check that object exists
    {					// If not we will fail
    if(warn)				//   If warnings desired we output
      {					//   message(s) as desired
      FMerror(20, 1);			//     Rotation not present
      FMerror(21,Gdec(i),1);		//     Cannot access object i
      }
    return false;			//   Rotation i is not present
    }
  return true;				// AOK, rotation is present
  }

// ____________________________________________________________________________
// iii           REFERENCE FRAME MAP PARAMETER SET SETUP FUNCTIONS 
// ____________________________________________________________________________

/* These functions try and set a coordinate frame from parameters found in
   a particular parameter set.  The initial and final axes names are input as

                 FrmName  (2) : AxesName - Initial Coordinate Axes
                 FrmName  (2) : AxesName - Final Coordinate Axes

   and individual obect rotations are defined by the Euler angles that take the
   final frame (associated with this coordinate frame) into the initial frame.

              FrmRot(#)  (3) : ( alpha, beta, gamma) - Framef Into Framei

   where the # is used to indicate the object that uses the rotation.

           Input        FrmMap  : Frame Mapping (this)
                        pset    : A parameter set
                        i       : Initial frame index
                        f       : Final frame index
                        pfx     : Prefix on rotation parameters
                        warn    : Warning level
                                    0 - no warnings
                                    1 - warnings
                                   >1 - fatal warnings
       Output           TF      : Frame Mapping is set
                                  from parameters in pset                    */

bool FrameMap::SetFrmMap(const ParameterSet& pset, int i, int f, int pfx, int warn)
  {
  ParameterSet  subpset;                        // Copy parameter set
  if(pfx != -1) subpset = pset.strip(pfx);      // to glean out those for
  else          subpset = pset;                 // the specified prefix only
  _EAs.clear();					// Remove all Euler angles
  _EA = EAzero;					// Clear generic Euler angles
  if(!SetNames(subpset, i, f, warn?true:false))	// Try to read frame names
    {						//   If unsuccessful then
    if(warn) FMerror(50, 1);			//   warn if desired and
    return false;				//   return failure
    }
  if(!SetNAxes(subpset, i, f, warn?true:false))	// Try to read # axes sets
    {						//   If unsuccessful then
    if(warn) FMerror(51, 1);			//   warn if desired and
    return false;				//   return failure
    }
  
  if(SetEAngles(subpset,i,f,false)) return true;// Done if one EAngles set
  if(SetEAngVec(subpset,i,f,warn?true:false))  return true;// Done if all EAngles read
  if(warn) FMerror(52, 1);			// Warn if desired and
  return false;					// return failure
  }

bool FrameMap::SetNames(const ParameterSet& pset, int i, int f, bool warn)
  {
  string pnb("FrmName(");			// Base frame name
  string pname = pnb + Gdec(i) + ")";		// Initial frame parameter
  string pstate;				// Temp string for parsing
  ParameterSet::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);			// Seek parameter in pset
  if(item == pset.end())                        // If not found, warn if
    {                                           // desired, return fail
    if(warn) FMerror(102,pname,1);
    return false;
    }
  (*item).parse(pname,Framei,pstate);		// If found, parse/set Framei
  pname = pnb + Gdec(f) + ")";			// Final frame parameter
  item = pset.seek(pname);			// Seek parameter in pset
  if(item == pset.end())                        // If not found, warn if
    {                                           // desired, return fail
    if(warn) FMerror(102,pname,1);
    return false;
    }
  (*item).parse(pname,Framef,pstate);		// If found, parse/set Framef
  return true;
  }

bool FrameMap::SetNAxes(const ParameterSet& pset, int i, int f, bool warn)
  {
  string pname = "FrmAxes("			// Reference frame axes
               + Gdec(i) + "," + Gdec(f) + ")";
  string pstate;				// Temp string for parsing
  ParameterSet::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);			// Seek parameter in pset
  if(item == pset.end())                        // If not found, warn if
    {                                           // desired, return fail
    if(warn) FMerror(102,pname,1);
    return false;
    }
  (*item).parse(pname,NAx,pstate);		// If found, parse/set # axes
  return true;
  }

bool FrameMap::SetEAngles(const ParameterSet& pset, int i, int f, bool warn)
  {
  string pname = "FrmRotAll("			// Reference frame axes
               + Gdec(i) + "," + Gdec(f) + ")";
  string pstate;				// Temp string for parsing
  ParameterSet::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);			// Seek parameter in pset
  if(item == pset.end())                        // If not found, warn if
    {                                           // desired, return fail
    if(warn) FMerror(102,pname,1);
    return false;
    }
  coord ea = coord(*item); 			// Set coordinate from pset
  _EA.alpha(ea.x()*DEG2RAD); 			// Set our alpha value
  _EA.beta(ea.y()*DEG2RAD);			// Set our beta value
  _EA.gamma(ea.z()*DEG2RAD);			// Set our gamma value
  return true;
  }

bool FrameMap::SetEAngVec(const ParameterSet& pset, int i, int f, bool warn)
  {
  string pnbase = "FrmRot("			// Reference frame axes
                 + Gdec(i) + "," + Gdec(f) + ",";
  string pstate;				// Temp string for parsing
  string pname;
  ParameterSet::const_iterator item;		// A pix into parameter list
  coord ea;
  EAngles EA;
  for(int k=0; k<NAx; k++)
    {
    pname = pnbase + Gdec(k) + ")"; 
    item = pset.seek(pname);			// Seek parameter in pset
    if(item == pset.end())			// If not found, warn if
      {						// desired, return fail
      if(warn) FMerror(102,pname,1);
      return false;
      }
    coord ea = coord(*item); 			// Set coordinate from pset
    EA.alpha(ea.x()*DEG2RAD); 			// Set our alpha value
    EA.beta( ea.y()*DEG2RAD);			// Set our beta value
    EA.gamma(ea.z()*DEG2RAD);			// Set our gamma value
    _EAs.push_back(EA);				// Save this set
    }
  return true;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                   REFERENCE FRAME MAP CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/* Each coordinate frame contains the name of the initial frame, the name of
   the final frame, and an array of Euler angles defined over (unspecified)
   objects that rotate from the final frame into the initial frame.          */

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

FrameMap::FrameMap()
  {
  Framei = "Initial";				// Set initial frame name
  Framef = "Final";				// Set final frame name
  _EA    = EAzero; 				// Set default Euler angles
  NAx    = 0;					// Set no axes sets
  }

// ----------------------------------------------------------------------------
//                    Constructors Using The Same Rotation
// ----------------------------------------------------------------------------

FrameMap::FrameMap(const string& Fi, const string& Ff, const EAngles& EA, int N)
  {
  Framei = Fi;				// Set initial frame name
  Framef = Ff;				// Set final frame name
  NAx    = N;				// Set number of axes sets
  _EA    = EA;				// Set Euler angles, all axes sets
  }

// ----------------------------------------------------------------------------
//                    Constructors Using Multiple Rotations
// ----------------------------------------------------------------------------

FrameMap::FrameMap(const string& Fi, const string& Ff, 
                                                  const vector<EAngles>& EAvec)
  {
  Framei = Fi;				// Set initial frame name
  Framef = Ff;				// Set final frame name
  NAx    = EAvec.size();		// Set number of axes sets 
  _EAs   = EAvec;			// Set rotations from Af into Ai
  _EA    = EAzero; 			// Set default Euler angles
  }

// ----------------------------------------------------------------------------
//                    Construction Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions will construct the coordinate frame from parameters found
   in a specified GAMMA parameter set. This allows the coordinate frame to be
   generated from a series of individual rotations as specified in an 
   external ASCII file.                                                      */

FrameMap::FrameMap(const ParameterSet& pset, int i, int f, int pfx, int warn)
  { 
  bool TF = SetFrmMap(pset,i,f,pfx,warn?1:0);	// Try to construct
  if(!TF && warn)				// If construction fails &
    { 						// warnings desired 
    FMerror(2, 1);				//   Problems in construction
    if(warn>1) FMfatal(22,Gdec(pfx));		//   Cant construct rotation
    else       FMerror(22,Gdec(pfx),1);
    }
  }

// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

FrameMap& FrameMap::operator=(const FrameMap& FrmMap) 
  {
  if(this == &FrmMap) return *this;	// Nothing if self-assignment
  Framei = FrmMap.Framei;			// Copy name of initial frame
  Framef = FrmMap.Framef;			// Copy name of final frame
  _EAs   = FrmMap._EAs;			// Copy vector of Euler angles
  _EA    = FrmMap._EA;			// Copy single Euler angles set
  NAx    = FrmMap.NAx;			// Copy number of axes sets
  return *this;				// We are not same as FrmMap
  }

FrameMap::~FrameMap () { }

// ____________________________________________________________________________
// B                     REFERENCE FRAME MAP ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions allow users to either get or set individual rotations in
   the vector of rotations. They also allow for obtaining or setting of the
   coordinate axes names.                                                    */

string  FrameMap::InitialFrame() const           { return Framei; }
string  FrameMap::FinalFrame()   const           { return Framef; }

void    FrameMap::InitialFrame(const string& Ai) { Framei = Ai; }
void    FrameMap::FinalFrame(const   string& Af) { Framef = Af; }
 
EAngles FrameMap::EA(int i) const              { return i<0?_EA:_EAs[i]; }
void    FrameMap::EA(const EAngles& EA, int i) 
  {
  if(i<0)           { _EA  = EA; _EAs.clear(); return; }
  if(_EA != EAzero) 
    { 
    _EAs = vector<EAngles>(NAx, _EA);
    _EAs[i] = EA;
    _EA = EAzero;
    }
  }

int FrameMap::size() const { return NAx; }

// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//      Functions To Make A Parameter Set From A Reference Frame Mapping
// ----------------------------------------------------------------------------

/* These functions will construct a parameter set from a reference frame map.
   Individual rotations will be placed into the parameter set as Euler 
   angle in degrees.                                                         */

FrameMap::operator ParameterSet( ) const
   { ParameterSet pset; pset += *this; return pset; }	

void operator+= (ParameterSet& pset, const FrameMap& FrmMap)
  { FrmMap.PSetAdd(pset, 0, 1); }

void FrameMap::PSetAdd(ParameterSet& pset, int i, int f, int pfx) const
  {
  int nr = size();				// Number of objects
  if(!nr) return;				// Nothing if no objects
  string prefx;                                 // Parameter prefix
  if(pfx != -1)                                 // Only use prefix if pfx
    prefx = string("[")+Gdec(pfx)+string("]");	// is NOT -1
  string    pname = prefx + "FrmName";		// Parameter name
  string    pvalue = Framei;			// Parameter value
  string    pstate("Reference Frame Name");	// Parameter statement
  SinglePar par(pname,pvalue,pstate);		// Parameter for initial axes
  pset.push_back(par);				// Add parameter to set
  pname  = prefx + "FrmName";			// Parameter name
  pvalue = Framef;				// Parameter value
  par = SinglePar(pname,pvalue,pstate);		// Parameter for final state
  pset.push_back(par);				// Add parameter to set
  pname  = prefx + "FrmAxes("			// Parameter name
         + Gdec(i) + "," + Gdec(f) + ")";
  pstate = "Numer of Mapped Axes Sets";		// Parameter statement
  par = SinglePar(pname,NAx,pstate);		// Parameter for # axes sets
  pset.push_back(par);				// Add parameter to set
  string sfx = "(" + Gdec(i) 			// Angle suffix name
             + "," + Gdec(f) + ")";
  
  if(_EA != EAzero)
    {
    pstate = Framef+" Into "+Framei+" (Deg)";	// Parameter statement
    pname  = prefx + "FrmRotAll" + sfx;		// Parameter name
    par    = _EA.param(pname, pstate);		// Parameter for Euler angles
    pset.push_back(par);			// Add parameter to set
    return;
    }

  pstate = Framef+" Into "+Framei+" (Degrees)";	// Parameter statement
  string pbi = "FrmRot(" + Gdec(i) 		// Angle suffix name
             + "," + Gdec(f) + ",";
  for(int l=0; l<nr; l++)			// Loop over individual
    {						// objects
    pname = prefx + pbi + Gdec(l) + ")";	//   Parameter name ith object
    par   = _EAs[l].param(pname,pstate,true);	//   Add parameter to set
    pset.push_back(par);			//   Add parameter to set
    } 
  }

// ----------------------------------------------------------------------------
//      Functions To Make A ASCII Parameter File From A Frame Mapping
// ----------------------------------------------------------------------------

/* These functions will construct an ASCII parameter file from a coordinate 
   frame.  Individual rotations will be placed into the file parameters as
   Euler angle in degrees.                                                   */

void FrameMap::write(const string& filename, int i, int f) const
   {
   ofstream ofstr(filename.c_str());	// Open filename for input
   if(!ofstr.good())                    // If file bad then exit
     {
     FMerror(1, filename);		// Filename problems
     FMfatal(40);			// Cannot write to file
     }
   ofstr.close();			// Close it now
   ParameterSet pset;			// Declare a parameter set
   PSetAdd(pset,i,f);			// Add ourself into parameter set
   pset.write(filename);		// Write parameter set to filename
   return;
   }

// ____________________________________________________________________________
// D                      REFERENCE FRAME MAP INPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions are used to specify the coordinate frame from either a 
   working GAMMA parameter set or an external ASCII file (in parameter set
   format).  The functions do NOT allow for the suffix (#), since it is 
   used/reserved for specific object rotation indices. The functions do allow
   for the prefix [#] so that multiple coordinate frames may be defined in the
   same file by switching the prefix indices.

           Input                FrmMap    : Frame Mapping
                                file    : Input file name
                                pset    : Parameter set
                                argv    : Vector of argc arguments
                                argn    : Argument index
                                pfx     : Coordinate frame index (default -1->none)
           Output               none    : Coordinate frame is read in
                                          from parameters in file filename
                                file    : Name of file used to set frame    */

bool FrameMap::read(const string& filename, int i, int j, int pfx, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename, warn?1:0))	// Read in pset from file
    {
    FMerror(1, filename,1);		// Filename problems
    if(warn > 1) FMfatal(41);		// Cant read from file, fatal 
    else         FMerror(41);		// or non-fatal warning
    return false;
    }
  if(!read(pset,i,j,pfx,warn?1:0))	// User overloaded function
    {
    FMerror(1, filename,1);		// Filename problems
    if(warn > 1) FMfatal(42);		// Cannot read file parameters
    else         FMerror(42);		// or non-fatal one
    return false;
    }
  return true;
  }

bool FrameMap::read(const ParameterSet& pset, int i, int f, int pfx, int warn)
  { 
  bool TF = SetFrmMap(pset,i,f,pfx,warn?1:0);	// Try & set ourself up
  if(!TF)                                       // If we didn't handle
    {                                           // setting ourself properly
    if(warn)                                    // then we'll issue some
      {                                         // warnings if desired
                   FMerror(33, 1);		//   Cant read parameters
      if(warn > 1) FMfatal(34);			//   Can't set from parameters
      else         FMerror(34,1);		//   or a warning issued
      }
    return false; 				// Return that we failed
    }
  return TF;
  }

string FrameMap::ask_read(int argc, char* argv[], int argn,
                                                         int i, int f, int pfx)
  {
  string filename;				// Name of parameter file
  query_parameter(argc, argv, argn,		// Get filename from command
  "\n\tFrame Mapping Filename? ", filename);	// Or ask for it
  read(filename, i, f, pfx); 			// Read rotation from filename
  return filename;
  }

string FrameMap::ask_read(int argc, char* argv[], int argn, 
                                      const string& def, int i, int f, int pfx)
  {
  string msg = "\n\tFrame Mapping Filename ["	// Query we will ask if
             + def + "]? ";                     // it is needed
  string filename = def;                        // Name of spin system file
  ask_set(argc,argv,argn,msg,filename);         // Or ask for it
  read(filename, i, f, pfx);			// Read system from filename
  return filename;                              // Return filename
  }

// ____________________________________________________________________________
// E                  REFERENCE FRAME MAP FORMATTED OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions will output information concerning the coordiante frame
   to any specified output stream.

           Input                FrmMap    : Frame Mapping (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Sum rotation only
                                           !0 - Full composite rotation
           Output               none    : Frame Mapping information
                                          placed into the output stream      */


ostream& FrameMap::print(ostream& ostr, int fflag) const
  {
  int nrot = size();
  if(!nrot)
    {
    string hdr("Empty Frame Mapping\n");
    string spacer = string(40-hdr.length()/2, ' ');
    ostr << "\n\n" << spacer << hdr << "\n";
    return ostr;
    }

/*       If Only One Set Of Euler Angles To Relate Frames Output Appears As

                           Frame Mapping: Framef --> Framei

                      Axes                      Euler Angles
                     Count               alpha     beta      gamma
                   --------            --------  --------  --------         
                     _NAx             _EA.alpha  _EA.beta  _EA.gamma         */

  string hdr = "Frame Mapping: " + Framei + " ===> " + Framef;
  string spacer = string((40-hdr.length()/2), ' ');
  ostr << "\n" << spacer << hdr << "\n" ;
  string spc1("   ");
  string line1 = "  Axes  " + spc1 + "       Euler Angles";
  string line2 = " Count  " + spc1 + " alpha     beta     gamma";
  string line3 = "--------" + spc1 + "--------  --------  --------";
  spacer = string((40-line3.length()/2), ' ');
  string spc2("  ");
  if(_EA != EAzero)
    {
    ostr << "\n" << spacer << line1;
    ostr << "\n" << spacer << line2;
    ostr << "\n" << spacer << line3;
    ostr << "\n   " << spacer << Gdec(NAx) << "       ";
    ostr << Gform("%8.4f", _EA.alpha()*RAD2DEG) << spc2
         << Gform("%8.4f", _EA.beta()*RAD2DEG)  << spc2
         << Gform("%8.4f", _EA.gamma()*RAD2DEG);
    return ostr;
    }

/*         If Multiple Axes Sets (Multiple Euler Angles) Output Appears As

                         Frame Mapping: Framef --> Framei

  Axes            Euler Angles          Axes          Euler Angles
 Index     alpha    beta     gamma     Index    alpha     beta    gamma
-------- -------- -------- --------   -------- -------- -------- --------         
   0     _EAs[0].alpha,beta,gamma        1     _EAs[1].alpha,beta,gamma
   2     _EAs[2].alpha,beta,gamma        3     _EAs[3].alpha,beta,gamma
   .               ....                  .               ....
   .               ....                  .               ....                */

  ostr << "\n Axes            Euler Angles               Axes           Euler Angles";
  ostr << "\n Index      alpha      beta     gamma       Index     alpha      beta     gamma";
  ostr << "\n"
         << "--------   --------  --------  --------    --------  --------  --------  --------";
  bool twocols = true; 
  for(int i=0; i<NAx; i++)
    {
    if(twocols == true) { twocols = false; ostr << "\n"; }
    else                { twocols = true; }
    ostr << "   " << Gdec(i) << "       ";
    ostr << Gform("%8.4f", _EAs[i].alpha()*RAD2DEG) << spc2
         << Gform("%8.4f", _EAs[i].beta()*RAD2DEG)  << spc2
         << Gform("%8.4f", _EAs[i].gamma()*RAD2DEG) << "   ";
    }
  ostr << "\n\n";
  return ostr;
  }

ostream& operator<< (ostream& out, const FrameMap& FrmMap) { return FrmMap.print(out); }

#endif							// FrameMap.cc
