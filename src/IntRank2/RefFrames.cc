/* RefFrames.cc *************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Reference Frames 			Implementation		**
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
** A variable of type RefFrame defines a transformation between a       **
** vector of coordinate axes in one frame into a vector of coordinate   **
** axes in another frame. A variable of type RefFrames is a vector of   **
** multiple RefFrame objects.                                           **
**                                                                      **
*************************************************************************/

#ifndef   RefFrames_cc_			// Is file already included?
#  define RefFrames_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <IntRank2/RefFrames.h>		// Include interface definition
#include <Basics/Gutils.h>		// Include GAMMA std errors
#include <Basics/Gconstants.h>		// Include GAMMA constants
#include <Level2/EAngles.h>		// Include Euler angles
#include <Basics/StringCut.h>		// Include form & dec functions
#ifdef _MSC_VER				// If we are using MSVC++
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                CLASS REFERENCE FRAMES ERROR HANDLING
// ____________________________________________________________________________

/*       Input                CFrm    : Reference Frame (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void RefFrames::RFSerror(int eidx, int noret) const
  {
  string hdr("Reference Frames");
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

volatile void RefFrames::RFSfatal(int eidx) const
  {
  RFSerror(eidx, 1);				// Output error message
  if(eidx) RFSerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }


/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Funciton pname
      5     Please Use Class Member Funciton pname
   default  Unknown Error - pname                                           */

void RefFrames::RFSerror(int eidx, const string& pname, int noret) const 
  {
  string hdr("Reference Frames");
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

volatile void RefFrames::RFSfatal(int eidx, const string& pname) const
  {
  RFSerror(eidx, 1);				// Output error message
  if(eidx) RFSfatal(0);				// State that its fatal
  exit(-1);					// Abort the program
  }

// ____________________________________________________________________________
// ii                REFERENCE FRAMES CHECKING FUNCTIONS
// ____________________________________________________________________________

bool RefFrames::ChkIdx(int i, bool warn) const
  {
  if(i<0 || i>=int(size()))		// Check that object exists
    {					// If not we will fail
    if(warn)				//   If warnings desired we output
      {					//   message(s) as desired
      RFSerror(20, 1);			//     Rotation not present
      RFSerror(21,Gdec(i),1);		//     Cannot access object i
      }
    return false;			//   Rotation i is not present
    }
  return true;				// AOK, rotation is present
  }

// ____________________________________________________________________________
// iii           REFERENCE FRAMES PARAMETER SET SETUP FUNCTIONS 
// ____________________________________________________________________________

/* These functions try and set a coordinate frame from parameters found in
   a particular parameter set.  The initial and final axes names are input as

                 RefFramesi  (2) : AxesName - Initial Coordinate Axes
                 RefFramesf  (2) : AxesName - Final Coordinate Axes

   and individual obect rotations are defined by the Euler angles that take the
   final frame (associated with this coordinate frame) into the initial frame.

              FrmRot(#)  (3) : ( alpha, beta, gamma) - Framef Into Framei

   where the # is used to indicate the object that uses the rotation.

           Input        CFrm    : Reference Frame (this)
                        pset    : A parameter set
                        pfx     : Prefix on rotation parameters
                        warn    : Warning level
                                    0 - no warnings
                                    1 - warnings
                                   >1 - fatal warnings
       Output           TF      : Reference Frame is set
                                  from parameters in pset                    */

bool RefFrames::SetRefFrm(const ParameterSet& pset, int pfx, int warn)
  {
  ParameterSet  subpset;                        // Copy parameter set
  if(pfx != -1) subpset = pset.strip(pfx);      // to glean out those for
  else          subpset = pset;                 // the specified prefix only
  _EAs.clear();					// Remove all Euler angles
  _EA = EAzero;					// Clear generic Euler angles
  if(!SetNames(subpset, warn))			// Try to read frame names
    {						//   If unsuccessful then
    if(warn) RFSerror(50, 1);			//   warn if desired and
    return false;				//   return failure
    }
  if(!SetNAxes(subpset, warn))			// Try to read # axes sets
    {						//   If unsuccessful then
    if(warn) RFSerror(51, 1);			//   warn if desired and
    return false;				//   return failure
    }
  
  if(SetEAngles(subpset, false)) return true;	// Done if one EAngles set
  if(SetEAngVec(subpset, warn))  return true;	// Done if all EAngles read
  if(warn) RFSerror(52, 1);			// Warn if desired and
  return false;					// return failure
  }

bool RefFrames::SetNames(const ParameterSet& pset, bool warn)
  {
  string pname("RefFramesi");			// Initial frame name
  string pstate;				// Temp string for parsing
  list<SinglePar>::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);			// Seek parameter in pset
  if(item == pset.end())                        // If not found, warn if
    {                                           // desired, return fail
    if(warn) RFSerror(102,pname,1);
    return false;
    }
  (*item).parse(pname,Framei,pstate);		// If found, parse/set Framei
  pname = string("RefFramesf");			// Final frame name
  item = pset.seek(pname);			// Seek parameter in pset
  if(item == pset.end())                        // If not found, warn if
    {                                           // desired, return fail
    if(warn) RFSerror(102,pname,1);
    return false;
    }
  (*item).parse(pname,Framef,pstate);		// If found, parse/set Framei
  return true;
  }

bool RefFrames::SetNAxes(const ParameterSet& pset, bool warn)
  {
  string pname("FrmAxes");			// Reference frame axes
  string pstate;				// Temp string for parsing
  list<SinglePar>::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);			// Seek parameter in pset
  if(item == pset.end())                        // If not found, warn if
    {                                           // desired, return fail
    if(warn) RFSerror(102,pname,1);
    return false;
    }
  (*item).parse(pname,NAx,pstate);		// If found, parse/set # axes
  return true;
  }

bool RefFrames::SetEAngles(const ParameterSet& pset, bool warn)
  {
  string pname("FrmRotAll");			// Reference frame axes
  string pstate;				// Temp string for parsing
  list<SinglePar>::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);			// Seek parameter in pset
  if(item == pset.end())                        // If not found, warn if
    {                                           // desired, return fail
    if(warn) RFSerror(102,pname,1);
    return false;
    }
  coord ea = coord(*item); 			// Set coordinate from pset
  _EA.alpha(ea.x()*DEG2RAD); 			// Set our alpha value
  _EA.beta(ea.y()*DEG2RAD);			// Set our beta value
  _EA.gamma(ea.z()*DEG2RAD);			// Set our gamma value
  return true;
  }


bool RefFrames::SetEAngVec(const ParameterSet& pset, bool warn)
  {
  string pnbase("FrmRot(");			// Reference frame axes
  string pstate;				// Temp string for parsing
  string pname;
  list<SinglePar>::const_iterator item;		// A pix into parameter list
  coord ea;
  EAngles EA;
  for(int i=0; i<NAx; i++)
    {
    pname = pnbase + Gdec(i) + ")"; 
    item = pset.seek(pname);			// Seek parameter in pset
    if(item == pset.end())			// If not found, warn if
      {						// desired, return fail
      if(warn) RFSerror(102,pname,1);
      return false;
      }
    coord ea = coord(*item); 			// Set coordinate from pset
    EA.alpha(ea.x()*DEG2RAD); 			// Set our alpha value
    EA.beta(ea.y()*DEG2RAD);			// Set our beta value
    EA.gamma(ea.z()*DEG2RAD);			// Set our gamma value
    _EAs.push_back(EA);				// Save this set
    }
  return true;
  }



// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                   REFERENCE FRAMES CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/* Each coordinate frame contains the name of the initial frame, the name of
   the final frame, and an array of Euler angles defined over (unspecified)
   objects that rotate from the final frame into the initial frame.          */

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

RefFrames::RefFrames() {};

// ----------------------------------------------------------------------------
//                    Constructors Using The Same Rotation
// ----------------------------------------------------------------------------

RefFrames::RefFrames(const string& Fi, const string& Ff, 
                                                      const EAngles& EA, int N)
  {
  Framei = Fi;				// Set initial frame name
  Framef = Ff;				// Set final frame name
  NAx    = N;				// Set number of axes sets
  _EA    = EA;				// Set Euler angles, all axes sets
  }

// ----------------------------------------------------------------------------
//                    Construction Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions will construct the vector of reference frame mappings from
   parameters found in a specified GAMMA parameter set. This allows the 
   vector of frame mappings to be generated from a series of individual 
   frame maps as specified in an external ASCII file.                        */

RefFrames::RefFrames(const ParameterSet& pset, int warn)
  { 
  bool TF = SetRefFrm(pset, pfx, warn?1:0);	// Try to construct
  if(!TF && warn)				// If construction fails &
    { 						// warnings desired 
    RFSerror(2, 1);				//   Problems in construction
    if(warn>1) RFSfatal(22,Gdec(pfx));		//   Cant construct rotation
    else       RFSerror(22,Gdec(pfx),1);
    }
  }

// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

RefFrames& RefFrames::operator=(const RefFrames& CFrm) 
  {
  if(this == &CFrm) return *this;	// Nothing if self-assignment
  Framei = CFrm.Framei;			// Copy name of initial frame
  Framef = CFrm.Framef;			// Copy name of final frame
  _EAs   = CFrm._EAs;			// Copy vector of Euler angles
  _EA    = CFrm._EA;			// Copy single Euler angles set
  NAx    = CFrm.NAx;			// Copy number of axes sets
  return *this;				// We are not same as CFrm
  }

RefFrames::~RefFrames () { }

// ____________________________________________________________________________
// B                     REFERENCE FRAMES ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions allow users to either get or set individual rotations in
   the vector of rotations. They also allow for obtaining or setting of the
   coordinate axes names.                                                    */

string  RefFrames::InitialAxes() const           { return Framei; }
string  RefFrames::FinalAxes()   const           { return Framef; }

void    RefFrames::InitialAxes(const string& Ai) { Framei = Ai; }
void    RefFrames::FinalAxes(const   string& Af) { Framef = Af; }
 
EAngles RefFrames::EA(int i) const              { return i<0?_EA:_EAs[i]; }
void    RefFrames::EA(const EAngles& EA, int i) 
  {
  if(i<0)           { _EA  = EA; _EAs.clear(); return; }
  if(_EA != EAzero) 
    { 
    _EAs = vector<EAngles>(NAx, _EA);
    _EAs[i] = EA;
    _EA = EAzero;
    }
  }

int RefFrames::size() const { return NAx; }

// ____________________________________________________________________________
// C                   PARAMETER & PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//          Functions To Make A Parameter Set From A Reference Frame
// ----------------------------------------------------------------------------

/* These functions will construct a parameter set from a coordinate frame.
   Individual rotations will be placed into the parameter set as Euler 
   angle in degrees.                                                         */

RefFrames::operator ParameterSet( ) const
   { ParameterSet pset; pset += *this; return pset; }	

void operator+= (ParameterSet& pset, const RefFrames& CFrm)
  { CFrm.PSetAdd(pset); }

void RefFrames::PSetAdd(ParameterSet& pset, int pfx) const
  {
  int nr = size();				// Number of objects
  if(!nr) return;				// Nothing if no objects
  string prefx;                                 // Parameter prefix
  if(pfx != -1)                                 // Only use prefix if pfx
    prefx = string("[")+Gdec(pfx)+string("]");	// is NOT -1
  string    pname = prefx + "RefFramesi";	// Parameter name
  string    pvalue = Framei;			// Parameter value
  string    pstate("Inital Reference Frame");	// Parameter statement
  SinglePar par(pname,pvalue,pstate);		// Parameter for initial axes
  pset.push_back(par);				// Add parameter to set
  pname  = prefx + "RefFramesf";			// Parameter name
  pvalue = Framef;				// Parameter value
  pstate = "Final Refernce Frame";		// Parameter statement
  par = SinglePar(pname,pvalue,pstate);		// Parameter for final state
  pset.push_back(par);				// Add parameter to set
  pname  = prefx + "FrmAxes";			// Parameter name
  pstate = "Numer of Frame Axes Sets";		// Parameter statement
  par = SinglePar(pname,NAx,pstate);		// Parameter for # axes sets
  pset.push_back(par);				// Add parameter to set
  
  if(_EA != EAzero)
    {
    pstate = Framef+" Into "+Framei+" (Deg)";	// Parameter statement
    pname  = prefx + "FrmRotAll";		// Parameter name
    par    = _EA.param(pname, pstate);		// Parameter for Euler angles
    pset.push_back(par);			// Add parameter to set
    return;
    }

  pstate = Framef+" Into "+Framei+" (Degrees)";	// Parameter statement
  string pbi("FrmRot(");			// Parameter base name
  coord ABG;					// Workign coordinate

  for(int i=0; i<nr; i++)			// Loop over individual
    {						// objects
    pname = prefx + pbi + Gdec(i) + ")";	//   Parameter name ith object
    ABG.x(_EAs[i].alpha()*RAD2DEG);		//   We use a coordinate for
    ABG.y(_EAs[i].beta()*RAD2DEG);		//   output and set the angles
    ABG.z(_EAs[i].gamma()*RAD2DEG);		//   to be in degrees
    par = ABG.param(pname, pstate);		//   Parameter for Euler angles
    pset.push_back(par);			//   Add parameter to set
    } 
  }

// ----------------------------------------------------------------------------
//      Functions To Make A ASCII Parameter File From A Reference Frame
// ----------------------------------------------------------------------------

/* These functions will construct an ASCII parameter file from a coordinate 
   frame.  Individual rotations will be placed into the file parameters as
   Euler angle in degrees.                                                   */

void RefFrames::write(const string& filename, int pfx) const
   {
   ofstream ofstr(filename.c_str());	// Open filename for input
   if(!ofstr.good())                    // If file bad then exit
     {
     RFSerror(1, filename);		// Filename problems
     RFSfatal(40);			// Cannot write to file
     }
   ofstr.close();			// Close it now
   ParameterSet pset;			// Declare a parameter set
   PSetAdd(pset,pfx);			// Add ourself into parameter set
   pset.write(filename);		// Write parameter set to filename
   return;
   }

// ____________________________________________________________________________
// D                      REFERENCE FRAMES INPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions are used to specify the coordinate frame from either a 
   working GAMMA parameter set or an external ASCII file (in parameter set
   format).  The functions do NOT allow for the suffix (#), since it is 
   used/reserved for specific object rotation indices. The functions do allow
   for the prefix [#] so that multiple coordinate frames may be defined in the
   same file by switching the prefix indices.

           Input                CFrm    : Reference Frame
                                file    : Input file name
                                pset    : Parameter set
                                argv    : Vector of argc arguments
                                argn    : Argument index
                                pfx     : Coordinate frame index (default -1->none)
           Output               none    : Coordinate frame is read in
                                          from parameters in file filename
                                file    : Name of file used to set frame    */

bool RefFrames::read(const string& filename, int pfx, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename, warn?1:0))	// Read in pset from file
    {
    RFSerror(1, filename,1);		// Filename problems
    if(warn > 1) RFSfatal(41);		// Cant read from file, fatal 
    else         RFSerror(41);		// or non-fatal warning
    return false;
    }
  if(!read(pset,pfx,warn?1:0))		// User overloaded function
    {
    RFSerror(1, filename,1);		// Filename problems
    if(warn > 1) RFSfatal(42);		// Cannot read file parameters
    else         RFSerror(42);		// or non-fatal one
    return false;
    }
  return true;
  }

bool RefFrames::read(const ParameterSet& pset, int pfx, int warn)
  { 
  bool TF = SetRefFrm(pset,pfx,warn?1:0);	// Try & set ourself up
  if(!TF)                                       // If we didn't handle
    {                                           // setting ourself properly
    if(warn)                                    // then we'll issue some
      {                                         // warnings if desired
                   RFSerror(33, 1);		//   Cant read parameters
      if(warn > 1) RFSfatal(34);		//   Can't set from parameters
      else         RFSerror(34,1);		//   or a warning issued
      }
    return false; 				// Return that we failed
    }
  return TF;
  }

string RefFrames::ask_read(int argc, char* argv[], int argn, int idx)
  {
  string filename;				// Name of parameter file
  query_parameter(argc, argv, argn,		// Get filename from command
  "\n\tReference Frame Filename? ", filename); // Or ask for it
  read(filename, idx); 				// Read rotation from filename
  return filename;
  }

string RefFrames::ask_read(int argc, char* argv[], int argn, const string& def, int idx)
  {
  string msg = "\n\tReference Frame Filename ["// Query we will ask if
             + def + "]? ";                     // it is needed
  string filename = def;                        // Name of spin system file
  ask_set(argc,argv,argn,msg,filename);         // Or ask for it
  read(filename, idx);                          // Read system from filename
  return filename;                              // Return filename
  }

// ____________________________________________________________________________
// E                  REFERENCE FRAMES FORMATTED OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions will output information concerning the coordiante frame
   to any specified output stream.

           Input                CFrm    : Reference Frame (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Sum rotation only
                                           !0 - Full composite rotation
           Output               none    : Reference Frame information
                                          placed into the output stream      */


ostream& RefFrames::print(ostream& ostr, int fflag) const
  {
  int nrot = size();
  if(!nrot)
    {
    string hdr("Empty Reference Frame\n");
    string spacer = string(40-hdr.length()/2, ' ');
    ostr << "\n\n" << spacer << hdr << "\n";
    return ostr;
    }

/*       If Only One Set Of Euler Angles To Relate Frames Output Appears As

                         Reference Frame: Framef --> Framei

                      Axes                      Euler Angles
                     Count               alpha     beta      gamma
                   --------            --------  --------  --------         
                     _NAx             _EA.alpha  _EA.beta  _EA.gamma         */

  string hdr = "Reference Frame: " + Framei + " ===> " + Framef;
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

                        Reference Frame: Framef --> Framei

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

ostream& operator<< (ostream& out, const RefFrames& CFrm) { return CFrm.print(out); }

#endif							// RefFrames.cc
