/* Powder.cc ****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Powder  				Implementation		**
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
** A variable of type Powder defines an array of Euler angles and a     **
** corresponding array of scaling factors. Each set of Euler angles     **
** defines a rotation that takes the desired axes (LAB) into initial    **
** axes (PAS). A corrsponding scaling factor defines the probability    **
** of finding that particular orientation relative to the other         **
** orientations.                                                        **
**                                                                      **
** Although the Euler angles used herein rotate (zum beispiel) LAB axes **
** into PAS axes by definition, GAMMA uses these Euler angles to take   **
** interactions from the PAS into the LAB frame. As is mentioned in     **
** countless places in GAMMA, this "inverse" rotation usage allows one  **
** to relate spherical coordinate angles (theta,phi) to the Euler       **
** angles (beta, alpha).                                                **
**                                                                      **
** Thus, assume that some intial axes are the PAS of some interaction.  **
** Then the Euler angles {alpha,beta,0} would be the set GAMMA uses to  **
** rotate into the laboratory axes so that the PAS z-axis will be at    **
** angle beta down form LAB +z axis and at angle alpha over from LAB    **
** +x axis. This is used for convenience even though {alpha, beta, 0}   **
** actually define the inverse rotation...... AOK? Ganz klar.           **
**                                                                      **
*************************************************************************/

#ifndef _Powder_cc_				// Is file already included?
#define _Powder_cc_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)	// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

//#include <IntRank2/Powder.h>		// Include interface definition
#include "Powder.h"			// Include interface definition
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
// i                CLASS POWDER ERROR HANDLING
// ____________________________________________________________________________

/*       Input                Pow    : Powder (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void Powder::POWerror(int eidx, int noret) const
  {
  string hdr("Powder");
  switch(eidx)
    {
    case  0: GAMMAerror(hdr,"Program Aborting.....",       noret);break; // (0)
    case  2: GAMMAerror(hdr,"Problems During Construction",noret);break; // (2)
    case 11: GAMMAerror(hdr,"Mismatched Dimensions",       noret);break; //(11)
 //   case 11: GAMMAerror(hdr,"Cannot Set Any Rotations",    noret);break; //(11)
 //   case 12: GAMMAerror(hdr,"Insufficient Parameters",     noret);break; //(12)
 //   case 20: GAMMAerror(hdr,"Accessed Object Is Absent",   noret);break; //(20)
 //   case 33: GAMMAerror(hdr,"Cannot Read Parameters",      noret);break; //(33)
 //   case 34: GAMMAerror(hdr,"Cannot Set From Parameters",  noret);break; //(34)
 //   case 40: GAMMAerror(hdr,"Cannot Write To File",        noret);break; //(40)
 //   case 41: GAMMAerror(hdr,"Cannot Read From File",       noret);break; //(41)
 //   case 42: GAMMAerror(hdr,"Cannot Read File Parameters", noret);break; //(42)
    }
  }

volatile void Powder::POWfatal(int eidx) const
  {
  POWerror(eidx, 1);				// Output error message
  if(eidx) POWerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }


/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Funciton pname
      5     Please Use Class Member Funciton pname
   default  Unknown Error - pname                                           */

void Powder::POWerror(int eidx, const string& pname, int noret) const 
  {
  string hdr("Powder");
  string msg;
  switch(eidx)
    {
    case  20: msg = string("Euler Angle Set Count Is ")  + pname; 	// (20)
              GAMMAerror(hdr, msg, noret); break;
    case  21: msg = string("Probabilty  Count Is ")  + pname; 		// (21)
              GAMMAerror(hdr, msg, noret); break;
//    case  21: msg = string("Can't Access Object ")     + pname; 	// (21)
//              GAMMAerror(hdr, msg, noret); break;
    case 101: msg = string("Can't Find Parameters For ") + pname;	//(101)
              GAMMAerror(hdr, msg, noret); break;
    default:  GAMMAerror(hdr, eidx, pname, noret);         break;
    }
  }

volatile void Powder::POWfatal(int eidx, const string& pname) const
  {
  POWerror(eidx, 1);				// Output error message
  if(eidx) POWfatal(0);				// State that its fatal
   GAMMAfatal();				// Clean exit from program
  }

// ____________________________________________________________________________
// ii                       POWDER CHECKING FUNCTIONS
// ____________________________________________________________________________


/*
bool Powder::ChkIdx(int i, bool warn) const
  {
  if(i<0 || i>=int(size()))		// Check that object exists
    {					// If not we will fail
    if(warn)				//   If warnings desired we output
      {					//   message(s) as desired
      POWerror(20, 1);			//     Rotation not present
      POWerror(21,Gdec(i),1);		//     Cannot access object i
      }
    return false;			//   Rotation i is not present
    }
  return true;				// AOK, rotation is present
  }
*/

// ____________________________________________________________________________
// iii           POWDER PARAMETER SET SETUP FUNCTIONS 
// ____________________________________________________________________________

/* These functions try and set a powder from parameters found in
   a particular parameter set.  The initial and final axes names are input as

                 CoordFrmi  (2) : AxesName - Initial Coordinate Axes
                 CoordFrmf  (2) : AxesName - Final Coordinate Axes

   and individual obect rotations are defined by the Euler angles that take the
   final frame (associated with this powder) into the initial frame.

              FrmRot(#)  (3) : ( alpha, beta, gamma) - Axesf Into Axesi

   where the # is used to indicate the object that uses the rotation.

           Input        Pow    : Powder (this)
                        pset    : A parameter set
                        pfx     : Prefix on rotation parameters
                        warn    : Warning level
                                    0 - no warnings
                                    1 - warnings
                                   >1 - fatal warnings
       Output           TF      : Powder is set
                                  from parameters in pset                    */

/*
bool Powder::SetCoordFrm(const ParameterSet& pset, int pfx, int warn)
  {
  ParameterSet  subpset;                        // Copy parameter set
  if(pfx != -1) subpset = pset.strip(pfx);      // to glean out those for
  else          subpset = pset;                 // the specified prefix only
  EAs.clear();					// Remove all Euler angles
  int i=0;					// Temporary rotation index
  while(SetRotation(subpset, i)) { i++; }	// Read successive rotations
  if(!size())					// If we found some, set the
    {						// as desired, perhaps stop
    POWerror(11, 1);				//   Cant set any rotations
    if(warn>1) POWfatal(12);			//   Insufficient parameters
    else       POWerror(12, 1);
    }
  return false;
  }


bool Powder::SetAxes(const ParameterSet& pset, bool warn)
  {
  string pname("CoordFrmi");			// Initial axes name
  string pname("CoordFrmi");			// Initial axes name
  }

bool Powder::SetRotation(const ParameterSet& pset, int idx)
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
  
bool Powder::GetEulerAngles(const ParameterSet& pset, EAngles& EA, int idx)
  { return EA.read(pset, idx); }
*/


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                   POWDER CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/* Each powder contains an array of Euler angles and a cooresponding array
   of probabilities.                                                         */

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

Powder::Powder() { }

// ----------------------------------------------------------------------------
//                       Constructors Using Two Vectors
// ----------------------------------------------------------------------------

Powder::Powder(const vector<EAngles>& EAs, const vector<double>& Ps,
                                                            const string& Name)
  {
  if(_EAs.size() != _Ps.size())		// Insures arrays sized compatible
    {					// If bad, issue warnings & quit
    POWerror(11, 1);			//    Dimension mismatch
    POWerror(20, Gdec(_EAs.size()), 1);	//    Number of Euler angles
    POWerror(21, Gdec(_Ps.size()),  1);	//    Number of probabilities
    POWfatal(2);			//    Error during construction
    }
  _EAs = EAs;				// Copy all Euler angle sets
  _Ps  = Ps;				// Copy all probabilties
  _PN  = Name;				// Set the powder name 
  _A3  = false;				// Assume only two angles
  int NO = _EAs.size();			// Number of orientations
  for(int i=0; i<NO && !_A3; i++)	// Scan input angles and look for
    {					// any having non-zero gamma values
    if(fabs(_EAs[i].gamma()) > 1.e-12)
      _A3 = true;
    }
  }

// ----------------------------------------------------------------------------
//                    Construction Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions will construct the powder from parameters found
   in a specified GAMMA parameter set. This allows the powder to be
   generated from a series of individual rotations as specified in an 
   external ASCII file.                                                      */

/*
Powder::Powder(const ParameterSet& pset, int pfx, int warn)
  { 
  bool TF = SetCoordFrm(pset, pfx, warn?1:0);	// Try to construct
  if(!TF && warn)				// If construction fails &
    { 						// warnings desired 
    POWerror(2, 1);				//   Problems in construction
    if(warn>1) POWfatal(22,Gdec(pfx));		//   Cant construct rotation
    else       POWerror(22,Gdec(pfx),1);
    }
  }
*/

// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

Powder& Powder::operator=(const Powder& Pow) 
  {
  if(this == &Pow) return *this;	// Nothing if self-assignment
  _EAs = Pow._EAs;			// Copy all Euler angle sets
  _Ps  = Pow._Ps;			// Copy all probabilties
  _PN  = Pow._PN;			// Copy the powder name
  _A3  = Pow._A3;			// Copy flag for three angles
  return *this;				// We are not same as Pow
  }

Powder::~Powder () { }


// ____________________________________________________________________________
// B                         POWDER ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions allow users to either get or set individual rotations in
   the vector of rotations. They also allow for obtaining or setting of the
   probabilites.                                                             */
 
EAngles Powder::EA(int i) const               { return _EAs[i]; }
void    Powder::EA(const EAngles& EA, int i)  { _EAs[i] = EA;
                                                if(EA.alpha()) _A3=true; }

double  Powder::Pop(int i) const              { return _Ps[i]; }
void    Powder::Pop(double& P, int i)         { _Ps[i] = P;    }

string  Powder::Name() const                  { return _PN; }
void    Powder::Name(const string& PN)        { _PN = PN;   }

int     Powder::size() const                  { return _EAs.size(); }

// ____________________________________________________________________________
// C                        PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//             Functions To Make A Parameter Set From A Powder
// ----------------------------------------------------------------------------

/* These functions will construct a parameter set from a powder.
   Individual rotations will be placed into the parameter set as Euler 
   angle in degrees.                                                         */

/*
Powder::operator ParameterSet( ) const
   { ParameterSet pset; pset += *this; return pset; }	

void operator+= (ParameterSet& pset, const Powder& Pow)
  { Pow.PSetAdd(pset); }

void Powder::PSetAdd(ParameterSet& pset, int pfx) const
  {
  int nr = size();				// Number of objects
  if(!nr) return;				// Nothing if no objects
  string prefx;                                 // Parameter prefix
  if(pfx != -1)                                 // Only use prefix if pfx
    prefx = string("[")+Gdec(pfx)+string("]");	// is NOT -1
  string    pname = prefx + "CoordFrmi";	// Parameter name
  string    pvalue = Axesi;			// Parameter value
  string    pstate("Inital Coordinate Axes");	// Parameter statement
  SinglePar par(pname,pvalue,pstate);		// Parameter for initial axes
  pset.push_back(par);				// Add parameter to set
  pname  = prefx + "CoordFrmf";			// Parameter name
  pvalue = Axesf;				// Parameter value
  pstate("Final Coordinate Axes");		// Parameter statement
  par = SinglePar(pname,pvalue,pstate);		// Parameter for final state
  pset.push_back(par);				// Add parameter to set
  
  pstate = Axesf+" Into "+Axesi+" (Degrees)";	// Paraemter statement
  string pbi("FrmRot(");			// Parameter base name
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
*/

// ----------------------------------------------------------------------------
//      Functions To Make A ASCII Parameter File From A Powder
// ----------------------------------------------------------------------------

/* These functions will construct an ASCII parameter file from a coordinate 
   frame.  Individual rotations will be placed into the file parameters as
   Euler angle in degrees.                                                   */

/*
void Powder::write(const string& filename, int pfx) const
   {
   ofstream ofstr(filename.c_str());	// Open filename for input
   if(!ofstr.good())                    // If file bad then exit
     {
     POWerror(1, filename);		// Filename problems
     POWfatal(40);			// Cannot write to file
     }
   ofstr.close();			// Close it now
   ParameterSet pset;			// Declare a parameter set
   PSetAdd(pset,pfx);			// Add ourself into parameter set
   pset.write(filename);		// Write parameter set to filename
   return;
   }
*/

// ____________________________________________________________________________
// D                      POWDER INPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions are used to specify the powder from either a 
   working GAMMA parameter set or an external ASCII file (in parameter set
   format).  The functions do NOT allow for the suffix (#), since it is 
   used/reserved for specific object rotation indices. The functions do allow
   for the prefix [#] so that multiple powders may be defined in the
   same file by switching the prefix indices.

           Input                Pow    : Powder
                                file    : Input file name
                                pset    : Parameter set
                                argv    : Vector of argc arguments
                                argn    : Argument index
                                pfx     : Coordinate frame index (default -1->none)
           Output               none    : Coordinate frame is read in
                                          from parameters in file filename
                                file    : Name of file used to set frame    */

/*
bool Powder::read(const string& filename, int pfx, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename, warn?1:0))	// Read in pset from file
    {
    POWerror(1, filename,1);		// Filename problems
    if(warn > 1) POWfatal(41);		// Cant read from file, fatal 
    else         POWerror(41);		// or non-fatal warning
    return false;
    }
  if(!read(pset,pfx,warn?1:0))		// User overloaded function
    {
    POWerror(1, filename,1);		// Filename problems
    if(warn > 1) POWfatal(42);		// Cannot read file parameters
    else         POWerror(42);		// or non-fatal one
    return false;
    }
  return true;
  }

bool Powder::read(const ParameterSet& pset, int pfx, int warn)
  { 
  bool TF = SetCoordFrm(pset,pfx,warn?1:0);	// Try & set ourself up
  if(!TF)                                       // If we didn't handle
    {                                           // setting ourself properly
    if(warn)                                    // then we'll issue some
      {                                         // warnings if desired
                   POWerror(33, 1);		//   Cant read parameters
      if(warn > 1) POWfatal(34);		//   Can't set from parameters
      else         POWerror(34,1);		//   or a warning issued
      }
    return false; 				// Return that we failed
    }
  return TF;
  }

string Powder::ask_read(int argc, char* argv[], int argn, int idx)
  {
  string filename;				// Name of parameter file
  query_parameter(argc, argv, argn,		// Get filename from command
  "\n\tPowder Filename? ", filename); // Or ask for it
  read(filename, idx); 				// Read rotation from filename
  return filename;
  }

string Powder::ask_read(int argc, char* argv[], int argn, const string& def, int idx)
  {
  string msg = "\n\tPowder Filename ["// Query we will ask if
             + def + "]? ";                     // it is needed
  string filename = def;                        // Name of spin system file
  ask_set(argc,argv,argn,msg,filename);         // Or ask for it
  read(filename, idx);                          // Read system from filename
  return filename;                              // Return filename
  }
*/

// ____________________________________________________________________________
// E                  POWDER FORMATTED OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions will output information concerning the powder to any
   specified output stream.

           Input                Pow    : Powder (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Sum rotation only
                                           !0 - Full composite rotation
           Output               none    : Powder information
                                          placed into the output stream      */


ostream& Powder::print(ostream& ostr, int fflag) const
  {
  int norient = size();
  if(!norient)
    {
    string hdr("Empty Powder\n");
    string spacer = string(40-hdr.length()/2, ' ');
    ostr << "\n\n" << spacer << hdr << "\n";
    return ostr;
    }

/*          Output Some Printed Lines Which Will Look Like The Following

                                 Powder PowderName

  Index   Scale    Theta     Phi       Index   Theta     Phi       Index   Theta     Phi
 ------- -------- -------- -------   -------- --------   ------- -------- --------
   ...     ...      ...      ...     ...      ...        ...     ...      ...  
   ...     ...      ...      ...     ...      ...        ...     ...      ...  
   ...     ...      ...      ...     ...      ...        ...     ...      ...  
   ...     ...      ...      ...     ...      ...        ...     ...      ...    */

//                                This Outputs The Header

  string hdr    = _PN + " Powder";
  string spacer = string((40-hdr.length()/2), ' ');
  ostr << "\n" << spacer << hdr;

  hdr    = "(" + Gdec(norient) + " Orientations";
  if(fflag == norient) hdr += ")";
  else                 hdr += ", " + Gdec(fflag) + " Output)"; 
  spacer = string((40-hdr.length()/2), ' ');
  ostr << "\n" << spacer << hdr << "\n" ;

//                                 This Sizes The Index Column

  norient = min(norient, fflag);
  string csp(" ");
  string I("Index");				// String for Index column header
  string IUL;					// String for Index column underline
  int ilen = Gdec(norient).length();		// Length of largest index
  int len  = I.length();			// Length of Index column header
  int il1=0, il2=0;
  string sil1, sil2;
  if(len >= ilen)				// If column header larger than number
    {						// then just center number under header
    IUL = string(len, '-');
    il1 = (len-ilen)/2;
    if(il1) sil1 = string(il1, ' ');
    il2 = len - il1 - ilen;
    if(il2>0) sil2 = string(il2, ' ');
    }
  else
    {
    IUL = string(ilen, '-');
    I   = string((ilen-len)/2, ' ') + I;
    I  += string(ilen-I.length(), ' ');
    }

//                               This Set Up The Other Columns

  string THE(" Theta  ");
  string PHI("  Phi   ");
  string ALP(" Alpha  ");
  string BET(" Beta   ");
  string GAM(" Gamma  ");
  string AUL("--------");
  string SCL(" Scale ");
  string SUL("-------");

//			          This Writes Column Headers
//              Two Columns If 3 Angle Powder, Three Columns If 2 Angle Powder

  int nc = 2;
  if(!_A3)
    {
    nc = 2;
    ostr << "\n" << I   << csp << SCL << csp << THE << csp << PHI << csp
         <<  " " << I   << csp << SCL << csp << THE << csp << PHI << csp
         <<  " " << I   << csp << SCL << csp << THE << csp << PHI;
    ostr << "\n" << IUL << csp << SUL << csp << AUL << csp << AUL << csp
         <<  " " << IUL << csp << SUL << csp << AUL << csp << AUL << csp
         <<  " " << IUL << csp << SUL << csp << AUL << csp << AUL;
    }
  else
    {
    nc = 1;
    ostr << "\n" << I   << csp << SCL << csp << ALP << csp << BET << csp << GAM << csp
         <<  " " << I   << csp << SCL << csp << ALP << csp << BET << csp << GAM;
    ostr << "\n" << IUL << csp << SUL << csp << AUL << csp << AUL << csp << AUL << csp
         <<  " " << IUL << csp << SUL << csp << AUL << csp << AUL << csp << AUL;
    }

  int newcol = nc;
  for(int i=0; i<norient; i++)
    {
    if(newcol==nc) { newcol=0;  ostr << "\n"; }
    else           { newcol++;  ostr <<  " "; }
   
    ostr << sil1 << Gdec(i, ilen) << sil2 << csp;
    ostr << Gform("%7.5f", _Ps[i]) << csp;
    if(!_A3)
      {
      ostr << Gform("%8.4f", _EAs[i].beta()*RAD2DEG)  << csp;
      ostr << Gform("%8.4f", _EAs[i].alpha()*RAD2DEG) << csp;
      }
    else
      {
      ostr << Gform("%8.4f", _EAs[i].alpha()*RAD2DEG) << csp;
      ostr << Gform("%8.4f", _EAs[i].beta()*RAD2DEG)  << csp;
      ostr << Gform("%8.4f", _EAs[i].gamma()*RAD2DEG) << csp;
      }
    }
  ostr << "\n\n";
  return ostr;
  }

ostream& operator<< (ostream& out, const Powder& Pow) { return Pow.print(out); }

// ____________________________________________________________________________
// F                             KNOWN POWDERS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//                              Spherical Powder
//-----------------------------------------------------------------------------

/* This function generates the simplest powder, just linear angle increments
   over the entire sphere. There will be Nbeta beta angles taken down from the
   +z axis spanning the range of (0, PI) or from (0, PI/2]. At each of these
   angles there will be Nalpha alpha angles taken over from the +x axis. These 
   second angles span the range of [0, 2PI). Since the probablity of having
   the orientation with beta = 0 and beta = PI is zero, these angles will NOT
   be included in the Euler angle set. The flag uhonly dictates if the beta
   range stays in the upper-hemisphere or not, the former allowed only if the
   number of alpha angles is set to 1. Having Nalpha=1 and uhonly=true should
   only be used when there is no asymmetry in the powder. Lastly, since the
   number of alpha angles is constant regardless of the angle beta, there is
   a scaling factor of sin(beta) applicable to each orientation at beta.
   This will clearly be a very direct, abeit computationally slow, means of
   generating an effective set of orientations overa a powder.               */

void Powder::Spherical(int Nbeta, int Nalpha, bool uhonly)
  {
  _EAs.clear();					// Remove existing angles
  _Ps.clear();					// Remove existing probabilties
  _PN = string("Spherical");			// Name powder spherical
  _A3 =false;					// It is two angle powder
  if(Nalpha) uhonly = false;			// Alpha average full sphere
  double bmax   = uhonly?90.0:180.0;		// Maximum beta to compute
  int    nbspan = uhonly?Nbeta:Nbeta+1;		// Number beta points to span 
  double binc   = bmax/double(nbspan); 		// Beta angle increment
  double beta   = binc;				// Start at first beta (not 0)
  double ainc   = 360.0/double(Nalpha); 
  double alph   = 0;				// Start at first alpha (0) 
  double sf;					// Scaling factor
  int i,j;
  for(i=0; i<Nbeta; i++, beta+=binc)		// Loop over beta angles (theta)
    {
    sf = sin(beta*DEG2RAD);			//   Use this  probability
    if(uhonly && beta==90.0) sf/=2.0;		//   Half if at 90 & half-sphere
    _EAs.push_back(EAngles(0,beta,0,true));	//   Store these angles
    _Ps.push_back(sin(beta*DEG2RAD));		//   Store this  probability
    for(j=1; j<Nalpha; j++)			//   Loop alpha angles (phi)
      {
      alph = double(j)*ainc;			//     The angle alpha
      _EAs.push_back(EAngles(alph,beta,0,true));//     Store these angles
      _Ps.push_back(sin(beta*DEG2RAD));		//     Store this  probability
      }
    }
  }

/* As is done above, this function also generates a simplest powder using
   linear angle increments over the entire sphere. There will be Nbeta beta
   angles taken down from the +z axis spanning the range of (0, PI) or from
   (0, PI/2]. At each of these angles there will be sin(beta)*Nalpha alpha
   angles taken over from the +x axis. These second angles span the range 
   of [0, 2PI). Since the probablity of having the orientation with 
   beta = 0 and beta = PI is zero, these angles will NOT be included in the
   powder array of Euler angles. The flag uhonly dictates if the beta range
   stays in the upper-hemisphere or not, the former allowed only if the
   number of alpha angles is set to 1. Having Nalpha=1 and uhonly=true should
   only be used when there is no asymmetry in the powder. Lastly, since the
   number of alpha angles is adjusted depending upon the angle beta, the
   scaling factor of sin(beta) is automatically included. Hence each 
   orientation will have the same probability.                              */

void Powder::Spherical2(int Nbeta, int Nalpha, bool uhonly)
  {
  _EAs.clear();					// Remove existing angles
  _Ps.clear();					// Remove existing probabilties
  _PN = string("Spherical2");			// Name powder spherical
  _A3 =false;					// It is two angle powder
  if(Nalpha) uhonly = false;			// Alpha average full sphere
  double bmax   = uhonly?90.0:180.0;		// Maximum beta to compute
  int    nbspan = uhonly?Nbeta:Nbeta+1;		// Number beta points to span 
  double binc   = bmax/double(nbspan); 		// Beta angle increment
  double beta   = binc;				// Start at first beta (not 0)
  double sf     = 1.0;				// Orientation scaling
  double ainc;  				// Alpha angle increment

  double alph, beta;				// Orientation Euler angles
  double asteps;				// Alpha angles at beta
  for(int j=1, i=0; j<=Nbeta; j++)		// Loop over beta angles
    {
    beta = double(j)*binc;			//   Set beta angle (deg)
    if(uhonly && beta==90.0) sf/=2.0;	        //   Half @ 90 & half-sphere
    asteps = NAlpha*sin(beta*DEG2RAD);		//   # alpha angles this beta
    ainc   = double(360.0)/double(asteps);	//   Alpha angle increment (deg)
    for(i=0; i<asteps; i++)			//   Loop over alpha angles
      {
      alph = ainc*double(i);			//     The angle alpha
      _EAs.push_back(EAngles(alph,beta,0,true));//     Store these angles
      _Ps.push_back(sf);			//     Store this  probability
      }
    }
  }

//-----------------------------------------------------------------------------
//                           Cheng Two Angle Powder
//-----------------------------------------------------------------------------

/* These functions set the powder according the the article by Cheng et. al.:
                                                                          
"Investigations of a nonrandom numerical method for multidimensional integration",
         Vera B. Cheng, Henry H. Suzukawa Jr. and Max Wolfsberg,
                      J.Chem.Phys, 59, 3992-9, (1973)                           

   The array pa contains values related to the powder orientation angle beta,
   or equivalently the angle theta. Given a powder quality value, PAQ, there
   will be N = pa[PAQ] crystallite orientations set in powder. For the kth 
   orientation, the angle beta will be given by PI*(k/N)..... this is
   just evenly incremented down from the z-axis within [0, PI) as in a
   spherical powder. The orientation probablity will be sin(beta), again as
   in a spherical powder. Of course the first value, k=0, will never be used
   since sin(theta) = 0.  The 2nd orientation angle, phi (or equivalently
   alpha), for the kth orientation is given by a somewhat more cryptic formula.
   For this powder, alpha = 2*PI*{(pb[PAQ]*k) % N}/N. This odd incrementation
   is of course the reason the Cheng powder differs from a spherical one.    */

void Powder::Cheng(int PAQ)
  {
  _EAs.clear();					// Remove existing angles
  _Ps.clear();					// Remove existing probabilties
  _PN = string("Spherical");			// Name powder spherical
  _A3 =false;					// It is two angle powder
  if(PAQ > 22) PAQ = 22;			// Insure max quality of 22
  if(PAQ < 0)  PAQ = abs(PAQ);			// Insure quality is positive

  int pa[] = {  2,3,5,8,13,21,34,55,89,144,233,377,610,987,1597,
                        2584,4181,6765,10946,17711,28657,43368,75025};
  int pb[] = {  1,1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,
                        1597,2584,4181,6765,10946,17711,28657};
  int N       = pa[PAQ];
  int M       = pb[PAQ];
  double binc = 180.0/double(N);		// Beta increment
  double ainc = 360.0*double(M)/double(N);	// Alpha increment
  double alpha, beta;
  double sf;
  for(int i=1; i<N; i++)			// Loop over orientations
    {
    beta  =  binc * double(i);			//   Angle beta  = theta
    alpha =  ainc * (i%N);		//   Angle alpha = phi
    sf    = sin(beta*DEG2RAD);			//   Scaling factor
    _EAs.push_back(EAngles(alpha,beta,0,true));	//   Store these angles
    _Ps.push_back(sin(beta*DEG2RAD));		//   Store this  probability
    }
  }

/*
void Powder::Cheng(int PAQ)
  {
  int value1[] = {     2,     10,    20,    30,     40,
                      50,    100,   150,   200,    300,
                     400,    500,   600,   700,    800,
                     900,   1000,  1100,  1200,   1300,
                    1400,   1500,  1600,  1700,   1800,
                    1900,   2000,  2100,  2200,   2300,
                    2400,   2500, 10000, 50000, 100000};
  int value2[] = {     1,      3,     3,     7,      3,
                       5,     15,    35,    55,     89,
                     127,     97,   103,   145,    189,
                     233,    313,   523,   573,    447,
                     159,    139,   205,   321,    291,
                     671,    297,   395,   697,    527,
                     549,    363,  3189,  9027,  27205};
  int value3[] = {     1,      5,     7,    11,     15,
                      13,     47,    63,    81,    137,
                     187,    229,   265,   223,    257,
                     355,    477,   391,   181,    191,
                     553,    621,   551,   789,    481,
                     829,    479,   993,   887,    827,
                     841,    917,  4713, 14857,  38057};
  }
*/


#endif							// Powder.cc
