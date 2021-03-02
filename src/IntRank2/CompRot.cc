/* CompRot.cc ***************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Composite Rotation 			Implementation		**
**                                                                      **
**      Scott Smith                                                     **
**      Copyright (c) 2001                                              **
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
** A variable of type rotation maintains a list of successive rotations **
** in terms of both Euler angles and Quaternions. The class contains    **
** two vectors of equal length, one containing a list of Euler angles   **
** and one containing a list of corresponding Quaternions. Access       **
** functions allow users to get or set any of the rotations in the      **
** vector, and to get summed rotations over any of the rotations.       **
**                                                                      **
*************************************************************************/

#ifndef   CompRot_cc_			// Is file already included?
#  define CompRot_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <IntRank2/CompRot.h>		// Include interface definition
#include <Basics/Gutils.h>		// Include GAMMA std errors
#include <Basics/Gconstants.h>		// Include GAMMA constants
#include <Level2/EAngles.h>		// Include Euler angles
#include <Level2/Quaternion.h>		// Include quaterions
#include <Basics/StringCut.h>		// Include form & dec functions
#if defined(_MSC_VER) || defined(__SUNPRO_CC)				// If we are using MSVC++
 #pragma warning (disable : 4786)       //   Kill STL namelength warnings
#endif

using std::string;			// Using libstdc++ strings
using std::vector;			// Using libstdc++ STL vectors
using std::ofstream;			// Using libstdc++ output file streams
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                CLASS COMPOSITE ROTATION ERROR HANDLING
// ____________________________________________________________________________

/*       Input                ROT     : Composite Rotation (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void CompRot::ROTerror(int eidx, int noret) const
  {
  string hdr("Composite Rotation");
  switch(eidx)
    {
    case  0: GAMMAerror(hdr,"Program Aborting.....",       noret);break; // (0)
    case  2: GAMMAerror(hdr,"Problems During Construction",noret);break; // (2)
    case 11: GAMMAerror(hdr,"Cannot Set Any Rotations",    noret);break; //(11)
    case 12: GAMMAerror(hdr,"Insufficient Parameters",     noret);break; //(12)
    case 20: GAMMAerror(hdr,"Accessed Rotation Absent",    noret);break; //(20)
    case 33: GAMMAerror(hdr,"Cannot Read Parameters",      noret);break; //(33)
    case 34: GAMMAerror(hdr,"Cannot Set From Parameters",  noret);break; //(34)
    case 40: GAMMAerror(hdr,"Cannot Write To File",        noret);break; //(40)
    case 41: GAMMAerror(hdr,"Cannot Read From File",       noret);break; //(41)
    case 42: GAMMAerror(hdr,"Cannot Read File Parameters", noret);break; //(42)
    case 50: GAMMAerror(hdr,"Bad Range In Summed Rotation",noret);break; //(50)
    case 51: GAMMAerror(hdr,"Cannot Make Summed Rotation", noret);break; //(51)
    }
  }

volatile void CompRot::ROTfatal(int eidx) const
  {
  ROTerror(eidx, 1);				// Output error message
  if(eidx) ROTerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }


/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Funciton pname
      5     Please Use Class Member Funciton pname
   default  Unknown Error - pname                                           */

void CompRot::ROTerror(int eidx, const string& pname, int noret) const 
  {
  string hdr("Composite Rotation");
  string msg;
  switch(eidx)
    {
    case  21: msg = string("Can't Access Rotation ")     + pname; 	// (21)
              GAMMAerror(hdr, msg, noret); break;
    case 101: msg = string("Can't Find Parameters For ") + pname;	//(101)
              GAMMAerror(hdr, msg, noret); break;
    default:  GAMMAerror(hdr, eidx, pname, noret);         break;
    }
  }

volatile void CompRot::ROTfatal(int eidx, const string& pname) const
  {
  ROTerror(eidx, 1);				// Output error message
  if(eidx) ROTfatal(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                COMPOSITE ROTATION CHECKING FUNCTIONS
// ____________________________________________________________________________


bool CompRot::ChkIdx(int i, int warn) const
  {
  if(i<0 || i>=int(Qs.size()))		// Check that rotation exists
    {					// If not we will fail
    if(warn)				//   If warnings desired we output
      {					//   message(s) and perhapse quit
                 ROTerror(20, 1);	//     Rotation not present
      if(warn>1) ROTfatal(21,Gdec(i));	//     Cannot access rotation i
      else       ROTerror(21,Gdec(i),1);
      }
    return false;			//   Rotation i is not present
    }
  return true;				// AOK, rotation is present
  }

bool CompRot::ChkRange(int i, int nr, int warn) const
  {
  bool TFnr = (nr > 0);			// Insure nr is positive
  bool TFi = ChkIdx(i,warn?1:0);	// Insure rotation i exists
  bool TFj = ChkIdx(i+nr-1,warn?1:0);	// Insure rotation j exists
  if(!TFi || !TFj || !TFnr)		// If either does not exist
    {					// then we have problems
    if(warn)				// If warnings desired we output
      {					// message(s) and perhaps quit
      		 ROTerror(50,1);	//   Bad range summed rotation
      if(warn>1) ROTfatal(51); 		//   Cannot make summed rotaiton
      else       ROTerror(51,1);
      }
    return false;			//   Rotation range not present
    }
  return true;				// AOK, rotation is present
  }

// ____________________________________________________________________________
// iii              COMPOSITE ROTATION AUXILIARY FUNCTIONS
// ____________________________________________________________________________


void CompRot::SetSum()
  {
  int nr = (int)Qs.size();		// Number of rotations
  if(!nr) 				// If no rotations then the
    { 					// summed quaterion is {0,0,0,1}
    sumQ  = quatern();			// & summed Euler angles are { 0,0,0 }
    sumEA = EAngles(coord0);
    }
  else					// If ther are rotations
    {
    sumQ = Qs[0];			//   Start with 1st quaternion
    for(int i=1; i<nr; i++)		//   Loop over all other rotations
      sumQ *= Qs[i];			//   and add them in to get the
    sumEA = sumQ.EA();			//   summed quaternion. The summed
    }					//   Euler angles computed from sum
  }

// ____________________________________________________________________________
// iv            COMPOSITE ROTATION PARAMETER SET SETUP FUNCTIONS 
// ____________________________________________________________________________

/* These functions try and set a composite rotation from parameters found in
   a particular parameter set.  Each rotation can be specified in one of three
   ways:
             1.) As Euler Angles - EAngles: (alpha,beta,gamma) EAngles(#)
             2.) As a Quaternion - quatern: (A,B,C,D)          Quaternion(#)
             3.) As three angles - double:  alpha              EAalpha(#)
                                   double:  beta               EAbeta(#)
                                   double   gamma              EAgamma(#)

   where the # is used to indicate the rotation index. The first rotation will
   have index 0 and the last rotation will be that which is present in 
   successive order from 0.                                                  */

bool CompRot::SetCmpRot(const ParameterSet& pset, int pfx, int warn)
  {
  ParameterSet  subpset;                        // Copy parameter set
  if(pfx != -1) subpset = pset.strip(pfx);      // to glean out those for
  else          subpset = pset;                 // the specified prefix only
  Qs.clear();					// Remove all quaterions
  EAs.clear();					// Remove all Euler angles
  int i=0;					// Temporary rotation index
  while(SetRotation(subpset, i)) { i++; }	// Read successive rotations
  if(Qs.size())					// If we found some, set the
    { SetSum(); return true; }			// summed rotation & return
  else						// If none read, issue warnings
    {						// as desired, perhaps stop
    ROTerror(11, 1);				//   Cant set any rotations
    if(warn>1) ROTfatal(12);			//   Insufficient parameters
    else       ROTerror(12, 1);
    }
  return false;
  }


bool CompRot::SetRotation(const ParameterSet& pset, int idx)
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
  
bool CompRot::GetEulerAngles(const ParameterSet& pset, EAngles& EA, int idx)
  { return EA.read(pset, idx); }

bool CompRot::GetQuaternion(const  ParameterSet& pset, quatern&  Q, int idx)
  { return Q.read(pset, idx); }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                   COMPOSITE ROTATION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________


/* Rotation construction. Since each rotation may contain a series of single
   rotations we allow use arrays for their construction. Thus the rotation may
   be constructed empty (no rotations), with a single rotation (either using
   Euler angles or a Quaternion), or with a vector of rotations (again
   either using Euler angles or Quaternions).                                */

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

CompRot::CompRot() : EAs(), Qs(), sumEA(coord0), sumQ() {}

// ----------------------------------------------------------------------------
//                   Constructors Using Single Rotations
// ----------------------------------------------------------------------------

CompRot::CompRot(const EAngles& EA) : EAs(1,EA),Qs(1,quatern(EA)) { SetSum(); }
CompRot::CompRot(const quatern& Q)  : EAs(1,Q.EA()),Qs(1,Q)       { SetSum(); }
CompRot::CompRot(double A, double B, double G)
         : EAs(1, EAngles(A,B,G)), Qs(1, quatern(EAs[0]))         { SetSum(); }

// ----------------------------------------------------------------------------
//                  Constructors Using Multiple Rotations
// ----------------------------------------------------------------------------

CompRot::CompRot(const vector<EAngles>& EAvec)
  {
  EAs = EAvec;
  int nr = int(EAs.size());
  Qs  = vector<quatern>(nr);
  for(int i=0; i<nr; i++)
    Qs[i] = quatern(EAs[i]);
  SetSum();
  }

CompRot::CompRot(const vector<quatern>& Qvec)
  {
  Qs = Qvec;
  int nr = int(Qs.size());
  EAs = vector<EAngles>(nr);
  for(int i=0; i<nr; i++)
    EAs[i] = Qs[i].EA();
  SetSum();
  }

// ----------------------------------------------------------------------------
//                    Construction Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions will construct the rotation from parameters found in a
   specified GAMMA parameter set. This allows the rotation to be generated 
   from individual rotations as specified in an external ASCII file.         */

CompRot::CompRot(const ParameterSet& pset, int idx, int warn)
  { 
  bool TF = SetCmpRot(pset, idx, warn?1:0);	// Try to construct
  if(!TF && warn)				// If construction fails &
    { 						// warnings desired 
    ROTerror(2, 1);				//   Problems in construction
    if(warn>1) ROTfatal(22,Gdec(idx));		//   Cant construct rotation
    else       ROTerror(22,Gdec(idx),1);
    }
  }


// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

void CompRot::operator= (const CompRot& ROT) 
  {
  EAs   = ROT.EAs;			// Copy vector of Euler angles
  Qs    = ROT.Qs;			// Copy vector of Quaternions
  sumEA = ROT.sumEA;			// Copy summed rotation Euler angles
  sumQ  = ROT.sumQ;			// Copy summed rotation Quaternion
  }

CompRot::~CompRot () { }


// ____________________________________________________________________________
// B                  INDIVIDUAL COMPOSITE ROTATION ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions allow users to either get or set individual rotations in
   the vector of rotations.                                                  */
 
EAngles CompRot::EA(int    i) const { return EAs[i];           }
quatern CompRot::Q(int     i) const { return Qs[i];            }
double  CompRot::alpha(int i) const { return (EAs[i]).alpha(); }
double  CompRot::beta(int  i) const { return (EAs[i]).beta();  }
double  CompRot::gamma(int i) const { return (EAs[i]).gamma(); }

void CompRot::EA(const  EAngles& ea, int i) { EAs[i]     = ea;   SetSum(); }
void CompRot::Q(const   quatern&  q, int i) { Qs[i]      = q;    SetSum(); }
void CompRot::alpha(double        A, int i) { (EAs[i]).alpha(A); SetSum(); }
void CompRot::beta(double         B, int i) { (EAs[i]).beta(B);  SetSum(); }
void CompRot::gamma(double        G, int i) { (EAs[i]).gamma(G); SetSum(); }

// ____________________________________________________________________________
// C                    SUMMED COMPOSITE ROTATION ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These allow users to access composite rotations, either over the entire
   vector or any sub-vector. Since composite rotations may depend upon
   multiple single rotations they may not be set, only obtained.             */
 
EAngles CompRot::EA()    const { return sumEA;         }
quatern CompRot::Q()     const { return sumQ;          }
double  CompRot::alpha() const { return sumEA.alpha(); }
double  CompRot::beta()  const { return sumEA.beta();  }
double  CompRot::gamma() const { return sumEA.gamma(); }

EAngles CompRot::EA(int    i, int nr) { return (Q(i,nr)).EA();     }
double  CompRot::alpha(int i, int nr) { return (EA(i,nr)).alpha(); }
double  CompRot::beta(int  i, int nr) { return (EA(i,nr)).beta();  }
double  CompRot::gamma(int i, int nr) { return (EA(i,nr)).gamma(); }
quatern CompRot::Q(int     i, int nr) 
  {
  ChkRange(i, nr, 2);			// Insure range is proper
  quatern Q = Qs[i];			// Start with rotation i
  for(int k=1; k<nr; k++) Q *= Qs[i+k];	// Add in successive rotations
  return Q;				// This is the summed rotation
  }

// ____________________________________________________________________________
// D                  COMPOSITE ROTATION PARSING FUNCTIONS
// ____________________________________________________________________________

/* These function allow users to obtain any part of the composite rotation   */

CompRot CompRot::operator() (int i, int nr) const
  {
  ChkRange(i, nr, 2);			// Insure range is proper
  CompRot CR(nr);			// Make an empty rotation
  for(int j=0, k=0; k<nr; k++, j++)	// Loop over rotations desired
    {
    CR.EAs[j] = EAs[k];			//   Clip Euler angle list
    CR.Qs[j]  = Qs[k];			//   Clip Quaternion list
    }
  CR.SetSum();				// Set new summed rotation
  return CR;
  }

void CompRot::operator+= (const EAngles& EA)
  {
  EAs.push_back(EA);
  Qs.push_back(quatern(EA));
  sumQ *= Qs[Qs.size()-1];
  sumEA = sumQ.EA(); 
  }

void CompRot::operator+= (const quatern& Q)
  {
  Qs.push_back(Q);
  EAs.push_back(Q.EA());
  sumQ *= Q;
  sumEA = sumQ.EA(); 
  }

// ____________________________________________________________________________
// E                        PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//            Functions To Make A Parameter Set From A Rotation
// ----------------------------------------------------------------------------

/* These functions will construct a parameter set from a composite rotation.
   Because an individual rotation may be specified either in terms of three
   Euler angles or as a single Quaterion we must choose one or the other for
   output output.  By default, these will be done as 3 Euler angles in the
   form of a GAMMA coordinate.                                               */

CompRot::operator ParameterSet( ) const
   { ParameterSet pset; pset += *this; return pset; }	

void operator+= (ParameterSet& pset, const CompRot& ROT)
  { ROT.PSetAdd(pset); }

void CompRot::PSetAdd(ParameterSet& pset, int pfx) const
  {
  string prefx;                                 // Parameter prefix
  if(pfx != -1)                                 // Only use prefix if pfx
    prefx = string("[")+Gdec(pfx)+string("]");	// is NOT -1
  SinglePar par;				// Temporary single parameter
  string pstate("Euler Angles (deg)");		// Parameter statement
  string pst = prefx + string("EAngles(");	// Parameter name base
  string pend(")");				// End of parameter name
  string pname;					// For full parameter name
  coord ABG;					// Temporary coordinate
  int nr = (int)Qs.size();			// Number of rotations
  for(int i=0; i<nr; i++)			// Loop over individual
    {						// rotations
    pname = pst + Gdec(i) + pend;		//   Parameter name ith rot.
    ABG.x(EAs[i].alpha()*RAD2DEG);		//   We use a coordinate for
    ABG.y(EAs[i].beta()*RAD2DEG);		//   output and set the angles
    ABG.z(EAs[i].gamma()*RAD2DEG);		//   to be in degrees
    par = ABG.param(pname, pstate);		//   Parameter for Euler angles
    pset.push_back(par);			//   Add parameter to set
    } 
  }

void CompRot::write(const string &filename, int pfx) const
   {
   ofstream ofstr(filename.c_str());	// Open filename for input
   if(!ofstr.good())                    // If file bad then exit
     {
     ROTerror(1, filename);		// Filename problems
     ROTfatal(40);			// Cannot write to file
     }
   ofstr.close();			// Close it now
   ParameterSet pset;			// Declare a parameter set
   PSetAdd(pset,pfx);			// Add ourself into parameter set
   pset.write(filename);		// Write parameter set to filename
   return;
   }

// ____________________________________________________________________________
// H                     COMPOSITE ROTATION INPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions are used to specify the rotation from either a working
   GAMMA parameter set or an external ASCII file (in parameter set format).
   The functions do NOT allow for the suffix (#), since it is used/reserved
   specific rotation indices. The functions do allow for the prefix [#] so that
   multiple composite rotations may be defined in the same file by switching
   the prefix indices. 

           Input                ROT     : Composite Rotation
                                file    : Input file name
                                pset    : Parameter set
                                argv    : Vecotr of argc arguments
                                argn    : Argument index
                                pfx     : Rotation prefix (default -1->none)
           Output               none    : Composite Rotation is read in
                                          from parameters in file filename
                                file    : Name of file used to set rotation  */

bool CompRot::read(const string &filename, int pfx, int warn)
   {
   ParameterSet pset;			// Declare a parameter set
   if(!pset.read(filename, warn?1:0))	// Read in pset from file
     {
     ROTerror(1, filename,1);		// Filename problems
     if(warn > 1) ROTfatal(41);		// Cant read from file, fatal 
     else         ROTerror(41);		// or non-fatal warning
     return false;
     }
   if(!read(pset,pfx,warn?1:0))		// User overloaded function
     {
     ROTerror(1, filename,1);		// Filename problems
     if(warn > 1) ROTfatal(42);		// Cannot read file parameters
     else         ROTerror(42);		// or non-fatal one
     return false;
     }
   return true;
   }

bool CompRot::read(const ParameterSet& pset, int pfx, int warn)
  { 
  bool TF = SetCmpRot(pset,pfx,warn?1:0);	// Try & set ourself up
  if(!TF)                                       // If we didn't handle
    {                                           // setting ourself properly
    if(warn)                                    // then we'll issue some
      {                                         // warnings if desired
                   ROTerror(33, 1);		//   Cant read parameters
      if(warn > 1) ROTfatal(34);		//   Can't set from parameters
      else         ROTerror(34,1);		//   or a warning issued
      }
    return false; 				// Return that we failed
    }
  return TF;
  }

string CompRot::ask_read(int argc, char* argv[], int argn, int idx)
  {
  string filename;				// Name of parameter file
  query_parameter(argc, argv, argn,		// Get filename from command
       "\n\tRotation filename? ", filename); 	// Or ask for it
  read(filename, idx); 				// Read rotation from filename
  return filename;
  }

// ____________________________________________________________________________
// J                  COMPOSITE ROTATION FORMATTED OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions will output information concerning the Rotation to any 
   specified output stream.

           Input                ROT     : Composite Rotation (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Sum rotation only
                                           !0 - Full composite rotation
           Output               none    : Composite Rotation placed
                                          into the output stream            */

ostream& CompRot::print(ostream& ostr, int fflag) const
  {
  if(!Qs.size())
    {
    ostr << "\t\n\tEmpty Composite Rotation\n";
    return ostr;
    }
  int nrot = (int)Qs.size();
  if(nrot == 1 || !fflag)
    ostr << "\n\t" << sumEA << "\t" << sumQ;
  else
    {
    ostr << "\n                                Composite Rotation";
    ostr << "\nRotation           Euler Angles"
         << "                         Quaternions";
    ostr << "\n Index     "
         << " alpha      beta      gamma    "
         << "    A         B         C         D";
    ostr << "\n"
         << "--------   --------  --------  --------   "
         << "--------  --------  --------  --------";
    for(int i=0; i<nrot; i++)
      {
      ostr << "\n   " << Gdec(i) << "       ";
      ostr << Gform("%8.4f", EAs[i].alpha()) << "  "
           << Gform("%8.4f", EAs[i].beta())  << "  "
           << Gform("%8.4f", EAs[i].gamma()) << "   ";

      ostr << Gform("%8.5f", Qs[i].A()) << "  "
           << Gform("%8.5f", Qs[i].B()) << "  "
           << Gform("%8.5f", Qs[i].C()) << "  "
           << Gform("%8.5f", Qs[i].D());
      }
//                      Last Row Is Summed Rotation
    ostr << "\n" << "  sum      ";
    ostr << Gform("%8.4f", sumEA.alpha()) << "  "
         << Gform("%8.4f", sumEA.beta())  << "  "
         << Gform("%8.4f", sumEA.gamma()) << "   ";
    ostr << Gform("%8.5f", sumQ.A()) << "  "
         << Gform("%8.5f", sumQ.B()) << "  "
         << Gform("%8.5f", sumQ.C()) << "  "
         << Gform("%8.5f", sumQ.D());
    }
  ostr << "\n\n";
  return ostr;
  }

ostream& operator<< (ostream& out, const CompRot& ROT) { return ROT.print(out); }

#endif							// CompRot.cc
