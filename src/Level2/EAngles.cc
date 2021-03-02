/* EAngles.cc ***************************************************-*-c++-*-
**									**
** 	                         G A M M A				**
**									**
**	Euler Angles		                    Implementation      **
**						 			**
**	Copyright (c) 2001					 	**
**	Scott Smith				 			**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                    		  	**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**                                                                      **
** This file contains a set of three Euler angles {alpha, beta, gamma}  **
** that are used in general rotations.  The class structure exists in   **
** order to maintain the angles within specific ranges as well as to    **
** facilite I/O and composite rotations.                                **
**                                                                      **
** The euler angle ranges maintained in this class are:                 **
**                                                                      **
**          alpha = gamma = [0,360]       beta = [0,180]                **
**                                                                      **
** This causes no loss in generality because rotations may still cover  **
** all of three dimensional space.                                      **
**                                                                      **
** Also note that all angles are maintained in units of radians herein. **
** Only Input/Output allows for degrees if desired.  The one exception	**
** is the constructor that takes a coordinate ABG in which it is	**
** assumed the angles are also specified in degrees.			**
**                                                                      **
*************************************************************************/

#ifndef _EAngles_cc_			// Is file already included?
#define _EAngles_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// this is the implementation
#endif

#if defined(_MSC_VER)			// If using MSVC++ then we
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif 

#include <Level2/EAngles.h>		// Include the header file
#include <Basics/Gconstants.h>		// Need constant PI
#include <Basics/Gutils.h>		// Need GAMMA error messages
#include <Level2/Quaternion.h>		// Include Quaternions
#include <Level1/coord.h>		// Include coordiantes 
#include <Basics/StringCut.h>		// Include Gdec, Gform (dec & form)
#include <stdlib.h>

using std::string;			// Using libstdc++ strings
using std::cout;			// Using libstdc++ standard output
using std::list;			// Using libstdc++ STL lists
using std::ofstream;			// Using libstdc++ output file streams
using std::ostream;			// Using libstdc++ output streams

const  EAngles EAzero(0,0,0);		// Zero Euler angles (constant)
double EAngles::AngCutoff = 1.e-10;	// Zero cutoff for angle (radians)

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE Functions ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                   Class Euler Angles Error Handling
// ____________________________________________________________________________

/*      Input                   EAx     : Euler Angles (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
        Output                  none    : Error message output
                                          Program execution stopped if fatal */
 
void EAngles::EAerror(int eidx, int noret) const
  {
  string hdr("Euler Angles");
  string msg;
  switch (eidx)
    {
    case 33: msg = string("Not Specified From Given Parameters");      // (33)
             GAMMAerror(hdr,msg,noret); break;
    case 34: msg = string("Cannot Properly Set From Param. Set");      // (34)
             GAMMAerror(hdr,msg,noret); break;
    case 40: msg = string("Can't Write To Parameter Set File");        // (40)
             GAMMAerror(hdr,msg,noret); break;
    case 41: msg = string("Can't Read Parameters From File");          // (41)
             GAMMAerror(hdr,msg,noret); break;
    case 42: msg = string("Can't Set From External File");             // (42)
             GAMMAerror(hdr,msg,noret); break;
    default: GAMMAerror(hdr,eidx,noret); break;
    }
  }
 
void EAngles::EAerror(int eidx, const string& pname, int noret) const
  {
  string hdr("Euler Angles");
  string msg;
  switch(eidx)
    {
    case  40: msg = string("Can't Write To File ")       + pname;       // (40)
              GAMMAerror(hdr, msg, noret); break;
    default:  GAMMAerror(hdr, eidx, pname, noret);         break;
    }
  }

volatile void EAngles::EAfatal(int eidx) const
  {
  EAerror(eidx, 1);
  if(eidx) EAerror(0);
  cout << "\n";
  exit(-1);
  }

// ____________________________________________________________________________
// ii            Class Euler Angles Private Facilitator Functions
// ____________________________________________________________________________

/* These are checking functions for Euler angles. They insure that the angles
   always reside within the class specified ranges:

                 alpha = gamma = [0,360)       beta = [0,180]
   
   Some care must be taken when the user inputs one or more angles that lie
   outside of these ranges.  For alpha and gamma one may simply add/subtract
   n*2*PI until the angles are within their specified ranges. For angle beta,
   since it does not span [0, 360), this does not entirely work! We can start
   by setting beta to be within [-PI, PI] & if it lies within [0, PI] all
   is OK. However, if beta lies within [-PI, 0) then the actual beta rotation
   angle should be |beta|.... but this does not account for the proper
   rotation of the x & y axes, their direction of rotation is reversed. To
   compensate we must add an additional PI rotation to both alpha & gamma.   */

void EAngles::SetAngles(double alpha, double beta, double gamma, bool deg)
  {
  if(deg)					// If angles input in degrees
    {						// convert them to radians
    alpha *= DEG2RAD;				// ALL Euler angles are kept
    beta  *= DEG2RAD;				// in radians (internally)
    gamma *= DEG2RAD;
    }
  _alpha = fmod(alpha, PIx2);			// Insure alpha = [0, 360)
  while(_alpha <  0.0)  { _alpha += PIx2; }	//   If < 0,    add n*2*PI
  while(_alpha >= PIx2) { _alpha -= PIx2; }	//   If > 2*PI, sub n*2*PI
  if(fabs(_alpha)      < AngCutoff) _alpha=0;	//   If ~0,    make exactly 0
  if(fabs(_alpha-PIx2) < AngCutoff) _alpha=0;	//   If ~2*PI, make exactly 0

  _gamma = fmod(gamma, PIx2);			// Insure beta = [0, 360)
  while(_gamma < 0.0)   { _gamma += PIx2; }	//   If < 0,    add n*2*PI
  while(_gamma >= PIx2) { _gamma -= PIx2; }	//   If > 2*PI, sub n*2*PI
  if(fabs(_gamma)      < AngCutoff) _gamma=0;	//   If ~0,    make exactly 0
  if(fabs(_gamma-PIx2) < AngCutoff) _gamma=0;	//   If ~2*PI, make exactly 0 

//  _beta = fmod(beta, PIx2);			// Insure beta = [0, 360)
  _beta = beta;
  while(_beta < -PI)   { _beta += PIx2; }	//   If < -PI, add n*2*PI
  while(_beta >  PI)   { _beta -= PIx2; }	//   If >  PI, add n*2*PI
  if(fabs(_beta)      < AngCutoff) _beta=0;	//   If ~0, make exactly 0
  if(fabs(_beta-PI)   < AngCutoff) _beta=PI;	//   If ~ PI, make exactly PI
  if(_beta >= 0) return;			// If beta = [0, PI] done
  _beta = -_beta;				// If beta = [-PI,0) 
  if(fabs(_beta)      < AngCutoff) _beta=0;	//   If ~0, make exactly 0
  if(fabs(_beta-PI)   < AngCutoff) _beta=PI;	//   If ~ PI, make exactly PI
  _alpha += PI;
  _alpha = fmod(_alpha, PIx2);			// Insure alpha = [0, 360)
  if(fabs(_alpha)      < AngCutoff) _alpha=0;	//   If ~0,    make exactly 0
  if(fabs(_alpha-PIx2) < AngCutoff) _alpha=0;	//   If ~2*PI, make exactly 0
  _gamma += PI;
  _gamma = fmod(_gamma, PIx2);			// Insure alpha = [0, 360)
  if(fabs(_gamma)      < AngCutoff) _gamma=0;	//   If ~0,    make exactly 0
  if(fabs(_gamma-PIx2) < AngCutoff) _gamma=0;	//   If ~2*PI, make exactly 0
  }

void EAngles::SetAlpha(double alpha, bool deg)
  {
//  if(deg) _alpha = fmod(alpha, 360.0) * DEG2RAD;
//  else    _alpha = fmod(alpha, PIx2);
  _alpha = fmod(alpha, PIx2);
  while(_alpha <  0.0)  { _alpha += PIx2; }
  while(_alpha >= PIx2) { _alpha -= PIx2; }
  if(fabs(_alpha)      < AngCutoff) _alpha=0;
  if(fabs(_alpha-PIx2) < AngCutoff) _alpha=0;
  }

void EAngles::SetBeta(double beta, bool deg)
  {
//  if(deg) _beta = fmod(beta, 180.0) * DEG2RAD;
//  else    _beta = fmod(beta, PIx2);
  _beta = fmod(beta, PIx2);
  while(_beta < 0.0)  { _beta += PI; }
  while(_beta > PI)   { _beta -= PI; }
  if(fabs(_beta)      < AngCutoff) _beta=0;
  if(fabs(_beta-PI)   < AngCutoff) _beta=PI;
  }

void EAngles::SetGamma(double gamma, bool deg)
  {
//  if(deg) _gamma = fmod(gamma, 360.0) * DEG2RAD;
//  else    _gamma = fmod(gamma, PIx2);
  _gamma = fmod(gamma, PIx2);
  while(_gamma < 0.0)   { _gamma += PIx2; }
  while(_gamma >= PIx2) { _gamma -= PIx2; }
  if(fabs(_gamma)      < AngCutoff) _gamma=0;
  if(fabs(_gamma-PIx2) < AngCutoff) _gamma=0;
  }

// ____________________________________________________________________________
// iii          Class Euler Angles Private Parameter Set Functions
// ____________________________________________________________________________

/* These functions try and set up Euler angles from parameters found in
   a particular parameter set.  Each object of type EAngles can be specified
   in one of two ways:

             1.) As Euler Angles - coord:   (alpha,beta,gamma) EAngles(#)
             2.) As three angles - double:  alpha              EAalpha(#)
                                   double:  beta               EAbeta(#)
                                   double   gamma              EAgamma(#)

   where the # is used to indicate the Euler angle index. It is assumed that
   the angles are specified in degrees, although the names may have "Rad"
   appended before "(#)" if the values are specified in radians. Again, we
   strictly enforce the angle ranges of alpha = gamma = [0,360]  and 
   beta = [0,180].                                                           */

bool EAngles::SetEAngles(const ParameterSet& pset, int idx, bool warn)
  {
  bool    TF = SetEASet(pset,   idx, warn);	// Try reading as coordinate
  if(!TF) TF = Set3Angles(pset, idx, warn);	// Try reading as 3 angles
/* sosi
  if(!TF && warn)
    {
    }
*/
  return TF;
  }

bool EAngles::SetEASet(const ParameterSet& pset, int idx, bool warn)
  {
  int nn=2;                                     // # allowed parameter names
  string pnames[2]={ "EAngles", "EAnglesRad" };	// Possible parameter names
  string app = string("(") + Gdec(idx) + ")";	// Parameter index append
  string pname;					// Full parameter name
  ParameterSet::const_iterator item;         // Pix in parameter list
  coord EA;					// Euler angles as coordinate
  bool deg = true;				// Flag if in degrees
  for(int i=0; i<nn; i++)                       // Loop possible parameters
    {
    pname = pnames[i];                          // Consruct parameter name
    if(idx != -1) pname += app;			// Use suffix if not -1 index
    item  = pset.seek(pname);                   // Pix in parameter list
    if(item != pset.end())                      // If the parameter is found
      {
      EA = coord(*item); 			//   Set coordinate from pset
      SetAlpha(deg?EA.x()*DEG2RAD:EA.x());	//   Set our alpha value
      SetBeta(deg?EA.z()*DEG2RAD:EA.y());	//   Set our beta value
      SetGamma(deg?EA.z()*DEG2RAD:EA.z());	//   Set our gamma value
      return true;				//   Return we were successful
      }
    deg = false;				// On 2nd we look in radians
    }
/*
sosi
  if(warn)
    {
    }
*/
  return false;					// Return unsuccessful
  }

bool EAngles::Set3Angles(const ParameterSet& pset, int idx, bool warn)
  {
  string pname;					// Full parameter name
  string pstate;				// Temp parameter comment
  string app;					// Name appendix: (idx)
  if(idx != -1)					// Set appendix only if
    app += string("(") + Gdec(idx) + ")";	// idx is not -1 (default)
  ParameterSet::const_iterator item;         // Pix in parameter list
  _alpha = 0;					// Zero current alpha value
  _beta  = 0;					// Zero current beta value
  _gamma = 0;					// Zero current gamma value

//			First Look For Angle Alpha                      

  bool TFA = false;				// Flag if we find alpha
  double A;					// Temp value for alpha
  pname = string("EAalpha") + app;		// Set parameter name
  item  = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// If the parameter is found
    {						// then set alpha in degrees
    (*item).parse(pname,A,pstate);		//   Get alpha from pset
    SetAlpha(A*DEG2RAD); 			//   Set our alpha value
    TFA = true;					//   Flag we found it
    }
  else						// If parameter not found look
    {						// for alpha in radians
    pname = string("EAalphaRad") + app;		//   Set parameter name
    item  = pset.seek(pname);			//   Pix in parameter list
    if(item != pset.end())			//   If the parameter is found
      {						//   then set alpha in radians
      (*item).parse(pname,A,pstate);		//     Get alpha from pset
      SetAlpha(A); 				//     Set our alpha value
      TFA = true;				//     Flag we found it
      }
    }

//			 Next Look For Angle Beta 

  bool TFB = false;				// Flag if we find beta
  double B;					// Temp value for beta
  pname = string("EAbeta") + app;		// Set parameter name
  item  = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// If the parameter is found
    {						// then set beta in degrees
    (*item).parse(pname,B,pstate);		//   Get beta from pset
    SetBeta(B*DEG2RAD); 			//   Set our beta value
    TFB = true;					//   Flag we found it
    }
  else						// If parameter not found look
    {						// for beta in radians
    pname = string("EAbetaRad") + app;		//   Set parameter name
    item  = pset.seek(pname);			//   Pix in parameter list
    if(item != pset.end())			//   If the parameter is found
      {						//   then set beta in radians
      (*item).parse(pname,B,pstate);		//     Get beta from pset
      SetBeta(B); 				//     Set our beta value
      TFB = true;				//     Flag we found it
      }
    }

//			 Lastly Look For Angle Gamma 

  bool TFG = false;				// Flag if we find gamma
  double G;					// Temp value for gamma
  pname = string("EAgamma") + app;		// Set parameter name
  item  = pset.seek(pname);			// Pix in parameter list
  if(item != pset.end())			// If the parameter is found
    {						// then set gamma in degrees
    (*item).parse(pname,G,pstate);		//   Get gamma from pset
    SetGamma(G*DEG2RAD); 			//   Set our gamma value
    TFG = true;					//   Flag we found it
    }
  else						// If parameter not found look
    {						// for gamma in radians
    pname = string("EAgammaRad") + app;		//   Set parameter name
    item  = pset.seek(pname);			//   Pix in parameter list
    if(item != pset.end())			//   If the parameter is found
      {						//   then set gamma in radians
      (*item).parse(pname,G,pstate);		//     Get gamma from pset
      SetGamma(G); 				//     Set our gamma value
      TFG = true;				//     Flag we found it
      }
    }

  if(TFA || TFB || TFG) return true;		// If any found, flag OK.
/* sosi
  if(warn)
    {
    }
*/
  return false;
  }

// ----------------------------------------------------------------------------
// ---------------------------- Public Functions ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A               Class Euler Angles Constructors/Destructor
// ____________________________________________________________________________

/* These constructors set up a new set of Euler angles. All need specification
   of the three components { alpha, beta, gamma }, either explicitly or by
   default. Of course, some of the angles may be zero.  Keep in mind that we
   restrict the angle ranges to alpha = gamma = [0,360] & beta = [0,180].
 
        Input Arguments                     Result
       ----------------       ----------------------------------------------
               -              Zero Euler Angles {alpha=0, beta=0, gamma=0}
       alpha,beta,gamma       Euler Angles with these 3 angles
              EA	      New Euler Angles indentical to those in EA    */

EAngles::EAngles ( ) { _alpha=0; _beta=0; _gamma=0; }
EAngles::EAngles(double alpha, double beta, double gamma, bool deg)
  { SetAngles(alpha, beta, gamma, deg); }

EAngles::EAngles(const coord& EA, bool deg)
  { SetAngles(EA.x(), EA.y(), EA.z(), deg); }

EAngles::EAngles(const EAngles& EA)
  {
  _alpha = EA._alpha;			// Copy Euler angle alpha (radians)
  _beta  = EA._beta;			// Copy Euler angle beta  (radians)
  _gamma = EA._gamma;			// Copy Euler angle gamma (radians)
  }

EAngles::~EAngles () { }		// No destruction needed here

EAngles& EAngles::operator= (const EAngles& EA)
  {
  if(this == &EA) return *this;		// Do nothing if already equal
  _alpha = EA._alpha;			// Copy Euler angle alpha (radians)
  _beta  = EA._beta;			// Copy Euler angle beta  (radians)
  _gamma = EA._gamma;			// Copy Euler angle gamma (radians)
  return *this;
  }

// ____________________________________________________________________________
// B                     Euler Angles Access Functions
// ____________________________________________________________________________

/* These allow users to directly access the individual Euler angles. Note 
   that all angles are returned and set assuming units of radians.           */

double EAngles::alpha() const { return _alpha; }
double EAngles::beta()  const { return _beta;  }
double EAngles::gamma() const { return _gamma; }

void EAngles::alpha(double A) { SetAlpha(A); }
void EAngles::beta(double  B) { SetBeta(B);  }
void EAngles::gamma(double G) { SetGamma(G); }

// ____________________________________________________________________________
// C           Class Euler Angles Composite Rotation Functions
// ____________________________________________________________________________

/* These functions handle composite rotations, i.e. generation of a set of
   Euler angles that represents two successive Euler angle rotations.
   Generation of summed rotation angles is done through use of Quaterions.
   Although GAMMA implicitly takes care of any problems, Euler angles define
   a rotation of one coordinate system into another. Quaternions define the
   rotation of data points into another, i.e. the opposite rotation. For this
   reason, the quaternion is generated from the inverse set of Euler angles. */

EAngles EAngles::operator*  (const EAngles& EA) const { return composite(EA); }
EAngles&    EAngles::operator*= (const EAngles& EA)       { *this= composite(EA); return *this; }
EAngles&    EAngles::operator&= (const EAngles& EA)       { *this= EA.composite(*this); return *this;}
EAngles EAngles::composite  (const EAngles& EA) const
  {
  quatern Q1(*this);			// Form Quaterion from us
  quatern Q2(EA);			// Form Quaternion from EA
  quatern Qcmp = Q1*Q2;			// Form Summed Quaternion
  return  Qcmp.EA();			// Return composite Euler Angles
  }

// ____________________________________________________________________________
// D          Class Euler Angles Parameter & Parameter Set Functions
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//                       Single Parameter Functions
//-----------------------------------------------------------------------------

        // Input               EA    : A set of Euler angles (this)
        //                     pname : A parameter name
        // Return              par   : A GAMMA parameter of type coordinate
        //                             with the name pname

SinglePar EAngles::param(const string& pname, bool deg) const
  {
  string  pstate("Euler Angles");		// Parameter statement
  if(deg) pstate += string(" (deg)");		// Adjust name if in degrees
  else    pstate += string(" (radians)");	// or in radians
  return param(pname, pstate, deg);		// Use overload function
  }

SinglePar EAngles::param(const string& pn, const string& ps, bool deg) const
  {
  double  sf = 1.0;				// Scaling for rad. vs. deg.
  if(deg) sf = RAD2DEG;				// Set for angles in degrees
  coord EA(_alpha*sf, _beta*sf, _gamma*sf);	// Use class coord format
  return EA.param(pn,ps); 			// Return parameter via coord
  }

// ----------------------------------------------------------------------------
//            Functions To Make A Parameter Set From Euler Angles
// ----------------------------------------------------------------------------

/* These will 1.) construct a parameter set containing an Euler angle set.
              2.) add an Euler angle set to an existing parameter set.
              3.) write an Euler angle set to an ASCII file in pset format.  */

EAngles::operator ParameterSet( ) const
   { ParameterSet pset; pset += *this; return pset; }

void operator+= (ParameterSet& pset, const EAngles& EA)
  { EA.PSetAdd(pset); }

void EAngles::PSetAdd(ParameterSet& pset, int idx, bool deg) const
  {
  string pname("EAngles");			// Set parameter name
  string pstate("Euler Angles");		// Parameter statement
  double sf = 1.0;				// Adjustment for radians/deg
  if(deg) pstate += string(" (deg)");		// Adjust these if they are
  else 						// to be output in radians
    {						// rather than degrees
    pstate += string(" (radians)");
    pname  += string("Rad");
    }
  string sufx;                                  // Parameter suffix
  if(idx != -1)                                 // Only use suffix if idx
    sufx = string("(")+Gdec(idx)+string(")");	// is NOT -1
  pname += sufx;				// Add suffix to parameter name
  SinglePar par;                                // Temporary single parameter
  if(deg) sf = RAD2DEG;				// Set for angles in degrees
  coord EA(_alpha*sf, _beta*sf, _gamma*sf);	// Use class coord format
  par = EA.param(pname, pstate);		// Parameter for Euler angles
  pset.push_back(par);				// Add parameter to set
  }

void EAngles::write(const string &filename, int idx, bool deg) const
   {
   ofstream ofstr(filename.c_str());		// Open filename for input
   if(!ofstr.good())				// If file bad then exit
     {
     EAerror(1, filename);			// Filename problems
     EAfatal(40);				// Cannot write to file
     }
   ofstr.close();				// Close it now
   ParameterSet pset;				// Declare a parameter set
   PSetAdd(pset,idx,deg);			// Add ourself into param. set
   pset.write(filename);			// Write param. set to filename
   return;
   }

// ----------------------------------------------------------------------------
//           Functions To Make Euler Angles From A Parameter Set
// ----------------------------------------------------------------------------

/* These will 1.) set Euler angles from values in a parameter set.
              2.) set Euler angles from parameters in an ASCII file.         */

bool EAngles::read(const string &filename, int idx, int warn)
   {
   ParameterSet pset;			// Declare a parameter set
   if(!pset.read(filename, warn?1:0))	// Read in pset from file
     {
     EAerror(1, filename,1);		// Filename problems
     if(warn > 1) EAfatal(41);		// Cant read from file, fatal
     else         EAerror(41);		// or non-fatal warning
     return false;
     }
   if(!read(pset,warn?1:0))		// User overloaded function
     {
     EAerror(1, filename,1);		// Filename problems
     if(warn > 1) EAfatal(42);		// Cannot read file parameters
     else         EAerror(42);		// or non-fatal one
     return false;
     }
   return true;
   }

bool EAngles::read(const ParameterSet& pset, int idx, int warn)
  {
  bool TF = SetEAngles(pset,idx,warn?1:0);	// Try & set ourself up
  if(!TF)					// If we didn't handle
    {						// setting ourself properly
    if(warn)					// then we'll issue some
      {						// warnings if desired
                   EAerror(33, 1);		//   Cant read parameters
      if(warn > 1) EAfatal(34);			//   Can't set from parameters
      else         EAerror(34,1);		//   or a warning issued
      }
    return false;				// Return that we failed
    }
  return TF;
  }


// ____________________________________________________________________________
// E                     Class Euler Angles I/O Functions
// ____________________________________________________________________________

// ------------------------ ASCII Output Functions ----------------------------

/*              Input           EA   : Euler Angles (this)
                                ostr : Output ASCII file stream
				deg  : Flag if output angles in degrees
			        hdr  : Flag if header is output
                Return          void : EA is sent to the output stream      */

ostream& EAngles::print(ostream& ostr, bool deg, bool hdr) const
{
  if(hdr)
  {
    if(deg) 
      ostr << "Euler Angles (deg): ";
    else    
      ostr << "Euler Angles (rad): ";
  }
  double sf = deg?RAD2DEG:1.0;
  string fmt("%8.4f");
  if(!deg) fmt = string("%7.5f");
  ostr << "Alpha = " << Gform(fmt.c_str(), _alpha*sf) << "  "
       << "Beta = "  << Gform(fmt.c_str(), _beta*sf)  << "  "
       << "Gamma = " << Gform(fmt.c_str(), _gamma*sf);
  return ostr;
}

ostream& operator << (ostream& ostr, const EAngles& EA)
  { return EA.print(ostr); }

// ____________________________________________________________________________
// F              Class Euler Angle Container Support Functions
// ____________________________________________________________________________

/* These allow the user to make STL lists and vectors from Euler Angles. If
   not present some compilers complain about code such as list<EAngles>

         Function                            Purpose
         ========   ===========================================================
          equal     Compares two sets of Euler angles to see if they produce
                    the same rotation. Equality can occur even if the Euler
                    angles are not the same, so long as the rotations defined
                    are equivalent (produce the same rotation matrix)
                    Note that operator == is NOT the same as equal. The former
                    does an exact comparison of angles, not rotations.
          inverse   Returns a set of Euler angles for the inverse rotation.
            RMx     Returns a 3x3 rotation matrix that will take one coordinate
                    system into another as defined by the Euler angles. The
                    flag inv allows for generation of the inverse rotation
                    matrix.                                                */


void EAngles::SetCutoff(double co)
  { AngCutoff = (co == -1.0)?1.e-10:fabs(co); }

bool EAngles::operator== (const EAngles& EA) const
  {
  if(this == &EA)                          return true;	 // Same set!
  if(fabs(_alpha - EA._alpha) > AngCutoff) return false; // Alpha mismatch
  if(fabs(_beta  - EA._beta)  > AngCutoff) return false; // Beta  mismatch
  if(fabs(_gamma - EA._gamma) > AngCutoff) return false; // Gamma mismatch
  return true;
  }

bool EAngles::operator!= (const EAngles& EA) const
  {
  if(this == &EA)                          return false; // Same set!
  if(fabs(_alpha - EA._alpha) > AngCutoff) return true;	 // Alpha mismatch
  if(fabs(_beta  - EA._beta)  > AngCutoff) return true;  // Beta  mismatch
  if(fabs(_gamma - EA._gamma) > AngCutoff) return true;  // Gamma mismatch
  return false;
  }

bool EAngles::operator<  (const EAngles& EA) const
  {
  if(_alpha < EA._alpha) return true;
  if(_alpha > EA._alpha) return false;
  if(_beta  < EA._beta)  return true;
  if(_beta  > EA._beta)  return false;
  if(_gamma < EA._gamma) return true;
  if(_gamma > EA._gamma) return false;
  return false;
  }

bool EAngles::operator>  (const EAngles& EA) const
  {
  if(_alpha > EA._alpha) return true;
  if(_alpha < EA._alpha) return false;
  if(_beta  > EA._beta)  return true;
  if(_beta  < EA._beta)  return false;
  if(_gamma > EA._gamma) return true;
  if(_gamma < EA._gamma) return false;
  return false;
  }

// ____________________________________________________________________________
// F                Class Euler Angle Auxiliary Functions
// ____________________________________________________________________________

/*       Function                            Purpose
         ========   ===========================================================
          equal     Compares two sets of Euler angles to see if they produce
                    the same rotation. Equality can occur even if the Euler
                    angles are not the same, so long as the rotations defined
                    are equivalent (produce the same rotation matrix)
                    Note that operator == is NOT the same as equal. The former
                    does an exact comparison of angles, not rotations.
          inverse   Returns a set of Euler angles for the inverse rotation.
            RMx     Returns a 3x3 rotation matrix that will take one coordinate
                    system into another as defined by the Euler angles. The
                    flag inv allows for generation of the inverse rotation
                    matrix.                                                  */

bool EAngles::equal(const EAngles& EA, double CUTOFF) const
  {
  if(*this == EA) return true;			// Check if same angles
  return (Rmx()-EA.Rmx()).is_zero(CUTOFF);	// Check rotations
  }

EAngles EAngles::inverse() const { quatern Q(*this, true); return Q.EA(); }

matrix EAngles::RMx(bool inv) const
  {
  matrix mx(3,3);
  double ca = cos(_alpha);
  double sa = sin(_alpha);
  double cb = cos(_beta);
  double sb = sin(_beta);
  double cg = cos(_gamma);
  double sg = sin(_gamma);	
  if(inv)
    {
    mx.put( ca*cb*cg - sa*sg,0,0);//    [                                    ]
    mx.put(-ca*cb*sg - sa*cg,0,1);//    |  CaCbCg-SaSg  -CaCbSg-SaCg    CaSb |
    mx.put( ca*sb,           0,2);//    |                                    |
    mx.put(sa*cb*cg+ca*sg,   1,0);//  -1|                                    |
    mx.put(-sa*cb*sg+ca*cg,  1,1);// R =|  SaCbCg+CaSg  -SaCbSg+CaCg    SaSb |
    mx.put(sa*sb,            1,2);//    |                                    |
    mx.put(-sb*cg,           2,0);//    |                                    |
    mx.put(sb*sg,            2,1);//    |     -SbCg         SbSg         Cb  |
    mx.put(cb,               2,2);//    [                                    ]
    }
  else
    {
    mx.put(ca*cb*cg - sa*sg, 0,0);//    [                                    ]
    mx.put(sa*cb*cg + ca*sg, 0,1);//    |  CaCbCg-SaSg   SaCbCg+CaSg   -SbCg |
    mx.put(-sb*cg,           0,2);//    |                                    |
    mx.put(-(sa*cg+ca*sg*cb),1,0);//    |                                    |
    mx.put(-sa*cb*sg+ca*cg,  1,1);// R =| -CaCbSg-SaCg  -SaCbSg+CaCg    SbSg |
    mx.put(sb*sg,            1,2);//    |                                    |
    mx.put(ca*sb,            2,0);//    |                                    |
    mx.put(sa*sb,            2,1);//    |      CaSb         SaSb         Cb  |
    mx.put(cb,               2,2);//    [                                    ]
    }
  return mx;
  }

matrix EAngles::Rmx()    const { return RMx(); }		// Deprecated
matrix EAngles::invRmx() const { return RMx(true); }		// Deprecated

#endif								// EAngles.cc
