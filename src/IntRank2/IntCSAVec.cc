/* IntCSAVec.cc *************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Shift Anisotropy Interactions Vector     Implementation		**
**                                                                      **
**      Scott Smith                                                     **
**      Copyright (c) 1996                                              **
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
** This class maintains a vector of rank2 SA interactions (from class   **
** IntCSA).  The class allows users to manipulate an array of such      **
** interactions simultaneously.                                         **
**                                                                      **
*************************************************************************/

#ifndef   IntCSAVec_cc_			// Is file already included?
#  define IntCSAVec_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <IntRank2/IntCSAVec.h>		// Include header file
#include <IntRank2/IntCSA.h>		// Include CSA interactions
#include <Basics/Gutils.h>		// Includes query knowledge
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/StringCut.h>		// Include Gdec and Gform
#if defined(_MSC_VER)			// If we are not using MSVC++
 #pragma warning (disable : 4786)       //   Kill STL namelength warnings
#endif

using std::cout;			// Using libstdc++ standard output
using std::string;			// Using libstdc++ strings
using std::vector;			// Using libstdc++ STL vectors
using std::ofstream;			// Using libstdc++ output file streams
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------  
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                 CLASS CSA INTERACTION VECTOR ERROR HANDLING
// ____________________________________________________________________________
 

void IntCSAVec::ICVerror(int eidx, int noret) const

	// Input		ICV	: CSA interaction vector (this)
        // 			eidx	: Error index
        //                      noret	: Flag for return (0=linefeed)
        // Output               none	: Error message

  {
  cout << "\nClass IntCSAVec: ";
  switch(eidx)
    {
    case 0:								// (0)
      cout << "Program Aborting.....";
      break;
    case 1:								// (1)
      cout << "Cannot Access CSA Interaction Between Spins";
      break;
    case 2:								// (2)
      cout << "Error During Vector Construction";
      break;
    case 3:								// (3)
      cout << "Cannot Construct Vector From Paramter Set";
      break;
    case 12:								// (12)
      cout << "Warning - Setting Asymmetry of a Zero Tensor";
      break;
    case 13:								// (13)
      cout << "Attempted CSA Access of Spin with Itself";
      break;
    case 14:								// (14)
      cout << "Sorry, CSA Tensor Operation Not Allowed Yet";
      break;
    case 15:								// (15)
      cout << "Electron Associated With Shift Anisotropy Interaction Disallowed";
      break;
    case 20:								// (20)
      cout << "Can't Write CSA Interaction Vector to Parameter File";
      break;
    case 21:								// (21)
      cout << "Can't Read CSA Interaction Vector from Parameter File";
      break;
    case 22:                                                            // (22)
      cout << "Problems Writing CSA Interaction Vector to Output FileStream";
      break;
    case 23:                                                            // (23)
      cout << "Cannot Output CSA Interaction Vector Parameters";
      break;
    case 24:                                                            // (24)
      cout << "Cannot Access Requested CSA Interaction";
      break;
    case 25:                                                            // (25)
      cout << "Can't Find Any CSA Interactions Amongst Parameters";
      break;
    default: 
      cout << "Unknown error";
      break;
    }
  if(!noret) cout << ".\n";
  }


volatile void IntCSAVec::ICVfatal(int eidx) const
     
	// Input		ICV	: CSA interaction vector (this)
        // 			eidx	: Error index
     	// Output		none	: Error message output
	//				  Program execution stopped
     
  {
  ICVerror(eidx, 1);			// First output the error
  if(eidx) ICVerror(0);			// Now output it fatal
  GAMMAfatal();					// Clean exit from program
  }


void IntCSAVec::ICVerror(int eidx, const string& pname, int noret) const

	// Input		ICV	: CSA interaction vector (this)
        // 			eidx	: Error index
        //                      pname	: Additional string for error
        //                      noret	: Flag for return (0=linefeed)
        // Output               none	: Error message

  {
  cout << "\nClass IntCSAVec: ";
  switch(eidx)
    {
    case 100:							// (100)
      cout << "Can't Read Parameter " << pname;
      break;
    case 101:							// (101)
      cout << "Problems with File " << pname;
      break;
    case 102:							// (102)
      cout << "Construct From Negative # Of Spins (" << pname << ")";
      break;
    case 120:							// (120)
      cout << "CSA Access Of Index " << pname << " Out Of Bounds";
      break;
    case 121:							// (121)
      cout << "Cannot Access CSA Interaction " << pname;
      break;
    case 130:							// (130)
      cout << "Parameter " << pname << " Is The Culprit!\n";
      break;
    default: 
      cout << "Unknown error";
      break;
    }
  if(!noret) cout << ".\n";
  }
 
// ____________________________________________________________________________
// ii               CSA INTERACTIONS VECTOR SETUP FUNCTIONS
// ____________________________________________________________________________
 
/* These functions determine specific aspects of a shift anisotropy 
   interactions vector from a specified parameter set. As the various
   interaction parameters tend to interact, the functions MUST be private
   because their misuse could produce an inconsistent interactions vector.

   The goal here is quite simple. We need to determine a series of shift
   anisotropy interactions from parameters found in an input GAMMA parameter
   set. Each interaction will be defined by unique parameters in the set that
   specify the values: { Iqn, PPM, CSA,eta, alpha, beta, gamma } Complexity
   arises because we need find a series of such interactions and because 
   the interactions allow for a variety of parameters for their definition.

   In this case, since shift anisotropy interactions are singly indexed, there
   is a one to one correspondence between the interaction index and the vector
   index. We simply assume that the first interaction has an index of 0 and
   that the rest are a series spanning 0, 1, 2, .... until the first not
   found.
	   Input		ICV	: CSA interaction vector (this)
                                pset    : A parameter set
                                idx     : Parameter prefix index
                                warn    : Flag for warnings
           Output               none    : Shift anisotropy interactions are
                                          set from parameters in pset

   We Use Only Parameters With Prefix [#] (unles idx=-1). As such we off this
   prefix from all parameter names prior to looking for the interactions.
   This is done before any calls to class IntCSA which doesn't use prefixes. */

bool IntCSAVec::setCIV(const ParameterSet& pset, int idx, bool warn)
  {
  clear();					// Remove any existing intacts.
  ParameterSet subpset;                         // Copy parameter set
  if(idx != -1) subpset = pset.strip(idx);      // to glean out those for
  else          subpset = pset;                 // the specified prefix index
  IntCSA C;					// Empty CSA interaction
  int k = 0; 					//   Start with 0th interaction
  while(C.read(subpset,k, false))		//   If we are able to read it
    {						//   then we keep looking
    push_back(C);				//     Put into the vector
    k++;					//     Index to next interaction
    }
  if(size()) return true;
  if(warn) ICVerror(25, 1);
  return false;
  }

// ____________________________________________________________________________
// iii         SHIFT ANISOTROPY INTERACTIONS VECTOR CHECKING FUNCTIONS
// ____________________________________________________________________________

/* These functions insure vector integrity and that access is not beyond
   the vector boundaries

	   Input		ICV	: CSA interaction vector (this)
	   			spin	: A CSA index
                                warn    : Warning level
           Output               TF      : Returns TRUE if spin is a valid
                                          interaction index, FALSE if not    */

bool IntCSAVec::CheckCI(int spin, int warn) const
  {
  if(spin>=int(size()) || spin<0)
    {
    if(warn)
      {
      ICVerror(120, Gdec(spin), 1);
      ICVerror(121, Gdec(spin), 1);
      }
    if(warn>1) ICVfatal(24);
    return false;
    }
  return true;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// A              CSA INTERACTION VECTOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

/* These functions are inherited from the ANSI standard "vector" class in the
   C++ library libstdc++.  I am listing them here so that I & other users
   don't have to keep looking them up all the time.                          */

IntCSAVec::IntCSAVec() {}
/*
IntCSAVec::IntCSAVec(int N)			Vector w/ N Interactions
IntCSAVec::IntCSAVec(int N, const IntCSA& CI)	Vector w/ N pars
IntCSAVec::IntCSAVec(const IntCSAVec& CSVec)	Vector copy of CSVec
IntCSAVec::IntCSAVec assign(N) 			Assign N element
IntCSAVec::~IntCSAVec()				Destructor of Vector         */

// ----------------------------------------------------------------------------
//                       Constructors Using Parameter Sets 
// ----------------------------------------------------------------------------

/* These functions allow users to fill up the shift anisotropy interactions
   vector from a GAMMA parameter set (external ASCII file). Variants enable
   higher classes such as spin systems to generate the vector when they have
   a knowledge of spin isotope types.

	   Input		ICV	: CSA interaction vector (this)
                                pset    : Parameter set
                                idx	: Index for CSA vector
	  				  i.e [idx] is parameter name prefix
                                warn    : Warning level
                                                0 - no warnings
                                                1 - warn if e- is isotope
                                                2 - fatal if e- is isotope
           Note                         : If idx is negative then no
                                          parameter indicies are used
	   Note				: Since this has no idea about
	  				  the number of interactions it will
	  				  count through CSA(i) until there
	  				  is none found starting with i=0     */
 
IntCSAVec::IntCSAVec(const ParameterSet& pset, int idx, int warn)
  {
  if(!setCIV(pset, idx, warn?true:false))
  if(warn)
    {
    ICVerror(3, 1);
    if(warn>1) ICVfatal(2);
    else       ICVerror(2);
    }
  }

IntCSAVec::IntCSAVec(const vector<Isotope>& Isos, 
                                         const ParameterSet& pset, int warn)
  {
  int ni = Isos.size();			// Number of interactions				
  IntCSA C, C0;				// An empty CSA interaction
  Isotope I;				// A specific spin isotope
  for(int i=0; i<ni; i++)		// Loop proposed interactions
    {
    C = C0;				//   Insure interaction empty
    I = Isos[i];			//   Get spin i isotope
    if(I.electron())			//   If spin i is an electron
      {					//   we cannot allow CSA interaction
      if(warn)				//	If warnings desired 
        {				//	issue them (and maybe fail!)
        ICVerror(15, 1);
        if(warn>1) ICVfatal(2);
        }
      }
    else C = IntCSA(I,pset,i,warn);
    push_back(C);			//   Set ith CSA interaction 
    }
  }

// ----------------------------------------------------------------------------
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------

void IntCSAVec::operator= (const IntCSAVec &ICV) { }
     IntCSAVec::~IntCSAVec() { }

// ____________________________________________________________________________ 
// B                CSA INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//                   Generic CSA Value Access Functions
// ---------------------------------------------------------------------------- 
        
        // Input                ICV     : CSA interaction vector
        //                      spin	: CSA index
        //                      val     : CSA interaction value
        //                                <=0: Isotropic shift (PPM)
	//				    1: CSA of spin (PPM)
        //                                  2: Asymmetry [0, 1] 
	//				    3: Alpha (degrees)
	//				    4: Beta (degrees)
	//				    5: Gamma (degrees)
        //                                  6: Theta (degrees, down from +z)
        //                                  7: Phi   (degrees, over from +z)
        // Output               none    : Get/Set CSA interaction value
        //                                for spicified CSA spin
 
void IntCSAVec::CValue(int spin, double val, int type)

  {
  if(!CheckCI(spin)) ICVfatal(1);		// Check CSA exists
  switch(type)
    {
    default:
    case 0: ((*this)[spin]).PPM(val);   break;	// Isotropic shift   (PPM)
    case 1: ((*this)[spin]).CSA(val);   break;	// Shift aniosotropy (PPM)
    case 2: ((*this)[spin]).eta(val);   break;	// Shfit asymmetry   [0,1]
    case 3: ((*this)[spin]).alpha(val); break;	// Here for alpha
    case 4: ((*this)[spin]).beta(val);  break;	// Here for beta
    case 5: ((*this)[spin]).gamma(val); break;	// Here for gamma
    case 6: ((*this)[spin]).theta(val); break;	// Here for theta
    case 7: ((*this)[spin]).phi(val);   break;	// Here for phi
//    case 4: ((*this)[spin]).delzz(val); break;	// Here for delzz
    }
  }


double IntCSAVec::CValue(int spin, int type) const

  {
  if(!CheckCI(spin)) ICVfatal(1);		// Check CSA exists
  double rval;
  switch(type)
    {
    default:
    case 0: rval = ((*this)[spin]).CSA();   break;	// Here for CSA
    case 1: rval = ((*this)[spin]).eta();   break;	// Here for asymmetry
    case 2: rval = ((*this)[spin]).theta(); break;	// Here for theta
    case 3: rval = ((*this)[spin]).phi();   break;	// Here for phi
//    case 4: rval = ((*this)[spin]).delzz(); break;	// Here for delzz
    }
  return rval;
  }


// ---------------------------------------------------------------------------- 
//                                 CSA Values
// ---------------------------------------------------------------------------- 

	// Input		ICV	: CSA interaction vector
	// 			spin	: CSA index
	// 			dz	: CSA coupling constant (Hertz)
	// Output		none	: Get/Set CSA spin coupling
	// Note				: Defined in class IntCSA as equal to
	//			          the CSA tensor delzz value

void   IntCSAVec::CSA(int spin, double dz)  { CValue(spin, dz, 0); }
double IntCSAVec::CSA(int spin) const       { return CValue(spin, 0); }
void   IntCSAVec::delz(int spin, double dz) { CValue(spin, dz, 0); }
double IntCSAVec::delz(int spin) const      { return CValue(spin, 0); }

// ---------------------------------------------------------------------------- 
//                        CSA Asymmetry Values
// ---------------------------------------------------------------------------- 

	// Input		ICV	: CSA interaction vector
	// 			spin	: CSA index
	// 			deta	: CSA interaction asymmetry
	// Output		none	: Get/Set CSA spin asymmetry
	// Note				: Defined in class IntCSA between [0,1]
	// Note				: Very unusual if nonzero!

void   IntCSAVec::eta(int spin, double ceta) { CValue(spin, ceta, 1); }
double IntCSAVec::eta(int spin) const        { return CValue(spin, 1); }

 
// ---------------------------------------------------------------------------- 
//               CSA Theta Orientation (Down From +z Axis)
// ---------------------------------------------------------------------------- 

	// Input		ICV	: CSA interaction vector
	// 			spin	: CSA index
	// 			dtheta	: CSA interaction angle (deg)
	// Output		none	: Get/Set CSA spin theta angle
	// Note				: Defined class IntCSA between [0,180]

void   IntCSAVec::theta(int spin, double dtheta) { CValue(spin, dtheta, 2); }
double IntCSAVec::theta(int spin) const          { return CValue(spin, 2); }

// ---------------------------------------------------------------------------- 
//               CSA Phi Orientation (Over From +x Axis)
// ---------------------------------------------------------------------------- 

	// Input		ICV	: CSA interaction vector
	// 			spin	: CSA index
	// 			dphi	: CSA interaction angle (deg)
	// Output		none	: Get/Set CSA spin phi angle
	// Note				: Defined in IntCSA between [0,360]

void   IntCSAVec::phi(int spin, double dphi)   { CValue(spin, dphi, 3); }
double IntCSAVec::phi(int spin) const          { return CValue(spin, 3); }
 
// ---------------------------------------------------------------------------- 
//                        Full CSA Interaction
// ---------------------------------------------------------------------------- 


IntCSA& IntCSAVec::operator() (int i) { return (*this)[i]; }
 
	// Input		ICV   : CSA interaction vector (this)
        //                      i     : A CSA interaction index
        // Ouput                DI    : The i'th CSA interaction in ICV
        // Note                       : Returns a reference to the interaction


IntCSA IntCSAVec::get(int spin) const

	// Input		ICV	: CSA interaction vector
	// 			spin	: CSA index
	// Output		DI	: Return rank 2 CSA interaction

  {
  if(!CheckCI(spin)) ICVfatal(1);		// Check CSA exists
  return (*this)[spin];				// Return requested IntCSA
  }
 
// ---------------------------------------------------------------------------- 
//                 Other CSA Interaction Vector Info
// ---------------------------------------------------------------------------- 


// int IntCSAVec::size() const			INHERITED
bool   IntCSAVec::nonzero() const
  {
  for(unsigned i=0; i<size(); i++)
    if((*this)[i].CSA()) return true;
  return false;
  }

// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//     Functions To Make A Parameter Set From A CSA Interaction Vector
// ----------------------------------------------------------------------------

// Note that the preferred means of specifying CSA interactions when there
// is no knowledge of spin isotope types nor spin coordinates are used for
// filling the parameter set.  For base parameters of individual interactions
// see class IntCSA. Any additional parameters for the interaction vector will
// be defined in the += function below (there currently aren't any).  

IntCSAVec::operator ParameterSet( ) const

	// Input		ICV   : CSA interaction vector
        //                      pset  : Parameter set
        // Output               pset  : Parameter set with only
        //                              CSA interaction parameters
 
  { ParameterSet pset; pset += *this; return pset; }		// Add in dynamic system parameters


void operator+= (ParameterSet& pset, const IntCSAVec &ICV)

	// Input		ICV	: CSA interaction vector
	//  			pset	: Parameter set
	// Output		pset	: Parameter set with only
	//			          CSA interaction parameters added

  { ICV.PSetAdd(pset); }

         
void IntCSAVec::PSetAdd(ParameterSet& pset, int idx) const
 
	// Input		ICV	: CSA interaction vector
        //                      pset    : Parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        // Output               void    : Vector parameters are
        //                                are added ot the parameter set
        //                                with interaction index idx
        // Note                         : The parameter names & types
        //                                output here MUST match those used
        //                                in setting the vector up
        //                                from parameters sets

  {
  for(unsigned i=0; i<size(); i++)	// Add each CSA interaction
    (*this)[i].PSetAdd(pset, i, idx);        // to the parameter set
  }

// ----------------------------------------------------------------------------
// Functions To Output CSA Interaction Vector To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

	// Input		ICV   	: CSA interaction vector (this)
	//			filename: Output file name
        //                      ofstr   : Output file stream
        //                      idx	: Parameter index value used for
	//				  prefix [#] in output names
        //                      warn    : Warning level
	// Output		none 	: CSA interaction vector is written
	//				  as a parameter set to file filename
        //                                written as a parameter set to
        //                                or to output filestream ofstr

bool IntCSAVec::write(const string &filename, int idx, int warn) const
  {
  if(!size()) return 1;
  ofstream ofstr(filename.c_str());	// Open filename for input
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!write(ofstr, idx, w2))		// If file bad then exit
    {
    ICVerror(101, filename, 1);		// Filename problems
    if(warn>1) ICVfatal(20);		// Fatal error
    return false;
    }
  ofstr.close();                        // Close it now
  return true;
  }
 
bool IntCSAVec::write(ofstream& ofstr, int idx, int warn) const
  {
  if(!size()) return 1;
  if(!ofstr.good())			// If file bad then exit
    {   
    if(warn) ICVerror(22);		//      Problems with file
    if(warn > 1) ICVfatal(23);		//      It's a fatal error
    return false;
    }   
  ParameterSet pset;			// Declare a parameter set
  for(unsigned i=0; i<size(); i++)	// Add each CSA interaction
    (*this)[i].PSetAdd(pset, i, idx);	// to the parameter set
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!pset.write(ofstr, w2))            // Use parameter set to write
    {
    if(warn)
      {
      ICVerror(22, 1);			// Problems writing to filestream
      if(warn>1) ICVfatal(23);		// Fatal error
      }
    return false;
    }  
  return true;
  }  


// ____________________________________________________________________________ 
// D               CSA INTERACTION VECTOR INPUT FUNCTIONS
// ____________________________________________________________________________ 

// ----------------------------------------------------------------------------
//        Direct Read of Vector From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

	// Input		ICV	: CSA interaction vector (this)
	// 			filename: Input filename
        //                      idx     : Parameter index value used for
	//				  prefix [#] in output names
	// Output		none	: CSA interaction vector filled
	//				  with parameters read from file
	// Input		ICV	: CSA interaction vector (this)
        //                      pset    : Parameter set
        //                      idx     : Parameter index value used for
        //                                prefix [#] in input names
	// Output		none	: CSA interaction vector filled
	//				  with parameters read from pset


bool IntCSAVec::read(const string &filename, int idx, int warn)
  {
  ParameterSet pset;                    // Declare a parameter set
  if(!pset.read(filename, warn?1:0))    // Read in pset from file
     {
     if(warn)
       {
       ICVerror(101, filename, 1);      //      Problems with file
       if(warn > 1) ICVfatal(21);       //      Its a fatal error
       else         ICVerror(21);       //      or a warning issued
       }
     return false;
     }
  return read(pset, idx, warn);
  }

bool IntCSAVec::read(const ParameterSet& pset, int idx, int warn)
  {
  bool TF;
  TF = setCIV(pset, idx, warn?true:false);
  if(!TF)
    {
    if(warn)
      {
      if(warn > 1) ICVfatal(21);        //      Its a fatal error
      else         ICVerror(21);        //      or a warning issued
      }
    return false;
    }
  return TF;
  }


// ----------------------------------------------------------------------------
//       Interactive Read of CSA Interaction Vector From An ASCII File
// ----------------------------------------------------------------------------


string IntCSAVec::ask_read(int argc, char* argv[], int argn)

	// Input		ICV    : CSA interaction vector (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
        // Output               void    : The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the spin system is read
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : The spin system is modifed (filled)

  {
  string filename;				// Name of spin system file  
  query_parameter(argc, argv, argn,		// Get filename from command
   "\n\tCSA Interaction Vector Filename? ", 
                                      filename);
  read(filename);		           	// Read system from filename
  return filename;
  }

// ____________________________________________________________________________ 
// E                 CSA INTERACTION VECTOR OUTPUT FUNCTIONS
// ____________________________________________________________________________ 

	// Input		ICV	: CSA interaction vector (this)
	// 			ostr	: Output stream
	//			full    : Flag for long vs short output
	// Output		non	: CSA interaction vector parameters
	//			          sent to the output stream

ostream& IntCSAVec::print(ostream& ostr, int full) const
  {
  string hdr="CSA Interactions Vector";		// Title for output
  string ctr;					// Spacer to center output 
  int ns = size();				// Get number of interactions
  if(!ns) 					// Exit if no interactions
    { 
    hdr += string(": Empty");
    ctr=string(40-hdr.length()/2, ' ');		// Set spacer to center output 
    ostr << "\n" << ctr << hdr;			// Output the header
    return ostr;
    }
  ctr = string(40-hdr.length()/2, ' ');		// Set spacer to center output 
  ostr << "\n" << ctr << hdr;			// Output the header
  if(!nonzero())				// If all interactions zero
    {						// just out the number we have
    hdr = Gdec(ns)				//	Make a new header
        + string(" Zero Interactions");
    ctr=string(40-hdr.length()/2, ' ');		// 	Spacer to center output 
    ostr << "\n\n" << ctr << hdr;		//	Output no nonzero ones
    return ostr;				//	Exit
    }
  cout << "\n";
  for(int j=0; j<ns; j++)			// Output CSA interactions
    {						// in order
    hdr = string("Interaction ") + Gdec(j);	//	Make a new header
    ctr=string(40-hdr.length()/2, ' ');		// 	Spacer to center output 
    ostr << "\n" << ctr << hdr;			//	Output a header
    ((*this)[j]).print(ostr, full);			//	Output the interaction
    ostr << "\n";				// 	Add a line spacer
    }
  return ostr;
  }

ostream& operator<< (ostream& out, const IntCSAVec& ICV)
  { return ICV.print(out); }


#endif
