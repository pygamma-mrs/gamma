/* IntGVec.cc ***************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Electron G Interactions Vector     Implementation		**
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
** This class maintains a vector of rank 2 electron G interactions (as  **
** defined in class IntG).  The class allows users, or more importantly **
** spin systems, to manipulate an array of such interactions either in	**
** concert or individually.						**
**                                                                      **
*************************************************************************/

#ifndef   IntGVec_cc_			// Is file already included?
#  define IntGVec_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#if defined(_MSC_VER)			// If using MSVC++ then we
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif 

#include <IntRank2/IntGVec.h>		// Include header file
#include <IntRank2/IntG.h>		// Include Electron G interactions
#include <Basics/Gutils.h>		// Include query parameter
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/Isotope.h>		// Include spin isotopes
#include <Basics/StringCut.h>		// Include Gdec and Gform

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ STL lists
using std::vector;			// Using libstdc++ STL vectors
using std::ofstream;			// Using libstdc++ output file streams
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------  
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                 ELECTRON G INTERACTION VECTOR ERROR HANDLING
// ____________________________________________________________________________
 
/*      Input                   IGV     : Electron G interaction vector (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : Added error message for output
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void IntGVec::IGVerror(int eidx, int noret) const
  {
  string hdr("Class IntGVec");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr,"Program Aborting.....",        noret); break;// (0)
    case 1: GAMMAerror(hdr,"Can't Access Spin G Interact.",noret); break;// (1)
    case 2: GAMMAerror(hdr,"Problems During Construction", noret); break;// (2)
    case 3: GAMMAerror(hdr,"Can't Construct From PSet",    noret); break;// (3)
    case 12:GAMMAerror(hdr,"Set Asymmetry Of Zero Interac",noret); break;// (12)
    case 14:GAMMAerror(hdr,"Sorry, Operation Inactive..." ,noret); break;// (14)
    case 15:GAMMAerror(hdr,"Nucleon G Interact Disallowed",noret); break;// (15)
    case 20:GAMMAerror(hdr,"Can't Write To Parameter File",noret); break;// (20)
    case 21:GAMMAerror(hdr,"Cant Read From Parameter File",noret); break;// (21)
    case 22:GAMMAerror(hdr,"Can't Write To Output Stream", noret); break;// (22)
    case 23:GAMMAerror(hdr,"Can't Ouptut Vector Params.",  noret); break;// (23)
    case 24:GAMMAerror(hdr,"Can't Access Requested Cmpnt.",noret); break;// (24)
    default:GAMMAerror(hdr,"Unknown Error",                noret); break;
    }
  }


volatile void IntGVec::IGVfatal(int eidx) const
  {
  IGVerror(eidx, 1);			// First output the error
  if(eidx) IGVerror(0);			// Now output it fatal
  GAMMAfatal();					// Clean exit from program
  }

/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Function pname
      5     Please Use Class Member Funciton pname
   default  Unknown Error - pname                                           */

void IntGVec::IGVerror(int eidx, const string& pname, int noret) const
  {
  string hdr("Class IntGVec");
  string msg;
  switch(eidx)
    {
    case 102:                                                   // (102)
      msg = string("Construct From Negative # Of Spins (")
          + pname + string(")");
      GAMMAerror(hdr, msg, noret); break;
    case 120:                                                   // (120)
      msg = string("Interaction Access Of Index ")
          + pname + string(" Out Of Bounds");
      GAMMAerror(hdr, msg, noret); break;
    case 121:                                                   // (121)
      msg = string("Cannot Access Electron G Interaction ")
          + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 130:                                                   // (130)
      msg = string("Parameter ")
          + pname + string(" Is The Culprit!");
      GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;
    }
  }
 
// ____________________________________________________________________________
// iii            ELECTRON G INTERACTIONS VECTOR SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions determine specific aspects of an electron g interactions
   vector from a specified parameter set. As the various interaction parameters
   tend to interact, the functions MUST be private because their misuse could
   produce an inconsistent interactions vector.

   The goal here is quite simple. We need to determine a series of electron g
   interactions from parameters found in an input GAMMA parameter set. Each 
   interaction will be defined by unique parameters in the set that specify the
   values: { Iqn, g, gA, eta, alpha, beta, gamma }. Complexity arises because
   we need find a series of such interactions and because the interactions 
   allow for a wide variety of parameters for their definition.

   In this case, since electron G interactions are singly indexed, there is a
   one to one correspondence between the interaction index and the vector
   index. We simply assume that the first interaction has an index of 0 and
   that the rest are a series spanning 0, 1, 2, .... until the first not
   found.
           Input                IGV     : G interactions vector (this)
                                pset    : A parameter set
                                idx     : Parameter prefix index
                                warn    : Flag for warnings
           Output               none    : Electron G interactions are
                                          set from parameters in pset

   We Use Only Parameters With Prefix [#] (unles idx=-1). As such we off this
   prefix from all parameter names prior to looking for the interactions.
   This is done before any calls to class IntG which doesn't use prefixes.  */

bool IntGVec::setGIV(const ParameterSet& pset, int idx, bool warn)
  {
  clear();                                      // Remove any existing intacts.
  ParameterSet subpset;                         // Copy parameter set
  if(idx != -1) subpset = pset.strip(idx);      // to glean out those for
  else          subpset = pset;                 // the specified prefix index
  IntG G;					// Empty G interaction
  int k = 0;                                    //   Start with 0th interaction
  while(G.read(subpset, k, false))		//   If we are able to read it
    {                                           //   then we keep looking
    push_back(G);                               //     Put into the vector
    k++;                                        //     Index to next interaction
    }
  if(size()) return true;
  if(warn) IGVerror(25, 1);
  return false;
  }

// ____________________________________________________________________________
// iii            ELECTRON G INTERACTIONS VECTOR CHECKING FUNCTIONS
// ____________________________________________________________________________
 
/* These functions insure vector integrity and that access is not beyond
   the vector boundaries

	   Input		IGV	: Electron G interaction vector (this)
	   			spin	: A G index
	  			warn    : Warning level
	   Output		TF	: Returns TRUE if spin is a valid
	  			  	  interaction index, FALSE if not    */

bool IntGVec::check_spin(int spin, int warn) const
  {
  if(spin>=int(size()) || spin<0)
    {
    if(warn)
      {
      IGVerror(120, Gdec(spin), 1);
      IGVerror(121, Gdec(spin), 1);
      }
    if(warn>1) IGVfatal(24);
    return false;
    }
  return true;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// A              ELECTRON G INTERACTION VECTOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 

// ----------------------------------------------------------------------------
//                     Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

/* These functions are inherited from the ANSI standard "vector" class in the 
   C++ library libstdc++.  I am listing them here so that I & other users don't
   have to keep looking them up all the time.
*/

   IntGVec::IntGVec() {}			// Empty Interaction Vector
/*
   IntGVec(int N)				Vector w/ N Interactions
   IntGVec(int N, const IntG& GI)		Vector w/ N pars
   IntGVec(const IntGVec& GVec)			Vector copy of GVec
   IntGVec assign(N)				Assign N element
   ~IntGVec()					Destructor of Vector         */

// ----------------------------------------------------------------------------
//                       Constructors Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions allow users to fill up the electron G interactions vector
   from a GAMMA parameter set (external ASCII file). Variants enable higher
   classes such as spin systems to generate the vector when they have a 
   knowledge of spin isotope types.

	   Input		IGV	: Electron G interaction vector (this)
                                pset    : Parameter set
                                idx     : Index for G vector
                                          i.e [idx] is parameter name prefix
                                warn    : Warning level
                                                0 - no warnings
                                                1 - warn if spin isn't e- 
                                                2 - fatal if spin isn't e- 
           Note                         : If idx is negative then no
                                          parameter indicies are used
           Note                         : Since this has no idea about
                                          the number of interactions it will
                                          count through G(i) until there
                                          is none found starting with i=0     */


IntGVec::IntGVec(const ParameterSet& pset, int idx, int warn)
  {
  if(!setGIV(pset, idx, warn?true:false))
  if(warn)
    {
    IGVerror(3, 1);
    if(warn>1) IGVfatal(2);
    else       IGVerror(2);
    }
  }

IntGVec::IntGVec(const vector<Isotope>& Isos,
                                         const ParameterSet& pset, int warn)
  {
  int ni = Isos.size();                 // Number of interactions
  IntG G, G0;				// An empty G interaction
  Isotope I;                            // A specific spin isotope
  for(int i=0; i<ni; i++)               // Loop proposed interactions
    {
    G = G0;                             //   Insure interaction empty
    I = Isos[i];                        //   Get spin i isotope
    if(!I.electron())			//   If spin i is not an electron
      {                                 //   we cannot allow G interaction
      if(warn)                          //      If warnings desired
        {                               //      issue them (and maybe fail!)
        IGVerror(15, 1);
        if(warn>1) IGVfatal(2);
        }
      }
    else G = IntG(I,pset,i,warn);
    push_back(G);                       //   Set ith G interaction
    }
  }

// ----------------------------------------------------------------------------
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------


//void IntGVec::operator= (const IntGVec &IGV)

	// Input		IGV	: Electron G interaction vector (this)
	// Output		none	: Electron G interaction vector is
	//				  constructed equivalent to sys

//  {
//  if(nspins) delete [] Cs;		// Delete any spin. interactions
//  nspins = IGV.nspins;			// Set number of Gs
//  if(nspins)
//    {
//    Cs = new IntG[nspins];		// New array of spin. interactions
//    for(int i=0;i<nspins;i++) 		// Copy spin. interactions
//      Cs[i] = IGV.Cs[i];
//    }
//  else Cs = NULL;
//  }


//IntGVec::~IntGVec () {}

	// Input		IGV	: Electron G interaction vector (this)
	// Output		none	: System IGV is destructed

//  {
//  if(Cs) delete [] Cs;			// Destroy any Electron G interactions
//  Cs = NULL;				// Insure this is indeed nothing
//  } 					// Rest destructs itself


// ____________________________________________________________________________ 
// B                ELECTRON G INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//                   Generic G Value Access Functions
// ---------------------------------------------------------------------------- 
        
        // Input                IGV     : Electron G interaction vector
        //                      spin	: G index
        //                      val     : Electron G interaction value
        //                                <=0: G of spin (PPM)
        //                                  1: Asymmetry [0, 1] 
        //                                  2: Theta (degrees, down from +z)
        //                                  3: Phi   (degrees, over from +z)
	//				    4: The delzz value (PPM)
        // Output               none    : Get/Set Electron G interaction value
        //                                for spicified G spin
 
void IntGVec::CValue(int spin, double val, int type)

  {
  if(!check_spin(spin)) IGVfatal(1);		// Check G exists
  switch(type)
    {
    default:
//    case 0: (*this)[spin].G(val);     break;	// Here for G
    case 1: (*this)[spin].eta(val);   break;	// Here for asymmetry
    case 2: (*this)[spin].theta(val); break;	// Here for theta
    case 3: (*this)[spin].phi(val);   break;	// Here for phi
//    case 4: (*this)[spin].delzz(val);  break;	// Here for delzz
    }
  }


double IntGVec::CValue(int spin, int type) const

  {
  if(!check_spin(spin)) IGVfatal(1);		// Check G exists
  double rval;
  switch(type)
    {
    default:
//    case 0: rval = (*this)[spin].G();     break;	// Here for G
    case 1: rval = (*this)[spin].eta();   break;	// Here for asymmetry
    case 2: rval = (*this)[spin].theta(); break;	// Here for theta
    case 3: rval = (*this)[spin].phi();   break;	// Here for phi
    case 4: rval = (*this)[spin].delzz();  break;	// Here for delzz
    }
  return rval;
  }


// ---------------------------------------------------------------------------- 
//                                 G Values
// ---------------------------------------------------------------------------- 

	// Input		IGV	: Electron G interaction vector
	// 			spin	: G index
	// 			dz	: G coupling constant (Hertz)
	// Output		none	: Get/Set G spin coupling
	// Note				: Defined in class IntG as equal to
	//			          the G tensor delzz value

//void   IntGVec::G(int spin, double dz)  { CValue(spin, dz, 0); }
//double IntGVec::G(int spin) const       { return CValue(spin, 0); }
void   IntGVec::delz(int spin, double dz) { CValue(spin, dz, 0); }
double IntGVec::delz(int spin) const      { return CValue(spin, 0); }

// ---------------------------------------------------------------------------- 
//                        G Asymmetry Values
// ---------------------------------------------------------------------------- 

	// Input		IGV	: Electron G interaction vector
	// 			spin	: G index
	// 			deta	: Electron G interaction asymmetry
	// Output		none	: Get/Set G spin asymmetry
	// Note				: Defined in class IntG between [0,1]
	// Note				: Very unusual if nonzero!

void   IntGVec::eta(int spin, double ceta) { CValue(spin, ceta, 1); }
double IntGVec::eta(int spin) const        { return CValue(spin, 1); }

 
// ---------------------------------------------------------------------------- 
//               G Theta Orientation (Down From +z Axis)
// ---------------------------------------------------------------------------- 

	// Input		IGV	: Electron G interaction vector
	// 			spin	: G index
	// 			dtheta	: Electron G interaction angle (deg)
	// Output		none	: Get/Set G spin theta angle
	// Note				: Defined class IntG between [0,180]

void   IntGVec::theta(int spin, double dtheta) { CValue(spin, dtheta, 2); }
double IntGVec::theta(int spin) const          { return CValue(spin, 2); }

// ---------------------------------------------------------------------------- 
//               G Phi Orientation (Over From +x Axis)
// ---------------------------------------------------------------------------- 

	// Input		IGV	: Electron G interaction vector
	// 			spin	: G index
	// 			dphi	: Electron G interaction angle (deg)
	// Output		none	: Get/Set G spin phi angle
	// Note				: Defined in IntG between [0,360]

void   IntGVec::phi(int spin, double dphi)   { CValue(spin, dphi, 3); }
double IntGVec::phi(int spin) const          { return CValue(spin, 3); }
 
// ---------------------------------------------------------------------------- 
//                        Full Electron G Interaction
// ---------------------------------------------------------------------------- 


IntG& IntGVec::operator() (int i) { return (*this)[i]; }
 
	// Input		IGV   : Electron G interaction vector (this)
        //                      i     : A Electron G interaction index
        // Ouput                DI    : The i'th Electron G interaction in IGV
        // Note                       : Returns a reference to the interaction


IntG IntGVec::get(int spin) const

	// Input		IGV	: Electron G interaction vector
	// 			spin	: G index
	// Output		DI	: Return rank 2 Electron G interaction

  {
  if(!check_spin(spin)) IGVfatal(1);		// Check G exists
  return (*this)[spin];				// Return requested IntG
  }
 
// ---------------------------------------------------------------------------- 
//                 Other Electron G Interaction Vector Info
// ---------------------------------------------------------------------------- 


//int IntGVec::size() const { return nspins; }

	// Input		IGV   : Electron G interaction vector
	// Output		nd    : Number of interactions in vector


int IntGVec::nonzero() const

	// Input		IGV   : Electron G interaction vector
	// Output		TF    : True if any interactions with a
	//				finite G value

   {
   for(int i=0; i<int(size()); i++)
     if((*this)[i].gxx()) return 1;
   return 0;
   }
// sosi - need better test above

// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//    Functions To Make A Parameter Set From A Electron G Interaction Vector
// ----------------------------------------------------------------------------

/* Note that the preferred means of specifying a single electron G interaction
   when there is no knowledge of the spin isotope type is used for filling the
   parameter set.  For base parameters of individual interactions see the class
   IntG. Any additional parameters for the interaction vector will be defined
   in the += function below (there currently aren't any).  Also note that we
   provide a means of writing the interaction(s) when the isotope type is known
   as this will be the case when a spin system controls the output.  In this
   latter case the parameters used are slightly different.

	   Input		IGV	: Electron G interaction vector
                                pset    : Parameter set
                                idx     : Parameter index value used for
                                          prefix [#] in output names
           Output               void    : Vector parameters are
                                          are added ot the parameter set
                                          with interaction index idx
           Note                         : The parameter names & types
                                          output here MUST match those used
                                          in setting the vector up
                                          from parameters sets               */

IntGVec::operator ParameterSet( ) const
  { ParameterSet pset; pset += *this; return pset; }

void operator+= (ParameterSet& pset, const IntGVec &IGV)
  { IGV.PSetAdd(pset); }

void IntGVec::PSetAdd(ParameterSet& pset, int idx) const
  {
  for(int i=0; i<int(size()); i++)
    (*this)[i].PSetAdd(pset, i, idx);
  }

// ----------------------------------------------------------------------------
//     Output Electron G Interaction Vector To ASCII From A Parmeter Set
// ----------------------------------------------------------------------------

int IntGVec::write(const string &filename, int idx, int warn) const

	// Input		IGV   	: Electron G interaction vector (this)
	//			filename: Output file name
        //                      idx	: Parameter index value used for
	//				  prefix [#] in output names
        //                      warn    : Warning level
	// Output		none 	: Electron G interaction vector is written
	//				  as a parameter set to file filename

  {
  if(!size()) return 1;
  ofstream ofstr(filename.c_str());	// Open filename for input
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!write(ofstr, idx, w2))		// If file bad then exit
    {
    IGVerror(1, filename, 1);		// Filename problems
    if(warn>1) IGVfatal(20);		// Fatal error
    return 0;
    }
  ofstr.close();                        // Close it now
  return 1;
  }


int IntGVec::write(ofstream& ofstr, int idx, int warn) const

	// Input		IGV   	: Electron G interaction vector (this)
        //                      ofstr   : Output file stream
        //                      idx     : Parameter index value used for
	//				  prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Electron G interaction vector is
        //                                written as a parameter set to
        //                                output filestream
 
  {
  if(!size()) return 1;
  if(!ofstr.good())			// If file bad then exit
    {   
    if(warn) IGVerror(22);		//      Problems with file
    if(warn > 1) IGVfatal(23);	//      It's a fatal error
    return 0;
    }   
  ParameterSet pset;			// Declare a parameter set
  for(int i=0; i<int(size()); i++)	// Add each Electron G interaction
    ((*this)[i]).PSetAdd(pset, i, idx);	// to the parameter set
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!pset.write(ofstr, w2))            // Use parameter set to write
    {
    if(warn)
      {
      IGVerror(22, 1);			// Problems writing to filestream
      if(warn>1) IGVfatal(23);	// Fatal error
      }
    return 0;
    }  
  return 1;
  }  


// ____________________________________________________________________________ 
// D               ELECTRON G INTERACTION VECTOR INPUT FUNCTIONS
// ____________________________________________________________________________ 

// ----------------------------------------------------------------------------
//        Direct Read of Vector From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

/* These functions allow users to fill the vector of interactions either from
   parameters found in an external ASCII file or contained in a GAMMA parameter
   set.
	   Input		IGV	: Electron G interaction vector (this)
	   			filename: Input filename
                             or pset    : Parameter set
                                idx     : Parameter index value used for
	  				  prefix [#] in parameter names
                                warn    : Warning output level
                                               0 = no warnings
                                               1 = warnings
                                              >1 = fatal warnings
	   Output		none	: Electron G interaction vector filled
	  				  with parameters read from file     */

bool IntGVec::read(const string &filename, int idx, int warn)
  {
  ParameterSet pset;                    // Declare a parameter set
  if(!pset.read(filename, warn?1:0))    // Try and read in pset
    {                                   // If we cannot read the file then
    if(warn)                            // we'll issue warnings as desired
      {
      IGVerror(1, filename, 1); 	//      Problems with file
      if(warn>1) IGVfatal(21);		//      This is a fatal problem
      else       IGVerror(21, 1);	//      Or maybe it aint so bad
      }
    return false;
    }
  return read(pset, idx, warn);
  }

bool IntGVec::read(const ParameterSet& pset, int idx, int warn) 
  {
  bool TF = setGIV(pset, idx, warn?1:0);	// Use overload to read
  if(!TF)                                       // If setIGvec didn't handle
    {                                           // the vector read from pset
    if(warn)                                    // then we'll issue errors
      {
      IGVerror(8, 1);				//    Problems with pset
      if(warn>1) IGVfatal(21);			//    This is a fatal problem
      if(warn>1) IGVerror(21);			//    Or maybe it isn't so bad..
      }
    return false;
    }
  return TF;
  }


// ----------------------------------------------------------------------------
//       Interactive Read of Electron G Interaction Vector From An ASCII File
// ----------------------------------------------------------------------------


string IntGVec::ask_read(int argc, char* argv[], int argn)

	// Input		IGV    : Electron G interaction vector (this)
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
   "\n\tElectron G Interaction Vector Filename? ", 
                                      filename);
  read(filename);		           	// Read system from filename
  return filename;
  }

// ____________________________________________________________________________ 
// E                 ELECTRON G INTERACTION VECTOR OUTPUT FUNCTIONS
// ____________________________________________________________________________ 

/* These functions provied ASCII output of the interaction vector into any
   output stream, including C++ standard output.

	   Input		IGV	: Electron G interaction vector (this)
	   			ostr	: Output stream
	  			full    : Flag for long vs short output
	   Output		non	: Electron G interaction vector
					  parameters sent to output stream  */


ostream& IntGVec::print(ostream& ostr, int full) const
  {
  string hdr="Electron G Interactions Vector";	// Title for output
  string ctr;					// Spacer to center output 
  if(!size()) 					// Exit if no interactions
    { 
    hdr += string(": NULL");
    ctr=string(40-hdr.length()/2, ' ');		// Set spacer to center output 
    ostr << "\n" << ctr << hdr;			// Output the header
    return ostr;
    }
  ctr = string(40-hdr.length()/2, ' ');		// Set spacer to center output 
  ostr << "\n" << ctr << hdr;			// Output the header
  if(!nonzero())				// If all interactions zero
    {						// just out the number we have
    hdr = Gdec(size())				//	Make a new header
        + string(" Zero Interactions");
    ctr=string(40-hdr.length()/2, ' ');		// 	Spacer to center output 
    ostr << "\n\n" << ctr << hdr;		//	Output no nonzero ones
    return ostr;				//	Exit
    }
  for(int j=0; j<int(size()); j++)		// Output Electron G interactions
    {						// in order
    hdr = string("Interaction ") + Gdec(j);	//	Make a new header
    ctr=string(40-hdr.length()/2, ' ');		// 	Spacer to center output 
    ostr << "\n" << ctr << hdr;			//	Output a header
    ((*this)[j]).print(ostr, full?true:false);		//	Output the interaction
    ostr << "\n";				// 	Add a line spacer
    }
  return ostr;
  }


ostream& operator<< (ostream& out, const IntGVec& IGV)
  { return IGV.print(out); }


#endif
