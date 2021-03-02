/* IntQuadVec.cc ************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Quadrupolar Interactions Vector     Implementation		**
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
** This class maintains a vector of rank 2 Quadrupolar interactions (as **
** defined in class IntQuad). The class lets users, or more importantly **
** spin systems, to manipulate an array of such interactions either in	**
** concert or individually.						**
**                                                                      **
*************************************************************************/

#ifndef   IntQuadVec_cc_		// Is file already included?
#  define IntQuadVec_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#if defined(_MSC_VER)			// If we are using MSVC++
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif

#include <IntRank2/IntQuadVec.h>	// Include header file
#include <IntRank2/IntQuad.h>		// Include Quadrupolar interactions
#include <Basics/Gutils.h>		// Include query parameter
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/Isotope.h>		// Include spin isotopes
#include <Basics/StringCut.h>		// Include Gform and Gdec
#include <stdlib.h>

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ STL lists
using std::vector;			// Using libstdc++ STL vectors
using std::ofstream;			// Using libstdc++ output file streams
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------  
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                  QUADRUPOLAR INTERACTION VECTOR ERROR HANDLING
// ____________________________________________________________________________
 
/*      Input                   IQV     : Quadrupolar interaction vector (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : Added error message for output
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void IntQuadVec::IQVerror(int eidx, int noret) const
  {
  string hdr("Class IntQuadVec");
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


volatile void IntQuadVec::IQVfatality(int eidx) const
  {
  IQVerror(eidx, 1);			// First output the error
  if(eidx) IQVerror(0);			// Now output it fatal
  GAMMAfatal();					// Clean exit from program
  }

/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Function pname
      5     Please Use Class Member Funciton pname
   default  Unknown Error - pname                                           */

void IntQuadVec::IQVerror(int eidx, const string& pname, int noret) const
  {
  string hdr("Class IntQuadVec");
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
      msg = string("Cannot Access Quadrupolar Interaction ")
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
// ii              QUADRUPOLAR INTERACTION VECTOR SETUP FUNCTIONS
// ____________________________________________________________________________
 
	// Input		IQV	: Quadrupolar interaction vector (this)
	// 			spin	: A G index
	//			warn    : Warning level
	// Output		TF	: Returns TRUE if spin is a valid
	//			  	  interaction index, FALSE if not

bool IntQuadVec::check_spin(int spin, int warn) const
  {
  if(spin>=int(size()) || spin<0)
    {
    if(warn)
      {
      IQVerror(120, Gdec(spin), 1);
      IQVerror(121, Gdec(spin), 1);
      }
    if(warn>1) IQVfatality(24);
    return false;
    }
  return true;
  }

// ---------------------------------------------------------------------------- 
//    Functions To Set Quadrupolar Interactions From Single Index Parameters
// ---------------------------------------------------------------------------- 

int IntQuadVec::getNInts(const ParameterSet& pset, int idx) const

	// Input		IQV	: Quadrupolar interaction vector (this)
        //                      pset	: A parameter set
	//			idx     : Parameter prefix index
        // Output               nd	: The number of Quadrupolar interactions
        //                                referenced by single indices in pset
	// Note				: This function looks for G(#) where
	//				  # = { 0, 1, ...., nd-1 }.

  {
  string pnameb("gxx(");			// Parameter base name
  string pnamee(")");				// Parameter ending
  string pfx;					// Parameter name prefix
  string sfx;					// Parameter name suffix
  if(idx >= 0)					// Set prefix if desired 
    {
    pfx = string("[")+Gdec(idx)+string("]");
    pnameb = pfx + pnameb;
    }
  int nd = 0;					// Number of Gs
  string pname = pnameb + Gdec(nd) + pnamee;	// Parameter name this G
  ParameterSet::const_iterator item;		// A pix into parameter list
  item = pset.seek(pname);			// Pix in pset for parameter
  while(item!=pset.end() && nd<500)		// See how many are in pset
    {
    nd ++;					//	Increment G count
    pname = pnameb + Gdec(nd) + pnamee;		// 	Adjust name next G
    item = pset.seek(pname);			// 	Pix pset for parameter
    }
  return nd;
  }


int IntQuadVec::setGs(const ParameterSet& pset, int idx, int warn)

	// Input		IQV	: Quadrupolar interaction vector (this)
	// 			pset	: A parameter set
	//			idx     : Parameter prefix index
	//			warn    : Warning flag
	//				      0 = no warnings
	//				     !0 = warnings
	// Output		none	: System Quadrupolar interactions are
	//			  	  set from parameters in pset
	// Note				: Assumed that the number of Gs
	//				  is ALREADY set and space allocated
	// Note				: This uses only parameters that
	//				  have the prefix "[idx]" if idx!=-1
 
/* This function uses the Quadrupolar interaction constructors directly.
   In turn, they use a single interaction (or spin) index.  They'll attempt
   to read, for each interaction in the vector, a set of parameters

                   { [Iso(i)/CI], G, Ctheta, Cphi, Ceta }

   where the shielding anisotropy is both read in and maintained in PPM units.
   G is related to the spatial tensor delzz via
 
      ^          3                                     1 
     / \ sigma = - del   = sigma  - sigma   = sigma  - - [ sigma  + sigma  ]
     ---         2    zz        ||        |        zz  2        xx       yy
                                         ---
   An Example Quadrupolar Interaction Definition Are The ASCII Lines

   CI(1)        (1) : 1.5                  - Spin I quantum number        
   G(1)         (1) : 40.0                 - Shift anisotropy (PPM)
   Gphi(1)      (1) : 45.0                 - Tensor orientation (deg)
   Gtheta(1)    (1) : 45.0                 - Tensor orientation (deg)
   Ceta(1)      (1) : 0.0                  - Tensor asymmetry [0,1]          */

//   Use Only Parameters With Prefix [#], Clip Off Prefix From Names
//   Since We'll Use IntQuad Class Directly Which Doesn't Use Prefixes
 
  {
  ParameterSet subpset;				// Copy parameter set
  if(idx != -1) subpset = pset.strip(idx);	// to glean out those for
  else          subpset = pset;			// the specified prefix index
  for(int i=0; i<int(size()); i++)		// Loop interactions and set 'em
    (*this)[i] = IntQuad(subpset,i,1);		// using the base class
  return 1;
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// A               QUADRUPOLAR INTERACTION VECTOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

/* These functions are inherited from the ANSI standard "vector" class in the 
   C++ library libstdc++.  I am listing them here so that I & other users don't
   have to keep looking them up all the time.
*/

   IntQuadVec::IntQuadVec() {};			// Empty Interaction Vector
/*
   IntQuadVec(int N)				Vector w/ N Interactions
   IntQuadVec(int N, const IntQuad& GI)		Vector w/ N pars
   IntQuadVec(const IntQuadVec& GVec)			Vector copy of GVec
   IntQuadVec assign(N)				Assign N element
   ~IntQuadVec()					Destructor of Vector         */

/*
IntQuadVec::IntQuadVec(const ParameterSet& pset, int indx, int warn)

	// Input		IQV	: Quadrupolar interaction vector (this)
        //                      pset    : Parameter set
        //                      idx	: Index for G vector
	//				  i.e [idx] is parameter name prefix
        //                      warn    : Warning level
        // Note                         : If idx is negative then no
        //                                parameter indicies are used
	// Note				: Since this has no idea about
	//				  the number of interactions it will
	//				  count through G(i) until there
	//				  is none found  starting with i=0
 
  {
//	      Attempt To Determine The Number of G Interactions

  int ns = getNInts(pset, indx);	// Get # of Quadrupolar interactions
  if(ns < 1)				// If we can't find any then
    {					// we'll have to abort the mission
    IQVerror(3,1);			// 	Attempt with - # spins
    IQVerror(3);			//	Cant construct from pset
    IQVfatality(2);			//	Error during constructions
    }
cout << "\n\tFound " << ns << " G Interactions...";
cout.flush();
//  *this = IntQuadVec(ns);		// Set up vector for ns Gs

//	   Attempt To Read { CI, G, Ctheta, Cphi, Ceta }

  setGs(pset, indx);			// Now set the Gs
  return;
  }


IntQuadVec::IntQuadVec(int ns, const ParameterSet& pset, int indx, int warn)

	// Input		IQV	: Quadrupolar interaction vector (this)
	//			ns	: Number of spins in vector
        //                      pset    : Parameter set
        //                      idx	: Index for G vector
	//				  i.e [idx] is parameter name prefix
        //                      warn    : Warning level
        // Note                         : If idx is negative then no
        //                                parameter indicies are used
 
  {
  if(ns<1)				// If number of spins negative
    {					// we can't construct a vector
    IQVerror(102,Gdec(ns),1);
    IQVfatality(2);
    }
//  *this = IntQuadVec(ns);		// Set up vector for ns Gs

//	   Attempt To Read { CI, G, Ctheta, Cphi, Ceta }

  setGs(pset, indx);		// Now set the Gs
  }
*/


// ----------------------------------------------------------------------------
//          This Constructor Supports Generation From A Spin System
// ----------------------------------------------------------------------------

/* In this case the vector of interactions is associated with a list of isotope
   types (Isos). The interaction vector length will match the vector of spin
   types. For each electron spin type a G interaction will be generated from
   parameters in the parameter set.
 
	   Input		IQV	: Quadrupolar interaction vector (this)
                                Isos    : Array of isotope/spin types
                                pset    : Parameter set
                                warn    : Warning level
                                                0 - no warnings
                                                1 - warn if spin is nucleon
	                                        2 - fatal if spin is nucleon
           Output               none    : Quadrupolar interaction vector
                                          constructed with spin labels
                                          in Isos and parameters in pset
           Note                         : We don't allow users to associate
                                          a nucleon with a G interaction. For
                                          each nucleus a NULL interaction
                                          is stored in the vector.           */
 
IntQuadVec::IntQuadVec(const vector<Isotope>& Isos,
                                           const ParameterSet& pset, int warn)
  {
  int ni = Isos.size();			// # possible interactions
  for(int i=0; i<ni; i++)		// Loop proposed interactions
    {
    if(Isos[i].electron())		//      Insure not an electron
      {					//	if so, don't make G
      if(warn)				//	interaction, but warn if 
        {				//	desired (or fail!)
        IQVerror(15, 1);
        if(warn>1) IQVfatality(2);
        }
      }
// sosik
//    push_back(IntQuad(Isos[i],pset,i,warn));// Nucleus, try & construct
    }
  }

// ----------------------------------------------------------------------------
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------


//void IntQuadVec::operator= (const IntQuadVec &IQV)

	// Input		IQV	: Quadrupolar interaction vector (this)
	// Output		none	: Quadrupolar interaction vector is
	//				  constructed equivalent to sys

//  {
//  if(nspins) delete [] Cs;		// Delete any spin. interactions
//  nspins = IQV.nspins;			// Set number of Gs
//  if(nspins)
//    {
//    Cs = new IntQuad[nspins];		// New array of spin. interactions
//    for(int i=0;i<nspins;i++) 		// Copy spin. interactions
//      Cs[i] = IQV.Cs[i];
//    }
//  else Cs = NULL;
//  }


//IntQuadVec::~IntQuadVec () {}

	// Input		IQV	: Quadrupolar interaction vector (this)
	// Output		none	: System IQV is destructed

//  {
//  if(Cs) delete [] Cs;			// Destroy any Quadrupolar interactions
//  Cs = NULL;				// Insure this is indeed nothing
//  } 					// Rest destructs itself


// ____________________________________________________________________________ 
// B                 QUADRUPOLAR INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//                   Generic G Value Access Functions
// ---------------------------------------------------------------------------- 
        
        // Input                IQV     : Quadrupolar interaction vector
        //                      spin	: G index
        //                      val     : Quadrupolar interaction value
        //                                <=0: G of spin (PPM)
        //                                  1: Asymmetry [0, 1] 
        //                                  2: Theta (degrees, down from +z)
        //                                  3: Phi   (degrees, over from +z)
	//				    4: The delzz value (PPM)
        // Output               none    : Get/Set Quadrupolar interaction value
        //                                for spicified G spin
 
void IntQuadVec::QValue(int spin, double val, int type)

  {
  if(!check_spin(spin)) IQVfatality(1);		// Check G exists
  switch(type)
    {
    default:
    case 0: (*this)[spin].QCC(val);  break;	// Here for qcc
    case 1: (*this)[spin].eta(val);   break;	// Here for asymmetry
    case 2: (*this)[spin].alpha(val); break;	// Here for alpha
    case 3: (*this)[spin].beta(val); break;	// Here for beta
    case 4: (*this)[spin].gamma(val); break;	// Here for gamma
    case 5: (*this)[spin].theta(val); break;	// Here for theta
    case 6: (*this)[spin].phi(val);   break;	// Here for phi
    }
  }

double IntQuadVec::QValue(int spin, int type) const

  {
  if(!check_spin(spin)) IQVfatality(1);		// Check G exists
  double rval;
  switch(type)
    {
    default:
    case 0: rval = (*this)[spin].QCC();   break;	// Here for qcc
    case 1: rval = (*this)[spin].eta();   break;	// Here for asymmetry
    case 2: rval = (*this)[spin].alpha(); break;	// Here for alpha
    case 3: rval = (*this)[spin].beta();  break;	// Here for beta
    case 4: rval = (*this)[spin].gamma(); break;	// Here for gamma
    case 5: rval = (*this)[spin].theta(); break;	// Here for theta
    case 6: rval = (*this)[spin].phi();   break;	// Here for phi
    }
  return rval;
  }

// ----------------------------------------------------------------------------
//                         Quadrupolar Coupling Constants
// ----------------------------------------------------------------------------

        // Input                IQV     : Quadrupolar interaction vector
        //                      spin    : Quadrupolar index
        //                      qcc     : Quadrupolar coupling constant (Hz)
        // Output               none    : Get/Set spin quadrupolar coupling
        // Note                         : QCC is identical to delzz for this
        //                                interaction and the anisotropy is
        //                                1.5 times the Quadrupolar coupling

void   IntQuadVec::QCC(int  spin, double qcc) { QValue(spin, qcc,      0); }
void   IntQuadVec::NQCC(int spin, double qcc) { QValue(spin, qcc,      0); }
void   IntQuadVec::delz(int spin, double qcc) { QValue(spin, qcc,      0); }
void   IntQuadVec::QA(int   spin, double qa)  { QValue(spin, qa*2./3., 0); }

double IntQuadVec::QCC(int  spin) const       { return QValue(spin, 0);     }
double IntQuadVec::NQCC(int spin) const       { return QValue(spin, 0);     }
double IntQuadVec::delz(int spin) const       { return QValue(spin, 0);     }
double IntQuadVec::QA(int   spin) const       { return QValue(spin, 0)*1.5; }

// ---------------------------------------------------------------------------- 
//                          Quadrupolar Asymmetry Values
// ---------------------------------------------------------------------------- 

	// Input		IQV	: Quadrupolar interaction vector
	// 			spin	: Interaction or spin index
	// 			qeta	: Quadrupolar interaction asymmetry
	// Output		none	: Get/Set Quadrupolar asymmetry
	// Note				: Defined to be between [0,1]

void   IntQuadVec::eta(int spin, double ceta) { QValue(spin, ceta, 1); }
double IntQuadVec::eta(int spin) const        { return QValue(spin, 1); }

//-----------------------------------------------------------------------------
//                     Quadrupolar Orientation Angle Access
//-----------------------------------------------------------------------------

/* There are 3 angles which orient the interaction spatial tensor relative to
   the interactions own principal axis system (PAS).  These are the set of
   Euler angles { alpha, beta, gamma }. They corrspond to rotating the
   coordinate axes first about the z-axis by alpha followed by a rotation about
   the new x-axis by beta and lastly about the new z-axis by gamma. For a
   symmetric tensor the last rotation is of no consequnce and set to zero, and
   for powder averages on the first two angles are used to sum over all
   spatial orientations. In these two cases the angles alpha and beta are one
   and the same as the spherical coordinate angles phi and theta respectively.
   Theta is the angle down from the PAS z-axis & the angle phi which is over
   from the PAS x-axis. The ranges of the angles (internally imposed on
   GAMMA's Euler angles) are alpha,gamma,phi = [0,360] and beta,theta=[0,180].
   These functions allow users to both obtain and set these angles for any
   intraction. Setting any of the angles will effectively reorient the
   spatial tensor (interaction).                                             */

double  IntQuadVec::alpha(int spin)       const { return QValue(spin, 2); }
double  IntQuadVec::beta(int spin)        const { return QValue(spin, 3); }
double  IntQuadVec::gamma(int spin)       const { return QValue(spin, 4); }
double  IntQuadVec::theta(int spin)       const { return QValue(spin, 5); }
double  IntQuadVec::phi(int spin)         const { return QValue(spin, 6); }
//EAngles IntQuadVec::orientation(int spin) const { return QValue(spin, 7); }

void IntQuadVec::alpha(int spin,double A) { QValue(spin, A , 2); }
void IntQuadVec::beta(int  spin,double B) { QValue(spin, B , 3); }
void IntQuadVec::gamma(int spin,double G) { QValue(spin, G , 4); }
void IntQuadVec::theta(int spin,double T) { QValue(spin, T , 5); }
void IntQuadVec::phi(int   spin,double P) { QValue(spin, P , 6); }

//void IntQuadVec::orientation(int spin, const EAngles& EA);
//void IntQuadVec::orientation(int spin,
//                            double A, double B, double G, bool deg=false);   */









// ---------------------------------------------------------------------------- 
//                        Full Quadrupolar Interaction
// ---------------------------------------------------------------------------- 

	// Input		IQV   : Quadrupolar interaction vector (this)
        //                      i     : A Quadrupolar interaction index
        // Ouput                DI    : The i'th Quadrupolar interaction in IQV
        // Note                       : Returns a reference to the interaction


      IntQuad& IntQuadVec::operator() (int i) { return (*this)[i]; }
const IntQuad& IntQuadVec::getcref(int     i) const { return (*this)[i]; }
      IntQuad  IntQuadVec::get(int         i) const
  {
  if(!check_spin(i)) IQVfatality(1);		// Check G exists
  return (*this)[i];				// Return requested IntQuad
  }
 
// ---------------------------------------------------------------------------- 
//                 Other Quadrupolar Interaction Vector Info
// ---------------------------------------------------------------------------- 


//int IntQuadVec::size() const { return nspins; }

	// Input		IQV   : Quadrupolar interaction vector
	// Output		nd    : Number of interactions in vector


int IntQuadVec::nonzero() const

	// Input		IQV   : Quadrupolar interaction vector
	// Output		TF    : True if any interactions with a
	//				finite G value

   {
   for(int i=0; i<int(size()); i++)
     if((*this)[i].qxx()) return 1;
   return 0;
   }
// sosi - need better test above

// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//    Functions To Make A Parameter Set From A Quadrupolar Interaction Vector
// ----------------------------------------------------------------------------

/* Note that the preferred means of specifying a single Quadrupolar interaction
   when there is no knowledge of the spin isotope type is used for filling the
   parameter set.  For base parameters of individual interactions see the class
   IntQuad. Any additional parameters for the interaction vector will be defined
   in the += function below (there currently aren't any).  Also note that we
   provide a means of writing the interaction(s) when the isotope type is known
   as this will be the case when a spin system controls the output.  In this
   latter case the parameters used are slightly different.

	   Input		IQV	: Quadrupolar interaction vector
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

IntQuadVec::operator ParameterSet( ) const
  { ParameterSet pset; pset += *this; return pset; }

void operator+= (ParameterSet& pset, const IntQuadVec &IQV)
  { IQV.PSetAdd(pset); }

void IntQuadVec::PSetAdd(ParameterSet& pset, int idx) const
  {
  for(int i=0; i<int(size()); i++)
    (*this)[i].PSetAdd(pset, i, idx);
  }


// ---------------------------------------------------------------------------- 
//     Functions To Make A Quadrupolar Interaction Vector From A Parameter Set
// ---------------------------------------------------------------------------- 
	// Input		IQV	: Quadrupolar interaction vector (this)
	// 			pset	: A parameter set
	// Output		none	: Quadrupolar interaction vector filled
	//				  with parameters in pset
	// Note                         : We just use a member function
	//                                which is needed to allow for
	//                                [#] prefixed parameter names

        // Input                IQV     : Quadrupolar interaction vector (this)
        //                      pset    : A parameter set
	//                      idx	: Parameter index value used for
	//                                prefix [#] in input names
	//                      warn	: Warning output level
	//                                      0 = no warnings
	//                                      1 = warnings
	//                                     >1 = fatal warnings
        // Output               TF	: Vector of interactions is filled
        //                                with parameters in pset
        // Note                         : This uses the assignment from pset
        //                                It exists to support prefix indices


void IntQuadVec::operator= (const ParameterSet& pset) { setIQVec(pset); }
bool IntQuadVec::setIQVec(const ParameterSet& pset, int idx, int warn)
  {
//      Use Only Parameters With Prefix [#], So Clip [#] From Names First

  ParameterSet  subpset;                        // Working parameter set
  if(idx != -1) subpset=pset.strip(idx);        // Get only params with [#]
  else          subpset=pset;                   // Or use full pset

//  First Attempt To Determine The Number Of G Interactions In Parameter Set

  bool TF=true;				// Track if we read OK
  int nd = getNInts(pset, idx);		// Get # of Quadrupolar interactions
  if(nd<1)				// If we can't find any then
    {					// we'll have to abort the mission
    IQVerror(3,1);
    IQVfatality(2);
    }

//  Now Attempt To Read nd Quadrupolar Interactions From Values In pset
//  Each IntQuad As { gxx,gyy,gzz,Gtheta,Gphi } or { Giso,Gdelz,Geta,Gtheta,Gphi }

  clear();				// Remove any existing interactions
//  for(int i=0; i<nd; i++)		// Loop interactions and set 'em
//    push_back(IntQuad(subpset,i,1));	// using the base class
  IntQuad G;
  for(int i=0; i<nd; i++)		// Loop interactions and set 'em
    {
    if(G.read(subpset,i))
      push_back(G);
    }
  return TF;
  }


// ----------------------------------------------------------------------------
//     Output Quadrupolar Interaction Vector To ASCII From A Parmeter Set
// ----------------------------------------------------------------------------

int IntQuadVec::write(const string &filename, int idx, int warn) const

	// Input		IQV   	: Quadrupolar interaction vector (this)
	//			filename: Output file name
        //                      idx	: Parameter index value used for
	//				  prefix [#] in output names
        //                      warn    : Warning level
	// Output		none 	: Quadrupolar interaction vector is written
	//				  as a parameter set to file filename

  {
  if(!size()) return 1;
  ofstream ofstr(filename.c_str());	// Open filename for input
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!write(ofstr, idx, w2))		// If file bad then exit
    {
    IQVerror(1, filename, 1);		// Filename problems
    if(warn>1) IQVfatality(20);		// Fatal error
    return 0;
    }
  ofstr.close();                        // Close it now
  return 1;
  }


int IntQuadVec::write(ofstream& ofstr, int idx, int warn) const

	// Input		IQV   	: Quadrupolar interaction vector (this)
        //                      ofstr   : Output file stream
        //                      idx     : Parameter index value used for
	//				  prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Quadrupolar interaction vector is
        //                                written as a parameter set to
        //                                output filestream
 
  {
  if(!size()) return 1;
  if(!ofstr.good())			// If file bad then exit
    {   
    if(warn) IQVerror(22);		//      Problems with file
    if(warn > 1) IQVfatality(23);	//      It's a fatal error
    return 0;
    }   
  ParameterSet pset;			// Declare a parameter set
  for(int i=0; i<int(size()); i++)	// Add each Quadrupolar interaction
    ((*this)[i]).PSetAdd(pset, i, idx);	// to the parameter set
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!pset.write(ofstr, w2))            // Use parameter set to write
    {
    if(warn)
      {
      IQVerror(22, 1);			// Problems writing to filestream
      if(warn>1) IQVfatality(23);	// Fatal error
      }
    return 0;
    }  
  return 1;
  }  


// ____________________________________________________________________________ 
// D                QUADRUPOLAR INTERACTION VECTOR INPUT FUNCTIONS
// ____________________________________________________________________________ 

// ----------------------------------------------------------------------------
//        Direct Read of Vector From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

/* These functions allow users to fill the vector of interactions either from
   parameters found in an external ASCII file or contained in a GAMMA parameter
   set.
	   Input		IQV	: Quadrupolar interaction vector (this)
	   			filename: Input filename
                             or pset    : Parameter set
                                idx     : Parameter index value used for
	  				  prefix [#] in parameter names
                                warn    : Warning output level
                                               0 = no warnings
                                               1 = warnings
                                              >1 = fatal warnings
	   Output		none	: Quadrupolar interaction vector filled
	  				  with parameters read from file     */

bool IntQuadVec::read(const string &filename, int idx, int warn)
  {
  ParameterSet pset;                    // Declare a parameter set
  if(!pset.read(filename, warn?1:0))    // Try and read in pset
    {                                   // If we cannot read the file then
    if(warn)                            // we'll issue warnings as desired
      {
      IQVerror(1, filename, 1); 	//      Problems with file
      if(warn>1) IQVfatality(21);	//      This is a fatal problem
      else       IQVerror(21, 1);	//      Or maybe it aint so bad
      }
    return false;
    }
  return read(pset, idx, warn);
  }

bool IntQuadVec::read(const ParameterSet& pset, int idx, int warn) 
  {
  bool TF = setIQVec(pset, idx, warn?1:0);	// Use overload to read
  if(!TF)                                       // If setIGvec didn't handle
    {                                           // the vector read from pset
    if(warn)                                    // then we'll issue errors
      {
      IQVerror(8, 1);				//    Problems with pset
      if(warn>1) IQVfatality(21);		//    This is a fatal problem
      if(warn>1) IQVerror(21);			//    Or maybe it isn't so bad..
      }
    return false;
    }
  return TF;
  }


// ----------------------------------------------------------------------------
//       Interactive Read of Quadrupolar Interaction Vector From An ASCII File
// ----------------------------------------------------------------------------


string IntQuadVec::ask_read(int argc, char* argv[], int argn)

	// Input		IQV    : Quadrupolar interaction vector (this)
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
   "\n\tQuadrupolar Interaction Vector Filename? ", 
                                      filename);
  read(filename);		           	// Read system from filename
  return filename;
  }

// ____________________________________________________________________________ 
// E                  QUADRUPOLAR INTERACTION VECTOR OUTPUT FUNCTIONS
// ____________________________________________________________________________ 

/* These functions provied ASCII output of the interaction vector into any
   output stream, including C++ standard output.

	   Input		IQV	: Quadrupolar interaction vector (this)
	   			ostr	: Output stream
	  			full    : Flag for long vs short output
	   Output		non	: Quadrupolar interaction vector
					  parameters sent to output stream  */


ostream& IntQuadVec::print(ostream& ostr, int full) const
  {
  string hdr="Quadrupolar Interactions Vector";	// Title for output
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
  for(int j=0; j<int(size()); j++)		// Output Quadrupolar interactions
    {						// in order
    hdr = string("Interaction ") + Gdec(j);	//	Make a new header
    ctr=string(40-hdr.length()/2, ' ');		// 	Spacer to center output 
    ostr << "\n" << ctr << hdr;			//	Output a header
    ((*this)[j]).print(ostr, full);		//	Output the interaction
    ostr << "\n";				// 	Add a line spacer
    }
  return ostr;
  }


ostream& operator<< (ostream& out, const IntQuadVec& IQV)
  { return IQV.print(out); }


#endif
