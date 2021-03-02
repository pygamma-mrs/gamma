/* IntHFVec.cc **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Hyperfine Interactions Vector     Implementation		**
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
** This class maintains a vector of rank 2 electron-nucleon hyperfine   **
** interactions (as defined in class IntHF). The class allows users, or **
** more importantly spin systems, to manipulate an array of such        **
** interactions either in concert or individually.                      **
**                                                                      **
*************************************************************************/

#ifndef   IntHFVec_cc_			// Is file already included?
#  define IntHFVec_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#if defined(_MSC_VER)			// If using MSVC++ then we
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif 

#include <IntRank2/IntHFVec.h>		// Include header file
#include <IntRank2/IntHF.h>		// Include Hyperfine interactions
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
// i                 HYPERFINE INTERACTION VECTOR ERROR HANDLING
// ____________________________________________________________________________
 
/*      Input                   IHFV     : Hyperfine interaction vector (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : Added error message for output
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void IntHFVec::IHFVerror(int eidx, int noret) const
  {
  string hdr("Class IntHFVec");
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
    case 30:GAMMAerror(hdr,"Can't Find Number of Spins",   noret); break;// (30)
    default:GAMMAerror(hdr,"Unknown Error",                noret); break;
    }
  }


volatile void IntHFVec::IHFVfatality(int eidx) const
  {
  IHFVerror(eidx, 1);				// First output the error
  if(eidx) IHFVerror(0);			// Now output it fatal
  GAMMAfatal();					// Clean exit from program
  }

/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Function pname
      5     Please Use Class Member Funciton pname
   default  Unknown Error - pname                                           */

void IntHFVec::IHFVerror(int eidx, const string& pname, int noret) const
  {
  string hdr("Class IntHFVec");
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
      msg = string("Cannot Access Hyperfine Interaction ")
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
// ii             HYPERFINE INTERACTION VECTOR SETUP FUNCTIONS
// ____________________________________________________________________________
 
	// Input		IHFV	: Hyperfine interaction vector (this)
	// 			spin	: A G index
	//			warn    : Warning level
	// Output		TF	: Returns TRUE if spin is a valid
	//			  	  interaction index, FALSE if not

bool IntHFVec::check_spin(int spin, int warn) const
  {
  if(spin>=int(size()) || spin<0)
    {
    if(warn)
      {
      IHFVerror(120, Gdec(spin), 1);
      IHFVerror(121, Gdec(spin), 1);
      }
    if(warn>1) IHFVfatality(24);
    return false;
    }
  return true;
  }

// ---------------------------------------------------------------------------- 
//     Functions To Set Hyperfine Interactions From Single Index Parameters
// ---------------------------------------------------------------------------- 

/* Hyperfine interactions may be specified using a single interaction index
   (as opposed to spin pair indices). The functions in this section attempt to
   read a vector of N such interactions that are index sequentially from [0,N).
   The value of N, if not set, weill be determined by the first interaction 
   that i s not fully specified in the parameter set.

   The function NInts tries to determine the maximum vector length, i.e. the
   number of interactions that should be created. This is done by looking for

	// Input		IHFV	: Hyperfine interaction vector (this)
        //                      pset	: A parameter set
	//			idx     : Parameter prefix index, [#]
                                warn    : Warning level
                                            0 - no warnings
                                            1 - warn if no spins specified
                                            2 - fatal if no spins specified
        // Output               void	: Interaction vector filled from
        //                                pset single indexed parameters    */

void IntHFVec::setHFs(const ParameterSet& pset, int idx, int warn)
  {
  clear();					// Remove any current
  ParameterSet subpset;				// Copy parameter set
  if(idx != -1) subpset = pset.strip(idx);	// to glean out those for
  else          subpset = pset;			// the specified prefix index
  IntHF HFtmp;					// Working HF interaction
  int i=0;					// Interaction index
  while(HFtmp.read(subpset, i++, 0))		// Fill vector until we cannot
    push_back(HFtmp);				// read in HF interaction i
  }


// ----------------------------------------------------------------------------
//   Functions To Set Hyperfine Interactions From Two Spin Indexed Parameters
// ----------------------------------------------------------------------------

/* Hyperfine interactions can be specified using two spin indices, i.e. spin
   pair indices (as opposed to using a single interaction index). These
   functions attempt to read all such interactions spanning N spins.  The
   vector vector will contain a hyperfine interaction for each unique spin
   pair, although it may be NULL if the interaction is not present in the
   parameter set.

	   Input		IHFV	: Hyperfine interaction vector (this)
                                pset    : A parameter set
                                indx    : Parameter prefix index ([#])
                                warn    : Warning level
                                            0 - no warnings
                                            1 - warn if no spins specified
                                            2 - fatal if no spins specified
           Output               ns      : # of spins specified in pset       */

int IntHFVec::getNSpins(const ParameterSet& pset, int indx, int warn) const
  {
  string prefix, pname;
  if(indx != -1)
    prefix = string("[")
           + Gdec(indx) + string("]");
  pname = prefix + string("NSpins");
  SinglePar par(pname);
  ParameterSet::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);                      // Try and get parameter
  if(item != pset.end())                        // Retrieve the number of spins
    {
    int npts;
    string pstate;
    (*item).parse(pname,npts,pstate);
    return npts;
    }
  else if(warn)
    {
    IHFVerror(30, 1);				// Can't find number of spins
    IHFVerror(130, pname, 1);			// Parameter pname is problem
    if(warn>1) IHFVfatality(3);			// Problems during construction
    else       IHFVerror(3);
    }
  return 0;
  }

void IntHFVec::setHFs(int ns, const ParameterSet& pset, int idx, int warn)
  {
  clear();					// Remove any current
  ParameterSet subpset;				// Copy parameter set
  if(idx != -1) subpset = pset.strip(idx);	// to glean out those for
  else          subpset = pset;			// the specified prefix index
  IntHF HFtmp;					// Working HF interaction
  int i, j;					// Spin indices
  for(i=0; i<ns-1; i++)				// Loop over all possible spin
    for(j=i+1; j<ns; j++)			// pairs, allow for NULL 
      push_back(IntHF(subpset, i, j, 0));	// interactions in this case
  }

// ----------------------------------------------------------------------------
//  Functions To Set Complete Hyperfine Interaction Vector From Parameter Set
// ----------------------------------------------------------------------------

/* Hyperfine interactions can be specified using either two spin indices or a
   single interaction index. For a vector of hyperfine interactions to be 
   generated from a parameter set we must first determine how many interactions
   the vector will contain. Once that is done we can loop over the all the
   interactions and have each read itself from the parameter set (class IntHF) 

           Input                IHFV    : Hyperfine interaction vector (this)
                                pset    : A parameter set
                                indx    : Parameter prefix index ([#])
                                warn    : Warning level
                                            0 - no warnings
                                            1 - warn if no spins specified
                                            2 - fatal if no spins specified
           Output               TF      : True if successulf in setting
                                          the vector from parameters in pset */

bool IntHFVec::setIHFVec(const ParameterSet& pset, int idx, int warn)
  {
  clear();				// Remove any existing interactions

//     First, Attempt To Read Hyperfine Interactions Using 2 Spin Indices

  int ns = getNSpins(pset, idx, 0);     // Try and read # of spins
//  int i, j, k;				// Spin indices
  int i, j;				// Spin indices
  if(ns)                                // If # spins was set, try &
    {                                   // read all HF interactions
    for(i=0; i<ns-1; i++)		//   Loop over the spin pairs
      for(j=i+1; j<ns; j++)
        push_back(IntHF(i,j,pset,0));
    return true;
    }

//      Next, Attempt To Read Hyperfine Interactions Using Single Indices

/*
  int nd = GetNInts(pset, idx);         // Get # of dipole interactions
  if(nd<1)                              // If we can't find any then
    {                                   // we'll have to abort the mission
    IDVerror(3,1);
    IDVfatality(2);
    }
  *this = IntDipVec(nd);                // Construct vector with nd dipoles
  setDs(pset, idx);                     // Now set the dipoles
*/
  return false;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// A              HYPERFINE INTERACTION VECTOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

/* These functions are inherited from the ANSI standard "vector" class in the 
   C++ library libstdc++.  I am listing them here so that I & other users don't
   have to keep looking them up all the time.
*/

IntHFVec::IntHFVec() { }
IntHFVec::IntHFVec(const IntHFVec& HFVec) { *this=HFVec; }

/*

	// Input		IHFV	: Hyperfine interaction vector (this)
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
 
IntHFVec::IntHFVec(const ParameterSet& pset, int indx, int warn)
  {
//	      Attempt To Determine The Number of G Interactions

  int ns = getNInts(pset, indx);	// Get # of Hyperfine interactions
  if(ns < 1)				// If we can't find any then
    {					// we'll have to abort the mission
    IHFVerror(3,1);			// 	Attempt with - # spins
    IHFVerror(3);			//	Cant construct from pset
    IHFVfatality(2);			//	Error during constructions
    }
cout << "\n\tFound " << ns << " G Interactions...";
cout.flush();
//  *this = IntHFVec(ns);		// Set up vector for ns Gs

//	   Attempt To Read { CI, G, Ctheta, Cphi, Ceta }

  setGs(pset, indx);			// Now set the Gs
  return;
  }


IntHFVec::IntHFVec(int ns, const ParameterSet& pset, int indx, int warn)

	// Input		IHFV	: Hyperfine interaction vector (this)
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
    IHFVerror(102,Gdec(ns),1);
    IHFVfatality(2);
    }
//  *this = IntHFVec(ns);		// Set up vector for ns Gs

//	   Attempt To Read { CI, G, Ctheta, Cphi, Ceta }

  setGs(pset, indx);		// Now set the Gs
  }
*/


// ----------------------------------------------------------------------------
//          This Constructor Supports Generation From A Spin System
// ----------------------------------------------------------------------------

/* In this case the vector of interactions is associated with a list of isotope
   types (Isos). The interaction vector length will be set appropriate to the
   # of spins in the vector.  This will be the number of unique spin pairs
   available.  For each nucleon-electron spin pair a hyperfine interaction will
   be generated from parameters in the parameter set.
 
	   Input		IHFV	: Hyperfine interaction vector (this)
                                Isos    : Array of isotope/spin types
                                pset    : Parameter set
                                warn    : Warning level
                                                0 - no warnings
                                                1 - warn if spin is nucleon
	                                        2 - fatal if spin is nucleon
           Output               none    : Hyperfine interaction vector
                                          constructed with spin labels
                                          in Isos and parameters in pset
           Note                         : We don't allow users to associate a
                                          nucleon-nucleon or electron-electron
                                          spin pairing with a HF interaction.
                                          For each improper spin pairing the
                                          vector contains a NULL interaction */

// ----------------------------------------------------------------------------
//                    Constructors Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions allow users to fill up hyperfine interactions vector from
   a GAMMA parameter set (external ASCII file). Variants enable higher classes
   such as spin systems to generate the vector when they have a knowledge of
   spin isotope types and hyperfine interaction parameters (parameter set).

	   Input		IHFV	: Hyperfine interaction vector (this)
        //                      pset    : Parameter set
        //                      idx     : Index for hyperfine vector
        //                                i.e [idx] is parameter name prefix
        //                      warn    : Warning level (applicable to use
        //                                of isotopes & parameter set)
        // Note                         : If idx is negative then no
        //                                parameter prefix will be used     */

 
//IntHFVec::IntHFVec(const vector<Isotope>& Is,const ParameterSet& pset,int wrn)
//  {
//  int i, j;				// Spin indices, interaction index
//int i;
//  int ns = Is.size();			// No. of spin isotopes
//  int ni = 0;				// No. HF interactions possible
//  for(i=1; i<ns; i++) ni+=i;		// Set # possible interactions
// sosik
/*
  for(i=0; i<ns-1; i++)			// Loop proposed interactions
    for(j=i+1; j<ns; j++)		//      (loop spin pairs)
      push_back(IntHF(Is[i],Is[j],pset,i,j,wrn));
*/
//  }


/*
IntHFVec::IntHFVec(const ParameterSet& pset, int idx, int warn)
  {
  if(!setDIV(pset, idx, warn?true:false))
    if(warn)
      {
      IDVerror(3, 1);
      if(warn>1) IDVfatal(2);
      else       IDVerror(2);
      }
  }
*/

IntHFVec::IntHFVec(const ParameterSet& pset,
                              const vector<Isotope>& Isos, int idx, int warn)
  {
  int ns = Isos.size();                         // Number of spins
  if(!ns) return;                               // Bail if no spins
  ParameterSet subpset;                         // Copy parameter set
  if(idx != -1) subpset = pset.strip(idx);      // to glean out those for
  else          subpset = pset;                 // the specified prefix index
  Isotope I,S;                                  // Spin isotopes
  int i, j;                                     // Individual spin indices
  IntHF H, H0;					// Working (NULL) interactions
  for(i=0; i<ns-1; i++)                         // Set hyperfine interactions
    for(j=i+1; j<ns; j++)                       // per each unique spin pair
      {
      H = H0;                                   //   Insure empty interaction
      if(I.nepair(S))				//   If e-/nucleus pair we read
        H.read(subpset,i,j,0);                  //   in the interaction, else 0
      push_back(H);                             //   Put interaction in vector
      }
  }

// ----------------------------------------------------------------------------
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------

/* Both assignment and destruction are handled by the base vector class.     */

void IntHFVec::operator= (const IntHFVec& HFV) { vector<IntHF>::operator=(HFV); }
     IntHFVec::~IntHFVec()                     { }

// ____________________________________________________________________________ 
// B                HYPERFINE INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//             Generic Hyperfine Interaction Value Access Functions
// ---------------------------------------------------------------------------- 
        
        // Input                IHFV	: Hyperfine interaction vector
        //                      spin	: G index
        //                      val     : Hyperfine interaction value
        //                                <=0: G of spin (PPM)
        //                                  1: Asymmetry [0, 1] 
        //                                  2: Theta (degrees, down from +z)
        //                                  3: Phi   (degrees, over from +z)
	//				    4: The delzz value (PPM)
        // Output               none    : Get/Set Hyperfine interaction value
        //                                for spicified G spin
 
void IntHFVec::CValue(int spin, double val, int type)

  {
  if(!check_spin(spin)) IHFVfatality(1);		// Check exists
  switch(type)
    {
    default:
    case 0: (*this)[spin].A(val);     break;	// Here for A
    case 1: (*this)[spin].eta(val);   break;	// Here for asymmetry
    case 2: (*this)[spin].theta(val); break;	// Here for theta
    case 3: (*this)[spin].phi(val);   break;	// Here for phi
//    case 4: (*this)[spin].delz(val);  break;	// Here for delzz
    }
  }


double IntHFVec::CValue(int spin, int type) const

  {
  if(!check_spin(spin)) IHFVfatality(1);		// Check G exists
  double rval;
  switch(type)
    {
    default:
    case 0: rval = (*this)[spin].A();     break;	// Here for A
    case 1: rval = (*this)[spin].eta();   break;	// Here for asymmetry
    case 2: rval = (*this)[spin].theta(); break;	// Here for theta
    case 3: rval = (*this)[spin].phi();   break;	// Here for phi
//    case 4: rval = (*this)[spin].delz();  break;	// Here for delzz
    }
  return rval;
  }


// ---------------------------------------------------------------------------- 
//                                 G Values
// ---------------------------------------------------------------------------- 

	// Input		IHFV	: Hyperfine interaction vector
	// 			spin	: G index
	// 			dz	: G coupling constant (Hertz)
	// Output		none	: Get/Set G spin coupling
	// Note				: Defined in class IntHF as equal to
	//			          the G tensor delzz value

//void   IntHFVec::G(int spin, double dz)  { CValue(spin, dz, 0); }
//double IntHFVec::G(int spin) const       { return CValue(spin, 0); }
void   IntHFVec::delz(int spin, double dz) { CValue(spin, dz, 0); }
double IntHFVec::delz(int spin) const      { return CValue(spin, 0); }

// ---------------------------------------------------------------------------- 
//                        G Asymmetry Values
// ---------------------------------------------------------------------------- 

	// Input		IHFV	: Hyperfine interaction vector
	// 			spin	: G index
	// 			deta	: Hyperfine interaction asymmetry
	// Output		none	: Get/Set G spin asymmetry
	// Note				: Defined in class IntG between [0,1]
	// Note				: Very unusual if nonzero!

void   IntHFVec::eta(int spin, double ceta) { CValue(spin, ceta, 1); }
double IntHFVec::eta(int spin) const        { return CValue(spin, 1); }

 
// ---------------------------------------------------------------------------- 
//               G Theta Orientation (Down From +z Axis)
// ---------------------------------------------------------------------------- 

	// Input		IHFV	: Hyperfine interaction vector
	// 			spin	: G index
	// 			dtheta	: Hyperfine interaction angle (deg)
	// Output		none	: Get/Set G spin theta angle
	// Note				: Defined class IntG between [0,180]

void   IntHFVec::theta(int spin, double dtheta) { CValue(spin, dtheta, 2); }
double IntHFVec::theta(int spin) const          { return CValue(spin, 2); }

// ---------------------------------------------------------------------------- 
//               G Phi Orientation (Over From +x Axis)
// ---------------------------------------------------------------------------- 

	// Input		IHFV	: Hyperfine interaction vector
	// 			spin	: G index
	// 			dphi	: Hyperfine interaction angle (deg)
	// Output		none	: Get/Set G spin phi angle
	// Note				: Defined in IntG between [0,360]

void   IntHFVec::phi(int spin, double dphi)   { CValue(spin, dphi, 3); }
double IntHFVec::phi(int spin) const          { return CValue(spin, 3); }
 
// ---------------------------------------------------------------------------- 
//                        Full Hyperfine Interaction
// ---------------------------------------------------------------------------- 


IntHF& IntHFVec::operator() (int i) { return (*this)[i]; }
 
	// Input		IHFV   : Hyperfine interaction vector (this)
        //                      i     : A Hyperfine interaction index
        // Ouput                DI    : The i'th Hyperfine interaction in IHFV
        // Note                       : Returns a reference to the interaction


IntHF IntHFVec::get(int hfi) const

	// Input		IHFV	: Hyperfine interaction vector
	// 			hif	: Hyperfine interaction index
	// Output		DI	: Return rank 2 Hyperfine interaction

  {
  if(!check_spin(hfi)) IHFVfatality(1);		// Check HF interaction exists
  return (*this)[hfi];				// Return requested IntHF
  }
 
// ---------------------------------------------------------------------------- 
//                 Other Hyperfine Interaction Vector Info
// ---------------------------------------------------------------------------- 


//int IntHFVec::size() const { return nspins; }

	// Input		IHFV   : Hyperfine interaction vector
	// Output		nd    : Number of interactions in vector


bool IntHFVec::nonzero() const
   {
//   for(unsigned i=0; i<size(); i++)
//     if((*this)[i].delz()) return true;
   return false;
   }

// sosi - need better test above

// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//    Functions To Make A Parameter Set From A Hyperfine Interaction Vector
// ----------------------------------------------------------------------------

/* This class has no implicit knowledge of the electron an nuclear spins. The
   preferred means of specifying a hyperfine G interaction when there is no 
   implicit knowledge of the spin types is that which is used for filling up
   the parameter set to be written. The base parameters of individual Hyperfine
   interactions in that case are 

                { GI(#), gxx(#), gyy(#), gcc(#), GTheta(#), GPhi(#) }.

   There are functions that provide a means of writing the interaction(s) when
   the spin types are known, as this will be the case when a spin system
   controls the output. In this latter case the parameters used when writting
   then interaction are slightly different.

   For base parameters of individual interactions see the class IntHF. Any
   additional parameters for the interaction vector will be defined in the +=
   function below (there currently aren't any).

          Input                GI      : G factor interaction
                               idx     : Interaction index (default -1)
                               pfx     : Interaction 2nd indx (def -1)
          Output               pset    : Parameter set with only
                            OR TF      : Return is FALSE if proper parameters
                                         for the interaction were not found
                                         electron G interaction parameters
          Note                         : The parameter names & types
                                         output here MUST match those used
                                         in setting the vector up
                                         from parameters sets 

    Function                                 Purpose
  ------------         -------------------------------------------------------
  ParameterSet         Convert interaction into a parameter set
  operator +=          Adds interaction to existing parameter set (friend)
  PSetAdd              Adds interaction to existing parameter set, this
                       allows for an index and a prefix in the parameters   */

IntHFVec::operator ParameterSet( ) const
  { ParameterSet pset; pset += *this; return pset; }

void operator+= (ParameterSet& pset, const IntHFVec &IHFV)
  { IHFV.PSetAdd(pset); }

void IntHFVec::PSetAdd(ParameterSet& pset, int idx) const
  {
  for(int i=0; i<int(size()); i++)
    (*this)[i].PSetAdd(pset, i, idx);
  }


// ---------------------------------------------------------------------------- 
//     Functions To Make A Hyperfine Interaction Vector From A Parameter Set
// ---------------------------------------------------------------------------- 
	// Input		IHFV	: Hyperfine interaction vector (this)
	// 			pset	: A parameter set
	// Output		none	: Hyperfine interaction vector filled
	//				  with parameters in pset
	// Note                         : We just use a member function
	//                                which is needed to allow for
	//                                [#] prefixed parameter names

        // Input                IHFV     : Hyperfine interaction vector (this)
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


void IntHFVec::operator= (const ParameterSet& pset) { setIHFVec(pset); }
/*
bool IntHFVec::setIHFVec(const ParameterSet& pset, int idx, int warn)
  {
//      Use Only Parameters With Prefix [#], So Clip [#] From Names First

  ParameterSet  subpset;                        // Working parameter set
  if(idx != -1) subpset=pset.strip(idx);        // Get only params with [#]
  else          subpset=pset;                   // Or use full pset

//  First Attempt To Determine The Number Of G Interactions In Parameter Set

  bool TF=true;				// Track if we read OK
  int nd = getNInts(pset, idx);		// Get # of Hyperfine interactions
  if(nd<1)				// If we can't find any then
    {					// we'll have to abort the mission
    IHFVerror(3,1);
    IHFVfatality(2);
    }

//  Now Attempt To Read nd Hyperfine Interactions From Values In pset
//  Each IntG As { gxx,gyy,gzz,Gtheta,Gphi } or { Giso,Gdelz,Geta,Gtheta,Gphi }

  clear();				// Remove any existing interactions
//  for(int i=0; i<nd; i++)		// Loop interactions and set 'em
//    push_back(IntHF(subpset,i,1));	// using the base class
  IntHF G;
  for(int i=0; i<nd; i++)		// Loop interactions and set 'em
    {
    if(G.read(subpset,i));
      push_back(G);
    }
  return TF;
  }
*/


// ----------------------------------------------------------------------------
//     Output Hyperfine Interaction Vector To ASCII From A Parmeter Set
// ----------------------------------------------------------------------------

int IntHFVec::write(const string &filename, int idx, int warn) const

	// Input		IHFV   	: Hyperfine interaction vector (this)
	//			filename: Output file name
        //                      idx	: Parameter index value used for
	//				  prefix [#] in output names
        //                      warn    : Warning level
	// Output		none 	: Hyperfine interaction vector is written
	//				  as a parameter set to file filename

  {
  if(!size()) return 1;
  ofstream ofstr(filename.c_str());	// Open filename for input
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!write(ofstr, idx, w2))		// If file bad then exit
    {
    IHFVerror(1, filename, 1);		// Filename problems
    if(warn>1) IHFVfatality(20);		// Fatal error
    return 0;
    }
  ofstr.close();                        // Close it now
  return 1;
  }


int IntHFVec::write(ofstream& ofstr, int idx, int warn) const

	// Input		IHFV   	: Hyperfine interaction vector (this)
        //                      ofstr   : Output file stream
        //                      idx     : Parameter index value used for
	//				  prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Hyperfine interaction vector is
        //                                written as a parameter set to
        //                                output filestream
 
  {
  if(!size()) return 1;
  if(!ofstr.good())			// If file bad then exit
    {   
    if(warn) IHFVerror(22);		//      Problems with file
    if(warn > 1) IHFVfatality(23);	//      It's a fatal error
    return 0;
    }   
  ParameterSet pset;			// Declare a parameter set
  for(int i=0; i<int(size()); i++)	// Add each Hyperfine interaction
    ((*this)[i]).PSetAdd(pset, i, idx);	// to the parameter set
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!pset.write(ofstr, w2))            // Use parameter set to write
    {
    if(warn)
      {
      IHFVerror(22, 1);			// Problems writing to filestream
      if(warn>1) IHFVfatality(23);	// Fatal error
      }
    return 0;
    }  
  return 1;
  }  


// ____________________________________________________________________________ 
// D               HYPERFINE INTERACTION VECTOR INPUT FUNCTIONS
// ____________________________________________________________________________ 

// ----------------------------------------------------------------------------
//        Direct Read of Vector From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

/* These functions allow users to fill the vector of interactions either from
   parameters found in an external ASCII file or contained in a GAMMA parameter
   set.
	   Input		IHFV	: Hyperfine interaction vector (this)
	   			filename: Input filename
                             or pset    : Parameter set
                                idx     : Parameter index value used for
	  				  prefix [#] in parameter names
                                warn    : Warning output level
                                               0 = no warnings
                                               1 = warnings
                                              >1 = fatal warnings
	   Output		none	: Hyperfine interaction vector filled
	  				  with parameters read from file     */

bool IntHFVec::read(const string &filename, int idx, int warn)
  {
  ParameterSet pset;                    // Declare a parameter set
  if(!pset.read(filename, warn?1:0))    // Try and read in pset
    {                                   // If we cannot read the file then
    if(warn)                            // we'll issue warnings as desired
      {
      IHFVerror(1, filename, 1); 	//      Problems with file
      if(warn>1) IHFVfatality(21);	//      This is a fatal problem
      else       IHFVerror(21, 1);	//      Or maybe it aint so bad
      }
    return false;
    }
  return read(pset, idx, warn);
  }

bool IntHFVec::read(const ParameterSet& pset, int idx, int warn) 
  {
  bool TF = setIHFVec(pset, idx, warn?1:0);	// Use overload to read
  if(!TF)                                       // If setIGvec didn't handle
    {                                           // the vector read from pset
    if(warn)                                    // then we'll issue errors
      {
      IHFVerror(8, 1);				//    Problems with pset
      if(warn>1) IHFVfatality(21);		//    This is a fatal problem
      if(warn>1) IHFVerror(21);			//    Or maybe it isn't so bad..
      }
    return false;
    }
  return TF;
  }


// ----------------------------------------------------------------------------
//       Interactive Read of Hyperfine Interaction Vector From An ASCII File
// ----------------------------------------------------------------------------


string IntHFVec::ask_read(int argc, char* argv[], int argn)

	// Input		IHFV    : Hyperfine interaction vector (this)
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
   "\n\tHyperfine Interaction Vector Filename? ", 
                                      filename);
  read(filename);		           	// Read system from filename
  return filename;
  }

// ____________________________________________________________________________ 
// E                 HYPERFINE INTERACTION VECTOR OUTPUT FUNCTIONS
// ____________________________________________________________________________ 

/* These functions provied ASCII output of the interaction vector into any
   output stream, including C++ standard output.

	   Input		IHFV	: Hyperfine interaction vector (this)
	   			ostr	: Output stream
	  			full    : Flag for long vs short output
	   Output		non	: Hyperfine interaction vector
					  parameters sent to output stream  */


ostream& IntHFVec::print(ostream& ostr, int full) const
  {
  string hdr="Hyperfine Interactions Vector";	// Title for output
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
  for(int j=0; j<int(size()); j++)		// Output Hyperfine interactions
    {						// in order
    hdr = string("Interaction ") + Gdec(j);	//	Make a new header
    ctr=string(40-hdr.length()/2, ' ');		// 	Spacer to center output 
    ostr << "\n" << ctr << hdr;			//	Output a header
    ((*this)[j]).print(ostr, full);		//	Output the interaction
    ostr << "\n";				// 	Add a line spacer
    }
  return ostr;
  }


ostream& operator<< (ostream& out, const IntHFVec& IHFV)
  { return IHFV.print(out); }

#endif						// IntHFVec.cc
