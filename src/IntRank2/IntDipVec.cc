/* IntDipVec.cc *************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Dipolar Interactions Vector 		Implementation		**
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
** This class maintains a vector of rank2 dipolar interactions (from	**
** class IntRank2).  The class allows users to manipulate an array of	**
** such interactions simultaneously.					**
**                                                                      **
*************************************************************************/

#ifndef   IntDipVec_cc_			// Is file already included?
#  define IntDipVec_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <IntRank2/IntDipVec.h>		// Include header file
#include <IntRank2/IntDip.h>		// Include dipolar interactions
#include <Basics/Gutils.h>		// Includes query knowledge
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Level1/coord.h>		// Include single coordinates
#include <Level1/coord_vec.h>		// Include coordinate vectors
#include <Basics/StringCut.h>		//   Done in StringCut
#include <stdlib.h>

using std::cout;			// Using libstdc++ standard output
using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ STL lists
using std::vector;			// Using libstdc++ STL vectors
using std::ofstream;			// Using libstdc++ output file streams
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------  
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                 CLASS DIPOLE INTERACTION VECTOR ERROR HANDLING
// ____________________________________________________________________________
 


	// Input		IDV	: Dipole interaction vector (this)
        // 			eidx	: Error index
        //                      noret	: Flag for return (0=linefeed)
        // Output               none	: Error message

void IntDipVec::IDVerror(int eidx, int noret) const
  {
  cout << "\nClass IntDipVec: ";
  switch(eidx)
    {
    case 0:								// (0)
      cout << "Program Aborting.....";
      break;
    case 1:								// (1)
      cout << "Cannot Access Dipolar Interaction Between Spins";
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
      cout << "Attempted Dipole Access of Spin with Itself";
      break;
    case 14:								// (14)
      cout << "Sorry, Dipolar Tensor Operation Not Allowed Yet";
      break;
    case 15:								// (15)
      cout << "Electron - Nucleus Spin Pair Disallowed";
      break;
    case 19:								// (19)
      cout << "Can't Fill Dipolar Interaction Vector from Parameter Set";
      break;
    case 20:								// (20)
      cout << "Can't Write Dipolar Interaction Vector to Parameter File";
      break;
    case 21:								// (21)
      cout << "Can't Read Dipolar Interaction Vector from Parameter File";
      break;
    case 22:                                                            // (22)
      cout << "Problems Writing Dipolar Interaction Vector to Output FileStream";
      break;
    case 23:                                                            // (23)
      cout << "Cannot Output Dipolar Interaction Vector Parameters";
      break;
    case 24:                                                            // (24)
      cout << "Cannot Access Requested Dipolar Interaction";
      break;
    case 25:                                                            // (25)
      cout << "Problems Setting Individual Dipolar Interactions";
      break;
    case 30:                                                            // (30)
      cout << "Cannot Determine The Number of Spins in Parameter File";
      break;
    case 31:                                                            // (31)
      cout << "Setting Spin To Default Isotope Type";
      break;
    case 32:                                                            // (32)
      cout << "Setting Interaction Between Spins To Zero";
      break;
    default: 
      cout << "Unknown error";
      break;
    }
  if(!noret) cout << ".\n";
  }


     
	// Input		IDV	: Dipole interaction vector (this)
        // 			eidx	: Error index
     	// Output		none	: Error message output
	//				  Program execution stopped
     
volatile void IntDipVec::IDVfatal(int eidx) const
  {
  IDVerror(eidx, 1);			// First output the error
  if(eidx) IDVerror(0);			// Now output it fatal
  cout << "\n";				// Keep the screen nice
  exit(-1);				// Abort the program
  }



	// Input		IDV	: Dipole interaction vector (this)
        // 			eidx	: Error index
        //                      pname	: Additional string for error
        //                      noret	: Flag for return (0=linefeed)
        // Output               none	: Error message

void IntDipVec::IDVerror(int eidx, const string& pname, int noret) const
  {
  cout << "\nClass IntDipVec: ";
  switch(eidx)
    {
    case 100:							// (100)
      cout << "Can't Read Parameter " << pname;
      break;
    case 101:							// (101)
      cout << "Problems with File " << pname;
      break;
    case 120:							// (120)
      cout << "Dipole Access Of Index " << pname << " Out Of Bounds";
      break;
    case 121:							// (121)
      cout << "Cannot Access Dipolar Interaction " << pname;
      break;
    case 122:							// (122)
      cout << "Cannot Access Dipolar Interaction Between Spins " << pname;
      break;
    case 123:							// (123)
      cout << "Cant Set Interaction Between Spin Pairs: " << pname;
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
     
volatile void IntDipVec::IDVfatal(int eidx, const string& pn) const
  {
  IDVerror(eidx, pn, 1);		// First output the error
  if(eidx) IDVerror(0);			// Now output it fatal
  cout << "\n";				// Keep the screen nice
  exit(-1);				// Abort the program
  }

// ____________________________________________________________________________
// ii               DIPOLE INTERACTION VECTOR SETUP FUNCTIONS
// ____________________________________________________________________________
 
/* These functions determine specific aspects of a dipolar interactions vector
   from a specified parameter set. As the various interaction parameters tend
   to interact, the functions MUST be private because their misuse could
   produce an inconsistent interactions vector.

   The goal here is quite simple. We need to determine a series of dipolar
   interactions from parameters found in an input GAMMA parameter set. Each
   interaction will be defined by unique parameters in the set that specify
   the values: { Iqn, Sqn, DCC ,eta, alpha, beta, gamma } Complexity arises 
   because we need find a series of such interactions and because dipolar
   interactions allow for a variety of parameters for their definition.

   Consider the kth dipolar interaction in the vector. It may have been 
   defined by parameters indexed with a single interaction number or it may
   have been defined by a pair of spin indices. If the former is true it is
   clear that the kth interaction is defined by the parameters using index
   k. If the latter is used what spin indices correspond to the kth dipolar
   interaction? That is a problem unless we know the total number of spins
   that are specified, and that is what we will demand. In order to use
   spin indexed parameters we MUST know the total number of spins defined.   */
   
// ---------------------------------------------------------------------------- 
//    Functions To Set Number Of Spins (For Using Spin Indexed Parameters)
// ---------------------------------------------------------------------------- 

/* As mentioned in the previous paragraphs, without knowledge of the total
   number of spins we cannot translate two spin indices into a the single
   interaction index the sets the position of the interaction in the vector.
   We allow several means of setting the total number of spins. These are as
   follows:

	1.) NSpins	The parameter NSpins will be assumed that number
        2.) Iso         A series of parameters Iso(i),   i=0,1,2,.. sets it
        3.) Dqn         A series of parameters Dqn(i),   i=0,1,2,.. sets it
	4.) DCC         A series of parameters DCC(i,j), i=0,1,2,.. sets it
	5.) Dxx         A series of parameters Dxx(i,j), i=0,1,2,.. sets it

   Note this does NOT allow for mixing up dipolar interaction designations.
   That is, one may not use Iso(i) & Iso(j) for one spin pair then switch to
   Dqn(i) & Dqn(j) for another spin pair. Stay consistent when using spin
   indices for dipolar designations.                                         */

bool IntDipVec::getNS(const ParameterSet& pset, int idx, 
                                                     int& ns, bool warn) const
  {
  if(getNSpins(pset,idx,ns,false))		// If parameter NSpins known
    return true;				// we are done
  ns = 0;					// We know of no spins yet
  int nns = 0;					// Working spin count
  getNIsos(pset,idx,nns,false);			// Try for ns using Iso params
  ns = gmax(nns, ns);				// Number of spins is max this
  getNqns(pset,idx,nns,false);			// Try for ns using Dqn params
  ns = gmax(nns, ns);				// Number of spins is max this
  getNdccs(pset,idx,nns,false);			// Try for ns using DCC params
  ns = gmax(nns, ns);				// Number of spins is max this
  getNdxxs(pset,idx,nns,false);			// Try for ns using Dxx params
  ns = gmax(nns, ns);				// Number of spins is max this
  if(ns) return true;
  return false;
  }

bool IntDipVec::getNSpins(const ParameterSet& pset, int indx, 
                                                     int& ns, bool warn) const
  {
  string prefix;				// Parameter prefix (if any)
  if(indx != -1)				// Set parameter name prefix
    prefix = string("[")			// [indx] if indx != -1
           + Gdec(indx) + string("]");
  string pname = prefix + string("NSpins");	// Set parameter name
  ParameterSet::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Try and get parameter
  if(item != pset.end()) 			// Retrieve the number of spins
    {
    string pstate;
    (*item).parse(pname,ns,pstate);		//   We found it, set ns
    return true;				//   Return we succeded
    }
  if(warn) IDVerror(30, 1);			// Can't determine # spins
  ns = 0;					// Insure no spins known
  return false;					// Signal our failure
  }
 
bool IntDipVec::getNIsos(const ParameterSet& pset,int pfx,
                                                     int& ni, bool warn) const
  {
  string prefix;				// Parameter prefix (if any)
  if(pfx != -1)					// Set parameter name prefix
    prefix = string("[")			// [pfx] if pfx != -1
           + Gdec(pfx) + string("]");
  string pbase1("Iso(");			// Parameter name base
  string pbase2("Ion(");			// Parameter name base
  string pend(")");				// Parameter name end
  string pname;					// For full parameter name
  bool OK = true;				// Flag if we know parameter
  ni = 0;					// Number of isotopes found
  ParameterSet::const_iterator item;         // A pix into parameter list
  while(OK)
    {
    pname = prefix + pbase1 + Gdec(ni) + pend; 	//   Set parameter name
    item = pset.seek(pname);			//   Search for parameter
    if(item != pset.end()) ni++;		//   Continue if found
    else
      {
      pname = prefix + pbase2 + Gdec(ni) + pend;//   Set parameter name
      if(item != pset.end()) ni++;		//   Continue if found
      else                   OK = false;	//   Quit if not
      }
    }
  if(ni) return true;				// If any found successful
  if(warn) IDVerror(30, 1);			// Can't determine # spins
  ni = 0;					// Insure no spins known
  return false;					// Else signal we failed
  }
 
bool IntDipVec::getNqns(const ParameterSet& pset,int indx,
                                                     int& nq, bool warn) const
  {
  string prefix;				// Parameter prefix (if any)
  if(indx != -1)				// Set parameter name prefix
    prefix = string("[")			// [indx] if indx != -1
           + Gdec(indx) + string("]");
  string pbase("Dqn(");				// Parameter name base
  string pend(")");				// Parameter name end
  string pname;					// For full parameter name
  bool OK = true;				// Flag if we know parameter
  nq = 0;					// No. of quantum #'s found
  ParameterSet::const_iterator item;         // A pix into parameter list
  while(OK)
    {
    pname = prefix + pbase + Gdec(nq) + pend; 	//   Set parameter name
    item = pset.seek(pname);			//   Search for parameter
    if(item != pset.end()) nq++;		//   Continue if found
    else                   OK = false;		//   Quit if not
    }
  if(nq) return true;				// If any found successful
  if(warn) IDVerror(30, 1);			// Can't determine # spins
  nq = 0;					// Insure no spins known
  return false;					// Else signal we failed
  }

bool IntDipVec::getNdccs(const ParameterSet& pset,int indx,
                                                     int& nd, bool warn) const
  {
  string prefix;				// Parameter prefix (if any)
  if(indx != -1)				// Set parameter name prefix
    prefix = string("[")			// [indx] if indx != -1
           + Gdec(indx) + string("]");
  string pbase("DCC(0,");			// Parameter name base
  string pend(")");				// Parameter name end
  string pname;					// For full parameter name
  bool OK = true;				// Flag if we know parameter
  nd = 1;					// No. of quantum #'s found
  ParameterSet::const_iterator item;         // A pix into parameter list
  while(OK)
    {
    pname = prefix + pbase + Gdec(nd) + pend; 	//   Set parameter name
    item = pset.seek(pname);			//   Search for parameter
    if(item != pset.end()) nd++;		//   Continue if found
    else                   OK = false;		//   Quit if not
    }
  if(nd>1) return true;				// If any found successful
  if(warn) IDVerror(30, 1);			// Can't determine # spins
  nd = 0;					// Insure no spins known
  return false;					// Else signal we failed
  }

bool IntDipVec::getNdxxs(const ParameterSet& pset,int indx,
                                                     int& nd, bool warn) const
  {
  string prefix;				// Parameter prefix (if any)
  if(indx != -1)				// Set parameter name prefix
    prefix = string("[")			// [indx] if indx != -1
           + Gdec(indx) + string("]");
  string pbase("Dxx(0,");			// Parameter name base
  string pend(")");				// Parameter name end
  string pname;					// For full parameter name
  bool OK = true;				// Flag if we know parameter
  nd = 1;					// No. of quantum #'s found
  ParameterSet::const_iterator item;		// A pix into parameter list
  while(OK)
    {
    pname = prefix + pbase + Gdec(nd) + pend; 	//   Set parameter name
    item = pset.seek(pname);			//   Search for parameter
    if(item != pset.end()) nd++;		//   Continue if found
    else                   OK = false;		//   Quit if not
    }
  if(nd>1) return true;				// If any found successful
  if(warn) IDVerror(30, 1);			// Can't determine # spins
  nd = 0;					// Insure no spins known
  return false;					// Else signal we failed
  }

// ----------------------------------------------------------------------------
//         Functions To Set The Entire Dipolar Interactions Vector
// ----------------------------------------------------------------------------

/* We allow class IntDip to read inidividual interactions from the parameter
   set. This simplifies our task down to dealing with the interaction indexing.
   If a single index is used per interaction this is simple because then the
   interaction index is one in the same as the vector index. If two indices
   (i.e. spin indices) are used per interactions we have a bit more trouble.
   In that case we MUST know the total number of spins we are dealing with in
   order to convert indices i & j into the vector index k. So, first we look
   for the total number of spins set in the parameters.  If found we will allow
   spin indexing and interaction indexing, if not we only allow interaction
   indexing.

	   Input		IDV	: Dipole interaction vector (this)
	   			pset	: A parameter set
	  			idx     : Parameter prefix index
				warn    : Flag for warnings
	   Output		none	: System dipolar interactions are
	  			  	  set from parameters in pset

   We Use Only Parameters With Prefix [#] (unles idx=-1). As such we off this
   prefix from all parameter names prior to looking for the interacitons. 
   This is done before any calls to class IntDip which doesn't use prefixes. */
 
bool IntDipVec::setDIV(const ParameterSet& pset, int idx, bool warn)
  {
  clear();                                      // Remove existing intacts.
  ParameterSet  subpset;			// Copy parameter set
  if(idx != -1) subpset = pset.strip(idx);	// to glean out those for
  else          subpset = pset;			// the specified prefix index

//     Try To Fill Vector Using Spin Pairs When # Of Spins Is Specified
//             (Or We Can Figure Out The # Of Spins Specified)

  int ns;					// Number of spins to track
  IntDip D, D0;					// Empty dipolar interaction
  string udef;					// Undefined interactions
  if(getNS(pset,idx,ns,0))			// If we know the # of spins
    {						// we allow spin indexing
    int i,j;					//   Two spin indices
    for(i=0; i<ns-1; i++)			//   Loop over i spins
      {
      for(j=i+1; j<ns; j++)			//   Loop over j spins
        {					//   Try and get D(i,j)
        D = D0;					//     Begin with empty D
        if(!D.read(subpset,i,j, false))		//     If unsuccessful read,
          {					//     store spin indices
          if(udef.length()) udef += string(", ");	//     of failed interact.
          udef += string("(") + Gdec(i)
               +  string(",") + Gdec(j)
               + string(")");
          }					//    Either way, store int
        push_back(D);				//    into the vector (may
        }					//    be empty)
      }
    if(udef.length() && warn)
      {
      if(warn)
        {
        IDVerror(25, 1);
        IDVerror(123, udef, 1);
        }
      return false;
      }
    return true;
    }

//      Try To Fill Vector With Either Spin Pair Or Interaction Indices

    int k = 0;					//   Start with 0th D
    while(D.read(subpset,k))			//   If we are able to read it
      {						//   then we keep looking
      push_back(D);				//     Put into the vector
      k++;					//     Index to next interaction
      }
   
  if(!size()) return false;
  return true;
  }

// ____________________________________________________________________________
// iii             DIPOLE INTERACTIONS VECTOR CHECKING FUNCTIONS
// ____________________________________________________________________________

/* These functions insure vector integrity and that access is not beyond
   the vector boundaries

           Input                IDV     : Dipole interaction vector (this)
                                dip     : A dipole index
                                warn    : Warning level
           Output               TF      : Returns TRUE if dip is a valid
                                          dipole index, FALSE if not         */

bool IntDipVec::CheckDI(int dip, int warn) const
  {
  if(dip>=int(size()) || dip<0)
    {
    if(warn)
      {
      IDVerror(120, Gdec(dip), 1);
      IDVerror(121, Gdec(dip), 1);
      }
    if(warn>1) IDVfatal(24);
    return 0;
    }
  return 1;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// A              DIPOLE INTERACTION VECTOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

/* These functions should be inherited from the ANSI standard "vector" class
   in the C++ library libstdc++? I am listing them here so that I & other
   users don't have to keep looking them up all the time and because the
   inheritance doesn't seem to always work.                                  */

IntDipVec::IntDipVec(int N) 
  { 
  if(N>0) 
    {
    IntDip DI;
    for(int i=0; i<N; i++) push_back(DI); 
    }
  }
IntDipVec::IntDipVec(const IntDipVec &IDV1) { *this = IDV1; }

// ----------------------------------------------------------------------------
//             Constructors Using Isotopes And Spin Coordinates
// ----------------------------------------------------------------------------

/* A pair of spin isotopes and associated spin coordinates defines a set of
   dipolar interactions. Problems will occur if any of the spin pairs mix up
   electrons and nuclei, since such a pairing is disallowed. A flag is added to
   let the user ignore such pairings (leaving an empty dipolar interaction),
   ignore these pairings but just leave the iteractions empty, or fail entirely
   if such a pairing is encountered. Note that this type of constructor is very
   much tailored for use in spin systems containing mutiple dipole-dipole
   interactions                                                              */

	// Input		IDV	: Dipole interaction vector (this)
	//			Isos	: Array of isotope types
	//			cvec	: Array of spin coordinates
        //                      warn    : Warning level
	//				    0 - no warnings
	//				    1 - warn if e-/nucleus pair
	//				    2 - fatal if e-/nucleus pair
	// Output		none	: Dipole interaction vector constructed
	//				  with isotope labels in Isos and spin
	//				  coordinates in cvec.
	// Note				: Vector size set by cvec size,
	//				  Array Isos must contain >= this size
	// Note				: We disallow dipolar interactions that
	//				  involve electron-nucleus pairs!  They
	//				  will trigger a NULL interaction in 
	//				  the vector & issue cause warnings if
	//				  warn is set nonzero

IntDipVec::IntDipVec(const vector<Isotope>& Isos,
                                               const coord_vec& cvec, int warn)
  {
  int nc = cvec.size();				// Number of spin coordinates
  int ni = Isos.size();				// Number of isotope types
  int ns = gmin(nc,ni);				// Number spins is smallest 
  Isotope ISI, ISS;				// For spin isotopes
  coord   ptI, ptS;				// For spin coordinates
  int i, j, k;					// Dummy variables
  IntDip DEmpty;				// An empty dipolar interaction
  for(i=0, k=0; i<ns-1; i++)			// Set dipolar interactions
    {						// First loop all I spins
    ptI = cvec.get(i);				//   Spin I coordinate
    ISI = Isos[i];				//   Spin I isotope
    for(j=i+1; j<ns; j++, k++)			//   Now loop all J>I spins
      {
      ptS = cvec.get(j);			//     Spin S coordinate
      ISS = Isos[j];				//     Spin S isotope
      if(!ISI.nepair(ISS))			//     Insure not an electron
        push_back(IntDip(ISI,ISS,ptI,ptS));	//     & set {I,S} interaction
      else
        {					//     If so, we will NOT set
        if(warn)				//     the interaction. If warn
          {					//     we will issue warnings &
          IDVerror(15, 1);			//     cease to run if desired.
          if(warn>1) IDVfatal(2);
          }
        push_back(DEmpty);
        }
      }
    }
  }

// ----------------------------------------------------------------------------
//                    Constructors Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions allow users to fill up the dipolar interactions vector from
   a GAMMA parameter set (external ASCII file). Variants enable higher classes
   such as spin systems to generate the vector when they have a knowledge of
   spin isotope types (and coordinates in the dipolar interaction case).

	// Input		IDV	: Dipole interaction vector (this)
        //                      pset    : Parameter set
        //                      idx	: Index for dipolar vector
	//				  i.e [idx] is parameter name prefix
        //                      warn    : Warning level (applicable to use
	//				  of isotopes & coordinates)
        // Note                         : If idx is negative then no
        //                                parameter prefix will be used     */
 
IntDipVec::IntDipVec(const ParameterSet& pset, int idx, int warn)
  {
  if(!setDIV(pset, idx, warn?true:false))
    if(warn)
      {
      IDVerror(3, 1);
      if(warn>1) IDVfatal(2);
      else       IDVerror(2);
      }
  }

IntDipVec::IntDipVec(const ParameterSet& pset, 
                              const vector<Isotope>& Isos, int idx, int warn)
  {
  int ns = Isos.size();				// Number of spins
  if(!ns) return;				// Bail if no dipoles
  ParameterSet subpset;				// Copy parameter set
  if(idx != -1) subpset = pset.strip(idx);	// to glean out those for
  else          subpset = pset;			// the specified prefix index
  Isotope I,S;					// Spin isotopes
  int i, j;					// Individual spin indices
  IntDip D, D0;					// Working (NULL) interactions
  for(i=0; i<ns-1; i++)				// Set dipolar interactions
    { 						// per each unique spin pair
    I = Isos[i];				//   Isotope type of 1st spin	
    for(j=i+1; j<ns; j++)			//   Loop over 2nd spins
      {
      S = Isos[j];				//   Isotope type of 2nd spin	
      D = D0;					//   Insure empty interaciotn
      if(!I.nepair(S))				//   If not e-/nucleus pair
        D.read(subpset,i,j,0);			//   we read in the interaction
      push_back(D);				//   Put interaction in vector
      }
    }
  if(!size() && warn)
    {
    IDVerror(3, 1);
    if(warn>1) IDVfatal(2);
    else       IDVerror(2);
    }
  }

// ----------------------------------------------------------------------------
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------

/* Both assignment and destruction are handled by the base vector class.     */

void IntDipVec::operator= (const IntDipVec& IDV) { vector<IntDip>::operator=(IDV); }
     IntDipVec::~IntDipVec () {} 	

// ____________________________________________________________________________ 
// B                DIPOLAR INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//                   Generic Dipolar Value Access Functions
// ---------------------------------------------------------------------------- 
        
        // Input                IDV     : Dipole interaction vector
        //                      dip     : Dipole index
        //                      val     : Dipolar interaction value
        //                                <=0: DCC coupling of spin (Hertz)
        //                                  1: Asymmetry [0, 1] 
        //                                  2: Theta (degrees, down from +z)
        //                                  3: Phi   (degrees, over from +z)
        // Output               none    : Get/Set dipolar interaction value
        //                                for spicified dipole dip
 
void IntDipVec::DValue(int dip, double val, int type)
  {
  if(!CheckDI(dip)) IDVfatal(1);		// Check dipole exists
  switch(type)
    {
    default:
    case 0: ((*this)[dip]).DCC(val);   break;	// Here for DCC
    case 1: ((*this)[dip]).eta(val);   break;	// Here for asymmetry
    case 2: ((*this)[dip]).theta(val); break;	// Here for theta
    case 3: ((*this)[dip]).phi(val);   break;	// Here for phi
    }
  }

double IntDipVec::DValue(int dip, int type) const
  {
  if(!CheckDI(dip)) 				// Insure dipole exists
    { IDVerror(1, 1); IDVfatal(120, Gdec(dip)); }
  double rval;
  switch(type)
    {
    default:
    case 0: rval = ((*this)[dip]).DCC();   break;	// Here for DCC
    case 1: rval = ((*this)[dip]).eta();   break;	// Here for asymmetry
    case 2: rval = ((*this)[dip]).theta(); break;	// Here for theta
    case 3: rval = ((*this)[dip]).phi();   break;	// Here for phi
    }
  return rval;
  }


// ---------------------------------------------------------------------------- 
//                         Dipolar Coupling Constants
// ---------------------------------------------------------------------------- 

	// Input		IDV	: Dipole interaction vector
	// 			dip	: Dipole index
	// 			dcc	: Dipolar coupling constant (Hertz)
	// Output		none	: Get/Set dipole dip coupling
	// Note				: Defined in class IntDip as equal to
	//			          the dipolar tensor delzz value

void   IntDipVec::DCC(int   dip, double dcc) { DValue(dip, dcc, 0); }
void   IntDipVec::Ddelz(int dip, double dcc) { DValue(dip, dcc, 0); }

double IntDipVec::DCC(int   dip) const       { return DValue(dip, 0); }
double IntDipVec::Ddelz(int dip) const       { return DValue(dip, 0); }

// ---------------------------------------------------------------------------- 
//                        Dipolar Asymmetry Values
// ---------------------------------------------------------------------------- 

	// Input		IDV	: Dipole interaction vector
	// 			dip	: Dipole index
	// 			deta	: Dipolar interaction asymmetry
	// Output		none	: Get/Set dipole dip asymmetry
	// Note				: Defined in class IntDip between [0,1]
	// Note				: Very unusual if nonzero!

void   IntDipVec::Deta(int dip, double deta) { DValue(dip, deta, 1); }
double IntDipVec::Deta(int dip) const        { return DValue(dip, 1); }

 
// ---------------------------------------------------------------------------- 
//               Dipolar Theta Orientation (Down From +z Axis)
// ---------------------------------------------------------------------------- 

	// Input		IDV	: Dipole interaction vector
	// 			dip	: Dipole index
	// 			dtheta	: Dipolar interaction angle (deg)
	// Output		none	: Get/Set dipole dip theta angle
	// Note				: Defined class IntDip between [0,180]

void   IntDipVec::Dtheta(int dip, double dtheta) { DValue(dip, dtheta, 2); }
double IntDipVec::Dtheta(int dip) const          { return DValue(dip, 2); }

// ---------------------------------------------------------------------------- 
//               Dipolar Phi Orientation (Over From +x Axis)
// ---------------------------------------------------------------------------- 

	// Input		IDV	: Dipole interaction vector
	// 			dip	: Dipole index
	// 			dphi	: Dipolar interaction angle (deg)
	// Output		none	: Get/Set dipole dip phi angle
	// Note				: Defined in IntDip between [0,360]

void   IntDipVec::Dphi(int dip, double dphi)   { DValue(dip, dphi, 3); }
double IntDipVec::Dphi(int dip) const          { return DValue(dip, 3); }
 
// ---------------------------------------------------------------------------- 
//                        Full Dipolar Interaction
// ---------------------------------------------------------------------------- 

/*	   Input		IDV	: Dipole interaction vector
	   			dip	: Dipole index
	   Output		DI	: Return rank 2 dipolar interaction
	   Note			        : These functions/operators return
					  either a reference, constant ref,
					  or copy of the dipolar interaction */

      IntDip& IntDipVec::operator() (int dip)       { return (*this)[dip]; }
const IntDip& IntDipVec::getcref(int     dip) const { return (*this)[dip]; }
      IntDip  IntDipVec::get(int         dip) const
  {
  if(!CheckDI(dip)) IDVfatal(1);		// Check dipole exists
  return (*this)[dip];				// Return requested IntDip
  }
 
// ---------------------------------------------------------------------------- 
//                 Other Dipolar Interaction Vector Info
// ---------------------------------------------------------------------------- 


double IntDipVec::Izval(int dip) const { return (*this)[dip].Izval(); }
double IntDipVec::Szval(int dip) const { return (*this)[dip].Izval(); }
 
        // Input                IDV   : Dipole interaction vector
        // Output               qn    : Quantum number of I or S (0.5, 1.5,..)


// int IntDipVec::size() const 					INHERITED
bool IntDipVec::nonzero() const
   {
   int nd = size();
   for(int i=0; i<nd; i++)
     if((*this)[i].DCC()) return true;
   return false;
   }

// ____________________________________________________________________________
// C                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//     Functions To Make A Parameter Set From A Dipolar Interaction Vector
// ----------------------------------------------------------------------------

// Note that the preferred means of specifying dipolar interactions when there
// is no knowledge of spin isotope types nor spin coordinates are used for
// filling the parameter set.  For base parameters of individual interactions
// see class IntDip. Any additional parameters for the interaction vector will
// be defined in the += function below (there currently aren't any).  


	// Input		IDV   : Dipole interaction vector
        //                      pset  : Parameter set
        // Output               pset  : Parameter set with only
        //                              dipolar interaction parameters
	// Output		pset	: Parameter set with only
	//			          dipolar interaction parameters added
 
IntDipVec::operator ParameterSet( ) const
  { ParameterSet pset; pset += *this; return pset; }		// Add in vector parameters

void operator+= (ParameterSet& pset, const IntDipVec &IDV)
  { IDV.PSetAdd(pset); }


         
 
	// Input		IDV	: Dipole interaction vector
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

void IntDipVec::PSetAdd(ParameterSet& pset, int idx) const
  {
  int nd = size();
  for(int i=0; i<nd; i++)	 	// Add each dipole interaction
    ((*this)[i]).PSetAdd(pset, i, idx);	// to the parameter set
  }

// ----------------------------------------------------------------------------
// Functions To Output Dipolar Interaction Vector To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

int IntDipVec::write(const string &filename, int idx, int warn) const

	// Input		IDV   	: Dipole interaction vector (this)
	//			filename: Output file name
        //                      idx	: Parameter index value used for
	//				  prefix [#] in output names
        //                      warn    : Warning level
	// Output		none 	: Dipole interaction vector is written
	//				  as a parameter set to file filename

  {
  if(size()) return 1;
  ofstream ofstr(filename.c_str());	// Open filename for input
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!write(ofstr, idx, w2))		// If file bad then exit
    {
    IDVerror(101, filename, 1);		// Filename problems
    if(warn>1) IDVfatal(20);		// Fatal error
    return 0;
    }
  ofstr.close();                        // Close it now
  return 1;
  }


int IntDipVec::write(ofstream& ofstr, int idx, int warn) const

	// Input		IDV   	: Dipole interaction vector (this)
        //                      ofstr   : Output file stream
        //                      idx     : Parameter index value used for
	//				  prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Dipolar interaction vector is
        //                                written as a parameter set to
        //                                output filestream
 
  {
  if(!size()) return 1;
  if(!ofstr.good())			// If file bad then exit
    {   
    if(warn) IDVerror(22);		//      Problems with file
    if(warn > 1) IDVfatal(23);	//      It's a fatal error
    return 0;
    }   
  ParameterSet pset;			// Declare a parameter set
  int nd = size();
  for(int i=0; i<nd; i++)		// Add each dipole interaction
    ((*this)[i]).PSetAdd(pset, i, idx);	// to the parameter set
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!pset.write(ofstr, w2))            // Use parameter set to write
    {
    if(warn)
      {
      IDVerror(22, 1);			// Problems writing to filestream
      if(warn>1) IDVfatal(23);	// Fatal error
      }
    return 0;
    }  
  return 1;
  }  


// ____________________________________________________________________________ 
// D               DIPOLAR INTERACTION VECTOR INPUT FUNCTIONS
// ____________________________________________________________________________ 

// ----------------------------------------------------------------------------
//        Direct Read of Vector From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

/* These next two read functions attempt to read in a series of dipolar 
   interactions whose parameter names are indexed with the prefix [#]. The
   interactions MUST be sequentially labeled, the first one not defined will
   be signal an end of the search. The individual interactions may be indexed
   using two spin indices, suffix (#1, #,2), or interaction indices, suffix 
   (#).

	   Input		IDV	: Dipole interaction vector (this)
	   			filename: Input filename
                                pset    : Parameter set
                                idx     : Parameter index value used for
	  				  prefix [#] in output names
                                warn    : Warning output label
                                           0 = no warnings
                                           1 = warnings
                                          >1 = fatal warnings
	   Output		none	: Dipole interaction vector filled
	  				  with parameters read from file
	  				  or with parameters read from pset
                                          Return true if interaction read   */

bool IntDipVec::read(const string& filename, int idx, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename, warn?1:0))	// Read in pset from file
     {   
     if(warn)
       {
       IDVerror(101, filename, 1);	//      Problems with file 
       if(warn > 1) IDVfatal(21);	//      Its a fatal error 
       else         IDVerror(21);	//	or a warning issued
       }
     return false;
     }
  if(!read(pset, idx, warn))
     {
     if(warn)
       {
       IDVerror(101, filename, 1);	//      Problems with file 
       if(warn > 1) IDVfatal(21);	//      Its a fatal error 
       else         IDVerror(21);	//	or a warning issued
       }
     return false;
     }
  return true;
  }

bool IntDipVec::read(const ParameterSet& pset, int idx, int warn)
  {
  bool TF;
  TF = setDIV(pset, idx, warn?true:false);
  if(!TF)   
    {
    if(warn)
      {
      if(warn > 1) IDVfatal(19);	//      Its a fatal error 
      else         IDVerror(19);	//	or a warning issued
      }
    return false;
    }
  return TF;
  }


// ----------------------------------------------------------------------------
//       Interactive Read of Dipole Interaction Vector From An ASCII File
// ----------------------------------------------------------------------------


string IntDipVec::ask_read(int argc, char* argv[], int argn)

	// Input		IDV    : Dipole interaction vector (this)
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
   "\n\tDipolar Interaction Vector Filename? ", 
                                      filename);
  read(filename);		           	// Read system from filename
  return filename;
  }

// ____________________________________________________________________________ 
// E               DIPOLE INTERACTION VECTOR OUTPUT FUNCTIONS
// ____________________________________________________________________________ 



	// Input		IDV	: Dipole interaction vector (this)
	// 			ostr	: Output stream
	//			full    : Flag for long vs short output
	// Output		non	: Dipole interaction vector parameters
	//			          sent to the output stream

ostream& IntDipVec::print(ostream& ostr, bool full) const
  {
  string hdr="Dipolar Interactions Vector";	// Title for output
  string ctr;					// Spacer to center output 
  if(!size()) 					// Exit if no interactions
    { 
    hdr += string(": NULL");
    ctr=string(40-hdr.length()/2, ' ');		// Set spacer to center output 
    ostr << "\n" << ctr << hdr << "\n";		// Output the header
    return ostr;
    }
  ctr = string(40-hdr.length()/2, ' ');		// Set spacer to center output 
  ostr << "\n" << ctr << hdr << "\n";		// Output the header
  int nd = size();
  if(!nonzero())				// If all interactions zero
    {						// just out the number we have
    hdr = Gdec(nd)				//	Make a new header
        + string(" Zero Interactions");
    ctr=string(40-hdr.length()/2, ' ');		// 	Spacer to center output 
    ostr << "\n\n" << ctr << hdr;		//	Output no nonzero ones
    return ostr;				//	Exit
    }
  for(int j=0; j<nd; j++)			// Output dipolar interactions
    {						// in order
    hdr = string("Interaction ") + Gdec(j);	//	Make a new header
    if(!((*this)[j].DCC())) hdr += ": Empty";
    ctr=string(40-hdr.length()/2, ' ');		// 	Spacer to center output 
    ostr << "\n" << ctr << hdr;			//	Output a header
    if(((*this)[j].DCC()))
      ((*this)[j]).print(ostr, full, false);	//	Output the interaction
    ostr << "\n";				// 	Add a line spacer
    }
  return ostr;
  }


ostream& operator<< (ostream& out, const IntDipVec& IDV)

	// Input		IDV	: Dipole interaction vector
	// 			out	: Output stream
	// Output		none	: Dipole interaction vector
	//				  parameters sent to output stream

  { return IDV.print(out); }


#endif
