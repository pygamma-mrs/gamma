/* IntRank2.cc **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Generic Rank 2 Interaction 		  Implementation	**
**                                                                      **
**      Scott Smith                                                     **
**      Copyright (c) 1997                                              **
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
** This class is the base class for irreducible rank 2 interactions in	**
** GAMMA. Primarily, it is a generic spin or spin-spin interaction that	**
** is composed of 3 elements: { Xi, A, T } where Xi is an interaction	**
** constant (for overall scaling(, A an irreducible rank 2 spatial	**
** tensor (class IntRank2A), and T an irreducible rank 2 spin tensor. 	**
**                                                                      **
** It also contains linked lists of irreducible rank 2 spin tensors, 	**
** applicable to all the interactions common to magnetic resonance.     **
** These lists allow GAMMA to avoid repeat generation of spin tensors,  **
** a likely occurance in the treatment of multi-spin systems because    **
** the spin tensors depend only on the Iz values of spins involved and  **
** the interaction type (in the single spin or spin pair Hilbert space) **
**                                                                      **
** Note that the use of linked lists doesn't avoid copying repeated     **
** spin tensors, but that is taken care of via referencing in GAMMA's   **
** matrix classes.                                                      **
**                                                                      **
** Finally, the class keeps track of a single set of Euler angles that	**
** specify the interaction orientation with respect to the laboratory   **
** frame.								**
**                                                                      **
*************************************************************************/

#ifndef   IntRank2_cc_			// Is file already included?
#  define IntRank2_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <IntRank2/IntRank2.h>		// Include interface defninition
#include <IntRank2/IntRank2A.h>		// Include interaction space tensors
#include <IntRank2/IntRank2T.h>		// Include interaction spin tensors
#include <list>				// Include libstdc++ STL lists
#include <Basics/Isotope.h>		// Know about isotopes
#include <Basics/Gutils.h>		// Know GAMMA error messaging
#include <Basics/Gconstants.h>		// Know DEG2RAD constant
#include <Basics/StringCut.h>		// StringCut as Gdec and Gform
#if defined(_MSC_VER) || defined(__SUNPRO_CC)				// If we are using MSVC++
 #pragma warning (disable : 4786)       //   Kill STL name length warnings
#endif	

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ STL lists
using std::vector;			// Using libstdc++ STL vectors
using std::ostream;			// Using libstdc++ output streams
using std::cout;			// Using libstdc++ standard output

typedef list<IntRank2T> IntRank2TList;	// Simplify using spin tensor lists

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                   RANK 2 INTERACTION ERROR HANDLING
// ____________________________________________________________________________
 

/*       Input                R2I     : Rank 2 interaction (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pn      : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void IntRank2::IR2error(int eidx, int noret) const
  {
  string hdr("Rank 2 Interaction");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr,"Program Aborting.....",        noret); break;// (0)
    case 1: GAMMAerror(hdr,"Problems During Construction", noret); break;// (1)
    case 8: GAMMAerror(hdr,"Theta Outside Range [0,180]",  noret); break;// (8)
    case 9: GAMMAerror(hdr,"Phi Outside Range [0,360]",    noret); break;// (9)
    case 10:GAMMAerror(hdr,"Asymmetry Outside Range [0,1]",noret); break;// (10)
    case 15:GAMMAerror(hdr,"Setting Asymmetry To Zero",    noret); break;// (15)
    case 20:GAMMAerror(hdr,"Setting Alpha Angle To Zero",  noret); break;// (20)
    case 21:GAMMAerror(hdr,"Setting Beta Angle To Zero",   noret); break;// (21)
    case 22:GAMMAerror(hdr,"Setting Gamma Angle To Zero",  noret); break;// (21)
    case 23:GAMMAerror(hdr,"Setting Theta Angle To Zero",  noret); break;// (23)
    case 24:GAMMAerror(hdr,"Setting Phi Angle To Zero",    noret); break;// (24)
    case 30:GAMMAerror(hdr,"Cannot Set External Field",    noret); break;// (30)
    case 31:GAMMAerror(hdr,"Cannot Set Larmor Frequency",  noret); break;// (31)
    case 32:GAMMAerror(hdr,"Cannot Set Base Frequency",    noret); break;// (32)
    default:GAMMAerror(hdr,"Unknown Error",                noret); break;
    }
  }  
 

volatile void IntRank2::IR2fatal(int eidx) const
  {
  IR2error(eidx, eidx);				// Output error message
  if(eidx) IR2error(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }
 
/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Funciton pname
      5     Please Use Class Member Funciton pname
   default  Unknown Error - pname                                           */


void IntRank2::IR2error(int eidx, const string& pname, int noret) const
  {
  string hdr("Rank 2 Interaction");
  string msg;
  switch(eidx)
    {
    case 110:                                                   // (110)
      msg = string("Unknown Isotope Type Used - ")  + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 111:                                                   // (111)
      msg = string("Cannot Parse Parameter ") + pname;
      GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;
    }
  }


// ____________________________________________________________________________
// ii                INTERACTION SPIN TENSOR SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions dictate whether interaction spin tensors are generated from
   scratch (and then stored in a linked list) or just copied from an existing
   spin tensor previously stored in a linked list.  Note that these have
   nothing to do with the interaction constant, spatial tensor, or orientation.

   These are PRIVATE because the value(s) of Ival and Sval of the interaction
   (IntRank2T) MUST be set prior to the use of these functions & their 
   misuse would be a bad thing.                                              */

// ----------------------------------------------------------------------------
//             Set Up Spin Tensor For A Single Spin Interaction 
// ----------------------------------------------------------------------------
 
/* We set the 5 spherical spin tensor components applicable to interactions
   involving a single spin interacting with a static magnetic field. Two such
   cases would be shift anisotropy (NMR) & g electron (EPR) interactions. Since
   this class maintains "scaled" spin tensors, all single spin - field
   interactions will share the same spin tensors. The 5 spin tensor components
   (arrays) will only differ between single spin interactions if the spin
   quantum value differs (or if the field orientation is set differently). 
   Typically, the static magnetic field is taken to be a normalized vector
   pointing along +z.  This class maintains a static linked list of such
   spin tensors to avoid regenerating them during computations involving
   multi-spin systems.
                                          +
                                       m  |
                            T    = (-1)  T
                             2,m          2,-m

                   1/2
                [4]                      1
          T   = |-| * I         T    = - - I           T    = 0
           2,0  [6]    z         2,1     2  +           2,2

                                  SINGLE SPIN - FIELD INTERACTIONS
           Input        R2I     : Rank 2 interaction (this)
           Output       none    : Interaction Spin Tensor spherical
                                  spin components are generated
	   Note			: Dont Forget To Set Ival & Sval 1st
           Note                 : For G interactions Ival is usually 2
                                  (Ival = 2*0.5 + 1, where I = 1/2 e-)
           Note                 : SA & G spin tensors are stored in the
                                  linked list SPFlist.  Thus, herein
                                  the spin tensor components are either
                                  constructed & stored, or copied from
                                  the spin tensor in the linked list.        */

void IntRank2::setSPF()
  {
  IntRank2TList::iterator item, prev;		// Item in spin tensor list
  item = SPFlist.begin();			// Start of SA & G tensor list
  int ihs;                                      // For list Hilbert spaces
  while(item != SPFlist.end())			// Loop through list entries
    {                                           // & try and find match
    ihs = (*item).IV();				// Entry Hilbert space
    if(ihs == Ival)				// If match just copy
      { 					// the spin tensor and get
      IntRank2T::operator= (*item);		// back to business
      return;
      }
    else if(ihs < Ival)                         // Not matched, goes in later
      {                                         // down in the list.  Lets
      prev = item;				// store where we're at then
      item++;					// and move to next list entry
      }
    else                                        // Here we must make a new
      {                                         // spin tensor.  So construct
      IntRank2T::setSPF();			// single spin interaction
      SPFlist.insert(item, *this);		// just one compared to
      return;
      }
    }
  IntRank2T::setSPF();                          // We will only reach this if
  SPFlist.push_back(*this);			// no list or at list end! So
  return;                                       // new interaction & append
  } 

// ----------------------------------------------------------------------------
//            Set Up Single Spin Interaciton Spin Tensor Components
// ----------------------------------------------------------------------------

/* Here we set the 5 spherical spin tensor components applicable to single spin
   interactions that are "interacting with themselves" in the sense that their
   own state modifies their environment which they interact. An example of such
   an interaciton would bea a quadrupolar interactions. However since this 
   class maintains "scaled" spin tensors, all single spin "self" interactions
   will share the same spin tensors. The five arrays will only differ between
   single spin "self" interactions if the spin quantum value differs.  This 
   class maintains a static linked list of such spin tensors to avoid
   regenerating them during computations involving multi-spin systems.

                                             +
                                          m  |
                               T    = (-1)  T
                                2,m          2,-m
           1/2
        [1]      2                      1                          1  2
  T   = |-| * [3I - I(I+1)]     T   = - -(I I + I I )       T    = - I
   2,0  [6]      z               2,1    2  + z   z +         2,2   2  +

                                            SINGLE SPIN INTERACTION 
          Input                R2I     : Rank 2 interaction (this)
          Output               none    : Interaction Spin Tensor spherical
                                         spin components are generated
	  Note			       : Dont Forget To Set Ival & Sval 1st
          Note                         : QUAD spin tensors are stored in the
                                         linked list SPQlist.  Thus, herein
                                         the spin tensor components are either
                                         constructed & stored, or copied from
                                         the spin tensor in the linked list.
          Note                         : For QUAD interactions Ival is always 
                                         greater than 2 (Ival=2*I+1, I>=1)   */
 
void IntRank2::setSPQ()
  {
  IntRank2TList::iterator item, prev;		// Item in spin tensor list
  item = SPQlist.begin();			// Start of Quad tensor list
  int ihs;                                      // For list Hilbert spaces
  while(item != SPQlist.end())			// Loop through list entries
    {                                           // & try and find match
    ihs = (*item).IV();				// Entry Hilbert space
    if(ihs == Ival)                             // If match just copy
      {                                         // the spin tensor and get
      IntRank2T::operator= (*item);		// back to business
      return;
      }  
    else if(ihs < Ival)                         // If not match, find where
      {                                         // it should be entered
      prev = item;
      item++;
      }
    else                                        // Here we must make a new
      {                                         // spin tensor.  So construct
      IntRank2T::setSPQ();			// single spin interaction
      SPQlist.insert(item, *this);		// and add it to the list
      return;
      }
    }
  IntRank2T::setSPQ();				// We will only reach this if
  SPQlist.push_back(*this);			// no list or at list end! So
  return;                                       // new interaction & append
  }  

// ----------------------------------------------------------------------------
//            Set Up Spin Pair Interaciton Spin Tensor Components
// ----------------------------------------------------------------------------

/* Here we set the 5 spherical spin tensor components applicable to spin-spin
   interactions. These can be akin to hyperfine interactions (elecron-nucleon)
   or dipolar interacions (electron-electron and nucleon-nucleon).  Since this
   class maintains "scaled" spin tensors, all spin-spin interactions share the
   same spin tensors. The five arrays will only differ between interactions if
   the spin quantum values differ.  This class maintains a static linked list
   of such spin tensors to avoid regenerating them during computations 
   involving multi-spin systems.
                                              +
                                           m  |
                                T    = (-1)  T
                                 2,m          2,-m
           1/2
        [1]                              1                            1
  T   = |-| * [3I S - I.S]       T   = - -(I S + I S )         T    = - I S
   2,0  [6]      z z              2,1    2  + z   z +           2,2   2  + +

                                             SPIN PAIR INTERACTIONS
          Input                R2I     : Rank 2 interaction (this)
	  Note			       : Dont Forget To Set Ival & Sval 1st
          Output               none    : Interaction Spin Tensor spherical
                                         spin components are generated       */

void IntRank2::setSPSP()
  {                                                                            
  IntRank2TList::iterator item, prev;		// Item in spin tensor list
  item = SPSPlist.begin();			// Start of Quad tensor list
  while(item != SPSPlist.end())			// Loop through list entries
    {                                           // and look for match
    if((*item).equalto(Ival,Sval))		// If we do find a match
      {                                         // just copy it and return
      IntRank2T::operator= (*item);
      return;
      }
    else if((*item).lessthan(Ival,Sval)) 	// If we don't find a match
      {                                         // we need to continue looking
      prev = item;				// until the list ends
      item++;
      }
    else                                        // Here we must make a new
      {                                         // entry.  So construct two
      IntRank2T::setSPSP();			// spin interaction spin tensor
      SPSPlist.insert(item, *this);		// list
      return;
      }
    }
  IntRank2T::setSPSP();				// We will only reach this if
  SPSPlist.push_back(*this);			// no list or at list end! So
  return;                                       // new interaction & append
  }  

// ____________________________________________________________________________
// iii            INTERACTION PARAMETER SET PARSING FUNCTIONS
// ____________________________________________________________________________

/* These functions glean information describing a generic rank 2 interaction
   from parameters in a GAMMA parameter set. Since base classes IntRank2A and
   IntRank2T handle irreducible rank 2 spatial and spin tensors respecitvely,
   we need only define the functions that deal with the reducible aspects
   of the rank 2 spatial components. This is done to support classes based on
   us! Even though we are irreducible, our derived interactions may utilize
   reducible rank 2 spatial tensors. We will also handle any additional
   parameter recognition that might be generically useful (e.g. Field).      */

// ----------------------------------------------------------------------------
//        Functions To Read Reducible Rank 2 Spherical Spatial Tensors
// ----------------------------------------------------------------------------

bool IntRank2::getAiAzAe(const ParameterSet& pset, const string& A,
           coord& Aize, int i, int j, bool warni, bool warnz, bool warne) const
  {
  double aiso  = 0;				// For isotropic value
  double adelz = 0; 				// For anisotropic value
  double aeta  = 0;				// For asymmetry value
  bool TFI = getAiso(pset,  A,aiso, i,j,warni);	// Try for Aiso first
  bool TFZ = getAaniso(pset,A,adelz,i,j,warnz);	// Try for Adelzz next
  bool TFE = getAeta(pset,  A,aeta, i,j,warne);	// Try for Aeta third
  Aize = coord(aiso, adelz, aeta);		// Set return coordinate
  return TFI&TFZ&TFE;				// Return how we did
  }

bool IntRank2::getAiso(const ParameterSet& pset, const string& A,
                             double& Aiso, int idxI, int idxS, bool warn) const
  {
  string Nidx("");				// Parameter name suffix
  if(idxI >= 0)					// Suffix only if idxI > -1
    {						//   Add idxI index to name
    Nidx  = string("(") + Gdec(idxI);	        //   as (#
    if(idxS > 0) 				//   If idxS > 0 add it to
      Nidx += string(",")+Gdec(idxS); 		//   to suffix as ,# (#,#
    Nidx += string(")");			//   Finish suffix (#) or (#,#)
    }
  string pname = A + "iso" + Nidx;		// String for parameter name
  string pstate;				// Temp string for parsing
  ParameterSet::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Seek parameter in pset
  if(item == pset.end()) 			// If not found, warn if
    { 						// desired, return fail
    if(warn) IR2error(111,pname,1);		// after zeroing the value
    Aiso = 0;					
    return false;
    }
  (*item).parse(pname,Aiso,pstate);		// If found, parse/set Aiso 
  return true;					// Return we found it OK
  }


bool IntRank2::getAaniso(const ParameterSet& pset, const string& A,
                            double& Adelz, int idxI, int idxS, bool warn) const
  {
  string Nidx("");				// Parameter name suffix
  if(idxI >= 0)					// Suffix only if idxI > -1
    {						//   Add idxI index to name
    Nidx  = string("(") + Gdec(idxI);	        //   as (#
    if(idxS > 0) 				//   If idxS > 0 add it to
      Nidx += string(",")+Gdec(idxS); 		//   to suffix as ,# (#,#
    Nidx += string(")");			//   Finish suffix (#) or (#,#)
    }
  string pname = A + "delz" + Nidx;		// String for parameter name
  string pstate;				// Temp string for parsing
  ParameterSet::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);			// Seek parameter in pset
  if(item != pset.end()) 			// If found, we simply parse
    {						// out the value & return
    (*item).parse(pname,Adelz,pstate);
    return true;
    }
  pname = A + "delzz" + Nidx;			// String for parameter name
  item = pset.seek(pname);			// Seek parameter in pset
  if(item != pset.end()) 			// If found, we simply parse
    {						// out the value & return
    (*item).parse(pname,Adelz,pstate);
    return true;
    }
  pname = A + "A" + Nidx;			// String for parameter name
  item = pset.seek(pname);			// Seek parameter in pset
  if(item != pset.end()) 			// If found, we parse out
    {						// the value, scale to delzz
    (*item).parse(pname,Adelz,pstate); 		// & then return
    Adelz *= 2.0/3.0;				//    anis = 2/3 delzz
    return true;
    }
  if(warn)					// Try as we might, we cannot
    {						// find this parameter in pset
    pname = A + "delz" + Nidx;			//   1st parameter name
    IR2error(111,pname,1);			//   We cannot find this one
    pname = A + "delzz" + Nidx;			//   2nd parameter name
    IR2error(111,pname,1);			//   We cannot find this one
    pname = A + "A" + Nidx;			//   3rd parameter name
    IR2error(111,pname,1);			//   We cannot find this one
    }
  return false;					// Return we failed to find it
  }

// ----------------------------------------------------------------------------
//     Functions To Read External Field Strengths, Larmor/Base Frequencies
// ----------------------------------------------------------------------------

/* These may be required for setting up field dependent interactions such
   as shift anisotropy and electron G interactions. The field is returned in
   Gauss and the only two parameters are Field and FieldT. The Larmor frequency
   (NMR) is specified by the lone parameter Omega which will be in MHz. The
   base frequency (EPR) is specified by the lone parameter GOmega which will
   be in in GHz. 

   Allowance are made so that all of thesse parameters can take a single 
   index, but normally they should not have one.                            */

bool IntRank2::getField(const ParameterSet& pset, const string& A,
                                         double& Bo, int idx, bool warn) const
  {
  Bo = 0;					// Default field strength
  string pstate; 	                        // Temp string for parsing
  ParameterSet::const_iterator item;         // A pix into parameter list

  string Nidx = "";                             // Name addition per index
  if(idx >= 0) Nidx                             // If index exists then set up
    += string("(") + Gdec(idx) + string(")");	// the parameter name to append

  string pname = A + "Field" + Nidx;		// Try for field strength (G)
  item = pset.seek(pname);                      // Pix in param. list for Field
  if(item != pset.end())                        // If Field has been defined
    {                                           // we'll parse the info and
    (*item).parse(pname,Bo,pstate);		// set field appropriately
    Bo = fabs(Bo);				// (this is in Gauss)
    return true;                                // Exit, we found Field
    }

  string pname2 = A + "FieldT" + Nidx; 		// Maybe field set in Tesla?
  item = pset.seek(pname2);			// Pix in param. list for FieldT
  if(item != pset.end())                        // If Field has been defined
    {                                           // we'll parse the info and
    (*item).parse(pname2,Bo,pstate);		// set field appropriately
    Bo *= 1.e4;					// (stored in Gauss)
    Bo = fabs(Bo);				// Must be positive
    return true;                                // Exit, we found Field
    }

  if(warn)                                      // If it hasn't been found
    { 						// issue warnings if needed
    IR2error(111, pname, 1);			//   Can't parse "Field"
    IR2error(40, 1);				//   Can't set external field
    }
  return false;
  }

bool IntRank2::getOmega(const ParameterSet& pset, const string& A,
                                         double& Om, int idx, bool warn) const
  {
  Om = 0;					// Default frequency
  string pstate; 	                        // Temp string for parsing
  ParameterSet::const_iterator item;         // A pix into parameter list

  string Nidx = "";                             // Name addition per index
  if(idx >= 0) Nidx                             // If index exists then set up
    += string("(") + Gdec(idx) + string(")");	// the parameter name to append

  string pname = A + "Omega" + Nidx;		// Try for frequency (MHz)
  item = pset.seek(pname);                      // Pix in param. list for freq.
  if(item != pset.end())                        // If Om has been defined
    {                                           // we'll parse the info and
    (*item).parse(pname,Om,pstate);		// set frequency appropriately
    return true;                                // Exit, we found frequency
    }

  if(warn)                                      // If it hasn't been found
    { 						// issue warnings if needed
    IR2error(111, pname, 1);			//   Can't parse "Omega"
    IR2error(31, 1);				//   Can't set Larmor frequency
    }
  return false;
  }

bool IntRank2::getGOmega(const ParameterSet& pset, const string& A,
                                          double& Om, int idx, bool warn) const
  {
  Om = 0;					// Default frequency
  string pstate; 	                        // Temp string for parsing
  ParameterSet::const_iterator item;         // A pix into parameter list

  string Nidx = "";                             // Name addition per index
  if(idx >= 0) Nidx                             // If index exists then set up
    += string("(") + Gdec(idx) + string(")");	// the parameter name to append

  string pname = A + "GOmega" + Nidx;		// Try for frequency (GHz)
  item = pset.seek(pname);                      // Pix in param. list for freq.
  if(item != pset.end())                        // If Om has been defined
    {                                           // we'll parse the info and
    (*item).parse(pname,Om,pstate);		// set frequency appropriately
    return true;                                // Exit, we found frequency
    }

  if(warn)                                      // If it hasn't been found
    { 						// issue warnings if needed
    IR2error(111, pname, 1);			//   Can't parse "GOmega"
    IR2error(32, 1);				//   Can't set base frequency
    }
  return false;
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A               RANK 2 INTERACTION CONSTRUCTORS AND DESTRUCTOR
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

/* Remember, we don't often DO NOT use IntRank2T constructor directly because
   we can avoid repreat building of the spin tensor using our list of spin
   tensors. Of course, for that we must know what type the spin tensor is... */

IntRank2::IntRank2() : IntRank2A(), IntRank2T() { _XI = 0.0; }

IntRank2::IntRank2(const IntRank2 &R2I) : IntRank2A(R2I), IntRank2T(R2I)
  { _XI = R2I._XI; }
 
// ---------------------------------------------------------------------------- 
//   Direct Single Spin & Spin - Field Constructors Taking Lots 'O Parameters
// ---------------------------------------------------------------------------- 

/* We set up an interaction involving a single spin interacting with a static
   magnetic field or with its own local field. Cases of the former are shift
   anisotropy (NMR) & electron G (EPR) interactions. The latter type would
   include quadrupolar interactions. As this class maintains "scaled" spatial
   & spin tensors, all single spin interactions will share the same spin 
   tensors. The five arrays will only differ between single spin interactions
   if the spin quantum value differs (or if field orientation is set different
   in spin-field interactions). Typically, the static magnetic field is taken
   to be a normalized vector pointing along +z.  Similarly, all interations
   share the same formulation of spatial tensor components.  
 
        Input                R2I     : Rank 2 interaction (this)
                             IsoI    : Spin isotope type 
                         or  qn      : Quantum number (e.g. 1.5, 0.5,...)
                             X       : Interaciton constant
                             ttype   : Spin tensor type
                                        true  = Spin Field (SA, G)
                                        false = Spin Self (Quad)
                             warn    : Warning output label
        Output               none    : Interaction constructed 
	Note                         : We do NOT use the base class (IntRank2T)
                                       directly herein because we can often
                                       avoid their construction by using a
                                       static linked list of interacions
                                            (SPQlist and SPFlist)            */

IntRank2::IntRank2(const string& IsoI, double X, 
                     double E, const EAngles& EA, bool ttype) :IntRank2A(E, EA)
  {
  Isotope II(IsoI);                     // Form isotope of type IsoI
  double mz = II.qn();			// Get spin quantum value (0.5,1.5,...)
  Ival = int(2.0*mz+1);			// Set Ival to 2I+1 (IntRank2T)
  Sval = 0;				// Set Sval to 0    (IntRank2T)
  _XI = X;				// Set the interaction constant
  (ttype)?setSPF():setSPQ(); 		// Set spin-field or spin-self
  }

IntRank2::IntRank2(const Isotope& IsoI, double X,
                     double E, const EAngles& EA, bool ttype) :IntRank2A(E, EA)
  {
  double mz = IsoI.qn();		// Get spin quantum value (0.5,1.5,...)
  Ival = int(2.0*mz+1);			// Set Ival to 2I+1 (IntRank2T)
  Sval = 0;				// Set Sval to 0    (IntRank2T)
  _XI = X;				// Set the interaction constant
  (ttype)?setSPF():setSPQ(); 		// Set spin-field or spin-self
  }

IntRank2::IntRank2(const double qn, double X,
                     double E, const EAngles& EA, bool ttype) :IntRank2A(E, EA)
  {
  Ival = int(2.0*qn+1);			// Set Ival to 2I+1 (IntRank2T)
  Sval = 0;				// Set Sval to 0    (IntRank2T)
  _XI = X;				// Set the interaction constant
  (ttype)?setSPF():setSPQ(); 		// Set spin-field or spin-self
  } 					// Functions setS* are member!

IntRank2::IntRank2(const string& IsoI, double X,
     const coord& AxAyAz, const EAngles& EA, bool ttype) :IntRank2A(AxAyAz, EA)
  {
  Isotope II(IsoI);                     // Form isotope of type IsoI
  double mz = II.qn();			// Get spin quantum value (0.5,1.5,...)
  Ival = int(2.0*mz+1);			// Set Ival to 2I+1 (IntRank2T)
  Sval = 0;				// Set Sval to 0    (IntRank2T)
  _XI = X;				// Set the interaction constant
  (ttype)?setSPF():setSPQ(); 		// Set spin-field or spin-self
  }

IntRank2::IntRank2(const Isotope& IsoI, double X,
     const coord& AxAyAz, const EAngles& EA, bool ttype) :IntRank2A(AxAyAz, EA)
  {
  double mz = IsoI.qn();		// Get spin quantum value (0.5,1.5,...)
  Ival = int(2.0*mz+1);			// Set Ival to 2I+1 (IntRank2T)
  Sval = 0;				// Set Sval to 0    (IntRank2T)
  _XI = X;				// Set the interaction constant
  (ttype)?setSPF():setSPQ(); 		// Set spin-field or spin-self
  }

IntRank2::IntRank2(double qn, double X,
     const coord& AxAyAz, const EAngles& EA, bool ttype) :IntRank2A(AxAyAz, EA)
  {
  Ival = int(2.0*qn+1);			// Set Ival to 2I+1 (IntRank2T)
  Sval = 0;				// Set Sval to 0    (IntRank2T)
  _XI = X;				// Set the interaction constant
  (ttype)?setSPF():setSPQ(); 		// Set spin-field or spin-self
  }
 
// ---------------------------------------------------------------------------- 
//          Direct Two Spin Constructors That Need Lots 'O Parameters
// ---------------------------------------------------------------------------- 
 
/* These constructors apply to spin-spin interactions (dipolar & hyperfine).
   Valid interaction types are specified in IntRank2T, these are currently
   HF and DIP. We limit the interaction types to spin pair interactions, other
   types will produce fatal errors. Remember, we typically DO NOT want to use
   the spin tensor (IntRank2T) constructors, opting instead to try and obtain
   the tensor arrays from our linked list.

        Input                R2I     : Rank 2 interaction (this)
                             IsoI    : Spin I isotope type
                             IsoS    : Spin S isotope type
                       or
                             Iqn     : Spin I quantum value (0.5, 1.5,...)
                             Sqn     : Spin S quantum value
                             delzz   : Tensor delzz value (Hz)
                             eta     : Tensor asymmetry value (default 0)
			     X	     : Interaction type                      */

IntRank2::IntRank2(const string& IsoI, const string& IsoS,
                                      double Xi, double eta, const EAngles& EA)
  { *this = IntRank2(Isotope(IsoI),Isotope(IsoS), Xi, eta, EA); }

IntRank2::IntRank2(const Isotope& IsI,const Isotope& IsS,
                                      double Xi, double eta, const EAngles& EA)
  {
  double iz = IsI.qn();			// Get Iz value of IsoI
  double sz = IsS.qn();			// Get Iz value of IsoS
  *this = IntRank2(iz,sz,Xi,eta,EA);	// Use overloaded constructor
  }

IntRank2::IntRank2(double Iqn, double Sqn,
                                      double Xi, double eta, const EAngles& EA)
         : IntRank2A(eta, EA)
  {
  Ival = int(2.0*Iqn + 1);		// Set the I Hilbert space
  Sval = int(2.0*Sqn + 1);		// Set the S Hilbert space
  _XI = Xi;				// Set the interaction constant
  setSPSP();				// Set up DIP spin tensor
  }
/*
IntRank2::IntRank2(const string& IsoI, const string& IsoS,
                             double Xi, const coord& AxAyAz, const EAngles& EA)
  { 
	*this = IntRank2(IsoI, IsoS, Xi, AxAyAz, EA); 
	// Currently this would cause infinite recursion.
	}
*/

IntRank2::IntRank2(const Isotope& IsoI, const Isotope& IsoS,
                             double Xi, const coord& AxAyAz, const EAngles& EA)
  {
  double iz = IsoI.qn();		// Get Iz value of IsoI
  double sz = IsoS.qn();		// Get Iz value of IsoS
  *this = IntRank2(iz,sz,Xi,AxAyAz,EA);	// Use overloaded constructor
  }

IntRank2::IntRank2(double Iqn, double Sqn, 
                             double Xi, const coord& AxAyAz, const EAngles& EA)
         : IntRank2A(AxAyAz, EA)
  {
  Ival = int(2.0*Iqn + 1);		// Set the I Hilbert space
  Sval = int(2.0*Sqn + 1);		// Set the S Hilbert space
  _XI  = Xi;				// Set the interaction constant
  setSPSP();				// Set up DIP spin tensor
  }

// ----------------------------------------------------------------------------
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------

void IntRank2::operator= (const IntRank2 &R2I1)
  {
  IntRank2A::operator=(R2I1);		// Copy the space tensor
  IntRank2T::operator=(R2I1);		// Copy the spin tensor
  _XI   = R2I1._XI;			// Copy interaction constant
  }

IntRank2::~IntRank2 () {}
 
// ____________________________________________________________________________
// B                     SPATIAL TENSOR COMPONENT ACCESS
// ____________________________________________________________________________
 
//-----------------------------------------------------------------------------
//                             Anisotropy Access
//-----------------------------------------------------------------------------

/*                  ^
   The anisotropy, /_\ A, relates to the delzz value of A. Because GAMMA uses
   irreducible rank 2 spatial tensors as a base class, both the anisotropy and
   the delzz value are constant for all spatial tensors.

                       1/2
                  [ 5 ]                  ^      3
          del   = |---]  = 0.51503      /_\ A = - del   = 0.77255
             zz   [6PI]                         2    zz                      */


// static double IntRank2A::delzz();				INHERITED
// static double IntRank2A::delA();				INHERITED

//-----------------------------------------------------------------------------
//			       Asymmetry Access
//-----------------------------------------------------------------------------

// Note that eta spans [0, 1] and is defined to be (Axx-Ayy)/Azz where
// |Azz| >= |Ayy| >= |Axx| in the treatment contained herein.
 
// double IntRank2A::eta( ) const;				INHERITED
// void   IntRank2A::eta(double Eta);				INHERITED
 
//-----------------------------------------------------------------------------
//		      Spherical Tensor Component Access
//-----------------------------------------------------------------------------

/* These allow one to access the irreducible spherical elements at a specific
   orientation. Most functions are derived from the base class IntRank2A. The
   exception are the PAS functions relative to those taking no angles. Since
   IntRank2A only exists in the PAS, here we reset the default function to
   use our only orientation and relegate the functions that take no arguments
   in IntRank2A to be end in PAS.  Note that these all use GAMMA scaling
   which sets A2m to be rank 2 spherical harmonics at all orientations if
   there is no asymmetry.                                                    */

/*
complex IntRank2::A20PAS(  ) const;				INHERITED
complex IntRank2::A21PAS(  ) const;				INHERITED
complex IntRank2::A2m1PAS( ) const;				INHERITED
complex IntRank2::A22PAS(  ) const;				INHERITED
complex IntRank2::A2m2PAS( ) const;				INHERITED

complex IntRank2::A20(  ) const;				INHERITED
complex IntRank2::A21(  ) const;				INHERITED
complex IntRank2::A2m1( ) const;				INHERITED
complex IntRank2::A22(  ) const;				INHERITED
complex IntRank2::A2m2( ) const;				INHERITED

complex IntRank2A::A20(double  A, double B, double G) const;	INHERITED
complex IntRank2A::A21(double  A, double B, double G) const;	INHERITED
complex IntRank2A::A2m1(double A, double B, double G) const;	INHERITED
complex IntRank2A::A22(double  A, double B, double G) const;	INHERITED
complex IntRank2A::A2m2(double A, double B, double G) const;	INHERITED

complex IntRank2A::A20(const  EAngles& EA) const;               INHERITED
complex IntRank2A::A21(const  EAngles& EA) const;               INHERITED
complex IntRank2A::A2m1(const EAngles& EA) const;               INHERITED
complex IntRank2A::A22(const  EAngles& EA) const;               INHERITED
complex IntRank2A::A2m2(const EAngles& EA) const;               INHERITED

complex IntRank2A::A2m(int m)                            const; INHERITED
complex IntRank2A::A2m(int m,double A,double B,double G) const; INHERITED
complex IntRank2A::A2m(int m, const EAngles& EA)         const; INHERITED

IR2ASph IntRank2A::SphCmp()                            const;   INHERITED
IR2ASph IntRank2A::SphCmp(double A,double B, double G) const;   INHERITED
IR2ASph IntRank2A::SphCmp(const EAngles& EA)           const;   INHERITED    */
  
//-----------------------------------------------------------------------------
//		           Cartesian Tensor Component Access
//-----------------------------------------------------------------------------

/*
  double IntRank2A::Axx() const;				INHERITED
  double IntRank2A::Axx(double theta, double phi=0) const;	INHERITED
  double IntRank2A::Ayy() const;				INHERITED
  double IntRank2A::Ayy(double theta, double phi=0) const;	INHERITED
  double IntRank2A::Azz() const;				INHERITED
  double IntRank2A::Azz(double theta, double phi=0) const;	INHERITED
  double IntRank2A::Ayx() const;				INHERITED
  double IntRank2A::Axy() const;				INHERITED
  double IntRank2A::Ayx(double theta, double phi=0) const;	INHERITED
  double IntRank2A::Axy(double theta, double phi=0) const;	INHERITED
  double IntRank2A::Azx() const;				INHERITED
  double IntRank2A::Axz() const;				INHERITED
  double IntRank2A::Azx(double theta, double phi=0) const;	INHERITED
  double IntRank2A::Axz(double theta, double phi=0) const;	INHERITED
  double IntRank2A::Azy( ) const;				INHERITED
  double IntRank2A::Ayz( ) const;				INHERITED
  double IntRank2A::Azy(double theta, double phi=0) const;	INHERITED
  double IntRank2A::Ayz(double theta, double phi=0) const;	INHERITED
  row_vector IntRank2A::CartComps() const;			INHERITED
  row_vector IntRank2A::CartComps(double theta, double phi=0) const;        */

  
//-----------------------------------------------------------------------------
//                         Orientation Angle Access
//-----------------------------------------------------------------------------

/* There are 3 angles which orient the spatial tensor relative to the inter-
   actions own principal axis system (PAS).  These are the set of Euler angles
   { alpha, beta, gamma }. They corrspond to rotating the coordinate axes first
   about the z-axis by alpha followed by a rotation about the new x-axis by
   beta and lastly about the new z-axis by gamma. For a symmetric tensor the
   last rotation is of no consequnce and set to zero, and for powder averages
   on the first two angles are used to sum over all spactial orientations.
   In these two cases the angles alpha and beta are one and the same as the
   spherical coordinate angles phi and theta respectively. Theta is the angle
   down from the PAS z-axis & the angle phi which is over from the PAS x-axis.
   The ranges of the angles (internally imposed on GAMMA's Euler angles) are
   alpha,gamma,phi = [0,360] and beta,theta=[0,180].  These functions allow
   users to both obtain and set these angles. Setting any of the angles will
   effictively reorient the spatial tensor.                                  */

/*
double  IntRank2A::alpha()       const;                         INHERITED
double  IntRank2A::beta()        const;                         INHERITED
double  IntRank2A::gamma()       const;                         INHERITED
double  IntRank2A::phi()         const;                         INHERITED
double  IntRank2A::theta()       const;                         INHERITED
EAngles IntRank2A::orientation() const;                         INHERITED

void IntRank2A::alpha(double A);                                INHERITED
void IntRank2A::beta(double  B);                                INHERITED
void IntRank2A::gamma(double G);                                INHERITED
void IntRank2A::phi(double   P);                                INHERITED
void IntRank2A::theta(double T);                                INHERITED
void IntRank2A::orientation(const EAngles& EA);                 INHERITED
void IntRank2A::orientation(double A, double B, double G, bool deg=false);   */
 
//-----------------------------------------------------------------------------
//                            Auxiliary Functions
//-----------------------------------------------------------------------------

/*  Function                             Purpose
   -----------   --------------------------------------------------------------
   AisoDelzEta   Given 3 PAS Cartesian components of a rank 2 spatial tensor
                 the function will determine the isotropic, anisotropic, and
                 asymmetry terms.  These are returned in a single coordinate.
   SortAxAyAz    Give three Cartesian components Axx, Ayy, and Azz in any
                 arbitrary they will are sorted so that |Azz| >= |Ayy| >= |Axx|
   CheckEta      Insures that an input asymmetry is in the range [0,1]
   PAS           True if tensor is in its PAS, false if not.
   Symmetric     True if tensor is symmetric (eta is zero), false if not
   CartMx        Returns a 3x3 matrix representing the Cartesian tensor
   Spherical     True if there is no scaling factor (XI interaction constant)
   Isotropic     True if there is no scaling factor (XI interaction constant)

// static coord  IntRank2::AisoDelzEta(const coord& AxAyAz);
// static void   IntRank2::SortAxAyAz(double& x, double& y, double& z);
//        bool   IntRank2::CheckEta(double eta, bool warn=true) const;
//        int    IntRank2::Symmetric( ) const;
//        matrix IntRank2::CartMx(double scale=1.0) const;                  */
 
bool   IntRank2::Spherical( )     const { return _XI?false:true; }
bool   IntRank2::Isotropic( )     const { return _XI?false:true; }
matrix IntRank2::CartMx(bool scf) const
  { return scf?IntRank2A::CartMx(_XI):IntRank2A::CartMx();}

// ____________________________________________________________________________
// C                       SPIN TENSOR COMPONENT ACCESS
// ____________________________________________________________________________

/* These functions allow users to access values associated with the interaction
   spin tensor. Since this class is derived from the base class IntRank2T, 
   all such functionality will be inherited from that class.                 */
 
//-----------------------------------------------------------------------------
//                          Spin Quantum Value Access
//-----------------------------------------------------------------------------

/*
  double IntRank2T::Izval()          const;			INHERITED
  double IntRank2T::Szval()          const;			INHERITED
  int    IntRank2T::IV()             const;			INHERITED
  int    IntRank2T::SV()             const;			INHERITED
  int    IntRank2T::HS()             const;			INHERITED
  double IntRank2T::qn(bool i=false) const;			INHERITED    */

//-----------------------------------------------------------------------------
//                      Spin Tensor Component Access
//-----------------------------------------------------------------------------

/*
  matrix IntRank2T::T20()           const;      // Return T2,0  INHERITED
  matrix IntRank2T::T21()           const;      // Return T2,1  INHERITED
  matrix IntRank2T::T2m1()          const;      // Return T2,-1 INHERITED
  matrix IntRank2T::T22()           const;      // Return T2,2  INHERITED
  matrix IntRank2T::T2m2()          const;      // Return T2,-2 INHERITED
  matrix IntRank2T::T2m(int m)      const;      // Return T2,m  INHERITED
  matrix IntRank2T::Tcomp(int comp) const;      // Return T2,m  INHERITED 

  matrix IntRank2T::T20(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T21(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m1(     const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T22(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m2(     const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m(int m,const vector<int>& HSs,int i,int j=-1) const;  */

//-----------------------------------------------------------------------------
//                      Spin Tensor Output Functions 
//-----------------------------------------------------------------------------

/*string* IntRank2T::TStrings(int M) const;			INHERITED    */


// ____________________________________________________________________________
// D                   RANK 2 INTERACTION CONSTANT ACCESS
// ____________________________________________________________________________
 
        // Input                R2I     : Rank 2 interaction (this)
        // Output               delzz   : Interaction delzz value
        // Note                         : Get/Set interaction delzz value
  
//double IntRank2::delz() const     { return DELZZ; }
//double IntRank2::delzz() const    { return DELZZ; }
//void   IntRank2::delz(double dz)  { DELZZ = dz; }
//void   IntRank2::delzz(double dz) { DELZZ = dz; }

double IntRank2::xi()            const { return _XI; }	// In radians/sec
void   IntRank2::xi(double xval)       { _XI = xval; }	// In radians/sec

// ____________________________________________________________________________
// E                          OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//     Functions That Generate Strings to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

string IntRank2::XiString() const
  {
  string XS("Xi Value");
  string SU;
  double xval = _XI;
       if(fabs(_XI) > 1.e6) { SU=string(" (x 10^-6):"); xval*=1.e-6; }
  else if(fabs(_XI) > 1.e3) { SU=string(" (x 10^-3):"); xval*=1.e-3; }
  else if(fabs(_XI) > 1.0)  { SU=string(":          ");              }
  else if(fabs(_XI) == 0.0) { SU=string(":          ");              }
  else if(fabs(_XI) < 1.e3) { SU=string(" (x 10^3): "); xval*=1.e3;  }
  else if(fabs(_XI) < 1.e6) { SU=string(" (x 10^6): "); xval*=1.e6;  }
  XS += SU + Gform("%10.2f", xval) + string(" /sec");
  return XS;
  }

//-----------------------------------------------------------------------------
//  Functions That Generate String Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

        // Input                R2I     : Rank 2 interaction (this)
        // Output               CSS     : Pointer to array of 5 strings
        // Note                         : The String array must be deleted
        //                                outside of this routine!
	// Note				: This outputs the interaction
	//				  values in generic fashion which is
	//				  not done in class IntRank2A because
	//				  we include spin quantum values

/*               Spin Quantum Number:       I
                 delzz:			xxxxx.xx ey
                 Asymmetry:                 x.xx                             */

vector<string> IntRank2::IR2AStrings() const
  {
  vector<string> SphStrings;
  if(Sval) SphStrings.push_back(StringIS());	// Spin Quantum Number(s)
  else     SphStrings.push_back(StringI());
  SphStrings.push_back(XiString());		// Interaction Constant
  SphStrings.push_back(AsymmetryString());	// Asymmetry
  return SphStrings;
  }
 

// string* IntRank2T::TStrings(int M) const	INHERITED

/*   				[ x.x, x.x, x.x]
  		        T     = [ x.x, x.x, x.x]
  		    	 2,m	[ x.x, x.x, x.x]
  				[ x.x, x.x, x.x]

  where M = { 0, 1, ..., 4 } ===> m = { 0, 1, -1, 2, -2 }                    */


ostream& IntRank2::printAT(ostream& ostr, IST_type ISTT) const
         
        // Input                R2I     : Rank 2 interaction (this)
        //                      ostr	: Output stream
	//			ISTT    : Interaction type
        // Output               none    : Interaction spatial tensor
        //                                sent to the output stream
	// Note				: This serves as a base print function
	//				  that the other IntName interaction
	//				  derived classes use.

/*	    Prints The Spherical Components, Spatial And Spin
       These Will Be Printed Lines Which Will Look Like The Following 
  		      (Repeated For All 5 m Values)
  
   						[ x.x, x.x, x.x]
  		A    = x.xxx		T     = [ x.x, x.x, x.x]
  		 2,m			 2,m	[ x.x, x.x, x.x]
  						[ x.x, x.x, x.x]             */

  {
  if(Izval() < 0.5)				// Just exit if nothing
    {
    string hdr("Empty ");
    switch(ISTT)
      {
      case SA:
        hdr += string("Shift Anisotropy");
        break;
      case G:
        hdr += string("Electron G");
        break;
      case QUAD:
        hdr += string("Quadrupolar");
        break;
      case DIP:
        hdr += string("Dipolar");
        break;
      case HF:
        hdr += string("Hyperfine");
        break;
      default:
        hdr += string("Rank 2");
        break;
      }
    hdr += string(" Interaction");
    int hl = hdr.length();
    string spacer = string(40-hl/2, ' ');
    ostr << "\n\n" << spacer << hdr << "\n";
    return ostr;
    }

  vector<string> DTS;				// For T2m strings
  vector<string> ASP = SphA2Strings();		// Get A2m strings
  int i, ll, hs = HS();				// Spin pair Hilbert space
  string spacer("  ");				// Space Between A & T
  string Asp = string(ASP[0].length(), ' ');	// Space for missing A part
  string Cen;					// Space to center on 80 cols 
  for(int m=0; m<5; m++)			// Loop 5 m components
    {
    DTS = TStrings(m);				// Get T2m array as strings
    ostr << "\n";				// Add a line break
    ll = Asp.length()				// Line lengths this m value
       + spacer.length() + DTS[0].length();
    if(80-ll>0) Cen=string((80-ll)/2, ' ');	// Set spacer to center lines
    else        Cen="";
    for(i=0; i<hs; i++)
      {
      ostr << "\n" << Cen;
      if(i == hs/2-1)    ostr << ASP[2*m];
      else if(i == hs/2) ostr << ASP[2*m+1];
      else               ostr << Asp;
      ostr << spacer << DTS[i];
      }
    }						// Do next m component
  return ostr;
  }


//-----------------------------------------------------------------------------
//   Functions That Generate Ouput Of The Rank 2 Shift Anistropy Interaction
//-----------------------------------------------------------------------------


// sosi seems to be problems here? cout << print leave some garbage....
ostream& IntRank2::print(ostream& ostr, int fflag) const

        // Input                R2I     : Rank 2 interaction (this)
        //                      ostr    : Output stream
        //                      fflag   : Format flag
        //                                  0 - Basic Parameters
        //                                 !0 - Full output
        // Output               none    : Interaction parameters
        //                                placed into the output stream
        // Note                         : This does NOT use the base class
        //                                virtual overload because we write
        //                                out two spin quantum values here?

  {
  if(Ival < 2.0)                             // Just exit if nothing
    {
    ostr << "\n\t\tEmpty Generic Irreducible Rank 2 Interaction\n";
    return ostr;
    }
  string hdr="Irreducible Rank 2 Interaction";	// Use this header
  vector<string> SS= IR2AStrings();		// Array of info strings
  IntRank2A::print(ostr, hdr, SS);		// Print spatial tensor
  if(fflag) printAT(ostr);                      // Print sphercial A & T
  ostr << "\n\n";                               // Add some linefeeds
  return ostr;
  }          

    
ostream& operator<< (ostream& out,const IntRank2& R2I) {return R2I.print(out);}
 
        // Input                out     : Output stream;
        // 			R2I     : Rank 2 interaction
        // Output                       : Modifies output stream
 
// ____________________________________________________________________________
// F                    SPIN TENSOR LINKED LIST ACCESS
// ____________________________________________________________________________

ostream& IntRank2::printISLList(ostream& ostr, const list<IntRank2T>& ISL)
  {
  string title;
  int nst = ISL.size();
  if(!nst)
    {
    title = string("List Is Currently Empty");
    ostr << "\n\n" << string((80-title.length())/2, ' ') << title;
    return ostr;
    }

  title = string("(") + Gdec(nst) + string(" Stored Spin Tensor");
  if(nst>1) title += "s)";
  else      title += ")";
  ostr << string((80-title.length())/2, ' ') << title << "\n";

  IntRank2TList::const_iterator item;		// Item in spin tensor list
  item = ISL.begin();				// Get 1st spin tensor
  int count = 0;				// Keep count as we output
  bool singsp = true;				// Assume single spin
  if((*item).Szval()) singsp = false;		// Set spin pair if Sz value

  if(singsp)
    {
    string spr(20, ' ');
    ostr << "\n" << spr << "Spin Tensor   Spin Quantum  Hilbert Space";
    ostr << "\n" << spr << "   Index         Number       Dimension  ";
    ostr << "\n" << spr << "-----------   ------------  -------------";
    while(item != ISL.end())			// Loop over all in the list
      {
      count++;
      ostr << "\n";
      ostr << spr      << Gdec(count,6)                  << string(8, ' ');
      ostr << "    "  << Gform("%4.1f", (*item).Izval()) << string(7, ' ');
      ostr << "    "  << Gdec((*item).HS(), 2);
      item++;
      }
    }
  else
    {
    while(item != ISL.end())			// Loop over all in the list
      {
      count++;
      ostr << "\n       " << Gdec(count,3) << ". ";
      ostr << "I = "   << Gform("%4.1f", (*item).Izval());
      ostr << ", S = "   << Gform("%4.1f", (*item).Szval())  << ", ";
      ostr << ", IV = "  << (*item).IV()     << ", ";
      ostr << ", SV = "  << (*item).SV()     << ", ";
      ostr << ", HS = "  << (*item).HS();
      item++;
      }
    }
  return ostr;
  }

ostream& IntRank2::printList(ostream& ostr, bool fflag)
  {
  ostr << "\n\t\t# of Stored G/Shift Anisotropy Spin Tensors: " << SPFlist.size();
  ostr << "\n\t\t# of Stored Quadrupolar Spin Tensors:        " << SPQlist.size();
  ostr << "\n\t\t# of Stored Hyperfine/Dipolar Tensors:       " << SPSPlist.size();
  if(fflag)
    {
    ostr << "\n\n";
    string hdr(" Rank 2 Interaction Spin Tensor List");
    if(SPFlist.size())
      {
      string title("CSA/Electron G");
      title += hdr;
      ostr << string((80-title.length())/2, ' ')
           << title << "\n";
      printISLList(ostr, SPFlist);
      }
    if(SPQlist.size())
      {
      string title("Quadrupolar");
      title += hdr;
      ostr << string((80-title.length())/2, ' ')
           << title << "\n";
      printISLList(ostr, SPQlist);
      }
    if(SPSPlist.size()) 
      {
      string title("Dipolar/Hyperfine");
      title += hdr;
      ostr << string((80-title.length())/2, ' ')
           << title << "\n";
      printISLList(ostr, SPSPlist);
      }
    ostr << "\n\n";
    }
  return ostr;
  }


ostream& IntRank2::printList(ostream& ostr, IST_type X, int fflag)
  {
  string hdr("");			// For header if any
  string spc("");			// Spacer to center header
  string lf("");			// Linefeed after header
  int    nst=0;				// Number of spin tensors
  if(fflag > 1)				// Set header if desired
    {
    switch(X)
      {
      case G:
      case SA:   hdr = string("CSA/Electron G");    break;
      case QUAD: hdr = string("Quadrupolar");       break;
      case HF:
      case DIP:  hdr = string("Dipolar/Hyperfine"); break;
      default:   hdr = string("Type Is Unknown");   break;
      }
    hdr += string(" Rank 2 Interaction Spin Tensor List");
    spc =  string((80-hdr.length())/2, ' ');
    lf  =  string("\n");
    }
  else if(!fflag)
    {
    hdr = string("\n\t\t# of Stored ");
    switch(X)
      {
      case G:    hdr += string("Electron G");       nst = SPFlist.size();  break;
      case SA:   hdr += string("Shift Anisotropy"); nst = SPFlist.size();  break;
      case QUAD: hdr += string("Quadrupolar");      nst = SPQlist.size();  break;
      case HF:   hdr += string("Hyperfine");        nst = SPSPlist.size(); break;
      case DIP:  hdr += string("Dipolar");          nst = SPSPlist.size(); break;
      default:   hdr += string("Unknown");          nst = 0;               break;
      }
    hdr += string("Spin Tensors: ");
    }
  switch(X)
    {
    case G:
    case SA:
      if(fflag) { ostr << spc << hdr << lf; printISLList(ostr, SPFlist); }
      else      { ostr << hdr << nst; }
      break;
    case QUAD:
      if(fflag) { ostr << spc << hdr << lf; printISLList(ostr, SPQlist); }
      else      { ostr << hdr << nst; }
      break;
    case HF:
    case DIP:
      if(fflag) { ostr << spc << hdr << lf; printISLList(ostr, SPSPlist); }
      else      { ostr << hdr << nst; }
      break;
    default:
    case UNK:
      ostr << "\n\t\tNo Clue What This Spin Tensor Type Is, Sorry.";
      break;
    }
  return ostr;
  }


// ____________________________________________________________________________
// G                RANK 2 INTERACTION COMPARISON FUNCTIONS
// ____________________________________________________________________________


int IntRank2::operator==(const IntRank2 &R2I2) const

        // Input                R2I     : Rank 2 interaction (this)
        //                      R2I2	: Another rank 2 interaction
        // Output               T/F     : TRUE if R2I2 and this are equal

  {
  double cutoff = 1.e-9;
if(IntRank2T::operator!=(R2I2)) 
  {
  cout << "\n\tFailed Spin Tensor Comparison.";
  return 0;
  }
if(IntRank2A::operator!=(R2I2))
  {
  cout << "\n\tFailed Space Tensor Comparison.";
 return 0;
  }
if(fabs(_XI-R2I2._XI) > cutoff)
  {
  cout << "\n\tFailed Xi Comparison.";
 return 0;
  }
/*
  if(IntRank2T::operator!=(R2I2)) return 0;
  if(IntRank2A::operator!=(R2I2)) return 0;
  if(fabs(I-R2I2.I) > cutoff) return 0;
  if(fabs(S-R2I2.S) > cutoff) return 0;
  if(fabs(_XI-R2I2._XI) > cutoff) return 0;
*/
  return 1;
  }


int IntRank2::operator!=(const IntRank2 &R2I2) const

        // Input                R2I     : Rank 2 interaction (this)
        //                      R2I2	: Another rank 2 interaction
        // Output               T/F     : FALSE if IR2A2 and this are equal

  {
  if(IntRank2T::operator!=(R2I2)) return 1;
  if(IntRank2A::operator!=(R2I2)) return 1;
  double cutoff = 1.e-11;
  if(fabs(_XI-R2I2._XI) > cutoff) return 1;
  return 0;
  }


// ____________________________________________________________________________
// H                RANK 2 INTERACTION HAMILTONIAN FUNCTIONS
// ____________________________________________________________________________

/* This section returns irreducible rank 2 interaction Hamiltonians. Because
   this class uses standardized spatial and spin tensors the Hamiltonians may
   all be derived from the same formula.

                        -2,2
                         ---           m                              
  H(alpha,beta,gamma) =  \   Xi  * (-1)  * A   (alpha,beta,gamma) * T    (i,j)
                         /                  2,m                      2,-m
                         ---
                          m
   
   The Hamiltonians will be returned in the single spin or spin pair Hilbert
   space of the interaction, unless a composite Hilbert space and some spin
   indices are supplie. All Hamiltonians will returned in the product basis 
   as simple matrices. Their units will be radians/sec.                      */


// ----------------------------------------------------------------------------
//                  First Order Interaction Hamiltonian
//       These Are SECULAR (Rotationally Invariant About Bo Field Axis)
//     Applicable When The Interaction Is A Small Perturbation To Zeeman
// ----------------------------------------------------------------------------

/* Note that for spin pair interactions this Hamiltonian is NOT invariant in a
   multiple rotating frame, but becomes time dependent! For example, in a
   heteronuclear dipolar interaction where the return is to be in the rotating
   frame of both I and S then H0 is time dependent from the I+S- and the I-S+
   terms in T20. If one chooses to work in the laboratory frame, i.e. add H0 to
   the lab-frame Zeeman Hamiltonian there will not be problems. Neither will
   there be problems if the interaction is single spin (electron G, CSA, Quad).
   But, to work in multiple rotating frames you must NOT use the flip-flop
   terms that occur in spin pair interactions. The latter assume a high-field
   limit! Individual (derived) spin pair interactions will likely supply the
   user with the ability to drop the flip-flop terms, that is NOT done here.

   The secular part of the interaction Hamiltonian is that which contains only
   those components which commute with z axis rotations.  Here we have

       [1]                                                           (0)
      H   (alpha,beta,gamma) = Xi * A   (alpha,beta,gamma) * T    = H
                                     2,0                      2,0

           Input                R2I     : Rank 2 interaction
                                alpha   : Euler angle (radians)
                                beta    : Euler angle (radians)
                                gamma   : Euler angle (radians)
                                EA      : Euler angles (radians)
                                HSs     : Array of spin Hilbert spaces
                                i       : First spin index (in HSs)
                                j       : Seond spin index (in HSs)
           Output               H       : Matrix for interacton Hamiltonian
           Note                         : No HSs: return in the spin Hilbert
                                          space of dimension (2I+1)[*(2S+1)]
                                          With HSs: return in composite spin
                                          Hilbert space                      */

matrix IntRank2::H0() const                    { return (_XI*A20())      * T0; }
matrix IntRank2::H0(double A, double B, double G) const
                                               { return (_XI*A20(A,B,G)) * T0; }
matrix IntRank2::H0(const EAngles& EA) const   { return (_XI*A20(EA))    * T0; }


matrix IntRank2::H0(const std::vector<int>& HSs, int i)			const
  { return CompSpace(H0(),HSs,i); }		// Return composite space H

matrix IntRank2::H0(const std::vector<int>& HSs, int i,
                                            double A, double B, double G) const
  { return CompSpace(H0(A,B,G),HSs,i); }	// Return composite space H

matrix IntRank2::H0(const std::vector<int>& HSs, int i, const EAngles& EA) const
  { return CompSpace(H0(EA),HSs,i); }		// Return composite space H


matrix IntRank2::H0(const std::vector<int>& HSs, int i, int j)            const
  { return CompSpace(H0(),HSs,i,j); }		// Return composite space H

matrix IntRank2::H0(const std::vector<int>& HSs, int i, int j,
                                            double A, double B, double G) const
  { return CompSpace(H0(A,B,G),HSs,i,j); }	// Return composite space H

matrix IntRank2::H0(const std::vector<int>& HSs, int i, int j,
                                                       const EAngles& EA) const
  { return CompSpace(H0(EA),HSs,i,j); }		// Return composite space H


// ----------------------------------------------------------------------------
//                          Full Interaction Hamiltonians
//                           These Use No Approximations! 
//   They Are Correct ln the Lab. Frame, But Probably NOT in a Rotating Frame
// ----------------------------------------------------------------------------

matrix IntRank2::H() const
  {
  IR2ASph As = SphCmp();			// Get spherical components
  matrix Hmx =  As.A20()*T0;			// Sum over the 5 components
         Hmx -= As.A2m1()*T1 + As.A21()*Tm1;	// using generic formulation
         Hmx += As.A2m2()*T2 + As.A22()*Tm2;
  Hmx.set_type(h_matrix_type);			// Insure it is Hermitian
  return _XI*Hmx;				// Scale and return 
  }						// Spin/SpinPair Hilbert Space

matrix IntRank2::H(double alpha, double beta, double gamma) const
  {
  IR2ASph As = SphCmp(alpha, beta, gamma);	// Get spherical components
  matrix Hmx =  As.A20()*T0;			// Sum over the 5 components
         Hmx -= As.A2m1()*T1 + As.A21()*Tm1;	// using generic formulation
         Hmx += As.A2m2()*T2 + As.A22()*Tm2;
  Hmx.set_type(h_matrix_type);			// Insure it is Hermitian
  return _XI*Hmx;				// Scale and return 
  }						// Spin/SpinPair Hilbert Space

matrix IntRank2::H(const EAngles& EA) const
  {
  IR2ASph As = SphCmp(EA);			// Get spherical components
  matrix Hmx =  As.A20()*T0;			// Sum over the 5 components
         Hmx -= As.A2m1()*T1 + As.A21()*Tm1;	// using generic formulation
         Hmx += As.A2m2()*T2 + As.A22()*Tm2;
  Hmx.set_type(h_matrix_type);			// Insure it is Hermitian
  return _XI*Hmx;				// Scale and return 
  }						// Spin/SpinPair Hilbert Space

matrix IntRank2::H(const std::vector<int>& HSs, int i, int j)            const
  { return CompSpace(H(),HSs,i,j); }		// Return composite space H

matrix IntRank2::H(const std::vector<int>& HSs, int i, int j,
                                           double A, double B, double G) const
  { return CompSpace(H(A,B,G),HSs,i,j); }	// Return composite space H

matrix IntRank2::H(const std::vector<int>& HSs, int i, int j,
                                                      const EAngles& EA) const
  { return CompSpace(H(EA),HSs,i,j); }		// Return composite space H


#endif							// IntRank2.cc


/* Final Comments:

    1. This class serves two purposes.  First, it is the base class for GAMMA's
       irreducible rank 2 spin interactions, hopefully making those classes
       quite simple in design.  Second it tracks three spin tensor linked lists
       soas to minimize storage of spin tensors.  The linked lists are based on
       the spin tensor class (IntRank2T) and thus CANNOT reside in spin tensor
       class itself (although that would perhaps be more logical.)
    2. This class (as IntRank2T) also has knowledge of isotopes so that we can
       have constructors that are intuitive in that they take spin types.
    3. Since we've taken the time to maintain spin tensor linked lists herein,
       all derived interaction classes derived should use the constructors of
       this class!                                                           */
 

