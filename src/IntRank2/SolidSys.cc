/* SolidSys.cc **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Solid Spin System                           Implementation	**
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
** A solid spin system is a collection of spin isotopes each of which   **
** has associated tensor properties (shift, quadrupole) and a relative  **
** position in 3-dimensional space (dipole).  Thus the solid system     **
** corresponds to a single molecule (or spin network, or crystallite)   **
** with a set orientation.                                              **
**                                                                      **
*************************************************************************/

#ifndef   Solid_sys_cc_			// Is file already included?
#  define Solid_sys_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#if defined(_MSC_VER)			// If we are using MSVC++
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif

#include <IntRank2/SolidSys.h>		// Include interface definition
#include <HSLib/SpinSystem.h>		// Includes spin system knowledge
#include <Basics/Gutils.h>		// Includes query knowledge
#include <Basics/StringCut.h>		// Includes Gdec knowledge
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/Isotope.h>		// Include spin isotopes
#include <IntRank2/IntDip.h>		// Include dipolar interactions
#include <IntRank2/IntCSA.h>		// Include shift anisotropy interacts
#include <IntRank2/IntQuad.h>		// Include quadrupolar interactions
#include <IntRank2/IntDipVec.h>		// Include dipolar interaction vectors
#include <IntRank2/IntQuadVec.h>	// Include quad. interaction vectors
#include <IntRank2/IntCSAVec.h>		// Include CSA interaction vectors
#include <Level1/coord.h>		// Include single coordinates
#include <Level1/coord_vec.h>		// Include coordinate vectors

using std::string;			// Using libstdc++ strings
using std::vector;			// Using libstdc++ STL vectors
using std::ofstream;			// Using libstdc++ output file streams
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                 CLASS SOLID SPIN SYSTEM ERROR HANDLING
// ____________________________________________________________________________
 
/*       Input 		      ssys    : Solid spin system (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void solid_sys::ssys_error(int eidx, int noret) const
  {
  string hdr("Class SolidSys");
  string m;
  switch(eidx)
    {
    case 0: GAMMAerror(hdr,"Program Aborting.....",        noret); break;// (0)
    case 1: GAMMAerror(hdr,"Can't Access Dipolar Interact",noret); break;// (1)
    case 2: GAMMAerror(hdr,"Electron Coordinate Specified",noret); break;// (2)
    case 3: GAMMAerror(hdr,"Spin Coord. List Incomplete",  noret); break;// (3)
    case 4: GAMMAerror(hdr,"If 1 Set, Set All Spin Coords",noret); break;// (4)
    case 11:GAMMAerror(hdr,"No Quad. Coupling For I=1/2",  noret); break;// (11)
    case 12:GAMMAerror(hdr,"Null Interact. Asymmetry Set", noret); break;// (12)
    case 13:GAMMAerror(hdr,"Access Of Self-Paired Dipole", noret); break;// (13)
    case 14:GAMMAerror(hdr,"Sorry, Dipolar Op Not Allowed",noret); break;// (14)
    case 20:GAMMAerror(hdr,"Can't Write To Parameter File",noret); break;// (20)
    case 21:GAMMAerror(hdr,"Can't Read From Param. File",  noret); break;// (21)
    case 22:GAMMAerror(hdr,"Can't Output To File Stream",  noret); break;// (22)
    case 23:GAMMAerror(hdr,"Cannot Output Parameters",     noret); break;// (23)
    case 24:GAMMAerror(hdr,"Cannot Set Up The System",     noret); break;// (24)
    case 30:GAMMAerror(hdr,"Can't Modify DCC=0 Dipole",    noret); break;// (30)
    case 31:GAMMAerror(hdr,"Can't Modify Electron Dipole", noret); break;// (31)
    case 50:m = string("Cannot Set Spin Coordinates, Need At Least ")
              + Gdec(spins()) + string(" Coordinates In Input Vector");
            GAMMAerror(hdr,m,                              noret); break;// (50)
    default:GAMMAerror(hdr,"Unknown Error",                noret); break;
    }
  }

volatile void solid_sys::ssys_fatal(int eidx) const
  {
  ssys_error(eidx, 1);			// First output the error
  if(eidx) ssys_error(0);		// Now output its fatal
  GAMMAfatal();					// Clean exit from program
  }


/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Funciton pname
      5     Please Use Class Member Funciton pname
   default  Unknown Error - pname                                           */


	// Input		ssys	: Solid spin system (this)
        // 			eidx	: Error index
        //                      pname	: Additional string for error
        //                      noret	: Flag for return (0=linefeed)
        // Output               none	: Error message

void solid_sys::ssys_error(int eidx, const string& pname, int noret) const
  {
  string hdr("Class SolidSys");
  string msg;
  switch(eidx)
    {
    case 102:							// (102)
      msg = string("Can't Get System, Index ")
          + pname + string(" From Parameter Set");
      GAMMAerror(hdr, msg, noret); break;
    case 103:							// (103)
      msg = string("Ignoring Coordinate Specification Of Spin ") 
          + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 104:							// (104)
      msg = string("Cannot Find Coordinate Specification For Spin ") 
          + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 105:							// (105)
      msg = string("Cannot Get/Set Coordinate Of Spin ") 
          + pname + string(", Outside Range [0,")
          + Gdec(spins()-1) + string("]");
      GAMMAerror(hdr, msg, noret); break;
    case 106:							// (106)
      msg = string("Cannot Set Coordinate Of An Electron, Spin ") 
          + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 107:							// (107)
      msg = string("Cannot Set Coordinate Of Spin ") 
          + pname
           + string(", Coordinates Do Not Exist For Any Spins");
      GAMMAerror(hdr, msg, noret); break;
    case 130:							// (130)
      msg = string("Parameter ") 
          + pname + string(" Is The Culprit!\n");
      GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;
    }
  }

volatile void solid_sys::ssys_fatal(int eidx, const string& pname) const
  {
  ssys_error(eidx, pname, 1);		// First output the error
  if(eidx) ssys_error(0);		// Now output its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                      SOLID SYSTEM SETUP FUNCTIONS
// ____________________________________________________________________________
 
// ---------------------------------------------------------------------------- 
//           Functions To Set Coordinates and Their Existence Flags
// ---------------------------------------------------------------------------- 

/* This spin system contains interactions defined over the system spins and
   spin pairs. These are contained in vectors and very simple functions are
   used deal with them. The one exception is dipolar interactions because we
   allow them to be set using spin coordinates.  The spin coordinates are 
   handled independently.... hence these functions.

   Normally, the coordinate vector could be read in directly using the class
   coord_vec.  But here we do this independently because we don't mind if
   some spins have no coordinates assigned, such as electrons.  So, here
   we'll read the coordinates in individually and just keep track of the 
   coordinates which haven't been specified.

   Note that in order to read the coordinates we at least need to know the 
   number of spins in our system.  Thus we keep this function private and 
   prevent it being called inappropriately.                                  */

	// Input		ssys	: Solid spin system (this)
	// 			pset	: A parameter set
        // Output               T/F	: True if spin coordiantes are found in
 	//				  the pset & they're read into ssys.
	// Note				: Assumes array cflags allocated


	// Input		ssys	: Solid spin system (this)
        // Output               none	: If there are spins in the system
 	//				  this zeros the coordinate existence
        // 				  flags.  If the flags don't exist,
	//				  then the array is constructed.  If
	// 				  no spins, the array is set to NULL.

void solid_sys::zero_cindx()
  {
  int ns = spins();				// Get number of spins
  cflags = vector<int>(ns,0);			// Fill with ns zeros
  }

int solid_sys::setCoords(const ParameterSet& pset)
  {
  int ccount = 0;			// Count coordinates read in
  int ns = spins();			// Number of system spins			
  SCoords = coord_vec(ns);		// Vector for coordinates
  coord pt;				// A working single point
  for(int i=0; i<ns; i++)		// Try and retrieve the
    { 					// coordinates for each spin
    if(pt.read(pset, i, 0))		// If coordinate is found
      {
      SCoords.put(pt,i);		//    Store the coordinate
      cflags[i] = 1;			//    Flag that we know it
      ccount++;				//    Track how many we know 
      }
    else cflags[i] = 0; 		// We can't find spin coord
    } 					// so flag we don't know it
  return ccount;
  }
 
// ---------------------------------------------------------------------------- 
//	             Functions To Set Dipolar Interactions
// ---------------------------------------------------------------------------- 

/* For this spin system we (should) know the number of spins and their isotope
   types.  The dipolar interactions are all to be placed into the dipolar
   interaction vector, Dvec.  This must be read in using spin (NOT interaction)
   indices.  For example, { Iso(i), coord(i), Iso(j), coord(j) } for all spin
   pairs.  However, not everyone likes to use spin coordinates for designating
   a dipolar interaction, so we must allow for their setting without the use
   of coordinates, but using { Isotopes, spin indices, + other info }

	   Input		ssys	: Solid spin system (this)
	   			pset	: A parameter set
	  			ccount  : Count of specified coordinates
	   Output		none	: System dipolar interactions are
	  			  	  set from parameters in pset
	   Note				: Assumes that
					    1.) System isotopes are set
	  				    2.) Any spin coordinates have
	  					been set (setCoords)
	  				    3.) All dipolar interactions are
	  				        specified by spin indices
	  				        (input formats can vary)     */
 
void solid_sys::setDs(const ParameterSet& pset, int ccount)
  {
  vector<Isotope> Isos = IsoVec();		// Get vector system isotopes
  if(ccount && cflags.size())			// If we know spin coordinates
    Dvec = IntDipVec(Isos, SCoords, 0);		// set interactions using them
  else						// If no spin coordinates
    Dvec = IntDipVec(pset,Isos,-1,0);		// set interactions using Isos
  }

// ---------------------------------------------------------------------------- 
//	          Functions To Set Shift Anisotropy Interactions
// ---------------------------------------------------------------------------- 

/* Each Spin That Isn't An Electron May Have A Shielding Anisotropy Interaction
   Set.  We Don't Particularly Care If There Isn't A Shielding Anisotropy 
   Interaction Defined For Any Nuclear Spin, The Vector Of Such Interactions
   Will Simply Have An Empty Value Set If No Anisotropy Has Been Specified.

	   Input		ssys	: Solid spin system (this)
	   			pset	: A parameter set
	   Output		none	: System shift anisotropy interactions
	  				  are set from parameters in pset
	   Note				: It isn't mandatory that there be
	  				  any shift interactions
	   Note				: Isotopes MUST be set prior to here.
	  				  Isotopes used to avoid confusion in
	  				  IntCSAVec over possible e- spins   */

void solid_sys::setCs(const ParameterSet& pset)
  {
  vector<Isotope> Isos = IsoVec();		// Get vector system isotopes
  Cvec = IntCSAVec(Isos, pset, 0);		// Use class constructor
  }

// ---------------------------------------------------------------------------- 
//	          Functions To Set Quadrupolar Interactions
// ---------------------------------------------------------------------------- 

/* We Must Enforce That Only Spins With I>0.5 Have A Quadrupolar Interactions
   Defined.  We Don't Care If There Isn't A Quadrupolar Interaction Defined
   For A Spin Even If It Does Have I>0.5				     */


	// Input		ssys	: Solid spin system (this)
	// 			pset	: A parameter set
        //                      warn    : Warning output label
        //                                 0 = no warnings
        //                                 1 = warnings
        //                                >1 = fatal warnings
	// Output		none	: System quadrupolar interactions are
	//			  	  set from parameters in pset
	// Note				: Assumes that all quadrupolar
	//				  interactions are specified by
	//				  spin index (although their input
	//				  form can vary)

void solid_sys::setQs(const ParameterSet& pset)
  {
  vector<Isotope> Isos = IsoVec();		// Get vector system isotopes
  Qvec = IntQuadVec(Isos, pset, 0);	// Use class constructor
  }

// ---------------------------------------------------------------------------- 
//	            Functions To Set Electron G Interactions
// ---------------------------------------------------------------------------- 

/* Each Spin That Is An Electron May Have An Electron G Tensor Interaction Set.
   We Don't Particularly Care If There Isn't A G Interaction Defined For Any
   Particular Electron Spin, The Vector Of Such Interactions Will Simply Have
   An Zero Values If No G Anisotropy Has Been Specified.

	   Input		ssys	: Solid spin system (this)
	   			pset	: A parameter set
	   Output		none	: System electron G interactions
	  				  are set from parameters in pset
	   Note				: It isn't mandatory that there be
	  				  any electron G interactions
	   Note				: Isotopes MUST be set prior to here.
	  				  Isotopes used to avoid confusion in
	  				  IntG over possible nuclear spins   */

void solid_sys::setGs(const ParameterSet& pset)
  {
  int i, ns = spins();			// Number of spins in system
  vector<Isotope> Isos = IsoVec();	// Get vector system isotopes
  Gvec = IntGVec(Isos, pset, 0);	// Use IntGVec class constructor
  for(i=0; i<ns; i++)			// Must set the isotropic g factors
    if(Isos[i].electron()) 		// as they may be set by {gxx,gyy,gzz}
      gfactor(i, Gvec[i].iso());
  }

// ----------------------------------------------------------------------------
//                  Functions To Set Hyperfine Interactions
// ----------------------------------------------------------------------------

/* Each Spin Pair That Involves An Electron And A Nucleon May Have A Hyperfine
   Interaction Set. We Don't Particularly Care If There Isn't A HF Interaction
   Defined For Any Particular Spin Pair, The Vector Of Such Interactions Will
   Simply Have An Zero Values If No HF Anisotropy Has Been Specified.

           Input                ssys    : Solid spin system (this)
                                pset    : A parameter set
           Output               none    : System hyperfine interactions
                                          are set from parameters in pset
           Note                         : It isn't mandatory that there be
                                          any hyperfine interactions
           Note                         : Isotopes MUST be set prior to here.
                                          Isotopes used to avoid confusion
                                          over electron vs nuclear spins     */

void solid_sys::setHFs(const ParameterSet& pset)
  {
  int i, j, k, ns = spins();		// Number of spins in system
  HFvec = IntHFVec(pset, IsoVec(), 0);	// Set interactions vector
  for(i=0, k=0; i<ns-1; i++)		// Must set the isotropic A factors
    for(j=i+1; j<ns; j++, k++) 		// as they may be set by {Axx,Ayy,Azz}
      if(nepair(i,j))			// but this is done only if we have 
        A(i,j,HFvec[k].iso());		// an electron-nucleon pair
  }
 
// ----------------------------------------------------------------------------
//                  Functions To Set The Full Spin System
// ----------------------------------------------------------------------------

	// Input		ssys	: Solid spin system (this)
	// 			pset	: A parameter set
        //                      idx     : Parameter index value
        //                      warn    : Warning output label
        //                                 0 = no warnings
        //                                 1 = warnings
        //                                >1 = fatal warnings
	// Output		none	: The entire system is set
	//			  	  from parameters in pset
	// Note				: This uses the assignment from pset
	//				  It exists to support prefix indices
	// Note				: Since solid_sys doesn't exactly
	//				  match base class spin_system's
	//				  assignment from pset, herein we	
	//				  go step by step.

int solid_sys::setSsys(const ParameterSet& pset, int idx, int warn)
  {
//   Use Only Parameters With Prefix [#], Clip Off Prefix From Names

  ParameterSet subpset;				// Copy parameter set
  if(idx != -1) subpset = pset.strip(idx);	// to glean out those for
  else          subpset = pset;			// the specified prefix index
  int ns = getSpins(subpset, 2);		// Get # of spins (spin_sys)
  *this = solid_sys(ns);			// (Re)Set system for ns spins

//             Set Values For Base Spin System (class spin_sys)

  int w=warnings();			// Store current warnings level
  setIs(subpset);			// Set isotope types (spin_sys)
  setName(subpset);			// Read in the system name (spin_sys)

// sosi
//         Set Values For Isotropic Spin System (class spin_system)
//          (Some Of These May Be Overridden By Tensor Quantities)

  setOm(subpset);		 	// Get spect. frequency (spin_system)
  setShifts(subpset); 			// Read in isotropic shifts
  warnings(1);				// Set to warn if J & e- only
  setJs(subpset); 			// Read in scalar couplings
  warnings(w);				// Reset original warnings level

//              Set Values For Spin Coordinates If They Exist
//	       Set Dipolar Interactions Using Spin Isotope Types 
// (Dipolar Interactions Will Use Coordinates Too If They Exist For All Nuclei)

  int cc = setCoords(subpset);		// Read in spin coordinates
  setDs(subpset, cc);			// Read in dipolar interactions

//         Set Values For Rank 2 Spin & Spin-Spin Interactions

  setCs(subpset);			// Read shift anisotropy interactions
  setQs(subpset);			// Read in quadrupolar interactions
  setGs(subpset);			// Read electron G interactions
  setHFs(subpset);			// Read hyperfine coupling interactions
// sosik - need to figure out this return...
int TF = 1;
  return TF;
  } 

// ____________________________________________________________________________
// iii                 SOLID SYSTEM ADMINISTRATION FUNCTIONS
// ____________________________________________________________________________
 
 
void solid_sys::ResetSOps(int spin)

        // Input                ssys    : Solid spin system (this)
        // Output               none    : Resets any spin operators for
        //                                specified spin
 
  {
spin++;
  return;
  }



// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// A              SOLID SPIN SYSTEM CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 



	// Input		ssys  : Solid spin system (this)
	// 			ssys1 : Solid spin system
	// Output		none  : Solid spin system constructed
	//		                equivalent to ssys1

solid_sys::solid_sys(int ns)                 : spin_system(ns) {zero_cindx();} 
solid_sys::solid_sys(const solid_sys &ssys1) : spin_system(ssys1)
  {
  Dvec  = ssys1.Dvec;			// Copy any dipolar interactions
  Cvec  = ssys1.Cvec;			// Copy any SA interactions
  Qvec  = ssys1.Qvec;			// Copy any quadrupolar interactions
  Gvec  = ssys1.Gvec;			// Copy any G  interactions
  HFvec = ssys1.HFvec;			// Copy any HF interactions
  cflags = ssys1.cflags;		// Coordinate flags
  SCoords = ssys1.SCoords;		// Copy spin coordinates
  }


void solid_sys::operator= (const solid_sys& ssys)
  {
  spin_system::operator=(ssys);		// Equate the spin system part
  Dvec  = ssys.Dvec;			// Copy any dipolar interactions
  Cvec  = ssys.Cvec;			// Copy any CSA interactions
  Qvec  = ssys.Qvec;			// Copy any Quadrupolar interacs
  Gvec  = ssys.Gvec;			// Copy any G   interactions
  HFvec = ssys.HFvec;			// Copy any HF  interactions
  cflags  = ssys.cflags;		// Coordinate flags
  SCoords = ssys.SCoords;		// Copy spin coordinates
  }

solid_sys::~solid_sys () { }	

// ____________________________________________________________________________ 
// B             SPIN COORDINATES VECTOR ACCESS & MANIPULATIONS
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//                          Get Spin Coordinates
// ---------------------------------------------------------------------------- 
 
coord solid_sys::getCoord(int i) const
 
        // Input                ssys    : Solid spin system
        //                      i       : spin index 
        // Output               coord   : Spin i coordinate
	// Note				: Non-existent coordinates
	//				  return as point (0,0,0)

  { 
  if(i<0 || i>spins()) 			// Insure valid spin index
    ssys_fatal(105, Gdec(i));
  if(i>SCoords.size()) return coord0;
  return SCoords(i);
  }
 

coord_vec solid_sys::getCoords() const { return SCoords; }

	// Input		ssys	: Solid spin system
	// Output		cvec	: Spin coordinates
	// Note				: Non-existent coordinates
	//				  return an empty coordinate vector


// ---------------------------------------------------------------------------- 
//                          Set Spin Coordinates
// ---------------------------------------------------------------------------- 
 
// sosiz - 1.) still need to deal with when is coord vector 0 and when are the
//             coordinate flags necessary?
//         2.) need to see if the linked list of interactions is working here

void solid_sys::setCoord(int i, coord& pt)
 
        // Input                ssys    : Solid spin system
        //                      i       : spin index 
	//			pt	: spin coordinate
        // Output               none    : Spin i coordinate set to pt
	// Note				: Coordinate can only be set if there
	//				  already exists a complete set!
 
  {
  string I, S;
  if(i<0 || i>spins()) 			// Insure valid spin index
    ssys_fatal(105, Gdec(i));
  if(!SCoords.size()) 			// Insure coordinates exist
    ssys_fatal(107, Gdec(i));
  SCoords.put(pt, i);			// Set the spin coordinate
  int j,k,l,ns=spins();			// Get the number of spins
  for(j=0,l=0; j<ns-1; j++)		// Loop over all spin pairs
    {
    I = symbol(j);
    for(k=j+1; k<ns; k++,l++)		// and see which involve the
      {
      S = symbol(k);
      if(k==i || j==i)			// new spin.  Check if the spin
        {				// is involved in any dipolar
        if(Dvec.DCC(l))			// interaction and, if so,
          Dvec(l) = IntDip(I,S,SCoords(j),SCoords(k));
        }
      }
    }
  }

void solid_sys::setCoords(const coord_vec& cvec)

	// Input		ssys	: Solid spin system
	// 			cvec	: Spin coordinates
	// Output		none    : Spin coordinates set to cvec
	// Note				: The size of cvec must be at least as
	//				  big as the number of system spins
	// Note				: The number of spins and the isotope
	//				  types should be know prior to this

  {
  int i, ns=spins();
  if(cvec.size() < ns)			// First insure there is a 
    ssys_fatal(50);			// coordinate for each spin
  if(cvec.size() > ns)			// If more coordinates, copy
    {					// coordinates 1 by 1 into
    SCoords = coord_vec(ns);		// the system's vector
    for(i=0; i<ns; i++)
      SCoords.put(cvec(i), i);
    } 					// For proper number of
  else SCoords = cvec;			// coordinates, just copy
  Dvec = IntDipVec(IsoVec(),SCoords,0);	// set dipolar interactions using em
  }

// ____________________________________________________________________________ 
// C                DIPOLAR INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________ 

/* The functions in this section provide access to the dipolar interactions in
   the system.

           Input                ssys    : Solid spin system
                                sI      : Spin I of spin system [0, nspins-1]
                                sS      : Spin S of spin system [0, nspins-1]
                                            1: Asymmetry [0, 1]
                                            2: Theta (degrees, down from +z)
                                            3: Phi   (degrees, over from +z)
           Output               none    : Get/Set dipolar interaction value
                                          for spicified spin pair

    Function  Arguments                        Result
    ========  =========  ======================================================
      DCC      i,j,dcc   Sets the dipolar coupling constant to dcc (Hz)
      DCC        i,j     Returns the dipolar coupling constant to dcc (Hz)
     Ddelz     i,j,dcc   Sets the dipolar coupling constant to dcc (Hz)
     Ddelz       i,j     Returns the dipolar coupling constant to dcc (Hz)
      DEta     i,j,dcc   Sets the dipolar interaction asymmetry [0, 1]
      DEta       i,j     Returns the dipolar interaction asymmetry [0,1]
     DTheta    i,j,the   Sets the dipolar interaction theta orientation (deg)
     DTheta      i,j     Returns the dipolar interaction theta orientat. (deg)
      DPhi     i,j,phi   Sets the dipolar interaction phi orientation (deg)
      DPhi       i,j     Returns the dipolar interaction phi orientat. (deg)
     IntDip      i,j     Returns the rank 2 dipolar interaction between i & j
    IntDipVec            Returns a vector of all system dipolar interactions

   Note that DCC = delzz for dipolar interactions and that setting DCC will
   kill any coordinates associated between two spins. The dipolar interaction
   asymmetry is typically zero. The theta orientation is that down from the
   lab frame to the PAS +z axis and is restricted to [0,180]. The phi
   orientation is from the lab +x axis to the PAX +x axis and it is restricted
   to [0, 360].  Note that any changes to theta and phi will kill any spin
   coordinates.
*/


// ----------------------------------------------------------------------------
//                   Generic Dipolar Value Access Functions
// ----------------------------------------------------------------------------

/*                                        <=0: DCC coupling of spin (Hertz)
                                            1: Asymmetry [0, 1]
                                            2: Theta (degrees, down from +z)
                                            3: Phi   (degrees, over from +z)
           Output               none    : Get/Set dipolar interaction value
                                          for spicified spin pair            */

void solid_sys::DValue(int sI, int sS, double val, int type)

  {
  if(!check_spins(sI, sS)) ssys_fatal(1); 	// Check that spins exist
  if(symbol(sI)=="e-" || symbol(sS)=="e-")	// Neither spin can be e-
    ssys_fatal(31);
  if(type != 1)					// Adjusting anything but eta 
    { 						// forces coordinate zeroing
//    if(cflags) delete[] cflags;			// or we risk an inconsistent
    SCoords = coord_vec(spins());		// set of dipolar interactions
    }
  int id = pairidx(sI, sS);			// Get dipole index
  if(!Dvec.DCC(id) && type>0)			// Can't set these values
    ssys_fatal(30);				// unless DCC is non-zero
  if(type<=0 && Dvec.Izval(id)<0.5)		// If setting DCC & no dip.int.
    { 						// construct a brand new one
    Dvec(id) = IntDip(qn(sI), qn(sS), val);
    return;
    }
  switch(type)
    {
    default:
    case 0: Dvec.DCC(id,    val); break;        // Here for DCC
    case 1: Dvec.Deta(id,   val); break;        // Here for asymmetry
    case 2: Dvec.Dtheta(id, val); break;        // Here for theta
    case 3: Dvec.Dphi(id,   val); break;        // Here for phi
    }
  }
     
double solid_sys::DValue(int sI, int sS, int type) const
  {
  if(!check_spins(sI, sS)) ssys_fatal(1);	// Insure spins exist.
  int id = pairidx(sI, sS);			// Get dipole index
  double rval;
  switch(type)
    {
    default:
    case 0: rval = Dvec.DCC(id);    break;	// Here for DCC
    case 1: rval = Dvec.Deta(id);   break;	// Here for asymmetry
    case 2: rval = Dvec.Dtheta(id); break;	// Here for theta
    case 3: rval = Dvec.Dphi(id);   break;	// Here for phi
    }
  return rval;
  }  

//                         Dipolar Coupling Constants
 
void   solid_sys::DCC(int   sI, int sS, double dcc) { DValue(sI, sS, dcc, 0); }
double solid_sys::DCC(int   sI, int sS) const       { return DValue(sI,sS,0); }
void   solid_sys::Ddelz(int sI, int sS, double dcc) { DValue(sI, sS, dcc, 0); }
double solid_sys::Ddelz(int sI, int sS) const       { return DValue(sI,sS,0); }

//                        Dipolar Asymmetry Values
 
void   solid_sys::Deta(int sI, int sS, double DE) { DValue(sI, sS, DE, 1);}
double solid_sys::Deta(int sI, int sS) const      { return DValue(sI, sS, 1); }

//               Dipolar Theta Orientation (Down From +z Axis)
 
void   solid_sys::Dtheta(int sI, int sS, double dthe) {DValue(sI,sS,dthe,2);}
double solid_sys::Dtheta(int sI, int sS) const        {return DValue(sI,sS,2);}
 
//               Dipolar Phi Orientation (Over From +x Axis)
 
void   solid_sys::Dphi(int sI, int sS, double dphi) { DValue(sI, sS, dphi, 2);}
double solid_sys::Dphi(int sI, int sS) const        { return DValue(sI,sS,2); }

//                     Dipolar Spin Tensor Components

matrix solid_sys::DTcomp(int sI, int sS, int m) const
  {
  if(!check_spins(sI, sS)) ssys_fatal(1); 	// First check pair exists
  int di = pairidx(sI, sS);			// Get dipole index
// sosi
//  return DipVec[di].T2m(HSvect(), sI, sS, m);
return matrix(di,di);
  }
 
//                         Full Dipolar Interaction

const IntDip& solid_sys::getDipInt(int sI, int sS) const
  {
  if(!check_spins(sI, sS)) ssys_fatal(1); 	// First check pair exists
  return Dvec.getcref(pairidx(sI,sS));		// Return requested IntDip
  }

const IntDip&    solid_sys::getDipInt(int dip) const { return Dvec.getcref(dip); }	
const IntDipVec& solid_sys::getDipVec()        const { return Dvec; }

// ____________________________________________________________________________ 
// C            SHIFT ANISOTROPY INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________ 

// sosi - things to do here....
//	1. This must be associated with all changes to the spin quantum number.
//	   That implies that it should be affected by isotope, element, etc. 
//	2. These should check that a.) I>1/2 and b.) Interaction exists
 
	
/******************************************************************************
*******************************************************************************
**									     **
**  Input     ssys    -	Solid spin system                                    **
**	      spin    -	Spin of spin system [0, nspins)                      **
**	      csa     - Chemical shift anisotropy (PPM)                      **
**	      Cdelz   - Chemical shift anisotropy delzz in Hz                **
**            Ceta    - Chemical shift anisotropy asymmetry (unitless)       **
**            Cthe    - Chemical shift anisotropy orientation (degrees)      **
**            Cphi    - Chemical shift anisotropy orientation (degrees)      **
**									     **
**  Output    void    - Sets value specified                                 **
**            double  - Returns value requested (Hz, unitless, or degrees)   **
**            IntCSA  - Returns CSA interaction				     **
**									     **
**  Function  QCC     - Get/Set Quadrupolar coupling constants	             **
**            Qdelz   - Get/Set Quadrupolar spatial tensor delzz values      **
**            Qeta    - Get/Set Quadrupolar asymmetry values                 **
**            Qtheta  - Get/Set Quadrupolar orientation angle [0,180]        **
**            Qphi    - Get/Set Quadrupolar orientation angle [0,360]        **
**            Qint    - Get Quadrupolar interaction                          **
**									     **
*******************************************************************************
******************************************************************************/

// ----------------------------------------------------------------------------
//                   Generic CSA Value Access Functions
// ----------------------------------------------------------------------------

	// Input		ssys	: Solid spin system
	// 			sI	: Spin I of spin system [0, nspins-1]
        //                      val     : CSA interaction value
        //                                <=0: CSA of spin (PPM)
        //                                  1: Asymmetry [0, 1]
        //                                  2: Theta (degrees, down from +z)
        //                                  3: Phi   (degrees, over from +z)
        // Output               none    : Get/Set CSA interaction value
        //                                for specified spin

void solid_sys::CValue(int sI, double val, int type)

  {
  if(!check_spin(sI))				// First check that spin
    ssys_fatal(1);				// exists.
  switch(type)
    {
    default:
    case 0: Cvec.CSA(sI, val);   break;		// Here for CSA
    case 1: Cvec.eta(sI, val);   break;		// Here for asymmetry
    case 2: Cvec.theta(sI, val); break;		// Here for theta
    case 3: Cvec.phi(sI, val);   break;		// Here for phi
//    case 4: Cvec.delz(sI, val);  break;		// Here for delzz
    }
  }
     
 
double solid_sys::CValue(int sI, int type) const
 
  {
  if(!check_spin(sI))				// First check that spin
    ssys_fatal(1);				// exists.
  double rval=0;;
  switch(type)
    {
    default:
    case 0: rval = Cvec.CSA(sI);   break;	// Here for CSA
    case 1: rval = Cvec.eta(sI);   break;	// Here for asymmetry
    case 2: rval = Cvec.theta(sI); break;	// Here for theta
    case 3: rval = Cvec.phi(sI);   break;	// Here for phi
//    case 4: rval = Cvec.delz(sI);  break;	// Here for delzz
    }
  return rval;
  }  


//                            CSA & delzz Values

	// Input		ssys	: Solid spin system
	// 			sI	: Spin I of spin system [0, nspins-1]
	// 			CSA	: Shift anisotropy of spin (CSA)
	// Output		none	: Get/Set CSA of spin I
	// Note				: Defined in class IntCSA as equal to
	//			          1.5 times tensor delzz value
	// Input		ssys	: Solid spin system
	// 			sI	: Spin I of spin system [0, nspins-1]
	// 			CE	: CSA asymmetry of spin
	// Output		none	: Get/Set CSA asymmetry
	// Note				: Defined class IntRank2 as [0,1]
 
	// Input		ssys	: Solid spin system
	// 			sI	: Spin I of spin system [0, nspins-1]
	// 			ctheta	: Get/Set spin angle (deg)
	// Output		none	: Set angle
	// Note				: Defined class IntRank2A as [0,180]
	// Input		ssys	: Solid spin system
	// 			sI	: Spin I of spin system [0, nspins-1]
	// 			dphi	: CSA angle of spin (deg)
	// Output		none	: Get/Set angle
	// Note				: Defined in class IntRank2A as [0,360]
	// Input		ssys	: Solid spin system
	// 			sI	: Spin I of spin system [0, nspins-1]
	// Output		CI	: Return rank 2 CSA interaction


 
void   solid_sys::CSA(int   sI, double cs) { CValue(sI, cs, 0); }
double solid_sys::CSA(int   sI) const      { return CValue(sI, 0); }
void   solid_sys::Cdelz(int sI, double dz) { CValue(sI, dz, 4); }
double solid_sys::Cdelz(int sI) const      { return CValue(sI, 4); }

//                            CSA Asymmetry Values
 
void   solid_sys::Ceta(int sI, double CE) { CValue(sI, CE, 1); }
double solid_sys::Ceta(int sI) const      { return CValue(sI, 1); }

//                   CSA Theta Orientation (Down From +z Axis)

void solid_sys::Ctheta(int   sI, double ctheta) { CValue(sI,ctheta,2); }
double solid_sys::Ctheta(int sI) const          { return CValue(sI,2); }
 
//                    CSA Phi Orientation (Over From +x Axis)

void   solid_sys::Cphi(int sI, double cphi) { CValue(sI, cphi, 2); }
double solid_sys::Cphi(int sI) const        { return CValue(sI,2); }
 
//                             Full CSA Interaction

IntCSA solid_sys::getCSAInt(int sI) const
  {
  if(!check_spin(sI)) ssys_fatal(1);		// Check that spin exists.
  return Cvec.get(sI);				// Return requested IntCSA
  }

// ____________________________________________________________________________
// E             QUADRUPOLAR INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________

/* The functions in this section provide access to the quadrupolar interactions
   in the system.

           Input                ssys    : Solid spin system
                                sI      : Spin I of spin system [0, nspins-1]
                                val     : Quadrupolar interaction value
                                          <=0: QCC of spin
                                            1: Asymmetry [0, 1]
                                            2: Theta (degrees, down from +z)
                                            3: Phi   (degrees, over from +z)
           Output               none    : Get/Set quadrupolar interaction
                                          value for specified spin
  									    
                     1/2                                                       
                  [2]                      1                                   
            T   = |-| * I         T    = - - I           T    = 0              
             2,0  [3]    z         2,1     2  +           2,2                  
  									       
                                   2                      3*NQCC	       
             DELZZ = NQCC = QCC = e qQ             w   = --------              
                                                    Q    2I(2I-1)	       
  									       	
                      1/2                  1/2    del               1/2        
           Q     [6*pi]     NQCC      [6*pi]         zz         [pi]	       
         xi    = |----| * --------- = |----|  * --------- = W * |--|           
                 [ 5  ]   2I*(2I-1)   [ 5  ]    2I*(2I-1)    Q  [15]           
  									       

    Function  Arguments                        Result
    ========  =========  ======================================================
      QCC       i,qcc    Sets the quadrupolar coupling constant to qcc (Hz)
      QCC         i      Returns the quadrupolar coupling constant to qcc (Hz)
     Qdelz      i,qcc    Sets the quadrupolar coupling constant to qcc (Hz)
     Qdelz        i      Returns the quadrupolar coupling constant to qcc (Hz)
      QEta      i,qcc    Sets the quadrupolar interaction asymmetry [0, 1]
      QEta        i      Returns the quadrupolar interaction asymmetry [0,1]
     QTheta     i,the    Sets the quad. interaction theta orientation (deg)
     QTheta       i      Returns the quad. interaction theta orientat. (deg)
      QPhi      i,phi    Sets the quadrupolar interaction phi orientation (deg)
      QPhi        i      Returns the quad. interaction phi orientat. (deg)
     IntQuad      i      Returns the rank 2 quadrupolar interaction for spin i
    IntQuadVec           Returns a vector of all system quad. interactions

   Note that the quadrupolar coupling (QCC) is defined to be identical to the
   interaction delzz value.  The theta orientation is that down from the
   lab frame to the PAS +z axis and is restricted to [0,180]. The phi
   orientation is from the lab +x axis to the PAX +x axis and it is
   restricted to [0, 360].                                                   */

// ----------------------------------------------------------------------------
//                   Generic Quadrupolar Value Access Functions
// ----------------------------------------------------------------------------

/*                                        <=0: DCC coupling of spin (Hertz)
                                            1: Asymmetry [0, 1]
                                            2: Theta (degrees, down from +z)
                                            3: Phi   (degrees, over from +z)
           Output               none    : Get/Set dipolar interaction value   */
  
//                   Generic Quadrupolar Value Access Functions

void solid_sys::QValue(int sI, double val, int type)
  {
  if(!check_spin(sI))				// First check that spin
    ssys_fatal(1);				// exists.
  switch(type)
    {
    default:
//    case 0: Qvec.delz(sI,  val); break;		// Here for delzz
    case 1: Qvec.eta(sI,   val); break;		// Here for asymmetry
    case 2: Qvec.theta(sI, val); break;		// Here for theta
    case 3: Qvec.phi(sI,   val); break;		// Here for phi
//    case 4: Qvec.delz(sI,  val); break;		// Here for delzz
    }
  }
     
 
double solid_sys::QValue(int sI, int type) const
 
  {
  if(!check_spin(sI))				// First check that spin
    ssys_fatal(1);				// exists.
  double rval=0;;
  switch(type)
    {
    default:
//   case 0: rval = Qvec.delz(sI);  break;	// Here for delzz
    case 1: rval = Qvec.eta(sI);   break;	// Here for asymmetry
    case 2: rval = Qvec.theta(sI); break;	// Here for theta
    case 3: rval = Qvec.phi(sI);   break;	// Here for phi
    }
  return rval;
  }  

//                          Quadrupolar Coupling Constants 

void   solid_sys::QCC(int   spin, double qcc)   { QValue(spin, qcc,   0); }
double solid_sys::QCC(int   spin) const         { return QValue(spin, 0); }
void   solid_sys::Qdelz(int spin, double delzz) { QValue(spin, delzz, 0); }
double solid_sys::Qdelz(int spin) const	        { return QValue(spin, 0); } 

//                      Quadrupolar Asymmetry Values
 
void   solid_sys::Qeta(int spin, double QE) { QValue(spin, QE,    1); }
double solid_sys::Qeta(int spin) const      { return QValue(spin, 1); }

//             Quadrupolar Theta Orientation (Down From +z Axis)
 
void   solid_sys::Qtheta(int spin, double qthe) {QValue(spin,qthe,  2);}
double solid_sys::Qtheta(int spin) const        {return QValue(spin,2);}
 
//             Quadrupolar Phi Orientation (Over From +x Axis)
 
void   solid_sys::Qphi(int spin, double qphi) {QValue(spin,qphi,  3);}
double solid_sys::Qphi(int spin) const        {return QValue(spin,3);}
 
//                 Quadrupolar Spin Tensor Components

matrix solid_sys::QTcomp(int spin, int m) const
  {
  if(!check_spin(spin)) ssys_fatal(1); 	// First check exists
// sosi
//  return QuadVec[di].T2m(HSvect(), spin, m);
return matrix(5,5);
  }
 
//                        Full Quadrupolar Interaction
         
const IntQuad& solid_sys::getQuadInt(int spin) const
  {
  if(!check_spin(spin)) ssys_fatal(1); 	// First check spin exists
  return Qvec.getcref(spin);			// Return requested IntQuad
  }
const IntQuadVec& solid_sys::getQuadVec()         const { return Qvec; }

// ____________________________________________________________________________
// F              ELECTRON G INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________

/* The functions in this section provide access to the electron G interactions
   in the system.

           Input                ssys    : Solid spin system
                                sI      : Spin I of spin system [0, nspins-1]
                                val     : Electron G interaction value
                                          <=0: DCC coupling of spin (Hertz)
                                            1: Asymmetry [0, 1]
                                            2: Theta (degrees, down from +z)
                                            3: Phi   (degrees, over from +z)
           Output               none    : Get/Set dipolar interaction value
                                          for spicified spin pair

    Function  Arguments                        Result
    ========  =========  ======================================================
      GCC      i,dcc     Sets the electron G coupling constant to dcc (Hz)
      GCC        i       Returns the electron G coupling constant to dcc (Hz)
     Gdelz     i,dcc     Sets the electron G coupling constant to dcc (Hz)
     Gdelz       i       Returns the electron G coupling constant to dcc (Hz)
      GEta     i,dcc     Sets the electron G interaction asymmetry [0, 1]
      GEta       i       Returns the electron G interaction asymmetry [0,1]
     GTheta    i,the     Sets the electron G interaction theta orientation (deg)
     GTheta      i       Returns the electron G interaction theta orientat. (deg)
      GPhi     i,phi     Sets the electron G interaction phi orientation (deg)
      GPhi       i       Returns the electron G interaction phi orientat. (deg)
      IntG       i       Returns the rank 2 electron G interaction for spin i
     IntGVec             Returns a vector of all system electron G interactions
*/

void solid_sys::GValue(int sI, double val, int type)

  {
  if(!check_spin(sI)) ssys_fatal(1);	 	// Check that spin exists
//  if(!Gvec[sI].delz() && type>0)		// Can't set these values
//    ssys_fatal(30);				// unless delzz is non-zero
  if(type<=0 && Gvec[sI].Izval()<0.5)		// If setting delzz & no G int.
    { 						// construct a brand new one
    Gvec(sI) = IntG(qn(sI), val);
    return;
    }
  switch(type)
    {
    default:
//    case 0: Gvec[sI].delz(val);  break;        // Here for delz
    case 1: Gvec[sI].eta(val);   break;        // Here for asymmetry
    case 2: Gvec[sI].theta(val); break;        // Here for theta
    case 3: Gvec[sI].phi(val);   break;        // Here for phi
    }
  }
     
double solid_sys::GValue(int sI, int type) const
  {
  if(!check_spin(sI)) ssys_fatal(1);		// Check that spin exists
  double rval;
  switch(type)
    {
    default:
//    case 0: rval = Gvec[sI].delz();  break;	// Here for DCC
    case 1: rval = Gvec[sI].eta();   break;	// Here for asymmetry
    case 2: rval = Gvec[sI].theta(); break;	// Here for theta
    case 3: rval = Gvec[sI].phi();   break;	// Here for phi
    }
  return rval;
  }  

//                         Dipolar Coupling Constants
 
void   solid_sys::Gdelz(int   sI, double gdzz) { GValue(sI, gdzz, 0); }

//                           Electron G Asymmetry Values

void   solid_sys::Geta(int sI, double ge) { GValue(sI, ge, 1); }
double solid_sys::Geta(int sI) const      { return GValue(sI, 1); }

//               Electron G Theta Orientation (Down From +z Axis)
 
void   solid_sys::Gtheta(int sI, double gthe) {GValue(sI,gthe,2);}
double solid_sys::Gtheta(int sI) const        {return GValue(sI,2);}
 
//               Electron G Phi Orientation (Over From +x Axis)
 
void   solid_sys::Gphi(int sI, double gphi) { GValue(sI, gphi, 2);}
double solid_sys::Gphi(int sI) const        { return GValue(sI,2); }

//                     Electron G Spin Tensor Components

matrix solid_sys::GTcomp(int sI, int m) const
  {
  if(!check_spin(sI)) ssys_fatal(1); 	// First check spin exists
  return (Gvec[sI]).T2m(m, HSvect(), sI);
  }
 
//                         Full Electron G Interaction

IntG solid_sys::getGInt(int sI) const
  {
  if(!check_spin(sI)) ssys_fatal(1); 	// First check spin exists
  return Gvec.get(sI);				// Return requested IntG
  }

IntGVec solid_sys::getGVec() const { return Gvec; }

// ____________________________________________________________________________
// G                HYPERFINE INTERACTION ACCESS & MANIPULATIONS
// ____________________________________________________________________________

/* The functions in this section provide access to the hyperfine interactions 
   in the system.

           Input                ssys    : Solid spin system
                                sI      : Spin I of spin system [0, nspins-1]
                                sS      : Spin S of spin system [0, nspins-1]
                                            1: Asymmetry [0, 1]
                                            2: Theta (degrees, down from +z)
                                            3: Phi   (degrees, over from +z)
           Output               none    : Get/Set hyperfine interaction value
                                          for specified spin pair

    Function  Arguments                        Result
    ========  =========  ======================================================
      DCC      i,j,dcc   Sets the hyperfine coupling constant to dcc (Hz)
      DCC        i,j     Returns the hyperfine coupling constant to dcc (Hz)
     HFdelz    i,j,dcc   Sets the hyperfine coupling constant to dcc (Hz)
     HFdelz      i,j     Returns the hyperfine coupling constant to dcc (Hz)
      HFEta    i,j,dcc   Sets the hyperfine interaction asymmetry [0, 1]
      HFEta      i,j     Returns the hyperfine interaction asymmetry [0,1]
     HFTheta   i,j,the   Sets the hyperfine interaction theta orientation (deg)
     HFTheta     i,j     Returns the hyperfine interaction theta orientat. (deg)
      HFPhi    i,j,phi   Sets the hyperfine interaction phi orientation (deg)
      HFPhi      i,j     Returns the hyperfine interaction phi orientat. (deg)
     IntHF       i,j     Returns the rank 2 hyperfine interaction between i & j
    IntHFVec             Returns a vector of all system hyperfine interactions

   Note that DCC = delzz for dipolar interactions and that setting DCC will
   kill any coordinates associated between two spins. The dipolar interaction
   asymmetry is typically zero. The theta orientation is that down from the
   lab frame to the PAS +z axis and is restricted to [0,180]. The phi
   orientation is from the lab +x axis to the PAX +x axis and it is restricted
   to [0, 360].  Note that any changes to theta and phi will kill any spin
   coordinates.
*/


//                   Generic Hyperfine Value Access Functions

void solid_sys::HFValue(int sI, int sS, double val, int type)

  {
  if(!check_spins(sI, sS)) ssys_fatal(1); 	// Check that spins exist
  int id = sI*spins()+sS;			// Get dipole index
  if(!Dvec.DCC(id) && type>0)			// Can't set these values
    ssys_fatal(30);				// unless DCC is non-zero
  if(type<=0 && Dvec.Izval(id)<0.5)		// If setting DCC & no dip.int.
    { 						// construct a brand new one
    Dvec(id) = IntDip(qn(sI), qn(sS), val);
    return;
    }
  switch(type)
    {
    default:
//    case 0: HFvec[id].delz(val); break; 	// Here for DCC
    case 1: HFvec[id].eta(val);   break;        // Here for asymmetry
    case 2: HFvec[id].theta(val); break;        // Here for theta
    case 3: HFvec[id].phi(val);   break;        // Here for phi
    }
  }
     
double solid_sys::HFValue(int sI, int sS, int type) const
  {
  if(!check_spins(sI, sS))			// First check that two spins
    ssys_fatal(1);				// exist.
  int id = sI*spins()+sS;			// Set dipole index
  double rval;
  switch(type)
    {
    default:
//    case 0: rval = HFvec[id].delz();  break;	// Here for DCC
    case 1: rval = HFvec[id].eta();   break;	// Here for asymmetry
    case 2: rval = HFvec[id].theta(); break;	// Here for theta
    case 3: rval = HFvec[id].phi();   break;	// Here for phi
    }
  return rval;
  }  

//                         Hyperfine Coupling Constants
 
//void   solid_sys::HFCC(int   sI, int sS, double dcc) { HFValue(sI, sS, dcc, 0); }
//double solid_sys::DCC(int   sI, int sS) const       { return HFValue(sI,sS,0); }
void   solid_sys::HFdelz(int sI, int sS, double dcc) { HFValue(sI, sS, dcc, 0); }
double solid_sys::HFdelz(int sI, int sS) const       { return HFValue(sI,sS,0); }

//                        Hyperfine Asymmetry Values
 
void   solid_sys::HFeta(int sI, int sS, double hfe) { HFValue(sI, sS, hfe, 1);}
double solid_sys::HFeta(int sI, int sS) const      { return HFValue(sI, sS, 1); }

//               Hyperfine Theta Orientation (Down From +z Axis)
 
void   solid_sys::HFtheta(int sI, int sS, double hthe) {HFValue(sI,sS,hthe,2);}
double solid_sys::HFtheta(int sI, int sS) const        {return HFValue(sI,sS,2);}
 
//               Hyperfine Phi Orientation (Over From +x Axis)
 
void   solid_sys::HFphi(int sI, int sS, double hphi) { HFValue(sI, sS, hphi, 2);}
double solid_sys::HFphi(int sI, int sS) const        { return HFValue(sI,sS,2); }

//                     Hyperfine Spin Tensor Components

matrix solid_sys::HFTcomp(int sI, int sS, int m) const
  {
  if(!check_spins(sI, sS)) ssys_fatal(1); 	// First check pair exists
  int hi = pairidx(sI, sS);			// Get spin pair index
  return HFvec[hi].T2m(m, HSvect(), sI, sS);
  }
 
//                         Full Hyperfine Interaction

IntHF solid_sys::getHFInt(int sI, int sS) const
  {
  if(!check_spins(sI, sS)) ssys_fatal(1); 	// First check pair exists
  return HFvec.get(pairidx(sI,sS));		// Return requested IntDip
  }

// ____________________________________________________________________________ 
// H                    ROTATION AND ORIENTATION FUNCTIONS
// ____________________________________________________________________________ 

/* These functions are used to rotate/orient either the spin system or a 
   specific interaction.  Note that orienting the spin system demands that
   ALL interactions are reoriented.                                          */

// ---------------------------------------------------------------------------- 
//	                 Functions To Orient The Spin System
// ---------------------------------------------------------------------------- 

/* The spin system as a whole is typically reoriented for two reasons. First
   for a powder average where system orientation with respect to the applied
   magnetic field will span all possibilities.  Second when the sample is 
   undergoing MAS the system is first oriented at the magic angle then rotated
   about that angle at a specified frequency.  There are of course other 
   possibilities, such as oriented systems                                   */

/*
void solid_sys::orient(double theta, double phi)
  {
  SCoords.rotate(ostr);			// Rotate spin coordinates
  Dvec.rotate(ostr); 			// Rotate dipolar interactions 
  Cvec.rotate(ostr);			// Rotate CSA interactions
  Qvec.rotate(ostr); 			// Rotate quad interactions
  Gvec.rotate(ostr); 			// Rotate e- G interactions
  HFvec.rotate(ostr); 			// Rotate hyperfine interacts.
  return ostr;
  }
*/

// ____________________________________________________________________________ 
// I                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________ 

// ---------------------------------------------------------------------------- 
//	  Functions To Make A Parameter Set From A Solid Spin System
// ---------------------------------------------------------------------------- 

/* Note that the preferred means of specifying system parameters are used for
   filling up the parameter set that will be returned regardless of whether
   the system was originally generated from parameter variants.  The perferred
   parameters for the base classes (class spin_system & class spin_sys) are
   defined in the associated class code. Additional default parameters for this
   class can be found in the rank 2 interaction vector classes.

	   Input		ssys	: Solid spin system
	    			pset	: Parameter set
                                          prefix [#] in output names
           Output               void	: System parameters are
                                          are added ot the parameter set
                                          with interaction index idx
           Note                         : The parameter names & types
                                          output here WILL match those used
                                          in setting the spin system 
                                          from parameters sets               */

solid_sys::operator ParameterSet( ) const
  { ParameterSet pset; pset += *this; return pset; }

void operator+= (ParameterSet& pset, const solid_sys &ssys)
  { ssys.PSetAdd(pset); }

void solid_sys::PSetAdd(ParameterSet& pset, int idx) const
  {
  spin_system::PSetAdd(pset, idx);	// Add isotropic system parameters
  int nc = cflags.size();		// Number of coordinates we know
  if(nc) SCoords.PSetAdd(pset, idx);// Next, add in coordinates existing
  else   Dvec.PSetAdd(pset, idx);	// else add in dipolar interactions
// sosix

  string pname;				// Now add unique solid_sys parameters
  string pval;				// by setting single parameter name, 
  string pstate;			// value and statement as strings
  SinglePar par;			// then adding the parameter to pset
//  int i, ns = ssys.spins();		// This is the number of spins




/*
  if(ssys.Cs)				// Next output any exisiting shift
    {					// anisotropy interactions (per spin)
    pstate = "Shift Anisotropy Interaction (PPM)";
    for(i=0; i<ns; i++)
      {
      if(ssys.Cdelz(i))
        {
        pname = string("AC(");
        pname += Gdec(i);
        pname += string(")");
        pval  = string("(")  + Gform("%f7.3", ssys.Cdelz(i));
        if(ssys.Ceta(i)) pval += string(", ") + Gform("%f5.3", ssys.Ceta(i));
        else             pval += string(", 0.0");
        if(ssys.Ctheta(i)) pval += string(", ") + Gform("%f7.3", ssys.Ctheta(i));
        else               pval += string(", 0.0");
        if(ssys.Cphi(i)) pval += string(", ") + Gform("%f7.3", ssys.Cphi(i));
        else             pval += string(", 0.0");
        pval += string(")");
        par = SinglePar(pname, pval, pstate);
        pset.add(par);
        }
      }
    }
  if(ssys.Qs)				// Now output any existiong quadrupolar
    {					// interactions (per spin)
    string pstate0 = "Quadrupolar Interaction ";
    double delzz;
    for(i=0; i<ns; i++)
      {
      delzz = ssys.Qdelz(i);
      if(delzz > 1.e6)
        {
        pstate = pstate0 + string("(MHz)");
        pname = string("AQMHz(");
        pname += Gdec(i);
        pname += string(")");
        delzz /= 1.e6;
        }
      else if(delzz > 1.e3)
        {
        pstate = pstate0 + string("(KHz)");
        pname = string("AQ(");
        pname += Gdec(i);
        pname += string(")");
        delzz /= 1.e3;
        }
      else
        {
        pstate = pstate0 + string("(Hz)");
        pname = string("AQHz(");
        pname += Gdec(i);
        pname += string(")");
        }
      if(delzz)
        {
        pval  = string("(")  + Gform("%f7.3", delzz);
        if(ssys.Qeta(i)) pval += string(", ") + Gform("%f5.3", ssys.Qeta(i));
        else             pval += string(", 0.0");
        if(ssys.Qtheta(i)) pval += string(", ") + Gform("%f7.3", ssys.Qtheta(i));
        else               pval += string(", 0.0");
        if(ssys.Qphi(i)) pval += string(", ") + Gform("%f7.3", ssys.Qphi(i));
        else             pval += string(", 0.0");
        pval += string(")");
        par = SinglePar(pname, pval, pstate);
        pset.add(par);
        }
      }
    }
*/



  return;
  } 


// ---------------------------------------------------------------------------- 
//	  Functions To Make A Solid Spin System From A Parameter Set
// ---------------------------------------------------------------------------- 


void solid_sys::operator= (const ParameterSet& pset) { setSsys(pset); }

	// Input		ssys	: Solid spin system (this)
	// 			pset	: A parameter set
	// Output		none	: Solid spin system filled with
	//			          parameters in pset
	// Note				: We just use a member function
	//				  which is needed to allow for
	//				  [#] prefixed parameter names




// ----------------------------------------------------------------------------
//    Functions To Output Solid State System To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

int solid_sys::write(const string &fileout, int idx, int warn) const

	// Input		ssys	: Solid spin system (this)
        //                      fileout	: Output file name
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Solid state system is written
        //                                as a parameter set to file fileout

  {
  if(!spins()) return 1;		// Nothing if no spins
  ofstream ofstr(fileout.c_str());	// Open filename for input
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!write(ofstr, idx, w2))            // If file bad then exit
    {
    if(warn==1) ssys_error(1,fileout);	// Problems with file
    else if(warn)
      {
      ssys_error(1,fileout,1);		// Problems with file
      ssys_fatal(20);		// Fatal error
      }
    return 0;
    }
  ofstr.close();                        // Close it now
  return 1;
  }


int solid_sys::write(ofstream& ofstr, int idx, int warn) const

	// Input		ssys	: Solid spin system (this)
        //                      ofstr   : Output file stream
        //                      idx     : Parameter index value used for
        //                                prefix [#] in output names
        //                      warn    : Warning level
        // Output               none    : Solid state system is written as a
        //                                parameter set to output filestream
	// Note				: This depends on function PSetAdd!

  {                                                         
  if(!spins()) return 1;		// Nothing if no spins
  if(!ofstr.good())                     // If file bad then exit
    {
    if(warn > 1)
      {
      ssys_error(22, 1);		//      Problems with file
      ssys_fatal(23);		//      It's a fatal error
      }
    if(warn) ssys_error(22);		//      Problems with file
    return 0;
    }
  ParameterSet pset;                 // Declare a parameter set
  PSetAdd(pset, idx);			// Add system to parameter set
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!pset.write(ofstr, w2))            // Use parameter set to write
    {
    if(warn)
      {
      ssys_error(22, 1);		// Problems writing to filestream
      if(warn>1) ssys_fatal(23);	// Fatal error
      }
    return 0;
    }  
  return 1;
  }  

 
// ____________________________________________________________________________
// F               SOLID STATE SPIN SYSTEM INPUT FUNCTIONS
// ____________________________________________________________________________
 
// ----------------------------------------------------------------------------
//   Direct Read of Solid State System From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

/* The following read functions utilize a single spin system index, for system
   being read.  They'll will try to read the spin system parameters that are
   prefixed by [#] where #=idx.  The default value, idx=-1, flags that there 
   is no prefix on the spin system paramters.

	   Input		ssys	: Solid spin system (this)
	   			filename: Input filename
                           or   pset    : Parameter set
                                idx     : Parameter index value used for
                                          prefix [#] in input names
                                warn    : Warning output label
                                           0 = no warnings
                                           1 = warnings
                                          >1 = fatal warnings
	   Output		TF	: Solid spin system filled with
	  				  parameters read from file or pset
	  				  Returns true if read properly      */

int solid_sys::read(const string &filename, int idx, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename, warn?1:0))	// Read in pset from file
    {					// If we cannot read the file
    if(warn)				// then issue warnings as desired
      {
      ssys_error(1, filename, 1);	//	Problems with file
      if(warn > 1) ssys_fatal(21);	//	Its a fatal error
      else         ssys_error(21);	//	or a warning issued
      }
    return 0;
    }
  return read(pset, idx, warn);		// Fill ssys with parameters
  }
                                                                               
int solid_sys::read(const ParameterSet& pset, int idx, int warn)
  { 
  int TF = setSsys(pset, idx, warn?1:0);	// Use overload to read 
  if(!TF)					// If setSsys didn't handle
    {						// the system read from pset
    if(warn)					// then we'll issue some
      {						// warnings if desired
      if(warn > 1) ssys_fatal(102, Gdec(idx));// Fatal error
      else         ssys_error(102, Gdec(idx));	// or a warning issued
      }
    return 0;
    }
  return TF;
  }
 
// ----------------------------------------------------------------------------
//          Interactive Read of Solid State System From An ASCII File
// ----------------------------------------------------------------------------

/* This function will ask the user for the ASCII (GAMMA parameter set) file
   containing the spin system parameters. There are simply too many to ask for
   individually.

           Input                ssys    : Solid spin system (this)
                                argc    : Number of arguments
                                argv    : Vector of argc arguments
                                argn    : Argument index
           Output               void    : The parameter argn of array argc
                                          is used to supply a filename
                                          from which the spin system is read
                                          If the argument argn is not in argv,
                                          the user is asked to supply a filename
           Note                         : File must be ASCII containing
                                          recognized sys parameters
           Note                         : Spin system is modifed (filled)    */

string solid_sys::ask_read(int argc, char* argv[], int argn)
  {
  string filename;				// Name of spin system file  
  query_parameter(argc, argv, argn,		// Get filename from command
       "\n\tSpin system filename? ", filename); // line or ask for them
  read(filename);		           	// Read system from filename
  return filename;
  }

string solid_sys::ask_read(int argc, char* argv[], int argn, const string& def)
  {
  string msg = "\n\tSpin system filename ["     // Query we will ask if
             + def + "]? ";                     // it is needed
  string filename = def;                        // Name of spin system file
  ask_set(argc,argv,argn,msg,filename);         // Or ask for it
  read(filename);                               // Read system from filename
  return filename;                              // Return filename
  }


// ____________________________________________________________________________ 
//                          STANDARD I/O FUNCTIONS
// ____________________________________________________________________________ 

/* These functions output either part or the entire solid state spin system.
   Many of these functions will print out lists of the interactions that are
   defined in the spin system as stored in the interaction vectors.

   Function                                 Output
   --------  ------------------------------------------------------------------
   printPs   Nuclear coordinates (for dipolar interactions).  This function is
             used instead of class vector output because some spins may have no
             coordinates. Unit flag: 1=Angstroms 2=nmeters 3=meters x=none
   printDs   Dipolar interactions sent to output stream
   printCs   Shift anisotropy interactions sent to output stream
   printQs   Quadrupolar interactions sent to output stream
   print     Entire spin system sent to output stream
     <<      Standard output of entire spin system                           */
 
// ----------------------------------------------------------------------------
//                    Print Nuclear Coordinates If Present
// ----------------------------------------------------------------------------

ostream& solid_sys::printPs(ostream& ostr, int units) const
  {
// sosi - need to set units!
units=1;

  if(!SCoords.size()) return ostr; 		// Exit if no coordinates
  int i, j, npts=spins();
  int cexists =0;
  for(i=0; i<npts && !cexists; i++)
    cexists += cflags[i];
  if(!cexists) return ostr;			// Exit if no coordinates
  double fact=1.0;
  ostr << "\nCoordinates";
  switch(units)
    {
    case 1:  fact = 1.e10; ostr << "\n(Angstroms)"; break;
    case 2:  fact = 1.e9;  ostr << "\n (nmeters)";  break;
    case 3:  fact = 1.0;   ostr << "\n (meters)";   break;
    default: fact = 1.0;                            break;
    }
  ostr << "\nSpin     :";
  for(i=0; i<spins(); i++) ostr << Gdec("%10d",i);
  for(i=0; i<3; i++)			// Loop x, y, z values
    {
    switch(i)
      {
      case 0: ostr << "\nX        :";
              for(j=0; j<npts; j++)
                {
                if(!cflags[j]) ostr << "  --------";
                else           ostr << Gform("%10.2f", fact*SCoords.x(j));
                }
              break;
      case 1: ostr << "\nY        :";
              for(j=0; j<npts; j++)
                {
                if(!cflags[j]) ostr << "  --------";
                else           ostr << Gform("%10.2f", fact*SCoords.y(j));
                }
              break;
      case 2: ostr << "\nZ        :";
              for(j=0; j<npts; j++)
                {
                if(!cflags[j]) ostr << "  --------";
                else          ostr << Gform("%10.2f", fact*SCoords.z(j));
                }
              break;
      }
    }
  ostr << "\n";
  return ostr;
  }
 
// ----------------------------------------------------------------------------
//                        Print Dipolar Interactions
// ----------------------------------------------------------------------------

ostream& solid_sys::printDs(ostream& ostr) const
  {
//		Just Exit If No Dipolar Interactions Present

  if(!Dvec.size())    return ostr;	// Bail if no interactions at all
  if(!Dvec.nonzero()) return ostr;	// Bail if no non-zero interactions

//	     Test That There Are Non-Zero Dipolar Interactions

  ostr << "\nDipolar Interactions";		// Start with the header
  int i,j;					// Temp spin indices
  int k,l;					// Temp interaction indices
  ostr << "\nSpin     :"; 			// Write row of spin labels
  for(i=0; i<spins(); i++)			// for column indices
    ostr << Gdec("%10d",i);
  int ns   = spins();				// Number of system spins
  int nout = 4;					// # output values (rows)
  int doeta=0;					// Count each non-zero cmpt
  for(k=0; k<spinpairs(); k++)			// Loop dipoles, view asymmetry
    if(Dvec.Deta(k)) doeta++; 			// Normally 0 & NOT printed
  if(doeta) nout++;				// Add a row if eta printed
  string blanks10("          ");                // Filler for non-existing #s
  string sthdrs[5] = { "DCC (kHz)",		// Dipolar tensor component
                       "R   (A)  ", "the (deg)",// headers
                       "phi (deg)", "eta      " };
  string althdr("DCC (G)  ");			// Alternate header e- e-
  bool   altflg;				// Alternate header flag
  double cmpval = 0;				// Tensor component value
  int    cmp;					// Component index
  string cmpstr;				// Tensor component header
 
  matrix DST = SCoords.distances(1);		// Get internuclear distances
  for(i=0, l=0; i<ns-1; i++)			// Begin with 1st spin index
    {						// and loop over spin pairs
    ostr <<"\nSpin "<< Gdec(i,2);
    for(cmp=0; cmp<nout; cmp++)
      {
      cmpstr = sthdrs[cmp];			//   Set value header
      altflg = false;				//   Assume not alternate
      if(!cmp && electron(i))			//   For e- with e- we use an
        { cmpstr=althdr; altflg=true; }		//   use alternate DCC header
      ostr << "\n" << cmpstr << ":";		//   Output value header
      for(j=0; j<=i; j++) ostr << blanks10; 	//   Skip pairs not output
      for(k=l; j<ns; j++, k++)			//   Loop over spins that are
        {					//   paired with spin i
        if(nepair(i, j)) ostr << "   -------";	//   Skip e- nucleon pairs
        else
          {
          if(altflg)      cmpval=(Dvec.DCC(k))*HZ2GAUSS;
          else if(cmp==0) cmpval=(Dvec.DCC(k))*1.e-3;
          else if(cmp==1) cmpval=DST.getRe(i,j);
          else if(cmp==2) cmpval=Dvec.Dtheta(k);
          else if(cmp==3) cmpval=Dvec.Dphi(k);
          else if(cmp==4) cmpval=Dvec.Deta(k);
          ostr << Gform("%10.2f", cmpval);
          }
        }
      }
    l += ns-i-1;				// Update interaction count
    if(i != ns-1) ostr << "\n";			// Space to next spin
    }
  return ostr;
  }

 
// ----------------------------------------------------------------------------
//                    Print Shift Anisotropy Interactions
// ----------------------------------------------------------------------------

ostream& solid_sys::printCs(ostream& ostr) const
  {
  if(!Cvec.size()) return ostr;			// Bail if no CSA interactions
  int i, someCs=0, ns=spins();			// Number of spins in system
  int *csatype;
  csatype = new int[ns];				// Flags for each spin
  for(i=0; i<ns; i++)
    {
    if(electron(i)) csatype[i] =-1;		//	This if e- spin
    else if(Cvec.CSA(i))			//	This if CSA value
      {
      csatype[i] = 1;				//	This if CSA value
      someCs++;
      }
    else                  csatype[i] = 0;	//	This if 0 CSA value
    }
  if(!someCs) return ostr;			// Exit if no shift anisotropy
  ostr << "\nShift Anisotropy Interactions";	// Heres the header
  ostr << "\nSpin     :";
  for(i=0; i<spins(); i++) ostr << Gdec("%10d",i);
  ostr << "\nCSA (PPM):";			// Write CSA's for each spin
  for(i=0; i<ns; i++)				// Done horizontally
    {
    switch(csatype[i])
      {
      case -1: ostr << "   -------"; break;	// 	Don't output e- CSA
      case  0: ostr << "      0.00"; break;	// 	Spins with CSA=0 set 0
      case 1:					//	Other, output CSA value
      default:					//	(output in degrees)
        ostr << Gform("%10.2f",Cvec.CSA(i));
      }
    }   
  ostr << "\neta [0,1]:";			// Write CSA asymmetry
  for(i=0; i<ns; i++)
    {
    switch(csatype[i])
      {
      case -1: ostr << "   -------"; break;	// 	Don't output e- CSA
      case  0: ostr << "      0.00"; break;	// 	Spins with CSA=0 set 0
      case 1:					//	Other, output CSA value
      default:					//	(output is unitless)
        if(Cvec.eta(i))
          ostr << Gform("%10.2f",Cvec.eta(i));
        else ostr << "      0.00";
      }   
    }   
  ostr << "\nthe (deg):";			// Write angle down from z
  for(i=0; i<ns; i++)
    {
    switch(csatype[i])
      {
      case -1: ostr << "   -------"; break;	// 	Don't output e- CSA
      case  0: ostr << "      0.00"; break;	// 	Spins with CSA=0 set 0
      case 1:					//	Other, output CSA value
      default:					//	(output in degrees)
        ostr << Gform("%10.2f", Cvec.theta(i));
      }
    }   
  ostr << "\nphi (deg):";			// Write angle over from x
  for(i=0; i<ns; i++)
    {
    switch(csatype[i])
      {
      case -1: ostr << "   -------"; break;	// 	Don't output e- CSA
      case  0: ostr << "      0.00"; break;	// 	Spins with CSA=0 set 0
      case 1:					//	Other, output CSA value
      default:					//	(output in degrees)
        ostr << Gform("%10.2f", Cvec.phi(i));
      }
    }
  delete [] csatype;
  return ostr;
  }

// ----------------------------------------------------------------------------
//                       Print Quadrupolar Interactions
// ----------------------------------------------------------------------------

ostream& solid_sys::printQs(ostream& ostr) const
  {
  if(!Qvec.size()) return ostr;			// Bail if no quad interactions
  int i, someQs=0, ns=spins();			// Number of spins in system
  int *quadtype;
  quadtype = new int[ns];			// Flags for each spin
  for(i=0; i<ns; i++)				// Check if quad. interactions
    {						// are present & what type
    if(qn(i) < 1.0) quadtype[i] = -1; 		// Type -1 => no quad. moment
    else if(Qvec[i].QCC())			// Type  1 => non-zero QCC
      {						// Type  0 => zero QCC
      someQs++;
      quadtype[i] = 1;
      }
    else quadtype[i] = 0;
    }
  if(someQs)					// If quad. interaction present
    {						// then print them ostr
    ostr << "\nQuadrupolar Interactions";	// Heres the header
    ostr << "\nSpin     :";
    for(i=0; i<spins(); i++) ostr << Gdec("%10d",i);
    ostr << "\nQCC (kHz):";			// Write QCC's for each spin
    for(i=0; i<ns; i++)				// Done horizontally
      {
      switch(quadtype[i])
        {
        case -1: ostr << "   -------"; break; 	//	Spins with I<1 put to -
        case  0: ostr << "      0.00"; break;	// 	Spins with QCC=0 set 0
        case 1:					//	Other, output QCC value
        default:				//	(output in kHz)
          ostr<<Gform("%10.2f",Qvec[i].QCC()*1e-3);
        }
      }   
    ostr << "\nwQ  (kHz):";			// Write quadrupolar frequency
    for(i=0; i<ns; i++)				// Done horizontally, each spin
      {
      switch(quadtype[i])
        {
        case -1: ostr << "   -------"; break; 	//	Spins with I<1 put to -
        case  0: ostr << "      0.00"; break;	// 	Spins with QCC=0 set 0
        case 1:					//	Other, output QCC value
        default:				//	(output in kHz)
          ostr<<Gform("%10.2f",Qvec[i].wQ()*1e-3);
        }
      }   
    ostr << "\neta [0,1]:";			// Write quad. asymmetry
    for(i=0; i<ns; i++)
      {
      switch(quadtype[i])
        {
        case -1: ostr << "   -------"; break; 	//	Spins with I<1 put to -
        case  0: ostr << "      0.00"; break;	// 	Spins with QCC=0 set 0
        case 1:					//	Other, output QCC value
        default:				//	(output in kHz)
          if(Qvec[i].eta())
            ostr << Gform("%10.2f",Qvec[i].eta());
          else ostr << "      0.00";
        }
      }   
    ostr << "\nthe (deg):";			// Write angle down from z
    for(i=0; i<ns; i++)
      {
      switch(quadtype[i])
        {
        case -1: ostr << "   -------"; break; 	//	Spins with I<1 put to -
        case  0: ostr << "      0.00"; break;	// 	Spins with QCC=0 set 0
        case 1:					//	Other, output QCC value
        default:				//	(output in kHz)
          ostr << Gform("%10.2f", Qvec[i].theta());
        }
      }   
    ostr << "\nphi (deg):";			// Write angle over from x
    for(i=0; i<ns; i++)
      {
      switch(quadtype[i])
        {
        case -1: ostr << "   -------"; break; 	//	Spins with I<1 put to -
        case  0: ostr << "      0.00"; break;	// 	Spins with QCC=0 set 0
        case 1:					//	Other, output QCC value
        default:				//	(output in kHz)
          ostr << Gform("%10.2f", Qvec[i].phi());
        }
      }   
    }
  delete [] quadtype;
  return ostr;
  }
 
// ----------------------------------------------------------------------------
//                     Print The Electron G Interactions
// ----------------------------------------------------------------------------

/* This prints out a listing of the electron G interactions in the system.
   Typically there will be such an interaction for each spin although those
   associated with a nucleon will be NULL.  Here we print interactions
   only for electrons that have a non-zero anisotropic component. The values
   are output in Cartesian terms by default.  The input flag allows for this
   (pf=0), spherical terms (pf>0), or both (pf<0)                            */

ostream& solid_sys::printGs(ostream& ostr, int pf) const
  {
  if(!Gvec.size()) return ostr;			// Bail if no G interactions
  int i, someGs=0, ns=spins();			// Number of spins in system
  int *etype;					// Flags if electron/nucleon
  etype = new int[ns];				// Allocate flags each spin
  for(i=0; i<ns; i++)				// Loop over the spins and see
    {						// if electrons or nucleons
    if(!spin_sys::isotope(i).electron()) 	//   This if spin is nucleon
      etype[i] =-1;				//   & set flag not to print
    else if(Gvec[i].iso())			//   This if G delzz value is
      { 					//   non-zero. Int this case
      etype[i] = 1; 				//   we flag it is printed
      someGs++; 				//   and count the # we print
      }
    else					//   This if G delzz is zero
      etype[i] = 0;				//   Flag we don't print it
    }
  if(!someGs) return ostr;			// Exit if no G interactions
  ostr << "\nElectron G Interactions";		// Heres the header
  ostr << "\nSpin     :";
  for(i=0; i<spins(); i++) ostr << Gdec("%10d",i);

//             These Are The Spatial Tensor Cartesian Components

  if(pf <=0)
    {
    ostr << "\nGxx      :";			// Write gxx's for each spin
    for(i=0; i<ns; i++)				// Done horizontally
      {
      switch(etype[i])
        {
        case -1: ostr << "   -------"; break;	//   Don't output nucleon gxx
        case  0: ostr << "   -------"; break;	//   Electron delzz=0 no output
        case  1:				//   Other, output gxx value
        default:				//	(This is unitless)
          ostr << Gform("%10.5f",Gvec[i].gxx());
        }
      }   
    ostr << "\nGyy      :";			// Write gyy's for each spin
    for(i=0; i<ns; i++)				// Done horizontally
      {
      switch(etype[i])
        {
        case -1: ostr << "   -------"; break;	//   Don't output nucleon gyy
        case  0: ostr << "   -------"; break;	//   Electron delzz=0 no output
        case  1:				//   Other, output gyy value
        default:				//	(This is unitless)
          ostr << Gform("%10.5f",Gvec[i].gyy());
        }
      }   
    ostr << "\nGzz      :";			// Write gzz's for each spin
    for(i=0; i<ns; i++)				// Done horizontally
      {
      switch(etype[i])
        {
        case -1: ostr << "   -------"; break;	//   Don't output nucleon gxx
        case  0: ostr << "   -------"; break;	//   Electron delzz=0 no ouptut
        case  1:				//   Other, output gzz value
        default:				//	(This is unitless)
          ostr << Gform("%10.5f",Gvec[i].gzz());
        }
      }   
    }   
//               These Are The Spatial Tensor Spherical Components

  if(pf)
    {
    ostr << "\ndelzz    :";			// Write G interact. anisotropy
    for(i=0; i<ns; i++)
      {
      switch(etype[i])
        {
        case -1: ostr << "   -------"; break;	//  Don't output nucleon eta
        case  0: ostr << "      0.00"; break;	//  Electron delzz=0 set 0
        case 1:					//  Other, output eta value
        default:				//    (output is unitless)
ostr << "dont' know";
//sosik
//            ostr << Gform("%10.2f",Gvec[i].delz());
        }   
      }   
    ostr << "\neta [0,1]:";			// Write G interact. asymmetry
    for(i=0; i<ns; i++)
      {
      switch(etype[i])
        {
        case -1: ostr << "   -------"; break;	//  Don't output nucleon eta
        case  0: ostr << "      0.00"; break;	//  Electron delzz=0 set 0
        case 1:					//  Other, output eta value
        default:				//    (output is unitless)
          if(Gvec.eta(i))
            ostr << Gform("%10.2f",Gvec[i].eta());
          else ostr << "      0.00";
        }   
      }   
    }   
//     These Spatial Tensor Values Apply To Both Catesian & Spherical Output

  ostr << "\nthe (deg):";			// Write angle down from z
  for(i=0; i<ns; i++)
    {
    switch(etype[i])
      {
      case -1: ostr << "   -------"; break;	//   Don't output nucleon theta
      case  0: ostr << "   -------"; break;	//   If delzz=0 don't output
      case 1:					//   Other, output theta 
      default:					//   value (degrees)
        ostr << Gform("%10.2f", Gvec[i].theta());
      }
    }   
  ostr << "\nphi (deg):";			// Write angle over from x
  for(i=0; i<ns; i++)
    {
    switch(etype[i])
      {
      case -1: ostr << "   -------"; break;	//   Don't output nucleon theta
      case  0: ostr << "   -------"; break;	//   If delzz=0 don't output
      case 1:					//   Other, output phi
      default:					//   value (degrees)
        ostr << Gform("%10.2f", Gvec[i].phi());
      }
    }
  delete [] etype;
  return ostr;
  }
 
// ----------------------------------------------------------------------------
//                     Print The Hyperfine Interactions
// ----------------------------------------------------------------------------

ostream& solid_sys::printHFs(ostream& ostr, int pf) const
  {
  if(!HFvec.size())   return ostr;	// Bail if no HF interactions
  if(!HFvec.nonzero()) return ostr;	// Bail if no non-zero interactions

//      We Have Hyperfine Interactions So We Print Them By Spin Pair

  ostr << "\nHyperfine Interactions";		// Start with the header
  int i,j;					// Temp spin indices
  int k, l;					// Temp interaction indices
  ostr << "\nSpin     :"; 			// Write row of spin labels
  for(i=0; i<spins(); i++)			// for column indices
    ostr << Gdec("%10d",i);
  int ns = spins();                             // Number of system spins
  string blanks10("          ");                // Filler for non-existing #s
  string sthdrs[5]={ "Axx (G)  ", "Ayy (G)  ",	// Hyperfine tensor component
                     "Azz (G)  ", 		// headers
                     "the (deg)", "phi (deg)"};	
  double cmpval=0;                              // Tensor component index
  for(i=0, l=0; i<ns-1; i++)			// Begin with 1st spin index
    {						// and loop over spin pairs
    ostr <<"\nSpin "<< Gdec(i,2);
    for(int cmp=0; cmp<5; cmp++)		// Loop components output
      {						// Makes 5 lines with a column
      ostr << "\n" << sthdrs[cmp] << ":";
      for(j=0; j<=i; j++) ostr << blanks10;	//   Loop spins paired j<=i
      for(k=l; j<ns; j++, k++)			//   Loop spins paired with i
        {					//   interaction k
        if(!nepair(i,j))
          ostr << "   -------";
        else
          {
          if(cmp==0)      cmpval=HFvec[k].hxx();
          else if(cmp==1) cmpval=HFvec[k].hyy();
          else if(cmp==2) cmpval=HFvec[k].hzz();
          else if(cmp==3) cmpval=HFvec[k].theta();
          else if(cmp==4) cmpval=HFvec[k].phi();
          ostr << Gform("%10.2f", cmpval);
          }
        }
      }
    ostr << "\n";
    l += ns-i-1;
    }
  return ostr;
  }

// ----------------------------------------------------------------------------
//                  Print The Solid State Spin System
// ----------------------------------------------------------------------------

ostream& solid_sys::print(ostream& ostr) const
  {
  spin_system::print(ostr);			// 1st write spin_system
  if(!spins()) return ostr;			// Exit now if no spins
  printPs(ostr);				// Print spin coordinates
  printDs(ostr); ostr << "\n";			// Print dipolar interactions 
  printCs(ostr); ostr << "\n";			// Print CSA interactions
  printQs(ostr); ostr << "\n";			// Print quad interactions
  printGs(ostr); ostr << "\n";			// Print e- G interactions
  printHFs(ostr); 				// Print hyperfine interacts.
  return ostr;
  }

ostream& operator<< (ostream& out, const solid_sys& ssys)
  { return ssys.print(out); }

// ____________________________________________________________________________
// I                         BASE SPIN SYSTEM DEALINGS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                These Are All Inherited As Is From Class spin_sys
// ----------------------------------------------------------------------------

/*
inline int    spins() const;			Get number of spins
inline int    spinpairs() const;		Get number of spin pairs
inline int    HS() const;			Get system Hilbert space             
inline int    HS(int spin) const;		Get spin Hilbert space
inline const  Isotope& isotope(int) const;	Get spin isotope(eg 23Na)
inline double weight(int) const;                Get spin atomic weight
inline string symbol(int spin) const;		Get spin symbol(eg 23Na)
inline double qn(int spin) const;		Get spin I value (eg .5) 
inline double qn() const;			Get system Fz value
inline string element(int spin) const;		Get spin type (eg Carbon)
inline string momentum(int spin) const;		Get spin momentum(eg 1/2)
inline string momentum() const;			Get system momentum -1 -1
inline double gamma(int spin) const;		Get gyromag. ratio s  T
inline double gamma(const string& iso) const;	Get gyromag. ratio of iso
inline row_vector qState(int state) const;	Get vector of state Fz
inline matrix     qStates() const;		Get array of spin Iz
inline double     qnState(int state) const;
inline col_vector qnStates() const;
inline row_vector qnDist() const;		Get vector of state dist.
inline row_vector CoherDist() const;		Get vector of coher. dist.
inline int    homonuclear() const;		Test homonuclear
inline int    heteronuclear() const;		Test heteronuclear
       int    spinhalf() const;			Test all spin 1/2
       int    electrons() const;		Test for electrons
       int    nepair(int i, int j) const;	Test for nucleus/e- pair
       int    enpair(int i, int j) const;	Test for e-/nucleus pair
       int    pairidx(int i, int j) const;	Get spin pair index
       bool   electron(int i) const;		Test if spin i electron
       bool   nucleon(int i) const;		Test if spin i nucleon
inline int    isotopes() const;			Get isotope count
inline string isotopes(int idx) const;		Get isotope spin count
inline int    isotopes(const string& I) const;	Test for isotope type
inline void   flags(int TF);			Set spin flags
inline void   flag(int spin, int TF);		Set spin flag
inline void   flag(const string& isoin, int TF);Set isotope flags
inline void   flag(const Isotope& isoin,int TF);Set isotope flags
inline int    flag(int spin) const;		Set spin flag
void get_flags(int* intarray) const;		Get array of spin flags
void set_flags(int* intarray, int TF) const;	Set spin flags
void set_flag(int* intarray, int spin, int TF) const;
void set_flag(int* intarray, const string& isoin, int TF) const;
void set_flag(int* intarray, const Isotope& isoin, int TF) const;
inline void name(const string& name);		Set system name
inline const  string& name(int i=-1) const;	Get system name
void warnings(int warnf);			Set warning level
int warnings() const;				Get warning level
string IsoDefault();				Get default isotope
void IsoDefault(const string& DI);		Set default isotope
int  getSpins(const ParameterSet& pset);	Read number spins
void setName(const ParameterSet& pset);	Read system name
void setIs(const ParameterSet& pset);	Read system isotopes         */

// ----------------------------------------------------------------------------
//    These Are All ReDefined Beyond Their Functionality In Class spin_sys
// ----------------------------------------------------------------------------

/*
virutal void     isotope(int, const Isotope&);	Set isotope type 
virtual void     isotope(int, const string&);	Set isotope type
virtual void     PSetAdd(ParameterSet& pset, int idx=-1) const;
virtual int      write(const string &filename, int ix=-1, int wn=2) const;
virtual int      write(ofstream& ofstr, int idx=-1, int warn=2) const;
virtual void     read(const string &filename);	Read spin system
virtual string   ask_read(int argc, char* argv[], int argn);
virtual ostream& print(ostream& out) const;	Print spin system           */

// ----------------------------------------------------------------------------
//          These Are The Class spin_sys Redefined Functions
// ----------------------------------------------------------------------------
 
        // Input                ssys    : Solid spin system (this)
        //                      Iso     : Isotope type
        //                      symbol  : Isotope type (e.g. 1H, 13C,...)
        // Output               none    : Solid spin system spin isotope
        //                                type is switched to Iso


void solid_sys::isotope(int spin, const string& symbol)

  {
  spin_sys::isotope(spin, symbol);		// Use base function
  if(symbol == "e-")
    {
    }
  else
    {
    }
// sosiz
  ResetSOps(spin);				// Reset the spin operators
  }


void solid_sys::isotope(int spin, const Isotope& Iso)

        // Input                ssys    : Solid spin system (this)
        //                      Iso     : Isotope type
        // Output               none    : Solid spin system spin isotope
        //                                type is switched to Iso

  {
  spin_sys::isotope(spin, Iso);			// Use base function
  ResetSOps(spin);				// Reset the spin operators
  }

const Isotope& solid_sys::isotope(int spin) const
  { return spin_sys::isotope(spin); }

#endif						// solid_sys.cc
