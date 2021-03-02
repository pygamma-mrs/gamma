/* IntDip.cc ****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Dipole-Dipole Interaction 	         Implementation		**
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
** A variable of type IntDip represents a dipole-dipole interaction	**
** defined for a particular nuclear spin pair. The interaction is 	**
** maintained in the form of spatial and spin tensors which are stored 	**
** in irreducible spherical format. 					**
**                                                                      **
** 1.) rank = 2         2.) symmetric           3.) eta = [0, 1]        **
**                                                                      **
** The tensor will internally be stored in irreducible spherical form   **
** at any chosen orientation.  The tensor principal axes will also be   **
** maintiained as well as a set of Euler angles that indicate the PAS   **
** orientation relative to some common axes.                            **
**                                                                      **
** Although not internal to this class, the interaction will blend with **
** a rank 2 spin tensor for the formation of an oriented dipolar    	**
** Hamiltonian.                                                         **
**                                                                      **
** In addition to the functions contained in the base spatial tensor    **
** class, this IntDip provides functions for building up the tensor	**
** and accessing the tensor elements from a "dipolar" standpoint. 	**
**                                                                      **
** The following defintions are used herein (Auv GAMMA normlized):	**
**                                                                      **
** 1.) PAS: |Azz| >= |Ayy| >= |Axx|					**
**                                                                      **
** 2.) PAS: Azz=2C, eta=(Axx-Ayy)/Azz, Axx=C(eta-1), Ayy=-C(1+eta)	**
**            								**
**                                                                      **
*************************************************************************/

#ifndef   IntDip_cc_			// Is file already included?
#  define IntDip_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <IntRank2/IntDip.h>		// Include interface definition
#include <Basics/Gconstants.h>		// Include PI, HBAR, & other constants
#include <Basics/Isotope.h>		// Include isotopes
#include <Basics/ParamSet.h>		// Include parameter sets
#include <HSLib/SpinOpSng.h>		// Include 1 spin operators
#include <Matrix/row_vector.h>		// Include row_vectors
#include <Matrix/matrix.h>		// Include matrices
#include <Level1/coord.h>		// Include coordinates
#include <Basics/Gutils.h>		// Include query parameter 
#include <string>
#include <fstream>
#include <Basics/StringCut.h>		// Using Gdec functions
#if defined(_MSC_VER) || defined(__SUNPRO_CC)				// If we are not using MSVC++
 #pragma warning (disable : 4786)       //   Kill STL namelength warnings
#endif

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ STL lists
using std::vector;			// Using libstdc++ STL vectors
using std::ofstream;			// Using libstdc++ output file streams
using std::ostream;			// Using libstdc++ output streams

const double MU   = 1.e-7;		// [mu/(4*pi)] J-sec -C  -m

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i             CLASS DIPOLE DIPOLE INTERACTION ERROR HANDLING
// ____________________________________________________________________________
 

/*       Input                D	      : Dipolar interaction (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void IntDip::IDerror(int eidx, int noret) const
  {
  string hdr("Dipolar Interaction");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr,"Program Aborting.....",        noret); break;// (0)
    case 2: GAMMAerror(hdr,"Problems During Construction", noret); break;// (2)
    case 5: GAMMAerror(hdr,"DCC Must Be >0 This Spin Pair",noret); break;// (5)
    case 6: GAMMAerror(hdr,"Using Absolute Value of DCC",  noret); break;// (6)
    case 7: GAMMAerror(hdr,"DCC Must Be <0 This Spin Pair",noret); break;// (7)
    case 8: GAMMAerror(hdr,"Negating Value of Input DCC",  noret); break;// (8)
    case 11:GAMMAerror(hdr,"I,S Must Be Multiples Of 1/2", noret); break;// (11)
    case 13:GAMMAerror(hdr,"Can't Construct From Par File",noret); break;// (13)
    case 14:GAMMAerror(hdr,"Info Troubles In Param File",  noret); break;// (14)
    case 16:GAMMAerror(hdr,"Setting Spin @ Default Iz=1/2",noret); break;// (16)
    case 17:GAMMAerror(hdr,"Setting Spin Type To Proton",  noret); break;// (17)
    case 20:GAMMAerror(hdr,"Can't Write To Parameter File",noret); break;// (20)
    case 21:GAMMAerror(hdr,"Can't Read From Param. File",  noret); break;// (21)
    case 22:GAMMAerror(hdr,"Can't Write To Output FileStr",noret); break;// (22)
    case 23:GAMMAerror(hdr,"Cannot Output Parameters",     noret); break;// (23)
    case 44:GAMMAerror(hdr,"Can't Read Coords. In PSet",   noret); break;// (44)
    case 45:GAMMAerror(hdr,"Same Coordinate For 2 Spins!", noret); break;// (45)
    case 46:GAMMAerror(hdr,"Found Only 1 Coord. In ParSet",noret); break;// (46)
    case 47:GAMMAerror(hdr,"Can't Determine Spin Isotopes",noret); break;// (47)
    case 49:GAMMAerror(hdr,"Improper Spin Designation",    noret); break;// (49)
    case 50:GAMMAerror(hdr,"Distance-Isotope Mismatch...", noret); break;// (50)
    case 51:GAMMAerror(hdr,"Can't Find Spin I Quantum #",  noret); break;// (51)
    case 52:GAMMAerror(hdr,"Can't Find Spin I Quantum #",  noret); break;// (52)
    case 60:GAMMAerror(hdr,"Use A Hyperfine Interaction?", noret); break;// (60)
    default:GAMMAerror(hdr,"Unknown Error",                noret); break;
    }
  }

volatile void IntDip::IDfatal(int eidx) const
  {
  IDerror(eidx, eidx);				// Output error message
  if(eidx) IDerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }


/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Function pname
      5     Please Use Class Member Function pname
   default  Unknown Error - pname                                           */

void IntDip::IDerror(int eidx, const string& pname, int noret) const
  {
  string hdr("Dipolar Interaction");
  string msg;
  switch(eidx)
    {
    case 44:                                                   // (44)
      msg = string("Problems Setting Interaction ")
          + string("Index ") + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 45:                                                   // (45)
      msg = string("Problems Setting Interaction ")
          + string("Between Spins ") + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 46:                                                   // (46)
      msg = string("Cannot Glean Interaction From Parameter ")
          + string("File, Spins ") + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 103:                                                   // (103)
      msg = string("Odd Asymmetry Value Of ")
          + pname + string(" Specified?");
      GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;
    }
  }

volatile void IntDip::IDfatal(int eidx, const string& pname) const
  {
  IDerror(eidx, pname, eidx);			// Output error message
  if(eidx) IDerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii    DIPOLE DIPOLE INTERACTION SETUP FUNCTIONS USING PARAMETER SETS
// ____________________________________________________________________________

/* These functions set up specific aspects of a dipolar interaction.  Since
   the various interaction parameters interact, the functions MUST be private
   because their misuse could produce an inconsistent interaction.
   
   The goal here is quite simple.  We must determine the following set of
   values for each dipolar interaction: { Iqn,Sqn,DCC,eta,alpha,beta,gamma }
   Complexity arises because we allow that a variety of parameters may be used
   to define these values. Additionally we allow that some parameters may be
   zero (e.g. eta) and that defaults may automatically be used. But, the end
   result remains the same, we seek the set of values that define a dipolar
   interaction in GAMMA.                                                     */

// ----------------------------------------------------------------------------
//	                 Complete Dipolar Interaction 
// ----------------------------------------------------------------------------

/* These functions will try and get all the parameters required to define a
   dipolar interaction: { Iqn,Sqn,DCC,eta,alpha,beta,gamma }. Although the 1st
   function is the generic means to accomplishing this, we will immediately
   branch into two categories: 1.) Defining the interaction with a single
   interaction index & 2.) Defining the interaction with two spin indices.

           Input                D       : Dipolar interaction (this)
                                pset    : A parameter set
                                idxI    : Index of first spin or interaction
                                idxS    : Index of first spin or unused (-1)
	  			warn 	: Warning level
	  				   0 - no warnings
	  				   1 - warnings
	  				  >1 - fatal warnings
           Output               TF	: Dipolar interaction is set 
                                          from parameters in pset
	  					TF = 0 couldn't read it
	  				        TF = 1 read with iso & coord
	  				        TF = 2 read with iso & DCC   */

bool IntDip::getDI(const ParameterSet& pset, 
         double& Iqn, double& Sqn, double& dcc, double& eta, EAngles& EA,
                                                  int idxI, int idxS, int warn)
  {
  if(idxS == -1) return getDI1(pset, Iqn, Sqn, dcc, eta, EA, idxI,       warn);
  else           return getDI2(pset, Iqn, Sqn, dcc, eta, EA, idxI, idxS, warn);
  }

bool IntDip::getDI1(const ParameterSet& pset, 
         double& Iqn, double& Sqn, double& dcc, double& eta, EAngles& EA,
                                                            int idxI, int warn)
  {				
//                First Determine The Two Spin Quantum Numbers
//     ( Default Is I=S=1/2, This Sets Dipolar Spin Tensor Hilbert Space )

  bool TFI = getIqn(pset, "Iqn", Iqn, idxI, 0);		// Try to get Iqn using
  if(!TFI) TFI = getIqn(pset, "DIqn", Iqn, idxI, 0);	// name Iqn or DIqn
  bool TFS = getIqn(pset, "Sqn", Sqn, idxI, 0);		// Try to get Sqn using
  if(!TFS) TFS = getIqn(pset, "DSqn", Sqn, idxI, 0);	// name Sqn or DSqn
  if(!TFI)						// If no Iqn found we
    {							// set a default value
    Iqn = 0.5;						// of 1/2 & warn
    if(warn) IDerror(51, 1);
    }
  if(!TFS)						// If no Iqn found we
    {							// set a default value
    Sqn = 0.5; 						// of 1/2 & warn
    if(warn) IDerror(52, 1);
    }
  
//       Try To Directly Read Dipolar Coupling, Asymmetry, Orientation
//       Dipolar Interaction Via { Iqn, Sqn, DCC(i), Deta(i), DEAs(i) }

//  1.) If DCC has been specified, these parameters will be used
//  2.) We don't mind that eta is not set, default will be zero
//  3.) Iqn & Sqn were set in the 1st section of this function
//  4.) We don't mind that no orientation is set, default is PAS
//  5.) Orientation set with either an Euler angle set or 3 individual angles

  if(getDCC(pset, dcc, idxI, -1, 0))			// 1.) Try and read in
    {							//     DCC, eta, Orient
    string Pbase("D");
    getAeta(pset, Pbase, eta, idxI, -1, 0);
    getOrientation(pset,Pbase,EA,idxI,-1,0);
    return true;
    }
  
//      Try To Directly Read Dipolar Cartesian Spatial Tensor Components
//     Dipolar Interaction Via { Iqn, Sqn, Dxx, Dxy, Dxz, Dyy, Dyz, Dzz }

//  1.) Duv values are used to set the dipolar coupling & asymmetry, DCC & eta
//  2.) If any Duv off-diagonals (u != v) set, D array sets orientation also
//  3.) If no Duv off-diagonals, orientation set with specified Euler angles 
//  4.) The input Duv values are taken to be in kHz
//  5.) Orientation set with either an Euler angle set or 3 individual angles

  coord DiDzDe;					// For Diso, Dzz, Deta
  if(getACart(pset,"D",DiDzDe,EA,idxI,-1,0))	// Try & use Cart. components
    {
    dcc = DiDzDe.y()*1.e3;			//  Set dipolar coupling 
    eta = DiDzDe.z();				//  Set the eta value
    return true;
    }
  
  return true;
  }

bool IntDip::getDI2(const ParameterSet& pset, 
         double& Iqn, double& Sqn, double& dcc, double& eta, EAngles& EA,
                                                  int idxI, int idxS, int warn)
  {
//                First Determine The Two Spin Quantum Numbers
//     ( Default Is I=S=1/2, This Sets Dipolar Spin Tensor Hilbert Space )

  string II, IS;				// String for isotope names
  Isotope ISI, ISS;				// Isotopes for spins
  string pn("D");				// Parameter name base
  bool TFI = false;				// Flag if we know isotopes
  if(getIsos(pset,idxI,idxS,II,IS,0))		// 1st try for isotope names
    {						// If successful, check them
    if(!SpinCheck(II,IS)) return false;		//    Insure valid isotopes
    ISI = Isotope(II);				//    An isotope of type II
    ISS = Isotope(IS);				//    An isotope of type IS
    if(!SpinCheck(ISI,ISS,TFI,0)) return false;	//    Disallow nucleus/e- pair
    TFI = true;					//    We know both isotopes
    Iqn = ISI.qn();				//    Know spin I quantum #
    Sqn = ISS.qn();				//    Know spin S quantum #
    }
  else if(!getIqns(pset,pn,Iqn,Sqn,idxI,idxS,0))// 2nd try for spin quant. #s
    { Iqn = 0.5; Sqn = 0.5; }			// Use default if not able
  
//            Try To Read Spin Coordinates If We Know Isotope Types
//       Dipolar Interaction Via { Iso(i), Iso(j), Coord(i), Coord(j) }

  if(TFI)					// If we know the isotopes
    if(getDI2(pset, dcc, eta, EA,		// we allow for spin coords.
                    ISI, ISS, idxI, idxS, warn))// The combination will set
      return true;				// DCC & orientation, eta=0
  
//        Try To Directly Read Dipolar Coupling, Asymmetry, Orientation
//     Dipolar Interaction Via { Iqn, Sqn, DCC(i,j), Deta(i,j), DEAs(i,j) }

//  1.) If DCC has been specified, these parameters will be used
//  2.) We don't mind that eta is not set, default will be zero
//  3.) Iqn & Sqn were set in the 1st section of this function
//  4.) We don't mind that no orientation is set, default is PAS
//  5.) Orientation set with either an Euler angle set or 3 individual angles

  if(getDCC(pset, dcc, idxI, idxS, 0))		// If DCC has been specified
    {						// we look for eta & the
    string Pbase("D");				// orientation
    getAeta(pset, Pbase, eta, idxI, idxS, 0);
    getOrientation(pset,Pbase,EA,idxI,idxS,0);
    return true;
    }
  
//      Try To Directly Read Dipolar Cartesian Spatial Tensor Components
//     Dipolar Interaction Via { Iqn, Sqn, Dxx, Dxy, Dxz, Dyy, Dyz, Dzz }

//  1.) Duv values are used to set the dipolar coupling & asymmetry, DCC & eta
//  2.) If any Duv off-diagonals (u != v) set, D array sets orientation also
//  3.) If no Duv off-diagonals, orientation set with specified Euler angles 
//  4.) The input Duv values are taken to be in kHz
//  5.) We don't mind that no orientation is set, default is PAS
//  6.) If used, Euler angles by either a set or 3 individual angles

  coord DiDzDe;					// For Diso, Dzz, Deta
  if(getACart(pset,"D",DiDzDe,EA,idxI,idxS,0))	// Try & use Cart. components
    {
    dcc = DiDzDe.y()*1.e3;			//  Set dipolar coupling 
    eta = DiDzDe.z();				//  Set the eta value
    return true;
    }
  
  return false;					// Tried everything, but we
  }						// cant read the interaction

//                  { Iso(i), Iso(j), Coord(i), Coord(j) }

bool IntDip::getDI2(const ParameterSet& pset, 
          double& dcc, double& eta, EAngles& EA,
          const Isotope& ISI, const Isotope& ISS, int idxI, int idxS, int warn)
  {
  coord ptI, ptS;				// If we know isotope types &
  if(getCoords(pset,ptI,ptS,idxI,idxS,0))	// can read spin coordinates,
    { 						// can set DCC & orientation
    dcc = CheckDCC(ISI, ISS, ptI, ptS);		// Calculate dipolar coupling
    eta = 0.0;					// No asymmetry in this case
    double phi = ptI.phi(ptS);			// Euler angle alpha 
    double the = ptI.theta(ptS);		// Euler angle beta
    EA = EAngles(phi,the,0.0);			// Euler angles set (orient)
    return true;				// Done, know full interaction
    }
  return false;					// Here if we failed to get
  }						// both spin coordinates

// ----------------------------------------------------------------------------
//         Get Dipolar Coupling Constant Value From A Parameter Set
// ----------------------------------------------------------------------------

/*                                    gamma * gamma * hbar
                        DELZZ = DCC =      I       S
                                      ---------------------
                                                3
                                               r

           Input                D       : Dipolar interaction (this)
                                pset    : A parameter set
                                DCC     : Dipolar coupling constant (Hz)
                                idxI    : Index of first spin or interaction
                                idxS    : Index of first spin
                                warn    : Warning level
                                           0 - no warnings
                                           1 - warnings
                                          >1 - fatal warnings
           Output               TF      : Return is true if the dipolar
                                          coupling constant has been found
                                          from parameters in pset. The value of
                                          DCC is set in Hz prior to exit.
	   Note				: There are two alternative ways that
					  parameters may define the DCC value.
					  1.) Through direct specification of
					      the dipolar Cartesian spatial
                                              tensor { Duv }
                                          2.) Though use of isotope labels and
                                              coordinates. This provides the
                                              interspin distance and the gyro-
                                              magnetic ratios needed for DCC

                         Currently Allowed DCC Parameters

          DCC,    DCC(#),    DCC(#,#)      - Dipolar Coupling in kHz
          DCCkHz, DCCkHz(#), DCCkHz(#,#)   - Dipolar Coupling in kHz
          DCCKHz, DCCKHz(#), DCCKHz(#,#)   - Dipolar Coupling in kHz
          DCCHz,  DCCHz(#),  DCCHz(#,#)    - Dipolar Coupling in Hz          */

bool IntDip::getDCC(const ParameterSet& pset, double& dcc, 
                                                 int idxI, int idxS, bool warn)
  {
  dcc = 0;					// Zero the value for starters
  string Nidx(""); 				// Parameter name suffix
  if(idxI != -1)				// Use suffix if idxI not -1
    {						// If using an interaction
    Nidx += string("(") + Gdec(idxI);		// index (idxS = -1) the suffix
    if(idxS != -1)				// will be (idxI) whereas if
      Nidx += string(",") + Gdec(idxS);		// using two spin indices the
    Nidx += string(")");			// suffix will be (idxI,idxS)
    }
  string Dnames[4] = { "DCC",    "DCCkHz", 	// Dipolar interaction param.
                       "DCCKHz", "DCCHz" }; 	// base names to set delzz.
  string pname, pstate;				// Strings for parameter
  ParameterSet::const_iterator item;		// A pix into parameter list
  SinglePar par;				// Single parameter for use
  for(int i=0; i<4; i++)			// Loop over possible delzz
    {                                           // parameters names
    pname = Dnames[i] +  Nidx;                  //      Parameter name
    item = pset.seek(pname);			//      Seek parameter in pset
    if(item != pset.end())			//      If it's been found
      {                                         //      parse the parameter
      (*item).parse(pname,dcc,pstate);		//      info and set delzz
      switch(i)                                 //      based on input type
        {
        case 0:                                 //      Here if DCC (kHz)
        case 1:					//	Here if DCCkHz
        case 2:					//	Here if DCCKHz
        default:
        dcc *= 1.e3;	
          break;
        case 3:                                 //      Here if DCCHz
          break;
        }
      return true;				//      Flag we've found this
      }
    }
  return false;
  }





// ********** 	     Set Up/Read The Whole Dipolar Interaction       **********


        // Input                D       : Dipolar interaction (this)
        //                      pset    : A parameter set
        //                      idxI    : Index of first  spin or interaction
        //                      idxS    : Index of second spin or unused
	//			warn 	: Warning level
	//				   0 - no warnings
	//				   1 - warnings
	//				  >1 - fatal warnings
        // Output               TF	: Dipolar interaction is set 
        //                                from parameters in pset
	//					TF = 0 couldn't read it
	//				        TF = 1 read with iso & coord
	//				        TF = 2 read with iso & DCC
// sosiz - still need to work out warnings here
//         the int TF part is not working, is it needed?

int IntDip::setDI(const ParameterSet& pset, int idxI, int idxS, int warn)
  {
  double Iz, Sz;
  double dcc, eta;
  EAngles EA;
  if(getDI(pset, Iz, Sz, dcc, eta, EA, idxI, idxS, warn))
    {
    _DCC = dcc;						// Set dipolar coupling
    double X = xi();					// Get interaction con.
    IntRank2::operator=(IntRank2(Iz,Sz,X,eta,EA));	// Generic interaction
    setT20wh();						// Added T2 component
    return 1;
    }
  return 0;
  }
 
// ----------------------------------------------------------------------------
//          Set Dipolar Interaction Spatial and Spin Tensor Components
// ----------------------------------------------------------------------------
 
/* ----------------------------------------------------------------------------
   These bypass the generalized GAMMA spatial (class space_T) and spin 
   (class spin_T) tensors, instead relying on the irreducible rank 2 
   base classes IntRank2A and IntRank2T which are more suitable for the
   treatment of dipole-dipole interactions.  Both are scaled such that they
   are interaction independent and adhere to the relative scaling between
   components of true spherical tensors.  The rank 2 spatial tensors are
   scaled such that, when rotated they are normalized rank 2 spherical
   harmonics in a symmetric case (eta=0), usually true in dipolar interactions.
   --------------------------------------------------------------------------*/

// void IntDip::setAs()		INHERITED	Set spatial tensor components

// void IntDip::setTs()         INHERITED	Done in spin tensor base class

        // Input                D	: Dipolar interaction (this)
        // Output               none    : Dipolar interaction spin
        //                                components are generated
	// Note				: T   = T0, T   = T1, T   = Tm1
	//				   2,0       2,1       2,-1
	//
	//     				        T   = T2 , T    = Tm 
	//				         2,2        2,-2
	// Note				: These should match the rank 2
	//				  spin tensors in class spin_T
	// Note				: An overall scaling factor could be
	//				  used on (all) these & should account
	//				  for any variations with literature.

//                                          +
//			                 m  |
// This Rule Appplies:        T    = (-1)  A
//			       2,m          2,-m

void IntDip::setT20wh()

        // Input                D	: Dipolar interaction (this)
        // Output               none    : Dipolar interaction weak
        //                                heteronuclear spin tensor 
	//				  components is generated
	//
	//     				       T   | = T20wh 
	//				        2,0|     
	//                                          weak heteronuclear
	//
	// Note				: This should be called during
	//				  all constructors after IntRank2T
	//				  has been generated

//          1/2                                1/2             1/2
//       [1]                      weak      [1]             [4]
// T   = |-| * [3I S - I.S]   ------------> |-| * [2I S ] = |-| I S
//  2,0  [6]      z z         heteronuclear [6]      z z    [6]  z z

  {                                                                             
  matrix IE = Ie(Ival);                         // The operator Ie
  matrix SE = Ie(Sval);                         // The operator Se
  matrix IZ = tensor_product(Iz(Ival), SE);     // The operator Iz
  matrix SZ = tensor_product(IE, Iz(Sval));     // The operator Sz
  T20wh  =  (2.*IZ*SZ)/sqrt(6.);		// T20  = (2IzSz)/sqrt(6)
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A          DIPOLE DIPOLE INTERACTION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ---------------------------------------------------------------------------- 
//                  Simple Constructors That Don't Do Much
// ---------------------------------------------------------------------------- 

IntDip::IntDip()                 : IntRank2()   { _DCC = 0; }
IntDip::IntDip(const IntDip &D1) : IntRank2(D1)
  {
  _DCC = D1._DCC;			// Copy dipolar coupling
  T20wh = D1.T20wh;			// Copy T20 weak heteronuclear
  }

// ----------------------------------------------------------------------------
//              Direct Constructors Using Spherical Components
// ----------------------------------------------------------------------------

/* Here we need to know the quantum numbers of the two spins involved, the
   dipolar coupling, and the interaction asymmetry (usually 0). Shoud a non-
   zero asymmetry actually be needed it is unitless and in the range [0,1]

           Input                D	: Dipolar interaction (this)
	  			II	: Spin I isotope type
	  			IS	: Spin S isotope type
	  			Iqn	: Spin I quantum number (e.g. 1.5)
	  			Sqn	: Spin S quantum number (e.g. 0.5)
                                DCC     : Tensor dipolar coupling value (Hz)
                                eta     : Tensor asymmetry value (default 0)
	   Output		none    : Dipolar interaction constructed
	   Note				: A fatal error will result if
	  				  either IsoI & IsoS are a mix of an
	  				  electron and a nucleon
	   Note				: A non-fatal error will occur if
	  				  the sign of DCC doesn't match the
	  				  sign of gamma(IsoI)*gamma(IsoS)
	   Note				: Here DCC=delzz                     */

IntDip::IntDip(const string& II, const string& IS,
                                     double DCC, double eta, const EAngles& EA)
  {
  if(!SpinCheck(II,IS,1)) IDfatal(2); 		// Insure spin types valid
  Isotope III(II);				// Make isotope for II
  Isotope IIS(IS);				// Make isotope for IS
  if(!SpinCheck(II,IS,false,true))		// Disallow e-/nucleus pair
    { IDerror(60, 1); IDfatal(2); }		// if so its a fatal error
  _DCC = CheckDCC(II, IS, DCC);			// Set proper DCC
  double Iz = III.qn();                  	// Get Iz value of II
  double Sz = IIS.qn();          	        // Get Iz value of IS
  double X = xi();				// Get interaction constant
  IntRank2::operator=(IntRank2(Iz,Sz,X,eta,EA));// Use generic interaction
  setT20wh();					// Create added T2 component
  }

IntDip::IntDip(const Isotope& II, const Isotope& IS,
                                     double DCC, double eta, const EAngles& EA)
  { 
  if(!SpinCheck(II,IS,false,true))		// Disallow e-/nucleus pair
    { IDerror(60, 1); IDfatal(2); }		// if so its a fatal error
  _DCC = CheckDCC(II, IS, DCC);			// Set proper DCC
  double Iz = II.qn();                  	// Get Iz value of II
  double Sz = IS.qn();          	        // Get Iz value of IS
  double X = xi();				// Get interaction constant
  IntRank2::operator=(IntRank2(Iz,Sz,X,eta,EA));// Use generic interaction
  setT20wh(); 					// Create added T2 component
  }

IntDip::IntDip(double Iz, double Sz, double DCC, double eta, const EAngles& EA)
  {
  if(!SpinCheck(Iz,Sz,true)) { IDfatal(2); }	// Insure Iz & Sz valid
  _DCC = DCC;					// Store our dipolar coupling
  double X = xi();				// Get interaction constant
  IntRank2::operator=(IntRank2(Iz,Sz,X,eta,EA));// Use generic interaction
  setT20wh(); 					// Create added T2 component
  }

// ----------------------------------------------------------------------------
//                  Direct Constructors Using Spin Coordinates
// ----------------------------------------------------------------------------

/* Here we need to know the quantum numbers of the two spins involved and spin
   coordinates from which we can determine the dipolar coupling. The 
   interaction asymmetry will automatically be set to 0. 

           Input                D       : Dipolar interaction (this)
                                II	: Spin I isotope type
                                IS	: Spin S isotope type
                                pt1     : Coordinate of spin I (Angstroms)
                                pt2     : Coordinate of spin S (Angstroms)
           Output               none    : Dipolar interaction constructed
           Note                         : A fatal error will result if
                                          either II & IS are a mix of an
                                          electron and a nucleon
           Note                         : A non-fatal error will occur if
                                          the sign of DCC doesn't match the
                                          sign of gamma(II)*gamma(IS)
           Note                         : For dipolar interactions DCC=delzz */

IntDip::IntDip(const string& II, const string& IS, 
                                            const coord& pt1, const coord& pt2)
  {
  if(!SpinCheck(II,IS,1)) IDfatal(2); 		// Insure spin types valid
  Isotope III(II);				// Make isotope for II
  Isotope IIS(IS);				// Make isotope for IS
  if(!SpinCheck(III,IIS,0,1))			// Insure not e-/nucleon pair
    { IDerror(60, 1); IDfatal(2); }		// Disallow e-/nucleus pairing
  _DCC = CheckDCC(III, IIS, pt1, pt2);		// Set proper DCC
  double Iz = III.qn();				// Get Iz value of II
  double Sz = IIS.qn();				// Get Iz value of IS
  double X = xi();				// Get interaction constant
  double E = 0.0;				// No asymmetry
  EAngles EA(pt1.phi(pt2), pt1.theta(pt2), 0.0);
  IntRank2::operator=(IntRank2(Iz,Sz,X,E,EA));	// Use IntRank2 constructor
  setT20wh();					// Create added T2 component
  }

IntDip::IntDip(const Isotope& II, const Isotope& IS, 
                                            const coord& pt1, const coord& pt2)
  {
  if(!SpinCheck(II,IS,0,1))			// Insure not e-/nucleon pair
    { IDerror(60, 1); IDfatal(2); }		// Disallow e-/nucleus pairing
  _DCC = CheckDCC(II, IS, pt1, pt2);		// Set proper DCC
  double Iz = II.qn();				// Get Iz value of II
  double Sz = IS.qn();				// Get Iz value of IS
  double X = xi();				// Get interaction constant
  double E = 0.0;				// No asymmetry
  EAngles EA = EAzero;
  IntRank2::operator=(IntRank2(Iz,Sz,X,E,EA));	// Use IntRank2 constructor
  setT20wh();					// Create added T2 component
  }

// ----------------------------------------------------------------------------
//            Direct Constructors That Use Cartesian Parameters
// ----------------------------------------------------------------------------

/* Here we try to use the set { Iqn, Sqn, Dxx, Dyy, Dzz, DEAngles }. The spin
   quantum numbers sets the dimension of the interaction spin Hilbert space.
   The values of Dxx, Dyy, and Dzz are in kHZ and should sum to zero. They 
   are used to determine the dipolar coupling DCC (kHz), and (if any) the
   asymmetry ETA ([0,1]). The interaction orientation is taken from Euler 
   angles DEAngles. We will allow for default settings of the orientation and 
   asymmetry so that those arguments need not be input during construction.
   Their default values are 0.                                               */

IntDip::IntDip(const string& II, const string& IS,
                                        const coord& DxDyDz, const EAngles& EA)
  {
  if(!SpinCheck(II,IS,1)) IDfatal(2); 		// Insure spin types valid
  Isotope III(II);				// Make isotope for II
  Isotope IIS(IS);				// Make isotope for IS
  if(!SpinCheck(II,IS,0,1))			// Disallow e-/nucleus pair
    { IDerror(60, 1); IDfatal(2); }		// if so its a fatal error
  coord DiDzDe = AisoDelzEta(DxDyDz);           // Switch to spherical values
  _DCC = CheckDCC(II, IS, DiDzDe.y());		// Set proper DCC
  double eta = DiDzDe.z();			// Get the asymmetry
  double Iz  = III.qn();                  	// Get Iz value of II
  double Sz  = IIS.qn();          	        // Get Iz value of IS
  double X   = xi();				// Get interaction constant
  IntRank2::operator=(IntRank2(Iz,Sz,X,eta,EA));// Use generic interaction
  setT20wh();					// Create added T2 component
  }

IntDip::IntDip(const Isotope& II, const Isotope& IS,
                                        const coord& DxDyDz, const EAngles& EA)
  { 
  if(!SpinCheck(II,IS,0,1))			// Disallow e-/nucleus pair
    { IDerror(60, 1); IDfatal(2); }		// if so its a fatal error
  coord DiDzDe = AisoDelzEta(DxDyDz);           // Switch to spherical values
  _DCC = CheckDCC(II, IS, DiDzDe.y());		// Set proper DCC
  double eta = DiDzDe.z();			// Get the asymmetry
  double Iz = II.qn();                  	// Get Iz value of II
  double Sz = IS.qn();          	        // Get Iz value of IS
  double X = xi();				// Get interaction constant
  IntRank2::operator=(IntRank2(Iz,Sz,X,eta,EA));// Use generic interaction
  setT20wh(); 					// Create added T2 component
  }

IntDip::IntDip(double Iz, double Sz, const coord& DxDyDz, const EAngles& EA)
  {
  if(!SpinCheck(Iz,Sz,1)) { IDfatal(2); } 	// Insure Iz & Sz valid
  coord DiDzDe = AisoDelzEta(DxDyDz);           // Switch to spherical values
  _DCC = DiDzDe.y();				// Set proper DCC
  double eta = DiDzDe.z();			// Get the asymmetry
  double X = xi();				// Get interaction constant
  IntRank2::operator=(IntRank2(Iz,Sz,X,eta,EA));// Use generic interaction
  setT20wh(); 					// Create added T2 component
  }

// ----------------------------------------------------------------------------
//                     Constructors Using Parameter Sets
// ----------------------------------------------------------------------------

/* These constructors attempt to set the dipolar interaction from parameters
   in a given Parameter set. Either a single interaction index may be supplied
   or two spin indices. An additional variation (in support of spin systems)
   takes two spin isotopes. Note that the read functions pretty much provide
   the same abilities as these functions.

   When using an single interaction index the dipolar interaction will be 
   ignorant with regards to any particular spins that are involved, knowing
   only their I & S values. Thus we are forced to use Iqn & Sqn for the 
   parameters that set the spin quantum numbers. Obviously, parameters Iso &
   Coord cannot be allowed with interaction indexing. Hence, interaction index
   parameter sets are

     { Iqn(i), Sqn(i), DCC(i), Deta(i), DEAngles(i) }
     { Iqn(i), Sqn(i), DCC(i), Deta(i), Dalpha(i), Dbeta(i), Dgamma(i) }
     { Iqn(i), Sqn(i), DCC(i), Deta(i), Dxx(i), Dyy(i), Dzz(i), DEAngles(i) }
     { Iqn(i), Sqn(i), DCC(i), Deta(i), Dalpha(i), Dbeta(i), Dgamma(i) }

   If the parameters Iqn and Sqn are not present in the file/parameter set,
   DI=DS=1/2 is used by default. If Deta is not present (& this is normally 
   the case) then Deta is assumed to be zero. Similarly, the default Euler
   angles are zero which will set the interaction in its own PAS.
 
   When using two spin indices the dipolar interaction can find spin isotope
   types in the paramter set. The spin quantum numbers may then be set from
   either Iso(i) and Iso(j) or from Dqn(i) and Dqn(j). If the parameters for
   setting the spin quantum numbers are not present then DI=DS=1/2 is used by
   default. With the spin indexing scheme we may include the use of spin 
   coordinates to set the strength & orientation.  Possible parameter sets are

          { Iso(i), Iso(j), Coord(i), Coord(j) }
          { Iso(i), Iso(j), DCC(i,j), Deta(i,j), DEAngles(i,j) }
          { Iso(i), Iso(j), DxDyDz(i,j), DEAngles(i,j) }
          { Dqn(i), Dqn(j), DCC(i,j), Deta(i,j), DEAngles(i,j) }
          { Dqn(i), Dqn(j), DxDyDz(i,j), DEAngles(i,j) }
 
   Individual Cartesian components Duv(i,j) may also be used in lieu of the
   dipolar coupling, asymmetry, & orientation (if off-diagonals are specified).

   When spin isotopes and coordinates are used the dipolar coupling 
   constant (DCC) is given by

                                [                     2 ]  / 3
                  DELZZ = DCC = | gamma * gamma * hbar  | / r
                                [      I       S        ]/

           Input                D       : Dipolar interaction
                                pset    : Parameter set
                                idxI    : 1st spin or interaction index
                                idxS    : 2nd spin interaction index
                                warn    : Warning output label
                                           0 = no warnings
                                           1 = warnings
                                          >1 = fatal warnings
           Output               none    : Dipolar interaction constructed
                                          from parameters in pset            */

IntDip::IntDip(const ParameterSet& pset, int idxI, int idxS, int warn)
  {
  if(!setDI(pset, idxI, idxS, warn?1:0))
    {
    if(warn)
      {
      IDerror(2, 1);			// Can't do construction
      if(warn > 1) IDfatal(13);		// Fatal error
      else         IDerror(13);		// or a warning issued
      }
    }					// Note that setDI will also
  }					// create T20wh as needed

/*
IntDip::IntDip(const ParameterSet& pset, double Iz, double Sz,
                                               int idxI, int idxS, int warn)
  {
  if(!SpinCheck(Iz,warn?true:false)		// Insure reasonable Iz & Sz
  || !SpinCheck(Sz,warn?true:false)) IDfatal(11);
  bool TF;
  TF =  getDCC(pset, _DCC, idx, -1, warn);	// Get DCC value (in Hz)
  double X = xi();				// Get interacion constant
  string pname("Deta");				// Parameter name for asym.
  double ET;
  getAeta(pset, pname, ET, idx, warn);		// Get asymmetry value (Deta)
  IntRank2::operator=(IntRank2(Iz,Sz,X,ET));	// Use IntRank2 constructor
  setT20wh();					// Create added T2 component
  }

IntDip::IntDip(int idxI, int idxS, ParameterSet& pset, int warn)
  {
  if(!setDI(pset, idxI, idxS, warn?1:0))
    {
    if(warn)
      {
      IDerror(2, 1);			// Can't do construction
      if(warn > 1) IDfatal(13);		// Fatal error
      else         IDerror(13);		// or a warning issued
      }
    }
  }
*/


// ----------------------------------------------------------------------------
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------

void IntDip::operator= (const IntDip &D1) 
  {
  IntRank2::operator=(D1);		// Copy the rank 2 tensor
  T20wh = D1.T20wh; 			// Also copy weak heteronuclear T20
  _DCC  = D1._DCC;			// Copy dipolar coupling constant
  }

IntDip::~IntDip () { }


// ____________________________________________________________________________
// B                     SPATIAL TENSOR COMPONENT ACCESS
// ____________________________________________________________________________
 
//-----------------------------------------------------------------------------
//	             Anisotropy Access (Dipolar Coupling)
//-----------------------------------------------------------------------------
 
/* By definition, the dipolar coupling constant is identical to the spatial
   tensor delzz value in GAMMA and this is 2/3 the tensor anisotorpy. The 
   relationship is formally given by

                  mu
                 --- * hbar * gamma  * gamma               1/2
                 4pi               i        j     1   [ 5  ]      D   2       D
 delzz = DCC   = ____________________________ = - - * |----|  * Xi  = - DA = w
            ij                 3                  2   [6*pi]          3
                           r
                            ij

   The functions here allow users to get and/or set the dipolar coupling
   constant as well as calculate it based on the formula above.
 
           Input                D	: Dipolar interaction
	  			IsoI	: Spin I isotope type (1H, 2H, ...)
	  			IsoS	: Spin S isotope type (1H, 2H, ...)
           Return               Rij	: Dipolar distance, in meters
	   Note				: The sign of DCC depends on the
	  				  product of gyromagnetic ratios
	  				  of II && SS
           Note                         : The value of delzz is equivalent to
                                          the dipolar coupling constant DCC
           Note                         : Some functions are static & do not
                                          need an instance of IntDip to work
                                                    -1      -2
           Note				: 1T = 1 J-C  -sec-m                 */


double IntDip::CheckDCC(const Isotope& II, const Isotope& IS, double dcc)
  {
  double pm = II.gamma()*IS.gamma();	// Form gammaI x gammaS
  if(dcc<0 && pm>0)			// We don't allow negative DCC values
    {					// for this isotope combination
    IDerror(5, 1);			//	Negative DCC forbidden
    IDerror(6);				//	Using absolute value
    return fabs(dcc);
    }
  else if(dcc>0 && pm<0)		// We don't allow positive DCC values
    {					// for this isotope combination
    IDerror(7, 1);			//	Positive DCC forbidden
    IDerror(8);				//	Using negated value
    return -dcc;
    }
  return dcc;
  }

double IntDip::CheckDCC(const Isotope& II, const Isotope& IS,
                                            const coord& pt1, const coord& pt2)
  {
  double R = pt1.Rad(pt2);
  return DCC(II,IS,R);			// Return the dipolar coupling
  }

double IntDip::DCC()          const { return _DCC; }	// In Hz
void   IntDip::DCC(double dz)       { _DCC = dz;   }	// In Hz

double IntDip::DCC(const string& IsoI, const string& IsoS, double Rij)
  { return DCC(Isotope(IsoI), Isotope(IsoS), Rij); }

double IntDip::DCC(const Isotope& IsoI, const Isotope& IsoS, double Rij)
  {
  double gI = IsoI.gamma();			// Get gamma of I (1/T-sec)
  double gS = IsoS.gamma();			// Get gamma of S (1/T-sec)
  double Rcubed = Rij*Rij*Rij;			// Get r^3 in meters
  return MU*HBAR*gI*gS/(Rcubed*PIx2); 		// This is DCC in Hz
  }

double IntDip::DCC(const string& II, const string& SS,
                                          const coord& ptII, const coord& ptSS)
  { return DCC(Isotope(II), Isotope(SS), ptII.Rad(ptSS)); }

double IntDip::DCC(const Isotope& II, const Isotope&     SS,
                                          const coord& ptII, const coord& ptSS)
  { return DCC(II, SS, ptII.Rad(ptSS)); }

double IntDip::W2DCC(const string& I, const string& S, double W)
  { return W2DCC(Isotope(I), Isotope(S), W); }

double IntDip::W2DCC(const Isotope& IsoI, const Isotope& IsoS, double W)
  {
  if(IsoI == IsoS) return W/3.0;		// For homonuclear, this
  else             return W/2.0;		// For heteronuclear, this
  }

double IntDip::DCC2W(const string& I, const string& S, double D)
  { return DCC2W(Isotope(I), Isotope(S), D); }

double IntDip::DCC2W(const Isotope& IsoI, const Isotope& IsoS, double D)
  {
  if(IsoI == IsoS) return D*3.0;		// For homonuclear, this
  else             return D*2.0;		// For heteronuclear, this
  }

/*  Note that the default base class IntRank2A also provides the GAMMA
    normalized anisotropy and delzz values. These are constants given by

                       1/2
                  [ 5 ]                  ^      3
          del   = |---]  = 0.51503      /_\ A = - del   = 0.77255
             zz   [6PI]                         2    zz

    Addiontal functions that are standard for all rank 2 interactions are
    found below.                                                            */

// static double IntRank2A::delzz();                            INHERITED
// static double IntRank2A::delA();                             INHERITED

double IntDip::aniso() const       { return 1.5*_DCC; }
void   IntDip::aniso(double aiso)  { _DCC = aiso/1.5; }
double IntDip::DA()    const       { return 1.5*_DCC; }
void   IntDip::DA(double aiso)     { _DCC = aiso/1.5; }

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

// These allow one to access the irreducible spherical elements at a specific   
// orientation without rotating the entire tensor.  In these functions theta is
// the angle down from the PAS z axis and phi the angle over from PAS x axis.
// Note that these all use GAMMA scaling which sets A2m to be rank 2 spherical
// harmonics at all oreintations if there is no asymmetry.
 
/*
  complex IntRank2A::A0() const;				INHERITED
  complex IntRank2A::A0(double theta, double phi) const;	INHERITED
  complex IntRank2A::A20( ) const;				INHERITED
  complex IntRank2A::A20(double theta, double phi) const;	INHERITED
  complex IntRank2A::A1() const;				INHERITED
  complex IntRank2A::A1(double theta, double phi) const;	INHERITED
  complex IntRank2A::A21( ) const;				INHERITED
  complex IntRank2A::A21(double theta, double phi) const;	INHERITED
  complex IntRank2A::Am1() const;				INHERITED
  complex IntRank2A::Am1(double theta,double phi) const;	INHERITED
  complex IntRank2A::A2m1( ) const;				INHERITED
  complex IntRank2A::A2m1(double theta, double phi) const;	INHERITED
  complex IntRank2A::A2() const;				INHERITED
  complex IntRank2A::A2(double theta, double phi) const;	INHERITED
  complex IntRank2A::A22( ) const;				INHERITED
  complex IntRank2A::A22(double theta, double phi) const;	INHERITED
  complex IntRank2A::Am2() const;				INHERITED
  complex IntRank2A::Am2(double theta,double phi) const;	INHERITED
  complex IntRank2A::A2m2( ) const;				INHERITED
  complex IntRank2A::A2m2(double theta, double phi) const;	INHERITED
  complex IntRank2A::Acomp(int comp) const;			INHERITED    */
  
//-----------------------------------------------------------------------------
//	       Cartesian Tensor Component Access - Normalized
//-----------------------------------------------------------------------------

// These allow one to access the irreducible Cartesian elements at a specific   
// orientation without rotating the entire tensor.  In these functions theta is
// the angle down from the PAS z axis and phi the angle over from PAS x axis.
// Note that these all use GAMMA scaling which sets A2m to be rank 2 spherical
// harmonics at all oreintations if there is no asymmetry.

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
//	      Cartesian Tensor Component Access - Not Normalized
//-----------------------------------------------------------------------------

double IntDip::dxx() const { return _DCC*RT6PIO5*Axx(); }
double IntDip::dyy() const { return _DCC*RT6PIO5*Ayy(); }
double IntDip::dzz() const { return _DCC*RT6PIO5*Azz(); }
double IntDip::dxy() const { return _DCC*RT6PIO5*Axy(); }
double IntDip::dyx() const { return _DCC*RT6PIO5*Ayx(); }
double IntDip::dxz() const { return _DCC*RT6PIO5*Axz(); }
double IntDip::dzx() const { return _DCC*RT6PIO5*Azx(); }
double IntDip::dyz() const { return _DCC*RT6PIO5*Ayz(); }
double IntDip::dzy() const { return _DCC*RT6PIO5*Azy(); }
  
double IntDip::dxx(double T, double P) const {return _DCC*RT6PIO5*Axx(P,T,0);}
double IntDip::dyy(double T, double P) const {return _DCC*RT6PIO5*Ayy(P,T,0);}
double IntDip::dzz(double T, double P) const {return _DCC*RT6PIO5*Azz(P,T,0);}
double IntDip::dyx(double T, double P) const {return _DCC*RT6PIO5*Ayx(P,T,0);}
double IntDip::dxy(double T, double P) const {return _DCC*RT6PIO5*Axy(P,T,0);}
double IntDip::dzx(double T, double P) const {return _DCC*RT6PIO5*Azx(P,T,0);}
double IntDip::dzy(double T, double P) const {return _DCC*RT6PIO5*Azy(P,T,0);}
double IntDip::dxz(double T, double P) const {return _DCC*RT6PIO5*Axz(P,T,0);}
double IntDip::dyz(double T, double P) const {return _DCC*RT6PIO5*Ayz(P,T,0);}

double IntDip::dxx(double A, double B, double G) const 
                                           { return _DCC*RT6PIO5*Axx(A,B,G); }
double IntDip::dyy(double A, double B, double G) const
                                           { return _DCC*RT6PIO5*Ayy(A,B,G); }
double IntDip::dzz(double A, double B, double G) const
                                           { return _DCC*RT6PIO5*Azz(A,B,G); }
double IntDip::dyx(double A, double B, double G) const
                                           { return _DCC*RT6PIO5*Ayx(A,B,G); }
double IntDip::dxy(double A, double B, double G) const
                                           { return _DCC*RT6PIO5*Axy(A,B,G); }
double IntDip::dzx(double A, double B, double G) const
                                           { return _DCC*RT6PIO5*Azx(A,B,G); }
double IntDip::dzy(double A, double B, double G) const
                                           { return _DCC*RT6PIO5*Azy(A,B,G); }
double IntDip::dxz(double A, double B, double G) const
                                           { return _DCC*RT6PIO5*Axz(A,B,G); }
double IntDip::dyz(double A, double B, double G) const
                                           { return _DCC*RT6PIO5*Ayz(A,B,G); }
  
double IntDip::dxx(const EAngles& EA) const { return _DCC*RT6PIO5*Axx(EA); }
double IntDip::dyy(const EAngles& EA) const { return _DCC*RT6PIO5*Ayy(EA); }
double IntDip::dzz(const EAngles& EA) const { return _DCC*RT6PIO5*Azz(EA); }
double IntDip::dyx(const EAngles& EA) const { return _DCC*RT6PIO5*Ayx(EA); }
double IntDip::dxy(const EAngles& EA) const { return _DCC*RT6PIO5*Axy(EA); }
double IntDip::dzx(const EAngles& EA) const { return _DCC*RT6PIO5*Azx(EA); }
double IntDip::dzy(const EAngles& EA) const { return _DCC*RT6PIO5*Azy(EA); }
double IntDip::dxz(const EAngles& EA) const { return _DCC*RT6PIO5*Axz(EA); }
double IntDip::dyz(const EAngles& EA) const { return _DCC*RT6PIO5*Ayz(EA); }

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
//        matrix IntRank2::CartMx(double scale=1.0) const
//        bool   IntRank2::Spherical( )
//        bool   IntRank2::Isotropic( )
//        matrix IntRank2::CartMx(bool scf)                                 */
 

matrix IntDip::Dmx() const
  {
  matrix DMX(3,3,h_matrix_type);
  DMX.put(  dxx(), 0, 0);
  DMX.put_h(dxy(), 0, 1);
  DMX.put_h(dxz(), 0, 2);
  DMX.put(  dyy(), 1, 1);
  DMX.put_h(dyz(), 1, 2);
  DMX.put(  dzz(), 2, 2);
  return DMX;
  }

// ____________________________________________________________________________
// C                       SPIN TENSOR COMPONENT ACCESS
// ____________________________________________________________________________
 
//-----------------------------------------------------------------------------
//                         Spin Quantum Value Access
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

/*     The commented member functions are INHERITED from class IntRank2T     */

// ----------------------------------------------------------------------------
//        Single Spin Or Spin Pair Hilbert Space Spin Tensor Components
// ----------------------------------------------------------------------------

//matrix IntRank2T::T20()           const;                // Return T2,0
//matrix IntRank2T::T21()           const;                // Return T2,1
//matrix IntRank2T::T2m1()          const;                // Return T2,-1
//matrix IntRank2T::T22()           const;                // Return T2,2
//matrix IntRank2T::T2m2()          const;                // Return T2,-2
//matrix IntRank2T::T2m(int m)      const;                // Return T2,m
//matrix IntRank2T::Tcomp(int comp) const;                // Return T2,m

// ----------------------------------------------------------------------------
//               Composite Hilbert Space Spin Tensor Components
// ----------------------------------------------------------------------------

/* Here HSs are single spin Hilbert space dimensions and the spins involved
   are indexed by i and j.                                                   */

//matrix IntRank2T::T20(       const vector<int>& HSs, int i, int j) const;
//matrix IntRank2T::T21(       const vector<int>& HSs, int i, int j) const;
//matrix IntRank2T::T2m1(      const vector<int>& HSs, int i, int j) const;
//matrix IntRank2T::T22(       const vector<int>& HSs, int i, int j) const;
//matrix IntRank2T::T2m2(      const vector<int>& HSs, int i, int j) const;
//matrix IntRank2T::T2m(int m, const vector<int>& HSs, int i, int j) const;

// ----------------------------------------------------------------------------
//      Additional Function To Handle "Perturbing" Dipolar Interactions
// ----------------------------------------------------------------------------

/* These are the components of T20 that are invariant under rotations about the
   z-axis in the case that the spin pair is heteronuclear. They facilitate
   construction of dipolar Hamiltonians in the rotationg frame when a dipolar
   interaction very weak relative to the Zeeman interaction. See T20wh and
   the function H0.                                                          */ 

matrix IntDip::T20het() const { return T20wh; }
matrix IntDip::T20het(const vector<int>& HSs, int i, int j) const
                             { return blow_up(T20wh, HSs, i, j); }

// ____________________________________________________________________________
// D          DIPOLE DIPOLE INTERACTION CONSTANT ACCESS FUNCTIONS
// ____________________________________________________________________________

/* The dipolar interaction constant is what sets the overall interaction
   strength and allows GAMMA to track different interaction types using
   "normalized" spatial and spin tensor (classes IntRank2A and IntRank2T). For
   dipolar interactions the interaction constant is defined as (in radians/sec)

                     1/2
               [6*pi]     mu
          -2 * |----|   * --- * hbar * gamma  * gamma               1/2
     D         [ 5  ]     4pi               i        j        [6*pi]
   xi   = ____________________________________________ = -2 * |----|  * DCC
     ij                       3                               [ 5  ]
                             r
                              ij

   Since the dipolar coupling constant is maintained in Hz internally whereas
   the interaction constant is kept in rad/sec we need a 2*Pi scaling.
 
	   Input		D	: Dipolar interaction
           Return               xi      : Dipolar. interaction constant
           Note                           Value is in radians/sec
          				            -1      -2
           Note				: 1T = 1 J-C  -sec-m
                                                                             */

double IntDip::xi(bool Hz) const { return Hz?-RT6PIO5*_DCC/PI:-2.0*RT6PIO5*_DCC; }

// ____________________________________________________________________________
// E             DIPOLE DIPOLE INTERACTION AUXILIARY FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                         Internuclear Distances
// ----------------------------------------------------------------------------

/* We cannot know the internuclear distance between the spins in our dipolar
   interaction because we don't know what kind of spins are involved.  We can
   however calculate the distance if we use the interaction dipolar coupling
   along with user specification of spin types (to get gyromagnetic ratios)  */ 

 
	// Input		D	: Dipolar interaction (this)
	//			II      : Spin I isotope type (1H, 2H, ...)
	//			SS      : Spin S isotope type (1H, 2H, ...)
	//			gchk	: Check if Iz & Sz Match II & SS 
	//					0 = don't bother
	//					1 = make sure
        // Return               R	: Dipolar distance, Angstroms
        //                                          -1      -2
        // Note				: 1T = 1 J-C  -sec-m
	// Note				: If we know the isotopes of the
	//				  two spins (input) and the dipolar
	//				  coupling (in class) we can determine
	//				  the internuclear distance.  This DOES 
	//				  check that II && SS have spin Iz
	//				  values matching those of D! 


/*                     [ mu                           ]1/3
                       | --- * hbar * gamma  * gamma  |
                       | 4pi               i        j |        [  NUM  ]
                 r   = | ____________________________ |     =  | ----- |
                  ij   |                              |        [ DENOM ]
                       |           DCC                |
                       [              ij              ]
                                                                  2 -1 -2
  The gyromagnetic ratio's herein will be ~1.e8/T-sec, or 1.e8 C m J  s.
  The product of hbar (~1.e-34 J-sec) and mu/4pi (~1.e-7 J-sec-C-m) with
  the gamma's from each spin will produce a value which is roughly

           -7              -34          16 -2 -2 4 -4    -25 3 -2 5 2
   NUM ~ 10  J-sec-C-m * 10   J-sec * 10  J  C  m s  = 10   m s  C m
  
  whereas the denominator will be ~ 1.e5/sec.  So, the returned values
  from this function can be estimated as around

                 [   -30 3 ]1/3     -10
           r   ~ | 10   m  |    = 10   m = 1 A
            ij   [         ]

  Thus the numerator above is roughly ~1.e-25                                */


double IntDip::R(const string& II, const string& SS, bool gchk) const
  {
  if(!_DCC) return _DCC;			// Quit if no interaction
  Isotope IsoI(II);				// Get I isotope
  Isotope IsoS(SS);				// Get S isotope
  if(gchk)					// If desired, make sure that
    {						// both Iz & Sz for II & SS
    if(IsoI.qn() != Izval()			// match interaction values
                         || IsoS.qn()!=Szval())
      IDerror(50);
    }
  double gI = IsoI.gamma();			// Get gamma of I (1/T-sec)
  double gS = IsoS.gamma();			// Get gamma of S (1/T-sec)
  double brack = MU*HBAR*gI*gS/(_DCC*PIx2); 	// For _DCC in Hz
  return pow(brack, 1.0/3.0)*1.e10;		// Return in Angstroms
  }

// ----------------------------------------------------------------------------
//                      Transition Splitting Frequencies
// ----------------------------------------------------------------------------

/* The dipolar frequency is the observed splitting between transitions caused
   from dipolar coupling of a spin pair.  For the spectrum of spin I coupled
   to spin S, there are 2*S+1 transitions (from the Dipolar Hamiltonian).  The
   values in these function are determined under the condition that the dipolar
   coupling is weak relative to the Zeeman interaction (high field approx.,
   first order terms only).  The splitting is orientationally dependent, the
   default value is when the tensor is oriented with it's principal axes (PAS)
   coinciding with the static field of the spectrometer.  If the tensor isn't
   PAS-field aligned, the splitting will vary with orientation according to

                       1            2                2                        
      W (theta,phi)  = - W   [ 3*cos (theta) - eta*sin (theta)cos(2*phi) ]
       D               2  D,o

   where theta is the angle down from the z axis and phi the angle over from
   the x axis of the internuclear vector.  Alternatively, if the Zeeman terms
   aint much stronger the splittings won't be equally spaced (2nd order terms).

   We CANNOT KNOW the splitting frequency between the spins in our dipolar
   interaction because we don't know what kind of spins are involved.  We can
   however calculate the frequency if we use the interaction dipolar coupling
   along with user specification of spin types (to get gyromagnetic ratios)

                Dipolar Splittings From Coupling to I=1/2 Nucleus

            HETERONUCLEAR                           HOMONUCLEAR
     ->    ->           ->  ->               ->    ->           ->   ->
     B  || r           B  | r                B  || r            B  | r
      o     ij          o -  ij               o     ij           o -  ij
 
   |         |       |         |           |         |         |  3     |
   |<-2*DCC->|       |<- DCC ->|           |<-3*DCC->|         |<--DCC->|
   |         |       |         |           |         |         |  2     |
 __|_________|__   __|_________|___     ___|_________|___  ____|________|___
 -DCC       DCC   -DCC/2     DCC/2      -3DCC/2    3DCC/2  -3DCC/4   3DCC/4


   Thus the splitting frequency will typically be some multiple of DCC. For
   the homonuclear case oriented along the magnetic field +z it will be 3*DCC
   and in the heteronuclear case it is 2*DCC.  The values change at different
   orientations.  

   Again, keep in mind that this class has NO knowledge of HOMO- vs. HETERO- 
   nuclear spin pairs.  Thus it cannot distinguish between the two.  That must
   be done in higher level classes or via user input (like these functions).  */
 
        // Input                D	: Dipolar interaction
	//			II      : Spin I isotope type (1H, 2H, ...)
	//			SS      : Spin S isotope type (1H, 2H, ...)
	//			gchk	: Check if Iz & Sz Match II & SS 
	//					0 = don't bother
	//					1 = make sure
        // Return               wD      : Dipolar frequency (Hz)
	// Note				: If we know the isotopes of the
	//				  two spins (input) and the dipolar
	//				  coupling (in class) we can determine
	//				  the high field splitting frequency.

double IntDip::W(const string& II, const string& SS, bool gchk) const
  {
  if(gchk)					// If desired, make sure that
    {						// both Iz & Sz for II & SS
    Isotope IsoI(II);				// Get I isotope
    Isotope IsoS(SS);				// Get S isotope
    if(IsoI.qn() != Izval()			// match interaction values
                         || IsoS.qn()!=Szval())
      IDerror(51);
    }
  if(II == SS) return 3.0*_DCC;			// For homonuclear, this
  else         return 2.0*_DCC;			// For heteronuclear, this
  }						// In Hz, as is _DCC stored

void IntDip::W(double w, const string& II, const string& SS, int gchk)
  {
  if(gchk)					// If desired, make sure that
    {						// both Iz & Sz for II & SS
    Isotope IsoI(II);				// Get I isotope
    Isotope IsoS(SS);				// Get S isotope
    if(IsoI.qn() != Izval()			// match interaction values
                         || IsoS.qn()!=Szval())
      IDerror(51);
    }
  if(II == SS) _DCC = w/3.0;			// Set _DCC homonuclear spins
  else         _DCC = w/2.0;			// Set _DCC heteronuclear spins
  }


// ----------------------------------------------------------------------------
//                 Spin & Interaction Index Checking Functions
// ----------------------------------------------------------------------------

bool IntDip::DCheck(const string& II, const string& IS, int warn) const
  {
  if(!Isotope::known(II))			// Make sure that string II
    {						// is a valid isotope in GAMMA
    if(warn)
      {
      ISTerror(110, II, 1);			//   Unknown isotope type
      if(warn > 1) ISTfatal(53);		//   Bad spin type
      else         ISTerror(53,1);
      }
    return false;
    }
  if(!Isotope::known(IS))			// Make sure that string IS
    {						// is a valid isotope in GAMMA
    if(warn)
      {
      ISTerror(110, IS, 1);			//   Unknown isotope type
      if(warn > 1) ISTfatal(53);		//   Bad spin type
      else         ISTerror(53,1);
      }
    return false;
    }
  return true;
  }

 
 
        // Input                D       : Dipolar interaction
        //                      Iqn,Sqn : Spin quantum values of I&S
        // Output               TF      : True if both I & S have spin 
        //                                quantum values that are non-zero 
        //                                positive multipes of 1/2 

int IntDip::Dspincheck(double Iqn, double Sqn) const
  {
  int TF = 1;                           // Assume everyone is OK
  int twoI = int(2.0*Iqn);              // We must insist that both I
  int twoS = int(2.0*Sqn);              // and S have spin angular momentum
  if(2.0*Iqn - double(twoI)) TF=0;      // which are non-zero integer
  else if(2.0*Sqn - double(twoS)) TF=0; // multiples of 1/2
  if(Iqn<=0.0 || Sqn<=0.0) TF=0;
  return TF;
  }

 
        // Input                D       : Dipolar interaction
	//			idxI       : Spin isotope type
	//			idxS       : Spin isotope type
        // Output               TF      : True if both I & S have spin 
        //                                indices that are appropriate

int IntDip::Dindexcheck(int idxI, int idxS) const
  {
  int TF=1;
  if(idxI == idxS) TF=0;
  else if(idxI<0 || idxS<0) TF=-1;
  return TF;
  }

// ____________________________________________________________________________
// F                        PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//      Functions To Make A Parameter Set From A Dipolar Interaction
// ----------------------------------------------------------------------------

// Note that the preferred means of specifying a dipolar interaction when there
// is no knowledge of spin isotope types nor spin coordinates is used for 
// filling the parameter set.  The base parameters of individual interactions 
// in that case are { DI(#), DS(#), DCC(#), DTheta(#), DPhi(#), DEta(#) }.
// If there is not eta value, as is typical for most dipolar interactions, the
// DEta parameter will not be written.

	// Input		D	: Dipolar interaction
	//  			pset	: Parameter set
        //                      idx     : Interaction index (default -1)
	//			pfx	: Interaction 2nd indx (def -1)
        // Output               void    : Interaction parameters are
        //                                are added ot the parameter set
        //                                with interaction index idx
	// Note				: The parameter names & types
	//				  output here MUST match those used
	//				  in setting dipolar interactions
	//				  from parameters sets
 
	// Output		pset	: Parameter set with the dipolar
	//			          interaction parameters added
	// Note				: Note that PSetAdd is used for
	//				  convenience here.


IntDip::operator ParameterSet( ) const
  { ParameterSet pset; pset += (*this);  return pset;}

void operator+= (ParameterSet& pset, const IntDip& D) { D.PSetAdd(pset); } 

void IntDip::PSetAdd(ParameterSet& pset, int idx, int pfx) const
  {
  string suffx;					// Parameter suffix
  if(idx != -1)					// Only use suffix if idx
    suffx = string("(")+Gdec(idx)+string(")");	// is NOT -1
  string prefx;					// Parameter prefix
  if(pfx != -1)					// Only use prefix if pfx
    prefx = string("[")+Gdec(pfx)+string("]");	// is NOT -1

  string pname  = string("IR2eta") + suffx;	// Add asymmetry parameter
  string pstate = string("Spin I Quantum Number");
  double pdatad = Izval();
  SinglePar par = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);

  pname  = prefx + string("DS") + suffx;	// Add DS spin quantum number
  pstate = string("Spin S Quantum Number");
  pdatad = Szval();
  par    = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);

  if(_DCC>999 || !_DCC)
    {
    pname  = prefx + string("DCC") + suffx;	// Add dip. coupling constant
    pstate = string("Dipolar Coupling Constant (kHz)");
    pdatad = _DCC*1.e-3;
    }
  else
    {
    pname  = prefx + string("DCCHz") + suffx;	// Add dip. coupling constant
    pstate = string("Dipolar Coupling Constant (Hz)");
    pdatad = _DCC;
    }
  par    = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);

  if(eta())
    {
    pname  = prefx + string("Deta") + suffx;	// Add asymmetry
    pstate = string("Dipolar Asymmetry)");
    pdatad = eta();
    par    = SinglePar(pname, pdatad, pstate);
    pset.push_back(par);
    }
  }

// ----------------------------------------------------------------------------
//    Functions To Output Dipolar Interaction To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

int IntDip::write(const string &filename, int idx, int pfx, int warn) const

	// Input		AD	: Dipolar interaction (this)
	//			filename: Output file name
	//			idx	: Interaction index (default -1)
	//			pfx	: Interaction 2nd indx (def -1)
	//			warn	: Warning level
	// Output		none 	: Dipolar interaction is
	//				  written as a parameter set to
	//				  file filename

  {
  ofstream ofstr(filename.c_str());	// Open filename for input
  if(!write(ofstr, idx, pfx, warn?1:0))	// If file bad then exit
    {
    IDerror(40, filename, 1);		// Filename problems
    if(warn>1) IDfatal(20);		// Fatal error
    return 0;
    }
  ofstr.close();			// Close it now
  return 1;
  }


int IntDip::write(ofstream& ofstr, int idx, int pfx, int warn) const

	// Input		AD	: Dipolar interaction (this)
	//			ofstr	: Output file stream
	//			idx	: Interaction index (default -1)
	//			pfx	: Interaction 2nd indx (def -1)
	//			warn	: Warning level
	// Output		none 	: Dipolar interaction is
	//				  written as a parameter set to
	//				  file filename

  {
  ParameterSet pset;			// Declare a parameter set
  PSetAdd(pset, idx, pfx);		// Add in interaction parameters
  if(!pset.write(ofstr, warn?1:0))	// Use parameter set to write
    {
    if(warn)
      {
      IDerror(22, 1);			// Problems writing to filestream
      if (warn>1) IDfatal(23);	// Fatal error
      }
    return 0;
    }
  return 1;
  }


// ____________________________________________________________________________
// H                   DIPOLAR INTERACTION INPUT FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//          Functions To Read ASCII File Using Two Spin Indicies
// ----------------------------------------------------------------------------

/* These next two read functions utilize two spin indices, for the spin pair
   involved in the dipolar interaction.  They'll try to read the parameter set

                   { Iso(i), Iso(j), Coord(i), Coord(j) }
   or
       { Iso(i), Iso(j), DCC(i,j), Dtheta(i,j), Dphi(i,j), Deta(i,j) }

   The functions do NOT allow for the suffix (#), since it is used by the spin
   index anyway. The functions do NOT allow for the prefix [#] either, because
   multiple dipole interactions can be defined in the same file by switching
   spin pair indices.  Multiple sets of interactions can be read using dipole
   interaction vectors (see class IntDipVec.)  The asymmetry, eta, is
   assumed to be zero (normally the case with dipolar interactions)          */

	// Input		D	: Dipolar interaction
	//			filename: Output file name
	//			idxI	: Index for spin I
	//			idxS	: Index for spin S
        //                      warn    : Warning output label
        //                                 0 = no warnings
        //                                 1 = warnings
        //                                >1 = fatal warnings
	// Output		TF      : Dipolar interaction is read in
	//				  from parameters in file filename
	//					TF = 0 couldn't read it
	//				        TF = 1 read with iso & coord
	//				        TF = 2 read with iso & DCC
	// Note				: Use of spin indices allow for
	//				  a variety of acceptable parameters
	// Input		D	: Dipolar interaction
	//			pset    : Parameter set
	//			IdxI	: Index for spin I
	//			IdxS	: Index for spin S
        //                      warn    : Warning output label
        //                                 0 = no warnings
        //                                 1 = warnings
        //                                >1 = fatal warnings
	// Output		TF 	: Dipolar interaction is read in
	//				  from parameters in pset
	//				  Return is true if interaction read
	//					TF = 0 couldn't read it
	//				        TF = 1 read with iso & coord
	//				        TF = 2 read with iso & DCC
	// Note				: Use of spin indices allow for
	//				  a variety of acceptable parameters


bool IntDip::read(const string &filename, int idxI, int idxS, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename, warn?1:0))	// Read in pset from file
    {					// If we cannot read the file at all
    if(warn)				// then issue warnings as desired
      {
      IDerror(40, filename, 1);		// 	Filename problems
      if(warn > 1) IDfatal(21);		// 	Fatal error
      else         IDerror(21);		//      or a warning issued
      }
    return false;			//  Return that we failed!
    }
  return read(pset, idxI, idxS, warn);	// User overloaded function
  }

bool IntDip::read(const ParameterSet& pset, int idxI, int idxS, int warn)
  {
  bool TF = setDI(pset, idxI, idxS, warn)?true:false;	// Use overload to read 
  if(!TF)						// If getDI didn't handle
    {							// setting interaction proper
    if(warn)						// then we'll issue some
      {							// warnings if desired
      string val;
      int eidx =44;
      if(idxI == -1)      { val = string(" None"); }
      else if(idxS == -1) { val = Gdec(idxI);      }
      else
        { val = Gdec(idxI) + string(" & ")
              + Gdec(idxS);
          eidx = 45;
        }
      if(warn > 1) IDfatal(eidx, val);		// Fatal error
      else         IDerror(eidx, val);		// or a warning issued
      }
    return false;				// Return that we failed
    }
  return TF;
  }

// ----------------------------------------------------------------------------
//       Functions To Read ASCII File Using One Interaction Index
// ----------------------------------------------------------------------------

/* These next two read functions will try an read the parameter set

                   { DI, DS, DCC, Dtheta, Dphi, Deta }

   Since there no spin indices given, parameters Iso and Coord are not
   allowed.  The functions do allow for the suffix (#) but it is not used
   in the default function call.  The functions do NOT allow for the prefix
   [#] since multiple interactions can be read using the suffix and multiple
   sets of interactions can be read using dipole interaction vectors (see
   class IntDipVec.)  If the parameters DI and DS are not present in the
   file/parameter set, DI=DS=1/2 is used by default.  If Deta is not
   present (& this is normally the case) then Deta is assumed to be zero.    */

//  { *(this) = IntDip(pset, idx); return 1; }	// Use the constructor
// sosi - this must be fixed to return TF properly




/* These next two read functions will try an read the parameter set 

                      { DCC, Dtheta, Dphi, Deta }

   Since there no spin indices given, parameters Iso and Coord are not 
   allowed.  Furthermore, the spin quantum numbers of the two spins
   involded in the interaction are taken to be 1/2, i.e. I=S=1/2. These
   functions do allow for the suffix (#) and the prefix [#], but they are
   not used in the default function call.  If Deta is not present
   (and this is normally the case) then Deta is assumed to be zero.	     */



	// Input		D	: Dipolar interaction
	//			filename: Output file name
	//			I       : Quantum I value (0.5, 1, 3.5, ...)
	//			S       : Quantum S value (0.5, 1, 3.5, ...)
	//			idx	: Interaction index (default -1->none)
	// Output		none 	: Dipolar interaction is read in
	//				  from parameters in file filename
	// Note				: Since no spin indices are used only
	//				  a few parameters are acceptable

void IntDip::read(const string &filename, double I, double S, int idx)
   {
   ParameterSet pset;                // Declare a parameter set
   if(!pset.read(filename, 1))		// Read in pset from file
     {
     IDerror(40, filename,1);		// Filename problems
     IDfatal(21);			// Fatal error
     }
   read(pset, I, S, idx);	 	// User overloaded function
   return;
   }



	// Input		D	: Dipolar interaction
	//			pset    : Parameter set
	//			Iqn	: Quantum I value (0.5, 1, ...)
	//			Sqn	: Quantum S value (0.5, 1, ...)
	//			idx	: Interaction index (default -1->none)
	// Output		none 	: Dipolar interaction is read in
	//				  from parameters in pset
	// Note				: Since no spin indices are used only
	//				  a few parameters are acceptable

void IntDip::read(const ParameterSet& pset, double Iqn, double Sqn, int idx)
  { 
  if(!Dspincheck(Iqn, Sqn))		// Insure reasonable Iqn & Sqn values
    IDfatal(11);
  *(this) = IntDip(Iqn, Sqn, pset, idx);// Use the constructor
  }


// ----------------------------------------------------------------------------
//        Interactive Ask/Read ASCII File Using A Single Interaction Index
// ----------------------------------------------------------------------------

        // Input                D	: Dipolar interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
	//			idx	: Interaction index
        // Output               String  : The parameter argn of array argc is
        //                                used to supply a filename from which
        //                                the dipolar interaction is read

string IntDip::ask_read(int argc, char* argv[], int argn, int idx)
  {
  string filename;				// Name of parameter file
  query_parameter(argc, argv, argn,		// Get filename from command
    "\n\tDipolar Interaction filename? ",	// Or ask for it
                                      filename);
  read(filename, idx); 				// Read system from filename
  return filename;
  }

string IntDip::ask_read(int argc, char* argv[], int argn, const string& def, int idx)
  {
  string msg = "\n\tDipolar Interaction filename ["     // Query we will ask if
             + def + "]? ";             	        // it is needed
  string filename = def;				// Name of spin system file
  ask_set(argc,argv,argn,msg,filename);			// Or ask for it
  read(filename, idx);					// Read system from filename
  return filename;					// Return filename
  }


// ----------------------------------------------------------------------------
//          Interactive Ask/Read ASCII File Using Two Spin Indicies
// ----------------------------------------------------------------------------


        // Input                D	: Dipolar interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
	//			idxI	: Index for spin I
	//			idxS	: Index for spin S
        // Output               String  : The parameter argn of array argc is
        //                                used to supply a filename from which
        //                                the dipolar interaction is read
        //                                If argument argn is not in argv, the
        //                                user is asked to supply a filename
        //                                The set filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing known IntDip parameters
        // Note                         : The interaction D is modifed (filled)

string IntDip::ask_read(int argc, char* argv[], int argn, int idxI, int idxS)
  {
  string filename;				// Name of parameter file
  query_parameter(argc, argv, argn,		// Get filename from command
    "\n\tDipolar Interaction filename? ",	// Or ask for it
                                      filename);
  read(filename, idxI, idxS); 			// Read system from filename
  return filename;
  }

        // Input                D	: Dipolar interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
	//			Iqn	: Quantum I value (0.5, 1, ...)
	//			Sqn	: Quantum S value (0.5, 1, ...)
	//			idx	: Interaction index (default -1->none)
        // Output               String  : The parameter argn of array argc is
        //                                used to supply a filename from which
        //                                the dipolar interaction is read
        //                                If argument argn is not in argv, the
        //                                user is asked to supply a filename
        //                                The set filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing known IntDip parameters
        // Note                         : The interaction D is modifed (filled)

string IntDip::ask_read(int argc, char* argv[], int argn,
                                              double Iqn, double Sqn, int idx)
  {
  string filename;				// Name of parameter file
  query_parameter(argc, argv, argn,		// Get filename from command
    "\n\tDipolar Interaction filename? ",	// Or ask for it
                                      filename);
  read(filename, Iqn, Sqn, idx); 		// Read system from filename
  return filename;
  }

// ----------------------------------------------------------------------------
//               Interactive Ask For All Kinds Interaction Info
// ----------------------------------------------------------------------------


        // Input                D	: Dipolar interaction (this)
	//			argc    : Number of arguments
	//			argv    : Array of arguments
	//			qn      : Query value
	//			DI	: Spin quantum number
	//			DS	: Spin quantum number
	//			Dcc	: Dipolar. coupling constant (Hz)
	//			Deta    : Dipolar. asymmetry
	//			Dtheta  : Dipolar. orientation angle
	//			Dphi	: Dipolar. orientation angle
        //                      Dflag   : Flag is eta is requested
        // Output               none    : The values of qn, I, Dcc, Deta,
	//				  Dtheta, and Dphi are set herein
	// Note				: This is INTERACTIVE!
	// Note				: This does NOT construct D!

void IntDip::ask(int argc, char* argv[], int& qn, double& DI, double& DS,
                                          double& Dcc, double& Deta, int Dflag)
  {
  askI(argc,   argv, qn++, DI);				// Ask for I type
  askS(argc,   argv, qn++, DS);				// Ask for S type
  askDCC(argc, argv, qn++, Dcc);			// Ask for DCC
//  query_parameter(argc, argv, qn++,			// Read in theta angle
//  "\n\tOrientation Down From PAS z [0, 180]? ", Dtheta);
//  query_parameter(argc, argv, qn++,			// Read in phi angle
//   "\n\tOrientation Over From PAS x [0, 360]? ", Dphi);
  if(Dflag)
    query_parameter(argc, argv, qn++,			// Read in the coupling
               "\n\tDipolar Asymmetry [0, 1]? ", Deta);
  else Deta = 0;
  }

void IntDip::askI(int argc, char* argv[], int qn, double& DI)
  {
  query_parameter(argc, argv, qn,			// Read in the I value
    "\n\tI Spin Quantum Value (0.5, 1, 1.5, ..)? ", DI);
  }

void IntDip::askS(int argc, char* argv[], int qn, double& DS)
  {
  query_parameter(argc, argv, qn,			// Read in the I value
    "\n\tS Spin Quantum Value (0.5, 1, 1.5, ..)? ", DS);
  }

void IntDip::askDCC(int argc, char* argv[], int qn, double& Dcc)
  {
  query_parameter(argc, argv, qn,			// Ask for frequency
                    "\n\tDipolar Coupling(kHz)? ", Dcc);// if not in argv
  Dcc *= 1.e3;						// Put this in Hz
  }


        // Input                D	: Dipolar interaction (this)
	//			argc    : Number of arguments
	//			argv    : Array of arguments
	//			qn      : Query value
        //                      Dflag   : Flag is eta is requested
        // Output               none    : D is set interactively
	// Note				: This is INTERACTIVE!

void IntDip::askset(int argc, char* argv[], int& qn, int Dflag)
  {
  double DI, DS, Dcc, DE=0;			// Temp. dip. values
  ask(argc,argv,qn,DI,DS,Dcc,DE,Dflag);		// Use the ask function
  *(this) = IntDip(DI,DS,Dcc,DE);		// Use assignment
  }

        // Input                D       : Dipolar interaction (this)
        //                      Dflag   : Flag is eta is requested
        // Output               none    : D is set interactively
        // Note                         : This is INTERACTIVE (by force)!
 
void IntDip::askset(int Dflag)
  {
  int qn = 1000;				// Temp query index
  int argc = 1;					// Temp # of argments
  char* argv[1];				// Temp array of args
  askset(argc, argv, qn, Dflag);		// Use overload
  }
 

// ____________________________________________________________________________
// I                  DIPOLAR INTERACTION OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//  Functions That Generate String Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------
 
/* string* IntRank2T::TStrings(int M) const;                     INHERITED

                                [ x.x, x.x, x.x]
                        T     = [ x.x, x.x, x.x]
                         2,m    [ x.x, x.x, x.x]
                                [ x.x, x.x, x.x]

  where M = { 0, 1, ..., 4 } ===> m = { 0, 1, -1, 2, -2 }                    */


vector<string> IntDip::CartAStrings(const string& CSForm) const
  {
  int nstr = 6;				// Size of string vector
  vector<string> Cartstrings(nstr);     // Vector of nstr strings
  double sf = 1.0;			// Scaling factor
  string sunit("(Hz)");			// Output units
  double maxDuv = Dmx().maxRe();	// Largest Cartesian element
  if(fabs(maxDuv) >= 1.e6)
    {
    sf = 1.e-6;
    sunit = string("(MHz)");
    }
  else if(fabs(maxDuv) >= 1.e3)
    {
    sf = 1.e-3;
    sunit = string("(kHz)");
    }
  Cartstrings[0] = "[D  , D  , D  ]";
  Cartstrings[1] = "[ xx   xy   xz]"
                 + string("   ")
                 + "["  + Gform(CSForm.c_str(), sf*dxx())
                 + ", " + Gform(CSForm.c_str(), sf*dxy())
                 + ", " + Gform(CSForm.c_str(), sf*dxz()) + "]";
  Cartstrings[2] = "[D  , D  , D  ]"
                 + string(" = ")
                 + "["  + Gform(CSForm.c_str(), sf*dyx())
                 + ", " + Gform(CSForm.c_str(), sf*dyy())
                 + ", " + Gform(CSForm.c_str(), sf*dyz()) + "]";
  Cartstrings[3] = "[ yx   yy   yz]"
                 + string("   ")
                 + "["  + Gform(CSForm.c_str(), sf*dzx())
                 + ", " + Gform(CSForm.c_str(), sf*dzy())
                 + ", " + Gform(CSForm.c_str(), sf*dzz()) + "]";
  Cartstrings[4] = "[D  , D  , D  ]";
  Cartstrings[5] = "[ zx   zy   zz]               " + sunit;
  int i, l, maxl=0;
  for(i=0; i<nstr; i++)
    {
    l = (Cartstrings[i]).length();
    if(maxl < l) maxl = l;
    }
  for(i=0; i<nstr; i++)
    {
    l = (Cartstrings[i]).length();
    if(maxl > l)
      Cartstrings[i] += string(maxl-l, ' ');
    }
  return Cartstrings;
  }

vector<string> IntDip::InfoStrings() const
  {
  vector<string> SphStrings;
  SphStrings.push_back(StringIS());
  string sDCC("Dipolar Coupling:  ");
  string DU(" Hz  ");
  double sf = 1.0;
  if(_DCC > 999.0) { sf=1.e-3; DU = string(" kHz "); }
  sDCC += Gform("%10.2f", _DCC*sf) + DU;
  SphStrings.push_back(sDCC);
  if(!PAS())
    {
    SphStrings.push_back(AlphaString(_EAs.alpha()));
    SphStrings.push_back(BetaString( _EAs.beta() ));
    SphStrings.push_back(GammaString(_EAs.gamma()));
    }
  else
    SphStrings.push_back(PASString());
  SphStrings.push_back(XiString());
  if(ETA) SphStrings[2] = AsymmetryString();
  return SphStrings;
  }      

vector<string> IntDip::DipAStrings() const
 
        // Input                D	: Dipolar interaction (this)
        // Output               CSS     : Pointer to array of 5 strings
        // Note                         : The String array must be deleted
        //                                outside of this routine!

/*               Spin Quantum Numbers:      I   , S
                 Coupling Constant:     xxxxx.xx xHz
                 Asymmetry:                 x.xx
                 Down From PAS z-Axis:    xxx.xx Degrees
                 Over From PAS x-Axis:    xxx.xx Degrees                     */

  {
  vector<string> SphStrings;
  SphStrings.push_back(StringIS());
  string sDCC("DCC:               ");
  string DU(" Hz  ");
  double sf = 1.0;
  if(_DCC > 999.0) { sf=1.e-3; DU = string(" kHz "); }
  sDCC += Gform("%10.2f", _DCC*sf) + DU;
  SphStrings.push_back(sDCC);
  SphStrings.push_back(XiString());
  if(ETA) SphStrings[2] = AsymmetryString();
  return SphStrings;
  }      
 
 
//-----------------------------------------------------------------------------
//     Functions That Generate Ouput Of The Rank 2 Dipolar Interaction
//-----------------------------------------------------------------------------

/* These functions will output information concerning the dipolar interaction
   to any output stream.

           Input                D	: Dipolar interaction (this)
                                ostr	: Output stream
	  			fflag   : Format flag
	  				    0 - Basic Parameters
	  				   !0 - Full output
           Output               none    : Dipolar interaction parameters 
                                          placed into the output stream
	   Note				: This does NOT use the base class
	  				  virtual overload because we write
	  				  out two spin quantum values here

  Function
  ========  ============================================
   print    
   printAT  Spherical components written to output stream
    << 
*/

ostream& IntDip::print(ostream& ostr, bool fflag, bool hdr) const
  {
  if(Izval() < 0.5)				// Just exit if nothing
    {
    string hdr("Empty Dipolar Interaction");
    string spacer = string(40-hdr.length()/2, ' ');
    ostr << "\n\n" << spacer << hdr << "\n";
    return ostr;
    }

/*          Output Some Printed Lines Which Will Look Like The Following

                               Dipolar Interaction

                                          [D  , D  , D  ]
   Spin Quantum Numbers:      x.xxxxxx    [ xx   xy   xz]   [ x.x, x.x, x.x]
   DCC (kHz):                 x.xxxxxx    [D  , D  , D  ] = [ x.x, x.x, x.x]
   Xi Value (rad/sec):        x.xxxxxx    [ yx   yy   yz]   [ x.x, x.x, x.x]
                                          [D  , D  , D  ]
                                          [ zx   zy   zz]                    */

  vector<string> Istrs = InfoStrings();		// Array of info strings
  vector<string> Astrs = CartAStrings("%7.2f");	// Array of A strings
  if(hdr)
    {
    string hdr = "Dipolar Interaction";		// Use this header
    string Spacer((40-hdr.length()/2), ' ');
    ostr << "\n\n" << Spacer << hdr << "\n";
    }
  IntRank2A::print(ostr, Astrs, Istrs);
  if(!fflag) return ostr;

/*          Now Output The Spherical Components, Spatial And Spin
       These Will Be Printed Lines Which Will Look Like The Following
                        (Repeated For All 5 m Values)

                                             [ x.xx  x.xx  x.xx  x.xx]
             A    = x.xxx            T     = [ x.xx  x.xx  x.xx  x.xx]
              2,m                     2,m    [ x.xx  x.xx  x.xx  x.xx]
                                             [ x.xx  x.xx  x.xx  x.xx]      */



  printAT(ostr);				// Print sphercial A & T
  ostr << "\n\n";				// Add some linefeeds
  return ostr;
  }


ostream& IntDip::printAT(ostream& ostr) const

/*	    Prints The Spherical Components, Spatial And Spin
       These Will Be Printed Lines Which Will Look Like The Following 
  		      (Repeated For All 5 m Values)
  
   						[ x.x, x.x, x.x]
  		A    = x.xxx		T     = [ x.x, x.x, x.x]
  		 2,m			 2,m	[ x.x, x.x, x.x]
  						[ x.x, x.x, x.x]             */

  { return IntRank2::printAT(ostr, DIP); }
 

ostream& operator<< (ostream& out, const IntDip& D) { return D.print(out); }

        // Input                out	: Output stream;
        //                      D	: Dipolar tensor to write
        // Output			: Modifies output stream

 
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
//                                  LEGACY
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
 
 

ostream& IntDip::printSpherical(ostream& ostr)
         
        // Input                D	: Dipolar spatial tensor (this)
        //                      ostr	: Output stream
        // Output               none    : Rank 2 spatial tensor parameters
        //                                sent to the output stream
	// Note				: Uses base class virtual overload

  {
  string hdr = "Dipolar Rank2 Spherical Spatial Tensor Components";  
  string Spacer = string((40-hdr.length()/2), ' '); 
  ostr << "\n" << Spacer << hdr;
  IntRank2A::printSpherical(ostr, 0);
  return ostr;
  }

 
ostream& IntDip::printCartesian(ostream& ostr)

        // Input                D	: Dipolar spatial tensor (this)
        //                      ostr	: Output stream
        // Output               none    : Rank 2 spatial tensor parameters
        //                                sent to the output stream
	// Note				: Uses base class virtual overload

  {
  string hdr = "Dipolar Rank2 Cartesian Spatial Tensor Components";  
  string Spacer = string((40-hdr.length()/2), ' '); 
  ostr << "\n" << Spacer << hdr;
  IntRank2A::printCartesian(ostr);
//  IntRank2A::printCartesian(ostr, 0);
// sosi
  return ostr;
  }

//ostream& IntDip::printCartesian(ostream& ostr, double theta, double phi)
         
        // Input                D	: Dipolar spatial tensor (this)
        //                      ostr	: Output stream
        //                      theta   : Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               none    : Rank 2 spatial tensor parameters
        //                                sent to the output stream
	// Note				: Uses base class virtual overload

/*     Axx = (A22+A2m2)/2 - A20/sqrt(6)      Axy = -i*(A22-A21)/2 = Ayx
       Ayy = -(A22+A2m2)/2 - A20/sqrt(6)     Axz = -(A21-A2m1)/2  = Azx
       Azz = sqrt(2/3)*A20                   Ayz = i*(A21-A2m1)/2 = Azy      */  

/*
  {
  string hdr = "Dipolar Rank2 Cartesian Spatial Tensor Components";  
  string Spacer = string((40-hdr.length()/2), ' '); 
  ostr << "\n" << Spacer << hdr;
  IntRank2A::printCartesian(ostr, theta, phi, 0);
  return ostr;
  }
*/

 
// ____________________________________________________________________________
// L                    DIPOLE DIPOLE HAMILTONIAN FUNCTIONS
// ____________________________________________________________________________

/* This section returns irreducible rank 2 dipolar interaction Hamiltonians.
   Because this class uses standardized spatial and spin tensors the
   Hamiltonians may all be derived from the same formula.

                        -2,2
                         ---           m
  H(alpha,beta,gamma) =  \   Xi  * (-1)  * A   (alpha,beta,gamma) * T    (i,j)
                         /                  2,m                      2,-m
                         ---
                          m

   The Hamiltonians will be returned in the spin pair Hilbert space of the
   interaction unless a composite Hilbert space and some spin index is
   supplied. All Hamiltonians will returned in the product basis as simple
   matrices. Their units will be radians/sec.                                */

// ----------------------------------------------------------------------------
//               First Order Dipolar Interaction Hamiltonian
//       These Are SECULAR (Rotationally Invariant About Bo Field Axis)
//  Applicable When The Dipolar Interaction Is A Small Perturbation To Zeeman
// ----------------------------------------------------------------------------

/* Note that H0 is NOT invariant in a multiple rotating frame, but becomes time 
   dependent!  For example, in a heteronuclear system where the return H0 is to
   be in the rotating frame of both I and S then H0 is time dependent from the
   I+S- and the I-S+ terms in T20.  If you work in the laboratory frame, i.e.
   add H0 to the lab-frame Zeeman Hamiltonian you will have not problems. But,
   to work in multiple rotating frames you must NOT use the flip-flop terms.
   The latter assumes a high-field limit!
 
   The secular part of the dipolar Hamiltonian is that which contains only
   those components which commute with z axis rotations.  Here we have
 
                [1]               D   D               D      (0)
               H  (theta,phi) = Xi * A (theta,phi) * T    = H
                D                     2,0             2,0    D
 
	   Input		D	: Dipolar interaction
				alpha   : Euler angle (radians)
				beta    : Euler angle (radians)
				gamma   : Euler angle (radians)
				EA      : Euler angles (radians)
				HSs     : Array of spin Hilbert spaces
				i       : First spin index (in HSs)
				j       : Seond spin index (in HSs)
				wh      : Flag for weak heteronuclear
           Output               H0	: The secular part of the dipolar
                                          Hamiltonian (default basis, Hz)
           Note				: Also called the 1st order dipolar
                                          interaction (perturbation theory)
           Note                      	: Rotationally invariant about z
           Note                         : No HSs: return in the spin Hilbert
                                          space of dimension (2I+1)[*(2S+1)]
                                          With HSs: return in composite spin
                                          Hilbert space                      */

// sosi - should enforce wh = true when isotope values differ....
matrix IntDip::H0(bool wh) const
  { return wh?(_XI*A20())*T20wh:IntRank2::H0(); }

matrix IntDip::H0(double A, double B, double G, bool wh) const
  { return wh?(_XI*A20(A,B,G))*T20wh:IntRank2::H0(A,B,G); }

matrix IntDip::H0(const EAngles& EA, bool wh) const
  { return wh?(_XI*A20(EA))*T20wh:IntRank2::H0(EA); }


matrix IntDip::H0(const vector<int>& HSs, int i, int j, bool wh) const
  { return wh?(_XI*A20())*T20het(HSs,i,j):(_XI*A20())*T20(HSs,i,j); }
//  { return wh?(_XI*A20())*T20het(HSs,i,j):IntRank2::H0(HSs,i,j); }

matrix IntDip::H0(const vector<int>& HSs, int i, int j, 
                  double A, double B, double G, bool wh) const
  { return wh?(_XI*A20(A,B,G))*T20het(HSs,i,j):IntRank2::H0(HSs,i,j,A,B,G); }

matrix IntDip::H0(const vector<int>& HSs, int i, int j, 
                             const EAngles& EA, bool wh) const
  { return wh?(_XI*A20(EA))*T20het(HSs,i,j):IntRank2::H0(HSs,i,j,EA); }

// ----------------------------------------------------------------------------
//                 Full Dipolar Interaction Hamiltonians
//       These Have No Approximations & Are Independent of Field Strength
//     (They Are Correct in The Laboratory Frame, NOT in a Rotating Frame) 
// ----------------------------------------------------------------------------

/* Here we make no approximations. As such, these Hamiltonians are not those
   used in the rotating frame. Here we simply sum over the 5 spatial and spin
   tensor components to produce the Hamiltonian. 

                         -2,2
                          ---   D       m    D                  D
        H (theta, phi) =  \   Xi  * (-1)  * A   (theta, phi) * T    (i,j)
         D                /                  2,m                2,-m
                          ---
                           m 
 
	   Input		D	: Dipolar interaction
				alpha   : Euler angle (radians)
				beta    : Euler angle (radians)
				gamma   : Euler angle (radians)
				EA      : Euler angles (radians)
				HSs     : Array of spin Hilbert spaces
				i       : First spin index (in HSs)
				j       : Seond spin index (in HSs)
           Output               H       : Matrix for dipolar Hamiltonian
           Note                         : No HSs: return in the spin Hilbert
                                          space of dimension (2I+1)[*(2S+1)]
                                          With HSs: return in composite spin
                                          Hilbert space                      */

// matrix IntRank2::H( ) const                                        INHERITED
// matrix IntRank2::H(double alpha, double beta, double gamma) const  INHERITED
// matrix IntRank2::H(const EAngles& EA) const                        INHERITED

// matrix IntDip::H(const vector<int>& HSs, int i, int j) const       INHERITED
// matrix IntDip::H(const vector<int>& HSs, int i, int j,             INHERITED
//                     double alpha, double beta, double gamma) const
// matrix IntDip::H(const vector<int>& HSs, int i, int j,             INHERITED
//                                           const EAngles& EA) const

// ____________________________________________________________________________
// M                 DIPOLE DIPOLE HAMILTONIAN FRIEND FUNCTIONS
// ____________________________________________________________________________

// These Don't Depend on Class IntDip Variables At All.  The Functions Just
// Mimic The Class' Hamiltonian Generation.


matrix HD0(double qn, double wDo, double eta, double theta, double phi)
 
	// Input		qn	: Duantum number (1, 1.5, 2.5,...)
	//			wDo     : PAS Dipolar frequency
	//			eta     : Dipolar asymmetry [0,1]
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               H0	: The secular part of the dipolar
        //                                Hamiltonian (default basis, Hz)
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note				: Also called the 1st order dipolar
        //                                interaction (perturbation theory)
        // Note                      	: Rotationally invariant about z
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1
 
//  The secular part of the dipolar Hamiltonian is that which returns
//  only those components which commute with z axis rotations.  Here we have
 
//            [1]             1                   2             (0)
//           H  (theta,phi) = - w (the,phi) * [3*I - I(I+1)] = H
//            D               6  D                z             D
 
// where 

//                      [ 1      2               1        2                   ]
// w (theta,phi) = W    | - [3cos (theta) - 1] + - eta sin (theta)*cos(2*phi) |
//  D               D,o [ 2                      2                            ]

// and                                                                          
//                                    3*DCC
//                            w    = --------
//                             D,o   2I(2I-1)
 
  {
  int Ival = int(2.*qn + 1);			// For 1 spin SOp functions
  matrix IE = Ie(Ival);				// The identity operator
  matrix IZ = Iz(Ival);				// The Fz operator
  double Ctheta = cos(theta*DEG2RAD);		// Need cos(theta)
  double Stheta = sin(theta*DEG2RAD);		// Need sin(theta)
  double C2phi = cos(2.0*phi*DEG2RAD);		// Need cos(2phi)
  double wD = wDo * ((3.0*Ctheta*Ctheta-1.0)	// "Oriented" wD
            + eta*Stheta*Stheta*C2phi);
  return (wD/12.0)*(3*(IZ*IZ) - (qn*(qn+1))*IE);// 1st Order Hamiltonian
  }  


matrix HD1(double Om, double qn, double wDo, double eta,
                                                      double theta, double phi)
 
        // Input                Om	: Field Strength (Larmor in Hz)
	// 			qn	: Duantum number (1, 1.5, 2.5,...)
	//			wDo     : PAS Dipolar frequency
	//			eta     : Dipolar asymmetry [0,1]
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HD1	: The 2nd order secular part of the
        //                                dipolar Hamiltonian (default
        //                                basis, Hz) for the interaction
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note				: Also called the 1st order dipolar
        //                                interaction (perturbation theory)
        // Note                      	: Rotationally invariant about z
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1
	// Note				: This will be zero in PAS if no eta
	//				  and small if Om >> DCC

//  [2]              -1   [ 1    ]2 [                          2     2   
// H   (theta,phi) = -- * | - w  |  | 2*V V  (theta,phi)*I *(4*I - 8I  - 1)
//  D                Om   [ 6  D ]  [    1 -1             z          z
//
//                                                              2    2      ]
//                                  + 2*V V  (theta,phi)*I *(2*I - 2I  - 1) |
//                                       2 -2             z          z      ]
// where
//
//                   1 [ [    2    2                               ]    4
// 2V V  (the,phi) = - | | eta *cos (2*phi) - 6*eta*cos(2*phi) + 9 | cos (the)
//   1 -1            2 [ [                                         ]
//
//                  [      2    2                                2   ]   2
//                - | 2*eta *cos (2*phi) - 6*eta*cos(2*phi) - eta + 9|cos (the)
//                  [                                                ]
//
//                  [    2      2             ] ]
//                + | eta * [cos (2*phi) - 1] | |
//                  [                         ] ]
//
// and
//                   3 [[  1    2    2          1                  3]    4
//  V V  (the,phi) = - || --*eta *cos (2*phi) - -*eta*cos(2*phi) + -| cos (the)
//   2 -2            4 [[ 12                    2                  4]
//
//                   [ -1    2    2          1    2   3 ]    2
//                 + | --*eta *cos (2*phi) + -*eta  - - | cos (theta)
//                   [  6                    3        2 ]
//
//                   [  1    2    2          1                  3 ] ]
//                 + | --*eta *cos (2*phi) + -*eta*cos(2*phi) + - | |
//                   [ 12                    2                  4 ] ]
 
//  in accordance with the article by P.P. Man "Dipolar Interactions" in
//  the Encyclopedia of Magnetic Resonance by Grant and Harris, Vol 6, Ped-Rel,
//  page 3840, Eq. (19), coupled with page 3841, Eq. (32).

  {
  double C2phi = cos(2.0*phi*DEG2RAD);		// cos(2.0*phi)
  double C2phisq = C2phi*C2phi;			// cos(2*phi)*cos(2*phi)
  double Ctheta = cos(theta*DEG2RAD);		// cos(theta)
  double Cthetasq = Ctheta*Ctheta;		// cos(theta)*cos(theta)
  double Ctheta4 = Cthetasq*Cthetasq;		// (cos(theta))^4
  double etasq = eta*eta;			// eta*eta
  double etasqC2phisq = etasq*C2phisq;		// (eta^2)*(cos(2*phi)^2)
  double etaC2phi = eta*C2phi;			// eta*cos(2*phi)

  double V1Vm11 = (etasqC2phisq - 6.0*etaC2phi + 9.0)*Ctheta4;
  double V1Vm12 = (2.0*etasqC2phisq - 6.0*etaC2phi + - etasq + 9.0)*Cthetasq;
  double V1Vm13 = etasq * (C2phisq - 1.0);
  double twoV1Vm1 = 0.5*(V1Vm11 - V1Vm12 + V1Vm13);

  double V2Vm21 = (etasqC2phisq/12.0 - 0.5*etaC2phi + 0.75)*Ctheta4;
  double V2Vm22 = (-etasqC2phisq/6.0 + etasq/3.0 - 1.5)*Cthetasq;
  double V2Vm23 = (etasqC2phisq/12.0 + 0.5*etaC2phi + 0.75);
  double twoV2Vm2 = 1.5*(V2Vm21 + V2Vm22 + V2Vm23);

  int Ival = int(2.*qn + 1);			// For 1 spin SOp functions
  matrix IZ = Iz(Ival);				// The operator Iz
  matrix IE = Ie(Ival);				// The operator E
  matrix Hmx =  twoV1Vm1*IZ*(4.0*qn*(qn+1)*IE - 8.0*IZ*IZ - IE);
  Hmx       += twoV2Vm2*IZ*(2.0*qn*(qn+1)*IE - 2.0*IZ*IZ - IE);
  Hmx       *= (-wDo*wDo/(36.0*Om));
  return Hmx;
  }
 

#endif						// IntDip.cc
