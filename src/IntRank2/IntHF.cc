/* IntHF.cc *****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Hyperfine Interaction 		         Implementation		**
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
** A variable of type IntHF represents a hyperfine interaction defined	**
** for a particular electron-nucleon spin pair. The interaction is 	**
** maintained in the form of spatial and spin tensors which are stored 	**
** in irreducible spherical format. 					**
**                                                                      **
** The following quantities are inherited/maintained in the base class  **
** IntRank2A (via IntRank2):                                            **
**                                                                      **
**  Asph[5] - Rank 2 irreducible spherical spatial tensor components.   **
**  ETA     - Rank 2 spatial tensor asymmetry value, range [0, 1].      **
**  THETA   - Spatial orientation angle (down from PAS z-axis)          **
**  PHI     - Spatial orientation angle (over from PAS x-axis)          **
**                                                                      **
** Similarly, the following quantities are inherited/maintained in the  **
** base class IntRank2T (via IntRank2):                                 **
**                                                                      **
**  Ival    - Hilbert space of 1st spin involved (typically 2 for e-)   **
**  Sval    - Hilbert space of 2nd spin involved                        **
**  T{0-2}  _ Spin tensor components T20, T21, and T22                  **
**  Tm{1,2} _ Spin tensor components T2-1 and T2-2                      **
**                                                                      **
** Additionally, this class and IntRank2 will maintain:                 **
**                                                                      **
**  I	     - First spin quantum number 				**
**  S        - Second spin quantum number				**
**  DELZZ    - Rank 2 spatial tensor delz value (Gauss).                **
**  AISO     - Isotropic hyperfine coupling (Gauss).                    **
**                                                                      **
** The tensor will internally be stored in irreducible spherical form   **
** at any chosen orientation.  The tensor principal axes will also be   **
** maintiained as well as a set of Euler angles that indicate the PAS   **
** orientation relative to some common axes.                            **
**                                                                      **
** In addition to the functions contained in the base spatial tensor    **
** class, this IntHF provides functions for building up the tensor      **
** accessing the tensor elements from a hyperfine standpoint.  This     **
** includes functions for reorienting the tensor, making Hamiltonians,  **
** and obtaining electron resonance conditions.                         **
**                                                                      **
** The following defintions are used herein (Auv normlized by delzz):   **
**                                                                      **
** 1.) PAS: |Azz| >= |Ayy| >= |Axx|                                     **
** 2.) PAS: Azz=1, eta=(Ayy-AXX)/2, Axx=(1+eta)/2, Ayy=(eta-1)/2        **
** 2.) PAS: Azz=2C, eta=(Axx-Ayy)/Azz, Axx=C(eta-1), Ayy=-C(1+eta)	**
**                                                                      **
** Individual spin tensors are stored in a linked list, handled by the  **
** class IntSTLList.  This prevents recalculation of the 5 arrays that  **
** are associated with such tensors when the spin(s) involved share the **
** same Iz & Sz values.  However, the arrays are still copied into new  **
** (equivalent) spin tensors, but that is done by reference in the      **
** GAMMA matrix class.                                                  **
**                                                                      **
*************************************************************************/

#ifndef   IntHF_cc_			// Is file already included?
#  define IntHF_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <IntRank2/IntHF.h>		// Include interface definition
#include <Basics/Gconstants.h>		// Include PI and other constants
#include <Basics/Isotope.h>		// Include isotopes
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/Gutils.h>		// Include parameter queries
#include <HSLib/SpinOpSng.h>		// Include 1 spin operators
#include <Matrix/row_vector.h>		// Include row_vectors
#include <Matrix/matrix.h>		// Include matrices
#include <Level1/coord.h>		// Include coordinates
#include <HSLib/SpinSys.h>		// Include base spin systems
#include <Basics/StringCut.h>
#if defined(_MSC_VER) || defined(__SUNPRO_CC)				// If we are using MSVC++
 #pragma warning (disable : 4786)       //   Kill STL namelength warnings
#endif

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ STL lists
using std::vector;			// Using libstdc++ STL vectors
using std::ofstream;			// Using libstdc++ output file streams
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i               CLASS HYPERFINE INTERACTION ERROR HANDLING
// ____________________________________________________________________________

/*       Input                HF      : Hyperfine interaction (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void IntHF::IHFerror(int eidx, int noret) const
  {
  string hdr("Hyperfine Interaction");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr,"Program Aborting.....",        noret); break;// (0)
    case 1: GAMMAerror(hdr,"Construction Spin Pair Bad",   noret); break;// (1)
    case 2: GAMMAerror(hdr,"Problems During Construction", noret); break;// (2)
    case 3: GAMMAerror(hdr,"Problems During Assignment.",  noret); break;// (3)
    case 4: GAMMAerror(hdr,"Demands Electron-Nucleon Pair",noret); break;// (4)
    case 8: GAMMAerror(hdr,"Theta (z Down) Beyond [0,180]",noret); break;// (8)
    case 9: GAMMAerror(hdr,"Phi (x Over) Outside [0, 360]",noret); break;// (9)
    case 10:GAMMAerror(hdr,"Asymmetry (eta) Beyond [0, 1]",noret); break;// (10)
    case 12:GAMMAerror(hdr,"Cant Set Hyperfine Anisotropy",noret); break;// (12)
    case 13:GAMMAerror(hdr,"Cannot Set Hyperfine Coupling",noret); break;// (13)
    case 14:GAMMAerror(hdr,"Cartesian Cmpt. Order Wrong",  noret); break;// (14)
    case 15:GAMMAerror(hdr,"Use |Azz| >= |Ayy| >= |Axx|",  noret); break;// (15)
    case 16:GAMMAerror(hdr,"Setting Spin Type To Proton",  noret); break;// (16)
    case 17:GAMMAerror(hdr,"Setting Spin Type To Electron",noret); break;// (17)
    case 20:GAMMAerror(hdr,"Can't Write To Parameter File",noret); break;// (20)
    case 21:GAMMAerror(hdr,"Cant Read From Parameter File",noret); break;// (21)
    case 22:GAMMAerror(hdr,"Insufficient File Parameters", noret); break;// (22)
    case 23:GAMMAerror(hdr,"Insufficient PSet Parameters", noret); break;// (23)
    case 25:GAMMAerror(hdr,"Trouble Setting Spin Tensor",  noret); break;// (25)
    case 26:GAMMAerror(hdr,"Using Default Quantum Numbers",noret); break;// (26)
    case 27:GAMMAerror(hdr,"Sorry, Setting I = S = 1/2",   noret); break;// (27)
    case 53:GAMMAerror(hdr,"Cannot Output Parameters",     noret); break;// (53)
    case 60:GAMMAerror(hdr,"Use A Dipolar Interaction?",   noret); break;// (60)
    }
  }

volatile void IntHF::IHFfatal(int eidx) const
  {
  IHFerror(eidx, 1);				// Output error message
  if(eidx) IHFerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }


/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Funciton pname
      5     Please Use Class Member Funciton pname
   default  Unknown Error - pname                                           */

void IntHF::IHFerror(int eidx, const string& pname, int noret) const 
  {
  string hdr("Hyperfine Interaction");
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
    case 101:                                                   // (101)
      msg = string("Can't Find Interaction Parameters For ")
          + pname;
      GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;
    }
  }

volatile void IntHF::IHFfatal(int eidx, const string& pname) const
  {
  IHFerror(eidx, pname, 1);			// Output error message
  if(eidx) IHFerror(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }


// ____________________________________________________________________________
// ii                HYPERFINE INTERACTION SETUP FUNCTIONS
// ____________________________________________________________________________

/* These by-pass the generalized GAMMA spatial-spin tensors (in class spin_T)
   and just produces the rank 2 space-spin tensor used for the hyperfine
   interaction directly.  They are scaled such that they are interaction
   independent and adhere to the relative scaling between components of true
   spherical tensors.  The rank 2 spatial tensors (used by this class) are
   scaled such that, when rotated they are normalized rank 2 spherical
   harmonics in a symmetric case (eta=0).                                    */

void IntHF::setT20wh()

        // Input                HF      : Hyperfine interaction (this)
        // Output               none    : Hyperfine interaction weak
        //                                heteronuclear spin tensor
        //                                components is generated
        //
        //                                     T   | = T20wh
        //                                      2,0|
        //                                          weak heteronuclear
        //
        // Note                         : This should be called during
        //                                all constructors after IntRank2T
        //                                has been generated
	// Note				: Unless we are near zero field
	//				  this component is will be valid

//          1/2                                1/2             1/2
//       [1]                      weak      [1]             [4]
// T   = |-| * [3I S - I.S]   ------------> |-| * [2I S ] = |-| I S
//  2,0  [6]      z z         heteronuclear [6]      z z    [6]  z z

  {
  matrix IE = Ie(Ival);                         // The operator Ie
  matrix SE = Ie(Sval);                         // The operator Se
  matrix IZ = tensor_product(Iz(Ival), SE);     // The operator Iz
  matrix SZ = tensor_product(IE, Iz(Sval));     // The operator Sz
  T20wh  =  (2.*IZ*SZ)/sqrt(6.);                // T20  = (2IzSz)/sqrt(6)
  }

// ____________________________________________________________________________
// iii        HYPERFINE INTERACTION PARAMETER SET SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions are used to parse a parameter set and fill up an electron
   HF interaction. They just route into others that read specific parts of the
   interaction. They oversee the process & quit if things go awry.  Since
   these are somewhat inter-related, we keep these functions private to avoid
   misuse.

   The goal here is quite simple. We must determine the following set of values
   for each hyperfine interaction: { Iqn,Sqn,A,Adelz,eta,alpha,beta,gamma }
   Complexity arises because we allow that a variety of parameters to be used
   to define these values. Additionally we allow that some parameters may be
   zero (e.g. eta) and that defaults may automatically be used. But, the end
   result remains the same, we seek the set of values that define a
   hyperfine interaction in GAMMA.                                           */

// ----------------------------------------------------------------------------
//                        Complete Hyperfine Interaction
// ----------------------------------------------------------------------------

/* These functions will try and get all the parameters required to define a
   hyperfine interaction: { Iqn,Sqn,A,Adelz,eta,alpha,beta,gamma }. Although
   the 1st function is the generic means to accomplishing this, we will
   immediately branch into two categories: 1.) Defining the interaction with a
   single interaction index & 2.) Defining the interaction with two spin
   indices.

           Input        HF      : Hyperfine interaction (this)
                        pset    : A parameter set
                        Iqn     : 1st spin quantum number
                        Sqn     : 2nd spin quantum number
                        hfc	: Isotropic hyperfine coupling (Gauss)
                        hfa     : Hyperfine anisotropy value (Gauss)
                        aeta    : Hyperfine asymmetry value [0,1]
                        EA      : Euler angles for orientation (radians)
                        idxI    : 1st spin or interaction index
                        idxS    : 2nd spin index
                        warn    : Warning level
                                    0 - no warnings
                                    1 - warnings
                                   >1 - fatal warnings
       Output           TF      : True if hyperfine interaction defined
                                  Argument values are set                    */

bool IntHF::getHFI(const ParameterSet& pset, double& Iqn, double& Sqn,
         double& hfc, double& hfa, double& aeta, EAngles& EA,
                                         int iI, int iS, bool warn) const
  {
  if(iS == -1) return getHFI1(pset,Iqn,Sqn,hfc,hfa,aeta,EA,iI,   warn);
  else         return getHFI2(pset,Iqn,Sqn,hfc,hfa,aeta,EA,iI,iS,warn);
  }

bool IntHF::getHFI1(const ParameterSet& pset, double& Iqn, double& Sqn,
         double& hfc, double& hfa, double& eta, EAngles& EA,
                                                    int idxI, bool warn) const
  {
  string pn("A");				// Parameter name base

//                First Determine The Two Spin Quantum Numbers
//    ( Default Is I=S=1/2, This Sets Hyperfine Spin Tensor Hilbert Space )

  bool TFI = getIqn(pset, "Iqn", Iqn, idxI, 0);         // Try to get Iqn
  if(!TFI) Iqn = 0.5;					// If failed, Iqn=1/2
  bool TFS = getIqn(pset, "Sqn", Sqn, idxI, 0);         // Try to get Sqn
  if(!TFS) Sqn = 0.5;					// If failed, Sqn=1/2

// Try To Directly Read Hyperfine Coupling, Anisotropy, Asymmetry, Orientation
//    Hypefine Interaction Via { Iqn, Sqn, A(i), AA(i), Aeta(i), AEAs(i) }

//  1.) If A has been specified, these parameters will be used
//  2.) We don't mind that eta is not set, default will be zero
//  3.) Iqn & Sqn were set in the 1st section of this function
//  4.) We don't mind that no orientation is set, default is PAS
//  5.) Orientation set with either an Euler angle set or 3 individual angles

  if(getHFC(pset, hfc, idxI, -1, 0))
    {                                                   //     DCC, eta, Orient
    getHFA(pset, hfa, idxI, -1, 0);			//   Loof for anisotropy
    getAeta(pset, pn, eta, idxI, -1, 0);		//   Look for asymmetry
    getOrientation(pset,pn,EA,idxI,-1,0);		//   Look for orient.
    return true;
    }

//      Try To Directly Read Hypefine Cartesian Spatial Tensor Components
//     Hypefine Interaction Via { Iqn, Sqn, Axx, Axy, Axz, Ayy, Ayz, Azz }

//  1.) Auv set the hyperfine coupling, anisotropy, & asymmetry, A, AA, & eta
//  2.) If any Auv off-diagonals (u != v) set, A array sets orientation also
//  3.) If no Auv off-diagonals, orientation set with specified Euler angles
//  4.) The input Auv values are taken to be in Gauss
//  5.) We don't mind that no orientation is set, default is PAS
//  6.) If used, Euler angles by either a set or 3 individual angles

  coord AiAzAe;                                 // For Aiso, Azz, Aeta
  if(getACart(pset,"A",AiAzAe,EA,idxI,-1,0))	// Try & use Cart. components
    {
    hfc = AiAzAe.x();                           //  Set hyperfine coupling
    hfa = AiAzAe.y();                           //  Set anisotropy value
    eta = AiAzAe.z();                           //  Set the eta value
    return true;
    }

  return false;
  }

bool IntHF::getHFI2(const ParameterSet& pset, double& Iqn, double& Sqn,
         double& hfc, double& hfa, double& eta, EAngles& EA,
                                          int idxI, int idxS, bool warn) const
  {
//                First Determine The Two Spin Quantum Numbers
//    ( Default Is I=S=1/2, This Sets Hyperfine Spin Tensor Hilbert Space )

  string II, IS;                                // String for isotope names
  Isotope ISI, ISS;                             // Isotopes for spins
  string pn("A");                               // Parameter name base
  bool TFI = false;                             // Flag if we know isotopes
  if(getIsos(pset,idxI,idxS,II,IS,0))           // 1st try for isotope names
    {                                           // If successful, check them
    if(!SpinCheck(II,IS)) return false;         //    Insure valid isotopes
    ISI = Isotope(II);                          //    An isotope of type II
    ISS = Isotope(IS);                          //    An isotope of type IS
    if(!SpinCheck(ISI,ISS,!TFI,0)) return false;//    Allow nucleus/e- pair
    TFI = true;                                 //    We know both isotopes
    Iqn = ISI.qn();                             //    Know spin I quantum #
    Sqn = ISS.qn();                             //    Know spin S quantum #
    }
  else if(!getIqns(pset,pn,Iqn,Sqn,idxI,idxS,0))// 2nd try for spin quant. #s
    { Iqn = 0.5; Sqn = 0.5; }                   // Use default if not able

//  Try To Directly Read Hypefine Coupling, Anisotropy, Asymmetry, Orientation
// Hypefine Interaction Via { Iqn, Sqn, A(i,j), AA(i,j), Aeta(i,j), DEAs(i,j) }

//  1.) If A has been specified, these parameters will be used
//  2.) We don't mind that eta is not set, default will be zero
//  3.) Iqn & Sqn were set in the 1st section of this function
//  4.) We don't mind that no orientation is set, default is PAS
//  5.) Orientation set with either an Euler angle set or 3 individual angles

  if(getHFC(pset, hfc, idxI, idxS, 0))
    {
    getHFA(pset, hfa, idxI, idxS, 0);			//   Loof for anisotropy
    getAeta(pset, pn, eta, idxI, idxS, 0);		//   Look for asymmetry
    getOrientation(pset,pn,EA,idxI,idxS,0);		//   Look for orient.
    return true;
    }

//      Try To Directly Read Hypefine Cartesian Spatial Tensor Components
//     Hypefine Interaction Via { Iqn, Sqn, Axx, Axy, Axz, Ayy, Ayz, Azz }

//  1.) Auv set the hyperfine coupling, anisotropy, & asymmetry, A, AA, & eta
//  2.) If any Auv off-diagonals (u != v) set, A array sets orientation also
//  3.) If no Auv off-diagonals, orientation set with specified Euler angles
//  4.) The input Auv values are taken to be in Gauss
//  5.) We don't mind that no orientation is set, default is PAS
//  6.) If used, Euler angles by either a set or 3 individual angles

  coord AiAzAe;                                 // For Aiso, Azz, Aeta
  if(getACart(pset,"A",AiAzAe,EA,idxI,idxS,0))  // Try & use Cart. components
    {
    hfc = AiAzAe.x();				//  Set hyperfine coupling
    hfa = AiAzAe.y(); 				//  Set anisotropy value
    eta = AiAzAe.z();                           //  Set the eta value
    return true;
    }
  return false;
  }

// ----------------------------------------------------------------------------
//         Get Isotropic Hypefine Coupling Value From A Parameter Set
// ----------------------------------------------------------------------------

/*	Input           HFI	: Hyperfine interaction (this)
			pset	: A parameter set
                        idxI	: 1st spin or interaction index
                        idxS	: 2nd spin index
                        hfc     : Hypefine coupling (Gauss)
                        warn    : Warning output flag
	Output          TF	: True if hyperfine coupling obtained
                                  constant from parameters in pset
	Note			: This WILL NOT alter the interaction
	Note			: Parameters are A, A(#), A(#,#)            */

bool IntHF::getHFC(const ParameterSet& pset, double& hfc,
                                          int idxI, int idxS, bool warn) const
  {
  ParameterSet::const_iterator item;         // A pix into parameter list
  string Nidx = "";                             // Name addition per index
  if(idxI != -1) 				// If index exists then set up
    {						// a proper parameter name
    Nidx += string("(") + Gdec(idxI);		// suffix
    if(idxS != -1)
      Nidx += string(", ") + Gdec(idxS); 
    Nidx += string(")");
    }

//          Look For Isotropic Hypefine coupling Using Parameter A

  string pname =  string("A") + Nidx;		// Hyperfine coupling (G) 
  item = pset.seek(pname);                      // Seek parameter in pset
  if(item != pset.end())                        // If parameter was found
    {                                           // parse the parameter
    string pstate;                              // Temp String for parameter
    (*item).parse(pname,hfc,pstate); 		// Glean out parameter value
    return true;                                // Return the coupling value
    }

  if(warn)					// If it hasn't been found
    {                                           // we can make an interaction
    IHFerror(2, pname, 1);			// Can't find A(idxI,idxS)
    IHFerror(13, 1);				// Can't set hyperfine coupling
    }
  return false;                                 // Could not read HF coupling
  }

// ----------------------------------------------------------------------------
//            Get Hypefine Anisotropy Value From A Parameter Set
// ----------------------------------------------------------------------------

/*	Input           HFI	: Hyperfine interaction (this)
			pset	: A parameter set
                        idxI	: 1st spin or interaction index
                        idxS	: 2nd spin index
                        hfa     : Hypefine anisotropy (Gauss)
                        warn    : Warning output flag
	Output          TF	: True if hyperfine coupling obtained
                                  constant from parameters in pset
	Note			: This WILL NOT alter the interaction
	Note			: Parameters are AA, AA(#), AA(#,#)         */

bool IntHF::getHFA(const ParameterSet& pset, double& hfa,
                                          int idxI, int idxS, bool warn) const
  {
  ParameterSet::const_iterator item;         // A pix into parameter list
  string Nidx = "";                             // Name addition per index
  if(idxI != -1) 				// If index exists then set up
    {						// a proper parameter name
    Nidx += string("(") + Gdec(idxI);		// suffix
    if(idxS != -1)
      Nidx += string(", ") + Gdec(idxS); 
    Nidx += string(")");
    }

//              Look For Hypefine Anisotropy Using Parameter AA

  string pname =  string("AA") + Nidx;		// Hyperfine coupling (G) 
  item = pset.seek(pname);                      // Seek parameter in pset
  if(item != pset.end())                        // If parameter was found
    {                                           // parse the parameter
    string pstate;                              // Temp String for parameter
    (*item).parse(pname,hfa,pstate); 		// Glean out parameter value
    return true;                                // Return the anisotropy
    }

  if(warn)					// If it hasn't been found
    {                                           // we can make an interaction
    IHFerror(2, pname, 1);			// Can't find A(idxI,idxS)
    IHFerror(12, 1);				// Can't set anisotropy
    }
  return false;                                 // Could not read anisotropy
  }

// ----------------------------------------------------------------------------
//                        Complete Hyperfile Interaction
// ----------------------------------------------------------------------------

/* This function employs all of the "get*" functions in the above sections
   to parse a parameter set for all of the values needed to define a hyperfine
   interaction, namely { Iqn,Sqn,A,Adelz,eta,alpha,beta,gamma }. If the
   interaction definition is found, we set the interaction or return false.  */

bool IntHF::setHFI(const ParameterSet& pset, int idxI, int idxS, int warn)
  {
  double  Iz;					// 1st spin quantum number
  double  Sz;					// 2nd spin quantum number
  double  hfc;                                  // Hyperfine coupling   (Gauss)
  double  hfa;                                  // Hyperfine anisotropy (Gauss)
  double  E;					// Hypefine asymmetry     [0,1]
  EAngles EA;                                   // Our Euler angles (radians)
  if(getHFI(pset, Iz, Sz, hfc, hfa, E,		// Try and get parameter vals
                         EA, idxI, idxS, warn?true:false))
    {                                           // If successful, set interact
    AISO     = hfc;				//   Set isotropic value (Gauss)
    DELZZ    = hfa/1.5;				//   Set PAS delzz value (Gauss)
    double X = xi();                            //   Set interaction constant
    IntRank2::operator=(IntRank2(Iz,Sz,X,E,EA));//   Set rank 2 interaction
    return true;                                //   Return we are successful
    }
  return false;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A          HYPERFINE INTERACTION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/* Hyperfine interaction construction. For full construction users need to
   define the 1.) rank 2 spatial tensor, 2.) rank 2 spin tensors, and 
   3.) interaction strength.  The 1st and 3rd can be supplied in spherical 
   terms as { Aiso, Azz, eta } or in Cartesian terms as { hxx, hyy, hzz }.
   For the spin tensor we need the two spin quantum value involved.

   GAMMA's rank 2 interactions are always stored in their PAS and as normalized
   irreducible spherical components. Scaling of the interaciton is tied up in
   the isotropic value (i.e. the A value).  Since Aiso = (1/3)[hxx+hyy+hzz]
   if the interaction is specified using these Cartesian components this is
   automatically set. If spherical components are used the A value must be 
   input.  The irreducible anisotropic rank 2 part of the spatial tensor will
   either be set by {Azz , eta} or from the reducible Cartesian rank 2 PAS
   components { hxx, hyy, hzz }.                                             */

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

IntHF::IntHF() : IntRank2()
  {
  AISO  = 0.0;					// No isotropic component
  DELZZ = 0.0;					// No anisotropic component
  }

IntHF::IntHF(const IntHF &HF1) : IntRank2(HF1) 
  {
  AISO  = HF1.AISO;				// Copy isotropic component
  DELZZ = HF1.DELZZ;				// Copy anisotropic component
  T20wh = HF1.T20wh;				// Copy T20 weak heteronuclear
  }

// ----------------------------------------------------------------------------
//              Direct Constructors Using PAS Cartesian Components
// ----------------------------------------------------------------------------

/* For this we need to know the 3 (reducible) PAS spatial tensor components,
   { Axx, Ayy, Azz } as well as the two spin quantum numbers & an orientation.
   It is assumed that the input Auv are in Gauss. Note that the three Cartesian
   values { Axx, Ayy, Azz } are equivalent to the three spherical values
   { Aiso, Adelzz, Aeta }.  Both specify the reducible rank 2 hyperfine spatial
   tensor in its principal axes.

        Input      HF      : Hyperfine interaction (this)
                   AxAyAz  : PAS HF Cartesian cmpts (Gauss)
                   IsoI    : Spin I isotope type
                   IsoS    : Spin S isotope type
        Output     none    : Hyperfine interaction constructed
                             from Axx, Ayy, Azz contained in
                             the input coordinate (in Gauss)                 */

IntHF::IntHF(const string& II,const string& IS,
                                          const coord& Axyz, const EAngles& EA)
  {
  if(!SpinCheck(II,IS)) IHFfatal(2);		// Insure spin types valid
  Isotope III(II);				// Make isotope for II
  Isotope IIS(IS);				// Make isotope for IS
  if(!SpinCheck(II,IS,true)) 			// Insure nucleus/e- pair
   { IHFerror(60, 1); IHFfatal(2); }		// Disallow e-/e- n/n pairing
  double Iz = III.qn();				// Get Iz value of II
  double Sz = IIS.qn();				// Get Iz value of IS
  coord  ADE = IntRank2A::AisoDelzEta(Axyz);	// Get Aiso, DELZZ, Eta
  AISO  = ADE.x();				// Set isotropic value 
  DELZZ = ADE.y();				// Set delzz value
  double X = xi();				// Get interaction constant
  double E = ADE.z();				// Get asymmetry (eta)
  IntRank2::operator=(IntRank2(Iz,Sz,X,E,EA));	// Use generic interaction
  setT20wh();					// Create added T2 component
  }

IntHF::IntHF(const Isotope& II,const Isotope& IS,
                                          const coord& Axyz, const EAngles& EA)
  {
  if(!SpinCheck(II,IS,true)) 			// Insure nucleus/e- pair
   { IHFerror(60, 1); IHFfatal(2); }		// Disallow e-/e- n/n pairing
  double Iz = II.qn();				// Get Iz value of II
  double Sz = IS.qn();				// Get Iz value of IS
  coord  ADE = IntRank2A::AisoDelzEta(Axyz);	// Get Aiso, DELZZ, Eta
  AISO  = ADE.x();				// Set isotropic value 
  DELZZ = ADE.y();				// Set delzz value
  double X = xi();				// Get interaction constant
  double E = ADE.z();				// Get asymmetry (eta)
  IntRank2::operator=(IntRank2(Iz,Sz,X,E,EA));	// Use generic interaction
  setT20wh();					// Create added T2 component
  }

IntHF::IntHF(double Iz, double Sz, const coord& Axyz, const EAngles& EA)
  {
  if(!SpinCheck(Iz,Sz,1)) { IHFfatal(2); }	// Insure Iz & Sz valid (n*1/2)
  coord ADE = IntRank2A::AisoDelzEta(Axyz);	// Get Aiso, DELZZ, Eta
  AISO  = ADE.x();				// Set isotropic value 
  DELZZ = ADE.y();				// Set delzz value
  double X = xi();				// Get interaction constant
  double E = ADE.z();				// Get asymmetry (eta)
  IntRank2::operator=(IntRank2(Iz,Sz,X,E, EA));	// Use generic interaction
  setT20wh();					// Create added T2 component
  }

// ----------------------------------------------------------------------------
//               Direct Constructors Using Spherical Components
// ----------------------------------------------------------------------------

/* For this we will need to know the three PAS spatial tensor components,
   { Aiso, Adelzz, Aeta } and the two spin quantum numbmers. It is assumed 
   that the input of both Aiso and Adelzz are in Gauss. Aeta will be
   restricted to reside between [0,1]

        Input       HF      : Hyperfine interaction (this)
                    IsoI    : Spin I isotope type
                    IsoS    : Spin S isotope type
        or          Iqn     : Spin I quantum number
                    Sqn     : Spin S quantum number
                    Aiso    : Spatial tensor Aiso value (Gauss)
                    Adelzz  : Spatial tensor delzz value (Gauss)
                    Aeta    : Spatial tensor asymmetry value (usually 0)
        Output      none    : Hyperfine interaction constructed              */


IntHF::IntHF(const string& II,const string& IS,
                    double Aiso, double Adelzz, double Ae, const EAngles& EA)
  {
  if(!SpinCheck(II,IS,1)) IHFfatal(2);		// Insure spin types valid
  Isotope III(II);				// Make isotope for II
  Isotope IIS(IS);				// Make isotope for IS
  *this = IntHF(III,IIS,Aiso,Adelzz,Ae,EA);	// User constructor overload
  }

IntHF::IntHF(const Isotope& II,const Isotope& IS,
                    double Aiso, double Adelzz, double Ae, const EAngles& EA)
  {
  if(!SpinCheck(II,IS,1,1))			// Insure e-/nucleus pair
   { IHFerror(60, 1); IHFfatal(2); }		// Disallow e-/e- n/n pairing
  double Iz = II.qn();				// Get Iz value of II
  double Sz = IS.qn();				// Get Iz value of IS
  AISO  = Aiso;					// Set isotropic value 
  DELZZ = Adelzz;				// Set delzz value
  double X = xi();				// Get interaction constant
  IntRank2::operator=(IntRank2(Iz,Sz,X,Ae,EA));	// Use generic interaction
  setT20wh();					// Create added T2 component
  }

IntHF::IntHF(double Iz, double Sz,
                    double Aiso, double Adelzz, double Ae, const EAngles& EA)
  {
  if(!SpinCheck(Iz,Sz,1)) { IHFfatal(2); }	// Insure Iz & Sz valid
  AISO  = Aiso;					// Set isotropic value 
  DELZZ = Adelzz;				// Set delzz value
  double X = xi();				// Get interaction constant
  IntRank2::operator=(IntRank2(Iz,Sz,X,Ae,EA));	// Use generic interaction
  setT20wh();					// Create added T2 component
  }

// ----------------------------------------------------------------------------
//                     Construction Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions will construct the interaction from parameters taken from a
   specified GAMMA parameter set. This is most useful when working with multi-
   spin systems that are setting up many interactions over the system using a
   single parameter file (or parameters read from an external ASCII file. The
   constructor using spin isotope types in conjunction with a parameter set is
   expressly for use in multispin systems (see IntHFVec & SolidSys)

       Input                HF      : Hyperfine interaction
                            pset    : Parameter set
                            idx     : Interaction index (default -1->none)
                            IdxI    : Index for spin I
                            IdxS    : Index for spin S
                            warn    : Flag to warn if no interaction found
       Output               none    : HF interaction constructed
                                      for spins with index idxI & idxS
                                      whose parameters are in pset
       Output               or      : HF interaction constructed
                                      for spins with quantum numbers qI & qS
                                      and parameters in pset                */
 
IntHF::IntHF(const ParameterSet& pset, int idxI, int idxS, int warn)
  {
  if(!setHFI(pset, idxI, idxS, warn?1:0))	// Set using parse function
    {
    if(warn)					// If failure, warn if desired
      {						// and maybe even abort
      IHFerror(23, 1);				//   Insufficient pset
      if(warn > 1) IHFfatal(2);			//   Problems in construction
      else         IHFerror(2);
      }
    }
  }


/*
IntHF::IntHF(const Isotope& II, const Isotope& IS, 
                        const ParameterSet& pset, int idxI, int idxS, int warn)
  {
  if(!nepair(II, IS))			// Insure electron-nucleon
    {					// If two nuclei or two electrons
    if(warn)				// issue warnings, maybe die
      {
      IHFerror(2, 1);			//   Problems during construction
      IHFerror(4, 1);			//   Need electron-nucleon pair
      if(warn > 1) IHFfatal(1);	//   Bad intraction spin pair
      }
    AISO = 0.0;
    return;
    }
  I    = II.qn();			//   Set I spin quantum value
  S    = IS.qn();			//   Set S spin quantum value
  Ival = II.HS();			//   Set I spin Hilbert space
  Sval = IS.HS();			//   Set S spin Hilbert space
  setHF();				//   Set up our spin tensor

  if(!setCartComp(pset,idxI,idxS))	// Try & read Cartesian components
    if(!setSphComp(pset,idxI,idxS))	// try & read spherical components
      AISO = 0;
    else
      {
      setAs();				// Set A2,m based on THETA,
      AISO = iso();			// Set the hyperfine isotropic value
      }
  }



IntHF::IntHF(const ParameterSet& pset, int idxI, int idxS, int warn)
  { getHFI(pset, idxI, idxS, warn?1:0); }

IntHF::IntHF(int idxI, int idxS, const ParameterSet& pset, int warn)
  { getHFI(pset, idxI, idxS, warn?1:0); }
 
IntHF::IntHF(const string& II, const string&IS, 
                        const ParameterSet& pset, int idxI, int idxS, int warn)
  { *this = IntHF(Isotope(II), Isotope(IS), pset, idxI, idxS, warn); }

*/

// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

void IntHF::operator= (const IntHF &HF1) 
  {
  IntRank2::operator=(HF1);		// Copy generic rank 2 interaction
  AISO  = HF1.AISO;			// Copy the hyperfine coupling 
  DELZZ = HF1.DELZZ;			// Copy the hyperfine anisotropy
  T20wh = HF1.T20wh;			// Copy weak heteronuclear T20
  }

IntHF::~IntHF () { }

// ____________________________________________________________________________
// B          HYPERFINE INTERACTION SPATIAL TENSOR ACCESS FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                      Spatial Tensor Isotropic Value
// ----------------------------------------------------------------------------

/* Here we only allow the hyperfine coupling to be set in Gauss. The coupling
   is define as
                                  1 [               ]
                           A    = _ | h  + h  + h   |
                            iso   3 [  xx   yy   zz ]                        */

double IntHF::iso()            const { return AISO; }
void   IntHF::iso(double aiso)       { AISO = aiso; }
double IntHF::A()              const { return AISO; }
void   IntHF::A(double aiso)         { AISO = aiso; }

// ----------------------------------------------------------------------------
//                      Spatial Tensor Anisotropy Value
// ----------------------------------------------------------------------------

/* The hyperfine anisotropy is

                  ^             1                 3
                 / \ A  = h   - - [ h   + h   ] = - del
                 ---       zz   2    xx    yy     2    zz

                           del   = h   - A
                                zz    zz    iso

   Note that the default base class IntRank2A also provides the GAMMA
   normalized anisotropy and delzz values. These are constants given by

                       1/2
                  [ 5 ]                  ^      3
          del   = |---]  = 0.51503      /_\ A = - del   = 0.77255
             zz   [6PI]                         2    zz 

                Input           HF      : Hyperfine interaction (this)
                                delzz   : The HF anisotropy value (Gauss)
                Output          delzz   : HF tensor delzz value             */


// static double IntRank2A::delzz();                            INHERITED
// static double IntRank2A::delA();                             INHERITED

double IntHF::aniso()            const { return 1.5*DELZZ; }
void   IntHF::aniso(double aiso)       { DELZZ = aiso/1.5; }
double IntHF::AA()               const { return 1.5*DELZZ; }
void   IntHF::AA(double aiso)          { DELZZ = aiso/1.5; }


// ----------------------------------------------------------------------------
//                      Spatial Tensor Asymmetry Value
// ----------------------------------------------------------------------------

// Note that eta spans [0, 1] and is defined to be (Axx-Ayy)/Azz where
// |Azz| >= |Ayy| >= |Axx| in the treatment contained herein.  These functions
// are inherited from the base class IntRank2A

//double IntHF::eta( ) const             INHERITED       Get the asymmetry
//void   IntHF::eta(double HFeta)        INHERITED       Set the asymmetry

//-----------------------------------------------------------------------------
//                    Spherical Tensor Component Access
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
complex IntRank2::A20PAS(  ) const;                             INHERITED
complex IntRank2::A21PAS(  ) const;                             INHERITED
complex IntRank2::A2m1PAS( ) const;                             INHERITED
complex IntRank2::A22PAS(  ) const;                             INHERITED
complex IntRank2::A2m2PAS( ) const;                             INHERITED

complex IntRank2::A20()  const;                                 INHERITED
complex IntRank2::A21()  const;                                 INHERITED
complex IntRank2::A2m1() const;                                 INHERITED
complex IntRank2::A22()  const;                                 INHERITED
complex IntRank2::A2m2() const;                                 INHERITED

complex IntRank2A::A20(double  A, double B, double G) const;    INHERITED
complex IntRank2A::A21(double  A, double B, double G) const;    INHERITED
complex IntRank2A::A2m1(double A, double B, double G) const;    INHERITED
complex IntRank2A::A22(double  A, double B, double G) const;    INHERITED
complex IntRank2A::A2m2(double A, double B, double G) const;    INHERITED

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
//                      Cartesian Tensor Component Access
//-----------------------------------------------------------------------------

/*
  double IntRank2A::Axx() const;                                INHERITED
  double IntRank2A::Ayy() const;                                INHERITED
  double IntRank2A::Azz() const;                                INHERITED
  double IntRank2A::Ayx() const;                                INHERITED
  double IntRank2A::Axy() const;                                INHERITED
  double IntRank2A::Azx() const;                                INHERITED
  double IntRank2A::Axz() const;                                INHERITED
  double IntRank2A::Azy() const;                                INHERITED
  double IntRank2A::Ayz() const;                                INHERITED


  double IntRank2A::Axx(double A, double B, double G) const;    INHERITED
  double IntRank2A::Ayy(double A, double B, double G) const;    INHERITED
  double IntRank2A::Azz(double A, double B, double G) const;    INHERITED
  double IntRank2A::Ayx(double A, double B, double G) const;    INHERITED
  double IntRank2A::Axy(double A, double B, double G) const;    INHERITED
  double IntRank2A::Azx(double A, double B, double G) const;    INHERITED
  double IntRank2A::Axz(double A, double B, double G) const;    INHERITED
  double IntRank2A::Azy(double A, double B, double G) const;    INHERITED
  double IntRank2A::Ayz(double A, double B, double G) const;    INHERITED

  double IntRank2A::Axx(const EAngles& EA) const;               INHERITED
  double IntRank2A::Ayy(const EAngles& EA) const;               INHERITED
  double IntRank2A::Azz(const EAngles& EA) const;               INHERITED
  double IntRank2A::Ayx(const EAngles& EA) const;               INHERITED
  double IntRank2A::Axy(const EAngles& EA) const;               INHERITED
  double IntRank2A::Azx(const EAngles& EA) const;               INHERITED
  double IntRank2A::Axz(const EAngles& EA) const;               INHERITED
  double IntRank2A::Azy(const EAngles& EA) const;               INHERITED
  double IntRank2A::Ayz(const EAngles& EA) const;               INHERITED

  row_vector IntRank2A::CartComps() const;                      INHERITED
  row_vector IntRank2A::CartComps(double theta, double phi=0) const;        */
 
// ----------------------------------------------------------------------------
//         Un-normalized Cartesian Components Of Hyperfine Tensor H
//                                                                 uv
// ----------------------------------------------------------------------------

/* These allow one to access the Cartesian elements of the full (unscaled)
   HF tensor at a specific orientation without rotating the entire tensor.  In
   these functions theta is the angle down from the lab. frame z-axis and phi
   the angle over from the lab. frame x-axis.  If no angles are specified, then
   the current orientation is used. The component output units will be Gauss.

                                    1/2
                            [ 6*PI ]
                      h   = | ---- |    * del   * A   + h
                       uv   [  5   ]         zz    uv    iso                 */


double IntHF::hxx() const { return AISO + DELZZ*RT6PIO5*Axx(); }
double IntHF::hyy() const { return AISO + DELZZ*RT6PIO5*Ayy(); }
double IntHF::hzz() const { return AISO + DELZZ*RT6PIO5*Azz(); }
double IntHF::hxy() const { return        DELZZ*RT6PIO5*Axy(); }
double IntHF::hyx() const { return        DELZZ*RT6PIO5*Ayx(); }
double IntHF::hxz() const { return        DELZZ*RT6PIO5*Axz(); }
double IntHF::hzx() const { return        DELZZ*RT6PIO5*Azx(); }
double IntHF::hyz() const { return        DELZZ*RT6PIO5*Ayz(); }
double IntHF::hzy() const { return        DELZZ*RT6PIO5*Azy(); }

double IntHF::hxx(double alpha, double beta, double gamma) const
                { return AISO + DELZZ*RT6PIO5*Axx(alpha,beta,gamma); }
double IntHF::hyy(double alpha, double beta, double gamma) const
                { return AISO + DELZZ*RT6PIO5*Ayy(alpha,beta,gamma); }
double IntHF::hzz(double alpha, double beta, double gamma) const
                { return AISO + DELZZ*RT6PIO5*Azz(alpha,beta,gamma); }
double IntHF::hyx(double alpha, double beta, double gamma) const
                { return        DELZZ*RT6PIO5*Ayx(alpha,beta,gamma); }
double IntHF::hxy(double alpha, double beta, double gamma) const
                { return        DELZZ*RT6PIO5*Axy(alpha,beta,gamma); }
double IntHF::hzx(double alpha, double beta, double gamma) const
                { return        DELZZ*RT6PIO5*Azx(alpha,beta,gamma); }
double IntHF::hzy(double alpha, double beta, double gamma) const
                { return        DELZZ*RT6PIO5*Azy(alpha,beta,gamma); }
double IntHF::hxz(double alpha, double beta, double gamma) const
                { return        DELZZ*RT6PIO5*Axz(alpha,beta,gamma); }
double IntHF::hyz(double alpha, double beta, double gamma) const
                { return        DELZZ*RT6PIO5*Ayz(alpha,beta,gamma); }

double IntHF::hxx(const EAngles& EA) const
                { return AISO + DELZZ*RT6PIO5*Axx(EA); }
double IntHF::hyy(const EAngles& EA) const
                { return AISO + DELZZ*RT6PIO5*Ayy(EA); }
double IntHF::hzz(const EAngles& EA) const
                { return AISO + DELZZ*RT6PIO5*Azz(EA); }
double IntHF::hyx(const EAngles& EA) const
                { return        DELZZ*RT6PIO5*Ayx(EA); }
double IntHF::hxy(const EAngles& EA) const
                { return        DELZZ*RT6PIO5*Axy(EA); }
double IntHF::hzx(const EAngles& EA) const
                { return        DELZZ*RT6PIO5*Azx(EA); }
double IntHF::hzy(const EAngles& EA) const
                { return        DELZZ*RT6PIO5*Azy(EA); }
double IntHF::hxz(const EAngles& EA) const
                { return        DELZZ*RT6PIO5*Axz(EA); }
double IntHF::hyz(const EAngles& EA) const
                { return        DELZZ*RT6PIO5*Ayz(EA); }

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
//        matrix IntRank2::CartMx(double scale=1.0) const
//        bool   IntRank2::Spherical( )
//        bool   IntRank2::Isotropic( )
//        matrix IntRank2::CartMx(bool scf)                                 */


/* This function also returns a scaled Cartesian spatial tensor. That tensor
   values may either be

        1.) GAMMA normalized Auv - Done With IntHF::CartMx()
        2.) Typical Suv values   - Done With IntHF::Smx(true);
        3.) Shown in lab frame   - Done With This Function

   For case 3.) the values are related to the GAMMA normalized (Auv) and
   typically presented values (guv) according to

                            1/2
                    [ 6*PI ]
             HF   = | ---- |    * del   * A   + Kdel    * HF
               uv   [  5   ]         zz    uv       u,v     iso

   where Kdel is a Kronecker delta function.                                 */

matrix IntHF::Amx() const
  {
  matrix HAmx = (DELZZ*RT6PIO5)*CartMx(1);		// Scaled A matrix
  matrix dmx(3,3,complex1,d_matrix_type);               // AIso matrix
  return HAmx + AISO*dmx;                               // The A matrix
  }

// ____________________________________________________________________________
// C                       SPIN TENSOR COMPONENT ACCESS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//                         Spin Quantum Value Access
//-----------------------------------------------------------------------------

/*
  double IntRank2T::Izval()          const;                     INHERITED
  double IntRank2T::Szval()          const;                     INHERITED
  int    IntRank2T::IV()             const;                     INHERITED
  int    IntRank2T::SV()             const;                     INHERITED
  int    IntRank2T::HS()             const;                     INHERITED
  double IntRank2T::qn(bool i=false) const;                     INHERITED    */

//-----------------------------------------------------------------------------
//                      Spin Tensor Component Access
//-----------------------------------------------------------------------------

/*     The commented member functions are INHERITED from class IntRank2T     */

// ----------------------------------------------------------------------------
//        Single Spin Or Spin Pair Hilbert Space Spin Tensor Components
// ----------------------------------------------------------------------------
/*
  matrix IntRank2T::T20()           const;      // Return T2,0  INHERITED
  matrix IntRank2T::T21()           const;      // Return T2,1  INHERITED
  matrix IntRank2T::T2m1()          const;      // Return T2,-1 INHERITED
  matrix IntRank2T::T22()           const;      // Return T2,2  INHERITED
  matrix IntRank2T::T2m2()          const;      // Return T2,-2 INHERITED
  matrix IntRank2T::T2m(int m)      const;      // Return T2,m  INHERITED
  matrix IntRank2T::Tcomp(int comp) const;      // comp=[-2,2]  INHERITED

// ----------------------------------------------------------------------------
//               Composite Hilbert Space Spin Tensor Components
// ----------------------------------------------------------------------------

  matrix IntRank2T::T20(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T21(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m1(     const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T22(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m2(     const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m(int m,const vector<int>& HSs,int i,int j=-1) const;  */

// ----------------------------------------------------------------------------
//      Additional Function To Handle "Perturbing" Hyperfine Interactions
// ----------------------------------------------------------------------------

/* These are the components of T20 that are invariant under rotations about the
   z-axis in the case that the spin pair is heteronuclear. Of course, for 
   hyperfine interacitons they are by nature heteronuclear. These funcitons
   facilitate construction of hypefine Hamiltonians in the rotationg frame 
   when the interaction very weak relative to the Zeeman interaction, and 
   that is most often the case. See T20wh and the function H0.               */


matrix IntHF::T20het() const { return T20wh; }
matrix IntHF::T20het(const vector<int>& HSs, int i, int j) const
                             { return blow_up(T20wh, HSs, i, j); }

// ____________________________________________________________________________
// D                  INTERACTION CONSTANT ACCESS FUNCTIONS
// ____________________________________________________________________________

/* The interaction constant is a scaling factor which is used to express all
   irreducible rank 2 interactions in a consistent fashion (class IntRank2)
   In fact, the interaction constant is maintained in the base class IntRank2
   as member variable _XI in radians/sec. This class simply relates the value
   of _XI to commonly used values that scale hyperfine interactions.

   We generate the interaction constant in radians/sec from the hyperfine
   delzz value (DELZZ) by direct conversion.  Since, for hyperfine coupling,
   DELZZ is maintained in Gauss, we convert from Gauss to rad/sec.


       HF          1 Hz        2.0*Pi rad
     xi   =  = ------------- * ---------- * DELZZ G = _XI (rad/sec)
       i       0.714567e-6 G     1 cycle

        // Input                HF      : Hyperfine interaction
        // Return               xi      : Interaction constant (rad/sec)     */

double IntHF::xi( ) const { return GAUSS2HZ*HZ2RAD*DELZZ; }

// ____________________________________________________________________________
// E                HYPEFINE INTERACTION AUXILIARY FUNCTIONS
// ____________________________________________________________________________

// None are defined as of yet

// ____________________________________________________________________________
// F                        PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//      Functions To Make A Parameter Set From A Hypefine Interaction
// ----------------------------------------------------------------------------

IntHF::operator ParameterSet( ) const
   { ParameterSet pset; pset += *this; return pset; }	

void operator+= (ParameterSet& pset, const IntHF &HF)
  { HF.PSetAdd(pset); }

void IntHF::PSetAdd(ParameterSet& pset, int idx, int pfx) const
  {
  string suffx;                                 // Parameter suffix
  if(idx != -1)                                 // Only use suffix if idx
    suffx = string("(")+Gdec(idx)+string(")");  // is NOT -1
  string prefx;                                 // Parameter prefix
  if(pfx != -1)                                 // Only use prefix if pfx
    prefx = string("[")+Gdec(pfx)+string("]");  // is NOT -1

  string pname  = prefx +string("Iqn") + suffx;  // Add I spin quantum number
  string pstate = string("Spin I Quantum Number");
  double pdatad = Izval();
  SinglePar par = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);

  pname  = prefx +string("Sqn") + suffx;	// Add S spin quantum number
  pstate = string("Spin S Quantum Number");
  pdatad = Szval();
  par = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);

  pname  = prefx + string("hxx") + suffx;       // Cartesian x component
  pstate = string("Cartesian X Component)");
  pdatad = hxx();
  par    = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);

  pname  = prefx + string("hyy") + suffx;       // Cartesian x component
  pstate = string("Cartesian Y Component)");
  pdatad = hyy();
  par    = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);

  pname  = prefx + string("hzz") + suffx;       // Cartesian x component
  pstate = string("Cartesian Z Component)");
  pdatad = hzz();
  par    = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);

  pname  = prefx + string("hEAngles") + suffx;  // Add orientation (degrees)
  pstate = string("Hyperfine Euler Angles (deg)");
  double a = _EAs.alpha() * RAD2DEG;
  double b = _EAs.beta()  * RAD2DEG;
  double g = _EAs.gamma() * RAD2DEG;
  coord EA(a, b, g);
  par = EA.param(pname, pstate);
  pset.push_back(par);
  } 

// ----------------------------------------------------------------------------
//  Functions To Output Hypefine Interaction To ASCII File (via Parameter Set)
// ----------------------------------------------------------------------------

/* These functions write parameters defining a hyperfine interaction into a
   specified ASCII file.

           Input                HF      : Hyperfine interaction
                                fn      : Output file name
                                ofstr   : Output file stream
                                idx     : Interaction index (default -1)
                                pfx     : Interaction 2nd indx (def -1)
                                wrn     : Warning level
           Output               none    : Dipolar interaction parameters
                                          written in parameter set format to
                                          file filename or ostream ofstr     */


int IntHF::write(const string &filename, int idx, int pfx, int warn) const
  {
  ofstream ofstr(filename.c_str());     // Open filename for input
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!write(ofstr, idx, pfx, w2))       // If file bad then exit
    {
    IHFerror(1, filename, 1);		// Filename problems
    if(warn>1) IHFfatal(20);		// Fatal error
    return 0;
    }
  ofstr.close();                        // Close it now
  return 1;
  }

int IntHF::write(ofstream& ofstr, int idx, int pfx, int warn) const
  {
  ParameterSet pset;                    // Declare a parameter set
  PSetAdd(pset, idx, pfx);              // Add in interaction parameters
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!pset.write(ofstr, w2))            // Use parameter set to write
    {                                   // out the interaction parameters
    if(warn)
      {
      IHFerror(52, 1);			// Problems writing to filestream
      if (warn>1) IHFfatal(53);		// Fatal error
      }
    return 0;
    }
  return 1;
  }

// ____________________________________________________________________________
// G                   HYPERFINE INTERACTION INPUT FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                     Functions To Read From An ASCII File 
// ----------------------------------------------------------------------------

/* These next two read functions utilize either two spin indices or a single
   interaction index. They'll try to read the hyperfine interaction form 
   parameters found either in an ASCII file or in a GAMMA parameter set. These
   functions do NOT allow for the prefix [#] because multiple hyperfine 
   interactions can be defined in the same file by switching either the spin 
   pair indices or the interaction index.  Multiple sets of interactions can be
   read using hyperfine interaction vectors (see class IntHFVec.)                 

	// Input		HF	: Hyperfine interaction
	//			filename: Output file name
	//			pset    : Parameter set
	//			idxI	: Index for spin I or interaction
	//			idxS	: Index for spin S (default -1)
        //                      warn    : Warning output label
        //                                 0 = no warnings
        //                                 1 = warnings
        //                                >1 = fatal warnings
	// Output		none 	: Interaction read from parameters
	//				  in file or pset                    */

bool IntHF::read(const string &filename, int idxI, int idxS, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename, warn?1:0))	// Read in pset from file
    {                                   // If we cannot read the file 
    if(warn)                            // then issue warnings as desired
      {
      IHFerror(1, filename, 1);		//      Filename problems
      if(warn > 1) IHFfatal(21);	//      Fatal error
      else         IHFerror(21);	//      or a warning issued
      }
    return false;                       //  Return that we failed!
    }
  return read(pset, idxI, idxS, warn);  // User overloaded function
  }


bool IntHF::read(const ParameterSet& pset, int idxI, int idxS, int warn)
  { 
  bool TF = setHFI(pset, idxI, idxS, warn?1:0);	// User overload to read
  if(!TF)                                       // If getHFI didn't handle
    {                                           // setting interaction proper
    if(warn)                                    // then we'll issue some
      {                                         // warnings if desired
      string val;
      int eidx = 44;
      if(idxI == -1)      { val = string(" None"); }
      else if(idxS == -1) { val = Gdec(idxI); }
      else
        { val = Gdec(idxI) + string(" & ")
              + Gdec(idxS);
          eidx = 45;
        }
                   IHFerror(23, 1);		//   Insufficient parameters
      if(warn > 1) IHFfatal(eidx, val);	// Fatal error
       else        IHFerror(eidx, val);	// or a warning issued
      }
    return false; 				// Return that we failed
    }
  return TF;
  }

// ____________________________________________________________________________
// H                    HYPERFINE HAMILTONIAN FUNCTIONS
// ____________________________________________________________________________

/* This section returns irreducible rank 2 hyperfine interaction Hamiltonians.
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
//               First Order Hyperfine Interaction Hamiltonian
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
//  These Are For When The Hyperfine Interaction Is A Perturbation To Zeeman
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

       [1]                      HF  HF                       HF     (0)
      H  (alpha,beta,gamma) = Xi * A   (alpha,beta,gamma) * T    = H
       HF                           2,0                      2,0    HF

           Input                HF      : Hyperfine interaction
                                alpha   : Euler angle (radians)
                                beta    : Euler angle (radians)
                                gamma   : Euler angle (radians)
                                EA      : Euler angles (radians)
                                HSs     : Array of spin Hilbert spaces
                                i       : First spin index (in HSs)
                                j       : Seond spin index (in HSs)
                                wh      : Flag for weak heteronuclear
           Output               H0      : The secular part of the hyperfine
                                          Hamiltonian (default basis, Hz)
           Note                         : Also called the 1st order Hyperfine
                                          interaction (perturbation theory)
           Note                         : Rotationally invariant about z
                                          when the wh flag is set true
           Note                         : No HSs: return in the spin Hilbert
                                          space of dimension (2I+1)[*(2S+1)]
                                          With HSs: return in composite spin
                                          Hilbert space                      */

matrix IntHF::H0(bool wh) const
  { return wh?(_XI*A20())*T20wh:IntRank2::H0(); }

matrix IntHF::H0(double A, double B, double G, bool wh) const
  { return wh?(_XI*A20(A,B,G))*T20wh:IntRank2::H0(A,B,G); }

matrix IntHF::H0(const EAngles& EA, bool wh) const
  { return wh?(_XI*A20(EA))*T20wh:IntRank2::H0(EA); }


matrix IntHF::H0(const vector<int>& HSs, int i, int j, bool wh) const
  { return wh?(_XI*A20())*T20het(HSs,i,j):IntRank2::H0(HSs,i,j); }

matrix IntHF::H0(const vector<int>& HSs, int i, int j,
                  double A, double B, double G, bool wh) const
  { return wh?(_XI*A20(A,B,G))*T20het(HSs,i,j):IntRank2::H0(HSs,i,j,A,B,G); }

matrix IntHF::H0(const vector<int>& HSs, int i, int j,
                             const EAngles& EA, bool wh) const
  { return wh?(_XI*A20(EA))*T20het(HSs,i,j):IntRank2::H0(HSs,i,j,EA); }


// ----------------------------------------------------------------------------
//                 Full Hyperfine Interaction Hamiltonians
//       These Have No Approximations & Are Independent of Field Strength
//     (They Are Correct in The Laboratory Frame, NOT in a Rotating Frame)
// ----------------------------------------------------------------------------

/* Here we make no approximations. As such, these Hamiltonians are not those
   used in the rotating frame. Here we simply sum over the 5 spatial and spin
   tensor components to produce the Hamiltonian.

                           -2,2
                            ---   HF      m    HF                         HF
  H (alpha, beta, gamma) =  \   Xi  * (-1)  * A   (alpha, beta, gamma) * T
   HF                       /                  2,m                        2,-m
                            ---
                             m                                               */

// matrix IntRank2::H( ) const                                        INHERITED
// matrix IntRank2::H(double alpha, double beta, double gamma) const  INHERITED
// matrix IntRank2::H(const EAngles& EA) const                        INHERITED

// matrix IntDip::H(const vector<int>& HSs, int i, int j) const       INHERITED
// matrix IntDip::H(const vector<int>& HSs, int i, int j,             INHERITED
//                     double alpha, double beta, double gamma) const
// matrix IntDip::H(const vector<int>& HSs, int i, int j,             INHERITED
//                                           const EAngles& EA) const

// ____________________________________________________________________________
// J                 HYPERFINE INTERACTION OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//  Functions That Generate String Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

/*  Function   Arguments                     Output
    ========   =========   ===================================================
    TStrings       m       String array for the mth component of spin tensor T
    GAStrings              String array for various interaction values

            HFStrings()                              TStrings(m)

     Hyperfine Coupling:   xxxx.xx Gauss          [ x.x, x.x, x.x]
     Hyperfine Anisotropy: xxxx.xx Gauss   T    = [ x.x, x.x, x.x]
     Hyperfine Asymmetry:  xxxx.xx Gauss    2,m   [ x.x, x.x, x.x]
     Down From PAS z-Axis:    x.xx Deg.           
     Over From PAS x-Axis:    x.xx Deg.
     Electron I Value:        I               m = [0,4] => {0,1,-1,2,-2}
     Nucleus  I Value:        S                                              */
 
vector<string> IntHF::CartAStrings(const string& CSForm) const
  {
  int nstr = 6;
  vector<string> Cartstrings(nstr);     // Vector of nstr strings
  Cartstrings[0] = "[a  , a  , a  ]";
  Cartstrings[1] = "[ xx   xy   xz]"
                 + string("   ")
                 + "["  + Gform(CSForm.c_str(), hxx())
                 + ", " + Gform(CSForm.c_str(), hxy())
                 + ", " + Gform(CSForm.c_str(), hxz()) + "]";
  Cartstrings[2] = "[a  , a  , a  ]"
                 + string(" = ")
                 + "["  + Gform(CSForm.c_str(), hyx())
                 + ", " + Gform(CSForm.c_str(), hyy())
                 + ", " + Gform(CSForm.c_str(), hyz()) + "]";
  Cartstrings[3] = "[ yx   yy   yz]"
                 + string("   ")
                 + "["  + Gform(CSForm.c_str(), hzx())
                 + ", " + Gform(CSForm.c_str(), hzy())
                 + ", " + Gform(CSForm.c_str(), hzz()) + "]";
  Cartstrings[4] = "[a  , a  , a  ]";
  Cartstrings[5] = "[ zx   zy   zz]";
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


vector<string> IntHF::InfoStrings() const
  {
  vector<string> SphStrings(5);				// Array of info strings
  int k=0;
  SphStrings[k++] = StringIS();
  SphStrings[k++] = string("Hyperfine Coupling:")
                  + Gform("%10.2f", AISO) + string(" G   ");
  SphStrings[k++] = string("Hyperfine Anisotropy:")
                  + Gform("%8.2f", DELZZ) + string(" G   ");
  SphStrings[k++] = string("Hyperfine Asymmetry:    ")
                  + Gform("%10.7f", ETA);
  SphStrings[k++] = XiString();
  return SphStrings;
  }

vector<string> IntHF::SphAStrings() const
  {
  vector<string> Sphstrings(4);
  int k=0;
  Sphstrings[k++] = StringIS();
  Sphstrings[k++] = string("Hyperfine Coupling:")
                  + Gform("%10.2f", AISO) + string(" Gauss");
  Sphstrings[k++] = string("Hyperfine Anisotropy:")
                  + Gform("%8.2f", DELZZ) + string(" Gauss");
  Sphstrings[k++] = string("Hyperfine Asymmetry:    ")
                  + Gform("%10.7f", ETA) + string(" ");
  return Sphstrings;
  }

//-----------------------------------------------------------------------------
//     Functions That Generate Ouput Of The Rank 2 Hyperfine Interaction
//-----------------------------------------------------------------------------

/* These functions will output information concerning the Hyperfine interaction
   to any output stream.

           Input                HF      : Hyperfine interaction (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Basic Parameters
                                           !0 - Full output
                                nrm     : Flag if GAMMA normalized output
           Output               none    : Hyperfine interaction placed
                                          into the output stream            */

ostream& IntHF::print(ostream& ostr, int fflag) const
  {
  if(Izval() < 0.5)
    {
    string hdr("Empty Hyperfine Interaction");
    string spacer = string(40-hdr.length()/2, ' ');
    ostr << "\n\n" << spacer << hdr << "\n";
    return ostr;
    }

/*  	    Output Some Printed Lines Which Will Look Like The Following 

  			       Hyperfine Interaction
  
   Hyperfine Coupling:    xxxxx.xx Gauss     [h  , h  , h  ] 
   Nucleus I Value:           I              [ xx   xy   xz]   [ x.x, x.x, x.x]
   Hyperfine Anisotropy:  xxxxx.xx Gauss     [h  , h  , h  ] = [ x.x, x.x, x.x]
   Hyperfine Asymmetry:       x.xx           [ yx   yy   yz]   [ x.x, x.x, x.x]
   Down From PAS z-Axis:    xxx.xx Degrees   [h  , h  , h  ]
   Over From PAS x-Axis:    xxx.xx Degrees   [ zx   zy   zz]                 */

  vector<string> Istrs = InfoStrings();			// Information strings
  vector<string> Astrs = CartAStrings("%6.3f");		// Cartesian A strings
  string hdr = "Hyperfine Interaction";			// Use this header
  string Spacer((40-hdr.length()/2), ' ');		// Spacer to center 
  ostr << "\n\n" << Spacer << hdr << "\n";		// Output centered hdr
  IntRank2A::print(ostr, Astrs, Istrs);
  if(!fflag) return ostr;

/*	    Now Output The Spherical Components, Spatial And Spin
       These Will Be Printed Lines Which Will Look Like The Following 
  		      (Repeated For All 5 m Values)
  
   						[ x.x, x.x, x.x]
  		A    = x.xxx		T     = [ x.x, x.x, x.x]
  		 2,m			 2,m	[ x.x, x.x, x.x]
  						[ x.x, x.x, x.x]             */

  printAT(ostr, HF);
  ostr << "\n\n";
  return ostr;
  }


ostream& operator<< (ostream& out, IntHF& HF) { return HF.print(out); }
ostream& IntHF::printSpherical(ostream& ostr)
         
        // Input                HF	: Hyperfine spatial tensor (this)
        //                      ostr	: Output stream
        // Output               none    : Rank 2 spatial tensor parameters
        //                                sent to the output stream
	// Note				: Uses base class virtual overload

  {
  ostr << "\t\n\t         Hyperfine Spatial Tensor";
  IntRank2A::printSpherical(ostr, 0);
  return ostr;
  }

 
/*
ostream& IntHF::printCartesian(ostream& ostr)

        // Input                HF	: Hyperfine spatial tensor (this)
        //                      ostr	: Output stream
        // Output               none    : Rank 2 spatial tensor parameters
        //                                sent to the output stream
	// Note				: Uses base class virtual overload

  {
  ostr << "\t\n\t         Hyperfine Spatial Tensor";
  IntRank2A::printCartesian(ostr, 0);
  return ostr;
  }

ostream& IntHF::printCartesian(ostream& ostr, double theta, double phi)
         
        // Input                HF	: Hyperfine spatial tensor (this)
        //                      ostr	: Output stream
        //                      theta   : Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               none    : Rank 2 spatial tensor parameters
        //                                sent to the output stream
	// Note				: Uses base class virtual overload

// Axx = (A22+A2m2)/2 - A20/sqrt(6)      Axy = -i*(A22-A21)/2 = Ayx
// Ayy = -(A22+A2m2)/2 - A20/sqrt(6)     Axz = -(A21-A2m1)/2  = Azx
// Azz = sqrt(2/3)*A20                   Ayz = i*(A21-A2m1)/2 = Azy

  {
  ostr << "\t\n\t         Hyperfine Spatial Tensor";
  IntRank2A::printCartesian(ostr, theta, phi, 0);
  return ostr;
  }
*/
 



















// ____________________________________________________________________________
// D                      HYPERFINE FREQUENCY FUNCTIONS
// ____________________________________________________________________________
 
// This frequency will be the splitting between the 2*I transitions contained
// in a Hyperfine Hamiltonian.  This will be when the hyperfine coupling is weak
// relative to the Zeeman interaction (high field approximation, first order
// terms only) if the tensor is oriented in it's principal axes (PAS).  If the
// tensor isn't aligned in it's PAS, the splitting will vary with orientation 
// according to

//                     1           2                2                        
//    W (theta,phi)  = - W  [ 3*cos (theta) - eta*sin (theta)cos(2*phi) ]
//     D               2  D

// where theta is the angle down from the PAS z axis and phi the angle over
// from the PAS x axis.  Alternatively, if the Zeeman terms aint much stronger
// the splittings won't be equally spaced at all (second order terms).
 
// Also, keep in mind that the Euler angles {phi,theta,gamma} which are kept
// with the tensor are used to relate the tensor PAS to some (common) set of
// coordinate axes.  They are not the phi and theta used in the above
// formula (unless you which the tensor aligned in the common axis system)


//double IntHF::wD( ) const { return DELZZ; }

        // Input                HF	: Hyperfine interaction
        // Return               wD      : Hyperfine frequency (Hz)

  


//void IntHF::wHF(double W)

        // Input                HF	: Hyperfine interaction
	//			W       : Hyperfine frequency (Hz)
        // Return               void	: The Hyperfine frequency 
	//			          of the interaction is set
	// Note				: The interaction I value
	//				  must be set prior to this 

//  {
//  DELZZ = W;
// sosi
//  }




// ____________________________________________________________________________
// L                    HYPERFINE HAMILTONIAN FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//               First Order Hyperfine Interaction Hamiltonian
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
//  These Are For When The Hyperfine Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------
 
/*
matrix IntHF::H0( ) const
 
	// Input		HF	: Hyperfine interaction
        // Output               H0	: The secular part of the hyperfine
        //                                Hamiltonian (default basis, Hz)
        //                                hyperfine interaction
        // Note				: Also called the 1st order hyperfine
        //                                interaction (perturbation theory)
        // Note                      	: Rotationally invariant about z
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1
 
//  The secular part of the hyperfine Hamiltonian is that which returns
//  only those components which commute with z axis rotations.  Here we have
 
//              [1]               D   D               D      (0)
//             H  (theta,phi) = Xi * A (theta,phi) * T    = H
//              D                     2,0             2,0    D
 
  { return (xi()*Acomp(0))*T20(); }


matrix IntHF::H0(double theta, double phi) const
 
	// Input		HF	: Hyperfine interaction
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               H0	: The secular part of the hyperfine
        //                                Hamiltonian (default basis, Hz)
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note				: Also called the 1st order hyperfine
        //                                interaction (perturbation theory)
        // Note                      	: Rotationally invariant about z
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1
 
//  The secular part of the hyperfine Hamiltonian is that which returns
//  only those components which commute with z axis rotations.  Here we have
 
//              [1]               HF   HF               HF     (0)
//             H  (theta,phi) = Xi  * A  (theta,phi) * T    = H
//              HF                     2,0              2,0    HF
 
  { return (xi()*A20(theta, phi))*T20(); }
*/

 
// ----------------------------------------------------------------------------
//              Second Order Hyperfine Interaction Hamiltonians
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
//  These Are For When The Hyperfine Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------

 
// ----------------------------------------------------------------------------
//      Summed First & Second Order Hyperfine Interaction Hamiltonians
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
//  These Are For When The Hyperfine Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------
 

//matrix IntHF::Hw(double Om) const

        // Input                sys  : Spin system
        //                      wD   : Hyperfine frequency (Hz)
	//			i    : Spin index
        // Output               HDw  : The secular part of the hyperfine
        //                             Hamiltonian (default basis, Hz)
	//			       for the spin i
	// Note			     : No asymmetry is considered here
	// Note			     : This is the sum of the 1st & 2nd order
	//			       hyperfine interactions (pert. theory)
	// Note			     : This is rotationally invariant about z

//  This function returns only the secular part of the second order hyperfine
//  Hamiltonian (from perturbation theory).  Note that this still assumes that
//  the hyperfine interaction is a perturbation to to the Zeeman interaction.

//                                (0)    (1)      [1]    [2]
//                       H    =  H    + H     =  H    + H
//                        D       D      D        D      D

//  { return H0() + H1(Om); }

 
//matrix IntHF::Hw(double Om, double theta, double phi) const
 
        // Input                Om	: Field Strength (Larmor in Hz)
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HHF	: The 2nd order secular part of the
        //                                hyperfine Hamiltonian (default
        //                                basis, Hz) for the interaction
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note				: Also called the 1st order hyperfine
        //                                interaction (perturbation theory)
        // Note                      	: Rotationally invariant about z
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1

//  { return H0(theta, phi) + H1(Om, theta, phi); }


// ----------------------------------------------------------------------------
//                 Full Hyperfine Interaction Hamiltonians
//       These Have No Approximations & Are Independent of Field Strength
//     (They Are Correct in The Laboratory Frame, NOT in a Rotating Frame) 
// ----------------------------------------------------------------------------

/*                              HF      m    HF                HF
              H (theta,phi) = Xi  * (-1)  * A   (theta,phi) * T
               HF                            2,m               2,-m          */
 
/*
matrix IntHF::H( ) const
 
	// Input		HF	: Hyperfine interaction
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1

  {
  matrix Hmx = Acomp(0)*T20();			// First set the m=0 terms
  if(norm(Acomp(1)))				// Then add in m=+/-1 terms
    Hmx -= (Acomp(1)*T2m1()+Acomp(2)*T21());	// if they are present
  if(norm(Acomp(3)))				// Then add in m=+/-2 terms
    Hmx += (Acomp(3)*T2m2()+Acomp(4)*T22());	// if they are present
  return xi()*Hmx;
  }


// ____________________________________________________________________________
// N                 HYPERFINE HAMILTONIAN FRIEND FUNCTIONS
// ____________________________________________________________________________

// These Don't Depend on Class IntHF Variables At All.  The Functions Just
// Mimic The Class' Hamiltonian Generation.


matrix HF0(double qn, double wDo, double eta, double theta, double phi)
 
	// Input		qn	: Quantum number (1, 1.5, 2.5,...)
	//			wDo     : PAS Hyperfine frequency
	//			eta     : Hyperfine asymmetry [0,1]
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               H0	: The secular part of the hyperfine
        //                                Hamiltonian (default basis, Hz)
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note				: Also called the 1st order hyperfine
        //                                interaction (perturbation theory)
        // Note                      	: Rotationally invariant about z
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1
 
//  The secular part of the hyperfine Hamiltonian is that which returns
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

*/

/*
matrix HF1(double Om, double qn, double wDo, double eta,
                                                      double theta, double phi)
 
        // Input                Om	: Field Strength (Larmor in Hz)
	// 			qn	: Quantum number (1, 1.5, 2.5,...)
	//			wDo     : PAS Hyperfine frequency
	//			eta     : Hyperfine asymmetry [0,1]
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HD1	: The 2nd order secular part of the
        //                                hyperfine Hamiltonian (default
        //                                basis, Hz) for the interaction
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note				: Also called the 1st order hyperfine
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
 
//  in accordance with the article by P.P. Man "Hyperfine Interactions" in
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
*/
 





// ____________________________________________________________________________
// P                          STANDARD I/O FUNCTIONS
// ____________________________________________________________________________

/*
string IntHF::ask_read(int argc, char* argv[], int argn, int idxI,int idxS)

        // Input                HF	: Hyperfine interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
	//			idxI	: Index for spin I
	//			idxS	: Index for spin S
        // Output               String  : The parameter argn of array argc is
        //                                used to supply a filename from which
        //                                the hyperfine interaction is read
        //                                If argument argn is not in argv, the
        //                                user is asked to supply a filename
        //                                The set filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing known IntHF parameters
        // Note                         : The interaction D is modifed (filled)

  {
  string filename;				// Name of parameter file
  query_parameter(argc, argv, argn,		// Get filename from command
    "\n\tHyperfine Interaction filename? ",	// Or ask for it
                                      filename);
  read(filename, idxI, idxS); 			// Read system from filename
  return filename;
  }



string IntHF::ask_read(int argc, char* argv[], int argn,
                                               double Iqn,double Sqn, int idx)

        // Input                HF	: Hyperfine interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
	//			Iqn	: Quantum I value (0.5, 1, ...)
	//			Sqn	: Quantum S value (0.5, 1, ...)
	//			idx	: Interaction index (default -1->none)
        // Output               String  : The parameter argn of array argc is
        //                                used to supply a filename from which
        //                                the hyperfine interaction is read
        //                                If argument argn is not in argv, the
        //                                user is asked to supply a filename
        //                                The set filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing known IntHF parameters
        // Note                         : The interaction D is modifed (filled)

  {
  string filename;				// Name of parameter file
  query_parameter(argc, argv, argn,		// Get filename from command
    "\n\tHyperfine Interaction filename? ",	// Or ask for it
                                      filename);
  read(filename, Iqn, Sqn, idx); 		// Read system from filename
  return filename;
  }


void IntHF::ask(int argc, char* argv[], int& qn, double& DI, double& DS,
          double& Dnqcc, double& Deta, double& Dtheta, double& Dphi, int Dflag)

        // Input                HF	: Hyperfine interaction (this)
	//			argc    : Number of arguments
	//			argv    : Array of arguments
	//			qn      : Query value
	//			DI	: Spin quantum number
	//			DS	: Spin quantum number
	//			Dnqcc   : Hyperfine. coupling constant (Hz)
	//			Deta    : Hyperfine. asymmetry
	//			Dtheta  : Hyperfine. orientation angle
	//			Dphi	: Hyperfine. orientation angle
        //                      Dflag   : Flag is DCC or wD requested
        // Output               none    : The values of qn, I, Dnqcc, Deta,
	//				  Dtheta, and Dphi are set herein
	// Note				: This is INTERACTIVE!

  {
  query_parameter(argc, argv, qn++,			// Read in the I value
      "\n\tI Spin Quantum Value (0.5, 1, 1.5, ..)? ", DI);
  query_parameter(argc, argv, qn++,			// Read in the I value
      "\n\tS Spin Quantum Value (0.5, 1, 1.5, ..)? ", DS);
  if(Dflag)
    {
    query_parameter(argc, argv, qn++,			// Read the frequency
       "\n\tHyperfine Frequency(kHz)? ", Dnqcc);
    Dnqcc *= 2.0*I*(2.0*I-1.0)/3.0;			// Set quad. coupling
    }
  else
    {
    query_parameter(argc, argv, qn++,			// Read in the coupling
       "\n\tHyperfine Coupling (kHz)? ", Dnqcc);
    }
  Dnqcc *= 1.e3;					// Put this in Hz
  query_parameter(argc, argv, qn++,			// Read in the coupling
       "\n\tHyperfine Asymmetry [0, 1]? ", Deta);
  query_parameter(argc, argv, qn++,			// Read in theta angle
  "\n\tOrientation Down From PAS z [1, 180]? ", Dtheta);
  query_parameter(argc, argv, qn++,			// Read in phi angle
   "\n\tOrientation Over From PAS x [0, 360]? ", Dphi);
  }


void IntHF::askset(int argc, char* argv[], int& qn, int Dflag)

        // Input                HF	: Hyperfine interaction (this)
	//			argc    : Number of arguments
	//			argv    : Array of arguments
	//			qn      : Query value
        //                      Dflag   : Flag is DCC or wD requested
        // Output               none    : D is set interactively
	// Note				: This is INTERACTIVE!

  {
  double DI, DS, Dnqcc, Deta, Dtheta, Dphi;
  ask(argc,argv,qn,DI,DS, Dnqcc,Deta,Dtheta, Dphi, Dflag);	// Use the ask function
  *(this) = IntHF(DI,DS, Dnqcc,Dtheta,Dphi,Deta);	// Use assignment
  }

 
void IntHF::askset(int Dflag)
 
        // Input                HF	: Hyperfine interaction (this)
        //                      Dflag   : Flag is DCC or wD requested
        // Output               none    : D is set interactively
        // Note                         : This is INTERACTIVE!
 
  {
  int qn = 1000;
  int argc = 1;
  char* argv[1];
  double I, S, dcc, eta, theta, phi;
  ask(argc,argv,qn,I,S,dcc,eta,theta, phi,Dflag);	// Use the ask function
  *(this) = IntHF(I,S,dcc,theta,phi,eta);         // Use assignment
  }

*/

#endif							// IntHF.cc
