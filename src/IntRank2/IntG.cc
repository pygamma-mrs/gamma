/* IntG.cc ******************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Electron G Interaction 	                 Implementation		**
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
** A variable of type IntG represents an electron G interaction defined **
** for a particular electron spin (or electron in a spin system).  The	**
** interaction is maintained in the form of spatial and	spin tensors	**
** which are stored in irreducible spherical format. In addition their  **
** is an isotropic value and an interaction constant that sets the 	**
** overall scaling on the irreducible rank 2 part.			**
**                                                                      **
** The following quantities are inherited/maintained in the base class	**
** IntRank2A (via IntRank2, normalized GAMMA scaling):			**
**                                                                      **
**  ETA     - Rank 2 spatial tensor asymmetry value, range [0, 1].	**
**                                                                      **
** Similarly, the following quantities are inherited/maintained in the	**
** base class IntRank2T (via IntRank2):					**
**                                                                      **
**  Ival    - Hilbert space of the spin involved (typically 2 for e-)	**
**  T{0-2}  _ Spin tensor components T20, T21, and T22			**
**  Tm{1,2} _ Spin tensor components T2-1 and T2-2			**
**                                                                      **
** The generic irreducible rank 2 interaction class IntRank2 will keep	**
** track of:								**
**                                                                      **
**  I       - Spin quantum number of electron involved                  **
**  _XI     - Interaction strength (radians/sec)			**
**                                                                      **
** Lastly, this class will itself maintain:				**
**                                                                      **
**  GISO    - Isotropic electron g-factor (unitless).			**
**  DELZZ   - Rank 2 spatial tensor delz value (unitless).		**
**  BoFIELD - Base external field strength (Gauss)			**
**                                                                      **
** The tensor will internally be stored in irreducible spherical form   **
** in the tensor PAS.							**
**                                                                      **
** In addition to the functions contained in the base spatial tensor    **
** class, this IntG provides functions for building up the tensor and	**
** accessing the tensor elements from a "g-factor" standpoint.  This	**
** includes functions for reorienting the tensor, making Hamiltonians,	**
** and obtaining electron resonance conditions.				**
**                                                                      **
** The following defintions are used herein (Auv normlized by delzz):	**
**                                                                      **
** 1.) PAS: |Azz| >= |Ayy| >= |Axx|                                     **
**                                                                      **
** 2.) PAS: Azz=2C, eta=(Axx-Ayy)/Azz, Axx=C(eta-1), Ayy=-C(1+eta)      **
**                                                                      **
** Individual spin tensors are stored in a linked list, handled by the  **
** class IntSTLList.  This prevents recalculation of the 5 arrays that  **
** are associated with such tensors when the spin(s) involved share the **
** same Iz value(s).  However, the arrays are still copied into new     **
** (equivalent) spin tensors, but that is done by reference in the      **
** GAMMA matrix class.                                                  **
**                                                                      **
*************************************************************************/

#ifndef   IntG_cc_			// Is file already included?
#  define IntG_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <IntRank2/IntG.h>		// Include interface definition
#include <IntRank2/IntRank2.h>		// Include base class
#include <IntRank2/IntRank2A.h>		// Include base class base class
#include <IntRank2/IntRank2T.h>		// Include base class base class
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/Gutils.h>		// Include parameter queries
#include <HSLib/SpinOpSng.h>		// Include single spin operators
#include <Level1/coord.h>		// Include coordinates
#include <Matrix/row_vector.h>		// Include row_vectors
#include <Matrix/matrix.h>		// Include matrices
#include <HSLib/SpinSys.h>		// Include base spin systems
#include <Basics/StringCut.h>
#include <string>
#include <list>				// Know STL lists
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
// i             CLASS G FACTOR INTERACTION ERROR HANDLING
// ____________________________________________________________________________
 

/*       Input                G       : Electron G interaction (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void IntG::IGerror(int eidx, int noret) const
  {
  string hdr("Electron G Interaction");
  switch(eidx)
    {
     case 0: GAMMAerror(hdr,"Program Aborting.....",        noret); break;// (0)
     case 1: GAMMAerror(hdr,"Construction From Rank!=2.",   noret); break;// (1)
     case 2: GAMMAerror(hdr,"Problems During Construction.",noret); break;// (2)
     case 3: GAMMAerror(hdr,"Problems During Assignment.",  noret); break;// (3)
     case 8: GAMMAerror(hdr,"Theta (z Down) Beyond [0,180]",noret); break;// (8)
     case 9: GAMMAerror(hdr,"Phi (x Over) Outside [0, 360]",noret); break;// (9)
     case 10:GAMMAerror(hdr,"Asymmetry (eta) Beyond [0, 1]",noret); break;// (10)
     case 12:GAMMAerror(hdr,"Set Asymmetry On Zero Tensor", noret); break;// (12)
     case 13:GAMMAerror(hdr,"Cannot Construct From File",   noret); break;// (13)
     case 14:GAMMAerror(hdr,"Cannot Set Shift Anisotropy",  noret); break;// (14)
     case 15:GAMMAerror(hdr,"Setting Asymmetry to Zero",    noret); break;// (15)
     case 16:GAMMAerror(hdr,"Setting Default I=1/2 Value",  noret); break;// (16)
     case 18:GAMMAerror(hdr,"Cannot Alter Interaction",     noret); break;// (18)
     case 19:GAMMAerror(hdr,"Can't Write To Parameter File",noret); break;// (19)
     case 20:GAMMAerror(hdr,"Can't Write To Parameter File",noret); break;// (20)
     case 21:GAMMAerror(hdr,"Cant Read From Parameter File",noret); break;// (21)
     case 22:GAMMAerror(hdr,"Insufficient File Parameters", noret); break;// (22)
     case 23:GAMMAerror(hdr,"Insufficient PSet Parameters", noret); break;// (23)
     case 24:GAMMAerror(hdr,"Insufficient Spatial Params",  noret); break;// (24)
     case 25:GAMMAerror(hdr,"Spin Is Not An Electron",      noret); break;// (25)
     case 30:GAMMAerror(hdr,"Quantum # < 1/2 Specified?",   noret); break;// (30)
     case 31:GAMMAerror(hdr,"Quantum # Not Multiple Of 1/2",noret); break;// (31)
     case 42:GAMMAerror(hdr,"G_T Use, Parameter Not Type 4",noret); break;// (42)
     case 43:GAMMAerror(hdr,"G_T Use, Must Be Rank 2",      noret); break;// (43)
     case 44:GAMMAerror(hdr,"G_T Use, Ignore Isotropic Val",noret); break;// (44)
     case 48:GAMMAerror(hdr,"Improper Spin Quantum Number", noret); break;// (48)
     case 49:GAMMAerror(hdr,"Improper Spin Designation",    noret); break;// (49)
     case 50:GAMMAerror(hdr,"Invalid Component, m=[-2,2]",  noret); break;// (50)
     case 51:GAMMAerror(hdr,"Cannot Write To File",   	    noret); break;// (51)
     case 52:GAMMAerror(hdr,"Cannot Write To FileStream",   noret); break;// (52)
     case 53:GAMMAerror(hdr,"Cannot Output Parameters",     noret); break;// (53)
     case 60:GAMMAerror(hdr,"Cannot Set From Parameters",   noret); break;// (60)
     default:GAMMAerror(hdr,"Unknown Error",                noret); break;
    }
  }


volatile void IntG::IGfatal(int eidx) const
  {
  IGerror(eidx,1);				// Output error message
  if(eidx) IGerror(0, 1);			// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Function pname
      5     Please Use Class Member Function pname
   default  Unknown Error - pname                                           */

void IntG::IGerror(int eidx, const string& pname, int noret) const 
  {
  string hdr("Electron G Interaction");
  string msg;
  switch(eidx)
    {
    case 101:                                                   // (101)
      msg = string("Can't Find Interaction Parameters For ")
          + pname;
      GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;
    }
  }

volatile void IntG::IGfatal(int eidx, const string& pname) const
  {
  IGerror(eidx, pname, eidx);                   // Output error message
  if(eidx) IGerror(0);                          // State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________
// ii                  G FACTOR INTERACTION SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions set up specific aspects of an electron G interaction.
   As the various interaction parameters interact, the functions MUST be
   private because their misuse could produce an inconsistent interaction.

   The goal here is quite simple. We must determine the following set of values
   for each G interaction: { Iqn, g, gA, eta, alpha, beta, gamma }
   Complexity arises because we allow that a variety of parameters may be used
   to define these values. Additionally we allow that some parameters may be
   zero (e.g. eta) and that defaults may automatically be used. But, the end
   result remains the same, we seek the set of values that define an electron
   G interaction in GAMMA.                                                   */

// ----------------------------------------------------------------------------
//                        Complete Electron G Interaction
// ----------------------------------------------------------------------------

/* These functions will try and get all the parameters required to define an
   electron G interaction: { Iqn, g, gA, eta, alpha, beta, gamma }. We will
   also see if we can figure out the external magnetic field strength (which is
   akin to knowing the base EPR frequency). To support G interactions in spin
   systems, which typically know spin isotope types and their spin quantum
   numbers, there is a separation between determining the spin quantum
   number and the rest of the interaction.                                   */

bool IntG::getGI(const ParameterSet& pset,
         double& Iqn, double& g, double& gA, double& eta, EAngles& EA,
                                         double& Bo, int idx, bool warn) const
  {
//                      Get The Spin Quantum Number
//                    (And Perhaps The Isotope Type)

  string pb("Iqn");                             // Parameter base name
  string II;                                    // String for isotope name
  Isotope ISI;                                  // Isotope for spin
  bool TFI = false;                             // Flag if we know isotope
  if(getIso(pset,II,idx,0))                     // 1. Try for isotope name
    {                                           //    If successful, check it
    if(!SpinCheck(II)) return false;            //    Insure valid isotope
    ISI = Isotope(II);                          //    An isotope of type II
    if(!SpinCheck(ISI,true,warn)) return false;	//    Disallow nuclear spin
    TFI = true;                                 //    We know the isotope
    Iqn = ISI.qn();                             //    Know spin I quantum #
    }
  else
    {
    if(getIqn(pset,pb,Iqn,idx,false))           // 2. Try for spin quant. #
      { if(!SpinCheck(Iqn)) return false; }     //    Insure valid  qn
    else Iqn = 0.5;                             // 3. Use default qn of 1/2
    }

//                      Get The External Field Strength
//                               Not Mandatory!

    getField(pset, Bo, idx, 0);

// Try To Directly Read Electron G Cartesian Spatial Tensor Components
//        Interaction Via { Iqn, Gxx, Gxy, Gxz, Gyy, Gyz, Gzz }

//  1.) Guv values set the g-factor, anisotropy & asymmetry: {GISO, DELZ, ETA}
//  2.) If any Guv off-diagonals (u != v) set, G array sets orientation also
//  3.) If no Guv off-diagonals, orientation set with specified Euler angles
//  4.) The input Guv values are taken to be unitless
//  5.) We don't mind that no orientation is set, default is PAS
//  6.) If used, Euler angles by either a set or 3 individual angles
//  7.) Cannot mix parameters with g and G for this read

  coord GiGzGe;                                 // For Giso, Gzz, Geta
  if(getACart(pset,"G",GiGzGe,EA,idx,-1,0))     // Try & use Cart. components
    {
    g   = GiGzGe.x();				//  Set the g factor
    gA  = GiGzGe.y();				//  Set the delzz value
    eta = GiGzGe.z();                           //  Set the eta value
    if(EA == EAzero)				//  Allow use of g with G
      getOrientation(pset,"g",EA,idx,-1,false); //  in parameter names
    return true;
    }
  if(getACart(pset,"g",GiGzGe,EA,idx,-1,0))     // Try & use Cart. components
    {
    g   = GiGzGe.x();				//  Set the g factor
    gA  = GiGzGe.y();				//  Set the delzz value
    eta = GiGzGe.z();                           //  Set the eta value
    if(EA == EAzero)				//  Allow use of g with G
      getOrientation(pset,"G",EA,idx,-1,false); //  in parameter names
    return true;
    }

//    Try To Directly Read g-Factor , Anisotropy, Asymmetry, & Orientation
//                Interaction Via { Iqn, g, gA, Geta, GEAngles }

//  1.) If g or G has been specified, these parameters will be used
//  2.) We don't mind that eta is not set, default will be zero
//  3.) The spin quantum number Iqn was set earlier in this function
//  4.) We don't mind that no orientation is set, default is PAS
//  5.) Orientation set with either an Euler angle set or 3 individual angles
//  6.) Can mix parameters with g and G for this read

  string Pbase1("G");				// Base for parameter names
  string Pbase2("g");				// Another base name
  if(getGIso(pset, g, idx, false)) 		// If g-factor specified
    {						// look for gA, eta, orient
    getGA(pset, gA, idx, false);		//   Try for anisotropy
    if(!getAeta(pset,"g",eta,idx,-1,false))	//   Try for asymmetry
        getAeta(pset,"G",eta,idx,-1,false); 	//   using Geta or geta
    if(!getOrientation(pset,"g",EA,idx,-1,false))//  Try for orientation
        getOrientation(pset,"G",EA,idx,-1,false);//  using GEAngles or
    return true;				//   gEAngles or [G/g]alpha
    }						//   etc. etc.

/*
  if(getCSA(pset, csa, idx, false))             //   Try for anisotropy
    {                                           //   (when no PPM set)
    ppm = 0;                                    //   No isotropic shift
    getAeta(pset, Pbase, eta, idx,-1,false);    //   Try for asymmetry
    getOrientation(pset,Pbase,EA,idx,-1,false); //   Try for orientation
    return true;
    }
*/

 if(warn)                                      // Try as we might, cannot get
   IGerror(60, 1);                             // can't set from parameters
 return false;
 }

bool IntG::getGI(const ParameterSet& pset, const Isotope& ISI,
                    double& g, double& gA, double& eta, EAngles& EA,
                                         double& Bo, int idx, bool warn) const
  {
//              Check For Valid Isotope Type, Get Isotope Name

  string II = ISI.symbol();                     // String for isotope name
  bool TFI = true;                              // Flag if we know isotope
  if(!SpinCheck(ISI,TFI,true)) return false;	// Disallow nuclear spin

  if(warn)                                      // Try as we might, cannot get
    {                                           // interaction from parameters
    IGerror(50, 1);                             // Can't set from parameters
    IGerror(51, 1);                             // Parameter set insufficient
    }
  return false;
 }

// ----------------------------------------------------------------------------
//            Get Isotropic G Value (G-Factor) From A Parameter Set
// ----------------------------------------------------------------------------

/*         Input                SA      : CSA interaction (this)
                                pset    : A parameter set
           			g       : The g factor (unitless)
                                idx     : Index value
                                warn    : Warning output flag
           Output               TF      : True if g factor is set
                                          from parameters in pset
           Note                         : This WILL NOT alter the interaction
           Note                         : Parameters are G, G(#), g, g(#)
           Note                         : The value is unitless              */

bool IntG::getGIso(const ParameterSet& pset,double& g,int idx,bool warn) const
  {
  string Nidx = "";                             // Name addition per index
  if(idx >= 0)					// If index exists then set up
    Nidx += string("(")+Gdec(idx)+string(")");	// the parameter name to append
  string Gnames[2] = { "g", "G" }; 		// G parameter
  string pname, pstate;                         // Strings for parameter
  ParameterSet::const_iterator item;         // A pix into parameter list
  for(int i=0; i<2; i++) 			// Loop over possible g
    {                                           // parameters names
    pname = Gnames[i] +  Nidx;                  //      Parameter name
    item = pset.seek(pname);			//      Seek parameter in pset
    if(item != pset.end()) 			//      If it's been found
      {                                         //      parse the parameter
      (*item).parse(pname,g,pstate);	//      info and set giso
      return true;
      }
    }
  if(warn)                                      // If it hasn't been found
    {                                           // we can make an interaction
    IGerror(2, "g" + Nidx, 1);			// Can't find g(idx)
    IGerror(13, 1);                             // Can't set g-factor
    }
  g = 0;
  return false;
  }  

// ----------------------------------------------------------------------------
//                Get Anisotropic G Value From A Parameter Set
// ----------------------------------------------------------------------------

/*         Input                SA      : CSA interaction (this)
                                pset    : A parameter set
           			gA	: The g anisotropy (unitless)
                                idx     : Index value
                                warn    : Warning output flag
           Output               TF	: True if g anisotropy (unitless)
                                          set from parameters in pset
           Note                         : This WILL NOT alter the interaction
           Note                         : Parameters are GA,GA(#),gA,gA(#)   */

bool IntG::getGA(const ParameterSet& pset, double& gA,int idx, bool warn) const
  {
  string Nidx = "";                             // Name addition per index
  if(idx >= 0)					// If index exists then set up
    Nidx += string("(")+Gdec(idx)+string(")");	// the parameter name to append
  string Gnames[2] = { "ga", "GA" }; 		// G parameter
  string pname, pstate;                         // strings for parameter
  ParameterSet::const_iterator item;         // A pix into parameter list
  for(int i=0; i<2; i++) 			// Loop over possible g
    {                                           // parameters names
    pname = Gnames[i] +  Nidx;                  //   Parameter name
    item = pset.seek(pname);			//   Seek parameter in pset
    if(item != pset.end())			//   If it's been found
      {                                         //   parse the parameter
      (*item).parse(pname,gA,pstate);		//   info and set delzz
      return true;
      }
    }
  if(warn)                                      // If it hasn't been found
    {                                           // we can make an interaction
    IGerror(2, "gA" + Nidx, 1);			// Can't find g(idx)
    IGerror(14, 1);                             // Can't set g anisotropy
    }
  gA = 0;
  return false;
  }

// ----------------------------------------------------------------------------
//          Functions To Read The Applied Field Or Base Frequency
// ----------------------------------------------------------------------------

/*         Input                SA      : CSA interaction (this)
                                pset    : A parameter set
           			Bo	: The field strength (Gauss)
                                idx     : Index value
                                warn    : Warning output flag
           Output               TF	: True if Bo (Gauss)
                                          set from parameters in pset
           Note                         : This WILL NOT alter the interaction
           Note                         : Allowd parameters are listed below

			Field,  FieldT, gField(#), gFieldT(#)
                        GOmega, Omega,  GField(#), GFieldT(#)
                                        gOmega(#), GOmega(#)                 */

bool IntG::getField(const ParameterSet& pset, double& Bo,
                                                      int idx, bool warn) const
  {
  if(IntRank2::getField(pset,"g",Bo,idx,0) 	// Try for gField(#)/gFieldT(#)
  || IntRank2::getField(pset,"G",Bo,idx,0)	// Try for GField(#)/GFieldT(#)
  || IntRank2::getField(pset,"", Bo, -1,0)) 	// Try for  Field   / FieldT
    return true;

  if(getGOmega(pset,"g",Bo,idx,0)		// Try for gGOmega(#)
  || getGOmega(pset,"G",Bo,idx,0)		// Try for GGOmega(#)
  || getGOmega(pset,"", Bo, -1,0))		// Try for GOmega, if found
    {						// convert from freq. to field
    Bo = Hvalue(fabs(Bo), GFREE);		// Gauss <-- Hvalue(GHz, none)
    return true;
    } 

  if(getOmega(pset, "", Bo,-1, 0))		// Try for Omega, if found
    {						// convert from freq. to field
    Bo= fabs(Bo)*1.e6*HZ2RAD/GAMMA1H;		// Gauss <-- MHz proton based
    return true;
    }
  if(warn)                                      // If it hasn't been found
    {                                           // we can make an interaction
    IGerror(2, "gField", 1);			// Can't find gField(idx)
    IGerror(14, 1);                             // Can't set field strength
    }
  Bo = 0;					// Cannot set field, so
  return false;					// return that we failed
  }

// ----------------------------------------------------------------------------
//                        Complete Electron G Interaction
// ----------------------------------------------------------------------------

/* This function employs all of the "get*" functions in the above sections
   to parse a parameter set for all of the values needed to define an electron
   G interaction, namely { Iqn, g, gA, eta, alpha, beta, gamma }. In addition
   a static external magnetic field strength (akin to knowing a base frequency)
   is attempted to be determined. If the interaction definition is found, we
   set the interaction or return false. To support G interactions in spin
   systems, which typically know spin isotope types and their spin quantum
   numbers, there is a special function that takes a spin isotope for the
   spin quantum designation.                                                 */

bool IntG::setGI(const ParameterSet& pset, int idx, bool warn)
  {
// sosi
  double  Iqn;                                  // Our spin quantum number
  double  g;					// Our isotropic val (g factor)
  double  delz;					// Our anisotropic value
  double  eta;                                  // Our asymmetry        [0,1]
  EAngles EA;                                   // Our Euler angles (radians)
  double  Bo;                                   // Our field strengty (Gauss)
  if(getGI(pset,Iqn,g,delz,eta,EA,Bo,idx,warn))	// Try and get parameter vals
    {                                           // If successful, set interact
    GISO     = g;				//   Set isotropic value
    DELZZ    = delz;                            //   Set anisotropic value
    BoFIELD  = Bo;                              //   Set field strength  (G)
    double X = xi();                            //   Set interaction constant
    IntRank2::operator=(IntRank2(Iqn,X,eta,EA));//   Set rank 2 interaction
    return true;                                //   Return we are successful
    }
  if(warn) IGerror(24, 1);			// Insufficient parameters
  return false;
  }

bool IntG::setGI(const Isotope& II,const ParameterSet& pset,int idx,bool warn)
  {
  double  g;					// Our isotropic val (g factor)
  double  delz;					// Our anisotropic value
  double  eta;                                  // Our asymmetry        [0,1]
  EAngles EA;                                   // Our Euler angles (radians)
  double  Bo;                                   // Our field strengty (Gauss)
  if(getGI(pset,II,g,delz,eta,EA,Bo,idx,warn))	// Try and get parameter vals
    {                                           // If successful, set interact
    GISO     = g;				//   Set isotropic value
    DELZZ    = delz;                            //   Set anisotropic value
    BoFIELD  = Bo;                              //   Set field strength  (G)
    double X = xi();                            //   Set interaction constant
    IntRank2::operator=(IntRank2(II,X,eta,EA)); //   Set rank 2 interaction
    return true;                                //   Return we are successful
    }
  if(warn) IGerror(25, 1);
  return false;
  }

// ____________________________________________________________________________
// iii           ELECTRON G INTERACTION SPECIAL SPIN TENSORS
// ____________________________________________________________________________

/* These by-pass the generalized GAMMA spatial-spin tensors (in class spin_T)
   and just produces the rank 2 space-spin tensor used for the electron G 
   interaction directly.  They are scaled such that they are interaction
   independent and adhere to the relative scaling between components of true
   spherical tensors.  The rank 2 spatial tensors (used by this class) are
   scaled such that, when rotated they are normalized rank 2 spherical
   harmonics in a symmetric case (eta=0).                                    */

/*
void IntG::setTs(coord& B)

        // Input                GI	: G factor interaction (this)
	//			B       : Oriented, normalized field vector
        // Output               none    : G factor interaction spherical-
        //                                spin components are generated
	// Note				: No check is made to see if the
	//				  Tsph array has be made!

  {
  double Bnorm = B.Rad();			// Norm of field vector
  complex Bz(B.z()/Bnorm);			// Z axis component of B
  complex Bp = (B.x()/Bnorm, B.y()/Bnorm);	// B+ = Bx + i*By
  complex Bm = (B.x()/Bnorm, -B.y()/Bnorm); 	// B- = Bx - i*By
  int Ival = 2;					// Hilbert space of e-
  matrix IM = Im(Ival);				// The operator I-
  matrix IP = Ip(Ival);				// The operator I+
  matrix IZ = Iz(Ival);				// The operator Iz
  Tsph[0]= (2./sqrt(6.))*IZ;			// T20  = 2*Iz/sqrt(6)
  Tsph[1]= (-0.5*Bz)*IP; 			// T2m1 =-1/2(Bz * I+)
  Tsph[2]= (0.5*Bz)*IM;				// T2m1 = 1/2(Bz * I-)
  if(norm(Bp))
    {
    Tsph[0] += (-0.5*Bm/sqrt(6.))*IM		// T20 -= (B-I- + B+I+) 
            -  (0.5*Bp/sqrt(6.))*IP;		//      / 2.0*sqrt(6) 
    Tsph[1] -= (0.5*Bp)*IZ;			// T21  -= 1/2(B+ * Iz) 
    Tsph[2] += (0.5*Bm)*IZ;			// T2-1 += 1/2(B- * Iz)
    Tsph[3] =  (0.5*Bp)*IP;			// T22   = 1/2(B+ * I+)
    Tsph[4] =  (0.5*Bm)*IM;			// T2m2  = 1/2(B- * I-)
    }
  }
*/

/* For this interaction the values { BoFIELD, DELZZ, _XI } are inter-related.
   If any one of them is set then the others may need to be updated. These
   functions will do this updating.                                         

    Function                             Purpose
   -----------   --------------------------------------------------------------
      setXi      Sets the value of _XI (in IntRank2) for the current values of
                 DELZZ and BoFIELD
      setBo      Sets the value of BoFIELD then adjusts the value of _XI
                 based on this value & DELZZ. Bo assumed in Gauss           */

void IntG::setXi() { _XI = RT6PIO5*BoFIELD*DELZZ*1.3996e6; }
void IntG::setBo(double Bo) { BoFIELD = fabs(Bo);  setXi(); }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

/* Electron g interaction construction. For full construction users need to
   define the 1.) irreducible rank 2 spatial tensor, 2.) irreducible rank 2
   spin tensors, 3.) interaction strength. Since GAMMA uses normalized spatial
   and spin tensors, the first two are relatively simple. For the spatial 
   tensor we need either { eta } or { gxx, gyy, gzz } and for the spin tensor
   we need only the spin quantum value involved, normally 0.5 for an electron.
   The interaction strength is a bit more finicky since it is tied up in the
   isotropic value of the reducible rank 2 spatial G tensor (the g-factor), 
   the anisotropic value of the irreducible part (related to delzz), as well as
   to the external field strength that the electron is in.

   Since g = (1/3)[gxx+gyy+gzz] and delzz = gzz - giso, if the interaction is
   specified using Cartesian components the isotropic and delzz values are 
   automatically set.  If not then both values should be set independently.
   The isotropic value can always be taken as that of the free electron 
   (gfree = 2.00231928) and the delzz value set to some small positive value.
   The true anisotropy, delg = gzz - 1/2(gxx+gyy), is equal to 3/2 delzz.
   Lastly, to obtain the final interaction strength - used in constructing
   Hamiltonians, we need to have either a field strength specification
   or a resonance frequency specification. These may also be done after the
   construction has taken place.                                             */

// ____________________________________________________________________________
// A          G FACTOR INTERACTION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

IntG::IntG() : IntRank2()   
  { 
  GISO    = 0.0;			// No isotropic g-factor
  DELZZ   = 0.0;			// No anisotropy
  BoFIELD = 0.0;			// No external static field
  }

IntG::IntG(const IntG &G1) : IntRank2(G1)
  {
  GISO    = G1.GISO;			// Copy the isotropic g-factor
  DELZZ   = G1.DELZZ;			// Copy the anisotropy
  BoFIELD = G1.BoFIELD;			// Copy the external field strength
  }

// ----------------------------------------------------------------------------
//            Direct Constructors That Use Spherical Parameters
// ----------------------------------------------------------------------------

/* Here we try to use the set { Iqn, g, GA, Geta, GEAs , Bo }. The spin
   quantum number sets the dimension of the interaction spin Hilbert space.
   g is the isotropic g value or g-factor, gA is the g anisotropy, and geta
   the g asymmetry ([0,1]). The interaction orientation is taken from Euler
   angles gEA and the external static field strength is specified by Bo.
   We will allow for default settings of the asymmetry, orientation, and field
   strength so that those arguments need not be input during construction.
   Their default values are all 0.

        Input      GI      : G factor interaction (this)
                   g       : Isotropic G value (unitless)
                   gA      : G anisotropy      (unitless)
                   eta	   : G asymmetry value [0,1]
                   EA      : G orientation Euler angles (radians)
                   Bo      : External field strength    (Gauss)
        Output     none    : GI interaction constructed                      */

IntG::IntG(const string& IsoI,
                 double g, double gA, double eta, const EAngles& EA, double Bo) 
  {
  if(!SpinCheck(IsoI)) IGfatal(2);              // Insure we know this isotope
  Isotope II(IsoI);                             // Get isotope for this type
  if(!SpinCheck(II, true)) IGfatal(2);		// Insure we are are electron
  GISO      = g;				// Set the isotropic value
  DELZZ     = gA/1.5;				// Set the g anisotropy
  BoFIELD   = Bo;				// Set the applied field (G)
  double Iz = II.qn();                          // Get Iz value of spin
  double X  = xi();                             // Get interaction constant
  IntRank2::operator=(IntRank2(Iz,X,eta,EA));   // Use generic interaction`
  }

IntG::IntG(const Isotope& IsoI,
                 double g, double gA, double eta, const EAngles& EA, double Bo) 
  {
  if(!SpinCheck(IsoI, true)) IGfatal(2);	// Insure we are are electron
  GISO      = g;				// Set the isotropic value
  DELZZ     = gA/1.5;				// Set the g anisotropy
  BoFIELD   = Bo;				// Set the applied field (G)
  double Iz = IsoI.qn(); 			// Get Iz value of spin
  double X  = xi();                             // Get interaction constant
  IntRank2::operator=(IntRank2(Iz,X,eta,EA));   // Use generic interaction`
  }

IntG::IntG(double Iqn,
                 double g, double gA, double eta, const EAngles& EA, double Bo) 
  {
  if(!SpinCheck(Iqn, false)) IGfatal(2);	// Insure we are m*1/2
  GISO      = g;				// Set the isotropic value
  DELZZ     = gA/1.5;				// Set the g anisotropy
  BoFIELD   = Bo;				// Set the applied field (G)
  double X  = xi();                             // Get interaction constant
  IntRank2::operator=(IntRank2(Iqn,X,eta,EA));  // Use generic interaction`
  }

// ----------------------------------------------------------------------------
//            Direct Constructors That Use Cartesian Parameters
// ----------------------------------------------------------------------------

/* Here we try to use the set { Iqn, Gxx, Gyy, Gzz, GEAs, GOm }. The spin
   quantum number sets the dimension of the interaction spin Hilbert space.
   The values of Gxx, Gyy, and Gzz are unitless and used to determine the
   isotropic g factro GIso, G anisotropy GA, and G asymmetry ETA ([0,1]). The
   interaction orientation is taken from Euler angles GEAs and the external
   static field strength specified by the spin Larmor frequency GOm (Hz). We
   will allow for default settings of the orientation and Larmor frequency so
   that those arguments need not be input during construction.  Their
   default values are 0.

        Input      GI      : G factor interaction (this)
                   Gcart   : Cartesian tensor values (unitless)
                   EA      : G orientation Euler angles (radians)
                   Bo      : External field strength    (Gauss)
        Output     none    : GI interaction constructed 
	Note		   : We insist |gzz| >= |gyy| >= |gxx|               */

IntG::IntG(const string& IsoI, const coord& GxGyGz, const EAngles& EA, double Bo)
  {
  if(!SpinCheck(IsoI)) IGfatal(2);              // Insure we know this isotope
  Isotope II(IsoI);                             // Get isotope for this type
  if(!SpinCheck(II, true)) IGfatal(2);		// Insure we are an electron
  coord GiGzGe = AisoDelzEta(GxGyGz);           // Switch to spherical values
  GISO      = GiGzGe.x();                       // Set isotropic g-factor
  DELZZ     = GiGzGe.y();			// Set anisotropic delzz
  BoFIELD   = Bo;				// Set the applied field
  double Iz = II.qn();                          // Get Iz value of spin
  double E  = GiGzGe.z();                       // Get asymmetry value
  double X  = xi();                             // Get interaction constant
  IntRank2::operator=(IntRank2(Iz,X,E,EA));     // Use generic interaction
  }

IntG::IntG(const Isotope& II, const coord& GxGyGz, const EAngles& EA, double Bo)
  {
  if(!SpinCheck(II, true)) IGfatal(2);		// Insure we are an electron
  coord GiGzGe = AisoDelzEta(GxGyGz);           // Switch to spherical values
  GISO      = GiGzGe.x();                       // Set isotropic g-factor
  DELZZ     = GiGzGe.y();			// Set anisotropic delzz
  BoFIELD   = Bo;				// Set the applied field
  double Iz = II.qn();                          // Get Iz value of spin
  double E  = GiGzGe.z();                       // Get asymmetry value
  double X  = xi();                             // Get interaction constant
  IntRank2::operator=(IntRank2(Iz,X,E,EA));     // Use generic interaction
  }

IntG::IntG(double Iqn, const coord& GxGyGz, const EAngles& EA, double Bo)
  {
  if(!SpinCheck(Iqn, false)) IGfatal(2);	// Insure we are m*1/2
  coord GiGzGe = AisoDelzEta(GxGyGz);           // Switch to spherical values
  GISO      = GiGzGe.x();                       // Set isotropic g-factor
  DELZZ     = GiGzGe.y();			// Set anisotropic delzz
  BoFIELD   = Bo;				// Set the applied field
  double E  = GiGzGe.z();                       // Get asymmetry value
  double X  = xi();                             // Get interaction constant
  IntRank2::operator=(IntRank2(Iqn,X,E,EA));	// Use generic interaction
  }

IntG::IntG(const coord& GxGyGz, const EAngles& EA, double Bo)
  {
  coord GiGzGe = AisoDelzEta(GxGyGz);           // Switch to spherical values
  GISO      = GiGzGe.x();                       // Set isotropic g-factor
  DELZZ     = GiGzGe.y();			// Set anisotropic delzz
  BoFIELD   = Bo;				// Set the applied field
  double E  = GiGzGe.z();                       // Get asymmetry value
  double X  = xi();                             // Get interaction constant
  IntRank2::operator=(IntRank2(0.5,X,E,EA));	// Use generic interaction
  }

// ----------------------------------------------------------------------------
//                    Construction Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions construct the interaction from parameters in a specified
   GAMMA parameter set. This is most useful when working with multispin systems
   that are setting up many interactions over the system using a single 
   parameter file (or parameters read from an external ASCII file.

        Input                GI      : G factor interaction
                             pset    : Parameter set
                             II      : Spin isotope type
                             idx     : Interaction index (default -1->none)
                             warn    : Flag to warn if no interaction found
        Output               none    : GI interaction constructed
                                       for spin with quantum number qn
                                       and parameters in pset
                                or   : GI interaction constructed
                                       for spin with index idxI
                                       and parameters in pset 
	Note			     : If no isotope has been set then
				       we default to I=-1/2 (e-)
        Note                         : Constructions that explicitly use
        			       isotope labels are for support of
        			       use by spin systems                   */

IntG::IntG(const ParameterSet& pset, int idx, int warn)
     : IntRank2()
  {
  if(!setGI(pset,idx,warn?true:false) && warn)             // Try & set interaction
    {                                           // and warn if troubles
    if(warn > 1) IGfatal(2);
    else         IGerror(2,1);
    }
  }  

IntG::IntG(const Isotope& II, const ParameterSet& pset, int idx, int warn)
  {
  if(!setGI(II, pset,idx,warn?true:false) && warn)		// Try & set interaction
    {                                           // and warn if troubles
    if(warn > 1) IGfatal(2);
    else         IGerror(2,1);
    }
  }  

// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

void IntG::operator= (const IntG &GI1) 
  {
  IntRank2::operator=(GI1);		// Copy the rank2 interaction
  GISO    = GI1.GISO;			// Copy the isotropic value
  DELZZ   = GI1.DELZZ;			// Copy the anisotropic delzz
  BoFIELD = GI1.BoFIELD;		// Copy the applied field 
  }

IntG::~IntG () {}			// Nothing to delete here

// ____________________________________________________________________________
// B     G INTERACTION SPATIAL TENSOR COMPONENT ACCESS FUNCTIONS
//
//                       { Giso, del  , eta, theta, phi }
//                                  zz
// ____________________________________________________________________________
 
// ----------------------------------------------------------------------------
//                      Spatial Tensor Isotropic Value
// ----------------------------------------------------------------------------

/* By definition                  1 [               ]
                           g    = - | g  + g  + g   |
                            iso   3 [  xx   yy   zz ]

   This value is almost always a bit larger that the free electron G value
   which is 2.00231928.                                                      */
 
double IntG::iso()            const { return GISO; }
double IntG::g()              const { return GISO; }
void   IntG::iso(double giso)       { GISO = giso; } 
void   IntG::g(double   giso)       { GISO = giso; } 

// ----------------------------------------------------------------------------
//                      Spatial Tensor Anisotropy Value
// ----------------------------------------------------------------------------

/* The g anisotropy is

               ^      3                         1
              / \ g = - del   = g  - g   = g  - - [ g  + g  ]
              ---     2    zz    ||   |     zz  2    xx   yy
                                     ---

                             del   = g   - g
                                zz    zz    iso

    Since the G tensor is stored unitless, the returned anisotropy will
    be unitless as well. Note that the default base class IntRank2A also
    provides the GAMMA normalized anisotropy and delzz values. These are
    constants given by

                       1/2
                  [ 5 ]                  ^      3
          del   = |---]  = 0.51503      /_\ A = - del   = 0.77255
             zz   [6PI]                         2    zz                      */

// static double IntRank2A::delzz();                            INHERITED
// static double IntRank2A::delA();                             INHERITED
 
double IntG::aniso() const    { return 1.5*DELZZ; }
double IntG::gA()    const    { return 1.5*DELZZ; }
void   IntG::aniso(double ga) { DELZZ = ga/1.5; } 
void   IntG::gA(double    ga) { DELZZ = ga/1.5; }

//double IntG::delz()   const { return DELZZ; }
//double IntG::delzz()  const { return DELZZ; }
//void   IntG::delz(double  dz) { DELZZ = dz; } 
//void   IntG::delzz(double dz) { DELZZ = dz; }

double IntG::gdelz()         const { return DELZZ; }
void   IntG::gdelz(double g)       { DELZZ = g; }
 

// ----------------------------------------------------------------------------
//                      Spatial Tensor Asymmetry Value
// ----------------------------------------------------------------------------
 
// Note that eta spans [0, 1] and is defined to be (Axx-Ayy)/Azz where
// |Azz| >= |Ayy| >= |Axx| in the treatment contained herein.  These functions
// are inherited from the base class IntRank2A

//double IntG::eta( ) const		INHERITED	Get the asymmetry
//void   IntG::eta(double Qeta)		INHERITED	Set the asymmetry

//-----------------------------------------------------------------------------
//                Spherical Tensor Component Access (Normalized)
//-----------------------------------------------------------------------------

/* These allow one to access the irreducible spherical elements at a specific
   orientation. Most functions are derived from the base class IntRank2A. The
   exception are the PAS functions relative to those taking no angles. Since
   IntRank2A only exists in the PAS, here we reset the default function to
   use our only orientation and relegate the functions that take no arguments
   in IntRank2A to be end in PAS.  Note that these all use GAMMA scaling
   which sets A2m to be rank 2 spherical harmonics at all orientations if
   there is no asymmetry. 

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
//               Cartesian Tensor Component Access (Normalized)
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
//             Un-normalized Cartesian Components Of ESR G Tensor G
//                                                                 uv
// ----------------------------------------------------------------------------
 
/* These allow one to access the Cartesian elements of the full (unscaled)
   G tensor at a specific orientation without rotating the entire tensor.  In
   these functions theta is the angle down from the lab. frame z-axis and phi
   the angle over from the lab. frame x-axis.  If no angles are specified, then
   the current orientation is used.

                                    1/2
                            [ 6*PI ]
                      g   = | ---- |    * del   * A   + g
                       uv   [  5   ]         zz    uv    iso                 */

double IntG::gxx() const { return GISO + DELZZ*RT6PIO5*Axx(); }
double IntG::gyy() const { return GISO + DELZZ*RT6PIO5*Ayy(); }
double IntG::gzz() const { return GISO + DELZZ*RT6PIO5*Azz(); }
double IntG::gxy() const { return        DELZZ*RT6PIO5*Axy(); }
double IntG::gyx() const { return        DELZZ*RT6PIO5*Ayx(); }
double IntG::gxz() const { return        DELZZ*RT6PIO5*Axz(); }
double IntG::gzx() const { return        DELZZ*RT6PIO5*Azx(); }
double IntG::gyz() const { return        DELZZ*RT6PIO5*Ayz(); }
double IntG::gzy() const { return        DELZZ*RT6PIO5*Azy(); }

double IntG::gxx(double alpha, double beta, double gamma) const
                         { return GISO + DELZZ*RT6PIO5*Axx(alpha,beta,gamma); }
double IntG::gyy(double alpha, double beta, double gamma) const
                         { return GISO + DELZZ*RT6PIO5*Ayy(alpha,beta,gamma); }
double IntG::gzz(double alpha, double beta, double gamma) const
                         { return GISO + DELZZ*RT6PIO5*Azz(alpha,beta,gamma); }
double IntG::gyx(double alpha, double beta, double gamma) const
                         { return        DELZZ*RT6PIO5*Ayx(alpha,beta,gamma); }
double IntG::gxy(double alpha, double beta, double gamma) const
                         { return        DELZZ*RT6PIO5*Axy(alpha,beta,gamma); }
double IntG::gzx(double alpha, double beta, double gamma)
                        const { return   DELZZ*RT6PIO5*Azx(alpha,beta,gamma); }
double IntG::gzy(double alpha, double beta, double gamma)
                        const { return   DELZZ*RT6PIO5*Azy(alpha,beta,gamma); }
double IntG::gxz(double alpha, double beta, double gamma)
                        const { return   DELZZ*RT6PIO5*Axz(alpha,beta,gamma); }
double IntG::gyz(double alpha, double beta, double gamma)
                        const { return   DELZZ*RT6PIO5*Ayz(alpha,beta,gamma); }

double IntG::gxx(const EAngles& EA) const {return GISO+DELZZ*RT6PIO5*Axx(EA);}
double IntG::gyy(const EAngles& EA) const {return GISO+DELZZ*RT6PIO5*Ayy(EA);}
double IntG::gzz(const EAngles& EA) const {return GISO+DELZZ*RT6PIO5*Azz(EA);}
double IntG::gyx(const EAngles& EA) const {return      DELZZ*RT6PIO5*Ayx(EA);}
double IntG::gxy(const EAngles& EA) const {return      DELZZ*RT6PIO5*Axy(EA);}
double IntG::gzx(const EAngles& EA) const {return      DELZZ*RT6PIO5*Azx(EA);}
double IntG::gzy(const EAngles& EA) const {return      DELZZ*RT6PIO5*Azy(EA);}
double IntG::gxz(const EAngles& EA) const {return      DELZZ*RT6PIO5*Axz(EA);}
double IntG::gyz(const EAngles& EA) const {return      DELZZ*RT6PIO5*Ayz(EA);}

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
   effictively reorient the spatial tensor. 

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

	1.) GAMMA normalized Auv - Done With IntG::mx()
	2.) Typical guv values   - Done With IntG::mx(true);
	3.) Shown in lab frame   - Done With This Function

   For case 3.) the values are related to the GAMMA normalized (Auv) and
   typically presented values (guv) according to

                            1/2
                    [ 6*PI ]
              G   = | ---- |    * del   * A   + Kdel    * g
               uv   [  5   ]         zz    uv       u,v    iso

   where Kdel is a Kronecker delta function.                                 */

matrix IntG::Gmx() const
  {
  matrix SAmx = RT6PIO5*CartMx(1);			// Scaled A matrix
  matrix dmx(3,3,complex1,d_matrix_type);		// GIso matrix
  return SAmx + GISO*dmx;				// The G matrix
  }

double IntG::Field()         const { return BoFIELD; }
double IntG::GOmega()        const { return Omvalue(BoFIELD, GFREE); }
double IntG::Frequency()     const { return Omvalue(BoFIELD, GISO);  }
void   IntG::Field(double      Bo) { setBo(Bo); }
void   IntG::Frequency(double GOm) { setBo(Hvalue(GOm,  GISO)); }
 
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
  int    IntRank2T::HS()             const;                     INHERITED    */


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
  matrix IntRank2T::Tcomp(int comp) const;      // comp=[-2,2]  INHERITED

  matrix IntRank2T::T20(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T21(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m1(     const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T22(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m2(     const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m(int m,const vector<int>& HSs,int i,int j=-1) const;  */

// ____________________________________________________________________________
// D          SHIFT ANISOTROPY INTERACTION CONSTANT ACCESS FUNCTIONS
// ____________________________________________________________________________

/* The interaction constants is a scaling factor which is used to express all
   irreducible rank 2 interactions in a consistent fashion (class IntRank2)
   In fact, the interaction constant is maintained in the base class IntRank2
   as member variable _XI in radians/sec. This class simply relates the value
   of _XI to commonly used values that scale electron G interactions.

   Electron G interaction constants are field dependent and defined as

                         1/2                         1/2
               G   [6*pi]   beta * H           [6*pi]    GOm
             Xi  = |----| * -------- * del   = |----|  * ----- * del
                   [ 5  ]      h          zz   [ 5  ]    g          zz
                                                          iso
  
   where Xi and GOm are in Hz and H is in Gauss. GOm in this case is the
   frequency at (isotropic) resonance. For typical EPR fields the value of
   Xi is quite large. Below we use the fact that beta = 9.2741x10^-21 erg/G
   and h = 6.6262x10^-27 erg-s/cycle, so beta/h = 1.3996x10^6 Hz/G.  Note that
   since this factor affects the anisotropic part of the g interaction if
   the anisotropy is zero, then delzz is zero and xi is zero....            */

double IntG::xiOm(double GOm) const { return RT6PIO5*GOm*1.e9*DELZZ/GISO; }
double IntG::xiBo(double Bo)  const { return RT6PIO5*Bo*DELZZ*BOHRMAG/PLANCK; }

double IntG::xi() const                       // Overwrites inherited fct
  { return BoFIELD?RT6PIO5*BoFIELD*DELZZ*BOHRMAG/PLANCK:0.0; }

void IntG::xi(double X) const
  {
  IGerror(60,1);                                // Cannot set Xi
  IGerror(61,1);                                // Xi depends on Omega & DELZZ
  IGfatal(62);                                  // Use function xiOm or xiBo
  X = 0.0;                                      // to avoid complaints
  }                                             // from the compiler

// ____________________________________________________________________________
// E                ELECTRON G INTERACTION AUXILIARY FUNCTIONS
// ____________________________________________________________________________

/* Note that the isotropic component, herein GISO, sets the base frequency for
   resonance as referenced to a free electron.  Given either a spectrometer
   frequency or a static field strength this base resonance frequency can be
   calculated.
 
              iso   GISO*Beta*H    GISO*(9.2741x10^-21 erg/G)*H
             W    = -----------  = ---------------------------- = GOm
              G          h            6.6262x10^-27 erg/Hz
 
   If we use a typical GISO ~2.003 & H ~3000 Gauss then W ~10 GHz.
   Note that, as is the commonly used convention in ESR (unlike NMR), these
   formulae return frequencies in the LABORATORY frame.  Substituting GISO
   with GISO-GFREE would yeild the frequencies in the rotating frame of a
   free electron.                                                            */

// ----------------------------------------------------------------------------
//      Spectrometer Frequency, Resonance Field, And g-Factor Conversions
// ----------------------------------------------------------------------------

/* These functions allow for conversions between field, frequency, & g factor.
   Employing the previous equation, if we have any two values we know the 3rd.

                   g*Beta*H              h*Om              h*Om
              Om = --------         H = ------        g = ------
                      h                 g*Beta            H*Beta

           Input                Om      : Spectrometer frequency (GHz)
                                H       : The resonance field (Gauss)
                                g       : The electron g factor
           Output               g       : The electron g factor
                                H       : The resonance field (Gauss)
                                Om      : The spectrometer frequency (GHz) 
           Note				: GAMMA defined constants are used
					  g*Beta = Bohr Magneton  -> BOHRMAG
					  h      = Plancks Const. -> PLANK   */

double IntG::gvalue(double Om, double H)
  {
  double PCOVERBM = PLANCK*1.e4/BOHRMAG;	// h/beta in Gauss/Hz
  return Om*1.e9*PCOVERBM/H;			// Return g (unitless)
  }

double IntG::Hvalue(double Om, double g)
  {
  double PCOVERBM = PLANCK*1.e4/BOHRMAG;		// h/beta in Gauss/Hz
  return Om*1.e9*PCOVERBM/g;			// Return H in Gauss
  }

double IntG::Omvalue(double H, double g)
  {
  double PCOVERBM = PLANCK*1.e4/BOHRMAG;		// h/beta in Gauss/Hz
  return g*H/PCOVERBM;				// Return Om (GHz)
  }


// ----------------------------------------------------------------------------
//                    Functions To Get Effective G Factors
// ----------------------------------------------------------------------------

// Note that many other "g" factor functions have been previously defined.

//    hfl           1       [     2                   2                   ]
//   g    =  g    + - del   | 3cos (theta)-1 + eta*sin (theta)*cos(2*phi) |
//    eff     iso   2    zz [                                             ]


double IntG::geff_hfl() const { return GISO + DELZZ; }

double IntG::geff_hfl(double theta, double phi) const

  {
  double therad = theta*DEG2RAD;
  double Ctheta = cos(therad);
  double GA = 3*Ctheta*Ctheta-1; 
  if(eta())
    {
    double phirad = phi*DEG2RAD;
    double Stheta = sin(therad);
    double C2phi = cos(2*phirad);
    GA += eta()*Stheta*Stheta*C2phi;
    }
  return GISO + 0.5*DELZZ*GA;
  }
 
// ----------------------------------------------------------------------------
//                      Functions To Get Resonance Frequency
//	    As a Function of g(theta,phi) and Applied H Field Strength
// ----------------------------------------------------------------------------
 
//                                        -21
//              beta*g   *H(G)   9.2741x10   erg/G * g   * H(G)   g    * H(G)
//                    eff                             eff          eff
//   GOm(GHz) = -------------- = ------------------------------ = -----------
//                    h                        -18                 GHZ2GAUSS
//                                    6.6262x10   erg/GHz
 
        // Input                GI	: G factor interaction
	//			H       : Static field (Gauss)
        // Return               void	: The field at resonance (Gauss)
	// Note				: A symmetric tensor is assumed
	//				  in the para & perp functions
	// Note				: When the field Ho ~ 3 kG and
	//				  g ~ 2.003 then GOm ~ 10 GHz

double IntG::GOm_iso()  const { return  GISO*BoFIELD/GHZ2GAUSS; }
double IntG::GOm_para() const { return gzz()*BoFIELD/GHZ2GAUSS; }
double IntG::GOm_perp() const { return gxx()*BoFIELD/GHZ2GAUSS; }
double IntG::GOm_xx()   const { return gxx()*BoFIELD/GHZ2GAUSS; }
double IntG::GOm_yy()   const { return gyy()*BoFIELD/GHZ2GAUSS; }
double IntG::GOm_zz()   const { return gzz()*BoFIELD/GHZ2GAUSS; }

double IntG::GOm_iso(double Ho)  const { return  GISO*Ho/GHZ2GAUSS; }
double IntG::GOm_para(double Ho) const { return gzz()*Ho/GHZ2GAUSS; }
double IntG::GOm_perp(double Ho) const { return gxx()*Ho/GHZ2GAUSS; }
double IntG::GOm_xx(double Ho)   const { return gxx()*Ho/GHZ2GAUSS; }
double IntG::GOm_yy(double Ho)   const { return gyy()*Ho/GHZ2GAUSS; }
double IntG::GOm_zz(double Ho)   const { return gzz()*Ho/GHZ2GAUSS; }

// ----------------------------------------------------------------------------
//                    Functions To Get Resonance Field Strength
//		As a Function of g(theta,phi) and RF-Field Frequency
// ----------------------------------------------------------------------------
 
/*                                  -27
              h * GOm      6.6262x10   erg/Hz * GOm(Hz)   GHZ2GAUSS*GOm(GHz)
        H = -----------  = ---------------------------- = -----------------
            beta * g                -21                          g
                    eff    9.2741x10   erg/G * geff               eff
 
           Input                GI	: G factor interaction
	  			GOm     : Excitation frequency (GHz)
           Return               void	: The field at resonance (Gauss)
	   Note				: A symmetric tensor is assumed
	  				  in the para & perp functions
	   Note				: When the field GOm ~ 10 GHz and
	  				  g ~ 2.003 then H ~ 3 kGauss        */

double IntG::field_iso(double  GOm) const { return GHZ2GAUSS*GOm/GISO; }
double IntG::field_para(double GOm) const { return GHZ2GAUSS*GOm/gzz(); }
double IntG::field_perp(double GOm) const { return GHZ2GAUSS*GOm/gxx(); }
double IntG::field_xx(double   GOm) const { return GHZ2GAUSS*GOm/gxx(); }
double IntG::field_yy(double   GOm) const { return GHZ2GAUSS*GOm/gyy(); }
double IntG::field_zz(double   GOm) const { return GHZ2GAUSS*GOm/gzz(); }

// ____________________________________________________________________________
// F                        PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//      Functions To Make A Parameter Set From An Electron G Interaction
// ----------------------------------------------------------------------------

/* This class has no implicit knowledge of the electron spin. For a typical
   electron (say in nitroxide) I=1/2 whereas in a molecular ion (say Mo3++)
   higher spin states will occur from crystal field splittings. The preferred
   means of specifying an electron G interaction when there is no knowledge of
   the electron spin isotope types is that which is used for filling up the
   parameter set.  The base parameters of individual interactions in that case
   are { Gqn(#), gxx(#), gyy(#), gcc(#), galpha(#), gbeta(#), ggamma(#) }.

          Input                GI      : G factor interaction
                               idx     : Interaction index (default -1)
          Output               pset    : Parameter set with only
                            OR TF      : Return is FALSE if proper parameters
       					 for the interaction were not found
                                         electron G interaction parameters

    Function                                 Purpose
  ------------         -------------------------------------------------------
  ParameterSet         Convert interaction into a parameter set
  operator +=          Adds interaction to existing parameter set (friend)
  PSetAdd              Adds interaction to existing parameter set, this
                       allows for an index and a prefix in the parameters   */

IntG::operator ParameterSet( ) const
  { ParameterSet pset; pset += *this; return pset; }

void operator+= (ParameterSet& pset, const IntG &GI)
  { GI.PSetAdd(pset); }

void IntG::PSetAdd(ParameterSet& pset, int idx, int pfx) const
  {
  string suffx;					// Parameter suffix
  if(idx != -1)					// Only use suffix if idx
    suffx = string("(")+Gdec(idx)+string(")");	// is NOT -1
  string prefx;					// Parameter prefix
  if(pfx != -1) 				// Only use prefix if pfx
    prefx = string("[")+Gdec(pfx)+string("]");	// is NOT -1

  string pname  = prefx +string("GI") + suffx;	// Add GI spin quantum number
  string pstate = string("Spin Quantum Number");// Values like 0.5, 1.0 ...
  double pdatad = Izval();
  SinglePar par = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);

  pname  = prefx + string("gxx") + suffx;	// Cartesian x component
  pstate = string("Cartesian X Component)");
  pdatad = gxx();
  par    = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);

  pname  = prefx + string("gyy") + suffx;	// Cartesian x component
  pstate = string("Cartesian Y Component)");
  pdatad = gxx();
  par    = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);

  pname  = prefx + string("gzz") + suffx;	// Cartesian x component
  pstate = string("Cartesian Z Component)");
  pdatad = gxx();
  par    = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);
  } 

// ----------------------------------------------------------------------------
//  Functions To Output Electron G Interaction To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

/* These functions write the G interaction into an ASCII file in GAMMA
   Parameter Set format.  That is, the resulting ASCII file may be read into
   any GAMMA program to create a G interaction identical to that written.

        // Input                GI      : G factor interaction (this)
        //                      filename: Output file name
        //                  or  ofstr   : Output file stream
        //                      idx     : Interaction index (default -1)
        //                      pfx     : Interaction 2nd indx (def -1)
        //                      warn    : Warning level
        // Output               none    : Interaction is written as a
        //                                parameter set to file filename
        //                                or into output file stream ofstr   */

bool IntG::write(const string &filename, int idx, int warn) const
   {
   ofstream ofstr(filename.c_str());	// Open filename for input
   if(!ofstr.good())                    // If file bad then exit
     {
     IGerror(1, filename);		// Filename problems
     IGfatal(20);			// Fatal error
     }
   bool TF = write(ofstr,idx,warn?1:0);	// Use overload to write
   if(TF && warn)			// If write failed & warnings desired
    {					// take the appropriate action
    IGerror(51, 1);			// Problems writing to file
    if (warn>1) IGfatal(53);		// Fatal error
    return false;
    }
   ofstr.close();
   return TF;
   }

bool IntG::write(ofstream& ofstr, int idx, int warn) const
  {
  ParameterSet pset;                    // Declare a parameter set
  PSetAdd(pset, idx);			// Add in interaction parameters
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!pset.write(ofstr, w2))            // Use parameter set to write
    {                                   // out the interaction parameters
    if(warn)
      {
      IGerror(52, 1);                   // Problems writing to filestream
      if (warn>1) IGfatal(53);		// Fatal error
      }
    return false;
    }
  return true;
  }

// ____________________________________________________________________________
// H                  ELECTRON G INTERACTION INPUT FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//      Direct Read of G Intraction From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

/* These functions allow users to specify G interactions from parameters in an
   external ASCII file or from parameters in a GAMMA parameter set. These two
   methods are nearly equivalent since any ASCII file will be immediately read
   into a parameter set and the parameter set subsequently used to set the
   interaction.

       Input                GI      : G factor interaction
                            filename: Output file name
                            pset    : Parameter set
                            idx     : Interaction index (default -1->none)
                            warn    : Warning output label
                                       0 = no warnings
                                       1 = warnings
                                      >1 = fatal warnings
       Output               none    : G factor interaction is read in
                                      from parameters in file filename
                                      or those in the parameter set          */

bool IntG::read(const string &filename, int idx, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename, warn?1:0))    // Read in pset from file
    {                                   // If we cannot read the file at all
    if(warn)                            // then issue warnings as desired
      {
      IGerror(1, filename, 1);		//      Filename problems
      if(warn > 1) IGfatal(21);      //      Fatal error
      else         IGerror(21);         //      or a warning issued
      }
    return false;			//  Return that we failed!
    }
  if(!read(pset, idx, warn?1:0)) 	// Use overloaded function
    {                                   // If we cannot read the file at all
    if(warn)                            // then issue warnings as desired
      {
      IGerror(1, filename, 1);		//      Filename problems
      if(warn > 1) IGfatal(22);      //      Fatal error
      else         IGerror(22);         //      or a warning issued
      }
    return false;			//  Return that we failed!
    }
  return true;				//  Return that all OK
  }

bool IntG::read(const ParameterSet& pset, int idx, int warn)
  {
  bool TF = setGI(pset, idx, warn?1:0);
  if(!TF)                                       // If setGI didn't handle
    {                                           // setting interaction proper
    if(warn)                                    // then we'll issue some
      {                                         // warnings if desired
      if(warn > 1) IGfatal(23);			// Fatal error
       else        IGerror(23,1); 		// or a warning issued
      }
    return false; 				// Return that we failed
    }
  return TF;
  }

// ----------------------------------------------------------------------------
//  Interactive Reading of G Intractions From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

        // Input                GI      : G factor interaction
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
        //                      idx     : Interaction index
        // Output               string  : The parameter argn of array argc is
        //                                used to supply a filename from which
        //                                the interaction is read

string IntG::ask_read(int argc, char* argv[], int argn, int idx)
    {
    string filename;					// Name of system file  
    query_parameter(argc, argv, argn,   	        // Get from command
      "\n\tElectron G interaction filename? ",filename);// Or ask for it
    read(filename, idx);				// Read from filename
    return filename;
    }

string IntG::ask_read(int argc, char* argv[], int argn, const string& def, int idx)
  {
  string msg = "\n\tElectron G Interaction filename ["	// Query we will ask if
             + def + "]? ";                             // it is needed
  string filename = def;                                // Name of spin system file
  ask_set(argc,argv,argn,msg,filename);                 // Or ask for it
  read(filename, idx);                                  // Read system from filename
  return filename;                                      // Return filename
  }

// ____________________________________________________________________________
// I               ELECTRON G INTERACTION OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//   Functions That Generate Simple Strings To Simplify & Modularize Printing
//-----------------------------------------------------------------------------

string IntG::GfactorString() const
  {
  string Sg = string("G factor (isotropic):   ")
            + Gform("%10.7f", GISO);
  return Sg;
  }

string IntG::GAString() const
  {
  string Ssft("G Anisotropy:           ");
  Ssft += Gform("%10.7f", DELZZ);
  return Ssft;
  }

string IntG::GFieldString() const
  {
  string Ssft("Base Field (kG):        ");
  Ssft += Gform("%10.7f", BoFIELD/1000.);
  return Ssft;
  }

string IntG::GFrequencyString() const
  {
  string Ssft("Base Frequency (GHz):   ");
  Ssft += Gform("%10.7f", Omvalue(BoFIELD,GFREE)*1.e-9);
  return Ssft;
  }

//-----------------------------------------------------------------------------
//  Functions That Generate String Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

/* These functions return spatial tensor information in string format. This is
   done to facilitate printing, in particular printing of G spatial tensors
   from within rank 2 interactions.  In this case we make a list of information
   thats usually displayed to the left of a 3x3 Cartesian matrix rep. of G.

                       Spin Quantum Number:       I
                       Isotropic G Factor     xxxxx.xx 
                       G Anisotrpy    :       xxxxx.xx 
                       Asymmetry:                 x.xx
                       Euler Angle Alpha:       xxx.xx Degrees
                       Euler Angle Beta:        xxx.xx Degrees
                       Euler Angle Gamma:       xxx.xx Degrees               */


vector<string> IntG::InfoStrings() const
  {
  vector<string> SphStrings;
  SphStrings.push_back(StringI());                      // Spin quantum #
  SphStrings.push_back(GfactorString());		// G factor
  SphStrings.push_back(GAString());			// G anisotropy
  SphStrings.push_back(AsymmetryString());              // Asymmetry
  if(PAS())                                             // PAS orientation
    SphStrings.push_back(PASString());                  // if in PAS
  else                                                  // Non PAS orientation
    {                                                   // we use 3 angles
    SphStrings.push_back(AlphaString());                // Alpha
    SphStrings.push_back(BetaString());                 // Beta
    SphStrings.push_back(GammaString());                // Gamma
    }
  SphStrings.push_back(GFieldString());			// Base Field
  if(BoFIELD)
    {
    SphStrings.push_back(GFrequencyString());		// Base Frequency
    SphStrings.push_back(XiString());			// Interaction const
    }
  return SphStrings;
  }


//-----------------------------------------------------------------------------
//  Functions That Generate String Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

/*                                     TStrings
                                       ========

                         [ x.x, x.x, x.x]
                 T     = [ x.x, x.x, x.x]    m = { 0, 1, -1, 2, -2 }     
                  2,m    [ x.x, x.x, x.x]
                         [ x.x, x.x, x.x]

                                     CartAStrings
                                     ============

                                [ xx   xy   xz]   [ x.x, x.x, x.x]
                  CartAStrings: [g  , g  , g  ] = [ x.x, x.x, x.x]
                                [ yx   yy   yz]   [ x.x, x.x, x.x]
                                [g  , g  , g  ]
                                [ zx   zy   zz]                              */


// string* IntRank2T::TStrings(int M) const;                     INHERITED

vector<string> IntG::CartAStrings(const string& CSForm) const
  {
  int nstr = 6; 
  vector<string> Cartstrings(nstr);	// Vector of nstr strings
  Cartstrings[0] = "[g  , g  , g  ]";
  Cartstrings[1] = "[ xx   xy   xz]"
                 + string("   ")
                 + "["  + Gform(CSForm.c_str(), gxx())
                 + ", " + Gform(CSForm.c_str(), gxy())
                 + ", " + Gform(CSForm.c_str(), gxz()) + "]";
  Cartstrings[2] = "[g  , g  , g  ]"
                 + string(" = ")
                 + "["  + Gform(CSForm.c_str(), gyx())
                 + ", " + Gform(CSForm.c_str(), gyy())
                 + ", " + Gform(CSForm.c_str(), gyz()) + "]";
  Cartstrings[3] = "[ yx   yy   yz]"
                 + string("   ")
                 + "["  + Gform(CSForm.c_str(), gzx())
                 + ", " + Gform(CSForm.c_str(), gzy())
                 + ", " + Gform(CSForm.c_str(), gzz()) + "]";
  Cartstrings[4] = "[g  , g  , g  ]";
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

//-----------------------------------------------------------------------------
//       Functions That Generate Ouput Of The Rank 2 G Interaction
//-----------------------------------------------------------------------------

/* These functions will output information concerning the G interaction to any
   output stream.

           Input                GI	: Electron G interaction (this)
                                ostr	: Output stream
	  			fflag   : Format flag
	  				    false - Basic Parameters
	  				    true  - Full output
           Output               none    : G spatial tensor parameters 
                                          placed into the output stream      */

ostream& IntG::print(ostream& ostr, bool fflag) const
  {
  if(!GISO && !_XI)
    {
    string hdr = "Empty Electron G Interaction";
    string Spacer((40-hdr.length()/2), ' ');
    ostr << "\n\n" << Spacer << hdr << "\n";
    return ostr;
    }

/*	    Output Some Printed Lines Which Will Look Like The Following 

  			    Electron G Interaction
  
   Isotropic G value:         x.xxxxxx    [g  , g  , g  ] 
   G Anisotropy:              x.xxxxxx    [ xx   xy   xz]   [ x.x, x.x, x.x]
   G Asymmetry:               x.xxxxxx    [g  , g  , g  ] = [ x.x, x.x, x.x]
   Base Field (kG):           x.xxxxxx    [ yx   yy   yz]   [ x.x, x.x, x.x]
   Base Frequency (GHz):      x.xxxxxx    [g  , g  , g  ]
   Xi Value:                  x.xxxxxx    [ zx   zy   zz]                    */

  vector<string> Gstrs = InfoStrings();		// Array of info strings
  vector<string> Astrs = CartAStrings("%6.3f");	// Array of A (G) strings
  string hdr = "Electron G Interaction";	// Use this header
  string Spacer((40-hdr.length()/2), ' ');	// To center header
  ostr << "\n\n" << Spacer << hdr << "\n";	// Output centered header
  IntRank2A::print(ostr, Astrs, Gstrs);		// Print G interaction
  if(!fflag) return ostr;			// Exit if no more output

/*	    Now Output The Spherical Components, Spatial And Spin
       These Will Be Printed Lines Which Will Look Like The Following 
                        (Repeated For All 5 m Values)
  
  		A    = x.xxx		T     = [ x.x, x.x]
  		 2,m			 2,m	[ x.x, x.x]                  */

  printAT(ostr, G);
  ostr << "\n\n";
  return ostr;
  }

ostream& operator<< (ostream& out, const IntG& GI) { return GI.print(out); }

ostream& IntG::STList(ostream& ostr, int fflag)
  {
  string hdr("Electron G Spin Tensors List");
  int hl = hdr.length();
  string spacer = string(40-hl/2, ' ');
  ostr << "\n" << spacer << hdr << "\n";
  IntRank2::printList(ostr, SA, 1);
  return ostr;
  }

// ____________________________________________________________________________
// J               ELECTRON G INTERACTION HAMILTONIAN FUNCTIONS
// ____________________________________________________________________________

/* This section returns a variety of electron G interaction Hamiltonians. For
   irreducible rank 2 electron G interaction Hamiltonians, because this class
   uses standardized spatial & spin tensors they may all be derived from the
   same formula.

                        -2,2
                         ---           m
  H(alpha,beta,gamma) =  \   Xi  * (-1)  * A   (alpha,beta,gamma) * T    (i,j)
                         /                  2,m                      2,-m
                         ---
                          m

   The Hamiltonians will be returned in the single spin Hilbert space of the
   interaction, unless a composite Hilbert space and some spin index is
   supplied. All Hamiltonians will returned in the product basis as simple
   matrices. Their units will be radians/sec.                                */

// ----------------------------------------------------------------------------
//			     Isotropic G Hamiltonian
// ----------------------------------------------------------------------------

/* We don't use spherical tensors for this part, it doesn't rotate anyway!
   I've built a "rotating frame" allowance into Hiso!  This is very NMR-like
   but may never be used by ESR puritans.
 
            iso   beta * H              beta * H          G
           H    = -------- * g   * S  = -------- * A   * T    = GOm * S
            G         h       iso   z       h       0,0   0,0          z

   Or in the free e- rotating frame:

         iso   beta * H              beta * H          
        H    = -------- * g   * S  = -------- * (g   - g ) * S  = GOm * S
         G         h       eff   z       h        iso   e     z          z


	   Input		GI	: G factor interaction (this)
	  			Om	: Spectrometer field (GHz)
	  			rotflg  : Rotating frame flag
           Output               Hiso	: The isotropic G Hamiltonian
           Note                      	: Rotationally invariant about z
	   Note				: This will return in the spin Hilbert
	  				  space of dimension 2S+1
	   Note				: If rotflg set, Hiso is referenced 
	  				  to free e- rotating frame          */

matrix IntG::Hiso(bool rotfrm) const
  {
  double fact = Frequency();			// Electron frequency (Hz)
  if(rotfrm) fact *= (1.0 - GISO/GFREE);	// Switch to rotating frame
  return (fact * Iz(Ival));			// Return Ham. in Hz
  }

matrix IntG::Hiso(vector<int> HSs, int i, bool rotfrm) const
  {
  double fact = Frequency();			// Electron frequency (Hz)
  if(rotfrm) fact *= (1.0 - GISO/GFREE);	// Switch to rotating frame
  return fact * blow_up(Iz(Ival), HSs, i); 	// Return Ham. in Hz
  }

// ----------------------------------------------------------------------------
//               First Order G factor Interaction Hamiltonian
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
//  These Are For When The G factor Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------

/*  The secular part of the G Hamiltonian contains only those Hamiltonian
    components which explicitly commute with z axis rotations.  Here we have
 
       [1]                      G                            G      (0)
      H  (alpha,beta,gamma) = Xi * A   (alpha,beta,gamma) * T    = H
       G                            2,0                      2,0    G
 
   Identical to the Anisotropic G Hamiltonian in a high field limit!
   Note that this Hamiltonian must be ADDED to the isotropic G Hamiltonian.

           Input                GI      : G factor interaction
           Output               H0      : Secular part of the g anisotropy
                                          Hamiltonian (default basis, Hz)
           Note                         : Also called the 1st order GI
                                          interaction (perturbation theory)
           Note                         : Rotationally invariant about z
           Note                         : This will be referenced to the lab
                                          frame in accordance with standard
                                          ESR definitions of the G tensor.   */


//matrix IntG::H0()                             const;                INHERITED
//matrix IntG::H0(double A, double B, double G) const;                INHERITED
//matrix IntG::H0(const EAngles& EA)            const;                INHERITED

//matrix IntG::H0(const vector<int>& HSs, int i)                    const;  IHT
//matrix IntG::H0(const vector<int>& HSs, int i,
//                                    double A, double B, double G) const;  IHT
//matrix IntG::H0(const vector<int>& HSs, int i, const EAngles& EA) const;  IHT

// ----------------------------------------------------------------------------
//                   Full Electron G Interaction Hamiltonians
//                        These Have No Approximations
//     (They Are Correct in The Laboratory Frame, NOT in a Rotating Frame)
// ----------------------------------------------------------------------------

/* Here we make no approximations. As such, these Hamiltonians are not those
   used in the rotating frame. Here we simply sum over the 5 spatial and spin
   tensor components to produce the Hamiltonian.

                         -2,2
                          ---   G       m    G                          G
 H (alpha, beta, gamma) = \   Xi  * (-1)  * A   (alpha, beta, gamma) * T    (i)
  G                       /                  2,m                        2,-m
                           ---
                            m

           Input                GI      : G factor interaction
                                alpha   : Euler angle (radians)
                                beta    : Euler angle (radians)
                                gamma   : Euler angle (radians)
                                EA      : Euler angles (radians)
                                HSs     : Array of spin Hilbert spaces
                                i       : Spin index (in HSs)
           Output               H       : Matrix for G Hamiltonian
           Note                         : No HSs: return in the spin Hilbert
                                          space of dimension (2I+1)
                                          With HSs: return in composite spin
                                          Hilbert space                      */

// matrix IntG::H( ) const                                            INHERITED
// matrix IntG::H(double alpha, double beta, double gamma) const      INHERITED
// matrix IntG::H(const EAngles& EA) const                            INHERITED

// matrix IntG::H(const vector<int>& HSs, int i) const                INHERITED
// matrix IntG::H(const vector<int>& HSs, int i,                      INHERITED
//                     double alpha, double beta, double gamma) const
// matrix IntG::H(const vector<int>& HSs, int i,                      INHERITED
//                                           const EAngles& EA) const












matrix IntG::Hiso(double GOm, int rotflg) const	// Take GOm in GHz!
  {
  double fact = GOm*1.e9;			// Electron frequency(Hz)
  if(rotflg) fact *= (1.0 - GISO/GFREE);	// Switch to rotating frame
  return (fact * Iz(Ival));			// Return Ham. in Hz
  }

matrix IntG::Hiso(vector<int> HSs, int i, double GOm, int rotflg) const
  {
  double fact = GOm*1.e9;			// Electron frequency(Hz)
  if(rotflg) fact *= (1.0 - GISO/GFREE);	// Switch to rotating frame
  return fact * blow_up(Iz(Ival), HSs, i); 	// Return Ham. in Hz
  }

matrix IntG::Hiso_lab(double GOm) const
  { return Hiso(GOm); }

matrix IntG::Hiso_lab(vector<int> HSs, int i, double GOm) const
  { return Hiso(HSs, i, GOm); }
 

// ----------------------------------------------------------------------------
//               First Order G factor Interaction Hamiltonian
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
//  These Are For When The G factor Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------

/*  The secular part of the G Hamiltonian contains only those Hamiltonian
    components which explicitly commute with z axis rotations.  Here we have
 
    [1]             G                   G      (0)
   H  (the,phi) = Xi * A   (the,phi) * T    = H
    G                   2,0             2,0    G
 
                  1              [     2                   2               ]
                = - beta*H*del   | 3cos (the) - 1 + eta*sin (the)cos(2*phi)| Iz
                  2           zz [                                         ]
 
                        del
                           zz    [     2                   2               ]
                = GOm * ------ * | 3cos (the) - 1 + eta*sin (the)cos(2*phi)| Iz
                        2*g      [                                         ]
                           iso
 
   Identical to the Anisotropic G Hamiltonian in a high field limit!
   Note that this Hamiltonian must be ADDED to the isotropic G Hamiltonian.

           Input                GI      : G factor interaction
                                Om      : Spectrometer field (GHz)
           Output               H0      : Secular part of the g anisotropy
                                          Hamiltonian (default basis, Hz)
           Note                         : Also called the 1st order GI
                                          interaction (perturbation theory)
           Note                         : Rotationally invariant about z
           Note                         : H0 spin Hilbert space dim. is 2
           Note                         : This will be referenced to the lab
                                          frame in accordance with standard
                                          ESR definitions of the G tensor.   */

matrix IntG::H0X(double Om) const 
  { return (xiOm(Om)*Acomp(0))*Tcomp(0); }

matrix IntG::H0X(double Om, double theta, double phi) const
  { return (xiOm(Om)*A20(theta, phi, 0))*Tcomp(0); }

matrix IntG::H0X(vector<int>HSs, int i, double Om) const 
  { return (xiOm(Om)*Acomp(0))*blow_up(Tcomp(0), HSs, i); }

matrix IntG::H0X(vector<int>HSs, int i, double Om, double theta, double phi) const
  { return (xiOm(Om)*A20(theta, phi, 0))*blow_up(Tcomp(0), HSs, i); }

// These next two generate the Hamiltonian more explicitly.  The results
// should match the those of the above functions.
 
matrix IntG::H0Direct(double GOm) const
  { return (GOm*1.e9*DELZZ/(2.0*GISO))*Iz(Ival); }


matrix IntG::H0Direct(double GOm, double Theta, double Phi) const
 
  {
  double therad = Theta*DEG2RAD;
  double Ctheta = cos(therad);
  double GA = 3*Ctheta*Ctheta-1; 
  if(eta())
    {
    double phirad = Phi*DEG2RAD;
    double Stheta = sin(therad);
    double C2phi = cos(2*phirad);
    GA += eta()*Stheta*Stheta*C2phi;
    }
  return (GOm*1.e9*DELZZ/(2.0*GISO))*Iz(Ival);
  }

// ----------------------------------------------------------------------------
//                Full G Interaction Anisotropy Hamiltonians
//       These Have No Approximations & Are Independent of Field Strength
//     (They Are Correct in The Laboratory Frame, NOT in a Rotating Frame) 
// ----------------------------------------------------------------------------

/*                      ---
    GA              G   \       m    G             G 
   H  (the,phi) = Xi  * /   (-1)  * A (the,phi) * T
                        ---          2,m           2,-m   
                         m

                      [       2                   2               
                = K * | [ 3cos (the) - 1 + eta*sin (the)cos(2*phi) ] S 
                      [                                               z
                                                                              ]
                + sin(the)*[ cos(the)(3-eta*cos(2*phi))S + eta*sin(2*phi)*S ] |
                                                        x                  y  ]
   where
                                                       del
                            1                             zz
                        K = - beta * H * del   = GOm * ------
                            2               zz          2*g
                                                           iso
 
 Note: There are no m = +/-2 spin components in rank 2 G factor interactions
       when the applied field is along the z-axis in the laboratory frame
     (the lab frame is the frame to which angles theta and phi are referenced)
*/
 

matrix IntG::HADirect(double GOm) const

  { return (GOm*1.e9*DELZZ/(GISO))*Iz(Ival); }


matrix IntG::HADirect(double GOm, double Theta, double Phi) const

  {
  double therad = Theta*DEG2RAD;
  double phirad = Phi*DEG2RAD;
  double Ctheta = cos(therad);
  double Stheta = sin(therad);
  double C2phi = cos(2*phirad);
  double KZ = 3*Ctheta*Ctheta-1; 
  double S2phi = sin(2.0*phirad);
  if(eta()) KZ += eta()*Stheta*Stheta*C2phi;
  double KX = Stheta*Ctheta*(3.0-eta()*C2phi);
  double KY = Stheta*eta()*S2phi;
  matrix H = KZ*Iz(Ival) + KX*Ix(Ival) + KY*Iy(Ival); 
  return (GOm*1.e9*DELZZ/(2.0*GISO))*H;
  }

// These next two use the spherical tensor components explicitly.  The results
// should match the those of the above two functions.


 
	// Input		GI	: G factor interaction
	//			GOm	: Spectrometer field (GHz)
	// Output		H	: The G Hamiltonian
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2
	// Note				: This does NOT sum over the m=+/-2
	//				  components (i.e. field along +z)
	// Note                         : Xi depends on the field herein 

	// Input		GI	: G factor interaction
	//			GOm	: Spectrometer field (GHz)
	// Output		H	: The G Hamiltonian
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2
	// Note				: This does NOT sum over the m=+/-2
	//				  components (i.e. field along +z)
	// Note                         : Xi depends on the field herein 

matrix IntG::HA(double GOm) const
  {
  matrix Hmx = Acomp(0)*Tcomp(0);		// First set the m=0 terms
  if(norm(Acomp(1)))				// Then add in m=+/-1 terms
    Hmx -= (Acomp(1)*T2m1()+Acomp(-1)*T21());	// if they are present
  return xiOm(GOm)*Hmx;
  }

matrix IntG::HA(double GOm, double Theta, double Phi) const
  {
  matrix Hmx = A20(Theta,Phi,0)*T20();		// First set the m=0 terms
  complex z  = A21(Theta, Phi,0);			// Get A2,1 component
  if(norm(z))					// Then add in m=+/-1 terms
    Hmx -= (z*T2m1()+A2m1(Theta,Phi,0)*T21());	// if they are present
  return xiOm(GOm)*Hmx;
  }

matrix IntG::HA(vector<int> HSs, int i, double GOm) const
  {
  matrix Hmx = Acomp(0)*blow_up(Tcomp(0),HSs,i);// First set the m=0 terms
  if(norm(Acomp(1)))				// Then add in m=+/-1 terms
    Hmx -= (Acomp(1)*T2m1()+Acomp(-1)*T21());	// if they are present
  return xiOm(GOm)*Hmx;
  }

matrix IntG::HA(vector<int> HSs, int i, double GOm, double T, double P) const
  {
  matrix Hmx = A20(T,P, 0)*blow_up(T20(),HSs,i);// First set the m=0 terms
  complex z  = A21(T, P, 0);			// Get A2,1 component
  if(norm(z))					// Then add in m=+/-1 terms
    {
    Hmx -= z*blow_up(T2m1(),HSs,i);
    Hmx -= A2m1(T,P,0)*blow_up(T21(),HSs,i);	// if they are present
    }
  return xiOm(GOm)*Hmx;
  }


// ----------------------------------------------------------------------------
//                      Full G Interaction Hamiltonians
//       These Have No Approximations & Are Independent of Field Strength
//     (They Are Correct in The Laboratory Frame, NOT in a Rotating Frame)
// ----------------------------------------------------------------------------

//                      iso   ---     m     G                      G
//     H (theta,phi) = H    + \   (-1)  * Xi  * A   (theta,phi) * T
//      G               G     /                  2,m               2,-m
//                            ---
//
//                      iso    A
//                   = H    + H  (theta, phi)
//                      G      G 

        // Input                GI      : G factor interaction
        //                      Om      : Spectrometer frequency (GHz)
        //                      Theta   : Orientation down from +z (degrees)
        //                      Phi     : Orientation over from +x (degrees)
        // Output               H       : The G Hamiltonian (in Hz)
        // Note                         : This will return in the spin Hilbert
        //                                space of dimension 2
        // Note                         : This does NOT sum over the m=+/-2
        //                                components (i.e. field along +z)
        // Note                         : Xi depends on the frequency herein

 
matrix IntG::HDirect(double Om) const { return Hiso(Om) + HADirect(Om); }
matrix IntG::HDirect(double Om, double Theta, double Phi) const
  { return Hiso(Om) + HADirect(Om, Theta, Phi); }
 
// These next two use the spherical tensor components explicitly.  The results
// should match the those of the above two functions.
 
matrix IntG::H(double Om) const { return Hiso(Om) + HA(Om); }
matrix IntG::H(double Om, double Theta, double Phi) const
  { return Hiso(Om) + HA(Om, Theta, Phi); }

matrix IntG::H(vector<int> HSs, int i, double Om)
 const { return Hiso(HSs, i, Om) + HA(HSs, i, Om); }
matrix IntG::H(vector<int> HSs, int i, double Om, double Theta, double Phi) const
  { return Hiso(HSs, i, Om) + HA(HSs, i, Om, Theta, Phi); }

// ____________________________________________________________________________
// I              G FACTOR INTERACTION STANDARD INPUT FUNCTIONS
// ____________________________________________________________________________


//void IntG::ask(int argc, char* argv[], int& qn, 
//         double& Qnqcc, double& Qeta, double& Qtheta, double& Qphi, int Qflag)

        // Input                Q	: G factor interaction (this)
	//			argc    : Number of arguments
	//			argv    : Array of arguments
	//			qn      : Query value
	//			QI	: Spin quantum number
	//			Qnqcc   : G. coupling constant (Hz)
	//			Qeta    : G. asymmetry
	//			Qtheta  : G. orientation angle
	//			Qphi	: G. orientation angle
        //                      Qflag   : Flag is QCC or wQ requested
        // Output               none    : The values of qn, I, Qnqcc, Qeta,
	//				  Qtheta, and Qphi are set herein
	// Note				: This is INTERACTIVE!

//  {
//  query_parameter(argc, argv, qn++,			// Read in the I value
//      "\n\tSpin Quantum Value (1, 1.5, ..)? ", QI);
//  if(Qflag)
//    {
//    query_parameter(argc, argv, qn++,			// Read the frequency
//       "\n\tG Frequency(kHz)? ", Qnqcc);
//    Qnqcc *= 2.0*I*(2.0*I-1.0)/3.0;			// Set quad. coupling
//    }
//  else
//    {
//    query_parameter(argc, argv, qn++,			// Read in the coupling
//       "\n\tG Coupling (kHz)? ", Qnqcc);
//    }
//  Qnqcc *= 1.e3;					// Put this in Hz
//  query_parameter(argc, argv, qn++,			// Read in the coupling
//       "\n\tG Asymmetry [0, 1]? ", Qeta);
//  query_parameter(argc, argv, qn++,			// Read in theta angle
//  "\n\tOrientation Down From PAS z [1, 180]? ", Qtheta);
//  query_parameter(argc, argv, qn++,			// Read in phi angle
//   "\n\tOrientation Over From PAS x [0, 360]? ", Qphi);
//  }


//void IntG::askset(int argc, char* argv[])

        // Input                GI	: G factor interaction (this)
	//			argc    : Number of arguments
	//			argv    : Array of arguments
        // Output               none    : GI is set interactively
	// Note				: This is INTERACTIVE!

//  {
//  double GII=0, GIdelz=0, GIeta=0, GItheta=0, GIphi=0;
//  ask(argc,argv,qn,GII,GIdelz,GIeta,GItheta, GIphi, GIflag);
//  *(this) = IntG(GIdelz,GIeta,GItheta,GIphi);
//  }

 
//void IntG::askset(int GIflag)
 
        // Input                GI	: G factor interaction (this)
        //                      GIflag	: Flag is delzz or G requested
        // Output               none    : GI is set interactively
        // Note                         : This is INTERACTIVE!
 
//  {
//  int qn = 1000;
//  int argc = 1;
//  char* argv[1];
//  double GII, GIdelz, GIeta, GItheta, GIphi;
//  ask(argc,argv,qn,GII,GIdelz,GIeta,GItheta, GIphi,GIflag);
//  *(this) = IntG(GII,GIdelz,GIeta,GItheta,GIphi);
//  }

// sosik



// ____________________________________________________________________________
// J                 ELECTRON G INTERACTION OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//  Functions That Generate String Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

/*  Function   Arguments                     Output
    ========   =========   ===================================================
    TStrings       m       String array for the mth component of spin tensor T
    GAStrings              String array for various interaction values

            GAStrings()                              TStrings(m)

     Isotropic G value:      x.xxxxxxx                [ x.x, x.x, x.x]
     G Anisotropy:           x.xxxxxxx        T     = [ x.x, x.x, x.x]
     G Asymmetry:            x.xxxxxxx         2,m    [ x.x, x.x, x.x]
     Down From PAS z-Axis:   x.xx Deg.                [ x.x, x.x, x.x]
     Over From PAS x-Axis:   x.xx Deg.   
     Applied Field (kG):     x.xxxxxxx
                                              m = [0,4] => {0,1,-1,2,-2}     */


         
        // Input                GI	: G spatial tensor (this)
        //                      ostr	: Output stream
        // Output               none    : G spatial tensor parameters 
	//				  sent to output stream
	// Note				: Uses base class virtual overload

ostream& IntG::printSpherical(ostream& ostr)
  {
  ostr << "\t\n\t    Electron G Spatial Tensor";
  IntRank2A::printSpherical(ostr, 0);
  return ostr; 
  }

 
   
        // Input                GI	: G spatial tensor (this)
        //                      ostr    : Output stream
        // Output               none    : G spatial tensor parameters 
	//				  sent to output stream
	// Note				: Uses base class virtual overload

// Axx = (A22+A2m2)/2 - A20/sqrt(6)      Axy = -i*(A22-A21)/2 = Ayx
// Ayy = -(A22+A2m2)/2 - A20/sqrt(6)     Axz = -(A21-A2m1)/2  = Azx
// Azz = sqrt(2/3)*A20                   Ayz = i*(A21-A2m1)/2 = Azy
 
ostream& IntG::printCartesian(ostream& ostr, const string& CSF, int tpf)
  {
  ostr << "\t\n\t    Electron G Spatial Tensor";
  IntRank2A::printCartesian(ostr, CSF, 0);
  return ostr;
  }

 
/*
ostream& IntG::printCartesian(ostream& ostr, double phi, double theta)
         
        // Input                GI	: G spatial tensor (this)
        //                      ostr	: Output stream
        //                      phi     : Orientation angle (degrees)
        //                      theta   : Orientation angle (degrees)
        // Output               none    : G spatial tensor parameters 
	//				  sent to output stream
	// Note				: Uses base class virtual overload

// Axx = (A22+A2m2)/2 - A20/sqrt(6)      Axy = -i*(A22-A21)/2 = Ayx
// Ayy = -(A22+A2m2)/2 - A20/sqrt(6)     Axz = -(A21-A2m1)/2  = Azx
// Azz = sqrt(2/3)*A20                   Ayz = i*(A21-A2m1)/2 = Azy
 
  {
  ostr << "\t\n\t    Electron G Spatial Tensor";
  IntRank2A::printCartesian(ostr, theta, phi, 0);
  return ostr; 
  }
*/

 
         
        // Input                GI	: G interaction (this)
        //                      ostr	: Output stream
        // Output               none    : Full G Cartiesian tensor
	//				  sent to output stream

//                     [g  , g  , g  ] 
//                     [ xx   xy   xz]   [ x.x, x.x, x.x]
//                     [g  , g  , g  ] = [ x.x, x.x, x.x]
//                     [ yx   yy   yz]   [ x.x, x.x, x.x]
//                     [g  , g  , g  ]
//                     [ zx   zy   zz]
//

ostream& IntG::printCartG(ostream& ostr, int tflag) const
  {
  if(tflag)
    ostr << "\n\n" << string(30, ' ')
         << "Electron G Cartesian Spatial Tensor\n";
  vector<string> Astrs = CartAStrings("%6.3f");
  string newgline = string("\n") + string(20, ' ');
  for(int i=0; i<int(Astrs.size()); i++)
    ostr << newgline << Astrs[i];
  return ostr;
  }

#endif								// IntG.cc
