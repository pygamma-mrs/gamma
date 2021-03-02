/* IntCSA.cc ****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Chemical Shift Anisotropy Interaction 	Implementation		**
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
** A variable of type IntCSA represents a chemical shift anisotropy     **
** interaction defined for a particular nuclear spin (or spin in a spin **
** system).  The interaction is maintained in the form of spatial and   **
** spin tensors which are stored in irreducible spherical format.       **
**                                                                      **
** 1.) rank = 2         2.) symmetric           3.) eta = [0, 1]        **
**                                                                      **
** The tensor will internally be stored in irreducible spherical form   **
** at any chosen orientation.  The tensor principal axes will also be   **
** maintiained as well as a set of Euler angles that indicate the PAS   **
** orientation relative to some common axes.                            **
**                                                                      **
** Although not internal to this class, the interaction will blend with **
** a rank 2 spin tensor for the formation of an oriented quadrupolar    **
** Hamiltonian.                                                         **
**                                                                      **
** In addition to the functions contained in the base spatial tensor    **
** class, this IntCSA provides functions for building up the tensor	**
** and accessing the tensor elements from a "shielding" standpoint.	**
**                                                                      **
** The following defintions are used herein (Auv GAMMA normlized):	**
**                                                                      **
** 1.) PAS: |Azz| >= |Ayy| >= |Axx|					**
**                                                                      **
** 2.) PAS: Azz=2C, eta=(Axx-Ayy)/Azz, Axx=C(eta-1), Ayy=-C(1+eta)	**
**            								**
*************************************************************************/

#ifndef   IntCSA_cc_			// Is file already included?
#  define IntCSA_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <IntRank2/IntCSA.h>		// Include interface definition
#include <Basics/Gconstants.h>		// Include PI and other constants
#include <Basics/Gutils.h>		// Include parameter queries, errors
#include <Basics/Isotope.h>		// Include isotopes
#include <Basics/ParamSet.h>		// Include parameter sets
#include <HSLib/SpinOpSng.h>		// Include 1 spin operators
#include <Level1/coord.h>		// Include coordinates
#include <Matrix/row_vector.h>		// Include row_vectors
#include <Matrix/matrix.h>		// Include matrices
#include <HSLib/SpinSys.h>		// Include base spin systems
#include <Basics/StringCut.h>		// Include Gdec and Gform
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
// i             CLASS SHIFT ANISOTROPY INTERACTION ERROR HANDLING
// ____________________________________________________________________________
 
/*       Input                SA      : CSA interaction (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void IntCSA::ICerror(int eidx, int noret) const
  {
  string hdr("CSA Interaction");
  switch(eidx)
    {
    case 0: GAMMAerror(hdr,"Program Aborting.....",        noret); break;// (0)
    case 2: GAMMAerror(hdr,"Problems During Construction.",noret); break;// (2)
    case 8: GAMMAerror(hdr,"Theta (z Down) Beyond [0,180]",noret); break;// (8)
    case 9: GAMMAerror(hdr,"Phi (x Over) Outside [0, 360]",noret); break;// (9)
    case 10:GAMMAerror(hdr,"Asymmetry (eta) Beyond [0, 1]",noret); break;// (10)
    case 13:GAMMAerror(hdr,"Cannot Set Isotropic Shift",   noret); break;// (13)
    case 14:GAMMAerror(hdr,"Cannot Set Shift Anisotropy",  noret); break;// (14)
    case 15:GAMMAerror(hdr,"Setting Asymmetry to Zero",    noret); break;// (15)
    case 16:GAMMAerror(hdr,"Setting Default I=1/2 Value",  noret); break;// (16)
    case 18:GAMMAerror(hdr,"Electrons Have No CSA Int.",   noret); break;// (18)
    case 20:GAMMAerror(hdr,"Can't Write To Parameter File",noret); break;// (20)
    case 50:GAMMAerror(hdr,"Cant Determine From Params",   noret); break;// (50)
    case 51:GAMMAerror(hdr,"Parameter Set Insufficient",   noret); break;// (51)
    case 53:GAMMAerror(hdr,"Cannot Output Parameters",     noret); break;// (53)
    case 60:GAMMAerror(hdr,"Cant Set Interaction Constant",noret); break;// (60)
    case 61:GAMMAerror(hdr,"Xi Depends On Bo and CSA",     noret); break;// (61)
    case 62:GAMMAerror(hdr,"Use Function xiOm or CSA",     noret); break;// (62)
    default:GAMMAerror(hdr,"Unknown Error",                noret); break;
    }
  }

volatile void IntCSA::ICfatal(int eidx) const
  {
  ICerror(eidx,1);				// Output error message
  if(eidx) ICerror(0,1);			// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname 
      4     Use Of Deprecated Function pname
      5     Please Use Class Member Function pname
   default  Unknown Error - pname                                           */

void IntCSA::ICerror(int eidx, const string& pname, int noret) const 
  {
  string hdr("CSA Interaction");
  string msg;
  switch(eidx)
    {
    case 44:                                                   // (44)
      msg = string("Problems Setting Interaction, Index")
          + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 102:							// (102)
      msg = string("Odd Spin Quantum Value Of ") 
          + pname + string(" Specified?");
      GAMMAerror(hdr, msg, noret); break;
    case 103:							// (103)
      msg = string("Odd Asymmetry Value Of ") 
          + pname + string(" Specified?");
      GAMMAerror(hdr, msg, noret); break;
    case 110:							// (110)
      msg = string("Unknown Isotope Type Used - ")  + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 111:							// (111)
      msg = string("Not Defined For Electrons - ")  + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 130:							// (130)
      msg = string("Parameter ")
          + pname + string(" Is The Culprit!");
      GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;
    }
  }

volatile void IntCSA::ICfatal(int eidx, const string& pname, int noret) const 
  {
  ICerror(eidx, pname,1);			// Output error message
  if(eidx) ICerror(0,1);			// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }


// ____________________________________________________________________________
// ii                SHIFT ANISOTROPY INTERACTION SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions set up specific aspects of a shift anisotropy interaction.
   As the various interaction parameters interact, the functions MUST be
   private because their misuse could produce an inconsistent interaction.

   The goal here is quite simple. We must determine the following set of values
   for each shift anisotropy interaction: { Iqn,PPM,CSA,eta,alpha,beta,gamma }
   Complexity arises because we allow that a variety of parameters may be used
   to define these values. Additionally we allow that some parameters may be
   zero (e.g. eta) and that defaults may automatically be used. But, the end
   result remains the same, we seek the set of values that define a shift
   anisotropy interaction in GAMMA.                                          */

// ----------------------------------------------------------------------------
//                     Complete Shift Anisotropy Interaction
// ----------------------------------------------------------------------------

/* These functions will try and get all the parameters required to define a
   shift anisotropy interaction: { Iqn,PPM,CSA,eta,alpha,beta,gamma }. We 
   also see if we can figure out the spin Larmor frequency (which is akin
   to knowing the static magnetic field). To support CSA interactions in spin
   systems, which typically know spin isotope types and their spin quantum
   numbers, there is a separation between determining the spin quantum 
   number and the rest of the interaction.                                   */

bool IntCSA::getCI(const ParameterSet& pset,
         double& Iqn, double& ppm, double& csa, double& eta, EAngles& EA,
                                         double& Om, int idx, bool warn) const
  {
//                      Get The Spin Quantum Number
//                    (And Perhaps The Isotope Type)

  string pb("Iqn");				// Parameter base name
  string II; 					// String for isotope name
  Isotope ISI; 					// Isotope for spin
  bool TFI = false;                             // Flag if we know isotope
  if(getIso(pset,II,idx,0))			// 1. Try for isotope name
    {                                           //    If successful, check it
    if(!SpinCheck(II)) return false;		//    Insure valid isotope
    ISI = Isotope(II);                          //    An isotope of type II
    if(!SpinCheck(ISI,false,0)) return false;	//    Disallow electron spin
    TFI = true;                                 //    We know the isotope
    Iqn = ISI.qn();                             //    Know spin I quantum #
    }
  else
    {
    if(getIqn(pset,pb,Iqn,idx,false))		// 2. Try for spin quant. #
      { if(!SpinCheck(Iqn)) return false; }	//    Insure valid  qn
    else Iqn = 0.5;  				// 3. Use default qn of 1/2
    }

//                      Get The Spin Larmor Frequency
//     Not Mandatory! If Shift Specified In Hz We Will Need It Though

  if(TFI) getOm(pset,Om,II,idx,0);		// If we know isotope type
  else    getOmega(pset,"S",Om,idx,0);		// If we dont know isotope

// Try To Directly Read Shift Anisotropy Cartesian Spatial Tensor Components
//        Interaction Via { Iqn, Sxx, Sxy, Sxz, Syy, Syz, Szz }

//  1.) Suv values set the shift, anisotropy & asymmetry: {SISO, SCSA, ETA} 
//  2.) If any Suv off-diagonals (u != v) set, S array sets orientation also
//  3.) If no Suv off-diagonals, orientation set with specified Euler angles
//  4.) The input Suv values are taken to be in PPM
//  5.) We don't mind that no orientation is set, default is PAS
//  6.) If used, Euler angles by either a set or 3 individual angles

  coord SiSzSe;                                 // For Siso, Szz, Seta
  if(getACart(pset,"S",SiSzSe,EA,idx,-1,0))	// Try & use Cart. components
    {
    ppm = SiSzSe.x();				//  Set the shift value 
    csa = 1.5*SiSzSe.y();			//  Set the anisotropy value
    eta = SiSzSe.z();                           //  Set the eta value
    return true;
    }

//      Try To Directly Read Shift, Anisotropy, Asymmetry, & Orientation
//             Interaction Via { Iqn, PPM, CSA, Seta, SEAngles }

//  1.) If PPM has been specified, these parameters will be used
//  2.) If we know the Larmor frequency (field strength) v can also be used
//  3.) We don't mind that eta is not set, default will be zero
//  4.) Iqn & Sqn were set in the 1st section of this function
//  5.) We don't mind that no orientation is set, default is PAS
//  6.) Orientation set with either an Euler angle set or 3 individual angles

  string Pbase("S");				// Base for parameter names
  if(getPPM(pset, ppm, idx, Om, false))		// If PPM has been specified
    {                                           // look for CSA, eta, orient
    getCSA(pset, csa, idx, false);		//   Try for anisotropy
    getAeta(pset, Pbase, eta, idx, -1, false);	//   Try for asymmetry
    getOrientation(pset,Pbase,EA,idx,-1,false);	//   Try for orientation
    return true;
    }

  if(getCSA(pset, csa, idx, false))		//   Try for anisotropy
    {						//   (when no PPM set)
    ppm = 0;					//   No isotropic shift
    getAeta(pset, Pbase, eta, idx,-1,false);	//   Try for asymmetry
    getOrientation(pset,Pbase,EA,idx,-1,false);	//   Try for orientation
    return true;
    }

  if(warn)					// Try as we might, cannot get
    {						// interaction from parameters 
    ICerror(50, 1);				// Can't set from parameters
    ICerror(51, 1);				// Parameter set insufficient
    }
  return false;
  }

bool IntCSA::getCI(const ParameterSet& pset, const Isotope& ISI,
                   double& ppm, double& csa, double& eta, EAngles& EA,
                                         double& Om, int idx, bool warn) const
  {
//              Check For Valid Isotope Type, Get Isotope Name

  string II = ISI.symbol();			// String for isotope name
  bool TFI = true;				// Flag if we know isotope
  if(!SpinCheck(ISI,TFI,0)) return false;	// Disallow electron spin

//                      Get The Spin Larmor Frequency
//     Not Mandatory! If Shift Specified In Hz We Will Need It Though

  getOm(pset,Om,II,idx,0);			// If we know isotope type
   
// Try To Directly Read Shift Anisotropy Cartesian Spatial Tensor Components
//        Interaction Via { Iso, Sxx, Sxy, Sxz, Syy, Syz, Szz }

//  1.) Suv values set the shift, anisotropy & asymmetry: {SISO, SCSA, ETA} 
//  2.) If any Suv off-diagonals (u != v) set, S array sets orientation also
//  3.) If no Suv off-diagonals, orientation set with specified Euler angles
//  4.) The input Suv values are taken to be in PPM
//  5.) We don't mind that no orientation is set, default is PAS
//  6.) If used, Euler angles by either a set or 3 individual angles

  coord SiSzSe;                                 // For Siso, Szz, Seta
  if(getACart(pset,"S",SiSzSe,EA,idx,-1,0))	// Try & use Cart. components
    {
    ppm = SiSzSe.x();				//  Set the shift value 
    csa = 1.5*SiSzSe.y();			//  Set the anisotropy value
    eta = SiSzSe.z();                           //  Set the eta value
    return true;
    }

//      Try To Directly Read Shift, Anisotropy, Asymmetry, & Orientation
//             Interaction Via { Iso, PPM, CSA, Seta, SEAngles }

//  1.) If PPM has been specified, these parameters will be used
//  2.) If we know the Larmor frequency (field strength) v can also be used
//  3.) We don't mind that eta is not set, default will be zero
//  4.) Iqn & Sqn were set in the 1st section of this function
//  5.) We don't mind that no orientation is set, default is PAS
//  6.) Orientation set with either an Euler angle set or 3 individual angles

  if(getPPM(pset, ppm, idx, Om, 0))		// If PPM has been specified
    {                                           // look for CSA, eta, orient
    string Pbase("S");                          //   Base for parmeter names
    getCSA(pset, csa, idx, 0);			//   Try for anisotropy
    getAeta(pset, Pbase, eta, idx, 0);		//   Try for asymmetry
    getOrientation(pset,Pbase,EA,idx,0);	//   Try for orientation
    return true;
    }

  if(warn)					// Try as we might, cannot get
    {						// interaction from parameters 
    ICerror(50, 1);				// Can't set from parameters
    ICerror(51, 1);				// Parameter set insufficient
    }
  return false;
  }

// ----------------------------------------------------------------------------
//            Get Isotropic Chemical Shift Value From A Parameter Set
// ----------------------------------------------------------------------------

/*	   Input		SA	: CSA interaction (this)
                                pset    : A parameter set
                                idx     : Index value
                                warn    : Warning output flag
           Output               PPM	: Chemical shift anisotropy (PPM)
                                          constant from parameters in pset
           Note                         : This WILL NOT alter the interaction
	   Note				: Parameters are PPM, PPM(#)
	   Note				: Only Accepted In PPM Here         */
              
bool IntCSA::getPPM(const ParameterSet& pset,
                              double& ppm, int idx, double Om, bool warn) const
  {
  ParameterSet::const_iterator item;         // A pix into parameter list
  string Nidx = "";                             // Name addition per index
  if(idx >= 0) Nidx                             // If index exists then set up
    += string("(") + Gdec(idx) + string(")"); 	// the parameter name to append

//              Look For Isotropic Shift Using Parameter PPM

  string pname =  string("PPM") + Nidx; 	// PPM parameter
  item = pset.seek(pname);			// Seek parameter in pset
  if(item != pset.end())			// If parameter was found
    {						// parse the parameter
    string pstate;				// Temp String for parameter
    (*item).parse(pname,ppm,pstate);		// Glean out parameter value
    return true;				// Return the shift value
    }

//              Look For Isotropic Shift Using Parameter v
//            (Only Allowed If We Know Our Larmor Frequency) 

  if(Om)					// Omega must not be zero
    {
    string pnamev =  string("v") + Nidx;	 	// v parameter
    item = pset.seek(pnamev);			// Seek parameter in pset
    if(item != pset.end())			// If parameter was found
      {						// parse the parameter
      string pstate;				// Temp String for parameter
      (*item).parse(pname,ppm,pstate);		// Glean out parameter value
      ppm /= Om;				// Switch from Hz to ppm
      return true;				// Return the shift value
      }
    }

  if(warn)					// If it hasn't been found
    {						// we can make an interaction 
    ICerror(2, pname, 1);			// Can't find PPM(idx)
    ICerror(13, 1);				// Can't set shift
    }
  return false;					// Could not read PPM
  }  
              
// ----------------------------------------------------------------------------
//             Get Shielding Anisotropy Value From A Parameter Set
// ----------------------------------------------------------------------------

/*	   Input		SA	: CSA interaction (this)
                                pset    : A parameter set
                                idx     : Index value
                                warn    : Warning output flag
           Output               CSA	: Chemical shift anisotropy (PPM)
                                          constant from parameters in pset
           Note                         : This WILL NOT alter the interaction
	   Note				: Parameters are CSA, CSA(#)

          Try & Read SA Interaction Chemical Shift Anisotropy { CSA }                                 
                                                                            
     ^          3                                     1                      
    / \ sigma = - del   = sigma  - sigma   = sigma  - - [ sigma  + sigma  ] 
    ---         2    zz        ||        |        zz  2        xx       yy  
                                        ---                                  */
              
bool IntCSA::getCSA(const ParameterSet& pset,
                                         double& csa, int idx, bool warn) const
  {
  ParameterSet::const_iterator item;         // A pix into parameter list
  string Nidx = "";                             // Name addition per index
  if(idx >= 0) Nidx                             // If index exists then set up
    += string("(") + Gdec(idx) + string(")"); 	// the parameter name to append
  string pname =  string("CSA") + Nidx; 	// CSA parameter
  item = pset.seek(pname);			// Seek parameter in pset
  if(item != pset.end())			// If parameter was found
    {						// parse the parameter
    string pstate;				// Temp String for parameter
    (*item).parse(pname,csa,pstate);		// Glean out parameter value
    return true;				// Return the CSA value
    }
  if(warn)					// If it hasn't been found
    {						// we can make an interaction 
    ICerror(2, pname, 1);			// Can't find CSA(idx)
    ICerror(14, 1);				// Can't set CSA
    }
  return false;					// Could not read CSA
  }  

// ----------------------------------------------------------------------------
//                       Get Spin Larmor Frequency (MHz)
// ----------------------------------------------------------------------------

/* Shift anisotropy interactions are field dependent. That is, their strength
   directly results from an externally applied magnetic field. The stronger
   the field the stronger the interaction. However, our interaction does NOT
   store the magnetic field strength. Rather, it stores a spin Larmor frequency
   because the Larmor frequency combines the field strength with the spin
   gyromagnetic ratio. Such a combination is needed because both the field and
   gyromagnetic ratio are needed to properly scale the interaction. Because we
   only know the interaction spin quantum number we have no intrinsic knowledge
   of any gyromagnetic ratio. We cannot convert between field and frequency
   without it. So, we track the Larmor frequency and spin quantum number rather
   than the field and the gyromagnetic ratio (and the spin quantum number).

   We can use the field strength and gyromagnetic ratio to set our Larmor
   frequency. That is done from the parameters "Iso" and "Field". Alternatively
   we can just set the Larmor frequency directly with the parameter "Omega".
   If we know the isotope type we will allow useage of "Field" to set the
   Larmor frequency. Else, we will look for "Omega" directly.

  	   Input		SA	: CSA interaction (this)
                                pset    : A parameter set
                                idx     : Index value
                                warn    : Warning output flag
           Output               Om 	: Spin Larmor frequency (MHz)
                                          from parameters in pset
           Note                         : Does NOT alter the interaction
	   Note				: Parameters are CSA, CSA(#)         */

bool IntCSA::getOm(const ParameterSet& pset,
                       double& Om, const string& II, int idx, bool warn) const
  {
  double Bo;					// For field strength
  string pnb("S");				// Parameter base name

  if(getField(pset,pnb,Bo,idx,false))		// If field (SField, SFieldT)
    {						// we convert to our Larmor
    Om = Bo*GAMMA1H*1.e-6/HZ2RAD;		// frequency (Gauss -> MHz)
    return true;
    }
  else if(getField(pset,"",Bo,-1,false))	// If field (Field, FieldT)
    {						// we convert to our Larmor
    Om = Bo*GAMMA1H*1.e-6/HZ2RAD;		// frequency (Gauss -> MHz)
    return true;
    }
  else if(getOmega(pset,pnb,Om,idx,false))	// If freq. (SOmega) specified
    { 						// it is our Larmor frequency
    return true;
    }
  else if(getOmega(pset,"",Om,-1,false))	// If freq. (Omega) specified
    {						// it is for 1H, our Larmor is
    Om *= Isotope(II).gamma()/GAMMA1H;		// scaled from proton Larmor
    return true;
    }
  else if(getGOmega(pset,"",Om,-1,false))	// If freq. (GOmega) specified
    {						// it is for e-, our Larmor is
    Isotope E("e-"); 				// scaled from electron freq.
    Om = Isotope(II).gamma()/fabs(E.gamma());
    Om *= 1.e3;
    return true;
    }
  Om = 0.0;
  return false;					// We cannot find any Larmor
  }						// frequency for ourself

// ----------------------------------------------------------------------------
//                     Complete Shift Anisotropy Interaction
// ----------------------------------------------------------------------------

/* This function employs all of the "get*" functions in the above sections
   to parse a parameter set for all of the values needed to define a shift
   anisotropy interaction, namely { Iqn,PPM,CSA,eta,alpha,beta,gamma }. In
   addition the spin Larmor frequency (akin to knowing the static magnetic
   field and gyromagnetic ration) is attempted to be determined. If the 
   interaction definition is found, we set the interaction or return false.
   To support CSA interactions in spin systems, which typically know spin
   isotope types and their spin quantum numbers, there is a special function
   that takes a spin isotope for the spin quantum designation.               */

bool IntCSA::setCI(const ParameterSet& pset, int idx, int warn)
  {
  double  Iqn;					// Our spin quantum number
  double  ppm;					// Our isotropic shift  (PPM)
  double  csa;					// Our shift anisotropy (PPM)
  double  eta;					// Our asymmetry        [0,1]
  EAngles EA;					// Our Euler angles (radians)
  double  Om;					// Our Larmor frequency (MHz)
  if(getCI(pset,Iqn,ppm,csa,eta,EA,Om,idx,warn?true:false))// Try and get parameter vals
    {						// If successful, set interact
    SISO     = ppm;				//   Set isotropic val   (PPM)
    SCSA     = csa;				//   Set anisotropic val (PPM)
    SOMEGA   = Om;				//   Set Larmor freq.    (MHz)
    double X = xi();				//   Set interaction constant
    IntRank2::operator=(IntRank2(Iqn,X,eta,EA));//   Set rank 2 interaction
    return true;				//   Return we are successful
    }
  return false;
  }

bool IntCSA::setCI(const Isotope& II,const ParameterSet& pset,int idx,int warn)
  {
  double  ppm;					// Our isotropic shift  (PPM)
  double  csa;					// Our shift anisotropy (PPM)
  double  eta;					// Our asymmetry        [0,1]
  EAngles EA;					// Our Euler angles (radians)
  double  Om;					// Our Larmor frequency (MHz)
  if(getCI(pset,II,ppm,csa,eta,EA,Om,idx,warn?true:false))	// Try and get parameter vals
    {						// If successful, set interact
    SISO     = ppm;				//   Set isotropic val   (PPM)
    SCSA     = csa;				//   Set anisotropic val (PPM)
    SOMEGA   = Om;				//   Set Larmor freq.    (MHz)
    double X = xi();				//   Set interaction constant
    IntRank2::operator=(IntRank2(II,X,eta,EA));	//   Set rank 2 interaction
    return true;				//   Return we are successful
    }
  return false;
  }

// ----------------------------------------------------------------------------
//      Set Shift Anisotropy Interaction Spatial and Spin Tensor Components
// ----------------------------------------------------------------------------

/* These bypass the generalized GAMMA spatial (class space_T) and spin
   (class spin_T) tensors, instead relying on the irreducible rank 2
   base classes IntRank2A and IntRank2T which are more suitable for the
   treatment of shift anisotropy interactions.  Both are scaled such that they
   are interaction independent and adhere to the relative scaling between
   components of true spherical tensors.  The rank 2 spatial tensors are
   scaled such that, when rotated they are normalized rank 2 spherical
   harmonics in a symmetric case (eta=0).                                    */

// void IntCSA::setTs()		INHERITED       Set spin tensor components
// void IntCSA::setAs()		INHERITED	Set spatial tensor components

/*
// sosi - this is not active.  Is it even needed?

void IntCSA::setTs(coord& B)

        // Input                SA	: CSA interaction (this)
	//			B       : Oriented, normalized field vector
        // Output               none    : CSA interaction spherical-
        //                                spin components are generated
	// Note				: No check is made to see if the
	//				  Tsph array has be made!

  {
  double Bnorm = B.Rad();			// Norm of field vector
  complex Bz(B.z()/Bnorm);			// Z axis component of B
  complex Bp = (B.x()/Bnorm, B.y()/Bnorm);	// B+ = Bx + i*By
  complex Bm = (B.x()/Bnorm, -B.y()/Bnorm); 	// B- = Bx - i*By
  int Ival = int(2.*I + 1);			// For 1 spin SOp functions
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

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A          SHIFT ANISOTROPY INTERACTION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                   Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

IntCSA::IntCSA() : IntRank2()
  { 
  SISO   = 0.0;				// No shift isotropy (PPM)
  SCSA   = 0.0;				// No shift anisotropy (PPM)
  SOMEGA = 0.0; 			// No Larmor frequency (Hz)
  }

IntCSA::IntCSA(const IntCSA &S) : IntRank2(S) 
  {
  SISO   = S.SISO;			// Copy the isotropy  (PPM)
  SCSA   = S.SCSA;			// Copy the anisotropy (PPM)
  SOMEGA = S.SOMEGA;			// Copy the Larmor frequency (Hz)
  }

// ----------------------------------------------------------------------------
//            Direct Constructors That Use Spherical Parameters
// ----------------------------------------------------------------------------

/* Here we try to use the set { Iqn, SIso, SCSA, Seta, SEAngles , SOm }. The
   spin quantum number sets the dimension of the interaction spin Hilbert space
   SIso is the isotropic chemical shift (PPM), SCSA is the chemical shift
   anisotropy (PPM), and Seta the shift asymmetry ([0,1]). The interaction
   orientation is taken from Euler angles SEAngles & the external static field
   strength is specified by the spin Larmor frequency SOm (Hz). We will
   allow for default settings of the asymmetry, orientation, and Larmor
   frequency so that those arguments need not be input during construction.
   Their default values are all 0.                                           */

IntCSA::IntCSA(const string& IsoI,
           double Siso, double Scsa, double eta, const EAngles& EA, double Om)
  {
  if(!SpinCheck(IsoI)) ICfatal(2);		// Insure we know this isotope
  Isotope II(IsoI);				// Get isotope for this type
  if(!SpinCheck(II, false)) ICfatal(2);		// Insure we are not electron
  SISO      = Siso;				// Set isotropic value (PPM)
  SCSA      = Scsa;				// Set anisotropic value (PPM)
  SOMEGA    = Om;				// Set Larmor frequency (Hz)
  double Iz = II.qn();				// Get Iz value of spin			
  double X  = xi();				// Get interaction constant
  IntRank2::operator=(IntRank2(Iz,X,eta,EA));	// Use generic interaction
  }

IntCSA::IntCSA(const Isotope& IsoI,
           double Siso, double Scsa, double eta, const EAngles& EA, double Om)
  {
  if(!SpinCheck(IsoI, false)) ICfatal(2);	// Insure we are not electron
  SISO      = Siso;				// Set isotropic value (PPM)
  SCSA      = Scsa;				// Set anisotropic value (PPM)
  SOMEGA    = Om;				// Set Larmor frequency (Hz)
  double Iz = IsoI.qn();			// Get Iz value of spin			
  double X  = xi();				// Get interaction constant
  IntRank2::operator=(IntRank2(Iz,X,eta,EA));	// Use generic interaction
  }

IntCSA::IntCSA(double Iqn, 
           double Siso, double Scsa, double eta, const EAngles& EA, double Om)
  {
  if(!SpinCheck(Iqn)) ICfatal(2); 		// Insure valid quantum #
  SISO     = Siso;				// Store isotropic shift (PPM) 
  SCSA     = Scsa;				// Set anisotropic value (PPM)
  SOMEGA   = Om;				// Set Larmor frequency (Hz)
  double X = xi();				// Get interaction constant
  IntRank2::operator=(IntRank2(Iqn,X,eta,EA));	// Use generic interaction
  }

// ----------------------------------------------------------------------------
//            Direct Constructors That Use Cartesian Parameters
// ----------------------------------------------------------------------------

/* Here we try to use the set { Iqn, Sxx, Syy, Szz, SEAngles, SOm }. The spin
   quantum number sets the dimension of the interaction spin Hilbert space. 
   The values of Sxx, Syy, and Szz are in PPM and used to determine the 
   isotropic chemical shift SIso (PPM), chemical shift anisotropy SCSA (PPM),
   and shift asymmetry ETA ([0,1]). The interaction orientation is taken from
   Euler angles SEAngles and the external static field strength specified by 
   the spin Larmor frequency SOm (Hz). We will allow for default settings of
   the orientation and Larmor frequency so that those arguments need not be
   input during construction.  Their default values are 0.                   */

IntCSA::IntCSA(const string& IsoI, const coord& SxSySz,
                                                  const EAngles& EA, double Om)
  {
  if(!SpinCheck(IsoI)) ICfatal(2);		// Insure we know this isotope
  Isotope II(IsoI);				// Get isotope for this type
  if(!SpinCheck(II, false)) ICfatal(2);		// Insure we are not electron
  coord SiSzSe = AisoDelzEta(SxSySz);		// Switch to spherical values
  SISO      = SiSzSe.x();			// Set isotropic value (PPM)
  SCSA      = 1.5*SiSzSe.y();			// Set anisotropic value (PPM)
  SOMEGA    = Om;				// Set Larmor frequency (Hz)
  double Iz = II.qn();				// Get Iz value of spin			
  double E  = SiSzSe.z();			// Get asymmetry value
  double X  = xi();				// Get interaction constant
  IntRank2::operator=(IntRank2(Iz,X,E,EA));	// Use generic interaction
  }

IntCSA::IntCSA(const Isotope& IsoI, const coord& SxSySz,
                                                  const EAngles& EA, double Om)
  {
  if(!SpinCheck(IsoI, false)) ICfatal(2);	// Insure we are not electron
  coord SiSzSe = AisoDelzEta(SxSySz);		// Switch to spherical values
  SISO      = SiSzSe.x();			// Set isotropic value (PPM)
  SCSA      = 1.5*SiSzSe.y();			// Set anisotropic value (PPM)
  SOMEGA    = Om;				// Set Larmor frequency (Hz)
  double Iz = IsoI.qn();			// Get Iz value of spin			
  double E  = SiSzSe.z();			// Get asymmetry value
  double X  = xi();				// Get interaction constant
  IntRank2::operator=(IntRank2(Iz,X,E,EA));	// Use generic interaction
  }

IntCSA::IntCSA(double Iqn, const coord& SxSySz, const EAngles& EA, double Om)
  {
  if(!SpinCheck(Iqn)) ICfatal(2); 		// Insure valid quantum #
  coord SiSzSe = AisoDelzEta(SxSySz);		// Switch to spherical values
  SISO      = SiSzSe.x();			// Set isotropic value (PPM)
  SCSA      = 1.5*SiSzSe.y();			// Set anisotropic value (PPM)
  SOMEGA    = Om;				// Set Larmor frequency (Hz)
  double E  = SiSzSe.z();			// Get asymmetry value
  double X  = xi();				// Get interaction constant
  IntRank2::operator=(IntRank2(Iqn,X,E,EA));	// Use generic interaction
  }

// ----------------------------------------------------------------------------
//       Constructors Using Parameter Sets & A Spin/Interaction Index
// ----------------------------------------------------------------------------
 
/* These constructors use a single interaction (or spin) index. They'll try to
   read the parameters

             { [Iso(i)/Iqn(i)], CSA, Ceta, Calpha, Cbeta, Cgamma }

   where the shielding anisotropy is related to the spatial tensor delzz via

      ^          3                                     1 
     / \ sigma = - del   = sigma  - sigma   = sigma  - - [ sigma  + sigma  ]
     ---         2    zz        ||        |        zz  2        xx       yy
                                         ---

 NOTE: This interaction maintains itself in PPM units!  That is somewhat odd
       in light of the fact that the externally applied field strength may
       not have been specified nor do we have knowledge of spin isotope types.
       But, we do track the spin Larmor frequency which effectively takes care
       of this and allows a conversion from PPM to Hz.
 
   When idx!=-1 the parameter Iso(idx) will be used preferentially over the
   cooresponding parameter Iqn(idx) to set the interaction spin quantum value.
   The functions do allow for the suffix (#) but its not used in the default
   function call.  The functions do NOT allow for the prefix [#] since multiple
   interactions can be read using the suffix and multiple sets of interactions
   can be read using shift anisotropy interaction vectors (see class IntCSAVec)
   If neither parameter Iso or Iqn are present in the file/parameter set, a
   spin quantum number of 1/2 is used by default. Many parameters will be
   assumed zero if not present (e.g. asymmetry & orientation)                */

        // Input                SA	: CSA interaction (this)
        //                      pset    : Parameter set
	//			II	: Isotope type
        //                      idx     : Interaction index (default -1->none)
        //                      warn    : Warning output label
        //                                 0 = no warnings
        //                                 1 = warnings
        //                                >1 = fatal warnings
        // Output               none    : Shift anisotropy interaction
	//				  constructed from parameters in pset
	// Note				: Constructions that explicitly use
	//				  isotope labels are for support of
	//				  use by spin systems
 
IntCSA::IntCSA(const ParameterSet& pset, int idx, int warn)
  {
  if(!setCI(pset,idx,warn) && warn)		// Try & set interaction
    {						// and warn if troubles
    if(warn > 1) ICfatal(2); 
    else         ICerror(2,1);
    }
  }

IntCSA::IntCSA(const Isotope& II, const ParameterSet& pset, int idx, int warn)
  {
  if(!setCI(II,pset,idx,warn) && warn)		// Try & set interaction
    {						// and warn if troubles
    if(warn > 1) ICfatal(2); 
    else         ICerror(2,1);
    }
  }
                                                                                
// ----------------------------------------------------------------------------
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------

/* Both assignment and destruction are handled by the base class IntRank2.   */

void IntCSA::operator= (const IntCSA &S)  
  { 
  IntRank2::operator=(S); 		// Copyh irreducible rank 2 part
  SISO   = S.SISO;			// Copy the isotropy  (PPM)
  SCSA   = S.SCSA;			// Copy the anisotropy (PPM)
  SOMEGA = S.SOMEGA;			// Copy the Larmor frequency (Hz)
  }

IntCSA::~IntCSA () {}

// ____________________________________________________________________________
// B                     SPATIAL TENSOR COMPONENT ACCESS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//                              Isotropy Access
//-----------------------------------------------------------------------------
 
/* Here we only allow the chemical shift to be set in PPM. No assumptions are
   made regarding the external field strength, gyromagnetic ration. or spin
   Larmor frequency.                                                         */                

double IntCSA::iso() const     { return SISO; }
void   IntCSA::iso(double ppm) { SISO = ppm; }
double IntCSA::PPM() const     { return SISO; }
void   IntCSA::PPM(double ppm) { SISO = ppm; }

//-----------------------------------------------------------------------------
//                             Anisotropy Access
//-----------------------------------------------------------------------------
 
/*  The shielding anisotropy is related to the spatial tensor delzz via

      ^          3                                     1 
     / \ sigma = - del   = sigma  - sigma   = sigma  - - [ sigma  + sigma  ]
     ---         2    zz        ||        |        zz  2        xx       yy
                                         ---
    
    Since the shift tensor is stored in PPM units the returned anisotropy will
    be in PPM units as well. Note that the default base class IntRank2A also
    provides the GAMMA normalized anisotropy and delzz values. These are
    constants given by

                       1/2
                  [ 5 ]                  ^      3
          del   = |---]  = 0.51503      /_\ A = - del   = 0.77255
             zz   [6PI]                         2    zz                      */
 

// static double IntRank2A::delzz();                            INHERITED
// static double IntRank2A::delA();                             INHERITED

double IntCSA::aniso() const     { return SCSA; }
void   IntCSA::aniso(double csa) { SCSA = csa; }
double IntCSA::CSA() const       { return SCSA; }
void   IntCSA::CSA(double csa)   { SCSA = csa; }
 
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

complex IntRank2::A20()  const;					INHERITED
complex IntRank2::A21()  const;					INHERITED
complex IntRank2::A2m1() const;					INHERITED
complex IntRank2::A22()  const;					INHERITED
complex IntRank2::A2m2() const;					INHERITED

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
  double IntRank2A::Ayy() const;				INHERITED
  double IntRank2A::Azz() const;				INHERITED
  double IntRank2A::Ayx() const;				INHERITED
  double IntRank2A::Axy() const;				INHERITED
  double IntRank2A::Azx() const;				INHERITED
  double IntRank2A::Axz() const;				INHERITED
  double IntRank2A::Azy() const;				INHERITED
  double IntRank2A::Ayz() const;				INHERITED

  double IntRank2A::Axx(double A, double B, double G) const;	INHERITED
  double IntRank2A::Ayy(double A, double B, double G) const;	INHERITED
  double IntRank2A::Azz(double A, double B, double G) const;	INHERITED
  double IntRank2A::Ayx(double A, double B, double G) const;	INHERITED
  double IntRank2A::Axy(double A, double B, double G) const;	INHERITED
  double IntRank2A::Azx(double A, double B, double G) const;	INHERITED
  double IntRank2A::Axz(double A, double B, double G) const;	INHERITED
  double IntRank2A::Azy(double A, double B, double G) const;	INHERITED
  double IntRank2A::Ayz(double A, double B, double G) const;	INHERITED

  double IntRank2A::Axx(const EAngles& EA) const;		INHERITED
  double IntRank2A::Ayy(const EAngles& EA) const;		INHERITED
  double IntRank2A::Azz(const EAngles& EA) const;		INHERITED
  double IntRank2A::Ayx(const EAngles& EA) const;		INHERITED
  double IntRank2A::Axy(const EAngles& EA) const;		INHERITED
  double IntRank2A::Azx(const EAngles& EA) const;		INHERITED
  double IntRank2A::Axz(const EAngles& EA) const;		INHERITED
  double IntRank2A::Azy(const EAngles& EA) const;		INHERITED
  double IntRank2A::Ayz(const EAngles& EA) const;		INHERITED

  row_vector IntRank2A::CartComps() const;			INHERITED
  row_vector IntRank2A::CartComps(double theta, double phi=0) const;        */

//-----------------------------------------------------------------------------
//               Cartesian Tensor Component Access (Not Normalized)
//-----------------------------------------------------------------------------

/* These allow one to access the Cartesian elements of the full (unscaled)
   shift anisotropy spatial tensor at any orientation.

                                    1/2
                            [ 6*PI ]
                      S   = | ---- |    * del   * A   + S
                       uv   [  5   ]         zz    uv    iso                 */


double IntCSA::Sxx() const { return SISO + SCSA*RT6PIO5*Axx()/1.5; }
double IntCSA::Syy() const { return SISO + SCSA*RT6PIO5*Ayy()/1.5; }
double IntCSA::Szz() const { return SISO + SCSA*RT6PIO5*Azz()/1.5; }
double IntCSA::Sxy() const { return        SCSA*RT6PIO5*Axy()/1.5; }
double IntCSA::Syx() const { return        SCSA*RT6PIO5*Ayx()/1.5; }
double IntCSA::Sxz() const { return        SCSA*RT6PIO5*Axz()/1.5; }
double IntCSA::Szx() const { return        SCSA*RT6PIO5*Azx()/1.5; }
double IntCSA::Syz() const { return        SCSA*RT6PIO5*Ayz()/1.5; }
double IntCSA::Szy() const { return        SCSA*RT6PIO5*Azy()/1.5; }


double IntCSA::Sxx(double alpha, double beta, double gamma) const
                    { return SISO + SCSA*RT6PIO5*Axx(alpha,beta,gamma)/1.5; }
double IntCSA::Syy(double alpha, double beta, double gamma) const
                    { return SISO + SCSA*RT6PIO5*Ayy(alpha,beta,gamma)/1.5; }
double IntCSA::Szz(double alpha, double beta, double gamma) const
                    { return SISO + SCSA*RT6PIO5*Azz(alpha,beta,gamma)/1.5; }
double IntCSA::Syx(double alpha, double beta, double gamma) const
                    { return        SCSA*RT6PIO5*Ayx(alpha,beta,gamma)/1.5; }
double IntCSA::Sxy(double alpha, double beta, double gamma) const
                    { return        SCSA*RT6PIO5*Axy(alpha,beta,gamma)/1.5; }
double IntCSA::Szx(double alpha, double beta, double gamma)
                   const { return   SCSA*RT6PIO5*Azx(alpha,beta,gamma)/1.5; }
double IntCSA::Szy(double alpha, double beta, double gamma)
                   const { return   SCSA*RT6PIO5*Azy(alpha,beta,gamma)/1.5; }
double IntCSA::Sxz(double alpha, double beta, double gamma)
                   const { return   SCSA*RT6PIO5*Axz(alpha,beta,gamma)/1.5; }
double IntCSA::Syz(double alpha, double beta, double gamma)
                   const { return   SCSA*RT6PIO5*Ayz(alpha,beta,gamma)/1.5; }

double IntCSA::Sxx(const EAngles& EA) const {return SISO+SCSA*RT6PIO5*Axx(EA)/1.5;}
double IntCSA::Syy(const EAngles& EA) const {return SISO+SCSA*RT6PIO5*Ayy(EA)/1.5;}
double IntCSA::Szz(const EAngles& EA) const {return SISO+SCSA*RT6PIO5*Azz(EA)/1.5;}
double IntCSA::Syx(const EAngles& EA) const {return      SCSA*RT6PIO5*Ayx(EA)/1.5;}
double IntCSA::Sxy(const EAngles& EA) const {return      SCSA*RT6PIO5*Axy(EA)/1.5;}
double IntCSA::Szx(const EAngles& EA) const {return      SCSA*RT6PIO5*Azx(EA)/1.5;}
double IntCSA::Szy(const EAngles& EA) const {return      SCSA*RT6PIO5*Azy(EA)/1.5;}
double IntCSA::Sxz(const EAngles& EA) const {return      SCSA*RT6PIO5*Axz(EA)/1.5;}
double IntCSA::Syz(const EAngles& EA) const {return      SCSA*RT6PIO5*Ayz(EA)/1.5;}

  
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
double  IntRank2A::alpha()       const;				INHERITED
double  IntRank2A::beta()        const;				INHERITED
double  IntRank2A::gamma()       const;				INHERITED
double  IntRank2A::phi()         const;				INHERITED
double  IntRank2A::theta()       const;				INHERITED
EAngles IntRank2A::orientation() const;				INHERITED

void IntRank2A::alpha(double A);				INHERITED
void IntRank2A::beta(double  B);				INHERITED
void IntRank2A::gamma(double G);				INHERITED
void IntRank2A::phi(double   P);				INHERITED
void IntRank2A::theta(double T);				INHERITED
void IntRank2A::orientation(const EAngles& EA);			INHERITED
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

        1.) GAMMA normalized Auv - Done With IntCSA::CartMx()
        2.) Typical Suv values   - Done With IntCSA::Smx(true);
        3.) Shown in lab frame   - Done With This Function

   For case 3.) the values are related to the GAMMA normalized (Auv) and
   typically presented values (guv) according to

                            1/2
                    [ 6*PI ]
              S   = | ---- |    * del   * A   + Kdel    * S
               uv   [  5   ]         zz    uv       u,v    iso

   where Kdel is a Kronecker delta function.                                 */

matrix IntCSA::Smx() const
  {
  matrix SAmx = (SCSA*RT6PIO5/1.5)*CartMx(1);		// Scaled A matrix
  matrix dmx(3,3,complex1,d_matrix_type);               // GIso matrix
  return SAmx + SISO*dmx;                               // The G matrix
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
   of _XI to commonly used values that scale shift anisotropy interactions.

   Shift anisotropy interaction constants are returned in units of radian/sec.
   That presents a bit of a problem since these interactions are maintained in
   units of PPM.  The relationship is

                1/2                    1/2                 1/2
     SA   [6*pi]                 [6*pi]              [8*pi]          ^
   xi   = |----| gamma*B *del  = |----| Omega*del  = |----| Omega * / \ sigma
          [ 5  ]        o    zz  [ 5  ]          zz  [ 15 ]         --- 

   where gamma is the spin gyromagnetic ratio, B an applied field strength and
   Omega a spin Larmor frequency.  To obtain the interaction constant in 
   angular frequency units we must convert using the Larmor frequency of the
   spin involved in the interaction. 
 
           Input                SA	: Shift anisotropy interaction
	  			IsoI	: Spin isotope type
	  			Om      : Spin Larmor frequency (Hz)
           Return               xi      : CSA. interaction constant
           Note                           Value is in radians/sec            */
 
double IntCSA::xiOm(const string& IsoI, double Om) const
  { 
  Isotope I(IsoI);				// Construct isotope for IsoI
  Isotope H("1H");				// This is a proton
  double IOm = Om*I.gamma()/H.gamma();		// Larmor frequency I spins
  return xiOm(IOm);				// Use overloaded function
  }						// Note: Om is Larmor of IsoI

double IntCSA::xiOm(double Om) const		// Here Om is the Larmor
  { return PIx2*sqrt(8.*PI/15.)*Om*1.e-6*SCSA; }// frequency of our spin (Hz)

double IntCSA::xi() const			// Overwrites inherited fct
  { return SOMEGA?xiOm(SOMEGA):0.0; }

void IntCSA::xi(double X) const
  {
  ICerror(60,1);				// Cannot set Xi
  ICerror(61,1);				// Xi depends on Omega & CSA
  ICfatal(62);					// Use function xiOm or CSA
  X = 0.0;					// Use to avoid complaints
  }						// from the compiler
 
// ____________________________________________________________________________
// E             SHIFT ANISOTROPY INTERACTION AUXILIARY FUNCTIONS
// ____________________________________________________________________________

 
// ____________________________________________________________________________
// F                        PARAMETER SET FUNCTIONS
// ____________________________________________________________________________
 
// ----------------------------------------------------------------------------
//  Functions To Make A Parameter Set From A Shielding Anisotropy Interaction
// ----------------------------------------------------------------------------

/* This class has no implicit knowledge of the nucleus type. The interaction
   may involve a 1H (I=1/2) or a 2H (I=1) or any GAMMA recognized nuclear spin.
   The preferred means of specifying a shielding anisotropy interaction when
   there is no knowledge of the nulcear spin isotope type is that which is used
   for filling up the parameter set.  The base parameters of individual
   interactions in that case are { CI(#), CSA(#), CTheta(#), CPhi(#), CEta(#) }.

          Input                SA      : Shift anisotropy interaction
                               idx     : Interaction index (default -1)
                               pfx     : Interaction 2nd indx (def -1)
          Output               pset    : Parameter set with only shielding
                                         anisotropy interaction parameters

    Function                                 Purpose
  ------------         -------------------------------------------------------
  ParameterSet         Convert interaction into a parameter set
  operator +=          Adds interaction to existing parameter set
  PSetAdd              Adds interaction to existing parameter set, this
                       allows for an index and a prefix in the parameters   */

IntCSA::operator ParameterSet( ) const
  { ParameterSet pset; pset += *this; return pset; } 

void operator+= (ParameterSet& pset, const IntCSA &SA) { SA.PSetAdd(pset); }

void IntCSA::PSetAdd(ParameterSet& pset, int idx, int pfx) const
  {                                                            
  string suffx;                                 // Parameter suffix
  if(idx != -1)                                 // Only use suffix if idx
    suffx = string("(")+Gdec(idx)+string(")");	// is NOT -1
  string prefx;                                 // Parameter prefix
  if(pfx != -1)                                 // Only use prefix if pfx
    prefx = string("[")+Gdec(pfx)+string("]");	// is NOT -1
 
  string pname  = prefx +string("Iqn") + suffx;	// Add spin quantum number
  string pstate = string("Spin Quantum Number");
  double pdatad = Izval();
  SinglePar par = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);
 
  pname  = prefx + string("PPM") + suffx;	// Add chemical shift (PPM)
  pstate = string("Chemical Shift (PPM)");
  pdatad = SISO;
  par    = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);
 
  pname  = prefx + string("CSA") + suffx;	// Add shift anisotropy (PPM)
  pstate = string("Shift Anisotropy (PPM)");
  pdatad = SCSA;
  par    = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);
 
  pname  = prefx + string("Seta") + suffx;	// Add asymmetry
  pstate = string("Shift Asymmetry [0,1]");
  pdatad = ETA;
  par    = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);

  pname  = prefx + string("SEAngles") + suffx;	// Add orientation (degrees)
  pstate = string("Shift Euler Angles (deg)");
  double a = _EAs.alpha() * RAD2DEG;
  double b = _EAs.beta()  * RAD2DEG;
  double g = _EAs.gamma() * RAD2DEG;
  coord EA(a, b, g);
  par = EA.param(pname, pstate);
  pset.push_back(par);

  if(SOMEGA)
    {
    pname  = prefx + string("Omega") + suffx;	// Add Larmor frequency (MHz)
    pstate = string("Larmor Frequency (MHz)");
    pdatad = SOMEGA;
    par    = SinglePar(pname, pdatad, pstate);
    pset.push_back(par);
    }
  }

// ----------------------------------------------------------------------------
// Functions To Output Shift Anisot. Interaction To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

/* These functions write the CSA interaction into an ASCII file in GAMMA
   Parameter Set format.  That is, the resulting ASCII file may be read into
   any GAMMA program to create a CSA interaction identical to that written. 

        // Input                SA	: Shift anisotropy interaction
        //                      filename: Output file name
        //                  or  ofstr   : Output file stream
        //                      idx     : Interaction index (default -1)
        //                      pfx     : Interaction 2nd indx (def -1)
        //                      warn    : Warning level
        // Output               none    : Interaction is written as a
        //                                parameter set to file filename
	//				  or into output file stream ofstr   */


int IntCSA::write(const string &filename, int idx, int pfx, int warn) const

  {
  ofstream ofstr(filename.c_str());	// Open filename for input
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!write(ofstr, idx, pfx, w2))       // If file bad then exit
    {
    ICerror(1, filename, 1);		// Filename problems
    if(warn>1) ICfatal(20);		// Fatal error
    return 0;
    }
  ofstr.close();                        // Close it now
  return 1;
  }
 
int IntCSA::write(ofstream& ofstr, int idx, int pfx, int warn) const
  {
  ParameterSet pset;			// Declare a parameter set
  PSetAdd(pset, idx, pfx);              // Add in interaction parameters
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!pset.write(ofstr, w2))            // Use parameter set to write
    {					// out the interaction parameters
    if(warn)
      {
      ICerror(52, 1);                   // Problems writing to filestream
      if (warn>1) ICfatal(53);       // Fatal error
      }
    return 0;
    }  
  return 1;
  }  

// ____________________________________________________________________________
// G              SHIFT ANISOTROPY INTERACTION INPUT FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//      Direct Read of SA Intraction From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

/* These next two read functions utilize a single spin/interaction index, for
   the spin involved in the interaction.  They'll try to read the parameter set

           { [Iso(i)/Iqn(i)] , CSA, Ceta, Calpha, Cbeta, Cgamma, Omega }

   The functions do allow for the suffix (#), but do NOT allow for the prefix
   [#] because multiple shift anisotropy interactions can be defined in the
   same file by switching spin/interaction indices.  Users can define multiple
   sets of interactions in the same file by using SA interaction vectors
   (see class IntCSAVec.)                                                    */
 
	// Input		SA	: CSA interaction (this)
        //                      filename: Input filename
        //                      pset    : Parameter set
        //                      idx     : Parameter index value used for
        //                                suffix (#) in input names
        //                      warn    : Warning output label
        //                                 0 = no warnings
        //                                 1 = warnings
        //                                >1 = fatal warnings
        // Output               none    : CSA interaction filled with
        //                                parameters read from file
	//				  or from parameter set
	// Note				: This function supercedes
	//				  that in IntRank2A

bool IntCSA::read(const string &filename, int idx, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename,warn?1:0))     // Read in pset from file
    {
    if(warn)
      {
      ICerror(1, filename, 1);		// Filename problems
      if(warn>1) ICfatal(21);		//      Its a fatal error
      else       ICerror(21,1);		//      or perhaps not
      }
    return false;
    }
  return read(pset, idx, warn);		// Use overload function
  }


bool IntCSA::read(const ParameterSet& pset, int idx, int warn)

  {
  bool TF = setCI(pset, idx, warn?1:0);
  if(!TF && warn)
    {
    string val;
    if(idx == -1) val = string(" None"); 
    else          val = Gdec(idx);
    if(warn > 1) ICfatal(44, val);         // Fatal error
    else         ICerror(44, val);         // or a warning issued
    }
  return TF;
  }
 

// ----------------------------------------------------------------------------
//  Interactive Reading of SA Intractions From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

	// Input		SA	: CSA interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
        //                      idx     : Interaction index
        // Output               string  : The parameter argn of array argc is
        //                                used to supply a filename from which
        //                                the interaction is read
 
string IntCSA::ask_read(int argc, char* argv[], int argn, int idx)
  {
  string filename;                              // Name of parameter file
  query_parameter(argc, argv, argn,             // Get filename from command
  "\n\tShift Anisotropy Interaction filename? ",// Or ask for it
                                     filename);
  read(filename, idx);                          // Read system from filename
  return filename;
  }

// ---------------------------------------------------------------------------- 
//               Interactive Ask For All Kinds Interaction Info 
// ---------------------------------------------------------------------------- 

void IntCSA::ask(int argc, char* argv[], int& qn, double& CI,
            double& Csa, double& Ceta, double& Ctheta, double& Cphi, int Cflag)
 
        // Input                SA      : CSA interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Array of arguments
        //                      qn      : Query value
        //                      CI      : Spin quantum number
        //                      Csa     : Shift anisotropy (PPM)
        //                      Ceta    : SA asymmetry
        //                      Ctheta  : SA orientation angle
        //                      Cphi    : SA orientation angle
        //                      Cflag   : Flag is eta is requested
        // Output               none    : The values of qn, I, Csa, Ceta,
        //                                Ctheta, and Cphi are set herein
        // Note                         : This is INTERACTIVE!
        // Note                         : This does NOT construct SA!

  {                                                                   
  query_parameter(argc, argv, qn++,                     // Read in the I value
    "\n\tI Spin Quantum Value (0.5, 1, 1.5, ..)? ", CI);
  query_parameter(argc, argv, qn++,                     // Read the anisotropy
                   "\n\tShift Anisotropy(PPM)? ", Csa);
  query_parameter(argc, argv, qn++,                     // Read in theta angle
  "\n\tOrientation Down From PAS z [0, 180]? ", Ctheta);
  query_parameter(argc, argv, qn++,                     // Read in phi angle
   "\n\tOrientation Over From PAS x [0, 360]? ", Cphi);
  if(Cflag)
    query_parameter(argc, argv, qn++,                   // Read in asymmetry
           "\n\tInteraction Asymmetry [0, 1]? ", Ceta);
  else Ceta = 0;
  }

  
void IntCSA::askset(int argc, char* argv[], int& qn, int SAflag)

        // Input                SA	: CSA interaction (this)
	//			argc    : Number of arguments
	//			argv    : Array of arguments
	//			qn      : Query value
        //                      SAflag	: Flag if asymmetry CSA requested
        // Output               none    : SA is set interactively
	// Note				: This is INTERACTIVE!

  {
  double SAI=0, SAdelz=0, SAeta=0, SAtheta=0, SAphi=0;
  ask(argc,argv,qn,SAI,SAdelz,SAeta,SAtheta, SAphi, SAflag);
  *(this) = IntCSA(SAI,SAdelz,SAeta,SAtheta,SAphi);
  }

 
void IntCSA::askset(int SAflag)
 
        // Input                SA	: CSA interaction (this)
        //                      SAflag	: Flag if asymmetry CSA requested
        // Output               none    : SA is set interactively
        // Note                         : This is INTERACTIVE!
 
  {
  int qn = 1000;
  int argc = 1;
  char* argv[1];
  double SAI, SAdelz, SAeta, SAtheta, SAphi;
  ask(argc,argv,qn,SAI,SAdelz,SAeta,SAtheta, SAphi,SAflag);
  *(this) = IntCSA(SAI,SAdelz,SAeta,SAtheta,SAphi);
  }


// ____________________________________________________________________________
// I                          OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//   Functions That Generate Simple Strings To Simplify & Modularize Printing
//-----------------------------------------------------------------------------

string IntCSA::ShiftString() const
  {
  string Ssft("Isotropic Shift  (PPM) ");
  Ssft += Gform("%11.7f", SISO);
  return Ssft;
  }

string IntCSA::CSAString() const
  {
  string Ssft("Shift Anisotropy (PPM)");
  Ssft += Gform("%12.7f", SCSA);
  return Ssft;
  }

string IntCSA::LarmorString() const
  {
  string Ssft("Larmor Frequency (MHz)");
  Ssft += Gform("%12.7f", SOMEGA);
  return Ssft;
  }


//-----------------------------------------------------------------------------
// Functions To Generate Information Strings To Simplify & Modularize Printing
//-----------------------------------------------------------------------------

/* These functions return spatial tensor information in string format. This is
   done to facilitate printing, in particular printing of CSA spatial tensors
   from within rank 2 interactions.  In this case we make a list of information
   thats usually displayed to the left of a 3x3 Cartesian matrix rep. of S.

                       Spin Quantum Number:       I
    	               Isotropic shift        xxxxx.xx PPM
                       Shift Anisotrpy:       xxxxx.xx PPM
                       Asymmetry:                 x.xx
                       Euler Angle Alpha:       xxx.xx Degrees
                       Euler Angle Beta:        xxx.xx Degrees
                       Euler Angle Gamma:       xxx.xx Degrees               */

vector<string> IntCSA::InfoStrings() const
 {
  vector<string> SphStrings;
  SphStrings.push_back(StringI());			// Spin quantum #
  SphStrings.push_back(ShiftString());			// Isotropic shift
  SphStrings.push_back(CSAString());			// Shift anisotropy
  SphStrings.push_back(AsymmetryString());		// Asymmetry
  if(PAS())						// PAS orientation
    SphStrings.push_back(PASString());			// if in PAS
  else							// Non PAS orientation
    {							// we use 3 angles
    SphStrings.push_back(AlphaString());		// Alpha
    SphStrings.push_back(BetaString());			// Beta
    SphStrings.push_back(GammaString());		// Gamma
    }
  SphStrings.push_back(LarmorString());
  SphStrings.push_back(XiString());
  return SphStrings;
  }

//-----------------------------------------------------------------------------
//  Functions That Generate String Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------


/* String* IntRank2T::TStrings(int M) const;                     INHERITED    

                                [ x.x, x.x, x.x]
                        T     = [ x.x, x.x, x.x]
                         2,m    [ x.x, x.x, x.x]
                                [ x.x, x.x, x.x]

  where M = { 0, 1, ..., 4 } ===> m = { 0, 1, -1, 2, -2 }                    */  
     
//-----------------------------------------------------------------------------
//   Functions That Generate Ouput Of The Rank 2 Shift Anistropy Interaction
//-----------------------------------------------------------------------------

        // Input                SA      : CSA interaction (this)
        //                      ostr    : Output stream
        //                      fflag   : Format flag
        //                                  0 - Basic Parameters
        //                                 !0 - Full output
        // Output               none    : CSA spatial tensor parameters 
        //                                placed into the output stream
        // Note                         : This does NOT use the base class
        //                                virtual overload because we write
        //                                out two spin quantum values here?
 
ostream& IntCSA::print(ostream& ostr, int fflag) const
  {
  if(Izval() < 0.5)                             // Just exit if nothing
    {
    string hdr("Empty Shift Anisotropy Interaction");
    int hl = hdr.length();
    string spacer = string(40-hl/2, ' ');
    ostr << "\n\n" << spacer << hdr << "\n";
    return ostr;
    }

/*          Output Some Printed Lines Which Will Look Like The Following

                               Shift Anisotropy Interaction

                                          [S  , S  , S  ]
   Spin Quantum Number:       x.xxxxxx    [ xx   xy   xz]   [ x.x, x.x, x.x]
   SCC (kHz):                 x.xxxxxx    [S  , S  , S  ] = [ x.x, x.x, x.x]
   Xi Value (rad/sec):        x.xxxxxx    [ yx   yy   yz]   [ x.x, x.x, x.x]
                                          [S  , S  , S  ]
                                          [ zx   zy   zz]                    */

  vector<string> SS = InfoStrings(); 		// Array of info strings
  string hdr="Shift Anisotropy Interaction";	// Use this header
  IntRank2A::print(ostr, hdr, SS);		// Print spatial tensor
  if(!fflag) return ostr;

/*          Now Output The Spherical Components, Spatial And Spin
       These Will Be Printed Lines Which Will Look Like The Following
                        (Repeated For All 5 m Values)

                                             [ x.xx  x.xx  x.xx  x.xx]
             A    = x.xxx            T     = [ x.xx  x.xx  x.xx  x.xx]
              2,m                     2,m    [ x.xx  x.xx  x.xx  x.xx]
                                             [ x.xx  x.xx  x.xx  x.xx]      */

  printAT(ostr); 				// Print sphercial A & T
  ostr << "\n\n";                               // Add some linefeeds
  return ostr;
  }

ostream& operator<< (ostream& out, const IntCSA& SA) { return SA.print(out); }

 
ostream& IntCSA::printAT(ostream& ostr) const

        // Input                SA      : CSA interaction (this)
        //                      ostr    : Output stream
        // Output               none    : CSA spatial tensor parameters 
        //                                placed into the output stream
        // Note                         : Uses base class virtual overload
 
 
/*          Prints The Spherical Components, Spatial And Spin
       These Will Be Printed Lines Which Will Look Like The Following
                      (Repeated For All 5 m Values)

                                                [ x.x, x.x, x.x]
                A    = x.xxx            T     = [ x.x, x.x, x.x]
                 2,m                     2,m    [ x.x, x.x, x.x]
                                                [ x.x, x.x, x.x]             */  

  { return IntRank2::printAT(ostr, SA); }


ostream& IntCSA::STList(ostream& ostr, int fflag)
  {
  string hdr("Shift Anisotropy Spin Tensors List");
  int hl = hdr.length();
  string spacer = string(40-hl/2, ' ');
  ostr << "\n" << spacer << hdr << "\n";
  IntRank2::printList(ostr, SA, 1);
  return ostr;
  }


// ____________________________________________________________________________
// M                    SHIFT ANISOTROPY HAMILTONIAN FUNCTIONS
// ____________________________________________________________________________

/* This section returns irreducible rank 2 shift anisotropy interaction
   Hamiltonians. Because this class uses standardized spatial and spin tensors
   the Hamiltonians may all be derived from the same formula.

                        -2,2
                         ---           m
  H(alpha,beta,gamma) =  \   Xi  * (-1)  * A   (alpha,beta,gamma) * T    (i,j)
                         /                  2,m                      2,-m
                         ---
                          m

   The Hamiltonians will be returned in the single spin Hilbert space of the
   interaction, unless a composite Hilbert space and some spin index is 
   supplied. All Hamiltonians will returned in the product basis as simple 
   matrices. Their units will be radians/sec.				     */

// ----------------------------------------------------------------------------
//               First Order Shift Anisotropy Interaction Hamiltonian
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
// Applicable When The Shift Anisotropy Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------

/*  The secular part of the shift anisotropy Hamiltonian is that which returns
    only those components which commute with z axis rotations.  Here we have

     [1]                       SA    SA                       SA     (0)
    H   (alpha,beta,gamma) = Xi   * A   (alpha,beta,gamma) * T    = H
     SA                              2,0                      2,0    SA


    Note that the interaciton constant is field dependent! The overall strength
    of the interaction very much depends upon the Larmor frequency (external
    field strength) and if this is not set the return Hamiltonian will be 0.

           Input                SA      : Shift anisotropy interaction
                                alpha   : Euler angle (radians)
                                beta    : Euler angle (radians)
                                gamma   : Euler angle (radians)
                                EA      : Euler angles (radians)
                                HSs     : Array of spin Hilbert spaces
                                i       : Spin index (in HSs)
           Output               H0      : The secular part of the shift
                                          anisotropy Hamiltonian
                                             (default basis, Hz)
           Note                         : Also called the 1st order SA
                                          interaction (perturbation theory)
           Note                         : Rotationally invariant about z
           Note                         : No HSs: return in the spin Hilbert
                                          space of dimension (2I+1)
                                          With HSs: return in composite
                                          spin space                         */

matrix IntCSA::H0() const
  { return IntRank2::H0(); }
matrix IntCSA::H0(double A, double B, double G) const
  { return IntRank2::H0(A,B,G); }
matrix IntCSA::H0(const EAngles& EA) const
  { return IntRank2::H0(EA); }

matrix IntCSA::H0(const vector<int>& HSs, int i) const
  { return IntRank2::H0(HSs,i); }
matrix IntCSA::H0(const vector<int>& HSs, int i, double A, double B, double G) const
  { return IntRank2::H0(HSs,i,A,B,G); }
matrix IntCSA::H0(const vector<int>& HSs, int i, const EAngles& EA) const
  { return IntRank2::H0(HSs,i,EA); }

// ----------------------------------------------------------------------------
//                      Full CSA Interaction Hamiltonians
//       These Have No Approximations & Are Independent of Field Strength
//     (They Are Correct in The Laboratory Frame, NOT in a Rotating Frame) 
// ----------------------------------------------------------------------------

/* Here we make no approximations. As such, these Hamiltonians are not those
   used in the rotating frame. Here we simply sum over the 5 spatial and spin
   tensor components to produce the Hamiltonian.

                         -2,2
                          ---   SA      m    SA                       SA
 H  (alpha,beta,gamma) =  \   Xi  * (-1)  * A   (alpha,beta,gamma) * T    (i,j)
  SA                      /                  2,m                      2,-m
                          ---
                           m

	   Input		SA	: CSA interaction
                                alpha   : Euler angle (radians)
                                beta    : Euler angle (radians)
                                gamma   : Euler angle (radians)
                                EA      : Euler angles (radians)
                                HSs     : Array of spin Hilbert spaces
                                i       : First spin index (in HSs)
                                j       : Seond spin index (in HSs)
           Output               H       : Matrix for dipolar Hamiltonian
           Note                         : No HSs: return in the spin Hilbert
                                          space of dimension (2I+1)
                                          With HSs: return in composite spin
                                          Hilbert space 

   => Users should set the external field strenth prior to this call.
   => There are no m = +/-2 spin components in rank 2 CSA interactions.      
   => We have our own functions to convert from radians/sec to cycles/sec    */

// matrix IntRank2::H( ) const                                        INHERITED
// matrix IntRank2::H(double alpha, double beta, double gamma) const  INHERITED
// matrix IntRank2::H(const EAngles& EA) const                        INHERITED

// matrix IntRank2::H(const vector<int>& HSs, int i) const            INHERITED
// matrix IntRank2::H(const vector<int>& HSs, int i,                  INHERITED
//                     double alpha, double beta, double gamma) const
// matrix IntRank2::H(const vector<int>& HSs, int i,                  INHERITED
//                                           const EAngles& EA) const














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


 
ostream& IntCSA::printSpherical(ostream& ostr)
         
        // Input                SA	: CSA spatial tensor (this)
        //                      ostr	: Output stream
        // Output               none    : CSA spatial tensor parameters 
	//				  sent to output stream
	// Note				: Uses base class virtual overload

  {
  ostr << "\t\n\t         CSA Spatial Tensor";
  IntRank2A::printSpherical(ostr, 0);
  return ostr; 
  }

 
ostream& IntCSA::printCartesian(ostream& ostr)
   
        // Input                SA	: CSA spatial tensor (this)
        //                      ostr    : Output stream
        // Output               none    : CSA spatial tensor parameters 
	//				  sent to output stream
	// Note				: Uses base class virtual overload

// Axx = (A22+A2m2)/2 - A20/sqrt(6)      Axy = -i*(A22-A21)/2 = Ayx
// Ayy = -(A22+A2m2)/2 - A20/sqrt(6)     Axz = -(A21-A2m1)/2  = Azx
// Azz = sqrt(2/3)*A20                   Ayz = i*(A21-A2m1)/2 = Azy
 
  {
  ostr << "\t\n\t         CSA Spatial Tensor";
//  IntRank2A::printCartesian(ostr, 0);
  return ostr;
  }

 
ostream& IntCSA::printCartesian(ostream& ostr, double phi, double theta)
         
        // Input                SA	: CSA spatial tensor (this)
        //                      ostr	: Output stream
        //                      phi     : Orientation angle (degrees)
        //                      theta   : Orientation angle (degrees)
        // Output               none    : CSA spatial tensor parameters 
	//				  sent to output stream
	// Note				: Uses base class virtual overload

// Axx = (A22+A2m2)/2 - A20/sqrt(6)      Axy = -i*(A22-A21)/2 = Ayx
// Ayy = -(A22+A2m2)/2 - A20/sqrt(6)     Axz = -(A21-A2m1)/2  = Azx
// Azz = sqrt(2/3)*A20                   Ayz = i*(A21-A2m1)/2 = Azy
 
  {
  ostr << "\t\n\t    CSA Spatial Tensor";
//  IntRank2A::printCartesian(ostr, theta, phi, 0);
  return ostr; 
  }









// ____________________________________________________________________________
// D                      SHIFT ANISOTROPY FREQUENCY FUNCTIONS
// ____________________________________________________________________________
 
// Note that the default CSA frequency is defined to be related to the
// nuclear shift anisotropy constant (NQCC) by the relationship
 
//                                              2
//                               3*NQCC       3e qQ
//                        W   = --------  = --------
//                         Q    2I(2I-1)    2I(2I-1)
 
// According to the above formula, this frequency will be the splitting between
// the 2*I transitions contained in a CSA Hamiltonian.  This will be
// observed at zero field or when the shift anisotropy is weak relative to
// the Zeeman interaction (high field approximation, first order terms only)
// if the tensor is oriented in it's principal axes (PAS).  If the tensor isn't
// aligned in it's PAS, the splitting will vary with orientation according to

//                     1           2                2                        
//    W (theta,phi)  = - W  [ 3*cos (theta) - eta*sin (theta)cos(2*phi) ]
//     Q               2  Q

// where theta is the angle down from the PAS z axis and phi the angle over
// from the PAS x axis.  Alternatively, if the Zeeman terms aint much stronger
// the splittings won't be equally spaced at all (second order terms).
 
// Also, keep in mind that the Euler angles {phi,theta,gamma} which are kept
// with the tensor are used to relate the tensor PAS to some (common) set of
// coordinate axes.  They are not the phi and theta used in the above
// formula (unless you which the tensor aligned in the common axis system)


double IntCSA::wC( ) const

        // Input                SA	: CSA interaction
        // Return               wQ      : CSA frequency (Hz)

  { return SCSA/1.5; }


void IntCSA::wC(double W)

        // Input                Q	: CSA interaction
	//			W       : CSA frequency (Hz)
        // Return               void	: The CSA frequency 
	//			          of the interaction is set
	// Note				: The interaction I value
	//				  must be set prior to this 

  {
  SCSA = W;
// sosi
  }






// ____________________________________________________________________________
// N                 SHIFT ANISOTROPY HAMILTONIAN FRIEND FUNCTIONS
// ____________________________________________________________________________

// These Don't Depend on Class IntCSA Variables At All.  The Functions Just
// Mimic The Class' Hamiltonian Generation.


matrix HC0(double qn, double wQo, double eta, double theta, double phi)
 
	// Input		qn	: Quantum number (1, 1.5, 2.5,...)
	//			wQo     : PAS CSA frequency
	//			eta     : CSA asymmetry [0,1]
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               H0	: The secular part of the quadrupolar
        //                                Hamiltonian (default basis, Hz)
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note				: Also called the 1st order quadrupolar
        //                                interaction (perturbation theory)
        // Note                      	: Rotationally invariant about z
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1
 
//  The secular part of the quadrupolar Hamiltonian is that which returns
//  only those components which commute with z axis rotations.  Here we have
 
//            [1]             1                   2             (0)
//           H  (theta,phi) = - w (the,phi) * [3*I - I(I+1)] = H
//            Q               6  Q                z             Q
 
// where 

//                      [ 1      2               1        2                   ]
// w (theta,phi) = W    | - [3cos (theta) - 1] + - eta sin (theta)*cos(2*phi) |
//  Q               Q,o [ 2                      2                            ]

// and                                                                          
//                                    3*QCC
//                            w    = --------
//                             Q,o   2I(2I-1)
 
  {
  int Ival = int(2.*qn + 1);			// For 1 spin SOp functions
  matrix IE = Ie(Ival);				// The identity operator
  matrix IZ = Iz(Ival);				// The Fz operator
  double Ctheta = cos(theta*DEG2RAD);		// Need cos(theta)
  double Stheta = sin(theta*DEG2RAD);		// Need sin(theta)
  double C2phi = cos(2.0*phi*DEG2RAD);		// Need cos(2phi)
  double wQ = wQo * ((3.0*Ctheta*Ctheta-1.0)	// "Oriented" wQ
            + eta*Stheta*Stheta*C2phi);
  return (wQ/12.0)*(3*(IZ*IZ) - (qn*(qn+1))*IE);// 1st Order Hamiltonian
  }  


matrix HC1(double Om, double qn, double wQo, double eta,
                                                      double theta, double phi)
 
        // Input                Om	: Field Strength (Larmor in Hz)
	// 			qn	: Quantum number (1, 1.5, 2.5,...)
	//			wQo     : PAS CSA frequency
	//			eta     : CSA asymmetry [0,1]
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HC1	: The 2nd order secular part of the
        //                                quadrupolar Hamiltonian (default
        //                                basis, Hz) for the interaction
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note				: Also called the 1st order quadrupolar
        //                                interaction (perturbation theory)
        // Note                      	: Rotationally invariant about z
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1
	// Note				: This will be zero in PAS if no eta
	//				  and small if Om >> QCC

//  [2]              -1   [ 1    ]2 [                          2     2   
// H   (theta,phi) = -- * | - w  |  | 2*V V  (theta,phi)*I *(4*I - 8I  - 1)
//  Q                Om   [ 6  Q ]  [    1 -1             z          z
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
 
//  in accordance with the article by P.P. Man "CSA Interactions" in
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
  Hmx       *= (-wQo*wQo/(36.0*Om));
  return Hmx;
  }
 
// sosi





#endif						// IntCSA.cc
