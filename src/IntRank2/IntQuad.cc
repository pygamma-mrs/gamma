/* IntQuad.cc ***************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Spatial Quadrupolar Tensor 		Implementation		**
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
** A variable of type IntQuad represents a quadrupolar interaction      **
** defined for a nuclear spin with a specific quantum number I.		**
** The interaction contains the essence of a scaled rank 2 spatial      **
** tensor which is restricted to have the following properties:         **
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
** class, this IntQuad provides functions for building up the tensor    **
** and accessing the tensor elements from a "quadrupolar" standpoint.   **
**                                                                      **
** The following defintions are used herein (Auv GAMMA normlized):	**
**                                                                      **
** 1.) PAS: |Azz| >= |Ayy| >= |Axx|					**
**                                                                      **
** 2.) PAS: Azz=2C, eta=(Axx-Ayy)/Azz, Axx=C(eta-1), Ayy=-C(1+eta)	**
**            								**
**                                    1/2				**          
**                           [   5   ]					**
**     where             C = | ----- |					**
**                           [ 24*PI ]					**
**            								**
**                                   2					**
** 3.) Quadrupolar Coupling: NQCC = e qQ				**
**                                                 2			**
**                                   3*NQCC      3e qQ			**
** 4.) Quadrupolar Frequency: w   = --------  = --------		**
**                             Q    2I(2I-1)    2I(2I-1)		**
**                                                                      **
*************************************************************************/

#ifndef   IntQuad_cc_			// Is file already included?
#  define IntQuad_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <IntRank2/IntQuad.h>		// Include interface definition
#include <Basics/Gconstants.h>		// Include PI and other constants
#include <Basics/SinglePar.h>		// Include parameter sets
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/Isotope.h>		// Include isotopes
#include <Basics/Gutils.h>		// Include parameter queries
#include <Matrix/matrix.h>		// Include matrices
#include <Matrix/row_vector.h>		// Include row_vectors
#include <HSLib/SpinOpSng.h>		// Include 1 spin operators
#include <HSLib/SpinSys.h>		// Include base spin systems
#include <Basics/StringCut.h>
#include <stdlib.h>
#include <vector>			// Inlcude libstdc++ STL vectors
#include <string>			// Inlcude libstdc++ STL strings
#if defined(_MSC_VER) || defined(__SUNPRO_CC)				// If we are using MSVC++
 #pragma warning (disable : 4786)       //   Kill STL namelength warnings
#endif

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ STL lists
using std::cout;			// Using libstdc++ standard output
using std::vector;			// Using libstdc++ STL vectors
using std::ofstream;			// Using libstdc++ output file streams
using std::ostream;			// Using libstdc++ output streams

const double xRT6PIO5 = sqrt(6.*PI/5.);		// Constant sqrt[6*PI/5]
const double xRT5O4PI = sqrt(5./(4.*PI)); 	// Constant sqrt[5/(4*PI)];
const double xRT5O24PI = sqrt(5./(24.*PI)); 	// Constant sqrt[5/(24*PI)];
const double xRT5O96PI = 0.5*xRT5O24PI; 	// Constant sqrt[5/(96*PI)];

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i             CLASS QUADRUPOLAR INTERACTION ERROR HANDLING
// ____________________________________________________________________________
 
/*       Input 		      Q	      : Quadrupolar interaction (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void IntQuad::Qerror(int eidx, int noret) const
  {
  string hdr("Quadrupolar Interaction");
  switch (eidx)
    {
    case 0: GAMMAerror(hdr,"Program Aborting.....",        noret); break;// (0)
    case 1: GAMMAerror(hdr,"Construction From Rank!=2.",   noret); break;// (1)
    case 2: GAMMAerror(hdr,"Problems During Construction.",noret); break;// (2)
    case 3: GAMMAerror(hdr,"Problems During Assignment.",  noret); break;// (3)
    case 8: GAMMAerror(hdr,"Theta (z Down) Beyond [0,180]",noret); break;// (8)
    case 9: GAMMAerror(hdr,"Phi (x Over) Outside [0, 360]",noret); break;// (9)
    case 10:GAMMAerror(hdr,"Asymmetry (eta) Beyond [0, 1]",noret); break;// (10)
    case 11:GAMMAerror(hdr,"I<=1/2 Is Not Allowed",        noret); break;// (11)
    case 12:GAMMAerror(hdr,"Set Asymmetry On Zero Tensor", noret); break;// (12)
    case 13:GAMMAerror(hdr,"Cannot Construct From File",   noret); break;// (13)
    case 18:GAMMAerror(hdr,"Cannot Alter Interaction",     noret); break;// (18)
    case 19:GAMMAerror(hdr,"No I Value Specified",         noret); break;// (19)
    case 20:GAMMAerror(hdr,"Can't Write To Parameter File",noret); break;// (20)
    case 21:GAMMAerror(hdr,"Cant Read From Parameter File",noret); break;// (21)
    case 25:GAMMAerror(hdr,"Spin Is An Electron",          noret); break;// (25)
    case 40:GAMMAerror(hdr,"Deprecated Parameter Use",     noret); break;// (40)
    case 41:GAMMAerror(hdr,"Please Use AQ Rather Than Q_T",noret); break;// (41)
    case 42:GAMMAerror(hdr,"Q_T Use, Parameter Not Type 4",noret); break;// (42)
    case 43:GAMMAerror(hdr,"Q_T Use, Must Be Rank 2",      noret); break;// (43)
    case 44:GAMMAerror(hdr,"Q_T Use, Ignore Isotropic Val",noret); break;// (44)
    case 50:GAMMAerror(hdr,"Invalid Component, m=[-2,2]",  noret); break;// (50)
    default:GAMMAerror(hdr,"Unknown Error",                noret); break;
    }
  }

//    case 15:GAMMAerror(hdr,"Setting Asymmetry to Zero",    noret); break;// (15)
//    case 16:GAMMAerror(hdr,"Setting Default I=1/2 Value",  noret); break;// (16)
//    case 19:GAMMAerror(hdr,"Can't Write To Parameter File",noret); break;// (19)
//    case 22:GAMMAerror(hdr,"Insufficient File Parameters", noret); break;// (22)
//    case 23:GAMMAerror(hdr,"Insufficient PSet Parameters", noret); break;// (23)
//    case 24:GAMMAerror(hdr,"Insufficient Spatial Params",  noret); break;// (24)
//    case 52:GAMMAerror(hdr,"Cannot Write To FileStream",   noret); break;// (52)
//    case 53:GAMMAerror(hdr,"Cannot Output Parameters",     noret); break;// (53)
//19      cout << "Error Setting Quadrupolar Frequency, No I Value Specified";
//40      cout << "You've Used A Generic Spatial Tensor Parameter, Please Update";



volatile void IntQuad::Qfatal(int eidx) const
  {
  Qerror(eidx);				// Output the error message
  if(eidx) Qerror(0, 1);		// State that its fatal
  GAMMAfatal();				// Clean exit from program
  }

/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Function pname
      5     Please Use Class Member Function pname
   default  Unknown Error - pname                                           */

void IntQuad::Qerror(int eidx, const string& pname, int noret) const
  {
  string hdr("Quadrupolar Interaction");
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

// ____________________________________________________________________________
// ii                QUADRUPOLAR INTERACTION SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions set up specific aspects of a quadrupolar interaction.
   As the various interaction parameters interact, the functions MUST be
   private because their misuse could produce an inconsistent interaction.

   The goal here is quite simple. We must determine the following set of values
   for each quadrupolar interaction: { Iqn, QCC, eta, alpha, beta, gamma }
   Complexity arises because we allow that a variety of parameters may be used
   to define these values. Additionally we allow that some parameters may be
   zero (e.g. eta) and that defaults may automatically be used. But, the end
   result remains the same, we seek the set of values that define a 
   quadrupolar interaction in GAMMA.                                         */

// ----------------------------------------------------------------------------
//                      Complete Quadrupolar Interaction
// ----------------------------------------------------------------------------

/* These functions will try and get all the parameters required to define a
   quadrupolar interaction: { Iqn, QCC, eta, alpha, beta, gamma }.           */

bool IntQuad::getQI(const ParameterSet& pset,
         double& Iqn, double& qcc, double& eta, EAngles& EA,
                                                      int idx, bool warn) const
  {
//                      Get The Spin Quantum Number
//                    (And Perhaps The Isotope Type)

  string II;                                    // String for isotope name
  Isotope ISI;                                  // Isotope for spin
  bool TFI = false;                             // Flag if we know isotope
  if(getIso(pset,II,idx,0))                     // 1. Try for isotope name
    {                                           //    If successful, check it
    if(!SpinCheck(II)) return false;            //    Insure valid isotope
    ISI = Isotope(II);                          //    An isotope of type II
    if(!SpinCheck(ISI,TFI,0)) return false;     //    Disallow electron spin
    TFI = true;                                 //    We know the isotope
    Iqn = ISI.qn();                             //    Know spin I quantum #
    if(!SpinCheck(Iqn,TFI,0)) return false;	//    Insure Iz > 1/2
    }
  else if(!getIqn(pset,"",Iqn,idx,0))           // 2. Try for spin quant. #
  { if(!SpinCheck(Iqn,true,0)) return false;	//    Insure valid  qn
  }
  else { Iqn = 1.0; }                           // 3. Use default qn of 1.0

//  Try To Directly Read Quadrupolar Cartesian Spatial Tensor Components
//         Interaction Via { Iqn, Qxx, Qxy, Qxz, Qyy, Qyz, Qzz }

//  1.) Quv values set the anisotropy & asymmetry: {QCC, ETA}
//  2.) If any Quv off-diagonals (u != v) set, Q array sets orientation also
//  3.) If no Quv off-diagonals, orientation set with specified Euler angles
//  4.) The input Quv values are taken to be in kHz
//  5.) We don't mind that no orientation is set, default is PAS
//  6.) If used, Euler angles by either a set or 3 individual angles

  coord QiQzQe;                                 // For Qiso, Qzz, Qeta
  if(getACart(pset,"Q",QiQzQe,EA,idx,-1,0))     // Try & use Cart. components
    {
    qcc = QiQzQe.y()*1.e3;			//  Set the anisotropy value
    eta = QiQzQe.z();                           //  Set the eta value
    return true;
    }
  else if(getACart(pset,"q",QiQzQe,EA,idx,-1,0))
    {
    qcc = QiQzQe.y()*1.e3;			//  Set the anisotropy value
    eta = QiQzQe.z();                           //  Set the eta value
    return true;
    }

//     Try To Directly Read Quadruplar Coupling, Asymmetry, & Orientation
//               Interaction Via { Iqn, QCC, Qeta, QEAngles }

//  1.) If QCC has been specified, these parameters will be used
//  3.) We don't mind that eta is not set, default will be zero
//  4.) Iqn was set in the 1st section of this function
//  5.) We don't mind that no orientation is set, default is PAS
//  6.) Orientation set with either an Euler angle set or 3 individual angles

  if(getQCC(pset, qcc, idx, 0))			// If QCC has been specified
    {                                           // look for eta & orientation
    string Pbase("Q");                          //   Base for parmeter names
    if(!getAeta(pset, Pbase, eta, idx, 0))	//   Try for asymmetry
      getAeta(pset,"q",eta,idx,0);
    if(!getOrientation(pset,Pbase,EA,idx,0))	//   Try for orientation
      getOrientation(pset,"q",EA,idx,0);
    return true;
    }

  if(warn)                                      // Try as we might, cannot get
    {                                           // interaction from parameters
    Qerror(50, 1);				// Can't set from parameters
    Qerror(51, 1);				// Parameter set insufficient
    }
  return false;
  }

// ----------------------------------------------------------------------------
//             Get Quadrupolar Coupling Value From A Parameter Set
// ----------------------------------------------------------------------------

/* This function attempts to determine the quadrupolar coupling constant from
   parameters in an input parameter set. Currently allowed parameters for QCC
   are:
                 QCC,    QCCkHz,  QCCKHz,  QCCHz,  QCCMHz,
                 WQ,     WQkHz,   WQKHz,   WQHz,   WQMHz,
                 NQCC,   NQCCkHz, NQCCKHz, NQCCHz, NQCCMHz,

           Input                Q       : Quadrupolar interaction (this)
                                pset    : A parameter set
                                idx     : Index value
				qcc     : Quadrupolar coupling
                                warn    : Warning output flag
           Output               TF	: True if quadrupolar coupling found
                                          from parameters in pset
           Note                         : Interaction is NOT altered         */

bool IntQuad::getQCC(const ParameterSet& pset, double& qcc, 
                                                    int idx, bool warn) const
  {
  string Nidx = "";                             // Name addition per index
  if(idx >= 0) Nidx                             // If index exists then set up
    += string("(") + Gdec(idx) + string(")");	// the parameter name to append
  string Qns[15]={ "WQ",      "WQkHz",  "WQKHz",// Quad. interaction parameter
                   "WQHz",    "WQMHz",  "QCC",  // base names to set delzz.
                   "QCCkHz",  "QCCKHz",         // These are order by priority
                   "QCCHz",   "QCCMHz", "NQCC", // in which they will be used
                   "NQCCkHz", "NQCCKHz",        // if encountered in the pset
                   "NQCCHz",  "NQCCMHz" };
  double I = Izval();				// Spin quantum number
  double Ifact = 2.*I*(2.*I-1.)/3.0;            // Scaling if frequency input
  string pname, pstate;                         // Strings for parameter
  ParameterSet::const_iterator item;         // A pix into parameter list
  SinglePar par;
  for(int i=0; i<15; i++)                       // Loop over possible delzz
    {                                           // parameters names
    pname = Qns[i] +  Nidx;                     //      Parameter name
    item = pset.seek(pname);                    //      Seek parameter in pset
    if(item != pset.end())                      //      If it's been found
      {                                         //      parse the parameter
      (*item).parse(pname,qcc,pstate);		//      info and set delzz
      switch(i)                                 //      based on input type
        {
        case  0: case  1: case 2:               //      Here if WQ in kHz
        qcc = qcc*Ifact*1.e3; break;
        case  3:                                //      Here if WQ in Hz
        qcc = qcc*Ifact;      break;
        case  4:                                //      Here if WQ in MHz
        qcc = qcc*Ifact*1.e6; break;
        case  5: case  6: case 7:
        case 10: case 11: case 12:              //      Here if QCC in kHz
        qcc = qcc*1.e3;       break;
        case  8: case 13:                       //      Here if QCC in Hz
        break;
        case  9: case 14:                       //      Here if QCC in MHz
        qcc = qcc*1.e6;       break;
        }
      return true;
      }
    }
  qcc = 0.0;
  return false;
  }

// ----------------------------------------------------------------------------
//                      Complete Quadrupolar Interaction
// ----------------------------------------------------------------------------

/* This function employs all of the "get*" functions in the above sections
   to parse a parameter set for all of the values needed to define a
   quadrupolar interaction, namely { Iqn,QCC,eta,alpha,beta,gamma }. If the
   interaction definition is found, we set the interaction or return false.  */

bool IntQuad::setQI(const ParameterSet& pset, int idx, bool warn)
  {
  double  Iqn;                                  // Our spin quantum number
  double  qcc;                                  // Our quad. coupling    (Hz)
  double  eta;                                  // Our asymmetry        [0,1]
  EAngles EA;                                   // Our Euler angles (radians)
  if(getQI(pset,Iqn,qcc,eta,EA,idx,warn))	// Try and get parameter vals
    {                                           // If successful, set interact
    _QCC = qcc; 	                            //   Set quad. coupling  (Hz)
// sosik
// not sure about setting xi without quantum number...
    double X = xi();                            //   Set interaction constant
    IntRank2::operator=(IntRank2(Iqn,X,eta,EA));//   Set rank 2 interaction
    return true;                                //   Return we are successful
    }
  return false;
  }


// ____________________________________________________________________________
// iii            QUADRUPOLAR INTERACTION CHECKING FUNCTIONS
// ____________________________________________________________________________

 
bool IntQuad::checkIHS(int eidx, int warn)
  {
  if(Ival < 3) 				// Insure quadrupolar interactions
    {					// have I value > 1/2 (HS > 2)
    if(warn)
      {
      Qerror(eidx, 1);			//   If not issue warning eidx
      if(warn > 1) Qfatal(11);	//   Stope execution if needed
      return false;
      }
    }
  return true;
  }


bool IntQuad::checkI(int eidx, int warn)
  {
//sosik
//  if(I<=0.5) 				// Insure quadrupolar interactions
if(Izval() <= 0.5)
    {					// have I value > 1/2 (HS > 2)
    if(warn)
      {
      Qerror(eidx, 1);
      if(warn > 1) Qfatal(11);
      return false;
      }
    }
  return true;
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

/* Quadrupolar interaction construction. For full construction users need to
   define the 1.) irreducible rank 2 spatial tensor, 2.) irreducible rank 2
   spin tensors, and 3.) interaction strength.  Since GAMMA uses normalized
   spatial and spin tensors, the first two are relatively simple. For the
   spatial tensor we need { eta, theta, phi } or { qxx, qyy, qzz, theta, phi }
   and for the spin tensor we need only the spin quantum value involved,
   often 1.0.  The scaling is more variable since it is defined in different
   ways, all of which relate to the anisotropic value of the irreducible part,
   delzz. Variations allow fo the strength to be set using a quadrupolar
   coupling constant or a quadrupolar frequency.  As these are tied into the
   quantum number of the spin involved we have to do some bookkeeping.       */

// ____________________________________________________________________________
// A          QUADRUPOLAR INTERACTION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

IntQuad::IntQuad()                  : IntRank2()   { _QCC = 0.0;     }
IntQuad::IntQuad(const IntQuad &Q1) : IntRank2(Q1) { _QCC = Q1._QCC; }

// ----------------------------------------------------------------------------
//              Direct Constructors Using Spherical Components
// ----------------------------------------------------------------------------

/* Here we need to know the quantum number of the spins involved, the
   quadrupolar coupling, and the interaction asymmetry. the asymmetry
   is unitless and in the range [0,1]

           Input                Q	: Quadrupolar interaction (this)
                                II      : Spin I isotope type
                                Iqn     : Spin I quantum number (e.g. 1.5)
                                QCC     : Quadrupolar coupling value (Hz)
                                eta     : Tensor asymmetry value (default 0)
           Output               none    : Quadrupolar interaction constructed
           Note                         : A fatal error will result if
                                          the spin quantum number isn't > 1/2
           Note                         : Here QCC=delzz                     */


/* Here we need to know the delzz value or else the PAS cartesian spatial 
   Quadrupolar tensor values {qxx, qyy, qzz }. If the interaction is constructed
   using Cartesian components we insist |qzz| >= |qyy| >= |qxx| or else the 
   asymmetry (eta) will not match that set in spherical terms. The other 
   parameters are optional and happy to be left at zero.  These are the
   orientation angles (theta & phi, degrees) & the asymmetry (eta, unitless)
   Unless explicitly specified as arguments, it is assumed that the spin on
   the nucleus is 1.

           Input                Q	: Quadrupolar interaction (this)
	  			qn	: Quantum number
	  			delzz   : Tensor delzz value (Hz)
	  			eta     : Tensor asymmetry value
	  			theta	: Tensor orientation from PAS (deg)
	  			phi	: Tensor orientation from PAS (deg)
	   Output		none    : Quadrupolar interaction constructed */

IntQuad::IntQuad(const string& IsoI, double qcc, double eta, const EAngles& EA)
  {
  if(!SpinCheck(IsoI)) Qfatal(2); 		// Insure we know this isotope
  Isotope II(IsoI);				// Get isotope for this type
  if(!SpinCheck(II, true)) Qfatal(2);		// Insure we are not electron
  _QCC      = qcc*1.e3;				// Set anisotropic value (PPM)
  double Iz = II.qn();                          // Get Iz value of spin
  if(!SpinCheck(Iz, true, true)) Qfatal(2);	// Insure Iz > 1/2
  double X  = xi();                             // Get interaction constant
  bool   F  = false;				// Flag we are spin-self
  IntRank2::operator=(IntRank2(Iz,X,eta,EA,F));	// Use generic interaction
  }

IntQuad::IntQuad(const Isotope& II, double qcc, double eta, const EAngles& EA)
  {
  if(!SpinCheck(II, true)) Qfatal(2);		// Insure we are not electron
  _QCC      = qcc*1.e3;				// Set anisotropic value (PPM)
  double Iz = II.qn();                          // Get Iz value of spin
  if(!SpinCheck(Iz, true, true)) Qfatal(2);	// Insure Iz > 1/2
  double X  = xi();                             // Get interaction constant
  bool   F  = false;				// Flag we are spin-self
  IntRank2::operator=(IntRank2(Iz,X,eta,EA,F));	// Use generic interaction
  }

IntQuad::IntQuad(double Iz, double qcc, double eta, const EAngles& EA)
  {
  if(!SpinCheck(Iz, true, true)) Qfatal(2);	// Insure Iz > 1/2
  _QCC      = qcc*1.e3;				// Set anisotropic value (PPM)
  double X  = xi();                             // Get interaction constant
  bool   F  = false;				// Flag we are spin-self
  IntRank2::operator=(IntRank2(Iz,X,eta,EA,F));	// Use generic interaction
  }

// ----------------------------------------------------------------------------
//              Direct Constructors Using Cartesian Components
// ----------------------------------------------------------------------------

/* Here we try to use the set { Iqn, Qxx, Qyy, Qzz, QEAngles }. The spin
   quantum number sets the dimension of the interaction spin Hilbert space.
   The values of Qxx, Qyy, and Qzz are in PPM and used to determine the
   quadrupolar coupling QCC (kHz) and asymmetry ETA ([0,1]). The interaction
   orientation is taken from Euler angles QEAngles. We will allow for default
   settings of the asymmetry and orientation so that those arguments need 
   not be input during construction.  Their default values are 0.     
   We will insist |Qzz| >= |Qyy| >= |Qxx| or else the asymmetry (eta) will not
   match that set using spherical terms. 

           Input                Q	: Quadrupolar interaction (this)
                                II      : Spin I isotope type
                                Iqn     : Spin I quantum number (e.g. 1.5)
                                QCC     : Quadrupolar coupling value (Hz)
                                eta     : Tensor asymmetry value (default 0)
           Output               none    : Quadrupolar interaction constructed
           Note                         : A fatal error will result if
                                          the spin quantum number isn't > 1/2
           Note                         : Here QCC=delzz                     */

IntQuad::IntQuad(const string& IsoI, const coord& QxQyQz, const EAngles& EA)
  {
  if(!SpinCheck(IsoI)) Qfatal(2); 		// Insure we know this isotope
  Isotope II(IsoI);				// Get isotope for this type
  if(!SpinCheck(II, true)) Qfatal(2);		// Insure we are not electron
  double Iz = II.qn();                          // Get Iz value of spin
  if(!SpinCheck(Iz, true, true)) Qfatal(2);	// Insure Iz > 1/2
  coord QiQzQe = AisoDelzEta(QxQyQz);           // Qwitch to spherical values
  _QCC      = QiQzQe.y()*1.e3;			// Set quadrupolar coupling (Hz)
  double E  = QiQzQe.z();                       // Get asymmetry value
  double X  = xi();                             // Get interaction constant
  bool   F  = false;				// Flag we are spin-self
  IntRank2::operator=(IntRank2(Iz,X,E,EA,F));	// Use generic interaction
  }

IntQuad::IntQuad(const Isotope& II, const coord& QxQyQz, const EAngles& EA)
  {
  if(!SpinCheck(II, true)) Qfatal(2);		// Insure we are not electron
  double Iz = II.qn();                          // Get Iz value of spin
  if(!SpinCheck(Iz, true, true)) Qfatal(2);	// Insure Iz > 1/2
  coord QiQzQe = AisoDelzEta(QxQyQz);           // Qwitch to spherical values
  _QCC      = QiQzQe.y()*1.e3;			// Set quadrupolar coupling (Hz)
  double E  = QiQzQe.z();                       // Get asymmetry value
  double X  = xi();                             // Get interaction constant
  bool   F  = false;				// Flag we are spin-self
  IntRank2::operator=(IntRank2(Iz,X,E,EA,F));	// Use generic interaction
  }

IntQuad::IntQuad(double Iz, const coord& QxQyQz, const EAngles& EA)
  {
  if(!SpinCheck(Iz, true, true)) Qfatal(2);	// Insure Iz > 1/2
  coord QiQzQe = AisoDelzEta(QxQyQz);           // Qwitch to spherical values
  _QCC      = QiQzQe.y()*1.e3;			// Set quadrupolar coupling (Hz)
  double E  = QiQzQe.z();                       // Get asymmetry value
  double X  = xi();                             // Get interaction constant
  bool   F  = false;				// Flag we are spin-self
  IntRank2::operator=(IntRank2(Iz,X,E,EA,F));	// Use generic interaction
  }

// ----------------------------------------------------------------------------
//                    Construction Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions will construct the interaction from parameters in a
   specified GAMMA parameter set.  This is most useful when working with 
   multispin systems that set up many quadrupolar interactions over the
   full system using a single parameter file (or parameters read from a
   single external ASCII file.)

        Input                QI      : Quadrupolar interaction
                             pset    : Parameter set
                             idx     : Interaction index (default -1->none)
                             idxI    : Index for the spin
                             warn    : Flag to warn if no interaction found
        Output               none    : QI interaction constructed
                                       for spin with quantum number qn
                                       and parameters in pset
                                or   : QI interaction constructed
                                       for spin with index idxI
                                       and parameters in pset                */

IntQuad::IntQuad(const ParameterSet& pset, int idx, int warn)
       : IntRank2()
  { 
  if(!setQI(pset, idx, warn?true:false))
    {
    cout << "\n\tTrouble setting interaction";
    exit(-1);
    }
  }

//IntQuad::IntQuad(const std::string& II, const ParameterSet& pset, int idx, int warn)
//   { (*this) = IntQuad(Isotope(II), pset, idx, warn); }

/*
IntQuad::IntQuad(const Isotope& II, const ParameterSet& pset, int idx, int warn)
       : IntRank2()
  { 
  Ival = int(2*II.qn())+1;		// Set I HS value
// sosik
//  I = II.qn();
  bool TF = checkIHS(11, warn?1:0);	// Insure I Hilbert space is OK
  if(!TF)				// If I < 1/2 cannot construct
    {					// any Quad interaction
    if(warn)				//   If warnings to be issued
      {
      if(warn<2) Qerror(2, 1);		//     State I < 1/2 NFG
      else       Qfatal(2); 		// in this parameter set
      }
    return;
    }
  TF = II.electron();			// See if spinis an electron
  if(TF)				// which cannot be involved in
    {					// any Quad interaction
    if(warn)				//   If warnings to be issued
      {
      Qerror(25, 1);			//     State spin is e-
      if(warn>1) Qfatal(2); 		// in this parameter set
      }
    return;
    }
// sosik
//  setQuad();				// Set spin tensor components
//  TF = setSphComp(pset, idx);  		// Try & read spherical components
//  if(!TF)                               // If we can't do spherical ones
//    TF = setCartComp(pset, idx);        // try & read Cartesian components
  }
*/
// sosixx

/*
IntQuad::IntQuad(double qn, const ParameterSet& pset, int idx)
  {
  int twoI = int(2.0*qn); 		// We must insist that I is an integer
  I = double(twoI)/2.0; 		// multiple of 1/2 and that I>1/2
  checkI(2, 11);			// Insure I is OK
  Tsph = new matrix[5];			// For spherical spin componenets
  setTs();				// Set 5 spherical spin components

  SinglePar par;			// A single parameter
  string pname;				// For parameter name
  string pstate;			// For parameter statement
  string sval;				// For parameter value
  string Nidx = "";			// Name addition per index
  if(idx >= 0) Nidx 			// If index exists then set up
    += string("(")+Gdec(idx)+string(")");// the parameter name to append

//---------------------------------------------------------------------------- 
//                Try & Read AQ Quadrupolar Coupling Spatial Tensor	     -
//  			   { delzz, eta, phi, theta }                        -
//									     -
// AQ	(2) : ( QCC, eta, theta, phi )		- Quadrupolar Interaction    -
//									     -
// --------------------------------------------------------------------------- 

  int found = setAsph(pset, idx);	// Look for IntRank2A spatial tensor

// --------------- Try & Read Q: PAS Spatial Quadrupolar Tensor --------------
 
//  int Qflag = 0;			// Flag if interaction found
//  pname = string("AQ") + Nidx;	// Quad. interaction parameter name
//  par = SinglePar(pname);		// Parameter for this name
//  if(pset.seek(par))			//      If parameter is in pset
//    {					//      get it as a spatial tensor
//    item = pset.seek(par);		//      Here's where it is
//    space_T X(pset(item));		//      Set as quad. spatial tensor
//    Qflag++;				//      Flag that we know it
//    }					// Note QCC here read in kHz!
// if(found) return;

  found = setAsphGen(pset, idx);	// Look for generic space tensor
    {
    Asph = new complex[5];		// For spherical space componenets
    setAs();				// Set 5 spherical space components
    return;				// Now we are done
    }
  pname = string("AQ") + Nidx;		// Hard as we look, there seems to
  Qerror(101, pname);			// be no Quad interaction defined
  Qfatal(0); 			// in this parameter set
  }
*/
//IntQuad::IntQuad(int idxI,         const ParameterSet& pset,             int warn=2);
//IntQuad::IntQuad(const string& II, const ParameterSet& pset, int idx=-1, int warn=2);
//IntQuad::IntQuad(double qn,        const ParameterSet& pset, int idx=-1, int warn=2);

 
// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

void IntQuad::operator= (const IntQuad &Q1)  { IntRank2::operator= (Q1); }
     IntQuad::~IntQuad() { }

// ____________________________________________________________________________
// B                     SPATIAL TENSOR COMPONENT ACCESS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//            Spatial Tensor Anisotropy Value (Quadrupolar Coupling)
// ----------------------------------------------------------------------------

/* Since quadrupolar interactions have no isotropic component, the interaction
   strength is entirely defined by the quadrupolar coupling (or delzz value).
   That would be handled entirely by class IntRank2 were it not for the
   different definitions in the literature and related parameters in common
   use.  So, we add in some functions here to make life easy for everybody.
   For this interaction the delzz value is the quadrupolar coupling, QCC,
   and it relates to the quadruplar frequency w  as shown below.
                                               Q
                                                                  2
                              2                      3*NQCC     3e qQ
        DELZZ = NQCC = QCC = e qQ             w   = -------- = --------
                                               Q    2I(2I-1)   2I(2I-1)

   Accordingly, the quadrupolar frequency will be the splitting between
   the 2*I transitions contained in a Quadrupolar Hamiltonian.  This will be
   observed at zero field or when the quadrupolar coupling is weak relative to
   the Zeeman interaction (high field approximation, first order terms only)
   if the tensor is oriented in it's principal axes (PAS).  If the tensor isn't
   aligned in it's PAS, the splitting will vary with orientation according to

                       1           2                2
      W (theta,phi)  = - W  [ 3*cos (theta) - eta*sin (theta)cos(2*phi) ]
       Q               2  Q

  where theta is the angle down from the PAS z axis and phi the angle over
  from the PAS x axis.  Alternatively, if the Zeeman terms aint much stronger
  the splittings won't be equally spaced at all (second order terms).

  Also, keep in mind that the Euler angles {phi,theta,gamma} which are kept
  with the tensor are used to relate the tensor PAS to some (common) set of
  coordinate axes.  They are not the phi and theta used in the above
  formula (unless you which the tensor aligned in the common axis system)

          Input                Q       : Quadrupolar interaction
                               delzz   : Nuclear quad. tensor delzz (Hz)
                               QCC     : Nuclear quad. coupling (Hz)
                               wQ      : Quadrupolar frequency (Hz)
          Output               delzz   : Nuclear quad. tensor delzz
                                         value in (Hz)
          Note                         : This is equivalent to the
                                         quadrupolar coupling

   Again, for this interaction the delzz value is equivalent to the Quadrupolar
   coupling constant. Also, note that since the quadrupolar frequency depends
   upon the spin quantum number it is not defined unless the spin is set.

   Note that the default base class IntRank2A also provides the GAMMA 
   normalized anisotropy and delzz values. These are constants given by

                       1/2
                  [ 5 ]                  ^      3
          del   = |---]  = 0.51503      /_\ A = - del   = 0.77255
             zz   [6PI]                         2    zz                      */

// static double IntRank2A::delzz();                            INHERITED
// static double IntRank2A::delA();                             INHERITED

double IntQuad::QCC()   const { return _QCC; }
double IntQuad::NQCC()  const { return _QCC; }
//double IntQuad::wQ( )   const { return (!I)?0.:3.*_QCC/(2.*I*(2.*I - 1.)); }
// sosk
double IntQuad::wQ( )   const
 { double I=Izval(); return (!I)?0.:3.*_QCC/(2.*I*(2.*I - 1.)); }
 
void   IntQuad::QCC(double   dz) { _QCC = dz; }
// sosik need to reset xi also
void   IntQuad::NQCC(double  dz) { _QCC = dz; }
// sosik need to reset xi also
void   IntQuad::wQ(double W)
  {
// sosik
  double I = Izval();
  if(!I) { Qerror(19, 1); Qfatal(18); }	// Fatal error
  double Ifact = 2.0*I*(2.0*I - 1.0);   // I based denomenator
  _QCC = W*Ifact/3.0;
// sosi need to reset xi also
  }

// ----------------------------------------------------------------------------
//                      Spatial Tensor Asymmetry Value
// ----------------------------------------------------------------------------

// Note that eta spans [0, 1] and is defined to be (Axx-Ayy)/Azz where
// |Azz| >= |Ayy| >= |Axx| in the treatment contained herein.  These functions
// are inherited from the base class IntRank2A

//double IntQuad::eta( ) const             INHERITED       Get the asymmetry
//void   IntQuad::eta(double Qeta)         INHERITED       Set the asymmetry

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
  row_vector IntRank2A::CartComps(double theta, double phi=0) const;         */

// ----------------------------------------------------------------------------
//   Un-normalized Cartesian Components Of Quadrupolar Spatial Tensor Q
//                                                                     uv
// ----------------------------------------------------------------------------

/* These allow one to access the Cartesian elements of the full (unscaled) Q
   tensor at a specific orientation without rotating the entire tensor.  In
   these functions theta is the angle down from the lab. frame z-axis and phi
   the angle over from the lab. frame x-axis.  If no angles are specified, then
   the current orientation is used.  Typically, these elements will be on the
   order of kHz. The relate to GAMMA's normlized Cartesian values according to
   (where Kdel is a Kronecker delta function)

                            1/2                         1/2
                    [ 6*PI ]                    [ 6*PI ]
              q   = | ---- |   * del   * A   =  | ---- |   * QCC * A
               uv   [  5   ]        zz    uv    [  5   ]            uv       */

double IntQuad::qxx() const { return _QCC*RT6PIO5*Axx(); }
double IntQuad::qyy() const { return _QCC*RT6PIO5*Ayy(); }
double IntQuad::qzz() const { return _QCC*RT6PIO5*Azz(); }
double IntQuad::qxy() const { return _QCC*RT6PIO5*Axy(); }
double IntQuad::qyx() const { return _QCC*RT6PIO5*Ayx(); }
double IntQuad::qxz() const { return _QCC*RT6PIO5*Axz(); }
double IntQuad::qzx() const { return _QCC*RT6PIO5*Azx(); }
double IntQuad::qyz() const { return _QCC*RT6PIO5*Ayz(); }
double IntQuad::qzy() const { return _QCC*RT6PIO5*Azy(); }

double IntQuad::qxx(double alpha, double beta, double gamma) const
                { return _QCC*RT6PIO5*Axx(alpha,beta,gamma); }
double IntQuad::qyy(double alpha, double beta, double gamma) const
                { return _QCC*RT6PIO5*Ayy(alpha,beta,gamma); }
double IntQuad::qzz(double alpha, double beta, double gamma) const
                { return _QCC*RT6PIO5*Azz(alpha,beta,gamma); }
double IntQuad::qyx(double alpha, double beta, double gamma) const
                { return _QCC*RT6PIO5*Ayx(alpha,beta,gamma); }
double IntQuad::qxy(double alpha, double beta, double gamma) const
                { return _QCC*RT6PIO5*Axy(alpha,beta,gamma); }
double IntQuad::qzx(double alpha, double beta, double gamma) const
                { return _QCC*RT6PIO5*Azx(alpha,beta,gamma); }
double IntQuad::qzy(double alpha, double beta, double gamma) const
                { return _QCC*RT6PIO5*Azy(alpha,beta,gamma); }
double IntQuad::qxz(double alpha, double beta, double gamma) const
                { return _QCC*RT6PIO5*Axz(alpha,beta,gamma); }
double IntQuad::qyz(double alpha, double beta, double gamma) const
                { return _QCC*RT6PIO5*Ayz(alpha,beta,gamma); }

double IntQuad::qxx(const EAngles& EA) const { return _QCC*RT6PIO5*Axx(EA); }
double IntQuad::qyy(const EAngles& EA) const { return _QCC*RT6PIO5*Ayy(EA); }
double IntQuad::qzz(const EAngles& EA) const { return _QCC*RT6PIO5*Azz(EA); }
double IntQuad::qxy(const EAngles& EA) const { return _QCC*RT6PIO5*Axy(EA); }
double IntQuad::qyx(const EAngles& EA) const { return _QCC*RT6PIO5*Ayx(EA); }
double IntQuad::qxz(const EAngles& EA) const { return _QCC*RT6PIO5*Axz(EA); }
double IntQuad::qzx(const EAngles& EA) const { return _QCC*RT6PIO5*Azx(EA); }
double IntQuad::qyz(const EAngles& EA) const { return _QCC*RT6PIO5*Ayz(EA); }
double IntQuad::qzy(const EAngles& EA) const { return _QCC*RT6PIO5*Azy(EA); }

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
   effectively reorient the spatial tensor.                                  */

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


// ____________________________________________________________________________
// C                       SPIN TENSOR COMPONENT ACCESS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//                         Spin Quantum Value Access
//-----------------------------------------------------------------------------

/*
  double IntRank2T::Izval()          const;                     INHERITED
  int    IntRank2T::IV()             const;                     INHERITED
  int    IntRank2T::HS()             const;                     INHERITED
  double IntRank2T::qn(bool i=false) const;                     INHERITED    */

//-----------------------------------------------------------------------------
//                      Spin Tensor Component Access
//-----------------------------------------------------------------------------

/* There are five spin tensor components for the quadrupolar interaction.

                                              +
                                           m  |
                                T    = (-1)  T
                                 2,m          2,-m

   T  = Tsph[0]; T  = Tsph[1]; T   = Tsph[2]; T  = Tsph[3]; T   = Tsph[4]
    2,0           2,1           2,-1           2,2           2,-2

            1/2
         [1]      2                         1                             1  2
   T   = |-| * [3I - I(I+1)]        T   = - -(I I + I I )          T    = - I
    2,0  [6]      z                  2,1    2  + z   z +            2,2   2  +

   These are provided by the base class IntRank2T.                           */

// ----------------------------------------------------------------------------
//        Single Spin Or Spin Pair Hilbert Space Spin Tensor Components
// ----------------------------------------------------------------------------

/* Here HSs are single spin Hilbert space dimensions and the spin involved
   is indexed by i.                                                      
  
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
//    Additional Function To Handle "Perturbing" Quadrupolar Interactions
// ----------------------------------------------------------------------------

/* These are the blends of the of T2(+/-)1 and T2(+/-)2 components that result
   from 1st order perturbation theory. The Hamiltonian built using these is an 
   adjustment term applicable when treating a quadrupolar interaction as a
   small perturbation to the Zeeman interaction. See the function H1. 

                        2    2                                2     2   
     T T   = I * [ 4 * I - 8I  - 1 ]       T T   = I * [ 2 * I  - 2I  - 1 ]
      1 -1    z         z                   2 -2    z         z              */


matrix IntQuad::T21m1() const
  {
double I = Izval();
int Ival = int(2.*Izval() + 1);
// sosik
//  int Ival  = int(2.*I + 1);			// For 1 spin SOp functions
  matrix IZ = Iz(Ival);				// The operator Iz
  matrix IE = Ie(Ival);				// The operator E
  matrix T  = IZ*(4.*I*(I+1)*IE - 8.*IZ*IZ - IE);
  return T;
  }

matrix IntQuad::T22m2() const
  {
  double I = Izval();
  int Ival = int(2.*Izval() + 1);
  matrix IZ = Iz(Ival);				// The operator Iz
  matrix IE = Ie(Ival);				// The operator E
  matrix T  = IZ*(2.*I*(I+1)*IE - 2.*IZ*IZ - IE);
  return T;
  }

matrix IntQuad::T21m1(const vector<int>& HSs,int i) const
  { return blow_up(T21m1(), HSs, i); }

matrix IntQuad::T22m2(const vector<int>& HSs,int i) const
  { return blow_up(T22m2(), HSs, i); }

// ____________________________________________________________________________
// D           QUADRUPOLAR INTERACTION CONSTANT ACCESS FUNCTIONS
// ____________________________________________________________________________

/* The quadrupolar interaction constant is what sets the overall interaction
   strength and allows GAMMA to track different interaction types using
   "normalized" spatial and spin tensor (classes IntRank2A and IntRank2T). For
   quadrupolar interactions the interaction constant is defined as

                       1/2                  1/2   del               1/2
           Q     [6*pi]     NQCC      [6*pi]         zz         [pi]
         xi    = |----| * --------- = |----|  * --------- = W * |--|
                 [ 5  ]   2I*(2I-1)   [ 5  ]    2I*(2I-1)    Q  [15]

    These values are (always) stored in radians/sec and will combine with the
    GAMMA spatial and spin tensors to form Hamiltonians for this and other
    interactions.
                            ---
                            \     Q      m    Q               Q
            H (theta,phi) = /   Xi * (-1)  * A (theta,phi) * T
             Q              ---               2,m             2,-m
                             m                                               */

// sosik
//double IntQuad::xi( ) const { return _QCC*xRT6PIO5/(2.*I*(2.*I-1.)); }
double IntQuad::xi( ) const
  {
  double I = Izval();
  return _QCC*xRT6PIO5/(2.*I*(2.*I-1.));
  }

// ____________________________________________________________________________
// E               QUADRUPOLAR INTERACTION AUXILIARY FUNCTIONS
// ____________________________________________________________________________

/* The quadrupolar frequency, as defined in GAMMA, relates the the nuclear
   quadrupolar coupling (NQCC) as
                                                 2
                                  3*NQCC       3e qQ
[1]                        W   = --------  = --------
                            Q    2I(2I-1)    2I(2I-1) 

   According to formula [1], the quadrupolar frequency will be the splitting 
   between the 2*I transitions contained in a Quadrupolar Hamiltonian. This 
   will be observed either at zero field or when the quadrupolar coupling is
   weak relative to the Zeeman interaction (high field approximation, 1st order
   Hamiltonian terms only) if the quadrupolar tensor orientation (PAS) is
   aligned with the external field axis. If the tensor isn't aligned this way,
   the splitting will vary with orientation according to

                       1           2                2                        
[2]   W (theta,phi)  = - W  [ 3*cos (theta) - eta*sin (theta)cos(2*phi) ]
       Q               2  Q

   where theta is the angle down from the PAS z axis and phi the angle over
   from the PAS x axis.  Alternatively, if the quadrupolar Hamiltonian is
   competetive with the Zeeman Hamiltonian, i.e. the quadrupolar coupling is
   on the same scale as the spin Larmor frequency, then the ovserved splittings
   won't be equally spaced at all (second order terms).
 
   Also, keep in mind that the Euler angles {phi,theta,gamma} which are kept
   with the tensor are used to relate the tensor PAS to some (common) set of
   coordinate axes.  They are not the phi and theta used in the above
   formula (unless you which the tensor aligned in the common axis system)   */

// ----------------------------------------------------------------------------
//                 Quadrupolar Frequency Conversion Functions
// ----------------------------------------------------------------------------

/* These are just some handy friend functions that allow quick conversion
   between the quadurupolar frequency and quadrupolar coupling.              */

double IntQuad::wQ2QCC(double QwQ,  double I) { return QwQ*2.*I*(2.*I-1.)/3.0; }
double IntQuad::QCC2wQ(double QQCC, double I) { return QQCC*3/(2.*I*(2.*I-1.)); }

// ----------------------------------------------------------------------------
//                       1st Order Quadrupolar Frequency
// ----------------------------------------------------------------------------
/* This is the splitting between transtions when the quadupolar interaction is
   much weaker than the Zeeman interaction.  These will only apply when that
   condition is applicable (but can be corrected by second order shifts). We
   employ equation [2] above.  Functions in this section allow users to obtain
   frequency values individually or produce multiple values in order to
   facilitate simple powder averaging.

           Input                Q       : Quadrupolar interaction
                                theta   : Orientation angle (degrees)
                                phi     : Orientation angle (degrees)
                                Ntheta  : Number of orientation increments
                                Nphi    : Number of orientation increments
           Return               wQ      : Quadrupolar frequency (Hz)
                                          defined for the tensor in
                                          it's current orientation
           Return               mx      : Array of components for building
                                          the 1st order quad. frequency
                                          over the sphere.
           Note                         : Theta spans [0, 180]
           Note                         : Phi spans [0, 360)
           Note                         : If no angles input the result is
                                          for the current tensor orientation */

double IntQuad::wQoriented() const { return wQ0(); } 
double IntQuad::wQ0() const
  {
// sosik
//  if(!I || !_QCC) return 0;			// No value if no quad
  if(!Izval() || !_QCC) return 0;		// No value if no quad
double THETA = _EAs.beta();
double PHI   = _EAs.alpha();
  if(!THETA) return wQ();			// Thats it along z of PAS
  double Ctheta = cos(THETA);			// Need cos(theta)
  if(!ETA) return wQ()*(3.*Ctheta*Ctheta-1.);
  double Stheta = sin(THETA);			// Need sin(theta)
  double C2phi = cos(2.*PHI); 
  return wQ()*0.5*(3.*Ctheta*Ctheta - 1. + ETA*Stheta*Stheta*C2phi);
  }


double IntQuad::wQoriented(double t, double p) const { return wQ0(t,p);}
double IntQuad::wQ0(double theta, double phi) const
  {
if(!Izval()) return 0;		// No value if no quad
// sosik
//  if(!I) return 0;                              // No value if no quad
  double thetarad = theta*DEG2RAD;		// Theta in radians
  double twophirad = 2.0*phi*DEG2RAD;		// Phi in radians
  double Ctheta = cos(thetarad);		// cos(theta)
  double Stheta = sin(thetarad);		// sin(theta)
  double C2phi = cos(twophirad);		// cos(2*phi)
  return 0.5*wQ()*(3.*Ctheta*Ctheta - 1. + ETA*Stheta*Stheta*C2phi);
  }

matrix IntQuad::wQoriented(int Nt, int Np) const { return wQ0(Nt,Np); }
matrix IntQuad::wQ0(int Ntheta, int Nphi) const
  {
  int nc = Ntheta;				// Number of columns to track
  if(Nphi > nc) nc = Nphi;			// (use max of Ntheta or Nphi)
  matrix mx(4, nc, complex0);			// Array for values
  double theinc = PI/double(Ntheta-1);		// Angle increment (radians)
  double theta = 0;				// This will be the angle
  double Ctheta, Stheta;
  int j=0;
  for(j=0; j<Ntheta; j++, theta += theinc)	// Loop required increments
    {
    Ctheta = cos(theta);
    Stheta = sin(theta);
    mx.put(1.5*Ctheta*Ctheta-0.5, 0, j);	// 	1/2*(3*cos(t)cos(t)-1)
    mx.put(0.5*ETA*Stheta*Stheta, 1, j);	// 	1/2*eta*sin(t)sin(t)
    mx.put(Stheta, 2, j);			// 	Row of sin(thetas)
    }
  if(!ETA) return mx;				// If symmetric, forget phi
  double phiinc2 = 4.0*PI/double(Nphi);		// 2x Angle increment (radians)
  double phi2 = 0;				// This will be the 2x angle
  for(j=0; j<Nphi; j++, phi2+= phiinc2)		// Loop required increments
    mx.put(cos(phi2),3,j);			//	cos(2*phi)
  return mx;
  }

// ----------------------------------------------------------------------------
//                     2nd Order Central Transition Shifts
// ----------------------------------------------------------------------------

/* These functions apply only to spins have I = m * 1/2 where m is odd & > 1!
   The value(s) returned are from a 2nd order perturbation treatment. 

                      2
                   - w
       (2)               Q   [          3 ]   [               4 
[3]   w (theta,phi) = ---- * | I(I+1) - - | * | A(phi,eta)*cos (theta)]
       Q              6*Om   [          4 ]   [                     z
        -1/2,1/2                      
                                                        2                     ]
                                        - B(phi,eta)*cos (theta) - C(phi,eta) |
                                                                              ]

                -27   9                      3                     2
   A(phi,eta) = --- + - * eta * cos(2*phi) - - [ eta * cos(2*phi) ]
                 8    4                      8

                30   1      2                      3                     2
   B(phi,eta) = -- - - * eta  - 2*eta*cos(2*phi) + - [ eta * cos(2*phi) ]
                 8   2                             4

                -3   1      2   1                  3                     2
   C(phi,eta) = -- + - * eta  - -*eta*cos(2*phi) - - [ eta * cos(2*phi) ]
                 8   3          4                  8                         */
 
        // Input                Q	: Quadrupolar interaction
        // Input                Om	: Field Strength (Larmor in Hz)
        // Return               wQ      : Shift in the central transition
	//				  frequency due to 2nd order effects
	//				  in an I=n*1/2 n=3,5,7,.... spin
	// Note				: The return is zero if I not proper


double IntQuad::wQcentral(double Om) const
  {
// sosik
double I = Izval();
double THETA = _EAs.beta();
double PHI   = _EAs.alpha();
  if(int(2*I)%2) return 0.0;				// Insure I = 3/2, 5/2, ....
  double Wo = 3.0*_QCC/(2.*I*(2.*I-1.));		// Base Quad. frequency
  double Ifact = I*(I+1) - 0.75; 			// Part of the prefactor
  double prefact = -Wo*Wo*Ifact/Om;			// Majority of the prefactor
  double Ctheta = cos(THETA); 				// Here is cosine(theta)
  double Cthetasq = Ctheta*Ctheta;			// Here is cos(theta*cos(theta)
  if(!ETA)						// For symmetric interaction,
    return (prefact/16.)*(1.-Cthetasq)*(9.*Cthetasq-1.);// use a simpler formula
  double thov8 = 3.0/8.0;
  double C2phi = cos(2.0*PHI);
  double etaC2phi = ETA*C2phi;
  double eC2phisq = etaC2phi*etaC2phi;
  double etasq = ETA*ETA;
  double A = -9.0*thov8 + 2.25*etaC2phi - thov8*eC2phisq;
  double B = 10.0*thov8 - 0.5*etasq - 2.0*etaC2phi + 0.75*eC2phisq;
  double C = -thov8 + (1.0/3.0)*etasq - 0.25*etaC2phi - thov8*eC2phisq;
  return (prefact/6.0)*(A*Cthetasq*Cthetasq + B*Cthetasq + C);
  }

matrix IntQuad::wQcentral(int Ntheta, int Nphi)
  {
  int nc = Ntheta;				// Number of columns to track
  if(Nphi > nc) nc = Nphi;			// (use max of Ntheta or Nphi)
  matrix mx(5, nc, complex0);			// Array for values
// sosik
double I=Izval();
  if(!(int(2*I)%2)) return mx;			// Insure I = 3/2, 5/2, ....
  double theinc = PI/double(Ntheta-1);		// Angle increment (radians)
  double theta = 0;				// This will be the angle
  int j=0;
  for(j=0; j<Ntheta; j++, theta += theinc)	// Loop required increments
    {
    mx.put(cos(theta), 0, j);			// 	Row of cos(thetas)
    mx.put(sin(theta), 1, j);			// 	Row of sin(thetas)
    }
  if(!ETA) return mx;				// If symmetric, forget phi
  double phiinc2 = 4.0*PI/double(Nphi);		// 2x Angle increment (radians)
  double phi2 = 0;				// This will be the 2x angle
  double etasq = ETA*ETA;			// Need eta*eta
  double m3o8 = -3.0/8.0;
  double Ao = -27.0/8.0;			// Constant part of A
  double Bo = 30.0/8.0 - 0.5*etasq;		// Constant part of B
  double Co = m3o8 + etasq/3.0;			// Constant part of C
  double C2phi, etaC2phi, eC2phisq;
  for(j=0; j<Nphi; j++, phi2+= phiinc2)		// Loop required increments
    {
    C2phi = cos(phi2);				//	cos(2*phi)
    etaC2phi = ETA*C2phi;			//	eta*cos(2*phi)
    eC2phisq = etaC2phi*etaC2phi;		//	[eta*cos(2phi)]^2
    mx.put(Ao+2.25*etaC2phi+m3o8*eC2phisq,2,j);	//	Row of A components
    mx.put(Bo-2.*etaC2phi+0.75*eC2phisq,3,j);	//	Row of B components
    mx.put(Co-0.25*etaC2phi+m3o8*eC2phisq,4,j);	//	Row of C components
    }
  return mx;
  }


// ---------------------- 2nd Order Quadrupolar Frequency ---------------------

// This is the splitting between transtions when the quadupolar interaction is
// weaker but not much much weaker than the Zeeman interaction.  These will
// only apply when that condition is applicable and should be added to the
// first order frequency (shifts)

double IntQuad::wQ1(double Om, double m) const
 
        // Input                Q	: Quadrupolar interaction
        //			Om	: Field Strength (Larmor in Hz)
	//			m       : Spin angular momentum z quantum #
        // Return               wQ      : 2nd order quadrupolar frequency
	//				  shift (Hz) for the transition m-1,m

//  [2]              -2   [ Xi ]2 [
// w   (theta,phi) = -- * | -- |  | A A  (theta,phi)*[24m(m-1) - 4I(I+1) + 9]
//  Q                Om   [ 2  ]  [  1 -1
//   m-1,m                      
//                                1                                           ]
//                              + -*A A  (theta,phi)*[12m(m-1) - 4I(I+1) + 6] |
//                                2  2 -2                                     ]

  {
  double mmm1 = m*(m-1);
// sosik - use of A2x can be improved
double I = Izval();
  double fourIIp1 = 4.0*I*(I+1);
  double A1Am1 = Re(A21()*A2m1());
  double A2Am2 = Re(A22()*A2m2());
  double W =  A1Am1*(24.0*mmm1 - fourIIp1 + 9.0);
  double Xi = xi();
  W += 0.5*A2Am2*(12.0*mmm1 - fourIIp1 + 6.0);
  W *= (-Xi*Xi/(2.*Om));
  return W;
  }

double IntQuad::wQ1(double Om, double m, double theta, double phi) const

        // Input                Q	: Quadrupolar interaction
        //			Om	: Field Strength (Larmor in Hz)
	//			m       : Spin angular momentum z quantum #
	//			theta   : Orientation angle (degrees) 
	//			phi	: Orientation angle (degrees) 
        // Return               wQ      : 2nd order quadrupolar frequency
	//				  shift (Hz) for the transition m-1,m

  {
  double mmm1 = m*(m-1);
// sosik
double I = Izval();
  double fourIIp1 = 4.0*I*(I+1);
  complex Asph1 = A21(phi, theta, 0.0);
  complex Asph2 = A22(phi, theta, 0.0);
  complex A1Am1 = -Asph1*conj(Asph1);
  complex A2Am2 = Asph2*conj(Asph2);
  double W =  Re(A1Am1)*(24.0*mmm1 - fourIIp1 + 9.0);
  double Xi = xi();
  W += 0.5*Re(A2Am2)*(12.0*mmm1 - fourIIp1 + 6.0);
  W *= (-Xi*Xi/(2.*Om));
  return W;
  }


matrix IntQuad::wQ1(int Ntheta, int Nphi)
 
        // Input                Q	: Quadrupolar interaction 
        //                      Ntheta	: Number of orientation increments
        //                      Nphi	: Number of orientation increments
        // Return               mx 	: Array of components for building
	//				  second order quad. frequency shifts 
	// Note				: Theta spans [0, 180]
	// Note				: Phi spans [0, 360)


//   A   (theta , phi ) = <0|mx|i> + Re(<2|mx|j>)*Re(<4|mx|i> + i*Im(<2|mx|i>)
//    2,1      i     j

//   A   (theta , phi ) = <1|mx|i>
//    2,2      i     j
//                      + Re(<3|mx|j>)*Re(<5|mx|i> + i*Im(<3|mx|j>)*Im(<5|mx|i>

  {
  int nc = Ntheta;				// Number of columns to track
  if(Nphi > nc) nc = Nphi;
  matrix mx(6, nc, complex0);			// Array for values
  double theinc = PI/double(Ntheta-1);		// Angle increment (radians)
  double theta = 0.0;				// This will be the angle
  double fac1 = 3.0*xRT5O24PI;			// A21 prefactor
  double fac2 = 3.0*xRT5O96PI;;			// A22 prefactor
  double Ctheta, Stheta;
  int j=0;
  for(j=0; j<Ntheta; j++)			// Loop required increments
    {
    Ctheta = cos(theta);
    Stheta = sin(theta);
    mx.put(fac1*Ctheta*Stheta, 0, j);		//	Row of A21 parts
    mx.put(fac2*Stheta*Stheta, 1, j);		//	Row of A22 parts
    mx.put(Stheta, 4, j); 
    mx.put(Ctheta, 5, j); 
    theta += theinc; 
    }
  double phiinc2 = 4.0*PI/double(Nphi);		// 2x Angle increment (radians)
  double phi2 = 0;				// This will be the 2x angle
  double facE1 = ETA*xRT5O24PI;			// A21 prefactor
  double facE2 = ETA*xRT5O96PI;			// A22 prefactor
  double C2phi, S2phi;
  for(j=0; j<Nphi; j++)				// Loop required increments
    {
    C2phi = cos(phi2);
    S2phi = sin(phi2);
    mx.put(facE1*complex(-C2phi,S2phi), 2, j);	//	Row of A21 parts
    mx.put(facE2*complex(C2phi,-2.*S2phi),3,j);	//	Row of A22 parts
    phi2 += phiinc2; 
    }
  return mx;
  }







// ____________________________________________________________________________
// F                        PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//      Functions To Make A Parameter Set From A Quadrupolar Interaction
//-----------------------------------------------------------------------------

/* This class has no implicit knowledge of the nuclear spin type. The preferred
   means of specifying a quadrupolar interaction when there is no knowledge of
   the spin isotope types is that which is used for filling up the parameter 
   set in the functions below. The base parameters of individual interactions
   in that case are { QI(#), QCC(#), Qeta(#), QTheta(#), QPhi(#) }.

          Input                QI      : Quadrupolar interaction
                               idx     : Interaction index (default -1)
                               pfx     : Interaction 2nd indx (def -1)
          Output               pset    : Parameter set with only QI
                            OR TF      : Return is FALSE if proper parameters
                                         for the interaction were not found
                                         Quadrupolar interaction parameters

    Function                                 Purpose
  ------------         -------------------------------------------------------
  ParameterSet         Convert interaction into a parameter set
  operator +=          Adds interaction to existing parameter set (friend)
  PSetAdd              Adds interaction to existing parameter set, this
                       allows for an index and a prefix in the parameters   */

IntQuad::operator ParameterSet( ) const
   { ParameterSet pset; pset += *this; return pset; }

void operator+= (ParameterSet& pset, const IntQuad &Q)
   { Q.PSetAdd(pset); }

void IntQuad::PSetAdd(ParameterSet& pset, int idx, int pfx) const
  {
  string suffx;                                // Parameter suffix
  if(idx != -1)                                // Only use suffix if idx
    suffx = string("(")+Gdec(idx)+string(")"); // is NOT -1
  string prefx;                                // Parameter prefix
  if(pfx != -1)                                // Only use prefix if pfx
    prefx = string("[")+Gdec(pfx)+string("]"); // is NOT -1

  string pname  = prefx +string("QI") + suffx;  // Add QI spin quantum number
  string pstate = string("Spin Quantum Number");
//  double pdatad = I;
// sosik
  double pdatad = Izval();
  SinglePar par = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);

  pname  = prefx + string("QCC") + suffx;	// Quadrupolar coupling
  pstate = string("Quadrupolar Coupling (kHz)");
  pdatad = QCC()*1.e-3;
  par    = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);
  
  pname  = prefx + string("Qeta") + suffx;	// Asymmetry
  pstate = string("Quadrupolar Asymmetry");
  pdatad = ETA;
  par    = SinglePar(pname, pdatad, pstate);
  pset.push_back(par);

  pname  = prefx + string("QEAngles") + suffx;  // Add orientation (degrees)
  pstate = string("Quadrupolar Euler Angles (deg)");
  double a = _EAs.alpha() * RAD2DEG;
  double b = _EAs.beta()  * RAD2DEG;
  double g = _EAs.gamma() * RAD2DEG;
  coord EA(a, b, g);
  par = EA.param(pname, pstate);
  pset.push_back(par);
  }

// ----------------------------------------------------------------------------
//  Functions To Output Quadrupolar Interaction To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

        // Input                SA      : Shift anisotropy interaction
        //                      fn      : Output file name
        //                      ofstr   : Output file stream
        //                      idx     : Interaction index (default -1)
        //                      pfx     : Interaction 2nd indx (def -1)
        //                      warn    : Warning level
        // Output               none    : Shift anisotropy interaction is
        //                                written as a parameter set to
        //                                file filename or output stream ofstr

int IntQuad::write(const string &filename, int idx, int pfx, int warn) const

  {
  ofstream ofstr(filename.c_str());     // Open filename for input
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!write(ofstr, idx, pfx, w2))       // If file bad then exit
    {
    Qerror(1, filename, 1);            // Filename problems
    if(warn>1) Qfatal(20);             // Fatal error
    return 0;
    }
  ofstr.close();                        // Close it now
  return 1;
  }

int IntQuad::write(ofstream& ofstr, int idx, int pfx, int warn) const
  {
  ParameterSet pset;                    // Declare a parameter set
  PSetAdd(pset, idx, pfx);              // Add in interaction parameters
  int w2 = 0;                           // Added warning level
  if(warn) w2 = 1;
  if(!pset.write(ofstr, w2))            // Use parameter set to write
    {                                   // out the interaction parameters
    if(warn)
      {
      Qerror(52, 1);                   // Problems writing to filestream
      if (warn>1) Qfatal(53);		// Fatal error
      }
    return 0;
    }
  return 1;
  }


// ____________________________________________________________________________
// G                 QUADRUPOLAR INTERACTION INPUT FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//    Direct Read of Q Interaction From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

/* These next two read functions utilize a single spin/interaction index, for
   the spin involved in the interaction.  They'll try to read the parameter set

                  { [Iso(i)/QI] , QCC, Qtheta, Qphi, Qeta }

   The functions do allow for the suffix (#), but do NOT allow for the prefix
   [#] because multiple quadrupolar interactions can be defined in the same
   file by switching spin/interaction indices.  Users can defined multiple
   sets of interactions in the same file by using Quad. interaction vectors
   (see class IntQuadVec.)

           Input                Q       : Quadrupolar interaction
                                filename: Output file name
                                pset    : Parameter set
                                idx     : Spin index
                                warn    : Warning output label
                                             0 = no warnings
                                             1 = warnings
                                             >1 = fatal warnings
           Output               none    : Quadrupolar factor interaction is
                                          read in from parameters in file
                                          filename or those in parameter set */

bool IntQuad::read(const string& filename, int idx, int warn)
  { 
  ParameterSet pset;                	// Declare a parameter set
  if(!pset.read(filename, warn?1:0))	// Read in pset from file
    {                                   // If we cannot read the file at all
    if(warn)                            // then issue warnings as desired
      {
      Qerror(1, filename, 1);          //      Filename problems
      if(warn > 1) Qfatal(21);      //      Fatal error
      else         Qerror(21);         //      or a warning issued
      }
    return false;                       //  Return that we failed!
    }
  if(!read(pset, idx, warn?1:0))        // Use overloaded function
    {                                   // If we cannot read the file at all
    if(warn)                            // then issue warnings as desired
      {
      Qerror(1, filename, 1);          //      Filename problems
      if(warn > 1) Qfatal(22);      //      Fatal error
      else         Qerror(22);         //      or a warning issued
      }
    return false;                       //  Return that we failed!
    }
  return true;                          //  Return that all OK
  }

bool IntQuad::read(const ParameterSet& pset,    int idx, int warn)
  {
  bool TF = setQI(pset, idx, warn?1:0);
  if(!TF)                                       // If getQI didn't handle
    {                                           // setting interaction proper
    if(warn)                                    // then we'll issue some
      {                                         // warnings if desired
      if(warn > 1) Qfatal(23);              // Fatal error
       else        Qerror(23,1);               // or a warning issued
      }
    return 0;                                   // Return that we failed
    }
  return TF;
  }

// ____________________________________________________________________________
// H                    QUADRUPOLAR HAMILTONIAN FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//               First Order Quadrupolar Interaction Hamiltonian
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
//  These Are For When The Quadrupolar Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------

/*  The secular part of the quadrupolar Hamiltonian is that which returns
    only those components which commute with z axis rotations.  Here we have

     [1]                       Q    Q                        Q      (0)
    H   (alpha,beta,gamma) = Xi  * A   (alpha,beta,gamma) * T    = H
     Q                              2,0                      2,0    Q

    The overall strength of the interaction very much depends upon
    interaction constant which is related to the quadrupolar coupling.

           Input                Q       : Quadrupolar interaction
                                alpha   : Euler angle (radians)
                                beta    : Euler angle (radians)
                                gamma   : Euler angle (radians)
                                EA      : Euler angles (radians)
                                HSs     : Array of spin Hilbert spaces
                                i       : Spin index (in HSs)
           Output               H0      : The secular part of the
                                          quadrupolar Hamiltonian
                                             (default basis, Hz)
           Note                         : Also called the 1st order Q
                                          interaction (perturbation theory)
           Note                         : Rotationally invariant about z
           Note                         : No HSs: return in the spin Hilbert
                                          space of dimension (2I+1)
                                          With HSs: return in composite
                                          spin space                         */

matrix IntQuad::H0() const
  { return IntRank2::H0(); }
matrix IntQuad::H0(double A, double B, double G) const
  { return IntRank2::H0(A,B,G); }
matrix IntQuad::H0(const EAngles& EA) const
  { return IntRank2::H0(EA); }

matrix IntQuad::H0(const vector<int>& HSs, int i) const
  { return IntRank2::H0(HSs,i); }
matrix IntQuad::H0(const vector<int>& HSs, int i, double A, double B, double G) const
  { return IntRank2::H0(HSs,i,A,B,G); }
matrix IntQuad::H0(const vector<int>& HSs, int i, const EAngles& EA) const
  { return IntRank2::H0(HSs,i,EA); }

// ----------------------------------------------------------------------------
//              Second Order Quadrupolar Interaction Hamiltonians
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
//  These Are For When The Quadrupolar Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------

/* When the quadrupolar Hamiltonian is weak relative to the Zeeman Hamiltonian
   but not much much weaker one may apply perturbation theory to determine a
   correction to the zeroth order Hamiltonian. This correction term is

    [2]          -1   [ 1    ]2 [                          2     2
   H   (A,B,G) = -- * | - w  |  | 2*A A  (A,B,G)*Iz*(4*I - 8Iz - 1)
    Q            Om   [ 3  Q ]  [    1 -1

                                                                2     2     ]
                                        + 2*A A  (A,B,G)*Iz*(2*I - 2Iz - 1) |
                                         2 -2                               ]

           Input                Om      : Field Strength (Larmor in Hz)
                                alpha   : Euler angle (radians)
                                beta    : Euler angle (radians)
                                gamma   : Euler angle (radians)
           Output               H1      : The 2nd order secular part of the
                                          quadrupolar Hamiltonian (default
                                          basis, Hz) for the interaction
           Note                         : Also called the 1st order quadrupolar
                                          interaction (perturbation theory)
           Note                         : This will be zero in PAS if no eta
                                          and small if Om >> QCC              */

matrix IntQuad::H1(double Om) const
  {
// sosi can use A2x better
  complex twoA1Am1 = 2.0*A21()*A2m1();
  complex twoA2Am2 = 2.0*A22()*A2m2();
  double Xi  = xi();
  matrix Hmx =  twoA1Am1*T21m1();
  Hmx       += twoA2Am2*T22m2();
  Hmx       *= (-Xi*Xi/(4.*Om));
  return Hmx;
  }
 
matrix IntQuad::H1(double Om, double alpha, double beta, double gamma) const
  {
  complex Asph1 = A21(alpha, beta, gamma);
  complex Asph2 = A22(alpha, beta, gamma);
  complex twoA1Am1 = -2.0*Asph1*conj(Asph1);
  complex twoA2Am2 =  2.0*Asph2*conj(Asph2);
  double Xi  = xi();
  matrix Hmx =  twoA1Am1*T21m1();
  Hmx       += twoA2Am2*T22m2();
  Hmx       *= (-Xi*Xi/(4.0*Om));
  return Hmx;
  }
 
matrix IntQuad::H1(double Om, const EAngles& EA) const
  {
  complex Asph1 = A21(EA);
  complex Asph2 = A22(EA);
  complex twoA1Am1 = -2.0*Asph1*conj(Asph1);
  complex twoA2Am2 =  2.0*Asph2*conj(Asph2);
  double Xi  = xi();
  matrix Hmx =  twoA1Am1*T21m1();
  Hmx       += twoA2Am2*T22m2();
  Hmx       *= (-Xi*Xi/(4.0*Om));
  return Hmx;
  }

matrix IntQuad::H1(vector<int>HSs, int i, double Om) const
  { return blow_up(H1(Om), HSs, i); }
 
matrix IntQuad::H1(vector<int>HSs, int i, double Om,
                                          double A, double B, double G) const
  { return blow_up(H1(Om, A, B, G), HSs, i); }

matrix IntQuad::H1(vector<int>HSs, int i, double Om, const EAngles& EA) const
  { return blow_up(H1(Om, EA), HSs, i); }
 
 
 
// ----------------------------------------------------------------------------
//      Summed First & Second Order Quadrupolar Interaction Hamiltonians
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
//  These Are For When The Quadrupolar Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------
 
/*  This function returns only the secular part of the second order quadrupolar
    Hamiltonian (from perturbation theory).  Note that this still assumes that
    the quadrupolar interaction is a perturbation to to the Zeeman interaction.

                                  (0)    (1)      [1]    [2]
                         H    =  H    + H     =  H    + H
                          Q       Q      Q        Q      Q
  
        // Input                Om	: Field Strength (Larmor in Hz)
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HQ	: The 2nd order secular part of the
        //                                quadrupolar Hamiltonian (default
        //                                basis, Hz) for the interaction
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note				: Also called the 1st order quadrupolar
        //                                interaction (perturbation theory)
        // Note                      	: Rotationally invariant about z
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1 unless the
	//				  composite Hilbert space is given    */

matrix IntQuad::Hw(double Om) const
  { return H0() + H1(Om); }

matrix IntQuad::Hw(double Om, double theta, double phi) const
  { return H0(phi, theta, 0.0) + H1(Om, phi, theta, 0.0); }

matrix IntQuad::Hw(vector<int>HSs, int i, double Om) const
  { return H0(HSs, i) + H1(HSs, i, Om); }

matrix IntQuad::Hw(vector<int>HSs, int i, double Om, double theta, double phi) const
  { return H0(HSs, i, phi, theta, 0.0) + H1(HSs, i, Om, phi, theta, 0.0); }

// ----------------------------------------------------------------------------
//                 Full Quadrupolar Interaction Hamiltonians
//       These Have No Approximations & Are Independent of Field Strength
//     (They Are Correct in The Laboratory Frame, NOT in a Rotating Frame) 
// ----------------------------------------------------------------------------

/* Here we make no approximations. As such, these Hamiltonians are not those
   used in the rotating frame. Here we simply sum over the 5 spatial and spin
   tensor components to produce the Hamiltonian. Recall that the 5 spatial
   tensor components, indexed by the spin angular momentum component
   m = {-2,-1,0,1,2}, map into the Asph array by the following:

                        {  0->0, 1->1, -1->2, 2->2, -2->4 }

                       -2,2
                        ---   Q       m    Q                  Q
      H (theta, phi) =  \   Xi  * (-1)  * A   (theta, phi) * T    
       Q                /                  2,m                2,-m
                        ---
                         m


           Input                Q       : Quadrupolar interaction
           Note                         : This will return in the spin Hilbert
                                          space of dimension (2I+1)          */

matrix IntQuad::H( ) const
  {
  matrix Hmx =  Acomp(0)*T0;
         Hmx -= Acomp(1)*T1 + Acomp(2)*Tm1;
         Hmx += Acomp(3)*T2 + Acomp(4)*Tm2;
  return xi()*Hmx;
  }

matrix IntQuad::H(double theta, double phi) const
  {
  matrix Hmx  = A20( phi,theta,0)*T0;
         Hmx -= A2m1(phi,theta,0)*T1 + A21(phi,theta,0)*Tm1;
         Hmx += A2m2(phi,theta,0)*T2 + A22(phi,theta,0)*Tm2;
  return xi()*Hmx;
  }

matrix IntQuad::H(const vector<int>& HSs, int i) const
  {
  matrix Hmx  = Acomp(0)*T20(HSs,i);
         Hmx -= Acomp(1)*T21(HSs,i) + Acomp(2)*T2m1(HSs,i);
         Hmx += Acomp(3)*T22(HSs,i) + Acomp(4)*T2m2(HSs,i);
  return xi()*Hmx;
  }

matrix IntQuad::H(const vector<int>& HSs, int i, double T, double P) const
  {
  matrix Hmx  = A20( P,T,0)*T20(HSs,i);
         Hmx -= A2m1(P,T,0)*T21(HSs,i) + A21(P,T,0)*T2m1(HSs,i);
         Hmx += A2m2(P,T,0)*T22(HSs,i) + A22(P,T,0)*T2m2(HSs,i);
  return xi()*Hmx;
  }


//matrix IntQuad::H(double theta, double phi) const
 
	// Input		Q	: Quadrupolar interaction
	// Note				: This will return in the spin Hilbert
	//				  space of dimension 2I+1

//                  Q
//                 w  [      2                     2                      2  2
// H (theta,phi) = _  | [3cos (theta) - 1 + eta*sin (theta)*cos(2*phi)][3I -I ]
//  Q              4  [                                                   z  

//                                                                    
//     + sin(theta)[3cos(theta) - eta*cos(theta)cos(2*phi)]*[I I + I I + I ]
//                                                            + z   z +

//                           2                               2           2   2
//                    + [3sin (theta) + eta*cos(2*phi)*(1+cos (theta))][I + I ]
//                                                                       +   -

//                                                                     2   2  ]
//                                     + i*eta*sin(2*phi)*cos(theta)*[I - I ] |
//                                                                     +   -  ]

//  {

// ----------------------------------------------------------------------------
//
//		  	Heres One Way to Do This			
//
//complex Asph20 = A0(theta,phi); 
//complex Asph21 = A1(theta,phi); 
//complex Asph22 = A2(theta,phi); 
//matrix Hmx = Asph20*Tsph[0];
//if(norm(Asph21))
//  Hmx -= Asph21*Tsph[2] - conj(Asph21)*Tsph[1];
//if(norm(Asph22))
//  Hmx += Asph22*Tsph[4] + conj(Asph22)*Tsph[3];
//return Xi*Hmx;
//
// ----------------------------------------------------------------------------
//
//		  	Heres A Second Way to Do This			

//  complex Asph20 = A0(theta,phi); 
//  complex Asph21 = A1(theta,phi); 
//  complex Asph22 = A2(theta,phi); 
//  int Ival = int(2.*I + 1);			// For 1 spin SOp functions
//  matrix IZ = Iz(Ival);				// The operator Iz
//  matrix IP = Ip(Ival);				// The operator I+
//  matrix IM = Im(Ival);				// The operator I-
//  matrix IE = Ie(Ival);				// The operator E
//  matrix IPsq = IP*IP;				// The I+I+ operator
//  matrix IMsq = IM*IM;				// The I-I- operator
//  matrix IX = Ix(Ival);				// The operator Ix
//  matrix IY = Iy(Ival);				// The operator Iy
//  matrix IXIZpIZIX = IX*IZ+IZ*IX;
//  matrix IYIZpIZIY = IY*IZ+IZ*IY;
//  matrix Hmx = Asph20/sqrt(6.0)*(3.*IZ*IZ-(I*(I+1))*IE);
//  Hmx -= Re(Asph21)*IXIZpIZIX + Im(Asph21)*IYIZpIZIY;
//  Hmx += (.5*Re(Asph22))*(IPsq+IMsq) + (complexi*.5*Im(Asph22))*(IMsq-IPsq);
//  return xi()*Hmx;

//
// ----------------------------------------------------------------------------
//
//		  	Heres A Third Way to Do This			
//
//  double Ctheta = cos(theta*DEG2RAD);		// cos(theta)
//  double Stheta = sin(theta*DEG2RAD);		// sin(theta)
//  double C2phi = cos(2.0*phi*DEG2RAD);	// cos(2.0*phi)
//  double S2phi = sin(2.0*phi*DEG2RAD);	// sin(2.0*phi)
//  double Cthetasq = Ctheta*Ctheta;		// cos(theta)*cos(theta)
//  double Sthetasq = Stheta*Stheta;		// sin(theta)*sin(theta)
//  matrix Hmx = (3.*Cthetasq-1.+ETA*Sthetasq*C2phi)*(3.*IZ*IZ-(I*(I+1))*IE);
//  Hmx += (1.5*Sthetasq + 0.5*ETA*C2phi*(1.+Cthetasq))*(IPsq+IMsq);
//  Hmx += (0.5*ETA*S2phi*Ctheta)*(IPsq-IMsq);
//  Hmx *= (wQ()/12.0);
//  return Hmx;
//  }


// ____________________________________________________________________________
// G      QUADRUPOLAR TENSOR SPHERICAL COMPONENTS FOR POWDER AVERAGES
// ____________________________________________________________________________
 
//matrix IntQuad::PTcomp(int comp)

	// Input		Q	: Quadrupolar interaction
	//			comp	: Spin operator component [0, 2]
   	// Output		T       : A spin operator (matrix) for the
	//				  quadrupolar Hamiltonian under
	//				  perturbation theory.
 
//                           1/2
//                          [1]      2                         1                             1  2
// First Order (comp=0) = |-| * [3I - I(I+1)]        T   = - -(I I + I I )          T    = - I
//                        2,0  [6]      z                  2,1    2  + z   z +            2,2   2  +
 



// ____________________________________________________________________________
// N                 QUADRUPOLAR HAMILTONIAN FRIEND FUNCTIONS
// ____________________________________________________________________________

// These Don't Depend on Class IntQuad Variables At All.  The Functions Just
// Mimic The Class' Hamiltonian Generation.

/*

matrix HQ0(double qn, double wQo, double eta, double theta, double phi)
 
	// Input		qn	: Quantum number (1, 1.5, 2.5,...)
	//			wQo     : PAS Quadrupolar frequency
	//			eta     : Quadrupolar asymmetry [0,1]
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


matrix HQ1(double Om, double qn, double wQo, double eta,
                                                      double theta, double phi)
 
        // Input                Om	: Field Strength (Larmor in Hz)
	// 			qn	: Quantum number (1, 1.5, 2.5,...)
	//			wQo     : PAS Quadrupolar frequency
	//			eta     : Quadrupolar asymmetry [0,1]
        //                      theta	: Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HQ1	: The 2nd order secular part of the
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
 
//  in accordance with the article by P.P. Man "Quadrupolar Interactions" in
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
*/ 

// ____________________________________________________________________________
// O                        PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

int IntQuad::set(const ParameterSet& pset, double I, int idx)

	// Input		Q	: Quadrupolar interaction
	//  			pset	: Parameter set
	//			I       : Quantum I value (0.5, 1, 3.5, ...)
	//			idx	: Tensor index (default -1)
	// Output		TF      : Q is set for a spin with quantum
	//				  number I by the parameters in the
	//			          parameter set pset.  Return is FALSE
	//				  if proper parameters are not found
	// Note				: Normally, the += function would
	//				  suffice, but here we support the 
	//				  ability to specify an index in
	//				  the parameter name(s) used

// sosix - must be careful about this being public! If so it must check that I
//         is valid (>1/2).  Also, one must make sure the arrays are there too.

// There are several ways to define a quadrupolar tensor for a specific spin.
// The most favorable (looked for initially) is in terms of a tensor specified
// by PAS (principal axis) components: { Aiso, delzz, eta, phi, theta, gamma }
// However Aiso = gamma = 0 in all quadrupolar tensors so these remain unused.
// An example of the parameter ASCII lines would be (removing the starting //)

//AQ	(4) : 2					- Quadrupolar Tensor (kHz)
//            ( 0.00,  1500.0,  0.25)
//            (  0.00,   0.00,   0.00)

// where the quadrupolar coupling is 1.5 MHz and the asymmetry is 0.25.
// Note - default units (for delzz) of the AQ parameter (quad tensor) are KHz!

// The next definition scanned for will be the tensor in terms of Cartesian
// components.  In that event there must be 9 componentsmin thetype
// Finally, a third option is available when one only needs a primitive tensor.
// In this case one simply defines the Quadrupolar coupling constant and an

  {
//sosi - needs work?
I += 1;
  SinglePar par;				// A single parameter
  ParameterSet::const_iterator item;         // A pix into parameter list
  string pname;					// For parameter name
  string pstate;				// For parameter statement
  string sval;					// For parameter value
  string Nidx = "";				// Name addition per index 
  if(idx >= 0) 					// If index exists then set up
    Nidx += string("(")+Gdec(idx)+string(")");	// the parameter name append

// --------------- Try & Read AQ: PAS Spatial Quadrupolar Tensor -------------

  int Qflag = 0;				// Flag if tensor found
  SinglePar Qpar;				// Quad. tensor parameter 
  string Qname = string("AQ") + Nidx; 		// Quad. tensor parameter name
  Qpar = SinglePar(Qname);			//	Parameter for this name
// sosi this is off now
//  if(pset.seek(Qpar))				// 	If parameter is in pset
//    {						//	then copy it into sys
//    item = pset.seek(Qname);			//	Here's where it is
//    *(this) = space_T(*item);			// 	Set quad. tensor
//    Qflag++;					//	Flag that we know it
//    }	 					// Note QCC here read in kHz!
  
// -------------- Try & Read Q_T: PAS Spatial Quadrupolar Tensor -------------
//                  (Old Parameter Name, Marked For Deletion!)

  if(!Qflag)					// Next try another parameter
    {						// & seek Quad spatial tensor
    Qname = string("Q_T") +  Nidx;	 	// Quad. tensor parameter name
    item = pset.seek(Qname);			//	Here's where it is
    if(item != pset.end())				// 	If parameter is in pset
      {						//	then copy it into sys
cout << "\nHere it is\n" << (*item);
//space_T X(pset(*item));
//cout << "\n\n\tRead Tensor: \n" << X;
//*(this) = X;
cout << "\n\n\tAssigned Tensor: \n" << *(this);
//      *(this) = space_T(*item);		// 	Set quad. tensor
      Qerror(40, 1);
      Qerror(41);
      Qflag++;					//	Flag that we know it
      } 					// Note QCC here read in kHz!
    }
cout << "\n\t\tQflag is " << Qflag;

// --------- Try & Read AQCart: Cartesian Spatial Quadrupolar Tensor ---------

  if(!Qflag)					// Next try another parameter
    {						// & seek Quad spatial tensor
    Qname = string("AQCart") +  Nidx; 		// Quad. tensor parameter name
    item = pset.seek(Qname);			//	Here's where it is
    if(item != pset.end())			// 	If parameter is in pset
      {						//	then copy it into sys
      Qflag++;					//	Flag that we know it
      cout << "\nClass IntQuad: Don't Yet "
           << "Know How To Read Cartesian " 
           << "Tensors!";
      }
    }

// --- Try & Read {QCC,Qeta,Qphi,Qtheta,Qchi }: PAS Spatial Quad. Tensor ----
//             (QCC Can Also Be QCCk for kHz or QCCM for MHz)

// sosi this is still not done
  if(!Qflag)					// Next try another parameter
    {						// & seek Quad spatial tensors
    Qname = string("QCC") +  Nidx;	 	// Quad. coupling param. name
cout << "\n\t Looking for " << Qname;
    item = pset.seek(Qname);			//	Here's where it is
    if(item != pset.end())			// 	If parameter is in pset
      {						//	then copy it into sys
//      (*item).parse(pname,pdatad,pstate);
      cout << "\nClass IntQuad: Don't Yet "
           << "Know How To Handle This!"; 
      }
      Qflag = 1;				//	Flag that we know it
    }
   return Qflag;
   } 


// ----------------------------------------------------------------------------
//      Functions To Make A Quadrupolar Interaction From A Parameter Set
// ----------------------------------------------------------------------------
 

int IntQuad::setAsph(const ParameterSet& pset, int idx)
 
	// Input		Q	: Quadrupolar interaction (this)
	// 			pset	: A parameter set
	//			idx	: Index value
        // Output               none    : Quadrupolar interaction
	//				  set from parameters in pset
 
/*---------------------------------------------------------------------------- 
**                Try & Read AQ Quadrupolar Coupling Spatial Tensor	    **
**  			   { delzz, eta, phi, theta }                       **
**									    **
** AQ	(2) : ( QCC, eta, theta, phi )		- Quadrupolar Interaction   **
**									    **
*****************************************************************************/ 
 
// Look for AQ parameter. Currently allowed parameters are the following:
//
//  AQ,    AQ(#) 	- Quadrupolar Coupling Specified in kHz 
//  AQkHz, AQkHz(#)	- Quadrupolar Coupling Specified in kHz
//  AQKHz, AQKHz(#)	- Quadrupolar Coupling Specified in kHz
//  AQHz,  AQHz(#)	- Quadrupolar Coupling Specified in Hz
//  AQMHz, AQMHz(#)	- Quadrupolar Coupling Specified in MHz

  {
  int found=0; 
  string Nidx = "";				// Name addition per index
  if(idx >= 0) Nidx 				// If index exists then set up
    += string("(")+Gdec(idx)+string(")");	// the parameter name to append
  string Qnames[5] = { "AQ", "AQkHz", "AQKHz",	// Quad. interaction parameter
		       "AQMHz", "AQHz" };
  string pname, pstate;				// Strings for parameter
  string val1, val2;
  ParameterSet::const_iterator item;         // A pix into parameter list
  SinglePar par;
  for(int i=0; i<5 && !found; i++) 		// Loop over possible spatial
    {						// tensor parameters names
    pname = Qnames[i] +  Nidx;			// 	Parameter name
    item = pset.seek(pname); 			//      Seek parameter in pset
    if(item != pset.end())			//	If it's been found
      {						//	parse the parameter
// sosi - The Pset parsing is goofing up String parameters if they have
//        spaces in them.......  Is that correct?
/*
sosi - This uses older style Regex, so it needs updating
cout << pset;
cout << (*item);
cout.flush();
      val1 = par.data();			//	2nd get the values
      cut(val1,"([ \t]*");			//	Remove "(" & blanks
cout << val1;
cout.flush();
      val2 = cut(val1,RXdouble);		//	Cut out the QCC value
      _QCC = atof(val2);			//	Convert to double
      cut(val1,"[ \t]*,[ \t]*");		//	Remove blanks,blanks
cout << val1;
cout.flush();
      val2 = cut(val1,RXdouble);		//	Cut out the ETA value
      ETA = atof(val2);				//	Convert to double
      cut(val1,"[ \t]*,[ \t]*");		//	Remove blanks,blanks
      val2 = cut(val1,RXdouble);		//	Cut out the THETA value
      THETA = atof(val2)*DEG2RAD;		// 	Convert to double
      cut(val1,"[ \t]*,[ \t]*");		//	Remove blanks,blanks
      val2 = cut(val1,RXdouble);		//	Cut out the PHI value
      PHI = atof(val2)*DEG2RAD;			// 	Convert to double
      switch(i)					//	based on input type
        {
        case 0:					//	Here if QCC in kHz
        case 1:
        case 2:
        default:
	  _QCC *= 1.e3; break;
        case 3:					// 	Here if QCC in MHz
	  _QCC *= 1.e6; break;
          break;
        case 4:					// 	Here if QCC in Hz
          break;
        }
      found++;
*/
      }
    }
  return found;
  }
 

/*---------------------------------------------------------------------------- 
**                Try & Read Q_T GENERIC (QUAD) Spatial Tensor		    **
**			  (This Should Be Avoided!!!!)			    **
**									    **
** Q_T	(4) : 2					- Quadrupolar Tensor (kHz)  **	
**            ( x.xx,  DCC(KHz),  eta)					    **
**            ( theta, phi,       x.xx)					    **
**									    **
*****************************************************************************/ 

// The parameter Q_T is a quadrupolar spatial tensor set as a generic rank 
// 2 (not irreducible) spatial tensor.  In order to read this easily I would
// have to include class space_T, the very class I am bypassing with the use
// of rank 2 interactions!  So, I will allow use of Q_T to set up IntQuad BUT
// I'll issue a warning AND (ugh!) I'll have redefine the parameter parse.


int IntQuad::setAsphGen(const ParameterSet& pset, int idx)
 
	// Input		Q	: Quadrupolar interaction (this)
	// 			pset	: A parameter set
	//			idx	: Index value
        // Output               none    : Quadrupolar interaction
	//				  set from parameters in pset
// sosi - CAREFUL! There may still be an inconsistency in the mapping of the
//        Euler angles alpha,beta,gamma into the two spherical angles theta,phi

  {
  int found=0; 
  string Nidx = "";				// Name addition per index
  if(idx >= 0) Nidx 				// If index exists then set up
    += string("(")+Gdec(idx)+string(")");	// the parameter name to append
  string Qnames = "Q_T"; 			// Quad. interaction parameter
  string pname, pstate;				// Strings for parameter
  string val1, val2;
  ParameterSet::const_iterator item;         // A pix into parameter list
  SinglePar par;
  pname = Qnames +  Nidx;			// Parameter name
  par = SinglePar(pname);			// Parameter for name
  item = pset.seek(pname); 			// Seek parameter in pset
  if(item != pset.end())				// If it's been found
    {						// parse the parameter
    Qerror(40,1);				// 	Warn outdated parameter
    Qerror(41);					// 	Suggest AQ not Q_T
/*
sosi - uses older style Regex so need to be updated
    par = pset(item); 				// 	1st get the parameter
    if(par.type() != 4) Qfatal(42); 		// 	Type must be space_T
    val1 = par.data();				// 	Get the parameter data
    val2 = cut(val1,RXint);			// 	Cut integer - rank
    if(atoi(val2) != 2) Qfatal(43);		// 	Must be rank 2 tensor
    cut(val1,"[ \t]*,[ \t]*");			// 	Cut blanks or commas
    cut(val1,"([ \t]*");			// 	Cut any ( and blanks
    val2 = cut(val1,RXdouble);			// 	Cut double value - Qiso
    if(atof(val2) != 0) Qerror(44);		//	Ignore isotropic value
    cut(val1,"[ \t]*,[ \t]*"); 			// 	Cut blanks or commas
    val2 = cut(val1,RXdouble);			//	Cut out the QCC value
    _QCC = atof(val2)*1.e3;			//	Convert to double (kHz)
    cut(val1,"[ \t]*,[ \t]*");			//	Remove blanks,blanks
    val2 = cut(val1,RXdouble);			//	Cut out the ETA value
    ETA = atof(val2);				//	Convert to double
    cut(val1,"[ \t]*)");			// 	Cut any blanks or )
    cut(val1,"[ \t]*,[ \t]*");			// 	Cut any blanks or commas
    cut(val1,"([ \t]*"); 			//      Cut any ( or blanks
    val2 = cut(val1,RXdouble);			//	Cut out the THETA value
    THETA = atof(val2)*DEG2RAD;			// 	Convert to double (deg)
    cut(val1,"[ \t]*,[ \t]*");			//	Remove blanks,blanks
    val2 = cut(val1,RXdouble);			//	Cut out the PHI value
    PHI = atof(val2)*DEG2RAD;			// 	Convert to double (deg)
*/
    found++;
    }
  return found;
  }

// ----------------------------------------------------------------------------
//        Interactive Read of Q Intraction From An ASCII File Or Input
// ----------------------------------------------------------------------------

string IntQuad::ask_read(int argc, char* argv[], int argn, int idx)
 
        // Input                Q       : Quadrupolar interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
        //                      idx     : Index for spin
        // Output               String  : The parameter argn of array argc is
        //                                used to supply a filename from which
        //                                the quadrupolar interaction is read
        //                                If argument argn is not in argv, the
        //                                user is asked to supply a filename
        //                                The set filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing known IntQuad parameters
        // Note                         : The interaction Q is modifed (filled)
	// Note				: Since no spin type is given, this will
	//				  also use the parameter Iso(#) to
	//				  the value of the spin quantum number

  {             
  string filename;                              // Name of parameter file
  query_parameter(argc, argv, argn,             // Get filename from command
    "\n\tQuadrupolar Interaction filename? ",	// Or ask for it
                                      filename);
  read(filename, idx);				// Read system from filename
  return filename;
  }


void IntQuad::ask(int argc, char* argv[], int& qn, double& QI,
          double& Qnqcc, double& Qeta, double& Qtheta, double& Qphi, int Qflag)

        // Input                Q	: Quadrupolar interaction (this)
	//			argc    : Number of arguments
	//			argv    : Array of arguments
	//			qn      : Query value
	//			QI	: Spin quantum number
	//			Qnqcc   : Quad. coupling constant (Hz)
	//			Qeta    : Quad. asymmetry
	//			Qtheta  : Quad. orientation angle
	//			Qphi	: Quad. orientation angle
        //                      Qflag   : Flag is QCC or wQ requested
        // Output               none    : The values of qn, I, Qnqcc, Qeta,
	//				  Qtheta, and Qphi are set herein
	// Note				: This is INTERACTIVE!

  {
  query_parameter(argc, argv, qn++,			// Read in the I value
      "\n\tSpin Quantum Value (1, 1.5, ..)? ", QI);
  if(Qflag)
    {
    query_parameter(argc, argv, qn++,			// Read the frequency
       "\n\tQuadrupolar Frequency(kHz)? ", Qnqcc);
// sosik
double I = Izval();
    Qnqcc *= 2.0*I*(2.0*I-1.0)/3.0;			// Set quad. coupling
    }
  else
    {
    query_parameter(argc, argv, qn++,			// Read in the coupling
       "\n\tQuadrupolar Coupling (kHz)? ", Qnqcc);
    }
  Qnqcc *= 1.e3;					// Put this in Hz
  query_parameter(argc, argv, qn++,			// Read in the coupling
       "\n\tQuadrupolar Asymmetry [0, 1]? ", Qeta);
  query_parameter(argc, argv, qn++,			// Read in theta angle
  "\n\tOrientation Down From PAS z [1, 180]? ", Qtheta);
  query_parameter(argc, argv, qn++,			// Read in phi angle
   "\n\tOrientation Over From PAS x [0, 360]? ", Qphi);
  }


void IntQuad::askset(int argc, char* argv[], int& qn, int Qflag)

        // Input                Q	: Quadrupolar interaction (this)
	//			argc    : Number of arguments
	//			argv    : Array of arguments
	//			qn      : Query value
        //                      Qflag   : Flag is QCC or wQ requested
        // Output               none    : Q is set interactively
	// Note				: This is INTERACTIVE!

  {
  double QI, Qnqcc, Qeta, Qtheta, Qphi;
  ask(argc,argv,qn,QI,Qnqcc,Qeta,Qtheta, Qphi, Qflag);	// Use the ask function
// sosik
//  *(this) = IntQuad(QI,Qnqcc,Qeta,Qtheta,Qphi); 	// Use assignment
  }

 
void IntQuad::askset(int Qflag)
 
        // Input                Q       : Quadrupolar interaction (this)
        //                      Qflag   : Flag is QCC or wQ requested
        // Output               none    : Q is set interactively
        // Note                         : This is INTERACTIVE!
 
  {
  int qn = 1000;
  int argc = 1;
  char* argv[1];
  double QI, Qnqcc, Qeta, Qtheta, Qphi;
  ask(argc,argv,qn,QI,Qnqcc,Qeta,Qtheta, Qphi,Qflag);	// Use the ask function
// sosik
//  *(this) = IntQuad(QI,Qnqcc,Qeta,Qtheta,Qphi);         // Use assignment
  }
 


// ____________________________________________________________________________
// J                 QUADRUPOLAR INTERACTION OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//  Functions That Generate String Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

/*  Function   Arguments                     Output
    ========   =========   ===================================================
    TStrings       m       String array for the mth component of spin tensor T
    QAStrings              String array for various interaction values

            QAStrings()                              TStrings(m)

     Isotropic G value:      x.xxxxxxx                [ x.x, x.x, x.x]
     G Anisotropy:           x.xxxxxxx        T     = [ x.x, x.x, x.x]
     G Asymmetry:            x.xxxxxxx         2,m    [ x.x, x.x, x.x]
     Down From PAS z-Axis:   x.xx Deg.                [ x.x, x.x, x.x]
     Over From PAS x-Axis:   x.xx Deg.
     Applied Field (kG):     x.xxxxxxx
                                              m = [0,4] => {0,1,-1,2,-2}     */

// string* IntRank2T::TStrings(int M) const;                     INHERITED

vector<string> IntQuad::CartAStrings(const string& CSForm) const
  {
  vector<string> Cartstrings(6);        	// Vector of 6 strings
  vector<string> Cstrs = CAStrings("q");
  double Vxx = qxx()*1.e-3;;
  double Vxy = qxy()*1.e-3;;
  double Vxz = qxz()*1.e-3;;
  double Vyy = qyy()*1.e-3;;
  double Vyz = qyz()*1.e-3;;
  double Vzz = qzz()*1.e-3;;
  Cartstrings[0] = Cstrs[0];
  Cartstrings[1] = Cstrs[1] + string("   ") + AuvString(Vxx, Vxy, Vxz, CSForm);
  Cartstrings[2] = Cstrs[2] + string(" = ") + AuvString(Vxy, Vyy, Vyz, CSForm);
  Cartstrings[3] = Cstrs[3] + string(" = ") + AuvString(Vxz, Vyz, Vzz, CSForm);
  Cartstrings[4] = Cstrs[4];
  Cartstrings[5] = Cstrs[5];

  int nstr = 6;					// Number of strings used
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

vector<string> IntQuad::SphAStrings() const
  {
  vector<string> Sphstrings(6);
  int k=0;
  Sphstrings[k++] = string("Frequency (kHz):   ")
                  + Gform("%10.2f", wQ()*1.e-3)
                  + string("     "); 
  Sphstrings[k++] = StringI();
  Sphstrings[k++] = string("Nuclear QCC (kHz): ")
                  + Gform("%10.2f", _QCC*1.e-3)
                  + string("     "); 
  Sphstrings[k++] = AsymmetryString();
  Sphstrings[k++] = ThetaString();
  Sphstrings[k++] = PhiString();
  return Sphstrings;
  }

//-----------------------------------------------------------------------------
//    Functions That Generate Ouput Of The Rank 2 Quadrupolar Interaction
//-----------------------------------------------------------------------------

/* These functions will output information concerning the Quadrupolar 
   interaction to any output stream.

        // Input                Q	: Quadrupolar interaction (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Basic Parameters
                                           !0 - Full output
                                nrm     : Flag if GAMMA normalized output
           Output               none    : Quad interaction parameters
                                          placed into the output stream      */

//ostream& IntQuad::print(ostream& ostr, int fflag, int nrm) const
ostream& IntQuad::print(ostream& ostr, int fflag) const
  {
  if(!_QCC)
    {
    ostr << "\n\n\t\tEmpty Quadrupolar Interaction\n";
    return ostr;
    }

//	    Output Some Printed Lines Which Will Look Like The Following 

//			       Quadrupolar Interaction
//
// Frequency:             xxxxx.xx xHz       [A  , A  , A  ] 
// Spin Quantum Number:       I              [ xx   xy   xz]   [ x.x, x.x, x.x]
// Coupling Constant:     xxxxx.xx xHz       [A  , A  , A  ] = [ x.x, x.x, x.x]
// Asymmetry:                 x.xx           [ yx   yy   yz]   [ x.x, x.x, x.x]
// Down From PAS z-Axis:    xxx.xx Degrees   [A  , A  , A  ]
// Over From PAS x-Axis:    xxx.xx Degrees   [ zx   zy   zz]
//

  vector<string> Astrs = CartAStrings("%6.3f");
  vector<string> Qstrs = SphAStrings();
  ostr << "\t\n\t\t\t\tQuadrupolar Interaction\n";
  int k = 0;
  int j = 0;
  ostr << "\n" << Qstrs[j++] << "   " << Astrs[k++];
  ostr << "\n" << Qstrs[j++] << "   " << Astrs[k++];
  ostr << "\n" << Qstrs[j++] << "   " << Astrs[k++];
  ostr << "\n" << Qstrs[j++] << "   " << Astrs[k++];
  ostr << "\n" << Qstrs[j++] << "   " << Astrs[k++];
  ostr << "\n" << Qstrs[j++] << "   " << Astrs[k++];
  if(!fflag) return ostr;

/*          Now Output The Spherical Components, Spatial And Spin
       These Will Be Printed Lines Which Will Look Like The Following
                        (Repeated For All 5 m Values)

                         Xi = x.xxx GHz = x.xxx G

                A    = x.xxx            T     = [ x.x, x.x]
                 2,m                     2,m    [ x.x, x.x]                  */

  ostr << "\n\n" << string(24, ' ')
       << " Xi = "  << _XI/1.e9
       << " GHz = " << _XI/1.3996e6 << " G";
  printAT(ostr, G);
  ostr << "\n\n";
  return ostr;
  }

ostream& operator<< (ostream& out, const IntQuad& Q) { return Q.print(out); }

//-----------------------------------------------------------------------------
//  Functions That Generate Ouput Of Cartesian and Spherical & Cartesian A
//-----------------------------------------------------------------------------

/*
ostream& IntQuad::printSpherical(ostream& ostr)
         

  {
  ostr << "\t\n\t         Quadrupolar Spatial Tensor";
  if(PAS()) ostr << "\n(Oriented Along The Principal Axes)";
  else
    {
    ostr << "\n\t     "
         << Gform("%7.2f", THETA*RAD2DEG)
         << " Degrees Down From PAS z";
    if(PHI)
      ostr << "\n\t     "
         << Gform("%7.2f", PHI*RAD2DEG)
         << " Degrees Over From PAS x";
     }
  ostr << "\n\n\tA    = " << Asph[0];
  ostr <<   "\n\t 2,0";
  ostr << "\n\n\tA    = " << Asph[1];
  ostr <<   "\n\t 2,1";
  ostr << "\n\n\tA    = " << Asph[2];
  ostr <<   "\n\t 2,-1";
  ostr << "\n\n\tA    = " << Asph[3];
  ostr <<   "\n\t 2,2";
  ostr << "\n\n\tA    = " << Asph[4];
  ostr <<   "\n\t 2,-2\n";
  return ostr; 
  }

 
ostream& IntQuad::printCartesian(ostream& ostr)
   
        // Input                Q       : Quad spatial tensor (this)
        //                      ostr    : Output stream
        // Output               none    : Quad spatial tensor parameters set to
        //                                output stream

// Axx = (A22+A2m2)/2 - A20/sqrt(6)      Axy = -i*(A22-A21)/2 = Ayx
// Ayy = -(A22+A2m2)/2 - A20/sqrt(6)     Axz = -(A21-A2m1)/2  = Azx
// Azz = sqrt(2/3)*A20                   Ayz = i*(A21-A2m1)/2 = Azy
 
  {
  double Vxx, Vyy, Vzz;
  double Vxy, Vxz, Vyz;
  Vxx = Re((Asph[4]+Asph[3])/2. - Asph[0]/sqrt(6.));
  Vyy = Re(-(Asph[4]+Asph[3])/2. - Asph[0]/sqrt(6.));
  Vzz = Re(sqrt(2.0/3.0)*Asph[0]); 
  Vxy = Re(complex(0,-0.5) * (Asph[4] - Asph[3]));
  Vxz = Re(- 0.5 * (Asph[1] - Asph[2]));
  Vyz = Re(complex(0,0.5) * (Asph[1] + Asph[2]));
  ostr << "\t\n\t         Quadrupolar Spatial Tensor";
  if(PAS()) ostr << "\n(Oriented Along The Principal Axes)";
  else
    {
    ostr << "\n\t     "
         << Gform("%7.2f", THETA*RAD2DEG)
         << " Degrees Down From PAS z";
    if(PHI)
      ostr << "\n\t     "
         << Gform("%7.2f", PHI*RAD2DEG)
         << " Degrees Over From PAS x";
     }
  ostr << "\n\n\t[A  , A  , A  ]   ";
  ostr << "\n\t[ xx   xy   xz]   ";
  ostr << "["  << Gform("%6.3f", Vxx)
       << ", " << Gform("%6.3f", Vxy)
       << ", " << Gform("%6.3f", Vxz) << "]";
  ostr << "\n\t[A  , A  , A  ] = "
       << "["  << Gform("%6.3f", Vxy)
       << ", " << Gform("%6.3f", Vyy)
       << ", " << Gform("%6.3f", Vyz) << "]";
  ostr << "\n\t[ yx   yy   yz]   "
       << "["  << Gform("%6.3f", Vxz)
       << ", " << Gform("%6.3f", Vyz)
       << ", " << Gform("%6.3f", Vzz) << "]";
  ostr << "\n\t[A  , A  , A  ]   "
       << "\n\t[ zx   zy   zz]   \n";
  return ostr;
  }

 
ostream& IntQuad::printCartesian(ostream& ostr, double phi, double theta)
         
        // Input                Q	: Quad spatial tensor (this)
        //                      ostr	: Output stream
        //                      phi     : Orientation angle (degrees)
        //                      theta   : Orientation angle (degrees)
        // Output               none	: Quad spatial tensor parameters set to
        //                                output stream

// Axx = (A22+A2m2)/2 - A20/sqrt(6)      Axy = -i*(A22-A21)/2 = Ayx
// Ayy = -(A22+A2m2)/2 - A20/sqrt(6)     Axz = -(A21-A2m1)/2  = Azx
// Azz = sqrt(2/3)*A20                   Ayz = i*(A21-A2m1)/2 = Azy
 
  {
  double Vxx, Vyy, Vzz;
  double Vxy, Vxz, Vyz;
  Vxx = Re((Asph[4]+Asph[3])/2. - Asph[0]/sqrt(6.));
  Vyy = Re(-(Asph[4]+Asph[3])/2. - Asph[0]/sqrt(6.));
  Vzz = Re(sqrt(2.0/3.0)*Asph[0]);
  Vxy = Re(complex(0,-0.5) * (Asph[4] - Asph[3]));
  Vxz = Re(- 0.5 * (Asph[1] - Asph[2]));
  Vyz = Re(complex(0,0.5) * (Asph[1] + Asph[2]));
  double Vxxr, Vyyr, Vzzr;
  double Vxyr, Vxzr, Vyzr;
  if(phi || theta)
    {
    Vxxr = Axx(theta, phi);
    Vyyr = Ayy(theta, phi);
    Vzzr = Azz(theta, phi);
    Vxyr = Axy(theta, phi);
    Vxzr = Axz(theta, phi);
    Vyzr = Ayz(theta, phi);
    }
  ostr << "\t\n\t    Quadrupolar Spatial Tensor";
  if(!phi && !theta) ostr << "\n(Oriented in Its PAS)";
  else
    {
    ostr << "\n(Oriented"
         << Gform("%7.2f", theta) << " Degrees Down From PAS z)";
    if(phi)
      ostr << " And Over "
         << Gform("%7.2f", phi) << " Degrees From PAS x";
     }
  ostr << "\n\n\t[A  , A  , A  ]   ";
  ostr << "\n\t[ xx   xy   xz]   ";
  ostr << "["  << Gform("%6.3f", Vxx)
       << ", " << Gform("%6.3f", Vxy)
       << ", " << Gform("%6.3f", Vxz) << "]";
  ostr << "\n\t[A  , A  , A  ] = "
       << "["  << Gform("%6.3f", Vxy)
       << ", " << Gform("%6.3f", Vyy)
       << ", " << Gform("%6.3f", Vyz) << "]";
  ostr << "\n\t[ yx   yy   yz]   "
       << "["  << Gform("%6.3f", Vxz)
       << ", " << Gform("%6.3f", Vyz)
       << ", " << Gform("%6.3f", Vzz) << "]";
  ostr << "\n\t[A  , A  , A  ]   "
       << "\n\t[ zx   zy   zz]   \n";
  return ostr; 
  }



ostream& IntQuad::printCartesian(ostream& ostr)

        // Input                QI      : Quadrupolar interaction (this)
        //                      ostr    : Output stream
        // Output               none    : Q spatial tensor parameters
        //                                sent to output stream
        // Note                         : Uses base class virtual overload

// Axx = (A22+A2m2)/2 - A20/sqrt(6)      Axy = -i*(A22-A21)/2 = Ayx
// Ayy = -(A22+A2m2)/2 - A20/sqrt(6)     Axz = -(A21-A2m1)/2  = Azx
// Azz = sqrt(2/3)*A20                   Ayz = i*(A21-A2m1)/2 = Azy

  {
  ostr << "\n\tQuadrupolar Interaction Cartesian Spatial Tensor\n";
  IntRank2A::printCartesian(ostr, 0);
  return ostr;
  }

ostream& IntQuad::printCartesian(ostream& ostr, double phi, double theta)

        // Input                QI      : Quadrupolar interaction (this)
        //                      ostr    : Output stream
        //                      phi     : Orientation angle (degrees)
        //                      theta   : Orientation angle (degrees)
        // Output               none    : Q spatial tensor parameters
        //                                sent to output stream
        // Note                         : Uses base class virtual overload

// Axx = (A22+A2m2)/2 - A20/sqrt(6)      Axy = -i*(A22-A21)/2 = Ayx
// Ayy = -(A22+A2m2)/2 - A20/sqrt(6)     Axz = -(A21-A2m1)/2  = Azx
// Azz = sqrt(2/3)*A20                   Ayz = i*(A21-A2m1)/2 = Azy

  {
  ostr << "\n\tQuadrupolar Interaction Cartesian Spatial Tensor\n";
  IntRank2A::printCartesian(ostr, theta, phi, 0);
  return ostr;
  }

ostream& IntQuad::printCartQ(ostream& ostr, int tflag) const

        // Input                QI      : Quadrupolar interaction (this)
        //                      ostr    : Output stream
        // Output               none    : Full Quadrupolar Cartiesian
	//				  spatial tensor sent to output stream

//                     [q  , q  , q  ]
//                     [ xx   xy   xz]   [ x.x, x.x, x.x]
//                     [q  , q  , q  ] = [ x.x, x.x, x.x]
//                     [ yx   yy   yz]   [ x.x, x.x, x.x]
//                     [q  , q  , q  ]
//                     [ zx   zy   zz]
//

  {
  if(tflag)
    ostr << "\n\n" << string(30, ' ')
         << "Quadrupolar Interaction Cartesian Spatial Tensor\n";
  vector<string> Astrs = CartAStrings("%6.3f");
  string newgline = string("\n") + string(20, ' ');
  for(int i=0; i<int(Astrs.size()); i++)
    ostr << newgline << Astrs[i];
  return ostr;
  }
*/
#endif							// IntQuad.cc
