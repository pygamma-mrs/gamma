/* IntRank2A.cc *************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Rank 2 Spatial Tensor			Implementation		**
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
** Class IntRank2A represents a very primitive rank 2 spatial tensor.	**
** Herein the tensor is restructed to have the following properties:	**
**                                                                      **
** 1.) rank=2      2.) traceless     3.) symmetric     4.) eta = [0, 1] **
**                                                                      **
** It is stored in the PAS with a GAMMA normalized delzz value. This	**
** is done because the tensor provided by this class is interaction	**
** independent, the scaling used has been chosen so that the spherical	**
** components are directly associated with normalized spherical 	**
** harmonics.  As such, this leaves the asymmetry, eta, as the only 	**
** value used to entirely  describe the	spatial tensor!			**
**									**
** Note: No rank zero terms (isotropic components) are maintained!	**
**       No rank one terms (antisymmetric components) are maintained!	**
**									**
** The asymmetry, eta, is independent of overall tensor scaling.	**
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
** Spatial tensors are used to describe the spatial part of GAMMA	**
** rank 2 interactions: See Class IntRank2.  These full interactions	**
** will be scaled to represent the interaction strength.		**
**            								**
*************************************************************************/

#ifndef IntRank2A_cc_			// Is file already included?
#  define IntRank2A_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <IntRank2/IntRank2A.h>		// Include interface definition
#include <IntRank2/CartMx2A.h>		// Include A conversion to PAS
#include <Basics/Gutils.h>		// Include GAMMA standared errors
#include <Basics/Gconstants.h>		// Include PI and other constants
#include <Basics/Isotope.h>		// Include isotopes
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Matrix/complex.h>		// Include complex numbers
#include <Matrix/row_vector.h>		// Include row_vectors
#include <Matrix/matrix.h>		// Include matrices
#include <Basics/Gutils.h>		// Include query_parameter & errors
#include <Level1/coord.h>		// Include coordinates
#include <Basics/StringCut.h>		// Include GAMMA string parsing

#if defined(_MSC_VER)			// If we are using MSVC++
 #pragma warning (disable : 4786)       //   Kill STL namelength warnings
#endif

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ STL list
using std::ofstream;			// Using libstdc++ output file streams
using std::vector;			// Using libstdc++ STL vectors
using std::ostream;			// Using libstdc++ output streams

const double RT6PIO5 = sqrt(6.*PI/5.); 		// Constant sqrt[6*PI/5]
const double RT5O4PI = sqrt(5./(4.*PI)); 	// Constant sqrt[5/(4*PI)];
const double RT5O6PI = sqrt(5./(6.*PI)); 	// Constant sqrt[5/(6*PI)];
const double RT5O24PI = sqrt(5./(24.*PI)); 	// Constant sqrt[5/(24*PI)];
const double RT5O96PI = 0.5*RT5O24PI; 		// Constant sqrt[5/(96*PI)];

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i           RANK 2 INTERACTION SPATIAL TENSOR ERROR HANDLING
// ____________________________________________________________________________

/*       Input                IR2A    : Rank 2 spatial tensor (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

void IntRank2A::IR2Aerror(int eidx, int noret) const
  {
  string hdr("Rank 2 Spatial Tensor");
  switch(eidx)
    {
    case 2: GAMMAerror(hdr,"Problems During Construction", noret); break;// (2)
    case 3: GAMMAerror(hdr,"Problems During Assignmant",   noret); break;// (3)
    case 4: GAMMAerror(hdr,"Nonzero Isotropic Cmp., Zero!",noret); break;// (4)
    case 5: GAMMAerror(hdr,"Isotropic (l=0)Comp., Zeroed!",noret); break;// (5)
    case 6: GAMMAerror(hdr,"AntiSymmetric (l=1) Zeroing!", noret); break;// (6)
    case 8: GAMMAerror(hdr,"Theta (z Down) Beyond [0,180]",noret); break;// (8)
    case 9: GAMMAerror(hdr,"Phi (x Over) Outside [0, 360]",noret); break;// (9)
    case 10:GAMMAerror(hdr,"Asymmetry (eta) Beyond [0, 1]",noret); break;// (10)
    case 11:GAMMAerror(hdr,"Setting Asymmetry to Zero",    noret); break;// (11)
    case 12:GAMMAerror(hdr,"Asymmetry Set In Zero Tensor!",noret); break;// (12)
    case 13:GAMMAerror(hdr,"Cannot Construct From File",   noret); break;// (13)
    case 14:GAMMAerror(hdr,"Cartesian component order bad",noret); break;// (14)
    case 15:GAMMAerror(hdr,"Ordering: |Azz|>=|Ayy|>=|Axx|",noret); break;// (15)
    case 18:GAMMAerror(hdr,"Cannot Alter Interaction",     noret); break;// (18)
    case 19:GAMMAerror(hdr,"Cannot Set Freq., No I Value", noret); break;// (19)
    case 20:GAMMAerror(hdr,"Cannot Write To Param. File",  noret); break;// (20)
    case 21:GAMMAerror(hdr,"Cannot Read From Param. File", noret); break;// (21)
    case 22:GAMMAerror(hdr,"Can't Write To Output FileStr",noret); break;// (22)
    case 23:GAMMAerror(hdr,"CAnnot Ouptut Parameters",     noret); break;// (23)
    case 25:GAMMAerror(hdr,"Invalid Cartesian Array",      noret); break;// (25)
    case 26:GAMMAerror(hdr,"Array Must Be 3x3 Symmetric",  noret); break;// (26)
    case 27:GAMMAerror(hdr,"Array Must Be Traceless",      noret); break;// (27)
    case 40:GAMMAerror(hdr,"Use Of Outdated Parameter!",   noret); break;// (40)
    case 41:GAMMAerror(hdr,"Use AQ Over Q_T Please",       noret); break;// (41)
    case 50:GAMMAerror(hdr,"Component Index Beyond [-2,2]",noret); break;// (50)
    case 51:GAMMAerror(hdr,"Cant Find Orientation Angles", noret); break;// (51)
    case 52:GAMMAerror(hdr,"Unable to Set Orientation",    noret); break;// (52)
    case 53:GAMMAerror(hdr,"Setting Euler Angles To Zero", noret); break;// (53)

    case 1: GAMMAerror(hdr,"Construct With Rank Not 2!",   noret); break;// (1)
    case 92:GAMMAerror(hdr,"Setting Theta To Zero Degrees",noret); break;// (92)
    case 93:GAMMAerror(hdr,"Setting Phi To Zero Degrees",  noret); break;// (93)

    default:GAMMAerror(hdr, eidx, noret);                          break;
    }
  }

volatile void IntRank2A::IR2Afatal(int eidx) const
  {
  IR2Aerror(eidx, 1);				// Output error message
  if(eidx) IR2Aerror(0);			// Write that its fatal
  GAMMAfatal();					// Clean exit from program
  }

/* Err Ind.     Default Message            Err Ind.     Default Message
   -------- -----------------------------  -------- ---------------------------
      1     Problems With File pname          2     Cannot Read Parameter pname
      3     Invalid Use of Function pname
      4     Use Of Deprecated Funciton pname
      5     Please Use Class Member Funciton pname
   default  Unknown Error - pname                                           */


void IntRank2A::IR2Aerror(int eidx, const string pname, int noret) const
  {
  string hdr("Rank 2 Spatial Tensor");
  string msg;
  switch(eidx)
    {
    case 101:                                                   // (101)
      msg = string("Can't Find Parameters For ")
          + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 111:                                                   // (111)
      msg = string("Cannot Parse Parameter ") + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 113:                                                   // (113)
      msg = string("Can't Find Coordinate For Spin ")
          + pname;
      GAMMAerror(hdr, msg, noret); break;
    case 200:                                                   // (200)
      msg = string("Odd Asymmetry Value Of ")
          + pname + string(" Specified?");
      break;
    case 201:                                                   // (201)
      msg = string("Odd Orientation (+z Down) Angle Of ")
          + pname + string(" Specified?");
      GAMMAerror(hdr, msg, noret); break;
    case 202:                                                   // (202)
      msg = string("Odd Orientation (+x Over) Angle Of ")
          + pname + string(" Specified?");
      GAMMAerror(hdr, msg, noret); break;
    case 203:                                                   // (203)
      msg = string("Cannot Set Asymmetry To ") + pname;
      GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;
    }
  }

// ____________________________________________________________________________
// ii         RANK 2 SPATIAL TENSOR FROM PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

/* In this section we specify how a rank 2 spatial tensor relates to a GAMMA
   parameter file.  That is, these functions determine how a spatial tensor may
   be read/set from an external ASCII file. These are private because some
   parameters are inter-related and hence must be set in a concerted manner. */

// ----------------------------------------------------------------------------
//         Functions To Read The Full Irreducible Rank 2 Spatial Tensor
// ----------------------------------------------------------------------------

// sosi - looks like a logic problem in this function... when is this used?
//        it doesn't appear to be used anywhere........ Remove this then or
//        fix & use in read function if all IntRank classes do not use it.
//        I have made getAX to replace this one, and it is now used internally.
//        The getA function was NOT used.  So, when getA dies rename all getAX
//        to get A and be done with it.

bool IntRank2A::getA(const ParameterSet& pset, const string& A,
                    double& Aeta, EAngles& EA, int iI, int iS, bool warn) const
  {
  coord Aod;					// For Cartesian A off-diags
  coord Axyz;					// For Cartesian A diagonals
  int nod = getAoffdiag(pset,A,Aod,iI,iS,warn);	// Try and read A off-diags
  int nd  = getAxAyAz(pset,A,Axyz,iI,iS,warn);	// Try and read A diagonals
  if(nod)					// If we read some off-diags
    {						// we try to parse for PAS
    matrix Amx(3,3, h_matrix_type);		//   Array for Cartesian mx
    Amx.put(Axyz.x(),  0, 0);			//   Set Axx value in mx
    Amx.put(Axyz.y(),  1, 1);			//   Set Ayy value in mx
    Amx.put(Axyz.z(),  2, 2);			//   Set Azz value in mx
    Amx.put_h(Aod.x(), 0, 1);			//   Set Axy value in mx 
    Amx.put_h(Aod.y(), 0, 2);			//   Set Axz value in mx
    Amx.put_h(Aod.z(), 1, 2);			//   Set Ayz value in mx
    CartMx2A Convert(Amx);			//   Try & covert into PAS
    Aeta = Convert.Eta();			//   Get the asymmetry
    EA   = Convert.EulerAngles();		//   Get the Euler angles
    }

  if(nd)					// If {Axx,Ayy,Azz} present
    {						// parse eta, try orientation
    coord Aide = AisoDelzEta(Axyz);		//   Parse for eta value
    Aeta = Aide.z();				//   Set the return eta value
    getOrientation(pset,A,EA,iI,iS,warn);	//   Try & set Euler angles
    return true;				//   We have found A values
    }
    
  getAeta(pset, A, Aeta, iI, iS, warn);		// Try for eta value directly
  getOrientation(pset,A,EA,iI,iS,warn);		// Try & set Euler angles
  return true;
  }

bool IntRank2A::getAX(const ParameterSet& pset, const string& A,
                 double& Aeta, EAngles& EA, int idxI, int idxS, int warn) const
  {

//                 First Try Using Any Cartesian Parameters

  coord Aod;						// Cartesian A offdiags
  coord Axyz;						// Cartesian A diags
  bool wf = warn>1?true:false;
  int nod = getAoffdiag(pset,A,Aod, idxI,idxS,wf);	// Read A off-diags
  int nd  = getAxAyAz(  pset,A,Axyz,idxI,idxS,wf);	// Read A diagonals
  if(nod)						// If any offdiags read
    {							// try to parse for PAS
    matrix Amx(3,3, h_matrix_type);			//   Cartesian mx
    Amx.put(Axyz.x(),  0, 0);				//   Set Axx in mx
    Amx.put(Axyz.y(),  1, 1);				//   Set Ayy in mx
    Amx.put(Axyz.z(),  2, 2);				//   Set Azz in mx
    Amx.put_h(Aod.x(), 0, 1);				//   Set Axy in mx 
    Amx.put_h(Aod.y(), 0, 2);				//   Set Axz in mx
    Amx.put_h(Aod.z(), 1, 2);				//   Set Ayz in mx
    CartMx2A Convert(Amx);				//   Covert into PAS
    Aeta = Convert.Eta();				//   Get the asymmetry
    EA = Convert.EulerAngles();				//   Get the Euler angles
    return true;					//   We are then done
    }

  if(nd)						// If any offdiags read
    {							// but no off-diagonals
    coord Aide = AisoDelzEta(Axyz);			//   Parse for eta
    Aeta = Aide.z();					//   Set eta value
    getOrientation(pset,A,EA,idxI,idxS,wf);		//   Try for Euler angs
    return true;					//   We are then done
    }

//        No Cartesian Parameters, So We Try For Spherical Parameters

  bool TF = getAeta(pset,A,Aeta,idxI,idxS,wf); 	// Try & read eta value
  if(!TF || !CheckEta(Aeta,wf))				// It must exist & be in
    { 							// the range from [0,1]
    if(warn) IR2Aerror(11,1);				//   Setting eta to 0
    Aeta = 0.0;						//   Set default 0
    return false;					//   Return we failed
    }
  getOrientation(pset,A,EA,idxI,idxS,wf); 		// Last try for orient.
  return true;
  }


bool IntRank2A::getACart(const ParameterSet& pset, const string& A,
                    coord& Aize, EAngles& EA, int iI, int iS, bool warn) const
  {
  coord Aod;					// For Cartesian A off-diags
  coord Axyz;					// For Cartesian A diagonals
  int nod = getAoffdiag(pset,A,Aod,iI,iS,warn);	// Try and read A off-diags
  int nd  = getAxAyAz(pset,A,Axyz,iI,iS,warn);	// Try and read A diagonals
  if(nod)					// If we read some off-diags
    {						// we try to parse for PAS
    matrix Amx(3,3, h_matrix_type);		//   Array for Cartesian mx
    Amx.put(Axyz.x(),  0, 0);			//   Set Axx value in mx
    Amx.put(Axyz.y(),  1, 1);			//   Set Ayy value in mx
    Amx.put(Axyz.z(),  2, 2);			//   Set Azz value in mx
    Amx.put_h(Aod.x(), 0, 1);			//   Set Axy value in mx 
    Amx.put_h(Aod.y(), 0, 2);			//   Set Axz value in mx
    Amx.put_h(Aod.z(), 1, 2);			//   Set Ayz value in mx
    CartMx2A Convert(Amx);			//   Try & covert into PAS
    Aize.x(Convert.Aiso());			//   Set isotropy
    Aize.y(Convert.delzz());			//   Set PAS delzz
    Aize.z(Convert.Eta());			//   Get the asymmetry
    EA   = Convert.EulerAngles();		//   Get the Euler angles
    return true;				//   We succeeded
    }
  if(nd)					// If {Axx,Ayy,Azz} present
    {						// parse eta, try orientation
    Aize = AisoDelzEta(Axyz);			//   Get {Aiso, Adelzz, Aeta}
    getOrientation(pset,A,EA,iI,iS,warn);	//   Try & set Euler angles
    return true;				//   We have found A values
    }
    
  return false;					// Just could not get A from
  }						// Cartesian components

// ----------------------------------------------------------------------------
//                     Functions To Read The Asymmetry (eta)
// ----------------------------------------------------------------------------

/* These functions are used to obtain the tensor asymmetry value, eta, from
   a parameter set. Note that the parameter read, "Pbase" is generic so that
   derived classes may supply their own parameter acceptable parameter names.  
   The default parameter name base of this class is "IR2", so IR2eta is one
   possible parameter name. The eta values must be unitless and restricted 
   between the range [0, 1]. Usually we're happy if no eta value is specified,
   it's just set to zero then, but this is controlled by the warnings flag.

   There are a few ways to set the asymmetry value. The easiest way is directly
   with use of a double value, parameter name "IR2eta". Since this is a 
   spherical tensor component, we must also allow for a Cartesian tensor 
   specification of eta. This is tied up in the PAS values {Axx, Ayy, Azz}.
   Lastly, if the user chooses to input a Cartesian tensor that is NOT in the
   PAS, i.e. { Axx, Axy, Axz, Ayy, Ayz, Azz }, the asymmetry and orientation
   are both tied into the same set of values so they must both be determined
   simultaneously.  Worse, this often demands an iterative proceedure to
   determine the Euler angles.                                               */

bool IntRank2A::getAeta(const ParameterSet& pset, const string& A,
                            double& Aeta, int idxI, int idxS, bool warn) const
  {
  Aeta = 0.0;					// Begin with no asymmetry
  string Nidx("");                              // Parameter name suffix
  if(idxI >= 0)                                 // Suffix only if idxI > -1
    {                                           //   Add idxI index to name
    Nidx  = string("(") + Gdec(idxI);            //   as (#
    if(idxS > 0)                                //   If idxS > 0 add it to
      Nidx += string(",")+Gdec(idxS);            //   to suffix as ,# (#,#
    Nidx += string(")");                        //   Finish suffix (#) or (#,#)
    }
  string pname = A + "eta" + Nidx;		// String for parameter name
  string pstate;                                // Temp string for parsing
  ParameterSet::const_iterator item;         // A pix into parameter list
  item = pset.seek(pname);                      // Seek parameter in pset
  if(item == pset.end())                        // If not found, warn if
    {                                           // desired, return fail
    if(warn) IR2Aerror(111,pname,1);
    return false;
    }
  (*item).parse(pname,Aeta,pstate);		// If found, parse/set Aeta
  return true;                                  // Return we found it OK
  }

// ----------------------------------------------------------------------------
//         Functions To Read Interaction Orientation (Euler Angles)
// ----------------------------------------------------------------------------

/* Each interaction stores one set of Euler Angles {alpha, beta, gamma} that
   describe its orientation with respect to its principal axis system (PAS).
   Typically these angles are those that relate the PAS to the laboratory
   frame, but they may easily be used to relate the PAS to some other frame
   (molecular, diffusion, ......). However, GAMMA contains composite
   rotations well suited for employing any number of reference frames. We are
   usually quite happy if no orientation is specified and we then just set the
   tensor in its PAS.  However, this is controlled by the warnings flag.

   All parameters will be prefixed with Pbase so that derived classes may use
   their own naming convention.  The default Pbase is "IR2" for this class.
   The easiest way to specify the orientation (after simply not setting it)
   is through use of a set of Euler angles { alpha, beta, gamma }. This may be
   done either with individual angles (e.g. IR2alpha, IR2beta, IR2gamma) or
   with a single Euler angle set (e.g. IR2EAngles). Another useful way, albeit
   more complicated, is to specify the 5 Cartesian tensor components directly,
   { Axx, Axy, Axz, Ayy, Ayz, Azz }. In that case the orientation and 
   asymmetry are both tied into the same set of values so they must both be
   determined at the same time.  Worse, this often demands an iterative
   proceedure to determine the Euler angles.                                 */

bool IntRank2A::getOrientation(const ParameterSet& pset, const string& Pbase,
                              EAngles& EA, int idxI, int idxS, bool warn) const
  {
  string Nidx = "";                             // Name addition per index
  if(idxI >= 0)                                 // Suffix only if idxI > -1
    {                                           //   Add idxI index to name
    Nidx  = string("(") + Gdec(idxI);            //   as (#
    if(idxS > 0)                                //   If idxS > 0 add it to
       Nidx += string(",")+Gdec(idxS);           //   to suffix as ,# (#,#
    Nidx += string(")");                        //   Finish suffix (#) or (#,#)
    }
  double pdatad;                                // Double for angle value
  string pstate;                                // Temp string for parsing

//         Try And Read A Set Of Euler Angles { Alpha, Beta, Gamma }
//                    ( These Will Externally Be in Degrees )

  string pname = Pbase + "EAngles" + Nidx;	// Euler angles parameter name
  coord EAtmp;					// Temp coordinate for read
  ParameterSet::const_iterator item;         // Pix in parameter list
  item  = pset.seek(pname);			// Get the pix for parameter
  if(item != pset.end())			// If the parameter is found
    {						// then we glean the EAngles
    EAtmp = coord(*item);			//   Read as a coordnate (deg)
    EA.alpha(EAtmp.x() * DEG2RAD);		//   Set Euler angles, we
    EA.beta(EAtmp.y()  * DEG2RAD);		//   convert degrees to 
    EA.gamma(EAtmp.z() * DEG2RAD); 		//   radians for storage
    return true;
    }

//         Try And Read 3 Individual Euler Angles { Alpha, Beta, Gamma }
//                    ( These Will Externally Be in Degrees )

  double alpha=0, beta=0, gamma=0;              // Temp Euler angles
  pname = Pbase + "Alpha" + Nidx; 		// Parameter name for alpha
  if(pset.getDouble(pname, pdatad))             // If it's been found
    {                                           // we look for beta and gamma
    alpha = pdatad*DEG2RAD;                     //   Keep the alpha value
    pname = Pbase + "Beta" + Nidx;              //   Parameter name for beta
    if(pset.getDouble(pname, pdatad))           //   Look for beta value
      beta = pdatad*DEG2RAD;                    //   Keep beta value if found
    pname = Pbase + "Gamma" + Nidx;             //   Parameter name for gamma
    if(pset.getDouble(pname, pdatad))           //   Look for gamma value
      gamma = pdatad*DEG2RAD;                   //   Keep gamma value if found
    EA = EAngles(alpha,beta,gamma);             //   Set Euler angles
    return true;
    }

//               Try And Read Euler Angles { Theta, Phi, --- }
//                    ( These Will Externally Be in Degrees )

  pname = Pbase + "Theta" + Nidx;               // Parameter name for theta
  if(pset.getDouble(pname, pdatad))             // If it's been found
    {                                           // we look for phi also
    beta = pdatad*DEG2RAD;			//   Keep theta value if found
    pname = Pbase + "Phi" + Nidx;               //   Parameter name for phi
    if(pset.getDouble(pname, pdatad))           //   Look for phi value
      alpha = pdatad*DEG2RAD;			//   Keep phi value if found
    EA = EAngles(alpha,beta,gamma);             //   Set Euler angles
    return true;
    }

//         Try And Read 3 Individual Euler Angles { alpha, beta, gamma }
//                    ( These Will Externally Be in Degrees )

  pname = Pbase + "alpha" + Nidx; 		// Parameter name for alpha
  if(pset.getDouble(pname, pdatad))             // If it's been found
    {                                           // we look for beta and gamma
    alpha = pdatad*DEG2RAD;                     //   Keep the alpha value
    pname = Pbase + "beta" + Nidx;              //   Parameter name for beta
    if(pset.getDouble(pname, pdatad))           //   Look for beta value
      beta = pdatad*DEG2RAD;                    //   Keep beta value if found
    pname = Pbase + "gamma" + Nidx;             //   Parameter name for gamma
    if(pset.getDouble(pname, pdatad))           //   Look for gamma value
      gamma = pdatad*DEG2RAD;                   //   Keep gamma value if found
    EA = EAngles(alpha,beta,gamma);             //   Set Euler angles
    return true;
    }

//                      Issue Warnings If No Angle Found
//                         (Only If Warning Flag Set)

  if(warn)					// Warn if desired & cannot
    {						// find any orientation angles
    IR2Aerror(51,1);				//    Cant find orient. angles
    IR2Aerror(52,1);				//    Unable to set orientation
    }
  return false;                                 // Could not find orientation
  }

// ----------------------------------------------------------------------------
//                 Functions To Read The Cartesian Components
// ----------------------------------------------------------------------------

/* Users may find it convenient to define a spatial tensor with Cartesian
   components.  This can either be the three PAS components { Axx, Ayy, Azz }
   or the five components that represent the tensor in some arbitrary 
   orientation { Axx, Axy, Axz, Ayy, Ayz, Azz }. If the former is used (i.e.
   Axy, Axz, and Ayz are not specified) the we can directly set the asymmetry
   fro the three values { Axx, Ayy, Azz } and look elsewhere for any 
   specification of the tensor orientation. If the latter is used then we have
   a more complicated situation on our hands because the orientation and the
   asymmetry are both tied into the same set of values.  In that case they must
   both be determined at the same time from Axx, Axy, Axz, Ayy, Ayz, Azz }. 
   In a general casel this will demand an iterative proceedure to determine
   the Euler angles.

   A few more items worthy of mention. First, all parameters will be prefixed
   with Pbase so that derived classes may use their own naming convention.
   The default Pbase here is "IR2A" for this class which means typical
   parameter names would be IR2Axx, IR2Axy, etc. Second, the functions do not
   discriminate between reducible and irreducible components. They do not 
   care if Axx + Ayy + Azz = 0 or not, thus allowing for use in derived
   classes that may have an isotropic component.                             */

int IntRank2A::getAxAyAz(const ParameterSet& pset, const string& A,
                              coord& Axyz, int idxI, int idxS, bool warn) const
  {
  string Nidx("");				// Parameter name suffix
  if(idxI >= 0)					// Suffix only if idxI > -1
    {						//   Add idxI index to name
    Nidx  = string("(") + Gdec(idxI);	        //   as (#
    if(idxS > 0) 				//   If idxS > 0 add it to
      Nidx += string(",")+Gdec(idxS); 		//   to suffix as ,# (#,#
    Nidx += string(")");			//   Finish suffix (#) or (#,#)
    }

  int cnt = 0;					// Count the values found
  double x=0,y=0,z=0;                           // Temp storage of cmpts
  string pstate;				// Temp string for parameter
  ParameterSet::const_iterator item;         // A pix into parameter list
 
  string pn = A + string("xx") + Nidx;		// Parameter name for Axx
  item = pset.seek(pn);				// Seek Axx parameter in pset
  if(item == pset.end()) 			// If not found, warn if
    { if(warn) IR2Aerror(111,pn,1); } 		// desired
  else						// If found, get Axx value
    { (*item).parse(pn,x,pstate); cnt++; } 	// set it to x and count it
 
  pn = A + string("yy") + Nidx;			// Parameter name for Ayy
  item = pset.seek(pn);				// Seek Ayy parameter in pset
  if(item == pset.end()) 			// If not found, warn if
    { if(warn) IR2Aerror(111,pn,1); } 		// desired
  else						// If found, get Ayy value
    { (*item).parse(pn,y,pstate); cnt++; } 	// set it to y and count it
 
  pn = A + string("zz") + Nidx;			// Parameter name for Azz
  item = pset.seek(pn);				// Seek Azz parameter in pset
  if(item == pset.end()) 			// If not found, warn if
    { if(warn) IR2Aerror(111,pn,1); }		// desired
  else						// If found, get Azz value
    { (*item).parse(pn,z,pstate); cnt++; } 	// set it to z and count it

  if(!cnt)					// If we have not found any
    {						// Auu values, try for all 3
    pn = A + "x" + A + "y" + A + "z" + Nidx;	// using coord AxAyAz
    item = pset.seek(pn);			//   Seek AxAyAz in pset
    if(item == pset.end()) 			//   If not found, warn if
      { if(warn) IR2Aerror(111,pn,1); } 	//   desired
    else					//   If found, get 3 Auu
      {						//   directly as coordinate 
      Axyz = coord(*item);			//   Return we found all 3
      return 3;
      }
    }
  Axyz = coord(x,y,z);				// Put results into Axyz
  return cnt;					// Set cnt of {Axx,Ayy,Azz}
  }                                             //   (these are unsorted!)

int IntRank2A::getAoffdiag(const ParameterSet& pset, const string& A,
                              coord& Aod, int idxI, int idxS, bool warn) const
  {
  string Nidx("");				// Parameter name suffix
  if(idxI >= 0)					// Suffix only if idxI > -1
    {						//   Add idxI index to name
    Nidx  = string("(") + Gdec(idxI);	        //   as (#
    if(idxS > 0) 				//   If idxS > 0 add it to
      Nidx += string(",")+Gdec(idxS); 		//   to suffix as ,# (#,#
    Nidx += string(")");			//   Finish suffix (#) or (#,#)
    }

  int cnt = 0;					// Count values found
  double xy=0,xz=0,yz=0;			// Temp storage of cmpts
  string pstate;				// Temp string for parameter
  ParameterSet::const_iterator item;         // A pix into parameter list
 
  string pn = A + string("xy") + Nidx;		// Parameter name for Axy
  item = pset.seek(pn);				// Seek Axy parameter in pset
  if(item == pset.end()) 			// If not found, warn if
    { if(warn) IR2Aerror(111,pn,1); }		// desired
  else						// If found, get Axy value
    { (*item).parse(pn,xy,pstate); cnt++; } 	// set it to xy and count it
 
  pn = A + string("xz") + Nidx;			// Parameter name for Axz
  item = pset.seek(pn);				// Seek Axz parameter in pset
  if(item == pset.end()) 			// If not found, warn if
    { if(warn) IR2Aerror(111,pn,1); }		// desired
  else						// If found, get Axz value
    { (*item).parse(pn,xz,pstate); cnt++; } 	// set it to xz and count it
 
  pn = A + string("yz") + Nidx;			// Parameter name for Ayz
  item = pset.seek(pn);				// Seek Ayz parameter in pset
  if(item == pset.end()) 			// If not found, warn if
    { if(warn) IR2Aerror(111,pn,1); }		// desired
  else						// If found, get Ayz value
    { (*item).parse(pn,yz,pstate); cnt++; } 	// set it to yz and count it

  Aod = coord(xy,xz,yz);			// Put results into Aod
  return cnt;					// Set cnt of {Axy,Axz,Ayz}
  }

// ----------------------------------------------------------------------------
//                 Functions To Read Int Two Spin Coordinates
// ----------------------------------------------------------------------------

/* Spin pair interactions, such as dipolar and hyperfine interactions, may have
   their orientation (and perhaps their strength) set by a pair of "spin" 
   coordinates. Here we read in two coordinates, one for each spin index and
   return them. Failure will occur if either spin coordinate is missing.     */

bool IntRank2A::getCoords(const ParameterSet& pset,
                   coord& ptI, coord& ptS, int idxI, int idxS, bool warn) const
  {
  int p1 = ptI.read(pset, idxI, warn);		// Try to get spin I coordinate
  int p2 = ptS.read(pset, idxS, warn);		// Try to get spin S coordinate
  if(p1 && p2) return true;			// OK if both have been found
  if(warn)					// If not both found we will
    {						// issue warnings if desired
    if(!p1) IR2Aerror(113, Gdec(idxI), 1);
    if(!p2) IR2Aerror(113, Gdec(idxS), 1);
    }
  return false;                                 // Could not find orientation
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A             RANK 2 INTERACTION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

IntRank2A::IntRank2A()                      { ETA=0.0; }
IntRank2A::IntRank2A(const IntRank2A &IR2A) { ETA=IR2A.ETA; _EAs=IR2A._EAs; }

// ----------------------------------------------------------------------------
//         Direct Constructors Using Cartesian Spatial Tensor Components
// ----------------------------------------------------------------------------
  
        // Input                IR2A	: Rank 2 spatial tensor (this)
	//			AxAyAz  : Cartesian PAS values
	// Output		none    : Rank 2 spatial tensor constructed
	// Note				: Theta PAS z down, phi PAS x over
	// Note				: Insist |Azz| >= |Ayy| >= |Axx|
	// Note				: Dont care if Auv are irreducible here

IntRank2A::IntRank2A(const coord& AxAyAz, const EAngles& EA)
  {
  double x    = AxAyAz.x();			// Input Axx
  double y    = AxAyAz.y();			// Input Ayy
  double z    = AxAyAz.z();			// Input Azz
  SortAxAyAz(x, y, z);				// Keep |Azz| >= |Ayy| >= |Axx|
  ETA   = (x-y)/(-x-y);				// Set asymmetry  
  _EAs = EA;					// Set orientation angles
  }
 
// ----------------------------------------------------------------------------
//         Direct Constructors Using Spherical Spatial Tensor Components
// ----------------------------------------------------------------------------
 
IntRank2A::IntRank2A(double eta, const EAngles& EA)
  {
  if(!CheckEta(eta)) IR2Aerror(2);		// Insure eta range [0,1]
  ETA = eta;					// Set the rank 2 asymmetry
  _EAs = EA;					// Set orientation angles
  }

// ----------------------------------------------------------------------------
//                       Constructors Using Parameter Sets
// ----------------------------------------------------------------------------

/* In this case the constructon is performed entirely from parameters found in
   the specified parameter set.  The indices allows users to read in multiple
   spatial tensors by setting either (#) or (#,#) suffix on parameter names.
   If the first # is default (idxI=-1) not suffix will be used.  If the second
   # is default (idxS=-1) then only a single (#) suffix is used, i.e. the #
   is an interaction index. When both #'s are used and the suffix is (#,#) 
   the numbers are loosely associated with spin indices.

   The parameter set is used to obtain only 2 pieces of information. The
   asymmetry value (eta) and the spatial tensor orientation (alpha,beta,gamma).
   Allowances are made so that the asymmetry can be set in either spherical
   or Cartesian terms and so that the orientation can be set using either
   individual angles or as on set of Euler angles.

           Input                IR2A    : Rank 2 spatial tensor
                                pset    : Parameter set
                                idxI    : Interaction index (default -1->none)
                                idxS    : Interaction index (default -1->none)
           Output               none    : Rank 2 spatial tensor constructed
                                          from parameters in pset             */

IntRank2A::IntRank2A(ParameterSet& pset, int idxI, int idxS, int warn)
  {
  string pnb("IR2");					// Base name for eta
  bool TF = getAX(pset,pnb,ETA,_EAs,idxI,idxS,1);	// Try & read tensor
  if(!TF)						// If unsuccessful,
    {							// make sure we are zero
    ETA  = 0.0;						// & take some action
    _EAs = EAzero;
    if(warn)
      {
      if(warn > 1)  IR2Afatal(21);			// action ala warn flag
      else          IR2Aerror(53);                        //   Setting EAs to 0
      }
    }
  
// sosi -this is old code that did not read cartesian spatial tensors
/*
  if(!getAeta(pset,pnb,ETA,idxI,idxS,warn?true:false)	// Try & read eta value
  || !CheckEta(ETA,warn?true:false))			// it must be in [0,1]
    {							// If unable take some
    if(warn > 2)  IR2Afatal(2);				// action ala warn flag
    else if(warn) IR2Aerror(11);    			//   Setting eta to 0
    ETA = 0.0;						//   Set default 0
    }
  if(!getOrientation(pset,pnb,_EAs,idxI,idxS,warn?true:false))	// Try & read orient.
    {
    if(warn > 2)  IR2Afatal(2);				// action ala warn flag
    else if(warn) IR2Aerror(53);    			//   Setting EAs to 0
    _EAs = EAzero;
    }
*/
  }
   
// ----------------------------------------------------------------------------
//                         Assignment and Destruction
// ----------------------------------------------------------------------------

void IntRank2A::operator= (const IntRank2A &IR2A) {ETA=IR2A.ETA; _EAs=IR2A._EAs;}
     IntRank2A::~IntRank2A()                      { }

// ____________________________________________________________________________
// B                RANK 2 ANISOTROPY & ASYMMETRY FUNCTIONS
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


double IntRank2A::delzz() { return     RT5O6PI; }
double IntRank2A::delA()  { return 1.5*RT5O6PI; }

//-----------------------------------------------------------------------------
//                             Asymmetry Access
//-----------------------------------------------------------------------------

/* These functions allow the user to get or set the spatial tensor asymmetry.
   In GAMMA we use the convention that eta spans [0, 1] and it is defined to be
   (Axx-Ayy)/Azz where |Azz| >= |Axx| >= |Ayy| in the treatment herein.

	 Input		IR2A	: Rank 2 spatial tensor
			Eta	: Rank 2 tensor asymmetry
	 Output		eta 	: Return asymmetry of rank 2
				  spatial tensor
                    OR  void    : Set the asymmetry of the spatial tensor    */

double IntRank2A::eta( ) const { return ETA; }
void   IntRank2A::eta(double Eta)
  {
  if(!CheckEta(Eta)) IR2Afatal(10); 		// Insure range is [0, 1]
  ETA = Eta;					// Set the eta value 
  }

// ____________________________________________________________________________
// C            RANK 2 TENSOR SPHERICAL COMPONENT ACCESS FUNCTIONS
// ____________________________________________________________________________
 
/* These functions allow one to access irreducible spherical elements at 
   specific orientations without rotating the entire tensor. In these functions
   theta is the angle down from the PAS z axis and phi the angle over from PAS
   x axis. Note that these all use GAMMA scaling which sets A2m to be rank 2 
   spherical harmonics at all orientations if there is no asymmetry.

           Input                IR2A	: Rank 2 spatial tensor 
                                theta   : Orientation angle (degrees)
                                phi     : Orientation angle (degrees)
           Return               A2m     : Rank 2 spherical element (Hz)
                                          for orientation {phi, theta}  
                                          from the tensor PAS                */

// ----------------------------------------------------------------------------
//                         m = 0 Component : A2,0 
// ----------------------------------------------------------------------------

/*                        1/2
                   [  5  ]          2                     2 
 A   (theta,phi) = |-----| * [ 3*cos (theta) - 1 + eta*sin (theta)cos(2*phi) ]
  2,0              [16*PI]
                                                  1/2                 1/2
         |                                 [  5  ]                 [3]
 A  (T,P)|      = Y   (T, P)   A   (0,0) = |-----|     A   (T,P) = |-| A  (T,P)
  2,0    |         2,0          2,0        [ 4*PI]      2,0        [2]  zz
          eta=0                                                              */

complex IntRank2A::A20PAS()       { return RT5O4PI; }
complex IntRank2A::A20()    const { return A20(_EAs); }

complex IntRank2A::A20(double A, double B, double G) const
  {
  double Cbeta   = cos(B);			// cos(beta)
  double Sbeta   = sin(B);			// sin(beta)
  double C2gamma = cos(2.0*G);	 		// cos(2*gamma)
  return sqrt(5./(16.*PI))*(3.*Cbeta*Cbeta - 1. + ETA*Sbeta*Sbeta*C2gamma);
  }

complex IntRank2A::A20(const EAngles& EA) const
  {
  double B=EA.beta();				// Angle beta  (radians)
  double G=EA.gamma();				// Angle gamma (radians)
  double Cbeta   = cos(B);			// cos(beta)
  double Sbeta   = sin(B);			// sin(beta)
  double C2gamma = cos(2.0*G);	 		// cos(2*gamma)
  return sqrt(5./(16.*PI))*(3.*Cbeta*Cbeta - 1. + ETA*Sbeta*Sbeta*C2gamma);
  }

// ----------------------------------------------------------------------------
//                 m = 1 Component : A2,1  [ A2,1 = -(A2,-1)* ]
// ----------------------------------------------------------------------------

/*                         1/2
                    [  5  ]              [
 A   (theta,phi)  = |-----| * sin(theta)*| 3*cos(theta)
  2,1               [24*PI]              [
                                                                            ]
                                - eta*[cos(theta)*cos(2*phi) - i*sin(2*phi)]|
                                                                            ]
          |                                                   
 A   (T,P)|    = Y   (T,P)   A   (0,0) = 0   A   (T,P) = -A  (T,P) - iA  (T,P)
  2,1     |       2,1         2,1             2,1          xz          yz
           eta=0                                                              */

complex IntRank2A::A21PAS()       { return 0.0; }
complex IntRank2A::A21()    const { return A21(_EAs);; }

complex IntRank2A::A21(double A, double B, double G) const
  {
  double  twoG    = 2.0*G;			// 2*gamma
  complex Ealpha  = exp(complexi*A);		// exp(i*alpha)
  double  Cbeta   = cos(B);			// cos(theta)
  double  Sbeta   = sin(B);			// sin(theta)
  double  S2beta  = sin(2.0*B);			// sin(theta)
  double  C2gamma = cos(twoG); 		  	// cos(2*gamma)
  double  S2gamma = sin(twoG);   		// sin(2*gamma)
  complex term    = complex(Cbeta*C2gamma,S2gamma);
  return -RT5O96PI*Ealpha*(3.*S2beta-(2.*ETA*Sbeta)*term);
  }

complex IntRank2A::A21(const EAngles& EA) const
  {
  double  A       = EA.alpha();			// Angle beta  (radians)
  double  B       = EA.beta();			// Angle beta  (radians)
  double  G       = EA.gamma();			// Angle gamma (radians)
  double  twoG    = 2.0*G;			// 2*gamma
  complex Ealpha  = exp(complexi*A);		// exp(i*alpha)
  double  Cbeta   = cos(B);			// cos(theta)
  double  Sbeta   = sin(B);			// sin(theta)
  double  S2beta  = sin(2.0*B);			// sin(theta)
  double  C2gamma = cos(twoG); 		  	// cos(2*gamma)
  double  S2gamma = sin(twoG);   		// sin(2*gamma)
  complex term    = complex(Cbeta*C2gamma,S2gamma);
  return -RT5O96PI*Ealpha*(3.*S2beta-(2.*ETA*Sbeta)*term);
  }

// ----------------------------------------------------------------------------
//               m = -1 Component : A2,-1  [ A2,-1 = (-A2,1)* ]
// ----------------------------------------------------------------------------

/*                          1/2
                     [  5  ]              [
 A    (theta,phi)  = |-----| * sin(theta)*| 3*cos(theta)
  2,-1               [24*PI]              [
                                                                            ]
                                - eta*[cos(theta)*cos(2*phi) + i*sin(2*phi)]|
                                                                            ]
          |                                                   
A    (T,P)|    = Y    (T,P)  A    (0,0)=0  A    (T,P) = A  (T,P) - iA  (T,P)
 2,-1     |       2,-1        2,-1          2,-1         xz          yz
           eta=0                                                             */
//  complex term    = complex(C2gamma*(1+Cbeta*Cbeta),2*S2gamma*Cbeta);

complex IntRank2A::A2m1PAS()       { return 0.0; }
complex IntRank2A::A2m1()    const { return A2m1(_EAs); }

complex IntRank2A::A2m1(double A, double B, double G) const
  {
  double  twoG    = 2.0*G;			// 2*gamma
  complex Ealpha  = exp(-complexi*A);		// exp(-i*alpha)
  double  Cbeta   = cos(B);			// cos(theta)
  double  Sbeta   = sin(B);			// sin(theta)
  double  S2beta  = sin(2.0*B);			// sin(theta)
  double  C2gamma = cos(twoG); 		  	// cos(2*gamma)
  double  S2gamma = sin(twoG);   		// sin(2*gamma)
  complex term    = complex(Cbeta*C2gamma,-S2gamma);
  return  RT5O96PI*Ealpha*(3.*S2beta-(2.*ETA*Sbeta)*term);
  }

complex IntRank2A::A2m1(const EAngles& EA) const
  {
  double  A       = EA.alpha();			// Angle beta  (radians)
  double  B       = EA.beta();			// Angle beta  (radians)
  double  G       = EA.gamma();			// Angle gamma (radians)
  double  twoG    = 2.0*G;			// 2*gamma
  complex Ealpha  = exp(-complexi*A);		// exp(-i*alpha)
  double  Cbeta   = cos(B);			// cos(theta)
  double  Sbeta   = sin(B);			// sin(theta)
  double  S2beta  = sin(2.0*B);			// sin(theta)
  double  C2gamma = cos(twoG); 		  	// cos(2*gamma)
  double  S2gamma = sin(twoG);   		// sin(2*gamma)
  complex term    = complex(Cbeta*C2gamma,-S2gamma);
  return  RT5O96PI*Ealpha*(3.*S2beta-(2.*ETA*Sbeta)*term);
  }

// ----------------------------------------------------------------------------
//               m = 2 Component : A2,2  [ A2,2 = (A2,-2)* ]
// ----------------------------------------------------------------------------
 
/*                        1/2
                   [  5  ]   [      2                  2 
 A   (theta,phi) = |-----|   [ 3*sin (theta) + eta*[cos (theta)+1]*cos(2*phi)
  2,2              [96*PI]   [
                                                                              ]
                                               - 2i*eta*cos(theta)*sin(2*phi) ]
                                                                              ]
                                                             1/2
                |                                     [  5  ]                   
       A   (T,P)|     = Y   (T,P)         A   (0,0) = |-----| * eta
        2,2     |        2,2               2,2        [24*PI]

                         1 [                     ]
             A   (T,P) = - | A  (T,P) - A  (T,P) | + i A  (T,P)
              2,2        2 [  xx         yy      ]      xy                   */

complex IntRank2A::A22PAS()       { return RT5O24PI;  }
complex IntRank2A::A22()    const { return A22(_EAs); }

complex IntRank2A::A22(double A, double B, double G) const
  {
  double  twoG    = 2.0*G;			// 2*gamma
  complex E2alpha = exp(2*complexi*A);		// exp(i*2*alpha)
  double  Cbeta   = cos(B);			// cos(theta)
  double  Sbeta   = sin(B);			// sin(theta)
  double  C2gamma = cos(twoG); 		  	// cos(2*gamma)
  double  S2gamma = sin(twoG);   		// sin(2*gamma)
  complex term    = complex(C2gamma*(1+Cbeta*Cbeta),2*S2gamma*Cbeta);
  return  RT5O96PI*E2alpha*(3.*Sbeta*Sbeta+ETA*term);
  }

complex IntRank2A::A22(const EAngles& EA) const
  {
  double  A       = EA.alpha();			// Angle beta  (radians)
  double  B       = EA.beta();			// Angle beta  (radians)
  double  G       = EA.gamma();			// Angle gamma (radians)
  double  twoG    = 2.0*G;			// 2*gamma
  complex E2alpha = exp(2*complexi*A);		// exp(i*2*alpha)
  double  Cbeta   = cos(B);			// cos(theta)
  double  Sbeta   = sin(B);			// sin(theta)
  double  C2gamma = cos(twoG); 		  	// cos(2*gamma)
  double  S2gamma = sin(twoG);   		// sin(2*gamma)
  complex term    = complex(C2gamma*(1+Cbeta*Cbeta),2*S2gamma*Cbeta);
  return  RT5O96PI*E2alpha*(3.*Sbeta*Sbeta+ETA*term);
  }


// ----------------------------------------------------------------------------
//              m = -2 Component : A2,-2  [ A2,-2 = (A2,2)* ]
// ----------------------------------------------------------------------------
 
/*                        1/2
                   [  5  ]   [      2                  2 
 A   (theta,phi) = |-----|   [ 3*sin (theta) + eta*[cos (theta)+1]*cos(2*phi)
  2,-2             [96*PI]   [
                                                                              ]
                                               + 2i*eta*cos(theta)*sin(2*phi) ]
                                                                              ]
                                                              1/2
                |                                      [  5  ]
      A    (T,P)|     = Y    (T,P)        A    (0,0) = |-----| * eta
       2,-2     |        2,-2              2,-2        [24*PI]

                          1 [                     ]
             A    (T,P) = - | A  (T,P) - A  (T,P) | - i A  (T,P)
              2,-2        2 [  xx         yy      ]      xy                   */

complex IntRank2A::A2m2PAS()       { return RT5O24PI;   }
complex IntRank2A::A2m2()    const { return A2m2(_EAs); }

complex IntRank2A::A2m2(double A, double B, double G) const
  {
  double  twoG    = 2.0*G;			// 2*gamma
  complex E2alpha = exp(-2*complexi*A);		// exp(-i*2*alpha)
  double  Cbeta   = cos(B);			// cos(theta)
  double  Sbeta   = sin(B);			// sin(theta)
  double  C2gamma = cos(twoG); 		  	// cos(2*gamma)
  double  S2gamma = sin(twoG);   		// sin(2*gamma)
  complex term    = complex(C2gamma*(1+Cbeta*Cbeta),-2*S2gamma*Cbeta);
  return  RT5O96PI*E2alpha*(3.*Sbeta*Sbeta+ETA*term);
  }

complex IntRank2A::A2m2(const EAngles& EA) const
  {
  double  A       = EA.alpha();			// Angle beta  (radians)
  double  B       = EA.beta();			// Angle beta  (radians)
  double  G       = EA.gamma();			// Angle gamma (radians)
  double  twoG    = 2.0*G;			// 2*gamma
  complex E2alpha = exp(-2*complexi*A);		// exp(-i*2*alpha)
  double  Cbeta   = cos(B);			// cos(theta)
  double  Sbeta   = sin(B);			// sin(theta)
  double  C2gamma = cos(twoG); 		  	// cos(2*gamma)
  double  S2gamma = sin(twoG);   		// sin(2*gamma)
  complex term    = complex(C2gamma*(1+Cbeta*Cbeta),-2*S2gamma*Cbeta);
  return  RT5O96PI*E2alpha*(3.*Sbeta*Sbeta+ETA*term);
  }

// ----------------------------------------------------------------------------
//                    m Component : A2,m  m = [-2,2]
// ----------------------------------------------------------------------------

complex IntRank2A::A2m(int m) const
  {
  switch(m)
    {
    case  0: return A20();  break;
    case  1: return A21();  break;
    case  2: return A22();  break;
    case -1: return A2m1(); break;
    case -2: return A2m2(); break;
    }
  return complex0;
  }

complex IntRank2A::A2m(int m, double alpha, double beta, double gamma) const
  {
  switch(m)
    {
    case  0: return A20(alpha,  beta, gamma); break;
    case  1: return A21(alpha,  beta, gamma); break;
    case  2: return A22(alpha,  beta, gamma); break;
    case -1: return A2m1(alpha, beta, gamma); break;
    case -2: return A2m2(alpha, beta, gamma); break;
    }
  return complex0;
  }

complex IntRank2A::A2m(int m, const EAngles& EA) const
  {
  switch(m)
    {
    case  0: return A20(EA);  break;
    case  1: return A21(EA);  break;
    case  2: return A22(EA);  break;
    case -1: return A2m1(EA); break;
    case -2: return A2m2(EA); break;
    }
  return complex0;
  }

// ----------------------------------------------------------------------------
//                         All 5 Spherical Components
// ----------------------------------------------------------------------------

/* These function generate and return (the equivalent of) all 5 irreducible
   rank 2 spatial tensor spherical components.                               */

IR2ASph IntRank2A::SphCmpPAS() const
  {
  IR2ASph SC;				// Spherical components
  SC._A20 = RT5O4PI; 			// A20(PAS) = sqrt[5/(4*PI)]
  SC._A21 = 0;				// A21(PAS) = 0.0
  SC._A22 = RT5O24PI*ETA;		// A22(PAS) = sqrt[5/(24*PI)] * eta
  return SC;
  }

IR2ASph IntRank2A::SphCmp() const 
  { return SphCmp(_EAs); }

IR2ASph IntRank2A::SphCmp(double A, double B, double G) const
  { EAngles EA(A,B,G); return SphCmp(EA); }
  
IR2ASph IntRank2A::SphCmp(const EAngles& EA) const
  {
  IR2ASph SC;					// Spherical components

  double A = EA.alpha();			// Angle alpha (radians)
  double B = EA.beta();				// Angle beta  (radians)
  double G = EA.gamma();			// Angle gamma (radians)

  complex Ealpha  = exp(complexi*A);		// exp(i*alpha)
  complex E2alpha = exp(2*complexi*A);		// exp(i*2*alpha)

  double  Cbeta   = cos(B);			// cos(beta)
  double  Sbeta   = sin(B);			// sin(beta)
  double  Cbetasq = Cbeta*Cbeta;		// cos(beta)*cos(beta)
  double  Sbetasq = Sbeta*Sbeta;		// sin(beta)*sin(beta)
  double  S2beta  = sin(2.0*B);			// sin(2*beta)

  double  twoG    = 2.0*G;			// 2*gamma
  double  C2gamma = cos(twoG);	 		// cos(2*gamma)
  double  S2gamma = sin(twoG);   		// sin(2*gamma)

  double RT5O16PI = sqrt(5./(16.*PI));
  complex z1(Cbeta*C2gamma,S2gamma);
  complex z2(C2gamma*(1+Cbetasq),2*S2gamma*Cbeta);

  SC._A20 =  RT5O16PI*(3.*Cbetasq - 1. + ETA*Sbetasq*C2gamma);
  SC._A21 = -RT5O96PI*Ealpha*(3.*S2beta-(2.*ETA*Sbeta)*z1);
  SC._A22 =  RT5O96PI*E2alpha*(3.*Sbetasq + ETA*z2);

  return SC;
  }

/* This function returns all 5 spherical PAS components but it uses
   the angular momentum indexing scheme m = {-2,-1,0,1,2}.                   */

complex IntRank2A::AcompPAS(int comp) const
  {
  switch(comp)
    {
    case  0: return RT5O4PI;      break;
    case  1: 
    case -1: return 0;            break;
    case  2: 
    case -2: return RT5O24PI*ETA; break;
    default:
      {
      IR2Aerror(50);
      IR2Afatal(0);
      }
    }
  return complex0; 
  }

complex IntRank2A::Acomp(int comp) const
  {
  switch(comp)
    {
    case  0: return A20(_EAs);  break;
    case  1: return A21(_EAs);  break;
    case -1: return A2m1(_EAs); break;
    case  2: return A22(_EAs);  break;
    case -2: return A2m2(_EAs); break;
    default:
      {
      IR2Aerror(50);
      IR2Afatal(0);
      }
    }
  return complex0;
  }


// ____________________________________________________________________________
// D            RANK 2 TENSOR CARTESIAN COMPONENT ACCESS FUNCTIONS
// ____________________________________________________________________________
 
/* These functions allow one to access the Cartesian elements at a specific 
   orientation without rotating the entire tensor.  In these functions theta is
   the angle down from the PAS z axis and phi the angle over from the PAS x 
   axis. Note that these all use GAMMA scaling which is based on the spherical
   components, A2m, being normalized to rank 2 spherical harmonics at all
   orientations if there is no asymmetry.

           Input                IR2A	: Rank 2 spatial tensor 
                                theta   : Orientation angle (degrees)
                                phi     : Orientation angle (degrees)
           Return               Auv     : Rank 2 cartesian component (Hz)
                                          for orientation {phi, theta}       */
 
// ----------------------------------------------------------------------------
//                           xx Component : Axx
// ----------------------------------------------------------------------------

/*                        1/2
                   [  5  ]          2                     2 
  A  (theta,phi) = |-----| * [ 3*sin (theta) - 1 + eta*cos (theta)cos(2*phi) ]
   xx              [24*PI]


                    1 [                        ]   [1] 1/2
         A  (T,P) = - | A   (T,P) + A    (T,P) | - |_|     A   (T,P)
          xx        2 [  2,2         2,-2      ]   [6]      2,0              */

double IntRank2A::AxxPAS() const { return -RT5O24PI*(1.0-ETA); }

double IntRank2A::Axx()    const
  { return Axx(_EAs.alpha(), _EAs.beta(), _EAs.gamma()); }

double IntRank2A::Axx(const EAngles& EA) const
  { return Axx(EA.alpha(), EA.beta(), EA.gamma()); }

double IntRank2A::Axx(double alpha, double beta, double gamma) const
  {
  double twoalpha = 2.0*alpha;
  double C2alpha = cos(twoalpha);		// cos(2.0*alpha)
  double S2alpha = sin(twoalpha);		// sin(2.0*alpha)

  double Cbeta   = cos(beta);			// cos(beta)
  double Sbeta   = sin(beta);			// sin(beta)
  double Csqbeta = Cbeta*Cbeta;			// cos(beta)*cos(beta)
  double Ssqbeta = Sbeta*Sbeta;			// sin(beta)*sin(beta)

  double twogamma = 2.0*gamma;
  double C2gamma  = cos(twogamma);		// cos(2.0*gamma)
  double S2gamma  = sin(twogamma);		// sin(2.0*gamma)

  double term1 = C2alpha*(3*Ssqbeta + ETA*(1.0+Csqbeta)*C2gamma);
  double term2 = -2.0*ETA*S2alpha*Cbeta*S2gamma;
  double term3 =  3.0*Csqbeta - 1.0 + ETA*Ssqbeta*C2gamma;
  return RT5O96PI*(term1+term2-term3);
  }

// ----------------------------------------------------------------------------
//                           yy Component : Ayy
// ----------------------------------------------------------------------------

/*                                     1/2
                                [  5  ] 
             A  (theta,phi) = - |-----| * [ 1 + eta*cos(2*phi) ]
              yy                [24*PI]
                            

                       1 [                   ]   [1]1/2
          A  (T,P) = - - | A   (T,P) + A     | - |-|    A   (T,P)
           yy          2 [  2,2         2,-2 ]   [6]     2,0                 */

double IntRank2A::AyyPAS() const { return -RT5O24PI*(1.0+ETA); }

double IntRank2A::Ayy()    const
  { return Ayy(_EAs.alpha(), _EAs.beta(), _EAs.gamma()); }

double IntRank2A::Ayy(const EAngles& EA) const
  { return Ayy(EA.alpha(), EA.beta(), EA.gamma()); }

double IntRank2A::Ayy(double alpha, double beta, double gamma) const
  {
  double twoalpha = 2.0*alpha;
  double C2alpha = cos(twoalpha);		// cos(2.0*alpha)
  double S2alpha = sin(twoalpha);		// sin(2.0*alpha)

  double Cbeta   = cos(beta);			// cos(beta)
  double Sbeta   = sin(beta);			// sin(beta)
  double Csqbeta = Cbeta*Cbeta;			// cos(beta)*cos(beta)
  double Ssqbeta = Sbeta*Sbeta;			// sin(beta)*sin(beta)

  double twogamma = 2.0*gamma;
  double C2gamma  = cos(twogamma);		// cos(2.0*gamma)
  double S2gamma  = sin(twogamma);		// sin(2.0*gamma)

  double term1 = -C2alpha*(3*Ssqbeta + ETA*(1.0+Csqbeta)*C2gamma);
  double term2 =  2.0*ETA*S2alpha*Cbeta*S2gamma;
  double term3 = -3.0*Csqbeta + 1.0 - ETA*Ssqbeta*C2gamma;
  return RT5O96PI*(term1+term2+term3);
  }

// ----------------------------------------------------------------------------
//                           zz Component : Azz
// ----------------------------------------------------------------------------

/*                         1/2
                    [  5  ]          2                     2
   A  (theta,phi) = |-----| * [ 3*cos (theta) - 1 + eta*sin (theta)cos(2*phi) ]
    zz              [24*PI]

                                     [ 3 ]1/2
                        A  (T,P) = - | - |    A   (T,P)
                         zz          [ 2 ]     2,0                           */
 
double IntRank2A::AzzPAS() const { return RT5O6PI; }

double IntRank2A::Azz()    const
  { return Azz(_EAs.alpha(), _EAs.beta(), _EAs.gamma()); }

double IntRank2A::Azz(const EAngles& EA) const
  { return Azz(EA.alpha(), EA.beta(), EA.gamma()); }

double IntRank2A::Azz(double alpha, double beta, double gamma) const
  {
  double Cbeta   = cos(beta);			// cos(beta)
  double Sbeta   = sin(beta);			// sin(beta)
  double Csqbeta = Cbeta*Cbeta;			// cos(beta)*cos(beta)
  double Ssqbeta = Sbeta*Sbeta;			// sin(beta)*sin(beta)

  double twogamma = 2.0*gamma;
  double C2gamma  = cos(twogamma);		// cos(2.0*gamma)

  return RT5O24PI*(3.0*Csqbeta-1+ETA*Ssqbeta*C2gamma);
  }

// ----------------------------------------------------------------------------
//                     xy & yx Components : Axy = Ayx
// ----------------------------------------------------------------------------

/*                              1/2
                         [  5  ]    
      A  (theta,phi) = - |-----| * [ eta*cos(theta)sin(2*phi) ] = A
       xy                [24*PI]                                   yx        

                                  1 [                       ]
                    A  (T,P) = -i - | A   (T,P) - A   (T,P) |
                     xy           2 [  2,2         2,1      ]                */


double IntRank2A::AyxPAS() const { return 0.0; }
double IntRank2A::AxyPAS() const { return 0.0; }

double IntRank2A::Ayx()    const
  { return Axy(_EAs.alpha(), _EAs.beta(), _EAs.gamma()); }
double IntRank2A::Axy()    const
  { return Axy(_EAs.alpha(), _EAs.beta(), _EAs.gamma()); }

double IntRank2A::Ayx(const EAngles& EA) const
  { return Axy(EA.alpha(), EA.beta(), EA.gamma()); }
double IntRank2A::Axy(const EAngles& EA) const
  { return Axy(EA.alpha(), EA.beta(), EA.gamma()); }

double IntRank2A::Ayx(double alpha, double beta, double gamma) const
  { return Axy(alpha, beta, gamma); }
double IntRank2A::Axy(double alpha, double beta, double gamma) const
  {
  double twoalpha = 2.0*alpha;
  double S2alpha = sin(twoalpha);		// sin(2.0*alpha)
  double C2alpha = cos(twoalpha);		// sin(2.0*alpha)

  double Cbeta   = cos(beta);			// cos(beta)
  double Sbeta   = sin(beta);			// sin(beta)
  double Csqbeta = Cbeta*Cbeta;			// cos(beta)*cos(beta)
  double Ssqbeta = Sbeta*Sbeta;			// sin(beta)*sin(beta)

  double twogamma = 2.0*gamma;
  double C2gamma  = cos(twogamma);		// cos(2.0*gamma)
  double S2gamma  = sin(twogamma);		// sin(2.0*gamma)

  double term1 = S2alpha*(3*Ssqbeta + ETA*(1.0+Csqbeta)*C2gamma);
  double term2 = 2.0*ETA*C2alpha*Cbeta*S2gamma;
  return RT5O96PI*(term1+term2);
  }

// ----------------------------------------------------------------------------
//                     zx & xz Components : Azx = Axz
// ----------------------------------------------------------------------------

/*                         1/2
                    [  5  ]    
 A  (theta,phi) = - |-----| * sin(theta)*cos(theta)*[3-eta*cos(2*phi)] = A
  xz                [24*PI]                                               zx

                                 1 [                        ]
                    A  (T,P) = - - | A   (T,P) - A    (T,P) |
                     xz          2 [  2,1         2,-1      ]                */


double IntRank2A::AzxPAS() const { return 0.0; }
double IntRank2A::AxzPAS() const { return 0.0; }

double IntRank2A::Azx()    const
  { return Axz(_EAs.alpha(), _EAs.beta(), _EAs.gamma()); }
double IntRank2A::Axz()    const
  { return Axz(_EAs.alpha(), _EAs.beta(), _EAs.gamma()); }

double IntRank2A::Azx(const EAngles& EA) const
  { return Axz(EA.alpha(), EA.beta(), EA.gamma()); }
double IntRank2A::Axz(const EAngles& EA) const
  { return Axz(EA.alpha(), EA.beta(), EA.gamma()); }

double IntRank2A::Azx(double alpha, double beta, double gamma) const
  { return Axz(alpha, beta, gamma); }
double IntRank2A::Axz(double alpha, double beta, double gamma) const
  {
  double Calpha  = cos(alpha);			// cos(alpha)
  double Salpha  = sin(alpha);			// sin(alpha)

  double Cbeta    = cos(beta);			// cos(beta)
  double Sbeta    = sin(beta);			// sin(beta)
  double S2beta   = sin(2.0*beta);		// sin(2.0*beta)

  double twogamma = 2.0*gamma;
  double C2gamma  = cos(twogamma);		// cos(2.0*gamma)
  double S2gamma  = sin(twogamma);		// sin(2.0*gamma)

  double term1 = Calpha*(3.*S2beta - 2.*ETA*Sbeta*Cbeta*C2gamma);
  double term2 = 2.0*ETA*Salpha*Sbeta*S2gamma;
  return RT5O96PI*(term1+term2);
  }

// ----------------------------------------------------------------------------
//                     zy & yz Components : Azy = Ayz
// ----------------------------------------------------------------------------
 
/*                           1/2
                      [  5  ]    
   A  (theta,phi) = - |-----| * eta*sin(theta)*sin(2*phi)] = A
    yz                [24*PI]                                 yz

                                 1 [                        ]
                    A  (T,P) = i - | A   (T,P) - A    (T,P) |
                     yz          2 [  2,1         2,-1      ]                */


double IntRank2A::AzyPAS( ) const { return 0.0; }
double IntRank2A::AyzPAS( ) const { return 0.0; }

double IntRank2A::Azy()    const
  { return Ayz(_EAs.alpha(), _EAs.beta(), _EAs.gamma()); }
double IntRank2A::Ayz()    const
  { return Ayz(_EAs.alpha(), _EAs.beta(), _EAs.gamma()); }

double IntRank2A::Azy(const EAngles& EA) const
  { return Ayz(EA.alpha(), EA.beta(), EA.gamma()); }
double IntRank2A::Ayz(const EAngles& EA) const
  { return Ayz(EA.alpha(), EA.beta(), EA.gamma()); }

double IntRank2A::Azy(double alpha, double beta, double gamma) const
  { return Ayz(alpha, beta, gamma); }
double IntRank2A::Ayz(double alpha, double beta, double gamma) const
  {
  double Calpha   = cos(alpha);			// cos(alpha)
  double Salpha   = sin(alpha);			// sin(alpha)

  double Cbeta    = cos(beta);			// cos(beta)
  double Sbeta    = sin(beta);			// sin(beta)
  double S2beta   = sin(2.0*beta);		// sin(2.0*beta)

  double twogamma = 2.0*gamma;
  double C2gamma  = cos(twogamma);		// cos(2.0*gamma)
  double S2gamma  = sin(twogamma);		// sin(2.0*gamma)

  double term1 = Salpha*(3.*S2beta - 2.*ETA*Sbeta*Cbeta*C2gamma);
  double term2 = 2.0*ETA*Calpha*Sbeta*S2gamma;
  return RT5O96PI*(term1-term2);
  }

// ----------------------------------------------------------------------------
//            All 9 Cartesian Components : Auv with u,v={x,y,z}
// ----------------------------------------------------------------------------

IR2ACart IntRank2A::CartCmpPAS() const { return CartCmp(0,0,0); }

IR2ACart IntRank2A::CartCmp() const { return CartCmp(_EAs); } 

IR2ACart IntRank2A::CartCmp(double A, double B, double G) const
  { EAngles EA(A,B,G); return CartCmp(EA); }

IR2ACart IntRank2A::CartCmp(const EAngles& EA) const
  {
  IR2ACart CC;					// Cartesian components

  double alpha = EA.alpha();			// Angle alpha (radians)
  double beta  = EA.beta();			// Angle beta  (radians)
  double gamma = EA.gamma();			// Angle gamma (radians)

  double twoalpha = 2.0*alpha;
  double Calpha   = cos(alpha);			// cos(alpha)
  double Salpha   = sin(alpha);			// sin(alpha)
  double C2alpha  = cos(twoalpha);		// cos(2.0*alpha)
  double S2alpha  = sin(twoalpha);		// sin(2.0*alpha)

  double Cbeta    = cos(beta);			// cos(beta)
  double Sbeta    = sin(beta);			// sin(beta)
  double S2beta   = sin(2.0*beta);		// sin(2.0*beta)
  double Csqbeta  = Cbeta*Cbeta;		// cos(beta)*cos(beta)
  double Ssqbeta  = Sbeta*Sbeta;		// sin(beta)*sin(beta)

  double twogamma = 2.0*gamma;
  double C2gamma  = cos(twogamma);		// cos(2.0*gamma)
  double S2gamma  = sin(twogamma);		// sin(2.0*gamma)

  double term1, term2, term3;

  term1   = C2alpha*(3*Ssqbeta + ETA*C2gamma*(1.0+Csqbeta));
  term2   = -2.0*ETA*S2alpha*S2gamma*Cbeta;
  term3   = 3.0*Csqbeta - 1.0 + ETA*Ssqbeta*C2gamma;
  CC._Axx =  RT5O96PI*(term1+term2-term3);
  CC._Ayy = -RT5O96PI*(term1+term2+term3);

  term1   = S2alpha*(3*Ssqbeta + ETA*(1.0+Csqbeta)*C2gamma);
  term2   = 2.0*ETA*C2alpha*Cbeta*S2gamma;
  CC._Axy = RT5O96PI*(term1+term2);

  term1   = Calpha*(3.*S2beta - 2.*ETA*Sbeta*Cbeta*C2gamma);
  term2   = 2.0*ETA*Salpha*Sbeta*S2gamma;
  CC._Axz = RT5O96PI*(term1 + term2);

  term1   = Salpha*((3.*S2beta) - 2.*ETA*Sbeta*Cbeta*C2gamma);
  term2   = 2.*ETA*Calpha*Sbeta*S2gamma;
  CC._Ayz = RT5O96PI*(term1-term2);

  return CC;
  }

// ----------------------------------------------------------------------------
//           Cartesian Components : { Axx, Ayy, Azz, Axy, Axz, Ayz }
// ----------------------------------------------------------------------------


row_vector IntRank2A::CartCompsPAS() const
  {
  complex V0  = A20PAS();		// Get A20
  complex V1  = A21PAS();		// Get A21
  complex Vm1 = A2m1PAS();		// Get A2-1
  complex V2  = A22PAS();		// Get A22
  complex Vm2 = A2m2PAS();		// Get A2-2
  row_vector vx(6);
  vx.put(((Vm2+V2)/2.-V0/sqrt(6.)),0);	// Axx = (A22+A2m2)/2 - A20/sqrt(6)
  vx.put((-(Vm2+V2)/2.-V0/sqrt(6.)),1);	// Ayy = -(A22+A2m2)/2 - A20/sqrt(6)
  vx.put((sqrt(2.0/3.0)*V0), 2);       	// Azz = sqrt(2/3)*A20
  vx.put((complex(0,-0.5)*(V2-Vm2)),3);	// Axy = -i*(A22-A21)/2 = Ayx
  vx.put((-0.5 * (V1 - Vm1)), 4);	// Axz = -(A21-A2m1)/2  = Azx
  vx.put((complex(0,0.5)*(V1+Vm1)),5);	// Ayz = i*(A21-A2m1)/2 = Azy
  return vx;
  }

row_vector IntRank2A::CartComps() const { return CartComps(_EAs); }

row_vector IntRank2A::CartComps(const EAngles& EA) const
  {
  row_vector vx(5);
  vx.put(Axx(EA), 0);		// Axx = (A22+A2m2)/2 - A20/sqrt(6)
  vx.put(Ayy(EA), 1);		// Axx = (A22+A2m2)/2 - A20/sqrt(6)
  vx.put(Azz(EA), 2);		// Axx = (A22+A2m2)/2 - A20/sqrt(6)
  vx.put(Axy(EA), 3);		// Axx = (A22+A2m2)/2 - A20/sqrt(6)
  vx.put(Axz(EA), 4);		// Axx = (A22+A2m2)/2 - A20/sqrt(6)
  vx.put(Ayz(EA), 5);		// Axx = (A22+A2m2)/2 - A20/sqrt(6)
  return vx;
  }

// ____________________________________________________________________________
// E             RANK 2 INTERACTION ORIENTATION ACCESS FUNCTIONS
// ____________________________________________________________________________

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

double  IntRank2A::alpha()       const { return _EAs.alpha(); }
double  IntRank2A::beta()        const { return _EAs.beta();  }
double  IntRank2A::gamma()       const { return _EAs.gamma(); }
double  IntRank2A::phi()         const { return _EAs.alpha(); }
double  IntRank2A::theta()       const { return _EAs.beta();  }
double  IntRank2A::chi()         const { return _EAs.gamma();  }
EAngles IntRank2A::orientation() const { return _EAs;         }

void IntRank2A::alpha(double A)                { _EAs.alpha(A); }
void IntRank2A::beta(double  B)                { _EAs.beta(B);  }
void IntRank2A::gamma(double G)                { _EAs.gamma(G); }
void IntRank2A::phi(double   P)                { _EAs.alpha(P); }
void IntRank2A::theta(double T)                { _EAs.beta(T);  }
void IntRank2A::chi(double   C)                { _EAs.gamma(C); }
void IntRank2A::orientation(const EAngles& EA) { _EAs = EA; }
void IntRank2A::orientation(double A, double B, double G, bool deg)
  { _EAs = EAngles(A,B,G,deg); }

void IntRank2A::rotate(const EAngles& EA) { _EAs = EA*_EAs; }

// ____________________________________________________________________________
// F                RANK 2 INTERACTION AUXILIARY FUNCTIONS
// ____________________________________________________________________________

/* Given the three PAS Cartesian components of an irreducible rank 2 spatial
   tensor function AxAyAz will determine the isotropic, anisotropic, and 
   asymmetry terms.  These will be returned in a single coordinate.          */

coord IntRank2A::AisoDelzEta(const coord& AxAyAz)	// Static, not constant
  {
  double xx   = AxAyAz.x();			// Get Axx (PAS)
  double yy   = AxAyAz.y();			// Get Ayy (PAS)
  double zz   = AxAyAz.z();			// Get Azz (PAS)
  double AisoVal = (1.0/3.0)*(xx+yy+zz);	// Get isotropic component
  if(fabs(AisoVal) < 1.e-12) AisoVal = 0.0;	// In sure zero if below delcut
  xx -= AisoVal;				// Remove isotropic component
  yy -= AisoVal;				// before sorting elements
  zz -= AisoVal;
  SortAxAyAz(xx, yy, zz);			// Keep |Azz| >= |Ayy| >= |Axx|
  double EtaVal = (xx-yy)/(zz);			// Calculate the asymmetry
  if(fabs(EtaVal) < 1.e-7)      EtaVal = 0.0;	// -> Small values set zero
  if(fabs(EtaVal - 1.) < 1.e-7) EtaVal = 1.0;	// -> Set exactly 1 if close
  double PASDelzz = zz - 0.5*(xx+yy);		// Calculate delzz value
  PASDelzz *= (2./3.);                          //  THIS IS NOT GAMMA SCALED
  if(fabs(PASDelzz) < 1.e-12) PASDelzz = 0.0;   //  (Small values set zero)
  if(EtaVal == 1) PASDelzz = fabs(PASDelzz);    // Use + delzz if eta is 1
  return coord(AisoVal,PASDelzz,EtaVal); 	//   (where Ayy = - Azz)
  }

void IntRank2A::SortAxAyAz(double& Ax, double& Ay, double& Az)
  {
  double tmp = 0;				// Temp for element sorting
  if(fabs(Ax) > fabs(Az))			// Make sure we abide by rule
    { tmp=Az; Az=Ax; Ax=tmp; }
  if(fabs(Ay) > fabs(Az))			//  |Azz| >= |Ayy| >= |Axx|
    { tmp=Az; Az=Ay; Ay=tmp; }
  if(fabs(Ax) > fabs(Ay))
    { tmp=Ay; Ay=Ax; Ax=tmp; }
  }

bool IntRank2A::CheckEta(double eta, bool warn) const
  {
  if(eta<0. || eta>1.)				// Insure ETA range [0, 1]
    {						// If not the case we will
    if(warn)					// issue warnings as needed
      {
      IR2Aerror(10, 1);				// Assymetry out of range
      IR2Aerror(203,Gform("%6.3f", eta),1);	// Cant set to eta
      }
    return false;
    }
  return true;
  }

bool IntRank2A::PAS( )       const { return _EAs == EAzero; }
bool IntRank2A::Symmetric( ) const { return ETA?false:true; }

matrix IntRank2A::CartMx(double scale) const
  {
  matrix mx;					// For 3x3 matrix
  if(PAS())					// For a tensor in PAS only the
    {						// diagonal elements be nonzero
    mx = matrix(3, 3, complex0, d_matrix_type);	//   Set array as diagonal
    if(Symmetric())				//   If symmetric tensor, eta=0
      {						//   and we know the normalized
      mx.put( 0.5,0,0);				//   elements so they are set
      mx.put(-0.5,1,1);				//   directly
      mx.put( 1.0,2,2);
      }
    else					//   If asymmetric tensor we
      {						//   calculate the elements as
      mx.put(Axx(),0,0);			//   they are filled up
      mx.put(Ayy(),1,1);
      mx.put(Azz(),2,2);
      }
    }
  else						// For a tensor not in its PAS
    {						// the array will be symmetric
    mx = matrix(3, 3, complex0, h_matrix_type);	//   Set the array symmetric
    mx.put(  Axx(),0,0);			//   Calculate and fill up the
    mx.put_h(Axy(),0,1);			//   array in symmetric fashion
    mx.put_h(Axz(),0,2);
    mx.put(  Ayy(),1,1);
    mx.put_h(Ayz(),1,2);
    mx.put(  Azz(),2,2);
    }
  mx *= scale;					// Lastly, scale the normalized
  return mx;					// array as we desire
  }
 
// ____________________________________________________________________________
// H                        PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

/* These functions provide an interface between a rank 2 spatial tensor and
   an external ASCII file in GAMMA parameter set format. Thus, such tensors
   may be specified in external files and easily read into a program. Similarly
   spatial tensors may be quicky stored in ASCII files for future (re)use.   */

// ----------------------------------------------------------------------------
//  Functions To Make A Parameter Set From A Rank 2 Irreducible Spatial Tensor
// ----------------------------------------------------------------------------

/* These functions export information defining an irreducible rank 2 spatial
   tensor to a GAMMA parameter set. There are two items to deal with: the
   spatial tensor asymmetry and its orientation. The easiest output would be
   a double for the eta value and a set of Euler angles for the orientation.

   Unfortunately, such parameters are often not what users would desire in
   expressing the spatial tensor of particular interactions. Furthermore, one
   may desire to use spin or interaction indices so that the parameters may
   be distinguishable from other rank 2 spatial tensors defined in the same
   parameter set. Such differences are mostly accomodated by making these
   functions virtual and allowing specific interactions to deal with any
   common variations. Thus we output only the eta and Euler angles. For the
   direct conversion and parameter set operators no indices will be used.

           Input                IR2A    : Rank 2 spatial tensor (this)
                                filename: Output file name
                                ofstr   : Output file stream
                                iI      : Interaction, 1st spin index (def. -1)
                                iS      :              2nd spin index (def. -1)
                                warn    : Warning level
           Output               none    : Spatial tensor is written as a
                                          parameter set to file filename
                                          or output file stream ofstr        */

IntRank2A::operator ParameterSet( ) const
  { ParameterSet pset; pset += *this; return pset; }

void operator+= (ParameterSet& pset, const IntRank2A &IR2A)
  { IR2A.PSetAdd(pset); }

void IntRank2A::PSetAdd(ParameterSet& pset, int iI, int iS, int warn) const
  {
  string suffx;                                 // Parameter name suffix
  if(iI != -1)					// Only use suffix if iI != -1
    {
    suffx = string("(") + Gdec(iI);		//   Add interaction/spin index
    if(iS != -1) suffx += string(", ")+Gdec(iS); //   2nd spin index if iS != -1
    suffx += string(")");			//   End of suffix
    }

  string pname  = string("IR2eta") + suffx;	// Add asymmetry parameter
  string pstate = string("Asymmetry [0,1]");
  SinglePar par = SinglePar(pname, ETA, pstate);
  pset.push_back(par);

  pname  = string("IR2EAs") + suffx;		// Add Euler angles (orient.) 
  pstate = string("Orientation Euler Angles ")
         + string(" (degrees)");
  double a = _EAs.alpha() * RAD2DEG;
  double b = _EAs.beta()  * RAD2DEG;
  double g = _EAs.gamma() * RAD2DEG;
  coord EA(a, b, g);
  par = EA.param(pname, pstate);                // Parameter for Euler angles
  pset.push_back(par);                          // Add parameter to set
  }

// ----------------------------------------------------------------------------
// Functions To Output Irreducible Spatial Tensor To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

        // Input                IR2A	: Rank 2 spatial tensor (this)
        //                      filename: Output file name
        //                      ofstr   : Output file stream
        //                      iI	: Interaction, 1st spin index (def. -1)
        //                      iS      :              2nd spin index (def. -1)
        //                      warn    : Warning level
        // Output               none    : Spatial tensor is written as a
	//				  parameter set to file filename
	//				  or output file stream ofstr

bool IntRank2A::write(const string &filename, int iI, int iS, int warn) const
  {
  ofstream ofstr(filename.c_str());	// Open filename for input
  if(!ofstr.good())			// If file bad then exit
    {
    IR2Aerror(1, filename);		// Filename problems
    if(warn>1) IR2Afatal(20);		// Fatal error
    else       IR2Aerror(20,1);		// Cannot write to file
    return false;
    }
  if(!write(ofstr, iI, iS, warn?1:0))	// Use overload to write parameters
    {					// If we cannot write these then
    IR2Aerror(40, filename, 1);		// Filename problems
    if(warn>1) IR2Afatal(22);		// Fatal error
    else       IR2Aerror(22,1);		// Cannot write to output stream
    return false;
    }
  ofstr.close();			// Close the file stream
  return true;				// We wrote parameters successfully
  }

bool IntRank2A::write(ofstream& ofstr, int iI, int iS, int warn) const
  {
  ParameterSet pset;                    // Declare a parameter set
  PSetAdd(pset, iI, iS); 		// Add spatial tensor parameters
  if(!pset.write(ofstr, warn?1:0))      // Use parameter set to write
    {
    if(warn)
      {
      if (warn>1) IR2Afatal(23);        // Fatal error
      else        IR2Aerror(23, 1);	// Cannot output parameters
      }
    return false;
    }
  return true;
  }

// ----------------------------------------------------------------------------
//  Functions To Read A Rank 2 Irreducible Spatial Tensor From An ASCII File
// ----------------------------------------------------------------------------

	// Input		IR2A	: Rank 2 spatial tensor
	//			filename: Output file name
	//			idx	: Interaction index (default -1->none)
        //                      warn    : Warning output label
        //                                 0 = no warnings
        //                                 1 = warnings
        //                                >1 = fatal warnings
	// Output		none 	: Rank 2 spatial tensor is read in
	//				  from parameters in file filename

bool IntRank2A::read(const string& filename, int idxI, int idxS, int warn)
  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(filename,warn?1:0))	// Read in pset from file
    {					// If we cannot read parameters
    if(warn)				// then issue warnings or quit?
      {
      IR2Aerror(1, filename, 1);	// Filename problems
      if(warn>1) IR2Afatal(21);		// Fatal error
      else       IR2Aerror(21,1);	// Can't read from parameter file
      }
    return false;
    }
  return read(pset, idxI, idxS, warn); 	// User overloaded function
  }

bool IntRank2A::read(ParameterSet& pset, int idxI, int idxS, int warn)
  {
  string pnb("IR2"); 					// Base name for params
  bool TF = getAX(pset,pnb,ETA,_EAs,idxI,idxS,1);	// Try & read tensor
  if(!TF)						// If unsuccessful,
    {							// make sure we are zero
    ETA  = 0.0;						// & take some action
    _EAs = EAzero;
    if(warn)
      {
      if(warn > 1)  IR2Afatal(21);			// action ala warn flag
      else          IR2Aerror(53);                        //   Setting EAs to 0
      }
    }
  return TF ;
/* This is the old read function that failed to look for cartesian input.
   Delete this once assured that no other IntRank2 classes use it.
  
  if(!getAeta(pset,pnb,ETA,idxI,idxS,warn?true:false) 	// Try & read eta value
  || !CheckEta(ETA,warn?true:false))			// it must be in [0,1]
    {                                                   // If unable take some
    if(warn > 2)  IR2Afatal(21);			// action ala warn flag
    else if(warn) IR2Aerror(11);                        //   Setting eta to 0
    ETA = 0.0;                                          //   Set default 0
    TF = false;
    }
  if(!getOrientation(pset,pnb,_EAs,idxI,idxS,warn?true:false))     // Try & read orient.
    {							// If unable take some
    if(warn > 2)  IR2Afatal(21);			// action ala warn flag
    else if(warn) IR2Aerror(53);                        //   Setting EAs to 0
    _EAs = EAzero;
    TF = false;
    }
*/
  return TF;
  }

// ____________________________________________________________________________
// I                      INTERACTIVE INPUT FUNCTIONS
// ____________________________________________________________________________

string IntRank2A::ask_read(int argc, char* argv[], int& argn, int na)
  {
  string filename;                              // Name of spin system file
  int idxI = -1;				// Interaction or spin index
  int idxS = -1;				// Second spin index
  
  query_parameter(argc, argv, argn,             // Get filename from command
       "\n\tSpin system filename? ", filename); // Or ask for it
  if(na == 1)
    {
    argn++;					// Move to next argument
    query_parameter(argc, argv, argn, 		// Get tensor index
         "\n\tSpatial tensor index? ", idxI);	// Or ask for it
    }
  else if(na == 2)
    {
    argn++;					// Move to next argument
    query_parameter(argc, argv, argn, 		// Get 1st spin index
         "\n\tFirst spin index? ", idxI);	// Or ask for it
    argn++;					// Move to next argument
    query_parameter(argc, argv, argn, 		// Get tensor index
         "\n\tSecond spin index? ", idxS);	// Or ask for it
    }
  read(filename, idxI, idxS);			// Read tensor from filename
  return filename;                              // Give back the filename
  }

/*
void IntRank2A::ask(int argc, char* argv[], int& argq,
                                      double& eta, double& theta, double& phi)

        // Input                IR2A	: Rank 2 spatial tensor (this)
	//			argc    : Number of arguments
	//			argv    : Array of arguments
        //                      argq	: Query value
	//			eta	: Rank 2. asymmetry
	//			theta	: Rank 2. orientation angle
	//			phi	: Rank 2. orientation angle
        // Output               none    : The values of eta, theta, phi
	//				  are set herein
	// Note				: This is INTERACTIVE!

  {
  query_parameter(argc, argv, argq++,			// Read in the coupling
       "\n\tInteraction Asymmetry [0, 1]? ", eta);
  query_parameter(argc, argv, argq++,			// Read in theta angle
  "\n\tOrientation Down From PAS z [1, 180]? ", theta);
  query_parameter(argc, argv, argq++,			// Read in phi angle
   "\n\tOrientation Over From PAS x [0, 360]? ", phi);
  }


void IntRank2A::askset(int argc, char* argv[], int& argq)

        // Input                IR2A	: Rank 2 spatial tensor (this)
	//			argc    : Number of arguments
	//			argv    : Array of arguments
        //                      argq	: Query value
        // Output               none    : IR2A is set interactively
	// Note				: This is INTERACTIVE!

  {
  double eta, theta, phi;
  ask(argc,argv,argq,eta,theta, phi);		// Use the ask function
  *(this) = IntRank2A(eta,theta,phi); 		// Use assignment
  }

 
void IntRank2A::askset( )
 
        // Input                IR2A	: Rank 2 spatial tensor (this)
        // Output               none    : IR2A is set interactively
        // Note                         : This is INTERACTIVE!
 
  {
  int argc = 1;
  char* argv[1];
  int argq = 1000;
  double eta, theta, phi;
  ask(argc,argv,argq,eta,theta,phi);		// Use the ask function
  *(this) = IntRank2A(eta,theta,phi); 	        // Use assignment
  }
*/
 

// ____________________________________________________________________________
// J                             OUTPUT FUNCTIONS
// ____________________________________________________________________________

//----------------------------------------------------------------------------- 
//   Functions That Generate Simple Strings To Simplify & Modularize Printing
//----------------------------------------------------------------------------- 

string IntRank2A::AsymmetryString() const
  { 
  string Seta("Asymmetry:              ");
  Seta += Gform("%10.7f", ETA);
  return Seta;
  }

string IntRank2A::PASString() const
  { 
  string SPAS("Orientation:             ");
         SPAS += "P.A.S.   ";
  return SPAS;
  }

string IntRank2A::AngleString(double angdeg) const
  {
  string SAng = Gform("%7.2f", angdeg) + string(" Deg.");
  return SAng;
  }

string IntRank2A::ThetaString() const
  { return ThetaString(_EAs.beta()*RAD2DEG); }
string IntRank2A::ThetaString(double thedeg) const
  { 
  string Stheta("Down From PAS z-Axis: ");
  Stheta += AngleString(thedeg);
  return Stheta;
  }

string IntRank2A::PhiString() const
  { return PhiString(_EAs.alpha()*RAD2DEG); }
string IntRank2A::PhiString(double phideg) const
  { 
  string Sphi("Over From PAS x-Axis: ");
  Sphi += AngleString(phideg);
  return Sphi;
  }

string IntRank2A::AlphaString() const
  { return AlphaString(_EAs.alpha()); }
string IntRank2A::AlphaString(double alpha) const
  { 
  string Salpha("Euler Angle Alpha:    ");
  Salpha += AngleString(alpha*RAD2DEG);
  return Salpha;
  }

string IntRank2A::BetaString() const
  { return BetaString(_EAs.beta()); }
string IntRank2A::BetaString(double beta) const
  { 
  string Sbeta("Euler Angle Beta:     ");
  Sbeta += AngleString(beta*RAD2DEG);
  return Sbeta;
  }

string IntRank2A::GammaString() const
  { return GammaString(_EAs.gamma()); }
string IntRank2A::GammaString(double gamma) const
  { 
  string Sgamma("Euler Angle Gamma:    ");
  Sgamma += AngleString(gamma*RAD2DEG);
  return Sgamma;
  }

string IntRank2A::AuvString(double Aux, double Auy, double Auz,
                                                   const string& CSForm) const
  { 
  string SAuv;
  SAuv += string("[")   + Gform(CSForm.c_str(), Aux);
  SAuv += string(", ")  + Gform(CSForm.c_str(), Auy);
  SAuv += string(", ")  + Gform(CSForm.c_str(), Auz) + string("]");
  return SAuv;
  } 

//----------------------------------------------------------------------------- 
// Functions That Generate Cartesian Strings To Simplify & Modularize Printing
//----------------------------------------------------------------------------- 

/*       Function CAString                 Function CartAStrings

          [A  , A  , A  ]           [A  , A  , A  ]
          [ xx   xy   xz]           [ xx   xy   xz]   [ x.x, x.x, x.x]
          [A  , A  , A  ]           [A  , A  , A  ] = [ x.x, x.x, x.x]
          [ yx   yy   yz]           [ yx   yy   yz]   [ x.x, x.x, x.x]
          [A  , A  , A  ]           [A  , A  , A  ]
          [ zx   zy   zz]           [ zx   zy   zz]                          */

vector<string> IntRank2A::CAStrings(const string& A) const
  {
  int ns = 6;
  int k  = 0;
  vector<string> Cartstrings(ns);			// Vector of 6 strings
  Cartstrings[k++] = string("[") + A + string("  , ") 	// Set first line
                                 + A + string("  , ")
                                 + A + string("  ]");
  Cartstrings[k++] = "[ xx   xy   xz]";			// Set second line
  Cartstrings[k++] = Cartstrings[0];			// Set third line
  Cartstrings[k++] = "[ yx   yy   yz]";			// Set fourth line
  Cartstrings[k++] = Cartstrings[0];			// Set fifth line
  Cartstrings[k++] = "[ zx   yy   yz]";			// Set sixth line
  return Cartstrings;
  }

vector<string> IntRank2A::CartAStrings(const string& CSForm) const
  {
  IR2ACart CCmps = CartCmp();		// Get Cartesian components
  return CartAStrings(CCmps, CSForm);	// Use function overload
  }

vector<string> IntRank2A::CartAStrings(const EAngles& EA,
                                                    const string& CSForm) const
  {
  IR2ACart CCmps = CartCmp(EA);		// Get Cartesian components
  return CartAStrings(CCmps, CSForm);	// Use function overload
  }

vector<string> IntRank2A::CartAStrings(const IR2ACart& CCmps,
                                                    const string& CSForm) const
  {
//                  Get The Cartesian Component Values

  double Vxx     = CCmps.Axx();		// Axx =  (A22+A2m2)/2 - A20/sqrt(6)
  double Vyy     = CCmps.Ayy();		// Ayy = -(A22+A2m2)/2 - A20/sqrt(6)
  double Vzz     = CCmps.Azz(); 	// Azz = sqrt(2/3)*A20
  double Vxy     = CCmps.Axy();		// Axy = -i*(A22-A2m2)/2 = Ayx
  double Vxz     = CCmps.Axz();		// Axz = -1*(A21-A2m1)/2 = Azx
  double Vyz     = CCmps.Ayz();		// Ayz =  i*(A21-A2m1)/2 = Azy

//           Set Strings For Output Of Cartesian Components

  int ns = 6;
  vector<string> Cartstrings(ns);	// Vector of 6 strings
  vector<string> Cstrs = CAStrings();	// Get array of base strings
  Cartstrings[0] = Cstrs[0];
  Cartstrings[1] = Cstrs[1] + string("   ") + AuvString(Vxx, Vxy, Vxz, CSForm);
  Cartstrings[2] = Cstrs[2] + string(" = ") + AuvString(Vxy, Vyy, Vyz, CSForm);
  Cartstrings[3] = Cstrs[3] + string(" = ") + AuvString(Vxz, Vyz, Vzz, CSForm);
  Cartstrings[4] = Cstrs[4];
  Cartstrings[5] = Cstrs[5];

//                    Standardize The Lengths Of Said Strings

  int i, l, maxl=0;
  for(i=0; i<ns; i++)
    {
    l = (Cartstrings[i]).length();
    if(maxl < l) maxl = l;
    }
  for(i=0; i<ns; i++)
    {
    l = (Cartstrings[i]).length();
    if(maxl > l)
      Cartstrings[i] += string(maxl-l, ' ');
    }
  return Cartstrings;
  }

//----------------------------------------------------------------------------- 
// Functions To Generate Information Strings To Simplify & Modularize Printing
//----------------------------------------------------------------------------- 

/* These functions return spatial tensor information in string format. This is
   done to facilitate printing, in particular printing of spatial tensors from
   within rank 2 interactions.  In this case we make a list of information 
   thats usually displayed to the left of a 3x3 Cartesian matrix rep. of A.

                         Asymmetry:                 x.xxxx
                         Euler Angle Alpha:       xxx.xxxx
                         Euler Angle Beta:        xxx.xxxx 
                         Euler Angle Gamma:       xxx.xxxx
                         Euler Angle Theta:       xxx.xxxx
                         Euler Angle Phi:         xxx.xxxx
                         Orientation:             PAS                        */

//  tmp = string("Asymmetry:            ")
//      + Gform("%7.2f", ETA) + string("        ");
//  InfoStrings.push_back(tmp);

vector<string> IntRank2A::InfoAStrings() const
  {
  vector<string> InfoStrings;
  string tmp;
  InfoStrings.push_back(AsymmetryString());
  if(!PAS())
    {
    InfoStrings.push_back(AlphaString(_EAs.alpha()));
    InfoStrings.push_back(BetaString( _EAs.beta() ));
    InfoStrings.push_back(GammaString(_EAs.gamma()));
    }
  else
    InfoStrings.push_back(PASString());
  return InfoStrings;
  }

vector<string> IntRank2A::InfoAStrings(const EAngles& EA) const
  {
  vector<string> InfoStrings(4);
  InfoStrings[0] = AsymmetryString();
  InfoStrings[1] = AlphaString(EA.alpha());
  InfoStrings[2] = BetaString(EA.beta());
  InfoStrings[3] = GammaString(EA.gamma());
  return InfoStrings;
  }

//----------------------------------------------------------------------------- 
// Functions That Generate Spherical Strings To Simplify & Modularize Printing
//----------------------------------------------------------------------------- 

/* These functions return spatial tensor components in string format.  This is
   done to facilitate printing, in particular printing of spatial tensors from
   within rank 2 interactions.

           Input                IR2A    : Rank 2 spatial tensor (this)
                                CSF     : Float print format (e.g. "%6.3f")
           Output               CSS     : Pointer to array of strings


     SphA2Strings               A    = xx.xxx
       Output                    2,0
                                A    = xx.xxx
                                 2,1
                                A    = xx.xxx
                                 2,-1
                                A    = xx.xxx
                                 2,2
                                A    = xx.xxx
                                 2,-2                                        */


vector<string> IntRank2A::SphA2Strings() const
  {
  vector<string> Sphstrings(10);
  string Aeq("A    = ");
  string Mcm(" 2,");
  string Top;
  string Bottom;
  string mlabs[5]={"0 ","1 ","-1","2 ","-2"};	// Labels for m
  int    mvals[5]={0,1,-1,2,-2};		// Values for m
  for(int m=0; m<5; m++)			// Loop 5 m components
    {
    Top    = Aeq;
    Bottom = Mcm + mlabs[m];
    if(fabs(Im(Acomp(mvals[m]))) > 1.e-6) 	// This for a complex value
      Top += string("(") 
          += string(Gform("%6.3f", Re(Acomp(mvals[m]))))
          += string(" ")
          += string(Gform("%6.3f", Im(Acomp(mvals[m]))))
          += string(")");
    else
      Top += Gform("%6.3f", Re(Acomp(mvals[m])))// This for a real value
          +  string(9, ' ');
    Bottom += string(Top.length()-Bottom.length(), ' ');
    Sphstrings[2*m]   = Top;
    Sphstrings[2*m+1] = Bottom;
    }
  return Sphstrings;
  }

//-----------------------------------------------------------------------------
//          Functions To Print The Tensor In Cartesian Format
//-----------------------------------------------------------------------------

/*	  Output Some Printed Lines Which Will Look Like The Following 

          Irreducible Rank 2 Cartesian Spatial Tensor Components
       (Orientation: Alpha = xxx.xx, Beta = xxx.xx, Gamma = xxx.xx)

                       [A  , A  , A  ] 
                       [ xx   xy   xz]   [ x.x, x.x, x.x]
                       [A  , A  , A  ] = [ x.x, x.x, x.x]
                       [ yx   yy   yz]   [ x.x, x.x, x.x]
                       [A  , A  , A  ]
                       [ zx   zy   zz]

           Input                IR2A	: Rank 2 spatial tensor (this)
                                ostr	: Output stream
				CSF     : Output format of A element
                                EA      : Euler angles for orientation
                                CAS     : Array of strings containing above
	  			tpf     : Title print flag (default 1 = print)
           Output               none    : Rank 2 spatial tensor parameters
                                          sent to the output stream
	   Note				: These are Cartesian components     */

ostream& IntRank2A::printCartesian(ostream& ostr, const string& CSF, int tpf)
  {
  string Spacer;
  string hdr;
  if(tpf)
    {
    hdr = "Irreducible Rank 2 Cartesian Spatial Tensor Components";
    Spacer = string((40-hdr.length()/2), ' ');
    ostr << "\n" << Spacer << hdr;
    if(tpf > 2)
      {
      hdr = "(Oriented Along The Principal Axes)";
      Spacer = string((40-hdr.length()/2), ' ');
      ostr << "\n" << Spacer << hdr;
      }
    ostr << "\n";
    }
  ostr << "\n";
  vector<string> CAS = CartAStrings(CSF);		// Get strings
  return printCartesian(ostr, CAS, tpf);		// Use overload
  }

ostream& IntRank2A::printCartesian(ostream& ostr, 
                                 const EAngles& EA, const string& CSF, int tpf)
  {
  string Spacer;
  string hdr;
  if(tpf)
    {
    hdr = "Rank 2 Cartesian Spatial Tensor Components";
    Spacer = string((40-hdr.length()/2), ' ');
    ostr << "\n" << Spacer << hdr;
    if(tpf > 1)
      {
      hdr = string("( Orientation: ")
         += string("Alpha = ") + Gform("%6.2f", EA.alpha()*RAD2DEG) + ", "
         += string("Beta = ")  + Gform("%6.2f", EA.beta()*RAD2DEG) + ", "
         += string("Gamma = ") + Gform("%6.2f", EA.gamma()*RAD2DEG) + " )";
      Spacer = string((40-hdr.length()/2), ' ');
      ostr << "\n" << Spacer << hdr;
      }
    ostr << "\n";
    }
  ostr << "\n";
  vector<string> CAS = CartAStrings(EA, CSF);		// Get strings
  return printCartesian(ostr, CAS, tpf);		// Use overload
  }

ostream& IntRank2A::printCartesian(ostream& ostr, 
                                            const vector<string>& CAS, int tpf)
  {
  int llength = CAS[0].length();		// Length of input lines
  string Lstart = string("\n");			// Line start
  if(llength < 80)				// Add space if center
    Lstart += string((80-llength)/2, ' ');	// on 80 cols possible
  ostr << Lstart << CAS[0];			// Print 1st line
  ostr << Lstart << CAS[1];			// Print 2nd line
  ostr << Lstart << CAS[2];			// Print 3rd line
  ostr << Lstart << CAS[3];			// Print 4th line
  ostr << Lstart << CAS[4];			// Print 5th line
  ostr << Lstart << CAS[5];			// Print 6th line
  ostr << "\n";
  return ostr;
  }

//-----------------------------------------------------------------------------
//            Here Are The Typically Utilized Output Functions
//-----------------------------------------------------------------------------

ostream& IntRank2A::print(ostream& ostr, vector<string>& CAS, vector<string>& IAS) const
//                vector<string>& CAS, vector<string>& IAS, int fflag, int tpf)
  {
  string Spacer("   ");					// Spacer in middle
  string blank = string(IAS[0].length(), ' ');		// Empty IAS for spacing
  int llength = IAS[0].length() + Spacer.length()	// Length of input lines
              + CAS[0].length();
  string Lstart = string("\n");				// Line start
  if(llength < 80)					// Add space if center
    Lstart += string((80-llength)/2, ' ');		// on 80 cols possible
  int ninfo = IAS.size();				// Number of info lines
  if(ninfo <= 6)					// Less lines on left
    {
    int nb = (6 - ninfo)/2;				// Lines with only right
    int i,j;
    for(j=0; j<nb; j++)
      ostr << Lstart << blank  << Spacer << CAS[j];
    for(i=0; i<ninfo; i++, j++)
      ostr << Lstart << IAS[i] << Spacer << CAS[j];
    for(; j<6; j++)
      ostr << Lstart << blank  << Spacer << CAS[j];
    ostr << "\n";
    }
  else
    {
    int nb = (ninfo-6)/2;
    int i,j;
    for(i=0; i<nb; i++)
      ostr << Lstart << IAS[i];
    for(j=0; j<6; j++, i++)
      ostr << Lstart << IAS[i] << Spacer << CAS[j];
    for(; i<ninfo; i++)
      ostr << Lstart << IAS[i];
    ostr << "\n";
    }
  return ostr;
  }


ostream& IntRank2A::print(ostream& ostr, int fflag, int tpf)

        // Input                IR2A	: Rank 2 spatial tensor (this)
        //                      ostr	: Output stream
	//			fflag   : Format flag
	//				    0 - Basic Parameters
	//				   !0 - Full output
	//			tpf     : Title print flag (default 1 = print)
        // Output               none    : Rank 2 spatial tensor parameters 
        //                                placed into the output stream

  {
/*	    Output Some Printed Lines Which Will Look Like The Following 

                        Irreducible Rank 2 Spatial Tensor
  
                                            [A  , A  , A  ] 
                                            [ xx   xy   xz]   [ x.x, x.x, x.x]
  Asymmetry:                 x.xxxx         [A  , A  , A  ] = [ x.x, x.x, x.x]
  Euler Angle Alpha:       xxx.xxxx         [ yx   yy   yz]   [ x.x, x.x, x.x]
  Euler Angle Beta:        xxx.xxxx         [A  , A  , A  ]
  Euler Angle Gamma:       xxx.xxxx         [ zx   zy   zz]                  */

  if(tpf) 						// Output header as
    {							// specified
    string hdr = "Irreducible Rank 2 Spatial Tensor";
    int hlen   = hdr.length();
    string hsp = string(40 - hlen/2, ' ');
    if(tpf) ostr << "\n" << hsp << hdr << "\n";
    }

  vector<string> CAS = CartAStrings(string("%6.3f"));
  vector<string> IAS = InfoAStrings();
  print(ostr, CAS, IAS);

/*   Output Some Printed Lines Which Will Look Something Like The Following 

  A    = xx.xxx   A    = xx.xxx   A    = xx.xxx   A    = xx.xxx   A    = xx.xxx
   2,0             2,1             2,-1            2,2             2,-2
                                                                             */
  printSpherical(ostr, -1);
  return ostr;
  }

ostream& operator<< (ostream& out, IntRank2A& IR2A) { return IR2A.print(out); }


//----------------------------------------------------------------------------- 
//     Functions That Generate Ouput Of The Spherical Rank 2 Tensor
//----------------------------------------------------------------------------- 

ostream& IntRank2A::print(ostream& ostr, const string& hdr,
                                               const vector<string>& SAS) const

        // Input                IR2A	: Rank 2 spatial tensor (this)
        //                      ostr	: Output stream
	//			hdr	: Header for output
	//			SAS     : Array of strings to be printed
	//				  to the left of the A array(s)
        // Output               none    : Rank 2 spatial tensor parameters 
	// Note				: For pretty output all the strings in
	//				  SAS should have the same length

  {
/*	    Output Some Printed Lines Which Will Look Like The Following 

  			       Some Header Which Was Input
  
  SAS info line 1          xxx.xx unit      [A  , A  , A  ] 
  SAS info line 2            x.xx           [ xx   xy   xz]   [ x.x, x.x, x.x]
         .                    .             [A  , A  , A  ] = [ x.x, x.x, x.x]
         .                    .             [ yx   yy   yz]   [ x.x, x.x, x.x]
         .                    .             [A  , A  , A  ]
         .                    .             [ zx   zy   zz]                  */

  int N = int(SAS.size());
  vector<string> CAS = CartAStrings(string("%6.3f"));	// Get A elements array
  string Spacer("   ");					// Between CAS & SAS
  int llength = SAS[0].length() + Spacer.length()	// Length of input lines
              + CAS[0].length();
  string Lstart = string("\n");				// Line start
  if(llength < 80)					// Add space if center
    Lstart += string((80-llength)/2, ' ');		// on 80 cols possible
  ostr << "\n"; 					// Output the header
  if(hdr.length() < 80)					// Add space if center
    ostr << string((80-hdr.length())/2, ' ');		// on 80 cols possible
  ostr << hdr << "\n";
  int i, j, k, nlines = N;				// Lines to print out
  if(nlines <= 6)					// This if printing
    {							// only 6 lines
    k = (6-N)/2;					// First SAS on row k 
    string blank = string(SAS[0].length(), ' ');	// Empty SAS for spacing
    for(i=0, j=0; i<6; i++)
      {
      if(i<k || j>=N)
        ostr << Lstart << blank << Spacer << CAS[i];
      else
        {
        ostr << Lstart << SAS[j] << Spacer << CAS[i];
        j++;
        }
      }
    }
  else							// This is printing
    {							// N lines
    k = (N-6)/2;					// First CAS on row k 
    for(i=0, j=0; i<N; i++)
      {
      if(i<k || j>=6)
        ostr << Lstart << SAS[i];
      else
        {
        ostr << Lstart << SAS[i] << Spacer << CAS[j];
        j++;
        }
      }
    }
  return ostr;
  }

 
ostream& IntRank2A::printSpherical(ostream& ostr, int tpf)
         
        // Input                IR2A	: Rank 2 spatial tensor (this)
        //                      ostr	: Output stream
	//			tpf     : Title print flag (default 1 = print)
        // Output               none    : Rank 2 spatial tensor parameters
        //                                sent to the output stream

  {
//	                     Output An Optional Header

  string Spacer;
  string hdr;
  if(tpf)
    {
    hdr = "Irreducible Rank 2 Spherical Spatial Tensor Components"; 
    Spacer = string((40-hdr.length()/2), ' ');
    ostr << "\n" << Spacer << hdr;
    if(tpf > 2)
      {
      hdr = "(Oriented Along The Principal Axes)";
      Spacer = string((40-hdr.length()/2), ' ');
      ostr << "\n" << Spacer << hdr;
      }
    }
  ostr << "\n\n";
  

/*   Output Some Printed Lines Which Will Look Something Like The Following 

  A    = xx.xxx   A    = xx.xxx   A    = xx.xxx   A    = xx.xxx   A    = xx.xxx
   2,0             2,1             2,-1            2,2             2,-2
                                                                             */
  vector<string> SBS = SphA2Strings();		// Get our spherical strings
  int i;
  int l1=0, l2=0;
  for(i=0; i<3; i++) l1 +=  SBS[2*i].length();
  for(i=3; i<5; i++) l2 +=  SBS[2*i].length();
  if(l1+l2+12 > 80)				// Use two lines here
    {
    Spacer = string((80 - l1 - 6)/2, ' ');
    ostr << Spacer;
    for(i=0; i<3; i++) ostr << SBS[2*i]   << "   ";
    ostr << "\n";  
    ostr << Spacer;
    for(i=0; i<3; i++) ostr << SBS[2*i+1] << "   ";
    ostr << "\n\n";  
    Spacer = string((80 - l2 - 6)/2, ' ');
    ostr << Spacer;
    for(i=3; i<5; i++) ostr << SBS[2*i]   << "   ";
    ostr << "\n";  
    ostr << Spacer;
    for(i=3; i<5; i++) ostr << SBS[2*i+1] << "   ";
    }
  else						// Use one line here
    {
    ostr << string((80 - l1 - l2 - 12)/2, ' ');
    for(i=0; i<5; i++) ostr << SBS[2*i] << "   ";
    ostr << "\n";  
    ostr << string((80 - l1 - l2 - 12)/2, ' ');
    for(i=0; i<5; i++) ostr << SBS[2*i+1] << "   ";
    }
  ostr << "\n";
  return ostr; 
  }
 
// ____________________________________________________________________________
// K         SPATIAL TENSOR SPHERICAL COMPONENTS FOR POWDER AVERAGES
// ____________________________________________________________________________
 
// ----------------- Parts Which Are Not ETA & Phi Dependent ------------------

 
        // Input                IR2A	: Rank 2 spatial tensor (this)
        //                      Ntheta  : Number of orientation increments
        // Return               A20	: Partial spherical tensor component
	//				  for orientation theta from the 
        //                                tensor PAS over a range of theta
	// Note				: This is sqrt[5/(4*PI)] in PAS
	//				: Add to Bsp0*sin(theta)*sin(theta)
	//				  in order to produce A20 in full
	// Note				: Theta spans [0, 180]

/*                                       1/2
                                  [  5  ]          2 
                 A  (theta,phi) = |-----| * [ 3*cos (theta) - 1 ]
                  20              [16*PI]                                    */


row_vector IntRank2A::Asp20(int Ntheta) const
  {
  row_vector vx(Ntheta);			// Array for values
  double fac = sqrt(5.0/(16.0*PI));		// A20 prefactor
  double fac3 = 3.0*fac;			// This is needed
  double theinc = PI/double(Ntheta-1);		// Angle increment (radians)
  double theta = 0;				// This will be the angle
  for(int i=0; i<Ntheta; i++)			// Loop required increments
    {
    vx.put(fac3*cos(theta)*cos(theta)-fac, i);
    theta += theinc; 
    }
  return vx;
  }

row_vector IntRank2A::Asp21(int Ntheta) const

        // Input                IR2A	: Rank 2 spatial tensor (this)
        //                      theta   : Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Return               A  	: Partial spherical tensor component
        //                       2,1      for orientation theta from the
        //                                tensor PAS over a range of theta
	//				: Add 2 Re(Bpow21)*cos(the))+iIm(Bpow21)
	//				  in order to produce A21 in full
	// Note				: Theta spans [0, 180]
 
/*                            1/2
                       [  5  ]              
         A   (theta) = |-----| * 3.0 * sin(theta) * cos(theta)
          2,1          [24*PI]                                               */
 
  {
  row_vector vx(Ntheta);			// Array for values
  double fac = 3.0*RT5O24PI;			// A21 prefactor
  double theinc = PI/double(Ntheta-1);		// Angle increment (radians)
  double theta = 0;				// This will be the angle
  for(int i=0; i<Ntheta; i++)			// Loop required increments
    {
    vx.put(fac*cos(theta)*sin(theta), i);
    theta += theinc; 
    }
  return vx;
  }


row_vector IntRank2A::Asp22(int Ntheta ) const
 
        // Input                IR2A	: Rank 2 spatial tensor (this)
        //                      theta   : Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Return               A  	: Partial spherical tensor component
        //                       2,2      for orientation theta from the
        //                                tensor PAS over a range of theta
	// Note				: This should be -1/2 in PAS
	//				: Add Re(Bsp22)*(1+cos(theta)cos(theta))
	//				   + iIm(Bsp22B)*cos(theta) in order
	//				  to produce A22 in full
	// Note				: Theta spans [0, 180]
 
/*                                        1/2
                                  3 [  3 ]      2
                   A   A(theta) = - |----|   sin (theta) 
                    2,2           2 [4*PI]                                   */

  {
  row_vector vx(Ntheta);			// Array for values
  double fac = sqrt(15.0/(32.0*PI));		// A22 prefactor
  double theinc = PI/double(Ntheta-1);		// Angle increment (radians)
  double theta = 0;				// This will be the angle
  for(int i=0; i<Ntheta; i++)			// Loop required increments
    {
    vx.put(fac*sin(theta)*sin(theta), i);
    theta += theinc; 
    }
  return vx;
  }


matrix IntRank2A::Asp2s(int Ntheta) const
 
        // Input                IR2A	: Rank 2 spatial tensor (this)
        //                      Ntheta  : Number of orientation increments
        // Return               A  	: Partial spherical tensor component
        //                       2,m      for orientation theta from the
        //                                tensor PAS over a range of theta
	// Note				: Theta spans [0, 180]

  {
  matrix mx(6, Ntheta);				// Array for values
  double theinc = PI/double(Ntheta-1);		// Angle increment (radians)
  double theta = 0;				// This will be the angle
  double fac0 = sqrt(5.0/(16.0*PI));		// A20 prefactor
  double fac03 = 3.0*fac0;			// This is needed
  double fac1 = 3.0*RT5O24PI;			// A21 prefactor
  double fac2 = sqrt(15.0/(32.0*PI));		// A22 prefactor
  double Ctheta, Stheta;
  for(int j=0; j<Ntheta; j++)			// Loop required increments
    {
    Ctheta = cos(theta);
    Stheta = sin(theta);
    mx.put(fac03*Ctheta*Ctheta-fac0, 0, j);	// 	Row of A20 parts
    mx.put(fac1*Ctheta*Stheta, 1, j);		//	Row of A21 parts
    mx.put(-mx.get(1,j), 2, j);			//	Row of A2-1 parts
    mx.put(fac2*Stheta*Stheta, 3, j);		//	Row of A22 parts
    mx.put(mx.get(3,j), 4, j);			//	Row of A2-2 parts
    mx.put(Stheta, 5, j);			//	Row of sin(theta)
    theta += theinc; 
    }
  return mx;
  }

// ____________________________________________________________________________
// L                  SPATIAL TENSOR COMPARISON FUNCTIONS
// ____________________________________________________________________________

/* These functions allow users to check if two spatial tensor are equivalent.
   Since eta is the only value that defines our irreducible rank 2 normalized
   spatial tensors, they will be equivalent if the asymmetry values match.

           Input                IR2A    : Rank 2 spatial tensor (this)
                                IR2A2   : Another rank 2 spatial tensor
           Output               T/F     : TRUE if IR2A2 and this are equal    */
 
bool IntRank2A::operator==(const IntRank2A &IR2A2) const
  { return (!(fabs(ETA-IR2A2.ETA) > 1.e-9)); }
 
bool IntRank2A::operator!=(const IntRank2A &IR2A2) const
  { return (fabs(ETA-IR2A2.ETA) > 1.e-9); }


#endif							// IntRank2A.cc
