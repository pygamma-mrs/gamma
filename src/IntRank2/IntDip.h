/* IntDip.h *****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Dipole Dipole Interaction 	  	   Interface		**
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
** A variable of type IntDip represents a dipole-dipole interaction 	**
** defined for a particular nuclear spin pair.  The interaction is 	**
** maintained in the form of spatial and spin tensors which are stored 	**
** in irreducible spherical format.					**
**                                                                      **
** 1.) rank = 2         2.) symmetric           3.) eta = [0, 1]        **
**                                                                      **
** The tensor will internally be stored in irreducible spherical form   **
** at any chosen orientation.  The tensor principal axes will also be   **
** maintiained as well as a set of Euler angles that indicate the PAS   **
** orientation relative to some common axes.  Normally, eta is zero	**
** for dipolar interactions (since the interaction is symmetric about	**
** the internuclear vector)						**
**                                                                      **
** The interaction blends a rank 2 spin tensor (for a spin pair)	**
** with a rank 2 spatial tensor and a dipolar interaction constant.	**
** This provides the ability to generate oriented dipolar Hamiltonians.	**
**                                                                      **
** In addition to the functions contained in the base spatial tensor    **
** class, this IntDip provides functions for building up the tensor   	**
** and accessing the tensor elements from a "dipolar" standpoint.   	**
**                                                                      **
** The following defintions are used herein (Auv normlized by delzz):	**
**                                                                      **
** 1.) PAS: |Azz| >= |Ayy| >= |Axx|					**
** 2.) PAS: Azz=1, eta=(Ayy-AXX)/2, Axx=(1+eta)/2, Ayy=(eta-1)/2        **
**                                                                      **
**                                                                      **
*************************************************************************/

#ifndef   IntDip_h_			// Is file already included?
#  define IntDip_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface
#endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Level1/coord.h>		// Include coordinates
#include <Level2/EAngles.h>		// Include Euler angles
#include <IntRank2/IntRank2A.h>		// Include base class headers
#include <IntRank2/IntRank2T.h>		// Include base class headers
#include <IntRank2/IntRank2.h>		// Include base class headers
#include <Matrix/matrix.h>		// Include GAMMA matrices
#include <Matrix/row_vector.h>		// Include GAMMA row vectors
#include <string>			// Include stdlibc++ strings
#include <fstream>			// Include stdlibc++ file streams
#include <list>				// Include stdlibc++ STL lists


//forward declarations
class IntDip;
MSVCDLL matrix HD0(double qn, double wCo, double eta=0, 
                                                 double theta=0, double phi=0);
MSVCDLL matrix HD1(double Om, double qn, double wCo, double eta,
                                             double theta=0.0, double phi=0.0);
class IntDip: public IntRank2
  {
  matrix T20wh;				// Spherical T20 weak heteronuclear
  double _DCC;				// Dipolar coupling constant (DELZZ)

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i              DIPOLE DIPOLE INTERACTION ERROR HANDLING
// ____________________________________________________________________________

/*       Input                D       : Dipolar interaction (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

         void IDerror(int eidx,                           int noret=0) const;
         void IDerror(int eidx, const std::string& pname, int noret=0) const;
volatile void IDfatal(int eidx)                                        const;
volatile void IDfatal(int eidx, const std::string& pname)              const;

// ____________________________________________________________________________
// ii     DIPOLE DIPOLE INTERACTION SETUP FUNCTIONS USING PARAMETER SETS
// ____________________________________________________________________________

/* These functions set up specific aspects of a dipolar interaction.  Since
   the various interaction parameters interact, the functions MUST be private
   because their misuse could produce an inconsistent interaction

   The goal here is quite simple.  We must determine the following set of
   values for each dipolar interaction: { Iqn,Sqn,DCC,eta,alpha,beta,gamma }
   Complexity arises because we allow that a variety of parameters may be used
   to define these values. Additionally we allow that some parameters may be
   zero (e.g. eta) and that defaults may automatically be used. But, the end
   result remains the same, we seek the set of values that define a dipolar
   interaction in GAMMA.                                                     */

// ----------------------------------------------------------------------------
//                       Complete Dipolar Interaction
// ----------------------------------------------------------------------------

/* These functions will try and get all the parameters required to define a
/  dipolar interaction: { Iqn,Sqn,DCC,eta,alpha,beta,gamma }. Although the 1st
/  function is the generic means to accomplishing this, we will immediately
/  branch into two categories: 1.) Defining the interaction with a single
/  interaction index & 2.) Defining the interaction with two spin indices.
                                                 int idxI, int idxS, int warn)

          Input                D       : Dipolar interaction (this)
                               pset    : A parameter set
                               idxI    : Index of first spin or interaction
                               idxS    : Index of first spin
                               warn    : Warning level
                                          0 - no warnings
                                          1 - warnings
                                         >1 - fatal warnings
          Output               TF      : Dipolar interaction is set
           Output               TF      : Dipolar interaction is set
                                                TF = 0 couldn't read it
                                                TF = 1 read with iso & coord
                                                TF = 2 read with iso & DCC   */

bool getDI(const ParameterSet& pset,
         double& Iqn, double& Sqn, double& dcc, double& eta, EAngles& EA,
                                                 int idxI, int idxS, int warn);
bool getDI1(const ParameterSet& pset,
         double& Iqn, double& Sqn, double& dcc, double& eta, EAngles& EA,
                                                           int idxS, int warn);
bool getDI2(const ParameterSet& pset,
         double& Iqn, double& Sqn, double& dcc, double& eta, EAngles& EA,
                                                 int idxI, int idxS, int warn);

bool getDI2(const ParameterSet& pset,
         double& dcc, double& eta, EAngles& EA,
         const Isotope& ISI, const Isotope& ISS, int idxI, int idxS, int warn);

// ----------------------------------------------------------------------------
//        Dipolar Interaction Setup Functions Using Spin Pair Indices
// ----------------------------------------------------------------------------
 
// **********        Set Up/Read This Spin Pair Isotope Types        **********
 
/*
int getIsos(const ParameterSet& pset, int idxI, int idxS,
                             std::string& IsoI, std::string& IsoS, int warn=1);
*/
 
        // Input                D       : Dipolar interaction (this)
        //                      pset    : A parameter set
        //                      idxI    : Index of first spin
        //                      idxS    : Index of first spin
        //                      IsoI    : String for spin I isotope type
        //                      IsoS    : String for spin S isotope type
        //                      warn    : Warning level
        //                                 0 - no warnings
        //                                 1 - warnings
        // Output               TF      : Values of IsoI & IsoS set from pset
        //                                Returns is true if both set 


// ----------------------------------------------------------------------------
//         Get Dipolar Coupling Constant Value From A Parameter Set
// ----------------------------------------------------------------------------
 
/*                     Try & Read D Dipolar Coupling Constant
 
                                      gamma * gamma * hbar
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
					  DCC is set in Hz prior to exit if true.

                         Currently Allowed DCC Parameters

          DCC,    DCC(#),    DCC(#,#)      - Dipolar Coupling in kHz
          DCCkHz, DCCkHz(#), DCCkHz(#,#)   - Dipolar Coupling in kHz
          DCCKHz, DCCKHz(#), DCCKHz(#,#)   - Dipolar Coupling in kHz
          DCCHz,  DCCHz(#),  DCCHz(#,#)    - Dipolar Coupling in Hz          */  

bool getDCC(const ParameterSet& pset, double& dcc,
                                                 int idxI, int idxS, bool warn);




// **********        Set Up/Read The Whole Dipolar Interaction       **********


int setDI(const ParameterSet& pset, int idxI, int idxS, int warn=2);

        // Input                D       : Dipolar interaction (this)
        //                      pset    : A parameter set
        //                      idxI    : Index of first spin
        //                      idxS    : Index of first spin
        //                      warn    : Warning level
        //                                 0 - no warnings
        //                                 1 - warnings
        //                                >1 - fatal warnings
        // Output               DCC     : Dipolar interaction is set
        //                                from parameters in pset


// ----------------------------------------------------------------------------
//          Set Dipolar Interaction Spatial and Spin Tensor Components
// ----------------------------------------------------------------------------

// void setAs();	INHERITED	 Set spatial tensor components
// void setTs();	INHERITED	 Done in IntRank2T base class
//                                            +
//                                         m  |
//                              T    = (-1)  T
//                               2,m          2,-m

//         T  = T0;   T  = T1;   T   = Tm1;   T  = T2;   T   = Tm2
//          2,0        2,1        2,-1         2,2        2,-2

//                   1/2
//                [4]                      1
//          T   = |-| * I         T    = - - I           T    = 0
//           2,0  [6]    z         2,1     2  +           2,2


void setT20wh();

        // Input                D       : Dipolar interaction (this)
        // Output               none    : Dipolar interaction weak
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

//          1/2                                1/2             1/2
//       [1]                      weak      [1]             [4]
// T   = |-| * [3I S - I.S]   ------------> |-| * [2I S ] = |-| I S
//  2,0  [6]      z z         heteronuclear [6]      z z    [6]  z z

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A          DIPOLE DIPOLE INTERACTION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

MSVCDLC IntDip();
MSVCDLC IntDip(const IntDip &D1);

// ----------------------------------------------------------------------------
//              Direct Constructors Using Spherical Components
// ----------------------------------------------------------------------------

/* Here we need to know the quantum numbers of the two spins involved, the
   dipolar coupling, and the interaction asymmetry (usually 0). Shoud a non-
   zero asymmetry actually be needed it is unitless and in the range [0,1]   */

        // Input                D       : Dipolar interaction (this)
        //                      IsoI    : Spin I isotope type
        //                      IsoS    : Spin S isotope type
	//			Iqn	: Spin I quantum number
	//			Sqn	: Spin S quantum number
        //                      DCC     : Dipolar coupling (Hz)
        //                      eta     : Tensor asymmetry value (default 0)
        // Output               none    : Dipolar interaction constructed
        // Note                         : A fatal error will result if
        //                                either IsoI & IsoS are a mix of an
        //                                electron and a nucleon
        // Note                         : A non-fatal error will occur if
        //                                the sign of DCC doesn't match the
        //                                sign of gamma(IsoI)*gamma(IsoS)
        // Note                         : For dipolar interactions DCC=delzz
 
MSVCDLC IntDip(const std::string& IsoI, const std::string& IsoS,
                         double DCC, double eta=0.0, const EAngles& EA=EAzero);

MSVCDLC IntDip(const Isotope&     IsoI, const Isotope&     IsoS,
                         double DCC, double eta=0.0, const EAngles& EA=EAzero);

MSVCDLC IntDip(double              Iqn, double              Sqn,
                         double DCC, double eta=0.0, const EAngles& EA=EAzero);

// ----------------------------------------------------------------------------
//              Direct Constructors Using Spin Coordinates
// ----------------------------------------------------------------------------

/* Here we need to know the quantum numbers of the two spins involved and spin
   coordinates from which we can determine the dipolar coupling. The
   interaction asymmetry will automatically be set to 0.

           Input                D       : Dipolar interaction (this)
                                II      : Spin I isotope type
                                IS      : Spin S isotope type
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
 
MSVCDLC IntDip(const std::string& IsoI, const std::string& IsoS,
                                          const coord& pt1, const coord& pt2);
MSVCDLC IntDip(const Isotope&     IsoI, const Isotope&     IsoS,
                                          const coord& pt1, const coord& pt2);

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

MSVCDLC IntDip(const std::string& II, const std::string& IS,
                                const coord& DxDyDz, const EAngles& EA=EAzero);
MSVCDLC IntDip(const Isotope& II, const Isotope& IS,
                                const coord& DxDyDz, const EAngles& EA=EAzero);
MSVCDLC IntDip(double Iz, double Sz, 
                                const coord& DxDyDz, const EAngles& EA=EAzero);

// ----------------------------------------------------------------------------
//                     Constructors Using Parameter Sets
// ----------------------------------------------------------------------------

/* These constructors attempt to set the dipolar interaction from parameters
   in a given Parameter set. Either a single interaction index may be supplied
   or two spin indices. An additional variation (in support of spin systems)
   takes two spin isotopes.

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

MSVCDLC IntDip(const ParameterSet& pset, int idxI, int idxS=-1, int warn=2);
//IntDip(const ParameterSet& pset, double Iz, double Is,
//                                         int idxI, int idxS=-1, int warn=2);


/*
IntDip(double Iqn, double Sqn, const ParameterSet& pset,
                                                       int idx=-1, int warn=2);
*/
 

// ----------------------------------------------------------------------------
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------

MSVCDLL void operator= (const IntDip &D1);
MSVCDLC ~IntDip();

// ____________________________________________________________________________
// B                     SPATIAL TENSOR COMPONENT ACCESS
// ____________________________________________________________________________
 
//-----------------------------------------------------------------------------
//                   Anisotropy Access (Dipolar Coupling)
//-----------------------------------------------------------------------------

/* By definition, the dipolar coupling constant is identical to the spatial
   tensor delzz value in GAMMA and this is 2/3 the tensor anisotropy. The
   relationship is formally given by

                  mu
                 --- * hbar * gamma  * gamma               1/2
                 4pi               i        j     1   [ 5  ]      D   2       D
 delzz = DCC   = ____________________________ = - - * |----|  * Xi  = - DA = w
            ij                 3                  2   [6*pi]          3
                           r
                            ij

   The above is in units of radians/sec. One should divide by 2*PI if values
   in Hz are desired.
   
   The functions here allow users to get and/or set the dipolar coupling
   constant as well as calculate it based on the formula above.

           Input                D       : Dipolar interaction
                                IsoI    : Spin I isotope type (1H, 2H, ...)
                                IsoS    : Spin S isotope type (1H, 2H, ...)
           Return               Rij     : Dipolar distance, in meters
           Note                         : The sign of DCC depends on the
                                          product of gyromagnetic ratios
                                          of II && SS
           Note                         : The value of delzz is equivalent to
                                          the dipolar coupling constant DCC
           Note                         : Some functions are static & do not
                                          need an instance of IntDip to work
                                                    -1      -2
           Note                         : 1T = 1 J-C  -sec-m                 */


MSVCDLL double CheckDCC(const Isotope& II, const Isotope& IS, double dcc);
MSVCDLL double CheckDCC(const Isotope& II, const Isotope& IS,
                                           const coord& pt1, const coord& pt2);

MSVCDLL static double DCC(const std::string& I, const std::string& S, double R);
MSVCDLL static double DCC(const Isotope&     I, const Isotope&     S, double R);
MSVCDLL static double DCC(const std::string& I, const std::string& S,
                                           const coord& ptI, const coord& ptS);
MSVCDLL static double DCC(const Isotope&     I, const Isotope&     S,
                                           const coord& ptI, const coord& ptS);
 
MSVCDLL double DCC() const;
MSVCDLL void   DCC(double dz);
                            
MSVCDLL static double W2DCC(const std::string& I, const std::string& S, double W);
MSVCDLL static double DCC2W(const std::string& I, const std::string& S, double D);
MSVCDLL static double W2DCC(const Isotope&     I, const Isotope&     S, double W);
MSVCDLL static double DCC2W(const Isotope&     I, const Isotope&     S, double D);

/*  Note that the default base class IntRank2A also provides the GAMMA
    normalized anisotropy and delzz values. These are constants given by

                       1/2
                  [ 5 ]                  ^      3
          del   = |---]  = 0.51503      /_\ A = - del   = 0.77255
             zz   [6PI]                         2    zz

    Addiontal functions that are standard for all rank 2 interacitons are
    found below.                                                            */

// static double IntRank2A::delzz();                            INHERITED
// static double IntRank2A::delA();                             INHERITED

MSVCDLL double aniso() const;
MSVCDLL void   aniso(double aiso);
MSVCDLL double DA()    const;
MSVCDLL void   DA(double aiso);

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
//	         Normalized Cartesian Tensor Component Access
//-----------------------------------------------------------------------------

/*
  double IntRank2A::AxxPAS() const;				INHERITED
  double IntRank2A::AyyPAS() const;				INHERITED
  double IntRank2A::AzzPAS() const;				INHERITED
  double IntRank2A::AyxPAS() const;				INHERITED
  double IntRank2A::AxyPAS() const;				INHERITED
  double IntRank2A::AzxPAS() const;				INHERITED
  double IntRank2A::AxzPAS() const;				INHERITED
  double IntRank2A::AzyPAS() const;				INHERITED
  double IntRank2A::AyzPAS() const;				INHERITED

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

  row_vector IntRank2A::CartCompsEA() const;			INHERITED
  row_vector IntRank2A::CartComps()   const;			INHERITED
  row_vector IntRank2A::CartComps(double theta, double phi=0) const;        */

//-----------------------------------------------------------------------------
//            Cartesian Tensor Component Access - Not Normalized
//-----------------------------------------------------------------------------

MSVCDLL double dxx() const;
MSVCDLL double dyy() const;
MSVCDLL double dzz() const;
MSVCDLL double dxy() const;
MSVCDLL double dyx() const;
MSVCDLL double dxz() const;
MSVCDLL double dzx() const;
MSVCDLL double dyz() const;
MSVCDLL double dzy() const;

MSVCDLL double dxx(double T, double P) const;
MSVCDLL double dyy(double T, double P) const;
MSVCDLL double dzz(double T, double P) const;
MSVCDLL double dyx(double T, double P) const;
MSVCDLL double dxy(double T, double P) const;
MSVCDLL double dzx(double T, double P) const;
MSVCDLL double dzy(double T, double P) const;
MSVCDLL double dxz(double T, double P) const;
MSVCDLL double dyz(double T, double P) const;

MSVCDLL double dxx(double A, double B, double G) const;
MSVCDLL double dyy(double A, double B, double G) const;
MSVCDLL double dzz(double A, double B, double G) const;
MSVCDLL double dyx(double A, double B, double G) const;
MSVCDLL double dxy(double A, double B, double G) const;
MSVCDLL double dzx(double A, double B, double G) const;
MSVCDLL double dzy(double A, double B, double G) const;
MSVCDLL double dxz(double A, double B, double G) const;
MSVCDLL double dyz(double A, double B, double G) const;

MSVCDLL double dxx(const EAngles& EA) const;
MSVCDLL double dyy(const EAngles& EA) const;
MSVCDLL double dzz(const EAngles& EA) const;
MSVCDLL double dyx(const EAngles& EA) const;
MSVCDLL double dxy(const EAngles& EA) const;
MSVCDLL double dzx(const EAngles& EA) const;
MSVCDLL double dzy(const EAngles& EA) const;
MSVCDLL double dxz(const EAngles& EA) const;
MSVCDLL double dyz(const EAngles& EA) const;

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
  double  IntRank2A::alpha()       const;			INHERITED
  double  IntRank2A::beta()        const;			INHERITED
  double  IntRank2A::gamma()       const;			INHERITED
  double  IntRank2A::phi()         const;			INHERITED
  double  IntRank2A::theta()       const;			INHERITED
  EAngles IntRank2A::orientation() const;			INHERITED

  void IntRank2A::alpha(double A);				INHERITED
  void IntRank2A::beta(double  B);				INHERITED
  void IntRank2A::gamma(double G);				INHERITED
  void IntRank2A::phi(double   P);				INHERITED
  void IntRank2A::theta(double T);				INHERITED
  void IntRank2A::orientation(const EAngles& EA);		INHERITED    */


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

        1.) GAMMA normalized Auv - Done With CartMx()
        2.) Typical Duv values   - Done With Smx(true);
        3.) Shown in lab frame   - Done With This Function

   For case 3.) the values are related to the GAMMA normalized (Auv) and
   typically presented values (Duv) according to

                                1/2
                        [ 6*PI ]
                  D   = | ---- |    * del   * A   + Kdel 
                   uv   [  5   ]         zz    uv       u,v 

   where Kdel is a Kronecker delta function.                                 */

MSVCDLL matrix Dmx() const;
 
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

MSVCDLL matrix T20het() const;
MSVCDLL matrix T20het(const std::vector<int>& HSs, int i, int j) const;
 
        // Input                D       : Dipolar interaction
        // Output               T       : Spherical spin tensor component
        //                       2,0      for m=0 under high-field assuming
        //                                the spins are heteronuclear
 
// ____________________________________________________________________________
// D             DIPOLAR INTERACTION CONSTANT ACCESS FUNCTIONS
// ____________________________________________________________________________

/* The dipolar interaction constant is what sets the overall interaction
   strength and allows GAMMA to track different interaction types using
   "normalized" spatial and spin tensor (classes IntRank2A and IntRank2T). For
   dipolar interactions the interaction constant is defined as (in radians/sec)

                          1/2
                    [6*pi]     mu        2
               -2 * |----|   * --- * hbar * gamma  * gamma
          D         [ 5  ]     4pi               i        j
        xi   = ____________________________________________ = del   = DCC
          ij                       3                             zz
                                  r                                          */


MSVCDLL double xi(bool Hz=false) const;
 
        // Input                D	: Dipolar interaction
        // Return               xi      : D. interaction constant
        // Note                           Value is in radians/sec
 
// ____________________________________________________________________________
// E               DIPOLAR INTERACTION AUXILIARY FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                     Internuclear Distance Functions
// ----------------------------------------------------------------------------

MSVCDLL double R(const std::string& II, const std::string& SS, bool gchk=1) const;

        // Input                D       : Dipolar interaction (this)
        //                      II      : Spin I isotope type (1H, 2H, ...)
        //                      SS      : Spin S isotope type (1H, 2H, ...)
	//			gchk	: Check spin types
        // Return               R       : Dipolar distance, Angstroms
        //                                          -1      -2
        // Note                         : 1T = 1 J-C  -sec-m


/*                     [ mu                           ]1/3
                       | --- * hbar * gamma  * gamma  |
                       | 4pi               i        j |        [  NUM  ]
                 r   = | ____________________________ |     =  | ----- |
                  ij   |                              |        [ DENOM ]
                       |           DCC                |  
                       [             ij               ]                      */

// ----------------------------------------------------------------------------
//                    Transition Splitting Frequencies
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
 

MSVCDLL double W(const std::string& II, const std::string& SS, bool gchk=true) const;

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

    
MSVCDLL void W(double w, const std::string& II, const std::string& SS, int gchk=1);
 
        // Input                D       : Dipolar interaction
        //                      w       : Dipolar frequency (Hz)
        //                      II      : Spin I isotope type (1H, 2H, ...)
        //                      SS      : Spin S isotope type (1H, 2H, ...)
        //                      gchk    : Check if Iz & Sz Match II & SS
        //                                      0 = don't bother
        //                                      1 = make sure
        // Return               void    : The Dipolar frequency (as DELZZ)
        //                                of the interaction is set

// ----------------------------------------------------------------------------
//                 Spin & Interaction Index Checking Functions
// ----------------------------------------------------------------------------

MSVCDLL bool DCheck(const std::string& I, const std::string& S, int warn=2) const;

MSVCDLL int Dspincheck(double Iqn, double Sqn) const;

	// Input		D	: Dipolar interaction
	//			Iqn,Sqn : Spin quantum values of I&S
	// Output		TF      : True if both I & S have spin
	//				  quantum values that are non-zero
	//				  positive multipes of 1/2


MSVCDLL int Dindexcheck(int Iidx, int Sidx) const;

        // Input                D       : Dipolar interaction
        //                      idxI       : Spin isotope type
        //                      idxS       : Spin isotope type
        // Output               TF      : True if both I & S have spin
        //                                idices that are appropriate
 

// ____________________________________________________________________________
// F                        PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//        Functions To Make A Parameter Set From A Dipolar Interaction
// ----------------------------------------------------------------------------

/* Note that the preferred means of specifying a dipolar interaction when there
   is no prior knowledge of spin isotope types nor spin coordinates is used for
   filling the parameter set.  The base parameters of individual interactions
   in that case are { DI(#), DS(#), DCC(#), DTheta(#), DPhi(#), DEta(#) }.  If
   there is no eta value, as is typical for most dipolar interactions, the
   DEta parameter will not be included.


           Function				Result
         ============  ========================================================
         ParameterSet  Converts Dipolar Interaction Into A Parameter Set
         operator+=    Adds Dipolar Interaction Parameters into Parameter Set
         PSetAdd       Similar to above but allows for indices (prefix/suffix)

        Input                D	: Dipolar interaction (this)
                             pset	: GAMMA parameter set
                             idx     : Interaction index (default -1)
                             pfx     : Interaction 2nd indx (def -1)         */

MSVCDLL              operator ParameterSet( ) const;
MSVCDLL friend void  operator+= (ParameterSet& pset, const IntDip &D);
MSVCDLL        void  PSetAdd(ParameterSet& pset, int idx=-1, int pfx=-1) const;

// ----------------------------------------------------------------------------
// Functions To Output Dipolar Interaction To An ASCII File (via Parameter Set)
// ----------------------------------------------------------------------------

/* These functions write parameters defining a dipolar interaction into a
   specified ASCII file. 
   is no prior knowledge of spin isotope types nor spin coordinates is used for
   filling the parameter set.  The base parameters of individual interactions
   in that case are { DI(#), DS(#), DCC(#), DTheta(#), DPhi(#), DEta(#) }.  If
   there is no eta value, as is typical for most dipolar interactions, the
   DEta parameter will not be included.

           Input                D	: Dipolar interaction (this)
                                fn	: Output file name
                                ofstr   : Output file stream
                                idx     : Interaction index (default -1)
                                pfx     : Interaction 2nd indx (def -1)
                                wrn	: Warning level
           Output               none    : Dipolar interaction parameters
                                          written in parameter set format to
                                          file filename or ostream ofstr     */

MSVCDLL int write(const std::string &fn,int idx=-1,int pfx=-1,int wrn=2) const;
MSVCDLL int write(std::ofstream& ofstr, int idx=-1,int pfx=-1,int wrn=2) const;
 
// ____________________________________________________________________________
// G                   DIPOLAR INTERACTION INPUT FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//          Functions To Read ASCII File Using Two Spin Indicies
// ----------------------------------------------------------------------------
 
/* These read functions utilize either two spin indices, one interaction
   index, or no interaction index. If the first index, idxI, is left as the
   default (-1) it will be assumed that no interaction index is used. If the
   value of idxI is not -1 but the value of idxS is the default of -1 then we
   assume that idxI is the interaction index. If neither idxI or idxS is left
   at their default value we use each as a spin index. 

   If the spin pair is specifiled for the dipolar interaction the we'll try to
   read the parameter set
 
                    { Iso(i), Iso(j), Coord(i), Coord(j) }
   
   If an interaction index is specified for the dipolar interaction the we'll
   try to read the parameter set
 
                           { DI, DS, DCC, DEAs, Deta }
 
   The functions do NOT allow for any additional suffix, since here the suffix
   is used by the spin or interaction index. The functions do NOT allow for 
   the prefix [#] either, because multiple dipole interactions can be defined 
   in the same file by switching spin pair or interaction indices.  Multiple
   sets of dipolar interactions can be read using dipole interaction vectors
   (see class IntDipVec.)  The asymmetry, eta, is assumed to be zero (normally
   the case with dipolar interactions) 

   When using the interaction index the parameters Iso and Coord are not
   allowed as these are intrinsically associated with single spins. If the 
   parameters DI and DS are not present in the file/parameter set, DI=DS=1/2 
   are used by default. 

           Input                D       : Dipolar interaction
                                filename: Output file name
                                pset    : Parameter set
                                IdxI    : Index for spin I or interaction
                                IdxS    : Index for spin S
	  			warn    : Warning output label
	  				   0 = no warnings
	  				   1 = warnings
	  				  >1 = fatal warnings
           Output               TF      : Dipolar interaction is read in
                                          from parameters in file filename
	  				  Return is true if read properly
                                                TF = 0 couldn't read it
                                                TF = 1 read with iso & coord
                                                TF = 2 read with iso & DCC   */

MSVCDLL bool read(const std::string &filename,int iI=-1,int iS=-1,int warn=2);
MSVCDLL bool read(const ParameterSet& pset,   int iI=-1,int iS=-1,int warn=2);
 
MSVCDLL void scan(const std::string &filename,int iI=-1,int iS=-1,int warn=2);
MSVCDLL void scan(const ParameterSet& pset,   int iI=-1,int iS=-1,int warn=2);

// ----------------------------------------------------------------------------
//       Functions To Read ASCII File Using One Interaction Index
// ----------------------------------------------------------------------------
 
//int read(const std::string &filename, int idx=-1, int warn=2);
//void read(const std::string &filename, int idx=-1);
 
        // Input                D       : Dipolar interaction
        //                      filename: Output file name
        //                      idx     : Interaction index (default -1->none)
	//			warn    : Warning output label
	//				   0 = no warnings
	//				   1 = warnings
	//				  >1 = fatal warnings
        // Output               none    : Dipolar interaction is read in
        //                                from parameters in file filename
        // Note                         : Since no spin indices are used only
        //                                a few parameters are acceptable
 

//int read(const ParameterSet& pset, int idx=-1, int warn=2);
//void read(const ParameterSet& pset, int idx=-1);
 
        // Input                D       : Dipolar interaction
        //                      pset    : Parameter set
        //                      idx     : Interaction index (default -1->none)
	//			warn    : Warning output label
	//				   0 = no warnings
	//				   1 = warnings
	//				  >1 = fatal warnings
        // Output               none    : Dipolar interaction is read in
        //                                from parameters in pset
        // Note                         : Since no spin indices are used only
        //                                a few parameters are acceptable
 
/* These next two read functions will try an read the parameter set
 
                      { DCC, Dtheta, Dphi, Deta }
 
   Since there no spin indices given, parameters Iso and Coord are not
   allowed.  Furthermore, the spin quantum numbers of the two spins
   involded in the interaction are taken to be 1/2, i.e. I=S=1/2. These
   functions do allow for the suffix (#) and the prefix [#], but they are
   not used in the default function call.  If Deta is not present
   (and this is normally the case) then Deta is assumed to be zero.          */

 
MSVCDLL void read(const std::string &filename, double Iqn, double Sqn, int idx=-1);
 
        // Input                D      : Dipolar interaction
        //                      filename: Output file name
        //                      Iqn     : Quantum I value (0.5, 1, ...)
        //                      Sqn     : Quantum S value (0.5, 1, ...)
        //                      idx     : Interaction index (default -1->none)
        // Output               none    : Dipolar interaction is read in
        //                                from parameters in file filename


MSVCDLL void read(const ParameterSet& pset, double Iqn, double Sqn, int idx=-1);
 
        // Input                D       : Dipolar interaction
        //                      pset    : Parameter set
        //                      Iqn     : Quantum I value (0.5, 1, ...)
        //                      Sqn     : Quantum S value (0.5, 1, ...)
        //                      idx     : Interaction index (default -1->none)
        // Output               none    : Dipolar interaction is read in
        //                                from parameters in pset

// ----------------------------------------------------------------------------
//        Interactive Ask/Read ASCII File Using A Single Interaction Index
// ----------------------------------------------------------------------------
 
MSVCDLL std::string ask_read(int argc, char* argv[], int argn,                    int idx=-1);
MSVCDLL std::string ask_read(int argc, char* argv[], int argn, const std::string& def, int idx=-1);
MSVCDLL void askI(int   argc, char* argv[], int qn, double& DI);
MSVCDLL void askS(int   argc, char* argv[], int qn, double& DS);
MSVCDLL void askDCC(int argc, char* argv[], int qn, double& Dcc);
 
        // Input                D       : Dipolar interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
        //                      idx     : Interaction index
        // Output               String  : The parameter argn of array argc is
        //                                used to supply a filename from which
        //                                the dipolar interaction is read
 
// ----------------------------------------------------------------------------
//          Interactive Ask/Read ASCII File Using Two Spin Indicies
// ----------------------------------------------------------------------------
 
MSVCDLL std::string ask_read(int argc, char* argv[], int argn, int iI,int iS);

        // Input                D       : Dipolar interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
        //                      iI	: Index for spin I
        //                      iS 	: Index for spin S
        // Output               String  : The parameter argn of array argc is
        //                                used to supply a filename from which
        //                                the dipolar interaction is read
        //                                If argument argn is not in argv, the
        //                                user is asked to supply a filename
        //                                The set filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing known IntDip parameters
        // Note                         : The interaction D is modifed (filled)


MSVCDLL std::string ask_read(int argc, char* argv[], int argn,
                                        double Iqn, double Sqn, int idx=-1);

        // Input                D       : Dipolar interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
        //                      Iqn     : Quantum I value (0.5, 1, ...)
        //                      Sqn     : Quantum S value (0.5, 1, ...)
        //                      idx     : Interaction index (default -1->none)
        // Output               String  : The parameter argn of array argc is
        //                                used to supply a filename from which
        //                                the dipolar interaction is read
        //                                If argument argn is not in argv, the
        //                                user is asked to supply a filename
        //                                The set filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing known IntDip parameters
        // Note                         : The interaction D is modifed (filled)

 
// ----------------------------------------------------------------------------
//               Interactive Ask For All Kinds Interaction Info
// ----------------------------------------------------------------------------
 
MSVCDLL void ask(int argc, char* argv[], int& argq, double& Iqn, double& Sqn,
                                     double& Cnqcc, double& Ceta, int Cflag=0);

        // Input                D       : Dipolar interaction (this)
        //                      argc    : Number of arguments 
        //                      argv    : Array of arguments 
        //                      argq 	: Query value 
        //                      Iqn	: Spin I quantum number
        //                      Sqn	: Spin S quantum number
        //                      Cnqcc   : D. coupling constant (Hz)
        //                      Ceta    : D. asymmetry 
        //                      theta	: D. orientation angle
        //                      Cphi    : D. orientation angle
	//			Cflag   : Flag is CCC or wDrequested
        // Output               none    : The values of qn, I, Cnqcc, Ceta,
        //                                theta, and Cphi are set herein 
        // Note                         : This is INTERACTIVE! 

 
MSVCDLL void askset(int argc, char* argv[], int& qn, int Cflag=0);
 
        // Input                D      : Dipolar interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Array of arguments
        //                      qn      : Query value
	//			Cflag   : Flag is CCC or wDrequested
        // Output               none    : Dis set interactively
        // Note                         : This is INTERACTIVE!
 

MSVCDLL void askset(int Cflag=0);

        // Input                D      : Dipolar interaction (this)
        // Output               none    : Dis set interactively
	//			Cflag   : Flag is CCC or wDrequested
        // Note                         : This is INTERACTIVE!

// ____________________________________________________________________________
// I                          OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//  Functions That Generate String Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

MSVCDLL std::vector<std::string> CartAStrings(const std::string& CSForm) const;
MSVCDLL std::vector<std::string> InfoStrings()                           const;
 
MSVCDLL std::vector<std::string> DipAStrings() const;
 
        // Input                D       : Dipolar interaction (this)
        // Output               CSS     : Pointer to array of 5 strings
        // Note                         : The String array must be deleted
        //                                outside of this routine!
 
/*               Spin Quantum Numbers:      I   , S
                 Coupling Constant:     xxxxx.xx xHz
                 Asymmetry:                 x.xx
                 Down From PAS z-Axis:    xxx.xx Degrees
                 Over From PAS x-Axis:    xxx.xx Degrees                     */  


MSVCDLL std::string* DipTStrings(int M) const;
 
        // Input                D       : Dipolar interaction (this)
        //                      M       : Ang. momentum component [0,4]
        // Output               TSS     : Pointer to array of hs strings
        //                                where hs is the spin pair
        //                                Hilbert space dimension
        // Note                         : The String array must be deleted
        //                                outside of this routine!
 
/*                              [ x.x, x.x, x.x]
                        T     = [ x.x, x.x, x.x]
                         2,m    [ x.x, x.x, x.x]
                                [ x.x, x.x, x.x]
 
  where M = { 0, 1, ..., 4 } ===> m = { 0, 1, -1, 2, -2 }                    */


//-----------------------------------------------------------------------------
//     Functions That Generate Ouput Of The Rank 2 Dipolar Interaction
//-----------------------------------------------------------------------------

MSVCDLL std::ostream& print(std::ostream& out, bool fflag=false, bool hdr=true) const;

        // Input                D      : Dipolar interaction (this)
        //                      ostr    : Output stream
        //                      fflag   : Output flag, true=full output
        // Output               none    : D spatial tensor parameters
        //                                placed into the output stream

/*               Spin Quantum Numbers:      I   , S
                 Coupling Constant:     xxxxx.xx xHz
                 Asymmetry:                 x.xx
                 Down From PAS z-Axis:    xxx.xx Degrees
                 Over From PAS x-Axis:    xxx.xx Degrees                     */

 
MSVCDLL std::ostream& printAT(std::ostream& ostr) const;

        // Input                D       : Dipolar spatial tensor (this)
        //                      ostr    : Output stream
        // Output               none    : Dipolar interaction spatial tensor
        //                                sent to the output stream
        // Note                         : Uses base class virtual overload


/*           Print The Spherical Components, Spatial And Spin
       These Will Be Printed Lines Which Will Look Like The Following
                      (Repeated For All 5 m Values)

                                                [ x.x, x.x, x.x]
                A    = x.xxx            T     = [ x.x, x.x, x.x]
                 2,m                     2,m    [ x.x, x.x, x.x]
                                                [ x.x, x.x, x.x]             */


MSVCDLL friend std::ostream& operator<< (std::ostream& out, const IntDip& D);

        // Input                out	: Output stream;
        //			D	: Dipolar tensor to write
        // Output			: Modifies output stream


// ____________________________________________________________________________
// L                    DIPOLE-DIPOLE HAMILTONIAN FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//               First Order Dipolar Interaction Hamiltonian
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
//  These Are For When The Dipolar Interaction Is A Perturbation To Zeeman
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

           Input                D       : Dipolar interaction
                                alpha   : Euler angle (radians)
                                beta    : Euler angle (radians)
                                gamma   : Euler angle (radians)
                                EA      : Euler angles (radians)
                                HSs     : Array of spin Hilbert spaces
                                i       : First spin index (in HSs)
                                j       : Seond spin index (in HSs)
                                wh      : Flag for weak heteronuclear
           Output               H0      : The secular part of the dipolar
                                          Hamiltonian (default basis, Hz)
           Note                         : Also called the 1st order dipolar
                                          interaction (perturbation theory)
           Note                         : Rotationally invariant about z
           Note                         : No HSs: return in the spin Hilbert
                                          space of dimension (2I+1)[*(2S+1)]
                                          With HSs: return in composite spin
                                          Hilbert space                      */

MSVCDLL matrix H0(                              bool wh=false) const;
MSVCDLL matrix H0(double A, double B, double G, bool wh=false) const;
MSVCDLL matrix H0(const EAngles& EA,            bool wh=false) const;

MSVCDLL matrix H0(const std::vector<int>& HSs, int i, int j, bool wh=false) const;
MSVCDLL matrix H0(const std::vector<int>& HSs, int i, int j,
                          double A, double B, double G, bool wh=false) const;
MSVCDLL matrix H0(const std::vector<int>& HSs, int i, int j,
                                     const EAngles& EA, bool wh=false) const;

// ----------------------------------------------------------------------------
//                 Full Dipolar Interaction Hamiltonians
//       These Have No Approximations & Are Independent of Field Strength
// ----------------------------------------------------------------------------

/* Here we make no approximations. As such, these Hamiltonians are not those
   used in the rotating frame. Here we simply sum over the 5 spatial and spin
   tensor components to produce the Hamiltonian. 

                        -2,2
                         ---   D       m    D                        D
 H (alpha,beta,gamma) =  \   Xi  * (-1)  * A   (alpha,beta,gamma) * T    (i,j)
  D                      /                  2,m                      2,-m
                         ---
                          m

           Input                D       : Dipolar interaction
				alpha	: Euler angle (radians)
				beta    : Euler angle (radians)
				gamma   : Euler angle (radians)
			 	EA	: Euler angles (radians)
				HSs	: Array of spin Hilbert spaces
				i	: First spin index (in HSs)
				j	: Seond spin index (in HSs)
	   Output		H	: Matrix for dipolar Hamiltonian
           Note                         : No HSs: return in the spin Hilbert
                                          space of dimension (2I+1)*(2S+1)
 					  With HSs: return in composite spin
					  Hilbert space                      */
 
// matrix IntRank2::H( ) const                                        INHERITED
// matrix IntRank2::H(double alpha, double beta, double gamma) const  INHERITED
// matrix IntRank2::H(const EAngles& EA) const                        INHERITED

// matrix IntRank2::H(const vector<int>& HSs, int i, int j) const     INHERITED
// matrix IntRank2::H(const vector<int>& HSs, int i, int j,           INHERITED
//                     double alpha, double beta, double gamma) const
// matrix IntRank2::H(const vector<int>& HSs, int i, int j,           INHERITED
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


MSVCDLL std::ostream& printSpherical(std::ostream& ostr);
MSVCDLL std::ostream& printCartesian(std::ostream& ostr);
MSVCDLL std::ostream& printCartesian(std::ostream& ostr, double theta, double phi=0);
         
// Axx = (A22+A2m2)/2 - A20/sqrt(6)      Axy = -i*(A22-A21)/2 = Ayx
// Ayy = -(A22+A2m2)/2 - A20/sqrt(6)     Axz = -(A21-A2m1)/2  = Azx
// Azz = sqrt(2/3)*A20                   Ayz = i*(A21-A2m1)/2 = Azy

        // Input                D	: D spatial tensor (this)
        //                      ostr	: Output stream
        //                      phi     : Orientation angle (degrees)
        //                      theta   : Orientation angle (degrees)
        // Output               none	: D spatial tensor parameters set to
        //                                output stream




// ____________________________________________________________________________
// M                 DIPOLE-DIPOLE HAMILTONIAN FRIEND FUNCTIONS
// ____________________________________________________________________________

// These Don't Depend on Class IntDip Variables At All.  The Functions Just
// Mimic The Class' Hamiltonian Generation.


MSVCDLL friend matrix HD0(double qn, double wCo, double eta, 
                                                 double theta, double phi);
 
        // Input                qn      : Quantum number (1, 1.5, 2.5,...)
        //                      wCo     : PAS Dipolar frequency
        //                      eta     : Dipolar asymmetry [0,1]
        //                      theta   : Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               H0      : The secular part of the Dipolar
        //                                Hamiltonian (default basis, Hz)
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note                         : Also called the 1st order Dipolar
        //                                interaction (perturbation theory)
        // Note                         : Rotationally invariant about z
        // Note                         : This will return in the spin Hilbert
        //                                space of dimension 2I+1

//  The secular part of the Dipolar Hamiltonian is that which returns
//  only those components which commute with z axis rotations.  Here we have

//            [1]             1                   2             (0)          
//           H  (theta,phi) = - w (the,phi) * [3*I - I(I+1)] = H
//            D              6  D               z             C

// where                                                          
 
//                      [ 1      2               1        2                   ]
// w (theta,phi) = W    | - [3cos (theta) - 1] + - eta sin (theta)*cos(2*phi) |
//  D               C,o [ 2                      2                            ]
 
// and
//                                    3*CCC
//                            w    = --------
//                             C,o   2I(2I-1)



MSVCDLL friend matrix HD1(double Om, double qn, double wCo, double eta,
                                             double theta, double phi);
                                                        
        // Input                Om      : Field Strength (Larmor in Hz)
	//                      qn      : Quantum number (1, 1.5, 2.5,...)
        //                      wCo     : PAS Dipolar frequency
        //                      eta     : Dipolar asymmetry [0,1]
        //                      theta   : Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HD1     : The 2nd order secular part of the
        //                                Dipolar Hamiltonian (default
        //                                basis, Hz) for the interaction
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note                         : Also called the 1st order Dipolar
        //                                interaction (perturbation theory)
        // Note                         : Rotationally invariant about z
        // Note                         : This will return in the spin Hilbert
        //                                space of dimension 2I+1


  };

#endif						// IntDip.h
