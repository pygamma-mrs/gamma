/* IntHF.h ******************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Hyperfine Interaction 		  	   Interface		**
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
** for a particular electron-nucleon spin pair.  The interaction is 	**
** maintained in the form of spatial and spin tensors which are stored 	**
** in irreducible spherical format.					**
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
**  I        - First spin quantum number                                **
**  S        - Second spin quantum number                               **
**  DELZZ    - Rank 2 spatial tensor delz value (Gauss).		**
**  AISO     - Isotropic hyperfine coupling (Gauss).			**
**  XI       - Hyperfine interaction constant (Hz)			**
**                                                                      **
** The tensor will internally be stored in irreducible spherical form   **
** at any chosen orientation.  The tensor principal axes will also be   **
** maintiained as well as a set of Euler angles that indicate the PAS   **
** orientation relative to some common axes.                            **
**                                                                      **
** In addition to the functions contained in the base spatial tensor    **
** class, this IntHF provides functions for building up the tensor   	**
** accessing the tensor elements from a hyperfine standpoint.  This     **
** includes functions for reorienting the tensor, making Hamiltonians,  **
** and obtaining electron resonance conditions.                         **
**                                                                      **
** The following defintions are used herein (Auv normlized by delzz):	**
**                                                                      **
** 1.) PAS: |Azz| >= |Ayy| >= |Axx|					**
** 2.) PAS: Azz=1, eta=(Ayy-AXX)/2, Axx=(1+eta)/2, Ayy=(eta-1)/2        **
**                                                                      **
** Individual spin tensors are stored in a linked list, handled by the  **
** class IntSTLList.  This prevents recalculation of the 5 arrays that  **
** are associated with such tensors when the spin(s) involved share the **
** same Iz & Sz values.  However, the arrays are still copied into new  **
** (equivalent) spin tensors, but that is done by reference in the      **
** GAMMA matrix class.                                                  **
**                                                                      **
*************************************************************************/

#ifndef   IntHF_h_			// Is file already included?
#  define IntHF_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface
#  endif

class complex;				// Know about complex numbers
class matrix;				// Know about matrices
class row_vector;			// Know about row vectors
//class ostream;				// Know about output streams
class spin_sys;				// Know about basic spin systems
class coord;				// Know about coordinates
class IntRank2A;			// Know about base class

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Basics/Isotope.h>		// Include spin isotopes
#include <IntRank2/IntRank2.h>		// Include base class headers
#include <IntRank2/IntRank2A.h>		// Include base class headers
#include <IntRank2/IntRank2T.h>		// Include base class headers
#include <string>			// Include libstdc++ strings
#include <vector>			// Include libstdc++ STL vectors

class IntHF: public IntRank2
  {
  double AISO;				// Hyperfine isotropic value (Gauss)
  double DELZZ;				// Hyperfine anisotropic value (Gauss)
  matrix T20wh;                         // Spherical T20 weak heteronuclear

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                    HYPERFINE INTERACTION ERROR HANDLING
// ____________________________________________________________________________

/*       Input 		      HF      : Hyperfine interaction (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

         void IHFerror(int eidx,                           int noret=0) const;
         void IHFerror(int eidx, const std::string& pname, int noret=0) const;
volatile void IHFfatal(int eidx)                                        const;
volatile void IHFfatal(int eidx, const std::string& pname)              const;

// ____________________________________________________________________________
// ii            HYPERFINE INTERACTION SETUP FUNCTIONS
// ____________________________________________________________________________

/* These by-pass the generalized GAMMA spatial-spin tensors (in class spin_T)
   and just produces the rank 2 space-spin tensor used for the hyperfine
   interaction directly.  They are scaled such that they are interaction
   independent and adhere to the relative scaling between components of true
   spherical tensors.  The rank 2 spatial tensors (used by this class) are
   scaled such that, when rotated they are normalized rank 2 spherical
   harmonics in a symmetric case (eta=0).                                    */

// void setAs();         INHERITED  Set spatial tensor comps (IntRank2)
// void setTs();         INHERITED  Set spin    tensor comps (IntRank2)


//void setTs();

        // Input                HF	: Hyperfine interaction (this)
        // Output               none    : Hyperfine interaction spherical
        //                                spin components are generated
        // Note                         : No check is made to see if the
        //                                Tsph array has be made!

/*                                            +
                                           m  |
                                T    = (-1)  T
                                 2,m          2,-m

   T  = Tsph[0]; T  = Tsph[1]; T   = Tsph[2]; T  = Tsph[3]; T   = Tsph[4]
    2,0           2,1           2,-1           2,2           2,-2

                     1/2
                  [4]                      1
            T   = |-| * I         T    = - - I           T    = 0
             2,0  [6]    z         2,1     2  +           2,2                */


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
   hyperfine interaction in GAMMA.                                          */

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
			hfc     : Isotropic hyperfine coupling (Gauss)
			hfa     : Hyperfine anisotropy value (Gauss)
			et a    : Hyperfine asymmetry value [0,1]
                        EA      : Euler angles for orientation (radians)
			idxI    : 1st spin or interaction index
                        idxS    : 2nd spin index
                        warn    : Warning level
                                    0 - no warnings
                                    1 - warnings
                                   >1 - fatal warnings
       Output           TF      : True if hyperfine interaction defined
                                  Argument values are set		     */

bool getHFI(const ParameterSet& pset, double& Iqn, double& Sqn,
         double& hfc, double& hfa, double& eta, EAngles& EA,
                                     int idxI, int idxS, bool warn=true) const;

bool getHFI1(const ParameterSet& pset, double& Iqn, double& Sqn,
         double& hfc, double& hfa, double& eta, EAngles& EA,
                                               int idxI, bool warn=true) const;

bool getHFI2(const ParameterSet& pset, double& Iqn, double& Sqn,
         double& hfc, double& hfa, double& eta, EAngles& EA,
                                     int idxI, int idxS, bool warn=true) const;


// ----------------------------------------------------------------------------
//         Get Isotropic Hypefine Coupling Value From A Parameter Set
// ----------------------------------------------------------------------------

/*      Input           HFI     : Hyperfine interaction (this)
                        pset    : A parameter set
                        idxI    : 1st spin or interaction index
                        idxS    : 2nd spin index
                        hfc     : Hypefine coupling (Gauss)
                        warn    : Warning output flag
        Output          TF      : True if hyperfine coupling obtained
                                  constant from parameters in pset
        Note                    : This WILL NOT alter the interaction
        Note                    : Parameters are A, A(#), A(#,#)             */

bool getHFC(const ParameterSet& pset, double& hfc,
                                  int idxI, int idxS=-1, bool warn=true) const;


// ----------------------------------------------------------------------------
//            Get Hypefine Anisotropy Value From A Parameter Set
// ----------------------------------------------------------------------------

/*      Input           HFI     : Hyperfine interaction (this)
                        pset    : A parameter set
                        idxI    : 1st spin or interaction index
                        idxS    : 2nd spin index
                        hfa     : Hypefine anisotropy (Gauss)
                        warn    : Warning output flag
        Output          TF      : True if hyperfine coupling obtained
                                  constant from parameters in pset
        Note                    : This WILL NOT alter the interaction
        Note                    : Parameters are AA, AA(#), AA(#,#)         */

bool getHFA(const ParameterSet& pset, double& hfa,
                                  int idxI, int idxS=-1, bool warn=true) const;

// ----------------------------------------------------------------------------
//                        Complete Hyperfile Interaction
// ----------------------------------------------------------------------------

/* This function employs all of the "get*" functions in the above sections
   to parse a parameter set for all of the values needed to define a hyperfine
   interaction, namely { Iqn,Sqn,A,Adelz,eta,alpha,beta,gamma }. If the
   interaction definition is found, we set the interaction or return false.  */

bool setHFI(const ParameterSet& ps,int idxI=-1,int idxS=-1,int warn=1);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A            HYPERFINE INTERACTION CONSTRUCTION, DESTRUCTION
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

MSVCDLC IntHF();
MSVCDLC IntHF(const IntHF &HF1);

// ----------------------------------------------------------------------------
//              Direct Constructors Using Cartesian Components
// ----------------------------------------------------------------------------

/*      Input      HF      : Hyperfine interaction (this)
                   II      : Spin I isotope type
                   IS      : Spin S isotope type
                   Iz      : Spin I quantum number
                   Sz      : Spin S quantum number
                   Axyz    : PAS HF Cartesian cmpts (Gauss)
		   EA      : Euler angles for orientation (radians)
        Output     none    : Hyperfine interaction constructed
                             from Axx, Ayy, Azz contained in
                             the input coordinate (in Gauss)                 */

MSVCDLC IntHF(const std::string& II, const std::string& IS,
                                  const coord& Axyz, const EAngles& EA=EAzero);
MSVCDLC IntHF(const Isotope&     II, const Isotope&     IS,
                                  const coord& Axyz, const EAngles& EA=EAzero);
MSVCDLC IntHF(      double       Iz,       double       Sz,
                                  const coord& Axyz, const EAngles& EA=EAzero);

// ----------------------------------------------------------------------------
//               Direct Constructors Using Spherical Components
// ----------------------------------------------------------------------------

/*      Input       HF	    : Hyperfine interaction (this)
                    IsoI    : Spin I isotope type
                    IsoS    : Spin S isotope type
        or          Iqn	    : Spin I quantum number
		    Sqn	    : Spin S quantum number
                    Aiso    : Spatial tensor Aiso value (Gauss)
                    Adelzz  : Spatial tensor delzz value (Gauss)
                    Aeta    : Spatial tensor asymmetry value [0,1]
        Output      none    : Hyperfine interaction constructed              */

MSVCDLC IntHF(const std::string& II, const std::string& IS,
        double Aiso, double Adelzz, double Aeta=0.0, const EAngles& EA=EAzero);
MSVCDLC IntHF(const Isotope&     II, const Isotope&     IS,
        double Aiso, double Adelzz, double Aeta=0.0, const EAngles& EA=EAzero);
MSVCDLC IntHF(      double       Iz,       double       Sz,
        double Aiso, double Adelzz, double Aeta=0.0, const EAngles& EA=EAzero);

// ----------------------------------------------------------------------------
//                    Construction Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions will construct the interaction from parameters in a
   specified GAMMA parameter set.  This is most useful when working with multi-
   spin systems that are setting up many interactions over the system using a
   single parameter file (or parameters read from an external ASCII file.

        Input                HF	     : Hyperfine interaction
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

MSVCDLC IntHF(const ParameterSet& pset, int idxI, int idxS=-1, int warn=2);

/*
IntHF(const ParameterSet& pset, int idxI, int idxS, int warn=2);
IntHF(int idxI, int idxS, const ParameterSet& pset, int warn=2);

IntHF(const std::string&  II, const std::string&  IS,
                    const ParameterSet& pset, int idxI, int IdxS, int warn=2);
IntHF(const  Isotope&     II, const Isotope&      IS,
                    const ParameterSet& pset, int idxI, int IdxS, int warn=2);
*/

// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

MSVCDLL void operator= (const IntHF &HF1);
MSVCDLC      ~IntHF();

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

MSVCDLL double iso() const;
MSVCDLL void   iso(double aiso);
MSVCDLL double A() const;
MSVCDLL void   A(double aiso);

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

MSVCDLL double aniso() const;
MSVCDLL void   aniso(double aiso);
MSVCDLL double AA()    const;
MSVCDLL void   AA(double aiso);

// ----------------------------------------------------------------------------
//                      Spatial Tensor Asymmetry Value
// ----------------------------------------------------------------------------

// Note that eta spans [0, 1] and is defined to be (Axx-Ayy)/Azz where
// |Azz| >= |Ayy| >= |Axx| in the treatment contained herein.  These functions
// are inherited from the base class IntRank2A

//MSVCDLL double eta( ) const;            INHERITED       Get the asymmetry
//void   eta(double HFeta);       INHERITED       Set the asymmetry

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
//                       Cartesian Tensor Component Access
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
   the current orientation is used.

                                    1/2
                            [ 6*PI ]
                      h   = | ---- |    * del   * A   + h
                       uv   [  5   ]         zz    uv    iso                 */

MSVCDLL double hxx() const;
MSVCDLL double hyy() const;
MSVCDLL double hzz() const;
MSVCDLL double hxy() const;
MSVCDLL double hyx() const;
MSVCDLL double hxz() const;
MSVCDLL double hzx() const;
MSVCDLL double hyz() const;
MSVCDLL double hzy() const;

MSVCDLL double hxx(double alpha, double beta, double gamma) const;
MSVCDLL double hyy(double alpha, double beta, double gamma) const;
MSVCDLL double hzz(double alpha, double beta, double gamma) const;
MSVCDLL double hyx(double alpha, double beta, double gamma) const;
MSVCDLL double hxy(double alpha, double beta, double gamma) const;
MSVCDLL double hzx(double alpha, double beta, double gamma) const;
MSVCDLL double hzy(double alpha, double beta, double gamma) const;
MSVCDLL double hxz(double alpha, double beta, double gamma) const;
MSVCDLL double hyz(double alpha, double beta, double gamma) const;

MSVCDLL double hxx(const EAngles& EA) const;
MSVCDLL double hyy(const EAngles& EA) const;
MSVCDLL double hzz(const EAngles& EA) const;
MSVCDLL double hxy(const EAngles& EA) const;
MSVCDLL double hyx(const EAngles& EA) const;
MSVCDLL double hxz(const EAngles& EA) const;
MSVCDLL double hzx(const EAngles& EA) const;
MSVCDLL double hyz(const EAngles& EA) const;
MSVCDLL double hzy(const EAngles& EA) const;

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

        1.) GAMMA normalized Auv - Done With CartMx()
        2.) Typical Suv values   - Done With Smx(true);
        3.) Shown in lab frame   - Done With This Function

   For case 3.) the values are related to the GAMMA normalized (Auv) and
   typically presented values (guv) according to

                            1/2
                    [ 6*PI ]
             HF   = | ---- |    * del   * A   + Kdel    * HF
               uv   [  5   ]         zz    uv       u,v     iso

   where Kdel is a Kronecker delta function.                                 */

MSVCDLL matrix Amx() const;

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
//          Spin Tensor Component Access (Spin Pair Hilbert Space)
//-----------------------------------------------------------------------------

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

MSVCDLL matrix T20het() const;
MSVCDLL matrix T20het(const std::vector<int>& HSs, int i, int j) const;

// ____________________________________________________________________________
// D             HYPERFINE INTERACTION CONSTANT ACCESS FUNCTIONS
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


MSVCDLL double xi( ) const;

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

/* IntHF has no implicit knowledge of the spins involved in the interaction,
   other than their spin quantum numbers. The preferred means of specifying a
   hyperfine interaction when there is no knowledge of the spin types types
   is that which is used for filling up the parameter set to be written. The
   base parameters of individual interactions in that case are

            { AI(#), AS(#), Axx(#), Ayy(#), Acc(#), AEAngles(#)}.

           Input                HF      : Hyperfine interaction
                                idx     : Interaction index (default -1)
                                pfx     : Interaction 2nd indx (def -1)
           Output               pset    : Parameter set with only
                                          electron G interaction parameters

     Function                                 Purpose
   ------------         -------------------------------------------------------
   ParameterSet         Convert interaction into a parameter set
   operator +=          Adds interaction to existing parameter set
   PSetAdd              Adds interaction to existing parameter set, this
                        allows for an index and a prefix in the parameters   */


MSVCDLL              operator ParameterSet( ) const;
MSVCDLL friend void operator+= (ParameterSet& pset, const IntHF &HF);
MSVCDLL void PSetAdd(ParameterSet& pset, int idx=-1, int pfx=-1) const;

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
           Output               none    : Hyperfine interaction parameters
                                          written in parameter set format to
                                          file filename or ostream ofstr     */

MSVCDLL int write(const std::string &fn,int idx=-1,int pfx=-1,int wrn=2) const;
MSVCDLL int write(std::ofstream& ofstr, int idx=-1,int pfx=-1,int wrn=2) const;

// ____________________________________________________________________________
// G                   HYPERFINE INTERACTION INPUT FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                    Functions To Read From An ASCII File
// ----------------------------------------------------------------------------

/* These next two read functions utilize either two spin indices or a single
   interaction index. They'll try to read the hyperfine interaction form
   parameters found either in an ASCII file or in a GAMMA parameter set. These
   functions do NOT allow for the prefix [#] because multiple hyperfine
   interactions can be defined in the same file by switching either the spin
   pair indices or the interaction index.  Multiple sets of interactions can be
   read using hyperfine interaction vectors (see class IntHFVec.)

        // Input                HF      : Hyperfine interaction
        //                      filename: Output file name
        //                      pset    : Parameter set
        //                      idxI    : Index for spin I or interaction
        //                      idxS    : Index for spin S (default -1)
        //                      warn    : Warning output label
        //                                 0 = no warnings
        //                                 1 = warnings
        //                                >1 = fatal warnings
        // Output               none    : Interaction read from parameters
        //                                in file or pset                    */

MSVCDLL bool read(const std::string &filename,  int idxI=-1,int IdxS=-1,int warn=2);
MSVCDLL bool read(const ParameterSet& pset,int idxI=-1,int idxS=-1,int warn=2);

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

           Input                HF	: Hyperfine interaction
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

MSVCDLL matrix H0(bool wh) const;
MSVCDLL matrix H0(double A, double B, double G, bool wh) const;
MSVCDLL matrix H0(const EAngles& EA, bool wh) const;


MSVCDLL matrix H0(const std::vector<int>& HSs, int i, int j, bool wh) const;
MSVCDLL matrix H0(const std::vector<int>& HSs, int i, int j,
                  double A, double B, double G, bool wh) const;
MSVCDLL matrix H0(const std::vector<int>& HSs, int i, int j,
                             const EAngles& EA, bool wh) const;

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
                             m

           Input                HF	: Hyperfine interaction
                                alpha   : Euler angle (radians)
                                beta    : Euler angle (radians)
                                gamma   : Euler angle (radians)
                                EA      : Euler angles (radians)
                                HSs     : Array of spin Hilbert spaces
                                i       : First spin index (in HSs)
                                j       : Seond spin index (in HSs)
           Output               H       : Matrix for hyperfine Hamiltonian
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
// I                 HYPERFINE INTERACTION OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//  Functions That Generate String Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

/*  Function   Arguments                     Output
    ========   =========   ===================================================
    TStrings       m       String array for the mth component of spin tensor T
    GAStrings              String array for various interaction values

            GAStrings()                              TStrings(m)

     Hyperfine Coupling:   xxxx.xx Gauss          [ x.x, x.x, x.x]
     Hyperfine Anisotropy: xxxx.xx Gauss   T    = [ x.x, x.x, x.x]
     Hyperfine Asymmetry:  xxxx.xx Gauss    2,m   [ x.x, x.x, x.x]
     Down From PAS z-Axis:    x.xx Deg.
     Over From PAS x-Axis:    x.xx Deg.
     Electron I Value:        I               m = [0,4] => {0,1,-1,2,-2}
     Nucleus I Value:         S                                              */

//vector<string> CartAStrings(const std::string& CSForm) const;
//vector<string> SphAStrings()                           const;
MSVCDLL std::vector<std::string> CartAStrings(const std::string& CSForm) const;
MSVCDLL std::vector<std::string> InfoStrings()                           const;
MSVCDLL std::vector<std::string> SphAStrings()                           const;


//-----------------------------------------------------------------------------
//     Functions That Generate Ouput Of The Rank 2 Hyperfine Interactions
//-----------------------------------------------------------------------------

/* These functions will output information concerning the hyperfine interaction
   to any output stream.

           Input                HFI	: Hyperfine interaction (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Basic Parameters
                                           !0 - Full output
                                nrm     : Flag if GAMMA normalized output
           Output               none    : Hyperfine interaction information
                                          placed into the output stream      */

MSVCDLL        std::ostream& print(std::ostream& out, int fflag=-1) const;
MSVCDLL friend std::ostream& operator<<  (std::ostream& out, IntHF& HF);


MSVCDLL std::ostream& printSpherical(std::ostream& ostr);
MSVCDLL std::ostream& printCartesian(std::ostream& ostr);
MSVCDLL std::ostream& printCartesian(std::ostream& ostr, double theta, double phi=0);
         
// Axx = (A22+A2m2)/2 - A20/sqrt(6)      Axy = -i*(A22-A21)/2 = Ayx
// Ayy = -(A22+A2m2)/2 - A20/sqrt(6)     Axz = -(A21-A2m1)/2  = Azx
// Azz = sqrt(2/3)*A20                   Ayz = i*(A21-A2m1)/2 = Azy

        // Input                HF	: Hyperfine spatial tensor (this)
        //                      ostr	: Output stream
        //                      phi     : Orientation angle (degrees)
        //                      theta   : Orientation angle (degrees)
        // Output               none	: Hyperfine spatial tensor parameters set to
        //                                output stream













// ____________________________________________________________________________
//                   HYPERFINE HAMILTONIAN FRIEND FUNCTIONS
// ____________________________________________________________________________

// These Don't Depend on Class IntHF Variables At All.  The Functions Just
// Mimic The Class' Hamiltonian Generation.


//friend matrix HHF0(double qn, double wCo, double eta=0, 
//                                                 double theta=0, double phi=0);
 
        // Input                qn      : Quantum number (1, 1.5, 2.5,...)
        //                      wCo     : PAS Hyperfine frequency
        //                      eta     : Hyperfine asymmetry [0,1]
        //                      theta   : Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               H0      : The secular part of the Hyperfine
        //                                Hamiltonian (default basis, Hz)
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note                         : Also called the 1st order Hyperfine
        //                                interaction (perturbation theory)
        // Note                         : Rotationally invariant about z
        // Note                         : This will return in the spin Hilbert
        //                                space of dimension 2I+1

//  The secular part of the Hyperfine Hamiltonian is that which returns
//  only those components which commute with z axis rotations.  Here we have

//            [1]             1                   2             (0)          
//           H  (theta,phi) = - w (the,phi) * [3*I - I(I+1)] = H
//            HF              6  HF               z             HF

// where                                                          
 
//                      [ 1      2               1        2                   ]  
// w (theta,phi) = W    | - [3cos (theta) - 1] + - eta sin (theta)*cos(2*phi) |
//  D              C,o [ 2                      2                            ]  
 
// and
//                                    3*CCC
//                            w    = --------
//                             C,o   2I(2I-1)



//friend matrix HHF1(double Om, double qn, double wCo, double eta,
//                                              double theta=0.0, double phi=0.0);
                                                        
        // Input                Om      : Field Strength (Larmor in Hz)         
        //                      qn      : Quantum number (1, 1.5, 2.5,...)
        //                      wCo     : PAS Hyperfine frequency
        //                      eta     : Hyperfine asymmetry [0,1]
        //                      theta   : Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HD1     : The 2nd order secular part of the
        //                                Hyperfine Hamiltonian (default
        //                                basis, Hz) for the interaction
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note                         : Also called the 1st order Hyperfine
        //                                interaction (perturbation theory)
        // Note                         : Rotationally invariant about z
        // Note                         : This will return in the spin Hilbert
        //                                space of dimension 2I+1




//friend void operator+= (ParameterSet& pset, const IntHF &HF);
//       void PSetAdd(ParameterSet& pset, int idx=-1, int pfx=-1) const;

//void write(const string &filename, int idx=-1);

	// Input		AHF	: Hyperfine interaction (this)
	//			filename: Output file name
	//			idx	: Tensor index (default -1)
	// Output		none 	: Hyperfine interaction is
	//				  written as a parameter set to
	//				  file filename

 
// ____________________________________________________________________________
//                            STANDARD I/O FUNCTIONS
// ____________________________________________________________________________


//string ask_read(int argc, char* argv[], int argn, int idxI,int idxS);

        // Input                HF	: Hyperfine interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
        //                      idxI    : Index for spin I
        //                      idxS    : Index for spin S
        // Output               String  : The parameter argn of array argc is
        //                                used to supply a filename from which
        //                                the hyperfine interaction is read
        //                                If argument argn is not in argv, the
        //                                user is asked to supply a filename
        //                                The set filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing known IntHF parameters
        // Note                         : The interaction D is modifed (filled)

//string ask_read(int argc, char* argv[], int argn,
//                                   double Iqn=0.5, double Sqn=0.5, int idx=-1);

        // Input                HF	: Hyperfine interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vecotr of argc arguments
        //                      argn    : Argument index
        //                      Iqn     : Quantum I value (0.5, 1, ...)
        //                      Sqn     : Quantum S value (0.5, 1, ...)
        //                      idx     : Interaction index (default -1->none)
        // Output               String  : The parameter argn of array argc is
        //                                used to supply a filename from which
        //                                the hyperfine interaction is read
        //                                If argument argn is not in argv, the
        //                                user is asked to supply a filename
        //                                The set filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing known IntHF parameters
        // Note                         : The interaction D is modifed (filled)

 
//void ask(int argc, char* argv[], int& argq, double& Iqn, double& Sqn,
//       double& Cnqcc, double& Ceta, double& theta, double& Cphi, int Cflag=0);

        // Input                HF	: Hyperfine interaction (this)
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

 
//void askset(int argc, char* argv[], int& qn, int Cflag=0);
 
        // Input                D      : Hyperfine interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Array of arguments
        //                      qn      : Query value
	//			Cflag   : Flag is CCC or wDrequested
        // Output               none    : Dis set interactively
        // Note                         : This is INTERACTIVE!
 

//void askset(int Cflag=0);

        // Input                D      : Hyperfine interaction (this)
        // Output               none    : Dis set interactively
	//			Cflag   : Flag is CCC or wDrequested
        // Note                         : This is INTERACTIVE!
  };

#endif							// IntHF.h
