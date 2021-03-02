/* IntRank2.h ***************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Generic Rank 2 Interaction		    Interface		**
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
** This class is the base class for rank 2 interactions in GAMMA. First	**
** and foremost, it serves as a generic irreducible rank 2 tensor. It	**
** also contains linked lists of irreducible rank 2 spin tensors that	**
** are applicable to the interactions common to magnetic resonance. 	**
** These lists allow GAMMA to avoid repeat generation of spin tensors, 	**
** a likely occurance in the treatment of multi-spin systems because 	**
** the spin tensors depend only on the Iz values of spins involved and	**
** the interaction type (in the single spin or spin pair Hilbert space)	**
** Note that the use of linked lists doesn't avoid copying repeated	**
** spin tensors, but that is taken care of via referencing in GAMMA's	**
** matrix classes. 							**
**                                                                      **
*************************************************************************/

#ifndef   IntRank2_h_			// Is file already included?
#  define IntRank2_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Basics/Isotope.h>		// Include spin isotope types
#include <IntRank2/IntRank2T.h>		// Include interaction spin tensors
#include <IntRank2/IntRank2A.h>		// Include interaction space tensors
#include <list>				// Include stdlibc++ STL lists
#include <string>			// Include stdlibc++ strings

static std::list<IntRank2T> SPFlist;	// List of spin-field interact. tensors
static std::list<IntRank2T> SPSPlist;	// List of spin-spin  interact. tensors
static std::list<IntRank2T> SPQlist;	// List of spin-self  interact. tensors

class IntRank2: public IntRank2A, public IntRank2T
  {
  double  _XI;				// Interaction constant (1/sec)
//  static std::list<IntRank2T> SPFlist;	// List of spin-field interact. tensors
//  static std::list<IntRank2T> SPSPlist;	// List of spin-spin  interact. tensors
//  static std::list<IntRank2T> SPQlist;	// List of spin-self  interact. tensors

  friend class IntDip;			// Full access to dipolar interactions
  friend class IntCSA;			// Full access to shift anisot. intacs.
  friend class IntQuad;			// Full access to quadrupolar interact.
  friend class IntG;			// Full access to electron G interacts.
  friend class IntHF;			// Full access to hyperfine interacts.
 
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
                                                                                
         void IR2error(int eidx,                        int noret=0) const;
volatile void IR2fatal(int eidx)                                     const;
         void IR2error(int eidx, const std::string& pn, int noret=0) const;

// ____________________________________________________________________________
// ii               INTERACTION SPIN TENSOR SETUP FUNCTIONS
// ____________________________________________________________________________

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
           Note                 : For G interactions Ival is usually 2
                                  (Ival = 2*0.5 + 1, where I = 1/2 e-)
           Note                 : SA & G spin tensors are stored in the
                                  linked list SPFlist.  Thus, herein
                                  the spin tensor components are either
                                  constructed & stored, or copied from
                                  the spin tensor in the linked list.        */

void setSPF();

// ----------------------------------------------------------------------------
//            Set Up Single Spin Interaction Spin Tensor Components
// ----------------------------------------------------------------------------

/* Here we set the 5 spherical spin tensor components applicable to single spin
   interactions. This currently includes only quadrupolar interactions. However
   since this class maintains "scaled" spin tensors, all single spin
   interactions would share the same spin tensors. The five arrays will only
   differ between single spin interactions if the spin quantum value differs.
   This class maintains a static linked list of such spin tensors to avoid
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
          Note                         : QUAD spin tensors are stored in the
                                         linked list SPQlist.  Thus, herein
                                         the spin tensor components are either
                                         constructed & stored, or copied from
                                         the spin tensor in the linked list.
          Note                         : For QUAD interactions Ival is always
                                         greater than 2 (Ival=2*I+1, I>=1)   */

void setSPQ();
 
// ----------------------------------------------------------------------------
//            Set Up Spin Pair Interaction Spin Tensor Components
// ----------------------------------------------------------------------------

/* Here we set the 5 spherical spin tensor components applicable to spin-spin
   interactions. These currently include hyperfine (elecron-nucleon) and
   dipolar (electron-electron and nucleon-nucleon) interactions.  Since this
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
          Output               none    : Interaction Spin Tensor spherical
                                         spin components are generated
          Note                         : Hyperfine interactions will typically
                                         have either Ival and/or Sval of 2
                                         (Ival = 2*0.5 + 1, where I = 1/2 e-)
          Note                         : Hyperfine interactions will always
                                         have either Ival and/or Sval of 2
                                         (Ival = 2*0.5 + 1, where I = 1/2 e-)*/

void setSPSP();

// ____________________________________________________________________________
// iii            INTERACTION PARAMETER SET PARSING FUNCTIONS
// ____________________________________________________________________________

/* These functions glean information describing a generic rank 2 interaction
   from parameters in a GAMMA parameter set. Since base classes IntRank2A and
   IntRank2T handle irreducible rank 2 spatial and spin tensors respecitvely,
   we need only define the functions that deal with the reducible aspects
   of the rank 2 spatatial components.  This is done to support classes
   based on us! Even though we are irreducible, our derived interactions
   may utilize reducible rank 2 spatial tensors.                             */

// ----------------------------------------------------------------------------
//        Functions To Read Reducible Rank 2 Spherical Spatial Tensors
// ----------------------------------------------------------------------------

/* A reducible rank 2 spherical spatial tensor will be a contain 9 components.
   If we neglect the rank 1 component there are only 5 components. If
   the interaction is then represented in its PAS there are only 3 components
   needed to characterize the spatial tensor { Aiso, delzz, eta } We look for
   parameters that will yeild these three values. ALL THREE VALUES MUST BE
   SPECIFIED or we return that we failed. But in this case such a failure does
   not mean we did not obtain enough information to characterize the spatial
   tensor. For example, some interactions have no isotropic component so it
   would be silly to even try and set its value (quadrupolar & dipolar 
   interactions have no isotropic value). Furthermore, many interactions have
   either no asymmetry or their asymmetry is neglected.  Again, we do not want
   to force declaration of a zero asymmetry value when it is not needed. For
   example dipolar interactions normally have no asymmetry and it is rare that
   calculations involving hyperfine intereactons include its asymmetry. Finally,
   interactions that have an isotropic term may have negligible anisotropy and
   hence no delzz term. Again, we don't want lack of the anisotropy to produce
   any errors. So, some of these functions have additional warning flags so
   that derived classes can turn then on and off appropriately. Any values not
   found will be zeroed so that they may be checked after the funciton call.

   Currently allowed parameters are:           Aiso,   Aiso(#),   Aiso(#,#)
 (A is the base parameter name input)          Adelz,  Adelz(#),  Adelz(#,#)
                                               Adelzz, Adelzz(#), Adelzz(#,#)
   By definition the anisotropy is             AA,     AA(#),     AA(#,#)
   larger than the delzz value by 	       Aeta,   Aeta(#),   Aeta(#,#)
   a factor of 1.5.                                                          */

bool getAiAzAe(const ParameterSet& pset, const std::string& A,
                    coord& Aize,  int idxI, int idxS=-1, 
                     bool warni=true, bool warnz=true, bool warne=false) const;

bool getAiso(const ParameterSet& pset, const std::string& A,
                    double& Aiso, int idxI, int idxS=-1, bool warn=true) const;

bool getAaniso(const ParameterSet& pset, const std::string& A,
                  double& Aaniso, int idxI, int idxS=-1, bool warn=true) const;


// ----------------------------------------------------------------------------
//     Functions To Read External Field Strengths, Larmor/Base Frequencies
// ----------------------------------------------------------------------------

/* These may be required for setting up field dependent interactions such
   as shift anisotropy and electron G interactions. The field is returned in
   Gauss and the only two parameters are Field and FieldT. These can take
   a single index, but normally they should not have one.                    */


bool getField(const  ParameterSet& pset, const std::string& A, double& Bo,
                                             int idx=-1, bool warn=true) const;
bool getOmega(const  ParameterSet& pset, const std::string& A, double& Om,
                                             int idx=-1, bool warn=true) const;
bool getGOmega(const ParameterSet& pset, const std::string& A, double& Om,
                                             int idx=-1, bool warn=true) const;

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
public:
 
// ____________________________________________________________________________
// A               RANK 2 INTERACTION CONSTRUCTORS AND DESTRUCTOR
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

MSVCDLC IntRank2();
MSVCDLC IntRank2(const IntRank2 &R2I);

// ----------------------------------------------------------------------------
//   Direct Single Spin & Spin - Field Constructors Taking Lots 'O Parameters
// ----------------------------------------------------------------------------

/* We set up an interaction involving a single spin interacting with a static
   magnetic field or with its own local field. Two cases of the former are
   anisotropy (NMR) & g electron (EPR) interactions. The latter type would
   include quadrupolar interactions Since this class maintains "scaled" spatial
   & spin tensors, all single spin interactions will share the same spin
   tensors. The five arrays will only differ between single spin interactions
   if the spin quantum value differs (or if field orientation is set different
   in spin-field interactions). Typically, the static magnetic field is taken
   to be a normalized vector pointing along +z.  Similarly, all interations
   share the same formulation of spatial tensor components.  These will
   differ only when orientation differs.  Note that allowed interaction types
   are defined in class IntRank2T.  In this function the allowed types are
   currently G, CSA, and QUAD.

        Input                R2I     : Rank 2 interaction (this)
                             IsI     : Spin isotope type
                             dz      : Tensor delzz value (Hz)
                         or
                             qn      : Quantum number (e.g. 1.5, 0.5,...)
                             E       : Tensor asymmetry value
                             T       : Spin tensor type
                                        true  = Spin Field (SA, G)
                                        false = Spin Self (Quad)
        Output               none    : Interaction constructed               */


MSVCDLC IntRank2(const std::string& IsI, double chi,
                          double eta=0, const EAngles& EA=EAzero, bool T=true);
MSVCDLC IntRank2(const Isotope&     IsI, double chi,
                          double eta=0, const EAngles& EA=EAzero, bool T=true);
MSVCDLC IntRank2(const double       Isi, double chi,
                          double eta=0, const EAngles& EA=EAzero, bool T=true);

MSVCDLC IntRank2(const std::string& IsI, double chi,
                   const coord& AxAyAz, const EAngles& EA=EAzero, bool T=true);
MSVCDLC IntRank2(const Isotope&     IsI, double chi,
                   const coord& AxAyAz, const EAngles& EA=EAzero, bool T=true);
MSVCDLC IntRank2(double             IsI, double chi,
                   const coord& AxAyAz, const EAngles& EA=EAzero, bool T=true);

// ----------------------------------------------------------------------------
//          Direct Two Spin Constructors That Need Lots 'O Parameters
// ----------------------------------------------------------------------------

/* These constructors apply to spin-spin interactions (dipolar & hyperfine).
   Valid interaction types are specified in IntRank2T, these are currently
   HF and DIP. We limit the interaction types to spin pair interactions, other
   types will produce fatal errors.

        Input                R2I     : Rank 2 interaction (this)
                    -        IsoI    : Spin I isotope type
                   |         IsoS    : Spin S isotope type
                   |   or
                   |         Iqn     : Spin I quantum value (0.5, 1.5,...)
                   |_        Sqn     : Spin S quantum value
                             Xi      : Interaction strength (1/sec)
                    -        eta     : Spatial asymmetry value (default 0)
                   |   or
                   |_        AxAyAz  : Tensor PAS Cartesian components
                             EA      : Euler angles (orientation)            */


MSVCDLC IntRank2(const std::string&  IsoI, const std::string& IsoS,
                          double Xi, double eta=0, const EAngles& EA = EAzero);
MSVCDLC IntRank2(const Isotope&      IsoI, const Isotope&     IsoS,
                          double Xi, double eta=0, const EAngles& EA = EAzero);
MSVCDLC IntRank2(      double        Iqn,        double       Sqn,
                          double Xi, double eta=0, const EAngles& EA = EAzero);
MSVCDLC IntRank2(const std::string&  IsoI, const std::string&      IsoS,
                   double Xi, const coord& AxAyAz, const EAngles& EA = EAzero);
MSVCDLC IntRank2(const Isotope&      IsoI, const Isotope&     IsoS,
                   double Xi, const coord& AxAyAz, const EAngles& EA = EAzero);
MSVCDLC IntRank2(      double        Iqn,        double       Sqn,
                   double Xi, const coord& AxAyAz, const EAngles& EA = EAzero);


// ---------------------------------------------------------------------------- 
//               Here Be The Assignment Operator & Destructor
// ----------------------------------------------------------------------------

MSVCDLL void    operator= (const IntRank2 &R2I1);
MSVCDLC virtual ~IntRank2 ();

// ____________________________________________________________________________
// B                     SPATIAL TENSOR COMPONENT ACCESS
// ____________________________________________________________________________
 
//-----------------------------------------------------------------------------
//			       Anisotropy Access
//-----------------------------------------------------------------------------

/*                  ^
   The anisotropy, /_\ A, relates to the delzz value of A. Because GAMMA uses
   irreducible rank 2 spatial tensors as a base class, both the anisotropy and
   the delzz value are constant for all spatial tensors.

                                  1/2          
                             [ 5 ]                  ^      3
   Class IntRank2A:  del   = |---]  = 0.51503      /_\ A = - del   = 0.77255
                        zz   [6PI]                         2    zz

                                                    ^      
   Class IntRank2:   del   = Xi * 0.51503          /_\ A = Xi * 0.77255
                        zz                             

 
static double IntRank2A::delzz();				INHERITED
static double IntRank2A::delA();				INHERITED    */
 
//-----------------------------------------------------------------------------
//			       Asymmetry Access
//-----------------------------------------------------------------------------

/* Note that the asymmetry, eta, spans [0, 1] and is defined in GAMMA to be
   (Axx-Ayy)/Azz where |Azz| >= |Ayy| >= |Axx|. Although the definition is for
   normalized spatial components it will hold regardless of scaling.
 
   double IntRank2A::eta( ) const;				INHERITED
   void   IntRank2A::eta(double Eta);				INHERITED    */
 
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
   there is no asymmetry.						     */
 
/*
complex A20PAS(  ) const; 			// INHERITED
complex A21PAS(  ) const;			// INHERITED 
complex A2m1PAS( ) const; 			// INHERITED
complex A22PAS(  ) const;			// INHERITED
complex A2m2PAS( ) const;			// INHERITED
 
complex A20(  ) const; 			// INHERITED
complex A21(  ) const;			// INHERITED 
complex A2m1( ) const; 			// INHERITED
complex A22(  ) const;			// INHERITED
complex A2m2( ) const;			// INHERITED

complex IntRank2A::A20(double  alpha,double beta,double gamma) const; INHERITED
complex IntRank2A::A21(double  alpha,double beta,double gamma) const; INHERITED
complex IntRank2A::A2m1(double alpha,double beta,double gamma) const; INHERITED
complex IntRank2A::A22(double  alpha,double beta,double gamma) const; INHERITED
complex IntRank2A::A2m2(double alpha,double beta,double gamma) const; INHERITED

complex IntRank2A::A20(const  EAngles& EA) const;		INHERITED
complex IntRank2A::A21(const  EAngles& EA) const;		INHERITED
complex IntRank2A::A2m1(const EAngles& EA) const;		INHERITED
complex IntRank2A::A22(const  EAngles& EA) const;		INHERITED
complex IntRank2A::A2m2(const EAngles& EA) const;		INHERITED

complex IntRank2A::A2m(int m)                            const; INHERITED
complex IntRank2A::A2m(int m,double A,double B,double G) const; INHERITED
complex IntRank2A::A2m(int m, const EAngles& EA)         const; INHERITED

IR2ASph IntRank2A::SphCmp()                            const;   INHERITED
IR2ASph IntRank2A::SphCmp(double A,double B, double G) const;   INHERITED
IR2ASph IntRank2A::SphCmp(const EAngles& EA)           const;   INHERITED

complex IntRank2A::AcompPAS(int comp) const;			INHERITED 
complex IntRank2A::Acomp(int    comp) const;			INHERITED    */
  
//-----------------------------------------------------------------------------
//	       Cartesian Tensor Component Access, Normalized
//-----------------------------------------------------------------------------

/* These allow one to access the irreducible tensor rank 2 Cartesian elements
   at a specific orientation without rotating the entire tensor. In these
   functions theta is the angle down from the PAS z axis and phi the angle over
   from PAS x axis. Note that these all use GAMMA scaling which sets A2m to be 
   rank 2 spherical harmonics at all oreintations if there is no asymmetry. 
   Scaling on the Caretsian components is a direct result of this.

double IntRank2A::AxxPAS() const;		       		      INHERITED
double IntRank2A::AyyPAS() const;		       		      INHERITED
double IntRank2A::AzzPAS() const;		       		      INHERITED
double IntRank2A::AyxPAS() const;		       		      INHERITED
double IntRank2A::AxyPAS() const;		       		      INHERITED
double IntRank2A::AzxPAS() const;		       		      INHERITED
double IntRank2A::AxzPAS() const;		       		      INHERITED
double IntRank2A::AzyPAS() const;		       		      INHERITED
double IntRank2A::AyzPAS() const;		       		      INHERITED

double IntRank2A::Axx() const;		        		      INHERITED
double IntRank2A::Ayy() const;		        		      INHERITED
double IntRank2A::Azz() const;		        		      INHERITED
double IntRank2A::Ayx() const;		        		      INHERITED
double IntRank2A::Axy() const;		        		      INHERITED
double IntRank2A::Azx() const;		        		      INHERITED
double IntRank2A::Axz() const;		        		      INHERITED
double IntRank2A::Azy() const;		        		      INHERITED
double IntRank2A::Ayz() const;		        		      INHERITED

double IntRank2A::Axx(double alpha, double beta, double gamma) const; INHERITED
double IntRank2A::Ayy(double alpha, double beta, double gamma) const; INHERITED
double IntRank2A::Azz(double alpha, double beta, double gamma) const; INHERITED
double IntRank2A::Ayx(double alpha, double beta, double gamma) const; INHERITED
double IntRank2A::Axy(double alpha, double beta, double gamma) const; INHERITED
double IntRank2A::Azx(double alpha, double beta, double gamma) const; INHERITED
double IntRank2A::Axz(double alpha, double beta, double gamma) const; INHERITED
double IntRank2A::Azy(double alpha, double beta, double gamma) const; INHERITED
double IntRank2A::Ayz(double alpha, double beta, double gamma) const; INHERITED

double IntRank2A::Axx(const EAngles& EA) const;	                      INHERITED
double IntRank2A::Ayy(const EAngles& EA) const;	                      INHERITED
double IntRank2A::Azz(const EAngles& EA) const;	                      INHERITED
double IntRank2A::Ayx(const EAngles& EA) const;	                      INHERITED
double IntRank2A::Axy(const EAngles& EA) const;	                      INHERITED
double IntRank2A::Azx(const EAngles& EA) const;	                      INHERITED
double IntRank2A::Axz(const EAngles& EA) const;	                      INHERITED
double IntRank2A::Azy(const EAngles& EA) const;	                      INHERITED
double IntRank2A::Ayz(const EAngles& EA) const;	                      INHERITED

row_vector IntRank2A::CartComps() const;			   INHERITED
row_vector IntRank2A::CartComps(double A,double B,double G) const; INHERITED
row_vector IntRank2A::CartComps(const EAngles& EA)          const; INHERITED */

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

// static coord  AisoDelzEta(const coord& AxAyAz);
// static void   SortAxAyAz(double& x, double& y, double& z);
//        bool   CheckEta(double eta, bool warn=true) const;
//        int    Symmetric( ) const;
//        matrix CartMx(double scale=1.0) const;                  */
 
MSVCDLL bool   Spherical( )       const;
MSVCDLL bool   Isotropic( )       const;
MSVCDLL matrix CartMx(bool scale) const;

// ____________________________________________________________________________
// C                       SPIN TENSOR COMPONENT ACCESS
// ____________________________________________________________________________
 
/* These functions allow users to access values associated with the interaction
   spin tensor. Since this class is derived from the base class IntRank2T,
   all such functionality will be inherited from that class.                 */

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
  matrix IntRank2T::T20()           const;	// Return T2,0	INHERITED
  matrix IntRank2T::T21()           const;	// Return T2,1	INHERITED
  matrix IntRank2T::T2m1()          const;	// Return T2,-1	INHERITED
  matrix IntRank2T::T22()           const;	// Return T2,2	INHERITED
  matrix IntRank2T::T2m2()          const;	// Return T2,-2	INHERITED
  matrix IntRank2T::T2m(int m)      const;	// Return T2,m	INHERITED

  matrix IntRank2T::T20(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T21(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m1(     const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T22(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m2(     const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m(int m,const vector<int>& HSs,int i,int j=-1) const;  */

//-----------------------------------------------------------------------------
//                      Spin Tensor Output Functions 
//-----------------------------------------------------------------------------

/*String* IntRank2T::TStrings(int M) const;			INHERITED    */


// ____________________________________________________________________________
// D                     RANK 2 INTERACTION ACCESS
// ____________________________________________________________________________

        // Input                SA      : Shift anisotropy interaction
        // Output               delzz   : Spatial tensor delzz value
	//			xi      : Interaction constant
        // Note                         : Get/Set interaction delzz value
	//                                or interaction constant
 
//double delz() const;		// Get delzz value
//double delzz() const;		// Get delzz value
//void   delz(double dz);		// Set delzz value
//void   delzz(double dz);		// Set delzz value

MSVCDLL virtual double xi() const;		// Get xi value (Hz)
MSVCDLL virtual void   xi(double xval);	// Set xi value (Hz)
 
// ____________________________________________________________________________
// E                          OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//     Functions That Generate Strings to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

MSVCDLL std::string XiString() const;

//-----------------------------------------------------------------------------
//  Functions That Generate String Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

MSVCDLL virtual std::vector<std::string> IR2AStrings() const;

        // Input                R2I     : Rank 2 interaction (this)
        // Output               CSS     : Pointer to array of 5 strings
        // Note                         : The String array must be deleted
        //                                outside of this routine!

/*               Spin Quantum Number:       I
                 delzz:                 xxxxx.xx PPM
                 Asymmetry:                 x.xx
                 Down From PAS z-Axis:    xxx.xx Degrees
                 Over From PAS x-Axis:    xxx.xx Degrees                     */  


MSVCDLL virtual std::ostream& printAT(std::ostream& ostr, IST_type ISTT=UNK) const;

        // Input                R2I     : Rank 2 interaction (this)
        //                      ostr    : Output stream
        //                      ISTT    : Interaction type
        // Output               none    : Interaction spatial tensor
        //                                sent to the output stream
        // Note                         : This serves as a base print function
        //                                that the other IntName interaction 
        //                                derived classes use.  It must be                               
        //                                virutal. 

/*          Prints The Spherical Components, Spatial And Spin 
       These Will Be Printed Lines Which Will Look Like The Following 
                      (Repeated For All 5 m Values) 
 
                                                [ x.x, x.x, x.x] 
                A    = x.xxx            T     = [ x.x, x.x, x.x] 
                 2,m                     2,m    [ x.x, x.x, x.x] 
                                                [ x.x, x.x, x.x]             */   

 
//-----------------------------------------------------------------------------
//   Functions That Generate Ouput Of The Rank 2 Shift Anistropy Interaction
//-----------------------------------------------------------------------------
 
//  "print" is virtual as IntRank2 is a base class for other IntName classes

MSVCDLL virtual std::ostream& print(std::ostream& ostr, int fflag=0) const;
 
        // Input                R2I     : Rank 2 interaction (this)
        //                      ostr    : Output stream
        //                      fflag   : Format flag
        //                                  0 - Basic Parameters 
        //                                 !0 - Full output
        // Output               none    : CSA spatial tensor parameters
        //                                placed into the output stream
        // Note                         : This does NOT use the base class
        //                                virtual overload because we write
        //                                out two spin quantum values here?
 
    
MSVCDLL friend std::ostream& operator<< (std::ostream& out,const IntRank2& R2I);
 
        // Input                out     : Output stream;                           
        //                      R2I     : Rank 2 interaction 
        // Output                       : Modifies output stream


// ____________________________________________________________________________
// F                    SPIN TENSOR LINKED LIST ACCESS
// ____________________________________________________________________________

/* These functions allow users to examine the lists of interactions maintained
   by this class. Each new interaction will produce a spatial and spin tensor.
   Since many interactions share the same spin tensor, this class keeps them
   in a static list and just reuses ones already constructed.  Here one may
   view what the list(s) look like.                                        

           Input		ostr    : Output stream
                                R2I     : Rank 2 interaction (this)
           			ISL     : Rank 2 interaction list
           			X       : Rank 2 interaction type
				phdr    : Flag to print header info
                                fflag   : Format flag
                                            0 - Basic List Parameters
                                           !0 - Full List Output
	   Output		osr     : Output stream after having info
	   				  regarding ISL placed into it       */

MSVCDLL static std::ostream& printISLList(std::ostream& ostr,
                                              const std::list<IntRank2T>& ISL);
MSVCDLL static std::ostream& printList(std::ostream& ostr, bool fflag=false);
MSVCDLL static std::ostream& printList(std::ostream& ostr,IST_type X,int=0);

// ____________________________________________________________________________
// G                 RANK 2 INTERACTION COMPARISON FUNCTIONS
// ____________________________________________________________________________

/* Spin systems and other classes may wish to use linked lists of rank 2 
   interactions. Use of C++ STL list and vector classes can be used for such
   purposes if the ususal comparison functions are defined. Below are the 
   definitions.

           Input                R2I     : Rank 2 interaction (this)
                                IST1    : Another interaction spin tensor
           Output               T/F     : TRUE if R2I2 is 1st or == R2I
        ///F_list ==                    - Equality
        ///F_list !=                    - Inequality                         */
 
MSVCDLL int operator==(const IntRank2 &R2I2) const;
MSVCDLL int operator!=(const IntRank2 &R2I2) const;

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
   space of the interaction. These will returned in the product basis as
   simple matrices. Their units will be radians/sec.                         */

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

MSVCDLL virtual matrix H0() const;
MSVCDLL virtual matrix H0(double alpha, double beta, double gamma) const;
MSVCDLL virtual matrix H0(const EAngles& EA) const;

MSVCDLL virtual matrix H0(const std::vector<int>& HSs, int i)          const;
MSVCDLL virtual matrix H0(const std::vector<int>& HSs, int i,
                                           double A, double B, double G) const;
MSVCDLL virtual matrix H0(const std::vector<int>& HSs, int i, 
                                                      const EAngles& EA) const;

MSVCDLL virtual matrix H0(const std::vector<int>& HSs, int i, int j)   const;
MSVCDLL virtual matrix H0(const std::vector<int>& HSs, int i, int j,
                                           double A, double B, double G) const;
MSVCDLL virtual matrix H0(const std::vector<int>& HSs, int i, int j, 
                                                      const EAngles& EA) const;


// ----------------------------------------------------------------------------
//                          Full Interaction Hamiltonians
//                           These Use No Approximations!
//   They Are Correct ln the Lab. Frame, But Probably NOT in a Rotating Frame
// ----------------------------------------------------------------------------

/* These are virtual because they only return the irreducible rank 2 component
   of an interaction Hamiltonian. Dervived classes may be reducible rank 2
   interactions and will need to add the rank 0 and rank 1 components into 
   this in order to obtain the complete Hamiltonian. 

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

   
MSVCDLL virtual matrix H() const;
MSVCDLL virtual matrix H(double alpha, double beta, double gamma) const;
MSVCDLL virtual matrix H(const EAngles& EA) const;

MSVCDLL virtual matrix H(const std::vector<int>& HSs, int i, int j)   const;
MSVCDLL virtual matrix H(const std::vector<int>& HSs, int i, int j,
                                           double A, double B, double G) const;
MSVCDLL virtual matrix H(const std::vector<int>& HSs, int i, int j, 
                                                      const EAngles& EA) const;
};
 
#endif								// IntRank2.h
