/* IntG.h *******************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Electron G Interaction			Interface		**
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
**  I       - Spin quantum number of electron involved			**
**  XI      - Interaction strength (radians/sec)			**
**                                                                      **
** Lastly, this class will itself maintain:				**
**                                                                      **
**  GISO    - Isotropic electron g-factor (unitless).			**
**  DELZZ   - Rank 2 spatial tensor delzz value (unitless).		**
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

#ifndef   IntG_h_			// Is file already included?
#  define IntG_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface
#  endif

class complex;				// Know about complex numbers
class coord;				// Know about coordinates
class matrix;				// Know about matrices
class row_vector;			// Know about row vectors
class spin_sys;				// Know about basic spin systems
class IntRank2A;			// Know about base class 

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Basics/ParamSet.h>		// Include parameter sets
#include <IntRank2/IntRank2A.h>		// Include base class headers
#include <IntRank2/IntRank2T.h>		// Include base class headers
#include <IntRank2/IntRank2.h>		// Include base class headers

class IntG: public IntRank2
  {
  double GISO;				// G isotropic value (gxx+gyy+gzz)/3
  double DELZZ;				// G anisotropic value (gzz - giso)
  double BoFIELD;			// Static magnetic field (Gauss)

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                  G FACTOR INTERACTION ERROR HANDLING
// ____________________________________________________________________________
 
/*       Input                GI      : Electron G interaction (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

         void IGerror(int eidx,                           int noret=0) const;
         void IGerror(int eidx, const std::string& pname, int noret=0) const;
volatile void IGfatal(int eidx)                                        const;
volatile void IGfatal(int eidx, const std::string& pname)              const;

// ____________________________________________________________________________
// ii                G FACTOR INTERACTION SETUP FUNCTIONS
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

bool getGI(const ParameterSet& pset,
         double& Iqn, double& g, double& gA, double& eta, EAngles& EA,
                                 double& Bo, int idx=-1, bool warn=true) const;

bool getGI(const ParameterSet& pset, const Isotope& ISI,
                    double& g, double& gA, double& eta, EAngles& EA,
                                 double& Bo, int idx=-1, bool warn=true) const;

// ----------------------------------------------------------------------------
//            Get Isotropic G Value (G-Factor) From A Parameter Set
// ----------------------------------------------------------------------------

/*         Input                SA      : CSA interaction (this)
                                pset    : A parameter set
                                idx     : Index value
                                warn    : Warning output flag
           Output               g       : The g factor (unitless)
                                          from parameters in pset
           Note                         : This WILL NOT alter the interaction
           Note                         : Parameters are G, G(#), g, g(#)
           Note                         : The value is unitless              */

bool getGIso(const ParameterSet& pset, double& g,
                                             int idx=-1, bool warn=true) const;

// ----------------------------------------------------------------------------
//                Get Anisotropic G Value From A Parameter Set
// ----------------------------------------------------------------------------

/*         Input                SA      : CSA interaction (this)
                                pset    : A parameter set
                                gA      : The g anisotropy (unitless)
                                idx     : Index value
                                warn    : Warning output flag
           Output               TF      : True if g anisotropy (unitless)
                                          set from parameters in pset
           Note                         : This WILL NOT alter the interaction
           Note                         : Parameters are GA,GA(#),gA,gA(#)   */

bool getGA(const ParameterSet& pset, double& gA,
                                             int idx=-1, bool warn=true) const;

// ----------------------------------------------------------------------------
//          Functions To Read The Applied Field Or Base Frequency
// ----------------------------------------------------------------------------

/*         Input                SA      : CSA interaction (this)
                                pset    : A parameter set
                                Bo      : The field strength (Gauss)
                                idx     : Index value
                                warn    : Warning output flag
           Output               TF      : True if Bo (Gauss)
                                          set from parameters in pset
           Note                         : This WILL NOT alter the interaction
           Note                         : Allowd parameters are listed below

                        Field,  FieldT, gField(#), gFieldT(#)
                        GOmega, Omega,  GField(#), GFieldT(#)
                                        gOmega(#), GOmega(#)                 */

bool getField(const ParameterSet& pset, double& Bo,
                                             int idx=-1, bool warn=true) const;


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

bool setGI(const ParameterSet& pset,         int idx=-1, bool warn=true);
bool setGI(const Isotope& II, const ParameterSet& pset,
                                                   int idx=-1, bool warn=true);

// ____________________________________________________________________________
// iii           ELECTRON G INTERACTION SPECIAL SPIN TENSORS
// ____________________________________________________________________________

/* For this interaction the values { BoFIELD, DELZZ, _XI } are inter-related.
   If any one of them is set then the others may need to be updated. These
   functions will do this updating.

    Function                             Purpose
   -----------   --------------------------------------------------------------
      setXi      Sets the value of _XI (in IntRank2) for the current values of
                 DELZZ and BoFIELD
      setBo      Sets the value of BoFIELD then adjusts the value of _XI
                 based on this value & DELZZ. Bo assumed in Gauss           */

void setXi();
void setBo(double Bo);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A               G FACTOR INTERACTION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/* Electron g interaction construction. For full construction users need to
   define the 1.) irreducible rank 2 spatial tensor, 2.) irreducible rank 2
   spin tensors, and 3.) interaction strength.  Since GAMMA uses normalized
   spatial and spin tensors, the first two are relatively simple. For the 
   spatial tensor we need { eta, theta, phi } or { gxx, gyy, gzz, theta, phi }
   and for the spin tensor we need only the spin quantum value involved,
   normally 0.5.  However, the scaling is more variable since it is tied up in
   the isotropic value of the reducible rank 2 spatial G tensor (i.e. the
   g-factor) as well as in the external field strength that the electron is in.
   Since g = (1/3)[gxx+gyy+gzz], if the interaction is specified using these
   Cartesian components the isotropic values is automatically set.  If not then
   the isotropic value should be set independently or it will be taken to be
   that of the free electron (gfree = 2.00231928). The anisotropic part of the
   spatial tensor needs to be set using either a field strength specification
   or a resonance frequency specification. These may also be done after the
   construction has taken place.                                             */

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

MSVCDLC IntG();					// Null g interaction
MSVCDLC IntG(const IntG &GI1);			// Duplicate g interaction

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
  
        Input      GI	   : G factor interaction (this)
		   g 	   : Isotropic G value (unitless)
	           gA      : G anisotropy      (unitless) 
	           geta    : G asymmetry value [0,1]        
		   EA      : G orientation Euler angles (radians)
		   Bo      : External field strength    (Gauss)
	Output	   none    : GI interaction constructed                      */

MSVCDLC IntG(const  std::string& Iso, double g, double gA,
                        double geta=0,  const EAngles& EA=EAzero, double Bo=0);
MSVCDLC IntG(const  Isotope&     Iso, double g, double gA,
                        double geta=0,  const EAngles& EA=EAzero, double Bo=0);
MSVCDLC IntG(double Iqn,              double g, double gA,
                        double geta=0,  const EAngles& EA=EAzero, double Bo=0);

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
   default values are 0. Note that construction using individual {Gxx,Gyy,Gzz}
   values is not supported because of conflicts with spherical constructors.

        Input      GI	   : G factor interaction (this)
                   Gcart   : Cartesian tensor values (unitless)
		   EA      : G orientation Euler angles (radians)
		   Bo      : External field strength    (Gauss)
	Output	   none    : GI interaction constructed                      */

MSVCDLC IntG(const std::string& I, const coord& Gcart, 
                                         const EAngles& EA=EAzero, double B=0);
MSVCDLC IntG(const Isotope&     I, const coord& Gcart,
                                         const EAngles& EA=EAzero, double B=0);
MSVCDLC IntG(double             I, const coord& Gcart,
                                         const EAngles& EA=EAzero, double B=0);
MSVCDLC IntG(                      const coord& Gcart,
                                         const EAngles& EA=EAzero, double B=0);

// ----------------------------------------------------------------------------
//                    Construction Using Parameter Sets
// ----------------------------------------------------------------------------
 
/* These functions will construct the interaction from parameters in a
   specified GAMMA parameter set.  This is most useful when working with multi-
   spin systems that are setting up many interactions over the system using a
   single parameter file (or parameters read from an external ASCII file.  

        Input                GI      : G factor interaction
                             pset    : Parameter set
                             idx     : Interaction index (default -1->none)
                             idxI    : Index for the spin
                             warn    : Flag to warn if no interaction found
        Output               none    : GI interaction constructed
                                       for spin with quantum number qn
                                       and parameters in pset
                                or   : GI interaction constructed
                                       for spin with index idxI
                                       and parameters in pset                */


MSVCDLC IntG(                  const ParameterSet& pset, int idx=-1, int warn=2);
MSVCDLC IntG(const Isotope& II,const ParameterSet& pset, int idx=-1, int warn=2);
 
// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

MSVCDLL void operator= (const IntG &GI1);
MSVCDLC      ~IntG();

// ____________________________________________________________________________
// B                   SPATIAL TENSOR COMPONENT ACCESS
//
//                     { Giso, del  , eta, theta, phi }
//                                zz
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
// 			   Isotropic Value Access
// ----------------------------------------------------------------------------
 
/* The isotropic value is 
                                    1 [               ]  
                             g    = - | g  + g  + g   |
                              iso   3 [  xx   yy   zz ]

   &* should be a value roughly near that of the free electron (2.00231928). */
 
MSVCDLL double iso() const;
MSVCDLL double g()   const;
MSVCDLL void   iso(double giso);
MSVCDLL void   g(double   giso);

// ----------------------------------------------------------------------------
// 			Spatial Tensor Anisotropy Value
// ----------------------------------------------------------------------------

/* The g anisotropy is  

               ^      3                        1
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

MSVCDLL double aniso() const;
MSVCDLL double gA()    const;
MSVCDLL void   aniso(double ga);
MSVCDLL void   gA(double    ga);

MSVCDLL double gdelz() const;
MSVCDLL void   gdelz(double gdz);


// ----------------------------------------------------------------------------
// 			Spatial Tensor Asymmetry Value
// ----------------------------------------------------------------------------

// Note that eta spans [0, 1] and is defined to be (Axx-Ayy)/Azz where
// |Azz| >= |Ayy| >= |Axx| in the treatment contained herein.  These functions
// are inherited from the base class IntRank2A

//double eta( ) const; 		INHERITED	Get the asymmetry
//void   eta(double Geta);	INHERITED	Set the asymmetry

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
 
/* These allow one to access the Cartesian elements of the full (unscaled) G
   tensor at a specific orientation without rotating the entire tensor.  In
   these functions theta is the angle down from the lab. frame z-axis and phi
   the angle over from the lab. frame x-axis.  If no angles are specified, then
   the current orientation is used.  Typically, these elements will be on the
   order of 2 (since gfree is 2.00231928). The relate to GAMMA's normlized
   Cartesian values according to (where Kdel is a Kronecker delta function)
 
                            1/2
                    [ 6*PI ]
              g   = | ---- |    * del   * A   + Kdel    * g
               uv   [  5   ]         zz    uv       u,v    iso               */
 
MSVCDLL double gxx() const;
MSVCDLL double gyy() const;
MSVCDLL double gzz() const;
MSVCDLL double gyx() const;
MSVCDLL double gxy() const;
MSVCDLL double gzx() const;
MSVCDLL double gxz() const;
MSVCDLL double gzy() const;
MSVCDLL double gyz() const;

MSVCDLL double gxx(double alpha, double beta, double gamma) const;
MSVCDLL double gyy(double alpha, double beta, double gamma) const;
MSVCDLL double gzz(double alpha, double beta, double gamma) const;
MSVCDLL double gyx(double alpha, double beta, double gamma) const;
MSVCDLL double gxy(double alpha, double beta, double gamma) const;
MSVCDLL double gzx(double alpha, double beta, double gamma) const;
MSVCDLL double gxz(double alpha, double beta, double gamma) const;
MSVCDLL double gzy(double alpha, double beta, double gamma) const;
MSVCDLL double gyz(double alpha, double beta, double gamma) const;

MSVCDLL double gxx(const EAngles& EA) const;
MSVCDLL double gyy(const EAngles& EA) const;
MSVCDLL double gzz(const EAngles& EA) const;
MSVCDLL double gyx(const EAngles& EA) const;
MSVCDLL double gxy(const EAngles& EA) const;
MSVCDLL double gzx(const EAngles& EA) const;
MSVCDLL double gxz(const EAngles& EA) const;
MSVCDLL double gzy(const EAngles& EA) const;
MSVCDLL double gyz(const EAngles& EA) const;

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

        1.) GAMMA normalized Auv - Done With mx()
        2.) Typical guv values   - Done With mx(true);
        3.) Shown in lab frame   - Done With This Function

   For case 3.) the values are related to the GAMMA normalized (Auv) and
   typically presented values (guv) according to

                            1/2
                    [ 6*PI ]
              G   = | ---- |    * del   * A   + Kdel    * g
               uv   [  5   ]         zz    uv       u,v    iso

   where Kdel is a Kronecker delta function.                                 */

MSVCDLL matrix Gmx()           const;		// Cartesian G spatial tensor
MSVCDLL double Field()         const;		// Get Base Field Bo  (Gauss)
MSVCDLL double GOmega()        const;		// Get Base Frequency (GHz)
MSVCDLL double Frequency()     const;		// Get Om Frequency   (GHz)
MSVCDLL void   Field(double      Bo);		// Set Bo Field       (Gauss)
MSVCDLL void   Frequency(double GOm);		// Set Om Frequency   (GHz)

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

   Electron G interaciton constants are field dependent and defined as

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

MSVCDLL double xiOm(double Om) const;
MSVCDLL double xiBo(double Bo) const;

MSVCDLL double xi()         const;		// Overwrites inherited function	
MSVCDLL void   xi(double X) const;		// Overwrites inherited function

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
           			Om      : The spectrometer frequency (GHz)   */

MSVCDLL static double gvalue(double Om, double H);
MSVCDLL static double Hvalue(double Om, double g);
MSVCDLL static double Omvalue(double H, double g);

// ----------------------------------------------------------------------------
//                    Functions To Get Effective G Factors
// ----------------------------------------------------------------------------
 
//     hfl           1       [     2                   2                   ]
//    g    =  g    + - del   | 3cos (theta)-1 + eta*sin (theta)*cos(2*phi) |
//     eff     iso   2    zz [                                             ]

MSVCDLL double geff_hfl() const;
MSVCDLL double geff_hfl(double theta, double phi) const;


// ----------------------------------------------------------------------------
//                      Functions To Get Resonance Frequency
//          As a Function of g(theta,phi) and Applied H Field Strength
// ----------------------------------------------------------------------------
 
//                                        -21
//              beta*g   *H(G)   9.2741x10   erg/G * g   * H(G)   g    * H(G)
//                    eff                             eff          eff
//   GOm(GHz) = -------------- = ------------------------------ = -----------
//                    h                        -18                 GHZ2GAUSS
//                                    6.6262x10   erg/GHz
 
        // Input                GI      : G factor interaction
        //                      H       : Static field (Gauss)
        // Return               void    : The field at resonance (Gauss)
        // Note                         : A symmetric tensor is assumed
        //                                in the para & perp functions
        // Note                         : When the field Ho ~ 3 kG and
        //                                g ~ 2.003 then GOm ~ 10 GHz
 
MSVCDLL double GOm_iso()  const;
MSVCDLL double GOm_para() const;
MSVCDLL double GOm_perp() const;
MSVCDLL double GOm_xx()   const;
MSVCDLL double GOm_yy()   const;
MSVCDLL double GOm_zz()   const;

MSVCDLL double GOm_iso(double  Ho) const;
MSVCDLL double GOm_para(double Ho) const;
MSVCDLL double GOm_perp(double Ho) const;
MSVCDLL double GOm_xx(double   Ho) const;
MSVCDLL double GOm_yy(double   Ho) const;
MSVCDLL double GOm_zz(double   Ho) const;


// ----------------------------------------------------------------------------
//                    Functions To Get Resonance Field Strength
//              As a Function of g(theta,phi) and RF-Field Frequency
// ----------------------------------------------------------------------------

//                                  -27
//            h * GOm      6.6262x10   erg/Hz * GOm(Hz)   GHZ2GAUSS*GOm(GHz)
//      H = -----------  = ---------------------------- = -----------------
//          beta * g                -21                          g
//                  eff    9.2741x10   erg/G * geff               eff

        // Input                GI      : G factor interaction
        //                      GOm     : Excitation frequency (GHz)
        // Return               void    : The field at resonance (Gauss)
        // Note                         : A symmetric tensor is assumed
        //                                in the para & perp functions
        // Note                         : When the field GOm ~ 10 GHz and
        //                                g ~ 2.003 then H ~ 3 kGauss

MSVCDLL double field_iso(double GOm)  const;
MSVCDLL double field_para(double GOm) const;
MSVCDLL double field_perp(double GOm) const;
MSVCDLL double field_xx(double GOm)   const;
MSVCDLL double field_yy(double GOm)   const;
MSVCDLL double field_zz(double GOm)   const;

// ____________________________________________________________________________
// F                G INTERACTION PARAMETER SET FUNCTIONS
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
   are { GI(#), gxx(#), gyy(#), gcc(#), GTheta(#), GPhi(#), GEta(#) }.

	   Input		GI	: G factor interaction
                                idx     : Interaction index (default -1)
                                pfx     : Interaction 2nd indx (def -1)
	   Output		pset	: Parameter set with only
	  			          electron G interaction parameters 

     Function                                 Purpose
   ------------         -------------------------------------------------------
   ParameterSet         Convert interaction into a parameter set
   operator +=          Adds interaction to existing parameter set
   PSetAdd              Adds interaction to existing parameter set, this
                        allows for an index and a prefix in the parameters   */

MSVCDLL operator ParameterSet( ) const;
MSVCDLL friend void operator+= (ParameterSet& pset, const IntG &GI);
MSVCDLL void PSetAdd(ParameterSet& pset, int idx=-1, int pfx=-1) const;

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

MSVCDLL bool write(const std::string &filename, int idx=-1, int warn=2) const;
MSVCDLL bool write(      std::ofstream& ofstr,  int idx=-1, int warn=2) const;

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

MSVCDLL bool read(const std::string&  fname, int idx=-1, int warn=2);
MSVCDLL bool read(const ParameterSet& pset,  int idx=-1, int warn=2);

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

MSVCDLL std::string ask_read(int argc, char* argv[], int argn,                    int idx=-1);
MSVCDLL std::string ask_read(int argc, char* argv[], int argn, const std::string& def, int idx=-1);

// ____________________________________________________________________________
// I                 ELECTRON G INTERACTION OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//   Functions That Generate Simple Strings To Simplify & Modularize Printing
//-----------------------------------------------------------------------------

MSVCDLL std::string GfactorString()    const;
MSVCDLL std::string GAString()         const;
MSVCDLL std::string GFieldString()     const;
MSVCDLL std::string GFrequencyString() const;

//-----------------------------------------------------------------------------
//  Functions That Generate String Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

/* These functions return spatial tensor information in string format. This is
   done to facilitate printing, in particular printing of G spatial tensors
   from within rank 2 interactions.  In this case we make a list of information
   thats usually displayed to the left of a 3x3 Cartesian matrix rep. of G.

                                   InfoStrings
                                   ===========

                       Spin Quantum Number:       I
                       Isotropic G Factor     xxxxx.xx
                       G Anisotrpy    :       xxxxx.xx
                       Asymmetry:                 x.xx
                       Euler Angle Alpha:       xxx.xx Degrees
                       Euler Angle Beta:        xxx.xx Degrees
                       Euler Angle Gamma:       xxx.xx Degrees


                                       TStrings
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


MSVCDLL std::vector<std::string> InfoStrings() const;
MSVCDLL std::vector<std::string> CartAStrings(const std::string& CSF) const;
// string* IntRank2T::TStrings(int M) const;                     INHERITED

//-----------------------------------------------------------------------------
//       Functions That Generate Ouput Of The Rank 2 G Interaction
//-----------------------------------------------------------------------------

/* These functions will output information concerning the G interaction to any
   output stream.

           Input                GI      : Electron G interaction (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            false - Basic Parameters
                                            true  - Full output
           Output               none    : G spatial tensor parameters
                                          placed into the output stream      */

MSVCDLL        std::ostream& print(std::ostream& out, bool fflag=false) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& out, const IntG& G);
MSVCDLL        std::ostream& printSpherical(std::ostream& ostr);
MSVCDLL static std::ostream& STList(std::ostream& ostr, int fflag=0);

// ____________________________________________________________________________
// J                    G FACTOR HAMILTONIAN FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                   The Isotropic G Interaction Hamiltonian
//                   (Rotationally Invariant About Any Axis)
// ----------------------------------------------------------------------------
 
/*   We don't use spherical tensors for this part, it doesn't rotate anyway!
     I've built a "rotating frame" allowance into Hiso!  This is very NMR-like
     but may never be used by ESR puritans.
 
            iso   beta * H              beta * H          G
           H    = -------- * g   * S  = -------- * A   * T    = GOm * S
            G         h       iso   z       h       0,0   0,0          z 
 
   Or in the free e- rotating frame:
 
         iso   beta * H              beta * H          
        H    = -------- * g   * S  = -------- * (g   - g ) * S  = GOm * S
         G         h       eff   z       h        iso   e     z          z   */
 
MSVCDLL matrix Hiso(bool rotfrm=false) const;
MSVCDLL matrix Hiso(std::vector<int> HSs, int i, bool rotfrm=false) const;

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
	   Note				: Does NOT include any isotropic terms
           Note                         : This will be referenced to the lab
                                          frame in accordance with standard
                                          ESR definitions of the G tensor.   */


//matrix H0()                             const;                INHERITED
//matrix H0(double A, double B, double G) const;                INHERITED
//matrix H0(const EAngles& EA)            const;                INHERITED

//matrix H0(const vector<int>& HSs, int i)                    const;  IHT
//matrix H0(const vector<int>& HSs, int i,
//                                    double A, double B, double G) const;  IHT
//matrix H0(const vector<int>& HSs, int i, const EAngles& EA) const;  IHT


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
	   Note				: Does NOT include any isotropic terms
           Note                         : No HSs: return in the spin Hilbert
                                          space of dimension (2I+1)
                                          With HSs: return in composite spin
                                          Hilbert space                      */

// matrix H( ) const                                            INHERITED
// matrix H(double alpha, double beta, double gamma) const      INHERITED
// matrix H(const EAngles& EA) const                            INHERITED













MSVCDLL matrix Hiso(double GOm, int rotflg=0) const;
MSVCDLL matrix Hiso(std::vector<int> HSs, int i, double GOm, int rotflg=0) const;
MSVCDLL matrix Hiso_lab(double GOm) const;
MSVCDLL matrix Hiso_lab(std::vector<int> HSs, int i, double GOm) const;
 

// ----------------------------------------------------------------------------
//               First Order G factor Interaction Hamiltonian
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
//  These Are For When The G factor Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------

//  The secular part of the G Hamiltonian contains only those Hamiltonian
//  components which explicitly commute with z axis rotations.  Here we have
 
  
//  [1]              G                 G      (0) 
// H  (the,phi) = Xi * A (the,phi) * T    = H  
//  G                   2,0           2,0    G 
 
//                1              [     2                   2               ] 
//              = - beta*H*del   | 3cos (the) - 1 + eta*sin (the)cos(2*phi)| Iz
//                2           zz [                                         ]  
  
//                      del 
//                         zz    [     2                   2               ] 
//              = GOm * ------ * | 3cos (the) - 1 + eta*sin (the)cos(2*phi)| Iz
//                      2*g      [                                         ]  
//                         iso 
  
// Identical to the Anisotropic G Hamiltonian in a high field limit!
// Note that this Hamiltonian must be ADDED to the isotropic G Hamiltonian.

MSVCDLL matrix H0X(double Om) const;
MSVCDLL matrix H0X(double Om, double theta, double phi=0) const;
MSVCDLL matrix H0X(std::vector<int> HSs, int i, double Om) const;
MSVCDLL matrix H0X(std::vector<int> HSs, int i, double Om, double theta, double phi=0) const;
 
        // Input                GI      : G factor interaction
        //                      Om      : Spectrometer field (GHz)
        // Output               H0      : Secular part of the g anisotropy
        //                                Hamiltonian (default basis, Hz)
        // Note                         : Also called the 1st order GI
        //                                interaction (perturbation theory)
        // Note                         : Rotationally invariant about z
        // Note                         : H0 spin Hilbert space dim. is 2
        // Note                         : This will be referenced to the lab
        //                                frame in accordance with standard
        //                                ESR definitions of the G tensor.
 
// These next two use the spherical tensor components explicitly.  The results 
// should match the those of the above two functions. 

MSVCDLL matrix H0Direct(double Om) const;
MSVCDLL matrix H0Direct(double Om, double theta, double phi=0) const;

// ----------------------------------------------------------------------------
//                Full G Interaction Anisotropy Hamiltonians
//       These Have No Approximations & Are Independent of Field Strength
//     (They Are Correct in The Laboratory Frame, NOT in a Rotating Frame)
// ----------------------------------------------------------------------------
 
//                      ---
//  GA              G   \       m    G             G
// H  (the,phi) = Xi  * /   (-1)  * A (the,phi) * T
//                      ---          2,m           2,-m
//                       m
 
//                    [       2                   2
//              = K * | [ 3cos (the) - 1 + eta*sin (the)cos(2*phi) ] S
//                    [                                               z
//                                                                           ]
//              + sin(the)*[ cos(the)(3-eta*cos(2*phi)S + eta*sin(2*phi)*S ] |
//                                                     x                  y  ]
// where
//                                                     del
//                          1                             zz
//                      K = - beta * H * del   = GOm * ------
//                          2               zz          2*g
//                                                         iso

// Note: There are no m = +/-2 spin components in rank 2 G factor interactions
//       when the applied field is along the z-axis in the laboratory frame
//   (the lab frame is the frame to which angles theta and phi are referenced)

        // Input                GI      : G factor interaction(this)
        //                      GOm     : Spectrometer field (GHz)
	//			Theta   : Orientation down from +z (degrees)
	//			Phi	: Orientation over from +x (degrees)
        // Output               H       : The Anisotropic G Hamiltonian
        // Note                         : These will return in the spin Hilbert
        //                                space of dimension 2
        // Note                         : These does NOT sum over the m=+/-2
        //                                components (i.e. field along +z)
        // Note                         : Xi depends on the frequency herein

MSVCDLL matrix HA(double GOm) const;
MSVCDLL matrix HA(double GOm, double Theta, double Phi=0) const;
MSVCDLL matrix HA(std::vector<int> HSs, int i, double GOm) const;
MSVCDLL matrix HA(std::vector<int> HSs, int i, double GOm, double T, double P=0) const;

// These next two use the spherical tensor components explicitly.  The results
// should match the those of the above two functions.

MSVCDLL matrix HADirect(double GOm) const;
MSVCDLL matrix HADirect(double GOm, double Theta, double Phi) const;

// ----------------------------------------------------------------------------
//                      Full G Interaction Hamiltonians
//       These Have No Approximations & Are Independent of Field Strength
//     (They Are Correct in The Laboratory Frame, NOT in a Rotating Frame)
// ----------------------------------------------------------------------------
 
/*                      iso   ---     m     G                      G
       H (theta,phi) = H    + \   (-1)  * Xi  * A   (theta,phi) * T
        G               G     /                  2,m               2,-m
                              ---
  
                        iso    A
                     = H    + H  (theta, phi)
                        G      G
 
        // Input                GI      : G factor interaction
        //                      Om      : Spectrometer frequency (GHz)
	//			Theta   : Orientation down from +z (degrees)
	//			Phi	: Orientation over from +x (degrees)
        // Output               H       : The G Hamiltonian (in Hz)
        // Note                         : This will return in the spin Hilbert
        //                                space of dimension 2
        // Note                         : This does NOT sum over the m=+/-2
        //                                components (i.e. field along +z)
        // Note                         : Xi depends on the frequency herein */

MSVCDLL matrix H(double Om) const;
MSVCDLL matrix H(double Om, double Theta, double Phi=0) const;
MSVCDLL matrix H(std::vector<int> HSs, int i, double Om) const;
MSVCDLL matrix H(std::vector<int> HSs, int i, double Om, double Theta, double Phi=0) const;

// These next two generate the Hamiltonians explicitly.  The results
// should match the those of the above two functions.

MSVCDLL matrix HDirect(double Om) const;
MSVCDLL matrix HDirect(double Om, double Theta, double Phi=0) const;



// ____________________________________________________________________________
// G              G FACTOR INTERACTION AUXILIARY FRIEND FUNCTIONS
// ____________________________________________________________________________


// ____________________________________________________________________________
// I              G FACTOR INTERACTION STANDARD INPUT FUNCTIONS
// ____________________________________________________________________________

 
//void ask(int argc, char* argv[], int& qn, double& CI,
//       double& Cnqcc, double& Ceta, double& theta, double& Cphi, int Cflag=0);
 
        // Input                GI      : G factor interaction (this)
        //                      argc    : Number of arguments 
        //                      argv    : Array of arguments 
        //                      qn      : Cuery value 
        //                      CI	: Spin quantum number
        //                      Cnqcc   : GI. coupling constant (Hz)
        //                      Ceta    : GI. asymmetry 
        //                      theta  : GI. orientation angle
        //                      Cphi    : GI. orientation angle
	//			Cflag   : Flag is CCC or wGIrequested
        // Output               none    : The values of qn, I, Cnqcc, Ceta,
        //                                theta, and Cphi are set herein 
        // Note                         : This is INTERACTIVE! 

 
//void askset(int argc, char* argv[], int& qn, int Cflag=0);
 
        // Input                GI      : G factor interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Array of arguments
        //                      qn      : Cuery value
	//			Cflag   : Flag is CCC or wGIrequested
        // Output               none    : GIis set interactively
        // Note                         : This is INTERACTIVE!
 

//void askset(int Cflag=0);

        // Input                GI      : G factor interaction (this)
        // Output               none    : GIis set interactively
	//			Cflag   : Flag is CCC or wGIrequested
        // Note                         : This is INTERACTIVE!
 
// sosik

//-----------------------------------------------------------------------------
//          Functions To Print The Tensor In Cartesian Format
//-----------------------------------------------------------------------------

        // Input                GI      : Electron G interaction (this)
        //                      ostr    : Output stream
        //                      T       : Orientation angle theta (degrees)
        //                      P       : Orientation angle phi   (degrees)
        //                      F       : Title print flag (default 1 = print)
        //                      CSF     : Element output format
        // Output               none    : Rank 2 spatial tensor parameters
        //                                sent to the output stream

/*           Output Some Printed Lines Which Will Look Like The Following


                       [A  , A  , A  ]
      Output From      [ xx   xy   xz]   [ x.x, x.x, x.x]
     Third Overload    [A  , A  , A  ] = [ x.x, x.x, x.x]
                       [ yx   yy   yz]   [ x.x, x.x, x.x]
                       [A  , A  , A  ]
                       [ zx   zy   zz]                                       */
 

MSVCDLL virtual std::ostream& printCartesian(std::ostream& ostr,
                                    const std::string& CSF="%6.3f", int tpf=2);
         
//std::ostream& printCartesian(std::ostream& ostr, double theta, double phi=0);
         
        // Input                G	: Electron G interaction (this)
        //                      ostr	: Output stream
        //                      phi     : Orientation angle (degrees)
        //                      theta   : Orientation angle (degrees)
        // Output               none	: Electron G interaction parameters set to
        //                                output stream
 

MSVCDLL std::ostream& printCartG(std::ostream& ostr, int tflag=1) const;
         
        // Input                GI      : G spatial tensor (this)
        //                      ostr    : Output stream
	//			tflag   : Title print flag
        // Output               none    : Full G Cartiesian tensor
        //                                sent to output stream
 
//                     [g  , g  , g  ] 
//                     [ xx   xy   xz]   [ x.x, x.x, x.x]
//                     [g  , g  , g  ] = [ x.x, x.x, x.x]
//                     [ yx   yy   yz]   [ x.x, x.x, x.x]
//                     [g  , g  , g  ]
//                     [ zx   zy   zz]
//
  };

#endif							// IntG.h
