/* IntQuad.h ****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Spatial Quadrupolar Tensor 		      Interface		**
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
** defined for a particular nuclear spin (or spin in a spin system).    **
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
** The following defintions are used herein (Auv normlized by delzz):	**
**                                                                      **
** 1.) PAS: |Azz| >= |Axx| >= |Ayy|					**
** 2.) PAS: Azz=1, eta=(Ayy-AXX)/2, Axx=(1+eta)/2, Ayy=(eta-1)/2        **
**                                                                      **
**                                   2					**
** 3.) Quadrupolar Coupling: NQCC = e qQ				**
**                                                 2			**
**                                   3*NQCC      3e qQ			**
** 4.) Quadrupolar Frequency: w   = --------  = --------		**
**                             Q    2I(2I-1)    2I(2I-1)		**
**                                                                      **
*************************************************************************/

#ifndef   IntQuad_h_			// Is file already included?
#  define IntQuad_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface
#  endif

extern const double xRT6PIO5; 		// Needed constant sqrt[6*PI/5]
extern const double xRT5O4PI;		// Needed constant sqrt[5/(4*PI)];
extern const double xRT5O24PI;		// Needed constant sqrt[5/(24*PI)];

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Basics/ParamSet.h>		// Include parameter sets
#include <IntRank2/IntRank2.h>		// Include parameter sets
#include <Matrix/complex.h>		// Include complex numbers
#include <Matrix/matrix.h>		// Include matrices
#include <Matrix/row_vector.h>		// Include row vectors
#include <iostream>
#include <string>

class IntQuad: public IntRank2
  {
  double _QCC;				// Quadruplolar coupling

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i             CLASS QUADRUPOLAR INTERACTION ERROR HANDLING
// ____________________________________________________________________________
 
/*      Input                QI      : Quadrupolar interaction (this)
                             eidx    : Flag for error type
                             noret   : Flag for return (0=return)
                             pname   : String included in message
        Output               none    : Error message
                                       Program execution stopped if fatal   */

void          Qerror(int eidx, int noret=0) const;
void          Qerror(int eidx, const std::string& pname, int noret=0) const;
volatile void Qfatal(int eidx) const;

// ____________________________________________________________________________
// ii                 QUADRUPOLAR INTERACTION SETUP FUNCTIONS
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

bool getQI(const ParameterSet& pset,
         double& Iqn, double& qcc, double& eta, EAngles& EA,
                                              int idx=-1, bool warn=true) const;

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
           Output               TF      : True if quadrupolar coupling found
                                          from parameters in pset
           Note                         : Interaction is NOT altered         */

bool getQCC(const ParameterSet& pset, double& qcc,
                                             int idx=-1, bool warn=true) const;

// ----------------------------------------------------------------------------
//                      Complete Quadrupolar Interaction
// ----------------------------------------------------------------------------

/* This function employs all of the "get*" functions in the above sections
   to parse a parameter set for all of the values needed to define a 
   quadrupolar interaction, namely { Iqn,QCC,eta,alpha,beta,gamma }. If the
   interaction definition is found, we set the interaction or return false.  */

bool setQI(const ParameterSet& pset, int idx=-1, bool warn=true);

// ____________________________________________________________________________
// iii            QUADRUPOLAR INTERACTION CHECKING FUNCTIONS
// ____________________________________________________________________________


bool checkIHS(int eidx=0, int warn=0);
bool checkI(int   eidx=0, int warn=0);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A          QUADRUPOLAR INTERACTION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

MSVCDLC IntQuad();
MSVCDLC IntQuad(const IntQuad &Q1);

// ----------------------------------------------------------------------------
//              Direct Constructors Using Spherical Components
// ----------------------------------------------------------------------------

/* Here we need to know the quantum number of the spins involved, the
   quadrupolar coupling, and the interaction asymmetry. the asymmetry
   is unitless and in the range [0,1]

           Input                Q       : Quadrupolar interaction (this)
                                II      : Spin I isotope type
                                Iqn     : Spin I quantum number (e.g. 1.5)
                                QCC     : Quadrupolar coupling value (Hz)
                                eta     : Tensor asymmetry value (default 0)
           Output               none    : Quadrupolar interaction constructed
           Note                         : A fatal error will result if
                                          the spin quantum number isn't > 1/2
           Note                         : Here QCC=delzz                     */

MSVCDLC IntQuad(const std::string& II, double qcc,
                                       double eta=0, const EAngles& EA=EAzero);
MSVCDLC IntQuad(const     Isotope& II, double qcc,
                                       double eta=0, const EAngles& EA=EAzero);
MSVCDLC IntQuad(          double   qn, double qcc,
                                       double eta=0, const EAngles& EA=EAzero);

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

           Input                Q       : Quadrupolar interaction (this)
                                II      : Spin I isotope type
                                Iqn     : Spin I quantum number (e.g. 1.5)
                                QCC     : Quadrupolar coupling value (Hz)
                                eta     : Tensor asymmetry value (default 0)
           Output               none    : Quadrupolar interaction constructed
           Note                         : A fatal error will result if
                                          the spin quantum number isn't > 1/2
           Note                         : Here QCC=delzz                     */

MSVCDLC IntQuad(const std::string& IsoI,const coord& Qxyz,const EAngles& EA=EAzero);
MSVCDLC IntQuad(const Isotope&  II,const coord& Qxyz,const EAngles& EA=EAzero);
MSVCDLC IntQuad(double          Iz,const coord& Qxyz,const EAngles& EA=EAzero);

// ----------------------------------------------------------------------------
//                     Construction Using Parameter Sets
// ----------------------------------------------------------------------------

/* These functions will construct the interaction from parameters in a
   specified GAMMA parameter set.  This is most useful when working with multi-
   spin systems that are setting up many interactions over the system using a
   single parameter file (or parameters read from an external ASCII file.

        Input                QI      : Quadrupolar interaction
                             pset    : Parameter set
                             idx     : Interaction index (default -1->none)
                             idxI    : Index for the spin
                             warn    : Flag to warn if no interaction found
					     0 = no warnings
					     1 = warnings
					     >1 = fatal warnings
        Output               none    : QI interaction constructed
                                       for spin with quantum number qn
                                       and parameters in pset
                                or   : QI interaction constructed
                                       for spin with index idxI
                                       and parameters in pset                */


MSVCDLC IntQuad(const ParameterSet& pset, int idx=-1, int warn=2);
//IntQuad(const std::string& II,  const ParameterSet& pset,
//                                                       int idx=-1, int warn=2);
//IntQuad(const Isotope& II, const ParameterSet& pset,
//                                                       int idx=-1, int warn=2);
//IntQuad(int idxI,         const ParameterSet& pset,             int warn=2);
//IntQuad(double qn,        const ParameterSet& pset, int idx=-1, int warn=2);
 
// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

MSVCDLL void operator= (const IntQuad &Q1);
MSVCDLC      ~IntQuad();

// ____________________________________________________________________________
// B                   SPATIAL TENSOR COMPONENT ACCESS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                      Spatial Tensor Anisotropy Value
// ----------------------------------------------------------------------------

/* Since quadrupolar interactions have no isotropic component, the interaction
   strength is entirely defined by the quadrupolar coupling (or delzz value).
   For this interaction the delzz value is the quadrupolar coupling, QCC, 
   and it relates to the quadruplar frequency w  and anisotropy as shown below.             
                                               Q
                                                                  2
                          2      2  ^                 3*NQCC     3e qQ
    DELZZ = NQCC = QCC = e qQ =  - / \ Q      w   = -------- = -------- 
                                 3 ---         Q    2I(2I-1)   2I(2I-1)
 
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

           Input                Q	: Quadrupolar interaction
                                delzz   : Nuclear quad. tensor delzz (Hz)
                                QCC     : Nuclear quad. coupling (Hz)
                                wQ      : Quadrupolar frequency (Hz)
           Output		delzz   : Nuclear quad. tensor delzz
	  				  value in (Hz)
           Note                         : This is equivalent to the
                                          quadrupolar coupling

    Again, for this interaction the delzz value is equivalent to the 
    Quadrupolar coupling constant. Also, note that since the quadrupolar
    frequency depends upon the spin quantum number it is not defined
    unless the spin is set. 

    Note that the default base class IntRank2A also provides the GAMMA
    normalized anisotropy and delzz values. These are constants given by

                       1/2
                  [ 5 ]                  ^      3
          del   = |---]  = 0.51503      /_\ A = - del   = 0.77255
             zz   [6PI]                         2    zz                      */



// static double IntRank2A::delzz();                            INHERITED
// static double IntRank2A::delA();                             INHERITED

MSVCDLL double QCC()   const;			// Set quadrupolar coupling
MSVCDLL double NQCC()  const;			// Set quadrupolar coupling
MSVCDLL double wQ()    const;			// Set quadrupolar frequency
MSVCDLL void   QCC(double  dz);		// Set quadrupolar coupling
MSVCDLL void   NQCC(double dz);		// Set quadrupolar coupling
MSVCDLL void   wQ(double    W);		// Set quadrupolar frequency

//-----------------------------------------------------------------------------
//                             Asymmetry Access
//-----------------------------------------------------------------------------

// Note that eta spans [0, 1] and is defined to be (Axx-Ayy)/Azz where
// |Azz| >= |Ayy| >= |Axx| in the treatment contained herein.  These functions
// are inherited from the base class IntRank2A

// double IntRank2A::eta( ) const;                              INHERITED
// void   IntRank2A::eta(double Eta);                           INHERITED

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
  row_vector IntRank2A::CartComps(double theta, double phi=0) const;        */

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

MSVCDLL double qxx() const;
MSVCDLL double qyy() const;
MSVCDLL double qzz() const;
MSVCDLL double qyx() const;
MSVCDLL double qxy() const;
MSVCDLL double qzx() const;
MSVCDLL double qxz() const;
MSVCDLL double qzy() const;
MSVCDLL double qyz() const;

MSVCDLL double qxx(double alpha, double beta, double gamma) const;
MSVCDLL double qyy(double alpha, double beta, double gamma) const;
MSVCDLL double qzz(double alpha, double beta, double gamma) const;
MSVCDLL double qyx(double alpha, double beta, double gamma) const;
MSVCDLL double qxy(double alpha, double beta, double gamma) const;
MSVCDLL double qzx(double alpha, double beta, double gamma) const;
MSVCDLL double qxz(double alpha, double beta, double gamma) const;
MSVCDLL double qzy(double alpha, double beta, double gamma) const;
MSVCDLL double qyz(double alpha, double beta, double gamma) const;

MSVCDLL double qxx(const EAngles& EA) const;
MSVCDLL double qyy(const EAngles& EA) const;
MSVCDLL double qzz(const EAngles& EA) const;
MSVCDLL double qyx(const EAngles& EA) const;
MSVCDLL double qxy(const EAngles& EA) const;
MSVCDLL double qzx(const EAngles& EA) const;
MSVCDLL double qxz(const EAngles& EA) const;
MSVCDLL double qzy(const EAngles& EA) const;
MSVCDLL double qyz(const EAngles& EA) const;

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

MSVCDLL matrix T21m1() const;
MSVCDLL matrix T22m2() const;

MSVCDLL matrix T21m1(const std::vector<int>& HSs, int i) const;
MSVCDLL matrix T22m2(const std::vector<int>& HSs, int i) const;

// ____________________________________________________________________________
// D           QUADRUPOLAR INTERACTION CONSTANT ACCESS FUNCTIONS
// ____________________________________________________________________________
 
/* The quadrupolar interaction constant is what sets the overall interaction
   strength and allows GAMMA to track different interaction types using 
   "normalized" spatian and spin tensor (classes IntRank2A and IntRank2T). For
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

MSVCDLL double xi( ) const;
 
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
      W (theta,phi)  = - W  [ 3*cos (theta) - eta*sin (theta)cos(2*phi) ]
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

MSVCDLL static double wQ2QCC(double QwQ, double I);
MSVCDLL static double QCC2wQ(double QCC, double I);

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
                                theta	: Orientation angle (degrees)
                                phi	: Orientation angle (degrees)
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
	   Note				: If no angles input the result is
	  				  for the current tensor orientation */
 
MSVCDLL double wQoriented()                           const;
MSVCDLL double wQ0()                                  const;
MSVCDLL double wQoriented(double theta, double phi=0) const;
MSVCDLL double wQ0(double        theta, double phi=0) const;
MSVCDLL matrix wQoriented(int   Ntheta, int   Nphi)   const;
MSVCDLL matrix wQ0(int          Ntheta, int   Nphi)   const;

// ----------------------------------------------------------------------------
//                       2nd Order Quadrupolar Frequency 
// ----------------------------------------------------------------------------

/* This is the splitting between transtions when the quadupolar interaction is
   weaker but not much much weaker than the Zeeman interaction.  These will
   only apply when that condition is applicable and should be added to the
   first order frequency (shifts).

   [2]              -2   [ Xi ]2 [
  w   (theta,phi) = -- * | -- |  | A A  (theta,phi)*[24m(m-1) - 4I(I+1) + 9]
   Q                Om   [ 2  ]  [  1 -1
    m-1,m
                                  1                                           ]
                                + -*A A  (theta,phi)*[12m(m-1) - 4I(I+1) + 6] |
                                  2  2 -2                                     ]

   Functions are provided which 1.) return a single shift given an input field
   strength (Omega) and 2.) return an array of components to facilitate powder 
   averaging.

           Input                Q       : Quadrupolar interaction
                                Om      : Field Strength (Larmor in Hz)
                                m       : Spin angular momentum z quantum #
                                theta   : Orientation angle (degrees)
                                phi     : Orientation angle (degrees)
                                Ntheta  : Number of orientation increments
                                Nphi    : Number of orientation increments
           Return               wQ      : 2nd order quadrupolar frequency
                                          shift (Hz) for the transition m-1,m
                                mx      : Array of components for building
                                          second order quad. frequency shifts
           Note                         : Theta spans [0, 180]
           Note                         : Phi spans [0, 360)                 */
 
 
MSVCDLL double wQ1(double Om, double m, double theta, double phi=0) const;
MSVCDLL double wQ1(double Om, double m) const;
MSVCDLL matrix wQ1(int Ntheta, int Nphi);
 
// ----------------------------------------------------------------------------
//                       2nd Order Central Transition Shifts
// ----------------------------------------------------------------------------

/* These functions apply only to spins have I = m * 1/2 where m is odd & > 1!
   The value(s) returned are from application of 2nd order perturbation theory
   when the quadrupolar interaciton is weaker, but not much much weaker, than
   the Zeeman interaction.  These will only apply when that condition is 
   applicable and should be added to the first order frequency (shifts).
 
                    2  
                 - w
  (2)               Q   [          3 ]   [               4
 w (theta,phi) = ---- * | I(I+1) - - | * | A(phi,eta)*cos (theta)]
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
                 8   3          4                  8

   Functions are provided which 1.) return a single shift given an input field
   strength (Omega) and 2.) return an array of components to facilitate powder 
   averaging.

           Input                Q       : Quadrupolar interaction
           Input                Om      : Field Strength (Larmor in Hz)
                                Ntheta  : Number of orientation increments
                                Nphi    : Number of orientation increments
           Return               wQ      : Shift in the central transition
                                          frequency due to 2nd order effects
                                          in an I=n*1/2 n=3,5,7,.... spin
           Return               mx      : Array of components for building
                                          second order quad. frequency shifts
                                          of the central transition over
                                          the sphere.
           Note                         : The return is zero if I not proper
           Note                         : Theta spans [0, 180]               */

MSVCDLL double wQcentral(double Om) const;
MSVCDLL matrix wQcentral(int Ntheta, int Nphi);




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

MSVCDLL             operator    ParameterSet( ) const;
MSVCDLL friend void operator+= (ParameterSet& pset, const IntQuad &Q);
MSVCDLL        void PSetAdd(ParameterSet& pset, int idx=-1, int pfx=-1) const;


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

MSVCDLL int write(const std::string &fn,int idx=-1,int pfx=-1,int wn=2) const;
MSVCDLL int write(std::ofstream& ofstr, int idx=-1,int pfx=-1,int wn=2) const;

// ____________________________________________________________________________
// G                 QUADRUPOLAR INTERACTION INPUT FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//      Direct Read of Q Intraction From An ASCII File Or A Parameter Set
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
	   Output		none    : Quadrupolar factor interaction is
					  read in from parameters in file
					  filename or those in parameter set */

MSVCDLL bool read(const std::string &filename, int idx=-1, int warn=2);
MSVCDLL bool read(const ParameterSet& pset,    int idx=-1, int warn=2);

// ____________________________________________________________________________
// H                    QUADRUPOLAR HAMILTONIAN FUNCTIONS
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

   The Hamiltonians will be returned in the single spin Hilbert space of the
   interaction unless a composite Hilbert space and some spin index is
   supplied. All Hamiltonians will returned in the product basis as simple
   matrices. Their units will be Hz.                                         */

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

MSVCDLL matrix H0( ) const;
MSVCDLL matrix H0(double alpha, double phi, double gamma) const;
MSVCDLL matrix H0(const EAngles& EA) const;

MSVCDLL matrix H0(const std::vector<int>& HSs, int i) const;
MSVCDLL matrix H0(const std::vector<int>& HSs, int i,
                               double alpha, double beta, double gamma) const;
MSVCDLL matrix H0(const std::vector<int>& HSs, int i, const EAngles& EA) const;
 

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

MSVCDLL matrix H1(double Om)                                          const;
MSVCDLL matrix H1(double Om, double alpha, double beta, double gamma) const;
MSVCDLL matrix H1(double Om, const EAngles& EA)                       const;

MSVCDLL matrix H1(std::vector<int>HSs, int i, double Om)              const;
MSVCDLL matrix H1(std::vector<int>HSs, int i, double Om,
                              double alpha, double beta, double gamma) const;
MSVCDLL matrix H1(std::vector<int>HSs, int i, double Om,
                              const EAngles& EA)                       const;

 
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

        // Input                Om      : Field Strength (Larmor in Hz)
        //                      theta   : Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HQ      : The 2nd order secular part of the
        //                                quadrupolar Hamiltonian (default
        //                                basis, Hz) for the interaction
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note                         : Also called the 1st order quadrupolar
        //                                interaction (perturbation theory)
        // Note                         : Rotationally invariant about z
        // Note                         : This will return in the spin Hilbert
        //                                space of dimension 2I+1 unless the
        //                                composite Hilbert space is given    */

MSVCDLL matrix Hw(double Om) const;
MSVCDLL matrix Hw(double Om, double theta, double phi=0) const;
MSVCDLL matrix Hw(std::vector<int>HSs, int i, double Om) const;
MSVCDLL matrix Hw(std::vector<int>HSs, int i, double Om, double theta, double phi=0) const;

// ----------------------------------------------------------------------------
//                 Full Quadrupolar Interaction Hamiltonians
//       These Have No Approximations & Are Independent of Field Strength
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

MSVCDLL matrix H( ) const;
MSVCDLL matrix H(double theta, double phi=0) const;
MSVCDLL matrix H(const std::vector<int>& HSs, int i) const;
MSVCDLL matrix H(const std::vector<int>& HSs, int i, double T, double P=0) const;
 
// ____________________________________________________________________________
// N                 QUADRUPOLAR HAMILTONIAN FRIEND FUNCTIONS
// ____________________________________________________________________________

// These Don't Depend on Class IntQuad Variables At All.  The Functions Just
// Mimic The Class' Hamiltonian Generation.


//friend matrix HQ0(double qn, double wQo, double eta=0, 
//                                                 double theta=0, double phi=0);
 
        // Input                qn      : Quantum number (1, 1.5, 2.5,...)
        //                      wQo     : PAS Quadrupolar frequency
        //                      eta     : Quadrupolar asymmetry [0,1]
        //                      theta   : Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               H0      : The secular part of the quadrupolar
        //                                Hamiltonian (default basis, Hz)
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note                         : Also called the 1st order quadrupolar
        //                                interaction (perturbation theory)
        // Note                         : Rotationally invariant about z
        // Note                         : This will return in the spin Hilbert
        //                                space of dimension 2I+1

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



//friend matrix HQ1(double Om, double qn, double wQo, double eta,
//                                              double theta=0.0, double phi=0.0);
                                                        
        // Input                Om      : Field Strength (Larmor in Hz)         
        //                      qn      : Quantum number (1, 1.5, 2.5,...)
        //                      wQo     : PAS Quadrupolar frequency
        //                      eta     : Quadrupolar asymmetry [0,1]
        //                      theta   : Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HQ1     : The 2nd order secular part of the
        //                                quadrupolar Hamiltonian (default
        //                                basis, Hz) for the interaction
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note                         : Also called the 1st order quadrupolar
        //                                interaction (perturbation theory)
        // Note                         : Rotationally invariant about z
        // Note                         : This will return in the spin Hilbert
        //                                space of dimension 2I+1


// ----------------------------------------------------------------------------
//      Functions To Make A Parameter Set From A Quadrupolar Interaction
// ----------------------------------------------------------------------------


MSVCDLL int set(const ParameterSet& pset, double I, int idx=-1);

        // Input                Q       : Quadrupolar interaction
	//  			pset	: Parameter set
        //                      I       : Quantum I value (0.5, 1, 3.5, ...)
	//			idx	: Tensor index (default -1)
	// Output		TF      : AQ is set by the parameters
	//			          in the parameter set pset
	//				  Return is FALSE if parameter
	//				  was not found
	// Note				: Normally, the += function would
	//				  suffice, but here we support the 
	//				  ability to specify an index in
	//				  the parameter name(s) used


 
// ----------------------------------------------------------------------------
//      Functions To Make A Quadrupolar Interaction From A Parameter Set
// ----------------------------------------------------------------------------

 
MSVCDLL int setAsph(const ParameterSet& pset, int idx=-1);

        // Input                Q       : Quadrupolar interaction (this)
        //                      pset    : A parameter set
        //                      idx     : Index value
        // Output               none    : Quadrupolar interaction
        //                                set from parameters in pset

/*----------------------------------------------------------------------------
**                Try & Read AQ Quadrupolar Coupling Spatial Tensor         **
**                         { delzz, eta, phi, theta }                       **
**                                                                          **
** AQ   (2) : ( DCC, eta, theta, phi )          - Quadrupolar Interaction   **
**                                                                          **
*****************************************************************************/
 
// Look for AQ parameter. Currently allowed parameters are the following:
//
//  AQ,    AQ(#)        - Quadrupolar Coupling Specified in kHz
//  AQkHz, AQkHz(#)     - Quadrupolar Coupling Specified in kHz
//  AQKHz, AQKHz(#)     - Quadrupolar Coupling Specified in kHz
//  AQHz,  AQHz(#)      - Quadrupolar Coupling Specified in Hz
//  AQMHz, AQMHz(#)     - Quadrupolar Coupling Specified in MHz
 
 
/*----------------------------------------------------------------------------
**                Try & Read T_Q GENERIC (QUAD) Spatial Tensor              **
**                        (This Should Be Avoided!!!!)                      **
**                                                                          **
** T_Q  (4) : 2                                 - Quadrupolar Tensor (kHz)  **
**            ( x.xx,  DCC(KHz),  eta)                                      **
**            ( theta, phi,       x.xx)                                     **
**                                                                          **
*****************************************************************************/
 
// The parameter T_Q is a quadrupolar spatial tensor set as a generic rank
// 2 (not irreducible) spatial tensor.  In order to read this easily I would
// have to include class space_T, the very class I am bypassing with the use
// of rank 2 interactions!  So, I will allow use of T_Q to set up IntQuad BUT
// I'll issue a warning AND (ugh!) I'll have redefine the parameter parse.

 
MSVCDLL int setAsphGen(const ParameterSet& pset, int idx=-1);
 
        // Input                Q       : Quadrupolar interaction (this)
        //                      pset    : A parameter set
        //                      idx     : Index value
        // Output               none    : Quadrupolar interaction
        //                                set from parameters in pset



        // Input                Q       : Quadrupolar interaction (this)
        //                      pset    : A parameter set
        //                      idx     : Index value
        // Output               none    : Quadrupolar interaction spin
        //                                quantum number is set from
        //                                parameters in pset
 
// ____________________________________________________________________________
// H                 QUADRUPOLAR INTERACTION INPUT FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//           Interactive Input From ASCII File Or From Console
// ----------------------------------------------------------------------------

/* These functions will query the user for information on building up a
   quadrupolar interaction.  The simplest form will just ask for a file that
   has the interaction parameters.  More complex forms will request specific
   interaction values.
  
           Input                Q       : Quadrupolar interaction (this) 
                                argc    : Number of arguments
                                argv    : Vecotr of argc arguments 
                                argn    : Argument index 
                                idx     : Index for spin 
           Output               string  : The parameter argn of array argc is 
                                          used to supply a filename from which 
                                          the quadrupolar interaction is read 
                                          If argument argn is not in argv, the 
                                          user is asked to supply a filename 
                                          The set filename is returned 
           Note                         : The file should be an ASCII file 
                                          containing known IntQuad parameters 
           Note                         : The interaction Q is modifed (filled) 
           Note                         : Since no spin type is given, this will
                                          also use the parameter Iso(#) to
                                          the value of the spin quantum number
*/

 
MSVCDLL std::string ask_read(int argc, char* argv[], int argn, int idx=-1);
MSVCDLL void ask(int argc, char* argv[], int& qn, double& QI,
        double& Qnqcc, double& Qeta, double& theta, double& Qphi, int Qflag=0);
 
        //                      QI	: Spin quantum number
        //                      Qnqcc   : Quad. coupling constant (Hz)
        //                      Qeta    : Quad. asymmetry 
        //                      theta  : Quad. orientation angle
        //                      Qphi    : Quad. orientation angle
	//			Qflag   : Flag is QCC or wQ requested
        // Output               none    : The values of qn, I, Qnqcc, Qeta,
        //                                theta, and Qphi are set herein 
        // Note                         : This is INTERACTIVE! 

 
MSVCDLL void askset(int argc, char* argv[], int& qn, int Qflag=0);
MSVCDLL void askset(int Qflag=0);
 
        // Input                Q       : Quadrupolar interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Array of arguments
        //                      qn      : Query value
	//			Qflag   : Flag is QCC or wQ requested
        // Output               none    : Q is set interactively
        // Note                         : This is INTERACTIVE!
 
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

   Frequency:             xxxxx.xx xHz       [A  , A  , A  ]
   Spin Quantum Number:       I              [ xx   xy   xz]   [ x.x, x.x, x.x]
   Coupling Constant:     xxxxx.xx xHz       [A  , A  , A  ] = [ x.x, x.x, x.x]
   Asymmetry:                 x.xx           [ yx   yy   yz]   [ x.x, x.x, x.x]
   Down From PAS z-Axis:    xxx.xx Degrees   [A  , A  , A  ]
   Over From PAS x-Axis:    xxx.xx Degrees   [ zx   zy   zz]
  
                                              m = [0,4] => {0,1,-1,2,-2}     */

// string* IntRank2T::TStrings(int M) const;                     INHERITED
MSVCDLL std::vector<std::string> CartAStrings(const std::string& CSForm) const;
MSVCDLL std::vector<std::string> SphAStrings()                           const;


//-----------------------------------------------------------------------------
//    Functions That Generate Ouput Of The Rank 2 Quadrupolar Interaction
//-----------------------------------------------------------------------------

/* These functions will output information concerning the Quadrupolar
   interaction to any output stream.

        // Input                Q       : Quadrupolar interaction (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Basic Parameters
                                           !0 - Full output
                                nrm     : Flag if GAMMA normalized output
           Output               none    : Quad interaction parameters
                                          placed into the output stream      */

MSVCDLL        std::ostream& print(std::ostream& out, int fflag=-1) const;
MSVCDLL friend std::ostream& operator<<    (std::ostream& out, const IntQuad& Q);

//-----------------------------------------------------------------------------
//  Functions That Generate Ouput Of Cartesian and Spherical & Cartesian A
//-----------------------------------------------------------------------------

//std::ostream& printSpherical(std::ostream& ostr);
//std::ostream& printCartesian(std::ostream& ostr);
//std::ostream& printCartesian(std::ostream& ostr, double theta, double phi=0);
  };

#endif								// IntQuad.h
