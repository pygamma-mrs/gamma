/* IntCSA.h *****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Chemical Shift Anisotropy Interaction 	     Interface		**
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
** A variable of type IntCSA represents a chemical shift anisotropy 	**
** interaction defined for a particular nuclear spin (or spin in a spin	**
** system).  The interaction is maintained in the form of spatial and	**
** spin tensors which are stored in irreducible	spherical format.	**
**                                                                      **
** 1.) rank = 2         2.) symmetric           3.) eta = [0, 1]        **
**                                                                      **
** The tensor will internally be stored in irreducible spherical form   **
** at any chosen orientation.  The tensor principal axes will also be   **
** maintiained as well as a set of Euler angles that indicate the PAS   **
** orientation relative to some common axes.                            **
**                                                                      **
** Although not internal to this class, the interaction will blend with **
** a rank 2 spin tensor for the formation of an oriented CSA 		**
** Hamiltonian.                                                         **
**                                                                      **
** In addition to the functions contained in the base spatial tensor    **
** class, this IntCSA provides functions for building up the tensor and	**
** accessing the tensor elements from a "shift anisotropy" standpoint.	**
**                                                                      **
** The following defintions are used herein (Auv normlized by delzz):	**
**                                                                      **
*************************************************************************/

#ifndef   IntCSA_h_			// Is file already included?
#  define IntCSA_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Basics/ParamSet.h>		// Include parameter sets
#include <IntRank2/IntRank2A.h>		// Include base class headers
#include <IntRank2/IntRank2T.h>		// Include base class headers
#include <IntRank2/IntRank2.h>		// Include base class headers
#include <Matrix/matrix.h>		// Include GAMMA matrices
#include <Matrix/row_vector.h>		// Include GAMMA row vectors
#include <string>			// Include stdlibc++ strings
#include <fstream>			// Include stdlibc++ file streams

//forward declarations
class IntCSA;
MSVCDLL matrix HC0(double qn, double wCo, double eta=0, 
                                                 double theta=0, double phi=0);
MSVCDLL matrix HC1(double Om, double qn, double wCo, double eta,
                                              double theta=0.0, double phi=0.0);
class IntCSA: public IntRank2
  {
  double SISO;                          // Shift isotropy   (PPM)
  double SCSA;				// Shift anisotropy (PPM)
  double SOMEGA;			// Spin Larmor frequency (Hz)

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i              SHIFT ANISOTROPY INTERACTION ERROR HANDLING
// ____________________________________________________________________________
 
/* These functions handle errors for this class.  They route into the main
   error functons of GAMMA (Basics Module)

        // Input                SA	: CSA interaction (this)
        //			eidx	: Error index
        //                      noret   : Flag for return (0=return)
	//			pn	: String included in message
        // Output               none	: Error message
	// 				  Executaion stopped if fatal        */

         void ICerror(int eidx,                        int noret=0) const;
volatile void ICfatal(int eidx)                                     const;
         void ICerror(int eidx, const std::string& pn, int noret=0) const;
volatile void ICfatal(int eidx, const std::string& pn, int noret=0) const;

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

bool getCI(const ParameterSet& pset,
         double& Iqn, double& ppm, double& csa, double& eta, EAngles& EA,
                                    double& Om, int idx, bool warn=true) const;

bool getCI(const ParameterSet& pset, const Isotope& ISI,
                      double& ppm, double& csa, double& eta, EAngles& EA,
                                    double& Om, int idx, bool warn=true) const;

// ----------------------------------------------------------------------------
//            Get Isotropic Chemical Shift Value From A Parameter Set
// ----------------------------------------------------------------------------

/*         Input                SA      : CSA interaction (this)
                                pset    : A parameter set
                                idx     : Index value
				Om      : Larmor frequency (MHz)
                                warn    : Warning output flag
           Output               PPM     : Chemical shift anisotropy (PPM)
                                          constant from parameters in pset
           Note                         : This WILL NOT alter the interaction
           Note                         : Parameters are PPM, PPM(#), v, v(#)
           Note                         : We only accepted PPM here unless
           				  the Larmor frequency (Om) is known */

bool getPPM(const ParameterSet& pset, double& ppm,
                                int idx=-1, double Om=0, bool warn=true) const;
 
// ----------------------------------------------------------------------------
//           Get Shielding Anisotropy Value From A Parameter Set
// ----------------------------------------------------------------------------
 
 
bool getCSA(const ParameterSet& pset, double& csa,
                                             int idx=-1, bool warn=true) const;
 
        // Input                SA      : CSA interaction (this)
        //                      pset    : A parameter set
        //                      idx     : Index value
        //                      warn    : Warning output label
        //                                 0 = no warnings
        //                                 1 = warnings
        //                                >1 = fatal warnings
        // Output               CSA     : Chemical shift anisotropy (PPM)
        //                                value set from parameters in pset
   
// Look for "delzz" value. Currently allowed parameters are the following:
//
//  CSA,    CSA(#)      - Chemical Shift Anisotropy in PPM

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

           Input                SA      : CSA interaction (this)
                                pset    : A parameter set
                                idx     : Index value
                                warn    : Warning output flag
           Output               Om      : Spin Larmor frequency (MHz)
                                          from parameters in pset
           Note                         : This WILL NOT alter the interaction
           Note                         : Parameters are CSA, CSA(#)         */

bool getOm(const ParameterSet& pset,
               double& Om, const std::string& II, int idx=-1, bool warn=true) const;


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

bool setCI(                  const ParameterSet& p,int i=-1,int w=1);
bool setCI(const Isotope& II,const ParameterSet& p,int i=-1,int w=1);

// ----------------------------------------------------------------------------
//      Set Shift Anisotropy Interaction Spatial and Spin Tensor Components
// ----------------------------------------------------------------------------

// void setTs();        INHERITED
// void setAs();	INHERITED	 Set spatial tensor components

//void setTs(coord& B);
 
        // Input                SA	: CSA interaction (this)
        //                      B       : Oriented, normalized field vector
        // Output               none    : CSA interaction spherical-
        //                                spin components are generated
        // Note                         : No check is made to see if the
        //                                Tsph array has be made!

//                                            +
//                                         m  |
//                              T    = (-1)  T
//                               2,m          2,-m

// T  = Tsph[0]; T  = Tsph[1]; T   = Tsph[2]; T  = Tsph[3]; T   = Tsph[4]
//  2,0           2,1           2,-1           2,2           2,-2

//          1/2
//       [1]   [       1             ]          1                       1
// T   = |-|   |2I B - - (B I + B I )|  T   = - - [B I + B I ]   T    = - B I
//  2,0  [6]   [  z z  2   - +   + - ]   2,1    2   z +   + z     2,2   2  + +


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A          SHIFT ANISOTROPY INTERACTION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------
 
MSVCDLC IntCSA();
MSVCDLC IntCSA(const IntCSA &SA1);
 
// ----------------------------------------------------------------------------
//            Direct Constructors That Use Spherical Parameters
// ----------------------------------------------------------------------------

/* Here we try to use the set { Iqn, SIso, SCSA, Seta, SEAs , SOm }. The spin
   quantum number sets the dimension of the interaction spin Hilbert space.
   SIso is the isotropic chemical shift (PPM), SCSA is the chemical shift
   anisotropy (PPM), and Seta the shift asymmetry ([0,1]). The interaction
   orientation is taken from Euler angles SEAs and the external static field
   strength is specified by the spin Larmor frequency SOm (Hz). We will
   allow for default settings of the asymmetry, orientation, and Larmor
   frequency so that those arguments need not be input during construction.
   Their default values are all 0.                                           */

MSVCDLC IntCSA(const std::string& IsoI, double Siso, double SdelA,
                      double eta=0.0, const EAngles& EA=EAzero, double Om=0.0);
MSVCDLC IntCSA(const Isotope& IsoI,     double Siso, double SdelA,
                      double eta=0.0, const EAngles& EA=EAzero, double Om=0.0);
MSVCDLC IntCSA(double qn,               double Siso, double SdelA,
                      double eta=0.0, const EAngles& EA=EAzero, double Om=0.0);

// ----------------------------------------------------------------------------
//            Direct Constructors That Use Cartesian Parameters
// ----------------------------------------------------------------------------

/* Here we try to use the set { Iqn, Sxx, Syy, Szz, SEAs, SOm }. The spin
   quantum number sets the dimension of the interaction spin Hilbert space.
   The values of Sxx, Syy, and Szz are in PPM and used to determine the
   isotropic chemical shift SIso (PPM), chemical shift anisotropy SCSA (PPM),
   and shift asymmetry ETA ([0,1]). The interaction orientation is taken from
   Euler angles SEAs and the external static field strength specified by the
   spin Larmor frequency SOm (Hz). We will allow for default settings of the
   orientation and Larmor frequency so that those arguments need not be input
   during construction.  Their default values are 0.                         */

MSVCDLC IntCSA(const std::string& IsoI, const coord& SxSySz,
                                      const EAngles& EA=EAzero, double Om=0.0);
MSVCDLC IntCSA(const Isotope& IsoI,    const coord& SxSySz,
                                      const EAngles& EA=EAzero, double Om=0.0);
MSVCDLC IntCSA(double qn,              const coord& SxSySz,
                                      const EAngles& EA=EAzero, double Om=0.0);

// ----------------------------------------------------------------------------
//    Constructors Using Parameter Sets & Single Spin/Interaction Indices
// ----------------------------------------------------------------------------

/* These constructors use a single interaction (or spin) index and are ignorant
   about any particular spin type involved (other than its I value).  They'll
   try to read the parameters
 
                   { [Iso(i)/CI], CSA, Ctheta, Cphi, Ceta }
 
   where the shielding anisotropy is related to the spatial tensor delzz via

      ^          3                                     1
     / \ sigma = - del   = sigma  - sigma   = sigma  - - [ sigma  + sigma  ]
     ---         2    zz        ||        |        zz  2        xx       yy
                                         ---

   NOTE: This interaction maintains itself in PPM units!  That is somewhat odd
         in light of the fact that the class has no knowledge of applied field
         strengths nor knowledge of spin isotope types.  As a result, we cannot
         build CSA Hamiltonians in Hz units unless we know a conversion from
         PPM to Hz (i.e. an applied field strength.)

   An Example Shift Anisotropy Interaction Definition Are The ASCII Lines
 
   CI(1)        (1) : 1.5                  - Spin I quantum number
   CSA(1)       (1) : 40.0                 - Shift anisotropy (PPM)
   Cphi(1)      (1) : 45.0                 - Tensor orientation (deg)
   Ctheta(1)    (1) : 45.0                 - Tensor orientation (deg)
   Ceta(1)      (1) : 0.0                  - Tensor asymmetry [0,1]

   When idx!=-1 the parameter Iso(idx) will be used preferentially over the
   cooresponding parameter CI(idx) to set the interaction spin quantum value.
   The functions do allow for the suffix (#) but its not used in the default
   function call.  The functions do NOT allow for the prefix [#] since multiple
   interactions can be read using the suffix and multiple sets of interactions
   can be read using shift anisotropy interaction vectors (see class IntCSAVec)
   If the parameter CI is not present in the file/parameter set, CI=1/2 is used
   by default.  If Ceta is not present then Ceta is assumed to be zero.      */  
 
        // Input                SA      : CSA interaction (this)
        //                      pset    : Parameter set
        //                      I	: Isotope type
        //                      idx     : Interaction index (default -1->none)
        //                      warn    : Warning output label
        //                                 0 = no warnings
        //                                 1 = warnings
        //                                >1 = fatal warnings
        // Output               none    : Shift anisotropy interaction
        //                                constructed from parameters in pset
        // Note                         : Constructions that explicitly use
        //                                isotope labels are for support of
        //                                use by spin systems


MSVCDLC IntCSA(const                  ParameterSet& pset,int idx=-1,int wn=2);
MSVCDLC IntCSA(const Isotope& I,const ParameterSet& pset,int idx=-1,int wn=2);

// ----------------------------------------------------------------------------
//               Here Are The Assignment Operator & Destructor
// ----------------------------------------------------------------------------
 
MSVCDLL void operator= (const IntCSA &SA1);
MSVCDLC      ~IntCSA();

// ____________________________________________________________________________
// B                     SPATIAL TENSOR COMPONENT ACCESS
// ____________________________________________________________________________
 
//-----------------------------------------------------------------------------
//                              Isotropy Access
//-----------------------------------------------------------------------------

/* Here we only allow the chemical shift to be set in PPM. No assumptions are
   made regarding the external field strength, gyromagnetic ratio, or spin
   Larmor frequency.                                                         */

MSVCDLL double iso() const;
MSVCDLL void   iso(double ppm);
MSVCDLL double PPM() const;
MSVCDLL void   PPM(double ppm);

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

MSVCDLL double aniso() const;
MSVCDLL void   aniso(double csa);
MSVCDLL double CSA() const;
MSVCDLL void   CSA(double csa);

//-----------------------------------------------------------------------------
//			       Asymmetry Access
//-----------------------------------------------------------------------------

// Note that eta spans [0, 1] and is defined to be (Axx-Ayy)/Azz where
// |Azz| >= |Ayy| >= |Axx| in the treatment contained herein.
 
// double IntRank2A::eta( ) const;				INHERITED
// void   IntRank2A::eta(double Eta);				INHERITED
 
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

//-----------------------------------------------------------------------------
//               Cartesian Tensor Component Access (Not Normalized)
//-----------------------------------------------------------------------------

/* These allow one to access the Cartesian elements of the full (unscaled)
   shift anisotropy spatial tensor at any orientation.

                                    1/2
                            [ 6*PI ]
                      S   = | ---- |    * del   * A   + S
                       uv   [  5   ]         zz    uv    iso                 */


MSVCDLL double Sxx() const;
MSVCDLL double Syy() const;
MSVCDLL double Szz() const;
MSVCDLL double Sxy() const;
MSVCDLL double Syx() const;
MSVCDLL double Sxz() const;
MSVCDLL double Szx() const;
MSVCDLL double Syz() const;
MSVCDLL double Szy() const;

MSVCDLL double Sxx(double alpha, double beta, double gamma) const;
MSVCDLL double Syy(double alpha, double beta, double gamma) const;
MSVCDLL double Szz(double alpha, double beta, double gamma) const;
MSVCDLL double Syx(double alpha, double beta, double gamma) const;
MSVCDLL double Sxy(double alpha, double beta, double gamma) const;
MSVCDLL double Szx(double alpha, double beta, double gamma) const;
MSVCDLL double Szy(double alpha, double beta, double gamma) const;
MSVCDLL double Sxz(double alpha, double beta, double gamma) const;
MSVCDLL double Syz(double alpha, double beta, double gamma) const;

MSVCDLL double Sxx(const EAngles& EA) const;
MSVCDLL double Syy(const EAngles& EA) const;
MSVCDLL double Szz(const EAngles& EA) const;
MSVCDLL double Syx(const EAngles& EA) const;
MSVCDLL double Sxy(const EAngles& EA) const;
MSVCDLL double Szx(const EAngles& EA) const;
MSVCDLL double Szy(const EAngles& EA) const;
MSVCDLL double Sxz(const EAngles& EA) const;
MSVCDLL double Syz(const EAngles& EA) const;

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
              S   = | ---- |    * del   * A   + Kdel    * S
               uv   [  5   ]         zz    uv       u,v    iso

   where Kdel is a Kronecker delta function.                                 */

MSVCDLL matrix Smx() const;
 
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
  matrix IntRank2T::Tcomp(int comp) const;      // comp=[-2,2]	INHERITED 

  matrix IntRank2T::T20(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T21(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m1(     const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T22(      const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m2(     const vector<int>& HSs,int i,int j=-1) const;
  matrix IntRank2T::T2m(int m,const vector<int>& HSs,int i,int j=-1) const;  */
 
// ____________________________________________________________________________
// E          SHIFT ANISOTROPY INTERACTION CONSTANT ACCESS FUNCTIONS
// ____________________________________________________________________________
 
/*      (This is GAMMA Defined & Most Useful In Treatment of Liquids)

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

           Input                SA      : Shift anisotropy interaction
                                IsoI    : Spin isotope type
                                Om1H	: Proton field strength (Hz)
                                Om      : Larmor precession frequency (Hz)
           Return               xi      : CSA. interaction constant
           Note                           Value is in radians/sec            */


MSVCDLL double xiOm(const std::string& IsoI, double Om1H) const;
MSVCDLL double xiOm(double Om) const;

MSVCDLL double xi()          const;	// Overwrites inherited function
MSVCDLL void   xi(double x)  const;	// Overwrites inherited function

// ____________________________________________________________________________ 
// F                        PARAMETER SET FUNCTIONS 
// ____________________________________________________________________________
  
// ---------------------------------------------------------------------------- 
//  Functions To Make A Parameter Set From A Shielding Anisotropy Interaction 
// ---------------------------------------------------------------------------- 
  
/* This class has no implicit knowledge of the nuclear spin. The interaction
   may involve a 1H (I=1/2) or a 2H (I=1) or any GAMMA recognize nuclear spin.
   The preferred means of specifying a shielding anisotropy interaction when
   there is no knowledge of the nulcear spin isotope type is that which is used
   for filling up the parameter set.  The base parameters of individual
   interactions in that case are { CI(#), CSA(#), CTheta(#), CPhi(#), CEta(#) }.

           Input                SA	: Shift anisotropy interaction
                                idx     : Interaction index (default -1)
                                pfx     : Interaction 2nd indx (def -1)
           Output               pset    : Parameter set with only shielding
                                          anisotropy interaction parameters
 
     Function                                 Purpose
   ------------ 	-------------------------------------------------------
   ParameterSet		Convert interaction into a parameter set
   operator +=          Adds interaction to existing parameter set
   PSetAdd              Adds interaction to existing parameter set, this
                        allows for an index and a prefix in the parameters   */

MSVCDLL operator ParameterSet( ) const;
MSVCDLL friend void operator+= (ParameterSet& pset, const IntCSA &SA);
MSVCDLL        void PSetAdd(ParameterSet& pset, int idx=-1, int pfx=-1) const;
    
// ----------------------------------------------------------------------------
//    Functions To Output S.A. Interaction To ASCII From A Parameter Set
// ----------------------------------------------------------------------------
 
        // Input                SA	: Shift anisotropy interaction
        //                      fn	: Output file name
        //                      ofstr   : Output file stream
        //                      idx     : Interaction index (default -1)
        //                      pfx     : Interaction 2nd indx (def -1)
        //                      warn    : Warning level
        // Output               none    : Shift anisotropy interaction is
        //                                written as a parameter set to
        //                                file filename or output stream ofstr

MSVCDLL int write(const std::string &fn,int idx=-1,int pfx=-1,int wrn=2) const;
MSVCDLL int write(std::ofstream& ofstr, int idx=-1,int pfx=-1,int wrn=2) const;
 
// ____________________________________________________________________________
// G              SHIFT ANISOTROPY INTERACTION INPUT FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//      Direct Read of SA Intraction From An ASCII File Or A Parameter Set
// ----------------------------------------------------------------------------

/* These next two read functions utilize a single spin/interaction index, for
   the spin involved in the interaction.  They'll try to read the parameter set

                  { [Iso(i)/CI] , CSA, Ctheta, Cphi, Ceta }

   The functions do NOT allow for the suffix (#), since it is used by the spin
   index anyway. The functions do NOT allow for the prefix [#] either, because
   multiple shift anisotropy interactions can be defined in the same file by 
   switching spin index.  Multiple sets of interactions can be read using SA
   interaction vectors (see class IntCSAVec.)                                */
 

 
        // Input                SA      : CSA interaction (this)
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
	//				  or contained in pset

MSVCDLL bool read(const std::string &filename, int idx=-1, int warn=2);
MSVCDLL bool read(const ParameterSet& pset,    int idx=-1, int warn=2);

// ----------------------------------------------------------------------------
//                    Interactive Ask/Read ASCII File
// ----------------------------------------------------------------------------

MSVCDLL std::string ask_read(int argc, char* argv[], int argn, int idx=-1);

        // Input                SA      : CSA interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
        //                      idx     : Interaction index
        // Output               String  : The parameter argn of array argc is
        //                                used to supply a filename from which
        //                                the interaction is read

 
// ----------------------------------------------------------------------------
//               Interactive Ask For All Kinds Interaction Info
// ----------------------------------------------------------------------------

MSVCDLL void ask(int argc, char* argv[], int& qn, double& CI,
         double& Csa, double& Ceta, double& Ctheta, double& Cphi, int Cflag=1);

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


MSVCDLL void askset(int argc, char* argv[], int& qn, int SAflag=0);
 
        // Input                SA      : CSA interaction (this)
        //                      argc    : Number of arguments
        //                      argv    : Array of arguments
        //                      qn      : Query value
        //                      SAflag  : Flag if delzz or CSA requested
        // Output               none    : SA is set interactively
        // Note                         : This is INTERACTIVE!
 

MSVCDLL void askset(int SAflag=0);

        // Input                SA      : CSA interaction (this)
        //                      SAflag  : Flag is delzz or CSA requested
        // Output               none    : SA is set interactively
        // Note                         : This is INTERACTIVE!

// ____________________________________________________________________________
// I                          OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//   Functions That Generate Simple Strings To Simplify & Modularize Printing
//-----------------------------------------------------------------------------

MSVCDLL std::string ShiftString()  const;
MSVCDLL std::string CSAString()    const;
MSVCDLL std::string LarmorString() const;

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

MSVCDLL std::vector<std::string> InfoStrings() const;

//-----------------------------------------------------------------------------
//  Functions That Generate String Arrays to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

/* vector<string> IntRank2T::TStrings(int M) const;                INHERITED    
 
                                [ x.x, x.x, x.x] 
                        T     = [ x.x, x.x, x.x] 
                         2,m    [ x.x, x.x, x.x] 
                                [ x.x, x.x, x.x] 
 
  where M = { 0, 1, ..., 4 } ===> m = { 0, 1, -1, 2, -2 }                    */   
 


        // Input                SA      : CSA interaction (this)
        // Output               CSS     : Pointer to array of 5 strings
        // Note                         : The String array must be deleted
        //                                outside of this routine!

//-----------------------------------------------------------------------------
//     Functions That Generate Ouput Of The Shift Anistropy Interaction
//-----------------------------------------------------------------------------
 
  
  
        // Input                D       : Dipolar interaction (this) 
        //                      ostr    : Output stream 
        //                      fflag   : Format flag 
        //                                  0 - Basic Parameters
        //                                 !0 - Full output 
        // Output               none    : CSA spatial tensor parameters
        //                                placed into the output stream
        // Note                         : This does NOT use the base class 
        //                                virtual overload because we write 
        //                                out two spin quantum values here?

MSVCDLL        std::ostream& print(std::ostream& ostr, int fflag=0) const;
MSVCDLL friend std::ostream& operator<<   (std::ostream& osr,  const IntCSA& C);





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


        // Input                out      : Output stream;
        //                      C	: Shift anisotropy tensor to write
        // Output                        : Modifies output stream


MSVCDLL std::ostream& printAT(std::ostream& ostr) const;
MSVCDLL std::ostream& printSpherical(std::ostream& ostr);
MSVCDLL std::ostream& printCartesian(std::ostream& ostr);
MSVCDLL std::ostream& printCartesian(std::ostream& ostr, double theta, double phi=0);
MSVCDLL static std::ostream& STList(std::ostream& ostr, int fflag=0);


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



// ____________________________________________________________________________
// D                   SHIFT ANISOTROPY FREQUENCY FUNCTIONS
// ____________________________________________________________________________

// This frequency will be the splitting between the 2*I transitions contained   
// in a CSA Hamiltonian.  This will be when the shift anisotropy is weak
// relative to the Zeeman interaction (high field approximation, first order
// terms only) if the tensor is oriented in it's principal axes (PAS).  If the
// tensor isn't aligned in it's PAS, the splitting will vary with orientation
// according to
 
// Also, keep in mind that the Euler angles {phi,theta,gamma} which are kept
// with the tensor are used to relate the tensor PAS to some (common) set of
// coordinate axes.  They are not the phi and theta used in the above
// formula (unless you which the tensor aligned in the common axis system)


MSVCDLL double wC( ) const;

        // Input                SA	: Shift anisotropy interaction
        // Return               wC      : Shift anisotropy frequency (Hz)


MSVCDLL void wC(double W);

        // Input                SA	: Shift anisotropy interaction
        //                      W       : Shift anisotropy frequency (Hz)
        // Return               void    : The shift anisotropy frequency
        //                                of the interaction is set
        // Note                         : The interaction I value
        //                                must be set prior to this


// ____________________________________________________________________________
// L                    SHIFT ANISOTROPY HAMILTONIAN FUNCTIONS
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
   matrices. Their units will be Hz.					     */

// ----------------------------------------------------------------------------
//               First Order Shift Anisotropy Interaction Hamiltonian
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
// Applicable When The Shift Anisotropy Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------

/*  The secular part of the quadrupolar Hamiltonian is that which returns
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

MSVCDLL matrix H0()                             const;
MSVCDLL matrix H0(double A, double B, double G) const;
MSVCDLL matrix H0(const EAngles& EA)            const;

MSVCDLL matrix H0(const std::vector<int>& HSs, int i)                               const;
MSVCDLL matrix H0(const std::vector<int>& HSs, int i, double A, double B, double G) const;
MSVCDLL matrix H0(const std::vector<int>& HSs, int i, const EAngles& EA)            const;

// ----------------------------------------------------------------------------
//              Second Order Shift anisotropy Interaction Hamiltonians
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
// Applicable When The Shift anisotropy Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------

//matrix H1(double Om) const;
//matrix H1(double Om, double theta, double phi=0) const;

        // Input                Om      : Field Strength (Larmor in Hz)
        // Output               HC1     : The 2nd order secular part of the
        //                                quadrupolar Hamiltonian (default
        //                                basis, Hz) for the interaction
        // Note                         : Also called the 1st order quadrupolar
        //                                interaction (perturbation theory)
 
//  [2]              -1   [ 1    ]2 [                          2     2   
// H   (theta,phi) = -- * | - w  |  | 2*A A  (theta,phi)*Iz*(4*I - 8Iz - 1)
//  SA               Om   [ 3  SA]  [    1 -1
//
//                                                              2     2     ]
//                                  + 2*A A  (theta,phi)*Iz*(2*I - 2Iz - 1) |
//                                       2 -2                               ]

 
// ----------------------------------------------------------------------------
//      Summed First & Second Order Shift anisotropy Interaction Hamiltonians
//       These are SECULAR (Rotationally Invariant About Bo Field Axis)
//  These Are For When The Shift anisotropy Interaction Is A Perturbation To Zeeman
// ----------------------------------------------------------------------------

//matrix Hw(double Om) const;
//matrix Hw(double Om, double theta, double phi) const;

        // Input                sys  : Spin system
        //                      wSA  : Shift anisotropy frequency (Hz)
        //                      i    : Spin index
        // Output               HCw  : The secular part of the quadrupolar
        //                             Hamiltonian (default basis, Hz)
        //                             for the spin i
        // Note                      : No asymmetry is considered here
        // Note                      : This is the sum of the 1st & 2nd order
        //                             quadrupolar interactions (pert. theory)
        // Note                      : This is rotationally invariant about z
 
//  This function returns only the secular part of the second order quadrupolar
//  Hamiltonian (from perturbation theory).  Note that this still assumes that
//  the quadrupolar interaction is a perturbation to to the Zeeman interaction.
 
//                                (0)    (1)      [1]    [2]
//                       H    =  H    + H     =  H    + H
//                        SA      SA     SA       SA     C


// ----------------------------------------------------------------------------
//                 Full Shift anisotropy Interaction Hamiltonians
//       These Have No Approximations & Are Independent of Field Strength
// ----------------------------------------------------------------------------
 

//matrix H( ) const;
 
        // Input                SA      : Shift anisotropy interaction
        // Note                         : This will return in the spin Hilbert
        //                                space of dimension 2I+1

//                                m     SA   SA                C
//            H (theta,phi) = (-1)  * Xi  * A   (theta,phi) * T
//             SA                            2,m               2,-m


// ____________________________________________________________________________
//                   SHIFT ANISOTROPY HAMILTONIAN FRIEND FUNCTIONS
// ____________________________________________________________________________

// These Don't Depend on Class IntCSA Variables At All.  The Functions Just
// Mimic The Class' Hamiltonian Generation.


MSVCDLL friend matrix HC0(double qn, double wCo, double eta, 
                                                 double theta, double phi);
 
        // Input                qn      : Quantum number (1, 1.5, 2.5,...)
        //                      wCo     : PAS Shift anisotropy frequency
        //                      eta     : Shift anisotropy asymmetry [0,1]
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
//            SA              6  SA               z             C

// where                                                          
 
//                      [ 1      2               1        2                   ]  
// w (theta,phi) = W    | - [3cos (theta) - 1] + - eta sin (theta)*cos(2*phi) |
//  SA              C,o [ 2                      2                            ]  
 
// and
//                                    3*CCC
//                            w    = --------
//                             C,o   2I(2I-1)



MSVCDLL friend matrix HC1(double Om, double qn, double wCo, double eta,
                                              double theta, double phi);
                                                        
        // Input                Om      : Field Strength (Larmor in Hz)         
        //                      qn      : Quantum number (1, 1.5, 2.5,...)
        //                      wCo     : PAS Shift anisotropy frequency
        //                      eta     : Shift anisotropy asymmetry [0,1]
        //                      theta   : Orientation angle (degrees)
        //                      phi     : Orientation angle (degrees)
        // Output               HC1     : The 2nd order secular part of the
        //                                quadrupolar Hamiltonian (default
        //                                basis, Hz) for the interaction
        //                                for orientation {phi, theta}
        //                                from the tensor PAS
        // Note                         : Also called the 1st order quadrupolar
        //                                interaction (perturbation theory)
        // Note                         : Rotationally invariant about z
        // Note                         : This will return in the spin Hilbert
        //                                space of dimension 2I+1

  };

#endif								// IntCSA.h
