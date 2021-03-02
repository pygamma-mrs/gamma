/* IntRank2A.h **************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Irreducible Rank 2 Spatial Tensor	Interface		**
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
** Class IntRank2A represents a primitive rank 2 interaction. The class	**
** contains the spatial part of a GAMMA scaled rank 2 interaction which **
** is restricted to have the following 5 properties:                    **
**                                                                      **
**        1.) rank=2        2.) traceless        3.) symmetric 		**
**                                                                      **
**          4.) eta = [0, 1]     5.) delz = sqrt[5/(6*PI)]		**
**                                                                      **
** The interaction will be stored as a spatial tensor that is in        **
** irreducible spherical form in its principal axes (PAS). It can 	**
** produce any or all components set to a specified orientation or      **
** Euler rotation from the PAS.						**
**                                                                      **
** Note: No rank zero terms (isotropic components) are maintained!      **
**       No rank one temrms (asymmetric components) are maintained!     **
**                                                                      **
** The rank 2 spatial tensor provided by this class is interaction      **
** independent and scaled so that the spherical components are directly **
** associated to normalized spherical harmonics.  The class by-passes   **
** GAMMA's generic spatial tensors (class space_T) for much improved    **
** computational efficiency.                                            **
**                                                                      **
** The following defintions are used herein (Auv GAMMA normlized):      **
**                                                                      **
** 1.)               |Azz| >= |Ayy| >= |Axx|                            **
**                                                                      **
** 2.)    Azz=2C, eta=(Axx-Ayy)/Azz, Axx=C(eta-1), Ayy=-C(1+eta)        **
**                                                                      **
**                                    1/2                               **
**                           [   5   ]                                  **
**     where             C = | ----- |                                  **
**                           [ 24*PI ]                                  **
**                                                                      **
**  Since these components are normalized, only eta is required to 	**
**  specify the rank 2 irreducible spherical tensor. From this alone	**
**  all tensor components can be generated for any orientation.	        **
**                                                                      **
*************************************************************************/
 
#ifndef   IntRank2A_h_			// Is file already included?
#  define IntRank2A_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface
#  endif

extern const double RT6PIO5; 		// Needed constant sqrt[6*PI/5]
extern const double RT5O4PI;		// Needed constant sqrt[5/(4*PI)];
extern const double RT5O24PI;		// Needed constant sqrt[5/(24*PI)];

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Matrix/matrix.h>		// Include GAMMA matrices
#include <Matrix/complex.h>		// Include GAMMA complex numbers
#include <Matrix/row_vector.h>		// Include GAMMA row vectors
#include <Level1/coord.h>		// Include coordinates
#include <Level2/EAngles.h>		// Inlcude Euler angles
#include <IntRank2/IntRank2ACmp.h>	// Inlcude rank 2 spatial cmpnts
#include <string>			// Include libstdc++ strings
#include <fstream>			// Include libstdc++ file streams
#include <iostream>			// Include libstdc++ file streams
#include <vector>			// Include libstdc++ STL vectors

//forward declarations
class IntRank2A;
MSVCDLL void GetEtaDelzz(const matrix& A, 
            double& Aeta, double& Adelz, double ecut=1.e-5, double dcut=1.e-9);

class IntRank2A
  {
  double ETA;				// Spatial tensor eta [0, 1]
  EAngles _EAs;                         // Euler angles (orient, radians)

  friend class IntRank2;		// Allow IntRank2 Full Access
  friend class IntCSA;			// Allow IntCSA   Full Access
  friend class IntDip;			// Allow IntDip   Full Access
  friend class IntQuad;			// Allow IntQuad  Full Access
  friend class IntHF;			// Allow IntHF    Full Access
  friend class IntG;			// Allow IntG     Full Access
   

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

         void IR2Aerror(int eidx,                       int noret=0) const;
volatile void IR2Afatal(int eidx)                                    const;
         void IR2Aerror(int eidx, const std::string pn, int noret=0) const;

// ____________________________________________________________________________
// ii         RANK 2 SPATIAL TENSOR FROM PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

/* In this section we specify how a rank 2 spatial tensor relates to a GAMMA
   parameter file.  That is, these functions determine how a spatial tensor may
   be read/set from an external ASCII file. These are private because some
   parameters are inter-related and hence must be set in a concerted manner.

   Since this class embodies a normalized irreducible rank 2 spatial tensor,
   only the asymmetry parameter is required to define it. Granted, users will
   rarely use this function directly in favor of using either Cartesian or 
   spherical means to specify reducible rank 2 spatial tensor. We will allow
   for setting the eta value in thos eterms, however that is a job best suited
   for specific interactions rather than here so one may account for any
   isotropic terms, rank 1 terms, and scaling when reading in the tensor. 

   We also read in the tensor orientation. This is normally done using a set
   of Euler angles that relate the tensor PAS to the laboratory frame. But note
   that GAMMA has composite rotations which generalize reference frames. Thus
   one can readily us these angles to relate to some other frame and/or use
   composite rotations to switch between any number of coordinate axes.      */    

// ----------------------------------------------------------------------------
//         Functions To Read The Full Irreducible Rank 2 Spatial Tensor
// ----------------------------------------------------------------------------

/* These functions are used to read the entire irreducible rank 2 spatial
   tensor. To deal with the various parameters we allow for specifying the
   tensor the following logic hierarchy will be used:

   1.) Look for { Axy, Axz, Ayz }. If found, look for { Axx, Ayy, Azz } and
       then parse the 3x3 Cartesian array for both the asymmetry and the
       tensor orientation. 
   2.) Look for { Axx, Ayy, Azz }. If found look for the orientation.
   3.) Look for { eta }. If found look for the orientation.

   The hard part is deciding when we allow any parameters to not be set. For
   example, Azz is redundent with { Axx, Ayy } for an irreducible tensor in
   the PAS. Also users may just neglect setting parameters such as eta,
   Ayz, or EAngles because they want them left at zero.                      */

bool getA(const ParameterSet& pset, const std::string& A,
       double& Aeta, EAngles& EA, int idxI, int idxS=-1, bool warn=true) const;

bool getAX(const ParameterSet& pset, const std::string& A,
           double& Aeta, EAngles& EA, int idxI, int idxS=-1, int warn=2) const;

bool getACart(const ParameterSet& pset, const std::string& A,
        coord& Aize, EAngles& EA, int idxI, int idxS=-1, bool warn=true) const;

// ----------------------------------------------------------------------------
//                   Functions To Get The Asymmetry (eta)
// ----------------------------------------------------------------------------

/* These functions are used to obtain the tensor asymmetry value eta from
   a parameter set.  Note that the parameter read, "Pbase" is generic so that
   derived classes may supply their own parameter acceptable parameter names.
   The eta values are unitless and restricted between the range [0, 1],  and
   we're happy if no eta value is specified, just switch the warnings off if
   one will often just not set eta. If eta is not read it's just set to zero.

           Input                IR2A    : Rank 2 spatial tensor (this)
                                pset    : A parameter set
                                A       : Parameter base name
				Aeta    : For return eta value
                                idxI    : Index value of first spin
                                idxS    : Index value of 2nd spin
                                warn    : Warning output switch
           Output               void    : Tensor asymmetry, eta,
                                          obtained from parameters in pset   */


bool getAeta(const ParameterSet& pset, const std::string& A,
                    double& Aeta, int idxI, int idxS=-1, bool warn=true) const;


// ----------------------------------------------------------------------------
//         Functions To Read Interaction Orientation (Euler Angles)
// ----------------------------------------------------------------------------

/* Each interaction stores one set of Euler Angles {alpha, beta, gamma} that
   describe its orientation with respect to its principal axis system (PAS).
   Typically these angles are those that relate the PAS to the laboratory
   frame, but they may easily be used to relate the PAS to some other frame
   (molecular, diffusion, ......). However, GAMMA contains composite
   rotations well suited for employing any number of reference frames.       */

bool getOrientation(const ParameterSet& pset, const std::string& Pbase,
                  EAngles& EA, int idxI=-1, int idxS=-1, bool warn=true) const;

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

int getAxAyAz(const   ParameterSet& pset, const std::string& A,
                  coord& Axyz, int idxI=-1, int idxS=-1, bool warn=true) const;
int getAoffdiag(const ParameterSet& pset, const std::string& A,
                  coord& Aod,  int idxI=-1, int idxS=-1, bool warn=true) const;

// ----------------------------------------------------------------------------
//                 Functions To Read Int Two Spin Coordinates
// ----------------------------------------------------------------------------

/* Spin pair interactions, such as dipolar and hyperfine interactions, may have
   their orientation (and perhaps their strength) set by a pair of "spin"
   coordinates. Here we read in two coordinates, one for each spin index and
   return them. Failure will occur if either spin coordinate is missing.     */

bool getCoords(const ParameterSet& pset,
             coord& ptI, coord& ptS, int idxI, int idxS, bool warn=true) const;

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A             RANK 2 INTERACTION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

MSVCDLC IntRank2A();
MSVCDLC IntRank2A(const IntRank2A &IR2Ab);

// ----------------------------------------------------------------------------
//         Direct Constructors Using Cartesian Spatial Tensor Components
// ----------------------------------------------------------------------------

/* These constructors require values for Axx, Ayy, and Azz. The Cartesian
   values are used only to set the spatial tensor asymmetry.  Note that since
   we are dealing with an irreducible spatial tensor the three Cartesian values
   MUST sum up to zero (tensor is traceless), so any isotropic value will be
   automatically removed where Aiso = (Axx+Ayy+Azz)/3. 

           Input                IR2A     : Rank 2 spatial tensor (this)
                                AxAyAz  : Cartesian PAS values
			OR
                                pt1     : Component 1 position
                                pt2     : Component 2 position
                                eta     : Tensor asymmetry

           Output               none    : Rank 2 spatial tensor constructed
           Note                         : Theta PAS z down, phi PAS x over
           Note                         : Insist |Azz| >= |Ayy| >= |Axx|     */

MSVCDLC IntRank2A(const coord& AxAyAz, const EAngles& EA=EAzero);
 
// ----------------------------------------------------------------------------
//         Direct Constructors Using Spherical Spatial Tensor Components
// ----------------------------------------------------------------------------

/* Since we are dealing with a normalized irreducible tensor, the only value
   that defines it (in the PAS) will be the asymmetry, eta.  The delzz value
   is automatically set to the normalized value of sqrt[5/(6*PI)].           */

MSVCDLC IntRank2A(double eta, const EAngles& EA=EAzero);

// ----------------------------------------------------------------------------
//                       Constructors Using Parameter Sets
// ----------------------------------------------------------------------------

/* In this case the constructon is performed entirely from parameters found in
   the specified parameter set.  The index allows users to read in multiple
   spatial tensors by setting a (#) suffix on all parameter names, #=idx
 
           Input                IR2A	: Rank 2 spatial tensor
                                pset    : Parameter set
                                idx     : Interaction index (default -1->none)
           Output               none    : Rank 2 spatial tensor constructed
                                          from parameters in pset             */

MSVCDLC IntRank2A(ParameterSet& pset, int idxI=-1, int idxS=-1, int warn=2);

// ----------------------------------------------------------------------------
//                         Assignment and Destruction
// ----------------------------------------------------------------------------

/* These must be virtual because our plan is to derive a class for generic
   rank 2 interactions (IntRank2) using this as its base class.               */

MSVCDLL virtual void operator= (const IntRank2A &IR2Ab);
MSVCDLC virtual      ~IntRank2A();

// ____________________________________________________________________________
// B                RANK 2 ANISOTROPY & ASYMMETRY FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//                               Anisotropy Access
//-----------------------------------------------------------------------------

/*                  ^
   The anisotropy, /_\ A, relates to the delzz value of A. Because GAMMA uses
   irreducible rank 2 spatial tensors as a base class, both the anisotropy and
   the delzz value are constant for all spatial tensors.

                       1/2
                  [ 5 ]                  ^      3
          del   = |---]  = 0.51503      /_\ A = - del   = 0.77255
             zz   [6PI]                         2    zz                      */

MSVCDLL static double delzz();
MSVCDLL static double delA();

//-----------------------------------------------------------------------------
//                             Asymmetry Access
//-----------------------------------------------------------------------------

/* These functions provide access to the spatial tensor asymmetry.  Remember
   that since this class is a normalized irreducible rank 2 spatial tensor the
   only value that defines it (in the PAS) is the assymetry value, eta. We
   restrict eta to span [0, 1] by using the definition eta = (Axx-Ayy)/Azz
   with |Azz| >= |Ayy| >= |Axx| in the treatment contained herein.

	   Input		IR2A	: Rank 2 spatial tensor
	   			Eta	: Rank 2 spatial tensor asymmetry
	   Output		eta 	: Return asymmetry of spatial tensor
	   			void    : Asymmetry value is set to Eta      */

MSVCDLL double eta( ) const;
MSVCDLL void   eta(double Eta);

// ____________________________________________________________________________
// C            RANK 2 TENSOR SPHERICAL COMPONENT ACCESS FUNCTIONS
// ____________________________________________________________________________
 
/* These allow one to access the irreducible spherical elements of the tensor
   either in the PAS or at a specific orientation without rotating the entire
   tensor. In these functions theta is the angle down from the PAS z axis & 
   phi the angle over from PAS x axis. The angles alpha, beta and gamma are 
   the three Euler angles
 
           Input                IR2A	: Rank 2 spatial tensor 
                                alpha	: Orientation angle (radians)
                                beta	: Orientation angle (radians)
                                gamma	: Orientation angle (radians)
				EA      : Orientation angles (radians)
           Return               A       : Rank 2 spherical element (Hz)
                                 2,m      for orientation {theta, phi}  
                                          or {alpha, beta, gamma} from 
	  				  the tensor PAS                     */

MSVCDLL static complex A20PAS();
MSVCDLL static complex A21PAS();
MSVCDLL static complex A2m1PAS();
MSVCDLL static complex A22PAS();
MSVCDLL static complex A2m2PAS();

MSVCDLL complex A20()  const;
MSVCDLL complex A21()  const;
MSVCDLL complex A2m1() const;
MSVCDLL complex A22()  const;
MSVCDLL complex A2m2() const;

MSVCDLL complex A20(double  alpha, double beta, double gamma) const;
MSVCDLL complex A21(double  alpha, double beta, double gamma) const;
MSVCDLL complex A2m1(double alpha, double beta, double gamma) const;
MSVCDLL complex A22(double  alpha, double beta, double gamma) const;
MSVCDLL complex A2m2(double alpha, double beta, double gamma) const;

MSVCDLL complex A20(const  EAngles& EA) const;
MSVCDLL complex A21(const  EAngles& EA) const;
MSVCDLL complex A2m1(const EAngles& EA) const;
MSVCDLL complex A22(const  EAngles& EA) const;
MSVCDLL complex A2m2(const EAngles& EA) const;

MSVCDLL complex A2m(int m)                                          const;
MSVCDLL complex A2m(int m, double alpha, double beta, double gamma) const;
MSVCDLL complex A2m(int m, const EAngles& EA)                       const;
 

MSVCDLL IR2ASph SphCmpPAS()                                     const;
MSVCDLL IR2ASph SphCmp()                                        const;
MSVCDLL IR2ASph SphCmp(double alpha, double beta, double gamma) const;
MSVCDLL IR2ASph SphCmp(const EAngles& EA)                       const;

/* This function returns all 5 spherical PAS components but it uses
   the angular momentum indexing scheme m = {-2,-1,0,1,2}.                   */

MSVCDLL complex AcompPAS(int comp) const;
MSVCDLL complex Acomp(int    comp) const;

/*                         1/2
                    [  5  ]          2                     2
   A  (theta,phi) = |-----| * [ 3*cos (theta) - 1 + eta*sin (theta)cos(2*phi) ]
    20              [16*PI]
        

                             1/2                       1/2
                       [  5 ]                     [ 3 ] 
           A   (PAS) = [----]          A   (EA) = | - |  A  (PAS)
            2,0        [4*PI]           2,0       [ 2 ]   zz                  */



/*                           1/2                                     
                      [  5  ]              [
   A   (theta,phi)  = |-----| * sin(theta)*| 3*cos(theta)
    2,1               [24*PI]              [
                                                                            ]
                                - eta*[cos(theta)*cos(2*phi) - i*sin(2*phi)]|
                                                                            ]

                                                  
           A   (PAS) = 0           A   (EA) = - A  (PAS) + i A  (PAS)
            2,1                     2,1          xz           yz              */

/*                              1/2                                 
                         [  5  ]                [
   A    (theta,phi)  = - |-----| * sin(theta) * | 3*cos(theta)
    2,-1                 [24*PI]                [
                                                                              ]
                                 - eta*[cos(theta)*cos(2*phi) + i*sin(2*phi)] |
                                                                              ]
                                                  
          A    (PAS) = 0          A    (EA) =  A  (PAS) - i A  (PAS)
           2,-1                    2,-1         xz           yz              */

 
/*                          1/2
                     [  5  ]   [      2                  2
   A   (theta,phi) = |-----|   [ 3*sin (theta) + eta*[cos (theta)+1]*cos(2*phi)
    2,2              [96*PI]   [
                                                                              ]
                                               - 2i*eta*cos(theta)*sin(2*phi) ]
                                                                              ]

 
        // Note                         : This is sqrt[5/(24*PI)]eta in PAS
	//				: This should be (Axx-Ayy)/2-iAxy

                            1/2
                     [  5  ]   [     2                  2
  A    (theta,phi) = |-----|   [3*sin (theta) + eta*[cos (theta)+1]*cos(2*phi)
   2,-2              [96*PI]   [
                                                                           ]
                                            + 2i*eta*cos(theta)*sin(2*phi) ]
                                                                           ] */


// ____________________________________________________________________________
// D            RANK 2 TENSOR CARTESIAN COMPONENT ACCESS FUNCTIONS
// ____________________________________________________________________________
 
/* These allow one to access the Cartesian elements at a specific orientation
   without rotating the entire tensor.  In these functions theta is the angle
   down from the PAS z axis and phi the angle over from the PAS x axis.

        // Input                IR2A	: Rank 2 spatial tensor 
        //                      theta   : Orientation angle (radians)
        //                      phi     : Orientation angle (radians)
        //                      alpha	: Orientation angle (radians)
        //                      beta	: Orientation angle (radians)
        //                      gamma	: Orientation angle (radians)
        // Return               A       : Rank 2 Cartesian element
        //                       uv       for orientation {theta, phi}  
        //                                or {alpha, beta, gamma} from 
	//				  the tensor PAS                     */

MSVCDLL double AxxPAS() const;
MSVCDLL double AxyPAS() const;
MSVCDLL double AxzPAS() const;
MSVCDLL double AyxPAS() const;
MSVCDLL double AyyPAS() const;
MSVCDLL double AyzPAS() const;
MSVCDLL double AzxPAS() const;
MSVCDLL double AzyPAS() const;
MSVCDLL double AzzPAS() const;

MSVCDLL double Axx() const;
MSVCDLL double Axy() const;
MSVCDLL double Axz() const;
MSVCDLL double Ayx() const;
MSVCDLL double Ayy() const;
MSVCDLL double Ayz() const;
MSVCDLL double Azx() const;
MSVCDLL double Azy() const;
MSVCDLL double Azz() const;

MSVCDLL double Axx(double alpha, double beta, double gamma) const;
MSVCDLL double Axy(double alpha, double beta, double gamma) const;
MSVCDLL double Axz(double alpha, double beta, double gamma) const;
MSVCDLL double Ayx(double alpha, double beta, double gamma) const;
MSVCDLL double Ayy(double alpha, double beta, double gamma) const;
MSVCDLL double Ayz(double alpha, double beta, double gamma) const;
MSVCDLL double Azx(double alpha, double beta, double gamma) const;
MSVCDLL double Azy(double alpha, double beta, double gamma) const;
MSVCDLL double Azz(double alpha, double beta, double gamma) const;

MSVCDLL double Axx(const EAngles& EA) const;
MSVCDLL double Axy(const EAngles& EA) const;
MSVCDLL double Axz(const EAngles& EA) const;
MSVCDLL double Ayx(const EAngles& EA) const;
MSVCDLL double Ayy(const EAngles& EA) const;
MSVCDLL double Ayz(const EAngles& EA) const;
MSVCDLL double Azx(const EAngles& EA) const;
MSVCDLL double Azy(const EAngles& EA) const;
MSVCDLL double Azz(const EAngles& EA) const;

MSVCDLL row_vector CartCompsPAS()                          const;
MSVCDLL row_vector CartComps()                             const;
MSVCDLL row_vector CartComps(double A, double B, double G) const;
MSVCDLL row_vector CartComps(const EAngles& EA)            const;

MSVCDLL IR2ACart CartCmpPAS()                          const;
MSVCDLL IR2ACart CartCmp()                             const;
MSVCDLL IR2ACart CartCmp(double A, double B, double G) const;
MSVCDLL IR2ACart CartCmp(const EAngles& EA)            const;
 
/* ----------------------------------------------------------------------------          

                          1/2
                   [  5  ]          2                     2
  A  (theta,phi) = |-----| * [ 3*sin (theta) - 1 + eta*cos (theta)cos(2*phi) ]
   xx              [24*PI]
 
                     1/2
              [  5  ]                     1 [             ]    -1/2
   A  (PAS) = |-----|           A  (EA) = - | A   + A     | - 6     * A 
    xx        [24*PI]            xx       2 [  2,2   2,-2 ]            2,0  
                  
   ----------------------------------------------------------------------------          

                             1/2
                      [  5  ]
   A  (theta,phi) = - |-----| * [ 1 + eta*cos(2*phi) ]
    yy                [24*PI]

                       1/2
                [  5  ]                        1 [             ]    -1/2
   A  (PAS) = - |-----| (1 + eta)  A  (EA) = - - | A   + A     | - 6     * A 
    yy          [24*PI]             yy         2 [  2,2   2,-2 ]            2,0  
                  

   ----------------------------------------------------------------------------          

                           1/2
                    [  5  ]          2                     2
   A  (theta,phi) = |-----| * [ 3*cos (theta) - 1 + eta*sin (theta)cos(2*phi) ]
    zz              [24*PI]
 
                                  1/2                     1/2
                            [ 5  ]                     [3] 
                 A  (PAS) = |----|           A  (EA) = |-|  * A 
                  zz        [6*PI]            zz       [2]     2,0  
                  
   ----------------------------------------------------------------------------          

                             1/2                                             
                      [  5  ]
   A  (theta,phi) = - |-----| * [ eta*cos(theta)sin(2*phi) ] = A
    xy                [24*PI]                                   yx

                                            [i]    [                      ]
    A  (PAS) = 0                A  (EA) = - |-|  * | A   (EA) - A    (EA) |
     xy                          xy         [2]    [  2,2        2,-2     ] 


   ----------------------------------------------------------------------------          

                             1/2                           
                      [  5  ]
   A  (theta,phi) = - |-----| * sin(theta)*cos(theta)*[3-eta*cos(2*phi)] = A
    xz                [24*PI]                                               zx
 
                                            [1]    [                      ]
    A  (PAS) = 0                A  (EA) = - |-|  * | A   (EA) - A    (EA) |
     xz                          xz         [2]    [  2,1        2,-1     ] 


   ----------------------------------------------------------------------------          

                                              1/2                           
                                       [  5  ]
   A  (theta,phi) = A  (theta,phi) = - |-----| * eta*sin(theta)*sin(2*phi)]
    yz               zy                [24*PI]

                                            [1]    [                      ]
    A  (PAS) = 0                A  (EA) = i |-|  * | A   (EA) + A    (EA) |
     yz                          yz         [2]    [  2,1        2,-1     ]   */


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

MSVCDLL double  alpha()       const;
MSVCDLL double  beta()        const;
MSVCDLL double  gamma()       const;
MSVCDLL double  phi()         const;
MSVCDLL double  theta()       const;
MSVCDLL double  chi()         const;
MSVCDLL EAngles orientation() const;

MSVCDLL void alpha(double A);
MSVCDLL void beta(double  B);
MSVCDLL void gamma(double G);
MSVCDLL void phi(double   P);
MSVCDLL void theta(double T);
MSVCDLL void chi(double   C);
MSVCDLL void orientation(const EAngles& EA);
MSVCDLL void orientation(double A, double B, double G, bool deg=false);

MSVCDLL void rotate(const EAngles& EA);

// ____________________________________________________________________________
// E                RANK 2 SPATIAL TENSOR FROM CARTESIAN ARRAY
// ____________________________________________________________________________

/* These functions are used to set up an irreducible rank 2 spatial tensor
   from an input 3x3 Cartesian tensor. The input array should be real symmetric
   and traceless. The functions are designed to first directly obtain through
   diagonalization values of eta and delzz, then perform a Newton minimzation
   to obtain a set of Euler angles by that will rotate the PAS into the axes of
   the input array. Note that only the asymmetry is stored in the class object,
   both the delzz value and Euler Angles are returned to whatever requested
   generation of the spatial tensor from the 3x3 array.                      */

MSVCDLL friend void GetEtaDelzz(const matrix& A, 
            double& Aeta, double& Adelz, double ecut, double dcut);

MSVCDLL friend bool GetEAngles(const matrix& A, const EAngles& EA);

// ____________________________________________________________________________
// F                RANK 2 INTERACTION AUXILIARY FUNCTIONS
// ____________________________________________________________________________

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
   CartMx        Returns a 3x3 matrix representing the Cartesian tensor      */

MSVCDLL static coord  AisoDelzEta(const coord& AxAyAz);
MSVCDLL static void   SortAxAyAz(double& Ax, double& Ay, double& Az);
MSVCDLL        bool   CheckEta(double eta, bool warn=true) const;
MSVCDLL        bool   PAS()                                const;
MSVCDLL        bool   Symmetric()                          const;
MSVCDLL        matrix CartMx(double scale=1.0)             const;

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

MSVCDLL virtual      operator    ParameterSet( ) const;
MSVCDLL friend  void operator+= (ParameterSet& pset, const IntRank2A &IR2A);
MSVCDLL virtual void PSetAdd(ParameterSet& pset, 
                                       int iI=-1, int iS=-1, int warn=2) const;

// ----------------------------------------------------------------------------
// Functions To Output Irreducible Spatial Tensor To ASCII From A Parameter Set
// ----------------------------------------------------------------------------

        // Input                IR2A    : Rank 2 spatial tensor (this)
        //                      filename: Output file name
        //                      ofstr   : Output file stream
        //                      iI      : Interaction, 1st spin index (def. -1)
        //                      iS      :              2nd spin index (def. -1)
        //                      warn    : Warning level
        // Output               none    : Spatial tensor is written as a
        //                                parameter set to file filename
        //                                or output file stream ofstr

MSVCDLL bool write(const std::string& filename,
                                       int iI=-1, int iS=-1, int warn=2) const;
MSVCDLL bool write(std::ofstream& ofstr,
                                       int iI=-1, int iS=-1, int warn=2) const;

// ----------------------------------------------------------------------------
//  Functions To Read A Rank 2 Irreducible Spatial Tensor From An ASCII File
// ----------------------------------------------------------------------------

/* These functions allow users to specify rank 2 interaciton spatial tensor
   from parameters in an external ASCII file or from parameters in a GAMMA
   parameter set. These two methods are nearly equivalent since any ASCII file
   will be immediately read into a parameter set and the parameter set 
   subsequently used to set the interaction.

	Input		    IR2A    : Rank 2 spatial tensor (this)
                            filename: Output file name
                            pset    : Parameter set
                            idx     : Interaction index (default -1->none)
                            warn    : Warning output label
                                       0 = no warnings
                                       1 = warnings
                                      >1 = fatal warnings
	Output              TF	    : True if rank 2 spatial tensor is read
				      in from parameters in file filename
                                      or those in the parameter set
	Note                        : Classes derived from IntRank2(A) will
	   			      usually have their own read function   */

MSVCDLL virtual bool read(const std::string &filename,
                                         int idxI=-1, int idxS=-1, int warn=2);
MSVCDLL virtual bool read(ParameterSet& pset,
                                         int idxI=-1, int idxS=-1, int warn=2);
 
// ____________________________________________________________________________
// I                          STANDARD I/O FUNCTIONS
// ____________________________________________________________________________

MSVCDLL std::string ask_read(int argc, char* argv[], int& argq, int na=0);
 
/*
virtual void ask(int argc, char* argv[], int& argq,
                                     double& eta, double& theta, double& phi);
 
        // Input                IR2A	: Rank 2 spatial tensor (this)
        //                      argc    : Number of arguments 
        //                      argv    : Array of arguments 
        //                      argq	: IR2Auery value 
        //                      eta	: Rank 2. asymmetry 
        //                      theta	: Rank 2. orientation angle
        //                      phi	: Rank 2. orientation angle
        // Output               none    : The values of eta, theta, and
        //                                phi are set herein 
        // Note                         : This is INTERACTIVE! 

 
virtual void askset(int argc, char* argv[], int& argq);
 
        // Input                IR2A	: Rank 2 spatial tensor (this)
        //                      argc    : Number of arguments
        //                      argv    : Array of arguments
        //                      argq	: IR2Auery value
        // Output               none    : IR2A is set interactively
        // Note                         : This is INTERACTIVE!
 

virtual void askset();

        // Input                IR2A	: Rank 2 spatial tensor (this)
        // Output               none    : IR2A is set interactively
        // Note                         : This is INTERACTIVE!
*/

// ____________________________________________________________________________
// J                             OUTPUT FUNCTIONS
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
// Functions That Generate Single Strings to Simplify and Modularize Printing
//-----------------------------------------------------------------------------

/*         Input                IR2A    : Rank 2 spatial tensor (this)
                                CSF     : Float print format (e.g. "%6.3f")

          AsymmetryString          Asymmetry:                 x.xx
          ThetaString              Down From PAS z-Axis:    xxx.xx Degrees
          PhiString                Over From PAS x-Axis:    xxx.xx Degrees   */

MSVCDLL std::string AsymmetryString()          const;
MSVCDLL std::string PASString()                const;
MSVCDLL std::string AngleString(double angdeg) const;
MSVCDLL std::string ThetaString()              const;
MSVCDLL std::string ThetaString(double thedeg) const;
MSVCDLL std::string PhiString()                const;
MSVCDLL std::string PhiString(double   phideg) const;
MSVCDLL std::string AlphaString()              const;
MSVCDLL std::string AlphaString(double alpha)  const;
MSVCDLL std::string BetaString()               const;
MSVCDLL std::string BetaString(double  beta)   const;
MSVCDLL std::string GammaString()              const;
MSVCDLL std::string GammaString(double gamma)  const;
MSVCDLL std::string AuvString(double Aux, double Auy, double Auz,
                                                 const std::string& CSF) const;

//-----------------------------------------------------------------------------
// Functions That Generate Cartesian Strings To Simplify & Modularize Printing
//-----------------------------------------------------------------------------

/*         Input                IR2A    : Rank 2 spatial tensor (this)
                                CSF     : Float print format (e.g. "%6.3f")
				EA      : Euler angles for orientation
				A       : String for element label (e.g. A)

         Function CAString                 Function CartAStrings

          [A  , A  , A  ]           [A  , A  , A  ]
          [ xx   xy   xz]           [ xx   xy   xz]   [ x.x, x.x, x.x]
          [A  , A  , A  ]           [A  , A  , A  ] = [ x.x, x.x, x.x]
          [ yx   yy   yz]           [ yx   yy   yz]   [ x.x, x.x, x.x]
          [A  , A  , A  ]           [A  , A  , A  ]
          [ zx   zy   zz]           [ zx   zy   zz]                          */


MSVCDLL std::vector<std::string> CAStrings(const std::string& A="A")  const;
MSVCDLL std::vector<std::string> CartAStrings(const 
                                               std::string& CSF="%6.3f") const;
MSVCDLL std::vector<std::string> CartAStrings(
                  const EAngles& EA,     const std::string& CSF="%6.3f") const;
MSVCDLL std::vector<std::string> CartAStrings(
                  const IR2ACart& CCmps, const std::string& CSF="%6.3f") const;


//-----------------------------------------------------------------------------
// Functions That Generate Information Strings To Simplify & Modularize Printing
//-----------------------------------------------------------------------------

/* These functions return spatial tensor information in string format. This is
   done to facilitate printing, in particular printing of spatial tensors from
   within rank 2 interactions.  In this case we make a list of information
   thats usually displayed to the left of a 3x3 Cartesian matrix rep. of A.  */

MSVCDLL std::vector<std::string> InfoAStrings() const;
MSVCDLL std::vector<std::string> InfoAStrings(const EAngles& EA) const;

//-----------------------------------------------------------------------------
// Functions That Generate Spherical Strings to Simplify & Modularize Printing
//-----------------------------------------------------------------------------

/*         Input                IR2A    : Rank 2 spatial tensor (this)
                                CSF     : Float print format (e.g. "%6.3f")
				EA      : Euler angles for orientation
				A       : String for element label (e.g. A)

                                A    = xx.xxx
                                 2,0
                                A    = xx.xxx
                                 2,1
                                A    = xx.xxx
                                 2,-1
                                A    = xx.xxx
                                 2,2 
                                A    = xx.xxx
                                 2,-2                                        */

//virtual std::vector<std::string> SphAStrings()                const;
MSVCDLL         std::vector<std::string> SphA2Strings()               const;

//-----------------------------------------------------------------------------
//          Functions To Print The Tensor In Cartesian Format
//-----------------------------------------------------------------------------

        // Input                IR2A    : Rank 2 spatial tensor (this)
        //                      ostr    : Output stream
        //                      T       : Orientation angle theta (degrees)
        //                      P       : Orientation angle phi   (degrees)
        //                      F       : Title print flag (default 1 = print)
	//                      CSF     : Element output format
        // Output               none    : Rank 2 spatial tensor parameters
        //                                sent to the output stream

/*          Output Some Printed Lines Which Will Look Like The Following


                       [A  , A  , A  ]
      Output From      [ xx   xy   xz]   [ x.x, x.x, x.x]
     Third Overload    [A  , A  , A  ] = [ x.x, x.x, x.x]
                       [ yx   yy   yz]   [ x.x, x.x, x.x]
                       [A  , A  , A  ]
                       [ zx   zy   zz]                                 */

MSVCDLL virtual std::ostream& printCartesian(std::ostream& ostr,
                                    const std::string& CSF="%6.3f", int tpf=2);
MSVCDLL virtual std::ostream& printCartesian(std::ostream& ostr, 
                 const EAngles& EA, const std::string& CSF="%6.3f", int tpf=2);
MSVCDLL virtual std::ostream& printCartesian(std::ostream& ostr, 
                               const std::vector<std::string>& CAS, int tpf=2);

//-----------------------------------------------------------------------------
//            Here Are The Typically Utilized Output Functions
//-----------------------------------------------------------------------------

        // Input                IR2A	: Rank 2 spatial tensor (this)
        //                      ostr    : Output stream
        //                      fflag   : Format flag
        //                                  0 - Basic Parameters
        //                                 !0 - Full output
	//			tpf	: Title print flag (1=print)
        // Output               none    : Rank 2 spatial tensor parameters
        //                                placed into the output stream
        // Input                out	: Output stream;
        //                      IR2A	: Rank 2 tensor to write
        // Output			: Modifies output stream

MSVCDLL std::ostream& print(std::ostream& ostr, 
           std::vector<std::string>& CAS, std::vector<std::string>& IAS) const;
MSVCDLL virtual std::ostream& print(std::ostream& out,
                                                      int fflag=-1, int tpf=1);
MSVCDLL friend std::ostream& operator<< (std::ostream& out, IntRank2A& IR2A);


//-----------------------------------------------------------------------------
//     Other Functions That Generate Ouput Of The Spherical Rank 2 Tensor
//-----------------------------------------------------------------------------

        // Input                IR2A     : Rank 2 spatial tensor (this)
        //                      ostr    : Output stream
        //                      hdr     : Header for output
        //                      SS      : Array of strings to be printed
        //                                to the left of the A array(s)
        //                      phi     : Orientation angle (degrees)
        //                      theta   : Orientation angle (degrees)
	//			tpf	: Title print flag (1=print)
        // Output               none    : Rank 2 spatial tensor parameters
        // Note                         : For pretty output all the strings in
        //                                SAS should have the same length

MSVCDLL virtual std::ostream& print(std::ostream& ostr, 
             const std::string& hdr, const std::vector<std::string>& SS) const;
MSVCDLL virtual std::ostream& printSpherical(std::ostream& ostr, int tpf=1);
         
// ____________________________________________________________________________
// K         SPATIAL TENSOR SPHERICAL COMPONENTS FOR POWDER AVERAGES
// ____________________________________________________________________________
 
/* These functions are used to facilitate powder averaging in some cases. They
   return vectors of pre-calcualted values for the tensor oriented at evenly
   space angles.

           Input                IR2A	: Rank 2 spatial tensor (this)
                                Ntheta  : Number of orientation increments
           Return               A2m	: Partial spherical tensor component
	  				  for orientation theta from the 
                                          tensor PAS over a range of theta's
	   Note				: Asp2m values need to be added to some
					  scaled Bsp2m values to generate the
				      	  full A2m(theta,ph) spatial tensor
					  component values at any orientation
	   Note				: Theta spans [0, 180]               */
   
// ----------------- Parts Which Are Not ETA & Phi Dependent ------------------

MSVCDLL row_vector Asp20(int Ntheta)  const;
MSVCDLL row_vector Asp21(int Ntheta)  const;
MSVCDLL row_vector Asp22(int Ntheta ) const;
MSVCDLL matrix     Asp2s(int Ntheta)  const;
 
/*                       1/2
                  [  5  ]          2 		   Values MUST be added to
 A  (theta,phi) = |-----| * [ 3*cos (theta) - 1 ]  Bsp20*sin(theta)*sin(theta)
  20              [16*PI]	`		   to produce A20(theta,phi)

 
                      1/2
               [  5  ]              			Values MUST be added
 A   (theta) = |-----| * 3.0 * sin(theta) * cos(theta)  to values from Bsp21 
  2,1          [24*PI]              			to get A21(the,phi)
 							+ Re[Bpow21)*cos(the)]
							+ iIm(Bpow21)

                        1/2
                3 [  3 ]      2
 A   A(theta) = - |----|   sin (theta) 	    Values MUST be added to Bsp22
  2,2           2 [4*PI]   		    to produce full A22(theta,phi)
	  				    + Re(Bsp22)*(1+cos(the)cos(the))
	  				    + iIm(Bsp22B)*cos(the)       */

// ____________________________________________________________________________
// L               SPATIAL TENSOR COMPARISON FUNCTIONS
// ____________________________________________________________________________
 
/* These functions allow users to check if two spatial tensor are equivalent.
   Since eta is the only value that defines our irreducible rank 2 normalized
   spatial tensors, they will be equivalent if the asymmetry values match.

           Input                IR2A	: Rank 2 spatial tensor (this)
           			IR2A2	: Another rank 2 spatial tensor
           Output               T/F	: TRUE if IR2A2 and this are equal
          /F_list ==			- Equality
          /F_list !=			- Inequality                          */


MSVCDLL virtual bool operator==(const IntRank2A &IR2A2) const;
MSVCDLL virtual bool operator!=(const IntRank2A &IR2A2) const;

  };

#endif						// IntRank2A.h
