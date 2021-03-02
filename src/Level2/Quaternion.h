/* Quaternion.h *************************************************-*-c++-*-
**									**
** 	                         G A M M A				**
**									**
**	Quaternion		                    Interface     	**
**						 			**
**	Copyright (c) 1991, 1999				 	**
**	Scott Smith				 			**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fur physikalische Chemie		 			**
**	8092 Zurich / Switzerland		 			**
**                                                    		  	**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**  									**
** This file contains class quatern, a data type which facilitates the	**
** use of quaternions.							**
**  									**
*************************************************************************/

///Chapter Class Quatern
///Section Overview
///Body    The class Quatern is a data type for a quaternion.
///Section Available Quaternion Functions

#ifndef   quatern_h_				// Is file already included?
#  define quatern_h_ 1				// If no, then remebmer it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// this is the interface
#  endif

# include <GamGen.h>				// Know MSVCDLL (__declspec)
# include <Level1/coord.h>			// Include 3D coordinates
# include <Level2/EAngles.h>			// Inlcude Euler angles
# include <Matrix/matrix.h>			// Include GAMMA matrices

class quatern
  {
  double        AQ, BQ, CQ, DQ;			// Four quaternion elements
  static int    QRange;				// Euler angle default range 
  static bool   SinPos, CosPos, TanPos, SCTF;	// Flags asin/acos/atan range
  static double ElemCutoff;                     // Zero cutoff for element
  typedef quatern quartern;			// Old name still OK (deprec.)

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                    CLASS QUATERNION ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   Qrt     : A quaternion (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

         void Qerror(int eidx,                          int noret=0) const;
         void Qerror(int eidx,   const std::string& pn, int noret=0) const;
volatile void Qfatal(int eidx=0)                                     const;
volatile void Qfatal(int eidx,   const std::string& pn)              const;

// ____________________________________________________________________________
// ii              PRIVATE FACILITATOR FUNCTIONS FOR QUATERNIONS
// ____________________________________________________________________________


/* These are testing functions for Quaternions.  There are two types of tests.
   The first deals with the quaternion norm which MUST always be one.  The 2nd
   deal with the 4x4 rotation matrix associated with the Quaternion. The
   array must not only be 4x4 but also orthonormal.                          */

bool CheckNorm(int warn=2) const;
bool CheckNorm(double A,double B,double C,double D,int warn=2) const;
bool CheckNorm(const matrix& Rotmx, bool warn=true) const; 

// ____________________________________________________________________________
// iii         PRIVATE PARAMETER PARSING FUNCTIONS FOR QUATERNIONS
// ____________________________________________________________________________

/* These allow for Quaternions to be set from parameters found in GAMMA
   parameter sets.  The Quaternion is specified by Quatern(#) where # is the
   index of the Quaternion. The parameter is assumed in string format and the
   data in the form (A, B, C, D) which will be parsed. Once parsed, the
   values will be checked for consistency.                                   */

bool SetQuatern(const ParameterSet& pset, int idx=-1, int warn=2);

bool SetSinPos() const;
bool SetCosPos() const;
bool SetTanPos() const;

double  FindBeta()                             const;
double  FindAlpha()                            const;
double  FindAlpha(double beta)                 const;
double  FindGamma()                            const;
double  FindGamma(double dbeta, double dalpha) const;
EAngles FindEAs()                              const;
double  GetAngle(double sinval, double cosval) const;

// ____________________________________________________________________________
//                QUATERNION COMMON INTERNAL FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

 public:

// ____________________________________________________________________________
// A                 CLASS QUATERNION CONSTRUCTORS/DETRUCTOR
// ____________________________________________________________________________

///Center Quaternion Algebraic
///F_list quatern		- Constructor
///F_list =			- Assignment

/* These constructors set up a new quaternion. All demand specification of the
   four quaternion components, either explicitly or by default. These may be
   set using three Euler angles as well. It is important to note that, unlike
   EAngles which always maintain angles in radians, a coordinate (ABG) is
   assumed to have the three angles specified in degrees.
 
        Input Arguments                     Result
        ---------------       ------------------------------------------
               -              Zero Rotation Quaternion {a=0, b=0, g=0}
           ABG,int	      Quaternion representing ABG or its inverse
            EA,int            Quaternion representing EA or its inverse
            Qrt,int           Quaternion equal to Qrt or its inverse
          QA,QB,QC,QD         Quaternion having these 4 components
                              (Insures A**2 + B**2 + C**2 + D**2 = 1)
 
   Note that the relationship between the Quaternion herein and the set of
   three Euler angles is in accordance with Spiess, JMR, 61, 356 (85). In       
   particular see Equation [3] of that article. Also see the article 
   Ernst, JMR, 63, 133 (85). The norm of the 4 components is always 1.      */

MSVCDLC quatern();
MSVCDLC quatern(const coord&   ABG, bool inv=false);
MSVCDLC quatern(const EAngles& EA,  bool inv=false);
MSVCDLC quatern(const quatern& Qrt, bool inv=false);
MSVCDLC quatern(const ParameterSet& pset, int idx=-1, int warn=2);
MSVCDLC quatern(double QA, double QB, double QC, double QD, bool inv=false);

// ----------------------------------------------------------------------------
//                  Assignment Operators & Destructor
// ----------------------------------------------------------------------------

// Remember, we assume coord is in degrees whereas EAngles is always in radians

MSVCDLC          ~quatern();
MSVCDLL quatern& operator= (const quatern& QRT);
MSVCDLL quatern& operator= (const coord&   EA);
MSVCDLL quatern& operator= (const EAngles& EA);

// ____________________________________________________________________________
// B                     QUATERNION ACCESS FUNCTIONS
// ____________________________________________________________________________

///Center Access Functions
///F_list A			- Access to Quaternion value A
///F_list B			- Access to Quaternion value B
///F_list C			- Access to Quaternion value C
///F_list D			- Access to Quaternion value D

MSVCDLL double A() const;
MSVCDLL double B() const;
MSVCDLL double C() const;
MSVCDLL double D() const;

// ____________________________________________________________________________
// C                 QUATERNION TO EULER ANGLE FUNCTIONS
//         Q == {A, B, C, D} ==> {alpha, beta, gamma} == EAngles
// ____________________________________________________________________________

///Center Euler Angle Functions
///F_list beta                  - Access to Euler angle beta        (radians)
///F_list alpha                 - Access to Euler angle alpha       (radians)
///F_list gamma                 - Access to Euler angle gamma       (radians)
///F_list EA                    - Access to Euler angles (EAngles in radians)
///F_list ABG			- Access to Euler angles (coord   in degrees)

/* These functions allow users to obtain any or all of the Euler angles 
   associated with a given Quaternion. The angle(s) is/are returned in
   units of radians.  These functions are implemented from the article

                     Spiess, JMR, 61, 356 (1985)
  
   The functions perform the inverse of eq. [7] in the above article.

   A second article is also useful:    Ernst, JMR, 63, 133 (85)

   The angle alpha can be determined more rapidly if angle beta is known
   and likewise angle gamma is generated more rapidly if both angles beta
   and alpha are known.  The corresponding functions allow the user to supply
   the beta (and alpha) values in order to facilitate the calculation.
   Note that there is sometimes an ambiguity concerning the returned angle
   that arises from use of the arc sine and arc cosine functions.  To resolve
   any conflicts, the angles may be generated from mulitple formulae.        */

MSVCDLL double  alpha() const;
MSVCDLL double  beta()  const;
MSVCDLL double  gamma() const;
MSVCDLL EAngles EA()    const;
MSVCDLL coord   ABG()   const;

// ____________________________________________________________________________
// D               CLASS QUATERNION COMPOSITE ROTATION FUNCTIONS
// ____________________________________________________________________________
     
///Center Composite Rotations
///F_list composite             - Composite rotations
///F_list RMx                   - Rotation matrix (4x4)

/* These functions allow users to deal with Quaternion rotations & composite
   rotations.  Given two quaternions, Qrt1 & Qrt2 (or the equivalent rotation
   in terms of Euler angles), one may obtain the quaternion for the composite
   rotation (i.e. Qrt1 followed by Qrt2).  The same may be accomplished by
   use of a Quaternion based rotation matrix (4x4).  These functions are based
   on an article by
                          Spiess, JMR, 61, 356 (1985)                         */

MSVCDLL quatern operator*  (const quatern& Q) const;
MSVCDLL quatern&    operator*= (const quatern& Q);
MSVCDLL quatern&    operator&= (const quatern& Q);
MSVCDLL friend  quatern  operator*  (const matrix&  Rmx,  const quatern& Q);

MSVCDLL quatern composite(const quatern& Q, bool rev=false) const;	// Deprecated
MSVCDLL matrix RotMx() const; 						// Deprecated
MSVCDLL matrix RMx() const;

MSVCDLL friend  quatern composite(const EAngles& EA1,  const EAngles& EA2);
MSVCDLL friend  quatern composite(const coord&   EA1,  const coord&   EA2);
MSVCDLL friend  quatern composite(const quatern& Qrt1, const quatern& Qrt2);	// Deprecated
MSVCDLL friend  quatern composite(const matrix&  Rmx,  const quatern& Qrt1);	// Deprecated

// ____________________________________________________________________________
// E                  CLASS QUATERNION AUXILIARY FUNCTIONS
// ____________________________________________________________________________

///Center General Functions
///F_list norm                  - Sum of squares
///F_list inverse		- Inverse quaternion
 
	// Input		Qrt	: Quaternion (this)
	// Output		Qnorm	: Computes a "norm" for the input
	//				  Quaternion by the formula
	//				      A**2+B**2+C**2+D**2
	// Note				: All valid quaternions should
	//				  have a norm of 1

MSVCDLL double  norm()    const;
MSVCDLL quatern inverse() const;

// ____________________________________________________________________________
// F                    CLASS QUATERNION I/O FUNCTIONS
// ____________________________________________________________________________

///Center I/O Functions
///F_list print               - Print to output filestream
///F_list <<                  - Standard output

MSVCDLL std::ostream& print(std::ostream& ostr, bool nf=true, bool hdr=true) const; 
MSVCDLL friend std::ostream& operator <<    (std::ostream& ostr, const quatern& Quar);

	// Input		Quar  : Quaternion
        //                      ostr  : Output stream
        // Output               none  : Quaternion is sent
	//				to the output stream

// ____________________________________________________________________________
// G              Class Quaternion Container Support Functions
// ____________________________________________________________________________

MSVCDLL bool operator== (const quatern& Quar) const;
MSVCDLL bool operator!= (const quatern& Quar) const;
MSVCDLL bool operator<  (const quatern& Quar) const;
MSVCDLL bool operator>  (const quatern& Quar) const;

// ____________________________________________________________________________
// H                  Class Quaternion Range Functions
// ____________________________________________________________________________
 
/* The Euler angles associated with a Quaternion can in principle have any 
   value. However, most users prefer to limit the range over which the three
   angles are defined.  These functions provide a means by which the Euler
   angles input to and returned from this class can be limited to reside
   within a specified range.  The range depends upon the flag "range". Users
   can associate values of range with implied Euler angle limits by adjusting
   the functions below.  Zero is reserved for no imposed limits.  As long as
   all 3D space can be sampled by the Euler angle range, this results in no
   loss in gerenality (multiple Euler angles sets produce the same rotation).
   They only affect of range will be in the Euler angle values output.  The
   quarternions and the associated rotation(s) will remain the same.        */

//       int  range() const; 
//static void range(int r);

// ____________________________________________________________________________
// I            Class Quaternion Parameter & Parameter Set Functions
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//                          Single Parameter Functions
//-----------------------------------------------------------------------------

        // Input                Qrt     : Quaternion (this)
        //                      pname   : A parameter name
        //                      pstate  : A parameter statement
        // Return               par     : A GAMMA parameter of type Quaternion
        //                                with default name and comment unless
        //                                specified

MSVCDLL SinglePar param()                                                    const;
MSVCDLL SinglePar param(const std::string& pname)                            const;
MSVCDLL SinglePar param(const std::string& pname, const std::string& pstate) const;

//-----------------------------------------------------------------------------
//                       Parameter Set From Quaternion
//-----------------------------------------------------------------------------

MSVCDLL operator ParameterSet( ) const;
MSVCDLL friend void operator+= (ParameterSet& pset, const quatern& Qrt);
MSVCDLL bool PSetAdd(ParameterSet& pset, int idx=-1, int pfx=-1) const;

//-----------------------------------------------------------------------------
//                       Parameter Set File From Quaternion
//-----------------------------------------------------------------------------

        // Input                Qrt     : Quaternion (this)
        //                      fo      : Output file name
        //                      of      : Output file stream
        //                      sfx     : Quaternion index  (default -1)
        //                      pfx     : Quaternion prefix (default -1)
        //                      warn    : Warning level
        // Output               none    : Quaterion is written as a parameter
        //                                to file or output file stream

MSVCDLL bool write(const std::string& fo,int idx=-1,int pfx=-1,int warn=2) const;
MSVCDLL bool write(    std::ofstream& of,int idx=-1,int pfx=-1,int warn=2) const;

//-----------------------------------------------------------------------------
//                         Quaternion From Parameter Set
//-----------------------------------------------------------------------------

        // Input                Qrt     : Quaternion (this)
        //                      filein  : An input (ASCII) file
        //                      pset    : A parameter set
        //                      indx    : Quaterion index
        //                      warn    : Warning level
        //                                  0 = no warnings
        //                                  1 = non-fatal warnings
        //                                  2 = fatal warnings
        // Output               none    : Quaterion filled with
        //                                values specified in filein
        //                                or in pset

MSVCDLL bool read(const std::string&  filein, int indx=-1, int warn=2);
MSVCDLL bool read(const ParameterSet& pset,   int indx=-1, int warn=2);

// ____________________________________________________________________________
// J             INVERSE TRIGONOMETRIC ANGLE DISCERNMENT FUNCTIONS
// ____________________________________________________________________________

/* These function help converting a Quaternion into Euler angles. Such
   conversions are more complicated than one would like because of the
   multi-valued nature of the (inverse) trigonometric functions. Because of
   this I've tried to simplyfy the procedure by 1st determininng how the
   computer sees the inverse functions asin & acos.

   The inverse sine and cosine functions in C/C++, asin and acos, will restrict
   their output values to be "principal values", always residing within a range
   over whcih sine and cosine are single valued. On some systems this is set
   to be the range of angles [-90, 90] whereas on some systems it is set to be
   the range [0, 180]. For example, on my Sun SPARC system  asin() returns in
   the range [-pi/2, pi/2] radians and acos() returns in the range [0,pi]
   radians.  The same is true for my SGI Irix systems. However, on my PC/Cygwin
   systems the ranges on both are [-pi/2, pi/2]. These functions just test a
   few values to discern what the standard is on the computer running this.  */

MSVCDLL static bool ASinPos();
MSVCDLL static bool ACosPos();
MSVCDLL static bool ATanPos();

// ----------------------------------------------------------------------------
//                 For Outputting Tests During Conversions
// ----------------------------------------------------------------------------

MSVCDLL void ShowConversion() const;

// ----------------------------------------------------------------------------
//               Check For Valid Quaternion Rotation Matrix
// ----------------------------------------------------------------------------

MSVCDLL static bool ValidRMx(const matrix& R, bool msgs=true);

  };

extern  quatern composite(const EAngles&,  const EAngles&);
extern  quatern composite(const coord&,  const coord&);
extern  quatern composite(const quatern&, const quatern&);	// Deprecated
extern  quatern composite(const matrix&,  const quatern&);	// Deprecated

//typedef quatern quartern;				// Either name OK
//int QRange = 1;				// a=g=[0,360],b=[0,pi]

#endif							// Quaternion.h

