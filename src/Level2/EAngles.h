/* EAngles.h *************************************************-*-c++-*-
**									**
** 	                         G A M M A				**
**									**
**	Euler Angles		                    Interface     	**
**						 			**
**	Copyright (c) 2001					 	**
**	Scott Smith				 			**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                    		  	**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**  									**
** This file contains a set of three Euler angles {alpha, beta, gamma}  **
** that are used in general rotations.  The class structure exists in   **
** order to maintain the angles within specific ranges as well as to    **
** facilite I/O and composite rotations.				**
**  									**
** The euler angle ranges maintained in this class are:			** 
**  									**
**          alpha = gamma = [0,360]       beta = [0,180]		**
**  									**
** This causes no loss in generality because rotations may still cover 	**
** all of three dimensional space.					**
**  									**
*************************************************************************/

///Chapter Class EAngles 
///Section Overview
///Body    The class EAngles is a data type for a three Euler angles.
///Section Available Euler Angle Functions

#ifndef   EAngles_h_				// Is file already included?
#  define EAngles_h_ 1				// If no, then remebmer it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// this is the interface
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <Basics/ParamSet.h>			// Include parameter sets
#include <Level1/coord.h>			// Include 3D coordinates
#include <string>				// Include libstdc++ strings


class EAngles
  {
         double _alpha;				// Euler Angle alpha (radians)
         double _beta;				// Euler Angle beta  (radians)
         double _gamma;				// Euler Angle gamma (radians)
  static double AngCutoff;			// Zero cutoff for angle

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE Functions ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                   Class Euler Angles Error Handling
// ____________________________________________________________________________

/*      Input                   EA      : Euler Angles (this)
                                eidx    : Error index
                                nrt     : Flag for linefeed (0=linefeed)
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

         void EAerror(int eidx,                        int n=0) const;
         void EAerror(int eidx, const std::string& pn, int n=0) const;
volatile void EAfatal(int eidx=0)                               const;

// ____________________________________________________________________________
// ii            Class Euler Angles Private Facilitator Functions 
// ____________________________________________________________________________

/* These are chekcing functions for Euler angles. They insure that the angles
   always reside within the class specified ranges:

               alpha = gamma = [0,360]       beta = [0,180]		     */

void SetAngles(double alpha, double beta, double gamma, bool d=false);
void SetAlpha(double  alpha, bool deg=false);	// Input in radians/deg
void SetBeta(double   beta,  bool deg=false);	// Input in radians/deg
void SetGamma(double  gamma, bool deg=false);	// Input in radians/deg

// ____________________________________________________________________________
// iii          Class Euler Angles Private Parameter Set Functions
// ____________________________________________________________________________

/* These functions try and set up Euler angles from parameters found in
   a particular parameter set.  Each object of type EAngles can be specified
   in one of two ways:

             1.) As Euler Angles - coord:   (alpha,beta,gamma) EAngles(#)
             2.) As three angles - double:  alpha              EAalpha(#)
                                   double:  beta               EAbeta(#)
                                   double   gamma              EAgamma(#)

   where the # is used to indicate the Euler angle index. It is assumed that
   the angles are specified in degrees, although the names may have "Rad"
   appended before "(#)" if the values are specified in radians. Again, we
   strictly enforce the angle ranges of alpha = gamma = [0,360]  and
   beta = [0,180].                                                           */

bool SetEAngles(const ParameterSet& pset, int idx=-1, bool warn=true);
bool SetEASet(const   ParameterSet& pset, int idx=-1, bool warn=true);
bool Set3Angles(const ParameterSet& pset, int idx=-1, bool warn=true);

// ----------------------------------------------------------------------------
// ---------------------------- Public Functions ------------------------------
// ----------------------------------------------------------------------------

 public:

// ____________________________________________________________________________
// A               Class Euler Angles Constructors/Destructor
// ____________________________________________________________________________

/* These constructors set up a new set of Euler angles. All need specification
   of the three components { alpha, beta, gamma }, either explicitly or by
   default. Of course, some of the angles may be zero.  Keep in mind that we
   restrict the angle ranges to alpha = gamma = [0,360] & beta = [0,180].

        Input Arguments                     Result
       ----------------       ----------------------------------------------
               -              Zero Euler Angles {alpha=0, beta=0, gamma=0}
              EA              New Euler Angles indentical to those in EA
       alpha,beta,gamma       Euler Angles with these 3 angles
                              (Angle ranges are strictly maintained)        */

MSVCDLC EAngles();
MSVCDLC EAngles(double alpha, double beta=0, double gamma=0, bool deg=false);
MSVCDLC EAngles(const coord&   EA, bool deg=true);
MSVCDLC EAngles(const EAngles& EA);
MSVCDLC ~EAngles();

MSVCDLL EAngles& operator= (const EAngles& EA);
//EAngles& operator= (const coord&   ABG);

// ____________________________________________________________________________
// B                     Euler Angles ACCESS Functions
// ____________________________________________________________________________

///Center Access Functions
///F_list alpha			- Access to Euler angle alpha
///F_list beta			- Access to Euler angle beta
///F_list gamma			- Access to Euler angle gamma

MSVCDLL double alpha() const;			// Get alpha (radians)
MSVCDLL double beta( ) const;			// Get beta  (radians)
MSVCDLL double gamma() const;			// Get gamma (radians)

MSVCDLL void alpha(double A);			// Set alpha (radians)
MSVCDLL void beta(double  B);			// Set beta  (radians)
MSVCDLL void gamma(double G);			// Set gamma (radians)

// ____________________________________________________________________________
// C            Class Euler Angles Composite Rotation Functions
// ____________________________________________________________________________

/* These functions handle composite rotations, i.e. generation of a set of
   Euler angles that represents two successive Euler angle rotations.
   Generation of summed rotation angles is done through use of Quaterions.   */

MSVCDLL EAngles operator*  (const EAngles& EA) const;
MSVCDLL EAngles&    operator*= (const EAngles& EA);
MSVCDLL EAngles&    operator&= (const EAngles& EA);
MSVCDLL EAngles composite  (const EAngles& EA) const;

// ____________________________________________________________________________
// D             Class Euler Angles Parameter & Parameter Set Functions
// ____________________________________________________________________________

//-----------------------------------------------------------------------------
//                       Single Parameter Functions
//-----------------------------------------------------------------------------

        // Input               EA    : A set of Euler angles (this)
        //                     pname : A parameter name
        // Return              par   : A GAMMA parameter of type coordinate
        //                             with the name pname

MSVCDLL SinglePar param(const std::string& pn,                      bool deg=true) const;
MSVCDLL SinglePar param(const std::string& pn,const std::string& ps,bool deg=true) const;

// ----------------------------------------------------------------------------
/* These will 1.) construct a parameter set containing an Euler angle set.
              2.) add an Euler angle set to an existing parameter set.
              3.) write an Euler angle set to an ASCII file in pset format.  */

MSVCDLL               operator ParameterSet( ) const;
MSVCDLL friend void   operator+= (ParameterSet& pset, const EAngles& EA);
MSVCDLL void PSetAdd(ParameterSet& pset,        int idx=-1, bool deg=true) const;
MSVCDLL void write(const std::string &filename, int idx=-1, bool deg=true) const;

// ----------------------------------------------------------------------------
//           Functions To Make Euler Angles From A Parameter Set
// ----------------------------------------------------------------------------

/* These will 1.) set Euler angles from values in a parameter set.
              2.) set Euler angles from parameters in an ASCII file.         */

MSVCDLL bool read(const std::string &filename, int idx=-1, int warn=2);
MSVCDLL bool read(const ParameterSet& pset,    int idx=-1, int warn=2);

// ____________________________________________________________________________
// E                    Class Euler Angles I/O Functions
// ____________________________________________________________________________

///Center I/O Functions
///F_list print               - Print to output filestream
///F_list <<                  - Standard output


MSVCDLL std::ostream& print(std::ostream& ostr, bool deg=true, bool hdr=true) const; 
MSVCDLL friend std::ostream& operator <<    (std::ostream& ostr, const EAngles& EA);

	// Input		EA    : Euler angles 
        //                      ostr  : Output stream
        // Output               none  : Euler angle set is sent
	//				to the output stream

// ____________________________________________________________________________
// F              Class Euler Angle Container Support Functions
// ____________________________________________________________________________

/* These allow the user to make STL lists and vectors from Euler Angles. If
   not present some compilers complain about code such as list<EAngles>      */

MSVCDLL static void SetCutoff(double co=-1);
MSVCDLL        bool operator== (const EAngles& EA) const;
MSVCDLL        bool operator!= (const EAngles& EA) const;
MSVCDLL        bool operator<  (const EAngles& EA) const;
MSVCDLL        bool operator>  (const EAngles& EA) const;

// ____________________________________________________________________________
// F                Class Euler Angle Auxiliary Functions
// ____________________________________________________________________________

/*       Function                            Purpose
         ========   ===========================================================
          equal     Compares two sets of Euler angles to see if they produce
                    the same rotation. Equality can occur even if the Euler
                    angles are not the same, so long as the rotations defined
                    are equivalent (produce the same rotation matrix)
	            Note that operator == is NOT the same as equal. The former
                    does an exact comparison of angles, not rotations.
          inverse   Returns a set of Euler angles for the inverse rotation.
            RMx     Returns a 3x3 rotation matrix that will take one coordinate      
                    system into another as defined by the Euler angles. The
                    flag inv allows for generation of the inverse rotation
                    matrix.                                                  */

MSVCDLL bool    equal(const EAngles& EA, double CUTOFF=1.e-10) const;
MSVCDLL EAngles inverse()                                      const;
MSVCDLL matrix  RMx(bool inv=false)                            const;

MSVCDLL matrix Rmx()    const;				// DEPRECATED
MSVCDLL matrix invRmx() const;				// DEPRECATED

  };

/*****************************************************************************/
/*****************************************************************************/
/*                       Class Euler Angle Constants                         */
/*****************************************************************************/
/*****************************************************************************/

extern const EAngles EAzero;		                // Zero Euler angles

#endif								// EAngles.h

