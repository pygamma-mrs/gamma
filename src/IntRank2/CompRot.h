/* CompRot.h ****************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      Composite Rotation				Interface	**
**                                                                      **
**      Scott Smith                                                     **
**      Copyright (c) 2001                                              **
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
** A variable of type rotation maintains a list of successive rotations **
** in terms of both Euler angles and Quaternions. The class contains 	**
** two vectors of equal length, one containing a list of Euler angles   **
** and one containing a list of corresponding Quaternions. Access	**
** functions allow users to get or set any of the rotations in the	**
** vector, and to get summed rotations over any of the rotations.	**
**                                                                      **
*************************************************************************/

#ifndef   CompRot_h_			// Is file already included?
#  define CompRot_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface 			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Basics/ParamSet.h>		// Include parameter sets
#include <Level1/coord.h>		// Include coordinates
#include <Level2/EAngles.h>		// Include Euler angles
#include <Level2/Quaternion.h>		// Include Quaternions
#include <string>			// Include stdlibc++ strings

class CompRot
  {
  std::vector<EAngles> EAs;		// Vector of Euler angles
  std::vector<quatern> Qs;		// Vector of Quaternions
  EAngles              sumEA;		// Summed rotation Euler angles
  quatern              sumQ;		// Summed rotation Quaternion

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                  COMPOSITE ROTATION ERROR HANDLING
// ____________________________________________________________________________

/*       Input 		      ROT     : Composite Rotation (this)
                              eidx    : Flag for error type
                              noret   : Flag for return (0=return)
                              pname   : String included in message
         Output               none    : Error message
                                        Program execution stopped if fatal   */

         void ROTerror(int eidx,                           int noret=0) const;
         void ROTerror(int eidx, const std::string& pname, int noret=0) const;
volatile void ROTfatal(int eidx)                                        const;
volatile void ROTfatal(int eidx, const std::string& pname)              const;
 
// ____________________________________________________________________________
// ii                COMPOSITE ROTATION CHECKING FUNCTIONS
// ____________________________________________________________________________

bool ChkIdx(int   i,         int warn=2) const;// Check index i OK
bool ChkRange(int i, int nr, int warn=2) const;// Check range [i,i+nr)

// ____________________________________________________________________________
// iii              COMPOSITE ROTATION AUXILIARY FUNCTIONS
// ____________________________________________________________________________

void SetSum();				// Calculate summed rotation

// ____________________________________________________________________________
// iv          COMPOSITE ROTATION PARAMETER SET SETUP FUNCTIONS
// ____________________________________________________________________________

/* These functions are used to parse a parameter set and fill up a rotation.

           Input        ROT     : Composite Rotation (this)
                        pset    : A parameter set
                        pfx     : Prefix on rotation parameters
                        warn    : Warning level
                                    0 - no warnings
                                    1 - warnings
                                   >1 - fatal warnings
       Output           TF      : Composite Rotation is set
                                  from parameters in pset		     */

bool SetCmpRot(const      ParameterSet& pset,int pfx=-1,  int warn=2);
bool SetRotation(const    ParameterSet& pset,             int idx=-1);
bool GetEulerAngles(const ParameterSet& pset,EAngles& EA, int idx=-1);
bool GetQuaternion(const  ParameterSet& pset,quatern&  Q, int idx=-1);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A                 COMPOSITE ROTATION CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/* Rotation construction. Since each rotation may contain a series of single
   rotations we allow use arrays for their construction. Thus the rotation may
   be constructed empty (no rotations), with a single rotation (either using
   Euler angles or a Quaternion), or with a vector of rotations (again
   either using Euler angles or Quaternions).                                */

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

MSVCDLC CompRot();

// ----------------------------------------------------------------------------
//                   Constructors Using Single Rotations
// ----------------------------------------------------------------------------

MSVCDLC CompRot(double alpha, double beta=0, double gamma=0);
MSVCDLC CompRot(const  EAngles& EA);
MSVCDLC CompRot(const  quatern& Q);
 
// ----------------------------------------------------------------------------
//                  Constructors Using Multiple Rotations
// ----------------------------------------------------------------------------

MSVCDLC CompRot(const std::vector<EAngles>& EAvec);
MSVCDLC CompRot(const std::vector<quatern>& Qvec);

// ----------------------------------------------------------------------------
//                    Construction Using Parameter Sets
// ----------------------------------------------------------------------------

MSVCDLC CompRot(const ParameterSet& pset, int idx=-1, int warn=2);

// ----------------------------------------------------------------------------
//                          Assignment and Destruction
// ----------------------------------------------------------------------------

MSVCDLL void operator= (const CompRot& ROT);
MSVCDLC      ~CompRot();

// ____________________________________________________________________________
// B                  INDIVIDUAL COMPOSITE ROTATION ACCESS FUNCTIONS
// ____________________________________________________________________________

MSVCDLL EAngles EA(int    i) const;
MSVCDLL quatern Q(int     i) const;
MSVCDLL double  alpha(int i) const;
MSVCDLL double  beta(int  i) const;
MSVCDLL double  gamma(int i) const;

MSVCDLL void EA(const  EAngles& ea, int i);
MSVCDLL void Q(const   quatern&  q, int i);
MSVCDLL void alpha(double        A, int i);
MSVCDLL void beta(double         B, int i);
MSVCDLL void gamma(double        G, int i);

// ____________________________________________________________________________
// C                    SUMMED COMPOSITE ROTATION ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These allow users to access composite rotations, either over the entire
   vector or any sub-vector. Since composite rotations may depend upon 
   multiple single rotations they may not be set, only obtained.             */

MSVCDLL EAngles EA()    const;
MSVCDLL quatern Q()     const;
MSVCDLL double  alpha() const;
MSVCDLL double  beta()  const;
MSVCDLL double  gamma() const;

MSVCDLL EAngles EA(int    i, int j);
MSVCDLL coord   ABG(int   i, int j);
MSVCDLL quatern Q(int     i, int j);
MSVCDLL double  alpha(int i, int j);
MSVCDLL double  beta(int  i, int j);
MSVCDLL double  gamma(int i, int j);

// ____________________________________________________________________________
// D                  COMPOSITE ROTATION PARSING FUNCTIONS
// ____________________________________________________________________________

/* These function allow users to obtain any part of the composite rotation   */

MSVCDLL CompRot operator() (int i, int nr) const;
MSVCDLL void    operator+= (const EAngles& EA);
MSVCDLL void    operator+= (const quatern& Q);

// ____________________________________________________________________________
// E                         PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

/* This class has no implicit knowledge of how an individual rotation has been
   expressed. Threre are three valid ways to express an individual rotation:

                   1.) As Euler Angles - EAngles: (alpha,beta,gamma) EA(#)
                   2.) As a Quaternion - quatern: (A,B,C,D)          QRT(#)
                   3.) As three angles - double:  alpha              alpha(#)
                                         double:  beta               beta(#)
                                         double   gamma              gamma(#)

   However, we must choose one for output of/storage in a parameter set. We
   choose the simplest to store, method 1.) above.

	   Input		AROT	: Composite Rotation (this)
                                pfx     : Rotation prefix (def -1)
	  			filename: Output file name
           Output               pset    : Parameter set with only
                                          composite rotation parameters
				void 	: Composite Rotation is written as
					  a parameter set to file filename

     Function                                 Purpose
   ------------         -------------------------------------------------------
   ParameterSet         Convert composite rotation into a parameter set
   operator +=          Adds composite rotation to existing parameter set
   PSetAdd              Adds composite rotation to existing parameter set, this
                        allows for a prefix to be used in the parameters.
   write		Writes rotation to ASCII file as parameter set       */

MSVCDLL             operator ParameterSet( ) const;
MSVCDLL friend void operator+= (ParameterSet& pset, const CompRot &ROT);
MSVCDLL        void PSetAdd(ParameterSet& pset,   int pfx=-1) const;
MSVCDLL        void write(const std::string &filename, int pfx=-1) const;

// ____________________________________________________________________________
// H                      COMPOSITE ROTATION INPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions are used to specify the rotation from either a working
   GAMMA parameter set or an external ASCII file (in parameter set format).
   The functions do NOT allow for the suffix (#), since it is used/reserved
   specific rotation indices. The functions do allow for the prefix [#] so that
   multiple composite rotations may be defined in the same file by switching
   the prefix indices.

           Input                ROT     : Composite Rotation
                                file    : Input file name
                                pset    : Parameter set
                                argv    : Vecotr of argc arguments
                                argn    : Argument index
                                idx     : Composite Rotation index (default -1->none)
           Output               none    : Composite Rotation is read in
                                          from parameters in file filename
                                file    : Name of file used to set rotation  */

MSVCDLL bool   read(const std::string& file,    int idx=-1, int warn=2);
MSVCDLL bool   read(const ParameterSet& pset,   int idx=-1, int warn=2);
MSVCDLL std::string ask_read(int argc, char* argv[], int argn,   int idx=-1);
 
// ____________________________________________________________________________
// J                       COMPOSITE ROTATION OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions will output information concerning the rotation to any
   specified output stream.

           Input                ROT	: Composite Rotation (this)
                                ostr    : Output stream
                                fflag   : Format flag
                                            0 - Sum rotation only
                                           !0 - Full composite rotation
           Output               none    : Composite Rotation information
                                          placed into the output stream      */

MSVCDLL        std::ostream& print(std::ostream& out, int fflag=-1) const;
MSVCDLL friend std::ostream& operator<<  (std::ostream& out, const CompRot& ROT);

  };

#endif							// CompRot.h
