/* CubicSys.h ***************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**	CubicSys 		           	Interface Definition	**
**							 	 	**
**	Copyright (c) 2000						**
**	Scott A. Smith				 			**
**  National High Magnetic Field Laboratory                           	**
**  1800 E. Paul Dirac Drive                                          	**
**  Box 4005                                                          	**
**  Tallahassee, FL 32306                                             	**
**						 			**
**  $Header: $
**						 			**
*************************************************************************/

/*************************************************************************
**									**
** Description								**
**							 		**
** Class CubicSys represents a magnetic ion in a crystal lattice with	**
** cubic symmetry. The system will be associated with a specific	**
** ionized element (e.g. Ce3+) exhibiting intrinsic spin (e.g. J=7/2).	**
** This system will be used in generating Hamiltonians and spectra in   **
** EPR computations. The class also manages a list of all possible ions	**
** intrinsically supported.						**
**									**
*************************************************************************/

#ifndef   CubicSys_h_				// Is file already included?
#  define CubicSys_h_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// If it is the GNU compiler
#    pragma interface				// then this is the interface
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <ESRLib/MagIon.h>			// Include Magnetic Ions
#include <fstream>				// Include file stream handling

class CubicSys
  {
  static MagIon* CubicIons;			// Array of known magnetic ions
  static int NIons;				// Number of known magnetic ions
  int Ion;					// A magnetic ion index (in list)


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                      CLASS CUBIC SYSTEM ERROR HANDLING
// ____________________________________________________________________________

/*      Input               CuSys   : CubicSys (this)
                            eidx    : Error index
                            noret   : Flag for linefeed (0=linefeed)
                            pname   : string in message
        Output              none    : Error message output
                                      Program execution stopped if fatal     */

void          CubicSys::CubicError(int eidx, int noret=0) const;
void          CubicSys::CubicError(int eidx, const std::string& pname, int noret=0) const;
volatile void CubicSys::CubicFatal(int eidx, const std::string& pname) const;

// ____________________________________________________________________________
// ii               CLASS CUBIC SYSTEM DEALINGS WITH ION LISTS
// ____________________________________________________________________________

/* This function sets up a static list of ions that are available for being
   a cubic system.  The list is shared by all cubic systems delcared in a GAMMA
   program. A non-zero value of the ion number (NIso) indicates that the list
   has been initialized so that one may avoid repeated initializations. The
   ion array, CubicIons, is first allocated and set to contain empty ions. Then
   we loop over the available ions and fill up CubicIons with appropriate
   values.                                                                   */

void CubicSys::SetIonList();


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________
// A                    CUBIC SYSTEM CONSTRUCTION, ASSIGNMENT
// ____________________________________________________________________________

/* These constructors set up a new cubic system.  There are few ways to make
   a Cubic System:

         Input Arguments                   Resulting System
      ---------------------    -------------------------------------------
               -               Empty Cubic System, No J (J=0)
            symbol             Cubic System With Specified Ion (e.g. Ce3+)

    The destructor here does nothing except delete the ion index because the
    ions list, although using some small amount of memory, must remain until
    all Cubic Systems in a  program are deleted.  This will not be done
    until program completion.                                                */

MSVCDLC      CubicSys::CubicSys();
MSVCDLC      CubicSys::CubicSys(const  CubicSys& I);
MSVCDLC      CubicSys::CubicSys(const  std::string& symbol);
MSVCDLL void CubicSys::operator=(const CubicSys& I);
MSVCDLC      CubicSys::~CubicSys();

// ____________________________________________________________________________
// B                       CUBIC SYSTEM ACCESS FUNCTIONS
// ____________________________________________________________________________

/*    Function        Return      Units                Example(s)
      --------     -------------  -----   -------------------------------
         qn        mz   (double)   hbar   4.0, 5.5, 9.5....
         HS        2I+1 (int)      none   10, 11, 16, .....
      momentum     mz   (string)   none   4, 5/2, 9/2, .....
       symbol           (string)   none   Ce3+, Pr3+, Nd3+, ....
       name             (string)   none   Cerium, Praseodymium, Ytterbium, ...
      element           (string)   none   Ce, Pr, Nd, Sm, ......
       number           (int)      none   58<-Ce, 59<-Pr, 60<-Nd, .....
       mass             (int)      amu    133<-Ce, 141<-Pr, 144<-Nd, ....
      weight            (double)   g/m    140.12<-Ce, 140.9077<-Pr, ....

	   Input                CuSys   : CubicSys (this)
	   Note			  		  		: All functions use CubicIons (ion data)  */

MSVCDLL       double       CubicSys::qn()       const;
MSVCDLL       int          CubicSys::HS()       const;
MSVCDLL       int          CubicSys::charge()   const;
MSVCDLL       std::string  CubicSys::momentum() const;
MSVCDLL const std::string& CubicSys::symbol()   const;
MSVCDLL const std::string& CubicSys::name()     const;
MSVCDLL const std::string& CubicSys::element()  const;
MSVCDLL       int          CubicSys::number()   const;
MSVCDLL       int          CubicSys::mass()     const;
MSVCDLL       double       CubicSys::weight()   const;
MSVCDLL       double       CubicSys::gJ()       const;
MSVCDLL       double       CubicSys::beta()     const;
MSVCDLL       double       CubicSys::gamma()    const;
MSVCDLL       double       CubicSys::F4()       const;
MSVCDLL       double       CubicSys::F6()       const;

// ____________________________________________________________________________
// C                        CUBIC SYSTEM I/O FUNCTIONS
// ____________________________________________________________________________

/* These functions allow users to write cubic systems in formatted ASCII to an
   output stream.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
     print     ostream     Writes cubic system in ASCII to output stream
      <<       ostream     Writes cubic system in ASCII to output stream     */

MSVCDLL        std::ostream& CubicSys::print(std::ostream& ostr) const;
MSVCDLL friend std::ostream& operator<<     (std::ostream& ostr, const CubicSys& CS);

// ____________________________________________________________________________
// D                         CUBIC SYSTEM LIST FUNCTIONS
// ____________________________________________________________________________

/* This class contains a static list of possible magnetic ions that can be the
   basis of a cubic system. The functions herein allow for manipulation of the
   magnetic ion list.

     Function        Output                       Description
   ------------  -------------  ----------------------------------------
       seek        ion index    Index in list of magnetic ion input
      exists       ion symbol   True if ion specified is present in list
       known       ion symbol   True if ion specified is present in list

   The seek function will return -1 if the ion is not found in the list.     */

MSVCDLL        int  CubicSys::seek(const   MagIon& MI);
MSVCDLL        bool CubicSys::exists(const std::string& symbol);
MSVCDLL static bool CubicSys::known(const  std::string& symbol);

// ____________________________________________________________________________
// E                     CubicSys Container Support Functions
// ____________________________________________________________________________

/* Aside for providing basic tests as to whether two systems are equvalent or
   not, these operators are necessary if any STL container classes are to
   be used based on cubic systems (e.g. list<CubicSys> or vector<CubicSys>)  */

MSVCDLL bool CubicSys::operator==(const CubicSys& CuSys) const;
MSVCDLL bool CubicSys::operator!=(const CubicSys& CuSys) const;
MSVCDLL bool CubicSys::operator<(const  CubicSys& CuSys) const;
MSVCDLL bool CubicSys::operator>(const  CubicSys& CuSys) const;

  };

#endif								// CubicSys.h

