/* CubicIon.h *****************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**	Magnetic Ion 		           	Interface Definition	**
**							 	 	**
**  Copyright (c) 2000							**
**  Scott A. Smith				 			**
**  National High Magnetic Field Laboratory                           	**
**  1800 E. Paul Dirac Drive                                          	**
**  Box 4005                                                          	**
**  Tallahassee, FL 32306                                             	**
**						 			**
**  $Header: $
**						 			**
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
** Class CubicIon represents a magnetic ion in a crystal lattice with   **
** cubic symmetry. The cubic ion will be associated with a specific     **
** ionized element (e.g. Ce3+) exhibiting intrinsic spin (e.g. J=7/2).  **
** These may be used in spin systems and to generate Hamiltonians and   **
** spectra in EPR computations. The class also manages a list of all    **
** possible cubic ions intrinsically supported.                         **
**                                                                      **
*************************************************************************/

#ifndef   CubicIon_h_				// Is file already included?
# define CubicIon_h_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)			// If it is the GNU compiler
#    pragma interface				// then this is the interface
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <ESRLib/CubicIonData.h>		// Include Cubic Ion Data
#include <Basics/Isotope.h>			// Inlcude GAMMA Isotopes
#include <fstream>				// Include file stream handling
#include <vector>				// Include libstdc++ STL vectors

class CubicIon : public Isotope
  {
  static std::vector<CubicIonData> CubicIons;	// Pointer to cubic ions list
         int			   Ion;		// An ion index

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                    CLASS CUBIC ION ERROR HANDLING
// ____________________________________________________________________________

/*      Input               CuSys   : CubicIon (this)
                            eidx    : Error index
                            noret   : Flag for linefeed (0=linefeed)
                            pname   : string in message
        Output              none    : Error message output
                                      Program execution stopped if fatal     */

void          CIError(int eidx,                           int noret=0) const;
void          CIError(int eidx, const std::string& pname, int noret=0) const;
volatile void CIFatal(int eidx, const std::string& pname) const;

// ____________________________________________________________________________
// ii           CLASS CUBIC ION DEALINGS WITH CUBIC ION LISTS
// ____________________________________________________________________________

/* This function sets up a static list of ions that are available for being
   a cubic system.  The list is shared by all cubic systems delcared in a GAMMA
   program. A non-zero value of the ion number (NIso) indicates that the list
   has been initialized so that one may avoid repeated initializations. The
   ion array, CubicIons, is first allocated and set to contain empty ions. Then
   we loop over the available ions and fill up CubicIons with appropriate
   values.                                                                   */

void AddCubicIons();

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________
// A                    CUBIC ION CONSTRUCTION, ASSIGNMENT
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

MSVCDLC      CubicIon();
MSVCDLC      CubicIon(const       CubicIon& CI);
MSVCDLC      CubicIon(const  std::string&   CI);
MSVCDLL void operator=(const      CubicIon& CI);
MSVCDLC      virtual ~CubicIon();

// ____________________________________________________________________________
// B                       CUBIC ION ACCESS FUNCTIONS
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

	   Input                CuSys   : CubicIon (this)
	   Note			  		  		: All functions use CubicIons (ion data)  */

//       double       qn()       const;		// INHERITED
//       int          HS()       const;		// INHERITED
//       std::string  momentum() const;		// INHERITED
// const std::string& symbol()   const;		// INHERITED
// const std::string& name()     const;		// INHERITED
// const std::string& element()  const;		// INHERITED
//       int          number()   const;		// INHERITED
//       int          mass()     const;		// INHERITED
//       double       weight()   const;		// INHERITED

MSVCDLL int    charge() const;
MSVCDLL double gJ()     const;
MSVCDLL double beta()   const;
MSVCDLL double gamma()  const;
MSVCDLL double F4()     const;
MSVCDLL double F6()     const;

// ____________________________________________________________________________
// C                        CUBIC ION I/O FUNCTIONS
// ____________________________________________________________________________

/* These functions allow users to write cubic systems in formatted ASCII to an
   output stream.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
     print     ostream     Writes cubic system in ASCII to output stream
      <<       ostream     Writes cubic system in ASCII to output stream     */

MSVCDLL        std::ostream& print(std::ostream& ostr) const;
MSVCDLL friend std::ostream& operator<<     (std::ostream& ostr, const CubicIon& CS);

// ____________________________________________________________________________
// D                         CUBIC ION LIST FUNCTIONS
// ____________________________________________________________________________

/* This class contains a static list of possible cubic ions that can be the
   basis of a cubic system. The functions herein allow for manipulation of the
   cubic ion list.

     Function        Output                       Description
   ------------  -------------  ----------------------------------------
       seek        ion index    Index in list of cubic ion input
      exists       ion symbol   True if ion specified is present in list
       known       ion symbol   True if ion specified is present in list
    initialize       void       Adds cubic ions to GAMMA isotopes list

   The seek function will return -1 if the ion is not found in the list.     */

MSVCDLL        int  seek(const   CubicIonData& CI);
MSVCDLL        bool exists(const std::string&  symbol);
MSVCDLL static bool known(const  std::string&  symbol);
MSVCDLL static void initialize();
MSVCDLL static int  size();
MSVCDLL static void PrintList(std::ostream& ostr, bool hdr=true);

// ____________________________________________________________________________
// E                     CubicIon Container Support Functions
// ____________________________________________________________________________

/* Aside for providing basic tests as to whether two systems are equvalent or
   not, these operators are necessary if any STL container classes are to
   be used based on cubic systems (e.g. list<CubicIon> or vector<CubicIon>)  */

MSVCDLL bool operator==(const CubicIon& CuSys) const;
MSVCDLL bool operator!=(const CubicIon& CuSys) const;
MSVCDLL bool operator<(const  CubicIon& CuSys) const;
MSVCDLL bool operator>(const  CubicIon& CuSys) const;

};

#endif								// CubicIon.h

