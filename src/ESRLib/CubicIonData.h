/* CubicIonData.h ***********************************************-*-c++-*-
**									**
** 	                      	  G A M M A				**
**									**
**	Cubic Ion Data                           Interface		**
**								 	**
**	Copyright (c) 2000					 	**
**	Scott A. Smith				 			**
**      National High Cubic Field Laboratory				**
**      1800 E. Paul Dirac Drive					**
**      Box 4005							**
**      Tallahassee, FL 32306						**
**						 			**
**      $Header:$
**								 	**
*************************************************************************/

/*************************************************************************
**							 		**
**  Description						 		**
**							 		**
** Class CubicIonData tracks informattion concerning a single ion that	**
** has spin angular momentum resulting primarily from its cubically	**
** symmetrhic (octahedral) environment. Each object of this type 	**
** contains the specifics of a single cubic ion for use in electron	**
** magnetic resonance computations.  Some access functions & output	**
** constants are also provided.  					**
**                                                                      **
** Note that users will rarely (if ever) deal with this file.  GAMMA	**
** contains a linked list of most known cubic ions, each element of the **
** list being an object of type CubicIonData. An individual Cubic Ion	**
** will point to an entry of this list (and hence the appropriate data)	**
** as well as to an entry into the more basic class IsotopeDataList.	**
**									**
*************************************************************************/

#ifndef   CubicIonData_h_		// Is this file already included?
#  define CubicIonData_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Basics/IsotopeData.h>		// Inlcude GAMMA IsotopeData class
#include <string>			// Include libstdc++ strings

class CubicIonData
  {
  friend class CubicIon;		// Allow Cubic Ions full access here
  std::string _symbol;          	// Isotope symbol (1H, 2H, 14N, ...)
  int         _charge;			// Ion charge (3, 6, -2, ...)
  double      _gJ;			// G scaling 1+J(J+1)+S(S+1)+L(L+1)/2J(J+1)
  double      _beta;
  double      _gamma;

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A              Cubic Ion Constructors, Assignment,Destructors
// ____________________________________________________________________________

/* These constructors set up a new cubic ion.  There are a few ways to make
   a cubic ion:

         Input Arguments                   Resulting System
      ---------------------    -------------------------------------------
               -               Empty cubic ion
              MI               Duplicate cubic ion
       HS,SY,NA,EL,NU,MA,WT    Cubic ion with properties set by input

    The constructor taking a full set of data argumnents allow users to
    make up their own cubic ions "on the fly" if necessary. The constructor
    which demands only a symbol value exists for convenience in other
    applications but may later cause troubles if other values left zero.,   */

MSVCDLC      CubicIonData();
MSVCDLC      CubicIonData(const CubicIonData& CID);
MSVCDLC      CubicIonData(const std::string&  CID);
MSVCDLL void operator= (  const CubicIonData& CID);
MSVCDLC      ~CubicIonData();

// ____________________________________________________________________________
// B                        Cubic Ion Access Functions
// ____________________________________________________________________________

/*    Function        Return      Units                Example(s)
      --------     -------------  -----   -------------------------------
       charge          (int)      none    3, -2, 5, .....
         gJ          (double)     none    
        beta         (double)     none    
       gamma         (double)     none                                       */

MSVCDLL std::string symbol() const;			// Stored
MSVCDLL      int    charge() const;			// Stored
MSVCDLL      double gJ()     const;			// Stored
MSVCDLL      double beta()   const;			// Stored
MSVCDLL      double gamma()  const;			// Stored

// ____________________________________________________________________________
// C                         Cubic Ion I/O Functions
// ____________________________________________________________________________

/* These functions allow users to write cubic ions in formatted ASCII to an
   output stream.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
     print     ostream     Writes cubic system in ASCII to output stream
      <<       ostream     Writes cubic system in ASCII to output stream     */

MSVCDLL std::ostream& print(std::ostream& ostr, int lf=1) const;
MSVCDLL friend  std::ostream& operator<< (std::ostream& ostr, const CubicIonData& CID);
};

#endif							// CubicIonData.h 
