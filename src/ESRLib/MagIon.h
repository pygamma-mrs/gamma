/* MagIon.h *****************************************************-*-c++-*-
**									**
** 	                      	  G A M M A				**
**									**
**	Magnetic Ion                                Interface		**
**								 	**
**	Copyright (c) 2000					 	**
**	Scott A. Smith				 			**
**      National High Magnetic Field Laboratory				**
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
**  Class MagIon embodies a single ion having spin angular momentum     **
**  that results primarily from its environment. Each object of this  	**
**  type contains all the required data corresponding to a single ion	**
**  for use in electron magnetic resonance computations.  Some access   **
**  functions & output constants are provided.                          **
**                                                                      **
**  Note that users will rarely (if ever) deal with this file.  GAMMA   **
**  contains a linked list of most known ions with spins, each element  **
**  of the list being an object of type MagIon.  Higher classes and     **
**  GAMMA based programs (such as CubicSys) use a linked list of ions.  **
**                                                                      **
**									**
*************************************************************************/

#ifndef   MagIon_h_			// Is this file already included?
#  define MagIon_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <string>			// Include string manipulations

class MagIon
  {
  friend class CubicSys;		// Allow Cubic Systems full access here
  int         _HS;			// Spin 2*I+1 values (11 if I=5/2)
  int         _charge;			// Ion charge (3, 6, -2, ...)
  std::string _symbol;			// Ion symbol (Ce3+, Pr3+, Nd3+, ...)
  std::string _name;			// Ion name (Cerium, , ..)
  std::string _element;			// Element type (Ce, Pr, ...)
  int         _number;			// Periodic number (58 for Ce, 68 for Er)
  int         _mass;			// Average atomic mass in amu (140 for Ce)
  double      _weight;			// Average atomic weight (in grams/mole)
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
// A              Magnetic Ion Constructors, Assignment,Destructors
// ____________________________________________________________________________

/* These constructors set up a new magnetic ion.  There are a few ways to make
   a magnetic ion:

         Input Arguments                   Resulting System
      ---------------------    -------------------------------------------
               -               Empty magnetic ion
              MI               Duplicate magnetic ion
       HS,SY,NA,EL,NU,MA,WT    Magnetic ion with properties set by input

    The constructor taking a full set of data argumnents allow users to
    make up their own magnetic ions "on the fly" if necessary. The constructor
    which demands only a symbol value exists for convenience in other
    applications but may later cause troubles if other values left zero.,   */

MSVCDLC      MagIon::MagIon();
MSVCDLC      MagIon::MagIon(const MagIon& MI);
MSVCDLC      MagIon::MagIon(int HS, const std::string& SY, int CH, const std::string& NA,
                              const std::string& EL, int NU, int MA, double WT);
MSVCDLC      MagIon::MagIon(const std::string& SYM);
MSVCDLL void MagIon::operator= (const MagIon& ID1);
MSVCDLC      MagIon::~MagIon();

// ____________________________________________________________________________
// B                        Magnetic Ion Access Functions
// ____________________________________________________________________________

/*    Function        Return      Units                Example(s)
      --------     -------------  -----   -------------------------------
         qn        mz   (double)   hbar   4.0, 5.5, 9.5....
         HS        2I+1 (int)      none   10, 11, 16, .....
       charge            (int)     none   3, -2, 5, .....
      momentum     mz   (string)   none   4, 5/2, 9/2, .....
       symbol           (string)   none   Ce3+, Pr3+, Nd3+, ....
       name             (string)   none   Cerium, Praseodymium, Ytterbium, ...
      element           (string)   none   Ce, Pr, Nd, Sm, ......
       number           (int)      none   58<-Ce, 59<-Pr, 60<-Nd, .....
       mass             (int)      amu    140<-Ce, 141<-Pr, 144<-Nd, ....
      weight            (double)   g/m    140.12<-Ce, 140.9077<-Pr, ....
      g-scaling			(double)   none   6/7<-Ce3+, 4/5<-Pr3+,....         */

MSVCDLL        double       MagIon::qn()	   const;		// Calculated
MSVCDLL        int          MagIon::HS()	   const;		// Stored
MSVCDLL        int          MagIon::charge()   const;       // Stored
MSVCDLL 	   std::string  MagIon::momentum() const;		// Calculated
MSVCDLL const  std::string& MagIon::symbol()   const;		// Stored
MSVCDLL const  std::string& MagIon::name()	   const;		// Stored
MSVCDLL const  std::string& MagIon::element()  const;		// Stored
MSVCDLL        int          MagIon::number()   const;		// Stored
MSVCDLL        int          MagIon::mass()	   const;		// Stored
MSVCDLL        double       MagIon::weight()   const;		// Stored
MSVCDLL        double       MagIon::gJ()       const;		// Stored
MSVCDLL        double       MagIon::beta()     const;		// Stored
MSVCDLL        double		MagIon::gamma()    const;		// Stored

// ____________________________________________________________________________
// C                         Magnetic Ion I/O Functions
// ____________________________________________________________________________

/* These functions allow users to write magnetic ions in formatted ASCII to an
   output stream.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
     print     ostream     Writes cubic system in ASCII to output stream
      <<       ostream     Writes cubic system in ASCII to output stream     */

MSVCDLL std::ostream& MagIon::print(std::ostream& ostr, int lf=1) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const MagIon& ID);
};

#endif															// MagIon.h

