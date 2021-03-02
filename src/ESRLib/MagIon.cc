/* MagIon.cc ****************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Magnetic Ion                                  Implementation	**
**                                                                      **
**  Copyright (c) 2000							**
**  Scott A. Smith							**
**  National High Magnetic Field Laboratory                           	**
**  1800 E. Paul Dirac Drive                                          	**
**  Box 4005                                                          	**
**  Tallahassee, FL 32306                                             	**
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/
     
/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
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
*************************************************************************/

#ifndef   MagIon_cc_			// Is this file already included?
#  define MagIon_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <ESRLib/MagIon.h>		// Include the interface
#include <Basics/StringCut.h>		// Include string cutting

using std::ostream;			// Using libstdc++ output streams
using std::string;			// Using libstdc++ strings

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A               Magnetic Ion Constructors, Destructor, Assignment
// ____________________________________________________________________________
 
/* These constructors set up a new magnetic ion.  There are a few ways to make
   a magnetic ion:
 
         Input Arguments                   Resulting System
      ---------------------    -------------------------------------------
               -               Empty magnetic ion
               MI              Duplicate magnetic ion
       HS,SY,NA,EL,NU,MA,WT    Magnetic ion with properties set by input

    The destructor here does nothing except delete the ion index because the
    ions list, although using some small amount of memory, must remain until
    all Cubic Systems in a  program are deleted.  This will not be done
    until program completion.                                                */
    
MagIon::MagIon()
  { _HS=0; _symbol=""; _charge=0; _name=""; _element=""; _number=0; _mass=0; _weight=0; }

MagIon::MagIon(const MagIon& MI)
  {
  _HS      = MI._HS;
  _charge  = MI._charge;
  _symbol  = MI._symbol;
  _name    = MI._name;
  _element = MI._element;
  _number  = MI._number;
  _mass    = MI._mass;
  _weight  = MI._weight;
  }

MagIon::MagIon(int HS, const string& SY, int CH, const string& NA,
			 const string& EL, int NU, int MA, double WT )
  {
  _HS      =  HS;
  _symbol  =  SY;
  _name    =  NA;
  _element =  EL;
  _number  =  NU;
  _mass    =  MA;
  _weight  =  WT;
  }

MagIon::MagIon(const string& symb)
 { _HS=0; _symbol=symb; _charge=0; _name=""; _element=""; _number=0; _mass=0; _weight=0; }

 void MagIon::operator= (const MagIon& MI)
  {
  _HS          = MI._HS;
  _charge      = MI._charge;
  _symbol      = MI._symbol;
  _name        = MI._name;
  _element     = MI._element;
  _number      = MI._number;
  _mass        = MI._mass;
  _weight      = MI._weight;
  }

MagIon::~MagIon() {}
 
// ____________________________________________________________________________
// B                        Magnetic Ion Access Functions
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
       mass             (int)      amu    140<-Ce, 141<-Pr, 144<-Nd, ....
      weight            (double)   g/m    140.12<-Ce, 140.9077<-Pr, ....     */

      double  MagIon::qn()      const { return (_HS ? (_HS-1)/2.0 : 0); }
      int     MagIon::HS()      const { return _HS; }
      int     MagIon::charge()  const { return _charge; }
const string& MagIon::symbol()  const { return _symbol; }
const string& MagIon::name()    const { return _name; }
const string& MagIon::element() const { return _element; }
      int     MagIon::number()  const { return _number; }
      int     MagIon::mass()    const { return _mass; }
      double  MagIon::weight()  const { return _weight; }
      double  MagIon::gJ()      const { return _gJ; }
      double  MagIon::beta()    const { return _beta; }
      double  MagIon::gamma()   const { return _gamma; }

string MagIon::momentum() const
  {
  if(!_HS) return string("0"); 		// If no _HS, then it's zero
  else if(_HS%2) return Gdec((_HS-1)/2);// If odd _HS, return fraction
  else return string(Gdec(_HS-1))+"/2";	// If even _HS, return whole int
  }

// ____________________________________________________________________________
// C                       Magnetic Ion I/O Functions
// ____________________________________________________________________________

/* These functions allow users to write magnetic ions in formatted ASCII to an
   output stream.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
     print     ostream     Writes cubic system in ASCII to output stream
      <<       ostream     Writes cubic system in ASCII to output stream     */

ostream& MagIon::print(ostream& ostr, int lf) const
  {
  ostr << "\nIon          " << _symbol
       << "\n Name        " << _name
       << "\n Charge      " << _charge
       << "\n Element     " << _element
       << "\n Number      " << _number
       << "\n Mass        " << _mass
       << "\n Weight      " << _weight
       << "\n Spin        " << qn();
  if(lf) ostr << "\n";
  return ostr;
  }

ostream& operator<< (ostream& ostr, const MagIon& MI) { return MI.print(ostr); }

#endif							// MagIon.cc
