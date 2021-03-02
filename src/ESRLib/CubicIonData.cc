/* CubicIonData.cc **********************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**  Cubic Ion Data				    Implementation	**
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
** Class CubicIonData tracks informattion concerning a single ion that  **
** has spin angular momentum resulting primarily from its cubically     **
** symmetrhic (octahedral) environment. Each object of this type        **
** contains the specifics of a single cubic ion for use in electron     **
** magnetic resonance computations.  Some access functions & output     **
** constants are also provided.                                         **
**                                                                      **
** Note that users will rarely (if ever) deal with this file.  GAMMA    **
** contains a linked list of most known cubic ions, each element of the **
** list being an object of type CubicIonData. An individual Cubic Ion   **
** will point to an entry of this list (and hence the appropriate data) **
** as well as to an entry into the more basic class IsotopeDataList.    **
**                                                                      **
*************************************************************************/

#ifndef   CubicIonData_cc_		// Is this file already included?
#  define CubicIonData_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <ESRLib/CubicIonData.h>	// Include the interface
#include <Basics/StringCut.h>		// Include string cutting 

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A               Cubic Ion Constructors, Destructor, Assignment
// ____________________________________________________________________________
 
/* These constructors set up a new cubic ion.  There are a few ways to make
   a cubic ion:
 
         Input Arguments                   Resulting System
      ---------------------    -------------------------------------------
               -               Empty cubic ion
               CID              Duplicate cubic ion
       HS,SY,NA,EL,NU,MA,WT    Magnetic ion with properties set by input

    The destructor here does nothing except delete the ion index because the
    ions list, although using some small amount of memory, must remain until
    all Cubic Systems in a  program are deleted.  This will not be done
    until program completion.                                                */
    
CubicIonData::CubicIonData() 
  {
  _symbol = "";
  _charge = 0;
  _gJ     = 0;
  _beta   = 0;
  _gamma  = 0;
  }

CubicIonData::CubicIonData(const CubicIonData& CID) 
  {
  _symbol = CID._symbol;
  _charge = CID._charge;
  _gJ     = CID._gJ;
  _beta   = CID._beta;
  _gamma  = CID._gamma;
  }

CubicIonData::CubicIonData(const std::string& CID) 
  { 
  _symbol = CID;
  _charge = 0;
  _gJ     = 0;
  _beta   = 0;
  _gamma  = 0;
  }

void CubicIonData::operator= (const CubicIonData& CID)
  {
  _symbol = CID._symbol;
  _charge = CID._charge;
  _gJ     = CID._gJ;
  _beta   = CID._beta;
  _gamma  = CID._gamma;
  }

CubicIonData::~CubicIonData() {}
 
// ____________________________________________________________________________
// B                        Cubic Ion Access Functions
// ____________________________________________________________________________
 
/*    Function        Return      Units                Example(s)
      --------     -------------  -----   -------------------------------
       charge          (int)      none    3, -2, 5, .....
         gJ          (double)     none
        beta         (double)     none
       gamma         (double)     none                                       */

std::string  CubicIonData::symbol()  const { return _symbol; }
int          CubicIonData::charge()  const { return _charge; }
double       CubicIonData::gJ()      const { return _gJ;     }
double       CubicIonData::beta()    const { return _beta;   }
double       CubicIonData::gamma()   const { return _gamma;  }

// ____________________________________________________________________________
// C                       Cubic Ion I/O Functions
// ____________________________________________________________________________

/* These functions allow users to write cubic ions in formatted ASCII to an
   output stream.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
     print     ostream     Writes cubic system in ASCII to output stream
      <<       ostream     Writes cubic system in ASCII to output stream     */

std::ostream& CubicIonData::print(std::ostream& ostr, int lf) const
  {
  ostr << "\nCubic Ion     " << _symbol
       << "\n Charge       " << _charge
       << "\n gJ Value     " << _gJ
       << "\n Beta         " << _beta
       << "\n Gammma       " << _gamma;
  if(lf) ostr << "\n";
  return ostr;
  }

std::ostream& operator<< (std::ostream& ostr, const CubicIonData& CID) 
  { return CID.print(ostr); }

#endif							// CubicIonData.cc
