/* CubicSys.cc **************************************************-*-c++-*-
**									**
** 	                   G A M M A 					**
**									**
**	CubicSys 			          Implementation	**
**						 			**
**  Copyright (c) 2000							**
**  Scott A. Smith							**
**  National High Magnetic Field Laboratory                           	**
**  1800 E. Paul Dirac Drive                                          	**
**  Box 4005                                                          	**
**  Tallahassee, FL 32306                                             	**
**							 		**
**      $Header: $
**							 		**
*************************************************************************/

/*************************************************************************
**						 			**
**  Description							 	**
**						 			**
** Class CubicSys represents a magnetic ion in a crystal lattice with	**
** cubic symmetry. The system will be associated with a specific	**
** ionized element (e.g. Ce3+) exhibiting intrinsic spin (e.g. J=7/2).	**
** This system will be used in generating Hamiltonians and spectra in   **
** EPR computations. The class also manages a list of all possible ions	**
** intrinsically supported.						**
**									**
*************************************************************************/
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

#ifndef   CubicSys_cc_			// Is the file already included?
#  define CubicSys_cc_ 1		// If no, then remember it 
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <ESRLib/CubicSys.h>		// Include the interface
#include <ESRLib/MagIon.h>		// Include list magnetic ions
#include <Basics/StringCut.h>		// Include string cutting
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <Basics/Gconstants.h>		// Include GAMMA constants
#include <string>			// Include string class
#include <fstream>			// Include I/O streams

int     CubicSys::NIons     = 0;	// Initialize to No Ions
MagIon* CubicSys::CubicIons = NULL;	// Initalize to NULL list

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                     CLASS CUBIC SYSTEM ERROR HANDLING
// ____________________________________________________________________________
 
/*      Input               CuSys   : CubicSys (this)
                            eidx    : Error index
                            noret   : Flag for linefeed (0=linefeed)
                            pname   : string in message
        Output              none    : Error message output
                                      Program execution stopped if fatal

The following error messages use the defaults set in the Gutils package
 
               Case                          Error Message
 
 NO PNAME       (0)                     Program Aborting.....
    		(9)                     Problems During Construction
                default                 Unknown Error 
 WITH PNAME     (1)                     Problems With File PNAME
                (2)                     Cannot Read Parameter PNAME
    		(3)  			Invalid Use Of Function PNAME
                default                 Unknown Error - PNAME                */

void CubicSys::CubicError(int eidx, int noret) const
  {
  std::string hdr("CubicSys");
  switch(eidx)
    {
    case 0:  GAMMAerror(hdr, 0, noret);    break;	// Program Aborting (0)
    default: GAMMAerror(hdr, eidx, noret); break;	// Unknown Error    (-1)
    }
  }

void CubicSys::CubicError(int eidx, const std::string& pname, int noret) const
  {
  string hdr("CubicSys");
  string msg;
  switch(eidx)
    {
    case 2: msg=string("Attempted Access of Unknown CubicSys Ion ")+pname;// (2)
      break;
    default: msg=string("Unknown Error - ") + pname; break;
    }
  GAMMAerror(hdr, msg, noret);
  }

volatile void CubicSys::CubicFatal(int eidx, const string& pname) const
  {                                                                 
  CubicError(eidx, pname, eidx);			// Output error message
  if(eidx) CubicError(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }
 
// ____________________________________________________________________________
// ii              CLASS CUBIC SYSTEM DEALINGS WITH ION LISTS
// ____________________________________________________________________________

/* This function sets up a static list of ions that are available for being
   a cubic system.  The list is shared by all cubic systems delcared in a GAMMA
   program. A non-zero value of the ion number (NIons) indicates that the list
   has been initialized so that one may avoid repeated initializations. The
   ion array, CubicIons, is first allocated and set to contain empty ions. Then
   we loop over the available magnetic ions and fill up CubicIons with
   appropriate values.                                                       */

void CubicSys::SetIonList()
  {
  NIons = 11;					// Set # magnetic ions we know 
  int i=0;					// Temp looping variable
  CubicIons = new MagIon[11];			// Allocate array for 'em
  MagIon* MI = CubicIons;			// Pointer to 1st magnetic ion
  for(; i<NIons; i++, MI++)			// Loop over all magnetic ions
    MI = new MagIon(); 				// and set them to default

/*        An array of known magnetic ion spin Hilbert spaces (2*I+1)         */

  int Spins[11] = {  6,  9, 10,  9,  6, 13, 16, 17, 16, 13, 8 };
  
/*                        An array of magnetic ion charges                  */
   
  int Charges[11] = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

/*     An array of known magnetic ion atomic numbers (58=Ce, 59=Pr, ...)     */

  int Numbers[11] = { 58, 59, 60, 61, 62, 65, 66, 67, 68, 69, 70 };

/*   An array of known magnetic ion atomic masses  (133=Ce, 141=Pr, ...)     */

  int Masses[11] = { 140, 141, 144, 145, 150, 159, 162, 165, 168, 169, 173 };

/*  Array of known magnetic ion ave. at. weights (Ce = 140.12 g/mole, ...)   */

  double Weights[11] = { 140.12,   140.9077, 144.24,   145.0,  150.36, 
                         158.9245, 160.50,   164.9304, 167.26, 168.9342,
                         173.04 };

/*                  An array of known magnetic ion names                     */

  string Names[11] = { "Cerium",  "Praseodymium", "Neodymium", 
                                                  "Promethium", "Samarium",
                       "Terbium", "Dysprosium",   "Holmium",    "Erbium",
                                                  "Thulium",    "Ytterbium" };

/*                      Array of magnetic ion elements                       */

  string Elements[11] = { "Ce", "Pr", "Nd", "Pm", "Sm",
                          "Tb", "Dy", "Ho", "Er", "Tm", "Yb" };
  
  double gJValues[11] = { 6.0/7.0, 4.0/5.0, 8.0/11.0, 3.0/5.0, 2.0/7.0,
                          3.0/2, 4.0/3.0, 5.0/4.0,  6.0/5.0, 7.0/6.0, 8.0/7.0 };
  
  double betaValues[11] = { 6.3492e-3, -7.3462e-4, -2.9111e-4,  4.0755e-4,
                            2.5012e-3,  1.2244e-4, -5.9200e-5, -3.3300e-5,
                            4.4400e-5,  1.6325e-4,  1.7316e-3 };
                            
  double gammaValues[11] = { 0.0,        6.0994e-5, -3.7988e-5,  6.0781e-5,
                             0.0,       -1.1212e-6,  1.0350e-6, -1.2937e-6,
                             2.0699e-6, -5.6061e-6,  1.4800e-4 };

  for(i=0, MI=CubicIons; i<NIons; i++, MI++)	// Loop all magnetic ions
    {
    MI->_HS          = Spins[i]; 		// Set spin Hilbert space
    MI->_number      = Numbers[i]; 		// Set the atomic number
    MI->_charge      = Charges[i];		// Set the ion charges
    MI->_mass        = Masses[i]; 		// Set the atomic mass
    MI->_weight      = Weights[i]; 		// Set the atomic weights
    MI->_name        = Names[i];		// Set the magnetic ion name 
    MI->_element     = Elements[i];		// Set the magnetic ion element
    MI->_gJ          = gJValues[i];		// Set the g scaling factor
    MI->_beta        = betaValues[i];           // Set the beta value
    MI->_gamma       = gammaValues[i];          // Set the gamma value
    MI->_symbol      = MI->_element 
                     + Gdec(MI->_charge);
    if(MI->_charge > 0) MI->_symbol += "+";
    if(MI->_charge < 0) MI->_symbol += "-";
    }
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A                     CUBIC SYSTEM CONSTRUCTION, ASSIGNMENT
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
 

CubicSys::CubicSys()
  {
  MagIon Cerium("Ce3+");		// Set mock default ion, Ce3+
  if(!CubicIons) SetIonList(); 		// If ions list is empty, must set
  Ion = seek(Cerium);			// Get index in ions list for Ce
  if(Ion < 0) CubicFatal(2, "Ce3+"); 	// If no Ce3+, FATAL ERROR!
  }

CubicSys::CubicSys(const string& symbol)
  {
  MagIon i(symbol);			// Form data for this magnetic ion 
  if(!CubicIons) SetIonList(); 		// If ions list is empty, must set
  Ion = seek(i);			// Get index into magnetic ions list
  if(Ion < 0) CubicFatal(2, symbol); 	// Fatal error, unknown magnetic ion
  
  }
  
     CubicSys::CubicSys(const CubicSys& I) : Ion(I.Ion) {}
void CubicSys::operator= (const CubicSys& CuSys) { Ion = CuSys.Ion; }
     CubicSys::~CubicSys() {};
 
// ____________________________________________________________________________
// B                        CUBIC SYSTEM ACCESS FUNCTIONS
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
      weight            (double)   g/m    140.12<-Ce, 140.9077<-Pr, ....

	Input		CuSys   : CubicSys (this)
	Note			: All functions use CubicIons (ion data)     */
	   
       double CubicSys::qn()       const { return (CubicIons[Ion]).qn(); }
       int    CubicSys::HS()       const { return (CubicIons[Ion]).HS(); }
       int    CubicSys::charge()   const { return (CubicIons[Ion]).charge(); }
       string CubicSys::momentum() const { return (CubicIons[Ion]).momentum();}
const string& CubicSys::symbol()   const { return (CubicIons[Ion]).symbol(); }
const string& CubicSys::name()     const { return (CubicIons[Ion]).name(); }
const string& CubicSys::element()  const { return (CubicIons[Ion]).element(); }
      int     CubicSys::number()   const { return (CubicIons[Ion]).number(); }
      int     CubicSys::mass()     const { return (CubicIons[Ion]).mass(); }
      double  CubicSys::weight()   const { return (CubicIons[Ion]).weight(); }
      double  CubicSys::gJ()       const { return (CubicIons[Ion]).gJ(); }
      double  CubicSys::beta()     const { return (CubicIons[Ion]).beta(); }
      double  CubicSys::gamma()    const { return (CubicIons[Ion]).gamma(); }
      
      double  CubicSys::F4() const
        {
        switch(HS())
          {
          case  5: return 12; break;
          case  6: return 60; break;			// Ce3+, Sm3+
          case  7: return 15; break;
          case  8:                                      // Yb3+
          case  9:					// Pr3+, Pm3+
          case 10: return 60; break;			// Nd3+
          case 11: return 14; break;
          case 12:
          case 13:					// Tb3+, Tm3+
          case 14:
          case 15:
          case 16:					// Er3+, Dy3+
          case 17: return 60; break;			// Ho3+
          default: return 0;
          }
        }

      double  CubicSys::F6() const
        {
        switch(HS())
          {
          case  5:
          case  6: return    0; break;			// Ce3+, Sm3+
          case  7: return  180; break;
          case  8:  					// Yb3+
          case  9: return 1260; break;			// Pr3+, Pm3+
          case 10: return 2520; break;			// Nd3+
          case 11: return 1260; break;
          case 12: return 3780; break;
          case 13:					// Tb3+, Tm3+
          case 14: return 7560; break;
          case 15: return 3780; break;
          case 16:					// Er3+, Dy3+
          case 17: return 13860; break;			// Ho3+
          default: return 0; 
          }
        }
 
// ____________________________________________________________________________
// C                        CUBIC SYSTEM I/O FUNCTIONS
// ____________________________________________________________________________
 
/* These functions allow users to write cubic systems in formatted ASCII to an
   output stream.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
     print     ostream     Writes cubic system in ASCII to output stream
      <<       ostream     Writes cubic system in ASCII to output stream     */

ostream& CubicSys::print(ostream& ostr) const
  {
  CubicIons[Ion].print(ostr, 0);
  ostr << "\n More Cubic System Info Here       \n";
  return ostr; 
  }
 
ostream& operator<< (ostream& ostr, const CubicSys& I) { return I.print(ostr); }
 
// ____________________________________________________________________________
// D                           CUBIC SYSTEM LIST FUNCTIONS
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
   
int CubicSys::seek(const MagIon& I)
  {
  if(!CubicIons) return -1;			// If no ions list we fail
  for(int i=0; i<NIons; i++)
    if(I.symbol() == CubicIons[i].symbol()) return i;
  return -1;
  }

bool CubicSys::exists(const string& symbol)
  {
  MagIon I(symbol);			// Get magnetic ion data for symbol
  if(!CubicIons) SetIonList(); 		// If CubicIons is empty, must
  if(seek(I)>=0) return true;
  return false;
  }

bool CubicSys::known(const string& symbol)
  { CubicSys X; return X.exists(symbol); }

// ____________________________________________________________________________
// E                    CubicSys Container Support Functions
// ____________________________________________________________________________

/* Aside for providing basic tests as to whether two systems are equvalent or
   not, these operators are necessary if any STL container classes are to
   be used based on cubic systems (e.g. list<CubicSys> or vector<CubicSys>)  */
   
bool CubicSys::operator== (const CubicSys& CS) const { return (Ion==CS.Ion);  }
bool CubicSys::operator!= (const CubicSys& CS) const { return (Ion!=CS.Ion);  }
bool CubicSys::operator<(const   CubicSys& CS) const { return (HS()<CS.HS()); }
bool CubicSys::operator>(const   CubicSys& CS) const { return (HS()>CS.HS()); }

#endif								// CubicSys.cc
