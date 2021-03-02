/* CubicIon.cc **************************************************-*-c++-*-
**									**
** 	                          G A M M A 				**
**									**
**	  Cubic Ion 			          Implementation	**
**						 			**
**  Copyright (c) 2000							**
**  Scott A. Smith							**
**  National High Magnetic Field Laboratory                           	**
**  1800 E. Paul Dirac Drive                                          	**
**  Box 4005                                                          	**
**  Tallahassee, FL 32306                                             	**
**							 		**
**  $Header: $
**							 		**
*************************************************************************/

/*************************************************************************
**						 			**
**  Description							 	**
**						 			**
** Class CubicIon represents a magnetic ion in a crystal lattice with	**
** cubic symmetry. The cubic ion will be associated with a specific	**
** ionized element (e.g. Ce3+) exhibiting intrinsic spin (e.g. J=7/2).	**
** This will be used in spin systems and to generate Hamiltonians and   **
** spectra in EPR computations. The class also manages a list of all    **
** possible cubic ions intrinsically supported.  			**
**									**
*************************************************************************/
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

#ifndef   CubicIon_cc_			// Is the file already included?
#  define CubicIon_cc_ 1		// If no, then remember it 
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <ESRLib/CubicIon.h>		// Include the interface
#include <ESRLib/CubicIonData.h>	// Include cubic ion data
#include <Basics/StringCut.h>		// Include string cutting
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <Basics/Gconstants.h>		// Include GAMMA constants
#include <string>			// Include string class
#include <fstream>			// Include I/O streams

#if defined(_MSC_VER) || defined(__SUNPRO_CC)	// If using MSVC++ then we
 #define min(a,b) (((a)<(b)) ? (a):(b))		// we need to define max & min
 #define max(a,b) (((a)>(b)) ? (a):(b))
#endif

/* Note that although CubicIons is a construct of class CubicIon, it must
   still be initialized because it is static. This helps out in during
   library linking (at least in MSVC++ it does                               */

std::vector<CubicIonData> CubicIon::CubicIons;		// Set empty list

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                     CLASS CUBIC ION ERROR HANDLING
// ____________________________________________________________________________
 
/*      Input               CuSys   : CubicIon (this)
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

void CubicIon::CIError(int eidx, int noret) const
  {
  std::string hdr("CubicIon");
  switch(eidx)
    {
    case 0:  GAMMAerror(hdr, 0, noret);    break;	// Program Aborting (0)
    default: GAMMAerror(hdr, eidx, noret); break;	// Unknown Error    (-1)
    }
  }

void CubicIon::CIError(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("CubicIon");
  std::string msg;
  switch(eidx)
    {
    case 2: msg=std::string("Attempted Access of Unknown CubicIon Ion ")+pname;// (2)
      break;
    default: msg=std::string("Unknown Error - ") + pname; break;
    }
  GAMMAerror(hdr, msg, noret);
  }

volatile void CubicIon::CIFatal(int eidx, const std::string& pname) const
  {                                                                 
  CIError(eidx, pname, eidx);			// Output error message
  if(eidx) CIError(0);				// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }
 
// ____________________________________________________________________________
// ii              CLASS CUBIC ION DEALINGS WITH CUBIC ION LISTS
// ____________________________________________________________________________

/* This function performs 2 tasks, related to the fact that each cubic ion is
   made of two parts: 1.) A GAMMA spin isotope and 2.) Added cubic ion info.
   Each GAMMA isotope, and hence each cubic ion, is simply an index in a large
   list of isotopes known to GAMMA. Since cubic ions are special isotopes we
   must add all the ones we wish to have in GAMMA to this isotopes list.  So,
   our first task is to add all cubic ions into the Isotopes list so that they
   are known to all of the platform. In addition to spin isotope information,
   cubic ions know about charge, octahedral symmetry, J factors, etc. This
   information is stored in a list of all know cubic ions in GAMMA, CubicIons.
   Each cubic ion also has a index in this CubicIons list to indicate its
   particular properties. Therefore, the second task we perform is to set up
   the list of cubic ion information for all GAMMA, CubicIons. 

   All of this needs to be done only once in the course of a simulation, so
   we keep a flag "_added" to know if we have added cubic ions to the isotopes
   list and built the cubic ions data list.  Once done, each cubic ion in a
   program simply has the two indices (one for each list) and this sets
   up all constant information appropriate for the ion.

   A final note. The first declaration of a cubic ion should invoke this
   function and thereby set up the cubic ions list automatically. If this is
   done and the Isotopes list has not been initialized this will also be done
   automatcally when the first cubic ion is added to the Isotopes list.      */


void CubicIon::AddCubicIons()			// Add cubic ions to GAMMA
  {						// IsotopesList for global use
  if(CubicIons.size()) return;			// Nothing if already added!
  int NIons = 11;				// Set # cubic ions we know 

/*        An array of known cubic ion spin Hilbert spaces (2*I+1)         */

  int Spins[11] = {  6,  9, 10,  9,  6, 13, 16, 17, 16, 13, 8 };
  
/*                        An array of cubic ion charges                  */
   
  int Charges[11] = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

/*     An array of known cubic ion atomic numbers (58=Ce, 59=Pr, ...)     */

  int Numbers[11] = { 58, 59, 60, 61, 62, 65, 66, 67, 68, 69, 70 };

/*   An array of known cubic ion atomic masses  (133=Ce, 141=Pr, ...)     */

  int Masses[11] = { 140, 141, 144, 145, 150, 159, 162, 165, 168, 169, 173 };

/*  Array of known cubic ion ave. at. weights (Ce = 140.12 g/mole, ...)   */

  double Weights[11] = { 140.12,   140.9077, 144.24,   145.0,  150.36, 
                         158.9245, 160.50,   164.9304, 167.26, 168.9342,
                         173.04 };

/*                  An array of known cubic ion names                     */

  std::string Names[11] = { "Cerium",  "Praseodymium", "Neodymium", 
                                                  "Promethium", "Samarium",
                       "Terbium", "Dysprosium",   "Holmium",    "Erbium",
                                                  "Thulium",    "Ytterbium" };

/*                      Array of cubic ion elements                       */

  std::string Elements[11] = { "Ce", "Pr", "Nd", "Pm", "Sm",
                          "Tb", "Dy", "Ho", "Er", "Tm", "Yb" };
  
  double gJValues[11] = { 6.0/7.0, 4.0/5.0, 8.0/11.0, 3.0/5.0, 2.0/7.0,
                          3.0/2, 4.0/3.0, 5.0/4.0,  6.0/5.0, 7.0/6.0, 8.0/7.0 };
  
  double betaValues[11] = { 6.3492e-3, -7.3462e-4, -2.9111e-4,  4.0755e-4,
                            2.5012e-3,  1.2244e-4, -5.9200e-5, -3.3300e-5,
                            4.4400e-5,  1.6325e-4,  1.7316e-3 };
                            
  double gammaValues[11] = { 0.0,        6.0994e-5, -3.7988e-5,  6.0781e-5,
                             0.0,       -1.1212e-6,  1.0350e-6, -1.2937e-6,
                             2.0699e-6, -5.6061e-6,  1.4800e-4 };

  Isotope Iele("e-");				// Default electron spin
  bool   elect;
  int    hs, atnum, mass,    charg;
  std::string symb,  name,    ele;
  double weight,  recept,  relfreq;
  double gJval,   betaval, gammaval;
  for(int i=0; i<NIons; i++)			// Loop all cubic ions
    {
    hs          = Spins[i]; 			// Get spin Hilbert space
    atnum       = Numbers[i]; 			// Get the atomic number
    mass        = Masses[i]; 			// Get the atomic mass
    weight      = Weights[i]; 			// Get the atomic weights
    name        = Names[i];			// Get the cubic ion name 
    ele         = Elements[i];			// Get the cubic ion element
    recept      = Iele.receptivity();		// Set the receptivity as e-
    relfreq     = Iele.relative_frequency();	// Set the relative freq. as e-
    elect       = true;				// Set it is an electron
    charg       = Charges[i];			// Get the ion charge 
    gJval       = gJValues[i];			// Get the g scaling factor
    betaval     = betaValues[i];		// Get the beta value
    gammaval    = gammaValues[i];		// Get the gamma value
    symb        = ele + Gdec(charg);		// Set the symbol to reflect
    if(charg > 0) symb += "+";			// the charge on the ion
    if(charg < 0) symb += "-";

    IsotopeData ID(hs,     symb,    name,  ele,	// Construct an IsotopeData
                   atnum,  mass,    weight, 
                   recept, relfreq, elect);
    Isotope::AddIsotope(ID);			// Add to GAMMA isotope list
    CubicIonData CD;				// Cubic ion data
    CD._symbol = symb;
    CD._charge = charg;				// Set the ion charge
    CD._gJ     = gJval;				// Set the g scaling factor
    CD._beta   = betaval;			// Set the beta value
    CD._gamma  = gammaval;			// Set the gamma value
    CubicIons.push_back(CD);			// Add ion to cubic ion list
    }
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A                    CUBIC ION CONSTRUCTION, ASSIGNMENT
// ____________________________________________________________________________

/* These constructors set up a new cubic ion.  There are few ways to make
   a new cubic ion:
 
         Input Arguments                   Resulting System
      ---------------------    -------------------------------------------
               -               Empty Magnetic Ion, No J (J=0)
            symbol             Magnetic Ion With Specified Ion (e.g. Ce3+)

    The destructor here does nothing except delete the ion index because the
    ions list, although using some small amount of memory, must remain until
    all Magnetic Ions in a  program are deleted.  This will not be done
    until program completion.                                                */
 

CubicIon::CubicIon()
  {
  if(!CubicIons.size()) AddCubicIons();	// Insure cubic ions are in lists
  CubicIonData CID("Ce3+");		// What we look for in ions list
  IsotopeData  ID("Ce3+");		// What we look for in isotopes list
  Iso = Isotope::seek(ID);		// Look for Ce3+ in isotopes list
  Ion = seek(CID);			// Look for Ce3+ in cubic ions list
  if(Ion < 0) CIFatal(2, "Ce3+"); 	// If no Ce3+, FATAL ERROR!
  if(Iso < 0) CIFatal(2, "Ce3+"); 	// If no Ce3+, FATAL ERROR!
  }					// Note We Use Ce3+ As Default Ion

CubicIon::CubicIon(const CubicIon& CI)
  {
  Iso = CI.Iso;				// Copy index in isotope list
  Ion = CI.Ion;				// Copy index in cubic ions list
  }

CubicIon::CubicIon(const std::string& CI)
  {
  if(!CubicIons.size()) AddCubicIons();	// Insure cubic ions are in list
  Iso = Isotope::seek(CI);		// Get index into isotopes   list
  Ion = seek(CI);			// Get index into cubic ions list
  if(Ion < 0) CIFatal(2, CI);	 	// If no CI, FATAL ERROR!
  if(Iso < 0) CIFatal(2, CI); 		// If no CI, FATAL ERROR!
  }
  
void CubicIon::operator= (const CubicIon& CI)
  { 
  Iso = CI.Iso;
  Ion = CI.Ion;
  }

CubicIon::~CubicIon() {};

 
// ____________________________________________________________________________
// B                        CUBIC ION ACCESS FUNCTIONS
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

	Input		CuSys   : CubicIon (this)
	Note			: All functions use CubicIons (ion data)     */
	   

//       double       CubicIon::qn()       const;               // INHERITED
//       int          CubicIon::HS()       const;               // INHERITED
//       std::string  CubicIon::momentum() const;               // INHERITED
// const std::string& CubicIon::symbol()   const;               // INHERITED
// const std::string& CubicIon::name()     const;               // INHERITED
// const std::string& CubicIon::element()  const;               // INHERITED
//       int          CubicIon::number()   const;               // INHERITED
//       int          CubicIon::mass()     const;               // INHERITED
//       double       CubicIon::weight()   const;               // INHERITED

int    CubicIon::charge()   const { return (CubicIons[Ion]).charge(); }
double CubicIon::gJ()       const { return (CubicIons[Ion]).gJ();     }
double CubicIon::beta()     const { return (CubicIons[Ion]).beta();   }
double CubicIon::gamma()    const { return (CubicIons[Ion]).gamma();  }
      
double  CubicIon::F4() const
  {
  switch(HS())
    {
    case  5: return 12; break;
    case  6: return 60; break;			// Ce3+, Sm3+
    case  7: return 15; break;
    case  8:					// Yb3+
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

double  CubicIon::F6() const
  {
  switch(HS())
    {
    case  5:
    case  6: return    0; break;		// Ce3+, Sm3+
    case  7: return  180; break;
    case  8:  					// Yb3+
    case  9: return 1260; break;		// Pr3+, Pm3+
    case 10: return 2520; break;		// Nd3+
    case 11: return 1260; break;
    case 12: return 3780; break;
    case 13:					// Tb3+, Tm3+
    case 14: return 7560; break;
    case 15: return 3780; break;
    case 16:					// Er3+, Dy3+
    case 17: return 13860; break;		// Ho3+
    default: return 0; 
    }
  }

// ____________________________________________________________________________
// C                        CUBIC ION I/O FUNCTIONS
// ____________________________________________________________________________
 
/* These functions allow users to write cubic ions in formatted ASCII to an
   output stream.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
     print     ostream     Writes cubic ion in ASCII to output stream
      <<       ostream     Writes cubic ion in ASCII to output stream     */

std::ostream& CubicIon::print(std::ostream& ostr) const
  {
  (CubicIons[Ion]).print(ostr, 0);
  ostr << "\n";
  (Isotopes[Iso]).print(ostr, 0,0);
  return ostr; 
  }
 
std::ostream& operator<< (std::ostream& ostr, const CubicIon& I) { return I.print(ostr); }

// ____________________________________________________________________________
// D                           CUBIC ION LIST FUNCTIONS
// ____________________________________________________________________________

/* This class contains a static list of possible cubic ions that can be the
   basis of a cubic ion. The functions herein allow for manipulation of the
   cubic ion list.

     Function        Output                       Description
   ------------  -------------  ----------------------------------------
       seek        ion index    Index in list of cubic ion input
      exists       ion symbol   True if ion specified is present in list
       known       ion symbol   True if ion specified is present in list
       
   The seek function will return -1 if the ion is not found in the list.     */
  
int CubicIon::seek(const CubicIonData& I)
  {
  int NIons = CubicIons.size();			// Number of ions stored
  if(!NIons) return -1;				// If no ions list we fail
  for(int i=0; i<NIons; i++)
    if(I.symbol() == (CubicIons[i]).symbol()) return i;
  return -1;
  }

bool CubicIon::exists(const std::string& symbol)
  {
  CubicIonData I(symbol);		// Get cubic ion data for symbol
  if(!CubicIons.size()) AddCubicIons();	// If CubicIons is empty, must
  if(seek(I)>=0) return true;
  return false;
  }

bool CubicIon::known(const std::string& symbol)
  { CubicIon X; return X.exists(symbol); }

void CubicIon::initialize()
  { CubicIon X; if(!CubicIons.size()) X.AddCubicIons(); }

int CubicIon::size() { return int(CubicIons.size()); }


void CubicIon::PrintList(std::ostream& ostr, bool hdr)
  {
  CubicIon X("Yb3+");				// Temporary cubic ion
  int Nions = X.size();				// Get size of ion list
  if(hdr)
    {
    std::string H("Currently Known Cubic Ions");
    ostr << CenterString(H);
    H = "(" + Gdec(X.size()) + " Available)";
    ostr << CenterString(H) << "\n";
    }

//                     Set Up Output Column Alignment

  int nl = Gdec(Nions).length() + 1;	// Length of printed index
  std::string nspc(nl+1, ' ');		// Spacer for index

  std::string sy;			// String for symbol
  int syl, syi, syf;			// For symbol alignment

  std::string csp("  ");		// Space between each columns
  int cspl = csp.length();		// Length of space

  std::string icol("     ");		// Space between isotopes

  int nml = 0;				// For Name column width
  int nmi, nmf;				// For Name column spacing
  int i;				// Dummy index
  CubicIonData CD;			// Temporary Cubic ion data
  CubicIon CI;				// Temporaryi Cubic Ion
  std::string symb, cname;

  for(i=0; i<Nions; i++)			// Find longest name
    {
    CD = CubicIons[i];
    symb = CD.symbol();
    CI = CubicIon(symb);
    cname = CI.name();
    nml = gmax(nml, int(cname.length()));
    }

  std::string nms("Name");		// For Name column header
  nmi = (nml-4)/2;
  nmf = nml - nmi - 4;
  if(nmi > 0) nms = std::string(nmi, ' ') + nms;
  if(nmf > 0) nms += std::string(nmf, ' ');

  int ncols = 2;			// Number of columns
  int next  = 0;			// For line end
  int cw  = nl + 1   + cspl;
  cw += 6            + cspl;
  cw += nms.length() + cspl;
  cw *= ncols;
  cw += (ncols-1)*icol.length();
  std::string cen(40-cw/2, ' ');

//			Write Column Headers

  ostr << "\n" << cen;
  for(i=0; i<ncols; i++)
    {
    ostr << "Indx"   << csp;		// Index
    ostr << "Symbol" << csp;		// Symbol
    ostr << nms      << csp;		// Name
    ostr << icol;			// Next Ion
    }

  ostr << "\n" << cen;
  for(i=0; i<ncols; i++)
    {
    ostr << std::string(nl+1, '=')         << std::string(cspl, ' ');
    ostr << std::string(6, '=')            << std::string(cspl, ' ');
    ostr << std::string(nms.length(), '=') << std::string(cspl, ' ');
    ostr << icol;
    }

  for(i=0; i<Nions; i++)
    {
    CD = CubicIons[i];
    sy = CD.symbol();
    CI = CubicIon(sy);

    syl = sy.length();
    syi = (6-syl)/2;
    syf = 6-syi-syl;
    if(syi>0) sy = std::string(syi, ' ') + sy;
    if(syf>0) sy += std::string(syf, ' ');
  
    nms = CI.name();
    nms += std::string(nml - nms.length(), ' ');
    if(!next) ostr << "\n" << cen;
    ostr << Gdec(i+1, nl) << "." << csp;		// Print index
    ostr << sy                   << csp;		// Print symbol
    ostr << nms                  << csp;
    next++;
    if(next >= ncols) next = 0;
    else              ostr << icol;
    }
  }

// ____________________________________________________________________________
// E                    CubicIon Container Support Functions
// ____________________________________________________________________________

/* Aside for providing basic tests as to whether two systems are equvalent or
   not, these operators are necessary if any STL container classes are to
   be used based on cubic ions (e.g. list<CubicIon> or vector<CubicIon>)  */
   
bool CubicIon::operator== (const CubicIon& CS) const { return (Ion==CS.Ion);  }
bool CubicIon::operator!= (const CubicIon& CS) const { return (Ion!=CS.Ion);  }
bool CubicIon::operator<(const   CubicIon& CS) const { return (HS()<CS.HS()); }
bool CubicIon::operator>(const   CubicIon& CS) const { return (HS()>CS.HS()); }

#endif								// CubicIon.cc
