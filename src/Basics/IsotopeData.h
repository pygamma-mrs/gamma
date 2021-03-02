/* IsotopeData.h ************************************************-*-c++-*-
**									**
** 	                          G A M M A				**
**									**
**	IsotopeData                                    Interface	**
**								 	**
**	Copyright (c) 1990, 1999			 		**
**	Tilo Levante, Scott A. Smith				 	**
**	Eidgenoessische Technische Hochschule			 	**
**	Labor fur physikalische Chemie				 	**
**	8092 Zurich / Switzerland				 	**
**								 	**
**      $Header: $
**								 	**
*************************************************************************/
     
/*************************************************************************
**							 		**
**  Description						 		**
**							 		**
**  Class IsotopeData embodies a single spin isotope.  Each object of	**
**  this type contains all the required data corresponding to a single	**
**  Isotope for use in magnetic resonance computation.  Some access	**
**  functions & output constants are provided.				**
**									**
**  Note that users will rarely (if ever) deal with this file.  GAMMA	**
**  contains a linked list of most known spin isotopes, each element of	**
**  the list being an object of type IsotopeData.  Higher classes and	**
**  GAMMA based programs use the linked list and class Isotope (a 	**
**  pointer into the isotopes list) rather than objects of this type.	**
**									**
*************************************************************************/
     
#ifndef   _IsotopeData_h_		// Is this file already included?
#  define _IsotopeData_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <string>			// Include libstdc++ strings
#include <vector>			// Include libstdc++ STL vectors
#include <GamGen.h>			// Know MSVCDLL (__declspec)

class IsotopeData
  {
  friend class Isotope;		// Allow Isotopes full access here
  int         _HS;		// Spin 2*I+1 values (2 if I=1/2)
  std::string _symbol;		// Isotope symbol (1H, 2H, 14N, ...)
  std::string _name;		// Isotope name (Hydrogen, Carbon, ..)
  std::string _element;		// Element type (H, He, ...)
  int         _number;		// Periodic number (1 for H, 3 for Li)
  int         _mass;		// Atomic mass in amu (protons + neutrons)
  double      _weight;		// Atomic weight (in grams/mole)
  double      _receptivity;	// Receptivity
  double      _relfreq;		// Relative frequency
  bool        _iselectron;	// Flag if electron (true) or nucleus (false)

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// A                    SINGLE ISOTOPE CONSTRUCTORS
// ____________________________________________________________________________


MSVCDLC IsotopeData();
MSVCDLC IsotopeData(const IsotopeData& data);
  
MSVCDLC IsotopeData(int HS_, const std::string& symb_, const std::string& name_, 
                        std::string element_, int number_, int mass_, double weight_,
                         double recept_, double rel_freq_, bool is_electron_=false);
 
        // Input              HS_ : Int for Isotope Hilbert space (2, 3, 4)
        //                  symb_ : String for Isotope symbol     (1H, 2H, 14N)
        //                  name_ : String for Isotope name       (Hydrogen)
        //               element_ : String for Isotope element    (H, Li, C)
        //                number_ : Int for the Isotope number	  (1@H, 3@Li)
        //                  mass_ : Int for the Isotope mass	  (1@H, 13@C)
        //                weight_ : Double Isotope weight	  (in g/mol)
        //                recept_ : Double Isotope receptivity	  ()
        //              rel_freq_ : Double Isotope rel. freq.	  (in MHz)
        // Output            this : Constructs IsotopeData from
        //                          all information which each isotope contains
 	// Note			  : There is NO check to insure that
	//			    the input parameters are consistent
	//			    Allows users to "make" pseudo isotopes!


MSVCDLC IsotopeData(const std::string& symbol_);
 
        // Input          symbol_ : A string of an Isotope symbol
        // Output            this : Constructs IsotopeData for the
        //                          isotope type described by symbol_
 	// Note			  : This ONLY sets the symbol, nothing
 	//			    else and will produce problems if the
 	//			    other data parameters aren't filled.

 
MSVCDLL IsotopeData& operator= (const IsotopeData& ID1);

MSVCDLC      ~IsotopeData();

// ____________________________________________________________________________ 
// B                    SINGLE ISOTOPE ACCESS FUNCTIONS
// ____________________________________________________________________________

/*
        Function  Return    Value                 Example(s)
        --------  -------   ------   ------------------------------------------
           qn     double     mz      0.5, 1.0, 1.5, ...... in units of hbar
           HS     integer   2*mz+1   2<-1/2, 3<-1,......  spin Hilbert space
        momentum  string     mz      1/2, 1, 3/2, ......  in units of hbar
         symbol   string             1H, 2H, 19F, 13C, ......
          name    string             Hydrogen, Lithium, Carbon, ....
         element  string             H, Li, F, C, ....
	 number   integer   at.#     1<-H, 3<-Li, 6<-C, ....
	  mass    integer   at.mass  1<-1H, 2<-2H, 13<-13C, .... in amu units
	 weight   double    at.wt.   1.00866<-1H, 7.016<-7Li, ... in grams/mole
         recept.  double             5680, 1540, .....
         rel.frq. double             400.13, 155.503, .. in MHz (1H based)   */

       MSVCDLL double       qn()	    const;		// Calculated
       MSVCDLL int          HS()	    const;		// Stored
       MSVCDLL std::string  momentum() const;		// Calculated
const  MSVCDLL std::string& symbol()   const;		// Stored
const  MSVCDLL std::string& name()	    const;		// Stored
const  MSVCDLL std::string& element()  const;		// Stored
       MSVCDLL int          number()   const;		// Stored
       MSVCDLL int          mass()	    const;		// Stored
       MSVCDLL double       weight()   const;		// Stored
       MSVCDLL double       recept()   const;		// Stored
       MSVCDLL bool         electron() const;		// Stored
       MSVCDLL double       rel_freq() const;		// Stored
 
// ____________________________________________________________________________
// C                    SINGLE ISOTOPE I/O FUNCTIONS
// ____________________________________________________________________________
 
        // Input             ostr : An output stream ostr
        //                     ID : An IsotopeData (this)
	//		       lf : Line feed flag
        // Output            ostr : The output stream, modified to 
        //                          contain the information about ID
 
MSVCDLL std::vector<std::string> printStrings(bool hdr=true) const;
MSVCDLL std::ostream& print(std::ostream& ostr, int lf=1, bool hdr=true) const;
MSVCDLL friend  std::ostream& operator<< (std::ostream& ostr, const IsotopeData& ID);


};						// End IsotopeData Class


#endif							// IsotopeData.h

