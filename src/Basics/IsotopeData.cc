/* IsotopeData.cc ***********************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      IsotopeData                                 Implementation	**
**                                                                      **
**      Copyright (c) 1990, 1999                                        **
**      Tilo Levante, Scott A. Smith                                    **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fur physikalische Chemie                                  **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/
     
/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  Class IsotopeData embodies a single spin isotope.  Each object of   **
**  this type contains all the required data corresponding to a single  **
**  Isotope for use in magnetic resonance computation.  Some access     **
**  functions & output constants are provided.                          **
**                                                                      **
**  Note that users will rarely (if ever) deal with this file.  GAMMA   **
**  contains a linked list of most known spin isotopes, each element of **
**  the list being an object of type IsotopeData.  Higher classes and   **
**  GAMMA based programs use the linked list and class Isotope (a       **
**  pointer into the isotopes list) rather than objects of this type.   **
**                                                                      **
*************************************************************************/

#ifndef   IsotopeData_cc_		// Is this file already included?
#  define IsotopeData_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <Basics/IsotopeData.h>		// Include the interface
#include <Basics/StringCut.h>		// Include string cutting
#include <iostream>			// Include libstdc++ iostreams
#include <vector>			// Include llibstdc++ STL vectors

using std::string;			// Using libstdc++ strings

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                      SINGLE ISOTOPE CONSTRUCTORS
// ____________________________________________________________________________

/* In these constructors the individual values are best set in the
   same order as they exist in the class definition.                         */
 
        // Input              HS_ : Int for Isotope Hilbert space (2, 3, 4)
        //                  symb_ : String for Isotope symbol     (1H, 2H, 14N)
        //                  name_ : String for Isotope name       (Hydrogen)
        //                  elem_ : String for Isotope element    (H, Li, C)
        //                  numb_ : Int for the Isotope number	  (1@H, 3@Li)
        //                  mass_ : Int for the Isotope mass	  (1@H, 13@C)
        //                  wght_ : Double Isotope weight	  (in g/mol)
        //                  rcpt_ : Double Isotope receptivity	  ()
        //                  rfrq_ : Double Isotope rel. freq.	  ()
	//		    elec_ : bool flag for electron/nucleus
        // Output            this : Constructs IsotopeData from
        //                          all information each isotope contains
 	// Note			  : Construction with string for symbol will
	//			    ONLY sets the isotope symbol, nothing
 	//			    else and will produce problems if the
 	//			    other data parameters aren't filled.


IsotopeData::IsotopeData()
  {
  _HS=0;
  _number=0;
  _mass=0;
  _weight=0;
  _receptivity=0;
  _relfreq=0;
  _iselectron=false;
  }

IsotopeData::IsotopeData(const IsotopeData& ID)
            : _HS(ID._HS),
              _symbol(ID._symbol),
              _name(ID._name),
              _element(ID._element),
              _number(ID._number),
              _mass(ID._mass),
              _weight(ID._weight),
              _receptivity(ID._receptivity),
              _relfreq(ID._relfreq),
              _iselectron(ID._iselectron) {}

IsotopeData::IsotopeData(int HS_, const string& symb_, const string& name_,
			 string elem_, int numb_,    int mass_, double wght_,
                         double rcpt_, double rfrq_, bool elec_)
            : _HS(HS_),
              _symbol(symb_),
              _name(name_),
              _element(elem_),
              _number(numb_),
              _mass(mass_),
              _weight(wght_),
              _receptivity(rcpt_),
              _relfreq(rfrq_),
              _iselectron(elec_) {}

IsotopeData::IsotopeData(const string& symb) : _symbol(symb)
  { 
  _HS          = 0;
  _name        = string("Unknown");
  _element     = string("XXX");
  _number      = 0;
  _mass        = 0;
  _weight      = 0;
  _receptivity = 0;
  _relfreq     = 0;
  _iselectron  = false;
  }

IsotopeData& IsotopeData::operator= (const IsotopeData& ID1)
  {
  _HS          = ID1._HS;
  _symbol      = ID1._symbol;
  _name        = ID1._name;
  _element     = ID1._element;
  _number      = ID1._number;
  _mass        = ID1._mass;
  _weight      = ID1._weight;
  _receptivity = ID1._receptivity;
  _relfreq     = ID1._relfreq;
  _iselectron  = ID1._iselectron;

  return (*this);
  }

IsotopeData::~IsotopeData() {}
 
// ____________________________________________________________________________
// B                    SINGLE ISOTOPE ACCESS FUNCTIONS
// ____________________________________________________________________________
 
/*      Function  Return    Value                 Example(s)
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

      double       IsotopeData::qn()       const { return (_HS ? (_HS-1)/2.0 : 0); }
      int          IsotopeData::HS()       const { return _HS; }
const string& IsotopeData::symbol()   const { return _symbol; }
const string& IsotopeData::name()     const { return _name; }
const string& IsotopeData::element()  const { return _element; }
      int          IsotopeData::number()   const { return _number; }
      int          IsotopeData::mass()     const { return _mass; }
      double       IsotopeData::weight()   const { return _weight; }
      bool         IsotopeData::electron() const { return _iselectron; }
      double       IsotopeData::recept()   const { return _receptivity; }
      double       IsotopeData::rel_freq() const { return _relfreq; }

string IsotopeData::momentum() const
  {
  if(!_HS) return string("0"); 		// If no _HS, then it's zero
  else if(_HS%2) return Gdec(int((_HS-1)/2));	// If odd _HS, return fraction
  else return string(Gdec(int(_HS-1)))+ string("/2");	// If even _HS, return whole int
  }

// ____________________________________________________________________________
// C                    SINGLE ISOTOPE I/O FUNCTIONS
// ____________________________________________________________________________

/* These send basic information about the isotope to the output stream.*/

        // Input               ID : An IsotopeData (this)
        //                   ostr : An output stream ostr
	//		       lf : Line feed flag
	//		      hdr : Flag if header output
        // Output            ostr : The output stream, modified to 
        //                          contain the information about ID

std::vector<string> IsotopeData::printStrings(bool hdr) const
  {
  std::vector<string> PStrings;
  if(hdr) 
    PStrings.push_back(string("Isotope       ") + _symbol);
  if(_name.length())
    PStrings.push_back(string(" Name         ") +  _name);
  else
    PStrings.push_back(string(" Name         Unknown"));
  if(_element.length())
    PStrings.push_back(string(" Element      ") +  _element);
  else
    PStrings.push_back(string(" Element      Unknown"));
    PStrings.push_back(string(" Number       ") +  Gdec(_number));
    PStrings.push_back(string(" Mass         ") +  Gdec(_mass) + string(" amu"));
    string Stmp(" Weight       ");
    string Sfrm("%8.4f");
    if(_weight < 100) Sfrm = string("%7.4f");
    if(_weight < 10)  Sfrm = string("%6.4f");
    PStrings.push_back( Stmp +  Gform(Sfrm, _weight) + string(" g/m"));
    PStrings.push_back(string(" Spin         ") +  Gform("%3.1f", qn()) + string(" (hbar)"));
    string X = " Type         ";
    if(_iselectron) X += "electron";
    else            X += "nucleus";
    PStrings.push_back(X);
  return PStrings;
  }

std::ostream& IsotopeData::print(std::ostream& ostr, int lf, bool hdr) const
  {
  std::vector<string> PStrings = printStrings(hdr);
  for(unsigned i=0; i<PStrings.size(); i++)
    ostr << PStrings[i] << std::endl;
  if(lf) ostr << std::endl;
  return ostr;
  }

std::ostream& operator<< (std::ostream& ostr, const IsotopeData& ID)
  { return ID.print(ostr); }


#endif							// IsotopeData.cc
