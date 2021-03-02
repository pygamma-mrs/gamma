/* IsotopeData.i */
// Swig interface file

%{
#include "Basics/IsotopeData.h"
%}

%include "std_string.i"
%include "std_vector.i"

%feature("autodoc", "1" );


%rename(__assign__) IsotopeData::operator=;


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


public:

IsotopeData();
IsotopeData(const IsotopeData& data);
  
IsotopeData(int HS_, const std::string& symb_, const std::string& name_, 
                        std::string element_, int number_, int mass_, double weight_,
                         double recept_, double rel_freq_, bool is_electron_=false);

IsotopeData(const std::string& symbol_);
 
IsotopeData& operator= (const IsotopeData& ID1);

~IsotopeData();

       double       qn()	    const;		// Calculated
       int          HS()	    const;		// Stored
       std::string  momentum() const;		// Calculated
const  std::string& symbol()   const;		// Stored
const  std::string& name()	    const;		// Stored
const  std::string& element()  const;		// Stored
       int          number()   const;		// Stored
       int          mass()	    const;		// Stored
       double       weight()   const;		// Stored
       double       recept()   const;		// Stored
       bool         electron() const;		// Stored
       double       rel_freq() const;		// Stored

std::vector<std::string> printStrings(bool hdr=true) const;

//std::ostream& print(std::ostream& ostr, int lf=1, bool hdr=true) const;
//friend  std::ostream& operator<< (std::ostream& ostr, const IsotopeData& ID);


};						// End IsotopeData Class


%extend IsotopeData {
%pythoncode %{

def __str__(self):
    """Prints out isotope data"""

    sss = ""
    sss += "qn: "          + str(self.qn()) + "\n"
    sss += "HS: "          + str(self.HS()) + "\n"
    sss += "momentum: "    + str(self.momentum()) + "\n"
    sss += "symbol: "      + str(self.symbol()) + "\n"
    sss += "name: "        + str(self.name()) + "\n"
    sss += "element: "     + str(self.element()) + "\n"
    sss += "number: "      + str(self.number()) + "\n"
    sss += "mass: "        + str(self.mass()) + "\n"
    sss += "weight: "      + str(self.weight()) + "\n"
    sss += "receptivity: " + str(self.recept()) + "\n"
    sss += "electron: "    + str(self.electron()) + "\n"
    sss += "rel freq: "    + str(self.rel_freq()) + "\n"

    return (sss)


%}
};

