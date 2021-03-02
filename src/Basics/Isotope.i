/* Isotope.i */
// Swig interface file

%{
#include "Basics/Isotope.h"
%}

%include "std_string.i"
%include "std_vector.i"

%feature("autodoc", "1" );

%rename(__eq__)  Isotope::operator== const;
%rename(__ne__)  Isotope::operator!= const;
%rename(__lt__)  Isotope::operator<  const;
%rename(__gt__)  Isotope::operator>  const;

class Isotope
{

protected:	// Give derived classes access

  static std::vector<IsotopeData> Isotopes;	// Pointer to isotopes list
  static double                   RELFRQ1H;	// Relative 1H Larmor frequency
         int                      Iso;		// An isotope index

         void Isoerror(int eidx,                           int noret=0) const;
         void Isoerror(int eidx, const std::string& pname, int noret=0) const;
volatile void Isofatal(int eidx)                                        const;
volatile void Isofatal(int eidx, const std::string& pname)              const;
 
void set_Isotope_list();
 
void SetRel1HF();

bool SetIsotope(const ParameterSet& pset, int idx=-1, bool warn=true);



public:


Isotope();

Isotope(const    Isotope& I);
Isotope(const    std::string&  I);

// Isotope& operator= (const      Isotope& I);

virtual  ~Isotope();

       double       qn()                 const;
       int          HS()                 const;
       std::string  momentum()           const;
const  std::string& symbol()             const;
const  std::string& name()               const;
const  std::string& element()            const;
       int          number()             const;
       int          mass()               const;
       double       weight()             const;
       double       gamma()              const;
       double       receptivity()        const;
       double       relative_frequency() const;
       bool         electron()           const;

virtual  bool read(const std::string& filename, int idx=-1, int warn=2);
virtual  bool read(const ParameterSet& pset,    int idx=-1, int warn=2);


std::vector<std::string> printStrings(bool hdr=true) const;

//virtual  std::ostream& print(std::ostream& ostr) const;

//friend   std::ostream& operator<< (std::ostream&    ostr, const Isotope& I);

virtual  int                      seek(const   IsotopeData& ID);
virtual  bool                     exists(const std::string& symbol);
static   bool                     known(const  std::string& symbol);
static   int                      size();
static   std::vector<std::string> PrintListStrings();

//static   void                     PrintList(std::ostream& ostr, bool hdr=true);

static  bool AddIsotope(const IsotopeData& ID, int warn=2);

virtual  bool operator== (const Isotope& I) const;
virtual  bool operator!= (const Isotope& I) const;
virtual  bool operator<  (const Isotope& I) const;
virtual  bool operator>  (const Isotope& I) const;

 bool nepair(const Isotope& S) const;
 bool enpair(const Isotope& S) const;
 bool eepair(const Isotope& S) const;
 bool nnpair(const Isotope& S) const;

};

%extend Isotope{
%pythoncode %{

def __str__(self):
    """Prints out Isotope"""

    sss = ""
    for v in self.printStrings():
        sss += str(v) + '\n'

    return (sss)

%}
};

