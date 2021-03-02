// SpinSystem.i
// Swig interface file.

%{
#include "HSLib/SpinSystem.h"
%}


%include "std_string.i"
%include "std_vector.i"
%include "HSLib/SpinSys.i"

%feature("autodoc", "1" );
//%feature("director" );

%rename(__assign__) spin_system::operator=;


class spin_system: public spin_sys
{

public:

spin_system(int spins=0);
spin_system(const spin_system &sys);

virtual       ~spin_system ();
spin_system&  operator= (const spin_system &sys);

virtual void   shifts(double shift=0);
virtual void   shift(int, double);
virtual double shift(int)          const;

double maxShift()                  const;
double maxShift(const std::string& Iso) const;
double minShift()                  const;
double minShift(const std::string& Iso) const;
double medianShift()               const;
double lab_shift(int)              const;	// Typically ~10^8 !

virtual void   offsetShifts(double OF,int i=0); 
virtual void   offsetShifts(double OF,const std::string& Iso); 
virtual void   PPM (int, double);

double PPM (int) const;


double gfactor(int    spin) const;
void   gfactor(int    spin, double g);
double eshift(int     spin) const;
double lab_eshift(int spin) const;
double efield(int     spin) const;
double efield_lab(int spin) const;


virtual void   Js(double Jval=0);
virtual void   J(int, int, double);
virtual void   J(double, int, int);
        double J(int, int) const;


virtual void   As(double Aval=0);
virtual void   A(int, int, double);
virtual void   A(double, int, int);
        double A(int, int)   const;
        double AHz(int, int) const;

void   Omega(double Om);				// Om : 1H spect. freq.  (MHz)
void   Omega(double Om, const std::string& iso);	// Om : iso spect. freq. (MHz)
double Omega(int spin=-1) const;			// Return Larmor of spin (MHz)
double Omega(const std::string& iso) const;	// Return Larmor of iso  (MHz)
double Bo() const;					// Return Field (Gauss)
void   OmegaAdjust(double Om);			// Reset 1H Larmor, PPM static
void   FieldAdjust(double B);			// Reset Field Strength

        
void   spectrometer_frequency(double freq);
double spectrometer_frequency() const;


void spflags(int TF);

void spflag(int spin1, int spin2, int TF);
 
int spflag(int spin1, int spin2) const;

double center(int spin=0);

double Nyquist(int spin,               double fact, double lwhh) const;
double Nyquist(const std::string& iso, double fact, double lwhh) const;
double Nyquist(const Isotope& iso,     double fact, double lwhh) const;

 
//operator ParameterSet( ) const;


//friend void operator+= (ParameterSet& pset, const spin_system &ss);


virtual void PSetAdd(ParameterSet& pset, int idx=-1) const;

void setJs(const ParameterSet& pset);

void setAs(const ParameterSet& pset);

void setShifts(const ParameterSet& pset);

void setGs(const ParameterSet& pset);

//virtual void operator= (const ParameterSet& pset);

virtual int write(const std::string &filename, int idx=-1, int warn=2) const;
//virtual int write(std::ofstream& ofstr,        int idx=-1, int warn=2) const; 

virtual int read(const std::string& fn,    int idx=-1, int warn=2);
virtual int read(const ParameterSet& pset, int idx=-1, int warn=2);

virtual std::string ask_read(int argc, char* argv[], int argn);
virtual std::string ask_read(int argc, char* argv[], int argn,
                                                    const std::string& def);

//std::ostream& printvs(std::ostream& ostr) const;
//std::ostream& printGs(std::ostream& ostr) const;
//std::ostream& printJs(std::ostream& ostr) const;
//std::ostream& printAs(std::ostream& ostr) const;
//std::ostream& printO(std::ostream&  ostr) const;


//virtual std::ostream& print(std::ostream& out, bool hdr=true) const;
//friend  std::ostream& operator<<(std::ostream& out, const spin_system& sys);


virtual std::vector<std::string> SYSStrings(int w1=10,int w2=12,int w3=1) const;
std::vector<std::string> VStrings(int   colwd=12, int digs=2) const;
std::vector<std::string> PPMStrings(int colwd=12, int digs=2) const;
std::vector<std::string> GFStrings(int  colwd=12, int digs=2) const;
std::vector<std::string> BeStrings(int  colwd=12, int digs=2) const;
std::vector<std::string> JStrings(int   colwd=12, int digs=2) const;
std::vector<std::string> AStrings(int   colwd=12, int digs=2) const;
std::vector<std::string> OmStrings(int  colwd=12, int digs=2) const;

/*
protected:
bool setOm(const ParameterSet& pset);
bool checkNotE(int spin, int warn=1) const;
bool checkNotN(int spin, int warn=1) const;
*/

};


%extend spin_system {
%pythoncode %{

def __str__(self):
    """Prints out spin system"""

    sss = ""
    for v in self.SYSStrings():
        sss += str(v) + '\n'

    return (sss)


def __repr__(self):
    return self.__str__()    


%}
};
