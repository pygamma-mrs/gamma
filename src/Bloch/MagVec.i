/* MagVec.i */

// Swig interface file.

%{
#include "Bloch/MagVec.h"
%}

%include "std_string.i"
%include "std_vector.i"

%feature("autodoc", "1" );

#ifdef SWIGPYTHON
%rename(__add__)  MagVec::operator+ const;
%rename(__iadd__) MagVec::operator+=;
%rename(__sub__)  MagVec::operator- const;
%rename(__isub__) MagVec::operator-=;
#endif

class MagVec : public col_vector
{

friend class BlochSys;

void MVerror(int ei,                        int nr=0) const;
void MVerror(int ei, const std::string& pn, int nr=0) const;

volatile void MVfatal(int ei) const;
volatile void MVfatal(int ei, const std::string& pn) const;

bool CheckNorms(const  std::vector<double>&  Ns, bool warn=true) const;
bool CheckRange(int cmp,                    bool warn=true) const;

bool SetVector(const   ParameterSet& P, int pfx=-1, bool W=true);
bool GetNVects(const   ParameterSet& P, int&     N, bool W=true) const;
bool SetSubVects(const ParameterSet& P, int      N, bool W=true);
bool GetCoord(const  ParameterSet& P, 
                                     coord& pt, int idx=-1, bool W=true) const;
bool GetMxMyMz(const   ParameterSet& P, 
            double& Mx, double& My, double& Mz, int idx=-1, bool W=true) const;

///bool GetVect(const ParameterSet& pset, int i, double& v, Isotope& I,
///             double& R1, double& R2, coord& Pt, int& Sp, bool warn=true) const;


public:


MagVec(int   N=0);
MagVec(const MagVec& MV);
MagVec(const col_vector& CV);

MagVec(double Mx, double My, double Mz);
MagVec(const coord& M);
MagVec(double Mx1,double My1,double Mz1,double Mx2,double My2,double Mz2);
MagVec(const coord& M1, const coord& M2);
MagVec(const std::vector<coord>& Ms);

///        ~MagVec();
///MagVec& operator= (const MagVec& MV);
///MagVec operator-  ()                 const;

MagVec  operator+  (const MagVec& M1) const;
MagVec& operator+= (const MagVec& M1);
MagVec  operator-  (const MagVec& M1) const;
MagVec& operator-= (const MagVec& M1);

///int size()   const;			// Inherited
int NComps() const;

double Mx(int cmp) const;
double My(int cmp) const;
double Mz(int cmp) const;

void   Mx(int cmp, double mx);
void   My(int cmp, double my);
void   Mz(int cmp, double mz);

double x(int     cmp=0) const;
double y(int     cmp=0) const;
double z(int     cmp=0) const;
double norm(int  cmp=0) const;
double theta(int cmp=0) const;
double phi(int   cmp=0) const;

std::vector<double> Norms() const;
void                Norms(const std::vector<double>& Ns);

double Norm(int i) const;
void   Norm(double nv, int i);

// The next line gave an error in swig.
//operator ParameterSet( ) const;
//friend void operator+= (ParameterSet& pset, const MagVec& MV);
bool PSetAdd(ParameterSet& pset, int pfx=-1)   const;

bool write(const std::string& filename, int pfx=-1, int warn=2) const;
bool write(std::ofstream& ofstr,        int pfx=-1, int warn=2) const;

bool read(const std::string& fn,    int idx=-1, int warn=2);
bool read(const ParameterSet& pset, int idx=-1, int warn=2);

std::string ask_read(int argc, char* argv[], int argn);
std::string ask_read(int argc, char* argv[], int argn,
                                   const std::string& def); 

std::vector<std::string> printStrings() const;
//std::ostream&            print(std::ostream& out, int np=20) const;
//friend std::ostream& operator<<(std::ostream& out, const MagVec& M);

MagVec Mx()    const;
MagVec My()    const;
MagVec Mz()    const;

static MagVec MxVec(int NC);
static MagVec MyVec(int NC);
static MagVec MzVec(int NC);

MagVec MxVec() const;
MagVec MyVec() const;
MagVec MzVec() const;

};

