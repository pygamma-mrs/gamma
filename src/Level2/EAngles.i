/* EAngles.i */

%{
#include "Level2/EAngles.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

%rename(__assign__) EAngles::operator=;

class EAngles
{
double _alpha;
double _beta;
double _gamma;
static double AngCutoff;

void EAerror(int eidx,                        int n=0) const;
void EAerror(int eidx, const std::string& pn, int n=0) const;
volatile void EAfatal(int eidx=0)                               const;

void SetAngles(double alpha, double beta, double gamma, bool d=false);
void SetAlpha(double  alpha, bool deg=false);
void SetBeta(double   beta,  bool deg=false);
void SetGamma(double  gamma, bool deg=false);

bool SetEAngles(const ParameterSet& pset, int idx=-1, bool warn=true);
bool SetEASet(const   ParameterSet& pset, int idx=-1, bool warn=true);
bool Set3Angles(const ParameterSet& pset, int idx=-1, bool warn=true);

public:

EAngles();
EAngles(double alpha, double beta=0, double gamma=0, bool deg=false);
EAngles(const coord&   EA, bool deg=true);
EAngles(const EAngles& EA);
~EAngles();

EAngles& operator= (const EAngles& EA);

double alpha() const;
double beta( ) const;
double gamma() const;

void alpha(double A);
void beta(double  B);
void gamma(double G);

EAngles operator*  (const EAngles& EA) const;
EAngles&    operator*= (const EAngles& EA);
EAngles&    operator&= (const EAngles& EA);
EAngles composite  (const EAngles& EA) const;

SinglePar param(const std::string& pn,                      bool deg=true) const;
SinglePar param(const std::string& pn,const std::string& ps,bool deg=true) const;

//operator ParameterSet( ) const;
void PSetAdd(ParameterSet& pset,        int idx=-1, bool deg=true) const;
void write(const std::string &filename, int idx=-1, bool deg=true) const;

bool read(const std::string &filename, int idx=-1, int warn=2);
bool read(const ParameterSet& pset,    int idx=-1, int warn=2);

static void SetCutoff(double co=-1);
bool operator== (const EAngles& EA) const;
bool operator!= (const EAngles& EA) const;
bool operator<  (const EAngles& EA) const;
bool operator>  (const EAngles& EA) const;

bool    equal(const EAngles& EA, double CUTOFF=1.e-10) const;
EAngles inverse()                                      const;
matrix  RMx(bool inv=false)                            const;

matrix Rmx()    const;
matrix invRmx() const;

};
