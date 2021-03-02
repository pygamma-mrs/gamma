/* Quaternion.i */

%{
#include "Level2/Quaternion.h"
%}

%feature("autodoc", "1" );

%rename(__assign__) quatern::operator=;

class quatern
{
double        AQ, BQ, CQ, DQ;			
static int    QRange;				
static bool   SinPos, CosPos, TanPos, SCTF;	
static double ElemCutoff;                    
typedef quatern quartern;			

void Qerror(int eidx,                          int noret=0) const;
void Qerror(int eidx,   const std::string& pn, int noret=0) const;
volatile void Qfatal(int eidx=0)                                     const;
volatile void Qfatal(int eidx,   const std::string& pn)              const;

bool CheckNorm(int warn=2) const;
bool CheckNorm(double A,double B,double C,double D,int warn=2) const;
bool CheckNorm(const matrix& Rotmx, bool warn=true) const; 

bool SetQuatern(const ParameterSet& pset, int idx=-1, int warn=2);

bool SetSinPos() const;
bool SetCosPos() const;
bool SetTanPos() const;

double  FindBeta()                             const;
double  FindAlpha()                            const;
double  FindAlpha(double beta)                 const;
double  FindGamma()                            const;
double  FindGamma(double dbeta, double dalpha) const;
EAngles FindEAs()                              const;
double  GetAngle(double sinval, double cosval) const;

public:

quatern();
quatern(const coord&   ABG, bool inv=false);
quatern(const EAngles& EA,  bool inv=false);
quatern(const quatern& Qrt, bool inv=false);
quatern(const ParameterSet& pset, int idx=-1, int warn=2);
quatern(double QA, double QB, double QC, double QD, bool inv=false);


~quatern();
quatern& operator= (const quatern& QRT);
quatern& operator= (const coord&   EA);
quatern& operator= (const EAngles& EA);

double A() const;
double B() const;
double C() const;
double D() const;

double  alpha() const;
double  beta()  const;
double  gamma() const;
EAngles EA()    const;
coord   ABG()   const;


quatern operator*  (const quatern& Q) const;
quatern&    operator*= (const quatern& Q);
quatern&    operator&= (const quatern& Q);

quatern composite(const quatern& Q, bool rev=false) const;
matrix RotMx() const; 						
matrix RMx() const;

double  norm()    const;
quatern inverse() const;

bool operator== (const quatern& Quar) const;
bool operator!= (const quatern& Quar) const;
bool operator<  (const quatern& Quar) const;
bool operator>  (const quatern& Quar) const;

SinglePar param()    const;
SinglePar param(const std::string& pname)  const;
SinglePar param(const std::string& pname, const std::string& pstate) const;

//operator ParameterSet( ) const;
bool PSetAdd(ParameterSet& pset, int idx=-1, int pfx=-1) const;

bool write(const std::string& fo,int idx=-1,int pfx=-1,int warn=2) const;
bool write(    std::ofstream& of,int idx=-1,int pfx=-1,int warn=2) const;

bool read(const std::string&  filein, int indx=-1, int warn=2);
bool read(const ParameterSet& pset,   int indx=-1, int warn=2);

static bool ASinPos();
static bool ACosPos();
static bool ATanPos();

void ShowConversion() const;

static bool ValidRMx(const matrix& R, bool msgs=true);

};

extern  quatern composite(const EAngles&,  const EAngles&);
extern  quatern composite(const coord&,  const coord&);
extern  quatern composite(const quatern&, const quatern&);
extern  quatern composite(const matrix&,  const quatern&);
