/* sys_gradz.i */

%{
#include "Gradients/sys_gradz.h"
%}


%include "std_vector.i"

%feature("autodoc", "1" );

class sys_gradz: public spin_system
{
int                 _NSS;	
double              dBodm;	
double              efflen;	
std::vector<double> Ps;	

void SysGZerr (int  eidx,                        int nr=0) const;
void SysGZerr(int   eidx, const std::string& st, int nr=0) const;
volatile void SysGZfatal(int eidx) const;

virtual int setSsys(const ParameterSet& pset,int idx=-1,int wrn=2);

public:

sys_gradz(int spins=0);
sys_gradz(const sys_gradz& sys);
~sys_gradz();
sys_gradz& operator= (const sys_gradz& sys);


void NSS(int nss);
int  NSS() const;

void   BoGrad(double bgrad);
double BoGrad()             const;
double GradVal(double dist) const;

void   SysLen(double len);
double SysLen()         const;
double SysDist(int nss) const;

spin_system SubSys(int nss)                    const;
double      SubSysShift(int nss, int spin)     const;
double      SubSysShift(double dist, int spin) const;
double      SubSysPPM(int nss, int spin)       const;
double      SubSysPPM(double dist, int spin)   const;

//operator ParameterSet( ) const;

void PSetAdd(ParameterSet& pset, int idx=-1) const;

void setSubSys(const ParameterSet& pset);

void setBoGrad(const ParameterSet& pset);

void setLength(const ParameterSet& pset);

void operator= (const ParameterSet& pset);

virtual int write(const std::string &filename, int idx=-1, int warn=2) const;
virtual int write(std::ofstream& ofstr,        int idx=-1, int warn=2) const;

virtual int read(const std::string& filename, int idx=-1, int warn=2);
virtual int read(const ParameterSet& pset, int idx=-1, int warn=2);

std::string ask_read(int argc, char* argv[], int argn);
std::string ask_read(int argc, char* argv[], int argn,
                                                     const std::string& def);
};

