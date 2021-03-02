/* MultiSys.i */

%{
#include "MultiSys/MultiSys.h"
%}

%include "std_vector.i"
%include "std_string.i"

%feature("autodoc", "1" );

%rename(__assign__) multi_sys::operator=;

class multi_sys

{
  std::string              _SysName;	
  std::vector<double>      _Pops;	
  std::vector<sys_dynamic> _Comps;	
  std::vector<ExchProc>	   _Exs;	

private:

void MSYSerror(int eidx,                           int noret=0) const;
void MSYSerror(int eidx, const std::string& pname, int noret=0) const;
volatile void MSYSfatal(int eidx)                                        const;
volatile void MSYSfatal(int eidx, const std::string& pname)              const;

bool getNComps(const ParameterSet& pset,int& ncmps,bool warn=true);
bool getMSName(const ParameterSet& pset, 
                                            std::string& name,bool warn=true);

bool getFName(const  ParameterSet& pset,
                                   std::string& name, int idx, bool warn=true);

bool getFNames(const ParameterSet& pset,
                              std::vector<std::string>& names, bool warn=true);

bool getComp(const ParameterSet& pset, int idx,
                                             sys_dynamic& cmp, bool warn=true);

bool getComps(const ParameterSet& pset, int ncmps,
                               std::vector<sys_dynamic>& cmps, bool warn=true);

bool getPop(const ParameterSet& pset, int idx,
                                            double& pop, bool warn=true) const;

bool getPops(const ParameterSet& pset, int ncmps,
                              std::vector<double>& pops, bool warn=true) const;

int  getNex(const ParameterSet& pset) const;

bool getProcesses(const ParameterSet& pset, 
                            std::vector<ExchProc>& procs, bool warn=true) const;

bool getMsys(const ParameterSet& pset,
            std::string& name, std::vector<sys_dynamic>& cmps,
       std::vector<double>& pops, std::vector<ExchProc>& procs, bool warn=true);

bool setMsys(const ParameterSet& pset,              bool warn=true);

bool CheckNComps(int                nc, bool warn=true) const;
bool CheckRange(unsigned             n, bool warn=true) const;
bool CheckProc(int                  ip, bool warn=true) const;
bool CheckProcs(                        bool warn=true) const;
bool CheckProc(const ExchProc& XP,      bool warn=true) const;
bool CheckField(const spin_system& sys, bool warn=true) const;

public:

multi_sys();

multi_sys(const multi_sys& msys);

multi_sys(double pop1, sys_dynamic &sys1, 
                               double pop2, sys_dynamic &sys2, double krate=0);
multi_sys& operator= (const multi_sys& msys);

~multi_sys();
void               name(const std::string& sysname);
const std::string& name() const;

void   pop(int icomp, double npop);
double pop(int icomp) const;
double popmin() const;
double popmax() const;

int                NComps() const;
void               Comp(int icomp, const sys_dynamic& sys);
const sys_dynamic& Comp(int icomp) const;
void               AddComp(const sys_dynamic& sys, double pop=0);

void CheckComp(unsigned n) const;

int       NExProcs() const;
const ExchProc& ExProc(int iex) const;
void      ExProc(const ExchProc& pr, int iex);
double    Kex(int iex) const;
void      Kex(double K, int iex);
int       NCompsLHS(int iex) const;
int       NCompsRHS(int iex) const;

bool homonuclear(int   cmp=-1) const;
bool heteronuclear(int cmp=-1) const;

int HS(int comp=-1) const;
int LS(int comp=-1) const;
std::vector<int> HSs()   const;
std::vector<int> LSs()   const;

const std::string symbol(int comp, int spin) const;

void Omega(double freq);
void Omega(double freq, const std::string& iso);

double Omega( ) const;
double Omega(const std::string& iso) const;                             

void write(std::string& filename, std::string basename = "comp");

bool read(const std::string&  filename, int warn=2);
bool read(const ParameterSet& pset,     int warn=2);

std::string ask_read(int argc, char* argv[], int argn);
std::string ask_read(int argc, char* argv[], int argn, const std::string& def);

std::vector<std::string> SpinMapStrs(int exp) const;
std::vector<std::string> LHSStrs() const;
std::vector<std::string> RHSStrs() const;
std::vector<std::string> EXPStrs() const;

};
