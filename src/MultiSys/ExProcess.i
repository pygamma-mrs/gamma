/* ExProcess.i */

%{
#include "MultiSys/ExProcess.h"
%}

%include "std_vector.i"
%include "std_string.i"

%feature("autodoc", "1" );

%rename(__assign__) ExchProc::operator=;

class ExchProc
  {
public:

  double               KRate;	
  std::vector<int>     LHSComps;
  std::vector<int>     RHSComps;
  std::vector<SpinMap> SpinMaps;


void XPerror(int eidx,                           int noret=0) const;
void XPerror(int eidx, const std::string& pname, int noret=0) const;
volatile void XPfatal(int eidx)                               const;
volatile void XPfatal(int eidx, const std::string& pname)     const;

bool getExch(const ParameterSet& pset, int idx,
                                          std::string& exch, bool warn=true) const;
bool parseExch(std::string& Exval,
               std::vector<int>& lhs, std::vector<int>& rhs, bool warn=true) const;
bool getComps(const ParameterSet& pset, int idx,
               std::vector<int>& lhs, std::vector<int>& rhs, bool warn=true) const;


bool getRate(const ParameterSet& pset, int idx,
                                               double& rate, bool warn=true) const;

bool getMappings(const ParameterSet& pset, int idx,
                                std::vector<SpinMap>& smaps, bool warn=true) const;

bool getXP(const ParameterSet& pset, double& rate,
  std::vector<int>& lhsc, std::vector<int>& rhsc, std::vector<SpinMap>& smaps,
                                               int idx, bool warn=true) const;
bool setXP(const ParameterSet& pset, int idx, bool warn=true);

bool CheckLHS(int comp, bool warn=true) const;
bool CheckRHS(int comp, bool warn=true) const;

ExchProc();
ExchProc(const ExchProc& proc);

ExchProc(const std::string& PROC, double Kex=0, int maxcomp=20);

ExchProc(const ParameterSet& pset, int ip=-1, int warn=2);

ExchProc& operator=(const ExchProc& pr);
~ExchProc();

ExchProc(int N_lhs, int N_rhs);

void intra_default(int ic1, int ic2, int nspins, double k);

double Kex() const;
void   Kex(double k);

int LHSComp(int comp) const;
int RHSComp(int comp) const;

int NCompsLHS() const;
int NCompsRHS() const;

bool mixes(int     comp, int comp1) const;
bool CompInLHS(int comp)            const;
bool CompInRHS(int comp)            const;
bool involves(int  comp, int lr=0)  const;

     int      NSpinMaps() const;
const SpinMap& SMap(int i) const;
     bool     SMap(int comp1, int sp1, int& comp2, int& sp2) const;
     void     add_pair(SpinMap);
     bool     mapped(int comp1, int s1, int comp2, int s2) const;
     bool     mapped(int comp1,         int comp2) const;

void mapping(const std::string& spair);
bool read(const std::string&  filename, int idx=-1, int warn=2);
bool read(const ParameterSet& pset,     int idx=-1, int warn=2);

static char Label(int i);
std::string LHSStr()     const;
std::string RHSStr()     const;

std::vector<std::string> SpinMapStrs() const;

std::ostream& print(std::ostream& ostr, int full=0) const;


};


