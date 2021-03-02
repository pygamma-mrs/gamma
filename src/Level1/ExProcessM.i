/* ExProcessM.h */
// Swig interface file.

%{
#include "Level1/ExProcessM.h"
%}

%include "std_string.i"
%include "std_vector.i"

%feature("autodoc", "1" );

%rename(__assign__) ExchProcM::operator=;

class ExchProcM
{

public:

double           KRate;		// Exchange rate (1/sec)
std::vector<int> Spins; 	// Spins involved in the exchange

void XPerror(int eidx,                           int noret=0) const;
void XPerror(int eidx, const std::string& pname, int noret=0) const;
volatile void XPfatal(int eidx)                               const;
volatile void XPfatal(int eidx, const std::string& pname)     const;

bool getExch(const ParameterSet& pset, int idx,
                     std::string& exch, bool warn=true) const;

bool parseExch(std::string& Exval,
                     std::vector<int>& sps, bool warn=true) const;

bool getComps(const ParameterSet& pset, int idx,
                     std::vector<int>& sps, bool warn=true) const;


bool getRate(const ParameterSet& pset, int idx,
                     double& rate, bool warn=true) const;

bool getXP(const ParameterSet& pset, double& rate,
                     std::vector<int>& sps, int idx, bool warn=true) const;

bool setXP(const ParameterSet& pset, int idx, bool warn=true);

bool CCheck(int comp, bool warn=true) const;
bool FCheck(          bool warn=true) const;

ExchProcM();
ExchProcM(const ExchProcM& proc);

ExchProcM(const ParameterSet& pset, int ip=-1, int warn=2);

ExchProcM& operator= (const ExchProcM& pr);

~ExchProcM();

double Kex() const;
void   Kex(double k);

int NComps()       const;
int NSpins()       const;
int Comp(int comp) const;

bool mixes(int    i, int j) const;
bool involves(int i)        const;

bool read(const std::string&  filename, int idx=-1, int warn=2);
bool read(const ParameterSet& pset,     int idx=-1, int warn=2);

std::string ExchStr() const;

//std::ostream& print(std::ostream& ostr, int full=0) const;
//friend std::ostream& operator<< (std::ostream& ostr, const ExchProcM& pro);

};


