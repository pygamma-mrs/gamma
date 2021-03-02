/* BlochSys.i */
// Swig interface file.

%{
#include "Bloch/BlochSys.h"
%}

%include "std_string.i"
%include "std_vector.i"

%feature("autodoc", "1" );

%rename(__assign__) BlochSys::operator=;


class BlochSys
{

std::vector <double>  Offsets;	// Mag. vector offsets
std::vector <Isotope> isotopes;	// Mag. vector associated isotope
std::vector <double>  R1rates;	// Mag. vector R1 rates
std::vector <double>  R2rates;	// Mag. vector R2 rates
std::vector <double>  Krates;		// Mag. vector exchange rates
std::vector <int>     Spins;		// Mag. vector associated spin
MagVec                _M;		// Magnetization vector

         void BSerror(int ei,                        int nr=0) const;
         void BSerror(int ei, const std::string& pn, int nr=0) const;
volatile void BSfatal(int ei) const;
volatile void BSfatal(int ei, const std::string& pn) const;

bool CheckR1s(const    std::vector<double>& R1s, bool warn=true) const;
bool CheckR2s(const    std::vector<double>& R2s, bool warn=true) const;
bool CheckIsos(const   std::vector<Isotope>& Is, bool warn=true) const;
bool CheckKs(const     std::vector<double>&  Ks, bool warn=true) const;
bool CheckSpins(int    ns1, int ns2,        bool warn=true) const;
bool CheckNorms(const  std::vector<double>&  Ns, bool warn=true) const;
bool CheckCoords(const coord_vec&       Ms, bool warn=true) const;

bool SetSystem(const ParameterSet& pset, int idx=-1, bool warn=true);

bool GetNSpins(const ParameterSet& pset,int& ns,bool wn=true) const;
bool GetNVects(const ParameterSet& pset,int& nm,bool wn=true) const;
bool SetVects(const  ParameterSet& pset,int   N,bool wn=true);

bool GetVect(const ParameterSet& pset, int i, double& v, Isotope& I,
                        double& R1, double& R2, int& Sp, bool warn=true) const;

bool GetW(const   ParameterSet& pset, int i, 
                                              double& v, bool warn=true) const;
bool GetIso(const ParameterSet& pset, int i,
                                             Isotope& I, bool warn=true) const;
bool GetR1(const  ParameterSet& pset, int i, 
                                             double& R1, bool warn=true) const;
bool GetR2(const  ParameterSet& pset, int i, 
                                             double& R2, bool warn=true) const;
bool GetSp(const  ParameterSet& pset, int i,
                                                int& Sp, bool warn=true) const;

bool SetExchange(const ParameterSet& pset, int nm, bool warn=true);

public:

BlochSys(int spins=0);
BlochSys(const BlochSys& sys);

BlochSys(double w, double R1, double R2);

BlochSys(const std::vector<double>& SH, 
               const std::vector<double>& R1s, const std::vector<double>& R2s);

BlochSys(const std::vector<double>& SH, const std::vector<Isotope>& Is,
               const std::vector<double>& R1s, const std::vector<double>& R2s);

BlochSys(const std::vector<double>& SH, 
               const std::vector<double>& R1s, const std::vector<double>& R2s,
                                                const std::vector<double>& Ks);

BlochSys(const std::vector<double>& SH, const std::vector<Isotope>& Is,
               const std::vector<double>& R1s, const std::vector<double>& R2s,
                                                const std::vector<double>& Ks);

BlochSys(const spin_system& sys, const RBasic& Rs);
BlochSys(const spin_system& sys, const matrix& Ks);
BlochSys(const spin_system& sys, const RBasic& Rs, const matrix& Ks);

BlochSys(const TTable1D& TT, const std::string& Iso=DEFISO);

~BlochSys();
BlochSys& operator= (const BlochSys& sys);

int NIso()         const;
int IsoMaxLength() const;

int NSpins()       const;

double R1(int i) const;
double T1(int i) const;
double R2(int i) const;
double T2(int i) const;
double LW(int i) const;

double MaxExchange() const;

std::vector<double> Norms() const;
void                Norms(const std::vector<double>& Ns);

double Norm(int i) const;
void   Norm(double nv, int i);

matrix H()                                       const;
matrix H(double gamB1, double w=0, double phi=0) const;

matrix B()                                       const;
matrix B(double gamB1, double w=0, double phi=0) const;
matrix R()                                       const;
matrix K()                                       const;

      MagVec  Meq() const;
const MagVec& Mo() const;
      MagVec  Mx() const;
      MagVec  My() const;
      MagVec  Mz() const;
      MagVec  Mss(const matrix& L, const matrix& R) const;
      MagVec  Mss(const matrix& L, const matrix& R,
                                                  const col_vector& Meq) const;

row_vector DetectMu(                 int u) const;
row_vector DetectMu(int k,           int u) const;
row_vector DetectMu(const std::string& I, int u) const;

row_vector DetectMx() const;
row_vector DetectMy() const;
row_vector DetectMz() const;

row_vector DetectMx(int i) const;
row_vector DetectMy(int i) const;
row_vector DetectMz(int i) const;

row_vector DetectMx(const std::string& I) const;
row_vector DetectMy(const std::string& I) const;
row_vector DetectMz(const std::string& I) const;

int size() const;

/*
virtual int write(const std::string &filename, int idx=-1, int warn=2) const;
virtual int write(std::ofstream& ofstr, int idx=-1, int warn=2) const; 
*/


bool read(const std::string& fn,    int idx=-1, int warn=2);
bool read(const ParameterSet& pset, int idx=-1, int warn=2);

std::string ask_read(int argc, char* argv[], int argn);
std::string ask_read(int argc, char* argv[], int argn,
                                                       const std::string& def); 

///std::ostream& print(std::ostream& out) const;
///friend std::ostream& operator<<(std::ostream& out, const BlochSys& sys);

};

