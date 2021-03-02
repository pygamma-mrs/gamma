/* sys_dynamic.i */
// Swig interface file.

%{
#include "LSLib/sys_dynamic.h"
%}

%include "std_vector.i"
%include "std_string.i"

%feature("autodoc", "1" );

%rename(__assign__) sys_dynamic::operator=;


class sys_dynamic: public spin_system, public coord_vec 
{

public:

sys_dynamic();
sys_dynamic(int spins);
sys_dynamic(const sys_dynamic &dsys);

sys_dynamic& operator= (const  sys_dynamic &dsys);

virtual      ~sys_dynamic ();


void   shifts(double shift=0);
void   shift(double shift, int spin);
void   shift(int spin, double shift);
double shift (int spin) const;
void   offsetShifts(double offset, int spin); 
void   offsetShifts(double offset, const std::string& iso); 

void PPM (int, double);

double     delz(int spin) const;
void       delz(int spin, double delzz);
double     Ceta(int spin) const; 
void       Ceta(int spin, double ceta); 
space_T    TC(int   spin) const;
void       TC(const space_T& A, int spin);
row_vector xiC_vector()   const; 
double     xiC(int spin)  const;
bool       CSA()          const;

void coords(const coord_vec& cvec, double cutoff=5.e-10);

bool Coord( ) const;

double  DCC(int   spin1, int spin2) const;
void    DCC(int   spin1, int spin2, double nu);
double  Ddelz(int spin1, int spin2) const;
void    Ddelz(int spin1, int spin2, double delzz);
double  Deta(int  spin1, int spin2) const;
void    Deta(int  spin1, int spin2, double Deta);
space_T AD(int    spin1, int spin2) const;
space_T AD(int dip) const;
int     dipoles() const;
int     dipole(int spin1, int spin2) const;
matrix  xiD_matrix() const;
bool    Dip( ) const;

double     QCC(int   spin) const;
void       QCC(int   spin, double nu);
double     Qdelz(int spin) const;
void       Qdelz(int spin, double delzz);
double     Qeta(int  spin) const;
void       Qeta(int  spin, double Qeta);

space_T    TQ(int    spin) const;
void       TQ(const space_T& A, int spin);

row_vector xiQ_vector( );
double     xiQ(int spin);
bool       Quad() const;

double     TR(int spin)  const;
double     tauR()        const;
row_vector xiR_vector( ) const;
double     xiR(int spin) const;

//operator ParameterSet( );

//friend void operator+= (ParameterSet& pset, sys_dynamic &dsys);

int setCoords(const ParameterSet& pset, int mand=0);

void setDip( );

void SetCSA(const ParameterSet& pset);

void setQuad(const ParameterSet& pset);

void setRand(const ParameterSet& pset);

void setTaus(const ParameterSet& pset, int mand=0);

bool setKs(const ParameterSet& pset, bool warn=true);

//void operator= (const ParameterSet& pset);

virtual void write(const std::string &filename);

virtual int read(const std::string &filename, int idx=-1, int warn=2);
virtual int read(const ParameterSet& pset,    int idx=-1, int warn=2);

virtual std::string ask_read(int argc, char* argv[], int argn);
virtual std::string ask_read(int argc, char* argv[], int argn, const std::string& def);

coord  taus() const;
double taux() const;
double tauy() const;
double tauz() const;

void taux(double tau);
void tauy(double tau);
void tauz(double tau);

double Kex(int p) const;
void   Kex(double K, int p);

//matrix Kex() const;

void   Kex_zero();
void   Kex(int i, int j, double K);
void   Kex(int N, int* Is, double K);

const std::vector<ExchProcM>& MExProcs() const;

std::vector<std::string> PtStrings(int w1=10, int w2=12, int digs=2) const;
std::vector<std::string> AQStrings(int w1=10, int w2=12, int digs=2) const;

//std::ostream& printAC(std::ostream& ostr) const;
//std::ostream& printAD(std::ostream& ostr) const;
//std::ostream& printAQ(std::ostream& ostr) const;
//std::ostream& printARDM(std::ostream& ostr) const;
//std::ostream& printTaus(std::ostream& ostr) const;
//std::ostream& printEX(std::ostream& ostr) const;
//virtual std::ostream& print(std::ostream& out) const;
//friend  std::ostream& operator<< (std::ostream& out, const sys_dynamic& dsys);
//std::ostream& print_D(std::ostream& out, int full=0);

};

