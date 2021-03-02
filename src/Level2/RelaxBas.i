/* RelaxBas.i */

%{
#include "Level2/RelaxBas.h"
%}

%include "std_vector.i"

%feature("autodoc", "1" );

%rename(__assign__) RBasic::operator=;

class RBasic
{
std::vector<double> R1rates;
std::vector<double> R2rates;
matrix              FZxNS;
col_vector          Fzs;	
std::vector<double> Cinf;	
std::vector<double> Csig;	
std::vector<gen_op> Izis;	
gen_op              H0;		
gen_op              Det;	
gen_op              SigInf;	
matrix              R2mx;	
 

void RBasErr(int   eidx,                           int noret=0) const;
void RBasErr(int   eidx, const std::string& pname, int noret=0) const;
void RBasFatal(int eidx)                                        const;
void RBasFatal(int eidx, const std::string& pname)              const;

bool GetNPoints(const ParameterSet& pset, int& ns, bool warn=true);
 
void CheckR1(double r1);
void CheckR2(double r2);
void CheckT1(double t1);
void CheckT2(double t2);
void CheckLW(double lw);

bool CheckSpin(double  i, bool warn=true) const;
bool CheckSpins(double N, bool warn=true) const;
bool CheckHS(double    H, bool warn=true) const;

bool SetIzs(const    spin_sys& sys,   bool warn=true);
bool SetCinfs(const  gen_op&   sigma, bool warn=true);
bool SetCsigs(const  gen_op&   sigma, bool warn=true);
bool TestLong(double cutoff=1.e-9);

void SetR2Mx();
void ZeroR2Mx();

public:

RBasic();
RBasic(const row_vector& vx);
RBasic(const RBasic& RB);
~RBasic();

RBasic operator= (const RBasic& RB);

int    spins() 	     const;	
int    HS()          const;
double T1(int i)     const;
double T2(int i)     const;
double R1(int i)     const;
double R2(int i)     const;
double LW(int i)     const;
double RB(int i, int type) const;

void T1(double val, int i);	
void T2(double val, int i);	
void R1(double val, int i);	
void R2(double val, int i);	
void LW(double val, int i);	
void RB(double val, int i, int type);

std::vector<double> T1s() const;	
std::vector<double> T2s() const;	
std::vector<double> R1s() const;	
std::vector<double> R2s() const;	
std::vector<double> LWs() const;	
std::vector<double> RBRates(int type) const;

bool   SetSystem(const spin_sys& sys, int warn=2);
bool   SetH0(const     gen_op&     H, int warn=2);
bool   SetDet(const    gen_op&     D, int warn=2);
bool   SetSigInf(const gen_op&     S, int warn=2);
matrix R2Mx();
matrix R2LOp();
  

matrix R1LOp();

gen_op     SigmaT1(const  gen_op& sigma);
gen_op     SigmaT2(const  gen_op& sigma);
col_vector SigmaC(const   gen_op& sigma);
col_vector SigmaCEq(const gen_op& sigeq);
col_vector SigmaCEq();
matrix     RC();
matrix     HC(const gen_op& H);
matrix     HC();
gen_op     Sigma(const col_vector& sigmaC);

double ReadT2(const  ParameterSet& pset, int sp, int idx=-1, int pf=0);
double ReadT1(const  ParameterSet& pset, int sp, int idx=-1, int pf=0);
double ReadLW(const  ParameterSet& pset, int sp, int idx=-1, int pf=0);
double ReadR2(const  ParameterSet& pset, int sp, int idx=-1, int pf=0);
double ReadR1(const  ParameterSet& pset, int sp, int idx=-1, int pf=0);
double ReadPar(const ParameterSet& P,int I,int t,int idx=-1, int pf=0);

std::vector<double> ReadT2s(const ParameterSet& p,int N,int idx=-1,int pf=0);
std::vector<double> ReadT1s(const ParameterSet& p,int N,int idx=-1,int pf=0);
std::vector<double> ReadLWs(const ParameterSet& p,int N,int idx=-1,int pf=0);
std::vector<double> ReadR2s(const ParameterSet& p,int N,int idx=-1,int pf=0);
std::vector<double> ReadR1s(const ParameterSet& p,int N,int idx=-1,int pf=0);
std::vector<double> ReadPars(const ParameterSet& p,
                                          int N, int type, int idx=-1,int pf=0);
int read(const std::string &filename, int idx=-1, int warn=2);
int read(const ParameterSet& pset,    int idx=-1, int warn=2);
 
gen_op Evolve(const spin_sys& sys, const gen_op& sig0, double t);
gen_op Evolve(const gen_op& sigmap, double t);

void FID(const gen_op& sigmap, double td, row_vector& fid, int N=0);

void       FID(const spin_sys& sys, const gen_op& sigmap,
                                          double td, row_vector& fid, int N=0);
row_vector FID(const spin_sys& sys, const gen_op& sigmap,
                                          double td,                  int N);

row_vector FID(const gen_op& sigmap, double td, int N);

std::vector <double>FzCoeffs(const spin_sys& sys, const gen_op& sigma);
};
