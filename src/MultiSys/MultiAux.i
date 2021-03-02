/* MultiAux.h */

%{
#include "MultiSys/MultiAux.h"
%}

%include "std_vector.i"
%include "std_string.i"

%feature("autodoc", "1" );

%rename(__assign__) spin_pair::operator=;
%rename(__assign__) process::operator=;

class spin_pair
  {
  public:

    int sub1;	
    int sp1;	
    int sub2;	
    int sp2;	

inline spin_pair() {};

inline spin_pair(int, int, int, int);

spin_pair(const spin_pair& Sp);

spin_pair(const std::string& SSP);   

inline spin_pair& operator = (const spin_pair& Sp);

inline ~spin_pair() {};

int Sub1() const;
int Sub2() const;
int Spin1() const;
int Spin2() const;
inline void print() const;
};


class process
  {
  public:

  double krate;			
  int* comp_lhs;		
  int n_lhs;			
  int* comp_rhs; 		
  int n_rhs;			
  int npair;			
  spin_pair* link;		

void process::XPerror(int eidx,                           int noret=0) const;
void process::XPerror(int eidx, const std::string& pname, int noret=0) const;
volatile void process::XPfatal(int eidx)                               const;

bool process::getExch(const ParameterSet& pset, int idx,
                                           std::string& exch, bool warn=true) const;

bool process::parseExch(std::string& Exval,
                std::vector<int>& lhs, std::vector<int>& rhs, bool warn=true) const;

bool process::getRate(const ParameterSet& pset, int idx,
                                               double& rate, bool warn=true) const;

bool process::getXP(const ParameterSet& pset, int idx, bool warn=true) const;
bool process::setXP(const ParameterSet& pset, int idx, bool warn=true) const;
inline process::process();
inline process::process(int N_lhs, int N_rhs);
process::process(const process& proc);
process(std::string& PROC, double Kex=0, int maxcomp=20);

process(const ParameterSet& pset, int ip=-1, int warn=2);

process& process::operator=(const process& pr);
void intra_default(int ic1, int ic2, int nspins, double k);
inline double get_k() const;

inline void set_k(double k);

int lhsindex(int comp);      

int rhsindex(int comp);      

int mixes(int ic1, int ic2);

int involves(int ic1, int lr=0);

int pairs() const;

void add_pair(spin_pair);

spin_pair get_pair(int) const;

int mapped(int c1, int s1, int c2, int s2) const;

void mapping(const std::string& spair);

bool process::read(const std::string&  filename, int idx=-1, int warn=2);
bool process::read(const ParameterSet& pset,     int idx=-1, int warn=2);

};				

inline spin_pair::spin_pair(int subA, int spA, int subB, int spB)
  {
  sub1 = subA; 			
  sp1  = spA;			
  sub2 = subB;			
  sp2  = spB;			
  }

inline void spin_pair::operator = (const spin_pair& Sp)
  {
  sub1 = Sp.sub1; 		
  sp1  = Sp.sp1; 		
  sub2 = Sp.sub2; 		
  sp2  = Sp.sp2; 		
  } 

inline void spin_pair::print() const
  {
  std::cout << "Component " <<sub1 << " spin #" <<sp1
            <<" --- to --- component ";
  std::cout << sub2 << " spin #" << sp2 << "\n";
  }

inline process::process()

  {
  krate = 0.0;			
  comp_lhs = NULL;		
  n_lhs = 0;			
  comp_rhs = NULL;		
  n_rhs = 0;
  npair = 0;		
  link = NULL;		
  }

inline process::process(int N_lhs, int N_rhs)

  {
  krate = 0.0;		         
  comp_lhs = new int[N_lhs];	
  n_lhs = N_lhs ;		
  comp_rhs = new int[N_rhs];	
  n_rhs = N_rhs;		
  npair = 0;			
  link = NULL;			
  }

inline double process::get_k() const { return (krate); }

inline void process::set_k(double k) { krate = k; }
