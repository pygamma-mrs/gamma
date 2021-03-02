/* decomp.i */

%{
#include "Level2/BaseDecomp.h"
%}

%feature("autodoc", "1" );

%rename(__assign__) decomp::operator=;

%include "std_vector.i"
%include "std_string.i"

class decomp
{
    
std::string         dname;		
int		    _NS;		
int                 _LS;		
std::vector<gen_op> BaseOps;		
std::vector<std::string> BaseNames;	
std::vector<std::string> BaseAltNames;	
std::vector<std::string> BaseCoeffNames;
std::vector<int>    BaseCoherences;	
std::vector<int>    BaseSpins;		
row_vector          BaseVals;	
std::vector<double> BaseCoefficients;
double              thresh;		
int                 prif;		
  
private:

void ODerror(int eidx, int noret=0) const;
volatile void ODfatal(int eidx)              const;

std::vector<int> sub_indices(int ils, int ssdim, int nss);
void sub_indices(int indices[], int ils, int ssdim, int nss);

void spin3halves(const spin_sys& sys);
void product_operators(const spin_sys& sys);

bool ChkIndex(int  i, bool warn=true) const;
bool ChkSize(int  ls, bool warn=true) const;

public:

decomp( );
decomp(const decomp &dec1);
decomp(const spin_sys& sys);
~decomp ( );
decomp& operator= (const decomp &dec);

void decompose(const gen_op& Op);
int size() const;			// Base Liouville space
int LS()   const;			// Base Liouville space
int HS()   const;			// Base Hilbert   space

std::vector<std::string> Names()           const;
std::vector<std::string> Names(int     m)  const;
void                Name(const std::string& name);
std::string              Name(          )  const;
std::string              OpName(int    i)  const;
std::string              AltOpName(int i)  const;
int                 MaxOpNameLen()    const;
int                 MaxOpAltNameLen() const;

int Coherence(int i) const;
int MaxCoherence()   const;

gen_op Op(const std::string& Opname) const;
gen_op Op(int i) const;

row_vector values()      const;
row_vector values(int m) const;
complex    value(int  i) const;

double bcoefficient(int i) const;

int index(const std::string& Opname) const;

std::vector<int> SortBySpins() const;

};

extern void PB_dec(const spin_sys &, const gen_op &);
