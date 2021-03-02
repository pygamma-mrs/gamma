/* HSprop.i */
// Swig interface file.


%{
#include "HSLib/HSprop.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

%rename(__mul__)  HSprop::operator* const;
%rename(__imul__) HSprop::operator*=;
%rename(__iand__) complex::operator&=;

%rename(__eq__)  HSprop::operator== const;
%rename(__ne__)  HSprop::operator!= const;
%rename(__lt__)  HSprop::operator<  const;
%rename(__gt__)  HSprop::operator>  const;

%rename(__assign__) HSprop::operator=;


extern gen_op prop(const gen_op& ham, double time); 
extern void prop_ip(gen_op& U, double time);
extern gen_op evolve(const gen_op& sigma, const gen_op& ham, double time);
extern void evolve_ip(gen_op& sigma, const gen_op& ham, double time);
extern gen_op evolve(const gen_op& sigma, const gen_op& U);
extern void evolve_ip(gen_op& sigma, const gen_op& U);


class HSprop
{

public:


HSprop();					// Null constructor
HSprop(int HS);				// Identity constructor

HSprop(const gen_op& H, double tevol);
HSprop(const gen_op& H, double tevol, bool prop);
HSprop(const HSprop& U);

~HSprop();

HSprop& operator= (const HSprop& U1);

double time()   const;
double length() const;
int    dim()    const;
//matrix Mx()     const;
basis  Bs()     const;
int    HS()     const;
int    LS()     const;
gen_op Op()     const;
gen_op H()      const;

void SetEBR() const;
void SetBasis(const gen_op& Op);

gen_op evolve(const gen_op& Op) const;

HSprop operator *  (const HSprop& U) const;
HSprop &   operator *= (const HSprop& U);
HSprop &   operator &= (const HSprop& U);

HSprop sim_trans(const gen_op& Op);

void sim_trans_ip(const gen_op& Op);

HSprop Pow(int n) const;


// Next six lines were commented out in .h file 
// prior to making .i file.     DCT, 09/28/09
///friend gen_op prop(const gen_op& ham, const double time);
///friend void   prop_ip(   gen_op& ham, const double time);
///friend gen_op evolve(const gen_op& sigma, const gen_op& ham, double time);
///friend void   evolve_ip(   gen_op& sigma, const gen_op& ham, double time);
///friend gen_op evolve(const gen_op& sigma, const gen_op& U);
///friend void   evolve_ip(   gen_op& sigma, const gen_op& U);


bool operator== (const HSprop& U) const;
bool operator!= (const HSprop& U) const;
bool operator<  (const HSprop& U) const;
bool operator>  (const HSprop& U) const;
 
//std::ostream& print(std::ostream& ostr, int full=0) const;
//friend std::ostream &operator << (std::ostream &ostr, const HSprop& U);

};

