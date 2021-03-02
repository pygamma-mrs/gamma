/* LSprop.i *****************************************************-*-c++-*/

%{
#include "LSLib/LSprop.h"
%}

%feature("autodoc", "1" );

%rename(__assign__) LSprop::operator=;
%rename(__mul__)  LSprop::operator* const;
%rename(__imul__) LSprop::operator*=;
%rename(__iand__) LSprop::operator&=;

super_op R_prop(super_op& eLt, gen_op sigmaeq);
super_op R_prop(super_op& L, gen_op& sigmaeq, double t);

class LSprop
{
super_op GOp;				
double Gt;				
private:
void LSPerror(int eidx, int noret=0) const;
void volatile LSPfatal(int error) const;
public:
LSprop();				
LSprop(int LS);			
LSprop(const gen_op& H, double tevol);
LSprop(const gen_op& H, double tevol, bool prop);
LSprop(const HSprop& U);
LSprop(const super_op& L, double tevol);
LSprop(const super_op& L, const densop& sigma_ss, double tevol);
LSprop(const super_op& G);
LSprop(const LSprop& G);
~LSprop();
LSprop& operator= (const LSprop& G1);
double   time()   const;
double   length() const;
int      dim()    const;
int      HS()     const;
int      LS()     const;
super_op LOp()    const;
void     L(const super_op& LOp);
void     length(double t);

void SetEBR() const;
void SetBasis(const super_op& LOp);

gen_op evolve(const gen_op& Op);

LSprop operator *  (const LSprop& G) const;
LSprop &  operator *= (const LSprop& G);
LSprop &  operator &= (const LSprop& G);
};

