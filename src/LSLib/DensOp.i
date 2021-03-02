// DensOp.i

%{
#include "LSLib/DensOp.h"
%}

%feature("autodoc", "1" );

%rename(__assign__) densop::operator=;

gen_op SigmaEq(const spin_sys& sys);
gen_op SigmaSS(const spin_sys& sys, super_op& L, super_op& R, int wrn);
gen_op SigmaSS(super_op& L, super_op& R, gen_op& seq, int wrn);

class densop : public gen_op
{
double Sigmat;
 
private:
 
void SIGMAerror(int eidx, int noret=0) const;
 
void volatile SIGMAfatality(int error) const;
 
public:
 
densop();
         
densop(const spin_sys& sys);
 
densop(const spin_sys& sys, super_op& L, super_op& R);

densop(super_op& L, super_op& R, gen_op& sigmaeq);
 
densop(gen_op& Op, double tevol);

densop(const densop& Sigma);
         
~densop();

densop & operator= (const densop& Sigma1);

double length() const;
 
void SetTrace(double tr);

};
