/* LSacquire.i */

%{
#include "LSLib/LSacquire.h"
%}

%feature("autodoc", "1" );

class row_vector;
class gen_op;		
class super_op;

void acquire(gen_op& sig0,             gen_op& D,super_op& L,double td,int N,row_vector& fid,double CO=1.e-18);

void FID(gen_op&     sig0,             gen_op& D,super_op& L,double td,int N,row_vector& fid,double CO=1.e-18);

void acquire(gen_op& sig0,gen_op& sigf,gen_op& D,super_op& L,double td,int N,row_vector& fid,double CO=1.e-18);

void FID(gen_op&     sig0,gen_op& sigf,gen_op& D,super_op& L,double td,int N,row_vector& fid,double CO=1.e-18);

void acquire(gen_op& sig0,             gen_op& D,super_op& L,row_vector& fid,double td,int np=0,double CO=1.e-18);

void FID(gen_op&     sig0,             gen_op& D,super_op& L,row_vector& fid,double td,int np=0,double CO=1.e-18);

void acquire(gen_op& sig0,gen_op& sigf,gen_op& D,super_op& L,row_vector& fid,double td,int np=0,double CO=1.e-18);

void FID(gen_op&     sig0,gen_op& sigf,gen_op& D,super_op& L,row_vector& fid,double td,int np=0,double CO=1.e-18);
    
row_vector acquire(gen_op& sig0,             gen_op& D,super_op& L,double td,int np,double CO=1.e-18);

row_vector FID(gen_op&     sig0,             gen_op& D,super_op& L,double td,int np,double CO=1.e-18);

row_vector acquire(gen_op& sig0,gen_op& sigf,gen_op& D,super_op& L,double td,int np,double CO=1.e-18);

row_vector FID(gen_op&     sig0,gen_op& sigf,gen_op& D,super_op& L,double td,int np,double CO=1.e-18);

void acquire(gen_op& sig, gen_op& D, super_op& G, row_vector &data, int np=0);

void FID(gen_op&     sig, gen_op& D, super_op& G, row_vector &data, int np=0);

/// KY - The following is in LSacquire.h file but no there is no 
// corresponing function in LSacquire.cc
/// void acquire(gen_op& sig, gen_op& D, LSprop&   G, row_vector &data, int np=0);

void FID(gen_op&     sig, gen_op& D, LSprop&   G, row_vector &data, int np=0);
    
row_vector acquire(gen_op& sig, gen_op& D, super_op& G, int np);
row_vector FID(gen_op&     sig, gen_op& D, super_op& G, int np);
row_vector acquire(gen_op& sig, gen_op& D, LSprop&   G, int np);
row_vector FID(gen_op&     sig, gen_op& D, LSprop&   G, int np);

void FIDx(gen_op sigma, gen_op& sigma0, gen_op& det, super_op &L,
	   			   row_vector &fid, double dt, int np=0);
void FIDrot(gen_op sigma, gen_op& sigma0, gen_op& det, super_op &L, gen_op& Fz,
                double Wrf, double time, row_vector &fid, double dt, int np=0);
