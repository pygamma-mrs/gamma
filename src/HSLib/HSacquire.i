/* HSacquire.i */

%{
#include "HSLib/HSacquire.h"
%}

%feature("autodoc", "1" );

void acquire(gen_op& sig0,gen_op& D,gen_op& H,double td,int N,row_vector& fid, double CO=1.e-18);
void acquire(gen_op& sig0,gen_op& D,gen_op& U,          int N,row_vector& fid, double CO=1.e-18);
void acquire(gen_op& sig0,gen_op& D,HSprop& U,          int N,row_vector& fid, double CO=1.e-18);

void FID(gen_op&     sig0,gen_op& D,gen_op& H,double td,int N,row_vector& fid, double CO=1.e-18);
void FID(gen_op&     sig0,gen_op& D,gen_op& U,          int N,row_vector& fid, double CO=1.e-18);
void FID(gen_op&     sig0,gen_op& D,HSprop& U,          int N,row_vector& fid, double CO=1.e-18);

row_vector acquire(gen_op& sig0, gen_op& D, gen_op& H, double td, int N, double CO=1.e-18);
row_vector acquire(gen_op& sig0, gen_op& D, gen_op& U,            int N, double CO=1.e-18);
row_vector acquire(gen_op& sig0, gen_op& D, HSprop& U,            int N, double CO=1.e-18);

row_vector FID(gen_op&     sig0, gen_op& D, gen_op& H, double td, int N, double CO=1.e-18);
row_vector FID(gen_op&     sig0, gen_op& D, gen_op& U,            int N, double CO=1.e-18);
row_vector FID(gen_op&     sig0, gen_op& D, HSprop& U,            int N, double CO=1.e-18);





