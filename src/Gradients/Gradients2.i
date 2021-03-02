/* Gradients2.i */

%{
#include "Gradients/Gradients2.h"
%}

%include "std_vector.i"

%feature("autodoc", "1" );

void Hzgrad(const sys_gradz& sys, gen_op& H0, gen_op* H);
std::vector<gen_op> Hzgrad(const sys_gradz& sys, gen_op& H0);

void Props(int NSS, gen_op* Hs, double t, gen_op* Us);
std::vector<gen_op> Props(std::vector<gen_op>& Hs, double t);
