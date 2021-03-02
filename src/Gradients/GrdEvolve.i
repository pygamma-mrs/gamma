/* GrdEvolve.i */

%{
#include "Gradients/GrdEvolve.h"
%}

%include "std_vector.i"

%feature("autodoc", "1" );

std::vector<gen_op> evolve(gen_op& sigma0, std::vector<gen_op>& Hs, double t);
std::vector<gen_op> evolve(gen_op& sigma0, std::vector<gen_op>& Us);

std::vector<gen_op> evolve(gen_op& sigma0, std::vector<gen_op>& Hs, RBasic& R, double t);

std::vector<gen_op> evolve(std::vector<gen_op>& sigmas0, gen_op& H, double t);
std::vector<gen_op> evolve(std::vector<gen_op>& sigmas0, gen_op& U);

std::vector<gen_op> evolve(std::vector<gen_op>& sigma0, std::vector<gen_op>& H, double t);
std::vector<gen_op> evolve(std::vector<gen_op>& sigma0, std::vector<gen_op>& U);

std::vector<gen_op> evolve(std::vector<gen_op>& sigma0, std::vector<gen_op>& Hs, RBasic& R, double t);

