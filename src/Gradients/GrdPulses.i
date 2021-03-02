/* GrdPulses.i */

%{
#include "Gradients/GrdPulses.h"
%}

%include "std_vector.i"
%include "std_string.i"

%feature("autodoc", "1" );

std::vector<gen_op> Ixpulse_Us(const sys_gradz& sys,
                                                const std::string& I, double angle);

std::vector<gen_op> Iypulse_Us(const sys_gradz& sys, 
                                                const std::string& I, double angle);

std::vector<gen_op> Sxpuls_U(const sys_gradz& sys, std::vector<gen_op>& Hs,
                     const std::string& I, double offset, double tp, double angle);
std::vector<gen_op> Sypuls_U(const sys_gradz& sys, std::vector<gen_op>& Hs,
                     const std::string& I, double offset, double tp, double angle);

std::vector<gen_op> Gxpulse_U(const sys_gradz& sys, std::vector<gen_op>& Hs,
               const std::string& I, double tp, double angle, int N, double cutoff);

std::vector<gen_op> Gypulse_U(const sys_gradz& sys, std::vector<gen_op>& Hs,
               const std::string& I, double tp, double angle, int N, double cutoff);

std::vector<gen_op> SincPulseXUs(const sys_gradz& sys, std::vector<gen_op>& Hs,
                      const std::string& I, double tp, double angle, int N, int NN);

std::vector<gen_op> SincPulseYUs(const sys_gradz& sys, std::vector<gen_op>& Hs,
                      const std::string& I, double tp, double angle, int N, int NN);

