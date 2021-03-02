/* MutExch.h */
// Swig Interface file

%{
#include "Level2/MutExch.h"
%}

%include "std_vector.i"

%feature("autodoc", "1" );

//matrix Kex(const spin_sys& sys, const ExchProcM& XP);

super_op Kex(const spin_sys& sys,
                          const std::vector<ExchProcM>& XPs, const basis& Bs);

super_op Kex(const sys_dynamic& sys, const basis& Bs);

