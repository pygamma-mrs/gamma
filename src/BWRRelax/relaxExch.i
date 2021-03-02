/* relaxExch.i */

%{
#include "BWRRelax/relaxExch.h"
%}


class gen_op;	
class super_op;	
class sys_dynamic;

super_op Rex(const sys_dynamic& sys);

super_op Rex(const sys_dynamic& sys, gen_op Op);

