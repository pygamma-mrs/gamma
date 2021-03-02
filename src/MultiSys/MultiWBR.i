/* MultiWBR.i */

%{
#include "MultiSys/MultiWBR.h"
%}

%feature("autodoc", "1" );

super_op RQQ(const multi_sys& msys, gen_op& H, int type=0, int level=4);
super_op RCC(const multi_sys& msys, gen_op& H, int type=0, int level=4);
super_op RDD(const multi_sys& msys, gen_op& H, int type=0, int level=4);
super_op RCQ(const multi_sys& msys, gen_op& H, int level=4);
super_op RQC(const multi_sys& msys, gen_op& H, int level=4);
