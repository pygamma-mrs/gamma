/* MultiExch.i */

%{
#include "MultiSys/MultiExch.h"
%}

%feature("autodoc", "1" );

super_op Xnm( const    multi_sys& msys);
matrix   Xnmp(const    multi_sys& msys, int p);
Xnmpblk(const multi_sys& msys, const ExchProc& Pro, matrix& Xp,
            double K, int cmpI, int cmpJ, int Io, int Iend, int Jo, int Jend);
Xnmpelem(const multi_sys& msys, const ExchProc& Pro, const row_vector& braI, const row_vector& ketI, const row_vector& braJ, const row_vector& ketJ, int cmpI, int cmpJ, int& hsnorm);

super_op Xm(const multi_sys& msys);

void Xnm(std::ostream&     ostr, const multi_sys& sys);
void Xnmp(std::ostream&    ostr, const multi_sys& msys, int p);
void Xnmpdblk(std::ostream& ostr, const multi_sys& msys, double K, int Io, int Iend);
void Xnmpblk(std::ostream& ostr,  const multi_sys& msys, const ExchProc& Pro,
             double K, int cmpI, int cmpJ, int Io, int Iend, int Jo, int Jend);

super_op XXnm( const    multi_sys& msys);
