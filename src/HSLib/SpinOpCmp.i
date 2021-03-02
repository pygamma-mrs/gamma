/* SpinOpCmp.i */

%{
#include "HSLib/SpinOpCmp.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

void SOpCmperror(int eidx, int noret=0);
void volatile SOpCmpfatal(int eidx);

spin_op Iu(const   spin_sys &sys, int spin, int type);
spin_op Ie(const   spin_sys &sys, int spin);
spin_op Iz(const   spin_sys &sys, int spin);
spin_op Ix(const   spin_sys &sys, int spin);
spin_op Iy(const   spin_sys &sys, int spin);
spin_op Ip(const   spin_sys &sys, int spin);
spin_op Im(const   spin_sys &sys, int spin);
spin_op Ia(const   spin_sys &sys, int spin);
spin_op Ib(const   spin_sys &sys, int spin);
spin_op Ipol(const spin_sys &sys, double m, int spin);

spin_op Fe(const spin_sys &sys, int spin);
spin_op Fx(const spin_sys &sys, int spin);
spin_op Fy(const spin_sys &sys, int spin);
spin_op Fz(const spin_sys &sys, int spin);
spin_op Fp(const spin_sys &sys, int spin);
spin_op Fm(const spin_sys &sys, int spin);
spin_op Fa(const spin_sys &sys, int spin);
spin_op Fb(const spin_sys &sys, int spin);
spin_op Fpol(const spin_sys &sys, double m, int spin);

spin_op Fe(const spin_sys &sys);
spin_op Fx(const spin_sys &sys);
spin_op Fy(const spin_sys &sys);
spin_op Fz(const spin_sys &sys);
spin_op Fp(const spin_sys &sys);
spin_op Fm(const spin_sys &sys);
spin_op Fa(const spin_sys &sys);
spin_op Fb(const spin_sys &sys);
spin_op Fpol(const spin_sys &sys, double m);

spin_op Fe(const spin_sys &sys, const std::string& iso);
spin_op Fx(const spin_sys &sys, const std::string& iso);
spin_op Fy(const spin_sys &sys, const std::string& iso);
spin_op Fz(const spin_sys &sys, const std::string& iso);
spin_op Fp(const spin_sys &sys, const std::string& iso);
spin_op Fm(const spin_sys &sys, const std::string& iso);
spin_op Fa(const spin_sys &sys, const std::string& iso);
spin_op Fb(const spin_sys &sys, const std::string& iso);
spin_op Fpol(const spin_sys &sys, double m, const std::string& iso);

spin_op Fe(const spin_sys& sys, const flagvec& sflags);
spin_op Fx(const spin_sys& sys, const flagvec& sflags);
spin_op Fy(const spin_sys& sys, const flagvec& sflags);
spin_op Fz(const spin_sys& sys, const flagvec& sflags);
spin_op Fp(const spin_sys& sys, const flagvec& sflags);
spin_op Fm(const spin_sys& sys, const flagvec& sflags);
spin_op Fa(const spin_sys& sys, const flagvec& sflags);
spin_op Fb(const spin_sys& sys, const flagvec& sflags);
//spin_op Fpol(const spin_sys& sys, double m, const flagvec& sflags);
 
spin_op Fe_sp(const spin_sys& sys);
spin_op Fx_sp(const spin_sys& sys);
spin_op Fy_sp(const spin_sys& sys);
spin_op Fz_sp(const spin_sys& sys);
spin_op Fp_sp(const spin_sys& sys);
spin_op Fm_sp(const spin_sys& sys);
spin_op Fa_sp(const spin_sys& sys);
spin_op Fb_sp(const spin_sys& sys);
spin_op Fpol_sp(const spin_sys& sys, double m);

spin_op Faxis(const    spin_sys& sys, int spin,             char axis);
spin_op Faxis(const    spin_sys& sys, const std::string& I, char axis);
spin_op Faxis_sp(const spin_sys &sys,                       char axis);
spin_op Faxis(const    spin_sys &sys,                       char axis);
//spin_op Faxis(const    spin_sys &sys, const flagvec& flags, char axis);

spin_op Fpol_gen_new(const spin_sys& sys,                      double m);
spin_op Fpol_gen_new(const spin_sys& sys,const std::string& I, double m);
//spin_op Fpol_gen_new(const spin_sys& sys,const flagvec& flags, double m);
spin_op Fpol_gen(const     spin_sys& sys,                      double m);

spin_op Ipdt(const spin_sys &sys, std::string name);
spin_op Fpdt(const spin_sys &sys, std::string name);

