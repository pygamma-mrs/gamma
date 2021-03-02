/* MultiLib.i */


%{
#include "MultiSys/MultiLib.h"
%}

%include "std_vector.i"

%feature("autodoc", "1" );

//gen_op multize(gen_op  OpFct(const spin_system &sys), const multi_sys &msys); 

gen_op multize(gen_op& Op, const multi_sys &msys, int icomp);

//gen_op multize(spin_op op(const spin_sys&),             const multi_sys& msys);
//gen_op multize(spin_op op(const spin_sys&,  const std::string&),
//                                    const std::string&, const multi_sys& msys);
//gen_op multize(spin_op op(const spin_sys&, double),
//                             double beta, const multi_sys &msys, int icomp=-1);
//gen_op multize(gen_op op(const spin_sys&, const gen_op&, double),
//	                      gen_op&, double, const multi_sys&, int icomp=-1);
//gen_op multize(gen_op op(const spin_sys&, const gen_op&, int, double),
//                       gen_op&, int, double, const multi_sys&, int icomp=-1);
//gen_op multize(gen_op op(const spin_sys&,const gen_op&,const std::string&, double),
//               gen_op&, const std::string&, double, const multi_sys&, int icomp=-1);
//gen_op multize(gen_op op(const spin_sys&, int, double),
//                  int nspin, double beta, const multi_sys &msys, int icomp=-1);
//gen_op multize(gen_op op(const spin_sys&, const std::string& iso, double),
//     const std::string& iso, double beta, const multi_sys &msys, int icomp=-1);
//gen_op multize(gen_op op(const spin_sys&, double),
//                             double beta, const multi_sys &msys, int icomp=-1);
//gen_op multize(gen_op op(const spin_sys&, double, double),
//                 double phi, double beta, const multi_sys &msys, int icomp=-1);
 
super_op multize(super_op& SOp, const multi_sys &msys, int icomp);

//super_op multize(super_op SOpFct(const gen_op&), const gen_op& Op, 
//                                                const multi_sys& msys);
//super_op multize(super_op sop(const sys_dynamic&,gen_op&,int,int), gen_op&,
//                                                  int, int, const multi_sys&);
//super_op multize(super_op sop(const sys_dynamic&,gen_op&,int), gen_op&,
//                                                         int, const multi_sys&);

basis D_basis(const multi_sys& msys);
std::vector<double> qStateLS(const multi_sys& msys, int I);
row_vector LS_qState_bra(const multi_sys& msys, int i);
row_vector LS_qState_ket(const multi_sys& msys, int i);
