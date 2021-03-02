/* relaxanalyze.i */

%{
#include "BWRRelax/relaxanalyze.h"
%}


void RDDel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,int a, int aa, int b, int bb,int DFS=0, int Windex=-1, int Sindex=0);
 
void RSSel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2, int a, int aa, int b, int bb, int DFS=0, int Windex=-1, int Sindex=0);
 
void RDSel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2, int a, int aa, int b, int bb, int DFS=0, int Windex=-1, int Sindex=0);
 
void RSDel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2, int a, int aa, int b, int bb, int DFS=0, int Windex=-1, int Sindex=0);

void RRRel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2, int a, int aa, int b, int bb,int DFS=0, int Windex=-1, int Sindex=0);
 
void RQQel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2, int a, int aa, int b, int bb, int DFS=0, int Windex=-1, int Sindex=0);
 
void RQSel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2, int a, int aa, int b, int bb, int DFS=0, int Windex=-1, int Sindex=0);
 
void RSQel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2, int a, int aa, int b, int bb, int DFS=0, int Windex=-1, int Sindex=0);
 
void RQDel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2, int a, int aa, int b, int bb, int DFS=0, int Windex=-1, int Sindex=0);
 
void RDQel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2, int a, int aa, int b, int bb, int DFS=0, int Windex=-1, int Sindex=0);
 
void Rijkl_el(const sys_dynamic& sys, gen_op& Ho, int rank,  spin_T* T1, spin_T* T2, std::string& Mlabel, std::string* line1, std::string* line2, std::string* line3, int* signs, int& Jterms, int& Lterms, int a, int aa, int b, int bb, int DFS=0, int Windex=-1, int Sindex=0);

void Rij_el(const sys_dynamic& sys, gen_op& Ho, int rank, spin_T* T1, spin_T* T2, std::string& Mlabel, std::string* line1, std::string* line2, std::string* line3, int* signs, int& Jterms, int& Lterms, int a, int aa, int b, int bb, int DFS=0, int Windex=-1, int Sindex=0);
 
void Rijk_el(const sys_dynamic& sys, gen_op& Ho, int rank, spin_T* T1, spin_T* T2, std::string& Mlabel, std::string* line1, std::string* line2, std::string* line3, int* signs, int& Jterms, int& Lterms, int a, int aa, int b, int bb, int DFS=0, int Windex=-1, int Sindex=0);
 
void Rkij_el(const sys_dynamic& sys, gen_op& Ho, int rank, spin_T* T1, spin_T* T2, std::string& Mlabel, std::string* line1, std::string* line2, std::string* line3, int* signs, int& Jterms, int& Lterms, int a, int aa, int b, int bb,int DFS=0, int Windex=-1, int Sindex=0);
void Rel_12(int hs, int rank, gen_op* T1s, gen_op* T2s, int& Jterms, complex* strsJ, int* trnsJ, int& Lterms, complex* strsL, int* trnsL, std::string* wlabs, int a, int aa, int b, int bb, double cutoff=1.e-4);

void Rel_12(int hs, int rank, gen_op* T1s, gen_op* T2s, int& nterms, complex* strs, int* trns,int a, int aa, int b, int bb, double cutoff=1.e-4);

void Rel_12_condense(int hs, int ntermsi, int& nterms, complex* strs, int* trns, int anti=0, double cutoff=1.e-4);
 
void Rel_12_condense(int hs, int ntermi, int& nterms, complex* strs, int* trns, std::string* wlabs, int anti=0);
 
void Rel(int ntermi, int& nterms, int npairs, complex* strs, int* trns, std::string* wlabs, int* cont, std::string* spns, std::string* Jlbs, std::string& Mlabel, std::string* line1, std::string* line2, std::string* line3, int* signs);

void Spin_labels(std::string* Lbls, const spin_sys& sys, int index=0);

void W_labels(std::string* Wlabels, const spin_sys& sys, gen_op &Op, int index=-1);

void Elem_labels(std::string* Lbls, std::string& R, std::string& M, std::string& S,int a, int aa, int b, int bb, int la=0, int laa=0, int lb=0, int lbb=0);

void Rel_clean(gen_op* T1s, gen_op* T2s, int rank);

void Rel(std::ostream& ostr, int nterms, std::string* line1, std::string* line2,
 std::string* line3, int* signs, std::string* Elabel, int add=0, int ncols=4);

void sort(int* indx, matrix& mx, int k, int type=0, int colf=0);
