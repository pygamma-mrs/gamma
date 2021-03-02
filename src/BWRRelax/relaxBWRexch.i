/* relaxBWRexch.i */

%{
#include "BWRRelax/relaxBWRexch.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

%rename(__assign__) WBRExch::operator=;

void REXijkl(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);

void REXijk(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);

void REXkij(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);

void REXij(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);
 
void REXmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
               double w0, double w1, double w2, int l, int level, int autoc, int DFS=0, int het=0);

void REXmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
         double w0, double w1, double w2, int l, int level, int autoc, int DFS, int het,
                               gen_op& Fz11, double W11, gen_op& Fz12, double W12,
                               gen_op& Fz21, double W21, gen_op& Fz22, double W22);

void REXrfijkl(super_op &LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);

void REXrfijk(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);
void REXrfkij(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);

void REXrfij(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
         double Wrflab, matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level, int DFS=0);

void REXrfmumu(super_op& LOp, gen_op* T1s, gen_op* T2s, matrix* J12,
                  double* J, double* w, int rank, int level, int autoc, int DFS=0, int het=0);

void REX_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12);

complex REX_4(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
		 	         int rank, int a, int b, int aa, int bb);

void REX_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix& J12,
                               gen_op& Fz11, double W11, gen_op& Fz12, double W12,
                               gen_op& Fz21, double W21, gen_op& Fz22, double W22);

complex REX_4(int hs, gen_op* T1s, gen_op* T2s, matrix& J12,
                        int rank, int a, int b, int aa, int bb,
                               gen_op& Fz1, double W1, gen_op& Fz2, double W2);

void REX_3(super_op& LOp, double* w, int rank, gen_op* T1s, gen_op* T2s,
					 matrix& J12, double cutoff=1.e-2);

void REXrf_4(super_op& LOp, int rank, gen_op* T1s, gen_op* T2s, matrix* J12a);

complex REXrf_4(int hs, gen_op* T1s, gen_op* T2s, matrix* J12,
		 	         int rank, int a, int b, int aa, int bb);

void REXrf_3(super_op& LOp, double* w, int rank, gen_op* T1s,
                                         gen_op* T2s, matrix* J12, double cutoff=1.e-6);




class WBRExch 
  {
  int DD; 	      
  int CC; 	      
  int QQ; 	      
 
  int DDdfs; 	      
  int CCdfs; 	      
  int QQdfs; 	      
 
  int DC; 	      
  int DCdfs; 	      

  int DQ; 	      
  int DQdfs; 	      

  int QC; 	      
  int QCdfs; 	      

  int level; 	      
  int type; 	      


void WBRerror(int eidx, int noret=0) const;
 
void WBRerror(int eidx, const std::string& pname, int noret=0) const;

volatile void WBRfatality(int eidx) const;

volatile void WBRfatality(int eidx, const std::string& pname) const;

void assign(const ParameterSet& pset, int DF=1, int CF=1, int QF=1);

public:

WBRExch();

WBRExch(const WBRExch& WBRE);

~WBRExch();

WBRExch& operator= (const WBRExch &WBRE);

void Level(int i);

int Level() const;
 
void Type(int i);

int Type() const;

void Dip(int i=1);

void DipDFS(int i=1);

void DipCSA(int i=1);

void DipCSADFS(int i=1);

void DipQuad(int i=1);

void DipQuadDFS(int i=1);

void CSA(int i=1);

void CSADFS(int i=1);

void CSADip(int i=1);

void CSADipDFS(int i=1);

void CSAQuad(int i=1);

void CSAQuadDFS(int i=1);

void Quad(int i=1);

void QuadDFS(int i=1);

void QuadDip(int i=1);

void QuadDipDFSQuad(int i=1);

void QuadCSA(int i=1);

void QuadCSADFS(int i=1);

double LWhh(const sys_dynamic& sys, const std::string& Iso);

//operator ParameterSet( ) const;

void SetZero();
void SetLevel(const ParameterSet& pset);
void SetType(const  ParameterSet& pset);
void SetDip(const   ParameterSet& pset);
void SetSA(const    ParameterSet& pset);
void SetQuad(const  ParameterSet& pset);
void SetDCX(const   ParameterSet& pset);

void SetDQX(const ParameterSet& pset);

void SetQCX(const ParameterSet& pset);

WBRExch& operator= (const ParameterSet& pset);

void prepQuad(const sys_dynamic& sys, matrix& Xis, spin_T* Ts, space_T* As) const;
 
void read(const std::string& filename);

void read(const std::string& filename, const sys_dynamic& sys);

void ask_read(int argc, char* argv[], int argn);

void ask_read(int argc, char* argv[], int argn, const sys_dynamic& sys);

void ask(int argc, char* argv[], int& argn);

super_op REX(const sys_dynamic& sys, gen_op& Ho, int fext=0);

super_op REXrf(const sys_dynamic& sys, gen_op& Heff, double Wrf, int fext=0);

};
