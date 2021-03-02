/* acquire1D.h */

%{
#include "Level2/acquire1D.h"
%}

%include "HSLib/GenOp.i"
%include "std_string.i"

%feature("autodoc", "1" );

%rename(__assign__) acquire1D::operator=;


class acquire1D 
{
public:


acquire1D();				// Null constructor
acquire1D(const acquire1D& ACQ1);	// Self constructor
acquire1D(gen_op& det, gen_op& H);

acquire1D(gen_op& det, gen_op& H, double cutoff);

//acquire1D(gen_op& det, HSprop& U);
//acquire1D(gen_op& det, HSprop& U, double cutoff);

acquire1D(gen_op& D, super_op& L, gen_op& sigi, double cut=1.e-12);

//acquire1D(matrix& D, super_op& L, gen_op& sigi, double cut=1.e-12);

acquire1D(gen_op& D, super_op& L,               double cut=1.e-12);

//acquire1D(matrix& D, super_op& L,               double cut=1.e-12);
//acquire1D(gen_op& D, LSprop&   G,               double cut=1.e-12);

~acquire1D();
acquire1D& operator= (const acquire1D& ACQ1);

const super_op& L()	    const;
const gen_op&   D()      const;
//const matrix    S()      const;
//const matrix&   Sinv()   const;
const TTable1D& TTable() const;
void  Detector(const gen_op& detect);


row_vector T(const gen_op& sigmap, double tinc,      int npts);
row_vector T(const gen_op& sigmap, int npts,         double tinc);

void       T(const gen_op& sigmap, row_vector& data, double tinc);
row_vector F(const gen_op& sigmap, int npts, double Fst, double Ffi);

void F(const gen_op& sigmap, row_vector& data, double Fst, double Ffi);

row_vector FD(const gen_op& sigma, int npts, double Fst, double Ffi);
void       FD(const gen_op& sigma, row_vector& data, double Fst, double Ffi);


const TTable1D& table(const gen_op& sigmap);
const TTable1D& table() const;

const TTable1D table_snapshot(const gen_op& sigmap);
const TTable1D table_snapshot() const;

void offset(double     F,    int inHz=1);
void offset(double     F,    int tr, int inHz);
void FRscale(double    Fscf);
void FRscale(double    Fscf, int tr);
void Iscale(double     Iscf);
void Iscale(double     Iscf, int tr);
void broaden(double    LWR,          int inHz=1);
void broaden(double    LWR,  int tr, int inHz);
void resolution(double res);
void pcorrect(double   Wpivot, complex& P);
complex pcorrect(double& w0, double w1, int order=5);

double Wmax()  const;
double LWmax() const;
void   setSort(int sf);
void   setConv(int cf);
//void table(const gen_op& sigma, std::ostream& ostr);

//void table(std::ostream& ostr) const;

int ls()          const;		// Liouville space size
int size()        const;		// Acquisition size
int full_size()   const;		// Full acquisition size
int transitions() const;		// Number of transitions

//void parameters(matrix& mx, double& SW, double& LW,
//                                     double& dt, int N, int pf=0) const;
//std::ostream& print(std::ostream&  ostr) const;
//std::ostream& print(std::ostream&  ostr, gen_op& sigmap);
//std::ostream& printT(std::ostream& ostr, gen_op& sigmap,
//                                         double tinc, int npts=10, int P2P=0);
//friend std::ostream& operator << (std::ostream& ostr, acquire1D& ACQ);

void           write(const std::string& fn) const;
//std::ofstream& write(std::ofstream&     fp) const;

 
/// * These two were already commented out before creating .i file. *
///  void read(string& fn, gen_op& Op1);
///  void read(string& fn, basis& bs);

 
void read(const std::string& fn);
//std::ifstream& read(std::ifstream&     fp);

};

