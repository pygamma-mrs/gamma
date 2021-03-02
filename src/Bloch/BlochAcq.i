/* BlochAcq.i */
// Swig interface file.

%{
#include "Bloch/BlochAcq.h"
%}

%include "std_vector.i"

%feature("autodoc", "1" );

%rename(__assign__) BlochAcq::operator=;


class BlochAcq 
{

int BS;			// Bloch dimension

matrix L;			// System "Liouvillian" matrix
matrix Sm1;			// L matrix inverse eigenvectors

row_vector det;		// The detection operator
col_vector Minf;		// Infinite time magnetization vector
complex trinf;		// Trace at infinite time Tr(det*Minf)

int pos;			// Acquire1D dimensions (pos <= ls)
row_vector A;			// A array (detector equivalent)
row_vector B;			// B array (time propagator equivalent)
std::vector<int> I;		// Array for indexing
double DCUTOFF;		// Detection cut value

TTable1D TTab;		// Transitions table for 1D spectrum
col_vector sigmap;		// Prepared density operator
double Icut;			// Transitiion intensity cut


void ACQerror(int eidx, int noret=0) const;

volatile void ACQfatality(int eidx=0) const;

void create( );

void make_table(const col_vector& Sp);


public:

BlochAcq();
BlochAcq(const row_vector& det, const matrix& L, 
                                   const col_vector& Minf, double cut=1.e-12);
BlochAcq(const row_vector& det, const matrix& L, double cut=1.e-12);

BlochAcq(const BlochAcq& ACQ1);

~BlochAcq();

BlochAcq& operator= (const BlochAcq& ACQ1);

row_vector T(const col_vector& M0, int npts,         double tinc);
void       T(const col_vector& M0, row_vector& data, double tinc);

row_vector F(const MagVec& M0, int npts,         double Fst, double Ffi);
void       F(const MagVec& M0, row_vector& data, double Fst, double Ffi);

const TTable1D& table(const col_vector& M0);
const TTable1D& table() const;

///int ls() const;
///int size() const;

int full_size() const;


///int transitions() const;
///void parameters(matrix& mx, double& SW, double& LW,
///                                     double& dt, int N, int pf=0) const;

//       std::ostream& print(std::ostream& ostr,                         bool pf=false) const;
//       std::ostream& print(std::ostream& ostr, const col_vector& sigp, bool pf=false) const;
//friend std::ostream& operator << (std::ostream&    ostr, const BlochAcq &ACQ);
};
