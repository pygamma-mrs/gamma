/* Matrix.i */
// Swig interface file

%{
#include "Matrix/matrix.h"
%}

%include "std_vector.i"
%include "std_string.i"

%feature("autodoc", "1" );

%rename(__add__)  Matrix::operator+;
%rename(__iadd__) Matrix::operator+=;
%rename(__sub__)  Matrix::operator-;
%rename(__isub__) Matrix::operator-=;

%rename(__mul__)  Matrix::operator*;
%rename(__imul__) Matrix::operator*=;
%rename(__div__)  Matrix::operator/ const;
%rename(__idiv__) Matrix::operator/=;

%rename(__eq__)  Matrix::operator== const;
%rename(__ne__)  Matrix::operator!= const;
%rename(__lt__)  Matrix::operator<  const;
%rename(__gt__)  Matrix::operator>  const;


void enable_blockdiag();
void disable_blockdiag();

class matrix
{

public:

matrix( );
matrix(int i);

matrix(int i, int j);
matrix(int i, int j, matrix_type t);
matrix(int i, int j, matrix_type t, hermitian_type h);
matrix(int i, int j, const complex& z);
matrix(int i, int j, const complex& z, matrix_type t);
matrix(int i, int j, const complex& z, matrix_type t, hermitian_type h);

matrix(int i, int j, double d);
matrix(int i, int j, double d, matrix_type t);
matrix(int i, int j, double d, matrix_type t, hermitian_type h);

matrix(const matrix& mx);
~matrix ();

//matrix&  operator= (const matrix& mx);

//const complex& operator() (int i, int j) const;

complex& operator() (int i, int j);

complex& elem(int  i, int j);
complex  get(int   i, int j) const;
double   getRe(int i, int j) const;
double   getIm(int i, int j) const; 
void     put(const   complex& z, int i, int j);
void     put_h(const complex& z, int i, int j);

matrix   get_block(int row, int col, int nrows, int ncols);
void     put_block(int row, int col, const matrix& mx);


hermitian_type stored_hermitian( ) const;
hermitian_type check_hermitian(double       d = GMxCut);
void set_hermitian(hermitian_type h = _hermitian);
hermitian_type test_hermitian(double        d = GMxCut) const;

matrix_type stored_type( ) const;
matrix_type test_type(matrix_type t, double d=GMxCut) const;
void        set_type(matrix_type t);
matrix_type check_type(const matrix_type t, const double d = GMxCut);
std::string mxtype() const;
 
int cols( ) const;
int rows( ) const;
int refs( ) const;
int pts( )  const;
         
bool is_symmetric(const double d=GMxCut) const;
bool is_hermitian(const double d=GMxCut) const;
bool is_unitary(const   double d=GMxCut) const;
bool is_real(const      double d=GMxCut) const;
bool is_imaginary(const double d=GMxCut) const;
bool is_complex(const   double d=GMxCut) const;
bool is_zero(const      double d=GMxCut) const;
bool is_diagonal(const  double d=GMxCut) const;
bool is_square()                         const;
 
matrix operator + (const matrix& mx) const;
matrix operator - (const matrix& mx) const;
matrix operator * (const matrix& mx) const;
matrix operator * (const complex& z) const;
matrix operator * (      double   d) const;
matrix operator / (const matrix& mx) const;
matrix operator / (const complex& z) const;
matrix operator / (      double   d) const;

//friend matrix operator * (const complex& z,  const matrix& mx);
//friend matrix operator * (      double d,    const matrix& mx);
 
matrix& operator += (const matrix& mx1);
matrix& operator -= (const matrix& mx1);
matrix& operator *= (const matrix& mx);
matrix& operator *= (const complex& z);
matrix& operator *= (double d);
matrix& operator /= (const matrix& mx);
matrix& operator /= (const complex& z);
matrix& operator /= (double d);

/*
friend matrix  Re(const matrix& mx);
friend matrix  Im(const matrix& mx);
friend matrix  conj(const matrix& mx);
friend matrix  transpose(const matrix& mx);
friend matrix  adjoint(const matrix& mx);
friend complex trace(const matrix& mx);
friend void enable_blockdiag();
friend void disable_blockdiag();
*/

matrix  operator- () const;
matrix  Re()         const;
matrix  Im()         const;
matrix  conj()       const;   
matrix  transpose()  const; 
matrix  adjoint()    const;
matrix  exp()        const;
complex trace()      const;    

matrix  swaprows(int i, int j);
matrix  swapcols(int i, int j);
matrix  permute( int i, int j);

double  maxRe() const;
double  maxIm() const;
complex maxZ()  const;
double  minRe() const;
double  minIm() const;
complex minZ()  const;

//friend complex trace(const matrix& mx1, const matrix& mx2);
complex trace(const matrix& mx2) const;
//friend matrix adjoint_times(const matrix& mx, const matrix& mx1);
//friend matrix times_adjoint(const matrix& mx, const matrix& mx1);

//friend  complex         det(const matrix& mx);
complex         det() const;
//friend  int             rank(const matrix& mx);
//friend matrix FFT(const  matrix& mx);
//friend matrix IFFT(const matrix& mx);

matrix FFT()  const;
matrix IFFT() const;

//friend matrix tensor_product(const matrix& mx1, const matrix& mx2);

static void    Header(bool   hf);
static void    PrintRI(bool  pi);
static void    PrintAll(bool pa);
static void    PictDim(int   pd);
static void    PrintVal(bool pv);
static void    PrintCols(int cl);
static void    PrintRows(int rl);
//static MxPrint PrintFlags();

/*
std::ostream& printHdr(std::ostream& ostr) const;
std::ostream& print(std::ostream&    ostr) const;
std::ostream& picture(std::ostream&  ostr) const;
friend std::ostream&  operator << (std::ostream& ostr, const matrix& mx);
std::ofstream& write(std::ofstream& fp, int form=0) const;
friend std::ofstream& write(std::ofstream& F,  const matrix& mx);
std::ifstream& read(std::ifstream&  F);
*/

//friend std::istream&  operator >> (std::istream& istr, matrix& mx);
void           ask(const matrix_type t=n_matrix_type);

matrix resize(int i, int j);
matrix diagonal_form();

bool same_reference_as(const matrix& mx) const;

void status(int full=0) const;

std::vector<int> BlockDiag(matrix&    BD, std::vector<int> &U) const;
void        SymTriDiag(matrix&  HTD, matrix& U) const;
void        HermTriDiag(matrix& STD, matrix& U) const;
void        SymDiag(matrix&      SD, matrix& U) const;
void        Diagonalize(matrix&   D, matrix& U) const;

/*
friend void diag(const matrix& mx,     matrix&   D, matrix& U);
friend matrix inv(const matrix& mx);
friend matrix LU(matrix& mx, int* indx);
friend matrix LUinv(matrix& B, int *indx, matrix& LU);
*/

void                TestEigenSystem(int pf=1) const;
void                TestTransform(const matrix& T, const matrix& S, int pf=1) const;
std::vector<double> ColumnNorms() const;
std::vector<double> TestIdentity(complex& TotalDev) const;
matrix              TestUnitary(std::ostream& ostr) const;
matrix              TestUTransform(const matrix& T, const matrix& U) const;

bool operator== (const matrix& mx) const;
bool operator!= (const matrix& mx) const;
bool operator<  (const matrix& mx) const;
bool operator>  (const matrix& mx) const;

};


%extend matrix {

%pythonbegin %{
import numpy
%}

%pythoncode %{

def __str__(self):
    
    mmm = []
    rr = self.rows()
    cc = self.cols()
    
    if self.is_real():
        for r in range(rr):
            mmm.append([])
            for c in range(cc):
                ij = self.getRe(r,c)
                (mmm[r]).append( ij )
    else:
        for r in range(rr):
            mmm.append([])
            for c in range(cc):
                ij = self.get(r,c)
                (mmm[r]).append( ij.real()+ij.imag()*1j )
    
    mmm = numpy.array(mmm)
    mmm_str = mmm.__str__()
    del mmm
    return  mmm_str 
        

def toList( self ):
    
    mmm = []
    rr = self.rows()
    cc = self.cols()
    
    if self.is_real():    
        for r in range(rr):
            mmm.append([])
            for c in range(cc):
                ij = self.getRe(r,c)
                (mmm[r]).append( ij )
    else:   
        for r in range(rr):
            mmm.append([])
            for c in range(cc):
                ij = self.get(r,c)
                (mmm[r]).append( ij.real()+ij.imag()*1j )    
    return  mmm 


def __repr__( self ):
    
    mmm = []
    rr = self.rows()
    cc = self.cols()

    if self.is_real():
        for r in range(rr):
            mmm.append([])
            for c in range(cc):
                ij = self.getRe(r,c)
                (mmm[r]).append( ij )
    else:
        for r in range(rr):
            mmm.append([])
            for c in range(cc):
                ij = self.get(r,c)
                (mmm[r]).append( ij.real()+ij.imag()*1j )
    
    mmm = numpy.array(mmm)

    return mmm.__repr__()


def toNParray( self ):

    rr = self.rows()
    cc = self.cols()
    
    if self.is_real():
        mmm = numpy.zeros((rr,cc))
        for r in range(rr):
            for c in range(cc):
                ij = self.getRe(r,c)
                mmm[r][c] = ij
    else:
        mmm = numpy.zeros((rr,cc), dtype=numpy.complex128)
        for r in range(rr):
            for c in range(cc):
                ij = self.get(r,c)
                mmm[r][c] =  ij.real()+ij.imag()*1j
    return mmm


def __add__( self, val ):
    if isinstance( val, matrix ):
#        print "matrix:: add_matrix"
        return( self.add_matrix(val))

    if isinstance(val,gen_op):
#       print "matrix:: matix + genop"
       return( matrix_plus_genop(self,val))


def __sub__( self, val ):
    if isinstance( val, matrix ):
#        print "matrix:: sub_matrix"
        return( self.sub_matrix(val))

    if isinstance(val,gen_op):
#       print "matrix:: matrix - gen_op"
       return( matrix_minus_genop(self,val))

%}
};
