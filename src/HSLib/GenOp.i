// GenOp.i
// Swig interface file.

%{
#include "HSLib/GenOp.h"
%}

%include "std_vector.i"
%include "std_string.i"
%template(genoprep_vec) std::vector <genoprep>;

%feature("autodoc", "1" );


%rename(__iadd__) gen_op::operator+=;

%rename(__neg__)  gen_op::operator-() const;
%rename(__isub__) gen_op::operator-=;

//%rename(__mul__)  gen_op::operator* const;
%rename(__imul__) gen_op::operator*=;
%rename(__idiv__) gen_op::operator/=;

%rename(__iand__) gen_op::operator&=;

%rename(__eq__)  gen_op::operator== const;
%rename(__ne__)  gen_op::operator!= const;
%rename(__lt__)  gen_op::operator<  const;
%rename(__gt__)  gen_op::operator>  const;

//%rename(__add__)  gen_op::operator+ const;
%rename(add_genop)  gen_op::operator+ const;

//%rename(__sub__)  gen_op::operator- const;
%rename(sub_genop)  gen_op::operator- const;

//%rename(__mul__)  gen_op::operator* const;
%rename(times_genop)  gen_op::operator* const;


%rename(genop_minus_spinop)  operator- (const gen_op&, const spin_op& );
%rename(genop_minus_matrix)  operator- (const gen_op&, const matrix&  );
%rename(matrix_minus_genop)  operator- (const matrix&, const gen_op&  );

%rename(genop_plus_spinop) operator+ (const gen_op&, const spin_op& );
%rename(genop_plus_matrix ) operator+ (const gen_op&, const matrix&  );
%rename(matrix_plus_genop ) operator+ (const matrix&, const gen_op& );

%rename(genop_times_complex) operator * (const gen_op&, const complex&);
%rename(complex_times_genop) operator * (const complex&,  const gen_op&);

%rename(genop_times_double) operator * (const gen_op&, double );
%rename(double_times_genop) operator * (double, const gen_op&);

%rename(genop_times_matrix) operator * (const gen_op&, const matrix&);
%rename(matrix_times_genop)  operator * (const matrix&,  const gen_op&);

%rename(genop_divide_complex)operator / (const gen_op&, const complex&);
%rename(genop_divide_float)operator / (const gen_op&, double);

%rename(Trace)  trace(const gen_op&);
%rename(Det)    det(const gen_op&);

typedef std::vector <genoprep> genoprep_vec;

class gen_op : public genoprep_vec
{

public:

gen_op();
gen_op(const matrix&  mx);

gen_op(const spin_op& SOp);
gen_op(const matrix&  mx1, const matrix& mx2);
gen_op(const matrix&  mx, const basis&  bs);

gen_op(const gen_op&  Op1);
//gen_op(const std::vector<matrix>& mxc, const std::vector<matrix>& bsc);


// sosi - the constructor below should be completely replaced by the above one
//gen_op(matrix* mxc,   int nc,      matrix* bsc=NULL);
//gen_op(int N,         matrix* mxs, matrix* bss);

virtual ~gen_op( );			// Floquet needs this virtual



gen_op operator+  (const gen_op& Op)  const;
gen_op &   operator+= (const gen_op& Op);
gen_op operator-  (const gen_op& Op)  const;
gen_op operator-  ()                  const;
gen_op &   operator-= (const gen_op& Op);
gen_op operator*  (const gen_op& Op)  const;
gen_op &   operator*= (const gen_op &Op);
gen_op &   operator&= (const gen_op &Op);

//friend gen_op operator- (const gen_op& Op);

//void operator = (const gen_op &Op1);

// sosi - added for alpha until spin op --> tp_matrix

//void   operator = (const spin_op &SOp);

gen_op &   operator += (const spin_op &SOp);
gen_op &   operator -= (const spin_op &SOp);

friend gen_op operator + (const gen_op& Op1, const spin_op& SOp);
friend gen_op operator - (const gen_op &Op1, const spin_op &SOp);

friend gen_op operator + (const gen_op& Op1, const matrix& mx);
friend gen_op operator - (const gen_op& Op1, const matrix& mx);

friend gen_op operator + (const matrix& mx,  const gen_op& Op1);
friend gen_op operator - (const matrix& mx,  const gen_op& Op1);

friend gen_op operator * (const gen_op& Op1, const matrix& mx);
friend gen_op operator * (const matrix& mx,  const gen_op& Op1);

/*
*/

/*
void   operator += (const matrix& mx);	
void   operator -= (const matrix &mx);
void   operator *= (const matrix& mx);
void   operator &= (const matrix& mx);
void   operator  = (const matrix& mx);
*/

friend gen_op operator * (const gen_op &Op1, const complex &z);
friend gen_op operator * (const complex &z,  const gen_op &Op1);
friend gen_op operator * (const gen_op &Op1, double z);
friend gen_op operator * (double z,          const gen_op &Op1);



gen_op & operator *= (const complex &z);
gen_op & operator *= (double r);

friend gen_op operator / (const gen_op &Op1, const complex &z);
friend gen_op operator / (const gen_op& Op1, double r);

gen_op & operator /= (const complex& z);
gen_op & operator /= (double r);

complex det() const;
complex trace() const;
complex trace(const gen_op& Op2) const;
complex proj(const gen_op &Op2, int norm=0) const;

int dim()    const;		// Hilbert space
int HS()     const;		// Hilbert space
int LS()     const;		// Liouville space
int dim_LS() const;		// Liouville space
  
gen_op exp() const;

gen_op exp(const complex& t, double cutoff=1.e-12) const;

gen_op Pow(int power) const;

/*
friend gen_op tensor_product(const gen_op &Op1, const gen_op &Op2);
friend gen_op log(const gen_op &Op1);
*/

gen_op sim_trans(const gen_op &Op2) const;
  
void sim_trans_ip(const gen_op &Op1);

// gen_op adjoint();  // Disabled!!

row_vector eigvals() const;
void eigvals(double* vx, int sort=0) const;

friend complex trace(const gen_op &Op);
friend complex det(const gen_op &Op);

//***************************************************************************
//****************************************************************************
// Marked for deletion....
/*
friend complex proj(const gen_op &Op1, const gen_op &Op2, int norm=0);  //**
friend int dim(const gen_op &Op);  					//**
friend gen_op exp(const gen_op &Op1);					//**
friend gen_op pow(const gen_op &Op1, int power);			//**
friend gen_op sim_trans(const gen_op &Op1, const gen_op &Op2);		//**
friend gen_op adjoint(const gen_op &Op1);				//**
*/
//***************************************************************************
//***************************************************************************

//matrix get_mx( )     const;
//matrix get_matrix( ) const;
//void   put_mx(const matrix &mx);
//void   put_matrix(const matrix &mx);

basis get_bs( )    const;
basis get_basis( ) const;
void  put_bs(const basis &bs);
void  put_basis(const basis &bs);


complex operator() (int row, int col) const;
complex get(int row, int col)         const;
void    put(const complex &z, int row, int col);

std::string name() const;
void name(const   std::string& n);
void bsname(const std::string& bn);

int exists( ) const;

col_vector superket() const;
void desuperket(const col_vector& mS);

gen_op project_sub(int ic) const;

int test_EBR() const;
int test_DBR() const;
int in_EBR()   const;
int in_DBR()   const;
	
void set_DBR() const;
void set_EBR() const;
void Op_base(const gen_op &Op1, double cutoff=1.e-12) const; 
void Op_base(const basis  &bs) const;

void status(int pf=0) const;

void setOnlyWBR( );

void Op_priority (int pty);

void SetLimits(int limit) const;
  
int  OpCheck(const gen_op& Op1, int warn=2) const;
//bool OpCheck(const matrix& mx, int warn=2) const;

int  OpCheck(const basis&  bs, int warn=2) const;
//int  OpCheck(const matrix& mx1, const matrix& mx2, int warn=2) const;
//int  OpCheck(const matrix& mx, const basis&  bs, int warn=2) const;

int OpCheck(int row,    int col, int warn=2) const;
int LimCheck(int limit,          int warn=1) const;

//std::ostream& print(std::ostream& out, int full=0) const;
//friend std::ostream& operator << (std::ostream& ostr, const gen_op &Op);
//std::ostream& eigenvalues(std::ostream& ostr, int rc=0, int ncol=4) const;

void write(const std::string& fn) const;
//std::ofstream& write(std::ofstream& fp) const;

//void read(const std::string& fn, gen_op& Op1); NOW DEPRECATED!

void read(const std::string& fn, const basis& bs);
void read(const std::string& fn);
//std::ifstream& read(std::ifstream &fp);

double TestEigenSystem(int pf=0);

bool is_symmetric(const double d=GMxCut) const;
bool is_hermitian(const double d=GMxCut) const;
bool is_unitary(const   double d=GMxCut) const;
bool is_real(const      double d=GMxCut) const;
bool is_imaginary(const double d=GMxCut) const;
bool is_complex(const   double d=GMxCut) const;
bool is_zero(const      double d=GMxCut) const;
bool is_diagonal(const  double d=GMxCut) const;
bool is_square()                         const;


/* Aside for providing basic tests as to whether to operators are equvalent or
   not, these operators are necessary if any STL container classes are to
   be used based on operators (e.g. list<gen_op> or vector<gen_op>)          */  

bool operator== (const gen_op& Op) const;
bool operator!= (const gen_op& Op) const;
bool operator<  (const gen_op& Op) const;
bool operator>  (const gen_op& Op) const;
 
};


inline gen_op I_gen_op(const basis& bs)
// Input		bs   : Basis
// Return		Op   : Identity Operator with basis bs 
{ 
    return gen_op(matrix(bs.size(),bs.size(),i_matrix_type),bs); 
};


%extend gen_op {

%pythonbegin %{
import numpy
%}

%pythoncode %{

def __str__(self):
    

    rr = self.dim()
    cc = self.dim()
    
    if not self.is_real():
        mmm = numpy.zeros((rr,cc), dtype=numpy.complex128)
        for r in range(rr):
            for c in range(cc):
                ij = self.get(r,c)
                mmm[r][c] =  ij.real()+ij.imag()*1j
    else:
        mmm = numpy.zeros((rr,cc))
        for r in range(rr):
            for c in range(cc):
                ij = self.get(r,c)
                mmm[r][c] = ij.real()

    return mmm.__str__()


def assign( self, mmm ):
    pass


def toList( self ):
    mmm = []
    rr = self.dim()
    cc = self.dim()

    if not self.is_real():
        for r in range(rr):
            mmm.append([])
            for c in range(cc):
                ij = self.get(r,c)
                (mmm[r]).append( ij.real()+ij.imag()*1j )
    else:
        for r in range(rr):
            mmm.append([])
            for c in range(cc):
                ij = self.get(r,c)
                (mmm[r]).append( ij.real() )

    return  mmm


def __repr__( self ):

    mmm = []
    rr = self.dim()
    cc = self.dim()

    if not self.is_real():
        mmm = numpy.zeros((rr,cc), dtype=numpy.complex128)
        for r in range(rr):
            for c in range(cc):
                ij = self.get(r,c)
                mmm[r][c] =  ij.real()+ij.imag()*1j
    else:
        mmm = numpy.zeros((rr,cc))
        for r in range(rr):
            for c in range(cc):
                ij = self.get(r,c)
                mmm[r][c] = ij.real()

    return mmm.__repr__()


def toNParray( self ):

    rr = self.dim()
    cc = self.dim()

    if not self.is_real():
        mmm = numpy.zeros((rr,cc), dtype=numpy.complex128)
        for r in range(rr):
            for c in range(cc):
                ij = self.get(r,c)
                mmm[r][c] =  ij.real()+ij.imag()*1j
    else:
        mmm = numpy.zeros((rr,cc))
        for r in range(rr):
            for c in range(cc):
                ij = self.get(r,c)
                mmm[r][c] = ij.real()
    return mmm


def __mul__( self, val ):
    if isinstance( val, gen_op ):
        #print "gen_op:: times_genop"
        return( self.times_genop(val))

    if  isinstance( val, spin_op ):
        #print "gen_op:: times_genop"
        return( self.times_genop(gen_op(val)))

    if isinstance( val, complex ):
        #print "gen_op:: genp * complex"
        return( genop_times_complex( self, val ) )

    if isinstance( val, float ):
        #print "gen_op:: genOp * double"
        return( genop_times_double( self, val ) )

    if isinstance( val, int ):
        #print "gen_op:: genOp * double"
        return( genop_times_double( self, val ) )

    if isinstance( val, matrix ):
        #print "gen_op:: genop * matrix"
        return( genop_times_matrix( self, val ) )


def __rmul__( self, val ):

    if  isinstance( val, spin_op ):
        #print "rmul:: times_genop"
        return( self.times_genop(gen_op(val)))

    if isinstance( val, float ):
        #print "gen_op:: double * genop"
        return( double_times_genop(  val,self ) )

    if isinstance( val, int ):
        #print "gen_op:: double * genop"
        return( double_times_genop(  val,self ) )

    if isinstance( val, complex ):
        #print "gen_op:: complex * genop"
        return( complex_times_genop(  val,self ) )

    if isinstance( val, matrix ):
        #print "gen_op:: matrix * genop"
        return( matrix_times_genop( val, self ) )

def __add__( self, val ):
    if isinstance( val, gen_op ):
        #print "gen_op:: add_genop"
        return( self.add_genop(val))

    if isinstance( val, spin_op ):
        #print "gen_op:: genop + spinop"
        return( genop_plus_spinop( self, val ))

    if isinstance( val, matrix ):
        #print "gen_op:: genop + matrix"
        return( genop_plus_matrix( self, val ))

def __radd__( self, val ):

    if isinstance( val, matrix ):
        #print "gen_op:: matrix + genop"
        return( matrix_plus_genop( val,self ))

def __sub__( self, val ):
    if isinstance( val, gen_op ):
        #print "gen_op::  - genop"
        return( self.sub_genop(val))

    if isinstance( val, spin_op ):
        #print "gen_op:: genop - spinop"
        return( genop_minus_spinop(self,val))

    if isinstance( val, matrix ):
        #print "gen_op:: genop - matrix"
        return( genop_minus_matrix(self,val))

def __rsub__( self, val ):

    if isinstance( val, matrix ):
        #print "gen_op:: matrix - gen_op"
        return( matrix_minus_genop( val,self ))


def __div__( self, val ):

    if isinstance( val, complex ):
        #print "gen_op:: gen_op / complex"
        return( genop_divide_complex( self, val ))

    if isinstance( val, float ):
        #print "gen_op:: gen_op / float"
        return( genop_divide_float( self, val ))

%}
};






