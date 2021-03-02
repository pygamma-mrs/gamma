// row_vector.i
// Swig interface file.

//%rename(rowVector) row_vector;

%{
#include "Matrix/row_vector.h"
%}

%include "std_vector.i"
%include "std_string.i"

%feature("autodoc", "2" );


%rename(__add__)  row_vector::operator+;
%rename(__iadd__) row_vector::operator+=;
%rename(__sub__)  row_vector::operator-;
%rename(__isub__) row_vector::operator-=;

%rename(__mul__)  row_vector::operator*;
%rename(__imul__) row_vector::operator*=;


%rename(pyG_FFT)  FFT(const  row_vector&);
%rename(pyG_IFFT) IFFT(const row_vector&);

%rename(complex_times_rowvector) operator* (const complex&, const row_vector&);
%rename(double_times_rowvector)  operator* (double,         const row_vector&);
//class col_vector;

class row_vector : public matrix
{

//friend class col_vector;


public:


row_vector();
row_vector(int i);
row_vector(int i, const complex& z);
row_vector(int i, double d);
row_vector(const row_vector& rvec);
//row_vector(const matrix& mx);

~row_vector();

//void     operator = (const row_vector& rvec);

//void     operator = (const col_vector& cvec);
//void     operator = (const matrix& mx);

complex& operator() (int i);
complex  get(int i)   const;

double   getRe(int i) const;
double   getIm(int i) const;
void     put(const complex& z, int i);

int elements( ) const;
int size( )     const;

friend row_vector operator * (const complex& z, const row_vector& rvec);
friend row_vector operator * (double         d, const row_vector& rvec);

row_vector operator + (const row_vector& rvec) const;
//row_vector operator + (const matrix&     mx)   const;

row_vector operator - (const row_vector& rvec) const;
//row_vector operator - (const matrix&     mx)   const;

//complex    operator * (const row_vector& rvec) const; // this causes a c++ error.
complex    operator * (const col_vector& cvec) const;

//row_vector operator * (const matrix& mx)       const;

row_vector operator * (const complex& z)       const;
row_vector operator * (double d)               const;
//row_vector operator / (const complex& z)       const;
//row_vector operator / (double d)               const;

row_vector& operator += (const row_vector& rvec1);
//void operator += (const matrix& mx);

row_vector& operator -= (const row_vector& rvec1);
//void operator -= (const matrix& mx);

row_vector& operator *= (      double d);
row_vector& operator *= (const complex& z);
row_vector& operator /= (      double d);
row_vector& operator /= (const complex& z);

/*
friend col_vector adjoint(const     row_vector& rvec);
friend col_vector transpose(const   row_vector& rvec);
friend complex    trace(const       row_vector& rvec);
*/

row_vector differential() const;

double           norm()              const;
complex          sum()               const;
double           maxRe()             const;
double           maxIm()             const;
complex          maxZ()              const;
double           minRe()             const;
double           minIm()             const;
complex          minZ()              const;
int              max(int type=0)     const;
int              min(int type=0)     const;
void             flip();
complex          sum(int st, int ne) const;
void             zero();

std::vector<int> sort(int type=0)    const;

friend row_vector FFT(const  row_vector& rvec);
friend row_vector IFFT(const row_vector& rvec);
/*
friend row_vector FFT(const  row_vector& rvec);
friend row_vector IFFT(const row_vector& rvec);

friend int operator==(const row_vector& rvec1, const row_vector& rvec2);
friend int operator!=(const row_vector& rvec1, const row_vector& rvec2);

friend row_vector product(const row_vector& rvec1, const row_vector& rvec2);
friend row_vector product(const row_vector& rvec,        col_vector& cvec);
*/

row_vector product()                        const;
row_vector product(const row_vector& cvec2) const;
row_vector product(const col_vector& rvec)  const;

/*
friend double  scalar_product(const row_vector& rvec);
friend complex scalar_product(const row_vector& rvec1, const row_vector& rvec2);

friend complex scalar_product(const row_vector& rvec,  const col_vector& cvec);
*/

row_vector  FFT()  const;
row_vector IFFT() const;

double  scalar_product()                        const;
complex scalar_product(const row_vector& cvec2) const;
complex scalar_product(const col_vector& rvec)  const;

std::string              hdrString() const;

//std::vector<std::string> printStrings(const MxPrint& PFlgs) const;

/*
       std::ostream& printcols(std::ostream& ost, int nc=4,int rc=1,int ne=0) const;
       std::ostream& print(std::ostream& ostr, int full=0) const;
friend std::ostream& operator << (std::ostream& ostr, const row_vector& rvec);
friend std::istream& operator >> (std::istream& istr, row_vector& rvec);
*/

enum { SMIS = 0, SIEMENS, SIEMENS_NOHDR, PLAIN_ASCII, ASCII_MT_DEG, SVS };

static row_vector read_pulse (const std::string filename, const int PulseFmt);

void ask();

};

%extend row_vector {
    
%pythonbegin %{
import numpy
import math
%}

%pythoncode %{

@classmethod
def from_list(class_object, row_vector_list ):
    co = class_object(len(row_vector_list))
    for i,v in enumerate(row_vector_list):
        co.put(pygamma.complex(v.real, v.imag),i)
    return co


def __str__(self):

    ll = self.size()
    if ll == 0:
        rrr = numpy.array([])
    else:
        if not self.is_real():
            rrr = numpy.zeros(ll, dtype = numpy.complex128)
            for i in range(ll):
                rrr[i] = self.getRe(i)+1j*self.getIm(i)
        else:
            rrr = numpy.zeros(ll );

            for i in range(ll):
               rrr[i] = self.getRe(i)

    return rrr.__str__()



def __repr__(self):

    ll = self.size()
    if ll == 0:
        rrr = numpy.array([])
    else:
        if not self.is_real():
            rrr = numpy.zeros(ll, dtype = numpy.complex128)
            for i in range(ll):
                rrr[i] = self.getRe(i)+1j*self.getIm(i)
        else:
            rrr = numpy.zeros(ll );

            for i in range(ll):
               rrr[i] = self.getRe(i)

    return rrr.__repr__()


def toNParray(self):

    ll = self.size()

    if ll == 0:
        rrr = numpy.array([])
    else:
        if not self.is_real():
            rrr = numpy.zeros(ll, dtype = numpy.complex128)
            for i in range(ll):
                rrr[i] = self.getRe(i)+1j*self.getIm(i)
        else:
            rrr = numpy.zeros(ll);

            for i in range(ll):
               rrr[i] = self.getRe(i)
    return rrr


def __getitem__( self, sss ):
    """ Returns a new row_vector using the numpy slice notation
        by creating a temporary numpy array of the row_vector
    """
    
    if isinstance(sss,slice):
        nnn = self.toNParray()
        nnn = nnn[sss.start:sss.stop:sss.step]
        
        array_length = len(nnn)
        new_row_vector = row_vector(array_length)
        
        for i in range(array_length):
            new_row_vector.put( complex(nnn[i].real,nnn[i].imag), i)
            
        del nnn
        return( new_row_vector )
        
    elif isinstance( sss, int ):
        valr = self.getRe(sss)
        vali = self.getIm(sss)
        return( complex( valr, vali ))


def __setitem__( self, i, v ):
    """Sets the pygamma complex value, int or float v at position i in row vector
     
       rrr = pygamma.row_vector(10)
       rrr[3] = pygamma.complex( 7, 5 )
       rrr[2] = 1
       rrr[1] = 9.0
    
    """
    if isinstance(v,float) or isinstance(v, int ):
        self.put( complex(v,0), i )
    elif isinstance(v,complex):
        self.put( v, i )


def __len__(self):
    return( self.size())
    

def toList(self):
    mmm = []
    size = self.size()
  
    if not self.is_real():
        for i in range(size):
            ij = self.get(i)
            mmm.append( ij.real()+ij.imag()*1j )
    else:
        for i in range(size):
            ij = self.getRe(i)
            mmm.append( ij )
    
    return  mmm 


def Real(self):
    """Returns numpy array of Real part of pygamma row_vector"""
    size = self.size()
    sss = numpy.zeros( size, dtype=numpy.float32)
    for i in range(size):
        sss[i] = self.getRe(i) 
    return( sss )


def Imag(self):
    """Returns numpy array of imaginary part of pygamma row_vector"""
    size = self.size()
    sss =  numpy.zeros( size, dtype=numpy.float32)
    for i in range(size):
        sss[i] = self.getIm(i) 
    return( sss )


def fft_1D(self):
    """FFT of row vector returning a pygamma row_vector"""
    spec = self.FFT()
    size = spec.pts()

    sss = row_vector( size, complex( 0 ))
    for i in range(size):
        sss.put( spec.get(0,i), i )
    return( sss )




    
def __rmul__(self, val ):

    if isinstance( val, int ) or isinstance( val, float ):
        return( double_times_rowvector( val, self))

    if isinstance( val, complex ):
        return( complex_times_rowvector( val, self))



     

%}
};
