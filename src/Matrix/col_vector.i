/* col_vector.i */
// Swig Interface File


//%rename(colVector) col_vector;

%{
#include "Matrix/col_vector.h"
%}


%include "std_vector.i"
%include "std_string.i"

%feature("autodoc", "1" );
//%feature("director" );


%rename(__add__)  col_vector::operator+;
%rename(__iadd__) col_vector::operator+=;
%rename(__sub__)  col_vector::operator-;
%rename(__isub__) col_vector::operator-=;

%rename(__mul__)  col_vector::operator*;
%rename(__imul__) col_vector::operator*=;



//class col_vector;			 // Know that col_vector is a class

class col_vector : public matrix
{

//friend class col_vector;

public:


col_vector ( );
col_vector (int i);
col_vector (int i, double d);
col_vector (int i, const complex& z);
col_vector (const col_vector& cvec);

//col_vector (const matrix& mx);

~col_vector ();

//col_vector& operator = (const col_vector& cvec);
//col_vector& operator = (const col_vector& rvec);

//col_vector& operator = (const matrix& mx);

complex&    operator() (int i);
complex     get(int i) const;

double      getRe(int i) const;
double      getIm(int i) const;
void        put(const complex& z, int i);

int elements( ) const;
int size( )     const;

/*
friend col_vector operator * (const complex& z, const col_vector& cvec);
friend col_vector operator * (double         d, const col_vector& cvec);
friend col_vector operator * (const matrix& mx, const col_vector& cvec);
*/

col_vector operator + (const col_vector& cvec) const;
//col_vector operator + (const matrix&     mx)   const;

col_vector operator - (const col_vector& cvec) const;
//col_vector operator - (const matrix&     mx)   const;

//complex    operator * (const col_vector& cvec) const; // causes and error in C++ code.
//matrix     operator * (const col_vector& rvec) const;
//matrix     operator * (const matrix&     mx)   const;

col_vector operator * (const complex&    z)    const;
col_vector operator * (double            d)    const;

//col_vector operator / (const complex&    z)    const;
//col_vector operator / (double            d)    const;

col_vector& operator += (const col_vector& cvec1);
//void operator += (const matrix&        mx);

col_vector& operator -= (const col_vector& cvec1);
//void operator -= (const matrix&        mx);

col_vector& operator *= (const complex&        z);
col_vector& operator *= (      double          d);

col_vector & operator /= (const complex&        z);
col_vector & operator /= (const double          d);

// col_vector operator- ();                   // MATRIX INHERITED

/*
friend col_vector adjoint(const   col_vector& cvec);
friend col_vector transpose(const col_vector& cvec);
friend complex    trace(const     col_vector& cvec);
*/

col_vector adjoint()      const;
col_vector transpose()    const;
complex    trace()        const;
col_vector differential() const;

double           norm()          const;
complex          sum()           const;
double           maxRe()         const;
double           maxIm()         const;
complex          maxZ()          const;
double           minRe()         const;
double           minIm()         const;
complex          minZ()          const;
int              max(int type=0) const;
int              min(int type=0) const;
void             flip();
complex          sum(int st, int ne) const;
void             zero();

std::vector<int> sort(int type=0) const;

/*
friend col_vector FFT(const  col_vector& cvec);
friend col_vector IFFT(const col_vector& cvec);
friend col_vector product(const col_vector& cvec1, const col_vector& cvec2);
friend col_vector product(const col_vector& cvec,        col_vector& rvec);
*/

col_vector product()                        const;
col_vector product(const row_vector& rvec)  const;
col_vector product(const col_vector& cvec2) const;

/*
friend double  scalar_product(const col_vector& cvec);
friend complex scalar_product(const col_vector& cvec1, const col_vector& cvec2);
friend complex scalar_product(const col_vector& cvec,  const col_vector& rvec);
*/

double  scalar_product()                        const;

complex scalar_product(const col_vector& cvec2) const;
complex scalar_product(const row_vector& rvec)  const;

std::string hdrString() const;

//std::vector<std::string> printStrings(const MxPrint& PFlgs) const;

/*
std::ostream& printcols(std::ostream& ost, int nc=4,int rc=1,int ne=0) const;
std::ostream& print(std::ostream& ostr, int full=0) const;
friend std::ostream& operator << (std::ostream& ostr, const col_vector& cvec);
friend std::istream& operator >> (std::istream& istr, col_vector& cvec);
*/

void ask();

};


%extend col_vector {

%pythonbegin %{
import numpy
import pygamma
%}


%pythoncode %{



@classmethod
def from_list(class_object, col_vector_list ):
    co = class_object(len(col_vector_list))
    for i,v in enumerate(col_vector_list):
        co.put(pygamma.complex(v.real, v.imag),i)
    return co



		
def __str__(self):
    sss = "["
    
    if self.size() == 0:
        return( "[]" )
	
    
    for i in range( self.size()-1 ):
    
        rr = self.getRe(i)
        ii = self.getIm(i)

        sss = sss + "" + str(rr)
        if ii >= 0.0:
            sss = sss  + "+" + str(ii) + "j, "
        else:
            sss = sss + str(ii) + "j, "
	    
    i = self.size()-1

    rr = self.getRe(i)
    ii = self.getIm(i)

    sss = sss + " "+str(rr)
    if ii >= 0.0:
        sss = sss  + "+" + str(ii) + "j]"
    else:
        sss = sss + str(ii) + "j]"
   
    return(sss)




def __repr__(self):

    ll = self.size()

    if ll == 0:
        rrr = numpy.array([])
    else:

        if self.is_real():
            rrr = numpy.zeros(ll );
            for i in range(ll):
               rrr[i] = self.getRe(i)
        else:
            rrr = numpy.zeros(ll, dtype = numpy.complex128)
            for i in range(ll):
                rrr[i] = self.getRe(i)+1j*self.getIm(i)



    return "pygamma.col_vector({})".format(rrr.__repr__())   


def toNParray(self):

    ll = self.size()

    if ll == 0:
        rrr = numpy.array([])
    else:

        if self.is_real():
            rrr = numpy.zeros(ll);
            for i in range(ll):
               rrr[i] = self.getRe(i)
        else:
            rrr = numpy.zeros(ll, dtype = numpy.complex128)
            for i in range(ll):
                rrr[i] = self.getRe(i)+1j*self.getIm(i)

    return rrr


def __getitem__( self, sss ):
    """ Returns a new col_vector using the numpy slice notation
        by creating a temporary numpy array of the col_vector
    """
    
    
    if isinstance(sss,slice):
        nnn = self.toNParray()
        nnn = nnn[sss.start:sss.stop:sss.step]
        
        array_length = len(nnn)
        new_col_vector = col_vector(array_length)
        
        for i in range(array_length):
            new_col_vector.put( complex(nnn[i].real,nnn[i].imag), i)
            
        del nnn
        return( new_col_vector )
        
    elif isinstance( sss, int ):
        valr = self.getRe(sss)
        vali = self.getIm(sss)
        return( complex( valr, vali ))


def __setitem__( self, i, v ):
    """Sets the pygamma complex value, int or float v at position i in col vector
     
       rrr = pygamma.col_vector(10)
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
    
    

#def __getslice__( self, i,j ):
#    """Returns a copy of pygamma col_vector between i and j
#    
#       rrr = pygamma.col_vector
#       sss = rrr[i,j]
#    """
#
#    size = j-i
#    if size < 0:
#        size = -1*size
#
#    sss = col_vector( size, complex( 0 ))
#
#    for k in range(size):
#        sss.put( self.get(k+i), k )
#    return( sss )

        
            

def toList(self):
    mmm = []
    size = self.size()

    
    

    
    if self.is_real():    
        for i in range(size):
            ij = self.getRe(i)
            mmm.append( ij )
            
    else:
        for i in range(size):
            ij = self.get(i)
            mmm.append( ij.real()+ij.imag()*1j )
    
    return  mmm 


def Real(self):
    """Returns numpy array of Real part of pygamma col_vector"""
    size = self.size()
    sss = numpy.zeros( size, dtype=numpy.float32)
    for i in range(size):
        sss[i] = self.getRe(i) 
    return( sss )




def Imag(self):
    """Returns numpy array of imaginary part of pygamma col_vector"""
    size = self.size()
    sss =  numpy.zeros( size, dtype=numpy.float32)
    for i in range(size):
        sss[i] = self.getIm(i) 
    return( sss )





def fft_1D(self):
    """FFT of col vector returning a pygamma col_vector"""
    spec = self.FFT()
    size = spec.pts()

    sss = col_vector( size, complex( 0 ))
    for i in range(size):
        sss.put( spec.get(0,i), i )
    return( sss )

def __rmul__( self,other ):
	
	return( self * other )

      

%}
};

