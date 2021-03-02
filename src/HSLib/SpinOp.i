/* SpinOp.i */
// Swig interface file.

%{
#include "HSLib/SpinOp.h"
%}

%include "std_string.i"
%include "GenOp.i"

%feature("autodoc", "1" );

%rename(__add__)  spin_op::operator+;
%rename(__iadd__) spin_op::operator+=;
%rename(__isub__) spin_op::operator-=;
%rename(__neg__)  spin_op::operator- ();


%rename(spinop_times_double)  operator*  (const spin_op&,  double);
%rename(double_times_spinop)  operator*  (double,          const spin_op& );
%rename(spinop_times_complex) operator*  (const spin_op&,  const complex&);
%rename(complex_times_spinop) operator*  (const complex&,  const spin_op& );

%rename(spinop_times_spinop ) operator*  (const spin_op&, const spin_op&);

%rename(spinop_divide_complex) operator/  (const spin_op&, const complex&);
%rename(spinop_divide_double)  operator/  (const spin_op&, double);


%rename(__imul__) spin_op::operator*=;
%rename(__idiv__) spin_op::operator/=;

%rename(__assign__) spin_op::operator=;

%rename(Trace) trace(const spin_op&);

class spin_op 
{

public:

spin_op();
//spin_op(int spins, matrix* prmxs);
spin_op(const spin_op& SOp);

~spin_op();

spin_op& operator= (const spin_op& SOp);

spin_op operator-  ()                    const;
spin_op operator+  (const spin_op& SOp1) const;

spin_op & operator-= (const spin_op& SOp1);
spin_op & operator+= (const spin_op& SOp1);

//friend spin_op          operator-  (const spin_op& SOp1, const spin_op& SOp2);
friend spin_op          operator*  (const spin_op& SOp1, const spin_op& SOp2);

spin_op & operator*= (const spin_op& SOp1);

friend spin_op       operator*  (const spin_op& SOp, const complex& z);
friend spin_op       operator*  (const complex& z,   const spin_op& SOp);

friend spin_op       operator*  (const spin_op& SOp, double d);
friend spin_op       operator*  (double d,           const spin_op& SOp);



spin_op & operator*= (const complex& z);
spin_op & operator*= (double d);

friend spin_op       operator/  (const spin_op& SOp, const complex& z);
friend spin_op       operator/  (const spin_op& SOp, double d);

spin_op & operator/= (const complex& z);
spin_op & operator/= (double d);

matrix   get_mx() const;
//operator matrix() const;

spin_op exp()     const; 		// Exponentation
spin_op adjoint() const; 		// Adjoint
complex trace()   const; 		// Trace

//friend spin_op exp(const spin_op& SOp); 		// Exponentation
//friend spin_op adjoint(const spin_op &SOp); 		// Adjoint
friend complex trace(const spin_op &SOp); 		// Trace

int spins( )    const;				// Number of spins
int refs( )     const;				// Full space mx refs
int refs(int i) const;				// Sub-space mx refs
int HS( )       const;				// Full Hilbert space

//void       print(std::ostream& ostr, int full=0) const;
//friend std::ostream& operator<< (std::ostream& ostr, const spin_op& SOp);

void status(int full=1) const;

void FaxisStruct(char axis) const;

};

%extend spin_op {

%pythonbegin %{
import numpy
%}

%pythoncode %{

def __str__(self):

    mmm = self.get_mx()
    return mmm.__str__()
 

def assign( self, mmm ):
    pass

def toList( self ):
    mmm = self.get_mx()
    return mmm.toList()
    


def __repr__( self ):
    mmm = self.get_mx()
    return mmm.__repr__()



def toNParray( self ):
    mmm = self.get_mx()
    return mmm.toNParray()



def __mul__(self, val ):

    if isinstance( val, int ) or isinstance( val, float ):
        return( spinop_times_double(  self,val))

    if isinstance( val, complex ):
        return( spinop_times_complex( self,val))

    if isinstance( val, spin_op ):
        return( spinop_times_spinop( self, val ))

    if isinstance( val, gen_op ):
        return( val.times_genop(gen_op(self)))





def __rmul__(self, val ):

    if isinstance( val, int ) or isinstance( val, float ):
        return( double_times_spinop( val, self))

    if isinstance( val, complex ):
        return( complex_times_spinop( val, self))

    if isinstance( val, spin_op ):
        return( spinop_times_spinop( val, self ))



def __div__( self, val ):

    if isinstance( val, complex ):
        return( spinop_divide_complex( self, val ))


    if isinstance( val, float ):
        return( spinop_divide_double( self, val ))

%}
};


