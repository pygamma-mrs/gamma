/* complex.i */
// Swig interface file


%{
#include "Matrix/complex.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );
%feature("director" );

%rename(__add__)  complex::operator+;
%rename(__iadd__) complex::operator+=;
%rename(__sub__)  complex::operator-;
%rename(__isub__) complex::operator-=;
%rename(__pos__)  complex::operator+ ();
%rename(__neg__)  complex::operator- ();

%rename(__mul__)  complex::operator* const;
%rename(__imul__) complex::operator*=;
%rename(__div__)  complex::operator/ const;
%rename(__idiv__) complex::operator/=;
//%rename(__iand__) complex::operator&=;

%rename(__eq__)  complex::operator== const;
%rename(__ne__)  complex::operator!= const;
%rename(__lt__)  complex::operator<  const;
%rename(__gt__)  complex::operator>  const;


class complex
{

public:

friend inline void Swap(complex& z1, complex& z2);
              
complex();
complex(double r,   double i=0.0);
complex(const complex& z);

//complex& operator= (const complex& z);
//complex& operator= (double r);

inline        double& Relem();
inline        double& Ielem();
inline const  double& Rec();
inline const  double& Imc();

// For PyGAMMA
inline        double real();
inline        double imag();

/*
friend inline double& zRe(complex& z);
friend inline double& zIm(complex& z);

friend inline double  Re(const complex& z);
friend inline double  Im(const complex& z);

friend inline void set_real_part(complex& z, double r);
friend inline void Re(complex& z, double r);

friend inline void set_imaginary_part(complex& z, double r);
friend inline void Im(complex& z, double r);
*/
	
inline complex  operator+  ()                  const;
inline complex  operator+  (const complex& z)  const;
inline complex  operator+  (double r)          const;
inline complex& operator+= (const complex& z);
inline complex& operator+= (double r);
inline complex  operator-  ()                  const;
inline complex  operator-  (double r);
inline complex  operator-  (const complex& z)  const;
inline complex& operator-= (const complex& z);
inline complex& operator-= (double r);

inline complex  operator*  (const complex& z)  const;
inline complex  operator*  (double r)          const;
inline complex& operator*= (const complex& z);
inline complex& operator*= (double r);
inline complex  operator/  (const complex& z)  const;
inline complex  operator/  (double r)          const;
inline complex& operator/= (const complex& z);
inline complex& operator/= (double r);


/*
friend inline complex operator+ (double r, const complex& z);
friend inline complex operator- (double r, const complex& z);
friend inline complex operator* (double r, const complex& z);
friend inline complex operator/ (double r, const complex& z);

friend inline void add(complex& z, const complex& z1, const complex& z2);
friend inline void add(complex& z,       double    r, const complex& z1);
friend inline void add(complex& z, const complex& z1,       double    r);
friend inline void sub(complex& z, const complex& z1, const complex& z2);
friend inline void sub(complex& z,       double    r, const complex& z1);
friend inline void sub(complex& z, const complex& z1,       double    r);
friend inline void mul(complex& z, const complex& z1, const complex& z2);
friend inline void mul(complex& z,       double    r, const complex& z1);
friend inline void mul(complex& z, const complex& z1,       double    r);
friend inline void div(complex& z, const complex& z1, const complex& z2);
friend inline void div(complex& z,       double    r, const complex& z1);
friend inline void div(complex& z, const complex& z1,       double    r);
*/

complex conj() const;
complex conj_times(const complex& z) const;


/*
friend inline complex conj(const complex& z);
friend inline complex conj(const complex& z, const complex& z1);

friend inline double AbsNorm(const     complex& z);
friend inline double square_norm(const complex& z);
friend inline double norm(const        complex& z);
friend inline double phase(const       complex& z);
friend inline void   norm(complex& z,  double r);
friend inline void   phase(complex& z, double r);

friend inline complex sqrt(const complex& z);
*/

inline complex Zexp() const;

/*
friend inline complex exp(const complex& z);

friend inline complex log(const complex& z);

friend inline complex pow(const complex& z, const complex& z1);
*/

complex Zsin()   const;
complex Zcos()   const;
complex Ztan()   const;
complex Zasin()  const;
complex Zacos()  const;
complex Zatan()  const;
complex Zsinh()  const;
complex Zcosh()  const;
complex Ztanh()  const;
complex Zasinh() const;
complex Zacosh() const;
complex Zatanh() const;

/*
friend complex sin(const   complex& z);
friend complex cos(const   complex& z);
friend complex tan(const   complex& z);
friend complex asin(const  complex& z);
friend complex acos(const  complex& z);
friend complex atan(const  complex& z);
friend complex sinh(const  complex& z);
friend complex cosh(const  complex& z);
friend complex tanh(const  complex& z);
friend complex asinh(const complex& z);
friend complex acosh(const complex& z);
friend complex atanh(const complex& z);
*/

static void Reiphi(bool TF);
static void delim(bool  TF);
static void scinot(bool TF);
static void digits(int  digs);
static void dadp(int    adp);

static int         dlength();			// Length of printed d
static std::string dformat();			// Format of printed d
static int         zlength();			// Length of printed z
static bool        normphase();		// True if R,phi out

static void complex_getf(bool& phase, bool& math, bool& science,
			                  int &digits, int &digs_aft_dpoint );


//friend void complex_getf(std::string& form);

std::string   printString() const;
//       std::ostream& print(std::ostream& ostr) const;

//friend std::ostream& operator<< (std::ostream& ostr, const complex& z);
//friend std::istream& operator>> (std::istream& istr, complex& z);


void write(const std::string& fn);
//void write(std::ofstream& fp) const;

void read(const std::string& fn);
//void read(std::ifstream& fp);

inline bool operator== (const complex& z) const;
inline bool operator!= (const complex& z) const;
       bool operator<  (const complex& z) const;
       bool operator>  (const complex& z) const;

};


extern const complex complex0;		// z = 0 : (0,0)
extern const complex complex1;		// z = 1 : (1,0)
extern const complex complexi;		// z = i : (0,1)

/*
inline complex::complex() { }
inline complex::complex(double r, double i) :re(r),im(i) { }
inline complex::complex(const complex& z) :re(z.re),im(z.im) { }

inline complex& complex::operator= (const complex &z)
{ re = z.re; im = z.im; return *this; }

inline complex& complex::operator= (double r)
{ re = r; im = 0;  return *this; }

inline void Swap(complex& z1, complex& z2)
{
  double tmp;
  tmp = z1.re; z1.re = z2.re; z2.re = tmp;
  tmp = z1.im; z1.im = z2.im; z2.im = tmp;
}


inline        double& complex::Relem()     { return re; }
inline        double& complex::Ielem()     { return im; }
inline const  double& complex::Rec()       { return re; }
inline const  double& complex::Imc()       { return im; }
inline        double& zRe(complex& z)      { return z.re; }
inline        double& zIm(complex& z)      { return z.im; }

inline        double  Re(const complex& z) { return z.re; }
inline        double  Im(const complex& z) { return z.im; }


inline void Re(complex& z,                 double r) { z.re = r; }
inline void Im(complex& z,                 double r) { z.im = r; }
inline void set_real_part(complex& z,      double r) { z.re = r; }
inline void set_imaginary_part(complex& z, double r) { z.im = r; }


inline complex  complex::operator+  () const
               { return *this; }

inline complex  complex::operator+  (const complex& z) const
               { return complex(re+z.re, im+z.im); }

inline complex  complex::operator+  (double r) const 
               { return complex(r+re,im); }

inline complex& complex::operator+= (const complex& z)
               { re+=z.re; im+=z.im; return (*this); }

inline complex& complex::operator+= (double r) 
               { re += r; return (*this);}

inline complex  complex::operator-  () const       
               { return complex(-re, -im); }

inline complex  complex::operator-  (double r)     
               { return complex(re-r, im); }

inline complex  complex::operator-  (const complex& z) const
               { return complex(re-z.re, im-z.im); }

inline complex& complex::operator-= (const complex& z)
               { re -= z.re; im -= z.im; return (*this); }

inline complex& complex::operator-= (double r) 
               { re -= r;  return (*this); }

inline complex  complex::operator* (const complex& z) const
               { return complex(re*z.re - im*z.im, re*z.im + im*z.re); } 

inline complex  complex::operator* (double r) const
               { return complex(r*re, r*im); }

inline complex& complex::operator*= (const complex& z)
               { double r = re*z.re - im*z.im; im = re*z.im + im*z.re; re = r;  return (*this);}

inline complex& complex::operator*= (double r) 
{ 
  re *= r; im *= r;  return (*this); 
}

inline complex complex::operator/ (const complex& z) const
{
  double r = z.re*z.re + z.im*z.im;
  return complex((re*z.re + im*z.im)/r, (z.re*im - z.im*re)/r);
}

inline complex complex::operator/ (double r) const
{ 
  return complex(re/r, im/r); 
}

inline complex& complex::operator/= (const complex& z)
  { 
  double r = z.re*z.re+z.im*z.im;
  double Xre = (re*z.re+im*z.im)/r;
  im = (z.re*im-z.im*re)/r;
  re = Xre;
  return *this;
  }

inline complex& complex::operator/= (double r) { re /= r; im /= r; return *this; }

inline complex operator+ (double r, const complex& z)
{ 
  return complex(r+z.re, z.im); 
}

inline complex operator- (double r, const complex& z)
{ 
  return complex(r-z.re, -z.im); 
}

inline complex operator* (double r, const complex& z)
{ 
  return complex(r*z.re, r*z.im); 
} 

inline complex operator/ (double r, const complex& z)
{
  double d = z.re*z.re + z.im*z.im;
  return complex((r*z.re)/d, (-z.im*r)/d);
}

inline void add(complex& z, const complex& z1, const complex& z2)
{ 
  z.re = z1.re + z2.re;
  z.im = z1.im + z2.im;
}

inline void add(complex& z,       double   r,  const complex& z1)
{
  z.re = r + z1.re;
  z.im = z1.im;
}


inline void add(complex& z, const complex& z1,       double   r)
{
  z.re = z1.re + r;
  z.im = z1.im;
}


inline void sub(complex& z, const complex& z1, const complex& z2)

{ 
  z.re = z1.re - z2.re;
  z.im = z1.im - z2.im;
}


inline void sub(complex& z,       double   r, const complex& z1)
{
  z.re = r - z1.re;
  z.im = z1.im;
}


inline void sub(complex& z, const complex& z1,       double   r)
{
  z.re = z1.re - r;
  z.im = z1.im;
}


inline void mul(complex& z, const complex& z1, const complex& b)
{
  z.re = z1.re*b.re - z1.im*b.im;
  z.im = z1.re*b.im + z1.im*b.re;
}


inline void mul(complex& z,       double   r, const complex& z1)
{
  z.re = r * z1.re;
  z.im = r * z1.im;
}


inline void mul(complex& z, const complex& z1,       double   r)
{
  z.re = z1.re * r;
  z.im = z1.im * r;
}


inline void div(complex& z, const complex& z1, const complex& z2)
{
  double x = z2.re*z2.re + z2.im*z2.im;
  z.re = (z1.re*z2.re + z1.im*z2.im)/x;
  z.im = (z2.re*z1.im - z2.im*z1.re)/x;
}

  
inline void div(complex& z,       double   r, const complex& b)
{
  double x = b.re*b.re + b.im*b.im;
  z.re = (r * b.re)/x;
  z.im = (-b.im * r)/x;
}


inline void div(complex& z, const complex& z1,       double   r)
{
  z.re = z1.re/r;
  z.im = z1.im/r;
}


inline complex conj(const complex& z1)
{ return complex(z1.re , - z1.im ); }

inline complex conj(const complex& z1, const complex& z2)
{ return complex( z1.re*z2.re+z1.im*z2.im, z1.re*z2.im-z1.im*z2.re); }

inline double AbsNorm(const complex& z)

{ return fabs(z.re) + fabs(z.im); }


inline double square_norm(const complex& z)

{ return z.re*z.re + z.im*z.im; }


inline double norm(const complex& z)

{
  #ifdef _MSC_VER
  // Need to confirm this has no effect on results.
  // Could also wrap hypot() in #pragma waning(disable: ...) then enable:
  return _hypot(z.re, z.im); 
  #else
  return hypot(z.re, z.im);
  #endif  
}


inline double phase(const complex& z)
{
  if(z==0) return 0;
  return atan2(z.im, z.re);
} 

inline void norm(complex& z, double r)
{ 
  r /= sqrt(z.re*z.re + z.im*z.im );
  z.re *= r;
  z.im *= r;
}

inline void phase(complex& z, double r)
{
  double tmp = norm(z);
  z.re = cos(r) * tmp;
  z.im = sin(r) * tmp;
}

inline complex sqrt(const complex& z1)
{
  complex z;
  double r = norm(z1);
  if(r==0.0)
  {
	z = 0;
  }
  else
  {
    double p = phase(z1)/2;
    r = sqrt(r);
    z.re = r*cos(p);
    z.im = r*sin(p);
  }

  return z;
}

inline complex complex::Zexp() const
{
  complex z;
  double e = exp(re);
  z.re = e * cos(im);
  z.im = e * sin(im);
  return z; 
}


inline complex exp(const complex& z1)
{
  complex z;
  double e = exp(z1.re);
  z.re = e * cos(z1.im);
  z.im = e * sin(z1.im);
  return z; 
}

inline complex log(const complex& z)
{ return complex(log(norm(z)), phase(z)); }


inline complex pow(const complex& z, const complex& z1)
{  return complex(exp(z1 * log(z))); }

// ____________________________________________________________________________
// G               Complex Class Container Support Functions
// ____________________________________________________________________________

inline bool complex::operator== (const complex& z1) const
{ return ( (re==z1.re) && (im==z1.im) ); }

inline bool complex::operator!= (const complex& z1) const
{ return ((re!=z1.re) || (im!=z1.im)); }

*/


%extend complex {
%pythoncode %{

def __str__(self):
    sss = "("
    sss = sss + str( self.real() ) + "," + str(self.imag()) + ")\n"
    return( sss )


def __repr__(self):
    sss =""
    sss = "pygamma.complex( " + str( self.real() ) + "," + str(self.imag()) + ")"
    return sss
    
    
def __rmul__( self,other ):
	
	return( self * other )


%}
};
