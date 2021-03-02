/* complex.h ************************************************************-*-c++-*-
**										**
**	                               G A M M A				**
**                                                              		**
**	Complex numbers	                                    Interface		**
**										**
**	Copyright (c) 1990, 2000				                **
**	Tilo Levante, Scott A. Smith						**
**	Eidgenoessische Technische Hochschule					**
**	Labor fur physikalische Chemie						**
**	8092 Zurich / Switzerland		 				**
**										**
**      $Header: $
**										**
*********************************************************************************/

#ifndef   Gcomplex_h_			// Is this file already included?
#define   Gcomplex_h_ 			// If no, then remember it

#if defined(GAMPRAGMA)		// Using the GNU compiler?
#pragma interface			// Then this is the interface
#endif

/******************************************************************************
**                                                                           **
** Description                                                               **
**                                                                           **
**  The class complex defines complex numbers for C++ with the usual         **
**  algebraic operations, higher functions, onversions and IO routines.      **
**                                                                           **
******************************************************************************/

#include <iostream>			// Include input/output streams
#include <cmath>			// Include math.h in libstdc++
#include <string>			// Include strings (for _form)
#include <fstream>			// Inlcude filestreams (ofstream)
#include <GamGen.h>			// Know MSVCDLL (__declspec)

class complex
  {
  double re, im;		// Real, imaginary values
  static bool _phase,		// T=(norm,phase), F*={Re, Im}
              _math,		// T=delimiters, F*=no delimiters
              _science;		// T=5.0e5, F*=5000
  static int  _digits,		// Number of total digits
              _digs_aft_dpoint;	// Number of digits after point
  static std::string _form;	// Output format string
  static std::string _zzer;	// Output zero string

// ____________________________________________________________________________
// i                       CLASS COMPLEX ERROR HANDLING
// ____________________________________________________________________________
 
/* These functions are typical for all GAMMA classes & modules.  They output
   error messages by relaying through common functions found in Gutils.  They 
   all take an error index which triggers specific messages to be output.  An
   optional string may be included.  The flag noret signals whether to add
   a line feed at the end of the error message.                              */

void Zerror(int eidx, int noret=0) const;
void Zerror(int eidx, const std::string& pname, int noret=0) const;
volatile void Zfatality(int eidx) const;

// ____________________________________________________________________________
// ii                     Class Complex Ouptut Format
// ____________________________________________________________________________

void SetForm();
//void SetForm(const std::string& fmt);


public:

MSVCDLL friend inline void Swap(complex& z1, complex& z2);

        // Input           z1, z2 : Two complex numbers
        // Output                 : Exchanges the two numbers

// ____________________________________________________________________________
// A                    COMPLEX CONSTRUCTION, ASSIGNMENT
// ____________________________________________________________________________
              
MSVCDLC          complex();
MSVCDLC          complex(double r,   double i=0.0);
MSVCDLC          complex(const complex& z);
MSVCDLC complex& operator= (const complex& z);
MSVCDLC complex& operator= (double r);

// ____________________________________________________________________________
// B                   COMPLEX ELEMENT ACCESS FUNCTIONS
// ____________________________________________________________________________

/* Here are functions to get and set either the real or imaginary components
   of a complex number.  There are two flavors of this, those that deal with
   the component references directly and those that deal with copies of these.
   Using the references is of course more efficient.  Additionally, provided
   are constant and non-constant versions.  That should cover everything.    

  Function/Operator/Variable                   Return Ability
  -------------------------- --------------------------------------------------
   z.re = z.Relem() = zRe(z) THE real component (z.re needs friend of complex)
   z.im = z.Ielem() = zIm(z) THE imag component (z.im needs friend of complex)
       z.Rec(), z.Imc()      Unchangable real/imaginary component reference
       Re(z), Im(z)          Copy of the real/imaginary component
   Re(z), set_real_part(z)   Set the real component
 Im(z),set_imaginary_part(z) Set the imaginary component                     */

MSVCDLL inline        double& Relem();
MSVCDLL inline        double& Ielem();
MSVCDLL inline const  double& Rec();
MSVCDLL inline const  double& Imc();

MSVCDLL inline        double real();
MSVCDLL inline        double imag();

MSVCDLL friend inline double& zRe(complex& z);
MSVCDLL friend inline double& zIm(complex& z);
MSVCDLL friend inline double  Re(const complex& z);
MSVCDLL friend inline double  Im(const complex& z);

MSVCDLL friend inline void set_real_part(complex& z, double r);
MSVCDLL friend inline void Re(complex& z, double r);
MSVCDLL friend inline void set_imaginary_part(complex& z, double r);
MSVCDLL friend inline void Im(complex& z, double r);
	
// ____________________________________________________________________________
// C                  COMPLEX BASIC ALGEBRAIC PROPERTIES
// ____________________________________________________________________________

/* These functions allow for simple mathematical operations between two
   complex numbers and between a complex number and a double.

   Operator Arguments      Result        Operator Arguments       Result
   -------- --------- -----------------  -------- --------- -------------------
      +         z            z              +       z,z1           z+z1
      +        z,r           z+r            +       r,z            z+r
      +=       z,z1          z+z1           +=        r            z+r
      -         z            -z             -       z,z1           z-z1
      -        z,r           z-r            -       r,z            r-z
      -=       z,z1          z-z1           -=      z,r            z-r
      *        z,z1          z*z1           *       z,r            z*r
      *=       z,z1          z*z1           *=      z,r            z*r
      /        z,z1          z/z1           /       z,r            z/r
      /=       z,z1          z/z1           /=      z,r            z/r      */

MSVCDLL inline complex  operator+  ()                  const;
MSVCDLL inline complex  operator+  (const complex& z)  const;
MSVCDLL inline complex  operator+  (double r)          const;
MSVCDLL inline complex& operator+= (const complex& z);
MSVCDLL inline complex& operator+= (double r);
MSVCDLL inline complex  operator-  ()                  const;
MSVCDLL inline complex  operator-  (double r);
MSVCDLL inline complex  operator-  (const complex& z)  const;
MSVCDLL inline complex& operator-= (const complex& z);
MSVCDLL inline complex& operator-= (double r);
MSVCDLL inline complex  operator*  (const complex& z)  const;
MSVCDLL inline complex  operator*  (double r)          const;
MSVCDLL inline complex& operator*= (const complex& z);
MSVCDLL inline complex& operator*= (double r);
MSVCDLL inline complex  operator/  (const complex& z)  const;
MSVCDLL inline complex  operator/  (double r)          const;
MSVCDLL inline complex& operator/= (const complex& z);
MSVCDLL inline complex& operator/= (double r);

MSVCDLL friend inline complex operator+ (double r, const complex& z);
MSVCDLL friend inline complex operator- (double r, const complex& z);
MSVCDLL friend inline complex operator* (double r, const complex& z);
MSVCDLL friend inline complex operator/ (double r, const complex& z);

// ____________________________________________________________________________
// D              ADDITIONAL BASIC ALGEBRAIC FUNCTIONS
// ____________________________________________________________________________

// Note - The use of these can be faster than the function in the above
//        set when using three numbers which are already defined.  In
//        all cases, the first argument is the result.

        // Input        z, z1, z2 : Three complex number
	//    or
        //              z, z1, r  : Two complex numbers and a real number
	//    or
        //              z, r, z1  : Two complex numbers and a real number
	// Output 		z : Alters z to one of the following
	//			    depending upon the function and overload

	//				z = z1 {+,-,*,/} z2
	//				z = z1 {+,-,*,/} r
	//				z = r  {+,-,*,/} z1

MSVCDLL friend inline void add(complex& z, const complex& z1, const complex& z2);
MSVCDLL friend inline void add(complex& z,       double    r, const complex& z1);
MSVCDLL friend inline void add(complex& z, const complex& z1,       double    r);
MSVCDLL friend inline void sub(complex& z, const complex& z1, const complex& z2);
MSVCDLL friend inline void sub(complex& z,       double    r, const complex& z1);
MSVCDLL friend inline void sub(complex& z, const complex& z1,       double    r);
MSVCDLL friend inline void mul(complex& z, const complex& z1, const complex& z2);
MSVCDLL friend inline void mul(complex& z,       double    r, const complex& z1);
MSVCDLL friend inline void mul(complex& z, const complex& z1,       double    r);
MSVCDLL friend inline void div(complex& z, const complex& z1, const complex& z2);
MSVCDLL friend inline void div(complex& z,       double    r, const complex& z1);
MSVCDLL friend inline void div(complex& z, const complex& z1,       double    r);

// ____________________________________________________________________________
// E                     ADDITIONAL BASIC ALGEBRAIC FUNCTIONS
// ____________________________________________________________________________
 
// ----------------------------------------------------------------------------
//                            Member Functions
// ----------------------------------------------------------------------------

/* These must be distinguished from the friend functions in this section. The
   member functions, if specified, will reside in the "gam" namespace and then
   dominate function usage within this class. For example, conj can then be
   used only for namespace gam complex numbers. This is not a problem for conj,
   but it would be a problem for functions like exp and log because we need
   those functions from namespace std to be at our disposal within this class.
   Consequently, any functions with names that overlap those in std that we
   need herein will be declared as friend functions.                         */ 
  
/*  Function/Operator                         Result
    -----------------   ------------------------------------------------------
         conj           Complex conjugate of z: a+ib --> a-ib
      conj_times(z)     Complex conj(z)*z1:     a+ib,c+id --> ac+bd+i(ad-bc)  
        conj(z)         Complex conjugate of z: a+ib --> a-ib
       conj(z,z1)       Complex conj(z)*z1:     a+ib,c+id --> ac+bd+i(ad-bc)
       AbsNorm(z)       Absolute norm of z:     a+ib -> |a| + |b|
     square_norm(z)     Square norm of z:       a+ib -> a*a + b*b
        norm(z)         Norm of z:              a+ib -> [a*a + b*b]^1/2
       phase(z)         Phase of z:             a+ib ->
       norm(z,r)        Norm of z:              norm set to r
       phase(r)         Phase of z:             phase set to r               */

MSVCDLL complex conj() const;
MSVCDLL complex conj_times(const complex& z) const;

MSVCDLL friend inline complex conj(const complex& z);
MSVCDLL friend inline complex conj(const complex& z, const complex& z1);

MSVCDLL friend inline double AbsNorm(const     complex& z);
MSVCDLL friend inline double square_norm(const complex& z);
MSVCDLL friend inline double norm(const        complex& z);
MSVCDLL friend inline double phase(const       complex& z);
MSVCDLL friend inline void   norm(complex& z,  double r);
MSVCDLL friend inline void   phase(complex& z, double r);


MSVCDLL friend inline complex sqrt(const complex& z);

	// Input 		z : A complex number, unaltered
	// Output 	       z1 : Returns the square root of z, z1 = sqrt(z)


MSVCDLL inline complex Zexp() const;

	// Input 		z : A complex number, unaltered
	// Output 	       z1 : Returns the exponential of z, z1 = exp(z)

MSVCDLL friend inline complex exp(const complex& z);

	// Input 		z : A complex number, unaltered
	// Output 	       z1 : Returns the exponential of z, z1 = exp(z)


MSVCDLL friend inline complex log(const complex& z);

	// Input 		z : A complex number, unaltered
	// Output 	       z1 : Returns the natural log of z, z1 = ln(z)
	// Errors		  : Error for z = (0,0)


MSVCDLL friend inline complex pow(const complex& z, const complex& z1);

	// Input 	     z,z1 : Complex numbers, unaltered             z1
	// Output 	       z2 : Returns the z to the z1 power, z2 = (z)



// ____________________________________________________________________________
// F           HYPERBOLIC & INVERSE TRIGONOMETRIC FUNCTIONS
// ____________________________________________________________________________

MSVCDLL complex Zsin()   const;
MSVCDLL complex Zcos()   const;
MSVCDLL complex Ztan()   const;
MSVCDLL complex Zasin()  const;
MSVCDLL complex Zacos()  const;
MSVCDLL complex Zatan()  const;
MSVCDLL complex Zsinh()  const;
MSVCDLL complex Zcosh()  const;
MSVCDLL complex Ztanh()  const;
MSVCDLL complex Zasinh() const;
MSVCDLL complex Zacosh() const;
MSVCDLL complex Zatanh() const;

MSVCDLL friend complex sin(const   complex& z);
MSVCDLL friend complex cos(const   complex& z);
MSVCDLL friend complex tan(const   complex& z);
MSVCDLL friend complex asin(const  complex& z);
MSVCDLL friend complex acos(const  complex& z);
MSVCDLL friend complex atan(const  complex& z);
MSVCDLL friend complex sinh(const  complex& z);
MSVCDLL friend complex cosh(const  complex& z);
MSVCDLL friend complex tanh(const  complex& z);
MSVCDLL friend complex asinh(const complex& z);
MSVCDLL friend complex acosh(const complex& z);
MSVCDLL friend complex atanh(const complex& z);

// ____________________________________________________________________________
// G                      CLASS COMPLEX I/O FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                     Functions To Alter The Output Format
// ----------------------------------------------------------------------------

/* The format used to output complex numbers is set using five factors. These
   are
           Variable     Type                Function
       ===============  ====     ============================================
       _phase           bool      true: R*e(i*phi),   false: a+ib
       _math            bool      true: dellimiters,  false: no delimiters
       _science         bool      true: sci.notation, false: no sci. notation
       _digits          int       Number of total digits (double)
       _digs_aft_dpoint int       Number of digits after point (double)

   The _digits, _digs_aft_dpoint, and science are used in combination to make
   a format string, e.g. %8.4f. Here are how things may appear:

   _phase T, _math T, _science F: 5.0*exp(4.0i)  where |z|=5 & phase=4
   _phase T, _math F, _science F: [5.0,4.0]      where |z|=5 & phase=4
   _phase F, _math T, _science F: 5.0+4.0i       where a=5 and b=4
   _phase F, _math F, _science F: (5.0,4.0)      where a=5 and b=4 (default)
   
   _digits will adjust how many digits are output with each double value.
   _digs_aft_dpoint will adjust how many digits follow the decimal point
   _science will affect whether the ouptut format for each double to be
            either %(_digits._digs_aft_dpoint)f or %(_digits._digs_aft_dpoint)e
     
   The default value of the above three is 6,2,false setting the format to
   be %6.2f = _form.                                                            */

MSVCDLL static void Reiphi(bool TF);
MSVCDLL static void delim(bool  TF);
MSVCDLL static void scinot(bool TF);
MSVCDLL static void digits(int  digs);
MSVCDLL static void dadp(int    adp);

// ----------------------------------------------------------------------------
//                    Functions To Facilite Complex Output
// ----------------------------------------------------------------------------

MSVCDLL static int         dlength();			// Length of printed d
MSVCDLL static std::string dformat();			// Format of printed d
MSVCDLL static int         zlength();			// Length of printed z
MSVCDLL static bool        normphase();		// True if R,phi out


//static void form(std::string& form);
// friend void complex_setf(int phase, int math, int science,
//			                    int digits, int digs_aft_dpoint );

MSVCDLL static void complex_getf(bool& phase, bool& math, bool& science,
			                  int &digits, int &digs_aft_dpoint );

MSVCDLL friend void complex_getf(std::string& form);
 
        // Input                form  : A string for the output format 
        // Output                void : This function just returns the
	//				string _form, for testing purposes.



// ----------------------------------------------------------------------------
//               Functions To Actually Perform Complex Output
// ----------------------------------------------------------------------------

	// Input 		ostr : An output stream
	//			z    : A complex number, unaltered
	// Output 		ostr : The modified output stream containing
	//			       complex number z in the format set
	//			       in complex_setf function

MSVCDLL        std::string   printString() const;
MSVCDLL        std::ostream& print(std::ostream& ostr) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const complex& z);

	// Input 		istr : An input stream
	//			z    : A complex number
	// Output 		istr : The complex number z is set
	//			       from the input stream istr
	//			       & the modified input stream
	//			       returned.
	// Note			     : Complex number z is read in the
	//			       format 5 4 as 5 + 4i
	// Bug			     : Currently no error handling

MSVCDLL friend std::istream& operator>> (std::istream& istr, complex& z);


        // Input                z     : Complex number (this)
        //                      fn    : Filename
        //                      fp    : File stream (pointing at z spot)
        // Return               void  : The complex number z is written
        //                              to a file called filename, or
        //                              into file fp at current location
        // Note                       : The file format is BINARY

MSVCDLL void write(const std::string& fn);
MSVCDLL void write(std::ofstream& fp) const;

        // Input                z     : Complex number (this)
        //                      fn    : Filename
        //                      fp    : File stream (pointing at z spot)
        // Return               void  : The complex number z is read
        //                              from a file called filename
        //                              or from file fp at current location

MSVCDLL void read(const std::string& fn);
MSVCDLL void read(std::ifstream& fp);

// ____________________________________________________________________________
// H         Complex Class Container Support Functions
// ____________________________________________________________________________

MSVCDLL inline bool operator== (const complex& z) const;
MSVCDLL inline bool operator!= (const complex& z) const;
MSVCDLL        bool operator<  (const complex& z) const;
MSVCDLL        bool operator>  (const complex& z) const;

  };

/************************************************************************/
/************************************************************************/
/*                       CLASS COMPLEX CONSTANTS			*/
/************************************************************************/
/************************************************************************/

extern const MSVCDLL complex complex0;		// z = 0 : (0,0)
extern const MSVCDLL complex complex1;		// z = 1 : (1,0)
extern const MSVCDLL complex complexi;		// z = i : (0,1)


/************************************************************************/
/************************************************************************/
/*                  CLASS COMPLEX INLINE FUNCTIONS			*/
/************************************************************************/
/************************************************************************/

// ____________________________________________________________________________
// A               COMPLEX CONSTRUCTION, ASSIGNMENT
// ____________________________________________________________________________

inline complex::complex() { }
inline complex::complex(double r, double i) :re(r),im(i) { }
inline complex::complex(const complex& z) :re(z.re),im(z.im) { }
inline complex& complex::operator= (const complex &z)
  { re = z.re; im = z.im; return *this; }
inline complex& complex::operator= (double r)
  { re = r; im = 0;  return *this; }


        // Input           z1, z2 : Two complex numbers
        // Output                 : Exchanges the two numbers

inline void Swap(complex& z1, complex& z2)
  {
  double tmp;
  tmp = z1.re; z1.re = z2.re; z2.re = tmp;
  tmp = z1.im; z1.im = z2.im; z2.im = tmp;
  }

// ____________________________________________________________________________
// B                   COMPLEX ELEMENT ACCESS FUNCTIONS
// ____________________________________________________________________________

/* Here are functions to get and set either the real or imaginary components
   of a complex number.  There are two flavors of this, those that deal with
   the component references directly and those that deal with copies of these.
   Using the references is of course more efficient.  Additionally, provided
   are constant and non-constant versions.  That should cover everything.    

  Function/Operator/Variable                   Return Ability
  -------------------------- --------------------------------------------------
   z.re = z.Relem() = zRe(z) THE real component (z.re needs friend of complex)
   z.im = z.Ielem() = zIm(z) THE imag component (z.im needs friend of complex)
       z.Rec(), z.Imc()      Unchangable real/imaginary component reference
       Re(z), Im(z)          Copy of the real/imaginary component
   Re(z), set_real_part(z)   Set the real component
 Im(z),set_imaginary_part(z) Set the imaginary component                     */

inline        double& complex::Relem()     { return re; }
inline        double& complex::Ielem()     { return im; }
inline const  double& complex::Rec()       { return re; }
inline const  double& complex::Imc()       { return im; }

inline        double complex::real()       { return re; }
inline        double complex::imag()       { return im; }

inline        double& zRe(complex& z)      { return z.re; }
inline        double& zIm(complex& z)      { return z.im; }
inline        double  Re(const complex& z) { return z.re; }
inline        double  Im(const complex& z) { return z.im; }

inline void Re(complex& z,                 double r) { z.re = r; }
inline void Im(complex& z,                 double r) { z.im = r; }
inline void set_real_part(complex& z,      double r) { z.re = r; }
inline void set_imaginary_part(complex& z, double r) { z.im = r; }

// ____________________________________________________________________________
// C              COMPLEX BASIC ALGEBRAIC PROPERTIES
// ____________________________________________________________________________

/* These functions allow for simple mathematical operations between two
   complex numbers and between a complex number and a double.

   Operator Arguments      Result        Operator Arguments       Result
   -------- --------- -----------------  -------- --------- -------------------
      +         -             z             +        z1           z+z1
      +         r            z+r            +=       z1           z+z1
      +=        r            z+r            -         -            -z
      -      z1,mx2      z1-mx2           *       z,mx           z*mx
      *      z1,mx2      z1*mx2           *       mx,d           d*mx
      /      z1,mx2      z1*inv(mx2)      *       d,mx           d*mx
      /      z,z         (1/z)*mx          *       mx,d         (1/d)*mx
      /      z,mx         z*inv(mx)         *       d,mx         d*inv(mx)  */



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
                                          { re *= r; im *= r;  return (*this); }
inline complex complex::operator/ (const complex& z) const
  {
  double r = z.re*z.re + z.im*z.im;
  return complex((re*z.re + im*z.im)/r, (z.re*im - z.im*re)/r);
  }

inline complex complex::operator/ (double r) const
  { return complex(re/r, im/r); }
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
                                           { return complex(r+z.re, z.im); }
inline complex operator- (double r, const complex& z)
                                           { return complex(r-z.re, -z.im); }
inline complex operator* (double r, const complex& z)
                                           { return complex(r*z.re, r*z.im); } 

inline complex operator/ (double r, const complex& z)
  {
  double d = z.re*z.re + z.im*z.im;
  return complex((r*z.re)/d, (-z.im*r)/d);
  }

// ____________________________________________________________________________
//                COMPLEX BASIC ALGEBRAIC PROPERTIES
// ____________________________________________________________________________


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


// ____________________________________________________________________________
// D               ADDITIONAL BASIC ALGEBRAIC FUNCTIONS
// ____________________________________________________________________________


inline complex conj(const complex& z1)
  { return complex(z1.re , - z1.im ); }

        // Input                z1 : A complex number, unaltered        *
        // Output                z : The complex conjugate of z1, z = z1


inline complex conj(const complex& z1, const complex& z2)
  { return complex( z1.re*z2.re+z1.im*z2.im, z1.re*z2.im-z1.im*z2.re); }

        // Input           z1, z2 : Two complex numbers, unaltered
        // Output               z : The returns conj(z1)*z2;



inline double AbsNorm(const complex& z)

        // Input                z : A complex number, unaltered
        // Output               r : Absolute norm of z, r = |re| + |im|

  { return fabs(z.re) + fabs(z.im); }


inline double square_norm(const complex& z)

        // Input                z : A complex number, unaltered
	//						      2       2
        // Output               r : Square norm of z, r = |re|  + |im|

  { return z.re*z.re + z.im*z.im; }


inline double norm(const complex& z)

        // Input                z : A complex number, unaltered
	//						[     2       2 ]
        // Output               r : Norm of z, r = sqrt | |re|  + |im|  |
	//						[               ]

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
  if(r==0.0) z = 0;
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


#endif								// complex.h
