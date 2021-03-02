/* complex.cc ***************************************************-*-c++-*-
**                                                                      **
** 	                              G A M M A                         **
**                                                                      **
**	Complex numbers                                 Implementation  **
**                                                                      **
**	Copyright (c) 1990                                              **
**	Tilo Levante, Scott A. Smith                                    **
**	Eidgenoessische Technische Hochschule                           **
**	Labor fur physikalische Chemie                                  **
**	8092 Zurich / Switzerland                                       **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**   Description                                                        **
**                                                                      **
** The class complex defines complex numbers for C++ with the usual     **
** algebraic operations, functions, conversion,	and I/O routines.       **
**                                                                      **
*************************************************************************/

#ifndef   Gcomplex_cc_			// Is the file already included?
#  define Gcomplex_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <Matrix/MxModBas.h>		// Include Matrix module errors
#include <Matrix/complex.h>		// Include the interface
#include <iostream>			// Include C++ IO streams

using std::string;			// Using libstdc++ strings

					// These set the default output form
bool   complex::_phase=false,		// T={norm,phase}, F={Re,Im}
       complex::_math=false,		// T=delimiters, F=no delimtiters
       complex::_science=false;		// Scientific notation
int    complex::_digits=6,		// Total digits, Digits after decimal
       complex::_digs_aft_dpoint=2;
string complex::_form="%6.2f";		// Output format string
string complex::_zzer="   0  ";		// Zero component


/************************************************************************/
/************************************************************************/
/*                       CLASS COMPLEX CONSTANTS                        */
/************************************************************************/
/************************************************************************/

const complex complex0(0,0);		// z = 0 : (0,0)
const complex complex1(1,0);		// z = 1 : (1,0)
const complex complexi(0,1);		// z = i : (0,1)

// ____________________________________________________________________________
// i                     CLASS COMPLEX ERROR HANDLING
// ____________________________________________________________________________

        // Input                z	: Complex number (this)
        // 			eidx	: Error index
	//			noret	: Flag for linefeed (0=linefeed)
	//			pname  : string in message

void complex::Zerror(int eidx, int noret) const
  {
  string hdr("Complex Number");
  switch (eidx)
    {
    case 0: MxModError(hdr, 0, noret); break;	// Program Aborting        (0)
    case 3: MxModError(hdr, 3, noret); break;	// !Construct From Pset    (3)
    case 4: MxModError(hdr, 4, noret); break;	// !Construct From File    (4)
    case 5: MxModError(hdr, 5, noret); break;	// !Write To Param File    (5)
    case 6: MxModError(hdr, 6, noret); break;	// !Write To Filestream    (6)
    default: MxModError(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }

void complex::Zerror(int eidx, const string& pname, int noret) const
  {
  string hdr("Complex Number");
  string msg;
  switch(eidx)
    {
    case 1:  MxModError(hdr, 1, pname, noret);   break;	// File Problems   (1)
    default: MxModError(hdr, -11, pname, noret); break;	// Unknown Error   (-1)
    }
  }
 
volatile void complex::Zfatality(int eidx) const
  {  
  Zerror(eidx, 1);				// Normal non-fatal error
  if(eidx) Zerror(0);				// Program aborting error
  MxModFatal();					// Keep screen nice, exit
  }

// ____________________________________________________________________________
// ii                     Class Complex Ouptut Format
// ____________________________________________________________________________


void complex::SetForm()				// Sets _form from three values
  {						// _digits, _digs_aft_dpoint 
  char ef = (_science)?'e':'f';			// and _science
  if(_digs_aft_dpoint<0)
    _form = string("%") + MxModdec(_digits) + ef;
  else
    _form = string("%") + MxModdec(_digits)
          + string(".") + MxModdec(_digs_aft_dpoint)
          + ef;
  }

//void complex::SetForm(const string& fmt)	// From fmt, sets _form as
//  {						// well as the three values
//  }


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

/*               These Are All Inlined In The Header File

inline        double& complex::Relem();
inline        double& complex::Ielem();
inline const  double& complex::Rec();
inline const  double& complex::Imc();
friend inline double& zRe(complex& z);
friend inline double& zIm(complex& z);
friend inline double  Re(const complex& z);
friend inline double  Im(const complex& z);

friend inline void set_real_part(complex& z, double r);
friend inline void Re(complex& z, double r);
friend inline void set_imaginary_part(complex& z, double r);
friend inline void Im(complex& z, double r);                                 */

// ____________________________________________________________________________
// D                     ADDITIONAL BASIC ALGEBRAIC FUNCTIONS
// ____________________________________________________________________________
      
/*  Function/Operator                         Result
    -----------------   ------------------------------------------------------
         conj           Complex conjugate of z: a+ib --> a-ib
      conj_times        Complex conj(z)*z1:     a+ib,c+id --> ac+bd+i(ad-bc)

*/
  
complex complex::conj() const { return complex(re, -im); }

complex complex::conj_times(const complex& z) const
  { return complex( re*z.re+im*z.im, re*z.im-im*z.re); }

// ____________________________________________________________________________
// E                           TRIGONOMETRIC FUNCTIONS
// ____________________________________________________________________________

/* Note that both member and non-member forms of these functions are supplied.
   To avoid naming conflicts, the member function names use a Z at their start.
   The non-member names are more commonly invoked and these are the ones that
   contain the explicit code to produce the desired result.                   */

complex complex::Zsin()   const { return sin(*this);   }
complex complex::Zcos()   const { return cos(*this);   }
complex complex::Ztan()   const { return tan(*this);   }
complex complex::Zasin()  const { return asin(*this);  }
complex complex::Zacos()  const { return acos(*this);  }
complex complex::Zatan()  const { return atan(*this);  }
complex complex::Zsinh()  const { return sinh(*this);  }
complex complex::Zcosh()  const { return cosh(*this);  }
complex complex::Ztanh()  const { return tanh(*this);  }
complex complex::Zasinh() const { return asinh(*this); }
complex complex::Zacosh() const { return acosh(*this); }
complex complex::Zatanh() const { return atanh(*this); }

complex sin(const complex& z)
  {
  double a = exp(z.im);
  double b = 1/a;
  return complex(sin(z.re) * (a+b)/2, cos(z.re) * (a-b)/2); 
  }

complex cos(const complex& z)

  { 
  double a = exp( z.im );
  double b = 1/a;
  return complex(cos(z.re) * (a+b)/2, -sin(z.re) * (a-b)/2 ); 
  }

complex tan(const complex& z )
  { 
  double a = exp(2 * z.im );
  double b = 1/a;
  double d = cos(2 * z.re ) + (a+b)/2;
  return complex(sin(2*z.re) / d, (a-b)/(2*d)); 
  } 

// _________________________ Inverse Trigonometric ____________________________

complex asin(const complex& z)
  { return - complexi*log(complexi*z + sqrt (1-z*z)); }

complex acos(const complex& z)
  { return complex(0,-1)*log(z + sqrt (z*z-1)); }

complex atan(const complex& z)
  { return complex(0,-0.5)*log((1+complexi*z)/(1-complexi*z)); }

// ______________________________ Hyperbolics _________________________________


complex sinh(const complex& z)
  { return complex( sinh(z.re)*cos(z.im),cosh(z.re)*sin(z.im) ); }

complex cosh(const complex& z)
  { return complex( cosh(z.re)*cos(z.im),sinh(z.re)*sin(z.im) ); }

complex tanh(const complex& z)
  {
  double tmp = cos(2*z.re) + cosh(2*z.im);
  return complex( sin(2*z.re)/tmp, sinh(2*z.im)/tmp );
  } 

// __________________________ Inverse Hyperbolics _____________________________

complex asinh(const complex& z) { return log(z + sqrt(z*z + 1)); }
complex acosh(const complex& z) { return log(z + sqrt(z*z - 1)); }
complex atanh(const complex& z) { return log((1+z)/(1-z))/2; }

// ____________________________________________________________________________
// F                       CLASS COMPLEX I/O FUNCTIONS
// ____________________________________________________________________________

void complex::Reiphi(bool TF)  { _phase   = TF;    }
void complex::delim(bool  TF)  { _math    = !TF;    }
void complex::scinot(bool TF)  { _science = TF; complex z; z.SetForm(); }
void complex::digits(int digs) { if(digs>0) _digits = digs; complex z; z.SetForm(); }
void complex::dadp(int   adp)  { if(adp<_digits-2) _digs_aft_dpoint=adp; complex z; z.SetForm(); }

bool complex::normphase()      { return _phase;  }
int    complex::dlength()      { return _digits; }
string complex::dformat()      { return _form;   }

int    complex::zlength()
{
  if(_phase)
  {
    if(_math) 
      return 2*_digits+7;			// R*exp(PHIi)
    else      
      return 2*_digits+3;			// [R,PHI]
  }

  if(_math)   
    return 2*_digits+2;			    // a + bi
  else        
    return 2*_digits+3;			    // (a,b)
}

/*                          ASCII I/O Based Functions                        */


        // Input                phase : TRUE  == output as norm and phase
        //                              FALSE == output as real & imag.   (*)
        //                       math : TRUE  == output as 5+4i or 5*exp(4i)
        //                              FALSE == output as (5,4) or [5,4] (*)
        //                    science : TRUE  == output as 5.0e5
        //                              FALSE == output as 500000         (*)
        //                     digs : number of digits total            (6)
        //            digs_aft_dpoint : digits after decimal point        (2)
        //                              (<0 means as many as possible)
        // Output                void : This sets the static values of
        //                              _phase, _math, _science, _digits,
        //                              and _digs_aft_dpoint.  As a
        //                              result, the printed form of complex
        //                              numbers are specified.

/*
void complex_setf(bool phase, bool math, bool science,
                                               int digs, int digs_aft_dpoint)
  {
  complex z;
  (phase)?z.Reiphi():z.aib();
  z.delim(!math);
  z.scinot(science);
  z.digits(digs);
  z.dadp(digs_aft_dpoint);
  char ef = (science)?'e':'f';
  if(digs_aft_dpoint<0)
    complex::_form = string("%") + MxModdec(digs) + ef;
  else
    complex::_form = string("%") + MxModdec(digs)
                   + string(".") + MxModdec(digs_aft_dpoint)
                   + ef;
  complex::_zzer   = string(digs/2, ' ') + string("0")
                   + string(digs-digs/2-1, ' ');

  }
*/

        // Input                phase : TRUE  == output as norm and phase
        //                              FALSE == output as real and imaginary
        //                      math  : TRUE  == output as 5+4i or 5*exp(4i)
        //                              FALSE == output as (5,4) or [5,4]
        //                    science :TRUE  == output as 5.0e5
        //                              FALSE == output as 500000
        //                     digits : number of digits total
        //			dadp  : digits after decimal point 
        //                              (<0 means as many as possible)
        // Output                void : This function sets the input variables
        //                              to the current static values of
        //                              _phase, _math, _science, _digits,
        //                              and _digs_aft_dpoint.

void complex::complex_getf(bool& phase, bool& math, bool& science, int& digits, int& dadp)
  {
  phase   =complex::_phase;
  math    =complex::_math;
  science =complex::_science;
  digits  =complex::_digits;
  dadp    =complex::_digs_aft_dpoint;
  }

        // Input                ostr : An output stream
        //                      z    : A complex number, won't be changed
        // Output               ostr : The modified output stream containing
        //                             complex number z in the format set
        //                             in complex_setf function

string complex::printString() const
  {
  string PS;				// Our print string

//            Set Output String Using R and Phi
//   _math==TRUE: R*exp(i*phi)    _math==FALSE: [R,phi]

  if(complex::_phase)				// Here _math=T=_phase
    { 
    double p = phase(*this);			// Get complex phase
    double n = norm(*this);			// Get complex norm
    if(complex::_math)				// Here _math=TRUE=_phase
      { 					// Format: R*exp(i*phi)
      if(n==0) 					// Length: 2*digits+7
        PS = _zzer				//    Set norm (zero)
           + string(_digits+7, ' ');		//    Set phase (empty)
      else
	{
	PS = MxModform(_form.c_str(),n)		//    Set norm	
	   + "*exp(";				//    phase header
	if(p==0)				//    If phase is zero
         PS += _zzer;				//    we write exactly 0
        else					//    If phase not zero
         PS+=MxModform(_form.c_str(),p);	//   we write the number
        PS += "i)";				//    phase end
	}
      return PS;				// We are done if writing
      }						// as R*exp(i*phi)

    else					// Here _math=F, phase=T
      {						// Format: [R,phi]
               PS =  "[";			// Length: 2*digits+3
      if(n==0) PS += _zzer;
      else     PS += MxModform(complex::_form.c_str(),n);
               PS += ",";
      if(p==0) PS += _zzer;
      else     PS += MxModform(complex::_form.c_str(),p);
               PS += "]";
      return PS;				// We are done if
      }						// as [R, phi]
    }

//                Set Output String Using a and b
//         _math==TRUE: Re+iIm    _math==FALSE: (Re,Im)

  if(complex::_math) 				// Here _math=T, phase=F
    {						// Format: re + im i
    if(re==0) PS = _zzer;			// Length: 2*digits+2 
    else      PS =  MxModform(_form.c_str(),re);
    if(im>=0) PS += "+";
    else      PS += "-";
    if(im==0) PS += _zzer;
    else      PS += MxModform(_form.c_str(),fabs(im));
              PS += "i";
    }
  else						// Here _math = FALSE = _phase
    {						// Format: (re,im)
              PS =  "(";			// Length: 2*digits+3
    if(re==0) PS += _zzer;
    else      PS += MxModform(_form.c_str(),re);
              PS += ",";
    if(im==0) PS += _zzer;
    else      PS += MxModform(_form.c_str(),im);
              PS += ")"; 
    }
  return PS;
  }

std::ostream& complex::print(std::ostream& ostr) const 
  { ostr << printString(); return ostr; }

std::ostream& operator << (std::ostream& ostr, const complex& z)
  { return z.print(ostr); }

        // Input                istr : An input stream
        //                      z    : A complex number
        // Output               istr : The complex number z is set
        //                             from the input stream istr
        //                             & the modified input stream
        //                             returned.
        // Note                      : Complex number z is read in the
        //                             format 5 4 as 5 + 4i
        // Bug                       : Currently no error handling

std::istream& operator >> (std::istream& istr, complex& z) {return istr >> z.re >> z.im;}

/*                         Binary I/O Based Functions                        */



        // Input                z     : Complex number (this)
        //                      fn    : Filename
        //                      fp    : File stream (pointing at z spot)
        // Return               void  : The complex number z is written
        //                              to a file called filename or
        //                              into file fp at current location
        // Note                       : The file format is BINARY

void complex::write(const string& fn)
  {
  std::ofstream fp;					// Construct a file
  fp.open(fn.c_str(),std::ios::out|std::ios::binary);	// Open file
  if(!fp)					// Check if file is OK
    {
    Zerror(1, fn, 1);				// File Problems 
    Zfatality(6);				// Can't write to output file
    }
  write(fp);					// Write z, use overload								// This makes warnings on egcs! sosi
  fp.close();					// Close file
  }

void complex::write(std::ofstream& fp) const
  {                                                               
  fp.write((char *) &re, sizeof(double));	// Write real component
  fp.write((char *) &im, sizeof(double));	// Write imaginary component
  }

        // Input                z     : Complex number (this)
        //                      fn    : Filename
        //                      fp    : File stream (pointing at z spot)
        // Return               void  : The complex number z is read
        //                              from a file called filename or
        //                              from file fp at current location

void complex::read(const string& fn)
  {                                                                  
  std::ifstream fp;					// Construct a file
  fp.open(fn.c_str(),std::ios::in|std::ios::binary); 	// Open file
  if(!fp)					// Check if file is OK
    {
    Zerror(1, fn, 1);		 		// File problems 
    Zfatality(11);				// Can't input from file
    }
  read(fp);					// Read z , use overload
  fp.close();					// Close file
  }

void complex::read(std::ifstream& fp)
  {                                                                  
  if(!fp)					// Check if filestream is OK
    {
    Zerror(1, 1);		 		// Input filestream problems 
    Zfatality(11);				// Can't input from file
    }
  if(fp.eof())
    {
    Zerror(1, 1);		 		// Input filestream problems 
    Zerror(12, 1);		 		// End of file reached
    Zfatality(11);				// Can't input from file
    }
  fp.read((char *)&re,sizeof(double));		// Read real component
  if(fp.eof())
    {
    Zerror(1, 1);		 		// Input filestream problems 
    Zerror(12, 1);		 		// End of file reached
    Zfatality(11);				// Can't input from file
    }
  fp.read((char *)&im,sizeof(double));		// Read imaginary component
  }

// ____________________________________________________________________________
// G               Complex Class Container Support Functions
// ____________________________________________________________________________

// bool complex::operator==(const complex& z) const     // INLINED in header
// bool complex::operator!=(const complex& z) const     // INLINED in header
 
bool complex::operator<(const complex& z) const
  { return (norm(*this) < norm(z)); }

bool complex::operator>(const complex& z) const
  { return (norm(*this) > norm(z)); }

#endif                                                            // complex.cc

