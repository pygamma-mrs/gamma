/* row_vector.cc ************************************************-*-c++-*-
**  									**
**                                G A M M A                    		**
**                                            				**
**      Row Vectors                                     Implementation	**
**                                            				**
**      Copyright (c) 1990, 2002                              		**
**      Tilo Levante & Scott A. Smith					**
**      Eidgenoessische Technische Hochschule           		**
**      Labor fuer physikalische Chemie                 		**
**      8092 Zuerich / Switzerland                      		**
**                                                      		**
**      $Header: $
**                                                      		**
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  Class row_vector defines row vectors for C++ with the usual         **
**  algebraic operations, higher functions and I/O routines.            **
**                                                                      **
*************************************************************************/

#ifndef   Grow_vector_cc_		// Is file already included?
#  define Grow_vector_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// this is the implementation
#  endif

#include <Matrix/row_vector.h>		// Include the interface
#include <Matrix/MxModBas.h>		// Include Matrix module errors
#include <Matrix/complex.h>		// Include GAMMA complex numbers
#include <vector>			// Include libstdc++ STL vectors
#include <cmath>			// Inlcude HUGE_VAL
#include <stdlib.h>


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                     CLASS ROW VECTOR ERROR HANDLING
// ____________________________________________________________________________
 

/*      Input 			rvec	: A row_vector (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : string in message
        Output                  none    : Error message output
                                          Program execution stopped if fatal
 
  The following error messages use the defaults set in the Gutils package

                Case                          Error Message

                (0)                     Program Aborting.....
                (9)                     Problems During Construction
                default                 Unknown Error                        */

using namespace std;

// Added this to handle the naming ambiguity of "exp" 
// from inside the row_vector class. DCT 10/05/09.
complex complex_exp(const complex& ab)
{
	return exp(ab);
}


void row_vector::RVerror(int eidx, int noret) const
  {                                                                             
  std::string hdr("Row Vector");
  std::string msg;
  switch (eidx)
    {
    case 29: MxModError(hdr,"Operation Is Disallowed",     noret); break;// (29)
    case 30: MxModError(hdr,"Multiple Rows Disallowed",    noret); break;// (30)
    case 31: MxModError(hdr,"Cannot Assign From Col Vec",  noret); break;// (31)
    case 40: MxModError(hdr,"Cannot Perform Construction", noret); break;// (40)
    case 41: MxModError(hdr,"Matrix-Vector Size Mismatch", noret); break;// (41) 
    case 42: MxModError(hdr,"Vector-Vector Size Mismatch", noret); break;// (42) 
    case 43: MxModError(hdr,"Can't Multiply 2 Col Vectors",noret); break;// (43)
    case 44: MxModError(hdr,"Trace On Vector Not Allowed", noret); break;// (44)
    default: MxModError(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }  
 
/* The following error messages use the defaults set in the Gutils package
 
                Case                          Error Message
 
                (1)                     Problems With File PNAME
                (2)                     Cannot Read Parameter PNAME
                (3)                     Invalid Use Of Function PNAME
                default                 Unknown Error - PNAME                */  
 
void row_vector::RVerror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Row Vector");
  std::string msg;
  switch(eidx)
    {
    case 5: msg = std::string("Bad Use Of ") + pname + std::string(" Function ");
             MxModError(hdr, msg, noret);  break;                      // (5)
    case 20 :msg = std::string("Construction From ") + pname + std::string(" Matrix");
             MxModError(hdr, msg, noret);  break;                      // (20)
    case 21: msg = std::string("Problems Using ") + pname + std::string(" Function ");
             MxModError(hdr, msg, noret);  break;                      // (21)
    default: MxModError(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }                              
 
volatile void row_vector::RVfatality(int eidx) const
  {
  RVerror(eidx, 1);                             // Normal non-fatal error
  if(eidx) RVerror(0);                          // Program aborting error
  MxModFatal();                                 // Keep screen nice, exit
  }

volatile void row_vector::RVfatality(int eidx, const std::string& pname) const
  {
  RVerror(eidx, pname, 1);                      // Normal non-fatal error
  if(eidx) RVerror(0);                          // Program aborting error
  MxModFatal();                                 // Keep screen nice, exit
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                 ROW VECTOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/*    Arguments                         Constructed Row Vector
    -------------       -------------------------------------------------------
          -             An empty row vector (size = 0)
          nc            A row vector with nc elements (columns) uninitialized
        nc, z           A row vector with nc elements (columns) all set to z
         rvec           A row vector equivalent to the input vector rvec
          mx            A row vector equivalent to input matrix mx
                        (in this case mx MUST be of dimension 1xn)           */

row_vector::row_vector( )                        : matrix()              { }
row_vector::row_vector(int nc)                   : matrix(1,nc)          { }
row_vector::row_vector(int nc,       double   d) : matrix(1,nc,d)        { }
row_vector::row_vector(int nc, const complex& z) : matrix(1,nc,z)        { }
row_vector::row_vector(const row_vector& rvec)   : matrix((matrix&)rvec) { }
row_vector::row_vector(const matrix& mx)         : matrix(mx)
  {
  if(mx.rows() != 1)
    {
    std::string em = std::string(MxModdec(mx.rows()))
              + std::string("x")
              + std::string(MxModdec(mx.cols()));
    RVerror(20,em,1);			        // Construct multirow array
    RVfatality(40);				// Can't do construction
    }
  }

row_vector::~row_vector () { }

// ____________________________________________________________________________
// B                    ROW VECTOR ACCESS AND ASSIGNMENT
// ____________________________________________________________________________

/* Function Arguments      Result        Function Arguments    Result
   -------- --------- -----------------  -------- --------- -------------------
     ( )        i     The ith element      get        i     Copy of ith element
    getRe       i     Real(ith element)   getIm       i     Imag(ith element)
    put        z,i    Set z as ith elem   elements   ---    Number of eleemnts
    size       ---    Number of elements

      -     rv1,rv2 Returns rv1-rv2                                        

    Function/Operator                         Result
    -----------------   ------------------------------------------------------
          =             Current row vector set equal to input rvec
    (int,int)           Reference to element <i|hmx|j> (Potential Danger)
    get(int,int)        Copy of element <i|dmx|j>      (Safe)
    put(int,int)        Assigns element <i|dmx|j>      (Return FALSE if fails)
    put_h(int,int)      Assigns both <i|dmx|j> & <j|dmx|i> (Return TRUE/FALSE)
    get_block(r,c,R,C)  Returns dmx block of size RxC starting at <r|dmx|c>
    put_block(r,c,mx)   Places mx into dmx at position <r|dmx|c> (TRUE/FALSE)

    Again, note that operations which return TF will trigger class matrix to
    take appropriate action when they return FALSE.  This usually means that
    the operation cannot be handled by a diagonal array so the array type
    is changed to be more generic.  Then the operation is reattempted.      */

void row_vector::operator = (const row_vector& rvec)
  { matrix::operator= ( (matrix&)rvec ); }
void row_vector::operator = (const col_vector& cvec)
  { RVerror(31,1); RVfatality(29); put(cvec.get(0),0); }
void row_vector::operator = (const matrix& mx)
  {
  if(mx.rows() != 1)
    {
    std::string em = std::string(MxModdec(mx.rows()))
              + std::string("x")
              + std::string(MxModdec(mx.cols()));
    RVerror(20,em,1);			// Construct multirow array
    RVfatality(40);			// Can't do construction
    }
  matrix::operator= ( (matrix&)mx );
  }

// *********************** row vector access functions ************************

complex& row_vector::operator() (int i) { return this->matrix::operator()(0,i); };
complex  row_vector::get(int i)   const { return matrix::get(0,i); }
double   row_vector::getRe(int i) const { return matrix::getRe(0,i); }
double   row_vector::getIm(int i) const { return matrix::getIm(0,i); }
void     row_vector::put(const complex& z, int i) { matrix::put(z,0,i); }
int      row_vector::elements()   const { return matrix::cols(); }
int      row_vector::size()       const { return cols(); }

// ____________________________________________________________________________
// F                  ROW VECTOR BINARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________

/* Note: The operators which must be defined here, as opposed to those
         inherited from class matrix, are they which return row_vector.
         (or at least, return a type other than what class matrix does)
         Definition  herein also allows for more explicit error messages

   Note: 1.) rvec +/- cvec is not defined, only works when (1x1) +/- (1x1) so
                           return is left as 1x1 matrix, error matrix error  
  

   Operator Arguments      Result        Operator Arguments    Result
   -------- --------- -----------------  -------- --------- -------------------
      +     rv,rv1    Returns rv+rv1        +     rv,mx     Returns rv+mx
      -     rv,rv1    Returns rv-rv1        -     rv,mx     Returns rv-mx
      -        rv    Returns -rv            -=    rv,rv1  rv1 subt. from rv
     +=     rv,rv1  rv1 added to rv         *=    rv,rv1  rv mult into rv1
      +     rv,rv1    Returns rv+rv1         *    rv1,rv2 Returns rv1*rv2
      -     rv1,rv2 Returns rv1-rv2                                          */

row_vector operator * (const complex& z, const row_vector& rvec)
  { return row_vector(z * (const matrix &)rvec); }

row_vector operator * (double d, const row_vector& rvec)
  { return row_vector(d * (const matrix &)rvec); }

row_vector row_vector::operator + (const row_vector& rvec) const
  {
  if(size() != rvec.size())
    {
    RVerror(42,1);			// Matrix-Vector size trouble
    RVfatality(21,"addition");		// Problems in addition
    }
  return row_vector((const matrix&)(*this) + (const matrix&)rvec);
  }

row_vector row_vector::operator + (const matrix& mx) const
  {
  if(size() != mx.cols() || mx.rows() != 1)
    {
    RVerror(41,1);				// Vector vector size trouble
    RVfatality(21,"addition");			// Problems in addition
    }
  return row_vector((const matrix&)(*this) + mx);
  }

row_vector row_vector::operator - (const row_vector& rvec) const
  {
  if(size() != rvec.size())
    {
    RVerror(42,1);				// Matrix-Vector size trouble
    RVfatality(21,"subtraction");		// Problems in subtraction
    }
  return row_vector((const matrix &)(*this) - (const matrix &)rvec);
  }

row_vector row_vector::operator - (const matrix& mx) const
  {
  if(size() != mx.cols() || mx.rows() != 1)
    {
    RVerror(41,1);				// Vector vector size trouble
    RVfatality(21,"subtraction");		// Problems in subtraction
    }
  return row_vector((const matrix &)(*this) - mx);
  }

complex row_vector::operator * (const row_vector& rvec) const
  {
  RVerror(43, 1);                       // Cannot multiply 2 col vects
  std::string pname("Row Vector * Row Vector");
  RVfatality(21,pname);
  return complex0+rvec.get(0);
  }
 
row_vector row_vector::operator * (const matrix& mx) const
  { 
  if(size() != mx.rows())
    {
    RVerror(42,1);				// Matrix-Vector size trouble
    RVfatality(21,"multiplication");		// Problems in multiplication
    }
  return row_vector((const matrix &)(*this) * mx);
  }

complex row_vector::operator * (const col_vector& cvec) const
  { 
  if(size() != cvec.size())
    {
    RVerror(41,1);				// Vector vector size trouble
    RVfatality(21,"multiplication");		// Problems in multiplication
    }
  return (((const matrix &)*this * (const matrix&)cvec).get(0,0));
  }

row_vector row_vector::operator * (const complex& z) const
  { return row_vector((const matrix &)(*this)*z); }

row_vector row_vector::operator * (double d) const
  { return row_vector((const matrix &)(*this)*d); }

row_vector row_vector::operator / (const complex& z) const
  { return row_vector((const matrix &)*this/z); }

row_vector row_vector::operator / (double d) const
  { return row_vector((const matrix &)*this/d); }


// ____________________________________________________________________________
// G               ROW VECTOR UNARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________

/* Note: These operators, as they are unable to change the type of the
  	 result, can only be those which return row_vector. 

   Operator Arguments      Result        Operator Arguments    Result
   -------- --------- -----------------  -------- --------- -------------------
      +=    rv,rv1    rv1 added to rv       +=    rv,mx     mx added to rv
      *=    rv,z      rv1 mult. by z        /=    rv,z      rv divided by z  */

row_vector& row_vector::operator += (const row_vector& rvec)
{
  if(!rvec.size()) 
    return *this;
  if(!size())
  {
    matrix::operator= ( (const matrix&)rvec );
    return *this;
  }
  if(size() != rvec.size())
  {
    RVerror(42, 1);                             // Vec-Vec dimensions bad
    std::string pname("Row Vector += Row Vector");
    RVfatality(21,pname);
  }
  matrix::operator+= ((const matrix &)rvec);
  
  return *this;
}

void row_vector::operator += (const matrix& mx)
  {
  if(!mx.rows()) return;                        // Nothing if mx NULL
  if(size() != mx.cols() || mx.rows() != 1)
   {
   RVerror(41, 1);   	                  // Vec-Mx dimensions bad
   std::string pname("Row Vector += Matrix");
   RVfatality(21,pname);
   }
  matrix::operator+= (mx);
  }

row_vector& row_vector::operator -= (const row_vector& rvec)
{
  if(!rvec.size()) 
    return *this;			   // Nothing if empty rvec
  if(!size())                                // Set -rvec if we are empty
  { 
    (*this) = rvec.matrix::operator-();
    return *this; 
  }
  if(size() != rvec.size())
  {
    RVerror(42, 1);                     // Vec-Vec dimensions bad
    std::string pname("Row Vector -= Row Vector");
    RVfatality(21,pname);
  }
  matrix::operator-= ((const matrix &)rvec);
  return *this;
}

void row_vector::operator -= (const matrix& mx)
  {
  if(!mx.rows()) return;                        // Nothing if mx NULL
  if(size() != mx.cols() || mx.rows() != 1)
    {
    RVerror(41, 1);                     // Vec-Mx dimensions bad
    std::string pname("Row Vector -= Matrix");
    RVfatality(21,pname);
    }
  matrix::operator-= (mx);
  }

row_vector& row_vector::operator *= (      double   d) 
{ 
  matrix::operator*= (d); 
  return *this;
}
row_vector& row_vector::operator *= (const complex& z) 
{ 
  matrix::operator*= (z); 
  return *this;
}
row_vector& row_vector::operator /= (      double   d) 
{ 
  matrix::operator/= (d); 
  return *this;
}
row_vector& row_vector::operator /= (const complex& z) 
{ 
  matrix::operator/= (z); 
  return *this;
}

// ____________________________________________________________________________
// J                    ROW VECTOR SIMPLE UNARY FUNCTIONS
// ____________________________________________________________________________

/* Note: The negation operator must be defined herein.  The function conjugate
         is inherited, but adjoint and transpose return col_vectors so they
         are also here. Transpose included just for specific error messaging.*/


        // Input		rvec : A row vector
        // Output adjoint	cvec : A column vector of the adjoint of rvec
        //                                                t        * T
        //                                 <cvec| = |rvec>  = |rvec >
        // Output transpose     cvec : A column vector of the transpose of rvec
        //                                               T
        //                                <cvec| = |rvec>
        // Output		void : Trace is disallowed on vectors
 
// row_vector row_vector::operator- () const;		// MATRIX INHERITED

col_vector adjoint(const row_vector& rvec)
  { return col_vector(adjoint((const matrix &)rvec)); }

col_vector transpose(const row_vector& rvec)
  { return col_vector(transpose((const matrix &)rvec)); }
 
complex trace(const row_vector& rvec)
  {
  rvec.RVfatality(44);                  // Vec-Mx dimensions bad
  return rvec.get(0);                   // Compiler likes cvec to be used
  }

row_vector row_vector::differential() const
  {
  int nr = size();
  row_vector rvec(nr);
  for(int i=0; i<nr-1; i++)
    rvec.put((*this).get(i+1) - (*this).get(i), i);
  rvec.put((rvec)(nr-2), nr-1);
  return rvec;
  }

// ____________________________________________________________________________
// I                    ROW VECTOR SIMPLE BINARY FUNCTIONS
// ____________________________________________________________________________

/*  Function  Output                       Description
    --------  ------  ---------------------------------------------------------
      norm    double  Returns the vector norm, sq.rt. of sum of squares
      sum     complex Returns the sum ov all vector elements
     maxRe    double  Returns largest real value in the vector
     maxIm    double  Returns largest imaginary value in the vector
     maxZ     complex Returns largest complex (biggest norm) value in vector
     minRe    double  Returns smallest real value in the vector
     minIm    double  Returns smallest imaginary value in the vector
     minZ     complex Returns smallest complex (biggest norm) value in vector
      max      int    Returns index of largest (  real, imag, norm } element
      min      int    Returns index of smallest { real, imag, norm } element
     flip     void    Inverses the element order in the vector
     zero     void    Sets all elements to zero
     sort     vector  Returns integer array of sorted indices

  The vector norm is given by 

                     [ size-1          ]
                     |  ---            |
                     |  \            2 |            2                *
          Norm = sqrt|  /  |<rvec|i>|  |, |<rvec|i>| = <rvec|i>*<rvec |i>
                     |  ---            |
                     [  i=0            ]

  & the following definitions apply to the max, min, and sort functions:
       
                          type         max/min
                            0       norm (default)
                           <0         imaginary                         
                           >0           real                                 */

double row_vector::norm() const
  {
  int n = size();			// Number of vector elements
  double d=0;				// Initialize norm to zero
  for(int i=0; i<n; i++)		// Take sum of the squares
    d += square_norm(get(i));
  return sqrt(d);			// Now return the square root
  }

complex row_vector::sum() const 
  {
  complex z(0);
  int n = size();                       // Number of vector elements
  for(int i=0; i<n; i++) z += get(i);	// Sum them all up
  return z;
  }

complex row_vector::sum(int st, int ne) const 
  {
  complex z(0);
  for(int i=st; i<st+ne; i++)
    z += get(i);
  return z;
  }

double  row_vector::maxRe() const { return getRe(max(1));  }
double  row_vector::minRe() const { return getRe(min(1));  }
double  row_vector::maxIm() const { return getIm(max(-1)); }
double  row_vector::minIm() const { return getIm(min(-1)); }
complex row_vector::maxZ()  const { return get(max(0));    }
complex row_vector::minZ()  const { return get(min(0));    }

int row_vector::max(int type) const
  {
  int n    = size();				// Number of vector elements
  int imax = 0;					// Index of maximum
  double max=-HUGE_VAL, tmp=0;			// Begin with very small max
  if(type>1)      type=1;			// If type>=1, maximum of reals
  else if(type<0) type=2;			// If type<0,  maximum of imags
  int i;					// Index of maximum value
  switch(type)
    {
    case 0:					// Index of max. norm value 
      for(i=0; i<n; i++)			// Loop over all points
        {					// & compare values
        tmp = square_norm(get(i));		//   Re*Re + Im*Im at i
        if(tmp > max) { max = tmp; imax = i; }	//   If largest value
        }					//   set as largest, store i
      break;
    case 1:					// Index of max. real value
      for(i=0; i<n; i++)			// Loop over all points
        {					// & compare values
        tmp = getRe(i);				//    Real value at i
        if(tmp > max) { max=tmp; imax=i; }	//    If largest value
        }					//    set as largest, store i
      break;
    case 2:					// Index of max. imag value
      for(i=0; i<n; i++)			// Loop over all points
        {					// & compare values
        tmp = getIm(i);				//    Imaginary value at i
        if(tmp > max) { max = tmp; imax = i; }	//    If largest value
        }					//    set as largest, store i
      break;
    }
  return imax;					// Return index of maximum
  }

int row_vector::min(int type) const
  {
  int n    = size();				// Number of vector elements
  int imin = 0;					// Index of minimum
  double min=HUGE_VAL, tmp=0;			// Begin with very larg min
  if(type>1)      type=1;			// If type>=1, minimum of reals
  else if(type<0) type=2;			// If type<0,  minimum of imags
  int i;					// Index of minimum value
  switch(type)
    {
    case 0:					// Index of min. norm value 
      for(i=0; i<n; i++)			// Loop over all points
        {					// & compare values
        tmp = square_norm(get(i));		//   Re*Re + Im*Im at i
        if(tmp < min) { min=tmp; imin=i; }	//   If smallest value
        }					//   set as smallest, store i
      break;
    case 1:					// Index of min. real value
      for(i=0; i<n; i++)			// Loop over all points
        {					// & compare values
        tmp = getRe(i);				//    Real value at i
        if(tmp < min) { min=tmp; imin=i; }	//    If smallest value
        }					//    set as smallest, store i
      break;
    case 2:					// Index of min. imag value
      for(i=0; i<n; i++)			// Loop over all points
        {					// & compare values
        tmp = getIm(i);				//    Imaginary value at i
        if(tmp < min) { min=tmp; imin=i; }	//    If smallest value
        }					//    set as smallest, store i
      break;
    }
  return imin;
  }

void row_vector::flip()
  {
  int n = size();			// Number of vector elements
  complex ztmp(0);			// Temporary value
  int i, j;
  for(i=0, j=n-1; i<n/2; i++,j--)	// Loop over the 1st half
    {
    ztmp += get(i);			// Store ith element
    put(get(j), i);			// Set ith elem to (n-1-i)th element
    put(ztmp, j);			// Set (n-1-i)th elem to ith elem
    }
  }

void row_vector::zero()
  { for(int i=0; i<size(); i++) put(complex0, i); }

std::vector<int> row_vector::sort(int type) const
  {
  int      nvals = size();		// Number of values to be sorted
  double *vals;				// Array of values to be sorted
  vals = new double[nvals];		// Allocate the array of values
  std::vector<int> indx;			// Allocate vector to return
  if(type < 0)  type =-1;
  if(type > 0)  type = 1;		// Insure type = { -1, 0, 1 }

//		Copy Column Or Row That Will Be Sorting Key

  int i=0,j=0;				// Array indices we'll use
  for(i=0; i<nvals; i++)		// First fill array vals with values 
    { 					// we want to sort & initialize the
    indx.push_back(i); 			// vector of sorted indices we return
    switch(type)
      {
      default:
      case 1:  vals[i] = getRe(i);     break;	// Sorting reals
      case 0:  vals[i] = square_norm(get(i)); break;	// Sorting norms
      case -1: vals[i] = getIm(i);   break;	// Sorting imaginaries
      }
    }

//	   Sort Column Or Row And Put Sorted Indices Into Return Vector

  double maxval;			// Value for maximum
  int maxind;				// Value of index at maximum
  for(i=0; i<nvals; i++)		// Now begin the sorting
    {					// by looping over the values
    maxval = vals[i];			// Assume this is maximum
    maxind = i;				// This is it's index
    for(j=i+1; j<nvals; j++)		// Now compare with all others
      {
      if(vals[j] > maxval)		// 	If we find a larger one
        {				//	then we set it as the maximum
        maxval = vals[j];		//	--> Store the bigger value
        maxind = j;			//	--> Store the big value index
        }
      }
    vals[maxind] = vals[i];		// Move value we started with to max
    j = indx[maxind];			// Here is the index of largest value
    indx[maxind] = indx[i];		// Copy current index to largest spot
    indx[i] = j;			// Copy index of largest to current spot
    }
  delete [] vals;			// Delete the values
  return indx;				// Return vector of indices
  }

// ____________________________________________________________________________
// J                 ROW VECTOR COMPLEX UNARY FUNCTIONS
// ____________________________________________________________________________

// ******************************      FFT      *******************************

row_vector FFT(const  row_vector& rvec) { return FFT((const  matrix &)rvec); }
row_vector IFFT(const row_vector& rvec) { return IFFT((const matrix &)rvec); }
 
        // Input        rvec : A row vector
        // Output      rvec1 : A row vector containing the discrete Fourier
        //                     (inverse) transform of the values in rvec
        // Note              : The vector rvec should be of dimension
        //                     1 x (2^n)
// ____________________________________________________________________________
// K                    ROW VECTOR CONTAINER SUPPPORT
// ____________________________________________________________________________

        // Input      rvec1  : A row vector
        //            rvec2  : A row vector
        // Output ==    T/F  : Returns TRUE if cvec1 is equal to cvec2.

//bool row_vector::operator==(const row_vector& rvec) const	// INHERITED
//bool row_vector::operator!=(const row_vector& rvec) const	// INHERITED
//bool row_vector::operator> (const row_vector& rvec) const	// INHERITED
//bool row_vector::operator< (const row_vector& rvec) const	// INHERITED

// ____________________________________________________________________________
// J                ROW VECTOR COMPLEX BINARY FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                              Straight Product
// ----------------------------------------------------------------------------

        // Input              rvec   : A row vector
        // Output             rvec1  : A "product" of rvec & itself
        // Note                      : This is NOT a scalar product
        //                             rvec1(i) = rvec(i) * rvec(i)

row_vector row_vector::product() const
  {
  int dim = size();
  row_vector rvec3(dim);
  for(int i=dim-1; i>=0; i--)
    rvec3.put(get(i)*get(i), i);
  return rvec3;
  }

        // Input              rvec1  : A row vector
        //                    rvec2  : Second row vector
        // Output             rvec3  : A "product" of rvec1 & rvec2
        // Note                      : This is NOT a scalar product
        //                             rvec3(i) = rvec1(i) * rvec2(i)


row_vector product(const row_vector& rvec, const col_vector& cvec)
  { return rvec.product(cvec); }

row_vector product(const row_vector& rvec1, const row_vector& rvec2)
  { return rvec1.product(rvec2); }

row_vector row_vector::product(const col_vector& cvec) const
  {
  int dim = size();
  if(dim != cvec.size()) 
    {
    RVerror(42,1);                              // Vec-vec size trouble
    RVerror(21,"Row Vector (x) Col Vec");       // Problems in this
    std::string pname("product");
    RVfatality(5, pname);
    }
  row_vector rvec3(dim);
  for(int i=dim-1; i>=0; i--)
    rvec3.put(get(i)*cvec.get(i), i);
  return rvec3;
  }

row_vector row_vector::product(const row_vector& rvec) const
  {
  int dim = size();
  if(dim != rvec.size()) 			//Size mismatch
    {
    RVerror(42,1);                              // Vec-vec size trouble
    RVerror(21,"Row Vector (x) Row Vec");       // Problems in this
    std::string pname("product");
    RVfatality(5, pname);
    }
  row_vector rvec3(dim);
  for(int i=dim-1; i>=0; i--)
    rvec3.put(get(i)*rvec.get(i), i);
  return rvec3;
  }

// ----------------------------------------------------------------------------
//                              Scalar Product
// ----------------------------------------------------------------------------

        // Input              cvec1  : A column vector
        //               cvec2,rvec  : Second column vector, or row vector
        // Output             z      : Scalar product of adjoint(cvec) & cvec

        //                             z = <cvec'|cvec>,   <cvec'| = adjoint(|rvec>)
        //                               = <cvec1'|cvec2>, <cvec2'| = adjoint(|rvec>)
        //                               = <rvec|cvec>
        // Input              rvec  : A row vector
        // Output             z      : Scalar product of rvec & adjoint(rvec)

        //                             z = <rvec|rvec'>,   |rvec'> = adjoint(<rvec|)

double scalar_product(const row_vector& rvec)
  { return rvec.scalar_product(); }

double row_vector::scalar_product() const
  {
  int dim = size();
  double sp=0;
  for(int i=dim-1; i>=0; i--)
    sp += square_norm(get(i));
  return sp;
  }

        // Input              rvec1  : A row vector
        //                    rvec2  : Second row vector
        // Output             z      : Scalar product of rvec1 & adjoint(rvec2)

        //                             z = <rvec1|rvec2'>, |rvec2'> = adjoint(<rvec|)
        // Input              rvec  : A row vector
        //                    cvec   : A column vector
        // Output             z      : Scalar product of rvec & cvec

	//			       z = <rvec|cvec>


complex scalar_product(const row_vector& rvec1, const row_vector& rvec2)
  { return rvec1.scalar_product(rvec2); }

complex scalar_product(const row_vector& rvec, const col_vector& cvec)
  { return rvec.scalar_product(cvec); }

complex row_vector::scalar_product(const row_vector& rvec) const
  {
  int dim = size();
  if(dim != rvec.size())
    {
    RVerror(42, 1);                             // Vec-Vec dimensions
    RVerror(21,"Row Vector (x) Row Vec");       // Problems in this
    std::string pname("scalar_product");
    RVfatality(5, pname);
    }
  complex z(0);
  for(int i=dim-1; i>=0; i--)
    z += (get(i)).conj_times(rvec.get(i));
  return z;
  }

complex row_vector::scalar_product(const col_vector& cvec) const
  {
  int dim = size();
  if(dim != cvec.size())
    {
    RVerror(42, 1);                             // Vec-Vec dimensions bad
    RVerror(21,"Row Vector (x) Col Vec");       // Problems in this
    std::string pname("scalar_product");
    RVfatality(5, pname);
    }
  complex z(0);
  for(int i=dim-1; i>=0; i--)
    z += get(i)*cvec.get(i);
  return z;
  }


// ____________________________________________________________________________
// L                      ROW VECTOR I/O FUNCTIONS
// ____________________________________________________________________________

// ------------------------ ASCII Output Functions ----------------------------

/*              Input		rvec	: A row vector (this)
                                ostr	: Output ASCII file stream
                                full	: Flag for amount of output
                Return          void	: rvec is sent to the output stream  */

std::string row_vector::hdrString() const
  {
  std::string hdr;
  if(!size()) hdr = "Empty Row Vector";
  else
    {
    std::string vtype("");
    if(is_real())           vtype = " Real";
    else if(is_imaginary()) vtype = " Imaginary";
    hdr = std::string("1 x ") + MxModdec(size()) 
        + vtype + std::string(" Row Vector");
    }
  return hdr;
  }

/* The row vector is placed into an array of strings, each string having
   equal length. If the number of vector elements exceeds the number of
   elements we allow on a single printed line, then a new line is begun.
   The number of elements per line is set by PFlgs.VxCols.                  */

std::vector<std::string> row_vector::printStrings(const MxPrint& PFlgs) const
  {
  std::vector<std::string> PStrings;
  int npts = size();
  if(!npts) return PStrings;
 
  int ptype = 0;				// Assume complex elements
  if(PFlgs.MxRIPrnt)				// See if we want just reals
    {						// or just imags printed
    if(is_real())           ptype = 1;		//   Real elements output
    else if(is_imaginary()) ptype = 2; 		//   Imag elements output
    }

  int elen;					// A single element length
  switch(ptype)					// Get the element length from
    {						// class complex. Depend on 
    default:					// complex class settings
    case 0: elen = complex::zlength(); break;	//   Element length if complex
    case 1:					//   Element length if real
    case 2: elen = complex::dlength(); break;	//   Element length if imaginary
    }
  std::string ezer = std::string(elen, ' ');	// Invisible element

  std::string efmt = complex::dformat();	// Real/Imag element format
  std::string vline;				// String for printed line

  int ncol = PFlgs.VxCols;			// Columns in output
  if(ptype) ncol *=2;				// Twice as many if real/imag
  int col=0, row=0;				// Column index
  for(int pos=0; pos<npts; pos++ )
    {
    if(col == 0) vline = "";
    if(ptype==1)					// Output as 
      vline += MxModform(efmt.c_str(),getRe(pos)); 	// a real element
    else if(ptype==2)					// Output as
      vline += MxModform(efmt.c_str(),getIm(pos)); 	// an imaginary element
    else  vline += get(pos).printString();         	// Output as complex
    col++;
    if(col < ncol)  vline += " ";
    else          { PStrings.push_back(vline); col=0; row++; }
    }
  if(col)					// If we have not output all
    {						// columns yet
    if(row)					//    Fill with blanks if
      while(col < ncol)				//    previous rows output
       { vline += ezer + " "; col++; }
    PStrings.push_back(vline);			//    Output last row
    }
  return PStrings;
  }

std::ostream& row_vector::printcols(std::ostream& ostr, int ncol, int rc, int nelem) const
  {
  if(!size()) return ostr;
  std::string marg(10, ' ');
  std::string spc(3, ' ');
  int col=0;
  if(nelem < 1) nelem = size(); 
  for(int i=0; i<nelem; i++ )
    {
    if(col == 0) ostr  << marg;
    if(!rc)       ostr << getRe(i);
    else if(rc<1) ostr << getIm(i);
    else          ostr << get(i);
    col++;
    if(col >= ncol) { ostr << "\n"; col=0; }
    else              ostr << spc;
    }
  return ostr;
  }

std::ostream& row_vector::print(std::ostream& ostr, int full) const
  {
  if(!size()) { ostr << "Empty Row Vector\n"; return ostr; }
  ostr << "GAMMA " << "1 x " << size() << " Row Vector";
  if(full)
    {
    ostr << "\n\n";
    int ncol = 4;
    int col=0;
    std::string marg(10, ' ');
    for(int pos=0; pos<size(); pos++ )
      {
      if(col == 0) ostr << marg;
      ostr << get(pos) << "  ";
      col++;
      if(col >= ncol) { ostr << "\n"; col=0; }
      }
    }
  else
    {
    int i = max(1);
    ostr << "\n\tMaximum Real At Point " << i << ", " << getRe(i);
    i = max(-1);
    ostr << "\n\tMaximum Imag At Point " << i << ", " << getIm(i);
    i = max(0);
    ostr << "\n\tMaximum Norm At Point " << i << ", " << sqrt(square_norm(get(i)));
    ostr << "\n\tVector Integral Is    " << sum();
    ostr << "\n\tVector Norm Is        " << norm();
    }
  ostr << "\n";
  return ostr;
  }


std::ostream& operator << (std::ostream& ostr, const row_vector& rvec) 
  { return rvec.print(ostr,1); }

std::istream& operator >> (std::istream& istr, row_vector& rvec)
  {
  matrix t1;
  istr >> t1;
  rvec = row_vector(t1);
  return istr;
  }



// code from strproc.cc

/**
 * trim_left - trim all whitespace from left side of a string
 * @in: string to be processed. This string is not modified; a
 * modified copy of the string is returned instead.
 */

string row_vector::trim_left (const string& in)
{
  int n = 0;
  string::const_iterator p = in.begin();

  if (in.length() == 0) {
    string out = "";
    return out;
  }

  while (p != in.end()) {
    if (isspace (*p))
      n++;
    else
      break;
    p++;
  }
  string out (in, n, in.length());
  return out;
}

/**
 * trim_right - trim all whitespace from right side of a string
 * @in: string to be processed. This string is not modified; a
 * modified copy of the string is returned instead.
 */
string row_vector::trim_right (const string& in)
{
  int n = 0;
  string::const_reverse_iterator p = in.rbegin();

  while (p != in.rend()) {
    if (isspace (*p))
      n++;
    else
      break;
    p++;
  }
  string out (in, 0, in.length() - n);
  return out;
}

/**
 * trim_all - remove whitespace from left and right sides of a string
 * @in: the string to be processed
 */
string row_vector::trim_all (const string& in)
{
  if (in.length() == 0) {
    string out = "";
    return out;
  }

  string out = trim_left (in);
  out = trim_right (out);
  return out;
}

/**
 * isws - use a restricted set of whitespace characters to classify @c
 * @c: character to be analyzed
 * Description: isws() returns 1 if @c is a space (` ') or tab (`\t')
 * character, and 0 otherwise.
 */
int row_vector::isws (const char c)
{
  if (c == ' ' || c == '\t')
	return 1;
  else
	return 0;
}

/**
 * squeeze - remove @line's excess separation characters and leading whitespace
 * @line: the string to be processed
 * Description: squeeze() returns a copy of @line with exactly one space
 * character (` ') between words. Excess space characters and tabs (`\t')
 * are removed in the returned copy. squeeze() also removes leading
 * whitespace from @line.
 */
string row_vector::squeeze (string line)
{
  int nwhite = 0;
  string squeezed, whchars = " \t";

  if (!line.length())
    return "";

  int pos = line.find_first_not_of (whchars);

  for (unsigned int i = pos; i < line.length(); i++) {
    if (isws (line[i]) && nwhite < 1) {    /* Keep one whitespace char.  */
	  squeezed += " ";
	  nwhite++;
	}
	else if ((isws (line[i])) && nwhite >= 1) {
	  nwhite++;
	}
	else {
	  squeezed += line[i];
	  nwhite = 0;
	}
  }
  return squeezed;
}

/**
 * split - split a string into tokens delimited by @sep
 * @s: a string holding a line of text to be split into
 * fields separated by @sep.
 * @sep: the character used to separate fields in line.
 * @field: a vector of strings to hold the tokenized version of @s.
 * Return Value: The number of fields in @s.
 */
int row_vector::split (string s, char sep, vector<string>& field)
{
  unsigned int i, nfield;
  unsigned int j = 0;
  string fld;

  if ((i = s.length()) == 0)
    return 0;

  i = nfield = 0;
  do {
    if (i < s.length() && s[i] == sep) { /* skip sep  */
      ++i;
    }
    else {
      j = s.find_first_of( sep, i );
      if (j > s.length())
	    j = s.length();
      fld = string( s, i, j-i );
    }
    if (nfield >= field.size())
      field.push_back( fld );
    else
      field[nfield] = fld;
    nfield++;
    i = j + 1;
  } while (j < s.length());

  return nfield;
}



//------------- This section imported from readpulse.cc --------------

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef DEG2RAD
#define DEG2RAD (PI/180.0)
#endif


/*******************************************************************
 * ReadPulse - read pulses from various forms of ASCII files:
 *
 * @filename: file name of the pulse profile file
 *
 * @PulseFmt: enumerated type of the pulse, corresponding to the
 * type of hardware the pulse was written by/for.
 *
 *******************************************************************/

//int ReadPulse (const char *filename, int PulseFmt, row_vector& v)
MSVCDLL row_vector row_vector::read_pulse (const string filename, const int PulseFmt)
{
  ifstream fin (filename.c_str(), ios_base::in);
  if (!fin) {
    //eprintf ("could not open pulse file %s", filename);
  }

  /* Decide what kind of file we have.  */
  switch (PulseFmt) {
	  case SMIS:
		//weprintf ("Unknown pulse file type; will use ideal pulses");
		break;
	  case SIEMENS:  /* Siemens file with header  */
		return row_vector::ReadSiemens (fin);
		break;
	  case SIEMENS_NOHDR:   /* Siemens file without header  */
		return row_vector::ReadSiemens_Nohdr (fin);
		break;
	  case SVS:
		return row_vector::ReadSVS (fin);
		break;
	  case PLAIN_ASCII:
		return row_vector::Read_Plain_ASCII (fin);
		//eprintf ("cannot handle PLAIN_ASCII type");
		break;
	  case ASCII_MT_DEG:
		return row_vector::Read_ASCII_mT_Deg (fin);
		break;
	  default:
		//eprintf ("unknown pulse type");
		break;
  }

  fin.close();

  return row_vector();
}


/*
 * ReadSiemens - read a Siemens pulse profile
 *
 * @fin: file stream of the pulse file
*/
row_vector row_vector::ReadSiemens (ifstream& fin)
{
  vector<string> fields;
  vector<double> pts;
  int i = 0, count = 0, in_data = 0;
  string line;

  while (getline (fin, line)) 
  {

    if (line.find ("Entry_Values") < line.length())
      in_data = 1;

    if (in_data) 
    {
      fields.clear();
      line = row_vector::squeeze (line);
      line = row_vector::trim_all (line);
      i = row_vector::split (line, ' ', fields);
      for (unsigned int j = 0; j < fields.size(); j++) 
      {
		if (is_decimal (fields[j].c_str())) 
        {
		  pts.push_back (atof (fields[j].c_str()));
		  count++;
		}
      }

      if (line.find ("End_Entry") < line.length())
	    in_data = 0;
    }
  }

  if (count % 2 != 0) 
  {
    //eprintf ("unequal number of amplitudes and phases");
  }
  else 
  {
    int npts = pts.size() / 2;
    row_vector v (npts, 0.0);
    for (int i = 0; i < npts; i++) 
    {
      complex ij = complex_exp(-complexi * DEG2RAD * pts[i + npts]);
      v.put (pts[i] * ij, i);
    }
    return v;
  }

  row_vector nv;
  return nv;      /* never get here */
}

/**
 * ReadSiemens_Nohdr - read and process a Siemens pulse file without header
 * @fin: file stream of the open pulse file
 */
row_vector row_vector::ReadSiemens_Nohdr (ifstream& fin)
{
  int pos, numpts, count = 0;
  string line, cooked;
  vector<string> field;
  vector<double> allpts;

  while (fin) {
    getline (fin, line);
    if (line.find ("End_Entry:") < line.length())
      break;

    cooked = row_vector::squeeze (line);
    pos    = row_vector::split (cooked, ' ', field);

    for (unsigned int i = 0; i < field.size(); i++) {
      string snum = field[i];
      const char* cptr = snum.data();
      /* No way to know where the amplitudes end and the
       * phases begin.  */
      allpts.push_back (atof (cptr));
      count++;
    }
    /* Now erase the contents of field, or else they accumulate.  */
    field.clear();
  }

  numpts = count / 2;
  vector<double> amps (numpts), phases (numpts);
  row_vector rvec (numpts);
  for (int i = 0; i < numpts; i++) {
    amps[i]   = allpts[i*2];			//bjs-May2002
    phases[i] = allpts[i*2 + 1];
  }

  for (int i = 0; i < numpts; i++) {				//bjs-Aug2001
	rvec.put ( amps[i] * complex_exp (-complexi * DEG2RAD * phases[i]), i);
  }

  return rvec;

}

/**
 * ReadSVS - Read in and process a complete SVS file
 * @fin - file stream of the open pulse file
 */
row_vector row_vector::ReadSVS (ifstream& fin)
{
  int pos, numpts = 0, in_data = 0;
  string line, cooked;
  double apt;
  vector<double> pts;
  vector<string> field;
  vector<double> amps, phases;

  const char* SIEM_HDR = "Num_Points:";
  const char* SIEM_VAL = "Entry_Values:";

  while (getline (fin, line)) {
    cooked = row_vector::trim_all (line);
    cooked = row_vector::squeeze (line);

    if (cooked.find ("End_Entry:") < cooked.length()) {
      in_data = 0;
      break;
    }

    if (cooked.find (SIEM_HDR) < cooked.length()) {
      /* Extract number of points in file  */
      pos = row_vector::split (cooked, ' ', field);
      string snum = field[1];
      const char* cptr = snum.data();
      numpts = atoi (cptr);
    }

    if (cooked.find (SIEM_VAL) < cooked.length()) {
      pos = row_vector::split (cooked, ' ', field);
      /* Excise SIEM_VAL from cooked string.  */
      pos = cooked.find_first_of (" ");
      cooked = cooked.substr (pos + 1);
      /* Now erase the contents of field, or else they accumulate.  */
      field.clear();
      in_data = 1;
    }

    if (in_data) {
      while (fin >> apt) {
	pts.push_back (apt);
      }
    }

  }

  row_vector rvec ( numpts );
  for (unsigned int i = 0; i < pts.size()/2; i++) {
    rvec.put (complex(pts[i], 0.0), i);
  }

  return rvec;
}


/**
 * ReadPlain_ASCII - read and process a Plain ASCII pulse file (Ampl(Hz),Phase(degr))
 * @fin: file stream of the open pulse file
 */
row_vector row_vector::Read_Plain_ASCII (ifstream& fin)
{
  int pos, numpts, count = 0;
  string line, cooked;
  vector<string> field;
  vector<double> allpts;

  while (fin) {
    getline (fin, line);
    if (line.find ("End_Entry:") < line.length())
      break;

    cooked = row_vector::squeeze (line);
    pos    = row_vector::split (cooked, ' ', field);

    for (unsigned int i = 0; i < field.size(); i++) {
      string snum = field[i];
      const char* cptr = snum.data();
      /* No way to know where the amplitudes end and the
       * phases begin.  */
      allpts.push_back( atof(cptr) );
      count++;
    }
    /* Now erase the contents of field, or else they accumulate.  */
    field.clear();
  }

  numpts = count / 2;
  
  row_vector rvec (numpts);
  for (int i = 0; i < numpts; i++) {
    rvec.put( complex( allpts[i*2],allpts[i*2+1] ), i );
  }
  
  return rvec;

}


/**
 * ReadPlain_ASCII_mT_Deg - read and process a Plain ASCII pulse file (Ampl(mT),Phase(degr))
 * @fin: file stream of the open pulse file
 */
row_vector row_vector::Read_ASCII_mT_Deg (ifstream& fin)
{
  int pos, numpts, count = 0;
  string line, cooked;
  vector<string> field;
  vector<double> allpts;

  while (fin) {
    getline (fin, line);
    if (line.find ("End_Entry:") < line.length())
      break;

    cooked = row_vector::squeeze (line);
    pos    = row_vector::split (cooked, ' ', field);

    for (unsigned int i = 0; i < field.size(); i++) {
      string snum = field[i];
      const char* cptr = snum.data();
      /* No way to know where the amplitudes end and the
       * phases begin.  */
      allpts.push_back( atof(cptr) );
      count++;
    }
    /* Now erase the contents of field, or else they accumulate.  */
    field.clear();
  }

  numpts = count / 2;
  
  row_vector rvec (numpts);
  for (int i = 0; i < numpts; i++) {
    rvec.put( complex( allpts[i*2] * 42600.0, allpts[i*2+1] ), i );
  }
  
  return rvec;

}


/**
 * is_decimal - decides whether a string represents a real number
 * @cp: string containing the suspect number
 */
int row_vector::is_decimal (const char *cp)
{
  /*
   * KNOWN BUG: returns true for a string of blanks
   */
  if (*cp == '\0')
    return 0;

  for ( ; *cp != '\0'; cp++) {
    if ( isalpha(*cp) || iscntrl(*cp) ) {
      return 0;
    }
  }
  return 1;
}


//--------------


// Input                rvec : Row vector (this)
// Output               void : The function sends questions to
//                             standard output interactively asking
//                             the user to supply the information to
//                             specify the vector.  The vector rvec is modified.
// Note                      : Function is for INTERACTIVE programs

void row_vector::ask()
  {
  int dim;
  double dr,di;
  std::cout << "\n\tPlease Input the Number of Elements: ";
  std::cin >> dim; 
  row_vector rvec(dim);
  for(int i=0; i<dim; i++ )
    {
    std::cout << "\n\tPlease Input Real and Imaginary Value of rvec("
         << i << ") [re im]: ";
    std::cin >> dr >> di;
    rvec.put(complex(dr,di),i);
    }
  *(this) = rvec;
  return;
  }

#endif						      // row_vector.cc
