/* col_vector.cc ************************************************-*-c++-*-
**									**
**		                   G A M M A				**
**                                                                      **
**	Column Vectors 		                 Implementation		**
**                                                                      **
**	Copyright (c) 1990, 2000					**
**	Tilo Levante, Scott A. Smith					**
**	Eidgenoessische Technische Hochschule				**
**	Labor fuer physikalische Chemie					**
**	8092 Zuerich / Switzerland					**
**                                                                      **
**      $Header: $
**                                                                      **
**                                                                      **
*************************************************************************/
/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
**  Class col_vector defines column vectors for C++ with the usual      **
**  algebraic operations, higher functions and I/O routines.            **
**                                                                      **
*************************************************************************/

#ifndef   Gcol_vector_cc_			// Is file already included?
#  define Gcol_vector_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler
#    pragma implementation			// Then this is the implementation
#  endif

#include <Matrix/col_vector.h>			// Include class interface
#include <Matrix/complex.h>			// Inlcude GAMMA complex numbers
#include <Matrix/MxModBas.h>			// Include Matrix module errors
#include <vector>				// Inlcude libstdc++ STL vectors
#include <cmath>				// Inlcude HUGE_VAL_VAL

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                    CLASS COLUMN VECTOR ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   cvec    : A column vector (this)
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

void col_vector::CVerror(int eidx, int noret) const
  {
  std::string hdr("Column Vector");
  std::string msg;
  switch (eidx)
    {
    case 29: MxModError(hdr,"Operation Is Disallowed",     noret); break;// (29)
    case 30: MxModError(hdr,"Multiple Columns Disallowed", noret); break;// (30)
    case 31: MxModError(hdr,"Cannot Assign From Row Vec",  noret); break;// (31)
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

void col_vector::CVerror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Column Vector");
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

volatile void col_vector::CVfatality(int eidx) const
  {
  CVerror(eidx, 1);                             // Normal non-fatal error
  if(eidx) CVerror(0);                          // Program aborting error
  MxModFatal();					// Keep screen nice, exit
  }

volatile void col_vector::CVfatality(int eidx, const std::string& pname) const
  {
  CVerror(eidx, pname, 1);                      // Normal non-fatal error
  if(eidx) CVerror(0);                          // Program aborting error
  MxModFatal();					// Keep screen nice, exit
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                  COLUMN VECTOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/*    Arguments                       Constructed Column Vector
    -------------       -------------------------------------------------------
          -             An empty col vector (size = 0)
          nc            A col vector with nc elements (columns) uninitialized
        nc, z           A col vector with nc elements (columns) all set to z
         cvec           A col vector equivalent to the input vector cvec
          mx            A col vector equivalent to input matrix mx
                        (in this case mx MUST be of dimension nx1)           */

col_vector::col_vector( )                       : matrix()              {}
col_vector::col_vector(int i)                   : matrix(i, 1)          {}
col_vector::col_vector(int i, double d)         : matrix(i, 1, d)       {}
col_vector::col_vector(int i, const complex& z) : matrix(i, 1, z)       {}
col_vector::col_vector(const col_vector& cvec)  : matrix((matrix&)cvec) {}
col_vector::col_vector(const matrix& mx)        : matrix(mx)
  {
  if(mx.cols() != 1)
    {
    std::string em = std::string(MxModdec(mx.rows()))
              + std::string("x")
              + std::string(MxModdec(mx.cols()));
    CVerror(20,em,1);                           // Construct multirow array
    CVfatality(40);                             // Can't do construction
    }
  }

col_vector::~col_vector () {}
 
// ____________________________________________________________________________
// B                   COLUMN VECTOR ACCESS AND ASSIGNMENT
// ____________________________________________________________________________

/* Function Arguments      Result        Function Arguments    Result
   -------- --------- -----------------  -------- --------- -------------------
     ( )        i     The ith element      get        i     Copy of ith element
    getRe       i     Real(ith element)   getIm       i     Imag(ith element)
    put        z,i    Set z as ith elem   elements   ---    Number of eleemnts
    size       ---    Number of elements

      -     cv1,cv2 Returns cv1-cv2

    Function/Operator                         Result
    -----------------   ------------------------------------------------------
          =             Current col vector set equal to input cvec
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

col_vector& col_vector::operator = (const col_vector& cvec)
  { 
  if(this != &cvec) matrix::operator= ( (matrix&)cvec );
  return *this;
  }

col_vector& col_vector::operator = (const row_vector& rvec)
  { CVerror(31,1); CVfatality(29); put(rvec.get(0),0); return *this; }

col_vector& col_vector::operator = (const matrix& mx)
  {
  if(mx.cols() != 1)
    {
    std::string em = std::string(MxModdec(mx.rows()))
              + std::string("x")
              + std::string(MxModdec(mx.cols()));
    CVerror(20,em,1);                   // Construct multicol array
    CVfatality(40);                     // Can't do construction
    }
  matrix::operator= ( (matrix&)mx );
  return *this;
  }

// *********************** col vector access functions ************************

complex& col_vector::operator() (int i) { return this->matrix::operator()(i,0); };
complex  col_vector::get(int i)   const { return matrix::get(i,0); }
double   col_vector::getRe(int i) const { return matrix::getRe(i,0); }
double   col_vector::getIm(int i) const { return matrix::getIm(i,0); }
void     col_vector::put(const complex& z, int i) { matrix::put(z,i,0); }
int      col_vector::elements( )  const { return rows(); }
int      col_vector::size( )      const { return rows(); }

// __________________________________________________________________________________
// F                  COLUMN VECTOR BINARY ARITHMETIC FUNCTIONS
// __________________________________________________________________________________

/* The operators which must be defined here, as opposed to those inherited from
   class matrix, are they which return col_vector (or at least, return a type
   other than what class matrix does.)D efinitions herein also allow for more
   explicit error messages.

   Additions/subtractions between column and row vectors are not defined here,
   consequently the error messages from doing so stems from class matrix.
   The only possible additions is when both are length 1 anyway.

   There are type casting problems when class matrix takes over some operators
   (e.g. mx*cvec, mx+cvec, mx-cvec). In that situation the result is proper but
   the return is matrix rather than column vector as one might like. This occurs
   because if matrix is the first object in the binary operation the matrix-matrix
   operator will be used over the matrix-column vector friend functions. Herein we
   must use friend functions because we cannot define the operators for matrix.

   Operator Arguments      Result        Operator Arguments       Result
   -------- --------- -----------------  -------- --------- -------------------
      +      cv1,cv2      cv1+cv2           *       mx,cv       mx*cv == cv1
      +      cv,mx        cv+mx             *       cv,z           z*cv
      -      cv1,cv2      cv1-cv2           *       cv,d           d*cv
      -      cv,mx        cv-mx             *       z,cv           z*cv
      *      cv1,cv2      cv1*cv2 (err)     *       d,cv           d*cv
      *      cv,rv        cv*rv == mx       /       cv,z         (1/z)*cv
      *      cv,mx        cv*mx == mx1      /       cv,d         (1/d)*cv         */

col_vector operator * (const complex& z, const col_vector& cvec)
  { return col_vector(z * (const matrix&)cvec); }

col_vector operator * (double d, const col_vector& cvec)
  { return col_vector(d * (const matrix&)cvec); }

col_vector col_vector::operator + (const col_vector& cvec) const
  {
  if(size() != cvec.size())
    {
    CVerror(42,1);					// Vec-vec size trouble
    CVfatality(21,"addition");				// Problems in addition
    }
  return col_vector((const matrix&)(*this) + (const matrix&)cvec);
  }

col_vector col_vector::operator + (const matrix& mx) const
  {
  if(size() != mx.rows() || mx.cols() != 1)
    {
    CVerror(41, 1);					// Mx-Vec dimensions bad
    std::string pname("Column Vector + Matrix");	// Error message string
    CVfatality(21,pname);				// Error in cvec + mx
    }
  return col_vector((const matrix&)(*this) + mx);
  }

col_vector col_vector::operator - (const col_vector& cvec) const
  {
  if(size() != cvec.size())
    {
    CVerror(42, 1);					// Vec-Vec dimensions bad
    std::string pname("Column Vector - Column Vector");	// Error message string
    CVfatality(21,pname);				// Error in cvec - cvec
    }
  return col_vector((const matrix&)(*this) - (const matrix&)cvec);
  }

col_vector col_vector::operator - (const matrix& mx) const
  {
  if(size() != mx.rows() || mx.cols() != 1)
    {
    CVerror(41, 1);					// Mx-Vec dimensions bad
    std::string pname("Column Vector - Matrix");	// Error message string
    CVfatality(21,pname);				// Error in cvec - mx
    }
  return col_vector((const matrix&)(*this) - mx);
  }

complex col_vector::operator * (const col_vector& cvec) const
  {
  CVerror(43, 1);					// Cant multiply 2 col vects
  std::string pname("Column Vector * Column Vector");	// Error message string
  CVfatality(21,pname);					// Error in cvec * cvec
  return complex0+cvec.get(0);				// This never returns
  }							// Needed over matrix to get
							// proper error messaging

matrix col_vector::operator * (const row_vector& rvec) const	//  mx  =  cvec  * rvec
  { return ((const matrix&)(*this))*((const matrix&)rvec); }	//(m x n) (m x 1) (1 x n)

matrix col_vector::operator * (const matrix& mx) const	//  mx1   =   cvec  *  mx
  { return ((const matrix&)(*this)) * mx; }		//(m x n)   (m x 1)  (1 x n)

// sosi - The next operator looses out to row_vector*col_vector in GCC 3.0
//        which is currently a friend function of class row_vector. It also is
//        dominated by mx*mx member function in some instances.

col_vector operator * (const matrix& mx, const col_vector& cvec)
  {
  if(cvec.size() != mx.cols())
    {
    cvec.CVerror(41, 1);				// Mx-Vec dimensions bad
    std::string pname("Mx * Column Vector");
    cvec.CVfatality(21,pname);
    }							//  cvec1  =   mx   * cvec
  return col_vector((mx) * (const matrix&)cvec);	// (m x 1)   (m x n) (n x 1)
  }

col_vector col_vector::operator * (const complex& z) const
  { return col_vector(z * (const matrix&)(*this)); }
col_vector col_vector::operator * (double d) const
  { return col_vector(d * (const matrix&)(*this)); }

col_vector col_vector::operator / (const complex& z) const
  { return col_vector((const matrix&)(*this)/z); }
col_vector col_vector::operator / (double d) const
  { return col_vector((const matrix&)(*this)/d); }

// ____________________________________________________________________________
// G                  COLUMN VECTOR UNARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________
 
/* Note: These operators, as they are unable to change the type of the
         result, can only be those which return col_vector.

   Operator Arguments      Result        Operator Arguments    Result
   -------- --------- -----------------  -------- --------- -------------------
      +=    cv,cv1    cv1 added to cv       +=    cv,mx     mx added to cv
      *=    cv,z      cv1 mult. by z        /=    cv,z      cv divided by z  */
 
col_vector& col_vector::operator += (const col_vector& cvec)
{
  if(!cvec.size()) 
    return *this;				// Nothing if cvec NULL
  if(!size())						// If we are NULL, become cvec
  {
    matrix::operator= ( (const matrix&)cvec );
    return *this;
  }
  if(size() != cvec.size())
  {
    cvec.CVerror(42, 1);				// Vec-Vec dimensions bad
    std::string pname("Column Vector += Column Vector");// Error message string
    cvec.CVfatality(21,pname);				// Error in cvec += cvec
  }
  
  matrix::operator+= ((const matrix&)cvec);
  return *this;
}

void col_vector::operator += (const matrix& mx)
  {
  if(!mx.rows()) return;				// Nothing if mx NULL
  if(size() != mx.rows() || mx.cols() != 1)
    {
    CVerror(41, 1);					// Vec-Mx dimensions bad
    std::string pname("Column Vector += Matrix");	// Error message string
    CVfatality(21,pname);				// Error in cvec += mx
    }
  matrix::operator+= (mx);
  }

col_vector& col_vector::operator -= (const col_vector& cvec)
{
  if(!cvec.rows()) 
    return *this;				// Nothing if mx NULL
  if(!size())						// Set -cvec if we are empty
  {
    (*this) = cvec.matrix::operator-();
    return *this;
  }
  if(size() != cvec.size())
  {
    CVerror(42, 1);					// Vec-Vec dimensions bad
    std::string pname("Column Vector -= Column Vector");// Error message string
    CVfatality(21,pname);				// Error in cvec -= cvec
  }
  
  matrix::operator-= ((const matrix&)cvec);
  return *this;
}

void col_vector::operator -= (const matrix& mx)
  {
  if(!mx.rows()) return;				// Nothing if mx NULL
  if(size() != mx.rows() || mx.cols() != 1)
    {
    CVerror(41, 1);					// Vec-Mx dimensions bad
    std::string pname("Column Vector -= Matrix");	// Error message string
    CVfatality(21,pname);				// Error in cvec -= mx
    }
  matrix::operator-= (mx);
  }

col_vector& col_vector::operator *= (const complex& z) 
{ 
  matrix::operator*= (z);
  return *this;  
}
col_vector& col_vector::operator *= (      double d)   
{ 
  matrix::operator*= (d); 
  return *this;
}
col_vector & col_vector::operator /= (const complex& z) 
{ 
  matrix::operator/= (z); 
  return *this;
}

col_vector & col_vector::operator /= (const double d)   
{ 
  matrix::operator/= (d); 
  return *this;
}

// ____________________________________________________________________________
// J                  COLUMN VECTOR SIMPLE UNARY FUNCTIONS
// ____________________________________________________________________________

// Note: The negation operator must be defined herein.  The function conjugate
//       is inherited, but adjoint and transpose return row_vectors so they
//       are also here.  Transpose is included just for specific error message.

        // Input		cvec : A column vector
        // Output adjoint	rvec : A row vector of the adjoint of cvec
        //                                               t        * T
        //                                <rvec| = |cvec>  = |cvec >
        // Output transpose	rvec : A row vector of the transpose of cvec
        //                                                T
        //                                 <rvec| = |cvec>
        // Output trace	        void : Trace is disallowed on vectors
 

// col_vector col_vector::operator- ();                   // INHERITED

row_vector adjoint(const col_vector& cvec)
  { return row_vector(adjoint((const matrix&)cvec)); }

row_vector transpose(const col_vector& cvec)
  { return row_vector(transpose((const matrix&)cvec)); }
 
complex trace(const col_vector& cvec)
  {
  cvec.CVfatality(44);			// Vec-Mx dimensions bad
  return cvec.get(0);			// Compiler likes cvec to be used
  }

row_vector col_vector::adjoint()   const
  { return row_vector(matrix::adjoint()); }
row_vector col_vector::transpose() const
  { return row_vector(matrix::transpose()); }
complex    col_vector::trace()      const
  {
  CVfatality(44);			// Vec-Mx dimensions bad
  return get(0);
  }

col_vector col_vector::differential() const
  {
  int nr = size();
  col_vector cvec(nr);
  for(int i=0; i<nr-1; i++)
    cvec.put((*this).get(i+1) - (*this).get(i), i);
  cvec.put((cvec)(nr-2), nr-1);
  return cvec;
  }
 
// ____________________________________________________________________________
// I                    COLUMN VECTOR SIMPLE BINARY FUNCTIONS
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
                     |  \            2 |            2                 *
          Norm = sqrt|  /  |<i|cvec>|  |, |<i|cvec>| = <i|cvec>*<i|cvec>
                     |  ---            |
                     [  i=0            ]

  & the following definitions apply to the max, min, and sort functions:
       
                          type         max/min
                            0       norm (default)
                           <0         imaginary                         
                           >0           real                                 */

double col_vector::norm() const
  {
  int n=size();
  double d=0;
  for (int i=0; i<n; i++)
    d += square_norm(get(i));
  return sqrt(d);
  }

complex col_vector::sum() const
  {
  complex z(0);
  int n = size();
  for(int i=0; i<n; i++) z += get(i);
  return z;
  }

complex col_vector::sum(int st, int ne) const
  {
  complex z(0);
  for(int i=st; i<st+ne; i++)
    z += get(i);
  return z;
  }

double  col_vector::maxRe() const { return getRe(max(1)); }
double  col_vector::maxIm() const { return getRe(max(-1)); }
complex col_vector::maxZ()  const { return get(max(0)); }
double  col_vector::minRe() const { return getRe(min(1)); }
double  col_vector::minIm() const { return getRe(min(-1)); }
complex col_vector::minZ()  const { return get(min(0)); }

int col_vector::max(int type) const
  {
  int n = size();			// Number of vector elements
  int imax=0;
  double max=-HUGE_VAL, tmp=0;
  if(type>1) type=1;
  else if(type<0) type=2;
  int i;
  switch(type)
    {
    case 0:
      for(i=0; i<n; i++)		// Sum them all up
        {
        tmp = square_norm(get(i));
        if(tmp > max) { max = tmp; imax = i; }
        }
      break;
    case 1:
      for(i=0; i<n; i++)		// Sum them all up
        {
        tmp = getRe(i);
        if(tmp > max) { max = tmp; imax = i; }
        }
      break;
    case 2:
      for(i=0; i<n; i++)		// Sum them all up
        {
        tmp = getIm(i);
        if(tmp > max) { max = tmp; imax = i; }
        }
      break;
    }
  return imax;
  }

int col_vector::min(int type) const
  {
  int n = size();			// Number of vector elements
  int imin=0;
  double min=HUGE_VAL, tmp=0;
  if(type>1) type=1;
  else if(type<0) type=2;
  int i;
  switch(type)
    {
    case 0:
      for(i=0; i<n; i++)		// Sum them all up
        {
        tmp = square_norm(get(i));
        if(tmp < min) { min = tmp; imin = i; }
        }
      break;
    case 1:
      for(i=0; i<n; i++)		// Sum them all up
        {
        tmp = getRe(i);
        if(tmp < min) { min = tmp; imin = i; }
        }
      break;
    case 2:
      for(i=0; i<n; i++)		// Sum them all up
        {
        tmp = getIm(i);
        if(tmp < min) { min = tmp; imin = i; }
        }
      break;
    }
  return imin;
  }

void col_vector::flip()
  {
  int n = size();			// Number of vector elements
  complex ztmp(0);                    // Temporary value
  int i, j;
  for(i=0, j=n-1; i<n/2; i++,j--)	// Loop over the 1st half
    {
    ztmp += get(i);			// Store ith element
    put(get(j), i);			// Set ith elem to (n-1-i)th element
    put(ztmp, j);			// Set (n-1-i)th elem to ith elem
    }
  }

void col_vector::zero()
  { for(int i=0; i<size(); i++) put(complex0, i); }

std::vector<int> col_vector::sort(int type) const
  {
  int      nvals = size();		// Number of values to be sorted
  double *vals;				// Array of values to be sorted
  vals = new double[nvals];		// Allocate the array of values
  std::vector<int> indx;		// Allocate vector to return
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
// J               COLUMN VECTOR COMPLEX UNARY FUNCTIONS
// ____________________________________________________________________________

// ******************************      FFT      *******************************

col_vector FFT(const  col_vector& cvec) { return FFT((const  matrix&)cvec); }
col_vector IFFT(const col_vector& cvec) { return IFFT((const matrix&)cvec); }
 
        // Input        cvec : A column vector
        // Output      cvec1 : A column vector containing the discrete Fourier
        //                     (inverse) transform of the values in cvec
        // Note              : The vector cvec should be of dimension
        //                     1 x (2^n)
 
// ____________________________________________________________________________
// K                     COLUMN VECTOR CONTAINER SUPPPORT
// ____________________________________________________________________________
 
        // Input      cvec1  : A column vector
        //            cvec2  : A column vector
        // Output       T/F  : Returns TRUE if cvec1 is equal to cvec2.

//bool col_vector::operator==(const col_vector& cvec) const	// INHERITED
//bool col_vector::operator!=(const col_vector& cvec) const	// INHERITED
//bool col_vector::operator> (const col_vector& cvec) const	// INHERITED
//bool col_vector::operator< (const col_vector& cvec) const	// INHERITED

// ____________________________________________________________________________
// K                     COLUMN VECTOR VECTOR PRODUCT FUNCTIONS
// ____________________________________________________________________________
 
//*************************** straight product ********************************
 
        // Input              cvec   : A column vector
        // Output             cvec1  : A "product" of cvec with itself
        // Note                      : This is NOT a scalar product
        
        //                               cvec1(i) = cvec(i) * cvec(i)

col_vector col_vector::product() const
  {
  int dim = size();
  col_vector cvec1(dim);
  for(int i=dim-1; i>=0; i--)
    cvec1.put(get(i)*get(i), i);
  return cvec1;
  }
 
        // Input              cvec   : A column vector
        //                    rvec   : A row vector
        // Output             cvec1  : A "product" of cvec and rvec
        // Note                      : This is NOT a scalar product
        //                               cvec1(i) = cvec(i) * rvec(i)
        //                               cvec3(i) = cvec1(i) * cvec2(i)

col_vector product(const col_vector& cvec, const row_vector& rvec)
  { return cvec.product(rvec); }

col_vector product(const col_vector& cvec1, const col_vector& cvec2)
  { return cvec1.product(cvec2); }

col_vector col_vector::product(const row_vector& rvec) const
  {
  int dim = size();
  if(dim != rvec.size())
    {
    CVerror(42,1);				// Vec-vec size trouble
    CVerror(21,"Col Vector (x) Row Vec");	// Problems in this
    std::string pname("product");
    CVfatality(5, pname);
    }
  col_vector cvec1(dim);
  for(int i=dim-1; i>=0; i--)
    cvec1.put(get(i)*rvec.get(i), i);
  return cvec1;
  }

col_vector col_vector::product(const col_vector& cvec2) const
  {
  int dim = size();
  if(dim != cvec2.size())
    {
    CVerror(42,1);				// Vec-vec size trouble
    CVerror(21,"Col Vector (x) Col Vec");	// Problems in this
    std::string pname("product");
    CVfatality(5, pname);
    }
  col_vector cvec3(dim);
  for(int i=dim-1; i>=0; i--)
    cvec3.put(get(i)*cvec2.get(i), i);
  return cvec3;
  }

 
//************************* scalar product ******************************
 
        // Input              cvec  : A Column vector
        // Output             z      : Scalar product of adjoint(cvec) & cvec

        //                             z = <cvec'|cvec>,   <cvec'| = adjoint(|cvec>)

double scalar_product(const col_vector& cvec)
  { return cvec.scalar_product(); }

double col_vector::scalar_product() const
  {
  int dim = size();
  double sp=0;
  for(int i=dim-1; i>=0; i--)
    sp += square_norm(get(i));
  return sp;
  }

        // Input              cvec1  : A column vector
        //                    cvec2  : Second column vector
        // Output             z      : Scalar product of adjoint(cvec1) & cvec2
 
        //                             z = <cvec1'|cvec2>, <cvec1'> = adjoint(|cvec>)
        // Input              cvec   : A column vector
        //                    rvec   : A row vector
        // Output             z      : Scalar product of rvec & cvec
 
        //                             z = <rvec|cvec>


complex scalar_product(const col_vector& cvec1, const col_vector& cvec2)
  { return cvec1.scalar_product(cvec2); }

complex scalar_product(const col_vector& cvec, const row_vector& rvec)
  { return cvec.scalar_product(rvec); }

complex col_vector::scalar_product(const col_vector& cvec2) const
  {
  int dim = size();
  if(dim != cvec2.size())
    {
    CVerror(42, 1);				// Vec-Vec dimensions bad
    CVerror(21,"Col Vector (x) Col Vec");	// Problems in this
    std::string pname("scalar_product");
    CVfatality(5, pname);
    }
  complex z(0);
  for(int i=dim-1; i>=0; i--)
    z += (cvec2.get(i)).conj_times(get(i));
  return z;
  }

complex col_vector::scalar_product(const row_vector& rvec) const
  {
  int dim = size();
  if(dim != rvec.size())
    {
    CVerror(42, 1);				// Vec-Vec dimensions bad
    CVerror(21,"Col Vector (x) Row Vec");	// Problems in this
    std::string pname("scalar_product");
    CVfatality(5, pname);
    }
  complex z(0);
  for(int i=dim-1; i>=0; i--)
    z += get(i)*rvec.get(i);
  return z;
  }


// ____________________________________________________________________________
// L                      COLUMN VECTOR I/O FUNCTIONS
// ____________________________________________________________________________

// ------------------------ ASCII Output Functions ----------------------------

/*              Input           cvec    : A column vector (this)
                                ostr    : Output ASCII file stream
                                full    : Flag for amount of output
                Return          void    : cvec is sent to the output stream  */


std::string col_vector::hdrString() const
  {
  std::string hdr;
  if(!size()) hdr = "Empty Column Vector";
  else
    {
    std::string vtype("");
    if(is_real())           vtype = " Real";
    else if(is_imaginary()) vtype = " Imaginary";
    hdr = MxModdec(size()) + std::string(" x 1")
        + vtype + std::string(" Column Vector");
    }
  return hdr;
  }
/* The column vector is placed into an array of strings, each string having
   equal length. If the number of vector elements exceeds the number of rows
   we allow for simplistic printing, then we begin an new column of vector
   elements. The number of rows before additonal columns are begun is set
   by PFlgs.VxRows. The max # of columns per line is set by PFlgs.VxCols.    */

std::vector<std::string> col_vector::printStrings(const MxPrint& PFlgs) const
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

//                    Print If Single Element Per Row

  int nrow = PFlgs.VxRows;			// Rows in output
  if(npts <= nrow)				// If less points than rows
    {						// no need for multiple columns
    for(int pos=0; pos<npts; pos++ )
      {
      if(ptype==1)					// Output as 
        vline = MxModform(efmt.c_str(),getRe(pos)); 	// a real element
      else if(ptype==2)					// Output as
        vline = MxModform(efmt.c_str(),getIm(pos)); 	// an imaginary element
      else  vline = get(pos).printString();         	// Output as complex   
      PStrings.push_back(vline);
      }
    return PStrings;
    }

//               Print If Mulitple Elements Per Row

  int ncol = PFlgs.VxCols;			// Columns in output
  if(ptype) ncol *=2;				// Twice as many if real/imag
  int row=0, col=0;				// Column index
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

std::ostream& col_vector::printcols(std::ostream& ostr, int ncol, int rc, int nelem) const
  {
  if(!size()) return ostr;
  std::string marg(10, ' ');
  std::string spc(3, ' ');
  int col=0;
  if(nelem < 1) nelem = size();
  for(int i=0; i<nelem; i++ )
    {
    if(col == 0)  ostr << marg;
    if(!rc)       ostr << getRe(i);
    else if(rc<1) ostr << getIm(i);
    else          ostr << get(i);
    col++;
    if(col >= ncol) { ostr << "\n"; col=0; }
    else              ostr << spc;
    }
  return ostr;
  }

std::ostream& col_vector::print(std::ostream& ostr, int full) const
  {
  if(!size()) { ostr << "Empty Column Vector\n"; return ostr; }
  ostr << "GAMMA " << "1 x " << size() << " Column Vector";
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

std::ostream& operator << (std::ostream& ostr, const col_vector& cvec)
  { return cvec.print(ostr, 1); }

std::istream& operator >> (std::istream& istr, col_vector& cvec)
  {
  matrix t1;
  istr >> t1;
  cvec = col_vector(t1);
  return istr;
  }

        // Input                cvec : Column vector (this)
        // Output               void : The function sends questions to
        //                             standard output interactively asking
        //                             the user to supply the information to
        //                             specify the vector. Vector is modified.
        // Note                      : Function is for INTERACTIVE programs
 
void col_vector::ask()
  {
  int dim;
  double dr,di;
  std::cout << "\n\tPlease Input the Number of Elements: ";
  std::cin >> dim; 
  col_vector cvec(dim);
  for(int i=0; i<dim; i++ )
    {
    std::cout << "\n\tPlease Input Real and Imaginary Value of cvec("
         << i << ") [re im]: ";
    std::cin >> dr >> di;
    cvec.put(complex(dr,di),i);
    }
  *(this) = cvec;
  return;
  }

#endif						// col_vector.cc
