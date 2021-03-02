/* h_matrix.cc **************************************************-*-c++-*-
**								   	**
**                               G A M M A 				**
**								   	**
**	Hermitian Matrix                        Implementation		**
**								   	**
**	Copyright (c) 1993						**
**	Scott Smith							**
**	University of California, Santa Barbara				**
**	Department of Chemistry						**
**	Santa Barbara, CA, 93106, USA					**
**							   		**
**      $Header: $
**								   	**
*************************************************************************/

/*************************************************************************
**							   		**
** Description						   		**
**									**
** The class h_matrix defines Hermitian arrays for C++ with in the	**
** GAMMA matrix class hierarchy.  Storage of Hermitian matrices is	**
** limited to the upper-triangular elements as are computations 	**
** involving them whenever possible.  Defined herein are also the usual	**
** algebraic operations {+,-,*,/}, and I/O routines, and complex	**
** functions {adjoint, trace,...}					**
**							   		**
*************************************************************************/

#ifndef   Gh_matrix_cc_				// Is file already included?
#  define Gh_matrix_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#  endif

#include <GamGen.h>				// Include OS specific stuff
#include <Matrix/h_matrix.h>			// Include Hermitian array interface
#include <Matrix/_matrix.h>     	        // Include base matrix type	
#include <Matrix/n_matrix.h>			// Include normal array interface
#include <Matrix/d_matrix.h>			// Include diagonal array interface
#include <Matrix/i_matrix.h>			// Include identity array interface
#include <Matrix/MxModBas.h>			// Include Matrix module errors 
#include <fstream>				// Include libstdc++ filestreams
#include <string>
#include <cmath>				// Inlcude HUGE_VAL_VAL
#include <stdlib.h>

#ifdef _USING_BLAS_
  #if defined(__APPLE__)
    #include <Accelerate/Accelerate.h>
  #else  
    #ifdef _USING_GOTOBLAS2_
    extern "C"{
    #endif
    #include <cblas.h>
    #ifdef _USING_GOTOBLAS2_
    }
    #endif
  #endif
#endif

#ifdef _USING_SUNPERFLIB_
 #include "gsunperf.h"
#endif 

#ifdef _USING_LAPACK_
 #if defined(__APPLE__)
  #include <Accelerate/Accelerate.h>
 #else
  #include <clapack.h>
 #endif  
#endif

#ifdef _USING_GOTOBLAS2_
 #define _USING_LAPACK_
 #define __CLPK_doublecomplex double
extern "C"
{ extern void zgeev_(const char*, const char*, const int*, __CLPK_doublecomplex *, const int*,
                     __CLPK_doublecomplex*, __CLPK_doublecomplex*, const int*, __CLPK_doublecomplex*,
                     const int*, __CLPK_doublecomplex*, const int*, double*, int*); 
  extern void zheev_(const char*, const char*, const int*, __CLPK_doublecomplex*, const int*,
                     double*, __CLPK_doublecomplex*, const int*, double*, int*);
}
#endif

// ____________________________________________________________________________
// i                    CLASS H_MATRIX ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   mx      : A Hermitian matrix (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : string in message
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void h_matrix::HMxerror(int eidx, int nr) const
  {
  std::string CL="Hermitian Matrix";
  Mxerror(CL, eidx, nr);
  }
 
void h_matrix::HMxerror(int eidx, const std::string& PN, int nr) const
  {
  std::string CL="Hermitian Matrix";
  Mxerror(CL, eidx, PN, nr);
  }
 
volatile void h_matrix::HMxfatal(int eidx) const
  {
  std::string CL="Hermitian Matrix";
  Mxfatality(CL, eidx);
  }
 
volatile void h_matrix::HMxfatal(int eidx, const std::string& PN) const
  {
  std::string CL="Hermitian Matrix";
  Mxfatality(CL, eidx, PN);
  }

// ____________________________________________________________________________
// A               CLASS H_MATRIX CONSTRUCTORS AND DESTRUCTOR
// ____________________________________________________________________________

/*  Arguments                      Constructed Array
    ---------		-------------------------------------------------------
        -	 	 An empty array
     nr, nc		 An (nr x nc) array, uninitialized
    nr, nc, z		 An (nr x nc) array, sets <i|mx|i>=Re(z) & <i<j|mx|j>=z
       hmx		 A duplicate of Hermitian array hmx

   Note that all constructors invariably call the base class constructor
   (i.e. the one in _matrix) which sets the row & column sizes.  The added
   code should simply allocate the data array of proper size. Similarly, the
   destructor only destroys the data array, anything else is destroyed in 
   the base class.                                                           */

h_matrix::h_matrix( )              : _matrix()      { data=new complex[size]; }	
h_matrix::h_matrix(int nr, int nc) : _matrix(nr,nc)
  {
  if(nr != nc)
    {
    HMxerror(9, 1);				// Problems in constructor
    HMxfatal(2);				// Rectangular array
    }
  size = (nr*nr+nr)/2;				// Size for upper triangle
  data = new complex[size];			// Allocate the array space
  }

h_matrix::h_matrix(int nr, int nc, const complex& z) : _matrix(nr,nc)
  {
  if(nr != nc)
    {
    HMxerror(9, 1);				// Problems in constructor
    HMxfatal(2);				// Rectangular array
    }
  size = (nr*nr+nr)/2;				// Size for upper triangle
  data = new complex[size];			// Allocate data array
  for(int pos=0; pos < size; pos ++)		// Set all array eleemnts
    data[pos] = z;
  for(int k=0; k<nr; k++)			// Now insure diagonal ones
    data[k*nr-(k*(k-1))/2] = Re(z);		// are real......
  }

h_matrix::h_matrix(const h_matrix& mx) : _matrix(mx)
  {
  size = mx.size;
  data = new complex[size];
  for(int pos=0; pos<size; pos ++)
    data[pos] = mx.data[pos];
  }

h_matrix::~h_matrix () { delete [] data; }		// Delete matrix data


// ____________________________________________________________________________
// B                  CLASS H_MATRIX ACCESS AND ASSIGNMENT
// ____________________________________________________________________________

/* There are two flavors of element access operators herein. The one that might
   commonly be used, i.e. mx(i,j), is DIRECTLY accessing the element <i|mx|j>.
   That can save CPU time by avoiding value copying when altering or obtaining
   elements.  However its misuse can lead to trouble in this class because not
   all Hermitian matrix elements are stored and there are restrictions on the
   elements so that the array remains Hermitian.  Two examples of when this is
   trouble would be "hmx(0,0)=comlexi;" and "hmx(2,0)=7".  The former fails as
   diagonal elements MUST be real and the latter fails because lower triangle
   elements are not stored so they don't exist.

   In constrast, the "get" function returns copies of the element & the "put"
   function checks that the element is being properly set.  Although these
   are slower they are absolutely safe to use.  The put function returns TF
   and will be FALSE if an element is attempted to be set non-hermitian. Note
   that a FALSE put will trigger class matrix (that which external users use)
   to change the matrix type and then procede with the put.
 
   The (i,j) operator is fast but dangerous if one is not careful or doesn't
   know what they are doing.  So, the (i,j) operator is used inside this class
   and other _matrix derived classes for speed, but this is protected from
   use outside - i.e. protected from the typical user. The typical user deals
   with the generic matrix class (class matrix) whose (i,j) operator works
   in "safe" mode.
 
    Function/Operator                         Result
    -----------------   ------------------------------------------------------
          =             Current (Hermitian) array set equal to input hmx
    (int,int)           Reference to element <i|hmx|j> (Potential Danger)
    get(int,int)        Copy of element <i|hmx|j>      (Safe)
    put(int,int)        Assigns element <i|hmx|j>      (Return FALSE if fails)
    put_h(int,int)      Assigns both <i|hmx|j> & <j|hmx|i> (Return TRUE/FALSE)
    get_block(r,c,R,C)  Returns hmx block of size RxC starting at <r|hmx|c>
    put_block(r,c,mx)   Places mx into hmx at position <r|hmx|c> (TRUE/FALSE)
 
    Again, note that operations which return TF will trigger class matrix to
    take appropriate action when they return FALSE.  This usually means that
    the operation cannot be handled by a Hermitian array so the array type
    is changed to be more generic.  The the operation is reattempted.       */
 
void h_matrix::operator= (const h_matrix& mx)
  {
  if(this == &mx) return; 		// Do nothing on self-assign
  delete [] data;			// First delete the data
  _matrix::operator= (mx);		// Assign by generic function
  data = new complex[size];		// Construct new data array
  for(int pos=0; pos<size; pos ++)	// Copy the data array
    data[pos] = mx.data[pos];
  }

complex& h_matrix::operator() (int i, int j)
  { 
  if(i<=j) return data[i*cols_-(i*(i-1))/2+j-i];
  cval = conj(data[j*cols_-(j*(j-1))/2+i-j]); 
  return cval;
  }

complex h_matrix::get(int i, int j) const
  { return (i<=j)?data[i*cols_-(i*(i-1))/2+j-i]
            :conj(data[j*cols_-(j*(j-1))/2+i-j]); }

bool h_matrix::put(const complex& z, int i, int j)
  {
  if(i==j)					// If setting a diagonal
    {						// it must be real or it
    if (fabs(Im(z)) > 1e-12) 			// If not real or off-diagonal
      {
//
//JABE, 9-9-2008
//
//Previously, this return of a false was used. This can lead to many problems though
//and now an error message is given and the program will stop
//
//    return false;			// cannot stay Hermitian
//
      std::cerr << "\n\nGAMMA ERROR: diagonal matrix element is too big for Hermitian matrix (" << Re(z) <<","<< Im(z) << "). \n Discarding imaginary part ...\n\n";
//    exit(99);
      }

    data[i*cols_-(i*(i-1))/2] = Re(z);	// Set the element

    return true;
    }
  else  				// Return False for non-diagonal elements
    return false;
  }

bool h_matrix::put_h(const complex& z, int i, int j) 
  { 
  if(i==j)					// If setting a diagonal
    { 
    if (fabs(Im(z)) > 1e-12) 			// If not real or off-diagonal
      {
//
//JABE, 9-9-2008
//
//Previously, this return of a false was used. This can lead to many problems though
//and now an error message is given and the program will stop
//
//    return false;			// cannot stay Hermitian
//
      std::cerr << "\n\nGAMMA ERROR: diagonal matrix element is too big for Hermitian matrix (" << Re(z) <<","<< Im(z) << "). \n Discarding imaginary part ...\n\n";
//    exit(99);
      }

    data[i*cols_-(i*(i-1))/2] = Re(z);
    }
  else if(i<j)					// This for upper triangle
    data[i*cols_-(i*(i-1))/2+j-i] = z;
  else if(i>j)					// This for lower triangle
    data[j*cols_-(j*(j-1))/2+i-j] = conj(z);
  return true;
  }


_matrix* h_matrix::get_block(int row, int col, int nrows, int ncols)
  {
  int r = rows_;				// Get the dimension of hmx
  if(!row && !col)                              // If block starts at <0|hmx|0> and
    if((nrows==r) && (ncols==r))		// spans all of hmx, just return hmx
      return this;
  int i,j,k,cnew;
  int pos, posnew;				// Old, new element positions
  if(row==col && nrows==ncols)			// Reqest is square about diagonal
    {						// so return will be Hermitian!
    h_matrix* hmx = new h_matrix(nrows,ncols);
    posnew = 0;					// Start with first hmx element
    for(i=0,k=row; i<nrows; i++,k++)		// Loop over all the rows of hmx
      {
      pos =  k*r - (k*(k-1))/2;			// Old position at <row|mx|row>
      cnew = ncols-i;				// Want this many elements this row 
      for(j=0; j<cnew; j++,pos++,posnew++)	// Loop over each of cnew columns
        hmx->data[posnew] = data[pos];		// & fill up the new matrix
      }
    return hmx;
    }
  else						// Request is non-square or else 
    {						// off-diagonal, so then we must
    n_matrix* nmx = new n_matrix(nrows, ncols); // return a n_matrix
    int IJ, ij, ji, nextr;			// Access to nmx: <I|nmx|J> = nmx->data[IJ]
    int iend = row+nrows;			// Last column of hmx in block
    int jend = col+ncols;			// Last column of hmx in block
    for(IJ=0, i=row; i<iend; i++)		// Loop over both nmx & hmx rows in block
      {
      j = col;					// Set initial j index for this row
      if(col < i)				// See if row starts in lower triangle
        {
        nextr = r-i;				// Row jump amount for hmx on row i
        ji = col*r-(col*(col-1))/2+i-col;	// <j|hmx|i> = data[ji] -> <col|hmx|i>
        for( ; j<i; j++, IJ++, ji+=nextr)	// Loop over ith hmx row, j<i (low. tri.)
          ((n_matrix*)nmx)->data[IJ] 		// <I|nmx|J> = <i|hmx|j> = <j|hmx*|i>
                              = conj(data[ji]);	// Note: ji = ii at loop exit
        for(ij=ji; j<jend; j++, IJ++, ij++)	// Loop over ith hmx row, j>=i (up. tri.)
          ((n_matrix*)nmx)->data[IJ]=data[ij];	// <I|nmx|J> = <i|hmx|j>
        }
      else					// Here if row starts in upper triangle
        {
        ij = i*r - (i*(i-1))/2 + col - i;	// <i|hmx|j> = data[ij] -> <i|hmx|col>
        for( ; j<jend; j++,IJ++,ij++)		// Loop over ith hmx row, j>=i (up. tri.)
          ((n_matrix*)nmx)->data[IJ]=data[ij];	// <I|nmx|J> = <i|hmx|j>
        }
      }
    return nmx;
    }
  }
		      
bool h_matrix::put_block(int row, int col, _matrix* mx)
  {
  int ncols = mx->cols();			// Columns in input block
  int nrows = mx->rows();			// Rows in input block
  int c = cols();				// Columns in hmx
  int r = rows();				// Rows in hmx
  if((row+nrows > r) || (col+ncols > c))	// Insure mx dim. not exceeded
    {
    HMxerror(5, "put_block", 1);		// Bad use of function 
    HMxfatal(52);				// Array dimensions exceeded
    }
  if(row!=col || nrows!=ncols)			// Block MUST be square about 
    return false;				// diag. else hmx -> nmx convert
  int i, j, k;
  int pos, pmx, cblk;
  switch(mx->stored_type())
    {
    case h_matrix_type:				// Put Hermitian matrix into Hermitian
      for(i=0, k=row; i<nrows; i++,k++)		// Loop over all the rows of mx
        {
        pos = k*c-(k*(k-1))/2;	 		// Position of <k|hmx|k>
        pmx = i*ncols-(i*(i-1))/2; 		// Position of <i|mx|i>
        cblk = ncols-i;				// Columns of mx to copy from row i
        for(j=0;j<cblk;j++,pos++,pmx++)		// Loop over each of cnew columns of mx
          data[pos] = 				// & fill up the corresponding hmx elements
                   ((h_matrix*)mx)->data[pmx];
        }
      break;
    case i_matrix_type:				// Put identity matrix into Hermitian
      for(i=0, k=row; i<nrows; i++,k++)		// Loop over all the rows of mx
        {
        pos = k*c-(k*(k-1))/2;	 		// Position of <k|hmx|k>
        cblk = ncols-i;				// Columns of mx to copy from row i
        data[pos] = complex1;	 		// Set <k|hmx|k> = 1;
        pos++;					// Skip to <k|hmx|k+1>
        cblk--;					// Decrement number of columns to zero
        for(j=0;j<cblk;j++,pos++)		// Loop over each of cnew columns of mx
          data[pos] = complex0;	 		// (upper triangle) and zero new matrix
        }
      break;
    case n_matrix_type:				// Put normal matrix into Hermitian
    default:					// Put unknown matrix into Hermitian
      return false;					// demands hmx -> nmx type conversion
      break;
    case d_matrix_type:				// Put diagonal matrix into Hermitian
      if(!mx->is_real())
        return false;
      for(i=0, k=row; i<nrows; i++,k++)		// Loop over all the rows of mx
        {
        pos = k*c-(k*(k-1))/2;	 		// Position of <k|hmx|k>
        cblk = ncols-i;				// Columns of mx to copy from row i
        data[pos] = ((d_matrix*)mx)->data[i]; 	// Set <k|hmx|k> = <i|dmx|i>
        pos++;					// Skip to <k|hmx|k+1>
        cblk--;					// Decrement number of columns
        for(j=0;j<cblk;j++,pos++)		// Loop over each of cnew columns of mx
          data[pos] = complex0;	 		// (upper triangle) and zero new matrix
        }
      break;
    }
  return true;
  }

// ____________________________________________________________________________
// C                 CLASS H_MATRIX HERMITIAN_TYPE HANDLING
// ____________________________________________________________________________

/* These functions handle the checking of whether an array is of a Hermitian
   type.  Since this class embodies a Hermitian matrix, the tests are easy and
   dont actually used the input value d.  These are somewhat outdated now that
   GAMMA has a Hermitian array class and will probably be removed in favor of
   a TF is_hermitian function function or something.  More appropriate would be
   an stored_unitary function.....
 
     Function            Output                       Description
   ----------------  --------------  ------------------------------------------
   stored_hermitian  hermitian_type  Returns _hermitian=TRUE
   test_hermitian    hermitian_type  Returns _hermitian=TRUE (d is unused)   */

hermitian_type h_matrix::stored_hermitian( )      const { return _hermitian; }
hermitian_type h_matrix::test_hermitian(double d) const
                                                 { return _hermitian; d=0.0; }

// ____________________________________________________________________________
// D                 CLASS H_MATRIX MATRIX_TYPE HANDLING
// ____________________________________________________________________________

/* These functions handle the checking of what type of array we have.  These
   types directly correspond to the matrix classes (such as h_matrix) derived
   from class _matrix.

     Function            Output                       Description
   ----------------  --------------  ------------------------------------------
   set_type            --------      Found in class matrix
   test_type           --------      Found in class matrix 
   stored_type       h_matrix_type   Alway returns we are Hermitian
   test_hermitian    *_matrix_type   Type hmx could be within d
   mxtype            string          Returns the string "Hermitian" 
 
   The test type looks to see if hmx could be an array of the matrix type
   specified as an input argument.  If it can within "d" the return type 
   is equal to the input type.  If it cannot the return is h_matrix_type     */

matrix_type h_matrix::stored_type( ) const { return h_matrix_type; }
matrix_type h_matrix::test_type(matrix_type m, double d) const
  {
  int f;				// T/F flag if conversion possible
  int pos, cblk;
  int i,j;				// Column and row counters
  int ncols = cols_;
  switch(m)
    {
    case h_matrix_type:			// h_matrix can be stored h_matrix
      return h_matrix_type;
    case n_matrix_type:			// h_matrix can always be stored n_matrix
      return n_matrix_type;
    case d_matrix_type:			// See if can be stored d_matrix
      f = 1;				// Start with f being true
      for(i=0; i<rows_&&f; i++)		// Loop over all the rows of hmx
        {
        pos = i*cols_-(i*(i-1))/2+1; 	// Put pos 1 past diagonal element <i|hmx|i+1>
        cblk = ncols-i-1;		// Non-diagonal cols to examine from row i
        for(j=0; j<cblk&&f; j++, pos++)	// Loop over non-diagonal hmx elements row i
	  f = (norm(data[pos]) < d);	// and insure they are essentially zero
        }
      return f?d_matrix_type:h_matrix_type;
      break;
    case i_matrix_type:			// See if can be stored i_matrix
      f = 1;				// Start with f being true
      for(i=0; i<rows_&&f; i++)		// Check for non-zero off-diagonals
        {
        pos = i*cols_-(i*(i-1))/2; 	// Put pos at diagonal element <i|hmx|i>
	f = (norm(data[pos]-1) < d);	// and insure that the diagonal is 1
        pos++;				// Move to next element
        cblk = ncols-i-1;		// Non-diagonal cols to examine from row i
        for(j=0; j<cblk&&f; j++, pos++)	// Loop over non-diagonal hmx elements row i
	  f = (norm(data[pos]) < d);	// at insure they are essentially zero
        }
      return f?i_matrix_type:h_matrix_type;
      break;
    default:
      return h_matrix_type;
      break;
    }
  return h_matrix_type;
  }
 
std::string h_matrix::mxtype() const { return std::string("Hermitian"); }
std::string h_matrix::mxtype(bool pf) const 
  {
  if(!pf) return std::string("Hermitian");
  else if(is_real())      return std::string("Real Symmetric");
  else if(is_imaginary()) return std::string("Imaginary Symmetric");
  return std::string("Hermitian");
  }

// ____________________________________________________________________________
// E                  CLASS H_MATRIX VARIOUS MATRIX CHECKS
// ____________________________________________________________________________

/*   Function    Output                       Description
   ------------  ------  ------------------------------------------------------
   is_symmetric   bool   Always TRUE for this matrix type
   is_hermitian   bool   Always TRUE for this matrix type
   is_unitary     bool   TF if inv(hmx) == adjoint mx, CPU intensive
   is_real        bool   TF if Im(<i|hmx|j>) < d for all i,j (def. d GMxCut)
   is_imaginary   bool   TF if Re(<i|hmx|j>) < d for all i,j (def. d GMxCut)
   is_complex     bool   TF if is_real && is_imaginary       (def. d GMxCut)
   is_zero        bool   TF if ||<i|hmx|j>|| < d for all i,j (def. d GMxCut)
   is_diagonal    bool   TF if ||<i|hmx|j>|| < d for all i!=j(def. d GMxCut)
   is_square      bool   TF if rows_ == cols_                (TRUE)
   is_equal       bool   TF if ||<i|hmx-mx|j>||<d    all i,j (def. d GMxCut)

Note that the unitary check herein tests if the array inverse is equal to its
adjoint.  So, the routine multiplies the matrix by its adjoint and then looks
to see how well the result conforms to an identity matrix.  This can use up
lots of time....

            ---                         ---
            \                   *t      \
   z(i,j) = /    <i|hmx|k><k|hmx  |j> = /    <i|hmx|k><k|hmx|j> = del +/- d
            ---                         ---                          ij
             k                           k                                  */

bool h_matrix::is_symmetric(double d) const { return is_real(d);  }
bool h_matrix::is_hermitian(double d) const { return true; d = 7; }

bool h_matrix::is_unitary(double d) const
  {
  bool f = true;				// Start with assumed unitary
  int nr = rows_;				// Dimension of hmx
  complex z;					// Complex number used in sum
  int i,j,k;					// Indices used in sum
  int ii, ik, ki, jj, jk, kj;			// Index equivs. of i,j,k pairs
  for(i=0; i<nr && f; i++)	                // Sum over all the matrix rows
    {
    ii = i*nr-(i*(i-1))/2;			// Diag.Elem. data[ii]=<i|hmx|i>
    for(j=i; j<nr && f; j++)			// Sum up.triang. of prdt (i<=j)
      {
      z = 0;					// Set <i|pdt|j> = 0
      jj = j*nr-(j*(j-1))/2;			// data[jj]=<j|hmx|j>
      ki = i;					// data[ki]=<k|hmx|i>@<0|hmx|i>
      kj = j;					// data[kj]=<k|hmx|j>@<0|hmx|j>
      for(k=0; k<i; k++,ki+=nr-k,kj+=nr-k)	// Sum k: i>k,j>k (l.t. 1st, u.t 2nd)
        z += conj(data[ki])*data[kj];		// 	<i|pdt|j> += <i|hmx|k><k|hmx|j>
      for(ik=ii; k<j; k++, ik++, kj+=nr-k)	// Sum over k: i<=k, k<j  (u.t. 1st, l.t. 2nd)
        z += data[ik]*data[kj];			// 	<i|pdt|j> += <i|hmx|k><k|hmx|j>
      for(jk=jj; k<nr; k++,ik++,jk++)		// Sum over k: i<=k, k>=j (u.t. 1st, u.t. 2nd)
        z += data[ik]*conj(data[jk]);		// 	<i|pdt|j> += <i|hmx|k><k|hmx|j>
      if(i == j)
        {					// Diagonal must be 1
        f = ((fabs(Re(z)-1)) < d);		//	check Re(<i|pdt|i>) is 1
        if(f)					//	and if that is true
          f = (fabs(Im(z)) < d);		//	check Im(<i|pdt|i>) is 0
        }
      else                                      // Off diagonal must be zero
        f = (norm(z) < d);			//	just check the norm here
      }
    }
  return f;
  }

bool h_matrix::is_real(double d) const
  {
  bool f = true;
  for(int r=0; (r<size)&&f; r++ )
    f = (fabs(Im(data[r])) <= d);
  return f;
  }

bool h_matrix::is_imaginary(double d) const
  {
  bool f=true;					// Assume it is imaginary
  for(int pos=0; pos<size && f; pos++)          // Loop over all the elements
    f = (fabs(Re(data[pos])) <= d);             // Stop if any real component
  if(f) f = (!is_zero(d));                      // Make sure this is not the zero mx
  return f;
  }

bool h_matrix::is_complex(double d) const
  {
  bool f=true;
  if(is_imaginary(d)) f = 0; 		// False if pure imaginary
  else if(is_real(d)) f = 0; 		// False if pure real
  return f;
  }

bool h_matrix::is_zero(double d) const
  {
  bool f = true;
  for(int r=0; (r<size)&&f; r++)
    f = (norm(data[r]) < d);
  return f;
  }

bool h_matrix::is_diagonal(double d) const
  {
  bool f = true;				// Assume it is diagonal
  int nr = rows_;				// Dimension of hmx
  complex *hii = data;				// hii = <i|hmx|i> -> <0|hmx|0>
  complex *hij;					// hij = <i|hmx|j>
  complex *hend = hii + (nr*nr+nr)/2;		// Data end: <nr-1|hmx|nr>
  int hrow = nr;                                // Row length
  for(; hii<hend && f; )			// Loop diagonal elements
    {
    hij = hii+1;				//  <i|hmx|j> -> <i|hmx|i+1>
    hii += hrow;				//  <i|hmx|i> -> <i+1|hmx|i+1>
    hrow--;					//  Row length is less now
    for(; hij<hii && f; hij++)
      f = bool(norm(*hij) < d);
    }
  return f;
  }

bool h_matrix::is_tridiagonal(double d) const
  {
  bool f = true;				// Assume it is tri-diagonal
  int nr = rows_;				// Dimension of hmx
  complex *hii = data;				// hii = <i|hmx|i> -> <0|hmx|0>
  complex *hij;					// hij = <i|hmx|j>
  complex *hend = hii + (nr*nr+nr)/2;		// Data end: <nr-1|hmx|nr>
  int hrow = nr;                                // Row length
  for(; hii<hend && f; )			// Loop diagonal elements
    {
    hij = hii+2;				//  <i|hmx|j> -> <i|hmx|i+2>
    hii += hrow;				//  <i|hmx|i> -> <i+1|hmx|i+1>
    hrow--;					//  Row length is less now
    for(; hij<hii && f; hij++)
      f = bool(norm(*hij) < d);
    }
  return f;
  }

bool h_matrix::is_square( ) const { return true; }
 
bool h_matrix::is_equal(_matrix* mx, double d) const
  {
  if(cols_ != mx->cols()) return false;			// We first check that
  if(rows_ != mx->rows()) return false;			// mx has the same size
  bool flag = mx->is_hermitian(d);			// and is Hermtiian
  int r, c;						// If so, then we are
  for(r=0; r<rows_&&flag; r++)				// forced to compare all
    for(c=r; c<rows_&&flag; c++)			// elements
       flag=(get(r,c)==(*mx)(r,c));
  return flag;
  }

// ____________________________________________________________________________
// F                  CLASS H_MATRIX ARITHMETIC FUNCTIONS
// ____________________________________________________________________________


_matrix* h_matrix::add(_matrix* mx)

	// Input            hmx : A h_matrix (this)
	//  	             mx : A second matrix
	// Output           sum : The result of adding hmx and mx;
	//			  sum = hmx + mx
	// Note			: Uses internal structure of other matrix types

  {
  if((cols_!=mx->rows()) || (cols_!=mx->cols()))// Check that sizes proper
    {
    HMxerror(5, "add", 1);			// Bad use of function 
    HMxfatal(51);				// Array dimensions mismatched
    }
  switch(mx->stored_type())
    {
    case h_matrix_type:				// Add two h_matrices -> Hermitian
      {
      h_matrix* sum = new h_matrix(cols_,cols_);	// New Hermitian matrix for sum
      complex *mij = &((h_matrix*)mx)->data[size-1];	// <i|mx|j> set to last mx element
      complex *hij = &data[size-1];		    	// <i|hmx|j> set to last hmx element
      complex *sij = &(sum->data[size-1]);	        // <i|sum|j> set at last sum element
      for(; hij>=data; sij--,hij--,mij--)		// Perform the addition elementwise
	(*sij) = (*hij) + (*mij);			// from the last to the first
      return sum;
      }
      break;
    case d_matrix_type:				// Add d_matrix to h_matrix
      {
      if(!mx->is_real())			// If dmx complex, sum not Hermitian
        {						// so the returned matrix will be n_matrix
        n_matrix* sum = new n_matrix(cols_,cols_);	// Generate new n_matrix for sum
        complex *dii = ((d_matrix*)mx)->data;		// Element of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *hij = data;				// Element of hmx: hij = <i|hmx|j> -> <0|hmx|0>
        complex *sij = ((n_matrix*)sum)->data;	// Element of sum: sij = <i|sum|j> -> <0|sum|0>
        complex *sji;					// Element of sum: sji = <j|sum|i>
        int i,j;
        for(i=0; i<cols_; i++,sij+=i,dii++)		// Loop over rows i, dii=<i|dmx|i>, sij=<i|sum|i>
          {
          *sij = *dii + *hij;				// <i|sum|i> = <i|dmx|i> + <i|hmx|i>
          sji = sij+cols_;				// <j|sum|i> -> <i+1|sum|i>
          for(j=i+1,hij++,sij++; j<cols_;		// Loop over upper-triangle elements
                          j++, hij++,sij++,sji+=cols_)	// Start at hij=<i|hmx|i+1>, sij=<i|sum|i+1>
            {
            *sij = *hij;				// <i|sum|j> = <i|mx|j>, i>j
            *sji = conj(*hij);			// <j|sum|i> = <i|mx*|j>, i<j
            }
          }
        return sum;
	}
      else						// For real dmx, sum will be Hermitian 
        {
	h_matrix* sum = new h_matrix(*this);		// New h_matrix equal to this
        for(int i=0; i<cols_; i++)			// Just add in diagonals now
          sum->data[i*cols_-(i*(i-1))/2]
                           += ((d_matrix*)mx)->data[i];
        return sum;
	}
      }
      break;
    case i_matrix_type:				// Add i_matrix to h_matrix -> Hermitian
      {
      h_matrix* sum = new h_matrix(*this);	// New h_matrix equal to hmx
      for(int i=0; i<cols_; i++)		// Just add 1 to its diagonal
        sum->data[i*cols_-(i*(i-1))/2] += complex1;
      return sum;
      }
      break;
    case n_matrix_type:				// Add n_matrix to h_matrix->n_matrix
      {
      n_matrix* sum = new n_matrix(cols_,cols_);// Construct new empty matrix sum
      complex *hij = data;			// Element of hmx: hij = <i|hmx|j> -> <0|hmx|0>
      complex *nij = ((n_matrix*)mx)->data;	// Element of nmx: nij = <i|nmx|j> -> <0|nmx|0>
      complex *mij = ((n_matrix*)sum)->data;    // Element of sum: mij = <i|sum|j> -> <0|sum|0>
      complex *nji, *mji;                       // Elements nji=<j|nmx|i>, mji=<j|sum|i>
      int i,j; 
      for(i=0; i<rows_; i++, nij+=i, mij+=i)    // Loop over the rows
        {
        (*mij) = (*nij) + (*hij);               // <i|sum|j> = <i|nmx|j> + <i|hmx|j>, i=j
        mji = mij + cols_;                      // <j|sum|i> -> <i+1|sum|i>
        nji = nij + cols_;                      // <j|nmx|i> -> <i+1|nmx|i>
        mij++;                                  // <i|sum|j> -> <i|sum|i+1>
        nij++;                                  // <i|nmx|j> -> <i|nmx|i+1>
        hij++;                                  // <i|hmx|j> -> <i|hmx|i+1>
        for(j=i+1; j<cols_; j++,hij++,mij++,nij++)    // Loop over the columns
          {
          (*mij) = (*nij) + (*hij);             // <i|sum|j> = <i|nmx|j> + <i|hmx|j>, i<j
          (*mji) = (*nji) + conj(*hij);         // <j|sum|i> = <j|nmx|i> + <i|hmx*|j>, i<j
          mji += cols_;                         // <j|sum|i> -> <j+1|sum|i>
          nji += cols_;                         // <j|nmx|i> -> <j+1|nmx|i>
          }
        }
      return sum;
      }
      break;
    default:						// Add generic matrix to h_matrix
      {
      n_matrix* sum = new n_matrix(cols_,cols_);	// Sructure we must return n_matrix
      int j, pos = 0;					// and use access functions
      int pnmxup, pnmxlw;				// Upper, lower triangle positions of sum
      for(int i=0; i<cols_; i++)			// Loop over all rows
        {
        pnmxup = i*cols_;				// Position of diagonal <i|sum|i>
        pnmxlw = pnmxup;				// Position of diagonal <i|sum|i>
        for(j=i; j<cols_; j++, pnmxup++, pos++)	// Loop over upper triangle elements
          {
          sum->data[pnmxup] = (*mx)(i,j) + data[pos];	// <i|sum|j> = <i|mx|j>+<i|hmx|j>, i<=j
          if(i != j)
            {
            pnmxlw--;
            sum->data[pnmxlw] = 			// <j|sum|i> = <j|mx|i>+<i|hmx|j>*, i>j
                           (*mx)(j,i) + conj(data[pos]);
            }
          }
	}
      return sum;
      }
    break;
    }
  }


_matrix* h_matrix::subtract(_matrix* mx)

	// Input            hmx : A h_matrix (this)
	//  	             mx : A second matrix
	// Output           mx1 : The result of subtracting mx from hmx
	// Note			: Uses internal structure of other matrix types

  {
  if((cols_!=mx->rows()) || (cols_!=mx->cols())) 	// Check that sizes proper
    {
    HMxerror(5, "subtract", 1);			// Bad use of function 
    HMxerror(51);				// Array dimensions mismatched
    return mx;
    }
  switch(mx->stored_type())
      {
      case h_matrix_type:				// Subtract two h_matrices
	{
        h_matrix* sub = new h_matrix(cols_,cols_);	// New Hermitian matrix sub
        complex *mij = &((h_matrix*)mx)->data[size-1];	// Pointer to last mx element
        complex *hij = &data[size-1];		    	// Pointer to last this element
        complex *sij = &(sub->data[size-1]);	        // Pointer to last sub element
        for(; hij>=data; sij--,hij--,mij--)		// Perform the subtraction elementwise
	  (*sij) = (*hij) - (*mij);			// from the last to the first
        return sub;
        }
        break;
      case d_matrix_type:				// Subtract d_matrix from h_matrix
        {						// Hermitian if dmx is real, normal else 
        if(!mx->is_real())
          {
          n_matrix* sub = new n_matrix(cols_,cols_);	// Generate new n_matrix for sub
          complex *dii = ((d_matrix*)mx)->data;		// Element of dmx: dii = <i|dmx|i> -> <0|dmx|0>
          complex *hij = data;				// Element of hmx: hij = <i|hmx|j> -> <0|hmx|0>
          complex *sij = ((n_matrix*)sub)->data;	// Element of sub: sij = <i|sub|j> -> <0|sub|0>
          complex *sji;					// Element of sub: sji = <j|sub|i>
          int i,j;
          for(i=0; i<cols_; i++,sij+=i,dii++)		// Loop rows i, dii=<i|dmx|i>, sij=<i|sub|i>
            {
            *sij = *hij - *dii;				// <i|sub|i> = <i|hmx|i> - <i|dmx|i>
            sji = sij+cols_;				// <j|sub|i> -> <i+1|sub|i>
            for(j=i+1,hij++,sij++; j<cols_;		// Loop over upper-triangle elements
                          j++, hij++,sij++,sji+=cols_)	// Start at hij=<i|hmx|i+1>, sij=<i|sub|i+1>
              {
              *sij = *hij;				// <i|sub|j> = <i|mx|j>, i>j
              *sji = conj(*hij);			// <j|sub|i> = <i|mx*|j>, i<j
              }
            }
	  return sub;
	  }
        else 
          {
	  h_matrix* mx1 = new h_matrix(*this);		// New h_matrix equal to this
          for(int i=0; i<cols_; i++)			// Just add in diagonals now
            mx1->data[i*cols_-(i*(i-1))/2]
                           -= ((d_matrix*)mx)->data[i];
          return mx1;
	  }
        }
        break;
      case i_matrix_type:				// Subtract i_matrix from h_matrix
        { 						// -> still Hermitian
        h_matrix* dif = new h_matrix(*this);		// New h_matrix equal to hmx
        for(int i=0; i<cols_; i++)			// Just subtract 1 to its diagonal
          dif->data[i*cols_-(i*(i-1))/2] -= complex1;
        return dif;
        }
        break;
      case n_matrix_type:				// Subtract n_matrix from h_matrix->n_matrix
        {
        n_matrix* sub = new n_matrix(cols_,cols_);      // Construct new empty matrix sub
        complex *hij = data;				// Element of hmx: hij = <i|hmx|j> -> <0|hmx|0>
        complex *nij = ((n_matrix*)mx)->data;		// Element of nmx: nij = <i|nmx|j> -> <0|nmx|0>
        complex *mij = ((n_matrix*)sub)->data;          // Element of sub: mij = <i|sub|j> -> <0|sub|0>
        complex *nji, *mji;                             // Elements nji=<j|nmx|i>, mji=<j|sub|i>
        int i,j; 
        for(i=0; i<rows_; i++, nij+=i, mij+=i)          // Loop over the rows
          {
          (*mij) = (*hij) - (*nij);                     // <i|sub|j> = <i|hmx|j> - <i|nmx|j>, i=j
          mji = mij + cols_;                            // <j|sub|i> -> <i+1|sub|i>
          nji = nij + cols_;                            // <j|nmx|i> -> <i+1|nmx|i>
          mij++;                                        // <i|sub|j> -> <i|sub|i+1>
          nij++;                                        // <i|nmx|j> -> <i|nmx|i+1>
          hij++;                                        // <i|hmx|j> -> <i|hmx|i+1>
          for(j=i+1; j<cols_; j++,hij++,mij++,nij++)    // Loop over the columns
            {
            (*mij) = (*hij) - (*nij);                   // <i|sub|j> = <i|hmx|j> - <i|nmx|j>, i<j
            (*mji) = conj(*hij) - (*nji);		// <j|sub|i> = <i|hmx*|j> - <j|nmx|i>, i<j
            mji += cols_;                               // <j|sub|i> -> <j+1|sub|i>
            nji += cols_;                               // <j|nmx|i> -> <j+1|nmx|i>
            }
          }
        return sub;
        }
        break;
      default:						// Subtract generic matrix from h_matrix
	{ 
        n_matrix* mx1 = new n_matrix(cols_,cols_);	// We must return n_matrix
        int j, pos = 0;					// and are forced to use access functions
        int pnmxup, pnmxlw;				// Upper, lower triangle positions of mx1
        for(int i=0; i<cols_; i++)			// Loop over all rows
          {
          pnmxup = i*cols_;				// Position of diagonal <i|mx1|i>
          pnmxlw = pnmxup;				// Position of diagonal <i|mx1|i>
          for(j=i; j<cols_; j++, pnmxup++, pos++)	// Loop over upper triangle elements
            {
            mx1->data[pnmxup] = data[pos] - (*mx)(i,j);	// <i|mx1|j> = <i|hmx|j>-<i|mx|j>, i<=j
            if(i != j)
              {
              pnmxlw--;
              mx1->data[pnmxlw] = 			// <j|mx1|i> = <i|hmx|j>*-<j|mx|i>, i>j
                          conj(data[pos]) - (*mx)(j,i);
              }
            }
	  }
	return mx1;
	}
      }
    }


_matrix* h_matrix::multiply(_matrix* mx)

	// Input            hmx : A h_matrix (this)
	//  	             mx : A second matrix
	// Output           mx1 : The result of multiplying mx into hmx
	// Note			: This is hmx*mx, != mx*hmx unless they commute
	// Note			: Uses internal structure of other matrix types

/*                                ---
                      <i|pdt|j> = \   <i|hmx|k> <k|mx|j>
                                  /
                                  ---
*/

  {
  if(cols_ != mx->rows())			// Insure dimensions match up
    { 						// hmx(cols_Xcols_)*(mx->rows()Xmx->cols())
    HMxerror(5, "multiply", 1);			// Bad use of function 
    HMxfatal(51);				// Array dimensions mismatched
    }
  switch(mx->stored_type())
    {
    case i_matrix_type: return this; break;	// Nothing if I*hmx, return hmx
    case h_matrix_type:				// Multiply h_matrix into h_matrix
      { 
#ifdef _USING_BLAS_
      int A_rows = rows();
      int A_cols = cols();            
      int B_rows = mx->rows();
      int B_cols = mx->cols();
      int C_rows = A_rows;
      int C_cols = B_cols;
      n_matrix* pdt =	new n_matrix(rows_, cols_);
      if(C_rows * C_cols > 16)  //4*4 = 16
      {
        double alpha[2] = {0,0};
        double beta[2] = {0,0};
        alpha[0] = 1.0;
        beta[0]  = 0.0;
        n_matrix* hmxA =	new n_matrix(A_rows,A_cols);		// Create new matrix h_matrix
        n_matrix* hmxB =	new n_matrix(B_rows,B_cols);		// Create new matrix h_matrix
        this->convert(hmxA);				// Convert h_matrix mx into normal matrix hmx
        mx->convert(hmxB);				// Convert h_matrix mx into normal matrix hmx
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A_rows, B_cols, A_cols, 
                         alpha, (double *)hmxA->data, A_cols, (double *)hmxB->data, B_cols, beta, (double *)pdt->data, C_cols);
//      std::cerr << "BLAS: h_matrix * h_matrix\n";
        delete hmxA;
        delete hmxB;
      }
      else
      {
          complex *p00 = pdt->data;			// Start of pdt: p00 = <0|pdt|0>
          complex *h00 = data;			// Start of hmx: h00 = <0|hmx|0>
          complex *m00 = ((h_matrix*)mx)->data;	// Start of mx: m00 = <0|mx|0>
          complex *pij, *hik, *hki, *h0i;
          complex *m0j, *mkj, *mjk;
          int i,j,k;
          for(i=0,pij=p00,h0i=h00; i<rows_; i++,h0i++)	// Loop over the rows of product
            {
	    for(j=0, m0j=m00; j<=i; j++,pij++,m0j++)	// Loop over the lower triangle of pdt (j<=i)
              {
              (*pij) = complex0;				// Initialize element: <i|pdt|j> = 0
	      for(k=0,hki=h0i,mkj=m0j; k<j;		// Begin looping over inner index k
                          k++,hki+=cols_-k, mkj+=cols_-k)
                  (*pij) += conj(*hki) * (*mkj);		// <i|pdt|j> += <k|hmx*|i><k|mx|j>, k<i, k<j
	        for(mjk=mkj; k<i; k++,hki+=cols_-k,mjk++)
                  (*pij) += conj(*hki) * conj(*mjk);	// <i|pdt|j> += <k|hmx*|i><j|mx*|k>, k<i, k>=j
	        for(hik=hki; k<cols_; k++,hik++,mjk++)
                  (*pij) += (*hik) * conj(*mjk); 		// <i|pdt|j> += <i|hmx|k><j|mx*|k>, k>=i, k>=j
                }
	      for(; j<cols_; j++,pij++,m0j++)		// Loop over the upper triangle of pdt (j>i) 
                {
                (*pij) = complex0;				// Initialize element: <i|pdt|j> = 0
	        for(k=0,hki=h0i,mkj=m0j; k<i;		// Begin looping over inner index k
                          k++,hki+=cols_-k, mkj+=cols_-k)
                  (*pij) += conj(*hki) * (*mkj);		// <i|pdt|j> += <k|hmx*|i><k|mx|j>, k<i, k<j
	        for(hik=hki; k<j; k++,hik++,mkj+=cols_-k)
                  (*pij) += (*hik) * (*mkj); 		// <i|pdt|j> += <i|hmx|k><k|mx|j>, k>=i, k<j
	        for(mjk=mkj; k<cols_; k++,hik++,mjk++)
                  (*pij) += (*hik) * conj(*mjk); 		// <i|pdt|j> += <i|hmx|k><j|mx*|k>, k>=i, k>=j
                }
              }
//      std::cerr << "NO BLAS: h_matrix * h_matrix\n";
      }
#else
      n_matrix* pdt = new n_matrix(rows_,cols_);// Construct a new normal matrix
      complex *p00 = pdt->data;			// Start of pdt: p00 = <0|pdt|0>
      complex *h00 = data;			// Start of hmx: h00 = <0|hmx|0>
      complex *m00 = ((h_matrix*)mx)->data;	// Start of mx: m00 = <0|mx|0>
      complex *pij, *hik, *hki, *h0i;
      complex *m0j, *mkj, *mjk;
      int i,j,k;
      for(i=0,pij=p00,h0i=h00; i<rows_; i++,h0i++)	// Loop over the rows of product
        {
	for(j=0, m0j=m00; j<=i; j++,pij++,m0j++)	// Loop over the lower triangle of pdt (j<=i)
          {
          (*pij) = complex0;				// Initialize element: <i|pdt|j> = 0
	  for(k=0,hki=h0i,mkj=m0j; k<j;		// Begin looping over inner index k
                      k++,hki+=cols_-k, mkj+=cols_-k)
              (*pij) += conj(*hki) * (*mkj);		// <i|pdt|j> += <k|hmx*|i><k|mx|j>, k<i, k<j
	    for(mjk=mkj; k<i; k++,hki+=cols_-k,mjk++)
              (*pij) += conj(*hki) * conj(*mjk);	// <i|pdt|j> += <k|hmx*|i><j|mx*|k>, k<i, k>=j
	    for(hik=hki; k<cols_; k++,hik++,mjk++)
              (*pij) += (*hik) * conj(*mjk); 		// <i|pdt|j> += <i|hmx|k><j|mx*|k>, k>=i, k>=j
            }
	  for(; j<cols_; j++,pij++,m0j++)		// Loop over the upper triangle of pdt (j>i) 
            {
            (*pij) = complex0;				// Initialize element: <i|pdt|j> = 0
	    for(k=0,hki=h0i,mkj=m0j; k<i;		// Begin looping over inner index k
                      k++,hki+=cols_-k, mkj+=cols_-k)
              (*pij) += conj(*hki) * (*mkj);		// <i|pdt|j> += <k|hmx*|i><k|mx|j>, k<i, k<j
	    for(hik=hki; k<j; k++,hik++,mkj+=cols_-k)
              (*pij) += (*hik) * (*mkj); 		// <i|pdt|j> += <i|hmx|k><k|mx|j>, k>=i, k<j
	    for(mjk=mkj; k<cols_; k++,hik++,mjk++)
              (*pij) += (*hik) * conj(*mjk); 		// <i|pdt|j> += <i|hmx|k><j|mx*|k>, k>=i, k>=j
            }
          }
#endif
	return pdt;
	}
        break;
    case n_matrix_type:				// Multiply h_matrix into an n_matrix
	{ 
#ifdef _USING_BLAS_
        int A_rows = rows();
        int A_cols = cols();            
        int B_cols = mx->cols();
        int C_rows = A_rows;
        int C_cols = B_cols;
        n_matrix* pdt =	new n_matrix(C_rows, C_cols);
        if(C_rows * C_cols > 16)  //4*4 = 16
        {
          double alpha[2] = {0,0};
          double beta[2] = {0,0};
          alpha[0] = 1.0;
          beta[0]  = 0.0;
          n_matrix* hmxA =	new n_matrix(A_rows,A_cols);		// Create new matrix h_matrix
          this->convert(hmxA);				// Convert h_matrix mx into normal matrix hmx
	  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A_rows, B_cols, A_cols, 
                           alpha, (double *)hmxA->data, A_cols, (double *)((n_matrix*)mx)->data, B_cols, beta, (double *)pdt->data, C_cols);
//        std::cerr << "BLAS: h_matrix * n_matrix\n";
	  delete hmxA;
        }
        else
        {
	  int c = mx->cols();				// Columns of product matrix
	  complex *p00 = pdt->data;			// Start of pdt: p00 = <0|pdt|0>
	  complex *h00 = data;				// Start of hmx: h00 = <0|hmx|0>
	  complex *n00 = ((n_matrix*)mx)->data;		// Start of nmx: n00 = <0|nmx|0>
	  complex *n10 = n00 + c;				// End of 1st nmx row: n10 = <1|nmx|0>
	  complex *pij, *hik, *hki, *h0i, *n0j, *nkj;
          int i, k;
	  for(i=0,pij=p00,h0i=h00; i<rows_; i++,h0i++)
	    for(n0j=n00; n0j<n10; pij++,n0j++)		// Effective loop over j
              {
              (*pij) = complex0;				// Initialize <i|pdt|j> to zero
	      for(k=0,nkj=n0j,hki=h0i; k<i;		// Begin looping over k
                                k++,nkj+=c,hki+=cols_-k)
                (*pij) += conj(*hki) * (*nkj);		// <i|pdt|j> += <k|hmx*|i><k|nmx|j>, i>k
	      for(hik=hki; k<cols_; k++,hik++,nkj+=c)
                (*pij) += (*hik) * (*nkj);		// <i|pdt|j> += <i|hmx|k><k|nmx|j>, i<=k
            }
//        std::cerr << "NON BLAS: h_matrix * n_matrix\n";
	}
#else
	int c = mx->cols();				// Columns of product matrix
	n_matrix* pdt =	new n_matrix(rows_,c);		// Construct a new normal matrix
	complex *p00 = pdt->data;			// Start of pdt: p00 = <0|pdt|0>
	complex *h00 = data;				// Start of hmx: h00 = <0|hmx|0>
	complex *n00 = ((n_matrix*)mx)->data;		// Start of nmx: n00 = <0|nmx|0>
	complex *n10 = n00 + c;				// End of 1st nmx row: n10 = <1|nmx|0>
	complex *pij, *hik, *hki, *h0i, *n0j, *nkj;
        int i, k;
	for(i=0,pij=p00,h0i=h00; i<rows_; i++,h0i++)
	  for(n0j=n00; n0j<n10; pij++,n0j++)		// Effective loop over j
            {
            (*pij) = complex0;				// Initialize <i|pdt|j> to zero
	    for(k=0,nkj=n0j,hki=h0i; k<i;		// Begin looping over k
                              k++,nkj+=c,hki+=cols_-k)
              (*pij) += conj(*hki) * (*nkj);		// <i|pdt|j> += <k|hmx*|i><k|nmx|j>, i>k
	    for(hik=hki; k<cols_; k++,hik++,nkj+=c)
              (*pij) += (*hik) * (*nkj);		// <i|pdt|j> += <i|hmx|k><k|nmx|j>, i<=k
            }
#endif
	return pdt;
	}
	break;
    case d_matrix_type:				// Multiply h_matrix into d_matrix
	{
	n_matrix* pdt =	new n_matrix(cols_,cols_); 	// Construct a new normal matrix
	complex *p00 = pdt->data; 			// Start of pdt: p00 = <0|pdt|0>
	complex *h00 = data;				// Start of hmx: h00 = <0|hmx|0>
	complex *d00 = ((d_matrix*)mx)->data;		// Start of dmx: d00 = <0|dmx|0>
        complex *dend = d00 + cols_;			// End of dmx: dend = <cols_|dmx|cols_+1>
	complex *pij, *h0i, *hij, *hji, *djj, *dii;	// These will be actual matrix elements
        int hrow;					// hji -> hji+hrow: <j|hmx|i> -> <j+1|hmx|i>
	for(dii=d00,h0i=h00,pij=p00; dii<dend;		// Effective loop over i
                                           h0i++,dii++)
          {
	  for(djj=d00,hji=h0i,hrow=cols_; djj<dii;	// Effective loop over j
			  pij++,djj++,hrow--,hji+=hrow)
            (*pij) = conj(*hji) * (*djj);		// <i|pdt|j> = <j|hmx*|i><j|dmx|j>, i>j
	  for(hij=hji; djj<dend; pij++,hij++,djj++)
            (*pij) = (*hij) * (*djj);			// <i|pdt|j> = <i|hmx|j><j|dmx|j>, i<=j
          }
//      std::cerr << "NON BLAS: h_matrix * d_matrix\n";
	return pdt;
	}
	break;
    default:						// Mult. generic mx into h_matrix
      { 
#ifdef _USING_BLAS_
        int A_rows = rows();
        int A_cols = cols();            
        int B_cols = mx->cols();
        int C_rows = A_rows;
        int C_cols = B_cols;
        n_matrix* pdt =	new n_matrix(C_rows, C_cols);
        if(C_rows * C_cols > 16)  //4*4 = 16
        {
          double alpha[2] = {0,0};
          double beta[2] = {0,0};
          alpha[0] = 1.0;
          beta[0]  = 0.0;
          n_matrix* hmxA =	new n_matrix(A_rows,A_cols);		// Create new matrix h_matrix
          this->convert(hmxA);				// Convert h_matrix mx into normal matrix hmx
	  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A_rows, B_cols, A_cols, 
                           alpha, (double *)hmxA->data, A_cols, (double *)((n_matrix*)mx)->data, B_cols, beta, (double *)pdt->data, C_cols);
//        std::cerr << "BLAS: h_matrix * unknown\n";
	  delete hmxA;
        }
        else
        {
          int c = mx->cols();				// Column dimension of product
          complex *p00 = pdt->data;			// Start of pdt: p00 = <0|pdt|0>
          complex *h00 = data;				// Start of hmx: h00 = <0|hmx|0>
          complex *pij, *hik, *hki, *h0i;
          int i,j,k;
          for(i=0,pij=p00,h0i=h00; i<rows_; i++,h0i++)	// Loop over the rows of product
            {
	    for(j=0; j<c; j++,pij++)			// Loop over the columns of product
              {
              (*pij) = complex0;			// Initialize element <i|pdt|j> to zero
	      for(k=0,hki=h0i; k<i; k++,hki+=cols_-k)	// Begin loop over inner index k
                  (*pij) += conj(*hki) * (*mx)(k,j);	// <i|pdt|j> += <k|hmx*|i><k|mx|j>, i>k
	      for(hik=hki; k<cols_; k++, hik++)
                  (*pij) += (*hik) * (*mx)(k,j); 		// <i|pdt|j> += <i|hmx|k><k|mx|j>, i<=k
	      }
//          std::cerr << "NON BLAS: h_matrix * unknown\n";
	  }
#else
      int c = mx->cols();				// Column dimension of product
      n_matrix* pdt =	new n_matrix(rows_,c);		// Construct a new normal matrix
      complex *p00 = pdt->data;			// Start of pdt: p00 = <0|pdt|0>
      complex *h00 = data;				// Start of hmx: h00 = <0|hmx|0>
      complex *pij, *hik, *hki, *h0i;
      int i,j,k;
      for(i=0,pij=p00,h0i=h00; i<rows_; i++,h0i++)	// Loop over the rows of product
        {
	for(j=0; j<c; j++,pij++)			// Loop over the columns of product
          {
          (*pij) = complex0;			// Initialize element <i|pdt|j> to zero
	  for(k=0,hki=h0i; k<i; k++,hki+=cols_-k)	// Begin loop over inner index k
              (*pij) += conj(*hki) * (*mx)(k,j);	// <i|pdt|j> += <k|hmx*|i><k|mx|j>, i>k
	  for(hik=hki; k<cols_; k++, hik++)
              (*pij) += (*hik) * (*mx)(k,j); 		// <i|pdt|j> += <i|hmx|k><k|mx|j>, i<=k
	  }
#endif
	return pdt;
	}
      }
    }
  return mx;						// Should not reach here
  }


_matrix* h_matrix::multiply(const complex& z)

	// Input            hmx : A h_matrix (this)
	//  	             z  : A complex number
	// Output           mx1 : The result of multiplying z*hmx

  {
  if(z == complex1)				// If z real, return Hermitian
    return this;
  if(Im(z) == 0)				// If z real, return Hermitian
    {						// else return an n_matrix
    h_matrix* hmx = new h_matrix(rows(),cols());
    for(int i=0; i<size; i++)
      hmx->data[i] = z*data[i];
    return hmx;
    }
  n_matrix* mx = new n_matrix(rows(),cols());
  int pnmxup, pnmxlw;				// Upper, lower triangle positions of mx
  int j, pos=0;
  for(int i=0; i<rows(); i++)			// Loop over all rows
    {
    pnmxup = i*rows()+i;			// Position of diagonal <i|mx|i>
    pnmxlw = pnmxup;				// Position of diagonal <i|mx|i>
    for(j=i; j<cols(); j++, pnmxup++, pos++)	// Loop over upper triangle elements
      {
      mx->data[pnmxup] = z*data[pos]; 		// <i|mx|j> = z*<i|hmx|j>, i<=j
      if(i != j)
        {
        pnmxlw += cols();
        mx->data[pnmxlw] = z*conj(data[pos]); 	// <j|mx|i> = z*<i|hmx*|j>*, i>j
        }
      }
    }
  return mx;
  }


      
        // Input                hmx  : Input Hermitian matrix (this)
        //                       mx  : Second matrix                         -1
        // Output               pdt  : New matrix which is the product hmx*mx
        // Note                      : This uses the function inv of mx
 
/*                    ---                -1                    -1
          <i|pdt|j> = \   <i|hmx|k> <k|mx  |j> = <i|hmx|i><i|mx  |j>
                      /
                      ---                                                    */
_matrix* h_matrix::divide(_matrix* mx)
  {
  if(cols_ != mx->cols())			// Insure dimensions match
    { 
    HMxerror(5, "divide", 1);			// Bad use of function 
    HMxfatal(51);				// Array dimensions mismatched
    }
  switch(mx->stored_type())
    {
    case i_matrix_type:				// Divide i_matrix into h_matrix
      return this;				// No change, just return hmx
      break;
    case d_matrix_type:				// Divide h_matrix by d_matrix
      return multiply((d_matrix*)mx->inv());	// Return hmx*inv(dmx)
      break;
    case h_matrix_type:				// Divide h_matrix by h_matrix
      return multiply((h_matrix*)mx->inv());	// Return hmx*inv(hmx)
      break;
    case n_matrix_type:				// Divide h_matrix by n_matrix
      return multiply((n_matrix*)mx->inv());	// Return hmx*inv(nmx)
      break;
    default:					// Divide d_matrix by generic mx
      break;
    }
  HMxerror(25, "divide", 1);			// Bad use of function 
  HMxfatal(23);					// Array dimensions mismatched
  return mx;
  } 


_matrix* h_matrix::divide(const complex& z)

	// Input            hmx : A h_matrix (this)
	//  	             z  : A complex number
	// Output           mx1 : The result of multiplying (1/z)*hmx
	// Note			: Uses internal structure of other matrix types

  {
  if(z == complex1) return this; 			// If z=1, return hmx
  if(z == complex0)					// If z=0, problems
    {
    HMxerror(18, 1);					//   Division by zero
    HMxerror(6, " mx/z", 1);				//   Matrix division
    HMxfatal(3, "divide");				//   Fail in divide function
    }
  complex z1 = 1/z;
  if(Im(z) == 0)					// If z real, return Hermitian
    {
    h_matrix* pdt = new h_matrix(rows_,cols_);
    complex *pij = pdt->data;	 			// Start of pdt: <i|pdt|j> -> <0|pdt|0>
    complex *hij = data;				// Start of hmx: <i|hmx|j> -> <0|hmx|0>
    complex *hend = hij + (cols_*cols_+cols_)/2;	// End of data in hmx: <cols_|hmx|cols_+1>
    for(; hij<hend; pij++,hij++)
      *pij = z1 * (*hij);				// <i|pdt|j> = (1/z) * <i|hmx|j>
    return pdt;
    }
  n_matrix* nmx = new n_matrix(rows_,cols_); 		// If z complex, return an n_matrix
  complex *nij = nmx->data;	 			// Start of nmx: <i|nmx|j> -> <0|nmx|0>
  complex *hij = data;					// Start of hmx: <i|hmx|j> -> <0|hmx|0>
  complex *hend = hij + (cols_*cols_+cols_)/2;		// End of data in hmx: <cols_|hmx|cols_+1>
  complex *nend = nij + cols_*cols_;			// End of data in nmx: <cols_|nmx|cols_+1>
  complex *nji;
  for(int i=0; hij<hend; i++, nij+=i)			// Loop over index i (rows of hmx)
    {							// <i|hmx|j>=<i|hmx|i>, <i|nmx|j>=<i|nmx|i>
    (*nij) = z1 * (*hij);				// <i|nmx|i> = (1/z) * <i|hmx|i> 
    nji = nij + cols_;					// <j|nmx|i> -> <i+1|nmx|i> (for j=i+1)
    nij++;						// <i|nmx|j> -> <i|nmx|i+1> (for j=i+1)
    hij++;						// <i|hmx|j> -> <i|hmx|i+1> (for j=i+1)
    for(; nji<nend; hij++,nij++,nji+=cols_)		// Effective loop over j (via nji), j>i
     {
     (*nij) = z1 * (*hij);				// <i|nmx|j> = (1/z) * <i|hmx|j>, j>i 
     (*nji) = z1 * conj(*hij);				// <j|nmx|i> = (1/z) * <i|hmx*|j>, j>i 
     }
   }
  return nmx;
  }
    
// ____________________________________________________________________________
// G              CLASS H_MATRIX UNARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________


	// Input            hmx : A h_matrix (this)
	//		     mx : A second matrix
	// Output           mx1 : The result of adding hmx & mx
	// Note			: This function avoids any matrix
	//			  copying if it can help it

_matrix* h_matrix::add_two(_matrix* mx)
  {
  if((mx->stored_type()) == h_matrix_type)	 // Add h_matrix into h_matrix
    { 
    if(rows_ != mx->rows())			// Insure matrix dimensions match
      {
      HMxerror(17, 1);				// Dimension problems
      HMxerror(6, "mx1 + mx2", 1);		// Bad use of function 
      HMxfatal(3, "add_two");			// Array dimensions mismatched
      }
    else
      {
      complex *a=&((h_matrix*)mx)->data[size-1];// Pointer to last element of mx
      complex *b = &data[size-1];		// Pointer to last element of hmx
      for(; b>=data; b--,a--)			// Perform addition elementwise
        (*a) += (*b);
      return mx;
      }
    }
  return add(mx); 				// If add hmx & !hmx, type changes.
  } 						// Needs copying, just use add function


	// Input            hmx : A h_matrix (this)
	//		     mx : A second matrix
	// Output           mx1 : The result of subtracting hmx from mx
	// Note			: This function avoids any matrix
	//			  copying if it can help it

_matrix* h_matrix::subtract_two(_matrix* mx)
  {
  if ((mx->stored_type()) == h_matrix_type)	 // Subtract h_matrix from h_matrix
    { 
    if(rows() != mx->rows())
      {
      HMxerror(17, 1);				// Dimension problems
      HMxerror(6, "mx1 - mx2", 1);		// Bad use of function 
      HMxfatal(3, "subtract_two");		// Array dimensions mismatched
      }
    else
      {
      complex *a=&((h_matrix*)mx)->data[size-1];// Pointer to last element of mx
      complex *b=&data[size-1];			// Pointer to last element of this
      for(; b>=data; b--,a--)			// Subtraction elementwise
        (*a) -= (*b);
      return mx;
      }
    }
  return mx->subtract(this);			// matrix type changes, result has to be
  } 						// copied anyway -> use subtract function


_matrix* h_matrix::multiply_two(_matrix* mx) 	// In this case, mx cannot be used for the
						// for the output.  If mx=h_matrix then its
  { return mx->multiply(this); }		// not allowed and if not h_matrix, type must
						// by modified (usually to n_matrix)!


_matrix* h_matrix::multiply_two(const complex &z)

  {
  if(z == complex1) return this;                   // If z is 1, don't do anything
  if(Im(z) == 0)                                // If z real, stay Hermitian
    {
    for(int i=0; i<size; i++)
      data[i] = z*data[i];
    return this;
    } 
  n_matrix* mx = new n_matrix(rows(),cols());   // If z has an imaginary component
  int pnmxup, pnmxlw;                           // Upper, lower triangle positions of mx
  int j, pos=0;
  for(int i=0; i<rows(); i++)                   // Loop over all rows
    { 
    pnmxup = i*rows()+i;                        // Position of diagonal <i|mx|i>
    pnmxlw = pnmxup;                            // Position of diagonal <i|mx|i>
    for(j=i; j<cols(); j++, pnmxup++, pos++)    // Loop over upper triangle elements
      {  
      mx->data[pnmxup] = z*data[pos];           // <i|mx|j> = z*<i|hmx|j>, i<=j
      if(i != j)
        {
        pnmxlw += cols();
        mx->data[pnmxlw] = z*conj(data[pos]);   // <j|mx|i> = z*<i|hmx*|j>*, i>j
        }
      }
    }
  return mx;
  } 

_matrix* h_matrix::divide_two(_matrix* mx)      { return mx->divide(this);  }
_matrix* h_matrix::divide_two(const complex &z) { return multiply_two(1/z); }


// ____________________________________________________________________________
// H                CLASS H_MATRIX SIMPLE UNARY FUNCTIONS
// ____________________________________________________________________________
 
/* These functions perform simple simple mathematical operations on a Hermitian
   matrix. Note that for Hermitian arrays the adjoint function does nothing.
 
    Function          Result             Function           Result
   ---------- -----------------------    --------  ------------------------
   negate     <i|hmx|j> -> -<i|hmx|j>       RE      <i|hmx|j>->Re(<i|hmx|j>)
   conjugate  <i|hmx|j> -> <i|hmx*|j>       IM      <i|hmx|j>->Im(<i|hmx|j>)
   transpose  <i|hmx|j> -> <j|hmx|i>      adjoint   <i|hmx|j>-><j|hmx*|i>   */


_matrix* h_matrix::negate()
  {
  h_matrix* hmx1 = new h_matrix(cols_, rows_);
  for(int i=0; i<size; i++) hmx1->data[i] = -data[i];
  return hmx1;
  }

_matrix* h_matrix::RE()
  {
  if(is_real()) return this;
  h_matrix* hmx1 = new h_matrix(cols_, rows_);
  for(int i=0; i<size; i++) hmx1->data[i] = Re(data[i]);
  return hmx1;
  }
 
 _matrix* h_matrix::IM()
   {
   if(is_imaginary()) return this;
   h_matrix* hmx1 = new h_matrix(cols_, rows_);
   for(int i=0; i<size; i++) hmx1->data[i] = Im(data[i]);
   return hmx1;
   }
 
_matrix* h_matrix::conjugate()
  {
  h_matrix* hmx1 = new h_matrix(cols_, rows_);
  for(int i=0; i<size; i++) hmx1->data[i] = conj(data[i]);
  return hmx1;
  }
  
_matrix* h_matrix::transpose()
  {
  h_matrix* hmx1 = new h_matrix(cols_, rows_);
  for(int i=0; i<size; i++) hmx1->data[i] = conj(data[i]);
  return hmx1;
  }

_matrix* h_matrix::adjoint() { return this; }
  
_matrix* h_matrix::mxexp()
  {
  _matrix* dmx;                                 // For eigenvalues
  _matrix* emx;                                 // For eigenvectors
  diag(dmx, emx);                               // Diagonalize mx to dmx & emx
  complex *d00 = ((d_matrix*)dmx)->data;        // Start of dmx: d00 = <0|dmx|0>
  complex *dii, *dend = d00 + cols_;            // Indexing on dmx diagonal
  for(dii=d00; dii<dend; dii++)                 // Effective loop over i
    (*dii) = exp(*dii);                         // Set <i|dmx|i> to exp(<i|dmx|i>
  _matrix* emxi  = emx->inv();                  // Get inverse of emx
  _matrix* pdt1  = dmx->multiply(emxi);
  _matrix* expmx = emx->multiply(pdt1);
  delete dmx;                                   // Clean up the dmx matrix
  delete emx;                                   // Clean up the emx matrix
  delete pdt1;                                  // Clean up the pdt1 matrix
  delete emxi;                                  // Clean up the emxi matrix
  return expmx;
  }

complex h_matrix::trace()

	// Input            hmx : A h_matrix (this)
	// Output           z   : The trace of hmx

  {
  complex z(0);
  for(int i=0; i<rows_; i++)
    z += data[i*cols_-(i*(i-1))/2];
  return z;
  }

/*     Function    Output                       Description
     ------------  -------  ---------------------------------------------------
       swaprows      T/F    Swaps two rows of the matrix, true if successful 
       swapcols      T/F    Swaps two columns of the matrix, true if successful 
        permute    matrix   Permutes rows & columns of the matrix
        maxRe      double   Returns largest real value in the array
        maxIm      double   Returns largest imaginary value in the array
        maxZ       complex  Returns largest complex (norm) value in array
        minRe      double   Returns smallest real value in the array
        minIm      double   Returns smallest imaginary value in the array
        minZ       complex  Returns smallest complex (norm) value in array   */


_matrix* h_matrix::swaprows(int I, int J)
  {
  n_matrix* nsr = NMX();                // Copy hmx to nmx 
  int I0 = I*cols_;                     // Index of <I|mx,mx'|0>
  int J0 = J*cols_;                     // Index of <J|mx,mx'|0>
  for(int k=0; k<cols_; k++)            // Swap rows I & J
    {
    nsr->data[I0+k] = get(J,k);		// Set <I|nmx|k> = <J|hmx|k>
    nsr->data[J0+k] = get(I,k);		// Set <J|nmx|k> = <I|hmx|k>
    }
  return nsr;
  }

_matrix* h_matrix::swapcols(int I, int J)
  {
  n_matrix* nsc = NMX();                // Copy hmx to nmx 
  for(int k=0; k<cols_; k++)            // Swap rows I & J
    {
    nsc->data[k*cols_+I] = get(k, J);    // Set <k|nmx|I> = <k|hmx|J>
    nsc->data[k*cols_+J] = get(k, I);    // Set <k|nmx|J> = <k|hmx|I>
    }
  return nsc;
  }

typedef n_matrix* NPtr;
typedef h_matrix* HPtr;

_matrix* h_matrix::permute(int I, int J)
  {
  h_matrix* phmx = new h_matrix(*this);			// A copy of hmx to use
  int i,j;
  for(j=I; j<cols_; j++) (*phmx).put_h(get(J,j),I,j);	// Set new row I from J
  for(j=J; j<cols_; j++) (*phmx).put_h(get(I,j),J,j);	// Set new row J from I
  for(i=0; i<=I; i++)    (*phmx).put_h(get(i,J),i,I);	// Set new col I from J
  for(i=0; i<=J; i++)    (*phmx).put_h(get(i,I),i,J);	// Set new col J from I
  (*phmx).put_h(get(J,J), I, I);			// New <i|hmx|i> was <j|hmx|j>
  (*phmx).put_h(get(I,I), J, J);			// New <j|hmx|j> was <i|hmx|i>
  (*phmx).put_h(get(J,I), I, J);			// New <j|hmx|j> was <i|hmx|i>
  return phmx;
  }

void h_matrix::permute_ip(int I, int J)
  {
//  complex zIJ = get(I,J);			// <I|hmx|J>
  complex zJI = get(J,I);			// <J|hmx|I>
  complex zII = get(I,I);			// <I|hmx|I>
  complex zJJ = get(I,I);			// <J|hmx|J>
  int i,j;
  complex *rowI, *rowJ;				// For matrix rows
  rowI = new complex[cols_];			// Allocate row
  rowJ = new complex[cols_];			// Allocate row
  for(j=0; j<cols_; j++)			// Fill rows
    {
    rowI[j] = get(I,j);				// <I|hmx|j>
    rowJ[j] = get(J,j);				// <J|hmx|j>
    }
  complex *colI, *colJ;				// For matrix cols
  colI = new complex[rows_];			// Allocate col
  colJ = new complex[rows_];			// Allocate col
  for(i=0; i<cols_; i++)			// Fill cols
    {
    colI[i] = get(i, I);			// <i|hmx|I>
    colJ[i] = get(i, J);			// <i|hmx|J>
    }
  for(j=I; j<cols_; j++) put_h(rowJ[j],I,j);	// Set new row I from J
  for(j=J; j<cols_; j++) put_h(rowI[j],J,j);	// Set new row J from I
  for(i=0; i<=I; i++)    put_h(colJ[i],i,I);	// Set new col I from J
  for(i=0; i<=J; i++)    put_h(colI[i],i,J);	// Set new col J from I
  put(zJJ, I, I);				// New <i|hmx|i> was <j|hmx|j>
  put(zII, J, J);				// New <j|hmx|j> was <i|hmx|i>
  put_h(zJI, I, J);				// New <j|hmx|j> was <i|hmx|i>
  delete [] rowI;
  delete [] rowJ;
  delete [] colI;
  delete [] colJ;
  return;
  }
    

/*
  NPtr nmx  = NMX(); 			// Copy hmx to nmx 
  NPtr pnmx = (NPtr)nmx->permute(I,J);	// Apply the permutaiton
  delete nmx;				// Don't need this now
  HPtr phmx = pnmx->HMX();		// Reset it to hmx
  delete pnmx;				// Don't need this now
  return phmx;
// NMX Row Swap
  n_matrix* nsr = new n_matrix(*this);          // Make new n_matrix
  int I0 = I*cols_;                             // Index of <I|mx,mx'|0>
  int J0 = J*cols_;                             // Index of <J|mx,mx'|0>
  for(int k=0; k<cols_; k++)                    // Swap rows I & J
    {
    nsr->data[I0+k] = data[J0+k];               // Set <I|mx'|k> = <J|mx|k>
    nsr->data[J0+k] = data[I0+k];               // Set <J|mx'|k> = <I|mx|k>
    }
*/

// sosiz - the above routine is a kludge, but it worked OK. A permutation 
//         on hmx will produce hmx but the code to do so is more difficult
//          within h_matrix.  So I just spawned it to normal matrix then
//          cast it back to hermitian.  This uses extra arrays ......
//         I've since fixed the routines to just stay Hermitian and the
//         above can be deleted when I'm convinced they work fine.

double h_matrix::maxRe() const
  { 
  double maxval=-HUGE_VAL;
  for(int i=0; i<size; i++) { maxval = gmax(Re(data[i]), maxval); }
  return (size)?maxval:0;  
  }

double h_matrix::maxIm() const
  { 
  double maxval=-HUGE_VAL;
  for(int i=0; i<size; i++)
    {
    maxval = gmax(Im(data[i]), maxval);
    maxval = gmax(-Im(data[i]), maxval);
    }
  return (size)?maxval:0;  
  }

complex h_matrix::maxZ()  const
  { 
  double maxval=-HUGE_VAL;
  complex maxz;
  for(int i=0; i<size; i++)
    {
    if(norm(data[i]) > maxval)
      {
      maxval = norm(data[i]);
      maxz = data[i];
      }
    }
  return (size)?maxz:complex0;
  }

double h_matrix::minRe() const
  { 
  double minval=HUGE_VAL;
  for(int i=0; i<size; i++) minval = gmin(Re(data[i]), minval);
  return (size)?minval:0;
  }
 
double h_matrix::minIm() const
  {
  double minval=HUGE_VAL;
  for(int i=0; i<size; i++)
    {
    minval = gmin(Im(data[i]), minval);
    minval = gmin(-Im(data[i]), minval);
    }
  return (size)?minval:0;
  }
 
complex h_matrix::minZ()  const
  {
  double minval=HUGE_VAL;
  complex minz;
  for(int i=0; i<size; i++)
    {
    if(norm(data[i]) < minval)
      {
      minval = norm(data[i]);
      minz = data[i];
      }
    }
  return (size)?minz:complex0;
  }

// ____________________________________________________________________________
// I                  CLASS H_MATRIX SIMPLE BINARY FUNCTIONS
// ____________________________________________________________________________



	// Input            hmx : A h_matrix (this)
	//		    mx  : A second matrix
	// Output           z   : The trace of hmx*mx
        // Note                 : This is faster than taking the matrix
        //                        product and then using unary trace!
        // Note                 : Uses internal structure of other matrix types
 
/*                ---                 --- ---
              z = \   <i|hmx*mx|i>  = \   \   <i|hmx|j><j|mx|i>
                  /                   /   /
                  ---                 --- ---
                   i                   i   j                                 */

complex h_matrix::trace(_matrix* mx)
  {
  if((rows_!=mx->cols()) || cols_!=mx->rows())		// Insure dimensions proper
    {
    HMxerror(17, 1);					// Dimension problems
    HMxerror(6, "trace(mx1 * mx2)", 1);			// Bad use of function 
    HMxfatal(3, "trace");				// Array dimensions mismatched
    }
  else
    {
    complex z(0);					// Initialize the result at zero
    switch(mx->stored_type())
      {
      case n_matrix_type: 				// Tr{h_matrix * n_matrix}    
	{
	complex *n00 = ((n_matrix*)mx)->data;		// Start of nmx: n00 = <0|nmx|0>
	complex *h00 = data;				// Start of hmx: h00 = <0|hmx|0>
        complex *nend = n00 + rows_*cols_;		// End of nmx: <rows_|nmx|cols_+1>
        complex *nji,*nii,*hij,*hji;			// These are the matrix elements
        int hrow;					// <j|hmx|i>+hrow -> <j+1|hmx|i>
        complex *h0i = h00;				// <0|hmx|i> -> <0|hmx|0>
        complex *n0i = n00;				// <j|nmx|0> -> <0|nmx|0>
	for(nii=n00; nii<nend; n0i++,h0i++,nii+=cols_+1)// Effective loop over i (via nii)
          {
          hji = h0i;					// <j|hmx|i> -> <0|hmx|i>
          nji = n0i;					// <j|nmx|i> -> <0|nmx|i>
	  hrow = cols_;					// <j|hmx|i>+hrow -> <j+1|hmx|i>
	  for( ; nji<nii; nji+=cols_,hrow--,hji+=hrow) 	// Begin effective loop over j, i>j
	    z += conj(*hji) * (*nji);			// Tr{hmx*mx} += <j|hmx*|i><j|mx|i>, i>j
	  for(hij=hji; nji<nend; hij++,nji+=cols_)	// Continue effective loop over j, i<=j
	    z += (*hij) * (*nji);			 // Tr{hmx*mx} += <i|hmx*|j><j|mx|i>, i<=j
          }
	return z;
	}
	break;
      case d_matrix_type:	 			// Tr{h_matrix * d_matrix}    
	{
        complex *h00 = data;				// Start of hmx: h00 = <0|hmx|0>
	complex *d00 = ((d_matrix*)mx)->data;		// Start of dmx: d00 = <0|dmx|0>
        complex *dend = d00 + cols_;			// End of dmx: <rows_|dmx|cols+1>
        complex *hii, *dii;				// These are the matrix elements
        int hrow = cols_;				// Row increment: <i|hmx|i> + hrow
	for(hii=h00,dii=d00; dii<dend;			//                  -> <i+1|hmx|i+1>
                                hii+=hrow,dii++,hrow--)
	  z += (*hii) * (*dii);				// Tr{hmx*dmx} += <i|hmx|i><i|dmx|i>
	return z;
	}
	break;
      case i_matrix_type: 				// Tr{h_matrix * i_matrix}    
        { return trace(); }				// Same as Tr{hmx}
        break;
      case h_matrix_type:                 
	{
        complex *h00 = data;				// Start of hmx: h00 = <0|hmx|0>
	complex *m00 = ((h_matrix*)mx)->data;		// Start of mx: m00 = <0|dmx|0>
        complex *hij,*hji, *mij, *mji;			// These are the matrix elements
        complex *h0i = h00;				// <0|hmx|i> -> <0|hmx|0>
        complex *m0i = m00;				// <0|mx|i> -> <0|mx|0>
        complex *mii = m00;				// <i|mx|i> -> <0|mx|0>
        complex *m11 = m00+cols_;			// End of 1st mx row <1|mx|1>
        int rinc;
	for(; m0i<m11; h0i++,m0i++)			// Effective loop over i (m0i)
          {
          hji = h0i;					// <j|hmx|i> = <0|hmx|i>
          mji = m0i;					// <j|mx|i> = <0|mx|i>
          rinc = cols_;					// Matrix row incrementation
	  for(; mji<mii; rinc--,hji+=rinc,mji+=rinc)	// Loop over j, j<i
	    z += conj(*hji) * (*mji);			// Tr{hmx*mx} += <j|hmx*|i><j|mx|i>
          mii += rinc;					// <i|mx|i> -> <i+1|mx|i+1>
	  for(hij=hji,mij=mji; mij<mii; hij++,mij++)	// Continue loop over j, j>=i
	    z += (*hij) * conj(*mij);			// Tr{hmx*mx} += <i|hmx|j><i|mx*|j>
          }
	return z;
	}
        break;
      default:                 
	{
        complex *h00 = data;				// Start of hmx: h00 = <0|hmx|0>
        complex *hij,*hji,*h0i;				// These are the matrix elements
        int i,j;
	for(i=0,h0i=h00; i<rows_; i++,h0i++)		// Loop over the index i
          {
	  for(j=0,hji=h0i; j<i; j++,hji+=cols_-j)	// Begin loop over index j, j<i
	    z += conj(*hji) * (*mx)(j,i);		// Tr{hmx*mx} += <j|hmx*|i><j|mx|i>, j<i
	  for(hij=hji; j<cols_; j++, hij++)		// Continue loop over index j, j>=i
	    z += (*hij) * (*mx)(j,i);			// Tr{hmx*mx} += <i|hmx|j><j|mx|i>, j>=i
          }
	return z;
	}
        break;
      }
    }
  return complex0;
  }



	// Input            hmx : A h_matrix (this)
	//		    mx  : A second matrix
	// Output           mx1 : The product of adjoint(hmx)*mx
        // Note                 : This is not much faster than taking the matrix
        //                        adjoint and then the product.
        // Note                 : Since, by definition, hmx is self adjoint,
	//			  this is equivalent to hmx*mx
        // Note                 : Uses internal structure of other matrix types
 
//                  ---                       ---
//     <i|pdt|j>  = \   <k|hmx*|i><k|mx|j> =  \   <i|hmx|k><k|mx|j>
//                  / 			      /
//                  ---                       ---
//                   k                         k

_matrix* h_matrix::adjoint_times(_matrix* mx)
  {
  if(rows_ != mx->rows())			// Insure dimensions are proper
    {
    HMxerror(17, 1);				// Dimension problems
    HMxerror(6, "adjoint(mx1) * mx2", 1);	// Bad use of function 
    HMxfatal(3, "adjoint_times");		// Array dimensions mismatched
    }
  else
    {
    switch(mx->stored_type())
      {
      case i_matrix_type:			// Adjoint[hmx]*imx = hmx
        { return this; }			// No change, just return hmx
        break;
      case d_matrix_type:			// Adjoint[hmx]*dmx = hmx*dmx
        { return multiply((d_matrix*)mx); }	// Return hmx*dmx
        break;
      case h_matrix_type:			// Adjoint[hmx]*hmx1 = hmx*hmx1
        { return multiply((h_matrix*)mx); }	// Return hmx*hmx
        break;
      case n_matrix_type:			// Adjoint[hmx]*nmx = hmx*nmx
        { return multiply((n_matrix*)mx); }	// Return hmx*nmx
        break;
      default:					// Adjoint[hmx]*unknown mx = hmx*mx
        break; 					// Return hmx*mx
      }
    }
  return multiply(mx);
  }



	// Input            hmx : A h_matrix (this)
	//		    mx  : A second matrix
	// Output           mx1 : The product of hmx * adjoint(mx)
        // Note                 : This is faster than taking the matrix
        //                        adjoint and then the product.
 
/*                                  ---
                       <i|pdt|j>  = \   <i|hmx|k><j|mx*|k>
                                    /
                                    ---
                                     k                                       */

_matrix* h_matrix::times_adjoint(_matrix* mx)
  {
  if(cols_ != mx->cols())				// Insure matrix dimensions proper
    {
    HMxerror(17, 1);					// Dimension problems
    HMxerror(6, "mx1*adjoint(mx2)", 1);			// Bad use of function 
    HMxfatal(3, "times_adjoint");			// Array dimensions mismatched
    return mx;
    }
  else
    switch (mx->stored_type())
      {
      case n_matrix_type: 				// hmx * adjoint(nmx)
	{ 
#ifdef _USING_BLAS_
        int A_rows = rows();
        int A_cols = cols();            
        int B_cols = mx->cols();
        int C_rows = A_rows;
        int C_cols = B_cols;
        n_matrix* pdt =	new n_matrix(C_rows, C_cols);
        if(C_rows * C_cols > 16)  //4*4 = 16
        {
          double alpha[2] = {0,0};
          double beta[2] = {0,0};
          alpha[0] = 1.0;
          beta[0]  = 0.0;
          n_matrix* hmxA =	new n_matrix(A_rows,A_cols);		// Create new matrix h_matrix
          this->convert(hmxA);				// Convert h_matrix mx into normal matrix hmx
	  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, A_rows, B_cols, A_cols, 
                           alpha, (double *)hmxA->data, A_cols, (double *)((n_matrix*)mx)->data, B_cols, beta, (double *)pdt->data, C_cols);
//        std::cerr << "BLAS: h_matrix * n_matrix\n";
	  delete hmxA;
        }
        else
        {
	int c = mx->rows();				// Columns of product matrix
	complex *p00 = pdt->data;			// Start of data in pdt: <0|pdt|0>
	complex *h00 = data;				// Start of data in hmx: <0|hmx|0>
	complex *n00 = ((n_matrix*)mx)->data;		// Start of data in nmx: <0|nmx|0>
	complex *pij,*hik,*hki,*h0i,*njk,*nj0,*nji;	// These are the matrix elements
	complex z;					// intermediate storage
        int i,j,hrow;
	for(pij=p00,h0i=h00,i=0; i<rows_; i++,h0i++)	// Effective loop over i
	  for(j=0,nj0=n00; j<c; j++,pij++)		// Effective loop over j
            {
            z = 0;					// Initialize the element at zero
            hrow = cols_;				// Initialize hmx row length
            hki = h0i;					// <k|hmx|i> = <0|hmx|i>
            njk = nj0;					// <j|nmx|k> = <j|nmx|0>
            nji = nj0+i;				// <j|nmx|i> = <j|nmx|i>
	    for(; njk<nji; njk++,hrow--,hki+=hrow)	// Begin loop over k, k<i (via njk)
	      z += conj(*hki) * conj(*njk);		// <i|pdt|j> += <k|hmx*|i><j|nmx*|k>, k<i
            nj0 += rows_;				// <j|nmx|0> -> <j+1|nmx|0>
	    for(hik=hki; njk<nj0; hik++,njk++)		// Continue loop over k, k>=i (via njk)
	      z += conj(*njk, *hik);			// <i|pdt|j> += <i|hmx|k><j|nmx*|k>, k>=i
	    *pij = z;
            }
         }
#else
	int c = mx->rows();				// Columns of product matrix
	n_matrix* pdt = new n_matrix(rows_,c);		// Create  a new n_matrix for result
	complex *p00 = pdt->data;			// Start of data in pdt: <0|pdt|0>
	complex *h00 = data;				// Start of data in hmx: <0|hmx|0>
	complex *n00 = ((n_matrix*)mx)->data;		// Start of data in nmx: <0|nmx|0>
	complex *pij,*hik,*hki,*h0i,*njk,*nj0,*nji;	// These are the matrix elements
	complex z;					// intermediate storage
        int i,j,hrow;
	for(pij=p00,h0i=h00,i=0; i<rows_; i++,h0i++)	// Effective loop over i
	  for(j=0,nj0=n00; j<c; j++,pij++)		// Effective loop over j
            {
            z = 0;					// Initialize the element at zero
            hrow = cols_;				// Initialize hmx row length
            hki = h0i;					// <k|hmx|i> = <0|hmx|i>
            njk = nj0;					// <j|nmx|k> = <j|nmx|0>
            nji = nj0+i;				// <j|nmx|i> = <j|nmx|i>
	    for(; njk<nji; njk++,hrow--,hki+=hrow)	// Begin loop over k, k<i (via njk)
	      z += conj(*hki) * conj(*njk);		// <i|pdt|j> += <k|hmx*|i><j|nmx*|k>, k<i
            nj0 += rows_;				// <j|nmx|0> -> <j+1|nmx|0>
	    for(hik=hki; njk<nj0; hik++,njk++)		// Continue loop over k, k>=i (via njk)
	      z += conj(*njk, *hik);			// <i|pdt|j> += <i|hmx|k><j|nmx*|k>, k>=i
	    *pij = z;
            }
#endif
	return pdt;
	}
	break;
      case d_matrix_type:				// hmx * adjoint(dmx)
	{
	n_matrix* pdt =	new n_matrix(cols_,cols_); 	// Construct a new normal matrix
	complex *p00 = pdt->data; 			// Start of pdt: p00 = <0|pdt|0>
	complex *h00 = data;				// Start of hmx: h00 = <0|hmx|0>
	complex *d00 = ((d_matrix*)mx)->data;		// Start of dmx: d00 = <0|dmx|0>
        complex *dend = d00 + cols_;			// End of dmx: dend = <cols_|dmx|cols_+1>
	complex *pij, *h0i, *hij, *hji, *djj, *dii;	// These will be actual matrix elements
        int hrow;					// hji -> hji+hrow: <j|hmx|i> -> <j+1|hmx|i>
	for(dii=d00,h0i=h00,pij=p00; dii<dend;		// Effective loop over i
                                           h0i++,dii++)
          {
	  for(djj=d00,hji=h0i,hrow=cols_; djj<dii;	// Effective loop over j
			  pij++,djj++,hrow--,hji+=hrow)
            (*pij) = conj(*hji) * conj(*djj);		// <i|pdt|j> = <j|hmx*|i><j|dmx*|j>, i>j
	  for(hij=hji; djj<dend; pij++,hij++,djj++)
            (*pij) = conj(*djj, *hij);			// <i|pdt|j> = <i|hmx|j><j|dmx*|j>, i<=j
	  }
	return pdt;
	}
	break;
      case h_matrix_type:				// hmx * adjoint(hmx1)
	return multiply((h_matrix*)mx);			// Same as hmx*hmx1
	break;
      case i_matrix_type:				// hmx * adjoint(imx)
	return this;					// Just return hmx
	break;
      default:						// hmx * adjoint(generic mx)
	{ 
#ifdef _USING_BLAS_
        int A_rows = rows();
        int A_cols = cols();            
        int B_cols = mx->cols();
        int C_rows = A_rows;
        int C_cols = B_cols;
        n_matrix* pdt =	new n_matrix(C_rows, C_cols);
        if(C_rows * C_cols > 16)  //4*4 = 16
        {
          double alpha[2] = {0,0};
          double beta[2] = {0,0};
          alpha[0] = 1.0;
          beta[0]  = 0.0;
          n_matrix* hmxA =	new n_matrix(A_rows,A_cols);		// Create new matrix h_matrix
          this->convert(hmxA);				// Convert h_matrix mx into normal matrix hmx
	  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, A_rows, B_cols, A_cols, 
                           alpha, (double *)hmxA->data, A_cols, (double *)((n_matrix*)mx)->data, B_cols, beta, (double *)pdt->data, C_cols);
//        std::cerr << "BLAS: h_matrix * unknown\n";
	  delete hmxA;
        }
        else
        {
	int c = mx->rows();				// Columns of matrix product
	complex *p00 = pdt->data; 			// Start of pdt: p00 = <0|pdt|0>
	complex *h00 = data;				// Start of hmx: h00 = <0|hmx|0>
	complex *pij, *hik, *hki, *h0i;			// These will be actual matrix elements
        int i, j, k;
	for(pij=p00,h0i=h00,i=0; i<rows_; i++,h0i++)
	  for(j=0; j<c; j++,pij++)
            {
	    for(k=0,hki=h0i; k<i; k++, hki+=cols_-k)	// Loop over internal index, i>k
	      *pij += conj((*mx)(j,k))*conj(*hki);	// <i|ptd|j> = <k|hmx*|i><j|mx*|k>, i>k
	    for(hik=hki; k<cols_; k++,hik++)		// Continut loop over internal index, i<=k
	      *pij += conj((*mx)(j,k), *hik);	 	// <i|pdt|j> = <i|hmx|k><j|mx*|k>, i<=k
            }
	}
#else
	int c = mx->rows();				// Columns of matrix product
	n_matrix* pdt = new n_matrix(rows_,c,complex0);	// Construct new n_matrix for result
	complex *p00 = pdt->data; 			// Start of pdt: p00 = <0|pdt|0>
	complex *h00 = data;				// Start of hmx: h00 = <0|hmx|0>
	complex *pij, *hik, *hki, *h0i;			// These will be actual matrix elements
        int i, j, k;
	for(pij=p00,h0i=h00,i=0; i<rows_; i++,h0i++)
	  for(j=0; j<c; j++,pij++)
            {
	    for(k=0,hki=h0i; k<i; k++, hki+=cols_-k)	// Loop over internal index, i>k
	      *pij += conj((*mx)(j,k))*conj(*hki);	// <i|ptd|j> = <k|hmx*|i><j|mx*|k>, i>k
	    for(hik=hki; k<cols_; k++,hik++)		// Continut loop over internal index, i<=k
	      *pij += conj((*mx)(j,k), *hik);	 	// <i|pdt|j> = <i|hmx|k><j|mx*|k>, i<=k
            }
#endif
	return pdt;
	}
      }
  }

 
// ____________________________________________________________________________
// J                 CLASS H_MATRIX COMPLEX UNARY FUNCTIONS
// ____________________________________________________________________________


complex h_matrix::det()
  {
  HMxerror(24, 1);				// Function not implemented
  HMxfatal(7, "Determinant");			// Function not yet implemented
  return complex0;
  }

int h_matrix::rank()
  {
  HMxerror(24, 1);				// Function not implemented
  HMxfatal(7, "Rank");				// Function not yet implemented
  return 0;
  }

 
// ____________________________________________________________________________
//                      CLASS H_MATRIX COMPLEX FUNCTIONS
// ____________________________________________________________________________

// ************************** tensor product ****************************

_matrix* h_matrix::tensor_product(_matrix* mx) 

        // Input            hmx : An h_matrix (this)
        //                   mx : A second matrix
        // Output           pdt : A matrix which is the tensor product
        //                        of the two input matrices

        //                           pdt        =   hmx (x) mx

        //                       (m*o x m*p)       (mxm)   (oxp)

        //                    <i*o+k|pdt|j*p+l> = <i|hmx|j><k|mx|l>

  {
  switch (mx->stored_type())
    {
    case n_matrix_type:				// Tensor Product: h_matrix (x) n_matrix
      {
      int hd = cols_;				// Dimension of hmx
      int cn = mx->cols();			// Columns of mx
      int rn = mx->rows();			// Rows of mx
      int rp = hd*rn;				// Columns in tensor product
      int cp = hd*cn;				// Columns in tensor product
      n_matrix* pdt = new n_matrix(rp,cp);	// New n_matix, pdt, much larger
      complex *pikjl = pdt->data;		// Start of data in pdt: <ik|pdt|jl>-><00|pdt|00>
      complex *h00 = data;			// Start of data in hmx: <0|hmx|0>
      complex *n00 = ((n_matrix*)mx)->data; 	// Start of data in nmx: <0|nmx|0>
      complex *nend = n00 + rn*cn;		// End of data in nmx: <0|nmx|0>
      complex *nkl, *nk0, *nkp10;		// Elements <k|nmx|l>, <k|nmx|0>, <k+1|nmx|0>
      complex *hij, *hji, *h0i;			// elements <i|hmx|j>, <j|hmx|i>, <0|hmx|i>
      complex *hii, *hxii;			// Elements <i|hmx|i> and <i|hmx|i> (temp)
      complex *h11 = h00 + hd;			// <1|hmx|1>
      int hrowi=hd;				// Used to track length of hmx rows, index i
      int hrowj;				// Used to track length of hmx rows, index j
      for(h0i=h00,hii=h00; h0i<h11;
                      h0i++,hii+=hrowi,hrowi--)
        for(nk0=n00; nk0<nend; nk0+=cn)		// Effective loop over k (nmx rows, via nk0)
          {
          hxii = hii;				// <i|hmx|i> (temp) -> <i|hmx|i> (stable)
          nkp10 = nk0 + cn;			// <k+1|nmx|0> = <k|nmx|0> + cn
          hrowj = hd;				// Length of nmx row, start row 0
          for(hji=h0i;hji<hxii;			// Effective loop over j (hmx columns, via hji), j<i
                             hrowj--,hji+=hrowj)
            for(nkl=nk0;nkl<nkp10;pikjl++,nkl++)// Effective loop over l (nmx columns, via nkl)
              *pikjl = conj(*hji) * (*nkl);	// <ik|pdt|jl> = <j|hmx*|i><k|nmx|l>, j<i
          hxii = hji + hrowj;			// <i|hmx|i> -> <i+1|hmx|i+1>
          for(hij=hji; hij<hxii; hij++)		// Continue effective looop over j, j>=i
            for(nkl=nk0;nkl<nkp10;pikjl++,nkl++)// Effective loop over l (nmx columns, via nkl)
              *pikjl = (*hij) * (*nkl);		// <ik|pdt|jl> = <i|hmx|j><k|nmx|l>, j>=i
          }
      return pdt;
      }
      break;
    case h_matrix_type:				// Tensor Product: h_matrix (x) h_matrix
      {
      int hd = cols_;				// Dimension of hmx
      int md = mx->cols();			// Dimension of mx
      int rd = hd*md;				// Dimenson of tensor product
      n_matrix* pdt = new n_matrix(rd,rd);	// New n_matix, pdt, much larger
      complex *p0000 = pdt->data;		// Start of data in pdt: <00|pdt|00>
      complex *h00 = data;			// Start of data in hmx: <0|hmx|0>
      complex *m00 = ((h_matrix*)mx)->data; 	// Start of data in mx: <0|mx|0>
      complex *pikjl;
      complex *mkl, *mlk, *m0k;
      complex *hij, *hji, *h0i;
      int i=0,j,k,l;
      for(pikjl=p0000,h0i=h00; i<hd; i++,h0i++)	// Loop over index i (rows of hmx)
        for(k=0,m0k=m00; k<md; k++,m0k++)	// Loop over the rows of nmx
          {
          for(j=0,hji=h0i; j<i; j++,hji+=hd-j)	// Loop over the columns of hmx, j<i
            {
            for(l=0,mlk=m0k; l<k; l++,pikjl++,mlk+=md-l)// Loop over the columns of nmx, l<k
              *pikjl = conj(*hji) * conj(*mlk);
            for(mkl=mlk; l<md; l++,pikjl++,mkl++)	// Loop over the columns of nmx
              *pikjl = conj(*hji) * (*mkl);
            }
          for(hij=hji; j<hd; j++,hij++)			// Loop over the columns of hmx, j>=i
            {
            for(l=0,mlk=m0k; l<k; l++,pikjl++,mlk+=md-l)// Loop over the columns of nmx, l<k
              *pikjl = (*hij) * conj(*mlk);
            for(mkl=mlk; l<md; l++,pikjl++,mkl++)	// Loop over the columns of nmx
              *pikjl = (*hij) * (*mkl);
            }
          }
// sosiz - This needs to return Hermitian... yes it works out that way.
//         But the math above would be a big pain so I am just doing it
//         cheap and dirty below.  I did work the above out once.. where is it?
      h_matrix* hpdt = pdt->HMX();
      delete pdt;
      return hpdt;
//      return pdt;
      }
      break;
    case d_matrix_type:				// Tensor Product h_matrix (x) d_matrix
      {
      int hd = cols_;				// Dimension of hmx
      int dd = mx->cols();			// Dimension of dmx
      int pd = hd*dd;				// Dimension of product
      complex *h00 = data;			// Start of data in hmx: <0|hmx|0>
      complex *d00 = ((d_matrix*)mx)->data; 	// Start of data in dmx: <0|dmx|0>
      complex *hend = h00 + (hd*hd+hd)/2;       // End of data in hmx: <hd|hmx|hd+1>
      complex *dend = d00 + dd;			// End of data in dmx: <dd|dmx|dd+1>
      complex *dkk;				// Element <k|dmx|k>
      if(mx->is_real())				// If dmx is real, the tensor product
        {					// will be Hermitian
        h_matrix* pdt = 			// New h_matrix for tensor product
                    new h_matrix(pd,pd,complex0);
        complex *pi0i0 = pdt->data;		// Start of data in pdt: <i0|pdt|i0>-><00|pdt|00>
        complex *pikjk, *pikik;			// Elements <ik|pdt|jk> & <ik|pdt|ik> 
        complex *hii, *hij, *zii;		// Elements <i|hmx|i>, <i|hmx|j>, & <i|hmx|i> (temp)
        int prow = pd;				// Length pdt row i, set at i=0
        int hrow = hd;				// Length hmx row i, set at i=0
        for(hii=h00; hii<hend;			// Effective loop over i (hmx rows, via hii)
                  hii+=hrow,pi0i0=pikik,hrow--)	// pi0i0 = <i0|pdt|i0> & <i|hmx|i>+hrow-><i+1|hmx|i+1>
          {
          for(pikik=pi0i0,dkk=d00; dkk<dend;	// Effective loop over k (dmx rows, via dkk)
                      dkk++,pikik+=prow,prow--)	// Loop over index k (dmx rows)
            {
            hij = hii;				// <i|hmx|j> = <i|hmx|i>
            zii = hii + hrow;			// <i|hmx|i> -> <i+1|hmx|i+1>
            for(pikjk=pikik; hij<zii;		// Effective loop over j (hmx columns), j>=i
                               hij++,pikjk+=dd)
              *pikjk = (*hij) * (*dkk); 	// <i*dd+k|pdt|j*dd+l> = del   <i|hmx|j><k|dmx|k>
            }					//                          k,l
          }
        return pdt;
        }
      else					// If dmx is complex, the tensor product
        {					// cannot be Hermitian, it must be n_matrix
        n_matrix* pdt =				// New n_matrix for tensor product
                    new n_matrix(pd,pd,complex0);
        complex *pi000 = pdt->data;		// Start of data in pdt: <i0|pdt|00> -> <00|pdt|00>
        complex *pend = pi000 + pd*pd;		// End of data in pdt: <pd|pdt|pd+1>
        complex *pikjk, *pik0k;			// Elements <ik|pdt|jk> & <ik|pdt|0k> 
        complex *hij, *hji, *h0i;		// Elements <i|hmx|j>, <j|hmx|i>, & <0|hmx|i>
        int i=0,j;
        for(h0i=h00; pi000<pend;		// Effective loop over i (hmx rows, via h0i)
                       i++,pi000+=dd*pd,h0i++)
          for(dkk=d00,pik0k=pi000; dkk<dend;	// Effective loop over k (dmx rows, via dkk)
                            pik0k+=pd+1,dkk++)
            {
            for(j=0,hji=h0i,pikjk=pik0k; j<i;	// Loop over index j (hmx columns), j<i
                      j++,hji+=hd-j,pikjk+=dd)
              *pikjk = conj(*hji)*(*dkk);
            for(hij=hji; j<hd;			// Continue loop over j (hmx columns), j>=i
                          j++,hij++,pikjk+=dd)
              *pikjk = (*hij) * (*dkk);		// <ik|pdt|jk> = <i*dd+k|pdt|j*dd+k> = <i|hmx|j><k|dmx|k>
            }
        return pdt;
        }
      }
      break;
    case i_matrix_type:				// Tensor Product h_matrix (x) i_matrix
      {
      int hd = cols_;				// Dimension of hmx
      int id = mx->cols();			// Dimension of imx
      int pd = hd*id;				// Dimension of product
      h_matrix* pdt=new h_matrix(pd,pd,complex0);	// Construct a new Hermitian array for product
      complex *p0000 = pdt->data;		// Start of data in pdt: <00|pdt|00>
      complex *h00 = data;			// Start of data in hmx: <0|hmx|0>
      complex *hend = h00 + (hd*hd+hd)/2;       // End of data in hmx: <hd|hmx|hd+1>
      complex *pikjk, *pikik;			// Elements <ik|pdt|jk> & <ik|pdt|ik> 
      complex *hij;
      complex *zii, *hii;
      int prow = pd;				// Length pdt row i, set at i=0
      int hrow = hd;				// Length hmx row i, set at i=0
      int k;
      complex *pi0i0 = p0000;

      for(hii=h00; hii<hend;			// Effective loop over i (hmx rows, via hii)
              hii+=hrow,pi0i0=pikik,hrow--)	// pi0i0 = <i0|pdt|i0> & <i|hmx|i>+hrow-><i+1|hmx|i+1>
        {
        k = 0;
        for(pikik=pi0i0; k<id; k++,pikik+=prow,prow--)	// Loop over index k (dmx rows)
          {
          hij = hii;
          zii = hii + hrow;				// <i|hmx|i> -> <i+1|hmx|i+1>
          for(pikjk=pikik; hij<zii; hij++,pikjk+=id)	// Effective loop over j (hmx columns), j>=i
            *pikjk = *hij; 				// <i*id+k|pdt|j*id+l> = del   * <i|hmx|j>
          }						//                          k,l
        }
      return pdt;
      }
      break;      
    default:					// Tensor Product h_matrix (x) generic matrix 
      { 
      int hd = cols_;				// Dimension of hmx
      int rg = mx->rows();			// Rows in gmx
      int cg = mx->cols();			// Columns in gmx
      n_matrix* pdt = new n_matrix(hd*rg,hd*cg);
      complex *pikjl = pdt->data;		// Start of data in pdt: <00|pdt|00>
      int i,j,k,l;
      for(i=0; i<hd; i++)
        for(k=0; k<rg; k++)
          for(j=0; j<hd; j++)
            for(l=0; l<cg; l++,pikjl++)
              *pikjl = (*this)(i,j)*(*mx)(k,l);
      return pdt;
      }
    }
  
  }

// ____________________________________________________________________________
// J                 CLASS H_MATRIX COMPLEX UNARY FUNCTIONS
// ____________________________________________________________________________

// Note:  The alogorithm for FFT is not found here.  The Hermitian matrix is
//        converted to a normal matrix and FFT happens in class n_matrix.




// ____________________________________________________________________________
// N                CLASS H_MATRIX DIAGONALIZATION FUNCTIONS
// ____________________________________________________________________________

/*************************************************************************
**                                                                      **
** The routine cred reduces a complex Hermitian matrix to a complex     **
** Hermitian tridiagonal form using the Householder algorithm. During	**
** the course of the transformation, the matrix becomes non-Hermitian.  **
** Thus, we initially convert the array to a normal matrix, perform the **
** tranformation, then convert back.					**
**                                                                      **
*************************************************************************/

// sosi - It would be nice NOT to have to convert hmx->nmx->hmx but this
//        would require a more complicated storage scheme during the
//        left transformation step, when the array is temporarily not
//        Hermitian. It would also be nice to keep U as hmx rather than
//        than nmx. Also, should change the data type of U to reflect that
//        it is Hermitian? But then that would mess up the next step...

void h_matrix::cred(_matrix* &U)
 
// Input        hmx     : A Hermitian matrix (this)
//              U       : Matrix (pointer) to contain eigenvectors
// Output       void    : Both matrices, hmx & U are altered
//                        Array hmx will be set tri-diagonal
//			              Array U will be Hermitian & Unitary
// Note			: On input U should be NULL.  It is
//			      set to a new array herein
 
{

  //	         Set Up A Complex Vector For Transformations

  complex *W;					// Transformation vector
  W = new complex[rows_];
  complex *wend = W + rows_;			// End of vector W + 1
  complex *wl=W, *wlp1, *wi, *wj;		// To access W elements

  //		       Set Up An Array Of Eigenvectors

  U = new n_matrix(rows_,rows_,complex0);	// Alocate eigenvector mx
  complex *U00 = ((n_matrix*)U)->data;		// U00 = <0|U|0>
  complex *Uend = U00 + rows_*rows_;		// End of array U
  complex *Ui0, *Uij, *Uil;			// To access U elements
  complex *Uii = U00;				// <i|U|i> <- <0|U|0>
  for(; Uii<Uend; Uii+=rows_+1)			// Initialize U to I matrix
    *Uii = complex1;

  //	            Set Up A Working Copy Of Input Array
  //		(Wish We Could Keep This Hermitian & Avoid It!)

  n_matrix* nmx = NMX();			// Copy hmx to nmx for now
  complex *n00  = nmx->data;			// n00  = <0|nmx|0>
  complex *nend = n00+rows_*rows_;		// nend = <rows_|nmx|rows_>+1
  complex *ni0, *nij, *nji, *nnil, *nli;		// To access nmx elements
  complex *nl0, *nlend, *nll, *nlp1l; 		// To access nmx elements

  //                      Householder Reduction Section

  /* We Loop Through Columns Of The Input Hermitian Matrix. First, All Elements
     That Are Not Tridiagonal (in the lth Row) Are Checked To See If They Are
     Non-Zero. If So, Then We Do Some Math For The Transformation. Note That On
     Successive Columns We Use The Updated Array With Alterations From The 
     Previous Column Adjustment In Place.                                      */

  double s, sw;					// Vector norms
  double csf;					// Column scaling factor
  for(int l=0; l<rows_-2; l++, wl++)		// Loop over l
  {

    //	Check For Non-Zero Tridiagonal Elements Lower Part of Column l
 
    csf = 0.0;					// Column scaling factor
    nll = n00 + l*rows_ + l; 			// Diagonal this l: <l|nmx|l>
    nnil = nll + rows_ + rows_;			// Start at <l+2|nmx|l> 
    for(; nnil<nend; nnil+=rows_)			// Get scaling factor from col
      csf += AbsNorm(*nnil);			// elements not tri-diagonal

    /*            Next Generate The Vector W For Column l If Needed

    This section assumes that  1. nll=<l|mx|l>   2. nend=<rows_|mx|rows_>+1
                               3. wl = <l|W>                                   */

    if(csf > 0.0)				// Only need to act if this
    {						// is a non-triagonal column
      wlp1 = wl+1;				// This is <l+1|W>
      nlp1l = nll + rows_;			// This is <l+1|nmx|l> 
      csf += AbsNorm(*nlp1l); 			// Adjust scaling factor
      sw = 0.0;
      nnil = nlp1l + rows_;			// Start at <l+2|nmx|l>
      for(; nnil<nend; nnil+=rows_)		// Loop down column l using
        sw += square_norm(*nnil/csf);		// lower non-tridiagonals
  
      s = sw + square_norm(*nlp1l/csf);
      s = csf*sqrt(s);
      if(*nlp1l != 0)				// <l+1|W> = 
      {				  	//                     1+s 
        mul(*wlp1, *nlp1l, 1.0+s/norm(*nlp1l));	//    <l+1|mx|l> * ------------
      }					//                 |<l+1|mx|l>|
      else 
      {
        *wlp1 = -s;				// <l+1|W> = -s
      }

      s = sw + square_norm(*wlp1/csf);
      s = sqrt(2/s)/csf;
      nnil = nlp1l + rows_;			// Start at <l+2|nmx|l>
      wi = wlp1 + 1;				// Begin with <l+2|W>
      for(; nnil<nend; wi++,nnil+=rows_)		// Loop column l, lower non-
        mul(*wi, *nnil, s);			// tridiags: <i|W>=s*<i|nmx|l>

      (*wlp1) *= s; 				// <l+1|W> *= s
         
      //                    Transform Matrix From The Left
      //                    This Will Alter The Input Array
 
/* For the current column (l), the diagonal values should already be set.
   The element directly below (tri-diagonal) is determined, as are others
   in below that (all zero).  In the transform, the elements in the lower
   right block are altered to an intermediate state.  This is depicted
   below for the transformation on the 1st and 3rd columns of a 5x5 array
   of <j|hmx|i> (with i>=l, j>l) Are Adjusted.  Note That The Vector W Is
   left untouched in these loops.
 
       t h h h h     t h h h h       t t 0 0 0     t t 0 0 0
       h h h h h     t x x x x       t t t 0 0     t t t 0 0
       h h h h h --> 0 x x x x       0 t t x x --> 0 t t x x
       h h h h h     0 x x x x       0 0 x x x     0 0 t y y
       h h h h h     0 x x x x       0 0 x x x     0 0 0 y y
                 l=0                           l=2
 
   No diagonal elements of the array are altered, nor are any of elements
   above the row with this same column index. Only the lower triangular
   elements of the active column (l) are set to their final values.  The
   rest (in a square block to the lower right) are altered to intermediate
   values corrected in the next transformation (From The Right).  This
   intermediate result will likely not be Hermitian as was the input array
 
   The most direct code to perform this transformation is show below, but this
   is not used in favor of code which is faster (albeit more complicated)

            complex ss;
            int j;
            for(i=l; i<rows_; i++)
              {
              ss = 0.0;
              for(j=l+1; j<rows_; j++) ss += conj(W[j],(*nmx)(j,i));
              for(j=l+1; j<rows_; j++) (*nmx)(j,i) -= ss*W[j];
              }

  This section assumes that  1. n00  = <0|mx|0>     2. wl   = <l|W>
                             3. wend = <rows_|W>+1  4. wlp1 = <l+1|W>        */


      nl0 = n00 + l*rows_;			// Position at <l|nmx|0>
      nlend = nl0 + rows_;			// End of row l
      nli = nl0 + l;				// Start <l|nmx|i>=<l|nmx|l>
      for(; nli<nlend; nli++)			// Loop over nmx row
      {
        wj = wlp1;				// Start at <l+1|W>
        nji = nli+rows_;			// Start at <l+1|nmx|i>
        register double ssi=0, ssr=0;		// Real & imag scaling
        for(; wj<wend; wj++,nji+=rows_)		// Loop down nmx column & add
        {					// <j|W*><j|nmxi> to ss
          register double tr = Re(*nji);	// Re(<j|nmx|i>
          register double ti = Im(*nji);	// Im(<j|nmx|i>
          register double wr = Re(*wj);		// Re(<j|W>) 		
          register double wi = Im(*wj);		// Im(<j|W>)
          ssr += wr*tr + wi*ti;
          ssi += wr*ti - wi*tr;
        }
        wj = wlp1;				// Start at <l+1|W>
        nji = nli+rows_;			// Start at <l+1|nmx|i>
        for(; wj<wend; wj++,nji+=rows_)		// Loop down nmx column & subt.
        {					// ss*<j|W> from <j|nmx|i>
          register double wr = Re(*wj);		// Re(<j|W>)
          register double wi = Im(*wj);		// Im(<j|W>)
          Re(*nji, Re(*nji)-ssr*wr+ssi*wi);	// Modify Re(j|nmx|i>)
          Im(*nji, Im(*nji)-ssr*wi-ssi*wr);	// Modify Im(j|nmx|i>)
        }
      }
             
//                    Transform Matrix From The Right
//                    This Will Alter The Input Array

/* For the current column (l), the vertical elements are already set.
   The transformation adjusts all of the elements to the right of the
   diagonal in the row with the same index as the column. Only the first
   (tri-diagonal) will be non-zero.  In addition, the block to the
   lower right of the diagonal will be altered to an intermediate state,
   except for the next diagonal which will be properly set.  This is depicted
   below for the transformation on the 1st and 3rd columns of a 5x5 array
 
       t h h h h     t t 0 0 0       t t 0 0 0     t t 0 0 0
       t x x x x     t t y y y       t t t 0 0     t t t 0 0
       0 x x x x --> 0 y y y y       0 t t x x --> 0 t t t 0
       0 x x x x     0 y y y y       0 0 t x x     0 0 t t y
       0 x x x x     0 y y y y       0 0 0 x x     0 0 0 y y
                 l=0                           l=2
 
   Intermediate results (array states prior to full tri-diagonalization)
   following the right transformation should be Hermitian.
 
   The most direct code to perform this transformation is show below, but this
   is not used in favor of code which is faster (albeit more complicated)

            complex ss;
            int j;
            for(i=l; i<rows_; i++)                                         
              {
              ss = 0.0;
              for(j=l+1; j<rows_; j++) ss += W[j]*(*nmx)(i,j);
              for(j=l+1; j<rows_; j++) (*nmx)(i,j) -= conj(W[j],ss);
              }

  This section assumes that  1. n00  = <0|mx|0>     2. wl   =<l|W>
                             3. wend = <rows_|W>+1  4. wlp1 = <l+1|W>        */

      ni0 = n00 + l*rows_;				// <i|nmx|0>=<l|nmx|0>
      nnil = ni0 + l;					// <i|nmx|l>=<l|nmx|l>
      for(nnil=ni0+l; ni0<nend; ni0+=rows_,nnil=ni0+l)	// Loop rows >= l
      {
        register double ssi=0, ssr=0;			// Real & imag scaling
        wj = wlp1;					// Start at <l+1|W>
        nij = nnil+1;					// Start at <i|nmx|l+1>
        for(; wj<wend; wj++,nij++) 			// Loop over nmx row i
        {						// and sum <j|W><i|nmx|j>
          register double wr = Re(*wj);			// Re(<j|W>)
          register double wi = Im(*wj);			// Im(<j|W>)
          register double tr = Re(*nij);		// Re(<i|nmx|j>)
          register double ti = Im(*nij);		// Im(<i|nmx|j>)
          ssr += wr*tr - wi*ti;
          ssi += wr*ti + wi*tr;
        }
        wj = wlp1;					// Start at <l+1|W>
        nij = nnil+1;					// Start at <i|nmx|l+1>
        for(; wj<wend; wj++,nij++) 			// Loop over nmx row i & subt.
        {						// ss*<j|W*> from <i|nmx|j>
          register double wr = Re(*wj);			// Re(<j|W>)
          register double wi = Im(*wj);			// Im(<j|W>)
          Re(*nij, Re(*nij)-ssr*wr-ssi*wi);		// New Re(<i|nmx|j>)
          Im(*nij, Im(*nij)+ssr*wi-ssi*wr);		// New Im(<i|nmx|j>)
        }
      }
             
//              Accumulate Transformations
//              (Maintains U As Both Hermitian And Unitary)
 
/* We Are Through Altering The Input Array For This Column Adjustment.  Now
   We Are Altering The Transformation Array U Using The Vector w. The Vector
   w Is Left Unaltered.  For Current Column "l", We Will Adjust All Array U
   Elements To The Right Of The Column.  This Is Depicted Below For A 5x5
   Array On The l=2 (3rd) Column, Where u Are Completed U Elements And x,y
   Represent Elements Not Quite Finished.
 
                           u u u x x    u u u u y
                           u u u x x    u u u u y
                           u u u x x -> u u u u y
                           u u u x x    u u u u y
                           u u u x x    u u u u y
 
   The most direct code to perform this transformation is show below, but this
   is not used in favor of code which is faster (albeit more complicated)
 
                    for(i=0; i<n; i++)
                       { 
                       ss = 0.0;
                       for(j=l+1; j<n; j++) ss += W[j]*U(i,j);
                       for(j=l+1; j<n; j++) U(i,j) -= conj(W[j],ss);
                       }

  This section assumes that  1. U00=<0|U|0>        2. Uend=<rows_|U|rows_>+1
                             3. wend = <rows_|W>+1 4. wlp1 = <l+1|W>         */

      for(Ui0=U00; Ui0<Uend; Ui0+=rows_)		// Loop all U rows
      {
        Uil = Ui0 + l;					// <i|U|l>
        register double ssi=0, ssr=0;			// Real & imag scaling
        for(wj=wlp1, Uij=Uil+1; wj<wend; wj++,Uij++)	// Loop columns l+1 on
        {
          register double wr = Re(*wj);			// Re(<j|W>)
          register double wi = Im(*wj);			// Im(<j|W>)
          register double tr = Re(*Uij);		// Re(<i|U|j>)
          register double ti = Im(*Uij);        	// Im(<i|U|j>)
          ssr += wr*tr - wi*ti;                 	// Add up real scaling
          ssi += wr*ti + wi*tr;				// Add up imag scaling
        }
        for(wj=wlp1,Uij=Uil+1; wj<wend; wj++,Uij++)	// Loop from j=l+1
        {
          register double wr = Re(*wj);			// Re(<j|W>)
          register double wi = Im(*wj);			// Im(<j|W>)
          Re(*Uij, Re(*Uij)-ssr*wr-ssi*wi);		// Set Re(<i|U|j>)
          Im(*Uij, Im(*Uij)+ssr*wi-ssi*wr);		// Set Im(<i|U|j>)
        }
      }
    }
  }    

  h_matrix* htmp = nmx->HMX();		// htmp = tri-diag-herm
  delete [] data;					// Delete our own data
  data = htmp->data;				// Set our data to htmp's
  htmp->data = NULL;				// Set htmp data empty
  delete htmp;						// Now delete htmp
  delete nmx;						// Remove working nmx
  delete [] W;						// Remove temp array
}
 


/*************************************************************************
**									**
** The routine rred will convert a general complex Hermitian tri-	**
** diagonal array into a real symmetric tri-diagonal array.  This 	**
** result may be subsequently used for a complete diagonalization.	**
** The tranformation is diagonal and unitary.				**
**                                                                      **
**      Input           hmx     : Input hermitian array (this)		**
**                                Assumed to be tri-diagonal		**
**                      mxev    : Any previously accumulated            **
**                                transformations in a normal array.    **
**                                If no previous transformations this   **
**                                should be an identity matrix (normal) **
**      Output          hmx     : Real symmetric tri-diagonal array     **
**                      mxev    : All accumulated transformations.      **
**                                                                      **
**  The routine loops over each of the diagonal elements of the tri-    **
**  diagonal array starting with the 2nd column (index 1).  It then     **
**  looks at the element directly left which may be complex.  If found  **
**  this is then made real (as is its symmetric partner), the next two  **
**  off-diagonals set appropriatly, & the transformation matrix fixed.  **
**  Below is a depiction of what the input tri-diagonal array is doing  **
**  during the looping over a 5x5 array.                                **
**                                                                      **
**  x x . . .    x ? . . .    x r . . .    x r . . .    x r . . .       **
**    x x . .      X x . .      x A x .    r x ? . .    r x r x .       **
**      x x . ->     x x . ->   B x x . -> . B X x . -> . r x C . ->    **
**        x x          x x          x x          x x        D   x       **
**          x            x            x            x            x       **
**                                                                      **
**   Initial    Check Row 1     Adjust    Check Row 2     Adjust        **
**                                                                      **
**  Above the x's are just the original hermitian tri-diagonal elements **
**  and the big X's are the diagonals used in the loop.  I used r's to  **
**  depict altered elements which are real and A,B,C,D to depict those  **
**  modified but not yet to the final form.  Note that elements on the  **
**  lower diagonal (B, D, ...) are used as storage locations.           **
**									**
*************************************************************************/

void h_matrix::rred(_matrix* (&U), int newU)

        // Input        hmx     : A tri-diagonal Hermitian array (this)
        //              U	: Previously accumulated transitions
	//		newU	: Flag whether to begin with U=I
	// Output	void    : The two input arrays will be modified.
	//			  The array hmx will return real symmetric 
	//			  tridiagonal and the array U will be
	//			  adjusted to contain the additional	
	//			  transformation from hmx in to hmx out.
	// Note			: This routine assumes that U is an
	//			  n_matrix and USES ITS INTERNAL STRUCTURE
	//			  It also assumes U is square and its
	//			  dimension matches that of hmx

  {
  int i;
  double normod;			// For norm of off-diagonal
  complex sf;				// Transformation scaling factor
  complex *hii = data + rows_;		// Diagonal element: hii = <1|hmx|1>
  complex *hup = data + 1;		// Elem. above diag: hup = <0|hmx|1>
  complex *hrt = hii + 1;		// Elem. right diag: hrt = <1|hmx|2>
  if(newU)				// Allocate array U if from scratch
    U=new n_matrix(rows_,rows_,complex0);
  complex *U00 = ((n_matrix*)U)->data;	// U00 = <0|U|0>
  complex *Uend = U00 + rows_*rows_;	// Uend = <rows_|U|rows_> + 1
  if(newU)				// If U is to begin from scratch
    {					// we must set it to I
    complex *Uii = U00;			// <i|U|i> <- <0|U|0>
    for(; Uii<Uend; Uii+=rows_+1)	// Initialize U to I matrix
    *Uii = complex1;
    }
  complex *Uji;				// Uji  = <j|U|i>
  int rl = rows_ - 1;			// Length of row containing hii
  for(i=1; i<rows_; i++)		// Loop rows of hmx (after row 0)
    {
    normod = norm(*hup);		// If off-diag. non-zero <i|hmx|i-1>
    if(normod!=0) 			// then we transform it to be real
      {
      sf = (*hup)/normod;		// 	Scaling factor
      (*hup) = normod;			// 	Set off-diagonal to new real
      if(i < rows_-1) (*hrt) *= sf;	// 	Adjust the next off diagonal
      for(Uji=U00+i;Uji<Uend;Uji+=rows_)// 	Accumulate this transform
        *Uji = conj(sf, *Uji);
      }
    hup = hrt;				// Right becomes up on next row 
    hii += rl;				// Jump to next row's diagonal
    hrt = hii + 1;			// New element right of diagonal hii
    rl--;				// New row length
    }
  }


/*************************************************************************
**									**
** The function sign is a C++ implementation of the standard FORTRAN	**
**  library function.							**
**									**
*************************************************************************/

inline double sign(double a,double b)
  { return double((b<0.0)?-fabs(a):fabs(a)); }

/*************************************************************************
**									**
** The routine tqli uses the QL algorithm for the diagonalization of	**
** a real symmetric tridiagonal matrix.					**
**									**
** This function uses the "sign" function!				**
**									**
*************************************************************************/

//sosi- this function appears to be loosing memory!!!!!

void h_matrix::tqli(_matrix* (&U), _matrix* (&D), int newU)

        // Input        hmx     : A tri-diagonal real symmetric array (this)
        //              U	: Pointer to mx to become eigenvectors mx
        //              D	: Pointer to mx to become diagonal mx
        //              newU    : Flag whether to begin with U=I
        // Output       void    : The matrices D and U are set to a diagonal
        //                        matrix of nmx eigenvalues and a matrix
        //                        the matrix of D eigenvectors respectively
	// Note			: The input array hmx is used only to supply
	//			  the diagonals & tri-diagonal elements
	// Note			: The pointer D must NOT point to any
	//			  allocated memory.  Memory will be allocated
	//			  herein and D set to be d_matrix* (diagonal)

  {
  if(!rows_) return;			// Do nothing for 0x0 array

//	       Copy Input Real Symmetric Tri-Diagonal Elements
// Diagonals Go Into Diagonal Matrix D, Off-Diagonals Into Vector e[n]

  D = new d_matrix(rows_,rows_);	// Allocate space for eigenvalues
  double* odvals;
  odvals = new double[rows_];
//  double odvals[rows_];			// Array for off-diagonals
  complex *hii = data; 			// hii = <i|hmx|i> -> <0|hmx|0>
  if(newU)				// Allocate array U if from scratch
    U=new n_matrix(rows_,rows_,complex0);
  complex *u00=((n_matrix*)U)->data;	// u00 = <0|U|0>
  complex *uend = u00 + rows_*rows_;	// End of U: after <rows_|U|rows>
  if(newU)				// If U is to begin from scratch
    {					// we must set it to I
    complex *uii = u00;			// <i|U|i> <- <0|U|0>
    for(; uii<uend; uii+=rows_+1)	// Initialize U to I matrix
    *uii = complex1;
    }
  complex *dst = ((d_matrix*)D)->data;	// dst = <0|D|0>
  complex *dii = dst;			// dii = <i|D|i> -> <0|D|0>
  int rl = rows_;			// Length of hmx row ii
  for(int ii=0; ii<rows_-1; ii++, dii++)
    {
    *dii = *hii;			// Store diagonal of row ii
    odvals[ii] = zRe(*(hii+1));		// Store off-diagonal of row ii
    hii += rl;				// Skip hii to next diagonal
    rl--;				// Decrement row length
    }
  odvals[rows_-1] = 0.0;		// Last one set to zero
  *dii = *hii;				// Last one <i|D|i> = <i|hmx|i>

  int i,l,m,iter;
  double dd,s,r,p,g,f,c,b,x;
  complex sxmki,cxmki,sxmkip1,cxmkip1;	// Complex temporaries

  complex *dmm, *dll, *diip1;		// Needed to get D elements
  complex *u0i, *uki, *ukip1;		// Needed to get U elements
  iter = 0;				// Zero iteration count
  for(l=0, dll=dst; l<rows_; l++,dll++)	// Loop for eigenvalues
    {
    do
      {
      m = l;				// Start with row m
      dmm = dll; 			// Pointer to <m|D|m>
      while(1)	 			// Look for an eigenvalue
        {
        if(m == rows_-1) break;		// all eigenvalues found
        dd = fabs(zRe(*dmm)) + fabs(zRe(*(dmm+1)));
        if(fabs(odvals[m])+dd == dd) break;// A new eigenvalue found
        m++;
        dmm++;				// Adjust to row <m|D|m>
        }
      if(m != l)
        {
        if(iter > 10*rows_) 		// Bail if too many iterations!
          {
          HMxerror(19, 1);		//   Problems during diag.
          HMxerror(20, 1);		//   Too Many Iterations
          HMxerror(55, 1);		//   Unlikely Error
          HMxerror(53, 1);		//   Please Report Problem
          HMxfatal(54);			//   GAMMA WWW site
          }
        iter++;				// Increment iteration count

// 				Form The Shift

//                                        1/2
//     [ <l+1|hmx|l+1> - <l|hmx|l>       ]
// r = | ------------------------- + 1.0 |       
//     [       2*<l|hmx|l+1>             ]

        g = (Re(*(dll+1))-Re(*dll))/(2.0*odvals[l]);
        r = sqrt(g*g + 1.0);
        x = odvals[l]/(g+sign(r,g));
        g = Re(*dmm) - Re(*dll) + x;

//                           This is The QL Step

        s = 1.0;					// Initialize s
        c = 1.0;					// Initialize c
        p = 0.0;					// Initialize p
        for(i=m-1,dii=dmm-1,diip1=dmm; i>=l; i--,dii--,diip1--)
          {
          f = s*odvals[i];
          b = c*odvals[i];
          if(fabs(f) >= fabs(g))
            {
            c = g/f;
            r = sqrt(c*c+1.0);
            odvals[i+1] = f*r;
            s = 1.0/r;
            c = c*s;
            }
          else
            {
            s = f/g;
            r = sqrt(s*s + 1.0);
            odvals[i+1] = g*r;
            c = 1.0/r;
            s = s*c;
            }
          g = Re(*diip1) - p;
          r = (Re(*dii)-g)*s + 2.0*c*b;
          p = s*r;
          Re(*diip1, g+p);			// <i+1|dmx|i+1> += g+p
          g = c*r - b;

//          Accumulate Transforms In the Eigenvector Array U
//   This Will Adjust Columns i and i+1 For This Interation in the Loop

          u0i = u00 + i;			// Set <0|U|i>
          uki = u0i;				// Set <k|mkev|i> = <0|U|i>
          ukip1 = uki+1;			// Set <k|mkev|i+1> = <0|U|i+1>
          for(uki=u0i, ukip1=uki+1; uki<uend;
                       uki+=rows_,ukip1+=rows_)
            {
	    mul(sxmki, s, *uki);		// sxmki = s*<k|U|i> (Need A Copy!)
	    mul(cxmki, c, *uki);		// cxmki = c*<k|U|i>
            mul(sxmkip1, s, *ukip1);		// sxmkip1 = s*<k|U|i+1>
            mul(cxmkip1, c, *ukip1);		// cxmkip1 = c*<k|U|i+1>
            sub(*uki, cxmki, sxmkip1);		// <k|U|i> = c*<k|U|i> - s*<k|U|i+1>
//            add(*mkip1, sxmki, cxmkip1);
// sosi - This is a little glitch.  GCC won't pick up on class complex "add" function
//        instead it thinks it's a bad call to h_matrix add......  I don't even know
//	  why the :: in front fixes things....
            ::add(*ukip1, sxmki, cxmkip1);
            }
          }
        Re((*dll), Re(*dll)-p);			// <l|dmx|l> -= p
        odvals[l] = g;				// Set adjusted <l|hmx|l+1>
        odvals[m] = 0.0;			// We've zeroed <m|hmx|m+1>
        }
      } while (m!=l);				// Loop till all eigenvalues found
    }
  delete [] odvals;
  }


/*************************************************************************
*									**
** The routine diag combines the routines cred,	rred and tqli (using	**
** the function sign) to diagonalize a complex Hermitian matrix.  	**
**									**
** During the course of the transformations, exclusively unitary	**
** (orthogonal) transformations are used. Errors come mainly from	**
** rounding. The resulting eigenbase of a Hermitian array deviates 	**
** from unitarity only negligibly.					**
**									**
** This is THE Entry Point to Hermitian Array Diagonalizations!		**
** This Is The Workhorse Behind Many GAMMA Calculations............!	**
** This Routine Needs to Be Working Very Well......................!	**
**									**
** Note That There Is Currently One Unique Feature Of This Function	**
** That Sets It Apart From Most Other GAMMA Matrix Functions. It Sets	**
** The Pointers mxd and mxev To Point To A Diagonal And A Normal Array	**
** Respectively.  In Turn, The Function Calls Diagonalization Routines	**
** That EXCLUSIVELY TAKE d_matrix and n_matrix!				**
**									**
*************************************************************************/

void h_matrix::diag(_matrix* (&D), _matrix* (&U))

// Input        hmx     : h_matrix (this)
//              D	    : Pointer to mx to become diagonal mx
//              U    	: Pointer to mx to become eigenvectors mx
//
// Output       void    : The matrices D and U are set to the diagonal
//                          matrix of hmx eigenvalues and the matrix
//                          of dmx eigenvectors respectively.
//
// Note 1               : Both U and D MUST be input pointing to
//                          NULL matrix. Their data will be allocated in
//                          cred and tqli respectively. Their reference
//                          counting MUST be handled externally.
// Note 2               : D will be a d_matrix, U n_matrix unitary

{

#if defined(_USING_SUNPERFLIB_) || defined(_USING_LAPACK_)

    int nrows = this->rows();
    int ncols = this->cols();
    if(nrows < 32)  // then do it the old fashioned gamma way.
    {
        h_matrix H(*this);			    // Make a copy of hmx to use
        H.cred(U);				        // To Hermitian tridiagonal form
        H.rred(U);				        // To real symm. tridiagonal form
        H.tqli(U, D);				    // To real diagonal form
        ((n_matrix*)U)->unitary = true;	// This should be unitary now
//	std::cerr << "h_matrix diag: using gamma code\n";
    }
    else
    {
	n_matrix* hmx =	new n_matrix(nrows, ncols);
	this->convert(hmx);
  	n_matrix* hmx1 = (n_matrix *)hmx->transpose();
#ifdef _USING_SUNPERFLIB_
        char jobz = 'V';  // Calculate "eigenvectors" AND "eigenvalues".
        char uplo = 'U';  // Upper triangular...
        int  N    = nrows;
        int  lda  = nrows;
        double *w_eig = new double[N];
        int info = -55555;
        zheev(jobz, uplo, N, (doublecomplex *)hmx1->data, lda, w_eig, &info);
#endif
#ifdef _USING_LAPACK_
        char jobz = 'V';  // Calculate "eigenvectors" AND "eigenvalues".
        char uplo = 'U';  // Upper triangular...
#if defined(__LP64__) /* this is needed on MacOS CLAPACK to make types match*/
        int  N    = nrows;
        int  lda  = nrows;
        int info = -55555;
        int lwork = 2*N + 10;
#else
        long int  N    = nrows;
        long int  lda  = nrows;
        long int info = -55555;
        long int lwork = 2*N + 10;
#endif
        double *w_eig = new double[N];
        complex *work= new complex [2*lwork];
        double *rwork = new double [3*N-2];
        zheev_(&jobz, &uplo, &N, (__CLPK_doublecomplex *) hmx1->data, &lda, w_eig, (__CLPK_doublecomplex *) work, &lwork, rwork, &info);
#endif
        if(info == 0)
        { // all is well.
        }
        else if (info > 0)
        { std::cerr << "\nDiagonalization failed to converge\n";
        }
        else if(info == -55555)
        { std::cerr << "\nReturn value, 'info', does not appear to have been set\n";
        }
        else
        { std::cerr << "\ninfo = " << info << "\n";
        }

        // Create the matrixes D and U.
        D = new d_matrix(nrows, ncols);
        U = new n_matrix(nrows, ncols);

        // Copy results back into D and U.
	complex *hmxp = hmx1->data;
        for(int i=0; i<nrows; i++)
        { D->put( w_eig[i], i, i);
          for(int j=0; j<ncols; j++)
          { U->put(hmxp[j*nrows+i], i, j);
          }
        }
        ((n_matrix*)U)->unitary = true;	// This should be unitary now
	delete hmx;
	delete hmx1;
	delete [] w_eig;
#ifdef _USING_LAPACK_
	delete [] work;
	delete [] rwork;
#endif
//	std::cerr << "h_matrix diag: using sunperflib code\n";
    }
#else
    h_matrix H(*this);			    // Make a copy of hmx to use
    H.cred(U);				        // To Hermitian tridiagonal form
    H.rred(U);				        // To real symm. tridiagonal form
    H.tqli(U, D);				    // To real diagonal form
    ((n_matrix*)U)->unitary = true;	// This should be unitary now
#endif
    return;				            // Now hmx = U*D*adj(U)
}


/*                         Diagonalization Routines

  Note That The Returned Diagonal & Eigenvectors Arrays Must Have Their 
  Referencing Set External To This Class.

           Input        hmx     : A Hermitian matrix (this)
                        mxd     : Pointer to mx to become diagonal mx
                        mxev    : Pointer to mx to become eigenvector mx
           Output       void    : The matrices mxd and mxev are set to
                                  the Hermitian matrix of dmx eigenvalues and
                                  the matrix of dmx eigenvectors respectively
           Note                 : For h_matrix, mxd should be real and 
                                  mxev both unitary and Hermitian
           Note                 : The reference count to dmx must be twice
                                  incremented external by the call origin.   */

// sosi - this is not optimized for h_matrix at all

//#include <sys/time.h>
//#include <sys/resource.h>

std::vector<int> h_matrix::BlockDiag(_matrix* (&BD), std::vector<int> &U)
  {
//cout << "This is h_matrix block diag\n";
//struct rusage me;
//getrusage(0, & me);
//cout << "block-diag routine: point 01: " << me.ru_utime.tv_sec+me.ru_utime.tv_usec/1.0e6 << " seconds\n";
//cout.flush();
  int nr = rows_;				// Matrix dimension
  //int count = 0;				// Count permutations
  BD = new h_matrix(*this);			// Start with copy

  int i,nze=0,bs=0,t5,t6;
  for(i=0; i<nr; i++) { U.push_back(i); }	// Set both as unpermuted
  std::vector<int> blkdims;
  for(t5=0;t5<nr;++t5)
  { if(nze==t5)
      ++nze;
    for(t6=nze;t6<nr;++t6)
    { if((*BD).get(t5,t6) != complex0 || (*BD).get(t6,t5) != complex0)
      { if(t6 != nze)
        { complex za, zb;
          int z1;
	  int t7;
//	  int k1, k2;
///       getrusage(0, & me);
///       std::cout << "block-diag routine: before perm: " << me.ru_utime.tv_sec+me.ru_utime.tv_usec/1.0e6 << " seconds\n";
//	  std::cout << "permuting "<< nze << " and " << t6 << "\n";
//	  for(k1=0;k1<nr;++k1)
//	  { for(k2=0;k2<nr;++k2)
//	       std::cout << (*BD).get(k1,k2) << "  ";
//	    std::cout << "\n";
//	  }
//        std::cout.flush();
	  for(t7=0;t7<nr;++t7)
	  { if(t7!=t6 && t7!=nze)
	    { za=(*BD).get(t7,t6);
	      zb=(*BD).get(t7,nze);
	      (*BD).put_h(zb,t7,t6);
	      (*BD).put_h(za,t7,nze);
	    }
	  }
	  za=(*BD).get(t6,t6);
	  zb=(*BD).get(nze,nze);
	  (*BD).put(zb,t6,t6);
	  (*BD).put(za,nze,nze);
//	  for(k1=0;k1<nr;++k1)
//	  { for(k2=0;k2<nr;++k2)
//	       std::cout << (*BD).get(k1,k2) << "  ";
//	    std::cout << "\n";
//	  }
	  z1=U[t6];
	  U[t6]=U[nze];
	  U[nze]=z1;
//        getrusage(0, & me);
//        std::cout << "block-diag routine: after  perm: " << me.ru_utime.tv_sec+me.ru_utime.tv_usec/1.0e6 << " seconds\n";
//        std::cout.flush();
	}
	++nze;
      }
    }
    if(nze==(t5+1) || t5==(nr-1))
    { blkdims.push_back(t5+1-bs);
      bs=t5+1;
    }
  }
//getrusage(0, & me);
//std::cout << "block-diag routine: point 02: " << me.ru_utime.tv_sec+me.ru_utime.tv_usec/1.0e6 << " seconds\n";
//std::cout.flush();
  return blkdims;
  }

void h_matrix::HermTriDiag(_matrix* (&HTD), _matrix* (&U))
  {
  if(is_tridiagonal()) 			// If we are already tri-diagonal
    {					// we don't need to do anything
    HTD = this;				//   Tri-diagonal form is us
    U   = new i_matrix(rows_, rows_);	//   Transformation matrix is I
    return;
    }
  h_matrix* hmx = new h_matrix(*this);	// Make a copy of hmx to use
  hmx->cred(U);				// To Hermitian tridiagonal form
  HTD = hmx;				// Assign pointer (dont del hmx data)
  }
 
void h_matrix::SymTriDiag(_matrix* (&STD), _matrix* (&U))
  {
  if(is_real()) 			// If already real tri-diagonal
    {					// we don't need to do anything
    STD = this;				//   In symmetric tri-diag. form
    U   = new i_matrix(rows_, rows_);	//   Transformation matrix is I
    return;
    }
  h_matrix* hmx = new h_matrix(*this);  // Make a copy of hmx to use
  hmx->rred(U, 1);			// To symmetric tridiagonal form
  STD = hmx;				// Assign pointer (dont del hmx data)
  }
 
void h_matrix::SymDiag(_matrix* (&D), _matrix* (&U)) { tqli(U, D, 1); }
 

// ____________________________________________________________________________
// O                    CLASS H_MATRIX INVERSION FUNCTIONS
// ____________________________________________________________________________

/* The inverse of a Hermitian matrix can be obtained in several ways.  For
   example (using H=Hermitian & U=unitary)

                    -1                t   -1       -1    t   -1
               H * H   = I = U * D * U * H   ==>  H   = U * D  * U

   But herein, we will use the LU decomposition to obtain the inverse. In this
   case we have (using H=Hermitian, L=lower triang., & U=upper triang.)

                               -1               -1
                          H * H  = I = L * U * H
  
   where we 1st solve L*Y=I for Y, then U*inv(H)=Y for inv(H) taking advantage
   of the special formats of L & U.

  The LU & LU inverse routines used here are derived from "Numerical Recipies",
  W.H.  Press, B.P.Flannery, S.A.Teukolsky, & W.T.Vetterling, Cambridge 
  University Press, 1989.  See pages 31 through 36. The code was extensively
  adapted for C++ and GAMMA arrays.

  The LU algorithm takes as input a Hermitian A and an integer array indx.
  It returns returns the LU decomposition of A' where A' is A with any needed
  row permutations. The latter (row permutations) are stored in the integer
  array indx (indx[0] < 0 flags no row changes).

  The LU decomposition formulates the equation (depicted for 4x4 A below) 

                                A = LU

       [A11 A12 A13 A14]   [ 1   0   0   0 ] [U11 U12 U12 U14]
       |A21 A22 A23 A24]   [L21  1   0   0 ] [ 0  U22 U23 U24|
       |A31 A32 A33 A34] = [L31 L32  1   0 ] [ 0   0  U33 U34|
       [A41 A42 A43 A44]   [L41 L42 L43  1 ] [ 0   0   0  U44]

  Neglecting implicit <i|L|i>, L & U are stored in in a single array.

                           [U11 U12 U12 U14]
                           [L21 U22 U23 U24|
                           [L31 L32 U33 U34|
                           [L41 L42 L43 U44]

  We have dropped the diagonal elements of L since we know they are always 1.
  The algorithm uses Crouts method in determining the elements of L and U
  which are stable under pivoting (row interchanges).

  The LU inverse algorithm takes an input matrix A' in its LU decomposition
  format where A' is some original matrix A with any needed row permutations
  to attain the LU form. The row permutations used to relate A to A' are stored
  in the integer array indx.  The function then proceeds to solve the following
  equation for |x> given a |b>:
                                 A|x> = |b>

  This is accomplished by considering (reformulating) the problem as

                    A|x> = (LU)|x> = L(U|x>) = L|y> = |b>

  using the LU decomposition. This is first solved for the vector |y>, where
  L|y> = |b>, easily accomplished since L is lower-triangular.  Subsequently
  the equation U|x> = |y> is solved for the desired |x>, again easily done
  since U is upper=triangular.  By repeating this process for many |x> and |b>
  pairs of columns, the inverse of the original array can be obtained.
                                    -1
        A [|xo> |x1> ... |xn>] = A A  =  [|bo> |b1> .... |bn>] = I 


    Function  Arguments               Result                      Notes
  ---------- --------- ----------------------------------- -------------------
     inv        ---    inverse(hmx) where inv(hmx)*hmx = I
     LU        int*    LU == hmx + row permutations        One Array L\U Form
     LUinv    int*, B  Solutionh X where LU*X = AX = B     Inverse if B=I


        A [|xo> |x1> ... |xn>] = A A  =  [|bo> |b1> .... |bn>] = I          */

_matrix* h_matrix::inv()
  {
// sosi - this will use up memory
  int hd = rows_;			// Matrix dimension
  i_matrix* I = new i_matrix(hd,hd);	// Construct an Identity matrix
  int *indx;
  indx = new int[hd];
//  int indx[hd];				// Array for any row permutations
  n_matrix* hLU = (n_matrix*)LU(indx);	// Perform LU decomposition of hmx
  _matrix* hinvmx = I->LUinv(indx, hLU);
  delete [] indx;
  return hinvmx; 
//  return I->LUinv(indx, hLU);		// Return the inverse via LU*hinv=I;
  }

_matrix* h_matrix::LU(int *indx)
  {
  n_matrix* nLU = new n_matrix(rows_,rows_);	// New n_matrix=this, start nLU
  convert(nLU);					// Copy hmx into nLU to start
  ((n_matrix*)nLU)->LU_decomp(indx);		// LU Decmp. of row permuted nmx
  return nLU;					// Return Ainv, a new n_matrix
  }

_matrix* h_matrix::LUinv(int *indx, _matrix* LU)
  {
  if(!LU->is_square())
    {
    HMxerror(17, 1);                                    //   Dimension problems
    HMxerror(6, " LU Inversion On Rectangular Mx", 1);  //   LU inversion multiply problems
    HMxfatal(3, "LUinv");                               //   Fail in times_adjoint function
    }
  if(rows_ != LU->rows())                               // Insure LU & D dimensions match
    {
    HMxerror(17, 1);                                    //   Dimension problems
    HMxerror(6, " LU Inversion LU|X>=|H> Mismatch", 1); //   LU inversion multiply problems
    HMxfatal(3, "LUinv");                               //   Fail in times_adjoint function
    }
  switch(LU->stored_type())
    {
    case i_matrix_type:				// LU is an identity matrix
      if(indx[0] < 0)				// If no row permutations, LU = A = I
        return this;				// so then A*X=B -> I*X=B; X=B
      break;
    case d_matrix_type:				// LU is a diagonal matrix
      if(indx[0] < 0)				// If no row permutations, LU = A = dmx 
        {
        d_matrix* LUinv = (d_matrix*)LU->inv();	// First get inv(LU) = inv(dmx)
        return LUinv->multiply(this);		// Then as D*X=B -> dmx*X=B; X=inv(dmx)*B;
        }
      break;
    case h_matrix_type:				// LU is a Hermitian matrix
    case n_matrix_type:				// LU is a normal matrix
    default:
      break;
    }
  n_matrix* X = new n_matrix(rows_,cols_);	// Constuct a new n_matrix
  convert(X);					// Convert B to X as n_matrix
  return X->LUinv(indx,LU);			// Use n_matrix algorithm
  }

// ____________________________________________________________________________
// P                      CLASS H_MATRIX I/O FUNCTIONS
// ____________________________________________________________________________

/* These functions perform both ASCII and Binary input/output operations on a
   Hermitian matrix.
 
    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   print      ostream     Writes array elemnts in ASCII to output stream
   write      ofstream    Writes array elemnts in Binary to output file stream
   read       ifstream    Read in array from Binary input file stream
   readASC    istream     Read in array from an input stream
   ask          ----      Interactively requests Hermitian matrix from user
   
   The binary format of h_matrix has data only, not size information.  The size
   is taken care of in classes matrix/_matrix and should done prior to use of
   these functions.  The data ordering is Re(<i|hmx|j>, Im(<i|hmx|j> columns
   then rows (i.e. row by row.)  If the flag=0 then only the diagonals are
   output (GAMMA format) whereas a non-zero flag will force all elements
   (diagonals and zero off-diagonals) to be output.                          */


std::vector<std::string> h_matrix::printStrings(const MxPrint& PFlgs) const
  {
  std::vector<std::string> PStrings;		// What to return

  int ptype = 0;                                // Complex elements output
  if(PFlgs.MxRIPrnt)                            // If we want just reals
    {                                           // or just imags, check
    if(is_real())           ptype = 1;          // Real elements output
    else if(is_imaginary()) ptype = 2;          // Imag elements output
    }

  int    elen;                                  // A single element length
  switch(ptype)                                 // Get the element length
    {                                           // from class complex. This 
    default:                                    // depends on the output
    case 0: elen = complex::zlength(); break;   // format
    case 1:
    case 2: elen = complex::dlength(); break;
    }
  std::string ezer("");				// String for "hidden" element
  if(!PFlgs.MxAll)				// If not printing all
    ezer = std::string(elen, ' ');		// elems, this is invisible
  std::string pline;
 
  std::string efmt = complex::dformat();        // Real/Imag element format
  int i,j, pos;
  for(i=0,pos=0; i<rows_; i++ )                 // Loop array rows
    {
    pline = std::string("");			// Empty row string

    if(!PFlgs.MxAll)				// If not printing all elems
      for(j=0; j<i; j++) pline += ezer;		// put blanks left of diag
    else if(ptype==1)                           // Output left of diag
      for(j=0; j<i; j++) 			// as real elements
        pline += MxModform(efmt.c_str(),Re(get(i,j)));
    else if(ptype==2)                           // Output left of diag
      for(j=0; j<i; j++) 			// as imaginary elements
        pline += MxModform(efmt.c_str(),Im(get(i,j)));
    else                                        // Output left of diag
      for(j=0; j<i; j++)			// as complex elements
        pline +=  get(i,j).printString();

    if(ptype==1)                                // Output rest as
      for(j=i; j<cols_; j++, pos++) 		// real elements
        pline += MxModform(efmt.c_str(),Re(data[pos]));
    else if(ptype==2)                           // Output rest as
      for(j=i; j<cols_; j++, pos++) 		// imaginary elements
        pline += MxModform(efmt.c_str(),Im(data[pos]));
    else 					// Output rest as
      for(j=i; j<cols_; j++, pos++) 		// complex elements
        pline +=  data[pos].printString();
    PStrings.push_back(pline);
    }
  return PStrings;
  }

std::vector<std::string> h_matrix::pictureStrings(const MxPrint& PFlgs) const
  {
  if(PFlgs.MxAll) { return printStrings(PFlgs); }	// If all elems
  std::vector<std::string> PStrings;			// What to return
  return PStrings;
  }

void h_matrix::print(std::ostream& ostr, const MxPrint & PF) const
  {
  int ptype = 0;                                // Complex elements output
  if(PF.MxRIPrnt)                               // If we want just reals
    {                                           // or just imags, check
    if(is_real())           ptype = 1;          // Real elements output
    else if(is_imaginary()) ptype = 2;          // Imag elements output
    }

  int    elen;                                  // A single element length
  switch(ptype)                                 // Get the element length
    {                                           // from class complex. This 
    default:                                    // depends on the output
    case 0: elen = complex::zlength(); break;   // format
    case 1:
    case 2: elen = complex::dlength(); break;
    }
  std::string ezer("");				// String for "hidden" element
  if(!PF.MxAll)				        // If not printing all
    ezer = std::string(elen, ' ');		// elems, this is invisible
//    ezer = string(elen+1, ' ');		// elems, this is invisible
 
  std::string efmt = complex::dformat();        // Real/Imag element format
  int clen = 40 - ((elen+1)*rows_-1)/2;         // Space to center 1 line
  std::string sp("");                           // Spacer to center a line
  if(clen>0) sp = std::string(clen, ' ');       // Set spacer for centering
  int i,j, pos;
  for(i=0,pos=0; i<rows_; i++ )                 // Loop array rows
    {
    ostr << sp;                                 // Space to center
    if(!PF.MxAll)				// If not printing all elems
      for(j=0; j<i; j++) ostr << ezer;          // put blanks left of diag
    else if(ptype==1)                           // Output left of diag
      for(j=0; j<i; j++) 			// as real elements
        ostr << MxModform(efmt.c_str(),Re(get(i,j)));
    else if(ptype==2)                           // Output left of diag
      for(j=0; j<i; j++) 			// as imaginary elements
        ostr << MxModform(efmt.c_str(),Im(get(i,j)));
    else                                        // Output left of diag
      for(j=0; j<i; j++) ostr << get(i,j); 	// as complex elements

    if(ptype==1)                                // Output rest as
      for(j=i; j<cols_; j++, pos++) 		// real elements
        ostr << MxModform(efmt.c_str(),Re(data[pos]));
    else if(ptype==2)                           // Output rest as
      for(j=i; j<cols_; j++, pos++) 		// imaginary elements
        ostr << MxModform(efmt.c_str(),Im(data[pos]));
    else 					// Output rest as
      for(j=i; j<cols_; j++, pos++) 		// complex elements
        ostr << data[pos];
    ostr << "\n";                               // Skip to next line
    }
  }

void h_matrix::picture(std::ostream& ostr, const MxPrint & PFlgs) const
  {
  int rlen = 2*rows_-1;				// Length of 1 row
  std::string sp("");				// Spacer to center row
  int len = 40-rlen/2;				// Space to center row
  if(len>0) sp = std::string(len, ' ');		// Set spacer to center
  int i, j, pos;				// Looping indices
  for(i=0, pos=0; i<rows_; i++ )		// Loop array rows
    {
    ostr << sp;					// Center row
    if(PFlgs.MxAll)				// If showing all elems
      for(j=0; j<i; j++)			// then draw in the
        if(norm(get(i,j))) ostr << "x ";	// those on lower triang.
        else               ostr << "0 ";
    else					// If not showing all
      for(j=0; j<i; j++) ostr << "  "; 		// then thse are empty
    for(j=i; j<cols_; j++)			// Now fill in those
      {
      if(norm(data[pos])) ostr << "x ";
      else                ostr << "0 ";
      pos++;					// Increment element
      }
    ostr << "\n";				// Skip to next row
    }
  }

        // Input            hmx : An h_matrix (this)
        //                   fp : A file stream
        //                 form : Flag for format in which to write hmx data
        //                          0 = GAMMA format (default); hmx data only
        //                          !0 = All matrix elements including 0's 
        // Note                 : For form!=0 data is written in rows, each pt
        //                        as two floats, real followed by imaginary
        //                        Re(<0|hmx|0>), Im(<0|hmx|0>, Re(0|hmx|1>, ...
        // Note                 : For form==0, only upper-triang. elements are
	//			  written, each pt as two floats, real then imag
	//			  Re(<0|hmx|0>),Im(<0|hmx|0>),Re(<1|hmx|1>),...

void h_matrix::write(std::ofstream& fp, int form)
  {
	// *** changed float to double.
  double dr,di;
  int i,j,pos,cpos;
  for(i=0, pos=0; i<rows_; i++)
    {
    for(j=0; j<i&&form; j++)
      {
      cpos = j*cols_-(j*(j-1))/2+i-j;
      dr = Re(data[cpos]);
      di = -Im(data[cpos]);
      fp.write((char*)&dr, sizeof(double));
      fp.write((char*)&di, sizeof(double));
      }
    for(j=i; j<cols_; j++,pos++)
      {
      dr = Re(data[pos]);
      di = Im(data[pos]);
      fp.write((char*)&dr, sizeof(double));
      fp.write((char*)&di, sizeof(double));
      }
    }
  }
 
 
void h_matrix::read(std::ifstream &fp)
 
        // Input            hmx : A h_matrix (this)
        //                   fp : A file stream
        // Output          void : Matrix is filled with data in file fp
        // Note                 : Read in binary format, data only
	//			  i.e. only the upper-triangular elements
        //                        Data ordering: Re(<i|mx|i>), Im(<i|mx|i>)
 
  {
  double dr,di;
  for(int pos=0; pos<size; pos++ )
    {
    fp.read((char*)&dr, sizeof(double));
    fp.read((char*)&di, sizeof(double));
    data[pos] = complex(dr,di);
    }
  }
     

void h_matrix::readASC(std::istream& istr)
  
        // Input            hmx : A h_matrix (this)
        //                 istr : An input stream
        // Output          void : The function fills the matrix hmx
        //                        from data in the input stream.  The
        //                        input begins with the row and column
        //                        dimensions (which should specify dmx.
	//                        Matrix nmx is modified.
        // Note                 : Function is for INTERACTIVE programs
// sosi - this is never accessable.  Class matrix always routes through n_matrix

  {
  int i,j;
  istr >> i >> j;
  resize(i,j);
  for(int pos=0; pos<size; pos++)
    istr >> data[pos]; 
  }


void h_matrix::ask( )

        // Input            hmx : A h_matrix (this)
        // Output          void : The function sends questions to
        //                        standard output interactively asking
        //                        the user to supply the matrix elements
        //                        to specify hmx.  Matrix hmx is modified.
        // Note                 : Function is for INTERACTIVE programs

  {
  float dr,di;
  for(int i=0, pos=0; i<rows_; i++ )
    {    
    std::cout << "\n\tPlease Input Value of <"
         << i << "|mx|" << i << "> [re]: ";
    std::cin  >> dr;
    data[pos] = complex(dr,0);
    pos++;
    for(int j=i+1; j<cols_; j++,pos++)
      {
      std::cout << "\n\tPlease Input Real and Imaginary Value of <"
         << i << "|mx|" << j << "> [re im]: ";
      std::cin >> dr >> di;
      data[pos] = complex(dr,di);
      }
    }
  }  

// ____________________________________________________________________________
// M                 CLASS H_MATRIX MISCELLANEOUS FUNCTIONS
// ____________________________________________________________________________


void h_matrix::resize(int i, int j)

	// Input          hmx   : An h_matrix (this)
	//		  i,j   : Row, column dimensions
	// Output         none  : Size of hmx readjusted to ixj
	// Note			: For hmx, failure if i!=j
	// Note			: Destroys data if not the same size

  {
  if(i != j)                    // Insure resized to square array
    {
    HMxerror(17, 1);                            //   Dimension problems
    HMxerror(6, " Matrix Resize", 1);           //   Matrix multiply problems
    HMxfatal(3, "resize");                      //   Fail in times_adjoint function
    }
  if((i*i+i)/2 != size)		// See if the new size is different
    { 				// If it is different, then
    _matrix::resize(i,j);	//	set new matrix sizes
    delete [] data;		//	delete the current data
    size = (i*i+i)/2;		// 	set the new data size
    data = new complex[size];	//	delcare a new data array
    }
  }


_matrix* h_matrix::copy() { return new h_matrix(*this); }

        // Input        hmx     : An h_matrix (this)
        // Output       mx      : A copy of the input h_matrix
        // Note                 : This allocates more memory
 


void h_matrix::convert(_matrix* mx)	

	// Input          hmx   : An h_matrix (this)
	//		  mx    : Matrix of a particular type
        // Output         mx    : Elements of hmx are placed into mx
        //                        in a fashion so as to maintain the mx
        //                        matrix type.  The result is the "best"
        //                        matrix conversion of hmx will reside in mx
        // Note                 : Since a hermitian matrix cannot always be converted
        //                        an identity or diagonal matrix, these conversions
        //                        may result in loss of data.
        // Note                 : Any pre-existing elements in mx will vanish
	// Note			: Used internal structure of other matrix types

// sosi - there must be a reason for resize use here, but damned if I recall what it is.

  {				
  switch (mx->stored_type())
    {
    case h_matrix_type:				// Matrix mx is Hermitian,
      (*(h_matrix *)mx) = (*this);		// no conversion necessary
      break;
    case n_matrix_type:				// Matrix mx is Normal, convert hmx to nmx
      {
      int hd = rows_;				// Dimension of hmx
      mx->resize(hd,hd);			// hmx into Normal mx
      complex *nij = ((n_matrix*)mx)->data;	// Start of data in nmx: <i|nmx|j>-><0|dmx|0>
      complex *hij = data;			// Start of data in hmx: <i|hmx|j>-><0|hmx|0>
      complex *hend = hij + (hd*hd+hd)/2;       // End of data in hmx: <rows_|hmx|cols_+1>
      complex *nend = nij + hd*hd;		// End of data in nmx: <rows_|hmx|cols_+1>
      complex *nji, *nii;			// Elements <j|nmx|i>, <i|nmx|i>
      for(nii=nij; hij<hend;nii+=hd+1)
        {
        *nii = *hij;				// <i|nmx|i> = <i|hmx|i>, i=j
        nji = nii + hd;				// <j|nmx|i> -> <i+1|nmx|i>
        nij = nii+1;				// <i|nmx|j> -> <i|nmx|i+1>
        hij++;					// <i|hmx|j> -> <i|hmx|i+1>
        for(; nji<nend; nij++,hij++,nji+=hd)
          {
          (*nij) = (*hij);			// <i|nmx|j> = <i|hmx|j>, i<j
          (*nji) = conj(*hij);			// <j|nmx|i> = <i|hmx*|j>, i<j
          }
        }
      }
      break;
    case d_matrix_type:				// Matrix mx is Diagonal, convert
      {						// hmx into Diagonal mx
      mx->resize(rows_,rows_);			// Resize, too bad if size mismatch!
      complex *dii = ((d_matrix*)mx)->data;	// Start of data in dmx: <i|dmx|i>-><0|dmx|0>
      complex *hii = data;			// Start of data in hmx: <i|hmx|i>-><0|hmx|0>
      complex *dend = dii + rows_;		// End of data in dmx: <rows_|dmx|rows+1>
      int hrow = rows_;				// Length of hmx row, start at row 0
      for(; dii<dend; dii++,hii+=hrow,hrow--)	// Copy <i|hmx|i> into <i|dmx|i>
        *dii = (*hii);				// & forget about hmx off-diagonals
      }
      break;
    case i_matrix_type:				// Matrix mx is Identity. Just use
      mx->resize(rows_, cols_);			// hmx size and none of its elements
      break;
    default:					// Matrix mx is of unknown type, so
      mx->resize(rows_,cols_);			// must copy over all of hmx elements
      complex *hij = data;			// Start of data in hmx: <i|hmx|j>-><0|hmx|0>
      int i, j=0;
      for(i=0; i<rows_; i++)
        {
	(*mx)(i,j) = (*hij);			// <i|mx|i> = <i|hmx|i>
        hij++;					// <i|hmx|i> -> <i|hmx|i+1>
	for(j=i+1; j<cols_; j++,hij++)
          {
	  (*mx)(i,j) = (*hij);			// <i|mx|j> = <i|hmx|j> i<j
	  (*mx)(j,i) = conj(*hij);		// <j|mx|i> = <i|hmx*|j> i<j
          }
        }
    }
  }


// ------------------------ New Conversion Routines ---------------------------
//                  (Been Having Troubles With The Old One)

/*      Input           hmx     : A Hermitian array (this)
        Output          mx      : A pointer to a newly created copy of
                                  hmx that has a differing type.
        Note                    : The conversion process MAY ALLOCATE
                                  new memory to hold the elements.
        Note                    : The conversion MAY LOOSE ELEMENTS
                                  if the new type doesn't hold all hmx's     */  
 


i_matrix* h_matrix::IMX()
  { 
  i_matrix* imx = new i_matrix(cols_, cols_);	// New i_matrix
  return imx;					// That's all we need
  }


d_matrix* h_matrix::DMX()
  {
  int id = rows_;				// Dimension of hmx
  d_matrix* dmx = new d_matrix(id,id);		// New empty d_matrix
  complex *hii = data;				// Start at <i|hmx|i>=<0|hmx|0>
  complex *dii = ((d_matrix*)dmx)->data;	// Start at <i|dmx|i>=<0|dmx|0>
  complex *hend = hii + (id*id+id)/2;		// Set hend as <id|hmx|id> + 1
  int hrow = id;				// The hmx row increment
  for(; hii<hend ;hii+=hrow,hrow--,dii++)	// Loop over hmx diagonals
    (*dii) = (*hii); 				// <i|dmx|i> = <i|hmx|i>
  return dmx;					// Return the diagonal array
  }

h_matrix* h_matrix::HMX() { return this; } 	// Just return myself

n_matrix* h_matrix::NMX()
  {
  int id = rows_;				// Dimension of hmx
  n_matrix* nmx = new n_matrix(id,id); 		// New empty n_matrix
  complex *nij = ((n_matrix*)nmx)->data;        // Start at <i|nmx|i>=<0|nmx|0>
  complex *hij = data;                          // Start at <i|hmx|i>=<0|hmx|0>
  complex *nji;					// Element nji=<j|nmx|i>
  int i,j;					// Added indexing variables
  for(i=0; i<id; i++, nij+=i) 			// Loop over matrix diagonals
    {
    (*nij) = (*hij);				// <i|nmx|j> = <i|hmx|j>, i=j
    nji = nij + cols_;				// <j|nmx|i> -> <i+1|nmx|i>
    nij++;					// <i|nmx|j> -> <i|nmx|i+1>
    hij++;					// <i|hmx|j> -> <i|hmx|i+1>
    for(j=i+1; j<cols_; j++,hij++,nij++)	// Loop columns from diagonal
      {
      (*nij) = (*hij);				// <i|nmx|j> = <i|hmx|j>, i<j
      (*nji) = conj(*hij);			// <j|nmx|i> = <i|hmx*|j>, i<j
      nji += cols_;				// <j|nmx|i> -> <j+1|nmx|i>
      }
    }
  return nmx;					// Now return the n_matrix
  }  

#endif						// h_matrix.cc
