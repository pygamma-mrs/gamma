/* n_matrix.cc **************************************************-*-c++-*-
**								   	**
**                             G A M M A		   		**
**									**
**	General Matrix	                          Implementation	**
**							   		**
**	Copyright (c) 1990, 1997			   		**
**	Tilo Levante, Scott A. Smith			  	 	**
**	Eidgenoessische Technische Hochschule		   		**
**	Labor fuer physikalische Chemie		   			**
**	8092 Zuerich / Switzerland		   			**
**						   			**
**      $Header: $
**						   			**
*************************************************************************/

/*************************************************************************
**									**
** Description								**
**									**
** The class n_matrix defines normal complex matrices for C++ with the	**
** usual operations +, -, *, / as well as in/output routines and the	**
** more commonly used matrix functions.					**
**									**
*************************************************************************/

#ifndef   Gn_matrix_cc_			// Is file already included?
#  define Gn_matrix_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <GamGen.h>			// Include OS specific stuff
#include <Matrix/n_matrix.h>		// Include the interface
#include <Matrix/d_matrix.h>		// Know about diagonal arrays
#include <Matrix/i_matrix.h>		// Know about identity arrays
#include <Matrix/h_matrix.h>		// Know about hermitian arrays
#include <Matrix/MxModBas.h>            // Include Matrix module errors 
#include <stdlib.h>
#include <fstream>                      // Include libstdc++ filestreams
#include <vector>			// Include libstdc++ STL vectors
#include <cmath>			// Inlcude HUGE_VAL_VAL


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
// i                  CLASS N_MATRIX CHECKING & ERROR HANDLING
// ____________________________________________________________________________

/*      Input               mx      : A normal complex matrix (this)
                            eidx    : Error index
                            noret   : Flag for linefeed (0=linefeed)
                            pname   : string in message
        Output              none    : Error message output
                                      Program execution stopped if fatal

  Note that these functions just route into the error functions in the base
  class _matrix.                                                             */

void n_matrix::NMxerror(int eidx, int nr) const
  {
  std::string CL="Complex Matrix";
  Mxerror(CL, eidx, nr);
  }

void n_matrix::NMxerror(int eidx, const std::string& PN, int nr) const
  {
  std::string CL="Complex Matrix";
  Mxerror(CL, eidx, PN, nr);
  }
 
volatile void n_matrix::NMxfatal(int eidx) const
  {
  std::string CL="Complex Matrix";
  Mxfatality(CL, eidx);
  }

volatile void n_matrix::NMxfatal(int eidx, const std::string& PN) const
  {
  std::string CL="Complex Matrix";
  Mxfatality(CL, eidx, PN);
  }

// ____________________________________________________________________________
// ii                        CLASS N_MATRIX CHECKING
// ____________________________________________________________________________
 
bool n_matrix::CheckDims(_matrix* mx, int warn)
  {
  if(cols_!=mx->rows() || cols_!=mx->cols())    // Insure mx dimension match
    {
    if(warn>0)
      {
      NMxerror(51,1);                           // Array dimensions mismatched
      NMxerror(31,1);                           // Row->Row Col-Col mismatch
      if(warn>1) NMxfatal(81);                  // Cannot continue on
      }
    return false;
    }
  return true;
  }
 
bool n_matrix::CheckDim(_matrix* mx, int mul, int warn)
  {
  bool TF;
  switch(mul)
    {
    default:
    case 0: TF = (cols_ == mx->cols()); break;
    case 1: TF = (cols_ == mx->rows()); break;
    case 2: TF = (rows_ == mx->rows()); break;
    case 3: TF = (rows_ == mx->cols()); break;
    }
  if(!TF)
    {
    if(warn>0)
      {
      NMxerror(51,1);                           // Array dimensions mismatched
      NMxerror(30,1);                           // Row->Col Col-Row mismatch
      if(warn>1) NMxfatal(81);                  // Cannot continue on
      }
    return false;
    }  
  return true;
  }  

// ____________________________________________________________________________
// A              CLASS N_MATRIX CONSTRUCTORS AND DESTRUCTOR
// ____________________________________________________________________________

/*  Arguments                              Constructed Array
    ---------           -------------------------------------------------------
        -                An empty array
     nr, nc              An (nr x nc) array, uninitialized
    nr, nc, z            An (nr x nc) array, sets <i|mx|i>=Re(z) & <i<j|mx|j>=z
       nmx               A duplicate of complex array nmx
 
   Note that all constructors invariably call the base class constructor
   (i.e. the one in _matrix) which sets the row & column sizes.  The added
   code should simply allocate the data array of proper size. Similarly, the
   destructor only destroys the data array, anything else is destroyed in
   the base class.                                                           */

n_matrix::n_matrix( )            : _matrix()
  {
  data = new complex[size]; 
  unitary = false;
  }

n_matrix::n_matrix(int i, int j) : _matrix(i,j)	
  {
  size = i*j;				// Compute the number of elements
  data = new complex[size];		// Allocate the data array
  unitary = false;
  }

n_matrix::n_matrix(int i, int j, const complex& z) : _matrix(i,j)			
  {
  size = i*j;				// Compute the number of elements
  data = new complex[size];		// Allocate the data array
  for(int pos=0; pos<size; pos ++)	// Set all elements to value z
    data[pos] = z;
  unitary = false;
  }

n_matrix::n_matrix(const n_matrix& mx) : _matrix(mx)
  {
  size = mx.size;			// Copy the matrix size
  data = new complex[size];		// Construct a new data array
  for(int pos=0; pos<size; pos++)	// Copy all the data points
    data[pos] = mx.data[pos];
  unitary = mx.unitary;
  }

n_matrix::~n_matrix () { delete [] data; }

// ____________________________________________________________________________
// B                CLASS N_MATRIX ACCESS AND ASSIGNMENT
// ____________________________________________________________________________

/* There are two flavors of element access operators herein. The one that might
   commonly be used, i.e. mx(i,j), is DIRECTLY accessing the element <i|mx|j>.
   That can save CPU time by avoiding value copying when altering or obtaining
   elements.  Its misuse can lead to trouble in some matrix classes because not
   all matrix elements are stored.  However, for normal arrays all elements are
   stored so there are no  restrictions on elements or element access.

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
    (int,int)           Reference to element <i|nmx|j> (Potential Danger)
    get(int,int)        Copy of element <i|nmx|j>      (Safe)
    put(int,int)        Assigns element <i|nmx|j>      (Return FALSE if fails)
    put_h(int,int)      Assigns both <i|nmx|j> & <j|nmx|i> (Return TRUE/FALSE)
    get_block(r,c,R,C)  Returns nmx block of size RxC starting at <r|nmx|c>
    put_block(r,c,mx)   Places mx into nmx at position <r|hmx|c> (TRUE/FALSE)

    Again, note that operations which return TF will trigger class matrix to  
    take appropriate action when they return FALSE.  This usually means that
    the operation cannot be handled by a normal array, but that is unlikely.
    Were it true, the operation is reattempted.                             */

void n_matrix::operator= (const n_matrix& mx)
  {
  if(this == &mx) return; 			// Do nothing on self-assign
  delete [] data;				// First delete the data
  _matrix::operator= (mx);			// Assign by generic function
  data = new complex[size];			// Construct new data array
  for(int pos=0; pos<size; pos ++)		// Copy the data array
    data[pos] = mx.data[pos];
  unitary = mx.unitary;				// Copy our unitary flag
  }

complex& n_matrix::operator() (int i, int j)  { return data[i*cols_+j]; }
complex  n_matrix::get(int i, int j) const    { return data[i*cols_+j]; }

bool n_matrix::put(const complex& z, int i, int j)
  {
  data[i*cols()+j] = z;				// Set the matrix element to z
  unitary = false;				// Breaks up unitary
  return true;					// Flag that this worked o.k.
  }

bool n_matrix::put_h(const complex& z, int i, int j) 
  {
  data[j*cols()+i] = conj(z);			// Set transpose element to z*
  data[i*cols()+j] = z;				// Set the element to z
  unitary = false;				// Breaks up unitary
  return true;					// Flag that this worked o.k.
  }

_matrix* n_matrix::get_block(int row, int col, int nrows, int ncols)
  {
  if(!row && !col)				// If block start @ <0|nmx|0> &
    if((nrows==rows()) && (ncols==cols()))	// spans all of nmx, return nmx
      return this;
  n_matrix* submx = new n_matrix(nrows, ncols);	// Constuct a new n_matrix
  int src = row*cols() + col;			// Pos.1st blk emt <row|nmx|col>
  complex* nIJ = &data[src];			// 1st blk emt:nIJ=<row|nmx|col>
  complex* mij = &(submx->data[0]);		// 1st ret emt:mij=<0|submx|0>
  int step = cols()-ncols; 			// Amt. to skip each row of nmx
  for(int i=0; i<nrows; i++, nIJ+=step)		// Loop all nrows of block
    for(int j=0; j<ncols; j++, nIJ++, mij++)	// Loop all ncols of block
      *mij = *nIJ;				// Set <i|submx|j> = <I|nmx|J>
  return submx;
  }

bool n_matrix::put_block(int row, int col, _matrix* mx)
  {
  int ncols = mx->cols();				// Get the number of columns to change
  int nrows = mx->rows();				// Get the number of rows to change
  int c = cols();					// Get the current number of columns
  if((row+nrows > rows())||(col+ncols > c))		// Check that put does not go beyond
    {
    NMxerror(20,1);					//    Block put exceeds dimensions
    NMxfatal(5,"put_block");				//    Bad use of put_block function
    return false;
    }
  else
    {
    int i, j;
    complex* nIJ = &data[row*c + col];			// Element of nmx: <I|nmx|J> -> <row|nmx|col>
    int step = c - ncols;	 			// Amount to skip nmx, each row
    switch(mx->stored_type())
      {
      case n_matrix_type:				// Put nmx into nmx
	{
	complex* mij = ((n_matrix*)mx)->data;	 	// Element of input mx: <i|mx|j> -> <0|mx|0>
	for(i=0; i<nrows; i++, nIJ+=step)	 	// Loop through all nrows to be modified
	  for(j=0; j<ncols; j++, mij++, nIJ++)		// Loop through all ncols to be modified
	    *nIJ = *mij;				// <I|nmx|J> = <i+row|nmx|j+col> = <i|mx|j>
	}
	break;
      case h_matrix_type:				// Put hmx into nmx
        {
        complex* nJI;					// For element of nmx of hmx adjoint element
	complex* mij = ((h_matrix*)mx)->data;	 	// Element of input hmx: <i|hmx|j> -> <0|hmx|0>
	for(i=0; i<nrows; i++, nIJ+=(step+i))	 	// Loop through all nrows to be modified
          {						// First set the nmx element of <i|hmx|i> 
	  *nIJ = *mij;					// <I|nmx|J> = <i+row|nmx|i+col> = <i|hmx|i>
          nJI = nIJ + c; 				// Set <J|nmx|I> = <i+1+row|i+col> 
          nIJ++;					// Set <I|nmx|J> = <i+row|nmx|i+1+col>
          mij++;					// Increment <i|hmx|i> to <i|hmx|i+1>
	  for(j=i+1; j<ncols; j++,mij++,nIJ++,nJI+=c)	// Loop through all ncols to be modified
            { 						// Set elements for <i|hmx|j> & <j|hmx|i> 
	    *nIJ = *mij;				// <I|nmx|J> = <i+row|nmx|j+col> = <i|hmx|j>
            *nJI = conj(*mij);				// <J'|nmx|I'> = <j+row|nmx|i+col> = <j|hmx|i>
	    }
	  }
        }
	break;
      case d_matrix_type:				// Put dmx into nmx
	{
	complex *dii = ((d_matrix*)mx)->data;		// Element of dmx: <i|dmx|i> -> <0|dmx|0>
	for(i=0; i<nrows; i++, nIJ+=step, dii++)	// Loop through all nrows to be modified
          {
	  for(j=0; j<ncols; j++, nIJ++)			// Loop through all ncols to be modified
	    *nIJ = 0; 					// <I|nmx|J> = <i+row|nmx|j+col> = <i|mx|j> = 0
	  data[c*(i+row)+i+col] = *dii; 		// <I|nmx|J> = <i+row|nmx|i+col> = <i|mx|i>
          }
	}						// use dmx access so don't rely on its structure
	break;
      case i_matrix_type:				// Put imx into nmx
	{
	for(i=0; i<nrows; i++, nIJ+=step) 		// Loop through all nrows to be modified
          {
	  for(j=0; j<ncols; j++, nIJ++)			// Loop through all ncols to be modified
	    *nIJ = 0; 					// <I|nmx|J> = <i+row|nmx|j+col> = <i|mx|j> = 0
	  data[c*(i+row)+i+col] = 1; 			// <I|nmx|J> = <i+row|nmx|i+col> = <i|mx|i> = 1
          }
	}
	break;
      default:						// Put unknown matrix into nmx
	{
	for(i=0; i<nrows; i++, nIJ+=step) 		// Loop through all nrows to be modified
	  for(j=0; j<ncols; j++, nIJ++)			// Loop through all ncols to be modified
	    *nIJ = (*mx)(i,j); 				// <I|nmx|J> = <i+row|nmx|j+col> = <i|mx|j>
	}						// Used mx access function since unknown type! 
      }
    unitary = false;					// Breaks up unitary
    return true;
    }
  }

// ____________________________________________________________________________
// C                 CLASS N_MATRIX HERMITIAN_TYPE HANDLING
// ____________________________________________________________________________

/* These functions handle the checking of whether an array is of a Hermitian
   type.  Since GAMMA has a Hermitian matrix type (h_matrix) it is easy to
   know the whether the array is stored Hermitian or not. The other funtion
   will test if the array is Hermitian to with d, where the test performed is
   norm(<i|mx|j> - <j|mx|i>*) < d for all elements. Of couse the array must
   also be square and the diagonals real values. These functions are a bit
   outdated now that GAMMA has a Hermitian array class and will probably be
   removed in favor of a TF is_hermitian function function or something.  More
   appropriate would be an stored_unitary function.....
 
     Function            Output                       Description
   ----------------  --------------  ------------------------------------------
   stored_hermitian  hermitian_type  Returns _hermitian=FALSE
   test_hermitian    hermitian_type  Returns _hermitian=T/F (within d)       */
 
// Note:  Functions set_hermitian & check_hermitian are found in class matrix!

hermitian_type n_matrix::stored_hermitian( ) const { return non_hermitian; }
hermitian_type n_matrix::test_hermitian(double d) const
  {
  int i, j, posji, posij=0;
  int f = (rows_==cols_);		// Insure matrix is square
  for(i=0; (i<rows_)&&f; i++,posij+=i)	// Then loop over the upper triangle
    {					// Start with nmx[posij] = <i|nmx|i>
    f = (fabs(Im(data[posij])) < d);	// Make sure the diagonal is real
    posji = posij + cols_;		// Set position nmx[posji] = <i+1|nmx|i>
    posij++;				// Set position nmx[posij] = <i|nmx|i+1>
    for(j=i+1; (j<rows_)&&f;		// Loop over upper to lower elements
             j++,posij++,posji+=cols_)	// nmx[posij] = <i|nmx|j> j>i
      f = (norm(data[posij]		// nmx[posji] = <j|nmx|i> i<j
              -conj(data[posji])) < d);
    }
  return f?_hermitian:non_hermitian;
  }

// ____________________________________________________________________________
// D                  CLASS N_MATRIX MATRIX_TYPE HANDLING
// ____________________________________________________________________________

/* These functions handle the checking of what type of array we have.  These
   types directly correspond to the matrix classes (such as n_matrix) derived
   from class _matrix.
 
     Function            Output                       Description
   ----------------  --------------  ------------------------------------------
   set_type            --------      Found in class matrix
   test_type           --------      Found in class matrix
   stored_type       n_matrix_type   Always returns we are normal complex
   test_hermitian    *_matrix_type   Type nmx could be within d
   mxtype            string          Returns the string "Complex"

   The test type looks to see if nmx could be an array of the matrix type
   specified as an input argument.  If it can within "d" the return type
   is equal to the input type.  If it cannot the return is n_matrix_type     */  
 
matrix_type n_matrix::stored_type( )  const { return n_matrix_type; }
std::string      n_matrix::mxtype()        const { return std::string("Complex"); }
std::string      n_matrix::mxtype(bool pf) const 
  {
  if(!pf) return std::string("Complex");
  else if(is_real())      return std::string("Real");
  else if(is_imaginary()) return std::string("Imaginary");
  return std::string("Complex");
  }

matrix_type n_matrix::test_type(matrix_type m, double d) const

        // Input                nmx  : Matrix (this)
        //                      t    : Matrix type
        //                      d    : Double precision number
        // Output               type : Tests nmx to see if it could be stored
        //                             as type t within the precison d.  It
        //                             if can, type t is returned. If it can't
        //                             the stored type is returned.
	// Note			     : Default value of d is GMxCut

  {
  int f;					// Test flag: True if conversion possible
  int i,j,pos;					// Column and row counters
  switch(m)
    {
    case h_matrix_type :			// See if n_matrix could be h_matrix
      f = (_hermitian == test_hermitian(d));	// Just use test_hermitian function
      return f?h_matrix_type:n_matrix_type;	// Return h_matrix or n_matrix
      break;
    case d_matrix_type :			// See if n_matrix could be d_matrix
      f = (rows_ == cols_);			// 1. Matrix must be square
      for(i=0, pos=0; (i<rows_) && f; i++)	// 2. All off diagonals must have a
	for(j=0; (j<cols_) && f; j++,pos++)	//    magnitude < d
	  if(i != j)
	    f = (norm(data[pos]) < d);
      return f?d_matrix_type:n_matrix_type;	// Return d_matrix or n_matrix
      break;
    case i_matrix_type :			// See if n_matrix could be i_matrix
      f = (rows_ == cols_);			// 1. Matrix must be square
      for(i=0, pos=0; (i<rows_) && f; i++)	// 2. All off diagonals must have a
	for(j=0; (j<cols_) && f; j++, pos++)	//    magnitude < d
	  if(i != j)				// 3. Add diagonals must have 
	    f = (norm(data[pos]) < d);		//    (magnitude - 1) < d
	  else
	    f = (norm(data[pos]-complex1) < d);
      return f?i_matrix_type:n_matrix_type;	// Return i_matrix or n_matrix
      break;
    default:					// See if any matrix can be stored
      return n_matrix_type;			// as n_matrix. Of course, always true
      break;
    }
  return n_matrix_type;
  }


// ____________________________________________________________________________
// E                  CLASS N_MATRIX VARIOUS MATRIX CHECKS
// ____________________________________________________________________________

/*   Function    Output                       Description
   ------------  ------  ------------------------------------------------------
   is_symmetric   bool   TF if <i|nmx|j>-<j|nmx*|i> < d,     (def. d GMxCut)
   is_hermitian   bool   TF if <i|nmx|j>-<j|nmx*|i> < d,     (def. d GMxCut)
   is_unitary     bool   TF if inv(nmx) == adjoint(nmx), CPU intensive
   is_real        bool   TF if Im(<i|nmx|j>) < d for all i,j (def. d GMxCut)
   is_imaginary   bool   TF if Re(<i|nmx|j>) < d for all i,j (def. d GMxCut)
   is_complex     bool   TF if is_real && is_imaginary       (def. d GMxCut)
   is_zero        bool   TF if ||<i|nmx|j>|| < d for all i,j (def. d GMxCut)
   is_diagonal    bool   TF if ||<i|nmx|j>|| < d for all i!=j(def. d GMxCut)
   is_square      bool   TF if rows_ == cols_
   is_equal       bool   TF if ||<i|nmx-mx|j>||<d    all i,j (def. d GMxCut)

The unitary check herein tests if the array inverse is equal to its adjoint.
The routine multiplies the matrix by its adjoint and then looks to see how
well the result conforms to an identity matrix. This can use lots of time....
                            ---
                            \                   *t 
                   z(i,j) = /    <i|nmx|k><k|nmx  |j> = del
                            ---                            ij
                             k                                              */

bool n_matrix::is_symmetric(double d) const
  {
  if(cols_ != rows_) return false;
  bool flag=true;
  int r, c;
  for(r=rows_-1; (r>=0)&&flag; r--)
    for(c=cols_-1; (c>=0)&&flag; c--)
      flag=(fabs(norm((*this).get(r,c)-(*this).get(c,r))) < d);
  return flag;
  }

bool n_matrix::is_hermitian(double d) const
  {
  if(cols_ != rows_) return false;
  bool flag=true;
  int r, c;
  for(r=rows_-1; (r>=0)&&flag; r--)
    for(c=cols_-1; (c>=0)&&flag; c--)
      flag=(fabs(norm((*this).get(r,c)-conj((*this).get(c,r)))) < d);
  return flag;
  }

bool n_matrix::is_unitary(double d) const
  {
  if(unitary) return unitary;
  bool f = true;
  int nr = rows_;
  f = (rows_ == cols_);				// Unitary must be square
  complex z;
  int i,j,k,pi0,pik,pjk;
  for(i=0,pi0=0;  i<nr && f; i++,pi0+=rows_)	// Sum over all the rows
    {
    for(j=0,pjk=0; j<nr && f; j++)		// Sum over all of the columns
      {
      z = 0;					// Set <i|pdt|j> = 0
      for(k=0,pik=pi0; k<nr; k++, pik++,pjk++)	// Sum over k, do the addition
        z += data[pik]*conj(data[pjk]);		// <i|pdt|j>+=<i|mx|k><j|mx*|k>
      if(i == j)				// Diagonal must be 1
        {
        f = ((fabs(Re(z)-1)) < d);
        if(f) f = (fabs(Im(z)) < d);
        }
      else f = (norm(z) < d); 			// Off diagonal must be zero
      }
    }
  return f;
  }

bool n_matrix::is_real(double d) const
  {
  bool f=1;
  for(int pos=0; pos<size && f; pos++)		// Loop over all the elements
    f = (fabs(Im(data[pos])) <= d);		// Stop if any element has an
  return f;					// imaginary component
  }

bool n_matrix::is_imaginary(double d) const
  {
  bool f=true;
  for(int pos=0; pos<size && f; pos++)		// Loop over all the elements
    f = (fabs(Re(data[pos])) <= d);		// Stop if any real component
  if(f)	f = (!is_zero(d)); 			// Insure not the zero mx
  return f;
  }

bool n_matrix::is_complex(double d) const
  {
  if(is_imaginary(d)) return false; 		// False if pure imaginary
  else if(is_real(d)) return false; 		// False if pure real
  return true;
  }

bool n_matrix::is_zero(double d) const
  {
  bool f=true;					// Assume we are zero
  for(int pos=0; pos<size && f; pos++)		// Loop over all the elements
    f = (norm(data[pos]) <= d);			// Stop if any non-zero elemt
  return f;
  }

bool n_matrix::is_diagonal(double d) const
  {
  if(!is_square()) return false;		// Not diagonal if not square
  bool f=true;					// Assume we are diagonal
  int cd = cols_;				// This is number of columns
  complex *ni0 = data;				// ni0 -> <0|nmx|0>
  complex *nend = ni0 + size;			// nend = <rows_-1|nmx|cols_>
  complex *nii;					// nii = <i|nmx|i>
  complex *nij;					// nij = <i|nmx|j>
  for(nii=ni0; nii<nend && f; ni0+=cd,nii+=cd+1)// Loop rows of nmx
    {
    for(nij=ni0; nij<nii && f; nij++)		// Check lower triangle
      f = (norm(*nij) < d);			// elements for non-zero
    for(nij++; nij<ni0+cd && f; nij++)		// Check upper triangle
      f = (norm(*nij) < d);			// elements for non-zero
    }
  return f;
  }

bool n_matrix::is_square() const { return (rows_ == cols_); }
 
bool n_matrix::is_equal(_matrix* mx, double d) const
  {
  if(cols_ != mx->cols())                                return false;
  if(rows_ != mx->rows())                                return false;
   bool flag=true;
   int r, c;
   for(r=rows_-1; (r>=0)&&flag; r--)
     for(c=cols_-1; (c>=0)&&flag; c--)
       flag=((*this).get(r,c)==(*mx)(r,c));
   return flag;
   }

// ____________________________________________________________________________
// F              CLASS N_MATRIX BINARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________

_matrix* n_matrix::add(_matrix* mx)

        // Input                nmx  : Input normal matrix (this)
	//			 mx  : Second matrix
        // Output               mx1  : New matrix which is the addition of
	//			       the two input arrays; mx1 = nmx + mx
	// Note			     : Uses internal structure of other matrix types

  {
  if((rows_!=mx->rows()) || (cols_!=mx->cols()))// Insure sizes compatible 
    {
    NMxerror(51,1);				//    Mismatched dimensions
    NMxfatal(20);				//    Cannot perform addition
    return mx;
    }
  else
    switch(mx->stored_type())
      {
      case n_matrix_type:				// Add two n_matrices
	{ 
	n_matrix* sum = new n_matrix(rows_,cols_);	// Construct new empty matrix sum
	complex *mij = &((n_matrix*)mx)->data[size-1];	// Element of mx:  mij = <i|mx|j>  -> <rows_|mx|cols_>
	complex *nij = &data[size-1];		    	// Element of nmx: nij = <i|nmx|j> -> <rows_|nmx|cols_>
	complex *sij = &(sum->data[size-1]);	        // Element of sum: sij = <i|sum|j> -> <rows_|sum|cols_>
	for(; nij>=data; sij--,nij--,mij--)		// Perform the addition element-wise
	  (*sij) = (*mij)+(*nij);			// <i|sum|j> = <i|nmx|j> + <i|mx|j>
	return sum;
	}
	break;
      case h_matrix_type:				// Add h_matrix to n_matrix
	{
  	n_matrix* sum = new n_matrix(rows_,cols_);	// Construct new empty matrix sum
        complex *hij = ((h_matrix*)mx)->data;		// Element of hmx: hij = <i|hmx|j> -> <0|hmx|0>
        complex *nij = data;				// Element of nmx: nij = <i|nmx|j> -> <0|nmx|0>
        complex *mij = ((n_matrix*)sum)->data;		// Element of sum: mij = <i|sum|j> -> <0|sum|0>
        complex *nji, *mji;				// Elements nji=<j|nmx|i>, mji=<j|sum|i>
        int i,j; 
        for(i=0; i<rows_; i++, nij+=i, mij+=i)		// Loop over the rows
          {
    	  (*mij) = (*nij) + (*hij);			// <i|sum|j> = <i|nmx|j> + <i|hmx|j>, i=j
          mji = mij + cols_;				// <j|sum|i> -> <i+1|sum|i>
          nji = nij + cols_;				// <j|nmx|i> -> <i+1|nmx|i>
          mij++;					// <i|sum|j> -> <i|sum|i+1>
          nij++;					// <i|nmx|j> -> <i|nmx|i+1>
          hij++;					// <i|hmx|j> -> <i|hmx|i+1>
  	  for(j=i+1; j<cols_; j++,hij++,mij++,nij++)	// Loop over the columns
            {
    	    (*mij) = (*nij) + (*hij);			// <i|sum|j> = <i|nmx|j> + <i|hmx|j>, i<j
  	    (*mji) = (*nji) + conj(*hij);		// <j|sum|i> = <j|nmx|i> + <i|hmx*|j>, i<j
            mji += cols_;				// <j|sum|i> -> <j+1|sum|i>
            nji += cols_;				// <j|nmx|i> -> <j+1|nmx|i>
            }
          } 
	return sum;
	}
	break;
      case d_matrix_type:				// Add d_matrix to n_matrix
	{
	n_matrix* sum = new n_matrix(*this);		// New n_matrix equal to this
	for(int i=0; i<rows_; i++)			// Deal with only the diagonals
	  sum->data[i*rows_+i] +=			// <i|sum|j> = <i|nmx|j>+del   <i|dmx|j>
                             ((d_matrix*)mx)->data[i];	//                          i,j
        sum->unitary = false;
	return sum;
	}
	break;
      case i_matrix_type:				// Add i_matrix to n_matrix
	{
	n_matrix* sum = new n_matrix(*this);		// New n_matrix equal to this
	for(int i=0; i<rows_; i++)			// Add 1 to its diagonal
	  sum->data[i*rows_+i] += complex1;		// <i|sum|j> = <i|nmx|j>+del   <i|imx|j>
        sum->unitary = false;
	return sum;					//                          i,j
	}
	break;
      default:						// Add generic matrix to n_matrix
	{ 
	n_matrix* sum = new n_matrix(rows_,cols_);
	int pos = 0;
	for(int i=0; i<rows_; i++)
	  for(int j=0; j<cols_; j++,pos++)
	    sum->data[pos] = data[pos]+(*mx)(i,j);
	return sum;
	}
      }
    }

_matrix* n_matrix::subtract(_matrix* mx)

        // Input                nmx  : Input normal matrix (this)
	//			 mx  : Second matrix
        // Output               sub  : New matrix which is the subtraction of
	//			       the two input arrays; sub = nmx - mx
	// Note			     : Uses internal structure of other matrix types

  {
  if((rows_!=mx->rows()) || (cols_!=mx->cols()))
    {
    NMxerror(51,1);				//    Mismatched dimensions
    NMxfatal(21);				//    Cannot perform subtraction
    return mx;
    }
  else
    switch(mx->stored_type())
      {
      case n_matrix_type:				// Subtract n_matrix from n_matrix
	{
	n_matrix* sub = new n_matrix(rows_,cols_);	// New empty n_matrix
	complex *mij = &((n_matrix*)mx)->data[size-1];	// Pointer to last mx element
	complex *nij = &data[size-1];		    	// Pointer to last nmx element
	complex *sij = &(sub->data[size-1]);	        // pointer to last sub element
	for(; nij>=data; sij--,nij--,mij--)		// Perform subtraction element-wise
	  (*sij) = (*nij)-(*mij);
	return sub;
	}
	break;
      case h_matrix_type:				// Add h_matrix to n_matrix
	{
  	n_matrix* sub = new n_matrix(rows_,cols_);	// Construct new empty matrix sub
        complex *hij = ((h_matrix*)mx)->data;		// Element of hmx: hij = <i|hmx|j> -> <0|hmx|0>
        complex *nij = data;				// Element of nmx: nij = <i|nmx|j> -> <0|nmx|0>
        complex *mij = ((n_matrix*)sub)->data;		// Element of sub: mij = <i|sub|j> -> <0|sub|0>
        complex *nji, *mji;				// Elements nji=<j|nmx|i>, mji=<j|sub|i>
        int i,j; 
        for(i=0; i<rows_; i++, nij+=i, mij+=i)		// Loop over the rows
          {
    	  (*mij) = (*nij) - (*hij);			// <i|sub|j> = <i|nmx|j> - <i|hmx|j>, i=j
          mji = mij + cols_;				// <j|sub|i> -> <i+1|sub|i>
          nji = nij + cols_;				// <j|nmx|i> -> <i+1|nmx|i>
          mij++;					// <i|sub|j> -> <i|sub|i+1>
          nij++;					// <i|nmx|j> -> <i|nmx|i+1>
          hij++;					// <i|hmx|j> -> <i|hmx|i+1>
  	  for(j=i+1; j<cols_; j++,hij++,mij++,nij++)	// Loop over the columns
            {
    	    (*mij) = (*nij) - (*hij);			// <i|sub|j> = <i|nmx|j> - <i|hmx|j>, i<j
  	    (*mji) = (*nji) - conj(*hij);		// <j|sub|i> = <j|nmx|i> - <i|hmx*|j>, i<j
            mji += cols_;				// <j|sub|i> -> <j+1|sub|i>
            nji += cols_;				// <j|nmx|i> -> <j+1|nmx|i>
            }
          } 
	return sub;
	}
	break;
      case d_matrix_type:				// Subtract d_matrix from n_matrix
	{
	n_matrix* sub = new n_matrix(*this);		// New n_matrix equal to this
	for(int i=0; i<rows_; i++)			// Deal with only the diagonals
	  sub->data[i*rows_+i] -=			// <i|sub|j> = <i|nmx|j>-del   <i|dmx|j>
                             ((d_matrix*)mx)->data[i];	//                          i,j
        sub->unitary = false;
	return sub;
	}
	break;
      case i_matrix_type:				// Subtract i_matrix from n_matrix
	{
	n_matrix* sub = new n_matrix(*this);		// New matrix equal to this
	for(int i=0; i<rows_; i++)			// Subtract 1 from its diagonal
	  sub->data[i*rows_+i] -= complex1;		// <i|sub|j> = <i|nmx|j>+del   <i|imx|j>
        sub->unitary = false;
	return sub;
	}
	break;
      default:						// Subtract generic matrix from n_matrix
	{ 
	n_matrix* sub = new n_matrix(rows_,cols_);
	int pos = 0;
	for(int i=0; i<rows_; i++)
	  for(int j=0; j<cols_; j++,pos++)
	    sub->data[pos] = data[pos]-(*mx)(i,j);
	return sub;
	}
      }
    }


_matrix* n_matrix::multiply(_matrix* mx)

// Input                nmx  : Input normal matrix (this)
//			 mx  : Second matrix
// Output               pdt  : New matrix which is the product nmx*mx
// Note                      : Uses internal structure of other matrix types

//
//                            ---
//                <i|pdt|j> = \   <i|nmx|k> <k|mx|j>
//                            /
//                            ---
//		               k
{
  if(cols() != mx->rows())				// Insure nmx & mx dimension match
  {
    NMxerror(51,1);					//    Mismatched dimensions
    NMxfatal(22);					//    Cannot perform multiplication
    return mx;
  }
  else
  {
    switch(mx->stored_type())
    {
      case n_matrix_type:				// Multiply n_matrix into n_matrix
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
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A_rows, B_cols, A_cols, 
                             alpha, (double *)data, A_cols, (double *)((n_matrix*)mx)->data, B_cols, beta, (double *)pdt->data, C_cols);
//              std::cerr << "BLAS: n_matrix * n_matrix\n";
            }
            else
            {
                int r = rows();					// Rows of product matrix (& of nmx)
                int c = mx->cols();				// Columns of product matrix (& of mx)
                int l = cols();					// Inner loop dimension (nmx cols & mx rows)
                complex *p00 = pdt->data;			// Start of data in pdt matrix
                complex *n00 = data;				// Start of data in nmx matrix
                complex *m00 = ((n_matrix*)mx)->data;	 	// Start of data in mx matrix
                complex *pend = p00 + r*c;			// End of data in pdt matrix
                complex *m0c = m00 + c;				// End of first row data in mx matrix
                complex *mend = m00 + c*l;			// End of data in mx matrix
                complex *nik = n00;				// Index for <i|nmx|k> set to <0|nmx|0>
                complex *pij, *ni0, *m0j, *mkj;
                complex z;					// Intermediate storage for <i|pdt|j>

                for(pij=p00, ni0=n00; pij<pend; ni0=nik)	// Effective sum over index i
                {
                  for(m0j=m00; m0j<m0c; m0j++, pij++)		// Effective sum over index j
                  {						// Effective sum over index k 
                    for(z=0,nik=ni0,mkj=m0j; mkj<mend; nik++,mkj+=c)
                      z += (*nik) * (*mkj);			// pij += <i|nmx|k><k|mx|j>
                      
                    *pij = z;
                  }
                }
//              std::cerr << "NON BLAS: n_matrix * n_matrix\n";

            }
#else
        int r = rows();					// Rows of product matrix (& of nmx)
        int c = mx->cols();				// Columns of product matrix (& of mx)
        int l = cols();					// Inner loop dimension (nmx cols & mx rows)
        n_matrix* pdt =	new n_matrix(r,c);		// Create new matrix for product
        complex *p00 = pdt->data;			// Start of data in pdt matrix
        complex *n00 = data;				// Start of data in nmx matrix
        complex *m00 = ((n_matrix*)mx)->data;	 	// Start of data in mx matrix
        complex *pend = p00 + r*c;			// End of data in pdt matrix
        complex *m0c = m00 + c;				// End of first row data in mx matrix
        complex *mend = m00 + c*l;			// End of data in mx matrix
        complex *nik = n00;				// Index for <i|nmx|k> set to <0|nmx|0>
        complex *pij, *ni0, *m0j, *mkj;
        complex z;					// Intermediate storage for <i|pdt|j>

        for(pij=p00, ni0=n00; pij<pend; ni0=nik)	// Effective sum over index i
        {
          for(m0j=m00; m0j<m0c; m0j++, pij++)		// Effective sum over index j
          {						// Effective sum over index k 
            for(z=0,nik=ni0,mkj=m0j; mkj<mend; nik++,mkj+=c)
              z += (*nik) * (*mkj);			// pij += <i|nmx|k><k|mx|j>
              
            *pij = z;
          }
        }
#endif
        return pdt;
        break;
      }
      
      case h_matrix_type:				// Multiply n_matrix into h_matrix
      {
#ifdef _USING_BLAS_
         int A_rows = rows();
         int A_cols = cols();            
         int B_rows = mx->rows();
         int B_cols = mx->cols();
         int C_rows = A_rows;
         int C_cols = B_cols;

         int r = rows();					// Rows of product matrix (& of nmx)
         int c = mx->cols();				// Columns of product matrix (& of mx)
         n_matrix* pdt =	new n_matrix(C_rows, C_cols);

         if(C_rows * C_cols > 16)  //4*4 = 16
         {
            double alpha[2] = {0,0};
            double beta[2] = {0,0};
            alpha[0] = 1.0;
            beta[0]  = 0.0;
            n_matrix* hmx =	new n_matrix(B_rows,B_cols);		// Create new matrix h_matrix
            mx->convert(hmx);				// Convert h_matrix mx into normal matrix hmx
	    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A_rows, B_cols, A_cols, 
                             alpha, (double *)data, A_cols, (double *)hmx->data, B_cols, beta, (double *)pdt->data, C_cols);
//          std::cerr << "BLAS: n_matrix * h_matrix\n";
	    delete hmx;
         }
         else
         { int l = cols();					// Inner loop dimension (nmx cols & mx rows)
           n_matrix* hmx =	new n_matrix(c,c);		// Create new matrix h_matrix
           mx->convert(hmx);				// Convert h_matrix mx into normal matrix hmx
           complex *p00 = pdt->data;			// Start of data in pdt matrix
           complex *n00 = data;				// Start of data in nmx matrix
           // sosi- looks like should be h_matrix below, fix & retest
           complex *h00 = ((n_matrix*)hmx)->data; 		// Start of data in hmx matrix
           complex *pend = p00 + r*c;			// End of data in pdt matrix
           complex *h0c = h00 + c;				// End of first row data in mx matrix
           // sosi- looks like next line is wrong also
           complex *hend = h00 + c*l;			// End of data in mx matrix
           complex *nik = n00;				// Index for <i|nmx|k> set to <0|nmx|0>
           complex *pij, *ni0, *h0j, *hkj;
           complex z;					// Intermediate storage for <i|pdt|j>
           for(pij=p00, ni0=n00; pij<pend; ni0=nik)	// Effective sum over index i
           {
             for(h0j=h00; h0j<h0c; h0j++, pij++)		// Effective sum over index j
             {						// Effective sum over index k 
              for(z=0,nik=ni0,hkj=h0j; hkj<hend; nik++,hkj+=c)
                 z += (*nik) * (*hkj);			// pij += <i|nmx|k><k|hmx|j>
               // sosi- looks incorrect? How does this work at all?
               *pij = z;
               }
           }
           delete hmx;					// Get rid of temporary hmx as n_matrix
//         std::cerr << "NON BLAS: n_matrix * h_matrix\n";
	 }
#else
        int r = rows();					// Rows of product matrix (& of nmx)
        int c = mx->cols();				// Columns of product matrix (& of mx)
        n_matrix* hmx =	new n_matrix(c,c);		// Create new matrix h_matrix
        mx->convert(hmx);				// Convert h_matrix mx into normal matrix hmx
        int l = cols();					// Inner loop dimension (nmx cols & mx rows)
        n_matrix* pdt =	new n_matrix(r,c);		// Create new matrix for product
        complex *p00 = pdt->data;			// Start of data in pdt matrix
        complex *n00 = data;				// Start of data in nmx matrix
        // sosi- looks like should be h_matrix below, fix & retest
        complex *h00 = ((n_matrix*)hmx)->data; 		// Start of data in hmx matrix
        complex *pend = p00 + r*c;			// End of data in pdt matrix
        complex *h0c = h00 + c;				// End of first row data in mx matrix
        // sosi- looks like next line is wrong also
        complex *hend = h00 + c*l;			// End of data in mx matrix
        complex *nik = n00;				// Index for <i|nmx|k> set to <0|nmx|0>
        complex *pij, *ni0, *h0j, *hkj;
        complex z;					// Intermediate storage for <i|pdt|j>
        for(pij=p00, ni0=n00; pij<pend; ni0=nik)	// Effective sum over index i
        {
          for(h0j=h00; h0j<h0c; h0j++, pij++)		// Effective sum over index j
            {						// Effective sum over index k 
            for(z=0,nik=ni0,hkj=h0j; hkj<hend; nik++,hkj+=c)
              z += (*nik) * (*hkj);			// pij += <i|nmx|k><k|hmx|j>
            // sosi- looks incorrect? How does this work at all?
            *pij = z;
            }
        }
        delete hmx;	  				// Get rid of temporary hmx as n_matrix
#endif
        return pdt;
        break;
	  }
      
      case d_matrix_type:				// Multiply n_matrix into d_matrix
	  {
        int c = cols();					// Rows of product matrix
        int r = rows();					// Cols of product matrix
        n_matrix* pdt = new n_matrix(r,c); 		// Create new normal matrix
        complex *p00 = pdt->data; 			// Start of data of pdt
        complex *n00 = data;				// Start of data of nmx
        complex *d00 = ((d_matrix*)mx)->data;		// Start of data of dmx
        complex *pend = p00 + r*c;			// End of data of pdt
        complex *dend = d00 + c;			// End of data in dmx
        complex *pij, *nij, *dij, *ni0;			// Matrix elements & <i|nmx|0>
        for(pij=p00, ni0=n00; pij<pend; ni0+=c)		// Effective loop over i
          for(dij=d00,nij=ni0; dij<dend; nij++,pij++,dij++)	// Effective loop over j
            *pij = *nij * *dij;				// <i|pdt|j> = <i|nmx|j><j|dmx|j>
//      std::cerr << "NON BLAS: n_matrix * d_matrix\n";
        return pdt;
        break;
      }
      
      case i_matrix_type:				// Multiply i_matrix into n_matrix
      {
        return this;					// Return nmx  
        break;
      }
      
      default:						// Multiply generic matrix into n_matrix
      { 						// Cannot use mx internal structure here
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
          cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A_rows, B_cols, A_cols, 
                             alpha, (double *)data, A_cols, (double *)((n_matrix*)mx)->data, B_cols, beta, (double *)pdt->data, C_cols);
//        std::cerr << "BLAS: n_matrix * unknown\n";
        }
        else
        {
         // Duplicate of "original" code from n_matrix.multiply() as of 05/11/10.
          n_matrix* pdt =					// First construct ne normal matrix
          new n_matrix(rows(),mx->cols(),complex0);
  
          int pos = 0;					// Position of 1st element in pdt
          int pos1 = 0;					// Position of 1st element in nmx
          for(int i=0;  i<pdt->rows(); i++, pos1+=cols())	// Loop over the rows of nmx (&pdt) 
            for(int j=0; j<pdt->cols(); j++,pos++)	// Loop over the columns mx (&pdt)
              for(int k=0; k<cols(); k++)			// Loop over the inner index
                    pdt->data[pos] += this->data[pos1+k]*(*mx)(k,j);                  
//      std::cerr << "NON BLAS: n_matrix * unknown\n";
        }
#else
        n_matrix* pdt =					// First construct ne normal matrix
        new n_matrix(rows(),mx->cols(),complex0);
        int pos = 0;					// Position of 1st element in pdt
        int pos1 = 0;					// Position of 1st element in nmx
        for(int i=0;  i<pdt->rows(); i++, pos1+=cols())	// Loop over the rows of nmx (&pdt) 
          for(int j=0; j<pdt->cols(); j++,pos++)	// Loop over the columns mx (&pdt)
            for(int k=0; k<cols(); k++)			// Loop over the inner index
                  pdt->data[pos] += this->data[pos1+k]*(*mx)(k,j);                  
#endif
        return pdt;
	  }
    }
  }

  // This shouldn't be needed as it should never get here, 
  // but it get's rid of a warning message.
  return mx;
}


_matrix* n_matrix::multiply(const complex& z)

        // Input                nmx  : Input normal matrix (this)
	//			z    : Complex number
        // Output               pdt  : New matrix which is the product z*nmx

  {
  if(z == complex1)					// Check if z is 1
    return this;					// If so, return nmx
  else if(z == complex0)				// Check if z is 0
    {							// If so, return zero matrix
    n_matrix* mx = new n_matrix(rows_,cols_,complex0);
    return mx;
    }
  n_matrix* mx = new n_matrix(rows_,cols_);		// If not, construct new matix
  for(int i=0; i<size; i++)				// Loop over all the elements
    mx->data[i] = z*data[i];				// <i|mx|j> = z * <i|nmx|j>
  return mx;
  }


_matrix* n_matrix::divide(_matrix* mx)

        // Input                nmx  : Input normal matrix (this)
	//			 mx  : Second matrix                         -1
        // Output               pdt  : New matrix which is the product nmx*mx
        // Note                      : This uses the function inv of mx

//                            ---                -1
//                <i|pdt|j> = \   <i|nmx|k> <k|mx  |j>
//                            /
//                            ---
//		               k
  {
  if(cols() != mx->rows())				// Insure nmx & mx dimension match
    {
    NMxerror(51,1);					//    Mismatched dimensions
    NMxfatal(23);					//    Cannot perform division
    return mx;
    }
  else
    switch(mx->stored_type())
      {
      case i_matrix_type:				// Divide i_matrix into n_matrix
	return this;					// Return nmx  
	break;
      case d_matrix_type:				// Divide n_matrix by d_matrix
	{
	int r = cols();					// Rows of product matrix
	int c = rows();					// Cols of product matrix
	n_matrix* pdt = new n_matrix(r,c); 		// Create new normal matrix
	complex *d00 = ((d_matrix*)mx)->data;		// Start of data of dmx
	complex *dend = d00 + c;			// End of data in dmx
        complex* dinv = new complex[cols()];		// Array for dmx inverse
        complex *djj, *dinvjj = dinv;
        for(djj=d00; djj<dend; djj++, dinvjj++)		// Loop over all columns &
          {						// set <j|dinv|j> = 1/<j|dmx|j>
          if(*djj == complex0)
            {
            NMxerror(27,1);				//    Inverse doesn't exist
            NMxfatal(23);				//    Cannot perform division
            }
          *dinvjj = 1.0/(*djj); 
          }
	complex *p00 = pdt->data; 			// Start of data of pdt
	complex *n00 = data;				// Start of data of nmx
	complex *pend = p00 + r*c;			// End of data of pdt
	complex *pij, *nij, *ni0;			// Matrix elements & <i|nmx|0>
	dend = dinv + c;				// End of data in dinv
	for(pij=p00, ni0=n00; pij<pend; ni0+=c)		// Effective loop over i
	  for(djj=dinv,nij=ni0; djj<dend; nij++,pij++,djj++)	// Effective loop over j
	    *pij = *nij * *djj;				// <i|pdt|j> = <i|nmx|j><j|inv(dmx)|j>
        delete dinv;
	return pdt;
	}
      case h_matrix_type:                               // Divide n_matrix by h_matrix
        return multiply((h_matrix*)mx->inv());          // Return nmx*inv(hmx)
        break;
      case n_matrix_type:                               // Divide n_matrix by n_matrix
        return multiply((n_matrix*)mx->inv());          // Return nmx*inv(nmx)
        break;
      default:						// Multiply n_matrix into n_matrix
        break;
      }
  NMxerror(25,"div",1);					// Division not fully implemented
  NMxfatal(23);						// Unable to do division
  return mx;
  }


_matrix* n_matrix::divide(const complex& z)

        // Input                nmx  : Input normal matrix (this)
	//			z    : Complex number
        // Output               pdt  : New matrix: (1/z)*nmx
	// Note			     : Handles code mx/z

  {
  if(z == complex1) return this; 		// If z is 1 return nmx
  else if(z == complex0)			// If z is 0, we have trouble
    {
    NMxerror(13,1);				//    Divisions by zero
    NMxfatal(23);				//    Cannot perform division
    }
  n_matrix* mx = new n_matrix(rows_,cols_);	// If not, construct new matix
  complex z1 = 1/z;				// Calculate 1/z
  for(int i=0; i<size; i++)			// Loop over all the elements
    mx->data[i] = z1*data[i];			// <i|mx|j> = (1/z) * <i|nmx|j>
  return mx;
  }


// ____________________________________________________________________________
// G              CLASS N_MATRIX UNARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________


_matrix* n_matrix::add_two(_matrix* mx)

        // Input                nmx  : Input normal matrix (this)
	//			mx   : Second matrix
        // Output               mx   : Matrix mx is added into nmx
	// Note			     : This is the implementation of mx += nmx

  {
  if((mx->stored_type()) == n_matrix_type)	 // Add two n_matrices
    { 
    if((rows_!=mx->rows()) || (cols_!=mx->cols()))
      {
      NMxerror(51,1);				//    Mismatched dimensions
      NMxfatal(20);				//    Cannot perform addition
      }
    complex *a=&((n_matrix*)mx)->data[size-1];	// Pointer to last element of mx
    complex *b = &data[size-1];			// Pointer to last element of this
    for(; b>=data ;b--,a--)			// Perform addition elementwise
      (*a) += (*b);
    ((n_matrix*)mx)->unitary = false;		// This is probably not unitary now
    return mx;
    }
  else						// Add n_matrix to non-n_matrix, since
    return add(mx); 				// matrix type changes, result has to be
  } 						// copied anyway -> use add function


_matrix* n_matrix::subtract_two(_matrix* mx)

  {
  if ((mx->stored_type()) == n_matrix_type)	 // Subtract n_matrix from n_matrix
    { 
    if((rows()!=mx->rows()) || (cols()!=mx->cols()))
      {
      NMxerror(51,1);				//    Mismatched dimensions
      NMxfatal(21);				//    Cannot perform subtraction
      }
    complex *a=&((n_matrix*)mx)->data[size-1];	// Pointer to last element of mx
    complex *b=&data[size-1];			// Pointer to last element of this
    for(; b>=data; b--,a--)			// Subtraction elementwise
      (*a)-=(*b);
    ((n_matrix*)mx)->unitary = false;		// This is probably not unitary now
    return mx;
    }
  else						// Subtract n_matrix to non-n_matrix, since
    return mx->subtract(this);			// matrix type changes, result has to be
  } 						// copied anyway -> use subtract function


_matrix* n_matrix::multiply_two(_matrix* mx) 	// In this case, mx cannot be used for the
						// for the output.  If mx=n_matrix then its
  { return mx->multiply(this); }		// not allowed and if not n_matrix, type must
						// by modified (usually to n_matrix)!



        // Input                nmx  : Input normal matrix (this)
	//			z    : Complex number
        // Output               nmx  : This matrix scaled by the complex number
	// Note			     : This is the implementation of mx *= z

_matrix* n_matrix::multiply_two(const complex &z)
  {
  for(int i=0; i<size; i++)
    data[i] = z*data[i];
  if(z != complex1) unitary=false;
  return this;
  }


_matrix* n_matrix::divide_two(_matrix* mx) { return mx->divide(this); }

_matrix* n_matrix::divide_two(const complex &z)

        // Input                nmx  : Input normal matrix (this)
	//			z    : Complex number
        // Output               nmx  : This matrix divided by the complex number
	// Note			     : This is the implementation of mx /= z

  {
  complex z1 = 1/z;
  for(int i=0; i<size; i++)
    data[i] = z1*data[i];
  if(z != complex1) unitary=false;
  return this;
  }

// ___________________________________________________________________________
// H              CLASS N_MATRIX SIMPLE UNARY FUNCTIONS
// ___________________________________________________________________________

/* These functions perform simple simple mathematical operations on a Hermitian
   matrix. Note that for Hermitian arrays the adjoint function does nothing.

    Function          Result             Function           Result
   ---------- -----------------------    --------  ------------------------
   negate     <i|nmx|j> -> -<i|nmx|j>       RE      <i|nmx|j>->Re(<i|nmx|j>)
   conjugate  <i|nmx|j> -> <i|nmx*|j>       IM      <i|nmx|j>->Im(<i|nmx|j>)
   transpose  <i|nmx|j> -> <j|nmx|i>      adjoint   <i|nmx|j>-><j|nmx*|i>   */

_matrix* n_matrix::negate()
  {
  n_matrix* mx = new n_matrix(rows_, cols_);	// Construct a new n_matrix
  for(int pos=0; pos<size; pos++)		// Fill elements with negated values
    mx->data[pos] = -data[pos];			// <i|mx|j> = -<i|nmx|j>
  mx->unitary = unitary;
  return mx;
  }
  
_matrix* n_matrix::RE()
  {
  if(is_real()) return this; 			// Return this if already real
  n_matrix* mx = new n_matrix(rows_, cols_);	// Construct a new n_matrix
  for(int pos=0; pos<size; pos++)		// Fill elements with real values
      mx->data[pos] = Re(data[pos]);
  mx->unitary = false;				// No longer unitary
  return mx;
  }
  
_matrix* n_matrix::IM()
  {
  if(is_imaginary()) return this; 		// Return this if already imaginary
  n_matrix* mx = new n_matrix(rows_, cols_);	// Construct a new n_matrix
  for(int pos=0; pos<size; pos++)		// Fill elements with imaginary part
    mx->data[pos] = Im(data[pos]);
  mx->unitary = false;				// No longer unitary
  return mx;
  }
  
_matrix* n_matrix::conjugate()
  {
  if(is_real()) return this; 			// Return this if just real
  n_matrix* mx = new n_matrix(rows_, cols_);	// Construct a new n_matrix
  for(int pos=0; pos<size; pos++)		// Fill elements with conjugated values
    mx->data[pos] = conj(data[pos]);		// <i|mx|j> = <i|nmx*|j>
  mx->unitary = unitary;
  return mx;
  }

_matrix* n_matrix::transpose()
  {
  n_matrix* tmx = new n_matrix(cols_, rows_);	// Construct a new n_matrix
  complex *n00 = data;				// Start of nmx: n00 = <0|nmx|0>
  complex *t00 = tmx->data;	 		// Start of tmx: t00 = <0|tmx|0>
  complex *t10 = t00 + rows_;	 		// 2nd tmx col : t10 = <1|tmx|0>
  complex *tend = t00 + cols_*rows_;		// End of data in tmx matrix
  complex *nij, *t0i, *tji;
  for(t0i=t00,nij=n00; t0i<t10; t0i++)		// Effective loop over i
    for(tji=t0i; tji<tend; tji+=rows_,nij++)	// Effective loop over j
      (*tji) = (*nij);				// <j|tmx|i> = <j|nmx|i>
  tmx->unitary = unitary;			// If nmx unitary, tmx unitary
  return tmx;
  }

_matrix* n_matrix::adjoint()
  {
  n_matrix* amx = new n_matrix(cols_, rows_);	// Construct a new n_matrix
  complex *n00 = data;				// Start of nmx: n00 = <0|nmx|0>
  complex *a00 = amx->data;	 		// Start of amx: a00 = <0|amx|0>
  complex *a10 = a00 + rows_;	 		// 2nd amx col : a10 = <1|amx|0>
  complex *aend = a00 + cols_*rows_;		// End of data in amx matrix
  complex *nij, *a0i, *aji;
  for(a0i=a00,nij=n00; a0i<a10; a0i++)		// Effective loop over i
    for(aji=a0i; aji<aend; aji+=rows_,nij++)	// Effective loop over j
      (*aji) = conj(*nij);			// <j|amx|i> = <j|nmx*|i>
  amx->unitary = unitary;			// If nmx unitary, amx is too
  return amx;
  }

/*  This curretly uses an eigensystem method for generating the exponential.

                   exp(A) = exp[E*D*inv(E)] = E*exp(D)*inv(E)

    Of course, this will fail when the eigensystem is ill-defined and there
    are other faster ways to approximate the exponential matrix. Hopefully
    later we can find the time to implement some of them too.                */

_matrix* n_matrix::mxexp()
  {
  _matrix* dmx;					// For eigenvalues
  _matrix* emx;					// For eigenvectors
  diag(dmx, emx);				// Diagonalize mx to dmx & emx 
  complex *d00 = ((d_matrix*)dmx)->data;	// Start of dmx: d00 = <0|dmx|0>
  complex *dii, *dend = d00 + cols_;		// Indexing on dmx diagonal
  for(dii=d00; dii<dend; dii++)			// Effective loop over i
    (*dii) = exp(*dii);				// Set <i|dmx|i> to exp(<i|dmx|i>
//  _matrix* pdt1  = dmx->multiply(emx);		// Perform exp(D)*E
//  _matrix* emxi  = emx->inv();			// Get inverse of emx
//  _matrix* expmx = emxi->multiply(pdt1);	// Perform inv(E)*exp(D)*E

  _matrix* emxi  = emx->inv();			// Get inverse of emx
  _matrix* pdt1  = dmx->multiply(emxi);
  _matrix* expmx = emx->multiply(pdt1);
  delete dmx;                                   // Clean up the dmx matrix
  delete emx;                                   // Clean up the emx matrix
  delete pdt1;					// Clean up the pdt1 matrix
  delete emxi;					// Clean up the emxi matrix
  return expmx;
  }


complex n_matrix::trace()
  {
  if(rows_ != cols_)
     {
     NMxerror(14,1);			//    Bad rectangular array use
     NMxfatal(24);			//    Cannot take trace
     }
  complex z(0);				// Initialize to zero
  for(int i=0; i<rows_; i++)		// Sum over the diagonal elements
    z += data[i*cols_+i];
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
 
_matrix* n_matrix::swaprows(int I, int J)
  {
  n_matrix* nsr = new n_matrix(*this);          // Make new n_matrix
  int I0 = I*cols_;				// Index of <I|mx,mx'|0>
  int J0 = J*cols_;				// Index of <J|mx,mx'|0>
  for(int k=0; k<cols_; k++)			// Swap rows I & J
    {
    nsr->data[I0+k] = data[J0+k];       	// Set <I|mx'|k> = <J|mx|k>
    nsr->data[J0+k] = data[I0+k];               // Set <J|mx'|k> = <I|mx|k>
    }
  return nsr;
  }

_matrix* n_matrix::swapcols(int I, int J)
  {
  n_matrix* nsc = new n_matrix(*this);  	// Make new n_matrix
  for(int k=0; k<rows_; k++)			// Swap columns I & J
    {
    nsc->data[k*cols_+I] = data[k*cols_+J];	// Set  <k|mx|I> = <k|mx|J>
    nsc->data[k*cols_+J] = data[k*cols_+I];	// Set  <k|mx|J> = <k|mx|I> orig
    }
  return nsc;
  }

_matrix* n_matrix::permute(int I, int J)
  {
  complex z;					// Temp for element swap
  n_matrix* pmx = (n_matrix*)swaprows(I,J);	// First swap rows I & J
  for(int k=0; k<rows_; k++)			// Now swap columns I & J
    {
    z = pmx->data[k*cols_+I];			// Save <k|nsr|I>
    pmx->data[k*cols_+I] = pmx->data[k*cols_+J];// Set  <k|mx|I> = <k|nsr|J>
    pmx->data[k*cols_+J] = z;			// Set  <k|mx|J> = <k|nsr|I>
    }
  return pmx;
  }

double n_matrix::maxRe() const
  {
  double maxval=-HUGE_VAL;
  for(int i=0; i<size; i++) maxval = gmax(Re(data[i]), maxval);
  return (size)?maxval:0;
  }

double n_matrix::maxIm() const
  {
  double maxval=-HUGE_VAL;
  for(int i=0; i<size; i++) maxval = gmax(Im(data[i]), maxval);
  return (size)?maxval:0;
  }
 
complex n_matrix::maxZ()  const
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

double n_matrix::minRe() const
  {
  double minval=HUGE_VAL;
  for(int i=0; i<size; i++) minval = gmin(Re(data[i]), minval);
  return (size)?minval:0;
  }
 
double n_matrix::minIm() const
  {
  double minval=HUGE_VAL;
  for(int i=0; i<size; i++) minval = gmin(Im(data[i]), minval);
  return (size)?minval:0;
  }
 
complex n_matrix::minZ() const
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

// ______________________________________________________________________
// I              CLASS N_MATRIX SIMPLE BINARY FUNCTIONS
// ______________________________________________________________________


complex n_matrix::trace(_matrix* mx)

        // Input            nmx : An n_matrix (this)
        //                   mx : A second matrix
        // Output            z  : A complex value which is the trace of
        //                        the product of two input arrays
        //                         z = Tr{nmx*mx} = Tr{mx*nmx}
        // Note                 : This is faster than taking the matrix
        //                        product and then using unary trace!
        // Note                 : Uses internal structure of other matrix types

//       ---                 --- ---
//   z = \   <i|nmx*mx|i>  = \   \   <i|nmx|j><j|mx|i>
//       /                   /   /
//       ---                 --- ---
//        i                   i   j

  {
  if(!CheckDim(mx,1,1) || !CheckDim(mx,3,1))	// Insure result is square
    NMxfatal(24);				// If not, cannot do trace
  complex z(0);				// Initialize the result at zero
  switch(mx->stored_type())
    {
    case n_matrix_type: 			// Tr{n_matrix * n_matrix}    
	{
        complex *n00 = data;			// Start of nmx: n00 = <0|nmx|0>
        complex *m00 = ((n_matrix*)mx)->data;	// Start of mx: m00 = <0|mx|0>
        complex *nz = n00 + rows_*cols_; 	// End of nmx: nz = <rows_|nmx|cols_+1>
        complex *mz = m00 + rows_*cols_; 	// End of mx: mz = <cols_|mx|rows_+1>
        complex *nij, *mji, *m0i;		// These are the matrix elements
	for(nij=n00,m0i=m00; nij<nz; m0i++)	// Effective loop over i
	  for(mji=m0i; mji<mz; nij++,mji+=rows_)// Effective loop over j
	    z += (*nij) * (*mji); 		// Tr{nmx*mx} += <i|nmx|j><j|mx|i>
	return z;
	}
	break;
    case h_matrix_type: 			// Tr{n_matrix * h_matrix}    
	{
        complex *n00 = data;			// Start of nmx: n00 = <0|nmx|0>
        complex *h00 = ((h_matrix*)mx)->data;	// Start of hmx: h00 = <0|hmx|0>
        complex *nij, *hij, *hji, *h0i;		// These are the matrix elements
        int i=0, j;
        hij = h00;                              // Keeps the compiler happy
        z = complex0;                           // Keeps the compiler happy
	for(h0i=h00,nij=n00; i<rows_; i++,h0i++)// Loop rows of nmx, cols of hmx
          {					// Split loop of nmx cows, hmx rows
	  for(hji=h0i,j=0; j<=i;		// to be over upper triangle of hmx
                j++,nij++,hij=hji,hji+=cols_-j)	// (i.e. hmx when j<=i, hmx* j>i)
	    z += (*nij) * (*hji);		// Tr{nmx*hmx} += <i|nmx|j><j|hmx|i>, j<i
          for(hij++; j<cols_; j++,nij++,hij++)
	    z += (*nij) * conj(*hij); 		// Tr{nmx*hmx} += <i|nmx|j><i|hmx*|j>, j>=i
          }
	return z;
	}
	break;
    case d_matrix_type: 			// Tr{n_matrix * d_matrix}    
	{
        complex *n00 = data;			// Start of nmx: n00 = <0|nmx|0>
        complex *d00 = ((d_matrix*)mx)->data;	// Start of dmx: d00 = <0|dmx|0>
        complex *dend = d00 + rows_;		// End of dmx: dend = <rows_|dmx|rows+1>
        complex *nii, *dii;			// These are the matrix elements
        int incd = cols_+1;			// <i|nii|i>+incd -> <i+1|nii|i+1>
	for(nii=n00,dii=d00; dii<dend;
                             nii+=incd, dii++)
          z += (*nii) * (*dii);			// Tr{nmx*dmx} += <i|nmx|i><i|dmx|i>
	return z;
	}
	break;
    case i_matrix_type: 			// Tr{n_matrix * i_matrix}    
	{
        for(int i=0; i<rows_; i++)		// Sum over the diagonal elements
          z += data[i*cols_+i];			// Tr{nmx*imx} += <i|nmx|i><i|imx|i>
        return z; 				//             += <i|nmx|i>
	}
	break;
    default:                 			// Tr{nmx * unknown matrix}
	{
        complex *n00 = data;			// Start of nmx: n00 = <0|nmx|0>
        complex *nij;				// This is matrix element <i|nmx|j>
        int i,j;
	for(i=0, nij=n00; i<rows_; i++)
	  for(j=0; j<cols_; j++, nij++)
	    z += (*nij) * (*mx)(j,i); 		// Tr{nmx*mx} += <i|nmx|j><j|mx|i>
	return z;
	}
    }
  }


_matrix* n_matrix::adjoint_times(_matrix* mx)
 
        // Input            nmx : An n_matrix (this)
        //                   mx : A second matrix
        // Output           mx1 : A matrix which is product of the
        //                        adjoint of n_matrix and mx
        //                                    T *
        //                         mx1 = [(nmx ) ] * mx
        // Note                 : This is faster than taking the adjoint
        //                        of nmx and then the product!  Use if the
        //                        adjoint of nmx is not needed. 
	// Note 		: The algorithm used is identical to
	// 			  multiply, but the used indices of nmx
	//			  are switched and conj(*a,*b) is used
	//			  where conj(z1,z2) = conj(z1)*z2!
	// Note		        : This uses knowledge of the internal matrix
	//			  structure of other classes

/*                                ---
                     <i|pdt|j>  = \   <k|nmx*|i><k|mx|j>
                                  /
                                  ---
                                   k                                         */

  {
  if(rows_ != mx->rows())				// Insure matrix dimensions match 
    {
    NMxerror(51,1);					//    Mismatched dimensions
    NMxfatal(5,"adjoint_times");			//    Bad use of function
    return mx;
    }
  else
    switch(mx->stored_type())
      {
      case n_matrix_type:				// Multiply n_matrix adjoint into n_matrix
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
		cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, A_rows, B_cols, A_cols, 
                             alpha, (double *)data, A_cols, (double *)((n_matrix*)mx)->data, B_cols, beta, (double *)pdt->data, C_cols);
//              std::cerr << "BLAS: n_matrix * n_matrix\n";
            }
            else
            {
	      int r = cols_;					// Rows of product matrix = cols of nmx 
	      int c = mx->cols();				// Columns of created matrix = cols of mx
	      int l = rows_;					// Dimension of inner loop
	      complex *p00 = pdt->data;			// Start of data of pdt
	      complex *n00 = data;				// Start of data of nmx
	      complex *m00 = ((n_matrix*)mx)->data;		// Start of data of mx
	      complex *pend = p00 + r*(mx->cols());		// End of data in pdt
	      complex *mend = m00 + c*l;			// End of data in mx
	      complex *m0r = m00 + c;				// End of first row in mx
	      complex *pij, *nki=n00, *mkj, *m0k, *n0i;
	      complex z;					// Intermediate storage of <i|pdt|j>
	      for(pij=p00, n0i=n00; pij<pend; n0i++)		// Effective loop over i
	        for(m0k=m00; m0k<m0r; m0k++, pij++)		// Effective loop over j
	          {
	          for(z=0,mkj=m0k,nki=n0i;			// Effective loop over k
                            mkj<mend; nki+=r, mkj+=c)
	            z += conj(*nki,*mkj);			// z += <k|nmx*|i><k|mx|j>
	          *pij = z;
	          }
	    }
#else
	int r = cols_;					// Rows of product matrix = cols of nmx 
	int c = mx->cols();				// Columns of created matrix = cols of mx
	int l = rows_;					// Dimension of inner loop
	n_matrix* pdt = new n_matrix(r,c); 		// Create new n_matrix for product
	complex *p00 = pdt->data;			// Start of data of pdt
	complex *n00 = data;				// Start of data of nmx
	complex *m00 = ((n_matrix*)mx)->data;		// Start of data of mx
	complex *pend = p00 + r*(mx->cols());		// End of data in pdt
	complex *mend = m00 + c*l;			// End of data in mx
	complex *m0r = m00 + c;				// End of first row in mx
	complex *pij, *nki=n00, *mkj, *m0k, *n0i;
	complex z;					// Intermediate storage of <i|pdt|j>
	for(pij=p00, n0i=n00; pij<pend; n0i++)		// Effective loop over i
	  for(m0k=m00; m0k<m0r; m0k++, pij++)		// Effective loop over j
	    {
	    for(z=0,mkj=m0k,nki=n0i;			// Effective loop over k
                      mkj<mend; nki+=r, mkj+=c)
	      z += conj(*nki,*mkj);			// z += <k|nmx*|i><k|mx|j>
	    *pij = z;
	    }
#endif
	  return pdt;
	}
	break;
      case d_matrix_type:				// Multiply n_matrix adjoint into d_matrix
	{
	n_matrix* pdt = new n_matrix(cols_,rows_);	// Create new n_matrix
	complex *p00 = pdt->data;			// Start of data in pdt: <0|pdt|0>
	complex *n00 = data;				// Start of data in nmx: <0|nmx|0>
	complex *d00 = ((d_matrix*)mx)->data;		// Start of data in dmx: <0|dmx|0>
	complex *pend = p00 + cols_*rows_;		// End of data in pdt: <cols_|pdt|rows_+1>
	complex *dend = d00 + rows_;			// End of data in dmx: <rows_|dmx|rows_+1>
	complex *pij, *nji, *djj, *n0i;			// These are the matrix elements
	for(pij=p00, n0i=n00; pij<pend; n0i++)		// Effective loop over i
	  for(djj=d00,nji=n0i;				// Effective loop over j
                     djj<dend; nji+=cols_,pij++,djj++)
	    *pij = conj(*nji, *djj);			// <i|pdt|j> = <j|nmx*|i> <j|dmx|j>
	return pdt;
	}
	break;
      case i_matrix_type:				// Multiply n_matrix adjoint into i_matrix
	return adjoint();				// This is just the adjoint
	break;
      case h_matrix_type:				// Multiply n_matrix adjoint into h_matrix
	{ 
#ifdef _USING_BLAS_
         int A_rows = rows();
         int A_cols = cols();            
         int B_rows = mx->rows();
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
            n_matrix* hmx =	new n_matrix(B_rows,B_cols);		// Create new matrix h_matrix
            mx->convert(hmx);				// Convert h_matrix mx into normal matrix hmx
	    cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, A_rows, B_cols, A_cols, 
                             alpha, (double *)data, A_cols, (double *)hmx->data, B_cols, beta, (double *)pdt->data, C_cols);
//          std::cerr << "BLAS: n_matrix * h_matrix\n";
	    delete hmx;
         }
         else
         {
	   complex *p00 = pdt->data;			// Start of data in pdt: <0|pdt|0>
	   complex *n00 = data;				// Start of data in nmx: <0|nmx|0>
	   complex *n10 = n00 + cols_;			// End of 1st nmx row: <1|nmx|0>
	   complex *h00 = ((h_matrix*)mx)->data;		// Start of data in hmx: <0|hmx|0>
	   complex *h11 = h00 + rows_;			// End of 1st hmx row: <1|hmx|1>
	   complex *pij,*nki,*n0i,*hjk,*hkj,*h0j,*hjj;	// These are the matrix elements
           int hrow;
	   for(pij=p00,n0i=n00; n0i<n10; n0i++)		// Effective loop over i (via n0i)
	     for(h0j=h00,hjj=h00; h0j<h11; pij++,h0j++)	// Effective loop over j (via h0j)
               {
               *pij = complex0;				// Initialize pij: <i|pdt|j> = 0
               nki = n0i;					// Set <k|nmx|i> = <0|nmx|i>
               hkj = h0j;					// Set <k|hmx|j> = <0|hmx|j>
               hrow = rows_;				// Row length of hmx
	       for( ; hkj<hjj; nki+=cols_,hrow--,hkj+=hrow)// Loop over the internal index, k<j (via hkj) 
	         (*pij) += conj(*nki, *hkj);		// <i|pdt|j> += <k|nmx*|i><k|hmx|j>
               hjj += hrow; 				// <j|hmx|j> -> <j+1|hmx|j+1>
	       for(hjk=hkj; hjk<hjj; nki+=cols_,hjk++)	// Continue loop over the internal index, k>=j 
	         (*pij) += conj(*nki) * conj(*hjk);	// <i|pdt|j> += <k|nmx*|i><j|hmx*|k>
               }
         }
#else
	n_matrix* pdt = new n_matrix(cols_,rows_);	// Create a new n_matrix, initialized to zero
	complex *p00 = pdt->data;			// Start of data in pdt: <0|pdt|0>
	complex *n00 = data;				// Start of data in nmx: <0|nmx|0>
	complex *n10 = n00 + cols_;			// End of 1st nmx row: <1|nmx|0>
	complex *h00 = ((h_matrix*)mx)->data;		// Start of data in hmx: <0|hmx|0>
	complex *h11 = h00 + rows_;			// End of 1st hmx row: <1|hmx|1>
	complex *pij,*nki,*n0i,*hjk,*hkj,*h0j,*hjj;	// These are the matrix elements
        int hrow;
	for(pij=p00,n0i=n00; n0i<n10; n0i++)		// Effective loop over i (via n0i)
	  for(h0j=h00,hjj=h00; h0j<h11; pij++,h0j++)	// Effective loop over j (via h0j)
            {
            *pij = complex0;				// Initialize pij: <i|pdt|j> = 0
            nki = n0i;					// Set <k|nmx|i> = <0|nmx|i>
            hkj = h0j;					// Set <k|hmx|j> = <0|hmx|j>
            hrow = rows_;				// Row length of hmx
	    for( ; hkj<hjj; nki+=cols_,hrow--,hkj+=hrow)// Loop over the internal index, k<j (via hkj) 
	      (*pij) += conj(*nki, *hkj);		// <i|pdt|j> += <k|nmx*|i><k|hmx|j>
            hjj += hrow; 				// <j|hmx|j> -> <j+1|hmx|j+1>
	    for(hjk=hkj; hjk<hjj; nki+=cols_,hjk++)	// Continue loop over the internal index, k>=j 
	      (*pij) += conj(*nki) * conj(*hjk);	// <i|pdt|j> += <k|nmx*|i><j|hmx*|k>
            }
#endif
	return pdt;
        }
        break;
      default:						// Multiply n_matrix adjoint into generic matrix
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
          cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans, A_rows, B_cols, A_cols, 
                             alpha, (double *)data, A_cols, (double *)((n_matrix*)mx)->data, B_cols, beta, (double *)pdt->data, C_cols);
//        std::cerr << "BLAS: n_matrix * unknown\n";
        }
        else
        {
  	  int c = mx->cols();				// Columns of created matrix = cols of mx
	  complex *p00 = pdt->data;			// Start of data in pdt: <0|pdt|0>
	  complex *n00 = data;				// Start of data in nmx: <0|nmx|0>
	  complex *pij, *nki, *n0i;			// These are the matrix elements
          int i,j,k;
	  for(i=0,pij=p00,n0i=n00; i<cols_; i++,n0i++)	// Loop over all the rows of pdt
	    for(j=0; j<c; j++,pij++)			// Loop over all the columns of pdt
	      for(k=0,nki=n0i; k<rows_; k++,nki+=cols_)	// Loop over the internal index 
	        (*pij) += conj(*nki, (*mx)(k,j));		// <i|pdt|j> += <k|nmx*|i><k|mx|j>

        }
#else
	int c = mx->cols();				// Columns of created matrix = cols of mx
	n_matrix* pdt = new n_matrix(cols_,c,complex0);	// Create a new n_matrix, initialized to zero
	complex *p00 = pdt->data;			// Start of data in pdt: <0|pdt|0>
	complex *n00 = data;				// Start of data in nmx: <0|nmx|0>
	complex *pij, *nki, *n0i;			// These are the matrix elements
        int i,j,k;
	for(i=0,pij=p00,n0i=n00; i<cols_; i++,n0i++)	// Loop over all the rows of pdt
	  for(j=0; j<c; j++,pij++)			// Loop over all the columns of pdt
	    for(k=0,nki=n0i; k<rows_; k++,nki+=cols_)	// Loop over the internal index 
	      (*pij) += conj(*nki, (*mx)(k,j));		// <i|pdt|j> += <k|nmx*|i><k|mx|j>
#endif
	return pdt;
	}
      }
  }


_matrix* n_matrix::times_adjoint(_matrix* mx)

        // Input            nmx : An n_matrix (this)
        //                   mx : A second matrix
        // Output           mx1 : A matrix which is product of the
        //                        n_matrix and the adjoint of mx
        //                                        T *
        //                         mx1 = nmx* [(mx ) ]
        // Note                 : This is faster than taking the adjoint
        //                        of mx and then the product!  Use if the
        //                        adjoint of mx is not needed.
	// Note 		: The algorithm used here is identical to
	// 			  multiply except the indices of mx are
	//			  switched and conj(*a,*b) is used!
	// Note		        : This uses knowledge of the internal matrix
	//			  structure of other classes

/*                               ---
                    <i|pdt|j>  = \   <i|nmx|k><j|mx*|k>
                                 /   
                                 ---
                                  k                                          */   
  {
  if(cols_ != mx->cols())			// Insure dimensions match
    {
    NMxerror(51,1);				//    Mismatched dimensions
    NMxfatal(5,"times_adjoint");		//    Bad use of function
    return mx;
    }
  else
    switch(mx->stored_type())
      {
      case n_matrix_type:		  	// Multiply n_matrix into n_matrix adjoint
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
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, A_rows, B_cols, A_cols, 
                             alpha, (double *)data, A_cols, (double *)((n_matrix*)mx)->data, B_cols, beta, (double *)pdt->data, C_cols);
//              std::cerr << "BLAS: n_matrix * n_matrix\n";
            }
            else
            {
	      int c = mx->rows();			// Columns of product matrix
	      int l = cols_;				// Dimension of inner loop
	      complex *p00 = pdt->data;		// Start of data in pdt: <0|pdt|0>
	      complex *n00 = data;			// Start of data in nmx: <0|nmx|0>
	      complex *m00 = ((n_matrix*)mx)->data;	// Start of data in mx: <0|mx|0>
	      complex *pend = p00 + rows_*c;		// End of data in pdt: <rows_|pdt|c+1>
	      complex *mend = m00 + c*cols_;		// End of data in mx: <c|mx|cols_+1>
	      complex *pij, *nik=n00, *ni0, *mjk, *mj0;
	      complex z;				// Intermediate storage
	      for(pij=p00, ni0=n00; pij<pend; ni0+=l)	// Effective loop over i
	        for(mj0=m00; mj0<mend; mj0+=l, pij++)	// Effective loop over j
                {
	          for(z=0, nik=ni0+l-1, mjk=mj0+l-1;	// Effective loop over k
                               nik>=ni0; nik--,mjk--)
	            z += conj(*mjk, *nik);		// <i|pij|j> += <i|nmx|k><j|mx*|k>
	          *pij = z;
	        }
	    }
#else
	int c = mx->rows();			// Columns of product matrix
	int l = cols_;				// Dimension of inner loop
	n_matrix* pdt = new n_matrix(rows_,c);	// Create a new n_matrix for result
	complex *p00 = pdt->data;		// Start of data in pdt: <0|pdt|0>
	complex *n00 = data;			// Start of data in nmx: <0|nmx|0>
	complex *m00 = ((n_matrix*)mx)->data;	// Start of data in mx: <0|mx|0>
	complex *pend = p00 + rows_*c;		// End of data in pdt: <rows_|pdt|c+1>
	complex *mend = m00 + c*cols_;		// End of data in mx: <c|mx|cols_+1>
	complex *pij, *nik=n00, *ni0, *mjk, *mj0;
	complex z;				// Intermediate storage
	for(pij=p00, ni0=n00; pij<pend; ni0+=l)	// Effective loop over i
	  for(mj0=m00; mj0<mend; mj0+=l, pij++)	// Effective loop over j
	    {
	    for(z=0, nik=ni0+l-1, mjk=mj0+l-1;	// Effective loop over k
                         nik>=ni0; nik--,mjk--)
	      z += conj(*mjk, *nik);		// <i|pij|j> += <i|nmx|k><j|mx*|k>
	    *pij = z;
	    }
#endif
	return pdt;
	}
	break;
      case d_matrix_type:				// Multiply n_matrix into d_matrix adjoint
	{
	int r = rows_;					// Rows of product matrix
	int c = cols_;					// Columns of product matrix
	n_matrix* pdt = new n_matrix(r,c);		// Create new normal matrix
	complex *p00 = pdt->data;			// Start of data in pdt
	complex *n00 = data;				// Start of data in nmx
	complex *m00 = ((d_matrix*)mx)->data;		// Start of data in mx
	complex *pend = p00 + r*c;			// End of data in pdt
	complex *mend = m00 + c;			// End of data in mx (uses d_matrix structure)
	complex *pij, *nij, *ni0, *mjj;
	for(pij=p00, ni0=n00; pij<pend; ni0+=c)		// Effective loop over i
	  for(mjj=m00,nij=ni0; mjj<mend;		// Effective loop over j
                           nij++, pij++, mjj++)
	    *pij = conj(*mjj,*nij);			// <i|pdt|j> = <i|nmx|k><j|mx*|k>
	return pdt;
	}
	break;
      case i_matrix_type:				// Multiply n_matrix into i_matrix adjoint
	return this;					// No change so return original matrix
	break;
      case h_matrix_type:				// Multiply n_matrix into h_matrix adjoint
	return multiply((h_matrix*)mx);			// Since hmx is self adjoint, just use multiply
	break;
      default:						// Multiply n_matrix into generic matrix adjoint
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
          cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, A_rows, B_cols, A_cols, 
                             alpha, (double *)data, A_cols, (double *)((n_matrix*)mx)->data, B_cols, beta, (double *)pdt->data, C_cols);
//        std::cerr << "BLAS: n_matrix * unknown\n";
        }
        else
        {
	int pos = 0;					// Index for pdt, set to 0 = <0|pdt|0>
	for(int i=0; i<rows_; i++)			// Loop over the rows of pdt
	  for(int j=0; j<pdt->cols(); j++,pos++)	// Loop over the columns of pdt
	    for(int k=0; k<cols_; k++)			// Loop over the index k
	      pdt->data[pos] += 			// <i|pdt|j> += <i|nmx|k><j|mx*|k>
                 conj((*mx)(j,k),(*this)(i,k));

	}
#else
	n_matrix* pdt =					// Construct new n_matrix, filled with 0's
          new n_matrix(rows_,mx->cols(),complex0);
	int pos = 0;					// Index for pdt, set to 0 = <0|pdt|0>
	for(int i=0; i<rows_; i++)			// Loop over the rows of pdt
	  for(int j=0; j<pdt->cols(); j++,pos++)	// Loop over the columns of pdt
	    for(int k=0; k<cols_; k++)			// Loop over the index k
	      pdt->data[pos] += 			// <i|pdt|j> += <i|nmx|k><j|mx*|k>
                 conj((*mx)(j,k),(*this)(i,k));
#endif
	return pdt;
	}
      }
  }


// ___________________________________________________________________________
// J                  CLASS N_MATRIX COMPLEX UNARY FUNCTIONS
// ___________________________________________________________________________

/*  Function/Operator                         Result
    -----------------   ------------------------------------------------------
         det            Determinant of the array (NOT IMPLEMENTED)
         rank           Rank of the array (NOT IMPLEMENTED)                 */

complex n_matrix::det()  { NMxfatal(25,"det");  return complex0; }
int     n_matrix::rank() { NMxfatal(25,"rank"); return 0; }


// ********************** Discrete Fourier Transform ******************** 

void n_matrix::FFT(int isign, bool comp)

        // Input          nmx   : An n_matrix (this)
        //                isign : Flag for transform or inverse transform
        //                                1=FFT, -1=IFFT
	//		  comp  : FFT compatibility flag.  If true the
	//			  transform is compatible with both Brigham
	//			  and MATLAB in that exp(iwt) produces a
	//			  positive peak at w upon transformation
	//			  If not true this will come at -w.
        // Output         nmx . : Modifies this to contain either its
        //                        discrete Fourier transform, or its
        //                        inverse discrete Fourier transform
	// Note			: The FFT alogorithm is from Numerical Recipies,
	// 			  Chapter 12, Section 2
	// Note			: The matrix is herein assumed to be a linear
	//			  array & transformed as a whole.  It does not
	//			  perform an FFT row-by-row (or column-by-column)
	// Note			: The matrix dimension must be 2*nn where nn is
	//			  an integer power of 2, check before entering!

  {
  complex w, wp, temp;
  double theta, sto2;
  int i,j,m,mmax,istep,ii;
  if(comp)				// Some FFT's are defined differently
    for(i=0; i<size; i++)		// and to account for the differences
      data[i] = conj(data[i]);		// we need to change sign of imaginaries
  if(isign<0)				// This for the inverse FFT  
    for(i=0,j=size/2; j<size; i++,j++)	// we must swap the the data
      {					// about the center point
      temp = data[j];			// [0,1,2,....N/2,....N-1]
      data[j] = data[i];		// becomes
      data[i] = temp;			// [N/2,...N-1,0,1,...N/2-1]   
      }
  for(j=0,i=0; i<size; i++)		// Bit reversal section
    {
    if(j>i)
      {
      temp = data[j];
      data[j] = data[i];
      data[i] = temp;
      }
    m = size/2;
    while((m>1) && (j>=m))
      {
      j = j-m;
      m = m/2;
      }
    j=j+m;
    }
  mmax=1;
  while (size>mmax)			// Begin the Danielson-Lanczos section
    {					// Outer loop executes log (nn) times
    istep = 2*mmax;			//                            2
    theta = PI/(isign*mmax);
    sto2 = sin(theta/2);
    wp = complex(-2*sto2*sto2, sin(theta));
//    wp = complex(-2*sqr(sin(theta/2)), sin(theta));	// sqr in libg++ builtin!
    w=1;
    for(ii=0; ii<mmax; ii++)
      {
      for(i=ii; i<size; i+=istep)
	{
	j = i+mmax;			// This is the Danielson Lanczos formula
	temp = w*data[j];
	data[j] = data[i] - temp;
	data[i] += temp;
	}
      w+=w*wp;
      }
    mmax *= 2;
    }
  if(isign>0)  
    for(i=0,j=size/2; j<size; i++,j++)	// Bit reversal section
      {
      temp = data[j];
      data[j] = data[i];
      data[i] = temp;
      }
  unitary = false;
  }


// ____________________________________________________________________________
// K            CLASS N_MATRIX COMPLEX BINARY (&TRINARY) FUNCTIONS
// ____________________________________________________________________________

// ***************************** Tensor Product *******************************

_matrix* n_matrix::tensor_product(_matrix* mx) 

        // Input            nmx : An n_matrix (this)
	//		     mx : A second matrix
        // Output           pdt : A matrix which is the tensor product
        //                        of the two input matrices

	//			     pdt        =   nmx (x) mx

	//			 (m*o x n*p)       (mxn)   (oxp)

	//		      <i*o+k|pdt|j*p+l> = <i|nmx|j><k|mx|l>
  {
  switch (mx->stored_type())
    {
    case n_matrix_type:				// Tensor Product n_matrix (x) n_matrix
      {
      int ca = cols_;				// Columns of nmx
      int ra = rows_;				// Rows of nmx
      int cb = mx->cols();			// Columns of mx
      int rb = mx->rows();			// Rows of mx
      n_matrix* pdt = new n_matrix(ra*rb,ca*cb);// New n_matix, pdt, much larger
      complex *p0000 = pdt->data;		// Start of data in pdt
      complex *n00 = data;			// Start of data in nmx
      complex *m00 = ((n_matrix*)mx)->data; 	// Start of data in mx
      complex *nend = n00 + rows_*cols_;	// End of data in nmx
      complex *mend = m00 + rb*cb;		// End of data in mx
      complex *nirowend = NULL;			// Index to track row end in nmx
      complex *mkrowend;			// Index to track row end in mx
      complex *pikjl,*nij,*mkl,*ni0,*mk0;
      for(pikjl=p0000, ni0=n00; ni0<nend;	// Effective loop over i
                                 ni0=nirowend)
	{
	nirowend = ni0 + cols_;			// Determine where row i of nmx ends 
	for(mk0=m00; mk0<mend; mk0=mkrowend)	// Effective loop over k, <k|mx|*>
	  {
	  mkrowend = mk0 + cb;			// Determine where row k of mx ends
	  for(nij=ni0; nij<nirowend; nij++)	// Effective loop over j, <i|nmx|j>
	    for(mkl=mk0; mkl<mkrowend; 		// Effective loop over l, <k|mx|l>
	                        mkl++,pikjl++)
              *pikjl = (*nij) * (*mkl);		// <ik|pdt|jl> = <i*cb+k|pdt|j*cb+l> = <i|nmx|j><k|mx|l>
	  }
	}
      return pdt;
      }
      break;
    case h_matrix_type:				// Tensor Product n_matrix (x) h_matrix
      { 
      int hd = mx->cols();			// Column (and row) dimension of hmx
      int pr = rows_*hd;			// Rows in tensor product
      int pc = cols_*hd;			// Columns in tensor product
      n_matrix* pdt = new n_matrix(pr, pc); 	// Construct new matrix for product
      complex *pikjk = pdt->data;		// Start of pdt: <ik|pdt|jl> -> <00|pdt|00>
      complex *n00 = data;			// Start of data in nmx
      complex *h00 = ((h_matrix*)mx)->data; 	// Start of hmx: <0|hmx|0>
      complex *pend = pikjk + pr*pc;		// End of data in pdt: <pr|pdt|pc+1>
      complex *nij;
      complex *h0k, *hkl, *hlk, *hkk, *HKK;
      complex *ni0 = n00;			// Start of nmx: <i|nmx|0> -> <0|nmx|0>
      complex *nip10 = ni0 + cols_;		// <i+1|nmx|0> -> <1|nmx|0>
      complex *h11 = h00 + hd;			// <1|hmx|1>
      int k, hrowl;
      for(; pikjk<pend; ni0+=cols_,nip10+=cols_)// Effective loop over i (nmx rows)
	for(k=0,h0k=h00; h0k<h11; k++,h0k++)	// Loop through all the rows of hmx
//	for(h0k=h00,HKK=h00; h0k<h11;		// Loop through all the rows of hmx
//	              h0k++,HKK+=hrowk,hrowk--)	// Loop through all the rows of hmx
          {
          HKK = h00 + k*hd-(k*(k-1))/2;		// Set <k|hmx|k>
	  for(nij=ni0; nij<nip10; nij++)	// Effective loop over j (nmx columns)
            {
            hrowl = hd;				// Row length of hmx, start at row 0
	    for(hlk=h0k,hkk=HKK; hlk<hkk;	// Effective loop over l (hmx columns), l<k
                    pikjk++,hrowl--,hlk+=hrowl)	// where hrowl is length of row l of hmx
	      *pikjk =(*nij) * conj(*hlk);	// <i*hd+k|pdt|j*hd+l> = <i|nmx|j>*<l|hmx*|k> 
            hkk += hrowl;			// <k|hmx|k> --> <k+1|hmx|k+1>
            for(hkl=hlk; hkl<hkk; pikjk++,hkl++)// Continue loop over l (hmx cols), k<=l
	      *pikjk = (*nij) * (*hkl);		// <i*hd+k|pdt|j*hd+l> = <i|nmx|j>*<k|hmx|l> 
            }
          }
      return pdt;
      }
      break;
    case d_matrix_type:				// Tensor Product n_matrix (x) d_matrix
      {
      int ca = cols_;				// Columns of nmx
      int ra = rows_;				// Rows of nmx
      int cb = mx->cols();			// Columns (and rows) of mx
      int cm = ca*cb;				// Columns in tensor product
      int rm = ra*cb;				// Rows in tensor product
      n_matrix* pdt =				// Construct new matrix for product
                   new n_matrix(rm,cm,complex0);// which has all elements of 0
      complex *p0000 = pdt->data;		// Start of data in pdt
      complex *n00 = data;			// Start of data in nmx
      complex *d00 = ((d_matrix*)mx)->data; 	// Start of data in dmx
      complex *nend = n00 + rows_*cols_;	// End of data in nmx
      complex *nrow = NULL;			// Dummy index for looping over 1 row in nmx
      complex *dend = d00 + cb;			// End of data in dmx
      complex *pikjk,*nij,*dkk,*pikj0,*pik00;
      for(pik00=p0000,nij=n00; nij<nend; 	// Effective i loop of <i|nmx|j>, pik00 = <i*cb+k|pdt|0>
                                  pik00+=cm*cb)
        for(nrow=nij+ca,pikj0=pik00; nij<nrow;	// Effective j loop of <i|nmx|j>, pikj0 = <i*cb+k|pdt|j*cb>
                               nij++,pikj0+=cb)
          for(dkk=d00,pikjk=pikj0; dkk<dend; 	// Effective loop over k, k = [0, cb)
                             dkk++,pikjk+=cm+1)
            *pikjk = (*nij) * (*dkk); 		// <ik|pdt|jl> = <i*cb+k|pdt|j*cb+l> = <i|nmx|j><k|dmx|l>
      return pdt;				//             = del   * <i|nmx|j><k|dmx|k>
      }						//                  k,l
      break;
    case i_matrix_type:				// Tensor Product n_matrix (x) i_matrix
      {
      int ca = cols_;				// Columns of nmx
      int cb = mx->cols();			// Columns (and rows) of mx
      int cm = cols_*cb;			// Columns in tensor product
      int rm = rows_*cb;			// Rows in tensor product
      n_matrix* pdt =				// Construct new matrix for product
                   new n_matrix(rm,cm,complex0);// which has all elements of 0
      complex *p0000 = pdt->data;		// Start of data in pdt
      complex *n00 = data;			// Start of data in nmx
      complex *nend = n00 + rows_*cols_;	// End of data in nmx
      complex *nrow = NULL;			// Dummy index for looping over 1 row in nmx
      complex *pikjk, *pikj0;			// Elements <i*cb+k|pdt|j*cb+k>, <i*cb+k|pdt|j*cb>
      complex *nij, *pik00;			// Elements <i|nmx|j> and <i*cb+k|pdt|0>
      int k;
      for(pik00=p0000,nij=n00; nij<nend; 	// Effective i loop of <i|nmx|j>, pik00 = <i*cb+k|pdt|0>
                                  pik00+=cm*cb)
        for(nrow=nij+ca,pikj0=pik00; nij<nrow;	// Effective j loop of <i|nmx|j>, pikj0 = <i*cb+k|pdt|j*cb>
                               nij++,pikj0+=cb)
          for(k=0,pikjk=pikj0; k<cb; 		// Effective loop over k, k = [0, cb)
          		       k++,pikjk+=cm+1)
            *pikjk = *nij; 			// <ik|pdt|jl> = <i*cb+k|pdt|j*cb+l> = <i|nmx|j><k|mx|l>
      return pdt;				//             = del   * <i|nmx|j> = del   <ik|pdt|jk>
      }						//                  k,l			k,l
      break;      
    default:					// Tensor Product n_matrix (x) generic matrix
      { 
      int rg = mx->rows();			// Rows of gmx
      int cg = mx->cols();			// Columns of gmx
      int rp = rows_*rg;			// Rows in tensor product
      int cp = cols_*cg;			// Columns in tensor product
      n_matrix* pdt = new n_matrix(rp,cp);	// Construct new matrix for product
      complex *p0 = pdt->data;			// Start of data in pdt: <00|pdt|00>
      complex *n00 = data;			// Start of data in nmx: <0|nmx|0>
      complex *nend = n00 + rows_*cols_;	// End of data in nmx: <rows_|nmx|cols+1>
      complex *pikjl, *nij, *ni0;		// <i*rg+k|pdt|j*cg+k>, <i|nmx|j>, <i|nmx|0>
      int k,l;
      for(pikjl=p0,ni0=n00;ni0<nend; ni0+=cols_)// Effective loop over i (rows of nmx)
	for(k=0; k<rg; k++)			// Loop over index k (rows of mx)
	  for(nij=ni0; nij<ni0+cols_; nij++)	// Effective loop over j (columns of nmx)
	    for(l=0; l<cg; l++,pikjl++)		// Loop over index l (columns of mx)
	      *pikjl = (*nij) * (*mx)(k,l);	// <i*rb+k|pdt|j*cb+l> = <i|nmx|j>*<k|mx|l> 
      return pdt;
      }
    }
  }


// ____________________________________________________________________________
// L                      CLASS N_MATRIX I/O FUNCTIONS
// ____________________________________________________________________________

/* These functions perform both ASCII and Binary input/output operations on a
   complex matrix.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   print      ostream     Writes array elements in ASCII to output stream
   picture    ostream     Writes array elements pictorically to output stream
   write      ofstream    Writes array elements in Binary to output file stream
   read       ifstream    Read in array from Binary input file stream
   readASC    istream     Read in array from an input stream
   ask          ----      Interactively requests complex matrix from user 
    
   The binary format of n_matrix has data only, not size information.  The size 
   is taken care of in classes matrix/_matrix and should done prior to use of 
   these functions.  The data ordering is Re(<i|hmx|j>, Im(<i|hmx|j> columns 
   then rows (i.e. row by row.)                                              */ 

std::vector<std::string> n_matrix::printStrings(const MxPrint& PFlgs) const
  {
  std::vector<std::string> PStrings;
  int ptype = 0;                                // Complex elements output
  if(PFlgs.MxRIPrnt)                               // If we want just reals
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
  std::string pline;
  std::string efmt = complex::dformat();        // Real/Imag element format
  int i,j,pos;
  for(i=0,pos=0; i<rows_; i++ )                 // Loop array rows
    {
    pline = std::string("");
    if(ptype==1)				// Output real elems only
      for(j=0; j<cols_; j++) 
        pline += MxModform(efmt.c_str(),Re(data[pos++]));
    else if(ptype==2)                           // Output imag elems only
      for(j=0; j<cols_; j++) 
        pline += MxModform(efmt.c_str(),Im(data[pos++]));
    else                                        // Output complex elems
      for(j=0; j<cols_; j++)
        pline += (data[pos++]).printString();
    PStrings.push_back(pline);				// Go to next line
    }
  return PStrings;
  }

std::vector<std::string> n_matrix::pictureStrings(const MxPrint& PFlgs) const
  {
  std::vector<std::string> PStrings;
//  int rlen = 2*rows_-1;                         // Length of 1 row
  std::string pline;
  int i, j, pos;                                // Looping indices
  for(i=0, pos=0; i<rows_; i++ )                // Loop array rows
    {
    pline = "";
    for(j=0; j<cols_; j++, pos++)		// Now fill in cols
      {
      if(norm(data[pos])) pline += "x";
      else                pline += "0";
      if(j+1 < cols_)     pline += " ";
      }
    PStrings.push_back(pline);
    }
  return PStrings;
  }

void n_matrix::print(std::ostream& ostr, const MxPrint& PF) const
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
 
  std::string efmt = complex::dformat();             // Real/Imag element format
  int clen = 40 - ((elen+1)*cols_-1)/2;         // Space to center 1 line
  std::string sp("");                                // Spacer to center a line
  if(clen>0) sp = std::string(clen, ' ');            // Set spacer for centering
  int i,j,pos;
  for(i=0,pos=0; i<rows_; i++ )                 // Loop array rows
    {
    ostr << sp;                                 // Space to center
    if(ptype==1)				// Output real elements
      for(j=0; j<cols_; j++) 
        ostr << MxModform(efmt.c_str(),Re(data[pos++]));
    else if(ptype==2)                           // Output imag elements
      for(j=0; j<cols_; j++) 
        ostr << MxModform(efmt.c_str(),Im(data[pos++]));
    else                                        // Output complex elems
      for(j=0; j<cols_; j++) ostr << data[pos++];
    ostr << "\n";				// Go to next line
    }
  }

void n_matrix::picture(std::ostream& ostr, const MxPrint& PF) const
  {
 int rlen = 2*rows_-1;                         // Length of 1 row
  std::string sp("");                                // Spacer to center row
  int len = 40-rlen/2;                          // Space to center row
  if(len>0) sp = std::string(len, ' ');              // Set spacer to center
  int i, j, pos;                                // Looping indices
  for(i=0, pos=0; i<rows_; i++ )                // Loop array rows
    {
    ostr << sp;                                 // Center row
    for(j=0; j<cols_; j++, pos++)		// Now fill in cols
      {
      if(norm(data[pos])) ostr << "x ";
      else                ostr << "0 ";
      }
    ostr << "\n";                               // Skip to next row
    }
  }

void n_matrix::write(std::ofstream &fp, int form)
  {
	// *** converted from float d to double d 
	// and from sizeof(float) to sizeof(double)
  double d;
  for(int i=0, pos=0; i<rows(); i++)	// First loop over the matrix rows
    for(int j=0; j<cols(); j++,pos++)	// Now loop over the matrix columns
      {
      d = Re(data[pos]);
      fp.write((char*)&d, sizeof(double));
      d = Im(data[pos]);
      fp.write((char*)&d, sizeof(double));
      }
  return;
  form = 0;				// Compiler likes form to be used
  }

void n_matrix::read(std::ifstream &fp)
  {
  float dr,di;
  for(int i=0, pos=0; i<rows(); i++ )
    for(int j=0; j<cols(); j++,pos++)
      {
      fp.read((char*)&dr, sizeof(float));
      fp.read((char*)&di, sizeof(float));
      data[pos] = complex(dr,di);
      }
  }

void n_matrix::readASC(std::istream& istr)
  {
  int i,j;
  istr >> i >> j;
  resize(i,j);
  for(int pos=0; pos<size; pos++)
    istr >> data[pos]; 
  }


void n_matrix::ask( )
  {
  float dr,di;
  for(int i=0, pos=0; i<rows(); i++ )
    for(int j=0; j<cols(); j++,pos++)
      {
      std::cout << "\n\tPlease Input Real and Imaginary Value of <"
           << i << "|mx|" << j << "> [re im]: ";
      std::cin >> dr >> di;
      data[pos] = complex(dr,di);
      }
  }


// ____________________________________________________________________________
// M               CLASS N_MATRIX MISCELLANEOUS FUNCTIONS
// ____________________________________________________________________________


void n_matrix::resize(int i, int j)

        // Input          nmx   : A n_matrix (this)
        //                i,j   : Row, column dimensions
        // Output         none  : Size of nmx readjusted to ixj
        // Note                 : Destroys data if not the same size,
        //                        no alterations if same size

  {
  _matrix::resize(i,j);			// Reset row & col stored values
  if((i*j) != size)			// Reset size and re-allocate
    {					// data if # elements changes 
    delete [] data;
    size = i*j;
    data = new complex[size];
    }
  }


_matrix* n_matrix::copy()

	// Input                this : Normal matrix (_matrix)
	// Output               mx   : Pointer to normal matrix (_matrix*)
	// Note                      : This returns a new n_matrix on exit.
	//                             Unlike class matrix where assignment
	//			       involves only referencing, a new copy
	//			       of the matrix data is done here

  { return new n_matrix(*this); }


void n_matrix::convert(_matrix* mx)
 
        // Input          nmx   : A n_matrix (this)
        //                mx    : A matrix of any type
        // Output         mx    : Elements of nmx are placed into mx
        //                        in a fashion so as to maintain the mx
        //                        matrix type.  The result is the "best"
        //                        matrix conversion of nmx will reside in mx
        // Note                 : Since a normal matrix cannot always be converted to
        //                        other matrix types, many conversions may result
	//			  in loss of data.
        // Note                 : Any existing elements in mx will vanish
        // Note                 : Uses internal structure of other matrix types
        // Note                 : The use of "resize" allocates a fresh data array

  {
  switch (mx->stored_type())
    {
    case n_matrix_type:				// Conversion of nmx into normal matrix
      (*(n_matrix *)mx) = (*this);		// is unnecessary
      break;
    case h_matrix_type:				// Conversion of nmx into hermitian matrix
      {
      if(rows_ != cols_)			// Insure that nmx is square
        {
        NMxerror(14,1);				//    Bad rectangular array use
        NMxfatal(70);				//    Cannot convert to Hermitian
        }
      mx->resize(rows_,cols_);			// Now set Hermitian matrix to right size
      int nd = rows_;                           // Dimension of square nmx
      complex *nij = data;			// Start of data in nmx: <i|nmx|j>-><0|nmx|0>
      complex *hij = ((h_matrix*)mx)->data;     // Start of data in hmx: <i|hmx|j>-><0|hmx|0>
      complex *hend = hij + (nd*nd+nd)/2;       // End of data in hmx: <nd|hmx|nd+1>
      complex *nip10 = nij + cols_;		// Element <i+1|nmx|0>
      int i = 0;				// Row index; <i|nmx|0> + i -> <i|nmx|i>
      for( ; hij<hend; i++,nip10+=cols_,nij+=i)	// Copy only <i|nmx|j> elements with i<=j
        {					// and thus, forget any lower triangle nmx elements
	(*hij) = Re(*nij);			// <i|hmx|i> = Re(<i|nmx|i>) j>i
        for(hij++,nij++; nij<nip10;hij++,nij++)
	  (*hij) = (*nij);			// <i|hmx|j> = <i|nmx|j> j>i
        }
      }
      break;
    case d_matrix_type:				// Conversion of nmx into diagonal matrix
      {
      if(rows_ != cols_)			// Insure that nmx is square
        {
        NMxerror(14,1);				//    Bad rectangular array use
        NMxfatal(71);				//    Cannot convert to diagonal
        }
      mx->resize(rows_,cols_);			// Now set diagonal matrix to right size
      complex *nii = data;			// Start of data in nmx: <i|nmx|j>-><0|nmx|0>
      complex *dii = ((d_matrix*)mx)->data;	// Start of data in dmx: <i|dmx|i>-><0|dmx|0>
      complex *dend = dii + rows_;		// End of data in dmx: <rows_|nmx|rows_+1>
      int nrow = rows_+1;			// Amount to move <i|nmx|i> to <i+1|nmx|i+1>
      for(; dii<dend; dii++,nii+=nrow)		// Loop over all the elements of dmx
	(*dii) = (*nii);			// <i|dmx|i> = <i|nmx|i>, forget other <i|nmx|j>
      }
      break;
    case i_matrix_type:				// Conversion of nmx into identity matrix
      if(rows_ != cols_)			// Insure that nmx is square
        {
        NMxerror(14,1);				//    Bad rectangular array use
        NMxfatal(72);				//    Cannot convert to identity
        }
      mx->resize(rows_, cols_);			// Set identity matrix size, loose <i|nmx|j>
      break;
    default:					// Conversion of nmx into an unknown matrix
      mx->resize(rows_,cols_);			// Set/Check matrix size
      complex *nij = data;			// Element <i|nmx|j> -> <0|nmx|0>
      int i,j;
      for(i=0; i<rows_; i++)			// and then copy all the elements <i|nmx|j> 
	for(j=0; j<cols_; j++,nij++)		// to the elements of <i|mx|j>
	  (*mx)(i,j) = (*nij);			// <i|mx|j> = <i|nmx|j>
    }
  }

//
// ------------------------ New Conversion Routines ---------------------------
//                  (Been Having Troubles With The Old One)

 
i_matrix* n_matrix::IMX() { return new i_matrix(cols_, cols_); }
 
        // Input                nmx     : A n_matrix (this)
        // Output               imx     : A new i_matrix matrix whose
        //                                elements are the same as nmx
        // Note                         : This will allocate NEW memory
        // Note                         : This will result in a total
        //                                loss of nmx elements
 
 
 
        // Input                nmx     : A n_matrix (this)
        // Output               dmx     : A new d_matrix matrix whose
        //                                elements are the same as nmx
        // Note                         : This will allocate NEW memory
        // Note                         : Array nmx MUST be square or fatal
	//				  error will occur herein
d_matrix* n_matrix::DMX()
  {                                             // We then loose nmx off diagonals
  if(rows_ != cols_)				// Array nmx MUST be square
    {
    NMxerror(14,1);				//    Bad rectangular array use
    NMxfatal(71);				//    Cannot convert to Hermitian
    }
  d_matrix* dmx = new d_matrix(rows_,rows_);	// Start with an empty d_matrix
  complex *nii = data;				// Start of data in nmx: <i|nmx|i>-><0|nmx|0>
  complex *dii = ((d_matrix*)dmx)->data;	// Start of data in dmx: <i|dmx|i>-><0|dmx|0>
  complex *dend = dii + rows_;			// End of data in dmx: <rows_|nmx|rows_+1>
  for(; dii<dend; dii++,nii+=rows_+1)		// Loop over all the elements of dmx
    (*dii) = (*nii);				// <i|dmx|i> = <i|nmx|i>, forget other <i|nmx|j>
  return dmx;
  }
 
 
                             
        // Input                nmx     : A n_matrix (this)
        // Output               hmx     : A new h_matrix matrix whose
        //                                elements are the same as nmx
        // Note                         : This will allocate NEW memory
        // Note                         : This may result in a partial
        //                                loss of nmx "non-hermitian" elements

h_matrix* n_matrix::HMX()
  {
  if(rows_ != cols_) 				// Array nmx MUST be square
    {
    NMxerror(14,1);				//    Bad rectangular array use
    NMxfatal(70);				//    Cannot convert to Hermitian
    }
  h_matrix* hmx = new h_matrix(rows_,rows_);	// Start with an empty h_matrix
  complex *nij = data;				// Start of data in nmx: <i|nmx|j>-><0|nmx|0>
  complex *hij = ((h_matrix*)hmx)->data;	// Start of data in hmx: <i|hmx|j>-><0|hmx|0>
  int i,j;
  for(i=0; i<rows_; i++, nij+=i)
    {
    (*hij) = zRe(*nij);				// <i|hmx|i> = Re<i|nmx|i>
    hij++;
    nij++;
    for(j=i+1; j<rows_; j++, hij++, nij++)
      (*hij) = (*nij);				// <i|hmx|j> = <i|nmx|j>, upper triangle
    }
  return hmx;
  }


n_matrix* n_matrix::NMX() { return this; }

        // Input                nmx     : A n_matrix (this)
        // Output               nmx     : A new n_matrix matrix whose
        //                                elements are the same as nmx
        // Note                         : This will NOT allocate NEW memory


 
// ____________________________________________________________________________
// N                CLASS N_MATRIX DIAGONALIZATION FUNCTIONS
// ____________________________________________________________________________

// Note: Due to the length of the diagonalization algorithm and associated
//       functions, they reside in the file n_matrix_diag.cc.  They are listed
//       below but commented out.

// void n_matrix::diag(_matrix *&mx, _matrix *&mx1)
// void n_matrix::cred(n_matrix& z)
// void n_matrix::rred(n_matrix& w)
// void n_matrix::tqli(n_matrix& z, d_matrix& a)
// complex* n_matrix::corth(int low, int igh)
// int n_matrix::comqr3(int low,int igh, complex *ort, d_matrix& w,
//				                 n_matrix& z, int flag)

// ____________________________________________________________________________
// O                    CLASS N_MATRIX INVERSION FUNCTIONS
// ____________________________________________________________________________

// Note: Due to the length of the inversion algorithm and associated functions,
//       some reside in the file n_matrix_inv.cc.  These are listed below but
//       commented out.

_matrix* n_matrix::inv()

        // Input            nmx : An n_matrix (this)
	// Output          ninv : Inverse of nmx as defined by
	//			         nmx * ninv = I
	// Note			: The algorithm solves for the inverse
	//			  using an LU decomposition of hmx as per
	//			      nmx * ninv = LU * ninv = I
	//			  by first generating I and LU then back
	//			  substituting to produce ninv

  {
  if(unitary) return adjoint();			// This is real new (sosi)
  int nd = rows_;				// Matrix dimension
  if(nd != cols_)
    {
    NMxerror(14,1);				//    Bad rectangular array use
    NMxfatal(25);				//    Cannot take inverse
    }
  int *indx;
  indx = new int[nd];
  n_matrix* nLU = new n_matrix(*this);          // Copy input array for LU decomp
  nLU->LU_decomp(indx);                         // LU decomposition of nLU
  n_matrix* I = new n_matrix(nd,nd,complex0);   // Constuct a new n_matrix
  for(int i=0; i<nd; i++) I->put(complex1,i,i); // Make it the identity matrix
  _matrix* mxinv = I->LUinv(indx,nLU);          // Back solve for inverse
  delete nLU;                                   // Clean up the nLU matrix
  delete I;                                     // Clean up the I matrix
  delete [] indx;
  return mxinv;                                 // Return the inverse of nmx
  }


        // Input            nmx : An n_matrix (this)
	//		    indx: Row permutation array
	// Output           nLU : LU decomposition of a row permutation
	//			  of the input matrix nmx, nmx', the row
	//			  exchanges recorded in array indx
	//			  nLU(low.tri.) * nLU(upp.tri.) = nmx'
	//			         L      *       U       = nmx'
	// Note			: The returned array has L on it's lower
	//			  triangle and U on its upper triangle,
	//			  the latter containing the diagonals
	//			  All <i|L|i> = 1, and are absent in nLU

_matrix* n_matrix::LU(int *indx)
  {
  n_matrix* nLU = new n_matrix(*this);	// New n_matrix equal to this, start nLU
  nLU->LU_decomp(indx);			// LU Decomposition of row permuted nmx
  nLU->unitary = false;			// No longer unitary if ever was
  return nLU;				// Return Ainv, a new n_matrix
  }


_matrix* n_matrix::LUinv(int *indx, _matrix* LU)

        // Input            B   : An n_matrix (this) of results
	//		    LU  : An LU decomposition of matrix A
	//			  (or row permutation of matrix A)
	//		    indx: Row permutation array
	// Output           X   : Solution to the equation: A|X> = |B>
	//			  where LU is the LU decomposition of A
	//			  or row permutation of A
	// Note		        : Matrix LU must be square
	// Note		        : Uses internal structure of other matrix types

  {
  int dlu = LU->rows();				// Dimension of LU matrix
  if(LU->cols() != dlu) 			// Insure LU is square
    {
    NMxerror(14,1);				//    Bad rectangular array use
    NMxfatal(26);				//    Problems with LU decomposition
    }
  if(dlu != rows_) 				// Dimension check: nLU(ixj)*x(jxk) = b(jxk)
    {
    NMxerror(51,1);				//    Array dimension mismatch
    NMxfatal(80);				//    Problems with LU decomposition
    }
  n_matrix* X = new n_matrix(*this);		// New n_matrix for X, X=B
  int rswap = (indx[0] >= 0);			// Flag if row swapping needed
  int i, ip;
  switch(LU->stored_type())
    {
    case i_matrix_type:				// Identity LU Matrix - LU*X = I*X = B
      {
      if(!rswap)				// If no row-permutations, then
        return this;				// LU*X = I*X = X = B , return B
      complex z;				// Temporary complex for row swapping				
      complex* x0k = X->data;			// Start of data in X: <0|X|k> -> <0|X|0>
      complex *xik;				// Element of X: <i|X|k>
      for(int k=0; k<cols_; k++,x0k++)		// Loop over k (X,B columns) 
        {
//	    First Solve For Column k of Y: L|Y> = I|Y> = |B> For |Y>
//		        (Just Perform Any Row Permutations!)

        for(xik=x0k,i=0;i<dlu;i++,xik+=cols_)	// Works forwards, |y(0)> is solved first
          { 					// because L is lower triangular
          ip = indx[i];				// Row permutations exist, so row is
          z = *xik; 				// swapped with row ip
          *xik = X->data[ip*cols_+k];		// So <i|B|k> is actually <ip|B|k>
          X->data[ip*cols_+k] = z;		// Now swap counterpart: <ip|B|k> <- <i|B|k>
          }

//	     Now Solve For Column k of X: U|X> = I|X> = |Y> For |X>
//   				(Don't Do Anything!)
        }
      return X;					// Return the results X
      }
      break;
    case d_matrix_type:				// Diagonal LU Matrix - LU*X = D*X = B
      {
      complex z;				// Temporary complex for row swapping				
      complex *lu00 = ((d_matrix*)LU)->data; 	// Start of data in LU matrix <0|LU|0>
      complex* luend = lu00 + dlu;		// End of data in LU: <rows_|LU|rows_+1>
      complex* x0k = X->data;			// Start with element <0|X|k> -> <0|X|0>
      complex *xik, *luii;			// Actual Elements: <i|X|k> and <i|LU|i>
      for(int k=0; k<cols_; k++,x0k++)		// Loop over k (X,B columns) 
        {					// and solve for the columns of X.

//	    First Solve For Column k of Y: L|Y> = D|Y> = |B> For |Y>
//		        (Just Perform Any Row Permutations!)

        xik = x0k;				// Set <i|X|k> to <0|X|k>, i.e. i=0
        for(i=0; i<dlu; i++,xik+=cols_)		// Works forwards, |y(0)> is solved first
          { 					// because L is lower triangular
          if(rswap)				// Swap if row permutations exist
            {					// (indx[0]<0 -> no permutations)
            ip = indx[i];			// Row permutations, row i is swapped with row ip
            z = *xik;
            *xik = X->data[ip*cols_+k];		// So <i|B|k> is actually <ip|B|k>
            X->data[ip*cols_+k] = z;		// Now swap counterpart: <ip|B|k> <- <i|B|k>
            }
          else
            *xik = X->data[i*cols_+k];		// No permutations, <i|Y|k> = <i|B|k>
          }

//	     Now Solve For Column k of X: U|X> = D|X = |Y> For |X>

//	               [	   --- 	             ]
//	          1    |           \ 		     |	  <i|Y|k>
//   <i|X|k> = ------- | <i|Y|k> - /  <i|U|j><j|X|k> | =  --------
//	       <i|U|i> |           --- 		     |    <i|LU|i>
//	               [           j>i 		     ]

        xik = x0k;				// Set <i|X|k> to <0|X|k>, i.e. i=0
	luii = lu00;				// Set <i|LU|i> to <0|LU|0>
        for(; luii<luend; xik+=cols_,luii++)	// Loop over all rows of column k of X
          *xik /= (*luii);			// This is now <i|X|k>
        }
      return X;					// Return the results X
      }
      break;
    case h_matrix_type:				// Hermitian LU Matrix (Extremely Unlikely!)
      {
      complex *lu00 = ((h_matrix*)LU)->data;	// Start of data in LU matrix
      complex *luend = lu00 + (dlu*dlu+dlu)/2;	// End of data in LU: <dlu|LU|dlu+1>
      complex *luii, *luij;
      complex* x00 = X->data;			// Start of data in X: <0|X|k>
      complex* xend = x00 + rows_*cols_;	// End of data in x: <rows_|LU|cols_+1>
      complex *xjk = x00;			// Element of X: <j|X|k> -> <0|X|0>
      complex *xik, *x0k;			// Element of X: <i|X|k>
      complex *luji, *lu0i;
      int lurow = dlu+1;			// Amount to move <i|LU|i> to <i+1|LU|i+1> 
      int lrow = dlu;
      complex sum = 0;
      x0k = x00;
      int ii;
      for(int k=0; k<cols_; k++,sum=0,x0k++)	// Loop over k (X,B columns)
        {					// and solve for column X|k> row-wise.

/*	        First Solve For Column k of Y: L|Y> = |B> For |Y>

  		         --- 	                       ---
  		         \ 			       \
     <i|Y|k> = <i|B|k> - /  <i|L|j><j|Y|k> = <i|B|k> - /  <j|LU*|i><j|Y|k> 
  		         --- 		               ---
  		         j<i 			       j<i                   */

        xik = x0k;				// Set <i|X|k> to <0|X|k>, i.e. i=0
        i = 0;
        for(ii=0,lu0i=lu00; i<dlu;	 	// Works forwards, <0|Y|k> is solved first
                        i++,xik+=cols_,lu0i++)	// because L is lower triangular
          {
          if(rswap)				// Must swap if there has been row permutations
            {					// (negative indx[0] indicates no permutations)
            ip = indx[i];			// Row permutations, row i is swapped with row ip
            sum = X->data[ip*cols_+k];		// Start with actual <i|B|k>, <ip|B|k> for sum
            X->data[ip*cols_+k] = *xik;		// Now swap row for next iteration: <ip|B|k> <- <i|B|k>
            }
          else
            sum = X->data[i*cols_+k];		// Start with actual <i|B|k>, <ip|B|k> for sum
          if(ii)
            {
            lrow = dlu;				// Amount to move from <j|LU|i> -> <j+1|LU|i>
            for(xjk=x0k,luji=lu0i; xjk<xik;	// Effective loop over j (LU columns), j<i
                 xjk+=cols_,lrow--,luji+=lrow)	//  (lower-tri of LU == upper-tri of LU*)
              sum -= conj(*luji) * (*xjk);	// sum -= <j|LU*|i>*<j|B|k>, i>j
            }
          else if(norm(sum))
            ii=1;
          *xik = sum;				// Now set <i|Y|k> to sum (stored in X)
          }
    
//	     Now Solve For Column k of X: U|X> = |Y> For |X>

//	               [	   --- 	             ]            [           ---                ]
//	          1    |           \ 		     |	    1     |           \                  |
//   <i|X|k> = ------- | <i|Y|k> - /  <i|U|j><j|X|k> | = -------- | <i|Y|k> - /  <i|LU|j><j|Y|k> |
//	       <i|U|i> |           --- 		     |   <i|LU|i> |           ---                |
//	               [           j>i 		     ]	          [           j<i                ]

        luii = luend -1;			// Set <i|LU|i> for last i value, <dlu-1|LU|dlu-1>
        lurow =	1;				// Amount to move from <i|LU|i> -> <i-1|LU|i-1>
        xik = xend - (cols_-k);			// Set <i|X|k> for last i value, <dlu-1|X|k>
        for(; xik>=x00;				// Effective loop over i (LU or U rows)
               lurow++,luii-=lurow,xik-=cols_)	// Work backwards because U is upper-triangular
          {					// Element  <dlu-1||X|k> is solved first
          xjk = xik + cols_;			// Begin with <j|X|k> = <i+1|X|k> (i.e. j=i+1) 
          luij = luii+1;			// Begin with <i|LU|j>=<i|LU|i+1> (i.e. j=i+1)
          for(; xjk<xend; luij++,xjk+=cols_) 	// Effective loop over j (LU cols; X,B,Y rows)
            *xik -= (*luij) * (*xjk);		// sum -= <i|LU|j>*<j|Y|k>, j>i so only LU upper-tri
          *xik /= (*luii);			// This is now <i|X|k>
	  }
        }
      return X;					// Return the results X
      }
      break;
    case n_matrix_type:				// Normal LU Matrix: LU*X = N*X = B
      {
      complex *lui0, *luii, *luij;
      complex *lu00 = ((n_matrix*)LU)->data; 	// Start of data in LU matrix
      complex* luend = lu00 + dlu*dlu;		// End of data in LU: <dlu|LU|dlu+1>
      complex* x00 = X->data;			// Start of data in X: <0|X|k>
      complex* xend = x00 + rows_*cols_;	// End of data in x: <rows_|LU|cols_+1>
      complex *xjk = x00;			// Element of X: <j|X|k> -> <0|X|0>
      complex *xik;				// Element of X: <i|X|k>
      int lurow = dlu+1;			// Amount to move <i|LU|i> to <i+1|LU|i+1> 
      complex sum = 0;
      complex *x0k = x00;
      int ii;
      for(int k=0; k<cols_; k++,sum=0,x0k++)	// Loop over k (X,B columns)
        {					// and solve for column X|k> row-wise.

/*	    First Solve For Column k of Y: L|Y> = |B> For |Y>

  		         --- 	                       ---
  		         \ 			       \
     <i|Y|k> = <i|B|k> - /  <i|L|j><j|Y|k> = <i|B|k> - /  <i|LU|j><j|Y|k> 
  		         --- 		               ---
  		         j<i 			       j<i                   */

        xik = x0k;				// Set <i|X|k> to <0|X|k>, i.e. i=0
        i = 0;
        for(ii=0,lui0=lu00; i<dlu;	 	// Works forwards, <0|Y|k> is solved first
                    i++,xik+=cols_,lui0+=cols_)	// because L is lower triangular
          {
          if(indx[0] > -1)			// Must swap if there has been row permutations
            {					// (negative indx[0] indicates no permutations)
            ip = indx[i];			// Row permutations, row i is swapped with row ip
            sum = X->data[ip*cols_+k];		// Start with actual <i|B|k>, <ip|B|k> for sum
            X->data[ip*cols_+k] = *xik;		// Now swap row for next iteration: <ip|B|k> <- <i|B|k>
            }
          else
            sum = X->data[i*cols_+k];		// Start with actual <i|B|k>, <ip|B|k> for sum
          if(ii)
            for(xjk=x0k,luij=lui0; xjk<xik;	// Effective loop over j (LU columns), j<i
                             xjk+=cols_,luij++) //   (lower triangle elements of LU)
              sum -= (*luij) * (*xjk);		// sum -= <i|LU|j>*<j|B|k>, i>j
          else if(norm(sum))
            ii=1;
          *xik = sum;				// Now set <i|Y|k> to sum (stored in X)
          }
    
//	     Now Solve For Column k of X: U|X> = |Y> For |X>

//	               [	   --- 	             ]            [           ---                ]
//	          1    |           \ 		     |	    1     |           \                  |
//   <i|X|k> = ------- | <i|Y|k> - /  <i|U|j><j|X|k> | = -------- | <i|Y|k> - /  <i|LU|j><j|Y|k> |
//	       <i|U|i> |           --- 		     |   <i|LU|i> |           ---                |
//	               [           j>i 		     ]	          [           j<i                ]

        luii = luend -1;			// Set <i|LU|i> for last i value, <dlu-1|LU|dlu-1>
        xik = xend - (cols_-k);			// Set <i|X|k> for last i value, <dlu-1|X|k>
        for(; xik>=x00; luii-=lurow,xik-=cols_)	// Effective loop over i (LU rows). Work backwards,
          {					// as U is upper-tri; <dlu-1||X|k> is solved first
          xjk = xik + cols_;			// Begin with <j|X|k> = <i+1|X|k> (i.e. j=i+1) 
          luij = luii+1;			// Begin with <i|LU|j>=<i|LU|i+1> (i.e. j=i+1)
          for(; xjk<xend; luij++,xjk+=cols_) 	// Effective loop over j (LU cols; X,B,Y rows)
            *xik -= (*luij) * (*xjk);		// sum -= <i|LU|j>*<j|Y|k>, j>i so only LU upper-tri
          *xik /= (*luii);			// This is now <i|X|k>
	  }
        }
      }
      return X;					// Return the results X
      break;
    default:					// LU is of an Unknown Matrix Type
      {
      complex* x00 = X->data;			// Start of data in X: <0|X|k>
      complex* xend = x00 + rows_*cols_;	// End of data in x: <rows_|LU|cols_+1>
      complex *xjk = x00;			// Element of X: <j|X|k> -> <0|X|0>
      complex *xik;				// Element of X: <i|X|k>
      complex sum = 0;
      int j, ii;
      complex *x0k = x00;
      for(int k=0; k<cols_; k++,sum=0,x0k++)	// Loop over k (X,B columns)
        {					// and solve for column X|k> row-wise.

/*	    First Solve For Column k of Y: L|Y> = |B> For |Y>

  		         --- 	                       ---
  		         \ 			       \
     <i|Y|k> = <i|B|k> - /  <i|L|j><j|Y|k> = <i|B|k> - /  <i|LU|j><j|Y|k> 
  		         --- 		               ---
  		         j<i 			       j<i                   */

        xik = x0k;				// Set <i|X|k> to <0|X|k>, i.e. i=0
        for(ii=0,i=0; i<dlu; i++,xik+=cols_) 	// Works forwards, <0|Y|k> is solved first
          { 					// because L is lower triangular
          if(rswap)				// Must swap if there has been row permutations
            {					// (negative indx[0] indicates no permutations)
            ip = indx[i];			// Row permutations, row i is swapped with row ip
            sum = X->data[ip*cols_+k];		// Start with actual <i|B|k>, <ip|B|k> for sum
            X->data[ip*cols_+k] = *xik;		// Now swap row for next iteration: <ip|B|k> <- <i|B|k>
            }
          else
            sum = X->data[i*cols_+k];		// Start with actual <i|B|k>, <ip|B|k> for sum
          if(ii)
            for(j=0,xjk=x0k;j<i;j++,xjk+=cols_) // loop over j (LU columns), j<i so LU low-triangle
              sum -= (*LU)(i,j) * (*xjk);	// sum -= <i|LU|j>*<j|B|k>, i>j
          else if(norm(sum))
            ii=1;
          *xik = sum;				// Now set <i|Y|k> to sum (stored in X)
          }
    
/*	     Now Solve For Column k of X: U|X> = |Y> For |X>

                  [	      --- 	        ]            [           ---                ]
             1    |           \ 		|       1    |           \                  |
<i|X|k> = ------- | <i|Y|k> - /  <i|U|j><j|X|k> | = -------- | <i|Y|k> - /  <i|LU|j><j|Y|k> |
          <i|U|i> |           --- 		|   <i|LU|i> |           ---                |
                  [           j>i 		]	     [           j<i                ] */

        xik = xend - (cols_-k);			// Set <i|X|k> for last i value, <dlu-1|X|k>
        for(i=dlu-1; xik>=x00; i--,xik-=cols_)	// Loop over i (LU rows). Work backwards,
          {					// as U is upper-tri; <dlu-1||X|k> is solved first
          xjk = xik + cols_;			// Begin with <j|X|k> = <i+1|X|k> (i.e. j=i+1) 
          for(j=i+1; j<dlu; j++,xjk+=cols_) 	// Loop over j (LU cols; X,B,Y rows)
            *xik -= (*LU)(i,j) * (*xjk);	// sum -= <i|LU|j>*<j|Y|k>, j>i so only LU upper-tri
	  *xik /= (*LU)(i,i); 			// This is now <i|X|k>
	  }
        }
      }
      return X;					// Return the results X
      break;
    }
  return X;					// Return the results X
  }


// ____________________________________________________________________________
// AA                 CLASS N_MATRIX DIAGONALIZATION FUNCTIONS
// ____________________________________________________________________________

/* These routines have the ability to diagonalize a general complex matrix. 

      Function
      ========    =============================================================
        cred      Reduces a complex Hermitian matrix to complex	Hermitian tri-
                  diagonal form using the Housholder algorithm.	
        rred      Converts a general complex Hermitian tri-diagonal array into
                  a real symmetric tri-diagonal array.
        tqli      Uses the QL algorithm to reduce an input tri-diagonal 
                  Hermitian array to a symmetric tridiagonal matrix.	
*/

	// Input	nmx	: A normal matrix (complex Hermitian!)
	//		U	: A second normal matrix
	// Output	void    : Both matrices, nmx & U are altered
	//			  Now nmx will be tri-diagonal and U
	//			  will be the transfomation array
	// Note			: It is assumed that array U is zeroed!

void n_matrix::cred(n_matrix& U)
  {
  int i,j,l,dim=rows();
  double s;			// norm of vector a(l+1,l) to a(dim-1,l)
  double sw;			// norm of vector a(l+2,l) to a(dim-1,l)
  double t;
  double csf;			// Column scaling factor
  complex* w;
  w = new complex[dim];

//                 Begin With Transformation Array As Identity

  for(i=0;i<dim;i++) U(i,i)=1.0;
 
//                      Householder Reduction Section
 
/* We Loop Through Columns Of The Input Hermitian Matrix. First, All Elements
   That Are Not Tridiagonal (in the lth Column) Are Checked To See If They Are
   Non-Zero. If So, Then We Do Some Math For The Transformation. Note That On
   Successive Columns We Use The Updated Array With Alterations From The 
   Previous Column Adjustment In Place.                                      */
   
  for(l=0; l<dim-2; l++) 		// Loop rows of input array
    {
 
//      Check For Non-Zero Tridiagonal Elements Lower Part of Column l

    csf = 0.0;				// This will be column scaling factor
    for(i=l+2; i<dim; i++) 		// Get scaling factor, it checks col
      csf += AbsNorm((*this)(i,l));	// l for !0 elements off tri-diagonal

/*            Next Generate The Vector W For Column l If Needed
 
  This section assumes that  1. nll=<l|mx|l>   2. nend=<rows_|mx|rows_>+1
                             3. wl = <l|W>                                   */
 
    if(csf > 0.0)                       // Do nothing if row "tridiagonal"
      {                                 // else we must transform.
      csf += AbsNorm((*this)(l+1,l));	// Adjust scaling factor (off-diag)
      sw = 0.0;
      for(i=l+2; i<dim; i++) 
	sw += square_norm((*this)(i,l)/csf);
      s = sw + square_norm((*this)(l+1,l)/csf);
      s = csf*sqrt(s);
      if((*this)(l+1,l)!=0)			// <l+1|W> =
	{					//                    1+s
	t = norm((*this)(l+1,l));		//    <l+1|mx|l> = ------------
	mul(w[l+1],(*this)(l+1,l),(1.0+s/t));   //                 |<l+1|mx|l>|
	}
      else w[l+1] = -s;				// <l+1|W> = -s
      s = sw+square_norm(w[l+1]/csf);
      s = sqrt(2/s)/csf;
	  
      for(i=l+2; i<dim; i++) 			// Get scaled column l
	mul(w[i],(*this)(i,l),s);		//  w[i] = s*<i|nmx|l>
      w[l+1] = w[l+1]*s;

//                        Transform Matrix From The Left
//		         This Will Alter The Input Array

/* For the current column (l), the diagonal values should already be set.
   The element directly below (tri-diagonal) is determined, as are others
   in below that (all zero).  In the transform, the elements in the lower
   right block are altered to an intermediate state.  This is depicted
   below for the transformation on the 1st and 3rd columns of a 5x5 array
   of <j|hmx|i> (with i>=l, j>l) Are Adjusted.  Note That The Vector w Is
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
 
               for(i=l; i<dim; i++)
                 {
                 ss = 0.0;
	         for(j=l+1; j<dim; j++) ss += conj(w[j],(*this)(j,i));
	         for(j=l+1; j<dim; j++) (*this)(j,i) -= ss*w[j];
                 }                                                           */

	for(i=l; i<dim; i++) 
	  {
	  register double ssi=0, ssr=0;
	  for(j=l+1; j<dim; j++)		// This loop calculates values
	    {					// ssi & ssr for this i
	    complex& t = (*this)(j,i);		// Set t as <j|mx|i>
	    register double wr = Re(w[j]);
	    register double wi = Im(w[j]);
	    register double tr = Re(t);
	    register double ti = Im(t);
	    ssr += wr*tr+wi*ti;
	    ssi += wr*ti-wi*tr;
	    }
	  for(j=l+1; j<dim; j++) 		// Loop adjusts the matrix elements
	    {					// below row l, column i & right
	    complex& t = (*this)(j,i);		// Set t as <j|mx|i>
	    register double wr = Re(w[j]);
	    register double wi = Im(w[j]);
	    Re(t, Re(t) -ssr*wr+ssi*wi);	// This will subtract off
	    Im(t, Im(t) -ssr*wi-ssi*wr);	// parts from element <j|mx|i>
	    }
	  }

//                        Transform Matrix From The Right
//			  This Will Alter The Input Array
 
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
 
               for(i=l; i<dim; i++)
                 {
                 ss = 0.0;
	         for(j=l+1; j<dim; j++) ss += w[j]*(*this)(i,j);
	         for(j=l+1; j<dim; j++) (*this)(i,j) -= conj(w[j],ss);
                 }                                                           */


	  for(i=l; i<dim; i++) 
	    {
	    register double ssi=0, ssr=0;
	    for(j=l+1; j<dim; j++) 		// This loop calculates the
	      {					// values ssi and ssr, this i
	      complex& t = (*this)(i,j);	// Set t = <i|mx|j> = mij
	      register double wr = Re(w[j]);
 	      register double wi = Im(w[j]);
	      register double tr = Re(t);
	      register double ti = Im(t);	// Now just add wj*mij to ss
	      ssr += wr*tr-wi*ti;		// ssr+=Re(wj)Re(mij)-Im(wj)Im(mij)
	      ssi += wr*ti+wi*tr;		// ssi+=Re(wj)Im(mij)+Im(wj)Re(mij)
	      }
	    for(j=l+1; j<dim; j++) 		// This loop adjusts
	      {					
	      complex& t = (*this)(i,j);	// Set t = <i|mx|j>
	      register double wr = Re(w[j]);
	      register double wi = Im(w[j]);
	      Re(t,Re(t)-ssr*wr-ssi*wi);	// Alters <i|mx|j> by
	      Im(t,Im(t)+ssr*wi-ssi*wr);	// subtractiong conj(wj)*ss
	      }
	    }
             
//                        Accumulate Transformations
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

                    for(i=0; i<dim; i++) 
	               {
	               ss = 0.0;
	               for(j=l+1; j<dim; j++) ss += w[j]*U(i,j);
	               for(j=l+1; j<dim; j++) U(i,j) -= conj(w[j],ss);
	               }                                                     */
 
       for(i=0; i<dim; i++) 			// Loop all rows of U
	{					// We first need to determine
	register double ssi=0, ssr=0;		// ssi & ssr for this row from
	for(j=l+1; j<dim; j++)			// upper triangle row elements
	  { 					// of U.
	  complex& t = U(i,j);			//   Get t = <i|U|j>
	  register double wr = Re(w[j]);	//   Re<j|w>
	  register double wi = Im(w[j]);	//   Im<j|w>
	  register double tr = Re(t);		//   Re<i|U|j>
	  register double ti = Im(t);		//   Im<i|U|j>
	  ssr += wr*tr-wi*ti;			//   Adjust ssr column j
	  ssi += wr*ti+wi*tr;			//   Adjust ssi column j
	  }
	for(j=l+1; j<dim; j++) 			// Having ssr & ssj for row i,
	  {					// adjust U row i upper triang. 
	  complex& t = U(i,j);			//   Get t = <i|U|j>
	  register double wr = Re(w[j]);	//   Re<j|w>
	  register double wi = Im(w[j]);	//   Im<j|w>
	  Re(t, Re(t)-ssr*wr-ssi*wi);		//   Adjust Re(<i|U|j>)
	  Im(t, Im(t)+ssr*wi-ssi*wr);		//   Adjust Im(<i|U|j>)
	  }
	}
      }
    }
  delete [] w;
  }



/*************************************************************************
**                                                                      **
** The routine rred will convert a general complex Hermitian tri-       **
** diagonal array into a real symmetric tri-diagonal array.  This       **
** result may be subsequently used for a complete diagonalization.      **
** The tranformation is diagonal and unitary.                           **
**                                                                      **
**                                                                      **
**	Input		nmx	: Input normal array (this)		**	
**				  Assumed Hermitian & tri-diagonal	**
**			mxev	: Any previously accumulated 		**
**				  transformations in a normal array.	**
**				  If no previous transformations this	**
**				  should be an identity matrix (normal)	**
** 	Output		nmx	: Real symmetric tri-diagonal array	**
**			mxev	: All accumulated transformations.	**	 
**                                                                      **
**  The routine loops over each of the diagonal elements of the tri-	**
**  diagonal array starting with the 2nd column (index 1).  It then	**
**  looks at the element directly left which may be complex.  If found **
**  this is then made real (as is its symmetric partner), the next two	**
**  off-diagonals set appropriatly, & the transformation matrix fixed.	**
**  Below is a depiction of what the input tri-diagonal array is doing  **
**  during the looping over a 5x5 array.				**
**                                                                      **
**  x x . . .    x x . . .    x r . . .	   x r . . .    x r . . .	**
**  x x x . .    ? X x . .    r x A . .	   r x A . .    r x r . .	**
**  . x x x . -> . x x x . -> . B x x .	-> . ? X x . -> . r x C . ->	**
**  . . x x x    . . x x x    . . x x x	   . . x x x    . . D x x	**
**  . . . x x    . . . x x    . . . x x	   . . . x x    . . . x x	**
**                                                                      **
**   Initial    Check Row 1     Adjust    Check Row 2     Adjust	**
**                                                                      **
**  x r . . .    x r . . .    x r . . .	   x r . . .   			**
**  r x r . .    r X r . .    r x r . .	   r x r . .   			**
**  . r x C . -> . r x r . -> . r x r .	-> . r X r . 			**
**  . . ? X x    . . r x E    . . r x E	   . . r x r   			**
**  . . . x x    . . . F x    . . . ? X	   . . . r x   			**
**                                                                      **
** Check Row 3     Adjust     Check Last     Adjust    			**
**                                                                      **
**  Above the x's are just the original hermitian tri-diagonal elements **
**  and the big X's are the diagonals used in the loop.  I used r's to	**
**  depict altered elements which are real (the norm of the element at  **
**  the preceding step) and A,B,C,D to depict those modified but not	**
**  yet to their final form.  Note that these appear in Hermitian pairs	**
**  { A = conj(B), C = conj(D), etc.} and that it is the elements on 	**
**  the lower diagonal (B, D, F..) that are initially modified from the	**
**  originals and then subsequently used for the next transformation.	**
**  It doesn't appear the uppers, the conjugates of the lowers, are 	**
**  used at all (A, C, E,..) even though they are set? Perhaps this is	**
**  to compensate for an input array not fully Hermitian?		**
**                                                                      **
*************************************************************************/

/* sosi - There are a few ways this can be sped up I think.
    1.) Replace the line nmx(l+1,l) = conj(ss,(*this)(l+1,l));
        with the line    nmx(l+1,l) = conj(ss, nmx(l+1,l))
    2.) Forget about even setting nmx(l,l+1)*=ss since we don't seem
        to ever use it and we definetely overwrite it if it is zero 
    3.) Switch to direct element access rather than recomputing the
        element index here.
*/
  
void n_matrix::rred(n_matrix& mxev, int newU)
  {
  n_matrix &nmx = *this;			// A handle on the array
  int i, l, nr=rows();				// These are just a set
  double normod;				// of working variables
  complex sf;
  for(l=1; l<nr; l++)				// Loop nmx rows
    {
    normod = norm(nmx(l,l-1));			// Get norm of <l|nmx|l-1> 
    if(normod!=0)				// (direcly left of <l|nmx|l>)
      {						// Act only if its non-zero
      sf = conj(nmx(l,l-1))/normod;		// <l|nmx|l-1>/||<l|nmx|l-1>||
      nmx(l,l-1) = normod;			// Set <l|nmx|l-1>=normod (left)
      nmx(l-1,l) = normod;			// Set <l-1|nmx|l>=normod (above)
      if(l < nr-1)				// Apply transform if 
        {					// still two rows from end
        nmx(l,l+1) *= sf;			// Mult. <l|nmxl+1> by sf (right)
        nmx(l+1,l) = conj(sf,(*this)(l+1,l));   // <l+1|nmx|l> = conj(s)*<l+1|this|l>
        }					// (below, using original below?)
      for(i=0; i<nr; i++) 			// Accumulate transform in
        mxev(i,l) = conj(sf,mxev(i,l));		// the mxev array now
      }
    }
  }



/********************************************************
*							*
*	The function sign is the C++ implementation	*
*	of the standard FORTRAN library function.	*
*							*
********************************************************/
inline double sign(double a,double b)
  { return double((b<0.0)?-fabs(a):fabs(a)); }


/*************************************************************************
**                                                                      **
** The routine tqli uses the QL algorithm to reduce an input tri-	**
** diagonal Hermitian array to a symmetric tridiagonal matrix.		**
**                                                                      **
**	Input:	nmx	: A symmetric tridiagonal matrix (this)		**
**		z	: Array of previous accumulated transformations	**
*	Output: a	: The tri-diagonalized matrix result		**
*		z	: Assumulated transformations.			**
**                                                                      **
*************************************************************************/

void n_matrix::tqli(n_matrix& z, d_matrix& a)
  {
  int i,k,l,m,iter,n=rows();
  double dd;
  double s,r,p,g,f,c,b;
  double x;
  double *e;
  e = new double[n];
  complex xx;
  if(n>0)
    {
    for(i=0;i<n-1;i++) {	// copy to workspace
      e[i]=zRe((*this)(i,i+1));
    };
    convert(&a);
    e[n-1] = 0.0;
    iter=0;
    for(l=0;l<n;l++)
      {						// Loop for eigenvalues
      do {
        m=l;
        while (1) {		// look for an eigenvalue
          if (m==n-1) break;	// all eigenvalues found
          dd=fabs(zRe(a(m,m)))+fabs(zRe(a(m+1,m+1)));
          if (fabs(e[m])+dd==dd) break;	// a new e.value found
          m=m+1;
        };
        if (m!=l) {
          if (iter>10*n)
            {
            NMxerror(29,1);			// Too many iterations
            NMxfatal(28);			// Unable to diagonalize
            }					// Almost surely bad input
          iter=iter+1;
          // form shift
          g=(Re(a(l+1,l+1))-Re(a(l,l)))/(2.0*e[l]);
          r=sqrt(g*g+1.0);
          x=e[l];
          x=x/(g+sign(r,g));
          g=Re(a(m,m))-Re(a(l,l))+x;
          s=1.0;
          c=1.0;
          p=0.0;
          // QL step
          for(i=m-1; i>=l; i--) {
            f=s*e[i];
            b=c*e[i];
            if (fabs(f)>=fabs(g)) {
              c=g/f;
              r=sqrt(c*c+1.0);
              e[i+1]=f*r;
              s=1.0/r;
              c=c*s;
            } else {
              s=f/g;
              r=sqrt(s*s+1.0);
              e[i+1]=g*r;
              c=1.0/r;
              s=s*c;
            };
            g=Re(a(i+1,i+1))-p;
            r=(Re(a(i,i))-g)*s+2.0*c*b;
            p=s*r;
            Re(a(i+1,i+1),g+p);
            g=c*r-b;
            // accumulate transforms
            for (k=0;k<n;k++) {
              complex& zi =z(k,i);
              complex& zi1=z(k,i+1);
	      mul (xx, s, zi);
              sub( zi, c*zi, s*zi1);
              ::add( zi1, xx, c*zi1);
            };
          };
          Re(a(l,l),Re(a(l,l))-p);
          e[l]=g;
          e[m]=0.0;
        };
      } while (m!=l);		// all eigenvalues found
    };
    }
  delete [] e;
  }


/*
     this subroutine is a translation of a complex analogue of
     the algol procedure orthes, num. math. 12, 349-368(1968)
     by martin and wilkinson.
     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).

     given a complex general matrix, this subroutine
     reduces a submatrix situated in rows and columns
     low through igh to upper hessenberg form by
     unitary similarity transformations.

     on input
        low and igh are integers determined by the balancing
          subroutine  cbal.  if  cbal  has not been used,
          set low=1, igh=n.

        a contains the complex input matrix.

     on output

        a of the hessenberg matrix.  information
          about the unitary transformations used in the reduction
          is stored in the remaining triangles under the
          hessenberg matrix.

        ort contains further information about the
          transformations.  only elements low through igh are used.

*/


complex* n_matrix::corth  (int low, int igh)

{
  int i,j,m;
  double f,g,h,sc;
  complex z;
  n_matrix& a=*this;
  int n=a.rows();
  complex *ort;
  ort=new complex[n];

  for (m=low+1; m<igh; m++) 
    {
      h=0;
      ort[m]=0;
      sc=0;
      for (i=m; i<=igh; i++) 
	sc+=AbsNorm( a(i,m-1) );

      if (sc>0.0) 
	{
	  for (i=igh;i>=m;i--) 
	    {
	      ort[i]= a(i,m-1) /sc;
	      h+=square_norm(ort[i]);
	    };
	  g=sqrt(h);
	  f=norm(ort[m]);
	  if (f==0.0)
	    {
	      ort[m]=g;
	      Re(a(m,m-1),sc);
	    }
	  else 
	    {
	      h+=f*g;
	      g/=f;
	      ort[m]*=(1.0+g);
	    };

	  for (j=m;j<n;j++)
	    {
	      z=0.0;
	      for (i=igh;i>=m;--i) 
		z += conj(ort[i],a(i,j));
	      z/=h;
	      for (i=m;i<=igh;i++) 
		a(i,j)-=z*ort[i];
	    };

	  for (i=0;i<=igh;i++) 
	    {
	      z=0.0;
	      for (j=igh;j>=m;--j) 
		z += ort[j]*a(i,j);
	      z/=h;
	      for (j=m;j<=igh;j++) 
		a(i,j)-= conj(ort[j],z);
	    };

	  ort[m]*=sc;
	  a(m,m-1) *= -g;
	};
    };
  return ort;
};

/*
     this subroutine is a translation of a unitary analogue of the
     algol procedure  comlr2, num. math. 16, 181-204(1970) by peters
     and wilkinson.
     handbook for auto. comp., vol.ii-linear algebra, 372-395(1971).
     the unitary analogue substitutes the qr algorithm of francis
     (comp. jour. 4, 332-345(1962)) for the lr algorithm.

     this subroutine finds the eigenvalues and eigenvectors
     of a complex upper hessenberg matrix by the qr
     method.  the eigenvectors of a complex general matrix
     can also be found if  corth  has been used to reduce
     this general matrix to hessenberg form.

     on input

        low and igh are integers determined by the balancing
          subroutine  cbal.  if  cbal  has not been used,
          set low=1, igh=n.

        ort contains information about the unitary trans-
          formations used in the reduction by  corth, if performed.
          only elements low through igh are used.  if the eigenvectors
          of the hessenberg matrix are desired, set ortr(j) and
          orti(j) to 0.0d0 for these elements.

        h contains the complex upper hessenberg matrix.
          their lower triangles below the subdiagonal contain further
          information about the transformations which were used in the
          reduction by  corth, if performed.  if the eigenvectors of
          the hessenberg matrix are desired, these elements may be
          arbitrary.

     on output

        ort, and the upper hessenberg portions of h
          have been destroyed.

        w contains the eigenvalues.  if an error
          exit is made, the eigenvalues should be correct
          for indices ierr+1,...,n.

        z comtains the eigenvectors.  the eigenvectors
          are unnormalized.  if an error exit is made, none of
          the eigenvectors has been found.

        ierr is set to
          zero       for normal return,
          j          if the limit of 30*n iterations is exhausted
                     while the j-th eigenvalue is being sought.

*/


int n_matrix::comqr3(int low,int igh,
		     complex *ort, d_matrix& w, n_matrix& z,
		     int flag)
{
  n_matrix &a=*this;
  int n=a.rows();
  int i,j,k,l=0,m,en,ii,jj,ll,nn,ip1,its,enm1,iend,lp1;
  double u,s,ss,nor,machep;
  complex sc,x,y,t,zz;
#define EPS	1.0e-13
  machep=EPS;
  /* unit matrix */
  for (i=0;i<n;i++) 
    for (j=0;j<n;j++) 
      z(i,j)=0;

  for (i=0;i<n;i++) 
    z(i,i)=1.0;

  /* accumulate transformations */
  iend=igh-low-1;
  if (iend>=0) 
    {
      for (ii=1;ii<=iend;ii++) 
	{
	  i=igh-ii;
	  if (ort[i-1]!=0.0) 
	    {
	      if (a(i-1,i-2)!=0.0) 
		{
		  nor=Re(conj(a(i-1,i-2),ort[i-1]));
		  ip1=i+1;
		  for (k=ip1;k<=igh;k++) 
		    ort[k-1]=a(k-1,i-2);
		  for (j=i;j<=igh;j++) 
		    {
		      sc=0.0;
		      for (k=i;k<=igh;k++) 
			sc += conj(ort[k-1],z(k-1,j-1));
		      sc/=nor;
		      for (k=i;k<=igh;k++) 
			z(k-1,j-1) +=sc*ort[k-1];
		    };
		};
	    };
	};
      /* real subdiagonal elements */
      l=low+1;
      for (i=l;i<=igh;i++) 
	{
	  if (i+1<igh)
	    ll=i+1;
	  else 
	    ll=igh;
	  if (Im(a(i-1,i-2))!=0.0) 
	    {
	      nor=norm(a(i-1,i-2));
	      y=a(i-1,i-2)/nor;
	      a(i-1,i-2)=nor;
	      for (j=i;j<=n;j++) 
		a(i-1,j-1)=conj(y,a(i-1,j-1));
	      for (j=1;j<=ll;j++) 
		a(j-1,i-1)*=y;
	      for (j=low;j<=igh;j++) 
	        z(j-1,i-1)*=y;
	    };
	};
    };

  /* isolated roots */
  for (i=1;i<=n;i++)
    if ((i<low)||(i>igh))
      w(i-1,i-1)=a(i-1,i-1);

  en=igh;
  t=0.0;

  /* iterate for eigenvalues */
  while (en>=low) 
    {
      its=0;
      enm1=en-1;
      while (1) 
	{
	  for (ll=low;ll<=en;ll++) 
	    {
	      l=en+low-ll;
	      if (l==low) break;
	      s=AbsNorm(a(l-2,l-2))+AbsNorm(a(l-1,l-1));
	      if (s<1.0) s=1.0;
	      if (fabs(Re(a(l-1,l-2)))<=machep*s) break;
	    };
	  if (l==en) break;
	  if (its==30)
	    return en;
	  /* form shift */
	  if ((its==10)||(its==20)) 
	    sc=fabs(Re(a(en-1,enm1-1))) 
	      + fabs(Re(a(enm1-1,en-3)));
	  else 
	    {
	      sc=a(en-1,en-1);
	      x=a(enm1-1,en-1)*Re(a(en-1,enm1-1));
	      if (x!=0.0)
		{
		  y=(a(enm1-1,enm1-1)-sc)/2.0;
		  zz = sqrt(y*y+x);
		  if (Re(conj(y,zz))<0.0) 
		    zz = -zz;
		  zz = x/(y+zz);
		  sc-=zz;
		};
	    };
	  for (i=low;i<=en;i++) 
	    a(i-1,i-1)-=sc;
	  t+=sc;
	  its=its+1;
	  /* QR decomposition */
	  lp1=l+1;
	  for (i=lp1;i<=en;i++)
	    {
	      s=Re(a(i-1,i-2));
	      nor=sqrt( square_norm(a(i-2,i-2))+s*s );
	      x=a(i-2,i-2)/nor;
	      w(i-2,i-2)=x;
	      a(i-2,i-2)=nor;
	      s /= nor;
	      a(i-1,i-2)=complex(0,s);
	      for (j=i;j<=n;j++) 
		{
		  y=a(i-2,j-1);
		  zz=a(i-1,j-1);
		  a(i-2,j-1)= conj(x,y)+s*zz;
		  a(i-1,j-1)= x*zz - s*y;
		};
	    };
	  s=Im(a(en-1,en-1));
	  if (s!=0.0)
	    {
	      nor=norm(a(en-1,en-1));
	      sc=conj(a(en-1,en-1)/nor);
	      a(en-1,en-1)=nor;
	      if (en!=n) 
		for (j=en+1;j<=n;j++)
		  a(en-1,j-1)*=sc;
	    };
	  /* calculate RQ */
	  for (j=lp1;j<=en;j++) 
	    {
	      x=w(j-2,j-2);
	      ss=Im(a(j-1,j-2));
	      for (i=1;i<=j;i++) 
		{
		  y=Re(a(i-1,j-2));
		  zz=a((i-1),j-1);
		  if (i!=j) 
		    {
		      Im(y,Im(a((i-1),j-2)));
		      Im(a((i-1),j-2), 
			 Re(x)*Im(y)+Im(x)*Re(y)+ss*Im(zz));
		    };
		  Re(a((i-1),j-2), 
		     Re(x)*Re(y)-Im(x)*Im(y)+ss*Re(zz));
		  a((i-1),j-1)=conj(x,zz)-ss*y;
		};
	      for (i=low;i<=igh;i++) 
		{
		  y=z((i-1),j-2);
		  zz=z((i-1),j-1);
		  z((i-1),j-2)=x*y+ss*zz;
		  z((i-1),j-1)=conj(x,zz)-ss*y;
		};
	    };
	  sc=conj(sc);
	  if (s!=0.0) 
	    {
	      for (i=1;i<=en;i++) 
		a((i-1),en-1)*=sc;
	      for (i=low;i<=igh;i++) 
		z((i-1),en-1)*=sc;
	    };
	};
      /* a root found */
      a((en-1),en-1)+=t;
      w(en-1,en-1)=a((en-1),en-1);
      en=enm1;

    };


  if (flag) 
    {
      /* all roots found, backsubstitution of eigenvectors */
      nor=0.0;
      for (i=1;i<=n;i++) 
	for (j=i;j<=n;j++) 
	  nor+=AbsNorm(a((i-1),j-1));
      if ((n==1)||(nor==0.0)) return 0;
      for (nn=2;nn<=n;nn++) 
	{
	  en=n+2-nn;
	  x=w(en-1,en-1);
	  a(en-1,en-1)=1;
	  enm1=en-1;
	  for (ii=1;ii<=enm1;ii++) 
	    {
	      i=en-ii;
	      zz=a((i-1),en-1);
	      if (i!=enm1)
		{
		  ip1=i+1;
		  for (j=ip1;j<=enm1;j++) 
		    zz+= a((i-1),j-1)*a((j-1),en-1);
		};
	      y=x-w(i-1,i-1);
	      if ((fabs(Re(y))<EPS)&&(fabs(Im(y))<EPS)) 
		Re(y,machep*nor);
	      a((i-1),en-1)=zz/y;
	    };
	};
    };
  enm1=n-1;
  /* eigenvectors of isolated root */
  for (i=1;i<=enm1;i++) 
    {
      if ((i<low)||(i>igh)) 
	for (j=i+1;j<=n;j++) 
	  z((i-1),j-1)=a((i-1),j-1);
    };
  /* multiply for eigenbase */
  for (jj=low;jj<=enm1;jj++) 
    {
      j=n+low-jj;
      if (j-1<igh) m=j-1;
      else m=igh;
      for (i=low;i<=igh;i++) 
	{
	  zz=z((i-1),j-1);
	  for (k=low;k<=m;k++) 
	    zz+=z((i-1),k-1)*a((k-1),j-1);
	  z((i-1),j-1)=zz;
	};
    };
  if (flag) 
    {
      /* normalize vectors */
      for (i=0;i<n;i++) 
	{
	  u=0.0;
	  for (j=0;j<n;j++) 
	    u=u+AbsNorm(z(j,i));
	  s=0.0;
	  for (j=0;j<n;j++) 
	    s+=square_norm(z(j,i)/u);
	  s=u*sqrt(s);
	  for (j=0;j<n;j++) 
	    z(j,i)/=s;
	};
    };
  return 0;
};




/*************************************************************************
**									**
** The routine diag combines the routines cred,	rred and tqli (using	**
** the function sign) to diagonalize a complex Hermitian matrix.  It	**
** uses the routines corth and comqr3 to diagonalise general complex 	**
** matrices.								**
**									**
** During the course of the transformations, exclusively unitary	**
** (orthogonal) transformations are used. Errors come mainly from	**
** rounding. The resulting eigenbase of a Hermitian array deviates 	**
** from unitarity only negligibly.					**
**									**
** This is THE Entry Point to Normal Complex Array Diagonalizations!	**
** This Is The Workhorse Behind Most GAMMA Calculations............!	**
** This Routine Needs to Be Working Very Well......................!	**
**									**
** Note That There Is Currently One Unique Feature Of This Function	**
** That Sets It Apart From Most Other GAMMA Matrix Functions. It Sets	**
** The Pointers mxd and mxev To Point To A Diagonal And A Normal Array	**
** Respectively.  In Turn, The Function Calls Diagonalization Routines	**
** That EXCLUSIVELY TAKE d_matrix and n_matrix!				**
**									**
*************************************************************************/
// sosi - use of this diag routine is a big pain!  I think the reason is
//        so that the pointers to mxd and mxev keep their types properly
//	  but most other functions have no problem with it anyway......

void n_matrix::diag(_matrix* (&mxd), _matrix* (&mxev))

        // Input        nmx     : A normal complex array, n_matrix (this)
        //              mxd     : Pointer to mx to become diagonal mx
        //              mxev    : Pointer to mx to become eigenvectors mx
        // Output       void    : The matrices mxd and mxev are set to
        //                        the diagonal matrix of nmx eigenvalues and
        //                        the matrix of dmx eigenvectors respectively
        // Note                 : The reference count to nmx must be twice
        //                        incremented external by the call origin.

  {
  int nrows = this->rows();
  int ncols = this->cols();
  mxd  = new d_matrix(nrows,ncols);			// Alocate space for eigenvalues
  mxev = new n_matrix(nrows,ncols,complex0); 		// Alocate space for eigenvectors	
  if(test_hermitian())				// See if matrix is Hermitian
    {
#if defined(_USING_SUNPERFLIB_) || defined(_USING_LAPACK_)
    if(nrows < 16)  // then do it the old fashioned gamma way.
    {
      n_matrix mx2(*this);			// Workspace to diagonalization
      mx2.cred(*(n_matrix *)mxev);		// To Hermitian tridiagonal form
      mx2.rred(*(n_matrix *)mxev);		// To real symm. tridiagonal form
      mx2.tqli(*(n_matrix *)mxev, 		// To diagonal form
                      *(d_matrix *)mxd);	
      ((n_matrix *)mxev)->unitary = true;		// Eigenvector array is unitary
//    std::cerr << "nh_matrix diag: using gamma code\n";
    }
    else
    {
	n_matrix* hmx = (n_matrix *)this->transpose();
#ifdef _USING_SUNPERFLIB_
        char jobz = 'V';  // Calculate "eigenvectors" AND "eigenvalues".
        char uplo = 'U';  // Upper triangular...
        int  N    = nrows;
        int  lda  = nrows;
        double *w_eig = new double[N];
        int info = -55555;
        zheev(jobz, uplo, N, (doublecomplex *)hmx->data, lda, w_eig, &info);
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
        zheev_(&jobz, &uplo, &N, (__CLPK_doublecomplex *) hmx->data, &lda, w_eig, (__CLPK_doublecomplex *) work, &lwork, rwork, &info);
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
        // Copy results back into mxd and mxev.
	complex *hmxp = hmx->data;
        for(int i=0; i<nrows; i++)
        { ((d_matrix *)mxd)->put( w_eig[i], i, i);
          for(int j=0; j<ncols; j++)
          { ((n_matrix *)mxev)->put(hmxp[j*nrows+i], i, j);
          }
        }
        ((n_matrix *)mxev)->unitary = true;		// Eigenvector array is unitary
	delete [] w_eig;
	delete hmx;
#ifdef _USING_LAPACK_
	delete [] work;
	delete [] rwork;
#endif
//	std::cerr << "nh_matrix diag: using sunperflib code\n";
    }
#else
    n_matrix mx2(*this);			// Workspace to diagonalization
    mx2.cred(*(n_matrix *)mxev);		// To Hermitian tridiagonal form
    mx2.rred(*(n_matrix *)mxev);		// To real symm. tridiagonal form
    mx2.tqli(*(n_matrix *)mxev, 		// To diagonal form
                      *(d_matrix *)mxd);	
    ((n_matrix *)mxev)->unitary = true;		// Eigenvector array is unitary
#endif
    }
  else						// Matrix isnt Hermitian
    {
#if defined(_USING_SUNPERFLIB_) || defined(_USING_LAPACK_)
    if(nrows < 128)  // then do it the old fashioned gamma way.
    {
      n_matrix mx2(*this);			// Workspace to diagonalization
      complex *tmp = mx2.corth(0,nrows-1);		// To upper Hessenberg form
      int err = mx2.comqr3(1, nrows, tmp,		// To diagonal form
          *(d_matrix*)mxd,*(n_matrix*)mxev, 1);
      if(err != 0)
      { NMxerror(29,1);				// Too many iterations
        NMxfatal(28);				// Unable to diagonalize
      }						// Almost surely bad input
      delete [] tmp;				// Delete tmp array from corth
//    std::cerr << "n_matrix diag: using gamma code\n";
    }
    else
    {
	n_matrix* hmx = (n_matrix *)this->transpose();
	n_matrix* hmxr = new n_matrix(nrows, ncols);
#ifdef _USING_SUNPERFLIB_
        char jobvl = 'N';  // Calculate "eigenvectors" AND "eigenvalues".
        char jobvr = 'V';  // Upper triangular...
        int  N    = nrows;
        int  lda  = nrows;
        complex *w_eig = new complex[N];
        int info = -55555;
        zgeev(jobvl, jobvr, N, (doublecomplex *)hmx->data, lda, (doublecomplex *)w_eig, 
	      (doublecomplex *)hmxr->data,N,(doublecomplex *)hmxr->data,N,&info);
#endif
#ifdef _USING_LAPACK_
        char jobvl = 'N';  // Calculate "eigenvectors" AND "eigenvalues".
        char jobvr = 'V';  // Upper triangular...
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
        complex *w_eig = new complex[N];
        complex *work= new complex [2*lwork];
	double *rwork = new double [3*N-2];
        zgeev_(&jobvl, &jobvr, &N, (__CLPK_doublecomplex *)hmx->data, &lda, (__CLPK_doublecomplex *)w_eig,
	       (__CLPK_doublecomplex *)hmxr->data,&N,(__CLPK_doublecomplex *)hmxr->data,&N,(__CLPK_doublecomplex *) work, &lwork, rwork, &info);
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
        // Copy results back into mxd and mxev.
	complex *hmxp = hmxr->data;
        for(int i=0; i<nrows; i++)
        { mxd->put( w_eig[i], i, i);
          for(int j=0; j<ncols; j++)
          { mxev->put(hmxp[j*nrows+i], i, j);
          }
        }
	delete hmxr;
	delete hmx;
	delete []  w_eig;
#ifdef _USING_LAPACK_
	delete [] work;
	delete [] rwork;
#endif
//	std::cerr << "n_matrix diag: using sunperflib code\n";
    }
#else
    n_matrix mx2(*this);			// Workspace to diagonalization
    complex *tmp = mx2.corth(0,nrows-1);		// To upper Hessenberg form
    int err = mx2.comqr3(1, nrows, tmp,		// To diagonal form
        *(d_matrix*)mxd,*(n_matrix*)mxev, 1);
    if(err != 0)
       {
       NMxerror(29,1);				// Too many iterations
       NMxfatal(28);				// Unable to diagonalize
       }					// Almost surely bad input
    delete [] tmp;				// Delete tmp array from corth
#endif
    }
  }

/*                         Diagonalization Routines

  Note That The Returned Diagonal & Eigenvectors Arrays MUST Have Their
  Referencing Set External To This Class. Also, As These Enter The Routine  
  As Pointers Which Will Be Set To Point To the Newly Transformed Arrays,
  They Should NOT Point To Any Array (Use No Memory But The Pointer)
 
           Input        nmx     : A general complex matrix (this)
                        mxd     : Pointer to mx to become diagonal mx
                        mxev    : Pointer to mx to become eigenvector mx
           Output       void    : The matrices mxd and mxev are set to
                                  the Hermitian matrix of dmx eigenvalues and
                                  the matrix of dmx eigenvectors respectively
           Note                 : For h_matrix, mxd should be real and
                                  mxev both unitary and Hermitian
           Note                 : The reference count to dmx must be twice
                                  incremented external by the call origin.   */

//#include <sys/time.h>
//#include <sys/resource.h>

std::vector<int> n_matrix::BlockDiag(_matrix* (&BD), std::vector<int> &U)
  {
//std::cout << "This is n_matrix block diag\n";
//struct rusage me;
//getrusage(0, & me);
//std::cout << "block-diag routine: point 01: " << me.ru_utime.tv_sec+me.ru_utime.tv_usec/1.0e6 << " seconds\n";
//std::cout.flush();
  int nr = rows_;				// Matrix dimension
  //int count = 0;				// Count permutations
  BD = new n_matrix(*this);			// Start with copy

  int i,nze=0,bs=0,t5,t6;
  for(i=0; i<nr; i++) { U.push_back(i); }	// Set both as unpermuted
  std::vector<int> blkdims;
  for(t5=0;t5<nr;++t5)
  { if(nze==t5)
      ++nze;
    for(t6=nze;t6<nr;++t6)
    { if((*BD).get(t5,t6) != complex0 || (*BD).get(t6,t5) != complex0)
      { if(t6 != nze)
        { complex z;
          int z1;
	  int t7;
//	  int k1, k2;
//        getrusage(0, & me);
//        std::cout << "block-diag routine: before perm: " << me.ru_utime.tv_sec+me.ru_utime.tv_usec/1.0e6 << " seconds\n";
//	  std::cout << "permuting "<< nze << " and " << t6 << "\n";
//	  for(k1=0;k1<nr;++k1)
//	  { for(k2=0;k2<nr;++k2)
//	       std::cout << (*BD).get(k1,k2) << "  ";
//	    std::cout << "\n";
//	  }
//        std::cout.flush();
	  for(t7=0;t7<nr;++t7)
	  {  z=(*BD).get(t7,t6);
	     (*BD).put((*BD).get(t7,nze),t7,t6);
	     (*BD).put(z,t7,nze);
	  }
	  for(t7=0;t7<nr;++t7)
	  {  z=(*BD).get(t6,t7);
	     (*BD).put((*BD).get(nze,t7),t6,t7);
	     (*BD).put(z,nze,t7);
	  }
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

void n_matrix::HermTriDiag(_matrix* (&HTD), _matrix* (&U))
  {
  if(!test_hermitian())				// See if matrix is Hermitian
    {
    std::cout << "\n\tn_matrix: Cannot Form Hermitian Tri-Diagonal Form, Sorry";
    return;
    }
  n_matrix *mx2 = new n_matrix(*this);		// Workspace to diagonalization
  int n = rows_;				// Matrix dimension
  U = new n_matrix(n,n,complex0); 		// Alocate space for eigenvectors	
  mx2->cred(*(n_matrix *)U);			// To Hermitian tridiagonal form
  HTD = mx2;
  }

void n_matrix::SymTriDiag(_matrix* (&STD), _matrix* (&U))
  {
  if(!test_hermitian())				// See if matrix is Hermitian
    {
    std::cout << "\n\tn_matrix: Cannot Form Symmetric Tri-Diagonal Form, Sorry";
    return;
    }
  n_matrix *mx2 = new n_matrix(*this);		// Workspace to diagonalization
  int dim = rows_;				// Matrix dimension
  U = new n_matrix(dim,dim,complex0); 		// Alocate space for eigenvectors	
  for(int i=0;i<dim;i++) (*U)(i,i)=1.0;		// Begin with U=I here
  mx2->rred(*(n_matrix *)U);			// To Symmetric tridiagonal form
  STD = mx2;
  }

void n_matrix::SymDiag(_matrix* (&D), _matrix* (&U))
  {
  if(!test_hermitian())				// See if matrix is Hermitian
    {
    std::cout << "\n\tn_matrix: Cannot Form Diagonal Form, Sorry";
    return;
    }
//  n_matrix *mx2 = new n_matrix(*this);		// Workspace to diagonalization
  int dim = rows_;				// Matrix dimension
  U = new n_matrix(dim,dim,complex0); 		// Alocate space for eigenvectors	
  for(int i=0;i<dim;i++) (*U)(i,i)=1.0;		// Begin with U=I here
  D = new d_matrix(dim, dim);			// Alocate space for eigenvalues
  tqli(*(n_matrix *)U, *(d_matrix *)D);		// To diagonal form
  }


// ____________________________________________________________________________
// BB                 CLASS N_MATRIX INVERSION FUNCTIONS
// ____________________________________________________________________________

/* LU Decomposition *************************************-*-c++-*-
**							 	**
**  This routine is taken from "Numerical Recipies", W.H.Press,	**
**  B.P.Flannery, S.A.Teukolsky, & W.T.Vetterling, Cambridge	**
**  University Press, 1989.  See pages 31 through 36.  The code **
**  has been modified for GAMMA into C++.			**
**							 	**
**  The algorithm takes an input matrix A & overwrites it with  **
**  the LU decompositions of A' where A' is A with any needed	**
**  permutations. Any permutations used	are stored in the	**
**  integer array indx.  The function returns the integer d	**
**  which will be +1 if an even	number of permutations was used **
**  or -1 if an odd number of permutations was necessary.  Sub- **
**  sequently, d can be used in taking the determinant of the	**
**  input matrix A (or LU).					** 
**								**  
**  The LU decomposition formulates the equation		**
**								**  
**		     	        A = LU				**  
**								**  
**  where L is lower triangular and U is upper triangular.	**
**								**
**  [A11 A12 A13 A14]   [L11  0   0   0 ] [U11 U12 U12 U14]	**
**  |A21 A22 A23 A24]   [L21 L22  0   0 ] [ 0  U22 U23 U24|	**
**  |A31 A32 A33 A34] = [L31 L32 L33  0 ] [ 0   0  U33 U34|	**
**  [A41 A42 A43 A44]   [L41 L42 L43 L44] [ 0   0   0  U44]	**
**								**
**  Both L and U have elements on the diagonal, but those of L	**
**  are always set to be 1: <i|L|i> = 1.			**
**								**
**  [A11 A12 A13 A14]   [ 1   0   0   0 ] [U11 U12 U12 U14]	**
**  |A21 A22 A23 A24]   [L21  1   0   0 ] [ 0  U22 U23 U24|	**
**  |A31 A32 A33 A34] = [L31 L32  1   0 ] [ 0   0  U33 U34|	**
**  [A41 A42 A43 A44]   [L41 L42 L43  1 ] [ 0   0   0  U44]	**
**								**
**  Neglecting the diagonal of L, both L and U can be overlaid	**
**  for storage in a single array (here overwriting A).		**
**								**
**                      [U11 U12 U12 U14]			**
**  			[L21 U22 U23 U24|			**
**  			[L31 L32 U33 U34|			**
**  			[L41 L42 L43 U44]			**
**								**
**  Finally, the algorithm uses Crouts method to determine the	**
**  elements of L and U.  This sets the diagonal elements of L	**
**  to 1 and follows a specific order in computation of the L	**
**  and U elements which allows A to be overwritten as L and U	**
**  are determined.  For Crouts method to be stable, partial	**
**  pivoting is used (row interchanges).			**
**								**
*****************************************************************/


	// Input		A     : Input normal matrix (this)
	//			indx  : Integer array for permutations
        // Return		A     : LU Decomposition of input A
        // 			indx  : Index of Permutations
	// Note			      : Assumes A is square
	// Note			      : Input matrix is overwritten!

int n_matrix::LU_decomp(int* indx)
  {
  int n = rows();			// Get dimension of A
  const double tiny = 1.0e-20;		// Set this to a small number
//  double vv[n];				// Vector for row scaling factors
  double *vv;
  vv = new double[n];
  int d = 1;				// Permutation counter (1 even, -1 odd)

  double big, temp;
  int i=0;
  for(i=0; i<n; i++)			// For each row of A find the
    {					// the largest element for pivoting 
    big=0.0;				// and store corresponding scaling
    for(int j=0; j<n; j++)		// factor in array vv
      {
      temp = norm(get(i,j));
      if(temp > big)
        big = temp;
      }
    if(fabs(big) < 1.e-9)		// Cannot have a zero row in A!
//    if(big == 0.0)			// Cannot have a zero row in A!
      {
      std::cout << "\nClass n_matrix: Singular matrix input\n";
      std::cout << "\nClass n_matrix: Cannot invert the array\n";
      exit(-1);
      }
    vv[i]=1.0/big;
    }

  complex Uij, Lij;			// Element of U and L
  int imax=0;
  double dum;
  complex dumz;
  for(int j=0; j<n; j++)		// Begin going through each A column
    {

// First implement equation below to compute the values of U (except <i|U|i>)
//
//		       i-1
//		       ---
// <i|U|j> = <i|A|j> - \   <i|L|k><k|U|j>
// 		       /
//		       ---
//		       k=0
//
// This is equation 2.3.12 on page 34 of the referenced text and its use is
// restricted to terms above the diagonal.

    for(i=0; i<j; i++)
      {
      Uij = get(i,j);
      for(int k=0; k<i; k++)
        Uij -= get(i,k)*get(k,j);
      put(Uij,i,j);
      }

// Next implement the equation below to compute the values of U (i>j, i<n)
// For the case when i=j, this equation is the same as thh previous equation
// when the scaling factor is neglected (and this is done).
//
//			       j-1
//	        1	       ---
// <i|L|j> = ------- <i|A|j> - \   <i|L|k><k|U|j>
// 	     <j|U|j>	       /
//			       ---
//			       k=0

    big=0.0;				// Used to look for largest pivot
    for (i=j; i<n; i++)
      {
      Lij = get(i,j);
      for(int k=0; k<j; k++)
	Lij -= get(i,k)*get(k,j);
      put(Lij,i,j);
      dum = vv[i] * norm(Lij);		// For looking at the pivot
      if (dum >= big)			// Use it if better that big
        {
	big=dum;
	imax=i;
	}
      }

    if (j != imax)			// This part interchanges the
      {					// rows of A via permutation
      for(int k=0; k<n; k++)		// The sign of d is switched
        {				// and the permutation stored
        dumz = get(imax,k);		// in the array indx
        (*this)(imax,k) = get(j,k);
	put(dumz,j,k);
	}
      d = -d;
      vv[imax] = vv[j];
      }

    indx[j] = imax;
    if(norm(get(j,j)) == 0.0)		// If <j|A|j> is zero make it
      put(tiny,j,j);			// small instead!
    if (j != n-1)			// Divide by the pivot element
      {					// If <j|A|j> was zero it has 
      dumz = 1.0/get(j,j);		// since been set to tiny so it
      for(i=j+1; i<n; i++)		// won't blow up!
        (*this)(i,j) *= dumz;
      }
    }					// Return to treat next column
  delete [] vv;
  return d;
  }


/* LU Back Substitution *********************************-*-c++-*-
**							 	**
**  This routine is taken from "Numerical Recipies", W.H.Press,	**
**  B.P.Flannery, S.A.Teukolsky, & W.T.Vetterling, Cambridge	**
**  University Press, 1989.  See pages 31 through 38.  The code **
**  has been modified for GAMMA into C++.			**
**							 	**
**  The algorithm takes an input matrix A' in its LU decomp-	**
**  osition format where A' is some original matrix A with any	**
**  needed row permutations to attain the LU form. The row	**
**  permutations used to relate A to A' are stored in the	**
**  integer array indx.  The function the proceeds to solve	**
**								**  
**			A|x> = |b>				**
**								**  
**  for |x> by considering the problem as			** 
**								**  
**              A|x> = (LU)|x> = L(U|x>) |b> 			**
**								**  
**  and first solving for the vector |y> where |y> = U|x>	**
**								**  
**			L|y> = |b>				**
**								**
**  followed by solving for |x> by				**
**								**
**			U|x> = |y>				**
**								**
**  Due to the triagular nature of the arrays L and U it is	**
**  relatively easy to solve the vector equations involving	**
**  them.  The only further complexity in this algorithm is	**
**  that the input |b> must be first permuted to |b'> so that	**
**  its element ordering match that of A'.  Following the	**
**  solution in which |x'> is obtained it is then un-permuted	**
**  to match the original A ordering.  The first equation	**
**  to be solved (actually involving L' not L) appears as	**
**								**
** 		 [ 1   0   0   0 ] [y1] = [b1]			**
** 		 [L21  1   0   0 ] [y2] = |b2|			**
** 		 [L31 L32  1   0 ] [y3] = |b3|			**
** 		 [L41 L42 L43  1 ] [y4] = [b4]			**
**								**
**  because the diagonal elements of the matrix L are always 1.	**
**  The first element of |y> (actually |y'>) is given by	**
**								**
**			y1 = b1					**
**								**
**  and then subseqent elements of this vector are found by	**
**								**
**			      [	     ---	  ]		**
**		         1    |      \	 	  |		**
**	         yi = ------- | b  - /  <i|L|j>y  |		**
**		      <i|L|i> [  i   ---	j ]		**
**								**
*****************************************************************/


	// Input		ALU   : Input matrix A in LU form (this)
	//			indx  : Integer array for permutations
	//			b     : Input column vector b (normal matrix)
        // Return		b     : The vector |x> of A|x> = |b>
	// Note			      : The input matrix ALU & the vector
	//				indx are not modified by this routine

void n_matrix::LU_backsub(int* indx, n_matrix& b)
  {
  int n = rows();			// Get dimension of array
  int i=0, ii=0, ip=0;
  complex sum=0;

//		First Solve L|y> = |b> For |y>

  for(i=0; i<n; i++)			// Works forwards, |y(0)>
    {					// solved first because L
    ip = indx[i];			// is lower triangular
    sum = b.get(ip,0);			// A is input permuted! Get
    b(ip,0) = b.get(i,0);		// intial b's in correct order
    if (ii)
      for (int j=0; j<i; j++)
        sum -= get(i,j)*b.get(j,0);
    else if(norm(sum))
      ii=1;
    b.put(sum,i,0);
    }

//		Now Solve U|x> = |y> For |x>

  for(i=n-1; i>=0; i--)			// Works backwards, |x(n-1)>
    {					// solved first because U
    sum = b.get(i,0);			// is upper-triangular
    for(int j=i+1; j<n; j++)
      sum -= get(i,j)*b.get(j,0);
    b(i,0) = sum/get(i,i);		// This is now |x(i)>
    }
  return;
  }


/* Matrix Inversion *************************************-*-c++-*-
**							 	**
**  This routine inv combines the decomposition of an input	**
**  matrix into its LU form with repetitive backsubstitution	**
**  to form the matrix inverse.  The problem can be formulated	**
**  as								**
**		     	   A * Ainv = I				**  
**								**  
**  where A (input matrix) and I are known and the equation	**
**  solved for Ainv.  The LU decomposition of A is described	**
**  by								**
**		     	        A = LU				**  
**								**  
**  where L is lower triangular and U is upper triangular.	**
**								**
**  [A11 A12 A13 A14]   [L11  0   0   0 ] [U11 U12 U12 U14]	**
**  |A21 A22 A23 A24]   [L21 L22  0   0 ] [ 0  U22 U23 U24|	**
**  |A31 A32 A33 A34] = [L31 L32 L33  0 ] [ 0   0  U33 U34|	**
**  [A41 A42 A43 A44]   [L41 L42 L43 L44] [ 0   0   0  U44]	**
**								**
**  Both L and U have elements on the diagonal, but those of L	**
**  are always set to be 1: <i|L|i> = 1.			**
**								**
**  [A11 A12 A13 A14]   [ 1   0   0   0 ] [U11 U12 U12 U14]	**
**  |A21 A22 A23 A24]   [L21  1   0   0 ] [ 0  U22 U23 U24|	**
**  |A31 A32 A33 A34] = [L31 L32  1   0 ] [ 0   0  U33 U34|	**
**  [A41 A42 A43 A44]   [L41 L42 L43  1 ] [ 0   0   0  U44]	**
**								**
**  Substitution of the LU decomposition into the original	**
**  equation produces						**
**								**  
**			LU * Ainv = I				**
**								**  
**  or equivalently						**
**								**  
**	     L * A' = I     where     A' = U * Ainv		**
**								**  
**  Since L is lower triangular it is easy to backsolve for	**
**  the matrix A'.  Once A' has been found the equation of	**
**  importance becomes						**
**								**  
**	     		A' = U * Ainv				**
**								**
**  and again, because U is upper triangular, it is trivial to	**
**  backsolve  for the desired matrix Ainv.			**
**								**
**								**
*****************************************************************/

	// Input		A     : Input normal matrix (this)
        // Return		Ainv  : Inverse of input A
	// Note			      : A is unaltered in this routine
	// Note			      : Assumes A is square!
// sosi: 6/23/92 - THIS STILL HAS TROUBLE FLAGGING A SINGULAR MATRIX
//	           FOR EXAMPLE TRY matrix(10,10,5.0);

_matrix* n_matrix::xinv()

  {
  complex z(0,0);
  int nr = rows();			// Get dimension of A
  n_matrix* Ainv = new n_matrix(nr,nr); // Construct matrix for inverse
  int* indx;
  indx = new int[cols()];
  n_matrix ALU(*this);			// Copy matrix A for LU decomp
  ALU.LU_decomp(indx);			// LU Decomposition of A
  n_matrix Ii(nr, 1, z);		// Column vector |Ii>
  n_matrix Ainvi(nr, 1, z);		// Column vector of inverse |Ainvi>
  int i=0, j=0;
  for(i=0; i<nr; i++)			// Go through the columns of I,    -1
    {					// |Ii>, solve for the columns of A  .
    Ii.data[i] = 1;			// Form column of I, |Ii>
    Ainvi = Ii; 			//		           -1 
    ALU.LU_backsub(indx, Ainvi);	// Get column of Ainv: A|Ai  > = |Ii>
    for(j=0; j<nr; j++)			// Copy Ainv column to Ainv
      (*Ainv)(j,i) = Ainvi.data[j];
    Ii.data[i] = 0;			// Unform column of I, |Ii>
    }
  delete [] indx;
  return Ainv;				// Return Ainv, a new n_matrix
  }

#endif						    // n_matrix.cc
