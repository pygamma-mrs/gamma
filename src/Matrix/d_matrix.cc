/* d_matrix.cc **************************************************-*-c++-*-
**								   	**
**      	               G A M M A				**
**								   	**
**	Diagonal Matrix                         Implementation		**
**						   			**
**	Copyright � 1990, 1991, 1992				   	**
**	Tilo Levante & Scott A. Smith				   	**
**	Eidgenoessische Technische Hochschule				**
**	Labor fuer physikalische Chemie				   	**
**	8092 Zuerich / Switzerland		 		  	**
**								   	**
**      $Header: $
**								   	**
*************************************************************************/

/*************************************************************************
**								   	**
** Description							   	**
**								   	**
** The class d_matrix defines diagonal complex matrices for C++ with 	**
** the usual operations +, -, *, / and I/O routines.			**
**								   	**
*************************************************************************/
     
#ifndef   Gd_matrix_cc_				// Is file already included?
#  define Gd_matrix_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#  endif

#include <GamGen.h>				// Include OS specific stuff
#include <Matrix/d_matrix.h>			// Include the interface
#include <Matrix/_matrix.h>             	// Include the base matrix class
#include <Matrix/n_matrix.h>			// Know about normal matrices
#include <Matrix/i_matrix.h>			// Know about identity matrices
#include <Matrix/h_matrix.h>			// Know about Hermitian matrices
#include <Matrix/MxModBas.h>			// Include Matrix module errors
#include <fstream>				// Know libstdc++ filestreams
#include <cmath>				// Know max and min functions
#include <Basics/StringCut.h>

// ____________________________________________________________________________
// i                     CLASS D_MATRIX ERROR HANDLING
// ____________________________________________________________________________

/*      Input               mx      : A normal complex matrix (this)
                            eidx    : Error index
                            noret   : Flag for linefeed (0=linefeed)
                            pname   : string in message
        Output              none    : Error message output
                                      Program execution stopped if fatal

  Note that these functions just route into the error functions in the base
  class _matrix.                                                             */

void d_matrix::DMxerror(int eidx, int nr) const
  {
  std::string CL="Diagonal Matrix";
  Mxerror(CL, eidx, nr);
  }

void d_matrix::DMxerror(int eidx, const std::string& PN, int nr) const
  {
  std::string CL="Diagonal Matrix";
  Mxerror(CL, eidx, PN, nr);
  }
 
volatile void d_matrix::DMxfatal(int eidx) const
  {
  std::string CL="Diagonal Matrix";
  Mxfatality(CL, eidx);
  }
 
volatile void d_matrix::DMxfatal(int eidx, const std::string& PN) const
  {
  std::string CL="Diagonal Matrix";
  Mxfatality(CL, eidx, PN);
  }
     
// ____________________________________________________________________________
// ii                      CLASS D_MATRIX CHECKING
// ____________________________________________________________________________


bool d_matrix::CheckDims(_matrix* mx, int warn)
  {
  if(cols_!=mx->rows() || cols_!=mx->cols())	// Insure mx dimension match
    {
    if(warn>0)
      {
      DMxerror(51,1);				// Array dimensions mismatched
      DMxerror(31,1);				// Row->Row Col-Col mismatch
      if(warn>1) DMxfatal(81);			// Cannot continue on
      }
    return false;
    }
  return true;
  }

bool d_matrix::CheckDim(_matrix* mx, bool mul, int warn)
  {
  bool TF;
  if(mul) TF = (cols_ == mx->rows());		// Insure dim. proper for mult 
  else    TF = (cols_ == mx->cols());		// Insure dim. proper for div
  if(!TF)
    {
    if(warn>0)
      {
      DMxerror(51,1);				// Array dimensions mismatched
      DMxerror(30,1);				// Row->Col Col-Row mismatch
      if(warn>1) DMxfatal(81);			// Cannot continue on
      }
    return false;
    }
  return true;
  }

// ____________________________________________________________________________
// A                 CLASS D_MATRIX CONSTRUCTORS AND DESTRUCTOR
// ____________________________________________________________________________

/*  Arguments                             Constructed Array
              ---------         -------------------------------------
                  -	 	 An empty array
               nr, nc		 An (nr x nc) array, uninitialized
              nr, nc, z		 An (nr x nc) array, sets <i|mx|i>=z
                 dmx		 A duplicate of diagonal array dmx

   Note that all constructors invariably call the base class constructor
   (i.e. the one in _matrix) which sets the row & column sizes.  The added
   code should simply allocate the data array of proper size. Similarly, the
   destructor only destroys the data array, anything else is destroyed in 
   the base class.                                                           */

d_matrix::d_matrix( ) : _matrix() { data = new complex[size]; }
d_matrix::d_matrix(int i, int j) : _matrix(i,j)	
  {
  if(i!=j) 
    { 						// Insure dmx is square
    DMxerror(9, 1);                             // Problems in constructor
    DMxfatal(50);                               // Rectangular array
    }
  size = i;					// Set the data size
  data = new complex[size];			// Set up data storage
  }

d_matrix::d_matrix(int i, int j, const complex& z) : _matrix(i,j)
  {
  if(i!=j) 
    { 						// Insure dmx is square
    DMxerror(9, 1);                             // Problems in constructor
    DMxfatal(50);                               // Rectangular array
    }
  size = i;					// Set the data size
  data = new complex[size];			// Set up data storage
  for(int pos=0; pos<size; pos ++)		// Set diagonal to z
    data[pos] = z;
  }

d_matrix::d_matrix(const d_matrix& dmx) : _matrix(dmx)
  {
  size = dmx.size;				// Set sizes the same
  data = new complex[size];			// Set up data storage
  for(int pos=0; pos<size; pos ++)		// Copy dmx matrix data
    data[pos] = dmx.data[pos];
  real = dmx.real;
  imag = dmx.imag;
  }

d_matrix::~d_matrix() { delete [] data; }
 
// ____________________________________________________________________________
// B                  CLASS D_MATRIX ACCESS AND ASSIGNMENT
// ____________________________________________________________________________

/* There are two flavors of element access operators herein. The one that might
   commonly be used, i.e. mx(i,j), is DIRECTLY accessing the element <i|mx|j>.
   That can save CPU time by avoiding value copying when altering or obtaining
   elements.  However its misuse can lead to trouble in this class because not
   all diagonal matrix elements are stored and there are restrictions on the
   elements so that the array remains diagonal.  An example of when this is
   trouble would be "dmx(1,0)=z;".  This fails because off-diagonal
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
          =             Current (diagonal) array set equal to input dmx
    (int,int)           Reference to element <i|dmx|j> (Potential Danger)
    get(int,int)        Copy of element <i|dmx|j>      (Safe)
    put(int,int)        Assigns element <i|dmx|j>      (Return FALSE if fails)
    put_h(int,int)      Assigns both <i|dmx|j> & <j|dmx|i> (Return TRUE/FALSE)
    get_block(r,c,R,C)  Returns dmx block of size RxC starting at <r|dmx|c>
    put_block(r,c,mx)   Places mx into dmx at position <r|dmx|c> (TRUE/FALSE)
 
    Again, note that operations which return TF will trigger class matrix to
    take appropriate action when they return FALSE.  This usually means that
    the operation cannot be handled by a diagonal array so the array type
    is changed to be more generic.  The the operation is reattempted.       */

void d_matrix::operator= (const d_matrix& mx)
  {
  if(this == &mx) return;			// Do nothing if equal
  delete [] data;				// Delete original data
  _matrix::operator= (mx);			// Now construct new mx
  data = new complex[size];
  for(int pos=0; pos<size; pos ++)		// Set elements equal
    data[pos] = mx.data[pos];
  }

complex& d_matrix::operator() (int i, int j)
  { return (i==j)?data[i]:Zdzero; }		// Return <i|dmx|i> or 0

complex d_matrix::get(int i, int j) const
  { return (i==j)?data[i]:complex0; }		// Return <i|dmx|i> or 0

bool d_matrix::put(const complex& z, int i, int j)
  {
  if(i==j) data[i] = z;				// Don't set unless on diag
  return(i==j);					// Else handle in matrix
  }

bool d_matrix::put_h(const complex& z, int i, int j) 
  {
  if(i==j) data[i] = z;				// Don't set unless on diag
  return(i==j);					// Else handle in matrix
  }

_matrix* d_matrix::get_block(int row, int col, int nrows, int ncols)
  {
  if(row==col && nrows==ncols)			// If square block about diagonal
    {						// then return another diagonal array
    if(!row && nrows==cols()) return this;	// If block is whole matrix, return dmx
    d_matrix* dmx = new d_matrix(nrows, ncols);	// Else return a smaller d_matrix
    for(int i=0; i<nrows; i++)
      dmx->data[i] = data[i+row];
    return dmx;
    }
  else						// If block rectangular or not about diagonal
    {						// then must return a normal matrix
    int i;
    n_matrix* nmx = new n_matrix(nrows,		// Return a n_matrix if block is
		            ncols, complex0);	// rectangular or off-diagonal
    if(row == col)				// If diagonals of dmx and nmx
      { 					// line up, just copy diagonal
      int end = (ncols>nrows)?nrows:ncols;
      for(i=0; i<end; i++)			// Loop over rows or columns 
        ((n_matrix*)nmx)->data[i*ncols+i]	// <i|nmx|i> = <i+col|dmx|i+col>
	                         = data[i+col];
      }
    else					// Getting a block of no particularn 
      {
      int j, I, J;
      for(i=0; i<nrows; i++)
        for(j=0,I=i+row,J=col; j<ncols; j++,J++)
          if(I == J)				// <i|nmx|j> = <i+row|nmx|j+col>
            ((n_matrix*)nmx)->data[i*ncols+j]	//           = <I|dmx|J> = del   * <I|dmx|I>
	                             = data[I];	//                            I,J
      }
    return nmx;
    }
  }

bool d_matrix::put_block(int row, int col, _matrix* mx)
  {
  if(row != col) return false; 			// Can't handle off-diag. put
  matrix_type mt = mx->stored_type();		// Get the input array type
  if(mx->rows()==1 && mx->cols()==1)		// If put involves only 1 
    { data[row]=(*mx).get(0,0); return true; }  // one element we handle it
  if(mt == h_matrix_type) return false;		// Can't handle hmx block put
  if(mt == n_matrix_type) return false;		// Can't handle nmx block put
  if(row+mx->rows() > rows_)	 		// Too big, so it won't work!
    {
    DMxerror(5, "put_block", 1);                // Bad use of function
    DMxfatal(52);                               // Array dimensions exceeded
    }
  if(mt == i_matrix_type)			// We can deal with putting
    {						// imx on diagonal. But not if
    if(mx->rows() == rows()) return false;      // it spans dmx, need dmx->imx 
    for(int i=0; i<mx->rows(); i++)		// Else just set matrix elements
      data[row+i] = 1;				// to 1 (from imx diagonals)
    return true; 
    }
  if(mt == d_matrix_type)			// We can deal wiht putting
    {						// dmx on the diagonal. Here we
    for(int i=0; i<mx->rows(); i++)		// just set matrix elements
      data[row+i] = ((d_matrix*)mx)->data[i];	// to those of input mx
    return true; 
    }
  return false;
  }
     
// ____________________________________________________________________________
// C                  CLASS D_MATRIX HERMITIAN_TYPE HANDLING
// ____________________________________________________________________________
/* These functions handle the checking of whether an array is of a Hermitian 
   type.  Since GAMMA has a Hermitian matrix type (h_matrix) it is easy to 
   know the whether the array is stored Hermitian or not. The other functions 
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
 
hermitian_type d_matrix::stored_hermitian( ) const { return non_hermitian; }

hermitian_type d_matrix::test_hermitian(double d) const
  {
  int f = 1;
  for(int r=0; (r<size)&&f; r++ )
    f = (fabs(Im(data[r])) < d);
  return f?_hermitian:non_hermitian;
  }
 
// ____________________________________________________________________________
// D                     CLASS D_MATRIX MATRIX_TYPE HANDLING
// ____________________________________________________________________________
 
/* These functions handle the checking of what type of array we have.  These
   types directly correspond to the matrix classes (such as d_matrix) derived
   from class _matrix.

     Function            Output                       Description
   ----------------  --------------  ------------------------------------------
   set_type            --------      Found in class matrix
   test_type           --------      Found in class matrix
   stored_type       d_matrix_type   Always returns we are diagonal
   test_hermitian    *_matrix_type   Type dmx could be within d
   mxtype            string          Returns the string "Diagonal"
 
   The test type looks to see if dmx could be an array of the matrix type
   specified as an input argument.  If it can within "d" the return type
   is equal to the input type.  If it cannot the return is d_matrix_type     */  

// Note:  Functions set_type & check_type are found in class matrix, not here.

matrix_type d_matrix::stored_type ( ) const { return d_matrix_type; }
matrix_type d_matrix::test_type(matrix_type m, double d) const
  {
  int f, i;
  switch(m)
    {
    case n_matrix_type: return n_matrix_type; break;
    case d_matrix_type: return d_matrix_type; break;
    case i_matrix_type:			// See if it can be stored i_matrix
     f = 1;
     for(i=0; (i<rows_)&&f; i++)	// Check for non-unity diagonals
       {
       if(fabs(Re(data[i]-complex1)) < d)
         f = (fabs(Im(data[i])) < d);
       else f = 0;
       }
     return f?i_matrix_type:d_matrix_type;
     break;
   case h_matrix_type:			// See if it can be stored h_matrix
     f = 1;
     for(i=0; (i<rows_)&&f; i++)	// Check for imaginary components
       f = (Im(data[i]) < d);		// on the diagonal
     return f?h_matrix_type:d_matrix_type;
     break;
   default:
     return d_matrix_type;
     break;
    }
  }

std::string d_matrix::mxtype()        const { return std::string("Diagonal"); }
std::string d_matrix::mxtype(bool pf) const 
  {
  if(!pf) return std::string("Diagonal");
  else if(is_real())      return std::string("Real Diagonal");
  else if(is_imaginary()) return std::string("Imaginary Diagonal");
  return std::string("Diagonal");
  }

// ____________________________________________________________________________
// E                 CLASS D_MATRIX VARIOUS MATRIX CHECKS
// ____________________________________________________________________________
 
/*   Function    Output                       Description
   ------------  ------  ------------------------------------------------------
   is_unitary     bool   TF if inv(_mx) == adjoint mx, CPU intensive
   is_real        bool   TF if Im(<i|_mx|j>) < d for all i,j (def. d GMxCut)
   is_imaginary   bool   TF if Re(<i|_mx|j>) < d for all i,j (def. d GMxCut)
   is_complex     bool   TF if is_real && is_imaginary       (def. d GMxCut)
   is_zero        bool   TF if ||<i|_mx|j>|| < d for all i,j (def. d GMxCut)
   is_equal       bool   TF if ||<i|dmx-mx|j>||<d    all i,j (def. d GMxCut)
 
                                                 *
  Unitary test is:       z(i,i) = <i|dmx|i><i|dmx |i> = 1                    */

bool d_matrix::is_symmetric(double d) const { return true; }
bool d_matrix::is_hermitian(double d) const { return is_real(d); }

bool d_matrix::is_unitary(double d) const
  {
  bool f = true;
  complex z;
  for(int r=0; (r<size)&&f; r++ )
    {
    z = data[r]*conj(data[r]);
    if((fabs(Re(z)-1)) < d)
      f = (fabs(Im(z)) < d);
    else f = false;
    }
  return f;
  }
 
bool d_matrix::is_real(double d) const
  {
  bool f = true;
  for(int r=0; (r<size)&&f; r++ )
    f = (fabs(Im(data[r])) < d);
  return f;
  }

bool d_matrix::is_imaginary(double d) const
  {
  bool f = true;
  for(int r=0; (r<size)&&f; r++ )
    f = (fabs(Re(data[r])) < d);
  if(f) f = (!is_zero(d));		// Insure this isn't a zero mx
  return f;
  }

bool d_matrix::is_complex(double d) const
  {
  bool f=true;
  if(is_imaginary(d)) f = 0; 		// False if pure imaginary
  else if(is_real(d)) f = 0;		// False if pure real
  return f;
  }

bool d_matrix::is_zero(double d) const
  {
  bool f = true;
  for(int r=0; (r<size)&&f; r++ )
    f = (norm(data[r]) < d);
  return f;
  }

bool d_matrix::is_diagonal(double d) const { return true; }
bool d_matrix::is_square()           const { return true; }
 
bool d_matrix::is_equal(_matrix* mx, double d) const
  {
  if(cols_ != mx->cols())                                return false;
  if(rows_ != mx->rows())                                return false;
  matrix_type mst = (*mx).test_type(d_matrix_type, d);
  if(mst != d_matrix_type)                               return false;
  bool f = true;
  for(int r=0; (r<size)&&f; r++ )
    f = (norm(data[r]-(*mx)(r,r)) < d);
  return f;
  }

// ____________________________________________________________________________
// F                CLASS D_MATRIX BINARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________

/* These functions allow for simple mathematical operations between a diagonal
   matrix and another matrix of unspecified type. Note that these functions
   return a pointer to the resulting array. That array must have it referencing
   done in the main matrix class or there will be memory leakage.

   Function    Result                   Details
   --------   --------   ---------------------------------------------------
     add       dmx+mx    Adjusts to mx type, sums only over active elements
   subtract    dmx-mx    Adjusts to mx type, subtact only over active elems.
   multiply    dmx*mx    Reduced calculation using only nonzero dmx elements
   multiply    dmx*z     Fast calculation only over dmx stored elements

*/

_matrix* d_matrix::add(_matrix* mx)
  {
  if(!CheckDims(mx, 1)) DMxfatal(20);			// Quit if mismatch
  switch (mx->stored_type())
    {
    case n_matrix_type:				// Add a n_matrix to a d_matrix
      { 
	n_matrix* sum = new n_matrix(*(n_matrix*)mx);	// <i|sum|j> = <i|mx|j>
	for(int i=0; i<cols_; i++)			// <i|sum|i> += <i|dmx|i>
          ((n_matrix*)sum)->data[i*cols_+i] += data[i];
	return sum;
      }
      break;
    case h_matrix_type:				// Add a h_matrix to a d_matrix
      { 
      if(is_real())					// If dmx is real, result Hermitian
        {
	  h_matrix* sum = new h_matrix(*(h_matrix*)mx);	// <i|sum|j> = <i|mx|j>
	  for(int i=0; i<cols_; i++) 			// <i|sum|i> += <i|dmx|i>
	    sum->data[i*cols_-(i*(i-1))/2] += data[i];
	  return sum;
	  }
      else						// If dmx is complex, result is normal
        {
	n_matrix* sum = new n_matrix(cols_,cols_);	// Generate new n_matrix for sum
        complex *dii = data;				// Element of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *hij = ((h_matrix*)mx)->data;		// Element of hmx: hij = <i|hmx|j> -> <0|hmx|0>
        complex *sij = ((n_matrix*)sum)->data;		// Element of sum: sij = <i|sum|j> -> <0|sum|0>
        complex *sji;					// Element of sum: sji = <j|sum|i>
        int i,j;
	for(i=0; i<cols_; i++,sij+=i,dii++)		// Loop over rows i, dii=<i|dmx|i>, sij=<i|sum|i> 
          {
          *sij = *dii + *hij;				// <i|sum|i> = <i|dmx|i> + <i|hmx|i> 
          sji = sij + cols_;				// <j|sum|i> -> <i+1|sum|i>
	  for(j=i+1,hij++,sij++; j<cols_;		// Loop over upper-triangle elements
                        j++, hij++,sij++,sji+=cols_)	// Start at hij=<i|hmx|i+1>, sij=<i|sum|i+1>
            {
	    *sij = *hij;	 			// <i|sum|j> = <i|mx|j>, i<j
	    *sji = conj(*hij);	 			// <j|sum|i> = <i|mx*|j>, i<j
            }
          }
	return sum;
	}
      }
	break;
    case d_matrix_type:				// Add two d_matrices
      { 
	d_matrix* sum = new d_matrix(*(d_matrix*)mx); 	// <i|sum|j> = <i|mx|j>
	for(int i=0; i<cols_; i++)			// Add in the diagonals of dmx
	  sum->data[i]+=data[i];			// <i|sum|i> += <i|dmx|i>
	return sum;
      }
	break;
    case i_matrix_type:				// Add an i_matrix to a d_matrix
      { 
	d_matrix* sum = new d_matrix(*this);		// Generate a new d_matrix
	for(int i=0;  i<cols_; i++)			// <i|sum|j> = <i|dmx|j>
	  sum->data[i] += complex1;			// <i|sum|i> += 1
	return sum;
      }
	break;
    default:						// Add a generic matrix to a d_matrix
      { 
	n_matrix* sum = new n_matrix(cols_,cols_);	// Construct new normal matrix
	int pos = 0;					// Copy all the elements
	for(int i=0;  i<cols_; i++)
	  for(int j=0; j<cols_; j++,pos++)
	    if(i==j)
	      sum->data[pos] = data[i] + (*mx)(i,j);	// <i|sum|j> = <i|dmx|i>+<i|mx|j> i=j
	    else
	      sum->data[pos] = (*mx)(i,j);		// <i|sum|j> = <i|mx|j> i!=j
	return sum;
      }
    }
  }


_matrix* d_matrix::subtract(_matrix* mx)
  {
  if(!CheckDims(mx, 1)) DMxfatal(21);			// Quit if mismatch
  switch (mx->stored_type())
    {
    case n_matrix_type:				// Subtract n_matrix from a d_matrix
	{ 
	n_matrix* dif = (n_matrix *)mx->negate();	// <i|dif|j> = -<i|mx|j>
	for(int i=0; i<cols_; i++)			// <i|dif|i> += <i|dmx|i> 
          ((n_matrix*)dif)->data[i*cols_+i] += data[i];
	return dif;
	}
	break;
    case h_matrix_type:				// Subtract h_matrix from a d_matrix
	{ 
        if(is_real())					// If dmx is real, result Hermitian
          {
	  h_matrix* dif = (h_matrix *)mx->negate();	// <i|dif|j> = -<i|mx|j>
	  for(int i=0;  i<cols_; i++) 			// <i|dif|i> += <i|dmx|i>
	    dif->data[i*cols_-(i*(i-1))/2]+=data[i];
	  return dif;
	  }
        else						// If dmx is complex, result is normal
          {
	  n_matrix* dif = new n_matrix(cols_,cols_);	// Generate new n_matrix for dif
          complex *dii = data;				// dii = <i|dmx|i> -> <0|dmx|0>
          complex *hij = ((h_matrix*)mx)->data;		// hij = <i|hmx|j> -> <0|hmx|0>
          complex *sij = ((n_matrix*)dif)->data;	// sij = <i|dif|j> -> <0|dif|0>
          complex *sji;					// sji = <j|dif|i>
          int i,j;
	  for(i=0; i<cols_; i++,sij+=i,dii++)		// Loop over rows i, dii=<i|dmx|i>, sij=<i|dif|i> 
            {
            *sij = *dii - *hij;				// <i|dif|i> = <i|dmx|i> - <i|hmx|i> 
            sji = sij+cols_;				// <j|dif|i> -> <i+1|dif|i>
	    for(j=i+1,hij++,sij++; j<cols_;		// Loop over upper-triangle elements
                          j++, hij++,sij++,sji+=cols_)	// Start at hij=<i|hmx|i+1>, sij=<i|dif|i+1>
              {
	      *sij = -(*hij);	 			// <i|dif|j> = -<i|mx|j>, i<j
	      *sji = -conj(*hij);	 		// <j|dif|i> = -<i|mx*|j>, i<j
              }
            }
	  return dif;
	  }
	}
	break;
    case d_matrix_type:				// Subtract two d_matrices
	{ 
	d_matrix* dif = (d_matrix *)mx->negate();	// <i|dif|j> = -<i|mx|j>
	for(int i=0;  i<cols_; i++)			// <i|dif|i> += <i|dmx|i>
	  dif->data[i] += data[i];
	return dif;
	}
	break;
    case i_matrix_type:				// Subtract an i_matrix from a d_matrix
	{ 
	d_matrix* dif = new d_matrix(*this);		// <i|dif|j> = <i|dmx|j>
	for(int i=0;  i<cols_; i++)			// <i|dif|i> -= 1
	  dif->data[i] -= complex1;
	return dif;
	}
	break;
    default:						// Subtract generic matrix from d_matrix
	{ 
	n_matrix* dif = new n_matrix(cols_,cols_);	// Construct new normal matrix
	int pos = 0;
	for(int i=0;  i<cols_; i++)
	  for(int j=0; j<cols_; j++,pos++)
	    if(i==j)
	      dif->data[pos] = data[i]-(*mx)(i,j);	// <i|dif|j> = <i|dmx|i>-<i|mx|j>, i==j
	    else
              dif->data[pos] =- (*mx)(i,j);		// <i|dif|j> = -<i|mx|j>, i!=j
	  return dif;
        }
    }
  }   


_matrix* d_matrix::multiply(_matrix* mx)
 
        // Input                dmx  : Input diagonal matrix (this)
        //                       mx  : Second matrix
        // Output               pdt  : New matrix which is the product dmx*mx
        // Note                      : Uses internal structure of other matrix types
 
//                        ---
//            <i|pdt|j> = \   <i|dmx|k> <k|mx|j> = <i|dmx|i><i|mx|j> 
//                        /
//                        ---

  {
  if(!CheckDim(mx)) 					// Quit if mismatch
    { DMxerror(5,"multiply",1); DMxfatal(22); }
  switch (mx->stored_type())
      {
      case d_matrix_type:				// Multiply two d_matrices
	{ 
	d_matrix* pdt = new d_matrix(rows_,cols_);	// Set pdt to new d_matrix
	for(int i=0; i<rows_; i++)
	  pdt->data[i] = this->data[i]			// <i|pdt|j> = del  <i|dmx|i><i|mx|i>
                          * ((d_matrix*)mx)->data[i];	//                ij
	return pdt;
	}
	break;
      case i_matrix_type:				// Multiply d_matrix into i_matrix
	return this;					// Product is dmx, return it
	break;
      case h_matrix_type:				// Multiply d_matrix into h_matrix
	{ 
	n_matrix* pdt = new n_matrix(rows_,cols_);	// Construct result as n_matrix
        complex *p00 = pdt->data;			// Start of pdt: p00 = <0|pdt|0>
        complex *d00 = data;				// Start of dmx: d00 = <0|dmx|0>
        complex *h00 = ((h_matrix*)mx)->data;		// Start of hmx: h00 = <0|hmx|0>
        int i=0,j;
        complex *pij, *dii, *hji, *hij;
	for(i=0,pij=p00,dii=d00; i<rows_; i++,dii++)	// Effective loop over i
          {
          hji = h00+i;					// Set <j|hmx|i> to <0|hmx|i>
	  for(j=0; j<i; j++,pij++,hji+=cols_-j)
	    (*pij) = (*dii) * conj(*hji);		// <i|pdt|j> = <i|dmx|i><j|hmx*|i>, i>j
	  for(hij=hji; j<cols_; j++,pij++,hij++)
	    (*pij) = (*dii) * (*hij);			// <i|pdt|j> = <i|dmx|i><i|hmx|j>, i<=j
          }
	return pdt;
        break;
	}
      case n_matrix_type:	             		// Multiply d_matrix into n_matrix
	{ 
        int c = mx->cols();				// Get pdt column dimension
	n_matrix* pdt = new n_matrix(rows_,c);		// Construct result as n_matrix
        complex *p00 = pdt->data;			// Start of pdt: p00 = <0|pdt|0>
        complex *d00 = data;				// Start of dmx: d00 = <0|dmx|0>
        complex *n00 = ((n_matrix*)mx)->data;		// Start of nmx: h00 = <0|nmx|0>
        complex *pend = p00 + rows_*c;			// Start of nmx: h00 = <0|nmx|0>
        complex *pij, *dii, *nij;			// These will be matrix elements
        complex *nrow;					// Dummy index for end of nmx row
	for(pij=p00,nij=n00,dii=d00; pij<pend; dii++)	// Effective loop over i
	  for(nrow=pij+c; pij<nrow; pij++,nij++)	// Effective loop over j
	    (*pij) = (*dii) * (*nij);			// <i|pdt|j> = <i|dmx|i><i|nmx|j>
	return pdt;
        break;
	}
      default:						// Multiply d_matrix into generic matrix
	{ 
        int c = mx->cols();				// Get pdt column dimension
	n_matrix* pdt = new n_matrix(rows_,c);		// Construct result as n_matrix
        complex *p00 = pdt->data;			// Start of pdt: p00 = <0|pdt|0>
        complex *d00 = data;				// Start of dmx: d00 = <0|dmx|0>
        complex *pij, *dii;				// This will be matrix elements
        int i,j;
	for(i=0, pij=p00, dii=d00; i<rows_; i++,dii++)	// Fill all elements <i|pdt|j>
	  for(j=0; j<c; j++,pij++)
	    (*pij) = (*dii) * (*mx)(i,j);		// <i|pdt|j> = <i|dmx|i><i|mx|j>
	return pdt;
	}
    }
  }


_matrix* d_matrix::multiply(const complex& z)
 
        // Input                dmx  : Input diagonal matrix (this)
        //                      z    : Complex number
        // Output               pdt  : New matrix which is the product z*dmx

  {
  if(z == complex1) return this;		// For z=1 just return dmx
  else if(z == complex0)			// Check if z is 0
    { 						//   If so, return a zero dmx
    d_matrix* pdt =
              new d_matrix(rows_,cols_,complex0);
    return pdt;
    }
  else						// Else must compute result
    {						// The product is another d_matrix
    d_matrix* pdt = new d_matrix(rows_,cols_);
    for(int i=0; i<size; i++)
     pdt->data[i] = z*data[i];			// <i|pdt|i> = z*<i|dmx|i>
    return pdt;
    }
  }


_matrix* d_matrix::divide(_matrix* mx)
         
        // Input                dmx  : Input diagonal matrix (this)
        //                       mx  : Second matrix                         -1
        // Output               pdt  : New matrix which is the product dmx*mx
        // Note                      : This uses the function inv of mx
 
/*                    ---                -1                    -1
          <i|pdt|j> = \   <i|dmx|k> <k|mx  |j> = <i|dmx|i><i|mx  |j>
                      /
                      ---
*/ 

  {
  if(!CheckDim(mx,false,1)) DMxfatal(23);		// Quit if mismatch
  switch(mx->stored_type())
    {
    case i_matrix_type: return this; break;		// No change if dmx/I
    case d_matrix_type:					// Divide dmx by dmx
      {
      d_matrix* pdt = new d_matrix(cols_,cols_);	// Create new diagonal matrix
      for(int i=0; i<cols_; i++) 			// Loop over all columns
        pdt->data[i] = data[i]/((d_matrix*)mx)->data[i]; // <i|pdt|i> = <i|dmx|i>/<i|mx|i>
      return pdt;
      }  
      break;
    case h_matrix_type:                                 // Divide dmx by hmx
      return multiply((h_matrix*)mx->inv());		// Return dmx*inv(hmx)
      break;
    case n_matrix_type:                                 // Divide dmx by nmx
      return multiply((n_matrix*)mx->inv());		// Return dmx*inv(nmx)
      break;
    default: break; 					// Divide dmx by generic mx
    }
  DMxerror(25,"divide",1);                              // Division not fully implemented
  DMxfatal(23);                                         // Unable to do division
  return this;
  }


         
        // Input                dmx  : Input diagonal matrix (this)
        //                      z    : Complex number
        // Output               pdt  : New matrix which is the product (1/z)*dmx
        // Note                      : This handles divisions in code such as mx/z

_matrix* d_matrix::divide(const complex& z)
  {
  if(z == complex1) return this;  		// If z=1, just return dmx
  else if(z == complex0)			// If z=0, fatal error
    {
    DMxerror(18,1);                             //    Divisions by zero
    DMxfatal(23);                               //    Cannot perform division
    }
  complex z1 = 1/z;
  d_matrix* mx = new d_matrix(rows_,cols_);
  for (int i=0; i<size; i++)
    mx->data[i] = z1*data[i];
  return mx;
  }
 
 
// ____________________________________________________________________________
// G                 CLASS D_MATRIX UNARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________
 

_matrix* d_matrix::add_two(_matrix* mx)	// Implements mx += d_matrix
 
        // Input                dmx  : Input diagonal matrix (this)
        //                      mx   : Second matrix
        // Output               sum  : Matrix mx has dmx added into it
        // Note                      : This is the implementation of mx += dmx
	// Note			     : Uses internal structure of other matrix types

  {
  int dd = rows_;				// Dimension of dmx
  if(!CheckDims(mx, 1)) DMxfatal(20);		// Quit if mismatch
  switch(mx->stored_type())
    {
      case n_matrix_type:			// Add a d_matrix into an n_matrix
	{
        complex *dii = data;			// Start of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *nii = ((n_matrix*)mx)->data;	// Start of nmx: nii = <i|nmx|i> -> <0|nmx|0>
        complex *dend = dii+dd;			// End of dmx: <dd|dmx|dd+1>
        int nrow = dd+1;			// Amount to take <i|nmx|i> -> <i+1|nmx|i+1>
	for(; dii<dend; dii++,nii+=nrow)
	  *nii += *dii;				// <i|nmx|i> += <i|dmx|i>
	return mx;
	}
	break;
      case h_matrix_type:			// Add a d_matrix into an h_matrix
	{
        complex *dii = data;			// Start of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *dend = dii+dd;			// End of dmx: <dd|dmx|dd+1>
        if(is_real())				// If dmx is real, result is Hermitian
          {
          complex *hii = ((h_matrix*)mx)->data;	// Start of hmx: hii = <i|hmx|i> -> <0|hmx|0>
          int hrow = dd;			// Amount to take <i|hmx|i> -> <i+1|hmx|i+1>
	  for(; dii<dend; dii++,hii+=hrow,hrow--)
	    *hii += *dii;			// <i|hmx|i> += <i|dmx|i>
	  return mx;
          }
        else					// If dmx is complex, result is Normal
          {
	  n_matrix* mx1 = new n_matrix(dd,dd);	// Construct a new n_matrix for sum
          complex *dii = data;			// Start of dmx: dii = <i|dmx|i> -> <0|dmx|0>
          complex *mii = ((n_matrix*)mx1)->data;// Start of mx:  mii = <i|mx|i> -> <0|mx|0>
          complex *hij = ((h_matrix*)mx)->data;	// Start of hmx: hij = <i|hmx|j> -> <0|hmx|0>
          complex *mend = mii + dd*dd;		// End of mx: <dd|mx|dd+1>
          complex *hend = hij + (dd*dd+dd)/2;	// End of data in hmx: <dd|hmx|dd+1>
          complex *mij, *mji;
          int mrow = dd + 1;			// Increment to take <i|mx|i> -> <i+1|mx|i+1>
	  for( ; hij<hend; dii++,mii+=mrow)	// Effective loop over i (hmx rows), at this
            {					// point hij = <i|hmx|i>
            *mii = *hij + *dii;			// <i|mx|i> = <i|hmx|i> + <i|dmx|i>
            hij++;				// <i|hmx|j> = <i|hmx|i> -> <i|hmx|i+1>
            mji = mii + dd;			// <j|mx|i> -> <i+1|mx|i>
            mij = mii + 1;			// <i|mx|j> -> <i|mx|i+1>
	    for(; mji<mend; mij++,mji+=dd,hij++)// Effecive loop over j (mx cols), j>i
              {
              *mji = conj(*hij);		// <j|mx|i> = <i|hmx*|j>, i<j
              *mij = *hij;			// <i|mx|j> = <i|hmx|j>, i<j
              }
            }
	  return mx1;
	  }
	}
	break;
      case d_matrix_type:			// Add d_matrix into another d_matrix
	{ 
        complex *dii = data;			// Start of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *mii = ((d_matrix*)mx)->data;	// Start of mx: mii = <i|mx|i> -> <0|mx|0>
        complex *dend = dii+dd;			// End of dmx: <dd|dmx|dd+1>
	for(; dii<dend; dii++,mii++)
	  *mii += *dii;				// <i|mx|i> += <i|dmx|i>
	return mx;
	}
	break;
      case i_matrix_type:			// Add d_matrix into an i_matrix
	{ 					// Type must change to d_matrix then
	d_matrix* mx1=new d_matrix(dd,dd);	// so constuct a new one
        complex *dii = data;			// Start of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *mii = ((d_matrix*)mx1)->data;	// Start of mx: mii = <i|mx|i> -> <0|mx|0>
        complex *dend = dii+dd;			// End of dmx: <dd|dmx|dd+1>
	for(;  dii<dend; dii++,mii++)
	  *mii = *dii + complex1;		// <i|mx|i> = <i|dmx|i> + <i|imx|i> = <i|dmx|i> + 1
	return mx1;
	}
	break;
      default:					// Add a d_matrix into a generic matrix
	{ 
	n_matrix* mx1 = new n_matrix(dd,dd);
	int pos = 0;
        int i,j;
	for(i=0; i<dd; i++)
	  for(j=0; j<dd; j++,pos++)
	    if(i == j)
	      mx1->data[pos]=data[i]+(*mx)(i,j);// <i|mx1|i> = <i|dmx|i> + <i|mx|i>
	    else
	      mx1->data[pos]=(*mx)(i,j);	// <i|mx1|j> = <i|mx|j>, j != i
	return mx1;
	}
      }
  }


_matrix* d_matrix::subtract_two(_matrix* mx)
 
        // Input                dmx  : Input diagonal matrix (this)
        //                      mx   : Second matrix
        // Output               mx   : Matrix mx has dmx subtracted from it
        // Note                      : This is the implementation of mx -= dmx
	// Note			     : Uses internal structure of other matrix types

  {
  int dd = rows_;				// Dimension of dmx
  if(!CheckDims(mx, 1)) DMxfatal(21);		// Quit if mismatch
  switch(mx->stored_type())
      {
      case n_matrix_type:			// Add a d_matrix into an n_matrix
	{
        complex *dii = data;			// Start of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *nii = ((n_matrix*)mx)->data;	// Start of nmx: nii = <i|nmx|i> -> <0|nmx|0>
        complex *dend = dii+dd;			// End of dmx: <dd|dmx|dd+1>
        int nrow = dd+1;			// Amount to take <i|nmx|i> -> <i+1|nmx|i+1>
	for(; dii<dend; dii++,nii+=nrow)
	  *nii -= *dii;				// <i|nmx|i> -= <i|dmx|i>
	return mx;
	}
	break;
      case h_matrix_type:			// Add a d_matrix into an h_matrix
	{
        complex *dii = data;			// Start of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *dend = dii + dd;		// End of dmx: <dd|dmx|dd+1>
        if(is_real())				// If dmx is real, result is Hermitian
          {
          complex *hii = ((h_matrix*)mx)->data;	// Start of hmx: hii = <i|hmx|i> -> <0|hmx|0>
          int hrow = dd;			// Amount to take <i|hmx|i> -> <i+1|hmx|i+1>
	  for(; dii<dend; dii++,hii+=hrow,hrow--)
	    *hii -= *dii;			// <i|hmx|i> -= <i|dmx|i>
	  return mx;
          }
        else					// If dmx is complex, result is Normal
          {
	  n_matrix* mx1 = new n_matrix(dd,dd);	// Construct a new n_matrix for sum
          complex *dii = data;			// Start of dmx: dii = <i|dmx|i> -> <0|dmx|0>
          complex *mii = ((n_matrix*)mx1)->data;// Start of mx:  mii = <i|mx|i> -> <0|mx|0>
          complex *hij = ((h_matrix*)mx)->data;	// Start of hmx: hij = <i|hmx|j> -> <0|hmx|0>
          complex *mend = mii + dd*dd;		// End of mx: <dd|mx|dd+1>
          complex *hend = hij + (dd*dd+dd)/2;	// End of data in hmx: <dd|hmx|dd+1>
          complex *mij, *mji;
          int mrow = dd + 1;			// Increment to take <i|mx|i> -> <i+1|mx|i+1>
	  for( ; hij<hend; dii++,mii+=mrow)	// Effective loop over i (hmx rows), at this
            {					// point hij = <i|hmx|i>
            *mii = *hij - *dii;			// <i|mx|i> = <i|hmx|i> - <i|dmx|i>
            hij++;				// <i|hmx|j> = <i|hmx|i> -> <i|hmx|i+1>
            mji = mii + dd;			// <j|mx|i> -> <i+1|mx|i>
            mij = mii + 1;			// <i|mx|j> -> <i|mx|i+1>
	    for(; mji<mend; mij++,mji+=dd,hij++)// Effecive loop over j (mx cols), j>i
              {
              *mji = conj(*hij);		// <j|mx|i> = <i|hmx*|j>, i<j
              *mij = *hij;			// <i|mx|j> = <i|hmx|j>, i<j
              }
            }
	  return mx1;
	  }
	}
	break;
      case d_matrix_type:			// Add d_matrix into another d_matrix
	{ 
        complex *dii = data;			// Start of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *mii = ((d_matrix*)mx)->data;	// Start of mx: mii = <i|mx|i> -> <0|mx|0>
        complex *dend = dii+dd;			// End of dmx: <dd|dmx|dd+1>
	for(; dii<dend; dii++,mii++)
	  *mii -= *dii;				// <i|mx|i> -= <i|dmx|i>
	return mx;
	}
	break;
      case i_matrix_type:			// Add d_matrix into an i_matrix
	{ 					// Type must change to d_matrix then
	d_matrix* mx1 = new d_matrix(dd,dd);	// so constuct a new one
        complex *dii = data;			// Start of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *mii = ((d_matrix*)mx1)->data;	// Start of mx: mii = <i|mx|i> -> <0|mx|0>
        complex *dend = dii+dd;			// End of dmx: <dd|dmx|dd+1>
	for(;  dii<dend; dii++,mii++)
	  *mii = complex1 - (*dii);		// <i|mx|i> = <i|imx|i> - <i|dmx|i> = 1 - <i|dmx|i>
	return mx1;
	}
	break;
      default:					// Add a d_matrix into a generic matrix
	{ 
	n_matrix* nmx = new n_matrix(dd,dd);
        complex *dii = data;			// Start of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *nij = ((n_matrix*)nmx)->data;	// Start of nmx: nij = <i|nmx|j> -> <0|nmx|0>
        int i,j;
	for(i=0; i<dd; i++,dii++)
	  for(j=0; j<dd; j++,nij++)
	    if(i == j)
	      *nij = (*mx)(i,j) - (*dii);	// <i|nmx|i> = <i|mx|i> - <i|dmx|i>
	    else
	      *nij = (*mx)(i,j);		// <i|nmx|j> = <i|mx|j>, j != i
	return nmx;
	}
      }
  }

_matrix* d_matrix::multiply_two(_matrix* mx)	// Implements mx *= d_matrix
 
        // Input                dmx  : Input diagonal matrix (this)
        //                      mx   : Second matrix
        // Output               mx   : Matrix mx is multiplied into dmx
        // Note                      : This is the implementation of mx *= dmx
	// Note			     : Uses internal structure of other matrix types

  {
  int dd = rows_;					// Dimension of dmx
  if(!CheckDim(mx)) 					// Quit if mismatch
    { DMxerror(5,"multiply_two",1); DMxfatal(22); }
  switch(mx->stored_type())
    {
    case n_matrix_type:			// Multiply d_matrix by n_matrix
      {
      complex *nij = ((n_matrix*)mx)->data;	// Start of nmx: mij = <i|nmx|j> -> <0|nmx|0>
      complex *d00 = data;			// Start of dmx: d00 = <0|dmx|0>
      complex *nend = nij + mx->rows()*dd;	// End of nmx: nend = <mx->rows()|nmx|dd+1>
      complex *dend = d00 + dd;		// End of dmx: dend = <dd|dmx|dd+1>
      complex *djj;
      while(nij<nend)
	  for(djj=d00; djj<dend; nij++,djj++)
	    (*nij) *= (*djj);			// <i|nmx|j> = <i|nmx|j> * <j|dmx|j>
      return mx;
      }
      break;
    case h_matrix_type:			// Multiply d_matrix by h_matrix
	{					// This results in a n_matrix (leider)
	n_matrix* nmx = new n_matrix(dd,dd);	// Construct a new n_matrix for result
        complex *nij = ((n_matrix*)nmx)->data;	// Start of nmx:  nij = <i|nmx|j> -> <0|nmx|0>
        complex *h0i = ((h_matrix*)mx)->data;	// Start of hmx: <0|hmx|i> -> <0|hmx|0>
        complex *d00 = data;			// Start of dmx: d00 = <0|dmx|0>
        complex *dx = d00 + dd;			// End of dmx: dx = <dd|dmx|dd+1>
        complex *djj, *dii;			// Elements <j|dmx|j>, <i|dmx|i>
        complex *hij, *hji;			// Elements <i|hmx|j>, <j|hmx|i>
        int hrow;				// Length of hmx row j 
	for(dii=d00; dii<dx; h0i++,dii++)	// Effective loop over index i (vid dii)
          {
	  for(hrow=dd,djj=d00,hji=h0i; djj<dii;	// Effective loop over index j (via djj), j<i
	          nij++,djj++,hrow--,hji+=hrow)	// starting at j=0, ending at j=i
	    *nij = conj(*hji) * (*djj);		// <i|nmx|j> = <j|hmx*|i>*<j|dmx|j>, j<i
	  for(hij=hji;djj<dx;nij++,djj++,hij++)	// Continue effective loop over j (via djj), j>=i
	    *nij = (*hij) * (*djj);		// <i|nmx|j> = <i|hmx|j>*<j|dmx|j>, j>=i
          }
	return nmx;
	}
        break;
    case d_matrix_type:			// Multiply d_matrix by d_matrix
	{ 
        complex *mii = ((d_matrix*)mx)->data;	// Start of mx:  mii = <i|mx|i> -> <0|mx|0>
        complex *dii = data;			// Start of dmx: d00 = <0|dmx|0>
        complex *dend = dii + dd;		// End of dmx: dend = <dd|dmx|dd+1>
	for(; dii<dend; dii++,mii++)
	  (*mii) *= (*dii);			// <i|mx|i> = <i|mx|i> * <i|dmx|i>
	return mx;
	}
	break;
    case i_matrix_type:			// Multiply d_matrix by d_matrix
        return this; 				// imx * dmx = dmx so return dmx
	break;
    default:					// Multiply generic matrix into d_matrix
	{ 
        int rd = mx->rows();			// Row dimension of mx 
	n_matrix* nmx = new n_matrix(rd,dd);	// Construct a new n_matrix for result
        complex *nij = ((n_matrix*)nmx)->data;	// Start of nmx:  nij = <i|nmx|j> -> <0|nmx|0>
        complex *d00 = data;			// Start of dmx: d00 = <0|dmx|0>
        complex *djj;				// Element <j|dmx|j>
        int i,j;
	for(i=0; i<rd; i++)			// Loop over the rows in nmx
	  for(j=0,djj=d00; j<dd; j++,nij++,djj++)
	    *nij = (*mx)(i,j) * (*djj);		// <i|nmx|j> = <i|mx|j> <j|dmx|j>
	return nmx;
      }
    }
  }


_matrix* d_matrix::multiply_two(const complex &z)	// Implements d_matrix *= z
 
        // Input                dmx  : Input diagonal matrix (this)
        //                      z    : Complex number
        // Output               mx   : Matrix dmx is multiplied by z
        // Note                      : This is the implementation of dmx *= z

  {
  if(z != complex1)
    {
    complex *dii = data;			// Start of dmx: d00 = <0|dmx|0>
    complex *dend = dii + cols_;		// End of dmx: dend = <dd|dmx|dd+1>
    for(; dii<dend; dii++)
      *dii *= z;
    }
  return this;
  }


_matrix* d_matrix::divide_two(_matrix* mx)		// Implements mx /= d_matrix
 
        // Input                dmx  : Input diagonal matrix (this)
        //                      mx   : Second matrix
        // Output               mx   : Matrix mx is multiplied into inverse of dmx
        // Note                      : This is the implementation of mx /= dmx
	// Note			     : Uses internal structure of other matrix types

// sosi - copy multiply_two and then use /=
  {
  if (rows()!=mx->cols())				// Perform a size check
    {
    DMxerror(17, 1);					//   Dimension problems
    DMxerror(6, " Matrix Division", 1);			//   Matrix division problems
    DMxfatal(3, "divide_two");
    return mx;
    }
  else
    switch (mx->stored_type())
      {
      case d_matrix_type:				// Divide two d_matrices
	{ 
	for(int j=0; j<rows(); j++)
	  ((d_matrix*)mx)->data[j] /= data[j];
	return mx;
	}
	break;
      case n_matrix_type:				// Divide a n_matrix by a d_matrix
	{
	for(int i=0, pos=0; i<mx->rows(); i++)
	  for(int j=0; j<rows(); j++,pos++)
	    ((n_matrix*)mx)->data[pos] /= data[j];
	return mx;
	break;
	}
      case i_matrix_type:				// Divide a i_matrix by a d_matrix
        {
	d_matrix* mx1 = new d_matrix(rows(),mx->cols());// Produce a new d_matrix
	for(int j=0; j<rows(); j++)
	  ((d_matrix*)mx1)->data[j] = 1.0/data[j];
	return mx1;
	break;
        }
      default:						// Divide a generic matrix by a d_matrix
	{ 
	n_matrix* mx1 = new n_matrix(rows(),mx->cols());// Produce a normal matrix
	int pos = 0;					// Position in mx1->data
	for(int i=0;  i<mx1->rows(); i++)
	  for(int j=0; j<mx1->cols(); j++,pos++)
	    mx1->data[pos] = (*mx)(i,j)/data[j];
// sosi - use of access function here
	return mx1;
	}
      }
  }


_matrix* d_matrix::divide_two( const complex &z)	// Implements d_matrix /= z
 
        // Input                dmx  : Input diagonal matrix (this)
        //                      z    : Complex number
        // Output               mx   : Matrix dmx is multiplied by (1/z)
        // Note                      : This is the implementation of dmx /= z

  {
  if(z != complex1)
    {
    complex z1 = 1/z;
    complex *dii = data;			// Start of dmx: d00 = <0|dmx|0>
    complex *dend = dii + cols_;		// End of dmx: dend = <dd|dmx|dd+1>
    for(; dii<dend; dii++)
      *dii *= z1;
    }
  return this;
  }
 
// ____________________________________________________________________________
// H                 CLASS D_MATRIX SIMPLE UNARY FUNCTIONS
// ____________________________________________________________________________

/* These functions perform simple simple mathematical operations on a Hermitian
   matrix. Note that for Hermitian arrays the adjoint function does nothing.
 
    Function          Result             Function           Result
   ---------- -----------------------    --------  ------------------------
   negate     <i|dmx|j> -> -<i|dmx|j>       RE      <i|dmx|j>->Re(<i|dmx|j>)
   conjugate  <i|dmx|j> -> <i|dmx*|j>       IM      <i|dmx|j>->Im(<i|dmx|j>)
   transpose  <i|dmx|j> -> <j|dmx|i>      adjoint   <i|dmx|j>-><j|dmx*|i>   */

_matrix* d_matrix::negate()
  {
  d_matrix* mx = new d_matrix(cols_, rows_);	// Construct a new d_matrix
  for(int i=0; i<rows_; i++)			// Loop over all the non-zero elements
    mx->data[i] = -data[i];			// <i|mx|i> = -<i|dmx|i>
  return mx;
  }

 _matrix* d_matrix::RE()
   {
   if(is_real()) return this;			// If real, no change
   d_matrix* mx = new d_matrix(cols_, rows_);	// If complex, construct new d_matrix
   for(int i=0; i<rows_; i++)			// Loop over all non_zero elements
     mx->data[i] = Re(data[i]);			// and just use the real parts
   return mx;
   }
  
 _matrix* d_matrix::IM()
   {
   if(is_imaginary()) return this;		// If imaginary, no change
   d_matrix* mx = new d_matrix(cols_, rows_);	// If complex, construct new d_matrix
   for(int i=0; i<rows_; i++)			// Loop over all non_zero elements
     mx->data[i] = Im(data[i]);			// and just use the real parts
   return mx;
   }

_matrix* d_matrix::conjugate()
  {
  if(is_real()) return this;			// If it is real, no change
  d_matrix* mx = new d_matrix(cols_, rows_);	// If complex, construct new d_matrix
  for(int i=0; i<rows_; i++)			// Loop over all non_zero elements
    mx->data[i] = conj(data[i]);		// <i|mx|i> = <i|dmx*|i>
  return mx;
  }

_matrix* d_matrix::transpose() { return this; }	// No change by transposition

_matrix* d_matrix::adjoint()
  {
  if(is_real())	return this;			// If it is real, no change
  d_matrix* mx = new d_matrix(cols_, rows_);	// If complex, construct new d_matrix
  for(int i=0; i<rows_; i++)			// Loop over all non_zero elements
    mx->data[i] = conj(data[i]);		// <i|mx|i> = <i|dmx*|i>
  return mx;
  }

_matrix* d_matrix::mxexp()
  {
  d_matrix* mx = new d_matrix(cols_, rows_);	// If complex, construct new d_matrix
  for(int i=0; i<rows_; i++)			// Loop over all non_zero elements
    mx->data[i] = exp(data[i]);			// <i|mx|i> = exp(<i|dmx*|i>)
  return mx;
  }

complex d_matrix::trace()
  {
  complex z(0);
  for(int i=0; i<rows_; i++) z += data[i];
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

_matrix* d_matrix::swaprows(int i, int j)
  {
  if(data[i]==data[j] && is_real())			// If elements equal
    {							// & real, Hermitian
    h_matrix* dsr = new h_matrix(rows_,rows_,complex0); // Make new h_matrix
    for(int k=0; k<rows_; k++) dsr->put(data[k],k,k);	// Make it dmx
    dsr->put(complex0,  i,i);                           // Do the row swap now,
    dsr->put(complex0,  j,j);                           // there are only 4
    dsr->put(data[i],   i,j);                           // elements to change
    return dsr;
    }
  n_matrix* dsr = new n_matrix(rows_,rows_,complex0);   // Make new h_matrix
  for(int k=0; k<rows_; k++) dsr->put(data[k],k,k);	// Make it dmx
  dsr->put(complex0,  i,i);				// Do the row swap now,
  dsr->put(complex0,  j,j);				// there are only 4
  dsr->put(data[i],   j,i);
  dsr->put(data[j],   i,j);
  return dsr;
  }

_matrix* d_matrix::swapcols(int i, int j)
  {
  if(data[i]==data[j] && is_real())			// If elements equal
    {							// & real, Hermitian
    h_matrix* dsc = new h_matrix(rows_,rows_,complex0); // Make new h_matrix
    for(int k=0; k<rows_; k++) dsc->put(data[k],k,k);	// Make it dmx
    dsc->put(complex0,  i,i);                           // Do the row swap now,
    dsc->put(complex0,  j,j);                           // there are only 4
    dsc->put(data[i],   i,j);                           // elements to change
    return dsc;
    }
  n_matrix* dsc = new n_matrix(rows_,rows_,complex0);   // Make new h_matrix
  for(int k=0; k<rows_; k++) dsc->put(data[k],k,k);	// Make it dmx
  dsc->put(complex0,  i,i);				// Do the row swap now,
  dsc->put(complex0,  j,j);				// there are only 4
  dsc->put(data[i],   i,j);
  dsc->put(data[j],   j,i);
  return dsc;
  }

_matrix* d_matrix::permute( int i, int j)
  {
  d_matrix* dp = new d_matrix(*this);			// Permutation  on
  dp->data[i] = data[j];				// dmx produces dmx
  dp->data[j] = data[i];
  return dp;
  }

double  d_matrix::maxRe() const
  { 
  double maxval=-HUGE;
  for(int i=0; i<rows_; i++) maxval = gmax(Re(data[i]), maxval);
  return (rows_)?maxval:0;  
  }

double  d_matrix::maxIm() const
  { 
  double maxval=-HUGE;
  for(int i=0; i<rows_; i++) maxval = gmax(Im(data[i]), maxval);
  return (rows_)?maxval:0;  
  }

complex d_matrix::maxZ()  const
  { 
  double maxval=-HUGE;
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

double  d_matrix::minRe() const
  { 
  double minval=HUGE;
  for(int i=0; i<rows_; i++) minval = gmin(Re(data[i]), minval);
  return (rows_)?minval:0;  
  }

double  d_matrix::minIm() const
  { 
  double minval=HUGE;
  for(int i=0; i<rows_; i++) minval = gmin(Im(data[i]), minval);
  return (rows_)?minval:0;  
  }

complex d_matrix::minZ()  const
  { 
  double minval=HUGE;
  complex minz;
  for(int i=0; i<rows_; i++)
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
// I                CLASS N_MATRIX SIMPLE BINARY FUNCTIONS
// ____________________________________________________________________________

 
 
        // Input            dmx : A d_matrix (this)
        //                   mx : A second matrix
        // Output            z  : A complex value which is the trace of
        //                        the product of two input arrays
        //                         z = Tr{dmx*mx} = Tr{mx*dmx}
        // Note                 : This is faster than taking the matrix
        //                        product and then using unary trace!
        // Note                 : Uses internal structure of other matrix types

//       ---                 ---
//   z = \   <i|dmx*mx|i>  = \   <i|dmx|i><i|mx|i>
//       /                   / 
//       ---                 ---
//        i                   i

complex d_matrix::trace(_matrix* mx)
  {
  if((rows_!=mx->cols())||cols_!=mx->rows())		// Insure dimensions match
    {
    DMxerror(17, 1);					//   Dimension problems
    DMxerror(6, " Trace on Rectangular Matrix", 1);	//   Matrix division problems
    DMxfatal(3, "trace");				//   Fail in trace function
    return complex0;
    }
  else
    {
    complex z(0);				// Initialize the result at zero
    switch(mx->stored_type())
      {
      case n_matrix_type:                       // Tr{d_matrix * n_matrix}
        {                  
        complex *d00 = data;			// Start of dmx: d00 = <0|dmx|0>
        complex *n00 = ((n_matrix*)mx)->data;	// Start of nmx: n00 = <0|nmx|0>
        complex *dz = d00 + cols_;		// End of dmx: dz = <rows_|dmx|cols_+1>
        complex *nii, *dii;			// These are the matrix elements
        for(nii=n00,dii=d00; dii<dz;
                           dii++,nii+=cols_+1)	
          z += (*dii) * (*nii);               // Tr{dmx*dmx} += <i|dmx|i><i|nmx|i>
        return z;
        }
        break;
      case h_matrix_type:                       // Tr{d_matrix * h_matrix}
        {
        complex *d00 = data;                    // Start of dmx: d00 = <0|dmx|0>
        complex *h00 = ((h_matrix*)mx)->data;   // Start of hmx: h00 = <0|hmx|0>
        complex *dz = d00 + rows_;		// End of dmx: dz = <rows_|dmx|rows+1>
        complex *dii, *hii;			// These are the matrix elements
        int hrow = cols_;			// hmx row increment: <i|hmx|i>+hrow
        for(hii=h00,dii=d00; dii<dz; 		//                        -> <i+1|hmx|i+1>
                        dii++,hii+=hrow,hrow--)
          z += (*dii) * (*hii);			// Tr{dmx*hmx} += <i|dmx|i><i|hmx|i>
        return z;
        }  
      case d_matrix_type:                       // Tr{d_matrix * d_matrix}
        {
        complex *d00 = data;			// Start of dmx: d00 = <0|dmx|0>
        complex *m00 = ((d_matrix*)mx)->data;	// Start of mx: m00 = <0|mx|0>
        complex *dz = d00 + rows_;		// End of dmx: dz = <rows_|dmx|rows+1>
        complex *dii, *mii;			// These are the matrix elements
        for(mii=m00,dii=d00; dii<dz; mii++,dii++)
          z += (*dii) * (*mii);                 // Tr{dmx*mx} += <i|dmx|i><i|mx|i>
        return z;
        }
        break;
      case i_matrix_type:                       // Tr{d_matrix * i_matrix}
        {
        for(int i=0; i<rows_; i++)		// Sum over the diagonal elements
          z += data[i];				// Tr{dmx*imx} += <i|dmx|i><i|imx|i>
        return z;				//             += <i|dmx|i>
        }
        break;
    default:					// Tr{d_matrix * unknown matrix}
        {
        complex *d00 = data;                    // Start of dmx: d00 = <0|dmx|0>
        complex *dii;                           // This is matrix element <i|dmx|i>
        int i;
        for(i=0, dii=d00; i<rows_; i++,dii++)
          z += (*dii) * (*mx)(i,i);           // Tr{dmx*mx} += <i|dmx|i><i|mx|i>
        return z;
        }
      }
    }
  }


_matrix* d_matrix::adjoint_times(_matrix* mx)

        // Input            dmx : A d_matrix (this)
        //                   mx : A second matrix
        // Output           mx1 : A matrix which is product of the
        //                        adjoint of d_matrix and mx
        //                                    T *
        //                         mx1 = [(dmx ) ] * mx
        // Note                 : This is faster than taking the adjoint
        //                        of nmx and then the product!  Use if the
        //                        adjoint of nmx is not needed.
        // Note                 : The algorithm used is identical to
        //                        multiply, but the used indices of nmx
        //                        are switched and conj(*a,*b) is used
        //                        where conj(z1,z2) = conj(z1)*z2!
        // Note                 : This uses knowledge of the internal matrix
        //                        structure of other classes
 
//                        ---
//            <i|pdt|j> = \   <k|dmx*|i> <k|mx|j> = <i|dmx*|i><i|mx|j> 
//                        /
//                        ---

  {
  if(rows_ != mx->rows())				// Insure matrix dimensions proper
    {
    DMxerror(17, 1);					//   Dimension problems
    DMxerror(6, " adjoint(mx1)*mx2", 1);		//   Matrix multiply problems
    DMxfatal(3, "adjoint_times");			//   Fail in adjoint_times function
    return mx;
    }
  else
    {
    switch(mx->stored_type())
      {
      case n_matrix_type:				// Adjoint of dmx * normal array
	{ 
        int c = mx->cols();				// Get pdt column dimension
	n_matrix* pdt = new n_matrix(rows_,c);		// Construct result as n_matrix
        complex *p00 = pdt->data;			// Start of pdt: p00 = <0|pdt|0>
        complex *d00 = data;				// Start of dmx: d00 = <0|dmx|0>
        complex *n00 = ((n_matrix*)mx)->data;		// Start of nmx: h00 = <0|nmx|0>
        complex *pend = p00 + rows_*c;			// Start of nmx: h00 = <0|nmx|0>
        complex *pij, *dii, *nij;			// These will be matrix elements
        complex *nrow;					// Dummy index for end of nmx row
	for(pij=p00,nij=n00,dii=d00; pij<pend; dii++)	// Effective loop over i
	  for(nrow=pij+c; pij<nrow; pij++,nij++)	// Effective loop over j
	    (*pij) = conj(*dii, *nij);			// <i|pdt|j> = <i|dmx*|i><i|nmx|j>
	return pdt;
        break;
	}
      case h_matrix_type:				// Multiply d_matrix into h_matrix
	{ 
	n_matrix* pdt = new n_matrix(rows_,cols_);	// Construct result as n_matrix
        complex *pij = pdt->data;			// Element of pdt: pij = <i|pdt|j> -> <0|pdt|0>
        complex *dii = data;				// Element of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *h00 = ((h_matrix*)mx)->data;		// Start of hmx: h00 = <0|hmx|0>
        int i=0,j;
        complex *hji, *hij;
	for(i=0; i<rows_; i++,dii++)			// Loop over the rows of the matrices
          {
          hji = h00+i;					// Set <j|hmx|i> to <0|hmx|i>
	  for(j=0; j<i; j++,pij++,hji+=cols_-j)		// Begin loop over j, i>j
	    (*pij) = conj(*dii) * conj(*hji);		// <i|pdt|j> = <i|dmx|i><j|hmx*|i>, i>j
	  for(hij=hji; j<cols_; j++,pij++,hij++)	// Continue loop over j, i<=j
	    (*pij) = conj(*dii, *hij);			// <i|pdt|j> = <i|dmx|i><i|hmx|j>, i<=j
          }
	return pdt;
        break;
	}
      case d_matrix_type:				// Adjoint of dmx * diagonal array
	{ 
	d_matrix* pdt = new d_matrix(rows_,rows_);	// New diagonal matrix for product
        complex *pii = pdt->data;			// Element of pdt: pii = <i|pdt|i> -> <0|pdt|0>
        complex *dii = data;				// Element of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *dend = dii + rows_;			// End of dmx: dend = <rows_|dmx|cols_+1>
        complex *mii = ((d_matrix*)mx)->data;		// Element of mx: mii = <i|mx|i> -> <0|mx|0>
	for( ; dii<dend; dii++,mii++,pii++)
          (*pii) = conj(*dii, *mii);	 		// <i|pdt|i> = <i|dmx*|i><i|mx|i>
	return pdt;
	}
	break;
      case i_matrix_type:				// Adjoint of dmx * identity array
	return adjoint();				// Just return the adjoint of dmx
	break;
      default:						// Adjoint of dmx * generic array
	{ 
        int c = mx->cols();				// Column dimension of product
	n_matrix* pdt = new n_matrix(rows_,c);
        complex *pij = pdt->data;			// Element of pdt: pij = <i|pdt|j> -> <0|pdt|0>
        complex *dii = data;				// Element of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        int i,j;
	for(i=0; i<rows_; dii++,i++)
	  for(j=0; j<c; j++,pij++)
	    (*pij) = conj(*dii, (*mx)(i,j)); 		// <i|pos|j> = <i|dmx*|i><i|mx|j>
	return pdt;
	}
      }
    }
  }


_matrix* d_matrix::times_adjoint(_matrix* mx)
 
        // Input            dmx : A d_matrix (this)
        //                   mx : A second matrix
        // Output           mx1 : A matrix which is product of the
        //                        n_matrix and the adjoint of mx
        //                                        T *
        //                         mx1 = nmx* [(mx ) ]
        // Note                 : This is faster than taking the adjoint
        //                        of mx and then the product!  Use if the
        //                        adjoint of mx is not needed.
        // Note                 : The algorithm used here is identical to
        //                        multiply except the indices of mx are
        //                        switched and conj(*a,*b) is used!
        // Note                 : This uses knowledge of the internal matrix
        //                        structure of other classes
 
//                        ---
//            <i|pdt|j> = \   <i|dmx|k> <j|mx*|k> = <i|dmx|i><j|mx|i> 
//                        /
//                        ---

  {
// sosi -switched intial size check as cols to cols vs. cols to rows
//       corrected for rectangular arrays and added h_matrix
  int r = rows_;
  if(rows_ != mx->cols())				// Insure matrix dimensions match
    {
    DMxerror(17, 1);					//   Dimension problems
    DMxerror(6, " mx1* adjoint(mx2)", 1);		//   Matrix multiply problems
    DMxfatal(3, "times_adjoint");			//   Fail in times_adjoint function
    return mx;
    }
  else
    switch (mx->stored_type())
      {
      case n_matrix_type:				// Multiply d_matrix into n_matrix adjoint
	{ 
        int c = mx->rows();				// Column dimension of product
	n_matrix* pdt = new n_matrix(r,c);
	int pos = 0;					// Position in mx1->data
	for(int i=0; i<r; i++)
	  for(int j=0; j<c; j++,pos++)
	    pdt->data[pos] = 				// <i|pdt|j> = <i|dmx|i>*<j|mx*|i> 
            data[i]*conj((*(n_matrix*)mx)(j,i));
// sosi - use of access function here
	return pdt;
	}
      case h_matrix_type:				// Multiply d_matrix into h_matrix adjoint
	return multiply((h_matrix*)mx);			// Since hmx is self adjoint, return dmx*hmx;
	break;
      case d_matrix_type:				// Multiply d_matrix into d_matrix adjoint
	{ 
	d_matrix* pdt = new d_matrix(rows_,rows_);	// New diagonal matrix for product
        complex *pii = pdt->data;			// Element of pdt: pii = <i|pdt|i> -> <0|pdt|0>
        complex *dii = data;				// Element of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *dend = dii + rows_;			// End of dmx: dend = <rows_|dmx|cols_+1>
        complex *mii = ((d_matrix*)mx)->data;		// Element of mx: mii = <i|mx|i> -> <0|mx|0>
	for( ; dii<dend; dii++,mii++,pii++)
          (*pii) = conj(*mii, *dii);	 		// <i|pdt|i> = <i|dmx|i><i|mx*|i>
	return pdt;
	}
	break;
      case i_matrix_type:				// Multiply d_matrix into i_matrix adjoint 
	return this;					// Just return dmx
	break;
      default:						// Multiply d_matrix into generic mx adjoint
	{ 
        int c = mx->cols();				// Column dimension of product
	n_matrix* pdt = new n_matrix(rows_,c);
        complex *pij = pdt->data;			// Element of pdt: pij = <i|pdt|j> -> <0|pdt|0>
        complex *dii = data;				// Element of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        int i,j;
	for(i=0; i<rows_; dii++,i++)
	  for(j=0; j<c; j++,pij++)
	    (*pij) = conj((*mx)(j,i), *dii); 		// <i|pos|j> = <i|dmx|i><j|mx*|i>
	return pdt;
	}
      }
  }

// __________________________________________________________________________________
// J                   CLASS D_MATRIX COMPLEX UNARY FUNCTIONS
// __________________________________________________________________________________

// Note:  Unlike class n_matrix, functions for inverse and diagonalization
//        are relatively simple and defined directly in the d_matrix.cc
// Note:  The alogorithm for FFT is not found here.  The diagonal matrix is
//        converted to a normal matrix and it's transformation occurs in
//        class n_matrix.

        // Input            dmx : A d_matrix (this)
        // Output          dinv : Returns the inverse of dmx as a
	//			  diagonal matrix

        //                               -1
        //                        I = dmx   * dmx = dinv * dmx

_matrix* d_matrix::inv()  
  {
  d_matrix* mx = new d_matrix(cols_,cols_);	// Construct a d_matrix for inverse
  complex *mii = mx->data;			// Element of mx: mii = <i|mx|i> -> <0|mx|0>
  complex *dii = data;				// Element of dmx: dii = <i|dmx|i> -> <0|dmx|0>
  complex *dend = dii + cols_;			// End of dmx & mx
  for(; dii<dend; mii++, dii++)
    (*mii) = 1/(*dii);
  return mx;
  }


_matrix* d_matrix::LU(int *indx)

        // Input            dmx : An d_matrix (this)
	//		    indx: Row permutation array
	// Output           dLU : LU decomposition of a row permutation
	//			  of the input matrix dmx, dmx', the row
	//			  exchanges recorded in array indx
	//			  dLU(low.tri.) * dLU(upp.tri.) = dmx'
	//			         L      *       U       = dmx'
  { 
  indx[0] = -1;				// Set indx to show no row permutations
  return this;
  }



        // Input            D   : A d_matrix (this) of results
	//		    LU  : An LU decomposition of matrix A
	//			  (or row permutation of matrix A)
	//		    indx: Row permutation array
	// Output           X   : Solution to the equation: A*X = LU*X = D
	//			  where LU is the LU decomposition of A
	//			  (with possible row permutation of A)
	// Note		        : Matrix LU must be square
	// Note		        : Uses internal structure of other matrix types

 _matrix* d_matrix::LUinv(int *indx, _matrix* LU)
  {
  if(!LU->is_square())
    {
    DMxerror(17, 1);					//   Dimension problems
    DMxerror(6, " LU Inversion On Rectangular Mx", 1);	//   LU inversion multiply problems
    DMxfatal(3, "LUinv");				//   Fail in times_adjoint function
    }
  if(rows_ != LU->rows())				// Insure LU & D dimensions match
    {
    DMxerror(17, 1);					//   Dimension problems
    DMxerror(6, " LU Inversion LU|X>=|D> Mismatch", 1);	//   LU inversion multiply problems
    DMxfatal(3, "LUinv");				//   Fail in times_adjoint function
    }
  switch(LU->stored_type())
    {
    case i_matrix_type:				// LU is an identity matrix
      if(indx[0] < 0)				// If no row permutations, LU = A = D
        return this;				// so then A*X=D -> I*X=D; X=D
      break;
    case d_matrix_type:				// LU is a diagonal matrix
      if(indx[0] < 0)				// If no row permutations, LU = A = dmx 
        return multiply(LU->inv());		// so then A*X=D -> D*X=D; X=inv(D)*D = D*inv(D);
      break;
    case h_matrix_type:				// LU is a Hermitian matrix
    case n_matrix_type:				// LU is a normal matrix
    default:
      break;
    }
  n_matrix* X = new n_matrix(rows_,cols_);	// Constuct a new n_matrix
  convert(X);					// Convert D to X as n_matrix
  _matrix* luinv = X->LUinv(indx,LU);		// Get the inverse
  delete X;					// Don't need this anymore
  return luinv;					// Use n_matrix algorithm
  }


// ____________________________________________________________________________
// N                CLASS N_MATRIX DIAGONALIZATION FUNCTIONS
// ____________________________________________________________________________

/*                             Diagonalization Routines

  This Don't Do Much Obviously. Returned Diagonal & Eigenvectors Arrays
  Must Have Their Referencing Set External To This Class.

           Input        dmx     : A diagonal matrix (this)  
                        mxd     : Pointer to mx to become diagonal mx
                        mxev    : Pointer to mx to become eigenvector mx
           Output       void    : The matrices mxd and mxev are set to
                                  the diagonal matrix of imx eigenvalues and
                                  the matrix of imx eigenvectors respectively
           Note                 : For d_matrix, both mx & mx1 are d_matrix
           Note                 : The reference count to dmx must be twice
                                  incremented external by the call origin.   */

std::vector<int> d_matrix::BlockDiag(_matrix* (&BF), std::vector<int> &U)
   {
   BF = this;
   int nr = rows_;				// Matrix dimension
   int i;
   for(i=0; i<nr; i++) { U.push_back(i); }	// Set both as unpermuted
   return std::vector<int>(1, cols_);
   }

void d_matrix::HermTriDiag(_matrix* (&HTD), _matrix* (&U))
  {
  HTD = this;
  U = new i_matrix(cols_, cols_);
  if(!is_real()) std::cout << "\n\tDiagonal Matrix: Cannot Form Hermitian TriDiagonal!";
  }

void d_matrix::SymTriDiag(_matrix* (&STD), _matrix* (&U))
  { STD = this; U = new i_matrix(cols_, cols_); }

void d_matrix::SymDiag(_matrix* (&D), _matrix* (&U))
  { D = this; U = new i_matrix(cols_, cols_); }

void d_matrix::diag(_matrix* (&D), _matrix* (&U))
  { D = this; U = new i_matrix(cols_, cols_); }
 
        // Input            dmx : A d_matrix (this)
        // Output             z : Value of the determinant; det{dmx}
 
//                            _____
//                      det =  | |  <i|dmx|i>
//                             | |

complex d_matrix::det()
  {
  complex z(1);					// Start with det = 1;
  complex *dii = data;				// dii = <i|dmx|i> -> <0|dmx|0>
  complex *dend = dii + cols_;			// End of dmx & mx
  for(; dii<dend; dii++) z *= (*dii);
  return z;
  }


int d_matrix::rank()
 
        // Input            dmx : A d_matrix (this)
        // Output             i : The value of the rank of dmx
	// Note			: For a diagonal matrix, the rank
	//			  is equivalent to the number of non-zero 
	//			  rows (or non-zero columns)
 
  {
  int r=0;					// Each non-zero column (or row)
  for(int i=0; i<rows_; i++)			// is linearly independent and adds
    if(data[i] != complex0) r++; 		// 1 to the rank
  if(r==0 && rows_)				// For a zero matrix of non-zero
    r++;					// dimension, set the rank to 1
  return r;
  }
 
 
// ____________________________________________________________________________
// K                 CLASS D_MATRIX COMPLEX BINARY FUNCTIONS
// ____________________________________________________________________________

// ***************************** tensor product *******************************
 
_matrix* d_matrix::tensor_product( _matrix* mx )

        // Input            dmx : A d_matrix (this)
        //                   mx : A second matrix
        // Output           pdt : A matrix which is the tensor product
        //                        of the two input matrices

        //                           pdt        =   dmx (x) mx

        //                       (m*n x m*o)       (mxm)   (nxo)

        //                    <i*n+k|pdt|j*o+l> = <i|dmx|j><k|mx|l>

  {
  switch (mx->stored_type())
    {
    case d_matrix_type:				// Tensor product dmx (x) dmx
      { 
      int dr = mx->rows();			// Rows of mx
      int dim = rows_*dr;			// Dimension of product
      d_matrix* pdt = new d_matrix(dim, dim);	// Tensor product is diagonal
      int pos=0;                                // Position in pdt: pdt[pos] = <i*dr+k|pdt|i*dr+k>
      for(int i=0; i<rows_; i++)		// Loop over the rows of dmx
        for(int k=0; k<dr; k++,pos++)		// Loop over the rows of mx
          pdt->data[pos] = 			// pdt[pos] = <i*dr+k|pdt|j*dr+l> = <i|dmx|j><k|mx|l>
              data[i]*((d_matrix*)mx)->data[k];	//          = del   del   <i|dmx|i><k|mx|k>
      return pdt;                               //               i,j   k,l
      }
      break;
    case i_matrix_type:				// Tensor product dmx (x) imx
      { 
      int ir = mx->rows();			// Rows of imx
      int dim = rows_*ir;			// Dimension of product
      d_matrix* pdt = new d_matrix(dim, dim);	// Tensor product is diagonal
      int pos=0;				// Position in pdt matrix
      for(int i=0; i<rows_; i++)		// Loop over the rows of dmx
        for(int k=0; k<ir; k++,pos++)		// Loop over the rows of imx
          pdt->data[pos] = data[i]; 		// pdt[pos] = <i*ir+k|pdt|j*ir+l> = <i|dmx|j><k|imx|l>
      return pdt; 				//          = del   del   <i|dmx|i>
      } 					//               i,j   k,l
      break;
    case n_matrix_type:				// Tensor Product dmx (x) nmx
      {
      int nc = mx->cols();			// Columns in nmx
      int nr = mx->rows();			// Rows in nmx
      int cm = cols_*nc;			// Columns in tensor product
      int rm = rows_*nr;			// Rows in tensor product
      n_matrix* pdt =				// Construct a new n_matrix for product
                   new n_matrix(rm,cm,complex0);
      complex *p0000 = pdt->data;		// Start of data in pdt: <00|pdt|00>
      complex *d00 = data;			// Start of data in dmx: <0|dmx|0>
      complex *n00 = ((n_matrix*)mx)->data;	// Start of data in nmx: <0|nmx|0>
      complex *dend = d00 + cols_;		// End of data in dmx: <rows_|dmx|rows_+1>
      complex *nend = n00 + nr*nc;		// End of data in nmx: <nr|nmx|nc+1>
      complex *nkp10 = NULL;			// Used to index <k+1|nmx|0>
      complex *pikil, *dii, *nkl;		// These are matrix elements
      complex *pi0i0;			
      for(pi0i0=p0000,dii=d00; dii<dend;	// Effective loop over i (via dii)
                        dii++,pi0i0+=cm*nr+nc)	// where pi0i0 = <i*nr|pdt|i*nc>
	for(nkl=n00,pikil=pi0i0; nkl<nend;	// Effective loop over k (via nkl)
                                 pikil+=cm-nc)	// where <ik|pdt|i0> -> <i(k+1)|pdt|i0> via pikil
	  for(nkp10=nkl+nc; nkl<nkp10;		// Effective loop over l (via nkl)
                                nkl++,pikil++)
	    *pikil = (*dii) * (*nkl); 		// <ik|pdt|jl> = <i*nr+k|pdt|j*nc+l> = <i|dmx|j><k|nmx|l>
      return pdt; 				//             = del   <i|dmx|i><k|nmx|l> = del   <ik|pdt|il>
      } 					//                  i,j                        i,j
      break;
    case h_matrix_type:				// Tensor product dmx (x) hmx
      { 
      int rh = mx->rows();			// Rows of Hermitian matrix
      int dim = rows_*rh;			// Dimension of tensor product
      if(is_real())				// See if dmx is real, if so
        {					// the tensor product is Hermitian
        h_matrix* pdt =				// Construct a new normal array
                new h_matrix(dim,dim,complex0);	// which is initialized to zero
        complex *p0000 = pdt->data;		// Start of data in pdt: <0,0|pdt|0,0>
        complex *h00 = ((h_matrix*)mx)->data;	// Start of data in hmx: <0|hmx|0>
        complex *d00 = data;			// Start of data in dmx: <0|dmx|0>
        complex *dend = d00 + rows_;		// End of data in dmx: <rows_|dmx|rows_+1>
        complex *hend = h00 + (rh*rh+rh)/2;	// End of data in hmx: <rh|hmx|rh+1>
        complex *pikil,*pikik,*hkl,*hkk,*dii;	// These are the matrix elements
        int prow = dim;				// Length of first row in pdt
        int hrow; 
        for(dii=d00,pikik=p0000;dii<dend;dii++)	// Effective loop over i (via dii)
          {
          hrow = rh;				// Length of 1st hmx row  (row k=0)
          hkk = h00;				// First hmx diagonal: <0|hmx|0>
          for(hkl=h00; hkl<hend; 		// Effective loop over index k (via hkl)
                           pikik+=prow, prow--)	// where prow is the length of row k of hmx
            {
            hkk += hrow;			// <k|hmx|k> -> <k+1|hmx|k+1>
            hrow--;				// Set length of the next hmx row
            for(pikil=pikik;hkl<hkk;		// Effective loop over l (via hkl), l = [k,rh)
                                 hkl++,pikil++)	// where indexes only over the upper triangle elements
              *pikil = (*dii) * (*hkl);		// <ik|pdt|jl> = <i*rh+k|pdt|j*rh+l> = <i|dmx|j><k|hmx|l>
            }					//             = del   <i|dmx|i><k|hmx|l>
          }					//                  i,j
        return pdt;
        } 
      else					// If dmx is complex, then the tensor
        {					// product is not Hermitian, its n_matrix
        n_matrix* pdt =				// Construct a new n_matrix for product
                new n_matrix(dim,dim,complex0);
        complex *p0000 = pdt->data;		// Start of data in pdt: <00|pdt|00>
        complex *d00 = data;			// Start of data in dmx: <0|dmx|0>
        complex *h00 = ((h_matrix*)mx)->data;	// Start of data in hmx: <0|hmx|0>
        complex *dend = d00 + cols_;		// End of data in dmx: <rows_|dmx|rows_+1>
        complex *hend = h00 + (rh*rh+rh)/2;	// End of data in hmx: <rh|hmx|rh+1>
        complex *pikil, *dii, *hkl, *hkk;	// These are matrix elements
        complex *pi0i0, *pikik, *pilik;			
        int hrow;
        for(pi0i0=p0000,dii=d00; dii<dend;	// Effective loop over i (via dii)
                       dii++,pi0i0+=dim*rh+rh)	// where pi0i0 = <i*rh|pdt|i*rh>
	  {
          hrow = rh;				// Length of row 0 in hmx
          hkk = h00;
	  for(hkl=h00,pikik=pi0i0; hkl<hend;	// Effective loop over k (via hkl)
                                 pikik+=dim+1)	// where <ik|pdt|ik> -> <i(k+1)|pdt|i(k+1)>
            {
            pikil = pikik;			// <ik|pdt|il> = <ik|pdt|ik>
	    *pikil = (*dii) * (*hkl);	 	// <ik|pdt|ik> = <i*rh+k|pdt|i*rh+k> = <i|dmx|i><k|hmx|k>
            pilik = pikik + dim;		// <il|pdt|ik> = <i*rh+k+1|pdt|i*rh+k>
            pikil++;				// <ik|pdt|il> = <i*rh+k|pdt|i*rh+k+1>
            hkl++;				// <k|hmx|l> = <k|hmx|k+1>
            hkk += hrow;			// <k|hmx|k> -> <k+1|hmx|k+1>
            hrow --;				// Length of hmx row k+1;
	    for(; hkl<hkk;			// Effective loop over l, l>k (via hkl)
	              pikil++,pilik+=dim,hkl++)
              {
	      *pikil = (*dii) * (*hkl); 	// <ik|pdt|il> = <i|dmx|i><k|hmx|l>, k<l
	      *pilik = (*dii) * conj(*hkl); 	// <il|pdt|ik> = <i|dmx|i><k|hmx*|l>, k<l
              }
            }
          }
        return pdt; 
        }
      break;
      }
    default:					// Tensor product dmx (x) generic matrix
      { 
      int gr = mx->rows();			// Rows in generic matrix
      int gc = mx->cols();			// Columns in generic matrix
      int pr = rows_*gr;			// Rows in product matrix
      int pc = cols_*gc;			// Columns in product matrix
      n_matrix* pdt = new n_matrix(pr, pc);	// Product is a new normal matrix
      int pos=0;				// Position in pdt matrix
      for(int i=0; i<rows_; i++)		// Loop through the rows of dmx
        for(int k=0; k<gr; k++)			// Loop over the rows of gmx
          for(int j=0; j<rows_; j++)		// Loop over the columns of dmx
	    for(int l=0; l<gc; l++,pos++)	// Loop over the columns of gmx
              if(i==j)
	 pdt->data[pos] = data[i]*(*mx)(k,l);	// pdt[pos] = <i*gr+k|pdt|j*gc+l>
      return pdt; 				//          = <i|dmx|j><k|gmx|l> = del   <k|gmx|l>
      } 					//                                    i,j
    }
  }


// ____________________________________________________________________________
// L                     CLASS D_MATRIX I/O FUNCTIONS
// ____________________________________________________________________________

/* These functions perform both ASCII and Binary input/output operations on a
   diagonal matrix.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   print      ostream     Writes array elemnts in ASCII to output stream
   write      ofstream    Writes array elemnts in Binary to output file stream
   read       ifstream    Read in array from Binary input file stream
   readASC    istream     Read in array from an input stream
   ask          ----      Interactively requests diagonal matrix from user
   
   The binary format of d_matrix has data only, not size information.  The size
   is taken care of in classes matrix/_matrix and should be written prior to
   use of these functions. The data ordering is Re(<i|dmx|j>, Im(<i|dmx|j>
   columns then rows (i.e. row by row.) If the flag=0 then only the diagonals
   are output (GAMMA format) whereas a non-zero flag will force all elements
   (diagonals and zero off-diagonals) to be output.                          */

std::vector<std::string> d_matrix::printStrings(const MxPrint& PFlgs) const
  {
  std::vector<std::string> PStrings;		// Strings we will return
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
  std::string ezer;				// Invisible or 0 element
  std::string b(" ");				// Space between elements
  if(!PFlgs.MxAll) ezer=std::string(elen, ' ');	// This if its invisible
  else             ezer=CenterString("0", elen);	// This if 0 is printed

  std::string efmt = complex::dformat();	// Real/Imag element format



  std::string pline;					// String for matrix row
  complex z;						// Temporary complex
  int i,j;						// Looping indices
  if(complex::normphase())				// Each element as (R,phi)
    {
    switch(ptype)
      {
      default:
      case 0:						//   Printing (R,phi)
        for(i=0; i<rows_; i++ )				//     Loop array rows
          {
          pline = "";					//       Zero row
          for(j=0; j<i; j++) pline += ezer + b;		//       Left of diag elems
          pline += (data[i]).printString();		//       This element
          for(j++; j<rows_; j++) pline += b + ezer;	//       Right of diag elems
          PStrings.push_back(pline);			//       Store this line
          }
        break;
      case 1:						//   Printing R only (real)
        for(i=0; i<rows_; i++ )				//     Loop array rows
          {
          pline = "";					//       Zero row
          for(j=0; j<i; j++)     pline += ezer + b;	//       Left of diag elems
          pline += MxModform(efmt.c_str(),Re(data[i]));	//       This element
          for(j++; j<rows_; j++) pline += b + ezer;	//       Right of diag elems
          PStrings.push_back(pline);			//       Store this line
          }
        break;
      case 2:						//   Print R, phi=0,90 (imag)
         for(i=0; i<rows_; i++ )			//     Loop array rows
          {
          pline = "";					//       Zero row
          for(j=0; j<i; j++) pline += ezer + b;		//       Left of diag elems
          pline += MxModform(efmt.c_str(),Im(data[i]));	//       This element
          for(j++; j<rows_; j++) pline += b + ezer;	//       Right of diag elems
          PStrings.push_back(pline);			//       Store this line
          }
        break;
      }
    }
  else						// Each element as (a, b)
    {
//    switch(ptype)
//      {
      for(i=0; i<rows_; i++ )			// Loop array rows
        {
        pline = "";
        z = data[i];				// Diagonal element
        for(j=0; j<i; j++) pline += ezer;	 	// Left of diag elems
        if(ptype==1)				// Output diagonal as 
          pline += MxModform(efmt.c_str(),Re(z)); 	// a real element
        else if(ptype==2)				// Output diagonal as
          pline += MxModform(efmt.c_str(),Im(z)); 	// an imaginary element
        else  pline += z.printString();         	// Output as complex
        for(j=i+1; j<cols_; j++)			// past the diagonal
          pline += ezer;				// if desired
        PStrings.push_back(pline);
        }
//      }
    }
  return PStrings;
  }

std::vector<std::string> d_matrix::pictureStrings(const MxPrint& PFlgs) const
  {
  if(PFlgs.MxAll) { return printStrings(PFlgs); }	// If all elems
  std::vector<std::string> PStrings;			// What to return
  return PStrings;
  }

void d_matrix::print(std::ostream& ostr, const MxPrint& PFlgs) const
  {
  int ptype = 0;				// Complex elements output
  if(PFlgs.MxRIPrnt)				// If we want just reals
    {						// or just imags, check
    if(is_real())           ptype = 1;		// Real elements output
    else if(is_imaginary()) ptype = 2; 		// Imag elements output
    }

  int    elen;					// A single element length
  switch(ptype)					// Get the element length
    {						// from class complex. This 
    default:					// depends on the output
    case 0: elen = complex::zlength(); break;	// format
    case 1:
    case 2: elen = complex::dlength(); break;
    }
  std::string ezer;					// Invisible or 0 element
  if(!PFlgs.MxAll) ezer = std::string(elen+1, ' ');	// This if its invisible
  else
    ezer = std::string(elen/2, ' ') + std::string("0")
         + std::string(elen-elen/2, ' ');

  std::string efmt = complex::dformat();		// Real/Imag element format
  int clen = 40 - ((elen+1)*rows_-1)/2;		// Space to center 1 line
  std::string sp("");				// Spacer to center a line
  if(clen>0) sp = std::string(clen, ' ');		// Set spacer for centering
  complex z;
  int i,j;
  for(i=0; i<rows_; i++ )			// Loop array rows
    {
    z = data[i];				// Diagonal element
    ostr << sp;					// Space to center
    for(j=0; j<i; j++) ostr << ezer;	 	// Left of diag elems
    if(ptype==1)				// Output diagonal as 
      ostr << MxModform(efmt.c_str(),Re(z)); 	// a real element
    else if(ptype==2)				// Output diagonal as
      ostr << MxModform(efmt.c_str(),Im(z)); 	// an imaginary element
    else ostr << z; 		        	// Output as complex
    if(PFlgs.MxAll)				// Output elements
      for(j=i+1; j<cols_; j++)			// past the diagonal
        ostr << ezer;				// if desired
    ostr << "\n";				// Skip to next line
    }
  }

void d_matrix::picture(std::ostream& ostr, const MxPrint & PFlgs) const
  {
//  char* V = NULL;
  int i,j;
/*
  if(PFlgs.MxPctVal)
    {
    double mx = norm(maxZ());			// Array maximum
    double mn = norm(minZ());			// Array minimum
    double del = mx-mn/rows_;			// Array range
    V = new char[rows_];
    if(del == 0) 
    }
*/
  int rlen = 2*rows_-1;				// Length of 1 row
  std::string sp("");				// Spacer to center row
  int len = 40-rlen/2;				// Length to center row
  if(len>0) sp = std::string(len, ' ');		// Set centering spacer
  for(i=0; i<rows_; i++)			// Loop over array rows
    {
    ostr << sp;					// Output centering spacer
    if(PFlgs.MxAll)				// If writing all elements
      for(j=0; j<i; j++) ostr << "0 ";		// put in values to diag.
    else					// If not, these are just
      for(j=0; j<i; j++) ostr << "  ";		// empty spaces
    if(!norm(data[j]))      ostr << "0 ";	// This is zero diagonal
//    else if(PFlgs.MxPctVal) ostr << V[j] << " ";// This is "relative" value
    else                    ostr << "x ";       // This is std filled diag
    if(PFlgs.MxAll)				// If writing all elements
      for(j=i+1; j<rows_; j++) ostr << "0 ";	// finish row past diagonal
    ostr << "\n";
    }
//  if(PFlgs.MxPctVal) delete [] V;
  }

void d_matrix::write(std::ofstream &fp, int form) const
  {
	// **** changed float to double
  double dr,di,dz=0.0;
  int j;
  for(int i=0; i<rows_; i++)			// Loop over all the rows
    {
    dr = Re(data[i]);				// Re(<i|dmx|i>)
    di = Im(data[i]);				// Im(<i|dmx|i>)
    if(!form)					// Just write diagonals if form
      {						// is GAMMA format (default)
      fp.write((char*)&dr, sizeof(double));
      fp.write((char*)&di, sizeof(double));
      }
    else					// Write all elements if form is
      {						// not GAMMA
      for(j=0; j<cols_; j++ )
        {  
        if(i==j)
          {
          fp.write((char*)&dr, sizeof(double));
          fp.write((char*)&di, sizeof(double));
          }
        else
          {
          fp.write((char*)&dz, sizeof(double));
          fp.write((char*)&dz, sizeof(double));
          }
        }
      }
    }
  return;
  }

void d_matrix::read(std::ifstream &fp)
  {
  float dr,di;
  for(int pos=0; pos<rows_; pos++ )
    {
    fp.read((char*)&dr, sizeof(float));
    fp.read((char*)&di, sizeof(float));
    data[pos] = complex(dr,di);
    }
  }

void d_matrix::readASC(std::istream& istr)
  {
  int i,j;
  istr >> i >> j;
  resize(i,j);
  for(int pos=0; pos<size; pos++ )
    istr >> data[pos]; 
  }
 
void d_matrix::ask( )
  {
  float dr,di;
  for(int i=0; i<rows_; i++ )
    {
    std::cout << "\n\tPlease Input Real and Imaginary Value of <"
         << i << "|mx|" << i << "> [re im]: ";
    std::cin >> dr >> di;
    data[i] = complex(dr,di);
    }
  }
 

// ____________________________________________________________________________
// M                CLASS D_MATRIX MISCELLANEOUS FUNCTIONS
// ____________________________________________________________________________

/* These are functions that don't fit into the other categories.
 
    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   resize      nr, nc     Resizes dmx to have nr rows & nc columns, type change
   copy         ----      Produces a copy of the dmx, allocates memory
   convert       mx       Matrix mx has its values set to those of dmx (no mem)
   IMX          ----      Identity matrix produced from dmx conversion (mem)
   DMX          ----      Diagonal matrix produced from dmx conversion (mem)
   HMX          ----      Hermitian matrix produced from dmx conversion (mem)
   NMX          ----      Complex matrix produced from dmx conversion (mem)
 
 The conversion functions will often drop elements, this is unavoidable. For
 example, one cannot generally convert a diagonal array to an identity array.
 There are no elements lost when using DMX or NMX.  The *MX functions
 will always allocate new memory for the converted matrix from dmx.         */

void d_matrix::resize(int i, int j)
  {
  if(i != j)			// Insure resized to square array
    {
    DMxerror(17, 1);				//   Dimension problems
    DMxerror(6, " Matrix Resize", 1);		//   Matrix multiply problems
    DMxfatal(3, "resize");			//   Fail in times_adjoint function
    }
  if(i != size)			// See if new size is any different
    { 				// If it is different, then
    _matrix::resize(i,j);	//	set the new col & row size
    delete [] data;		//	delete the exisiting data
    size = i;			//	set the new data size
    data = new complex[size];	//	declare a new data array
    }
  }

_matrix* d_matrix::copy() { return new d_matrix(*this); }

void d_matrix::convert(_matrix* mx)
  {
  switch(mx->stored_type())
    {
    case d_matrix_type:				// Conversion of diagonal to diagonal
      (*(d_matrix *)mx)=(*this);		// is unnecessary
      break;
    case i_matrix_type:				// Convert diagonal to identity matrix
      mx->resize(rows_, cols_);			// Set matrix size, loose all <i|dmx|j>
      break;
    case h_matrix_type:				// Convert diagonal to Hermitian
      {
      int dd = rows_;				// Dimension of dmx
      mx->resize(dd,dd);			// Switch/Check size of hmx
      complex *hii = ((h_matrix*)mx)->data;	// Start of data in hmx: <i|hmx|i>-><0|hmx|0>
      complex *dii = data;			// Start of data in dmx: <i|dmx|j>-><0|hmx|0>
      complex *dend = dii + dd;			// End of data in dmx: <dd|dmx|dd+1>
      complex *hend = hii + (dd*dd+dd)/2;	// End of data in hmx: <dd|hmx|dd+1>
      complex *hij=hii;				// <i|hmx|j> -> <0|hmx|0>
      int hrow = dd;				// Row increment hmx: <i|hmx|i> -> <i+1|hmx|i+1>
      for(; hij<hend; hij++)			// Loop over all elements of hmx
        (*hij) = complex0;			// <i|hmx|j> = 0
      for(; dii<dend; dii++,hii+=hrow,hrow--)	// Loop over all elements of dmx
	(*hii) = Re(*dii);			// <i|hmx|i> = Re(<i|dmx|i>)
      }
      break;
    case n_matrix_type:				// Convert diagonal to a normal matrix
      {
      int dd = rows_;				// Dimension of dmx
      mx->resize(dd,dd);			// Switch/Check the size of nmx
      complex *nii = ((n_matrix*)mx)->data;	// Start of data in nmx: <i|nmx|i>-><0|nmx|0>
      complex *nij = nii;			// Element <i|nmx|j>-><0|nmx|0>
      complex *dii = data;			// Start of data in dmx: <i|dmx|j>-><0|hmx|0>
      complex *nend = nii + dd*dd;		// End of data in nmx: <dd|nmx|dd+1>
      complex *dend = dii + dd;			// End of data in dmx: <dd|dmx|dd+1>
      int nrinc = dd+1;				// Row increment nmx: <i|nmx|i> -> <i+1|nmx|i+1>
      for(; nij<nend; nij++)			// Loop over all elements of nmx
        (*nij) = complex0;			// <i|nmx|j> = 0
      for(; dii<dend; dii++,nii+=nrinc)		// Loop over all elements of dmx
	(*nii) = (*dii);			// <i|nmx|i> = <i|dmx|i>
      }
      break;
    default:					// Convert diagonal matrix to unknown type
      {
      int dd = rows_;				// Dimension of dmx
      mx->resize(dd,dd);			// Switch/Check the size of mx 
      complex *dii = data;			// Start of data in dmx: <i|dmx|j>-><0|hmx|0>
      for(int i=0; i<dd; i++, dii++)		// Then set <i|mx|j> = <i|dmx|j>
        {
	for(int j=0; j<dd; j++)
	  (*mx)(i,j) = complex0;		// <i|mx|j> = 0
        (*mx)(i,i) = (*dii);			// <i|mx|i> = <i|dmx|i>
        }
      }
    }
  }

//
// ------------------------ New Conversion Routines ---------------------------
//                  (Been Having Troubles With The Old One)
 
 
i_matrix* d_matrix::IMX() { return new i_matrix(cols_, cols_); }
d_matrix* d_matrix::DMX() { return this; }

h_matrix* d_matrix::HMX() 
  {						// Only dmx reals will transfer!
  int id = rows_;				// Dimension of dmx
  h_matrix* hmx = new h_matrix(id,id,complex0);	// Start with a zero matrix
  complex *hii = ((h_matrix*)hmx)->data;	// Start of data in hmx: <i|hmx|i>-><0|hmx|0>
  complex *dii = data;				// Start of data in dmx: <i|dmx|i>-><0|dmx|0>
  complex *hend = hii + (id*id+id)/2;		// End of data in hmx: <id|hmx|id+1>
  int hrow = id;				// Row inc. hmx: <i|hmx|i> -> <i+1|hmx|i+1>
  for(; hii<hend ;hii+=hrow,hrow--,dii++)	// Loop over diagonal elements of hmx
    (*hii) = Re(*dii);				// <i|hmx|i> = Re(<i|dmx|i>
  return hmx;
  }
  
n_matrix* d_matrix::NMX() 
  {
  int id = rows_;				// Dimension of dmx
  n_matrix* nmx = new n_matrix(id,id,complex0);	// Start with a zero matrix
  complex *nii = ((n_matrix*)nmx)->data;	// Start of data in nmx: <i|nmx|i>-><0|nmx|0>
  complex *dii = data;				// Start of data in dmx: <i|dmx|i>-><0|dmx|0>
  complex *nend = nii + id*id;			// End of data in nmx: <id|nmx|id+1>
  int nrinc = id+1;				// Row inc. nmx: <i|nmx|i> -> <i+1|nmx|i+1>
  for(; nii<nend; nii+=nrinc,dii++)		// Loop over diagonal elements of nmx
    (*nii) = (*dii);				// <i|nmx|i> = <i|dmx|i>
  return nmx;
  }

#endif							// d_matrix.cc
