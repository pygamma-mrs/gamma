/* i_matrix.cc **************************************************-*-c++-*-
**								   	**
**                                G A M M A				**
**						   			**
**	Identity Matrix	                           Implementation	**
**						                    	**
**	Copyright (c) 1990, 1991, 1992		                    	**
**	Tilo Levante, Scott A. Smith		                    	**
**	Eidgenoessische Technische Hochschule	                    	**
**	Labor fuer physikalische Chemie		                    	**
**	8092 Zuerich / Switzerland		                    	**
**						                    	**
**      $Header: $
**						                    	**
*************************************************************************/

/*************************************************************************
**								   	**
**  Description				 	  			**
**						   			**
**  The class i_matrix defines identity matrices 			**
**  for C++ with the usual mathematical operations			**
**  +, -, *, / and in/output routines.					**
**						   			**
*************************************************************************/
     
#ifndef Gi_matrix_cc_			// Is file already included?
#  define Gi_matrix_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Matrix/i_matrix.h>		// Include the interface
#include <Matrix/_matrix.h>		// Include the base matrix
#include <Matrix/n_matrix.h>		// Include normal matrices
#include <Matrix/d_matrix.h>		// Include diagonal matrices
#include <Matrix/h_matrix.h>		// Include hermitian matrices
#include <fstream>			// Include libstdc++ filestreams
#include <string>			// Include libstdc++ strings


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                  CLASS I_MATRIX CHECKING & ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   mx      : An identity matrix (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : string in message
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void i_matrix::IMxerror(int eidx, int nr) const
  {
  std::string CL="Identity Matrix";
  Mxerror(CL, eidx, nr);
  }

void i_matrix::IMxerror(int eidx, const std::string& PN, int nr) const
  {
  std::string CL="Identity Matrix";
  Mxerror(CL, eidx, PN, nr);
  }
 
volatile void i_matrix::IMxfatal(int eidx) const
  {
  std::string CL="Identity Matrix";
  Mxfatality(CL, eidx);
  }
 
volatile void i_matrix::IMxfatal(int eidx, const std::string& PN) const
  {
  std::string CL="Identity Matrix";
  Mxfatality(CL, eidx, PN);
  }

// ____________________________________________________________________________
// ii                      CLASS I_MATRIX CHECKING
// ____________________________________________________________________________

bool i_matrix::CheckDims(_matrix* mx, int warn) const
  {
  if(cols_!=mx->rows() || cols_!=mx->cols())    // Insure mx dimension match
    {
    if(warn>0)
      {
      IMxerror(51,1);                           // Array dimensions mismatched
      IMxerror(31,1);                           // Row->Row Col-Col mismatch
      if(warn>1) IMxfatal(81);                  // Cannot continue on
      }
    return false;
    }
  return true;
  }

// ____________________________________________________________________________
// A                 CLASS I_MATRIX CONSTRUCTORS AND DESTRUCTOR
// ____________________________________________________________________________

/*  Arguments                              Constructed Array
    ---------           -------------------------------------------------------
        -                An empty array
     nr, nc              An (nr x nc) array, MUST have nr=nc
    nr, nc, z            An (nr x nc) array, MUST have nr=nc & z=1
       imx               A duplicate of Identity array imx
 
   Note that all constructors invariably call the base class constructor
   (i.e. the one in _matrix) which sets the row & column sizes.  The added
   code does only checking as there is no data array to allocate. Similarly,
   the destructor does nothing since there is no data array, anything else
   is destroyed in the base class.                                           */

i_matrix::i_matrix( )            : _matrix() { }
i_matrix::i_matrix(int i, int j) : _matrix(i,j)
  { if(i!=j) { IMxerror(1,1); IMxfatal(2); } }
i_matrix::i_matrix(int i, int j, const complex& z) : _matrix(i,j)
  {
  if(i!=j)           { IMxerror(1,1); IMxfatal(2);  }    // Insure mx square
  if (z != complex1) { IMxerror(1,1); IMxfatal(50); }    // Insure diago. 1
  }
i_matrix::i_matrix(const i_matrix& imx) : _matrix(imx) { }
i_matrix::~i_matrix () { }


// ____________________________________________________________________________
// B                    CLASS I_MATRIX ACCESS AND ASSIGNMENT
// ____________________________________________________________________________

/* There are two flavors of element access operators herein. The one that might
   commonly be used, i.e. mx(i,j), is DIRECTLY accessing the element <i|mx|j>.
   That can save CPU time by avoiding value copying when altering or obtaining
   elements.  However its misuse can lead to trouble in this class because no
   identity matrix elements are stored and there are restrictions on the
   elements so that the array remains Identity.  An example of when this is
   trouble would be "imx(0,0)=complexi;" and "imx(2,0)=7".  The former fails as
   diagonal elements MUST be 1 and the latter fails because off-diagonal
   elements are not stored so they don't exist.

   In constrast, the "get" function returns copies of the element & the "put"
   function checks that the element is being properly set.  Although these
   are slower they are absolutely safe to use.  The put function returns TF
   and will be FALSE if an element is attempted to be set non-identity. Note
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
          =             Current (Identity) array set equal to input imx
    (int,int)           Reference to element <i|hmx|j> (Potential Danger)
    get(int,int)        Copy of element <i|imx|j>      (Safe)
    put(int,int)        Assigns element <i|imx|j>      (Return FALSE if fails)
    put_h(int,int)      Assigns both <i|imx|j> & <j|imx|i> (Return TRUE/FALSE)
    get_block(r,c,R,C)  Returns imx block of size RxC starting at <r|imx|c>
    put_block(r,c,mx)   Places imx into hmx at position <r|imx|c> (TRUE/FALSE)
 
    Again, note that operations which return TF will trigger class matrix to
    take appropriate action when they return FALSE.  This usually means that
    the operation cannot be handled by a Hermitian array so the array type
    is changed to be more generic.  The the operation is reattempted.       */
 
void i_matrix::operator= (const i_matrix& mx)
  { if(this == &mx) return;  _matrix::operator= (mx); }
complex& i_matrix::operator() (int i, int j) { return (i==j)?Zione:Zizero; }
complex  i_matrix::get(int i, int j) const { return (i==j)?complex1:complex0; }
bool     i_matrix::put(const complex& z, int i, int j)
  { 
  if((i==j) && (z==complex1)) return true;					
  return false;					
  } 

bool i_matrix::put_h(const complex& z, int i, int j) 
 { if((i==j) && (z==complex1)) return true; return false; }


_matrix* i_matrix::get_block(int row, int col, int nrows, int ncols)
 
        // Input      imx : An i_matrix (this)
        //            row : Row position of block start
        //            col : Column position of block start
        //          nrows : Number of rows in the block
        //          ncols : Number of columns in the block
        // Output         : Returns a (nrows x ncols) matrix
        //                  containing the elements of imx
        // Note           : Matrix type output will vary
        //                  depending upon block requested
        // Note           : Uses internal structure of n_matrix
         
  {
  if(row==col && nrows==ncols)			// See if square block about diagonal
    {
    if(!row && ncols==cols_) return this;	// If block is whole matrix, return imx
    return new i_matrix(nrows, ncols);		// Else return a smaller i_matrix
    }
  else
    {
    int i;
    n_matrix* mx = new n_matrix(nrows,		// Return a n_matrix block is
		          ncols,complex0);	// rectangular or off-diagonal
    if(row == col)				// If diagonals of imx and nmx
      {						// line up, just copy diagonal	
      int end = (ncols>nrows)?nrows:ncols;
      for(i=0; i<end; i++)	 		// Loop over rows or columns
        ((n_matrix*)mx)->data[i*ncols+i]=complex1; // <i|nmx|i> = <i|imx|i> = 1
      }
    else
      {
      int I,J,j;
      for(i=0; i<nrows; i++)
        for(j=0,I=i+row,J=col; j<ncols; j++,J++)
          if(I == J) 				// <i|nmx|j> = <i+row|nmx|j+col>
            ((n_matrix*)mx)->data[i*ncols+j]	//           = <I|imx|J> = del   * 1
                                       = complex1; //                            I,J
      }
    return mx;
    }
  }


bool i_matrix::put_block(int row, int col, _matrix* mx)
 
        // Input      imx : An i_matrix (this)
        //            row : Row position of block start
        //            col : Column position of block start
        //             mx : A pointer to a matrix
        // Output      TF : Matrix imx is modified to contain the block
        //                  of data in matrix mx.  TRUE indicates all
        //                  O.K., FALSE signals class matrix to
        //                  handle this function directly
 
  {
  if(mx->stored_type() == i_matrix_type)        // Can put in i_matrix if input 
    if(row == col)                              // block on the diagonal and not
      {                                         // too big
      if(row+mx->rows() > rows_)               // This if too big!
        { 
        IMxerror(52,1); 
        IMxfatal(21,"put_block");
        }
      return true;					// Put O.K., but doesn't change mx 
      }
  return false;					// For other puts must do dmx->nmx 
  }


// ____________________________________________________________________________
// C                  CLASS I_MATRIX HERMITIAN_TYPE HANDLING
// ____________________________________________________________________________
 
// Note: For functions set_hermitian & check_hermitian see class matrix!
 
hermitian_type i_matrix::stored_hermitian( )      const { return _hermitian; }
hermitian_type i_matrix::test_hermitian(double d) const
                                                 { return _hermitian; d=0.0; }
 
// ____________________________________________________________________________
// D                 CLASS I_MATRIX MATRIX_TYPE HANDLING
// ____________________________________________________________________________

/* These functions handle the checking of what type of array we have.  These
   types directly correspond to the matrix classes (such as i_matrix) derived
   from class _matrix.

     Function            Output                       Description
   ----------------  --------------  ------------------------------------------
   set_type            --------      Found in class matrix
   test_type           --------      Found in class matrix 
   stored_type       i_matrix_type   Alway returns we are Hermitian
   test_hermitian    *_matrix_type   Type imx could be within d
   mxtype            string          Returns the string "Identity" 
 
   The test type looks to see if imx could be an array of the matrix type
   specified as an input argument.  If it can within "d" the return type 
   is equal to the input type.  If it cannot the return is i_matrix_type     */

matrix_type i_matrix::stored_type( ) const { return i_matrix_type; }
matrix_type i_matrix::test_type(const matrix_type m, const double d) const	
                                           { return m; double dtmp; dtmp=d; }
std::string i_matrix::mxtype()        const { return std::string("Identity"); }
std::string i_matrix::mxtype(bool pf) const { return std::string("Identity"); }

// ____________________________________________________________________________
// E                   CLASS I_MATRIX VARIOUS MATRIX CHECKS
// ____________________________________________________________________________

/*   Function    Output                       Description
   ------------  ------  ------------------------------------------------------
   is_symmetric   bool   TF if Im(<i|imx|j>) < d for all i,j (TRUE,  d unused)
   is_hermitian   bool   TF if Im(<i|imx|j>) < d for all i,j (TRUE,  d unused)
   is_unitary     bool   TF if inv(imx) == adjoint imx       (TRUE)
   is_real        bool   TF if Im(<i|imx|j>) < d for all i,j (TRUE,  d unused)
   is_imaginary   bool   TF if Re(<i|imx|j>) < d for all i,j (FALSE, d unused)
   is_complex     bool   TF if is_real && is_imaginary       (FALSE, d unused)
   is_zero        bool   TF if ||<i|imx|j>|| < d for all i,j (FALSE, d unused)
   is_diagonal    bool   TF if ||<i|imx|j>|| < d for all i!=j(TRUE,  d unused)
   is_square      bool   TF if rows_ == cols_                (TRUE)
   is_equal       bool   TF if ||<i|imx-mx|j>||<d    all i,j  

   These checks are consistent with other matrix classes.  The value d is
   not used and no real math is performed execpt in the "is_equal" function */
 
bool i_matrix::is_symmetric(double d) const { return true;  d = 0.0; } 
bool i_matrix::is_hermitian(double d) const { return true;  d = 0.0; } 
bool i_matrix::is_unitary(double d)   const { return true;  d = 0.0; } 
bool i_matrix::is_real(double d)      const { return true;  d = 0.0; } 
bool i_matrix::is_imaginary(double d) const { return false; d = 0.0; }
bool i_matrix::is_complex(double d)   const { return false; d = 0.0; }
bool i_matrix::is_zero(double d)      const { return false; d = 0.0; }
bool i_matrix::is_diagonal(double d)  const { return true;  d = 0.0; }
bool i_matrix::is_square( )           const { return true; }

bool i_matrix::is_equal(_matrix* mx, double d) const
  {
  if(cols_ != mx->cols())                                return false;
  if(rows_ != mx->rows())                                return false;
  return ((*mx).test_type(i_matrix_type, d) == i_matrix_type);
  }

// ____________________________________________________________________________
// F              CLASS I_MATRIX BINARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________
 
/* These functions allow for simple mathematical operations between a identity
   matrix and another matrix of unspecified type. Note that these functions
   return a pointer to the resulting array. That array must have it referencing
   done in the main matrix class or there will be memory leakage.

   Function    Result                   Details
   --------   --------   ---------------------------------------------------
     add       imx+mx    Adjusts to mx type, sums only over active elements
   subtract    imx-mx    Adjusts to mx type, subtact only over active elems.
   multiply    imx*mx    Reduced calculation using only nonzero imx elements
   multiply    imx*z     Fast calculation only over imx stored elements

*/


_matrix* i_matrix::add(_matrix* mx)
  {
  if(!CheckDims(mx, 1)) IMxfatal(20);		// Quit if mismatch
  switch (mx->stored_type())
    {
      case n_matrix_type:			// Add a n_matrix to an i_matrix
	{ 
	n_matrix* sum =				// Construct new n_matrix equal to mx
                 new n_matrix(*(n_matrix*)mx);	// <i|sum|j> = <i|mx|j>
	for(int i=0; i<cols_; i++)		// <i|sum|i> += <i|imx|i> = 1
          ((n_matrix*)sum)->data[i*cols_+i]
                                      += complex1;
	return sum;
	}
	break;
      case h_matrix_type:			// Add a h_matrix to an i_matrix
	{
	h_matrix* sum				// Construct a new h_matrix equal to mx
                = new h_matrix(*(h_matrix*)mx);
	for(int i=0; i<cols_; i++)		// <i|sum|j> = <i|hmx|j>
	  sum->data[i*cols_-(i*(i-1))/2]	// <i|sum|i> += <i|imx|i> = 1
                                      += complex1;
	return sum;
	}
      case d_matrix_type:			// Add a d_matrix to an i_matrix
	{	
	d_matrix* sum				// Construct new d_matrix equal to mx
                = new d_matrix(*(d_matrix*)mx);
	for(int i=0; i<cols_; i++)		// <i|sum|j> = <i|dmx|j>
	  sum->data[i] += complex1;		// <i|sum|i> += <i|imx|i> = 1
	return sum;
	}
	break;
      case i_matrix_type:			// Add two i_matrices
	return new d_matrix(cols_,		// Return is a diagonal array
                             cols_,complex(2));	// <i|sum|j> = del   * 2.0
	break;					//                i,j
      default:					// Add a generic matrix to an i_matrix
	{ 
	n_matrix* sum = new n_matrix(cols_,cols_);
	int pos = 0;
	for(int i=0; i<cols_; i++)
	  for(int j=0; j<cols_; j++,pos++)
	    if(i==j)
	      sum->data[pos]=(*mx)(i,j)+complex1;
	    else
	      sum->data[pos]=(*mx)(i,j);
	return sum;
	}
    }
  }


_matrix* i_matrix::subtract(_matrix* mx)
  {
  if(!CheckDims(mx, 1)) IMxfatal(21);		// Quit if mismatch
  switch(mx->stored_type())
    {
    case n_matrix_type:			// Subtract a n_matrix from an i_matrix
	{ 
        n_matrix* dif=(n_matrix *)mx->negate();	// <i|dif|j> = -<i|mx|j>
        for(int i=0; i<cols_; i++)		// <i|dif|i> += <i|imx|i> = 1 
          ((n_matrix*)dif)->data[i*cols_+i]
                                      += complex1;
	return dif;
	}
	break;
      case h_matrix_type:			// Add a h_matrix to an i_matrix
	{
        h_matrix* dif =				// <i|dif|j> = -<i|mx|j>
                      (h_matrix *)mx->negate();
        for(int i=0; i<cols_; i++)		// <i|dif|i> += <i|imx|i> = 1
          dif->data[i*cols_-(i*(i-1))/2]
                                      += complex1;
        return dif;
        }
        break;
      case d_matrix_type:			// Subtract a d_matrix from an i_matrix
	{ 
	d_matrix* mx1 = new d_matrix(cols_, cols_);
	for(int i=0; i<cols_; i++)
	  mx1->data[i] = 1.0 - ((d_matrix*)mx)->data[i];
	return mx1;
	}
	break;
      case i_matrix_type:			// Subtract an i_matrix from an i_matrix;
	return new d_matrix(cols_,cols_,complex0);// Left with a zero matrix
	break;
      default:					// Generic matrix subtracted from i_matrix
	{ 
	n_matrix* mx1 = new n_matrix(cols_,cols_);
	int pos = 0;
	for(int i=0; i<cols_; i++)
	  for(int j=0; j<cols_; j++,pos++)
	    if(i != j)
	      mx1->data[pos] =- (*mx)(i,j);
	    else
	      mx1->data[pos] = complex1 - (*mx)(i,j);
	return mx1;
	}
    }
  }


 
        // Input                imx  : Input identity matrix (this)
        //                       mx  : Second matrix
        // Output               pdt  : New matrix which is the product imx*mx
        // Note                      : The matrix mx is returned as multiplication
        //                             by the identity matrix leaves it unchanged.
        //                             Matrix dimensions are checked though.

_matrix* i_matrix::multiply(_matrix* mx)
  {
  if(cols() != mx->rows())
    {
    IMxerror(17, 1);                            //   Dimension problems
    IMxerror(6, " mx1 * mx2", 1);		//   Matrix multiply problems
    IMxfatal(3, "multiply");			//   Fail in times_adjoint function
    }
  return mx;					// I*mx = mx
  }

        // Input                imx  : Input identity matrix (this)
        //                      z    : Complex number
        // Output               pdt  : New matrix which is the product z*imx

_matrix* i_matrix::multiply(const complex& z)
  {
  if(z == complex1) return this;			// I*1 = I
  else return new d_matrix(rows(),cols(),z);	// I*z = z*I -> diagonal
  }


 
        // Input                imx  : Input identity matrix (this)
        //                       mx  : Second matrix                         -1
        // Output               pdt  : New matrix which is the product imx*mx
        // Note                      : Still under construction
        // Note                      : This uses the function inv
 
//                            ---                -1
//                <i|pdt|j> = \   <i|imx|k> <k|mx  |j>
//                            /
//                            ---

_matrix* i_matrix::divide(_matrix* mx)
  {
  if(cols_ != mx->rows())			// Test if size compatible
    {
    IMxerror(17, 1);                            //   Dimension problems
    IMxerror(6, " mx1 / mx2", 1);		//   Matrix divide problems
    IMxfatal(3, "divide");			//   Fail in divide function
    }
  else
    switch(mx->stored_type())
      {
      case i_matrix_type: return this; break;   // Divide imx into imx
      case d_matrix_type:                       // Divide imx by dmx
        return (d_matrix*)mx->inv(); break;     // Return inv(dmx)
      case h_matrix_type:                       // Divide imx by hmx
        return (h_matrix*)mx->inv(); break;     // Return inv(hmx)
      case n_matrix_type:                       // Divide imx by nmx
        return (n_matrix*)mx->inv(); break;     // Return imx*inv(nmx)
      default:                                  // Divide imx by generic mx
        return mx->inv();			// Return inv(mx)
        break;
      }   
    IMxerror(6, " mx1 / mx2", 1);		//   Matrix divide problems
    IMxerror(3, " divide", 1);			//   Matrix divide problems
    IMxfatal(25);				//   Not fully implemented
    return this;	                        // Compiler confused about return
  }


 
        // Input                imx  : Input identity matrix (this)
        //                      z    : Complex number
        // Output               pdt  : New matrix which is the product (1/z)*imx
        // Note                      : This handles divisions in code such as mx/z

_matrix* i_matrix::divide(const complex& z)
  {
  if(z == complex1) return this;		// No change: I/1 = I
  else if(z==complex0)
    {
    IMxerror(18, 1);				//   Divide by zero
    IMxerror(6, " mx1 / z", 1);			//   Matrix divide problems
    IMxfatal(3, " divide");			//   Matrix divide problems
    }
  return new d_matrix(rows(),cols(),1/z);	// I/z = (1/z)*I -> diagonal
  }
 
// ____________________________________________________________________________
// G                  CLASS I_MATRIX UNARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________
 

 
        // Input                imx  : Input identity matrix (this)
        //                      mx   : Second matrix
        // Output               sum  : Matrix mx is added into imx
        // Note                      : This is the implementation of mx += imx
        // Note                      : Uses internal structure of other matrix types

_matrix* i_matrix::add_two(_matrix* mx)		// Implementation of mx += I
  {
  int id = rows_;				// Dimension of imx
  if((id!=mx->rows()) || (id!=mx->cols()) )	// Check that array dimensions match 
    {
    IMxerror(17, 1);                            //   Dimension problems
    IMxerror(6, " mx1 + mx2", 1);		//   Matrix divide problems
    IMxfatal(3, "add_two");			//   Fail in divide function
    }
  else
    switch(mx->stored_type())
      {
      case n_matrix_type:			// Add an i_matrix into a n_matrix
	{
        complex *nii = ((n_matrix*)mx)->data;	// Start of nmx: nii = <i|nmx|i> -> <0|nmx|0>
        complex *nend = nii + id*id;		// End of nmx: <rows_|nmx|cols|+1>
        int nrow = id+1;			// Amount to take <i|nmx|i> -> <i+1|nmx|i+1>
	for(; nii<nend; nii+=nrow)
	  *nii += complex1;
	return mx;
	}
	break;
      case h_matrix_type:			// Add an i_matrix into an h_matrix
	{
        complex *hii = ((h_matrix*)mx)->data;	// Start of hmx: hii = <i|hmx|i> -> <0|hmx|0>
        complex *hend = hii + (id*id+id)/2;	// End of data in hmx: <id|hmx|id+1>
        int hrow = id;				// Amount to take <i|hmx|i> -> <i+1|hmx|i+1>
	for(; hii<hend; hii+=hrow, hrow--)
	  *hii += complex1;
	return mx;
	}
	break;
      case d_matrix_type:			// Add an i_matrix to a d_matrix
	{ 
        complex *dii = ((d_matrix*)mx)->data;	// Start of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *dend = dii + id;		// End of dmx: <cols_|dmx|cols_+1>
	for(; dii<dend; dii++)
	  *dii += complex1;
	return mx;
	}
	break;
      case i_matrix_type:			// Add an i_matrix to an i_matrix
	return new d_matrix(id,id,complex(2));	// Return d_matrix, <i|dmx|i> = 2
	break;
      default:					// Add an i_matrix to a generic matrix
	{ 
	n_matrix* mx1 = new n_matrix(id,id);	// Construct new nmx for result
	int pos = 0;
	for(int i=0; i<id; i++)
	  for(int j=0; j<id; j++,pos++)
	    if(i == j)
	      mx1->data[pos]=1+(*mx)(i,j);
	    else
	      mx1->data[pos]=(*mx)(i,j);
	return mx1;
	}
      }
  return mx;					// Always have some return
  }


 
        // Input                imx  : Input identity matrix (this)
        //                      mx   : Second matrix
        // Output               sum  : Matrix mx has imx subtracted from it
        // Note                      : This is the implementation of mx -= imx
        // Note                      : Uses internal structure of other matrix types

_matrix* i_matrix::subtract_two(_matrix* mx)	// Implementation of mx -= I
  {
  int id = rows_;				// Dimension of imx
  if((id!=mx->rows()) || (id!=mx->cols()) )	// Check that array dimensions match 
    {
    IMxerror(17, 1);                            //   Dimension problems
    IMxerror(6, " mx1 - mx2", 1);		//   Matrix divide problems
    IMxfatal(3, "subtract_two");		//   Fail in divide function
    }
  else
    switch(mx->stored_type())
      {
      case n_matrix_type:			// Subtract an i_matrix from a n_matrix
	{
        complex *nii = ((n_matrix*)mx)->data;	// Start of nmx: nii = <i|nmx|i> -> <0|nmx|0>
        complex *nend = nii + id*id;		// End of nmx: <rows_|nmx|cols|+1>
        int nrow = id+1;			// Amount to take <i|nmx|i> -> <i+1|nmx|i+1>
	for(; nii<nend; nii+=nrow)
	  *nii -= complex1;
	return mx;
	}
	break;
      case h_matrix_type:			// Subtract an i_matrix from an h_matrix
	{
        complex *hii = ((h_matrix*)mx)->data;	// Start of hmx: hii = <i|hmx|i> -> <0|hmx|0>
        complex *hend = hii + (id*id+id)/2;	// End of data in hmx: <id|hmx|id+1>
        int hrow = id;				// Amount to take <i|hmx|i> -> <i+1|hmx|i+1>
	for(; hii<hend; hii+=hrow, hrow--)
	  *hii -= complex1;
	return mx;
	}
	break;
      case d_matrix_type:			// Subtract an i_matrix from a d_matrix
	{ 
        complex *dii = ((d_matrix*)mx)->data;	// Start of dmx: dii = <i|dmx|i> -> <0|dmx|0>
        complex *dend = dii + id;		// End of dmx: <cols_|dmx|cols_+1>
	for(; dii<dend; dii++)
	  *dii -= complex1;
	return mx;
	}
	break;
      case i_matrix_type:			// Subtract an i_matrix from an i_matrix
	return new d_matrix(id, id, complex0);	// Return d_matrix, <i|dmx|i> = 0
	break;
      default:					// Subtract an i_matrix from a generic matrix
	{ 
	n_matrix* mx1 = new n_matrix(id,id);	// Construct new nmx for result
	int pos = 0;
	for(int i=0; i<id; i++)
	  for(int j=0; j<id; j++,pos++)
	    if(i == j)
	      mx1->data[pos]=(*mx)(i,j)-complex1;
	    else
	      mx1->data[pos]=(*mx)(i,j);
	return mx1;
	}
      }
  return mx;					// Always have some return
  }


 
        // Input                imx  : Input identity matrix (this)
        //                      mx   : Second matrix
        // Output               sum  : Matrix mx is multiplied into imx
        // Note                      : This is the implementation of mx *= imx
        // Note                      : Since mx*imx = mx, mx is returned

_matrix* i_matrix::multiply_two(_matrix* mx)
  {
  if(rows_ != mx->cols()) 			// Check size compatability
    {
    IMxerror(17, 1);                            //   Dimension problems
    IMxerror(6, " mx1 * mx2", 1);		//   Matrix divide problems
    IMxfatal(3, "multiply_two");		//   Fail in divide function
    }
  return mx;					// mx*I = mx
  }


 
        // Input                imx  : Input identity matrix (this)
        //                      z    : Complex number
        // Output               mx   : Matrix imx is multiplied by z
        // Note                      : This implements imx*=z (but not z*=imx)

_matrix* i_matrix::multiply_two(const complex &z)
  {
  if(z != complex1)
    return new d_matrix(rows_,cols_,z);		// I*z -> diagonal if z not 1
  return this;					// I*1 = I
  }


_matrix* i_matrix::divide_two(_matrix* mx) { return mx; }	 // mx/I = mx
_matrix* i_matrix::divide_two(const complex &z)
 
        // Input                imx  : Input identity matrix (this)
        //                      z    : Complex number
        // Output               mx   : Matrix imx is multiplied by (1/z)
        // Note                      : This implements imx/=z (but not z/=imx)

  {
  if(z != complex1)
    return new d_matrix(rows_,cols_, 1/z);	// I/z -> diagonal
  return this;					// I/1 = I
  }

// ____________________________________________________________________________
// H                 CLASS I_MATRIX SIMPLE UNARY FUNCTIONS
// ____________________________________________________________________________
/* These functions perform simple simple mathematical operations on an Identity
   matrix. Note that for Identity arrays many of these do nothing.

    Function          Result             Function           Result
   ---------- -----------------------    --------  ------------------------
   negate     <i|imx|j> -> -<i|imx|j>       RE      <i|imx|j>->Re(<i|imx|j>)
   conjugate  <i|imx|j> -> <i|imx*|j>       IM      <i|imx|j>->Im(<i|imx|j>)
   transpose  <i|imx|j> -> <j|imx|i>      adjoint   <i|imx|j>-><j|imx*|i>    */
 
_matrix* i_matrix::RE()        { return this; }		    // imx is real
_matrix* i_matrix::conjugate() { return this; }		    // imx self adjoint
_matrix* i_matrix::transpose() { return this; }		    // imx self adjoint
_matrix* i_matrix::adjoint()   { return this; }		    // imx self adjoint
_matrix* i_matrix::negate()    { return new d_matrix(rows_,cols_,-complex1); } 
_matrix* i_matrix::IM()        { return new d_matrix(rows_,cols_, complex0); }
_matrix* i_matrix::mxexp()     { return new d_matrix(rows_,cols_, exp(1.0)); }
complex  i_matrix::trace()     { return complex(rows_); }

/*     Function    Output                       Description
     ------------  -------  ---------------------------------------------------
       swaprows      T/F    Swaps two rows of the matrix
       swapcols      T/F    Swaps two columns of the matrix
        permute    matrix   Permutes rows & columns of the matrix
        maxRe      double   Returns largest real value in the array
        maxIm      double   Returns largest imaginary value in the array
        maxZ       complex  Returns largest complex (norm) value in array
        minRe      double   Returns smallest real value in the array
        minIm      double   Returns smallest imaginary value in the array
        minZ       complex  Returns smallest complex (norm) value in array   */

_matrix* i_matrix::swaprows(int i, int j)
  { 
  h_matrix* Isr = new h_matrix(rows_,rows_,complex0);	// Make new h_matrix
  for(int k=0; k<rows_; k++) Isr->put(complex1,k,k); 	// Make it imx
  Isr->put(complex0,  i,i);				// Do the row swap now,
  Isr->put(complex0,  j,j);				// there are only 4
  Isr->put_h(complex1,i,j);				// elements to change
  return Isr;
  }
_matrix* i_matrix::swapcols(int i, int j) { return swaprows(i,j); }
_matrix* i_matrix::permute( int i, int j) { return this; }
double   i_matrix::maxRe()   const { return (rows_)?1:0;                  }
double   i_matrix::maxIm()   const { return 0;                            }
complex  i_matrix::maxZ()    const { return (rows_)?complex1:complex0;    }
double   i_matrix::minRe()   const { return (rows_==1)?1:0;               }
double   i_matrix::minIm()   const { return 0;                            }
complex  i_matrix::minZ()    const { return (rows_==1)?complex1:complex0; }

// ____________________________________________________________________________
// I                 CLASS I_MATRIX SIMPLE BINARY FUNCTIONS
// ____________________________________________________________________________


        // Input            imx : An i_matrix (this)
        //                   mx : A second matrix
        // Output            z  : A complex value wish is the trace of
        //                        the product of two input arrays
        //                         z = Tr{imx*mx} = Tr{mx}
        // Note                 : This is faster than taking the matrix
        //                        product and then using unary trace!

complex i_matrix::trace(_matrix* mx)
  {
  if((rows_ !=mx->cols()) || cols_!=mx->rows())	// Insure dimensions proper
    {
    IMxerror(17, 1);                            //   Dimension problems
    IMxerror(6, " trace(mx1 * mx2)", 1);	//   Matrix divide problems
    IMxfatal(3, "trace");			//   Fail in divide function
    }
  return mx->trace();				// Tr{I*mx} = Tr{mx}
  }

_matrix* i_matrix::adjoint_times(_matrix* mx) { return mx; }
_matrix* i_matrix::times_adjoint(_matrix* mx) { return mx->adjoint(); }

        // Input            imx : An i_matrix (this)
        //                   mx : A second matrix
        // Output           mx1 : adjoint_times produces
        //                                    T *
        //                         mx1 = [(imx ) ] * mx = mx
        //                                        T *      T*
        //                         mx1 = nmx* [(mx ) ] = mx

// ____________________________________________________________________________
// J                  CLASS I_MATRIX COMPLEX UNARY FUNCTIONS
// ____________________________________________________________________________
 
// Note:  The alogorithm for FFT is not found here.  The identity matrix is
//        converted to normal matrix then transformed in class n_matrix.

complex i_matrix::det()  { return complex1; }
int     i_matrix::rank() { return rows_;    }	
 
        // Input            imx : An i_matrix (this)
        // Output             z : Value of the determinant; det{imx}
        // Output             i : Rank of matrix imx
        // Note                 : For i_matrix, det=1, rank=dim


// ____________________________________________________________________________
// K                 CLASS I_MATRIX COMPLEX BINARY FUNCTIONS
// ____________________________________________________________________________
  
// ***************************** tensor product *******************************



        // Input            imx : An i_matrix (this)
        //                   mx : A second matrix
        // Output           pdt : A matrix which is the tensor product
        //                        of the two input matrices
	// Note			: Uses internal structure of other matrix types 

        //                           pdt        =   imx (x) mx

        //                       (m*n x n*o)       (mxn)   (nxo)

        //                    <i*n+k|pdt|j*o+l> = <i|imx|j><k|mx|l>

_matrix* i_matrix::tensor_product (_matrix* mx)
  {
  switch (mx->stored_type())
    {
    case d_matrix_type:				// Tensor product imx (x) dmx
      { 
      int r = rows_;				// Get dimension of imx
      int dr = mx->rows();			// Get dimension of dmx
      d_matrix* pdt = new d_matrix(r*dr, r*dr);	// Consruct a new dmx
      int pos=0;				// Position in pdt: pdt[pos] = <i*dr+k|pdt|i*dr+k>
      for(int i=0; i<r; i++)			// Loop over the rows of imx
	for(int k=0; k<dr; k++,pos++)		// Loop over the rows of dmx
	  pdt->data[pos] = 			// pdt[pos] = <i*dr+k|pdt|j*dr+l> = <i|imx|j><k|dmx|l>
	             ((d_matrix*)mx)->data[k]; 	//          = del   del   <i|imx|i><k|nmx|k> = del       <k|nmx|k>
      return pdt;				//               i,j   k,l                        i,j;k,l
      }
      break;

   case i_matrix_type:				// Tensor product imx (x) dmx
     {						// Return another (larger) imx
     int dim = rows_*mx->rows();		// Get the dimension of the product
     return new i_matrix(dim,dim);
     }
     break;

    case n_matrix_type:				// Tensor product imx (x) nmx
      {
      int ca = cols_;				// Columns of imx
      int cb = mx->cols();			// Columns of nmx
      int rb = mx->rows();			// Rows of nmx
      int cm = cols_*cb;			// Columns in tensor product
      int rm = cols_*rb;			// Rows in tensor product
      n_matrix* pdt = 				// Construct a new normal array
                new n_matrix(rm,cm,complex0);	// which is initialized to zero
      complex *p0000 = pdt->data;		// Start of data in pdt
      complex *n00 = ((n_matrix*)mx)->data;	// Start of data in nmx
      complex *nend = n00 + cb*rb;		// End of data in  nmx
      complex *nrow = NULL;			// Dummy index: loop over row of nmx
      complex *pikil,*nkl,*pi0i0;
      int i;
      for(pi0i0=p0000, i=0; i<ca;		// Effective loop over i, i = [0, ca)
                         i++, pi0i0+=cm*rb+cb)
	for(nkl=n00,pikil=pi0i0; nkl<nend;	// Effective loop over k, k = [0, rb)
                                 pikil+=cm-cb)
	  for(nrow=nkl+cb; nkl<nrow;		// Effective loop over l, l = [0, cb)
                               nkl++,pikil++)
	    *pikil = *nkl;			// <ik|pdt|jl> = <i*rb+k|pdt|j*cb+l> = <i|imx|j><k|nmx|l>
      return pdt;				//             = del   <i|imx|i><k|nmx|l> = del   <k|nmx|l>
      }						//                  i,j                        i,j
      break;
    case h_matrix_type:				// Tensor product imx (x) hmx
      { 
      int r = rows_;				// Get dimension of imx
      int rh = mx->rows();			// Rows of Hermitian matrix
      int dim = r*rh;
      h_matrix* pdt = 				// Construct a new normal array
               new h_matrix(dim,dim,complex0);	// which is initialized to zero
      complex *p0000 = pdt->data;		// Start of data in pdt: <0,0|pdt|0,0>
      complex *h00 = ((h_matrix*)mx)->data;	// Start of data in hmx: <0|hmx|0>
      complex *hend = h00 + (rh*rh+rh)/2;	// End of data in hmx: <rh|hmx|rh+1>
      complex *pikil, *pikik, *hkl, *hkk;	// These are the matrix elements
      int prow = dim;				// Length of first row in pdt
      int i,hrow;
      for(i=0,pikik=p0000; i<r; i++)		// Loop over the rows of imx
        {
        hrow = rh; 				// Length of 1st hmx row  (row k=0)
        hkk = h00;				// First hmx diagonal: <0|hmx|0>
        for(hkl=h00;hkl<hend;pikik+=prow,prow--)// Effective loop over index k
          {
          hkk += hrow;				// <k|hmx|k> -> <k+1|hmx|k+1>
          hrow--;				// Set length of the next hmx row
	  for(pikil=pikik;hkl<hkk;hkl++,pikil++)// Effective loop over l (via hkl), l = [k,rh)
	    *pikil = *hkl;			// <ik|pdt|jl> = <i*rh+k|pdt|j*rh+l> = <i|imx|j><k|hmx|l>
          } 					//             = del   <i|imx|i><k|hmx|l> = del   <k|hmx|l>
        } 					//                  i,j                        i,j
      return pdt;
      }	
      break;
    default:					// Tensor product imx (x) generic matrix
      { 
      int r = rows_;				// Get dimension of imx
      int rg = mx->rows();			// Rows of generic matrix
      int cg = mx->cols();			// Columns of generic matrix
      n_matrix* pdt = 				// Construct a new normal array
               new n_matrix(r*rg,r*cg,complex0);	// which is initialized to zero
      int pos=0;				// Position in pdt: pdt[pos] = <i*rg+k|pdt|j*cg+l>
      for(int i=0; i<r; i++)
        for(int k=0; k<rg; k++)
	  for(int j=0; j<r; j++)
	    for(int l=0; l<cg; l++,pos++)
              if(i==j)
	         pdt->data[pos] = (*mx)(k,l);	// <ik|pdt|jl> = <i*rb+k|pdt|j*cb+l> = <i|imx|j><k|gmx|l>
      return pdt;				//             = del   <i|imx|i><k|gmx|l> = del   <k|gmx|l>
      }						//                  i,j                        i,j
    }
  }
 
// ____________________________________________________________________________
// L                     CLASS I_MATRIX I/O FUNCTIONS
// ____________________________________________________________________________

/* These functions perform both ASCII and Binary input/output operations on an
   identity matrix.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   print      ostream     Writes array elements in ASCII to output stream
   picture    ostream     Writes array elements pictorically to output stream
   write      ofstream    Writes array elements in Binary to output file stream
   read       ifstream    Read in array from Binary input file stream
   readASC    istream     Read in array from an input stream
   ask          ----      Interactively requests complex matrix from user

   The binary format of i_matrix data only, not size information.  The size
   is taken care of in classes matrix/_matrix and should done prior to use of
   these functions.  The data ordering is Re(<i|imx|j>, Im(<i|imx|j> columns
   then rows (i.e. row by row.)                                              */

std::vector<std::string> i_matrix::printStrings(const MxPrint& PFlgs) const
  {
  std::vector<std::string> PStrings;
  if(!PFlgs.MxAll) 			// No unstored elems printed
    {					// So just print out 1
    PStrings.push_back(std::string("1"));
    return PStrings;
    }
  int i,j;				// Element indicies
  std::string pline;			// String for array row
  std::string z("0");
  std::string one("1");
  std::string b(" ");
  for(i=0; i<rows_; i++ )		// Loop rows
    {
    pline = "";
    for(j=0; j<i; j++)
      pline += z + b;
    pline += one;
    for(j++; j<rows_; j++)
      pline += b + z;
    PStrings.push_back(pline);
    }
  return PStrings;
  }

std::vector<std::string> i_matrix::pictureStrings(const MxPrint& PFlgs) const
  {
  if(PFlgs.MxAll) { return printStrings(PFlgs); }	// If all elems
  std::vector<std::string> PStrings;			// What to return
  std::string pline;
  std::string blnk(" ");
  std::string one("1");
  int i, j;
  for(i=0; i<rows_; i++)				// Loop over all rows
    {
    pline = "";
    for(j=0; j<rows_; j++)
      {
      if(j==i) pline += one;
      else     pline += blnk;
      if(j+1 < rows_) pline += blnk;
       }
    PStrings.push_back(pline);
    }
  return PStrings;
  }


void i_matrix::print(std::ostream& ostr, const MxPrint& PFlgs) const
  {
  if(!PFlgs.MxAll) 			// No unstored elems printed
    {					// So just print out 1
    ostr << std::string(39, ' ') << 1;
    return;
    }
  int rlen = 2*rows_-1;			// Length of 1 array row
  std::string sp("");			// Spacer to center row
  int len = 40-rlen/2;			// Spacing to center row
  if(len>0) sp = std::string(len, ' ');	// Set spacer to center
  int i,j;				// Element indicies
  for(i=0; i<rows_; i++ )		// Loop rows
    {
    ostr << sp;
    for(j=0; j<i; j++)       ostr << "0 ";
                             ostr << "1 ";
    for(j=i+1; j<cols_; j++) ostr << "0 ";
    ostr << "\n";
    }
  }

void i_matrix::picture(std::ostream& ostr, const MxPrint& PFlgs) const
  {
  if(PFlgs.MxAll) { print(ostr,PFlgs); return; }
  int rlen = 2*rows_-1;
  std::string sp("");
  int len = 40-rlen/2;
  if(len>0) sp = std::string(len, ' ');
  for(int i=0; i<rows_; i++ )
    {
    ostr << sp;
    for(int j=0; j<i; j++) ostr << "  ";
      ostr << "1 ";
    ostr << "\n";
    }
  }

void i_matrix::write(std::ofstream &fstr, int form) const

        // Input            imx : An i_matrix (this)
        //                 fstr : A file stream
        //                 form : Flag for format in which to write dmx data
        //                               0 = GAMMA format (default); no data
        //                              !0 = All matrix elements including 1's & 0's 
        // Output          fstr : File modified to contain imx data in binary
        // Note                 : For form!=0 the data is written in rows, each point
        //                        as two floats, real followed by imaginary
        //                        Re(<0|imx|0>), Im(<0|imx|0>, Re(0|imx|1>, ...
        //                        where Re(<i|imx|j>) = Im(<i|imx|j>) = 0 when i!=j 
	//			  and Re(<i|imx|j>) =1, Im(i|imx|j>) = 0 when i==j
        // Note                 : For form==0, no elements are written, imx contains
        //                        no explicit data - the matrix values are implicit

  {						// If GAMMA binary form, do nothing
  if(form)					// If not GAMMA binary form, then
    {						// write out all of th elements
    float d0 = 0.0;
    float d1 = 1.0;
    for(int i=0; i<rows_; i++ )
      for(int j=0; j<cols_; j++ )
        {
        if(i==j)
          fstr.write((char*)&d1, sizeof(float));	//	 Diagonal real part
        else
          fstr.write((char*)&d0, sizeof(float));	//	 Off-diagonal real part
        fstr.write((char*)&d0, sizeof(float));		//	 All elements imaginary part
        }
    }
  return;
  }



        // Input            imx : An i_matrix (this)
        //                 fstr : File
        // Output          fstr : File unmodified!
        // Note                 : This is the matching function to
        //                        i_matrix::write!  The write can only
        //                          1.) write a GAMMA i_matrix, which
	//				contains no information so no
	//				write is performed to the file
	//			    2.) write an i_matrix as if all of 
	//				its elements were present in memory
	//			  For read of imx written by 1.), we read
	//			  nothing.  Read of imx written by 2.) is
	//			  not possible since all elements are present,
	//			  it must be read in as a normal matrix!

void i_matrix::read(std::ifstream& fstr)
  {
  if(!fstr.is_open()) 			// Compiler likes fstr to be used
    {
    IMxerror(1, 1);			//   Problems with input filestream
    IMxerror(26, 1); 			//   Binary filestream closed
    IMxfatal(3, "read");		//   Fail in read function
    }
  } 


void i_matrix::readASC(std::istream& istr)

// sosi - this makes no sense to me.  Input i & j for the dimensions
//        of imx.  Then resize forces i==j and i=j=size.  So, unless
//        the imx is already set up it cannot be "read" and if already
//        set up why bother to "read" it?
  {
  int i,j;
  istr >> i >> j;				// Input rows & cols
  resize(i,j);
  }


// ____________________________________________________________________________
// M                 CLASS I_MATRIX MISCELLANEOUS FUNCTIONS
// ____________________________________________________________________________
/* These are functions that don't fit into the othere categories.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   resize      nr, nc     Resizes imx to have nr rows & nc columns, type change
   copy         ----      Produces a copy of the imx, allocates memory
   convert       mx       Matrix mx has its values set to those of imx (no mem)
   IMX          ----      Identity matrix produced from imx conversion (mem)
   DMX          ----      Diagonal matrix produced from imx conversion (mem)
   HMX          ----      Identity matrix produced from imx conversion (no mem)
   NMX          ----      Complex matrix produced from imx conversion  (mem)
   
   The conversion functions will often drop elements, this is unavoidable. For
   example, one cannot generally convert an HMX array to a identity array.
   There are no elements lost when using HMX or NMX.  The *MX functions
   will always allocate new memory for the converted matrix from hmx.      */ 

void i_matrix::resize(int i, int j)
  {
  if(i != j)			// Insure resized to square array
    {
    IMxerror(17, 1);                            //   Dimension problems
    IMxerror(6, " Matrix Resize", 1);           //   Matrix resize problems
    IMxfatal(3, "resize");                      //   Fail in resize function
    }
  if(i != size)
    {
    _matrix::resize(i,j);	// Set the new col & row size
    size = i;			// Reset the total size
    }
  }


_matrix* i_matrix::copy() { return new i_matrix(*this); }
void i_matrix::convert(_matrix* mx)
  {
  switch (mx->stored_type())
    {
    case i_matrix_type:				// Conversion of imx to identity is
      (*(i_matrix *)mx) = (*this);		// not necessary
      break;
    case d_matrix_type:				// Convert imx to a diagonal matrix
      {
      int id = rows_;				// Dimension of imx
      mx->resize(id,id);			// Set/Check matrix size of dmx
      complex *dii = ((d_matrix*)mx)->data;     // Start of data in hmx: <i|hmx|i>-><0|hmx|0>
      complex *dend = dii + id;			// End of data in dmx: <dd|dmx|dd+1>
      for(; dii<dend; dii++)			// Loop through all elements of dmx
        (*dii) = complex1;
      }
      break;
    case h_matrix_type:				// Conversion of imx to Hermitian matrix
      {
      int id = rows_;				// Dimension of imx
      mx->resize(id,id);			// Set/Check matrix size of hmx
      complex *hii = ((h_matrix*)mx)->data;     // Start of data in hmx: <i|hmx|i>-><0|hmx|0>
      complex *hend = hii + (id*id+id)/2;       // End of data in hmx: <id|hmx|id+1>
      complex *hij=hii;                         // <i|hmx|j> -> <0|hmx|0>
      int hrow = id;                            // Row increment hmx: <i|hmx|i> -> <i+1|hmx|i+1>
      for(; hij<hend; hij++)			// Loop over all elements of hmx
        (*hij) = complex0;			// <i|hmx|j> = 0
      for(; hii<hend ;hii+=hrow,hrow--)		// Loop over all elements of imx
        (*hii) = complex1;				// <i|hmx|i> = 1
      }
      break;
    case n_matrix_type:				// Conversion of imx to normal matrix
      {
      int id = rows_;				// Dimension of imx
      mx->resize(id,id);			// Set/Check matrix size of nmx
      complex *nii = ((n_matrix*)mx)->data;	// Start of data in nmx: <i|nmx|i>-><0|nmx|0>
      complex *nij = nii;			// Element <i|nmx|j>-><0|nmx|0>
      complex *nend = nii + id*id;		// End of data in nmx: <id|nmx|id+1>
      int nrinc = id+1;				// Row increment nmx: <i|nmx|i> -> <i+1|nmx|i+1>
      for(; nij<nend; nij++)			// Loop over all elements of nmx
        (*nij) = complex0;			// <i|nmx|j> = 0
      for(; nii<nend; nii+=nrinc)		// Loop over all elements of imx
        (*nii) = complex1;				// <i|nmx|i> = 1
      }
      break;
    default:					// Conversion of imx to unknown matrix
      {
      int id = rows_;                           // Dimension of imx
      mx->resize(id,id);                        // Switch/Check the size of mx
      for(int i=0; i<id; i++)			// Then set <i|mx|j> = <i|imx|j>
        {
        for(int j=0; j<id; j++)
          (*mx)(i,j) = complex0;			// <i|mx|j> = 0
        (*mx)(i,i) = complex1;			// <i|mx|i> = <i|imx|i>
        }
      }
    }
  }

 
//
// ------------------------ New Conversion Routines ---------------------------
//                  (Been Having Troubles With The Old One)


i_matrix* i_matrix::IMX() { return this; }
d_matrix* i_matrix::DMX() { return new d_matrix(cols_,cols_,complex1); }
h_matrix* i_matrix::HMX()
  {
  int id = rows_;				// Dimension of imx
  h_matrix* hmx = new h_matrix(id,id,complex0);	// Start with a zero matrix
  complex *hii = ((h_matrix*)hmx)->data;	// Start of data in hmx: <i|hmx|i>-><0|hmx|0>
  complex *hend = hii + (id*id+id)/2;		// End of data in hmx: <id|hmx|id+1>
  int hrow = id;				// Row inc. hmx: <i|hmx|i> -> <i+1|hmx|i+1>
  for(; hii<hend ;hii+=hrow,hrow--)		// Loop over diagonal elements of hmx
  (*hii) = complex1;				// <i|hmx|i> = 1
  return hmx;
  }


n_matrix* i_matrix::NMX()
  {
  int id = rows_;				// Dimension of imx
  n_matrix* nmx = new n_matrix(id,id,complex0);	// Start with a zero matrix
  complex *nii = ((n_matrix*)nmx)->data;	// Start of data in nmx: <i|nmx|i>-><0|nmx|0>
  complex *nend = nii + id*id;			// End of data in nmx: <id|nmx|id+1>
  int nrinc = id+1;				// Row inc. nmx: <i|nmx|i> -> <i+1|nmx|i+1>
  for(; nii<nend; nii+=nrinc)			// Loop over diagonal elements of nmx
    (*nii) = complex1;				// <i|nmx|i> = 1
  return nmx;
  }

// ____________________________________________________________________________
// N                CLASS I_MATRIX DIAGONALIZATION FUNCTIONS
// ____________________________________________________________________________

/* These functions reflect individual steps used in general complex matrix
   diagonalizations. Since how a diagonalization is performed depends upon the
   underlying matrix type.  Of course, in the case of identity matrices, such
   routines don't do a damn thing. They must exist but will just return imx
   back.

   Function   Arguments                         Result
 ------------ --------- ------------------------------------------------------
  BlockDiag    mx, mx   Get blocked array form:          mx = U * BD  * inv(U)
 SymTriDiag    mx, mx   Get symmetric tridiagonal form:  mx = U * STD * inv(U)
 HermTriDiag   mx, mx   Get Hermitian tridiagonal array: mx = U * HTD * inv(U)
   SymDiag     mx, mx   Diagonalize symmetric array:     mx = U * SD  * inv(U)
 Diagonalize   mx, mx   Diagonalize array:               mx = U *  D  * inv(U)
    diag       mx, mx   Diagonalize array (non-member):  mx = U *  D  * inv(U)

 In all cases, the original array will remain unchanged and both the array U
 and "diagonalized" matrix will return as (pointers to) identity matrices.
 The input arrays should be NULL entering this function, and there referencing
 done upon return.  Since the return here is invariably the same matrix, all
 of this amounts to just incrementing the "this" imx reference count by 2.   */

std::vector<int> i_matrix::BlockDiag(_matrix*    (&BD), std::vector<int> &U)
  { BD = this;
  int nr = rows_;				// Matrix dimension
  int i;
  for(i=0; i<nr; i++) { U.push_back(i); }	// Set both as unpermuted
  return std::vector<int>(1, rows_); }
void        i_matrix::HermTriDiag(_matrix* (&HTD), _matrix* (&U))
  { HTD = this; U = this; }
void        i_matrix::SymTriDiag(_matrix*  (&STD), _matrix* (&U))
  { STD = this; U = this; }
void        i_matrix::SymDiag(_matrix*       (&D), _matrix* (&U))
  { D = this;   U = this; }
void        i_matrix::diag(_matrix *         (&D), _matrix * (&U))
  { D = this;   U = this; }

// ____________________________________________________________________________
// O                    CLASS I_MATRIX INVERSION FUNCTIONS
// ____________________________________________________________________________
 

_matrix* i_matrix::inv() { return this; }
	
        // Input            imx : An i_matrix (this)
        // Output           imx : Returns imx as an identity matrix is
        //                        it's own inverse
        //                             -1
        //                        I = I  * I = I * I

 _matrix* i_matrix::LU(int *indx) { indx[0] = -1; return this; }

        // Input            imx : An i_matrix (this)
	//		    indx: Row permutation array
	// Output           iLU : LU decomposition of a row permutation
	//			  of the input matrix imx, imx', the row
	//			  exchanges recorded in array indx
	//			  iLU(low.tri.) * iLU(upp.tri.) = imx'
	//			         L      *       U       = imx'

 _matrix* i_matrix::LUinv(int *indx, _matrix* LU)

        // Input            I   : An i_matrix (this) of results
	//		    LU  : An LU decomposition of matrix A
	//			  (or row permutation of matrix A)
	//		    indx: Row permutation array
	// Output           X   : Solution to the equation: A*X = I
	//			  where LU is the LU decomposition of A
	//			  or row permutation of A
	// Note		        : Matrix LU must be square
	// Note		        : If indx=NULL, X will be inv(A)
	// Note		        : Uses internal structure of other matrix types

  {
  if(!(LU->is_square())) 			// Insure LU is a square array
    {
    IMxerror(17, 1);                                    //   Dimension problems
    IMxerror(6, "LU Inversion On Rectangular Mx", 1);   //   LU inversion multiply problems
    IMxfatal(3, "LUinv");                               //   Fail in times_adjoint function
    }
  if(rows_ != LU->rows())			// Insure LU and B dimensions match
    { 						//     LU(mxm)*X(mxm) = B(mxm)
    IMxerror(17, 1);                                    //   Dimension problems
    IMxerror(6, "LU Inversion LU|X>=|I> Mismatch", 1); //   LU inversion multiply problems
    IMxfatal(3, "LUinv");                               //   Fail in times_adjoint function
    }
  switch(LU->stored_type())
    {
    case i_matrix_type:				// LU is an identity matrix
      if(indx[0] < 0)				// If no row permutations, LU = A = I
        return this;				// so then A*X=I -> I*X=I; X=I
      break;
    case d_matrix_type:				// LU is a diagonal matrix
      if(indx[0] < 0)				// If no row permutations, LU = A = dmx 
        return LU->inv();			// so then A*X=I -> dmx*X=I; X=inv(dmx);
      break;
    case h_matrix_type:				// LU is a Hermitian matrix
    case n_matrix_type:				// LU is a normal matrix
    default:
      break;
    }
// sosi - This had a memory leak! So, I put a kludge in here to patch it.
//        Now, I just build an identity matrix called I as n_matrix by hand, use
//        the n_matrix routine LUinv, and delete I before returning the inverse.
//        Obviously this should be fixed to use a real I matrix as it will be much
//        faster!  12/21/96
 
  int nd = rows_;
  n_matrix* I = new n_matrix(nd,nd,complex0);   // Constuct a new n_matrix
  for(int i=0; i<nd; i++) I->put(complex1,i,i); // Make it the identity matrix
  _matrix* luinv = I->LUinv(indx,LU);
  delete I;
  return luinv;
  }

#endif								// i_matrix.cc
