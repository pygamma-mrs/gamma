/* _matrix.h ****************************************************-*-c++-*-
**									**
**                              G A M M A				**
**									**
**	Prototype Matrix                    Interface Definition	**
**									**
**	Copyright (c) 1990, 1991, 1992					**
**	Tilo Levante & Scott A. Smith					**
**	Eidgenoessische Technische Hochschule				**
**	Labor fuer physikalische Chemie					**
**	8092 Zuerich / Switzerland					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** Description								**
**									**
** The class _matrix defines a "dummy matrix" in C++ for GAMMA.  This 	**
** serves as a proto-type and front end to all classes which handle the	**
** specific matrix types.  The externally used matrix class, matrix,	**
** routes through this 	dummy class in finding specialized routines.	**
**									**
*************************************************************************/

#ifndef   G_matrix_h_			// Is file already included?
#  define G_matrix_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Include system specifics
#include <iostream>			// Inlcude libstdc++ IO streams
#include <Matrix/complex.h>		// Include complex numbers
#include <string>			// Include libstdc++ strings
#include <vector>			// Include libstdc++ STL vectors

const double GMxCut = 1.e-12;		// General machine cutoff

struct MxPrint 				// Matrix ASCII output info
  {
  bool MxHdr;				// Flag if header is output
  bool MxRIPrnt;			// Flag if print real vs imag
  bool MxAll;				// Flag if all elemts output
  int  MxPctDim;			// When to use pictorial
  bool MxPctVal;			// Flag if pictorial qualif.
  int  VxCols;				// Columns in row vector print
  int  VxRows;				// Rows in column vector print

  MxPrint(bool A, bool B, bool C, int D, bool E, int F, int G)
    { MxHdr=A; MxRIPrnt=B; MxAll=C; MxPctDim=D; MxPctVal=E; VxCols=F; VxRows=G; } 
  };

const MxPrint MxPrDefs(true,true,true,5,true,4,30);

enum matrix_type {  _matrix_type,	// base class
		    n_matrix_type,	// normal matrix class
		    d_matrix_type,	// diagonal matrix class
		    i_matrix_type,	// identity matrix class
		    h_matrix_type	// hermitian matrix class
		 };

/* The matrix_type above is used for each derived class of _matrix. It
   is this flag which will often be used to determine how a derived
   matrix behaved.  The "enum" simply assigns a value to each type.
   So, _matrix_type==0, n_matrix_type==1,.... and so on.   The ordering
   isn't particularly important because all tests are to with the names
   rather than the number values.                                       */

enum hermitian_type { non_hermitian=0, _hermitian=1 };

/* The hermitian_type above is also used for each derived class of _matrix.
   Although GAMMA has a complex Hermitian matrix class, some matrices that
   are not stored as Hermitian may indeed be Hermitian (e.g. an I matrix).
   In that case, hermitian_type flag can be set to indicate that fact.  */

class n_matrix;				// Know about derived normal arrays
class d_matrix;				// Know about derived diagonal arrays
class i_matrix;				// Know about derived identity arrays
class h_matrix;				// Know about derived Hermitian arrays

class _matrix {
                int rows_,		// Number of rows in array
		    cols_,		// Number of columns in array
	            size,		// Number of stored elements
	            references_;	// Number of references to this matrix
					// will be used by class matrix


friend class n_matrix;			// Allow normal arrays full access
friend class d_matrix;			// Allow diagonal arrays full access
friend class i_matrix;			// Allow identity arrays full access
friend class h_matrix;			// Allow Hermitian arrays full access

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                      CLASS _MATRIX ERROR HANDLING
// ____________________________________________________________________________
 
/*      Input                   mx      : A _matrix (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : string in message
        Output                  none    : Error message output
                                          Program execution stopped if fatal */
 
         void _MxError(                     int eidx, int nr=0) const;
         void _MxError(int eidx,const std::string& PN, int n=0) const;
volatile void _MxFatal(                  int eidx)              const;
volatile void _MxFatal(int eidx, const std::string& PN)         const;

         void Mxerror(const std::string& CL, int eidx, int n=0) const;
         void Mxerror(const std::string& CL, int eidx,
                                        const std::string& PN, int nr=0) const;
volatile void Mxfatality(const std::string& CL,  int eidx)      const;
volatile void Mxfatality(const std::string& CL,  int eidx,
                                         const std::string& PN)          const;
public:
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A               CLASS _MATRIX CONSTRUCTORS AND DESTRUCTOR
// ____________________________________________________________________________
 
        _matrix( );
        _matrix(int i, int j);
        _matrix(const _matrix& m);
virtual ~_matrix();

// ____________________________________________________________________________
// B                   CLASS _MATRIX ACCESS AND ASSIGNMENT
// ____________________________________________________________________________
 

// ************************* assignment for _matrix ***************************

void operator= (const _matrix& mx);

        // Input                _mx	: A _matrix (this)
	//			 mx	: Another _matrix
        // Output       	void	: _mx is set equal to mx


// ************************ _matrix access functions **************************

/* These elements either get or set a paricular element.  They are designed
   with a global matrix structure in mind, that is the external global matrix
   class handles matrix manipulations. Thus the access operator (i,j) is here
   allowed. Normally this would be dangerous since elements many matrix types
   are not actually stored and thus cannot be set with this operator. However,
   we will BE CAREFUL within the _matrix derived classes to not ruin the 
   structure with misuse of the access function won't we?  Note that the put
   functions will return T/F depending upon their success.                   */

virtual complex& operator() (int i, int j);
virtual complex  get(int i, int j) const;
virtual bool     put(const complex& z, int i, int j);
virtual bool     put_h(const complex& z, int i, int j);

virtual _matrix* get_block(int row, int col, int nr, int nc);

        // Input                _mx	: A _matrix (this)
        //			z 	: A complex number
	//			row, col: Postion where the block starts
	//			nr, nc	: Size of the block
	// Output		_mx*	: Returns a block of dimension (rowxcol)
	// 				  taken from _mx starting at <row|_mx|col>


virtual bool put_block(int row, int col, _matrix* mx);

        // Input                _mx     : A _matrix (this)
        //                      row,col : Row and column indices
        //                      m       : Another _matrix
        // Output               TF      : True if _mx is properly modified
        // Note                         : Matrix m is put into array _mx
        //                                starting at position <row|_mx|col>

// ******************* _matrix details access functions ***********************

const int  rows();		     // Number of matrix rows
const int  cols();		     // Number of matrix columns
const int  refs();		     // Number of matrix references
const int  pts();		     // Number of matrix points
      int& references();	     // THE actual reference count

 
// ____________________________________________________________________________
// C                  CLASS _MATRIX HERMITIAN_TYPE HANDLING
// ____________________________________________________________________________

/* These functions handle the checking of whether an array is of a Hermitian
   type.  Since base matrices have no type, the tests always fail and don't
   actually use the input value d.  These are somewhat outdated now that
   GAMMA has a Hermitian array class and will probably be removed in favor of
   a TF is_hermitian function function or something.  More appropriate would be
   an stored_unitary function.....

     Function            Output                       Description
   ----------------  --------------  ------------------------------------------
   stored_hermitian  hermitian_type  Returns _hermitian=TRUE
   test_hermitian    hermitian_type  Returns _hermitian=TRUE (d is unused)   */


virtual hermitian_type stored_hermitian( ) const; 
virtual hermitian_type test_hermitian(double d=GMxCut) const;

// ____________________________________________________________________________
// D                   CLASS _MATRIX MATRIX_TYPE HANDLING
// ____________________________________________________________________________

/* These functions handle the checking of what type of array we have.  These
   types directly correspond to the matrix classes (such as i_matrix) derived
   from class _matrix.

     Function            Output                       Description
   ----------------  --------------  ------------------------------------------
   set_type            --------      Found in class matrix
   test_type           --------      Found in class matrix 
   stored_type        _matrix_type   Alway returns we are null
   test_hermitian    *_matrix_type   Type _mx could be within d
   mxtype            string          Returns the string "Null" 
 
   The test type looks to see if _mx could be an array of the matrix type
   specified as an input argument.  If it can within "d" the return type 
   is equal to the input type.  If it cannot the return is _matrix_type     */

// Note:  Functions set set_type & check_type are in class matrix, not here.

virtual matrix_type stored_type( ) const;
virtual matrix_type
                 test_type(matrix_type m=n_matrix_type, double d=GMxCut) const;
virtual std::string mxtype() const;
virtual std::string mxtype(bool pf) const;

// ____________________________________________________________________________
// E                   CLASS _MATRIX VARIOUS MATRIX CHECKS
// ____________________________________________________________________________


/*   Function    Output                       Description
   ------------  ------  ------------------------------------------------------
   is_unitary     int    TF if inv(_mx) == adjoint mx, CPU intensive
   is_real        int    TF if Im(<i|_mx|j>) < d for all i,j (def. d GMxCut)
   is_imaginary   int    TF if Re(<i|_mx|j>) < d for all i,j (def. d GMxCut)
   is_complex     int    TF if is_real && is_imaginary       (def. d GMxCut)
   is_zero        int    TF if ||<i|_mx|j>|| < d for all i,j (def.d GMxCut)
   is_equal       int    TF if ||<i|_mx-mx|j>||<d    all i,j (def.d GMxCut)  */

virtual bool is_symmetric(double          d = GMxCut) const;
virtual bool is_hermitian(double          d = GMxCut) const;
virtual bool is_unitary(double            d = GMxCut) const;
virtual bool is_real(double               d = GMxCut) const;
virtual bool is_imaginary(double          d = GMxCut) const;
virtual bool is_complex(double            d = GMxCut) const;
virtual bool is_zero(double               d = GMxCut) const;
virtual bool is_diagonal(double           d = GMxCut) const;
virtual bool is_square()                              const;
virtual bool is_equal(_matrix* mx, double d = GMxCut) const;

// ____________________________________________________________________________
// F                CLASS _MATRIX BINARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________
 
/* All Arithmetic Function Are Handled By Derived Classes of _matrix.  Thus,
   Access of These ALWAYS Produces A Fatal Error.  Users CANNOT Perform Any
   Arithmetic Operations On NULL Matrices Which Have No Dimension Or Data.   */
 
// --------------------- Matrix With Matrix Operations ------------------------
 
virtual _matrix* add(_matrix* mx);
virtual _matrix* subtract(_matrix* mx);
virtual _matrix* multiply(_matrix* mx);
virtual _matrix* divide(_matrix* mx);

        // Input                _mx     : A _matrix (this)
        //                      mx	: Pointer to another matrix
        // Ouput                mx1	: Pointer to the {addition, difference,
	//				  product, division } array of the two
	//				  input matrices
	// Note				: Order is _mx-mx, _mx*mx, _mx/mx
 
// --------------------- Matrix With Scalar Operations ------------------------

virtual _matrix* multiply(const complex& z);
virtual _matrix* divide(const complex& z);

        // Input                _mx     : A _matrix (this)
	//			z	: A complex number
        // Ouput                mx1	: Pointer to the array which is _mx
	//				  {multiplied, divided} by z

// ____________________________________________________________________________
// G                 CLASS _MATRIX UNARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________

/* These functions allow for simple mathematical operations between two         
   matrices and between a matrix and a scalar.
 
                    Function     Arguments          Result 
                    ------------ --------- -------------------------
                    add_two       _mx,mx   mx+_mx (mx += this)
                    subtract_two  _mx,mx   mx-_mx (mx -= this)
                    multiply_two  _mx,mx   mx*_mx (mx *= this)
                    multiply_two  _mx,z    z*_mx  (mx = z*this)
                    divide_two    _mx,mx   mx*_mx (mx *= inv(this))
                    divide_two    _mx,z    z*_mx  (mx = z*inv(this))        */

virtual _matrix* add_two ( _matrix* mx );
virtual _matrix* subtract_two( _matrix* mx );
virtual _matrix* multiply_two( _matrix* mx);
virtual _matrix* multiply_two( const complex &z );
virtual _matrix* divide_two( _matrix* mx );
virtual _matrix* divide_two( const complex &z );

// ____________________________________________________________________________
// H                  CLASS _MATRIX SIMPLE UNARY FUNCTIONS
// ____________________________________________________________________________
 
/* These functions perform simple simple mathematical operations on a Hermitian
   matrix. Note that for Hermitian arrays the adjoint function does nothing.
 
    Function          Result             Function           Result
   ---------- -----------------------    --------  ------------------------
   negate     <i|_mx|j> -> -<i|_mx|j>       RE      <i|_mx|j>->Re(<i|_mx|j>)
   conjugate  <i|_mx|j> -> <i|_mx*|j>       IM      <i|_mx|j>->Im(<i|_mx|j>)
   transpose  <i|_mx|j> -> <j|_mx|i>      adjoint   <i|_mx|j>-><j|_mx*|i>   */

virtual _matrix* negate();
virtual _matrix* RE();
virtual _matrix* IM();
virtual _matrix* conjugate();
virtual _matrix* transpose();
virtual _matrix* mxexp();
virtual _matrix* adjoint();
virtual complex  trace();

/*     Function    Output                       Description
     ------------  -------  ---------------------------------------------------
       swaprows      T/F    Swaps two rows of the matrix
       swapcols      T/F    Swaps two columns of the matrix
        permute    matrix   Permutes rows & columns of the matrix
        maxRe      double  Returns largest real value in the array
        maxIm      double  Returns largest imaginary value in the array
        maxZ       complex Returns largest complex (norm) value in array
        minRe      double  Returns smallest real value in the array
        minIm      double  Returns smallest imaginary value in the array
        minZ       complex Returns smallest complex (norm) value in array    */

virtual _matrix* swaprows(int i, int j);
virtual _matrix* swapcols(int i, int j);
virtual _matrix* permute(int  i, int j);
virtual double   maxRe()   const;
virtual double   maxIm()   const;
virtual complex  maxZ()    const;
virtual double   minRe()   const;
virtual double   minIm()   const;
virtual complex  minZ()    const;

// ____________________________________________________________________________
// I                  CLASS _MATRIX SIMPLE BINARY FUNCTIONS
// ____________________________________________________________________________


virtual complex  trace( _matrix *mx );
virtual _matrix* adjoint_times(_matrix* mx);
virtual _matrix* times_adjoint(_matrix* mx);

// ____________________________________________________________________________
// J                  CLASS _MATRIX COMPLEX UNARY FUNCTIONS
// ____________________________________________________________________________


virtual _matrix* inv();

        // Input        this: a matrix
        // Output       a pointer to a new matrix


virtual _matrix* LU(int* indx);

        // Input        this: a matrix
	//			indx : Row permutation array
	// Output               mx2  : LU decomposition of a row permutation
	//			       of the input matrix mx, mx', the row
	//			       exchanges recorded in array indx
	//			       mx2(low.tri.) * mx2(upp.tri.) = mx'
	//			              L      *       U       = mx'


virtual _matrix* LUinv(int *indx, _matrix* LU);

        // Input            B   : A matrix (this) of results
	//		    LU  : An LU decomposition of matrix A
	//			  (or row permutation of matrix A)
	//		    indx: Row permutation array
	// Output           X   : Solution to the equation: A*X = B
	//			  where LU is the LU decomposition of A
	//			  or row permutation of A
	// Note		        : Matrix LU must be square


// ____________________________________________________________________________
// N                CLASS _MATRIX DIAGONALIZATION FUNCTIONS
// ____________________________________________________________________________

/*                             Diagonalization Routines

           Input        _mx	: A null matrix (this)
                        D	: Pointer to mx to become diagonal mx 
                        U	: Pointer to mx to become eigenvector mx 
           Output       void    : The matrices D and U are set to 
                                  the diagonal matrix of _mx eigenvalues & 
                                  matrix of _mx eigenvectors respectively 
           Note                 : Just a placeholder, Disallowed             */
 
virtual std::vector<int> BlockDiag(_matrix*   (&BF),  std::vector<int> &U);
virtual void             HermTriDiag(_matrix* (&HTD), _matrix* (&U)); 
virtual void             SymTriDiag(_matrix*  (&STD), _matrix* (&U)); 
virtual void             SymDiag(_matrix*     (&D),   _matrix* (&U));

virtual void diag(_matrix* (&D),          _matrix* (&U));

virtual complex det();

        // Input        this: a matrix
        // Output       the determinat of the matrix


virtual int rank();

        // Input        this: a matrix
        // Output       the rank of the matrix


// ____________________________________________________________________________
// K                CLASS _MATRIX COMPLEX BINARY FUNCTIONS
// ____________________________________________________________________________

// ****************************** tensor product ******************************

virtual _matrix* tensor_product ( _matrix* mx );

        // Input        this, mx: two matrices
        // Output               : a matrix of size
        //                        (this->rows()*mx->rows(),
        //                         this->cols()*mx->cols()) with the elemts:
        //                        (i*mx->rows()+j,k*mx->cols()+l)=
        //                         this(i,k)*mx(j,l)


// ____________________________________________________________________________
// L                      CLASS _MATRIX I/O FUNCTIONS
// ____________________________________________________________________________

/* These functions perform both ASCII and Binary input/output operations on a
   _matrix.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   print      ostream     Writes array elements in ASCII to output stream
   picture    ostream     Writes array elements pictorically to output stream
   write      ofstream    Writes array elements in Binary to output file stream
   read       ifstream    Read in array from Binary input file stream
   readASC    istream     Read in array from an input stream
   ask          ----      Interactively requests complex matrix from user 
    
   The binary format of _matrix has no data, nor size information.  Thus, we
   should never output these things in binary format.                        */

virtual std::vector<std::string> printStrings(const   MxPrint& pf) const;
virtual std::vector<std::string> pictureStrings(const MxPrint& pf) const;

virtual void print(std::ostream&   ostr, const MxPrint& PFlgs) const;
virtual void picture(std::ostream& ostr, const MxPrint& PFlgs) const;
virtual void write(std::ofstream &fp, int form=0) const;
virtual void read(std::ifstream& fp);
virtual void readASC(std::istream& istr);
virtual void ask( );
 
// ____________________________________________________________________________
// M                   CLASS _MATRIX MISCELLANEOUS FUNCTIONS
// ____________________________________________________________________________

virtual void resize(int i, int j);

	// Input                _mx  	: A _matrix (this)
        //			i,j	: Row, column dimensions
        // Output		none	: Size of _mx readjusted to ixj
        // Note				: For _mx, general failure
         

virtual _matrix* copy();

	// Input                _mx  	: A _matrix (this)
        // Output               mx	: Pointer to _matrix
        // Note				: This returns a new _matrix on exit.
        //				  Unlike class matrix where assignment
        //				  involves only referencing, a new copy
        //				  of the matrix data is done here
 

virtual void convert ( _matrix* mx );
 
	// Input                _mx  	: A _matrix (this)
        //			mx	: A matrix of any type
        // Output		mx 	: Elements of _mx are placed into mx
        //				  in a fashion so as to maintain the mx
        //				  matrix type.  The result is the "best"
        //				  conversion of _mx into type of mx
 

//
// ------------------------ New Conversion Routines ---------------------------
//                  (Been Having Troubles With The Old One)

virtual i_matrix* IMX();
virtual d_matrix* DMX();
virtual h_matrix* HMX();
virtual n_matrix* NMX();

        // Input                _mx     : A _matrix (this)
        // Output               imx     : A new i/d/h/n_matrix matrix whose
        //                                elements are the same as _mx
        // Note                         : These will die


virtual _matrix* convert(matrix_type mxtype);
 
	// Input                _mx  	: A _matrix (this)
        //			mxtype	: A matrix type
        // Output		mx 	: A new matrix of type mxtype is made
	//				  whose elements are the same as _mx
        //				  if possible.

};

static       complex ZNULL(0,0);		   // Complex changable value

#endif
