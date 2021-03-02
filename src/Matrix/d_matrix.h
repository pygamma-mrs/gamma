/* d_matrix.h ***************************************************-*-c++-*-
**									**
**                               G A M M A				**
**									**
**	Diagonal Matrix                      Interface Definition	**
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
**                                                                      **
** Description                                                          **
**                                                                      **
** The class d_matrix defines diagonal complex matrices for C++ with    **
** the usual operations +, -, *, / and I/O routines.                    **
**                                                                      **
*************************************************************************/

#ifndef   Gd_matrix_h_			// Is file already included?
#  define Gd_matrix_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <Matrix/_matrix.h>		// Know base matrix type
#include <Matrix/complex.h>		 // Know GAMMA complex class
#include <string>			// Know libstdc++ strings
#include <vector>			// Know libstdc++ STL vectors

static complex Zdzero(0);		  // Complex public changable 0

class d_matrix : public _matrix 
  {
  complex* data;			// The actual stored data
  bool real;				// Flag if real array
  bool imag;				// Flag if imaginary array 
      
// ----------------------------------------------------------------------------
// ---------------------------- PRIVATE FUNCTIONS -----------------------------
// ----------------------------------------------------------------------------
// ____________________________________________________________________________
// i                      CLASS D_MATRIX ERROR HANDLING
// ____________________________________________________________________________
 
/*      Input                   dmx     : A diagonal matrix (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : string in message
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

         void DMxerror(int idx,                        int nr=0) const; 
         void DMxerror(int idx, const std::string& PN, int nr=0) const;
volatile void DMxfatal(int idx)                                  const;
volatile void DMxfatal(int idx, const std::string& PN)           const;

// ____________________________________________________________________________
// ii                     CLASS D_MATRIX CHECKING
// ____________________________________________________________________________

virtual bool CheckDims(_matrix* mx,                int warn=1);
virtual bool CheckDim(_matrix*  mx, bool mul=true, int warn=1);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

  friend class n_matrix;
  friend class i_matrix;
  friend class h_matrix;
  
// ____________________________________________________________________________
// A              CLASS D_MATRIX CONSTRUCTORS AND DESTRUCTOR
// ____________________________________________________________________________

/*  Arguments                      Constructed Array
    ---------           -------------------------------------------------------
        -                An empty array
     nr, nc              An (nr x nc) array, uninitialized
    nr, nc, z            An (nr x nc) array, sets <i|mx|i>=Re(z) & <i<j|mx|j>=z
       dmx               A duplicate of diagonal array dmx
 
   Note that all constructors invariably call the base class constructor
   (i.e. the one in _matrix) which sets the row & column sizes.  The added
   code should simply allocate the data array of proper size. Similarly, the
   destructor only destroys the data array, anything else is destroyed in
   the base class.                                                           */

        d_matrix( );
        d_matrix(int i, int j);
        d_matrix(int i, int j, const complex& z);
        d_matrix(const d_matrix& m);
virtual ~d_matrix();

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

void             operator= (const d_matrix& dmx);
virtual complex& operator() (int i, int j);
virtual complex  get(int i, int j) const;
virtual bool     put (const complex& z, int i, int j);
virtual bool     put_h(const complex& z, int i, int j);
virtual _matrix* get_block(int row, int col, int nrows, int ncols);
virtual bool     put_block(int row, int col, _matrix* mx );

// ____________________________________________________________________________
// C                 CLASS D_MATRIX HERMITIAN_TYPE HANDLING
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
 
// Note:  Functions set_hermitian & check_hermitian are found in class matrix

virtual hermitian_type stored_hermitian ( )              const;
virtual hermitian_type test_hermitian(double d = GMxCut) const;
 
// ____________________________________________________________________________
// D                  CLASS D_MATRIX MATRIX_TYPE HANDLING
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

// Note:  Functions set set_type & check_type are in class matrix, not here.
 

virtual matrix_type stored_type( ) const;
virtual matrix_type test_type(matrix_type m = n_matrix_type,
                                                        double d=GMxCut) const;
virtual std::string mxtype()        const;
virtual std::string mxtype(bool pf) const;

// ____________________________________________________________________________
// E                  CLASS D_MATRIX VARIOUS MATRIX CHECKS
// ____________________________________________________________________________

/*   Function    Output                       Description
   ------------  ------  ------------------------------------------------------
   is_symmetric   bool   TF if Im(<i|dmx|j>) < d for all i,j (def. d GMxCut)
   is_hermitian   bool   TF if Im(<i|dmx|j>) < d for all i,j (def. d GMxCut)
   is_unitary     bool   TF if inv(_mx) == adjoint mx, CPU intensive
   is_real        bool   TF if Im(<i|dmx|j>) < d for all i,j (def. d GMxCut)
   is_imaginary   bool   TF if Re(<i|dmx|j>) < d for all i,j (def. d GMxCut)
   is_complex     bool   TF if is_real && is_imaginary       (def. d GMxCut)
   is_zero        bool   TF if ||<i|dmx|j>|| < d for all i,j (def. d GMxCut)
   is_diagonal    bool   TF if ||<i|dmx|j>|| < d for all i!=j(TRUE,d GMxCut)
   is_square      bool   TF if rows_ == cols_                (TRUE)
   is_equal       bool   TF if ||<i|dmx-mx|j>||<d    all i,j (def. d GMxCut) */

virtual bool is_symmetric(double          d=GMxCut) const;
virtual bool is_hermitian(double          d=GMxCut) const;
virtual bool is_unitary(double            d=GMxCut) const;
virtual bool is_real(double               d=GMxCut) const;
virtual bool is_imaginary(double          d=GMxCut) const;
virtual bool is_complex(double            d=GMxCut) const;
virtual bool is_zero(double               d=GMxCut) const;
virtual bool is_diagonal(double           d=GMxCut) const;
virtual bool is_square( )                           const;
virtual bool is_equal(_matrix* mx, double d=GMxCut) const;

// ____________________________________________________________________________
// F              CLASS D_MATRIX BINARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________
 
/* These functions allow for simple mathematical operations between two         
   matrices and between a matrix and a scalar.
 
   Function Arguments      Result        Function Arguments       Result
   -------- --------- -----------------  -------- --------- -------------------
     add     dmx,mx        dmx+mx        multiply  dmx,mx         dmx*mx
   subtract  dmx,mx        dmx-mx        multiply  dmx,z          z*dmx
    divide   dmx,mx        dmx*inv(mx)    divide   dmx,z         (1/z)*mx
      /      z,mx         z*inv(mx)         *       d,mx         d*inv(mx)  */
 
virtual _matrix* add(_matrix* mx);
virtual _matrix* subtract(_matrix* mx);
virtual _matrix* multiply( _matrix* mx );
virtual _matrix* multiply(const complex& z);
virtual _matrix* divide(_matrix* mx);
virtual _matrix* divide(const complex& z);

// ____________________________________________________________________________
// G                 CLASS D_MATRIX UNARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________

/* These functions allow for simple mathematical operations between two         
   matrices and between a matrix and a scalar.
 
                    Function     Arguments          Result 
                    ------------ --------- -------------------------
                    add_two       dmx,mx   mx+dmx (mx += this)
                    subtract_two  dmx,mx   mx-dmx (mx -= this)
                    multiply_two  dmx,mx   mx*dmx (mx *= this)
                    multiply_two  dmx,z    z*dmx  (mx = z*this)
                    divide_two    dmx,mx   mx*dmx (mx *= inv(this))
                    divide_two    dmx,z    z*dmx  (mx = z*inv(this))        */

virtual _matrix* add_two(_matrix* mx);
virtual _matrix* subtract_two( _matrix* mx ); 
virtual _matrix* multiply_two( _matrix* mx );
virtual _matrix* multiply_two(const complex &z );
virtual _matrix* divide_two( _matrix* mx );
virtual _matrix* divide_two(const complex &z );

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


virtual _matrix* negate();
virtual _matrix* RE();
virtual _matrix* IM();
virtual _matrix* conjugate();
virtual _matrix* transpose();
virtual _matrix* adjoint();
virtual _matrix* mxexp();
virtual complex  trace();

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

virtual _matrix*         swaprows(int i, int j);
virtual _matrix*         swapcols(int i, int j);
virtual _matrix*         permute( int i, int j);
virtual double           maxRe()                            const;
virtual double           maxIm()                            const;
virtual complex          maxZ()                             const;
virtual double           minRe()                            const;
virtual double           minIm()                            const;
virtual complex          minZ()                             const;

// ____________________________________________________________________________
// I                CLASS D_MATRIX SIMPLE BINARY FUNCTIONS
// ____________________________________________________________________________


virtual complex trace(_matrix* mx);

        // Input            dmx : A d_matrix (this)
        //                   mx : A second matrix
        // Output            z  : A complex value which is the trace of
        //                        the product of two input arrays
        //                         z = Tr{dmx*mx} = Tr{mx*dmx}
        // Note                 : This is faster than taking the matrix
        //                        product and then using unary trace!

 
virtual _matrix* adjoint_times(_matrix* mx);

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


virtual _matrix* times_adjoint(_matrix* mx);

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
 
 
// ____________________________________________________________________________
// J                  CLASS D_MATRIX COMPLEX UNARY FUNCTIONS
// ____________________________________________________________________________
 

virtual _matrix* inv();
      
        // Input            dmx : A d_matrix (this)
        // Output          dinv : Returns the inverse of dmx as a
        //                        diagonal matrix
 
        //                               -1
        //                        I = dmx   * dmx = dinv * dmx


virtual _matrix* LU(int *indx);

        // Input            dmx : An d_matrix (this)
	//		    indx: Row permutation array
	// Output           dLU : LU decomposition of a row permutation
	//			  of the input matrix dmx, dmx', the row
	//			  exchanges recorded in array indx
	//			  dLU(low.tri.) * dLU(upp.tri.) = dmx'
	//			         L      *       U       = dmx'


virtual _matrix* LUinv(int *indx, _matrix* LU);

        // Input            B   : A d_matrix (this) of results
	//		    LU  : An LU decomposition of matrix A
	//			  (or row permutation of matrix A)
	//		    indx: Row permutation array
	// Output           X   : Solution to the equation: A*X = B
	//			  where LU is the LU decomposition of A
	//			  or row permutation of A
	// Note		        : Matrix LU must be square
	// Note		        : Uses internal structure of other matrix types


// ____________________________________________________________________________
// N                CLASS D_MATRIX DIAGONALIZATION FUNCTIONS
// ____________________________________________________________________________

/*                             Diagonalization Routines

  This Don't Do Much Obviously. Returned Diagonal & Eigenvectors Arrays
  Must Have Their Referencing Set External To This Class.

           Input        dmx     : A diagonal matrix (this)
                        mxd     : Pointer to mx to become diagonal mx
                        U	: Pointer to mx to become eigenvector mx
           Output       void    : The matrices mxd and mxev are set to
                                  the diagonal matrix of imx eigenvalues and
                                  the matrix of imx eigenvectors respectively
           Note                 : For d_matrix, both mx & mx1 are d_matrix
           Note                 : The reference count to dmx must be twice
                                  incremented external by the call origin.   */


virtual std::vector<int> BlockDiag(_matrix*    &BD, std::vector<int> &U);
virtual void             HermTriDiag(_matrix* &HTD, _matrix* &U);
virtual void             SymTriDiag(_matrix*  &STD, _matrix* &U);
virtual void             SymDiag(_matrix*       &D, _matrix* &U);

virtual void diag(_matrix* &D,                    _matrix* &U);
 
virtual complex det();
 
        // Input            dmx : A d_matrix (this)
        // Output             z : Value of the determinant; det{dmx}

//                            _____                                  
//                      det =  | |  <i|dmx|i>
//                             | |


  virtual int rank();

        // Input            dmx : A d_matrix (this)
        // Output             i : The value of the rank of dmx
        // Note                 : For a diagonal matrix, the rank
        //                        is equivalent to the number of non-zero
        //                        rows (or non-zero columns)

    
// ____________________________________________________________________________
// K              CLASS D_MATRIX COMPLEX BINARY FUNCTIONS
// ____________________________________________________________________________
 
// ***************************** tensor product *******************************


virtual _matrix* tensor_product(_matrix* mx);

        // Input            dmx : A d_matrix (this)
        //                   mx : A second matrix
        // Output           pdt : A matrix which is the tensor product
        //                        of the two input matrices
 
        //                           pdt        =   dmx (x) mx
 
        //                       (m*n x m*o)       (mxm)   (nxo)
 
        //                    <i*n+k|pdt|j*o+l> = <i|dmx|j><k|mx|l>
 

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
   (diagonals and zero off-diagonals) to be output.                         */

virtual std::vector<std::string> printStrings(const   MxPrint& pf) const;
virtual std::vector<std::string> pictureStrings(const MxPrint& pf) const;

virtual void print(std::ostream&   ostr, const MxPrint& pf) const; 
virtual void picture(std::ostream& ostr, const MxPrint& pf) const; 
virtual void write(std::ofstream &fp, int form=0) const;
virtual void read(std::ifstream &fp);
virtual void readASC(std::istream& istr);
virtual void ask( );

// ____________________________________________________________________________
// M                 CLASS D_MATRIX MISCELLANEOUS FUNCTIONS
// ____________________________________________________________________________
 
 
virtual void resize(int i, int j);
     
        // Input          dmx   : A d_matrix (this)
        //                i,j   : Row, column dimensions
        // Output         none  : Size of dmx readjusted to ixj
        // Note                 : For dmx, failure if i!=j
        // Note                 : Destroys data if not the same size,
        //                        no alterations if same size
 

virtual _matrix* copy();
 
        // Input                dmx  : Diagonal matrix (this)
        // Output               mx   : Pointer to diagonal matrix
        // Note                      : This returns a new d_matrix on exit.
        //                             Unlike class matrix where assignment
        //                             involves only referencing, a new copy
        //                             of the matrix data is done here


virtual void convert(_matrix* mx);

        // Input          dmx   : A d_matrix (this)
        //                mx    : A matrix of any type
        // Output         mx    : Elements of dmx are placed into mx
        //                        in a fashion so as to maintain the mx
        //                        matrix type.  The result is the "best"
        //                        matrix conversion of dmx will reside in mx
        // Note                 : Because this function uses direct element
        //                        access in assignment (e.g. mx(i,j) = z) it
        //                        has the ability to really mess up inherent
        //                        matrix structure (e.g. create a Hermitian
        //                        matrix with complex diagonals). So be careful
        //                        in adding code & test this function thoroughly!
        // Note                 : Since a diagonal matrix cannot be converted to
        //                        an identitiy matrix unless all of its elements are
        //                        1, this conversion results in loss of data.
        // Note                 : Any existing elements in mx will vanish
 

//
// ------------------------ New Conversion Routines ---------------------------
//                  (Been Having Troubles With The Old One)


virtual i_matrix* IMX();

        // Input                dmx     : A diagonal matrix (this)
        // Output               imx	: A new identity matrix whose
        //                                elements are the "same" as dmx
        // Note                         : New imx will have NO dmx elements
        // Note                         : This will allocate memory


virtual d_matrix* DMX();

        // Input                dmx     : A diagonal matrix (this)
        // Output               dmx     : A new d_matrix matrix whose
        //                                elements are the same as dmx
        // Note                         : This will NOT allocate memory
	//				: Return is just dmx


virtual h_matrix* HMX();

        // Input                dmx     : A diagonal matrix (this)
        // Output               nmx     : A new h_matrix matrix whose
        //                                elements are the same as dmx
        // Note                         : This will allocate memory
 
 
virtual n_matrix* NMX();
 
        // Input                dmx     : A diagonal matrix (this)
        // Output               nmx     : A new n_matrix matrix whose
        //                                elements are the same as dmx
        // Note                         : This will allocate memory
};
#endif					// d_matrix.cc
