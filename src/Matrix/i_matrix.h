/* i_matrix.h ***************************************************-*-c++-*-
**								   	**
**                                G A M M A				**
**						   			**
**	Identity Matrix	                 Interface Definition	        **
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

#ifndef   Gi_matrix_h_			// Is file already included?
#  define Gi_matrix_h_ 1		// If no, then remember that 
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <Matrix/_matrix.h>		 // Include the base matrix
#include <Matrix/complex.h>		 // Include libstdc++ complex
#include <Matrix/n_matrix.h>             // Include normal matrices
#include <Matrix/d_matrix.h>             // Include diagonal matrices
#include <Matrix/h_matrix.h>             // Include hermitian matrices
#include <vector>                        // Include libstdc++ STL vectors 

static complex Zizero(0);		  // Complex public changable 0
static complex Zione(1);		  // Complex public changable 1
	
class i_matrix : public _matrix
  {

// ____________________________________________________________________________
// i                      CLASS I_MATRIX ERROR HANDLING
// ____________________________________________________________________________
 
/*      Input                   imx     : An Identity matrix (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : string in message
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

         void IMxerror(int idx,                       int nr=0) const; 
         void IMxerror(int idx, const std::string& PN,int nr=0) const;
volatile void IMxfatal(int idx)                                 const;
volatile void IMxfatal(int idx, const std::string& PN)          const;

friend void volatile i_matrix_err(int error);

// ____________________________________________________________________________
// ii                      CLASS I_MATRIX CHECKING
// ____________________________________________________________________________

virtual bool CheckDims(_matrix* mx, int warn=1) const;

 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

  friend class d_matrix;
  friend class n_matrix;
  friend class h_matrix;
  
// ____________________________________________________________________________
// A                CLASS I_MATRIX CONSTRUCTORS AND DESTRUCTOR
// ____________________________________________________________________________
 
/*  Arguments                      Constructed Array
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

inline i_matrix( );
       i_matrix(int i, int j);
       i_matrix(int i, int j, const complex& z);
inline i_matrix(const i_matrix& m);
inline virtual ~i_matrix();

// ____________________________________________________________________________
// B                    CLASS I_MATRIX ACCESS AND ASSIGNMENT
// ____________________________________________________________________________

/* There are two flavors of element access operators herein. The one that might
   commonly be used, i.e. mx(i,j), is DIRECTLY accessing the element <i|mx|j>.
   That can save CPU time by avoiding value copying when altering or obtaining
   elements.  However its misuse can lead to trouble in this class because no
   identity matrix elements are stored and there are restrictions on the
   elements so that the array remains identity.  An example of when this is
   trouble would be "imx(1,1)=7;".  This fails because only diagonal
   elements may be set to one and no off-diagonal elements exist.
   
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
          =             Current (Hermitian) array set equal to input hmx
    (int,int)           Reference to element <i|imx|j> (Potential Danger)
    get(int,int)        Copy of element <i|imx|j>      (Safe)
    put(int,int)        Assigns element <i|imx|j>      (Return FALSE if fails)
    put_h(int,int)      Assigns both <i|imx|j> & <j|imx|i> (Return TRUE/FALSE)
    get_block(r,c,R,C)  Returns imx block of size RxC starting at <r|imx|c>
    put_block(r,c,mx)   Places mx into imx at position <r|imx|c> (TRUE/FALSE)

    Again, note that operations which return TF will trigger class matrix to
    take appropriate action when they return FALSE.  This usually means that
    the operation cannot be handled by a identity array so the array type
    is changed to be more generic.  The the operation is reattempted.       */

inline void             operator= (const i_matrix& mx);
virtual inline complex& operator() (int i, int j);
virtual inline complex  get(int i, int j) const;
virtual inline bool     put(const complex& z, int i, int j);
virtual inline bool     put_h(const complex& z, int i, int j);
virtual _matrix*        get_block(int r,int c, int nr, int nc);
virtual inline bool     put_block(int row, int col, _matrix* mx);

// ____________________________________________________________________________
// C                CLASS I_MATRIX HERMITIAN_TYPE HANDLING
// ____________________________________________________________________________

/* These functions handle the checking of whether an array is of a Hermitian
   type.  Since identity matrices are Hermitian, the tests are easy and don't
   actually use the input value d.  These are somewhat outdated now that
   GAMMA has a Hermitian array class and will probably be removed in favor of
   a TF is_hermitian function function or something.  More appropriate would be
   an stored_unitary function.....

     Function            Output                       Description
   ----------------  --------------  ------------------------------------------
   stored_hermitian  hermitian_type  Returns _hermitian=TRUE
   test_hermitian    hermitian_type  Returns _hermitian=TRUE (d is unused)   */


virtual inline hermitian_type stored_hermitian( ) const; 
virtual inline hermitian_type test_hermitian(double d = GMxCut) const;
 
// ____________________________________________________________________________
// D                  CLASS I_MATRIX MATRIX_TYPE HANDLING
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

// Note:  Functions set_type & check_type are found in class matrix, not here.

virtual inline matrix_type stored_type( ) const;
virtual inline matrix_type test_type(
                      matrix_type m = n_matrix_type, double d = GMxCut) const;
virtual std::string mxtype() const;
virtual std::string mxtype(bool pf) const;

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
   is_squre       bool   TF if rows_ == cols                 (TRUE)
   is_equal       bool   TF if ||<i|imx-mx|j>||<d    all i,j

   These checks are consistent with other matrix classes.  The value d is
   not used and no real math is performed execpt in the "is_equal" function */

virtual inline bool is_symmetric(double d = GMxCut) const;
virtual inline bool is_hermitian(double d = GMxCut) const;
virtual inline bool is_unitary(double   d = GMxCut) const;
virtual inline bool is_real(double      d = GMxCut) const;
virtual inline bool is_imaginary(double d = GMxCut) const;
virtual inline bool is_complex(double   d = GMxCut) const;
virtual inline bool is_zero(double      d = GMxCut) const;
virtual inline bool is_diagonal(double  d = GMxCut) const;
virtual        bool is_square( )                    const;

virtual        bool is_equal(_matrix* mx, double d = GMxCut) const;

// ____________________________________________________________________________
// F                CLASS I_MATRIX BINARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________
 
/* These functions allow for simple mathematical operations between two
   matrices and between a matrix and a scalar.
 
   Function Arguments      Result        Function Arguments       Result
   -------- --------- -----------------  -------- --------- -------------------
     add     imx,mx        imx+mx        multiply  imx,mx         imx*mx=mx
   subtract  imx,mx        imx-mx        multiply  imx,z          z*imx
    divide   imx,mx        imx*inv(mx)    divide   imx,z         (1/z)*imx
      /      z,mx         z*inv(mx)         *       d,mx         d*inv(imx)  */
 
virtual        _matrix* add(_matrix* mx);
virtual        _matrix* subtract(_matrix* mx);
virtual        _matrix* multiply(_matrix* mx);
virtual inline _matrix* multiply(const complex& z);
virtual inline _matrix* divide(_matrix* mx);
virtual inline _matrix* divide(const complex& z);

// ____________________________________________________________________________
// G                CLASS I_MATRIX UNARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________

/* These functions allow for simple mathematical operations between two         
   matrices and between a matrix and a scalar.

                    Function     Arguments          Result
                    ------------ --------- -------------------------
                    add_two       imx,mx   mx+imx (mx += this)
                    subtract_two  imx,mx   mx-imx (mx -= this)
                    multiply_two  imx,mx   mx*imx (mx *= this)
                    multiply_two  imx,z    z*imx  (mx = z*this)
                    divide_two    imx,mx   mx*imx (mx *= inv(this))
                    divide_two    imx,z    z*imx  (mx = z*inv(this))        */
 
virtual        _matrix* add_two(_matrix* mx);
virtual        _matrix* subtract_two(_matrix* mx);
virtual        _matrix* multiply_two(_matrix* mx);
virtual        _matrix* multiply_two(const complex &z);
virtual inline _matrix* divide_two(_matrix* mx);
virtual inline _matrix* divide_two(const complex &z);

// ____________________________________________________________________________
// H                  CLASS I_MATRIX SIMPLE UNARY FUNCTIONS
// ____________________________________________________________________________

/* These functions perform simple simple mathematical operations on an Identity
   matrix. Note that for Identity arrays many of these do nothing.
 
    Function          Result             Function           Result
   ---------- -----------------------    --------  ------------------------
   negate     <i|imx|j> -> -<i|imx|j>       RE      <i|imx|j>->Re(<i|imx|j>)
   conjugate  <i|imx|j> -> <i|imx*|j>       IM      <i|imx|j>->Im(<i|imx|j>)
   transpose  <i|imx|j> -> <j|imx|i>      adjoint   <i|imx|j>-><j|imx*|i>   */


virtual inline _matrix* negate();		// Return is -1 d_matrix
virtual inline _matrix* RE();			// Return is this
virtual inline _matrix* IM();			// Return is 0 d_matrix
virtual        _matrix* conjugate();		// Return is this
virtual inline _matrix* transpose();		// Return is this
virtual inline _matrix* adjoint();		// Return is this
virtual inline _matrix* mxexp();		// Return is this
virtual inline complex  trace();		// Return is size

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

virtual _matrix* swaprows(int i, int j);
virtual _matrix* swapcols(int i, int j);
virtual _matrix* permute( int i, int j);
virtual double   maxRe()   const;
virtual double   maxIm()   const;
virtual complex  maxZ()    const;
virtual double   minRe()   const;
virtual double   minIm()   const;
virtual complex  minZ()    const;

// ____________________________________________________________________________
// I                CLASS I_MATRIX SIMPLE BINARY FUNCTIONS
// ____________________________________________________________________________

        // Input            imx : An i_matrix (this)
        //                   mx : A second matrix
	// Outpue (trace)    z  : A complex value wish is the trace of
        //                        the product of two input arrays
        //                         z = Tr{imx*mx} = Tr{mx}
        // Output (adj_tim) mx1 : A matrix which is just mx as the
	//			  product of adjoint imx and mx is mx
        //                                    T *
        //                         mx1 = [(imx ) ] * mx = mx
        // Output (tim_adj) mx1 : A matrix which is the adjoing of mx
        //                                        T *      T*
        //                         mx1 = nmx* [(mx ) ] = mx
        // Note                 : imx.trace(mx) is faster than taking
	//			  trace(imx*mx), but not much.

virtual inline complex trace(_matrix *mx);
virtual       _matrix* adjoint_times(_matrix* mx);
virtual       _matrix* times_adjoint(_matrix* mx);

// ____________________________________________________________________________
// J              CLASS I_MATRIX COMPLEX UNARY FUNCTIONS
// ____________________________________________________________________________

// Note:  The alogorithm for FFT is not found here.  The identity matrix is
//        converted to a normal matrix and it's transformation occurs in
//	  class n_matrix.
 
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

virtual std::vector<int> BlockDiag(_matrix*    (&BD), std::vector<int> &U);
virtual void             HermTriDiag(_matrix* (&HTD), _matrix* (&U));
virtual void             SymTriDiag(_matrix*  (&STD), _matrix* (&U));
virtual void             SymDiag(_matrix*       (&D), _matrix* (&U));
virtual void             diag(_matrix*          (&D), _matrix* (&U));
 
// ____________________________________________________________________________
// O                    CLASS I_MATRIX INVERSION FUNCTIONS
// ____________________________________________________________________________

/* These functions parallel those found in other GAMMA matrix classes.  They
   do not need to do anything for identity matrices.

    Function  Arguments               Result                      Notes
  ---------- --------- ----------------------------------- --------------------
     inv        ---    inverse(imx) where inv(imx)*imx = I Inverse is this
     LU        int*    LU == imx + row permutations        LU is this
     LUinv    int*, A  Solution X where LU*X = AX = I      Inverse A if no indx

  Note that the LU function takes the LU decomposition of any square array, 
  and that is not likely to be an identiy array.                             */

virtual inline _matrix* inv();
virtual inline _matrix* LU(int* indx);
virtual        _matrix* LUinv(int *indx, _matrix* LU);

virtual inline complex det();

        // Input            imx : An i_matrix (this)
        // Output 	      z : Value the determinant; det{imx}
	// Note			: For i_matrix, z=1

virtual inline int rank();

        // Input            imx : An i_matrix (this)
        // Output 	      i : Rank of matrix imx
	// Note			: For i_matrix, rank = dim
 
// ____________________________________________________________________________
// K                CLASS I_MATRIX COMPLEX BINARY FUNCTIONS
// ____________________________________________________________________________

// ************************** tensor product ****************************

virtual _matrix* tensor_product(_matrix* mx);
 
        // Input            imx : An i_matrix (this)
        //                   mx : A second matrix
        // Output           pdt : A matrix which is the tensor product
        //                        of the two input matrices
 
        //                           pdt        =   imx (x) mx
 
        //                       (m*n x n*o)       (mxn)   (nxo)
 
        //                    <i*n+k|pdt|j*o+l> = <i|imx|j><k|mx|l>


// ____________________________________________________________________________
// L                      CLASS I_MATRIX I/O FUNCTIONS
// ____________________________________________________________________________

/* These functions perform both ASCII and Binary input/output operations on an
   Identity matrix.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   print      ostream     Writes array elemnts in ASCII to output stream
   write      ofstream    Writes array elemnts in Binary to output file stream
   read       ifstream    Read in array from Binary input file stream
   readASC    istream     Read in array from an input stream
   ask          ----      Interactively requests identity matrix from user
   
   The binary format of i_matrix has data, but no size information.  The size
   is taken care of in classes matrix/_matrix and should done prior to use of
   these functions.  The data ordering is Re(<i|hmx|j>, Im(<i|hmx|j> columns
   then rows (i.e. row by row.)                                             */

virtual std::vector<std::string> printStrings(const   MxPrint& pf) const;
virtual std::vector<std::string> pictureStrings(const MxPrint& pf) const;

virtual void print(std::ostream&   ostr, const MxPrint& pf) const;
virtual void picture(std::ostream& ostr, const MxPrint& pf) const;

virtual inline void write(std::ofstream& fp, int form=0) const;
virtual inline void read(std::ifstream& fp);
virtual inline void readASC(std::istream& istr);

// ______________________________________________________________________
// M              CLASS I_MATRIX MISCELLANEOUS FUNCTIONS
// ______________________________________________________________________


virtual inline void resize(int i, int j);

	// Input          *this : i_matrix
	//		    i,j : rows and cols values
        // Output 	  *this : resized matrix with i rows & j cols


virtual inline _matrix* copy();

	// Input          *this : i_matrix
        // Output 	     mx : A copy of the input i_matrix


virtual void convert(_matrix* mx);

	// Input          *this : i_matrix
	// 		     mx : pointer to a matrix
        // Output 	  *this : converts *this to a new matrix
	//			  having the same type as mx

//
// ------------------------ New Conversion Routines ---------------------------
//                  (Been Having Troubles With The Old One)

        // Input                imx     : A identity matrix (this)
        // Output               imx1    : A new i_matrix matrix
        // Output               dmx     : A new d_matrix matrix
        // Output               hmx     : A new h_matrix matrix
        // Output               nmx     : A new n_matrix matrix
        //                                Return matrix elements are the
	//                                same as were in imx
        // Note                         : This will allocate memory
        //	                        : except in the case of returning imx
	//                                which just returns itself

virtual i_matrix* IMX();
virtual d_matrix* DMX(); 
virtual h_matrix* HMX();
virtual n_matrix* NMX();

  };

#endif							// i_matrix.h
