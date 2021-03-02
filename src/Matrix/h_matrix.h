/* h_matrix.h ***************************************************-*-c++-*-
**								   	**
**                                G A M M A			   	**
**						   			**
**	Hermitian Matrix	            Interface definition   	**
**								   	**
**	Copyright (c) 1993						**
**	Scott Smith							**
**	University of California, Santa Barbara				**
**	Department of Chemistry						**
**	Santa Barbara, CA, 93106, USA					**
**								   	**
**      $Header: $
**								    	**
*************************************************************************/

/*************************************************************************
**								   	**
** Description							   	**
**						   			**
** The class h_matrix defines Hermitian	matrices for C++ with the 	**
** usual operations +, -, *, / and in/output routines.			**
**						   			**
*************************************************************************/

#ifndef   Gh_matrix_h_				// Is file already included
#  define Gh_matrix_h_ 1			// No, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the interface
#  endif

#include <Matrix/_matrix.h>
#include <Matrix/complex.h>
#include <vector>				// Include libstc++ STL vectors

class h_matrix : public _matrix
  {
	
  complex* data;				// Matrix elements
  complex cval;					// Temp for 1 element

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                      CLASS H_MATRIX ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   hmx	: A Hermitian matrix (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : string in message
        Output                  none    : Error message output
                                          Program execution stopped if fatal */
 
         void HMxerror(int idx,                        int nr=0) const;
         void HMxerror(int idx, const std::string& PN, int nr=0) const;
volatile void HMxfatal(int idx)                                  const;
volatile void HMxfatal(int idx, const std::string& PN)           const;

friend void volatile h_matrix_err(int error);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

    friend class n_matrix;
    friend class d_matrix;
    friend class i_matrix;

// ____________________________________________________________________________
// A                 CLASS H_MATRIX CONSTRUCTORS AND DESTRUCTOR
// ____________________________________________________________________________
 
/*  Arguments                      Constructed Array
    ---------           -------------------------------------------------------
        -                An empty array
     nr, nc              An (nr x nc) array, uninitialized
    nr, nc, z            An (nr x nc) array, sets <i|mx|i>=Re(z) & <i<j|mx|j>=z
       hmx               A duplicate of Hermitian array hmx
 
   Note that all constructors invariably call the base class constructor
   (i.e. the one in _matrix) which sets the row & column sizes.  The added
   code should simply allocate the data array of proper size. Similarly, the
   destructor only destroys the data array, anything else is destroyed in
   the base class.                                                           */

h_matrix( );
h_matrix(int i, int j);
h_matrix(int i, int j, const complex& z);
h_matrix(const h_matrix& m);
virtual ~h_matrix();

// ____________________________________________________________________________
// B                    CLASS H_MATRIX ACCESS AND ASSIGNMENT
// ____________________________________________________________________________

/* There are two flavors of element access operators herein. The one that might
   commonly be used, i.e. mx(i,j), is DIRECTLY accessing the element <i|mx|j>.
   That can save CPU time by avoiding value copying when altering or obtaining
   elements.  However its misuse can lead to trouble in this class because not
   all Hermitian matrix elements are stored and there are restrictions on the
   elements so that the array remains Hermitian.  Two examples of when this is
   trouble would be "hmx(0,0)=complexi;" and "hmx(2,0)=7". The former fails as
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

void             operator= (const h_matrix& mx);
virtual complex& operator() (int i, int j);
virtual complex  get(int i, int j) const;
virtual bool     put(const complex& z, int i, int j);
virtual bool     put_h(const complex& z, int i, int j);
virtual _matrix* get_block(int row, int col, int nrows, int ncols);
virtual bool     put_block(int row, int col, _matrix* mx);

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

virtual hermitian_type stored_hermitian( ) const; 
virtual hermitian_type test_hermitian(double d = GMxCut) const;

 
// ____________________________________________________________________________
// D                  CLASS H_MATRIX MATRIX_TYPE HANDLING
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
   mxtype            string          Returns the string "Hermitian"          */

virtual matrix_type stored_type( ) const;
virtual matrix_type test_type(matrix_type m=n_matrix_type,
                                                      double d=GMxCut) const;
virtual std::string mxtype()        const;
virtual std::string mxtype(bool pf) const;

// ____________________________________________________________________________
// E                  CLASS H_MATRIX VARIOUS MATRIX CHECKS
// ____________________________________________________________________________

/*   Function    Output                       Description
   ------------  ------  ------------------------------------------------------
   is_symmetric   bool   True if array is real
   is_hermitian   bool   Always true
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
to see how well the result conforms to an identity matrix.  This can use lots
of time..
            ---                         ---
            \                   *t      \
   z(i,j) = /    <i|hmx|k><k|hmx  |j> = /    <i|hmx|k><k|hmx|j> = del +/- d
            ---                         ---                          ij
             k                           k                                  */

virtual bool is_symmetric(double          d=GMxCut) const;
virtual bool is_hermitian(double          d=GMxCut) const;
virtual bool is_unitary(double            d=GMxCut) const;
virtual bool is_real(double               d=GMxCut) const;
virtual bool is_imaginary(double          d=GMxCut) const;
virtual bool is_complex(double            d=GMxCut) const;
virtual bool is_zero(double               d=GMxCut) const;
virtual bool is_diagonal(double           d=GMxCut) const;
virtual bool is_square()                            const;
virtual bool is_equal(_matrix* mx, double d=GMxCut) const;
bool is_tridiagonal(double d=GMxCut) const;
  
// ____________________________________________________________________________
// F                CLASS H_MATRIX BINARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________

/* These functions allow for simple mathematical operations between two         
   matrices and between a matrix and a scalar.
 
   Function Arguments      Result        Function Arguments       Result
   -------- --------- -----------------  -------- --------- -------------------
     add     hmx,mx        hmx+mx        multiply  hmx,mx         hmx*mx
   subtract  hmx,mx        hmx-mx        multiply  hmx,z          z*hmx
    divide   hmx,mx        hmx*inv(mx)    divide   hmx,z         (1/z)*mx
      /      z,mx         z*inv(mx)         *       d,mx         d*inv(mx)  */

virtual _matrix* add(_matrix* mx);
virtual _matrix* subtract(_matrix* mx);
virtual _matrix* multiply(_matrix* mx);
virtual _matrix* multiply(const complex& z);
virtual _matrix* divide(_matrix* mx);
virtual _matrix* divide(const complex& z);
 
// ____________________________________________________________________________
// G              CLASS H_MATRIX UNARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________
 
/* These functions allow for simple mathematical operations between two         
   matrices and between a matrix and a scalar.
 
                    Function     Arguments          Result 
                    ------------ --------- -------------------------
                    add_two       hmx,mx   mx+hmx (mx += this)
                    subtract_two  hmx,mx   mx-hmx (mx -= this)
                    multiply_two  hmx,mx   mx*hmx (mx *= this)
                    multiply_two  hmx,z    z*hmx  (mx = z*this)
                    divide_two    hmx,mx   mx*hmx (mx *= inv(this))
                    divide_two    hmx,z    z*hmx  (mx = z*inv(this))        */

virtual _matrix* add_two(_matrix* mx);
virtual _matrix* subtract_two(_matrix* mx);
virtual _matrix* multiply_two(_matrix* mx);
virtual _matrix* multiply_two(const complex &z);
virtual _matrix* divide_two(_matrix* mx);
virtual _matrix* divide_two(const complex &z);
 
// ____________________________________________________________________________
// H                CLASS H_MATRIX SIMPLE UNARY FUNCTIONS
// ____________________________________________________________________________

/* These functions perform simple simple mathematical operations on a Hermitian
   matrix. Note that for Hermitian arrays the adjoint function does nothing.

    Function          Result             Function           Result
   ---------- -----------------------    --------  ------------------------
   negate     <i|hmx|j> -> -<i|hmx|j>       RE      <i|hmx|j>->Re(<i|hmx|j>)
   conjugate  <i|hmx|j> -> <i|hmx*|j> 	    IM      <i|hmx|j>->Im(<i|hmx|j>)
   transpose  <i|hmx|j> -> <j|hmx|i>      adjoint   <i|hmx|j>-><j|hmx*|i>   */

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

virtual _matrix* swaprows(  int I, int J);
virtual _matrix* swapcols(  int I, int J);
virtual _matrix* permute(   int I, int J);
        void     permute_ip(int I, int J);
virtual double   maxRe()   const;
virtual double   maxIm()   const;
virtual complex  maxZ()    const;
virtual double   minRe()   const;
virtual double   minIm()   const;
virtual complex  minZ()    const;

// ____________________________________________________________________________
// I                 CLASS H_MATRIX SIMPLE BINARY FUNCTIONS
// ____________________________________________________________________________

virtual complex trace(_matrix* mx);
 
        // Input            hmx : A h_matrix (this)
        //                  mx  : A second matrix
        // Output           z   : The trace of hmx*mx


virtual _matrix* adjoint_times(_matrix* mx);

	// Input          *this : An h_matrix
	//                   mx : matrix
	// Output           mx1 : Returns a pointer to a new matrix
	//                                      adjoint
	//                        mx1  = (*this)        * mx = mx


virtual _matrix* times_adjoint(_matrix* mx);

	// Input	this,mx: two matricies
	// Output              : return this*adjoint(mx)
 
 
// ____________________________________________________________________________
// J                 CLASS H_MATRIX COMPLEX UNARY FUNCTIONS
// ____________________________________________________________________________
 

virtual complex det();

	// Input          *this : An h_matrix
        // Output             z : The determinant of the matrix


virtual int rank();

	// Input          *this : An h_matrix
        // Output             z : The rank of the matrix

 
 
// ______________________________________________________________________
// K         CLASS H_MATRIX COMPLEX BINARY (&TRINARY) FUNCTIONS
// ______________________________________________________________________
 
virtual _matrix* tensor_product(_matrix* mx);

	// Input          *this : An h_matrix
	// 		     mx : A matrix
	// Output               : A matrix of size
	//                        (this->rows()*mx->rows(),
	//                         this->cols()*mx->cols()) with the elements:
	//                        (i*mx->rows()+j,k*mx->cols()+l)=
	//                         this(i,k)*mx(j,l)

// ____________________________________________________________________________
// L                     CLASS H_MATRIX I/O FUNCTIONS
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
   then rows (i.e. row by row.)                                             */

virtual std::vector<std::string> printStrings(const   MxPrint& pf) const;
virtual std::vector<std::string> pictureStrings(const MxPrint& pf) const;

virtual void print(std::ostream&   ostr, const MxPrint& PF) const;
virtual void picture(std::ostream& ostr, const MxPrint& PF) const;
virtual        void write(std::ofstream &fp, int form=0);
virtual        void read(std::ifstream &fp);
virtual inline void readASC(std::istream& istr);
virtual        void ask( );
   
// ____________________________________________________________________________
// M                CLASS H_MATRIX MISCELLANEOUS FUNCTIONS
// ____________________________________________________________________________

/* These are functions that don't fit into the othere categories.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   resize      nr, nc     Resizes hmx to have nr rows & nc columns, type change
   copy         ----      Produces a copy of the hmx, allocates memory
   convert       mx       Matrix mx has its values set to those of hmx (no mem)
   IMX          ----      Identity matrix produced from hmx conversion (mem)
   DMX          ----      Diagonal matrix produced from hmx conversion (mem)
   HMX          ----      Hermitian matrix produced from hmx conversion (mem)
   NMX          ----      Complex matrix produced from hmx conversion (mem)
   
   The conversion functions will often drop elements, this is unavoidable. For
   example, one cannot generally convert a Hermitian array to a diagonal array.
   There are no elements lost when using HMX or NMX.  The *MX functions
   will always allocate new memory for the converted matrix from hmx.      */  

virtual void      resize(int i, int j);
virtual _matrix*  copy();
virtual void      convert(_matrix* mx);
virtual i_matrix* IMX();	// hmx -> imx (elements lost, m used)
virtual d_matrix* DMX();	// hmx -> dmx (off-diags lost, m used)
virtual h_matrix* HMX();	// hmx -> hmx (no loss, memory used)
virtual n_matrix* NMX();	// hmx -> nmx (no loss, memory used)

// ____________________________________________________________________________
// N                CLASS H_MATRIX DIAGONALIZATION FUNCTIONS
// ____________________________________________________________________________

/* Hermitian arrays have a few important properties worth noting here. These
   following points are directly from "Linear Algebra and its Applications" by
   Gilbert Strang, Academic Press Inc., 1980 (using H=hermitian & U=unitary):
   
   1. ALL Hermitian array eigenvalues are real (p. 222) & may be degenerate.
   2. Hermitian arrays may be diagonalized by a UNITARY tranformation (p. 225)
      Accordingly, there exists a complete set of orthonormal eigenvectors.
   3. Unitary arrays, U, have the property that inv(U) = adjoint(U) (p. 223).
   4. Since H is diagonlizable, exp(Ht) is never singular (p 207) and may be
      determined by U*exp(Dt)*adjoint(U) where D is the eigenvalue array of
      H. This works even in the case of degenerate eigenvalues.

   In GAMMA, the Hermitian array diagonalization will take place in three
   distinct steps:
   
           cred           -1  rred           -1   tqli            -1
       H -------> U *HTD*U   ------> U *RTD*U   -------> U * D * U
                   1      1           2      2

   where H=Hermitian, HTD=Hermitian Tridiagonal, RTD=Real Tridiagonal,
   D=Diagonal and U=Unitary.  The functions used to do this are as follows:

      cred	: Takes a Hermitian array (H=*this) & a NULL pointer to what
                  is to become a unitary transformation array U.  Modifies
                  (*this) to Hermitian tri-diagonal & U to unitary using a
                  Householder algorithm.
      rred	: Takes a tri-diagonal Hermitian array (H=*this) & pointer to
                  any previous unitary transformation array U or an identity
                  matrix (U=I).  Modifies (*this) to symmetric tri-diagonal
                  array (still h_matrix) & U to unitary using a Householder
                  algorithm.
      tqli	: Takes a tri-diagonal symmetric array (H=*this) and pointer
                  to any previous unitary transformation array U or an 
                  identity matrix (U=I) as well as a pointer to a NULL array
                  that is to contain the array eigenvalues, D.  Modifies D
                  to real diagonal & U to unitary using a QR decomposition 
                  algorithm.
      diag      : Applies cred, rred, and tqli successively.  It is essential
                  that that the input array (pointers) point to empty NULL
                  arrays as their memory will be allocated in the routines
                  (cred and tqli respectively for U & D.) 

  Note That The Returned Diagonal & Eigenvectors Arrays MUST Have Their
  Referencing Set External To This Class. Also, As These Enter The Routine
  As Pointers Which Will Be Set To Point To the Newly Transformed Arrays,
  They Should NOT Point To Any Array (Use No Memory But The Pointer)         */

        void             cred(_matrix* (&U));
        void             rred(_matrix* (&U),                int newU=0);
        void             tqli(_matrix* (&U), _matrix* (&D), int newU=0);
virtual void             diag(_matrix* (&D), _matrix* (&U));
virtual std::vector<int> BlockDiag(_matrix*    (&BD), std::vector<int> &U);
virtual void             HermTriDiag(_matrix* (&HTD), _matrix* (&U));
virtual void             SymTriDiag(_matrix*  (&STD), _matrix* (&U));
virtual void             SymDiag(_matrix*       (&D), _matrix* (&U));


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

virtual _matrix* inv();
virtual _matrix* LU(int *indx);
virtual _matrix* LUinv(int *indx, _matrix* LU);

}; 
#endif 								// h_matrix.h
