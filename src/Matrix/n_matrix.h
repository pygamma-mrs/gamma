/* n_matrix.h ***************************************************-*-c++-*-
**                                                                      **
**                             G A M M A                                **
**                                                                      **
**      General Matrix                            Interface             **
**                                                                      **
**      Copyright (c) 1990, 1997                                        **
**      Tilo Levante, Scott A. Smith                                    **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zuerich / Switzerland                                      **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The class n_matrix defines normal complex matrices for C++ with the  **
** usual operations +, -, *, / as well as in/output routines and the    **
** more commonly used matrix functions.                                 **
**                                                                      **
*************************************************************************/

#ifndef   Gn_matrix_h_			// Is the file already included?
#  define Gn_matrix_h_ 1		// If no, then include it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <Matrix/_matrix.h>		  // Include base matrix class
#include <Matrix/d_matrix.h>		  // Include diagonal arrays
#include <Matrix/i_matrix.h>		  // Include identity arrays
#include <vector>		  	  // Include libstdc++ STL vectors		

#ifndef PI				// Define PI if not done
#define PI 3.14159265358979323846	// Needed for FFT routine
#endif

class n_matrix : public _matrix {
	
complex* data;				  // Pointer to complex data array
bool     unitary;			  // Flag if unitary or not

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                      CLASS N_MATRIX ERROR HANDLING
// ____________________________________________________________________________
 
/*      Input                   nmx     : A complex matrix (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : string in message
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

         void NMxerror(int idx,                        int nr=0) const; 
         void NMxerror(int idx, const std::string& PN, int nr=0) const;
volatile void NMxfatal(int idx)                                  const;
volatile void NMxfatal(int idx, const std::string& PN)           const;

// ____________________________________________________________________________
// ii                        CLASS N_MATRIX CHECKING
// ____________________________________________________________________________
 
bool CheckDims(_matrix* mx,          int warn=1);
bool CheckDim(_matrix* mx,  int mul, int warn=1);

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

  friend class d_matrix;
  friend class i_matrix;
  friend class h_matrix;

// ____________________________________________________________________________
// A               CLASS N_MATRIX CONSTRUCTORS AND DESTRUCTOR
// ____________________________________________________________________________

/*  Arguments                      Constructed Array                            
    ---------           -------------------------------------------------------
        -                An empty array
     nr, nc              An (nr x nc) array
    nr, nc, z            An (nr x nc) array with all <i|nmx|j> = z
       nmx               A duplicate of complex array nmx
 
   Note that all constructors invariably call the base class constructor
   (i.e. the one in _matrix) which sets the row & column sizes.  The added
   code should simply allocate the data array of proper size. Similarly, the
   destructor only destroys the data array, anything else is destroyed in
   the base class.                                                           */
 
         n_matrix( );
         n_matrix(int i, int j);
         n_matrix(int i, int j, const complex& z);
         n_matrix(const n_matrix& m);
virtual  ~n_matrix();

// ____________________________________________________________________________
// B               CLASS N_MATRIX ACCESS AND ASSIGNMENT
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
 

void             operator= (const n_matrix& mx);
virtual complex& operator() (int i, int j);
virtual complex  get(int i, int j) const;
virtual bool     put(const complex& z, int i, int j);
virtual bool     put_h(const complex& z, int i, int j);
virtual _matrix* get_block(int row, int col, int nrows, int ncols);
virtual bool     put_block(int row, int col, _matrix* mx);

// ____________________________________________________________________________
// C               CLASS N_MATRIX HERMITIAN_TYPE HANDLING
// ____________________________________________________________________________

/* These functions handle the checking of whether an array is of a Hermitian
   type.  Since GAMMA has a Hermitian matrix type (h_matrix) it is easy to
   know the whether the array is stored Hermitian or not. The other funtion
   will test if the array is Hermitian to with d, where the test performed is
   norm(<i|mx|j> - <j|mx|i>*) < d for all elements. These functions are a bit
   outdated now that GAMMA has a Hermitian array class and will probably be
   removed in favor of a TF is_hermitian function function or something.  More
   appropriate would be an stored_unitary function.....
 
     Function            Output                       Description
   ----------------  --------------  ------------------------------------------
   stored_hermitian  hermitian_type  Returns _hermitian=FALSE
   test_hermitian    hermitian_type  Returns _hermitian=T/F (within d)       */
 
// Note:  Functions set_hermitian & check_hermitian are found in class matrix!

virtual hermitian_type stored_hermitian( ) const; 
virtual hermitian_type test_hermitian(double d = GMxCut) const;


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

// Note: Functions set_type & check_type are in class matrix, not here.

virtual matrix_type stored_type( ) const;
virtual matrix_type test_type(matrix_type m=n_matrix_type,double d=GMxCut) const;
virtual std::string mxtype()        const;
virtual std::string mxtype(bool pf) const;

// ____________________________________________________________________________
// E                 CLASS N_MATRIX VARIOUS MATRIX CHECKS
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

Note that the unitary check herein tests if the array inverse is equal to its
adjoint.  So, the routine multiplies the matrix by its adjoint and then looks
to see how well the result conforms to an identity matrix.  This can use up
lots of time....
                            ---
                            \                   *t 
                   z(i,j) = /    <i|nmx|k><k|nmx  |j> 
                            --- 
                             k                                              */

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
// F              CLASS N_MATRIX BINARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________

 
/* These functions allow for simple mathematical operations between two
   matrices and between a matrix and a scalar.

   Operator Arguments      Result        Operator Arguments       Result
   -------- --------- -----------------  -------- --------- -------------------
      +      mx1,mx2      mx1+mx2           *       mx,z           z*mx
      -      mx1,mx2      mx1-mx2           *       z,mx           z*mx
      *      mx1,mx2      mx1*mx2           *       mx,d           d*mx
      /      mx1,mx2      mx1*inv(mx2)      *       d,mx           d*mx
      /      mx,z         (1/z)*mx          *       mx,d         (1/d)*mx
      /      z,mx         z*inv(mx)         *       d,mx         d*inv(mx)  */


virtual _matrix* add(_matrix* mx);
 
        // Input                nmx  : Input normal matrix (this)
        //                       mx  : Second matrix
        // Output               mx1  : New matrix which is the addition of
        //                             the two input arrays; mx1 = nmx + mx
 

virtual _matrix* subtract(_matrix* mx);

        // Input                nmx  : Input normal matrix (this)
        //                       mx  : Second matrix
        // Output               mx1  : New matrix which is the subtraction of
        //                             the two input arrays; mx1 = nmx - mx


virtual _matrix* multiply(_matrix* mx);
 
        // Input                nmx  : Input normal matrix (this)
        //                       mx  : Second matrix
        // Output               pdt  : New matrix which is the product nmx*mx
        // Note                      : This uses the internal structure of
        //                             other _matrix derived classes!  Thus if
        //                             any implementation changes so must this!
 
//                            ---
//                <i|pdt|j> = \   <i|nmx|k> <k|mx|j>
//                            /
//                            ---
//   


virtual _matrix* multiply(const complex& z);

        // Input                nmx  : Input normal matrix (this)
        //                      z    : Complex number
        // Output               pdt  : New matrix which is the product z*nmx


virtual _matrix* divide(_matrix* mx);
 
        // Input                nmx  : Input normal matrix (this)
        //                       mx  : Second matrix                         -1
        // Output               pdt  : New matrix which is the product nmx*mx
        // Note                      : Still under construction
        // Note                      : This uses the function inv
 
//                            ---                -1
//                <i|pdt|j> = \   <i|nmx|k> <k|mx  |j>
//                            /
//                            ---
//  


virtual _matrix* divide(const complex& z);
 
        // Input                nmx  : Input normal matrix (this)
        //                      z    : Complex number
        // Output               pdt  : New matrix which is the product (1/z)*nmx
        // Note                      : This handles divisions in code such as mx/z
 

// ____________________________________________________________________________
// G                CLASS N_MATRIX UNARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________
 
/* These functions allow for simple mathematical operations between two
   matrices and between a matrix and a scalar.  These are member functions so
   below I use mx to represent *this.

   Operator Arguments      Result        Operator Arguments     Result
   -------- --------- -----------------  -------- --------- ----------------
      +=     mx,mx1       mx+=mx1           *=      mx,z         mx*=z
      -=     mx,mx1       mx-=mx1           *=      mx,d         mx*=d
      *=     mx,mx1       mx*=mx1           /=      mx,z         mx*=(1/z)
      /=     mx,mx1      mx*=inv(mx1)       /=      mx,d         mx*=(1/d)
      /      mx,z         (1/z)*mx          *       mx,d         (1/d)*mx
      /      z,mx         z*inv(mx1)        *       d,mx         d*inv(mx)  */

virtual _matrix* add_two(_matrix* mx);
virtual _matrix* subtract_two(_matrix* mx);
virtual _matrix* multiply_two(_matrix* mx);
virtual _matrix* multiply_two(const complex &z);
virtual _matrix* divide_two(_matrix* mx);
virtual _matrix* divide_two(const complex &z);

// ____________________________________________________________________________
// H                CLASS N_MATRIX SIMPLE UNARY FUNCTIONS
// ____________________________________________________________________________
   
/* These functions perform simple simple mathematical operations on a complex
   matrix.

    Function          Result             Function           Result           
   ---------- -----------------------    --------  ------------------------
   negate     <i|nmx|j> -> -<i|nmx|j>       RE      <i|nmx|j>->Re(<i|nmx|j>)
   conjugate  <i|nmx|j> -> <i|nmx*|j>       IM      <i|nmx|j>->Im(<i|nmx|j>)
   transpose  <i|nmx|j> -> <j|nmx|i>      adjoint   <i|nmx|j>-><j|nmx*|i>   */

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

virtual _matrix* swaprows(int i, int j);
virtual _matrix* swapcols(int i, int j);
virtual _matrix* permute( int i, int j);
virtual double   maxRe() const;
virtual double   maxIm() const;
virtual complex  maxZ()  const;
virtual double   minRe() const;
virtual double   minIm() const;
virtual complex  minZ()  const;

// ____________________________________________________________________________
// I                 CLASS N_MATRIX SIMPLE BINARY FUNCTIONS
// ____________________________________________________________________________


virtual complex trace(_matrix* mx);

	// Input            nmx : An n_matrix (this)
	//		     mx : A second matrix
	// Output            z  : A complex value wish is the trace of
	//                        the product of two input arrays
	//                         z = Tr{nmx*mx} = Tr{mx*nmx}
	// Note			: This is faster than taking the matrix	
	//			  product and then using unary trace!
 

virtual _matrix* adjoint_times(_matrix* mx);

	// Input            nmx : An n_matrix (this)
	//		     mx : A second matrix
	// Output           mx1 : A matrix which is product of the
	//                        adjoint of n_matrix and mx
	//				      T *
	//                         mx1 = [(nmx ) ] * mx
	// Note			: This is faster than taking the adjoint
	//			  of nmx and then the product!  Use if the
	//			  adjoint of nmx is not needed. 


virtual _matrix* times_adjoint(_matrix* mx);

	// Input            nmx : An n_matrix (this)
	//		     mx : A second matrix
	// Output           mx1 : A matrix which is product of the
	//                        n_matrix and the adjoint of mx
	//					  T *
	//                         mx1 = nmx* [(mx ) ]
	// Note			: This is faster than taking the adjoint
	//			  of mx and then the product!  Use if the
	//			  adjoint of mx is not needed. 


// ____________________________________________________________________________
// J                 CLASS N_MATRIX COMPLEX UNARY FUNCTIONS
// ____________________________________________________________________________

/*  Function/Operator                         Result
    -----------------   ------------------------------------------------------
         det            Determinant of the array (NOT IMPLEMENTED)
         rank           Rank of the array (NOT IMPLEMENTED)                 */

virtual complex det();
virtual int     rank();
  

// ********************** Discrete Fourier Transform ******************** 

void FFT(int isign, bool comp);

	// Input          *this : An n_matrix
        //                isign : Flag for transform or inverse transform
	//		                  1=FFT, -1=IFFT
	//		   comp	: Flag for FFT compatiblity with Brigham
	//			  and MATLAB (or original GAMMA)
        // Output          this : Modifies this to contain either its
	//			  discrete Fourier transform, or its
	//			  inverse discrete Fourier transform

// ____________________________________________________________________________
// K         CLASS N_MATRIX COMPLEX BINARY (&TRINARY) FUNCTIONS
// ____________________________________________________________________________


// ************************** tensor product ****************************

virtual _matrix* tensor_product(_matrix* mx);

        // Input            nmx : An n_matrix (this)
        //                   mx : A second matrix
        // Output           pdt : A matrix which is the tensor product
        //                        of the two input matrices

        //                           pdt        =   nmx (x) mx

        //                       (m*o x n*p)       (mxn)   (oxp)

        //                    <i*o+k|pdt|j*p+l> = <i|nmx|j><k|mx|l>


// ____________________________________________________________________________
// L                       CLASS N_MATRIX I/O FUNCTIONS
// ____________________________________________________________________________

/* These functions perform both ASCII and Binary input/output operations on a
   complex matrix.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   print      ostream     Writes array elemnts in ASCII to output stream
   write      ofstream    Writes array elemnts in Binary to output file stream
   read       ifstream    Read in array from Binary input file stream
   readASC    istream     Read in array from an input stream
   ask          ----      Interactively requests complex matrix from user
   
   The binary format of n_matrix has data only, not size information.  The size
   is taken care of in classes matrix/_matrix and should done prior to use of
   these functions.  The data ordering is Re(<i|hmx|j>, Im(<i|hmx|j> columns
   then rows (i.e. row by row.)                                             */

virtual std::vector<std::string> printStrings(const   MxPrint& pf) const;
virtual std::vector<std::string> pictureStrings(const MxPrint& pf) const;

virtual void print(std::ostream&   ostr, const MxPrint & PF) const;
virtual void picture(std::ostream& ostr, const MxPrint & PF) const;
virtual        void write(std::ofstream &fp, int form=0);
virtual        void read(std::ifstream &fp);
virtual inline void readASC(std::istream& istr);
virtual        void ask( );

/* Note that the value of form in writing a normal matrix has no consequence
   because all normal matrix elements are output.  This is not the case in
   many other matrix types                                                   */

// ____________________________________________________________________________
// M                CLASS N_MATRIX MISCELLANEOUS FUNCTIONS
// ____________________________________________________________________________

/* These are functions that don't fit into the othere categories.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   resize      nr, nc     Resizes nmx to have nr rows & nc columns, type change
   copy         ----      Produces a copy of the nmx, allocates memory
   convert       mx       Matrix mx has its values set to those of nmx (no mem)
   IMX          ----      Identity matrix produced from nmx conversion (mem)
   DMX          ----      Diagonal matrix produced from nmx conversion (mem)
   HMX          ----      Hermitian matrix produced from nmx conversion (mem)
   NMX          ----      Complex matrix produced from nmx conversion (mem)
   
   The conversion functions will often drop elements, this is unavoidable. For
   example, one cannot generally convert a normal array to a diagonal array.
   There are no elements lost when using HMX or NMX.  The *MX functions
   will always allocate new memory for the converted matrix from hmx.      */ 

virtual void      resize(int i, int j);
virtual _matrix*  copy();
virtual void      convert(_matrix* mx);
virtual i_matrix* IMX();
virtual d_matrix* DMX();
virtual h_matrix* HMX();                                                      
virtual n_matrix* NMX();
 

// ____________________________________________________________________________
// N                CLASS N_MATRIX DIAGONALIZATION FUNCTIONS
// ____________________________________________________________________________

/*                         Diagonalization Routines

  Note That The Two Matrices Which Are Given As Arguments MUST Have Their
  Referencing Set External To This Class. Also, As These Enter The Routine
  As Pointers Which Will Be Set To Point To the Newly Transformed Arrays,
  They Should NOT Point To Any Array (Use No Memory But The Pointer)

           Input        nmx     : A general complex matrix (this)    
                        mxd     : Pointer to mx to become "diagonalized"
                        mxev    : Pointer to mx to become eigenvector array
           Output       void    : The matrices mxd and mxev are set to
                                  the Hermitian matrix of dmx eigenvalues and
                                  the matrix of dmx eigenvectors respectively
           Note                 : For h_matrix, mxd should be real and
                                  mxev both unitary and Hermitian
           Note                 : The reference count to dmx must be twice
                                  incremented external by the call origin.   */
 
virtual std::vector<int> BlockDiag(_matrix*    (&BD), std::vector<int> &U);
virtual void             HermTriDiag(_matrix* (&HTD), _matrix* (&U));
virtual void             SymTriDiag(_matrix*  (&STD), _matrix* (&U));
virtual void             SymDiag(_matrix*       (&D), _matrix* (&U));

virtual void diag(_matrix* (&D),          _matrix* (&U));

        // Input        nmx     : n_matrix (this)
        //              mxd     : Pointer to mx to become diagonal mx
        //              mxev    : Pointer to mx to become eigenvectors mx
        // Output       void    : The matrices mxd and mxev are set to
        //                        the diagonal matrix of nmx eigenvalues and
        //                        the matrix of dmx eigenvectors respectively
        // Note                 : The reference count to nmx must be twice
        //                        incremented external by the call origin.

void     cred  (n_matrix& z);
void     rred  (n_matrix& w, int newU=0);
void     tqli  (n_matrix& z, d_matrix& a);
complex* corth (int low, int igh);
int      comqr3(int low,int igh,complex *ort,d_matrix& w,n_matrix& z,int flag);

// ____________________________________________________________________________
// O                   CLASS N_MATRIX INVERSION FUNCTIONS
// ____________________________________________________________________________

// sosi - these are the old inversion routines and will be deleted eventually

virtual _matrix* xinv();
             int LU_decomp(int* indx);
            void LU_backsub (int* indx, n_matrix& b);

	// Input            nmx : An n_matrix (this)
	// Output            mx : A new matrix which is the inverse of nmx
	//                                 -1
	//                         mx = nmx  ; mx * nmx = I

/* The inverse of a  matrix can be obtained in several ways.  Herein
   we will use the LU decomposition to obtain the inverse. In this
   case we have (using N=Normal Complex, L=lower triang., & U=upper triang.)

                               -1               -1
                          N * N  = I = L * U * N
  
   where we 1st solve L*Y=I for Y, then U*inv(N)=Y for inv(N) taking advantage
   of the special formats of L & U.

  The LU & LU inverse routines used here are derived from "Numerical Recipies",
  W.H.  Press, B.P.Flannery, S.A.Teukolsky, & W.T.Vetterling, Cambridge 
  University Press, 1989.  See pages 31 through 36. The code was extensively
  adapted for C++ and GAMMA arrays.

  The LU algorithm takes as input a normal array A and an integer array indx.
  It returns returns the LU decomposition of N' where N' is N with any needed
  row permutations. The latter (row permutations) are stored in the integer
  array indx (indx[0] < 0 flags no row changes).

  The LU decomposition formulates the equation (depicted for 4x4 N below) 

                                N = LU

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

  The LU inverse algorithm takes an input matrix N' in its LU decomposition
  format where N' is some original matrix N with any needed row permutations
  to attain the LU form. The row permutations used to relate N to N' are stored
  in the integer array indx.  The function then proceeds to solve the following
  equation for |x> given a |b>:
                                 N|x> = |b>

  This is accomplished by considering (reformulating) the problem as

                    N|x> = (LU)|x> = L(U|x>) = L|y> = |b>

  using the LU decomposition. This is first solved for the vector |y>, where
  L|y> = |b>, easily accomplished since L is lower-triangular.  Subsequently
  the equation U|x> = |y> is solved for the desired |x>, again easily done
  since U is upper=triangular.  By repeating this process for many |x> and |b>
  pairs of columns, the inverse of the original array can be obtained.
                                    -1
        N [|xo> |x1> ... |xn>] = N N  =  [|bo> |b1> .... |bn>] = I 


    Function  Arguments               Result                      Notes
  ---------- --------- ----------------------------------- -------------------
     inv        ---    inverse(nmx) where inv(nmx)*nmx = I
     LU        int*    LU == nmx + row permutations        One Array L\U Form
     LUinv    int*, B  Solutionh X where LU*X = NX = B     Inverse if B=I


        N [|xo> |x1> ... |xn>] = N N  =  [|bo> |b1> .... |bn>] = I          */


virtual _matrix* inv();
virtual _matrix* LU(int *indx);
virtual _matrix* LUinv(int *indx, _matrix* LU);

};
#endif								    // n_matrix.h
