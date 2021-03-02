/* matrix.h *****************************************************-*-c++-*-
**									**
**                                G A M M A				**
**									**
**	General Matrix                                  Interface	**
**									**
**	Copyright (c) 1990, 1991, 1992					**
**	Tilo Levante, Scott A. Smith					**
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
** The class matrix defines matrices for GAMMA in C++ with the usual    **
** operations +, -, *, /, and I/O routines.		                **
**                                                                      **
** It uses the following special matrix classes.                        **
**                                                                      **
**        _matrix: dummy matrix                                         **
**       n_matrix: normal matrix                                        **
**	 d_matrix: diagonal matrix                                      **
**       i_matrix: identity matrix                                      **
**       h_matrix: Hermitian matrix                                     **
**                                                                      **
*************************************************************************/

#ifndef   Gmatrix_h_			// Is file already included?
#  define Gmatrix_h_ 1			// If no, then remember it
#  if defined GAMPRAGMA			// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <Matrix/_matrix.h>		// Include base matrix class

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <string>			// Include libstdc++ strings
#include <iostream>			// Include C++ input/output
#include <fstream>			// Include file streams
#include <vector>			// Include libstdc++ STL vectors

MSVCDLL void enable_blockdiag();
MSVCDLL void disable_blockdiag();


class matrix
  {
  _matrix *m;				// Just pointer to stored matrix
  static bool    FFTComp;               // FFT compatibility (Brigham) 
  static bool    BlkDiag;		// Flag to block during diagonalize
  static MxPrint PrntFlgs;		// ASCII output controls

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                    CLASS MATRIX ERROR HANDLING
// ____________________________________________________________________________

/*	Input			mx      : A matrix (this)
                            eidx    : Error index
                            noret   : Flag for linefeed (0=linefeed)
                            pname   : string in message
	Output			none    : Error message output
        				  Program execution stopped if fatal     */

void Mxerror(int eidx, int noret=0) const;
void Mxerror(int eidx, const std::string& pname, int noret=0) const;
volatile void Mxfatality(int eidx) const;
volatile void Mxfatality(int eidx, const std::string& pname) const;

// ____________________________________________________________________________
// ii             CLASS MATRIX VIRTUAL MATRIX HANDLING
// ____________________________________________________________________________

/******************************************************************************
**                                                                           **
** Class matrix tries to reduce the copying of matrices if possible.         **
** In many instances only "virtual copies" of a matrix need be made          **
** rather than a physical (in memory) copy.  A virtual matrix is just        **
** a reference to a matrix (in memory). Production of a virtual copy         **
** then simply increments a counter of the number of references to the       **
** matrix data and no actual data copy is made.                              **
**                                                                           **
** Note that this class (matrix) handles virutal versus real copying. The    **
** matrix type classes (_matrix derived: i_matrix, d_matrix, h_matirx,....)  **
** do not deal with such issues, expecting this generic class to handle it.  **
** With that in mind one can avoid the possibility of having memory leaks.   **
**                                                                           **
******************************************************************************/

friend _matrix* virtual_copy(_matrix* mx);
friend void     virtual_delete(_matrix *mx);
friend _matrix* virtual_to_real_copy(_matrix *mx);

matrix(_matrix* m);

	// Input           m  : Pointer to _matrix
	// Output          mx : A matrix is constructed from m
	// Note                 This is not allowed outside the class


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:


// ____________________________________________________________________________
// A                  CLASS MATRIX CONSTRUCTORS AND DESTRUCTOR
// ____________________________________________________________________________
 
/*    Arguments                         Constructed Array
    -------------       -------------------------------------------------------
         -               An empty array
     nr,nc,mt,ht         An (nr x nc) array of type mt & hermitian type ht
    nr, nc, z, mt        An (nr x nc) array with all <i|mx|j> = z of type mt
         mx              A duplicate of array mx

   The posible values of nt are currently:

	            n_matrix_type == normal_matrix  (default)
	            h_matrix_type == Hermitian matrix
	            d_matrix_type == diagonal matrix
	            i_matrix_type == identity matrix

   The possible values of hermitian_type are

		_hermitian == matrix is Hermitian
		non_hermitian == matrix is not Hermitian (default)           */

MSVCDLC matrix( );
MSVCDLC matrix(int i);

MSVCDLC matrix(int i, int j);
MSVCDLC matrix(int i, int j,                   matrix_type t);
MSVCDLC matrix(int i, int j,                   matrix_type t, hermitian_type h);

MSVCDLC matrix(int i, int j, const complex& z);
MSVCDLC matrix(int i, int j, const complex& z, matrix_type t);
MSVCDLC matrix(int i, int j, const complex& z, matrix_type t, hermitian_type h);

MSVCDLC matrix(int i, int j, double d);
MSVCDLC matrix(int i, int j, double d,         matrix_type t);
MSVCDLC matrix(int i, int j, double d,         matrix_type t, hermitian_type h);

MSVCDLC matrix(const matrix& mx);
MSVCDLC ~matrix ();
 
// ____________________________________________________________________________
// B                    CLASS MATRIX ACCESS AND ASSIGNMENT
// ____________________________________________________________________________

/* There are two flavors of element access operators herein. The one that might
   commonly be used, i.e. mx(i,j), CANNOT directly access the element <i|mx|j>.
   To be absolutely clear: mx(i,j) is NOT a reference to the elemment <i|mx|j>!
   That is because many matrix types do not have all elements stored in memory,
   and one cannot get a reference to something that simply does not exist. Thus  
   having direct element access can lead to big time troubles.  As an example,
   consider what might happen with mx(2,3) = 7; when mx is (internally) stored
   as a diagonal array.  Although <2|mx|3> is 0 in such a case, since it is NOT
   in memory it cannot be referred to nor be set to any value. So, what must
   occur is something more complicated - namely the array can no longer be kept
   diagonal and memory allocated to handle non-diagonal elements.

   The best way to handle this?  There isn't one AFIK. What GAMMA will do is
   allow access to the constant element reference (get it but not set it) so
   it is safe. The "get" function returns copies of the element & the "put"
   function checks that the element is being properly set, so these two are
   always safe as well. Finally, for the user who needs to access the elements
   for speed (hopefully aware of what they are doing......) the elem function
   is provided.  A final WARNING: elem is a dangerous function and can not be
   used without some user thought processes taking part.  The code will always
   compile but the results may not be what one might blindly think.

    Function/Operator                         Result
    -----------------   ------------------------------------------------------
          =             Current array set equal to input mx
    (int,int)           Reference to element <i|mx|j>
    get(int,int)        Copy of element <i|mx|j>
    put(int,int)        Assigns element <i|mx|j>
    put_h(int,int)      Assigns both <i|mx|j> & <j|mx|i>
    get_block(r,c,R,C)  Returns mx block of size RxC starting at <r|hmx|c>
    put_block(r,c,mx1)  Places mx1 into mx at position <r|hmx|c>

    Note that the assignment operator here only makes a virtual copy of the
    array. The data is copied only if a change is made to the matrix. Block
    access (put/get_block) begins at element <r|mx|c> and ends with element
    will be <r+R|mx|c+C> where desired put/get array is R x C               */

MSVCDLL       matrix&  operator= (const matrix& mx);
MSVCDLL const complex& operator() (int i, int j) const;
MSVCDLL       complex& operator() (int i, int j);
MSVCDLL       complex& elem(int  i, int j);
MSVCDLL       complex  get(int   i, int j) const;
MSVCDLL       double   getRe(int i, int j) const;
MSVCDLL       double   getIm(int i, int j) const; 
MSVCDLL       void     put(const   complex& z, int i, int j);
MSVCDLL       void     put_h(const complex& z, int i, int j);

MSVCDLL       matrix   get_block(int row, int col, int nrows, int ncols);
MSVCDLL       void     put_block(int row, int col, const matrix& mx);

// ____________________________________________________________________________
// C                CLASS MATRIX HERMITIAN_TYPE HANDLING
// ____________________________________________________________________________

/*   Function      Output                       Description
 ----------------  ------  -------------------------------------------------
 stored_hermitian  h_type  Returns stored hermitian type
 check_hermitian   h_type  Returns stored hermitian type (with mx(ij)<d = 0)
 set_hermitian     void    Sets matrix hermitian, may result in data loss
 test_hermitian    h_type  Returns possible hermitian type (w/ mx(ij)<d = 0) */

MSVCDLL hermitian_type stored_hermitian( ) const;
MSVCDLL hermitian_type check_hermitian(double       d = GMxCut);
MSVCDLL void           set_hermitian(hermitian_type h = _hermitian);
MSVCDLL hermitian_type test_hermitian(double        d = GMxCut) const;

// ____________________________________________________________________________
// D                   CLASS MATRIX MATRIX_TYPE HANDLING
// ____________________________________________________________________________

/*   Function    Output                       Description
   -----------  -------  ------------------------------------------------------
   stored_type  mx_type  Returns current matrix type
   test_type    h_type   Returns possible hermitian type (w/ mx(ij)<d = 0)
   set_type     void     Set matrix to specified type, may result in data loss
   check_type   mx_type  Test if mx could be matrix type (with mx(ij)<d = 0)
   mxtype	  string   Return a string for the current matrix type       

    Note: test_type will return tested type if true, its own type if false
    Note: check_type will convert the arrray if no data loss will occur.
          If it is not possible, the current matrix type is returned.         */

MSVCDLL matrix_type stored_type( ) const;
MSVCDLL matrix_type test_type(matrix_type t, double d=GMxCut) const;
MSVCDLL void        set_type(matrix_type t);
MSVCDLL matrix_type check_type(const matrix_type t, const double d = GMxCut);
MSVCDLL std::string mxtype() const;

// ____________________________________________________________________________
// E             CLASS MATRIX VARIOUS MATRIX CHECKS & PARAMETERS
// ____________________________________________________________________________
 
/*  Function  Output                       Description
    --------  ------  ---------------------------------------------------------
      rows     int    Returns # of matrix rows       - handled by _matrix
      cols     int    Returns # of matrix columns    - handled by _matrix
      refs     int    Returns # of matrix references - handled by _matrix
      pts      int    Returns # of matrix elements   - handled by _matrix    */
 
MSVCDLL int cols( ) const;
MSVCDLL int rows( ) const;
MSVCDLL int refs( ) const;
MSVCDLL int pts( )  const;
 
/*   Function    Output                       Description
   ------------  ------  ------------------------------------------------------
   is_symmetric   bool    TF if <i|mx|j>-<j|mx|i>  < d       (def. d GMxCut)
   is_hermitian   bool    TF if <i|mx|j>-<j|mx|i>* < d       (def. d GMxCut)
   is_unitary     bool    TF if inv(mx) == adjoint mx, CPU intensive
   is_real        bool    TF if Im(<i|mx|j>) < d for all i,j (def. d GMxCut)
   is_imaginary   bool    TF if Re(<i|mx|j>) < d for all i,j (def. d GMxCut)
   is_complex     bool    TF if is_real && is_imaginary      (def. d GMxCut)
   is_zero        bool    TF if ||<i|mx|j>|| < d for all i,j (def. d GMxCut) 
   is_diagonal    bool    TF if ||<i|mx|j>|| < d for all i!=j(TRUE,d GMxCut)
   is_square      bool    TF if rows_ == cols_ 

Note that the unitary check herein tests if the array inverse is equal to its
adjoint.  So, the routine multiplies the matrix by its adjoint and then looks
to see how well the result conforms to an identity matrix.  This can use up
lots of time....
                          ---
                          \                 *t
                 z(i,j) = /    <i|mx|k><k|mx  |j> = del   +/- d
                          ---                          ij                   */
        
MSVCDLL bool is_symmetric(const double d=GMxCut) const;
MSVCDLL bool is_hermitian(const double d=GMxCut) const;
MSVCDLL bool is_unitary(const   double d=GMxCut) const;
MSVCDLL bool is_real(const      double d=GMxCut) const;
MSVCDLL bool is_imaginary(const double d=GMxCut) const;
MSVCDLL bool is_complex(const   double d=GMxCut) const;
MSVCDLL bool is_zero(const      double d=GMxCut) const;
MSVCDLL bool is_diagonal(const  double d=GMxCut) const;
MSVCDLL bool is_square()                         const;
 
// ____________________________________________________________________________
// F                 CLASS MATRIX BINARY ARITHMETIC FUNCTIONS
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

MSVCDLL matrix operator + (const matrix& mx) const;
MSVCDLL matrix operator - (const matrix& mx) const;
MSVCDLL matrix operator * (const matrix& mx) const;
MSVCDLL matrix operator * (const complex& z) const;
MSVCDLL matrix operator * (      double   d) const;
MSVCDLL matrix operator / (const matrix& mx) const;
MSVCDLL matrix operator / (const complex& z) const;
MSVCDLL matrix operator / (      double   d) const;

MSVCDLL friend matrix operator * (const complex& z,  const matrix& mx);
MSVCDLL friend matrix operator * (      double d,    const matrix& mx);
//MSVCDLL friend matrix operator / (const complex& z,  const matrix& mx);
//MSVCDLL friend matrix operator / (      double d,    const matrix& mx);

// ____________________________________________________________________________
// G                  CLASS MATRIX UNARY ARITHMETIC FUNCTIONS
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
 
MSVCDLL matrix& operator += (const matrix& mx1);
MSVCDLL matrix& operator -= (const matrix& mx1);
MSVCDLL matrix& operator *= (const matrix& mx);
MSVCDLL matrix& operator *= (const complex& z);
MSVCDLL matrix& operator *= (double d);
MSVCDLL matrix& operator /= (const matrix& mx);
MSVCDLL matrix& operator /= (const complex& z);
MSVCDLL matrix& operator /= (double d);

// ____________________________________________________________________________
// H                  CLASS MATRIX SIMPLE UNARY FUNCTIONS
// ____________________________________________________________________________

/*              Input          mx : Matrix (this)
                Output         mx1: A matrix derived from mx
                Note              : Operations performed _matrix derived class

    Function/Operator           Output          Description
    -----------------       --------------      -------------------------------
        operator-             mx1 = -mx         negation
           Re                 mx1 = Re{mx}      zero imaginaries
           Im                 mx1 = IM{mx}      zero reals
          conj                mx1 = mx*         conjugate   <i|mx1|j>=<i|mx*|j>
        transpose             mx1 = mxt         transpose   <i|mx1|j>=<j|mx|i>
        adjoint               mx1 = mxt*        adjoint     <i|mx1|j>=<j|mx*|i>
          exp                 mx1 = exp(mx)     exponential <i|mx1|j>=<j|mx*|i>
         trace                  z = Tr{mx}      trace       z = sum <i|mx|i> */

MSVCDLL friend matrix  Re(const matrix& mx);
MSVCDLL friend matrix  Im(const matrix& mx);
MSVCDLL friend matrix  conj(const matrix& mx);
MSVCDLL friend matrix  transpose(const matrix& mx);
MSVCDLL friend matrix  adjoint(const matrix& mx);
MSVCDLL friend complex trace(const matrix& mx);
MSVCDLL friend void enable_blockdiag();
MSVCDLL friend void disable_blockdiag();

MSVCDLL matrix  operator- () const;
MSVCDLL matrix  Re()         const;
MSVCDLL matrix  Im()         const;
MSVCDLL matrix  conj()       const;   
MSVCDLL matrix  transpose()  const; 
MSVCDLL matrix  adjoint()    const;
MSVCDLL matrix  exp()        const;
MSVCDLL complex trace()      const;    

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

MSVCDLL matrix  swaprows(int i, int j);
MSVCDLL matrix  swapcols(int i, int j);
MSVCDLL matrix  permute( int i, int j);
MSVCDLL double  maxRe() const;
MSVCDLL double  maxIm() const;
MSVCDLL complex maxZ()  const;
MSVCDLL double  minRe() const;
MSVCDLL double  minIm() const;
MSVCDLL complex minZ()  const;

// ____________________________________________________________________________
// I                CLASS MATRIX SIMPLE BINARY FUNCTIONS
// ____________________________________________________________________________

/*     Function    Output                       Description
    -------------  -------  ---------------------------------------------------
        trace      complex  Returns trace of matrix product: Tr{mx1*mx2}
    adjoint_times  matrix   Returns adjoint times array:     adjoint(mx)*mx1
    times_adjoint  matrix   Returns array times adjoint:     mx*adjoint(mx1) */

MSVCDLL friend complex trace(const matrix& mx1, const matrix& mx2);
MSVCDLL complex trace(const matrix& mx2) const;
MSVCDLL friend matrix adjoint_times(const matrix& mx, const matrix& mx1);
MSVCDLL friend matrix times_adjoint(const matrix& mx, const matrix& mx1);

// ____________________________________________________________________________
// J               CLASS MATRIX COMPLEX UNARY FUNCTIONS
// ____________________________________________________________________________

MSVCDLL friend  complex         det(const matrix& mx);
MSVCDLL         complex         det() const;
MSVCDLL friend  int             rank(const matrix& mx);

// ***************************      FFT      ****************************

/* These functions handle Fast Fourier Transforms (FFT) and inverse
   Fast Fourier Transforms (IFFT).

	 Input	      mx  : Matrix to be transformed
	 Output      mx1 : Transformed matrix
	 Note	          : The static matrix flag FFTComp sets whether
                          the transform is "folded" or not    n    n
	 Note            : Matrix dimension must be either (1x2 or 2 x1)*/

MSVCDLL friend matrix FFT(const  matrix& mx);
MSVCDLL friend matrix IFFT(const matrix& mx);

MSVCDLL matrix FFT()  const;
MSVCDLL matrix IFFT() const;

// ____________________________________________________________________________
// K            CLASS MATRIX COMPLEX BINARY (&TRINARY) FUNCTIONS
// ____________________________________________________________________________

//------------------------------------------------------------------------------
//                                Tensor Product
//------------------------------------------------------------------------------

	// Input     mx1,mx2 : Matrices
	// Output 	 mx3 : Matrix tensor product of mx1 & mx2
	//			   mx3 = mx1 (x) mx2
        //                        (mx1->rows()*mx2->rows(),
        //                         mx1->cols()*mx2->cols()) with the elemts:
        //                        (i*mx2->rows()+j,k*mx2->cols()+l)=
        //                         mx1(i,k)*mx2(j,l)

MSVCDLL friend matrix tensor_product(const matrix& mx1, const matrix& mx2);

// ____________________________________________________________________________
// L                     CLASS MATRIX I/O FUNCTIONS
// ____________________________________________________________________________

//------------------------------------------------------------------------------
//                      ASCII OUTPUT CONTROL FUNCTIONS
//------------------------------------------------------------------------------

/* These functions control how arrays will be output. Those with a "*" apply
   only to matrix numerical element output (non-pictorial).  Those with a +
   apply to pictorial output only. Control of numerical element output format
   is set in class complex.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
     Header     bool       Flags whether to write a header above the matrix
   *PrintRI     bool      *Flags whether to write only reals/imags if possible
    PrintAll    bool       Flags whether to write all elements (non-stored)
    PictDim     int        Dimension where output switches to pictorial rep.
   +PrintVal    bool      +Flags whether valuate elems. in pictorial output  */

MSVCDLL static void    Header(bool   hf);
MSVCDLL static void    PrintRI(bool  pi);
MSVCDLL static void    PrintAll(bool pa);
MSVCDLL static void    PictDim(int   pd);
MSVCDLL static void    PrintVal(bool pv);
MSVCDLL static void    PrintCols(int cl);
MSVCDLL static void    PrintRows(int rl);
MSVCDLL static MxPrint PrintFlags();

//-----------------------------------------------------------------------------
//                          ASCII OUTPUT FUNCTIONS
//-----------------------------------------------------------------------------

/* These functions allow users to write matrices in formatted ASCII to an 
   output stream.  They also allow users to obtain strings representing 
   array rows which allows for full flexibility in array output. 

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
    printHdr   ostream     Writes header for placement above the array  
     print     ostream     Writes array elemnts in ASCII to output stream
    picture    ostream     Writes array pictorially to output stream         */

MSVCDLL        std::ostream& printHdr(std::ostream& ostr) const;
MSVCDLL        std::ostream& print(std::ostream&    ostr) const;
MSVCDLL        std::ostream& picture(std::ostream&  ostr) const;
MSVCDLL friend std::ostream&  operator << (std::ostream& ostr, const matrix& mx);

//-----------------------------------------------------------------------------
//                            BINARY I/O FUNCTIONS
//-----------------------------------------------------------------------------

/* These functions perform both ASCII and Binary input/output operations on a
   matrices.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   write      ofstream    Writes array elemnts in Binary to output file stream
   read       ifstream    Read in array from Binary input file stream

   Binary I/O the format of matrices is partially dictated by the flag "form".
   By default (form=0) the output will be GAMMA format which is compact in that
   special matrix types dictate how many elements are output.  For form!=0
   all array elements are output even if they don't exist in memory.  For the
   former the array type, # rows, # columns, & data are all output.  For the
   latter only the data is output (generic array).  The data ordering is
   Re(<i|hmx|j>, Im(<i|hmx|j> columns then rows (i.e. row by row.)  Note that
   binary read & write functions MUST exactly match in format!               */  

MSVCDLL        std::ofstream& write(std::ofstream& fp, int form=0) const;
MSVCDLL friend std::ofstream& write(std::ofstream& F,  const matrix& mx);
MSVCDLL        std::ifstream& read(std::ifstream&  F);

//-----------------------------------------------------------------------------
//                            ASCII INPUT FUNCTIONS
//-----------------------------------------------------------------------------

/* These functions perform ASCII input of arrays.

    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   readASC    istream     Read in array from an input stream    
   ask          ----      Interactively requests matrix from user            */

MSVCDLL friend std::istream&  operator >> (std::istream& istr, matrix& mx);
MSVCDLL        void           ask(const matrix_type t=n_matrix_type);

// ____________________________________________________________________________
// M                   CLASS MATRIX MISCELLANEOUS FUNCTIONS
// ____________________________________________________________________________

/* These are functions that don't fit into the other categories.
 
    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   resize      nr, nc     Resizes mx to have nr rows & nc columns, nr*nc=const
   copy         ----      Produces a copy of the hmx, allocates memory
   convert       mx       Matrix mx has its values set to those of hmx (no mem)
 
   The conversion functions will often drop elements, this is unavoidable. For
   example, one cannot generally convert a Hermitian array to a diagonal array.
   There are no elements lost when using HMX or NMX.  The *MX functions
   will always allocate new memory for the converted matrix from hmx.      */

MSVCDLL matrix resize(int i, int j);
MSVCDLL matrix diagonal_form();

	// Input        this : Matrix
	// Output 	     : The matrix this is resized to be
	//		       larger with its elements on the diagonal

MSVCDLL bool same_reference_as(const matrix& mx) const;

	// Input        this : Matrix
	//		  mx : Matrix
        // Output            : Returns TRUE if both matricies
	//		       reference to the identical storage
	//		       location (& thus are identical arrays)

MSVCDLL void status(int full=0) const;

	// Input		mx    : Matrix (this)
	//			full  : Flag for amount of output
	// Output		void  : Outputs matrix status

// ____________________________________________________________________________
// N                CLASS MATRIX DIAGONALIZATION FUNCTIONS
// ____________________________________________________________________________

/* These functions reflect individual steps used in general complex matrix
   diagonalizations. Since how a diagonalization is performed depends upon the
   underlying matrix type, the generic functions will handle this process,
   namely "Diagonalize" or "diag".  The other functions allow the user to
   monitor individual steps taken in the general case. 

   Function   Arguments                         Result
 ------------ --------- ------------------------------------------------------
  BlockDiag    mx, mx   Get blocked array form:          mx = U * BD  * inv(U)
 SymTriDiag    mx, mx   Get symmetric tridiagonal form:  mx = U * STD * inv(U)
 HermTriDiag   mx, mx   Get Hermitian tridiagonal array: mx = U * HTD * inv(U)
   SymDiag     mx, mx   Diagonalize symmetric array:     mx = U * SD  * inv(U)
 Diagonalize   mx, mx   Diagonalize array:               mx = U *  D  * inv(U)
    diag       mx, mx   Diagonalize array (non-member):  mx = U *  D  * inv(U)

 In all cases, the original array will remain unchanged and the arguments will
 will be set to the new matrix and transformation array (U) that will convert
 between the new form and the original array. In the case of block diagonal-
 ization a vector of block sizes will be returned. When possible, specialized
 routines will be used based on the matrix type and the return arrays will 
 be of the best suited form (e.g diagonal, block, unitary, etc....)          */

MSVCDLL   std::vector<int> BlockDiag(matrix&    BD, std::vector<int> &U) const;
MSVCDLL        void        SymTriDiag(matrix&  HTD, matrix& U) const;
MSVCDLL        void        HermTriDiag(matrix& STD, matrix& U) const;
MSVCDLL        void        SymDiag(matrix&      SD, matrix& U) const;
MSVCDLL        void        Diagonalize(matrix&   D, matrix& U) const;

MSVCDLL friend void diag(const matrix& mx,     matrix&   D, matrix& U);
 
// ____________________________________________________________________________
// O                    CLASS MATRIX INVERSION FUNCTIONS
// ____________________________________________________________________________

/* Herein, we use the LU decomposition to obtain a matrix inverse. In this
   case we have (using A=matrix , L=lower triang., & U=upper triang.)

                               -1               -1
                          A * A  = I = L * U * A
  
   where we 1st solve L*Y=I for Y, then U*inv(A)=Y for inv(H) taking advantage
   of the special formats of L & U.

  The LU & LU inverse routines used here are derived from "Numerical Recipies",
  W.H.  Press, B.P.Flannery, S.A.Teukolsky, & W.T.Vetterling, Cambridge 
  University Press, 1989.  See pages 31 through 36. The code was extensively
  adapted for C++ and GAMMA arrays.

  The LU algorithm takes as input a matrix A and an integer array indx.
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
     inv        ---    inverse(mx) where inv(mx)*mx = I
     LU        int*    LU == mx + row permutations         One Array L\U Form
     LUinv    int*, B  Solutionh X where LU*X = AX = B     Inverse if B=I   */

MSVCDLL friend matrix inv(const matrix& mx);
MSVCDLL friend matrix LU(matrix& mx, int* indx);
MSVCDLL friend matrix LUinv(matrix& B, int *indx, matrix& LU);

// ____________________________________________________________________________
// P                     CLASS MATRIX TESTING FUNCTIONS
// ____________________________________________________________________________

/* Since GAMMA matrices and the matrix class hierarchy are quite complicated
   I have directly included some core testing functions to check that the
   diagonalization routines are working properly.

       Function    Arguments                           Result
   --------------- ---------   ------------------------------------------------
   TestEigenSystem    int      Diagonalize the array & check its eigensystem.
   TestTransform   mx,mx,int   Check if mx can be remade from its eigensystem.
   ColumNorms         ---      Will output column norms (for eigenvectors).
   TestIdentity     complex    Check if array is idenity with z (A*inv(A)=I)
   TestUniary        mx, mx    Constructs the product U*adjoint(U) & tests
   TestUTransform    mx,mx     Test unitary tranform U*T*adj(U)              */

MSVCDLL void                TestEigenSystem(int pf=1) const;
MSVCDLL void                TestTransform(const matrix& T, const matrix& S, int pf=1) const;
MSVCDLL std::vector<double> ColumnNorms() const;
MSVCDLL std::vector<double> TestIdentity(complex& TotalDev) const;
MSVCDLL matrix              TestUnitary(std::ostream& ostr) const;
MSVCDLL matrix              TestUTransform(const matrix& T, const matrix& U) const;

// ____________________________________________________________________________
// Q                   MATRIX CONTAINER SUPPORT FUNCTIONS
// ____________________________________________________________________________

/* Aside for providing basic tests as to whether to matrics are equvalent or
   not, these operators are necessary in any STL container classes are to
   be used based on matrices (e.g. list<matrix> or vector<matrix>)           */

MSVCDLL bool operator== (const matrix& mx) const;
MSVCDLL bool operator!= (const matrix& mx) const;
MSVCDLL bool operator<  (const matrix& mx) const;
MSVCDLL bool operator>  (const matrix& mx) const;

};

#endif 						                  // matrix.h
