/* _matrix.cc ***************************************************-*-c++-*-
**									**
**                               G A M M A				**
**									**
**	Prototype Matrix                           Implementation	**
**									**
**	Copyright (c) 1990, 1991, 1992					**
**	Tilo Levante, Scott A. Smith					**
**	Eidgenoessische Technische Hochschule				**
**	Labor fuer physikalische Chemie					**
**	8092 Zuerich / Switzerland					**
**									**
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The class _matrix defines a "dummy matrix" in C++ for GAMMA.  This   **
** serves as a proto-type and front end to all classes which handle the **
** specific matrix types.  The externally used matrix class, matrix,    **
** routes through this  dummy class in finding specialized routines.    **
**                                                                      **
*************************************************************************/

#ifndef   G_matrix_cc_			// Is file already included?
#  define G_matrix_cc_ 1		// If no, then include it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Matrix/_matrix.h>		// Include the interface
#include <Matrix/MxModBas.h>		// Include Matrix module errors
#include <fstream>                      // Include libstdc++ file streams
#include <Matrix/i_matrix.h>            // Include identity interface
#include <Matrix/d_matrix.h>		// Include diagonal matrices
#include <Matrix/h_matrix.h> 		// Include Hermitian matrices
#include <Matrix/n_matrix.h>		// Include normal matrices
#include <vector>			// Include libstdc++ vectors
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i                      CLASS _MATRIX ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   mx      : A _matrix (this)                      
				CL	: Matrix class name (_matrix derived)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : string in message
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void _matrix::_MxError(int eidx, int nr) const
  { std::string CL="Base Matrix"; Mxerror(CL, eidx, nr); }
void _matrix::_MxError(int eidx, const std::string& PN, int nr) const
  { std::string CL="Base Matrix"; Mxerror(CL, eidx, PN, nr); }

void _matrix::Mxerror(const std::string& CL, int eidx, int noret) const
  {
  std::string msg;
  switch (eidx)
    {
    case 1:  MxModError(CL,"Problems During Construction",noret); break;// (1)
    case 2:  MxModError(CL,"Rectangular Array Forbidden",noret);  break;// (2)
    case 3:  MxModError(CL,"Multiply Size Incompatibility",noret);break;// (3)
    case 4:  MxModError(CL,"Create Rectangular Identity Matrix?");break;// (4)
    case 6:  MxModError(CL,"Bad Internal Component Access",noret);break;// (6)
    case 9:  MxModError(CL,"Problems During Construction",noret); break;// (9)
    case 10: MxModError(CL,"Block Put Exceeds Matrix Dimension"); break;// (10)
    case 11: MxModError(CL,"Resize Attempt to Rectangular Array");break;// (11)
    case 12: MxModError(CL,"Division Not Fully Implemented");     break;// (12)
    case 13: MxModError(CL,"Division By Zero");                   break;// (13)
    case 14: MxModError(CL,"Bad Rectangular Array Use",noret);    break;// (14)
    case 20: MxModError(CL,"Unable To Perform Addition",noret);   break;// (20)
    case 21: MxModError(CL,"Cannot Perform Subtraction",noret);   break;// (21)
    case 22: MxModError(CL,"Unable To Do Multiplication",noret);  break;// (22)
    case 23: MxModError(CL,"Unable To Do Division",noret);        break;// (23)
    case 24: MxModError(CL,"Unable To Take Trace",noret);         break;// (24)
    case 25: MxModError(CL,"Unable To Take Inverse", noret);      break;// (25)
    case 26: MxModError(CL,"Problems In LU Decomposition",noret); break;// (26)
    case 27: MxModError(CL,"Array Inverse Does Not Exist",noret); break;// (27)
    case 28: MxModError(CL,"Unable To Do Diagonalization",noret); break;// (28)
    case 29: MxModError(CL,"Too Many Function Iterations",noret); break;// (29)
    case 30: MxModError(CL,"Row<->Col Col<->Row Mismatch",noret); break;// (30)
    case 31: MxModError(CL,"Row<->Row Col<->Col Mismatch",noret); break;// (31)
    case 42: MxModError(CL,"Problems During Assignment",noret);   break;// (42)
    case 50: MxModError(CL,"Array May Contain only 1's or 0's!"); break;// (50)
    case 51: MxModError(CL,"Mixing Mismatched Dimensions",noret); break;// (51)
    case 52: MxModError(CL,"Array Dimensions Exceeded",noret);    break;// (52)
    case 56: MxModError(CL,"Element Access Out Of Range",noret);  break;// (56)
    case 67: MxModError(CL,"Cannot Access Element",noret);        break;// (67)
    case 70: MxModError(CL,"Cannot Convert To Hermitian",noret);  break;// (70)
    case 71: MxModError(CL,"Cannot Convert To Diagnoal", noret);  break;// (71)
    case 72: MxModError(CL,"Cannot Convert To Identity", noret);  break;// (72)
    case 73: MxModError(CL,"Cannot Convert To Complex", noret);   break;// (73)
    case 80: MxModError(CL,"Dimen. Mismatch in A*X=LU*X=B",noret);break;// (80)
    case 81: MxModError(CL,"Cannot Continue Calculations", noret);break;// (81)
    default: MxModError(CL, eidx, noret);                         break;
    }
  }  

void _matrix::Mxerror(const std::string& CL, int eidx,
                                  const std::string& pname, int noret) const
  {                                                                             
  std::string msg;
  switch(eidx)
    {
    case 5: msg = std::string("Bad Use Of Function ");
             MxModError(CL,msg+pname,noret);  break;                  // (5)
    case 21:msg = std::string("Disallowed Use of Function ");
             MxModError(CL,msg+pname,noret);  break;                  // (21)
    case 22:msg = std::string("Disallowed Use of Operator ");
             MxModError(CL,msg+pname,noret);  break;                  // (22)
    case 23:msg = std::string("Disallowed Use of Array Access ");
             MxModError(CL,msg+pname,noret);  break;                  // (23)
    case 24:msg = std::string("Disallowed Use of Array I/O ");
             MxModError(CL,msg+pname,noret);  break;                  // (24)
    case 25:msg = std::string("Sorry, Function ") + pname 
                + std::string(" Not Fully Implemented");
             MxModError(CL,msg,noret);  break;                        // (25)
    default: MxModError(CL, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }  

volatile void _matrix::_MxFatal(int eidx) const
  { std::string CL="Base Matrix"; Mxfatality(CL, eidx); }

volatile void _matrix::_MxFatal(int eidx, const std::string& PN) const
  { std::string CL="Base Matrix"; Mxfatality(CL, eidx, PN); }   

volatile void _matrix::Mxfatality(const std::string& CL, int eidx) const
  {                                                                 
  Mxerror(CL, eidx, 1);			// Normal non-fatal error
  if(eidx) Mxerror(CL, 0);		// Program aborting error
  MxModFatal();	
  }

volatile void _matrix::Mxfatality(const std::string& CL, int eidx,
                                                 const std::string& PN) const
  {                                                                 
  Mxerror(CL, eidx, PN, 1);			// Normal non-fatal error
  if(eidx) Mxerror(CL, 0);			// Program aborting error
  MxModFatal();					// End program
  }


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A               CLASS _MATRIX CONSTRUCTORS AND DESTRUCTOR
// ____________________________________________________________________________

/*  Arguments                      Constructed Array
    ---------           -------------------------------------------------------
        -                An empty array
     nr, nc              An (nr x nc) array, uninitialized
       _mx               A duplicate of prototype array mx

   Note that all matrix classes that are derived from _matrix will invariably
   call these constructors when they are built, so that their row and column
   sizes are set. However, here we do not allocate any space since the purpose
   of the derive classes is to minimize the number of elements stored.       */  

_matrix::_matrix ( )
  {
  rows_      = 0;			// Default matrix is size 0x0 
  cols_      = 0;
  size       = 0;			// Nothing stored
  references_= 0;
  }

_matrix::_matrix(int i, int j)
  {
  rows_      = i;			// Size of the matrix	
  cols_      = j;
  size       = 0;			// Nothing stored
  references_= 0;
  }

_matrix::_matrix(const _matrix& m)
  {
  rows_      = m.rows_;			// Size of the matrix	
  cols_      = m.cols_;
  size       = m.size;			// Nothing stored
  references_= 0;
  }

_matrix::~_matrix () { }				// Nothing to destruct

 
// ____________________________________________________________________________
// B                   CLASS _MATRIX ACCESS AND ASSIGNMENT
// ____________________________________________________________________________

                                                                                
// ************************* assignment for _matrix ***************************

	// Input                _mx     : A _matrix (this)
        // Output       	m1	: A copy of the _matrix mx
	// Note				: Does not modify the reference count!

void _matrix::operator= (const _matrix& m)
  {
  if(this == &m) return;			// Ignore if self-assign
  rows_ 	= m.rows_;			// Set the row size
  cols_ 	= m.cols_;			// Set the column size
  size  	= m.size;			// Set the data size
  }


// ******************* _matrix element access functions ***********************

complex& _matrix::operator() (int i, int j)
             { _MxFatal(23,"Element Access mx(i,j)"); i=0; j=0; return ZNULL; }
 
complex _matrix::get(int i, int j) const
         { _MxFatal(23,"Element Access get(i,j)"); i=0; j=0; return complex0; }
 
bool _matrix::put(const complex& z, int i, int j)
   { _MxFatal(23,"Element Set put(i,j)"); i=0; j=int(Re(z)); return false; }

bool _matrix::put_h(const complex& z, int i, int j) 
 { _MxFatal(23,"Element Set put_h(i,j)"); i=0; j=int(Re(z)); return false; }

// ********************* _matrix block access functions ***********************

_matrix* _matrix::get_block(int row, int col, int nr, int nc)
  { _MxFatal(23,"Get Array Block"); row=0; col=0; nr=0; nc=0; return this; }

bool _matrix::put_block(int row, int col, _matrix* m)
       { _MxFatal(23,"Set Array Block"); row=0; col=0; m=NULL; return false; }

// ******************* _matrix details access functions ***********************

const int _matrix::rows()  { return rows_; }		// Handled directly
const int _matrix::cols()  { return cols_; }		// Handled directly
const int _matrix::refs()  { return references_; }	// Handled directly
const int _matrix::pts()   { return size; }		// Handled directly
int& _matrix::references() { return references_; }	// Handled directly

// ____________________________________________________________________________
// C        CLASS _MATRIX HERMITIAN_TYPE & MATRIX_TYPE HANDLING
// ____________________________________________________________________________


hermitian_type _matrix::stored_hermitian( )      const {return _hermitian; }
hermitian_type _matrix::test_hermitian(double d) const
                                                       {return _hermitian;d=0;}

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
   mxtype                string      Returns the string "Null" 
 
   The test type looks to see if _mx could be an array of the matrix type
   specified as an input argument.  If it can within "d" the return type 
   is equal to the input type.  If it cannot the return is _matrix_type     */

// Note:  Functions set set_type & check_type are in class matrix, not here.

matrix_type _matrix::stored_type( ) const             { return _matrix_type; }
matrix_type _matrix::test_type(matrix_type m, double d) const
                               { return _matrix_type; d=0.0; m=_matrix_type; }
std::string _matrix::mxtype() const                      { return std::string("Null"); }
std::string _matrix::mxtype(bool pf) const               { return std::string("Null"); }

// ____________________________________________________________________________
// E                   CLASS _MATRIX VARIOUS MATRIX CHECKS
// ____________________________________________________________________________

/* These functions return either TRUE(1) or FALSE(0).  We define _matrix as
   being non-unitary, real, zero, non-square.  Note that the setting of d
   in these function just placates the compiler which likes d to be used 

     Function    Output                       Description
   ------------  ------  ------------------------------------------------------
   is_symmetric   int    TF if Im(<i|_mx|j>) < d for all i,j (def. d GMxCut)
   is_hermitian   int    TF if Im(<i|_mx|j>) < d for all i,j (def. d GMxCut)
   is_unitary     int    TF if inv(_mx) == adjoint mx, CPU intensive
   is_real        int    TF if Im(<i|_mx|j>) < d for all i,j (def. d GMxCut)
   is_imaginary   int    TF if Re(<i|_mx|j>) < d for all i,j (def. d GMxCut)
   is_complex     int    TF if is_real && is_imaginary       (def. d GMxCut)
   is_zero        int    TF if ||<i|_mx|j>|| < d for all i,j (def.d GMxCut)
   is_diagonal    bool   TF if ||<i|_mx|j>|| < d for all i!=j(def. d GMxCut)
   is_square      bool   TF if rows_ == cols_
   is_equal       int    TF if ||<i|_mx-mx|j>||<d    all i,j (def.d GMxCut)  */

bool _matrix::is_symmetric(double d) const { return true;  d = 0.0; }
bool _matrix::is_hermitian(double d) const { return true;  d = 0.0; }
bool _matrix::is_unitary(double d)   const { return false; d = 0.0; }
bool _matrix::is_real(double d)      const { return true;  d = 0.0; }
bool _matrix::is_imaginary(double d) const { return false; d = 0.0; }
bool _matrix::is_complex(double d)   const { return false; d = 0.0; }
bool _matrix::is_zero(double d)      const { return true;  d = 0.0; }
bool _matrix::is_diagonal(double d)  const { return true;  d = 0.0; }
bool _matrix::is_square()            const { return false; }

bool _matrix::is_equal(_matrix* mx, double d) const
                              { return ((*mx).stored_type() == _matrix_type); }

// ____________________________________________________________________________
// F                CLASS _MATRIX BINARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________

// All Arithmetic Function Are Handled By Derived Classes of _matrix.  Thus,
// Access of These ALWAYS Produces A Fatal Error.  Users CANNOT Perform Any
// Arithmetic Operations On NULL Matrices Which Have No Dimension Or Data.

// --------------------- Matrix With Matrix Operations ------------------------

_matrix* _matrix::add(_matrix* mx)      { _MxFatal(22, "+"); return mx; }
_matrix* _matrix::subtract(_matrix* mx) { _MxFatal(22, "-"); return mx; }
_matrix* _matrix::multiply(_matrix* mx) { _MxFatal(22, "*"); return mx; }
_matrix* _matrix::divide(_matrix* mx)   { _MxFatal(22, "/"); return mx; }

// --------------------- Matrix With Scalar Operations ------------------------

_matrix* _matrix::multiply(const complex& z) 
                                       { _MxFatal(22, "*"); return this; }
_matrix* _matrix::divide(const   complex& z)
                                       { _MxFatal(22, "/"); return this; }

// ____________________________________________________________________________
// G                 CLASS _MATRIX UNARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________

// All Arithmetic Function Are Handled By Derived Classes of _matrix.  Thus,
// Access of These ALWAYS Produces A Fatal Error.  Users CANNOT Perform Any
// Arithmetic Operations On NULL Matrices Which Have No Dimension Or Data.
 
// --------------------- Matrix With Matrix Operations ------------------------

_matrix* _matrix::add_two(_matrix* mx)      { _MxFatal(22, "+"); return mx; }
_matrix* _matrix::subtract_two(_matrix* mx) { _MxFatal(22, "-"); return mx; }
_matrix* _matrix::multiply_two(_matrix* mx) { _MxFatal(22, "*"); return mx; }
_matrix* _matrix::divide_two(_matrix* mx)   { _MxFatal(22, "/"); return mx; }

// --------------------- Matrix With Scalar Operations ------------------------

_matrix* _matrix::multiply_two(const complex& z)
                                            { _MxFatal(22, "*"); return this; }
_matrix* _matrix::divide_two(const complex& z)
                                            { _MxFatal(22, "/"); return this; }

 
// ____________________________________________________________________________
// H                  CLASS _MATRIX SIMPLE UNARY FUNCTIONS
// ____________________________________________________________________________

/* All Arithmetic Function Are Handled By Derived Classes of _matrix.  Thus,
   Access of These ALWAYS Produces A Fatal Error.  Users CANNOT Perform Any
   Arithmetic Operations On NULL Matrices Which Have No Dimension Or Data.   */
 
_matrix* _matrix::negate()    {_MxFatal(5, "Negation");          return this; }
_matrix* _matrix::RE()        {_MxFatal(5, "Real");              return this; }
_matrix* _matrix::IM()        {_MxFatal(5, "Imaginary");         return this; }
_matrix* _matrix::conjugate() {_MxFatal(5, "Complex Conjugate"); return this; }
_matrix* _matrix::transpose() {_MxFatal(5, "Transpose");         return this; }
_matrix* _matrix::adjoint()   {_MxFatal(5, "Adjoint");           return this; }
_matrix* _matrix::mxexp()     {_MxFatal(5, "Exponential");       return this; }
complex  _matrix::trace()     {_MxFatal(5, "Trace");         return complex0; }

_matrix* _matrix::swaprows(int i,int j) {_MxFatal(5,"swaprows"); return this;}
_matrix* _matrix::swapcols(int i,int j) {_MxFatal(5,"swapcols"); return this;}
_matrix* _matrix::permute( int i,int j) {_MxFatal(5, "permute"); return this;}
double   _matrix::maxRe()   const { _MxFatal(5, "maxRe");   return 0;        }
double   _matrix::maxIm()   const { _MxFatal(5, "maxIm");   return 0;        }
complex  _matrix::maxZ()    const { _MxFatal(5, "maxZ");    return complex0; }
double   _matrix::minRe()   const { _MxFatal(5, "minRe");   return 0;        }
double   _matrix::minIm()   const { _MxFatal(5, "minIm");   return 0;        }
complex  _matrix::minZ()    const { _MxFatal(5, "minZ");    return complex0; }
 
// ____________________________________________________________________________
// I                  CLASS _MATRIX SIMPLE BINARY FUNCTIONS
// ____________________________________________________________________________

/* All Arithmetic Function Are Handled By Derived Classes of _matrix.  Thus,
   Access of These ALWAYS Produces A Fatal Error.  Users CANNOT Perform Any
   Arithmetic Operations On NULL Matrices Which Have No Dimension Or Data.   */
 
complex  _matrix::trace(_matrix *mx) 
                             { _MxFatal(5, "Trace");      return complex0; }
_matrix* _matrix::adjoint_times (_matrix* mx)
                             { _MxFatal(5, "Adjoint-Times"); return this;  }
_matrix* _matrix::times_adjoint(_matrix* mx)
                             { _MxFatal(5, "Times-Adjoint"); return this;  }

// ____________________________________________________________________________
// J                  CLASS _MATRIX COMPLEX UNARY FUNCTIONS
// ____________________________________________________________________________

/* All Arithmetic Functions Are Handled By Derived Classes of _matrix.  Thus,
   Access of These ALWAYS Produces A Fatal Error.  Users CANNOT Perform Any
   Arithmetic Operations On NULL Matrices (Which Have No Dimension Or Data.) */

_matrix* _matrix::inv() { _MxFatal(5, "Inverse");                return this; }
_matrix* _matrix::LU(int* indx) 
                   { _MxFatal(5, "LU Decomposition"); indx=NULL; return this; }
_matrix* _matrix::LUinv(int* indx, _matrix* LU)
                { _MxFatal(5, "LU Inverse"); indx=NULL; LU=NULL; return this; }

std::vector<int> _matrix::BlockDiag(_matrix* (&BF), std::vector<int> &U)
             { _MxFatal(5, "Block Diagonalization"); BF=NULL; U=std::vector<int>();
                                                        return std::vector<int>(); }
void        _matrix::HermTriDiag(_matrix* (&HTD), _matrix* (&U))
             { _MxFatal(5, "Hermitian TriDiagonalization"); HTD=NULL; U=NULL; }
void        _matrix::SymTriDiag(_matrix* (&STD), _matrix* (&U))
             { _MxFatal(5, "Symmetric TriDiagonalization"); STD=NULL; U=NULL; }
void        _matrix::SymDiag(_matrix* (&D), _matrix* (&U))
             { _MxFatal(5, "Diagonalization"); D=NULL; U=NULL; }
void        _matrix::diag(_matrix* (&D), _matrix* (&U))
             { _MxFatal(5, "Diagonalization"); D=NULL; U=NULL; }

complex     _matrix::det() { _MxFatal(5, "Determinant"); return complex0; }
int             _matrix::rank() { _MxFatal(5, "Rank");       return 0;        }

// ____________________________________________________________________________
// K                CLASS _MATRIX COMPLEX BINARY FUNCTIONS
// ____________________________________________________________________________
 

_matrix* _matrix::tensor_product(_matrix* mx)
                                  { _MxFatal(5, "Tensor Product"); return mx; }

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

std::vector<std::string> _matrix::printStrings(const   MxPrint& pf) const {return std::vector<std::string>(); }
std::vector<std::string> _matrix::pictureStrings(const MxPrint& pf) const {return std::vector<std::string>(); }

void _matrix::print(std::ostream&   ostr, const MxPrint& PFlgs) const {return;}
void _matrix::picture(std::ostream& ostr, const MxPrint& PFlgs) const {return;}

void _matrix::write(std::ofstream &fstr, int form) const
                         { _MxFatal(5, "Binary Write"); form=0; fstr.close(); }

void _matrix::read(std::ifstream &fstr)  
                                  { _MxFatal(5, "Binary Read"); fstr.close(); }

void _matrix::readASC(std::istream& istr)
                                 { _MxFatal(5, "ASCII Read"); int i; istr>>i; }

void _matrix::ask( )                        { _MxFatal(5, "Interactive Ask"); }


// ____________________________________________________________________________
// M                 CLASS _MATRIX MISCELLANEOUS FUNCTIONS
// ____________________________________________________________________________


void _matrix::resize(int i, int j) { rows_ = i; cols_ = j; }

_matrix* _matrix::copy() { return new _matrix(*this); }

void _matrix::convert(_matrix* mx) { _MxFatal(21,"convert"); mx=NULL; }

 
// ------------------------ New Conversion Routines ---------------------------
//                  (Been Having Troubles With The Old One)
 

i_matrix* _matrix::IMX() { _MxFatal(21,"IMX"); return new i_matrix(1,1); }
d_matrix* _matrix::DMX() { _MxFatal(21,"DMX"); return new d_matrix(1,1); }
h_matrix* _matrix::HMX() { _MxFatal(21,"HMX"); return new h_matrix(1,1); }
n_matrix* _matrix::NMX() { _MxFatal(21,"NMX"); return new n_matrix(1,1); }

        // Input                _mx     : A _matrix (this)
        //                      mxtype  : A matrix type
        // Output               mx      : A new matrix of type mxtype is made
        //                                whose elements are the same as _mx
    
_matrix* _matrix::convert(matrix_type mxtype)
  { _MxFatal(5, "convert"); mxtype=_matrix_type; return this; }

#endif							// _matrix.cc
