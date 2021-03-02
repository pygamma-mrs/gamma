/* matrix.cc ****************************************************-*-c++-*-
**                                                                      **
**                             G A M M A				**
**                                                                      **
**	Matrix		                           Implementation 	**
**                                                                      **
**      Copyright (c) 1990, 1991, 1992                                  **
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
** The class matrix defines matrices for GAMMA in C++ with the usual    **
** operations +, -, *, /, and I/O routines.                             **
**                                                                      **
** It uses the following special matrix classes.                        **
**                                                                      **
**       _matrix: dummy matrix                                          **
**      n_matrix: normal matrix                                         **
**      d_matrix: diagonal matrix                                       **
**      i_matrix: identity matrix                                       **
**      h_matrix: Hermitian matrix                                      **
**                                                                      **
*************************************************************************/

#ifndef   Gmatrix_cc_			// Is file already included?
#  define Gmatrix_cc_ 1			// No, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#  endif

#include <Matrix/matrix.h>		// Include the interface
#include <Matrix/_matrix.h>		// Include base matrix class
#include <Matrix/n_matrix.h>		// Include normal complex arrays
#include <Matrix/d_matrix.h>		// Include diagonal arrays
#include <Matrix/i_matrix.h>		// Include identity arrays
#include <Matrix/h_matrix.h>		// Include Hermitian arrays
#include <Matrix/MxModBas.h>		// Include Matrix module errors
#include <vector>                       // Include libstdc++ STL vectors
#include <iostream>			// Include libstdc++ io streams
#include <fstream>			// Include libstdc++ file streams
#include <sstream>

bool    matrix::FFTComp=true; 		// FFT compatibility (Brigham)
bool    matrix::BlkDiag=true; 		// Block diag. (Sure, Why Not)
MxPrint matrix::PrntFlgs(true,		// Set to print matrix header
                         true,		// Set to print real,imag,complex
                         false,		// Set to print only stored elements
                         30,		// Switch to pictorial at dim=30
                         false,		// Do not qualify pictorial output
                         4,		// Print 4 cols for row vector output
                         30);		// Print 30 rows for col vector output

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                    CLASS MATRIX ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   mx      : A marix (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : string in message
        Output                  none    : Error message output
                                          Program execution stopped if fatal

  The following error messages use the defaults set in the Gutils package

                Case                          Error Message                

                (0)                     Program Aborting.....
                (9)                     Problems During Construction
                default                 Unknown Error                                      

                (1)                     Problems With File PNAME
                (2)                     Cannot Read Parameter PNAME
                (3)                     Invalid Use Of Function PNAME
                default                 Unknown Error - PNAME                */

void matrix::Mxerror(int eidx, int noret) const 
  {
  std::string hdr("Matrix");
  std::string msg;
  switch (eidx)
    {
    case 6:  MxModError(hdr,"Bad Internal Component Access",noret);break;//(6)
    case 7:  MxModError(hdr,"Accessing Empty Operator",noret);    break;// (7)
    case 20: MxModError(hdr,"Unable To Perform Addition",noret);  break;// (20)
    case 28: MxModError(hdr,"Non-Square Diagonalization",noret);  break;// (28)
    case 30: MxModError(hdr,"Failed Requested Operation",noret);  break;// (30)
    case 40: MxModError(hdr,"FFT Still on Vectors Only", noret);  break;// (40)
    case 41: MxModError(hdr,"Base 2 FFT, Needs 2^N Pts", noret);  break;// (41)
    case 42: MxModError(hdr,"Base 2 IFFT, Needs 2^N Pts", noret); break;// (42)
    case 49: MxModError(hdr,"Cannot Set Hermitian Type",  noret); break;// (49)
    case 50: MxModError(hdr,"Cannot Set Matrix Type",     noret); break;// (50)
    case 51: MxModError(hdr,"Matrix Dimension Trouble",   noret); break;// (51)
    case 52: MxModError(hdr,"Cant Read From Binary File", noret); break;// (52)
    case 53: MxModError(hdr,"Element Size Mismatch",      noret); break;// (53)
    default: MxModError(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  }  

void matrix::Mxerror(int eidx, const std::string& pname, int noret) const
  {
  std::string hdr("Matrix");
  std::string msg;
  switch(eidx)
    {
    case 5: msg = std::string("Bad Use Of ") + pname + std::string(" Function ");
             MxModError(hdr,msg,noret);  break;                  // (5)
    case 6: msg = std::string("Problems Using ") + pname + std::string(" Function ");
             MxModError(hdr,msg,noret);  break;                  // (6)
    case 9 :msg = std::string("Construction With ") + pname + std::string("OpReps");
             MxModError(hdr,msg,noret);  break;                  // (9)
    case 52: msg = std::string("Dimension Is ") + pname + std::string("?");
             MxModError(hdr,msg,noret);  break;  
    default: MxModError(hdr, eidx, pname, noret); break;// Unknown Error  (-1)
    }
  }  

volatile void matrix::Mxfatality(int eidx) const
  {
  Mxerror(eidx, 1); 				// Normal non-fatal error
  if(eidx) Mxerror(0);				// Program aborting error
  MxModFatal();                                 // Keep screen nice, exit
  }
     
volatile void matrix::Mxfatality(int eidx, const std::string& pname) const
  {                                                                 
  Mxerror(eidx, pname, 1);			// Normal non-fatal error
  if(eidx) Mxerror(0); 				// Program aborting error
  MxModFatal();                                 // Keep screen nice, exit
  }

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
** matrix type classes (_matrix derived: i_matrix, d_matrix, h_matrix,....)  **
** do not deal with such issues, expecting this generic class to handle it.  **
** With that in mind one can avoid the possibility of having memory leaks.   **
**                                                                           **
******************************************************************************/

_matrix* virtual_copy(_matrix* mx)          { mx->references()++; return mx; }
void     virtual_delete(_matrix *mx)
                  { mx->references()--; if(mx->references() <= 0) delete mx; }

	// Input                mx   : Pointer to _matrix
	// Output               mx1  : Pointer to _matrix
	// Note                      : This makes a true copy of mx
	//			          (via mx->copy) to mx1, making
	//			           a new copy of the matrix data
	// Note                      : USE WHEN _matrix DATA TO BE ALTERED

_matrix* virtual_to_real_copy(_matrix *mx)
  {
  if(mx->references() > 1)         // If multiple references to mx
    {
    _matrix* mx1 =                 // Copy the data to mx1 and
       virtual_copy(mx->copy());   // Increment the reference count to mx1
    virtual_delete(mx);		// Decrement references to mx
    return mx1;			// Return the new matrix pointer, mx1
    }
  else return mx;                  // If 1 reference then just return
  }

matrix::matrix(_matrix* mx) { m = virtual_copy(mx); }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A             CLASS MATRIX CONSTRUCTORS AND DESTRUCTOR
// ____________________________________________________________________________

/*    Arguments                         Constructed Array
    -------------       -------------------------------------------------------
         -               An empty array
     nr,nc,mt,ht         An (nr x nc) array of type mt & hermitian type ht
    nr, nc, z, mt        An (nr x nc) array with all <i|mx|j> = z of type mt
         mx              A duplicate of array mx

        Possible nt Values                        Possible ht Values
  =================================   ========================================
  n_matrix_type == Normal (default)   _hermitian == Hermitian
  h_matrix_type == Hermitian          non_hermitian == Not Hermitian (default)
  d_matrix_type == Diagonal
  i_matrix_type == Identity                                                 

  Note that since an empty matrix will ALWAYS be present (if any matrices are
  are present), the reference count of any empty array will alway be the 
  number of empty arrays in a program + 1. See 1st constructor dummymatrix. */

matrix::matrix( )
  { 
  static _matrix* dummymatrix = NULL;		// Keep the first NULL matrix
  if(dummymatrix == NULL)                        // All sucessive NULL matrices
    dummymatrix = virtual_copy(new _matrix());	// will reference the 1st
  m = virtual_copy(dummymatrix);                 // Now set mx matrix to NULL
  }

matrix::matrix(int i)
  { m = virtual_copy(new n_matrix(i,1)); }

/*
matrix::matrix(int i, int j, matrix_type t, hermitian_type h)
  {
  switch(t) 
    {
    case d_matrix_type: m = virtual_copy(new d_matrix(i,j)); break;
    case i_matrix_type: m = virtual_copy(new i_matrix(i,j)); break;
    case h_matrix_type: m = virtual_copy(new h_matrix(i,j)); break;
    case n_matrix_type:
    default:            m = virtual_copy(new n_matrix(i,j));
    }
  return;
  hermitian_type x; x = h;		  // Compiler likes h to be used
  }
*/
matrix::matrix(int i, int j, const complex& z, matrix_type t, hermitian_type h)
  {
  switch(t)
    {
    case d_matrix_type: m = virtual_copy(new d_matrix(i,j,z));break;
    case i_matrix_type: m = virtual_copy(new i_matrix(i,j,z));break;
    case h_matrix_type: m = virtual_copy(new h_matrix(i,j,z));break;
    default:            m = virtual_copy(new n_matrix(i,j,z));
    }
  /*
  return;
  hermitian_type x;  x = h;		       // Compiler likes h to be used
  */
  }

matrix::matrix(int i, int j, const complex& z, matrix_type t)
  {
  switch(t)
    {
    case d_matrix_type: m = virtual_copy(new d_matrix(i,j,z));break;
    case i_matrix_type: m = virtual_copy(new i_matrix(i,j,z));break;
    case h_matrix_type: m = virtual_copy(new h_matrix(i,j,z));break;
    default:            m = virtual_copy(new n_matrix(i,j,z));
    }
  /*
  return;
  hermitian_type x;  x = h;		       // Compiler likes h to be used
  */
  }

matrix::matrix(int i, int j, double d, matrix_type t, hermitian_type h)
  {
  switch(t)
    {
    case d_matrix_type: m = virtual_copy(new d_matrix(i,j,d));break;
    case i_matrix_type: m = virtual_copy(new i_matrix(i,j,d));break;
    case h_matrix_type: m = virtual_copy(new h_matrix(i,j,d));break;
    default:            m = virtual_copy(new n_matrix(i,j,d));
    }
  }

matrix::matrix(int i, int j, double d, matrix_type t)
  {
  switch(t)
    {
    case d_matrix_type: m = virtual_copy(new d_matrix(i,j,d));break;
    case i_matrix_type: m = virtual_copy(new i_matrix(i,j,d));break;
    case h_matrix_type: m = virtual_copy(new h_matrix(i,j,d));break;
    default:            m = virtual_copy(new n_matrix(i,j,d));
    }
  }

matrix::matrix(int i, int j)
  { m = virtual_copy(new n_matrix(i,j,0)); }

matrix::matrix(int i, int j, const complex& z)
  { m = virtual_copy(new n_matrix(i,j,z)); } 

matrix::matrix(int i, int j,double d)
  { m = virtual_copy(new n_matrix(i,j,d)); } 

matrix::matrix(int i, int j, matrix_type t)
  {
  switch(t) 
    {
    case d_matrix_type: m = virtual_copy(new d_matrix(i,j)); break;
    case i_matrix_type: m = virtual_copy(new i_matrix(i,j)); break;
    case h_matrix_type: m = virtual_copy(new h_matrix(i,j)); break;
    case n_matrix_type:
    default:            m = virtual_copy(new n_matrix(i,j));
    }
  }

matrix::matrix(int i, int j, matrix_type t, hermitian_type h)
  {
  switch(t) 
    {
    case d_matrix_type: m = virtual_copy(new d_matrix(i,j)); break;
    case i_matrix_type: m = virtual_copy(new i_matrix(i,j)); break;
    case h_matrix_type: m = virtual_copy(new h_matrix(i,j)); break;
    case n_matrix_type:
    default:            m = virtual_copy(new n_matrix(i,j));
    }
  }

matrix::matrix(const matrix& mx) { m = virtual_copy(mx.m); }
matrix::~matrix()                { virtual_delete(m); }


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

matrix& matrix::operator= (const matrix& mx)
  {
  if(this == &mx) return *this;		// Nothing if self-assignment
  virtual_delete(m);			// Delete the matrix of this
  m = virtual_copy(mx.m);		// Reference this to mx matrix
  return *this;
  }

//const complex& matrix::operator() const (int i, int j) { return (*m)(i,j); }
complex& matrix::operator() (int i, int j) { return (*m)(i,j); }
complex& matrix::elem(int i, int j)        { return (*m)(i,j); }
//complex  matrix::get(int i, int j)   const { return m->get(i,j); }
//double   matrix::getRe(int i, int j) const { return (m->get(i,j)).Rec(); }
//double   matrix::getIm(int i, int j) const { return (m->get(i,j)).Imc(); }		
complex  matrix::get(int i, int j)   const { if (i>=rows() || j>=cols()) Mxerror(999,1); return m->get(i,j); }
double   matrix::getRe(int i, int j) const { if (i>=rows() || j>=cols()) Mxerror(999,1); return (m->get(i,j)).Rec(); }
double   matrix::getIm(int i, int j) const { if (i>=rows() || j>=cols()) Mxerror(999,1); return (m->get(i,j)).Imc(); }		

void matrix::put(const complex& z, int i, int j)
  {
  m = virtual_to_real_copy(m);		 // Make a new copy (this changes)
  if(!m->put(z,i,j))                       // If the _matrix class put fails
    {                                      // then the type must change
    switch((*m).stored_type())		 // depending on the starting type
      {
      case i_matrix_type:                   // Put value into identity matrix
        if(i == j) set_type(d_matrix_type); // Set z diagonal, imx->dmx
        else       set_type(n_matrix_type); // Set z offdiag,  imx->nmx
        break;

      case d_matrix_type:		        // Put value into diagonal matrix
      case h_matrix_type:                  // Put value into Hermitian matrix
                   set_type(n_matrix_type);// dmx->nmx or hmx->nmx
        break;
      case n_matrix_type:                  // Put value into normal matrix
        Mxerror(6, "put", 1);		 // Trouble with put function!
        Mxfatality(30);                    // Failed requested operation
        break;
      default:                             // Put value into generic matrix
        set_type(n_matrix_type);           // Switch it to a normal matrix 
        break;
      }
    if(!m->put(z,i,j))                     // Now use put with new matrix type
      {
      Mxerror(6, "put", 1);                // Trouble with put function!
      Mxfatality(30);                      // Failed requested operation
      }
    }
  }

void matrix::put_h(const complex& z, int i, int j)
  {
  m = virtual_to_real_copy(m);				// Make a new copy (mx changes)
  if(!m->put_h(z,i,j))					// See if handled by matrix class
    {
    switch((*m).stored_type())
      {
      case i_matrix_type:				// Value into i_matrix
        if(i==j) set_type(d_matrix_type);		// Put z diagonal, array to diag
        else     set_type(h_matrix_type);		// Else switch to Hermitian
        break;
      case d_matrix_type:				// Value into d_matrix
        if(m->is_real()) set_type(h_matrix_type);	// If dmx real become Hermitian
        else             set_type(n_matrix_type);	// Else switch to normal
        break;
      case h_matrix_type: set_type(n_matrix_type); break;  // Value int h_matrix
      case n_matrix_type: 				// Value into n_matrix should work!
        Mxerror(6, "Hermitian Element Put", 1);		// (i.e. handled in 2nd line above)
        Mxerror(56, 1);   		                // Check matrix indices
        Mxfatality(3, "put_h");				// This trouble in put_h function
        break;
      default:            set_type(n_matrix_type); break;  // Put value into generic matrix
      }
    if(!m->put_h(z,i,j))				// Give it another try
      {
      Mxerror(6, "Hermitian Element Put", 1);		// (i.e. handled in 2nd line above)
      Mxerror(56, 1);   		                // Check matrix indices
      Mxfatality(3, "put_h");				// This trouble in put_h function
      }							// Die if this doesn't work!
    }
  }


matrix matrix::get_block(int row, int col, int nrows, int ncols)
  { return matrix(m->get_block(row,col,nrows,ncols)); }

void matrix::put_block(int row, int col, const matrix& mx)
  {
  if(!row && !col)			// It is possible to avoid construction of
    {					// a new matrix if put starts at <0|mx|0> 
    if(this == &mx) return;		// One way is see if putting mx into itself
    else if(mx.rows()==rows() &&	// Another is if put mx replaces all of mx1 
                   mx.cols() == cols())
      {
      virtual_delete(m);		// Delete the matrix of mx1
      m = virtual_copy(mx.m);		// Reference matrix of mx1 to mx matrix
      return;
      }
    }
  if(stored_type() == i_matrix_type)	// Another way the matrix may not change is
    {
    if(mx.stored_type()==i_matrix_type)	// if an I matrix s put on an I matrix diagonal
      if(row == col)			// If so and the dimensions are O.K., the
        if(row+mx.rows() <= rows())	// matrix is unaffected by the block put
          return;
     }
  m = virtual_to_real_copy(m);		// Make a real copy, mx1 will change
  if(!m->put_block(row, col, mx.m))	// Hope _matrix class can do it
    {					// if not, type will have to change
    switch((*m).stored_type())
      {
      case i_matrix_type:		// Put block into identity matrix
        if(row == col)			// See if on the diagonal
          {
          if(mx.stored_type()		// If put diagonal block on diagonal
                      == d_matrix_type) // then identity becomes diagonal
            set_type(d_matrix_type);
          else if(mx.stored_type()	// If put Hermitian block on diagonal
                      == h_matrix_type) // then identity becomes Hermitian
            set_type(h_matrix_type);
          else				// For anything else, change identity
            set_type(n_matrix_type);	// to a normal matrix
          }
        else				// If put block off diagonal, identity
          set_type(n_matrix_type);	// matrix must become a normal matrix
        break;
      case d_matrix_type:		// Put block into diagonal matrix
        if(row == col && is_real())	// See if put on the diagonal and 
          {
          if(mx.stored_type()		// If put Hermitian block on diagonal
                      == h_matrix_type) // then a real diagonal becomes Hermitian
            set_type(h_matrix_type);
          else				// If put Non-Hermitian block, change diagonal
            set_type(n_matrix_type);	// to a normal matrix
          }
        else				// For anything else, change diagonal
          set_type(n_matrix_type);	// to a normal matrix
        break;
      default:				// Put value into generic 
        set_type(n_matrix_type);
        break;
      }
    if(!m->put_block(row, col, mx.m))	// Try put_block again
      {
      Mxerror(6, "Put Block In Matrix", 1);
      Mxfatality(3, "put_block");	// This trouble in put_block function
      }
    }
  }
      
// ____________________________________________________________________________
// C                CLASS MATRIX HERMITIAN_TYPE HANDLING
// ____________________________________________________________________________

/*   Function      Output                       Description
 ----------------  ------  -------------------------------------------------
 stored_hermitian  h_type  Returns stored hermitian type
 check_hermitian   h_type  Returns stored hermitian type (with mx(ij)<d = 0)
 set_hermitian     void    Sets matrix hermitian, may result in data loss
 test_hermitian    h_type  Returns possible hermitian type (w/ mx(ij)<d = 0) 

 Note that the default value of d (used for roundoff) is GMxCut which can 
 be set by the user.                                                         */

hermitian_type matrix::stored_hermitian( ) const 
  { return m->stored_hermitian(); }

hermitian_type matrix::test_hermitian(double d) const
  { return m->test_hermitian(d); }

// sosi - the hermitian_type isn't even used!
// sosi - should not copy here if input is hermitian
// why do not individual classes handle this change
// why isn't the unknown matrix type handled?

void matrix::set_hermitian(hermitian_type h)
  {
  m = virtual_to_real_copy(m);            // Make copy which we'll alter
  switch((*m).stored_type())
    {
    case d_matrix_type:			// Set diagonal to Hermitian
      {					// Just zero the imaginary parts
      for(int i=0; i<(*m).rows(); i++)	// on the diagonal & its Hermitian
	set_imaginary_part((*m)(i,i),0);
      break;
      }
    case n_matrix_type: 		// Set normal matrix to Hermitian 
      {					// Copy <i|mx|j> to <j|mx|i>* for i>j
      int i;
      for(i=0; i<(*m).rows(); i++)	// ant set <i|mx|i> to real
	for(int j=0; j<i; j++)
	  (*m)(i,j) = ((*m)(j,i)).conj();
      for(i=0; i<(*m).rows(); i++) 
	set_imaginary_part((*m)(i,i),0);
      break;
      }
    case i_matrix_type:			// Set identity or Hermitian to Hermitian
    case h_matrix_type:			// does nothing, they are Hermitian 
      break;
    default:				// Set unknown type Hermitian, don't know
      Mxerror(49, 1);			//    Cannot set  hermitian matrix type
      Mxerror(6, "set_hermitian", 1);	//    Trouble w/ set_hermitian fnct!
      Mxfatality(30);			//    Failed operation
    }
  return;
  hermitian_type x;
  x = h;				// Compiler likes h to be used
  }
      
hermitian_type matrix::check_hermitian(double d)
  {
  hermitian_type h = test_hermitian(d);	// First test if it can be Hermitian
  set_hermitian(h);                       // Set the type accordingly
  return h;                               // Return T/F depending on test
  }

// ____________________________________________________________________________
// D                    CLASS MATRIX MATRIX_TYPE HANDLING
// ____________________________________________________________________________

/*   Function    Output                       Description
   -----------  -------  ------------------------------------------------------
   stored_type  mx_type  Returns current matrix type, no checking is done
   test_type    mx_type  Returns mx_type if array can be that type to within d
                         This function will return its own type in test fails
   set_type     void     Set matrix to specified type, may result in data loss
			 This does matrix copying of (valid) elements
   check_type   mx_type  Test if mx could be matrix type (with mx(ij)<d = 0)
                         Then set the type to that specified if possible
   mxtype       string   Return a string for the current matrix type       

    Note: test_type will return tested type if true, its own type if false
    Note: check_type will convert the arrray if no data loss will occur.
          If it is not possible, the current matrix type is returned.         */

matrix_type matrix::stored_type() const             { return m->stored_type(); }
matrix_type matrix::test_type(matrix_type t, double d) const
                                                   { return m->test_type(t,d); }
void matrix::set_type(matrix_type t)
  { 
  if(m->stored_type() == t) return;		// Nothing if no type change
  _matrix *mx = NULL;				// Pointer to converted array
  switch(t)					// Convert *m into new array type
    {
    case n_matrix_type: mx = m->NMX(); break;	// These will return a pointer to
    case d_matrix_type: mx = m->DMX(); break;	// to a new (unreferenced) array of
    case i_matrix_type: mx = m->IMX(); break;	// the proper type, leaving the 
    case h_matrix_type: mx = m->HMX(); break;	// array m points to untouched
    default: 
      Mxerror(50, 1);				//    Cannot set matrix type
      Mxerror(6, "set_type", 1);		//    Trouble w/ set_type fnct!
      Mxfatality(30);				//    Failed operation
    }
  virtual_delete(m);				// Delete one reference to m
  m = virtual_copy(mx);				// Now set m to mx, (mx refs++)
  }

matrix_type matrix::check_type(const matrix_type t, const double d)
               { matrix_type tmp = test_type(t,d); set_type(tmp); return tmp; }
std::string matrix::mxtype() const                           { return m->mxtype(); }

// ____________________________________________________________________________
// E              CLASS MATRIX VARIOUS MATRIX CHECKS & PARAMETERS
// ____________________________________________________________________________

/*  Function  Output                       Description
    --------  ------  ---------------------------------------------------------
      rows     int    Returns # of matrix rows       - handled by _matrix
      cols     int    Returns # of matrix columns    - handled by _matrix
      refs     int    Returns # of matrix references - handled by _matrix
      pts      int    Returns # of matrix elements   - handled by _matrix    */

int matrix::rows( ) const { return m->rows(); }
int matrix::cols( ) const { return m->cols(); }
int matrix::refs( ) const { return m->refs(); }
int matrix::pts( )  const { return m->pts();  }

/*   Function    Output                       Description
   ------------  ------  ------------------------------------------------------
   is_symmetric   bool   TF if symmetric within d (default GMxCut)
   is_hermitian   bool   TF if hermitian within d (default GMxCut)
   is_unitary     bool   TF if unitary   within d (default GMxCut), much CPU
   is_real        bool   TF if real      within d (default GMxCut)
   is_imaginary   bool   TF if imaginary within d (default GMxCut)
   is_complex     bool   TF if complex   within d (default GMxCut)
   is_zero        bool   TF if zero      within d (default GMxCut) 
   is_diagonal    bool   TF if diagonal  within d (default GMxCut)
   is_square      bool   TF if square                                        */

bool matrix::is_symmetric(const double d) const { return m->is_symmetric(d); }
bool matrix::is_hermitian(const double d) const { return m->is_hermitian(d); }
bool matrix::is_unitary(const   double d) const { return m->is_unitary(d);   }
bool matrix::is_real(const      double d) const { return m->is_real(d);      }
bool matrix::is_imaginary(const double d) const { return m->is_imaginary(d); }
bool matrix::is_complex(const   double d) const { return m->is_complex(d);   }
bool matrix::is_zero(const      double d) const { return m->is_zero(d);      }
bool matrix::is_diagonal(const  double d) const { return m->is_diagonal(d);  }
bool matrix::is_square()                  const { return m->is_square();     }

// ____________________________________________________________________________
// F                 CLASS MATRIX BINARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________
 
/* These functions allow for simple mathematical operations between two
   matrices and between a matrix and a scalar. All of these operations are 
   done in the _matrix derived classes for faster throughput.  Note that the
   matrix reference counting is done by the constructor used in the return.

   Operator Arguments      Result        Operator Arguments       Result
   -------- --------- -----------------  -------- --------- -------------------
      +      mx1,mx2      mx1+mx2           *       mx,z           z*mx
      -      mx1,mx2      mx1-mx2           *       z,mx           z*mx
      *      mx1,mx2      mx1*mx2           *       mx,d           d*mx
      /      mx1,mx2      mx1*inv(mx2)      *       d,mx           d*mx
      /      mx,z         (1/z)*mx          *       mx,d         (1/d)*mx
      /      z,mx         z*inv(mx)         *       d,mx         d*inv(mx)  */

matrix matrix::operator+ (const matrix& mx) const
                                 { return matrix(m->add(mx.m)); }
matrix matrix::operator- (const matrix& mx) const
                                 { return matrix(m->subtract(mx.m)); }
matrix matrix::operator* (const matrix& mx) const
                                 { return matrix(m->multiply(mx.m)); }
matrix matrix::operator* (const complex& z) const
                                 { return matrix(m->multiply(z)); }
matrix matrix::operator* (      double d) const
                                 { return matrix(m->multiply(complex(d))); }
matrix matrix::operator/ (const matrix& mx) const
                                 { return matrix(m->divide(mx.m)); }
matrix matrix::operator/ (const complex& z) const
                                 { return matrix(m->divide(z)); }
matrix matrix::operator/ (        double d) const
                                 { return matrix(m->divide(complex(d))); }
matrix operator* (const complex& z, const matrix& mx)
                                 { return matrix(mx.m->multiply(z)); }
matrix operator* (        double d,   const matrix& mx)
                                 { return matrix(mx.m->multiply(complex(d))); }

// ____________________________________________________________________________
// G                 CLASS MATRIX UNARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________
 

// sosi - what happens in these next two when this is NULL?

matrix& matrix::operator += (const matrix& mx)	// Implementation of this += mx

  {
  if(!mx.rows()) return *this;			// Do nothing if mx is NULL
  if(!rows())					// If we are empty the result
    {						// is just mx
    virtual_delete(m);
    m = virtual_copy(mx.m);
    return *this;
    }
  m = virtual_to_real_copy(m);		// First make a real copy of this
  _matrix* tmp = mx.m->add_two(m);	// Route to _matrix class of mx!
  if(tmp != m)
    {
    virtual_delete(m);
    m = virtual_copy(tmp);
    }
  return *this;
  }

matrix& matrix::operator -= (const matrix& mx)	// Implement this -= mx

  {
  if(!mx.rows()) return *this;			// Do nothing if mx is NULL
  if(!rows())					// If we are empty the result
    {						// is just mx
    virtual_delete(m);
    m = virtual_copy(mx.m->negate());
    return *this;
    }
  m = virtual_to_real_copy(m);		// First make a real copy of this
  _matrix* tmp = mx.m->subtract_two(m);	// Route to _matrix class of mx!
  if(tmp != m)
    {
    virtual_delete(m);
    m = virtual_copy(tmp);
    }
  return *this;
  }


matrix& matrix::operator *= (const matrix& mx)	// Implements this *= mx

	// Input                mx0  : Input matrix (this)
	// 			mx   : Second matrix
	// Output		mx0  : Matrix mx0 is modified by multiplication
	//			       into mx

  {
  if(mx.m->stored_type()==i_matrix_type)// No change if mx is i_matrix
    if(cols() == mx.cols())
      return *this;
  m = virtual_to_real_copy(m);		// First make a real copy of this
  _matrix* tmp = mx.m->multiply_two(m);	// Route to _matrix class of mx!
  if(tmp != m)
    {
    virtual_delete(m);
    m = virtual_copy(tmp);
    }
  return *this;
  }

  
	// Input                mx   : Input matrix (this)
	// 			z    : Complex number
	// Output		mx   : Matrix mx is scaled by z
	//			       mx *= z

matrix& matrix::operator *= (const complex& z)	// Implements this *= z
  {
  if(z == complex1) return *this;	// If z is 1, no change in mx
  m = virtual_to_real_copy(m);		// First make a real copy of this
  _matrix* tmp = m->multiply_two(z);	// Route to _matrix class of this
  if(tmp != m)
    {
    virtual_delete(m);
    m = virtual_copy(tmp);
    }
  return *this;
  }


matrix& matrix::operator *= (double d)	// Implements this *= d
  
	// Input                mx   : Input matrix (this)
	// 			d    : A real number
	// Output		mx   : Matrix mx is scaled by d
	//			       mx *= d

  {
  if(d==1) return *this; 			// If d is 1, no change in mx
  m = virtual_to_real_copy(m);			// First make a real copy of this
  _matrix* tmp = m->multiply_two(complex(d));	// Route to _matrix class of this
  if(tmp != m)
    {
    virtual_delete(m);
    m = virtual_copy(tmp);
    }
  return *this;
  }


matrix& matrix::operator /= (const matrix& mx)	// Implements this /= mx

	// Input                mx0  : Input matrix (this)
	// 			mx   : Second matrix
	// Output		mx0  : Matrix mx0 is modified by division by
	//			       mx

  {
  if(mx.m->stored_type()==i_matrix_type)// No change if mx is i_matrix
    if(cols() == mx.cols())
      return *this;
  m = virtual_to_real_copy(m);		// First make a real copy of this
  _matrix* tmp = mx.m->divide_two(m);	// Route to _matrix class of mx!
  if(tmp != m)
    {
    virtual_delete(m);			// Dereference the original m
    m = virtual_copy(tmp);
    }
  return *this;
  }


matrix& matrix::operator /= (const complex& z)	// Implements this /= z

	// Input                mx   : Input matrix (this)
	// 			z    : Complex number
	// Output		mx   : Matrix mx is scaled by z
	//			       mx *= 1/z
  {
  if(z == complex1) return *this;	// If z is 1, no change in mx
  m = virtual_to_real_copy(m);		// First make a real copy of this
  _matrix* tmp = m->divide_two(z);	// Route to _matrix class of this
  if(tmp != m)
    {
    virtual_delete(m);
    m = virtual_copy(tmp);
    }
  return *this;
  }


matrix& matrix::operator /= (double d)	// Implements this /= d

	// Input                mx   : Input matrix (this)
	// 			d    : A real number
	// Output		mx   : Matrix mx is scaled by d
	//			       mx *= 1/d
  {
  if(d==1) return *this; 			// If d is 1, no change in mx
  m = virtual_to_real_copy(m);			// First make a real copy of this
  _matrix* tmp = m->divide_two(complex(d));	// Route to _matrix class of this
  if(tmp != m)
    {
    virtual_delete(m);
    m = virtual_copy(tmp);
    }
  return *this;
  }

// ____________________________________________________________________________
// H                  CLASS MATRIX SIMPLE UNARY FUNCTIONS
// ____________________________________________________________________________
 
/* 		Input          mx : Matrix (this)
        	Output         mx1: A matrix derived from mx
		Note 	     	  : Operations performed _matrix derived class

    Function/Operator	  Output 		    Description
    ------------------    -------------  --------------------------------------
	operator-          mx1 = -mx	  negation
           Re		   mx1 = Re{mx}	  zero imaginaries
           Im		   mx1 = IM{mx}   zero reals
          conj	           mx1 = mx* 	  conjugate <i|mx1|j>=<i|mx*|j>
        transpose 	   mx1 = mxt	  transpose <i|mx1|j>=<j|mx|i>
        adjoint	           mx1 = mxt*	  adjoint   <i|mx1|j>=<j|mx*|i>
	 exp               mx1 = exp(mx)  exponentiation 
	 trace		     z = Tr{mx}	  trace     z = sum <i|mx|i>         */

//matrix operator- (const matrix& mx) { return matrix(mx.m->negate());    }
matrix Re(const matrix& mx)         { return matrix(mx.m->RE());        }
matrix Im(const matrix& mx)         { return matrix(mx.m->IM());        }
matrix conj(const matrix& mx)       { return matrix(mx.m->conjugate()); }
matrix transpose(const matrix& mx)  { return matrix(mx.m->transpose()); }
matrix adjoint(const matrix& mx)    { return matrix(mx.m->adjoint());   }
complex trace(const matrix& mx)     { return mx.m->trace();             }

matrix  matrix::operator- () const  { return matrix(m->negate());    }
matrix  matrix::Re()         const  { return matrix(m->RE());        }
matrix  matrix::Im()         const  { return matrix(m->IM());        }
matrix  matrix::conj()       const  { return matrix(m->conjugate()); }
matrix  matrix::transpose()  const  { return matrix(m->transpose()); }
matrix  matrix::adjoint()    const  { return matrix(m->adjoint());   }
matrix  matrix::exp()        const  { return matrix(m->mxexp());     }
complex matrix::trace()      const  { return        m->trace();      }

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

//  sosiz - at the very least you must check the bounds on i,j
matrix matrix::swaprows(int i, int j) { return matrix(m->swaprows(i,j)); }
matrix matrix::swapcols(int i, int j) { return matrix(m->swapcols(i,j)); }
matrix matrix::permute( int i, int j) { return matrix(m->permute(i,j));  }

double  matrix::maxRe() const   { return m->maxRe(); }
double  matrix::maxIm() const   { return m->maxIm(); }
complex matrix::maxZ()  const   { return m->maxZ();  }
double  matrix::minRe() const   { return m->minRe(); }
double  matrix::minIm() const   { return m->minIm(); }
complex matrix::minZ()  const   { return m->minZ();  }

// ____________________________________________________________________________
// I                  CLASS MATRIX SIMPLE BINARY FUNCTIONS
// ____________________________________________________________________________

/* These next functions just call the _matrix class functions. The _matrix
   functions return a pointer to the _matrix containing the result and
   the matrix result is constructed by virtual copy during the return

      Function     Output                     Description
   -------------  ------  ----------------------------------------------------
      trace          z    Trace(mx1*mx2) done efficiently
   adjoint_times    mx    Adjoint(mx1)*mx2 done efficiently
   times_adjoint    mx    mx1*Adjoint(mx2) done efficiently                 */

complex trace(const matrix& mx, const matrix& mx1) {return mx.m->trace(mx1.m);}
complex matrix::trace(const matrix& mx1) const     {return m->trace(mx1.m);}
matrix  adjoint_times(const matrix& mx, const matrix& mx1)
  { return matrix(mx.m->adjoint_times(mx1.m)); }
matrix times_adjoint(const matrix& mx, const matrix& mx1)
  { return matrix(mx.m->times_adjoint(mx1.m)); }
void enable_blockdiag() {matrix::BlkDiag=true;}
void disable_blockdiag() {matrix::BlkDiag=false;}

// ____________________________________________________________________________
// J               CLASS MATRIX COMPLEX UNARY FUNCTIONS
// ____________________________________________________________________________

complex         det(const matrix& mx)  { return mx.m->det();  }
complex matrix::det() const            { return m->det();     }
int             rank(const matrix& mx) { return mx.m->rank(); }


// ******************************      FFT      *******************************

/* These functions handle Fast Fourier Transforms (FFT) and inverse
   Fast Fourier Transforms (IFFT).

	 Input	      mx  : Matrix to be transformed
	 Output      mx1 : Transformed matrix
	 Note	          : The static matrix flag FFTComp sets whether
                          the transform is "folded" or not    n    n
	 Note            : Matrix dimension must be either (1x2 or 2 x1)     */

//int IsBase2(double dim)
//  {
//  while(dim > 1) dim /= 2;
//  (dim == 1)?return 1:return 0;
//  }

matrix FFT(const matrix& mx) { return mx.FFT(); }
matrix IFFT(const matrix& mx) { return mx.IFFT(); }

matrix matrix::FFT() const
  {
  matrix mx1(rows(),cols(),n_matrix_type);
  double dim=0;					// Check that matrix is either
  int TF = 1;					// (1 x 2**n) or (2**n x 1)
  if(rows() == 1)      dim=double(cols());
  else if(cols() == 1) dim=double(rows());
  else
    {
    Mxerror(5,"FFT", 1);			// Bad use of FFT function
    Mxfatality(40);				// Only allowed on vectors
    }
  if(dim != 1 && TF)
    {
    while(dim > 1) dim /= 2;
    if(dim != 1)   TF = 0;
    }
  else TF = 0;
  if(!TF) 				// If here then we have not used
    {					// a base 2 dimension for FFT!
    Mxerror(5,"FFT", 1);		// Bad use of FFT function
    Mxfatality(41);			// Need base 2 # of points
    }
  m->convert(mx1.m);			// Copy data to new n_matrix
  ((n_matrix*)(mx1.m))->FFT(1,FFTComp);	// Use the n_matrix FFT algorithm
  return mx1;				// Return mx1;
  }

matrix matrix::IFFT() const
  {
  matrix mx1(rows(),cols(),n_matrix_type);	// Output array
  double dim=0;					// Dimension in 1D
  bool basetwo = true;				// Assume dimension OK	
  if(rows() == 1)      dim = double(cols());	// See if col. dim OK
  else if(cols() == 1) dim = double(rows());	// See if row dim OK
  else basetwo = false;				// Else not base 2 FFT
  if(dim != 1 && basetwo) 			// Check that dim is
    { 						// of base 2 (our FFT)
    while(dim > 1) dim /= 2; 			// Need dim of 1 x 2^n
    if(dim != 1) basetwo = false; 		// or 2^n x 1
    }
  else basetwo = false;
  if(!basetwo)
    {
    Mxerror(5,"IFFT", 1);			// Bad use of FFT function
    Mxfatality(40);				// Only allowed on vectors
    }
  m->convert(mx1.m);				// Copy data to new n_matrix
  ((n_matrix*)(mx1.m))->FFT(-1,FFTComp);	// Use the n_matrix FFT algorithm
  mx1 /= complex(mx1.rows()*mx1.cols());	// Scale back by matrix size
  return mx1;					// Return mx1;
  }

// ____________________________________________________________________________
// K            CLASS MATRIX COMPLEX BINARY (&TRINARY) FUNCTIONS
// ____________________________________________________________________________
 
 
// ***************************** tensor product *******************************
 

matrix tensor_product(const matrix& mx1, const matrix& mx2)
     { return matrix(mx1.m->tensor_product(mx2.m)); }

        // Input            mx1 : A matrix
        //                  mx2 : A second matrix
        // Output           pdt : A matrix which is the tensor product
        //                        of the two input matrices

        //                           pdt        =   mx1 (x) mx2

        //                       (m*o x n*p)       (mxn)   (oxp)

        //                    <i*o+k|pdt|j*p+l> = <i|mx1|j><k|mx2|l>

// ____________________________________________________________________________
// L                       CLASS MATRIX I/O FUNCTIONS
// ____________________________________________________________________________

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

void    matrix::Header(bool   hf) {          PrntFlgs.MxHdr    = hf; }
void    matrix::PrintRI(bool  pi) {          PrntFlgs.MxRIPrnt = pi; }
void    matrix::PrintAll(bool pa) {          PrntFlgs.MxAll    = pa; }
void    matrix::PictDim(int   pd) { if(pd>1) PrntFlgs.MxPctDim = pd; }
void    matrix::PrintVal(bool pv) {          PrntFlgs.MxPctVal = pv; }
void    matrix::PrintCols(int cl) { if(cl>0) PrntFlgs.VxCols   = cl; }
void    matrix::PrintRows(int rl) { if(rl>0) PrntFlgs.VxRows   = rl; }
MxPrint matrix::PrintFlags()      { return   PrntFlgs; }

//-----------------------------------------------------------------------------
//                          ASCII OUTPUT FUNCTIONS
//-----------------------------------------------------------------------------

/* These functions allow users to write matrices in formatted ASCII to an 
   output stream.  They also allow users to obtain strings representing 
   array rows which allows for full flexibility in array output.
 
    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
     print     ostream     Writes array elements in ASCII to output stream
      <<       ostream     Writes array elementes in ASCII to output stream
    picture    ostream     Writes array pictorially to output stream         */  


std::ostream& matrix::printHdr(std::ostream& ostr) const
  {
  std::string hdr = MxModdec(rows()) + std::string(" x ") + MxModdec(cols());
  hdr += std::string(" ") + m->mxtype(PrntFlgs.MxRIPrnt);
  hdr += std::string(" Matrix");
  int len = 40-hdr.length()/2;
  if(len>0) ostr << std::string(len, ' ');
  ostr << hdr;
  return ostr;
  }

std::ostream& matrix::print(std::ostream& ostr) const
  {
  if(PrntFlgs.MxHdr) { printHdr(ostr); ostr << "\n\n"; }
  if(rows() > PrntFlgs.MxPctDim)
    m->picture(ostr,PrntFlgs);
  else
    m->print(ostr, PrntFlgs);
  return ostr;
  }

std::ostream& matrix::picture(std::ostream& ostr) const
  {
  if(PrntFlgs.MxHdr) { printHdr(ostr); ostr << "\n\n"; }
  m->picture(ostr, PrntFlgs); 
  return ostr;
  }

std::ostream& operator << (std::ostream& ostr, const matrix& mx)
  { mx.print(ostr); return ostr; }
 
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

std::ofstream& matrix::write(std::ofstream &fp, int form) const
  {
  if(!form)
    {
    char c;				// Character to flag matrix type
    int ro = rows();                    // Number of matrix rows
    int co = cols();                    // Number of matrix columns
    switch(stored_type())
      {
      case i_matrix_type:		// Write an identity matrix
        c = 'i';			// Flag that it is identity
        fp.write((char*)&c, sizeof(char));
        fp.write((char*)&ro, sizeof(int));	// Write the number of rows
        break;
      case d_matrix_type:		// Write a diagonal matrix
        c = 'd';			// Flag that it is diagonal
        fp.write((char*)&c, sizeof(char));
        c = 'c';			// Flag that is is complex
        fp.write((char*)&c, sizeof(char));
        fp.write((char*)&ro, sizeof(int));	// Write the number of rows
        break;
      case h_matrix_type:		// Write a Hermitian matrix
        c = 'h';			// Flag that it is Hermitian
        fp.write((char*)&c, sizeof(char));
        fp.write((char*)&ro, sizeof(int));	// Write the number of rows
        break;
      case n_matrix_type:		// Write a normal matrix
        c = 'n';			// Flag that it is normal
        fp.write(&c, sizeof(char));
        c = 'c';			// Flag that is is complex
        fp.write((char*)&c, sizeof(char));
        fp.write((char*)&ro, sizeof(int));	// Write the number of rows
        fp.write((char*)&co, sizeof(int));	// Write the number of columns
        break;
      default: break;
      }
    }
  m->write(fp, form);			// Now use appropriate matrix class
  return fp;				// Return the (advanced) file
  }

std::ofstream& write(std::ofstream &a, const matrix& mx)
  { mx.m->write(a,1); return a; }
 
// --------------------------- Binary Input Functions -------------------------

std::ifstream& matrix::read(std::ifstream& fp)
  {
  int fploc = fp.tellg();			// Get the current file location
  char c1, c2;					// Characters for matrix type
  int r, c;					// Read matrix rows and columns
  int bigmx = 100000000;			// Used to check mx dimensions
  fp.read(&c1, sizeof(char));			// Attempt to read matrix type
  switch(c1)					// (Done if GAMMA binary format)
    {
    case 'n':					// GAMMA n_matrix
      fp.read((char*)&c2, sizeof(char));		// Continue to read matrix type
      fp.read((char*)&r, sizeof(int));			// Read the number of rows
      fp.read((char*)&c, sizeof(int));			// Read the number of columns 
      if(c2=='c' && r>0 && c>0 && r*c<bigmx)	// Check that all is reasonable 
        {
        virtual_delete(m);			// Delete the actual mx matrix
        m = virtual_copy(new n_matrix(r,c));	// Construct a new n_matrix
        m->read(fp);				// Read in the n_matrix
        return fp;				// Now stop here
        }
      break;
    case 'h':					// GAMMA h_matrix
      fp.read((char*)&r, sizeof(int));		// Read the number of rows
      if(r>0 && r*r<bigmx)			// Check dimensions reasonable 
        {
        virtual_delete(m);			// Delete the actual mx matrix
        m = virtual_copy(new h_matrix(r,r));	// Construct a new h_matrix
        m->read(fp);				// Read in the h_matrix
        return fp;				// Now stop here
        }
      break;
    case 'd':					// GAMMA d_matrix
      fp.read((char*)&c2, sizeof(char));	// Continue to read matrix type
      fp.read((char*)&r, sizeof(int));		// Read the number of rows
      if(c2=='c' && r>0 && r*r<bigmx)		// Check that all is reasonable
        {
        virtual_delete(m);			// Delete the actual mx matrix
        m = virtual_copy(new d_matrix(r,r));	// Construct a new d_matrix
        m->read(fp);				// Read in the d_matrix
        return fp;				// Now stop here
        }
      break;
    case 'i':					// GAMMA i_matrix
      fp.read((char*)&r, sizeof(int));		// Read the number of rows
      if(r>0 && r*r<bigmx)			// Check dimensions reasonable 
        {
        virtual_delete(m);			// Delete the actual mx matrix
        m = virtual_copy(new i_matrix(r,r));	// Construct a new i_matrix
        m->read(fp);				// Read in the i_matrix
        return fp;				// Now stop here
        }
      break;
    default:					// File doesnt point to GAMMA mx
      r = rows();				// Assume mx has proper rows
      c = cols();				// Assume mx has proper cols
      if(r<0 || c<0 || r*c>bigmx)		// Check that all is reasonable
        {
        Mxerror(50,1);				// Matrix dimension trouble
        std::string serr = MxModdec(r)
                         + std::string(" x ")
                         + MxModdec(c);
        Mxerror(52, serr, 1);
        Mxerror(6, "read", 1);			//    Trouble w/ read fnct!
        Mxfatality(30);				//    Failed operation
        }
      fp.seekg(fploc, std::ios::beg);		// Return to initial file spot
      matrix_type mt = stored_type();		// Store the current matrix type
      virtual_delete(m);			// Delete the actual mx matrix
      m = virtual_copy(new n_matrix(r,c));	// Construct a new n_matrix
      m->read(fp);				// Read in matrix as a n_matrix
      check_type(mt);				// If able, mx to original type
      return fp;				// Now stop here
      break;
    }
  Mxerror(6, "read", 1);			//    Trouble w/ read fnct!
  Mxerror(52,1);				//    Cant read from binary file
  Mxfatality(30);				//    Failed operation
  return fp;
  }

//-----------------------------------------------------------------------------
//                            ASCII INPUT FUNCTIONS
//-----------------------------------------------------------------------------
    
/* These functions perform ASCII input of arrays.
 
    Function  Arguments                           Result
   ---------- ---------   -----------------------------------------------------
   readASC    istream     Read in array from an input stream
   ask          ----      Interactively requests matrix from user            */  
 
std::istream& operator >> (std::istream& istr, matrix& mx)
  {
  virtual_delete(mx.m);				// Delete the actual matrix of mx
  mx.m = virtual_copy(new n_matrix(0,0));	// Set the matrix to a new n_matrix
  mx.m->readASC(istr);				// Read as a n_matrix ONLY
  return istr;					// Return the input stream
  }						// Only n_matrix::readASC used

void matrix::ask(const matrix_type t)
  {
  int r,c;
  switch(t) 
    {
    case _matrix_type:				// Fill up a _matrix
    case n_matrix_type:				// Fill up a n_matrix
    default:					// Fill up unknown matrix type
      std::cout << "\n\tPlease Input the Number of Rows: ";
      std::cin >> r;
      std::cout << "\n\tPlease Input the Number of Columns: ";
      std::cin >> c;
      virtual_delete(m);			// Delete the current matrix of mx
      m = virtual_copy(new n_matrix(r,c));	// Set matrix to new n_matrix 
      m->ask();					// Route to n_matrix::ask
      break;
    case d_matrix_type:
      std::cout << "\n\tPlease Input the Matrix Dimension ";
      std::cin >> r;
      virtual_delete(m);			// Delete the current matrix of mx
      m = virtual_copy(new d_matrix(r,r));	// Set matrix to new n_matrix 
      m->ask();					// Route to d_matrix::ask
      break;
    case i_matrix_type:
      virtual_delete(m);			// Delete the current matrix of mx
      std::cout << "\n\tPlease Input the Matrix Dimension ";
      std::cin >> r;
      m = virtual_copy(new i_matrix(r,r));	// Set matrix to new i_matrix 
      break;					// Thats it!
    case h_matrix_type:
      virtual_delete(m);			// Delete the current matrix of mx
      std::cout << "\n\tPlease Input the Matrix Dimension ";
      std::cin >> r;
      m = virtual_copy(new h_matrix(r,r));	// Set matrix to new h_matrix 
      m->ask();					// Route to h_matrix::ask
      break;
    }
  return;
  }


// ____________________________________________________________________________
// M                   CLASS MATRIX MISCELLANEOUS FUNCTIONS
// ____________________________________________________________________________


matrix matrix::resize(int i, int j)

	// Input                mx0  : Matrix (this)
	// 			i,j  : Row and column dimensions
	// Output               mx   : Matrix mx0 is resized to have the
	//			       dimension (i x j) and this is 
	//			       output as the new matrix mx
	// Note			     : Fails if # elements changes


// GNU COMPILER BUG:  FOR SOME REASON THIS ROUTINE IS INCONSISTENT
//                    IN MATRIX REFERENCE COUNTING!!  IT SEEMS THAT
//		      THE ADDITION OF 1 OUTPUT STATMENT TELLS THE
//		      COMPILER TO WORK CORRECTLY ?!.
// sosi - I haven't checked this again since gcc 2.7.#

  {
  int r = rows();			// Rows of this
  int c=cols();				// Cols of this
  if(r*c != i*j)			// Check that # elements same
    {
    Mxerror(50,1);			// Matrix dimension trouble
    Mxerror(52,1);			// Element size mismatch
    Mxerror(6, "resize", 1);		// Trouble w/ resize fnct!
    Mxfatality(30);			// Failed operation
    matrix mx(*this);
    return mx;
    }
  else
    if(stored_type() == n_matrix_type)	// Normal matrix conversion
      {
      matrix mx(*this);			// Make a copy of this
      mx.m = virtual_to_real_copy(mx.m);// Insure they are independent
      mx.m->resize(i,j);		// Resize the new matrix
      return mx;
      }
    else				// Other matrix conversion
      {
      matrix mx(i, j, n_matrix_type);	// Construct a new matrix
      int ii =i-1;
      int jj=j-1;
      for(r=rows()-1; r>=0; r--)	// Loop backwards over all
	for(c=cols()-1; c>=0; c--)	// elements of this
	  {
	  mx(ii,jj) = (*this)(r,c);
	  jj--;
	  if(jj<0)
	    {
	    jj = j-1;
	    ii--;
	    }
	  }
      return mx;
      }
    }

matrix matrix::diagonal_form()
  {
  matrix mx(rows()*cols(),rows()*cols(),d_matrix_type);
  int r = rows();
  int c = cols();
  int i = r*c;
  i--;
  for(r=rows()-1; r>=0; r--)
    for(c=cols()-1; c>=0; c--)
      {
      mx(i,i)=(*this)(r,c);
      i--;
      }
  return mx;			// return mx
  }

	// Input		mx1   : Matrix (this)
	//			mx    : A second matrix
	// Output		TF    : True if mx & mx point
	//				to the same _matrix

bool matrix::same_reference_as(const matrix& mx) const
  { return (mx.m==m); }

	// Input		mx    : Matrix (this)
	//			full  : Flag for amount of output
	// Output		void  : Outputs matrix status
	// Note			      : Currently there is no matrix
	//				type for unitary so this is 
	//				tested independently

void matrix::status(int full) const
  {
  if(i_matrix_type == stored_type())		// Output stored matrix type
    std::cout << "\n\tIdentity matrix";
  else if(d_matrix_type == stored_type())
    std::cout << "\n\tDiagonal matrix";
  else if(n_matrix_type == stored_type())
    std::cout << "\n\tNormal matrix";
  else if(h_matrix_type == stored_type())
    std::cout << "\n\tHermitian matrix";
  else
    std::cout << "\n\tUnknown or Null matrix type";
  std::cout << ", " << refs() << " Reference";	// Number of references
  if(refs() > 1) std::cout << "s";

  double ecutoff=1.e-5;
  if(full)
    {
    int mxcnt=0, mxtypes[5]={0,0,0,0,0};
    if(i_matrix_type != stored_type() &&
                          i_matrix_type == test_type(i_matrix_type))
      {
      mxtypes[0] = 1;
      mxcnt++;
      }
    if(d_matrix_type != stored_type() &&
                          d_matrix_type == test_type(d_matrix_type, ecutoff))
      {
      mxtypes[1] = 1;
      mxcnt++;
      }
    if(h_matrix_type != stored_type() &&
                          h_matrix_type == test_type(h_matrix_type, ecutoff))
      {
      mxtypes[2] = 1;
      mxcnt++;
      }
    if(n_matrix_type != stored_type() &&
                          n_matrix_type == test_type(n_matrix_type, ecutoff))
      {
      mxtypes[3] = 1;
      mxcnt++;
      }
    if(is_unitary())
      {
      mxtypes[4] = 1;
      mxcnt++;
      }
    std::cout << "\n\tOther Possible Storage Type"; 		// Output possible types
    if(!mxcnt) std::cout << "s: None";
    else
      {
      if(mxcnt > 1) std::cout << "s: ";
      else          std::cout << ": ";
      mxcnt = 0;
      if(mxtypes[0])
        {
        std::cout << "Identity";
        mxcnt++;
        }
      if(mxtypes[1])
        {
        if(mxcnt) std::cout << ", ";
        std::cout << "Diagonal";
        mxcnt++;
        }
      if(mxtypes[2])
        {
        if(mxcnt) std::cout << ", ";
        std::cout << "Hermitian";
        mxcnt++;
        }
      if(mxtypes[3])
        {
        if(mxcnt) std::cout << ", ";
        std::cout << "Normal";
        mxcnt++;
        }
      if(mxtypes[4])
        {
        if(mxcnt) std::cout << ", ";
        std::cout << "Unitary";
        }
      }
    }
  std::cout << "\n";
  return;
  }


matrix inv(const matrix &mx) { return matrix(mx.m->inv()); }

	// Input                mx   : Input matrix
	// Output               mx2  : Inverse of input matrix mx
	//			             mx * mx2 = I



/* LU ***********************************************************-*-c++-*-
**								 	**
**  This routine used here is taken from "Numerical Recipies", W.H.	**
**  Press, B.P.Flannery, S.A.Teukolsky, & W.T.Vetterling, Cambridge 	**
**  University Press, 1989.  See pages 31 through 36. The code has been	**
**  extensively adapted for C++ and GAMMA arrays.			**
**								 	**
**  The algorithm takes an input matrix A and an integer array indx &	**
**  returns the LU decomposition of A' where A' is A with any needed	**
**  row permutations, the latter are stored in the integer array indx	**
**  (indx[0] < 0 flags no row changes).					**
**									**  
**  The LU decomposition formulates the equation (depicted for 4x4 A)	**
**									**  
**		     	        A = LU					**  
**									**
**     [A11 A12 A13 A14]   [ 1   0   0   0 ] [U11 U12 U12 U14]		**
**     |A21 A22 A23 A24]   [L21  1   0   0 ] [ 0  U22 U23 U24|		**
**     |A31 A32 A33 A34] = [L31 L32  1   0 ] [ 0   0  U33 U34|		**
**     [A41 A42 A43 A44]   [L41 L42 L43  1 ] [ 0   0   0  U44]		**
**									**
**  Neglecting implicit <i|L|i>, L & U are stored in in a single array.	**
**									**
**                         [U11 U12 U12 U14]				**	
**  			   [L21 U22 U23 U24|				**
**  	   		   [L31 L32 U33 U34|				**
**  		   	   [L41 L42 L43 U44]				**
**									**
**  The algorithm uses Crouts method in determining the elements of L	**
**  U which is stable under pivoting (row interchanges).		**
**									**
*************************************************************************/

matrix LU(matrix &mx, int* indx)
    { return matrix(mx.m->LU(indx)); }

	// Input                mx   : Input matrix
	//			indx : Row permutation array
	// Output               mx2  : LU decomposition of a row permutation
	//			       of the input matrix mx, mx', the row
	//			       exchanges recorded in array indx
	//			       mx2(low.tri.) * mx2(upp.tri.) = mx'
	//			              L      *       U       = mx'


/* LUinv ********************************************************-*-c++-*-
**								 	**
**  This routine is taken from "Numerical Recipies", W.H.Press, B.P.	**
**  Flannery, S.A.Teukolsky, & W.T.Vetterling, Cambridge University	**
**  Press, 1989.  See pages 31 through 38.  The code has been largely	**
**  modified for C++ and GAMMA.						**
**								 	**
**  The algorithm takes an input matrix LU and in integer array indx,	**
**  where LU is the LU decomposition of A' and A' is a row permutation	**
**  of a matrix A.  The row permuations which relate A' and A are kept	**
**  in indx. The function the proceeds to solve				**
**									**  
**			       A * X = B				**
**									**  
**  for X in a two step process. Using P to indicate row permutations	** 
**  (contained in indx)  The equation is re-formulated as		**
**									**  
**               A * X = (PLU) * X = PL*(U*X) = PL*Y = B 		**
**									**  
**  with U*X = Y, equation     PL * Y = B      is readily solved as	**
**  L is lower triangular		       and P simple row inter-	**
**  changes.  Subsequently, the desired matrix X is obtained from 	**
**									**
**			       U * X = Y				**
**									**
**  which is also easily solved because U is upper-triangular.  Note	**
**  that all row permutations are taken care of in forming Y because	**
**  if A' = LU, A = PA', then A = PLU.					**
**									**
*************************************************************************/

matrix LUinv(matrix &B, int *indx, matrix& LU)

	// Input            B   : Input matrix of results
	//		    LU  : An LU decomposition of matrix A
	//			  (or row permutation of matrix A)
	//		    indx: Row permutation array
	// Output           X   : Solution to the equation: A*X = B
	//			  where LU is the LU decomposition of A
	//			  or row permutation of A
	// Note		        : Matrix LU must be square
	// Note		        : If indx = NULL, and B=I, X will be inv(A)

{ return matrix(B.m->LUinv(indx,LU.m)); }

 
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

std::vector<int> matrix::BlockDiag(matrix& BF, std::vector<int> &U) const
  {
  if(rows() != cols()) 			// Cant block non-square array 
    { 					// so this will be a fatal error
    Mxerror(5, "BlockDiag", 1);	        //   Bad use of BlockDiag function
    Mxfatality(28);			//   Diagonalization on non-square mx
    }
  virtual_delete(BF.m);			// Delete BF, it will be overwritten
  std::vector<int> blks;		// Array of blocks
  blks = m->BlockDiag(BF.m, U);	// Block diagonalize mx: make BF & U
  BF.m = virtual_copy(BF.m);		// Now copy the new matrices as they
  return blks;				// Return array of blocks 
  }

void matrix::HermTriDiag(matrix& HTD, matrix& U) const
  {
  if(rows() != cols()) 			// Cant diagonalize non-square array 
    { 					// so this will be a fatal error
    Mxerror(5, "HermTriDiag", 1);	//   Bad use of HermTriDiag function
    Mxfatality(28);			//   Diagonalization on non-square mx
    }
  virtual_delete(HTD.m);		// Delete HTD, it will be overwritten
  virtual_delete(U.m);			// Delete U, it will be overwritten
  m->HermTriDiag(HTD.m, U.m);		// Tridiagonalize mx: make HTD & U
  HTD.m = virtual_copy(HTD.m);		// Now copy the new matrices as they
  U.m = virtual_copy(U.m);		// have bypassed reference counting
  }

void matrix::SymTriDiag(matrix& STD, matrix& U) const
  {
  if(rows() != cols()) 			// Cant diagonalize non-square array 
    { 					// so this will be a fatal error
    Mxerror(5, "SymTriDiag", 1);	//   Bad use of HermTriDiag function
    Mxfatality(28);			//   Diagonalization on non-square mx
    }
  virtual_delete(STD.m);		// Delete STD, it will be overwritten
  virtual_delete(U.m);			// Delete U, it will be overwritten
  m->SymTriDiag(STD.m, U.m);		// Tridiagonalize mx: make STD, U
  STD.m = virtual_copy(STD.m);		// Now copy the new matrices as they
  U.m = virtual_copy(U.m);		// have bypassed reference counting
  }

void matrix::SymDiag(matrix& D, matrix& U) const
  {
  if(rows() != cols()) 			// Cant diagonalize non-square array 
    { 					// so this will be a fatal error
    Mxerror(5, "SymDiag", 1);		//   Bad use of SymDiag function
    Mxfatality(28);			//   Diagonalization on non-square mx
    }
  virtual_delete(D.m);			// Delete D, it will be overwritten
  virtual_delete(U.m);			// Delete U, it will be overwritten
  m->SymDiag(D.m, U.m);			// Diagonalize mx: make D, U
  D.m = virtual_copy(D.m);		// Now copy the new matrices as they
  U.m = virtual_copy(U.m);		// have bypassed reference counting
  }

void matrix::Diagonalize(matrix& D, matrix& U) const
  {
  if(rows() != cols()) 			// Cant diagonalize non-square array 
    { 					// so this will be a fatal error
    Mxerror(5, "Diagonalize", 1);	//   Bad use of Diagonalize function
    Mxfatality(28);			//   Diagonalization on non-square mx
    }
  virtual_delete(D.m);			// Delete D, it will be overwritten
  virtual_delete(U.m);			// Delete U, it will be overwritten
  m->diag(D.m, U.m);			// Diagonalize mx, make D & U
  D.m = virtual_copy(D.m);		// Now copy the new matrices as they
  U.m = virtual_copy(U.m);		// have bypassed reference counting
  }

void diag(const matrix& mx, matrix& D, matrix& S)
  {
  if(mx.rows() != mx.cols()) 		// Cant diagonalize non-square array 
    { 					// so this will be a fatal error
    mx.Mxerror(5, "diag", 1);		//   Bad use of diag function
    mx.Mxfatality(28);			//   Diagonalization on non-square mx
    }
  if(!matrix::BlkDiag)			// If not trying to first block array
    {
    virtual_delete(D.m);		//   Delete D, it'll be overwritten
    virtual_delete(S.m);		//   Delete S, it'll be overwritten
    mx.m->diag(D.m, S.m);		//   Diagonalize mx, make D & S
    D.m = virtual_copy(D.m);		//   Now copy new matrices as they've
    S.m = virtual_copy(S.m);		//   bypassed reference counting
    return;
    }
  std::vector<int> blks;		// Array of block dimensions
  matrix BD;				// Arrays: Block Form, Permutation
  std::vector<int> P;		// Array of block dimensions
  blks = mx.BlockDiag(BD, P);		// Block diagonalize mx: make BD & P
    int nblks = blks.size();		// Get the number of blocks
  if(nblks == 1)			// If there is only one block then
    {					// we are back to where we were
    virtual_delete(D.m);		//   Delete D, it'll be overwritten
    virtual_delete(S.m);		//   Delete S, it'll be overwritten
    mx.m->diag(D.m, S.m);		//   Diagonalize mx, make D & S
    D.m = virtual_copy(D.m);		//   Now copy new matrices as they've
    S.m = virtual_copy(S.m);		//   bypassed reference counting
    return;
    }
  int nr = mx.rows(); 			// There are blocks, work block-wise
  D = matrix(nr,nr,0.0,d_matrix_type);	// Initialize D for eigenvalues
  matrix S1 = matrix(nr,nr,1.0,d_matrix_type);	// Initialize S for eigenvectors
  int I=0;				// The block row,column index
  matrix subBLK;			// Array for working with sub-blocks
  matrix subD;
  matrix subS;
  int bd; 
//std::cout << "\n\n\tHere is The Original Array: " << mx;
//std::cout << "\n\n\tHere is The Blocked Array: " << BD;
//std::cout << "\n\n\tBlock Diag. On Mx Of Dimension " << nr;
//std::cout << "\n\tThere are " << nblks << " Blocks\n";
  for(int i=0; i<nblks; i++)		// Loop over the blocks
    {
    bd = blks[i];			//   Get block dimension
//std::cout << "\n\tDiagonalizing Block " << i+1
//     << " of Dimension " << bd; 
//std::cout << "\n\tStarting At Index " << I; 
//std::cout.flush();
    if(bd == 1)				//   If block size is only 1 then
      {					//   BD has a isolated diagonal
//std::cout << "\nHere is Block Diagonalized Elem: " << BD(I,I);
      D.put(BD(I,I),I,I);		//   Then eigenvector is 1 column
      I++;				//   and eigenvalue is in array BD
//std::cout << "\nHere is Resulting EigenValue: " << D(I,I);
//std::cout << "\nHere is Resulting EigenVector: " << complex1;
      }
    else				//   If block is larger than one we
      {					//   will clip it out and diagonalize
      virtual_delete(subD.m);		//   Delete D, it'll be overwritten
      virtual_delete(subS.m);		//   Delete S, it'll be overwritten
      subBLK=BD.get_block(I,I,bd,bd);	//   it explicitly
//std::cout << "\nHere is SubBlock To Diagonalize: " << subBLK;
      subBLK.m->diag(subD.m, subS.m);	//   Diagonalize mx, make D & S
      virtual_copy(subD.m);	//   Now copy new matrices as they've
      virtual_copy(subS.m);	//   bypassed reference counting
//std::cout << "\nHere is Resulting EigenValues: " << subD;
//std::cout << "\nHere is Resulting EigenVectors: " << subS;
      D.put_block(I,I,subD);		//   Put block of eigenvalues in
      S1.put_block(I,I,subS);		//   Put block of eigenvectors in
      I += bd;				//   Up index to next block start
      }
//std::cout << "\nHere is D: " << D;
//std::cout << "\nHere is S: " << S1;
//std::cout.flush();
    }
//getrusage(0, & me);
//std::cout << "diag routine: point 03: " << me.ru_utime.tv_sec+me.ru_utime.tv_usec/1.0e6 << " seconds\n";
//std::cout.flush();
//std::cout << D;
//std::cout << S1;
  S = matrix(nr,nr,0.0,d_matrix_type);	// Initialize S for eigenvectors
  for(int i=0;i<nr;++i)
  { for(int j=0;j<nr;++j)
    { S.put(S1.get(i,j),P[i],j);
    }
  }
//std::cout << S;
//getrusage(0, & me);
//std::cout << "diag routine: point 04: " << me.ru_utime.tv_sec+me.ru_utime.tv_usec/1.0e6 << " seconds\n";
//std::cout.flush();
  }

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
   TestUnitary       mx, mx    Constructs the product U*adjoint(U) & tests
   TestUTransform    mx,mx     Test unitary tranform U*T*adj(U)              */

// sosi  - note that this assumes that EV is unitary.....?!
void matrix::TestEigenSystem(int pf) const
  {
  matrix D, EV;				// Diagonalized mx & mx Eigenvalues
  diag(*this, D, EV);			// Get the diagonalized array
  matrix I   = EV*inv(EV);		// Calculate I = mx*inverse(mx)
  matrix mxp = EV*D*inv(EV); 	        // Calculate mxp = EV*D*inverse(EV)
//  matrix I = EV*EV.adjoint();		// Calculate I = mx*adj(mx)
//  matrix mxp = EV*D*EV.adjoint();         // Calculate mxp = EV*D*adjoint(EV)
  if(pf > 2)
    {
    std::cout << "\n\n" << std::string(30, ' ') << "Array Being Tested" << *(this);
    std::cout << "\n\n" << std::string(30, ' ') << "Array EigenVectors" << EV;
    std::cout << "\n\n" << std::string(30, ' ') << "Array EigenValues"  << D;
    if(pf > 3)
      {
      std::cout << "\n\n" << std::string(23, ' ')
           << "Eig.Vec. * inverse(Eig.Vec.) = I?" << I;
      std::cout << "\n\n" << std::string(18, ' ')
           << "Eig.Vec. * Eig.Val * inverse(Eig.Vec.) = Mx?" << mxp;
      }
    }

//               Fill Array ColNorms With Column Norm Values

  int i,j,dim=cols();
  double ColNorm, cn;			// For column norms
  double* ColNorms;			// Array of column norms
  ColNorms = new double[dim];		// Allocate column norms array
  double NormDevTot = 0;
  for(j=0; j<dim; j++)			// Loop over the columns
    {
    ColNorm = 0;			// Begin with norm of zero
    for(i=0; i<dim; i++)		// Loop over the columns
      {
      cn = norm(EV.get(i,j));		// Get norm of <i|eigenvector|j>
      ColNorm += cn*cn;			// Add to column sum
      }
    ColNorms[j] = sqrt(ColNorm);	// Store the column norm
    NormDevTot += fabs(ColNorms[j]-1.0); 
    }

//                      Calculate mx*inv(mx) = 1

  double* ColZeros;			// For summed column off-diags
  ColZeros = new double[dim];		// Allocate column off-diags array
  double ColZero;			// Value for a column 
  double TotalZero=0;
  double TotalIDevR=0, TotalIDevI=0;
  for(i=0; i<dim; i++)			// Loop over the rows
    {
    ColZero=0;				// Begin with no value
    for(j=0; j<dim; j++)		// Loop over the columns
      {
      if(i!=j) 
        ColZero += norm(I.get(i,j));
      else
        {
        TotalIDevR += fabs(I.getRe(i,i) - 1.0);
        TotalIDevI += fabs(I.getIm(i,i));
        }
      }
    ColZeros[i] = ColZero;
    TotalZero += ColZero;
    }

  if(pf > 1)
    {
    std::cout << "\n" << std::string(30, ' ') << "EigenVector Tests";
    std::string marg = std::string(3, ' ');
    std::string tab = std::string(3, ' ');
    std::string lp("(");
    std::string rp(")");
    std::string line1 = marg + "Column";
    std::string line2 = marg + " Index";
    std::string line3 = marg + "------";
    line1 +=       tab  + std::string("      Column Norm      "); 
    line2 +=       tab  + std::string("Value (Deviation x10^9)"); 
    line3 +=       tab  + std::string("-----------------------");
    line1 +=       tab  + std::string("              Mx*adjoint(Mx)=I");
    line2 +=       tab  + std::string("Off-diag Sum   Diagonal   (Deviation x10^9)");
    line3 +=       tab  + std::string("-------------------------------------------");
    std::cout << "\n\n" << line1;
    std::cout << "\n"   << line2;
    std::cout << "\n"   << line3;
    for(i=0; i<dim; i++)
      {
      line1 =  marg + MxModform("%4d", i) + std::string("  "); 
      line1 += tab  + MxModform("%7.3f", ColNorms[i]) + std::string(" ")
                    + lp + MxModform("%12.8f", (ColNorms[i]-1)*1.e9) + rp;
      line1 += tab  + MxModform("%10.3f", ColZeros[i]) + std::string(5,' ')
                    + MxModform("%8.4f", norm(I.get(i,i))) + std::string(5, ' ')
                    + MxModform("%12.8f", fabs(norm(I.get(i,i))-1)*1.e9);
      std::cout << "\n" << line1;
      }
    }
  std::cout << "\n";
  std::cout << "\n\tColumn Norm Total Deviation:           " << NormDevTot;
  std::cout << "\n\tInverse Total Deviation:               " 
       << TotalIDevR + TotalIDevI + TotalZero;

//                      Calculate EV*D*inverse(EV) = mx

  double DEV =0;
  for(i=0; i<dim; i++)				// Loop over the rows
    for(j=0; j<dim; j++)			// Loop over the columns
      DEV += norm(get(i,j)-mxp.get(i,j));

  std::cout << "\n\tMx-Ev*D*inverse(Ev) Total Deviation:   " << DEV;

  delete [] ColNorms;
  delete [] ColZeros;
  }


void matrix::TestTransform(const matrix& T, const matrix& S, int pf) const
  {
  matrix I = S.TestUnitary(std::cout);	// Test/Get I=S*adj(S)
  matrix mxp = TestUTransform(T, S); 	// Test/Get S*T*adjoint(S)
  if(pf > 2)
    {
    std::cout << "\n\n" << std::string(30, ' ') << "Array Being Tested\n" << *(this);
    std::cout << "\n\n" << std::string(30, ' ') << "Array Transform\n" << S;
    std::cout << "\n\n" << std::string(30, ' ') << "Array Transformed\n"  << T;
    if(pf > 3)
      {
      std::cout << "\n\n" << std::string(23, ' ') << "S * adj(S) = I?\n" << I;
      std::cout << "\n\n" << std::string(18, ' ') << "S * T * inv(S) = Mx?\n" << mxp;
      }
    }
  }


	// Input	mx	: A matrix (this)
	// Output	vect    : A vector of the column norms
	//			  of mx
	// Note			: I use this to look at unitary
 	//			  arrays resulting from Hermitian
	//			  matrix diagonalizations

std::vector<double> matrix::ColumnNorms() const
  {
  int i,j,dim=cols();
  double ColNorm, cn;                   // For column norms
  std::vector<double> ColNorms(dim);		// Allocate column norms array
  for(j=0; j<dim; j++)                  // Loop over the columns
    {
    ColNorm = 0;                        // Begin with norm of zero
    for(i=0; i<dim; i++)                // Loop over the columns
      {  
      cn = norm(get(i,j));		// Get norm of <i|U|j>
      ColNorm += cn*cn;                 // Add to column sum
      }  
    ColNorms[j] = sqrt(ColNorm);        // Store the column norm
    }
  return ColNorms;
  }


std::vector<double> matrix::TestIdentity(complex& TotalDev) const

	// Input	mx	: A matrix (this)
	//		TotalDev: A placeholder for more output
	// Output	vect    : A vector of the column sums over
	//			  off-diagonal elements
	//		TotalDev: Set to contain the total real
	//			  and imaginary deviation from 1
	//			  summed over the diagonal
	// Note			: I use this to look at unitary
 	//			  arrays resulting from Hermitian
	//			  matrix diagonalizations where
	//				I = U*adjoint(U)

  {
  int i,j,dim=cols();
  std::vector<double> ColZeros(dim);		// Allocate column norms array
  double ColZero;                       // Value for a column
  double TotalIDevR=0, TotalIDevI=0;
  for(i=0; i<dim; i++)                  // Loop over the rows
    {
    ColZero=0;                          // Begin with no value
    for(j=0; j<dim; j++)                // Loop over the columns
      {
      if(i!=j)
        ColZero += norm(get(i,j));
      else
        {
        TotalIDevR += fabs(getRe(i,i) - 1.0);
        TotalIDevI += fabs(getIm(i,i));
        }
      }
    ColZeros[i] = ColZero;
    }
  TotalDev = complex(TotalIDevR, TotalIDevI);
  return ColZeros;
  }


matrix matrix::TestUnitary(std::ostream& ostr) const

	// Input        mx      : A matrix (this)
        //              ostr    : An output stream
        // Output       Imx     : Product of mx*adjoint(mx)
        //                        Output stream ostr is modified
        //                        with info as to how close Imx is
        //                        to an indentity matrix

  {

//        Calculate Column Norms, Should Be 1 For Unitary Arrays

  int i, j, dim=cols();				// Get array dimension
  std::vector<double> ColNorms = ColumnNorms();	// Get column norms
  double NormDevTot = 0;			// For total norm deviation
  for(j=0; j<dim; j++)				// Loop over the columns
    NormDevTot += fabs(ColNorms[j]-1.0); 	// Get total norm deviation

//               Calculate mx*adjoint(mx) = mx*inv(mx) = 1
//           This Is Valid Only For Unitary Arrays Of Course

  matrix I = (*this)*adjoint();		// Calculate I = U*adj(U)
  complex TotalIDev;				// For deviation from I
  std::vector<double> ColZeros 			// Get total deviation from I
                   = I.TestIdentity(TotalIDev);
  double TotalZero=0;				// Total stray from zero
  for(i=0; i<dim; i++)				// Loop over the rows
    TotalZero += ColZeros[i];			// Get total stray from 0
  double TotalIDevR = zRe(TotalIDev);
  double TotalIDevI = zIm(TotalIDev);

  ostr << "\n" << std::string(30, ' ') << "Unitary Array Tests";
  std::string marg = std::string(3, ' ');
  std::string tab = std::string(3, ' ');
  std::string lp("(");
  std::string rp(")");
  std::string line1 = marg + "Column";
  std::string line2 = marg + " Index";
  std::string line3 = marg + "------";
  line1 +=       tab  + std::string("      Column Norm      "); 
  line2 +=       tab  + std::string("Value (Deviation x10^9)"); 
  line3 +=       tab  + std::string("-----------------------");
  line1 +=       tab  + std::string("              Mx*adjoint(Mx)=I");
  line2 +=       tab  + std::string("Off-diag Sum   Diagonal   (Deviation x10^9)");
  line3 +=       tab  + std::string("-------------------------------------------");
  ostr << "\n\n" << line1;
  ostr << "\n"   << line2;
  ostr << "\n"   << line3;
  for(i=0; i<dim; i++)
    {
    line1 =  marg + MxModform("%4d", i) + std::string("  "); 
    line1 += tab  + MxModform("%7.3f", ColNorms[i]) + std::string(" ")
                  + lp + MxModform("%12.8f", (ColNorms[i]-1)*1.e9) + rp;
    line1 += tab  + MxModform("%10.3f", ColZeros[i]) + std::string(5,' ')
                  + MxModform("%8.4f", norm(I.get(i,i))) + std::string(5, ' ')
                  + MxModform("%12.8f", fabs(norm(I.get(i,i))-1)*1.e9);
    ostr << "\n" << line1;
    }
  ostr << "\n\n\tColumn Norm Total Deviation:     " << NormDevTot;
  ostr << "\n\tAdjoint=Inverse Total Deviation: " 
       << TotalIDevR + TotalIDevI + TotalZero;
  return I;
  }


matrix matrix::TestUTransform(const matrix& T, const matrix& U) const

        // Input        mx      : A matrix (this)
	//		T	: Transformed mx
	//		U       : Transformation m atrix
        //              ostr    : An output stream
	//		mx'	: U*T*adj(U)

//           Calculate U*D*adjoint(U) = U*D*inv(U) = mx
//          This Is Valid Only For Hermitian Arrays Of Course
  {
  matrix mxp = U*T*U.adjoint();         // Calculate mxp = U*T*adjoint(U)
  int i,j,dim=cols();
  double DEV =0;
  for(i=0; i<dim; i++)                  // Loop over the rows
    for(j=0; j<dim; j++)                // Loop over the columns
      DEV += norm(get(i,j)-mxp.get(i,j));
  std::cout << "\n\tMx-U*D*adj(U) Total Deviation:   " << DEV;
  return mxp;
  }

// ____________________________________________________________________________
// Q                MATRIX CONTAINER SUPPORT FUNCTIONS
// ____________________________________________________________________________

/* Aside for providing basic tests as to whether to matrics are equvalent or
   not, these operators are necessary in any STL container classes are to
   be used based on matrices (e.g. list<matrix> or std::vector<matrix>)           */

bool matrix::operator== (const matrix& mx2) const
 {
 if(same_reference_as(mx2)) return true;		// Try for easy + test!
 if(rows()!=mx2.rows())     return false;		// Try for easy - test!
 if(cols()!=mx2.cols())	    return false;		// Try for easy - test!

//       Not Clearly Equal Or Unequal, We're Forced To Element Compare

 matrix_type mt1 =     (*m).stored_type();		// Get 1st array type
 matrix_type mt2 = (*(mx2.m)).stored_type();		// Get 2nd array type
 if(mt1 == i_matrix_type) return m->is_equal(mx2.m);	// Easy if array imx 
 if(mt2 == i_matrix_type) return mx2.m->is_equal(m);    // Easy if array imx
 if(mt1 == d_matrix_type) return m->is_equal(mx2.m);	// Easy if array dmx 
 if(mt2 == d_matrix_type) return mx2.m->is_equal(m);    // Easy if array dmx
 if(mt1 == h_matrix_type) return m->is_equal(mx2.m);	// Harder if array hmx 
 if(mt2 == h_matrix_type) return mx2.m->is_equal(m);    // Harder if array hmx
 if(mt1 == n_matrix_type) return m->is_equal(mx2.m);	// Hardest if array nmx 
 if(mt2 == n_matrix_type) return mx2.m->is_equal(m);    // Hardest if array nmx
 bool flag=true;					// If we just don't
 int r, c;						// know either matrix
 for(r=rows()-1; (r>=0)&&flag; r--)			// type we do a brute
   for(c=cols()-1; (c>=0)&&flag; c--)			// force comparison
     flag=(get(r,c)==mx2.get(r,c));
 return flag;
 }

bool matrix::operator!=(const matrix& mx) const {return !(*this==mx); }	
bool matrix::operator<(const matrix& mx) const  {return (rows()<mx.rows());}
bool matrix::operator>(const matrix& mx) const  {return (rows()>mx.rows());}

#endif						// matrix.cc
