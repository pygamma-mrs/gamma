/* col_vector.h *************************************************-*-c++-*-
**                                                                      **
**                                  G A M M A                           **
**                                                                      **
**	Column Vectors		     	                  Interface	**
**                                                                      **
**      Copyright (c) 1990, 2000                                        **
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
**  Description                                                         **
**                                                                      **
**  Class col_vector defines column vectors for C++ with the usual      **
**  algebraic operations, higher functions and I/O routines.            **
**                                                                      **
*************************************************************************/

#ifndef   Gcol_vector_h_		// Is file already included?
#  define Gcol_vector_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Include system specifics
#include <iostream>			// Include libstdc++ io streams
#include <vector>			// Include libstdc++ STL vectors
#include <Matrix/complex.h>		// Include GAMMA complex numbers
#include <Matrix/matrix.h>	        // Include GAMMA matrices
#include <Matrix/row_vector.h>		// Include GAMMA row vectors

class row_vector;			 // Know that row_vector is a class

class col_vector : public matrix
   {
   friend class row_vector;

private:

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS COLUMN VECTOR ERROR HANDLING
// ____________________________________________________________________________


/*      Input                   cvec    : A column vector (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
                                pname   : string in message
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

         void CVerror(int eidx, int noret=0) const;
         void CVerror(int eidx, const std::string& pname, int noret=0) const;
volatile void CVfatality(int eidx) const;
volatile void CVfatality(int eidx, const std::string& pname) const;

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:


// ____________________________________________________________________________
// A                  COLUMN VECTOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

///F_list col_vector            - Constructor
///F_list ~                     - Destructor

/*    Arguments                         Constructed Column Vector
    -------------       -------------------------------------------------------
          -             An empty column vector (size = 0)
          nc            A col vector with nc elements (columns) uninitialized
        nc, z           A col vector with nc elements (columns) all set to z
          rv            A col vector equivalent to the input vector rv
          mx            A col vector equivalent to input matrix mx
                        (in this case mx MUST be of dimension nx1)           */

MSVCDLC col_vector ( );
MSVCDLC col_vector (int i);
MSVCDLC col_vector (int i, double d);
MSVCDLC col_vector (int i, const complex& z);
MSVCDLC col_vector (const col_vector& cvec);
MSVCDLC col_vector (const matrix& mx);
MSVCDLC ~col_vector ();

// ____________________________________________________________________________
// B               COLUMN VECTOR ACCESS AND ASSIGNMENT
// ____________________________________________________________________________

/*  Function/Operator                         Result
    -----------------   ------------------------------------------------------
          =             Current column vector set equal to input cvec
        (int)           Reference to element <i|cvec>
       get(int)         Copy of element <i|cvec>
      getRe(int)        Copy of real part of element <i|cvec>
      getIm(int)        Copy of imaginary part of element <i|cvec>
       put(int)         Assigns element <i|cvec>                            */

MSVCDLL col_vector& operator = (const col_vector& cvec);
MSVCDLL col_vector& operator = (const row_vector& rvec);
MSVCDLL col_vector& operator = (const matrix& mx);
MSVCDLL complex&    operator() (int i);
MSVCDLL complex     get(int i) const;
MSVCDLL double      getRe(int i) const;
MSVCDLL double      getIm(int i) const;
MSVCDLL void        put(const complex& z, int i);

// ____________________________________________________________________________
// E             CLASS COLUMN VECTOR VARIOUS CHECKS & PARAMETERS
// ____________________________________________________________________________

/*  Function  Output                       Description
    --------  ------  ---------------------------------------------------------
    elements   int    Returns # of vector rows
      size     int    Returns # of vector rows
      rows     int    Returns # of vector rows       - handled by matrix
      cols     int    Returns # of vector columns    - handled by matrix
      refs     int    Returns # of vector references - handled by matrix
      pts      int    Returns # of vector elements   - handled by matrix     */

MSVCDLL int elements( ) const;
MSVCDLL int size( )     const;
/*
int matrix::rows();		     // Number of matrix rows       INHERITED
int matrix::cols();		     // Number of matrix columns    INHERITED
int matrix::refs();		     // Number of matrix references INHERITED
int matrix::pts();		     // Number of matrix points     INHERITED
*/

// ____________________________________________________________________________
// F                COLUMN VECTOR BINARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________

/* The operators which must be defined here, as opposed to those inherited from
   class matrix, are they which return col_vector (or at least, return a type
   other than what class matrix does.)D efinitions herein also allow for more
   explicit error messages.

   Additions/subtractions between column and row vectors are not defined here,
   consequently the error messages from doing so stems from class matrix.
   The only possible additions is when both are length 1 anyway.

   There are type casting problems when class matrix takes over some operators
   (e.g. mx*cvec, mx+cvec, mx-cvec). In that situation the result is proper but
   the return is matrix rather than column vector as one might like.

   Operator Arguments      Result        Operator Arguments       Result
   -------- --------- -----------------  -------- --------- -------------------
      +      cv1,cv2      cv1+cv2           *       mx,cv       mx*cv == cv1
      +      cv,mx        cv+mx             *       cv,z           z*cv
      -      cv1,cv2      cv1-cv2           *       cv,d           d*cv
      -      cv,mx        cv-mx             *       z,cv           z*cv
      *      cv1,cv2      cv1*cv2 (err)     *       d,cv           d*cv
      *      cv,rv        cv*rv == mx       /       cv,z         (1/z)*cv
      *      cv,mx        cv*mx == mx1      /       cv,d         (1/d)*cv    */

MSVCDLL friend col_vector operator * (const complex& z, const col_vector& cvec);
MSVCDLL friend col_vector operator * (double         d, const col_vector& cvec);
MSVCDLL friend col_vector operator * (const matrix& mx, const col_vector& cvec);

MSVCDLL col_vector operator + (const col_vector& cvec) const;
MSVCDLL col_vector operator + (const matrix&     mx)   const;
MSVCDLL col_vector operator - (const col_vector& cvec) const;
MSVCDLL col_vector operator - (const matrix&     mx)   const;
MSVCDLL complex    operator * (const col_vector& cvec) const;
MSVCDLL matrix     operator * (const row_vector& rvec) const;
MSVCDLL matrix     operator * (const matrix&     mx)   const;
MSVCDLL col_vector operator * (const complex&    z)    const;
MSVCDLL col_vector operator * (double            d)    const;
MSVCDLL col_vector operator / (const complex&    z)    const;
MSVCDLL col_vector operator / (double            d)    const;

// ____________________________________________________________________________
// G                COLUMN VECTOR UNARY ARITHMETIC FUNCTIONS
// ____________________________________________________________________________

// Note: These operators, as they are unable to change the type of the
//       result, can only be those which return col_vector.

MSVCDLL col_vector& operator += (const col_vector& cvec1);
MSVCDLL void operator += (const matrix&        mx);
MSVCDLL col_vector& operator -= (const col_vector& cvec1);
MSVCDLL void operator -= (const matrix&        mx);
MSVCDLL col_vector& operator *= (const complex&        z);
MSVCDLL col_vector& operator *= (      double          d);
MSVCDLL col_vector & operator /= (const complex&        z);
MSVCDLL col_vector & operator /= (const double          d);

// ____________________________________________________________________________
// J                   COLUMN VECTOR SIMPLE UNARY FUNCTIONS
// ____________________________________________________________________________

// Note: The negation operator must be defined herein.  The function conjugate
//       is inherited, but adjoint and transpose return col_vectors so they
//       are also here.  Transpose is included just for specific error message.

/*  Function  Output                       Description
    --------  ------  ---------------------------------------------------------
      -       cvec    Returns the vector mutiplied by -1 (negation) -<cvec|
    adjoint   rvec    Returns row vector where <i|cvec> = <i|cvec>*
   transpose  rvec    Returns row vector where <i|cvec> = <i|cvec>
     trace    ----    This is disallowed for row vectors (but MUST exist)    */

// col_vector operator- ();                   // MATRIX INHERITED

MSVCDLL friend row_vector adjoint(const   col_vector& cvec);
MSVCDLL friend row_vector transpose(const col_vector& cvec);
MSVCDLL friend complex    trace(const     col_vector& cvec);

MSVCDLL row_vector adjoint()      const;
MSVCDLL row_vector transpose()    const;
MSVCDLL complex    trace()        const;
MSVCDLL col_vector differential() const;

// ____________________________________________________________________________
// I                    COLUMN VECTOR SIMPLE BINARY FUNCTIONS
// ____________________________________________________________________________

/*  Function  Output                       Description
    --------  ------  ---------------------------------------------------------
      norm    double  Returns the vector norm, sq.rt. of sum of squares
      sum     complex Returns the sum ov all vector elements
     maxRe    double  Returns largest real value in the vector
     maxIm    double  Returns largest imaginary value in the vector
     maxZ     complex Returns largest complex (biggest norm) value in vector
     minRe    double  Returns smallest real value in the vector
     minIm    double  Returns smallest imaginary value in the vector
     minZ     complex Returns smallest complex (biggest norm) value in vector
      max      int    Returns index of largest (  real, imag, norm } element
      min      int    Returns index of smallest { real, imag, norm } element
     flip     void    Inverses the element order in the vector
     zero     void    Sets all elements to zero
     sort     vector  Integer vector of sorted indices

  The vector norm is given by

                     [ size-1          ]
                     |  ---            |
                     |  \            2 |            2                   *
          Norm = sqrt|  /  |<i|cvec>|  |, |<i|cvec>| = <i|cvec>*<i|cvec>
                     |  ---            |
                     [  i=0            ]

  & the following definitions apply to the max, min, and sort functions:

                          type         max/min
                            0       norm (default)
                           <0         imaginary
                           >0           real                                 */

MSVCDLL double           norm()          const;
MSVCDLL complex          sum()           const;
MSVCDLL double           maxRe()         const;
MSVCDLL double           maxIm()         const;
MSVCDLL complex          maxZ()          const;
MSVCDLL double           minRe()         const;
MSVCDLL double           minIm()         const;
MSVCDLL complex          minZ()          const;
MSVCDLL int              max(int type=0) const;
MSVCDLL int              min(int type=0) const;
MSVCDLL void             flip();
MSVCDLL complex          sum(int st, int ne) const;
MSVCDLL void             zero();
MSVCDLL std::vector<int> sort(int type=0) const;

// ____________________________________________________________________________
// J                  COLUMN VECTOR COMPLEX UNARY FUNCTIONS
// ____________________________________________________________________________

// ******************************      FFT      *******************************

MSVCDLL friend col_vector FFT(const  col_vector& cvec);
MSVCDLL friend col_vector IFFT(const col_vector& cvec);

        // Input        cvec : A column vector
        // Output      cvec1 : A column vector containing the discrete Fourier
        //                     transform of the values in cvec
        // Note              : The vector cvec should be of dimension
        //                     1 x (2^n)

// ____________________________________________________________________________
// K                     COLUMN VECTOR CONTAINER SUPPPORT
// ____________________________________________________________________________

        // Input      cvec1  : A column vector
        //            cvec2  : A column vector
        // Output ==    T/F  : Returns TRUE if cvec1 is equal to cvec2.
	// Output !=    T/F  : Returns TRUE if cvec1 is NOT equal to cvec2.
        // Output >     T/F  : Returns TRUE if cvec1 is equal to cvec2.
	// Output <     T/F  : Returns TRUE if cvec1 is NOT equal to cvec2.

//bool operator==(const col_vector& cvec) const	// INHERITED
//bool operator!=(const col_vector& cvec) const	// INHERITED
//bool operator> (const col_vector& cvec) const	// INHERITED
//bool operator< (const col_vector& cvec) const	// INHERITED

// ____________________________________________________________________________
// J                COLUMN VECTOR COMPLEX BINARY FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                              Straight Product
// ----------------------------------------------------------------------------

        // Input              cvec1  : A col vector
        //               cvec2,rvec  : Second column vector, or row vector
        // Output             rvec3  : A "product" of two vectors
        // Note                      : This is NOT a scalar product

        //                             cvec3(i) = cvec1(i) * cvec1(i)
        //                                      = cvec1(i) * cvec2(i)
        //                                      = cvec(i) * rvec(i)

MSVCDLL friend col_vector product(const col_vector& cvec1, const col_vector& cvec2);
MSVCDLL friend col_vector product(const col_vector& cvec,        row_vector& rvec);

MSVCDLL col_vector product()                        const;
MSVCDLL col_vector product(const row_vector& rvec)  const;
MSVCDLL col_vector product(const col_vector& cvec2) const;

// ----------------------------------------------------------------------------
//                              Scalar Product
// ----------------------------------------------------------------------------

        // Input              cvec1  : A column vector
        //               cvec2,rvec  : Second column vector, or row vector
        // Output             z      : Scalar product of adjoint(cvec) & cvec

        //                             z = <cvec'|cvec>,   <cvec'| = adjoint(|rvec>)
        //                               = <cvec1'|cvec2>, <cvec2'| = adjoint(|rvec>)
        //                               = <rvec|cvec>

MSVCDLL friend double  scalar_product(const col_vector& cvec);
MSVCDLL friend complex scalar_product(const col_vector& cvec1, const col_vector& cvec2);
MSVCDLL friend complex scalar_product(const col_vector& cvec,  const row_vector& rvec);

MSVCDLL double  scalar_product()                        const;
MSVCDLL complex scalar_product(const col_vector& cvec2) const;
MSVCDLL complex scalar_product(const row_vector& rvec)  const;

// ____________________________________________________________________________
// L                       COLUMN VECTOR I/O FUNCTIONS
// ____________________________________________________________________________

// --------------------- ASCII Input/Output Functions -------------------------

/*              Input           cvec    : A column vector (this)
                                ostr    : Output ASCII file stream
                                full    : Flag for amount of output
                Return          void    : cvec is sent to the output stream  */

MSVCDLL std::string              hdrString() const;
MSVCDLL std::vector<std::string> printStrings(const MxPrint& PFlgs) const;

MSVCDLL        std::ostream& printcols(std::ostream& ost, int nc=4,int rc=1,int ne=0) const;
MSVCDLL        std::ostream& print(std::ostream& ostr, int full=0) const;
MSVCDLL friend std::ostream& operator << (std::ostream& ostr, const col_vector& cvec);
MSVCDLL friend std::istream& operator >> (std::istream& istr, col_vector& cvec);

        // Input                rvec : Row vector (this)
        // Output               void : The function sends questions to
        //                             standard output interactively asking
        //                             the user to supply the information to
        //                             specify the vector.  rvec is modified.
        // Note                      : Function is for INTERACTIVE programs

MSVCDLL void ask();

};

#endif								// col_vector.h
