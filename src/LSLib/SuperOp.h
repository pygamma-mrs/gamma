/* SuperOp.h ****************************************************-*-c++-*-
**							         	**
**                                 G A M M A				**
**						         		**
**	Superoperators                                 Interface	**
**									**
**	Copyright (c) 1993						**
**	Scott Smith							**
**	Eidgenoessische Technische Hochschule				**
**	Labor fuer physikalische Chemie					**
**	8092 Zuerich / Switzerland					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
**  Description								**
**									**
**  Class superoperator defines a superoperator in Liouville space for	**
**  C++.  The class contains common algebraic operations (+, -, *, /),	**
**  more complex mathematical functions	 (exp) and I/O routines.	**
**									**
*************************************************************************/

///Chapter Class Superoperator (super_op)
///Section Overview
///Body    The class
///Section Available Superoperator Functions

#ifndef   Super_op_h_			// Is this file already included?
#  define Super_op_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Include system specifics
#include <Matrix/matrix.h>		// Know about matrices
#include <Matrix/complex.h>		// Know about GAMMA complex #s
#include <HSLib/Basis.h>		// Know about bases
#include <HSLib/GenOp.h>		// Know about general operators
#include <string>			// Know stdlibc++ string class

class super_op;  // Forward declaration so following functions will compile

// These declarations were added here,
// as the "friend" versions were not considered declarations
// but rather as defining the relationship to the super_op class.

MSVCDLL super_op left(const gen_op& Op); 	// LOp*Op1 = Op*Op1
MSVCDLL super_op left(const matrix& mx, const basis& bs);
MSVCDLL super_op left(const matrix& mx); 
MSVCDLL super_op right(const gen_op& Op); 	// LOp*Op1 = Op1*Op
MSVCDLL super_op right(const matrix& mx); 	// LOp*Op1 = Op1*mx
MSVCDLL super_op right(const matrix& mx, const basis& bs);
MSVCDLL super_op Hsuper(const gen_op& Heff);
 
MSVCDLL super_op U_transform(const gen_op& Op);
MSVCDLL super_op U_transform(const matrix& mx);

MSVCDLL super_op commutator(const gen_op& Op); 	// LOp*Op1 = [Op, Op1]
MSVCDLL super_op commutator(const matrix& mx);	// LOp*mx = [mx, mx1] 

MSVCDLL super_op d_commutator(const gen_op& Op, const complex& z = complex1);
MSVCDLL super_op d_commutator(const matrix& mx);
MSVCDLL super_op d_commutator(const gen_op& Op1, const gen_op& Op2);
MSVCDLL super_op d_commutator(const gen_op& Op1, const gen_op& Op2,
                                                             const complex& z);
MSVCDLL super_op d_commutator (const matrix& mx1, const matrix& mx2);



class super_op
 {

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------

  mutable matrix mx;			// LOp matrix (in Liouville space)
  mutable basis Hbs;			// LOp basis (in Hilbert space)
  mutable basis Lbs;			// LOp basis (in Liouville space)
  int HSp;				// Hilbert Space Dimension
  int LSp;				// Liouville Space Dimension


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     SUPEROPERATOR ERROR FUNCTIONS
// ____________________________________________________________________________

        // Input                LOp     : Superoperator (this)
        //                      eidx    : Error index
        //                      noret   : Flag for linefeed (0=linefeed)
        //                      pname   : string in message

         void LOperror(int eidx, int noret=0) const;
         void LOperror(int eidx, const std::string& pname, int noret=0) const;
volatile void LOpfatal(int eidx) const;
volatile void LOpfatal(int eidx, const std::string& pname) const;

// ____________________________________________________________________________
// ii           USEFUL INTERNAL FUNCTIONS FOR SUPER OPERATORS
// ____________________________________________________________________________


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A           CLASS SUPER OPERATOR SET CONSTRUCTORS/DESTRUCTOR
// ____________________________________________________________________________

///Center Basic Functions
// F_list super_op	     - Constructor

MSVCDLC super_op();
MSVCDLC super_op(const matrix& mx);
MSVCDLC super_op(const matrix& mx,  const matrix& bs);
MSVCDLC super_op(const matrix& mx,  const basis& bs);
MSVCDLC super_op(const std::vector<matrix>& mxc, const std::vector<matrix>& bsc);
MSVCDLC super_op(matrix* mxc, int nc, matrix* bsc=NULL);

        // Input                nc  : Number of components in multi-sys
        //                      mxc : SOp submatrices associated with comps
        //                      bsc : HS basis matrices associated with comps
        // Output               SOp : Super_op (this) containing
        //                            mxc as diagonal blocks of mx, and bsc
        //                            as diagonal blocks of Hbs
	// Note			    : Added in support of multi_sys!


MSVCDLC super_op(const super_op& LOp1);
MSVCDLC super_op(const gen_op& Op1, const gen_op& Op2);
MSVCDLC ~super_op( );

// ____________________________________________________________________________
// B      SUPER OPERATOR FUNCTIONS: SUPER OPERATOR WITH SUPER OPERATOR
// ____________________________________________________________________________

///Center Superoperator with Superoperator Functions
// F_list +		     - Addition of superoperators
// F_list +=		     - Unary addition of superoperators
// F_list -		     - Subtraction of superoperators
// F_list -		     - Negation of superoperators
// F_list -=		     - Unary subtraction of superoperators
// F_list *		     - Multiplication of superoperators
// F_list *=		     - Unary multiplication of superoperators
// F_list =		     - Assignment of superoperators

        // Input                LOp  : A super operator (this)
        //                      LOp1 : A super operator
	// Return (+)		LOpx : Superoperator which is the addition
        //	                       of the two input superoperators,
	//			               LOpx =  LOp + LOp
	//			       Result in basis of LOp
	// Return (+=)		LOp  : Superoperator (this) having the input
        //	                       superoperator LOp1 added to it
	//			               LOp = LOp + LOp1
	//			       Result in basis of LOp
	// Return (-)		LOpx : Superoperator which is the subtraction
        //	                       of the two input superoperators,
	//			               LOpx = LOp - LOp1
	//			       Result in basis of LOp
	// Return (-)		LOpx : Superoperator which is the negative
        //	                       of the input superoperator
	//                                       LOp = - LOp1
	// Return (-=)		LOp  : Superoperator (this) having the input
        //	                       superoperator LOp1 subtracted from it
	//			                 LOp = LOp - LOp1
	//			       Result in basis of LOp

MSVCDLL super_op operator +  (const super_op& LOp1) const;
MSVCDLL super_op& operator += (const super_op& LOp1);
MSVCDLL super_op operator -  (const super_op& LOp1) const;
MSVCDLL super_op operator -  () const;
MSVCDLL super_op& operator -= (const super_op& LOp1);

        // Input                LOp  : A super operator (this)
        //                      LOp1 : A super operator
	// Return (*)		LOpx : Superoperator which is the multiplication
        //	                       of the two input superoperators,
	//			             LOpx = LOp1 * LOp2
	//			       Result in basis of LOp1
	// Return (*=)		LOp  : Superoperator (this) multiplied into
	//			       input superoperator LOp1
	//			              LOp = LOp * LOp1
	//			       Result in basis of LOp
        // Return (&=)		LOp  : Super Operator (this) which
        //                             has been multiplied by LOp1
        //                                    LOp = LOp1 * LOp
	//			       Result in basis of LOp1
	// Return (=)		LOp1 : Superoperator which is a copy of
	//			       the input superoperator.

MSVCDLL super_op operator *  (const super_op& LOp1) const;
MSVCDLL super_op& operator *= (const super_op& LOp1);
MSVCDLL super_op& operator &= (const super_op& LOp1);
MSVCDLL void     operator  = (const super_op& LOp1);

// ____________________________________________________________________________
// C        SUPER OPERATOR FUNCTIONS: SUPER OPERATOR WITH OPERATOR 
// ____________________________________________________________________________

///Center Superoperator with Operator Functions
// F_list *		     - Superoperator times operator

MSVCDLL friend gen_op operator * (const super_op& LOp, const gen_op& Op1);
MSVCDLL friend gen_op operator * (const gen_op& Op1, const super_op& LOp);

	// Input		LOp  : Superoperator.
	// 			Op1  : General operator.
	// Return		LOp  : General operator which is the
	//			       multiplication of the input
	//			       superoperator and the operator,
	//			       Op =  LOp * Op1.
	// Note		    	     : Result in basis of LOp1
	// Note			     : Order matters, Op*LOp1 is undefined!


// ____________________________________________________________________________
// D         SUPER OPERATOR FUNCTIONS, SUPER OPERATOR WITH SCALAR
// ____________________________________________________________________________

///Center Superoperator with Scalar Functions
// F_list *		     - Superoperator multiplication by scalar 
// F_list *=		     - Unary LOp multiplication by scalar 
// F_list /		     - Superoperator division by scalar 
// F_list /		     - Unary superoperator division by scalar 

/* These functions allow for simple mathematical operations between a super-
   operators in Liouville space and an constant.  Functions exist that provide
   users with the ability to multiply and divide superoperators with constants.

   Operator Arguments      Result        Operator Arguments       Result
   -------- --------- -----------------  -------- --------- -------------------
      *      LOp, z        z*LOp             /      LOp,z       (1/z)*LOp
      *      z, LOp        z*LOp            /=      this,z      (1/z)* this
      *=     this, z       z*this                                            */

MSVCDLL friend super_op  operator *  (const super_op& LOp1, const complex& z);
MSVCDLL friend super_op  operator *  (const complex& z,     const super_op& LOp1);
MSVCDLL friend super_op  operator *  (const super_op& LOp1, double d);
MSVCDLL friend super_op  operator *  (double d,             const super_op& LOp1);
MSVCDLL super_op&   operator *= (const complex& z);
MSVCDLL super_op&   operator *= (double d);
MSVCDLL friend super_op  operator /  (const super_op& LOp1, const complex& z);
MSVCDLL friend super_op  operator /  (const super_op& LOp1, double d);
MSVCDLL super_op&   operator /= (const complex& z);
MSVCDLL super_op&   operator /= (double d);

// ____________________________________________________________________________
// E                    SUPER OPERATOR COMPLEX FUNCTIONS
// ____________________________________________________________________________

///Center Superoperator Complex Functions
// F_list left		      - Left translation superoperator
// F_list right		      - Right translation superoperator
// F_list commutator	      - Commutation superoperator
// F_list d_commutator	      - Double Commutator superoperator
// F_list U_transform	      - Unitary transformation superoperator
// F_list project 	      - Projection superoperator
// F_list exp 	     	      - Superoperator exponential

// ----------------------------------------------------------------------------
//                  Left And Right Transformation Superoperators
// ----------------------------------------------------------------------------

/* Left:   LOp*Op1 = Op*Op1           where             LOp = Op (X) E

                                                                      t
   Right:  LOp*Op1 = Op1*Op           where             LOp = E (X) Op       */

MSVCDLL friend super_op left(const gen_op& Op); 	// LOp*Op1 = Op*Op1
MSVCDLL friend super_op left(const matrix& mx, const basis& bs);
MSVCDLL friend super_op left(const matrix& mx); 
MSVCDLL friend super_op right(const gen_op& Op); 	// LOp*Op1 = Op1*Op
MSVCDLL friend super_op right(const matrix& mx); 	// LOp*Op1 = Op1*mx
MSVCDLL friend super_op right(const matrix& mx, const basis& bs);

// ----------------------------------------------------------------------------
//                          Commutation Superoperators
// ----------------------------------------------------------------------------

/*                           t
    LOp = Op (X) E - E (x) Op          where       LOp*Op1 = [Op, Op1]       */

MSVCDLL friend super_op commutator(const gen_op& Op); 	// LOp*Op1 = [Op, Op1]
MSVCDLL friend super_op commutator(const matrix& mx);	// LOp*mx = [mx, mx1] 

// ----------------------------------------------------------------------------
//                      Double Commutation Superoperators
// ----------------------------------------------------------------------------
 
// sosi - these still need to be made valid for composite Liouville space

/*                                       T           T           T     T 
          LOp = Op1*Op2 (X) E - Op1 X Op2 - Op2 X Op1 + E (X) Op1 * Op2
   where
                           LOp*Op3 = [Op1,[Op2, Op3]]

   These functions return LOp, which acts according to the equations below
   (respectively) where X is any arbitrary operator (or array) in Hilbert
   space. The other operators are they from which LOp is built. 

	1.) LOp*X =   [Op,  [Op,  X]]
	2.) LOp*X =   [mx,  [mx,  X]]
	3.) LOp*X =   [Op1, [Op2, X]]
	4.) LOp*X = z*[Op1, [Op2, X]]
	5.) LOp*X =   [mx1, [mx2, X]]                                       */
 
MSVCDLL friend super_op d_commutator(const gen_op& Op, const complex& z);
MSVCDLL friend super_op d_commutator(const matrix& mx);
MSVCDLL friend super_op d_commutator(const gen_op& Op1, const gen_op& Op2);
MSVCDLL friend super_op d_commutator(const gen_op& Op1, const gen_op& Op2,
                                                             const complex& z);
MSVCDLL friend super_op d_commutator (const matrix& mx1, const matrix& mx2);

/*                                       T           T          T     T
         LOp = (mx1)(mx2) X E - mx1 X mx2 - mx2 X mx1 + E X (mx1 )(mx2 )     */

// ----------------------------------------------------------------------------
//                    Unitary Transformation Superoperators
// ----------------------------------------------------------------------------
 
/*                     *                                              -1
      LOp = Op (X) Op             where        LOp*Op1 = Op * Op1 * Op       */
 
MSVCDLL friend super_op U_transform(const gen_op& Op);
MSVCDLL friend super_op U_transform(const matrix& mx);

// ----------------------------------------------------------------------------
//                           Projection Superoperators
// ----------------------------------------------------------------------------
 
MSVCDLL friend super_op project(const gen_op& Op);

	// Input		Op    : General operator
	// Output		LOp   : Projection superoperator
	//				LOp*Op1 = z*Op, w/ Op1 = z*Op + .....
	//					     t
	//				LOp = Op(x)Op


MSVCDLL friend super_op project(const matrix& mx);

	// Input		mx    : Matrix
	// Output		LOp   : Projection superoperator
	//				LOp*Op1 = z*mx, w/ Op1 = z*mx + .....
	//					     t
	//				LOp = mx(x)mx

// ----------------------------------------------------------------------------
//                          Exponential Superoperators
// ----------------------------------------------------------------------------

MSVCDLL super_op exp() const;

        // Input                LOp   : Superoperator (this)
        // Return               ExpLOp: Exponential of LOp
        //                              ExpLOp = exp(LOp)
        // Note                       : Exponential output in EBR of LOp
        // Note                       : L0p's EBR is generated herein


MSVCDLL super_op exp(const complex& t, double cutoff=1.e-12) const;

        // Input                LOp   : Superoperator (this)
	//			t     : Exponential factor
	//			cutoff: Exponential factor roundoff
        // Return               ExpLOp: Exponential of LOp
        //                              ExpLOp = exp(t*LOp)
        // Note                       : Exponential output in EBR of LOp
        // Note                       : L0p's EBR is generated herein
	// Note			      : Value of t is considered 0 if
	//				it's magnituded is less than cutoff

  
MSVCDLL friend super_op exp(const super_op& LOp1);

	// Input		LOp1  : Superoperator
        // Return		LOp   : Exponential of LOp1
	//				LOp = exp(LOp1)
        // Note			      : Computed in EBR of LOp1

  
MSVCDLL friend super_op exp(const super_op& LOp1, const complex& t);

	// Input		LOp1  : Superoperator
        // Return		LOp   : Exponential of LOp1
	//				LOp = exp(LOp1*t)
        // Note			      : Computed in EBR of LOp1
  

MSVCDLL friend super_op pow(const super_op& LOp, int power);

	// Input		LOp   : Superoperator (this)
	//			power : A power factor
        // Return		LOp1  : Initial superoperator taken
	//				to a particular power
        // Note			      : Computed in EBR of LOp

// ____________________________________________________________________________
// F                  SUPER OPERATOR BASIS MANIPULATIONS
// ____________________________________________________________________________

// F_list set_EBR		- Set superoperator eigenbasis
// F_list set_HBR		- Set superop to default Liouville space basis 
// F_list set_DBR	        - Set superop default Hilbert space basis 
// F_list LOp_base	        - Set superoperator basis

MSVCDLL void set_EBR() const;

	// Input		LOp	: General operator (this)
	// Output		none	: LOp set to EBR 
        // Note                         : Not a constant function, both
        //                                matrix and basis may change

	
MSVCDLL void set_HBR() const;

	// Input		LOp	: Superoperator (this)
	// Output		none	: LOp set to default LS basis
        // Note                         : Not a constant function, both
        //                                matrix and basis may change
 
         
MSVCDLL void set_DBR() const;
 
        // Input                LOp  : Superoperator (this)
        // Output               none : LOp set into a default Hilbert
        //                             space basis.


MSVCDLL void LOp_base(const super_op& LOp1) const;

	// Input		LOp1 : Superoperator.
	// 			LOp2 : Superoperator.
	// Return		none : Superoperator LOp2 put
	//			       into the basis of LOp1


MSVCDLL void LOp_Hbase(const super_op& LOp1, int warn=0) const;

	// Input		LOp  : Superoperator (this).
	//			LOp1 : Superoperator.
	//			warn : Flag if warnings desired (off)
	// Return		none : Superoperator LOp is put
	//			       into the Hilbert space basis of LOp1
	// Note			     : This function performs similarity
	//			       transformations in Liouville space.
	//			       It is likely to be computationally
	//			       intensive, avoid it if possible.

	// Input		LOp  : Superoperator.
	// 			Op   : General operator.
	// Return		none : General operator Op put
	//			       into the Hilbert space basis
	//			       of superoperator LOp

MSVCDLL void LOp_base(const    gen_op& Op) const;
MSVCDLL void SetHSBaseOf(const gen_op& Op) const;

// ____________________________________________________________________________
// G                   SUPER OPERATOR AUXILIARY FUNCTIONS
// ____________________________________________________________________________

///Center Superoperator Auxiliary Functions


MSVCDLL int HS()   const;			// Operator Hilbert space
MSVCDLL int size() const;			// Operator Liouville space
MSVCDLL int dim()  const;			// Operator Liouville space
MSVCDLL int LS()   const;			// Operator Liouville space

	// Input		LOp  : Superoperator (this).
	// Return		LSp  : Liouville space dimension of LOp
	// F_list size	  	     - Get superoperator Liouville dimension

  
MSVCDLL void eigenvalues(int nc = 4, int ri=0) const;

	// Input		LOp   : Superoperator (this)
	//			nc    : Number of columns to print in
	//			ri    : Flag for Real(0), Complex(1), Im(-1)
        // Return		None  : LOp eigenvalues sent to
	//				standard output
        // Note			      : Sets EBR of LOp
	// F_list eigenvalues	      - Output superoperator eigenvalues


// ____________________________________________________________________________
// H                 SUPER OPERATOR COMPONENT MANIPULATIONS
// ____________________________________________________________________________

///Center Superoperator Component Manipulation Functions

// -------------------------- Matrix Manipulations ----------------------------

MSVCDLL matrix Mx()     const;
MSVCDLL matrix get_mx() const;

        // Output               mx    : Matrix rep of LOp in Liouville space

MSVCDLL void put_mx(const matrix& mx1);

	// Input		LOp   : Superoperator (this)
        // 			mx1   : Matrix  in Liouville space
        // Output		none  : LOp altered to have mx as
	//				its representation
	// F_list put_mx              - Set superoperator matrix 

// -------------------------- Basis Manipulations -----------------------------

MSVCDLL basis Bs()       const;
MSVCDLL basis get_basis() const;

	// Input		LOp   : Superoperator (this)
        // Output               Hbs   : Hilbert space basis of LOp
	// Note			      : LOp static, basis is copied
	// F_list get_basis           - Get LOp Hilbert space basis


MSVCDLL void put_basis(const basis& Hbs);

	// Input		LOp   : Superoperator (this)
        // 			Hbs   : Basis in Hilbert space
        // Output		none  : LOp altered to have Hbs1 as
	//				its current basis without a
	//			        change in the LOp matrix
	// F_list put_basis           - Set superoperator basis

	// Input		LOp   : Superoperator (this)
        // Output               Lbs   : Hilbert space basis of LOp
	// Note			      : LOp static, basis is copied
	// F_list get_Lbasis          - Get LOp Liouville space basis

MSVCDLL basis LBs()        const;
MSVCDLL basis get_Lbasis() const;


MSVCDLL void put_Lbasis(const basis& Lbs);

	// Input		LOp   : Superoperator (this)
        // 			Lbs   : Basis in Liouville space
        // Output		none  : LOp altered to have Lbs as
	//				its current Liouville basis
	//			        without change in the LOp matrix
	// F_list put_Lbasis          - Set superoperator Liouville basis


// -------------------- Individual Element Manipulations ----------------------

MSVCDLL complex operator() (int row, int col) const;

	// Input		LOp   : Superoperator (this)
        //                      row   : Row index
        //                      col   : Column index
        // Output               z     : Value of <row|LOp|col>
	// Note			      : LOp remains static, so LOp(i,j) = z
	//				will not work, unlike the matrix analog
	// F_list ()                  - Get superoperator matrix element


MSVCDLL void put(int row, int col, const complex& z);

	// Input		LOp   : Superoperator (this)
        //                      row   : Row index
        //                      col   : Column index
	//			z     : Complex number
        // Output               none  : Value of <row|LOp|col> set to z
	// F_list put                 - Set superoperator matrix element


MSVCDLL complex get(int row, int col) const;

	// Input		LOp   : Superoperator (this)
        //                      row   : Row index
        //                      col   : Column index
        // Output               z     : Value of <row|LOp|col>
	// Note			      : LOp remains unchanged, element copy 	
	// F_list get                 - Get superoperator matrix element


// ____________________________________________________________________________
// I                        CLASS SUPER OPERATOR CHECKS
// ____________________________________________________________________________

/* This set of functions insures that the LOp/Op/mx/bs is properly square &/or
   of compatible dimension to mix with another LOp/Op/mx/bs.  All return TRUE
   if dimensions look OK.  The functions take a warning level flag which will
   be used to decide what happens under a FAIL.  If the flag is 0 then nothing
   will result.  If the flags is 1 a non-fatal error is output.  If the flag
   is >1 a fatal error is output and program execution stopped.              */  

MSVCDLL bool checkLOp(const super_op& LOp1,                 int warn=2) const;
MSVCDLL bool checkLOp(const gen_op& Op,                     int warn=2) const;
MSVCDLL bool checkLOp(const matrix& mx,                     int warn=2) const;
MSVCDLL bool checkLOp(const matrix& mx1, const matrix& mx2, int warn=2) const;
MSVCDLL bool checkLOp(const matrix& mx,  const basis& bs,   int warn=2) const;
MSVCDLL bool checkLOp(int row,           int col,           int warn=2) const;

// ____________________________________________________________________________
// J                       CLASS SUPER OPERATOR TESTS
// ____________________________________________________________________________

///Center Test Functions

MSVCDLL void status() const;

	// Input		LOp   : Superoperator (this)
	// Output		void  : Outputs superoperator status


MSVCDLL int operator == (const super_op& LOp1);

	// Input		LOp   : Superoperator (this)
        //                      LOp   : Superoperator
        // Output               T_F   : TRUE if LOp = LOp1


MSVCDLL int below(double d) const;

	// Input		LOp   : Superoperator (this)
        //                      d     : Number
        // Output               T_F   : TRUE if all norm(<i|LOp|j>) <= d


// ____________________________________________________________________________
// K                  CLASS SUPER OPERATOR I/O FUNCTIONS
// ____________________________________________________________________________

///Center Superoperator I/O Functions
// F_list print              - ASCII output of superoperator to filestream
// F_list <<                 - Standard output of superoperator

// -------------------------- ASCII Output Functions --------------------------
 
/*         Input                LOp  : Superoperator (this)
                                ostr : Output ASCII file stream
                                flag : Flag for amount of output
           Return               ostr : LOp is set to the output stream       */

MSVCDLL std::ostream& print(std::ostream& ostr, int flag=0) const;
MSVCDLL friend std::ostream& operator << (std::ostream& ostr, const super_op& LOp);

// -------------------------- Binary Output Functions -------------------------
 
/*         Input                LOp  : Superoperator (this)
                                fn   : Output binary filename
				fp   : File stream (pointing at LOp spot)
           Return               void : LOp is written in binary format to
				       the file or filestream.               */

MSVCDLL void           write(const std::string& fn)  const;
MSVCDLL std::ofstream& write(std::ofstream& fp) const;

// --------------------------- Binary Input Functions -------------------------
 
/*         Input                LOp  : Superoperator (this)
                                fn   : Input binary filename
				fp   : File stream (pointing at LOp spot)
				Op   : Operator used to set working HS basis
				LOp  : Superoperator used to set working bases
           Return               void : LOp is read in binary format to
				       the file or filestream.
	   Note                      : Superoperators read in this fashion will
				       NOT share the same bases even if they
				       did so when written.  This must be done
                                       explicitly if desired.                */

MSVCDLL void           read(const std::string& fn);
MSVCDLL void           read(const std::string& fn, const gen_op&   Op);
MSVCDLL void           read(const std::string& fn, const super_op& LOp1);
MSVCDLL std::ifstream& read(std::ifstream& fp);

// ____________________________________________________________________________
// L                 CLASS SUPEROPERATOR LEFTOVER FUNCTIONS
// ____________________________________________________________________________


MSVCDLL friend super_op Hsuper(const gen_op& Heff);

	// Input		Heff  : Effective Hamiltonian (Hz)
	// Output		LOp   : Hamiltonian commutation superoperator
	//
	//       <a,aa|LOp|b,bb> = del  * del     <b|H   |b> - <bb|H   |bb>
	//			      ab     aa,bb    eff           eff
	//
	// Note			      :	Lop is returned in angular frequency
	//				units


MSVCDLL friend super_op HsuperX(const gen_op& Heff);

        // Input                Heff  : Effective Hamiltonian (Hz)
        // Output               LOp   : Hamiltonian commutation superoperator
        //
        //       <a,aa|LOp|b,bb> = del  * del     <b|H   |b> - <bb|H   |bb>
        //                            ab     aa,bb    eff           eff
        //
        //
        // Note                       : LOp is returned in angular frequency
        //                              units

}; 

#endif 						// SuperOp.h

