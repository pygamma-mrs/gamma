/* SpinOp.h *****************************************************-*-c++-*-
**									**
**                               G A M M A       			**
**									**
**      Spin Operators                 			Interface	**
**									**
**      Copyright (c) 1991						**
**      Z.L. Madi, S. Smith, T. Levante					**
**      Eidgenoessische Technische Hochschule				**
**      Labor fur physikalische Chemie					**
**      8092 Zurich / Switzerland					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
**  Description								**
**									**
**  Class spin_op embodies quantum mechanical spin operators in either	**
**  a single spin or composite spin Hilbert space.  The selectivity is	**
**  defined by a list of spins typically specified by a spin system	**
**  (in class spin_sys, see SpinSys files).				**
**									**
**  The composite space matrix representation is in the natural basis	**
**  constructed from the direct products of single spin bases.		**
**									**
**  Spin operators may internally exists as a list of single spin	**
**  arrays whose tensor products form the spin operator, or as a single **
**  composite space array, or both. 					**
**								 	**
*************************************************************************/

#ifndef SpinOp_h_			// Is file already included?
#  define SpinOp_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/SpinSys.h>		// Know about base spin systems
#include <Matrix/matrix.h>		// Know about matrices

class spin_op 
  {
  int      nspins;			// Assoicated number of spins
  int*     Hspaces;			// Associated spin Hilbert spaces
  int*     spinflags;			// Associated spin flags
  matrix*  pr;				// Array of single spin matrices  (ptr)
  mutable  matrix mx;			// Composite Hilbert space matrix (ptr)
  static   matrix FSmx;			// Default full space matrix

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i                   CLASS SPIN OPERATOR ERROR HANDLING
// ____________________________________________________________________________

        // Input                SOp     : Spin operator (this)
        //                      eidx    : Error index
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error message
        //                                Program execution stopped (fatal)

         void SOperror(int    eidx, int noret=0) const;
volatile void SOpfatality(int eidx)              const;

// ____________________________________________________________________________ 
// ii                CLASS SPIN PRIVATE AUXILIARY FUNCTIONS
// ____________________________________________________________________________


void blow_up() const;

	// Input		SOp  : Spin operator (this)
	// Return		none : Spin operator SOp is produced and
	//			       stored in the full Hilbert space
	// Note		             : Does nothing if full Hilbert space
	//			       representation exists

void blow_up(matrix_type t) const;
        // Input                SOp     : Spin operator (this)
        // Return               void    : The full Hilbert space matrix
        //                                representaion of SOp is produced
        //                                if it does not already exist.
        // Note                         : A 1x1 I matrix is the default "mx" &
        //                                flags that the true mx isn't present
        // Note                         : Nothing is done if full mx exists
        // sosi                         : Hopefully the tensor_product function
        //                                takes into account matrix structure?
        //                                but this forces the issue! 
     
// ____________________________________________________________________________
// iii             CLASS SPIN OPERATOR CHECKING FUNCTIONS
// ____________________________________________________________________________
 
 
int checkSpin(int i, int warn=2) const;

        // Input                SOp     : Spin operator (this)
        //                      warn    : Warning level
        //                                      0 = no warning
        //                                      1 = non-fatal  warning
        //                                     >1 = fatal warning
        // Output               TF      : True if the spin exists


int checkSys(const spin_op& SOp1, int warn=1) const;

        // Input                SOp     : Spin operator (this)
        //                      SOp1    : Another spin operator
        //                      warn    : Warning level
        //                                      0 = no warning
        //                                      1 = non-fatal  warning
        //                                     >1 = fatal warning
        // Output               TF      : True if the systems associate with
        //                                the two operators match

// ____________________________________________________________________________
// iv             CLASS SPIN OPERATOR SETUP FUNCTIONS
// ____________________________________________________________________________

void CopySubSpaces(matrix* prmxs);
void CopySubSpaces(const spin_op& SOp);
void CopySpinFlags(matrix* prmxs);
void CopySpinFlags(const spin_op& SOp);
void ZeroSpinFlags();
void BlendSpinFlags(const spin_op& SOp);
void CopySubArrays(matrix* prmxs);
void CopySubArrays(const spin_op& SOp);
void DelSubArrays();

int CopyFullMx(const spin_op& SOp1, int warn=2) const;
 
        // Input                SOp     : Spin operator (this)
        //                      SOp1    : Another spin operator
        //                      warn    : How to handle errors
        //                                      0 = no warning
        //                                      1 = output warning
        //                                     >1 = fatal problem
        // Return               TF      : This routine copies the full
        //                                space array of SOp1 to SOp2
        //                                (no size checking!!) and returns
        //                                TRUE if copy is proper!

 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

public:

// ____________________________________________________________________________
// A              CLASS SPIN OPERATOR CONSTRUCTORS, DESTRUCTOR
// ____________________________________________________________________________

/* These constructors set up a new spin operator.  There are several ways to 
   make a spin operator as listed below:
 
             Input Arguments                Resulting Spin Operator
          ---------------------    -----------------------------------------
                   -               Empty Spin Operator
             spins, prmxs          SOp with N subspace & N subspace arrays
                  SOp              SOp that is a duplicate of the on input   */

MSVCDLC spin_op();
MSVCDLC spin_op(int spins, matrix* prmxs);
MSVCDLC spin_op(const spin_op& SOp);
MSVCDLC ~spin_op();
spin_op& operator= (const spin_op& SOp);
 
// ____________________________________________________________________________
// B                 SPIN OPERATOR - SPIN OPERATOR FUNCTIONS
// ____________________________________________________________________________

/* These functions allow for simple mathematical operations between two spin
   operators.  This includes addition, subtraction, multiplication. There is
   one unary function as well, negation.

   Operator Arguments      Result        Operator Arguments    Result
   -------- --------- -----------------  -------- --------- -------------------
      -        SOp    Returns -SOp          -=    SOp,SOp1  SOp1 subt. from SOp
      +     SOp,SOp1  Returns SOp+SOp1       *    SOp1,SOp2 Returns SOp1*SOp2
     +=     SOp,SOp1  SOp1 added to SOp     *=    SOp,SOp1  SOp mult into SOp1
      -     SOp1,SOp2 Returns SOp1-SOp2                                      */

MSVCDLL        spin_op operator-  ()                    const;
MSVCDLL        spin_op operator+  (const spin_op& SOp1) const;
MSVCDLL        spin_op &    operator+= (const spin_op& SOp1);
MSVCDLL friend spin_op          operator-  (const spin_op& SOp1, const spin_op& SOp2);
MSVCDLL        spin_op &    operator-= (const spin_op& SOp1);
MSVCDLL friend spin_op          operator*  (const spin_op& SOp1, const spin_op& SOp2);
MSVCDLL        spin_op &    operator*= (const spin_op& SOp1);
 
// ____________________________________________________________________________
// C                   SPIN OPERATOR - SCALAR FUNCTIONS
// ____________________________________________________________________________

/* These functions allow for two mathematical operations between a scalar &
   a spin operators, multiplication & division. 

 Operator Arguments      Result        Operator Arguments    Result
 -------- --------- -----------------  -------- --------- -------------------
    *      z,SOp    Returns z*SOp          *=    SOp,z    SOp multiplied by z
    *      SOp,z    Returns z*SOp          *=    SOp,d    SOp multiplied by d
    *      d,SOp    Returns d*SOp          /     SOp,z    Returns (1/z)*SOp
    *      SOp,d    Returns SOp*d          /     SOp,d    Returns (1/d)*SOp
    /=     SOp,d    SOp mult. by (1/d)     /=    SOp,z    SOp mult. by (1/z) */

MSVCDLL friend spin_op       operator*  (const spin_op& SOp, const complex& z);
MSVCDLL friend spin_op       operator*  (const spin_op& SOp, double d);
MSVCDLL friend spin_op       operator*  (const complex& z,   const spin_op& SOp);
MSVCDLL friend spin_op       operator*  (double d,           const spin_op& SOp);
MSVCDLL        spin_op & operator*= (const complex& z);
MSVCDLL        spin_op & operator*= (double d);
MSVCDLL friend spin_op       operator/  (const spin_op& SOp, const complex& z);
MSVCDLL friend spin_op       operator/  (const spin_op& SOp, double d);
MSVCDLL        spin_op & operator/= (const complex& z);
MSVCDLL        spin_op & operator/= (double d);

// ____________________________________________________________________________
// D                  SPIN OPERATOR - MATRIX FUNCTIONS
// ____________________________________________________________________________

MSVCDLL matrix   get_mx() const;
MSVCDLL operator matrix() const;

	// Input		SOp  : Spin operator (this)
	// Return		mx   : Matrix which is the
        //	                       equivalent of the input spin operator

// ____________________________________________________________________________
// E                  SPIN OPERATOR ASSOCIATED FUNCTIONS
// ____________________________________________________________________________

/* "useful" functions, probably unneeded when class spin_op is gone!         */

MSVCDLL spin_op exp()     const; 		// Exponentation
MSVCDLL spin_op adjoint() const; 		// Adjoint
MSVCDLL complex trace()   const; 		// Trace

MSVCDLL friend spin_op exp(const spin_op& SOp); 		// Exponentation
MSVCDLL friend spin_op adjoint(const spin_op &SOp); 		// Adjoint
MSVCDLL friend complex trace(const spin_op &SOp); 		// Trace

/*                 These are common access functions                         */

MSVCDLL int spins( )    const;				// Number of spins
MSVCDLL int refs( )     const;				// Full space mx refs
MSVCDLL int refs(int i) const;				// Sub-space mx refs
MSVCDLL int HS( )       const;				// Full Hilbert space

// ____________________________________________________________________________
// F                       SPIN OPERATOR I/O FUNCTIONS
// ____________________________________________________________________________

	// Input		SOp  : Spin operator (this)
	//			ostr : An output stream
	//			full : Flag for how much output
	// Return		void : SOp is placed into the output stream


MSVCDLL        void       print(std::ostream& ostr, int full=0) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const spin_op& SOp);
   
        // Input                SOp     : Spin operator (this)
        //                      full    : Flag for amount of output
        // Output               void    : Outputs SOp status

MSVCDLL void status(int full=1) const;
 
// ____________________________________________________________________________
// G                           SPIN OPERATOR KLUDGE
// ____________________________________________________________________________

// This function forces a Spin operator's matrix to take on a specific
// matrix structure.  Its a kludge because so far GAMMA's matrices don't know
// how to keep a proper matrix structure in the "blow up" function.  When this
// is fixed in the matrix classes, remove this function and all places where
// its used because it will be no longer necessary.

MSVCDLL void FaxisStruct(char axis) const;

};
 
#endif								// SpinOp.h
