/* basis.h ******************************************************-*-c++-*-
**                                                                      **
**                                 G A M M A                            **
**                                                                      **
**      Basis                                      Interface 		**
**                                                                      **
**      Copyright (c) 1991, 1992, 1997                                  **
**                                                                      **
**      Scott Smith                                                     **
**      Copyright (c) 1996                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      Tilo Levante                                                    **
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
** The class basis provides arrays that defines a basis in GAMMA for    **
** operators and superoperators.  These are static private structures   **
** that are difficult to alter.  The class currently supports composite **
** spaces, i.e. bases that span multiple sub-spaces.                    **
**                                                                      **
** Likely, each bais is	just a matrix residing in the Hilbert/Liouville	**
** space of an associated operator/superoperator.  However, the class	**
** will also handle composite spaces which are direct products, i.e.	**
** the basis has defined subspaces with specific dimensions. The number	**
** of the dimensions and their individual sizes are also maintained by	**
** this class.								**
**                                                                      **
*************************************************************************/

///Chapter Class Basis
///Section Overview
///Body    The class Basis handles basis arrays which
///        are used in transforming operators and
///	   superoperators between their different basis
///	   representations.
///Section Available Basis Functions

#ifndef   Gbasis_h_			// Is file already included?
#  define Gbasis_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/matrix.h>		// Knowledge of matrices
#include <vector>			// Know libstdc++ STL vectors

class basis : private matrix
  { 
  int  nc;				// Number of components
  int* ncd;		 		// Dimensions of components
  std::string bname;			// Basis name


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i                       CLASS BASIS ERROR HANDLING
// ____________________________________________________________________________

/*      Input		       bs      : Basis (this)
                               eidx    : Error index
                               noret   : Flag for linefeed (0=linefeed)
                               pname   : Added error message for output
       Output                  none    : Error message output
                                         Program execution stopped if fatal

  The following error messages use the defaults set in the Gutils package

               Case                          Error Message

               (0)                     Program Aborting.....
               (1)                     Problems With Input File Stream
               (2)                     Problems With Output File Stream
               default                 Unknown Error                        */

         void BSerror(int eidx,                           int noret=0) const;
         void BSerror(int eidx, const std::string& pname, int noret=0) const;
volatile void BSfatal(int eidx)                                        const;
volatile void BSfatal(int eidx, const std::string& pname)              const;

 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
public:
  
// ___________________________________________________________________________
// A                   BASIS CONSTRUCTION, DESTRUCTION
// ___________________________________________________________________________

///Center Basis Algebraic
/// F_list basis 	- Constructor
/// F_list ~	 	- Destructor

// ---------------------------------------------------------------------------
//                           Simple Constructors
// ---------------------------------------------------------------------------

/* The produce default bases, essentially identity matrices either in an 
   single space or in a composite space.                                    */

MSVCDLC basis();
MSVCDLC basis(int dim);
MSVCDLC basis(const std::vector<int> dims);

MSVCDLC basis(const matrix& mx, int nc=1, int* ncd=NULL);
MSVCDLC basis(const basis&  bs, const matrix& mx);
MSVCDLC basis(const basis&  bs);

// ---------------------------------------------------------------------------
//                         Destruction And Assignment
// ---------------------------------------------------------------------------

MSVCDLC        ~basis();
MSVCDLL basis& operator=(const basis& bs);
MSVCDLL basis& operator=(const matrix& mx);

// ___________________________________________________________________________
// B                       BASIC BASIS FUNCTIONS
// ___________________________________________________________________________

///Center Basic Basis Functions 

	// Input		: Basis (this)
	// Output		: Dimension of the basis
	/// F_list dim	 	- Basis dimension

MSVCDLL int         size() const;
MSVCDLL int         dim()  const;
MSVCDLL std::string name() const;
MSVCDLL void        name(const std::string& nm);

// ____________________________________________________________________________
// C            BASIC BASIS FUNCTIONS SUPPORTING SUBSPACES
// ____________________________________________________________________________
 

MSVCDLL int dim_LS() const;

	// Input		: Basis (this)
	// Output	   int  : Dimension of LS derived from this basis

MSVCDLL int sub_N() const;

	// Input		: Basis (this)
	// Output		: Number of multi_sys components in rep

MSVCDLL int sub_dim(int ic) const;

	// Input	       int  : Component of multi_sys
	// Output	       int  : Dimension associated with ic

MSVCDLL int sub_anchor(int ic) const;

	// Input	       int  : Component of multi_sys
	// Output	       int  : X,Y (X=Y) of the upper left corner of
	//			    : corresponding block 

MSVCDLL int sub_anchor_LS (int ic) const;

	// Input	       int  : Component of multi_sys
	// Output	       int  : X,Y (X=Y) of the upper left corner of
	//			    : corresponding block in LS

MSVCDLL int which_sub_LS(int i) const;

        // Input                 int  : Number of basis state 
        // Output                int  : Number of multi_sys component to which
	//			      : this state belongs
    

MSVCDLL friend basis defLSbasis(const basis& bs);

        // Input            bs  : A (Hilbert space) basis
        // Output           Lbs : A default Liouville space basis
        //                        which has a sub-space Liouville structure
        //                        that coincides with bs Hilbert sub-space(s)


MSVCDLL friend basis defbasis(const basis& bs);

        // Input            bs  : A (Hilbert space) basis
        // Output           hbs : A default basis which has a sub-space
	//			  structure that coincides with bs


// ____________________________________________________________________________
// D                    BASIS ARRAY & TRANSFORMATION ACCESS
// ____________________________________________________________________________

/// F_list U	 	- Basis conversion matrix
/// F_list convert 	- Matrix conversion into new basis
/// F_list convert_back	- Matrix conversion out of basis
 
/* These functions allow access to the basis transformation array and apply
   that transformation to an input matrix.

   Function   Arguments                           Result 
 ------------ --------- ------------------------------------------------------- 
      U         ----    Returns the basis array as a matrix
   convert       mx     Applies basis transformation (put mx in bs)
 convert_back    mx     Applies inverse basis transform (take mx out of bs)

                           -1             bs          *t
        convert:   mx1 = bs  * mx * bs =======> = (bs)  * mx * bs
                                       unitary

                                     -1   bs                   *t
   convert_back:   mx1 = bs * mx * bs   =======> bs * mx * (bs)
                                        unitary


The two "convert" functions below are at the moment very poor. This is because
GAMMA quickly looses track of whether an array (basis) is unitary... forcing
us to to a time-consuming check to see if it is.                             */
 
MSVCDLL matrix U()                            const;
MSVCDLL matrix get_matrix()                   const;
MSVCDLL matrix get_mx()                       const;
MSVCDLL matrix convert(const      matrix& mx) const;
MSVCDLL matrix convert_back(const matrix& mx) const; 

// ___________________________________________________________________________
// E                     MORE COMPLEX BASIS FUNCTIONS
// ___________________________________________________________________________

// ***************************** tensor product *******************************

// These functions are also "kludges" in the sense that better matrix handling
// would make them totally unnecessary.  A scan of the GAMMA sources should
// reveal that these functions are rarely used (I Hope)


MSVCDLL friend basis tensor_product(const basis& bs1, const basis& bs2);

        // Input            bs1 : A basis
        //                  bs2 : A second basis
        // Output           pdt : A basis which is the tensor product
        //                        of the two input bases
        // Note                 : This mirrors the handling of tensor
        //                        products by class matrix, but because
        //                        of private inheritance in class basis
        //                        it needs to be in this class too.

        //                           pdt        =   bs1 (x) bs2

        //                       (m*o x n*p)       (mxn)   (oxp)

        //                    <i*o+k|pdt|j*p+l> = <i|bs1|j><k|bs2|l>
  
// ___________________________________________________________________________
// F                    BASIS COMPARISON/TEST FUNCTIONS
// ___________________________________________________________________________

///Center Basis Comparison and Test Functions 
/// F_list ==		  - Basis comparison
/// F_list !=             - Basis comparison
/// F_list isDefualtBasis - Default basis check
/// F_list refs		  - Number of basis references

MSVCDLL bool operator==(const basis& bs2) const;
MSVCDLL bool operator!=(const basis& bs2) const;
MSVCDLL bool isDefaultBasis() const;


MSVCDLL int refs() const;

	// Input		: Basis (this)
	// Output	    int : Returns number of references to the
        //			  basis matrix
  

MSVCDLL bool check(const basis& bs1) const;

        // Input                : Bases to be compared, (this) and bs1
        // Outptut              : True if matrix dimensions, number of
        //                        components in the representation and
        //                        dimensions associated with components
        //                        coincide
  
// ___________________________________________________________________________
// G                         BASIS I/O FUNCTIONS
// ___________________________________________________________________________

///Center Basis Input Output Functions 
///F_list <<		- Basis standard output
///F_list write		- Basis output to binary

// ----------------------- ASCII Output Functions ----------------------------

/*              Input           bs   : Basis (this)
                                ostr : Output ASCII file stream
                                full : Flag for amount of output
                Return          void : bs is sent to the output stream       */  

MSVCDLL        std::ostream& print(std::ostream& ostr, int full=1) const;
MSVCDLL friend std::ostream& operator<<  (std::ostream& ostr, const basis& bs);

// ------------------------ Binary Output Functions ---------------------------
 
/*              Input           bs   : Basis (this)
                                fn   : Input binary filename
                                fp   : File stream (pointing at bs
                Return          void : Op is written to either the
                                       specified file or filestream.
                Note                 : Output format is partially set by
                                       class matrix (matrix typing)          */
  
MSVCDLL void           write(const std::string& fn) const;
MSVCDLL std::ofstream& write(std::ofstream& fp)     const;

// ------------------------ Binary Input Functions ----------------------------
 
/*              Input           bs   : Basis (this)
                                fn   : Input binary filename
                                fp   : File stream (pointing at bs)
                Return          void : bs is read in from either the
                                       specified file or filestream.         */  
MSVCDLL void           read(const std::string& fn);
MSVCDLL std::ifstream& read(std::ifstream&     fp);

    
// ____________________________________________________________________________
// H                         BASIS TESTING FUNCTIONS
// ____________________________________________________________________________

/* Since a Basis is fundamental to both operators and superoperators in GAMMA
   I have added these tests to provide a means to insure they are working
   correctly. It is paramount that we know how to properly convert an operator
   (or superoperator) between bases. In order for that to occur we MUST be able
   to obtain the inverse of the basis transformation array. In certains cases
   this is just the adjoint (if the basis is unitary) but often times we will
   not be so fortunate, such as when dealing with relaxation superoperator.
   A good test would then be U * inv(U) = 1. That will be proivded here.     */

MSVCDLL double TestBasis(int pf=0) const;

};

#endif						// basis.h
