/* GenOp.h ******************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      General Operator                               Interface        **
**                                                                      **
**      Copyright (c) 1990, 1991, 1992                                  **
**      Scott Smith                                                     **
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
**  Class gen_op defines general operators for GAMMA in C++.  It        **
**  contains all of all the algebraic operations (+, -, *, /), most     **
**  useful complex functions (exp, log) and input/output routines.      **
**                                                                      **
*************************************************************************/

#ifndef   _GenOp_h_			// Is file already included?
#  define _GenOp_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/matrix.h>		// Include matrices
#include <HSLib/GenOpRep.h>		// Include single representations
#include <HSLib/SpinOp.h> 		// Include until spin_op -> tp_matrix!
#include <vector>			// Include stdlibc++ STL vectors

//forward declaration
class gen_op;
MSVCDLL complex proj(const gen_op &Op1, const gen_op &Op2, int norm=0);  //**

class gen_op : public std::vector <genoprep>
{
    
// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------
 
  mutable genoprep*   DBR;		// Op in Default Basis Rep.	(DBR)
  mutable genoprep*   EBR;		// Op in EigenBasis Rep.	(EBR)
  mutable genoprep*   WBR;		// Op in Working Basis Rep.	(WBR)
          std::string OpName;		// Name of Operator     
  static  int         MaxReps;		// Maximum allowed number of reps
  static  int         DBPr;		// Default Basis Rep priority
  static  int         EBPr;		// EigenBasis    Rep priority

private:

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                 CLASS GENERAL OPERATOR ERROR HANDLING
// ____________________________________________________________________________

/*      Input               Op      : General operator (this)
                            eidx    : Error index
                            noret   : Flag for linefeed (0=linefeed)
                            pname   : string in message
        Output              none    : Error message output
                                      Program execution stopped if fatal     */

 
void GenOperror(int eidx, int noret=0) const;
void GenOperror(int eidx, const std::string& pname, int noret=0) const;
volatile void GenOpfatality(int eidx) const;
volatile void GenOpfatality(int eidx, const std::string& pname) const;

// ____________________________________________________________________________
// ii               GENERAL OPERATOR COMMON INTERNAL FUNCTIONS
// ____________________________________________________________________________

/* These functions help in manipulating the operator representations.  They     
   individually don't worry about Op integrity so they are PRIVATE and must
   be used with caution. Op integrity demands that if any representations
   exist then WBR should point to one of the, whereas if no representations
   exist WBR,EBR,&DBR should all be NULL.
 
 Function              Action                          Precautions
 -------- ---------------------------------   -----------------------------
  ZeroOp  Zeroes OpReps, Nulls All Pointers   This should always be fine
  setNULL Nulls All Pointers                  OpReps Can Exist, But NO WBR!
  AddRep  Adds New OpRep, Sets WBR EBR etc.   Pre-Set WBR,DBR,EBR as needed  */ 
 
      genoprep* Obegin();
const genoprep* ObeginC() const;
void ZeroOp();
void setNULL() const;
void AddRep(const gen_op& Op1);
void AddRep(const  genoprep& OpRep);
void AddRepM(const genoprep& OpRep) const;

int GetIndex(const genoprep* OpRep) const;
//void push_back(const genoprep& OpRep);
// ____________________________________________________________________________
// iii           INTERNAL OPERATOR REPRESENTATION MANIPULATIONS
// ____________________________________________________________________________

/* These functions look for a particular OpRep and will set the WBR to that
   representation if found.  They do not alter the operator in any other way.
   Note that these are constant functions because only the OpRep that WBR
   points to is changed, i.e. WBR may be set to point to a new OpRep but the
   OpRep isn't altered (except to maybe change a priority value).  This is
   possible because WBR was declared mutable in the class definition.

   Functions  Argument              Result                       Caution
   ---------  --------  --------------------------------  ---------------------
   check_DBR    Op1     If Op1 DBR=WBR, sets Op DBR=WBR   NULL WBR's,DBR's Bad?
   check_EBR    Op1     If Op1 WBR=EBR, sets Op EBR=WBR   NULL WBR's,EBR's Bad?
   check_DBR     -      If Op WBR tests DBR, Set DBR=WBR  WBR MUST NOT BE NULL
   check_EBR   cutoff   If Op WBR tests EBR, Set EBR=WBR  WBR MUST NOT BE NULL
 
   Why are these private? For one thing, it isn't mandatory that a DBR or
   EBR exists for any particular operator (same for WBR in a NULL operator)
   so accessing OpReps through these pointers could lead to big trouble. But
   also, these will over write any existing OpRep pointer without checking.
   Use caution! Mostly these are called after construction, assignment, or
   an operation that leaves the operator with only 1 OpRep, WBR.             */

void check_DBR(const gen_op &Op1)      const;
void check_EBR(const gen_op &Op1)      const;
void check_DBR( )                      const;
void check_EBR(double cutoff = 1.e-12) const;

// ********************* Finding Specific Representations ********************

/* Functions  Argument              Result                       Caution
   ---------  --------  --------------------------------  ---------------------
    FindRep      bs     Pointer to OpRep with basis bs     Return NULL !found
    SetRep       bs     Set WBR to OpRep with basis bs     Return 0 if unable
    LP_rep       -      Pointer to OpRep with low prior.         -----

 Note that LPR, the OpRep of lowest priority, is deleted in delete_rep. We hope
 that higher priority reps are stored at list top from IBR and lower priority
 reps at the bottom of the list. WBR is preferentially NOT chosen as LPR.    */

const genoprep* FindRep(const basis &bs) const;
int             SetRep(const basis &bs)  const;
const genoprep* LP_rep()                 const;

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A             CLASS GENERAL OPERATOR CONSTRUCTORS/DETRUCTOR
// ____________________________________________________________________________

/// Center General Operator Algebraic 
//  F_list = 		     - Constructor
//  F_list ~ 		     - Destructor 
 
/* These constructors set up a new operator.  There are several ways to make
   an operator as listed below:
 
             Input Arguments                Resulting Operator
             ---------------       -----------------------------------------
                   -               Empty Operator, No Representations
                  mx               Op With 1 Representation, mx in DBR
                  SOp              Op With 1 Representation, SOp's mx in DBR
                mx1,mx2            Op With 1 Representation, {mx1,mx2}
                 mx,bs             Op With 1 Representation, {mx,bs}
               *mx,n,*bs           Op With 1 Representation, {MXCmp,BSCmp} 
               N,*mx,*bs           Op With N Representations {mx[i],bs[i]} 
 
   Above, MXCmp & BSCmp are larger composite space arrays formed by taking
   smaller space arrays/bases and placing them on the diagonal of the
   larger array/basis to result in a block-diagonal form.  This supports
   use of multi-system spin systems in GAMMA.                                */
	
MSVCDLC gen_op();
MSVCDLC gen_op(const matrix&  mx);
MSVCDLC gen_op(const spin_op& SOp);
MSVCDLC gen_op(const matrix&  mx1,         const matrix& mx2);
MSVCDLC gen_op(const matrix&  mx,          const basis&  bs);
MSVCDLC gen_op(const gen_op&  Op1);
MSVCDLC gen_op(const std::vector<matrix>& mxc,  const std::vector<matrix>& bsc);
// sosi - the constructor below should be completely replaced by the above one
MSVCDLC gen_op(matrix* mxc,   int nc,      matrix* bsc=NULL);
MSVCDLC gen_op(int N,         matrix* mxs, matrix* bss);
MSVCDLC virtual ~gen_op( );			// Floquet needs this virtual

// ____________________________________________________________________________
// B               OPERATOR FUNCTIONS, OPERATOR WITH OPERATOR 
// ____________________________________________________________________________

/// Center Operator - Operator Functions 
//  F_list + 		     - Operator-Operator Addition
//  F_list += 		     - Operator-Operator Unary Addition
//  F_list - 		     - Operator-Operator Subtraction 
//  F_list - 		     - Operator Negation
//  F_list -= 		     - Operator-Operator Unary Subtraction 
//  F_list * 		     - Operator-Operator Multiplication
//  F_list *= 		     - Operator-Operator Unary Multiply
//  F_list &= 		     - Operator-Operator Reverse Unary Multiply
//  F_list = 		     - Operator-Operator Assignment

/*     Operator    Arguments                           Result 
   --------------- ---------   ------------------------------------------------ 
          +         Op1,Op2    Op1+Op2 in WBR of Op1
          +=        Op1,Op2    Op1 = Op1+Op2 in WBR of Op1
          -            Op1     -Op1 in WBR of Op1
          -         Op1,Op2    Op1-Op2 in WBR of Op1
          -=        Op1,Op2    Op1 = Op1-Op2 in WBR of Op1
          *         Op1,Op2    Op1*Op2 in WBR of Op1
          *=        Op1,Op2    Op1 = Op1*Op2 in WBR of Op1
          &=        Op1,Op2    Op1 = Op2*Op1 in WBR of Op1                   */

MSVCDLL gen_op operator+  (const gen_op& Op)  const;
MSVCDLL gen_op &   operator+= (const gen_op& Op);
MSVCDLL gen_op operator-  (const gen_op& Op)  const;
MSVCDLL gen_op operator-  ()                  const;
MSVCDLL gen_op &   operator-= (const gen_op& Op);
MSVCDLL gen_op operator*  (const gen_op& Op)  const;
MSVCDLL gen_op &   operator*= (const gen_op &Op);
MSVCDLL gen_op &   operator&= (const gen_op &Op);

//MSVCDLL friend gen_op operator- (const gen_op& Op);

MSVCDLL void operator = (const gen_op &Op1);

	// Input		Op1  : General operator.
	// 			Op   : General operator (this).
	// Return		Op   : Operator which is a copy of
	//			       the input operator, Op = Op1.
	// Note		             : Result EXCLUSIVELY in WBR of Op1

// ----------------------------------------------------------------------------
//            OPERATOR FUNCTIONS, OPERATOR WITH SPIN OPERATOR 
// ----------------------------------------------------------------------------

// sosi - added for alpha until spin op --> tp_matrix

MSVCDLL void   operator = (const spin_op &SOp);
MSVCDLL friend gen_op operator + (const gen_op &Op1, const spin_op &SOp);
MSVCDLL friend gen_op operator - (const gen_op &Op1, const spin_op &SOp);
MSVCDLL gen_op &   operator += (const spin_op &SOp);
MSVCDLL gen_op &   operator -= (const spin_op &SOp);

// ____________________________________________________________________________
// C                 OPERATOR FUNCTIONS, OPERATOR WITH MATRIX
// ____________________________________________________________________________

/// Center Operator - Matrix Functions 

/* These functions blend operators with matrices. The arrays must be square
   and match the dimension of the operartor with which it is interacting. All
   arrays are assumed equal to an operator in the default basis. As such, the
   result of matrix-operator interactions is always in the default basis.    */

// F_list + 		     - Operator-Matrix Addition
// F_list += 		     - Operator-Matrix Unary Addition
// F_list - 		     - Operator-Matrix Subtraction
// F_list -= 		     - Operator-Matrix Unary Subtraction
// F_list * 		     - Operator-Matrix Multiplication
// F_list *= 		     - Operator-Matrix Unary Multiplication
// F_list &= 		     - Operator-Matrix Reverse Unary Subtraction
// F_list = 		     - Operator-Matrix Assignment

/* These functions blend operators with matrices. The arrays must be square
   and match the dimension of the operartor with which it is interacting. All
   arrays are assumed equal to an operator in the default basis. As such, the
   result of matrix-operator interactions is always in the default basis. 

       Operator    Arguments                           Result 
   --------------- ---------   ------------------------------------------------ 
          +          Op,mx     Return Op1 = Op+mx            in DBR of Op
          +          mx,Op     Return Op1 = mx+Op            in DBR of Op
          +=         Op,mx     Set Op = Op+mx                in DBR of Op
          -          Op,mx     Return Op1 = Op-mx            in DBR of Op
          -          mx,Op     Return Op1 = mx-Op            in DBR of Op
          -=         Op,mx     Set Op = Op-mx                in DBR of Op
          *          Op,mx     Return Op1 = Op*mx            in DBR of Op
          *          mx,Op     Return Op1 = mx*Op            in DBR of Op
          *=         Op,mx     Set Op = Op*mx                in DBR of Op
          &=         Op,mx     Set Op = mx*Op                in DBR of Op
          =          Op,mx     Set Op = mx                   in DBR          */

MSVCDLL friend gen_op operator + (const gen_op& Op1, const matrix& mx);
MSVCDLL friend gen_op operator + (const matrix& mx,  const gen_op& Op1);
MSVCDLL friend gen_op operator - (const gen_op& Op1, const matrix& mx);
MSVCDLL friend gen_op operator - (const matrix& mx,  const gen_op& Op1);
MSVCDLL friend gen_op operator * (const gen_op& Op1, const matrix& mx);
MSVCDLL friend gen_op operator * (const matrix& mx,  const gen_op& Op1);

MSVCDLL        void   operator += (const matrix& mx);	
MSVCDLL        void   operator -= (const matrix &mx);
MSVCDLL        void   operator *= (const matrix& mx);
MSVCDLL        void   operator &= (const matrix& mx);
MSVCDLL        void   operator  = (const matrix& mx);

// ____________________________________________________________________________
// D                  OPERATOR FUNCTIONS, OPERATOR WITH SCALAR
// ____________________________________________________________________________

/// Center Operator - Scalar Functions 
// F_list * 		     - Operator-Scalar Multiplication
// F_list *= 		     - Operator-Scalar Unary Multiplication
// F_list / 		     - Operator-Scalar Division
// F_list *= 		     - Operator-Scalar Unary Division

MSVCDLL friend gen_op operator * (const gen_op &Op1, const complex &z);
MSVCDLL friend gen_op operator * (const complex &z,  const gen_op &Op1);
MSVCDLL friend gen_op operator * (const gen_op &Op1, double z);
MSVCDLL friend gen_op operator * (double z,          const gen_op &Op1);

	// Input		Op1  : General operator
        //                      z    : Complex number
	// Return	        Op   : General operator which is
	//			       the input operator multiplied by
        //	                       the complex number, z*Op1
	// Note			     : Result Op EXCLUSIVELY in WBR


MSVCDLL gen_op & operator *= (const complex &z);
MSVCDLL gen_op & operator *= (double r);

MSVCDLL friend gen_op operator / (const gen_op &Op1, const complex &z);


	// Input		Op1  : General operator
        //			z    : Complex number
	// Return	        Op   : General operator which is the
	//			       input operator Op1 multiplied
	//			       by the inverted complex number.
	//			       Op = (1/z) * Op1
	// Note			     : Result EXCLUSIVELY in WBR of Op1


MSVCDLL friend gen_op operator / (const gen_op& Op1, double r);

	// Input		Op1  : General operator
        //			z    : A real number
	// Return	        Op   : General operator which is the
	//			       input operator Op1 multiplied
	//			       by the inverted real number.
	//			       Op = (1/r) * Op1
	// Note			     : Result EXCLUSIVELY in WBR of Op1
     

MSVCDLL gen_op & operator /= (const complex& z);
MSVCDLL gen_op & operator /= (double r);

// ____________________________________________________________________________
// E                       COMPLEX OPERATOR FUNCTIONS
// ____________________________________________________________________________

/// Center Complex Operator Functions 
/// F_list det		     - Operator Determinant
/// F_list trace	     - Operator Trace
/// F_list proj		     - Operator Projection
/// F_list dim		     - Operator dimension
/// F_list exp		     - Operator exponential 


MSVCDLL complex det() const;

	// Input		Op   : General operator (this)
        // Return		z    : Complex number, determinant of Op
	// Note			     : Performed in WBR of Op


MSVCDLL complex trace() const;

	// Input		Op   : General operator (this)
        // Return		z    : Complex number, trace of Op
	// Note			     : Performed in WBR of Op


MSVCDLL complex trace(const gen_op& Op2) const;

	// Input		Op1  : General operator (this)
	// 			Op2  : General operator
        // Return		z    : Complex, trace of (Op1*Op2)
	// Note			     : Performed in WBR of Op1


MSVCDLL complex proj(const gen_op &Op2, int norm=0) const;

	// Input		Op1  : General operator (this)
	// 			Op2  : General operator
	// 			norm : Value of normalization
	//			       if zero it is computed
        // Return		z    : Complex, projection of Op1 on Op2
	//				      t		    t
	//				tr(Op2 * Op1)/tr(Op2  * Op2)
	// Note			     : Performed in WBR of Op2


MSVCDLL int dim()    const;		// Hilbert space
MSVCDLL int HS()     const;		// Hilbert space
MSVCDLL int LS()     const;		// Liouville space
MSVCDLL int dim_LS() const;		// Liouville space

	// Input		Op   : General operator (this)
        // Return		int  : Op Liouville space dimension
  
MSVCDLL gen_op exp() const;

	// Input		Op1  : General operator (this)
        // Return		Op   : exponential of Op1
	//			       Op = exp(Op1)
        // Note			     : Computed in EBR of Op1


MSVCDLL gen_op exp(const complex& t, double cutoff=1.e-12) const;

        // Input                Op    : Operator (this)
        //                      t     : Exponential factor
        //                      cutoff: Exponential factor roundoff
        // Return               ExpOp : Exponential of Op
        //                              ExpOp = exp(t*Op)
        // Note                       : Exponential output in EBR of Op
        // Note                       : Op's EBR is generated herein
        // Note                       : Value of t is considered 0 if
        //                              it's magnituded is less than cutoff


MSVCDLL gen_op Pow(int power) const;

        // Input                Op      : Operator (this)
        //                      power   : Exponential power
        // Return               Op^n    : Op taken to the nth power
        // Note                         : Output in EBR of Op
        // Note                         : Op's EBR is generated herein


MSVCDLL friend gen_op tensor_product(const gen_op &Op1, const gen_op &Op2);

        // Input                Op1  : General operator.
        //                      Op2  : General operator.
        // Return               Op   : Operator tensor product of the two input
        //                             operators, Op =  Op1 X Op2.
        // Note                      : Order matters - Op1XOp2 != Op2XOp1

 
MSVCDLL friend gen_op log(const gen_op &Op1);
 
        // Input                Op1  : General operator
        // Return               Op   : logarithms of Op1
        //                             Op = log(Op1)
        // Note                      : Computed in EBR of Op1


MSVCDLL gen_op sim_trans(const gen_op &Op2) const;

	// Input		Op1  : General operator (this)
	// 			Op2  : General operator
        // Return		Op   : General operator from similarity
	//			       transform of Op2 by Op1 
	// 			       Op = Op1 * Op2 * [Op1]
        // Note			     : Op EXCLUSIVELY in WBR of Op1
	// F_list sim_trans	     - Operator similarity transformation

  
MSVCDLL void sim_trans_ip(const gen_op &Op1);

	// Input		Op   : General operator (this)
	// 			Op1  : General operator
        // Return		None : Similarity transform of Op by Op1 
	//					            t
	// 			       Op = Op1 * Op * [Op1]
        // Note			     : Op EXCLUSIVELY in WBR of Op1
	// F_list sim_trans_ip	     - Operator similarity transformation

  
//  gen_op adjoint();

	// Input		Op1  : General operator(this)
        // Return		Op   : Hermitian adjoint of Op1
	//                          	          *t
	// 			       Op = [Op1]  (transpose & conjugate)
        // Note                      : EXCLUSIVELY in the WBR of Op
	// F_list adjoint	     - Operator Hermitian adjoint


MSVCDLL row_vector eigvals() const;
MSVCDLL void       eigvals(double* vx, int sort=0) const;

	// Input		Op   : General operator (this)
        // 			vx   : Vector of doubles 
	//			sort : Flag for eigenvalue sorting 
        // Return		none : Vector w filled with eigenvalues
	// F_list eigvals	     - Operator Eigenvalues


/***************************************************************************
****************************************************************************
  Now marked for deletion, these have been changed to member functions
  Must wait until matrix functions are also member - Sosi 2/21/92         */

MSVCDLL friend complex det(const gen_op &Op);					//**
MSVCDLL friend complex trace(const gen_op &Op);					//**
MSVCDLL friend complex trace(const gen_op &Op1, const gen_op &Op2);		//**
MSVCDLL friend complex proj(const gen_op, const gen_op, int);  //**
MSVCDLL friend int dim(const gen_op &Op);  					//**
MSVCDLL friend gen_op exp(const gen_op &Op1);					//**
MSVCDLL friend gen_op pow(const gen_op &Op1, int power);			//**
MSVCDLL friend gen_op sim_trans(const gen_op &Op1, const gen_op &Op2);		//**
MSVCDLL friend gen_op adjoint(const gen_op &Op1);				//**

/***************************************************************************
***************************************************************************/

// ____________________________________________________________________________
// F                      OPERATOR COMPONENT MANIPULATIONS
// ____________________________________________________________________________

/// Center Internal Access
//  F_list get_mx	      - Retrieve Operator Matrix
//  F_list put_mx	      - Assign Operator Matrix
//  F_list get_basis	      - Retrieve Operator Basis
//  F_list put_basis	      - Assign Operator Basis
//  F_list ()		      - Retrieve Operator Element
//  F_list put		      - Assign Operator Element
//  F_list get		      - Retrieve Operator Element

/*     Function    Arguments                           Result
   --------------- ---------   ------------------------------------------------
       get_mx         ---      Returns the matrix of current OpRep (WBR)
     get_matrix       ---      Returns the matrix of current OpRep (WBR)
       put_mx          mx      Sets WBR matrix to be mx, other OpReps deleted
     put_matrix        mx      Sets WBR matrix to be mx, other OpReps deleted
       get_bs         ---      Returns the basis of current OpRep (WBR)
     get_basis        ---      Returns the basis of current OpRep (WBR)
       put_bs          bs      Sets WBR basis to be bs, other OpReps deleted
     put_basis         bs      Sets WBR basis to be bs, other OpReps deleted
       (i,j)         int,int   Returns <i|Op|j> in current OpRep (WBR) 
        get          int, int  Returns <i|Op|j> in current OpRep (WBR) 
        put          int, int  Sets <i|Op|j> in WBR, others deleted

  Note that Op(i,j) returns a copy of the operator element, NOT the element.
  Thus, code such as "Op(i,j) = z;" will not work!                           */

MSVCDLL matrix get_mx( )     const;
MSVCDLL matrix get_matrix( ) const;
MSVCDLL void   put_mx(const matrix &mx);
MSVCDLL void   put_matrix(const matrix &mx);

MSVCDLL basis get_bs( )    const;
MSVCDLL basis get_basis( ) const;
MSVCDLL void  put_bs(const basis &bs);
MSVCDLL void  put_basis(const basis &bs);

MSVCDLL complex operator() (int row, int col) const;
MSVCDLL complex get(int row, int col)         const;
MSVCDLL void    put(const complex &z, int row, int col);

// ____________________________________________________________________________
// G                      OPERATOR AUXILIARY FUNCTIONS
// ____________________________________________________________________________

// F_list exists	      - Test for Operator Existence
// F_list name		      - Get/Set Operator Name
// F_list bsname	      - Get Operator Basis Name (WBR)

MSVCDLL std::string name() const;
MSVCDLL void name(const   std::string& n);
MSVCDLL void bsname(const std::string& bn);

MSVCDLL int exists( ) const;

	// Input		Op    : General operator (this)
	// Output		int   : TRUE if Op not NULL


MSVCDLL col_vector superket() const;

        // Input                Op    : General operator (this)
        // Output               vec   : Vector of alligned operator elements
	// Note			      : This function supports MultiSys


MSVCDLL void desuperket(const col_vector& mS);

        // Input                vec   : Superket iniitally associated with Op
        // Output               Op    : General operator
	// Note			      : This function supports MultiSys


MSVCDLL gen_op project_sub(int ic) const;

        // Input                Op    : General operator (this)
        //                      int   : Multi_sys component ic
        // Output                     : Reduced general operator for comp ic
	// Note			      : This function supports MultiSys

// ____________________________________________________________________________
// H                  OPERATOR REPRESENTATION MANIPULATIONS
// ____________________________________________________________________________

// F_list set_DBR	      - Set Operator to Default Basis
// F_list set_EBR	      - Set Operator to EigenBasis
// F_list set_DBR	      - Set Operator to Default Basis
// F_list setOnlyWBR	      - Zero All Operator Reps. but Current
// F_list test_EBR	      - Test for Operator Eigenbasis
// F_list status	      - Output current Operator status
// F_list Op_priority	      - Set Priority for Current Operator Rep.
// F_list set_limits	      - Insure Operator Reps. Under Limit

/* These functions are quick checks to see of a particular reperesentation
   exists and/or whether the current represenation is the working OpRep.  The
   functions check for EBR & DBR, an eigenbasis and the default basis OpReps.*/

MSVCDLL int test_EBR() const;
MSVCDLL int test_DBR() const;
MSVCDLL int in_EBR()   const;
MSVCDLL int in_DBR()   const;

/* These functions are will put an operator into a particular basis.  A new
   OpRep will be formed if the representation doesn't already exist.  The
   Op_base functions set "this" into the basis/Op.WBR input as an argument. 
   A cutoff can be specified that will "check" the validity of an EBR.       */
	
MSVCDLL void set_DBR() const;
MSVCDLL void set_EBR() const;
MSVCDLL void Op_base(const gen_op &Op1, double cutoff=1.e-12) const; 
MSVCDLL void Op_base(const basis  &bs) const;

MSVCDLL void status(int pf=0) const;

	// Input		Op    : General operator (this)
	// Output		void  : Outputs operator status


MSVCDLL void setOnlyWBR( );

	// Input		Op   : General operator (this)
	// Output		none : Deletes all reps but WBR of Op
	// Note			     : Does nothing if Op is NULL
	// Note			     : At finish nreps = 1 (temporarily)


// ******************* Representation Limits and Priorities *******************

MSVCDLL void Op_priority (int pty);

	// Input		Op   : General operator (this).
	// 			pty  : Priority value
	// Output		none : Assigns Op WBR priority
	// Note			     : Higher pty value implies rep
	//			       preferentially maintained

MSVCDLL void SetLimits(int limit) const;

	// Input		Op   : General operator (this).
	// 			limit: Representation limit number 
	// Output		none : Insures number of Op reps
	//			       does not exceed limit
// sosi *** How is this different from limit reps? 


// ____________________________________________________________________________
// I                  CLASS OPERATOR EQUALITY & INEQUALITY
// ____________________________________________________________________________

/* These functions check whether two operators are equivalent.  They are taken
   to be equal if their default basis matrix representations are the same. If
   class gen_op is protecting its data properly, there is no way to generate
   two operators having the same DBR and not being equal.  Not that this then
   assumes that the number of representations may be unequal even though the
   operators are taken as equal.

           Input                Op      : First operator (this) 
                                Op1     : Second operator
           Output               TF      : True if Op = Op1 (==)
                                          True if Op != Op2 (!=)
           sosi                         : Add == to genoprep someday?      */ 

// sosi - these are now in section M
//int operator== (const gen_op& Op1) const;
//int operator!= (const gen_op& Op1) const;

// ____________________________________________________________________________
// J                        CLASS OPERATOR CHECKS
// ____________________________________________________________________________

/* This set of functions insures that the Op/mx/bs is properly square and/or   
   of compatible dimension to mix with another Op/mx/bs.  All return TRUE if
   dimensions look OK.  The functions take a warning level flag which will
   be used to decide what happens under a FAIL.  If the flag is 0 then nothing
   will result.  If the flags is 1 a non-fatal error is output.  If the flag
   is >1 a fatal error is output and program execution stopped.              */
  
MSVCDLL int  OpCheck(const gen_op& Op1,                    int warn=2) const;
MSVCDLL bool OpCheck(const matrix& mx,                     int warn=2) const;
MSVCDLL int  OpCheck(const basis&  bs,                     int warn=2) const;
MSVCDLL int  OpCheck(const matrix& mx1, const matrix& mx2, int warn=2) const;
MSVCDLL int  OpCheck(const matrix& mx,  const basis&  bs,  int warn=2) const;

/* This set of functions insures that the Op limits are met, Op element
   access is within matrix representation bounds, etc.                       */

MSVCDLL int OpCheck(int row,    int col, int warn=2) const;
MSVCDLL int LimCheck(int limit,          int warn=1) const;


// ____________________________________________________________________________
// K                     CLASS OPERATOR I/O FUNCTIONS
// ____________________________________________________________________________

// ------------------------ ASCII Output Functions ----------------------------

/* 		Input		Op   : Operator (this)
				ostr : Output ASCII file stream
                                full : Flag for amount of output
		Return		void : Op is sent to the output stream 
		Note		     : Op WBR representation only	     */

MSVCDLL        std::ostream& print(std::ostream& out, int full=0) const;
MSVCDLL friend std::ostream& operator << (std::ostream& ostr, const gen_op &Op);
MSVCDLL std::ostream& eigenvalues(std::ostream& ostr, int rc=0, int ncol=4) const;

// ------------------------ Binary Output Functions ---------------------------

/* 		Input		Op   : Operator (this)
				fn   : Input binary filename
				fp   : File stream (pointing at Op)
		Return		void : Op is written to either the
				       specified file or filestream.
		Note		     : Only WBR representation is output     */

MSVCDLL void write(const std::string& fn) const;
MSVCDLL std::ofstream& write(std::ofstream& fp) const;

// ------------------------ Binary Input Functions ----------------------------

/* 		Input		Op   : Operator (this)
				fn   : Input binary filename
				fp   : File stream (pointing at Op)
				bs   : An optional basis in which to
				       force the operator into so that many
				       read operators may share one basis!
		Return		void : Op is read in from either the
				       specified file or filestream.
		Note		     : Read Op will reside in only 1 basis   */

//void read(const std::string& fn, gen_op& Op1);		NOW DEPRECATED!
MSVCDLL void read(const std::string& fn, const basis& bs);
MSVCDLL void read(const std::string& fn);
MSVCDLL std::ifstream& read(std::ifstream &fp);

// ____________________________________________________________________________
// L                     CLASS OPERATOR TESTING FUNCTIONS
// ____________________________________________________________________________
 
/* Since GAMMA operators take the brunt of the work in most MR simulations, it
   is essential that they are in proper working order. As such, I have included
   some means to test various essential aspects of them. In particular, two
   basic aspects are the ability to handle their diagonalization and multiple
   representation (OpRep) tracking. Some of these tests simply parallel the
   functionality provided in class matrix.

     Function    Arguments                           Result
 --------------- --------- ------------------------------------------------
 TestEigenSystem    int    Diagonalize the operator & check its eigensystem  */

MSVCDLL double TestEigenSystem(int pf=0);

// ----------------------------------------------------------------------------
//                     Matrix Functionality For GenOp
// ----------------------------------------------------------------------------

/* This allows GenOp to act like matrix when using these testing functions.

     Function    Output                       Description
   ------------  ------  ------------------------------------------------------
   is_symmetric   bool    TF if <i|mx|j>-<j|mx|i>  < d       (def. d GMxCut)
   is_hermitian   bool    TF if <i|mx|j>-<j|mx|i>* < d       (def. d GMxCut)
   is_unitary     bool    TF if inv(mx) == adjoint mx, CPU intensive
   is_real        bool    TF if Im(<i|mx|j>) < d for all i,j (def. d GMxCut)
   is_imaginary   bool    TF if Re(<i|mx|j>) < d for all i,j (def. d GMxCut)
   is_complex     bool    TF if is_real && is_imaginary      (def. d GMxCut)
   is_zero        bool    TF if ||<i|mx|j>|| < d for all i,j (def. d GMxCut)
   is_diagonal    bool    TF if ||<i|mx|j>|| < d for all i!=j(TRUE,d GMxCut)
   is_square      bool    TF if rows_ == cols_                               */

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
// M                 CLASS OPERATOR CONTAINER SUPPORT FUNCTIONS
// ____________________________________________________________________________
 
/* Aside for providing basic tests as to whether to operators are equvalent or
   not, these operators are necessary if any STL container classes are to
   be used based on operators (e.g. list<gen_op> or vector<gen_op>)          */  

MSVCDLL bool operator== (const gen_op& Op) const;
MSVCDLL bool operator!= (const gen_op& Op) const;
MSVCDLL bool operator<  (const gen_op& Op) const;
MSVCDLL bool operator>  (const gen_op& Op) const;
 
};


// ____________________________________________________________________________
//                   CLASS OPERATOR INLINE UTILITY FUNCTIONS
// ____________________________________________________________________________

inline gen_op I_gen_op(const basis& bs)

	// Input		bs   : Basis
	// Return		Op   : Identity Operator with basis bs 

  { return gen_op(matrix(bs.size(),bs.size(),i_matrix_type),bs); }
 
#endif 							// GenOp.h
