/* FloqOp.h *****************************************************-*-c++-*-
**                                                                      **
**                           G A M M A                                  **
**                                                                      **
**      Floquet Operators                           Interface		**
**                                                                      **
**      Copyright (c) 1992                                              **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      Scott A. Smith                                                  **
**      Copyright (c) 1994, 1997                                        **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** Class floq_op embodies Floquet operators in GAMMA.  The class sets	**
** defines the properties and allowed operations such operators such	**
** that they may be freely used in GAMMA based programs.		**
**                                                                      **
*************************************************************************/

#ifndef   Gfloq_op_h_			// Is file already included?
#  define Gfloq_op_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// is in the GNU compiler
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/matrix.h>		// Include knowledge of arrays
#include <HSLib/Basis.h>		// Include knowledge of bases
#include <HSLib/GenOp.h>		// Include knowledge of operators 

class floq_op : public gen_op

{

// ----------------------------------------------------------------------------
// ------------------------------- STRUCTURE ----------------------------------
// ----------------------------------------------------------------------------

private:
                                // can be set by the user directly
  double Om;    		// Frequency of Fourier Expansion
  int N ;			// 2*N+1: Dimension of `Photon-space'
  int hs;			// Dimension of Hilbert space

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
      
// ____________________________________________________________________________
// i                 CLASS FLOQUET OPERATOR ERROR HANDLING
// ____________________________________________________________________________

        // Input                FOp     : Floquet operator (this)
        //                      eidx    : Error index
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error message
        //                                Program execution stopped (fatal)

         void FOperror(int eidx, int noret=0) const;
volatile void FOpfatality(int eidx) const;

// ____________________________________________________________________________
// ii            FLOQUET OPERATOR COMPATIBILITY CHECKING
// ____________________________________________________________________________
 
        // Input                FOp     : Floquet operator (this)
        //                      FOp1    : A second floquet operator
	//			Fmx	: An array in Floquet space
	//			N1,N2   : Photon indices
	//			H1,H2	: Operator indices (Hilbert space)
        //                      warn    : Warning level
        //                                    0 = No warnings (default)
        //                                    1 = Non-fatal warnings
        //                                   >1 = Fatal warnings
        // Output (FOP1, warn)  TF      : True if the two operators match in
        //                                three areas -
        //                                     1. Phonon space sizes (N)
        //                                     2. Hilbert space sizes (hs)
        //                                     3. Fourier frequencies (Om)
        // Output (Fmx, warn)   TF      : True if the operator & matrix match
        //                                Floquet space dimensions
        // Output (N1,N2,H1,H2) TF      : True if the element <I|FOp|J> is
        //                                within Floquet space dimensions where
	//				  I=(N+N1)*hs+H1 & J=(N+N2)*hs+H2

int FOpCheck(const floq_op& FOp1,            int warn=0) const;
int FOpCheck(const matrix& Fmx,              int warn=0) const;
int FOpCheck(int N1, int N2, int H1, int H2, int warn=0) const;

// ____________________________________________________________________________
// iii            FLOQUET OPERATOR - HILBERT OPERATOR BASES
// ____________________________________________________________________________

void SetBasis(const gen_op& Op);

        // Input		FOp	: Floquet operator (this)
	//			Op	: A Hilbert space opertor
	// Output		void	: FOp has its basis set to the
	//				  tensor product E X OpBs
	// Note				: PRIVATE - no FOp/Op compatibility
	//				  checks are made herein


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A             CLASS FLOQUET OPERATOR CONSTRUCTORS/DESTRUCTOR
// ____________________________________________________________________________

public:
  
MSVCDLC floq_op( );
MSVCDLC floq_op(const floq_op& F);

	// Input                N_   : Floquet space dimension (2*N+1)
        //                      hs_  : Hilbert space dimension (2*hs+1)
        //                   omega_  : Rotation frequency
        // 	                mx   : Matrix
        //                      bs   : Basis
        // Output		Op   : Floquet operator (this) constructed
        //                             2*N+1*hs and k*wr(k=-N...N) on main
        //                             diagonal.
	//			       mx, no bs: with matrix mx as DBR
	//			       mx, bs:    with matrix mx & basis bs 

MSVCDLC floq_op(int N_, int hs_, double omega_);
MSVCDLC floq_op(int N_, int hs_, double omega_, const matrix& mx);
MSVCDLC floq_op(int N_, int hs_, double omega_, const matrix& mx, const basis& bs);

	// Input		FOp1 : Floquet operator.
	// 			FOp  : Floquet operator (this).
	// Return		FOp  : Operator which is a copy of
	//			       the input operator, FOp = FOp1.
	// Note		             : Result EXCLUSIVELY in WBR of FOp1
  
MSVCDLL floq_op& operator = (const floq_op& FOp1);
MSVCDLC          ~floq_op();

// ____________________________________________________________________________
// B                   FLOQUET OPERATOR ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions are used to access (but not set) internal dimensions and
   parameters regarding a Floquet opertor.  The following functions are set:

   hsdim  - returns base Hilbert space dimension
   phodim - returns photon space dimension
   phodim - returns Fourier expansion frequency
   dim    - returns full dimension of FloqOp (2N+1)*hs
   size   - returns full dimension of FloqOp (2N+1)*hs                       */

MSVCDLL int    hsdim()  const;
MSVCDLL int    phodim() const;
MSVCDLL double omega()  const;
MSVCDLL int    dim()    const;
MSVCDLL int    size()   const;

// ____________________________________________________________________________
// C   FLOQUET OPERATOR FUNCTIONS, FLOQUET OPERATOR WITH FLOQUET OPERATOR 
// ____________________________________________________________________________

/* These functions handle all dealings between Floquet operators and Floquet
   operators.

 Function   Arguments      Return       Function   Arguments      Return
 --------   ---------   -------------   --------   ---------   -------------
    +       FOp1,FOp2   FOp=FOp1+FOp2      *       FOp1,FOp2   FOp=FOp1*FOp2
    +=      FOp1        FOp=FOp+FOp1       *=      FOp1        FOp=FOp*FOp1
    -       FOp1,FOp2   FOp=FOp1-FOp2      &=      FOp1        FOp=FOp1*FOp 
    -=      FOp1        FOp=FOp-FOp1
    -       FOp1        FOp=-FOp1                                            */

MSVCDLL        floq_op operator +  (const floq_op &FOp2) const;
MSVCDLL        void    operator += (const floq_op &FOp1);
MSVCDLL        floq_op operator -  (const floq_op &FOp2) const;
MSVCDLL        floq_op operator -  ()                    const;
MSVCDLL        void    operator -= (const floq_op &FOp1);
MSVCDLL        floq_op operator *  (const floq_op &FOp2) const;
MSVCDLL        void    operator *= (const floq_op &FOp1);
MSVCDLL        void    operator &= (const floq_op &FOp1);

// ____________________________________________________________________________
// D        FLOQUET OPERATOR FUNCTIONS, FLOQUET OPERATOR WITH MATRIX
// ____________________________________________________________________________

/* These functions handle all dealings between Floquet operators and Floquet
   dimension matrices.

 Function   Arguments      Return       Function   Arguments      Return
 --------   ---------   -------------   --------   ---------   -------------
    +       FOp1,Fmx    FOp=FOp1+Fmx       *       FOp1,FMx    FOp=FOp1*FMx
    +       Fmx,FOp1    FOp=FOp1+Fmx       *       FMx,FOp1    FOp=FMx*FOp1
    +=      FOp1        FOp=FOp+FOp1       *=      FMx         FOp=FOp*FMx
    -       FOp1,Fmx    FOp=FOp1-Fmx       &=      FMx         FOp=FMx*FOp 
    -=      Fmx         FOp=FOp-Fmx        -       Fmx,FOp1    FOp=Fmx-FOp1  */

MSVCDLL        floq_op operator +  (const matrix&  Fmx)  const;
MSVCDLL friend floq_op          operator +  (const matrix&  Fmx,  const floq_op& FOp1);
MSVCDLL        void    operator += (const matrix&  Fmx);
MSVCDLL        floq_op operator -  (const matrix&  Fmx)  const;
MSVCDLL friend floq_op          operator -  (const matrix&  Fmx,  const floq_op& FOp1);
MSVCDLL        void    operator -= (const matrix&  Fmx);
MSVCDLL        floq_op operator *  (const matrix&  Fmx)  const;
MSVCDLL friend floq_op          operator *  (const matrix&  Fmx,  const floq_op& FOp1);
MSVCDLL        void    operator *= (const matrix&  Fmx);
MSVCDLL        void    operator &= (const matrix&  Fmx);

// ____________________________________________________________________________
// E     FLOQUET  OPERATOR FUNCTIONS, FLOQUET OPERATOR WITH SCALAR
// ____________________________________________________________________________

/* These functions handle all dealings between Floquet operators and scalar
   values (int, double, complex).

 Function   Arguments      Return       Function   Arguments      Return
 --------   ---------   ------------    --------   ---------  ----------------
    *        FOp1,z     FOp = z*FOp1       *=       this, z   void, z*this
    *        z,FOp1     FOp = z*FOp1       *=       this, d   void, d*this
    *        FOp1,d     FOp = d*FOp1        /       FOp1, z   FOp = (1/d)*this
    *        d,FOp1     FOp = d*FOp1       /=       this, z   void, (1/d)*this
                                                                             */
MSVCDLL        floq_op operator *  (const complex &z) const;
MSVCDLL friend floq_op          operator *  (const complex& z, const floq_op &FOp1);
MSVCDLL        floq_op operator *  (double d)         const;
MSVCDLL friend floq_op          operator *  (double d,         const floq_op &Op1);
MSVCDLL        void    operator *= (const complex& z);
MSVCDLL        void    operator *= (double d);
MSVCDLL        floq_op operator /  (const complex &z) const;
MSVCDLL        floq_op operator /  (double d)         const;
MSVCDLL        void    operator /= (const complex& z);
MSVCDLL        void    operator /= (double d);

// ____________________________________________________________________________
// F                  COMPLEX  FLOQUET OPERATOR FUNCTIONS
// ____________________________________________________________________________


MSVCDLL friend gen_op pho_trace(const floq_op& Op);

	// Input		Op   : Floquet operator
        // Output          gen_op    : Operator containing trace(s) over    
        //                             Photon space

MSVCDLL friend floq_op exp(const floq_op& Op1);
	
        // Input		Op1  : Floquet operator
        // Return		Op   : exponential of Op1
	//			       Op = exp(Op1)
        // Note			     : Computed in EBR of Op1

MSVCDLL floq_op exp();

        // Input                FlOp    : Floquet operator (*this)
        // Return               ExpFlOp : Exponential of Op1
        //                                   Op = exp(Op1)

MSVCDLL friend floq_op prop(floq_op& FLOQHAM, double time);
 
	// Input		Op1   : Floquet Hamiltonian
        //                      time  : Evolution time
        // Return		FP    : Floquet Propagator in DBR
	

MSVCDLL friend floq_op fprop(floq_op &FLOQHAM, double time);
  
	// Input		Op1   : Floquet Hamiltonian
        //                      time  : Evolution time
        // Return		FP    : Floquet Propagator in DBR
	
// ___________________________________________________________________________
// G               FLOQUET OPERATOR COMPONENT MANIPULATIONS
// ___________________________________________________________________________

// --------------------- Floquet Submatrix Manipulations ---------------------

MSVCDLL gen_op operator() (int N1, int N2) const;
MSVCDLL gen_op get_block  (int N1, int N2) const;
     
        // Input                FOp     : Floquet operator (this)
        //                      N1      : Row index of photon dimension
        //                      N2      : Column index of photon dimension
        // Output               Op      : Operator at <N1|FOp|N2>


MSVCDLL void put_block(const gen_op &Op1, int N1, int N2);

        // Input                FOp     : Floquet operator (this)
        //                      N1      : Row index of photon dimension
        //                      N2      : Column index of photon dimension
	//			Op	: General operator
        // Output               none	: Op will be put at <N1|FOp|N2> 



MSVCDLL void put_sdiag(const gen_op &Op, int sdn_);
        
	// Input                FOp     : Floquet operator (this)
	//                       Op 	: General operator
	//                      sdn_	: Number of side diagonals to be set
	// Output                FOp1	: Floquet operator with side diagonals
	//				: set to Op.  The number of side
	//				  diagonals set will be sdn_

  
// -------------------- Individual Element Manipulations ----------------------


MSVCDLL void put(const complex& z, int row, int col);
       
       // Input                      : Floq_operator (this)
       //                 row,col    : Indices 
       // Output            none     : Value of <row|Op|col> in WBR
       //                            : set to z
       // Note                       : floq_op will be changed
 

MSVCDLL complex get(int row, int col) const;

       // Input                      : Floq_operator (this)
       //                 row,col    : Indices 
       // Output            none     : Value of <row|Op|col> in WBR
       // Note                       : floq_op remains unchanged


MSVCDLL void put(const complex& z, int N1_, int N2_, int H1_ , int H2_ );
       
       // Input                      : Floq_operator (this)
       //                 N1_,N2_    : Indices of photon space:-N...N
       //                 H1_,H2_    : Indices of Hilbert space:0...hs
       // Output            none     : Value of <N2,H2|Op|H1,N1> will be
       //                            : set to z
       // Note                       : floq_op will be changed
 

MSVCDLL complex get(int N1_ , int N2_ , int H1_ , int H2_) const;


       // Input                      : Floq_operator (this)
       //                 N1_,N2_    : Indices of photon space:-N...N
       //                 H1_,H2_    : Indices of Hilbert space:0...hs
       // Output                z    : complex value of <N2_,H2_|Op|H1_,N1_>
       // Note                       : floq_op remains unchanged

// ____________________________________________________________________________
// H                OPERATOR REPRESENTATION MANIPULATIONS
// ____________________________________________________________________________

	// Input		Op    : Floquet operator (this)
	// Output (set_DBR)	none  : Op WBR set to DBR 
	// Output (set_EBR)	none : Op WBR set to EBR 
 
MSVCDLL void set_DBR() const;
MSVCDLL void set_EBR() const;

// ____________________________________________________________________________
// I                 CLASS FLOQUET OPERATOR I/O FUNCTION
// ____________________________________________________________________________

	// Input		FOp	: Floquet operator (this)
	//			ostr    : An output stream
	// Output		void	: Output stream altered by FOP 
	//				  details being placed into it

MSVCDLL void print(std::ostream& ostr, int full=1, std::string spacer="\t") const;
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const floq_op& FOp);

 // ___________________________________________________________________________ 
 // J             INTERNAL FLOQUET OPERATOR MANIPULATIONS
 // ___________________________________________________________________________       

        // Input                     : Floquet Operator (this)
        // Output add_omega  F_op    : Floquet operator with omegas added
        //                           : on main diagonal
        // Output sub_omega  F_op    : Floquet operator with omegas subtracted
        //                           : on main diagonal   

MSVCDLL void add_omega();           
MSVCDLL void sub_omega();

};

#endif						// FloqOp.h
