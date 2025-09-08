/* Floq2Op.h ********************************************-*-c++-*-
**                                                              **
**                            G A M M A                         **
**                                                              **
**      Floquet Operator                        Interface	**
**                                                              **
**      Copyright (c) 1992                                      **
**	Marc Baldus						**
**      Eidgenoessische Technische Hochschule                   **
**      Labor fuer physikalische Chemie                         **
**      8092 Zurich / Switzerland                               **
**                                                              **
**      $Header: $
**								**
*****************************************************************/

/*****************************************************************
**								**
** 	Description				    **
**						    **
**	The class floquet operator defines	    **
**	an operator for C++ with all the	    **
**	necessary operations (+, -, *, /)	    **
**	functions (exp) and input/output	    **
**	routines				    **
**						    **
*****************************************************/

#ifndef   GFloq2Op_h_			// Is file already included
#  define GFloq2Op_h_ 1			// no, then remember it
#  if defined(GAMPRAGMA)		// is it the GNU compiler
#    pragma interface			// this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/matrix.h>
#include <HSLib/Basis.h>
#include <HSLib/GenOp.h>

class floq2_op : /* private  */ public gen_op

{
    
// ----------------------------------------------------------------------
// ---------------------------- STRUCTURE -------------------------------
// ----------------------------------------------------------------------

private:
                                // can be set by the user directly
  double _omega1,_omega2;    	// Frequency of Fourier Expansion
  int N1,N2 ;			// 2*N+1: Dimension of `Photon-space'
  int hs;			// Dimension of Hilbert space



// ----------------------------------------------------------------------
// ------------------------ PRIVATE FUNCTIONS ---------------------------
// ----------------------------------------------------------------------
      
// ______________________________________________________________________
//                 CLASS FLOQUET OPERATOR ERROR HANDLING
// ______________________________________________________________________

        // Input                FOp     : Floquet operator (this)
        //                      eidx    : Error index
        //                      noret   : Flag for return (0=linefeed)
        // Output               none    : Error message
        //                                Program execution stopped (fatal)

//         void FOperror(int eidx, int noret=0) const;
//volatile void FOpfatal(int eidx) const;

friend void floq2_op_error (int error);
friend void volatile floq2_op_fatality (int error);

// ----------------------------------------------------------------------
// ------------------------- PUBLIC FUNCTIONS ---------------------------
// ----------------------------------------------------------------------

public:

// ______________________________________________________________________
// A             CLASS FLOQUET OPERATOR CONSTRUCTORS/DETRUCTOR
// ______________________________________________________________________

MSVCDLC floq2_op();
MSVCDLC floq2_op(const floq2_op& FOp);

        // Input                N_   : 2*N+1: Dimension of Floquet space
        //                      hs_  : 2*hs+1: Dimension of spin space
        //                   omega_  : Rotation frequency
        // Output                    : returns a new floq_op with Dimension
        //                             2*N+1*hs and k*wr(k=-N...N) on main
        //                             diagonal.

	// Input                N_   : Floquet space dimension
        //                      hs_  : Hilbert space dimension
        //                   omega_  : Rotation frequency
        // 	                mx   : Matrix
        // Output		Op   : Floquet operator (this) constructed
	//			       with matrix mx as DBR
        // Note			     : Assumes mx DBR of Op
 
MSVCDLC floq2_op (int N1_,int N2_,int hs_,double omega1_,double omega2_)
   :gen_op(matrix((2*N1_+1)*(2*N2_+1)*hs_,(2*N1_+1)*(2*N2_+1)*hs_,0,n_matrix_type))
      , _omega1(omega1_),_omega2(omega2_), N1(N1_),N2(N2_), hs(hs_) {};

	// Input                N    : First Floquet space dimension
	// Input                M    : Second Floquet space dimension
        //                      HS   : Hilbert space dimension
        //                     Om1   : First rotation frequency
        //                     Om2   : Second rotation frequency
        // 	                mx   : Matrix
        //                      bs   : Basis
        // Output		Op   : Floquet operator (this) constructed
	//			       with matrix mx and basis bs
        // Note			     : Assumes mx DBR of Op
 
MSVCDLC floq2_op(int N,int M,int HS,double Om1,double Om2,matrix& mx);
MSVCDLC floq2_op(int N,int M,int HS,double Om1,double Om2,matrix& mx, basis &bs);

MSVCDLC floq2_op& operator = (const floq2_op& FOp1);
MSVCDLC           ~floq2_op ();

// ____________________________________________________________________________
// B                   FLOQUET OPERATOR ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions are used to access (but not set) internal dimensions and
   parameters regarding a Floquet opertor.  The following functions are set:

   hsdim   - returns base Hilbert space dimension
   phodim1 - returns first photon space dimension
   phodim2 - returns second photon space dimension
   omega1  - returns first Fourier expansion frequency
   omega2  - returns second Fourier expansion frequency
   dim     - returns full dimension of FloqOp
   size    - returns full dimension of FloqOp	                      */

MSVCDLL int    dim()     const;
MSVCDLL int    size()    const;
MSVCDLL int    hsdim()   const;
MSVCDLL int    phodim1() const;
MSVCDLL int    phodim2() const;
MSVCDLL double omega1()  const;
MSVCDLL double omega2()  const;

// ______________________________________________________________________
// C FLOQUET OPERATOR FUNCTIONS, FLOQUET OPERATOR WITH FLOQUET OPERATOR 
// ______________________________________________________________________

/* These functions handle all dealings between Floquet operators and Floquet
   operators.

 Function   Arguments      Return       Function   Arguments      Return
 --------   ---------   -------------   --------   ---------   -------------
    +       FOp1,FOp2   FOp=FOp1+FOp2      *       FOp1,FOp2   FOp=FOp1*FOp2
    +=      FOp1        FOp=FOp+FOp1       *=      FOp1        FOp=FOp*FOp1
    -       FOp1,FOp2   FOp=FOp1-FOp2      &=      FOp1        FOp=FOp1*FOp 
    -=      FOp1        FOp=FOp-FOp1
    -       FOp1        FOp=-FOp1                                            */

MSVCDLL friend floq2_op operator + (floq2_op &Op1 , floq2_op &Op2);
MSVCDLL void operator += (floq2_op &Op1);
MSVCDLL friend floq2_op operator - (floq2_op &Op1, floq2_op &Op2);
MSVCDLL friend floq2_op operator - (floq2_op &Op1);
MSVCDLL void operator -= (floq2_op &Op1);
MSVCDLL friend floq2_op operator * (floq2_op &Op1, floq2_op &Op2);
MSVCDLL void operator *= (floq2_op &Op1);
MSVCDLL void operator &= (floq2_op& Op1);

// ______________________________________________________________________
//      FLOQUET OPERATOR FUNCTIONS, FLOQUET OPERATOR WITH MATRIX
// ______________________________________________________________________

/* These functions handle all dealings between Floquet operators and Floquet
   dimension matrices.

 Function   Arguments      Return       Function   Arguments      Return
 --------   ---------   -------------   --------   ---------   -------------
    +       FOp1,Fmx    FOp=FOp1+Fmx       *       FOp1,FMx    FOp=FOp1*FMx
    +       Fmx,FOp1    FOp=FOp1+Fmx       *       FMx,FOp1    FOp=FMx*FOp1
    +=      FOp1        FOp=FOp+FOp1       *=      FMx         FOp=FOp*FMx
    -       FOp1,Fmx    FOp=FOp1-Fmx       &=      FMx         FOp=FMx*FOp 
    -=      Fmx         FOp=FOp-Fmx        -       Fmx,FOp1    FOp=Fmx-FOp1  */

MSVCDLL friend floq2_op operator + (floq2_op &Op1, matrix &mx);
MSVCDLL friend floq2_op operator + (matrix &mx, floq2_op &Op1);
MSVCDLL void operator += (const matrix& mx);
MSVCDLL friend floq2_op operator - (floq2_op &Op1, matrix &mx);
MSVCDLL friend floq2_op operator - (matrix &mx, floq2_op &Op1);
MSVCDLL void operator -= (const matrix& mx);
MSVCDLL friend floq2_op operator * (floq2_op &Op1, matrix &mx);
MSVCDLL friend floq2_op operator * (matrix &mx, floq2_op &Op1);
MSVCDLL void operator *= (matrix &mx);
MSVCDLL void operator &= (matrix &mx);

// ______________________________________________________________________
//       FLOQUET OPERATOR FUNCTIONS, FLOQUET OPERATOR WITH SCALAR
// ______________________________________________________________________

/* These functions handle all dealings between Floquet operators and scalar
   values (int, double, complex).

 Function   Arguments      Return       Function   Arguments      Return
 --------   ---------   ------------    --------   ---------  ----------------
    *        FOp1,z     FOp = z*FOp1       *=       this, z   void, z*this
    *        z,FOp1     FOp = z*FOp1       *=       this, d   void, d*this
    *        FOp1,d     FOp = d*FOp1        /       FOp1, z   FOp = (1/d)*this
    *        d,FOp1     FOp = d*FOp1       /=       this, z   void, (1/d)*this
                                                                               */

/*
MSVCDLL        floq_op operator *  (const complex &z) const;
MSVCDLL friend floq_op operator *  (const complex& z, const floq_op &FOp1);
MSVCDLL        floq_op operator *  (double d)         const;
MSVCDLL friend floq_op operator *  (double d,         const floq_op &Op1);
MSVCDLL        void    operator *= (const complex& z);
MSVCDLL        void    operator *= (double d);
MSVCDLL        floq_op operator /  (const complex &z) const;
MSVCDLL        floq_op operator /  (double d)         const;
MSVCDLL        void    operator /= (const complex& z);
*/

MSVCDLL friend floq2_op operator * (floq2_op &Op1, complex &z);
MSVCDLL friend floq2_op operator * (complex &z, floq2_op &Op1);
MSVCDLL friend floq2_op operator * (floq2_op &Op1, double z);
MSVCDLL friend floq2_op operator * (double z, floq2_op &Op1);
MSVCDLL        void     operator *= (const complex& z);
MSVCDLL        void     operator *= (double d);
MSVCDLL        floq2_op operator /  (const complex& z) const;
MSVCDLL        floq2_op operator /  (      double   d) const;
MSVCDLL        void     operator /= (const complex& z);
MSVCDLL        void     operator /= (      double   d);

// ______________________________________________________________________
//               COMPLEX  FLOQUET OPERATOR FUNCTIONS
// ______________________________________________________________________

//  friend gen_op pho_trace (floq2_op &Op);

	// Input		Op   : Floquet operator
        // Output          gen_op    : Operator containing trace(s) over    
        //                             Photon space

MSVCDLL friend floq2_op exp(floq2_op &Op1);
	
        // Input		Op1  : Floquet operator
        // Return		Op   : exponential of Op1
	//			       Op = exp(Op1)
        // Note			     : Computed in EBR of Op1

MSVCDLL friend floq2_op expm(floq2_op &Op1);
	
        // Input		Op1  : Floquet operator
        // Return		Op   : exponential of Op1 using Pade
	//			       Op = exp(Op1)

/*  friend floq2_op prop (floq2_op &FLOQHAM , double &time);
 
	// Input		Op1   : Floquet Hamiltonian
        //                      time  : Evolution time
        // Return		FP    : Floquet Propagator in DBR
	

  friend floq_op fprop (floq_op &FLOQHAM , double &time);
  
	// Input		Op1   : Floquet Hamiltonian
        //                      time  : Evolution time
        // Return		FP    : Floquet Propagator in DBR
	
*/
  
// ______________________________________________________________________
//             FLOQUET OPERATOR COMPONENT MANIPULATIONS
// ______________________________________________________________________

// ------------------ Floquet Submatrix Manipulations -----------------

/*  gen_op operator() (int row, int col);

	// Input		Op    : Floquet operator (this)
        //                      row   : Row index
        //                      col   : Column index
        // Output             gen_op  : at position (row,col) of Floquet Matrix
     
  gen_op get_block (int row, int col);

	// Input		Op    : Floquet operator (this)
        //                      row   : Row index
        //                      col   : Column index
        // Output             gen_op  : at position (row,col) of Floquet Matrix			      	
*/

MSVCDLL void put_block (gen_op &Op1, int r1, int c1, int r2, int c2);

	// Input		Op    : Floquet operator (this)
        //                      row   : Row index
        //                      col   : Column index
	//			Op1   : General operator
        // Output               none  : Op1 will be put at <row|Op|col> 

MSVCDLL void put_sdiag (gen_op &Op1, int sdn_1, int sdn_2);
        
       // Input                 Op    : Floquet operator (this)
       //                       Op1   : General operator
       //                      sdn_   : Number of sidediagonal to be set
       // Output             floq_op  : Floquet operator with sidediagonal
       //                             : set

  
// ----------------- Individual Element Manipulations -------------------
/*
  void put (complex &z, int row, int col );
       
       // Input                      : Floq_operator (this)
       //                 row,col    : Indices 
       // Output            none     : Value of <row|Op|col> in WBR
       //                            : set to z
       // Note                       : floq_op will be changed
 
  complex get ( int row, int col);


       // Input                      : Floq_operator (this)
       //                 row,col    : Indices 
       // Output            none     : Value of <row|Op|col> in WBR
       // Note                       : floq_op remains unchanged

  void put (complex &z, int N1_, int N2_, int H1_ , int H2_ );
       
       // Input                      : Floq_operator (this)
       //                 N1_,N2_    : Indices of photon space:-N...N
       //                 H1_,H2_    : Indices of Hilbert space:0...hs
       // Output            none     : Value of <N2,H2|Op|H1,N1> will be
       //                            : set to z
       // Note                       : floq_op will be changed
 
  complex get ( int N1_ , int N2_ , int H1_ , int H2_ );


       // Input                      : Floq_operator (this)
       //                 N1_,N2_    : Indices of photon space:-N...N
       //                 H1_,H2_    : Indices of Hilbert space:0...hs
       // Output                z    : complex value of <N2_,H2_|Op|H1_,N1_>
       // Note                       : floq_op remains unchanged
*/
// ______________________________________________________________________
//                  OPERATOR REPRESENTATION MANIPULATIONS
// ______________________________________________________________________


// *************  Set Operator to a Specific Representation *************
	
MSVCDLL void set_DBR();
MSVCDLL void set_EBR();

// ______________________________________________________________________
//                   CLASS FLOQUET OPERATOR I/O FUNCTION
// ______________________________________________________________________

	// Input		ostr : Output stream
	// 			Op   : Floquet operator
	// Return		     : Output stream into which Op is put

MSVCDLL friend std::ostream& operator << (std::ostream& ostr, const floq2_op &Op);

 // ______________________________________________________________________       
 //               INTERNAL FLOQUET OPERATOR MANIPULATIONS
 // ______________________________________________________________________        

        // Input                     : Floquet Operator (this)
        // Output add_omegas F_op    : Floquet operator with omegas added
        //                           : on main diagonal
        // Output sub_omegas F_op    : Floquet operator with omegas subtracted
        //                           : on main diagonal   
            
MSVCDLL void add_omegas();
MSVCDLL void sub_omegas();

};  

#endif 						// Floq2Op.h
