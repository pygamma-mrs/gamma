/* LSAux.h ******************************************************-*-c++-*-
**									**
**	                         G A M M A				**
**				          				**
**	Auxiliary Liouville Space Functions                 Interface	**
**						        		**
**	Copyright (c) 1993		        			**
**	Scott A. Smith		        				**
**	University of California, Santa Barbara				**
**	Department of Chemistry		        			**
**	Santa Barbara CA. 93106 USA	        			**
**					        			**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** Description								**
**									**
** Useful functions in dealing with Liuoville space entities.		**
**									**
*************************************************************************/

///Chapter Liouville Space Auxiliary Functions
///Section Overview
///Body The ...
///Section Available Liouville Space Auxiliary Functions

#ifndef GLS_aux_h_			// Is this file already included?
#  define GLS_aux_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Include system specifics
#include <LSLib/SuperOp.h>		// Know about superoperators
#include <Matrix/matrix.h>		// Know about matrices

 
// ____________________________________________________________________________
// **************** General Liouville Space Auxiliary Functions ***************
// ____________________________________________________________________________

  
MSVCDLL void print(const super_op& LOp, double cutoff=1.e-6, int nc=4, int ri=0);

	// Input		LOp   : Superoperator
	//			nc    : Number of columns
	//			ri    : Flag Real(0), Complex(1), Im(-1)
        // Return		None  : LOp elements sent to standard I/O
        // Note			      : Sets EBR of LOp

  
MSVCDLL void eigenvalues(super_op& LOp, int sort=1, int nc=4, int ri=0);

	// Input		LOp   : Superoperator
	//			sort  : Flag for sorting elements
	//			nc    : Number of columns
	//			ri    : Flag Real(0), Complex(1), Im(-1)
        // Return		None  : LOp eigenvalues sent to
	//				standard output
        // Note			      : Sets EBR of LOp


// ______________________________________________________________________
// *********************** COHERENCE SORTING ****************************
// ______________________________________________________________________

//??? Note: These need to be refined.  See note in the .cc file

MSVCDLL matrix UOrderMQC(const spin_sys &sys);

	// Input		sys	: Spin system
	// Output		U	: Transformation matrix
	//				  for coherence ordering


MSVCDLL super_op OrderMQC(const super_op &LOp, const matrix &U);

	// Input		LOp	: Superoperator
	//			U	: Transformation matrix
	// Output		R	: Superoperator LOp which has
	//				  been coherence ordered


MSVCDLL super_op OrderMQC(super_op &LOp, spin_sys &sys);

	// Input		LOp	: Superoperator
	//			sys	: Spin system
	//			R	: Superoperator LOp which has
	//				  been coherence ordered according
	//				  to sys default basis order


//super_op Hsuper(gen_op& Heff);

	// Input		Heff  : Effective Hamiltonian (Hz)
	// Output		LOp   : Hamiltonian commutation superoperator
	//
	//       <a,aa|LOp|b,bb> = del  * del     <b|H   |b> - <bb|H   |bb>
	//			      ab     aa,bb    eff           eff
	//
	// Note			      :	Lop is returned in angular frequency
	//				units


//super_op HsuperX(gen_op& Heff);

        // Input                Heff  : Effective Hamiltonian (Hz)
        // Output               LOp   : Hamiltonian commutation superoperator
        //
        //       <a,aa|LOp|b,bb> = del  * del     <b|H   |b> - <bb|H   |bb>
        //                            ab     aa,bb    eff           eff
        //
        //
        // Note                       : LOp is returned in angular frequency
        //                              units
 


// ______________________________________________________________________
// ********* MATRIX FUNCTIONS THAT SHOULD BE IN CLASS MATRIX ************
// ______________________________________________________________________


MSVCDLL matrix solve_it(matrix& X, matrix& Uguess, matrix& b, int lim=10);

	// Input		X     : Iteration matrix
	//			Uguess: Initial guess for U
	//			b     : Inhomogeneous part
	//			lim   : Maximum number of iterations allowed
        // Return		Ui    : Solution obtained by iteration to
	//
	//				     |U > = X|U   > + |b>
 	//			               i       i-1


MSVCDLL matrix invert_it(matrix& X);

	// Input		X     : Matrix
        // Return		Y     : Inverse of X
	// Note			      : Assumes X is Hermitian


MSVCDLL int LU_decomp(matrix& A, int* indx);

	// Input		A     : Input matrix
	//			indx  : Integer array for permutations
        // Return		A     : LU Decomposition of input A
        // 			indx  : Index of Permutations
	// Note			      : Assumes A is real, non-singular
	//				and square


MSVCDLL void LU_backsub(matrix &ALU, int* indx, matrix& b);

	// Input		ALU   : Input matrix A in LU form
	//			indx  : Integer array for permutations
	//			b     : Input column vector b
        // Return		      :
	// Note			      :


MSVCDLL matrix LU_invert(matrix& A);

	// Input		A     : Input matrix
        // Return		Ainv  : Inverse of input A
	// Note			      : Assumes A is , non-singular
	//				and square

#endif 					// LSAux.h
