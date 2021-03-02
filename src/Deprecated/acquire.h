/* acquire.h ********************************************-*-c++-*-
**								**
** 	               G A M M A				**
**								**
**	Class Acquire	                     Interface		**
**								**
**	Copyright (c) 1991, 1992, 1993		 		**
**	Tilo Levante, Scott Smith		 		**
**	Eidgenoessische Technische Hochschule	 		**
**	Labor fur physikalische Chemie		 		**
**	8092 Zurich / Switzerland		 		**
**								**
**      $Header: $
**								**
*****************************************************************/

/*****************************************************************
**								**
**  This file contains the implementation of class acquire, a	**
**  data type which facilitates the repetitive computation of	**
**  expectation values.  Use of this class can drastically	**
**  reduce computation time when multiple expectation values	**
**  of the type <Op(t)> = Trace<Op*sigma(t)> are desired for	**
**  over evenly spaced time increments.				**
**								**
*****************************************************************/

///Chapter Class Acquire
///Section Overview
///Body    The class Acquire contains the computational core
///        necessary for determining <Op(t)> = Tr{Op*sigma(t)}
///        repeatedly over evenly spaced time increments.
///        Here OpD is an operator (for a physical quantity)
///	   and Op(t) is another operator evolving in time
///        under a static Hamiltonian.  Use of Acquire can
///        increase computationaly speed and reduce code
///	   complexity.
///Section Available Acquire Functions


#ifndef   Gacquire_h_		// Is this file already included?
#  define Gacquire_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)	// Using the GNU compiler?
#    pragma interface		// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/Basis.h>
#include <Matrix/complex.h>

class matrix;			// Knowledge of matrices
class row_vector;		// Knowledge of row vectors
class gen_op;			// Knowledge of operators
class super_op;			// Knowledge of superoperators

class acquire 
  {
  int pos;			// Acquire dimensiton (# summation points)
  complex *A;			// Complex array (detector equivalent)
  complex *B;			// Complex array (time prop. equivalent)
  int *I,*J;			// Arrays for indices
  basis bs;			// Hilbert space computation basis
  matrix Sm1;			// Liouville space inverse eigenvectors
  matrix siginf;		// Infinite time density matrix
  complex trinf;		// Trace at infinite time
  int LS;			// Liouville space dimension
// sosi added dwell time for alpha, is it needed?
  double tdw;


// ----------------------------------------------------------------------
// ------------------------ PRIVATE FUNCTIONS ---------------------------
// ----------------------------------------------------------------------

// ______________________________________________________________________
//                    CLASS ACQUIRE ERROR HANDLING
// ______________________________________________________________________


// sosi - the alpha version has this constructor:
// create(gen_op& det, gen_op& prop, double dw)
// there is an acquire() function with these arguments

  void create(gen_op &det, gen_op &prop, int debug=0);

	// Input	det   : detection operator (in trace)
	//	 	prop  : propagation operator (in trace)
	// Output	this  : An acquire, parameters to facilitate
	//			repetetive acquisition calculations
	// Note		      : Liouville acquire - ls=0, J!=0, Sm1=0


  void create(gen_op& det, super_op& eLt, double cutoff=1.e-12);

	// Input	det   : detection operator (in trace)
	//	 	eLt   : propagation superoperator (in trace)
	// Output	this  : An acquire, parameters to facilitate
	//			repetetive acquisition calculations
	// Note		      : eLt must be a superoperator which
	//			operates directly on the density matrix
	//			for the incrementation time
	// Note		      : Liouville acquire - ls!=0, J=NULL, Sm1!=0


  void create(gen_op& det, super_op& eLt, gen_op sigma, int debug=0);

	// Input	det   : detection operator (in trace)
	//	 	eLt   : propagation superoperator (in trace)
	//		sigma : infinite time density operator
	// Output	this  : An acquire, parameters to facilitate
	//			repetetive acquisition calculations
	// Note		      : eLt must be a superoperator which
	//			operates directly on the density matrix
	//			for the incrementation time
	// Note		      : Redfield acquire - ls!=0, J=NULL, Sm1!=0


// ----------------------------------------------------------------------
// ------------------------- PUBLIC FUNCTIONS ---------------------------
// ----------------------------------------------------------------------

 public:

// ______________________________________________________________________
//                CLASS ACQUIRE CONSTRUCTION, DESTRUCTION
// ______________________________________________________________________

///Center Acquire Algebraic

// ------------------------ NULL Constructor ----------------------------

MSVCDLC acquire ();

	// Input	none  :
	// Output	ACQ   : A NULL acquire (this)
	///F_list acquire     - Constructor

// --------------- Hilbert Space Treatment Constructors -----------------

MSVCDLC acquire(gen_op &det, gen_op &prop);

	// Input	det   : Detection operator
	//		prop  : A propagation operator	
	// Output	ACQ   : Acquire (this) is constructed


MSVCDLC acquire(matrix &det, gen_op &prop);

	// Input	det   : Detection matrix
	//		prop  : A propagation operator	
	// Output	ACQ   : Acquire (this) is constructed


MSVCDLC acquire(gen_op &det, gen_op &H, double dw);

	// Input	det   : Detection operator (trace computation)
	//	 	H     : Evolution Hamiltonian (trace computation)
	//		dw    : Dwell time between acquisition points
	// None		ACQ   : Acquire (this) is constructed


MSVCDLC acquire(matrix &det, gen_op &H, double dw);

	// Input	det   : Detection matrix (trace computation)
	//	 	H     : Evolution Hamiltonian (trace computation)
	//		dw    : Dwell time between acquisition points
	// None		ACQ   : Acquire (this) is constructed

// -------------- Liouville Space Treatment Constructors ----------------

MSVCDLC acquire(gen_op &det, super_op &L);

	// Input	det   : Detection operator
	//		L     : A superoperator for time increment
	// Output	ACQ   : Acquire (this) is constructed



MSVCDLC acquire(matrix &det, super_op &L);

	// Input	det   : Detection matrix
	//		L     : A superoperator for time increment
	// Output	ACQ   : Acquire (this) is constructed


MSVCDLC acquire(gen_op &det, super_op &R, double dw);

	// Input	det   : Detection operator (trace computation)
	//	 	R     : Evolution Superoperator (trace computation)
	//		dw    : Dwell time between acquisition points
	// None		ACQ   : Acquire (this) is constructed


MSVCDLC acquire(matrix &det, super_op &R, double dw);

	// Input	det   : Detection matrix (trace computation)
	//	 	R     : Evolution Superoperator (trace computation)
	//		dw    : Dwell time between acquisition points
	// None		ACQ   : Acquire (this) is constructed


// ------------------ Redfield Treatment Constructors ----------------

MSVCDLC acquire(gen_op &det, super_op &L, gen_op& sigma);

	// Input	det   : Detection operator
	//		L     : A superoperator for time increment
	//		sigma : Infinite time density matrix
	// Output	ACQ   : Acquire (this) is constructed


 acquire(matrix &det, super_op &L, gen_op& sigma);

	// Input	det   : Detection matrix
	//		sigma : Infinite time density matrix
	//		L     : A superoperator for time increment
	// Output	ACQ   : Acquire (this) is constructed


MSVCDLC acquire(gen_op &det, super_op &R, gen_op& sigma, double dw);

	// Input	det   : Detection operator (trace computation)
	//		sigma : Infinite time density matrix
	//	 	R     : Evolution Superoperator (trace computation)
	//		dw    : Dwell time between acquisition points
	// None		ACQ   : Acquire (this) is constructed


MSVCDLC acquire(matrix &det, super_op &R, gen_op& sigma, double dw);

	// Input	det   : Detection matrix (trace computation)
	//		sigma : Infinite time density matrix
	//	 	R     : Evolution Superoperator (trace computation)
	//		dw    : Dwell time between acquisition points
	// None		ACQ   : Acquire (this) is constructed


// ----------------------- Self Constructor -----------------------------

MSVCDLC acquire (const acquire& ACQ1);

	// Input	ACQ1  : Acquire
	// None		ACQ   : Acquire (this) constructed from ACQ1

// --------------------------- Destructor -------------------------------

MSVCDLC ~acquire();

	// Input	ACQ   : An acquire (this)
	// Output	none  : ACQ is destructed

// -------------------------- Assignment --------------------------------

MSVCDLL void operator = (const acquire& ACQ1);

	// Input	ACQ1  : Acquire
	// None		ACQ   : Acquire (this) copied from ACQ1
	///F_list =           - Assignment


// ______________________________________________________________________
//                    EXPECTATION VALUE CALCULATIONS
// ______________________________________________________________________

///Center Basic Functions

MSVCDLL void operator () (gen_op &sigma, row_vector& fid, int np=0);

	// Input	ACQ   : An acquire (this)
	// 		sigma : Density matrix (operator propagated)
	//		fid   :	Data row_vector containing the result
	//		np    : Number of points to compute
	// Output	None  : fid data row_vector filled
	///F_list ()          - Acquistion

// This routine fills the data block fid according to the formula
 
//               t                    -1              t                  n  -1
//  fid(n) = <det *S|exp(L    *n*td)|S  *simga> = <det *S|[exp(L    *td)] |S  *sigma>
//                        diag                                  diag
 
// where det is a detection operator, sigma the state of the system at the
// start of some series of evolution steps by time increment td, and the
// exponential Liouvillian L is
 
//                                                    -1  
//                         exp(L*t) = S*exp(L    *t)*S
//                                           diag
 
// In terms of class acquire, this is then
 
//                         n                                th              -1
//  fid(n) = Sum { A[i] * B [i] * C[i] } where C[i] is the i   element of |S  * sigma>
 
// Variation 1: Hilbert space treatment - L is already diagonal, S = inv(S) = I
// Variation 2: Redfield treatment      - It is not sigma which evolved but sigma-sigmainf
//					  thus component from this must be added and subtracted.
//					  (Strictly true only if <det|exp(L*t)|sigmainf> !=0)
 

MSVCDLL void operator () (gen_op& sigma, matrix& mx, double cutoff=1.e-10);

	// Input	ACQ   : An acquire (this)
	// 		sigma : Density matrix (operator propagated)
        //              mx    : A matrix for output data
        //              cutoff: An intensity cutoff, contributions with
        //                      norms below this value are not output
        // Output       void  : The matrix mx is output as a posx2 array.
	//			Each row corresponds to one of the pos components
	//			detected according to det, evolving according to
	//			exp(Lt), and beginning in a state described by sigma.
	//			Component i is a complex exponential oscillating
	//			at frequency Im(mx(i,0)) & decaying at rate Re(mx(i,0))
	//			with complex intensity mx(i,1)
        // Note               : No sorting is performed!  No components are combined!
        // Note               : If ACQ was formed with exp(L*td) with td != 1, then
	//		                     THIS ROUTINE WILL FAIL
 
// This routine fills the matrix mx with components of the time evolved signal
// fid(t) as defined by the formula
 
//                                t                 -1
//                   fid(t) = <det *S|exp(L    *t)|S  *simga>
//                                         diag
 
// where det is a detection operator, sigma the state of the system at the start
// of some evolution, and the exponential Liouvillian L which describes the evolution is
 
//                                                    -1  
//                        exp(L*t) = S*exp(L    *t)*S
//                                          diag
 
// In terms of class acquire, this is then
 
//                   fid(t) = Sum { A[i] * B[i]*t * C[i] }
//                             i
//                    th              -1
// where C[i] is the i   element of |S  * sigma>.  The matrix mx contains a row
// for each signal component, i.e. for each i.  Each row, or each component, has
// two parts: column 0 = B[i], column 1 = A[i]*C[i].
 
//               mx(i,0) = B[i]           mx(i,1) = A[i] * C[i]

//  Since B[i]*t are complex exponentials, column 0 contains exponential precession
// (imaginaries) and decay (real) rates.  Column 1 contains the complex exponential
// intensities (magnitudes and phases).  

// Note that an "analog" frequency domain signal can be obtained by Fourier transformation
// of the above equations.  The result would be that the values in the matrix then
// correspond to complex Lorentzians with frequencies and witdhs in column 0 and the
// intensities in column 1.
 
// Variation 1: Hilbert space treatment - L is already diagonal, S = inv(S) = I
// Variation 2: Redfield treatment      - It is not sigma which evolved but sigma-sigmainf
//					  thus component from this must be added and subtracted.
//					  (Strictly true only if <det|exp(L*t)|sigmainf> !=0)
 
// ______________________________________________________________________
//                  CLASS ACQUIRE AUXILIARY FUNCTIONS
// ______________________________________________________________________


MSVCDLL int ls() const;

	// Input	ACQ   : An acquire (this)
	// Output	ls    : Liouville space size


MSVCDLL int size();

	// Input	ACQ   : An acquire (this)
	// Output	size  : Size of the arrays A & B
	///F_list size        - Reduced space size


MSVCDLL int full_size();

	// Input	ACQ   : An acquire (this)
	// Output	size  : Full size of operators
	//			the Hilbert space squared
	///F_list full_size   - Full space size


// ____________________________________________________________________________
//                        CLASS ACQUIRE I/O FUNCTIONS
// ____________________________________________________________________________

// ----------------- ASCII Output Functions (For Debugging) -------------------

MSVCDLL std::ostream& print(std::ostream &ostr);

        // Input                ACQ   : Acquire
        //                      ostr  : Output stream
        // Output               none  : Acquire vector is sent
        //                              to the output stream


MSVCDLL std::ostream& print(std::ostream &ostr, gen_op& sigma);

        // Input                ACQ   : Acquire
        //                      ostr  : Output stream
        // Output               none  : Acquire vector is sent
        //                              to the output stream


MSVCDLL friend std::ostream &operator << (std::ostream &ostr, acquire &ACQ);

	// Input		ACQ   : Acquire
        //                      ostr  : Output stream
        // Output               none  : Acquire vector is sent
	//				to the output stream


MSVCDLL std::ostream& printT(std::ostream &ostr, gen_op& sigma, int np);

        // Input                ACQ   : Acquire
        //                      ostr  : Output stream
        //                      sigma : Density matrix (operator propagated)
        //                      np    : Number of points to compute
        // Output               ostr  : The output stream modified by
        //                                the acquisition time domain
        //                                calculation

// ______________________________________________________________________
//                  CLASS ACQUIRE OUTDATED FUNCTIONS
// ______________________________________________________________________


MSVCDLL void operator () (gen_op &sigma, gen_op &sigma0, super_op &R,
						 row_vector& fid, int np=0);

	// Input	ACQ   : An acquire (this)
	// 		sigma : Density matrix (operator propagated)
	// 		sigma0: Equilibrium matrix
	//		R     : Relaxation superoperator
	//		fid   :	Data row_vector containing the result
	//		np    : Number of points to compute
	// Output	None  : Acquisition "fid" data row_vector filled
	// Note		      : This routine should utilze a relaxation
	//			superoperator R which utilizes the
	//			secular approximation!
  

MSVCDLL void operator () (gen_op &sigma, gen_op &sigma0, super_op &R,
				row_vector& fid, double offset, int np=0);

	// Input	ACQ   : An acquire (this)
	// 		sigma : Density matrix (operator propagated)
	// 		sigma0: Equilibrium matrix
	//		R     : Relaxation superoperator
	//		fid   :	Data row_vector containing the result
	//		offset: RF field frequency offset
	//		np    : Number of points to compute
	// Output	None  : Acquisition "fid" data row_vector filled
	// Note		      : This routine should utilze a relaxation
	//			superoperator R which utilizes the
	//			secular approximation!

};

#endif						// acquire.h
