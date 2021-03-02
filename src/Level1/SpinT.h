/* SpinT.h **********************************************-*-c++-*-
**								**
** 	                   G A M M A				**
**								**
**	Spin Tensors	                       Interface   	**
**							 	**
**	Copyright (c) 1991, 1992			 	**
**	Scott Smith					        **
**	Eidgenoessische Technische Hochschule	 		**
**	Labor fur physikalische Chemie		 		**
**	8092 Zurich / Switzerland		 		**
**						 		**
**      $Header: $
**                						**
*****************************************************************/

/*****************************************************************
**								**
** Description:							**
**								**
** This file contains the definition and workings of spin	**
**  tensors.							**
**								**
*****************************************************************/

///Chapter Class Spin Tensor (spin_T)
///Section Overview
///Body The class Spin Tensor (spin_T) embodies spin tensors
///Section Available Spin Tensor Functions      

#ifndef   Gspin_T_h_			// Is file already included?
#  define Gspin_T_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Level1/coord.h>		// Include coordinates
#include <HSLib/GenOp.h>		// Include general operators
#include <HSLib/SpinOp.h>		// Include spin operators
#include <Level1/SpaceT.h>		// Include spatial tensors


// Forward declaration of spin_T class
// for the benefit of following functions.
class spin_T;

// Added function declarations to augment friend declartions
// that are included in the class.

MSVCDLL spin_T T1(const spin_sys &sys, int spin);
 
MSVCDLL spin_T T11(const spin_sys &sys, int spin);

MSVCDLL spin_op T1(const spin_sys &sys, int spin, int l, int m);

MSVCDLL spin_op T10(const spin_sys &sys, int spin, int m);

MSVCDLL spin_op T10(spin_op &Ie, int m);

MSVCDLL spin_op T11(const spin_sys &sys, int spin, int m);

MSVCDLL spin_op T11(spin_op &Im, spin_op &Iz, spin_op &Ip, int m);

MSVCDLL spin_T T2(const spin_sys &sys, int spin1, int spin2);

MSVCDLL spin_T T22wh(const spin_sys &sys, int spin1, int spin2);

MSVCDLL spin_T T22(const spin_sys &sys, int spin1, int spin2);

MSVCDLL spin_T T22(const spin_sys &sys, spin_op &Im1, spin_op &Iz1, spin_op &Ip1,
			           spin_op &Im2, spin_op &Iz2, spin_op &Ip2);

MSVCDLL spin_op T2(const spin_sys &sys, int spin1, int spin2, int l, int m);

MSVCDLL spin_op T2(spin_op &Im1, spin_op &Iz1, spin_op &Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int l, int m);

MSVCDLL spin_op T20(const spin_sys &sys, int spin1, int spin2, int m);

MSVCDLL spin_op T20(spin_op &Im1, spin_op &Iz1, spin_op& Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int m);

MSVCDLL spin_op T21(const spin_sys &sys, int spin1, int spin2, int m);

MSVCDLL spin_op T21(spin_op &Im1, spin_op &Iz1, spin_op& Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int m);

MSVCDLL spin_op T22(const spin_sys &sys, int spin1, int spin2, int m);

MSVCDLL spin_op T22(spin_op &Im1, spin_op &Iz1, spin_op& Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int m);

MSVCDLL spin_T T2(const spin_sys &sys, int spin, const coord &vect);
MSVCDLL spin_T T2(const spin_sys &sys, spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect);

MSVCDLL spin_T T2SS(const spin_sys &sys, int spin, const coord &vect, int rev=0);

MSVCDLL spin_T T2SS(const spin_sys &sys, spin_op &Im, spin_op &Iz,
				      spin_op &Ip, const coord &vect, int rev=0);

MSVCDLL spin_T T22(const spin_sys &sys, int spin, const coord &vect);
MSVCDLL spin_T T22(const spin_sys &sys, spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect);

MSVCDLL spin_T T22SSirr(const spin_sys &sys, int spin, const coord &vect, int rev=0);

MSVCDLL spin_T T22SSirr(const spin_sys &sys, spin_op &Im, spin_op &Iz,
					 spin_op &Ip, const coord &vect, int rev=0);

MSVCDLL spin_op T2(const spin_sys &sys, int spin, const coord &vect, int l, int m);
MSVCDLL spin_op T2(spin_op &Im, spin_op &Iz, spin_op &Ip,
 					   const coord &vect, int l, int m);
MSVCDLL spin_op T20(const spin_sys &sys, int spin, const coord &vect, int m);
MSVCDLL spin_op T20(spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect, int m);
MSVCDLL spin_op T21(const spin_sys &sys, int spin, const coord &vect, int m);
MSVCDLL spin_op T21(spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect, int m);
MSVCDLL spin_op T22(const spin_sys &sys, int spin, const coord &vect, int m);
MSVCDLL spin_op T22(spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect, int m);

MSVCDLL spin_op T2SS(const spin_sys &sys, int spin, const coord &vect, int l, int m, int rev=0);

MSVCDLL spin_op T2SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
 				   const coord &vect, int l, int m, int rev=0);

MSVCDLL spin_op T20SS(const spin_sys &sys, int spin, const coord &vect, int m, int rev=0);

MSVCDLL spin_op T20SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
						 const coord &vect, int m, int rev=0);

MSVCDLL spin_op T21SS(const spin_sys &sys, int spin, const coord &vect, int m, int rev=0);

MSVCDLL spin_op T21SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
					 const coord &vect, int m, int rev=0);

MSVCDLL spin_op T22SS(const spin_sys &sys, int spin, const coord &vect, int m, int rev=0);

MSVCDLL spin_op T22SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
						const coord &vect, int m, int rev=0);

MSVCDLL spin_op T_comp(spin_T &SphT, int L, int M);
// ?? Remove this, replaced by member function

MSVCDLL spin_T T_mult(spin_T &SphT1, spin_T &SphT2);

MSVCDLL spin_op T_mult(spin_T &SphT1, spin_T &SphT2, int L, int M);

MSVCDLL spin_T T_rot(spin_T &SphT1, double alpha, double beta, double gamma);
// ??? this function should be removed, replaced by the member function!

MSVCDLL spin_op T_rot(spin_T &SphT1, int l, int m,
			 double alpha, double beta, double gamma);
// ??? this function should be removed, replaced by the member function!

MSVCDLL spin_op T_prod(spin_T &SphT, space_T &SphA, int l, int m);
MSVCDLL spin_op T_prod(space_T &SphA, spin_T &SphT, int l, int m);

MSVCDLL spin_op T_prod(spin_T &SphT, space_T &SphA, int l);
MSVCDLL spin_op T_prod(space_T &SphA, spin_T &SphT, int l);

MSVCDLL spin_op T_prod(spin_T &SphT, space_T &SphA);
MSVCDLL spin_op T_prod(space_T &SphA, spin_T &SphT);

MSVCDLL double Clebsch_Gordan(int a, int b, int alpha, int beta, int c, int gamma);
MSVCDLL double Wigner_3j(int a, int b, int c, int alpha, int beta, int gamma);



class spin_T
  {

// ----------------------------------------------------------------------------
// ------------------------------ STRUCTURE -----------------------------------
// ----------------------------------------------------------------------------

  spin_sys* sys; 		// Pointer to a spin system
  int rank;			// Rank of spin tensor
  spin_op*** pr;		// pointer to pointer to SOp pointers

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                       SPIN TENSOR ERROR HANDLING
// ____________________________________________________________________________


void volatile spin_T_error(const int i);

	// Input		i    : Error Flag
	// Output		none : Error Message Output


void volatile spin_T_fatality (int error);

	// Input		none :
	// Output		none : Stops Execution & Error Message Output



// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

  public:

// ____________________________________________________________________________
// A                    SPIN TENSOR CONSTRUCTORS, DESTRUCTOR
// ____________________________________________________________________________

///Center Spin Tensor Algebraic

MSVCDLC spin_T();

	// Input		none :
	// Return		SphT : nullspin tensor
	///F_list spin_T	     - Spin tensor constructor


MSVCDLC spin_T(const spin_sys &sys);

	// Input		sys  : spin system
	// Return		SphT : spin tensor associated with
        //	                       the input spin system


MSVCDLC spin_T(const spin_T &SphT);

	// Input		SphT : spin tensor
	// Return		SphT1: spin tensor duplicate of input
        //	                       spin tensor


MSVCDLC spin_T(const spin_T &SphT, int l);

	// Input		SphT : spin tensor
	//                      l    : rank
	// Return		SphT1: spin tensor duplicate of input
        //	                       spin tensor


MSVCDLC ~spin_T();

	// Input		SphT : Sperical Tensor
	// Return		     : none, deletes current Tensor


// ____________________________________________________________________________
//                         SPIN TENSOR UNARY OPERATIONS
// ____________________________________________________________________________


MSVCDLL spin_T & operator = (const spin_T &SphT);

	// Input		SphT : spin tensor
	// Return		SphT1: spin tensor equivalent to the
        //	                       input spin tensor
	// Note		             : SphT and SphT1 must be associated
	//			       with the same spin system
	///F_list =		     - Spin tensor assignment 
 

// ____________________________________________________________________________
// B           GENERAL RANK 1 SPHERICAL SPIN TENSOR FUNCTIONS
// ____________________________________________________________________________

/// Rank 1 Spin Tensor Functions

// ********************  Full Rank One Treatment ************************


MSVCDLL friend spin_T T1(const spin_sys &sys, int spin);

	// Input		sys   : spin system
	// 			spin  : spin index
	// Output		SphT  : rank 1 Spin Tensor for spin specified
	///F_list T1		      - Spin tensor, rank 1


// *****************  Irreducible Rank One Spin Tensor ******************

 
MSVCDLL friend spin_T T11(const spin_sys &sys, int spin);

	// Input		sys   : spin system
	// 			spin  : spin index
	// Output		SphT  : rank 1 Irreducible Spin Tensor for
	//				spin specified
	///F_list T11		      - Spin tensor, rank 1 irreducible
 

// **********  Specific Irreducible Rank One Tensor Components **********


MSVCDLL friend spin_op T1(const spin_sys &sys, int spin, int l, int m);

	// Input		sys   : spin system
	// 			spin  : spin index
	//			l     : component index
	//			m     : component index
	// Output		SOp   : rank 1 Spin Tensor component Tlm
	//				for spins specified
	///F_list T1		      - Spin tensor, rank 1 irreducible component
 

MSVCDLL friend spin_op T10(const spin_sys &sys, int spin, int m);

	// Input		sys   : spin system
	// 			spin  : spin index
	//			m     : component index
	// Output		SOp   : rank 1 Spin Tensor component T00
	//				for spins specified
	///F_list T10		      - Spin tensor, rank 1 m = 0 irreducible component


MSVCDLL friend spin_op T10(spin_op &Ie, int m);

	// Input		Ie    : Identity matrix
	//			m     : Component index		      1
	// Output		SOp   : Rank 1 spin tensor component T
	//							      0m

MSVCDLL friend spin_op T11(const spin_sys &sys, int spin, int m);

	// Input		sys   : spin system
	// 			spin  : spin index
	//			m     : component index
	// Output		SOp   : rank 1 Spin Tensor component T1m
	//				for spins specified
	///F_list T10		      - Spin tensor, rank 1 m = 1 irreducible component


MSVCDLL friend spin_op T11(spin_op &Im, spin_op &Iz, spin_op &Ip, int m);

	// Input		Im    : Spin operator I-
	// 			Iz    : Spin operator Iz
	// 			Ip    : Spin operator I+
	//			m     : Component index		      1
	// Output		SOp   : Rank 1 spin tensor component T
	//							      1m
 

// ______________________________________________________________________
// C            GENERAL RANK 2 SPHERICAL SPIN TENSOR FUNCTIONS
// ______________________________________________________________________

/// Rank 2 Spin Tensor Functions

// ********************  Full Rank Two Treatment ************************


MSVCDLL friend spin_T T2(const spin_sys &sys, int spin1, int spin2);

	// Input		sys   : spin system
	// 			spin1 : spin index
	// 			spin2 : spin index
	// Output		SphT  : rank 2 Spin Tensor for spins specified
	///F_list T2		      - Spin tensor, rank 2 general
 

// *****************  Irreducible Rank Two Spin Tensor ******************


MSVCDLL friend spin_T T22wh(const spin_sys &sys, int spin1, int spin2);

        // Input                sys   : Spin system
        //                      spin1 : Spin index
        //                      spin2 : Spin index                      2
        // Output               SOp   : Irreducible rank 2 spin tensor T (ij)
        //                                                              2
        // Note                       : This routine will leave all but
        //                              the m=0 component empty if spin1 and
        //                              spin 2 are not of the same spin type
 

MSVCDLL friend spin_T T22(const spin_sys &sys, int spin1, int spin2);

	// Input		sys   : spin system
	// 			spin1 : spin index
	// 			spin2 : spin index
	// Output		SphT  : irreducible rank 2 Spin Tensor
	//				for spins specified
	///F_list T22		      - Spin tensor, rank 2 general irreducible


MSVCDLL friend spin_T T22(const spin_sys &sys, spin_op &Im1, spin_op &Iz1, spin_op &Ip1,
			           spin_op &Im2, spin_op &Iz2, spin_op &Ip2);

	// Input		sys   : Spin system
	// 			Im1   : I- spin 1
	// 			Iz1   : Iz spin 1
	// 			Ip1   : I+ spin 1
	// 			Im2   : I- spin 2
	// 			Iz2   : Iz spin 2
	// 			Ip2   : I+ spin 2			2
	// Output		SOp   : Irreducible rank 2 spin tensor T (12)
	//							        2


// **********  Specific Irreducible Rank Two Tensor Components **********


MSVCDLL friend spin_op T2(const spin_sys &sys, int spin1, int spin2, int l, int m);

	// Input		sys   : spin system
	// 			spin1 : spin index
	// 			spin2 : spin index
	//			l     : rank index
	//			m     : component index
	// Output		SOp   : rank 2 Spin Tensor component Tlm
	//				for spins specified
	///F_list T22		      - Spin tensor, rank components


MSVCDLL friend spin_op T2(spin_op &Im1, spin_op &Iz1, spin_op &Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int l, int m);

	// Input		sys   : spin system
	// 			Im1   : I- spin 1
	// 			Iz1   : Iz spin 1
	// 			Ip1   : I+ spin 1
	// 			Im2   : I- spin 2
	// 			Iz2   : Iz spin 2
	// 			Ip2   : I+ spin 2
	//			l     : component index
	//			m     : Component index               2
	// Output		SOp   : Rank 2 Spin Tensor component T  (12)
	//							      lm
 

MSVCDLL friend spin_op T20(const spin_sys &sys, int spin1, int spin2, int m);

	// Input		sys   : spin system
	// 			spin1 : spin index
	// 			spin2 : spin index
	//			m     : component index
	// Output		SOp   : rank 2 Spin Tensor component T00
	//				for spins specified
 

MSVCDLL friend spin_op T20(spin_op &Im1, spin_op &Iz1, spin_op& Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int m);

	// Input		Im1   : I- spin 1
	// 			Iz1   : Iz spin 1
	// 			Ip1   : I+ spin 1
	// 			Im2   : I- spin 2
	// 			Iz2   : Iz spin 2
	// 			Ip2   : I+ spin 2
	//			m     : Component index               2
	// Output		SOp   : Rank 2 Spin Tensor component T  (12)
	//							      0m
	///F_list T20		      - Spin tensor, rank components
 

MSVCDLL friend spin_op T21(const spin_sys &sys, int spin1, int spin2, int m);

	// Input		sys   : spin system
	// 			spin1 : spin index
	// 			spin2 : spin index
	//			m     : component index
	// Output		SOp   : rank 2 Spin Tensor component T1m
	//				for spins specified
	///F_list T21		      - Spin tensor, rank components


MSVCDLL friend spin_op T21(spin_op &Im1, spin_op &Iz1, spin_op& Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int m);

	// Input		Im1   : I- spin 1
	// 			Iz1   : Iz spin 1
	// 			Ip1   : I+ spin 1
	// 			Im2   : I- spin 2
	// 			Iz2   : Iz spin 2
	// 			Ip2   : I+ spin 2
	//			m     : Component index               2
	// Output		SOp   : Rank 2 Spin Tensor component T  (12)
	//							      1m
 

MSVCDLL friend spin_op T22(const spin_sys &sys, int spin1, int spin2, int m);

	// Input		sys   : spin system
	// 			spin1 : spin index
	// 			spin2 : spin index
	//			m     : component index
	// Output		SOp   : rank 2 Spin Tensor component T2m
	//				for spins specified
	///F_list T22		      - Spin tensor, rank components
 

MSVCDLL friend spin_op T22(spin_op &Im1, spin_op &Iz1, spin_op& Ip1,
		spin_op &Im2, spin_op &Iz2, spin_op &Ip2, int m);

	// Input		Im1   : I- spin 1
	// 			Iz1   : Iz spin 1
	// 			Ip1   : I+ spin 1
	// 			Im2   : I- spin 2
	// 			Iz2   : Iz spin 2
	// 			Ip2   : I+ spin 2
	//			m     : Component index               2
	// Output		SOp   : Rank 2 Spin Tensor component T  (12)
	//							      2m


// ______________________________________________________________________
// D         GENERAL RANK 2 SPHERICAL SPIN-SPACE TENSOR FUNCTIONS
// ______________________________________________________________________


// ********************  Full Rank Two Treatment ************************


MSVCDLL friend spin_T T2(const spin_sys &sys, int spin, const coord &vect);
MSVCDLL friend spin_T T2(const spin_sys &sys, spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect);

MSVCDLL friend spin_T T2SS(const spin_sys &sys, int spin, const coord &vect, int rev);

	// Input		sys   : spin system
	// 			spin  : spin index
	// 			vect  : Cartesian vector
	// 			rev   : Reverse tensor order
	// Output		SphT  : rank 2 "Spin" Tensor for spin specified


MSVCDLL friend spin_T T2SS(const spin_sys &sys, spin_op &Im, spin_op &Iz,
				      spin_op &Ip, const coord &vect, int rev);

	// Input                sys   : Spin system
   	//                      Im    : Spin operator I-
   	//                      Iz    : spin operator Iz
   	//                      Ip    : Spin operator I+
   	//                      vect  : Cartesian vector
	// 			rev   : Reverse tensor order
   	// Output               SphT  : rank 2 "Spin" Tensor for spin specified
 

// *****************  Irreducible Rank Two Spin Tensor ******************


MSVCDLL friend spin_T T22(const spin_sys &sys, int spin, const coord &vect);
MSVCDLL friend spin_T T22(const spin_sys &sys, spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect);

MSVCDLL friend spin_T T22SSirr(const spin_sys &sys, int spin, const coord &vect, int rev);

	// Input		sys   : spin system
	// 			spin  : spin index
	// 			vect  : Cartesian vector
	// 			rev   : Reverse tensor order
	// Output		SphT  : irreducible rank 2 "Spin" Tensor
	//				for spin specified


MSVCDLL friend spin_T T22SSirr(const spin_sys &sys, spin_op &Im, spin_op &Iz,
					 spin_op &Ip, const coord &vect, int rev);

  	// Input                sys   : Spin system
  	//                      Im    : Spin operator I-
  	//                      Iz    : spin operator Iz
  	//                      Ip    : Spin operator I+
  	//                      vect  : Cartesian vector
	// 			rev   : Reverse tensor order
  	// Output               SphT  : irreducible rank 2 "Spin" Tensor
  	//                              for spin specified


// *****  Specific Irreducible Rank Two Space-Spin Tensor Components ****


MSVCDLL friend spin_op T2(const spin_sys &sys, int spin, const coord &vect, int l, int m);
MSVCDLL friend spin_op T2(spin_op &Im, spin_op &Iz, spin_op &Ip,
 					   const coord &vect, int l, int m);
MSVCDLL friend spin_op T20(const spin_sys &sys, int spin, const coord &vect, int m);
MSVCDLL friend spin_op T20(spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect, int m);
MSVCDLL friend spin_op T21(const spin_sys &sys, int spin, const coord &vect, int m);
MSVCDLL friend spin_op T21(spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect, int m);
MSVCDLL friend spin_op T22(const spin_sys &sys, int spin, const coord &vect, int m);
MSVCDLL friend spin_op T22(spin_op &Im, spin_op &Iz, spin_op &Ip, const coord &vect, int m);

MSVCDLL friend spin_op T2SS(const spin_sys &sys, int spin, const coord &vect, int l, int m, int rev);

	// Input		sys   : spin system
	// 			spin  : spin index
	//			vect  : Cartesian vector
	//			l     : rank index
	//			m     : component index
	// 			rev   : Reverse tensor order
	// Output		SOp   : rank 2 "Spin" Tensor component Tlm
	//				for spin specified


MSVCDLL friend spin_op T2SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
 				   const coord &vect, int l, int m, int rev);

  	// Input		Im    : Spin operator I-
  	//                      Iz    : spin operator Iz
  	//                      Ip    : Spin operator I+
  	//                      vect  : Cartesian vector
	//			l     : rank index
	//			m     : component index
	// 			rev   : Reverse tensor order
	// Output		SOp   : rank 2 "Spin" Tensor component Tlm
	//				for spin specified


MSVCDLL friend spin_op T20SS(const spin_sys &sys, int spin, const coord &vect, int m, int rev);

	// Input		sys   : Spin system
	// 			spin  : Spin index
	//			vect  : Cartesian vector
	//			m     : Component index
	// 			rev   : Reverse tensor order (unused)
	// Output		SOp   : Rank 2 "Spin" Tensor component T0m
	//				for spin specified


MSVCDLL friend spin_op T20SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
						 const coord &vect, int m, int rev);

  	// Input		Im    : Spin operator I-
  	//                      Iz    : Spin operator Iz
  	//                      Ip    : Spin operator I+
	//			vect  : Cartesian vector
	//			m     : Component index
	// 			rev   : Reverse tensor order (unused)
	// Output		SOp   : Rank 2 "Spin" Tensor component T0m
	//				for spin specified


MSVCDLL friend spin_op T21SS(const spin_sys &sys, int spin, const coord &vect, int m, int rev);

	// Input		sys   : Spin system
	// 			spin  : Spin index
	//			vect  : Cartesian vector
	//			m     : Component index
	// 			rev   : Reverse tensor order
	// Output		SOp   : Rank 2 "Spin" Tensor component T1m
	//				for spin specified


MSVCDLL friend spin_op T21SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
					 const coord &vect, int m, int rev);

  	// Input		Im    : Spin operator I-
  	//                      Iz    : Spin operator Iz
  	//                      Ip    : Spin operator I+
	//			vect  : Cartesian vector
	//			m     : Component index
	// 			rev   : Reverse tensor order
	// Output		SOp   : Rank 2 "Spin" Tensor component T1m
	//				for spin specified


MSVCDLL friend spin_op T22SS(const spin_sys &sys, int spin, const coord &vect, int m, int rev);

	// Input		sys   : Spin system
	// 			spin  : Spin index
	//			vect  : Cartesian vector
	//			m     : Component index
	// 			rev   : Reverse tensor order (unused)
	// Output		SOp   : Rank 2 "Spin" Tensor component T2m
	//				for spin specified


MSVCDLL friend spin_op T22SS(spin_op &Im, spin_op &Iz, spin_op &Ip,
						const coord &vect, int m, int rev);

  	// Input		Im    : Spin operator I-
  	//                      Iz    : Spin operator Iz
  	//                      Ip    : Spin operator I+
	//			vect  : Cartesian vector
	//			m     : Component index
	// 			rev   : Reverse tensor order (unused)
	// Output		SOp   : Rank 2 "Spin" Tensor component T2m
	//				for spin specified


// ______________________________________________________________________
// E                      SPIN TENSOR FUNCTIONS
// ______________________________________________________________________


MSVCDLL friend spin_op T_comp(spin_T &SphT, int L, int M);
// ?? Remove this, replaced by member function

MSVCDLL spin_op component(int l, int m);

	// Input		SphT : Spin Tensor(this)
	// 			l    : momentum index
	// 			m    : momentum index
	// Output		SOp  : Spin Operator, the l,m
	//			     : component of Spin Tensor SphT
	///F_list T_comp	     - Spin tensor, rank components

// ______________________________________________________________________
// F                    SPIN TENSOR MULTIPLICATIONS
// ______________________________________________________________________
 

MSVCDLL friend spin_T T_mult(spin_T &SphT1, spin_T &SphT2);

	// Input		SphT1 : Spin Tensor
	// 			SphT2 : Spin Tensor
	// Output		SphT  : Spin Tensor which is the product
	//			        of the two input tensors
	//				SphT = SphT1 x SphT2
	///F_list T_mult	      - Spin tensor, rank components
 

MSVCDLL friend spin_op T_mult(spin_T &SphT1, spin_T &SphT2, int L, int M);

	// Input		SphT1 : Spin Tensor
	// 			SphT2 : Spin Tensor
	// Output		SOp   : Spin operator which is the L,M
	//				component of tensor SphT, the
	//			        product of the two input tensors
	//				SphT = SphT1 x SphT2
 
 
// ______________________________________________________________________
// G                     SPIN TENSOR ROTATIONS
// ______________________________________________________________________


MSVCDLL friend spin_T T_rot(spin_T &SphT1, double alpha, double beta, double gamma);
// ??? this function should be removed, replaced by the member function!

MSVCDLL friend spin_op T_rot(spin_T &SphT1, int l, int m,
			 double alpha, double beta, double gamma);
// ??? this function should be removed, replaced by the member function!


MSVCDLL spin_T rotate(double alpha, double beta, double gamma);

	// Input		SphT1 : Spin Tensor(this)
	// 			alpha : Euler Angle
	// 			beta  : Euler Angle
	// 			gamma : Euler Angle
	// Output		SphT  : Spin Tensor which is the input
	//			        tensor in the coordinate system
	//				rotated by input Euler angles
	///F_list rotate	      - Spin tensor rotation
 

MSVCDLL spin_T rotate(const coord &EA);

	// Input		SphT1 : Spin Tensor(this)
	// 			EA    : Set of Euler angles
	// Output		SphT  : Spin Tensor which is the input
	//			        tensor in the coordinate system
	//				rotated by input Euler angles


MSVCDLL spin_op rotate(int l, int m, double alpha, double beta, double gamma);

	// Input		SphT1 : Spin Tensor
	//			l     : Component index
	//			m     : Component index
	// 			alpha : Euler Angle
	// 			beta  : Euler Angle
	// 			gamma : Euler Angle
	// Output		Alm   : Spatial Tensor component of
	//			        the input tensor in the coordinate
	//				system rotated by input Euler angles


MSVCDLL spin_op rotate(int l, int m, const coord &EA);

	// Input		SphT1 : Spin Tensor(this)
	//			l     : Component index
	//			m     : Component index
	// 			EA    : Set of Euler angles
	// Output		Tlm   : Spatial Tensor component of
	//			        the input tensor in the coordinate
	//				system rotated by input Euler angles
 

MSVCDLL friend spin_op T_prod(spin_T &SphT, space_T &SphA, int l, int m);
MSVCDLL friend spin_op T_prod(space_T &SphA, spin_T &SphT, int l, int m);

	// Input		SphT  : Irreducible Spin Tensor
	// 			SphA  : Irreducible Spatial Tensor
	//			l     : Tensor Rank Index
	//			m     : Tensor Component Index
	// Output		SOp   : Spin operator which is the result
	//				of the scalar product of  the l,m 
	//			        component ot SphT with the l-m of SphA
	//				      m   
	//				  (-1)  * T   * A
	//				           lm    l-m
	//				or equivalently
	//				      m   
	//				  (-1)  * A   * T
	//				           lm    l-m
	///F_list rotate	      - Spin tensor product 
 

MSVCDLL friend spin_op T_prod(spin_T &SphT, space_T &SphA, int l);
MSVCDLL friend spin_op T_prod(space_T &SphA, spin_T &SphT, int l);

	// Input		SphT  : Irreducible Spin Tensor
	// 			SphA  : Irreducible Spatial Tensor
	//			l     : Tensor Rank Index
	// Output		SOp   : Spin operator which is the result
	//				of the scalar product of  the rank l 
	//			        components of SphT with SphA
	//
	//				  Sum  [     m              ]
	//				  over | (-1)  * T   * A    |
	//				   m   [          lm    l-m ]
	//				or equivalently
	//				  Sum  [     m              ]
	//				  over | (-1)  * A   * T    |
	//				   m   [          lm    l-m ]
 
 
MSVCDLL friend spin_op T_prod(spin_T &SphT, space_T &SphA);
MSVCDLL friend spin_op T_prod(space_T &SphA, spin_T &SphT);

	// Input		SphT  : Irreducible Spin Tensor
	// 			SphA  : Irreducible Spatial Tensor
	// Output		SOp   : Spin operator which is the result
	//				of the scalar product of  the rank l 
	//			        components of SphT with SphA
	//
	//			        Sum   Sum  [     m              ]
	//			        over  over | (-1)  * T   * A    |
	//			         l     m   [          lm    l-m ]
	//				or equivalently
	//			        Sum   Sum  [     m              ]
	//			        over  over | (-1)  * A   * T    |
	//			         l     m   [          lm    l-m ]
 

// ______________________________________________________________________
// I                     AUXILIARY SPIN TENSOR FUNCTIONS
// ______________________________________________________________________


MSVCDLL int Rank();

	// Input		SphT  : Spin tensor (this)
	// Output		r     : Tensor rank
	///F_list Rank                - Spin tensor rank


MSVCDLL friend double Clebsch_Gordan(int a, int b, int alpha, int beta, int c, int gamma);

	// Input		a     : momentum index
	// 			b     : momentum index
	// 			alpha : z-momentum index
	// 			beta  : z-momentum index
	// 			c     : momentum index
	// 			gamma : z-momentum index
	// Output		r     : Clebsch-Gordan coefficients
	// Note 		      : Only Integral (Non-Negative) J
	///F_list Clebsch_Gordan      - Clebsch-Gordan coefficients 


MSVCDLL friend double Wigner_3j(int a, int b, int c, int alpha, int beta, int gamma);

	// Input		a     : momentum index
	// 			b     : momentum index
	// 			c     : momentum index
	// 			alpha : z-momentum index
	// 			beta  : z-momentum index
	// 			gamma : z-momentum index
	// Output		r     : Wigner 3-j coefficients
	///F_list Wigner_3j	      - Wigner 3-j coefficients 
 

// ______________________________________________________________________
//                       SPIN TENSOR I/O FUNCTIONS
// ______________________________________________________________________


MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const spin_T &SphT);

	// Input		ostr : string
	// 			SphT : spin tensor
	// Return		     : stream, prints spin tensor components
	///F_list <<		     - standard output for spin tensors

  };

#endif						// SpinT.h
