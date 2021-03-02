/* LSprop.h *****************************************************-*-c++-*-
**									**
**	                        G A M M A				**
**                                                                      **
**	Liouville Space Propagator 		Interface		**
**                                                                      **
**      Copyright (c) 1998                                              **
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** The class LSprop defines a propagator in spin Liouville space. Such	**
** propagators are simply Liouville space operators which will evolve a	**
** density operator for a specific length of time in a particular set	**
** of rotating frames.							**
**									**
*************************************************************************/

///Chapter Class Liouville Space Propagator 
///Section Overview
///Body    None
///Section Available Spin System Functions

#ifndef   LSprop_h_			// Is file already included?
#  define LSprop_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// this is the interface
#  endif

#include <GamGen.h>			// Include system specifics
#include <LSLib/SuperOp.h>		// Know superoperator info
#include <LSLib/DensOp.h>		// Know density operator class
#include <HSLib/HSprop.h>		// Know Hilbert space propagators


MSVCDLL extern super_op R_prop(super_op& eLt, gen_op sigmaeq);
MSVCDLL extern super_op R_prop(super_op& L, gen_op& sigmaeq, double t);

MSVCDLL extern void evolve_ip(gen_op &sigma, super_op &GOp);
MSVCDLL extern gen_op evolve(const gen_op &sigma, super_op &LOp, const double time);
MSVCDLL extern gen_op evolve(gen_op &sigma, super_op &GOp);



class LSprop
  {
  super_op GOp;				// Propagator superoperator
  double Gt;				// Evolution time
//double *frames;			// Rotating frames

 
private:
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i              CLASS LIOUVILLE SPACE PROPAGATOR ERROR HANDLING
// ____________________________________________________________________________
 
        // Input                G	: LS propagator (this)
        //                      eidx    : Error flag
        //                      noret   : Return flag
        // Output               none    : Error Message Output
 
 
void          LSPerror(int eidx, int noret=0) const;
void volatile LSPfatal(int error) const;
                                                                
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
public:
 
// ____________________________________________________________________________
// A         CLASS LIOUVILLE SPACE PROPAGATOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

/* A NULL propagator is evident because it will have no Liouville space. There
   are propagators (e.g. ideal pulses) that do not evolve the system in time
   as they are infinitely short, so Gt=0 is no indicatation that the propagator
   is zero or identity. An identity propagator is determined IFF the propagator
   matrix is the identity matrix. Again, one must NOT use Ut to determine
   that. This has ramifications in other areas of this class!                */
  
MSVCDLC LSprop();				// Null constructor
MSVCDLC LSprop(int LS);				// Identity constructor


MSVCDLC LSprop(const gen_op& H, double tevol);
MSVCDLC LSprop(const gen_op& H, double tevol, bool prop);
MSVCDLC LSprop(const HSprop& U);

MSVCDLC LSprop(const super_op& L, double tevol);
MSVCDLC LSprop(const super_op& L, const densop& sigma_ss, double tevol);
 
        //                      L       : Active Liouvillian (1/sec)
        //                      sigma_ss: Steady state density operator 
        //                      tevol   : Evolution time (seconds)
        // Output               LSprop  : Propagator for the input Hamiltonian
 
/* The constructor below is used when it is desired to convert a superoperator
   that is already a propagator into a object of type LSprop. For example, this
   can occur if one creates a unitary transform superoperator (U_transform)
   from a general operator (gen_op) that is itself a Hilbert space propagator.
   zum beispeil: gen_op H; gen_op U = prop(H,t); super_op G = U_transform(U);
   Of course, there not many good reasons (perhaps building ideal pulse props)
   to do this because it breaks the ease of using propagators. Instead use code
   such as: gen_op H; HSprop U(H,t); LSprop(U); In the rare times where it is
   necessary to use the constructor below, make sure you specify the 
   propagator evolution time after construction or it will be left as zero.  */

MSVCDLC LSprop(const super_op& G);
  
MSVCDLC         LSprop(const LSprop& G);	// Self-construction
MSVCDLC         ~LSprop();			// Destruction
LSprop& operator= (const LSprop& G1);		// Assignment
 
// ____________________________________________________________________________
// B               LIOUVILLE SPACE PROPAGATOR ACCESS FUNCTIONS
// ____________________________________________________________________________

/*              Function  Output                  Purpose
                --------  ------  ------------------------------------
                  time    double  The length of G in seconds.
                 length   double  The length of G in seconds.
                  dim      int    The Liouville space dimension of G.
                  HS       int    The Hilbert space dimension of G.
                  LS       int    The Liouville space dimension of G.        */

MSVCDLL double   time()   const; 				// Evolve time
MSVCDLL double   length() const;				// Evolve time
MSVCDLL int      dim()    const; 				// LS dimension
MSVCDLL int      HS()     const; 				// HS dimension
MSVCDLL int      LS()     const;				// LS dimension
MSVCDLL super_op LOp()    const;				// Superoperator

MSVCDLL void     L(const super_op& LOp);
MSVCDLL void     length(double t);

// ____________________________________________________________________________
// C                        PROPAGATOR BASIS FUNCTIONS
// ____________________________________________________________________________

MSVCDLL void SetEBR() const;
MSVCDLL void SetBasis(const super_op& LOp);

        // Input                G	: A LS propagator (this)
        //                      LOp     : A superoperator
        // Output               G	: G put into current basis of LOp


 
// ____________________________________________________________________________
// D                   PROPAGATOR EVOLUTION FUNCTIONS 
// ____________________________________________________________________________ 
 
/************** Evolution Under This Superoperator Propagator ****************/

MSVCDLL gen_op evolve(const gen_op& Op);

        // Input                G       : A LS propagator (this)
        //                      Op      : An operator
        // Output               Op1     : Op evolved under prop G

/************ Evolution Under A Static Superoperator (Liouvillian) ***********/ 

MSVCDLL friend gen_op evolve(const gen_op &sigma, super_op &LOp, const double time);

	// Input		sigma : Op to be propagated (dens. mx.)
	//	 		LOp   : Propagation superoperator (rad/sec)
	//	 		time  : Evolution time (seconds)
	// Output		sigma1: Sigma evolved by LOp for time

/************** Evolution Under A Supplied Superoperator Propagator **********/

MSVCDLL friend gen_op evolve(gen_op &sigma, super_op &GOp);

	// Input		sigma : Op. to be propagated (density matrix)
	//	 		GOp   : Super operator propagator
	// Output		Op    : Operator, density matrix propagated
	// Note			      : The superoperator GOp is here assumed
	//				to be a propagation superoperator!


MSVCDLL friend void evolve_ip(gen_op &sigma, super_op &GOp);

	// Input		sigma : Operator propagated (density matrix)
	//	 		GOp   : Super operator propagator
	// Output		none  : Density matrix, sigma, is modified

 
// ____________________________________________________________________________
// C           PROPAGATOR FUNCTIONS, PROPAGATOR WITH PROPAGATOR
// ____________________________________________________________________________
 

 
        // Input                G1	: A LS Propagator.
        //                      G2	: A LS Propagator.
        // Return               G	: LS propagator product of the 2 input
        //                                propagators, G =  G1 * G2.
        // Note				: Order matters - G1*G2 != G2*G1
        // F_list *			- Propagator-Propagator Multiplication

 
 
        // Input                G	: LS propagator (this).
        //                      G1	: LS propagator.
        // Return               Op	: LS Propagator which is the input
        //				  propagator multiplied into G1
        //                                        G = G * G1
        // Note				: Order matters - G1*G2 != G2*G1
        // F_list *			- Prop-Prop Unary Multiplication
 
        // Input                G	: LS propagator (this).
        //                      G1	: LS propagator.
        // Return               G	: LS Propagator which is the input
        //				  propagator multiplied by G1
        //                                        G = G1 * G
        // Note				: Order matters - G1*G2 != G2*G1
        // F_list *			- Prop-Prop Reverse Unary Multiply
 
MSVCDLL LSprop operator *  (const LSprop& G) const;
MSVCDLL LSprop& operator *= (const LSprop& G);
MSVCDLL LSprop& operator &= (const LSprop& G);
 
// ____________________________________________________________________________
// D           PROPAGATOR FUNCTIONS, PROPAGATOR WITH SUPEROPERATOR
// ____________________________________________________________________________
 
 
MSVCDLL friend LSprop operator * (super_op& LOp, LSprop& G);
 
        // Input                LOp     : A superoperator
        //                      G       : Another superoperator
        // Return               G1      : Propagator that is the product of
        //                                the superoperator and propagator,
        //                                      G1 = LOp*G
        // Note                         : Order matters - LOp*G != G*LOp


// ____________________________________________________________________________
// Y               LIOUVILLE SPACE PROPAGATOR OUTPUT FUNCTIONS
// ____________________________________________________________________________
 
  
MSVCDLL friend void set_trace(gen_op& sigma, double tr);
  

MSVCDLL friend super_op R_prop(super_op& eLt, gen_op sigmaeq);

// Input		eLt	: Exponential Liouvillian to relax
//				  the density matrix
//			sigmaeq : Density matrix at equilibrium
//				  (or steady state density matrix)
//			Ho      : Operator in Hilbert space
// Output		LOp     : Superoperator which fullfills
//
//	   eLt * {sigma-sigmaeq} + sigmaeq  = LOp * sigma
//
// Note				: LOp is constructed in the current basis
//				  of sigmaeq


MSVCDLL friend super_op R_prop(super_op& L, gen_op& sigmaeq, double t);

// Input		eLt	: Liouvillian for evolution
//			sigmaeq : Density matrix at equilibrium
//				  (or steady state density matrix)
//			t       : Evolution time
// Output		LOp     : Superoperator which fullfills
//
//	   eLt * {sigma-sigmaeq} + sigmaeq  = LOp * sigma



// ____________________________________________________________________________
// Z               LIOUVILLE SPACE PROPAGATOR OUTPUT FUNCTIONS
// ____________________________________________________________________________
 
        // Input                G       : A LS propagator
        //                      ostr    : An output stream 
	//			full	: Flag for output amount
        // Output               ostr    : Output stream that has had
        //                                G written into it 

MSVCDLL std::ostream& print(std::ostream& ostr, int full=0) const;
MSVCDLL friend std::ostream &operator << (std::ostream &ostr, const LSprop& G);

};

#endif						// LSprop.h
