/* DensOp.h *****************************************************-*-c++-*-
**									**
**	                        G A M M A				**
**                                                                      **
**	Density Operator			Interface		**
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
** The class densop defines a density operator. This operator is by &	**
** large a general operator.  As such it normally resides in Hilbert	**
** space or a composite Hilbert space.  In additon to general operator	**
** activity, density operators track their time since "equilibrium" and	**
** "rotating frames".  They also allow themselves to to be evolved by	**
** GAMMA's Hilbert space density operators &/or Liouville space		**
** propagators.								**
**									**
*************************************************************************/

///Chapter Class Density Operator
///Section Overview
///Body    None
///Section Available Density Operator Functions

#ifndef   Densop_h_			// Is file already included?
#  define Densop_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// this is the interface
#  endif

#include <GamGen.h>			// Include system specifics
#include <HSLib/GenOp.h>		// Know operator information
#include <LSLib/SuperOp.h>		// Know about superoperator class
#include <HSLib/SpinSys.h>		// Know about spin system class

//forward declarations
class densop;

MSVCDLL gen_op SigmaEq(const spin_sys& sys);
MSVCDLL gen_op SigmaSS(const spin_sys& sys, super_op& L, super_op& R, int wrn=0);
MSVCDLL gen_op SigmaSS(super_op& L, super_op& R, gen_op& seq, int wrn=0);

class densop : public gen_op
  {
  double Sigmat;			// Evolution time
//double *frames;			// Rotating frames

 
private:
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// i              CLASS DENSITY OPERATOR ERROR HANDLING
// ____________________________________________________________________________
 
 
void SIGMAerror(int eidx, int noret=0) const;
 
        // Input                U	: Density operator (this)
        //                      eidx    : Error flag
        //                      noret   : Return flag
        // Output               none    : Error Message Output
 
 
void volatile SIGMAfatality(int error) const;
 
        // Input                U	: Density operator (this)
        //                      eidx    : Error flag
        // Output               none    : Stops execution & error Message

                                                                
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
public:
 
// ____________________________________________________________________________
// A         CLASS DENSITY OPERATOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 
  
MSVCDLC densop();
         
        // Input                none    :
        // Output               U       : A NULL density operator (this)
        ///F_list               densop  - Constructor

                                                      
MSVCDLC densop(const spin_sys& sys);
 
        // Input               sys	: A spin system
        // Output              Sigma	: The density matrix describing
        //				  the system under a high-temp
        //				  approximation.
        // Note				: The identity matrix is subtracted
        //				  out, leaving the operator traceless
        //				  (normally Tr{density matrix} = 1)
        // Note				: The density matrix is scaled by a
        //				  constant scaling factor.  This does
        //				  not change the relative intensities
        //				  of any observed transitions


MSVCDLC densop(const spin_sys& sys, super_op& L, super_op& R);

	// Input		sys   : Spin system
	// 			L     : Full Liouvillian
	//			R     : Relaxation/Exchange Superoperator
        // Return		Sigma : Steady state density operator
	// Note			      : Both R and L must be in the same
	//				units (rad/sec).  Be careful because
	//		 		usually Ho is generated in Hz!
	// Note			      : Because L is singular, the steady
	//				state density matrix must be
	//				determined via a reduced Liouville
	//				space.  In turn, this function is
	//				forced to use matrix routines rather
	//				than the (code simple) superoperator
	//				functions!
	// Note			      : Because this algorithm assumes the
	//				trace of all density matrices is 1
	//				we must first set the equilibrium matrix
 

MSVCDLC densop(super_op& L, super_op& R, gen_op& sigmaeq);
 
        // Input                L       : Full Liouvillian (rad/sec)
        //                      R       : Relaxation Superoperator (rad/sec)
        //                      sigmaeq : Equilibrium density operator
        // Return               sss     : Steady state density matrix
        // Note                         : Careful R & L units, usually Ho is Hz!
        // Note                         : L is singular, thus steady state must
        //                                be determined in a reduced Liouville
        //                                space.  Best done with matrix routines
        //                                not (code simple) superop. functions!
        // Note                         : Algorithm needs density matrix trace = 1
        //                              Function insures this then returns trace
        //                              back to the original value.
         
//      L|s  > = R|s  >    -------> [L - L|S><E|]|s  > = R|s  > - |L1>
//         ss       eq                             ss       eq


MSVCDLC densop(gen_op& Op, double tevol);

        //                      Op      : Active density operator
        //                      tevol   : Evolution time (seconds)
        // Output               densop  : Sigma set to Op at time tevol

  
MSVCDLC densop(const densop& Sigma);
         
        // Input                none    :
        // Output               Sigma	: A NULL density operator (this)
        ///F_list               densop  - Constructor

 
MSVCDLC ~densop();
         
        // Input                Sigma	: A density operator
        // Output               void    : Density operator destructed
 
MSVCDLL densop & operator= (const densop& Sigma1);
 
        // Input                Sigma       : A density operator (this)
        //                      Sigma1      : Another density operator
        // Output               void    : Sigma is set equal to Sigma1
 



// ____________________________________________________________________________
// B               DENSITY OPERATOR ACCESS FUNCTIONS
// ____________________________________________________________________________


MSVCDLL double length() const;

        // Input                Sigma	: A density operator (this)
        // Output               t	: Time the density operator spans

 
// ____________________________________________________________________________
// C           PROPAGATOR FUNCTIONS, PROPAGATOR WITH PROPAGATOR
// ____________________________________________________________________________
 
// ____________________________________________________________________________
// D                  DENSITY OPERATOR AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 
 
MSVCDLL void SetTrace(double tr);
 
        // Input                Sigma   : A density operator (this)
        //                      tr      : Value for Sigma trace 
        // Output               void    : Sets Sigma trace to be tr
 

// ____________________________________________________________________________
// X                 DENSITY OPERATOR FRIEND FUNCTIONS
// ____________________________________________________________________________

/*
   By definition

               1    [ -H ]   1 [          2    3        ]           -H
     sigmaeq = - exp|----] = - | 1 + X + X  + X  + .... | where X = --
               Z    [ kT ]   Z [                        ]           kT
 
   When kT >> -H, the identity matrix is neglected, and the factor ZkT removed,
 
             ~ 1 [     -H ]     -H
     sigmaeq = - | 1 + -- | --> --- --> -H
               Z [     kT ]     ZkT
 
   For the isotropic static Hamiltonian in NMR, neglecting chemical shifts and
   coupling constants (e.g. energy level populations are insignificantly
   affected by these contributions, so H is the Zeeman Hamiltonian)
 
                 ---                           ---  gamma
             ~   \                             \         i     homonuclear
     sigmaeq = - /   - hbar * gamma  * Iz  --> /    ------ Iz  -----------> Fz
                 ---               i     i     ---  gamma    i
                  i                             i        0
 
   where the operator has be again rescaled by hbar*gamma
                                                         0                   */  
 
MSVCDLL friend gen_op SigmaEq(const spin_sys& sys);
 
        // Input                sys     : A spin system
        // Output               Sigma   : An operator representing the
        //                                spin system sys at equilibrium
 



	// Input		sys   : Spin system
	// 			L     : Full Liouvillian
	//			R     : Relaxation Superoperator
        //                      seq   : Equilibrium density operator
        // Return		sss   : Steady state density matrix
	// Note			      : Both R and L must be in the same
	//				units (rad/sec).  Be careful because
	//		 		usually Ho is generated in Hz!
	// Note			      : Because L is singular, the steady
	//				state density matrix must be
	//				determined via a reduced Liouville
	//				space.  In turn, this function is
	//				forced to use matrix routines rather
	//				than the (code simple) superoperator
	//				functions!
	// Note			      : Because this algorithm assumes the
	//				trace of all density matrices is 1
	//				we must first set the equilibrium matrix
	//				to satisfy this condition them reset the
	//				final answer back to trace = 0 as is used
	//				the rest of the GAMMA routines
        // Note                       : L is singular, thus steady state must
        //                              be determined in a reduced Liouville
        //                              space.  Best done with matrix routines
        //                              not (code simple) superop. functions!
        // Note                       : Algorithm needs density matrix trace = 1
        //                              Function insures this then returns trace
        //                              back to the original value.

//	L|s  > = R|s  >    -------> [L - L|S><E|]|s  > = R|s  > - |L1>
//         ss       eq		                   ss       eq


MSVCDLL friend gen_op SigmaSS(const spin_sys& sys, super_op& L, super_op& R, int wrn);
MSVCDLL friend gen_op SigmaSS(super_op& L, super_op& R, gen_op& seq,         int wrn);


// ____________________________________________________________________________
// Z               DENSITY OPERATOR OUTPUT FUNCTIONS
// ____________________________________________________________________________
 
        // Input                Sigma	: A density operator
        //                      ostr    : An output stream 
	//			full	: Flag for output amount
        // Output               ostr    : Output stream that has had
        //                                Sigma written into it 

MSVCDLL std::ostream& print(std::ostream& ostr, int full=0) const;
MSVCDLL friend std::ostream &operator << (std::ostream &ostr, const densop& Sigma);

};
 
#endif							// densop.h
