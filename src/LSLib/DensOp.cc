/* DensOp.cc ****************************************************-*-c++-*-
**									**
**	                        G A M M A				**
**                                                                      **
**	Density Operator			Implementation		**
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
**                                                                      **
** Description                                                          **
**                                                                      **
** The class densop defines a density operator. This operator is by &   **
** large a general operator.  As such it normally resides in Hilbert    **
** space or a composite Hilbert space.  In additon to general operator  **
** activity, density operators track their time since "equilibrium" and **
** "rotating frames".  They also allow themselves to to be evolved by   **
** GAMMA's Hilbert space density operators &/or Liouville space 	**
** propagators.  							**
**                                                                      **
*************************************************************************/

#ifndef _densop_cc_				// Is file already included?
#  define _densop_cc_ 1				// If no, then remember it
#  if defined(GAMPRAGMA)	// Using the GNU compiler?
#    pragma implementation			// this is the implementation
#  endif

#include <LSLib/DensOp.h>			// Include the header
#include <LSLib/LSAux.h>			// Know about LU decomposition
#include <HSLib/GenOp.h>			// Must know operators
#include <LSLib/SuperOp.h>			// Must know superoperators
#include <HSLib/SpinOpCmp.h>			// Must know Fz operator
#include <string> 				// Must know about strings
#include <HSLib/SpinSys.h>			// Know about spin systems
#include <Basics/StringCut.h>			// Know about the Gdec functio
#include <stdlib.h>

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------ 
// ---------------------------------------------------------------------------- 
 
// ____________________________________________________________________________ 
// i                 CLASS DENSITY OPERATOR ERROR HANDLING
// ____________________________________________________________________________
  
  
void densop::SIGMAerror(int eidx, int noret) const
  
        // Input                Sigma	: Density operator (this)
        //                      eidx    : Error flag    
        //                      noret   : Return flag 
        // Output               none    : Error Message Output 
                                                               
  {
  std::cout << "\nClass Density Operator: ";
  switch(eidx)
    {
    case 0:                                                             // (0)
      std::cout << "Program Aborting....";
      break;
    case 1:                                                             // (1)
      std::cout << "Error During Construction";
      break;
    case 10:								// (10)
      std::cout << "Steady State Does Not Exist";
      break;
    default:
      std::cout << "Unknown Error (Number " << eidx << ")";
    }
  if(!noret) std::cout << ".\n";
  else       std::cout << ".";
  }

void volatile densop::SIGMAfatality(int eidx) const
  {                                                                       
  SIGMAerror(eidx);
  if(eidx) SIGMAerror(0);
  std::cout << "\n";
  exit(-1);
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------ 
// ---------------------------------------------------------------------------- 
  
// ____________________________________________________________________________ 
// A            CLASS DENSITY OPERATOR CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________ 
  
 
densop::densop() { Sigmat = 0; }
          
        // Input                none    :
        // Output               Sigma	: A NULL density operator (this) 
        ///F_list               densop  - Constructor 
 

densop::densop(const spin_sys& sys) : gen_op(SigmaEq(sys)) { Sigmat = 0; }

	// Input	       sys : A spin system
	// Output		Op : The density matrix describing
	//			     the system under a high-temp
	//			     approximation.
	// Note			   : The identity matrix is subtracted
	//			     out, leaving the operator traceless
	//			      (normally Tr{density matrix} = 1)
	// Note			   : The density matrix is scaled by a
	//			     constant scaling factor.  This does
	//			     not change the relative intensities
	//			     of any observed transitions

/*
   By definition

               1    [ -H ]   1 [          2    3        ]           -H
     sigmaeq = - exp|----| = - | 1 + X + X  + X  + .... | where X = --
               Z    [ kT ]   Z [                        ]           kT

   When kT >> ||H||, the powers of X become small very quickly. Under such
   circumstances we can drop all powers higher than 1. Furthermore, the
   1 is just along for the ride (under most density operator evolutions)
   so we can just remove it. Lastly, the factors Z & kT are just scaling
   contants which one may leave out since the global scaling is not
   terribly important. These may be removed too:


             ~ 1 [     -H ]     -H
     sigmaeq = - | 1 + -- | --> --- --> -H 
               Z [     kT ]     ZkT

   A few points:  1.) We are usually in the high temperature limit
                  2.) The 1 is needed in some relaxation treatments
                  3.) The value of Z (partition function) is simply the
                      sytsem spin Hilbert space.
                  4.) The factor Z should be used when comparing 
                      magnetization intensities in systems of different
                      Hilbert space dimensions.
	
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


        // Input                sys	: Spin system
	// Input 		L	: Full Liouvillian (rad/sec)
	//			R	: Relaxation Superoperator (rad/sec)
	//			sigmaeq	: Equilibrium density operator
        // Return		sss	: Steady state density matrix
	// Note				: Careful R & L units, usually Ho is Hz!
	// Note				: L is singular, thus steady state must
	//				  be determined in a reduced Liouville
	//				  space.  Best done with matrix routines
	//				  not (code simple) superop. functions!
	// Note				: Algorithm needs density matrix trace = 1
	//		 		  Function insures this then returns trace
	//			 	  back to the original value.

//	L|s  > = R|s  >    -------> [L - L|S><E|]|s  > = R|s  > - |L1>
//         ss       eq		                   ss       eq

densop::densop(const spin_sys& sys, super_op& L, super_op& R) 
       :gen_op(SigmaSS(sys, L, R))     { Sigmat = 0; }
 
densop::densop(super_op& L, super_op& R, gen_op& sigmaeq)
       :gen_op(SigmaSS(L, R, sigmaeq)) { Sigmat = 0; }


densop::densop(gen_op& Op, double tevol) : gen_op(Op) { Sigmat = tevol; }

        //                      Op	: An operator
        //                      tevol	: Evolution time (seconds)
        // Output               densop	: Sigma set to Op at time tevol

                                                                                
densop::densop(const densop& Sigma) : gen_op(Sigma) { Sigmat = Sigma.Sigmat; }
         
        // Input                none    :
        // Output               Sigma	: A NULL density operator (this)


densop::~densop() { }
         
        // Input                Sigma	: A density operator (this)
        // Output               void	: Density operator destructed
 


        // Input                Sigma	: A density operator (this)
        // 			Sigma1	: Another density operator
        // Output               void	: Sigma is set equal to Sigma1

densop & densop::operator= (const densop& Sigma1)
  {
  gen_op::operator=(Sigma1);
  Sigmat = Sigma1.Sigmat;
  return (*this);
  }

// ____________________________________________________________________________
// B                  DENSITY OPERATOR ACCESS FUNCTIONS
// ____________________________________________________________________________


double densop::length() const { return Sigmat; }

        // Input                Sigma	: A density operator (this)
        // Output               Sigmat	: Time of the density operator

// ____________________________________________________________________________
// D                  DENSITY OPERATOR AUXILIARY FUNCTIONS
// ____________________________________________________________________________



        // Input                Sigma	: A density operator (this)
	//			tr	: Value for Sigma trace
	// Output		void    : Sets Sigma trace to be tr

void densop::SetTrace(double tr)
  {
  double acttr = Re(trace(*this));	// Get current trace
  double difftr = acttr - tr;		// Get trace adjustment
  int hs = size();			// Get Sigma Hilbert space
  complex z;				// Temp value
  double fact;				// Temp value
  int i;				// Temp value
  if(difftr)				// Adjust all diagonal 
    {					// elements to fix trace
    fact = difftr/double(hs);
    for(i=0; i<hs; i++)
      {  
      z = get(i,i) - fact;
      put(z, i,i);
      }  
    }  
  }



// ____________________________________________________________________________
// C           PROPAGATOR FUNCTIONS, PROPAGATOR WITH PROPAGATOR
// ____________________________________________________________________________

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


        // Input                sys	: A spin system
        // Output               Sigma	: An operator representing the
	//			          spin system sys at equilibrium

gen_op SigmaEq(const spin_sys& sys)
  {
  gen_op Op;
  if(sys.homonuclear()) Op = Fz(sys);
  else
    {
    double gamma0, gammai;
    gamma0 = sys.gamma(0);
    Op = Iz(sys,0);
    for(int i=1; i<sys.spins(); i++)
      {
      gammai = sys.gamma(i);
      Op += (gammai/gamma0)*Iz(sys,i);
      }
    }
//  spin_op Z = Fz(sys);
//  Op = gen_op(Z.matrix());
//  Op = gen_op(Z.matrix(), sys.get_basis());
  return Op;
  } 


//gen_op sigma_eq(spin_sys& sys, basis& bs)
//return Op;

	// Input		bs : basis
	// Output		Op : high temperature equilibrium
	//			     density matrix in basis bs

//{
//  spin_op Z = Fz(sys);
//  Op = gen_op(Z.matrix(), bs);
//  return;
//} 

// replaced by new routine from Riccardo Balzan balzan.riccardo@parisdescartes.fr
// since this routine give different results if the operators are in different 
// bases. New routine seems to work well.
// 11/2014 MAER

//gen_op SigmaSS(const spin_sys& sys, super_op& L, super_op& R, int warn)
//
//	// Input		sys   : Spin system
//	// 			L     : Full Liouvillian (rad/sec)
//	//			R     : Relaxation Superoperator (rad/sec)
//        // Return		sss   : Steady state density matrix
//	// Note			      : Careful R & L units, usually Ho is Hz!
//	// Note			      : L is singular, thus steady state must
//	//				be determined in a reduced Liouville
//	//				space.  Best done with matrix routines
//	//				not (code simple) superop. functions!
//	// Note			      : Algorithm needs density matrix trace = 1
//	//				Function insures this then returns trace
//	//				back to the original value.
//
////	L|s  > = R|s  >    -------> [L - L|S><E|]|s  > = R|s  > - |L1>
////         ss       eq		                   ss       eq
//
//  {
//  gen_op seq;
//  double d = 1.e-9; 
//  R.LOp_base(L);			// Insure R in Liouville base of L
//  if(R.below(d))			// If R==0, |s  > = 0
//    {					//            ss
//    densop X;
//    if(warn)
//      {
//      densop X;
//      X.SIGMAerror(10);
//      }
//    return seq;
//    }
//  else if(L == R)			// If R==L, |s  > = s  >
//    return SigmaEq(sys);		//	      ss     eq	
//  else
//    {
//    basis Lbs = L.get_basis();		// Store the Hilbert basis of L
//    seq = SigmaEq(sys);			// Equilibrium density matrix
//    int hs = seq.dim();			// Get the Hilbert space size
//    int ls = hs*hs;			// Compute the Liouville space
//    matrix Emx(hs,hs,i_matrix_type);	// Temporary Identity matrix
//    matrix trace1 = Emx/complex(hs);	// Scale by dimension
//    seq += trace1;			// Set trace to 1 for this
//    seq.Op_base(Lbs);			// Put operator into L basis
//
////	Calculate Equation Right Hand Side of Steady State Equation
////
////			  '
////			|s  > = R|s  > - |L1>
////			  eq       eq
//
//    matrix L1, ssp, sub_ssp;
//    L1=(L.get_mx()).get_block(0,0,ls,1);// L1 matrix = 1st column of L
//    ssp=((R*seq).get_mx()).resize(ls,1);// Modified seq, basis of R
//    ssp -= L1;				// |ssp> = R|seq> - |L1>
//    sub_ssp = ssp.get_block(1,0,ls-1,1);// Reduced space seq (as matrix)
//
////	Calculate Equation Left Hand Side of Steady State Equation
////
////			X = L - L|S><E|
////
//
//    matrix Smx(hs, hs, 0.0);		// Temproray Null matrix
//    Smx.put(1.0, 0, 0);			// 1st element of Smx set to 1 
//    gen_op S(Smx, Lbs);			// Create the S operator
//    gen_op E(Emx, Lbs);			// Create the S operator
//    super_op SE(S,E);			// SE = |S><E|
//    super_op X = L;
//    X -= L*SE; 				// X = L - L|S><E|
//    matrix sub_X((X.get_mx()).		// Reduced space X
//	    get_block(1,1,ls-1,ls-1));
//
//// Form the LU Decomposition of the Reduced & Modified Liouville Matrix
//  
///*                      '                    -1  '
//  	     X|s  > = |s  >  ---->  |s  > = X  |s  >
//             - -ss     -eq           -ss    -   -eq
//  
//                       '                       -1  '
//           LU|s  > = |s  >  ---->  |s  > = (LU)  |s  >
//           -- -ss     -eq           -ss     --    -eq
//*/
//
//    matrix ALU = sub_X;
//    int *indx;
//    indx = new int[ls-1];
//    LU_decomp(ALU,indx);
//
///*		Calculate the Steady State Density Matrix
//  
//                         -1  '               trace
//  	     |s  > = (LU)  |s  >  ;   |s  > ------> |s  >
//              -ss     --   -eq         -ss            ss
//*/
//
//
//    matrix sub_sss = sub_ssp;		// Steady state in subspace
//    LU_backsub(ALU,indx,sub_sss);	// Back solve for steady state
//    matrix sssmx(hs, hs, 0.0);		// Reform the full steady state
//    int k=0;				// density matrix in Hilbert
//    for(int i=0; i<hs; i++)		// space.
//      for(int j=0; j<hs; j++)
//        if(i!=0 || j!=0)		// All but first element
//          {
//          sssmx.put(sub_sss.get(k,0),i,j);
//          k++;
//          }
//    sssmx.put(1.0-trace(sssmx),0,0);	// First element from trace
//    sssmx.resize(hs,hs);		// Back to a square array
//    sssmx -= trace1;			// Set trace back to 0
//    gen_op sss(sssmx, seq.get_basis()); // Put this into a gen. op.
//    delete [] indx;
//    return sss;
//    }
//  }
//
//
//gen_op SigmaSS(super_op& L, super_op& R, gen_op& seq, int warn)
//
//	// Input 		L     : Full Liouvillian (rad/sec)
//	//			R     : Relaxation Superoperator (rad/sec)
//	//			seq   : Equilibrium density operator
//        // Return		sss   : Steady state density matrix
//	// Note			      : Careful R & L units, usually Ho is Hz!
//	// Note			      : L is singular, thus steady state must
//	//				be determined in a reduced Liouville
//	//				space.  Best done with matrix routines
//	//				not (code simple) superop. functions!
//	// Note			      : Algorithm needs density matrix trace = 1
//	//				Function insures this then returns trace
//	//				back to the original value.
//
//	L|s  > = R|s  >    -------> [L - L|S><E|]|s  > = R|s  > - |L1>
//         ss       eq		                   ss       eq
//
//  {
//  double d = 1.e-9; 
//  R.LOp_base(L);			// Insure R in Liouville base of L
//  if(R.below(d))			// If R==0, |s  > = 0
//    {					//            ss
//    if(warn)
//      {
//      densop X;
//      X.SIGMAerror(10);
//      }
//    return gen_op();
//    }
//  else if(L == R)			// If R==L, |s  > = s  >
//    return seq;				//	      ss     eq	
//  else
//    {
//    basis Lbs = L.get_basis();		// Store the Hilbert basis of L
//    int hs = seq.dim();			// Get the Hilbert space size
//    int ls = hs*hs;			// Compute the Liouville space
//    matrix Emx(hs,hs,i_matrix_type);	// Temporary Identity matrix
//    matrix trace1 = Emx/complex(hs);	// Scale by dimension
//    seq += trace1;			// Set trace to 1 for this
//    seq.Op_base(Lbs);			// Put operator into L basis
//
////	Calculate Equation Right Hand Side of Steady State Equation
////
////			  '
////			|s  > = R|s  > - |L1>
////			  eq       eq
//
//    matrix L1, ssp, sub_ssp;
//    L1=(L.get_mx()).get_block(0,0,ls,1);// L1 matrix = 1st column of L
//    ssp=((R*seq).get_mx()).resize(ls,1);// Modified seq, basis of R
//    ssp -= L1;				// |ssp> = R|seq> - |L1>
//    sub_ssp = ssp.get_block(1,0,ls-1,1);// Reduced space seq (as matrix)
//
////	Calculate Equation Left Hand Side of Steady State Equation
////
////			X = L - L|S><E|
////
//
//    matrix Smx(hs, hs, 0.0);		// Temproray Null matrix
//    Smx.put(1.0, 0, 0);			// 1st element of Smx set to 1 
//    gen_op S(Smx, Lbs);			// Create the S operator
//    gen_op E(Emx, Lbs);			// Create the S operator
//    super_op SE(S,E);			// SE = |S><E|
//    super_op X = L;
//    X -= L*SE; 				// X = L - L|S><E|
//    matrix sub_X((X.get_mx()).		// Reduced space X
//	    get_block(1,1,ls-1,ls-1));
//
//// Form the LU Decomposition of the Reduced & Modified Liouville Matrix
////
////                      '                    -1  '
//	     X|s  > = |s  >  ---->  |s  > = X  |s  >
//           - -ss     -eq           -ss    -   -eq
//
//                     '                       -1  '
//         LU|s  > = |s  >  ---->  |s  > = (LU)  |s  >
//         -- -ss     -eq           -ss     --    -eq
//
//    matrix ALU = sub_X;
//    int *indx;
//    indx = new int[ls-1];
//    LU_decomp(ALU,indx);
//
////		Calculate the Steady State Density Matrix
////
////                       -1  '               trace
////	     |s  > = (LU)  |s  >  ;   |s  > ------> |s  >
////            -ss     --   -eq         -ss            ss
//
//
//    matrix sub_sss = sub_ssp;		// Steady state in subspace
//    LU_backsub(ALU,indx,sub_sss);	// Back solve for steady state
//    matrix sssmx(hs, hs, 0.0);		// Reform the full steady state
//    int k=0;				// density matrix in Hilbert
//    for(int i=0; i<hs; i++)		// space.
//      for(int j=0; j<hs; j++)
//        if(i!=0 || j!=0)		// All but first element
//          {
//          sssmx.put(sub_sss.get(k,0),i,j);
//          k++;
//          }
//    sssmx.put(1.0-trace(sssmx),0,0);	// First element from trace
//    sssmx.resize(hs,hs);		// Back to a square array
//    sssmx -= trace1;			// Set trace back to 0
//    gen_op sss(sssmx, seq.get_basis()); // Put this into a gen. op.
//    delete [] indx;
//    return sss;
//    }
//  }


gen_op SigmaSS(super_op& L, super_op& R, gen_op& s_eq, int warn){

  // Input 		L     : Full Liouvillian (rad/sec)
  //			R     : Relaxation Superoperator (rad/sec)
  //			seq   : Equilibrium density operator
  // Return		sss   : Steady state density matrix
  // Note		      : Careful R & L units, usually Ho is Hz!
  // Note                     : Direct calculation in Liouvillian eigenbasis

  //Note                      : Solving L|s_ss>=R|s_eq> for |s_ss>

  //Note                      : L is set to eigenbasis
  //Note                      : R is modified to same Hilbert space as L and default Liouvillian space
  //Note                      : s_eq is not modified (copied to seq)
  //Note                      : sss has the same hilbert space as s_eq

  const double dm = 1.e-9;    // Threshold to determine if a Liouvillian matrix is negligible
  const double dv = 1.e-9;    // Threshold to determine if an eigenvalue is null
  const int ls=L.LS();        //Liouville space dimension
  const int hs=L.HS();        //Hilbert space dimension   

//0) Quick dimensional check
  if((R.HS()!=hs)||(s_eq.HS()!=hs)){
    if(warn){				
      densop X;
//    X.SIGMAerror(10);
    }
    return s_eq; 
  }
//CASE 1: negligible relaxation; steady state is null
  if(R.below(dm)) {			// If R==0, |s  > = 0
    if(warn){				//               ss
      densop X;
//    X.SIGMAerror(10);
    }
    return gen_op(); 
  }
//CASE 2: negligible evolution; steady state is equilibrium
  if((L-R).below(dm)){			// If (L-R)==0, |s  > = s  >
    return s_eq;			//	             ss     eq	
  }
//CASE 3: calculation of steady state operators

  gen_op seq=s_eq;                        //copy of equilibrium density matrix

//1) Ensure all of relevant operators in correct basis
  L.set_EBR();                            //Set L in eigenbasis representation
  R.LOp_base(L); R.set_HBR();             //Set R to same Hilbert base of L and then to default Liouvillian representation
  seq.Op_base(L.get_basis());             //Set equilibrium density matrix to same Hilbert space as L
                                          //NOTE: In equations R'|s'_eq>=(B^-1)R|s_eq> so no need to calculate
                                          //|s'_eq> if R is in HBR
//2) Extraction of relevant matrices
  matrix mL=L.get_mx();                           //Diagonal Liouvillian
  matrix mR=R.get_mx();                           //Relaxation term
  matrix mB=L.get_Lbasis().get_mx();              //Transformation to Diagonal Liouvillian basis
  matrix mB_inv(ls,ls,i_matrix_type);mB_inv/=mB;  //Inverse transformation
  matrix vseq=seq.get_mx().resize(ls,1);          //Equilibrium density matrix in vector form
  matrix vI(ls, 1, complex(1,0));                 //Identity vector

//3) Performing calculations
  matrix vR=mB_inv*mR*vseq;                       //Right hand side of equation R'|s'_eq> = (B^-1)R|s_eq>
  matrix vL=mL*vI;                                //vector containing the eigenvalues on the diagonal of L'

//4) Determining vsss' (in L eigenbasis)
  matrix vsss(ls,1);                              //Steady state vector in Liouvillian eigenbasis
  for(int i=0; i<ls; i++){                        //for every eigenvalue
    if(norm(vL.get(i,0))<dv){                     //If it is null
      if(norm(vR.get(i,0))<dv){                   //And if also the right hand side term (R'|s_eq>) is null
	vsss.put(complex(0,0),i,0);               //The element of the vsss vector is null
      }
      else{                                       //Otherwise (null eigenvalue && non-null right hand side) a problem appeared
	if(warn){				  //Throws an error "steady state does not exist"
	  densop X;
//        X.SIGMAerror(10);  
	}                                         //the impossibility to determine the steady state (divergent evolution)
	return gen_op();                          //And returns a null generic_operator
      }
    }
    else{vsss.put(vR.get(i,0)/vL.get(i,0),i,0);}  //Otherwise, if the eigenvalue is non-null then calculates the relative steady state element
  }

//5) Arranging for gen_op output
  matrix msss=(mB*vsss).resize(hs,hs);   //calculatin |s_ss>=B|s'_ss> (default Liouvillian representation) in matrix form
  gen_op sss=seq;                        //creating output operator in same hilbert basis as all the previous calculations
  sss.put_mx(msss);                      //placing matrix in operator
  sss.Op_base(s_eq);                     //Setting output density matrix in input basis
  return sss;                            //output operator
}


gen_op SigmaSS(const spin_sys& sys, super_op& L, super_op& R, int warn)
{  gen_op seq = SigmaEq(sys);
  return SigmaSS(L, R, seq, warn);  
}


// ____________________________________________________________________________
// Y                 DENSITY OPERATOR AUXILIARY FUNCTIONS
// ____________________________________________________________________________

/*

ostream& densop::printMQC(const spin_system& sys, ostream& ostr) const

        // Input                Sigma   : A density operator (this)
	//			ostr	: An output stream
        // Output               ostr	: Output stream that has had
        //              		  Sigma written into it

  {
  int ns = sys.spins();			// Spins in the sytsem
  String labels[ns+1];			// Allocate space for labels
  for(int i=0; i<ns; i++)
    {
    switch(i)
      { 
      case 0:  label[i] = "ZQC"; break;
      case 1:  label[i] = "SQC"; break;
      case 2:  label[i] = "DQC"; break;
      case 3:  label[i] = "TQC"; break;
      default: label[i] = Gdec(i) + String("QC"); break;
      }
    }
  label[ns] = "POP";
  return ostr;
  }
*/


// ____________________________________________________________________________
// Z                 DENSITY OPERATOR OUTPUT FUNCTIONS
// ____________________________________________________________________________


std::ostream& densop::print(std::ostream& ostr, int full) const

        // Input                Sigma   : A density operator (this)
	//			ostr	: An output stream
	//			full	: Flag for amount of output
        // Output               ostr	: Output stream that has had
        //              		  Sigma written into it

  {
  ostr << "\nDensity Operator Time: " << Sigmat << " s\n";
  if(full >= 0) gen_op::print(ostr, full);
  return ostr;
  }


std::ostream &operator << (std::ostream &ostr, const densop &Sigma)

        // Input                Sigma	: A density operator
        //                      ostr    : Output stream
        // Output               none    : Density operator info is sent
        //                                to the output stream

  {
  Sigma.print(ostr);
  return ostr;
  }

#endif						// DensOp.cc
