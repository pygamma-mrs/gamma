/* MultiLOp.cc **************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Multiple System SuperOperators             Implementation 	**
**                                                                      **
**      Copyright (c) 1995                                              **
**      Nikolai Skrynnikov                                              **
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header:
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** This file of functions support multi_sys, the GAMMA class for        **
** handling mulitple spin systems.  The routines herein generally       **
** involve such a spin system and build up common superoperators.       **
** The functions mirror those defined in the Liouville space library 	**
** (see module LSLib) but will return operators that exist in the       **
** composite Liouville spin space spanning the components in multi_sys. **
** The Liouville space is a direct product of the systems involved. 	**
**                                                                      **
*************************************************************************/

#ifndef   MultiLOp_cc			// Is the file already included?
#  define MultiLOp_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <MultiSys/MultiLOp.h>		// Include our header file
#include <MultiSys/MultiLib.h>		// Include the multize functions
#include <MultiSys/MultiSys.h>		// Include multi_sys spin systems
#include <HSLib/GenOp.h>		// Include general operators
#include <HSLib/HSham.h>		// Include H.S. Hamiltonians
#include <LSLib/SuperOp.h>		// Include general superoperators
#include <LSLib/LSAux.h>		// Include function Hsuper
#include <HSLib/SpinSys.h>		// Include base spin ssytems


        // Input		Heff	: Effective Hamiltonian (Hz)
        // Output               LOp	: Hamiltonian commutation superoperator
	// Note				: This superoperator is returned
	//				  in the default basis
        // Note				: LOp is returned in angular frequency
        //                                units
	// Note				: The multize function below must
	//				  be defined for this to work

super_op Hsuper(const multi_sys& msys, const gen_op& Heff)
  { return multize(Hsuper,Heff,msys); }


	// Input		msys    : A multi_spin spin system
	// Output		Lo	: The isotropic Hamiltonian
	//				  commutation superoperator.
	// Note				: This superoperator is returned
	//				  in the default basis


super_op Lo(const multi_sys &msys)
  {
  int nc = msys.NComps();		// Number of components
  matrix *mxc;				// Array for superop matrices
  mxc = new matrix[nc];
  sys_dynamic sys;			// Scratch spin system
  gen_op Op;				// Scratch operator
  super_op LOp;				// Scratch superoperator
  for(int i=0; i<nc; i++)
    {
    sys = msys.Comp(i);			// The ith spin system
    Op = Ho(sys);			// Isotropic H for ith comp
    Op.set_DBR();			// Insure its in a default basis
    LOp = commutator(Op);		// Commutation superoperator of Op
    mxc[i] = LOp.get_mx();
    (mxc[i]).set_type(n_matrix_type);
    mxc[i] *= ((-complex1)*2*PI);
    }
  super_op LLOp(mxc, nc);			// implies default basis
  delete [] mxc;
  return LLOp;
  }



	// Input	H 	: A operator
	// Output	LOp     : A superoperator which will
	//			  diagonalize H


super_op U_LS(gen_op& H)
  {
  basis bs = H.get_basis();			// Get H's basis 
  int nc = bs.sub_N();				// Number of subspaces
  matrix *mxc;				// Storage for subspace arrays
  mxc = new matrix[nc];
  gen_op Op, op_U;
  matrix mx_sys;
  super_op sop_U; 
  for(int i=0; i<nc; i++)			// Loop over all components
    {
    Op = H.project_sub(i);			// H's ith comp. (ith subspace)
    Op.set_EBR();				// Set H's ith comp. in eigenbasis
    mx_sys = (Op.get_basis()).U();
//    mx_sys = (const matrix&)Op.get_basis();
    op_U = gen_op(mx_sys);  
    sop_U = U_transform(op_U);			// That's U x U
    mxc[i] = sop_U.get_mx();
    }
  super_op LOp(mxc, nc);			// implies default basis
  delete [] mxc;
  return LOp;
  }


super_op Uinv_LS(gen_op& H)

	// Input	H 	: A operator
	// Output	LOp     : Operator op for multi_sys msys

  {
  basis bs = H.get_basis();			// Get H's basis 
  int nc = bs.sub_N();				// Number of subspaces
  matrix *mxc;				// Storage for subspace arrays
  mxc = new matrix[nc];
  gen_op op_sys, op_U;
  matrix mx_sys;
  super_op sop_U;
  for (int i=0; i<nc; i++)
    {
    op_sys = H.project_sub(i);
    op_sys.set_EBR();                         // basis gives diagonalizing U
    mx_sys = (op_sys.get_basis()).U();
//    mx_sys = (const matrix&) op_sys.get_basis();
    mx_sys = inv(mx_sys);			// this is U(-1)
    op_U = gen_op(mx_sys);
    sop_U = U_transform(op_U);                    // That's U(-1) x U(-1)
    mxc[i] = sop_U.get_mx();
    }
  super_op LOp(mxc, nc);			// implies default basis
  delete [] mxc;
  return LOp;
}


super_op Op_Ebase(super_op& L, gen_op& H)
  {
  super_op Uh = U_LS(H);
  super_op Uh_1 = Uinv_LS(H);
// sosi - don't know about this one
//  super_op uLu = Uh_1*L*Uh; 				// Put into eigenbasis of H
  super_op uLu = Uh_1; 				// Put into eigenbasis of H
  uLu *= L;
  uLu *= Uh;
  return uLu;
  } 

#endif							// MultiLOp.cc
