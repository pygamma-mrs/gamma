/* MultiLib.cc **************************************************-*-c++-*-
**									**
**                                G A M M A 				**
**									**
**      Multiple System Library		           Implementation	**
**									**
**      Copyright (c) 1995						**
**      Nikolai Skrynnikov						**
**      Dr. Scott A. Smith						**
**      National High Magnetic Field Laboratory				**
**      1800 E. Paul Dirac Drive					**
**      Tallahassee Florida, 32306-4005					**
**									**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** This module of function supports the multi_sys, the GAMMA class 	**
** handling mulitple spin systems.  The routines herein	generally	**
** involve such a spin system and build up common operators, in this	**
** case in a direct product space of the systems involved.		**
**									**
*************************************************************************/

#ifndef MultiLib_cc			// Is the file already included?
#define MultiLib_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// Then this is the implementation
#endif

#include <HSLib/SpinSys.h>
#include <MultiSys/MultiLib.h>		// Include the header file
#include <MultiSys/MultiSys.h>		// Include multi_sys spin systems
//#include <MultiSys/MultiAux.h>		// Include multi_sys auxiliary
//#include <HSLib/HSauxil.h>		// Include common sigma_eq
//#include <HSLib/HSham.h>		// Include common Hamiltonians
//#include <HSLib/PulseI.h>		// Include ideal pulses
#include <HSLib/SpinOp.h>		// Include simple spin operators
//#include <HSLib/SpinOpCmp.h>		// Include composite spin operators
//#include <HSLib/SpinOpRot.h>		// Include rotation spin operators
#include <HSLib/GenOp.h>		// Inlcude general operators
//#include <BWRRelax/relaxDip.h>		// Include dipolar relaxation
#include <BWRRelax/relaxCSA.h>		// Include CSA relaxation
#include <BWRRelax/relaxQuad.h>		// Include quadrupole relaxation
#include <BWRRelax/relaxQCSA.h>		// Include CSA-Quad X-correlation
#include <stdlib.h>

using std::vector;			// Using libstdc++ STL vectors
using std::string;			// Using libstdc++ strings

// ____________________________________________________________________________
// A  Generic Functions For Making Primitive Operators in Composite Space
// ____________________________________________________________________________
// sosi - why does this use n_matrix_type?

/*
                   Function Form                             Example
 --------------------------------------------------- --------------------------
 gen_op Name(const spin_sys&)                        Ho(sys)
 gen_op& Op                                          gen_op H = Hcs(sys)

        Input      OpFct   : A function Returning Op (gen_op)
                   Op      : An operator
                   msys    : A multi_sys spin system
                   icomp   : System components affected (-1 = all)
        Output     MOp     : Operator op for multi_sys msys                  */

gen_op multize(gen_op op(const spin_system&), const multi_sys &msys)
  {
  int nc = msys.NComps();			// Number of components
  vector<matrix> mxc;				// Array of matrices
  vector<matrix> bsc;				// Array of bases
  sys_dynamic sys;				// Temp spin system
  gen_op op_sys;				// Temp operator
  for(int i=0; i<nc; i++)			// Loop over msys components
    {
    sys = msys.Comp(i);				//   Get component
    op_sys = op(sys);				//   Operator for component
    mxc.push_back(op_sys.get_mx());		//   Matrix for component
    (mxc[i]).set_type(n_matrix_type);		//   It is a normal matrix
    bsc.push_back((op_sys.get_basis()).U());
    }
  gen_op MultiOp(mxc,bsc);			// Make composiste operator
  return MultiOp;				// Return the operator
  }

gen_op multize(gen_op& Op, const multi_sys &msys, int icomp)
  {
// sosi - check that comp exists
  int nc = msys.NComps();			// Number of components
  vector<matrix> mxc;				// Array of matrices
  vector<matrix> bsc;				// Array of bases
  int hs;
  sys_dynamic sys;				// Temp spin system
  gen_op OpComp;				// Operator component
  for(int i=0; i<nc; i++)			// Loop over msys components
    {
    hs = msys.HS(i);				// Component Hilbert space
    if(i == icomp)				// If its the specified
      {						// component 
      if(hs != Op.dim())			// Insure proper dim
	{
// sosi - need proper error message here 
        exit(-1);
	}
      mxc.push_back(Op.get_mx());		// Store its matrix 
      bsc.push_back((Op.get_basis()).U());	// Store its basis
      }
    else					// If its not the specified
      {						// component
      matrix mxd(hs,hs,0.,d_matrix_type);	//   Temp diagonal mx
      matrix mxi(hs,hs,i_matrix_type);		//   Temp identity mx
      mxc.push_back(mxd);			//   Matrix rep is zero
      bsc.push_back(mxi);			//   Basis is identity
      }
    }
  return gen_op(mxc, bsc);
  }

// ____________________________________________________________________________
// B      General Operator <=== Spin Operator Function(spin system)
// ____________________________________________________________________________

/*                 (Often These Are Rotations Or Spin Operators)

                    Function Form                             Example
 --------------------------------------------------- --------------------------
 spin_op Name(const spin_sys&)                       Fx(sys)
 spin_op Name(const spin_sys&, const string&)        Fx(sys, 1H)
 spin_op Name(const spin_sys&, double)               Rz(sys, 90.0)

        Input      SOpFct  : A function Returning SOp (spin_op)
                   msys    : A multi_sys spin system
                   double  : Anything (often a pulse angle, e.g. 90.0)
                   string  : Anything (often an isotope name, e.g. "1H")
                   icomp   : System components affected (-1 = all)
        Output     MOp     : Operator op for multi_sys msys
                             that parallels SOpFct but resides in
                             the direct product Hilbert space of msys        */

gen_op multize(spin_op op(const spin_sys&), const multi_sys &msys)
  {
  int nc = msys.NComps();			// Number of components
  vector<matrix> mxc;				// Array of matrices
  vector<matrix> bsc;				// Array of bases
  gen_op Op;					// Scratch operator
  for(int i=0; i<nc; i++)			// Loop over spin components
    {
    Op = op(msys.Comp(i));			// Requested Op, component i
    mxc.push_back(Op.get_mx());			// Store its matrix 
    bsc.push_back((Op.get_basis()).U());	// Store its basis
    }
  return gen_op(mxc, bsc);
  }

gen_op multize(spin_op op(const spin_sys&, const string&),
                                     const string& line, const multi_sys &msys)
  {
  int nc = msys.NComps();			// Number of components
  vector<matrix> mxc;				// Array of matrices
  vector<matrix> bsc;				// Array of bases
  gen_op Op;					// Scratch operator
  for(int i=0; i<nc; i++)			// Loop over spin components
    {
    Op = op(msys.Comp(i), line);		// Requested Op, component i
    mxc.push_back(Op.get_mx());			// Store its matrix 
    bsc.push_back((Op.get_basis()).U());	// Store its basis
    }
  return gen_op(mxc, bsc);
  }

gen_op multize(spin_op op(const spin_sys&, double),
                                 double beta, const multi_sys &msys, int icomp)
  {
  int nc = msys.NComps();			// Number of components
  if(icomp != -1) msys.CheckComp(icomp);	// Insure component exists
  vector<matrix> mxc;				// Array of matrices
  vector<matrix> bsc;				// Array of bases
  sys_dynamic sys;				// System for single component
  gen_op op_sys, sigma_sys;
  int subHS=0;
  for(int i=0; i<nc; i++)			// Loop the system components
    {
    sys = msys.Comp(i);			//	Component i (sys dynamic)
//    sigma_sys = sigma.project_sub(i);		//	Density Op this component
    op_sys = gen_op(op(sys, beta));		//	Requested Op this component
    if((icomp == -1) || (i == icomp))		//	Add it if it was selected
      {
      mxc.push_back(op_sys.get_mx());		//	Matrix for operator, comp. i
      bsc.push_back((op_sys.get_basis()).U());	//	Basis for operator, comp. i
      mxc[i].set_type(n_matrix_type);		//	Insure matrix is n_matrix
      bsc[i].set_type(n_matrix_type);		//	Insure basis is n_matrix
      }
    else
      {
      subHS = sys.HS();				//	Hilbert space of component i
      mxc.push_back(matrix(subHS,subHS,i_matrix_type)); // 	Op matrix is I for this component
      bsc.push_back(mxc[i]);			//	Basis is I for this component
      mxc[i].set_type(n_matrix_type);
      bsc[i].set_type(n_matrix_type);
      }
    }
  return gen_op(mxc, bsc);			// Return in composite space
  }

// ____________________________________________________________________________
// C    General Operator <=== General Operator Function(spin_sys, gen_op)
// ____________________________________________________________________________

/*                (Often These Are Evolution Or Pulse Functions)

                    Function Form                             Example
 --------------------------------------------------- --------------------------
 Name(const spin_sys&,gen_op&,double)                Iypuls(sys, sigma, 90.0)
 Name(const spin_sys&,gen_op&,int, double)           Iypuls(sys, sigma, 1, 90.)
 Name(const spin_sys&,gen_op&,const string&, double) Iypuls(sys,sigma,"1H",90.)

        Input      OpFct   : A function Returning Op (gen_op)
                   msys    : A multi_sys spin system
                   icomp   : System components affected (-1 = all)
                   double  : Anything (often a pulse angle, e.g. 90.0)
                   int     : Anyting (often a spin index)
                   string  : Anything (often an isotope name, e.g. "1H")
        Output     MOp     : Operator op for multi_sys msys
                             that parallels OpFct but resides in
                             the direct product Hilbert space of msys        */

gen_op multize(gen_op op(const spin_sys&, const gen_op&, double),
                  gen_op& sigma, double beta, const multi_sys &msys, int icomp)
  {
  int nc = msys.NComps();			// Number of components
  if(icomp != -1) msys.CheckComp(icomp);	// all blocks; otherwise over 
  vector<matrix> mxc;				// Array of matrices
  vector<matrix> bsc;				// Array of bases
  sys_dynamic dsys;
  gen_op op_sys, sigma_sys;
  for(int i=0; i<nc; i++)			// Loop over components
    {
    dsys = msys.Comp(i);			// This is the ith component
    sigma_sys = sigma.project_sub(i);		// sigma in component i subspace
    op_sys = op(dsys, sigma_sys, beta);		// Call desired subspace function
    if((icomp == -1) || (i == icomp))
      {
      mxc.push_back(op_sys.get_mx());
      bsc.push_back((op_sys.get_basis()).U());
      mxc[i].set_type(n_matrix_type);
      bsc[i].set_type(n_matrix_type);
      }
    else
      {
      mxc.push_back(sigma_sys.get_mx());	// remains untouched
      bsc.push_back((sigma_sys.get_basis()).U());
      mxc[i].set_type(n_matrix_type);
      bsc[i].set_type(n_matrix_type);
      }
    }
  return gen_op(mxc, bsc);
  }

gen_op multize(gen_op op(const spin_sys&, const gen_op&, int, double),
       gen_op& sigma, int nspin, double beta, const multi_sys &msys, int icomp)
  {
  int nc = msys.NComps();			// Number of components
  if (icomp != -1)				// -1 to perform operation over
        msys.CheckComp(icomp);			// all blocks; otherwise over
                                                // icomp block
  vector<matrix> mxc;				// Array of matrices
  vector<matrix> bsc;				// Array of bases
  sys_dynamic sys;
  gen_op op_sys, sigma_sys;
  for(int i=0; i<nc; i++)
    {
    sys = msys.Comp(i);
    sigma_sys = sigma.project_sub(i);
    op_sys = op(sys, sigma_sys, nspin, beta);
    if ((icomp == -1) || (i == icomp))
      {
      mxc.push_back(op_sys.get_mx());
      bsc.push_back((op_sys.get_basis()).U());
      mxc[i].set_type(n_matrix_type);
      bsc[i].set_type(n_matrix_type);
      }
    else
      {
      mxc.push_back(sigma_sys.get_mx());	// remains untouched
      bsc.push_back((sigma_sys.get_basis()).U());
      mxc[i].set_type(n_matrix_type);
      bsc[i].set_type(n_matrix_type);
      }
    }
  return gen_op(mxc, bsc);
  }


gen_op multize(gen_op op(const spin_sys&,const gen_op&,const string&, double),
   gen_op& sigma,const string& iso,double beta,const multi_sys &msys,int icomp)
  {
  int nc = msys.NComps();			// Number of components
  if(icomp != -1)				// -1 to perform operation over
        msys.CheckComp(icomp);                 // all blocks; otherwise over
                                                // icomp block
  vector<matrix> mxc;				// Array of matrices
  vector<matrix> bsc;				// Array of bases
  sys_dynamic sys;
  gen_op op_sys, sigma_sys;
  for(int i=0; i<nc; i++)
    {
    sys = msys.Comp(i);
    sigma_sys = sigma.project_sub(i);
    op_sys = op(sys, sigma_sys, iso, beta);
    if ((icomp == -1) || (i == icomp))
      {
      mxc.push_back(op_sys.get_mx());
      bsc.push_back((op_sys.get_basis()).U());
      mxc[i].set_type(n_matrix_type);
      bsc[i].set_type(n_matrix_type);
      }
    else
      {
      mxc.push_back(sigma_sys.get_mx());    // remains untouched
      bsc.push_back((sigma_sys.get_basis()).U());
      mxc[i].set_type(n_matrix_type);
      bsc[i].set_type(n_matrix_type);
      }
    }
  return gen_op(mxc, bsc);
  }

// ____________________________________________________________________________
// D    General Operator <=== General Operator Function(spin_sys, values)
// ____________________________________________________________________________

/*          (Often These Are Propagator or Rotation Functions)

                    Function Form                             Example
 --------------------------------------------------- --------------------------
 Name(const spin_sys&, int, double)                  Iypuls_U(sys, 1, 90.0)
 Name(const spin_sys&, const string&, double)        Iypuls_U(sys, "13C", 90.0)
 Name(const spin_sys&, double)                       Iypuls_U(sys, 90.0)
 Name(const spin_sys&, double, double)               Ixypuls_U(sys, 90.0, 90.0)

        Input      OpFct   : A function Returning Op (gen_op)
                   msys    : A multi_sys spin system
                   icomp   : System components affected (-1 = all)
                   double  : Anything (often a pulse angle, e.g. 90.0)
                   int     : Anyting (often a spin index)
                   string  : Anything (often an isotope name, e.g. "1H")
        Output     MOp     : Operator op for multi_sys msys
                             that parallels OpFct but resides in
                             the direct product Hilbert space of msys        */

gen_op multize(gen_op op(const spin_sys&, int, double),
                    int nspin, double beta, const multi_sys &msys, int icomp)
  {
  int nc = msys.NComps();			// Number of components
  if(icomp != -1) msys.CheckComp(icomp);	// Insure component exists
  vector<matrix> mxc;				// Array of matrices
  vector<matrix> bsc;				// Array of bases
  sys_dynamic sys;				// System for single component
  gen_op op_sys, sigma_sys;
  int subHS=0;
  for(int i=0; i<nc; i++)			// Loop the system components
    {
    sys = msys.Comp(i);			//	Component i (sys dynamic)
//    sigma_sys = sigma.project_sub(i);		//	Density Op this component
    op_sys = op(sys, nspin, beta);		//	Requested Op this component
    if((icomp == -1) || (i == icomp))		//	Add it if it was selected
      {
      mxc.push_back(op_sys.get_mx());		//	Matrix for operator, comp. i
      bsc.push_back((op_sys.get_basis()).U());	//	Basis for operator, comp. i
      mxc[i].set_type(n_matrix_type);		//	Insure matrix is n_matrix
      bsc[i].set_type(n_matrix_type);		//	Insure basis is n_matrix
      }
    else
      {
      subHS = sys.HS();				//	Hilbert space of component i
      mxc.push_back(matrix(subHS,subHS,i_matrix_type)); // 	Op matrix is I for this component
      bsc.push_back(mxc[i]);				//	Basis is I for this component
      mxc[i].set_type(n_matrix_type);
      bsc[i].set_type(n_matrix_type);
      }
    }
  return gen_op(mxc, bsc);			// Return in composite space
  }

gen_op multize(gen_op op(const spin_sys&, const string& iso, double),
          const string& iso, double beta, const multi_sys &msys, int icomp)
  {
  int nc = msys.NComps();			// Number of components
  if(icomp != -1) msys.CheckComp(icomp);	// Insure component exists
  vector<matrix> mxc;				// Array of matrices
  vector<matrix> bsc;				// Array of bases
  sys_dynamic sys;				// System for single component
  gen_op op_sys, sigma_sys;
  int subHS=0;
  for(int i=0; i<nc; i++)			// Loop the system components
    {
    sys = msys.Comp(i);			//	Component i (sys dynamic)
//    sigma_sys = sigma.project_sub(i);		//	Density Op this component
    op_sys = op(sys, iso, beta);		//	Requested Op this component
    if((icomp == -1) || (i == icomp))		//	Add it if it was selected
      {
      mxc.push_back(op_sys.get_mx());			//	Matrix for operator, comp. i
      bsc.push_back((op_sys.get_basis()).U());	//	Basis for operator, comp. i
      mxc[i].set_type(n_matrix_type);		//	Insure matrix is n_matrix
      bsc[i].set_type(n_matrix_type);		//	Insure basis is n_matrix
      }
    else
      {
      subHS = sys.HS();				//	Hilbert space of component i
      mxc.push_back(matrix(subHS,subHS,i_matrix_type)); // 	Op matrix is I for this component
      bsc.push_back(mxc[i]);			//	Basis is I for this component
      mxc[i].set_type(n_matrix_type);
      bsc[i].set_type(n_matrix_type);
      }
    }
  return gen_op(mxc, bsc);			// Return in composite space
  }


gen_op multize(gen_op op(const spin_sys&, double),
                    double beta, const multi_sys &msys, int icomp)
  {
  int nc = msys.NComps();			// Number of components
  if(icomp != -1) msys.CheckComp(icomp);	// Insure component exists
  vector<matrix> mxc;				// Array of matrices
  vector<matrix> bsc;				// Array of bases
  sys_dynamic sys;				// System for single component
  gen_op op_sys, sigma_sys;
  int subHS=0;
  for(int i=0; i<nc; i++)			// Loop the system components
    {
    sys = msys.Comp(i);			//	Component i (sys dynamic)
//    sigma_sys = sigma.project_sub(i);		//	Density Op this component
    op_sys = op(sys, beta);			//	Requested Op this component
    if((icomp == -1) || (i == icomp))		//	Add it if it was selected
      {
      mxc.push_back(op_sys.get_mx());		//	Matrix for operator, comp. i
      bsc.push_back((op_sys.get_basis()).U());	//	Basis for operator, comp. i
      mxc[i].set_type(n_matrix_type);		//	Insure matrix is n_matrix
      bsc[i].set_type(n_matrix_type);		//	Insure basis is n_matrix
      }
    else
      {
      subHS = sys.HS();				//	Hilbert space of component i
      mxc.push_back(matrix(subHS,subHS,i_matrix_type)); // 	Op matrix is I for this component
      bsc.push_back(mxc[i]);			//	Basis is I for this component
      mxc[i].set_type(n_matrix_type);
      bsc[i].set_type(n_matrix_type);
      }
    }
  return gen_op(mxc, bsc);			// Return in composite space
  }

gen_op multize(gen_op op(const spin_sys&, double, double),
               double phi, double beta, const multi_sys &msys, int icomp)
  {
  int nc = msys.NComps();			// Number of components
  if(icomp != -1) msys.CheckComp(icomp);	// Insure component exists
  vector<matrix> mxc;				// Array of matrices
  vector<matrix> bsc;				// Array of bases
  sys_dynamic sys;				// System for single component
  gen_op op_sys, sigma_sys;
  int subHS=0;
  for(int i=0; i<nc; i++)			// Loop the system components
    {
    sys = msys.Comp(i);			//	Component i (sys dynamic)
//    sigma_sys = sigma.project_sub(i);		//	Density Op this component
    op_sys = op(sys, phi, beta);		//	Requested Op this component
    if((icomp == -1) || (i == icomp))		//	Add it if it was selected
      {
      mxc.push_back(op_sys.get_mx());		//	Matrix for operator, comp. i
      bsc.push_back((op_sys.get_basis()).U());	//	Basis for operator, comp. i
      mxc[i].set_type(n_matrix_type);		//	Insure matrix is n_matrix
      bsc[i].set_type(n_matrix_type);		//	Insure basis is n_matrix
      }
    else
      {
      subHS = sys.HS();				//	Hilbert space of component i
      mxc.push_back(matrix(subHS,subHS,i_matrix_type)); // 	Op matrix is I for this component
      bsc.push_back(mxc[i]);				//	Basis is I for this component
      mxc[i].set_type(n_matrix_type);
      bsc[i].set_type(n_matrix_type);
      }
    }
  return gen_op(mxc, bsc);			// Return in composite space
  }

// sosik






// ----------------------------------------------------------------------------
//	          Generic Functions Which Use Superoperators
//  This one takes a single superoperator of one component and expands
//  the superoperator matrix to the composite Liouville space dimensions
// ----------------------------------------------------------------------------



super_op multize(super_op& SOp, const multi_sys &msys, int icomp)
  {
// sosi - check that comp exists
  int nc = msys.NComps();			// Number of components
  vector<matrix> mxc;				// Array of matrices
  vector<matrix> bsc;				// Array of bases
  int ls, hs;
  sys_dynamic sys;				// Temp spin system
  gen_op OpComp;				// Operator component
  for(int i=0; i<nc; i++)			// Loop over msys components
    {
    ls = msys.LS(i);				// Component Liouville space
    hs = msys.HS(i);				// Component Hilbert space
    if(i == icomp)				// If its the specified
      {						// component 
      if(ls != SOp.dim())			// Insure proper dim
	{
// sosi - need proper error message here 
        exit(-1);
	}
      mxc.push_back(SOp.get_mx());		// Store its matrix 
      bsc.push_back((SOp.get_basis()).U());	// Store its basis
      }
    else					// If its not the specified
      {						// component
      matrix mxd(ls,ls,0.,d_matrix_type);	//   Temp diagonal mx
      matrix mxi(hs,hs,i_matrix_type);		//   Temp identity mx
      mxc.push_back(mxd);			//   Matrix rep is zero
      bsc.push_back(mxi);			//   Basis is identity
      }
    }
  return super_op(mxc, bsc);
  }



// ----------------------------------------------------------------------------
//	          Generic Functions Which Use Superoperators
//            (Often These Are Relaxation Or Exchange Functions)
// ----------------------------------------------------------------------------

super_op multize(super_op SOpFct(const gen_op&), const gen_op& Op, const multi_sys& msys)
// sosi - basis problems with this function?!

        // Input        SOpFct	: Function which returns a  superoperator
        //                        given an operator as its only argument
        //              Op	: An operator
        //              msys    : A multi_sys spin system 
        // Output       MSOp    : Superoperator SOp for multi_sys msys 
        //                        in the direct product Liouville space 
        // Note                 : This routine calls an arbitrary function 
        //                        of the following form 
        //                              super_op SOpFct(gen_op& Op) 
        //                        where Op is the only function arguments
	//			  and SopFct is the name of the function
	// Note			: There is no check herein to insure that
	//			  the composite space of Op is the same
	//			  as that which is the default of msys

  {
  int nc = msys.NComps();			// Number of components
  matrix *mxc;				// Array for superop matrices
  mxc = new matrix[nc];
  gen_op OpComp;				// Operator component
  super_op SOpComp;				// Superoperator component
  for(int i=0; i<nc; i++)			// Loop all msys components
    {
    OpComp = Op.project_sub(i);			// Get Op's ith component
    SOpComp = SOpFct(OpComp);			// Superoperator ith component
    mxc[i] = SOpComp.get_mx();			// Store SOpComp matrix rep.
    }
  super_op LOp(mxc, nc);			// implies default basis
  delete [] mxc;
  return LOp;
  }


super_op multize(super_op sop(const sys_dynamic&,gen_op&,int,int), gen_op& H,
                                     int type, int level, const multi_sys& msys)
// sosi - basis problems with this function?!   This is
//        the main relaxation superoperator interface.

        // Input        sop     : Function which returns a  
        //                        superoperator given a system, operator, 
        //                        and two integers
        //              H       : An operator (isotropic Hamiltonian) 
	//		type    : Computation type  
	//		level   : Computation level 
        //              msys    : A multi_sys spin system 
        // Output       MSOp    : Superoperator op for multi_sys msys 
        //                        in the direct product Liouville space 
        // Note                 : This routine calls an arbitrary function 
        //                        of the following form 
        //                        super_op Name(sys_dynamic&, gen_op&, int, int) 
        //                        where H, type and level will be the arguments    

  {
  int nc = msys.NComps();			// Number of components
  vector<matrix> mxc;				// Array of superop matrices
  vector<matrix> bsc;				// Array of superop HS bases
  sys_dynamic dsys;				// Temp dynamic spin system
  gen_op H_sys;					// Temp isotropic Hamiltonian
  super_op R;					// Temp (relaxation) superop
  for(int i=0; i<nc; i++)			// Loop all components (systems)
    {
    dsys = msys.Comp(i);			// Get component i
    H_sys = H.project_sub(i);			// H's ith comp. (i's H.space)
    R = sop(dsys, H_sys, type, level);		// Rfunct comp. i (i's L.space)
    mxc.push_back(R.get_mx());			// Store matrix rep this superop
    bsc.push_back((R.get_basis()).U());		// Store HS basis for this rep
    }
  return super_op(mxc, bsc);			// Implies default basis
  }


super_op multize(super_op sop(const sys_dynamic&,gen_op&,int), gen_op& H,
                                          int level, const multi_sys& msys)

        // Input        sop     : Function which returns an 
        //                        superoperator given a system, operator, 
        //                        and two integers
        //              H       : An operator (isotropic Hamiltonian) 
	//		level   : Computation level 
        //              msys    : A multi_sys spin system 
        // Output       MSOp    : Superoperator op for multi_sys msys 
        //                        in the direct product Liouville space 
        // Note                 : This routine calls an arbitrary function 
        //                        of the following form 
        //                        super_op Name(sys_dynamic&, gen_op&, int, int)
        //                        where H, type and level will be the arguments    

// sosi - basis problems with this function?!
  {
  int nc = msys.NComps();			// Number of components
  matrix *mxc;				// Array for superop matrices
  mxc = new matrix[nc];
  sys_dynamic dsys;				// Temp dynamic spin system
  gen_op H_sys;					// Temp isotropic Hamiltonian
  super_op R;					// Temp (relaxation) superop.
  for(int i=0; i<nc; i++)			// Loop all components (systems)
    {
    dsys = msys.Comp(i);			// Get component i
    H_sys = H.project_sub(i);			// H's ith comp. (i's H. space)
    R = sop(dsys, H_sys, level);		// R's ith comp. (i's L. space)
    mxc[i] = R.get_mx();			// Store matrix rep of superop
    }
  super_op LOp(mxc, nc);			// implies default basis
  delete [] mxc;
  return LOp;
  }

	// Input	msys	: A multi_sys spin system
	// Output	bs      : A default basis in the
	//			  composite Hilbert space
	//			  of msys

basis D_basis(const multi_sys& msys) { return basis(msys.HSs()); } 

// ____________________________________________________________________________
// E            Composite Liouville Space SuperBras And SuperKets
// ____________________________________________________________________________

/* A basis function, i, in Hilbert space corresponds to a set of spin quantum
   numbers, mz, for the spins involved in the system. For example, a single
   spin system of 3 spin 1/2 particles will have the 8 basis functions

  		        |+++>, |++->, |+-+>, ... |---> 
                          0      1      2          7 

   where + = mz of 1/2 and - = mz of -1/2. In a multiple spin system the same
   situation exists except that the # of basis functions is the sum of the
   individual spin systems basis function count, and the basis functions occur
   in blocks because there is no blending between the systems. For example, if
   we had a system of three components A, B, and C which had spin Hilbert 
   spaces of 6, 8, and 4 respectively we could write our basis functions as

        |A0>, |A1>, ... |A5>, |B0>, |B1>, ... |B7>, |C0>, |C1>, ... |C7>
         0     1         5     6     7         9     10    11        17

   Within each block the spin quantum values of the other components are not
   of consequence and may be ignored.

   We can (and need to) do the same thing in Liouville space. Again the basis
   functions correspond to individual spin quantum numbers but each is based on
   both a <bra| and |ket> of the Hilbert space basis functions. For our single
   system of 3 spin 1/2 particles we have

 |+++><+++|, |+++><++-|, ..., |+++><---|, |++-><+++|, |++-><++-|,... |---><---|
     0           1                7            8          9              63

   Note that in this basis product the bra is that changing indices fastest.
   And again, when we form a multipe spin system that is a composite of other
   systems the total number of basis functions is the addition of their 
   respective basis funtion.  For the system of three components A, B, and C 
   which had spin Hilbert spaces of 6, 8, and 4 respectively we get Liouville
   spaces of 36, 64, and 16 which sums to a total dimension of 116. We can use
   the simpler nomenclature to write them as

 |A0><A0|,|A0><A1|,...|A6><A6|,|B0><B0|,|B0><B1|,...|B7><B7|,|C0><C1|,...|C3><C3|
     0       1           35       36       37           99      100        115

   Notice that there is correspondance between the basis function index and the
   block (or component) with which the function is associated. 

   The functions supplied here are quite simple & based on the above reasoning.
   Given a basis function index these return the spin states for all spins in
   the associated component. For example consider again the multiple spin 
   system of ABC. In Hilbert space basis, function 9 has spin states associated
   with system B in its 7th basis function. In Liouville state the 37th function
   has spin states associated with the bra of system B in its 1st state and with
   the ket of system B in its 2nd state. Note that there is never any mixing of
   different components so their spin states are ignored.

		Input	I	: Basis function index (full LS)
  			msys	: Multiple spin system
  		Output	mk	: Bra/Ket of basis state I
  		Note		: For these functions basis ordering
  				  is as      |+++>, |++->, |+-+>, ... 
  				  instead of |+++>, |-++>, |+-+>, ...
  				  so the ordering is reversed!               */

/* ----------------------------------------------------------------------------
   Here we desire an array of spin mz values for Liouville space basis function
   I. This corresponds to 1 and only 1 component spin system so we simply need
   the that component's spin mz values for Liouville state I. There are 2 such
   arrays, one associated with the bra and one with the ket. For the bra array
   the Hilbert space functions just cycle through. That is, Liouville index I
   corresponds to index Io%HS in Hilbert space where HS is the spin Hilbert 
   of the component involved and Io is the index I relative to the start of 
   the block containing the components basis functions. The ket array has its
   Hilbert space indexing jumping around so that for Liouville index I we need
   Hilbert space index Io - HS*(Io%HS).

   Lets be clear on this by looking at an example. In the ABC system we have
   been using, suppose one wanted the values for basis function 98. That will
   correspond to component B which has its block beginning at 36 and a H.S.
   dimension of 8. The H.S.  basis function index used for the bra is simply
   (98-36)%8 = (62%8) = 6. Similarly for the ket we would use H.S. index of
   62 - 8*(62%8) = 62-8*6 = 62-56 = 7; Both of these numbers were indicated
   earlier for basis function 98 
                                          |B7><B6|
                                             98
*/


row_vector LS_qState_bra(const multi_sys& msys, int I)
  {
  basis Dbs          = D_basis(msys);		// Get default LS basis of msys
  int cmp            = Dbs.which_sub_LS(I);	// Component hosting state I
  int i              = I-Dbs.sub_anchor_LS(cmp);// Local state I index 
  sys_dynamic sys    = msys.Comp(cmp);		// Component for state i
  int         dim    = sys.HS();		// Component Hilbert space 
  int         nspins = sys.spins();		// No. spins in the component
  int         j      = i/dim;			// Mod. div., j bra in |j><k|
  row_vector mk1 = sys.qState(j);		// jth state spin mz values
  row_vector mk(nspins);			// Vector for LS basis ket
  for(int it=0; it<nspins; it++)		// Just loop spins (elements)
    mk.put(mk1.get(nspins-it-1),it); 		// and reverse the ordering
  return mk;					// This is the state vector
  } 

vector<double> qStateLS(const multi_sys& msys, int I)
  {
  basis Dbs          = D_basis(msys);		// Get default LS basis of msys
  int   ncmp	     = msys.NComps();		// Number of components
  int i, j, dim, nspins;
  sys_dynamic sys;
  vector<double> mk;
  row_vector mk1;
  for(int cmp=0; cmp<ncmp; cmp++)		// Loop system components
    {
    sys    = msys.Comp(cmp);			// Component for state I
    dim    = sys.HS();				// Component Hilbert space 
    i      = I - Dbs.sub_anchor_LS(cmp);	// Local state I index 
    nspins = sys.spins();			// No. spins in the component
    j      = i/dim;				// Mod. div., j bra in |j><k|
    mk1    = sys.qState(j);			// jth state spin mz values
    for(int it=0; it<nspins; it++)		// Just loop spins (elements)
      mk.push_back(mk1.getRe(nspins-it-1)); 	// and reverse the ordering
    }
  return mk;					// This is the state vector
  } 

row_vector LS_qState_ket(const multi_sys& msys, int I)
  {
  basis Dbs          = D_basis(msys);		// Get default LS basis of msys
  int cmp            = Dbs.which_sub_LS(I);	// Component hosting state I
  int i              = I-Dbs.sub_anchor_LS(cmp);// Local state I index 
  sys_dynamic sys    = msys.Comp(cmp);		// Component for state I/i
  int         dim    = sys.HS();		// Component Hilbert space 
  int         nspins = sys.spins();		// No. spins in the component
  int k = i - (dim*(i/dim));			// Mod. div., k ket in |j><k|
  row_vector mk1 = sys.qState(k);		// kth state spin mz values 
  row_vector mk(nspins);			// Vector for LS basis ket
  for (int it=0; it<nspins; it++)		// Just loop spins (elements)
    mk.put(mk1.get(nspins-it-1),it); 		// and reverse the ordering
  return mk;
  }

#endif							// MultiLib.cc
