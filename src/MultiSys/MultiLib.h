/* MultiLib.h ***************************************************-*-c++-*-
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      Multiple System Library		               Interface	**
**                                                                      **
**      Copyright (c) 1995                                              **
**      Nikolai Skrynnikov                                              **
**      Dr. Scott A. Smith                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** This module of function supports the multi_sys, the GAMMA class      **
** handling mulitple spin systems.  The routines herein generally       **
** involve such a spin system and build up common operators, in this    **
** case in a direct product space of the systems involved.              **
**                                                                      **
*************************************************************************/

#ifndef   MultiLib_h			// Is the file already included?
#  define MultiLib_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/GenOp.h>		// Knowledge of operators
#include <LSLib/SuperOp.h>		// Knowledge of superoperators
#include <Matrix/row_vector.h>		// Knowledge of row vectors
#include <MultiSys/MultiSys.h>		// Know multiple spin systems
#include <HSLib/SpinSystem.h>		// Know isotropic systems
#include <LSLib/sys_dynamic.h>		// Know anisotropic systems
#include <vector>			// Know libstdc++ vectors

 
// ____________________________________________________________________________
// A  Generic Functions For Making Primitive Operators in Composite Space
// ____________________________________________________________________________
/*
                    Function Form                             Example
 --------------------------------------------------- --------------------------
 gen_op Name(const spin_sys&)                        Ho(sys)
 gen_op& Op                                          gen_op H = Hcs(sys)

	Input	   OpFct   : A function Returning Op (gen_op)
                   Op      : An operator
                   msys    : A multi_sys spin system
        	   icomp   : System components affected (-1 = all)
	Output     MOp     : Operator op for multi_sys msys                  */

MSVCDLL gen_op multize(gen_op  OpFct(const spin_system &sys), const multi_sys &msys); 
MSVCDLL gen_op multize(gen_op& Op, const multi_sys &msys, int icomp);

// ____________________________________________________________________________
// B      General Operator <=== Spin Operator Function(spin system)
// ____________________________________________________________________________

/*                 (Often These Are Rotations Or Spin Operators)

                    Function Form                             Example
 --------------------------------------------------- --------------------------
 spin_op Name(const spin_sys&)                       Fx(sys)
 spin_op Name(const spin_sys&, const string&)        Fx(sys, 1H)
 spin_op Name(const spin_sys&, double)               Rz(sys, 90.0)

	Input	   SOpFct  : A function Returning SOp (spin_op)
                   msys    : A multi_sys spin system
		   double  : Anything (often a pulse angle, e.g. 90.0)
		   string  : Anything (often an isotope name, e.g. "1H")
        	   icomp   : System components affected (-1 = all)
	Output     MOp     : Operator op for multi_sys msys
			     that parallels SOpFct but resides in
        		     the direct product Hilbert space of msys        */

MSVCDLL gen_op multize(spin_op op(const spin_sys&),             const multi_sys& msys);
MSVCDLL gen_op multize(spin_op op(const spin_sys&,  const std::string&),
                                    const std::string&, const multi_sys& msys);
MSVCDLL gen_op multize(spin_op op(const spin_sys&, double),
                             double beta, const multi_sys &msys, int icomp=-1);

// ____________________________________________________________________________
// C    General Operator <=== General Operator Function(spin_sys, gen_op)
// ____________________________________________________________________________

/*                (Often These Are Evolution Or Pulse Functions)

                    Function Form                             Example
 --------------------------------------------------- --------------------------
 Name(const spin_sys&,gen_op&,double)                Iypuls(sys, sigma, 90.0)
 Name(const spin_sys&,gen_op&,int, double)           Iypuls(sys, sigma, 1, 90.)
 Name(const spin_sys&,gen_op&,const string&, double) Iypuls(sys,sigma,"1H",90.)

	Input	   OpFct   : A function Returning Op (gen_op)
                   msys    : A multi_sys spin system
        	   icomp   : System components affected (-1 = all)
		   double  : Anything (often a pulse angle, e.g. 90.0)
                   int     : Anyting (often a spin index)
		   string  : Anything (often an isotope name, e.g. "1H")
	Output     MOp     : Operator op for multi_sys msys
			     that parallels OpFct but resides in
        		     the direct product Hilbert space of msys        */

MSVCDLL gen_op multize(gen_op op(const spin_sys&, const gen_op&, double),
	                      gen_op&, double, const multi_sys&, int icomp=-1);
MSVCDLL gen_op multize(gen_op op(const spin_sys&, const gen_op&, int, double),
                         gen_op&, int, double, const multi_sys&, int icomp=-1);
MSVCDLL gen_op multize(gen_op op(const spin_sys&,const gen_op&,const std::string&, double),
          gen_op&, const std::string&, double, const multi_sys&, int icomp=-1);

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

	Input	   OpFct   : A function Returning Op (gen_op)
                   msys    : A multi_sys spin system
        	   icomp   : System components affected (-1 = all)
		   double  : Anything (often a pulse angle, e.g. 90.0)
                   int     : Anyting (often a spin index)
		   string  : Anything (often an isotope name, e.g. "1H")
	Output     MOp     : Operator op for multi_sys msys
			     that parallels OpFct but resides in
        		     the direct product Hilbert space of msys        */

MSVCDLL gen_op multize(gen_op op(const spin_sys&, int, double),
                  int nspin, double beta, const multi_sys &msys, int icomp=-1);
MSVCDLL gen_op multize(gen_op op(const spin_sys&, const std::string& iso, double),
     const std::string& iso, double beta, const multi_sys &msys, int icomp=-1);
MSVCDLL gen_op multize(gen_op op(const spin_sys&, double),
                             double beta, const multi_sys &msys, int icomp=-1);
MSVCDLL gen_op multize(gen_op op(const spin_sys&, double, double),
                 double phi, double beta, const multi_sys &msys, int icomp=-1);
 


// ----------------------------------------------------------------------------
//                  Generic Functions Which Use Superoperators
//              (Often These Are Relaxation Or Exchange Functions)
// ----------------------------------------------------------------------------

MSVCDLL super_op multize(super_op& SOp, const multi_sys &msys, int icomp);

MSVCDLL super_op multize(super_op SOpFct(const gen_op&), const gen_op& Op, const multi_sys& msys);
 
        // Input        SOpFct  : Function which returns a  superoperator
        //                        given an operator as its only argument
        //              Op      : An operator
        //              msys    : A multi_sys spin system
        // Output       MSOp    : Superoperator SOp for multi_sys msys
        //                        in the direct product Liouville space
        // Note                 : This routine calls an arbitrary function
        //                        of the following form 
        //                              super_op SOpFct(gen_op& Op) 
        //                        where Op is the only function arguments
        //                        and SopFct is the name of the function
        // Note                 : There is no check herein to insure that
        //                        the composite space of Op is the same
        //                        as that which is the default of msys 
 

MSVCDLL super_op multize(super_op sop(const sys_dynamic&,gen_op&,int,int), gen_op&,
                                                    int, int, const multi_sys&);

        // Input        sop     : Function which returns an  
        //                        superoperator given a system, operator,  
        //                        and two integers
        //              H       : An operator (isotropic Hamiltonian)  
        //              type    : Computation type  
        //              level   : Computation level 
        //              msys    : A multi_sys spin system  
        // Output       MSOp    : Superoperator op for multi_sys msys  
        //                        in the direct product Liouville space 
        // Note                 : This routine calls an arbitrary function  
        //                        of the following form  
        //                        super_op Name(sys_dynamic&, gen_op&, int, int) 
        //                        where H, type and level will be the arguments     
 

MSVCDLL super_op multize(super_op sop(const sys_dynamic&,gen_op&,int), gen_op&,
                                                         int, const multi_sys&);
 
        // Input        sop     : Function which returns an 
        //                        superoperator given a system, operator,
        //                        and an integer  
        //              H       : An operator (isotropic Hamiltonian)
        //              level   : Computation level 
        //              msys    : A multi_sys spin system 
        // Output       MSOp    : Superoperator op for multi_sys msys
        //                        in the direct product Liouville space
        // Note                 : This routine calls an arbitrary function
        //                        of the following form 
        //                        super_op Name(sys_dynamic&, gen_op&, int, int)
        //                        where H, type and level will be the arguments

// ____________________________________________________________________________
//                PRODUCTION OF MULTISPIN SYSTEM SUPEROPERATORS
// ____________________________________________________________________________
 
MSVCDLL basis D_basis(const multi_sys& msys);
MSVCDLL std::vector<double> qStateLS(const multi_sys& msys, int I);
MSVCDLL row_vector LS_qState_bra(const multi_sys& msys, int i);
MSVCDLL row_vector LS_qState_ket(const multi_sys& msys, int i);

#endif						// MultiLib.h
