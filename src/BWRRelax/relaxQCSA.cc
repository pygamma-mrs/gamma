/* relaxQCSA.cc **********************************************************
**									**
** 	                           G A M M A				**
**									**
**	NMR Quadrupole-CSA X-Correlation 	Implementation 		**
**                                                                      **
**      Copyright (c) 1997                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/


/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** Herein are functions dealing with quadrupole-CSA cross correlation   **
** effects in NMR as according to WBR relaxation theory.                **
**                                                                      **
*************************************************************************/

#ifndef _relax_QCSA_cc_			// Is this file already included?
#define _relax_QCSA_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation         		// This is the implementation
#endif


#include <BWRRelax/relaxQCSA.h>		// Include header description
#include <BWRRelax/relaxNMR.h>		// Include base WBR functions
#include <BWRRelax/relaxQuad.h>		// Include quadrupolar relaxation
#include <BWRRelax/relaxCSA.h>		// Include CSA relaxation
#include <BWRRelax/relaxRF.h>		// Include rf-field functions
#include <Level1/nmr_tensor.h>		// Include commont spin tensors
#include <LSLib/sys_dynamic.h>		// Include aniostropic systems
#include <stdlib.h>

// ______________________________________________________________________
// *********** Quadrupole & CSA Cross Relaxation Superoperators *********
// ______________________________________________________________________


super_op RQCX(const sys_dynamic& sys, gen_op& Ho, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : Quadrupole & CSA cross relaxation
	//			        superoperator
	// Note			      :	Computed in the eigenbasis of Ho


   {
   int hs = sys.HS();				// Total system Hilbert space
   int ls = hs*hs;				// Total system Liouville space
   Ho.set_EBR();				// Insure Ho is in its eigenbasis (Hz)
   matrix mx(ls, ls, 0.0);			// Construct zero superoperator
   super_op LOp(mx, Ho.get_basis());		// Put into Ho eigenbasis

//		    Prepare the Spectral Density Components

   double taus[5];				// Get 5 taus for spectral densities
   taust(taus, sys.taus());			// Set 5 taus
   double chi = chit(sys.taus());		// Get the system chi value
   double *w;				// Set up for system energy levels (LAB)
   w = new double[hs];
   gen_op Holab; 				// these will be in Hertz, ~10**8 in value
   if(abs(level) > 1)				// Needed for higher level computations
     {						// involving spectral density functions
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);				//	This is H in the lab frame
     Holab.Op_base(Ho, 1.e-7);			//	It must be in Ho eigenbasis
     if(!Holab.test_EBR())			//	and should recognize its EBR
       std::cout << "\n\tWarning relax_QCSA: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

//	      Get the Quadrupolar & CSA Cross Relaxation Superoperator

   RQCX(LOp, sys, Ho, w, taus, chi, level);	// Set LOp to Quad & CSA R superop
   delete [] w;
   return LOp;
   }


void RQCX(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                                   double* taus, double chi, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Quadrupole & CSA cross relaxation superop
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

  matrix xiCs = xiCSA(sys);			// Matrix of Xi's for CSA
  matrix xiQs = xiQ(sys);			// Matrix of Xi's for Quadrupoles

//		   Prepare the CSA & Quad Spin Tensors

  int ns = sys.spins();				// Total number of spins
  spin_T *Tcsa;					// CSA spin tensors, each spin
  spin_T *Tquad;				// Quad spin tensors, each spin
  Tcsa = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
  Tquad = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
  int i;
  for(i=0; i<ns; i++)				// Set CSA & Quad spin tensors
    {
    if(xiCs.getRe(i,i)) Tcsa[i] = T_CS2(sys,i);	
    if(xiQs.getRe(i,i)) Tquad[i] = T_Q(sys,i);
    }

//		Prepare the CSA & Quad Spatial Tensors

   space_T *Acsa; 				// CSA spatial tensors for each spin
   space_T *Aquad; 				// Quad spatial tensors for each spin
   Acsa = new space_T[ns];			// Compiler didn't like space_T A[ns]
   Aquad = new space_T[ns];			// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set CSA & Quad space tensors
     {
     if(xiCs.getRe(i,i)) Acsa[i] = sys.TC(i);
     if(xiQs.getRe(i,i)) Aquad[i] = sys.TQ(i);
     }

//	    Determine Quadrupole & CSA Cross Relaxation Superoperator

   Rij(LOp, sys, Ho, w, xiQs, xiCs, Aquad, Acsa,
                  Tquad, Tcsa, taus, chi, 0, level);
   if(level == 4)
     Rkij(LOp, sys, Ho, w, xiCs, xiQs, Acsa, Aquad,
                  Tcsa, Tquad, taus, chi, 0, level);
   else
     LOp *= 2.0;
   return;
   }


// ______________________________________________________________________
// ************ Quadrupole-CSA Relaxation Superoperators ****************
// ______________________________________________________________________


super_op RQC(const sys_dynamic& sys, gen_op& Ho, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : Quadrupole-CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


   {
   int hs = sys.HS();				// Total system Hilbert space
   int ls = hs*hs;				// Total system Liouville space
   Ho.set_EBR();				// Insure Ho is in its eigenbasis (Hz)
   matrix mx(ls, ls, 0.0);			// Construct zero superoperator
   super_op LOp(mx, Ho.get_basis());

//		    Prepare the Spectral Density Components

   double taus[5];				// Get 5 taus for spectral densities
   taust(taus, sys.taus());
   double chi = chit(sys.taus());		// Get the system chi value
   double *w;				// Set up for system energy levels (LAB)
   w = new double[hs];
   gen_op Holab; 				// these will be in Hertz, ~10**8 in value
   if(abs(level) > 1)				// Needed for higher level computations
     {						// involving spectral density functions
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);				//	This is H in the lab frame
     Holab.Op_base(Ho, 1.e-7);			//	It must be in Ho eigenbasis
     if(!Holab.test_EBR())			//	and should recognize its EBR
       std::cout << "\n\tWarning relax_QCSA: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

//	         Get the Quadrupolar-CSA Relaxation Superoperator

   RQC(LOp, sys, Ho, w, taus, chi, level);	// Set LOp to Quad-CSA R superop
   delete [] w;
   return LOp;
   }


void RQC(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Quad-CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

  matrix xiCs = xiCSA(sys);			// Matrix of Xi's for CSA
  matrix xiQs = xiQ(sys);			// Matrix of Xi's for Quadrupoles

//		   Prepare the CSA & Quad Spin Tensors

  int ns = sys.spins();				// Total number of spins
  spin_T *Tcsa;					// CSA spin tensors, each spin
  spin_T *Tquad;				// Quad spin tensors, each spin
  Tcsa = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
  Tquad = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
  int i;
  for(i=0; i<ns; i++)				// Set CSA & Quad spin tensors
    {
    if(xiCs.getRe(i,i)) Tcsa[i] = T_CS2(sys,i);	
    if(xiQs.getRe(i,i)) Tquad[i] = T_Q(sys,i);
    }

//		Prepare the CSA & Quad Spatial Tensors

   space_T *Acsa; 				// CSA spatial tensors for each spin
   space_T *Aquad; 				// Quad spatial tensors for each spin
   Acsa = new space_T[ns];			// Compiler didn't like space_T A[ns]
   Aquad = new space_T[ns];			// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set CSA & Quad space tensors
     {
     if(xiCs.getRe(i,i)) Acsa[i] = sys.TC(i);
     if(xiQs.getRe(i,i)) Aquad[i] = sys.TQ(i);
     }

//	    Determine Quadrupole-CSA Relaxation Superoperator

   Rij(LOp, sys, Ho, w, xiQs, xiCs, Aquad, Acsa,
                            Tquad, Tcsa, taus, chi, 0, level);
   return;
   }


// ------------ Quadrupole-CSA, With An Applied RF-Field ----------------


// ______________________________________________________________________
// ************ CSA-Quadrupole Relaxation Superoperators ****************
// ______________________________________________________________________


super_op RCQ(const sys_dynamic& sys, gen_op& Ho, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Quadrupole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


   {
   int hs = sys.HS();				// Total system Hilbert space
   int ls = hs*hs;				// Total system Liouville space
   Ho.set_EBR();				// Insure Ho is in its eigenbasis (Hz)
   matrix mx(ls, ls, 0.0);			// Construct zero superoperator
   super_op LOp(mx, Ho.get_basis());

//		    Prepare the Spectral Density Components

   double taus[5];				// Get 5 taus for spectral densities
   taust(taus, sys.taus());
   double chi = chit(sys.taus());		// Get the system chi value
   double *w;				// Set up for system energy levels (LAB)
   w = new double[hs];
   gen_op Holab; 				// these will be in Hertz, ~10**8 in value
   if(abs(level) > 1)				// Needed for higher level computations
     {						// involving spectral density functions
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);				//	This is H in the lab frame
     Holab.Op_base(Ho, 1.e-7);			//	It must be in Ho eigenbasis
     if(!Holab.test_EBR())			//	and should recognize its EBR
       std::cout << "\n\tWarning relax_QCSA: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

//	         Get the CSA-Quadrupolar Relaxation Superoperator

   RCQ(LOp, sys, Ho, w, taus, chi, level);	// Set LOp to CSA-Quad R superop
   delete [] w;
   return LOp;
   }


void RCQ(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Quadrupole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

  matrix xiCs = xiCSA(sys);			// Matrix of Xi's for CSA
  matrix xiQs = xiQ(sys);			// Matrix of Xi's for Quadrupoles

//		   Prepare the CSA & Quad Spin Tensors

  int ns = sys.spins();				// Total number of spins
  spin_T *Tcsa;					// CSA spin tensors, each spin
  spin_T *Tquad;				// Quad spin tensors, each spin
  Tcsa = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
  Tquad = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
  int i;
  for(i=0; i<ns; i++)				// Set CSA & Quad spin tensors
    {
    if(xiCs.getRe(i,i)) Tcsa[i] = T_CS2(sys,i);	
    if(xiQs.getRe(i,i)) Tquad[i] = T_Q(sys,i);
    }

//		Prepare the CSA & Quad Spatial Tensors

   space_T *Acsa; 				// CSA spatial tensors for each spin
   space_T *Aquad; 				// Quad spatial tensors for each spin
   Acsa = new space_T[ns];			// Compiler didn't like space_T A[ns]
   Aquad = new space_T[ns];			// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set CSA & Quad space tensors
     {
     if(xiCs.getRe(i,i)) Acsa[i] = sys.TC(i);
     if(xiQs.getRe(i,i)) Aquad[i] = sys.TQ(i);
     }

//	    Determine CSA-Quadrupole Relaxation Superoperator

    Rij(LOp, sys, Ho, w, xiCs, xiQs, Acsa, Aquad,
                  Tcsa, Tquad, taus, chi, 0, level);
   return;
   }

// ------------ CSA-Quadrupole, With An Applied RF-Field ----------------

super_op RCQrf(const sys_dynamic& sys, gen_op& Heff, double Wrflab, int level)

	// Input		sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian (Hertz)
	//			Wrflab: RF-Field frequency, lab frame (Hertz)
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Quadrupole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff


   {
   int hs = sys.HS();				// Total system Hilbert space
   int ls = hs*hs;				// Total system Liouville space
   Heff.set_EBR();				// Insure Heff is in its eigenbasis (Hz)
   matrix mx(ls, ls, 0.0);			// Construct zero superoperator
   super_op LOp(mx, Heff.get_basis());

//		    Prepare the Spectral Density Components

   double taus[5];				// Get 5 taus for spectral densities
   taust(taus, sys.taus());
   double chi = chit(sys.taus());		// Get the system chi value
   double *w;				// Set up for system energy levels (LAB)
   w = new double[hs];
   if(abs(level) > 1)				// Needed for higher level computations
     Heff.eigvals(w);

//	         Get the CSA-Quadrupolar Relaxation Superoperator

   RCQrf(LOp,sys,Heff,w,Wrflab,taus,chi,level);	// Set LOp to CSA-Quad R superop
   delete [] w;
   return LOp;
   }


void RCQrf(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			Wrflab: Field frequency lab frame (Hz)
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Quadrupole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff

   {
//		    Prepare the Interaction Constants

  matrix xiCs = xiCSA(sys);			// Matrix of Xi's for CSA
  matrix xiQs = xiQ(sys);			// Matrix of Xi's for Quadrupoles

//		   Prepare the CSA & Quad Spin Tensors

  int ns = sys.spins();				// Total number of spins
  spin_T *Tcsa;					// CSA spin tensors, each spin
  spin_T *Tquad;				// Quad spin tensors, each spin
  Tcsa = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
  Tquad = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
  int i;
  for(i=0; i<ns; i++)				// Set CSA & Quad spin tensors
    {
    if(xiCs.getRe(i,i)) Tcsa[i] = T_CS2(sys,i);	
    if(xiQs.getRe(i,i)) Tquad[i] = T_Q(sys,i);
    }

//		Prepare the CSA & Quad Spatial Tensors

   space_T *Acsa; 				// CSA spatial tensors for each spin
   space_T *Aquad; 				// Quad spatial tensors for each spin
   Acsa = new space_T[ns];			// Compiler didn't like space_T A[ns]
   Aquad = new space_T[ns];			// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set CSA & Quad space tensors
     {
     if(xiCs.getRe(i,i)) Acsa[i] = sys.TC(i);
     if(xiQs.getRe(i,i)) Aquad[i] = sys.TQ(i);
     }

//	    Determine CSA-Quadrupole Relaxation Superoperator

    Rrfkij(LOp, sys, Heff, w, Wrflab, xiCs, xiQs, Acsa, Aquad,
                  Tcsa, Tquad, taus, chi, 0, level);
   return;
   }


#endif 						// __RELAX_QCSA_CC_

