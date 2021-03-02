/* relaxDCSA.cc **********************************************************
**                                                                      **
**                                G A M M A                             **
**                                                                      **
**      NMR Dipole-CSA X-Correlation                 Implementation     **
**                                                                      **
**	Scott A. Smith						 	**
**      Copyright (c) 2002                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fuer physikalische Chemie				 	**
**	8092 Zurich / Switzerland				 	**
**									**
**	University of California, Santa Barbara				**
**	Department of Chemistry						**
**	Santa Barbara CA. 93106 USA					**
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** These are functions dealing with dipole-CSA cross correlation 	**
** effects in NMR as according to WBR relaxation theory.                **
**                                                                      **
*************************************************************************/

#ifndef _relax_DCSA_cc_		// Is this file already included?
#define _relax_DCSA_cc_ 1	// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation          // This is the implementation
#endif


#include <BWRRelax/relaxDCSA.h>		// Include header
#include <Level1/nmr_tensor.h>		// Inlcude common MR spin tensors
#include <stdlib.h>


// ______________________________________________________________________
// ************* Dipole & CSA Cross Relaxation Superoperators ***********
// ______________________________________________________________________


super_op RDCX(const sys_dynamic& sys, gen_op& Ho, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : Dipole & CSA cross relaxation
	//			        superoperator
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
   double *w;		   		// Set up for system energy levels (LAB)
   w = new double[hs];
   gen_op Holab; 				// these will be in Hertz, ~10**8 in value
   if(abs(level) > 1)				// Needed for higher level computations
     {						// involving spectral density functions
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);				//	This is H in the lab frame
     Holab.Op_base(Ho, 1.e-7);			//	It must be in Ho eigenbasis
     if(!Holab.test_EBR())			//	and should recognize its EBR
       std::cout << "\n\tWarning relax_DCSA: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

//	      Get the Dipolar & CSA Cross Relaxation Superoperator

   RDCX(LOp, sys, Ho, w, taus, chi, level);	// Set LOp to Dip & CSA R superop
   delete [] w;
   return LOp;
   }


void RDCX(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Dipole & CSA cross relaxation superop
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

  matrix xiCs = xiCSA(sys);			// Matrix of Xi's for CSA
  matrix xiDs = xiD(sys);			// Matrix of Xi's for Dipoles

//			Prepare the CSA Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *Tcsa;				// Spin tensors for each spin
   Tcsa = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)				// Set them all to CSA spin tensors
     if(Re(xiCs.get(i,i)))
       Tcsa[i] = T_CS2(sys,i);

//		      Prepare the Dipolar Spin Tensors

   int ndip = sys.dipoles();			// Total number of dipoles
   spin_T *Tdip;				// Spin tensors for each dipole
   Tdip = new spin_T[ndip];			// Compiler didn't like spin_T T[ndip]
   int ij=0;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xiDs.get(i,j)))
         {
         Tdip[ij] = T_D(sys,i,j);
         ij++;
	 }

//			Prepare the CSA Spatial Tensors

   space_T *Acsa; 				// Spatial tensors for each spin
   Acsa = new space_T[ns];			// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them all to CSA space tensors
     if(Re(xiCs.get(i,i)))
       Acsa[i] = sys.TC(i);

//	               Prepare the Dipolar Space Tensors

   space_T *Adip;
   Adip = new space_T[ndip];			// Compiler didn't like space_T A[ndip]
   ij=0;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xiDs.get(i,j)))
         {
//         Adip[ij] = sys.TD(i,j);
         ij++;
	 }

//	    Determine Dipole & CSA Cross Relaxation Superoperator

   Rijk(LOp, sys, Ho, w, xiDs, xiCs, Adip, Acsa,
                  Tdip, Tcsa, taus, chi, 0, level);
   if(level == 4)
     Rkij(LOp, sys, Ho, w, xiCs, xiDs, Acsa, Adip,
                  Tcsa, Tdip, taus, chi, 0, level);
   else
     LOp *= 2.0;
   return;
   }


// ______________________________________________________________________
// ************** Dipole-CSA Relaxation Superoperators ******************
// ______________________________________________________________________


super_op RDC(const sys_dynamic& sys, gen_op& Ho, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : Dipole-CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


   {
//   int ns = sys.spins();			// Total number of spins
//   int ndip = sys.dipoles();			// Total number of dipoles
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
       std::cout << "\n\tWarning relax_DCSA: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

//	         Get the Dipolar-CSA Relaxation Superoperator

   RDC(LOp, sys, Ho, w, taus, chi, level);	// Set LOp to Dip-CSA R superop
   delete [] w;
   return LOp;
   }


void RDC(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Dip-CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

  matrix xiCs = xiCSA(sys);			// Matrix of Xi's for CSA
  matrix xiDs = xiD(sys);			// Matrix of Xi's for Dipoles

//			Prepare the CSA Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *Tcsa;				// Spin tensors for each spin
   Tcsa = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)			// Set them all to CSA spin tensors
     if(Re(xiCs.get(i,i)))
       Tcsa[i] = T_CS2(sys,i);

//		      Prepare the Dipolar Spin Tensors

   int ndip = sys.dipoles();			// Total number of dipoles
   spin_T *Tdip;				// Spin tensors for each dipole
   Tdip = new spin_T[ndip];			// Compiler didn't like spin_T T[ndip]
   int ij=0;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xiDs.get(i,j)))
         {
         Tdip[ij] = T_D(sys,i,j);
         ij++;
	 }

//			Prepare the CSA Spatial Tensors

   space_T *Acsa; 				// Spatial tensors for each spin
   Acsa = new space_T[ns];			// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them all to CSA space tensors
     if(Re(xiCs.get(i,i)))
       Acsa[i] = sys.TC(i);

//	               Prepare the Dipolar Space Tensors

   space_T *Adip;
   Adip = new space_T[ndip];			// Compiler didn't like space_T A[ndip]
   ij=0;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xiDs.get(i,j)))
         {
//         Adip[ij] = sys.TD(i,j);
         ij++;
	 }

//	    Determine Dipole-CSA Relaxation Superoperator

    Rijk(LOp, sys, Ho, w, xiDs, xiCs, Adip, Acsa,
                  Tdip, Tcsa, taus, chi, 0, level);
   return;
   }


// -------------- Dipole-CSA, With An Applied RF-Field ------------------


// ______________________________________________________________________
// ************** CSA-Dipole Relaxation Superoperators ******************
// ______________________________________________________________________


super_op RCD(const sys_dynamic& sys, gen_op& Ho, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Dipole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho


   {
//   int ns = sys.spins();			// Total number of spins
//   int ndip = sys.dipoles();			// Total number of dipoles
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
       std::cout << "\n\tWarning relax_DCSA: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

//	         Get the CSA-Dipolar Relaxation Superoperator

   RCD(LOp, sys, Ho, w, taus, chi, level);	// Set LOp to CSA-Dip R superop
   delete [] w;
   return LOp;
   }


void RCD(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Dipole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

  matrix xiCs = xiCSA(sys);			// Matrix of Xi's for CSA
  matrix xiDs = xiD(sys);			// Matrix of Xi's for Dipoles

//			Prepare the CSA Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *Tcsa;				// Spin tensors for each spin
   Tcsa = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)				// Set them all to CSA spin tensors
     if(Re(xiCs.get(i,i)))
       Tcsa[i] = T_CS2(sys,i);

//		      Prepare the Dipolar Spin Tensors

   int ndip = sys.dipoles();			// Total number of dipoles
   spin_T *Tdip;				// Spin tensors for each dipole
   Tdip = new spin_T[ndip];			// Compiler didn't like spin_T T[ndip]
   int ij=0;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xiDs.get(i,j)))
         {
         Tdip[ij] = T_D(sys,i,j);
         ij++;
	 }

//			Prepare the CSA Spatial Tensors

   space_T *Acsa; 				// Spatial tensors for each spin
   Acsa = new space_T[ns];			// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them all to CSA space tensors
     if(Re(xiCs.get(i,i)))
       Acsa[i] = sys.TC(i);

//	               Prepare the Dipolar Space Tensors

   space_T *Adip;
   Adip = new space_T[ndip];			// Compiler didn't like space_T A[ndip]
   ij=0;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xiDs.get(i,j)))
         {
//         Adip[ij] = sys.TD(i,j);
         ij++;
	 }

//	    Determine CSA-Dipole Relaxation Superoperator

    Rkij(LOp, sys, Ho, w, xiCs, xiDs, Acsa, Adip,
                  Tcsa, Tdip, taus, chi, 0, level);
   return;
   }

// -------------- CSA-Dipole, With An Applied RF-Field ------------------

super_op RCDrf(const sys_dynamic& sys, gen_op& Heff, double Wrflab, int level)

	// Input		sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian (Hertz)
	//			Wrflab: RF-Field frequency, lab frame (Hertz)
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Dipole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff


   {
//   int ns = sys.spins();			// Total number of spins
//   int ndip = sys.dipoles();			// Total number of dipoles
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

//	         Get the CSA-Dipolar Relaxation Superoperator

   RCDrf(LOp,sys,Heff,w,Wrflab,taus,chi,level);	// Set LOp to CSA-Dip R superop
   delete [] w;
   return LOp;
   }


void RCDrf(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			Wrflab: Field frequency lab frame (Hz)
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : CSA-Dipole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff

   {
//		    Prepare the Interaction Constants

  matrix xiCs = xiCSA(sys);			// Matrix of Xi's for CSA
  matrix xiDs = xiD(sys);			// Matrix of Xi's for Dipoles

//			Prepare the CSA Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *Tcsa;				// Spin tensors for each spin
   Tcsa = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)				// Set them all to CSA spin tensors
     if(Re(xiCs.get(i,i)))
       Tcsa[i] = T_CS2(sys,i);

//		      Prepare the Dipolar Spin Tensors

   int ndip = sys.dipoles();			// Total number of dipoles
   spin_T *Tdip;				// Spin tensors for each dipole
   Tdip = new spin_T[ndip];			// Compiler didn't like spin_T T[ndip]
   int ij=0;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xiDs.get(i,j)))
         {
         Tdip[ij] = T_D(sys,i,j);
         ij++;
	 }

//			Prepare the CSA Spatial Tensors

   space_T *Acsa; 				// Spatial tensors for each spin
   Acsa = new space_T[ns];			// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them all to CSA space tensors
     if(Re(xiCs.get(i,i)))
       Acsa[i] = sys.TC(i);

//	               Prepare the Dipolar Space Tensors

   space_T *Adip;
   Adip = new space_T[ndip];			// Compiler didn't like space_T A[ndip]
   ij=0;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xiDs.get(i,j)))
         {
//         Adip[ij] = sys.TD(i,j);
         ij++;
	 }

//	    Determine CSA-Dipole Relaxation Superoperator

    Rrfkij(LOp, sys, Heff, w, Wrflab, xiCs, xiDs, Acsa, Adip,
                  Tcsa, Tdip, taus, chi, 0, level);
   return;
   }


#endif /* __RELAX_DCSA_CC__ */

