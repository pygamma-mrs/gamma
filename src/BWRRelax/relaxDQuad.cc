/* relaxDQuad.cc *********************************************************
**									**
** 	                          G A M M A				**
**									**
** NMR Dipole-Quadrupole X-Correlation Quad           Implementation  	**
**						 			**
**      Copyright (c) 1997                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**								 	**
**      $Header: $
**								 	**
*************************************************************************/


/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** These are functions dealing with dipole-quadrupole cross correlation	**
** effects in NMR as according to WBR relaxation theory.                **
**                                                                      **
*************************************************************************/

#ifndef _relax_DQuad_cc_			// Is this file already included?
#define _relax_DQuad_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <BWRRelax/relaxDQuad.h>			// Include the header info
#include <BWRRelax/relaxNMR.h>			// Include base WBR functions
#include <BWRRelax/relaxQuad.h>			// Include quadrupolar relaxation
#include <BWRRelax/relaxDip.h>			// Include dipolar relaxation
#include <BWRRelax/relaxRF.h>			// Include rf-field functions
#include <LSLib/sys_dynamic.h>
#include <Level1/nmr_tensor.h>

using std::cout;

// ______________________________________________________________________
// ************* Dipole & Quad Cross Relaxation Superoperators ***********
// ______________________________________________________________________



	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : Dipole & Quad cross relaxation
	//			        superoperator
	// Note			      :	Computed in the eigenbasis of Ho


super_op RDQX(sys_dynamic& sys, gen_op& Ho, int level)
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

   double *w;					// Set up for system energy levels (LAB)
   w = new double[hs];			//  compiler did not like double w[hs], as hs is not known at compile time,

//   double w[hs];				// Set up for system energy levels (LAB)
   
   gen_op Holab; 				// these will be in Hertz, ~10**8 in value
   if(abs(level) > 1)				// Needed for higher level computations
     {						// involving spectral density functions
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);				//	This is H in the lab frame
     Holab.Op_base(Ho, 1.e-7);			//	It must be in Ho eigenbasis
     if(!Holab.test_EBR())			//	and should recognize its EBR
       cout << "\n\tWarning relax_DQuad: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

//	      Get the Dipolar & Quad Cross Relaxation Superoperator

   RDQX(LOp, sys, Ho, w, taus, chi, level);	// Set LOp to Dip & Quad R superop
   return LOp;
   }


void RDQX(super_op& LOp, sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Dipole & Quad cross relaxation superop
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

  matrix xiQs = xiQ(sys);			// Matrix of Xi's for Quad
  matrix xiDs = xiD(sys);			// Matrix of Xi's for Dipoles

//			Prepare the Quad Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *Tquad;				// Spin tensors for each spin
   Tquad = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)				// Set them all to Quad spin tensors
     if(Re(xiQs.get(i,i)))
       Tquad[i] = T_Q(sys,i);

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

//			Prepare the Quad Spatial Tensors

   space_T *Acsa; 				// Spatial tensors for each spin
   Acsa = new space_T[ns];			// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them all to Quad space tensors
     if(Re(xiQs.get(i,i)))
       Acsa[i] = sys.TQ(i);

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

//	    Determine Dipole & Quad Cross Relaxation Superoperator

   Rijk(LOp, sys, Ho, w, xiDs, xiQs, Adip, Acsa,
                  Tdip, Tquad, taus, chi, 0, level);
   if(level == 4)
     Rkij(LOp, sys, Ho, w, xiQs, xiDs, Acsa, Adip,
                  Tquad, Tdip, taus, chi, 0, level);
   else
     LOp *= 2.0;
   return;
   }


// ______________________________________________________________________
// ************** Dipole-Quad Relaxation Superoperators ******************
// ______________________________________________________________________


super_op RDQ(sys_dynamic& sys, gen_op& Ho, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : Dipole-Quad relaxation superoperator
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

   double *w;					// Set up for system energy levels (LAB)
   w = new double[hs];			//  compiler did not like double w[hs];

// double w[hs];				// Set up for system energy levels (LAB)
   
   gen_op Holab; 				// these will be in Hertz, ~10**8 in value
   if(abs(level) > 1)				// Needed for higher level computations
     {						// involving spectral density functions
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);				//	This is H in the lab frame
     Holab.Op_base(Ho, 1.e-7);			//	It must be in Ho eigenbasis
     if(!Holab.test_EBR())			//	and should recognize its EBR
       cout << "\n\tWarning relax_DQuad: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

//	         Get the Dipolar-Quad Relaxation Superoperator

   RDQ(LOp, sys, Ho, w, taus, chi, level);	// Set LOp to Dip-Quad R superop
   return LOp;
   }


void RDQ(super_op& LOp, sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Dip-Quad relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

  matrix xiQs = xiQ(sys);			// Matrix of Xi's for Quad
  matrix xiDs = xiD(sys);			// Matrix of Xi's for Dipoles

//			Prepare the Quad Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *Tquad;				// Spin tensors for each spin
   Tquad = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)			// Set them all to Quad spin tensors
     if(Re(xiQs.get(i,i)))
       Tquad[i] = T_Q(sys,i);

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

//			Prepare the Quad Spatial Tensors

   space_T *Acsa; 				// Spatial tensors for each spin
   Acsa = new space_T[ns];			// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them all to Quad space tensors
     if(Re(xiQs.get(i,i)))
       Acsa[i] = sys.TQ(i);

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

//	    Determine Dipole-Quad Relaxation Superoperator

    Rijk(LOp, sys, Ho, w, xiDs, xiQs, Adip, Acsa,
                  Tdip, Tquad, taus, chi, 0, level);
   return;
   }


// -------------- Dipole-Quad, With An Applied RF-Field ------------------


// ______________________________________________________________________
// ************** Quad-Dipole Relaxation Superoperators ******************
// ______________________________________________________________________


super_op RQD(sys_dynamic& sys, gen_op& Ho, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			level : Relaxation treatment level
	// Output		LOp   : Quad-Dipole relaxation superoperator
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

   double *w;					// Set up for system energy levels (LAB)
   w = new double[hs];			//  compiler did not like double w[hs];

//   double w[hs];				// Set up for system energy levels (LAB)
  
   gen_op Holab; 				// these will be in Hertz, ~10**8 in value
   if(abs(level) > 1)				// Needed for higher level computations
     {						// involving spectral density functions
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);				//	This is H in the lab frame
     Holab.Op_base(Ho, 1.e-7);			//	It must be in Ho eigenbasis
     if(!Holab.test_EBR())			//	and should recognize its EBR
       cout << "\n\tWarning relax_DQuad: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

//	         Get the Quad-Dipolar Relaxation Superoperator

   RQD(LOp, sys, Ho, w, taus, chi, level);	// Set LOp to Quad-Dip R superop
   return LOp;
   }


void RQD(super_op& LOp, sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Quad-Dipole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

  matrix xiQs = xiQ(sys);			// Matrix of Xi's for Quad
  matrix xiDs = xiD(sys);			// Matrix of Xi's for Dipoles

//			Prepare the Quad Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *Tquad;				// Spin tensors for each spin
   Tquad = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)				// Set them all to Quad spin tensors
     if(Re(xiQs.get(i,i)))
       Tquad[i] = T_Q(sys,i);

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

//			Prepare the Quad Spatial Tensors

   space_T *Acsa; 				// Spatial tensors for each spin
   Acsa = new space_T[ns];			// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them all to Quad space tensors
     if(Re(xiQs.get(i,i)))
       Acsa[i] = sys.TQ(i);

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

//	    Determine Quad-Dipole Relaxation Superoperator

    Rkij(LOp, sys, Ho, w, xiQs, xiDs, Acsa, Adip,
                  Tquad, Tdip, taus, chi, 0, level);
   return;
   }

// -------------- Quad-Dipole, With An Applied RF-Field ------------------

super_op RQDrf(sys_dynamic& sys, gen_op& Heff, double Wrflab, int level)

	// Input		sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian (Hertz)
	//			Wrflab: RF-Field frequency, lab frame (Hertz)
	//			level : Relaxation treatment level
	// Output		LOp   : Quad-Dipole relaxation superoperator
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

   double *w;					// Set up for system energy levels (LAB)
   w = new double[hs];			//  compiler did not like double w[hs];

//   double w[hs];				// Set up for system energy levels (LAB)
   
   if(abs(level) > 1)				// Needed for higher level computations
     Heff.eigvals(w);

//	         Get the Quad-Dipolar Relaxation Superoperator

   RQDrf(LOp,sys,Heff,w,Wrflab,taus,chi,level);	// Set LOp to Quad-Dip R superop
   return LOp;
   }


void RQDrf(super_op& LOp, sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			Wrflab: Field frequency lab frame (Hz)
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	//			level : Relaxation treatment level
	// Output		LOp   : Quad-Dipole relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff

   {
//		    Prepare the Interaction Constants

  matrix xiQs = xiQ(sys);			// Matrix of Xi's for Quad
  matrix xiDs = xiD(sys);			// Matrix of Xi's for Dipoles

//			Prepare the Quad Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *Tquad;				// Spin tensors for each spin
   Tquad = new spin_T[ns];			// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)				// Set them all to Quad spin tensors
     if(Re(xiQs.get(i,i)))
       Tquad[i] = T_Q(sys,i);

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

//			Prepare the Quad Spatial Tensors

   space_T *Acsa; 				// Spatial tensors for each spin
   Acsa = new space_T[ns];			// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them all to Quad space tensors
     if(Re(xiQs.get(i,i)))
       Acsa[i] = sys.TQ(i);

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

//	    Determine Quad-Dipole Relaxation Superoperator

    Rrfkij(LOp, sys, Heff, w, Wrflab, xiQs, xiDs, Acsa, Adip,
                  Tquad, Tdip, taus, chi, 0, level);
   return;
   }


#endif 						// relaxDQuad.cc
