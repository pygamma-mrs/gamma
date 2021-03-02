/* relax_Dip.cc **************************************************
**                                                              **
**                        G A M M A                             **
**                                                              **
**      NMR Dipolar Relaxation Functions    Implementation      **
**                                                              **
**      Copyright (c) 1991, 1992, 1993                          **
**      Scott A. Smith                                          **
**                                                              **
**      Eidgenoessische Technische Hochschule                   **
**      Labor fuer physikalische Chemie                         **
**      8092 Zurich / Switzerland                               **
**                                                              **
**	University of California, Santa Barbara			**
**	Department of Chemistry					**
**	Santa Barbara CA. 93106 USA				**
**							 	**
**      $Header: $
**						 		**
*****************************************************************/

/*****************************************************************
**							 	**
**  Description						 	**
**							 	**
**  The functions herein provide easy access to many of the 	**
**  quantities associated with dipolar relaxation in isotropic	**
**  NMR systems.						**
**						 		**
*****************************************************************/

#ifndef _relax_Dip_cc_		// Is this file already included?
#define _relax_Dip_cc_ 1	// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#endif

#include <BWRRelax/relaxDip.h>		// Include the header
#include <Level1/SpinT.h>		 // Include spin tensors
#include <BWRRelax/relaxNMR.h>		// Include relaxation engine
#include <BWRRelax/relaxRF.h>		 // Include relaxation with rf
#include <HSLib/HSham.h>		// Include NMR Hamiltonians
#include <Level1/nmr_tensor.h>		// Include NMR tensors
#include <stdlib.h>

const double mu0d4pi = 1.0e-7;		// mu0/4p (J-sec C  m  )
const double pi2 = 6.283185307;		// 2*pi
const double hbar = 1.05459e-34;	// hbar (J-sec)


// ____________________________________________________________________________
// ***************** Dipole-Dipole Relaxation Superoperators ******************
// ____________________________________________________________________________


// ---------------------- Dipole-Dipole, No RF-Field --------------------------


void RDD(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
                          double* taus, double chi, int type, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note				Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

   matrix xis = xiD(sys);			// Dipolar interaction constants

//	                Prepare the Spin Tensors

   int ns = sys.spins();			// Total number of spins
   int ndip = sys.dipoles();			// Total number of dipoles
   spin_T *T;					// Spin tensors for each dipole
   T = new spin_T[ndip];			// Compiler didn't like spin_T T[ndip]
   int ij=0;
   int i=0;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xis.get(i,j)))
         {
         T[ij] = T_D(sys,i,j);
         ij++;
	 }

//	               Prepare the Space Tensors

   space_T *A;					// These 2 steps as compiler
   A = new space_T[ndip];			// didn't like space_T A[ndip]
   ij=0;					// Begin with first dipole
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute dipolar spin tensors
       if(Re(xis.get(i,j)))
         {
//         A[ij] = sys.TC(i,j);
         ij++;
	 }

//	           Dipolar Relaxation Superoperator 

   Rijkl(LOp,sys,Ho,w,xis,xis,A,A,T,T,taus,chi,type,level);
   }


super_op RDD(const sys_dynamic& sys, gen_op& Ho, int type, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note				Computed in the eigenbasis of Ho

   {
//   int ns = sys.spins();			// Total number of spins
//   int ndip = sys.dipoles();			// Total number of dipoles
   int hs = sys.HS();				// Total system Hilbert space
   int ls = hs*hs;				// Total system Liouville space
   Ho.set_EBR();				// Insure Ho is in its eigenbasis (Hz)
   matrix mx(ls, ls, complex0);			// Construct zero superoperator
   super_op LOp(mx, Ho.get_basis());

//		    Prepare the Spectral Density Components

   double taus[5];				// Get 5 taus for spectral densities
   taust(taus, sys.taus());
   double chi = chit(sys.taus());		// Get the system chi value
   double *w;				// Set up for system energy levels (LAB)
   w = new double[hs];
   if(abs(level) > 1)				// Needed for higher level computations
     {						// where spectral densities are evaluated
     gen_op HZ = Hz(sys); 
     HZ.Op_base(Ho, 1.e-7);
     if(!HZ.test_EBR())
       std::cout << "\n\tWarning relax_Dip: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     for(int k=0; k<hs; k++)
       w[k] = Re(HZ.get(k,k)) + Re(Ho.get(k,k));
     }

//	         Get the Pure Dipolar Relaxation Superoperator

   RDD(LOp,sys,Ho,w,taus,chi,type,level);	// Set LOp to Dip-Dip R superop
   delete [] w;
   return LOp;
   }


super_op RDD(const spin_system& sys, gen_op& Ho, double tau,
				 matrix& dist, int type, int level)

	// Input		sys   : Spin system
	//			Ho    : General operator
	// 			tau   : Correlation time (in seconds)
	// 			dist  : Distance matrix (in meters)
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	//				computed in the basis of Ho

   {
   int ns = sys.spins();			// Total number of spins
   int ndip = 0;				// Total number of dipoles
   int i;
   for(i=1; i<= ns; i++)
     ndip += i;
   int hs = sys.HS();				// Total system Hilbert space
   int ls = hs*hs;				// Total system Liouville space
   Ho.set_EBR();				// Insure Ho is in its eigenbasis (Hz)
   matrix mx(ls, ls, 0.0);			// Construct zero superoperator
   super_op LOp(mx, Ho.get_basis());

//		    Prepare the Spectral Density Components

   double *w;				// Set up for system energy levels (LAB)
   w = new double[hs];
   gen_op Holab; 				// these will be in Hertz, ~10**8 in value
   if(abs(level) > 1)				// Needed for higher level computations
     {						// Where spectral densities are evaluated
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);				//	This is H is the lab frame
     Holab.Op_base(Ho,1.e-7);			//	It must in the Ho eigenbasis
     if(!Holab.test_EBR())	 		//	and should know its an eigenbasis
       std::cout << "\n\tWarning relax_Dip: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

//		    Prepare the Interaction Constants

   matrix xis = xiD(sys, dist);

//	                Prepare the Spin Tensors

   spin_T *T;					// Spin tensors for each dipole
   T = new spin_T[ndip];			// Compiler didn't like spin_T T[ndip]
   int ij=0;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xis.get(i,j)))
         {
         T[ij] = T_D(sys,i,j);
         ij++;
	 }

//	           Dipolar Relaxation Superoperator 

   Rijkl(LOp,sys,Ho,w,xis,xis,T,T,tau,type,level);
   delete [] w;

   return LOp;
   }


super_op RDD_Jgen(const sys_dynamic& sys, gen_op& Ho, int type, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note				Computed in the eigenbasis of Ho

   {
   int ns = sys.spins();			// Total number of spins
   int ndip = sys.dipoles();			// Total number of dipoles
   int hs = sys.HS();				// Total system Hilbert space
   int ls = hs*hs;				// Total system Liouville space
   Ho.set_EBR();				// Insure Ho is in its eigenbasis (Hz)
   matrix mx(ls, ls, 0.0);			// Construct zero superoperator
   super_op LOp(mx, Ho.get_basis());

//			Prepare the Spin Tensors

   spin_T *T;					// Spin tensors for each dipole
   T = new spin_T[ndip];			// Compiler didn't like spin_T T[ndip]
   int ij=0;
   for (int i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for (int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       {
       T[ij] = T_D(sys,i,j);
       ij++;
       }
//
//		    Prepare the Interaction Constants
//		     and Spectral Density Components

   matrix xis = xiD(sys);			// Dipolar interaction constants
//   matrix xis = sys.xiD_matrix();		// Dipolar interaction constants
   double tau = sys.taux();			// Use only 1 tau (spherical top)
   double *w;				// Set up for system energy levels (LAB)
   w = new double[hs];
   gen_op Holab; 				// these will be in Hertz, ~10**8 in value
   if(abs(level) > 1)				// Needed for higher level computations
     {						// where spectral densities are evaluated
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);				//	This is H is the lab frame
     Holab.Op_base(Ho,1.e-7);			//	It must in the Ho eigenbasis
     if(!Holab.test_EBR())	 		//	and should know its an eigenbasis
       std::cout << "\n\tWarning relax_Dip: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

//	      Dipolar Relaxation Superoperator Computation

   Rijkl(LOp,sys,Ho,w,xis,xis, T,T,tau,type,level);
   delete [] w;
   return LOp;
   }

// ------------- Dipole-Dipole, With An Applied RF-Field ----------------

void RDDrf(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
                 double Wrflab, double* taus, double chi, int type, int level)

	// Input		sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note				Computed in the eigenbasis of Heff

   {
//		    Prepare the Interaction Constants

   matrix xis = xiD(sys);			// Dipolar interaction constants

//	                Prepare the Spin Tensors

   int ns = sys.spins();			// Total number of spins
   int ndip = sys.dipoles();			// Total number of dipoles
   spin_T *T;					// Spin tensors for each dipole
   T = new spin_T[ndip];			// Compiler didn't like spin_T T[ndip]
   int ij=0;
   int i;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xis.get(i,j)))
         {
         T[ij] = T_D(sys,i,j);
         ij++;
	 }

//	               Prepare the Space Tensors

   space_T *A;
   A = new space_T[ndip];			// Compiler didn't like space_T A[ndip]
   ij=0;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xis.get(i,j)))
         {
//         A[ij] = sys.TC(i,j);
         ij++;
	 }

//	           Dipolar Relaxation Superoperator 


   Rrfijkl(LOp, sys, Heff, w, Wrflab, xis, xis, A, A,	// Dipole-dipole interactions
                  T, T, taus, chi, type, level);
   return;
   }


super_op RDDrf(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type, int level)

	// Input		sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian (Hertz)
	//			Wrf   : RF-Field frequecy, lab frame (Hertz)
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff

   {
//   int ns = sys.spins();			// Total number of spins
//   int ndip = sys.dipoles();			// Total number of dipoles
   int hs = sys.HS();				// Total system Hilbert space
   int ls = hs*hs;				// Total system Liouville space
   Heff.set_EBR();				// Insure Heff is in its eigenbasis
   matrix mx(ls, ls, 0.0);			// Construct zero superoperator mx
   super_op LOp(mx, Heff.get_basis());		// Construct zero superoperator

//		    Prepare the Spectral Density Components

   double taus[5];				// Get 5 taus for spectral densities
   taust(taus, sys.taus());
   double chi = chit(sys.taus());		// Get the chi value
   double *w;				// Set up system energy levels
   w = new double[hs];
   matrix* J12 = NULL;				// Matrices of spectral densities
   if(abs(level) > 1)				// Needed for higher level computations
     {
     J12 = new matrix[5];	
     Heff.eigvals(w);				// Frequencies w in field rotating frame
     }

//	         Get the Pure Dipolar Relaxation Superoperator

   RDDrf(LOp,sys,Heff,w,Wrf,taus,chi,type,level);	// Set LOp to Dip-Dip R superop
   delete [] w;
   return LOp;
   }


// ---- Dipole-Dipole, No Applied RF-Field, Dynamic Frequency Shift  ----


void RDDds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
                          double* taus, double chi, int type, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note				Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

   matrix xis = xiD(sys);			// Dipolar interaction constants

//	                Prepare the Spin Tensors

   int ns = sys.spins();			// Total number of spins
   int ndip = sys.dipoles();			// Total number of dipoles
   spin_T *T;					// Spin tensors for each dipole
   T = new spin_T[ndip];			// Compiler didn't like spin_T T[ndip]
   int ij=0;
   int i;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xis.get(i,j)))
         {
         T[ij] = T_D(sys,i,j);
         ij++;
	 }

//	               Prepare the Space Tensors

   space_T *A;
   A = new space_T[ndip];			// Compiler didn't like space_T A[ndip]
   ij=0;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xis.get(i,j)))
         {
//         A[ij] = sys.TC(i,j);
         ij++;
	 }

//	           Dipolar Relaxation Superoperator 

   Rijklds(LOp, sys, Ho, w, xis, xis, A, A,	// Dipole-dipole interactions
                  T, T, taus, chi, type, level);
   return;
   }


super_op RDDds(const sys_dynamic& sys, gen_op& Ho, int type, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note				Computed in the eigenbasis of Ho

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
     {						// where spectral densities are evaluated
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);				//	This is H is the lab frame
     Holab.Op_base(Ho,1.e-7);			//	It must in the Ho eigenbasis
     if(!Holab.test_EBR())	 		//	and should know its an eigenbasis
       std::cout << "\n\tWarning relax_Dip: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

//	         Get the Pure Dipolar Relaxation Superoperator

   RDDds(LOp,sys,Ho,w,taus,chi,type,level);	// Set LOp to Dip-Dip R superop
   delete [] w;
   return LOp;
   }


// -- Dipole-Dipole, With An Applied RF-Field, Dynamic Frequency Shift  -


void RDDrfds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double* w,
                 double Wrflab, double* taus, double chi, int type, int level)

	// Input		sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note				Computed in the eigenbasis of Heff

   {
//		    Prepare the Interaction Constants

   matrix xis = xiD(sys);			// Dipolar interaction constants

//	                Prepare the Spin Tensors

   int ns = sys.spins();			// Total number of spins
   int ndip = sys.dipoles();			// Total number of dipoles
   spin_T *T;					// Spin tensors for each dipole
   T = new spin_T[ndip];			// Compiler didn't like spin_T T[ndip]
   int ij=0;
   int i;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xis.get(i,j)))
         {
         T[ij] = T_D(sys,i,j);
         ij++;
	 }

//	               Prepare the Space Tensors

   space_T *A;
   A = new space_T[ndip];			// Compiler didn't like space_T A[ndip]
   ij=0;
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     for(int j=i+1; j<ns; j++)			// Compute all dipolar spin tensors
       if(Re(xis.get(i,j)))
         {
//         A[ij] = sys.TC(i,j);
         ij++;
	 }

//	           Dipolar Relaxation Superoperator 


   Rrfijklds(LOp, sys, Heff, w, Wrflab, xis, xis, A, A,	// Dipole-dipole interactions
                  T, T, taus, chi, type, level);
   return;
   }


super_op RDDrfds(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type, int level)

	// Input		sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian (Hertz)
	//			Wrf   : RF-Field frequecy, lab frame (Hertz)
	// 			type  : Relaxation type to compute
	//				   0 = dipole-dipole AC & CC
	//				   + = dipole-dipole AC
	//				   - = dipole-dipole CC
	//			level : Relaxation treatment level
	// Output		LOp   : Dipolar relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff

   {
//   int ns = sys.spins();			// Total number of spins
//   int ndip = sys.dipoles();			// Total number of dipoles
   int hs = sys.HS();				// Total system Hilbert space
   int ls = hs*hs;				// Total system Liouville space
   Heff.set_EBR();				// Insure Heff is in its eigenbasis
   matrix mx(ls, ls, 0.0);			// Construct zero superoperator mx
   super_op LOp(mx, Heff.get_basis());		// Construct zero superoperator

//		    Prepare the Spectral Density Components

   double taus[5];				// Get 5 taus for spectral densities
   taust(taus, sys.taus());
   double chi = chit(sys.taus());		// Get the chi value
   double *w;				// Set up system energy levels
   w = new double[hs];
   matrix* J12 = NULL;				// Matrices of spectral densities
   if(abs(level) > 1)				// Needed for higher level computations
     {
     J12 = new matrix[5];	
     Heff.eigvals(w);				// Frequencies w in field rotating frame
     }

//	         Get the Pure Dipolar Relaxation Superoperator

   RDDrfds(LOp,sys,Heff,w,Wrf,taus,chi,type,level);	// Set LOp to Dip-Dip R superop
   delete [] w;
   return LOp;
   }


// ______________________________________________________________________
// *********** Dipole-Dipole Interaction & Coupling Constants ***********
// ______________________________________________________________________

matrix xiD(const sys_dynamic& dsys, double cutoff)


	// Input		cvec  : A coordinate vector(this)
	// 			cutoff: Cutoff value to zero xis
	// Return		dximx : A matrix of dipolar interaction
	//				constants (unit rad/sec)
	//					  -1      -2
	// Note			      : 1T = 1 J-C  -sec-m	
	//
	//		          1/2
	//	            [6*pi]      mu.
	//	       -2 * |----|    * --- * hbar * gamma  * gamma
	//	  D         [ 5  ]      4pi               i        j
	//	xi   = _____________________________________________
	//	  ij		         3
	//			        r
	//			         ij

  {
  matrix dximx(dsys.spins(), dsys.spins(), complex0);
  int npts = dsys.spins();		// Number of spins in the system
  double xiij = 0;
  for(int i=0; i<npts-1; i++)		// Loop over all spin pairs
    for(int j=i+1; j<npts; j++)
      {
      xiij = xiD(dsys, i, j);		// Get the dipolar xi for this pair
      if(fabs(xiij) < cutoff)
        xiij = 0.0;
      dximx.put(xiij, i, j);		// Store in the xi matrix
      dximx.put(xiij, j, i);
      }
  return dximx;
  }


double xiD(const sys_dynamic& dsys, int i, int j)

	// Input		dsys  : A dynamic spin system
	// 			i,j   : Spin indices
	// Return		xi    : Dipolar interaction constant
	//				for spin pair i,j (rad/sec)
	//					  -1      -2
	// Note			      : 1T = 1 J-C  -sec-m	
	//
	//		          1/2
	//	            [6*pi]      mu.
	//	       -2 * |----|    * --- * hbar * gamma  * gamma
	//	  D         [ 5  ]      4pi               i        j
	//	xi   = _____________________________________________
	//	  ij		         3
	//			        r
	//			         ij

  {
  const double a = 1.941625913;		// sqrt[6*pi/5]
  const double K = -2.0*a*mu0d4pi*hbar;
  double xiij = 0;
  if(i != j)
    {
    coord pti = dsys.get(i); 		// spin i coords      -1  -1
    double gi = dsys.gamma(i);		// gamma of spin i sec  -T
    coord ptj = dsys.get(j); 		// spin j coords      -1  -1
    double gj = dsys.gamma(j);		// gamma of spin j sec  -T
    double rij = Rad(pti, ptj);		// dist. between i & j (meters)
    double r3 = rij*rij*rij;
    xiij = K*gi*gj/r3;
    }
  return xiij;
  }


matrix xiD(const spin_system& sys, matrix& dist, int angs, double cutoff)


	// Input		sys   : Spin system
	// 		        dist  : Distances between spins
	// 			angs  : Flag for angstrom distances
	// 			cutoff: Cutoff value to zero xis
	// Return		xis   : Dipolar interaction constant
	//				matrix (unit rad/sec)
	//					  -1      -2
	// Note			      : 1T = 1 J-C  -sec-m	
	//
	//		          1/2
	//	            [6*pi]     mu. 
	//	       -2 * |----|   * --- * hbar * gamma  * gamma
	//	  D         [ 5  ]     4pi               i        j
	//	xi   = ____________________________________________
	//	  ij		         3
	//			        r
	//			         ij

  {
  matrix xis(sys.spins(), sys.spins(), complex0);
  int ns = sys.spins();
  double xiij = 0;
  double gi=0, gj=0, rij=0;
  for(int i=0; i<ns-1; i++)		// Loop over all spin pairs
    {
    gi = sys.gamma(i);
    for(int j=i+1; j<ns; j++)
      {
      gj = sys.gamma(j);
      rij = Re(dist.get(i,j));
      xiij = xiD(gi, gj, rij, angs); 
      if(fabs(xiij) < cutoff)
        xiij = 0.0;
      xis.put(xiij, i, j);
      xis.put(xiij, j, i);
      }
    }
  return xis;
  }


double xiD(double gi, double gj, double rij, int angs)

	// Input		gi    : Gamma of 1st spin
	// Input		gj    : Gamma of 2nd spin
	// Input		rij   : Distance between spins
	// 			angs  : Flag for angstrom distances
	// Return		xiij  : Dipolar interaction constant
	//				(unit rad/sec)
	//					  -1      -2
	// Note			      : 1T = 1 J-C  -sec-m	
	//
	//		          1/2
	//	            [6*pi]     mu. 
	//	       -2 * |----|   * --- * hbar * gamma  * gamma
	//	  D         [ 5  ]     4pi               i        j
	//	xi   = ____________________________________________
	//	  ij		         3
	//			        r
	//			         ij

  {
  const double a = 1.941625913;		// sqrt[6*pi/5]
  if(angs) rij *= 1.0e-10;		// Switch from Angs. to meters
  double r3 = rij*rij*rij;
  double xiij = -2.0*a*mu0d4pi*hbar*gi*gj/r3;
  return xiij;
  }


matrix DCC(const sys_dynamic& dsys)


	// Input		cvec  : A coordinate vector(this)
	// 			angs  : Flag for angstrom distances
	// Return		DCCmx : A matrix of dipolar coupling
	//				constants (unit rad/sec)
	//					  -1      -2
	// Note			      : 1T = 1 J-C  -sec-m	
	//
	//	       mu.
	//	       --- * hbar * gamma  * gamma
	//	       4pi               i        j
	//     DCC   = ____________________________
	//	  ij		  3
	//			 r
	//			  ij

  {
  matrix DCCmx(dsys.spins(), dsys.spins(), complex0);
  int npts = dsys.spins();		// Number of spins in the system
  double DCCij = 0;
  for(int i=0; i<npts-1; i++)		// Loop over all spin pairs
    for(int j=i+1; j<npts; j++)
      {
      DCCij = DCC(dsys, i, j);		// Get spin pair dipolar coupling const.
      DCCmx.put(DCCij, i, j);		// Store in the DCC matrix
      DCCmx.put(DCCij, j, i);
      }
  return DCCmx;
  }


 double DCC(const sys_dynamic& dsys, int i, int j)

	// Input		dsys  : A dynamic spin system
	// 			i,j   : Spin indices
	// Return		DCCij : Dipolar coupling constant
	//				for spin pair i,j (rad/sec)
	//					  -1      -2
	// Note			      : 1T = 1 J-C  -sec-m	
	//
	//
	//	       mu.
	//	       --- * hbar * gamma  * gamma
	//	       4pi               i        j
	//     DCC   = ____________________________
	//	  ij		  3
	//			 r
	//			  ij

  {
  const double K = mu0d4pi*hbar;
  double DCCij = 0;
  if(i != j)
    {
    coord pti = dsys.get(i); 		// spin i coords      -1  -1
    double gi = dsys.gamma(i);		// gamma of spin i sec  -T
    coord ptj = dsys.get(j); 		// spin j coords      -1  -1
    double gj = dsys.gamma(j);		// gamma of spin j sec  -T
    double rij = Rad(pti, ptj);		// dist. between i & j (meters)
    double r3 = rij*rij*rij;
    DCCij = K*gi*gj/r3;
    }
  return DCCij;
  }


 double DCC(double gi, double gj, double rij, int angs)

	// Input		gi    : Gamma of 1st spin
	// Input		gj    : Gamma of 2nd spin
	// Input		rij   : Distance between spins
	// 			angs  : Flag for angstrom distances
	// Return		DCCij : Dipolar coupling constant
	//				(unit rad/sec)
	//					  -1      -2
	// Note			      : 1T = 1 J-C  -sec-m	
	//
	//
	//	       mu.
	//	       --- * hbar * gamma  * gamma
	//	       4pi               i        j
	//     DCC   = ____________________________
	//	  ij		  3
	//			 r
	//			  ij

  {
  if(angs) rij *= 1.0e-10;		// Switch from Angs. to meters
  double r3 = rij*rij*rij;
  double DCCij = mu0d4pi*hbar*gi*gj/r3;
  return DCCij;
  }


// ____________________________________________________________________________
// ____________________________________________________________________________
// ____________________________________________________________________________
//                 Classical Two-Spin Dipolar Relaxation Functions
// ____________________________________________________________________________
// ____________________________________________________________________________

/* These functions are only "classical" in the sense that they apply a two-spin
   approximation.  Any particular spin will relax from a sum over all it's
   dipolar pairings with all other spins.  However, the equations implemented
   are indeed based firmly on quantum mechanics.  Unless indicated otherwise,
   they do assume that relaxation stems from random diffusive motion of spins
   moving as an isotropic top.						     */
 
// ____________________________________________________________________________
//                         Longitudinal Relaxation, T1
// ____________________________________________________________________________

/*
	      We'll Use The Following Two Formulae As Taken From
       Farrar & Becker, "Pulse & F.T. NMR", Academic Press, New York, 1971 
  
For Spin I Relaxed by Unlike Spin S:
------------------------------------
------------------------------------
  
 DDU   1     2   2      2        [ 1               3          3            ]
R   = --- = g * g * hbar * S(S+1)| -- J (w - w ) + - J (w ) + - J (w + w ) |
 1     T     I   S               [ 12  0  I   S    2  1  I    4  2  I   S  ]
  	   1
  
       2 2  _2        tau  [        2               6               12        ]
    = g g * h * S(S+1)____ | ________________ + __________ + ________________ |
       I S               6 |            2   2        2   2              2   2 |
                      15r  | 1 + (w -w ) tau    1 + w tau    1 + (w +w ) tau  |
                           [       I  S              I             I  S       ]
  
  
For Spin I Relaxed by Like Spin:
--------------------------------
--------------------------------
  
  
 DDL   1    3  4      2         [                ]
R   = --- = _ g * hbar * I(I+1) | J (w) + J (2w) |
 1     T    2                   [  1       2     ]
        1
  
           4      2         tau [     1              4      ]
    = 2 * g * hbar * I(I+1) ___ | __________ +  ___________ |
                              6 |      2   2          2   2 |
                            5r  [ 1 + w tau     1 + 4w tau  ]
               
  
Note That For SI Units (which GAMMA uses) We Need 2 mu /(4*pi) Factors Also
                                                      0

Multiplied Into These Equations.  Of Course, We All Know The Conversion

                                      Joule-second
                          1 Telsa = 1 ------------ 2
                                      Coulomb-meter                          */


row_vector R1_DD(const sys_dynamic& sys) 

	// Input	       sys : Spin system
	// Output		R1 : R1 values for spins in the system
	//			     via dipolar relaxation.
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  {
  int ns = sys.spins();			// Get spins in the system
  row_vector R1s(ns);			// Vector of R1 rates
  complex R1;
  for(int i=0; i<ns; i++)		// Loop through all the spins
    {
    R1 = complex(R1_DD(sys,i), 0);	// Get R1 for this spin
    R1s.put(R1,i);			// Store it in the vector 
    }     
  return R1s;
  }


double R1_DD(const sys_dynamic& sys, int i) 

	// Input	       sys : Spin system
	//		         i : Spin index
	// Output		R1 : R1 value for spin i dipolar relaxed by all
	//			     other spins in the system
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  {
  int ns = sys.spins();			// Get spins in the system
  double R1 = 0;			// Assume no relaxation
  for(int j=0; j<ns; j++)		// Loop over all spins
    if(j != i) R1 += R1_DD(sys, i, j);	// Apply two spin approximation
  return R1;
  }


double R1_DD(const sys_dynamic& sys, int i, int j) 

	// Input	       sys : Spin system
	//		         i : Spin index
	//			 j : Spin index
	// Output		R1 : R1 value for spin i dipolar relaxed
	//			     by spin j
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {						//              2 -2 -1
    double gammai = sys.gamma(i); 		//                             -1 -1
    double gammaj = sys.gamma(j);		// Gyromagnetic ratios i&j (sec  T  )
    double r = sys.distance(i,j);		// Distance between spins (meters)
    double fi, fj;
    fi = mu0d4pi*gammai*gammai*hbar/(r*r*r);	// Careful not to overload exponent
    fj = mu0d4pi*gammaj*gammaj*hbar/(r*r*r);
    double wi = sys.lab_shift(i)*pi2;		// Lab frequency spin i (rad/sec)
    double wj = sys.lab_shift(j)*pi2;		// Lab frequency spin j (rad/sec)
    double IjIjp1 = sys.qn(j)*(sys.qn(j)+1.0);	// Ij * (Ij+1)
    double tau = sys.taux();			// Use x-axis tau value only (sec) 
    double factor = fi*fj*IjIjp1*tau;		// Total pre-factor
    double term1, term2, term3;
    double wimwj, wipwj;
    double R1 = 0;

    if(sys.isotope(i) != sys.isotope(j))
      {
      factor *= (1.0/15.0);			// Total pre-factor
      wimwj = wi - wj;
      wipwj = wi + wj;
      term1 = 2.0/(1+(wimwj*wimwj*tau*tau));	// 2/[1+(wi-wj)^2*tau^2] 
      term2 = 6.0/(1+(wi*wi*tau*tau)); 		// 6/[1+wi^2*tau^2]
      term3 = 12.0/(1+(wipwj*wipwj*tau*tau)); 	// 12/[1+(wi+wj)^2*tau^2]
      }
    else
      {
      factor *= (2.0/5.0);			// Total pre-factor
      term1 = 1.0/(1+(wi*wi*tau*tau)); 
      term2 = 4.0/(1+(4.0*wi*wi*tau*tau)); 
      term3 = 0;
      }
    R1 = factor * (term1+term2+term3);
    return R1;
    }


double R1_DD_max(const sys_dynamic& sys, int i) 

        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            R1max : Maximum R1 value of all spins in
        //                           the system of isotope the same as the
        //                           input spin i
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.

    {
    int ns = sys.spins();               // Get spins in the system
    double R1, R1max=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(sys.isotope(i)==sys.isotope(j))// Only consider spins like spin i
        {
        R1 = R1_DD(sys, j);             // R1 spin j, 2 spin approximation
        if(R1 > R1max)                  // Take the largest R1 of these
          R1max = R1;
        }
    return R1max;
    }


double R1_DD_max(const sys_dynamic& sys, const std::string& Iso) 

        // Input               sys : Spin system
        //                     Iso : A string for an isotope type
        // Output            R1max : Maximum R1 value of all spins in
        //                           the system of isotope the same as Iso
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.

    {
    Isotope I(Iso);
    int ns = sys.spins();               // Get spins in the system
    double R1, R1max=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(I==sys.isotope(j))		// Only consider spins of type Iso
        {
        R1 = R1_DD(sys, j);             // R1 spin i, 2 spin approximation
        if(R1 > R1max)                  // Take the largest R1 of these
          R1max = R1;
        }
    return R1max;
    }


double R1_DD_max(const sys_dynamic& sys) 

	// Input	       sys : Spin system
	// Output	     R1max : Maximum R1 value of all spin spins in
	//			     the system under dipolar relaxation
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {
    int ns = sys.spins();		// Get spins in the system
    double R1, R1max = 0;
    for(int i=0; i<ns; i++)		// Loop over all spins
      {
      R1 = R1_DD(sys, i);		// Two spin approximation
      if(R1 > R1max)
        R1max = R1;
      }
    return R1max;
    }


row_vector T1_DD(const sys_dynamic& sys) 

	// Input	       sys : Spin system
	// Output		T1 : T1 values for spins in the system
	//			     via dipolar relaxation.
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {
    int ns = sys.spins();
    row_vector T1s(ns);			// Vector of T1 times
    row_vector R1s = R1_DD(sys);	// Vector of R1 rates
    double T1;
    for(int i=0; i<ns; i++)		// Loop through all the spins
      {
      T1 = 1.0/Re(R1s.get(i));		// Get T1 for this spin
      T1s.put(T1,i);			// Store it in the vector 
      }     
    return T1s;
    }


double T1_DD(const sys_dynamic& sys, int i) 

	// Input	       sys : Spin system
	//		         i : Spin index
	// Output		T1 : T1 value for spin i dipolar relaxed by all
	//			     other spins in the system
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    { return 1.0/R1_DD(sys,i); }


double T1_DD(const sys_dynamic& sys, int i, int j) 

	// Input	       sys : Spin system
	//		         i : Spin index
	//			 j : Spin index
	// Output		T1 : T1 value for spin i dipolar relaxed
	//			     by spin j
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    { return 1.0/R1_DD(sys,i,j); }


double T1_DD_max(const sys_dynamic& sys, int i) 

        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            T1max : Maximum T1 value of all spins in
        //                           the system of isotope the same as the
        //                           input spin i
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.

    {
    int ns = sys.spins();               // Get spins in the system
    double T1, T1max=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(sys.isotope(i)==sys.isotope(j))// Only consider spins like spin i
        {
        T1 = T1_DD(sys, j);             // T1 spin j, 2 spin approximation
        if(T1 > T1max)                  // Take the largest T1 of these
          T1max = T1;
        }
    return T1max;
    }


double T1_DD_max(const sys_dynamic& sys, const std::string& Iso) 

        // Input               sys : Spin system
        //                     Iso : A string for an isotope type
        // Output            T1max : Maximum T1 value of all spins in
        //                           the system of isotope the same as Iso
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.

    {
    Isotope I(Iso);
    int ns = sys.spins();               // Get spins in the system
    double T1, T1max=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(I==sys.isotope(j))		// Only consider spins of type Iso
        {
        T1 = T1_DD(sys, j);             // T1 spin i, 2 spin approximation
        if(T1 > T1max)                  // Take the largest T1 of these
          T1max = T1;
        }
    return T1max;
    }


double T1_DD_max(const sys_dynamic& sys) 

	// Input	       sys : Spin system
	// Output	     T1max : Maximum T1 value of all spin spins in
	//			     the system under dipolar relaxation
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {
    int ns = sys.spins();		// Get spins in the system
    double T1, T1max = 0;
    for(int i=0; i<ns; i++)		// Loop over all spins
      {
      T1 = T1_DD(sys, i);		// Two spin approximation
      if(T1 > T1max)
        T1max = T1;
      }
    return T1max;
    }


// ------------------- Transverse Relaxation, T2 ------------------------

//  Farrar & Becker, "Pulse & F.T. NMR", Academic Press, New York, 1971 
//
//  For Spin I Relaxed by Unlike, Uncoupled Spin S
//        Eq. 4.19 page 55 & Eq. 4.9 page 51
//
//  DDU   1     2   2      2         [ 1          1           3          3          3            ]
// R   = --- = g * g * hbar * S(S+1) | _ J (0) + __ J (w  ) + _ J (w ) + _ J (w ) + _ J (w + w ) |
//  2     T     I   S                [ 6  0      24  0  IS    4  1  I    2  1  S    8  2  I   S  ]
//	   2
//
//
//              2   2      2         tau  [          1            3            6               6          ]
//           = g * g * hbar * S(S+1) ____ | 4 + ___________ + __________ + __________ + _________________ |
//              I   S                   6 |          2    2        2   2        2   2               2   2 |
//                                   15r  |     1 + w  tau    1 + w tau    1 + w tau    1 + (w + w ) tau  |
//             			          [	     IS            I            S             I   S       ]
//
//
//  For Spin I Relaxed by Unlike, Coupled Spin S
//
//  DDJ   1     2   2      2         tau  [          1            3              3              6         ]
// R   = --- = g * g * hbar * S(S+1) ____ | 4 + ___________ + __________ + __________ + _________________ |
//  2     T     I   S                   6 |          2    2        2   2        2   2               2   2 |
//         2                         15r  |     1 + w  tau    1 + w tau    1 + w tau    1 + (w + w ) tau  |
//             				  [	     IS            I            S             I           ]
//             
//
//  For Spin I Relaxed by Like Spin
//  Eq. 4.16 page 55,  Eq. 4.9 page 51, and Eq. 4.21 page 56
//
//  DDL   1     4      2         [ 3         15          3         ]
// R   = --- = g * hbar * I(I+1) | _ J (0) + __ J (w ) + _ J (2w ) |
//  2     T                      [ 8  0       4  1  I    8  2   I  ]
//	   2
//
//              4      2         tau [           5            2     ]
//           = g * hbar * I(I+1) ___ | 3 + __________ + ___________ |
//                                 6 |          2   2         2   2 |
//                               5r  [     1 + w tau    1 + 4w tau  ]
//             
//
// FOR SI UNITS (WHICH GAMMA USES) WE NEED 2 mu /(4*pi) FACTORS ALSO IN FRONT!
//
//

row_vector R2_DD(const sys_dynamic& sys) 

	// Input	       sys : Spin system
	// Output		R2 : R2 values for spins in the system
	//			     via dipolar relaxation.
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {
    int ns = sys.spins();		// Get spins in the system
    row_vector R2s(ns);			// Vector of R2 times
    complex R2;
    for(int i=0; i<ns; i++)		// Loop through all the spins
      {
      R2 = complex(R2_DD(sys,i), 0);	// Get R2 for this spin
      R2s.put(R2,i);			// Store it in the vector 
      }     
    return R2s;
    }


double R2_DD(const sys_dynamic& sys, int i) 

	// Input	       sys : Spin system
	//		         i : Spin index
	// Output		R2 : R2 value for spin i dipolar relaxed by all
	//			     other spins in the system
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {
    int ns = sys.spins();		// Get spins in the system
    double R2 = 0;
    for(int j=0; j<ns; j++)		// Loop over all spins
      if(j != i)
        R2 += R2_DD(sys, i, j);		// Two spin approximation
    return R2;
    }


double R2_DD(const sys_dynamic& sys, int i, int j) 

	// Input	       sys : Spin system
	//		         i : Spin index
	//			 j : Spin index
	// Output		R2 : R2 value for spin i dipolar relaxed
	//			     by spin j
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {						//              2 -2 -1
    const double mu0d4pi = 1.0e-7;		// mu0/4p (J-sec C  m  )
    const double pi2 = 6.283185307;		// 2*pi
    const double hbar = 1.05459e-34;		// hbar (J-sec)
    double gammai = sys.gamma(i); 		//                             -1 -1
    double gammaj = sys.gamma(j);		// Gyromagnetic ratios i&j (sec  T  )
    double r = sys.distance(i,j);		// Distance between spins (meters)
    double fi, fj;
    fi = mu0d4pi*gammai*gammai*hbar/(r*r*r);	// Careful not to overload exponent
    fj = mu0d4pi*gammaj*gammaj*hbar/(r*r*r);
    double wi = sys.lab_shift(i)*pi2;		// Lab frequency spin i (rad/sec)
    double wj = sys.lab_shift(j)*pi2;		// Lab frequency spin j (rad/sec)
    double IjIjp1 = sys.qn(j)*(sys.qn(j)+1.0);	// Ij * (Ij+1)
    double tau = sys.taux();			// Use x-axis tau value only (sec) 
    double factor = fi*fj*IjIjp1*tau;		// Total pre-factor
    double term1, term2, term3;
    double term4=0, term5=0;
    double wimwj, wipwj;
    double R2 = 0;
    if(sys.isotope(i) == sys.isotope(j))	// "DDL", Like spins
      {
      factor /= 5.0;
      term1 = 3.0;
      term2 = 5.0/(1+wi*wi*tau*tau); 
      term3 = 2.0/(1+4.0*wi*wi*tau*tau); 
      }
    else					// "DDU", Unlike spins
      {
      factor /= 15.0;
      wimwj = wi - wj;
      wipwj = wi + wj;
      term1 = 4.0;
      term2 = 1.0/(1+wimwj*wimwj*tau*tau); 	// 1/[1+(wi-wj)^2*tau^2]
      term3 = 3.0/(1+wi*wi*tau*tau); 		// 3/[1+(wi)^2*tau^2]
      term4 = 6.0/(1+wj*wj*tau*tau); 		// 6/(1+(wj)^2*tau^2]
      term5 = 6.0/(1+wipwj*wipwj*tau*tau); 	// 6/(1+(wi+wj)^2*tau^2]
      if(sys.J(i,j))				// For "DDJ", spins coupled
        term4 /= 2.0;
      }
    R2 = factor*(term1+term2+term3+term4+term5);
    return R2;
    }


double R2_DD_max(const sys_dynamic& sys, int i)
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            R2max : Maximum R2 (under dipolar relaxation) value
        //                           of all spins in the system of isotope type 
        //                           the same as the input spin i
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    int ns = sys.spins();               // Get spins in the system
    double R2, R2max = 0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(sys.isotope(i)==sys.isotope(j))// Only consider spins like spin i
        {
        R2 = R2_DD(sys, j);             // Two spin approximation
        if(R2 > R2max)
          R2max = R2;
        }
    return R2max;
    }


double R2_DD_max(const sys_dynamic& sys, const std::string& Iso)
 
        // Input               sys : Spin system
        //                     Iso : A string for an isotope type
        // Output            R2max : Maximum R2 (under dipolar relaxation) value
        //                           of all spins in the system of isotope Iso
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    Isotope I(Iso);
    int ns = sys.spins();               // Get spins in the system
    double R2, R2max = 0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(I == sys.isotope(j))		// Only consider spins like Iso
        {
        R2 = R2_DD(sys, j);             // Two spin approximation
        if(R2 > R2max)
          R2max = R2;
        }
    return R2max;
    }
 

double R2_DD_max(const sys_dynamic& sys) 

	// Input	       sys : Spin system
	// Output	     R2max : Maximum R2 value of all spin spins in
	//			     the system under dipolar relaxation
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {
    int ns = sys.spins();		// Get spins in the system
    double R2, R2max = 0;
    for(int i=0; i<ns; i++)		// Loop over all spins
      {
      R2 = R2_DD(sys, i);		// Two spin approximation
      if(R2 > R2max)
        R2max = R2;
      }
    return R2max;
    }


row_vector T2_DD(const sys_dynamic& sys) 

	// Input	       sys : Spin system
	// Output		T2 : T2 values for spins in the system
	//			     via dipolar relaxation.
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {
    int ns = sys.spins();
    row_vector T2s(ns);			// Vector of T2 times
    row_vector R2s = R2_DD(sys);	// Vector of R2 rates
    double T2;
    for(int i=0; i<ns; i++)		// Loop through all the spins
      {
      T2 = 1.0/Re(R2s.get(i));		// Get T2 for this spin
      T2s.put(T2,i);			// Store it in the vector 
      }     
    return T2s;
    }


double T2_DD(const sys_dynamic& sys, int i) 

	// Input	       sys : Spin system
	//		         i : Spin index
	// Output		T2 : T2 value for spin i dipolar relaxed by all
	//			     other spins in the system
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    { return (1.0/R2_DD(sys,i)); }


double T2_DD(const sys_dynamic& sys, int i, int j) 

	// Input	       sys : Spin system
	//		         i : Spin index
	//			 j : Spin index
	// Output		T2 : T2 value for spin i dipolar relaxed
	//			     by spin j
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    { return (1.0/R2_DD(sys,i,j)); }


double T2_DD_max(const sys_dynamic& sys, int i)
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            T2max : Maximum T2 (under dipolar relaxation) value
        //                           of all spins in the system of isotope type 
        //                           the same as the input spin i
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    int ns = sys.spins();               // Get spins in the system
    double T2, T2max = 0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(sys.isotope(i)==sys.isotope(j))// Only consider spins like spin i
        {
        T2 = T2_DD(sys, j);             // Two spin approximation
        if(T2 > T2max)
          T2max = T2;
        }
    return T2max;
    }


double T2_DD_max(const sys_dynamic& sys, const std::string& Iso)
 
        // Input               sys : Spin system
        //                     Iso : A string for an isotope type
        // Output            T2max : Maximum T2 (under dipolar relaxation) value
        //                           of all spins in the system of isotope Iso
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    Isotope I(Iso);
    int ns = sys.spins();               // Get spins in the system
    double T2, T2max = 0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(I == sys.isotope(j))		// Only consider spins like Iso
        {
        T2 = T2_DD(sys, j);             // Two spin approximation
        if(T2 > T2max)
          T2max = T2;
        }
    return T2max;
    }


double T2_DD_max(const sys_dynamic& sys) 

	// Input	       sys : Spin system
	// Output	     T2max : Maximum T2 value of all spin spins in
	//			     the system under dipolar relaxation
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {
    int ns = sys.spins();		// Get spins in the system
    double T2, T2max = 0;
    for(int i=0; i<ns; i++)		// Loop over all spins
      {
      T2 = T2_DD(sys, i);		// Two spin approximation
      if(T2 > T2max)
        T2max = T2;
      }
    return T2max;
    }


// ----------------------- Dipolar Linewidths ---------------------------
//
//                                 DD
//                                R
//                           DD    2      1.0
//                         LW   = --- = -------
//                           hh    pi    DD 
//                                     	T  * pi
//                                       2
//

row_vector LWhh_DD(const sys_dynamic& sys) 

	// Input	       sys : Spin system
	// Output		LW : Expected linewidths at half-height
 	//			     for spins in the system via dipolar relaxation.
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {
    int ns = sys.spins();
    row_vector LWs(ns);			// Vector of linewidths
    row_vector R2s = R2_DD(sys);	// Vector of R2 rates
    double R2;
    for(int i=0; i<ns; i++)		// Loop through all the spins
      {
      R2 = Re(R2s.get(i));		// Get R2 for this spin
      LWs.put(R2/3.14159,i);		// Store it in the vector 
      }     
    return LWs;
    }


double LWhh_DD(const sys_dynamic& sys, int i) 

	// Input	       sys : Spin system
	//		         i : Spin index
	// Output		LW : Expected linewidth at half-height
	//			     for spin i dipolar relaxed by the system
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    { return (R2_DD(sys,i)/3.14159); }


double LWhh_DD(const sys_dynamic& sys, int i, int j) 

	// Input	       sys : Spin system
	//		         i : Spin index
	//			 j : Spin index
	// Output		LW : Expected linewidth at half-height
	//			     for spin i dipolar relaxed by spin j
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    { return (R2_DD(sys,i,j)/3.14159); }


double LWhh_DD_max(const sys_dynamic& sys, int i)

        // Input               sys : Spin system
        //                       i : A spin index       
        // Output               LW : Maximum linewidth at half-height
        //                           under dipolar relaxation over all
        //                           spins in the system of isotope the
        //                           same as the input spin i
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.

  { return (R2_DD_max(sys, i)/3.14159); }


double LWhh_DD_max(const sys_dynamic& sys, const std::string& Iso)

        // Input               sys : Spin system
        //                     Iso : A spin isotope label
        // Output               LW : Maximum linewidth at half-height
        //                           under dipolar relaxation over all
        //                           spins in the system of isotope Iso
        // Note                    : This routine assumes a two spin approximation
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.

  { return (R2_DD_max(sys, Iso)/3.14159); }


  double LWhh_DD_max(const sys_dynamic& sys) 

	// Input 	       sys : Spin system
	// Output		LW : Maximum linewidth at half-height
	//			     over all spins under dipolar relaxation

  { return (R2_DD_max(sys)/3.14159); }


// ----------------- Nuclear Overhauser Enhancements --------------------
//
//  NOE Enhancement For Spin Relaxed by Spin S

//                [    -1                 6         ]
//                | ___________ + _________________ |
//             g  |      2    2               2   2 |
//              I | 1 + w  tau    1 + (w + w ) tau  |
//                [	   IS           S   I       ]		   
// NOE = ___________________________________________________  = rho
//          [     1             3                6         ]       NOE
//          | ___________ + __________ + _________________ |
//       g  |      2    2        2   2               2   2 |
//        S | 1 + w  tau    1 + w tau    1 + (w + w ) tau  |
//          [	   IS            I             I   S       ]
//             
//             
//                       eta    = 1 + rho
//                          NOE          NOE
//
//

  double NOE(const sys_dynamic& sys, int i, int j, double eta) 

	// Input	       sys : Spin system
	//		         i : Spin index
	//			 j : Spin index
	//		       eta : flag for % enhancement
	// Output	       NOE : Expected NOE for for spin i due to
	//			     dipolar relaxation by spin j
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {
    double gammai = sys.gamma(i); 		//                             -1 -1
    double gammaj = sys.gamma(j);		// Gyromagnetic ratios i&j (sec  T  )
    double wi = sys.lab_shift(i)*pi2;		// Lab frequency spin i (rad/sec)
    double wj = sys.lab_shift(j)*pi2;		// Lab frequency spin j (rad/sec)
    double tau = sys.taux();			// Use x-axis tau value only (sec) 
    double wimwj = wi - wj;
    double wipwj = wi + wj;
    double term1 = 1.0/(1+wimwj*wimwj*tau*tau);
    double term2 = 1.0/(1+wi*wi*tau*tau);
    double term3 = 1.0/(1+wipwj*wipwj*tau*tau);
    double NUM = gammaj*(6.0*term3 - term1);
    double DEN = gammai*(term1+3.0*term2+6*term3);
    double NOE = NUM/DEN;
    if(eta)
      NOE += 1.0; 
    return NOE;
    }


// ----------- Multiple Quantum Transition Relaxation Rates -------------

//
//  ZQT - For Spin I Relaxed by Unlike, Coupled Spin S
//
//  DDJ   1     2   2      2         tau  [     2             3            3      ]
// R   = --- = g * g * hbar * S(S+1) ____ | ___________ + __________ + __________ |
//  2     T     I   S                   6 |      2    2        2   2        2   2 |
//         2                         15r  | 1 + w  tau    1 + w tau    1 + w tau  |
//             				  [	 IS            I            S     ]
//             
//
//  SQT - For Spin I Relaxed by Unlike, Coupled Spin S
//
//  DDJ   1     2   2      2         tau  [          1            3              3              6         ]
// R   = --- = g * g * hbar * S(S+1) ____ | 4 + ___________ + __________ + __________ + _________________ |
//  2     T     I   S                   6 |          2    2        2   2        2   2               2   2 |
//         2                         15r  |     1 + w  tau    1 + w tau    1 + w tau    1 + (w + w ) tau  |
//             				  [	     IS            I            S             I           ]
//             
//
//  DQT - For Spin I Relaxed by Unlike, Coupled Spin S
//
//  DDJ   1     2   2      2         tau  [     3              3              12        ]
// R   = --- = g * g * hbar * S(S+1) ____ | __________ + __________ + _________________ |
//  2     T     I   S                   6 |      2   2        2   2               2   2 |
//         2                         15r  | 1 + w tau    1 + w tau    1 + (w + w ) tau  |
//             				  [	 I            S             I   S	]
//             
//
// FOR SI UNITS (WHICH GAMMA USES) WE NEED 2 mu /(4*pi) FACTORS ALSO IN FRONT!
//
//

  row_vector R2_DDMQT(const sys_dynamic& sys, int MQC) 

	// Input	       sys : Spin system
	//		       MQC : Multiple quantum coherence
	// Output		R2 : R2 values for spins in the system
	//			     via dipolar relaxation.
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {
    int ns = sys.spins();		// Get spins in the system
    row_vector R2s(ns);			// Vector of R2 times
    complex R2;
    for(int i=0; i<ns; i++)		// Loop through all the spins
      {
      R2 = complex(R2_DDMQT(sys,MQC,i), 0);	// Get R2 for this spin
      R2s.put(R2,i);			// Store it in the vector 
      }     
    return R2s;
    }


  double R2_DDMQT(const sys_dynamic& sys, int MQC, int i) 

	// Input	       sys : Spin system
	//		       MQC : Multiple quantum coherence
	//		         i : Spin index
	// Output		R2 : R2 value for spin i dipolar relaxed by all
	//			     other spins in the system
	// Note			   : This routine assumes a two spin approximation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {
    int ns = sys.spins();		// Get spins in the system
    double R2 = 0;
    for(int j=0; j<ns; j++)		// Loop over all spins
      if(j != i)
        R2 += R2_DDMQT(sys, MQC, i, j);	// Two spin approximation
    return R2;
    }


  double R2_DDMQT(const sys_dynamic& sys, int MQC, int i, int j) 

	// Input	       sys : Spin system
	//		       MQC : Multiple quantum coherence
	//		         i : Spin index
	//			 j : Spin index
	// Output		R2 : R2 value for MQC involving spin i dipolar relaxed
	//			     by spin j
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

    {
    double R2 = 0;
    if(MQC == 1)
      R2 = R2_DD(sys, i, j);
    else if(sys.J(i,j) == 0 || abs(MQC) > 2)
      R2 = 0;
    else
      {
      double gammai = sys.gamma(i); 		//                             -1 -1
      double gammaj = sys.gamma(j);		// Gyromagnetic ratios i&j (sec  T  )
      double r = sys.distance(i,j);		// Distance between spins (meters)
      double fi, fj;
      fi = mu0d4pi*gammai*gammai*hbar/(r*r*r);	// Careful not to overload exponent
      fj = mu0d4pi*gammaj*gammaj*hbar/(r*r*r);
      double wi = sys.lab_shift(i)*pi2;		// Lab frequency spin i (rad/sec)
      double wj = sys.lab_shift(j)*pi2;		// Lab frequency spin j (rad/sec)
      double IjIjp1 = sys.qn(j)*(sys.qn(j)+1.0);// Ij * (Ij+1)
      double tau = sys.taux();			// Use x-axis tau value only (sec) 
      double factor = fi*fj*IjIjp1*tau/15.0;	// Total pre-factor
      double term1, term2, term3, term4;
      double wimwj, wipwj;
      term2 = 3.0/(1+wi*wi*tau*tau); 
      term3 = 3.0/(1+wj*wj*tau*tau); 
      if(MQC == 0)
        {
        wimwj = wi - wj;
        term1 = 2.0/(1+wimwj*wimwj*tau*tau); 
        term4 = 0.0;
        }
      else
        {
        wipwj = wi + wj;
        term1 = 0.0;
        term4 = 12.0/(1+wipwj*wipwj*tau*tau); 
        }  
      R2 = factor*(term1+term2+term3+term4);
      }  
    return R2;
    }



#endif /* __RELAX_DIP_CC__ */
