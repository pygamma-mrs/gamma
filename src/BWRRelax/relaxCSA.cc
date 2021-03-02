/* relaxCSA.cc ***********************************************************
**									**
**	                        G A M M A				**
**								 	**
**	CSA Relaxation                           Implementation		**
**							 		**
**	Copyright (c) 1991, 1992, 1993				 	**
**	Scott A. Smith						 	**
**									**
**	Eidgenoessische Technische Hochschule			 	**
**	Labor fuer physikalische Chemie				 	**
**	8092 Zurich / Switzerland			 		**
**									**
**	University of California, Santa Barbara				**
**	Department of Chemistry						**
**	Santa Barbara CA. 93106 USA					**
**								 	**
**      $Header: $
**									**
*************************************************************************/

/*************************************************************************
**									**
** Description							 	**
**									**
** The following functions provide easy access to several CSA		**
** relaxation related functions.					**
**									**
*************************************************************************/

#ifndef _relax_CSA_cc_			// Is this file already included?
#define _relax_CSA_cc_ 1		// If no, then remember that it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation         		// This is the implementation
#endif


#include <BWRRelax/relaxCSA.h>		// Include header file
#include <Level1/nmr_tensor.h>		// Include common spin tensors
#include <LSLib/sys_dynamic.h>		// Include anisotropic systems
#include <stdlib.h>

const double K = 1.941625913;		// sqrt[6*pi/5]
const double pi2 = 6.283185307;		// 2 * pi

// ____________________________________________________________________________
// ****************** CSA-CSA Relaxation Superoperators ***********************
// ____________________________________________________________________________

void RCC(super_op& LOp, const sys_dynamic& sys, gen_op& H, double*w,
                                 double* taus, double chi, int type, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			H     : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of H

   {
//		    Prepare the Interaction Constants

  matrix Xis = xiCSA(sys);			// Matrix of Xi's for CSA

//			Prepare the Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *T;					// Spin tensors for each spin
   T = new spin_T[ns];				// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)				// Set them all to CSA spin tensors
     if(Re(Xis.get(i,i))) T[i] = T_CS2(sys,i);

//			Prepare the Spatial Tensors

   space_T *A; 					// Spatial tensors for each spin
   A = new space_T[ns];				// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them all to CSA space tensors
     if(Re(Xis.get(i,i))) A[i] = sys.TC(i);

//		    Perform Computation Over All Spins

   Rij(LOp, sys, H, w, Xis, Xis, A, A,		// All CSA parameters
                  T, T, taus, chi, type, level);
   return;
   }


super_op RCC(const sys_dynamic& sys, gen_op& H, int type, int level)

	// Input		sys   : Dynamic spin system
	//			H     : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of H

   {
   int hs = sys.HS();				// Total system Hilbert space
   int ls = hs*hs;				// Total system Liouville space
   H.set_EBR();				        // Insure H is in its eigenbasis
   matrix mx(ls, ls, 0.0);			// Construct zero superoperator
   super_op LOp(mx, H.get_basis());

   double taus[5];				// Get 5 taus for spectral densities
   taust(taus, sys.taus());
   double chi = chit(sys.taus());		// Get the system chi value
   double *w;				// Set up for system energy levels (LAB)
   w = new double[hs];
   gen_op Hlab;				// these will be in Hertz, ~10**8 in value
   if(abs(level) > 1)				// Used for higher level calculations
     {						// where spectral densities are evaluated
     Hlab = Hcs_lab(sys);
     Hlab += HJ(sys);		//	This is H in the lab frame
     Hlab.Op_base(H, 1.e-7);			//	It must be in the H eigenbasis
     if(!Hlab.test_EBR())
       std::cout << "\n\tWarning relax_CSA: "
            << " Unable to Obtain Proper H(lab) Eigenbasis";
     Hlab.eigvals(w);
     }


   RCC(LOp,sys,H,w,taus,chi,type,level);	// Set LOp to CSA-CSA relaxation superop
   delete [] w;
   return LOp;
   }

// ---------------- CSA-CSA, With An Applied RF-Field -------------------

void RCCrf(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int type, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			Wrflab: Field freqency in the lab frame (Hz)
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

   matrix Xis = xiCSA(sys);			// Matrix of Xi's for CSA

//			Prepare the Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *T;					// Spin tensors for each spin
   T = new spin_T[ns];				// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)				// Set them all to CSA spin tensors
     if(Re(Xis.get(i,i))) T[i] = T_CS2(sys,i);

//			Prepare the Spatial Tensors

   space_T *A; 					// Spatial tensors for each spin
   A = new space_T[ns];				// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them all to CSA space tensors
     if(Re(Xis.get(i,i)))
       A[i] = sys.TC(i);

//		    Perform Computation Over All Spins

   Rrfij(LOp, sys, Heff, w, Wrflab, Xis, Xis, A, A,	// All CSA parameters
                  T, T, taus, chi, type, level);
   return;
   }


super_op RCCrf(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type, int level)

	// Input                sys   : Dynamic spin system
	//                      Heff  : Effective Hamiltonian (Hertz)
	//                      Wrf   : RF-Field frequecy, lab frame (Hertz)
	//                      phi   : Field phase angle (Degrees)
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff

   {
   int hs = sys.HS();				// Total system Hilbert space
   int ls = hs*hs;				// Total system Liouville space
   Heff.set_EBR();				// Insure Heff is in its eigenbasis
   matrix mx(ls, ls, 0.0);			// Construct zero superoperator
   super_op LOp(mx, Heff.get_basis());

   double taus[5];				// Get 5 taus for spectral densities
   taust(taus, sys.taus());
   double chi = chit(sys.taus());		// Get the system chi value
   double *w;				     // Set up for system energy levels (LAB)
   w = new double[hs];
   matrix* J12 = NULL;				// Matrices for spectral densities
   if(abs(level) > 1)				// Used for higher level calculations
     {
     J12 = new matrix[5];
     Heff.eigvals(w);
     }

   RCCrf(LOp,sys,Heff,w,Wrf,taus,chi,type,level);	// Set LOp to CSA-CSA relaxation superop
   delete [] w;
   return LOp;
   }


// ----- CSA-CSA, No Applied RF-Field, Dynamic Frequency Shift  ------


void RCCds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int type, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

   matrix Xi = xiCSA(sys);			// Matrix of Xi's for CSA

//			Prepare the Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *T;					// Spin tensors for each spin
   T = new spin_T[ns];				// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)			// Set them to quadrupolar spin tensors
     if(Re(Xi.get(i,i)))
       T[i] = T_CS2(sys,i);

//			Prepare the Spatial Tensors

   space_T *A; 					// Spatial tensors for each spin
   A = new space_T[ns];				// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them to quadrupolar space tensors
     if(Re(Xi.get(i,i)))
       A[i] = sys.TC(i);

//		    Perform Computation Over All Spin Pairs

   Rijds(LOp, sys, Ho, w, Xi, Xi, A, A,		// All CSA parameters
               T, T, taus, chi, type, level);
   return;
   }


super_op RCCds(const sys_dynamic& sys, gen_op& Ho, int type, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
   int hs = sys.HS();				// Total system Hilbert space
   int ls = hs*hs;				// Total system Liouville space
   Ho.set_EBR();				// Insure Ho is in its eigenbasis
   matrix mx(ls, ls, 0.0);			// Construct zero superoperator
   super_op LOp(mx, Ho.get_basis());

   double taus[5];				// Get 5 taus for spectral densities
   taust(taus, sys.taus());
   double chi = chit(sys.taus());		// Get the system chi value
   double *w;				// Set up for system energy levels (LAB)
   w = new double[hs];
   gen_op Holab;				// these will be in Hertz, ~10**8 in value
   if(abs(level) > 1)				// Used for higher level calculations
     {
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);
     Holab.Op_base(Ho, 1.e-7);			//	It must be in the Ho eigenbasis
     if(!Holab.test_EBR())
       std::cout << "\n\tWarning relax_CSA: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

   RCCds(LOp, sys, Ho, w, taus, chi, type, level);
   delete [] w;
   return LOp;
   }


// ---- CSA-CSA, With An Applied RF-Field, Dynamic Frequency Shift  -----


void RCCrfds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int type, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			Wrflab: Field freqency in the lab frame (Hz)
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

   matrix Xis = xiCSA(sys);			// Matrix of Xi's for CSA

//			Prepare the Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *T;					// Spin tensors for each spin
   T = new spin_T[ns];				// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)			// Set them all to CSA spin tensors
     if(Re(Xis.get(i,i)))
       T[i] = T_CS2(sys,i);

//			Prepare the Spatial Tensors

   space_T *A; 					// Spatial tensors for each spin
   A = new space_T[ns];				// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them all to CSA space tensors
     if(Re(Xis.get(i,i)))
       A[i] = sys.TC(i);

//		    Perform Computation Over All Spins

   Rrfijds(LOp, sys, Heff, w, Wrflab, Xis, Xis, A, A,	// All CSA parameters
                  T, T, taus, chi, type, level);
   return;
   }


  super_op RCCrfds(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type, int level)

	// Input                sys   : Dynamic spin system
	//                      Heff  : Effective Hamiltonian (Hertz)
	//                      Wrf   : RF-Field frequecy, lab frame (Hertz)
	//                      phi   : Field phase angle (Degrees)
	// 			type  : Relaxation type to compute
	//				   0 = CSA-CSA AC & CC
	//				   + = CSA-CSA AC
	//				   - = CSA-CSA CC
	//			level : Relaxation treatment level
	// Output		LOp   : CSA relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff

   {
   int hs = sys.HS();				// Total system Hilbert space
   int ls = hs*hs;				// Total system Liouville space
   Heff.set_EBR();				// Insure Heff is in its eigenbasis
   matrix mx(ls, ls, 0.0);			// Construct zero superoperator
   super_op LOp(mx, Heff.get_basis());

   double taus[5];				// Get 5 taus for spectral densities
   taust(taus, sys.taus());
   double chi = chit(sys.taus());		// Get the system chi value
   double *w;				// Set up for system energy levels (LAB)
   w = new double[hs];
   matrix* J12 = NULL;				// Matrices for spectral densities
   if(abs(level) > 1)				// Used for higher level calculations
     {
     J12 = new matrix[5];
     Heff.eigvals(w);
     }

   RCCrfds(LOp,sys,Heff,w,Wrf,taus,chi,type,level);	// Set LOp to CSA-CSA relaxation superop
   delete [] w;
   return LOp;
   }


// ______________________________________________________________________
// ********************* CLASSICAL CSA RELAXATION ***********************
// ______________________________________________________________________


// ------------------ Longitudinal Relaxation, T1 -----------------------

// Ferrar & Becker, "Pulse and Fourier Transform NMR", Academic Press, New York, 1971 
//                         Eq. 4.28 page 59 and Eq. 4.9 page 51.
//
//
//  I    1            2    2         2  tau [       2         ]
// R  = --- = (gamma ) (B ) (s  - s )   ___ | _______________ |
//  1    T          I    0    ||   |        |         2     2 |
//        1                   ||  ---   15  | 1 + (w ) (tau)  |
//             				    [	    I         ]
//
//
//                    2     I  2   tau  [       3         ]
//          = (Omega )  (del  )   _____ | _______________ |
//                  I       zz          |       2       2 |
//                                 10   | 1 + (w ) (tau)  |
//             				[       I         ]
//             
//                                3
// where  (s  - s ) =  /\ sigma = - del  , the first two are what
//          ||   |    /__\        2    zz
//          ||  ---
//
// is commonly called the CSA, the latter what is stored in GAMMA!
//


row_vector R1_CC(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		R1 : R1 values for all system spins relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons

  {
  int ns = sys.spins();				// Number of spins in the system
  row_vector R1s(ns);				// Vector of relaxation rates
  double R1=0;
  for(int i=0; i<ns; i++)
    {
    R1 = R1_CC(sys, i);				// Get the rate for spin i
    R1s.put(complex(R1,0),i);			// Store the relaxation rate
    }
  return R1s;
  }


double R1_CC(const sys_dynamic& sys, int i)

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		R1 : R1 value for spin i relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons

  {
  double W = sys.Omega(i);			// Larmor frequency (MHz)
  W *= pi2;					// Switch to angular frequency
  double wI = sys.lab_shift(i)*pi2;		// Shift of i in the lab frame (1/sec)
  double tau = sys.taux();			// Use x-axis tau value only (sec) 
  double delzz = sys.delz(i);			// Get anisotropy (PPM)
  double factor = W*W*delzz*delzz/10.0;		// Prefactor
  double term = 3.0*tau/(1+wI*wI*tau*tau); 	// Spectral density term
  double R1 = factor * term;
  return R1;
  }


double R1_CC_max(const sys_dynamic& sys, int i)
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            R1max : Maximum R1 value due to CSA relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    int ns = sys.spins();               // Number of spins in the system
    double R1, R1max=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(sys.isotope(i)==sys.isotope(j))// Only consider spins like spin i
        {
        R1 = R1_CC(sys, j);             // Get the rate for spin j
        if(R1 > R1max)                  // Take the largest R1 of these
          R1max = R1;
        }
    return R1max;
    }


double R1_CC_max(const sys_dynamic& sys, const std::string& Iso)
 
        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            R1max : Maximum R1 value due to CSA relaxation
        //                           over all spins in the system of isotope
        //                           type Iso
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    Isotope I(Iso);
    int ns = sys.spins();               // Number of spins in the system
    double R1, R1max=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(I == sys.isotope(j))		// Only consider spins like spin i
        {
        R1 = R1_CC(sys, j);             // Get the rate for spin j
        if(R1 > R1max)                  // Take the largest R1 of these
          R1max = R1;
        }
    return R1max;
    }


double R1_CC_max(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		R1 : Maximum of all R1 values of the system spins
	//			     relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  {
  int ns = sys.spins();				// Number of spins in the system
  double R1=0, R1max=0;
  for(int i=0; i<ns; i++)
    {
    R1 = R1_CC(sys, i);				// Get the rate for spin i
    if(R1max < R1)				// Set it to R1max if its big
      R1max = R1;
    }
  return R1max;
  }


  row_vector T1_CC(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output	       T1s : T1 values for all system spins relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons

  {
  int ns = sys.spins();
  row_vector T1s(ns);
  for(int i=0; i<ns; i++)
    T1s.put(complex(1.0/R1_CC(sys,i),0),i);
  return T1s;
  }


  double T1_CC(const sys_dynamic& sys, int i)

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		T1 : T1 value for spin i relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons

  { return (1.0/R1_CC(sys, i)); }


  double T1_CC_max(const sys_dynamic& sys, int i)
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            T1max : Maximum T1 value due to CSA relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    int ns = sys.spins();               // Number of spins in the system
    double T1, T1max=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(sys.isotope(i)==sys.isotope(j))// Only consider spins like spin i
        {
        T1 = T1_CC(sys, j);             // Get the rate for spin j
        if(T1 > T1max)                  // Take the largest T1 of these
          T1max = T1;
        }
    return T1max;
    }


double T1_CC_max(const sys_dynamic& sys, const std::string& Iso)
 
        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            T1max : Maximum T1 value due to CSA relaxation
        //                           over all spins in the system of isotope
        //                           type Iso
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    Isotope I(Iso);
    int ns = sys.spins();               // Number of spins in the system
    double T1, T1max=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(I == sys.isotope(j))		// Only consider spins like spin i
        {
        T1 = T1_CC(sys, j);             // Get the rate for spin j
        if(T1 > T1max)                  // Take the largest T1 of these
          T1max = T1;
        }
    return T1max;
    }


  double T1_CC_max(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		T1 : Maximum of all T1 values of the system spins
	//			     relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  {
  int ns = sys.spins();				// Number of spins in the system
  double T1=0, T1max=0;
  for(int i=0; i<ns; i++)
    {
    T1 = T1_CC(sys, i);				// Get the rate for spin i
    if(T1max < T1)				// Set it to T1max if its big
      T1max = T1;
    }
  return T1max;
  }


// ------------------- Transverse Relaxation, T2 ------------------------

// Ferrar & Becker, "Pulse and Fourier Transform NMR", Academic Press, New York, 1971 
//                               Eq. 4.29 page 59.
//
//
//  I    1            2    2         2  tau [       6              ]
// R  = --- = (gamma ) (H ) (s  - s )   ___ | _______________  + 8 |
//  2    T          I    0    ||   |        |         2     2      |
//        2                   ||  ---   90  | 1 + (w ) (tau)       |
//             				    [	    I              ]
//
//
//                    2     I  2   tau  [       3              ]
//          = (Omega )  (del  )   _____ | _______________  + 4 |
//                  I       zz          |       2       2      |
//                                 20   | 1 + (w ) (tau)       |
//             				[       I              ]
//             
//

  row_vector R2_CC(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		R2 : R2 values for system spins relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons

  {
  int ns = sys.spins();				// Number of spins in the system
  row_vector R2s(ns);				// Vector of relaxtion rates
  double R2=0;
  for(int i=0; i<ns; i++)
    {
    R2 = R2_CC(sys,i); 				// Calculate R2
    R2s.put(complex(R2,0),i);			// Store the relaxation rate
    }
  return R2s;
  }


  double R2_CC(const sys_dynamic& sys, int i)

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		R2 : R2 value for spin i relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons

  {
  double omega = sys.Omega(i);			// Larmor frequency (MHz)
  omega *= pi2;					// Switch to angular frequency
  double wi = sys.lab_shift(i);			// Shift of i in the lab frame (Hz)
  wi *= pi2;					// Switch to angular frequency
  double tau = sys.taux();			// Use x-axis tau value only (sec) 
  double delzz = sys.delz(i);			// Get anisotropy (PPM)

  double factor = omega*omega*delzz*delzz/20.0;	// Prefactor
  double term1 = 3.0/(1+wi*wi*tau*tau); 	// Spectral density term
  double term2 = 4.0;
  double R2 = factor * tau * (term1 + term2);
  return R2;
  }


  double R2_CC_max(const sys_dynamic& sys, int i)
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            R2max : Maximum R2 value due to CSA relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    int ns = sys.spins();               // Number of spins in the system
    double R2, R2max=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(sys.isotope(i)==sys.isotope(j))// Only consider spins like spin i
        {
        R2 = R2_CC(sys, j);             // Get the rate for spin j
        if(R2 > R2max)                  // Take the largest R2 of these
          R2max = R2;
        }
    return R2max;
    }


double R2_CC_max(const sys_dynamic& sys, const std::string& Iso)
 
        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            R2max : Maximum R2 value due to CSA relaxation
        //                           over all spins in the system of isotope
        //                           type Iso
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    Isotope I(Iso);
    int ns = sys.spins();               // Number of spins in the system
    double R2, R2max=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(I == sys.isotope(j))		// Only consider spins like spin i
        {
        R2 = R2_CC(sys, j);             // Get the rate for spin j
        if(R2 > R2max)                  // Take the largest R2 of these
          R2max = R2;
        }
    return R2max;
    }


double R2_CC_max(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		R2 : Maximum of all R2 values of the system spins
	//			     relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  {
  int ns = sys.spins();				// Number of spins in the system
  double R2=0, R2max=0;
  for(int i=0; i<ns; i++)
    {
    R2 = R2_CC(sys, i);				// Get the rate for spin i
    if(R2max < R2)				// Set it to R2max if its big
      R2max = R2;
    }
  return R2max;
  }


row_vector T2_CC(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		T2 : T2 values for system spins relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons

  {
  int ns = sys.spins();
  row_vector R2s = R2_CC(sys);
  row_vector T2s(ns);
  for(int i=0; i<ns; i++)
    T2s.put((1.0/R2s.get(i)),i);
  return T2s;
  }


  double T2_CC(const sys_dynamic& sys, int i)

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		T2 : T2 value for spin i relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the CSA interaction
	//			     constant as it is used for comparisons

  { return (1.0/R2_CC(sys, i)); }


  double T2_CC_max(const sys_dynamic& sys, int i)
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            T2max : Maximum T2 value due to CSA relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    int ns = sys.spins();               // Number of spins in the system
    double T2, T2max=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(sys.isotope(i)==sys.isotope(j))// Only consider spins like spin i
        {
        T2 = T2_CC(sys, j);             // Get the rate for spin j
        if(T2 > T2max)                  // Take the largest T2 of these
          T2max = T2;
        }
    return T2max;
    }


double T2_CC_max(const sys_dynamic& sys, const std::string& Iso)
 
        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            T2max : Maximum T2 value due to CSA relaxation
        //                           over all spins in the system of isotope
        //                           type Iso
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    Isotope I(Iso);
    int ns = sys.spins();               // Number of spins in the system
    double T2, T2max=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(I == sys.isotope(j))		// Only consider spins like spin i
        {
        T2 = T2_CC(sys, j);             // Get the rate for spin j
        if(T2 > T2max)                  // Take the largest T2 of these
          T2max = T2;
        }
    return T2max;
    }


  double T2_CC_max(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		T2 : Maximum of all T2 values of the system spins
	//			     relaxed by CSA
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  {
  int ns = sys.spins();				// Number of spins in the system
  double T2=0, T2max=0;
  for(int i=0; i<ns; i++)
    {
    T2 = T2_CC(sys, i);				// Get the rate for spin i
    if(T2max < T2)				// Set it to T2max if its big
      T2max = T2;
    }
  return T2max;
  }


// ------------------------- CSA Linewidths -----------------------------
//
//                                 CC
//                                R
//                           CC    2      1.0
//                         LW   = --- = -------
//                           hh    pi    CC 
//                                     	T  * pi
//                                       2
//

  row_vector LWhh_CC(const sys_dynamic& sys) 

	// Input 	       sys : Spin system
	// Output	       LWs : Expected linewidths at half-height
	//			     for system spins relaxed by CSA

// 		     LineWidth via CSA T2 Relaxation

//				LW = R /pi
//			              2	

  {
  row_vector LWs = R2_CC(sys); 
  LWs /= 3.14159;
  return LWs;
  }


  double LWhh_CC(const sys_dynamic& sys, int i) 

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		LW : Expected linewidth at half-height
	//			     for spin i relaxed by CSA

  {
  double R2 = R2_CC(sys, i); 
  return (R2/3.14159);
  }


  double LWhh_CC_max(const sys_dynamic& sys, int i)
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output          LWhhmax : Maximum LWhh value due to CSA relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    int ns = sys.spins();               // Number of spins in the system
    double LWhh, LWhhmax=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(sys.isotope(i)==sys.isotope(j))// Only consider spins like spin i
        {
        LWhh = LWhh_CC(sys, j);		// Get the linewidth for spin i
        if(LWhh > LWhhmax)		// Take the largest LWhh of these
          LWhhmax = LWhh;
        }
    return LWhhmax;
    }
 

double LWhh_CC_max(const sys_dynamic& sys, const std::string& Iso)
 
        // Input               sys : Spin system
        //                     Iso : A string for an Isotope type
        // Output          LWhhmax : Maximum LWhh value due to CSA relaxation
        //                           over all spins in the system of isotope
        //                           type as Iso
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    Isotope I(Iso);
    int ns = sys.spins();               // Number of spins in the system
    double LWhh, LWhhmax=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(I == sys.isotope(j))		// Only consider spins like spin i
        {
        LWhh = LWhh_CC(sys, j);		// Get the linewidth for spin i
        if(LWhh > LWhhmax)		// Take the largest LWhh of these
          LWhhmax = LWhh;
        }
    return LWhhmax;
    }
 

  double LWhh_CC_max(const sys_dynamic& sys) 

	// Input 	       sys : Spin system
	// Output		LW : Maximum expected linewidth at half-height
	//			     of all spins relaxed by CSA in the system

  {
  double LW=0, LWmax=0;
  for(int i=0; i<sys.spins(); i++)
    {
    LW = LWhh_CC(sys,i);
    if(LWmax < LW)
      LWmax = LW;
    }
  return LWmax;
  }


// ______________________________________________________________________
// *************** CSA-CSA Relaxation Auxiliary Functions ***************
// ______________________________________________________________________


 matrix xiCSA(const sys_dynamic& dsys)


	// Input		dsys  : A dynamic system
	// Return		xis   : A matrix of CSA interaction
	//				constants (rad/sec)
	// Note			      : The PPM units of delzz are
	//				counteracted by the MHz units
	//				on Omega
	//
	//		 1/2
	//   CSA   [6*pi]
	// xi    = |----| * gamma * B * del  (i) = K * Omega  * del  (i)
	//   i     [ 5  ]        i   0     zz               i      zz

  {
  matrix xis(dsys.spins(), dsys.spins(), complex0, d_matrix_type);
  double delzz, xii;
  double Omi;				// Spectrometer frequency (MHz)
  for(int i=0; i<dsys.spins(); i++)
    {
    Omi = dsys.Omega(i);		// Larmor frequency of spin i (MHz)
    delzz = dsys.delz(i);		// delzz of spin i (PPM)
    xii = K*pi2*Omi*delzz;		// 2*pi converts Hz to rad/sec
    xis.put(complex(xii,0.0), i, i);	// Store the Xi value for spin i
    }
  return xis;
  }


 double xiCSA(const sys_dynamic& dsys, int i)

	// Input		dsys  : A dynamic system
	//			i     : Spin index
	// Return		xii   : The gamma CSA interaction
	//				constant (rad/sec)
	// Note			      : The PPM units of delzz are
	//				counteracted by the MHz units
	//				on Omega
	//
	//		 1/2
	//   CSA   [6*pi]
	// xi    = |----| * gamma * B * del  (i) = K * Omega  * del  (i)
	//   i     [ 5  ]        i   0     zz               i      zz

  {
  double Omi = dsys.Omega(i);		// Frequency of spin i (MHz)
  double delzz = dsys.delz(i);		// delzz of spin i (PPM)
  double xii = K*pi2*Omi*delzz;		// 2*pi converts Hz to rad/sec
  return xii;
  }


 matrix xiCSA(const spin_system& sys, double* CSAs)


	// Input		dsys  : A dynamic system
	//			CSAs  : Vector of CSA values (PPM)
	// Return		xis   : A matrix of CSA interaction
	//				constants (rad/sec)
	// Note			      : The PPM units of CSA's are
	//				counteracted by the MHz units
	//				on Omega
	//
	//		 1/2
	//   CSA   [6*pi]               2  /\          2 
	// xi    = |----| * gamma * B * - /__\sigma  = _ * K * Omega  * CSA
	//   i     [ 5  ]        i   0  3          i   3            i      i

  {
  matrix xis(sys.spins(), sys.spins(), complex0, d_matrix_type);
  double csa=0, xii=0;
  double Omi;				// Spectrometer frequency (MHz)
  for(int i=0; i<sys.spins(); i++)
    {
    Omi = sys.Omega(i);			// Larmor frequency of spin i (MHz)
    csa = CSAs[i];			// CSA of spin i (PPM)
    xii = K*pi2*(2.0/3.0)*Omi*csa;	// 2*pi converts Hz to rad/sec
    xis.put(complex(xii,0.0), i, i);	// Store the Xi value for spin i
    }
  return xis;
  }


 double xiCSA(const spin_system& sys, int i, double csa)

	// Input		dsys  : A dynamic system
	//			i     : Spin index
	//			csa   : CSA value of spin i (PPM)
	// Return		xii   : The gamma CSA interaction
	//				constant (rad/sec)
	// Note			      : The PPM units of CSA is
	//				counteracted by the MHz units
	//				on Omega
	//
	//		 1/2
	//   CSA   [6*pi]               2  /\          2 
	// xi    = |----| * gamma * B * - /__\sigma  = _ * K * Omega  * CSA
	//   i     [ 5  ]        i   0  3          i   3            i      i

  {
  double Omi = sys.Omega(i);		// Frequency of spin i (MHz)
  double xii = K*pi2*(2.0/3.0)*Omi*csa;	// 2*pi converts Hz to rad/sec
  return xii;
  }


 row_vector CSA(const sys_dynamic& dsys)


	// Input		dsys  : A dynamic system
	// Return		csas  : A vector of CSA values in PPM
	//
	//   /\                2
	//  /__\sigma = CSA  = - del  (i)
	//           i     i   3    zz 

  {
  row_vector csas(dsys.spins());
  double csai;
  for(int i=0; i<dsys.spins(); i++)
    {
    csai = (2.0/3.0)*dsys.delz(i);	// Larmor frequency of spin i (MHz)
    csas.put(complex(csai,0), i);	// Store the Xi value for spin i
    }
  return csas;
  }


 double CSA(const sys_dynamic& dsys, int i)

	// Input		dsys  : A dynamic system
	// Return		csa   : CSA value of spin i in PPM
	//
	//   /\                2
	//  /__\sigma = CSA  = - del  (i)
	//           i     i   3    zz 

  { return (2.0/3.0)*dsys.delz(i); }


#endif /* __RELAX_CSA_CC__ */

