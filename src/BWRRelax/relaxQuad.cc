/* relax_Quad.cc *****************************************
**							**
**	              G A M M A				**
**						 	**
**	NMR Library, Quadrupolar Relaxation Functions	**
**						 	**
**	Implementation   				**
**						 	**
**	Copyright (c) 1993			 	**
**	Scott A. Smith				 	**
**							**
**	University of California, Santa Barbara		**
**	Department of Chemistry				**
**	Santa Barbara CA. 93106 USA			**
**						 	**
**      $Header:
**						 	**
** Modifications:					**
**						 	**
*********************************************************/

/*********************************************************
**						 	**
** 	Description				 	**
**						 	**
**	The following functions provide easy access	**
**	to several quadrupolar relaxation functions	**
**						 	**
*********************************************************/

#ifndef RelaxQuad_			// Is this file already included?
#define RelaxQuad_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// This is the implementation
#endif

#include <Basics/Gconstants.h>		// Know Pix2
#include <BWRRelax/relaxQuad.h>		// Include the file header
#include <BWRRelax/relaxNMR.h>		// Include generic relaxation
#include <Level1/nmr_tensor.h>		// Include common spin tensors
#include <LSLib/sys_dynamic.h>		// Include aniostropic systems
#include <stdlib.h>

const double K = 1.941625913;			// sqrt[6*pi/5]

// ______________________________________________________________________
// ********* Quadrupole-Quadrupole Relaxation Superoperators ************
// ______________________________________________________________________

// ----------------- Quad-Quad, No Applied RF-Field ------------------

void RQQ(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
                     double* taus, double chi, int type, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Q-Q AC & CC
	//				   + = Q-Q AC
	//				   - = Q-Q CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quadrupolar relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

   matrix Xi = xiQ(sys);			// Quadrupolar interaction constants

//			Prepare the Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *T;					// Spin tensors for each spin
   T = new spin_T[ns];				// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)				// Set them to quadrupolar spin tensors
     if(Re(Xi.get(i,i)))
       T[i] = T_Q(sys,i);

//			Prepare the Spatial Tensors

   space_T *A; 					// Spatial tensors for each spin
   A = new space_T[ns];				// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them to quadrupolar space tensors
     if(Re(Xi.get(i,i)))
       A[i] = sys.TQ(i);

//		    Perform Computation Over All Spin Pairs

   Rij(LOp, sys, Ho, w, Xi, Xi, A, A,		// All quadrupolar parameters
               T, T, taus, chi, type, level);
   return;
   }


super_op RQQ(const sys_dynamic& sys, gen_op& Ho, int type, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = Q-Q AC & CC
	//				   + = Q-Q AC
	//				   - = Q-Q CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quadrupolar relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//   int ns = sys.spins();			// Total number of spins
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
     {						// where spectral densities are evaluated
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);				//	This is H in the lab frame
     Holab.Op_base(Ho, 1.e-7);			//	It must be the Ho eigenbasis
     if(!Holab.test_EBR())
       std::cout << "\n\tWarning relax_Quad: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

   RQQ(LOp, sys, Ho, w, taus, chi, type, level);
   delete [] w;
   return LOp;
   }

// --------------- Quad-Quad, With An Applied RF-Field ------------------

void RQQrf(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int type, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			Wrflab: Field freqency in the lab frame (Hz)
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Quad-Quad AC & CC
	//				   + = Quad-Quad AC
	//				   - = Quad-Quad CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quad relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

   matrix Xis = xiQ(sys);			// Quadrupolar interaction constants

//			Prepare the Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *T;					// Spin tensors for each spin
   T = new spin_T[ns];				// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)				// Set them all to Quad spin tensors
     if(Re(Xis.get(i,i)))
       T[i] = T_Q(sys,i);

//			Prepare the Spatial Tensors

   space_T *A; 					// Spatial tensors for each spin
   A = new space_T[ns];				// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them all to Quad space tensors
     if(Re(Xis.get(i,i)))
       A[i] = sys.TQ(i);

//		    Perform Computation Over All Spins

   Rrfij(LOp, sys, Heff, w, Wrflab, Xis, Xis, A, A,	// All Quad parameters
                         T, T, taus, chi, type, level);
   return;
   }


super_op RQQrf(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type, int level)

	// Input                sys   : Dynamic spin system
	//                      Heff  : Effective Hamiltonian (Hertz)
	//                      Wrf   : RF-Field frequecy, lab frame (Hertz)
	//                      phi   : Field phase angle (Degrees)
	// 			type  : Relaxation type to compute
	//				   0 = Quad-Quad AC & CC
	//				   + = Quad-Quad AC
	//				   - = Quad-Quad CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quad relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff

   {
//   int ns = sys.spins();			// Total number of spins
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

   RQQrf(LOp,sys,Heff,w,Wrf,taus,chi,type,level);	// Set LOp to Quad-Quad relaxation superop
   delete [] w;
   return LOp;
   }

// ---- Quad-Quad, No Applied RF-Field, Dynamic Frequency Shift  -----


void RQQds(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int type, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Q-Q AC & CC
	//				   + = Q-Q AC
	//				   - = Q-Q CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quadrupolar relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

   matrix Xi = xiQ(sys);			// Quadrupolar interaction constants

//			Prepare the Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *T;					// Spin tensors for each spin
   T = new spin_T[ns];				// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)				// Set them to quadrupolar spin tensors
     if(Re(Xi.get(i,i)))
       T[i] = T_Q(sys,i);

//			Prepare the Spatial Tensors

   space_T *A; 					// Spatial tensors for each spin
   A = new space_T[ns];				// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them to quadrupolar space tensors
     if(Re(Xi.get(i,i)))
       A[i] = sys.TQ(i);

//		    Perform Computation Over All Spin Pairs

   Rijds(LOp, sys, Ho, w, Xi, Xi, A, A,		// All Quadrupolar parameters
               T, T, taus, chi, type, level);
   return;
   }


super_op RQQds(const sys_dynamic& sys, gen_op& Ho, int type, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = Q-Q AC & CC
	//				   + = Q-Q AC
	//				   - = Q-Q CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quadrupolar relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//   int ns = sys.spins();			// Total number of spins
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
     {						// where spectral densities are evaluated
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);				//	This is H in the lab frame
     Holab.Op_base(Ho, 1.e-7);			//	It must be the Ho eigenbasis
     if(!Holab.test_EBR())
       std::cout << "\n\tWarning relax_Quad: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

   RQQds(LOp, sys, Ho, w, taus, chi, type, level);
   delete [] w;
   return LOp;
   }


// -- Quad-Quad, With An Applied RF-Field, Dynamic Frequency Shift  --


void RQQrfds(super_op& LOp, const sys_dynamic& sys, gen_op& Heff, double*w,
                 double Wrflab, double* taus, double chi, int type, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Heff  : Effective Hamiltonian
	//			Wrflab: Field freqency in the lab frame (Hz)
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Quad-Quad AC & CC
	//				   + = Quad-Quad AC
	//				   - = Quad-Quad CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quad relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

   matrix Xis = xiQ(sys);			// Quadrupolar interaction constants

//			Prepare the Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *T;					// Spin tensors for each spin
   T = new spin_T[ns];				// Compiler didn't like spin_T T[ns]
   int i;
   for(i=0; i<ns; i++)			// Set them all to Quad spin tensors
     if(Re(Xis.get(i,i)))
       T[i] = T_Q(sys,i);

//			Prepare the Spatial Tensors

   space_T *A; 					// Spatial tensors for each spin
   A = new space_T[ns];				// Compiler didn't like space_T A[ns]
   for(i=0; i<ns; i++)				// Set them all to Quad space tensors
     if(Re(Xis.get(i,i)))
       A[i] = sys.TQ(i);

//		    Perform Computation Over All Spins

   Rrfijds(LOp, sys, Heff, w, Wrflab, Xis, Xis, A, A,	// All Quad parameters
                         T, T, taus, chi, type, level);
   return;
   }


super_op RQQrfds(const sys_dynamic& sys, gen_op& Heff, double Wrf, int type, int level)

	// Input                sys   : Dynamic spin system
	//                      Heff  : Effective Hamiltonian (Hertz)
	//                      Wrf   : RF-Field frequecy, lab frame (Hertz)
	//                      phi   : Field phase angle (Degrees)
	// 			type  : Relaxation type to compute
	//				   0 = Quad-Quad AC & CC
	//				   + = Quad-Quad AC
	//				   - = Quad-Quad CC
	//			level : Relaxation treatment level
	// Output		LOp   : Quad relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Heff

   {
//   int ns = sys.spins();			// Total number of spins
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

   RQQrfds(LOp,sys,Heff,w,Wrf,taus,chi,type,level);	// Set LOp to Quad-Quad relaxation superop
   delete [] w;
   return LOp;
   }


// ______________________________________________________________________
// **************** CLASSICAL QUADRUPOLAR RELAXATION ********************
// ______________________________________________________________________

// -------------------- Longitudinal Relaxation -------------------------

// Sudmeier, Anderson, and Frye, "Calculation of Nuclear Spin Relaxation Times",
//        Conc. Magn. Reson., 1990, 2, 197-212: page 201 equation [25]
//
//                                     2 
//       1      3(2I+3)      2 [    eta  ][   2*tau         8*tau     ]
// R  = --- = ----------- QCC  |1 + ---  || ---------- + ------------ |
//  1    T        2            [     3   ]|      2   2         2   2  |
//        1   400I (2I-1)                 [ 1 + w tau    1 + 4w tau   ]
//
// where QCC is the quadrupolar coupling constant in angular frequency units


row_vector R1_QQ(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		R1 : R1 values for all spins under quad. relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  {
  int ns = sys.spins();				// Number of spins in the system
  row_vector R1s(ns);				// Vector of relaxation rates
  double R1=0;
  for(int i=0; i<ns; i++)
    {
    R1 = R1_QQ(sys, i);				// Get the rate for spin i
    R1s.put(complex(R1,0),i);			// Store the relaxation rate
    }
  return R1s;
  }


double R1_QQ(const sys_dynamic& sys, int i) 

	// Input	       sys : Spin system
	// 	 	         i : Spin index
	// Output		R1 : R1 value for spin due to quadrupolar relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


  {
  double I = sys.qn(i);				// Get the spin quantum number
  double Num = 3.0*(2.0*I+3.0);			// Calculate the prefactor numerator
  double Den = 400.0*I*I*(2*I-1.0);		// Calculate the prefactor denominator
  double prefac = Num/Den;
  double tau = sys.taux();			// Use x-axis tau value only (sec) 
  double w = sys.lab_shift(i);			// Get the shift in MHz
  w *= PIx2;					// Switch to angular frequency
  double Q = sys.QCC(i);			// Get the quadrupolar coupling constant
  double eta = sys.Qeta(i);			// Get the quadrupolar asymmetry
  double etat = 1.0 + eta*eta/3.0;		// Set the eta term factor
  double term1 = 2.0*tau/(1+w*w*tau*tau); 	// Spectral density term
  double term2 = 8.0*tau/(1+4.0*w*w*tau*tau); 	// Spectral density term
  Q *= PIx2;
  double R1 = prefac*Q*Q*etat*(term1+term2);
  return R1;
  }


double R1_QQ_max(const sys_dynamic& sys, int i)
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            R1max : Maximum R1 value due to quad relaxation
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
        R1 = R1_QQ(sys, j);             // Get the rate for spin j
        if(R1 > R1max)                  // Take the largest R1 of these
          R1max = R1;
        }
    return R1max;
    }
 
 
double R1_QQ_max(const sys_dynamic& sys, const std::string& Iso)
 
        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            R1max : Maximum R1 value due to quad relaxation
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
      if(I == sys.isotope(j))           // Only consider spins like spin i
        {
        R1 = R1_QQ(sys, j);             // Get the rate for spin j
        if(R1 > R1max)                  // Take the largest R1 of these
          R1max = R1;
        }
    return R1max;
    }


double R1_QQ_max(const sys_dynamic& sys)
 
        // Input               sys : Spin system
        // Output               R1 : Maximum of all R1 values of the system spins
        //                           relaxed by quadrupolar mechanisms
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
  {
  int ns = sys.spins();                         // Number of spins in the system
  double R1=0, R1max=0;
  for(int i=0; i<ns; i++)
    {
    R1 = R1_QQ(sys, i);                         // Get the rate for spin i
    if(R1max < R1)                              // Set it to R1max if its big
      R1max = R1;
    }
  return R1max;
  }


row_vector T1_QQ(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		T1 : T1 values for all spins under quad. relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  {
  int ns = sys.spins();
  row_vector T1s(ns);
  for(int i=0; i<ns; i++)
    T1s.put(complex(1.0/R1_QQ(sys,i),0),i);
  return T1s;
  }


double T1_QQ(const sys_dynamic& sys, int i)

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		T1 : T1 value for spin i under quad. relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  { return (1.0/R1_QQ(sys, i)); }
 
 
double T1_QQ_max(const sys_dynamic& sys, int i)
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            T1max : Maximum T1 value due to quad relaxation
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
        T1 = T1_QQ(sys, j);             // Get the rate for spin j
        if(T1 > T1max)                  // Take the largest T1 of these
          T1max = T1;
        }
    return T1max;
    }
 
 
double T1_QQ_max(const sys_dynamic& sys, const std::string& Iso)
 
        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            T1max : Maximum T1 value due to quad relaxation
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
      if(I == sys.isotope(j))           // Only consider spins like spin i
        {
        T1 = T1_QQ(sys, j);             // Get the rate for spin j
        if(T1 > T1max)                  // Take the largest T1 of these
          T1max = T1;
        }
    return T1max;
    }

 
double T1_QQ_max(const sys_dynamic& sys)
 
        // Input               sys : Spin system
        // Output               T1 : Maximum of all T1 values of the system spins
        //                           relaxed by quadrupolar mechanisms
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
  {
  int ns = sys.spins();                         // Number of spins in the system
  double T1=0, T1max=0;
  for(int i=0; i<ns; i++)
    {
    T1 = T1_QQ(sys, i);                         // Get the rate for spin i
    if(T1max < T1)                              // Set it to T1max if its big
      T1max = T1;
    }
  return T1max;
  }
 

// --------------------- Transverse Relaxation --------------------------

// Sudmeier, Anderson, and Frye, "Calculation of Nuclear Spin Relaxation Times",
//        Conc. Magn. Reson., 1990, 2, 197-212: page 201 equation [26]
//
//                                     2 
//       1    3(2I+3)*tau    2 [    eta  ][          5            2       ]
// R  = --- = ----------- QCC  |1 + ---  || 3 + ---------- + ------------ |
//  2    T        2            [     3   ]|          2   2         2   2  |
//        2   400I (2I-1)                 [     1 + w tau    1 + 4w tau   ]
//
// where QCC is the quadrupolar coupling constant in angular frequency units


row_vector R2_QQ(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		R2 : R2 values for spins via quadrupolar relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  {
  int ns = sys.spins();				// Number of spins in the system
  row_vector R2s(ns);				// Vector of relaxtion rates
  double R2=0;
  for(int i=0; i<ns; i++)
    {
    R2 = R2_QQ(sys,i); 			// Calculate R2
    R2s.put(complex(R2,0),i);			// Store the relaxation rate
    }
  return R2s;
  }


double R2_QQ(const sys_dynamic& sys, int i) 

	// Input	       sys : Spin system
	// 	 	         i : Spin index
	// Output		R2 : R2 value for spin due to quadrupolar relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.


  {
  double I = sys.qn(i);				// Get the spin quantum number
  double Num = 3.0*(2.0*I+3.0);			// Calculate the prefactor numerator
  double Den = 400.0*I*I*(2.0*I-1.0);		// Calculate the prefactor denominator
  double prefac = Num/Den;
  double w = sys.lab_shift(i);			// Get the shift in MHz
  w *= PIx2;					// Switch to angular frequency
  double tau = sys.taux();			// Use x-axis tau value only (sec) 
  double Q = sys.QCC(i);			// Get the quadrupolar coupling constant
  double eta = sys.Qeta(i);			// Get the quadrupolar asymmetry
  Q *= PIx2;					// Switch to angular frequency
  double etat = 1.0 + eta*eta/3.0;		// Set the eta term factor
  double term1 = 3.0*tau;			// Spectral density term
  double term2 = 5.0*tau/(1+w*w*tau*tau); 	// Spectral density term
  double term3 = 2.0*tau/(1+4.0*w*w*tau*tau); 	// Spectral density term
  double R1 = prefac*Q*Q*etat*(term1+term2+term3);
  return R1;
  }
 
 
double R2_QQ_max(const sys_dynamic& sys, int i)
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            R2max : Maximum R2 value due to quad relaxation
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
        R2 = R2_QQ(sys, j);             // Get the rate for spin j
        if(R2 > R2max)                  // Take the largest R2 of these
          R2max = R2;
        }
    return R2max;
    }


double R2_QQ_max(const sys_dynamic& sys, const std::string& Iso)
 
        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            R2max : Maximum R2 value due to quad relaxation
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
      if(I == sys.isotope(j))           // Only consider spins like spin i
        {
        R2 = R2_QQ(sys, j);             // Get the rate for spin j
        if(R2 > R2max)                  // Take the largest R2 of these
          R2max = R2;
        }
    return R2max;
    }
 
 
double R2_QQ_max(const sys_dynamic& sys)
 
        // Input               sys : Spin system
        // Output               R2 : Maximum of all R2 values of the system spins
        //                           relaxed by quadrupolar mechanisms
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
  {
  int ns = sys.spins();                         // Number of spins in the system
  double R2=0, R2max=0;
  for(int i=0; i<ns; i++)
    {
    R2 = R2_QQ(sys, i);                         // Get the rate for spin i
    if(R2max < R2)                              // Set it to R2max if its big
      R2max = R2;
    }
  return R2max;
  }
 

row_vector T2_QQ(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		T2 : T2 values for spins via quadrupolar relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  {
  int ns = sys.spins();
  row_vector R2s = R2_QQ(sys);
  row_vector T2s(ns);
  for(int i=0; i<ns; i++)
    T2s.put((1.0/R2s.get(i)),i);
  return T2s;
  }


double T2_QQ(const sys_dynamic& sys, int i)

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		T2 : T2 value for spini via quadrupolar relaxation
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  { return (1.0/R2_QQ(sys, i)); }


double T2_QQ_max(const sys_dynamic& sys, int i)
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output            T2max : Maximum T2 value due to quad relaxation
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
        T2 = T2_QQ(sys, j);             // Get the rate for spin j
        if(T2 > T2max)                  // Take the largest T2 of these
          T2max = T2;
        }
    return T2max;
    }


double T2_QQ_max(const sys_dynamic& sys, const std::string& Iso)
 
        // Input               sys : Spin system
        //                     Iso : A string for an isotope
        // Output            T2max : Maximum T2 value due to quad relaxation
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
      if(I == sys.isotope(j))           // Only consider spins like spin i
        {
        T2 = T2_QQ(sys, j);             // Get the rate for spin j
        if(T2 > T2max)                  // Take the largest T2 of these
          T2max = T2;
        }
    return T2max;
    }

 
double T2_QQ_max(const sys_dynamic& sys)
 
        // Input               sys : Spin system
        // Output               T2 : Maximum of all T2 values of the system spins
        //                           relaxed by quadrupolar mechanisms
        // Note                    : This routine assumes a symmetric shift tensor
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
  {
  int ns = sys.spins();                         // Number of spins in the system
  double T2=0, T2max=0;
  for(int i=0; i<ns; i++)
    {
    T2 = T2_QQ(sys, i);                         // Get the rate for spin i
    if(T2max < T2)                              // Set it to T2max if its big
      T2max = T2;
    }
  return T2max;
  }
 

// -------------------- Quadrupolar Linewidths -------------------------

//                                 QQ
//                                R
//                           QQ    2      1.0
//                         LW   = --- = -------
//                           hh    pi    QQ 
//                                     	T  * pi
//                                       2
//

row_vector LWhh_QQ(const sys_dynamic& sys) 

	// Input 	       sys : Spin system
	// Output	       LWs : Expected linewidths at half-height
	//			     for system spins under quad relaxtion

  {
  row_vector LWs = R2_QQ(sys); 
  LWs /= 3.14159;
  return LWs;
  }


double LWhh_QQ(const sys_dynamic& sys, int i) 

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		LW : Expected linewidth at half-height
	//			     for spin i under quad relaxtion

  {
  double R2 = R2_QQ(sys, i); 
  return (R2/3.14159);
  }

 
double LWhh_QQ_max(const sys_dynamic& sys, int i)
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output          LWhhmax : Maximum LWhh value due to quad relaxation
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
        LWhh = LWhh_QQ(sys, j);         // Get the linewidth for spin i
        if(LWhh > LWhhmax)              // Take the largest LWhh of these
          LWhhmax = LWhh;
        }
    return LWhhmax;
    }


double LWhh_QQ_max(const sys_dynamic& sys, const std::string& Iso)
 
        // Input               sys : Spin system
        //                     Iso : A string for an Isotope type
        // Output          LWhhmax : Maximum LWhh value due to quad relaxation
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
      if(I == sys.isotope(j))           // Only consider spins like spin i
        {
        LWhh = LWhh_QQ(sys, j);         // Get the linewidth for spin i
        if(LWhh > LWhhmax)              // Take the largest LWhh of these
          LWhhmax = LWhh;
        }
    return LWhhmax;
    }
 

double LWhh_QQ_max(const sys_dynamic& sys) 

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		LW : Expected linewidth at half-height
	//			     for spin i under quad relaxtion

  {
  int ns = sys.spins();
  double R2max = 0, R2;
  for(int i=0; i<ns; i++)
    {
    R2 = R2_QQ(sys,i);
    if(R2 > R2max)
      R2max = R2;
    }
  return (R2max/3.14159);
  }


// ______________________________________________________________________
// ********** QUADRUPOLE-QUADRUPOLE RELAXATION AUXILIARY FUNCTIONS ******
// ______________________________________________________________________


matrix xiQ(const sys_dynamic& sys)


	// Input		dsys  : A dynamic system (this)
	// Return		dximx : A matrix of quadrupolar interaction
	//				constants (rad/sec)
	//
	//		 1/2   QCC   	        1/2  del  (i)
	//   Q     [6*pi]         i       [6*pi]        zz
	// xi    = |----| * ----------  = |----|  * ----------
	//   i     [ 5  ]   2I (2I -1)    [ 5  ]    2I (2I -1)
	//                    i   i                   i   i

  {
  matrix xivec(sys.spins(), sys.spins(), complex(0.0,0.0), d_matrix_type);
  double I=0.5, delzz, xii, Ifact;
  for(int i=0; i<sys.spins(); i++)
    {
    xii = 0;					// Set the Xi value to zero
    I = sys.qn(i);				// Get the spin quantum number
    if(I > 0.5)					// Only mess with I>1/2 spins
      {
      delzz = sys.QCC(i);			// delzz of spin i (Hz)
      Ifact = 2.0*I*(2.0*I - 1.0);		// I based denomenator
      xii = K*PIx2*delzz/Ifact;			// 2*pi converts Hz to rad/sec
      }
    xivec.put(xii, i, i);			// Store the Xi value for spin i
    }
  return xivec;
  }


double xiQ(const sys_dynamic& sys, int i)

	// Input		dsys  : A dynamic system (this)
	//			i     : Spin index
	// Return		xii   : Quadrupolar interaction constant (rad/sec)
	//
	//		 1/2   QCC   	        1/2  del  (i)
	//   Q     [6*pi]         i       [6*pi]        zz
	// xi    = |----| * ----------  = |----|  * ----------
	//   i     [ 5  ]   2I (2I -1)    [ 5  ]    2I (2I -1)
	//                    i   i                   i   i

  {
  double xii = 0;
  double I = sys.qn(i);
  double delzz=0;
  double Ifact=0;
  if(I > 0.5)
    {
    delzz = sys.QCC(i);				// delzz of spin i (Hz)
    Ifact = 2.0*I*(2.0*I - 1.0);		// I based denomenator
    xii = K*PIx2*delzz/Ifact;			// 2*pi converts Hz to rad/sec
    }
  return xii;
  }


#endif						// relaxQuad.cc


