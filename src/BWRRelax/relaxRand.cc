/* relaxRand.cc ******************************************
**							**
**	               G A M M A			**
**						 	**
**	Random Field Relaxation      Implementation 	**
**						 	**
**	Copyright (c) 1993			 	**
**	Scott A. Smith				 	**
**							**
**	University of California, Santa Barbara		**
**	Department of Chemistry				**
**	Santa Barbara CA. 93106 USA			**
**						 	**
**      $Header: $
**						 	**
*********************************************************/

/*********************************************************
**						 	**
** 	Description				 	**
**						 	**
**	The following functions provide easy access	**
**	to Random Field relaxation related functions	**
**						 	**
*********************************************************/

#ifndef _relax_Random_cc_	// Is this file already included?
#define _relax_Random_cc_ 1	// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation          // This is the implementation
#endif


#include <BWRRelax/relaxRand.h>		// Include the header file
#include <BWRRelax/relaxNMR.h>		// Generic relaxation  routines
#include <Level1/nmr_tensor.h>		// Include common spin tensors
#include <LSLib/sys_dynamic.h>		// Include aniostropic system
#include <stdlib.h>


// ______________________________________________________________________
// *************** Random Field Relaxation Superoperators ***************
// ______________________________________________________________________


void RRRx(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                                     double tau, int type, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			tau   : System random field correlation times
	// 			type  : Relaxation type to compute
	//				   0 = RDM-RDM AC & CC
	//				   + = RDM-RDM AC
	//				   - = RDM-RDM CC
	//			level : Relaxation treatment level
	// Output		LOp   : Random Field relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho
	// Note			      :	Uses a spherical top spectral density

   {
//		    Prepare the Interaction Constants

   matrix Xis = xiRDM(sys);			// Matrix of Xi's for RDM
   if(!tau)					// Force extreme narrowing if no
     {						// tau value has been specified
     level=0;
     tau = 1.e-15;				// This is the default value used
     }						// in the xi functions also!

//			Prepare the Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *T;					// Spin tensors for each spin
   T = new spin_T[ns];				// Compiler didn't like spin_T T[ns]
   for(int i=0; i<ns; i++)			// Set them all to RDM spin tensors
     if(Re(Xis.get(i,i)))
       T[i] = T_RF(sys,i);

//		    Perform Computation Over All Spins

   Rij(LOp, sys, Ho, w, Xis, Xis,		// All Random Field parameters
                  T, T, tau, 1, type, level);
   return;
   }


super_op RRRx(const sys_dynamic& sys, gen_op& Ho, int type, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = RDM-RDM AC & CC
	//				   + = RDM-RDM AC
	//				   - = RDM-RDM CC
	//			level : Relaxation treatment level
	// Output		LOp   : Random Field relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho
	// Note			      :	Uses a spherical top spectral density

   {
//   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
   int ls = hs*hs;				// Total system Liouville space
   Ho.set_EBR();				// Insure Ho is in its eigenbasis
   matrix mx(ls, ls, 0.0);			// Construct zero superoperator
   super_op LOp(mx, Ho.get_basis());

   double tau = sys.tauR();			// Get random field tau value
   double *w;		   		 	// Set up for system energy levels (LAB)
   w = new double[hs];
   gen_op Holab;				// these will be in Hertz, ~10**8 in value
   if(abs(level) > 1)				// Used for higher level calculations
     {						// involving spectral density functions
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);				//	This is H in the lab frame
     Holab.Op_base(Ho, 1.e-7);			//	It must be in the Ho eigenbasis
     if(!Holab.test_EBR())
       std::cout << "\n\tWarning relax_Rand: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

   RRRx(LOp,sys,Ho,w,tau,type,level);		// Set LOp to RDM-RDM relaxation superop
   delete [] w;
   return LOp;
   }


void RRR(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double*w,
                     double* taus, double chi, int type, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = RDM-RDM AC & CC
	//				   + = RDM-RDM AC
	//				   - = RDM-RDM CC
	//			level : Relaxation treatment level
	// Output		LOp   : Random Field relaxation superoperator
	// Note			      :	Computed in the eigenbasis of Ho

   {
//		    Prepare the Interaction Constants

  matrix Xis = xiRDM(sys);			// Matrix of Xi's for RDM

//			Prepare the Spin Tensors

   int ns = sys.spins();			// Total number of spins
   spin_T *T;					// Spin tensors for each spin
   T = new spin_T[ns];				// Compiler didn't like spin_T T[ns]
   for(int i=0; i<ns; i++)			// Set them all to RDM spin tensors
     if(Re(Xis.get(i,i)))
       T[i] = T_RF(sys,i);

//			Prepare the Spatial "Tensors"

// Due to poorly defined spatial variations in the random field,
// only a single value is present, no asymmetry or Euler angles.
// Thus rather than a spatial tensor, only a double is needed.  This
// value is include in the chi's and here A is a dummy tensor only.

   space_T *A; 					// Spatial tensors for each spin
   A = new space_T[ns];				// Compiler didn't like space_T A[ns]
//  for(i=0; i<ns; i++)				// Set them all to RDM space tensors
//    if(Re(Xis.get(i,i)))
//      A[i] = sys.TC(i);

//		    Perform Computation Over All Spins

//   Rij(LOp, sys, Ho, w, Xis, Xis, A, A,	// All Random Field parameters
//                  T, T, taus, chi, type, level);
   Rij_rdm(LOp, sys, Ho, w, Xis, Xis, A, A,	// All Random Field parameters
                  T, T, taus, chi, type, level);
   return;
   }


super_op RRR(const sys_dynamic& sys, gen_op& Ho, int type, int level)

	// Input		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	// 			type  : Relaxation type to compute
	//				   0 = RDM-RDM AC & CC
	//				   + = RDM-RDM AC
	//				   - = RDM-RDM CC
	//			level : Relaxation treatment level
	// Output		LOp   : Random Field relaxation superoperator
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
int ns = sys.spins();
taus[0] = sys.TR(ns);
if(taus[0] == 0)
 taus[0] = 1.0e-15;
   double *w;		   		 	// Set up for system energy levels (LAB)
   w = new double[hs];
   gen_op Holab;				// these will be in Hertz, ~10**8 in value
   if(abs(level) > 1)				// Used for higher level calculations
     {
     Holab = Hcs_lab(sys);
     Holab += HJ(sys);				//	This is H in the lab frame
     Holab.Op_base(Ho, 1.0e-7);			//	It must be in the Ho eigenbasis
     if(!Holab.test_EBR())
       std::cout << "\n\tWarning relax_Rand: "
            << " Unable to Obtain Proper Ho(lab) Eigenbasis";
     Holab.eigvals(w);
     }

   RRR(LOp,sys,Ho,w,taus,chi,type,level);	// Set LOp to RDM-RDM relaxation superop
   delete [] w;
   return LOp;
   }



void Rij_rdm(super_op& LOp, const sys_dynamic& sys, gen_op& Ho, double* w,
         matrix& xi1s, matrix& xi2s, space_T* A1, space_T* A2,
         spin_T* T1, spin_T* T2, double* taus, double chi, int type, int level)

	// Input                LOp   : Superoperator
	// 			sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			w     : Transition frequencies
	//			xi1s  : Interaction constant, mu1
	//			xi2s  : Interaction constant, mu2
	//			A1    : Spatial tensor, mu1 (dummy)
	//			A2    : Spatial tensor, mu2 (dummy)
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			taus  : Sytem 5 effective correlation times
	//			chi   : Spin system chi value
	// 			type  : Relaxation type to compute
	//				   0 = Auto & Cross Correlation
	//				   + = Auto Correlation Only
	//				   - = Cross Correlation Only
	//			level : Relaxation treatment level

/*	              Two Single Spin Mechanisms

                      --- ---             --- ---
                      \   \               \   \
               LOp += /   /    R        = /   /   R
  	              --- ---   mu1,mu2   --- ---  i,j
                      mu1 mu2              i   j                       */

   {
   int rank=1;
   double xi1, xi1xi2;			// For specific interaction constants
//   double xi1, xi2, xi1xi2;			// For specific interaction constants
   double c1s[5];				// Set up 5 coefficients interaction 1
//   double c2s[5];				// Set up 5 coefficients interaction 2
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   T1s = new gen_op[2*rank+1];
   double w0=0,w1=0,w2=0;			// Used for J's levels 0 & 1
//   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
   int m;					// z-component of ang. momentum
   double cutoff = 1e-12;
   for(int i=0; i<ns; i++)			// Sum over spins i (mu1)
     {
     xi1 = Re(xi1s.get(i,i));			// Get spin i (mu1) interaction constant
     if(fabs(xi1) > cutoff)			// Only add non-trivial contributions
       {
       for(m=-rank; m<=rank; m++)		// Put spin tensor for i (mu1) into a
         {					// vector of operators in basis of Ho
         T1s[m+rank] = gen_op(((T1[i]).component(rank,m)));
//         T1s[m+rank] = gen_op(((T1[i]).component(rank,m)).matrix());
//         T1s[m+rank] = gen_op(T1[i].component(rank,m),
//                               sys.get_basis());
         T1s[m+rank].Op_base(Ho);
         }
       for(int j=0; j<ns; j++)			// Sum over spins j
         if(i==j)				// Auto-correlation term only
           {
           xi1xi2 = xi1*xi1;			// xi(mu1,i)*xi(mu1,i)
           if(fabs(xi1xi2) > cutoff)		// Only add non-trivial contributions
             Rmumu_rdm(LOp, T1s, T1s, w, hs, taus,
                c1s,c1s,xi1xi2,w0,w1,w2,0);
// forces correlated and extreme narrowing!
           }
       }
     }						// Increment first spin
   return;
   if(A1[0].exists()) type=0;			// Compiler likes these used
   if(A2[0].exists()) level=0;			// Compiler likes these used
   if(T2 == NULL) chi = xi2s.rows();		// Compiler likes these used
   }


   void Rmumu_rdm(super_op& LOp, gen_op* T1s, gen_op* T2s, double* w, int hs,
              double* taus, double* c1s, double* c2s, double xi1xi2,
               double w0, double w1, double w2, int level, int autoc)

	// Input                LOp   : Superoperator
	//                      T1s   : Spin tensor components of mu1
	// 			T2s   : Spin tensor components of mu2
	//			w     : Transisiton frequencies
	//			hs    : Hilbert space size
	//			taus  : 5 effective correlation times
	//			c1s   : J coefficients of mu1
	//			c2s   : J coefficients of mu2
	//			xi1xi2: Interaction constant product 
	//			w0    : Zero quantum transition frequency
	//			w1    : Single quantum transition frequency
	//			w2    : Double quantum transition frequency
	//			level : Relaxation treatment level
	//			autoc : Flag for auto correlation vs cross correlation
	// Output		LOp   : Relaxation superoperator
	// Note			      :	Computed in the basis T1s and T2s

//                            LOp += R   (Level)
//                                    1,2  

   {
   int rank = 1;
   matrix J12;					// Used for J's levels > 1
   double J0, J1, J2;				// Used for J's levels 0 & 1
   switch(level)
     {
     case 4:					// Level 4 mu1-mu2: element by element
       J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spectral densities
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       R_4(LOp, rank, T1s, T2s, J12);
       break;
     case -4:					// Level 4 mu1-mu2: double commutator
       J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spectral densities
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       R_4s(LOp, rank, T1s, T2s, J12);
       break;
     case 3:					// Level 3 mu1-mu2: element by element
       J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spectral densities
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       R_3(LOp, w, rank, T1s, T2s, J12);
       break;
     case -3:					// Level 3 mu1-mu2: double commutator
       J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spectral densities
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       R_3s(LOp, w, rank, T1s, T2s, J12);
       break;
     case 2:					// Level 2 mu1-mu2: element by element
       J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spectral densities
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       R_2(LOp, rank, T1s, T2s, J12);
       break;
     case -2:					// Level 2 mu1-mu2: double commutator
       J12 = J_reduced(w,hs,taus,c1s,c2s,1);	// Get all reduced spectral densities
       J12 *= complex(xi1xi2);			// Scale by the interaction constants
       R_2s(LOp, rank, T1s, T2s, J12);
       break;
     case 1: 					// Level 1 mu1-mu2: double commutator
       J0 = J_reduced(taus, c1s, c1s, w0, 1);	// J(0)
       J1 = J_reduced(taus, c1s, c1s, w1, 1);	// J(w0)
       J2 = J_reduced(taus, c1s, c1s, w2, 1);	// J(2w0)
       if(autoc)
         R_AC_1(T1s, LOp, rank,
                     J0*xi1xi2,J1*xi1xi2,J2*xi1xi2);
       else
         R_CC_1(T1s, T2s, LOp, rank,
                   J0*xi1xi2,J1*xi1xi2,J2*xi1xi2);
       break;
     case 0:					// Level 0 mu1-mu2: element by element
       J0 = J_gen(taus[0],0.0);			// Need J(0) only, extreme narrowing
       J0 *= xi1xi2;
       if(fabs(J0) > 1.e-6)
         R_0(LOp,rank,T1s,T2s,complex(J0));
       break;
     default:					// Level 0 mu1-mu2: double commutator
       J0 = J_reduced(taus,c1s,c2s,0.0,1);	// Need J(0) only, extreme narrowing
       if(fabs(xi1xi2*J0) > 1.e-6)
       { if(autoc)
         {
           R_AC_0(T1s, LOp, rank, xi1xi2*J0);
         }
         else
           R_CC_0(T1s,T2s,LOp,rank,xi1xi2*J0);
       }
       break;
     }
   return;
   }

 
// ______________________________________________________________________
// ***************** CLASSICAL RANDOM FIELD RELAXATION ******************
// ______________________________________________________________________


// ------------------ Longitudinal Relaxation, T1 -----------------------


  row_vector R1_RR(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		R1 : R1 values for all system spins relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons

  {
  int ns = sys.spins();				// Number of spins in the system
  row_vector R1s(ns);				// Vector of relaxation rates
  double R1=0;
  for(int i=0; i<ns; i++)
    {
    R1 = R1_RR(sys, i);				// Get the rate for spin i
    R1s.put(complex(R1,0),i);			// Store the relaxation rate
    }
  return R1s;
  }


  double R1_RR(const sys_dynamic& sys, int i)

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		R1 : R1 value for spin i relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons

// sosi: this sets all R1/T1
  { return R2_RR(sys,i); }


  double R1_RR_max(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		R1 : Maximum of all R1 values of the system spins
	//			     relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  {
  int ns = sys.spins();				// Number of spins in the system
  double R1=0, R1max=0;
  for(int i=0; i<ns; i++)
    {
    R1 = R1_RR(sys, i);				// Get the rate for spin i
    if(R1max < R1)				// Set it to R1max if its big
      R1max = R1;
    }
  return R1max;
  }


  row_vector T1_RR(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output	       T1s : T1 values for all system spins relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons

  {
  int ns = sys.spins();
  row_vector T1s(ns);
  for(int i=0; i<ns; i++)
    T1s.put(complex(1.0/R1_RR(sys,i),0),i);
  return T1s;
  }


  double T1_RR(const sys_dynamic& sys, int i)

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		T1 : T1 value for spin i relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons

  { return (1.0/R1_RR(sys, i)); }


  double T1_RR_max(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		T1 : Maximum of all T1 values of the system spins
	//			     relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  {
  int ns = sys.spins();				// Number of spins in the system
  double T1=0, T1max=0;
  for(int i=0; i<ns; i++)
    {
    T1 = T1_RR(sys, i);				// Get the rate for spin i
    if(T1max < T1)				// Set it to T1max if its big
      T1max = T1;
    }
  return T1max;
  }


// ------------------- Transverse Relaxation, T2 ------------------------


  row_vector R2_RR(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		R2 : R2 values for system spins relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons

  {
  int ns = sys.spins();				// Number of spins in the system
  row_vector R2s(ns);				// Vector of relaxtion rates
  double R2=0;
  for(int i=0; i<ns; i++)
    {
    R2 = R2_RR(sys,i); 				// Calculate R2
    R2s.put(complex(R2,0),i);			// Store the relaxation rate
    }
  return R2s;
  }


  double R2_RR(const sys_dynamic& sys, int i)

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		R2 : R2 value for spin i relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons

// sosi: this sets all R2/T2

//				LW = R /pi
//			              2	
  { return 3.14*LWhh_RR(sys,i); }


  double R2_RR_max(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		R2 : Maximum of all R2 values of the system spins
	//			     relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  {
  int ns = sys.spins();				// Number of spins in the system
  double R2=0, R2max=0;
  for(int i=0; i<ns; i++)
    {
    R2 = R2_RR(sys, i);				// Get the rate for spin i
    if(R2max < R2)				// Set it to R2max if its big
      R2max = R2;
    }
  return R2max;
  }


  row_vector T2_RR(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		T2 : T2 values for system spins relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons

  {
  int ns = sys.spins();
  row_vector R2s = R2_RR(sys);
  row_vector T2s(ns);
  for(int i=0; i<ns; i++)
    T2s.put((1.0/R2s.get(i)),i);
  return T2s;
  }


  double T2_RR(const sys_dynamic& sys, int i)

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		T2 : T2 value for spin i relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.
	// Note			   : This routine does not use the R.F. interaction
	//			     constant as it is used for comparisons

  { return (1.0/R2_RR(sys, i)); }


  double T2_RR_max(const sys_dynamic& sys)

	// Input 	       sys : Spin system
	// Output		T2 : Maximum of all T2 values of the system spins
	//			     relaxed by R.F.
	// Note			   : This routine assumes a symmetric shift tensor
	// Note			   : Here, the motion is assumed to be diffusive and
	//  			     the system moving as an isotropic top.

  {
  int ns = sys.spins();				// Number of spins in the system
  double T2=0, T2max=0;
  for(int i=0; i<ns; i++)
    {
    T2 = T2_RR(sys, i);				// Get the rate for spin i
    if(T2max < T2)				// Set it to T2max if its big
      T2max = T2;
    }
  return T2max;
  }


// ------------------------- R.F. Linewidths -----------------------------
//
//                                 RR
//                                R
//                           RR    2      1.0
//                         LW   = --- = -------
//                           hh    pi    RR 
//                                     	T  * pi
//                                       2
//

  row_vector LWhh_RR(const sys_dynamic& sys) 

	// Input 	       sys : Spin system
	// Output	       LWs : Expected linewidths at half-height
	//			     for system spins relaxed by R.F.

// 		     LineWidth via R.F. T2 Relaxation

//				LW = R /pi
//			              2	

  {
  row_vector LWs = R2_RR(sys); 
  LWs /= 3.14159;
  return LWs;
  }


  double LWhh_RR(const sys_dynamic& sys, int i) 

	// Input 	       sys : Spin system
	// 			 i : Spin index
	// Output		LW : Expected linewidth at half-height
	//			     for spin i relaxed by R.F.

// sosi: this sets all linewidths

  { return sys.TR(i); }


  double LWhh_RR_max(const sys_dynamic& sys, int i)
 
        // Input               sys : Spin system
        //                       i : A spin index       
        // Output          LWhhmax : Max. LWhh value from Rand. Field relaxation
        //                           over all spins in the system of isotope
        //                           type the same as the input spin i
        // Note                    : Here, the motion is assumed to be diffusive and
        //                           the system moving as an isotropic top.
 
    {
    int ns = sys.spins();               // Number of spins in the system
    double LWhh, LWhhmax=0;
    for(int j=0; j<ns; j++)             // Loop over all spins
      if(sys.isotope(i)==sys.isotope(j))// Only consider spins like spin i
        {
        LWhh = LWhh_RR(sys, j);		// Get the linewidth for spin i
        if(LWhh > LWhhmax)		// Take the largest LWhh of these
          LWhhmax = LWhh;
        }
    return LWhhmax;
    }
 

  double LWhh_RR_max(const sys_dynamic& sys, const std::string& Iso)
 
        // Input               sys : Spin system
        //                     Iso : A string for an Isotope type
        // Output          LWhhmax : Max. LWhh value due to Rand.Field relaxation
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
        LWhh = LWhh_RR(sys, j);		// Get the linewidth for spin i
        if(LWhh > LWhhmax)		// Take the largest LWhh of these
          LWhhmax = LWhh;
        }
    return LWhhmax;
    }
 

  double LWhh_RR_max(const sys_dynamic& sys) 

	// Input 	       sys : Spin system
	// Output		LW : Maximum expected linewidth at half-height
	//			     of all spins relaxed by Rand.Field in the system

  {
  double LW=0, LWmax=0;
  for(int i=0; i<sys.spins(); i++)
    {
    LW = LWhh_RR(sys,i);
    if(LWmax < LW)
      LWmax = LW;
    }
  return LWmax;
  }

// ______________________________________________________________________
// ************* Random Field Relaxation Auxiliary Functions ************
// ______________________________________________________________________


 matrix xiRDM(const sys_dynamic& dsys)


	// Input		dsys  : A dynamic system
	// Return		xis   : A matrix of RDM interaction
	//				constants (rad/sec)
	// Note			      : The random field xi values are
	//				based on SQT linewidths

  {
  matrix xis(dsys.spins(), dsys.spins(), complex0, d_matrix_type);
  complex z;
  for(int i=0; i<dsys.spins(); i++)
    {
    z = complex(dsys.xiR(i));
    xis.put(z, i, i);		// Store the Xi value for spin i
    }
  return xis;
  }


 double xiRDM(const sys_dynamic& dsys, int i)

	// Input		dsys  : A dynamic system
	// Return		xis   : A matrix of RDM interaction
	//				constants (rad/sec)
	//
  { return dsys.xiR(i); }


#endif /* __RELAX_Random_CC__ */

