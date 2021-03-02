/* FloqAcq.cc ***************************************************-*-c++-*-
**									**
**	Acquisitons in Floquet Space		Implementation		**
**									**
**	Copyright (c) 1992					 	**
**	Tilo Levante             		 			**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fuer physikalische Chemie				 	**
**	8092 Zurich / Switzerland				 	**
**								 	**
**      $Header: 
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** Description							 	**
**								 	**
** This module supports use of Floquet operators in GAMMA simulations.	**
** In particular, these functions allow for facile calculations of FIDs	**
** (free induction decays). It must be used in conjunction with GAMMA	**
** Floquet operator classes.						**
**								 	**
*************************************************************************/

#ifndef   Floq_acq_cc_				// Is file already included?
#  define Floq_acq_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation			// The this is implementation
#  endif

#include <Floquet/FloqAcq.h>			// Include the header file
#include <Level2/TrnsTable1D.h>			// Include transition tables
#include <vector>				// Include libstdc++ vectors
#include <Basics/Gconstants.h>			// Include PI, etc.

// ____________________________________________________________________________
//		          Generic Free Induction Decays
// ____________________________________________________________________________


void FID(floq_op sig0,gen_op& D,floq_op& H,double dt,int N,row_vector& fid)


        // Input                sig0  : Operator propagated (initial dens. mx.)
        //                      D     : Detection operator in trace computation
        //                      H     : Hamiltonian for propagation (in Hertz)
        //                      dt    : Evolution time (seconds, per point)
        //                      N     : Number of points to generate
        //                      fid   : Data vector containing the result
        // Output                     : None, FID data vector filled
        // Note                       : Assumed fid is at least as large as N
        // Note                       : Propagator is unitless, so 2*PI factor
        //                              used to get radians in exponent

/* At acquisition point k, which corresponds to time t = k*td, we have
   when working in the eigenbasis of the propagator U = exp(-2*pi*i*H*td)

         ---               itm/Om
fid(k) = \    Trace { D * e       * sigma(t) }
         /             m
         ---
          m

         ---               i*k*td*m/Om    k                       k
       = \    Trace { D * e            * U (td) * sig0 * adjoint[U (td)] }
         /             m
         ---
          m

         ---                i*k*td*m/Om       k                           k
       = \    { <i|D |j> * e            * <j|U (td)|j><j|sig0|i>*conj(<i|U (td)|i>) }
         /          m
         ---
        i,j,m

         ---                                            i*td*m/Om k
       = \    { <i|D |j><j|sig0|i>*[<j|U|j>*conj(i|U|i>e         ]  }
         /          m
         ---
        i,j,m

         ---                   k      ---               k 
       = \    { A(i,j) * B(i,j)  }  = \    { A(p) * B(p)  }
         /                            /          
         ---                          ---
         i,j                           p                                     */

  {
  H.set_EBR();				// Put H into an eigenbase
  sig0.Op_base(H);			// Put sigma in H eigenbase
  complex z(0,-2.0*PI*dt);              // Exponential factor for H -> U
//  floq_op U = (z*H).exp();		// Set td evolution propagator
// sosi	- use of member function above would be better	
floq_op U = exp(z*H); // evolution propagator
  double omega = H.omega();		// Fourier expansion frequency
  int nf = H.phodim();			// Phonon space dimension
  int fs = H.size();			// Floquet space dimension
  int ls = (2*nf+1)*fs*fs;		// Full Liouville dimension

  complex *A = new complex[ls];		// Array for A(p)
  complex *B = new complex[ls];		// Array for B(p)

  int i,j,m,pos;
  double Pi2Omdt = dt*omega*2.0*PI;
  complex ph;
  for(pos=0, m=-nf; m<=nf; m++)		// Generate A & B arrays
    {
    floq_op FD(nf,H.hsdim(),omega);	//   For floquet detection op.
    FD.put_sdiag(D, m);			//   Fill detection operator
    FD.Op_base(H);			//   Put FD in H eigenbas
    ph = exp(complex(0, m*Pi2Omdt));	//   Exponential factor
    for(i=0; i<fs; i++)			//   Loop Floquet transitions
      {					//   (single Liouville space)
      for(j=0; j<fs; j++)		//   Only count non-zero
	{ 				//   contributors to the acquisiton
	A[pos]=FD.get(i,j)*sig0.get(j,i);
	B[pos]=conj(U.get(i,i),U.get(j,j))*ph;
	if((square_norm(A[pos])>1.0e-12)&& 
	   (square_norm(B[pos])>1.0e-12)) 
	  pos ++;
	}
      }
    }

  for(int k=0; k<N; k++)		// Loop desired FID points (times)
    {
    complex z(0);			// Begin with zero point  
    for(int p=0; p<pos; p++)		// Loop all point contributors
      {
      z += A[p];			//   Add contribution from p
      A[p] *= B[p];			//                      k+1
      } 				//   Adjust A so its A*B   (next pt)
    fid.put(z,k);				// Set point, move to next k (time)
    } 
  delete [] A;				// Delete complex array A
  delete [] B;				// Delete complex array B
  }


void FID(floq_op &sigma , gen_op& det, floq_op& ham, row_vector& fid)

	// Input	sigma : operator propagated (initial density mx)
	// 		det   : "detction operator" in trace computation
	//	 	ham   : "hamiltonian" for propagation (in Hertz)
	//		fid   :	data row_vector containing the result
	// Output	      : none, fid data row_vector filled

  {
int np = 0;				// Number of points in the fid
int dw = 0;				// Dwell time
//  fid.p_set::get("Dim", np);		// Number of points
//  fid.p_set::get("dwell",dw);		// Dwell time
  FID(sigma,det,ham,dw,np,fid);		// Call fully parameterized FID
  }



// direct calculation in frequency range


void spec(floq_op sigma, gen_op& det, floq_op& ham, 
          double minF, double maxF, int np, row_vector& spec)

	// Input		sigma : operator propagated (initial density mx)
	// 			det   : "detction operator" in trace computation
	//	 		ham   : "hamiltonian" for propagation (in Hertz)
	//	 		minF,maxF    : Frequency range in (in Hertz)
	//	 		np    : number of points to generate
	//			spec  :	data row_vector containing the result
	// Output		      : none, spec data row_vector filled
	// Note			      : it is assumed row_vector spec is at least
	//				as large as np

{
  spec = row_vector(np);
  ham.set_EBR();	// ham into its eigenbase
  sigma.Op_base(ham);	// density operator into eigenbase of ham

  //  \i = sqrt(-1)
  //  assuming t=0 at beginn of FID
  // spec   = Sum m trace ( det_m * e^\i 2 Pi t m\omega 
  //                        * e^(\i 2 Pi ham t) * sigma * e^(-\i 2 Pi ham t) )
  //        = Sum i,j,m ( <i|det_m|j>e^(\i t m\omega)
  //                      e^(\i 2 Pi <j|ham|j> t)*<j|sigma|i>*e^(-\i 2 Pi <i|ham|i> t)
  //        = Sum i,j,m ( <i|det_m|j><j|sigma|i> * 
  //                  e^(\i 2 Pi t (<j|ham|j>-<i|ham|i>+m \omega)) )
  //        = Sum i,j ( A(i,j) * e^(\i 2 Pi t B(i,j)) )

  double omega = ham.omega();
  double dF    = (maxF-minF)/(np-1);
  int nf = ham.phodim();			// Phonon space dimension
  int n  = ham.size();				// Floquet space dimension
  int i,j,m,pos;

  

  complex A,B;

  for (m=-nf; m<=nf; m++)
    {
      floq_op fdet(nf,ham.hsdim(),omega);
      fdet.put_sdiag (det, m);
//      fdet.put_block(det,0,m);
      fdet.Op_base(ham);
      double f = m*omega;
      for (i=0; i<n; i++)
	for (j=0; j<n; j++)
     
	  {
	    A = fdet.get(i,j)*sigma.get(j,i);
	    B = ham.get(j,j)-ham.get(i,i)-f;
	    pos = (int)((Re(B)-minF)/dF + 0.5);
	    if ((square_norm(A)>1.0e-12)&&
                (pos>=0)&&(pos<np)) 
	           // use elements only when they contribute to the trace
	      spec.put(spec(pos)+A, pos);
	  }
    }
}


void spec_maspowder(floq_op Fsig, gen_op& D, floq_op& H, 
                            double minF, double maxF, int np, row_vector& spec)

        // Input                Fsig  : Operator propagated (initial dens. mx.)
        //                      D     : Detection operator in trace computation
        //                      H     : Hamiltonian for propagation (in Hertz)
	//	 		minF  : Minimum frequency of output spectrum (Hz)
	//			maxF  : Maximum frequency of output spectrum (Hz)
	//	 		np    : Number of points to generate
	//			spec  :	Data vector for the result
	// Output		void  : Block spec is filled with a frequency
	//				domain spectrum from [minF, maxF] that
	//				is due to system Fsig evolution under
	//				Hamiltonian H & detected by D
	// sosi			      : What this has to do with MAS or Powders
	//				is beyond me.  I don't see anything in here
	//				that deals with either of these.  Also the
	//				cutoff should probably first check FDet
	//				then Fsig then pos.  Most fails will be on
	//				FDet which should be sparse (for F+/F-)

  {
  spec = row_vector(np, complex0);		// Set new row vector
  H.set_EBR();					// Put H into an eigenbase
  Fsig.Op_base(H);				// Put sigma in H eigenbase
  double Om = H.omega();			// Fourier expansion frequency
  int NP = H.phodim();				// Phonon space dimension
  int NF = H.size();				// Floquet space dimension
  floq_op FD(NP,H.hsdim(),Om);			// Empty Floquet detection op.
  FD.put_sdiag(D, 0);				// Put in D, fill detection op.
  FD.Op_base(H);				// Insure FD is in H eigenbase
  double dF    = (maxF-minF)/double(np-1);
  int i,j,pos;
  complex A,B;
  for(i=0; i<NF; i++)				// Loop the full Floquet space
    for(j=0; j<NF; j++)
      {
      A   = FD.get(i,j)*Fsig.get(j,i);		// <i|FD|j><j|FD|i>
      B   = H.get(j,j)-H.get(i,i);		// <j|H|j>-<i|H|i> = Wji
      pos = (int)((Re(B)-minF)/dF + 0.5);	// Wji position in spectrum
      if((square_norm(A)>1.0e-12)&&		// Ignore if outside block
	                   (pos>=0)&&(pos<np)) 	// or of limited intensity
      spec.put(spec(pos) + A, pos);		// Add contribution to spectrum
      }
  }


// ____________________________________________________________________________
// F                   DIRECT TRANSITION TABLE GENERATION
// ____________________________________________________________________________

/* These functions return a 1D transitions table given a prepared system
   (density operator) that will be evolved during the an acquisition.  There
   are a few other factors which will influence the table generated by these
   functions.

        ICUT    This is an intensity cutoff value.  Any transition whose
                intensity (norm) is below this value will be removed from the
                returned table.
*/

TTable1D table(const floq_op& Fsig, const gen_op& D, const floq_op& H) 

  {
  double DCUTOFF = 1.e-3;		// Set the detection cutoff level
  H.set_EBR();                          // Put H into an eigenbase
  Fsig.Op_base(H);			// Put sigma in H eigenbase
  double Om = H.omega();		// Get Fourier expansion frequency
  int NP = H.phodim();			// Phonon space dimension
  int NF = H.size();			// Floquet space dimension
  floq_op FD(NP,H.hsdim(),Om);		// Empty Floquet detection operator
  FD.put_sdiag(D, 0);			// Put in D, fill detection operator
  FD.Op_base(H);			// Insure FD is in H eigenbase

/* We need to know how many relevant transitions exist in order to build
   an nice transitions table.  To do that we'll calculate all the transition
   intensities and compare them to the detection cutoff, DCUTOFF.  Since we
   have to calculate them anyway, we'll store the intensities.  Typically
   there will be many less viable transitions than there are Floquet
   coherences (NF*NF), so it is deemed best just to dynamically store these 
   in order to save space at the expense of CPU use.                         */

  int i,j,pos;				// Looping indices
  std::vector<complex> Is;			// Array of intensities
  std::vector<complex> WR;			// Array of frequencies,rates
  complex z;				// Working intensity
  for(i=0, pos=0; i<NF; i++)		// Loop the full Floquet space
    for(j=0; j<NF; j++, pos++)
      {
      z = FD.get(i,j)*Fsig.get(j,i);	// Intensity = <i|D|j><j|sigma|i>
      if(square_norm(z) > DCUTOFF)	// If this value is above the detection
        {				// cutoff, then store the intensity
        Is.push_back(z);		// then calculate and store the
        z = H.get(j,j)-H.get(i,i);	// frequency: <j|H|j>-<i|H|i> = Wji
        WR.push_back(z);
        }
      }

/* Now we will construct and fill up an array of transitions.  This is done
   in a format for generation of a TrnsTable (transitions table).            */

  int ntr = Is.size();			// Total number of viable transitions
  matrix mx(ntr, 2);			// Allocate an array for transitions
  for(pos=0; pos<ntr; pos++)		// Loop through viable transitions
    {					// and put the intensity,frequency
    mx.put(Is[pos], pos, 1);		// info into the array
    mx.put(complex(0,Re(WR[pos])),pos,0);
    }
  return TTable1D(mx);
  }

void transitions(const floq_op& Fsig0, const gen_op& D, const floq_op& H, 
                                   double minF, double maxF, int np, row_vector& spec)
  {
  spec = row_vector(np);		// Set new row vector
  H.set_EBR();                          // Put H into an eigenbase
  Fsig0.Op_base(H);			// Put sigma in H eigenbase
  double Om = H.omega();		// Get Fourier expansion frequency
  int NP = H.phodim();			// Phonon space dimension
  int NF = H.size();			// Floquet space dimension
  floq_op FDet(NP,H.hsdim(),Om);	// Empty Floquet detection operator
  FDet.put_sdiag(D, 0);			// Put in D, fill detection operator
  FDet.Op_base(H);			// Insure FDet is in H eigenbase
//  double dF = (maxF-minF)/(np-1);	// Frequency increment between points
  int i,j,pos;
  complex A,B;
  for(i=0, pos=0; i<NF; i++)		// Loop the full Floquet space
    for(j=0; j<NF; j++, pos++)
      {
      A   = FDet.get(i,j)*Fsig0.get(j,i);// <i|FD|j><j|FD|i>
      B   = H.get(j,j)-H.get(i,i);	// <j|H|j>-<i|H|i> = Wji
      if((square_norm(A)>1.0e-12))
      std::cout << "\n\t\t" << pos << ".   " << A << "   " << B;
      }
  }



void FID_vega(floq_op sigma, gen_op& det, floq_op& ham, double dw, int np, row_vector& fid)

	// Input		sigma : operator propagated (initial density mx)
	// 			det   : "detction operator" in trace computation
	//	 		ham   : "hamiltonian" for propagation (in Hertz)
	//	 		dw    : evolution time (seconds, per point)
	//	 		np    : number of points to generate
	//			fid   :	data row_vector containing the result
	// Output		      : none, fid data row_vector filled
	// Note			      : it is assumed row_vector fid is at least
	//				as large as np

  {
  floq_op FP;		// propagator
  ham.set_EBR();	// ham into its eigenbase
  sigma.Op_base(ham);	// density operator into eigenbase of ham

  FP = exp( complex(0,-2.0*PI*dw)*ham ); // evolution propagator

  //  \i = sqrt(-1)
  //  assuming t=0 at beginn of FID
  // fid(k) = Sum m trace ( det_m * e^\itm\omega 
  //                        * (FP^k) * sigma * (adjoint(FP)^k) )
  //        = Sum i,j,m ( <i|det_m|j>e^(\i(k*dw)m\omega)
  //                      <j|FP|j>^k*<j|sigma|i>*conj(<i|FP|i>)^k )
  //        = Sum i,j,m ( <i|det_m|j><j|sigma|i> * 
  //                  (<j|FP|j>*conj(<i|FP|i>)*e^(\i dw m 2 PI \omega))^k )
  //        = Sum i,j ( A(i,j) * B(i,j)^k )

  double omega = ham.omega();
  int nf = ham.phodim();
  int n  = ham.size();
  int size = (2*nf+1)*n*n;
  int i,j,m,pos;

  

  complex *A = new complex[size];
  complex *B = new complex[size];

  for (pos=0, m=-nf; m<=nf; m++)
    {
      floq_op fdet(nf,ham.hsdim(),omega);
//      fdet.put_sdiag (det, m);
      fdet.put_block(det,0,m);
      fdet.Op_base(ham);
      complex ph=exp(complex(0,dw*m*omega*2*PI));
      for (i=0; i<n; i++)
	for (j=0; j<n; j++)
     
	  {
	    A[pos] = fdet.get(i,j)*sigma.get(j,i);
	    B[pos] = conj(FP.get(i,i),FP.get(j,j))*ph;
	    if ((square_norm(A[pos])>1.0e-12)&&
		(square_norm(B[pos])>1.0e-12)) 
	           // use elements only when they contribute to the trace
	      pos ++;
	  }
    }

  for (int k=0; k<np; k++)
    {
      complex z(0);
      for (int p=0; p<pos; p++)
	{
	  z += A[p];
	  A[p] *= B[p];
	}
      fid.put(z,k);
    } 
//fid.p_set::add("dwell", dw, "Dwell Time");
//fid.p_set::add("SW", 1.0/(2.0*dw), "Spectral Width");

  delete []A;
  delete []B;
}


void FID_vega( floq_op &sigma , gen_op& det, floq_op& ham, row_vector& fid)

	// Input	sigma : operator propagated (initial density mx)
	// 		det   : "detction operator" in trace computation
	//	 	ham   : "hamiltonian" for propagation (in Hertz)
	//		fid   :	data row_vector containing the result
	// Output	      : none, fid data row_vector filled

{
int np = 0;				// Number of points in the fid
int dw = 0;				// Dwell time
//  fid.p_set::get("Dim", np);		// Number of points
//  fid.p_set::get("dwell",dw);		// Dwell time
  FID(sigma,det,ham,dw,np,fid);		// Call fully parameterized FID
}



// direct calculation in frequency range


void spec_vega(floq_op sigma, gen_op& det, floq_op& ham, 
          double minF, double maxF, int np, row_vector& spec)

	// Input		sigma : operator propagated (initial density mx)
	// 			det   : "detction operator" in trace computation
	//	 		ham   : "hamiltonian" for propagation (in Hertz)
	//	 		minF,maxF    : Frequency range in (in Hertz)
	//	 		np    : number of points to generate
	//			spec  :	data row_vector containing the result
	// Output		      : none, spec data row_vector filled
	// Note			      : it is assumed row_vector spec is at least
	//				as large as np

{
  spec = row_vector(np);
  ham.set_EBR();	// ham into its eigenbase
  sigma.Op_base(ham);	// density operator into eigenbase of ham

  //  \i = sqrt(-1)
  //  assuming t=0 at beginn of FID
  // spec   = Sum m trace ( det_m * e^\i 2 Pi t m\omega 
  //                        * e^(\i 2 Pi ham t) * sigma * e^(-\i 2 Pi ham t) )
  //        = Sum i,j,m ( <i|det_m|j>e^(\i t m\omega)
  //                      e^(\i 2 Pi <j|ham|j> t)*<j|sigma|i>*e^(-\i 2 Pi <i|ham|i> t)
  //        = Sum i,j,m ( <i|det_m|j><j|sigma|i> * 
  //                  e^(\i 2 Pi t (<j|ham|j>-<i|ham|i>+m \omega)) )
  //        = Sum i,j ( A(i,j) * e^(\i 2 Pi t B(i,j)) )

  double omega = ham.omega();
  double dF    = (maxF-minF)/(np-1);
  int nf = ham.phodim();
  int n  = ham.size();
  int i,j,m,pos;

  

  complex A,B;

  for (m=-nf; m<=nf; m++)
    {
      floq_op fdet(nf,ham.hsdim(),omega);
//      fdet.put_sdiag (det, m);
      fdet.put_block(det,0,m);
      fdet.Op_base(ham);
      double f = m*omega;
      for (i=0; i<n; i++)
	for (j=0; j<n; j++)
     
	  {
	    A = fdet.get(i,j)*sigma.get(j,i);
	    B = ham.get(j,j)-ham.get(i,i)-f;
	    pos = (int)((Re(B)-minF)/dF + 0.5);
	    if ((square_norm(A)>1.0e-12)&&
                (pos>=0)&&(pos<np)) 
	           // use elements only when they contribute to the trace
	      spec.put(spec(pos) + A, pos);
	  }
    }
}


#endif 								// FloqAcq.cc

