/* LSacquire.cc *************************************************-*-c++-*-
**                                                                      **
**                            G A M M A                                 **
**                                                                      **
**      Liouville Space Acquisitions               Implementation       **
**                                                                      **
**      Copyright (c) 1990, 1991, 1992                                  **
**      Tilo Levante, Scott Smith                                       **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The GAMMA Platform Provides Functions for Simulation of Magnetic     **
** Resonance Experiments and Other Associated   Mathematical            **
** Capabilities.  The Set of Functions Herein Allow for the Simulation  **
** of Free Induction Decays(FID's) in Spin Liouville Space & Add	**
** Abilities Related To Simulated Data Acquistion.                      **
**                                                                      **
*************************************************************************/

#ifndef _LSacq_cc_			// Is file already included?
#define _LSacq_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#endif

#include <LSLib/LSacquire.h>		// Include the file header
#include <Matrix/row_vector.h>		// Knowledge of row vectors
#include <HSLib/GenOp.h>		// Knowledge of operators
#include <LSLib/SuperOp.h>		// Knowledge of superoperators
#include <HSLib/HSprop.h>		// Knowledge of nmr evolutions

// ____________________________________________________________________________
// A                Generic Liouville Space Acquisition Functions
// ____________________________________________________________________________

/* The mathematics for generating a series of expectation values that evolve
   through each step under an exponential of an input superoperator is sketched
   out below. Here L is the superoperator, sigma the initial density operator,
   siginf the density operator for the system after infinite time evolution, &
   D the operator representing the properties we want to sample. It is assumed
   that if the infinite time density operator is specified then the eponential
   of L acts on the difference density operator & that there is then a need for
   semi-classical difference density operator terms. If there is no infinite
   time density operator is supplied then the exponential of L acts directly on
   the density operator (same equations with siginf set to zero)
 
                            -Lt                      
         fid(k) = <adj(D)| e   k |sigma - siginf> + <adj(D)|siginf>

                             -Et    -1                  
         fid(k) = <adj(D)*S|e   k |S  * (sigma - siginf)> + <adj(D)|siginf>

                    -1
    where  L = S*E*S   with E being a diagonal superoperator containing the
    eigenvalues of L and S a matrix containing the the eigenvectors of L.
    Letting delsig = sigma-siginf, K=<adj(D)|siginf>, DAS = adj(D)*S we obtain
    (where E is a diagonal superoperator)

                                    -Et    -1                  
                     fid(k) = <DAS|e   k |S  * delsig> + K

    Expanding this, where the index l spans Liouville space, we obtain

                       ---
                       \              -Et        -1                  
              fid(k) = /   <DAS|l><l|e   k|l><l|S  * delsig> + K
                       ---
                        l

     Using td for the time between fid pts (expectation values) & setting the
     elements of the eigenvalue array to <l|E|l> = -iw  + R   we get the form
                                                      ll   ll
                          ---
                          \              -1          (-iw  + R  )k*t      
             fid(k) = K + /   <DAS|l><l|S  * delsig>e    ll   ll    d
                          ---
                           l
                                                           -1
     This simplifies a bit if we make a vector <A| = <DAS|l><l|S  >
     and restrict sum over non-zero <A| elements.

                                ---                    ---
                                \             k        \       k
                   fid(k) = K + /   <A|p><p|B>   = K + /   A  B 
                                ---                    ---  p  p
                                 p                      p
 
     The FID functions which take a superoperator and a time as arguments all
     implement the last form when generating a series of expectation values
     evenly incremented in time. The routines will immediately generate the 
     vectors A and B then quicky do the summation over p. This process is
     repeated to generate the acqusition points by taking powers of B elements.
     When there is not siginf specified then K will be zero and delsig = sigma.

	Input		sig  : Operator propagated (initial dens. mx)
	                sigf : Equilibrium/Steady-State density matrix
	   		D    : Detection operator in trace computation
	   		L    : Liouvillian for propagation (in rad/sec)
		        U    : Propagator for one increment
	   		td   : Evolution time (seconds, per point)
	   		N    : Number of points to generate
	  		data : Data vector containing the result
	Output		     : none, fid data vector filled
	Note		     : If no block size is specified it will
			       be taken as the vector size
	Note		     : If no vector is specified one will be
			       generated and returned
	Note		     : This routine is typically used where
       			       the Liouvillian is

	  		          L = iH + R     &     [H, R] != 0          */

//                      Generic Acquisiton Functions
//                    (acquire == FID, No Propagators)

// sosi - why is this so lame? It does not save much time in generating anything....
// and does not follow the prescription given above for doing the math either!

void FIDL(gen_op& sigma,gen_op& siginf,gen_op& det,super_op& L,row_vector& fid,double dt,int np,double CO)
  {
//sosi - still need to fix this
CO += 1.0;
  if(np <= 0) np = fid.size();		// Set up block size if not specified
//  complex K = trace(det*siginf);	// Find the K value <det|siginf>
  gen_op delsig = sigma - siginf;	// Set up the difference density matrix
  L.set_EBR();				// Put superoperator into its eigenbase
  L.LOp_base(delsig);			// Put delta sigma into LOp Hilbert basis
  L.LOp_base(siginf);			// Put siginf into LOp Hilbert basis
  L.LOp_base(det);			// Put det into LOp Hilbert basis
  basis S = L.get_Lbasis();		// Get superoperator basis
  matrix Sm1 = inv(S.U());		// Get superoperator basis inverse
  matrix mx  = delsig.get_mx();		// Get delta sigma as array
  matrix Sm1delsig =			// Get |inv(S)*delsigma> (superket)
            Sm1*mx.resize(L.size(), 1);
matrix mxd;

  double tk;				// Evolution time from initial point
  super_op emDtk;			// For L Exponential eigenvalues
  for(int k=0; k<np; k++)
    {
    tk = double(k)*dt;			// Get FID time at point k
    emDtk = exp(L, -tk);		// Compute exponential at that time
// sosi - there is an unexplained - sign taken out here. It should be in so
//        there must be one missing somewhere else.
//    emDtk = exp(L, -tk);		// Compute exponential at that time
//if(k==1) cout << emDtk.get_mx();
    mx = (S.U() * emDtk.get_mx()) * Sm1delsig;
    mx = mx.resize(sigma.dim(), sigma.dim());
    mx += siginf.get_mx();		// Reform full relaxed density matrix
mxd = adjoint(det.get_mx());
fid.put(trace(mxd,mx), k);
//    fid.put(trace(det.get_mx(), mx), k);// Compute FID point k
    }
  }

void FIDL(gen_op& sigma,               gen_op& det, super_op& L,row_vector& fid,double dt,int np,double CO)
  { gen_op siginf; FIDL(sigma,siginf,det,L,fid,dt,np,CO); }

row_vector FIDL(gen_op& sig,gen_op& sigf,gen_op& D,super_op& L,double dt,int np,double CO)
  {
  row_vector data(np);
  FIDL(sig, sigf, D, L, data, dt, np, CO);
  return data;
  }

row_vector FIDL(gen_op& sig,             gen_op& D,super_op& L,double dt,int np,double CO)
  { gen_op sigi; return FIDL(sig,sigi,D,L,dt,np,CO); }

// ____________________________________________________________________________
// B       Acquisition Functions Under Constant Liouvillian Evolution
// ____________________________________________________________________________

//          Functions Which Take A row_vector As An Argument
//                   (acquire == FID, No Propagators)

void acquire(gen_op& sig,             gen_op& D,super_op &L,double dt,int np,row_vector &data,double CO)
  { FIDL(sig,D,L,data,dt,np,CO); }
void FID(    gen_op& sig             ,gen_op& D,super_op &L,double dt,int np,row_vector &data,double CO)
  { FIDL(sig,D,L,data,dt,np,CO); }
void acquire(gen_op& sig,gen_op& sigi,gen_op& D,super_op &L,double dt,int np,row_vector &data,double CO)
  { FIDL(sig,sigi,D,L,data,dt,np,CO); }
void FID(    gen_op& sig,gen_op& sigi,gen_op& D,super_op &L,double dt,int np,row_vector &data,double CO)
  { FIDL(sig,sigi,D,L,data,dt,np,CO); }

void acquire(gen_op& sig,             gen_op& D,super_op &L,row_vector &data,double dt,int np,double CO)
  { FIDL(sig,D,L,data,dt,np,CO); }
void FID(    gen_op& sig             ,gen_op& D,super_op &L,row_vector &data,double dt,int np,double CO)
  { FIDL(sig,D,L,data,dt,np,CO); }
void acquire(gen_op& sig,gen_op& sigi,gen_op& D,super_op &L,row_vector &data,double dt,int np,double CO)
  { FIDL(sig,sigi,D,L,data,dt,np,CO); }
void FID(    gen_op& sig,gen_op& sigi,gen_op& D,super_op &L,row_vector &data,double dt,int np,double CO)
  { FIDL(sig,sigi,D,L,data,dt,np,CO); }

//                  Functions Which Return A row_vector 
//                    (acquire == FID, No Propagators)

row_vector acquire(gen_op& sig,             gen_op& D,super_op& L,double dt,int np,double CO)
  { return FIDL(sig,D,L,dt,np,CO); }
row_vector FID(    gen_op& sig,             gen_op& D,super_op& L,double dt,int np,double CO)
  { return FIDL(sig,D,L,dt,np,CO); }
row_vector acquire(gen_op& sig,gen_op& sigf,gen_op& D,super_op& L,double dt,int np,double CO)
  { return FIDL(sig,sigf,D,L,dt,np,CO); }
row_vector FID(    gen_op& sig,gen_op& sigf,gen_op& D,super_op& L,double dt,int np,double CO)
  { return FIDL(sig,sigf,D,L,dt,np,CO); }

// ____________________________________________________________________________
// C  Acquisition Functions Under Constant Superoperator Propagator Evolution
// ____________________________________________________________________________

/* The mathematics for generating a series of expectation values that evolve
   through each step under a superoperator propagator G is sketched out below.
   Here G is the propagator, sigma the initial density operator and D the
   operator representing the properties we want to sample. It is assumed that
   G acts directly on the density operator with no need for the semi-classical
   difference density operator terms.
                                            
                              k                     k  -1                  
           fid(k) = <adj(D)| G |sigma> = <adj(D)*S|E |S  * sigma>

                -1
   where  G = S*E*S and E is a diagonal superoperator. Letting DAS = adj(D)*S 
   we expand the equation where the index l spans Liouville space.

                          ---
                          \              k       -1                  
                 fid(k) = /   <DAS|l><l|E |l><l|S  * sigma>
                          ---
                           l

                                                           -1
   This simplifies a bit if we make a vector <A| = <DAS|l><l|S  *sigma>
   and restrict sum over non-zero <A| elements.

                             ---                ---
                             \             k    \       k
                    fid(k) = /   <A|p><p|E>   = /   A  E
                             ---                ---  p  p
                              p                  p
 
   The FID functions all implement the last form when generating a series
   of expectation values evenly incremented in time when give the propagator
   G as an argument.  The routines will immediately generate the vectors A and 
   E then quicky do the sum. The fid points are repeatly generated by taking 
   powers of E elements.

	Input		sig  : Operator propagated (initial dens. mx)
	                sigf : Equilibrium/Steady-State density matrix
	   		D    : Detection operator in trace computation
		        G    : Propagator for one increment
	   		N    : Number of points to generate
	  		data : Data vector containing the result
	Output		     : none, fid data vector filled
	Note		     : If no block size is specified it will
			       be taken as the vector size
	Note		     : If no vector is specified one will be
			       generated and returned
	Note		     : This routine is typically used where
       			       the superoperator propagator G acts a

	  		          |sigma(td) = G|sigma>                      */


//          Functions Which Take A row_vector As An Argument
//                           (acquire == FID)

void acquire(gen_op& sig, gen_op& D, super_op& G, row_vector& data, int np)
  { FID(sig,D,G,data,np); }

// sosi - I had to fix this, the others may still be bad!
void FID(gen_op& sigma,   gen_op& D, super_op& G, row_vector& data, int np)
  {
  if(np <= 0) np = data.size();		// Set up block size if not specified
  G.set_EBR();				// Put superoperator into its eigenbase
  G.LOp_base(sigma);			// Put sigma into LOp Hilbert basis
  G.LOp_base(D);			// Put D into LOp Hilbert basis
  basis S    = G.get_Lbasis();		// Get superoperator basis
  matrix Sm1 = inv(S.U());		// Get superoperator basis inverse
  matrix mx  = sigma.get_mx();		// Get sigma as a matrix
  int ls     = G.size();		// Get Liouville space dimension
  int hs     = sigma.dim();		// Get Hilbert space dimension
  matrix Sm1sig = Sm1*mx.resize(ls,1);	// Determine inv(S)*sigma as matrix
  matrix Gmx = G.get_mx();		// Get the G array (diagonal)
  matrix Umx(ls,ls,i_matrix_type);
  matrix Dmx = D.get_mx();		// Get detection operator array
  for(int k=0; k<np; k++)		// Loop over acquisition points
    {
    mx = (S.U()*Umx)*Sm1sig;		// Evolve sigma to this point
    mx = mx.resize(hs,hs);		// Set back into square array
    data.put(trace(Dmx, mx), k);	// Compute FID point
    Umx *= Gmx;				// Set Gmx for evolve to next point
    }
  }

void FID(gen_op& sigma,   gen_op& D, LSprop& G, row_vector& data, int np)
  { super_op GOP = G.LOp(); FID(sigma,D,GOP,data,np); }

//                  Functions Which Return A row_vector 
//                           (acquire == FID)

row_vector acquire(gen_op& sig, gen_op& D, super_op &G, int np)
  { return FID(sig, D, G, np); }

row_vector FID(gen_op&    sig, gen_op& D, super_op &G, int np)
  { row_vector data(np,complex0); FID(sig, D, G, data, np); return data; }

row_vector acquire(gen_op&    sig, gen_op& D, LSprop& G, int np)
  { return FID(sig, D, G, np); }

row_vector FID(gen_op&    sig, gen_op& D, LSprop& G, int np)
  { row_vector data(np,complex0); super_op GOP = G.LOp(); FID(sig, D, GOP, data, np); return data; }

// ____________________________________________________________________________
// D      Other Liouville Space Acquisition Related Function
// ____________________________________________________________________________

// Older that need updating 

	// Input                sigma : Initial density mx
	//                      sigma0: Equilibrium/Steady-State density matrix
	//                      det   : "detection operator" in trace comp.
	//                      L     : Liouvillian for relaxation & evolution
	//                      dt    : Evolution time (seconds, per point)
	//                      np    : Number of points to generate
	//                      fid   : Data row_vector containing the result
	// Output                     : None, fid data row_vector filled
	// Note			      : The relaxation superoperator should be
	//				a Liouvillian, i.e. it must contain the
	//				time evolution components.
	// Note                       : It is assumed row_vector fid is at least
	//                              as large as np
	// Note                       : Routine is used when Liouvillian is
	//		                   L = iH + R     &     [H, R] != 0

  // fid(k) = Trace { det *  sigma(tk) },  where tk = td*k
  //        = Trace { det *  [[exp(-L*tk) * (sigma - sigma0)] + sigma0] }
  //        = Trace { det *  [[exp(-L*tk) * delsigma] + sigma0] }
  //        = Trace { det *  [[exp(-L*dt) * exp(-L*dt(k-1)*delsigma] + sigma0] }
  //        = Trace { det *  [[exp(-L*dt) * delsigma(k-1)] + sigma0] }
  //        = trace { det * [[ S * exp(-D*dt) * S-1 * delsigma(k-1)] + sigma0] }

void FIDx(gen_op sigma, gen_op& sigma0, gen_op& det, super_op &L,
	   			             row_vector &fid, double dt, int np)
    {
    L.LOp_base(sigma);			// Put density op. in Hilbert basis of L
    L.LOp_base(sigma0);			// (Eq./St.State) dens. op. in Hilbert basis of L
    L.LOp_base(det);			// Put detect. op. inHilbert basis of L
    super_op eLdt = exp(L, -dt);	// Exponential superop. for dwell time
    eLdt.set_HBR();			// Put Exp. superop into default basis
    gen_op delsigma = sigma-sigma0;	// Difference density matrix
    gen_op sigmat = sigma;		// Density matrix at different times
    complex z;
    for(int i=0; i<np; i++)		// Sum over all FID points
      {
      if(i != 0)
        {
        delsigma = eLdt * delsigma;
        sigmat = delsigma + sigma0;
        }
      z = trace(sigmat, det);
      fid.put(z, i); 
      }
    }


void FIDrot(gen_op sigma, gen_op& sigma0, gen_op& det, super_op &L, gen_op& Fz,
	                double Wrf, double time, row_vector &fid, double dt, int np)
    
	// Input                sigma : Initial density mx
	//                      sigma0: Equilibrium/Steady-State density matrix
	//                      det   : The "detection operator" in trace computation
	//                      L     : Liouvillian for relaxation & evolution
	//                      Fz    : Appropriate Fz operator
	//                      Wrf   : Rotating frame frequency
	//                      fid   : Data row_vector containing the result
	//                      dt    : Evolution time (seconds, per point)
	//                      np    : Number of points to generate
	// Output                     : None, fid data row_vector filled
	// Note                       : It is assumed row_vector fid is at least
	//                              as large as np
	// Note                       : This routine is used where the Liouvillian is
	//		                   L = iH + R     &     [H, R] != 0

  // fid(k) = Trace { det *  sigma(tk) },  where tk = td*k
  //        = Trace { det *  [[exp(-L*tk) * (sigma - sigma0)] + sigma0] }
  //        = Trace { det *  [[exp(-L*tk) * delsigma] + sigma0] }
  //        = Trace { det *  [[exp(-L*dt) * exp(-L*dt(k-1) * delsigma] + sigma0] }
  //        = Trace { det *  [[exp(-L*dt) * delsigma(k-1)] + sigma0] }
  //        = trace { det * [[ S * exp(-D*dt) * S-1 * delsigma(k-1)] + sigma0] }

    {
    sigma = evolve(sigma, Fz, -Wrf*time);	// Put sigma in Wrf rotating frame
    L.LOp_base(sigma);				// Put sigma into Hilbert basis of L
    L.LOp_base(sigma0);				// (Eq./steady-state) dens.op. to H-basis of L
    L.LOp_base(det);				// Put detection operator into Hilbert basis of L
    super_op eLdt = exp(L, -dt);		// Form exponential superoperator for dwell time
    L.set_HBR();				// Put L into Liouvillian default basis
    gen_op delsigma = sigma - sigma0;		// Difference density matrix (Wrf rot. frame)
    gen_op sigmat = sigma;			// sigma matrix at different times
    complex z;
    for(int k=0; k<np; k++)			// Loop over all FID points
      {
      if(k != 0)
        {
        delsigma = eLdt * delsigma;		// delsigma(k) from delsigma(k-1)
        sigmat = delsigma + sigma0;		// sigma(k) in Wrf rotating frame
        time += dt;
        }
      sigmat = evolve(sigmat,Fz, Wrf*time);	// sigma(k) in initial rotating frame
      z = trace(sigmat, det);			// FID(k) from Tr{sigma(k)*det}
      fid.put(z, k); 
      }
    }


//void FID(gen_op& sigma, gen_op& sigma0, gen_op& det, gen_op& H,
//                        super_op &R, row_vector &fid, double dt, int np=0)

  // fid(k) = trace { det * sigma(tk) }
  // fid(k) = trace { det * exp(-iHtk) * [[exp(-Rtk) * (sigma - sigma0)] + sigma0] * exp(iHtk) }
  // fid(k) = trace { det * exp(-iHsuptk) * [[exp(-Rtk) * (sigma - sigma0)] + sigma0] }

	// Input		sigma : Operator propagated (initial density mx)
	//			sigma0: Operator at equilibrium (or steady state)
	// 			det   : Detection operator (in trace computation)
	//	 		ham   : Hamiltonian for propagation (in Hertz)
	//			R     : Relaxation superoperator (rad/sec)
	//			fid   :	Data row_vector containing the result
	//	 		dt    : Evolution time (seconds, per point)
	//	 		np    : Number of points to generate
	// Output		      : None, FID data row_vector filled
	// Note			      : It is assumed row_vector fid is at least
	//				as large as np
	// Note			      : The relaxation superoperator should not
	//				contain the time precession and should commute
	//				with the Hamiltonian.

  //{
  //super_op Hsup = Hsuper(H);		// Set up Hamiltonian superoperator (rad/sec)
  //FID(sigma,sigma0,det,Hsup,R,fid,dt,np)// Use the overloaded superop function
  //return;
  //}


//void FID(gen_op& sigma, gen_op& sigma0, gen_op& det, super_op& H,
                        //super_op &R, row_vector &fid, double dt, int np=0)

  // fid(k) = trace { det * sigma(tk) }
  // fid(k) = trace { det * exp(-iHtk) * [[exp(-Rtk) * (sigma - sigma0)] + sigma0] }
  //        = trace { det * exp(-iHtk) * [[ S * exp(-Dtk) * S-1 * (sigma - sigma0)] + sigma0] }

	// Input		sigma : Operator propagated (initial density mx)
	//			sigma0: Operator at equilibrium (or steady state)
	// 			det   : Detection operator (in trace computation)
	//	 		H     : Hamiltonian superoperator (rad/sec)
	//			R     : Relaxation superoperator (rad/sec)
	//			fid   :	Data row_vector containing the result
	//	 		dt    : Evolution time (seconds, per point)
	//	 		np    : Number of points to generate
	// Output		      : None, FID data row_vector filled
	// Note			      : It is assumed row_vector fid is at least
	//				as large as np
	// Note			      : The relaxation superoperator should not
	//				contain the time precession and should commute
	//				with the Hamiltonian superoperator.

  //{
  //if(np <= 0) np = fid.size();		// Set up the block size if not specified
  //gen_op delsig = sigma - sigma0;	// Set up the difference density matrix
  //L.set_EBR();				// Put superoperator into its eigenbase
  //H.set_EBR();				// Put hamiltonian into its eigenbase
  //L.LOp_base(delsig);			// Put delta sigma into LOp Hilbert basis
  //L.LOp_base(sigma0);			// Put sigma0 into LOp Hilbert basis
  //L.LOp_base(det);			// Put det into LOp Hilbert basis
  //matrix S = L.get_Lbasis();		// Get superoperator basis
  //matrix Sm1 = inv(S);			// Get superoperator basis inverse
  //matrix mx = delsig.get_mx();
  //matrix Sm1delsig = Sm1*mx.resize(L.size(), 1);
  //double tk;
  //super_op emDtk, emHtk;
  //for(int k=0; k<np; k++)
    //{
    //tk = double(k)*dt;
    //emHtk = exp(L, -tk);
    //emDtk = exp(L, -tk);
    //mx = (S * emDtk.get_mx()) * Sm1delsig;
    //mx = exp
    //mx = mx.resize(sigma.dim(), sigma.dim());
    //mx += sigma0.get_mx();
    //fid(k) = trace(det.get_mx(), mx);
    //}
  //}

#endif 						// LSacquire.cc
