/* LSacquire.h **************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Liouville MR Acquisitions                   Interface 		**
**                                                                      **
**      Copyright (c) 1990, 1991, 1992                                  **
**      Scott Smith                          				**
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
** Capabilities.  The Set of Functions HereinAllows for the Simulation  **
** of Free Induction Decays(FID's) in Spin Hilbert Space & Adds         **
** Abilities Related To Simulated Data Acquistion.                      **
**                                                                      **
*************************************************************************/

#ifndef   LSacq_h_			// Is file already included?
#  define LSacq_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler
#    pragma interface			// This is the interface
#  endif

class row_vector;			// Know GAMMA row vectors
class gen_op;				// Know GAMMA Hilbert space operators
class super_op;				// Know GAMMA Liouville superoperators

#include <GamGen.h>			// Include system specifics
#include <LSLib/LSprop.h>		// Include Liouville space propagators

// ____________________________________________________________________________
// A      Free Induction Decays Under Liouvillian Superoperator Evolution
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

//          Functions Which Take A row_vector As An Argument
//                           (acquire == FID)
//         The overloads here allow for a default value set for N

MSVCDLL void acquire(gen_op& sig0,             gen_op& D,super_op& L,double td,int N,row_vector& fid,double CO=1.e-18);
MSVCDLL void FID(gen_op&     sig0,             gen_op& D,super_op& L,double td,int N,row_vector& fid,double CO=1.e-18);
MSVCDLL void acquire(gen_op& sig0,gen_op& sigf,gen_op& D,super_op& L,double td,int N,row_vector& fid,double CO=1.e-18);
MSVCDLL void FID(gen_op&     sig0,gen_op& sigf,gen_op& D,super_op& L,double td,int N,row_vector& fid,double CO=1.e-18);

MSVCDLL void acquire(gen_op& sig0,             gen_op& D,super_op& L,row_vector& fid,double td,int np=0,double CO=1.e-18);
MSVCDLL void FID(gen_op&     sig0,             gen_op& D,super_op& L,row_vector& fid,double td,int np=0,double CO=1.e-18);
MSVCDLL void acquire(gen_op& sig0,gen_op& sigf,gen_op& D,super_op& L,row_vector& fid,double td,int np=0,double CO=1.e-18);
MSVCDLL void FID(gen_op&     sig0,gen_op& sigf,gen_op& D,super_op& L,row_vector& fid,double td,int np=0,double CO=1.e-18);
    
//                  Functions Which Return A row_vector 
//                           (acquire == FID)

MSVCDLL row_vector acquire(gen_op& sig0,             gen_op& D,super_op& L,double td,int np,double CO=1.e-18);
MSVCDLL row_vector FID(gen_op&     sig0,             gen_op& D,super_op& L,double td,int np,double CO=1.e-18);
MSVCDLL row_vector acquire(gen_op& sig0,gen_op& sigf,gen_op& D,super_op& L,double td,int np,double CO=1.e-18);
MSVCDLL row_vector FID(gen_op&     sig0,gen_op& sigf,gen_op& D,super_op& L,double td,int np,double CO=1.e-18);

// ____________________________________________________________________________
// B      Free Induction Decays Under Superoperator Propagator Evolution
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

MSVCDLL void acquire(gen_op& sig, gen_op& D, super_op& G, row_vector &data, int np=0);
MSVCDLL void FID(gen_op&     sig, gen_op& D, super_op& G, row_vector &data, int np=0);
MSVCDLL void acquire(gen_op& sig, gen_op& D, LSprop&   G, row_vector &data, int np=0);
MSVCDLL void FID(gen_op&     sig, gen_op& D, LSprop&   G, row_vector &data, int np=0);
    
//                  Functions Which Return A row_vector 
//                           (acquire == FID)

MSVCDLL row_vector acquire(gen_op& sig, gen_op& D, super_op& G, int np);
MSVCDLL row_vector FID(gen_op&     sig, gen_op& D, super_op& G, int np);
MSVCDLL row_vector acquire(gen_op& sig, gen_op& D, LSprop&   G, int np);
MSVCDLL row_vector FID(gen_op&     sig, gen_op& D, LSprop&   G, int np);

// ____________________________________________________________________________
// C      Other Liouville Space Acquisition Related Functions
// ____________________________________________________________________________

// Older that need updating 

MSVCDLL void FIDx(gen_op sigma, gen_op& sigma0, gen_op& det, super_op &L,
	   			   row_vector &fid, double dt, int np=0);
MSVCDLL void FIDrot(gen_op sigma, gen_op& sigma0, gen_op& det, super_op &L, gen_op& Fz,
                double Wrf, double time, row_vector &fid, double dt, int np=0);

#endif 						// LSacquire.h
