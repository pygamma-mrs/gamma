/* relaxJ.h *************************************-*-c++-*-
**							**
**	               G A M M A			**
**						 	**
**	Spectral Density Functions    Interface		**
**						 	**
**	Copyright (c) 1990, 1991, 1992, 1993	 	**
**	Scott Smith				 	**
**	Eidgenoessische Technische Hochschule	 	**
**	Labor fuer physikalische Chemie		 	**
**	8092 Zurich / Switzerland		 	**
**						 	**
**      $Header: $
**						 	**
*********************************************************/

/*********************************************************
**						 	**
** 	Description				 	**
**						 	**
**	These are a collection of routines which	**
**	compute spectral density functions or		**
**	quantities related to spectral density		**
**	functions.					**
**						 	**
*********************************************************/

#ifndef   Relax_J_h_			// Is file already included?
#  define Relax_J_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Level1/coord.h>		// Include coordinates
#include <Level1/SpaceT.h>		// Include spatial tensors
#include <HSLib/GenOp.h>		// Include general operators


/********************** Error handling *******************/


void J_error (int i);

	// Input		int  : Error Flag
	// Output		none : Error Message Output



void volatile J_fatality (int error);

	// Input		none :
	// Output		none : Stops Execution & Error Message Output


// ______________________________________________________________________
//                      General Spectral Density Function
// ______________________________________________________________________


MSVCDLL double J_gen(double tau, double omega, int hertz=0);

	// Input		tau  : correlation time (sec)
	//	 		omega: frequency, (Hz or rad/sec)
	//			rad  : flag for Hz or rad/sec
	// Output		J(w) : spectral density value at omega
	// 			       for isotropic rigid diffusion
	// Note		             : output units match tau	


MSVCDLL matrix J_gen(double tau, double* w, int hs, int hertz=0);

	// Input		tau  : correlation time (sec)
	//	 		w    : vector of frequenies, (Hz or 1/sec)
	//			hs   : Hilbert space size
	//			hertz: flag for Hz or rad/sec
	// Output		J    : matrix of spectral density values
	//			       at the transition frequencies between
	// 			       the values contained in the vector w
	// Note		             : output units match tau	


MSVCDLL matrix J_gen_shft(double tau, double* w, double shift, int hs, int hertz=0);

	// Input		tau  : correlation time (sec)
	//	 		w    : vector of frequenies, (Hz or 1/sec)
	//			shift: transition frequency shift
	//			hs   : Hilbert space size
	//			hertz: flag for Hz or rad/sec
	// Output		J    : matrix of spectral density values
	//			       at the transition frequencies between
	// 			       the values contained in the vector w
	//			       shifted by the amount shift
	// Note		             : output units match tau	

// ______________________________________________________________________
//                  Anisotropic Rotor Auxiliary Functions 
// ______________________________________________________________________
//

MSVCDLL void tausD(double* tau, double Dx, double Dy, double Dz);

	// Input		tau  : double vector for tau values
	// 			Dx   : x-axis rot. diffusion const. (rad/sec)
	// 			Dy   : y-axis rot. diffusion const. (rad/sec)
	// 			Dz   : z-axis rot. diffusion const. (rad/sec)
	// Output		tau  : five effective correlation times (sec)
	// 			       for isotropic rigid diffusion


MSVCDLL void tausD(double* tau, const coord& Ds);

	// Input		tau  : double vector for tau values
	// 			Ds   : rotational diffusion consts. (rad/sec)
	// Output		tau  : five effective correlation times (sec)
	// 			       for isotropic rigid diffusion


MSVCDLL double chiD(double Dx, double Dy, double Dz);

	// Input		Dx   : x-axis rot. diffusion const. (rad/sec)
	// 			Dy   : y-axis rot. diffusion const. (rad/sec)
	// 			Dz   : z-axis rot. diffusion const. (rad/sec)
	// Output		chi  : chi value for reduced spectral density
	// 			       function coefficient computations

MSVCDLL double chiD(const coord& Ds);

	// Input		Ds   : rotational diffusion consts. (rad/sec)
	// Output		chi  : chi value for reduced spectral density
	// 			       function coefficient computations


MSVCDLL void taust(double* tau, double taux, double tauy, double tauz);

	// Input		tau  : double vector for tau values
	// 			taux : x-axis correlation time (sec)
	// 			tauy : y-ayis correlation time (sec)
	// 			tauz : z-azis correlation time (sec)
	// Output		tau  : five effective correlation times (sec)
	// 			       for isotropic rigid diffusion


MSVCDLL void taust(double* tau, const coord& taus);

	// Input		tau  : double vector for tau values
	// 			taus : correlation times (sec)
	// Output		tau  : five effective correlation times (sec)
	// 			       for isotropic rigid diffusion


MSVCDLL double chit(double taux, double tauy, double tauz);

	// Input		taux : x-axis correlation time (sec)
	// 			tauy : y-ayis correlation time (sec)
	// 			tauz : z-azis correlation time (sec)
	// Output		chi  : chi value for reduced spectral density
	// 			       function coefficient computations


MSVCDLL double chit(const coord& taus);

	// Input		taus : correlation times (sec)
	// Output		chi  : chi value for reduced spectral density
	// 			       function coefficient computations


MSVCDLL void Jcoeffs(double* c, double alpha, double beta, double gamma,
						    double chi=0, double eta=0);

	// Input		c    : double vector for 5 coefficients
	// 			alpha: Euler angle (radians)
	// 			beta : Euler angle (radians)
	// 			gamma: Euler angle (radians)
	//			chi  : Ratio of rotational correlation times
	// 			eta  : Tensor value
	// Output		c    : five spectral density function coefficients
	// 			       for isotropic rigid diffusion


MSVCDLL void Jcoeffs(double* c, const coord& EA, double chi=0, double eta=0);

	// Input		c    : double vector for 5 coefficients
	// 			EA   : Euler angles (radians)
	//			chi  : Ratio of rotational correlation times
	// 			eta  : Tensor value
	// Output		c    : five spectral density function coefficients
	// 			       for isotropic rigid diffusion


// ______________________________________________________________________
//             Anisotropic Rotor Spectral Density Functions 
// ______________________________________________________________________
//

MSVCDLL double J_reduced(double *tau, double *c1, double *c2, double omega, int hertz=0);

	// Input		tau  : five correlation times (sec)
	// 			c1   : five coefficients of 1
	// 			c2   : five coefficients of 2
	//	 		omega: frequency, (Hz or rad/sec)
	//			rad  : flag for Hz or rad/sec
	// Output		J(w) : spectral density value at omega
	// 			       for isotropic rigid diffusion
	// Note		             : output units match tau	


MSVCDLL matrix J_reduced(gen_op& Op, double *tau, double *c1, double *c2, int hertz=0);

	// Input		Op   : general operator
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transition
	// 			       frequencies of Op
	// Note		             : output units match tau	


MSVCDLL matrix J_reduced(double* w, int size, double *tau, double *c1, double *c2, int hertz=0);

	// Input		w    : double vector of frequencies
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transitions
	// 			       between the frequencies in w
	// Note		             : output units match tau	


MSVCDLL matrix J_red_shft(double* w, double shift, int size, double *tau,
						 double *c1, double *c2, int hertz=0);

	// Input		w    : double vector of frequencies (energies)
	//			shift: transition frequency shift
	//			size : Hilbert space size
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transitions
	// 			       between the frequencies in w
	//			       shifted by the amount shift
	// Note		             : output units match tau	


MSVCDLL matrix J_reduced(row_vector w, double *tau, double *c1, double *c2, int hertz=0);

	// Input		w    : row vector of frequencies
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transitions
	// 			       between the frequencies in w
	// Note		             : output units match tau	


MSVCDLL matrix J_red_shft(row_vector w, double shift, double *tau,
					 double *c1, double *c2, int hertz=0);

	// Input		w    : row vector of frequencies (energies)
	//			shift: transition frequency shift
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transitions
	// 			       between the frequencies in w
	// Note		             : output units match tau	


MSVCDLL double J_reduced(const coord& taus, const coord& EA1, double eta1,
			 const coord& EA2, double eta2, double omega, int hertz=0);

	// Input		taus : correlation times (sec)
	//	 		EA1  : Euler angles mu1 (radians)
	//	 		eta1 : eta value mu1
	//	 		EA2  : Euler angles mu2 (radians)
	//	 		eta2 : eta value mu2
	//	 		omega: frequency, (Hz or rad/sec)
	//			rad  : flag for Hz or rad/sec
	// Output		J(w) : spectral density value at omega
	// 			       for isotropic rigid diffusion
	// Note		             : output units match tau	


MSVCDLL double J_reduced(const coord& taus, space_T& T1, space_T& T2, double omega, int hertz=0);

	// Input		taus : correlation times (sec)
	//			T1   : spatial tensor mu1
	//			T2   : spatial tensor mu2
	//	 		omega: frequency, (Hz or rad/sec)
	//			rad  : flag for Hz or rad/sec
	// Output		J(w) : spectral density value at omega
	// 			       for isotropic rigid diffusion
	// Note		             : output units match tau	


// ______________________________________________________________________
// Anisotropic Rotor Spectral Density Functions - Dynamic Frequency Shift 
// ______________________________________________________________________
//

MSVCDLL double Q_reduced(double *tau, double *c1, double *c2, double omega, int hertz=0);

	// Input		tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//	 		omega: frequency, (Hz or rad/sec)
	//			rad  : flag for Hz or rad/sec
	// Output		J(w) : spectral density value at omega
	// 			       for isotropic rigid diffusion
	// Note		             : output units match tau	



MSVCDLL matrix Q_reduced(gen_op& Op, double *tau, double *c1, double *c2, int hertz=0);

	// Input		Op   : general operator
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transition
	// 			       frequencies of Op
	// Note		             : output units match tau	


MSVCDLL matrix Q_reduced(double* w, int size, double *tau, double *c1, double *c2, int hertz=0);

	// Input		w    : double vector of frequencies (energies)
	//			size : Hilbert space size
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transitions
	// 			       between the frequencies in w
	// Note		             : output units match tau	


MSVCDLL matrix Q_red_shft(double* w, double shift, int size, double *tau,
						 double *c1, double *c2, int hertz=0);

	// Input		w    : double vector of frequencies (energies)
	//			shift: transition frequency shift
	//			size : Hilbert space size
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transitions
	// 			       between the frequencies in w
	// Note		             : output units match tau	


MSVCDLL matrix Q_reduced(row_vector w, double *tau, double *c1, double *c2, int hertz=0);

	// Input		w    : row vector of frequencies (energies)
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transitions
	// 			       between the frequencies in w
	// Note		             : output units match tau	


MSVCDLL matrix Q_red_shft(row_vector w, double shift, double *tau,
					 double *c1, double *c2, int hertz=0);

	// Input		w    : row vector of frequencies (energies)
	//			shift: transition frequency shift
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transitions
	// 			       between the frequencies in w
	// Note		             : output units match tau	


MSVCDLL double Q_reduced(const coord& taus, const coord& EA1, double eta1, 
                   const coord& EA2, double eta2, double omega, int hertz=0);

	// Input		taus : correlation times (sec)
	//	 		EA1  : Euler angles mu1 (radians)
	//	 		eta1 : eta value mu1
	//	 		EA2  : Euler angles mu2 (radians)
	//	 		eta2 : eta value mu2
	//	 		omega: frequency, (Hz or rad/sec)
	//			rad  : flag for Hz or rad/sec
	// Output		J(w) : spectral density value at omega
	// 			       for isotropic rigid diffusion
	// Note		             : output units match tau	


MSVCDLL double Q_reduced(const coord& taus, space_T& T1, space_T& T2, double omega, int hertz=0);

	// Input		taus : correlation times (sec)
	//			T1   : spatial tensor mu1
	//			T2   : spatial tensor mu2
	//	 		omega: frequency, (Hz or rad/sec)
	//			rad  : flag for Hz or rad/sec
	// Output		J(w) : spectral density value at omega
	// 			       for isotropic rigid diffusion
	// Note		             : output units match tau	

// ______________________________________________________________________
//                        Lipari and Szabo Functions
//             Reference : JACS, Vol. 104, #17, 4545-4570, 1982
//             Specifics : Equations (8) - (12) on page 4560
// ______________________________________________________________________
  
// ************************* Overall Isotropic **************************

MSVCDLL double J_LZ_iso(double S, double tauM, double taue, double omega);
  
	// Input		S    : general order parameter, [0,1]
	//	 		tauM : overall correlation time (sec)
	//	 		taue : effective correlation time (sec)
	// Output		J(w) : spectral density value at omega
	// Note		             : output units match tau	


// ************************ Overall Anisotropic *************************

MSVCDLL double J_LZ_aniso(double S, double A, double tau1,
				double tau2, double taue, double omega);
  
	// Input		S    : general order parameter, [0,1]
	// 			A    : general order parameter, [0,1]
	//	 		tau1 : overall correlation time (sec)
	//	 		tau2 : overall correlation time (sec)
	//	 		taue : effective correlation time (sec)
	// Output		J(w) : spectral density value at omega
	// Note		              : output units match tau	



#endif

