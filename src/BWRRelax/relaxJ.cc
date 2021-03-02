/* relaxJ.cc ************************************-*-c++-*-
**							**
**	              G A M M A				**
**						 	**
**	Spectral Density Functions   Implementation	**
**						 	**
**	Copyright (c) 1993			 	**
**	Scott Smith				 	**
**	Eidgenoessische Technische Hochschule	 	**
**	Labor fuer physikalische Chemie		 	**
**	8092 Zurich / Switzerland		 	**
**						 	**
**      $Header: $
**						 	**
** Modifications:					**
**						 	**
*********************************************************/

/*********************************************************
**						 	**
** 	Description				 	**
**						 	**
**	The GAMMA NMR Library Provides Functions	**
**	for the Simulation of Magnetic Resonance	**
**	Experiments and Associated Mathematical		**
**	Capabilities.  These are a collection of	**
**	routines whichcompute spectral density		**
**	functions or quantities related to spectral	**
**	density	functions.				**
**						 	**
*********************************************************/

#ifndef _relax_J_cc_		// Is file already included?
#define _relax_J_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation		// This is the implementation
#endif

#include <BWRRelax/relaxJ.h>
#include <stdlib.h>


// ______________________________________________________________________
//                Spectral Density Function Error Handling
// ______________________________________________________________________


void J_error (int i)

	// Input		int  : Error Flag
	// Output		none : Error Message Output

{
  std::cout << "\nSpectral Density Functions: ";
  switch (i)
    {
    case 0:
      std::cout << "\nUnknown error.";
      break;
    default:
      std::cout << "Unknown error.\n";
      break;
    }
}


void volatile J_fatality (int error)

	// Input		none :
	// Output		none : Stops Execution & Error Message Output

{
  J_error (error);
  std::cout << "\nSpectral Density Functions: Program Aborting.\n";
  exit(-1);
}


// ______________________________________________________________________
//                      General Spectral Density Function
// ______________________________________________________________________


  double J_gen(double tau, double omega, int hertz)

	// Input		tau  : correlation time (sec)
	//	 		omega: frequency, (Hz or rad/sec)
	//			rad  : flag for Hz or rad/sec
	// Output		J(w) : spectral density value at omega
	// 			       for isotropic rigid diffusion
	// Note		             : output units match tau	

  {
  double J;
  if(hertz) omega *= (6.283185307);	// Switch to radians per second if needed
  J = tau/(1+tau*tau*omega*omega);
  return J;
  }


  matrix J_gen(double tau, double* w, int hs, int hertz)


	// Input		tau  : correlation time (sec)
	//	 		w    : vector of frequenies, (Hz or 1/sec)
	//			hertz: flag for Hz or rad/sec
	// Output		J    : matrix of spectral density values
	//			       at the transition frequencies between
	// 			       the values contained in the vector w
	// Note		             : output units match tau	

  {
  matrix J(hs, hs, complex0); 
  double Jij=0;
  double omega=0;
  for(int i=0; i<hs; i++)
    for(int j=0; j<hs; j++)
      {
      omega = w[i]-w[j];		// Compute transition frequency
      if(hertz) omega *= (6.283185307);	// Switch to rad/sec if needed
      Jij=tau/(1+tau*tau*omega*omega);	// Calculate J value
      J.put(Jij,i,j);
      }
  return J;
  }


  matrix J_gen_shft(double tau, double* w, double shift, int hs, int hertz)


	// Input		tau  : correlation time (sec)
	//	 		w    : vector of frequenies, (Hz or 1/sec)
	//			shift: transition frequency shift
	//			hs   : Hilbert space size
	//			hertz: flag for Hz or rad/sec
	// Output		J    : matrix of spectral density values
	//			       at the transition frequencies between
	// 			       the values contained in the vector w
	// Note		             : output units match tau	

  {
  matrix J(hs, hs, complex0);
  double Jij=0;
  double omega=0;
  for(int i=0; i<hs; i++)
    for(int j=0; j<hs; j++)
      {
      omega = w[i] - w[j] + shift;	// Compute shifted transition frequency
      if(hertz) omega *= (6.283185307);	// Switch to rad/sec if needed
      Jij=tau/(1+tau*tau*omega*omega);	// Calculate J value
      J.put(Jij,i,j);
      }
  return J;
  }


// ______________________________________________________________________
//                  Anisotropic Rotor Auxiliary Functions 
// ______________________________________________________________________
//

  void tausD(double* tau, double Dx, double Dy, double Dz)

	// Input		tau  : double vector for tau values
	// 			Dx   : x-axis rot. diffusion const. (rad/sec)
	// 			Dy   : y-axis rot. diffusion const. (rad/sec)
	// 			Dz   : z-axis rot. diffusion const. (rad/sec)
	// Output		tau  : five effective correlation times (sec)
	// 			       for isotropic rigid diffusion

  {
  double Dp = 0.5*(Dx+Dy);
  double Dm = 0.5*(Dx-Dy);
  double Dtmp1 = Dz - Dp;
  double Dtmp2 = 2.0*sqrt((Dtmp1*Dtmp1) + 3.0*Dm*Dm);
  tau[0] = 1.0/(2.0*Dp + 4.0*Dz);		// tau -2 case
  tau[1] = 1.0/(5.0*Dp - 3.0*Dm + Dz);		// tau -1 case
  tau[2] = 1.0/(4.0*Dp + 2.0*Dz - Dtmp2);	// tau  0 case
  tau[3] = 1.0/(5.0*Dp + 3.0*Dm + Dz);		// tau +1 case
  tau[4] = 1.0/(4.0*Dp + 2.0*Dz + Dtmp2);	// tau +2 case
  return;
  }


  void tausD(double* tau, const coord& Ds)

	// Input		tau  : double vector for tau values
	// 			Ds   : rotational diffusion consts. (rad/sec)
	// Output		tau  : five effective correlation times (sec)
	// 			       for isotropic rigid diffusion

  {
  double Dp = 0.5*(Ds.x()+Ds.y());
  double Dm = 0.5*(Ds.x()-Ds.y());
  double Dtmp1 = Ds.z() - Dp;
  double Dtmp2 = 2.0*sqrt((Dtmp1*Dtmp1) + 3.0*Dm*Dm);
  tau[0] = 1.0/(2.0*Dp + 4.0*Ds.z());		// tau -2 case
  tau[1] = 1.0/(5.0*Dp - 3.0*Dm + Ds.z());	// tau -1 case
  tau[2] = 1.0/(4.0*Dp + 2.0*Ds.z() - Dtmp2);	// tau  0 case
  tau[3] = 1.0/(5.0*Dp + 3.0*Dm + Ds.z());	// tau +1 case
  tau[4] = 1.0/(4.0*Dp + 2.0*Ds.z() + Dtmp2);	// tau +2 case
  return;
  }


  double chiD(double Dx, double Dy, double Dz)

	// Input		Dx   : x-axis rot. diffusion const. (rad/sec)
	// 			Dy   : y-axis rot. diffusion const. (rad/sec)
	// 			Dz   : z-axis rot. diffusion const. (rad/sec)
	// Output		chi  : chi value for reduced spectral density
	// 			       function coefficient computations

  {
  double chi = 1.5707963;			// default is pi/2
  double Dp = 0.5*(Dx+Dy);
  double Dm = 0.5*(Dx-Dy);
  if(Dz != Dp)
    {
//    double arg = sqrt(3.0)*Dm/(Dz - Dm);
//    chi = 0.5 * atan(arg);
// sosi - the above seems to be giving incorrect answers.
//        DEFINATELY, the denominator needs Dp not Dm!
//        Also, the factor of 1/2 seems unnecessary and
//	  the sign is incorrect.  These latter two are not
//	  at the moment worked out in explicit detail.
//        I changed it but be wary!
//	  The next two lines are the patch 5/31/94
    double arg = sqrt(3.0)*Dm/(Dz - Dp);
    chi = -atan(arg);
    }
  return chi;
  }


  double chiD(const coord& Ds)

	// Input		Ds   : rotational diffusion consts. (rad/sec)
	// Output		chi  : chi value for reduced spectral density
	// 			       function coefficient computations

  { return chiD(Ds.x(), Ds.y(), Ds.z()); }


  void taust(double* tau, double taux, double tauy, double tauz)

	// Input		tau  : double vector for tau values
	// 			taux : x-axis correlation time (sec)
	// 			tauy : y-ayis correlation time (sec)
	// 			tauz : z-azis correlation time (sec)
	// Output		tau  : five effective correlation times (sec)
	// 			       for isotropic rigid diffusion

  {
  double Dx = 1.0/(6.0*taux);		// First compute the three rotational
  double Dy = 1.0/(6.0*tauy);		// diffusion constants
  double Dz = 1.0/(6.0*tauz);
  tausD(tau, Dx, Dy, Dz);
  return;
  }


  void taust(double* tau, const coord& taus)

	// Input		tau  : double vector for tau values
	// 			taus : correlation times (sec)
	// Output		tau  : five effective correlation times (sec)
	// 			       for isotropic rigid diffusion

  {
  double Dx = 1.0/(6.0*taus.x());	// First compute the three rotational
  double Dy = 1.0/(6.0*taus.y());	// diffusion constants
  double Dz = 1.0/(6.0*taus.z());
  tausD(tau, Dx, Dy, Dz);
  return;
  }


  double chit(double taux, double tauy, double tauz)

	// Input		taux : x-axis correlation time (sec)
	// 			tauy : y-ayis correlation time (sec)
	// 			tauz : z-azis correlation time (sec)
	// Output		chi  : chi value for reduced spectral density
	// 			       function coefficient computations

  {
  double Dx = 1.0/(6.0*taux);		// First compute the three rotational
  double Dy = 1.0/(6.0*tauy);		// diffusion constants
  double Dz = 1.0/(6.0*tauz);
  return chiD(Dx, Dy, Dz);
  }


  double chit(const coord& taus)

	// Input		taus : correlation times (sec)
	// Output		chi  : chi value for reduced spectral density
	// 			       function coefficient computations

  {
  double Dx = 1.0/(6.0*taus.x());	// First compute the three rotational
  double Dy = 1.0/(6.0*taus.y());	// diffusion constants
  double Dz = 1.0/(6.0*taus.z());
  return chiD(Dx, Dy, Dz);
  }


  void Jcoeffs(double* c, double alpha, double beta, double gamma,
						 double chi, double eta)

	// Input		c    : double vector for 5 coefficients
	// 			alpha: Euler angle (radians)
	// 			beta : Euler angle (radians)
	// 			gamma: Euler angle (radians)
	//			chi  : Rotational correlation time ratio (radians)
	// 			eta  : Tensor value (unitless)
	// Output		c    : five spectral density function coefficients
	// 			       for isotropic rigid diffusion

  {
  double sqrt3 = sqrt(3.0);
  const double K0 = 0.546274215;		// sqrt[15/(16*pi)]
  const double K1 = 1.092548431;		// sqrt[15/(4.0*pi)]
  const double K2 = 0.315391565;		// sqrt[5/(16.0*pi)]

  double S2g = sin(2.0*gamma);
  double C2g = cos(2.0*gamma);
  double S2a = sin(2.0*alpha);
  double C2a = cos(2.0*alpha);

  double l = sin(beta)*cos(alpha);
  double m = sin(beta)*sin(alpha);
  double n = cos(beta);

  double etapart = (1.0+n*n)*C2a*C2g - 2.0*S2a*S2g*n;
  double X = (l*l - m*m) + (eta/3.0)*etapart;
  double Y = (3.0*n*n -1.0) + (eta*(1.0-n*n)*C2g);

  etapart = (eta/3.0)*((1+n*n)*S2a*S2g + 2.0*n*C2a*S2g);
  c[0] = K0 * (2.0*l*m + etapart);			// c -2 case

  etapart = (eta/3.0)*(m*S2g - l*n*C2g);
  c[1] = K1 * (l*n + etapart);				// c -1 case

  c[2] = K2 * (sqrt3*X*sin(chi) - Y*cos(chi));		// c  0 case

  etapart = (eta/3.0)*(l*S2g + m*n*C2g);
  c[3] = K1 * (m*n - etapart);				// c +1 case

  c[4] = K2 * (sqrt3*X*cos(chi) + Y*sin(chi));		// c +2 case

  return;
  }


  void Jcoeffs(double* c, const coord& EA, double chi, double eta)

	// Input		c    : double vector for 5 coefficients
	// 			EA   : Euler angles (radians)
	//			chi  : Ratio of rotational correlation times
	// 			eta  : Tensor value
	// Output		c    : five spectral density function coefficients
	// 			       for isotropic rigid diffusion

  {
  Jcoeffs(c, EA.x(), EA.y(), EA.z(), chi, eta);
  return;
  }


// ______________________________________________________________________
//             Anisotropic Rotor Spectral Density Functions 
// ______________________________________________________________________
//

  double J_reduced(double *tau, double *c1, double *c2, double omega, int hertz)

	// Input		tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//	 		omega: frequency, (Hz or rad/sec)
	//			rad  : flag for Hz or rad/sec
	// Output		J(w) : spectral density value at omega
	// 			       for isotropic rigid diffusion
	// Note		             : output units match tau	

  {
  double J, Jtmp;
  if(hertz) omega *= (6.283185307);	// Switch to radians per second if needed
  J = 0;
  for(int i=0; i<5; i++)
    {
    Jtmp = tau[i]/(1+tau[i]*tau[i]*omega*omega);
    J += c1[i]*c2[i]*Jtmp;
    }
  return J/5.0;
  }


  matrix J_reduced(gen_op& Op, double *tau, double *c1, double *c2, int hertz)


	// Input		Op   : general operator
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transition
	// 			       frequencies of Op
	// Note		             : output units match tau	

  {
  matrix Jmx(Op.dim(), Op.dim(), complex0);
  Op.set_EBR();
  int nl = Op.dim(); 
  double w = 0;
  for(int i=0; i<nl; i++)
    for(int j=0; j<nl; j++)
      {
      w = Re(Op.get(i,i)-Op.get(j,j));
      Jmx.put(J_reduced(tau,c1,c2,w,hertz), i, j);
      }
  return Jmx;
  }


  matrix J_reduced(double* w, int size, double *tau, double *c1, double *c2, int hertz)


	// Input		w    : double vector of frequencies (energies)
	//			size : Hilbert space size
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transitions
	// 			       between the frequencies in w
	// Note		             : output units match tau	

  {
  matrix Jmx(size, size, complex0);
  double wtr = 0;
  for(int i=0; i<size; i++)
    for(int j=0; j<size; j++)
      {
      wtr = w[i] - w[j];
      Jmx.put(J_reduced(tau,c1,c2,wtr,hertz), i, j);
      }
  return Jmx;
  }


  matrix J_red_shft(double* w, double shift, int size, double *tau,
						 double *c1, double *c2, int hertz)


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

  {
  matrix Jmx(size, size, complex0);
  double wtr = 0;
  for(int i=0; i<size; i++)
    for(int j=0; j<size; j++)
      {
      wtr = w[i] - w[j] + shift;
      Jmx.put(J_reduced(tau,c1,c2,wtr,hertz), i, j);
      }
  return Jmx;
  }


  matrix J_reduced(row_vector w, double *tau, double *c1, double *c2, int hertz)


	// Input		w    : row vector of frequencies (energies)
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transitions
	// 			       between the frequencies in w
	// Note		             : output units match tau	

  {
  matrix Jmx(w.size(), w.size(), complex0);
  complex wtr = 0;
  for(int i=0; i<w.size(); i++)
    for(int j=0; j<w.size(); j++)
      {
      wtr = w.get(i)-w.get(j);
      Jmx.put(J_reduced(tau,c1,c2,Re(wtr),hertz), i, j);
      }
  return Jmx;
  }


matrix J_red_shft(row_vector w, double shift, double *tau,
					 double *c1, double *c2, int hertz)

	// Input		w    : row vector of frequencies (energies)
	//			shift: transition frequency shift
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transitions
	// 			       between the frequencies in w
	// Note		             : output units match tau	

  {
  matrix Jmx(w.size(), w.size(), complex0);
  complex wtr = 0;
  for(int i=0; i<w.size(); i++)
    for(int j=0; j<w.size(); j++)
      {
      wtr = w.get(i)-w.get(j)+shift;
      Jmx.put(J_reduced(tau,c1,c2,Re(wtr),hertz), i, j);
      }
  return Jmx;
  }


  double J_reduced(const coord& taus, const coord& EA1, double eta1,
                            const coord& EA2, double eta2, double omega, int hertz)

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

  {
  double tau[5], c1[5], c2[5];
  taust(tau, taus);			// get 5 effective taus
  double chi = chit(taus);		// get chi value
  Jcoeffs(c1, EA1, chi, eta1);		// get 5 coefficients, mu 1
  Jcoeffs(c2, EA2, chi, eta2);		// get 5 coefficients, mu 2
  return J_reduced(tau, c1, c2, omega, hertz);
  }


  double J_reduced(const coord& taus, space_T& T1, space_T& T2, double omega, int hertz)

	// Input		taus : correlation times (sec)
	//			T1   : spatial tensor mu1
	//			T2   : spatial tensor mu2
	//	 		omega: frequency, (Hz or rad/sec)
	//			rad  : flag for Hz or rad/sec
	// Output		J(w) : spectral density value at omega
	// 			       for isotropic rigid diffusion
	// Note		             : output units match tau	

  { return J_reduced(taus, T1.PASys_EA(), T1.eta(), T2.PASys_EA(), T2.eta(), omega, hertz); }
// **** At this point, still need to clearly define tensor angles & angle access


// ______________________________________________________________________
// Anisotropic Rotor Spectral Density Functions - Dynamic Frequency Shift 
// ______________________________________________________________________
//

  double Q_reduced(double *tau, double *c1, double *c2, double omega, int hertz)

	// Input		tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//	 		omega: frequency, (Hz or rad/sec)
	//			rad  : flag for Hz or rad/sec
	// Output		J(w) : spectral density value at omega
	// 			       for isotropic rigid diffusion
	// Note		             : output units match tau	

  {
  double J, Jtmp;
  if(hertz) omega *= (6.283185307);	// Switch to radians per second if needed
  J = 0;
  for(int i=0; i<5; i++)
    {
    Jtmp = tau[i]/(1+tau[i]*tau[i]*omega*omega);
    J += omega*tau[i]*c1[i]*c2[i]*Jtmp;
    }
  return J/5.0;
  }


  matrix Q_reduced(gen_op& Op, double *tau, double *c1, double *c2, int hertz)


	// Input		Op   : general operator
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transition
	// 			       frequencies of Op
	// Note		             : output units match tau	

  {
  matrix Jmx(Op.dim(), Op.dim(), complex0);
  Op.set_EBR();
  int nl = Op.dim(); 
  double w = 0;
  for(int i=0; i<nl; i++)
    for(int j=0; j<nl; j++)
      {
      w = Re(Op.get(i,i)-Op.get(j,j));
      Jmx.put(Q_reduced(tau,c1,c2,w,hertz), i, j);
      }
  return Jmx;
  }


  matrix Q_reduced(double* w, int size, double *tau, double *c1, double *c2, int hertz)


	// Input		w    : double vector of frequencies (energies)
	//			size : Hilbert space size
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transitions
	// 			       between the frequencies in w
	// Note		             : output units match tau	

  {
  matrix Jmx(size, size, complex0);
  double wtr = 0;
  for(int i=0; i<size; i++)
    for(int j=0; j<size; j++)
      {
      wtr = w[i] - w[j];
      Jmx.put(Q_reduced(tau,c1,c2,wtr,hertz), i, j);
      }
  return Jmx;
  }


  matrix Q_red_shft(double* w, double shift, int size, double *tau,
						 double *c1, double *c2, int hertz)


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

  {
  matrix Qmx(size, size, complex0);
  double wtr = 0;
  for(int i=0; i<size; i++)
    for(int j=0; j<size; j++)
      {
      wtr = w[i] - w[j] + shift;
      Qmx.put(Q_reduced(tau,c1,c2,wtr,hertz), i, j);
      }
  return Qmx;
  }


  matrix Q_reduced(row_vector w, double *tau, double *c1, double *c2, int hertz)


	// Input		w    : row vector of frequencies (energies)
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transitions
	// 			       between the frequencies in w
	// Note		             : output units match tau	

  {
  matrix Jmx(w.size(), w.size(), complex0);
  complex wtr = 0;
  for(int i=0; i<w.size(); i++)
    for(int j=0; j<w.size(); j++)
      {
      wtr = w.get(i)-w.get(j);
      Jmx.put(Q_reduced(tau,c1,c2,Re(wtr),hertz), i, j);
      }
  return Jmx;
  }


  matrix Q_red_shft(row_vector w, double shift, double *tau,
					 double *c1, double *c2, int hertz)


	// Input		w    : row vector of frequencies (energies)
	//			shift: transition frequency shift
	// 			tau  : five correlation times (sec)
	// 			c1   : five coefficients of interaction 1
	// 			c2   : five coefficients of interaction 2
	//			hertz: flag for Hz or rad/sec
	// Output		J(w) : spectral density values at all transitions
	// 			       between the frequencies in w
	// Note		             : output units match tau	

  {
  matrix Qmx(w.size(), w.size(), complex0);
  complex wtr = 0;
  for(int i=0; i<w.size(); i++)
    for(int j=0; j<w.size(); j++)
      {
      wtr = w.get(i)-w.get(j)+shift;
      Qmx.put(Q_reduced(tau,c1,c2,Re(wtr),hertz), i, j);
      }
  return Qmx;
  }


double Q_reduced(const coord& taus, const coord& EA1, double eta1,
                         const coord& EA2, double eta2, double omega, int hertz)

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

  {
  double tau[5], c1[5], c2[5];
  taust(tau, taus);			// get 5 effective taus
  double chi = chit(taus);		// get chi value
  Jcoeffs(c1, EA1, chi, eta1);		// get 5 coefficients, mu 1
  Jcoeffs(c2, EA2, chi, eta2);		// get 5 coefficients, mu 2
  return Q_reduced(tau, c1, c2, omega, hertz);
  }


  double Q_reduced(const coord& taus, space_T& T1, space_T& T2, double omega, int hertz)

	// Input		taus : correlation times (sec)
	//			T1   : spatial tensor mu1
	//			T2   : spatial tensor mu2
	//	 		omega: frequency, (Hz or rad/sec)
	//			rad  : flag for Hz or rad/sec
	// Output		J(w) : spectral density value at omega
	// 			       for isotropic rigid diffusion
	// Note		             : output units match tau	

  {
  return Q_reduced(taus, T1.PASys_EA(), T1.eta(), T2.PASys_EA(), T2.eta(), omega, hertz);
  }

// ______________________________________________________________________
//                        Lipari and Szabo Functions
//             Reference : JACS, Vol. 104, #17, 4545-4570, 1982
//             Specifics : Equations (8) - (12) on page 4560
// ______________________________________________________________________
//
// ************************* Overall Isotropic **************************

  double J_LZ_iso(double S, double tauM, double taue, double omega)
  
	// Input		S    : general order parameter, [0,1]
	//	 		tauM : overall correlation time (sec)
	//	 		taue : effective correlation time (sec)
	// Output		J(w) : spectral density value at omega
	// Note		             : output units match tau	

    {
    double J, tau;
    tau = (tauM*taue)/(tauM+taue); 
    J = (2.0/5.0)*(S*S*J_gen(tauM,omega) + (1-S*S)*J_gen(tau, omega));
    return J;
    }


// ************************ Overall Anisotropic *************************

  double J_LZ_aniso(double S, double A, double tau1,
				double tau2, double taue, double omega)
  
	// Input		S    : general order parameter, [0,1]
	// 			A    : general order parameter, [0,1]
	//	 		tau1 : overall correlation time (sec)
	//	 		tau2 : overall correlation time (sec)
	//	 		taue : effective correlation time (sec)
	// Output		J(w) : spectral density value at omega
	// Note		              : output units match tau	

    {
    double J, taup, taupp;
    taup = (tau1*taue)/(tau1+taue); 
    taupp = (tau2*taue)/(tau2+taue); 
    J = A*S*S*J_gen(tau1,omega) + (1-A)*S*S*J_gen(tau2, omega);
    J += A*(1-S*S)*J_gen(taup,omega) + (1-A)*(1-S*S)*J_gen(taupp, omega);
    J *= (2.0/5.0);
    return J;
    }


#endif

