/* HSham.cc ****************************************************-*-c++-*-*
**									**
**	                        G A M M A	 			**
**								 	**
**	NMR Hamiltonians                        Implementation 		**
**								 	**
**	Copyright (c) 1991, 1992, 1993		 			**
**	Scott Smith				 			**
**	Eidgenoessische Technische Hochschule	 			**
**	Labor fuer physikalische Chemie		 			**
**	8092 Zurich / Switzerland		 			**
**								 	**
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
** 	Description						 	**
**								 	**
**  This module of the GAMMA NMR Library provides functions		**
**  designed to facilitate construction of the more commonly 		**
**  encountered Hamiltonians in magnetic resonance.			**
**						 			**
*************************************************************************/
// sosi - Write a check that the secular approximation is valid for HJ!
//        perhaps one in which some scalar coupling flags are automagically set!
// sosi - One day add a basis to spin_sys so all of these operators share
//        the exact same basis. Zum beispiel, then one can use on return 
//  Op = gen_op(SOp.matrix(), sys.get_basis());

#ifndef _HSham_cc_			// Is file already included?
#define _HSham_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma implementation			// this is the implementation
#endif

#include <HSLib/HSham.h>		// Include the file headers
#include <string>			// Knowledge of strings 
#include <HSLib/SpinSystem.h>		// Knowldege of spin systems
#include <HSLib/GenOp.h>		// Knowldege of general operators
#include <HSLib/SpinOp.h>		// Knowledge of spin operators
#include <HSLib/SpinOpCmp.h>		// Knowledge of composite spin ops.
#include <HSLib/SpinOpRot.h>		// Knowledge of rotated spin ops.

// ____________________________________________________________________________
// A                 ISOTROPIC CHEMICAL SHIFT HAMILTONIANS
// ____________________________________________________________________________

/* These functions produce the isotropic component of the chemical shift
   Hamiltonian. This can be done either in the laboratory frame or in multiple
   rotating frames (one frame per isotope type.)  Here we take the static
   magnetic field to be pointing in the +z direction in the lab frame.
 
     Function                             Result
     --------  --------------------------------------------------------------
       Hcs     Isotropic shift Hamiltonian in m-rotating frame where m is the
               the number of unique isotopes in the spin system sys.
     Hcs_lab   Isotropic shift Hamiltonian in the laboratory frame.
 
   The resulting operators will be output in the product basis.  They will
   be diagonal.  System shifts are input relative to the reference frespec
   (or Omega for the isotope).  Spins with a positive gyromagnetic ratio will
   be deshielded from whatever resonates at frespec if they have a positve
   shift. The - sign in the formula below stems from E = -<u|B> where the
   field vector B points along positive z in GAMMA.  In the rotating frame
   the shifts (and H elements) will be on a kHz scale whereas those in the 
   laboratory frame will be in MHz with common NMR field strengths. 
 
                ---                                    ---
          H   = \   - v  * I    <=======>    H       = \   - V  * I 
           cs   /      i    z,i               cs_lab   /      i    z,i
                ---                                    ---                   */

gen_op Hcs(const spin_system& sys)
  {
  spin_op SOp;
  int ns = sys.spins();
  for(int i=0; i<ns; i++)
    {
    if(sys.shift(i) != 0.0 || !i)
      SOp -= sys.shift(i) * Iz(sys,i);
    }
  gen_op H(SOp);
  H.name("Shift Hamiltonian");
  return H;
  } 

gen_op Hcs_lab(const spin_system& sys)
  {
  int ns = sys.spins();
  spin_op SOp;
  for(int i=0; i<ns; i++)
    if(sys.lab_shift(i) != 0.0 || !i)
      SOp -= sys.lab_shift(i) * Iz(sys,i);
//  Op = gen_op(SOp.matrix(), sys.get_basis());
//  sosi: Can't use above without basis in spin sys
//  Op = gen_op(SOp.matrix());
  gen_op H(SOp);
  H.name("Shift Hamiltonian, Lab. Frame");
  return H;
  } 


// ____________________________________________________________________________
// B                  ISOTROPIC SCALAR COUPLING HAMILTONIANS
// ____________________________________________________________________________

/* These functions produce the isotropic component of the scalar coupling
   Hamiltonian. They come in a couple of flavors, depending upon the level of
   approximation you choose to make.

     Function                             Result
     --------  --------------------------------------------------------------
       HJw     Hamiltonian with all weak coupling, no matter what nuclei.
        HJ     Hamiltonian with all strong coupling, no matter what nuclei.
       HJwh    Hamiltonian with weak hetero-, strong homo-nuclear coupling

   The resulting operators will be output in the product basis.  They will be
   either diagonal (all weak coupling) or Hermitian. Note that it is the
   last one (HJwh) that is most often used in NMR. HJ is used for computations
   in the laboratory frame and HJw used when you want to assume a very
   strong field (and have your spectra look spiffy.)                         */

gen_op HJw(const spin_system& sys)
  {
  spin_op SOp;
  int ns = sys.spins();
  for(int i=0; i<ns-1; i++)
    for(int j=i+1; j<ns; j++)
      if(!sys.nepair(i,j) && fabs(sys.J(i,j)) > 1.e-5)
	SOp += sys.J(i,j) * Iz(sys,i) * Iz(sys,j);
  gen_op H(SOp);
  H.name("Weak Scalar Coupling Hamiltonian");
  return H;
  }

gen_op HJ(const spin_system& sys)
  {
  spin_op SOp;
  spin_op SOp1;
  int ns = sys.spins();
  for(int i=0; i<ns-1; i++)
    for(int j=i+1; j<ns; j++)
      if(!sys.nepair(i,j) && fabs(sys.J(i,j)) > 1.e-5)
	{
	SOp1  = Iz(sys,i) * Iz(sys,j);
	SOp1 += Iy(sys,i) * Iy(sys,j);
	SOp1 += Ix(sys,i) * Ix(sys,j);
	SOp  += sys.J(i,j) * SOp1;
	}
  gen_op H(SOp);
  H.name("Weak Scalar Coupling Hamiltonian");
  return H;
  }

gen_op HJwh(const spin_system& sys)
  {
  spin_op SOp;
  spin_op SOp1;
  int ns = sys.spins();
  for(int i=0; i<ns-1; i++)
    {
    for(int j=i+1; j<ns; j++)
      {
      if(!sys.nepair(i,j) && fabs(sys.J(i,j)) > 1.e-5)
	{
	SOp1 = Iz(sys,i) * Iz(sys,j);		// Secular part of F(i)*F(j)
	if(sys.isotope(i) == sys.isotope(j))
	  {
	  SOp1 += Iy(sys,i) * Iy(sys,j);	// Non-secular part of F(i)*F(j)
	  SOp1 += Ix(sys,i) * Ix(sys,j);
	  }
	SOp1 *= sys.J(i,j);
	SOp += SOp1;
	}
      }
    }
  gen_op H(SOp);
  H.name("Weak Heteronuclear Scalar Coupling Hamiltonian");
  return H;
  }



	// Input	       sys  : Spin system
	// 		       iso  : Isotope type
	// Output	        HJ  : Isotropic scalar coupling
	//			      Hamiltonian (default basis)
	// Note			      This gives a "decoupled" Hamiltonian
	//			        . Heteronuclear J's are weak
	//			        . J's involving iso are set 0
	// 

gen_op HJd(const spin_system& sys, const std::string& iso)
  {
  spin_op SOp;
  spin_op SOp1;
  int ns = sys.spins();
  for(int i=0; i<ns-1; i++)
    {
    for(int j=i+1; j<ns; j++)
      if(fabs(sys.J(i,j)) > 0.0001
           && sys.symbol(i) != iso  && sys.symbol(j) != iso)
	{
	SOp1 = Iz(sys,i) * Iz(sys,j);		// Secular part of F(i)*F(j)
	if(sys.isotope(i) == sys.isotope(j))
	  {
	  SOp1 += Iy(sys,i) * Iy(sys,j);	// Non-secular part of F(i)*F(j)
	  SOp1 += Ix(sys,i) * Ix(sys,j);
	  }
	SOp1 *= sys.J(i,j);
	SOp += SOp1;
	}
    }
//  Op = gen_op(SOp.matrix(), sys.get_basis());
//  Op = gen_op(SOp.matrix());
  return gen_op(SOp);
  }


// ____________________________________________________________________________
// C                          ISOTROPIC HAMILTONIANS
// ____________________________________________________________________________

/* These functions provide the isotropic high-resolution NMR Hamiltonian. This
   will include chemical shift and scalar coupling terms, but nothing else.
   Variations allow for weak coupling, weak hetero-nuclear coupling, and 
   the Hamiltonian in the laboratory frame.

   Function                          Returned Operator
   --------  ------------------------------------------------------------------
      Ho     Commonly used liquid NMR Hamiltonian, weak heteronuclear coupling
      How    Same as above but with forced weak scalar couplings for all
     Ho_lab  Isotropic Hamiltonian but referenced to the laboratory frame

   Note that both Ho and How are generated in one or more rotating frames. That
   is, they assume that there are m independent rotating frames where m is the
   number of unique isotopes in the spin system.  If the system has spins that
   resonate close together (as will occur with some isotope types, especially
   in low field simulations), these Hamiltonians will not be appropriate. In
   such case put all shifts into the same rotating frame or simply work in 
   laboratory frame. Also, remember you can always build any spin Hamiltonian
   you want by adding up spin operators multiplied by system parameters      */

gen_op Ho(const spin_system& sys)
  { gen_op H = Hcs(sys);     H += HJwh(sys); H.name("Isotropic Hamiltonian");              return H; }

gen_op How(const spin_system& sys)
  { gen_op H = Hcs(sys);     H += HJw(sys);  H.name("Isotopic Weak Coupling Hamiltonian"); return H; }

gen_op Ho_lab(const spin_system& sys)
  { gen_op H = Hcs_lab(sys); H += HJ(sys);   H.name("Lab. Frame Isotropic Hamiltonian");   return H; }


// ____________________________________________________________________________
// D                          ZEEMAN HAMILTONIAN
// ____________________________________________________________________________

/* These functions provide the Zeeman Hamiltonian in in the laboratory frame.
   Overloaded values allow one to obtain the Hamiltonian that applies to only
   one spin.
              ---                     ---                  ---
        H   = \   - gamma * B * I   = \   - Omega  * I   = \   H
         Z    /          i   0   zi   /          i    zi   /    Z,i
              ---                     ---                  ---
               i                        i                   i                */

gen_op Hz(const spin_system& sys)
  {
  spin_op SOp;
  int ns = sys.spins();
  for(int i=0; i<ns; i++)
    SOp -= sys.Omega(i) * 1.e+6 * Iz(sys,i);
//  Op = gen_op(SOp.matrix(), sys.get_basis());
//  Op = gen_op(SOp.matrix());
  gen_op H(SOp);
  H.name("Zeeman Hamiltonian");
  return H;
  } 

gen_op Hz(const spin_system& sys, const std::string& I)
  {
  spin_op SOp;
//  double x = -1.e+6 * sys.Omega(I);
//  SOp = x*Iz(sys,I);
  SOp += (-1.0*sys.Omega(I)*1.e+6) * Fz(sys,I);
//  Op = gen_op(SOp.matrix(), sys.get_basis());
//  Op = gen_op(SOp.matrix());
  return gen_op(SOp);
  } 


// ____________________________________________________________________________
// E                     RF-FIELD RELATED HAMILTONIANS
// ____________________________________________________________________________


gen_op H1(const spin_system& sys, const std::string& iso, double gamB1, double phi)

	// Input		sys    :	Spin system
	// 			gamB1 : Applied rf-field strength (Hertz!)
	// 			iso   :	Isotype type to which field is applied
	// Output		H1    :	Isotropic liquid Hamiltonian (Hertz)
	// Note			      :	H1 is output in the frame rotating about
	//			     	the z-axis at the rf-field frequency
	// Note			      : H1 acts only on spins of type iso
	// Note			      : gamB1 sets the field strength as to that
	//				used for a 10 usec 1H pulse of 90 deg.
	//				2.5e4(cyc/sec)*1e-6(sec)=.25cyc=90 deg.
	//
	//		H  = B  * gamma   *  F(phi,{i})
	//		 1    1        {i}    xy

  {
  spin_op SOp;
  SOp = gamB1*Fxy(sys, iso, phi);	// Determine H1 in Hertz 
//  Op = gen_op(SOp.matrix(), sys.get_basis());
//  Op = gen_op(SOp.matrix());
  return gen_op(SOp);
  } 


gen_op Heff(spin_sys& sys, gen_op &H0, const std::string& iso, double Wrf, double gamB1, double phi)

	// Input		sys  : Spin system
	// 		      	H0   : Current static Hamiltonian (Hertz)
	//		      	iso  : Isotope type being affected
	//		      	Wrf  : Frequency of the applied rf-field (Hz)
	// 			gamB1: Applied rf-field strength (Hertz)
	// 			phi  : Applied rf-field phase angle (degrees)
	// Output		Heff : Effective Hamiltonian (Hertz), defined in
	//     			       the rotating frame about the z-axis at Wrf
	// Note			     : H1 acts only on spins of isotope type iso
	// Note			2.5e4: Strength for 10 usec 1H pulse of 90 degrees
	// Note			     : Careful with the sign of Wrf! For a shift of
	//			       100 Hz need Wrf of -100 Hz for resonance!
	//
	//          H    = H  - gamma   * F          + w  * F
	//	     eff    0        {i}   xy,phi,{i}   rf   z,{i}

  {
  spin_op SOp;				// A working spin operator
  gen_op Op(H0);
//  basis bs = sys.get_basis();		// The default system product basis
  if(Wrf)
    {
    SOp = Wrf*Fz(sys, iso);		// Term for shift into the Wrf rot. frame
    Op += gen_op(SOp);			// Op is now Ho - Wrf*Fz: rot. frame shift
//    Op += gen_op(SOp.matrix());		// Op is now Ho - Wrf*Fz: rot. frame shift
//    Op += gen_op(SOp.matrix(), bs);	// Op is now Ho - Wrf*Fz: rot. frame shift
    }
  if(gamB1)
    {
    SOp = gamB1*Fxy(sys, iso, phi);	// Field Hamiltonian in rot. frame (negated)
    Op -= gen_op(SOp);			// Now add in the field component
//    Op -= gen_op(SOp.matrix());		// Now add in the field component
//    Op -= gen_op(SOp.matrix(), bs);	// Now add in the field component
    }
  return Op;
  } 


// ____________________________________________________________________________
// F                ISOTROPIC ELECTRON G FACTOR HAMILTONIANS
// ____________________________________________________________________________

/* These functions produce the isotropic component of the electron G tensor
   Hamiltonian. This can be done either in the laboratory frame or in the e-
   rotating frame (at the spectrometer base frequency which is assumed set to
   that used in EPR.)  Here we take the static magnetic field to be pointing in
   the +z direction in the lab frame.
 
     Function                             Result
     --------  --------------------------------------------------------------
       Hcs     Isotropic shift Hamiltonian in m-rotating frame where m is the
               the number of unique isotopes in the spin system sys.
     Hcs_lab   Isotropic shift Hamiltonian in the laboratory frame.

   The resulting operators will be output in the product basis.  They will
   be diagonal.  The rotating frame Hamiltonian is produced in that valid
   for the "free electron". G tensors are input unitless and electron 
   resonance frequencies are field dependent. This Hamiltonian has a negative 
   sign for consistency with the nuclear spin chemical shift Hamiltonian,
   likely this is NOT the sign convention ESR people like.  Can't please
   everybody though. In GAMMA we MUST INSIST THAT Hg_lab == Hg + Hz! All 
   electrons are "shielded" from the free electron (negative gamma, positive G).

                e-                                      e-
                ---                                     ---
          H   = \   - v   * S    <=======>    H       = \   - V    * I 
           cs   /      e,i   z,i               cs_lab   /      e,i    z,i
                ---                                     --- 
                 i                                       i                   */

gen_op Hg(const spin_system& sys)

  {
  if(!sys.Omega())
    std::cout << "\n\tWarning: Electron G Hamiltonian Request With No Field Set!";
  spin_op SOp; 					// Empty spin operator
  int ns = sys.spins();                         // Number of spins
  for(int i=0; i<ns; i++)                       // Loop spins
    if(sys.electron(i)) 			// Only sum over electrons
      SOp -= sys.eshift(i) * Iz(sys,i);
  gen_op H(SOp);
  H.name("Electron G Hamiltonian");
  return H;
  }

gen_op Hg_lab(const spin_system& sys)
  {                                 
  if(!sys.Omega())
    std::cout << "\n\tWarning: Electron G Hamiltonian Request With No Field Set!";
  spin_op SOp;  				// Empty spin operator
  int ns = sys.spins();                         // Number of spins
  for(int i=0; i<ns; i++)                       // Loop spins
    if(sys.electron(i)) 			// Only sum over electrons
      SOp -= sys.lab_eshift(i) * Iz(sys,i);
  gen_op H(SOp);
  H.name("Lab. Frame Electron G Hamiltonian");
  return H;
  }

// ____________________________________________________________________________
// G            ELECTRON HYPERFINE COUPLING RELATED HAMILTONIANS
// ____________________________________________________________________________
                                                                                
gen_op HAw(const spin_system& sys)
 
        // Input               sys  : Spin system
        // Output               HJw : Isotropic "weak" hyperfine coupling
        //                            Hamiltonian (default basis)
        // Note                       In this routine all A's are weak
 
  {
  spin_op SOp;					// Empty spin operators
  int ns = sys.spins();                         // Number of spins
  for(int i=0; i<ns-1; i++)                     // Loop spin pairs
    for(int j=i+1; j<ns; j++)
      if(sys.enpair(i,j) && fabs(sys.A(i,j)) > 1.e-4)
        SOp += sys.AHz(i,j)*Iz(sys,i)* Iz(sys,j);
  gen_op H(SOp);
  H.name("Hyperfile Coupling Hamiltonian");
  return H;
  }

// ____________________________________________________________________________
// H                   QUADRUPOLAR RELATED HAMILTONIANS
// ____________________________________________________________________________


gen_op HQsec(const spin_system& sys, double wQ, int i)

        // Input                sys  : Spin system
        //                      Q    : Quadrupolar coupling (Hz)
	//			i    : Spin index
        // Output               HQsec: The secular part of the quadrupolar
        //                             Hamiltonian (default basis, Hz)
	//			       for the spin i

//                                        2    2
//                          HQ  = wQ*[3*Iz  - I ]

  {
  spin_op SOp; 				// Start with empty spin operator
  double Ival = sys.qn(i); 		// Get the spin quantum number
  if(Ival<=0.5)				// Nothing if I = 1/2  
    return gen_op(SOp);			// Create a null operator
  SOp += 2.0*Iz(sys,i)*Iz(sys,i);	//      Add term 3IzIz
  SOp -= Ix(sys,i)*Ix(sys,i);		//      Subtract term IxIx
  SOp -= Iy(sys,i)*Iy(sys,i);		//      Subtract term IyIy
  return gen_op(wQ*SOp);		// Switch to an operator
  }

#endif 							// HSham.cc
