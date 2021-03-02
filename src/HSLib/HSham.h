/* HSham.h ******************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Hilbert Space MR Hamiltonians 		     Implementation	**
**                                                                      **
**      Copyright (c) 1991,1998                                         **
**      Scott Smith                                                     **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zurich / Switzerland                                       **
**                                                                      **
**      $Header: $
**                                                                      **
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The GAMMA MR Library Provides Functions for the Simulation of 	**
** Magnetic Resonance Experiments and Associated Mathematical		**
** Capabilities.  These functions supply many of the simple 		**
** Hamiltonians	commonly used in NMR.					**
**						 			**
*************************************************************************/

#ifndef   HSham_h_			// Is file already included?
#  define HSham_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <HSLib/SpinSystem.h>		// Knowledge of spin systems
#include <HSLib/GenOp.h>		// Knowledge of operators

// ____________________________________________________________________________
// A                  ISOTROPIC CHEMICAL SHIFT HAMILTONIANS
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

MSVCDLL gen_op Hcs(const spin_system &sys);
MSVCDLL gen_op Hcs_lab(const spin_system &sys);

// ____________________________________________________________________________
// B               ISOTROPIC SCALAR COUPLING HAMILTONIANS
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

MSVCDLL gen_op HJw(const  spin_system &sys);
MSVCDLL gen_op HJ(const   spin_system &sys);
MSVCDLL gen_op HJwh(const spin_system &sys);
MSVCDLL gen_op HJd(const spin_system& sys, const std::string& iso);

	// Input	       sys  : Spin system
	// 		       iso  : Isotope type
	// Output	        HJ  : Isotropic scalar coupling
	//			      Hamiltonian (default basis)
	// Note			      This gives a "decoupled" Hamiltonian
	//			        . Heteronuclear J's are weak
	//			        . J's involving iso are set 0
	// 


// ____________________________________________________________________________
// C                 ISOTROPIC HIGH RESOLUTION NMR HAMILTONIANS
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

MSVCDLL gen_op Ho(const     spin_system &ss);
MSVCDLL gen_op How(const    spin_system &ss);
MSVCDLL gen_op Ho_lab(const spin_system &ss);

// ____________________________________________________________________________
// D                           ZEEMAN HAMILTONIAN
// ____________________________________________________________________________

/* These functions provide the Zeeman Hamiltonian in the laboratory frame.
   Overloaded values allow one to obtain the Hamiltonian that applies to only
   one spin.
              ---                     ---                  ---
        H   = \   - gamma * B * I   = \   - Omega  * I   = \   H
         Z    /          i   0   zi   /          i    zi   /    Z,i
              ---	              ---                  ---
               i                        i                   i                */

MSVCDLL gen_op Hz(const spin_system& sys);
MSVCDLL gen_op Hz(const spin_system& sys, const std::string& I);

// ____________________________________________________________________________
// E                      RF-FIELD RELATED HAMILTONIANS
// ____________________________________________________________________________


MSVCDLL gen_op H1(const spin_system& sys, const std::string& iso,
                                    double gamB1=2.5e4, double phi=0.0);

	// Input		ss    :	Spin system
	// 			gamB1 : Applied rf-field strength (Hertz!)
	// 			iso   :	Isotype type to which field is applied
	// Output		H1    :	Isotropic liquid Hamiltonian (rad/sec)
	// Note			      :	H1 is output in the frame rotating about
	//			     	the z-axis at the rf-field frequency
	// Note			      : H1 acts only on spins of isotope type iso
	// Note			      : gamB1 sets the field strength equal to that
	//				used for a 10 usec proton pulse of 90 degrees
	//				2.5e4(cyc/sec)*1e-6(sec) = .25 cyc = 90 degrees
	//
	//		H  = B  * gamma   *  F(phi,{i})
	//		 1    1        {i}    xy


MSVCDLL gen_op Heff(spin_sys& sys, gen_op &H0, const std::string& iso,
                                  double Wrf=0, double gamB1=2.5e4, double phi=0.0);

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
	//
	//          H    = H  - gamma   * F          + w  * F
	//	     eff    0        {i}   xy,phi,{i}   rf   z,{i}

// ____________________________________________________________________________
// F              ISOTROPIC ELECTRON G FACTOR HAMILTONIANS
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
 
                e-
                ---                                    ---
          H   = \   - v   * S    <=======>    H       = \   - V    * I
           cs   /      e,i   z,i               cs_lab   /      e,i    z,i
                ---                                    ---                   */
 
MSVCDLL gen_op Hg(const     spin_system& sys);
MSVCDLL gen_op Hg_lab(const spin_system& sys);
 
// ____________________________________________________________________________
// G            ELECTRON HYPERFINE COUPLING RELATED HAMILTONIANS
// ____________________________________________________________________________

                                                                                
MSVCDLL gen_op HAw(const spin_system& sys);
 
        // Input               sys  : Spin system
        // Output               HJw : Isotropic "weak" hyperfine coupling
        //                            Hamiltonian (default basis)
        // Note                       In this routine all A's are weak
 
 
 
// ____________________________________________________________________________
// I                   QUADRUPOLAR RELATED HAMILTONIANS
// ____________________________________________________________________________

MSVCDLL gen_op HQsec(const spin_system& sys, double wQ, int i);

        // Input                sys  : Spin system
        //                      Q    : Quadrupolar coupling (Hz)
        //                      i    : Spin index
        // Output               HQsec: The secular part of the quadrupolar
        //                             Hamiltonian (default basis, Hz)
        //                             for the spin i

//                                       2    2
//                         HQ  = wQ*[3*Iz  - I ]

 
#endif 							//  HSham.h

