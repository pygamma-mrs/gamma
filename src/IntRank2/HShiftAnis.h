/* HShiftAnis.h *************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      Oriented Shift Anisotropy Hamiltonians 	     Interface		**
**                                                                      **
**      Scott Smith                                                     **
**      Copyright (c) 2001                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
** This module of the GAMMA MR platform provides Hamiltonians for       **
** shift anisotropy interations.  The functions are typically designed  **
** to take a spin system containing shift anisotropy interactions as an **
** input argument. Differing functions allow users to obtain the 	**
** Hamiltonian under various approximation levels, at any orientation,  **
** and for either all system shift anisotropy interations or just that  **
** for a particular spin. Note that the returned Hamiltonians will	**
** reside in the Hilbert space of the spin system (i.e. the composite	**
** spin Hilbert space.)  						**
**                                                                      **
*************************************************************************/

#ifndef   HCSA_h_			// Is file already included?
#  define HCSA_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <IntRank2/SolidSys.h>		// Knowledge of spin systems
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

MSVCDLL gen_op Hcs(const     spin_system &sys);
MSVCDLL gen_op Hcs_lab(const spin_system &sys);

//matrix IntCSA::H( ) const;
//matrix IntCSA::H(double theta, double phi) const;
//matrix IntCSA::H(const vector<int>& HSs, int i) const;
//matrix IntCSA::H(const vector<int>& HSs, int i, double T, double P) const;

#endif 							//  HShiftAnis.h

