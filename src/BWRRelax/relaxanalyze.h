/* relax_analyze.h *************************************-*-c++-*-*
**								**
**	                   G A M M A	 			**
**							 	**
**	NMR Relaxation Analysis             Interface 		**
**							 	**
**	Copyright (c) 1993			 		**
**	Scott Smith				 		**
**	University of Utah					**
**	Department of Chemistry					**
**      Salt Lake City, UT, 84112, USA                          **
**							 	**
**      $Header:
**						 		**
** Modifications:						**
**							 	**
******************************************************************/

/*****************************************************************
**							 	**
**  Description						 	**
**							 	**
**  This module of the GAMMA relaxation package provides	**
**  routines designed to facilitate the analysis of relaxation	**
**  in arbitrary spin systems.					**
**						 		**
*****************************************************************/

#ifndef   Relax_analyze_h_		// Is this file already included?
#  define Relax_analyze_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

# include <GamGen.h>			// Know MSVCDLL (__declspec)
# include <LSLib/SuperOp.h>		// Include superoperators
# include <LSLib/sys_dynamic.h>		// Include dynamic spin systems
# include <Level1/SpinT.h>		// Include spin tensors
# include <BWRRelax/relaxDip.h>		// Include Dipolar relaxation
# include <BWRRelax/relaxCSA.h>		// Include CSA relaxation
# include <BWRRelax/relaxQuad.h>		// Include Quadrupolar relaxation
# include <BWRRelax/relaxRand.h>		// Include Random Field relaxation
# include <BWRRelax/relaxDCSA.h>		// Include Dipolar-CSA relaxation
//# include <BWRRelax/relaxDQ.h>		// Include Dipolar-Quadrupolar relaxation
//# include <BWRRelax/relaxCSAQ.h>		// Include CSA-Quadrupolar relaxation

// ______________________________________________________________________
//                            WAVEFUNCTIONS
// ______________________________________________________________________

MSVCDLL void RDDel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
                              int a, int aa, int b, int bb,
                               int DFS=0, int Windex=-1, int Sindex=0);
 
//                                >i     l>k
//                DD         --- --- --- ---        DD
//         <a,aa|R  |b,bb> = \   \   \   \   <a,aa|R    |b,bb>
//                           /   /   /   /          ijkl
//                           --- --- --- ---
//                            i   j   k   l

MSVCDLL void RSSel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
                              int a, int aa, int b, int bb,
                               int DFS=0, int Windex=-1, int Sindex=0);
 
//                SASA         --- ---        SASA
//         <a,aa|R    |b,bb> = \   \   <a,aa|R    |b,bb>
//                             /   /          ij
//                             --- ---
//                              i   j

MSVCDLL void RDSel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
                              int a, int aa, int b, int bb,
                               int DFS=0, int Windex=-1, int Sindex=0);
 
//                                >i
//                DSA         --- --- ---        DSA
//         <a,aa|R   |b,bb> = \   \   \   <a,aa|R   |b,bb>
//                            /   /   /          ijk
//                            --- --- ---
//                             i   j   k
 

MSVCDLL void RSDel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
                              int a, int aa, int b, int bb,
                               int DFS=0, int Windex=-1, int Sindex=0);

//                                    >i
//                SAD         --- --- ---        SAD
//         <a,aa|R  |b,bb> = \   \   \    <a,aa|R   |b,bb>
//                           /   /   /           kij
//                           --- --- ---
//                            k   i   j

 
MSVCDLL void RRRel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
                              int a, int aa, int b, int bb,
                               int DFS=0, int Windex=-1, int Sindex=0);
 
//                RFRF         --- ---        RFRF
//         <a,aa|R    |b,bb> = \   \   <a,aa|R    |b,bb>
//                             /   /          ij
//                             --- ---
//                              i   j

 
MSVCDLL void RQQel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
                              int a, int aa, int b, int bb,
                               int DFS=0, int Windex=-1, int Sindex=0);
 
//                QQ         --- ---        QQ
//         <a,aa|R  |b,bb> = \   \   <a,aa|R  |b,bb>
//                           /   /          ij
//                           --- ---
//                            i   j
 

MSVCDLL void RQSel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
                              int a, int aa, int b, int bb,
                               int DFS=0, int Windex=-1, int Sindex=0);
 
//                QSA         --- ---        QSA
//         <a,aa|R   |b,bb> = \   \   <a,aa|R   |b,bb>
//                            /   /         ij
//                            --- ---
//                             i   j

MSVCDLL void RSQel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
                              int a, int aa, int b, int bb,
                               int DFS=0, int Windex=-1, int Sindex=0);
 
//                SAQ         --- ---        SAQ
//         <a,aa|R   |b,bb> = \   \   <a,aa|R   |b,bb>
//                            /   /          ij
//                            --- ---
//                             i   j


MSVCDLL void RQDel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
                              int a, int aa, int b, int bb,
                               int DFS=0, int Windex=-1, int Sindex=0);
 
//                QD         --- ---        QD
//         <a,aa|R  |b,bb> = \   \   <a,aa|R  |b,bb>
//                           /   /          ij
//                           --- ---
//                            i   j


MSVCDLL void RDQel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
                              int a, int aa, int b, int bb,
                               int DFS=0, int Windex=-1, int Sindex=0);
 
//                DQ         --- ---        DQ
//         <a,aa|R  |b,bb> = \   \   <a,aa|R  |b,bb>
//                           /   /          ij
//                           --- ---
//                            i   j
 
// ______________________________________________________________________
// *********** RELAXATION MATRIX WHOLE ELEMENT STRING FUNCTIONS *********
// ______________________________________________________________________


MSVCDLL void Rijkl_el(const sys_dynamic& sys, gen_op& Ho, int rank,  spin_T* T1, spin_T* T2,
                           std::string& Mlabel, 
   std::string* line1, std::string* line2, std::string* line3, int* signs, int& Jterms, int& Lterms,
                              int a, int aa, int b, int bb,
                               int DFS=0, int Windex=-1, int Sindex=0);

        // Input                sys   : Dynamic spin system
        //                      Ho    : Static Hamiltonian
        //                      rank  : Rank of the interactions
        //                      T1    : Spin tensor, mu1
        //                      T2    : Spin tensor, mu2
        //                      Mlabel: Label for two mechanisms
        //                      line1 : Strings for the upper line of component
        //                      line2 : Strings for the middle line of component
        //                      line3 : Strings for the lower line of component
        //                      signs : Flags for the signs of the components
        //                      a, aa : 1st transition indices
        //                      b, bb : 2nd transition indices
        //                      DFS   : Flag for dynamic frequency shifts
        //                                      > 0  J Terms Only
        //                                        0  Both J and L Terms (Default)
        //                                      < 0  L Terms Only
        //                      Windex: Flag for frequency indexing
        //                                 >1 = Always use spin indices
        //                                  0 = Use spin indices if needed
        //                                 -1 = Never use spin indices (default)
        //                                 -2 = Assume extreme narrowing
        //                                 -3 = Assume extreme broadening
        //                      Sindex: Flag for spin indexing
        //                                  0 = Use only if needed (default)
        //                                  1 = Always use spin indices
        //                                 <1 = Never use spin indices
        // Output               void  :
        // Note                       : This routine assumes that each of
        //                              the two interactions involve 2 spins
        // Note                       : The above implies the rank is 2
 
//                                --- --- --- ---
//                    12          \   \   \   \    12
//             <a,a'|R  |b,b'> =  /   /   /   /   R
//                                --- --- --- ---  ij,kl
//                                 i   j   k   l
 

MSVCDLL void Rij_el(const sys_dynamic& sys, gen_op& Ho, int rank,
                  spin_T* T1, spin_T* T2,
                           std::string& Mlabel, 
   std::string* line1, std::string* line2, std::string* line3, int* signs, int& Jterms, int& Lterms,
                              int a, int aa, int b, int bb,
                               int DFS=0, int Windex=-1, int Sindex=0);
 
        // Input                sys   : Dynamic spin system
        //                      Ho    : Static Hamiltonian
        //                      rank  : Rank of the interactions
        //                      T1    : Spin tensor, mu1
        //                      T2    : Spin tensor, mu2
        //                      Mlabel: Label for two mechanisms
        //                      line1 : Strings for the upper line of component
        //                      line2 : Strings for the middle line of component
        //                      line3 : Strings for the lower line of component
        //                      signs : Flags for the signs of the components
        //                      a, aa : 1st transition indices
        //                      b, bb : 2nd transition indices
        //                      DFS   : Flag for dynamic frequency shifts
        //                                      > 0  J Terms Only
        //                                        0  Both J and L Terms (Default)
        //                                      < 0  L Terms Only
        //                      Windex: Flag for frequency indexing
        //                                 >1 = Always use spin indices
        //                                  0 = Use spin indices if needed
        //                                 -1 = Never use spin indices (default)
        //                                 -2 = Assume extreme narrowing
        //                                 -3 = Assume extreme broadening
        //                      Sindex: Flag for spin indexing
        //                                  0 = Use only if needed (default)
        //                                  1 = Always use spin indices
        //                                 <1 = Never use spin indices
        // Output               void  :
        // Note                       : This routine assumes that each of
        //                              the two interactions involve 1 spin
        //                              e.g. SA-SA, Q-Q, Q-SA, SA-Q, RF-RF
    
//                       Two Single Spin Mechanisms
 
//                                        --- ---
//                            12          \   \    12
//                     <a,a'|R  |b,b'> =  /   /   R
//                                        --- ---  i,j
//                                         i   j
 
MSVCDLL void Rijk_el(const sys_dynamic& sys, gen_op& Ho, int rank,
                  spin_T* T1, spin_T* T2,
                           std::string& Mlabel,  
   std::string* line1, std::string* line2, std::string* line3, int* signs, int& Jterms, int& Lterms,
                              int a, int aa, int b, int bb,
                               int DFS=0, int Windex=-1, int Sindex=0);
 
        // Input                sys   : Dynamic spin system
        //                      Ho    : Static Hamiltonian
        //                      rank  : Rank of the interactions
        //                      T1    : Spin tensor, mu1
        //                      T2    : Spin tensor, mu2
        //                      Mlabel: Label for two mechanisms
        //                      line1 : Strings for the upper line of component
        //                      line2 : Strings for the middle line of component
        //                      line3 : Strings for the lower line of component
        //                      signs : Flags for the signs of the components
        //                      a, aa : 1st transition indices
        //                      b, bb : 2nd transition indices
        //                      DFS   : Flag for dynamic frequency shifts
        //                                      > 0  J Terms Only
        //                                        0  Both J and L Terms (Default)
        //                                      < 0  L Terms Only
        //                      Windex: Flag for frequency indexing
        //                                 >1 = Always use spin indices
        //                                  0 = Use spin indices if needed
        //                                 -1 = Never use spin indices (default)
        //                                 -2 = Assume extreme narrowing
        //                                 -3 = Assume extreme broadening
        //                      Sindex: Flag for spin indexing
        //                                  0 = Use only if needed (default)
        //                                  1 = Always use spin indices
        //                                 <1 = Never use spin indices
        // Output               void  :
        // Note                       : This routine assumes that the first
        //                              interaction involves two spins and
        //                              the second interaction involves 1 spin
        //                              e.g. D-SA, D-Q
 
//               Mechanism 1 = Spin-Pair; Mechanism 2 = Spin
 
//                                        --- --- ---
//                            12          \   \   \    12
//                     <a,a'|R  |b,b'> =  /   /   /   R
//                                        --- --- ---  ij,k
//                                         i   j   k
 

MSVCDLL void Rkij_el(const sys_dynamic& sys, gen_op& Ho, int rank,
                  spin_T* T1, spin_T* T2,
                           std::string& Mlabel, 
   std::string* line1, std::string* line2, std::string* line3, int* signs, int& Jterms, int& Lterms,
                              int a, int aa, int b, int bb,
                               int DFS=0, int Windex=-1, int Sindex=0);
 
        // Input                sys   : Dynamic spin system
        //                      Ho    : Static Hamiltonian
        //                      rank  : Rank of the interactions
        //                      T1    : Spin tensor, mu1
        //                      T2    : Spin tensor, mu2
        //                      Mlabel: Label for two mechanisms
        //                      line1 : Strings for the upper line of component
        //                      line2 : Strings for the middle line of component
        //                      line3 : Strings for the lower line of component
        //                      signs : Flags for the signs of the components
        //                      a, aa : 1st transition indices
        //                      b, bb : 2nd transition indices
        //                      DFS   : Flag for dynamic frequency shifts
        //                                      > 0  J Terms Only
        //                                        0  Both J and L Terms (Default)
        //                                      < 0  L Terms Only
        //                      Windex: Flag for frequency indexing
        //                                 >1 = Always use spin indices
        //                                  0 = Use spin indices if needed
        //                                 -1 = Never use spin indices (default)
        //                                 -2 = Assume extreme narrowing
        //                                 -3 = Assume extreme broadening
        //                      Sindex: Flag for spin indexing
        //                                  0 = Use only if needed (default)
        //                                  1 = Always use spin indices
        //                                 <1 = Never use spin indices
        // Output               void  :
        // Note                       : This routine assumes that the first
        //                              interaction involves one spin and
        //                              the second interaction involves 2 spins
        //                              e.g. SA-D, Q-D
 
//                  Mechanism 1 = Spin; Mechanism 2 = Spin-Pair
 
//                                        --- --- ---
//                            12          \   \   \    12
//                     <a,a'|R  |b,b'> =  /   /   /   R
//                                        --- --- ---  k,ij
//                                         k   i   j
 
// ______________________________________________________________________
// ************ NUMERICAL RELAXATION MATRIX ELEMENT FUNCTIONS ***********
// ______________________________________________________________________
 
// These functions fill numerical arrays which represent the terms which
// contribute to a specific relaxation matrix element from a specific
// pair of interaction types (relaxation mechanisms) and a specific set
// of involved spins (e.g. ij-kl, ij-k, i-j).  Terms from both the sym-
// metric and anti-symmetric (dynamic frequency shift) spectral density
// components are treated.
 
MSVCDLL void Rel_12(int hs, int rank, gen_op* T1s, gen_op* T2s,
          int& Jterms, complex* strsJ, int* trnsJ,
               int& Lterms, complex* strsL, int* trnsL,
        std::string* wlabs, int a, int aa, int b, int bb, double cutoff=1.e-4);

        // Input                hs    : Spin system Hilbert space
        //                      rank  : Rank of the two interactions
        //                      T1s   : Spin tensor components, 1st interaction
        //                      T2s   : Spin tensor components, 2nd interaction
        //                      Jterms: Number of non-zero J(w) terms in summation
        //                      strsJ : Array of interaction strengths for J(w)
        //                      trnsJ : Array of transition indices for J(w)
        //                      Lterms: Number of non-zero L(w) terms in summation
        //                      strsL : Array of interaction strengths for L(w)
        //                      trnsL : Array of transition indices for L(w)
        //                      Wlbls : String array of transition labels
        //                      a, aa : 1st transition indices
        //                      b, bb : 2nd transition indices
        //                      cutoff: Value at which interaction is taken to be zero
        // Output               void  : Arrays strsJ, trnsJ, strsL, and trnsL are filled.
        //                              The values of Jterms and Lterms is adjusted
        // Note                       : T1s, T2s are assumed in their proper bases
        // Note                       : If J(L)terms != 0 at start, terms are added to arrays
        // Note                       : The subscripts 1 and 2 are cryptic in this context
        //                              as they are often used simultaneously for indexing
        //                              both a mechanism and the spin(s) involved in relaxation
        //
        //                                          Mechanism     Rank  1    2
        //                                      ----------------- ---- ---- ----
        //                                       Dipolar           2    ij   kl
        //                                       Shift Anisotropy  2    i    j
        //                                       Quadrupolar       2    i    j
        //                                       Random Field      1    i    j
 
/* The purpose of this routine is determine all non-zero contributions to the relaxation
   matrix element <a,a'|R|b,b'> for the interactions 1 with 2 as specified by the two
   input spin tensors (via their components), T1s and T2s.  The contribution to <a,a'|R|b,b'>
   is given by
                      rank    ls
   <a,a'|R   |b,b'>   ---   [ ---
          1,2         \     | \                 m       m
   ---------------- = /     | /   delta     <a|T |g><b|T |g> J  (w  )
        Xi            ---   | ---      a',b'    1       2    ~12  gb
          1,2       m=-rank [  g
  
                                 m        m                       m        m
                           - <a|T |b><a'|T |b'> J  (w    ) - <b'|T |a'><b|T |a> J  (w  )
                                 1        2     ~12  b'a'         1        2    ~12  ab
 
                                                     ls
                                                     ---                                      ]
                                                     \               m        m               |
                                                  +  /   delta   <g|T |a'><g|T |b'> J  (w   ) |
                                                     ---      a,b    1        2     ~12  b'g  |
                                                      g                                       ]
  
   where the nomenclature is J  (w) = J  (w) + L  (w), and L  (w) = w*tau*J  (w).  L(w) are
                             ~12       12       12          12             12
  
   the antisymmetric dynamic frequency shift terms whereas J(w) are the normal symmetric
   spectral density functions. The indices where a,a',b, and b' indicate different eigenstates,
   and there are implicit spin labels on both the tensor components and on the spectral density
   functions due to the mechanism indices 1 and 2.  The expression is reduced in this function
   into a sum over nterms according to
  
                           < Jterms                        < Lterms
          <a,a'|R   |b,b'>   ---                             ---
                 1,2         \                               \
          ---------------- = /     strsJ[i] * J(w        ) + /     strsL[j] * L(w        )
              Xi             ---                 trnsJ[i]    ---                 trnsL[j]
                1,2          i=0                             j=0
  
   where trans[i] relates to the original eigenstates according to
  
                              w           w        trnsJ[i] = a*hs+a'
                               trnsJ[i] =  aa'           a' = tr%hs;
                                                         a  = (tr-a')/hs;
  
   There will be (2*rank+1)*[ls + 2 + ls] terms in the summation over m & g, however the spin tensor
   products act a a good filter so that most terms will be zero.  Any terms in which the spin tensor
   product has a magnitude below cutoff are herein considered zero and will not be included in the
   ouput arrays. Furthmore, many terms will add if they have the same spectral density functions.
   The end effect of this routine then, is the smallest possible amount of information needed to
   characterize the relaxation the relaxation matrix element for the pair of interactions (1,2) and
   the specified set of spins in terms of scaled spectral density functions.                       */
 

MSVCDLL void Rel_12(int hs, int rank, gen_op* T1s, gen_op* T2s,
                  int& nterms, complex* strs, int* trns,
		 	  int a, int aa, int b, int bb, double cutoff=1.e-4);

        // Input                hs    : Spin system Hilbert space
        //                      rank  : Rank of the two interactions
        //                      T1s   : Spin tensor components, 1st interaction
        //                      T2s   : Spin tensor components, 2nd interaction
        //                      strs  : Array of interaction strengths
        //                      trns  : Array of transition indices
        //                      a, b  : 1st transition indices
        //                      a, aa : 1st transition indices
        //                      b, bb : 2nd transition indices
        // Output               void  : Arrays strs and trns are filled, nterms value set
        // Note                       : T1s, T2s are assumed in their proper bases
        // Note                       : If nterms != 0 at start, terms are added to arrays
        // Note                       : The subscripts 1 and 2 are cryptic in this context
        //                              as they are often used simultaneously for indexing
        //                              both a mechanism and the spin(s) involved in relaxation
        //
        //                                          Mechanism     Rank  1    2
        //                                      ----------------- ---- ---- ----
        //                                       Dipolar           2    ij   kl
        //                                       Shift Anisotropy  2    i    j
        //                                       Quadrupolar       2    i    j
        //                                       Random Field      1    i    j
 
/* The purpose of this routine is determine all non-zero contributions to the relaxation
   matrix element <a,a'|R|b,b'> for the interactions 1 with 2 as specified by the two
   input spin tensors (via their components), T1s and T2s.  The contribution to <a,a'|R|b,b'>
   is given by
                      rank    ls
   <a,a'|R   |b,b'>   ---   [ ---
          1,2         \     | \                 m       m
   ---------------- = /     | /   delta     <a|T |g><b|T |g> J  (w  )
        Xi            ---   | ---      a',b'    1       2     12  gb
          1,2       m=-rank [  g
  
                                 m        m                       m        m
                           - <a|T |b><a'|T |b'> J  (w    ) - <b'|T |a'><b|T |a> J  (w  )
                                 1        2      12  b'a'         1        2     12  ab       
 
                                                     ls
                                                     ---                                      ]
                                                     \               m        m               |
                                                  +  /   delta   <g|T |a'><g|T |b'> J  (w   ) |
                                                     ---      a,b    1        2      12  b'g  |
                                                      g                                       ]
 
   where a,a',b, and b' indicate different eigenstates, and there are implicit spin labels on
   both the tensor components and on the spectral density functions due to the mechanism indices
   1 and 2.  The expression is reduced in this function into a sum over nterms according to
  
                                            < nterms
                           <a,a'|R   |b,b'>   ---
                                  1,2         \
                           ---------------- = /      strs[i] * J(w        )
                                Xi            ---                 trans[i]
                                  1,2         i=0
  
   where trans[i] relates to the original eigenstates according to
  
                              w           w         trns[i] = a*hs+a'
                               trans[i] =  aa'           a' = tr%hs;
                                                         a  = (tr-a')/hs;
  
   There will be (2*rank+1)*[ls + 2 + ls] terms in the summation over m & g, however the spin tensor
   products act a a good filter so that most terms will be zero.  Any terms in which the spin tensor
   product has a magnitude below cutoff are herein considered zero and will not be included in the
   ouput arrays.

   The end effect of this routine then, is to generate the information to construct the relaxation 
   matrix element for a pair of interactions (1,2) and specified set of spins in terms of scaled
   spectral density functions.  For example, for dipole-dipole relaxation in a CH system one might
   obtain

  
                     <1,0|R|1,0>          DD                 DD              DD
                     -----------  = -0.5 J    (WC-WH) + 0.2 J    (WH) + 0.3 J    (-WC+WH)
                       Xi Xi              CHCH               CHCH            CHCH
                        CH  CH
  
   Other routines in this module are provided to combine like terms and add these terms
   to terms for the same element from other spins and other mechanisms.                          */


MSVCDLL void Rel_12_condense(int hs, int ntermsi, int& nterms,
                            complex* strs, int* trns, int anti=0, double cutoff=1.e-4);
 
        // Input                hs    : Spin system Hilbert space
        //                      ntermi: Initial index to start with
        //                      nterms: Number of terms contributing to R element
        //                      strs  : Array of interaction strengths
        //                      trns  : Array of transition indices
        //                      anti  : Flag for J(w) versus L(w)
        //                                   0 : Assume J(-w) = J(w)  symmetric
        //                                  !0 : Assume L(-w) = -L(w) anti-symmetric
        //                      cutoff: Value at which interaction is taken to be zero
        // Output               void  : Arrays strs and trns are modified, nterms value reset
        //                              as terms are combined with the same J(w) values
        // Note                       : IT IS ASSUMED HEREIN THAT ALL COMPONENTS FROM
        //                              [ntermsi, nterms) STEM FROM THE SAME RELAXATION
        //                              MECHANISMS AND INVOLVE THE SAME SPINS
        // Note                       : This condenses according to a transition number
        //                              index only.  Transitions with the same frequency
        //                              but of a different index are not added herein
        //                              (for example WH1 will not combine with WH2)
 
//     The intent of this routine is to reduce the number of explicit terms contributing
//  to a particular relaxtion matrix element. For example, the following will be performed.
//
//       DD                 DD              DD                     DD                 DD
// -0.5 J    (WC-WH) + 0.2 J    (WH) + 0.3 J    (-WC+WH) --> -0.2 J    (WC-WH) + 0.2 J    (WH)
//       CHCH               CHCH            CHCH                   CHCH               CHCH
//
// Note that it is assumed that all input terms correspond to the same relaxation mechanisms,
// (in this case DD) and that they are all associated with the same set of spins (here CHCH).
//     Furthermore, it is assumed that the spectral densities are either symmetric, i.e.
//
//                         DD               DD                   DD
//                   c1 * J    (WH) + c2 * J    (-WH) = [c1+c2] J    (WH)
//                         CHCH             CHCH                 CHCH
//
//     or anti-symmetric, such as in the case of the terms for dynamic frequency shifts.
//
//                         DD               DD                   DD
//                   c1 * L    (WH) + c2 * L    (-WH) = [c1-c2] L    (WH)
//                         CHCH             CHCH                 CHCH
//
//      The inherent symmetry is indicated by proper setting of the input flag anti.
//       
// The spectral densities are input in a compact form, where the coefficients are contained
// in the array strs (interaction strengths) and the transitions are indicated by the array
//  trns (transition indices).  These arrays are allowed to contain other contributions to
//  the relaxation matrix element, either from other spins or from other mechanisms.  Thus,
//  the first term to consider is specified by the value ntermi and the number of terms to
//  consider is specified by the value of nterms.  The value of nterms may herein be reduced.

 
MSVCDLL void Rel_12_condense(int hs, int ntermi, int& nterms, complex* strs,
                                             int* trns, std::string* wlabs, int anti=0);
 
        // Input                hs    : Spin system Hilbert space
        //                      ntermi: Initial index to start with
        //                      nterms: Number of non-zero terms contributing
        //                              to R element
        //                      strs  : Array of interaction strengths
        //                      trns  : Array of transition indices
        //                      anti  : Flag for J(w) versus L(w)
        //                                   0 : Assume J(-w) = J(w)  symmetric
        //                                  !0 : Assume L(-w) = -L(w) anti-symmetric
        // Output               void  : Arrays strs and trns are modified as
        //                              all terms are set with positive w in J(w)
        // Note                       : IT IS ASSUMED HEREIN THAT ALL COMPONENTS FROM
        //                              [ntermsi, nterms) STEM FROM THE SAME RELAXATION
        //                              MECHANISMS AND INVOLVE THE SAME SPINS
        // Note                       : This may modifiy the transition index based on
        //                              the transition label in Wlbls.  All transitions
        //                              are set so that the first frequency is positive.
        //                              Then all with the same label are set to have the
        //                              same transition index.
 
// The intent of this routine is set all explicit terms contributing to a particular relaxtion
// matrix element so that the first eigenstate of the transition is positive. For example, the
// following will be performed.
//
//              DD              DD                     DD              DD
//         0.2 J    (WH) + 0.3 J    (-WC+WH) -->  0.2 J    (WH) - 0.3 J    (WC-WH)
//              CHCH            CHCH                  CHCH             CHCH
//
// Note that it is assumed that all input spectral densities are either symmetric, i.e.
//
//                                DD               DD
//                           c * J    (-WH) = c * J    (WH)
//                                CHCH             CHCH
//
// The inherent symmetry is indicated by proper setting of the input flag anti.
//
// The spectral densities are input in a compact form, where the coefficients are contained
// in the array strs (interaction strengths).  The transitions are indicated by the array
// trns (transition indices) and the array Wlbls (transition labels).  These arrays are
// allowed to contain other contributions to the relaxation matrix element, either from other
// spins, other mechanisms, or the opposite symmetry.  Thus, the first term to consider is
// specified by the value ntermi and the number of terms to consider is specified by the value
// of nterms.
 
// ______________________________________________________________________
// ************* STRING RELAXATION MATRIX ELEMENT FUNCTIONS *************
// ______________________________________________________________________

MSVCDLL void Rel(int ntermi, int& nterms, int npairs, complex* strs, int* trns,
		 std::string* wlabs, int* cont, std::string* spns, std::string* Jlbs,
		 std::string& Mlabel, std::string* line1, std::string* line2, std::string* line3, int* signs);

        // Input                ntermi: Number of terms currently in line1, line2, line3, signs
        //                      nterms: Number of terms contributing to <a,aa|R|b,bb>
        //                      npairs: Number spin parings in these terms
        //                      strs  : Array for interaction strengths for these terms
        //                      trns  : Array for transition indices for these terms
        //                      Wlbls : Array of transition labels for the system
        //                      cont  : Array of contributions per spin pairing
        //                      spns  : Array of spin labels per spin pairing
        //                      Jlbs  : Array of J labels per 1,2 pair
        //                      Mlabel: String label fo two mechanisms
        //                      cutoff: Value at which interaction is taken to be zero
        // Output               void  : The contribution to <a,a'|R|b,b'> from the
        //                              mechanisms Mlabel and spins spns are added to the arrays
        //                              line1, line2, line3, and signs starting at index ntermi.
        //                              The number of terms in these arrays, nterms, is updated
        // Note                       : ALL TERMS ARE ASSUMED FROM THE SAME RELAXATION MECHANISM

//  The purpose of this routine is to fill String arrays with a formal (non-numerical)
//  description of the components to a particular relaxation matrix element stemming from
//  the interaction of two relaxation mechanisms (DD, SASA, DSA, SAD, ....).  There may
//  many terms involved due to a loop over spin pairs or pairs of spin pairs as well as a
//  loop over transitions (spectral density functions).  Fore example, the contribution
//  from the shift anisotropy (SA) to <0,1|R|0,1> in a CH system is
  
//               SASA                SASA             SASA            SASA    
//         <0,1|R    |0,1> = + 0.50 J    (WC) + 0.67 J    (0) + 0.50 J    (WH)
//                                   CC               HH              HH      
//
//                                   SASA            SASA    
//                           + 0.67 L    (0) - 0.50 L    (WH)
//                                   HH              HH     
 
//  That is accomplished relative to the function arguments according to the scheme
//  presented in the expression:
//
//                            nterms
//                             ---                                      ---
//              Mlabel         \                      Mlabel            \             line1[i]
//       <a,a'|R      |b,b'> = /    strs[i] * Jlbls[j]       (w[i]) --> /     sign[i] line2[i]
//              1,2            ---                    spns[j]           ---           line3[i]
//                           i=termi                                  i=termi
//
//  where the index j is incremented (the J symbol used may be changed) for each new combination
//  of spins involved.  The number of unique spin combinations is set by npairs and the number of
//  terms per each npair value is cont[npair].  This allows the user to externally specify the
//  J symbol used according to the variety of standards found in NMR literature (J,L,K,Q,j,k,..)
//  and specify the spin labels to use (e.g. in dipolar relaxation ijij vs. ij; ijik vs. jik).
//  Because the user may desire to combine these terms with other terms for the same relaxation
//  matrix element (due to other interacions), the routine puts first term of the nterms input into
//  the string arrays sign, line1, line2, and line3, beginning at index i.

// ______________________________________________________________________
// *************** RELAXATION ANALYSIS AUXILIARY FUNCTIONS **************
// ______________________________________________________________________


MSVCDLL void Spin_labels(std::string* Lbls, const spin_sys& sys, int index=0);

        // Input                Lbls  : String a for spin labels
        //                      sys   : A spin system
        //                      index : Flag for inclusion of spin index
        //                                  0 = Use only if needed
        //                                  1 = Always use spin indices
        //                                 <1 = Never use spin indices
        // Output               void  : Fills the array of spin labels


// The intent of this routine is to                     DD
// construct labels for the spins which can       -0.5 J    (WC-WH)
// be used to produce subscripts in J labels -------->  CHCH

// In this routine, each spin of the system is assigned a unique spin label.
//  This will be as small as possible yet still maintiain it's uniqueness. 
// In the worst case, the label will be the spin's symbol and a spin index.
// For example: 14N0, 14N1, 14N2.  However, if there are no other isotopes of a
// particular spin type, the isotope label will be neglected and 14N0 will
// just be N0.  Next, if there are no other spins of a particular isotope type
// the spin index will be neglected and only N is used.  If the user is familiar
// with the spin system then there can be no question of the spin by the returned
// index no matter what their ordering.  One other tidbit, rather that use 1H
// 2H, and 3H for the hydrogen isotopes, here they are switched to H,D,T


MSVCDLL void W_labels(std::string* Wlabels, const spin_sys& sys, gen_op &Op, int index=-1);

	// Input 	Wlabels : An array of strings
	//       	sys	: A basic spin system
	// 		Op	: General operator (associated to sys)
	//			index : Flag for inclusion of spin index
	//				   >1 = Always use spin indices
	//				    0 = Use spin indices if needed
	//				   -1 = Never use spin indices (default)
	//				   -2 = Assume extreme narrowing
	//				   -3 = Assume extreme broadening
	//				    
	// Return	void    : The string array Wlabels is filled with system
	//			  transitions ordered from the eigenbasis
	//			  of operator Op in the format specified
	// Note			: Transition ordering depends on Op input basis!
	// Note			: For Extreme Narrowing All Labels are W
	// Note			: For Extreme Broadening Labels not Near 0 are X
	// Note			: The view of E.B. and E.N. should be reworked
	//			  because it is currently very in-efficient

// The intent of this routine is to                     DD
// construct labels for the transitions which      -0.5 J    (WC-WH)
// can be used in construction of spectral               CHCH   ^
// density functions ___________________________________________|

// In this routine, each transition in the system is assigned a transition label.
// Such labels are NOT unique because differences due to chemical shifts and
// scalar couplings are neglected. Thus, for example, all single quantum proton
// transitions will be called WH even though their frequencies may vary slightly
// due to shifts and couplings.  What the transitions are depends upon the current
// basis of the input operator Ho (labeled as such because it is usually a
// Hamiltonian and the basis it's eigenbasis).


MSVCDLL void Elem_labels(std::string* Lbls, std::string& R, std::string& M, std::string& S,
            int a, int aa, int b, int bb, int la=0, int laa=0, int lb=0, int lbb=0);

	// Input 	       Lbls   : Strings for relaxation matrix
	//			 	element labels
	//	               R      : Relaxation matrix label
	//	               M      : Mechanisms label
	//		       S      : Spins label
        //                     a, aa  : 1st transition indices
        //                     b, bb  : 2nd transition indices
        //                     la,laa : Desired string length of the
	//			        1st transition indices
        //                     laa,lbb: Desired string length of the
	//				2nd transition indices
	// Output		void  : Fills the string array with lines
	//				of the element label
	// Note			      : Set strings Mlabel and Slabel to
	//				"" if not to be printed.
	//			         
	// Note			      : If R is only blanks, then three blank
	//				strings of its length are output

// The intent of this routine is to           DD         <--- Lbl[0]
// construct labels for a particular   <a,a'|R    |b,b'> <--- Lbl[1]
// relaxation matrix element.  This           CHCH       <--- Lbl[2]
// one at the right, for example.

MSVCDLL void Rel_clean(gen_op* T1s, gen_op* T2s, int rank);

        // Input                T1s   : Array of general operators
        //                      T2s   : Array of general operators
        //                      rank  : Rank of a tensor
        // Output               void  : The arrays T1s and T2s contain
        //                              the components of two tensors
        //                              of rank "rank".  These are zeroed.
        // Note                       : This routine is necessary if 
        //                              memory is allocated by use of the
        //                              "new" function in setting up the
        //                              two tensor component arrays.

// ______________________________________________________________________
// *************** RELAXATION MATRIX OUTPUT STRING FUNCTIONS ************
// ______________________________________________________________________

 
MSVCDLL void Rel(std::ostream& ostr, int nterms, std::string* line1, std::string* line2,
         std::string* line3, int* signs, std::string* Elabel, int add=0, int ncols=4);


// ______________________________________________________________________
// ************ MISCELLANEOUS USEFUL FUNCTIONS FOR RELAXATION ***********
// ______________________________________________________________________

// sosi - some if not all of these functions will find an different home
//        in GAMMA some day 

// MSVCDLL void ask_relax(int argc, char* argv[], int& argn,
//                 super_op& R, const sys_dynamic& sys, gen_op& H, int pflag=1);

        // Input        argc    : Number of command line arguments
        //              argv    : Command line arguments
        //              argn    : Initial command line argument
        //                        for relaxation parameters
        //              R       : Relaxation superoperator
        //              sys     : Dynamic spin system
        //              H       : Isotropic Hamiltonian
	//		pflag   : Flag for printing computaions
        // Output       none    : Function is void.  The relaxation
        //                        matrix R has different effects added
        //                        in depending upon user requests


MSVCDLL void sort(int* indx, matrix& mx, int k, int type=0, int colf=0);

        // Input        indx  : An array of integers
	//		mx    : A matrix
 	//		k     : An index
	//		type  : Basis for the sorting
	//			 0 = sort real values (default)
	//			>0 = sort norms
	//		        <0 = sort imaginaries
	//		colf  : Row or column sort
	//			 0 = sort rows (default)
	//			!0 = sort columns
	// Output	void  : The integers in indx will be filled
	//			with the sorted indices
	// Note		      : The matrix is left unaltered


#endif										// RelaxAnalyze.h 

