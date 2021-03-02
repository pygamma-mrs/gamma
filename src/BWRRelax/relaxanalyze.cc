// sosi - Rijkl_el (etc.) one can use a function for spns[ijkl] to produce the
//        compact W&G nomenclature for the spin labels
// sosi - Rijkl_el (etc.) one can use a function for Jlbs[ijkl] to produce
//        any combination of labels for J
// sosi - Rel? : one could also track xi values according for each spin
//               set the same way the J's are done. That would also remove the
//	         problem of setting the precision

/* relax_analyze.cc ************************************-*-c++-*-*
**								**
**	                   G A M M A	 			**
**							 	**
**	NMR Relaxation Analysis             Implementation	**
**							 	**
**	Copyright (c) 1993			 		**
**	Scott Smith				 		**
**	University of Utah					**
**	Department of Chemistry					**
**      Salt Lake City, UT, 84112, USA                          **
**							 	**
*****************************************************************/

/*****************************************************************
**							 	**
**  Description						 	**
**							 	**
**  This module of the GAMMA relaxation package provides	**
**  routines designed to facilitate the analysis of relaxation	**
**  in arbitrary spin systems.					**
**						 		**
*****************************************************************/

#ifndef _relax_analyze_cc_	// Is this file already included?
#define _nmr_ham_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma implementation		// This is the implementation

#endif
#if defined(_MSC_VER)			// If using MSVC++ then we
 #pragma warning (disable : 4786)       // Kill STL namelength warnings
#endif 

# include <BWRRelax/relaxanalyze.h>
# include <Basics/StringCut.h>
#include <stdlib.h>

// ______________________________________________________________________
// ********* RELAXATION SUPEROPERATOR ELEMENT OUTPUT FUNCTIONS **********
// ______________________________________________________________________

// ------------- Dipole-Dipole Relaxation Matrix Elements ---------------

void RDDel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
	                      int a, int aa, int b, int bb,
                               int DFS, int Windex, int Sindex)

//                                >i     l>k
//                DD         --- --- --- ---        DD
//         <a,aa|R  |b,bb> = \   \   \   \   <a,aa|R    |b,bb>
//                           /   /   /   /          ijkl
//                           --- --- --- ---
//                            i   j   k   l

	// Input 		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			T1    : Spin tensor, mu1 for dipolar 
	//			T2    : Spin tensor, mu2
	//			a, b  : 1st transition indices
	//			aa, bb: 2nd transition indices
	//			DFS   : Flag for dynamic frequency shifts
	//					> 0  J Terms Only				    
	//					  0  Both J and L Terms (Default)
	//					< 0  L Terms Only
	//			Windex: Flag for frequency indexing
        //                                 >1 = Always use spin indices
        //                                  0 = Use spin indices if needed
        //                                 -1 = Never use spin indices (default)
        //                                 -2 = Assume extreme narrowing
        //                                 -3 = Assume extreme broadening
	//			Sindex: Flag for spin indexing
        //                                  0 = Use only if needed (default)
        //                                  1 = Always use spin indices
        //                                 <1 = Never use spin indices
	// Output		void  :
	// Note			      : This routine assumes that each of
	//				the two interactions are rank 2
	//				and that they both involve dipolar
	//				interactions
// sosi - then the input here only needt T1 for auto, not T1 and T2
//        and the same goes for all the auto correlated stuff

   {
   std::string line1[100];				// Array of strings for line 1
   std::string line2[100];				// Array of strings for line 2
   std::string line3[100];				// Array of strings for line 3
   int signs[100];				// Array of flags for the signs
   std::string Elabel[3];				// strings for element label, one each line
   std::string Rlabel = "R";				// Label for relaxation matrix
   std::string Mlabel = "DD";			// Label for relaxation mechanisms
   std::string Slabel = "";				// Label for spin involved (none here)
   Elem_labels(Elabel, Rlabel, Mlabel, Slabel,	// Get a label for the element
                                 a, aa, b, bb);
   int Jterms=0, Lterms=0;			// Start with no terms (J(w) or L(w)) 
   Rijkl_el(sys, Ho, 2, T1, T2, Mlabel,		// Get all the strings for the terms
    line1, line2, line3, signs, Jterms, Lterms,	// that contribute to this particular element
                  a,aa,b,bb,DFS,Windex,Sindex);
   Rel(std::cout, Jterms+Lterms, line1, line2,	// Now output the element in string format
                         line3, signs, Elabel);
   return;
   }


void RSSel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
	                      int a, int aa, int b, int bb,
                               int DFS, int Windex, int Sindex)

//                SASA         --- ---        SASA
//         <a,aa|R    |b,bb> = \   \   <a,aa|R    |b,bb>
//                             /   /          ij
//                             --- --- 
//                              i   j

	// Input 		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			a, b  : 1st transition indices
	//			aa, bb: 2nd transition indices
	//			DFS   : Flag for dynamic frequency shifts
	//					> 0  J Terms Only				    
	//					  0  Both J and L Terms (Default)
	//					< 0  L Terms Only
	//			Windex: Flag for frequency indexing
        //                                 >1 = Always use spin indices
        //                                  0 = Use spin indices if needed
        //                                 -1 = Never use spin indices (default)
        //                                 -2 = Assume extreme narrowing
        //                                 -3 = Assume extreme broadening
	//			Sindex: Flag for spin indexing
        //                                  0 = Use only if needed (default)
        //                                  1 = Always use spin indices
        //                                 <1 = Never use spin indices
	// Output		void  :
	// Note			      : This routine assumes that each of
	//				the two interactions are rank 2
	//				and that they both involve 1 spin
	//				e.g. SA-SA, Q-Q, Q-SA, SA-Q, RF-RF

   {
   std::string line1[100];				// Array of strings for line 1
   std::string line2[100];				// Array of strings for line 2
   std::string line3[100];				// Array of strings for line 3
   int signs[100];				// Array of flags for the signs
   std::string Elabel[3];				// strings for element label, one each line
   std::string Rlabel = "R";				// Label for relaxation matrix
   std::string Mlabel = "SASA";			// Label for relaxation mechanisms
   std::string Slabel = "";				// Label for spin involved (none here)
   Elem_labels(Elabel, Rlabel, Mlabel, Slabel,	// Get a label for the element
                                 a, aa, b, bb);
   int Jterms=0, Lterms=0;
   Rij_el(sys, Ho, 2, T1, T2, Mlabel, 
    line1, line2, line3, signs, Jterms, Lterms,
                  a,aa,b,bb,DFS,Windex,Sindex);
   Rel(std::cout, Jterms+Lterms, line1, line2,
                         line3, signs, Elabel);
   return;
   }


void RDSel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
	                      int a, int aa, int b, int bb,
                               int DFS, int Windex, int Sindex)

//                                >i
//                DSA         --- --- ---        DSA
//         <a,aa|R   |b,bb> = \   \   \   <a,aa|R   |b,bb>
//                            /   /   /          ijk
//                            --- --- --- 
//                             i   j   k

   {
   std::string line1[100];				// Array of strings for line 1
   std::string line2[100];				// Array of strings for line 2
   std::string line3[100];				// Array of strings for line 3
   int signs[100];				// Array of flags for the signs
   std::string Elabel[3];				// strings for element label, one each line
   std::string Rlabel = "R";				// Label for relaxation matrix
   std::string Mlabel = "DSA";			// Label for relaxation mechanisms
   std::string Slabel = "";				// Label for spin involved (none here)
   Elem_labels(Elabel, Rlabel, Mlabel, Slabel,	// Get a label for the element
                                 a, aa, b, bb);
   int Jterms=0, Lterms=0;
   Rijk_el(sys, Ho, 2, T1, T2, Mlabel, 
    line1, line2, line3, signs, Jterms, Lterms,
                  a,aa,b,bb,DFS,Windex,Sindex);
   Rel(std::cout, Jterms+Lterms, line1, line2,
                         line3, signs, Elabel);
   return;
   }


void RSDel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
	                      int a, int aa, int b, int bb,
                               int DFS, int Windex, int Sindex)

//                                    >i
//                SAD         --- --- ---        SAD
//         <a,aa|R  |b,bb> = \   \   \    <a,aa|R   |b,bb>
//                           /   /   /           kij
//                           --- --- ---
//                            k   i   j

   {
   std::string line1[100];				// Array of strings for line 1
   std::string line2[100];				// Array of strings for line 2
   std::string line3[100];				// Array of strings for line 3
   int signs[100];				// Array of flags for the signs
   std::string Elabel[3];				// strings for element label, one each line
   std::string Rlabel = "R";				// Label for relaxation matrix
   std::string Mlabel = "SAD";			// Label for relaxation mechanisms
   std::string Slabel = "";				// Label for spin involved (none here)
   Elem_labels(Elabel, Rlabel, Mlabel, Slabel,	// Get a label for the element
                                 a, aa, b, bb);
   int Jterms=0, Lterms=0;
   Rkij_el(sys, Ho, 2, T1, T2, Mlabel, 
    line1, line2, line3, signs, Jterms, Lterms,
                  a,aa,b,bb,DFS,Windex,Sindex);
   Rel(std::cout, Jterms+Lterms, line1, line2,
                         line3, signs, Elabel);
   return;
   }


void RRRel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
	                      int a, int aa, int b, int bb,
                               int DFS, int Windex, int Sindex)

//                RFRF         --- ---        RFRF
//         <a,aa|R    |b,bb> = \   \   <a,aa|R    |b,bb>
//                             /   /          ij
//                             --- ---
//                              i   j

   {
   std::string line1[100];				// Array of strings for line 1
   std::string line2[100];				// Array of strings for line 2
   std::string line3[100];				// Array of strings for line 3
   int signs[100];				// Array of flags for the signs
   std::string Elabel[3];				// strings for element label, one each line
   std::string Rlabel = "R";				// Label for relaxation matrix
   std::string Mlabel = "RFRF";			// Label for relaxation mechanisms
   std::string Slabel = "";				// Label for spin involved (none here)
   Elem_labels(Elabel, Rlabel, Mlabel, Slabel,	// Get a label for the element
                                 a, aa, b, bb);
   int Jterms=0, Lterms=0;
   Rij_el(sys, Ho, 1, T1, T2, Mlabel, 
    line1, line2, line3, signs, Jterms, Lterms,
                  a,aa,b,bb,DFS,Windex,Sindex);
   Rel(std::cout, Jterms+Lterms, line1, line2,
                         line3, signs, Elabel);
   return;
   }


void RQQel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
	                      int a, int aa, int b, int bb,
                               int DFS, int Windex, int Sindex)

//                QQ         --- ---        QQ
//         <a,aa|R  |b,bb> = \   \   <a,aa|R  |b,bb>
//                           /   /          ij
//                           --- --- 
//                            i   j

   {
   std::string line1[100];				// Array of strings for line 1
   std::string line2[100];				// Array of strings for line 2
   std::string line3[100];				// Array of strings for line 3
   int signs[100];				// Array of flags for the signs
   std::string Elabel[3];				// strings for element label, one each line
   std::string Rlabel = "R";				// Label for relaxation matrix
   std::string Mlabel = "QQ";			// Label for relaxation mechanisms
   std::string Slabel = "";				// Label for spin involved (none here)
   Elem_labels(Elabel, Rlabel, Mlabel, Slabel,	// Get a label for the element
                                 a, aa, b, bb);
   int Jterms=0, Lterms=0;
   Rij_el(sys, Ho, 2, T1, T2, Mlabel, 
    line1, line2, line3, signs, Jterms, Lterms,
                  a,aa,b,bb,DFS,Windex,Sindex);
   Rel(std::cout, Jterms+Lterms, line1, line2,
                         line3, signs, Elabel);
   return;
   }
	

void RQSel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
	                      int a, int aa, int b, int bb,
                               int DFS, int Windex, int Sindex)

//                QSA         --- ---        QSA
//         <a,aa|R   |b,bb> = \   \   <a,aa|R   |b,bb>
//                            /   /         ij
//                            --- ---
//                             i   j

   {
   std::string line1[100];				// Array of strings for line 1
   std::string line2[100];				// Array of strings for line 2
   std::string line3[100];				// Array of strings for line 3
   int signs[100];				// Array of flags for the signs
   std::string Elabel[3];				// strings for element label, one each line
   std::string Rlabel = "R";				// Label for relaxation matrix
   std::string Mlabel = "QSA";			// Label for relaxation mechanisms
   std::string Slabel = "";				// Label for spin involved (none here)
   Elem_labels(Elabel, Rlabel, Mlabel, Slabel,	// Get a label for the element
                                 a, aa, b, bb);
   int Jterms=0, Lterms=0;
   Rij_el(sys, Ho, 2, T1, T2, Mlabel, 
    line1, line2, line3, signs, Jterms, Lterms,
                  a,aa,b,bb,DFS,Windex,Sindex);
   Rel(std::cout, Jterms+Lterms, line1, line2,
                         line3, signs, Elabel);
   return;
   }


void RSQel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
	                      int a, int aa, int b, int bb,
                               int DFS, int Windex, int Sindex)

//                SAQ         --- ---        SAQ
//         <a,aa|R   |b,bb> = \   \   <a,aa|R   |b,bb>
//                            /   /          ij
//                            --- --- 
//                             i   j

   {
   std::string line1[100];			// Array of strings for line 1
   std::string line2[100];			// Array of strings for line 2
   std::string line3[100];			// Array of strings for line 3
   int signs[100];				// Array of flags for the signs
   std::string Elabel[3];			// strings for element label, one each line
   std::string Rlabel = "R";			// Label for relaxation matrix
   std::string Mlabel = "SAQ";			// Label for relaxation mechanisms
   std::string Slabel = "";			// Label for spin involved (none here)
   Elem_labels(Elabel, Rlabel, Mlabel, Slabel,	// Get a label for the element
                                 a, aa, b, bb);
   int Jterms=0, Lterms=0;
   Rij_el(sys, Ho, 2, T1, T2, Mlabel, 
    line1, line2, line3, signs, Jterms, Lterms,
                  a,aa,b,bb,DFS,Windex,Sindex);
   Rel(std::cout, Jterms+Lterms, line1, line2,
                         line3, signs, Elabel);
   return;
   }


void RQDel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
	                      int a, int aa, int b, int bb,
                               int DFS, int Windex, int Sindex)

//                QD         --- ---        QD
//         <a,aa|R  |b,bb> = \   \   <a,aa|R  |b,bb>
//                           /   /          ij
//                           --- ---
//                            i   j

   {
   std::string line1[100];				// Array of strings for line 1
   std::string line2[100];				// Array of strings for line 2
   std::string line3[100];				// Array of strings for line 3
   int signs[100];				// Array of flags for the signs
   std::string Elabel[3];				// strings for element label, one each line
   std::string Rlabel = "R";				// Label for relaxation matrix
   std::string Mlabel = "QD";			// Label for relaxation mechanisms
   std::string Slabel = "";				// Label for spin involved (none here)
   Elem_labels(Elabel, Rlabel, Mlabel, Slabel,	// Get a label for the element
                  a,aa,b,bb,DFS,Windex,Sindex);
   int Jterms=0, Lterms=0;
   Rkij_el(sys, Ho, 2, T1, T2, Mlabel, 
    line1, line2, line3, signs, Jterms, Lterms,
                                 a, aa, b, bb);
   Rel(std::cout, Jterms+Lterms, line1, line2,
                         line3, signs, Elabel);
   return;
   }


void RDQel(const sys_dynamic& sys, gen_op& Ho, spin_T* T1, spin_T* T2,
	                      int a, int aa, int b, int bb,
                               int DFS, int Windex, int Sindex)

//                DQ         --- ---        DQ
//         <a,aa|R  |b,bb> = \   \   <a,aa|R  |b,bb>
//                           /   /          ij
//                           --- ---
//                            i   j

   {
   std::string line1[100];				// Array of strings for line 1
   std::string line2[100];				// Array of strings for line 2
   std::string line3[100];				// Array of strings for line 3
   int signs[100];				// Array of flags for the signs
   std::string Elabel[3];				// strings for element label, one each line
   std::string Rlabel = "R";				// Label for relaxation matrix
   std::string Mlabel = "DQ";			// Label for relaxation mechanisms
   std::string Slabel = "";				// Label for spin involved (none here)
   Elem_labels(Elabel, Rlabel, Mlabel, Slabel,	// Get a label for the element
                                 a, aa, b, bb);
   int Jterms=0, Lterms=0;
   Rijk_el(sys, Ho, 2, T1, T2, Mlabel, 
    line1, line2, line3, signs, Jterms, Lterms,
                  a,aa,b,bb,DFS,Windex,Sindex);
   Rel(std::cout, Jterms+Lterms, line1, line2,
                         line3, signs, Elabel);
   return;
   }


// ______________________________________________________________________
// *********** RELAXATION MATRIX WHOLE ELEMENT STRING FUNCTIONS *********
// ______________________________________________________________________


  void Rijkl_el(const sys_dynamic& sys, gen_op& Ho,
            int rank, spin_T* T1, spin_T* T2, std::string& Mlabel,
         std::string* line1, std::string* line2, std::string* line3, int* signs,
              int& Jterms, int& Lterms,
                              int a, int aa, int b, int bb,
                               int DFS, int Windex, int Sindex)

	// Input 		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			rank  : Rank of the interactions
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			Mlabel: Label for two mechanisms
	//			line1 : strings for the upper line of component
	//			line2 : strings for the middle line of component
	//			line3 : strings for the lower line of component
	//			signs : Flags for the signs of the components
	//			a, aa : 1st transition indices
	//			b, bb : 2nd transition indices
	//			DFS   : Flag for dynamic frequency shifts
	//					> 0  J Terms Only				    
	//					  0  Both J and L Terms (Default)
	//					< 0  L Terms Only
	//			Windex: Flag for frequency indexing
        //                                 >1 = Always use spin indices
        //                                  0 = Use spin indices if needed
        //                                 -1 = Never use spin indices (default)
        //                                 -2 = Assume extreme narrowing
        //                                 -3 = Assume extreme broadening
	//			Sindex: Flag for spin indexing
        //                                  0 = Use only if needed (default)
        //                                  1 = Always use spin indices
        //                                 <1 = Never use spin indices
	// Output		void  :
	// Note			      : This routine assumes that each of
	//				the two interactions involve 2 spins
	// Note			      : The above implies the rank is 2

//                                --- --- --- ---
//                    12          \   \   \   \    12
//             <a,a'|R  |b,b'> =  /   /   /   /   R
//	                          --- --- --- ---  ij,kl
//                                 i   j   k   l

   {
std::string Llbl = "L";
std::string Jlbl = "J";
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[2*rank+1];
   T2s = new gen_op[2*rank+1];
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
   int m;					// z-component of ang. momentum
   std::string slabs;
   std::string *Wlbls;				// Array of transition labels
   Wlbls = new std::string[hs*hs];
   W_labels(Wlbls, sys, Ho, Windex);		// Fill up the transition labels
   int ndip = sys.dipoles();			// Total number of dipoles
   int ndpairs=ndip*ndip;			// Total number of dipole pairs
   int ij,kl;
   int maxt = 2*(2*rank+1)*(hs+1);		// Max. contributions per dip. pair

   complex *strsJ;			// Array for interaction strengths (J's)
   int *trnsJ;			// Array for transition indices (J's)
   int *contJ;				// Number of contributors per dipole pair
   std::string *Jlbs;			// Array for J labels (J, can use L,K, etc.)
   complex *strsL;			// Array for interaction strengths (L's)
   int *trnsL;			// Array for transition indices (L's)
   int *contL;				// Number of contributors per dipole pair
   std::string *Llbs;			// Array for L labels (L, can use L,Q, etc.)
   std::string *spns;			// Array for spin labels (ijkl, ijjk, etc.)

   strsJ = new complex[ndpairs*maxt];
   trnsJ = new int[ndpairs*maxt];	// Array for transition indices (J's)
   contJ = new int[ndpairs];		// Number of contributors per dipole pair
   Jlbs  = new std::string[ndpairs];		// Array for J labels (J, can use L,K, etc.)
   strsL = new complex[ndpairs*maxt];	// Array for interaction strengths (L's)
   trnsL = new int[ndpairs*maxt];		// Array for transition indices (L's)
   contL = new int[ndpairs];		// Number of contributors per dipole pair
   Llbs  = new std::string[ndpairs];		// Array for L labels (L, can use L,Q, etc.)
   spns  = new std::string[ndpairs];		// Array for spin labels (ijkl, ijjk, etc.)

   int Jterm0=0, Lterm0=0;			// Total terms contributing to <a,a'|R|b,b'>
   int ijkl = 0;				// Counter for dipole-dipole interactions
   ij=0;
   std::string *Slbls;				// strings for spin labels
   Slbls = new std::string[ns];
   Spin_labels(Slbls, sys, Sindex);		// Fill up spin labels
   for(int i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     {
     for(int j=i+1; j<ns; j++)
       {
       for(m=-rank; m<=rank; m++)		// Put spin tensor for i,j into a
         {					// vector of operators in basis of Ho
         T1s[m+rank] = T1[ij].component(rank,m);
         T1s[m+rank].Op_base(Ho);
         }
       kl = 0;					// dipole kl count to zero
       for(int k=0; k<ns-1; k++)		// Sum spin pairs k&l: dipole kl
         {
         for(int l=k+1; l<ns; l++)
           {
           slabs = Slbls[i] + Slbls[j]
                 + Slbls[k] + Slbls[l];
           if(ij == kl) 			// Auto-correlation term
             {					//                      12
             Rel_12(hs,rank,T1s,T1s, 		// Get terms for <a,a'|R    |b,b'>
                    Jterms,strsJ,trnsJ, 	//			ijkl
                    Lterms,strsL,trnsL,
                              Wlbls,a,aa,b,bb);
             }
           else if(ij != kl)			// Cross-correlation term
             {
             for(m=-rank; m<=rank; m++)		// Put spin tensor for k,l into a
               {				// vector of operators in basis of Ho
               T2s[m+rank] = T2[kl].component(rank,m);
               T2s[m+rank].Op_base(Ho);
               }				//                      12
             Rel_12(hs,rank,T1s,T2s, 		// Get terms for <a,a'|R    |b,b'>
                    Jterms,strsJ,trnsJ, 	//			ijkl
                    Lterms,strsL,trnsL,
                              Wlbls,a,aa,b,bb);
             }
           contJ[ijkl] = Jterms-Jterm0;		// Number of terms this dipole pair
           Jlbs[ijkl] = Jlbl;			// Store J label for these terms
           contL[ijkl] = Lterms-Lterm0;		// Number of terms this dipole pair
           Llbs[ijkl] = Llbl;			// Store L label for these terms
           spns[ijkl] = slabs;			// Store spin label for these terms
           kl++;				// Increment second dipole
           ijkl++;				// Increment total dipole pair count
           Jterm0 = Jterms;			// Update the starting J terms counter
           Lterm0 = Lterms;			// Update the starting L terms counter
           }
         }
       ij++;					// Increment first dipole
       }
     }
   if(Jterms && DFS >= 0)
     Rel(0, Jterms, ndpairs, strsJ, trnsJ,	// Put J(w) terms into string arrays
           Wlbls, contJ, spns, Jlbs, Mlabel, 	// line[1-3], and track signs
                   line1, line2, line3, signs);
   else
     Jterms = 0;
   if(Lterms && DFS <= 0)
     Rel(Jterms, Lterms, ndpairs, strsL, trnsL,	// Add L(w) terms to string arrays
           Wlbls, contL, spns, Llbs, Mlabel, 	// line[1-3], track signs also
                   line1, line2, line3, signs);
   else
     Lterms = 0;
   Rel_clean(T1s, T2s, rank);
   delete [] Wlbls;
   delete [] strsJ;
   delete [] trnsJ;
   delete [] contJ;	
   delete [] Jlbs;	
   delete [] strsL;	
   delete [] trnsL;	
   delete [] contL;	
   delete [] Llbs;		
   delete [] spns;
   delete [] Slbls;

   return;
   }


void Rij_el(const sys_dynamic& sys, gen_op& Ho,
            int rank, spin_T* T1, spin_T* T2, std::string& Mlabel,
        std::string* line1, std::string* line2, std::string* line3, int* signs,
              int& Jterms, int& Lterms,
                              int a, int aa, int b, int bb,
                               int DFS, int Windex, int Sindex)

	// Input 		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			rank  : Rank of the interactions
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			Mlabel: Label for two mechanisms
	//			line1 : strings for the upper line of component
	//			line2 : strings for the middle line of component
	//			line3 : strings for the lower line of component
	//			signs : Flags for the signs of the components
	//			a, aa : 1st transition indices
	//			b, bb : 2nd transition indices
	//			DFS   : Flag for dynamic frequency shifts
	//					> 0  J Terms Only				    
	//					  0  Both J and L Terms (Default)
	//					< 0  L Terms Only
	//			Windex: Flag for frequency indexing
        //                                 >1 = Always use spin indices
        //                                  0 = Use spin indices if needed
        //                                 -1 = Never use spin indices (default)
        //                                 -2 = Assume extreme narrowing
        //                                 -3 = Assume extreme broadening
	//			Sindex: Flag for spin indexing
        //                                  0 = Use only if needed (default)
        //                                  1 = Always use spin indices
        //                                 <1 = Never use spin indices
	// Output		void  :
	// Note			      : This routine assumes that each of
	//				the two interactions involve 1 spin
	//				e.g. SA-SA, Q-Q, Q-SA, SA-Q, RF-RF

//	                 Two Single Spin Mechanisms

//                                        --- ---
//                            12          \   \    12
//                     <a,a'|R  |b,b'> =  /   /   R
//	                                  --- ---  i,j
//                                         i   j

   {
Sindex=0;					// Keep compiler happy
std::string Llbl = "L";
std::string Jlbl = "J";
//   double cutoff = 1e-12;
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[2*rank+1];
   T2s = new gen_op[2*rank+1];
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
   int m;					// z-component of ang. momentum
   std::string slabs;
   std::string *Wlbls;				// Array of transition labels
   Wlbls = new std::string[hs*hs];
   W_labels(Wlbls, sys, Ho, Windex);		// Fill up the transition labels
   int npairs=ns*ns;				// Total number of spin pairs
   int maxt = 2*(2*rank+1)*(hs+1);		// Max. contributions per dip. pair

   complex *strsJ;			// Array for interaction strengths (J's)
   int *trnsJ;			// Array for transition indices (J's)
   int *contJ;				// Number of contributors per dipole pair
   std::string *Jlbs;			// Array for J labels (J, can use L,K, etc.)
   complex *strsL;			// Array for interaction strengths (L's)
   int *trnsL;			// Array for transition indices (L's)
   int *contL;				// Number of contributors per dipole pair
   std::string *Llbs;			// Array for L labels (L, can use L,Q, etc.)
   std::string *spns;			// Array for spin labels (ijkl, ijjk, etc.)

   strsJ = new complex[npairs*maxt];
   trnsJ = new int[npairs*maxt];	// Array for transition indices (J's)
   contJ = new int[npairs];		// Number of contributors per dipole pair
   Jlbs  = new std::string[npairs];		// Array for J labels (J, can use L,K, etc.)
   strsL = new complex[npairs*maxt];	// Array for interaction strengths (L's)
   trnsL = new int[npairs*maxt];		// Array for transition indices (L's)
   contL = new int[npairs];		// Number of contributors per dipole pair
   Llbs  = new std::string[npairs];		// Array for L labels (L, can use L,Q, etc.)
   spns  = new std::string[npairs];		// Array for spin labels (ijkl, ijjk, etc.)

   int Jterm0=0, Lterm0=0;			// Total terms contributing to <a,a'|R|b,b'>
   int i,j,ij=0;
   std::string *Slbls;				// strings for spin labels
   Slbls = new std::string[ns];
   Spin_labels(Slbls, sys);			// Fill up spin labels
   for(i=0; i<ns; i++)				// Sum over spins i (mu1)
     {
     for(m=-rank; m<=rank; m++)			// Put spin tensor for i (mu1) into a
       {					// vector of operators in basis of Ho
       T1s[m+rank] = T1[i].component(rank,m);
//                             sys.get_basis());
//       T1s[m+rank] = gen_op(T1[i].component(rank,m),
//                             sys.get_basis());
       T1s[m+rank].Op_base(Ho);
       }
     for(j=0; j<ns; j++)			// Sum over spins j
       {
       slabs = Slbls[i] + Slbls[j];
       if(i==j)					// Auto-correlation term
         {					//                      12
         Rel_12(hs,rank,T1s,T1s, 		// Get terms for <a,a'|R  |b,b'>
             Jterms,strsJ,trnsJ,	 	//			ij
                 Lterms,strsL,trnsL,		// in compact numerical format
                              Wlbls,a,aa,b,bb);
         }
       else if(i!=j)				// Cross-correlation term
         {
         for(m=-rank; m<=rank; m++)		// Put spin tensor for j into a
           {					// vector of operators in basis of Ho
           T2s[m+rank] = T2[j].component(rank,m);
//           T2s[m+rank] = gen_op(T2[j].component(rank,m),
//                                 sys.get_basis());
           T2s[m+rank].Op_base(Ho);
           }					//                      12
         Rel_12(hs,rank,T1s,T2s, 		// Get terms for <a,a'|R  |b,b'>
             Jterms,strsJ,trnsJ,	 	//			ij
                      Lterms,strsL,trnsL, 	// in compact numerical format 
                              Wlbls,a,aa,b,bb);
         }
       contJ[ij] = Jterms-Jterm0;		// Number of terms this spin pair
       contL[ij] = Lterms-Lterm0;		// Number of terms this spin pair
       Jlbs[ij] = Jlbl;				// Store J label for these terms
       Llbs[ij] = Llbl;				// Store L label for these terms
       spns[ij] = slabs;			// Store spin label for these terms
       ij++;					// Increment total spin pair count
       Jterm0 = Jterms;				// Update the starting counter
       Lterm0 = Lterms;				// Update the starting counter
       }
     }
   if(Jterms && DFS >= 0)
     Rel(0, Jterms, npairs, strsJ, trnsJ,	// Put J(w) terms into string arrays
           Wlbls, contJ, spns, Jlbs, Mlabel, 	// line[1-3], and track signs
                   line1, line2, line3, signs);
   else
     Jterms = 0;
   if(Lterms && DFS <= 0)
     Rel(Jterms, Lterms, npairs, strsL, trnsL,	// Add L(w) terms to string arrays
           Wlbls, contL, spns, Llbs, Mlabel, 	// line[1-3], track signs also
                   line1, line2, line3, signs);
   else
     Lterms = 0;
   Rel_clean(T1s, T2s, rank);
   delete [] Wlbls;
   delete [] strsJ;
   delete [] trnsJ;
   delete [] contJ;	
   delete [] Jlbs;	
   delete [] strsL;	
   delete [] trnsL;	
   delete [] contL;	
   delete [] Llbs;		
   delete [] spns;
   delete [] Slbls;
   return;
   }


void Rijk_el(const sys_dynamic& sys, gen_op& Ho,
           int rank, spin_T* T1, spin_T* T2, std::string& Mlabel,
        std::string* line1, std::string* line2, std::string* line3, int* signs,
                int& Jterms, int& Lterms,
                              int a, int aa, int b, int bb,
                               int DFS, int Windex, int Sindex)

	// Input 		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			rank  : Rank of the interactions
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			Mlabel: Label for two mechanisms
	//			line1 : strings for the upper line of component
	//			line2 : strings for the middle line of component
	//			line3 : strings for the lower line of component
	//			signs : Flags for the signs of the components
	//			a, aa : 1st transition indices
	//			b, bb : 2nd transition indices
	//			DFS   : Flag for dynamic frequency shifts
	//					> 0  J Terms Only				    
	//					  0  Both J and L Terms (Default)
	//					< 0  L Terms Only
	//			Windex: Flag for frequency indexing
        //                                 >1 = Always use spin indices
        //                                  0 = Use spin indices if needed
        //                                 -1 = Never use spin indices (default)
        //                                 -2 = Assume extreme narrowing
        //                                 -3 = Assume extreme broadening
	//			Sindex: Flag for spin indexing
        //                                  0 = Use only if needed (default)
        //                                  1 = Always use spin indices
        //                                 <1 = Never use spin indices
	// Output		void  :
	// Note			      : This routine assumes that the first
	//				interaction involves two spins and
	//				the second interaction involves 1 spin
	//				e.g. D-SA, D-Q

//	         Mechanism 1 = Spin-Pair; Mechanism 2 = Spin

//                                        --- --- ---
//                            12          \   \   \    12
//                     <a,a'|R  |b,b'> =  /   /   /   R
//	                                  --- --- ---  ij,k
//                                         i   j   k

   {
Sindex=0;					// Keep compiler happy
std::string Jlbl = "J";
std::string Llbl = "L";
//   double cutoff = 1.e-12;
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[2*rank+1];
   T2s = new gen_op[2*rank+1];
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
//   int ls = hs*hs;				// Total system Liouville space
   int m;					// z-component of ang. momentum
   std::string slabs;
   std::string *Wlbls;				// Array of transition labels
   Wlbls = new std::string[hs*hs];
   W_labels(Wlbls, sys, Ho, Windex);		// Fill up the transition labels
   int ndip = sys.dipoles();			// Total number of dipoles
   int npairs = ndip*ns;			// Total number of dipole-spin pairs
   int maxt = 2*(2*rank+1)*(hs+1);		// Max. contributions per dip. pair


   complex *strsJ;			// Array for interaction strengths (J's)
   int *trnsJ;			// Array for transition indices (J's)
   int *contJ;				// Number of contributors per dipole pair
   std::string *Jlbs;			// Array for J labels (J, can use L,K, etc.)
   complex *strsL;			// Array for interaction strengths (L's)
   int *trnsL;			// Array for transition indices (L's)
   int *contL;				// Number of contributors per dipole pair
   std::string *Llbs;			// Array for L labels (L, can use L,Q, etc.)
   std::string *spns;			// Array for spin labels (ijkl, ijjk, etc.)

   strsJ = new complex[npairs*maxt];
   trnsJ = new int[npairs*maxt];	// Array for transition indices (J's)
   contJ = new int[npairs];		// Number of contributors per dipole pair
   Jlbs  = new std::string[npairs];		// Array for J labels (J, can use L,K, etc.)
   strsL = new complex[npairs*maxt];	// Array for interaction strengths (L's)
   trnsL = new int[npairs*maxt];		// Array for transition indices (L's)
   contL = new int[npairs];		// Number of contributors per dipole pair
   Llbs  = new std::string[npairs];		// Array for L labels (L, can use L,Q, etc.)
   spns  = new std::string[npairs];		// Array for spin labels (ijkl, ijjk, etc.)

   int Jterm0=0, Lterm0=0;			// Total terms contributing to <a,a'|R|b,b'>
   int i,j,ij=0,ijk=0;
   std::string *Slbls;				// strings for spin labels
   Slbls = new std::string;
   Spin_labels(Slbls, sys);			// Fill up spin labels
   for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij
     {
     for(j=i+1; j<ns; j++,ij++)
       {
       for(m=-rank; m<=rank; m++)		// Put spin tensor for i,j into a
         {					// vector of operators in basis of Ho
         T1s[m+rank] = T1[ij].component(rank,m);
//         T1s[m+rank] = gen_op(T1[ij].component(rank,m),
//                          sys.get_basis());
         T1s[m+rank].Op_base(Ho);
         }
       for(int k=0; k<ns; k++)			// Sum over spins k
         {
         slabs = Slbls[i] + Slbls[j] + Slbls[k];
         for(m=-rank; m<=rank; m++)		// Put spin tensor for k into a
           {					// vector of operators in basis of Ho
           T2s[m+rank] = T2[k].component(rank,m);
//           T2s[m+rank] = gen_op(T2[k].component(rank,m),
//                                 sys.get_basis());
           T2s[m+rank].Op_base(Ho);
           }					//                      12
         Rel_12(hs,rank,T1s,T2s, 		// Get terms for <a,a'|R   |b,b'>
             Jterms,strsJ,trnsJ, 		//			ijk
                 Lterms,strsL,trnsL,		// in compact numerical format 
                              Wlbls,a,aa,b,bb);
         contJ[ijk] = Jterms-Jterm0;		// Number of terms this spin pair
         contL[ijk] = Lterms-Lterm0;		// Number of terms this spin pair
         Jlbs[ijk] = Jlbl;			// Store J label for these terms
         Llbs[ijk] = Llbl;			// Store L label for these terms
         spns[ijk] = slabs;			// Store spin label for these terms
         ijk++;					// Increment total dipole-spin pair count
         Jterm0 = Jterms;			// Update the starting counter
         Lterm0 = Lterms;			// Update the starting counter
         } 					// Increment spin k (mu2)
       }					// Increment pair ij (mu1)
     }
   if(Jterms && DFS >= 0)
     Rel(0, Jterms, npairs, strsJ, trnsJ,	// Put J(w) terms into string arrays
           Wlbls, contJ, spns, Jlbs, Mlabel, 	// line[1-3], and track signs
                   line1, line2, line3, signs);
   else
     Jterms = 0;
   if(Lterms && DFS <= 0)
     Rel(Jterms, Lterms, npairs, strsL, trnsL,	// Add L(w) terms to string arrays
           Wlbls, contL, spns, Llbs, Mlabel, 	// line[1-3], track signs also
                   line1, line2, line3, signs);
   else
     Lterms = 0;
   Rel_clean(T1s, T2s, rank);
   delete [] Wlbls;
   delete [] strsJ;
   delete [] trnsJ;
   delete [] contJ;	
   delete [] Jlbs;	
   delete [] strsL;	
   delete [] trnsL;	
   delete [] contL;	
   delete [] Llbs;		
   delete [] spns;
   delete [] Slbls;
   return;
   }


void Rkij_el(const sys_dynamic& sys, gen_op& Ho,
           int rank, spin_T* T1, spin_T* T2, std::string& Mlabel,
        std::string* line1, std::string* line2, std::string* line3, int* signs,
                    int& Jterms, int& Lterms,
                              int a, int aa, int b, int bb,
                               int DFS, int Windex, int Sindex)

	// Input 		sys   : Dynamic spin system
	//			Ho    : Static Hamiltonian
	//			rank  : Rank of the interactions
	//			T1    : Spin tensor, mu1
	//			T2    : Spin tensor, mu2
	//			Mlabel: Label for two mechanisms
	//			line1 : strings for the upper line of component
	//			line2 : strings for the middle line of component
	//			line3 : strings for the lower line of component
	//			signs : Flags for the signs of the components
	//			a, aa : 1st transition indices
	//			b, bb : 2nd transition indices
	//			DFS   : Flag for dynamic frequency shifts
	//					> 0  J Terms Only				    
	//					  0  Both J and L Terms (Default)
	//					< 0  L Terms Only
	//			Windex: Flag for frequency indexing
        //                                 >1 = Always use spin indices
        //                                  0 = Use spin indices if needed
        //                                 -1 = Never use spin indices (default)
        //                                 -2 = Assume extreme narrowing
        //                                 -3 = Assume extreme broadening
	//			Sindex: Flag for spin indexing
        //                                  0 = Use only if needed (default)
        //                                  1 = Always use spin indices
        //                                 <1 = Never use spin indices
	// Output		void  :
	// Note			      : This routine assumes that the first
	//				interaction involves one spin and
	//				the second interaction involves 2 spins
	//				e.g. SA-D, Q-D

//	            Mechanism 1 = Spin; Mechanism 2 = Spin-Pair

//                                        --- --- ---
//                            12          \   \   \    12
//                     <a,a'|R  |b,b'> =  /   /   /   R
//	                                  --- --- ---  k,ij
//                                         k   i   j

   {
Sindex=0;					// Keep compiler happy
std::string Jlbl = "J";
std::string Llbl = "L";
//   double cutoff = 1.0e-12;
   gen_op *T1s;					// Compiler didn't like gen_op T1s[5]
   gen_op *T2s;					// These are for spin tensor components
   T1s = new gen_op[2*rank+1];
   T2s = new gen_op[2*rank+1];
//   double w0=0,w1=0,w2=0,wi=0,wj=0;		// Used for J's levels 0 & 1
   int ns = sys.spins();			// Total number of spins
   int hs = sys.HS();				// Total system Hilbert space
//   int ls = hs*hs;				// Total system Liouville space
   int m;					// z-component of ang. momentum
   std::string slabs;
   std::string *Wlbls;
   Wlbls = new std::string[hs*hs];				// Array of transition labels
   W_labels(Wlbls, sys, Ho, Windex);		// Fill up the transition labels
   int ndip = sys.dipoles();			// Total number of spin pairs
   int npairs = ndip*ns;			// Total number of dipole-spin pairs
   int maxt = 2*(2*rank+1)*(hs+1);		// Max. contributions per pair

   complex *strsJ;			// Array for interaction strengths (J's)
   int *trnsJ;			// Array for transition indices (J's)
   int *contJ;				// Number of contributors per dipole pair
   std::string *Jlbs;			// Array for J labels (J, can use L,K, etc.)
   complex *strsL;			// Array for interaction strengths (L's)
   int *trnsL;			// Array for transition indices (L's)
   int *contL;				// Number of contributors per dipole pair
   std::string *Llbs;			// Array for L labels (L, can use L,Q, etc.)
   std::string *spns;			// Array for spin labels (ijkl, ijjk, etc.)

   strsJ = new complex[npairs*maxt];
   trnsJ = new int[npairs*maxt];	// Array for transition indices (J's)
   contJ = new int[npairs];		// Number of contributors per dipole pair
   Jlbs  = new std::string[npairs];		// Array for J labels (J, can use L,K, etc.)
   strsL = new complex[npairs*maxt];	// Array for interaction strengths (L's)
   trnsL = new int[npairs*maxt];		// Array for transition indices (L's)
   contL = new int[npairs];		// Number of contributors per dipole pair
   Llbs  = new std::string[npairs];		// Array for L labels (L, can use L,Q, etc.)
   spns  = new std::string[npairs];		// Array for spin labels (ijkl, ijjk, etc.)

   int Jterm0=0, Lterm0=0;			// Total terms contributing to <a,a'|R|b,b'>
   int i,j,ij,kij=0;
   std::string *Slbls;				// strings for spin labels
   Slbls = new std::string[ns];
   Spin_labels(Slbls, sys);			// Fill up spin labels
   for(int k=0; k<ns; k++)			// Sum over spins k (mu1)
     {
     for(m=-rank; m<=rank; m++)			// Put spin tensor for k (mu1) into a
       {					// vector of operators in basis of Ho
       T1s[m+rank] = T1[k].component(rank,m);
//       T1s[m+rank] = gen_op(T1[k].component(rank,m),
//                               sys.get_basis());
       T1s[m+rank].Op_base(Ho);
       }
     ij = 0;					// Set dipole count to zero
     for(i=0; i<ns-1; i++)			// Sum spin pairs i&j: dipole ij (mu2)
       {
       for(j=i+1; j<ns; j++)
         {
         slabs = Slbls[k] + Slbls[i] + Slbls[j];
         for(m=-rank; m<=rank; m++)		// Put spin tensor for i,j into a
           {					// vector of operators in basis of Ho
           T2s[m+rank] = T2[ij].component(rank,m);
//           T2s[m+rank] = gen_op(T2[ij].
//                          component(rank,m),
//                            sys.get_basis());
             T2s[m+rank].Op_base(Ho);
           }					//                      12
         Rel_12(hs,rank,T1s,T2s, 		// Get terms for <a,a'|R   |b,b'>
             Jterms,strsJ,trnsJ, 		//			kij
                 Lterms,strsL,trnsL,		// in compact numerical format 
                              Wlbls,a,aa,b,bb);
         contJ[kij] = Jterms-Jterm0;		// Number of terms this spin pair
         contL[kij] = Lterms-Lterm0;		// Number of terms this spin pair
         Jlbs[kij] = Jlbl;			// Store J label for these terms
         Llbs[kij] = Llbl;			// Store L label for these terms
         spns[kij] = slabs;			// Store spin label for these terms
         Jterm0 = Jterms;			// Update the starting counter
         Lterm0 = Lterms;			// Update the starting counter
         ij++;
         kij++;					// Increment total dipole-spin pair count
         }					// Increment second dipole (mu2)
       }
     }						// Increment spin (mu1)
   if(Jterms && DFS >= 0)
     Rel(0, Jterms, npairs, strsJ, trnsJ,	// Put J(w) terms into string arrays
           Wlbls, contJ, spns, Jlbs, Mlabel, 	// line[1-3], and track signs
                   line1, line2, line3, signs);
   else
     Jterms = 0;
   if(Lterms && DFS <= 0)
     Rel(Jterms, Lterms, npairs, strsL, trnsL,	// Add L(w) terms to string arrays
           Wlbls, contL, spns, Llbs, Mlabel, 	// line[1-3], track signs also
                   line1, line2, line3, signs);
   else
     Lterms = 0;
   Rel_clean(T1s, T2s, rank);

   delete [] Wlbls;
   delete [] strsJ;
   delete [] trnsJ;
   delete [] contJ;	
   delete [] Jlbs;	
   delete [] strsL;	
   delete [] trnsL;	
   delete [] contL;	
   delete [] Llbs;		
   delete [] spns;
   delete [] Slbls;

   return;
   }

// ______________________________________________________________________
// ************ NUMERICAL RELAXATION MATRIX ELEMENT FUNCTIONS ***********
// ______________________________________________________________________

// These functions fill numerical arrays which represent the terms which
// contribute to a specific relaxation matrix element from a specific
// pair of interaction types (relaxation mechanisms) and a specific set
// of involved spins (e.g. ij-kl, ij-k, i-j).  Terms from both the sym-
// metric and anti-symmetric (dynamic frequency shift) spectral density
// components are treated.

void Rel_12(int hs, int rank, gen_op* T1s, gen_op* T2s,
    int& Jterms, complex* strsJ, int* trnsJ,
      int& Lterms, complex* strsL, int* trnsL,
        std::string* Wlbls, int a, int aa, int b, int bb, double cutoff)

	// Input		hs    : Spin system Hilbert space
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			Jterms: Number of non-zero J(w) terms in summation
	//			strsJ : Array of interaction strengths for J(w)
	//			trnsJ : Array of transition indices for J(w)
	//			Lterms: Number of non-zero L(w) terms in summation
	//			strsL : Array of interaction strengths for L(w)
	//			trnsL : Array of transition indices for L(w)
	//			Wlbls : string array of transition labels
	//			a, aa : 1st transition indices
	//			b, bb : 2nd transition indices
	//			cutoff: Value at which interaction is taken to be zero
	// Output		void  : Arrays strsJ, trnsJ, strsL, and trnsL have
	//				information added to them starting at Jterms,Lterms
	//				The values of Jterms and Lterms is adjusted
	// Note			      : T1s, T2s are assumed in their proper bases
	// Note			      : If J(L)terms != 0 at start, terms are added to arrays
	// Note			      : The subscripts 1 and 2 are cryptic in this context
	// 			        as they are often used simultaneously for indexing
	//				both a mechanism and the spin(s) involved in relaxation
	//
	//					    Mechanism	  Rank	1    2
	//					----------------- ---- ---- ----
	//					 Dipolar	   2    ij   kl
	//					 Shift Anisotropy  2    i    j
	//					 Quadrupolar	   2    i    j
	//					 Random Field      1    i    j

/* The purpose of this routine is determine all non-zero contributions to the relaxation
   matrix element <a,a'|R|b,b'> for the interactions 1 with 2 as specified by the two
   input spin tensors (via their components), T1s and T2s.  The contribution to <a,a'|R|b,b'>
   is given by

                      rank    ls
   <a,a'|R   |b,b'>   ---   [ ---
          1,2         \     | \                 m       m 
   ---------------- = /     | /   delta     <a|T |g><b|T |g> J  (w  )
  	Xi            ---   | ---      a',b'    1       2    ~12  gb
          1,2	    m=-rank [  g
  
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
                1,2	     i=0                             j=0
  
   where trans[i] relates to the original eigenstates according to
  
                              w           w        trnsJ[i] = a*hs+a' 
                               trnsJ[i] =  aa'           a' = tr%hs;
                                                         a  = (tr-a')/hs;
  
   There will be (2*rank+1)*[ls + 2 + ls] terms in the summation over m & g, however the spin tensor
   products act a a good filter so that most terms will be zero.  Any terms in which the spin tensor
   product has a magnitude below cutoff are herein considered zero and will not be included in the
   ouput arrays. Furthermore, many terms will add if they have the same spectral density functions.
   The end effect of this routine then, is the smallest possible amount of information needed to
   characterize the relaxation the relaxation matrix element for the pair of interactions (1,2) and
   the specified set of spins in terms of scaled spectral density functions.                      */

  {
cutoff=0;					// Keep compiler happy
  int Jtermi = Jterms;				// There are already this many J(w) terms
  int Ltermi = Lterms;				// There are already this many L(w) terms
  Rel_12(hs,rank,T1s,T2s,			// Get terms for <a,a'|R  |b,b'> for J(w)
                 Jterms,strsJ,trnsJ,a,aa,b,bb);	//			12
  int trm;
  for(trm=Jtermi; trm<Jterms; trm++,Lterms++)	// Now copy these into arrays for the L(w)
    {						// terms.  They are the same at this point
    strsL[Lterms] = strsJ[trm];			// although J is symmetric but L is asymmetric
    trnsL[Lterms] = trnsJ[trm];
    }
  Rel_12_condense(hs, Jtermi, Jterms,		// Place all the terms so that have a postive
                         strsJ, trnsJ, Wlbls);	// w in J(w) by use of symmetry J(-w) = J(w)
  Rel_12_condense(hs, Jtermi, 			// Add together terms having the same J(w)
                        Jterms, strsJ, trnsJ);
  Rel_12_condense(hs, Ltermi, Lterms,		// Place all the terms so that they have a
                     strsL, trnsL, Wlbls, 1);	// positive w in L(w) by use of symmetry L(-w) = -L(w)
  Rel_12_condense(hs, Ltermi,			// Add together terms having the same L(w)
                    Lterms, strsL, trnsL, 1);
  return;
  }


  void Rel_12(int hs, int rank, gen_op* T1s, gen_op* T2s,
                    int& nterms, complex* strs, int* trns, 
		 	  int a, int aa, int b, int bb, double cutoff)

	// Input		hs    : Spin system Hilbert space
	//			rank  : Rank of the two interactions
	// 			T1s   : Spin tensor components, 1st interaction
	// 			T2s   : Spin tensor components, 2nd interaction
	//			strs  : Array of interaction strengths
	//			trns  : Array of transition indices
	//			a, aa  : 1st transition indices
	//			b, bb : 2nd transition indices
	//			cutoff: Value at which interaction is taken to be zero
	// Output		void  : Arrays strs and trns are filled, nterms value set
	// Note			      : T1s, T2s are assumed in their proper bases
	// Note			      : If nterms != 0 at start, terms are added to arrays
	// Note			      : The subscripts 1 and 2 are cryptic in this context
	// 			        as they are often used simultaneously for indexing
	//				both a mechanism and the spin(s) involved in relaxation
	//
	//					    Mechanism	  Rank	1    2
	//					----------------- ---- ---- ----
	//					 Dipolar	   2    ij   kl
	//					 Shift Anisotropy  2    i    j
	//					 Quadrupolar	   2    i    j
	//					 Random Field      1    i    j

/* The purpose of this routine is determine all non-zero contributions to the relaxation
   matrix element <a,a'|R|b,b'> for the interactions 1 with 2 as specified by the two
   input spin tensors (via their components), T1s and T2s.  The contribution to <a,a'|R|b,b'>
   is given by
                      rank    ls
   <a,a'|R   |b,b'>   ---   [ ---
          1,2         \     | \                 m       m 
   ---------------- = /     | /   delta     <a|T |g><b|T |g> J  (w  )
  	Xi            ---   | ---      a',b'    1       2     12  gb
          1,2	    m=-rank [  g
  
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
                                  1,2	      i=0
  
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
   to terms for the same element from other spins and other mechanisms.                       */
 
   {
//   int ntermsi = nterms;			// There are this many terms already
   int k=0, g=0;				// Indices for m-component & wavefunction 
   int tr = 0;					// This is the active transition
   complex str;					// This is the interaction strength
   for(int m=-rank; m<=rank; m++)
     { 						//                     m         m
     str = T1s[k].get(a,b)*T2s[k].get(aa,bb); 	// RII: -J(bb-aa)*<a|T1 |b><aa|T2 |bb>
     if(norm(str) > cutoff)
       {
       tr = bb*hs+aa;				//	Determine the transition index
       strs[nterms] = -str;			//	Store the coupling strength
       trns[nterms] = tr;			//	Store the J(W) transition label
       nterms++;				//	Increment the count of non-zero terms
       } 					//                     m         m
     str = T1s[k].get(bb,aa)*T2s[k].get(b,a); 	// RIII: -J(a-b)*<bb|T1 |aa><b|T2 |a>
     if(norm(str) > cutoff)
       {
       tr = a*hs+b;				//	Determine the transition index
       strs[nterms] = -str;			//	Store the coupling strength
       trns[nterms] = tr;			//	Store the J(W) transition label
       nterms++;				//	Increment the count of non-zero terms
       }
     for(g=0; g<hs; g++)
       {
       if(aa == bb)
         { 					//                    m        m
	 str = T1s[k].get(a,g)*T2s[k].get(b,g); // RI:    J(g-b)*<a|T1 |g><b|T2 |g>
         if(norm(str) > cutoff)
           {
           tr = g*hs+b;				//	Determine the transition index
           strs[nterms] = str;			//	Store the coupling strength
           trns[nterms] = tr;			//	Store the J(W) transition label
           nterms++;				//	Increment the count of non-zero terms
           }
         }
       if(a == b)
         { 					//                     m         m
         str=T1s[k].get(g,aa)*T2s[k].get(g,bb); // RIV:   J(bb-g)*<g|T1 |aa><g|T2 |bb>
         if(norm(str) > cutoff)
           {
           tr = bb*hs+g;			//	Determine the transition index
           strs[nterms] = str;			//	Store the coupling strength
           trns[nterms] = tr;			//	Store the J(W) transition label
           nterms++;				//	Increment the count of non-zero terms
           }
         }
       }
     k++;
     }
   return;
   }


void Rel_12_condense(int hs, int ntermi, int& nterms,
                      complex* strs, int* trns, int anti, double cutoff)

	// Input		hs    : Spin system Hilbert space
	//			ntermi: Initial index to start with
	//			nterms: Number of terms contributing to R element
	//			strs  : Array of interaction strengths
	//			trns  : Array of transition indices
	//			anti  : Flag for J(w) versus L(w)
	//				     0 : Assume J(-w) = J(w)  symmetric
	//				    !0 : Assume L(-w) = -L(w) anti-symmetric
	//			cutoff: Value at which interaction is taken to be zero
	// Output		void  : Arrays strs and trns are modified, nterms value reset
	//				as terms are combined with the same J(w) values
	// Note			      : IT IS ASSUMED HEREIN THAT ALL COMPONENTS FROM
	//				[ntermsi, nterms) STEM FROM THE SAME RELAXATION
	//				MECHANISMS AND INVOLVE THE SAME SPINS
	// Note			      : This condenses according to a transition number
	//				index only.  Transitions with the same frequency
	//				but of a different index are not added herein
	//				(for example WH1 will not combine with WH2)

//     The intent of this routine is to reduce the number of explicit terms contributing
//  to a particular relaxtion matrix element. For example, the following will be performed.
//
//       DD                 DD              DD                     DD                 DD
// -0.5 J    (WC-WH) + 0.2 J    (WH) + 0.3 J    (-WC+WH) --> -0.2 J    (WC-WH) + 0.2 J    (WH)
//       CHCH               CHCH            CHCH                   CHCH               CHHH
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

  {
   int a,b;				// Use for basis function indices
   int same, newnum=0;			// Flag for unique terms, term counter 
   int i,j,tr,trinv;
   for(i=ntermi; i<nterms; i++)		// Loop over all the contributing terms
     {
     tr = trns[i];			// Get the ith component's transition
     same = 0;				// Assume it is unique
     for(j=0; j<newnum && !same; j++)	// Check for match with previous one
       {
       if(tr == trns[j+ntermi])		// If it matches a previous transition
         {				// then flag that its the same and
         same = 1;			// just add it to the previous one
         strs[j+ntermi] += strs[i]; 
         }
       else				// Also check if it is the inverse of
         {				// a previous transition.  In this case
         a = tr%hs;			// it is added dependent upon the 
         b = (tr-a)/hs;			// symmetry of the spectral density
         trinv = a*hs+b;		// function
         if(trinv == trns[j+ntermi])
           {
           same = 1;
           if(anti)
             strs[j+ntermi] -= strs[i];
           else
             strs[j+ntermi] += strs[i];
           }
         }
       }
     if(!same)				// If unique, increment count of unique
       { 				// components of the element
       trns[newnum+ntermi] = tr;
       strs[newnum+ntermi] = strs[i];
       newnum++;
       }
     }
   nterms = ntermi+newnum;		// Now there are only this many terms

//    During Combination of Terms With The Same J(w) in the Previous Section
//    Some Terms May Have Cancelled Each Other Out.  Here, Contributions With
//                  No Appreciable Intensity Are Removed.

   for(i=ntermi; i<nterms; i++)		// Loop over all the current terms
     {
     if(norm(strs[i]) < cutoff)		// See if the intensity is now zero
       {				// If so, them mover all the following
       for(j=i+1; j<nterms; j++)	// terms back in the array one spot.
         {
         strs[j-1] = strs[j];
         trns[j-1] = trns[j];
         }
       nterms--;			// Reduce the number of terms by 1
       i--;				// Recheck same element, its different now
       }
     }
   return;
   }


void Rel_12_condense(int hs, int ntermi, int& nterms, complex* strs,
                                               int* trns, std::string* Wlbls, int anti)

	// Input		hs    : Spin system Hilbert space
	//			ntermi: Initial index to start with
	//			nterms: Number of non-zero terms contributing
	//				to R element
	//			strs  : Array of interaction strengths
	//			trns  : Array of transition indices
	//			anti  : Flag for J(w) versus L(w)
	//				     0 : Assume J(-w) = J(w)  symmetric
	//				    !0 : Assume L(-w) = -L(w) anti-symmetric
	// Output		void  : Arrays strs and trns are modified as
	//				all terms are set with positive w in J(w)
	// Note			      : IT IS ASSUMED HEREIN THAT ALL COMPONENTS FROM
	//				[ntermsi, nterms) STEM FROM THE SAME RELAXATION
	//				MECHANISMS AND INVOLVE THE SAME SPINS
	// Note			      : This may modifiy the transition index based on
	//				the transition label in Wlbls.  All transitions
	//				are set so that the first frequency is positive.
	//				Then all with the same label are set to have the
	//				same transition index.

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
// or anti-symmetric, such as in the case of the terms for dynamic frequency shifts.
//
//                                DD                DD
//                           c * L    (-WH) = -c * L    (WH)
//                                CHCH              CHCH
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

  {
  std::string trlab, neg = "-";
  int a, b, tr, trinv;
  int i;
  for(i=ntermi; i<nterms; i++)		// Loop over all the contributing terms
    {
    tr = trns[i];			// Get the ith component's transition
    trlab = Wlbls[tr];			// Get the transition label 
    if(!trlab.find("-"))		// See if the transition is negative
      {					// If so, set it to the opposite transition
      a = tr%hs;
      b = (tr-a)/hs;
      trinv = a*hs+b;
      trns[i] = trinv;
      if(anti)
        strs[i] *= -1.0;
      }
    }

//    Now, Look For Terms Which Have The Same Transition Labels.  If Two Have
//     The Same Transition Label Then Set Thier Transition Indices The Same

  int j, same;
  for(i=ntermi; i<nterms; i++)		// Loop over all the contributing terms
    {
    same = 0;
    tr = trns[i];			// Get the ith component's transition
    trlab = Wlbls[tr];			// Get the transition label 
    for(j=ntermi; j<nterms && !same; j++)
      if(trlab == Wlbls[trns[j]])
        {
        trns[i] = trns[j];
        same = 1;
        }
    }
  return; 
  }

// ______________________________________________________________________
// ************* STRING RELAXATION MATRIX ELEMENT FUNCTIONS *************
// ______________________________________________________________________


  void Rel(int ntermi, int& nterms, int npairs, complex* strs, int* trns,
             std::string* Wlbls, int* cont, std::string* spns, std::string* Jlbs,
               std::string& Mlabel, std::string* line1, std::string* line2, std::string* line3, int* signs)

	// Input 		ntermi: Number of terms currently in line1, line2, line3, signs
	//      		nterms: Number of terms contributing to <a,aa|R|b,bb>
	//			npairs: Number spin parings in these terms
	//			strs  : Array for interaction strengths for these terms
	//			trns  : Array for transition indices for these terms
	//			Wlbls : Array of transition labels for the system
	//			cont  : Array of contributions per spin pairing
	//			spns  : Array of spin labels per spin pairing
	//			Jlbs  : Array of J labels per 1,2 pair
	//			Mlabel: string label fo two mechanisms
	//			cutoff: Value at which interaction is taken to be zero
	// Output		void  : The contribution to <a,a'|R|b,b'> from the
	//				mechanisms Mlabel and spins spns are added to the arrays
	//				line1, line2, line3, and signs starting at index ntermi.
	//			        The number of terms in these arrays, nterms, is updated
	// Note			      : ALL TERMS ARE ASSUMED FROM THE SAME RELAXATION MECHANISM

//  The purpose of this routine is to fill string arrays with a formal (non-numerical)
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
//	        1,2            ---                    spns[j]           ---           line3[i]
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

    {

//int intlen = 4;					// Set intensity length
//int intdec = 2;					// Set intensity decimals
    int Mlen = Mlabel.length();				// Length of mechanism label
    std::string Ilabel, Jlabel, Wlabel;			// strings for I,J,W labels 
    int Ilen, Jlen, Wlen;				// Lengths of I,J,W labels
    std::string Mlab;					// string for spin labels
    std::string Slab;					// string for spin labels
    int Slen;						// Length of spin labels
    int dellen;						// Difference Slen-Mlen

//			 Here We Loop Over All Interaction Pairs 1,2
//       For Each Pair 1,2 There May Be Many Terms Contributing to <a,a'|R  |b,b'>
//	                  This Number is Stored in the Array cont         12

    int p12term = 0;
    int term = 0;
    double inten;
    for(int p12=0; p12<npairs && term<nterms; p12++)	// Loop over all the 1,2 pairs
      {							// (e.g. ij-kl, i-j, ij-k, k-ij)
      Jlabel = std::string(" ") + Jlbs[p12];			// Build the label for J: e.g. " J"
      Jlen = Jlabel.length();
      Slab = spns[p12];					// Build the label for spins: e.g. "CHNH"
      Slen = Slab.length();				// Determine whether mechanism (e.g. DSA)
      Mlab = Mlabel;					// Copy the mechanism label
      dellen = Slen-Mlen;				// label or spins lable is longer.  Then
      if(dellen>0)					// set it so both are the same length
        Mlab += std::string(dellen, ' ');			// so that they will align easier in string
      else if(dellen)					// format
        {
        dellen = abs(dellen);
        Slab += std::string(dellen, ' ');
        Slen += dellen;
        }
      for(p12term=0; p12term<cont[p12]; p12term++)	// Now loop over all the terms <a,a'|R|b,b'>
        {						//                                    12
        inten = Re(strs[term]);
//        Ilabel = std::string(dtoa(fabs(inten),'f',4,2));
        Ilabel = Gform("%4.2f", fabs(inten));		// Build a label for the intensity
        Ilen = Ilabel.length();
        Wlabel = std::string("(") + Wlbls[trns[term]]	// Build a label for the transition
                 + std::string(")");
        Wlen = Wlabel.length();
        if(inten < 0)					// Now store:	1.) The sign of the intensity 
          signs[ntermi+term] = 0; 			//		2.) A string for the mechanism (1st line)
        else						//		3.) A string for the J(w)      (2nd line)
          signs[ntermi+term] = 1; 			//		4.) A string for the spins     (3rd line)

        line1[ntermi+term] = std::string(Ilen+Jlen, ' ') + Mlab              + std::string(Wlen, ' ');	//      DSA
        line2[ntermi+term] =        Ilabel + Jlabel + std::string(Slen, ' ') + Wlabel;		// 0.5 J   (2WH)
        line3[ntermi+term] = std::string(Ilen+Jlen, ' ') + Slab              + std::string(Wlen, ' ');	//      CCH
        term++;
        }
      }
    return;
    }

// ______________________________________________________________________
// *************** RELAXATION ANALYSIS AUXILIARY FUNCTIONS **************
// ______________________________________________________________________


  void Spin_labels(std::string* Lbls, const spin_sys& sys, int index)

	// Input 		Lbls  : string a for spin labels
	//				Dimension is number of spins
	//			sys   : A spin system
	//			index : Flag for inclusion of spin index
	//				    0 = Use only if needed
	//				    1 = Always use spin indices
	//				   <1 = Never use spin indices
	// Output		void  : Fills the array of spin labels


// The intent of this routine is to                     DD
// construct labels for the spins which can       -0.5 J    (WC-WH)
// be used to produce subscripts in J labels -------->  CHCH  ^  ^
// as well as transition spin labels _________________________|__|

// In this routine, each spin of the system is assigned a unique spin label.
//  This will be as small as possible yet still maintiain it's uniqueness. 
// In the worst case, the label will be the spin's symbol and a spin index.
// For example: 14N0, 14N1, 14N2.  However, if there are no other isotopes of a
// particular spin type, the isotope label will be neglected and 14N0 will
// just be N0.  Next, if there are no other spins of a particular isotope type
// the spin index will be neglected and only N is used.  If the user is familiar
// with the spin system then there can be no question of the spin by the returned
// index no matter what their ordering.  One other tidbit, rather that use 1H
// 2H, and 3H for the hydrogen isotopes, here they are switched to H,D,T.

   {
   int ns = sys.spins();			// Number of spins in the system
   std::string *symbs;				// Allocate an array for spin symbols
   std::string *sanac;				// Array for element types (H,C,N,...)
   int *snumb;				// Array for isotope number(1,2,3,14,15)
   int *iso_unique;				// Array of isotope uniqueness flags
   int *spn_unique;				// Array of spin uniqueness flags

   symbs = new std::string[ns];				// Allocate an array for spin symbols
   sanac = new std::string[ns];				// Array for element types (H,C,N,...)
   snumb = new int[ns];				// Array for isotope number(1,2,3,14,15)
   iso_unique = new int[ns];				// Array of isotope uniqueness flags
   spn_unique = new int[ns];				// Array of spin uniqueness flags

   int i,j;
   for(i=0; i<ns; i++)				// Loop over all the spins, obtain the
     { 						// spin symbols and then parse them
     symbs[i] = sys.symbol(i);			//	This is the symbol    (14N)
     snumb[i] = (sys.isotope(i)).number();	//	This is the isotope # (14)
     sanac[i] = (sys.isotope(i)).element();	//	Rest is the element   (N)
     if(symbs[i] == "2H")			//	Switch 2H to 2D for deuterium
       {
       symbs[i] = std::string("2D");
       sanac[i] = std::string("D");
       }
     if(symbs[i] == "3H")			//	Switch 3H to 3T for tritium
       {
       symbs[i] = std::string("3T");
       sanac[i] = std::string("T");
       }
     }
   for(i=0; i<ns; i++)				// Now loop over all the spins and
     { 						// look for other spins of this type
     spn_unique[i]=1;				//	Assume there aren't any (no indexing)
     if(index>0) spn_unique[i]=0;		//	If desired, force spin indexing
     if(index==0)				//	If use only when necessary, then
       {					//	we must loop over all other spins
       for(j=0; j<ns && spn_unique[i]; j++)	//	and compare their types.
         if(i!=j && symbs[i]==symbs[j])		//	Set it unique flag to 0 if others
           spn_unique[i] = 0;
       }
     }
   for(i=0; i<ns; i++)				// Now loop over all the spins and
     { 						// look for other isotopes of this type
     iso_unique[i]=1;				//	Assume there aren't any
     for(j=0; j<ns && iso_unique[i]; j++)	//	Compare with all other spins
       if(i!=j && sanac[i]==sanac[j]		//	Set it unique flag to 0 if others
                      && snumb[i]!=snumb[j])
         iso_unique[i] = 0;
     }
   for(i=0; i<ns; i++)				// Now loop over all the spins and
     {						// assemble the unique spin label
     Lbls[i] = "";				//	Start with no label
     if(!iso_unique[i])				//	Add isotope number if other
       Lbls[i] += std::string(Gdec(snumb[i]));	//	isotopes of this atom type
     Lbls[i] += sanac[i];			//	Add the element designation
     if(!spn_unique[i])				//	Add the spin number if there
       Lbls[i] += std::string(Gdec(i));		//	are other spins of this type
     }

   delete [] symbs;				// Allocate an array for spin symbols
   delete [] sanac;				// Array for element types (H,C,N,...)
   delete [] snumb;				// Array for isotope number(1,2,3,14,15)
   delete [] iso_unique;				// Array of isotope uniqueness flags
   delete [] spn_unique;				// Array of spin uniqueness flags

   return;
   }


  void W_labels(std::string* Wlabels, const spin_sys& sys, gen_op &Op, int index)

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

  {
  if(index < -3) index = -1;			// Insure no odd index values input
  matrix B((Op.get_basis()).U());		// Get the eigenbasis matrix
  int hs = B.rows();				// Basis size, spin Hilbert space
  if(!hs) return;				// Exit if Null Operator
//  int ls = hs*hs;				// Liouville space, # of transitions
  matrix sps = sys.qStates();			// System pdt basis function spin states
  int *pbf;					// Find the largest product basis function
  pbf = new int[hs];
  double maxcoeff; 				// contributor to each of the wave functions
  int wf, bf;					// This is equivalent to finding the largest
  for(wf=0; wf<hs; wf++)			// number in each column of the basis B
    {
    maxcoeff = 0;
    pbf[wf] = 0;
    for(bf=0; bf<hs; bf++)
    if(norm(B(bf,wf)) > maxcoeff)
      {
      maxcoeff = norm(B(bf,wf));
      pbf[wf] = bf;				// pbf[wf] is the index of the product basis
      } 					// function that dominates wavefunction wf
    }
  int ns = sys.spins();				// Get the number of spins in the system
  int *flip;					// This will store indices of flipped spins
  double *trqn;				// This will store delta Iz of flipped spins
  std::string *Slbls;				// This will store the spin labels used in W
 
  flip = new int[ns];					// This will store indices of flipped spins
  trqn = new double[ns];				// This will store delta Iz of flipped spins
  Slbls = new std::string[ns];				// This will store the spin labels used in W

  Spin_labels(Slbls, sys, index);		// Construct the spin labels
  int flips, same;
  int cnt, tr=0;
  int i, j;
  std::string spin0, spin1, *spins;
  spins = new std::string[ns];
  std::string label;					// Temporary label for the transition
  int wfi, wfj;					// Transiiton Initial & final wavefunctions
  int pbi, pbj;					// Corresponding dominant basis functions
  double qtot = 0;				// Used for extremene broadening (sosi)
  for(wfi=0; wfi<hs; wfi++)
    {
    pbi = pbf[wfi];				// Largest pdt basis function of wfi
    for(wfj=0; wfj<hs; wfj++)
      {
      pbj = pbf[wfj];				// Largest pdt basis function of wfj
      flips = 0;				// Start with no spin flips in transition
      for(i=0; i<ns; i++)			// Loop through spins and find the 
        if(sps(pbi,i) != sps(pbj,i)) 		// ones that have flipped by comparing
          { 					// the spin states.  If different, store
          flip[flips] = i;			// the spins index as well as the del Iz
          trqn[flips]=Re(sps(pbi,i)-sps(pbj,i));// involved in the flip.  Then keep count
          flips++;				// of the number of spins flipped
          }
      if(!flips)				// If no spins flip, set the label to 0
        Wlabels[tr] = std::string("0");
      else if(flips == 1)			// If one just 1 spin flips, its an SQT
        {
        spins[0] = Slbls[flip[0]];		//	Get the label of the spin flipped
        label = "";				//	Start with no label
        if(trqn[0] < 0) label += "-";		//	Explicitly set the sign of the transition
        if(abs(int(trqn[0])) != 1)		//	If the del Iz is 1, don't write it
          label += Gdec(abs(int(trqn[0])));	//	i.e. write 3WH but not 1WH, thats just WH
        label += std::string("W") + spins[0];	//	Now add in the frequency part WX
        Wlabels[tr] = label;			//	Set this temporary label to the stored label
        if(index == -3) Wlabels[tr] = "X";	// Crude enforcement of extreme broadening (sosi)
        }
      else if(flips == 2)			// If two spins flipped, its either a DQT or a ZQT
        {
        spins[0] = Slbls[flip[0]];		//	Get the label of the 1st spin flipped
        spins[1] = Slbls[flip[0]];		//	Get the label of the 2nd spin flipped
        label = "";				//	Start with no label
        if(spins[0] != spins[1])		//	If different spin types, must use both
          {					//	spins in the label.
          if(trqn[0] < 0) label += "-";		//	Have explicit sign on 1st only if negative
          if(abs(int(trqn[0])) != 1)		//	If the del Iz is 1, don't write the value
            label += Gdec(abs(int(trqn[0]))); 	//	i.e. write 3WH but not 1WH, thats just WH
          label += std::string("W") + spins[0]; 	//	Now add in the frequency part WX of 1st spin
          if(trqn[1] < 0) label += " - ";	//	For the 2nd spin, use the sign as a separator
          else            label += " + ";
          if(abs(int(trqn[1])) != 1)		//	Still, if del Iz is 1 for this flip, dont
            label += Gdec(abs(int(trqn[1])));	//	bother to write it out.
          label += std::string("W") + spins[1];	//	Now add in the frequency part WX of 2nd spin
          Wlabels[tr] = label;
          if(index == -3 && trqn[0]+trqn[1]!=0) // Crude enforcement of extreme broadening (sosi)
                             Wlabels[tr] = "X";
          }
        else					//	If two spins of the same type flipped,
          {					//	just add them together.  They may cancel
          trqn[0] += trqn[1];			//	if homonuclear ZQC.  If so just set the
          if(int(trqn[0]) == 0)			//	label to 0
            Wlabels[tr] = std::string("0");
          else 					// 	Or they may add if homonuclear DQC.
            {					//	If so then start the label with a - if
            if(trqn[0] < 0) label += "-";	//	negative.  
            if(abs(int(trqn[0])) != 1) 		//	If the del Iz is 1, don't write it
              label += Gdec(abs(int(trqn[0]))); 	//	i.e. write 3WH but not 1WH, thats just WH
            label += std::string("W") + spins[0]; 	//	Now add in the frequency part WX
            Wlabels[tr] = label;		//	Set this temporary label to the stored label
            if(index == -3) Wlabels[tr] = "X";	// Crude enforcement of extreme broadening (sosi)
            }
          }
        }
      else					// If more that two spins have flipped, this
        { 					// is either a MQT with |M|>2 or a lower QT
        label = "";				// that has combination flips.  Start with no label.
        cnt = 0;				// Keep count of unique spin types involved
        for(i=0; i<flips; i++)			// Loop over all spins that are flipped
          {
          spin0 = Slbls[flip[i]];		// This is the label for the flipped spin i
          same = 0;				// Assume its label is different from all other spins
          for(j=0; j<i && !same; j++)		// Loop over all previous flipped spins
            {					// and insure that this one hasn't been
            spin1 = Slbls[flip[j]];		// accounted for yet
            if(spin0 == spin1)			// If it has been already accounted for, flag it
              same = 1;
            }
          if(!same)				// This will occur only if spin i is the first
            {					// of its type to be added to the label
            for(j=i+1; j<flips; j++)		// Loop over all next flipped spins
              {					// If they are the same type then just
              spin1 = Slbls[flip[j]];		// accounted for them here
              if(spin0 == spin1)
                trqn[i] += trqn[j];
              }
            if(trqn[i] != 0)			// Now add this spin type contribution to
              {					// the label (if it is not zero)
              qtot += trqn[i];			// For extreme broadening (sosi)
              if(!cnt)
                {
                if(trqn[i] < 0) label += "-";	// For the first type in the label, write
                }				// sign only if negative.  For subsequent
              else				// spin types, use the sign as a spacer
                {
                if(trqn[i] < 0) label += " - ";
                else            label += " + ";
                }
              if(abs(int(trqn[i])) != 1)	// If the spin type has a del Iz of 1
                label += Gdec(abs(int(trqn[i])));// don't bother putting the 1 in the label
                label += std::string("W") + spin0; 	// Now add the WX part, wher X is spin label
              cnt++;
              }
            }
          }
        Wlabels[tr] = label;			// Now set the temp label to the stored one
        if(index == -3 && qtot) 		// Crude enforcement of extreme broadening (sosi)
                         Wlabels[tr] = "X";
        }
      if(index == -2) Wlabels[tr] = "0";	// Crude enforcement of extreme narrowing (sosi)
      qtot = 0;					// Used for extremene broadening (sosi)
      tr++;
      }
    }
  delete [] pbf;
  delete [] flip;					// This will store indices of flipped spins
  delete [] trqn;				// This will store delta Iz of flipped spins
  delete [] Slbls;				// This will store the spin labels used in W
  delete [] spins;

  return;
  }


  void Elem_labels(std::string* Lbls, std::string& R, std::string& M, std::string& S,
            int a, int aa, int b, int bb, int la, int laa, int lb, int lbb)

	// Input 	       Lbls   : std::strings for relaxation matrix
	//			 	element labels
	//	               R      : Relaxation matrix label
	//	               M      : Mechanisms label
	//		       S      : Spins label
        //                     a, aa  : 1st transition indices
        //                     b, bb  : 2nd transition indices
        //                     la,laa : Desired string length of the
	//			        1st transition indices
        //                     lb,lbb : Desired string length of the
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

   {
   int Slen, Mlen, dellen;			// For string lengths
   if(R == std::string(R.length(), ' '))		// If blank R then just make all
     {						// three lines blanks.  This is handy
     Lbls[0] = R;				// For adding strings so element terms
     Lbls[1] = R;
     Lbls[2] = R;
     }
   else
     {
     Lbls[1] = std::string("<");			// Start of R:	<
     if(la) Lbls[1] += std::string(Gdec(a,la));	//		<a
     else   Lbls[1] += std::string(Gdec(a));
     Lbls[1] += std::string(",");			//		<a,
     if(laa) Lbls[1] += std::string(Gdec(aa,laa));	//		<a,aa
     else    Lbls[1] += std::string(Gdec(aa));
     Lbls[1] += std::string("|") + R;		//		<a,aa|R
     std::string end = std::string("|");			//		|
     if(lb) end += std::string(Gdec(b,lb));		//		|b
     else   end += std::string(Gdec(b));
     end += std::string(",");			//		|b,
     if(lbb) end += std::string(Gdec(bb,lbb));	//		|b,bb
     else    end += std::string(Gdec(bb));
     end += std::string(">");			//		|b,bb>
     Slen = S.length();				// Get the spins label length
     Mlen = M.length();				// Get the mechanisms label length
     dellen = Slen-Mlen;			// See which is longest
     if(dellen>0)				// Set it so both are the same length
       M += std::string(dellen, ' ');		// so that they will align easier in
     else if(dellen)				// string manipulations of the element
       {					// Both will be set to be length of
       dellen = abs(dellen);			// Slen at the end
       S += std::string(dellen, ' ');
       Slen += dellen;
       }
     Lbls[0] = std::string(Lbls[1].length(), ' ')   + M                 + std::string(end.length(), ' ');
     Lbls[2] = std::string(Lbls[1].length(), ' ')   + S                 + std::string(end.length(), ' ');
     Lbls[1] = Lbls[1]                         + std::string(Slen, ' ') + end;
     }
   return;
   }


  void Rel_clean(gen_op* T1s, gen_op* T2s, int rank)

	// Input 		T1s   : Array of general operators
	//			T2s   : Array of general operators
	//			rank  : Rank of a tensor
	// Output		void  : The arrays T1s and T2s contain
	//				the components of two tensors
	//				of rank "rank".  These are zeroed.
	// Note			      : This routine is necessary if 
	//				memory is allocated by use of the
	//				"new" function in setting up the
	//				two tensor component arrays.

   {
   gen_op Op;					// Use a NULL operator
   int cmps = 2*rank+1;				// There are this many components
   int i;
   if(T1s)					// If T1 tensor components exist
     {						// then loop through and set them
     for(i=0; i<cmps; i++)			// all to NULL operators
       T1s[i] = Op;
     T1s = NULL;
     }
   if(T2s)					// If T2 tensor components exist
     {						// then loop through and set them
     for(i=0; i<cmps; i++)			// all to NULL operators
       T2s[i] = Op;
     T2s = NULL;
     }
   return;
   }

// ______________________________________________________________________
// ********** RELAXATION MATRIX OUTPUT ELEMENT STRING FUNCTIONS *********
// ______________________________________________________________________


  void Rel(std::ostream& ostr, int nterms, std::string* line1, std::string* line2,
	   std::string* line3, int* signs, std::string* Elabel, int add, int ncols)

	// Input		ostr  : An output stream
	//      		nterms: Number of terms contributing to <a,aa|R|b,bb>
	//			line1 : strings for the upper line of the printed element
	//			line2 : strings for the middle line of the printed element
	//			line3 : strings for the lower line of the printed element
	//			signs : Flags for the signs of the components
	//			Elabel: Relaxation matrix element label
	//			add   : Flag for whether to add lines or not
	//				       add = 0  -> Print Element Label
	//				       add >= 0 -> Don't Print Label, Its Been Done
	//						   Label had length of add
	//			cols  : Number of columns to print out
	// Output		void  : The contribution to <a,a'|R|b,b'> from the
	//				mechanisms ss and spins sl are sent into
	//				the output stream ostr.
	// Note			      : T1s, T2s are assumed in their proper bases

// In this routine, the nterms components that contribute to a particular
// relaxation matrix element are sent into the output stream ostr.  The element
// can be over many relaxation mechanisms, one pair of mechanisms, or one pair of
// mechanism in which only particular spins are involved.  This is determined by
// the formation of the input string arrays which are generated with other routines
// in this module.  Also, the number of terms, or equivalently, how detailed the
// element is characterized is also determined by how the string arrays are constructed.

//                   --- [ ---
//        ss         \   | \                 m       m 
// <a,a'|R  |b,b'> = /   | /   delta     <a|T |g><b|T |g> J  (w  )
//	  sl         --- | ---      a',b'    1       2     sl  gb
//		      m  [  g
//
//                               m        m                       m        m
//                         - <a|T |b><a'|T |b'> J  (w    ) - <b'|T |a'><b|T |a> J  (w  )
//	     	                 1        2      sl  b'a'         1        2     sl  ab

//
//                                                   ---                                      ]
//                                                   \               m        m               |
//                                                +  /   delta   <g|T |a'><g|T |b'> J  (w   ) |
//	     	                                     ---      a,b    1        2      sl  b'g  |
//                                                    g                                       ]

    {

    std::string slab[3];					// Array of labels for the signs
    slab[0] = std::string(" - ");				// These will be sign labels
    slab[1] = std::string(" + ");
    slab[2] = std::string("   ");
    std::string sgnspace(3, ' ');				// Spaces the length of printed sign

    std::string elab[3];					// Array of labels for the matrix
    if(add)						// If terms just added to previously
      {							// printed ones, no label, just spaces
      elab[0] = std::string(abs(add)+2, ' ');		// Note that three extra spaces are added
      elab[1] = std::string(abs(add)+2, ' ');		// due to the idea that a " =" was used
      elab[2] = std::string(abs(add)+2, ' ');		// in an earlier printed lines
      }
    else
      {							// Label is to be added in output
      elab[0] = Elabel[0] + std::string(2, ' ');		// Tack on 2 spaces here because of " ="
      elab[1] = Elabel[1] + std::string(" = ");		// Add an " = " to end of input label
      elab[2] = Elabel[2] + std::string(2, ' ');		// Tack on 2 spaces here because of " ="
      }
    std::string espace(elab[1].length(), ' ');		// Space the length of the "label"

    if(!add && !nterms) 				// If the label is to be printed
      {							// yet there are no terms, then just state
      ostr << "\n" << elab[0];			 	// that it is indeed zero
      ostr << "\n" << elab[1] << std::string("0");
      ostr << "\n" << elab[2];
      }

    int term = 0;
//    int term, ntermi = 0;
    int i, j;
    for(i=0; i<nterms; i+=ncols)
      {
      ostr << "\n";					// Skip down to Line 1
      for(j=0,term=i; j<ncols && term<nterms; j++)	// Print out Line 1, column j 
        { 						// 	If 1st column and 1st term  
        if(j==0 && i==0 && !add) ostr << elab[0]; 	//	and label added, print the label
	  else if(j==0)  ostr << espace;		// 	line 1, else spaces if 1st column
        ostr << sgnspace << line1[term];		//	Now output spaces for sign and term
        term++;						//	Increment term count
        }
      ostr << "\n";					// Skip down to Line 2
      for(j=0,term=i; j<ncols && term<nterms; j++)	// Print out Line 2, column j
        { 						// 	If 1st column and 1st term  
        if(j==0 && i==0 && !add) ostr << elab[1];	// 	and not added, print the label
        else if(j==0) ostr << espace; 			// 	line 2, else spaces if 1st column
        if(signs[term])					// Now put in a + or - sign unless it
          { 						// is a positive sign following a label
          if(j==0 && i==0 && !add) ostr << slab[2];	//	Put out no sign
          else ostr << slab[1];				//	Put out a + sign
          }
        else ostr << slab[0];				// 	Put out a - sign
        ostr << line2[term];				//	Now output the term
        term++;						//	Increment term count
        }
      ostr << "\n";					// Skip down to Line 3
      for(j=0,term=i; j<ncols && term<nterms; j++)	// Print out Line 3, column j
        { 						// 	If 1st column and 1st term  
        if(j==0 && i==0 && !add) ostr << elab[2]; 	//	and label added, print the label
	  else if(j==0) ostr << espace;			// 	line 3, else spaces if 1st column
        ostr << sgnspace << line3[term];		//	Now output spaces for sign and term
        term++;						//	Increment term count
        }
      }
    return;
    }

// ______________________________________________________________________
// ************ MISCELLANEOUS USEFUL FUNCTIONS FOR RELAXATION ***********
// ______________________________________________________________________


// KY - Commenting out this function - not consistent with that in relaxBWR.*
//   void ask_relax(int argc, char* argv[], int& argn,
//                  super_op& R, const sys_dynamic& sys, gen_op& H, int pflag)
// 
//         // Input        argc    : Number of command line arguments
//         //              argv    : Command line arguments
//         //              argn    : Initial command line argument
//         //                        for relaxation parameters
//         //              R       : Relaxation superoperator
//         //              sys     : Dynamic spin system
//         //              H       : Isotropic Hamiltonian
// 	//		pflag   : Flag for printing computaions
//         // Output       none    : Function is void.  The relaxation
//         //                        matrix R has different effects added
//         //                        in depending upon user requests
// 
//   //    {
//   // int type=0;						// Set for cross relaxation
// // int level=4;					// Set for non-secular  
// //  std::string dip = "n", dipdfs = "n",			// Set for no dipolar
//     //         csa = "n", csadfs = "n",			// Set for no shift anisotropy
//   //     dipcsa = "n",				// Set for no dipole-shift anisotropy
// 	 //       quad = "n", quaddfs = "n",			// Set for no quadrupolar			
// 		  //   dipquad = "n", csaquad = "n",		// Set for no dipole-quad, csa-quad
// 		       // rdm = "n";					// Set for no random field				
// 
// //			   See About Dipolar Relaxation
// 
//     query_parameter(argc, argv, argn,                   // Dipolar relaxation?
//      "\n\tInclude Dipolar Relaxation (y/n)? ", dip);
//     argn++;
//     if(dip=="y")
//       {
//       query_parameter(argc, argv, argn,			// Dipolar DFS relaxation effects
//          "\n\tInclude Dipolar Dynamic Frequency Shifts (y/n)? ", dipdfs);
//       argn++;
//       }
// 
// //			     See About CSA Relaxation
// 
//     int SA = 0;						// Assume no CSA tensors
//     int i;
//     for(i=0; i<sys.spins() && !SA; i++)			// See if any CSA tensors
//       if((sys.TC(i)).exists())				// are present in sys
//         SA = 1;
//     if(SA)						// If there are CSA tensors
//       {							// see if CSA relaxation is
//       query_parameter(argc, argv, argn,			// desired
//             "\n\tInclude CSA Relaxation (y/n)? ", csa);
//       argn++;
//       if(csa=="y")
//         {
//         query_parameter(argc, argv, argn,		// CSA DFS relaxation effects
//           "\n\tInclude CSA Dynamic Frequency Shifts (y/n)? ", csadfs);
//         argn++;
//         }
//       if((dip=="y") && (csa=="y"))
//         {
//         query_parameter(argc, argv, argn,                   // Dipole-CSA relaxation
//              "\n\tInclude Dipole-CSA Relaxation (y/n)? ", dipcsa);
//         argn++;
//         }
//       }
// 
// //		     See About Quadrupolar Relaxation
// 
//     int Q = 0;						// Assume no Quadrupolar tensors
//     for(i=0; i<sys.spins() && !Q; i++)		// See if any Quadrupolar tensors
//       if((sys.TQ(i)).exists())				// are present in sys
//         Q = 1;
//     if(Q)						// If there are Quadrupolar tensors
//       {							// see if Quadrupolar relaxation is
//       query_parameter(argc, argv, argn,			// desired
//             "\n\tInclude Quadrupolar Relaxation (y/n)? ", quad);
//       argn++;
//       if(quad=="y")
//         {
//         query_parameter(argc, argv, argn,		// Quadrupolar DFS relaxation effects
//           "\n\tInclude Quadrupolar Dynamic Frequency Shifts (y/n)? ", quaddfs);
//         argn++;
//         }
//       if((dip=="y") && (quad=="y"))
//         {
//         query_parameter(argc, argv, argn,                   // Dipole-Quadrupolar relaxation
//              "\n\tInclude Dipolar-Quadrupolar Relaxation (y/n)? ", dipquad);
//         argn++;
//         }
//       if((csa=="y") && (quad=="y"))
//         {
//         query_parameter(argc, argv, argn,                   // CSA-Quadrupolar relaxation
//              "\n\tInclude CSA-Quadrupolar Relaxation (y/n)? ", csaquad);
//         argn++;
//         }
//       }
// 
// //		     See About Random Field Relaxation
// 
//     query_parameter(argc, argv, argn,                   // Random Field relaxation
//      "\n\tInclude Random Field Relaxation (y/n)? ", rdm);
//     argn++;
// 
// //		     Now Compute All the Relaxation Effects
// 
//     if(dip == "y")
//       {
//       if(pflag) std::cout << "\n\tComputing The Dipolar Relaxation Matrix";
//       R += RDD(sys, H, type, level);
//       }
//     complex icmplx(0,1);                                // z = 0 + 1i
//     if(dipdfs == "y")
//       {
//       if(pflag) std::cout << "\n\tComputing The Dipolar DFS Relaxation Matrix";
//       R += icmplx*RDDds(sys, H, type, level);
//       }
//     if(csa == "y")
//       {
//       if(pflag) std::cout << "\n\tComputing The CSA Relaxation Matrix";
//       R += RCC(sys, H, type, level);
//       }
//     if(csadfs == "y")
//       {
//       if(pflag) std::cout << "\n\tComputing The CSA DFS Relaxation Matrix";
//       R += icmplx*RCCds(sys, H, type, level);
//       }
//     if(dipcsa == "y")
//       {
//       if(pflag) std::cout << "\n\tComputing The Dipole-CSA Relaxation Matrix";
//       R += RDCX(sys, H, level);
//       }
//     if(quad == "y")
//       {
//       if(pflag) std::cout << "\n\tComputing The Quadrupolar Relaxation Matrix";
//       R += RQQ(sys, H, type, level);
//       }
//     if(quaddfs == "y")
//       {
//       if(pflag) std::cout << "\n\tComputing The Quadrupolar DFS Relaxation Matrix";
//       R += icmplx*RQQds(sys, H, type, level);
//       }
// // sosi - need to build these functions
// //    if(dipquad == "y")
// //      {
// //      if(pflag) std::cout << "\n\tComputing The Dipole-Quadrupolar Relaxation Matrix";
// //      R += RDQX(sys, H, level);
// //      }
// //    if(csaquad == "y")
// //      {
// //      if(pflag) std::cout << "\n\tComputing The CSA-Quadrupolar Relaxation Matrix";
// //      R += RCQX(sys, H, level);
// //      }
//     if(rdm == "y")
//       {
//       if(pflag) std::cout << "\n\tComputing The Random Field Relaxation Matrix";
//       R += complex(5*PI/2)*RRR(sys, H, type, level);
//       }
//     return;
//     }


  void sort(int* indx, matrix& mx, int k, int type, int colf)

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

  {
  int i=0,j=0;
  double *vals;			// Array of values to be sorted
  int nvals = mx.rows();	// Number of values to be sorted
  if(colf) nvals = mx.cols();	// Sort the columns based on values
  vals = new double[nvals];
  if(type > 0) type = 1;
  else if(type) type = 2;
  for(i=0; i<nvals; i++)	// Fill vals with values to be sorted
    {				// and intialize the indices
    indx[i] = i;
    switch(type)
      {
      case 0:			// Sort real values
      default:
        if(colf)
          vals[i] = Re(mx.get(k,i));
        else
          vals[i] = Re(mx.get(i,k));
        break;
      case 1:			// Sort norms
        if(colf)
          vals[i] = norm(mx.get(k,i));
        else
          vals[i] = norm(mx.get(i,k));
        break;
      case 2:			// Sort imaginary values
        if(colf)
          vals[i] = Im(mx.get(k,i));
        else
          vals[i] = Im(mx.get(i,k));
        break;
      }
    }
  double maxval;
  int maxind;
  for(i=0; i<nvals; i++)	// Now perform the sort
    {
    maxval = vals[i];		// Assume this is maximum
    maxind = i;			// This is it's index
    for(j=i+1; j<nvals; j++)	// Now compar with all others
      {
      if(vals[j] > maxval)
        {
        maxval = vals[j];	//	This is greater
        maxind = j;		//	This is its index
        }
      }
    vals[maxind] = vals[i];	// Swap the values
    j = indx[maxind];		// Swap the absolute indices
    indx[maxind] = indx[i];	// Now fix the stored indces
    indx[i] = j;
    }
  delete [] vals;		// Delete the values
  return;
  }


#endif /* __RELAX_ANALYZE_CC__ */ 
 
