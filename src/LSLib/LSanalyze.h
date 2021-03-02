/* LSanalyze.h *************************************************-*-c++-*-*
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      MR Liouville Space Analyze                   Interface 		**
**                                                                      **
**      Copyright (c) 1993                                              **
**      Scott Smith                                                     **
**      University of Utah                                              **
**      Department of Chemistry                                         **
**      Salt Lake City, UT, 84112, USA                                  **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The GAMMA Platform Provides Functions for Simulation of Magnetic     **
** Resonance Experiments and Other Associated Mathematical              **
** Capabilities.  The Set of Functions Herein Allows for the Analysis   **
** of Spin Liouville Spaces Associated with Composite Spin Systems.     **
** This module of the GAMMA MR Library provides functions designed to   **
** facilitate the analysis of the more common aspects of MR with        **
** respect to a quantum mechanical viewpoint.  As GAMMA provides full   **
** access to all its defined quantities, these routines serve just to   **
** sort information and present it in an orderly fashion.               **
**                                                                      **
*************************************************************************/

#ifndef   LSanalyze_h_			// Is file already included?
#  define LSanalyze_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// This is the interface
#  endif

#include <GamGen.h>			// Include system specifics
#include <HSLib/GenOp.h>		// Include Hilbert operators
#include <LSLib/SuperOp.h>		// Include Liouville operators
#include <Basics/Gconstants.h>		// Include constants
#include <string>			// Include libstdc++ strings
 
// ____________________________________________________________________________
// A                     WAVEFUNCTIONS / BASIS FUNCTIONS
// ____________________________________________________________________________

/* The wavefunctions will span the Liouville space of the superoperator LOp.
   Each wavefunction will be some linear combination of the default super-
   operator basis functions.  The default superoperator basis functions are
   inherent to each superoperator.  They are determined directly from the
   differences of Hilbert space wavefunctions, and these Hilbert space wave-
   functions are contained in the Hilbert space basis internally present in
   LOp.  Each of these Hilbert space wavefunctions are related to the GAMMA
   default basis functions (inherent to spin system sys).  If LOp is in its
   eigenbasis, its wavefunctions are eigenfunctions.  If the LOp Hilbert
   space basis an eigenbasis, then the LOp wavefunctions are linear com-
   binations of Hilbert space eigenfunction differences.  Each wavefunction
   corresponds to a coherence between two state of the system.               */

 
MSVCDLL void wf_labels(std::string* wflabels, const spin_sys& sys, super_op &LOp,
                       double cut=1.e-4, int type=-1, int pbf=-1, int pfz=0);
 
        // Input        wflabels: An array of strings
        //              sys     : A basic spin system
        //              LOp     : A superoperator (associated to sys)
        //              cut     : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              type    : Flag for type of transition labels
        //                         <0 - by Liouville number, [0, ls) (DEFAULT)
        //                          0 - Hilbert transitions, {1->2, 2->3,..}
        //                         >0 - Hilbert transitions, {a|aa>->b|ab>, ...}
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //                              INACTIVE IF TYPE <= 0
        //              pfz     : Flag for inclusion of total Fz in Hilbert basis
        //                              INACTIVE IF TYPE <= 0
        // Return       wflabels: An array of strings containg all system
        //                        wavefunctions in the current working basis
        //                        of operator LOp in the format specified
        // Note                 : Handles complex basis function coefficients!
        // Note                 : Wavefunction ordering depends on LOp input basis!


MSVCDLL void wf_labels(std::string* wflabels, int* filter, const spin_sys& sys, super_op &LOp,
                         double cut=1.e-4, int type=-1, int pbf=-1, int pfz=0);

        // Input        wflabels: An array of strings
        //              filter  : An array of integers, output filter
        //              sys     : A basic spin system
        //              LOp     : A superoperator (associated to sys)
        //              cut     : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              type    : Flag for type of transition labels
        //                         <0 - by Liouville number, [0, ls) (DEFAULT)
        //                          0 - Hilbert transitions, {1->2, 2->3,..}
        //                         >0 - Hilbert transitions, {a|aa>->b|ab>, ...}
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //                              INACTIVE IF TYPE <= 0
        //              pfz     : Flag for inclusion of total Fz in Hilbert basis
        //                              INACTIVE IF TYPE <= 0
        // Return       wflabels: An array of strings containg all system
        //                        wavefunctions in the current working basis
        //                        of operator LOp in the format specified
        // Note                 : Handles complex basis function coefficients!
        // Note                 : Wavefunction ordering depends on LOp input basis!


MSVCDLL void wf_labels(std::string* wflabels, const spin_sys& sys, const matrix &B,const matrix& HB, double cut=1.e-4, int type=-1, int pbf=-1, int pfz=0);

        // Input        wflabels: An array of strings
        //              sys     : A basic spin system
        //              LB      : A superoperator basis
        //              HB      : A Hilbert basis (associated to sys)
        //              cut     : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              type    : Flag for type of transition labels
        //                         <0 - by Liouville number, [0, ls) (DEFAULT)
        //                          0 - Hilbert transitions, {1->2, 2->3,..}
        //                         >0 - Hilbert transitions, {a|aa>->b|ab>, ...}
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //                              INACTIVE IF TYPE <= 0
        //              pfz     : Flag for inclusion of total Fz in Hilbert basis
        //                              INACTIVE IF TYPE <= 0
        // Return       wflabels: An array of strings containg all system
        //                        wavefunctions in the current working basis
        //                        of operator LOp in the format specified
        // Note                 : Handles complex basis function coefficients!
        // Note                 : Wavefunction ordering depends on LOp input basis!

 
MSVCDLL void wf_labels(std::string* wflabels, int* index, const spin_sys& sys, const matrix &B,const matrix& HB, double cut=1.e-4, int type=-1, int pbf=-1, int pfz=0);
 
        // Input        wflabels: An array of strings
        //              filter  : An array of integers, output filter
        //              sys     : A basic spin system
        //              LB      : A superoperator basis
        //              HB      : A Hilbert basis (associated to sys)
        //              cut     : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              type    : Flag for type of transition labels
        //                         <0 - by Liouville number, [0, ls) (DEFAULT)
        //                          0 - Hilbert transitions, {1->2, 2->3,..}
        //                         >0 - Hilbert transitions, {a|aa>->b|ab>, ...}
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //                              INACTIVE IF TYPE <= 0
        //              pfz     : Flag for inclusion of total Fz in Hilbert basis
        //                              INACTIVE IF TYPE <= 0
        // Return       wflabels: An array of strings containg all system
        //                        wavefunctions in the current working basis
        //                        of operator LOp in the format specified
        // Note                 : Handles complex basis function coefficients!
        // Note                 : Wavefunction ordering depends on LOp input basis!
 

// ____________________________________________________________________________
// B                                EIGENVALUES
// ____________________________________________________________________________

/* The eigenvalues will span the Liouville space of the superoperator LOp.
   Each eigenvalue will correspond to an element of LOp after LOp has been
   placed into it representation in an eigenbasis.  Each eigenvalue will then
   correspond to a particular eigenfunction of the superoperator.  If the
   superoperator is a commutation superoperator of the Hamiltonian, then the
   eigenvalues will be transition energies of the system.                    */  

 
MSVCDLL void ev_labels(std::string* evlabels, super_op& LOp, double cutoff=1.e-6);
 
        // Input        evlabels: An array of strings
        //              LOp     : A superoperator
        //              cutoff  : Eigenvalues of magnitude below this are
        //                        assumed to be zero
        // Return       evlabels: An array of strings containg all the
        //                        operator eigenvalues (in the eigenbasis)
        //                        of operator LOp.
        // Note                 : Handles complex eigenvalues!
        // Note                 : Eigenvalue ordering depends on LOp eigenbasis.
        //                        Either order the eigenbasis prior to entering
        //                        function or sort the values after output


MSVCDLL void ev_labels(std::string* evlabels, int* filter, super_op& LOp, double cutoff=1.e-6);

        // Input        evlabels: An array of strings
        //              filter  : Array of integers to filter values
        //              LOp     : A superoperator
        //              cutoff  : Eigenvalues of magnitude below this are
        //                        assumed to be zero
        // Return       evlabels: An array of strings containg all the
        //                        operator eigenvalues (in the eigenbasis)
        //                        of operator LOp.
        // Note                 : Handles complex eigenvalues!
        // Note                 : Eigenvalue ordering depends on LOp eigenbasis.
        //                        Either order the eigenbasis prior to entering
        //                        function or sort the values after output

 
// ____________________________________________________________________________
// C                             TRANSITIONS
// ____________________________________________________________________________

/* Transitions between eigenstates are associated with two different viewpoints 
   1.) It is taken as a jump from one eigenfunction to another eigenfunction,
       so there are two eigenfunctions associated with each transition.

   2)  The transition has a value assocated with it, the value being the
       difference between two eigenvalues.  There is only one transition value,
       the net of two eigenvalues.  If operator Op is the system Hamiltonian
       then the transition value is the transition frequency.  If operator Op
       is Fz, then the value is the net Fz associated with the transtion.

   Note that transition values in Hilbert space (above) are equivalent to 
   single eigenvalues of superoperators in Liouville space.  Similarly,
   transitions between two eigenstates in Hilbert space are associated with
   a single wavefunction in Liouville space.                                 */

// ____________________________________________________________________________
// D                       EIGENSYSTEM OUTPUT FUNCTIONS
// ____________________________________________________________________________
 
// ------------------ System Wavefunctions in Liouville Space -----------------
 
MSVCDLL void wavefunctions(std::ostream& ostr, const spin_sys& sys, super_op &LOp,
         double cutoff=1.e-4, int type=-1, int pbf=-1, int pfz=0, int title=1);
 
        // Input        ostr    : An output stream
        //              sys     : A spin system
        //              LOp     : A superoperator (associated to sys)
        //              cutoff  : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              type    : Flag for type of transition labels
        //                         <0 - by Liouville number, [0, ls) (DEFAULT)
        //                          0 - Hilbert transitions, {1->2, 2->3,..}
        //                         >0 - Hilbert transitions, {a|aa>->b|ab>, ...}
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //                              INACTIVE IF TYPE <= 0
        //              pfz     : Flag for inclusion of total Fz in Hilbert basis
        //                              INACTIVE IF TYPE <= 0
        //              title   : Flag for inclusion of a title
        // Return       void    : The system wavefunctions are sent into the
        //                      : supplied output stream ostr.


MSVCDLL void wavefunctions(std::ostream& ostr, int* filter, const spin_sys& sys, super_op &LOp,
              double cutoff=1.e-4, int type=-1, int pbf=-1, int pfz=0, int title=1);
 
        // Input        ostr    : An output stream
        //              filter  : An array of integers, output filter
        //              sys     : A spin system
        //              LOp     : A superoperator (associated to sys)
        //              cutoff  : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              type    : Flag for type of transition labels
        //                         <0 - by Liouville number, [0, ls) (DEFAULT)
        //                          0 - Hilbert transitions, {1->2, 2->3,..}
        //                         >0 - Hilbert transitions, {a|aa>->b|ab>, ...}
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //                              INACTIVE IF TYPE <= 0
        //              pfz     : Flag for inclusion of total Fz in Hilbert basis
        //                              INACTIVE IF TYPE <= 0
        //              title   : Flag for inclusion of a title
        // Return       void    : The system wavefunctions are sent into the
        //                      : supplied output stream ostr.
      

// -------------------- Eigensystem in Liouville Space ------------------------
 
MSVCDLL void eigensystem(std::ostream& ostr, const spin_sys& sys, super_op& LOp,
                    double cute=1.e-6, double cutc=1.e-4,
                            int type=-1, int pbf=-1, int pfz=0, int title=1);
 
        // Input        ostr    : An output stream
        //              sys     : A spin system
        //              LOp     : A superoperator
        //              cute    : Eigenvalues below this are assumed 0
        //              cutc    : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              type    : Flag for type of transition labels
        //                         <0 - by Liouville number, [0, ls) (DEFAULT)
        //                          0 - Hilbert transitions, {1->2, 2->3,..}
        //                         >0 - Hilbert transitions, {a|aa>->b|ab>, ...}
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //                              INACTIVE IF TYPE <= 0
        //              pfz     : Flag for inclusion of total Fz in Hilbert basis
        //                              INACTIVE IF TYPE <= 0
        //              title   : Flag for inclusion of a title
        // Return       none    : All eigenvalues and eigenfunctions
        //                        of operator LOp are sent into ostream ostr
 
// ____________________________________________________________________________
// E                             SELECTION RULES
// ____________________________________________________________________________
 
#endif 				               		// LSanalyze.h

