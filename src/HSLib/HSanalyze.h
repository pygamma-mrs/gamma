/* HSanalyze.h *************************************************-*-c++-*-*
**                                                                      **
**	                              G A M M A				  **
**                                                                      **
**	MR Hilbert Space Analyze                     Interface		  **
**                                                                      **
**	 Copyright (c) 1993                                              **
**	 Scott Smith					 		         **
**      University of Utah                      	                       **
**      Department of Chemistry                         	         **
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
** of Spin Hilbert Spaces Associated with Composite Spin Systems.       **
**                                                                      **
*************************************************************************/

#ifndef   HSanalyze_h_				// Is file already included?
#  define HSanalyze_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// This is the interface
#  endif

#include <GamGen.h>				// Know MSVCDLL (__declspec)
#include <HSLib/GenOp.h>                         // Include operators
#include <string>                                // Know about strings
const std::string alphabeta[7]                   // Spin state labels
                = {"a","b","g","d","e","w","x"};	

// ____________________________________________________________________________
// A                             WAVEFUNCTIONS
// ____________________________________________________________________________


/* The wavefunctions will span the Hilbert space of the operator Op. Each
   wavefunction will be some linear combination of the default basis functions
   GAMMA uses that are inherent to spin system sys.  The linear combination
   will depend upon the current Op basis.  If Op is in its eigenbasis, wave-
   functions are eigenfunctions.  Each wavefunction corresponds to a state
   of the system.                                                            */

 
MSVCDLL std::string qStatel(const spin_sys& sys, int bf);
 
        // Input                sys        : A basic spin system
        //                      bf         : A state (basis function) index
        // Output               bflabel    : A label for the input state in terms
        //                                   of individual spin states.  For example
        //                                   ab for up-down nn {I1=I2=1/2}, or gb
        //                                   for down-middle in {I1=1,I2=1/2}
        // Note                            : The string alphabeta is set to
        //                                   handle spins with maximum I=7/2
        // Note                            : Individual spin labels are as follows:
        //                                      I=1/2: 1/2 - a; -1/2 - b;
        //                                      I=1:    1  - a;   0  - b;  -1  - g;
        //                                      I=3/2: 3/2 - a;  1/2 - b; -1/2 - g; -3/2 - d;
 
 
MSVCDLL void wf_labels(std::string* wflabels, const spin_sys& sys, gen_op &Op,
                                    double cutoff=1.e-4, int pbf=1, int pfz=1);
 
        // Input        wflabels: An array of strings
        //              sys     : A basic spin system
        //              Op      : General operator (associated to sys)
        //              cutoff  : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT
        // Return       void    : The string array wflabels is filled with system
        //                        wavefunctions in the current working basis
        //                        of operator Op in the format specified
        // Note                 : Handles complex basis function coefficients!
        // Note                 : Wavefunction ordering depends on Op input basis!


MSVCDLL void wf_labels(std::string* wflabels, int* filter, const spin_sys& sys,
                     const gen_op &Op, double cut=1.e-4, int pbf=1, int pfz=1);

        // Input        wflabels: An array of strings
        //              filter  : An array of integers, output filter
        //              sys     : A basic spin system
        //              Op      : General operator (associated to sys)
        //              cut     : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT
        // Return       void    : The string array wflabels is filled with system
        //                        wavefunctions in the current working basis
        //                        of operator Op in the format specified
        // Note                 : Handles complex basis function coefficients!
        // Note                 : Wavefunction ordering depends on Op input basis!


MSVCDLL void wf_labels(std::string* wflabels, const spin_sys& sys, const matrix &B,
                                    double cutoff=1.e-4, int pbf=1, int pfz=1);
 
        // Input        wflabels: An array of strings
        //              sys     : A basic spin system
        //              B       : A basis (associated to sys)
        //              cutoff  : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT
        // Return       void    : The string array wflabels is filled with system
        //                        wavefunctions in the current working basis
        //                        of operator Op in the format specified
        // Note                 : Handles complex basis function coefficients!
        // Note                 : Wavefunction ordering depends on Op input basis!

 
MSVCDLL void wf_labels(std::string* wflabels, int* filter, const spin_sys& sys,
                      const matrix &B, double cut=1.e-4, int pbf=1, int pfz=1);
 
        // Input        wflabels: An array of strings
        //              filter  : An array of integers, output filter
        //              sys     : A basic spin system
        //              B       : A basis (associated to sys)
        //              cut     : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT
        // Return       void    : The string array wflabels is filled with system
        //                        wavefunctions in the current working basis
        //                        of operator Op in the format specified
        // Note                 : Handles complex basis function coefficients!
        // Note                 : Wavefunction ordering depends on Op input basis!


// ____________________________________________________________________________ 
// B                                EIGENVALUES 
// ____________________________________________________________________________ 
 
/* The eigenvalues will span the Liouville space of the superoperator LOp.
   Each eigenvalue will correspond to an element of LOp after LOp has been
   placed into it representation in an eigenbasis.  Each eigenvalue will then
   correspond to a particular eigenfunction of the superoperator.  If the 
   superoperator is a commutation superoperator of the Hamiltonian, then the
   eigenvalues will be transition energies of the system.                    */


MSVCDLL void ev_labels(std::string* evlabels, gen_op& Op, double cutoff=1.e-6);

        // Input        evlabels: An array of strings
        //              Op      : General operator
        //              cutoff  : Eigenvalues of magnitude below this are
        //                        assumed to be zero
        // Return       evlabels: An array of strings containg all the
        //                        operator eigenvalues (in the eigenbasis)
        //                        of operator Op.
        // Note                 : Handles complex eigenvalues!
        // Note                 : Eigenvalue ordering depends on Op eigenbasis.
        //                        Either order the eigenbasis prior to entering
        //                        function or sort the values after output
 

MSVCDLL void ev_labels(std::string* evlabels, int* filter, gen_op& Op, double cutoff=1.e-6);

        // Input        evlabels: An array of strings
        //              filter  : Array of integers to filter values
        //              Op      : General operator
        //              cutoff  : Eigenvalues of magnitude below this are
        //                        assumed to be zero
        // Return       evlabels: An array of strings containg all the
        //                        operator eigenvalues (in the eigenbasis)
        //                        of operator Op.
        // Note                 : Handles complex eigenvalues!
        // Note                 : Eigenvalue ordering depends on Op eigenbasis.
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


MSVCDLL void tref_labels(std::string* trlabels, const spin_sys& sys, gen_op &Op,
                           int type=0, double cut=1.e-4, int pbf=1, int pfz=0);

        // Input        trlabels: An array of strings
        //              sys     : A basic spin system
        //              Op      : General operator (associated to sys)
        //              type    : Flag for type of transition labels
        //                        0 = label states by number (DEFAULT)
        //                        1 = label states by eigenfunctions 
        //              cut     : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //                              INACTIVE IF TYPE = 0
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //                              INACTIVE IF TYPE = 0
        //              pfz     : Flag for inclusion of total Fz 
        //                              INACTIVE IF TYPE = 0
        // Return       void    : The string array trlabels is filled with system
        //                        transitions ordered from the eigenbasis
        //                        of operator Op in the format specified
        // Note                 : Transition ordering depends on Op input basis!


MSVCDLL void tref_labels(std::string* trlabels, const spin_sys& sys, const matrix &B,
                        int type=0, double cutoff=1.e-4, int pbf=1, int pfz=0);

        // Input        trlabels: An array of strings
        //              sys     : A basic spin system
        //              B       : Eigenbasis matrix (associated to sys)
        //              type    : Flag for type of transition labels
        //                        0 = label states by number (DEFAULT)
        //                        1 = label states by eigenfunctions 
        //              cutoff  : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //                              INACTIVE IF TYPE = 0
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //                              INACTIVE IF TYPE = 0
        //              pfz     : Flag for inclusion of total Fz
        //                              INACTIVE IF TYPE = 0
        // Return       void    : The string array trlabels is filled with system
        //                        transitions ordered from the eigenbasis
        //                        of operator Op in the format specified
        // Note                 : Transition ordering depends on Op input basis!
        // Note                 : There is no check, but B should be an eigenbasis
        //                        of some operator associated with sys


// ------------------------ Hilbert Transition Types --------------------------


MSVCDLL void tran_types(std::string* trtypes, const spin_sys& sys, gen_op &Op,
                                                 int type=0, double cut=1.e-4);

	// Input 	trtypes : An array of strings
	//       	sys	: A basic spin system
	// 		Op	: General operator (associated to sys)
	//		type    : Flag for type of transition labels
	//			  0 = label states by number (DEFAULT)
	//			  1 = label states by eigenfunctions 
	//		cut     : Basis functions with coefficients below
	//			  this value are not taken as part of a wf
	//			        INACTIVE IF TYPE = 0
	//		pbf	: Flag for default basis function format
	//			   <=0 - by number, [0, hs)
	//			   >=0 - by spin states, e.g. abaa (DEFAULT)
	//			        INACTIVE IF TYPE = 0
        //		pfz	: Flag for inclusion of total Fz 
	//			        INACTIVE IF TYPE = 0
	// Return	void    : The string array trtypes is filled with system
	//			  transitions ordered from the eigenbasis
	//			  of operator Op in the format specified
	// Note			: Transition ordering depends on Op input basis!


MSVCDLL void tran_types(std::string* trtypes, const spin_sys& sys, const matrix &B,
                                              int type=0, double cutoff=1.e-4);

	// Input 	trtypes: An array of strings
	//       	sys	: A basic spin system
	// 		B 	: Eigenbasis matrix (associated to sys)
	//		type    : Flag for type of transition labels
	//			  0 = label states by number (DEFAULT)
	//			  1 = label states by eigenfunctions 
	// Return	void    : The string array trtypes is filled with system
	//			  transitions ordered from the eigenbasis
	//			  of operator Op in the format specified
	// Note			: Transition ordering depends on Op input basis!
	// Note			: There is no check, but B should be an eigenbasis
	//			  of some operator associated with sys


// ----------------------- Hilbert Transition Energies ------------------------

MSVCDLL void trev_labels(std::string* trlabels, gen_op& Op, double cutoff=1.e-6);

        // Input        trlabels: An array of strings
        //              Op      : General operator
        //              cutoff  : Transitions of magnitude below this are
        //                        assumed to be zero
        // Return       trlabels: An array of strings containg all the
        //                        operator transition values (in the eigenbasis)
        //                        of operator Op.
        // Note                 : Handles complex transitions!
        // Note                 : Transition ordering depends on Op eigenbasis.
        //                        Either order the Op eigenbasis prior to entering
        //                        function or sort the values after output


// ____________________________________________________________________________
// D                       EIGENSYSTEM OUTPUT FUNCTIONS
// ____________________________________________________________________________


// ------------------ System Wavefunctions in Hilbert Space -------------------

MSVCDLL void wavefunctions(std::ostream& ostr, const spin_sys& sys, gen_op &Op,
                       double cutoff=1.e-4, int pbf=1, int pfz=1, int title=1);
 
        // Input        ostr    : An output stream
        //              sys     : A spin system  
        //              Op      : General operator (associated to sys)
        //              cutoff  : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              pbf     : Flag for default basis function format 
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT)
        //              title   : Flag for inclusion of a title
        // Return       void    : The system wavefunctions are sent into the
        //                      : supplied output stream ostr.


MSVCDLL void wavefunctions(std::ostream& ostr, int* filt, const spin_sys& sys, gen_op &Op,
                       double cutoff=1.e-4, int pbf=1, int pfz=1, int title=1);
 
        // Input        ostr    : An output stream
        //              filt    : An array of integers, output filter
        //              sys     : A spin system
        //              Op      : General operator (associated to sys)
        //              cutoff  : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT)
        //              title   : Flag for inclusion of a title
        // Return       void    : The system wavefunctions are sent into the
        //                      : supplied output stream ostr.
 

// -----------------------Eigensystem in Hilbert Space ------------------------

MSVCDLL void eigensystem(std::ostream& ostr, const spin_sys& sys, gen_op& Op,
      double cute=1.e-6, double cutc=1.e-4, int pbf=1, int pfz=1, int title=1);

        // Input        ostr    : An output stream
        //              sys     : A spin system
        //              Op      : A general operator
        //              cute    : Eigenvalues below this are assumed 0
        //              cutc    : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT)
        //              title   : Flag for inclusion of a title
        // Return       none    : All eigenvalues and eigenfunctions
        //                        of operator Op are sent into ostream ostr


MSVCDLL void eigensystem(std::ostream& ostr, int* filt, const spin_sys& sys, gen_op& Op,
      double cute=1.e-6, double cutc=1.e-4, int pbf=1, int pfz=1, int title=1);
 
        // Input        ostr    : An output stream
        //              filt    : An array of integers, output filter
        //              sys     : A spin system
        //              Op      : A general operator
        //              cute    : Eigenvalues below this are assumed 0
        //              cutc    : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT)
        //              title   : Flag for inclusion of a title
        // Return       none    : All eigenvalues and eigenfunctions
        //                        of operator Op are sent into ostream ostr


// --------------------- System Transitions in Hilbert Space ------------------

MSVCDLL void transitions(std::ostream& ostr, const spin_sys& sys, gen_op& Op, int type=0,
      double cutt=1.e-6, double cutc=1.e-4, int pbf=1, int pfz=1, int title=1);

        // Input        ostr    : An output stream
        //              sys     : A spin system
        //              Op      : A general operator
        //              type    : Flag for type of transition labels
        //                        0 = label states by number (DEFAULT)
        //                        1 = label states by eigenfunctions 
        //              cutt    : Transitions below this are assumed 0
        //              cutc    : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //                              INACTIVE IF TYPE = 0
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //                              INACTIVE IF TYPE = 0
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT)
        //                              INACTIVE IF TYPE = 0
        //              title   : Flag for inclusion of a title
        // Return       none    : All transitions and transition values
        //                        of operator Op are sent into ostream ostr


MSVCDLL void transitions(std::ostream& ostr, int* filt, const spin_sys& sys, gen_op& Op, int type=0,
         double cutt=1.e-6, double cutc=1.e-4, int pbf=1, int pfz=1, int title=1);

        // Input        ostr    : An output stream
        //              filt    : An array of integers, output filter
        //              sys     : A spin system
        //              Op      : A general operator
        //              type    : Flag for type of transition labels
        //                        0 = label states by number (DEFAULT)
        //                        1 = label states by eigenfunctions 
        //              cutt    : Transitions below this are assumed 0
        //              cutc    : Basis functions with coefficients below
        //                        this value are not taken as part of a wf
        //                              INACTIVE IF TYPE = 0
        //              pbf     : Flag for default basis function format
        //                         <=0 - by number, [0, hs)
        //                         >=0 - by spin states, e.g. abaa (DEFAULT)
        //                              INACTIVE IF TYPE = 0
        //              pfz     : Flag for inclusion of total Fz   (DEFAULT)
        //                              INACTIVE IF TYPE = 0
        //              title   : Flag for inclusion of a title
        // Return       none    : All transitions and transition values
        //                        of operator Op are sent into ostream ostr


// ____________________________________________________________________________
// E                             SELECTION RULES
// ____________________________________________________________________________
 
// -------------------- Selection Rules in Hilbert Space ----------------------
 
MSVCDLL void ev_select(int* select, gen_op &Op, double val1, int type=0,
                              double val2=0, double cutoff=1.e-4 , int reim=1);
 
        // Input        select  : An array of integers
        //              Op      : General operator (associated to sys)
        //              val1    : Value used in determining the selection
        //              type    : Type of selection applied
        //                         0 - Only levels having ev == val1 (DEFAULT)
        //                         1 - All levels having ev != val1
        //                         2 - All levels having ev >= val1
        //                         3 - All levels having ev <= val1
        //                         4 - All levels having val2 <= ev <= val1
        //              val2    : Value used in determining the selection,
        //                        unused unless type=4, (0 DEFAULT) 
        //              cutoff  : If fabs(ev-val1) < cutoff, assume equal
        //              reim    : Type of ev used
        //                         0 - Use the norm of the complex ev
        //                        >0 - Use the real component of the ev (DEFAULT)
        //                        <0 - Use the imaginary component of ev
        // Return       void    : The integer array is filled with 0's and 1's
        //                        depending on whether the system sys wavefunctions
        //                        satisfy the input selection rules.
        // Note                 : Selection ordering depends on Op input basis!

 
MSVCDLL void tr_select(int* select, gen_op &Op, double val1, int type=0,
                              double val2=0, double cutoff=1.e-4 , int reim=1);

        // Input        select  : An array of integers
        //              Op      : General operator (associated to sys)
        //              val1    : Value used in determining the selection
        //              type    : Type of selection applied
        //                         0 - Only transitions having dev == val1 (DEFAULT)
        //                         1 - All transitions having dev != val1
        //                         2 - All transitions having dev >= val1
        //                         3 - All transitions having dev <= val1
        //                         4 - All transitions having val2 <= dev <= val1
        //              val2    : Value used in determining the selection,
        //                        unused unless type=4, (0 DEFAULT) 
        //              cutoff  : If fabs(dev-val1) < cutoff, assume equal
        //              reim    : Type of dev used
        //                         0 - Use the norm of the complex dev
        //                        >0 - Use the real component of the dev (DEFAULT)
        //                        <0 - Use the imaginary component of dev
        // Return       void    : The integer array is filled with 0's and 1's
        //                        depending on whether the system sys wavefunctions
        //                        satisfy the input selection rules.
        // Note                 : Selection ordering depends on Op input basis!


#endif							                                  // HSanalyze.h

