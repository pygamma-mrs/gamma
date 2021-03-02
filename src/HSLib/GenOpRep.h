/* GenOpRep.h ***************************************************-*-c++-*-
**                                                                      **
**                               G A M M A                              **
**                                                                      **
**      General Operator Represenation              Interface           **
**                                                                      **
**      Copyright (c) 1990, 1991, 1992                                  **
**      Scott Smith                                                     **
**      Eidgenoessische Technische Hochschule                           **
**      Labor fuer physikalische Chemie                                 **
**      8092 Zuerich / Switzerland                                      **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
** The class genoprep defines a single operator representaion for a     **
** variable of type gen_op, i.e. a GAMMA general operator.  Each        **
** representation consists of a matrix, a basis, and a priority #.      **
**                                                                      **
*************************************************************************/

#ifndef _GenOpRep_h_ 			// Is file already included?
#define _GenOpRep_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA) 			// Using the GNU compiler?
#    pragma interface			// If yes, then we are the interface
#endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Matrix/matrix.h>		// Know about GAMMA matrices
#include <HSLib/Basis.h>		// Know about operator bases

class genoprep
  {
public:
 
         matrix RepMx;			// Op matrix representation
         basis  RepBs;			// Associated Op basis (array)
         int    RepPty;			// Representation priority
  static bool   BSPrnt;			// Flag to print basis array

// ____________________________________________________________________________
// i                CLASS OPERATOR REPRESENTATION ERROR HANDLING
// ____________________________________________________________________________

	// Input		OpRep	: Operator representation (this)
        //                      eidx    : Error index
        //                      noret   : Flag for linefeed (0=linefeed)
        // Output               void    : An error message is output   (error)
        //                                or program execution stopped (fatal)
 
         void OpReperror(int eidx, int noret=0) const;
volatile void OpRepfatal(int eidx)              const;
 
// ____________________________________________________________________________
// A            OPERATOR REPRESENTATION CONSTRUCTORS/DETRUCTOR
// ____________________________________________________________________________

MSVCDLC genoprep();
MSVCDLC genoprep(const genoprep& OpRep);
MSVCDLC genoprep(const matrix& mx, const basis& bs, int pty);
MSVCDLC ~genoprep();
genoprep& operator= (const genoprep& OpRep); 

// ____________________________________________________________________________
// J                    OPERATOR REPRESENTATION CHECKS
// ____________________________________________________________________________


	// Input		OpRep: A operator representation
	// 			warn : A warning level
	//			         0 = no warning
	//			         1 = non-fatal warning
	//			         2 = fatal warnings
	// Output		bool : True if OprRep matrix mx square.
        //			       & its dimension matches its basis

MSVCDLL bool OpRepCheck(int warn=2) const;

// ____________________________________________________________________________
// K                  CLASS OPERATOR REPRESENTATION I/O FUNCTIONS
// ____________________________________________________________________________
 
// ------------------------ ASCII Output Functions ----------------------------
 
/*              Input           OpRep: Operator representation (this)
                                ostr : Output ASCII file stream
                                full : Flag for amount of output
                Return          void : OpRep is sent to the output stream    */
 
MSVCDLL        std::ostream& print(std::ostream& ostr, int full=0) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr,const genoprep &OpRep);
 
// ------------------------ Binary Output Functions ---------------------------

/*              Input           OpRep: Operator representation (this)           
                                fn   : Input binary filename
                                fp   : File stream (pointing at OpRep spot)
                Return          void : OpRep is written to either the
                                       specified file or filestream.
                Note                 : Output format is partially set by
                                       class matrix (matrix typing)          */  

MSVCDLL void           write(const std::string& fn) const;
MSVCDLL std::ofstream& write(std::ofstream& fp)     const;

// ------------------------ Binary Input Functions ----------------------------

/*              Input           OpRep: Operator representation (this)
                                fn   : Input binary filename
                                fp   : File stream (pointing at OpRep)
                Return          void : OpRep is read in from either the
                                       specified file or filestream.         */  

MSVCDLL void           read(const std::string& fn);
MSVCDLL std::ifstream& read(std::ifstream& fp);

// ____________________________________________________________________________
// L        OPERATOR REPRESENTATION LIST/VECTOR SUPPORT FUNCTIONS
// ____________________________________________________________________________

MSVCDLL bool operator==(const genoprep& OpRep) const;
MSVCDLL bool operator!=(const genoprep& OpRep) const;
MSVCDLL bool operator<(const  genoprep& OpRep) const;
MSVCDLL bool operator>(const  genoprep& OpRep) const;
 
};

#endif							// GenOpRep.h
