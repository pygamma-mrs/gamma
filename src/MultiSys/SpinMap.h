/* SpinMap.h ****************************************************-*-c++-*-
**									**
**                              G A M M A				**
**									**
**      Spin Mapping				Interface		**
**									**
**      Copyright (c) 2001 						**
**      Dr. Scott A. Smith						**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**									**
**      $Header: $
**									**
**************************************************************************

**************************************************************************
**									**
** This file contains a simple auxiliary class to support mulitple	**
** spin systems and non-mutual exchange. Specifically it handles the	**
** mapping of a spin in one system onto the spin of another system. The	**
** class has no knowledge of of spins or spin systems, so it acutally	**
** maps an index of one component into an index of a second component.	** 
**									**
** It does that by tracing 4 integers.  The first two are a component	**
** and (spin) index, as are the latter two.  The component indices	**
** may be used to index particular spin systems in a variable of type	**
** multi_sys. The (spin) indices then would index particular spins in	**
** their respective systems (or multi_sys components).			**
**									**
*************************************************************************/

#ifndef   GSpinMap_h_ 			// Is this file already included?
#  define GSpinMap_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <iostream>			// Knowledge of cout & NULL
#include <string>			// Include libstdc++ string
#include <vector>			// Include libstdc++ STL vectors
#include <Basics/ParamSet.h>		// Inlcude parameter sets

class SpinMap
  {
  public:

  int sub1;					// First component
  int sp1;					// Spin of first component
  int sub2;					// Second component
  int sp2;					// Spin of second component

// ____________________________________________________________________________
// i                      CLASS SPIN MAP ERROR HANDLING
// ____________________________________________________________________________

        // Input                eidx    : Error flag
        //                      pn	: String included in error
        //                      nr	: Flag for return (0=return)
        // Output               none    : Output process error message
        //                                Program execution stop (fatal)

         void SMerror(int eidx,                        int nr=0) const;
volatile void SMfatal(int eidx)                                  const;
         void SMerror(int eidx, const std::string& pn, int nr=0) const;
//volatile void SMfatal(int eidx, const std::string& pname)     const;

// ____________________________________________________________________________
// ii                   CLASS SPIN MAP PARAMETER SET PARSING
// ____________________________________________________________________________

/* These functions allow for a spin map to be set from parameters in a
   specified parameter set.                                                  */

        // Input                SMap    : A spin map (this)
        //                      pset    : Input parameter set
        //                      idx     : Exchange process index
        //                      mdx     : Spin mapping index
        //                      warn    : Warning flag
        // Output               TF      : Spin map is filled with
        //                                parameters read from file

bool getSMStr(const ParameterSet& pset, int idx, int mdx,
                                         std::string& sm, bool warn=true) const;

bool getSM(const ParameterSet& pset, int idx, int mdx,
         int& comp1, int& spin1, int& comp2, int& spin2, bool warn=true) const;

bool setSM(const ParameterSet& pset,int idx,int mdx,bool warn=true);

//_________________________________________________________________________________
// iii                   CLASS SPIN MAP CHECKING FUNCTIONS
//_________________________________________________________________________________

bool Check(bool warn=true) const;
bool Check(int c1, int s1, int c2, int s2, bool warn=true) const;

//_____________________________________________________________________________
// A                           SPIN MAP CONSTRUCTORS
//_____________________________________________________________________________

        // Input                spair   : A spin pairing (this)
	//			c1,s1   : 1st component, spin involved
	//			c2,s2   : 2nd component, spin involved
        // Output               void	: A spin pairing is created

MSVCDLC      SpinMap();
MSVCDLC      SpinMap(int c1, int s1, int c2, int s2);
MSVCDLC      SpinMap(const SpinMap& SM);
MSVCDLC      SpinMap(const std::string& SM);   
MSVCDLL SpinMap& operator = (const SpinMap& SM);
MSVCDLC      ~SpinMap() {};

//_____________________________________________________________________________
// B	                  SPIN MAP ACCESS FUNCTIONS
//_____________________________________________________________________________

MSVCDLL int Sub1()  const;			// Component 1 index
MSVCDLL int Sub2()  const;			// Component 2 index
MSVCDLL int Spin1() const;			// Spin 1 index
MSVCDLL int Spin2() const;			// Spin 2 index

// ____________________________________________________________________________
// C                         SPIN MAP INPUT FUNCTIONS
// ____________________________________________________________________________

        // Input                spair   : A spin pairing (this)
        //                      filename: Input filename
        //                      pset    : Input parameter set
        //                      idx     : Exchange process index
	//			mdx	: Spin mapping index
        //                      warn    : Warning output level
        //                                      0 = no warnings
        //                                      1 = warnings
        //                                     >1 = fatal warnings
        // Output               TF      : Spin map is filled with
        //                                parameters read from file
        //                                TRUE if read is successful
        // Note                         : The file should be an ASCII file
        //                                containing recognized parameters

MSVCDLL bool read(const std::string&  filename,int idx,int mdx,int warn=2);
MSVCDLL bool read(const ParameterSet& pset,    int idx,int mdx,int warn=2);

//_____________________________________________________________________________
// D		                  SPIN MAP OUTPUT
//_____________________________________________________________________________

        // Input                spair   : A spin pairing (this)
        //                      ostr	: Output stream
        // Output               none	: Spin pairing writtn to output stream
        ///F_list print			- Write system to output stream
        ///F_list <<			- Standard Output

MSVCDLL        void          print() const;
MSVCDLL        std::ostream& print(std::ostream& ostr) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const SpinMap& Sp);

};

#endif								// SpinMap.h

