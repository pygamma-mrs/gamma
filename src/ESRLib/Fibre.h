/* Fibre.h ******************************************************-*-c++-*-
**                                                                      **
**                                 G A M M A                            **
**                                                                      **
**    Fibre						Interface       **
**                                                                      **
**    Copyright (c) 1999                                                **
**    Scott A. Smith, Alexandra King                                    **
**    National High Magnetic Field Laboratory                           **
**    1800 E. Paul Dirac Drive                                          **
**    Box 4005                                                          **
**    Tallahassee, FL 32306                                             **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
** Class Fibre embodies a myocin fibre containing nitroxide.		**
**                                                                      **
*************************************************************************/

#ifndef _Fibre_h_			// Is this file already included?
#define _Fibre_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#endif

#include <string>				// Include libstdc++ strings
#include <ESRLib/AngleSet.h>			// Include angle sets

class Fibre
  {
  AngleSet AlphaMol, BetaMol, GammaMol;		// Molecular frame angles
  AngleSet AlphaSam, BetaSam, GammaSam;
  AngleSet AlphaLab, BetaLab, GammaLab;		// Laboratory frame angles
  double AlphaKernel, BetaKernel;

private:

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

//____________________________________________________________________________
// i                    CLASS Fibre ERROR HANDLING
//____________________________________________________________________________


void Fibre::FibError(int eidx, int noret=0) const;
 
	// Input	eidx		: Error Flag
	//         	noret		: Flag for linefeed
	// Output	none		: Error Message Output

void Fibre::FibError(int eidx, const std::string& pname, int noret=0) const;
 
	// Input	*this	    	: A Fibre
     	// 		eidx    	: Flag for error type 
	//		pname   	: Additional output std::string
     	//		noret   	: Flag for return (0=return) 
     	// Output	none    	: Error message

volatile void Fibre::FibFatality(int eidx) const;
 
	// Input	*this		: A Fibre
    	// 		eidx    	: Error index 
     	// Output	none    	: Error message output
     	// 				  Program execution stopped 

volatile void Fibre::FibFatality(int eidx, const std::string& pname) const;
 
	// Input 	*this		: A Fibre
  	//      	eidx    	: Error index 
	//		pname   	: Additional output string
    	// Output	none    	: Error message output
     	//				  Program execution stopped 
// ____________________________________________________________________________
// ii                CLASS AngleSet PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

 
int Fibre::SetKernels(const ParameterSet& pset, int warn=1); 
  
        // Input                ASet    : An angle set (this)
        //                      pset    : A parameter set 
        //                      warn    : Warning level
        //                              0 = no warnings 
        //                              1 = non-fatal warnings
        //                              2 = fatal warnings 
        // Output               TF      : True if kernel values read/set
        //                                from parameters in pset    
        // Note                         : The kernel parameter names we'll seek 
        //                                are "" and ""

public:

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A                  CLASS Fibre CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 

Fibre::Fibre();

        // Input                Fib     : Fiber (this)
	// Output		none	: Fibre constructed (default)
 

Fibre::Fibre(const Fibre& Fib1);
 
        // Input                Fib     : Fiber (this)
        //                      Fib1    : Another Fiber
        // Output               Fib     : Fiber (this) constructed
        //                                identical to Fib1
  
 
void Fibre::operator= (const Fibre& Fib1);
 
        // Input                Fib     : Fiber (this)
        //                      Fib1    : Another Fiber
        // Output               none    : Fib set identical to Fib1 


// ____________________________________________________________________________
// C                      CLASS Fibre INPUT FUNCTIONS
// ____________________________________________________________________________

int Fibre::read(const std::string& filein, int indx=-1, int warn=0);

	// Input	*this	: A Fibre
	//		filein	: An input (ASCII) file
	//		indx	: Point index
	//		warn	: Warning output level
	// Output	none	: Fibre filled with parameters specified
	//			  in the file filein

int Fibre::read(const ParameterSet& pset, int indx=-1, int warn=0);

	// Input	*this	: A Fibre
	//  		pset	: A parameter set
	//		indx	: Fibre index
	// Output	none	: Fibre read from parameters in pset
	//			  with index indx.
	// Note			: If parameter not found it is set to 0
	//			  & error message is sent to std output.

// ____________________________________________________________________________
// D                    Class Fibre Output Functions
// ____________________________________________________________________________


std::ostream& Fibre::print(std::ostream& ostr) const;
 
        // Input             Fib        : A fiber (this)
        //                   ostr       : An output stream
        // Output            ostr       : The output stream, modified to
        //                                contain the information about
        //                                the fiber
 
friend std::ostream& operator<< (std::ostream& ostr, const Fibre& Fib);
 
        // Input             ostr       : An output stream ostr
        //                   Fib        : A fiber
        // Output            ostr       : The output stream, modified to 
        //                                contain the information about
        //                                the fiber

// ____________________________________________________________________________
//                       CLASS Angle_Set ACCESS FUNCTIONS
// ____________________________________________________________________________

  };

#endif							// Fibre.h
