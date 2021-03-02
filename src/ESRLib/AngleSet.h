/* AngleSet.h ***************************************************-*-c++-*-
**                                                                      **
**                                 G A M M A                            **
**                                                                      **
**    Angle Set			   	                Interface	**
**                                                                      **
**    Copyright (c) 1999						**
**    Scott A. Smith, Alexandra King					**
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
** Class AngleSet embodies a set of angles that define a range of 	**
** orientations allowed for an object.					**
**                                                                      **
*************************************************************************/

#ifndef _AngleSet_h_			// Is this file already included?
#define _AngleSet_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)				// Using the GNU compiler?
#    pragma interface 			// Then this is the interface
#endif

#include <string>			// Include libstdc++ strings
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets

class AngleSet			
  {
  double Angle;				//  Current angle in degrees
  double Average;			//  Average angle in degrees
  double Range;				//  Angle range in degrees 
  double StepAng;			//  Step size in degrees

private:


// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                     CLASS AngleSet ERROR HANDLING
// ____________________________________________________________________________


void AngleSet::ASetError(int eidx, int noret=0) const;
 
	// Input                ASet		: AngleSet (this)
	// 			eidx		: Error Flag
	//         		noret   	: Flag for linefeed
	// Output       	none    	: Error Message Output

void AngleSet::ASetError(int eidx, const std::string pname, int noret=0) const;
 
	// Input                ASet		: AngleSet (this)
     	//              	eidx    	: Flag for error type 
	//		 	pname   	: Additional output string
     	//              	noret   	: Flag for return (0=return) 
     	// Output 		none    	: Error message

volatile void AngleSet::ASetFatality(int eidx) const;
 
	// Input                ASet		: AngleSet (this)
    	//       		eidx    	: Error index 
     	// Output      		none    	: Error message output
     	//                      	    	  Program execution stopped

volatile void AngleSet::ASetFatality(int eidx, const std::string& pname) const;
 
	// Input                ASet		: AngleSet (this)
  	//      		eidx    	: Error index 
	//			pname   	: Additional output string
    	// Output     		none    	: Error message output
     	//                      	        Program execution stopped 

// ____________________________________________________________________________
//                   CLASS AngleSet PARAMETER SET FUNCTIONS
// ____________________________________________________________________________

/* These are facilitator functions which allow the AngleSet parameters to be
   input from a parameter set.

           Input                ASet    : An angle set (this)
                                pset    : A parameter set
                                name    : Parameter base name
           Output               TF      : True if angle set parameter is
	  				  read/set from parameter in pset   */


int AngleSet::SetAngle(const ParameterSet& pset,const std::string& name,int warn=0);
int AngleSet::SetAverage(const ParameterSet& pset,const std::string& N,int warn=0);
int AngleSet::SetRange(const ParameterSet& pset,const std::string& N,int warn=0);
int AngleSet::SetStepAng(const ParameterSet& pset,const std::string& N,int warn=0);


public:


// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                      CLASS AngleSet CONSTRUCTORS
// ____________________________________________________________________________


AngleSet::AngleSet(double Ang=0, double Ave=0, double Adel=0, double Astep=0);
	
        // Input                ASet	: AngleSet (this) 
	// 			Ang 	: Orientation Angle (degrees)
	// 			Ave	: Average angle (If angle range > 1)
	// 			Adel  	: Angle range [Ave-Adel,Ave+Adel]
	//			Astep	: Angle step size
	// Output		   	: AngleSet (this) constructed
	// Note 			: Default values set to 0

AngleSet::AngleSet(const AngleSet& ASet1); 
	
        // Input                ASet	: AngleSet (this) 
	// 			ASet1	: Another AngleSet
	// Output		Aset    : AngleSet (this) constructed
	// 	            	    	  identical to ASet1

void AngleSet::operator= (const AngleSet& ASet1);
	
        // Input                ASet	: AngleSet (this) 
	// 			ASet1  	: Another AngleSet
	// Output		none  	: Aset set identical to ASet1

// ____________________________________________________________________________
// B                     CLASS AngleSet ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions allow acces to the individual values within an AngleSet.   */

double AngleSet::GetX()    const;
double AngleSet::GetAv()   const;
double AngleSet::GetDel()  const;
double AngleSet::GetStep() const;

void AngleSet::SetX(double    Ang=0);		
void AngleSet::SetAv(double   Ave=0);	 
void AngleSet::SetDel(double  Adel=0);	
void AngleSet::SetStep(double Astep=1);

// ____________________________________________________________________________
// C                      CLASS AngleSet I/O FUNCTIONS
// ____________________________________________________________________________

/* These functions perform Angle Set output to an ASCII file stream

           Input                ASet            : AngleSet (this)
                                ostr            : Output stream
                                marg            : Margin used before output
           Output               none            : AngleSet parameters are
                                                  sent to the output stream  */

       std::ostream& AngleSet::printBase(std::ostream& ostr, const std::string& marg) const;
       std::ostream& AngleSet::print(std::ostream& ostr) const;
friend std::ostream& operator << (std::ostream& ostr, const AngleSet& ASet);

// ____________________________________________________________________________
// D                      CLASS AngleSet INPUT FUNCTIONS
// ____________________________________________________________________________
 

int AngleSet::read(const std::string& Fin,const std::string& Name,int indx=-1,int wrn=0);
 
        // Input        ASet    : AngleSet (this) 
        //              Fin     : An input (ASCII) file
        //              Name    : ASet parameter base name
        //              indx    : Parameter(s) prefix value
        //              wrn	: Warning output level
        //                      	0 = no warnings
        //                      	1 = non-fatal warnings
        //                      	2 = fatal warnings
        // Output       T/F     : True if ASet filled with values read
        //                        from the file Fin


int AngleSet::read(const ParameterSet& P,const std::string& N,int I=-1,int wrn=0);
 
        // Input        ASet    : AngleSet (this) 
        //		P	: A parameter set
        //              N	: ASet parameter base name
        //              I	: Parameter(s) prefix value
        //             	wrn	: 	Warning level
        //                      	0 = no warnings
        //                      	1 = non-fatal warnings
        //                      	2 = fatal warnings
        // Output       T/F     : True if ASet filled with values read
        //                        from the parameter set P


void AngleSet::ask_read(int argc, char* argv[], int& argn, int indx=-1);

	// Input				: AngleSet parameters
	//			argc		: Number of parameters
	//			argv		: Number of arguments
	//			argn		: Argument index
	//			indx		: AngleSet index
	// Output		void		: The parameter argn of array argc
	//					  is used to supply a filename
	//					  from which the AngleSet parameters
	//					  are read. If the argument argn is 
	//					  not in argv, the user is asked for a 
	//					  filename.
	// Note			 	: The file shoud be an ASCII file
	//					  containing recognized sys parameters
	// Note				: AngleSet parameters are modifed (filled)



  };

#endif						// AngleSet.h
