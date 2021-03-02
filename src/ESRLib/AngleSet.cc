/* AngleSet.cc **************************************************-*-c++-*-
**                                                                      **
**                              G A M M A                               **
**                                                                      **
**      AngleSet 				    Implementation	**
**                                                                      **
**      Copyright (c) 1999                                              **
**      Dr. Scott A. Smith and Allie King                               **              
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Box 4005                                                        **
**      Tallahassee, FL 32306                                           **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
** Class AngleSet embodies a set of angles that define a range of       **
** orientations allowed for an object.                                  **
**                                                                      **
*************************************************************************/

#ifndef   AngleSet_cc_			// Is this file already included?
#  define AngleSet_cc_ 1		// If no, then remember it
//#  if defined(GAMPRAGMA)			// Using the GNU compiler?
//#    pragma implementation		// This is the implementation
//#endif

#include <ESRLib/AngleSet.h>		// Include class header
#include <Basics/ParamSet.h>		// Include GAMMA parameter sets
#include <Basics/Gutils.h>		// Include GAMMA errors & queries
//#include <stream.h>
#include <iostream>

// ----------------------------------------------------------------------------
// ------------------ CLASS AngleSet PRIVATE FUNCTIONS ------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                    CLASS AngleSet ERROR HANDLING
// ____________________________________________________________________________


 
	// Input		eidx		: Error Flag
	//         		noret   	: Flag for linefeed
	// Output       	none    	: Error Message Output
 

/* The following error messages use the defaults set in the Gutils package

                Case                          Error Message

                (0)                     Program Aborting.....
                (9)                     Problems During Construction
                default                 Unknown Error                        */
 

void AngleSet::ASetError(int eidx, int noret) const
  {
  string hdr("Angle Set");
  string msg;
  switch(eidx)
    {
    case 5:						 	       	// (5)
      msg = string("Cannot Construct From Parameter Set");
      GAMMAerror(hdr, msg, noret); break;
    default: GAMMAerror(hdr, eidx, noret); break;// Usually Unknown Error  (-1)
    }
  if(!noret) cout << "\n";
  }
 
 void AngleSet::ASetError(int eidx, const string pname, int noret) const
 
	// Input        		    	: AngleSet (this)
     	//              	eidx    	: Flag for error type 
	//		 	pname   	: Additional output string
     	//              	noret   	: Flag for return (0=return) 
     	// Output 		none    	: Error message

  {    
  cout << "\nAngle Set: ";
  switch(eidx)
    {
    case 40:                                      		// (40)
      cout << "Problems with File " << pname;
      break;
    case 77:						        	// (77)
      cout << "Cannot Find AngleSet Parameter" << pname << " !?";
      break;
    case 78:						        	// (78)
      cout << "Setting Angle Set Parameter " << pname << " to 0";
      break;
    case 100:                                           	// (100)
      cout << "Cannot Read Parameter " << pname;
      break;
    case 101:                                             	// (101)
      cout << "Cannot Use Parameter File " << pname;
      break;
    case 130:                                           	// (130)
      cout << "Parameter " << pname << " Is The Culprit!\n";
      break;
    default:
      cout << "Unknown error";
      break;
    }
  if(!noret) cout << ".\n";
  }  

volatile void AngleSet::ASetFatality(int eidx) const
 
	// Input                    	: AngleSet (this)
    	//       		eidx	: Error index 
     	// Output      	none    	: Error message output
     	//                          	  Program execution stopped 

  {
  ASetError(eidx, eidx);		// Output error message
  if(eidx) ASetError(0);		// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }
 
volatile void AngleSet::ASetFatality(int eidx, const string& pname) const
 
	// Input  		 		: AngleSet (this)
  	//      		eidx    	: Error index 
	//			pname   	: Additional output string
    	// Output     	none    	: Error message output
     	//                              Program execution stopped 
  
  {
  ASetError(eidx, pname, eidx);				// Output error message
  if(eidx) ASetError(0);                          // State that its fatal
  GAMMAfatal();					// Clean exit from program
  }


// ____________________________________________________________________________
// ii                CLASS AngleSet PARAMETER SET FUNCTIONS
// ____________________________________________________________________________


int AngleSet::SetAngle(const ParameterSet& pset, const string& name, int warn)

	// Input		ASet	: An angle set (this)
	//  			pset	: A parameter set
	//			name	: Parameter base name 
	// Output		TF	: True if angle parameter is read/set
	//				  from parameter in pset 
	// Note				: The Angle parameter name we'll seek
	//				  is given by "name + Ang"

  {
  string pstate;				// Temp string
  string pname = name + string("Ang");		// Set angle parameter name
  list<SinglePar>::const_iterator item;		// Iterator into pset
  item = pset.seek(pname);			// Look for parameter
  if(item == pset.end())			// If we can't find parameter
    {						// then issue errors as needed
    if(warn)
      {
      ASetError(77, pname, 1);			// Can't find X param info
      if(warn > 1)  ASetFatality(5);
      else          ASetError(78, pname);	// Setting nameX to 0
      }
    Angle = 0;					// Set default angle as 0
    return 0;					// Flag we can't get it
    }
  (*item).parse(pname,Angle,pstate);		// Set the angle value
  return 1;					// Flag all is O.K.
  }


int AngleSet::SetAverage(const ParameterSet& pset,const string& N,int warn)
 
        // Input                ASet    : An angle set (this)
        //                      pset    : A parameter set
        //                      N       : Parameter base name
        // Output               TF      : True if average parameter is read/set
        //                                from parameter in pset
        // Note                         : The Average parameter name we'll seek
        //                                is given by "name + Ave"

  {
  string pstate;				// Temp string
  string pname = N + string("Ave");		// Set average parameter name
  list<SinglePar>::const_iterator item;		// Iterator into pset
  item = pset.seek(pname);			// Look for parameter
  if(item == pset.end())			// If we can't find parameter
    {						// then issue errors as needed
    if(warn)
      {
      ASetError(77, pname, 1);			// Can't find parametger pname
      if(warn > 1)  ASetFatality(5);		// Bail if fatal error
      else          ASetError(78, pname);	// Flag setting average to 0
      }
    Average = 0;				// Set default average as 0
    return 0;					// Flag we can't get it
    }
  (*item).parse(pname,Average,pstate);		// Set the average value
  return 1;					// Flag all is O.K.
  }


int AngleSet::SetRange(const ParameterSet& pset,const string& N,int warn)
 
        // Input                ASet    : An angle set (this)
        //                      pset    : A parameter set
        //                      N       : Parameter base name
        // Output               TF      : True if range parameter is read/set
        //                                from parameter in pset
        // Note                         : The Range parameter name we'll seek
        //                                is given by "N + Del"

  {
  string pstate;				// Temp string
  string pname = N + string("Del");		// Set range parameter name
  list<SinglePar>::const_iterator item;		// Iterator into pset
  item = pset.seek(pname);			// Look for parameter
  if(item == pset.end())			// If we can't find parameter
    {						// then issue errors as needed
    if(warn)
      {
      ASetError(77, pname, 1);			// Can't find parametger pname
      if(warn > 1)  ASetFatality(5);		// Bail if fatal error
      else          ASetError(78, pname);	// Flag setting range to 0
      }
    Range = 0;					// Set default range as 0
    return 0;					// Flag we can't get it
    }
  (*item).parse(pname,Range,pstate);		// Set the range value
  return 1;					// Flag all is O.K.
  }


int AngleSet::SetStepAng(const ParameterSet& pset,const string& N,int warn)
 
        // Input                ASet    : An angle set (this)
        //                      pset    : A parameter set
        //                      N       : Parameter base name
        // Output               TF      : True if steps parameter is read/set
        //                                from parameter in pset
        // Note                         : The Range parameter name we'll seek
        //                                is given by "N + Step"

  {
  string pstate;				// Temp string
  string pname = N + string("Step");		// Set steps parameter name
  list<SinglePar>::const_iterator item;		// Iterator into pset
  item = pset.seek(pname);			// Look for parameter
  if(item == pset.end())			// If we can't find parameter
    {						// then issue errors as needed
    if(warn)
      {
      ASetError(77, pname, 1);			// Can't find parametger pname
      if(warn > 1)  ASetFatality(5);		// Bail if fatal error
      else          ASetError(78, pname);	// Flag setting steps to 0
      }
    StepAng = 0;				// Set default steps as 0
    return 0;					// Flag we can't get it
    }
  (*item).parse(pname,StepAng,pstate);		// Set the steps value
  return 1;					// Flag all is O.K.
  }

// ____________________________________________________________________________
// A                       CLASS AngleSet CONSTRUCTORS
// ____________________________________________________________________________

	
AngleSet::AngleSet(double Ang, double Ave, double Adel, double Astep) 
        
        // Input                ASet    : AngleSet (this) 
        //                      Ang     : Orientation Angle (degrees)
        //                      Ave     : Average angle (If angle range > 1)
        //                      Adel    : Angle range [Ave-Adel,Ave+Adel]
        //                      Astep   : Angle step size
        // Output                       : AngleSet (this) constructed
        // Note                         : Default values set to 0

  {
  Angle   = Ang;
  Average = Ave;
  Range   = Adel; 
  StepAng = Astep;
  }

AngleSet::AngleSet(const AngleSet& ASet1) 
         
        // Input                ASet    : AngleSet (this)
        //                      ASet1   : Another AngleSet
        // Output               Aset    : AngleSet (this) constructed
        //                                identical to ASet1

  {
  Angle   = ASet1.Angle;
  Average = ASet1.Average;
  Range   = ASet1.Range;
  StepAng = ASet1.StepAng; 
  }

void AngleSet::operator= (const AngleSet& ASet1) 

        // Input                ASet    : AngleSet (this)
        //                      ASet1   : Another AngleSet
        // Output               none    : Aset set identical to ASet1

  {
  Angle   = ASet1.Angle;
  Average = ASet1.Average;
  Range   = ASet1.Range;
  StepAng = ASet1.StepAng; 
  }

// ____________________________________________________________________________
// B                    CLASS AngleSet ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These function allow acces to the individual values within an AngleSet.   */
 
double AngleSet::GetX()	   const	{ return Angle;   }
double AngleSet::GetAv()   const	{ return Average; }
double AngleSet::GetDel()  const 	{ return Range;	  }
double AngleSet::GetStep() const	{ return StepAng; }

void AngleSet::SetX(double    Ang)	{ Angle   = Ang;   } 
void AngleSet::SetAv(double   Ave)	{ Average = Ave;   } 
void AngleSet::SetDel(double  Adel)	{ Range   = Adel;  }
void AngleSet::SetStep(double Astep)	{ StepAng = Astep; } 


// ____________________________________________________________________________
// C                     CLASS AngleSet Output FUNCTIONS
// ____________________________________________________________________________

/* These functions perform Angle Set output to an ASCII file stream

           Input                ASet            : AngleSet (this)
	  			ostr		: Output stream
	  			marg		: Margin used before output
	   Output		none		: AngleSet parameters are
	  					  sent to the output stream  */

ostream& AngleSet::printBase(ostream& ostr, const string& marg) const
  {
  ostr << marg << "Current Angle   (degrees): " << Angle;
  ostr << marg << "Average Angle   (degrees): " << Average;
  ostr << marg << "Range Fram Ave. (degrees): " << Range << " (Ave +/-)";
  ostr << marg << "Angle Step Size (degrees): " << StepAng;
  return ostr;
  }


ostream& AngleSet::print(ostream& ostr) const
  {
  ostr << "\n\n\tAngleSet Parameters\n:";
  printBase(ostr, "\n\t");
  ostr << "\n\n";
  return ostr;
  }


ostream& operator << (ostream& ostr, const AngleSet& A)
  {
  A.print(ostr);
  return ostr;
  }


// ____________________________________________________________________________
// D                      CLASS AngleSet INPUT FUNCTIONS
// ____________________________________________________________________________
 
int AngleSet::read(const string& Fin, const string& Name, int indx, int warn)

        // Input 	ASet    : AngleSet (this) 
	//		Fin	: An input (ASCII) file
	//		Name	: ASet parameter base name
	//		indx	: Parameter(s) prefix value
	//		warn   	: Warning output level
	// Output	T/F	: True if ASet filled with values read
	//                        from the file Fin

  {
  ParameterSet pset;			// Declare a parameter set
  if(!pset.read(Fin, warn?1:0))		// Read in pset from file
    {					// If we can't read in the pset
    if(warn)				// we'll issue warning and/or quit
      {					// as desired
      ASetError(40, Fin, 1);		//   Problems with file Fin
      if(warn>1) ASetFatality(101,Fin);	//   Can't read from file Fin
      else       ASetError(101,1);
      }
    return 0;  
    }   
  return read(pset, Name, indx, warn);	// Use overload if pset OK
  }

int AngleSet::read(const ParameterSet& pset,const string& N,int indx,int warn)
 
        // Input        ASet    : AngleSet (this) 
        //              pset	: A parameter set
        //              N       : ASet parameter base name
        //              indx	: Parameter(s) prefix value
        //              warn	:       Warning level
        //                              0 = no warnings
        //                              1 = non-fatal warnings
        //                              2 = fatal warnings
        // Output       T/F     : True if ASet filled with values read
        //                        from the parameter set pset
  {
  int TF =  SetAngle(pset,   N, warn);
      TF *= SetAverage(pset, N, warn);
      TF *= SetRange(pset,   N, warn);
      TF *= SetStepAng(pset,   N, warn);
  return TF;
  }

#endif							// AngleSet.cc
