/* Fibre.cc *****************************************************-*-c++-*-
**                                                                      **
**                                 G A M M A                            **
**                                                                      **
**    Fibre                                         Implementation 	**
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
** Class Fibre embodies a myocin fibre containing nitroxide.            **
**                                                                      **
*************************************************************************/

#ifndef   Fibre_cc_			// Is this file already included?
#  define Fibre_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <Basics/ParamSet.h>			// Include GAMMA parameters
#include <Basics/Gutils.h>			// Include GAMMA queries,errors
#include <ESRLib/Fibre.h>			// Include class interface

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

//____________________________________________________________________________
// I                    CLASS Fibre ERROR HANDLING
//____________________________________________________________________________


void Fibre::FibError(int eidx, int noret) const
 
	// Input	eidx		: Error Flag
	//         	noret		: Flag for linefeed
	// Output	none		: Error Message Output
	// Note 			: Standard messages in Gutils used

/*  0: Program Aborting.......
    9: Problems During Construction                                         */
 
  {
  string hdr("Myocin Fiber");
  switch(eidx)
    {
    case 0:  GAMMAerror(hdr, 0, noret);  break;	// Program Aborting	   (0)
    default: GAMMAerror(hdr,eidx,noret); break;	// Use Default Message
    }
  if(!noret) cout << ".\n";
  }
 
void Fibre::FibError(int eidx, const string& pname, int noret) const
 
	// Input	*this	    	: a Fibre
     	// 		eidx    	: Flag for error type 
	//		pname   	: Additional output string
     	//		noret   	: Flag for return (0=return) 
     	// Output	none    	: Error message

/*  1: Problems With File pname
    2: Cannot Read Parameter pname                                          */

  {    
  string hdr("Myocin Fiber");
  string msg;
  switch(eidx)
    {
    case 78:					        	// (78)
      msg = string("Setting Fibre Parameter ") + pname + string(" to 0");
      GAMMAerror(hdr, msg, noret); 	break;
    case 101:                                             	// (101)
      msg = string("Can't Use Parameter File ") + pname;
      GAMMAerror(hdr, msg, noret); 	break;
    case 130:                                           	// (130)
      msg = string("Parameter ") + pname + string(" Is The Culprit!");
      GAMMAerror(hdr, msg, noret); 	break;
    default: GAMMAerror(hdr, eidx, pname, noret); break;
    }
  if(!noret) cout << ".\n";
  }  

volatile void Fibre::FibFatality(int eidx) const
 
	// Input	*this		: a Fibre
    	// 		eidx    	: Error index 
     	// Output	none    	: Error message output
     	// 				  Program execution stopped 

  {
  FibError(eidx, eidx);			// Output error message
  if(eidx) FibError(0);			// State that its fatal
  GAMMAfatal();					// Clean exit from program
  }
 
volatile void Fibre::FibFatality(int eidx, const string& pname) const
 
	// Input 	*this		: a Fibre
  	//      	eidx    	: Error index 
	//		pname   	: Additional output string
    	// Output	none    	: Error message output
     	//				  Program execution stopped 
  
  {
  FibError(eidx, pname, eidx);		// Output error message
  if(eidx) FibError(0);			// State that its fatal
  GAMMAfatal();					// Clean exit from programam
  }
 
// ____________________________________________________________________________
// ii                  CLASS Fiber PARAMETER SET FUNCTIONS
// ____________________________________________________________________________
 

int Fibre::SetKernels(const ParameterSet& pset, int warn)
 
        // Input                ASet    : An angle set (this)
        //                      pset    : A parameter set
        //			warn    : Warning level
        //                              0 = no warnings
        //                              1 = non-fatal warnings
        //                              2 = fatal warnings
        // Output               TF      : True if kernel values read/set
        //                                from parameters in pset
        // Note                         : The kernel parameter names we'll seek
        //                                are "" and ""

  {
  string pstate;                                // Temp string
  string pname("AlphaKernel");			// Set 1st kernel pname
  list<SinglePar>::const_iterator item;         // Iterator into pset
  item = pset.seek(pname);                      // Look for parameter
  if(item == pset.end())                        // If we can't find parameter
    {                                           // then issue errors as needed
    if(warn)
      {  
      FibError(77, pname, 1);                  // Can't find parametger pname
      if(warn > 1)  FibFatality(5);            // Bail if fatal error
      else          FibError(78, pname);       // Flag setting kernel to 1
      }  
    AlphaKernel=1;				// Set default kernel as 1
    return 0;                                   // Flag we can't get it
    }
  (*item).parse(pname,AlphaKernel,pstate);	// Set the 1st kernel value
  pname = string("BetaKernel");			// Set 2nd kernel pname
  item = pset.seek(pname);                      // Look for parameter
  if(item == pset.end())                        // If we can't find parameter
    {                                           // then issue errors as needed
    if(warn)
      {  
      FibError(77, pname, 1);                  // Can't find parametger pname
      if(warn > 1)  FibFatality(5);            // Bail if fatal error
      else          FibError(78, pname);       // Flag setting kernel to 1
      }  
    BetaKernel=1;				// Set default kernel as 1
    return 0;                                   // Flag we can't get it
    }
  (*item).parse(pname,BetaKernel,pstate);	// Set the 2nd kernel value
  return 1;                                     // Flag all is O.K.
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------
 
// ____________________________________________________________________________
// A                  CLASS Fibre CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

Fibre::Fibre()

  {
  AlphaKernel = 1;
  BetaKernel  = 1;
  }
 

Fibre::Fibre(const Fibre& Fib1)

        // Input                Fib	: Fiber (this)
        //                      Fib1	: Another Fiber
        // Output               Fib	: Fiber (this) constructed
        //                                identical to Fib1
 
  {
  AlphaMol    = Fib1.AlphaMol;
  BetaMol     = Fib1.BetaMol;
  GammaMol    = Fib1.GammaMol;
  AlphaSam    = Fib1.AlphaSam;
  BetaSam     = Fib1.BetaSam;
  GammaSam    = Fib1.GammaSam;
  AlphaLab    = Fib1.AlphaLab;
  BetaLab     = Fib1.BetaLab;
  GammaLab    = Fib1.GammaLab;
  AlphaKernel = Fib1.AlphaKernel;
  BetaKernel  = Fib1.BetaKernel;
  }
 

void Fibre::operator= (const Fibre& Fib1)
 
        // Input                Fib	: Fiber (this)
        //                      Fib1	: Another Fiber
        // Output               none    : Fib set identical to Fib1
 
  {
  AlphaMol    = Fib1.AlphaMol;
  BetaMol     = Fib1.BetaMol;
  GammaMol    = Fib1.GammaMol;
  AlphaSam    = Fib1.AlphaSam;
  BetaSam     = Fib1.BetaSam;
  GammaSam    = Fib1.GammaSam;
  AlphaLab    = Fib1.AlphaLab;
  BetaLab     = Fib1.BetaLab;
  GammaLab    = Fib1.GammaLab;
  AlphaKernel = Fib1.AlphaKernel;
  BetaKernel  = Fib1.BetaKernel;
  }

// ____________________________________________________________________________
// C                      CLASS Fibre INPUT FUNCTIONS
// ____________________________________________________________________________

int Fibre::read(const string& filein, int indx, int warn)

	// Input	*this	: A Fibre
	//		filein	: An input (ASCII) file
	//		indx	: Point index
	//		warn	: Warning output level
	// Output	none	: Fibre filled with parameters specified in
	//			  filein
  {
  ParameterSet pset;				// Declare a parameter set
  if(!pset.read(filein, warn?1:0))		// Read in pset from file
    {
    if(warn)
      {
      FibError(1, filein, 1);			// Problems with file filein
      if(warn>1) FibFatality(101,filein);	// Can't read from file filein
      else FibError(101,filein,1);
      }
    return 0;
    }
  return read(pset,indx,warn);
  }

int Fibre::read(const ParameterSet& pset, int indx, int warn)

	// Input	*this	: A Fibre
	//  		pset	: A parameter set
	//		indx	: Fibre index
	// Output	none	: Fibre read from parameters
 	//			  in pset with index indx.
	// Note			: If parameter not found it is set to 0 and
 	//			  error message is sent to standard output.

  {
  int TF = 1;
  TF *= AlphaMol.read(pset, "AlphaMol", indx, (warn)?1:0);
  TF *= BetaMol.read(pset,  "BetaMol",  indx, (warn)?1:0);
  TF *= GammaMol.read(pset, "GammaMol", indx, (warn)?1:0);
  TF *= AlphaSam.read(pset, "AlphaSam", indx, (warn)?1:0);
  TF *= BetaSam.read(pset,  "BetaSam",  indx, (warn)?1:0);
  TF *= GammaSam.read(pset, "GammaSam", indx, (warn)?1:0);
  TF *= AlphaLab.read(pset, "AlphaLab", indx, (warn)?1:0);
  TF *= BetaLab .read(pset, "BetaLab",  indx, (warn)?1:0);
  TF *= GammaLab.read(pset, "GammaLab", indx, (warn)?1:0);
  TF *= SetKernels(pset, (warn)?1:0);
  return TF;
  }

// ____________________________________________________________________________
// D                    Class Fibre Output Functions
// ____________________________________________________________________________


ostream& Fibre::print(ostream& ostr) const
 
        // Input             Fib	: A fiber (this)
        //                   ostr 	: An output stream
        // Output            ostr 	: The output stream, modified to
        //                          	  contain the information about
        //                          	  the fiber

  {
  ostr << "\n\tWe Gots Fibre";
  ostr << AlphaMol;
  ostr << BetaMol;
  ostr << GammaMol;
  ostr << AlphaSam;
  ostr << BetaSam;
  ostr << GammaSam;
  ostr << AlphaLab;
  ostr << BetaLab;
  ostr << GammaLab;
  ostr << AlphaKernel;
  ostr << BetaKernel;
  return ostr;
  }

 
ostream& operator<< (ostream& ostr, const Fibre& Fib) {return Fib.print(ostr);}

        // Input             ostr 	: An output stream ostr
        //                   Fib	: A fiber
        // Output            ostr 	: The output stream, modified to
        //                          	  contain the information about
        //                          	  the fiber

// ____________________________________________________________________________
//                       CLASS Angle_Set ACCESS FUNCTIONS
// ____________________________________________________________________________

/*
double Fibre::GetAlphaMolX()			{ return AlphaMolX; }			
void   Fibre::SetAlphaMolX(double A=0)	{ AlphaMolX = A; }

double Fibre::GetAlphaMolAv()			{ return AlphaMolAv; }
void   Fibre::SetAlphaMolAv(double A=0)	{ AlphaMolAv = A; }

double Fibre::GetAlphaMolDel()		{ return AlphaMolDel; }			
void   Fibre::SetAlphaMolDel(double A=0)	{ AlphaMolDel = A; }

double Fibre::GetAlphaMolStep()		{ return AlphaMolStep; }			
void   Fibre::SetAlphaMolStep(double A=0)	{ AlphaMolStep = A; }

double Fibre::GetBetaMolX()			{ return BetaMolX; }			
void   Fibre::SetBetaMolX(double A=0)	{ BetaMolX = A; }

double Fibre::GetBetaMolAv()			{ return BetaMolAv; }			
void   Fibre::SetBetaMolAv(double A=0)	{ BetaMolAv = A; }

double Fibre::GetBetaMolDel()			{ return BetaMolDel; }			
void   Fibre::SetBetaMolDel(double A=0)	{ BetaMolDel = A; }

double Fibre::GetBetaMolStep()		{ return BetaMolStep; }			
void   Fibre::SetBetaMolStep(double A=0)	{ BetaMolStep = A; }	

double Fibre::GetGammaMolX()			{ return GammaMolX; }			
void   Fibre::SetGammaMolX(double A=0)	{ GammaMolX = A; }

double Fibre::GetGammaMolAv()			{ return GammaMolAv; }			
void   Fibre::SetGammaMolAv(double A=0)	{ GammaMolAv = A; }

double Fibre::GetGammaMolDel()		{ return GammaMolDel; }			
void   Fibre::SetGammaMolDel(double A=0)	{ GammaMolDel = A; }

double Fibre::GetGammaMolStep()		{ return GammaMolStep; }			
void   Fibre::SetGammaMolStep(double A=0)	{ GammaMolStep = A; }	

double Fibre::GetAlphaSamX()			{ return AlphaSamX; }			
void   Fibre::SetAlphaSamX(double A=0)	{ AlphaSamX = A; }

double Fibre::GetAlphaSamAv()			{ return AlphaSamAv; }			
void   Fibre::SetAlphaSamAv(double A=0)	{ AlphaSamAv = A; }

double Fibre::GetAlphaSamDel()		{ return AlphaSamDel; }			
void   Fibre::SetAlphaSamDel(double A=0)	{ AlphaSamDel = A; }

double Fibre::GetAlphaSamStep()		{ return AlphaSamStep; }			
void   Fibre::SetAlphaSamStep(double A=0)	{ AlphaSamStep = A; }

double Fibre::GetBetaSamX()			{ return BetaSamX; }			
void   Fibre::SetBetaSamX(double A=0)	{ BetaSamX = A; }

double Fibre::GetBetaSamAv()			{ return BetaSamAv; }			
void   Fibre::SetBetaSamAv(double A=0)	{ BetaSamAv = A; }

double Fibre::GetBetaSamDel()			{ return BetaSamDel; }			
void   Fibre::SetBetaSamDel(double A=0)	{ BetaSamDel = A; }

double Fibre::GetBetaSamStep()		{ return BetaSamStep; }			
void   Fibre::SetBetaSamStep(double A=0)	{ BetaSamStep = A; }	

double Fibre::GetGammaSamX()			{ return GammaSamX; }			
void   Fibre::SetGammaSamX(double A=0)	{ GammaSamX = A; }

double Fibre::GetGammaSamAv()			{ return GammaSamAv; }			
void   Fibre::SetGammaSamAv(double A=0)	{ GammaSamAv = A; }

double Fibre::GetGammaSamDel()		{ return GammaSamDel; }			
void   Fibre::SetGammaSamDel(double A=0)	{ GammaSamDel = A; }

double Fibre::GetGammaSamStep()		{ return GammaSamStep; }			
void   Fibre::SetGammaSamStep(double A=0)	{ GammaSamStep = A; }

double Fibre::GetAlphaLabX()			{ return AlphaLabX; }			
void   Fibre::SetAlphaLabX(double A=0)	{ AlphaLabX = A; }

double Fibre::GetAlphaLabAv()			{ return AlphaLabAv; }			
void   Fibre::SetAlphaLabAv(double A=0)	{ AlphaLabAv = A; }

double Fibre::GetAlphaLabDel()		{ return AlphaLabDel; }			
void   Fibre::SetAlphaLabDel(double A=0)	{ AlphaLabDel = A; }

double Fibre::GetAlphaLabStep()		{ return AlphaLabStep; }			
void   Fibre::SetAlphaLabStep(double A=0)	{ AlphaLabStep = A; }

double Fibre::GetBetaLabX()			{ return BetaLabX; }			
void   Fibre::SetBetaLabX(double A=0)	{ BetaLabX = A; }

double Fibre::GetBetaLabAv()			{ return BetaLabAv; }			
void   Fibre::SetBetaLabAv(double A=0)	{ BetaLabAv = A; }

double Fibre::GetBetaLabDel()			{ return BetaLabDel; }			
void   Fibre::SetBetaLabDel(double A=0)	{ BetaLabDel = A; }

double Fibre::GetBetaLabStep()		{ return BetaLabStep; }			
void   Fibre::SetBetaLabStep(double A=0)	{ BetaLabStep = A; }	

double Fibre::GetGammaLabX()			{ return GammaLabX; }			
void   Fibre::SetGammaLabX(double A=0)	{ GammaLabX = A; }

double Fibre::GetGammaLabAv()			{ return GammaLabAv; }			
void   Fibre::SetGammaLabAv(double A=0)	{ GammaLabAv = A; }

double Fibre::GetGammaLabDel()		{ return GammaLabDel; }			
void   Fibre::SetGammaLabDel(double A=0)	{ GammaLabDel = A; }

double Fibre::GetGammaLabStep()		{ return GammaLabStep; }			
void   Fibre::SetGammaLabStep(double A=0)	{ GammaLabStep = A; }

double Fibre::GetAlphaKernel()		{ return AlphaKernel; }			
void   Fibre::SetAlphaKernel(double A=0)	{ AlphaKernel = A; }

double Fibre::GetBetaKernel()			{ return BetaKernel; }			
void   Fibre::SetBetaKernel(double A=0)	{ BetaKernel = A; }
*/

#endif								// Fibre.cc

