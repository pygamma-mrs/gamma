/* SinglePar.cc *************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**      Single Parameter                           Implementation	**
**									**
**      Copyright (c) 1994						**
**      Scott A. Smith							**
**									**
**      $Header: $
**									**
*************************************************************************/
     
/*************************************************************************
**									**
**  Description								**
**									**
**  Class SinglePar embodies a single parameter that may be used to	**
**  export/import values to/from ASCII files.  They also form the entry	**
**  type for GAMMA parameter lists which perform the same for multiple	**
**  parameters. Access functions are provided to get at the individual	**
**  entities in the parameter.						**
**									**
*************************************************************************/

#ifndef   SinglePar_cc_			// Is this file already included?
#  define SinglePar_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <Basics/SinglePar.h>		// Include the interface
#include <GamGen.h>			// Include OS specifics
#include <Basics/StringCut.h>		// Include strings & string cutting
#include <Basics/Gutils.h>		// Include GAMMA errors
#include <cstdlib>			// Include functions atoi, atof, etc.
#include <fstream>			// Include file streams
#include <string>			// Include libstdc++ strings
#include <iostream>			// Include input output streams (cout)
#include <cstring>

using std::string;			// Using libstdc++ strings
using std::vector;			// Using libstdc++ vectors

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i                    SINGLE PARAMETER ERROR HANDLING
// ____________________________________________________________________________

        // Input		par	: A parameter (this)
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
 
void SinglePar::SPerror(int eidx, int noret) const
  {
  string hdr("Isotope");
  string msg;
  switch(eidx)
    {
    case 1: msg = string("Typecasting Input To Integer Parameter");	// (1)
            GAMMAerror(hdr, eidx, noret); break;
    case 2: msg = string("Typecasting Input To Double Parameter");		// (2)
            GAMMAerror(hdr,msg,noret); break;
    case 3: msg = string("Typecasting Input To String Parameter");		// (3)
             GAMMAerror(hdr,msg,noret); break;
    default: GAMMAerror(hdr, eidx, noret); break;
    }
  }    

volatile void SinglePar::SPfatality(int eidx) const
  {                                                                 
  SPerror(eidx, 1);				// Output error message
  if(eidx) SPerror(0);                          // State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________ 
// ii                       TYPE MORPHING FUNCTIONS
// ____________________________________________________________________________

        // Input              par : A single parameter (this)
        //                   Type : The parameter type
        //                   Name : The parameter name
        //                   State: The parameter statement
        // Output            void : The parameter name is changed
        // Note                   : These are a little finicky as they
        //                          enforce some parameter types to have        
        //                          reasonable data values 

// Additional functions name & state are public and found in file section B

void SinglePar::type(int Type)
  {
  if(Type == ParType) return;		// Do nothing if no type change
  switch(Type)				// Handle change on casewise
    {
    case -1:				// A Null parameter
      ParName  = string("");	// so we need to have nothing else
      ParData  = ParName;
      ParState = ParName;
      ParType = Type;
      break;
    case 0:				// An integer parameter
      ParData = cutInt(ParData);	// Insure data is integer
      if(!ParData.length())
        ParData = string("0");
      ParType = Type;
      break;
    case 1:				// A double parameter
      ParData = cutDouble(ParData);	// Insure data is double
      if(!ParData.length())
        ParData = string("0.0");
      ParType = Type;
      break;
    default:				// We don't deal with morphing
      ParType = Type;			// another other more complicated
      break; 				// or user defined parameter types
    }
  }


void SinglePar::data(string Data)
  {
  switch(ParType)			// Handle change on casewise
    {
    case -1:				// Null parameter, no value allowed
      break;				// any value for this
    case 0:				// An integer parameter
      ParData = cutWhite(Data);		// Insure data is integer
      if(!ParData.length())
        ParData = string("0");
      break;
    case 1:				// A double parameter
      ParData = cutDouble(Data);	// Insure data is double
      if(!ParData.length())
        ParData = string("0.0");
      break;
    default:				// We don't deal with morphing
      ParData = Data;			// another other more complicated
      break; 				// or user defined parameter types
    }
  }
 
// ____________________________________________________________________________
// iii                SINGLE PARAMETER INPUT FUNCTIONS
// ____________________________________________________________________________

/* These Are Functions That Facilitate Input Of Parameters From ASCII Files.
   They Could Be Abused So They Are Kept Private Here.                       */

	// Input	    par : A single parameter (this)
	//		    line: String representing part of ASCII input line
	// Return	    stat: TRUE if coord parameter set from line
 
int SinglePar::setCoord(string& input)
  {
  cutBlksXBlks(input, "(");			// Cut out past ( 
  string scratch = cutDouble(input);		// Cut out a double
  if(scratch == string("")) return -1;			// Bail if no 1st ordinate
  ParData = string("( ") + scratch;		// Start data with ( #
  cutBlksXBlks(input, ",");			// Remove blanks & comma
  scratch = cutDouble(input);			// Cut out a double
  if(scratch == string("")) return -1;			// Bail if no 2nd ordinate
  ParData += string(", ") + scratch;		// Data now ( #, #
  cutBlksXBlks(input, ",");			// Remove blanks & comma
  scratch = cutDouble(input);			// Cut out a double
  if(scratch == string("")) return -1;			// Bail if no 3nd ordinate
  ParData += string(", ") + scratch 		// Data now ( #, #, # ) 
          +  string(" )");
  cutBlksXBlks(input, ")");			// Cut out past closing )
  cutBlksXBlks(input, "-");			// Cut remaining past -
  ParState = input;				// Rest is the comment
  return 1;
  }



	// Input	    par : A single parameter (this)
	//		    inp : Input file stream
	//		    line: String representing part of ASCII input line
	// Return	     TF	: TRUE if tensor parameter set from line
 
// sosi - still need to remove [#] on 2nd & 3rd lines when multi-sys?
int SinglePar::setTensor(std::ifstream& inp, string& input)
  {
  ParData = cutInt(input);			// Cut int, tensor rank
  if(ParData == string("")) return 0;			// Exit if no rank set
  int rank = atoi(ParData.c_str());		// This is the tensor rank
  if(rank < 0) return 0;			// No negative ranked tensors
  cutBlksXBlks(input, "-");			// Cut remaining past -
  ParState = input;				// Rest is the comment
  string scratch;				// To parse 2nd & 3rd lines
  int retval = 1;				// Default return value
  switch(rank)					//	Switch read rank based
    {
    case 0:					// Rank 0 spatial tensor
    case 1:					// Rank 1 spatial tensor
      retval = 0;
      break;
    case 2:					// Rank 2 spatial tensor
      {
      char buf[200];				// Need a character buffer
      inp.getline(buf, 200);			// Read next (2nd) line
      input = string(buf);			//   convert char to string
//      getline(inp, input);			// Read next (2nd) line
      if(inp.eof()) return -1;			//   Exit if no data
      cutWhite(input);				//   Remove starting blanks
      SinglePar P2;				//   Use temp parameter
      if(!P2.setCoord(input)) return -1;	//   Bail if no tensor info
      ParData += string(", ")+P2.data();	//   Add in next data
      inp.getline(buf, 200);			// Read next (2nd) line
      input = string(buf);			//   convert char to string
//      getline(inp, input);			// Read next (3rd) line
      SinglePar P3;				//   Use temp parameter
      if(!P3.setCoord(input)) return -1;	//   Bail if no tensor info
      ParData += string(", ")+P3.data();	//   Add in next data
      break;
      }
    default:
      retval = -1;
      break;
    }
  return retval;
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A              SINGLE PARAMETER CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

SinglePar::SinglePar() : ParType(-1) {}

SinglePar::SinglePar(const SinglePar& par)
          : ParName(par.ParName), ParType(par.ParType),
            ParData(par.ParData), ParState(par.ParState) {}

// ----------------------------------------------------------------------------
//                Constructors That Set Up The Parameter
// ----------------------------------------------------------------------------
 
SinglePar::SinglePar(const string& pn, int pt, const string& pd, const string& ps)
          : ParName(pn), ParType(pt), ParData(pd),             ParState(ps) {}

SinglePar::SinglePar(const string& pn,                   int pd, const string& ps)
          : ParName(pn), ParType(0),  ParData(Gdec(pd)),       ParState(ps) {}
  
SinglePar::SinglePar(const string& pn,                double pd, const string& ps)
          : ParName(pn), ParType(1),  ParData(Gform("%g",pd)), ParState(ps) {}
  
SinglePar::SinglePar(const string& pn,         const string& pd, const string& ps)
          : ParName(pn), ParType(2),  ParData(pd),             ParState(ps) {}

// ----------------------------------------------------------------------------
//                Odd Constructors, Destructor, Assignment
// ----------------------------------------------------------------------------

SinglePar::SinglePar(const string& pname) : ParName(pname) {}
 
SinglePar::~SinglePar() {}

SinglePar& SinglePar::operator=(const SinglePar& par)
  { 
  if(this == &par) return *this;	// Forget self assignment
  ParName  = par.ParName;		// Copy contents of par into us
  ParType  = par.ParType,
  ParData  = par.ParData;
  ParState = par.ParState;
  return *this;				// Return ourself
  }

// ____________________________________________________________________________
// B                  SINGLE PARAMETER ACCESS FUNCTIONS
// ____________________________________________________________________________

/*         Input               par : A single parameter (this)                  
           Output (name)      name : The name of parameter par
                  (type)      type : Integer designation of type
                  (data)      data : string containing parameter data
                  (state)    state : string containing parameter comment
	   Note			   : These are inline & reside in the header */

const string& SinglePar::name() const  { return ParName; }
      int     SinglePar::type() const  { return ParType; }
const string& SinglePar::data() const  { return ParData; }
const string& SinglePar::state() const { return ParState; } 
 

/*         Input              par : A single parameter (this)
                             Name : The parameter name
                             State: The parameter statement
           Output            void : The parameter name is changed           */

void SinglePar::name(const  string& N) { if(ParType>=0) ParName  = N; }
void SinglePar::state(const string& S) { if(ParType>=0) ParState = S; }

// Additional functions type & data are private and found in file section ii 
 
// ____________________________________________________________________________
// C                      SINGLE PARAMETER PARSING FUNCTIONS
// ____________________________________________________________________________

/* These functions take a generic parameter and attempt to split it apart 
   into components as indicated in the function argument list.  They return
   TF depending upon whether the parsing was sucessful.                      */

        // Input              par : A parameter (this)
        //                   name : string for parameter name
	//		    state : Parameter statement
	// Output	     none : Name and statement set from par   

void SinglePar::GetNS(string& name, string& state) const
  {
  name  = ParName;			// Set the parameter name
  state = ParState;			// Set the parameter comment
  }

        // Input              par : A parameter (this)
        //                   name : string for parameter name
	//		        i : Integer parameter value
	//		     flag : Warning flag if type cast needed
	//		    state : Parameter statement
        // Output              TF : Returns TRUE if the parameter is
	//			    properly parsed as an integer
	// Note			  : name, i and state are set

bool SinglePar::parse(string& name, int& i, string& state, int warn) const
  {
  GetNS(name, state);			// Set name and statement
  bool TF = true;			// Assume we can parse OK
  if(ParType != 0)			// Insure par is integer
    {					// If not we'll need to convert
    TF = false;
    if(warn) SPerror(1);
    }
  i = atoi(ParData.c_str());
  return TF;
  }

        // Input              par : A parameter (this)
        //                   name : string for parameter name
	//		        d : Double parameter value
	//		     flag : Warning flag if type cast needed
	//		    state : Parameter statement
        // Output              TF : Returns TRUE if the parameter is
	//			    properly parsed as a double
	// Note			  : name, d and state are set

bool SinglePar::parse(string& name, double& d, string& state, int warn) const
  {
  GetNS(name, state);			// Set name and statement
  bool TF = true;			// Assume we can parse OK
  if(ParType != 1)			// Insure parameter is double
    {					// If not we'll need to convert
    TF = false;
    if(warn) SPerror(2);
    }
  d = atof(ParData.c_str());
  return TF;
  }

        // Input              par : A parameter (this)
        //                   name : string for parameter name
	//		        d : Double parameter value
	//		     flag : Warning flag if type cast needed
	//		    state : Parameter statement
        // Output              TF : Returns TRUE if the parameter is
	//			    properly parsed as a double
	// Note			  : name, val and state are set

bool SinglePar::parse(string& name, string& val, string& state, int warn) const
  {
  GetNS(name, state);			// Set name and statement
  bool TF = true;
  if(ParType != 2)
    {
    TF = false;
    if(warn) SPerror(3);
    }
  val = ParData;
  return TF;
  }

        // Input              par : A parameter (this)
        //                   name : string for parameter name
	//		 dx,dy,dz : Double parameter values
	//		     flag : Warning flag if type cast needed
	//		    state : Parameter statement
        // Output              TF : Returns TRUE if the parameter is
	//			    properly parsed as a coordinate
	// Note			  : name, dx, dy, dz, and state are set

bool SinglePar::parse(string& name, double& dx, double& dy, double& dz,
                                       string& state, int warn) const
  {
  bool TF = true;
  GetNS(name, state);				// Set name and statement
  if(ParType != 3)
    {
    TF = false;
    if(warn) std::cout << "\nTypecasting";
    }
  string val1 = ParData;			// Set val to all of data
  cutParBlks(val1);				// Remove leading ( + blanks
  string val2 = cutDouble(val1);		// Clip a double from front
  dx = atof(val2.c_str());			// Set it to a double, dx
  cutBlksXBlks(val1, ",");			// Clip off blank , blank
  val2 = cutDouble(val1);			// Clip off another double
  dy = atof(val2.c_str());			// Set it to a double, dy
  cutBlksXBlks(val1, ",");			// Clip off blank , blank
  val2 = cutDouble(val1);			// Clip off the third ordinate
  dz = atof(val2.c_str());			// Set it to a double, dz
  return TF;					// Now we're done
  }


        // Input              par : A parameter (this)
        //                   name : string for parameter name
	//		     rank : An integer value
	//         diso,delz,deta : Double parameter values
	//       alpha,beta,gamma : Double parameter values
	//		     flag : Warning flag if type cast needed
	//		    state : Parameter statement
        // Output              TF : Returns TRUE if the parameter is
	//			    properly parsed as a tensor
	// Note			  : name, rank diso, delz, deta,
	//			    alpha, beta, gamma, and state are set

bool SinglePar::parse(string& name, int& rank, double& diso, double& delz,
             double& deta, double& alpha, double& beta, double& gamma,
                                          string& state, int warn) const
  {
  bool TF = true;
  name = ParName;
  state = ParState;
  if(ParType != 4)
    {
    TF = false;
    if(warn) std::cout << "\nTypecasting";
    }
  string val1, val2;
  val1 = ParData;
  cutParBlks(val1);
  val2 = cutInt(val1);
  rank = atoi(val2.c_str());
  cutParBlks(val1);
  val2 = cutDouble(val1);
  diso = atof(val2.c_str());
  cutBlksXBlks(val1, ",");
  val2 = cutDouble(val1);
  delz = atof(val2.c_str());
  cutBlksXBlks(val1, ",");
  val2 = cutDouble(val1);
  deta = atof(val2.c_str());

  cutWhite(val1);
  cutBlksXBlks(val1, ",");
  cutParBlks(val1);
  val2 = cutDouble(val1);
  alpha = atof(val2.c_str());
  cutBlksXBlks(val1, ",");
  val2 = cutDouble(val1);
  beta = atof(val2.c_str());
  cutBlksXBlks(val1, ",");
  val2 = cutDouble(val1);
  gamma = atof(val2.c_str());
  return TF;
  }


// ____________________________________________________________________________
// D                  SINGLE PARAMETER OUTPUT FUNCTIONS
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Write In Simple Format, As To Screen
// ----------------------------------------------------------------------------

/* These routines are used for writing out a single parameter into an output
   stream.  The format is set for easy viewing but is NOT compatibile with
   GAMMA parameter format.  That is, this is good for output but not I/O.

        // Input             par  : Parameter (this)
        // 		     ostr : An output stream
        // Output            ostr : The output stream, modified to 
        //                          contain the information about
	//			    the parameter                            */

vector<string> SinglePar::printStrings() const
  {
  vector<string> PStrings;
  if(ParType != -1)
    {
    PStrings.push_back(string("Parameter: ")  + ParName);
    string st;
    switch(ParType)
      {
      case 0:  st = string("Integer (0)");    break;
      case 1:  st = string("Double (1)");     break;
      case 2:  st = string("String (2)");     break;
      case 3:  st = string("Coordinate (3)"); break;
      case 4:  st = string("Tensor (4)");     break;
      default: st = string("Unknown (")
                  + Gdec(ParType) + string(")");
      }
    PStrings.push_back(string("  Type   - ") + st);
    PStrings.push_back(string("  Data   - ") + ParData);
    PStrings.push_back(string("  State  - ") + ParState);
    }
  else
    PStrings.push_back(string("Parameter: Empty"));
  return PStrings;
  }

std::ostream& SinglePar::print(std::ostream& ostr) const
  {
  ostr << "\nParameter: " << ParName 		// Output name (string)
       << "\n  Type   - " << ParType		// Output type (int)
       << "\n  Data   - " << ParData		// Output data (string)
       << "\n  State  - " << ParState		// Output state (string)
       << "\n";
  return ostr;
  }

std::ostream& operator<< (std::ostream& ostr, const SinglePar& par)
  { return par.print(ostr); }


// ----------------------------------------------------------------------------
//             Write In Standardized Format, As To Parameter Set
// ----------------------------------------------------------------------------

/* This writes the parameter into an ASCII output file stream in GAMMA
   Parameter Set format. This allows for the parameter to be re-read into
   a program from the ASCII file using simple commands.              
     
	   Input	    par : A single parameter (this)
	  		  ostr  : An output filestream
	                 namelen: Length of parameter name column (default 10)
	   Return	     TF	: TRUE if parameter sucessfully written
	   Note		   	: The parameter is written in GAMMA
	  			  parameter format
	   Note			: Any new parameter types added to
          			  to GAMMA must be put into the switch       */

bool SinglePar::write(std::ofstream& ostr, int namelen) const
  {
  bool TF = true;				// Assume write will work
  ostr << ParName;				// Write the name 

  int len = namelen - ParName.length();		// Name output length
  if(len <= 0) len=1;				// Insure positve #
  ostr << string(len, ' ');			// Add spaces so columnated
  ostr << " (" << ParType << ") : ";		// Write the type in () & :
  string val1, val2;				// Temp string variables
  switch(ParType)
    {
    case 0:					// Integer parameter
    case 1:					// Double parameter
    case 2:					// string parameter
    case 3:					// Coordinate parameter
    default:					// Default parameter
      ostr << ParData;
      ostr << "\t- " << ParState << "\n";
      break;
    case 4:					// Tensor parameter
      val1 = ParData;
      val2 = cutInt(val1);
      ostr << val2;
      ostr << "\t- " << ParState << "\n";
      cutBlksXBlks(val1, ",");
      val2 = cutParBlks(val1);
      ostr << "\t\t" << val2 << "\n";
      cutBlksXBlks(val1, ",");
      ostr << "\t\t" << val1 << "\n";
      break;
    }
  return TF;
  }

 
// ____________________________________________________________________________
// E                  SINGLE PARAMETER INPUT FUNCTIONS
// ____________________________________________________________________________
     
	// Input	    par : A single parameter (this)
	//		    inp : An input filestream
	// Return	     TF	:  1 if parameter read from fstream
	//			   0 if end of file reached
	//			  -1 if partial parameter read
	// Note		   	: The parameter is filled up
	//			  from the input filestream
	// Note			: Any new parameter types added to
        //			  to GAMMA must be put into the switch
	// Note			: string parameters can NOT contain blanks
	//			  in their data but they may have - in them

int SinglePar::read(std::ifstream& inp)
  {
  string input, scratch;			// Read string, scratch string
  char buf[200];				// Need a character buffer
#ifdef _MSC_VER
  strcpy_s(buf, "");
#else
  strcpy(buf, "");
#endif
  inp.getline(buf, 200, '\n');			// Read line, to \n or 200 chars
  if(inp.eof())					// If we encounter the end of
    {						// the file, we exit if there
    input = string(buf);			// is no data on that line. This
    if(input.length() < 2) return 0;		// avoids missing data on lines
    }						// that contain the file end

  input = string(buf);			// Convert to string
  ParName = cutString(input);			// Remove any initial string
  if(!ParName.length()) return -1;		// Can't read any parameter!
  if(input[0] != '(')   return -1;		// Begin looking for (#)
  input = input.substr(1);			// Clip ( from input start
  scratch = cutInt(input,0);
  if(!scratch.length()) return -1;		// Can't read any parameter!
  ParType = atoi(scratch.c_str());		// Clip integer from input
  if(input[0] != ')')   return -1;		// Look for end of (#)
  input = input.substr(1);			// Clip 0 from input start
  scratch = cutBlksXBlks(input, ":");		// Clip : from string
  if(!scratch.length()) return -1;		// Can't read any parameter!
  int retval = 0;				// Default return value  

  switch(ParType)
    {
    case 0:					// Data is integer
	ParData = cutInt(input);		//	Cut out an int
        cutBlksXBlks(input, "-");		//      Cut remaining past -	
        ParState = input;			//	Rest is the comment
        retval = 1;				//      We succeeded
	break;
      case 1:					// Data is double
	ParData = cutDouble(input);		//	Cut out a double
        cutBlksXBlks(input, "-");		//      Cut remaining past -
        ParState = input;			//	Rest is the comment
        retval = 1;				//      We succeeded
	break;
      case 2:					// Data is string or unknown
      default:
        ParData = cutString(input); 		//	Cut until we find " "
        cutBlksXBlks(input, "-");		//      Cut remaining past -	
        ParState = input;			//	Rest is the comment
        retval = 1;				//      We succeeded
	break;
      case 3:					// Data is coordinate
        retval = setCoord(input);		//	Use private function
	break;
      case 4:					// Data is spatial tensor
        retval = setTensor(inp, input);		//	Use private function
	break;
      }
  return retval;
  }

// ____________________________________________________________________________
// E         SINGLE PARAMETER LIST (Parameter Set) SUPPORT FUNCTIONS
// ____________________________________________________________________________

bool SinglePar::operator==(const SinglePar& par) const
  {
  if(ParName != par.ParName) return false;
  if(ParType != par.ParType) return false;
  if(ParData != par.ParData) return false;
  return true;
  }

bool SinglePar::operator!=(const SinglePar& par) const
  { return (!((*this) == par)); }
 
bool SinglePar::operator<(const SinglePar& par) const
  { return (ParName < par.ParName); }

bool SinglePar::operator>(const SinglePar& par) const
  { return (ParName > par.ParName); }

#endif

