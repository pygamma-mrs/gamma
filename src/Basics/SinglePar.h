/* SinglePar.h **************************************************-*-c++-*-
**									**
** 	                        G A M M A				**
**									**
**	Single Parameter                               Interface	**
**								 	**
**      Scott Smith                                                     **
**      Copyright (c) 1996                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**      $Header: $
**								 	**
*************************************************************************/

/*************************************************************************
**								 	**
**  Description							 	**
**								 	**
** Class SinglePar stores a single parameter used in storing parameter  **
** sets and the ASCII I/O of parameters (&classes) to/from files.	**
** Access functions are provided to get at the individual entities in	**
** the parameter.							**
**									**
*************************************************************************/
     
#ifndef   SinglePar_h_			// Is this file already included?
#define   SinglePar_h_  1		// If no, then remember it

#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know OS dependent stuff
#include <iostream>			// Include file streams
#include <fstream>			// Include {if,of}streams
#include <string>			// Include string class
#include <vector>			// Include libstdc++ STL vectors

class SinglePar
  {
  std::string ParName;			// Parameter name
  int         ParType;			// Parameter type (int, double, char,.)
  std::string ParData;			// Parameter data (as a string)
  std::string ParState;			// Parameter statement

public:
 

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                    SINGLE PARAMETER ERROR HANDLING
// ____________________________________________________________________________

        // Input                par     : A parameter (this)
        //                      eidx    : Flag for error type
        //                      noret   : Flag for return (0=return)
        // Output               none    : Error message
        //                                Program execution stopped

         void SPerror(int   eidx, int noret=0) const;
volatile void SPfatality(int eidx)             const;

// ____________________________________________________________________________
// ii                       TYPE MORPHING FUNCTIONS
// ____________________________________________________________________________

        // Input              par : A single parameter (this)
	//                   Type : The parameter type
	//                   Name : The parameter name
	//                   State: The parameter statement
        // Output            void : The parameter name is changed
	// Note			  : These are a little finicky as they
	//			    enforce some parameter types to have
	//			    reasonable data values
 
// Additional functions name & state are public and found in file section B

void type(int         Type);
void data(std::string Data);

// ____________________________________________________________________________
// iii                SINGLE PARAMETER INPUT FUNCTIONS
// ____________________________________________________________________________

/* These Are Functions That Facilitate Input Of Parameters From ASCII Files.
   They Could Be Abused So They Are Kept Private Here.                       */

int setCoord(std::string& input);

	// Input	    par : A single parameter (this)
	//		    line: String representing part of ASCII input line
	// Return	     TF	: TRUE if coord parameter set from line
 

int setTensor(std::ifstream& inp, std::string& input);

	// Input	    par : A single parameter (this)
	//		    inp : Input parameter stream
	//		    line: String representing part of ASCII input line
	// Return	     TF	: TRUE if tensor parameter set from line
 

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A              SINGLE PARAMETER CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________

// ----------------------------------------------------------------------------
//                  Simple Constructors That Don't Do Much
// ----------------------------------------------------------------------------

MSVCDLC SinglePar();				// Default constructor
MSVCDLC SinglePar(const SinglePar& par);	// Self constructor

// ----------------------------------------------------------------------------
//                Constructors That Set Up The Parameter
// ----------------------------------------------------------------------------
 
        // Input            pname : A string for the parameter name
        //                  pdata : A {int, double, string} for parameter data
        //                  ptype : An integer designating the parameter type
        //                 pstate : A string for a comment concerning parameter

MSVCDLC SinglePar(const std::string& pname, int ptype, const std::string& pdata, const std::string& pstate);
MSVCDLC SinglePar(const std::string& pname,                           int pdata, const std::string& pstate);
MSVCDLC SinglePar(const std::string& pname,                        double pdata, const std::string& pstate);
MSVCDLC SinglePar(const std::string& pname,            const std::string& pdata, const std::string& pstate);

// ----------------------------------------------------------------------------
//                Odd Constructors, Destructor, Assignment
// ----------------------------------------------------------------------------

        // Input            pname : A string for a parameter name
        // Output            this : Constructs a dummy SinglePar with the
        //                          name input
 	// Note			  : This ONLY sets the name, nothing else
  
                   SinglePar(const std::string& pname);
MSVCDLC            ~SinglePar();
MSVCDLL SinglePar& operator=(const SinglePar& par);

// ____________________________________________________________________________
// B                  SINGLE PARAMETER ACCESS FUNCTIONS
// ____________________________________________________________________________
 
        // Input               par : A single parameter (this)
        // Output (name)      name : The name of parameter par
	//        (type)      type : Integer designation of type
        //        (data)      data : string containing parameter data
        //        (state)    state : string containing parameter comment

  
MSVCDLL const std::string&  name()  const;		// Can't alter name here
MSVCDLL       int           type()  const;		// Can't alter type here
MSVCDLL const std::string&  data()  const;		// Can't alter data here
MSVCDLL const std::string&  state() const;		// Can't alter state here


        // Input              par : A single parameter (this)
	//                   Name : The parameter name
	//                   State: The parameter statement
        // Output            void : The parameter name is changed

// Additional functions type & data are private and found in file section ii

MSVCDLL void name(const  std::string& Name);		// This sets name to Name
MSVCDLL void state(const std::string& State);	// This sets state to State

// ____________________________________________________________________________
// C                      SINGLE PARAMETER PARSING FUNCTIONS
// ____________________________________________________________________________
 
/* These functions take a generic parameter and attempt to split it apart
   into components as indicated in the function argument list.  They return
   TF depending upon whether the parsing was sucessful.                      */
 
        // Input              par : A parameter (this)
        //                   name : string for parameter name
        //                  state : Parameter statement
        // Output            none : Name and statement set from par
 
MSVCDLL void GetNS(std::string& name, std::string& state) const;

        // Input              par : A parameter (this)
        //                   name : string for parameter name
        //                    val : Parameter value
        //               dx,dy,dz : Double parameter values
        //                   rank : An integer value
        //         diso,delz,deta : Double parameter values
        //       alpha,beta,gamma : Double parameter values
        //                   flag : Warning flag if type cast needed
        //                  state : Parameter statement
        // Output              TF : Returns TRUE if the parameter is
        //                          properly parsed
        // int,double,string	  : name, val and state are set
	// dx, dy, dz		  : name, dx, dy, dz, and state are set
	// rank,diso,delz,deta..  : name,rank,diso,delz,deta,alpha,beta,.. set
 

MSVCDLL bool parse(std::string& name, int&         val, std::string& state, int warn=0) const;
MSVCDLL bool parse(std::string& name, double&      val, std::string& state, int warn=0) const; 
MSVCDLL bool parse(std::string& name, std::string& val, std::string& state, int warn=0) const;

MSVCDLL bool parse(std::string& name, double& dx, double& dy, double& dz,
                                               std::string& state, int warn=0) const;

MSVCDLL bool parse(std::string& name, int& rank, double& diso, double& delz,
               double& deta, double& alpha, double& beta, double& gamma,
                                               std::string& state, int warn=0) const;

// ____________________________________________________________________________
// D                  SINGLE PARAMETER OUTPUT FUNCTIONS
// ----------------------------------------------------------------------------
 
// ----------------------------------------------------------------------------
//                   Write In Simple Format, As To Screen
// ----------------------------------------------------------------------------

        // Input             par  : Parameter (this)
        //                   ostr : An output stream
        // Output            ostr : The output stream, modified to 
        //                          contain the information about
        //                          the parameter

MSVCDLL        std::vector<std::string> printStrings() const;
MSVCDLL        std::ostream&            print(std::ostream& ostr) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const SinglePar& par);

// ----------------------------------------------------------------------------
//             Write In Stanardiszed Format, As To Parameter Set
// ----------------------------------------------------------------------------

MSVCDLL bool write(std::ofstream& ostr, int namelen=10) const;

        // Input            par : A single parameter (this)
        //                ostr  : An output filestream
        //               namelen: Length of parameter name column
        // Return            TF : TRUE if parameter sucessfully written
        // Note                 : The parameter is written in GAMMA
        //                        parameter format
        // Note                 : Any new parameter types added to
        //                        to GAMMA must be put into the switch


// ____________________________________________________________________________
// E                  SINGLE PARAMETER INPUT FUNCTIONS
// ____________________________________________________________________________

MSVCDLL int read(std::ifstream& inp);
     
        // Input            par : A single parameter (this)
        //                  inp : An input filestream
        // Return           T/F : Returns true if parameter read
        //                        from the input filestream
        // Note                 : Any new parameter types added to
        //                        to GAMMA must be put into this function

// ____________________________________________________________________________
// E         SINGLE PARAMETER LIST (Parameter Set) SUPPORT FUNCTIONS
// ____________________________________________________________________________

MSVCDLL bool operator==(const SinglePar& par) const;
MSVCDLL bool operator!=(const SinglePar& par) const;
MSVCDLL bool operator< (const SinglePar& par) const;
MSVCDLL bool operator> (const SinglePar& par) const;

};

/*****************************************************************************/
/*****************************************************************************/
/*                  CLASS SINGLE PARAMETER INTRINSIC TYPES                   */
/*****************************************************************************/
/*****************************************************************************/

/* All current GAMMA parameters of type SinglePar strore their data in
   string format.  The data string must be then parsed to produce values
   in accordance with the parameter.  Below details the parameters which
   are defined by GAMMA via type and how the data string exists.  Any
   class that uses variables of type SinglePar must parse via this format.

        Type    External    string Format	  Example Parse Code
      (ParType)    Type         (ParData)                   
      ------- ---------- ------------------- ----------------------------
         0    Integer           val	     int i = atoi(par.data())
         1    Double            val	     double d = atof(par.data())
         2    string	        val          string s = par.data()

         3    Coordinate     (x,y,z)         string dat = par.data();
                                             double x=0,y=0,z=0;
  		                             string Sval;
                                             cut(dat,"([ \t]*");
  				             Sval = cut(dat,RXdouble);
   				             x = atof(Sval);
   				             cut(dat,"[ \t]*,[ \t]*");
   				             Sval = cut(dat,RXdouble);
          				     y = atof(Sval);
             				     cut(dat,"[ \t]*,[ \t]*");
          				     Sval = cut(dat,RXdouble);
         				     z = atof(Sval);
          				     coord pt(x,y,z);

         4    Tensor	 i, (x,y,z), (a,b,c) string dat = par.data();
                                             int rank;
                                             double iso, delz, eta;
                                             double alpha, beta, gamma;
  		                             string Sval;
                                             Sval = cut(dat,RXint);
                                   	     rank = atoi(Sval)
                                             cut(dat,"([ \t]*");
  				             Sval = cut(dat,RXdouble);
   				             iso = atof(Sval);
   				             cut(dat,"[ \t]*,[ \t]*");
   				             Sval = cut(dat,RXdouble);
          				     delz = atof(Sval);
             				     cut(dat,"[ \t]*,[ \t]*");
          				     Sval = cut(dat,RXdouble);
         				     eta = atof(Sval);
   				             cut(dat,"[ \t]*,[ \t]*");
                                             cut(dat,"([ \t]*");
  				             Sval = cut(dat,RXdouble);
   				             alpha = atof(Sval);
   				             cut(dat,"[ \t]*,[ \t]*");
   				             Sval = cut(dat,RXdouble);
          				     beta= atof(Sval);
             				     cut(dat,"[ \t]*,[ \t]*");
          				     Sval = cut(dat,RXdouble);
         				     gamma = atof(Sval);            */


#endif								// SinglePar.h

