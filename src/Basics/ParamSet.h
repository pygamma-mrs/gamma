/* ParamSet.h ***************************************************-*-c++-*-
**                                                                      **
**                            G A M M A                                 **
**                                                                      **
**      Parameter Set                             	Interface       **
**                                                                      **
**      Scott A. Smith                                                  **
**      Copyright (c) 1997, 2002                                        **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/

/*************************************************************************
**                                                                      **
** Description                                                          **
**                                                                      **
** The class ParameterSet contains a linked list of GAMMA parameters,   **
** the latter embodied in class SinglePar.  Each entry in a Parameter   **
** Set is of type SinglePar, see class SinglePar.                       **
**                                                                      **
*************************************************************************/

#ifndef   ParameterSet_h			// Is this file included?
#  define ParameterSet_h 1			// If not, then remember it
#  if defined(GAMPRAGMA)			// Using the GNU compiler?
#    pragma interface				// Then this is the interface
#  endif

#include <GamGen.h>				// Know OS dependent stuff
#include <Basics/SinglePar.h>                   // Inlcude GAMMA parameters
#include <list>                                 // Include libstdc++ lists
#include <string>				// Include libstdc++ strings
#include <iostream>				// Include libstdc++ IO streams
#include <fstream>				// Include {if,of}streams

/* sosi - gcc 2.8.1 on Solaris 6 has trouble with public std::list<SinglePar>
          as the base class.  This is not true with 2.95.2 on Irix 6.5, 
          MSVC++ 5, and egcs-2.91.60 under CygWin.  Removal of std:: here is
          fine with gcc 2.8.1 though.  Soooo, what to do?  I'll type def it
          and hope for the best....... yes this seems OK with all C++ flavors.
          Here is the story then: I have replaced the simple class declaration
          below with the typedef and the declaration using the typedef.  When
          2.8.2 is ancient I'll just set things back.        June 22 2000

          Now it should be ancient as gcc 3.1 is out. I'll set it back but
          am leaving both forms commented out in case its needed...
                                                             Sept 20 2002   */
                    
/*
class ParameterSet : public std::list<SinglePar>	// Direct class declare 
							//       OR
typedef std::list<SinglePar> stdlistSP;			// Using typedef for
class ParameterSet : public stdlistSP			// indirect declare  */

const std::list<SinglePar> GamSParInit;
const std::vector<int> GamIntVecInit;

typedef std::list<SinglePar> stdlistSP;		// Using typedef on STL list

class ParameterSet : public stdlistSP
{
private:

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i                     CLASS PARAMETER SET ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   pset    : Parameter set (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
				pname   : Added error message for output
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

         void Perror(int eidx, int noret=0) const;
volatile void Pfatality(int eidx) const;
         void Perror(int eidx, const std::string& pname, int noret=0) const;

public:

// ____________________________________________________________________________
// A              PARAMETER SET CONSTRUCTORS/ASSIGNMENTS/DESTRUCTOR
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't 
   have to keep looking them up all the time.

   ParameterSet()				Empty Parameter Set
   ParameterSet(int N)				Parameter Set w/ N Parameters
   ParameterSet(int N, const SinglePar& par)	Parameter Set w/ N pars
   ParameterSet(const ParameterSet& pset)	Parameter Set copy of pset
   ParameterSet assign(N)			Assign N element
*/

// ____________________________________________________________________________
// B                PARAMETER SET ITERATORS & MEMBER ACCESS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't 
   have to keep looking them up all the time.

   list<SinglePar>::iterator ParameterSet::begin()	Pointer to 1st element
   list<SinglePar>::iterator ParameterSet::end()	Pointer to last element
   SinglePar                 ParameterSet::front()	First element
   SinglePar                 ParameterSet::back()	Last element
*/

// ____________________________________________________________________________
// C                PARAMETER SET LIST & QUEUE OPERATIONS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't 
   have to keep looking them up all the time.

   ParameterSet::push_back(const SinglePar& par)	Add par to list end
   ParameterSet::pop_back()				Remove par at list end
   ParameterSet::push_front(const SinglePar& par)	Add par to list start
   ParameterSet::pop_front(const SinglePar& par)	Remove par at list start

   ParameterSet::insert(iterator p, SinglePar& par)	Add par before p
   ParameterSet::erase(iterator p)			Remove par at p
   ParameterSet::clear()				Remove all list entries
*/

// ____________________________________________________________________________
// D                PARAMETER SET ADDITIONAL QUEUE OPERATIONS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't 
   have to keep looking them up all the time.

   int  ParameterSet::size()				Number of entries
   bool ParameterSet::empty()				TRUE if pset empty
   bool ParameterSet::operator==(ParameterSet pset)	TRUE if psets equal
   bool ParameterSet::operator!=(ParameterSet pset)	TRUE if psets not equal
*/

// ____________________________________________________________________________
// E                    PARAMETER SET AUXILIARY FUNCTIONS
// ____________________________________________________________________________
  
/* These functions are "list" type of functions that have been added to make
   the list of single parameters (i.e. parameter list) do simple things needed
   for ready access of single parameters.
  
    Function   Arguments                     Result 
   ----------  ---------  ----------------------------------------------------- 
    contains    string    Returns true/false if parameter with name is in pset
       "       SinglePar  Returns true/false if single parmeters is in pset
     seek       string    Returns iterator in pset for parameter with name
       "       SinglePar  Returns iterator in pset for single parameter
     strip       int      Returns parameter set with parameters named [#]name
   countpar   string,int  Counts contiguous parameters in pset of name name(#)

   Note that in the function "strip" the returned parameter set contains
   parameters whos names no longer are prefixed with [#].                   */

MSVCDLL int contains(const std::string& pname) const;
MSVCDLL int contains(const SinglePar& par)     const;

MSVCDLL stdlistSP::const_iterator seek(const std::string& pname) const;
MSVCDLL stdlistSP::const_iterator seek(const SinglePar&   par)   const;
 
MSVCDLL ParameterSet strip(int indx) const;
MSVCDLL int          countpar(const std::string& pnamein, int idx0=0);

// ____________________________________________________________________________
// F                     PARAMETER SET OUTPUT FUNCTIONS
// ____________________________________________________________________________
 
///F_list <<                     - Standard Output
///F_list write                  - ASCII parameter file

/* These are functions that output formatted information concerning the
   Parameter set to a specified output stream of file. The latter function,
   write, will produce an ASCII file that is self-readable by this class.

                Input           pset : Parameter set (this)
                                ostr : Output ASCII file stream
                                file : ASCII output file (making pset file)
				warn : Warning level flag for write fail
                Return          void : pset is sent to the output stream 
                                int  : T/F write to output pset file.        */

MSVCDLL std::vector<std::string> printStrings()           const;
       std::ostream&     print(std::ostream& out) const;
MSVCDLL friend std::ostream& operator<< (std::ostream& out, const ParameterSet& pset);

MSVCDLL bool write(const std::string& fileout, int warn=2) const;
MSVCDLL bool write(std::ofstream& ofstr,       int warn=2) const;

// ____________________________________________________________________________
// G                    PARAMETER SET INPUT FUNCTIONS
// ____________________________________________________________________________
 
        // Input                pset    : Parameter set (this)
        //                      filein  : Input filename
        //                      inp     : Input filestream
        //                      fflag   : Flag for fatal warnings.
	//			fflag   : Flag for fatal warnings.
	//				     0 : No error messages
	//				     1 : Non-fatal messages
	//				  else : Messages, stop execution
        // Output               int     : True if file filein or input stream
        //                                inp has been scanned for parameters.
        // Note                         : Parameter set filled with parameters
        //                                read from the input file filein
	//				  or input stream inp
        //

MSVCDLL bool read(const std::string& filein, int fflag=0);
MSVCDLL bool read(std::ifstream& inp,        int fflag=0);

// ____________________________________________________________________________
// H                    PARAMETER SET INTERACTIVE FUNCTIONS
// ____________________________________________________________________________

        // Input                pset    : Parameter set (this)
        //                      argc    : Number of arguments
        //                      argv    : Vector of argc arguments
        //                      argn    : Argument index
        // Output               filename: The parameter argn of array argc
        //                                is used to supply a filename
        //                                from which the parameter set is read
        //                                If the argument argn is not in argv,
        //                                the user is asked to supply a filename
        //                                The filename is returned
        // Note                         : The file should be an ASCII file
        //                                containing recognized sys parameters
        // Note                         : The parameter set is modifed (filled)

MSVCDLL std::string ask_read(int argc, char* argv[], int argn);

        // Input                pset    : Parameter set (this)
        //                      value	: Value of parameter desired
        // Output               TF      : True if parameter named "name" exists
        //                                in pset and sucessfully parsed as
        //                                a {string, int, double}. The value
	//				  will be set to parameter value

MSVCDLL bool getParameter(const std::string& name, std::string& value) const;
MSVCDLL bool getParameter(const std::string& name, int&         value) const;
MSVCDLL bool getParameter(const std::string& name, double&      value) const;

MSVCDLL bool getString(const std::string& name, std::string& value) const;
MSVCDLL bool getInt(const    std::string& name, int&         value) const;
MSVCDLL bool getDouble(const std::string& name, double&      value) const;

};

//   typedef std::multimap<std::string, std::string> Map;
//   typedef Map::iterator iterator;
   
#endif						// End ParamSet.h
