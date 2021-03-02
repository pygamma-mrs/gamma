/* GamTest.h ****************************************************-*-c++-*-
**									**
** 	                        G A M M A				**
**									**
**	GAMMA Test                                 Interface		**
**								 	**
**      Scott Smith                                                     **
**      Copyright (c) 2000                                              **
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**      email: gamma@magnet.fsu.edu                                     **
**      www: http://gamma.magnet.fsu.edu                                **
**                                                                      **
**      $Header: $
**                                                                      **
*************************************************************************/
     
/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
** Class GamTest embodies a test set for the GAMMA platform testing     **
** hierarchy. Such tests are specific to a particular set of GAMMA      **
** modules. A GAMMA test will always deal with set of modules, or       **
** equivalently, a set of source subdirectories from the main GAMMA     **
** source directory, $GAMMA/src/. Thus, this class allows for facile    **
** construction of bulk GAMMA testing. It does so by running module     **
** tests where each module test is associated with a subdirectory from  **
** the main GAMMA source directory.					**
**                                                                      **
*************************************************************************/

#ifndef   GamTest_h_			// Is this file already included?
#  define GamTest_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Testing/ModTest.h>		// Include class tests
#include <list>				// Include libstdc++ STL lists
#include <iostream>			// Include file streams
#include <fstream>			// Include {if,of}streams
#include <string>			// Include string class

/* sosi - gcc 2.8.1 on Solaris 6 has trouble with public std::list<ModTest>
          as the base class.  This is not true with 2.95.2 on Irix 6.5, 
          MSVC++ 5, and egcs-2.91.60 under CygWin.  Removal of std:: here is
          fine with gcc 2.8.1 though.  Soooo, what to do?  I'll type def it
          and hope for the best....... yes this seems OK with all C++ flavors.
          Here is the story then: I have replaced the simple class declaration
          below with the typedef and the declaration using the typedef.  When
          2.8.2 is ancient I'll just set things back.        June 22 2000    */

//class GamTest: public std::list<ModTest>

typedef std::list<ModTest> stdlistGT;

class GamTest: public stdlistGT
{
       int         TestStatus;               // Test status (True/False)
  std::string      TestName;                 // Test name
  std::string      TestDescription;          // Test description
       int         TestRunLevel;             // Test run level
  std::vector<int> TestResults;              // Array of test results

private:
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                      GAMMA TEST ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   GT      : A GAMMA test (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

         void GTerror(int eidx, int noret=0) const;
volatile void GTfatal(int eidx)              const;

// ____________________________________________________________________________
// ii                   GAMMA TEST INITIALIZATION FUNCTIONS
// ____________________________________________________________________________

/* This function will initialize the test results vector TestResults. This      
   should be done whenever testing is done if either 1.) The status flag is
   -1, indicating no tests have ever been run, or if the current number of
   tests does not match the size of the tests results vector.  The latter
   indicates that additional tests have been added, potentially with a
   different ordering, since the test results were generated previously.    */

void SetResults(int force=0);                 

// ____________________________________________________________________________
// iii                 GAMMA TEST MODULE TEST INDEXING FUNCTIONS
// ____________________________________________________________________________
 
/* These function deal with indicies of Module Test tests. They allow one to
   obtain the test Pix from either an interger or a string and they perform
   range checking to insure that the requested test exists in the section.   */
 
bool                          CheckIndex(int k, int w=1) const;
std::list<ModTest>::iterator       GetPixNC(int k);
std::list<ModTest>::const_iterator GetPix(int k)              const;
std::list<ModTest>::const_iterator GetPix(const std::string& N)    const;
int                           GetIndex(const std::string& N)  const;

public:

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A            GAMMA TEST CONSTRUCTION, DESTRUCTION, ASSIGNMENT
// ____________________________________________________________________________
 
/* Were it not for the additional variables defined in this class, all the
   constructors would just be inherited from the ANSI standard "list" class
   in the C++ library libstdc++. Since this is not the case, I am listing all
   such constructors (commented out) and will only (re)define the ones that
   I care about.

GamTest()                               Empty GAMMA Test
GamTest(int N)                          Module Test w/ N Tests
GamTest(int N, const ModTest& MT)       Module Test w/ N Module Tests
GamTest(const GamTest& GT)              GAMMA Test copy of GT
GamTest assign(N)                       Assign N element
~GamTest()                              Destructor of Module Test           */

MSVCDLC            GamTest();                   // Empty Module Test
MSVCDLC            GamTest(const  GamTest& GT); // Module Test copy of GT
MSVCDLL GamTest& operator= (const GamTest& GT); // Assignment operator
 
// ____________________________________________________________________________
// B                   GAMMA TEST ITERATORS & MEMBER ACCESS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't
   have to keep looking them up all the time.

list<SectTest>::iterator begin()      Pointer to 1st element
list<SectTest>::iterator end()        Pointer to last element
SectTest                 front()      First element
SectTest                 back()       Last element                  */

// ____________________________________________________________________________
// C                     GAMMA TEST LIST & QUEUE OPERATIONS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't 
   have to keep looking them up all the time.

push_back(const  SectTest& MT)         Add MT to list end
pop_back()                             Remove MT at list end
push_front(const SectTest& MT)         Add MT to list start
pop_front(const  SectTest& MT)         Remove MT at list start
insert(iterator p, SectTest& MT)       Add MT before p
erase(iterator p)                      Remove MT at p
clear()                                Remove all list entries      */

// ____________________________________________________________________________
// D                    GAMMA TEST ADDITIONAL QUEUE OPERATIONS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't 
   have to keep looking them up all the time.

int  size()                          Number of entries
bool empty()                         TRUE if pset empty
bool operator==(GamTest GT)          TRUE if psets equal
bool operator!=(GamTest GT)          TRUE if psets not equal        */

// ____________________________________________________________________________
// E                     GAMMA TEST LIST AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 
/* These functions are "list" type of functions that have been added to make 
   the list of Module Tests (i.e. GAMMA test) do simple things needed for ready 
   access of Moduel Tests.
   
    Function   Arguments                     Result 
   ----------  ---------  ----------------------------------------------------- 
    contains    string    Returns true/false if test of name is in GAMMA tests
       "       ModTest    Returns true/false if module test is in GAMMA tests
     seek       string    Returns iterator in GAMMA tests for test with name
       "       ModTest    Returns iterator in GAMMA tests for module test    */ 
 
MSVCDLL int contains(const std::string& N) const;
MSVCDLL int contains(const ModTest& MT)  const;

MSVCDLL std::list<ModTest>::const_iterator seek(const std::string& N) const;
MSVCDLL std::list<ModTest>::const_iterator seek(const ModTest& MT) const;
 
// ____________________________________________________________________________
// F                       GAMMA TEST ACCESS FUNCTIONS
// ____________________________________________________________________________
 
/* These functions allow access to simple details of individual Module Tests
   or all tests in one lump sum. They will perform NO TESTING (we want code
   a bit more sophisticated for that - see a later section).  These functions
   facilitate external programs in generating nicely formatted output with the
   GAMMA Test and Module Test results. The lower case names parallel those
   that are present in the Class, Section and Single Test classes.

     Function                             Result
  ---------------  ---------------------------------------------------------
      GetName      Returns string for name of GAMMA or Module Test
     GetNames      Returns string vector names of all Module Tests
  GetDescription   Returns string for description of GAMMA or Module Test
  GetDescriptions  Returns string vector descriptions of all Module Tests
     GetStatus     Returns int for status of or Module Test
    GetStatuses    Returns int vector for status of all Module Tests
    GetRunLevel    Returns int for run level of GAMMA or Module Test
   GetRunLevels    Returns int vector for run level of all Module Tests
    GetResults     Returns int vector for all current Module Test Results     */

MSVCDLL std::string    GetName()                    const;
MSVCDLL std::string    GetName(int k)               const;
//std::string    GetName(int k, int j)        const;
MSVCDLL std::vector<std::string> GetNames()                   const;
//vector<string> GetNames(int k)              const;
MSVCDLL std::string         GetDescription()             const;
MSVCDLL std::string         GetDescription(int k)        const;
//string         GetDescription(int k, int j) const;
MSVCDLL std::vector<std::string> GetDescriptions()            const;
//vector<string> GetDescriptions(int k)       const;
MSVCDLL int            GetStatus()                  const;
MSVCDLL int            GetStatus(int k)             const;
MSVCDLL std::vector<int>    GetStatuses()                const;
MSVCDLL int            GetRunLevel()                const;
MSVCDLL int            GetRunLevel(int k)           const;
MSVCDLL std::vector<int>    GetRunLevels()               const;
MSVCDLL std::vector<int>    GetResults();
MSVCDLL std::vector<int>    GetResults(int k);

/*         Input                GT : A GAMMA test (this)
                  (name)      name : The name of the test     (setting name)
                  (status)  status : Current test status      (setting status)
                  (descr)    descr : Current test description (setting describ)
                  (runlev)  runlev : Current test run level   (setting runlev)
           Output (name)      name : The name of the test     (getting name)
                  (status)  status : Current test status      (getting status)
                  (descr)    descr : Current test description (getting describ)
                  (runlev)  runlev : Current test run level   (getting runlev)
                              void : If setting interal value                */  
 
MSVCDLL const std::string& name()     const;                // Get name
MSVCDLL       int          status()   const;                // Get status
MSVCDLL const std::string& describe() const;                // Get descript
MSVCDLL       int          runlevel() const;                // Get runlev
MSVCDLL       void         name(const  std::string& Name);  // Set name
MSVCDLL       void         status(     int          Status);// Set status
MSVCDLL       void         describe(const std::string& D);  // Set decript
MSVCDLL       void         runlevel(   int          RLev);  // Set runlev
 
// ____________________________________________________________________________
// G                      GAMMA TEST TESTING FUNCTIONS
// ____________________________________________________________________________
 
/* These functions will invoke any or all of the Module Tests that are listed
   as part of this GAMMA Test. The results of the Module Tests will be stored
   in the vector TestResults and the summed results placed into TestStatus.
   The value of TestStatus (default -1) indicates the whether the tests have
   been run (as does the vector TestResults).  There are three basic means for
   testing:
             1.) Run All Of The Module Tests                   TestMods
             2.) Run All Module Tests Until A Failure Occurs   TestMod
             3.) Run A Test On A Specific Module Test          TestMod

   In each case the "force" flag dictates whether these tests should be run
   anew or alternatively reuse test results if previsouly run (status!=-1).  */
 
//-----------------------------------------------------------------------------
//                             Test All Modules
//-----------------------------------------------------------------------------

MSVCDLL int TestMods(std::ostream& ostr, int anew=0, int keepon=0);
 
//-----------------------------------------------------------------------------
//                           Test Specific Module
//-----------------------------------------------------------------------------
 
MSVCDLL int TestMod(std::ostream& ostr, const int tidx,     
                                                 int anew = 0, int keepon = 0);
MSVCDLL int TestMod(std::ostream& ostr, const std::string& tnam,
                                                 int anew = 0, int keepon = 0);

// ____________________________________________________________________________
// H                       GAMMA TEST OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions facilitate output of GAMMA tests in a nicely formatted
   fashion. Note that when test results are listed only the status flags are
   utilized. That is, NO TESTS ARE RUN.  So, any testing should be performed
   prior to calling test output functions.

   Header:  This function makes a nice title for the GAMMA test, below
   ======   which the PASS/FAIL results of the module tests should be listed.
            Module test PASS/FAIL listing is done in the Module Test
            class (ModTest::Result) and should conform to the format set
            herein.
   Results: This function loops through and spits out all the Module Test
   ResRec:  Whereas the above Results function tabulates module test results,
   ======   this function will go into each module and output results from
            individual class tests. Typically this is done only on modules
            that have failed their test. It will also have class tests call
            results from individual section tests if the function resursion
            level (nlevels) is set greater than 1.                           */

MSVCDLL std::ostream& Header(std::ostream& ostr)              const;
MSVCDLL std::ostream& Results(std::ostream& ostr, int goon=1) const;
MSVCDLL std::ostream& ResRec(std::ostream& ostr,  int goon=1, int nlevels=1);


// ____________________________________________________________________________
// I                       GAMMA TEST INTERACTIVE FUNCTIONS
// ____________________________________________________________________________
 
/* These functions facilitate the running of a GAMMA test. Recall that each 
   test runs through a series of module tests, ultimately the test returning
   TRUE/FALSE on completion. However there are several additional facets of
   GAMMA testing.
 
   1.) Reuse Last Module Test Results: If previously performed, we can just
                                       bypass a test and use its T/F status
   2.) Continue Through Failed Tests:  If a module test fails should testing
                                       stop or should one go on to the next
   3.) Recurse Through Class Tests:    If a module test fails should testing
                                       recurse into class tests for that
                                       module , displaying specific failures?
 
   These functions assume you want some output text from running the tests,
   otherwise you would just use the TestLevel function directly, right?      */  

MSVCDLL int AskRun(std::ostream& ostr);


// ____________________________________________________________________________
// I                       SECTION TEST AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 
MSVCDLL std::ostream& ListTests(std::ostream& ostr)  const;
MSVCDLL std::ostream& FinishTest(std::ostream& ostr) const;


};

#endif							// GamTest.h
