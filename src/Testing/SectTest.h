/* SectTest.h ***************************************************-*-c++-*-
**									**
** 	                        G A M M A				**
**									**
**	Section Test                                 Interface		**
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
** Class SectTest embodies a test set for the GAMMA platform testing    **
** hierarchy. Such tests are specific to a particular GAMMA section in  **
** a source code file. A section test will always deal with one section **
** of a source file (class or collection of associated functions) that  **
** lies within a single GAMMA Module. A section test may contain 	**
** multiple single tests.                                               **
**                                                                      **
*************************************************************************/

#ifndef   SectTest_h_			// Is this file already included?
#  define SectTest_h_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <Testing/SingleTest.h>		// Include single tests
#include <list>				// Include libstdc++ STL lists
#include <vector>			// Include libstdc++ STL vectors
#include <iostream>			// Include file streams
#include <fstream>			// Include {if,of}streams
#include <string>			// Include string class

/* sosi - gcc 2.8.1 on Solaris 6 has trouble with public std::list<SingleTest>
          as the base class.  This is not true with 2.95.2 on Irix 6.5, 
          MSVC++ 5, and egcs-2.91.60 under CygWin.  Removal of std:: here is
          fine with gcc 2.8.1 though.  Soooo, what to do?  I'll type def it
          and hope for the best....... yes this seems OK with all C++ flavors.
          Here is the story then: I have replaced the simple class declaration
          below with the typedef and the declaration using the typedef.  When
          2.8.2 is ancient I'll just set things back.        June 22 2000    */

//class SectTest: public std::list<SingleTest>

typedef std::list<SingleTest> stdlistST;

class SectTest: public stdlistST
{
       int         TestStatus;               // Test status (True/False)
  std::string      TestName;                 // Test name
  std::string      TestDescription;          // Test description
       int         TestRunLevel;	     // Test run level
  std::vector<int> TestResults;		     // Array of test results

private:
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                      SECTION TEST ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   ST      : A single test (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

         void STerror(int eidx, int noret=0) const;
volatile void STfatality(int eidx)           const;

// ____________________________________________________________________________
// ii                 SECTION TEST INITIALIZATION FUNCTIONS
// ____________________________________________________________________________
 
/* This function will initialize the test results vector TestResults. This 
   should be done whenever testing is done if either 1.) The status flag is
   -1, indicating no tests have ever been run, or if the current number of
   tests does not match the size of the tests results vector.  The latter 
   indicates that additional tests have been added, potentially with a
   different ordering, since the test results were generated previously.    */
 
void SetResults(int force=0);

// ____________________________________________________________________________
// iii               SECTION TEST SINGLE TEST INDEXING FUNCTIONS
// ____________________________________________________________________________

/* These function deal with indicies of Single Test tests. They allow one to
   obtain the test Pix from either an interger or a string and they perform
   range checking to insure that the requested test exists in the section.   */

bool                                  CheckIndex(int k, int w=1) const;
std::list<SingleTest>::iterator       GetPixNC(int k);
std::list<SingleTest>::const_iterator GetPix(int k)              const;
std::list<SingleTest>::const_iterator GetPix(const   std::string& N)  const;
int                                   GetIndex(const std::string& N)  const;

public:

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A            SECTION TEST CONSTRUCTION, DESTRUCTION, ASSIGNMENT
// ____________________________________________________________________________
 
/* Were it not for the additional variables defined in this class, all the
   constructors would just be inherited from the ANSI standard "list" class 
   in the C++ library libstdc++. Since this is not the case, I am listing all
   such constructors (commented out) and will only (re)define the ones that
   I care about.

SectTest()                               Empty Section Test
SectTest(int N)                          Section Test w/ N Tests
SectTest(int N, const SingleTest& PT)    Section Test w/ N PT
SectTest(const SectTest& PT)             Section Test copy of PT
SectTest assign(N)                       Assign N element
~SectTest()                              Destructor of Section Test          */

MSVCDLC            SectTest();                    // Empty Section Test
MSVCDLC            SectTest(const SectTest& PT);  // Section Test copy of PT
MSVCDLL SectTest& operator= (const SectTest& PT); // Assignment operator

// ____________________________________________________________________________
// B                   SECTION TEST ITERATORS & MEMBER ACCESS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't
   have to keep looking them up all the time.

list<SingleTest>::iterator begin()      Pointer to 1st element
list<SingleTest>::iterator end()        Pointer to last element
SingleTest                 front()      First element
SingleTest                 back()       Last element               */

// ____________________________________________________________________________
// C                     SECTION TEST LIST & QUEUE OPERATIONS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't 
   have to keep looking them up all the time.

push_back(const  SingleTest& ST)       Add ST to list end
pop_back()                             Remove ST at list end
push_front(const SingleTest& ST)       Add ST to list start
pop_front(const  SingleTest& ST)       Remove ST at list start

insert(iterator p, SingleTest& ST)     Add ST before p
erase(iterator p)                      Remove ST at p
clear()                                Remove all list entries     */

// ____________________________________________________________________________
// D                    SECTION TEST ADDITIONAL QUEUE OPERATIONS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't 
   have to keep looking them up all the time.

int  size()                            Number of entries
bool empty()                           TRUE if pset empty
bool operator==(SectTest PT)           TRUE if psets equal
bool operator!=(SectTest PT)           TRUE if psets not equal     */

// ____________________________________________________________________________
// E                     SECTION TEST LIST AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 
/* These functions are "list" type of functions that have been added to make 
   the list of single test (i.e. class test) do simple things needed for ready 
   access of single test. 
   
    Function   Arguments                     Result 
   ----------  ---------  ----------------------------------------------------- 
    contains    string    Returns true/false if test of name is in class tests
       "      SingleTest  Returns true/false if single test is in class tests
     seek       string    Returns iterator in class tests for test with name
       "      SingleTest  Returns iterator in class tests for single test    */ 
 
MSVCDLL int contains(const std::string& tname) const;
MSVCDLL int contains(const SingleTest& ST)     const;

MSVCDLL std::list<SingleTest>::const_iterator seek(const std::string& tname) const;
MSVCDLL std::list<SingleTest>::const_iterator seek(const SingleTest& ST) const;

// ____________________________________________________________________________
// F                       SECTION TEST ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions allow access to simple details of individual Single Tests
   or all tests in one lump sum. They will perform NO TESTING (we want code
   a bit more sophisticated for that - see a later section).  These functions
   facilitate external programs in generating nicely formatted output with the
   Section Test and Single Test results. The lower case names parallel those
   that are present in the Single Test class.
 
     Function                             Result
  ---------------  ---------------------------------------------------------
      GetName      Returns string for name of Section or Single Test
     GetNames      Returns string vector names of all Single Tests
  GetDescription   Returns string for description of Section or Single Test
  GetDescriptions  Returns string vector descriptions of all Single Tests
     GetStatus     Returns int for status of Section or Single Test
    GetStatuses    Returns int vector for status of all Single Tests
    GetRunLevel    Returns int for run level of Section or Single Test
   GetRunLevels    Returns int vector for run level of all Single Tests 
    GetResults     Returns int vector for all current Single Test Results   */


MSVCDLL std::string              GetName()             const;
MSVCDLL std::string              GetName(int k)        const;
MSVCDLL std::vector<std::string> GetNames()            const;
MSVCDLL std::string              GetDescription()      const;
MSVCDLL std::string              GetDescription(int k) const;
MSVCDLL std::vector<std::string> GetDescriptions()     const;
MSVCDLL int                      GetStatus()           const;
MSVCDLL int                      GetStatus(int k)      const;
MSVCDLL std::vector<int>         GetStatuses()         const;
MSVCDLL int                      GetRunLevel()         const;
MSVCDLL int                      GetRunLevel(int k)    const;
MSVCDLL std::vector<int>         GetRunLevels()        const;
MSVCDLL std::vector<int>         GetResults();


/*         Input                ST : A section test (this)
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
// G                      SECTION TEST TESTING FUNCTIONS
// ____________________________________________________________________________

/* These functions will invoke any or all of the Single Tests that are listed
   as part of this Section Test. The results of the Single Tests will be stored
   in the vector TestResults and the summed results placed into TestStatus.
   The value of TestStatus (default -1) indicates the whether the tests have
   been run (as does the vector TestResults).  There are three basic means for
   doing the testing:
                       1.) Run All Of The Single Tests 
                       2.) Run All Single Tests Until A Failure Occurs
                       3.) Run A Test On A Specific Single Test
   
   In each case the "force" flag dictates whether these tests should be run 
   anew or alternatively reuse test results if previsouly run (status!=-1).  */

//-----------------------------------------------------------------------------
//                             Run All Single Tests
//-----------------------------------------------------------------------------

MSVCDLL int TestSingles(std::ostream& ostr, int anew=0, int keepon=0);

//-----------------------------------------------------------------------------
//                  Run Specific Single Test (By Index Or Name)
//-----------------------------------------------------------------------------
 
MSVCDLL int TestSingle(std::ostream& ostr, int tidx, int anew=0);
MSVCDLL int TestSingle(std::ostream& ostr, const std::string& tnam, int a=0);

// ____________________________________________________________________________
// H                       SECTION TEST OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions facilitate output of section testing in a nicely formatted
   fashion. Note that when test results are listed only the status flags are
   utilized. That is, NO TESTS ARE RUN.  So, any testing should be performed
   prior to calling test output functions.

   Header:  This function makes a nice title for the section's single tests,
   ======   below which the PASS/FAIL results of the single tests should be
            listed. Single test PASS/FAIL listing is done in the SingleTest
            class (SingleTest::Result) and should conform to the format set
            herein.
   Result:  This function makes a simple single line output for listing under
   ======   ClassTest headers (ClassTest::Header).  The format of the output
            should conform to the format of that header.
   Results: This function loops through and spits out all the SingleTest
   =======  class results of the section (SingleTest::Result)
   ResRec:  Whereas the above Results function tabulates single test results,
   ======   this function will go into each single test output its detailed
            results. Typically this is done only on single tests that
            have failed.                                                    */

MSVCDLL std::ostream& Header(std::ostream&  ostr, const std::string& CN);
MSVCDLL std::ostream& Header(std::ostream&  ostr)             const;
MSVCDLL std::ostream& Result(std::ostream&  ostr)             const;
MSVCDLL std::ostream& Results(std::ostream& ostr, int goon=1) const;
MSVCDLL std::ostream& ResRec(std::ostream&  ostr, int keepon, int nl=1);

// ____________________________________________________________________________
// I                   SECTION TEST INTERACTIVE FUNCTIONS
// ____________________________________________________________________________

/* These functions facilitate the running a section test. Recall that each
   section test runs through a series of single tests, ultimately the test 
   returning TRUE/FALSE on completion. However there are several additional 
   facets of section testing.

   1.) Continue Through Failed Tests:   If a single test fails should testing
                                        stop or should one go on to the next
   2.) Recurse Through Single Tests:    If a single test fails should testing
                                        recurse into that test while displaying
                                        specific test outcome features?

      Recursion Levels:   0 = Output Single Test Pass/Fail Results 
                          1 = Output Single Test Full Details 

   These functions assume you want some output text from running the tests,
   otherwise you would just use the TestLevel function directly, right?     */

MSVCDLL int AskRun(std::ostream& ostr);

// ____________________________________________________________________________
// J                       SECTION TEST AUXILIARY FUNCTIONS
// ____________________________________________________________________________

MSVCDLL std::ostream& ListTests(std::ostream& ostr) const;
};

#endif							// SectTest.h

