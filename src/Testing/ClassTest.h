/* ClassTest.h *************************************************-*-c++-*-
**									**
** 	                        G A M M A				**
**									**
**	Class Test                                    Interface		**
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
** Class ClassTest embodies a test set for the GAMMA platform testing   **
** hierarchy. Such tests are specific to a particular GAMMA class or    **
** source code file. A class test will always deal with one source file **
** (class or collection of associated functions) within one GAMMA       **
** Module. That is, any one class (or file of functions) will be named  **
** $GAMMA/src/"module"/"class". Thus this class will be intrinsically   **
** associated with such a file and all its tests apply to it. This is   **
** done through a series of tests associated with each section of the   **
** class (or file of functions).                                        **
**                                                                      **
*************************************************************************/

#ifndef   ClassTest_h_			// Is this file already included?
#  define ClassTest_h_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma interface			// Then this is the interface
#  endif

#include <GamGen.h>			// Know MSVCDLL (__declspec)
#include <Testing/SectTest.h>		// Include section tests
#include <list>				// Include libstdc++ STL lists
#include <iostream>			// Include file streams
#include <fstream>			// Include {if,of}streams
#include <string>			// Include string class

/* sosi - gcc 2.8.1 on Solaris 6 has trouble with public std::list<SectTest>
          as the base class.  This is not true with 2.95.2 on Irix 6.5, 
          MSVC++ 5, and egcs-2.91.60 under CygWin.  Removal of std:: here is
          fine with gcc 2.8.1 though.  Soooo, what to do?  I'll type def it
          and hope for the best....... yes this seems OK with all C++ flavors.
          Here is the story then: I have replaced the simple class declaration
          below with the typedef and the declaration using the typedef.  When
          2.8.2 is ancient I'll just set things back.        June 22 2000    */

//class ClassTest: public std::list<SectTest>

typedef std::list<SectTest> stdlistCT;

class ClassTest: public stdlistCT
{
       int         TestStatus;               // Test status (True/False)
  std::string      TestName;                 // Test name
  std::string      TestDescription;          // Test description
       int         TestRunLevel;             // Test run level
  std::vector<int> TestResults;              // Array of test results
       bool        TestType;		     // Test type flag (class vs fcts)

private:
 
// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// i                      CLASS TEST ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   CT      : A class test (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

         void CTerror(int eidx, int noret=0) const;
volatile void CTfatality(int eidx)           const;

// ____________________________________________________________________________
// ii                   CLASS TEST INITIALIZATION FUNCTIONS
// ____________________________________________________________________________

/* This function will initialize the test results vector TestResults. This      
   should be done whenever testing is done if either 1.) The status flag is
   -1, indicating no tests have ever been run, or if the current number of
   tests does not match the size of the tests results vector.  The latter
   indicates that additional tests have been added, potentially with a
   different ordering, since the test results were generated previously.    */

void SetResults(int force=0);                 

// ____________________________________________________________________________
// iii               CLASS TEST SECTION TEST INDEXING FUNCTIONS
// ____________________________________________________________________________
 
/* These function deal with indicies of Section Test tests. They allow one to
   obtain the test Pix from either an interger or a string and they perform
   range checking to insure that the requested test exists in the section.   */
 
bool                                CheckIndex(int k, int w=1) const;
std::list<SectTest>::iterator       GetPixNC(int k);
std::list<SectTest>::const_iterator GetPix(int k)              const;
std::list<SectTest>::const_iterator GetPix(const   std::string& N)    const;
int                                 GetIndex(const std::string& N)  const;


public:

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A            CLASS TEST CONSTRUCTION, DESTRUCTION, ASSIGNMENT
// ____________________________________________________________________________
 
/* Were it not for the additional variables defined in this class, all the
   constructors would just be inherited from the ANSI standard "list" class
   in the C++ library libstdc++. Since this is not the case, I am listing all
   such constructors (commented out) and will only (re)define the ones that
   I care about.

ClassTest()                               Empty Class Test
ClassTest(int N)                          Class Test w/ N Tests
ClassTest(int N, const SectTest& ST)      Class Test w/ N Section Tests
ClassTest(const ClassTest& CT)            Class Test copy of CT
ClassTest assign(N)                       Assign N element
~ClassTest()                              Destructor of Class Test           */

MSVCDLC            ClassTest();                     // Empty Class Test
MSVCDLC            ClassTest(const  ClassTest& CT); // Class Test copy of CT
MSVCDLL ClassTest& operator= (const ClassTest& CT); // Assignment operator
 
// ____________________________________________________________________________
// B                   CLASS TEST ITERATORS & MEMBER ACCESS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't
   have to keep looking them up all the time.

list<SectTest>::iterator begin()      Pointer to 1st element
list<SectTest>::iterator end()        Pointer to last element
SectTest                 front()      First element
SectTest                 back()       Last element                */

// ____________________________________________________________________________
// C                     CLASS TEST LIST & QUEUE OPERATIONS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't 
   have to keep looking them up all the time.

push_back(const  SectTest& ST)         Add ST to list end
pop_back()                             Remove ST at list end
push_front(const SectTest& ST)         Add ST to list start
pop_front(const  SectTest& ST)         Remove ST at list start
insert(iterator p, SectTest& ST)       Add ST before p
erase(iterator p)                      Remove ST at p
clear()                                Remove all list entries    */

// ____________________________________________________________________________
// D                    CLASS TEST ADDITIONAL QUEUE OPERATIONS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't 
   have to keep looking them up all the time.

int  size()                            Number of entries
bool empty()                           TRUE if pset empty
bool operator==(ClassTest CT)          TRUE if psets equal
bool operator!=(ClassTest CT)          TRUE if psets not equal    */

// ____________________________________________________________________________
// E                     CLASS TEST LIST AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 
/* These functions are "list" type of functions that have been added to make 
   the list of section tests (i.e. class test) do simple things needed for
   ready access of section tests.
   
    Function   Arguments                     Result 
   ----------  ---------  ----------------------------------------------------- 
    contains    string    Returns true/false if test of name is in class tests
       "       SectTest   Returns true/false if section test is in class tests
     seek       string    Returns iterator in class tests for test with name
       "       SectTest   Returns iterator in class tests for section test   */ 
 
MSVCDLL int contains(const std::string& tname) const;
MSVCDLL int contains(const SectTest& ST)       const;

MSVCDLL std::list<SectTest>::const_iterator seek(const std::string& tname) const;
MSVCDLL std::list<SectTest>::const_iterator seek(const SectTest& ST) const;
 
// ____________________________________________________________________________
// F                       CLASS TEST ACCESS FUNCTIONS
// ____________________________________________________________________________
 
/* These functions allow access to simple details of individual Section Tests
   or all tests in one lump sum. They will perform NO TESTING (we want code
   a bit more sophisticated for that - see a later section).  These functions
   facilitate external programs in generating nicely formatted output with the
   Class Test and Section Test results. The lower case names parallel those
   that are present in the Single Test class.

     Function                             Result
  ---------------  ---------------------------------------------------------
      GetName      Returns string for name of Class or Section Test
     GetNames      Returns string vector names of all Section Tests
  GetDescription   Returns string for description of Class or Section Test
  GetDescriptions  Returns string vector descriptions of all Section Tests
     GetStatus     Returns int for status of Class or Section Test
    GetStatuses    Returns int vector for status of all Section Tests
    GetRunLevel    Returns int for run level of Class or Section Test
   GetRunLevels    Returns int vector for run level of all Section Tests
    GetResults     Returns int vector for all current Section Test Results   */

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

/*         Input                CT : A class test (this)
                  (name)      name : The name of the test     (setting name)
                  (status)  status : Current test status      (setting status)
                  (descr)    descr : Current test description (setting describ)
                  (runlev)  runlev : Current test run level   (setting runlev)
                  (T/F)       type : Current type type        (setting type)
           Output (name)      name : The name of the test     (getting name)
                  (status)  status : Current test status      (getting status)
                  (descr)    descr : Current test description (getting describ)
                  (runlev)  runlev : Current test run level   (getting runlev)
                  (T/F)       type : Current test type        (getting type)
                              void : If setting interal value                */  
 
MSVCDLL const std::string& name()     const;                // Get name
MSVCDLL       int          status()   const;                // Get status
MSVCDLL const std::string& describe() const;                // Get descript
MSVCDLL       int          runlevel() const;                // Get runlev
MSVCDLL       bool         type()     const;                // Get type
MSVCDLL       void         name(const  std::string& Name);  // Set name
MSVCDLL       void         status(     int          Status);// Set status
MSVCDLL       void         describe(const std::string& D);  // Set decript
MSVCDLL       void         runlevel(   int          RLev);  // Set runlev
MSVCDLL       void         type(       bool            T);  // Set type
 
// ____________________________________________________________________________
// G                      CLASS TEST TESTING FUNCTIONS
// ____________________________________________________________________________
 
/* These functions will invoke any or all of the Section Tests that are listed
   as part of this Class Test. The results of the Section Tests will be stored
   in the vector TestResults and the summed results placed into TestStatus.
   The value of TestStatus (default -1) indicates the whether the tests have
   been run (as does the vector TestResults).  There are three basic means for
   doing the testing:
                       1.) Run All Of The Section Tests
                       2.) Run All Section Tests Until A Failure Occurs
                       3.) Run A Test On A Specific Section Test

   In each case the "force" flag dictates whether these tests should be run
   anew or alternatively reuse test results if previsouly run (status!=-1).  */
 
//-----------------------------------------------------------------------------
//                             Tests All Sections
//-----------------------------------------------------------------------------

MSVCDLL int TestSects(std::ostream& ostr, int force=0, int keepon=0);
 
//-----------------------------------------------------------------------------
//                  Test Specific Section (By Index Or Name)
//-----------------------------------------------------------------------------
 
MSVCDLL int TestSect(std::ostream& ostr, int tidx, int frc=0, int keepon=0);
MSVCDLL int TestSect(std::ostream& ostr, const std::string& tnam, int frc=0, int k=0); 

// ____________________________________________________________________________
// H                       CLASS TEST OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions facilitate output of class testing in a nicely formatted
   fashion. Note that when test results are listed only the status flags are
   utilized. That is, NO TESTS ARE RUN.  So, any testing should be performed
   prior to calling test output functions.
   
   Header:  This function makes a nice title for the class's section tests,
   ======   below which the PASS/FAIL results of the section tests should be
            listed. Section test PASS/FAIL listing is done in the Section Test
            class (SectTest::Result) and should conform to the format set
            herein.
   Result:  This function makes a simple single line output for listing under
   ======   ModTest headers (ModTest::Header).  The format of the output
            should conform to the format of that header.
   Results: This function loops through and spits out all the Section Test 
   =======  results of the class (SectTest::Result)
   ResRec:  Whereas the above Results function tabulates section test results,
   ======   this function will go into each section and output results from
            individual (single) tests. Typically this is done only on sections
            that have failed their test. It will also have single tests
            output results detailed if the funciton recursion level
            (nlevels) is set greater than 1.                                */

MSVCDLL std::ostream& Header(std::ostream& ostr, const std::string& MN);
MSVCDLL std::ostream& Header(std::ostream& ostr)              const;
MSVCDLL std::ostream& Result(std::ostream& ostr)              const;
MSVCDLL std::ostream& Results(std::ostream& ostr, int goon=1) const;
MSVCDLL std::ostream& ResRec(std::ostream&  ostr, int keepon, int nl=1);
 
// ____________________________________________________________________________
// I                       CLASS TEST INTERACTIVE FUNCTIONS
// ____________________________________________________________________________

MSVCDLL int AskRun(std::ostream& ostr);
 
// ____________________________________________________________________________
// J                       SECTION TEST AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 
MSVCDLL std::ostream& ListTests(std::ostream& ostr) const;

};

#endif							// ClassTest.h

