/* ModTest.cc ***************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**      Module Test                                Implementation	**
**									**
**      Scott A. Smith							**
**      Copyright (c) 2000						**
**      National High Magnetic Field Laboratory                         **
**      1800 E. Paul Dirac Drive                                        **
**      Tallahassee Florida, 32306-4005                                 **
**      email: gamma@magnet.fsu.edu                                     **
**      www: http://gamma.magnet.fsu.edu                                **
**                                                                      **
**      $Header: $
**									**
*************************************************************************/
    
/*************************************************************************
**                                                                      **
**  Description                                                         **
**                                                                      **
** Class ModTest embodies a test set for the GAMMA platform testing     **
** hierarchy. Such tests are specific to a particular GAMMA module.     **
** A module test will always deal with a single subdirectory in the     **
** GAMMA source directory, e.g. $GAMMA/src/"module". Thus this class    **
** will be intrinsically be associated with one such subdirectory and   **
** test all of the code therein.  It does so by running class tests     ** 
** where each class test is associated with a file within the module    **
** subdirectory.                                                        **
**                                                                      **
*************************************************************************/

#ifndef   ModTest_cc_			// Is this file already included?
#  define ModTest_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <Testing/ModTest.h>		// Include the interface
#include <Testing/ClassTest.h>          // Include class tests 
#include <Testing/SectTest.h>           // Include section tests
#include <Testing/SingleTest.h>         // Include single tests
#include <Basics/Gutils.h>              // Need GAMMA error messages
#include <stdlib.h>
#include <fstream>			// Include file streams
#include <string>                       // Include libstdc++ strings
#include <iostream>                     // Include input output streams (cout)
#include <vector>			// Include libstdc++ STL vectors
#include <list>				// Include libstdc++ STL lists

using std::string;			// Using libstdc++ strings
using std::list;			// Using libstdc++ lists
using std::vector;			// Using libstdc++ vectors
using std::ostream;			// Using libstdc++ output streams
using std::cout;			// Using libstdc++ standard output
using std::cin;				// Using libstdc++ standard input

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i                        MODULE TEST ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   MT      : A section test (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void ModTest::MTerror(int eidx, int noret) const
  {
  string hdr("Class Test");
  string msg;
  switch (eidx)
    {
    case 11: msg = string("Accessed Section Test Out Of Range");
             GAMMAerror(hdr,msg,noret); break;                          //(11)
    case 12: msg = string("Continuing With No Test Name");
             GAMMAerror(hdr,msg,noret); break;                          //(12)
    case 13: msg = string("Continuing With No Test Description");
             GAMMAerror(hdr,msg,noret); break;                          //(13)
    case 14: msg = string("Continuing With Status Of -1");
             GAMMAerror(hdr,msg,noret); break;                          //(14)
    case 15: msg = string("Continuing With Run Level Of 0");
             GAMMAerror(hdr,msg,noret); break;                          //(15)
    case 16: msg = string("Stopping Section Testing");
             GAMMAerror(hdr,msg,noret); break;                          //(16)
    case 17: msg = string("Unable To Run Specified Section Test");
             GAMMAerror(hdr,msg,noret); break;                          //(17)
    case 18: msg = string("Section Test By Name Does Not Exist");
             GAMMAerror(hdr,msg,noret); break;                          //(18)
    default: GAMMAerror(hdr,eidx,noret); break;
    }
  }    

volatile void ModTest::MTfatal(int eidx) const
  {                                                                 
  MTerror(eidx, eidx);                          // Output error message
  if(eidx) MTerror(0);                          // State that its fatal
  GAMMAfatal();					// Clean exit from program
  }
 
// ____________________________________________________________________________
// ii                  MODULE TEST INITIALIZATION FUNCTIONS
// ____________________________________________________________________________
 
/* This function will initialize the test results vector TestResults. This
   should be done whenever testing is done if either 1.) The status flag is
   -1, indicating no tests have ever been run, or if the current number of
   tests does not match the size of the tests results vector.  The latter
   indicates that additional tests have been added, potentially with a
   different ordering, since the test results were generated previously.    */
 
void ModTest::SetResults(int force)
  {
  if(TestResults.size() != size())              // If the vector is the wrong
    {                                           // size we make a new on
    TestResults = vector<int>(size());          //   Set the vector anew
    TestStatus = -1;                            //   Flag we are untested
    }
  if(force || TestStatus==-1)                   // We restart results if forced
    {                                           // or we have not tested yet
    for(int i=0; i<int(size()); i++)            // Initialize all results to
      TestResults[i] = -1;                      // be -1 (i.e. untested)
    }
  }

// ____________________________________________________________________________
// iii               MODULE TEST SECTION TEST INDEXING FUNCTIONS
// ____________________________________________________________________________

/* These function deal with indicies of Class Test tests. They allow one to    
   obtain the test Pix from either an interger or a string and they perform
   range checking to insure that the requested test exists in the section.   */
 
bool ModTest::CheckIndex(int k, int warn) const
  {
  if(k<0 || k>=(int)size())                     // Insure test is in range
    {
    if(warn)
       {
       MTerror(11, 1);                          // Accessed element out of range
       if(warn>1) MTerror(16);                  // Stopping class testing
       } 
    return false;
    }
  return true;
  }  
 
std::list<ClassTest>::iterator ModTest::GetPixNC(int k)
  {
  list<ClassTest>::iterator i=begin();          // Begin at class tests start
  for(int j=0; j<k && i!=end(); j++) i++;       // Move to requested test
  return i;                                     // Output Pix into list
  }
 
std::list<ClassTest>::const_iterator ModTest::GetPix(int k) const
  {
  list<ClassTest>::const_iterator i=begin();    // Begin at class tests start
  for(int j=0; j<k && i!=end(); j++) i++;       // Move to requested test
  return i;                                     // Output Pix into list
  }
 
std::list<ClassTest>::const_iterator ModTest::GetPix(const string& N) const
  {
  list<ClassTest>::const_iterator i=begin();    // Begin at class tests start
  while(i != end())                             // and loop over all of its
    {                                           // class tests
    if(N == (*i).name()) return i;              //   Return i is name match
    i++;                                        //   Else go to the next Pix
    }
  return i;                                     // Output Pix into list
  }  

int ModTest::GetIndex(const string& N) const
  {
  int k=0;                                      // String index
  list<ClassTest>::const_iterator i=begin();    // Begin at class tests start
  while(i != end())                             // and loop over all of its
    {                                           // single tests
    if(N == (*i).name()) return k;              //   Return i is name match
    i++;                                        //   Else go to the next Pix
    k++;                                        //   and increment integer
    }
  if(i==end()) k=-1;                            // Return -1 if not in list
  return k;                                     // Output integer index
  }
 
// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                   MODULE TEST CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
                                                                               
/* Were it not for the additional variables defined in this class, all the
   constructors would just be inherited from the ANSI standard "list" class
   in the C++ library libstdc++. Since this is not the case, I am listing all
   such constructors (commented out) and will only (re)define the ones that
   I care about.

ModTest()                               Empty Module Test
ModTest(int N)                          Module Test w/ N Tests
ModTest(int N, const ClassTest& CT)     Module Test w/ N Class Tests
ModTest(const ModTest& MT)              Module Test copy of MT
ModTest assign(N)                       Assign N element
~ModTest()                              Destructor of Module Test          */

ModTest::ModTest() : stdlistMT()
  {
  TestStatus      = -1;                         // Flag as untested
  TestName        = string("");                 // Set no name
  TestDescription = string("");                 // Set no description
  TestRunLevel    = 0;                          // Set run level 0
  }
 
ModTest::ModTest(const ModTest& MT) : stdlistMT(MT)
  {
  TestStatus      = MT.TestStatus;              // Copy module test status
  TestName        = MT.TestName;                // Copy module test name
  TestDescription = MT.TestDescription;         // Copy module test descript
  TestRunLevel    = MT.TestRunLevel;            // Copy module test run level
  TestResults     = MT.TestResults;             // Copy module test results
  }
 
ModTest& ModTest::operator= (const ModTest& MT)
  {
  if(this == &MT) return *this;                 // Don't copy if same one
  stdlistMT::operator=(MT);                     // Copy Class Test list
  TestStatus      = MT.TestStatus;              // Copy module test status
  TestName        = MT.TestName;                // Copy module test name
  TestDescription = MT.TestDescription;         // Copy module test descript
  TestRunLevel    = MT.TestRunLevel;            // Copy module test run level
  TestResults     = MT.TestResults;             // Copy module test results
  return *this;
  }

// ____________________________________________________________________________
// B                     MODULE TEST ITERATORS & MEMBER ACCESS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't
   have to keep looking them up all the time.

list<ClassTest>::iterator ModTest::begin()      Pointer to 1st element
list<ClassTest>::iterator ModTest::end()        Pointer to last element
ClassTest                 ModTest::front()      First element
ClassTest                 ModTest::back()       Last element                */

// ____________________________________________________________________________
// C                     MODULE TEST LIST & QUEUE OPERATIONS
// ____________________________________________________________________________ 
 
/* These functions are inherited from the ANSI standard "list" class in the C++ 
   library libstdc++.  I am listing them here so that I & other users don't  
   have to keep looking them up all the time. 
 
ModTest::push_back(const  ClassTest& ST)       Add CT to list end
ModTest::pop_back()                            Remove CT at list end
ModTest::push_front(const ClassTest& ST)       Add CT to list start
ModTest::pop_front(const  ClassTest& ST)       Remove CT at list start
 
ModTest::insert(iterator p, ClassTest& ST)     Add CT before p
ModTest::erase(iterator p)                     Remove CT at p
ModTest::clear()                               Remove all list entries       */

// ____________________________________________________________________________ 
// D                    MODULE TEST ADDITIONAL QUEUE OPERATIONS
// ____________________________________________________________________________ 
 
/* These functions are inherited from the ANSI standard "list" class in the C++ 
   library libstdc++.  I am listing them here so that I & other users don't  
   have to keep looking them up all the time. 
 
int  ModTest::size()                            Number of entries
bool ModTest::empty()                           TRUE if pset empty
bool ModTest::operator==(ModTest MT)            TRUE if psets equal
bool ModTest::operator!=(ModTest MT)            TRUE if psets not equal      */

// ____________________________________________________________________________
// E                     MODULE TEST LIST AUXILIARY FUNCTIONS
// ____________________________________________________________________________

/* These functions are "list" type of functions that have been added to make
   the list of class tests (i.e. module test) do simple things needed for
   ready access of class tests.
  
    Function   Arguments                     Result
   ----------  ---------  -----------------------------------------------------
    contains    string    Returns true/false if test of name is in module tests
       "       ClassTest  Returns true/false if class test is in module tests
     seek       string    Returns iterator in module tests for test with name
       "       ClassTest  Returns iterator in module tests for class test    */

int ModTest::contains(const string& testname) const
  { return (seek(testname) != end()); }
int ModTest::contains(const ClassTest& CT) const
  { return contains(CT.name()); }

std::list<ClassTest>::const_iterator ModTest::seek(const string& tname) const
  {
  ClassTest CT;                                 // Working class test
  list<ClassTest>::const_iterator item=begin(); // Begin at input CT start
  while(item != end())                          // and loop over all of its
    {                                           // class tests
    CT = *item;                                 //   This is a class test
    if(tname == CT.name()) return item;         //   Return item is name match
    item++;                                     //   Go to the next Pix
    }
  return item;
  }

std::list<ClassTest>::const_iterator ModTest::seek(const ClassTest& CT) const
  { return seek(CT.name()); }

// ____________________________________________________________________________
// F                       MODULE TEST ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions allow access to simple details of individual Class Tests
   or all tests in one lump sum. They will perform NO TESTING (we want code
   a bit more sophisticated for that - see a later section).  These functions
   facilitate external programs in generating nicely formatted output with the
   Module Test and Class Test results. The lower case names parallel those
   that are present in the Section and Single Test classes.

     Function                             Result
  ---------------  ---------------------------------------------------------
      GetName      Returns string for name of Module or Class Test
     GetNames      Returns string vector names of all Class Tests
  GetDescription   Returns string for description of Module or Class Test
  GetDescriptions  Returns string vector descriptions of all Class Tests
     GetStatus     Returns int for status of Module or Class Test
    GetStatuses    Returns int vector for status of all Class Tests
    GetRunLevel    Returns int for run level of Module or Class Test
   GetRunLevels    Returns int vector for run level of all Class Tests
    GetResults     Returns int vector for all current Class Test Results     */
 
string ModTest::GetName() const
  { return TestName; }
 
string ModTest::GetName(int k) const
  {
  if(!CheckIndex(k,1))                          // Insure test is in range
    { MTerror(12); return string(""); }         // If not, go on w/ no range
  list<ClassTest>::const_iterator i=begin();    // Begin at class tests start
  for(int j=0; j<k; j++) i++;                   // Move to requested test
  return (*i).name();                           // Output test name
  }
 
vector<string> ModTest::GetNames() const
  {
  vector<string> Names;
  list<ClassTest>::const_iterator i=begin();    // Begin at class tests start
  while(i != end())                             // Loop until end of tests
    { Names.push_back((*i).name()); i++; }      // while putting names in vector
  return Names;
  } 

    
string ModTest::GetDescription() const
  { return TestDescription; }
    
string ModTest::GetDescription(int k) const
  {
  if(!CheckIndex(k,1))                          // Insure test is in range
    { MTerror(13); return string(""); }         // If not, go on w/ no descript
  list<ClassTest>::const_iterator i=begin();    // Begin at class tests start
  for(int j=0; j<k; j++) i++;                   // Move to requested test
  return (*i).describe();                       // Output test description
  }
 
vector<string> ModTest::GetDescriptions() const
  {
  vector<string> Descripts;
  list<ClassTest>::const_iterator i=begin();    // Begin at class tests start
  while(i != end())                             // Loop until end of tests
    { Descripts.push_back((*i).describe()); i++; }
  return Descripts;
  }
 
int ModTest::GetStatus() const
  { return TestStatus; }
 
int ModTest::GetStatus(int k) const
  {
  if(!CheckIndex(k,1))                          // Insure test is in range
    { MTerror(14); return -1; }                 // If not, go on w/ status -1
  list<ClassTest>::const_iterator i=begin();    // Begin at class tests start
  for(int j=0; j<k; j++) i++;                   // Move to requested test
  return (*i).status();                         // Output test status
  }
 
vector<int> ModTest::GetStatuses() const
  {
  vector<int> Stati;
  list<ClassTest>::const_iterator i=begin();     // Begin at class tests start
  while(i != end())                             // Loop until end of tests
    { Stati.push_back((*i).status()); i++; }
  return Stati;
  }

int ModTest::GetRunLevel() const
  { return TestRunLevel; }
 
int ModTest::GetRunLevel(int k) const
  {
  if(!CheckIndex(k,1))                          // Insure test is in range
    { MTerror(15); return 0; }                  // If not, go on w/ level 0
  list<ClassTest>::const_iterator i=begin();    // Begin at class tests start
  for(int j=0; j<k; j++) i++;                   // Move to requested test
  return (*i).runlevel();                       // Output test run level
  }
 
vector<int> ModTest::GetRunLevels() const
  {
  vector<int> RLevs;
  list<ClassTest>::const_iterator i=begin();    // Begin at class tests start
  while(i != end())                             // Loop until end of tests
    { RLevs.push_back((*i).runlevel()); i++; }
  return RLevs;
  } 
    
vector<int> ModTest::GetResults() { SetResults(); return TestResults; }
 
 
/*         Input                MT : A module test (this)
                  (name)      name : The name of the test     (setting name)
                  (status)  status : Current test status      (setting status)
                  (descr)    descr : Current test description (setting describ)
                  (runlev)  runlev : Current test run level   (setting runlev)
           Output (name)      name : The name of the test     (getting name)
                  (status)  status : Current test status      (getting status)
                  (descr)    descr : Current test description (getting describ)
                  (runlev)  runlev : Current test run level   (getting runlev)
                              void : If setting interal value                */  
 
const string& ModTest::name()     const  { return TestName; }
      int     ModTest::status()   const  { return TestStatus; }
const string& ModTest::describe() const  { return TestDescription; }
      int     ModTest::runlevel() const  { return TestRunLevel; }
 
void ModTest::name(const     string& N) { TestName        = N; }
void ModTest::status(        int     S) { TestStatus      = S; }
void ModTest::describe(const string& D) { TestDescription = D; }
void ModTest::runlevel(      int     R) { TestRunLevel    = R; }

// ____________________________________________________________________________
// G                        MODULE TEST TESTING FUNCTIONS
// ____________________________________________________________________________

/* These functions will invoke any or all of the Class Tests that are listed
   as part of this Class Test. Results of the Class Tests will be stored
   in the vector TestResults and the summed results placed into TestStatus.
   The value of TestStatus (default -1) indicates the whether the tests have
   been run (as does the vector TestResults).  There are three basic means for
   testing:
             1.) Run All Of The Class Tests			TestClasses
             2.) Run All Class Tests Until A Failure Occurs     TestClass
             3.) Run A Test On A Specific Class Test            TestClass
  
   In each case the "force" flag dictates whether these tests should be run
   anew or alternatively reuse test results if previsouly run (status!=-1).  */
 
//-----------------------------------------------------------------------------
//                            Tests All Classes
//-----------------------------------------------------------------------------
 
int ModTest::TestClasses(ostream& ostr, int anew, int keepon)
  {
  SetResults(anew);                             // Initialize results(if want)
  list<ClassTest>::iterator i=begin();          // Begin at class tests start
  int tidx = 0;                                 // This is plain test index
  TestStatus = 1;                               // Assume we test OK
  int goodtest = 1;                             // A flag if we fail
  while(i!=end() && goodtest)                   // Loop until end of tests
    {
    (*i).runlevel(TestRunLevel);                // Set test run level
    goodtest=((*i).TestSects(ostr,anew,keepon));//   Run class test
    TestResults[tidx++] = goodtest;             //   Store test result
    TestStatus *= abs(goodtest);                //   Sum result for section
    i++;                                        //   Increment to next test
    if(keepon) goodtest=1;                      //   Continue on failure
    }                                           ///     (if desired)
  return TestStatus;
  }
 
//-----------------------------------------------------------------------------
//                   Test Specific Class (By Index Or Name)
//-----------------------------------------------------------------------------
 
int ModTest::TestClass(ostream& ostr, const int tidx, int frc, int keepon)
  {
  if(!CheckIndex(tidx, 1))                      // Insure test tidx exists
    { MTfatal(17); }                            // Cannot run requested test
  list<ClassTest>::iterator item;               // A Pix into list of tests
  item = GetPixNC(tidx);                        // Get Pix for requested test
  return ((*item).TestSects(ostr,frc,keepon));	// Run test, return result
  }

int ModTest::TestClass(ostream& ostr, const string& tname, int frc, int keepon)
  {
  int tidx = GetIndex(tname);                   // Get text index
  if(!CheckIndex(tidx, 1))                      // Insure test tname exists
    { MTerror(18); MTfatal(17); }               // Cannot run requested test
  list<ClassTest>::iterator item;               // A Pix into list of tests
  item = GetPixNC(tidx);                        // Get Pix for requested test
  return ((*item).TestSects(ostr,frc,keepon)); 	// Run test, return result
  }
 
// ____________________________________________________________________________
// H                       MODULE TEST OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions facilitate output of module testing in a nicely formatted
   fashion. Note that when test results are listed only the status flags are
   utilized. That is, NO TESTS ARE RUN.  So, any testing should be performed
   prior to calling test output functions.
   
   Header:  This function makes a nice title for the module's class tests,
   ======   below which the PASS/FAIL results of the class tests should be
            listed. Class test PASS/FAIL listing is done in the Class Test
            class (ClassTest::Result) and should conform to the format set
            herein.
   Result:  This function makes a simple single line output for listing under
   ======   GAMMA Test headers (GamTest::Header).  The format of the output
            should conform to the format of that header.
   Results: This function loops through and spits out all the Class Test
   =======  results of the module (ClassTest::Result)
   ResRec:  Whereas the above Results function tabulates class test results,
   ======   this function will go into each class and output results from
            individual secton tests. Typically this is done only on classes
            that have failed their test. It will also have section tests
            rall results from single tests if the funciton recursion level
            (nlevels) is set greater than 1.                                */


ostream& ModTest::Header(ostream& ostr) const
  {
  string TName = TestName;
  if(!TName.length()) TName = string("Unknown");
  int len = 44 + TName.length();
  int ls = 42-len/2;
  string s1 = string(ls, ' ');
  ostr << "\n\n" << s1;
  ostr << "Testing GAMMA "<< TName << " Module Classes & Functions\n";
  ostr << s1 << string(13+TName.length()+28, '*');
  ostr << "\n\n";
  ostr << "  Class/File                        Description                        Outcome\n";
  ostr << " ============  ======================================================  =======\n";
  return ostr;
  }

ostream& ModTest::Result(ostream& ostr) const
  {        
  string CName = TestName;
  if(!CName.length()) CName = string("Unknown");
  string CDesc = TestDescription;
  if(!CDesc.length()) CDesc = string("Unknown");

  ostr << " " << CName;
  int len = CName.length();
  ostr << string(12-len, ' ');
  ostr << "  "
       << CDesc;
  len = CDesc.length();
  ostr << string(54-len, ' ');
  ostr << "    ";
  if(TestStatus>0)      ostr << "PASS";
  else if(!TestStatus)  ostr << "FAIL";
  else                  ostr << "----";
  ostr << "\n";
  return ostr;
  }

ostream& ModTest::Results(ostream& ostr, int keepon) const
  {
  list<ClassTest>::const_iterator item;         // Iterator in Class test list
  bool goon = true;				// Flag to continue output
  item = begin();                               // Begin at section tests start
  while(item != end() && goon)                  // Loop until end of tests &
    {						// output class test results
    (*item).Result(ostr);			//   Output class result
    if(!(*item).status() && !keepon) goon=false;//   Quit if fail & not go on
    item++;					//   Move to next class
    }
  return ostr;
  }

ostream& ModTest::ResRec(ostream& ostr, int keepon, int nlevels)
  {
  list<ClassTest>::iterator clpix;             // Pix in class test list
  clpix = begin();                             // Begin at 1st class test
  int clstat;                                  // Class test status
  while(clpix != end())                        // Loop until end of tests
    { 
    clstat = (*clpix).status();                //   Class test result
    if(!clstat)                                //   If class test failed
      {                                        //   then we will look at
      (*clpix).Header(ostr, TestName);         //   the class's section test
      (*clpix).Results(ostr, keepon);          //   results
      if(nlevels>1)                            //   If desired we can also
        (*clpix).ResRec(ostr,keepon,nlevels--);//   results spit out class 
      }                                        //   section results
    clpix++;
    if(!keepon  && !clstat) clpix = end();  //     Bail on fail if !keepon
    } 
  return ostr;
  }

// ____________________________________________________________________________
// I                      MODULE TEST INTERACTIVE FUNCTIONS
// ____________________________________________________________________________

/* These functions facilitate the running a module test. Recall that each 
   module test runs through a series of class tests, ultimately the test 
   returning TRUE/FALSE on completion. However there are several additional 
   facets of module testing.

   1.) Continue Through Failed Tests: If a class test fails should testing
                                      stop or continue into the next class
   3.) Recurse Through Section Tests: If a class test fails should testing
                                      recurse into section tests for that
                                      class, displaying specific failures?
   
      Recursion Levels:   0 = Output Class   Test Pass/Fail Results
                          1 = Output Section Test Pass/Fail Results
                          2 = Output Single  Test Pass/Fail Results
                          3 = Output Single  Test Full Details

   These functions assume you want some output text from running the tests,
   otherwise you would just use the TestLevel function directly, right?      */  

int ModTest::AskRun(ostream& ostr)
  {
  int keepon=0;
  string yn;
  cout << "\n\t" << "Proceed Through Failed Tests [y,n]? ";
  cin >> yn;
  if(yn == "y") keepon=1;
  int recurse = 1;
  cout << "\n\t" << "Failed Class Recursion Levels [0,3]? ";
  cin >> recurse;
 
  TestRunLevel = 0;                             // Set testing silent
  TestClasses(ostr, 0, keepon);                 // Run the class tests
  Header(ostr);                                 // Output module header
  Results(ostr, keepon);                        // Output module results
  if(!TestStatus && recurse>=1)                 // If any class test failed
    ResRec(ostr, keepon, recurse--);		// examine it in detail
  return TestStatus;
  }  

// ____________________________________________________________________________
// J                       MODULE TEST AUXILIARY FUNCTIONS
// ____________________________________________________________________________

ostream& ModTest::ListTests(ostream& ostr) const
  {
  list<ClassTest>::const_iterator i=begin();     // Begin at class tests start
   while(i != end())                             // Loop until end of tests
     {   ostr << (*i).describe() << "\n";
     i++;
     }
  return ostr;
  }

ostream& ModTest::FinishTest(ostream& ostr) const
  {
  ostr << "\n\n";
  if(TestStatus>=1)
    {
    ostr << string(15, ' ')
         << "GAMMA Module " << TestName << " Has Successfully Passed Its Tests";
    ostr << "\n" << string(15, ' ') << string(47+TestName.length(), '=');
    ostr << "\n" << string(15, ' ') << string(47+TestName.length(), '=');
    }
  else if(!TestStatus)
    {
    ostr << string(19, ' ')
         << "GAMMA Has Failed Testing. Please Report This!";
    ostr << "\n" << string(19, ' ') << string(45, '=');
    ostr << "\n" << string(19, ' ') << string(45, '=') << "\n\n";
    ostr << "\tIf you would be so kind, it's a simple matter to report these results.\n";
    ostr << "\tSimply \"Cut And Paste\" the output above into a GAMMA Bug Report Form\n";
    ostr << "\ton WWW or into an e-mail message. For the former, use the form at\n\n";
    ostr << "\t\t\thttp:://gamma.magnet.fsu.edu/info/bugrep/\n\n";
    ostr << "\tThis form will also be available at GAMMA mirror sites and in your local\n";
    ostr << "\tGAMMA documentation if it has been installed (Information:Bug Report).\n";
    ostr << "\tFor reporting via e-mail, send the information to\n\n";
    ostr << "\t\t\t\tgamma@magnet.fsu.edu\n\n";
    ostr << "\tIn either case, please include details about your system and C++ compiler\n";
    ostr << "\tand your e-mail address if you desire a response. The effort is appreciated.\n";
    ostr << "\n\n";
    }
  else
    ostr << string(19, ' ')
         << "GAMMA Module Testing Has Not Been Completed";
  ostr << "\n\n";
  return ostr;
  }



#endif							// ModTest.cc
