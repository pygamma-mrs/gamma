/* GamTest.cc ***************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**      GAMMA Test                                Implementation	**
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
** Class GamTest embodies a test set for the GAMMA platform testing     **
** hierarchy. Such tests are specific to a particular set of GAMMA      **
** modules. A GAMMA test will always deal with set of modules, or       **
** equivalently, a set of source subdirectories from the main GAMMA     **
** source directory, $GAMMA/src/. Thus, this class allows for facile    **
** construction of bulk GAMMA testing. It does so by running module     **
** tests where each module test is associated with a subdirectory from  **
** the main GAMMA source directory.                                     **
**                                                                      **
*************************************************************************/

#ifndef   GamTest_cc_			// Is this file already included?
#  define GamTest_cc_ 1			// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <Testing/GamTest.h>		// Include the interface
#include <Testing/ModTest.h> 	        // Include module tests 
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
//using std::cout;			// Using libstdc++ standard output
using std::cin;				// Using libstdc++ standard input

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i                        GAMMA TEST ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   GT      : A GAMMA test (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void GamTest::GTerror(int eidx, int noret) const
  {
  string hdr("GAMMA Test");
  string msg;
  switch (eidx)
    {
    case 11: msg = string("Accessed Module Test Out Of Range");
             GAMMAerror(hdr,msg,noret); break;                          //(11)
    case 12: msg = string("Continuing With No Test Name");
             GAMMAerror(hdr,msg,noret); break;                          //(12)
    case 13: msg = string("Continuing With No Test Description");
             GAMMAerror(hdr,msg,noret); break;                          //(13)
    case 14: msg = string("Continuing With Status Of -1");
             GAMMAerror(hdr,msg,noret); break;                          //(14)
    case 15: msg = string("Continuing With Run Level Of 0");
             GAMMAerror(hdr,msg,noret); break;                          //(15)
    case 16: msg = string("Stopping Module Testing");
             GAMMAerror(hdr,msg,noret); break;                          //(16)
    case 17: msg = string("Unable To Run Specified Module Test");
             GAMMAerror(hdr,msg,noret); break;                          //(17)
    case 18: msg = string("Module Test By Name Does Not Exist");
             GAMMAerror(hdr,msg,noret); break;                          //(18)
    case 19: msg = string("Cannot Output Module's Class Test Results");
             GAMMAerror(hdr,msg,noret); break;                          //(19)
    default: GAMMAerror(hdr,eidx,noret); break;
    }
  }    

volatile void GamTest::GTfatal(int eidx) const
  {                                                                 
  GTerror(eidx, eidx);                          // Output error message
  if(eidx) GTerror(0);                          // State that its fatal
  GAMMAfatal();					// Clean exit from program
  }
 
// ____________________________________________________________________________
// ii                  GAMMA TEST INITIALIZATION FUNCTIONS
// ____________________________________________________________________________
 
/* This function will initialize the test results vector TestResults. This
   should be done whenever testing is done if either 1.) The status flag is
   -1, indicating no tests have ever been run, or if the current number of
   tests does not match the size of the tests results vector.  The latter
   indicates that additional tests have been added, potentially with a
   different ordering, since the test results were generated previously.    */
 
void GamTest::SetResults(int force)
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
// iii               GAMMA TEST SECTION TEST INDEXING FUNCTIONS
// ____________________________________________________________________________

/* These function deal with indicies of Module Test tests. They allow one to    
   obtain the test Pix from either an interger or a string and they perform
   range checking to insure that the requested test exists in the section.   */
 
bool GamTest::CheckIndex(int k, int warn) const
  {
  if(k<0 || k>=(int)size())                     // Insure test is in range
    {
    if(warn)
       {
       GTerror(11, 1);                          // Accessed element out of range
       if(warn>1) GTerror(16);                  // Stopping module testing
       } 
    return false;
    }
  return true;
  }  
 
std::list<ModTest>::iterator GamTest::GetPixNC(int k)
  {
  list<ModTest>::iterator i=begin();            // Begin at module tests start
  for(int j=0; j<k && i!=end(); j++) i++;       // Move to requested test
  return i;                                     // Output Pix into list
  }
 
std::list<ModTest>::const_iterator GamTest::GetPix(int k) const
  {
  list<ModTest>::const_iterator i=begin();      // Begin at module tests start
  for(int j=0; j<k && i!=end(); j++) i++;       // Move to requested test
  return i;                                     // Output Pix into list
  }
 
std::list<ModTest>::const_iterator GamTest::GetPix(const string& N) const
  {
  list<ModTest>::const_iterator i=begin();      // Begin at module tests start
  while(i != end())                             // and loop over all of its
    {                                           // module tests
    if(N == (*i).name()) return i;              //   Return i is name match
    i++;                                        //   Else go to the next Pix
    }
  return i;                                     // Output Pix into list
  }  

int GamTest::GetIndex(const string& N) const
  {
  int k=0;                                      // String index
  list<ModTest>::const_iterator i=begin();      // Begin at module tests start
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
// A                   GAMMA TEST CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
                                                                               
/* Were it not for the additional variables defined in this class, all the
   constructors would just be inherited from the ANSI standard "list" class
   in the C++ library libstdc++. Since this is not the case, I am listing all
   such constructors (commented out) and will only (re)define the ones that
   I care about.

GamTest()                               Empty Module Test
GamTest(int N)                          Module Test w/ N Tests
GamTest(int N, const ModTest& MT)       Module Test w/ N Module Tests
GamTest(const GamTest& GT)              Module Test copy of GT
GamTest assign(N)                       Assign N element
~GamTest()                              Destructor of Module Test          */

GamTest::GamTest() : stdlistGT()
  {
  TestStatus      = -1;                         // Flag as untested
  TestName        = string("");                 // Set no name
  TestDescription = string("");                 // Set no description
  TestRunLevel    = 0;                          // Set run level 0
  }
 
GamTest::GamTest(const GamTest& GT) : stdlistGT(GT)
  {
  TestStatus      = GT.TestStatus;              // Copy GAMMA test status
  TestName        = GT.TestName;                // Copy GAMMA test name
  TestDescription = GT.TestDescription;         // Copy GAMMA test descript
  TestRunLevel    = GT.TestRunLevel;            // Copy GAMMA test run level
  TestResults     = GT.TestResults;             // Copy GAMMA test results
  }
 
GamTest& GamTest::operator= (const GamTest& GT)
  {
  if(this == &GT) return *this;                 // Don't copy if same one
  stdlistGT::operator=(GT);                     // Copy Module Test list
  TestStatus      = GT.TestStatus;              // Copy GAMMA test status
  TestName        = GT.TestName;                // Copy GAMMA test name
  TestDescription = GT.TestDescription;         // Copy GAMMA test descript
  TestRunLevel    = GT.TestRunLevel;            // Copy GAMMA test run level
  TestResults     = GT.TestResults;             // Copy GAMMA test results
  return *this;
  }

// ____________________________________________________________________________
// B                     GAMMA TEST ITERATORS & MEMBER ACCESS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't
   have to keep looking them up all the time.

std::list<ModTest>::iterator GamTest::begin()      Pointer to 1st element
std::list<ModTest>::iterator GamTest::end()        Pointer to last element
ModTest                 GamTest::front()      First element
ModTest                 GamTest::back()       Last element                   */

// ____________________________________________________________________________
// C                     GAMMA TEST LIST & QUEUE OPERATIONS
// ____________________________________________________________________________ 
 
/* These functions are inherited from the ANSI standard "list" class in the C++ 
   library libstdc++.  I am listing them here so that I & other users don't  
   have to keep looking them up all the time. 
 
GamTest::push_back(const  ModTest& MT)          Add MT to list end
GamTest::pop_back()                             Remove MT at list end
GamTest::push_front(const ModTest& MT)          Add MT to list start
GamTest::pop_front(const  ModTest& MT)          Remove MT at list start
 
GamTest::insert(iterator p, ModTest& MT)        Add MT before p
GamTest::erase(iterator p)                      Remove MT at p
GamTest::clear()                                Remove all list entries      */

// ____________________________________________________________________________ 
// D                    GAMMA TEST ADDITIONAL QUEUE OPERATIONS
// ____________________________________________________________________________ 
 
/* These functions are inherited from the ANSI standard "list" class in the C++ 
   library libstdc++.  I am listing them here so that I & other users don't  
   have to keep looking them up all the time. 
 
int  GamTest::size()                            Number of entries
bool GamTest::empty()                           TRUE if pset empty
bool GamTest::operator==(GamTest GT)            TRUE if psets equal
bool GamTest::operator!=(GamTest GT)            TRUE if psets not equal      */

// ____________________________________________________________________________
// E                     GAMMA TEST LIST AUXILIARY FUNCTIONS
// ____________________________________________________________________________

/* These functions are "list" type of functions that have been added to make
   the list of Module Tests (i.e. GAMMA test) do simple things needed for
   ready access of Module Tests.
  
    Function   Arguments                     Result
   ----------  ---------  -----------------------------------------------------
    contains    string    Returns true/false if test of name is in GAMMA tests
       "        ModTest   Returns true/false if module test is in GAMMA tests
     seek       string    Returns iterator in GAMMA tests for test with name
       "        ModTest   Returns iterator in GAMMA tests for module test    */

int GamTest::contains(const string& testname) const
  { return (seek(testname) != end()); }
int GamTest::contains(const ModTest& MT) const
  { return contains(MT.name()); }

std::list<ModTest>::const_iterator GamTest::seek(const string& tname) const
  {
  ModTest MT;                                   // Working module test
  list<ModTest>::const_iterator item=begin();   // Begin at input MT start
  while(item != end())                          // and loop over all of its
    {                                           // module tests
    MT = *item;                                 //   This is a module test
    if(tname == MT.name()) return item;         //   Return item is name match
    item++;                                     //   Go to the next Pix
    }
  return item;
  }

std::list<ModTest>::const_iterator GamTest::seek(const ModTest& MT) const
  { return seek(MT.name()); }

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
     GetStatus     Returns int for status of GAMMA or Module Test
    GetStatuses    Returns int vector for status of all Module Tests
    GetRunLevel    Returns int for run level of GAMMA or Module Test
   GetRunLevels    Returns int vector for run level of all Module Tests
    GetResults     Returns int vector for all current Module Test Results     */
 
string GamTest::GetName() const
  { return TestName; }
 
string GamTest::GetName(int k) const
  {
  if(!CheckIndex(k,1))                          // Insure test is in range
    { GTerror(12); return string(""); }         // If not, go on w/ no range
  list<ModTest>::const_iterator i=begin();      // Begin at module tests start
  for(int j=0; j<k; j++) i++;                   // Move to requested test
  return (*i).name();                           // Output test name
  }
 
vector<string> GamTest::GetNames() const
  {
  vector<string> Names;
  list<ModTest>::const_iterator i=begin();      // Begin at module tests start
  while(i != end())                             // Loop until end of tests
    { Names.push_back((*i).name()); i++; }      // while putting names in vector
  return Names;
  } 

    
string GamTest::GetDescription() const
  { return TestDescription; }
    
string GamTest::GetDescription(int k) const
  {
  if(!CheckIndex(k,1))                          // Insure test is in range
    { GTerror(13); return string(""); }         // If not, go on w/ no descript
  list<ModTest>::const_iterator i=begin();      // Begin at module tests start
  for(int j=0; j<k; j++) i++;                   // Move to requested test
  return (*i).describe();                       // Output test description
  }
 
vector<string> GamTest::GetDescriptions() const
  {
  vector<string> Descripts;
  list<ModTest>::const_iterator i=begin();      // Begin at module tests start
  while(i != end())                             // Loop until end of tests
    { Descripts.push_back((*i).describe()); i++; }
  return Descripts;
  }
 
int GamTest::GetStatus() const
  { return TestStatus; }
 
int GamTest::GetStatus(int k) const
  {
  if(!CheckIndex(k,1))                          // Insure test is in range
    { GTerror(14); return -1; }                 // If not, go on w/ status -1
  list<ModTest>::const_iterator i=begin();      // Begin at module tests start
  for(int j=0; j<k; j++) i++;                   // Move to requested test
  return (*i).status();                         // Output test status
  }
 
vector<int> GamTest::GetStatuses() const
  {
  vector<int> Stati;
  list<ModTest>::const_iterator i=begin();      // Begin at module tests start
  while(i != end())                             // Loop until end of tests
    { Stati.push_back((*i).status()); i++; }
  return Stati;
  }

int GamTest::GetRunLevel() const
  { return TestRunLevel; }
 
int GamTest::GetRunLevel(int k) const
  {
  if(!CheckIndex(k,1))                          // Insure test is in range
    { GTerror(15); return 0; }                  // If not, go on w/ level 0
  list<ModTest>::const_iterator i=begin();      // Begin at module tests start
  for(int j=0; j<k; j++) i++;                   // Move to requested test
  return (*i).runlevel();                       // Output test run level
  }
 
vector<int> GamTest::GetRunLevels() const
  {
  vector<int> RLevs;
  list<ModTest>::const_iterator i=begin();      // Begin at module tests start
  while(i != end())                             // Loop until end of tests
    { RLevs.push_back((*i).runlevel()); i++; }
  return RLevs;
  } 
    
vector<int> GamTest::GetResults() { SetResults(); return TestResults; }
 
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
 
const string& GamTest::name()     const  { return TestName; }
      int     GamTest::status()   const  { return TestStatus; }
const string& GamTest::describe() const  { return TestDescription; }
      int     GamTest::runlevel() const  { return TestRunLevel; }
 
void GamTest::name(const     string& N) { TestName        = N; }
void GamTest::status(        int     S) { TestStatus      = S; }
void GamTest::describe(const string& D) { TestDescription = D; }
void GamTest::runlevel(      int     R) { TestRunLevel    = R; }

// ____________________________________________________________________________
// G                        GAMMA TEST TESTING FUNCTIONS
// ____________________________________________________________________________

/* These functions will invoke any or all of the Module Tests that are listed
   as part of this GAMMA Test. Results of the Module Tests will be stored
   in the vector TestResults and the summed results placed into TestStatus.
   The value of TestStatus (default -1) indicates the whether the tests have
   been run (as does the vector TestResults).  There are three basic means for
   doing the testing:
                       1.) Run All Of The Module Tests
                       2.) Run All Module Tests Until A Failure Occurs
                       3.) Run A Test On A Specific Module Test
  
   In each case the "force" flag dictates whether these tests should be run
   anew or alternatively reuse test results if previsouly run (status!=-1).  */
 
//-----------------------------------------------------------------------------
//                              Test All Modules
//-----------------------------------------------------------------------------
 
int GamTest::TestMods(ostream& ostr, int anew, int keepon)
  {
  SetResults(anew);                             // Initialize results(if want)
  list<ModTest>::iterator i=begin();            // Begin at module tests start
  int tidx = 0;                                 // This is plain test index
  TestStatus = 1;                               // Assume we test OK
  int goodtest = 1;                             // A flag if we fail
  while(i!=end() && goodtest)                   // Loop until end of tests
    {
    (*i).runlevel(TestRunLevel);                // Set test run level
    goodtest=((*i).TestClasses(ostr,anew,keepon));//   Run module test
    TestResults[tidx++] = goodtest;             //   Store test result
    TestStatus *= abs(goodtest);                //   Sum result for GAMMA test
    i++;                                        //   Increment to next test
    if(keepon) goodtest=1;                      //   Continue on failure
    }                                           //     (if desired)
  return TestStatus;
  }

//-----------------------------------------------------------------------------
//                     Run Specific Module (By Name Or Index)
//-----------------------------------------------------------------------------
 
int GamTest::TestMod(ostream& ostr, int tidx, int anew, int keepon)
  {
  if(!CheckIndex(tidx, 1))                        // Insure test tidx exists
    { GTfatal(17); }                              // Cannot run requested test
  list<ModTest>::iterator item;                   // A Pix into list of tests
  item = GetPixNC(tidx);                          // Get Pix for requested test
  return ((*item).TestClasses(ostr,anew,keepon)); // Run test, return result
  }

int GamTest::TestMod(ostream& ostr, const string& tname, int anew, int  keepon)
  {
  int tidx = GetIndex(tname);                     // Get text index
  if(!CheckIndex(tidx, 1))                        // Insure test tname exists
    { GTerror(18); GTfatal(17); }                 // Cannot run requested test
  list<ModTest>::iterator item;                   // A Pix into list of tests
  item = GetPixNC(tidx);                          // Get Pix for requested test
  return ((*item).TestClasses(ostr,anew,keepon)); // Run test, return result
  }
 
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
   =======  results of the GAMMA test (ModTest::Result) 
   ResRec:  Whereas the above Results function tabulates module test results,
   ======   this function will go into each module and output results from
            individual class tests. Typically this is done only on modules
            that have failed their test. It will also have class tests call
            results from individual section tests if the function resursion
            level (nlevels) is set greater than 1.                           */

ostream& GamTest::Header(ostream& ostr) const
  {
  string TName = TestName;
  if(!TName.length()) TName = string("UNKNOWN BULK");
  int len = 6 + TName.length() + 8;
  int ls = 40-len/2;
  string s1 = string(ls, ' ');
  ostr << "\n\n" << s1 << "GAMMA "<< TName << " TESTING\n";
  int tlen = TestDescription.length();          // Length of test description
  if(tlen) 			                // If there is a description,
     {                                          // include it below the title
     ostr << string(38-tlen/2, ' ')             //        (centered)
          << "(" << TestDescription << ")\n";
     }
  ostr << " -----------------------------------------------------------------------------\n";
  ostr << " -----------------------------------------------------------------------------\n";
  ostr << "\n";
  ostr << "    Module                           Description                       Outcome\n";
  ostr << " ============  ======================================================  =======\n";
  return ostr;
  }

ostream& GamTest::Results(ostream& ostr, int keepon) const
  {
  list<ModTest>::const_iterator item;           // Iterator in Module test list
  bool goon = true;				// Flag to continue output
  item = begin();                               // Begin at section tests start
  while(item != end() && goon)                  // Loop until end of tests &
   {						// output module test results
   (*item).Result(ostr);			//   Output module result
   if(!(*item).status() && !keepon) goon=false;	//   Quit if fail & not go on
   item++; 					//   Move to next module
   }
  return ostr;
  }

ostream& GamTest::ResRec(ostream& ostr, int keepon, int nlevels)
  {
  list<ModTest>::iterator modpix;		// Pix in module test list
  modpix = begin();				// Begin at 1st module test
  int modstat;					// Module test status
  while(modpix != end())			// Loop until end of tests
     { 
     modstat = (*modpix).status();		//   Module test result
     if(!modstat)				//   If module test failed
        {					//   then we will look at
        (*modpix).Header(ostr);			//   the module's class test
        (*modpix).Results(ostr, keepon);	//   results
        if(nlevels>1)				//   If desired we can also
        (*modpix).ResRec(ostr,keepon,nlevels--);//   results spit out class 
        }					//   section results
     modpix++;
     if(!keepon  && !modstat) modpix = end();  //     Bail on fail if !keepon
     } 
  return ostr;
  }

// ____________________________________________________________________________
// I                       GAMMA TEST INTERACTIVE FUNCTIONS
// ____________________________________________________________________________
 
/* These functions facilitate the running of a GAMMA test. Recall that each 
   test runs through a series of module tests, ultimately the test returning
   TRUE/FALSE on completion. However there are several additional facets of
   GAMMA testing.
 
   1.) Continue Through Failed Tests:  If a module test fails should testing
                                       stop or continue on to the next module
   2.) Recurse Through Class Tests:    If a module test fails should testing
                                       recurse into class tests for that
                                       module , displaying specific failures?

      Recursion Levels:    0 = Output Module  Test Pass/Fail Results 
                           1 = Output Class   Test Pass/Fail Results 
                           2 = Output Section Test Pass/Fail Results 
                           3 = Output Single  Test Pass/Fail Results 
                           4 = Output Single  Test Full Details 
 
   These functions assume you want some output text from running the tests,
   otherwise you would just use the TestMods function directly, right?      */  

int GamTest::AskRun(ostream& ostr)
  {
  int keepon=0;
  string yn;
  ostr << "\n\t" << "Proceed Through Failed Tests [y,n]? ";
  cin >> yn;
  if(yn == "y") keepon=1;
  int recurse = 1;
  ostr << "\n\t" << "Failed Module Recursion Levels [0,4]? ";
  cin >> recurse;

  TestRunLevel = 0;                             // Set testing silent
  TestMods(ostr, 0, keepon);                    // Run the GAMMA test
  Header(ostr);                                 // Output test header
  Results(ostr, keepon);                        // Output test results
  if(!TestStatus && recurse>=1)                 // If any module test failed
     ResRec(ostr, keepon, recurse--);		// examine them in detail
  return TestStatus;
  }  

// ____________________________________________________________________________
// J                       GAMMA TEST AUXILIARY FUNCTIONS
// ____________________________________________________________________________

ostream& GamTest::ListTests(ostream& ostr) const
  {
  list<ModTest>::const_iterator i=begin();       // Begin at module tests start
   while(i != end())                             // Loop until end of tests
     {   ostr << (*i).describe() << "\n";
     i++;
     }
  return ostr;
  }

ostream& GamTest::FinishTest(ostream& ostr) const
  {
  ostr << "\n\n";
  if(TestStatus>=1)
    {
    ostr << string(20, ' ')
         << "GAMMA Has Successfully Passed Its Tests";
    ostr << "\n" << string(20, ' ') << string(39, '=');
    ostr << "\n" << string(20, ' ') << string(39, '=');
    }
  else if(TestStatus==0)
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
    {
    ostr << string(19, ' ')
         << "GAMMA Testing Has Not Been Completed!";
    }
  ostr << "\n\n";
  return ostr;
  }


#endif							// GamTest.cc
