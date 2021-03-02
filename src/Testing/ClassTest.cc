/* ClassTest.cc *************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**      Class Test                                Implementation	**
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
**									**
**  Description								**
**									**
** Class ClassTest embodies a test set for the GAMMA platform testing   **
** hierarchy. Such tests are specific to a particular GAMMA class or    **
** source code file. A class test will always deal with one source file **
** (class or collection of associated functions) within one GAMMA 	**
** Module. That is, any one class (or file of functions) will be named  **
** $GAMMA/src/"module"/"class". Thus this class will be intrinsically   **
** associated with such a file and all its tests apply to it. This is   **
** done through a series of tests associated with each section of the   **
** class (or file of functions).					**
**									**
*************************************************************************/

#ifndef   ClassTest_cc_			// Is this file already included?
#  define ClassTest_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <GamGen.h>			// Include OS specific stuff
#include <Testing/ClassTest.h>		// Include the interface
#include <Testing/SectTest.h>           // Include section tests
#include <Testing/SingleTest.h>         // Include single tests
#include <Basics/Gutils.h>              // Need GAMMA error messages
#include <stdlib.h>
#include <fstream>			// Include file streams
#include <string>                       // Include libstdc++ strings
#include <iostream>                     // Include input output streams (cout)

using std::string;			// Using libstdc++ strings

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i                        CLASS TEST ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   CT      : A class test (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void ClassTest::CTerror(int eidx, int noret) const
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

volatile void ClassTest::CTfatality(int eidx) const
  {                                                                 
  CTerror(eidx, eidx);                          // Output error message
  if(eidx) CTerror(0);                          // State that its fatal
  GAMMAfatal();					// Clean exit from program
  }
 
// ____________________________________________________________________________
// ii                  CLASS TEST INITIALIZATION FUNCTIONS
// ____________________________________________________________________________
 
/* This function will initialize the test results vector TestResults. This
   should be done whenever testing is done if either 1.) The status flag is
   -1, indicating no tests have ever been run, or if the current number of
   tests does not match the size of the tests results vector.  The latter
   indicates that additional tests have been added, potentially with a
   different ordering, since the test results were generated previously.    */
 
void ClassTest::SetResults(int force)
  {
  if(TestResults.size() != size())              // If the vector is the wrong
    {                                           // size we make a new on
    TestResults = std::vector<int>(size());          //   Set the vector anew
    TestStatus = -1;                            //   Flag we are untested
    }
  if(force || TestStatus==-1)                   // We restart results if forced
    {                                           // or we have not tested yet
    for(int i=0; i<int(size()); i++)            // Initialize all results to
      TestResults[i] = -1;                      // be -1 (i.e. untested)
    }
  }

// ____________________________________________________________________________
// iii               CLASS TEST SECTION TEST INDEXING FUNCTIONS
// ____________________________________________________________________________

/* These function deal with indicies of Section Test tests. They allow one to    
   obtain the test Pix from either an interger or a string and they perform
   range checking to insure that the requested test exists in the section.   */
 
bool ClassTest::CheckIndex(int k, int warn) const
  {
  if(k<0 || k>=(int)size())                     // Insure test is in range
    {
    if(warn)
       {
       CTerror(11, 1);                          // Accessed element out of range
       if(warn>1) CTerror(16);                  // Stopping section testing
       } 
    return false;
    }
  return true;
  }  
 
std::list<SectTest>::iterator ClassTest::GetPixNC(int k)
  {
  std::list<SectTest>::iterator i=begin();           // Begin at section tests start
  for(int j=0; j<k && i!=end(); j++) i++;       // Move to requested test
  return i;                                     // Output Pix into list
  }
 
std::list<SectTest>::const_iterator ClassTest::GetPix(int k) const
  {
  std::list<SectTest>::const_iterator i=begin();     // Begin at section tests start
  for(int j=0; j<k && i!=end(); j++) i++;       // Move to requested test
  return i;                                     // Output Pix into list
  }
 
std::list<SectTest>::const_iterator ClassTest::GetPix(const string& N) const
  {
  std::list<SectTest>::const_iterator i=begin();     // Begin at section tests start
  while(i != end())                             // and loop over all of its
    {                                           // section tests
    if(N == (*i).name()) return i;              //   Return i is name match
    i++;                                        //   Else go to the next Pix
    }
  return i;                                     // Output Pix into list
  }  

int ClassTest::GetIndex(const string& N) const
  {
  int k=0;                                      // String index
  std::list<SectTest>::const_iterator i=begin();     // Begin at section tests start
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
// A                   CLASS TEST CONSTRUCTION, DESTRUCTION
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

ClassTest::ClassTest() : stdlistCT()
  {
  TestStatus      = -1;                         // Flag as untested
  TestName        = string("");                 // Set no name
  TestDescription = string("");                 // Set no description
  TestRunLevel    = 0;                          // Set run level 0
  TestType        = true;			// Set type as class
  }
 
ClassTest::ClassTest(const ClassTest& CT) : stdlistCT(CT)
  {
  TestStatus      = CT.TestStatus;              // Copy class test status
  TestName        = CT.TestName;                // Copy class test name
  TestDescription = CT.TestDescription;         // Copy class test descript
  TestRunLevel    = CT.TestRunLevel;            // Copy class test run level
  TestResults     = CT.TestResults;             // Copy class test results
  TestType        = CT.TestType;		// Copy class test type
  }
 
ClassTest& ClassTest::operator= (const ClassTest& CT)
  {
  if(this == &CT) return *this;                 // Don't copy if same one
  stdlistCT::operator=(CT);                     // Copy Section Test list
  TestStatus      = CT.TestStatus;              // Copy class test status
  TestName        = CT.TestName;                // Copy class test name
  TestDescription = CT.TestDescription;         // Copy class test descript
  TestRunLevel    = CT.TestRunLevel;            // Copy class test run level
  TestResults     = CT.TestResults;             // Copy class test results
  TestType        = CT.TestType;		// Copy class test type
  return *this;
  }

// ____________________________________________________________________________
// B                     CLASS TEST ITERATORS & MEMBER ACCESS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't
   have to keep looking them up all the time.

std::list<SectTest>::iterator ClassTest::begin()      Pointer to 1st element
std::list<SectTest>::iterator ClassTest::end()        Pointer to last element
SectTest                 ClassTest::front()      First element
SectTest                 ClassTest::back()       Last element                */

// ____________________________________________________________________________
// C                     CLASS TEST LIST & QUEUE OPERATIONS
// ____________________________________________________________________________ 
 
/* These functions are inherited from the ANSI standard "list" class in the C++ 
   library libstdc++.  I am listing them here so that I & other users don't  
   have to keep looking them up all the time. 
 
ClassTest::push_back(const  SectTest& ST)        Add ST to list end
ClassTest::pop_back()                            Remove ST at list end
ClassTest::push_front(const SectTest& ST)        Add ST to list start
ClassTest::pop_front(const  SectTest& ST)        Remove ST at list start
 
ClassTest::insert(iterator p, SectTest& ST)      Add ST before p
ClassTest::erase(iterator p)                     Remove ST at p
ClassTest::clear()                               Remove all list entries     */


// ____________________________________________________________________________ 
// D                    CLASS TEST ADDITIONAL QUEUE OPERATIONS
// ____________________________________________________________________________ 
 
/* These functions are inherited from the ANSI standard "list" class in the C++ 
   library libstdc++.  I am listing them here so that I & other users don't  
   have to keep looking them up all the time. 
 
int  ClassTest::size()                            Number of entries
bool ClassTest::empty()                           TRUE if pset empty
bool ClassTest::operator==(ClassTest CT)          TRUE if psets equal
bool ClassTest::operator!=(ClassTest CT)          TRUE if psets not equal    */ 

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

int ClassTest::contains(const string& testname) const
  { return (seek(testname) != end()); }
int ClassTest::contains(const SectTest& ST) const
  { return contains(ST.name()); }

std::list<SectTest>::const_iterator ClassTest::seek(const string& tname) const
  {
  SectTest ST;                                  // Working section test
  std::list<SectTest>::const_iterator item=begin();  // Begin at input CT start
  while(item != end())                          // and loop over all of its
    {                                           // section tests
    ST = *item;                                 //   This is a section test
    if(tname == ST.name()) return item;         //   Return item is name match
    item++;                                     //   Go to the next Pix
    }
  return item;
  }

std::list<SectTest>::const_iterator ClassTest::seek(const SectTest& ST) const
  { return seek(ST.name()); }

// ____________________________________________________________________________
// F                       CLASS TEST ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions allow access to simple details of individual Section Tests
   or all tests in one lump sum. They will perform NO TESTING (we want code
   a bit more sophisticated for that - see a later section).  These functions
   facilitate external programs in generating nicely formatted output with the
   Class Test and Section Test results.  The lower case names parallel those
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

 
string ClassTest::GetName() const
  { return TestName; }
 
string ClassTest::GetName(int k) const
  {
  if(!CheckIndex(k,1))                          // Insure test is in range
    { CTerror(12); return string(""); }         // If not, go on w/ no range
  std::list<SectTest>::const_iterator i=begin();     // Begin at section tests start
  for(int j=0; j<k; j++) i++;                   // Move to requested test
  return (*i).name();                           // Output test name
  }
 
std::vector<string> ClassTest::GetNames() const
  {
  std::vector<string> Names;
  std::list<SectTest>::const_iterator i=begin();     // Begin at section tests start
  while(i != end())                             // Loop until end of tests
    { Names.push_back((*i).name()); i++; }      // while putting names in vector
  return Names;
  } 

    
string ClassTest::GetDescription() const
  { return TestDescription; }
    
string ClassTest::GetDescription(int k) const
  {
  if(!CheckIndex(k,1))                          // Insure test is in range
    { CTerror(13); return string(""); }         // If not, go on w/ no descript
  std::list<SectTest>::const_iterator i=begin();     // Begin at section tests start
  for(int j=0; j<k; j++) i++;                   // Move to requested test
  return (*i).describe();                       // Output test description
  }
 
std::vector<string> ClassTest::GetDescriptions() const
  {
  std::vector<string> Descripts;
  std::list<SectTest>::const_iterator i=begin();     // Begin at section tests start
  while(i != end())                             // Loop until end of tests
    { Descripts.push_back((*i).describe()); i++; }
  return Descripts;
  }
 
int ClassTest::GetStatus() const
  { return TestStatus; }
 
int ClassTest::GetStatus(int k) const
  {
  if(!CheckIndex(k,1))                          // Insure test is in range
    { CTerror(14); return -1; }                 // If not, go on w/ status -1
  std::list<SectTest>::const_iterator i=begin();     // Begin at section tests start
  for(int j=0; j<k; j++) i++;                   // Move to requested test
  return (*i).status();                         // Output test status
  }
 
std::vector<int> ClassTest::GetStatuses() const
  {
  std::vector<int> Stati;
  std::list<SectTest>::const_iterator i=begin();     // Begin at section tests start
  while(i != end())                             // Loop until end of tests
    { Stati.push_back((*i).status()); i++; }
  return Stati;
  }

int ClassTest::GetRunLevel() const
  { return TestRunLevel; }
 
int ClassTest::GetRunLevel(int k) const
  {
  if(!CheckIndex(k,1))                          // Insure test is in range
    { CTerror(15); return 0; }                  // If not, go on w/ level 0
  std::list<SectTest>::const_iterator i=begin();     // Begin at section tests start
  for(int j=0; j<k; j++) i++;                   // Move to requested test
  return (*i).runlevel();                       // Output test run level
  }
 
std::vector<int> ClassTest::GetRunLevels() const
  {
  std::vector<int> RLevs;
  std::list<SectTest>::const_iterator i=begin();     // Begin at section tests start
  while(i != end())                             // Loop until end of tests
    { RLevs.push_back((*i).runlevel()); i++; }
  return RLevs;
  } 
    
std::vector<int> ClassTest::GetResults() { SetResults(); return TestResults; }
 
 
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
 
const string& ClassTest::name()     const  { return TestName; }
      int     ClassTest::status()   const  { return TestStatus; }
const string& ClassTest::describe() const  { return TestDescription; }
      int     ClassTest::runlevel() const  { return TestRunLevel; }
      bool    ClassTest::type()     const  { return TestType; }
 
void ClassTest::name(const     string& N) { TestName        = N; }
void ClassTest::status(        int     S) { TestStatus      = S; }
void ClassTest::describe(const string& D) { TestDescription = D; }
void ClassTest::runlevel(      int     R) { TestRunLevel    = R; }
void ClassTest::type(           bool   T) { TestType        = T; }

// ____________________________________________________________________________
// G                        CLASS TEST TESTING FUNCTIONS
// ____________________________________________________________________________

/* These functions will invoke any or all of the Section Tests that are listed
   as part of this Section Test. Results of the Section Tests will be stored
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
//                           Tests All Sections
//-----------------------------------------------------------------------------
 
int ClassTest::TestSects(std::ostream& ostr, int frc, int keepon)
  {
  SetResults(frc);                              // Initialize results(if want)
  std::list<SectTest>::iterator i=begin();           // Begin at section tests start
  int tidx = 0;                                 // This is plain test index
  TestStatus = 1;                               // Assume we test OK
  int goodtest = 1;                             // A flag if we fail
  while(i!=end() && goodtest)                   // Loop until end of tests
    {
    (*i).runlevel(TestRunLevel);                // Set test run level
    goodtest=((*i).TestSingles(ostr,frc,keepon));//  Run section test
    TestResults[tidx++] = goodtest;             //   Store test result
    TestStatus *= abs(goodtest);                //   Sum result for section
    i++;                                        //   Increment to next test
    if(keepon) goodtest=1;                      //   Continue on failure
    }                                           ///     (if desired)
  return TestStatus;
  }
 
//-----------------------------------------------------------------------------
//                 Test Specific Sections (Via Index Or Name)
//-----------------------------------------------------------------------------
 
int ClassTest::TestSect(std::ostream& ostr, int tidx, int frc, int keepon)
  {
  if(!CheckIndex(tidx, 1))                      // Insure test tidx exists
    { CTfatality(17); }                         // Cannot run requested test
  std::list<SectTest>::iterator item;                // A Pix into list of tests
  item = GetPixNC(tidx);                        // Get Pix for requested test
  return ((*item).TestSingles(ostr,frc,keepon));// Run test, return result
  }

int ClassTest::TestSect(std::ostream& ostr, const string& tname, int frc, int kepon)
  {
  int tidx = GetIndex(tname);                   // Get text index
  if(!CheckIndex(tidx, 1))                      // Insure test tname exists
    { CTerror(18); CTfatality(17); }            // Cannot run requested test
  std::list<SectTest>::iterator item;                // A Pix into list of tests
  item = GetPixNC(tidx);                        // Get Pix for requested test
  return ((*item).TestSingles(ostr,frc,kepon));	// Run test, return result
  }
 
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
            output results detailed if the function recursion level
            (nlevels) is set greater than 1.                                */


std::ostream& ClassTest::Header(std::ostream& ostr, const string& ModName)
  {
  int Mlen = ModName.length();			// Module name
  if(!Mlen) return Header(ostr);		// If no module, use overload
  string CName = TestName;			// Get class name
  int Clen = CName.length();			// Get name length
  if(!Clen) CName = string("Unknown");		// If no name, use default
  int len = 21 + Mlen + 2 + Clen;		// Length of title
  (TestType)?len+=6:len+=10;			// Adjust for test type
  int ls = 40-len/2;
  string s1 = string(ls, ' ');
  ostr << "\n\n" << s1				// Output the title line
       << "Testing GAMMA Module "<< ModName
       << ", " << CName;
  if(TestType) ostr << " Class\n";		// Adjust title for test type
  else         ostr << " Functions\n";
  int tlen = TestDescription.length();          // Length of test description
  if(TestDescription.length())                  // If the is a description,
     {                                          // include it below the title
     ostr << string(39-tlen/2, ' ')	//        (centered)
          << "(" << TestDescription << ")";
     len = gmax(len,tlen+2);                     // Adjust length for underline
     s1 = string(40-len/2, ' ');                // Adjust initial spacer
     }
  ostr << "\n" << s1 << string(len, '=')        // Output the title underline
       << "\n" << s1 << string(len, '=');       // which is two rows of =====
  ostr << "\n\n";
  ostr << " Test Section                       Description                        Outcome\n";
  ostr << " ============  ======================================================  =======\n";
  return ostr;
  }

std::ostream& ClassTest::Header(std::ostream& ostr) const
  {
  string CName = TestName;			// Get class name
  int Clen = CName.length();			// Get name length
  if(!Clen) CName = string("Unknown");		// If no name, use default
  int len = 14 + Clen;				// Length of title
  (TestType)?len+=6:len+=10;			// Adjust for test type
  int ls = 40-len/2;                            // Total length of output line
  string s1 = string(ls, ' ');                  // Initial spacer
  ostr << "\n\n" << s1                          // Output the title line
       << "Testing GAMMA "<< CName;		// Adjusted for test type
  if(TestType) ostr << " Class\n";
  else         ostr << " Functions\n";
  int tlen = TestDescription.length();          // Length of test description
  if(tlen)                                      // If there is a description,
     {                                          // include it below the title
     ostr << string(39-tlen/2, ' ')	//        (centered)
          << "(" << TestDescription << ")";
     len = gmax(len,tlen+2);                     // Adjust length for underline
     s1 = string(40-len/2, ' ');                // Adjust initial spacer
     }
  ostr << "\n" << s1 << string(len, '=')        // Output the title underline
       << "\n" << s1 << string(len, '=');       // which is two rows of =====
  ostr << "\n\n";
  ostr << " Test Section                       Description                        Outcome\n";
  ostr << " ============  ======================================================  =======\n";
  return ostr;
  }

   
std::ostream& ClassTest::Result(std::ostream& ostr) const
  {                                                                           
  string CName = TestName;
  if(!CName.length()) CName = string("Unknown");
  string CDesc = TestDescription;
  if(!CDesc.length()) CDesc = string("Unknown");

  ostr << " " << CName;
  int len = CName.length();
  if(len<12) ostr << string(12-len, ' ');
  ostr << "  "
       << CDesc;
  len = CDesc.length();
  if(len<54) ostr << string(54-len, ' ');
  ostr << "   ";
  if(TestStatus>0)      ostr << " PASS";
  else if(!TestStatus)  ostr << " FAIL";
  else                  ostr << "NOTEST";
  ostr << "\n";
  return ostr;
  }
 
std::ostream& ClassTest::Results(std::ostream& ostr, int keepon) const
  {
  std::list<SectTest>::const_iterator item;          // Iterator in Section test list
  bool goon = true;				// Flag to continue output
  item = begin();                               // Begin at section tests start
  while(item != end() && goon)                  // Loop until end of tests &
    {                                           // output section test results
    (*item).Result(ostr);                       //   Output section result
    if(!(*item).status() && !keepon) goon=false;//   Quit if fail & not go on
    item++;                                     //   Move to next section
    }
  return ostr;
  }

std::ostream& ClassTest::ResRec(std::ostream& ostr, int keepon, int nlevels)
  {
  std::list<SectTest>::iterator stpix;              // Pix in section test list
  stpix = begin();                             // Begin at 1st section test
  int ststat;                                  // Section test status
  while(stpix != end())                        // Loop until end of tests
    { 
    ststat = (*stpix).status();                //   Section test result
    if(!ststat)                                //   If sect test failed
      {                                        //   then we will look at
      (*stpix).Header(ostr, TestName);         //   the section's single
      (*stpix).Results(ostr, keepon);          //   test results
      if(nlevels>1)                            //   If desired we can also
        (*stpix).ResRec(ostr,keepon,nlevels--);//   results spit out class 
      }                                        //   section results
    stpix++;
    if(!keepon  && !ststat) stpix = end();     // Bail on fail if !keepon
    } 
  return ostr;
  }
 
// ____________________________________________________________________________
// I                       CLASS TEST INTERACTIVE FUNCTIONS
// ____________________________________________________________________________

/* These functions facilitate the running a class test. Recall that each class
   test runs through a series of section tests, ultimately the test returning
   TRUE/FALSE on completion. However there are several additional facets of
   class testing.

   2.) Continue Through Failed Tests:   If a section test fails should testing
                                        stop or continue on to the section next
   3.) Recurse Through Single Tests:    If a section test fails should testing
                                        recurse into single tests for that
                                        section, displaying specific failures?
  
      Recursion Levels:   0 = Output Section Test Pass/Fail Results
                          1 = Output Single  Test Pass/Fail Results 
                          2 = Output Single  Test Full Details 

   These functions assume you want some output text from running the tests,
   otherwise you would just use the TestLevel function directly, right?      */

int ClassTest::AskRun(std::ostream& ostr)
  {
  int keepon=0;
  string yn;
  std::cout << "\n\t" << "Proceed Through Failed Tests [y,n]? ";
  std::cin >> yn;
  if(yn == "y") keepon=1;
  int recurse = 1;
  std::cout << "\n\t" << "Failed Section Recursion Levels [0,2]? ";
  std::cin >> recurse;

  TestRunLevel = 0;				// Set testing silent
  TestSects(ostr, keepon, recurse);		// Run the section tests
  Header(ostr);  				// Output class header
  Results(ostr);				// Output class results
  if(!TestStatus && recurse>=1)                 // If any section test failed
    ResRec(ostr, keepon, recurse--);            // examine it in detail
  return TestStatus;
  }

// ____________________________________________________________________________
// J                       CLASS TEST AUXILIARY FUNCTIONS
// ____________________________________________________________________________

std::ostream& ClassTest::ListTests(std::ostream& ostr) const
  {
  std::list<SectTest>::const_iterator i=begin();     // Begin at single tests start
   while(i != end())                            // Loop until end of tests
     {   ostr << (*i).describe() << "\n";
     i++;
     }
  return ostr;
  }

#endif							// ClassTest.cc
