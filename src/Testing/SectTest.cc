/* SectTest.cc **************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**      Section Test                                Implementation	**
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
** Class SectTest embodies a test set for the GAMMA platform testing    **
** hierarchy. Such tests are specific to a particular GAMMA section in  **
** a source code file. A section test will always deal with one section **
** of a source file (class or collection of associated functions) that  **
** lies within a single GAMMA Module. A section test may contain        **
** multiple single tests.                                               **
**                                                                      **
*************************************************************************/

#ifndef   SectTest_cc_			// Is this file already included?
#  define SectTest_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <GamGen.h>			// Include OS specific stuff
#include <Testing/SectTest.h>		// Include the interface
#include <Testing/SingleTest.h>		// Include single tests
#include <Basics/Gutils.h>              // Need GAMMA error messages
#include <stdlib.h>
#include <fstream>			// Include file streams
#include <string>                       // Include libstdc++ strings
#include <iostream>                     // Include input output streams (cout)

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

void SectTest::STerror(int eidx, int noret) const
  {
  std::string hdr("Section Test");
  std::string msg;
  switch (eidx)
    {
    case 11: msg = std::string("Accessed Single Test Out Of Range");
             GAMMAerror(hdr,msg,noret); break;                          //(11)
    case 12: msg = std::string("Continuing With No Test Name");
             GAMMAerror(hdr,msg,noret); break;                          //(12)
    case 13: msg = std::string("Continuing With No Test Description");
             GAMMAerror(hdr,msg,noret); break;                          //(13)
    case 14: msg = std::string("Continuing With Status Of -1");
             GAMMAerror(hdr,msg,noret); break;                          //(14)
    case 15: msg = std::string("Continuing With Run Level Of 0");
             GAMMAerror(hdr,msg,noret); break;                          //(15)
    case 16: msg = std::string("Stopping Section Testing");
             GAMMAerror(hdr,msg,noret); break;                          //(16)
    case 17: msg = std::string("Unable To Run Specified Single Test");
             GAMMAerror(hdr,msg,noret); break;                          //(17)
    case 18: msg = std::string("Single Test By Name Does Not Exist");
             GAMMAerror(hdr,msg,noret); break;                          //(18)
    default: GAMMAerror(hdr,eidx,noret); break;
    }
  }    

volatile void SectTest::STfatality(int eidx) const
  {                                                                 
  STerror(eidx, eidx);                          // Output error message
  if(eidx) STerror(0);                          // State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________ 
// ii                 SECTION TEST INITIALIZATION FUNCTIONS
// ____________________________________________________________________________

/* This function will initialize the test results vector TestResults. This 
   should be done whenever testing is done if either 1.) The status flag is
   -1, indicating no tests have ever been run, or if the current number of
   tests does not match the size of the tests results vector.  The latter 
   indicates that additional tests have been added, potentially with a
   different ordering, since the test results were generated previously.    */

void SectTest::SetResults(int force)
  {
  if(TestResults.size() != size())		// If the vector is the wrong
    {						// size we make a new on
    TestResults = std::vector<int>(size());		//   Set the vector anew 
    TestStatus = -1;				//   Flag we are untested
    }
  if(force || TestStatus==-1)			// We restart results if forced
    {						// or we have not tested yet
    for(int i=0; i<int(size()); i++)		// Initialize all results to
      TestResults[i] = -1;			// be -1 (i.e. untested)
    }
  }

// ____________________________________________________________________________
// iii               SECTION TEST SINGLE TEST INDEXING FUNCTIONS
// ____________________________________________________________________________
 
/* These function deal with indicies of Single Test tests. They allow one to 
   obtain the test Pix from either an interger or a string and they perform
   range checking to insure that the requested test exists in the section.   */

bool SectTest::CheckIndex(int k, int warn) const
  {
  if(k<0 || k>=(int)size())			// Insure test is in range
    {
    if(warn)
       {
       STerror(11, 1);				// Accessed element out of range
       if(warn>1) STerror(16);			// Stopping section testing
       }
    return false;
    }
  return true; 
  }

std::list<SingleTest>::iterator SectTest::GetPixNC(int k)
  {
  std::list<SingleTest>::iterator i=begin();         // Begin at single tests start
  for(int j=0; j<k && i!=end(); j++) i++;	// Move to requested test
  return i;					// Output Pix into list
  }

std::list<SingleTest>::const_iterator SectTest::GetPix(int k) const
  {
  std::list<SingleTest>::const_iterator i=begin();   // Begin at single tests start
  for(int j=0; j<k && i!=end(); j++) i++;	// Move to requested test
  return i;					// Output Pix into list
  }

std::list<SingleTest>::const_iterator SectTest::GetPix(const std::string& N) const
  {
  std::list<SingleTest>::const_iterator i=begin();  // Begin at single tests start
  while(i != end())                             // and loop over all of its
    {                                           // single tests
    if(N == (*i).name()) return i;              //   Return i is name match
    i++;                                        //   Else go to the next Pix
    }
  return i;					// Output Pix into list
  }

int SectTest::GetIndex(const std::string& N) const
  {
  int k=0;					// String index
  std::list<SingleTest>::const_iterator i=begin();  	// Begin at single tests start
  while(i != end())                             // and loop over all of its
    {                                           // single tests
    if(N == (*i).name()) return k;              //   Return i is name match
    i++;                                        //   Else go to the next Pix
    k++;					//   and increment integer
    }
  if(i==end()) k=-1;				// Return -1 if not in list
  return k;					// Output integer index
  }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                   SECTION TEST CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 
 
/* Were it not for the additional variables defined in this class, all the
   constructors would just be inherited from the ANSI standard "list" class 
   in the C++ library libstdc++. Since this is not the case, I am listing all
   such constructors (commented out) and will only (re)define the ones that
   I care about.

SectTest::SectTest()                               Empty Section Test
SectTest::SectTest(int N)                          Section Test w/ N Tests
SectTest::SectTest(int N, const SingleTest& ST)    Section Test w/ N ST
SectTest::SectTest(const SectTest& PT)             Section Test copy of PT
SectTest::SectTest assign(N)                       Assign N element
SectTest::~SectTest()                              Destructor                */

SectTest::SectTest() : stdlistST()
  {
  TestStatus      = -1; 			// Flag as untested
  TestName        = std::string(""); 		// Set no name
  TestDescription = std::string(""); 		// Set no description
  TestRunLevel    = 0; 				// Set run level 0
  }

SectTest::SectTest(const SectTest& PT) : stdlistST(PT)
  {
  TestStatus      = PT.TestStatus; 		// Copy section test status
  TestName        = PT.TestName; 		// Copy section test name
  TestDescription = PT.TestDescription; 	// Copy section test descript
  TestRunLevel    = PT.TestRunLevel; 		// Copy section test run level
  TestResults     = PT.TestResults; 		// Copy section test results
  }

SectTest& SectTest::operator= (const SectTest& PT)
  {
  if(this == &PT) return *this;			// Don't copy if same one
  stdlistST::operator=(PT);			// Copy Single Test list
  TestStatus      = PT.TestStatus; 		// Copy section test status
  TestName        = PT.TestName; 		// Copy section test name
  TestDescription = PT.TestDescription; 	// Copy section test descript
  TestRunLevel    = PT.TestRunLevel; 		// Copy section test run level
  TestResults     = PT.TestResults; 		// Copy section test results
  return *this;
  }

// ____________________________________________________________________________
// B                     SECTION TEST ITERATORS & MEMBER ACCESS
// ____________________________________________________________________________

/* These functions are inherited from the ANSI standard "list" class in the C++
   library libstdc++.  I am listing them here so that I & other users don't
   have to keep looking them up all the time.

std::list<SingleTest>::iterator SectTest::begin()      Pointer to 1st element
std::list<SingleTest>::iterator SectTest::end()        Pointer to last element
SingleTest                 SectTest::front()      First element
SingleTest                 SectTest::back()       Last element               */

// ____________________________________________________________________________
// C                     SECTION TEST LIST & QUEUE OPERATIONS
// ____________________________________________________________________________ 
 
/* These functions are inherited from the ANSI standard "list" class in the C++ 
   library libstdc++.  I am listing them here so that I & other users don't  
   have to keep looking them up all the time. 
 
SectTest::push_back(const    SingleTest& ST)     Add ST to list end
SectTest::pop_back()                             Remove ST at list end
SectTest::push_front(const   SingleTest& ST)     Add ST to list start
SectTest::pop_front(const    SingleTest& ST)     Remove ST at list start
SectTest::insert(iterator p, SingleTest& ST)     Add ST before p
SectTest::erase(iterator p)                      Remove ST at p
SectTest::clear()                                Remove all list entries    */ 

// ____________________________________________________________________________ 
// D                    SECTION TEST ADDITIONAL QUEUE OPERATIONS
// ____________________________________________________________________________ 
 
/* These functions are inherited from the ANSI standard "list" class in the C++ 
   library libstdc++.  I am listing them here so that I & other users don't  
   have to keep looking them up all the time. 
 
int  SectTest::size()                            Number of entries
bool SectTest::empty()                           TRUE if pset empty
bool SectTest::operator==(SectTest PT)          TRUE if psets equal
bool SectTest::operator!=(SectTest PT)          TRUE if psets not equal      */ 

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

int SectTest::contains(const std::string& testname) const
  { return (seek(testname) != end()); }
int SectTest::contains(const SingleTest& ST) const
  { return contains(ST.name()); }

std::list<SingleTest>::const_iterator SectTest::seek(const std::string& tname) const
  {
  SingleTest ST;                                // Working single test
  std::list<SingleTest>::const_iterator item=begin();// Begin at input PT start
  while(item != end())                          // and loop over all of its
    {                                           // single tests
    ST = *item;                                 //   This is a single test
    if(tname == ST.name()) return item;         //   Return item is name match
    item++;                                     //   Go to the next Pix
    }
  return item;
  }

std::list<SingleTest>::const_iterator SectTest::seek(const SingleTest& ST) const
  { return seek(ST.name()); }

// ____________________________________________________________________________
// F                    SECTION TEST ACCESS FUNCTIONS
// ____________________________________________________________________________

/* These functions allow access to simple details of individual Single Tests
   or all tests in one lump sum. They will perform NO TESTING (we want code
   a bit more sophisticated for that - see a later section).  These functions
   facilitate external programs in generating nicely formatted output with the
   Section Test and Single Test results.  The lower case names parallel those
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

std::string SectTest::GetName() const
  { return TestName; }

std::string SectTest::GetName(int k) const
  {
  if(!CheckIndex(k,1)) 				// Insure test is in range
    { STerror(12); return std::string(""); } 	// If not, go on w/ no range
  std::list<SingleTest>::const_iterator i=begin();   // Begin at single tests start
  for(int j=0; j<k; j++) i++;			// Move to requested test
  return (*i).name();				// Output test name
  }

std::vector<std::string> SectTest::GetNames() const
  {
  std::vector<std::string> Names;
  std::list<SingleTest>::const_iterator i=begin();   // Begin at single tests start
  while(i != end())				// Loop until end of tests 
    { Names.push_back((*i).name()); i++; }	// while putting names in vector
  return Names;
  }

std::string SectTest::GetDescription() const
  { return TestDescription; }

std::string SectTest::GetDescription(int k) const
  {
  if(!CheckIndex(k,1)) 				// Insure test is in range
    { STerror(13); return std::string(""); } 	// If not, go on w/ no descript
  std::list<SingleTest>::const_iterator i=begin();   // Begin at single tests start
  for(int j=0; j<k; j++) i++;			// Move to requested test
  return (*i).describe();			// Output test description
  }

std::vector<std::string> SectTest::GetDescriptions() const
  {
  std::vector<std::string> Descripts;
  std::list<SingleTest>::const_iterator i=begin();   // Begin at single tests start
  while(i != end())				// Loop until end of tests 
    { Descripts.push_back((*i).describe()); i++; }
  return Descripts;
  }

int SectTest::GetStatus() const
  { return TestStatus; }

int SectTest::GetStatus(int k) const
  {
  if(!CheckIndex(k,1)) 				// Insure test is in range
    { STerror(14); return -1; } 		// If not, go on w/ status -1
  std::list<SingleTest>::const_iterator i=begin();   // Begin at single tests start
  for(int j=0; j<k; j++) i++;			// Move to requested test
  return (*i).status();				// Output test status
  }

std::vector<int> SectTest::GetStatuses() const
  {
  std::vector<int> Stati;
  std::list<SingleTest>::const_iterator i=begin();   // Begin at single tests start
  while(i != end())				// Loop until end of tests 
    { Stati.push_back((*i).status()); i++; }
  return Stati;
  }

int SectTest::GetRunLevel() const
  { return TestRunLevel; }

int SectTest::GetRunLevel(int k) const
  {
  if(!CheckIndex(k,1)) 				// Insure test is in range
    { STerror(15); return 0; } 			// If not, go on w/ level 0
  std::list<SingleTest>::const_iterator i=begin();   // Begin at single tests start
  for(int j=0; j<k; j++) i++;			// Move to requested test
  return (*i).runlevel();			// Output test run level
  }

std::vector<int> SectTest::GetRunLevels() const
  {
  std::vector<int> RLevs;
  std::list<SingleTest>::const_iterator i=begin();   // Begin at single tests start
  while(i != end())				// Loop until end of tests 
    { RLevs.push_back((*i).runlevel()); i++; }
  return RLevs;
  }

std::vector<int> SectTest::GetResults() { SetResults(); return TestResults; }

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

const std::string& SectTest::name()     const  { return TestName; }
      int     SectTest::status()   const  { return TestStatus; }
const std::string& SectTest::describe() const  { return TestDescription; }
      int     SectTest::runlevel() const  { return TestRunLevel; }

void SectTest::name(const     std::string& N) { TestName        = N; }
void SectTest::status(        int     S) { TestStatus      = S; }
void SectTest::describe(const std::string& D) { TestDescription = D; }
void SectTest::runlevel(      int     R) { TestRunLevel    = R; }

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
//                           Run All Single Tests
//-----------------------------------------------------------------------------

int SectTest::TestSingles(std::ostream& ostr, int anew, int keepon)
  {
  SetResults(anew);				// Initialize results(if want)
  std::list<SingleTest>::iterator i=begin(); 	// Begin at single tests start
  int tidx = 0;					// This is plain test index
  TestStatus = 1;				// Assume we test OK
  int goodtest = 1;				// A flag if we fail
  while(i!=end() && goodtest)			// Loop until end of tests 
    {
    (*i).runlevel(TestRunLevel);		// Set test run level
    goodtest = ((*i).runtest(ostr, anew));	//   Run single test
    TestResults[tidx++] = goodtest;		//   Store test result
    TestStatus *= abs(goodtest);		//   Sum result for section 
    i++;					//   Increment to next test
    if(keepon) goodtest=1;			//   Continue on failure
    }						///     (if desired)
  return TestStatus;
  }

//-----------------------------------------------------------------------------
//                  Run Specific Single Test (By Index or Name)
//-----------------------------------------------------------------------------

int SectTest::TestSingle(std::ostream& ostr, int tidx, int anew)
  {
  if(!CheckIndex(tidx, 1)) 			// Insure test tidx exists
    { STfatality(17); }				// Cannot run requested test
  std::list<SingleTest>::iterator item;              // A Pix into list of tests
  item = GetPixNC(tidx);			// Get Pix for requested test
  return ((*item).runtest(ostr, anew));	// Run test, return result
  }

int SectTest::TestSingle(std::ostream& ostr, const std::string& tname, int anew)
  {
  int tidx = GetIndex(tname);			// Get text index
  if(!CheckIndex(tidx, 1)) 			// Insure test tname exists
    { STerror(18); STfatality(17); }		// Cannot run requested test
  std::list<SingleTest>::iterator item; 	        // A Pix into list of tests
  item = GetPixNC(tidx);			// Get Pix for requested test
  return ((*item).runtest(ostr, anew));	// Run test, return result
  }

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

std::ostream& SectTest::Header(std::ostream& ostr, const std::string& ClassName)
  {
  int Clen = ClassName.length();		// Class name length
  if(!Clen) return Header(ostr);		// Use overload if no class
  std::string SName = TestName;			// Copy the section name
  int Slen = SName.length();			// Length of section name
  if(!Slen) {SName=std::string("Unknown"); Slen=7;}	// If no section name, default
  if(Slen>12) Slen=12;				// No more than 12 chars here
  int len = 20 + Clen + 10 + Slen;		// Total length of title line
  int ls = 40-len/2;				// Total length of output line
  std::string s1 = std::string(ls?ls:1, ' ');		// Initial spacer
  ostr << "\n\n" << s1				// Output the title line
       << "Testing GAMMA Class "<< ClassName
       << ", Section " << SName << "\n";
  int tlen = TestDescription.length();		// Length of test description
  if(tlen > 54) tlen=54;			// No more that 54 chars here
  if(TestDescription.length())                  // If the is a description,
     {                                          // include it below the title
     ostr << std::string(39-tlen/2, ' ')             //        (centered)
          << "(" << TestDescription << ")";
     len = gmax(len,tlen+2);                     // Adjust length for underline
     s1 = std::string(40-len/2, ' ');		// Adjust initial spacer
     }
  ostr << "\n" << s1 << std::string(len, '=')        // Output the title underline 
       << "\n" << s1 << std::string(len, '=');       // which is two rows of =====
  ostr << "\n\n";
  ostr << "     Test                           Description                        Outcome\n";
  ostr << " ============  ======================================================  =======\n";
  return ostr;
  }

std::ostream& SectTest::Header(std::ostream& ostr) const
  {
  std::string SName = TestName;			// Copy the section name
  int Slen = SName.length();			// Length of section name
  if(!Slen) {SName=std::string("Unknown"); Slen=7;}  // If no section name, default
  if(Slen>12) Slen=12;				// No more than 12 chars here
  int len = 22 + Slen;				// Total length of title line
  int ls = 40-len/2;                            // Total length of output line
  std::string s1 = std::string(ls?ls:1, ' ');             // Initial spacer
  int tlen = TestDescription.length();		// Length of test description
  if(tlen > 54) tlen=54;			// No more that 54 chars here
  ostr << "\n\n" << s1                          // Output the title line
       << "Testing GAMMA Section " << SName << "\n";
  if(TestDescription.length())			// If the is a description,
     { 						// include it below the title
     ostr << std::string(39-tlen/2, ' ') 		//        (centered)
          << "(" << TestDescription << ")";
     len = gmax(len,tlen+2);			// Adjust length for underline
     s1 = std::string(40-len/2, ' ');		// Adjust initial spacer
     }
  ostr << "\n" << s1 << std::string(len, '=')        // Output the title underline 
       << "\n" << s1 << std::string(len, '=');       // which is two rows of =====
  ostr << "\n\n";
  ostr << "     Test                           Description                        Outcome\n";
  ostr << " ============  ======================================================  =======\n";
  return ostr;
  }


std::ostream& SectTest::Result(std::ostream& ostr) const
  {
  std::string CName = TestName;
  int nlen = CName.length();
  if(!nlen) CName = std::string("Unknown");
  if(nlen > 12) { nlen=12; CName=TestName.substr(0,12); }
  std::string CDesc = TestDescription;
  int dlen = CDesc.length();
  if(!dlen) CDesc = std::string("Unknown");
  if(dlen > 54) { dlen=54; CDesc=TestDescription.substr(0,54); }
 
  ostr << " " << CName;
  ostr << std::string(12-nlen, ' ');
  ostr << "  "
       << CDesc;
  ostr << std::string(54-dlen, ' ');
  ostr << "  ";
  if(TestStatus>0)      ostr << "  PASS";
  else if(!TestStatus)  ostr << "  FAIL";
  else                  ostr << " NOTEST";
  ostr << "\n";
  return ostr;
  }
 
std::ostream& SectTest::Results(std::ostream& ostr, int keepon) const
  {
  std::list<SingleTest>::const_iterator item;        // Iterator in Single test list
  bool goon = true;                             // Flag to continue output
  item = begin();                               // Begin at Single tests start
  while(item != end() && goon)                  // Loop until end of tests &
    {                                           // output single test results
    (*item).Result(ostr);                       //   Output single test result
    if(!(*item).status() && !keepon) goon=false;//   Quit if fail & not go on
    item++;                                     //   Move to next single test
    }
  return ostr;
  }

// sosiz
std::ostream& SectTest::ResRec(std::ostream& ostr, int keepon, int nlevels)
  {
  std::list<SingleTest>::iterator stpix;            // Pix in single test list
  stpix = begin();                             // Begin at 1st single test
  int ststat;                                  // Single test status
  while(stpix != end())                        // Loop until end of tests
    { 
    ststat = (*stpix).status();                //   Single test result
    if(!ststat)                                //   If single test failed
      {                                        //   then we will look at
      (*stpix).Header(ostr, TestName);         //   its problesm in detail
      (*stpix).runlevel(4);
      (*stpix).runtest(ostr,1);
      }                                        //   section results
    stpix++;
    if(!keepon  && !ststat) stpix = end();  //     Bail on fail if !keepon
    } 
  return ostr;
  nlevels++;
  }

// ____________________________________________________________________________
// I                   SECTION TEST INTERACTIVE FUNCTIONS
// ____________________________________________________________________________

/* These functions facilitate the running a section test. Recall that each
   section test runs through a series of single tests, ultimately the test 
   returning TRUE/FALSE on completion. However there are several additional 
   facets of section testing.

   1.) Continue Through Failed Tests:   If a single test fails should testing
                                        stop or continue on with single tests
   2.) Recurse Through Single Tests:    If a single test fails should testing
                                        recurse into that test while displaying
                                        specific test outcome features?

      Recursion Levels:   0 = Output Single Test Pass/Fail Results 
                          1 = Output Single Test Full Details 

   These functions assume you want some output text from running the tests,
   otherwise you would just use the TestLevel function directly, right?       */

int SectTest::AskRun(std::ostream& ostr)
  {
  int keepon=0;
  std::string yn;
  std::cout << "\n\t" << "Proceed Through Failed Tests [y,n]? ";
  std::cin >> yn;
  if(yn == "y") keepon=1;
  int recurse = 1;
  std::cout << "\n\t" << "Recurse Into Failed Single Tests [y,n]? ";
  std::cin >> recurse;

  TestRunLevel = 0;				// Set Testing Silent
  TestSingles(ostr, keepon, recurse);		// Run The Single Tests
  Header(ostr);					// Output Section Header
  Results(ostr);                                // Output Section Results
  if(!TestStatus && recurse>=1)                  // If any single test failed
    ResRec(ostr, keepon, recurse--);            // examine it in detail
  return TestStatus;
  }


// ____________________________________________________________________________
// J                       SECTION TEST AUXILIARY FUNCTIONS
// ____________________________________________________________________________
 
std::ostream& SectTest::ListTests(std::ostream& ostr) const
  {
  std::list<SingleTest>::const_iterator i=begin();   // Begin at single tests start
   while(i != end())				// Loop until end of tests 
     {   ostr << (*i).describe() << "\n";
     i++;
     }   
  return ostr;
  }

#endif							// SectTest.cc
