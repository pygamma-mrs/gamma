/* SingleTest.cc ************************************************-*-c++-*-
**									**
**                                 G A M M A				**
**									**
**      Single Test                                Implementation	**
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
** Class SingleTest embodies a single test in the GAMMA platform	**
** testing hierarchy. Such tests are direct, examining specific aspects **
** of the platform. Thus a single test will always deal with one        **
** section within one class within one Module.				**
**									**
*************************************************************************/

#ifndef _SingleTest_cc_			// Is this file already included?
#  define _SingleTest_cc_ 1		// If no, then remember it
#  if defined(GAMPRAGMA)		// Using the GNU compiler?
#    pragma implementation		// Then this is the implementation
#  endif

#include <Testing/SingleTest.h>		// Include the interface
#include <Basics/Gutils.h>              // Need GAMMA error messages
#include <fstream>			// Include file streams
#include <string>                       // Include libstdc++ strings
#include <iostream>                     // Include input output streams (cout)

using std::string;			// Using libstdc++ strings
using std::ostream;			// Using libstdc++ output streams

// ----------------------------------------------------------------------------
// --------------------------- PRIVATE FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________ 
// i                      SINGLE TEST ERROR HANDLING
// ____________________________________________________________________________

/*      Input                   ST      : A single test (this)
                                eidx    : Error index
                                noret   : Flag for linefeed (0=linefeed)
        Output                  none    : Error message output
                                          Program execution stopped if fatal */

void SingleTest::STerror(int eidx, int noret) const
  {
  string hdr("Single Test");
  string msg;
  switch (eidx)
    {
    case 11: msg = string("Quaternion Norm Deviates From 1");
             GAMMAerror(hdr,msg,noret); break;                          //(11)
    case 12: msg = string("Set {A,B,C,D} Forms an Invalid Quaternion");
             GAMMAerror(hdr,msg,noret); break;                          //(12)
    default: GAMMAerror(hdr,eidx,noret); break;
    }
  }    

volatile void SingleTest::STfatality(int eidx) const
  {                                                                 
  STerror(eidx, eidx);                          // Output error message
  if(eidx) STerror(0);                          // State that its fatal
  GAMMAfatal();					// Clean exit from program
  }

// ____________________________________________________________________________ 
// ii                       SINGLE TEST DEFAULTS
// ____________________________________________________________________________

/* Each test points to a function that takes as arguments an output stream and
   an integer and returns an integer. When there is no particular test function
   (such as during default contruction), then the test should just point to a
   default function which return true (the integer 1). This is that function.*/

int TrueFct(ostream& ostr, int RL) { return 1; ostr << ""; RL++; }

// ----------------------------------------------------------------------------
// ---------------------------- PUBLIC FUNCTIONS ------------------------------
// ----------------------------------------------------------------------------

// ____________________________________________________________________________
// A                   SINGLE TEST CONSTRUCTION, DESTRUCTION
// ____________________________________________________________________________
 
/* These constructors set up a new single test. 

        Input Arguments                     Result
        ---------------       -------------------------------------------------
               -              An empty test, just good to have, returns true
              ST              Makes an identical copy of the test
*/

SingleTest::SingleTest()
  { 
  TestStatus      = -1;					// Set the test status
  TestName        = string("Default");			// No default name
  TestDescription = string("No Description Specified");	// No default description
  TestTest        = &TrueFct;				// Default test function
  }

SingleTest::SingleTest(const SingleTest& ST)
  { 
  TestStatus      = ST.TestStatus; 
  TestName        = ST.TestName;
  TestDescription = ST.TestDescription;
  TestTest        = ST.TestTest;
  }

SingleTest::~SingleTest () { }
void SingleTest::operator= (const SingleTest &ST)
  { 
  if(this == &ST) return;		// Nothing if self-assignment
  TestStatus      = ST.TestStatus; 
  TestName        = ST.TestName;
  TestDescription = ST.TestDescription;
  TestTest        = ST.TestTest;
  }

// ____________________________________________________________________________
// B                       SINGLE TEST ACCESS FUNCTIONS
// ____________________________________________________________________________

/*         Input                ST : A single test (this)
                  (name)      name : The name of the test     (setting name)
                  (status)  status : Current test status      (setting status)        
                  (descr)    descr : Current test description (setting describ)        
                 (runlev)   runlev : Current test run level   (setting runlev)
		  (test)      test : Current test             (setting test)
           Output (name)      name : The name of the test     (getting name)
                  (status)  status : Current test status      (getting status)
                  (descr)    descr : Current test description (getting describ)
                  (runlev)  runlev : Current test run level   (getting runlev)
		  (test)      test : Result of test test      (getting test)
                              void : If setting interal value  
  
   Note that none of these function will run the test. They only allow users
   to access the current test parameters.                                    */

const string& SingleTest::name()      const  { return TestName; }
      int     SingleTest::status()    const  { return TestStatus; }
const string& SingleTest::describe()  const  { return TestDescription; }
      int     SingleTest::runlevel()  const  { return TestRunLevel; }

  typedef int GTestFct(std::ostream&, int);	// Simplify code (hopefully)
  GTestFct* SingleTest::test()     const  { return TestTest; } ;				// Test function (pointer to)

void SingleTest::name(const      string& N)   { TestName        = N;    }
void SingleTest::status(         int     S)   { TestStatus      = S;    }
void SingleTest::describe(const  string& D)   { TestDescription = D;    }
void SingleTest::runlevel(       int     R)   { TestRunLevel    = R;    }
void SingleTest::test(int (*T)(ostream&,int)) { TestTest        = T;    }

// ____________________________________________________________________________
// C                       SINGLE TEST TESTING FUNCTIONS
// ____________________________________________________________________________

/* Each test has a pointer to a function that is used for actual testing.  This
   test is not invoked unless the user requests that it be run.  The functions
   in this section allow for that. Note that once a test has been run test 
   status flag is set (default = -1 ==> untested) so that subsequent testing
   within a loop can readily be avoided. Also note that each test may produce
   output that is sent to the specified output stream. The amount of output
   often depends upon the current class runlevel (TestRunLevel)              */

int SingleTest::runtest(ostream& ostr, int force)
  {
  if(force)					// If forced, we run test even 
    TestStatus = (*TestTest)(ostr,TestRunLevel);// if it has been run before
  else if(TestStatus == -1)			// If not, we run test only if
    TestStatus = (*TestTest)(ostr,TestRunLevel);// it has not been run before
  return TestStatus;
  }

int SingleTest::RunLevel0(std::ostream& ostr, int force)
  { TestRunLevel=0; return runtest(ostr,force); }

int SingleTest::RunLevel1(std::ostream& ostr, int force)
  { TestRunLevel=1; return runtest(ostr,force); }

int SingleTest::RunLevel2(std::ostream& ostr, int force)
  { TestRunLevel=2; return runtest(ostr,force); }

int SingleTest::RunLevel3(std::ostream& ostr, int force)
  { TestRunLevel=3; return runtest(ostr,force); }

int SingleTest::RunLevel4(std::ostream& ostr, int force)
  { TestRunLevel=4; return runtest(ostr,force); }

// ____________________________________________________________________________
// D                       SINGLE TEST OUTPUT FUNCTIONS
// ____________________________________________________________________________

/*              Input           ST   : A single test (this)
                                ostr : Output ASCII file stream 
                                full : Flag for amount of output 
                Return          void : Test is sent to the output stream    */ 

ostream& SingleTest::print(ostream& ostr) const
  {
  ostr << "\nTest Name: " << TestName 			// Output name
       << "\n  Type        - "  << TestStatus		// Output status
       << "\n  RunLevel    - "  << TestRunLevel		// Output runlevel
       << "\n  Description - "  << TestDescription	// Output description
       << "\n";
  return ostr;
  }

ostream& operator<< (ostream& ostr, const SingleTest& par)
  { return par.print(ostr); }

// ____________________________________________________________________________
// E                  SINGLE TEST STL LIST SUPPORT FUNCTIONS
// ____________________________________________________________________________

/* These functions are MUST be defined if one desired to build C++ STL lists
   and vectors of single test, i.e. list<SingleTest> and vector<SingleTest>.
   They don't make a lot of sense and aren't really used, it is just that
   some compilers complain if they are not present (MSVC++ for one)          */

bool SingleTest::operator==(const SingleTest& ST) const
  {
  if(TestName   != ST.TestName)   return false;
  if(TestStatus != ST.TestStatus) return false;
  return true;
  }

bool SingleTest::operator!=(const SingleTest& ST) const
  { return (!((*this) == ST)); }
 
bool SingleTest::operator<(const SingleTest& ST) const
  { return (TestName < ST.TestName); }

bool SingleTest::operator>(const SingleTest& ST) const
  { return (TestName > ST.TestName); }


// ____________________________________________________________________________
// F                       SINGLE TEST OUTPUT FUNCTIONS
// ____________________________________________________________________________

/* These functions facilitate output of single tests in a nicely formatted
   fashion. Note that when a test result is listed only its status flag is
   utilized. That is, NO TEST IS RUN.  So, any testing should be performed
   prior to calling the result output function. If single test details are
   needed, this is done during the test run with the appropriate output flags
   set, not in these simple output functions.
   
   Header:  This function makes a nice title for the single test below which
   ======   the test specific results can be listed. This is typically used
            when a test has failed and specifics aspects of that failure are
            requested. In this instance the single test is re-run in verbose
            mode after the header is output.
   Result:  This function makes a simple single line output for listing under
   ======   SectTest headers (SectTest::Header).  The format of the output
            should conform to the format of that header.                     */

ostream& SingleTest::Header(ostream& ostr, const string& SectName) const
  {
  int Slen = SectName.length();			// Section name length
  if(!Slen) return Header(ostr);		// Use overload if no section
  string TName = TestName;			// Copy single test name
  int Tlen = TName.length();                    // Length of single test name
  if(!Tlen) {TName=string("Unknown"); Tlen=7;}  // If no test name, default
  int len = 8 + Slen + 8 + Tlen;                // Total length of title line
  int ls = 40-len/2;                            // Total length of output line
  string s1 = string(ls, ' ');                  // Initial spacer
  int tlen = TestDescription.length();          // Length of test description
  ostr << "\n\n"     << s1                      // Output the title line
       << "Section " << SectName
       << " * Test "  << TName << "\n";
  if(TestDescription.length())                  // If the is a description,
     {                                          // include it below the title
     ostr << string(40-(tlen+2)/2, ' ')         //        (centered)
          << "(" << TestDescription << ")";
     if(tlen+2>len)
       {
       len = tlen+2;                            // Adjust length for underline
       s1 = string(40-len/2, ' ');              // Adjust initial spacer
       }
     }
  ostr << "\n" << s1 << string(len, '=')        // Output the title underline 
       << "\n" << s1 << string(len, '=');       // which is two rows of =====
  ostr << "\n\n";
  return ostr;
  }

ostream& SingleTest::Header(ostream& ostr) const
  {
  string TName = TestName;                      // Copy the single test name
  int Tlen = TName.length();                    // Length of single test name
  if(!Tlen) {TName=string("Unknown"); Tlen=7;}  // If no test name, default
  int len = 12 + Tlen;                          // Total length of title line
  int ls = 40-len/2;                            // Total length of output line
  string s1 = string(ls, ' ');                  // Initial spacer
  int tlen = TestDescription.length();          // Length of test description
  ostr << "\n\n" << s1                          // Output the title line
       << "Single Test " << TName << "\n";
  if(TestDescription.length())                  // If the is a description,
     {                                          // include it below the title
     ostr << string(40-(tlen+2)/2, ' ')         //        (centered)
          << "(" << TestDescription << ")";
     if(tlen+2>len)
       {
       len = tlen+2;                            // Adjust length for underline
       s1 = string(40-len/2, ' ');              // Adjust initial spacer
       }
     }
  ostr << "\n" << s1 << string(len, '=')        // Output the title underline 
       << "\n" << s1 << string(len, '=');       // which is two rows of =====
  ostr << "\n\n";
  return ostr;
  }

ostream& SingleTest::Result(ostream& ostr) const
  {
  string CName = TestName;
  int nlen = CName.length();
  if(!nlen) CName = string("Unknown");
  if(nlen > 12) { nlen=12; CName = TestName.substr(0,12); }
  string CDesc = TestDescription;
  int dlen = CDesc.length();
  if(!dlen) CDesc = string("Unknown");
  if(dlen > 54) { dlen=54; CDesc = TestDescription.substr(0,54); }
 
  ostr << " " << CName;
  ostr << string(12-nlen, ' ');
  ostr << "  "
       << CDesc;
  ostr << string(54-dlen, ' ');
  ostr << "    ";
  if(TestStatus>0)      ostr << "PASS";
  else if(!TestStatus)  ostr << "FAIL";
  else                  ostr << "----";
  ostr << "\n";
  return ostr;
  }

#endif							// SingleTest.cc
