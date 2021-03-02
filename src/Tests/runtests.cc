// This file contains a set of gamma tests
#include "gamma.h"
#include "testsuite.h"

using namespace std;


int main(int argc, char* argv[])
{
  int out;
  string testname;


  if(argc < 2)
  {
    testname = "all";
  }
  else
  {
	testname = argv[1];
  }

  string sysfile ="gsh_test.sys";
  string pulse180file = "bjs180_1.txt";

  if (testname == "iso")
    out = GammaTest::iso_test();
  else if (testname == "fid")
    out = GammaTest::fid_test(sysfile);
  else if (testname == "spinecho")
    out = GammaTest::spinecho_test(sysfile);
  else if (testname == "spinecho_realpulse")
    out = GammaTest::spinecho_realpulse_test(sysfile, pulse180file);
  else if (testname == "press_realpulse")
    out = GammaTest::press_realpulses_test(sysfile, pulse180file); 
  else if (testname == "fid_exchange")
    out = GammaTest::fid_exchange_test(sysfile);
  else if (testname == "all")
  {
    out = GammaTest::iso_test();
    out = GammaTest::fid_test(sysfile);
    out = GammaTest::spinecho_test(sysfile);
    out = GammaTest::spinecho_realpulse_test(sysfile, pulse180file);
    out = GammaTest::press_realpulses_test(sysfile, pulse180file);
    out = GammaTest::fid_exchange_test(sysfile);
  }
  else 
  {
    cout << "Sorry, I don't know about a test named " << testname << endl;
  }

  return 0;
}
