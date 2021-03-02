// This file contains interfaces to a set of gamma tests
#include "gamma.h"
#include <string>

#ifndef _TEST_SUITE_H_2009_11_27
#define _TEST_SUITE_H_2009_11_27

class GammaTest

{

public:

static int iso_test();

static int fid_test(std::string & sysfile);

static int spinecho_test(std::string & sysfile);

static int spinecho_realpulse_test(std::string & sysfile, std::string & pulse180file);

static int press_realpulses_test(std::string & sysfile, std::string & pulse180file);

static int fid_exchange_test(std::string & sysfile);

};

#endif // TEST_SUITE_H_2009_11_27
