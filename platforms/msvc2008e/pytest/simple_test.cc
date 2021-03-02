#include <iostream>

// **** READ THIS ***

// You need to manually INSTALL pygamma in a suitable location 
// (e.g. C:\Python25\Lib\site-packages) on your machine
// before running this project.

// NOTE: simple_test.cc is just a place holder, and a file to "touch"
// to get the project to rebuild.

// The real work for this "pytests" project is done as a post-build
// event.  To view this , right click on the pytests project and select 
// properties, then go to:
// Configuration Properties-->Build Events-->Post Build Events
// and click on the property, "Command Line").

// The code for run_tests.py is in "...gamma/trunk/src/Tests". 
// Because of the directory where this code is run, it will open the file:
// gamma_pytests.suite and run the python code in the direcory
// .../gamma/test/src/pyTests

// Happy testing.

int main(int argc, char *argv[])
{
	  std::cout << "Greetings Earthlings!";
    return 0;
}