/* Gutils.i */
// Swig interface file

%{
#include "Basics/Gutils.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );


void query_parameter(int argc,char* argv[],
                                  int par,const std::string& Q, std::string& V);
void query_parameter(int argc,char* argv[],
                                  int par,const std::string& Q, double&      V);
void query_parameter(int argc,char* argv[],
                                  int par,const std::string& Q, int&         V);

bool ask_set(int argc,char* argv[],int par,const std::string& Q, std::string& V);
bool ask_set(int argc,char* argv[],int par,const std::string& Q, double&      V);
bool ask_set(int argc,char* argv[],int par,const std::string& Q, int&         V);
 
void GAMMAerror(const std::string& hdr, const std::string& msg,             int noret=0);
void GAMMAerror(const std::string& hdr, int eidx,                           int noret=0);
void GAMMAerror(const std::string& hdr, int eidx, const std::string& pname, int noret=0);
volatile void GAMMAfatal();

