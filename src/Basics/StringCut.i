/* StringCut.i */
// Swig interface file

%{
#include "Basics/StringCut.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

std::string cutWhite(std::string&     Sinp);
std::string cutString(std::string&    Sinp,                       bool xwhite=true);
std::string cutParBlks(std::string&   Sinp); 
std::string cutBlksXBlks(std::string& Sinp, const std::string& X, bool xwhite=true);
std::string cutDouble(std::string&    Sinp,                       bool xwhite=true);
std::string cutInt(std::string&       Sinp,                       bool xwhite=true);

std::string Gitoa(int i);
std::string Gdec(int i);
std::string Gdec2(long li);
std::string Gdec(const std::string& fmt, int i);
std::string Gdec(int i, int digs);

std::string Gform(const std::string& fmt, double d);
std::string Gform(const std::string& fmt, int i);

//bool Greadline(std::ifstream& s, std::string& Sout, char terminator='\n');

std::string   CenterString(const std::string& str, int width=80);

//std::ostream& CenterString(std::ostream& ostr, const std::string& str, int width=80);

