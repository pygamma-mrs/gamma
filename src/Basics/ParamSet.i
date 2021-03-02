// ParameterSet.i

%{
#include "Basics/ParamSet.h"
%}

%include "std_string.i"
%include "std_vector.i"
%include "std_list.i"
%template(stdlistSP) std::list<SinglePar>;


const std::list<SinglePar> GamSParInit;
const std::vector<int> GamIntVecInit;

typedef std::list<SinglePar> stdlistSP;		// Using typedef on STL list

class ParameterSet : public stdlistSP
{

//public:


int contains(const std::string& pname) const;
int contains(const SinglePar& par)     const;

stdlistSP::const_iterator seek(const std::string& pname) const;
stdlistSP::const_iterator seek(const SinglePar&   par)   const;
 
ParameterSet strip(int indx) const;
int countpar(const std::string& pnamein, int idx0=0);

std::vector<std::string> printStrings() const;

//std::ostream& print(std::ostream& out) const;
//friend std::ostream& operator<< (std::ostream& out, const ParameterSet& pset);

bool write(const std::string& fileout, int warn=2) const;
//bool write(std::ofstream& ofstr,       int warn=2) const;

bool read(const std::string& filein, int fflag=0);
//bool read(std::ifstream& inp,        int fflag=0);

std::string ask_read(int argc, char* argv[], int argn);

bool getParameter(const std::string& name, std::string& value) const;
bool getParameter(const std::string& name, int& value) const;
bool getParameter(const std::string& name, double& value) const;

bool getString(const std::string& name, std::string& value) const;
bool getInt(const    std::string& name, int& value) const;
bool getDouble(const std::string& name, double& value) const;

};

