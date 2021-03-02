/* SinglePar.i */
// Swig interface file

%module SinglePar

%{
#include "Basics/SinglePar.h"
%}

%include "std_string.i"
%include "std_vector.i"

%feature("autodoc", "1" );

%rename(__eq__)  SinglePar::operator== const;
%rename(__ne__)  SinglePar::operator!= const;
%rename(__lt__)  SinglePar::operator<  const;
%rename(__gt__)  SinglePar::operator>  const;

%rename(__assign__) SinglePar::operator=;


class SinglePar
{

public:
 
SinglePar();				        // Default constructor
SinglePar(const SinglePar& par);	// Self constructor

SinglePar(const std::string& pname, int pdata, const std::string& pstate);
SinglePar(const std::string& pname, double pdata, const std::string& pstate);
SinglePar(const std::string& pname, const std::string& pdata, const std::string& pstate);
SinglePar(const std::string& pname, int ptype, const std::string& pdata, const std::string& pstate);

SinglePar(const std::string& pname);
~SinglePar();
SinglePar& operator=(const SinglePar& par);


void SPerror(int eidx, int noret=0) const;

volatile void SPfatality(int eidx)    const;

void type(int Type);
void data(std::string Data);

int setCoord(std::string& input);

//int setTensor(std::ifstream& inp, std::string& input);
  
const std::string&  name()  const;		// Can't alter name here
const std::string&  data()  const;		// Can't alter data here
const std::string&  state() const;		// Can't alter state here

int type() const;		                // Can't alter type here

void name(const  std::string& Name);	// This sets name to Name
void state(const std::string& State);	// This sets state to State
void GetNS(std::string& name, std::string& state) const;
 
bool parse(std::string& name, int& val, std::string& state, int warn=0) const;

bool parse(std::string& name, double& val, std::string& state, int warn=0) const; 
bool parse(std::string& name, std::string& val, std::string& state, int warn=0) const;

bool parse(std::string& name, double& dx, double& dy, double& dz,
                                               std::string& state, int warn=0) const;

bool parse(std::string& name, int& rank, double& diso, double& delz,
                    double& deta, double& alpha, double& beta, double& gamma,
                                        std::string& state, int warn=0) const;


std::vector<std::string> printStrings() const;

/*
std::ostream& print(std::ostream& ostr) const;
friend std::ostream& operator<< (std::ostream& ostr, const SinglePar& par);
bool write(std::ofstream& ostr, int namelen=10) const;
int read(std::ifstream& inp);
*/

bool operator==(const SinglePar& par) const;
bool operator!=(const SinglePar& par) const;
bool operator< (const SinglePar& par) const;
bool operator> (const SinglePar& par) const;

};


/*
class SinglePar
{

public:
 
MSVCDLC SinglePar();				        // Default constructor
MSVCDLC SinglePar(const SinglePar& par);	// Self constructor

MSVCDLC SinglePar(const std::string& pname, int pdata, const std::string& pstate);
MSVCDLC SinglePar(const std::string& pname, double pdata, const std::string& pstate);
MSVCDLC SinglePar(const std::string& pname, const std::string& pdata, const std::string& pstate);
MSVCDLC SinglePar(const std::string& pname, int ptype, const std::string& pdata, const std::string& pstate);

SinglePar(const std::string& pname);
MSVCDLC ~SinglePar();
MSVCDLL SinglePar& operator=(const SinglePar& par);


void SPerror(int eidx, int noret=0) const;

volatile void SPfatality(int eidx)    const;

void type(int Type);
void data(std::string Data);

int setCoord(std::string& input);

//int setTensor(std::ifstream& inp, std::string& input);
  
MSVCDLL const std::string&  name()  const;		// Can't alter name here
MSVCDLL const std::string&  data()  const;		// Can't alter data here
MSVCDLL const std::string&  state() const;		// Can't alter state here

MSVCDLL int type() const;		                // Can't alter type here

MSVCDLL void name(const  std::string& Name);	// This sets name to Name
MSVCDLL void state(const std::string& State);	// This sets state to State
MSVCDLL void GetNS(std::string& name, std::string& state) const;
 
MSVCDLL bool parse(std::string& name, int& val, std::string& state, int warn=0) const;

MSVCDLL bool parse(std::string& name, double& val, std::string& state, int warn=0) const; 
MSVCDLL bool parse(std::string& name, std::string& val, std::string& state, int warn=0) const;

MSVCDLL bool parse(std::string& name, double& dx, double& dy, double& dz,
                                               std::string& state, int warn=0) const;

MSVCDLL bool parse(std::string& name, int& rank, double& diso, double& delz,
                    double& deta, double& alpha, double& beta, double& gamma,
                                        std::string& state, int warn=0) const;


MSVCDLL std::vector<std::string> printStrings() const;


//MSVCDLL std::ostream& print(std::ostream& ostr) const;
//MSVCDLL friend std::ostream& operator<< (std::ostream& ostr, const SinglePar& par);
//MSVCDLL bool write(std::ofstream& ostr, int namelen=10) const;
//MSVCDLL int read(std::ifstream& inp);

//MSVCDLL bool operator==(const SinglePar& par) const;
//MSVCDLL bool operator!=(const SinglePar& par) const;
//MSVCDLL bool operator< (const SinglePar& par) const;
//MSVCDLL bool operator> (const SinglePar& par) const;

};
*/
