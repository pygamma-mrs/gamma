/* GenOpRep.h */

%{
#include "HSLib/GenOpRep.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

%rename(__eq__)  genoprep::operator== const;
%rename(__ne__)  genoprep::operator!= const;
%rename(__lt__)  genoprep::operator<  const;
%rename(__gt__)  genoprep::operator>  const;

%rename(__assign__) genoprep::operator=;


class genoprep
{

public:
 
matrix RepMx;			// Op matrix representation
basis  RepBs;			// Associated Op basis (array)
int    RepPty;			// Representation priority

static bool   BSPrnt;			// Flag to print basis array
 
void OpReperror(int eidx, int noret=0) const;

volatile void OpRepfatal(int eidx) const;

genoprep();
genoprep(const genoprep& OpRep);
genoprep(const matrix& mx, const basis& bs, int pty);
~genoprep();

genoprep& operator= (const genoprep& OpRep); 

bool OpRepCheck(int warn=2) const;

//       std::ostream& print(std::ostream& ostr, int full=0) const;
//friend std::ostream& operator<< (std::ostream& ostr,const genoprep &OpRep);

void write(const std::string& fn) const;

//std::ofstream& write(std::ofstream& fp)     const;

void read(const std::string& fn);

//std::ifstream& read(std::ifstream& fp);

bool operator==(const genoprep& OpRep) const;
bool operator!=(const genoprep& OpRep) const;
bool operator<(const  genoprep& OpRep) const;
bool operator>(const  genoprep& OpRep) const;

};
