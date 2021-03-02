/* SpinMap.i */

%{
#include "MultiSys/SpinMap.h"
%}

%include "std_string.i"
%include "std_vector.i"

%feature("autodoc", "1" );

%rename(__assign__) SpinMap::operator=;

class SpinMap
  {
  public:

    int sub1;	
    int sp1;	
    int sub2;	
    int sp2;	

    void SMerror(int eidx,                        int nr=0) const;
    volatile void SMfatal(int eidx)                                  const;
    void SMerror(int eidx, const std::string& pn, int nr=0) const;

    bool getSMStr(const ParameterSet& pset, int idx, int mdx,
                                         std::string& sm, bool warn=true) const;

    bool getSM(const ParameterSet& pset, int idx, int mdx,
         int& comp1, int& spin1, int& comp2, int& spin2, bool warn=true) const;

    bool setSM(const ParameterSet& pset,int idx,int mdx,bool warn=true);

    bool Check(bool warn=true) const;
    bool Check(int c1, int s1, int c2, int s2, bool warn=true) const;

    SpinMap();
    SpinMap(int c1, int s1, int c2, int s2);
    SpinMap(const SpinMap& SM);
    SpinMap(const std::string& SM);   
    SpinMap& operator = (const SpinMap& SM);
    ~SpinMap() {};
    int Sub1()  const;			// Component 1 index
    int Sub2()  const;			// Component 2 index
    int Spin1() const;			// Spin 1 index
    int Spin2() const;			// Spin 2 index

    bool read(const std::string&  filename,int idx,int mdx,int warn=2);
    bool read(const ParameterSet& pset,    int idx,int mdx,int warn=2);
    void          print() const;

};


