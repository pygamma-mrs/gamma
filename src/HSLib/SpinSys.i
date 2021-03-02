/* SpinSys.i */
// Swig interface file.

%{
#include "HSLib/SpinSys.h"
%}

%include "std_vector.i"
%include "std_string.i"

%rename(__assign__) spin_sys::operator=;


class spin_sys
{

public:

    int check_spin(int  spin, int die=1) const;
    int check_spins(int spin1, int spin2, int die=1) const;

    spin_sys();
    spin_sys(int  spins);
    spin_sys(const  spin_sys& sys);

    ~spin_sys();


   virtual spin_sys& operator=(const spin_sys& sys);

	// int??  shouldn't this be changed to bool?
    //int operator==(const spin_sys& sys) const;
    //int operator!=(const spin_sys &sys) const;

    int spins()      const;
    int spinpairs()  const;
    int HS()         const;
    int HS(int spin) const;


    virtual void isotope(int, const std::string&);

    //virtual void  isotope(int, const Isotope&);
    //virtual const  Isotope&  isotope(int)  const;


    double weight(int)                    const;
    std::string symbol(int spin)          const;

    double      qn(int spin)              const;
    double      qn()                      const;

    std::string element(int spin)         const;
    std::string momentum(int spin)        const;
    std::string momentum()                const;

    double  gamma(int spin)               const;
    double  gamma(const std::string& iso) const;

    //const  std::vector<Isotope>& IsoVec() const;

    std::vector<int> HSvect() const;
    row_vector  qState(int state) const;

    //matrix qStates() const;

    double qnState(int state) const;

    col_vector qnStates() const;

    row_vector qnDist() const;

    row_vector CoherDist() const;

    bool   homonuclear()                  const;
	bool   heteronuclear()                const;
	bool   electron(int i)                const;
	bool   nucleon(int i)                 const;
	bool   spinhalf()                     const;
	int    electrons()                    const;
	int    nucleons()                     const;
	bool   nepair(int i,  int j)          const;
	bool   enpair(int i,  int j)          const;
	bool   eepair(int i,  int j)          const;
	bool   nnpair(int i,  int j)          const;
	int    pairidx(int i, int j)          const;
	int    isotopes()                     const;


	std::string isotopes(int idx)         const;
	bool  isotopes(const std::string& I)  const;


	void    SetFlag(int spin, bool TF);
	void    SetFlags(bool TF);
	void    SetFlags(const std::string& isoin, bool TF);
	//void    SetFlags(const Isotope& Iso, bool TF);
	bool    GetFlag(int i) const;

    /*
	flagvec GetFlags() const;
	flagvec GetFlags(bool TF) const;
	flagvec GetFlags(int spin, bool TF, bool DefTF=0) const;
	flagvec GetFlags(const std::string& isoin,bool TF,bool DTF=0) const;
	flagvec GetFlags(const Isotope& isoin,bool TF,bool DTF=0) const;
    */

	void  name(const std::string& name);
	const std::string& name(int i=-1) const;

	void warnings(int warnf);
	int  warnings() const;

	std::string IsoDefault();

	void IsoDefault(const std::string& DI);


	//operator ParameterSet( ) const;
	//friend  void operator+=(ParameterSet& pset, const spin_sys &ss);


	virtual void PSetAdd(ParameterSet& pset, int idx=-1) const;

	int getSpins(const ParameterSet& pset, int warn=0) const;

	void setName(const ParameterSet& pset);

	void setIs(const ParameterSet& pset);

	//void operator=(const ParameterSet& pset);

    //virtual int write(const std::string &filename, int ix=-1, int wn=2) const;
    //virtual int write(std::ofstream& ofstr, int idx=-1, int warn=2) const;

    //virtual int read(const std::string& filename, int idx=-1, int warn=2);
    //virtual int read(const ParameterSet& pset,    int idx=-1, int warn=2);


    virtual std::string ask_read(int argc, char* argv[], int argn);
    virtual std::string ask_read(int argc, char* argv[], int argn,
                                             const std::string& def);

    basis get_basis() const;

    //matrix BasisMap1() const;
    //matrix BasisMap2() const;
    //matrix TransitionMap2() const;

    //virtual std::ostream& print(std::ostream& out, bool hdr=true) const;

    //friend  std::ostream& operator<<(std::ostream& out, const spin_sys& sys);

    virtual std::vector<std::string> printstrings() const;
    virtual std::vector<std::string> SYSStrings(int w1=10,int w2=5,int w3=1) const;

    std::vector<std::string> SIStrings(int  colwd=10) const;
    std::vector<std::string> SYMStrings(int colwd=10) const;
    std::vector<std::string> SAMStrings(int colwd=10) const;
};

%extend spin_sys{
%pythoncode %{

def __str__(self):
    """Prints out spin sys"""

    sss = ""
    for v in self.printstrings():
        sss += str(v) + '\n'

    return (sss)


def __repr__(self):
    return self.__str__()


%}
};


