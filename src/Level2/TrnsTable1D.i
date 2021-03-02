/* TrnsTable1D.h */

%{
#include "Level2/TrnsTable1D.h"
%}

%include "std_vector.i"
%include "std_string.i"
%include "std_list.i"

%rename(__assign__) TTable1D::operator=;

namespace std {
   %template(StringVector) vector<string>;
   %template(IntVector)    vector<int>;
   %template(DoubleVector) vector<double>;
}


class TTable1D: private matrix 
{

public:

TTable1D();
//TTable1D(const matrix& mx);
//TTable1D(const matrix& mx, int warn);
TTable1D(const TTable1D& TTab1);

~TTable1D();

TTable1D& operator= (const TTable1D& TTab1); 
double center(bool wa=true);
void   offset(double F, int inHz=1);
void   offset(double F, int tr, int inHz);
void   FRscale(double Fscf);
void   FRscale(double Fscf, int tr);
void   BC(double res=-1.0);
 
void Iscale(double         Iscf);
void Iscale(double         Iscf, int tr);
void Iscale(const complex& Iscf);
void Iscale(const complex& Iscf, int tr);
void Iremove(double        dcut=-1.0);

void broaden(double LWR, int inHz=1);
void broaden(double LWR, int tr, int inHz);

void resolution(double res);

row_vector              T(int npts,         double tinc) const;
void                    T(row_vector& data, double tinc) const;
std::vector<row_vector> Ts(int npts,        double tinc) const;
std::vector<int>        TCutoffs(int npts,  double tinc) const;

row_vector              F(int npts,        double Fst,double Ffi) const;
void                    F(row_vector& data,double Fst,double Ffi) const;
std::vector<row_vector> Fs(int npts,       double Fst,double Ffi) const;
row_vector FD(int N,            double fstart, double fend) const;
void       FD(row_vector& data, double fstart, double fend) const;

complex pcorrect(double& w0,    double w1, int order);
void    pcorrect(double Wpivot, complex& P);

double     R2(int tr) const;
double     Fr(int tr) const;
complex    I(int  tr) const;
row_vector Tr(int tr) const;

bool   LineWidths()  const;
bool   Intensities() const;
bool   Phases()      const;
int    size()        const;
double FRmax()       const;
double FRmin()       const;
double Tdmin()       const;
double LWmax()       const;
double LWmin()       const;
double Imax()        const;
double Noisemax()    const;
//matrix mx()	         const;

//friend TTable1D sum(const TTable1D& TT1, const TTable1D& TT2, double res=1.e-6);

std::vector<int> Sort(int k, int type, int colf) const;
 
bool readPSet(const std::string& filein, int indx=-1, int warn=2);
bool readPSet(const ParameterSet& pset,  int indx,    int warn=2);

void setType(int     typ);
void setSort(int     sf);
void setConv(double  cf);
void setIcut(double  ct);
void setInorm(double inorm=0.0);
void setSN(double    S2N);
void setHprint(int   hp);
void setRprint(int   rp);
void setLWprint(int  lwp);
void setT2print(int  t2p);
void setPHprint(int  php);
void setFreqRev();
 
int    getType()    const;
int    getSort()    const;
double getConv()    const;
double getIcut()    const;
double getInorm()   const;
double getSN()      const;
int    getHprint()  const;
int    getRprint()  const;
int    getLWprint() const;
int    getT2print() const;
int    getPHprint() const;
bool   getFreqRev() const;

std::vector<std::string> printStrings()           const;
//std::ostream&            print(std::ostream& out) const;

//friend std::ostream& operator << (std::ostream &ostr, const TTable1D& TTab);
//std::ostream& status(std::ostream& ostr) const;

void           write(const std::string& fn) const;
//std::ofstream& write(std::ofstream&     fp) const;

void        dbwrite_old(const std::string& fileName, 
                        const std::string& compname,  // metabolite name			 
                        const double& lowppm, 
			            const double& highppm, 
			            const double& specfreq, 
			            const double& reffreq,
			            const int& loop,
						const std::vector<std::string> & headerLines) const;

void        dbwrite(     const std::string& fileName, 
						 const std::string& compname, 
						 const double& specfreq,
						 const int& numberspins,
						 const int& loop,
						 const std::vector<std::string> & header) const;						
                            
unsigned int  calc_spectra( std::vector<double> & freq,
							std::vector<double> & ampl,
							std::vector<double> & phase,
							double specfreq,
							int numberspins,                                       
							double freqtol = 0.001,
							double phasetol = 45.0,
							double lowppm = -1000.0,
							double highppm = 1000.0) const;
                            
unsigned int  calc_spectra2(std::vector<double> & freq,
							std::vector<double> & ampl,
							std::vector<double> & phase,
							double specfreq,
							int numberspins,                                       
							double freqtol = 0.001,
							double phasetol = 45.0,
							double lowppm = -1000.0,
							double highppm = 1000.0,
                            double normal = 0.0,
                            bool verbose = false) const;                            
 
void           read(const std::string& fn);
//std::ifstream& read(std::ifstream&     fp);
 
//std::ostream& printT(std::ostream& ostr, double tinc, int npts, int P2P=0)
 
//std::ostream& printF(std::ostream& ostr, int npts, double Fst, double Ffi);

// sosi  this is deprecated.... ever used?
//friend void offset(matrix& mx, double F, double LWR, int inHz=0);
};

