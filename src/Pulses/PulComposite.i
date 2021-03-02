/* PulComposite.i */
// Swig interface file.

%{
#include "Pulses/PulComposite.h"
%}

%include "std_string.i"

%rename(__assign__) PulComposite::operator=;

// Note: the comment "///" is used in this file to indicate that
// the commented out line has a method that was already commented
// out before creating this .i file.

class PulWaveform;

class PulComposite : public PulWaveform
{

public:

PulComposite();

PulComposite(const PulWaveform& pulwf,
                                  const spin_system& sys, const std::string& isoch);

PulComposite(const PulWaveform& pulwf,
           gen_op& H, gen_op& FX, gen_op& FY, gen_op& FZ, const std::string& isoch);

PulComposite(const PulWaveform& pulwf,
             const spin_system& sys, const super_op& LOp, const std::string& isoch);

PulComposite(const PulComposite& PT1);

virtual ~PulComposite();

PulComposite & operator = (const PulComposite& CPul1);

gen_op GetH(int i) const;

super_op L0(int i);
super_op GetL0(int i);

super_op Leff(int i);
super_op GetLeff(int i);

//densop SigSS(int i);

virtual HSprop GetU(int i=-1);

virtual HSprop GetU(int i, double td);

virtual HSprop GetU(double td);

virtual HSprop GetUsum(int i=-1);

virtual HSprop GetUmult(int N);

//virtual LSprop GetG(int i=-1);
//LSprop GetG(int i, double td);
//LSprop GetG(double td);
//virtual LSprop GetGsum(int i=-1);
//LSprop GetGmult(int N);



///int        steps()   const;
///std::string     name()    const;
///double     length()  const;

std::string     channel() const;

///row_vector values()  const;

///double  strength(int i) const;
///double  phase(int i)    const;
///double  length(int i)   const;
///complex value(int i)    const;
 
gen_op       FZ()        const;
super_op     ROp()       const;
//densop       SigEq()     const;
double       Precision() const;

///double steps(double td) const;			INHERITED
///int fullsteps(double td) const;			INHERITED
///double sumlength(int i) const;			INHERITED

///int gamB1const() const;
///int phaseconst() const;
///int timeconst() const;

virtual void scalegB1(double sf); 

void setRelax(const spin_system& sys, const super_op& LOp);

void FIDheader(int typ, int rlx=0) const;

void FIDpoint(int typ, int pt, int iWFs, int iSTs) const;

void FIDvalue(int typ, double td, const complex& z) const;

void FIDtell(double SW) const;

virtual double FIDsync(double& SW, int warn=0) const;

virtual int FIDtest(double td, int& nWFs, int& nSTPs, double& tr) const;

row_vector FIDsynchWF(int npts, int nWFs,
                                       gen_op &D, gen_op& sigmap, int track=0);

row_vector FIDsynchST(int npts, int nSTs,
                                       gen_op &D, gen_op& sigmap, int track=0);
 
row_vector FIDsynchFR(int npts, int nFRs,
                                       gen_op &D, gen_op& sigmap, int track=0);

virtual row_vector FID(int N, double td, gen_op &D, gen_op& sp, int track=0);

virtual row_vector FIDRsynchWF(int npts, int nWFs,
                                       gen_op &D, gen_op& sigmap, int track=0);

virtual row_vector FIDRsynchST(int npts, int nSTs,
                                       gen_op &D, gen_op& sigmap, int track=0);

virtual row_vector FIDRsynchFR(int npts, int nFRs,
                                       gen_op &D, gen_op& sigmap, int track=0);

virtual row_vector FIDR(int N, double td, gen_op &D, gen_op& sp, int track=0);

///void GP(int type, int split, int ends) const;         INHERITED

//virtual std::ostream& printEvolve(std::ostream &ostr, double td) const;

//virtual std::ostream& printFID(std::ostream &ostr, double td, int n) const;

//virtual std::ostream& printInfo(std::ostream &ostr) const;

//virtual std::ostream& print(std::ostream &ostr, int full=0) const;

//friend std::ostream &operator << (std::ostream &ostr, const PulComposite &CPul);

};

