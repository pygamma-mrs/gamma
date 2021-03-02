/* PulCycle.i */
// Swig interface file.

%{
#include "Pulses/PulCycle.h"
%}

%include "std_string.i"

%rename(__assign__) PulCycle::operator=;


class PulCycle : public PulComposite
{

  std::string      CYCname;		// Pulse cycle name
  int         CYCsteps;			// Pulse cycle # of steps
  row_vector  CYCvals;			// Pulse cycle step phase (deg)
  double      CYCtp;			// Pulse cycle length (sec)
  int*        Pindex;			// Pulse cycle propagator indexing
  int         Pcount;			// Pulse cycle propagator count
  int         CUcount;			// Pulse cycle propagator count
  HSprop*     CUsteps;			// Pulse cycle step propagators
  HSprop*     CUsums;			// Pulse cycle summed propagators
  int         CGcount;			// Pulse cycle super propagator count
  LSprop*     CGsteps;			// Pulse cycle step superprops
  LSprop*     CGsums;			// Pulse cycle summed superprops

private:

void CYCerror(int eidx, int noret=0) const;

void volatile CYCfatality(int error) const;

void deleteCIndxs();

void deleteCUprops();

void deleteCGprops();

void copyCIndxs(const PulCycle& CYC1);

void copyCUprops(const PulCycle& CYC1);

void copyCGprops(const PulCycle& CYC1);

void SetCIndxs( );

void SetCUs();

void SetCGs( );

void SetBasis(gen_op& Op);

public:

PulCycle();

PulCycle(const PulComposite& P, const row_vector& S, std::string N);

PulCycle(const PulCycle& PT1);

~PulCycle();

// void operator = (const PulCycle& CYC1);
// reworked to this:
PulCycle& operator= (const PulCycle& CYC1);

//gen_op GetH(PulComposite& PC, int i=-1) const

virtual HSprop GetCU(int i=-1);

LSprop GetCG(int i, int j, double td);

virtual HSprop GetCUsum(int i=-1);

HSprop GetCUsum(int i, int j);

virtual HSprop GetCUmult(int N);

LSprop GetCG(int i=-1);

LSprop GetCGsum(int i=-1);

LSprop GetCGsum(int i, int j);

LSprop GetCGmult(int N);

int          steps()		const;
int          WF_steps()		const;
std::string  name()		const;
std::string  WF_name()		const;
row_vector   values()		const;
row_vector   WF_values()	const;
double       length()		const;
double       WF_length()	const;

complex value(int i)  const;
double  phase(int i)  const;

double steps(double td) const;

double cycles(double td) const;

int fullcycles(double td=-1) const;

void scalegB1(double sf);

row_vector IvsT(int split, int ends, int N=1) const;

row_vector PvsT(int split, int ends, int N=1, double ph=0) const;

void GP(int ty=1, int spl=0, int ed=0, int N=1, double p=0) const;

void FM(int ty=1, int spl=0, int ed=0, int N=1, double p=0) const;

double FIDsync(double& SW) const;

int FIDtest(double td, int& nCYs, int& nWFs, int& nSTPs, double& tr) const;

row_vector FIDsynchCYC(int npts, int nCYs,
                                       gen_op &D, gen_op& sigmap, int track=0);

row_vector FIDWFsynch(int npts, int nWFs,
                                       gen_op &D, gen_op& sigmap, int track=0);

row_vector FIDSTsynch(int npts, int nSTs,
                                       gen_op &D, gen_op& sigmap, int track=0);

row_vector FID(int N, double td, gen_op &D, gen_op& sp, int F=0);

virtual row_vector FIDRsynchCYC(int npts, int nCYs,
                                       gen_op &D, gen_op& sigmap, int track=0);

virtual row_vector FIDRWFsynch(int npts, int nWFs,
                                       gen_op &D, gen_op& sigmap, int track=0);

virtual row_vector FIDRSTsynch(int npts, int nSTs,
                                       gen_op &D, gen_op& sigmap, int track=0);

row_vector FIDR(int N, double td, gen_op &D, gen_op& sp, int F=0);

//std::ostream& printEvolve(std::ostream &ostr, double td) const;
//virtual std::ostream& printFID(std::ostream &ostr, double td, int npts) const;
//std::ostream& printSteps(std::ostream &ostr) const;
//std::ostream& printInfo(std::ostream &ostr) const;
//std::ostream& printBase(std::ostream &ostr) const;
//std::ostream& print(std::ostream &ostr, int full=0) const;
//friend std::ostream &operator << (std::ostream &ostr, const PulCycle &CYC);

};

