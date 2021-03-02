/* PulTrain.h */
// Swig interface file.

%{
#include "Pulses/PulTrain.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

%rename(__assign__) PulTrain::operator=;



class PulTrain : public PulComposite
{

std::string Name;		// Pulse train name
PulCycle PTCyc;			// Pulse train cycle
PulTrainSCyc PTSCyc;	// Pulse train supercycle
int Cycles;				// Flag if cycle defined
int SCycles;		    // Flag if supercycle defined


private:

void PTerror(int eidx, int noret=0) const;

void volatile PTfatality(int error);


public:

PulTrain();
 
PulTrain(const PulComposite& CPul, std::string N="");

PulTrain(const PulComposite& CPul, const PulCycle& Cyc, std::string N="");

PulTrain(const PulComposite& CPul, const PulCycle& Cyc,
                                   const PulSupCycle& SCyc, std::string N="");

PulTrain(const PulTrain& PT1);

~PulTrain();

PulTrain& operator = (const PulTrain& PT1);

HSprop GetU(double td);

//std::ostream& info(std::ostream& ostr, double td) const;
//std::ostream& info(std::ostream& ostr, double td, int npts) const;

row_vector FID(int npts, double td, gen_op &D, gen_op& sigmap);
 
row_vector FIDR(int npts, double td, gen_op &D, gen_op& sigmap);

//std::ostream& printEvolve(std::ostream &ostr, double td, int full=0) const;
//std::ostream& printCycle(std::ostream &ostr, int full=0) const;
//std::ostream& printSCycle(std::ostream &ostr, int full=0) const;
//std::ostream& print(std::ostream &ostr, int full=0) const;
//friend std::ostream &operator << (std::ostream &ostr, const PulTrain &PT);
 
};

