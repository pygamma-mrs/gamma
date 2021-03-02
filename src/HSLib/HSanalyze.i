/* HSanalyze.i */

%{
#include "HSLib/HSanalyze.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

// Spin state labels
const std::string alphabeta[7]
                = {"a","b","g","d","e","w","x"};	
 
std::string qStatel(const spin_sys& sys, int bf);

void wf_labels(std::string* wflabels, const spin_sys& sys, gen_op &Op,
                                    double cutoff=1.e-4, int pbf=1, int pfz=1);

void wf_labels(std::string* wflabels, int* filter, const spin_sys& sys,
                     const gen_op &Op, double cut=1.e-4, int pbf=1, int pfz=1);

void wf_labels(std::string* wflabels, const spin_sys& sys, const matrix &B,
                                    double cutoff=1.e-4, int pbf=1, int pfz=1);

void wf_labels(std::string* wflabels, int* filter, const spin_sys& sys,
                      const matrix &B, double cut=1.e-4, int pbf=1, int pfz=1);

void ev_labels(std::string* evlabels, gen_op& Op, double cutoff=1.e-6);

void ev_labels(std::string* evlabels, int* filter, gen_op& Op, double cutoff=1.e-6);

void tref_labels(std::string* trlabels, const spin_sys& sys, gen_op &Op,
                           int type=0, double cut=1.e-4, int pbf=1, int pfz=0);

void tref_labels(std::string* trlabels, const spin_sys& sys, const matrix &B,
                        int type=0, double cutoff=1.e-4, int pbf=1, int pfz=0);

void tran_types(std::string* trtypes, const spin_sys& sys, gen_op &Op,
                                                 int type=0, double cut=1.e-4);

void tran_types(std::string* trtypes, const spin_sys& sys, const matrix &B,
                                              int type=0, double cutoff=1.e-4);

void trev_labels(std::string* trlabels, gen_op& Op, double cutoff=1.e-6);


// Comment out declarations referring to ostream...

//void wavefunctions(std::ostream& ostr, const spin_sys& sys, gen_op &Op,
//                       double cutoff=1.e-4, int pbf=1, int pfz=1, int title=1);

//void wavefunctions(std::ostream& ostr, int* filt, const spin_sys& sys, gen_op &Op,
//                       double cutoff=1.e-4, int pbf=1, int pfz=1, int title=1);

//void eigensystem(std::ostream& ostr, const spin_sys& sys, gen_op& Op,
//      double cute=1.e-6, double cutc=1.e-4, int pbf=1, int pfz=1, int title=1);

//void eigensystem(std::ostream& ostr, int* filt, const spin_sys& sys, gen_op& Op,
//      double cute=1.e-6, double cutc=1.e-4, int pbf=1, int pfz=1, int title=1);

//void transitions(std::ostream& ostr, const spin_sys& sys, gen_op& Op, int type=0,
//      double cutt=1.e-6, double cutc=1.e-4, int pbf=1, int pfz=1, int title=1);

//void transitions(std::ostream& ostr, int* filt, const spin_sys& sys, gen_op& Op, int type=0,
//         double cutt=1.e-6, double cutc=1.e-4, int pbf=1, int pfz=1, int title=1);
 
void ev_select(int* select, gen_op &Op, double val1, int type=0,
                              double val2=0, double cutoff=1.e-4 , int reim=1);
 
void tr_select(int* select, gen_op &Op, double val1, int type=0,
                              double val2=0, double cutoff=1.e-4 , int reim=1);

