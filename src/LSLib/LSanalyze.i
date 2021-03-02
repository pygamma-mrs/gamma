/* LSanalyze.i */

%{
#include "LSLib/LSanalyze.h"
%}

%include "std_string.i"

%feature("autodoc", "1" );

void wf_labels(std::string* wflabels, const spin_sys& sys, super_op &LOp, double cut=1.e-4, int type=-1, int pbf=-1, int pfz=0);

void wf_labels(std::string* wflabels, int* filter, const spin_sys& sys, super_op &LOp,double cut=1.e-4, int type=-1, int pbf=-1, int pfz=0);

void wf_labels(std::string* wflabels, const spin_sys& sys, const matrix &B, const matrix& HB, double cut=1.e-4, int type=-1, int pbf=-1, int pfz=0);

void wf_labels(std::string* wflabels, int* index, const spin_sys& sys, const matrix &B, const matrix& HB, double cut=1.e-4, int type=-1, int pbf=-1, int pfz=0);
 
void ev_labels(std::string* evlabels, super_op& LOp, double cutoff=1.e-6);
 
void ev_labels(std::string* evlabels, int* filter, super_op& LOp, double cutoff=1.e-6);

void wavefunctions(std::ostream& ostr, const spin_sys& sys, super_op &LOp,
         double cutoff=1.e-4, int type=-1, int pbf=-1, int pfz=0, int title=1);
 
void wavefunctions(std::ostream& ostr, int* filter, const spin_sys& sys, super_op &LOp,
              double cutoff=1.e-4, int type=-1, int pbf=-1, int pfz=0, int title=1);
 
void eigensystem(std::ostream& ostr, const spin_sys& sys, super_op& LOp,
                    double cute=1.e-6, double cutc=1.e-4,
                            int type=-1, int pbf=-1, int pfz=0, int title=1);
 
