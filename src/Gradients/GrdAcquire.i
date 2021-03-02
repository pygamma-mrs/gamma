/* GrdAcquire.i */
// Swig interface file.

%{
#include "Gradients/GrdAcquire.h"
%}

%include "std_vector.i"
%include "std_string.i"

%feature("autodoc", "1" );

void acquire(std::vector<gen_op>& sigmas, gen_op& D, std::vector<gen_op>& Hs,
                                       double td, int t2pts, row_vector& fid);
void acquire(std::vector<gen_op>& sigmas, gen_op& D, std::vector<gen_op>& Us,
                                                  int t2pts, row_vector& fid);
void FID(std::vector<gen_op>&     sigmas, gen_op& D, std::vector<gen_op>& Hs,
                                       double td, int t2pts, row_vector& fid);
void FID(std::vector<gen_op>&     sigmas, gen_op& D, std::vector<gen_op>& Us,
                                                  int t2pts, row_vector& fid);

row_vector acquire(std::vector<gen_op>& sigmas, gen_op& D, std::vector<gen_op>& Hs,
                                                         double td, int t2pts);
row_vector acquire(std::vector<gen_op>& sigmas, gen_op& D, std::vector<gen_op>& Us,
                                                                    int t2pts);
row_vector FID(std::vector<gen_op>&     sigmas, gen_op& D, std::vector<gen_op>& Hs,
                                                         double td, int t2pts);
row_vector FID(std::vector<gen_op>&     sigmas, gen_op& D, std::vector<gen_op>& Us,
                                                                    int t2pts);

row_vector acquire(RBasic& RB, std::vector<gen_op>& sigmas,
                          gen_op& D, std::vector<gen_op>& Hs, double td, int t2pts);
row_vector FID(RBasic&     RB, std::vector<gen_op>& sigmas,
                          gen_op& D, std::vector<gen_op>& Hs, double td, int t2pts);

row_vector acquire1DT(gen_op& D, std::vector<gen_op>& Hs,
                    std::vector<gen_op>& sigmas, int t2pts, double td, bool norm=1);

row_vector acquire1DT(gen_op& D, std::vector<gen_op>& Hs,
                            gen_op& sigmas, int t2pts, double td, bool norm=1);


complex detect(gen_op& D, std::vector<gen_op>& sigmas);

