/* SpinOpRot.i */
// Swig interface file.

%{
#include "HSLib/SpinOpRot.h"
%}

%include "std_vector.i"
%include "std_string.i"

%feature("autodoc", "1" );

spin_op Rx(const spin_sys& sys,    int spin,                double beta);
spin_op Rx(const spin_sys& sys,    const std::string& iso,  double beta);
spin_op Rx(const spin_sys& sys,                             double beta);
spin_op Rx(const spin_sys& sys,    const flagvec& flags,    double beta);
spin_op Rx_sp(const spin_sys& sys,                          double beta);

spin_op Ry(const spin_sys& sys,    int spin,               double beta);
spin_op Ry(const spin_sys& sys,    const std::string& iso, double beta);
spin_op Ry(const spin_sys& sys,                            double beta);
spin_op Ry(const spin_sys& sys,    const flagvec& flags,   double beta);
spin_op Ry_sp(const spin_sys& sys,                         double beta);

spin_op Rz(const spin_sys& sys,    int spin,                  double beta);
spin_op Rz(const spin_sys& sys,    const std::string& iso,         double beta);
spin_op Rz(const spin_sys& sys,                               double beta);
spin_op Rz(const spin_sys& sys,    const flagvec& flags, double beta);
spin_op Rz_sp(const spin_sys& sys,                            double beta);

spin_op Raxis(const spin_sys& sys,int I,             double beta, char axis);
spin_op Raxis(const spin_sys& sys,const std::string& iso, double beta, char axis);
spin_op Raxis(const spin_sys& sys,                   double beta, char axis);
spin_op Raxis(const spin_sys& sys,const flagvec& flags,double B,char A);
spin_op Raxis_sp(const spin_sys& sys,                double beta, char axis);

spin_op Rxy(const spin_sys& sys,int spin,              double phi,double beta);
spin_op Rxy(const spin_sys& sys,const std::string& iso,     double phi,double beta);
spin_op Rxy(const spin_sys& sys,                       double phi,double beta);
spin_op Rxy(const spin_sys& sys,const flagvec& flags,double phi,double B);
spin_op Rxy_sp(const spin_sys& sys,                    double phi,double beta);

spin_op Ryz(const spin_sys& sys, int spin,         double theta,double beta);
spin_op Ryz(const spin_sys& sys, const std::string& iso,double theta,double beta);
spin_op Ryz(const spin_sys& sys,                   double theta,double beta);
spin_op Ryz(const spin_sys& sys, const flagvec& flags, double T,double B);
spin_op Ryz_sp(const spin_sys& sys,                double theta,double beta);

spin_op Rzx(const spin_sys& sys, int spin,          double theta, double beta);
spin_op Rzx(const spin_sys& sys, const std::string& iso, double theta, double beta);
spin_op Rzx(const spin_sys& sys,                    double theta, double beta);
spin_op Rzx(const spin_sys& sys, const flagvec& flags, double P, double B);
spin_op Rzx_sp(const spin_sys& sys,                 double phi, double beta);

spin_op Rplane(const spin_sys& sys,int spin, double phi, double beta, char p);
spin_op Rplane(const spin_sys& sys,const std::string& iso,
                                             double phi, double beta, char p);
spin_op Rplane(const spin_sys& sys,          double phi, double beta, char p);
spin_op Rplane(const spin_sys& S,const flagvec& F,
                                             double phi, double beta, char p);
spin_op Rplane_sp(const spin_sys &sys,       double phi, double beta, char p);

spin_op Rxyz(const spin_sys& sys,int spin,double theta,double phi,double beta);
spin_op Rxyz(const spin_sys& sys,const std::string& iso,
                                                  double theta,double phi,double beta);
spin_op Rxyz(const spin_sys& sys,         double theta,double phi,double beta);
spin_op Rxyz(const spin_sys& sys,const flagvec& flags,
                                                  double theta,double phi,double beta);
spin_op Rxyz_sp(const spin_sys& sys,      double theta,double phi,double beta);

spin_op Rspace(const spin_sys & sys, const flagvec& flags,
                                        double theta, double phi, double beta);

spin_op R_Euler(const spin_sys &sys, int spin, double a, double b, double g);
spin_op R_Euler(const spin_sys &sys,const std::string& iso,
                                               double alpha,double beta,double gamma);
spin_op R_Euler(const spin_sys &sys,   double alpha,double beta,double gamma);
spin_op R_Euler_sp(const spin_sys &sys,double alpha,double beta,double gamma);

spin_op R_Euler_plane(const spin_sys &sys, const flagvec& flags,
                                     double alpha, double beta, double gamma); 

spin_op Ixy(const spin_sys &sys, int spin,                  double theta);
spin_op Fxy(const spin_sys &sys, int spin,                  double theta);
spin_op Fxy(const spin_sys &sys, const std::string& iso,    double theta);
spin_op Fxy(const spin_sys &sys,                            double theta);
spin_op Fxy(const spin_sys &sys, const flagvec& flags,      double theta);
spin_op Fxy_sp(const spin_sys &sys,                         double theta);

spin_op Ip(const spin_sys &sys, int spin,                  double theta);
spin_op Fp(const spin_sys &sys, int spin,                  double theta);
spin_op Fp(const spin_sys &sys, const std::string& iso,         double theta);
spin_op Fp(const spin_sys &sys, const flagvec& flags, double theta);
spin_op Fp(const spin_sys &sys,                            double theta);
spin_op Fp_sp(const spin_sys &sys,                         double theta);

spin_op Im(const spin_sys &sys, int spin,                  double theta);
spin_op Fm(const spin_sys &sys, int spin,                  double theta);
spin_op Fm(const spin_sys &sys, const std::string& iso,    double theta);
spin_op Fp(const spin_sys &sys, const flagvec& flags,      double theta);
spin_op Fm(const spin_sys &sys,                            double theta);
spin_op Fm_sp(const spin_sys &sys,                         double theta);

spin_op Fplane(const spin_sys&sys, double theta, char OPtype);

spin_op RotSpinOp(const spin_op& R, const spin_op& F);
