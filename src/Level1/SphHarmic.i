/* SphHarmic.i */

%{
#include "Level1/SphHarmic.h"
%}

%feature("autodoc", "1" );

void Y_error(int eidx, int noret=0);

void Y_error(int eidx, const std::string& pname, int noret=0);

void volatile Y_fatality(int error);

double Y00rad();

double Y10rad(double theta);

complex Y11rad(double theta, double phi);

complex Y1m1rad(double  theta, double  phi);

double Y20rad(double theta);

complex Y21rad(double theta, double phi);

complex Y2m1rad(double theta, double phi);

complex Y22rad(double theta, double phi);

complex Y2m2rad(double theta, double phi);

double Y30rad(double theta);

complex Y31rad(double  theta, double  phi);

complex Y3m1rad(double theta, double  phi);

complex Y32rad(double theta, double phi);

complex Y3m2rad(double theta, double phi);

complex Y33rad(double theta, double phi);

complex Y3m3rad(double theta, double phi);

complex Ylmrad(int l, int m, double  theta, double  phi);

double Y00();

double Y10(double theta);

complex Y11(double theta, double phi);

complex Y1m1(double  theta, double  phi);

double Y20(double  theta);

complex Y21(double  theta, double  phi);

complex Y2m1(double theta, double phi);

complex Y22(double theta, double phi);

complex Y2m2(double theta, double phi);

double Y30(double theta);

complex Y31(double theta, double phi);

complex Y3m1(double theta, double phi);

complex Y32(double theta, double phi);

complex Y3m2(double theta, double phi);

complex Y33(double theta, double phi);

complex Y3m3(double theta, double phi);

complex Ylm(int l, int m, double theta, double phi);

