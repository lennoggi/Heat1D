#include <cmath>
#include "include/Declare_functions.hh"
#include "Parameters.hh"

using namespace std;

// Gaussian
double gaussian(const double &x) {
    constexpr double norm   = 1./(SIGMA*sqrt(2.*M_PI));
    constexpr double expfac = 1./(2*SIGMA*SIGMA);
    const     double xs     = x - X0;
    return norm*exp(-expfac*xs*xs);
}
