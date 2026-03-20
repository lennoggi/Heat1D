#ifndef CHECK_PARAMETERS_HH
#define CHECK_PARAMETERS_HH

#include "../Parameters.hh"


static_assert(INITIAL_DATA == GAUSSIAN);

static_assert(ALPHA > 0.);

static_assert(NX > 1);
static_assert(L  > 0.);

static_assert(NT > 1, "Need at least two time steps for correct output handling");
static_assert(DT > 0.);

static_assert(SIGMA > 0.);
static_assert(X0 - 3.*SIGMA > 0. and X0 + 3*SIGMA < L - L/static_cast<double>(NX),
              "Initial function is quite out of bounds");

static_assert(OUT_EVERY > 0);
static_assert(VERBOSE or not VERBOSE);


#endif  // CHECK_PARAMETERS_HH
