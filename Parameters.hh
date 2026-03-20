#ifndef PARAMETERS_HH
#define PARAMETERS_HH


// ***** DON'T TOUCH *****
#define GAUSSIAN 0
// ***********************

// Initial data
#define INITIAL_DATA GAUSSIAN


// Diffusion constant
constexpr inline double ALPHA = 1.;

// Left and right time-independent function values (Dirichlet BCs)
constexpr inline double F0 = 0.;
constexpr inline double FL = 0.;


// Number of spatial grid points
constexpr inline size_t NX = 1000;

// Spatial bounds are [0, L - L/N]
constexpr inline double L = 1.;


// Number of time steps
constexpr inline size_t NT = 200;

// Time step
constexpr inline double DT = 1.e-04;


// Initial function center
constexpr inline double X0 = 0.3*L;

// Initial wavefunction half-width
constexpr inline double SIGMA = 0.03*L;


// Output frequency
constexpr inline size_t OUT_EVERY = 1;

// Output filename
#define FILENAME "Data.h5"

// Verbosity
constexpr inline bool VERBOSE = true;


#endif  // PARAMETERS_HH
