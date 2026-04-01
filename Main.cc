#include <cassert>
#include <array>
#include <sstream>
#include <iostream>

#include <hdf5.h>

#include "include/Declare_functions.hh"
#include "include/Check_parameters.hh"

#include "Parameters.hh"

using namespace std;


int main() {
    // Print some info
    #if (INITIAL_DATA == GAUSSIAN)
    cout << "Gaussian initial data"     << endl
         << "Initial center: " << X0    << endl
         << "Initial sigma:  " << SIGMA << endl
         << endl;
    #else
    #error "Invalid initial data"
    #endif


    // Build the initial data
    constexpr double dx = L/static_cast<double>(NX);
    array<double, NX> x, u, diag, rhs;

    for (auto i = decltype(NX){0}; i < NX; ++i) {
        x[i] = i*dx;  // NOTE: only used for output

        #if (INITIAL_DATA == GAUSSIAN)
        u[i] = gaussian(x[i]);
        #else
        #error "Invalid initial data"
        #endif
    }


    /* ------------------------
     * Initialize the HDF5 file
     * ------------------------ */
    const auto file_id = H5Fcreate(FILENAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    assert(file_id >= 0);

    constexpr hsize_t fdims = NX;
    const auto fspace_id    = H5Screate_simple(1, &fdims, nullptr);
    assert(fspace_id >= 0);


    // ***** Output frequency ****
    constexpr hsize_t one   = 1;
    const auto space_one_id = H5Screate_simple(1, &one, nullptr);
    assert(space_one_id >= 0);

    const auto out_freq_dset_id = H5Dcreate(file_id, "Output frequency", H5T_NATIVE_INT, space_one_id,
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(out_freq_dset_id >= 0);
    assert(H5Dwrite(out_freq_dset_id, H5T_NATIVE_INT, space_one_id, space_one_id, H5P_DEFAULT, &OUT_EVERY) >= 0);
    assert(H5Dclose(out_freq_dset_id) >= 0);
    assert(H5Sclose(space_one_id) >= 0);


    // ***** Position *****
    const auto x_dset_id = H5Dcreate(file_id, "Position",  H5T_NATIVE_DOUBLE, fspace_id,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(x_dset_id >= 0);
    assert(H5Dwrite(x_dset_id, H5T_NATIVE_DOUBLE, fspace_id, fspace_id, H5P_DEFAULT, x.data()) >= 0);
    assert(H5Dclose(x_dset_id) >= 0);


    // ***** Solution *****
    const auto u_group_id = H5Gcreate(file_id, "Solution",
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(u_group_id >= 0); 

    const auto u0_dset_id = H5Dcreate(u_group_id, "Iteration 0",
                                      H5T_NATIVE_DOUBLE, fspace_id,
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(u0_dset_id >= 0);
    assert(H5Dwrite(u0_dset_id, H5T_NATIVE_DOUBLE, fspace_id, fspace_id, H5P_DEFAULT, u.data()) >= 0);
    assert(H5Dclose(u0_dset_id) >= 0);



    /* ====================
     * Begin time evolution
     * ==================== */
    constexpr auto nxm1 = NX - 1;
    constexpr auto nxm2 = NX - 2;

    constexpr double r           = ALPHA*DT/(dx*dx);
    constexpr double one_plus_2r = 1. + 2.*r;
    constexpr double r_u0        = r*U0;
    constexpr double r_uL        = r*UL;

    constexpr hsize_t tdims = NT/OUT_EVERY;  // NOTE: integer division
    array<double, tdims> time;

    // n == 0 has been written to file already
    for (auto n = decltype(NT){1}; n < NT; ++n) {
        diag[0]     = one_plus_2r;
         rhs        = u;
         rhs[0]    += r_u0;
         rhs[nxm1] += r_uL;

        /* Thomas' tridiagonal algorithm to solve the linear system
         *   / b[0] c[0]  0    0   ... \ / u[0] \   / d[0] \
         *   | a[1] b[1] c[1]  0   ... | | u[1] |   | d[1] |
         *   |  0   a[2] b[2] c[2] ... | | u[2] |   | d[2] |
         *   |  0    0   a[3] b[3] ... | | u[3] | = | d[3] |
         *   |  .    .    .    .    .  | |  .   |   |  .   |
         *   |  .    .    .    .    .  | |  .   |   |  .   |
         *                                                                      */
        // Forward elimination
        for (auto i = decltype(NX){1}; i < NX; ++i) {
            const auto m = -r/diag[i-1];   // m = a[i]/b[i-1]
            diag[i]  = one_plus_2r + m*r;  // b[i] -= m*c[i-1]
             rhs[i] -= m*rhs[i-1];         // d[i] -= m*d[i-1]
        }

        // Backward solve
        u[nxm1] = rhs[nxm1]/diag[nxm1];
        for (auto i = decltype(NX){nxm2}; i > 0; --i) {
            u[i] = (rhs[i] + r*u[i+1])/diag[i];  // u[i] = (d[i] - c[i]*x[i+1])/b[i]
        }


        // Write data to file
        static_assert(sizeof(std::size_t) <= sizeof(long long),
                      "size_t too big for long long");
        const auto [quot, rem] = lldiv(static_cast<long long>(n), static_cast<long long>(OUT_EVERY));
        if (rem == 0) {
            time[quot] = n*DT;

            ostringstream iteration_ss;
            iteration_ss << "Iteration " << n;

            const auto un_dset_id = H5Dcreate(u_group_id, iteration_ss.str().c_str(),
                                              H5T_NATIVE_DOUBLE, fspace_id,
                                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            assert(un_dset_id >= 0);
            assert(H5Dwrite(un_dset_id, H5T_NATIVE_DOUBLE, fspace_id, fspace_id, H5P_DEFAULT, u.data()) >= 0);
            assert(H5Dclose(un_dset_id) >= 0);

            if constexpr (VERBOSE) {
                cout << "Iteration " << n << ": data written to file" << endl;
            }
        }

        if constexpr (VERBOSE) {
            cout << "Iteration " << n << ": done" << endl;
        }
    }


    // Write time to file 
    const auto tspace_id = H5Screate_simple(1, &tdims, nullptr);
    assert(tspace_id >= 0);

    const auto t_dset_id = H5Dcreate(file_id, "Time", H5T_NATIVE_DOUBLE, tspace_id,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(t_dset_id >= 0);

    assert(H5Dwrite(t_dset_id, H5T_NATIVE_DOUBLE, tspace_id, tspace_id, H5P_DEFAULT, time.data()) >= 0);

    assert(H5Dclose(t_dset_id) >= 0);
    assert(H5Sclose(tspace_id) >= 0);
 

    // Close up
    assert(H5Sclose(fspace_id) >= 0);
    assert(H5Gclose(u_group_id) >= 0);
    assert(H5Fclose(file_id) >= 0);

    cout << "Data written to '" << FILENAME << "'" << endl;

    return 0;
}
