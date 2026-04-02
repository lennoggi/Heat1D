# Heat1D
A code solving the 1D heat equation
$$
\begin{equation*}
    \frac{\partial u}{\partial t}\left(t, x\right) = \alpha\frac{\partial u}{\partial x^2}\left(t, x\right)\;,
\end{equation*}
$$
where $\alpha > 0$, with Dirichlet boundary conditions
$$
\begin{align*}
    u\left(t, 0\right) &= u_0\\
    u\left(t, L\right) &= u_L
\end{align*}
$$
by implicit time integration (backward Euler).

## Minimal requirements
- A C++ compiler supporting the `c++17` standard
- `cmake`
- The HDF5 library
- `python3` with `numpy`, `matplotlib`, and `h5py` to generate evolution snapshots
- `ffmpeg` to generate movies

## Usage
1. Tune the parameters in `Parameters.hh`
2. Compile:
   ```
   ./build.sh
   ```
3. Run:
   ```
   ./install/bin/Heat1D_exe
   ```
4. To plot evolution snapshots and make movies, see `Utils/README.md`
