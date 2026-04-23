# Heat1D
A suite of codes solving the heat equation
$$
\begin{equation*}
    \frac{\partial u}{\partial t}\left(t, x\right) = \alpha\nabla^2 u\left(t, x\right)\;,
\end{equation*}
$$
where $\alpha > 0$, in one or two dimensions with Dirichlet boundary conditions
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
   ./install/bin/<executable>
   ```
4. To plot evolution snapshots and make movies, see `Utils/README.md` within each subdirectory
