# 3-part RT - Steps 2 and 3

Follow the steps to compile and run the Step 1 of the 2-part RT code for optically thick analysis of NASA's
Carruthers Geocorona Observatory mission.

## Forward model, Python wrapper

### Inputs

1) **alt_file_LR.dat** - grid, LR - low resolution, 32 values
2) **orbitinfo.txt**
3) **rt_bkg.txt**
4) **thermo_bkg.txt**

### Run

1. Change the directory to the RT folder directory

```
cd /path_to_RT_code_folder
```

#### Using 'distutils'

2) Compile the forward and los models executable through Python.

```
python -m numpy.f2py -c subroutines_lyao.f corona.f global_parameters.f lyao_rt.f -m forward -llapack
python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_inv6.f -m los_inv6 --backend meson

```

These distutils builds are single-threaded (no OpenMP).

#### Using 'meson'

Numpy migrates to Meson build system: https://numpy.org/devdocs/f2py/buildtools/distutils-to-meson.html#f2py-meson-distutils

2) Compile the forward and los models

```
source lapack_flags.sh
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f corona.f global_parameters.f lyao_rt.f -m forward --backend meson --dep openmp "${F2PY_LAPACK[@]}"
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_inv6.f -m los_inv6 --backend meson
```

Build dependencies of the forward model (`lyao_rt.f`):

- **LAPACK is required**: the source-function matrix is factored with LAPACK `DGETRF`. Without it the module
  compiles but fails to import (undefined symbol `dgetrf_`). `lapack_flags.sh` sets the link options in
  `F2PY_LAPACK` (all `build_*.sh` scripts source it):
  - inside a conda env: the env's LAPACK (`-L$CONDA_PREFIX/lib -llapack`, OpenBLAS on conda-forge) plus an rpath;
  - otherwise the system LAPACK (`-llapack`): Accelerate on macOS, `liblapack-dev` or `libopenblas-dev` on
    Debian/Ubuntu.

  It doesn't use meson's `--dep lapack`, which only searches pkg-config and fails in conda envs (conda-forge's
  `liblapack` ships no `lapack.pc`).
- `--dep openmp` is optional and makes the forward model multithreaded (without it, it runs on one core).

The LOS model (`lyao_los.f`) can also be built with `--dep openmp`, as `build_1x.sh` does. Don't use it with the
large fixed-`MAXLOS` drivers (e.g. `driver_los_lyao_65536.f`, `driver_los_lyao_131072.f`) or the hres forward
model: OpenMP places their large local arrays (about 11 x MAXLOS x 8 bytes in the drivers) on the stack, which
overflows it.

The number of threads is `OMP_NUM_THREADS` (default: all cores). `rt_inversion.inversion_init` splits the cores
between its worker processes automatically. Its process pool uses the `spawn` start method on every platform: on
Linux the default `fork` hangs once the parent process has run the OpenMP code (GNU libgomp is not fork-safe), so
scripts that call it need the usual `if __name__ == "__main__":` guard.

On Debian/Ubuntu without conda, install the compilers and LAPACK with
`sudo apt install gfortran liblapack-dev python3-dev` (or `libopenblas-dev` for a faster LAPACK), plus
`pip install numpy meson ninja`.

Coveniently you can run one of the build shell scripts:  **build_1x.sh**,  **build_nfi.sh**,  **build_wfi.sh** or inversion **build_inv6.sh**

```
./build_nfi.sh
```

### Output represents various arrays

**lyao_los_1x.DAT** -  file with the model LOS radiance for geometry of orbitinfo_1x.txt
**lyao_los_NFI.DAT** -  file with the model LOS radiance for geometry of orbitinfo_NFI.txt
**lyao_los_WFI.DAT** -  file with the model LOS radiance for geometry of orbitinfo_WFI.txt

### Troubleshooting

The Fortran doesn't compile with either backend

```
conda install conda-forge::gfortran

or

brew install gcc
gfortran --version

```

Install meson and ninja
Install gpython (gcc): https://fortran-lang.org/learn/os_setup/install_gfortran/

Combination of installing/uninstalling gfortran with conda and brew (Mac) helped. Also linking

