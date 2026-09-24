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
python -m numpy.f2py -c subroutines_lyao.f corona.f global_parameters.f lyao_rt.f -m forward
python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_inv6.f -m los_inv6 --backend meson

```

#### Using 'meson'

Numpy migrates to Meson build system: https://numpy.org/devdocs/f2py/buildtools/distutils-to-meson.html#f2py-meson-distutils

2) Compile the forward and los models

```
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f corona.f global_parameters.f lyao_rt.f -m forward --backend meson
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_inv6.f -m los_inv6 --backend meson
```

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

