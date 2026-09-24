#!/bin/bash

# LAPACK link options for the forward model (lyao_rt.f calls DGETRF)
source "$(dirname "$0")/lapack_flags.sh"

# forward
#python -m numpy.f2py -c -m forward subroutines_lyao.f corona.f global_parameters.f lyao_rt.f
# los
#python -m numpy.f2py -c -m los_inv6 subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_inv6.f

# Migration to Meson build system: https://numpy.org/devdocs/f2py/buildtools/distutils-to-meson.html#f2py-meson-distutils
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f corona.f global_parameters.f lyao_rt.f -m forward --backend meson --dep openmp "${F2PY_LAPACK[@]}"
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f corona.f global_parameters_hres.f lyao_rt_hres.f -m forward_hres --backend meson
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_inv6.f -m los_6 --backend meson
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_128.f -m los_128 --backend meson
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_256.f -m los_256 --backend meson
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_512.f -m los_512 --backend meson
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_1024.f -m los_1024 --backend meson
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_16384.f -m los_16384 --backend meson
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_32768.f -m los_32768 --backend meson
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_65536.f -m los_65536 --backend meson
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_131072.f -m los_131072 --backend meson

mv forward* ../
mv los* ../