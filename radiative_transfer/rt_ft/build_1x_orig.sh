#!/bin/bash

# forward
#python -m numpy.f2py -c -m forward subroutines_lyao.f corona.f global_parameters.f lyao_rt.f
# los
#python -m numpy.f2py -c -m los_inv6 subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_inv6.f

# Migration to Meson build system: https://numpy.org/devdocs/f2py/buildtools/distutils-to-meson.html#f2py-meson-distutils
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f corona.f global_parameters.f lyao_rt.f -m forward --backend meson --dep openmp --dep lapack
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_512.f -m los_1x --backend meson

mv forward* ../
mv los* ../