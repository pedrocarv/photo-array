#!/bin/bash

# Compile a separate LOS executable for the desired viewing geometry (orbitinfo.txt file) and LOS execution.
# This version creates los.out with MAXLOS in the driver_los_lyao.f file with following values:
#     sequential    1x 1024   NFI 1024*1024   WFI 512*512
#     parallel 8    1x  128   NFI  128*1024   WFI 128*512
#     parallel 16   1x   64   NFI   64*1024   WFI  64*512
#     parallel 32   1x   32   NFI   32*1024   WFI  32*512
#     inversion     6 (INV6   MAXLOS = 6)
# It also reads orbitinfo.txt (no more _1x, NFI_parallel etc.)

# forward
#python -m numpy.f2py -c -m forward subroutines_lyao.f corona.f global_parameters.f lyao_rt.f
# los
#python -m numpy.f2py -c -m los_wfi subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_wfi.f

# Migration to Meson build system: https://numpy.org/devdocs/f2py/buildtools/distutils-to-meson.html#f2py-meson-distutils
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f corona.f global_parameters.f lyao_rt.f -m forward --backend meson
FC="gfortran" python -m numpy.f2py -c subroutines_lyao.f lyao_los.f corona.f global_parameters.f driver_los_lyao_wfi.f -m los_wfi --backend meson

mv forward* ../
mv los* ../
