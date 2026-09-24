#!/bin/bash

# -c:        Compile and assemble, but do not link. Compile to an object file, rather than producing a standalone program. This flag is useful if your program source code is split into multiple files. The object files produced by this command can later be linked together into a complete program.
# -std=f95:  Enforces strict compliance with the Fortran 95 standard. This is like -pedantic, except it generates errors instead of warnings.
#            The default value for std is gnu, which specifies a superset of the Fortran 95 standard that includes all of the extensions supported by GNU Fortran, although warnings will be given for obsolete extensions not recommended for use in new code. The legacy value is equivalent but without the warnings for obsolete extensions, and may be useful for old non-standard programs.
# -C:        Do not discard comments. All comments are passed through to the output file, except for comments in processed directives, which are deleted along with the directive. 
# -S:        Compile only; do not assemble or link.
# -o <file>: Place the output into <file>.

# Compile the forward model
gfortran subroutines_lyao.f corona.f global_parameters.f lyao_rt.f -o forward.out

# Compile a separate LOS executable for the desired viewing geometry (orbitinfo.txt file) and LOS execution.
# This version creates los.out with MAXLOS in the driver_los_lyao.f file with following values:
#     inversion     6 (INV6   MAXLOS = 6)
# It also reads orbitinfo.txt (no more orbitinfo_1x.txt, NFI_parallel etc.)
gfortran subroutines_lyao.f lyao_los.f corona.f driver_los_lyao_inv6.f global_parameters.f -o los.out

# mv forward.out ../
# mv los.out ../