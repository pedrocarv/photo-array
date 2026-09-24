# Sourced by the build_*.sh scripts: f2py options to link LAPACK, which the
# forward model (lyao_rt.f) needs for DGETRF.
#
# Links -llapack directly instead of meson's `--dep lapack`, whose only lookup
# method is pkg-config and conda-forge's liblapack ships no lapack.pc.
#   - Inside a conda env: the env's LAPACK (OpenBLAS on conda-forge), with an
#     rpath so the module finds it at import time.
#   - Otherwise: the system LAPACK (Accelerate on macOS; liblapack-dev or
#     libopenblas-dev on Debian/Ubuntu).
F2PY_LAPACK=(-llapack)
if [ -n "${CONDA_PREFIX:-}" ]; then
    F2PY_LAPACK=(-L"$CONDA_PREFIX/lib" -llapack)
    export LDFLAGS="${LDFLAGS:-} -Wl,-rpath,$CONDA_PREFIX/lib"
fi
