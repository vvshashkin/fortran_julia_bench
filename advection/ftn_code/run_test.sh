#!/bin/bash

echo "test with following setup:"
cat namelist

echo ""
ifort -fast -cpp params.f90 stvec_mod.f90 flux.f90 operator.f90 time_scheme.f90 main.f90 &&
echo "ifort basic variant:" &&
time ./a.out

echo ""
ifort -fast -cpp -DADV_OPT params.f90 stvec_mod.f90 flux.f90 operator.f90 time_scheme.f90 main.f90 &&
echo "ifort with optimized code for improved inlining/ipo:" &&
time ./a.out

echo ""
ifort -fast -cpp -DHARD_CODE params.f90 stvec_mod.f90 flux.f90 operator.f90 time_scheme.f90 main.f90 &&
echo "ifort with hardcoded up4 flux:" &&
time ./a.out

echo ""
gfortran -Ofast -cpp params.f90 stvec_mod.f90 flux.f90 operator.f90 time_scheme.f90 main.f90 &&
echo "gfortran basic variant:" &&
time ./a.out

echo ""
gfortran -Ofast -cpp -DADV_OPT params.f90 stvec_mod.f90 flux.f90 operator.f90 time_scheme.f90 main.f90 &&
echo "gfortran with optimized code for improved inlining/ipo:" &&
time ./a.out

echo ""
gfortran -Ofast -cpp -DHARD_CODE params.f90 stvec_mod.f90 flux.f90 operator.f90 time_scheme.f90 main.f90 &&
echo "gfortran with hardcoded up4 flux:" &&
time ./a.out
