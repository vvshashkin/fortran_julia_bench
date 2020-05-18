#!/bin/bash

echo "test with following setup:"
cat namelist

echo ""
echo "ifort basic variant:"
ifort -O3 -cpp params.f90 stvec_mod.f90 flux.f90 operator.f90 time_scheme.f90 main.f90 -ipo
time ./a.out
rm a.out

echo ""
echo "ifort with optimized code for improved inlining/ipo:"
ifort -O3 -cpp -DADV_OPT params.f90 stvec_mod.f90 flux.f90 operator.f90 time_scheme.f90 main.f90 -ipo
time ./a.out
rm a.out

echo ""
echo "gfortran basic variant:"
gfortran -O3 -cpp -flto params.f90 stvec_mod.f90 flux.f90 operator.f90 time_scheme.f90 main.f90
time ./a.out
rm a.out

echo ""
echo "gfortran with optimized code for improved inlining/ipo:"
gfortran -O3 -flto -cpp -DADV_OPT params.f90 stvec_mod.f90 flux.f90 operator.f90 time_scheme.f90 main.f90
time ./a.out
rm a.out
