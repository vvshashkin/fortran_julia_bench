#!/bin/bash

echo "test with following setup:"
cat namelist

echo ""
echo "ifort no-objects variant:"
ifort -O3 -ipo -cpp params.f90 operator.f90 time_scheme.f90 main.f90
time ./a.out
rm a.out

echo ""
echo "ifort no-objects variant with hard-coded up4 flux:"
ifort -O3 -ipo -cpp -DHARD_CODE params.f90 operator.f90 time_scheme.f90 main.f90
time ./a.out
rm a.out

echo ""
echo "gfortran no-objects variant:"
gfortran -O3 -cpp -flto params.f90  operator.f90 time_scheme.f90 main.f90
time ./a.out
rm a.out

echo ""
echo "gfortran no-objects variant with hard-coded up4 flux:"
gfortran -O3 -cpp -DHARD_CODE -flto params.f90  operator.f90 time_scheme.f90 main.f90
time ./a.out
rm a.out
