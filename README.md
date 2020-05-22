# fortran_julia_bench
The idea of this project is to compare the speed of computations for two popular high-performance computing languages
using the simple hydrodynamic-type problem.
Also to estimate the overhead introduced by using derived-types and type-bound procedures instead of simple
array-based fortran code.

## The test problem
The test problem is 2D passive tracer advection equation written in the flux-form:

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial q}{\partial t}=-\nabla(\vec{v}q),">

where <img src="https://render.githubusercontent.com/render/math?math=q"> is the tracer density, <img src="https://render.githubusercontent.com/render/math?math=\vec{v}"> is constant horizontal wind field.

The domain is square of 1x1 size with biperiodic boundary conditions.

The problem can be considered as the reduction of real hydrodynamics-solver code for benchmarking purposes, also one can think of it as a practical task for numerical-methods-in-hydrodynamics students.

## Solution method
The advection equation is discretized in space on the quadrilater grid with NxN cells. Finite-differences / finite-volume method gives the following equation for the evolution of tracer density in i,j-th cell:

<img src="https://render.githubusercontent.com/render/math?math=\frac{\partial q_{i,j}}{\partial t}=-\frac{(F_{i %2B 1/2,j}-F_{i-1/2,j} %2B F_{i,j %2B 1/2}-F_{i,j-1/2})}{\Delta x^2},">

where <img src="https://render.githubusercontent.com/render/math?math=F=\vec n\cdot\vec v q"> is a mass-flux at cell interface with normal <img src="https://render.githubusercontent.com/render/math?math=\vec n">. Very different flux schemes can be used. Currently we implemented 1-st and 4-th order accurate upstream fluxes and weno5 flux.

Time integration of the discrete problem can be performed by any suitable time-stepping scheme. RK4 scheme is curently used.

## Code structure
The main requirement for the code is to flexible, that means the following:
- the run-time choice of numerical schemes must be controled externally by setting the configuration parameters,
- new numerical schemes (eg. flux schemes) can be implemented easily.

To meet these requirements the following abstractions are separated:
|    Abstraction          |        Fortran implementation     |      Julia implementation       |
|-------------------------|-----------------------------------|---------------------------------|
|Domain parameters        | params derived type               |  params struct                  |
| Time integration scheme | derived type with step method     | function/callable struct q1 = step(q0, operator,dt)        |
| RHS of the problem      | operator_t derived type           | rhs=operator(q0)                |
| State vector            | derived type with basic arithmetic operations defined (extentions of stvec abstract class) | struct with basic arithmetic operations defined                |
| Flux                    | function compatible with clflux abstract function interface |  function                      |

The specific choice of time integration and flux schemes is controlled by namelist parameters (Fortran) and symbol-type arguments passed to the main function of Julia code. Abstract types and interfaces in Fortran code allows simple implementation and usage of different time-stepping schemes and fluxes. Duck-typing helps Julia to be flexible without explicitly defining interfaces.

The resulting Fortran code contains ~650 lines, Julia code is shorter (~240 lines). Non-flexible fortran code (array processing, no derived type and abstract interfaces has ~170 lines).

## Results
### General remarks
Julia run-time is measured by @btime. Fortran run time is measures by runing

```
time ./a.out
```

3-5 times and averaging the results (user+system).

Fortran compilation flags are `ifort -fast` (the same as `-O3 -ipo -no-prec-div -fp-model fast=2`) and `gfortran -Ofastgfortran`

@inbounds and @inline are placed here and there in the Julia code to increase speed.
```bash
$ ifort --version
ifort (IFORT) 19.1.0.166 20191121
Copyright (C) 1985-2019 Intel Corporation.  All rights reserved.

$ gfortran --version
GNU Fortran (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0
Copyright (C) 2017 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

$ julia --version
julia version 1.3.1
```

### 4-th order upstream flux test

Parameters: N=1000, tscheme=rk4_opt, flux=up4, nstep=100

#### Fortran

|Compiler                        | ifort | ifort  | ifort | gfortran | gfortran | gfortran  |
|--------------------------------|--------|--------|-------|---------|----------|-----------|
|Code version (see comment below)|   1    |    2   |   3   |      1  |     2    |      3    |
|Run time(s)                     | 6.2    | 3.9    | 2.6   | 14.1    | 14.1     |  3.0      |

Code versions:
  1. flux as function pointer in the operator type:
```fortran
call flux_conv(fout%p, fin%p, params, this%fluxfun, this%hw)
```
  2. flux is explicitly passed to the flux_conv subroutine:
```fortran
            if(this%flux_scheme_name == "up1") then
                call flux_conv(fout%p, fin%p, params, up1f, 1)
            else if(this%flux_scheme_name == "up4") then
                call flux_conv(fout%p, fin%p, params, up4f, 3)
            else if(this%flux_scheme_name == "weno5") then
                call flux_conv(fout%p, fin%p, params, weno5f, 3)
            end if
```
  3. hardcoded up4 flux

#### Non-object oriented fortran code

|Compiler                        |  ifort1 | ifort2 | gfortran1 | gfortran2|
|--------------------------------|---------|--------|-----------|----------|
|Code version (see comment below)|   1     |    2   |   1       |      2   |
|  Run time (s)                  | 2.5     | 3.8    | 2.8       | 12.0     |

Code versions:
  1. hard-coded flux
  2. flux as function which returns 1 value

#### Julia

|code verion       |     1   | 2        |     3    |
|------------------|---------|----------|----------|
|runtime (s)       | 3.3     | 3.3      |   ~50    |

Julia code versions:
  1. basic version
  2. hardcoded up4 flux
  3. using array sections instead of array views as argument to flux function (10x slow down)

### Weno5 test
N=1000, nstep=100, tscheme=k4_opt, flux=weno5
| ifort (function pointer)    | ifort (if/then/else) | gfortran (if/then/else) | julia (basic)  |
|-------------------------------|----------------------|-------------------------|----------------|
| 12.1      | 9.9               | 23.3      |  13.8   |

*) -fp-model fast -no-prec-div

### Results summary
Currently it looks like Julia runs indeed! The speed of computations in this problems is almost the same as for Fortran compiled with ifort. Actually, the considered problem tests the efficiency of small-function inlining and optimization. The performance critical code is the following:
```julia
    @inbounds for j=1:params.N
        for i=1:params.N+1
            flx[i,j] = fluxf(m,view(q,i:i+2m-1,j+m),params.u[i,j])
            #flx[i,j] = fluxf(m,q[i:i+2m-1,j+m],params.u[i,j])     #10x slow-down
        end
    end
```
```julia
@inline function up4(m::Int,q,u)
    za1 = 0.5+0.5sign(u)
    za2 = 1.0-za1
    return @inbounds u*(za1*(3.0q[m+1]+13.0q[m]-5.0q[m-1]+q[m-2])+
                        za2*(3.0q[m]+13.0q[m+1]-5.0q[m+2]+q[m+3]))/12.0
end
```
Fortran:
```fortran
        do j=1,N
            do i=0,N
                flx(i,j) = fluxf(hw,q(i-hw+1:i+hw,j),params%u(i,j))
            end do
        end do
```
```fortran
    real(kind=8) function up4f(hw, q, u) result(fl)
        integer(kind=4), intent(in) :: hw
        real(kind=8),    intent(in) :: q(-hw+1:hw)
        real(kind=8),    intent(in) :: u

        real(kind=8)             :: za1, za2

        za1 = .5_8+sign(.5_8,u)
        za2 = 1._8-za1
        fl = u*(za1*(3.0_8*q(1)+13.0_8*q(0)-5.0_8*q(-1)+q(-2))+ &
                za2*(3.0_8*q(0)+13.0_8*q(1)-5.0_8*q( 2)+q( 3)))/12.0_8
    end function up4f
```
Julia inliner is very good, it manages to generate code as fast as hard-coded version of up4 (and even faster). The ifort small-function code is 30-50% slower as compared to the hardcoded version.

Fortran old style non-object oriented code is not significantly faster than object-oriented one with ifort, unless the flux is implemented as derived type (it kills all interprocedural optimization). GNU fortran is only slightly slower than Intel fortran for simple code and completely fails optimization when faced with some problems like small-function inlining and/or all these derived types :(.

(Was not shown) Attempt to implement "natural" stvec arithmetics (`f1=a*f2+b*f3`) with fortran resulted in about 3 seconds slower code (up4 test) as compared to less elegant but robust `call f%lincomb3(f1,f2,a,b)`. Try tscheme="rk4" instead of "rk4_opt" to compare.

## Running the tests
To run the tests on your machine just clone the repository and run `./run_test.sh` in advection/ftn_code directory and `julia run_test.jl` in advection/julia_code. If you don't have ifort, please, modify run_test.sh file accordingly.
