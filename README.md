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

Fortran compilation flags are `ifort -O3 -ipo` and `gfortran -O3 -flto`

@inbounds and @inline are placed here and there in the Julia code to increase speed.

### 4-th order upstream flux test

Parameters: N=1000, tscheme=rk4_opt, flux=up4, nstep=100

#### Fortran

|Compiler                        | ifort | ifort  | ifort | gfortran | gfortran | gfortran  |
|--------------------------------|--------|--------|-------|---------|----------|-----------|
|Code version (see comment below)|   1    |    2   |   3   |      1  |     2    |      3    |
|Run time(s)                     | 6.9    | 4.2    | 3.2    | 16.4   | 14.5     |  5.0      |

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
|  Run time (s)                  | 3.0     | 4.5    | 5.8       | 15.0     |

Code versions:
  1. hard-coded flux
  2. flux as function which returns 1 value

#### Julia

|code verion       |     1   | 2        |     3    |
|------------------|---------|----------|----------|
|runtime (s)       | 4.0     | 4.5      |   ~50    |

Julia code versions:
  1. basic version
  2. hardcoded up4 flux (slower than function, how this is possible?)
  3. using array sections instead of array views as argument to flux function (10x slow down)

### Weno5 test
N=1000, nstep=100, tscheme=k4_opt, flux=weno5
| ifort (function pointer)    | ifort (if/then/else) | gfortran (if/then/else) | julia (basic)  |
|-------------------------------|----------------------|-------------------------|----------------|
| 14.3      | 12.7              | 77.2      |  14.6   |
|           | 12.1 (*see below)  |           |         |

*) -fp-model fast -no-prec-div
