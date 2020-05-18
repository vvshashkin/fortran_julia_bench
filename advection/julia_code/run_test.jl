using BenchmarkTools
include("advection_mod.jl")

advection.adv(nstep=1,tscheme=:rk4_opt)
println("setup: N=1000, nstep=10, tscheme=:rk4_opt, flux = :up4")
@btime advection.adv(N=1000,nstep=10,tscheme=:rk4_opt,flux=:up4)

println("setup: N=1000, nstep=10, tscheme=:rk4, flux = :up4")
@btime advection.adv(N=1000,nstep=10,tscheme=:rk4,flux=:up4)
