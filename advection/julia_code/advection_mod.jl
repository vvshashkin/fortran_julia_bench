module advection

using LoopVectorization

include("params.jl")
include("stvec.jl")
include("operator.jl")
include("time_scheme.jl")

function adv(;N=1000,dt=0.0005,Xmax=1.0,nstep=1,tscheme=:rk4,flux=:up4,type=Float64,
             u=(x,y)->1.0,
             v=(x,y)->1.0,
             q=(x,y)->exp(-100*((x-0.5Xmax)^2+(y-0.5Xmax)^2)/Xmax^2))
    params = params_t{type}(convert(type,Xmax),N,u,v)
    f1     = stvec_t{type}(N,initq(params, q))

    if(flux == :up1)
        fluxfun = up1
        halo_width = 1
    elseif(flux == :up4)
        fluxfun = up4
        halo_width = 3
    elseif(flux == :weno5)
        fluxfun = weno5
        halo_width = 3
    else
        println("Error: unknown flux scheme ",flux)
        return
    end

    if(tscheme==:rk4)
        oper = f-> adv_oper(f,params,halo_width,fluxfun)
        ts! = (f1,oper,dt) -> rk4!(f1,oper,dt)
    elseif(tscheme == :rk4_opt)
        oper = (f,p)-> adv_oper!(f,p,params,halo_width,fluxfun)
        ts! = (f,oper!,dt) -> init_rk4_opt!(f1)(f,oper!,dt)
    else
        println("Error: unknown time-steping scheme ",tscheme)
        return
    end

    for it=1:nstep
        ts!(f1,oper,dt)
    end
    return f1
end

function initq(params, q)
    Q = Array{Float64}(undef,params.N, params.N)
    for j=1:params.N
        y = (j-0.5)*params.dx
        for i=1:params.N
            x = (i-0.5)*params.dx
            Q[i,j] = q(x,y)
        end
    end
    return(Q)
end

end
