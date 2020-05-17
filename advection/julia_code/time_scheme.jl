function rk4!(f::stvec_t,oper,dt::Real)
    k1 = oper(f)
    y = f+0.5dt*k1
    k2 = oper(y)
    y = f+0.5dt*k2
    k3 = oper(y)
    y = f+dt*k3
    k4 = oper(y)

    copy!(f+(dt/6.0)*(k1+2.0k2+2.0k3+k4),f)
    return
end

struct rk4_opt!{T}
    y  :: stvec_t{T}
    k1 :: stvec_t{T}
    k2 :: stvec_t{T}
    k3 :: stvec_t{T}
    k4 :: stvec_t{T}
end

function init_rk4_opt!(f::stvec_t)
    y  = simil(f)
    k1 = simil(f)
    k2 = simil(f)
    k3 = simil(f)
    k4 = simil(f)
    return rk4_opt!(y,k1,k2,k3,k4)
end

function (ts::rk4_opt!)(f::stvec_t,oper!,dt::Real)
    oper!(ts.k1,f)
    lin_comb!(ts.y,f,ts.k1,1.0,0.5dt)
    oper!(ts.k2,ts.y)
    lin_comb!(ts.y,f,ts.k2,1.0,0.5dt)
    oper!(ts.k3,ts.y)
    lin_comb!(ts.y,f,ts.k3,1.0,dt)
    oper!(ts.k4,ts.y)

    lin_comb!(f,f,ts.k1,1.0,dt/6.0)
    lin_comb!(f,f,ts.k2,1.0,dt/3.0)
    lin_comb!(f,f,ts.k3,1.0,dt/3.0)
    lin_comb!(f,f,ts.k4,1.0,dt/6.0)
    return
end
