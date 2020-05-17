function adv_oper(fin::stvec_t, params::params_t,m::Int,fluxf)
    fout = simil(fin)
    adv_oper!(fout,fin,params,m,fluxf)
    return fout
end
function adv_oper!(fout::stvec_t, fin::stvec_t, params::params_t, m::Int, fluxf)

    q = periodic_bc(fin.q,params.N,m)

    flx = Array{eltype(fin.q)}(undef,params.N+1,params.N)
    fly = Array{eltype(fin.q)}(undef,params.N,params.N+1)

    for j=1:params.N
        for i=1:params.N+1
            flx[i,j] = fluxf(m,q[i:i+2m-1,j+m],params.u[i,j])
        end
    end
    for j=1:params.N+1
        for i=1:params.N
            fly[i,j] = fluxf(m,q[i+m,j:j+2m-1],params.v[i,j])
        end
    end

    for j=1:params.N
        for i=1:params.N
            fout.q[i,j] = -(flx[i+1,j]-flx[i,j]+fly[i,j+1]-fly[i,j])/params.dx
        end
    end
    return
end

function up1(m,q,u)
    za1 = 0.5+0.5sign(u)
    za2 = 1.0-za1
    return u*(za1*q[m]+za2*q[m+1])
end
function up4(m,q,u)
    za1 = 0.5+0.5sign(u)
    za2 = 1.0-za1
    return u*(za1*(3.0q[m+1]+13.0q[m]-5.0q[m-1]+q[m-2])+
              za2*(3.0q[m]+13.0q[m+1]-5.0q[m+2]+q[m+3]))/12.0
end

function periodic_bc(q,N,m)
    q1 = Array{eltype(q)}(undef,N+2m,N+2m)

    q1[m+1:m+N,m+1:m+N] .= q[1:N,1:N]
    q1[1:m,m+1:m+N] .= q[N-m+1:N,1:N]
    q1[N+m+1:N+2m,m+1:m+N]  .= q[1:m,1:N]
    q1[1:N+2m,1:m] .= q1[1:N+2m,N+1:N+m]
    q1[1:N+2m,N+m+1:N+2m] .= q1[1:N+2m,m+1:2m]

    return q1
end
