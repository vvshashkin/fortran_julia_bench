function adv_oper(fin::stvec_t, params::params_t,m::Int,fluxf)
    fout = simil(fin)
    adv_oper!(fout,fin,params,m,fluxf)
    return fout
end
function adv_oper!(fout::stvec_t, fin::stvec_t, params::params_t, m::Int, fluxf)

    #m=1  #simulation of final very optimized code
    q = periodic_bc(fin.q,params.N,m)

    flx = Array{eltype(fin.q)}(undef,params.N+1,params.N)
    fly = Array{eltype(fin.q)}(undef,params.N,params.N+1)

    @inbounds for j=1:params.N
        for i=1:params.N+1
            flx[i,j] = fluxf(m,view(q,i:i+2m-1,j+m),params.u[i,j])
            #hardcoded up4, slightly slower (WTF!?)
            #za1 = 0.5+0.5sign(params.u[i,j])
            #za2 = 1.0-za1
            #flx[i,j] = params.u[i,j]*(za1*(3.0q[i+m,j+m]+13.0q[i+m-1,j+m]-5.0q[i+m-2,j+m]+q[i+m-3,j+m])+
            #                          za2*(3.0q[i+m-1,j+m]+13.0q[i+m,j+m]-5.0q[i+m+1,j+m]+q[i+m+2,j+m]))/12.0

            #flx[i,j] = fluxf(m,q[i:i+2m-1,j+m],params.u[i,j])     #10x slow-down
        end
    end
    @inbounds for j=1:params.N+1
        for i=1:params.N
            fly[i,j] = fluxf(m,view(q,i+m,j:j+2m-1),params.v[i,j])

            #hardcoded up4, slightly slower (WTF!?)
            #za1 = 0.5+0.5sign(params.v[i,j])
            #za2 = 1.0-za1
            #fly[i,j] = params.v[i,j]*(za1*(3.0q[i+m,j+m]+13.0q[i+m,j+m-1]-5.0q[i+m,j+m-2]+q[i+m,j+m-3])+
            #                          za2*(3.0q[i+m,j+m-1]+13.0q[i+m,j+m]-5.0q[i+m,j+m+1]+q[i+m,j+m+2]))/12.0

            #fly[i,j] = fluxf(m,q[i+m,j:j+2m-1],params.v[i,j])     #10x slow-down
        end
    end

    @inbounds for j=1:params.N
        for i=1:params.N
            fout.q[i,j] = -(flx[i+1,j]-flx[i,j]+fly[i,j+1]-fly[i,j])/params.dx
        end
    end
    return
end

@inline function up1(m::Int,q,u)
    za1 = 0.5+0.5sign(u)
    za2 = 1.0-za1
    return @inbounds u*(za1*q[m]+za2*q[m+1])
end
@inline function up4(m::Int,q,u)
    za1 = 0.5+0.5sign(u)
    za2 = 1.0-za1
    return @inbounds u*(za1*(3.0q[m+1]+13.0q[m]-5.0q[m-1]+q[m-2])+
                        za2*(3.0q[m]+13.0q[m+1]-5.0q[m+2]+q[m+3]))/12.0
end

function periodic_bc(q,N,m)
    q1 = Array{eltype(q)}(undef,N+2m,N+2m)

    @inbounds q1[m+1:m+N,m+1:m+N] .= q[1:N,1:N]
    @inbounds q1[1:m,m+1:m+N] .= q[N-m+1:N,1:N]
    @inbounds q1[N+m+1:N+2m,m+1:m+N]  .= q[1:m,1:N]
    @inbounds q1[1:N+2m,1:m] .= q1[1:N+2m,N+1:N+m]
    @inbounds q1[1:N+2m,N+m+1:N+2m] .= q1[1:N+2m,m+1:2m]

    return q1
end
