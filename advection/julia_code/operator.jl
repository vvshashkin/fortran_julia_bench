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
@inline function weno5(m::Int,q,u)
  a1 = 0.5*(1+sign(u))
  a2 = 0.5*(1-sign(u))
  a = a1*q[m-2]+a2*q[m+3]
  b = a1*q[m-1]+a2*q[m+2]
  c = a1*q[m  ]+a2*q[m+1]
  d = a1*q[m+1]+a2*q[m]
  e = a1*q[m+2]+a2*q[m-1]

  q1 = (2.0a-7.0b+11.0c)/6.0
  q2 = (-b+5.0c+2.0d)/6.0
  q3 = (2.0c+5.0d-e)/6.0
  s1 = (13.0/12.0)*(a-2.0b+c)^2 + 0.25(a-4.0b+3.0c)^2
  s2 = (13.0/12.0)*(b-2.0c+d)^2 + 0.25(d-b)^2
  s3 = (13.0/12.0)*(c-2.0d+e)^2 + 0.25(3.0c-4.0d+e)^2
  eps = 1e-6
  w1 = 0.1/(eps+s1)^2
  w2 = 0.6/(eps+s2)^2
  w3 = 0.3/(eps+s3)^2
  return u*(w1*q1+w2*q2+w3*q3)/(w1+w2+w3)
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
