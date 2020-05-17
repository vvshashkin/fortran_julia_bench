struct params_t{T}
    N    :: Int
    Xmax :: T
    dx   :: T
    u    :: Array{T,2}
    v    :: Array{T,2}
    function params_t{T}(Xmax ::T,N::Int,u,v) where(T)
        U = Array{T,2}(undef,N+1,N) #velocity is defined at cell interfaces => N+1 in one of directions
        V = Array{T,2}(undef,N,N+1)

        dx = Xmax/N
        for j=1:N
            y = (j-0.5)*dx
            for i=1:N+1
                x = (i-1)*dx
                U[i,j] = u(x,y)
            end
        end
        for j=1:N+1
            y = (j-1)*dx
            for i=1:N
                x = (i-0.5)*dx
                V[i,j] = v(x,y)
            end
        end

        return new{T}(N,Xmax,dx,U,V)
    end
end
