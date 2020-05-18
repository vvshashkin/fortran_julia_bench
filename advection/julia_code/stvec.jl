struct stvec_t{T}
    N :: Int
    q :: Array{T,2}
end

Base.:+(f1::stvec_t{T},f2::stvec_t{T}) where(T) = stvec_t{T}(f1.N,f1.q .+ f2.q)

Base.:-(f1::stvec_t{T},f2::stvec_t{T}) where(T) = stvec_t{T}(f1.N,f1.q .- f2.q)

Base.:*(a::Real,f::stvec_t{T}) where(T) = stvec_t{T}(f.N,a.*f.q)
Base.:*(f::stvec_t{T},a::Real) where(T) = stvec_t{T}(f.N,a.*f.q)

Base.:/(f::stvec_t{T},a::Float64) where(T) = stvec_t{T}(f.N,f.q ./ a)

function lin_comb!(f0::stvec_t,f1::stvec_t, f2::stvec_t, a::Real, b::Real)
    f0.q .= a.*f1.q.+b.*f2.q
    return
end

function copy!(fin::stvec_t, fout::stvec_t)
    fout.q .= fin.q
    return
end

function simil(f::stvec_t{T}) where(T)
    return stvec_t{T}(f.N,similar(f.q))
end
