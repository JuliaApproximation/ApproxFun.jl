

## This file has methods to stably divide the jacobi singularity


function dirichlet_divide_singularity{T<:Number}(b::Bool,v::Vector{T})
    n=length(v)
    w=zeros(T,n-1)
    s=b?1:-1

    if n==1
        return w
    elseif n==2
        w[1]=-s*v[2]
        return w
    end

    w[n-1]=-2s*v[n]

    for k=n-2:-1:2
        @inbounds w[k]=-2s*(v[k+1] - .5w[k+1])
    end

    @inbounds w[1]=-s*(v[2] - .5w[2])


    w
end





# dirichletrange_divide_singularity divides by (1-x^2)
# where the coefficients are ChebyshevDirichlet{1,1}
# assume that v[1] == v[2] ==0.
# the output coefficients are Chebyshev


function dirichlet_divide_singularity{T<:Number}(v::Vector{T})
    n=length(v)
    if n ≤ 2    # assumes v[1]==v[2]==0 which is deleted
        return zeros(T,1)
    end

    w=zeros(T,n-2)



    w[n-2]=-4v[n]
    if n>3
        w[n-3]=-4v[n-1]
    end

    for k=n-4:-1:2
        @inbounds w[k]=-4*(v[k+2] - .25w[k+2])
    end

    w[1]=-2v[3]
    if n≥5
        w[1]+=0.5w[3]
    end

    w
end


divide_singularity(s::Bool,f::Fun)=Fun(dirichlet_divide_singularity(s,Fun(f,s?ChebyshevDirichlet{0,1}:ChebyshevDirichlet{1,0}).coefficients),Chebyshev(domain(f)))
divide_singularity(f::Fun)=Fun(dirichlet_divide_singularity(Fun(f,ChebyshevDirichlet{1,1}).coefficients),Chebyshev(domain(f)))

function divide_singularity(s::@compat(Tuple{Int,Int}),f::Fun)
    if s[1]>0 && s[2]>0
        divide_singularity((s[1]-1,s[2]-1),divide_singularity(f))
    elseif s[1]>0
        divide_singularity((s[1]-1,s[2]),divide_singularity(false,f))
    elseif s[2]>0
        divide_singularity((s[1],s[2]-1),divide_singularity(true,f))
    else
        @assert s[1]==s[2]==0
        f
    end
end


