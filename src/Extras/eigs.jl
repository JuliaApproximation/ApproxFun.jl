

Base.eigvals(A::BandedOperator,n::Integer)=eigvals(full(A[1:n,1:n]),full(Conversion(domainspace(A),rangespace(A))[1:n,1:n]))

function Base.eigvals{T<:Operator}(A::Vector{T},n::Int)
    nf = length(A)-1
    for k=1:nf
        @assert isa(A[k],Functional)
    end
    A = promotedomainspace(A,choosedomainspace(A))
    typ = mapreduce(eltype,promote_type,A)
    C = Conversion(domainspace(A[end]),rangespace(A[end]))

    A1,C1 = zeros(typ,n,n),zeros(typ,n,n)
    for k=1:nf
        A1[k,1:n] = A[k][1:n]
    end
    A1[1+nf:end,1:n] = A[end][1:n-nf,1:n]
    C1[1+nf:end,1:n] = C[1:n-nf,1:n]
    ev = sort!(eigvals(A1,C1),by=abs)
    ev[!isinf(ev)]
end
