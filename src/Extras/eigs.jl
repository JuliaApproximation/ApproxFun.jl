

for OP in (:(Base.eigvals),:(Base.eigs))
    @eval $OP(A::Operator,n::Integer;tolerance::Float64=100eps())=$OP([A],n;tolerance=tolerance)
end

Base.eigvals{T<:Operator}(A::Vector{T},n::Integer;tolerance::Float64=100eps())=eigs(A,n;tolerance=tolerance)[1]

function Base.eigs{T<:Operator}(A::Vector{T},n::Integer;tolerance::Float64=100eps())
    nf = length(A)-1
    for k=1:nf
        @assert isafunctional(A[k])
    end
    A = promotedomainspace(A,choosedomainspace(A))
    typ = mapreduce(eltype,promote_type,A)

    ds=domainspace(A[end])
    C = Conversion(ds,rangespace(A[end]))

    A1,C1 = zeros(typ,n,n),zeros(typ,n,n)
    for k=1:nf
        A1[k,1:n] = A[k][1:n]
    end
    A1[1+nf:end,1:n] = A[end][1:n-nf,1:n]
    C1[1+nf:end,1:n] = C[1:n-nf,1:n]

    λ,V=eig(A1,C1)

    pruneeigs(λ,V,ds,tolerance)
end

function pruneeigs(λ,V,ds,tolerance)
    retλ=eltype(λ)[]
    retV=Fun{typeof(ds),eltype(V)}[]
    n=length(λ)
    for k=1:n
        if slnorm(V,n-3:n,k)≤tolerance
            push!(retλ,λ[k])
            push!(retV,Fun(V[:,k],ds))
        end
    end
    retλ,retV
end
