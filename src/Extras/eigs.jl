


Base.eigvals(A::Operator,n::Integer;tolerance::Float64=100eps()) =
    eigs(A,n;tolerance=tolerance)[1]

function Base.eigs(A::Operator,n::Integer;tolerance::Float64=100eps())
    typ = eltype(A)

    ds=domainspace(A)
    C = Conversion(ds,rangespace(A))

    A1,C1 = zeros(typ,n,n),zeros(typ,n,n)
    A1[1:end,1:n] = A[1:n,1:n]
    C1[1:end,1:n] = C[1:n,1:n]

    λ,V=eig(A1,C1)

    pruneeigs(λ,V,ds,tolerance)
end

Base.eigvals(Bcs::Operator,A::Operator,n::Integer;tolerance::Float64=100eps()) =
    eigs(Bcs,A,n;tolerance=tolerance)[1]

function Base.eigs(Bcs_in::Operator,A_in::Operator,n::Integer;tolerance::Float64=100eps())
    Bcs, A = promotedomainspace([Bcs_in, A_in])

    nf = size(Bcs,1)
    @assert isfinite(nf)

    typ = promote_type(eltype(Bcs),eltype(A))

    ds=domainspace(A)
    C = Conversion(ds,rangespace(A))

    A1,C1 = zeros(typ,n,n),zeros(typ,n,n)
    A1[1:nf,1:n] = Bcs[1:nf,1:n]
    A1[nf+1:end,1:n] = A[1:n-nf,1:n]
    C1[nf+1:end,1:n] = C[1:n-nf,1:n]

    λ,V = eig(A1,C1)

    λ, V = pruneeigs(λ,V,ds,tolerance)
    p = sortperm(λ; lt=(x,y) -> isless(abs(x),abs(y)))
    λ[p], V[p]
end

function pruneeigs(λ,V,ds,tolerance)
    retλ=eltype(λ)[]
    retV=VFun{typeof(ds),eltype(V)}[]
    n=length(λ)
    for k=1:n
        if slnorm(V,n-3:n,k)≤tolerance
            push!(retλ,λ[k])
            push!(retV,Fun(ds,V[:,k]))
        end
    end
    retλ,retV
end
