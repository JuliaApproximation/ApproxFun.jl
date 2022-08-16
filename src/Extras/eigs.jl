


eigvals(A::Operator,n::Integer;tolerance::Float64=100eps()) =
    eigs(A,n;tolerance=tolerance)[1]

"""
    λ, V = eigs(A::Operator, n::Integer; tolerance::Float64=100eps())

Compute the eigenvalues and eigenvectors of the operator `A`. This is done in the following way:

* Truncate `A` into an n×n matrix `A₁`.
* Compute eigenvalues and eigenvectors of `A₁`.
* Filter out those eigenvectors of `A₁`, which are approximately eigenvectors
of `A` as well. The `tolerance` argument controls, which eigenvectors of the approximation are kept.
"""
function eigs(A::Operator,n::Integer;tolerance::Float64=100eps())
    typ = eltype(A)

    ds=domainspace(A)
    C = Conversion(ds,rangespace(A))

    A1,C1 = zeros(typ,n,n),zeros(typ,n,n)
    A1[1:end,1:n] = A[1:n,1:n]
    C1[1:end,1:n] = C[1:n,1:n]

    λ,V=eigen(A1,C1)

    pruneeigs(λ,V,ds,tolerance)
end

eigvals(Bcs::Operator,A::Operator,n::Integer;tolerance::Float64=100eps()) =
    eigs(Bcs,A,n;tolerance=tolerance)[1]

"""
    λ, V = eigs(BC::Operator, A::Operator, n::Integer; tolerance::Float64=100eps())

Compute `n` eigenvalues and eigenvectors of the operator `A`,
subject to the boundary conditions `BC`.

# Examples
```jldoctest
julia> #= We compute the spectrum of the second derivative,
          subject to zero boundary conditions.
          We solve this eigenvalue problem in the Chebyshev basis =#

julia> S = Chebyshev();

julia> D = Derivative(S, 2);

julia> BC = Dirichlet(S);

julia> λ, v = ApproxFun.eigs(BC, D, 100);

julia> λ[1:10] ≈ [-(n*pi/2)^2 for n in 1:10] # compare with the analytical result
true
```
"""
function eigs(Bcs_in::Operator,A_in::Operator,n::Integer;tolerance::Float64=100eps())
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

    λ,V = eigen(A1,C1)

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
