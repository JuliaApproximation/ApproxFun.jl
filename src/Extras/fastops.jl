###
# This file contains BLAS/BandedMatrix overrides for operators
# that depend on the structure of BandedMatrix
####




#####
# Conversions
#####

function BandedMatrix(S::SubOperator{T,ConcreteConversion{Chebyshev{DD,RR},Ultraspherical{Int,DD,RR},T},
                              Tuple{UnitRange{Int},UnitRange{Int}}}) where {T,DD,RR}
    # we can assume order is 1
    ret = BandedMatrix{eltype(S)}(undef, size(S), bandwidths(S))
    kr,jr = parentindices(S)
    dg = diagindshift(S)

    @assert -bandwidth(ret,1) ≤ dg ≤ bandwidth(ret,2)-2

    ret[band(dg)] .= 0.5
    ret[band(dg+1)] .= 0.0
    ret[band(dg+2)] .= -0.5

    # correct first entry
    if 1 in kr && 1 in jr
        ret[1,1] = 1.0
    end

    ret
end

function BandedMatrix(V::SubOperator{T,ConcreteConversion{Ultraspherical{LT,DD,RR},Ultraspherical{LT,DD,RR},T},
                                                                  Tuple{UnitRange{Int},UnitRange{Int}}}) where {T,LT,DD,RR}

    n,m = size(V)
    V_l, V_u = bandwidths(V)
    ret = BandedMatrix{eltype(V)}(undef, (n,m), (V_l,V_u))
    kr,jr = parentindices(V)
    dg = diagindshift(V)


    λ = order(rangespace(parent(V)))
    c = λ-one(T)

    # need to drop columns



    1-n ≤ dg ≤ m-1 && (ret[band(dg)] .= c./(jr[max(0,dg)+1:min(n+dg,m)] .- 2 .+ λ))
    1-n ≤ dg+1 ≤ m-1 && (ret[band(dg+1)] .= 0)
    1-n ≤ dg+2 ≤ m-1 && (ret[band(dg+2)] .= c./(2 .- λ .- jr[max(0,dg+2)+1:min(n+dg+2,m)]))

    ret
end


#####
# Derivatives
#####



function BandedMatrix(S::SubOperator{T,ConcreteDerivative{Chebyshev{DD,RR},K,T},
                                                     Tuple{UnitRange{Int},UnitRange{Int}}}) where {T,K,DD,RR}

    n,m = size(S)
    ret = BandedMatrix{eltype(S)}(undef, (n,m), bandwidths(S))
    kr,jr = parentindices(S)
    dg = diagindshift(S)

    D = parent(S)
    k = D.order
    d = domain(D)

    C=convert(T,pochhammer(one(T),k-1)/2*(4/(complexlength(d)))^k)


    # need to drop columns


    if 1-n ≤ dg+k ≤ m-1
        ret[band(dg+k)] .= C.*(jr[max(0,dg+k)+1:min(n+dg+k,m)] .- one(T))
    end

    ret
end


function BandedMatrix(S::SubOperator{T,ConcreteDerivative{Ultraspherical{LT,DD,RR},K,T},
                                                  Tuple{UnitRange{Int},UnitRange{Int}}}) where {T,K,DD,RR,LT}
    n,m = size(S)
    ret = BandedMatrix{eltype(S)}(undef, (n,m), bandwidths(S))
    kr,jr = parentindices(S)
    dg = diagindshift(S)

    D = parent(S)
    k = D.order
    λ = order(domainspace(D))
    d = domain(D)

    C = convert(T,pochhammer(one(T)*λ,k)*(4/(complexlength(d)))^k)
    ret[band(dg+k)] .= C

    ret
end

## These are both hacks that apparently work

function BandedMatrix(S::SubOperator{T,PP,Tuple{UnitRange{Int},UnitRange{Int}}}) where {T,PP<:PartialInverseOperator}
    kr,jr = parentindices(S)
    P = parent(S)
    #ret = BandedMatrix{eltype(S)}(undef, size(S), bandwidths(S))
    ret = BandedMatrix{eltype(S)}(undef, (last(kr),last(jr)), bandwidths(P))
    b = bandwidth(P, 2)
    #@assert first(kr) == first(jr) == 1

    @inbounds for j in 1:last(jr)
        kk = colrange(ret, j)
        if j in kk
            ret[j,j] = inv(P.cache[j,j])
        end
        for k in first(kk):min(last(kk),j-1)
            t = zero(T)
            for i = max(k,j-b-1):j-1
                t += ret[k,i]*P.cache[i,j]
            end
            ret[k,j] = -t/P.cache[j,j]
        end
    end

    ret[kr,jr]
end

function BandedMatrix(S::SubOperator{T,ConcreteConversion{QuotientSpace{SP,O,D,R},SP,T},Tuple{UnitRange{Int},UnitRange{Int}}}) where {SP,O,D,R,T}
    kr,jr = parentindices(S)
    C = parent(S)
    #ret = BandedMatrix{eltype(S)}(undef, size(S), bandwidths(S))
    ret = BandedMatrix{eltype(S)}(undef, (last(kr),last(jr)), bandwidths(C))
    #@assert first(kr) == first(jr) == 1

    sp = domainspace(C)
    F = sp.F
    A = F.factors
    n = size(A, 1)
    B = sp.bcs[1:n,1:last(jr)+n]
    x = sp.x
    @inbounds for j in 1:last(jr)
        kk = colrange(ret, j)
        if j in kk
            ret[j,j] = one(R)
        end
        for jj = 1:n, ii = 1:n
            A[ii,jj] = B[ii,j+jj]
        end
        for ii = 1:n
            x[ii] = -B[ii,j]
        end
        if norm(x) > 8*norm(A)*eps(R)
            mutable_lu!(F)
            ldiv!(F, x)
        end
        for k = first(kk)+1:last(kk)
            ret[k,j] = x[k-j]
        end
    end

    ret[kr,jr]
end
