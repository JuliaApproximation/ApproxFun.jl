###
# This file contains BLAS/BandedMatrix overrides for operators
# that depend on the structure of BandedMatrix
####



# default copy is to loop through
# override this for most operators.
function default_BandedMatrix(S::Operator)
    Y=BandedMatrix{eltype(S)}(undef, size(S), bandwidths(S))

    for j=1:size(S,2),k=colrange(Y,j)
        @inbounds inbands_setindex!(Y,S[k,j],k,j)
    end

    Y
end


# default copy is to loop through
# override this for most operators.
function default_RaggedMatrix(S::Operator)
    data=Array{eltype(S)}(0)
    cols=Array{Int}(size(S,2)+1)
    cols[1]=1
    for j=1:size(S,2)
        cs=colstop(S,j)
        K=cols[j]-1
        cols[j+1]=cs+cols[j]
        resize!(data,cols[j+1]-1)

        for k=1:cs
            data[K+k]=S[k,j]
        end
    end

    RaggedMatrix(data,cols,size(S,1))
end

function default_Matrix(S::Operator)
    n, m = size(S)
    if isinf(n) || isinf(m)
        error("Cannot convert $S to a Matrix")
    end

    eltype(S)[S[k,j] for k=1:n, j=1:m]
end




# The diagonal of the operator may not be the diagonal of the sub
# banded matrix, so the following calculates the row of the
# Banded matrix corresponding to the diagonal of the original operator


diagindshift(S,kr,jr) = first(kr)-first(jr)
diagindshift(S::SubOperator) = diagindshift(S,parentindices(S)[1],parentindices(S)[2])


#TODO: Remove
diagindrow(S,kr,jr) = bandwidth(S,2)+first(jr)-first(kr)+1
diagindrow(S::SubOperator) = diagindrow(S,parentindices(S)[1],parentindices(S)[2])



#####
# Conversions
#####

function convert(::Type{BandedMatrix},
               S::SubOperator{T,ConcreteConversion{Chebyshev{DD,RR},Ultraspherical{Int,DD,RR},T},
                              Tuple{UnitRange{Int},UnitRange{Int}}}) where {T,DD,RR}
    # we can assume order is 1
    ret = BandedMatrix{eltype(S)}(undef, size(S), bandwidths(S))
    kr,jr = parentindices(S)
    dg = diagindshift(S)

    @assert -bandwidth(ret,1) ≤ dg ≤ bandwidth(ret,2)-2

    ret[band(dg)] = 0.5
    ret[band(dg+1)] = 0.0
    ret[band(dg+2)] = -0.5

    # correct first entry
    if 1 in kr && 1 in jr
        ret[1,1] = 1.0
    end

    ret
end

function convert(::Type{BandedMatrix},V::SubOperator{T,ConcreteConversion{Ultraspherical{LT,DD,RR},Ultraspherical{LT,DD,RR},T},
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
    1-n ≤ dg+1 ≤ m-1 && (ret[band(dg+1)] = 0)
    1-n ≤ dg+2 ≤ m-1 && (ret[band(dg+2)] .= c./(2 .- λ .- jr[max(0,dg+2)+1:min(n+dg+2,m)]))

    ret
end


#####
# Derivatives
#####



function convert(::Type{BandedMatrix},S::SubOperator{T,ConcreteDerivative{Chebyshev{DD,RR},K,T},
                                                     Tuple{UnitRange{Int},UnitRange{Int}}}) where {T,K,DD,RR}

    n,m = size(S)
    ret = BandedMatrix{eltype(S)}(undef, (n,m), bandwidths(S))
    kr,jr = parentindices(S)
    dg = diagindshift(S)

    D = parent(S)
    k = D.order
    d = domain(D)
    C=T(pochhammer(one(T),k-1)/2*(4/(d.b-d.a))^k)

    # need to drop columns


    if 1-n ≤ dg+k ≤ m-1
        ret[band(dg+k)] .= C.*(jr[max(0,dg+k)+1:min(n+dg+k,m)] .- one(T))
    end

    ret
end


function convert(::Type{BandedMatrix},S::SubOperator{T,ConcreteDerivative{Ultraspherical{LT,DD,RR},K,T},
                                                  Tuple{UnitRange{Int},UnitRange{Int}}}) where {T,K,DD,RR,LT}
    n,m = size(S)
    ret = BandedMatrix{eltype(S)}(undef, (n,m), bandwidths(S))
    kr,jr = parentindices(S)
    dg = diagindshift(S)

    D = parent(S)
    k = D.order
    λ = order(domainspace(D))
    d = domain(D)

    C = T(pochhammer(one(T)*λ,k)*(4/(d.b-d.a))^k)
    ret[band(dg+k)] = C

    ret
end


####
# Toeplitz/Hankel
####


# αn,α0,αp give the constant for the negative, diagonal and positive
# entries.  The usual case is 1,2,1
function toeplitz_axpy!(αn,α0,αp,neg,pos,kr,jr,ret)
    dat=ret.data
    m=size(dat,2)

    dg=diagindrow(ret,kr,jr)

    # diagonal
    if dg ≥ 1
        α0p=α0*pos[1]
        @simd for j=1:m
            @inbounds dat[dg,j]+=α0p
        end
    end

    # positive entries
    for k=2:min(length(pos),dg)
        αpp=αp*pos[k]
        @simd for j=1:m
            @inbounds dat[dg-k+1,j]+=αpp
        end
    end

    # negative entries
    for k=1:min(length(neg),size(dat,1)-dg)
        αnn=αn*neg[k]
        @simd for j=1:m
            @inbounds dat[dg+k,j]+=αnn
        end
    end

    ret
end

# this routine is for when we want a symmetric toeplitz
function sym_toeplitz_axpy!(αn,α0,αp,cfs,kr,jr,ret)
    dg=diagindrow(ret,kr,jr)

    dat=ret.data
    m=size(dat,2)
    # diagonal
    if dg ≥ 1
        α0p=α0*cfs[1]
        @simd for j=1:m
            @inbounds dat[dg,j]+=α0p
        end
    end

    # positive entries
    for k=2:min(length(cfs),dg)
        αpp=αp*cfs[k]
        @simd for j=1:m
            @inbounds dat[dg-k+1,j]+=αpp
        end
    end

    # negative entries
    for k=2:min(length(cfs),size(dat,1)-dg+1)
        αnn=αn*cfs[k]
        @simd for j=1:m
            @inbounds dat[dg+k-1,j]+=αnn
        end
    end

    ret
end

toeplitz_axpy!(α,neg,pos,kr,jr,ret) =
    toeplitz_axpy!(α,α,α,neg,pos,kr,jr,ret)

sym_toeplitz_axpy!(α0,αp,cfs,kr,jr,ret) =
    sym_toeplitz_axpy!(αp,α0,αp,cfs,kr,jr,ret)


function hankel_axpy!(α,cfs,kr,jr,ret)
    dat=ret.data

    st=stride(dat,2)-2
    mink=first(kr)+first(jr)-1

    N=length(dat)

    # dg gives the row corresponding to the diagonal of the original operator
    dg=diagindrow(ret,kr,jr)

    # we need the entry where the first entry is written to
    # this is going to be a shift of the diagonal of the true operator
    dg1=dg+first(kr)-first(jr)

    for k=mink:min(length(cfs),size(dat,1)+mink-dg1)
        dk=k-mink+dg1
        nk=k-mink
        αc=α*cfs[k]
        @simd for j=dk:st:min(dk+nk*st,N)
            @inbounds dat[j]+=αc
        end
    end

    ret
end
