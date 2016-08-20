###
# This file contains BLAS/BandedMatrix overrides for operators
# that depend on the structure of BandedMatrix
####



# default copy is to loop through
# override this for most operators.
function default_bandedmatrix(S::Operator)
    l,u = bandwidth(S,1),bandwidth(S,2)
    Y=BandedMatrix(eltype(S),size(S,1),size(S,2),l,u)

    for j=1:size(S,2),k=colrange(Y,j)
        @inbounds Y.data[k-j+u+1,j]=S[k,j]
    end

    Y
end




# The diagonal of the operator may not be the diagonal of the sub
# banded matrix, so the following calculates the row of the
# Banded matrix corresponding to the diagonal of the original operator


diagindrow(S,kr,jr) =
    bandwidth(S,2)+first(jr)-first(kr)+1


diagindrow(S::SubOperator) =
    diagindrow(S,parentindexes(S)[1],parentindexes(S)[2])



#####
# Conversions
#####

function Base.convert{T,C<:Chebyshev,U<:Ultraspherical{1}}(::Type{BandedMatrix},
                        S::SubOperator{T,ConcreteConversion{C,U,T},Tuple{UnitRange{Int},UnitRange{Int}}})
    ret=bzeros(S)
    kr,jr=parentindexes(S)


    dat=ret.data
    dg1=diagindrow(S)
    dg3=dg1-2

    @assert 1 ≤ dg1 ≤ size(dat,1)
    @assert 1 ≤ dg3 ≤ size(dat,1)

    @simd for j=1:size(dat,2)
        @inbounds dat[dg1,j]=0.5
        @inbounds dat[dg3,j]=-0.5
    end

    if 1 in kr && 1 in jr
        dat[dg1,1]=1.0
    end

    ret
end

function Base.convert{T,m,λ,DD}(::Type{BandedMatrix},S::SubOperator{T,ConcreteConversion{Ultraspherical{m,DD},Ultraspherical{λ,DD},T},
                                                                              Tuple{UnitRange{Int},UnitRange{Int}}})
    ret=bzeros(S)

    kr,jr=parentindexes(S)
    dat=ret.data

    shft = first(jr)-1

    dg1=diagindrow(S)
    dg3=dg1-2

    c=λ-one(T)  # this supports big types
    @simd for j=1:size(dat,2)
        @inbounds dat[dg1,j]=c/(j+shft - 2 + λ)
    end

    @simd for j=max(1,3-shft):size(dat,2)
        @inbounds dat[dg3,j]=-c/(j+shft - 2 + λ)
    end

    ret
end


#####
# Derivatives
#####



function Base.convert{T,K,DD}(::Type{BandedMatrix},S::SubOperator{T,ConcreteDerivative{Chebyshev{DD},K,T},
                                                                Tuple{UnitRange{Int},UnitRange{Int}}})
    D=parent(S)
    m=D.order
    d=domain(D)

    ret=bzeros(S)
    u=bandwidth(ret,2)

    kr,jr=parentindexes(S)
    dat=ret.data

    shft = first(jr)-1

    dg=diagindrow(S)-m   # mth superdiagonal

    C=(.5pochhammer(one(T),m-1)*(4./(d.b-d.a)).^m)::T
    @simd for j=max(m+1-shft,1):size(dat,2)
        @inbounds dat[dg,j]=C*(j+shft-one(T))
    end

    ret
end


function Base.convert{T,K,DD,λ}(::Type{BandedMatrix},S::SubOperator{T,ConcreteDerivative{Ultraspherical{λ,DD},K,T},
                                                                Tuple{UnitRange{Int},UnitRange{Int}}})
    D=parent(S)
    m=D.order
    d=domain(D)

    ret=bzeros(S)
    u=bandwidth(ret,2)

    kr,jr=parentindexes(S)
    dat=ret.data

    shft = first(jr)-1

    dg=diagindrow(S)-m   # mth superdiagonal

    C=(pochhammer(one(T)*λ,m)*(4./(d.b-d.a)).^m)::T
    @simd for j=max(m+1-shft,1):size(dat,2)
        @inbounds dat[dg,j]=C
    end

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
