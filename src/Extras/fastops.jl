###
# This file contains BLAS/Base.copy overrides for operators
# that depend on the structure of BandedMatrix
####



# default copy is to loop through
# override this for most operators.
function default_copy(S::SubBandedMatrix)
    l,u = bandwidth(S,1),bandwidth(S,2)
    Y=BandedMatrix(eltype(S),size(S,1),size(S,2),l,u)

    for (k,j) in eachbandedindex(S)
        @inbounds Y.data[k-j+u+1,j]=S[k,j]
    end

    Y
end

Base.copy(S::SubBandedMatrix) = default_copy(S)




# The diagonal of the operator may not be the diagonal of the sub
# banded matrix, so the following calculates the row of the
# Banded matrix corresponding to the diagonal of the original operator

function diagindrow(S::SubBandedMatrix)
    u=bandwidth(S,2)
    kr,jr=parentindexes(S)
    u+first(jr)-first(kr)+1
end

#####
# Conversions
#####

function Base.copy{T,C<:Chebyshev,U<:Ultraspherical{1}}(S::SubBandedMatrix{T,ConcreteConversion{C,U,T},
                                                                              Tuple{UnitRange{Int},UnitRange{Int}}})
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

function Base.copy{T,m,λ,DD}(S::SubBandedMatrix{T,ConcreteConversion{Ultraspherical{m,DD},Ultraspherical{λ,DD},T},
                                                                              Tuple{UnitRange{Int},UnitRange{Int}}})
    ret=bzeros(S)

    kr,jr=parentindexes(S)
    dat=ret.data

    shft = jr[1]-1

    dg1=diagindrow(S)
    dg3=dg1-2

    c=λ-one(T)  # this supports big types
    @simd for j=1:size(dat,2)
        @inbounds dat[dg3,j]=c/(j+shft - 2 + λ)
    end

    @simd for j=max(3,3-shft):size(dat,2)
        @inbounds dat[dg1,j]=-c/(j+shft - 2 + λ)
    end

    ret
end


#####
# Derivatives
#####



function Base.copy{T,K,DD,λ}(S::SubBandedMatrix{T,ConcreteDerivative{Ultraspherical{λ,DD},K,T},
                                                                Tuple{UnitRange{Int},UnitRange{Int}}})
    D=parent(S)
    m=D.order
    d=domain(D)

    ret=bzeros(S)
    u=bandwidth(ret,2)

    kr,jr=parentindexes(S)
    dat=ret.data

    shft = jr[1]-1

    dg=diagindrow(S)-m   # mth superdiagonal

    if λ == 0
        C=(.5pochhammer(one(T),m-1)*(4./(d.b-d.a)).^m)::T
        @simd for j=max(m+1-shft,1):size(dat,2)
            @inbounds dat[dg,j]=C*(j+shft-one(T))
        end
    else
        C=(pochhammer(one(T)*λ,m)*(4./(d.b-d.a)).^m)::T
        @simd for j=max(m+1-shft,1):size(dat,2)
            @inbounds dat[dg,j]=C
        end
    end

    ret
end


####
# Toeplitz/Hankel
####



####
# Multiplication
####
