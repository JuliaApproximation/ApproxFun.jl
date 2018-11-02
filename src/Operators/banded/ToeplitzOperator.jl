export ToeplitzOperator, HankelOperator


mutable struct ToeplitzOperator{T<:Number} <: Operator{T}
    negative::Vector{T}
    nonnegative::Vector{T}
end


ToeplitzOperator(V::Vector{T},W::Vector{Q}) where {T<:Number,Q<:Number} =
    ToeplitzOperator{promote_type(T,Q)}(V,W)
ToeplitzOperator(V::AbstractVector,W::AbstractVector) =
    ToeplitzOperator(collect(V),collect(W))

convert(::Type{Operator{TT}},T::ToeplitzOperator) where {TT} =
    ToeplitzOperator(convert(Vector{TT},T.negative),convert(Vector{TT},T.nonnegative))

for op in (:(Base.real), :(Base.imag))
    @eval $op(T::ToeplitzOperator) = ToeplitzOperator($op(T.negative), $op(T.nonnegative))
end

function SymToeplitzOperator(V::Vector)
    W=V[2:end]
    V=copy(V)
    V[1]*=2
    ToeplitzOperator(W,V)
end

for OP in (:domainspace,:rangespace)
    @eval $OP(T::ToeplitzOperator) = ℓ⁰
end

getindex(T::ToeplitzOperator,k::Integer,j::Integer) =
    toeplitz_getindex(T.negative,T.nonnegative,k,j)

function toeplitz_getindex(negative::AbstractVector{T},nonnegative::AbstractVector{T},k::Integer,j::Integer) where T
    if 0<k-j≤length(negative)
        negative[k-j]
    elseif 0≤j-k≤length(nonnegative)-1
        nonnegative[j-k+1]
    else
        zero(T)
    end
end

function toeplitz_getindex(cfs::AbstractVector{T},k::Integer,j::Integer) where T
    if k==j && !isempty(cfs)
        2cfs[1]
    elseif 0<k-j≤length(cfs)-1
        cfs[k-j+1]
    elseif 0<j-k≤length(cfs)-1
        cfs[j-k+1]
    else
        zero(T)
    end
end

function BandedMatrix(S::SubOperator{T,ToeplitzOperator{T},Tuple{UnitRange{Int},UnitRange{Int}}}) where T
    ret = BandedMatrix(Zeros, S)

    kr,jr=parentindices(S)

    neg=parent(S).negative
    pos=parent(S).nonnegative

    toeplitz_axpy!(1.0,neg,pos,kr,jr,ret)
end



bandwidths(T::ToeplitzOperator)=(length(T.negative),length(T.nonnegative)-1)


# slice of a ToeplitzOPerator is a ToeplitzOperator

function Base.getindex(T::ToeplitzOperator,kr::InfRanges,jr::InfRanges)
    sh=first(jr)-first(kr)
    st=step(jr)
    @assert st==step(kr)
    if sh ≥0
        ToeplitzOperator([reverse!(T.nonnegative[1:sh]);T.negative],T.nonnegative[sh+1:st:end])
    else
        ToeplitzOperator(T.negative[-sh+1:st:end],[reverse!(T.negative[1:-sh]);T.nonnegative])
    end
end





## Hankel Operator


mutable struct HankelOperator{T<:Number} <: Operator{T}
    coefficients::Vector{T}
end

for OP in (:domainspace,:rangespace)
    @eval $OP(T::HankelOperator) = ℓ⁰
end

HankelOperator(V::AbstractVector)=HankelOperator(collect(V))

HankelOperator(f::Fun)=HankelOperator(f.coefficients)



@eval convert(::Type{Operator{TT}},T::HankelOperator) where {TT}=HankelOperator(convert(Vector{TT},T.coefficients))

function hankel_getindex(v::AbstractVector,k::Integer,j::Integer)
   if k+j-1 ≤ length(v)
        v[k+j-1]
    else
        zero(eltype(v))
    end
end

getindex(T::HankelOperator,k::Integer,j::Integer) =
    hankel_getindex(T.coefficients,k,j)


function BandedMatrix(S::SubOperator{T,HankelOperator{T},Tuple{UnitRange{Int},UnitRange{Int}}}) where T
    ret=BandedMatrix(Zeros, S)

    kr,jr=parentindices(S)
    cfs=parent(S).coefficients

    hankel_axpy!(1.0,cfs,kr,jr,ret)
end


bandwidths(T::HankelOperator) = (max(0,length(T.coefficients)-1),max(0,length(T.coefficients)-1))



## algebra


function Base.maximum(T::ToeplitzOperator)
    if isempty(T.negative)
        maximum(T.nonnegative)
    elseif isempty(T.nonnegative)
        maximum(T.negative)
    else
        max(maximum(T.negative),maximum(T.nonnegative))
    end
end

-(T::ToeplitzOperator)=ToeplitzOperator(-T.negative,-T.nonnegative)
*(c::Number,T::ToeplitzOperator)=ToeplitzOperator(c*T.negative,c*T.nonnegative)

-(H::HankelOperator)=HankelOperator(-H.coefficients)
*(c::Number,H::HankelOperator)=HankelOperator(c*H.coefficients)


## inv

function Base.inv(T::ToeplitzOperator)
    @assert length(T.nonnegative)==1
    ai=\(T,[1.0];maxlength=100000)
    ToeplitzOperator(ai[2:end],ai[1:1])
end

function Fun(T::ToeplitzOperator)
   if length(T.nonnegative)==1
      Fun(Taylor(),[T.nonnegative;T.negative])
    elseif length(T.negative)==0
        Fun(Hardy{false}(),T.nonnegative)
    else
        Fun(Laurent(Circle()),interlace(T.nonnegative,T.negative))
    end
end
