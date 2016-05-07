export ToeplitzOperator, HankelOperator, LaurentOperator






type ToeplitzOperator{T<:Number} <: BandedOperator{T}
    negative::Vector{T}
    nonnegative::Vector{T}
end

ToeplitzOperator{T<:Number,Q<:Number}(V::Vector{T},W::Vector{Q})=ToeplitzOperator{promote_type(T,Q)}(V,W)
function ToeplitzOperator(V::Vector)
    W=V[2:end]
    V=copy(V)
    V[1]*=2
    ToeplitzOperator(W,V)
end



getindex(T::ToeplitzOperator,k::Integer,j::Integer) =
    toeplitz_getindex(T.negative,T.nonnegative,k,j)

function toeplitz_getindex{T}(negative::AbstractVector{T},nonnegative::AbstractVector{T},k::Integer,j::Integer)
    if 0<k-j≤length(negative)
        negative[k-j]
    elseif 0≤j-k≤length(nonnegative)-1
        nonnegative[j-k+1]
    else
        zero(T)
    end
end

function toeplitz_getindex{T}(cfs::AbstractVector{T},k::Integer,j::Integer)
    if k==j
        2cfs[1]
    elseif 0<k-j≤length(cfs)-1
        cfs[k-j+1]
    elseif 0<j-k≤length(cfs)-1
        cfs[j-k+1]
    else
        zero(T)
    end
end


bandinds(T::ToeplitzOperator)=(-length(T.negative),length(T.nonnegative)-1)


# slice of a ToeplitzOPerator is a ToeplitzOperator

function Base.getindex(T::ToeplitzOperator,kr::AbstractCount,jr::AbstractCount)
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


type HankelOperator{T<:Number} <: BandedOperator{T}
    coefficients::Vector{T}
end

HankelOperator(f::Fun)=HankelOperator(f.coefficients)


for TYP in (:Operator,:BandedOperator)
    @eval Base.convert{TT}(::Type{$TYP{TT}},T::HankelOperator)=HankelOperator(convert(Vector{TT},T.coefficients))
end

function hankel_getindex(v::AbstractVector,k::Integer,j::Integer)
   if k+j-1 ≤ length(v)
        v[k+j-1]
    else
        zero(eltype(v))
    end
end

getindex(T::HankelOperator,k::Integer,j::Integer) =
    hankel_getindex(T.coefficients,k,j)



bandinds(T::HankelOperator)=(1-length(T.coefficients),length(T.coefficients)-1)



## Laurent Operator

type LaurentOperator{T<:Number} <: BandedOperator{T}
    negative::Vector{T}
    nonnegative::Vector{T}
end

for ATYP in (:Operator,:BandedOperator), TYP in(:ToeplitzOperator,:LaurentOperator)
    @eval Base.convert{TT}(::Type{$ATYP{TT}},T::$TYP)=$TYP(convert(Vector{TT},T.negative),
                                                                            convert(Vector{TT},T.nonnegative))
end




function laurent_getindex{T}(negative::AbstractVector{T},nonnegative::AbstractVector{T},k::Integer,j::Integer)
    # switch to double-infinite indices
    k=iseven(k)?-div(k,2):div(k-1,2)
    j=iseven(j)?-div(j,2):div(j-1,2)

    if 0<k-j≤length(negative)
        negative[k-j]
    elseif 0≤j-k≤length(nonnegative)-1
        nonnegative[j-k+1]
    else
        zero(T)
    end
end





getindex(T::LaurentOperator,k::Integer,j::Integer)=laurent_getindex(T.negative,T.nonnegative,k,j)

shiftbandinds(T::LaurentOperator)=-length(T.negative),length(T.nonnegative)-1
function bandinds(T::LaurentOperator)
    sbi=shiftbandinds(T)
    min(2sbi[1],-2sbi[end]),max(2sbi[end],-2sbi[1])
end


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

for TYP in (:ToeplitzOperator,:LaurentOperator)
    @eval begin
        -(T::$TYP)=$TYP(-T.negative,-T.nonnegative)
        *(c::Number,T::$TYP)=$TYP(c*T.negative,c*T.nonnegative)
    end
end

-(H::HankelOperator)=HankelOperator(-H.coefficients)
*(c::Number,H::HankelOperator)=HankelOperator(c*H.coefficients)


## inv

function Base.inv(T::ToeplitzOperator)
    @assert length(T.nonnegative)==1
    ai=linsolve(T,[1.0];maxlength=100000)
    ToeplitzOperator(ai[2:end],ai[1:1])
end

function Fun(T::ToeplitzOperator)
   if length(T.nonnegative)==1
      Fun([T.nonnegative;T.negative],Taylor())
    elseif length(T.negative)==0
        Fun(T.nonnegative,Hardy{false}())
    else
        Fun(interlace(T.nonnegative,T.negative),Laurent(Circle()))
    end
end
