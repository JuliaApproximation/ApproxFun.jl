# We need to implement some functionality for the ApproxFun constructor to work
real{T}(::Type{Dual{T}})=Dual{ApproxFun.real(T)}

for OP in (:plan_chebyshevtransform,:plan_ichebyshevtransform)
    @eval $OP{D<:Dual}(v::Vector{D})=$OP(realpart(v))
end
chebyshevtransform{D<:Dual}(v::Vector{D},plan...)=dual(chebyshevtransform(realpart(v),plan...),
                                                        chebyshevtransform(dualpart(v),plan...))

#TODO: Hardy{false}
for OP in (:plan_transform,:plan_itransform)
    for TYP in  (:Fourier,:Laurent,:SinSpace)
        @eval $OP{T<:Dual,D<:Domain}(S::$TYP{D},x::Vector{T})=$OP(S,realpart(x))
    end
end

for OP in (:transform,:itransform)
    for TYP in (:Fourier,:Laurent,:SinSpace)
        @eval $OP{T<:Dual,D<:Domain}(S::$TYP{D},x::Vector{T},plan)=dual($OP(S,realpart(x),plan),$OP(S,dualpart(x),plan))
    end
end

chop!(f::Fun,d::Dual)=chop!(f,realpart(d))
samplenorm{T<:Dual}(v::Vector{T})=norm(realpart(v))

