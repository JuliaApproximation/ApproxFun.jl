# We need to implement some functionality for the ApproxFun constructor to work
real{T}(::Type{Dual{T}})=Dual{ApproxFun.real(T)}
plan_chebyshevtransform{D<:Dual}(v::Vector{D})=ApproxFun.plan_chebyshevtransform(realpart(v))
chebyshevtransform{D<:Dual}(v::Vector{D},plan...)=dual(chebyshevtransform(realpart(v),plan...),chebyshevtransform(epsilon(v),plan...))
chop!(f::Fun,d::Dual)=chop!(f,realpart(d))
samplenorm{T<:Dual}(v::Vector{T})=norm(realpart(v))

