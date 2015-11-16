## CurveSpace


Space{S<:Fourier}(d::PeriodicCurve{S})=Fourier(d)
Space{S<:Laurent}(d::PeriodicCurve{S})=Laurent(d)

#TODO: Make type stable
Base.convert(::Type{Curve},f::Fun)=isa(domain(f),IntervalDomain)?IntervalCurve(f):PeriodicCurve(f)

# function evaluate{C<:Curve,TT}(f::AbstractVector,S::Space{TT,C},x::Number)
#     rts=roots(domain(S).curve-x)
#     @assert length(rts)==1
#     evaluate(Fun(f,setdomain(S,canonicaldomain(S))),first(rts))
# end



identity_fun{C<:Curve,TT}(d::Space{TT,C})=Fun(domain(d).curve.coefficients,
                                              setdomain(space(d.domain.curve),domain(d)))
