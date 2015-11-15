## CurveSpace

Space{S<:IntervalDomain}(d::Curve{S})=Chebyshev(d)
Space{S<:PeriodicDomain}(d::Curve{S})=Fourier(d)

# function evaluate{C<:Curve,TT}(f::AbstractVector,S::Space{TT,C},x::Number)
#     rts=roots(domain(S).curve-x)
#     @assert length(rts)==1
#     evaluate(Fun(f,setdomain(S,canonicaldomain(S))),first(rts))
# end


identity_fun{C<:Curve,TT}(d::Space{TT,C})=Fun(d.domain.curve.coefficients,d)
