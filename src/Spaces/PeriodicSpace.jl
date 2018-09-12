Fun(::typeof(identity), d::Circle) = Fun(Laurent(d),[d.center,0.,d.radius])

Space(d::PeriodicDomain) = Fourier(d)
Space(d::Circle) = Laurent(d)


## Evaluation

Evaluation(d::PeriodicDomain,x::Number,n...) = Evaluation(Laurent(d),complex(x),n...)

## Definite Integral

DefiniteIntegral(d::PeriodicDomain) = DefiniteIntegral(Laurent(d))
DefiniteLineIntegral(d::PeriodicDomain) = DefiniteLineIntegral(Laurent(d))

## Toeplitz
union_rule(A::Space{<:PeriodicInterval}, B::Space{<:IntervalOrSegment}) =
    union(Space(Interval(domain(A))), B)
