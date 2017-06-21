export continuity



Space(d::IntervalDomain) = Chebyshev(d)

identity_fun{T<:Number}(d::Segment{T}) = Fun(Chebyshev(d),T[(d.b+d.a)/2,(d.b-d.a)/2])


## Calculus



# the default domain space is higher to avoid negative ultraspherical spaces
Integral(d::IntervalDomain,n::Integer) = Integral(Ultraspherical(1,d),n)

for Func in (:DefiniteIntegral,:DefiniteLineIntegral)
    @eval begin
        #TODO: this may be misleading
        $Func(d::IntervalDomain) = $Func(JacobiWeight(-.5,-.5,Chebyshev(d)))
        function $Func(α::Number,β::Number,d::IntervalDomain)
            @assert α == β
            @assert round(Int,α+.5) == α+.5
            @assert round(Int,α+.5) >= 0
            $Func(JacobiWeight(α,β,Ultraspherical(round(Int,α+.5),d)))
        end
        $Func(α::Number,β::Number) = $Func(α,β,Interval())
    end
end
