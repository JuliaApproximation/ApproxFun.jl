export continuity



Space(d::IntervalOrSegment) = Chebyshev(d)
Space(d::FullSpace{<:Real}) = Chebyshev(Line())

Fun(::typeof(identity), d::Segment{T}) where {T<:Number} =
    Fun(Chebyshev(d), T[(rightendpoint(d)+leftendpoint(d))/2, complexlength(d)/2])


## Calculus



# the default domain space is higher to avoid negative ultraspherical spaces
Integral(d::IntervalOrSegment,n::Integer) = Integral(Ultraspherical(1,d),n)

for Func in (:DefiniteIntegral,:DefiniteLineIntegral)
    @eval begin
        #TODO: this may be misleading
        $Func(d::IntervalOrSegment) = $Func(JacobiWeight(-0.5,-0.5,Chebyshev(d)))
        function $Func(α::Number,β::Number,d::IntervalOrSegment)
            @assert α == β
            @assert round(Int,α+.5) == α+.5
            @assert round(Int,α+.5) >= 0
            $Func(JacobiWeight(α,β,Ultraspherical(round(Int,α+.5),d)))
        end
        $Func(α::Number,β::Number) = $Func(α,β,ChebyshevInterval())
    end
end
