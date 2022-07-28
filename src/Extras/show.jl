## Domains

function show(io::IO,d::Line)
    if d.center == angle(d) == 0 && d.α == d.β == -1.
        print(io,"ℝ")
    elseif  d.α == d.β == -1.
        print(io,"Line($(d.center),$(angle(d)))")
    else
        print(io,"Line($(d.center),$(angle(d)),$(d.α),$(d.β))")
    end
end

function show(io::IO,d::Ray)
    if d.orientation && angle(d)==0
        print(io,"【$(d.center),∞❫")
    elseif  d.orientation && angle(d)==1.0π
        print(io,"【$(d.center),-∞❫")
    elseif  d.orientation
        print(io,"【$(d.center),exp($(angle(d))im)∞❫")
    elseif !d.orientation  && angle(d)==0
        print(io,"❪∞,$(d.center)】")
    elseif !d.orientation && angle(d)==1.0π
        print(io,"❪-∞,$(d.center)】")
    else
        print(io,"❪exp($(angle(d))im)∞,$(d.center)】")
    end
end

show(io::IO,d::PeriodicSegment) = print(io,"【$(leftendpoint(d)),$(rightendpoint(d))❫")
show(io::IO,d::Circle) =
    print(io,(d.radius==1 ? "" : string(d.radius))*
                    (d.orientation ? "🕒" : "🕞")*
                    (d.center==0 ? "" : "+$(d.center)"))
show(io::IO,d::Point) = print(io,"Point($(d.x))")


## Spaces

for typ in ("Chebyshev","Fourier","Laurent","Taylor","SinSpace","CosSpace")
    TYP = Meta.parse(typ)
    @eval function show(io::IO,S::$TYP{D,R}) where {D,R}
        print(io,$typ*"(")
        show(io,domain(S))
        print(io,")")
    end
end



function show(io::IO, S::Ultraspherical)
    print(io,"Ultraspherical($(order(S)),")
    show(io,domain(S))
    print(io,")")
end

function show(io::IO,S::Jacobi)
    S.a == S.b == 0 ? print(io,"Legendre(") : print(io,"Jacobi($(S.b),$(S.a),")
    show(io,domain(S))
    print(io,")")
end


show(io::IO, S::Chebyshev{<:ChebyshevInterval}) = print(io, "Chebyshev()")
show(io::IO, S::Ultraspherical{<:Any,<:ChebyshevInterval}) = print(io, "Ultraspherical($(order(S)))")
show(io::IO, S::Jacobi{<:ChebyshevInterval}) =
    S.a == S.b == 0 ? print(io,"Legendre()") : print(io,"Jacobi($(S.b),$(S.a))")

show(io::IO, S::NormalizedPolynomialSpace) = (print(io, "Normalized"); show(io, S.space))

function show(io::IO,s::JacobiWeight)
    d=domain(s)
    #TODO: Get shift and weights right
    if s.α==s.β
        print(io,"(1-x^2)^$(s.α)[")
    elseif s.β==0
        print(io,"(1-x)^$(s.α)[")
    elseif s.α==0
        print(io,"(1+x)^$(s.β)[")
    else
        print(io,"(1+x)^$(s.β)*(1-x)^$(s.α)[")
    end

    show(io,s.space)
    print(io,"]")
end

function show(io::IO,s::LogWeight)
    d=domain(s)
    #TODO: Get shift and weights right
    if s.α==s.β
        print(io,"log((1-x^2)^$(s.α))[")
    else
        print(io,"log((1+x)^$(s.α)*(1-x)^$(s.β))[")
    end

    show(io,s.space)
    print(io,"]")
end
