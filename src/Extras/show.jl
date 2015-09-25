## Domains

Base.show(io::IO,d::Interval)=print(io,"【$(d.a),$(d.b)】")
function Base.show(io::IO,d::Line)
    if d.center == d.angle == 0 && d.α == d.β == -1.
        print(io,"❪-∞,∞❫")
    elseif  d.α == d.β == -1.
        print(io,"Line($(d.center),$(d.angle))")
    else
        print(io,"Line($(d.center),$(d.angle),$(d.α),$(d.β))")
    end
end
Base.show(io::IO,d::PeriodicInterval)=print(io,"【$(d.a),$(d.b)❫")
Base.show(io::IO,d::Circle)=print(io,(d.radius==1?"":string(d.radius))*"⨀"*(d.center==0?"":"+$(d.center)"))


## Spaces


for typ in ("Chebyshev","Fourier","Laurent")
    TYP=parse(typ)
    @eval function Base.show{D}(io::IO,S::$TYP{D})
        print(io,$typ*"(")
        show(io,domain(S))
        print(io,")")
    end
end



function Base.show(io::IO,s::JacobiWeight)
    d=domain(s)
    #TODO: Get shift and weights right
    if s.α==s.β
        print(io,"(1-x^2)^$(s.α)[")
    else
        print(io,"(1+x)^$(s.α)*(1-x)^$(s.β)[")
    end

    show(io,s.space)
    print(io,"]")
end

function Base.show(io::IO,s::MappedSpace)
    show(io,s.space)
    print(io,"↦")
    show(io,s.domain)
end

function Base.show(io::IO,s::VectorSpace)
    print(io,"[")
    show(io,s.space)
    print(io,"]_"*string(s.dimensions[1]))
end

function Base.show(io::IO,s::MatrixSpace)
    print(io,"[")
    show(io,s.space)
    print(io,"]_("*string(s.dimensions[1])*","*string(s.dimensions[2])*")")
end


function Base.show(io::IO,s::SumSpace)
    show(io,s.spaces[1])
    print(io,"⊕")
    show(io,s.spaces[2])
end

function Base.show(io::IO,s::TensorSpace)
    for i=1:d-1
        show(io,s.spaces[i])
        print(io,"⊗")
    end
    show(io,s.spaces[d])
end




## Fun

function Base.show(io::IO,f::Fun)
    print(io,"Fun(")
    show(io,f.coefficients)
    print(io,",")
    show(io,f.space)
    print(io,")")
end
