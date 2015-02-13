Base.show(io::IO,d::Interval)=print(io,"[$(d.a),$(d.b)]")
Base.show(io::IO,d::PeriodicInterval)=print(io,"[$(d.a),$(d.b))")

for typ in ("Chebyshev","Fourier","Laurent")
    TYP=parse(typ)
    @eval function Base.show(io::IO,S::$TYP)
        print(io,$typ*"(")
        show(io,domain(S))
        print(io,")")
    end
end

function Base.show{S}(io::IO,s::ArraySpace{S,1})
    show(io,s.dimensions[1])
    show(io,s.space)
end

function Base.show(io::IO,f::Fun)
    print(io,"Fun(")
    show(io,f.coefficients)
    print(io,",")
    show(io,f.space)
    print(io,")")
end


