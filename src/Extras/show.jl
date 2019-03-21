## Domains

show(io::IO,d::Segment) = print(io,"the segment [$(leftendpoint(d)),$(rightendpoint(d))]")
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

show(io::IO,::ConstantSpace{AnyDomain}) = print(io,"ConstantSpace")
show(io::IO,S::ConstantSpace) = print(io,"ConstantSpace($(domain(S)))")
show(io::IO,f::Fun{ConstantSpace{AnyDomain}}) =
    print(io,"$(convert(Number,f)) anywhere")

show(io::IO,f::Fun{ConstantSpace{DD,RR}}) where {DD,RR} =
    print(io,"$(convert(Number,f)) on $(domain(f))")

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

function show(io::IO,s::QuotientSpace)
    show(io,s.space)
    print(io," /\n")
    show(io,s.bcs;header=false)
end


function show(io::IO,ss::SumSpace)
    s = components(ss)
    show(io,s[1])
    for sp in s[2:end]
        print(io,"⊕")
        show(io,sp)
    end
end


function show(io::IO,ss::PiecewiseSpace)
    s = components(ss)
    show(io,s[1])
    for sp in s[2:end]
        print(io,"⨄")
        show(io,sp)
    end
end

summary(ss::ArraySpace) = string(Base.dims2string(length.(axes(ss))), " ArraySpace")
function show(io::IO,ss::ArraySpace;header::Bool=true)
    header && print(io,summary(ss)*":\n")
    show(io, ss.spaces)
end

function show(io::IO,s::TensorSpace)
    d = length(s.spaces)
    for i=1:d-1
        show(io,s.spaces[i])
        print(io,"⊗")
    end
    show(io,s.spaces[d])
end

function show(io::IO,s::SubSpace)
    print(io,s.space)
    print(io,"|")
    show(io,s.indexes)
end


## Fun

show(io::IO, ::MIME"text/plain", f::Fun) = show(io, f)

function show(io::IO, f::Fun)
    print(io,"Fun(")
    show(io,f.space)
    print(io,",")
    show(io,f.coefficients)
    print(io,")")
end

## MultivariateFun

show(io::IO, ::MIME"text/plain", f::MultivariateFun) = show(io, f)

function show(io::IO,L::LowRankFun)
    print(io,"LowRankFun on ",space(L)," of rank ",rank(L),".")
end

function show(io::IO,P::ProductFun)
    print(io,"ProductFun on ",space(P),".")
end



## Operator

summary(B::Operator) =
    string(typeof(B).name.name)*" : "*string(domainspace(B))*" → "*string(rangespace(B))

struct PrintShow
    str
end
Base.show(io::IO,N::PrintShow) = print(io,N.str)


function show(io::IO,B::Operator;header::Bool=true)
    header && println(io,summary(B))
    dsp=domainspace(B)

    if !isambiguous(domainspace(B)) && (eltype(B) <: Number)
        if isbanded(B) && isinf(size(B,1)) && isinf(size(B,2))
            BM=B[1:10,1:10]

            M=Matrix{Any}(undef,11,11)
            fill!(M,PrintShow("⋅"))
            for j = 1:size(BM,2),k = colrange(BM,j)
                M[k,j]=BM[k,j]
            end

            for k=max(1,11-bandwidth(B,2)):11
                M[k,end]=PrintShow("⋱")
            end
            for j=max(1,11-bandwidth(B,1)):10
                M[end,j]=PrintShow("⋱")
            end

            print_array(io, M)
        elseif isinf(size(B,1)) && isinf(size(B,2))
            BM=B[1:10,1:10]

            M=Matrix{Any}(undef,11,11)
            for k=1:10,j=1:10
                M[k,j]=BM[k,j]
            end

            M[1,end]=PrintShow("⋯")
            M[end,1]=PrintShow("⋮")

            for k=2:11
                M[k,end]=M[end,k]=PrintShow("⋱")
            end

            print_array(io, M)
        elseif isinf(size(B,1))
            BM=B[1:10,1:size(B,2)]

            M=Matrix{Any}(undef,11,size(B,2))
            for k=1:10,j=1:size(B,2)
                M[k,j]=BM[k,j]
            end
            for k=1:size(B,2)
                M[end,k]=PrintShow("⋮")
            end

            print_array(io, M)
        elseif isinf(size(B,2))
            BM=B[1:size(B,1),1:10]

            M=Matrix{Any}(undef,size(B,1),11)
            for k=1:size(B,1),j=1:10
                M[k,j]=BM[k,j]
            end
            for k=1:size(B,1)
                M[k,end]=PrintShow("⋯")
            end

            print_array(io, M)
        else
            print_array(io, AbstractMatrix(B)[1:size(B,1),1:size(B,2)])
        end
    end
end
