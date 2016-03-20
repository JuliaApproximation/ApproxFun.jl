## Domains

Base.show(io::IO,d::Interval)=print(io,"【$(d.a),$(d.b)】")
function Base.show(io::IO,d::Line)
    if d.center == angle(d) == 0 && d.α == d.β == -1.
        print(io,"❪-∞,∞❫")
    elseif  d.α == d.β == -1.
        print(io,"Line($(d.center),$(angle(d)))")
    else
        print(io,"Line($(d.center),$(angle(d)),$(d.α),$(d.β))")
    end
end

function Base.show(io::IO,d::Ray)
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

Base.show(io::IO,d::PeriodicInterval)=print(io,"【$(d.a),$(d.b)❫")
Base.show(io::IO,d::Circle)=print(io,(d.radius==1?"":string(d.radius))*(d.orientation?"🕒":"🕞")*(d.center==0?"":"+$(d.center)"))
Base.show(io::IO,d::Point)=print(io,"Point($(d.x))")


function Base.show(io::IO,s::UnionDomain)
    show(io,s[1])
    for d in s[2:end]
        print(io,"∪")
        show(io,d)
    end
end

## Spaces

Base.show(io::IO,::AnySpace)=print(io,"AnySpace")

Base.show(io::IO,::ConstantSpace{AnyDomain})=print(io,"ConstantSpace")
Base.show(io::IO,S::ConstantSpace)=print(io,"ConstantSpace($(domain(S)))")


for typ in ("Chebyshev","Fourier","Laurent")
    TYP=parse(typ)
    @eval function Base.show{D}(io::IO,S::$TYP{D})
        print(io,$typ*"(")
        show(io,domain(S))
        print(io,")")
    end
end

function Base.show{λ,D}(io::IO,S::Ultraspherical{λ,D})
    print(io,"Ultraspherical{$λ}(")
    show(io,domain(S))
    print(io,")")
end

function Base.show(io::IO,S::Jacobi)
    print(io,"Jacobi($(S.a),$(S.b),")
    show(io,domain(S))
    print(io,")")
end



function Base.show(io::IO,s::JacobiWeight)
    d=domain(s)
    #TODO: Get shift and weights right
    if s.α==s.β
        print(io,"(1-x^2)^$(s.α)[")
    elseif s.α==0
        print(io,"(1-x)^$(s.β)[")
    elseif s.β==0
        print(io,"(1+x)^$(s.α)[")
    else
        print(io,"(1+x)^$(s.α)*(1-x)^$(s.β)[")
    end

    show(io,s.space)
    print(io,"]")
end

function Base.show(io::IO,s::LogWeight)
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
    show(io,s[1])
    for sp in s[2:end]
        print(io,"⊕")
        show(io,sp)
    end
end

function Base.show(io::IO,s::TupleSpace)
    print(io,"⟨")
    show(io,s[1])
    for sp in s[2:end]
        print(io,",")
        show(io,sp)
    end
    print(io,"⟩")
end

function Base.show(io::IO,s::PiecewiseSpace)
    show(io,s[1])
    for sp in s[2:end]
        print(io,"⨄")
        show(io,sp)
    end
end

function Base.show(io::IO,s::TensorSpace)
    d = length(s)
    for i=1:d-1
        show(io,s[i])
        print(io,"⊗")
    end
    show(io,s[d])
end




## Fun

function Base.show(io::IO,f::Fun)
    print(io,"Fun(")
    show(io,f.coefficients)
    print(io,",")
    show(io,f.space)
    print(io,")")
end

## MultivariateFun

function Base.show(io::IO,L::LowRankFun)
    print(io,"LowRankFun on ",space(L)," of rank ",rank(L),".")
end

function Base.show(io::IO,P::ProductFun)
    print(io,"ProductFun on ",space(P))
end



## Operator

Base.summary(B::Operator)=string(typeof(B).name.name)*":"*string(domainspace(B))*"↦"*string(rangespace(B))


function Base.show(io::IO,B::BandedOperator;header::Bool=true)
    header && println(io,summary(B))

    BM=B[1:10,1:10]

    M=Array(Any,11,11)
    fill!(M,PrintShow(""))
    for kj=eachbandedindex(BM)
        M[kj]=BM[kj]
    end

    for k=max(1,11-bandinds(B,2)):11
        M[k,end]=PrintShow("⋱")
    end
    for j=max(1,11+bandinds(B,1)):10
        M[end,j]=PrintShow("⋱")
    end

    Base.showarray(io,M;header=false)
end

function Base.show(io::IO,F::Functional;header::Bool=true)
    header && println(io,summary(F))

    V=Array{Any}(1,11)
    copy!(V,F[1:10])
    V[end]=PrintShow("⋯")

    Base.showarray(io,V;header=false)
end

function Base.writemime{T<:Functional}(io::IO, ::MIME"text/plain", A::Vector{T};header::Bool=true)
    n = length(A)
    header && for k=1:n println(io,summary(A[k])) end
    M=Array{Any}(n,11)
    for k=1:n
        M[k,1:10] = A[k][1:10]
        M[k,end]=PrintShow("⋯")
    end
    Base.with_output_limit(()->Base.print_matrix(io, M))
end

function Base.writemime{T<:Operator}(io::IO, ::MIME"text/plain", A::Vector{T};header::Bool=true)
    nf = length(A)-1
    for k=1:nf @assert isa(A[k],Functional) end
    header && for k=1:nf+1 println(io,summary(A[k])) end
    M=Array{Any}(11,11)
    fill!(M,PrintShow(""))
    for k=1:nf
        M[k,1:10] = A[k][1:10]
        M[k,end]=PrintShow("⋯")
    end

    MM=Array{Any}(11-nf,11)
    fill!(MM,PrintShow(""))

    B = A[end]
    BM=B[1:10-nf,1:10]

    for kj=eachbandedindex(BM)
        MM[kj]=BM[kj]
    end

    for k=1+nf:10,j=1:10
        M[k,j] = MM[k-nf,j]
    end

    for k=max(1,11-bandinds(B,2)+nf):11
        M[k,end]=PrintShow("⋱")
    end
    for j=max(1,11+bandinds(B,1)-nf):10
        M[end,j]=PrintShow("⋱")
    end

    Base.with_output_limit(()->Base.print_matrix(io, M))
end
