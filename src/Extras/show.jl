

## Counts

Base.show(io::IO,c::UnitCount) = print(io,"$(c.start):∞")
Base.show(io::IO,c::Count) = print(io,"$(c.start):$(c.step):∞")


## Domains

Base.show(io::IO,d::Segment) = print(io,"【$(d.a),$(d.b)】")
function Base.show(io::IO,d::Line)
    if d.center == angle(d) == 0 && d.α == d.β == -1.
        print(io,"ℝ")
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

Base.show(io::IO,::ConstantSpace{AnyDomain}) = print(io,"ConstantSpace")
Base.show(io::IO,S::ConstantSpace) = print(io,"ConstantSpace($(domain(S)))")
Base.show(io::IO,f::Fun{ConstantSpace{AnyDomain}}) =
    print(io,"$(Number(f)) anywhere")

Base.show{DD}(io::IO,f::Fun{ConstantSpace{DD}}) =
    print(io,"$(Number(f)) on $(domain(f))")

for typ in ("Chebyshev","Fourier","Laurent","Taylor","SinSpace","CosSpace")
    TYP=parse(typ)
    @eval function Base.show{D}(io::IO,S::$TYP{D})
        print(io,$typ*"(")
        show(io,domain(S))
        print(io,")")
    end
end

function Base.show(io::IO,S::Ultraspherical)
    print(io,"Ultraspherical($(order(S)),")
    show(io,domain(S))
    print(io,")")
end

function Base.show(io::IO,S::Jacobi)
    S.a == S.b == 0 ? print(io,"Legendre(") : print(io,"Jacobi($(S.b),$(S.a),")
    show(io,domain(S))
    print(io,")")
end



function Base.show(io::IO,s::JacobiWeight)
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



function Base.show(io::IO,s::SumSpace)
    show(io,s[1])
    for sp in s[2:end]
        print(io,"⊕")
        show(io,sp)
    end
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

function Base.show(io::IO,s::SubSpace)
    print(io,s.space)
    print(io,"|")
    show(io,s.indexes)
end


## Fun

function Base.show(io::IO,f::Fun)
    print(io,"Fun(")
    show(io,f.space)
    print(io,",")
    show(io,f.coefficients)
    print(io,")")
end

## MultivariateFun

function Base.show(io::IO,L::LowRankFun)
    print(io,"LowRankFun on ",space(L)," of rank ",rank(L),".")
end

function Base.show(io::IO,P::ProductFun)
    print(io,"ProductFun on ",space(P),".")
end



## Operator

Base.summary(B::Operator) = string(typeof(B).name.name)*":"*string(domainspace(B))*"→"*string(rangespace(B))


function Base.show(io::IO,B::Operator;header::Bool=true)
    header && println(io,summary(B))
    dsp=domainspace(B)

    if !isambiguous(domainspace(B)) && eltype(B) <: Number
        if isbanded(B) && isinf(size(B,1)) && isinf(size(B,2))
            BM=B[1:10,1:10]

            M=Matrix{Any}(11,11)
            fill!(M,PrintShow(""))
            for j = 1:size(BM,2),k = colrange(BM,j)
                M[k,j]=BM[k,j]
            end

            for k=max(1,11-bandinds(B,2)):11
                M[k,end]=PrintShow("⋱")
            end
            for j=max(1,11+bandinds(B,1)):10
                M[end,j]=PrintShow("⋱")
            end

            showarray(io,M;header=false)
        elseif isinf(size(B,1)) && isinf(size(B,2))
            BM=B[1:10,1:10]

            M=Matrix{Any}(11,11)
            for k=1:10,j=1:10
                M[k,j]=BM[k,j]
            end

            M[1,end]=PrintShow("⋯")
            M[end,1]=PrintShow("⋮")

            for k=2:11
                M[k,end]=M[end,k]=PrintShow("⋱")
            end

            showarray(io,M;header=false)
        elseif isinf(size(B,1))
            BM=B[1:10,1:size(B,2)]

            M=Array(Any,11,size(B,2))
            for k=1:10,j=1:size(B,2)
                M[k,j]=BM[k,j]
            end
            for k=1:size(B,2)
                M[end,k]=PrintShow("⋮")
            end

            showarray(io,M;header=false)
        elseif isinf(size(B,2))
            BM=B[1:size(B,1),1:10]

            M=Matrix{Any}(size(B,1),11)
            for k=1:size(B,1),j=1:10
                M[k,j]=BM[k,j]
            end
            for k=1:size(B,1)
                M[k,end]=PrintShow("⋯")
            end

            showarray(io,M;header=false)
        else
            showarray(io,AbstractMatrix(B)[1:size(B,1),1:size(B,2)];header=false)
        end
    end
end


function Base.show{T<:Operator}(io::IO, ::MIME"text/plain", A::Vector{T}; header::Bool=true)
    nf = length(A)-1
    header && for k=1:nf+1 println(io,summary(A[k])) end
    if all(Ak -> isafunctional(Ak), A[1:nf]) && isbanded(A[end]) &&
            isinf(size(A[end],1)) && isinf(size(A[end],2)) && eltype(A[end]) <: Number &&
            all(Ak -> !isambiguous(domainspace(Ak)), A)
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

        for j = 1:size(BM,2), k = colrange(BM,j)
            MM[k,j]=BM[k,j]
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
end
