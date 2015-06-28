export ⊕

## SumSpace encodes a space that can be decoupled as f(x) = a(x) + b(x) where a is in S and b is in V

if VERSION≥v"0.4.0-dev"
    immutable SumSpace{SV,T,d} <: FunctionSpace{T,d}
        spaces::SV
        SumSpace(dom::Domain)=new(tuple(map(typ->typ(dom),SV.parameters)...))
        SumSpace(sp::Tuple)=new(sp)
    end
else
    immutable SumSpace{SV,T,d} <: FunctionSpace{T,d}
        spaces::SV
        SumSpace(dom::Domain)=new(map(typ->typ(dom),SV))
        SumSpace(sp::Tuple)=new(sp)
    end
end


SumSpace(sp::Tuple)=SumSpace{typeof(sp),mapreduce(basistype,promote_type,sp),ndims(first(sp))}(sp)
SumSpace(A::SumSpace,B::SumSpace)=SumSpace((A.spaces...,B.spaces...))
SumSpace(A::FunctionSpace,B::SumSpace)=SumSpace((A.spaces...,B.spaces...))
SumSpace(A::SumSpace,B::FunctionSpace)=SumSpace((A.spaces...,B))
SumSpace(A::FunctionSpace,B::FunctionSpace)=SumSpace((A,B))
SumSpace(sp::Array)=SumSpace(tuple(sp...))

canonicalspace(A::SumSpace)=SumSpace(sort([A.spaces...]))

⊕(A::FunctionSpace,B::FunctionSpace)=SumSpace(A,B)
⊕(f::Fun,g::Fun)=Fun(interlace(coefficients(f),coefficients(g)),space(f)⊕space(g))





Base.getindex(S::SumSpace,k)=S.spaces[k]

domain(A::SumSpace)=domain(A[1])  # TODO: this assumes all spaces have the same domain
setdomain(A::SumSpace,d::Domain)=SumSpace(map(sp->setdomain(sp,d),A.spaces))


spacescompatible(A::SumSpace,B::SumSpace)=all(map(spacescompatible,A.spaces,B.spaces))




function union_rule(A::SumSpace,B::FunctionSpace)
    if B in A.spaces
        A
    else
        error("Implement")
    end
end
function union_rule(A::SumSpace,B::SumSpace)
    @assert length(A.spaces)==length(B.spaces)==2
    if spacescompatible(A,B)
        A
    elseif spacescompatible(A[1],B[2]) && spacescompatible(B[1],A[2])
        A
    else
        error("Implement")
    end
end


function union_rule(B::SumSpace,A::ConstantSpace)
    for sp in B.spaces
        if union(A,sp)==sp
            return B
        end
    end

    NoSpace()
end


coefficients(cfs::Vector,A::SumSpace,B::SumSpace)=defaultcoefficients(cfs,A,B)



function coefficients(cfs::Vector,A::FunctionSpace,B::SumSpace)
    @assert length(B.spaces)==2
    if spacescompatible(A,B.spaces[1])
        interlace(cfs,[zero(eltype(cfs))])
    elseif spacescompatible(A,B.spaces[2])
        interlace([zero(eltype(cfs))],cfs)
    else
        defaultcoefficients(cfs,A,B)
    end
end




## routines

evaluate{D<:SumSpace,T}(f::Fun{D,T},x)=mapreduce(vf->evaluate(vf,x),+,vec(f))
for OP in (:differentiate,:integrate)
    @eval $OP{D<:SumSpace,T}(f::Fun{D,T})=$OP(vec(f,1))⊕$OP(vec(f,2))
end

# assume first domain has 1 as a basis element

function Base.ones(S::SumSpace)
    if union(ConstantSpace(),S.spaces[1])==S.spaces[1]
        ones(S[1])⊕zeros(S[2])
    else
        zeros(S[1])⊕ones(S[2])
    end
end

function Base.ones{T<:Number}(::Type{T},S::SumSpace)
    @assert length(S.spaces)==2
    if union(ConstantSpace(),S.spaces[1])==S.spaces[1]
        ones(T,S[1])⊕zeros(T,S[2])
    else
        zeros(T,S[1])⊕ones(T,S[2])
    end
end


# vec

Base.vec{D<:SumSpace,T}(f::Fun{D,T},k)=Fun(f.coefficients[k:length(space(f).spaces):end],space(f)[k])
Base.vec(S::SumSpace)=S.spaces
Base.vec{S<:SumSpace,T}(f::Fun{S,T})=Fun[vec(f,j) for j=1:length(space(f).spaces)]



## values

itransform(S::SumSpace,cfs)=Fun(cfs,S)[points(S,length(cfs))]


## SumSpace{ConstantSpace}
# this space is special


conversion_rule{V<:FunctionSpace}(SS::SumSpace{@compat(Tuple{ConstantSpace,V})},::V)=SS
Base.vec{V,TT,d,T}(f::Fun{SumSpace{@compat(Tuple{ConstantSpace,V}),TT,d},T},k)=k==1?Fun(f.coefficients[1],space(f)[1]):Fun(f.coefficients[2:end],space(f)[2])
Base.vec{V,TT,d,T}(f::Fun{SumSpace{@compat(Tuple{ConstantSpace,V}),TT,d},T})=Any[vec(f,1),vec(f,2)]

