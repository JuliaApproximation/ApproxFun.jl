export ⊕

## SumSpace encodes a space that can be decoupled as f(x) = a(x) + b(x) where a is in S and b is in V


abstract DirectSumSpace{SV,T,DD,d} <: Space{T,DD,d}




for TYP in (:SumSpace,:TupleSpace)
    @eval begin
        immutable $TYP{SV,T,DD,d} <: DirectSumSpace{SV,T,DD,d}
            spaces::SV
            $TYP(dom::Domain)=new(tuple(map(typ->typ(dom),SV.parameters)...))
            $TYP(sp::Tuple)=new(sp)
        end

        $TYP(sp::Tuple)=$TYP{typeof(sp),mapreduce(basistype,promote_type,sp),typeof(domain(first(sp))),ndims(first(sp))}(sp)

        $TYP(A::$TYP,B::$TYP)=$TYP(tuple(A.spaces...,B.spaces...))

        $TYP(A::Space,B::$TYP)=$TYP(tuple(A,B.spaces...))
        $TYP(A::$TYP,B::Space)=$TYP(tuple(A.spaces...,B))
        $TYP(A::Space...)=$TYP(A)
        $TYP(sp::Array)=$TYP(tuple(sp...))

        canonicalspace(A::$TYP)=$TYP(sort([A.spaces...]))

        spacescompatible(A::$TYP,B::$TYP)=length(A.spaces)==length(B.spaces)&&all(map(spacescompatible,A.spaces,B.spaces))

        setdomain(A::$TYP,d::Domain)=$TYP(map(sp->setdomain(sp,d),A.spaces))
    end
end






Base.getindex(S::DirectSumSpace,k)=S.spaces[k]

domain(A::DirectSumSpace)=domain(A.spaces[end])  # TODO: this assumes all spaces have the same domain
                                           #        we use end to avoid ConstantSpace





function spacescompatible(A::Tuple,B::Tuple)
    if length(A) != length(B)
        return false
    end
    #assumes domain doesn't impact sorting
    asort=sort([A...]);bsort=sort([B...])
    for k=1:length(asort)
        if !spacescompatible(asort[k],bsort[k])
            return false
        end
    end

    return true
end


function union_rule(A::SumSpace,B::SumSpace)
    @assert length(A.spaces)==length(B.spaces)==2
    if spacescompatible(A,B)
        A
    elseif spacescompatible(A.spaces,B.spaces)
        A≤B?A:B
    else
        #TODO: should it be attempted to union subspaces?
        SumSpace(union(A.spaces,B.spaces))
    end
end

function union_rule(A::SumSpace,B::Space)
    if B in A.spaces
        A
    else
        A⊕B
    end
end


function union_rule(B::SumSpace,A::ConstantSpace)
    for sp in B.spaces
        if isconvertible(A,sp)
            return B
        end
    end

    NoSpace()
end




## routines

evaluate{D<:SumSpace,T}(f::Fun{D,T},x)=mapreduce(vf->evaluate(vf,x),+,vec(f))
for OP in (:differentiate,:integrate)
    @eval $OP{D<:SumSpace,T}(f::Fun{D,T})=⊕(map($OP,vec(f))...)
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

Base.vec{D<:DirectSumSpace,T}(f::Fun{D,T},k)=Fun(f.coefficients[k:length(space(f).spaces):end],space(f)[k])
Base.vec(S::DirectSumSpace)=S.spaces
Base.vec{S<:DirectSumSpace,T}(f::Fun{S,T})=Fun[vec(f,j) for j=1:length(space(f).spaces)]



## values

itransform(S::SumSpace,cfs)=Fun(cfs,S)(points(S,length(cfs)))


## SumSpace{ConstantSpace}
# this space is special

union_rule{V}(SS::SumSpace{Tuple{ConstantSpace,V}},W::ConstantSpace)=SS
function union_rule{V}(SS::SumSpace{Tuple{ConstantSpace,V}},W::SumSpace)
    a=length(SS.spaces)==2?SS.spaces[2]:$TYP(SS.spaces[2:end])
    if isa(W.spaces[1],ConstantSpace)
        b=length(W.spaces)==2?W.spaces[2]:$TYP(W.spaces[2:end])
        $TYP(SS.spaces[1],union(a,b))
    else
        $TYP(SS.spaces[1],union(SS.spaces[2],W))
    end
end
union_rule{V}(SS::SumSpace{Tuple{ConstantSpace,V}},
                   W::Space)=SumSpace(SS.spaces[1],union(SS.spaces[2],W))

for TYP in (:SumSpace,:TupleSpace)
    @eval begin
        conversion_rule{V,W}(SS::$TYP{Tuple{ConstantSpace,V}},
                             TT::$TYP{Tuple{ConstantSpace,W}})=$TYP(SS.spaces[1],
                                                                             conversion_type(SS.spaces[2],TT.spaces[2]))


        Base.vec{V,TT,DD,d,T}(f::Fun{$TYP{Tuple{ConstantSpace,V},TT,DD,d},T},k)=k==1?Fun(f.coefficients[1],space(f)[1]):Fun(f.coefficients[2:end],space(f)[2])
        Base.vec{V,TT,DD,d,T}(f::Fun{$TYP{Tuple{ConstantSpace,V},TT,DD,d},T})=Any[vec(f,1),vec(f,2)]

        #TODO: fix
        function Base.vec{W,TT,DD,d,T}(f::Fun{$TYP{Tuple{ConstantSpace,ConstantSpace,W},TT,DD,d},T},k)
            if k≤2
                Fun(f.coefficients[k],space(f)[k])
            else
                @assert k==3
                Fun(f.coefficients[k:end],space(f)[k])
            end
        end

        function Base.vec{V,W,TT,DD,d,T}(f::Fun{$TYP{Tuple{ConstantSpace,V,W},TT,DD,d},T},k)
            if k==1
                Fun(f.coefficients[1],space(f)[1])
            else
                Fun(f.coefficients[k:2:end],space(f)[k])
            end
        end


        Base.vec{V,W,TT,DD,d,T}(f::Fun{$TYP{Tuple{ConstantSpace,V,W},TT,DD,d},T})=Any[vec(f,1),vec(f,2),vec(f,3)]
    end
end



#support tuple set
Base.start{SS<:DirectSumSpace}(f::Fun{SS})=start(vec(f))
Base.done{SS<:DirectSumSpace}(f::Fun{SS},k)=done(vec(f),k)
Base.next{SS<:DirectSumSpace}(f::Fun{SS},k)=next(vec(f),k)
