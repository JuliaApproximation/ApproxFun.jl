

## Space promotion for InterlaceOperator
# It's here because we need DirectSumSpace

for TYP in (:PiecewiseSpace,:ArraySpace)
    @eval begin
        function promotedomainspace(A::InterlaceOperator{T,2},sp::$TYP) where T
            if domainspace(A) == sp
                return A
            end
            @assert size(A.ops,2) == length(sp)
            InterlaceOperator([promotedomainspace(A.ops[k,j],sp[j]) for k=1:size(A.ops,1),j=1:size(A.ops,2)],$TYP)
        end
        function interlace_choosedomainspace(ops,rs::$TYP)
            @assert length(ops) == length(rs)
            # this ensures correct dispatch for unino
            sps = Array{Space}(
                filter(x->!isambiguous(x),map((op,s)->choosedomainspace(op,s),ops,rs)))
            if isempty(sps)
                UnsetSpace()
            else
                union(sps...)
            end
        end
    end
end




function continuity(sp::PiecewiseSpace,order::Integer)
    m=ncomponents(sp)
    B=zeros(Operator{prectype(sp)},m-1,m)

    for k=1:m-1
        B[k,k] = Evaluation(component(sp,k),last,order)
        B[k,k+1] = -Evaluation(component(sp,k+1),first,order)
    end

    InterlaceOperator(B,PiecewiseSpace,ArraySpace)
end

function continuity(sp::PiecewiseSpace,kr::UnitRange)
    @assert first(kr)==0
    m=ncomponents(sp)
    B=zeros(Operator{prectype(sp)},length(kr)*(m-1),m)
    for r in kr
        B[(m-1)*r+1:(m-1)*(r+1),:] = continuity(sp,r).ops
    end
    InterlaceOperator(B,PiecewiseSpace,ArraySpace)
end


continuity(d::UnionDomain,k) = continuity(Space(d),k)
continuity(d) = continuity(d,0)

blockdiag(A::PlusOperator) = mapreduce(blockdiag, +, A.ops)
blockdiag(A::TimesOperator) = mapreduce(blockdiag, .*, A.ops)

# TODO: general wrappers

Evaluation(S::SumSpace,x,order) =
    EvaluationWrapper(S,x,order,
        InterlaceOperator(RowVector(vnocat(map(s->Evaluation(s,x,order),components(S))...)),SumSpace))


ToeplitzOperator(G::Fun{MatrixSpace{S,DD,RR},V}) where {S,RR,V,DD} = interlace(map(ToeplitzOperator,Array(G)))

## Sum Space




## Conversion

function coefficients(v::AbstractVector,a::ArraySpace,b::ArraySpace)
    if a==b
        v
    else
        interlace(map((f,s)->Fun(f,s),Fun(a,v),b),b)
    end
end


# ArraySpace is straight forward

function Conversion(a::ArraySpace, b::ArraySpace)
    @assert size(a) == size(b)
    ConversionWrapper(InterlaceOperator(Diagonal(Conversion.(vec(a.spaces), vec(b.spaces))), a, b))
end


# Sum Space and PiecewiseSpace need to allow permutation of space orders
for TYP in (:SumSpace,:PiecewiseSpace)
    @eval function Conversion(S1::$TYP,S2::$TYP)
        v1 = collect(S1.spaces)
        v2 = collect(S2.spaces)

        sort1 = sort(collect(v1))
        sort2 = sort(collect(v2))

        T = promote_type(eltype(domain(S1)),eltype(domain(S2)))

        if any(s->!isinf(dimension(s)),v1) || any(s->!isinf(dimension(s)),v2)
            @assert length(S1.spaces) == length(S2.spaces) == 2
            if hasconversion(S1, S2.spaces[1])
                ops = Operator{T}[Conversion(S1.spaces[1], S2.spaces[1]) Conversion(S1.spaces[2], S2.spaces[1]);
                                ZeroOperator(T, S1.spaces[1], S2.spaces[2]) ZeroOperator(T, S1.spaces[2], S2.spaces[2])]
            elseif hasconversion(S1, S2.spaces[2])
                ops = Operator{T}[ZeroOperator(T, S1.spaces[1], S2.spaces[1]) ZeroOperator(T, S1.spaces[2], S2.spaces[1]);
                                  Conversion(S1.spaces[1], S2.spaces[2]) Conversion(S1.spaces[2], S2.spaces[2])]
            else
                error("Not implemented")
            end
            ConversionWrapper(InterlaceOperator(ops, S1, S2))
        elseif sort1 == sort2
            # swaps sumspace order
            ConversionWrapper(PermutationOperator{T}(perm(v1,v2),S1,S2))
        elseif all(map(hasconversion,v1,v2))
            # we can blocmk convert
            ConversionWrapper(SpaceOperator(
                InterlaceOperator(Diagonal([map(Conversion,v1,v2)...]),$TYP),
                S1,S2))
        elseif all(map(hasconversion,sort1,sort2))
            # we can blocmk convert
            P1 = PermutationOperator{T}(perm(v1,sort1),S1,$TYP(sort1))
            P2 = PermutationOperator{T}(perm(sort2,v2),$TYP(sort2),S2)
            ConversionWrapper(TimesOperator(
                [P2,InterlaceOperator(Diagonal([map(Conversion,sort1,sort2)...]),$TYP),P1]))
        elseif map(canonicalspace,S1.spaces) == map(canonicalspace,S2.spaces)
            error("Not implemented")
        else
            # try sorting canonicalspace
            csort1 = sort(collect(map(canonicalspace,v1)))
            csort2 = sort(collect(map(canonicalspace,v2)))
            if csort1 == csort2
                # we can block convert after permuting
                prm = perm(map(canonicalspace,v1),map(canonicalspace,v2))
                ds2 = $TYP(S1.spaces[perm])
                P = PermutationOperator{T}(prm, S1, ds2)
                ConversionWrapper(TimesOperator(Conversion(ds2,S2),P))
            elseif all(map(hasconversion,csort1,csort2))
                C1 = Conversion(S1,$TYP(csort1))
                C2 = Conversion($TYP(csort1),$TYP(csort2))
                C3 = Conversion($TYP(csort2),S2)

                ConversionWrapper(TimesOperator([C3,C2,C1]))
            else
                # we don't know how to convert so go to default
                defaultConversion(S1,S2)
            end
        end
    end
end

function hasconversion(a::SumSpace, b::Space)
    for n=1:length(a.spaces)
        if hasconversion(a.spaces[n],b) ≠ true
            return false
        end
    end
    return true
end

function Conversion(a::SumSpace, b::Space)
    if !hasconversion(a, b)
        throw(ArgumentError("Cannot convert $a to $b"))
    end

    m=zeros(Operator{promote_type(prectype(a), prectype(b))},1,length(a.spaces))
    for n=1:length(a.spaces)
        m[1,n]=Conversion(a.spaces[n],b)
    end
    return ConversionWrapper(InterlaceOperator(m, a, b, cache(interlacer(a)), cache(BlockInterlacer((repeated(1),))), (1-dimension(b),dimension(a)-1)))
end



for (OPrule,OP) in ((:conversion_rule,:conversion_type),(:maxspace_rule,:maxspace),
                        (:union_rule,:union))
    for TYP in (:SumSpace,:PiecewiseSpace)
        @eval function $OPrule(S1sp::$TYP,S2sp::$TYP)
            S1 = components(S1sp)
            S2 = components(S2sp)
            cs1,cs2=map(canonicalspace,S1),map(canonicalspace,S2)
            if length(S1) != length(S2)
                NoSpace()
            elseif canonicalspace(S1sp) == canonicalspace(S2sp)  # this sorts S1 and S2
                S1sp ≤ S2sp ? S1sp : S2sp  # choose smallest space by sorting
            elseif cs1 == cs2
                # we can just map down
                # $TYP(map($OP,S1.spaces,S2.spaces))
                # this is commented out due to Issue #13261
                newspaces = [$OP(S1[k],S2[k]) for k=1:length(S1)]
                if any(b->b==NoSpace(),newspaces)
                    NoSpace()
                else
                    $TYP(newspaces)
                end
            elseif sort(collect(cs1)) == sort(collect(cs2))
                # sort S1
                p=perm(cs1,cs2)
                $OP($TYP(S1[p]),S2)
            elseif length(S1) == length(S2) == 2  &&
                    $OP(S1[1],S2[1]) != NoSpace() &&
                    $OP(S1[2],S2[2]) != NoSpace()
                #TODO: general length
                $TYP($OP(S1[1],S2[1]),
                     $OP(S1[2],S2[2]))
            elseif length(S1) == length(S2) == 2  &&
                    $OP(S1[1],S2[2])!=NoSpace() &&
                    $OP(S1[2],S2[1])!=NoSpace()
                #TODO: general length
                $TYP($OP(S1[1],S2[2]),
                     $OP(S1[2],S2[1]))
            else
                NoSpace()
            end
        end
    end
end




## Derivative

#TODO: do in @calculus_operator?

for (Op,OpWrap) in ((:Derivative,:DerivativeWrapper),(:Integral,:IntegralWrapper))
    @eval begin
        $Op(S::PiecewiseSpace,k::Integer) =
            $OpWrap(InterlaceOperator(Diagonal([map(s->$Op(s,k),components(S))...]),PiecewiseSpace),k)
        function $Op(S::ArraySpace,k::Integer)
            ops = map(s->$Op(s,k),S)
            $OpWrap(InterlaceOperator(Diagonal(ops),S,ArraySpace(reshape(rangespace.(ops),size(S)))),k)
        end
    end
end

function Derivative(S::SumSpace,k::Integer)
    # we want to map before we decompose, as the map may introduce
    # mixed bases.
    if typeof(canonicaldomain(S))==typeof(domain(S))
        DerivativeWrapper(InterlaceOperator(Diagonal([map(s->Derivative(s,k),components(S))...]),SumSpace),k)
    else
        DefaultDerivative(S,k)
    end
end

choosedomainspace(M::CalculusOperator{UnsetSpace},sp::SumSpace)=mapreduce(s->choosedomainspace(M,s),union,sp.spaces)

## Multiplcation for Array*Vector

function Multiplication(f::Fun{MatrixSpace{S,DD,RR}},sp::VectorSpace{S2,DD2,RR2}) where {S,DD,RR,S2,DD2,RR2}
    @assert size(space(f),2)==length(sp)
    m=Array(f)
    MultiplicationWrapper(f,interlace(Operator{promote_type(cfstype(f),prectype(sp))}[Multiplication(m[k,j],sp[j]) for k=1:size(m,1),j=1:size(m,2)]))
end




## Multiply components

function Multiplication(f::Fun{PW},sp::PiecewiseSpace) where PW<:PiecewiseSpace
    p=perm(domain(f).domains,domain(sp).domains)  # sort f
    vf=components(f)[p]

    MultiplicationWrapper(f,InterlaceOperator(Diagonal([map(Multiplication,vf,sp.spaces)...]),PiecewiseSpace))
end

Multiplication(f::Fun{SumSpace{SV1,D,R1}},sp::SumSpace{SV2,D,R2}) where {SV1,SV2,D,R1,R2} =
    MultiplicationWrapper(f,mapreduce(g->Multiplication(g,sp),+,components(f)))
Multiplication(f::Fun,sp::SumSpace) =
    MultiplicationWrapper(f,InterlaceOperator(Diagonal([map(s->Multiplication(f,s),components(sp))...]),SumSpace))

Multiplication(f::Fun, sp::PiecewiseSpace) = MultiplicationWrapper(f, Multiplication(Fun(f,sp),sp))


# we override coefficienttimes to split the multiplication down to components as union may combine spaes

coefficienttimes(f::Fun{S1},g::Fun{S2}) where {S1<:SumSpace,S2<:SumSpace} = mapreduce(ff->ff*g,+,components(f))
coefficienttimes(f::Fun{S1},g::Fun) where {S1<:SumSpace} = mapreduce(ff->ff*g,+,components(f))
coefficienttimes(f::Fun,g::Fun{S2}) where {S2<:SumSpace} = mapreduce(gg->f*gg,+,components(g))


coefficienttimes(f::Fun{S1},g::Fun{S2}) where {S1<:PiecewiseSpace,S2<:PiecewiseSpace} =
    Fun(map(coefficienttimes,components(f),components(g)),PiecewiseSpace)


## Definite Integral

# This makes sure that the defaults from a given Domain are respected for the UnionDomain.

DefiniteIntegral(d::UnionDomain) =
    DefiniteIntegral(PiecewiseSpace(map(domainspace,map(DefiniteIntegral,d.domains))))
DefiniteLineIntegral(d::UnionDomain) =
    DefiniteLineIntegral(PiecewiseSpace(map(domainspace,map(DefiniteLineIntegral,d.domains))))


DefiniteIntegral(sp::PiecewiseSpace) =
    DefiniteIntegralWrapper(InterlaceOperator(hnocat(map(DefiniteIntegral,sp.spaces)...),sp,ConstantSpace(rangetype(sp))))
DefiniteLineIntegral(sp::PiecewiseSpace) =
    DefiniteLineIntegralWrapper(InterlaceOperator(hnocat(map(DefiniteLineIntegral,sp.spaces)...),sp,ConstantSpace(rangetype(sp))))

## TensorSpace of two PiecewiseSpaces

Base.getindex(d::TensorSpace{Tuple{PWS1,PWS2}},i::Integer,j::Integer) where {PWS1<:PiecewiseSpace,PWS2<:PiecewiseSpace} =
    d[1][i]⊗d[2][j]
Base.getindex(d::TensorSpace{Tuple{PWS1,PWS2}},i::AbstractRange,j::AbstractRange) where {PWS1<:PiecewiseSpace,PWS2<:PiecewiseSpace} =
    PiecewiseSpace(d[1][i])⊗PiecewiseSpace(d[2][j])

## ProductFun

##  Piecewise

function components(U::ProductFun{PS}) where PS<:PiecewiseSpace
    ps=space(U,1)
    sp2=space(U,2)
    m=length(ps)
    C=coefficients(U)
    [ProductFun(C[k:m:end,:],component(ps,k),sp2) for k=1:m]
end
