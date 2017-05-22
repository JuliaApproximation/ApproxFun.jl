

## Space promotion for InterlaceOperator
# It's here because we need DirectSumSpace

for TYP in (:PiecewiseSpace,:ArraySpace)
    @eval begin
        function promotedomainspace{T}(A::InterlaceOperator{T,2},sp::$TYP)
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







for op in (:dirichlet,:neumann,:continuity,:ivp)
    @eval $op(d::PiecewiseSpace,k...) = InterlaceOperator($op(d.spaces,k...),PiecewiseSpace,ArraySpace)
    @eval $op(d::UnionDomain,k...) = InterlaceOperator($op(d.domains,k...),PiecewiseSpace,ArraySpace)
end



Base.blkdiag(A::PlusOperator)=mapreduce(blkdiag,+,A.ops)
Base.blkdiag(A::TimesOperator)=mapreduce(blkdiag,.*,A.ops)

# TODO: general wrappers

Evaluation(S::SumSpace,x,order) =
    EvaluationWrapper(S,x,order,
        InterlaceOperator(hcat(map(s->Evaluation(s,x,order),components(S))...),SumSpace))


ToeplitzOperator{S,T,V,DD}(G::Fun{MatrixSpace{S,T,DD,1},V}) = interlace(map(ToeplitzOperator,Array(G)))

## Sum Space




## Conversion

function coefficients(v::AbstractVector,a::ArraySpace,b::ArraySpace)
    if a==b
        v
    else
        interlace(map((f,s)->Fun(f,s),Fun(a,v),b),b)
    end
end


# Sum Space and PiecewiseSpace need to allow permutation of space orders
for TYP in (:SumSpace,:PiecewiseSpace)
    @eval function Conversion(S1::$TYP,S2::$TYP)
        if any(s->!isinf(dimension(s)),S1.spaces) || any(s->!isinf(dimension(s)),S2.spaces)
            error("Need to implement finite dimensional case")
        elseif sort([S1.spaces...])==sort([S2.spaces...])
            # swaps sumspace order
            ConversionWrapper(SpaceOperator(
            PermutationOperator(promote_type(eltype(domain(S1)),eltype(domain(S2))),S1.spaces,S2.spaces),
                          S1,S2))
        elseif all(map(hasconversion,S1.spaces,S2.spaces))
            # we can blocmk convert
            ConversionWrapper(InterlaceOperator(Diagonal([map(Conversion,S1.spaces,S2.spaces)...]),$TYP))
        elseif map(canonicalspace,S1.spaces)==map(canonicalspace,S2.spaces)
            error("Not implemented")
        elseif sort([map(canonicalspace,S1.spaces)...])==sort([map(canonicalspace,S2.spaces)...])
            # we can block convert after permuting
            P=PermutationOperator(promote_type(eltype(domain(S1)),eltype(domain(S2))),
                                  map(canonicalspace,S1.spaces),
                                  map(canonicalspace,S2.spaces))
            ds2=$TYP(S1.spaces[P.perm])
            ConversionWrapper(TimesOperator(Conversion(ds2,S2),SpaceOperator(P,S1,ds2)))
        elseif all(map(hasconversion,reverse(S1.spaces),S2.spaces))
            # special case that comes up, especially for two spaces
            rS1=SumSpace(reverse(S1.spaces))
            ConversionWrapper(TimesOperator(
                Conversion(rS1,S2),
                SpaceOperator(PermutationOperator(length(S1.spaces):-1:1),
                            S1,rS1)))
        elseif all(map(hasconversion,sort([map(canonicalspace,S1.spaces)...]),sort([map(canonicalspace,S2.spaces)...])))
            #TODO: general case
            @assert length(S1.spaces)==2
            ds2=$TYP(S1.spaces[[2,1]])
            ConversionWrapper(TimesOperator(Conversion(ds2,S2),Conversion(S1,ds2)))

        else
            # we don't know how to convert so go to default
            defaultConversion(S1,S2)
        end
    end
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
            elseif sort([cs1...]) == sort([cs2...])
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

for TYP in (:PiecewiseSpace,:ArraySpace),(Op,OpWrap) in ((:Derivative,:DerivativeWrapper),
                                                         (:Integral,:IntegralWrapper))
    @eval $Op(S::$TYP,k::Integer) =
        $OpWrap(InterlaceOperator(Diagonal([map(s->$Op(s,k),components(S))...]),$TYP),k)
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

function Multiplication{S,T,DD,S2,T2,DD2,dim}(f::Fun{MatrixSpace{S,T,DD,dim}},sp::VectorSpace{S2,T2,DD2,dim})
    @assert size(space(f),2)==length(sp)
    m=Array(f)
    MultiplicationWrapper(f,interlace(Operator{promote_type(eltype(f),eltype(sp))}[Multiplication(m[k,j],sp[j]) for k=1:size(m,1),j=1:size(m,2)]))
end




## Multiply components

function Multiplication{PW<:PiecewiseSpace}(f::Fun{PW},sp::PiecewiseSpace)
    p=perm(domain(f).domains,domain(sp).domains)  # sort f
    vf=components(f)[p]

    MultiplicationWrapper(f,InterlaceOperator(Diagonal([map(Multiplication,vf,sp.spaces)...]),PiecewiseSpace))
end

Multiplication{SV1,SV2,T2,T1,D,d}(f::Fun{SumSpace{SV1,T1,D,d}},sp::SumSpace{SV2,T2,D,d}) =
    MultiplicationWrapper(f,mapreduce(g->Multiplication(g,sp),+,components(f)))
Multiplication(f::Fun,sp::SumSpace) =
    MultiplicationWrapper(f,InterlaceOperator(Diagonal([map(s->Multiplication(f,s),components(sp))...]),SumSpace))


# we override coefficienttimes to split the multiplication down to components as union may combine spaes

coefficienttimes{S1<:SumSpace,S2<:SumSpace}(f::Fun{S1},g::Fun{S2}) = mapreduce(ff->ff*g,+,components(f))
coefficienttimes{S1<:SumSpace}(f::Fun{S1},g::Fun) = mapreduce(ff->ff*g,+,components(f))
coefficienttimes{S2<:SumSpace}(f::Fun,g::Fun{S2}) = mapreduce(gg->f*gg,+,components(g))


coefficienttimes{S1<:PiecewiseSpace,S2<:PiecewiseSpace}(f::Fun{S1},g::Fun{S2})=depiece(map(coefficienttimes,components(f),components(g)))


## Definite Integral

# This makes sure that the defaults from a given Domain are respected for the UnionDomain.

DefiniteIntegral(d::UnionDomain) =
    DefiniteIntegral(PiecewiseSpace(map(domainspace,map(DefiniteIntegral,d.domains))))
DefiniteLineIntegral(d::UnionDomain) =
    DefiniteLineIntegral(PiecewiseSpace(map(domainspace,map(DefiniteLineIntegral,d.domains))))


DefiniteIntegral(sp::PiecewiseSpace) =
    DefiniteIntegralWrapper(InterlaceOperator(hcat(map(DefiniteIntegral,sp.spaces)...),sp,ConstantSpace()))
DefiniteLineIntegral(sp::PiecewiseSpace) =
    DefiniteLineIntegralWrapper(InterlaceOperator(hcat(map(DefiniteLineIntegral,sp.spaces)...),sp,ConstantSpace()))

## TensorSpace of two PiecewiseSpaces

Base.getindex{PWS1<:PiecewiseSpace,PWS2<:PiecewiseSpace}(d::TensorSpace{Tuple{PWS1,PWS2}},i::Integer,j::Integer) =
    d[1][i]⊗d[2][j]
Base.getindex{PWS1<:PiecewiseSpace,PWS2<:PiecewiseSpace}(d::TensorSpace{Tuple{PWS1,PWS2}},i::Range,j::Range) =
    PiecewiseSpace(d[1][i])⊗PiecewiseSpace(d[2][j])

## ProductFun

##  Piecewise

function components{PS<:PiecewiseSpace}(U::ProductFun{PS})
    ps=space(U,1)
    sp2=space(U,2)
    m=length(ps)
    C=coefficients(U)
    [ProductFun(C[k:m:end,:],component(ps,k),sp2) for k=1:m]
end
