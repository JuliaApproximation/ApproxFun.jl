

for op in (:dirichlet,:neumann,:continuity,:ivp)
    @eval $op(d::PiecewiseSpace,k...)=interlace($op(d.spaces,k...))
    @eval $op(d::UnionDomain,k...)=interlace($op(d.domains,k...))
end



## diag provides a way to convert between DiagonalInterlaceOperator and bacn
# function blkdiagm{B<:Operator}(v::Vector{B})
#     if spacescompatible(map(domainspace,v)) && spacescompatible(map(rangespace,v))
#         # arrayspace
#         DiagonalInterlaceOperator(v)
#     elseif domainscompatible(map(domainspace,v)) && domainscompatible(map(rangespace,v))
#         # tuplespace
#         DiagonalInterlaceOperator(v)
#     else
#         DiagonalPiecewiseOperator(v)
#     end
# end

Base.blkdiag(A::DiagonalInterlaceOperator)=A.ops
Base.blkdiag(A::PlusOperator)=mapreduce(blkdiag,+,A.ops)
Base.blkdiag(A::TimesOperator)=mapreduce(blkdiag,.*,A.ops)

# TODO: general wrappers
for TYP in (:DerivativeWrapper,:ConversionWrapper)
    @eval Base.blkdiag{DT<:DiagonalInterlaceOperator}(A::($TYP{DT}))=A.op.ops
end

Base.blkdiag{FT<:PiecewiseSpace,OT<:DiagonalInterlaceOperator}(A::MultiplicationWrapper{FT,OT})=A.op.ops


## Vector



immutable DiagonalArrayOperator{B<:BandedOperator,T<:Number} <: BandedOperator{T}
    op::B
    dimensions::Tuple{Vararg{Int}}
end

DiagonalArrayOperator{T}(op::BandedOperator{T},dms::Tuple{Vararg{Int}})=DiagonalArrayOperator{typeof(op),T}(op,dms)
#DiagonalArrayOperator{T}(op::BandedOperator{T},dms::Int)=DiagonalArrayOperator(op,(dms,))


function bandinds(D::DiagonalArrayOperator)
    bra,brb=bandinds(D.op)
    n=*(D.dimensions...)
    n*bra,n*brb
end


function addentries!(D::DiagonalArrayOperator,A,kr::Range,::Colon)
    n=*(D.dimensions...)
    for k=1:n
        stride_addentries!(D.op,k-n,k-n,n,n,A,kr)
    end
    A
end


for op in (:domainspace,:rangespace)
    @eval $op(D::DiagonalArrayOperator)=ArraySpace($op(D.op),D.dimensions)
end


Derivative(AS::ArraySpace,k::Integer)=DerivativeWrapper(DiagonalArrayOperator(Derivative(AS.space,k),size(AS)),k)

function conversion_rule(AS::ArraySpace,BS::ArraySpace)
    if size(AS)==size(BS)
        ArraySpace(conversion_type(AS.space,BS.space),size(AS))
    else
        NoSpace()
    end
end

for OP in (:maxspace,:conversion_type)
    @eval function $OP(AS::ArraySpace,BS::ArraySpace)
        if size(AS)==size(BS)
            ArraySpace($OP(AS.space,BS.space),size(AS))
        else
            NoSpace()
        end
    end
end


function Conversion(AS::ArraySpace,BS::ArraySpace)
    @assert size(AS)==size(BS)
    ConversionWrapper(DiagonalArrayOperator(Conversion(AS.space,BS.space),size(AS)))
end

ToeplitzOperator{S,T,V,DD}(G::Fun{MatrixSpace{S,T,DD,1},V})=interlace(map(ToeplitzOperator,mat(G)))
ToeplitzOperator{S,T,V,DD}(G::Fun{MatrixSpace{S,T,DD,1},V})=interlace(map(ToeplitzOperator,mat(G)))

## Sum Space




## Conversion


# TupleSpace maps down
function Conversion(a::TupleSpace,b::TupleSpace)
    m=findlast(s->dimension(s)==1,a.spaces)
    @assert all(s->dimension(s)==1,a.spaces[1:m-1]) &&
            a.spaces[1:m-1]==b.spaces[1:m-1]

    if m==0
        ConversionWrapper(DiagonalInterlaceOperator(map(Conversion,a.spaces,b.spaces),TupleSpace))
    elseif length(a)==m
        SpaceOperator(FiniteOperator(eye(m)),a,b)
    elseif length(a)==m+1
        ConversionWrapper(BlockOperator(eye(m),zeros(m,0),zeros(0,m),Conversion(a[end],b[end])))
    else
        ConversionWrapper(BlockOperator(eye(m),zeros(m,0),zeros(0,m),Conversion(TupleSpace(a[m+1:end]),TupleSpace(b[m+1:end]))))
    end
end





# Sum Space and PiecewiseSpace need to allow permutation of space orders
for TYP in (:SumSpace,:PiecewiseSpace)
    @eval function Conversion(S1::$TYP,S2::$TYP)
        if any(s->dimension(s)!=Inf,S1.spaces) || any(s->dimension(s)!=Inf,S2.spaces)
            error("Need to implement finite dimensional case")
        elseif sort([S1.spaces...])==sort([S2.spaces...])
            # swaps sumspace order
            ConversionWrapper(SpaceOperator(
            PermutationOperator(promote_type(eltype(domain(S1)),eltype(domain(S2))),S1.spaces,S2.spaces),
                          S1,S2))
        elseif all(map(hasconversion,S1.spaces,S2.spaces))
            # we can blocmk convert
            ConversionWrapper(DiagonalInterlaceOperator(map(Conversion,S1.spaces,S2.spaces),$TYP))
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
            TimesOperator(Conversion(ds2,S2),Conversion(S1,ds2))

        else
            # we don't know how to convert so go to default
            defaultconversion(S1,S2)
        end
    end
end




for (OPrule,OP) in ((:conversion_rule,:conversion_type),(:maxspace_rule,:maxspace))
    for TYP in (:SumSpace,:PiecewiseSpace)
        @eval function $OPrule(S1::$TYP,S2::$TYP)
            cs1,cs2=map(canonicalspace,S1.spaces),map(canonicalspace,S2.spaces)
            if length(S1.spaces)!=length(S2.spaces)
                NoSpace()
            elseif canonicalspace(S1)==canonicalspace(S2)  # this sorts S1 and S2
                S1 ≤ S2?S1:S2  # choose smallest space by sorting
            elseif cs1==cs2
                # we can just map down
                # $TYP(map($OP,S1.spaces,S2.spaces))
                # this is commented out due to Issue #13261
                newspaces=[$OP(S1[k],S2[k]) for k=1:length(S1.spaces)]
                if any(b->b==NoSpace(),newspaces)
                    NoSpace()
                else
                    $TYP(newspaces)
                end
            elseif sort([cs1...])== sort([cs2...])
                # sort S1
                p=perm(cs1,cs2)
                $OP($TYP(S1.spaces[p]),S2)
            elseif length(S1)==length(S2)==2  &&
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

    # TupleSpace doesn't allow reordering
    @eval function $OPrule(S1::TupleSpace,S2::TupleSpace)
        K=findfirst(s->!isa(s,ConstantSpace),S1)-1

        if length(S1)==length(S2)  &&
                all(s->isa(s,ConstantSpace),S2[1:K]) &&
                all(s->!isa(s,ConstantSpace),S2[K+1:end])
            newspaces=[$OP(S1[k],S2[k]) for k=K+1:length(S1)]

            if any(b->b==NoSpace(),newspaces)
                NoSpace()
            else
                TupleSpace(tuple(S1[1:K]...,newspaces...))
            end
        else
            NoSpace()
        end
    end

end




## Derivative

#TODO: do in @calculus_operator?

for TYP in (:SumSpace,:PiecewiseSpace,:TupleSpace),(Op,OpWrap) in ((:Derivative,:DerivativeWrapper),
                                                                   (:Integral,:IntegralWrapper))
    @eval $Op(S::$TYP,k)=$OpWrap(DiagonalInterlaceOperator(map(s->$Op(s,k),S.spaces),$TYP),k)
end




## Multiplcation for Array*Vector

function Multiplication{S,T,DD}(f::Fun{MatrixSpace{S,T,DD,1}},sp::VectorSpace)
    @assert size(space(f),2)==length(sp)
    m=mat(f)
    MultiplicationWrapper(f,interlace(BandedOperator{promote_type(eltype(f),eltype(sp))}[Multiplication(m[k,j],sp.space) for k=1:size(m,1),j=1:size(m,2)]))
end




## Multiply pieces

function Multiplication{PW<:PiecewiseSpace}(f::Fun{PW},sp::PiecewiseSpace)
    vf=vec(f)
    @assert length(vf)==length(sp)
    MultiplicationWrapper(f,DiagonalInterlaceOperator(map(Multiplication,vf,sp.spaces),PiecewiseSpace))
end


function bandinds{S<:SumSpace,SS<:SumSpace}(M::Multiplication{S,SS})
    a,b=vec(M.f)
    sp=domainspace(M)
    bandinds(Multiplication(a,sp)+Multiplication(b,sp))
end
function rangespace{S<:SumSpace,SS<:SumSpace}(M::Multiplication{S,SS})
    a,b=vec(M.f)
    sp=domainspace(M)
    rangespace(Multiplication(a,sp)+Multiplication(b,sp))
end
function addentries!{S<:SumSpace,SS<:SumSpace}(M::Multiplication{S,SS},A,k,::Colon)
    a,b=vec(M.f)
    sp=domainspace(M)
    addentries!(Multiplication(a,sp)+Multiplication(b,sp),A,k,:)
end




function bandinds{S,SS<:SumSpace}(M::Multiplication{S,SS})
    a,b=vec(domainspace(M))
    Ma=Multiplication(M.f,a)
    Mb=Multiplication(M.f,b)

    bandinds(DiagonalInterlaceOperator((Ma,Mb),SumSpace))
end
function rangespace{S,SS<:SumSpace}(M::Multiplication{S,SS})
    a,b=vec(domainspace(M))
    Ma=Multiplication(M.f,a)
    Mb=Multiplication(M.f,b)

    rangespace(Ma)⊕rangespace(Mb)
end
function addentries!{S,SS<:SumSpace}(M::Multiplication{S,SS},A,k,::Colon)
    a,b=vec(domainspace(M))
    Ma=Multiplication(M.f,a)
    Mb=Multiplication(M.f,b)

    addentries!(DiagonalInterlaceOperator((Ma,Mb),SumSpace),A,k,:)
end




## Definite Integral

# This makes sure that the defaults from a given Domain are respected for the UnionDomain.

DefiniteIntegral(d::UnionDomain) = DefiniteIntegral(PiecewiseSpace(map(domainspace,map(DefiniteIntegral,d.domains))))
DefiniteLineIntegral(d::UnionDomain) = DefiniteLineIntegral(PiecewiseSpace(map(domainspace,map(DefiniteLineIntegral,d.domains))))

####### This is a hack to get the Faraday Cage working.
function getindex{PWS<:PiecewiseSpace,T}(Σ::DefiniteLineIntegral{PWS,T},kr::Range)
    d = domain(Σ)
    n = length(d)
    promote_type(T,eltype(d))[k ≤ n? one(T) : zero(T) for k=kr]
end
datalength{PWS<:PiecewiseSpace,T}(Σ::DefiniteLineIntegral{PWS,T})=length(domain(Σ))
####### This is a hack to get the Faraday Cage working.

## TensorSpace of two PiecewiseSpaces

Base.getindex{PWS1<:PiecewiseSpace,PWS2<:PiecewiseSpace}(d::TensorSpace{Tuple{PWS1,PWS2}},i::Integer,j::Integer)=d[1][i]⊗d[2][j]
Base.getindex{PWS1<:PiecewiseSpace,PWS2<:PiecewiseSpace}(d::TensorSpace{Tuple{PWS1,PWS2}},i::Range,j::Range)=PiecewiseSpace(d[1][i])⊗PiecewiseSpace(d[2][j])

## ProductFun

##  Piecewise

function pieces{PS<:PiecewiseSpace}(U::ProductFun{PS})
    ps=space(U,1)
    sp2=space(U,2)
    m=length(ps)
    C=coefficients(U)
    [ProductFun(C[k:m:end,:],ps[k],sp2) for k=1:m]
end
