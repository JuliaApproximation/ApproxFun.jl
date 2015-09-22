export sumblkdiagm



immutable DiagonalPiecewiseOperator{T<:Number,B<:Operator} <: AbstractDiagonalInterlaceOperator{T,B}
    ops::Vector{B}
end

DiagonalPiecewiseOperator{B<:Operator}(v::Vector{B})=DiagonalPiecewiseOperator{mapreduce(eltype,promote_type,v),B}(v)
DiagonalPiecewiseOperator{TT<:Tuple}(v::Union{Vector{Any},TT})=DiagonalPiecewiseOperator(Operator{mapreduce(eltype,promote_type,v)}[v...;])



for op in (:domainspace,:rangespace)
    @eval $op(D::DiagonalPiecewiseOperator)=PiecewiseSpace(map($op,D.ops))
end




Derivative(sp::PiecewiseSpace)=DerivativeWrapper(DiagonalPiecewiseOperator(map(Derivative,sp.spaces)),1)
Derivative(sp::PiecewiseSpace,k::Integer)=DerivativeWrapper(DiagonalPiecewiseOperator(map(s->Derivative(s,k),sp.spaces)),k)


function Multiplication{PW<:PiecewiseSpace}(f::Fun{PW},sp::PiecewiseSpace)
    vf=vec(f)
    @assert length(vf)==length(sp)
    MultiplicationWrapper(f,DiagonalPiecewiseOperator(BandedOperator{promote_type(eltype(f),eltype(sp))}[Multiplication(vf[k],sp[k]) for k=1:length(vf)]))
end

function Conversion(f::PiecewiseSpace,g::PiecewiseSpace)
    @assert length(f)==length(g)
    ConversionWrapper(DiagonalPiecewiseOperator(BandedOperator{Float64}[Conversion(f[k],g[k]) for k=1:length(f)]))
end

for op in (:dirichlet,:neumann,:continuity,:ivp)
    @eval $op(d::PiecewiseSpace,k...)=interlace($op(d.spaces,k...))
    @eval $op(d::UnionDomain,k...)=interlace($op(d.domains,k...))
end



## diag provides a way to convert between DiagonalInterlaceOperator and bacn
function blkdiagm{B<:Operator}(v::Vector{B})
    if spacescompatible(map(domainspace,v)) && spacescompatible(map(rangespace,v))
        DiagonalInterlaceOperator(v)
    else
        DiagonalPiecewiseOperator(v)
    end
end

Base.blkdiag(A::AbstractDiagonalInterlaceOperator)=A.ops
Base.blkdiag(A::PlusOperator)=mapreduce(blkdiag,+,A.ops)
Base.blkdiag(A::TimesOperator)=mapreduce(blkdiag,.*,A.ops)

# TODO: general wrappers
for TYP in (:DerivativeWrapper,:ConversionWrapper)
    @eval Base.blkdiag{DT<:AbstractDiagonalInterlaceOperator}(A::($TYP{DT}))=A.op.ops
end

Base.blkdiag{FT<:PiecewiseSpace,OT<:AbstractDiagonalInterlaceOperator}(A::MultiplicationWrapper{FT,OT})=A.op.ops


## Vector



immutable DiagonalArrayOperator{B<:BandedOperator,T<:Number} <: BandedOperator{T}
    op::B
    dimensions::@compat(Tuple{Vararg{Int}})
end

DiagonalArrayOperator{T}(op::BandedOperator{T},dms::@compat(Tuple{Vararg{Int}}))=DiagonalArrayOperator{typeof(op),T}(op,dms)
#DiagonalArrayOperator{T}(op::BandedOperator{T},dms::Int)=DiagonalArrayOperator(op,(dms,))


function bandinds(D::DiagonalArrayOperator)
    bra,brb=bandinds(D.op)
    n=*(D.dimensions...)
    n*bra,n*brb
end


function addentries!(D::DiagonalArrayOperator,A,kr::Range)
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


immutable SumInterlaceOperator{T,B<:Operator} <: AbstractDiagonalInterlaceOperator{T,B}
    ops::Vector{B}
end

SumInterlaceOperator{B<:Operator}(v::Vector{B})=SumInterlaceOperator{mapreduce(eltype,promote_type,v),B}(v)
SumInterlaceOperator(v::Vector{Any})=SumInterlaceOperator(Operator{mapreduce(eltype,promote_type,v)}[v...])


Base.convert{BT<:SumInterlaceOperator}(::Type{BT},SI::BT)=SI
Base.convert{BT<:Operator}(::Type{BT},SI::SumInterlaceOperator)=SumInterlaceOperator(BandedOperator{eltype(BT)}[SI.ops...])

for op in (:domainspace,:rangespace)
    @eval $op(S::SumInterlaceOperator)=SumSpace($op(S.ops[1]),$op(S.ops[2]))
end

sumblkdiagm{B<:Operator}(v::Vector{B})=SumInterlaceOperator(v)



## Conversion



function Conversion(S1::SumSpace,S2::SumSpace)
    if sort([S1.spaces...])==sort([S2.spaces...])
        # swaps sumspace order
        ConversionWrapper(SpaceOperator(
        PermutationOperator(promote_type(eltype(domain(S1)),eltype(domain(S2))),S1.spaces,S2.spaces),
                      S1,S2))
    elseif all(map(hasconversion,S1.spaces,S2.spaces))
        # we can blocmk convert
        ConversionWrapper(sumblkdiagm([map(Conversion,S1.spaces,S2.spaces)...]))
    elseif map(canonicalspace,S1.spaces)==map(canonicalspace,S2.spaces)
        error("Not implemented")
    elseif sort([map(canonicalspace,S1.spaces)...])==sort([map(canonicalspace,S2.spaces)...])
        # we can block convert after permuting
        P=PermutationOperator(promote_type(eltype(domain(S1)),eltype(domain(S2))),
                              map(canonicalspace,S1.spaces),
                              map(canonicalspace,S2.spaces))
        ds2=SumSpace(S1.spaces[P.perm])
        ConversionWrapper(TimesOperator(Conversion(ds2,S2),SpaceOperator(P,S1,ds2)))
    elseif all(map(hasconversion,sort([map(canonicalspace,S1.spaces)...]),sort([map(canonicalspace,S2.spaces)...])))
        #TODO: general case
        @assert length(S1.spaces)==2
        ds2=SumSpace(S1.spaces[[2,1]])
        TimesOperator(Conversion(ds2,S2),Conversion(S1,ds2))
    else
        # we don't know how to convert so go to default
        defaultconversion(S1,S2)
    end
end




for TYP in (:SumSpace,:PiecewiseSpace),  (OPrule,OP) in ((:conversion_rule,:conversion_type),(:maxspace_rule,:maxspace))
    @eval begin
        function $OPrule(S1::$TYP,S2::$TYP)
            cs1,cs2=map(canonicalspace,S1.spaces),map(canonicalspace,S2.spaces)
            if length(S1.spaces)!=length(S2.spaces)
                NoSpace()
            elseif canonicalspace(S1)==canonicalspace(S2)  # this sorts S1 and S2
                S1 ≤ S2?S1:S2  # choose smallest space by sorting
            elseif cs1==cs2
                # we can just map down
                # $TYP(map($OP,S1.spaces,S2.spaces))
                # this is commented out due to Issue #13261
                newspaces=[$OP(S1.spaces[k],S2.spaces[k]) for k=1:length(S1.spaces)]
                if any(b->b==NoSpace(),newspaces)
                    NoSpace()
                else
                    $TYP(newspaces)
                end
            elseif sort([cs1...])== sort([cs2...])
                # sort S1
                p=perm(cs1,cs2)
                $OP($TYP(S1.spaces[p]),S2)
            elseif length(S1.spaces)==length(S2.spaces)==2  &&
                    $OP(S1.spaces[1],S2.spaces[2])!=NoSpace() &&
                    $OP(S1.spaces[2],S2.spaces[1])!=NoSpace()
                #TODO: general length
                $TYP($OP(S1.spaces[1],S2.spaces[2]),
                     $OP(S1.spaces[2],S2.spaces[1]))
            else
                NoSpace()
            end
        end
    end
end



## Derivative

#TODO: do in @calculus_operator?
Derivative(S::SumSpace,k)=DerivativeWrapper(sumblkdiagm([Derivative(S.spaces[1],k),Derivative(S.spaces[2],k)]),k)
Integral(S::SumSpace,k)=IntegralWrapper(sumblkdiagm([Integral(S.spaces[1],k),Integral(S.spaces[2],k)]),k)




## Multiplcation for Array*Vector

function Multiplication{S,T,DD}(f::Fun{MatrixSpace{S,T,DD,1}},sp::VectorSpace)
    @assert size(space(f),2)==length(sp)
    m=mat(f)
    MultiplicationWrapper(f,interlace(BandedOperator{promote_type(eltype(f),eltype(sp))}[Multiplication(m[k,j],sp.space) for k=1:size(m,1),j=1:size(m,2)]))
end



## Multiply pieces

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
function addentries!{S<:SumSpace,SS<:SumSpace}(M::Multiplication{S,SS},A,k)
    a,b=vec(M.f)
    sp=domainspace(M)
    addentries!(Multiplication(a,sp)+Multiplication(b,sp),A,k)
end




function bandinds{S,SS<:SumSpace}(M::Multiplication{S,SS})
    a,b=vec(domainspace(M))
    Ma=Multiplication(M.f,a)
    Mb=Multiplication(M.f,b)

    bandinds(DiagonalInterlaceOperator([Ma,Mb]))
end
function rangespace{S,SS<:SumSpace}(M::Multiplication{S,SS})
    a,b=vec(domainspace(M))
    Ma=Multiplication(M.f,a)
    Mb=Multiplication(M.f,b)

    rangespace(Ma)⊕rangespace(Mb)
end
function addentries!{S,SS<:SumSpace}(M::Multiplication{S,SS},A,k)
    a,b=vec(domainspace(M))
    Ma=Multiplication(M.f,a)
    Mb=Multiplication(M.f,b)

    addentries!(DiagonalInterlaceOperator([Ma,Mb]),A,k)
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

Base.getindex{PWS1<:PiecewiseSpace,PWS2<:PiecewiseSpace}(d::TensorSpace{@compat(Tuple{PWS1,PWS2})},i::Integer,j::Integer)=d[1][i]⊗d[2][j]
Base.getindex{PWS1<:PiecewiseSpace,PWS2<:PiecewiseSpace}(d::TensorSpace{@compat(Tuple{PWS1,PWS2})},i::Range,j::Range)=PiecewiseSpace(d[1][i])⊗PiecewiseSpace(d[2][j])

## ProductFun

##  Piecewise

function pieces{PS<:PiecewiseSpace}(U::ProductFun{PS})
    ps=space(U,1)
    sp2=space(U,2)
    m=length(ps)
    C=coefficients(U)
    [ProductFun(C[k:m:end,:],ps[k],sp2) for k=1:m]
end
