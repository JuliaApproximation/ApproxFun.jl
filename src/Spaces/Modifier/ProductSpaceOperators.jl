export sumblkdiagm



immutable DiagonalPiecewiseOperator{T<:Number,B<:Operator} <: AbstractDiagonalInterlaceOperator{T,B}
    ops::Vector{B}
end

DiagonalPiecewiseOperator{B<:Operator}(v::Vector{B})=DiagonalPiecewiseOperator{mapreduce(eltype,promote_type,v),B}(v)
DiagonalPiecewiseOperator(v::Vector{Any})=DiagonalPiecewiseOperator(Operator{mapreduce(eltype,promote_type,v)}[v...;])


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
    ConversionWrapper(DiagonalPiecewiseOperator([Conversion(f[k],g[k]) for k=1:length(f)]))
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
    dimensions::(Int...)
end

DiagonalArrayOperator{T}(op::BandedOperator{T},dms::(Int...))=DiagonalArrayOperator{typeof(op),T}(op,dms)
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

ToeplitzOperator{S,T,V}(G::Fun{ArraySpace{S,2,T,1},V})=interlace(map(ToeplitzOperator,mat(G)))
ToeplitzOperator{S,T,V}(G::Fun{ArraySpace{S,2,T,1},V})=interlace(map(ToeplitzOperator,mat(G)))

## Sum Space


immutable SumInterlaceOperator{T<:Number,B<:Operator} <: AbstractDiagonalInterlaceOperator{T,B}
    ops::Vector{B}
end

SumInterlaceOperator{B<:Operator}(v::Vector{B})=SumInterlaceOperator{mapreduce(eltype,promote_type,v),B}(v)
SumInterlaceOperator(v::Vector{Any})=SumInterlaceOperator(Operator{mapreduce(eltype,promote_type,v)}[v...])

Base.convert{BT<:Operator}(::Type{BT},SI::SumInterlaceOperator)=SumInterlaceOperator(BandedOperator{eltype(BT)}[SI.ops...])

for op in (:domainspace,:rangespace)
    @eval $op(S::SumInterlaceOperator)=SumSpace($op(S.ops[1]),$op(S.ops[2]))
end

sumblkdiagm{B<:Operator}(v::Vector{B})=SumInterlaceOperator(v)



## Conversion
# swaps sumspace order

immutable BiSwapOperator <: BandedOperator{Float64} end
bandinds(::BiSwapOperator)=-1,1

function addentries!(::BiSwapOperator,A,kr::Range)
    for k=kr
        if isodd(k)
            A[k,k+1] += 1
        else
            A[k,k-1] += 1
        end
    end

    A
end



function Conversion{A,B,T}(S1::SumSpace{A,B,T},S2::SumSpace{B,A,T})
    @assert S1.spaces[1]==S2.spaces[2] && S1.spaces[2]==S2.spaces[1]
    ConversionWrapper(SpaceOperator(BiSwapOperator(),S1,S2))
end


function conversion_type{A,B,T}(S1::SumSpace{A,B,T},S2::SumSpace{B,A,T})
    if S1.spaces[1]==S2.spaces[2] && S1.spaces[2]==S2.spaces[1]
        S1 #Arbitraty
    else
        NoSpace()
    end
end


## Derivative

#TODO: do in @calculus_operator?
Derivative(S::SumSpace,k::Integer)=DerivativeWrapper(sumblkdiagm([Derivative(S.spaces[1],k),Derivative(S.spaces[2],k)]),k)
Integral(S::SumSpace,k::Integer)=IntegralWrapper(sumblkdiagm([Integral(S.spaces[1],k),Integral(S.spaces[2],k)]),k)




## Multiplcation for Array*Vector

function Multiplication{S,T,Q}(f::Fun{ArraySpace{S,2,T,1}},sp::ArraySpace{Q,1})
    @assert size(space(f),2)==length(sp)
    m=mat(f)
    MultiplicationWrapper(f,interlace(BandedOperator{promote_type(eltype(f),eltype(sp))}[Multiplication(m[k,j],sp.space) for k=1:size(m,1),j=1:size(m,2)]))
end

