export sumblkdiagm



immutable DiagonalPiecewiseOperator{T<:Number,B<:Operator} <: AbstractDiagonalInterlaceOperator{T,B}
    ops::Vector{B}
end

DiagonalPiecewiseOperator{B<:Operator}(v::Vector{B})=DiagonalPiecewiseOperator{mapreduce(eltype,promote_type,v),B}(v)
DiagonalPiecewiseOperator(v::Vector{Any})=DiagonalPiecewiseOperator(Operator{mapreduce(eltype,promote_type,v)}[v...])


for op in (:domainspace,:rangespace)
    @eval $op(D::DiagonalPiecewiseOperator)=PiecewiseSpace(map($op,D.ops))
end




Derivative(sp::PiecewiseSpace)=DerivativeWrapper(DiagonalPiecewiseOperator(map(Derivative,sp.spaces)),1)
Derivative(sp::PiecewiseSpace,k::Integer)=DerivativeWrapper(DiagonalPiecewiseOperator(map(s->Derivative(s,k),sp.spaces)),k)


function Multiplication{PW<:PiecewiseSpace}(f::Fun{PW},sp::PiecewiseSpace)
    vf=vec(f)
    @assert length(vf)==length(sp)
    MultiplicationWrapper(f,DiagonalPiecewiseOperator([Multiplication(vf[k],sp[k]) for k=1:length(vf)]))
end

function Conversion(f::PiecewiseSpace,g::PiecewiseSpace)
    @assert length(f)==length(g)
    ConversionWrapper(DiagonalPiecewiseOperator(Operator[Conversion(f[k],g[k]) for k=1:length(f)]))
end

for op in (:dirichlet,:neumann,:ivp)
    @eval $op(d::PiecewiseSpace)=interlace($op(d.spaces))
    @eval $op(d::UnionDomain)=interlace($op(d.domains))
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
function Conversion(AS::ArraySpace,BS::ArraySpace)
    @assert size(AS)==size(BS)
    ConversionWrapper(DiagonalArrayOperator(Conversion(AS.space,BS.space),size(AS)))
end


## Sum Space


immutable SumInterlaceOperator{T<:Number,B<:Operator} <: AbstractDiagonalInterlaceOperator{T,B}
    ops::Vector{B}
end

SumInterlaceOperator{B<:Operator}(v::Vector{B})=SumInterlaceOperator{mapreduce(eltype,promote_type,v),B}(v)
SumInterlaceOperator(v::Vector{Any})=SumInterlaceOperator(Operator{mapreduce(eltype,promote_type,v)}[v...])

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



function Conversion{A,B,T,D}(S1::SumSpace{A,B,T,D},S2::SumSpace{B,A,T,D})
    @assert S1.spaces[1]==S2.spaces[2] && S1.spaces[2]==S2.spaces[1]
    ConversionWrapper(SpaceOperator(BiSwapOperator(),S1,S2))
end


function conversion_type{A,B,T,D}(S1::SumSpace{A,B,T,D},S2::SumSpace{B,A,T,D})
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

function Multiplication{S,T,D,Q}(f::Fun{ArraySpace{S,2,T,D}},sp::ArraySpace{Q,1})
    @assert size(space(f),2)==length(sp)
    m=mat(f)
    MultiplicationWrapper(f,interlace(Operator[Multiplication(m[k,j],sp.space) for k=1:size(m,1),j=1:size(m,2)]))
end

