export sumblkdiagm



immutable DiagonalPiecewiseOperator{T<:Number,B<:Operator} <: AbstractDiagonalInterlaceOperator{T,B}
    ops::Vector{B}
end

DiagonalPiecewiseOperator{B<:Operator}(v::Vector{B})=DiagonalPiecewiseOperator{mapreduce(eltype,promote_type,v),B}(v)
DiagonalPiecewiseOperator(v::Vector{Any})=DiagonalPiecewiseOperator(Operator{mapreduce(eltype,promote_type,v)}[v...])


for op in (:domainspace,:rangespace)
    @eval $op(D::DiagonalInterlaceOperator)=depiece(map($op,D.ops))
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

function addentries!(::BiSwapOperator,A::ShiftArray,kr::Range)
    for k=kr
        if isodd(k)
            A[k,1] += 1
        else
            A[k,-1] += 1
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


