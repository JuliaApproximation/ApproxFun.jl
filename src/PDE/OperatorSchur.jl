


toarray{T<:Operator}(B::Array{T},n)=Float64[    B[k][j] for  k=1:length(B),j=1:n];
function toarray{T<:Operator}(A::Vector{T},n::Integer,m::Integer)
    ret = zeros(n,m)

    nbc = isafunctional(A[end])?length(A):length(A)-1
    for k=1:nbc
        ret[k,:]=A[k][1:m]
    end

    if nbc < length(A)
        ret[nbc+1:end,:]=A[end][1:n-nbc,1:m]
    end

    ret
end

function pdetoarray(Byin,Lin,Min,ny::Integer)
    Yop=promotespaces([Lin,Min])


    By=toarray(Byin,ny)
    nbcy=length(Byin)

    Ly=Yop[1][1:ny-nbcy,1:ny]|>full
    My=Yop[2][1:ny-nbcy,1:ny]|>full

    By,Ly,My
end

function pdetoarray(Byin,L::Vector,ny::Integer)
    Yop=promotespaces(L)
    By=toarray(Byin,ny)
    nbcy=length(Byin)
    By,[Yk[1:ny-nbcy,1:ny] for Yk in Yop]
end


function cont_reduce_dofs( R,A::Array )
    if length(R) > 0
        for k = 1:size(R,1)
            A = A - A[:,k]*R[k:k,:]
        end
    end

    A
end


abstract AbstractOperatorSchur{BT,MT}

domainspace(S::AbstractOperatorSchur)=S.domainspace
rangespace(S::AbstractOperatorSchur)=S.rangespace
domain(S::AbstractOperatorSchur)=domain(domainspace(S))


type DiagonalOperatorSchur{MT<:Number} <:AbstractOperatorSchur{MT,MT}
    ops::Vector{Vector{MT}}


    domainspace::Space
    rangespace::Space
end

Base.eltype{MT}(::DiagonalOperatorSchur{MT})=MT

function DiagonalOperatorSchur{T1<:Number,T2<:Number}(R::Vector{T1},T::Vector{T2},d,r)
    PT=promote_type(T1,T2)
    A=Array(Vector{PT},2)
    A[1]=convert(Vector{PT},R)
    A[2]=convert(Vector{PT},T)
    DiagonalOperatorSchur(A,d,r)
end

function DiagonalOperatorSchur(L::Operator,M::Operator,n::Integer)
    Yop=promotespaces([L,M])
    DiagonalOperatorSchur(diag(Yop[1][1:n,1:n]|>full),diag(Yop[2][1:n,1:n]|>full),domainspace(Yop[1]),rangespace(Yop[2]))
end


function DiagonalOperatorSchur{O<:Operator}(L::Vector{O},n::Integer)
    Yop=promotespaces(L)

    ##TODO: type

    ops=Array(Vector{mapreduce(eltype,promote_type,L)},length(Yop))
    for k=1:length(L)
        ops[k]=diag(Yop[k][1:n,1:n])
    end

    DiagonalOperatorSchur(ops,domainspace(Yop[1]),rangespace(Yop[2]))
end

Base.size(S::DiagonalOperatorSchur,k)=length(S.ops[1])
Base.size(S::DiagonalOperatorSchur)=size(S,1),size(S,2)

getdiagonal(S::DiagonalOperatorSchur,k,j)=S.ops[j][k]

numbcs(::DiagonalOperatorSchur)=0


type OperatorSchur{BT<:Number,MT<:Number} <:AbstractOperatorSchur{BT,MT}
    bcP::Matrix{BT}  # permute columns of bcs
    bcQ::Matrix{BT}  # bc normalizer

    bcs::Matrix{BT} # normalized bcs

    # C == QRZ',  D == QTZ'
    # where C/D have degrees of freedom removed
    Q::Matrix{MT}
    Z::Matrix{MT}

    R::Matrix{MT}
    T::Matrix{MT}

    # L[:,1:k] and M[:,1:k]  so we know how the columns are killed
    Lcols::Matrix{MT}
    Mcols::Matrix{MT}

    domainspace::Space
    rangespace::Space
end

#make sure cols are same type as ops
OperatorSchur{M<:Number}(bcP,bcQ,bcs,Q::Array{M,2},Z::Array{M,2},R::Array{M,2},T::Array{M,2},Lcols::Array,Mcols::Array,ds,rs)=
    OperatorSchur(bcP,bcQ,bcs,Q,Z,R,T,convert(Array{M},Lcols),convert(Array{M},Mcols),ds,rs)

Base.size(S::OperatorSchur)=size(S.bcP)
Base.size(S::OperatorSchur,k)=size(S.bcP,k)

Base.eltype{BT,MT}(::OperatorSchur{BT,MT})=promote_type(BT,MT)

numbcs(S::OperatorSchur)=size(S.bcQ,1)

getdiagonal(S::OperatorSchur,k,j)=j==1?S.R[k,k]:S.T[k,k]



OperatorSchur{FT<:Operator}(B::Vector{FT},L::UniformScaling,M::Operator,n::Integer) =
    OperatorSchur(B,ConstantOperator(L),M,n)
OperatorSchur{FT<:Operator}(B::Vector{FT},L::Operator,M::UniformScaling,n::Integer) =
    OperatorSchur(B,L,ConstantOperator(M),n)




OperatorSchur(B,L::SparseMatrixCSC,M::SparseMatrixCSC,ds,rs)=OperatorSchur(B,full(L),full(M),ds,rs)
function OperatorSchur(B::Array,L::Array,M::Array,ds,rs)
    B,Q,L,M,P=regularize_bcs(B,L,M)

    K = size(B,1)

    Lcols=L[:,1:K];Mcols=M[:,1:K]

    L=cont_reduce_dofs(B,L)
    M=cont_reduce_dofs(B,M)

    C=L[:,K+1:end]
    D=M[:,K+1:end]
    CD=schurfact(C,D)
    Q2=CD[:left];Z2=CD[:right]
    R=CD[:S]; T=CD[:T]

    OperatorSchur(P,Q,B,Q2,Z2,R,T, Lcols,Mcols,ds,rs)
end


function OperatorSchur{FT<:Operator}(B::Vector{FT},L::Operator,M::Operator,n::Integer)
    L,M=promotespaces([L,M])
    OperatorSchur(pdetoarray(B,L,M,n)...,domainspace(L),rangespace(L))
end


#####
# StrideOperatorSchur
####

immutable StrideOperatorSchur{MT<:Number} <:AbstractOperatorSchur{Float64,MT}
    odd::OperatorSchur{Float64,MT}
    even::OperatorSchur{Float64,MT}
end


function StrideOperatorSchur(L,M,n)
    L,M=promotespaces([L,M])
    L1=DestrideOperator(L,-1,-1,2,2);M1=DestrideOperator(M,-1,-1,2,2)
    L2=DestrideOperator(L,0,0,2,2);M2=DestrideOperator(M,0,0,2,2)

    B=FillFunctional(2.0)
    O1=OperatorSchur([B],L1,M1,div(n,2))
    O2=OperatorSchur([B],L2,M2,div(n,2))

    StrideOperatorSchur(O1,O2)
end

for OP in (:domainspace,:rangespace)
    @eval $OP(S::StrideOperatorSchur)=$OP(S.odd)
end

## Decide which data structure


function Base.schurfact{FT<:Operator,O<:Operator}(B::Vector{FT},A::Vector{O},n::Integer)
    if isempty(B) && all(LM->bandinds(LM)==(0,0),A)
        DiagonalOperatorSchur(A,n)
    elseif length(A)==2
        L,M=promotespaces([A[1],A[2]])
#         if length(B)==2 &&
#                      gcd(stride(L),stride(M))==2 &&
#                      isa(B[1],ConcreteEvaluation{Ultraspherical{0},Bool,Int,Float64}) &&
#                      isa(B[2],ConcreteEvaluation{Ultraspherical{0},Bool,Int,Float64}) &&
#                      !B[1].x && B[2].x
#              StrideOperatorSchur(L,M,n)
#         else
            OperatorSchur(B,L,M,n)
#        end
    else
        error("Schur factorization unknown for more than 2 non-diagonal operators.")
    end
end
