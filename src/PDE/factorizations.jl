#############
#  PDEOperatorSchur
# represent a splitting rank 2 operator that's been discretized in the Y direction using
# the Schur decomposition
#
#  S is an AbstractOperatorSchur that encodes the schur decomposition
#############


abstract AbstractPDEOperatorSchur

immutable PDEOperatorSchur{OS<:AbstractOperatorSchur,LT<:Number,MT<:Number,ST<:Number,FT} <: AbstractPDEOperatorSchur
    Bx::FT
    Lx::Operator{LT}
    Mx::Operator{MT}
    S::OS

    indsBx::Vector{Int}
    indsBy::Vector{Int}

    Rdiags::Vector{SavedBandedOperator{ST}}
end

#TODO: Functional type?
Base.eltype{OS,LT,MT,ST,FT}(PO::PDEOperatorSchur{OS,LT,MT,ST,FT})=promote_type(eltype(PO.S),LT,MT,ST)


function PDEOperatorSchur{LT<:Number,MT<:Number,BT<:Number,ST<:Number}(Bx,Lx::Operator{LT},Mx::Operator{MT},S::AbstractOperatorSchur{BT,ST},indsBx,indsBy)
    Xops=promotespaces([Lx,Mx])
    Lx=SavedBandedOperator(Xops[1]);Mx=SavedBandedOperator(Xops[2])

    ny=size(S,1)
    nbcs=numbcs(S)
    RT=promote_type(LT,MT,BT,ST)
    Rdiags=Array(SavedBandedOperator{RT},ny)
    resizedata!(Lx,ny);resizedata!(Mx,ny)


    for k=1:ny-nbcs
        ##TODO: Do block case
        Rdiags[k]=SavedBandedOperator(PlusOperator(Operator{RT}[getdiagonal(S,k,1)*Lx,getdiagonal(S,k,2)*Mx]))
        resizedata!(Rdiags[k],ny)
    end

    if isempty(Bx)
        funcs=nothing
    elseif length(Bx)==1
        funcs=cache(Bx[1])
        resizedata!(funcs,:,ny)
    else
        funcs=cache(InterlaceOperator(Bx))
        resizedata!(funcs,:,ny)
    end


    PDEOperatorSchur(funcs,Lx,Mx,S,indsBx,indsBy,Rdiags)
end




PDEOperatorSchur(Bx::Vector,Lx::Operator,Mx::Operator,S::AbstractOperatorSchur) =
    PDEOperatorSchur(Bx,Lx,Mx,S,[1:length(Bx)],length(Bx)+[1:numbcs(S)])
PDEOperatorSchur(Bx::Vector,Lx::Operator,Mx::UniformScaling,S::AbstractOperatorSchur) =
    PDEOperatorSchur(Bx,Lx,ConstantOperator(Mx.λ),S)
PDEOperatorSchur(Bx::Vector,Lx::UniformScaling,Mx::Operator,S::AbstractOperatorSchur) =
    PDEOperatorSchur(Bx,ConstantOperator(Lx.λ),Mx,S)


function PDEOperatorSchur(Bx,By,A::Operator,ny::Integer,indsBx,indsBy)
    @assert iskronsumop(A)
    ops=sumops(A)
    @assert all(isproductop,ops)
    @assert length(ops)==2

    PDEOperatorSchur(Bx,dekron(ops[1],1),dekron(ops[2],1),schurfact(By,[dekron(ops[1],2),dekron(ops[2],2)],ny),indsBx,indsBy)
end




bcinds(S::PDEOperatorSchur,k)=k==1?S.indsBx:S.indsBy
numbcs(S::AbstractPDEOperatorSchur,k)=length(bcinds(S,k))
Base.length(S::PDEOperatorSchur)=size(S.S,1)


function PDEOperatorSchur(A::Vector,ny::Integer)
    indsBx,Bx=findfunctionals(A,1)
    indsBy,By=findfunctionals(A,2)

    LL=A[end]


    PDEOperatorSchur(Bx,By,LL,ny,indsBx,indsBy)
end



for op in (:domainspace,:rangespace)
    @eval begin
        $op(P::PDEOperatorSchur,k::Integer)=k==1?$op(P.Lx):$op(P.S)
        $op(L::PDEOperatorSchur)=$op(L,1)⊗$op(L,2)
    end
end


domain(P::PDEOperatorSchur,k::Integer)=k==1?domain(P.Lx):domain(P.S)
domain(P::PDEOperatorSchur)=domain(domainspace(P))




## Constructuor


Base.schurfact{LT<:Number,MT<:Number,BT<:Number,ST<:Number}(Bx,Lx::Operator{LT},Mx::Operator{MT},S::AbstractOperatorSchur{BT,ST},indsBx,indsBy) =
    PDEOperatorSchur(Bx,Lx,Mx,S,indsBx,indsBy)

function Base.schurfact(Bx,By,A::Operator,ny::Integer,indsBx,indsBy)
    @assert iskronsumop(A)
    opsx,opsy=simplifydekron(A)
    @assert length(opsx)==length(opsy)==2

    schurfact(Bx,
              opsx[1],opsx[2],
              schurfact(By,opsy,ny),
              indsBx,indsBy)
end



function Base.schurfact{BT<:Operator}(A::Vector{BT},S::ProductDomain,ny::Integer)
    indsBx,Bx=findfunctionals(A,1)
    indsBy,By=findfunctionals(A,2)

    LL=A[end]

    schurfact(Bx,By,LL,ny,indsBx,indsBy)
end



Base.schurfact{BT<:Operator}(A::Vector{BT},n::Integer)=schurfact(A,domain(A[end]),n)
Base.schurfact(A::Operator,n::Integer)=schurfact([A],n)





## discretize


function discretize{OT<:Operator}(A::Vector{OT},S...)
    if iskronsumop(A[end])&&length(simplifydekron(A[end])[1])==2
        schurfact(A,S...)
    elseif isdiagop(A[end],2)
        schurfact(A,S...)
    else
        kronfact(A,S...)
    end
end
discretize(A::Operator,n::Integer)=discretize([A],n)
