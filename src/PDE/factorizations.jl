#############
#  PDEOperatorSchur
# represent a splitting rank 2 operator that's been discretized in the Y direction using
# the Schur decomposition
#
#  S is an AbstractOperatorSchur that encodes the schur decomposition
#############


abstract AbstractPDEOperatorSchur

immutable PDEOperatorSchur{OS<:AbstractOperatorSchur,LT<:Number,MT<:Number,ST<:Number,FT<:Operator} <: AbstractPDEOperatorSchur
    Bx::Vector{FT}
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
        funcs=Array(SavedFunctional{Float64},0)
    else
        funcs=Array(SavedFunctional{mapreduce(eltype,promote_type,Bx)},length(Bx))
        for k=1:length(Bx)
            funcs[k]=SavedFunctional(Bx[k])
            resizedata!(funcs[k],ny)
        end
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



#############
# PDEStrideOperatorSchur
#############


immutable PDEStrideOperatorSchur  <: AbstractPDEOperatorSchur
    odd
    even
end


function PDEStrideOperatorSchur(Bx,Lx::Operator,Mx::Operator,S::StrideOperatorSchur,indsBx,indsBy)
    @assert indsBy==[3,4]
    Podd=PDEOperatorSchur(Bx,Lx,Mx,S.odd,indsBx,[indsBy[1]])
    Peven=PDEOperatorSchur(Bx,Lx,Mx,S.even,indsBx,[indsBy[1]])

    PDEStrideOperatorSchur(Podd,Peven)
end

###############
# Product
# Represents an operator on e.g. a Disk
#############

immutable PDEProductOperatorSchur{ST<:Number,
                                  FT<:Operator,
                                  DS<:AbstractProductSpace,
                                  RS<:AbstractProductSpace,
                                  S<:Space,
                                  V<:Space} <: AbstractPDEOperatorSchur
    Bx::Vector{Vector{FT}}
    Rdiags::Vector{SavedBandedOperator{ST}}
    domainspace::DS
    rangespace::RS
    indsBx::Vector{Int}
end

PDEProductOperatorSchur{ST,FT,S,V}(Bx::Vector{Vector{FT}},
                                   Rdiags::Vector{SavedBandedOperator{ST}},
                                   ds::AbstractProductSpace{Tuple{S,V}},rs,indsBx)=PDEProductOperatorSchur{ST,FT,
                                                                                                    typeof(ds),typeof(rs),
                                                                                                    S,V}(Bx,Rdiags,ds,rs,indsBx)

Base.eltype{ST}(::PDEProductOperatorSchur{ST})=ST

Base.length(S::PDEProductOperatorSchur)=length(S.Rdiags)

domain(P::PDEProductOperatorSchur)=domain(domainspace(P))

function PDEProductOperatorSchur{OT<:Operator}(A::Vector{OT},sp::AbstractProductSpace,nt::Integer)
    Bx=A[1:end-1]
    L=A[end]

    @assert isdiagop(L,2)

    BxV=Array(Vector{SavedFunctional{Float64}},nt)
    Rdiags=Array(SavedBandedOperator{Complex{Float64}},nt)
    for k=1:nt
        op=diagop(L,k)
        csp=columnspace(sp,k)
        Rdiags[k]=SavedBandedOperator(op)
        #TODO: Type
        BxV[k]=SavedFunctional{Float64}[SavedFunctional(diagop(Bxx,k)) for Bxx in Bx]
        resizedata!(Rdiags[k],2nt+100)
        for Bxx in BxV[k]
            resizedata!(Bxx,2nt+100)
        end
    end
    PDEProductOperatorSchur(BxV,Rdiags,sp,rangespace(L),1:length(Bx))
end

PDEProductOperatorSchur{OT<:Operator}(A::Vector{OT},sp::BivariateDomain,nt::Integer)=PDEProductOperatorSchur(A,Space(sp),nt)


##TODO: Larger Bx
bcinds(S::PDEProductOperatorSchur,k)=k==1?S.indsBx:[]
domainspace(S::PDEProductOperatorSchur)=S.domainspace
domainspace(S::PDEProductOperatorSchur,k)=S.domainspace[k]
rangespace(S::PDEProductOperatorSchur)=S.rangespace


## Constructuor


Base.schurfact{LT<:Number,MT<:Number,BT<:Number,ST<:Number}(Bx,Lx::Operator{LT},Mx::Operator{MT},S::AbstractOperatorSchur{BT,ST},indsBx,indsBy)=PDEOperatorSchur(Bx,Lx,Mx,S,indsBx,indsBy)
Base.schurfact{LT<:Number,MT<:Number,ST<:Number}(Bx,Lx::Operator{LT},Mx::Operator{MT},S::StrideOperatorSchur{ST},indsBx,indsBy)=PDEStrideOperatorSchur(Bx,Lx,Mx,S,indsBx,indsBy)

function Base.schurfact{T}(Bx,By,A::Operator{BandedMatrix{T}},ny::Integer,indsBx,indsBy)
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





Base.schurfact{BT<:Operator}(A::Vector{BT},S::BivariateDomain,n::Integer)=PDEProductOperatorSchur(A,S,n)


Base.schurfact{BT<:Operator}(A::Vector{BT},n::Integer)=schurfact(A,domain(A[end]),n)
Base.schurfact{T}(A::BivariateOperator{T},n::Integer)=schurfact([A],n)

function *(A::PDEProductOperatorSchur,F::ProductFun)
    ret=copy(F.coefficients)
    for k=1:length(ret)
        ret[k]=A.Rdiags[k]*ret[k]
    end
    ProductFun(ret,space(F))
end



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
discretize{T}(A::BivariateOperator{T},n::Integer)=discretize([A],n)
