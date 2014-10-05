
export lap,grad,timedirichlet, ⊗


type PDEOperator
    # 2x2 array 
    # [Lx Ly;
    # [Mx My] 
    # that represents
    # Lx⊗Ly + Mx⊗My
 
    ops::Array{Operator,2}
end



PDEOperator(A,B)=PDEOperator([A B])


domainspace(L::PDEOperator,j::Integer)=findmindomainspace(L.ops[:,j])
rangespace(L::PDEOperator,j::Integer)=findmaxrangespace(L.ops[:,j])

for op in (:domainspace,:rangespace)
    @eval begin
        $op(L::PDEOperator)=$op(L,1)⊗$op(L,2)
    end
end


function domain(LL::PDEOperator,j::Integer)
    for k=1:size(LL.ops,1)
        dx=domain(LL.ops[k,j]) 
        if dx != AnyDomain()
            return dx
        end
    end
    return AnyDomain()
end

function domain(LL::PDEOperator)
    @assert size(LL.ops,2)==2
    domain(LL,1)⊗domain(LL,2)
end


Base.transpose(LL::PDEOperator)=PDEOperator(LL.ops[:,end:-1:1])


## Space promotion


function promotedomainspace(P::PDEOperator,S::FunctionSpace,col::Integer)
    if col==1
        PDEOperator([promotedomainspace(P.ops[1,1],S) P.ops[1,2];
                     promotedomainspace(P.ops[2,1],S) P.ops[2,2]])    
    else
        @assert col==2
        PDEOperator([P.ops[1,1] promotedomainspace(P.ops[1,2],S);
                     P.ops[2,1] promotedomainspace(P.ops[2,2],S)])        
    end
end
promotedomainspace(P::PDEOperator,S::TensorSpace)=promotedomainspace(promotedomainspace(P,S[1],1),S[2],2)



## Algebra

⊗(A::Operator,B::Operator)=PDEOperator(A,B)
⊗(A::Operator,B::UniformScaling)=A⊗ConstantOperator(1.0B.λ)
⊗(A::UniformScaling,B::Operator)=ConstantOperator(1.0A.λ)⊗B

⊗(A::Fun,B::Fun)=Multiplication(A)⊗Multiplication(B)
⊗(A,B::Fun)=A⊗Multiplication(B)
⊗(A::Fun,B)=Multiplication(A)⊗B

⊗{T<:Operator}(A::Vector{T},B::Operator)=PDEOperator[PDEOperator(Ai,B) for Ai in A]
⊗{T<:Operator}(A::Operator,B::Vector{T})=PDEOperator[PDEOperator(A,Bi) for Bi in B]
⊗{T<:Operator}(A::Vector{T},B::UniformScaling)=A⊗ConstantOperator(1.0B.λ)
⊗{T<:Operator}(A::UniformScaling,B::Vector{T})=ConstantOperator(1.0A.λ)⊗B



function +(A::PDEOperator,B::PDEOperator)
    ret=copy(A.ops)
    for k=1:size(B.ops,1)
        ##TODO: might not be ordered
        if ret[k,1]==B.ops[k,1]
            ret[k,2]+=B.ops[k,2]
        elseif ret[k,2]==B.ops[k,2]
            ret[k,1]+=B.ops[k,1]            
        else
            ret=[ret;B.ops[k,:]]
        end
    end
    PDEOperator(ret)
end

function lap(d::ProductDomain)
    @assert length(d.domains)==2
    Dx=Derivative(d.domains[1])
    Dy=Derivative(d.domains[2])    
    Dx^2⊗I+I⊗Dy^2
end


function -(A::PDEOperator)
    ops=copy(A.ops)
    for k=1:size(ops,1)
        ops[k,1]=-ops[k,1]
    end
    PDEOperator(ops)
end

+(A::UniformScaling,B::PDEOperator)=B+ConstantOperator(1.0A.λ)⊗ConstantOperator(1.0)
+(B::PDEOperator,A::UniformScaling)=B+ConstantOperator(1.0A.λ)⊗ConstantOperator(1.0)
-(A::UniformScaling,B::PDEOperator)=-B+ConstantOperator(1.0A.λ)⊗ConstantOperator(1.0)
-(B::PDEOperator,A::UniformScaling)=B+ConstantOperator(-1.0A.λ)⊗ConstantOperator(1.0)

-(A::PDEOperator,B::PDEOperator)=A+(-B)

function *(c::Number,A::PDEOperator)
    ops=copy(A.ops)
    for k=1:size(ops,1)
        ops[k,1]=c*ops[k,1]
    end
    PDEOperator(ops)
end
*(A::PDEOperator,c::Number)=c*A
function *(A::PDEOperator,B::PDEOperator)
    # TODO: higher rank operators
    if A==B
        @assert size(A.ops,1)==size(B.ops,1)==2
        @assert size(A.ops,2)==size(B.ops,2)==2    
        PDEOperator([A.ops[1,1]^2 A.ops[1,2]^2;
                     2A.ops[1,1]*A.ops[2,1] A.ops[1,2]*A.ops[2,2];
                     A.ops[2,1]^2 A.ops[2,2]^2
                    ])    
    else
        @assert size(A.ops,1)==size(B.ops,1)==1
        @assert size(A.ops,2)==size(B.ops,2)==2
        PDEOperator([A.ops[1,1]*B.ops[1,1] A.ops[1,2]*B.ops[1,2]])
    end
end



##TODO how to determine whether x or y?
function *(a::Fun,A::PDEOperator)
    ops = copy(A.ops)
    for k=1:size(ops,1)
        ops[k,1]=a*ops[k,1]
    end
    PDEOperator(ops)
end

Base.diff(d::ProductDomain,k)=k==1?Base.diff(d.domains[1])⊗I:I⊗Base.diff(d.domains[2])  
grad(d::ProductDomain)=[Base.diff(d,k) for k=1:length(d.domains)]




function dirichlet(d::ProductDomain)
    @assert length(d.domains)==2
    Bx=dirichlet(d.domains[1])
    By=dirichlet(d.domains[2])
    [Bx⊗I,I⊗By]
end

function neumann(d::ProductDomain)
    @assert length(d.domains)==2
    Bx=neumann(d.domains[1])
    By=neumann(d.domains[2])
    [Bx⊗I,I⊗By]
end

function timedirichlet(d::ProductDomain)
    @assert length(d.domains)==2
    Bx=dirichlet(d.domains[1])
    Bt=dirichlet(d.domains[2])[1]
    [I⊗Bt,Bx⊗I]
end


function *{S,T}(L::PDEOperator,f::LowRankFun{S,T})
    @assert size(L.ops,2)==2
    @assert size(L.ops,1)==2    
    n=length(f.A)
    A=Array(Fun,2n)
    B=Array(Fun,2n)
    
    for k=1:n
        A[k]=L.ops[1,1]*f.A[k]
        B[k]=L.ops[1,2]*f.B[k]        
        A[k+n]=L.ops[2,1]*f.A[k]
        B[k+n]=L.ops[2,2]*f.B[k]        
    end
    
    LowRankFun(A,B)
end



## Schur PDEOperator
# represent an operator that's been discretized in the Y direction

abstract AbstractPDEOperatorSchur

type PDEOperatorSchur{OS<:AbstractOperatorSchur,LT<:Number,MT<:Number,ST<:Number,FT<:Functional} <: AbstractPDEOperatorSchur
    Bx::Vector{FT}
    Lx::Operator{LT}
    Mx::Operator{MT}
    S::OS
    
    indsBx::Vector{Int}
    indsBy::Vector{Int}
    
    Rdiags::Vector{SavedBandedOperator{ST}}
end


function PDEOperatorSchur{LT<:Number,MT<:Number,BT<:Number,ST<:Number}(Bx,Lx::Operator{LT},Mx::Operator{MT},S::AbstractOperatorSchur{BT,ST},indsBx,indsBy)
    ny=size(S,1)
    nbcs=numbcs(S)
    Rdiags=Array(SavedBandedOperator{promote_type(LT,MT,BT,ST)},ny)
    Xops=promotespaces([Lx,Mx])
    Lx=SavedBandedOperator(Xops[1]);Mx=SavedBandedOperator(Xops[2])
    
    resizedata!(Lx,ny);resizedata!(Mx,ny)
    
    
    for k=1:ny-nbcs
        ##TODO: Do block case
        Rdiags[k]=SavedBandedOperator(getdiagonal(S,k,1)*Lx + getdiagonal(S,k,2)*Mx)
        resizedata!(Rdiags[k],ny)
    end

    
    PDEOperatorSchur(Bx,Lx,Mx,S,indsBx,indsBy,Rdiags)
end


PDEOperatorSchur(Bx::Vector,Lx::Operator,Mx::Operator,S::AbstractOperatorSchur)=PDEOperatorSchur(Bx,Lx,Mx,S,[1:length(Bx)],length(Bx)+[1:numbcs(S)])
PDEOperatorSchur(Bx::Vector,Lx::Operator,Mx::UniformScaling,S::AbstractOperatorSchur)=PDEOperatorSchur(Bx,Lx,ConstantOperator(Mx.λ),S)
PDEOperatorSchur(Bx::Vector,Lx::UniformScaling,Mx::Operator,S::AbstractOperatorSchur)=PDEOperatorSchur(Bx,ConstantOperator(Lx.λ),Mx,S)


function PDEOperatorSchur(Bx,By,A::PDEOperator,ny::Integer,indsBx,indsBy)
   @assert size(A.ops)==(2,2)
   PDEOperatorSchur(Bx,A.ops[1,1],A.ops[2,1],OperatorSchur(By,A.ops[1,2],A.ops[2,2],ny),indsBx,indsBy)
end

isxfunctional(B::PDEOperator)=size(B.ops,1)==1&&size(B.ops,2)==2&&typeof(B.ops[1,1])<:Functional
isyfunctional(B::PDEOperator)=size(B.ops,1)==1&&size(B.ops,2)==2&&typeof(B.ops[1,2])<:Functional
ispdeop(B::PDEOperator)=!isxfunctional(B)&&!isyfunctional(B)


bcinds(S::PDEOperatorSchur,k)=k==1?S.indsBx:S.indsBy
numbcs(S::AbstractPDEOperatorSchur,k)=length(bcinds(S,k))
Base.length(S::PDEOperatorSchur)=size(S.S,1)



function PDEOperatorSchur{T<:PDEOperator}(A::Vector{T},ny::Integer)
    indsBx=find(isxfunctional,A)
    Bx=Functional[(@assert Ai.ops[1,2]==ConstantOperator{Float64}(1.0); Ai.ops[1,1]) for Ai in A[indsBx]]
    indsBy=find(isyfunctional,A)
    By=Functional[(@assert Ai.ops[1,1]==ConstantOperator{Float64}(1.0); Ai.ops[1,2]) for Ai in A[indsBy]]
    inds=find(ispdeop,A)
    @assert length(inds)==1&&inds[1]==length(A)
    
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





## Product

# Represents an operator on e.g. a Disk
type PDEProductOperatorSchur{ST<:Number,FT<:Functional} <: AbstractPDEOperatorSchur
    Bx::Vector{FT}
    Rdiags::Vector{SavedBandedOperator{ST}}
    domainspace::AbstractProductSpace
end

Base.length(S::PDEProductOperatorSchur)=length(S.Rdiags)

function PDEProductOperatorSchur{T<:PDEOperator}(A::Vector{T},sp::AbstractProductSpace,nt::Integer)
    indsBx=find(isxfunctional,A)
    Bx=Functional[(@assert Ai.ops[1,2]==ConstantOperator{Float64}(1.0); Ai.ops[1,1]) for Ai in A[indsBx]]
    
    @assert length(Bx)==1  #TODO: this is too restrictve
    
    indsBy=find(isyfunctional,A)
    @assert length(indsBy)==0
    inds=find(ispdeop,A)
    @assert length(inds)==1&&inds[1]==length(A)
    
    L=A[end]
    
    
    
    L=promotedomainspace(L,sp[2],2)
    S=OperatorSchur([],L.ops[1,2],L.ops[2,2],nt)
    @assert(isa(S,DiagonalOperatorSchur))
    Lx=L.ops[1,1];Mx=L.ops[2,1]
    #TODO: Space Type
    BxV=Array(SavedFunctional{Float64},nt)
    Rdiags=Array(SavedBandedOperator{Complex{Float64}},nt)    
    for k=1:nt
        csp=columnspace(sp,k)
        Rdiags[k]=SavedBandedOperator(promotedomainspace(getdiagonal(S,k,1)*Lx+getdiagonal(S,k,2)*Mx,csp))
        BxV[k]=SavedFunctional(promotedomainspace(Bx[1],csp))
        resizedata!(Rdiags[k],2nt+100)
        resizedata!(BxV[k],2nt+100)        
    end  
    PDEProductOperatorSchur(BxV,Rdiags,sp)
end


##TODO: Larger Bx
bcinds(S::PDEProductOperatorSchur,k)=k==1?[1]:[]
domainspace(S::PDEProductOperatorSchur)=S.domainspace
domainspace(S::PDEProductOperatorSchur,k)=S.domainspace[k]

# for op in (:domainspace,:rangespace)
#     @eval begin
#         $op(P::PDEProductOperatorSchur,k::Integer)=k==1?$op(P.Lx):$op(P.S)
#        $op(L::PDEProductOperatorSchur)=$op(L,1)⊗$op(L,2)
#     end
# end






type ProductRangeSpace <: BivariateFunctionSpace
    S::PDEProductOperatorSchur
end

rangespace(S::PDEProductOperatorSchur)=ProductRangeSpace(S)

function space(S::ProductRangeSpace,k)
    @assert k==2
    S.S.domainspace[2]
end 

function coefficients{S,V,SS,T}(f::ProductFun{S,V,SS,T},sp::ProductRangeSpace)
    @assert space(f,2)==space(sp,2)
    
    n=size(f,2)
    F=[coefficients(f.coefficients[k],rangespace(sp.S.Rdiags[k])) for k=1:n]
    m=mapreduce(length,max,F)
    ret=zeros(T,m,n)
    for k=1:n
        ret[1:length(F[k]),k]=F[k]
    end
    ret    
end



## Constructuor

Base.schurfact{T<:PDEOperator}(A::Vector{T},n::Integer)=PDEOperatorSchur(A,n)
Base.schurfact(A::PDEOperator,n::Integer)=schurfact([A],n)
Base.schurfact{T<:PDEOperator}(A::Vector{T},S::AbstractProductSpace,n::Integer)=PDEProductOperatorSchur(A,S,n)


