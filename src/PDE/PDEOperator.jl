
export lap,grad,timedirichlet


type PDEOperator
    ops::Array{Operator,2}
end



PDEOperator(A,B)=PDEOperator([A B])



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


⊗(A::Operator,B::Operator)=PDEOperator(A,B)
⊗(A::Operator,B::UniformScaling)=A⊗ConstantOperator(1.0B.λ)
⊗(A::UniformScaling,B::Operator)=ConstantOperator(1.0A.λ)⊗B

⊗(A::Fun,B::Fun)=MultiplicationOperator(A)⊗MultiplicationOperator(B)
⊗(A,B::Fun)=A⊗MultiplicationOperator(B)
⊗(A::Fun,B)=MultiplicationOperator(A)⊗B

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

function lap(d::TensorDomain)
    @assert length(d.domains)==2
    Dx=Base.diff(d.domains[1])
    Dy=Base.diff(d.domains[2])    
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
    @assert size(A.ops,1)==size(B.ops,1)==1
    @assert size(A.ops,2)==size(B.ops,2)==2
    PDEOperator([A.ops[1,1]*B.ops[1,1] A.ops[1,2]*B.ops[1,2]])
end

##TODO how to determine whether x or y?
function *(a::Fun,A::PDEOperator)
    ops = copy(A.ops)
    for k=1:size(ops,1)
        ops[k,1]=a*ops[k,1]
    end
    PDEOperator(ops)
end

Base.diff(d::TensorDomain,k)=k==1?Base.diff(d.domains[1])⊗I:I⊗Base.diff(d.domains[2])  
grad(d::TensorDomain)=[Base.diff(d,k) for k=1:length(d.domains)]




function dirichlet(d::TensorDomain)
    @assert length(d.domains)==2
    Bx=dirichlet(d.domains[1])
    By=dirichlet(d.domains[2])
    [Bx⊗I,I⊗By]
end

function neumann(d::TensorDomain)
    @assert length(d.domains)==2
    Bx=neumann(d.domains[1])
    By=neumann(d.domains[2])
    [Bx⊗I,I⊗By]
end

function timedirichlet(d::TensorDomain)
    @assert length(d.domains)==2
    Bx=dirichlet(d.domains[1])
    Bt=dirichlet(d.domains[2])[1]
    [I⊗Bt,Bx⊗I]
end


function *(L::PDEOperator,f::Fun2D)
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
    
    Fun2D(A,B)
end



## Schur PDEOperator
# represent an operator that's been discretized in the Y direction

type PDEOperatorSchur{LT<:Number,MT<:Number,OT<:Number,BT<:Number,ST<:Number,FT<:Functional}
    Bx::Vector{FT}
    Lx::Operator{LT}
    Mx::Operator{MT}
    S::OperatorSchur{OT,BT} 
    
    indsBx::Vector{Int}
    indsBy::Vector{Int}
    
    Rdiags::Vector{SavedBandedOperator{ST}}
end

function PDEOperatorSchur{LT<:Number,MT<:Number,BT<:Number,ST<:Number}(Bx,Lx::Operator{LT},Mx::Operator{MT},S::OperatorSchur{BT,ST},indsBx,indsBy)
    ny=size(S,1)
    nbcs=numbcs(S)
    Rdiags=Array(SavedBandedOperator{promote_type(LT,MT)},ny)
    Xops=promotespaces([Lx,Mx])
    Lx=SavedBandedOperator(Xops[1]);Mx=SavedBandedOperator(Xops[2])
    
    resizedata!(Lx,ny);resizedata!(Mx,ny)
    
    
    for k=1:ny-nbcs
        ##TODO: Do block case
        Rdiags[k]=SavedBandedOperator(S.R[k,k]*Lx + S.T[k,k]*Mx)
        resizedata!(Rdiags[k],ny)
    end

    
    PDEOperatorSchur(Bx,Lx,Mx,S,indsBx,indsBy,Rdiags)
end


PDEOperatorSchur(Bx::Vector,Lx::Operator,Mx::Operator,S::OperatorSchur)=PDEOperatorSchur(Bx,Lx,Mx,S,[1:length(Bx)],length(Bx)+[1:numbcs(S)])
PDEOperatorSchur(Bx::Vector,Lx::Operator,Mx::UniformScaling,S::OperatorSchur)=PDEOperatorSchur(Bx,Lx,ConstantOperator(Mx.λ),S)
PDEOperatorSchur(Bx::Vector,Lx::UniformScaling,Mx::Operator,S::OperatorSchur)=PDEOperatorSchur(Bx,ConstantOperator(Lx.λ),Mx,S)


function PDEOperatorSchur(Bx,By,A::PDEOperator,ny::Integer,indsBx,indsBy)
   @assert size(A.ops)==(2,2)
   PDEOperatorSchur(Bx,A.ops[1,1],A.ops[2,1],OperatorSchur(By,A.ops[1,2],A.ops[2,2],ny),indsBx,indsBy)
end

isxfunctional(B::PDEOperator)=size(B.ops,1)==1&&size(B.ops,2)==2&&typeof(B.ops[1,1])<:Functional
isyfunctional(B::PDEOperator)=size(B.ops,1)==1&&size(B.ops,2)==2&&typeof(B.ops[1,2])<:Functional
ispdeop(B::PDEOperator)=!isxfunctional(B)&&!isyfunctional(B)


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

Base.schurfact{T<:PDEOperator}(A::Vector{T},n::Integer)=PDEOperatorSchur(A,n)


domainspace(P::PDEOperatorSchur,k::Integer)=k==1?domainspace(P.Lx):domainspace(P.S)
rangespace(P::PDEOperatorSchur,k::Integer)=k==1?rangespace(P.Lx):rangespace(P.S)
domain(P::PDEOperatorSchur,k::Integer)=k==1?domain(P.Lx):domain(P.S)

