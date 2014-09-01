
export lap,grad


type PDEOperator
    ops::Array{Operator,2}
end



PDEOperator(A,B)=PDEOperator([A B])



function domain(LL::PDEOperator,j::Integer)
    for k=1:size(LL.ops,1)
        dx=domain(LL.ops[k,j]) 
        if dx != Any
            return dx
        end
    end
    return Any
end

function domain(LL)
    @assert size(LL.ops,2)==2
    domain(LL,1)⊗domain(LL,2)
end


⊗(A::Operator,B::Operator)=PDEOperator(A,B)
⊗(A::Operator,B::UniformScaling)=A⊗ConstantOperator(1.0B.λ)
⊗(A::UniformScaling,B::Operator)=ConstantOperator(1.0A.λ)⊗B

⊗{T<:Operator}(A::Vector{T},B::Operator)=PDEOperator[PDEOperator(Ai,B) for Ai in A]
⊗{T<:Operator}(A::Operator,B::Vector{T})=PDEOperator[PDEOperator(A,Bi) for Bi in B]
⊗{T<:Operator}(A::Vector{T},B::UniformScaling)=A⊗ConstantOperator(1.0B.λ)
⊗{T<:Operator}(A::UniformScaling,B::Vector{T})=ConstantOperator(1.0A.λ)⊗B



function +(A::PDEOperator,B::PDEOperator)
    ret=copy(A.ops)
    for k=1:size(B.ops,1),j=1:size(A.ops,1)
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
    Dx=diff(d.domains[1])
    Dy=diff(d.domains[2])    
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

function grad(d::TensorDomain)
    @assert length(d.domains)==2
    Dx=diff(d.domains[1])
    Dy=diff(d.domains[2])    
    [Dx⊗I,I⊗Dy]
end





function dirichlet(d::TensorDomain)
    @assert length(d.domains)==2
    Bx=dirichlet(d.domains[1])
    By=dirichlet(d.domains[2])
    [Bx⊗I,I⊗By]
end




## Schur PDEOperator
# represent an operator that's been discretized in the Y direction

type PDEOperatorSchur{T,FT<:Functional}
    Bx::Vector{FT}
    Lx::Operator{T}
    Mx::Operator{T}
    S::OperatorSchur{T} 
    
    indsBx::Vector{Int}
    indsBy::Vector{Int}
    
    Rdiags::Vector{SavedBandedOperator{T}}
end

function PDEOperatorSchur{T}(Bx,Lx::Operator{T},Mx::Operator{T},S::OperatorSchur{T},indsBx,indsBy)
    ny=size(S,1)
    nbcs=numbcs(S)
    Rdiags=Array(SavedBandedOperator{T},ny)
    Xops=promotespaces([Lx,Mx])
    Lx=SavedBandedOperator(Xops[1]);Mx=SavedBandedOperator(Xops[2])
    
    resizedata!(Lx,ny);resizedata!(Mx,ny)
    
    
    for k=1:ny-nbcs
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

