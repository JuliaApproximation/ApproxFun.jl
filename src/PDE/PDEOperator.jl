
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

