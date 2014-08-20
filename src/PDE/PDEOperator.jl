type PDEOperator
    ops::Array{Operator,2}
end



PDEOperator(A,B)=PDEOperator([A B])
⊗(A::Operator,B::Operator)=PDEOperator(A,B)
⊗(A::Operator,B::UniformScaling)=PDEOperator(A,ConstantOperator(1.0B.λ))
⊗(A::UniformScaling,B::Operator)=PDEOperator(ConstantOperator(1.0A.λ),B)
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

+(A::UniformScaling,B::PDEOperator)=B+ConstantOperator(1.0)⊗ConstantOperator(1.0)
+(B::PDEOperator,A::UniformScaling)=B+ConstantOperator(1.0)⊗ConstantOperator(1.0)

function grad(d::TensorDomain)
    @assert length(d.domains)==2
    Dx=diff(d.domains[1])
    Dy=diff(d.domains[2])    
    [Dx⊗I,I⊗Dy]
end