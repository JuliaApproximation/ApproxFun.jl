
## This allows "dekron"ing a PDE operator into two ODE operators
# used in schurfact

# productop is a product of two operator


typealias WrapperOps Union{ConversionWrapper,MultiplicationWrapper,DerivativeWrapper,IntegralWrapper,LaplacianWrapper}


isproductop(a) = iskronop(a)  # all kron ops are product ops


iskronop(::) = false
iskronop(::KroneckerOperator) = true


iskronop(A::Union{WrapperOps,SpaceOperator,ConstantTimesOperator})=iskronop(A.op)
iskronop(A::TimesOperator)=all(iskronop,A.ops)
iskronop{T,DS<:TensorSpace}(::ConstantOperator{T,DS}) = true
iskronop{T,S<:TensorSpace,V<:TensorSpace}(::ZeroOperator{T,S,V}) = true


iskronsumop(::)=false
iskronsumop(A::PlusOperator)=all(a->iskronop(a)||iskronsumop(a),A.ops)
iskronsumop(A::Union{WrapperOps,SpaceOperator,ConstantTimesOperator})=iskronsumop(A.op)


#dekron gives the operators that make up a productop

dekron(K::KroneckerOperator,k)=K.ops[k]

dekron(S::WrapperOps,k)=dekron(S.op,k)
dekron(S::TimesOperator,k)=TimesOperator(Operator{eltype(eltype(S))}[dekron(op,k) for op in S.ops])
dekron(S::SpaceOperator,k)=SpaceOperator(dekron(S.op,k),domainspace(S)[k],rangespace(S)[k])
#TODO: dekron(S::SpaceOperator,k)=SpaceOperator(dekron(S.op,k),domainspace(S)[k],rangespace(S)[k])
dekron(sp::ConstantTimesOperator,k)=k==1?sp.λ*dekron(sp.op,k):dekron(sp.op,k)
dekron{T,DS<:TensorSpace}(C::ConstantOperator{T,DS},k) =
    k==1?ConstantOperator(C.λ,C.space[1]):ConstantOperator(one(C.λ),C.space[k])
dekron{T,S<:TensorSpace,V<:TensorSpace}(C::ZeroOperator{T,S,V},k) =
    ZeroOperator(C.domainspace[k],C.rangespace[k])


dekron(K)=dekron(K,1),dekron(K,2)


sumops(A)=A.ops
sumops(A::Union{SpaceOperator,WrapperOps})=sumops(A.op)
sumops(A::ConstantTimesOperator)=Operator{eltype(A)}[A.λ*op for op in sumops(A.op)]
function dekron(T::Type,A,k,::Colon)
    ret=Array(Operator{T},0)
    if iskronop(A)
        push!(ret,dekron(A,k))
    elseif iskronsumop(A)
        for op in sumops(A)
            da=dekron(op,k,:)
            for a in da
                push!(ret,a)
            end
        end
    else
       error("dekron not defined unless iskronop or iskronsumop" )
    end
    ret
end

dekron(A,k,::Colon)=dekron(eltype(eltype(A)),A,k,:)

dekron(A,::Colon,::Colon)=dekron(A,1,:),dekron(A,2,:)

function reducekronsum!(opsx,opsy)
    @assert length(opsx)==length(opsy)
    k=1
    while k ≤ length(opsx)
        for j=length(opsx):-1:k+1
            if opsx[k]==opsx[j]
               opsy[k]+=opsy[j]
               deleteat!(opsx,j)
               deleteat!(opsy,j)
           end
         end
        k+=1
    end
    opsx,opsy
end

function simplifydekron(A)
    opsx,opsy=dekron(A,:,:)
    map!(simplify,opsx);map!(simplify,opsy)
    reducekronsum!(opsy,opsx)
    reducekronsum!(opsx,opsy)
end
