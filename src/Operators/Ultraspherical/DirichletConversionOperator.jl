export DirichletConversionOperator

## ConversionOperator

type DirichletConversionOperator <: BandedOperator{Float64}
    left::Int
    right::Int
end

DirichletConversionOperator(B::ChebyshevDirichletSpace)=DirichletConversionOperator(B.left,B.right)

domainspace(M::DirichletConversionOperator)=ChebyshevDirichletSpace(AnyDomain(),M.left,M.right)
rangespace(M::DirichletConversionOperator)=UltrasphericalSpace(0)

bandinds(C::DirichletConversionOperator)=0,(C.left+C.right)

function addentries!(C::DirichletConversionOperator,A::ShiftArray,kr::Range1)
    if C.left == 1 &&  C.right == 0
        toeplitz_addentries!(ShiftVector([1.,1.],1),A,kr)
    elseif  C.left == 0 &&  C.right == 1
        toeplitz_addentries!(ShiftVector([1.,-1.],1),A,kr)
    elseif C.left == C.right == 1
        toeplitz_addentries!(ShiftVector([1.,0.,-1.],1),A,kr)    
    elseif C.left == C.right == 2
        for k=kr
            A[k,0]=1
            A[k,4]=2*(k+1)/k-1
            if k>= 3
                A[k,2]=-2*(k-1)/(k-2)
            end
        end
        
        A
    else
        error("Higher order Dirichlet not implmented")
    end
end

function ConversionOperator(B::ChebyshevDirichletSpace,A::UltrasphericalSpace)
    if A.order == 0
        DirichletConversionOperator(B)
    else
        ConversionOperator(UltrasphericalSpace(0),A)*DirichletConversionOperator(B)
    end
end


