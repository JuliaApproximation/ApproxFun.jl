export DirichletConversionOperator

## ConversionOperator

type DirichletConversionOperator <: BandedOperator{Float64}
    left::Int
    right::Int
end

domainspace(M::DirichletConversionOperator)=UltrasphericalDirchletSpace(0,Any,1,0)
rangespace(M::DirichletConversionOperator)=UltrasphericalSpace(0)

bandinds(C::DirichletConversionOperator)=0,(C.left+C.right)

function addentries!(C::DirichletConversionOperator,A::ShiftArray,kr::Range1)
    if C.left == 1 &&  C.right == 0
        toeplitz_addentries!(ShiftVector([1.,1.],1),A,kr)
    elseif  C.left == 0 &&  C.right == 1
        toeplitz_addentries!(ShiftVector([1.,-1.],1),A,kr)
    elseif C.left == C.right == 1
        toeplitz_addentries!(ShiftVector([1.,0.,-1.],1),A,kr)    
    else
        error("Higher order Dirichlet not implmented")
    end
end




