# DirichletSpaces


immutable ChebyshevDirichlet{left,right} <: PolynomialSpace
    domain::Union(IntervalDomain,AnyDomain)
    ChebyshevDirichlet(d)=new(d)
    ChebyshevDirichlet()=new(Interval())    
end

spacescompatible{l,r}(a::ChebyshevDirichlet{l,r},b::ChebyshevDirichlet{l,r})=domainscompatible(a,b)

ChebyshevDirichlet()=ChebyshevDirichlet{1,1}()



canonicalspace(S::ChebyshevDirichlet)=Chebyshev(domain(S))


## Dirichlet Conversion

addentries!(C::Conversion{ChebyshevDirichlet{1,0},Chebyshev},A,kr::Range1)=toeplitz_addentries!(ShiftVector([1.,1.],1),A,kr)
addentries!(C::Conversion{ChebyshevDirichlet{0,1},Chebyshev},A,kr::Range1)=toeplitz_addentries!(ShiftVector([1.,-1.],1),A,kr)
function addentries!(C::Conversion{ChebyshevDirichlet{1,1},Chebyshev},A,kr::Range)
    A=toeplitz_addentries!(ShiftVector([1.,0.,-1.],1),A,kr)    
    if kr[1]==1
        A[1,2]+=-1
    end
    A
end
function addentries!(C::Conversion{ChebyshevDirichlet{2,2},Chebyshev},A,kr::Range)
    for k=kr
        A[k,k]=1
        A[k,k+4]=2*(k+1)/k-1
        if k>= 3
            A[k,k+2]=-2*(k-1)/(k-2)
        end
    end
    
    A
end
bandinds(::Conversion{ChebyshevDirichlet{1,0},Chebyshev})=0,1
bandinds(::Conversion{ChebyshevDirichlet{0,1},Chebyshev})=0,1
bandinds(::Conversion{ChebyshevDirichlet{1,1},Chebyshev})=0,2

conversion_rule(b::ChebyshevDirichlet,a::Chebyshev)=b

# return the space that has banded Conversion to the other
# function conversion_rule(a::ChebyshevDirichlet,b::Ultraspherical)
#     @assert domainscompatible(a,b)
#     
#     a
# end




## Evaluation Functional


function getindex(B::Evaluation{ChebyshevDirichlet{1,0},Bool},kr::Range)
    d = domain(B)
    
    if B.x == false && B.order == 0
        Float64[k==1?1.0:0.0 for k=kr]
    elseif B.x == true && B.order == 0
        Float64[k==1?1.0:2.0 for k=kr]
    else
        getindex(Evaluation(d,B.x,B.order)*Conversion(domainspace(B)),kr)
    end
end

function getindex(B::Evaluation{ChebyshevDirichlet{0,1},Bool},kr::Range)
    d = domain(B)
    
    if B.x == true && B.order == 0
        Float64[k==1?1.0:0.0 for k=kr]
    elseif B.x == false && B.order == 0
        Float64[k==1?1.0:-(-1)^k*2.0 for k=kr]        
    else
        getindex(Evaluation(d,B.x,B.order)*Conversion(domainspace(B)),kr)
    end
end

function getindex(B::Evaluation{ChebyshevDirichlet{1,1},Bool},kr::Range)
   tol = 200.*eps()
    d = domain(B)
    
    if B.x == false && B.order == 0
        Float64[k==1?1.0:(k==2?-2.0:0.0) for k=kr]
    elseif B.x == true && B.order == 0
        Float64[k==1?1.0:0.0 for k=kr]
    else
        getindex(Evaluation(d,B.x,B.order)*Conversion(domainspace(B)),kr)
    end
end

Evaluation(sp::ChebyshevDirichlet,x::Float64,ord::Integer)=EvaluationWrapper(sp,x,ord,Evaluation(domain(sp),x,ord)*Conversion(sp))
