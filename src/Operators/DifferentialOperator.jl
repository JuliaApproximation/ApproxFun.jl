export DifferentialOperator
export conversionmatrix,derivativematrix
export multiplicationmatrix

## DifferentialOperator constructors

type DifferentialOperator{T<:IFun} <: InfiniteOperator
    coefficients::Vector{T}
    domain::IntervalDomain
end

##TODO: ensure any funs match coefficients

DifferentialOperator(cfs::Vector)=DifferentialOperator(cfs,Interval())
DifferentialOperator(cfs,d::Vector)=DifferentialOperator(cfs,Interval(d))
DifferentialOperator{T<:Real}(cfs::Vector{T},d::IntervalDomain)=DifferentialOperator(IFun([1.],d)*cfs,d)
DifferentialOperator(cfs::Vector{Any},d::IntervalDomain)=DifferentialOperator(IFun([1.],d)*cfs,d)



Base.diff(d::IntervalDomain,n::Integer)=DifferentialOperator([zeros(n),1.],d)
Base.diff(d::IntervalDomain)=Base.diff(d,1)
Base.eye(d::IntervalDomain)=DifferentialOperator([1.],d)



##Information about operator



differentialorder(d::DifferentialOperator)=length(d.coefficients)-1

##TODO: Should band always include 0?

function subbandindex(d::DifferentialOperator)
    ret = 0
    
    for k = 1:length(d.coefficients)
        ret = min(ret,k - length(d.coefficients[k]))
    end
    
    ret
end

function superbandindex(d::DifferentialOperator)
    od=differentialorder(d)
  
    ret = 0
    
    for k = 1:length(d.coefficients)
        if (norm(d.coefficients[k].coefficients) > eps())
            ret = max(ret,2(od - k) + k + length(d.coefficients[k]))
        end
    end
    
    ret
end

bandrange(d::DifferentialOperator)=subbandindex(d):superbandindex(d)




function derivativematrix(μ::Integer,d::Domain,m::Integer)
    if(μ == 0)    
        speye(m,m)
    else
        ret = spzeros(m,m)

        for l=1:m-μ
            ret[l, l + μ] =  2.^(μ-1).*factorial(μ-1).*(2./(d.b-d.a)).^μ.*(μ+l-1)
        end
        
        ret
    end
end

function conversionmatrix(λ::Integer,m::Integer)
    ret = spzeros(m,m)

    if λ == 0
        for l=1:m
            ret[l, l] =  l==1? 1. : .5
        end
        
        for l=1:m-2
            ret[l, l + 2] =  -.5
        end
        
        ret    
    else
        for l=1:m
            ret[l, l] =  λ./(l - 1 + λ)
        end
        
        for l=1:m-2
            ret[l, l + 2] =  -λ./(l + λ + 1)
        end
        
        ret    
    end
end

function conversionmatrix(od::Range1,m::Integer)        
    ret = speye(m,m)
    
    for λ = od[1]:od[end]-1
       ret = conversionmatrix(λ,m)*ret 
    end
    
    ret
end

function derivativematrix(μ::Integer,od::Range1,d::Domain,m::Integer)
    @assert od[1] == 0  #TODO: more general
    
    conversionmatrix(μ:od[end],m)*derivativematrix(μ,d,m)
end


function Base.getindex(op::DifferentialOperator,k::Range1,j::Range1)
    m = max(k[end],j[end])  #TODO: smart indexing
    
    ret = spzeros(m,m)
    
    od = differentialorder(op)
    
    for l = 1:length(op.coefficients)
        ret += conversionmatrix(l-1:od,m)*multiplicationmatrix(l-1,op.coefficients[l],m)*derivativematrix(l-1,op.domain,m)
    end
    
    ret[k,j]
end



## Multiplication matrix



multiplicationmatrix(f::IFun,n::Integer)=.5(ToeplitzOperator(f)[1:n,1:n] + [spzeros(1,n),HankelOperator(f.coefficients[2:end])[1:n-1,1:n]])

umultiplicationmatrix(f::IFun,n::Integer)=.5(ToeplitzOperator(f)[1:n,1:n]  - HankelOperator(f.coefficients[3:end])[1:n,1:n])


multiplicationmatrix(a::Number,n::Integer)=a*speye(n)

multiplicationmatrix(l::Integer,a::Number,n)=a*speye(n)

function multiplicationmatrix(l::Integer,a::IFun,n)
    if length(a) == 1
        multiplicationmatrix(a.coefficients[1],n)
    elseif l==0
        multiplicationmatrix(a,n)
    elseif l==1
        umultiplicationmatrix(a,n)
    elseif 
        error("higher order multiplication not yet implemented")
    end
end

## ultraconv

ultraconversion(v::Vector{Float64},μ::Integer)=conversionmatrix(0:μ,length(v))*v
ultraiconversion(v::Vector{Float64},μ::Integer)=conversionmatrix(0:μ,length(v))\v
coefficients(f::IFun,m::Integer)=ultraconversion(f.coefficients,m)



## Addidition of Differential operators

function +(A::DifferentialOperator,B::DifferentialOperator)
    @assert A.domain == B.domain
    
    n = max(length(A.coefficients),length(B.coefficients))
    
    DifferentialOperator(pad(A.coefficients,n) + pad(B.coefficients,n),A.domain)
end

-(A::DifferentialOperator)=DifferentialOperator(-A.coefficients,A.domain)

function -(A::DifferentialOperator,B::DifferentialOperator)
    @assert A.domain == B.domain
    
    n = max(length(A.coefficients),length(B.coefficients))
    
    DifferentialOperator(pad(A.coefficients,n) - pad(B.coefficients,n),A.domain)
end

*(a::Number,B::DifferentialOperator)=DifferentialOperator(a*B.coefficients,B.domain)
*(a::IFun,B::DifferentialOperator)=DifferentialOperator(a*B.coefficients,B.domain)


function +(a::IFun, B::DifferentialOperator)
    @assert a.domain == B.domain

    cfs = deepcopy(B.coefficients)
    cfs[1] = cfs[1] + a
    DifferentialOperator(cfs,B.domain)
end

+(B::DifferentialOperator,a::IFun)=a+B
-(B::DifferentialOperator,a::IFun)=-a + B



