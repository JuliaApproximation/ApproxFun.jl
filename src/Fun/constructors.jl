Fun{T<:Union(Int64,Complex{Int64})}(coefs::Vector{T},d::FunctionSpace)=Fun(1.0*coefs,d)


function Fun(f::Function,d::FunctionSpace,n::Integer)
    pts=points(d,n)
    f1=f(pts[1])
    

    if isa(f1,Vector) && !isa(d,VectorFunctionSpace) 
        return Fun(f,VectorFunctionSpace(d,length(f1)),n)
    end
    
        
    T=typeof(f1)
        
    vals=T[f(x) for x in pts]
    Fun(transform(d,vals),d)
end

# the following is to avoid ambiguity
# Fun(f::Fun,d) should be equivalent to Fun(x->f[x],d)
#TODO: fall back to Fun(x->f[x],d) if conversion not implemented?
Fun(f::Fun,d::FunctionSpace)=Fun(coefficients(f,d),d)
Fun{T<:FunctionSpace}(f::Fun,::Type{T})=Fun(f,T(domain(f)))
Fun{T<:FunctionSpace}(c::Number,::Type{T})=Fun(c,T(AnyDomain()))



Fun{T<:FunctionSpace}(f,::Type{T})=Fun(f,T())
Fun{T<:FunctionSpace}(f,::Type{T},n::Integer)=Fun(f,T(),n)

Fun(f,d::Domain)=Fun(f,Space(d))
Fun(f,d::Domain,n)=Fun(f,Space(d),n)
Fun{T<:Domain}(f,::Type{T})=Fun(f,T())


Fun(c::Number)=Fun([c])

# We do zero special since zero exists even when one doesn'
Fun{T<:FunctionSpace}(c::Number,::Type{T})=c==0?zeros(T(AnyDomain())):c*ones(T(AnyDomain()))
Fun(c::Number,d::Domain)=c==0?zeros(d):c*ones(d)
Fun(c::Number,d::FunctionSpace)=c==0?zeros(d):c*ones(d)
Fun(c::Number,n::Integer)=Fun([c],n)

## List constructor

Fun{T<:Domain}(c::Number,dl::Vector{T})=Fun(c,UnionDomain(dl))
Fun{T<:Domain}(f,dl::Vector{T})=Fun(f,UnionDomain(dl))
Fun{T<:Domain}(f,dl::Vector{T},n::Integer)=Fun(f,UnionDomain(dl),n)

## Adaptive constructors

function randomFun(f::Function,d::IntervalDomain)
    @assert d == Interval()

    #TODO: implement other domains
    
    Fun(chebyshevtransform(randomadaptivebary(f)),d)
end


# function veczerocfsFun(f::Function,d::IntervalDomain)
#     #reuse function values
# 
#     tol = 200*eps()
# 
#     for logn = 4:20
#         cf = Fun(f, d, 2^logn + 1)
#         cfs=coefficients(cf)  ##TODO: general domain
#         
#         if norm(cfs[:,end-8:end],Inf) < tol*norm(cfs[:,1:8],Inf)
#             nrm=norm(cfs,Inf)
#             return map!(g->chop!(g,10eps()*nrm),cf)
#         end
#     end
#     
#     warn("Maximum length reached")
#     
#     Fun(f,d,2^21 + 1)
# end

function zerocfsFun(f::Function,d::FunctionSpace)
    #TODO: reuse function values?
    f0=f(first(domain(d)))

    if !isa(d,VectorFunctionSpace) && isa(f0,Vector)       
        return zerocfsFun(f,VectorFunctionSpace(d,length(f0)))
    end

    tol = 200*eps()

    r=rand(d)
    fr=f(r)

    for logn = 4:20
        cf = Fun(f, d, 2^logn + 1)
        absc=abs(cf.coefficients)
        maxabsc=maximum(absc)
        if maxabsc==0 && fr==0
            return(zeros(d))
        end
        
        # we allow for transformed coefficients being a different size
        ##TODO: how to do scaling for unnormalized bases like Jacobi?
        if length(cf) > 8 && maximum(absc[end-8:end]) < tol*maxabsc &&  norm(cf[r]-fr)<1E-4  
            return chop!(cf,10eps()*maxabsc)
        end
    end
    
    warn("Maximum length reached")
    
    Fun(f,d,2^21 + 1)
end




function abszerocfsFun(f::Function,d::FunctionSpace)
    #reuse function values

    tol = 200eps();

    for logn = 4:20
        cf = Fun(f, d, 2^logn + 1)
        
        if maximum(abs(cf.coefficients[end-8:end])) < tol
            return chop!(cf,10eps())
        end
    end
    
    warn("Maximum length reached")
    
    Fun(f,d,2^21 + 1)
end


function Fun(f::Function, d::FunctionSpace; method="zerocoefficients")
    if f==identity
        identity_fun(d)
    elseif f==zero # zero is always defined
        zeros(Float64,d)
    elseif f==one
        ones(Float64,d)
    elseif method == "zerocoefficients"
        zerocfsFun(f,d)
    elseif method == "abszerocoefficients"
        abszerocfsFun(f,d)
    else
        randomFun(f,d)    
    end
end
Fun(f::Function,d::Domain;opts...)=Fun(f,Space(d);opts...)





## Aliases


FFun(x,d::PeriodicDomain)=Fun(x,Laurent(d))
FFun(x,d::PeriodicDomain,n...)=Fun(x,Laurent(d),n...)
FFun{T<:Number}(x,d::Vector{T})=Fun(x,Laurent(d))
FFun{T<:Number}(x,d::Vector{T},n...)=Fun(x,Laurent(d),n...)
FFun(f,n::Integer)=Fun(f,Laurent(PeriodicInterval()),n)
FFun(f)=Fun(f,Laurent(PeriodicInterval()))



Fun(f::Function,n::Integer)=Fun(f,Interval(),n)
Fun{T<:Number}(f::Function,d::Vector{T},n::Integer)=Fun(f,Interval(d),n)
Fun{T<:Number}(cfs::Vector{T})=Fun(1.0*cfs,Interval())
Fun{T<:Number,M<:Number}(cfs::Vector{M},d::Vector{T})=Fun(1.0*cfs,Interval(d))
Fun{T<:Number}(f::Function,d::Vector{T})=Fun(f,Interval(d))



function Fun(cfs::Vector{Any})
    @assert isempty(cfs)
    Fun(Float64[])
end
function Fun(cfs::Vector{Any},s::FunctionSpace)
    @assert isempty(cfs)
    Fun(Float64[],s)
end