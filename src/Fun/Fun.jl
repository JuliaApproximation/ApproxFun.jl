

include("Domain.jl")
include("FunctionSpace.jl")


##  Constructors



##TODO: No zero length funs
type Fun{S<:FunctionSpace,T<:Union(Float64,Complex{Float64})} 
    coefficients::Vector{T}
    space::S
end

##Coefficient routines
#TODO: domainscompatible?
coefficients(f::Fun,msp::FunctionSpace)=spaceconversion(f.coefficients,space(f),msp)
coefficients{T<:FunctionSpace}(f::Fun,::Type{T})=coefficients(f,T(AnyDomain()))
canonicalcoefficients(f::Fun)=coefficients(f,canonicalspace(f.space))  
coefficients(f::Fun)=f.coefficients
##Convert routines


Base.convert{T<:Number,S<:DomainSpace}(::Type{Fun{S,T}},x::Number)=x*ones(T,S(AnyDomain()))
Base.convert{T<:Number,S<:FunctionSpace}(::Type{Fun{S,Complex{Float64}}},f::Fun{S,T})=Fun(convert(Vector{Complex{Float64}},f.coefficients),f.space)
Base.promote_rule{T<:Number,S<:FunctionSpace}(::Type{Fun{S,Complex{Float64}}},::Type{Fun{S,T}})=Fun{S,Complex{Float64}}
Base.promote_rule{T<:Number,IF<:Fun}(::Type{IF},::Type{T})=IF


Base.zero{T,S<:DomainSpace}(::Type{Fun{S,T}})=zeros(T,S(AnyDomain()))
Base.one{T,S<:DomainSpace}(::Type{Fun{S,T}})=ones(T,S(AnyDomain()))
for op in (:(Base.zeros),:(Base.ones))
    @eval ($op){S,T}(f::Fun{S,T})=$op(T,f.space)
end



## domain


## General routines

domain(f::Fun)=domain(f.space)
#domain(f::FFun)=f.domain
domain{T<:Fun}(v::Vector{T})=map(domain,v)


for op = (:tocanonical,:tocanonicalD,:fromcanonical,:fromcanonicalD)
    @eval ($op)(f::Fun,x)=($op)(domain(f),x)
end

space(f::Fun)=f.space
spacescompatible(f::Fun,g::Fun)=spacescompatible(space(f),space(g))



##Evaluation

evaluate{S,T}(f::Fun{S,T},x)=evaluate(Fun(f,domain(f)),x)  #Default is convert to Chebyshev
Base.getindex(f::Fun,x)=evaluate(f,x)

for op in (:(Base.first),:(Base.last))
    @eval $op{S,T}(f::Fun{S,T})=f[$op(domain(f))]
end




##Data routines
values(f::Fun,dat...)=itransform(f.space,f.coefficients,dat...) 
points(f::Fun)=points(f.space,length(f))
Base.length(f::Fun)=length(f.coefficients)



## Manipulate length


pad!(f::Fun,n::Integer)=pad!(f.coefficients,n)
pad(f::Fun,n::Integer)=Fun(pad(f.coefficients,n),f.space)


function chop!{S,T}(f::Fun{S,T},tol::Real)
    chop!(f.coefficients,tol)
    if length(f.coefficients) == 0
        f.coefficients = [zero(T)]
    end
    
    f
end
chop(f::Fun,tol)=chop!(Fun(copy(f.coefficients),f.space),tol)
chop!(f::Fun)=chop!(f,eps())


## Addition and multiplication




for op = (:+,:-)
    @eval begin
        function ($op)(f::Fun,g::Fun)
            if spacescompatible(f,g)
                n = max(length(f),length(g))
                f2 = pad(f,n); g2 = pad(g,n)
            
                Fun(($op)(f2.coefficients,g2.coefficients),domain(f)!=AnyDomain()?f.space:g.space)
            else 
                #TODO: is it better to convert to minspace?
                m=maxspace(f.space,g.space)
                $op(Fun(f,m),Fun(g,m)) # convert to same space
            end
        end

        ($op){N<:Number}(f::Fun,c::N)=$op(f,c*ones(f))
        ($op){N<:Number}(c::N,f::Fun)=$op(c*ones(f),f)    
    end
end 



fasttimes(f2,g2)=Fun(chebyshevtransform(values(f2).*values(g2)),domain(f2))




for op = (:*,:.*,:./,:/)
    @eval ($op)(f::Fun,c::Number) = Fun(($op)(f.coefficients,c),f.space)
end 

-(f::Fun)=Fun(-f.coefficients,f.space)
-(c::Number,f::Fun)=-(f-c)


for op = (:*,:.*,:+)
    @eval ($op)(c::Number,f::Fun)=($op)(f,c)
end




function .^(f::Fun,k::Integer)
    if k == 0
        1.
    elseif k > 0
        f.*f.^(k-1)
    else
        f./f.^(k+1)
    end
end



## Norm


Base.dot(f::Fun,g::Fun)=sum(conj(f).*g)
Base.norm(f::Fun)=real(sqrt(dot(f,f)))



## Mapped functions

import Base.imag, Base.real, Base.conj

for op = (:real,:imag,:conj) 
    @eval ($op){T,D<:DomainSpace{Float64}}(f::Fun{D,T}) = Fun(($op)(f.coefficients),f.space)
end

Base.abs2{S}(f::Fun{S,Float64})=f.^2
Base.abs2{S}(f::Fun{S,Complex{Float64}})=real(f).^2+imag(f).^2

##  integration

function Base.cumsum(f::Fun)
    cf = integrate(f)
    cf - first(cf)
end



function differentiate(f::Fun,k::Integer)
    @assert k >= 0
    (k==0)?f:differentiate(differentiate(f),k-1)
end

Base.diff(f::Fun,n...)=differentiate(f,n...)



==(f::Fun,g::Fun) =  (f.coefficients == g.coefficients && f.space == g.space)





## Overrideable


# multiplying 2 Funs, we assume this can be done by transform
# the parametrizations are to make it the broadest definition

function .*{T,N,S}(f::Fun{S,T},g::Fun{S,N})
    @assert domainscompatible(f,g)
    #TODO Coefficient space version
    n = length(f) + length(g) - 1
    f2 = pad(f,n); g2 = pad(g,n)
    
    sp=space(f)
    chop!(Fun(transform(sp,values(f2).*values(g2)),sp),10eps())
end

# When the spaces differ we promote and multiply
function .*{T,N,S,V}(f::Fun{S,T},g::Fun{V,N})
    sp=minspace(space(f),space(g))
    if sp==NoSpace() # see if a multiplication operator is implemented
        Multiplication(f,space(g))*g
    else
        Fun(f,sp).*Fun(g,sp)
    end
end


function Base.sum{S,T}(f::Fun{S,T})
    cf=integrate(f)
    last(cf) - first(cf)
end


integrate{D,T}(f::Fun{D,T})=integrate(Fun(f,domain(f)))



include("constructors.jl")

