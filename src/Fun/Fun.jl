

include("Domain.jl")
include("FunctionSpace.jl")


##  Constructors



##TODO: No zero length funs
type Fun{T<:Union(Float64,Complex{Float64}),S<:FunctionSpace} 
    coefficients::Vector{T}
    space::S
end

##Coefficient routines
#TODO: domainscompatible?
coefficients(f::Fun,msp::FunctionSpace)=spaceconversion(f.coefficients,space(f),msp)
coefficients{T<:FunctionSpace}(f::Fun,::Type{T})=coefficients(f,T(AnyDomain()))
coefficients(f::Fun)=coefficients(f,Space(domain(f)))  #TODO: coefficients() returns canonical coefficients.  is this confusing??

##Convert routines


Base.convert{T<:Number,S<:DomainSpace}(::Type{Fun{T,S}},x::Number)=x*ones(T,S(AnyDomain()))
Base.convert{T<:Number,S<:FunctionSpace}(::Type{Fun{Complex{Float64},S}},f::Fun{T,S})=Fun(convert(Vector{Complex{Float64}},f.coefficients),f.space)
Base.promote_rule{T<:Number,S<:FunctionSpace}(::Type{Fun{Complex{Float64},S}},::Type{Fun{T,S}})=Fun{Complex{Float64},S}
Base.promote_rule{T<:Number,IF<:Fun}(::Type{IF},::Type{T})=IF


Base.zero{T,S<:DomainSpace}(::Type{Fun{T,S}})=zeros(T,S(AnyDomain()))
Base.one{T,S<:DomainSpace}(::Type{Fun{T,S}})=ones(T,S(AnyDomain()))
for op in (:(Base.zeros),:(Base.ones))
    @eval ($op){T}(f::Fun{T})=$op(T,f.space)
end



## domain


## General routines

domain(f::Fun)=domain(f.space)
#domain(f::FFun)=f.domain
domain(::Number)=Any
domain{T<:Fun}(v::Vector{T})=map(domain,v)


for op = (:tocanonical,:tocanonicalD,:fromcanonical,:fromcanonicalD)
    @eval ($op)(f::Fun,x)=($op)(domain(f),x)
end

space(f::Fun)=f.space
spacescompatible(f::Fun,g::Fun)=spacescompatible(space(f),space(g))



##Evaluation

evaluate{T,S}(f::Fun{T,S},x)=evaluate(Fun(f,domain(f)),x)  #Default is convert to Chebyshev
Base.getindex(f::Fun,x)=evaluate(f,x)

for op in (:(Base.first),:(Base.last))
    @eval $op{T,S}(f::Fun{T,S})=f[$op(domain(f))]
end




##Data routines
values(f::Fun)=itransform(f.space,f.coefficients) 
points(f::Fun)=points(f.space,length(f))
Base.length(f::Fun)=length(f.coefficients)



## Manipulate length


pad!(f::Fun,n::Integer)=pad!(f.coefficients,n)
pad(f::Fun,n::Integer)=Fun(pad(f.coefficients,n),f.space)


function chop!{T}(f::Fun{T},tol::Real)
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
                $op(Fun(f,domain(f)),Fun(g,domain(g))) # convert to Chebyshev
            end
        end

        ($op){N<:Number,T<:Number}(f::Fun{T},c::N)=$op(f,c*ones(f))
        ($op){N<:Number,T<:Number}(c::N,f::Fun{T})=$op(c*ones(f),f)    
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

import Base.norm

norm(f::Fun)=real(sqrt(sum(f.*conj(f))))



## Mapped functions

import Base.imag, Base.real, Base.conj

for op = (:real,:imag,:conj) 
    @eval ($op){T,D<:IntervalDomainSpace}(f::Fun{T,D}) = Fun(($op)(f.coefficients),f.space)
end

Base.abs2(f::Fun{Float64})=f.^2
Base.abs2(f::Fun{Complex{Float64}})=real(f).^2+imag(f).^2

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

function .*{T,N,S}(f::Fun{T,S},g::Fun{N,S})
    @assert domainscompatible(f,g)
    #TODO Coefficient space version
    n = length(f) + length(g) - 1
    f2 = pad(f,n); g2 = pad(g,n)
    
    sp=space(f)
    chop!(Fun(transform(sp,values(f2).*values(g2)),sp),10eps())
end

# When the spaces differ we promote and multiply
function .*{T,N,S,V}(f::Fun{T,S},g::Fun{N,V})
    sp=minspace(space(f),space(g))
    Fun(f,sp).*Fun(g,sp)
end


function Base.sum{T}(f::Fun{T})
    cf=integrate(f)
    last(cf) - first(cf)
end


integrate{T,D}(f::Fun{T,D})=integrate(Fun(f,domain(f)))



include("constructors.jl")

