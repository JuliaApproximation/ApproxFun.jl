export Fun,evaluate,values,points
export coefficients,integrate,differentiate,domain,space,linesum

include("Domain.jl")
include("Space.jl")

##  Constructors


type Fun{S,T}
    coefficients::Vector{T}
    space::S
    function Fun(coeff::Vector{T},sp::S)
        @assert length(coeff)≤dimension(sp)
        new(coeff,sp)
    end
end

"""
    Fun(coefficients,space)

returns a Fun with coefficients in the space
"""
Fun(coeff::Vector,sp::Space)=Fun{typeof(sp),eltype(coeff)}(coeff,sp)
Fun{T<:Integer}(coeff::Vector{T},sp::Space)=Fun(1.0coeff,sp)

function Fun(v::Vector{Any},sp::Space)
    @assert isempty(v)
    Fun(Float64[],sp)
end

##Coefficient routines
#TODO: domainscompatible?

"""
    coefficients(fun,space)

returns the coefficients of a fun in a possibly different space
"""
function coefficients(f::Fun,msp::Space)
    #zero can always be converted
    if length(f)==1 && f.coefficients[1]==0
        f.coefficients
    else
        coefficients(f.coefficients,space(f),msp)
    end
end
coefficients{T<:Space}(f::Fun,::Type{T})=coefficients(f,T(domain(f)))
coefficients(f::Fun)=f.coefficients
coefficients(c::Number,sp::Space)=Fun(c,sp).coefficients


##Convert routines


Base.convert{T,S}(::Type{Fun{S,T}},f::Fun{S})=Fun(convert(Vector{T},f.coefficients),f.space)
Base.convert{T,S}(::Type{Fun{S,T}},f::Fun)=Fun(Fun(convert(Vector{T},f.coefficients),f.space),S(domain(f)))

Base.convert{T,S}(::Type{Fun{S,T}},x::Number)=x==0?zeros(T,S(AnyDomain())):x*ones(T,S(AnyDomain()))
Base.convert{S}(::Type{Fun{S}},x::Number)=x==0?zeros(S(AnyDomain())):x*ones(S(AnyDomain()))
Base.convert{IF<:Fun}(::Type{IF},x::Number)=Fun(x)
Base.promote_rule{T,V,S}(::Type{Fun{S,T}},::Type{Fun{S,V}})=Fun{S,promote_type(T,V)}


Base.zero(::Type{Fun})=Fun(0.)
Base.zero{T,S<:Space}(::Type{Fun{S,T}})=zeros(T,S(AnyDomain()))
Base.one{T,S<:Space}(::Type{Fun{S,T}})=ones(T,S(AnyDomain()))
for op in (:(Base.zeros),:(Base.ones))
    @eval ($op){S,T}(f::Fun{S,T})=$op(T,f.space)
end

Base.zero(f::Fun)=zeros(f)
Base.one(f::Fun)=ones(f)

Base.eltype{S,T}(::Fun{S,T})=T



## domain


## General routines


domain(f::Fun)=domain(f.space)
domain{T<:Fun}(v::Vector{T})=map(domain,v)


setdomain(f::Fun,d::Domain)=Fun(f.coefficients,setdomain(space(f),d))

for op in (:tocanonical,:tocanonicalD,:fromcanonical,:fromcanonicalD,:invfromcanonicalD)
    @eval $op(f::Fun,x...)=$op(domain(f),x...)
end

for op in (:tocanonical,:tocanonicalD)
    @eval $op(d::Domain)=$op(d,Fun(identity,d))
end
for op in (:fromcanonical,:fromcanonicalD,:invfromcanonicalD)
    @eval $op(d::Domain)=$op(d,Fun(identity,canonicaldomain(d)))
end


space(f::Fun)=f.space
spacescompatible(f::Fun,g::Fun)=spacescompatible(space(f),space(g))
canonicalspace(f::Fun)=canonicalspace(space(f))
canonicaldomain(f::Fun)=canonicaldomain(domain(f))


##Evaluation

function evaluate(f::AbstractVector,S::Space,x...)
    csp=canonicalspace(S)
    if spacescompatible(csp,S)
        error("Override evaluate for " * string(typeof(csp)))
    else
        evaluate(coefficients(f,S,csp),csp,x...)
    end
end

evaluate(f::Fun,x...)=evaluate(f.coefficients,f.space,x...)


Base.call(f::Fun,x...)=evaluate(f,x...)

for op in (:(Base.first),:(Base.last))
    @eval $op{S,T}(f::Fun{S,T})=f($op(domain(f)))
end




##Data routines


values(f::Fun,dat...)=itransform(f.space,f.coefficients,dat...)
points(f::Fun)=points(f.space,length(f))
Base.length(f::Fun)=length(f.coefficients)


function Base.stride(f::Fun)
    # Check only for stride 2 at the moment
    # as higher stride is very rare anyways
    M=maxabs(f.coefficients)
    for k=2:2:length(f)
        if abs(f.coefficients[k])>40*M*eps()
            return 1
        end
    end

    2
end



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
chop!(f::Fun)=chop!(f,eps(eltype(f.coefficients)))


## Addition and multiplication

for op in (:+,:-,:(.+),:(.-))
    @eval begin
        function $op(f::Fun,g::Fun)
            if spacescompatible(f,g)
                n = max(length(f),length(g))
                f2 = pad(f,n); g2 = pad(g,n)

                Fun(($op)(f2.coefficients,g2.coefficients),isambiguous(domain(f))?g.space:f.space)
            else
                m=union(f.space,g.space)
                if isa(m,NoSpace)
                    error("Cannot "*string($op)*" because no space is the union of "*string(typeof(f.space))*" and "*string(typeof(g.space)))
                end
                $op(Fun(f,m),Fun(g,m)) # convert to same space
            end
        end
        $op{S,T<:Number}(f::Fun{S,T},c::T)=c==0?f:$op(f,Fun(c))
        $op(f::Fun,c::Number)=$op(f,Fun(c))
        $op(f::Fun,c::UniformScaling)=$op(f,c.λ)
        $op(c::UniformScaling,f::Fun)=$op(c.λ,f)
    end
end

for op in (:+,:(.+))
    @eval $op(c::Number,f::Fun)=c==0?f:$op(Fun(c),f)
end

for op in (:-,:(.-))
    @eval $op(c::Number,f::Fun)=c==0?-f:$op(Fun(c),f)
end

# equivalent to Y+=a*X
axpy!(a,X::Fun,Y::Fun)=axpy!(a,coefficients(X,space(Y)),Y)
function axpy!(a,xcfs::Vector,Y::Fun)
    if a!=0
        n=length(Y); m=length(xcfs)

        if n≤m
            resize!(Y.coefficients,m)
            for k=1:n
                @inbounds Y.coefficients[k]+=a*xcfs[k]
            end
            for k=n+1:m
                @inbounds Y.coefficients[k]=a*xcfs[k]
            end
        else #X is smaller
            for k=1:m
                @inbounds Y.coefficients[k]+=a*xcfs[k]
            end
        end
    end

    Y
end



for op = (:*,:.*,:./,:/)
    @eval $op(f::Fun,c::Number) = Fun($op(f.coefficients,c),f.space)
end

-(f::Fun)=Fun(-f.coefficients,f.space)
-(c::Number,f::Fun)=-(f-c)


for op = (:*,:.*,:+,:(.+))
    @eval $op(c::Number,f::Fun)=$op(f,c)
end


function .^{S,T}(f::Fun{S,T},k::Integer)
    if k == 0
        ones(space(f))
    elseif k==1
        f
    elseif k > 1
        f.*f.^(k-1)
    else
        1./f.^(-k)
    end
end

Base.inv{S,T}(f::Fun{S,T})=1./f

# Integrals over two Funs, which are fast with the orthogonal weight.

export dotu, linedotu, linedot

dotu(f::Fun,g::Fun)=defaultdotu(f,g)
dotu(c::Number,g::Fun)=sum(c*g)
dotu(g::Fun,c::Number)=sum(g*c)

linedotu(f::Fun,g::Fun)=defaultlinedotu(f,g)
linedotu(c::Number,g::Fun)=linesum(c*g)
linedotu(g::Fun,c::Number)=linesum(g*c)

# Having fallbacks allow for the fast implementations.

defaultdotu(f::Fun,g::Fun)=sum(f.*g)
defaultlinedotu(f::Fun,g::Fun)=linesum(f.*g)

# Conjuations

Base.dot(f::Fun,g::Fun)=dotu(conj(f),g)
Base.dot(c::Number,g::Fun)=dotu(conj(c),g)
Base.dot(g::Fun,c::Number)=dotu(conj(g),c)

linedot(f::Fun,g::Fun)=linedotu(conj(f),g)
linedot(c::Number,g::Fun)=linedotu(conj(c),g)
linedot(g::Fun,c::Number)=linedotu(conj(g),c)

## Norm

Base.norm(f::Fun) = norm(f,2)

function Base.norm(f::Fun,p::Number)
    if p < 1
        return error("p should be 1 ≤ p ≤ ∞")
    elseif 1 ≤ p < Inf
        return abs(sum(abs2(f)^(p/2)))^(1/p)
    else
        return maxabs(f)
    end
end

function Base.norm(f::Fun,p::Int)
    if 1 ≤ p < Inf
        return iseven(p) ? abs(sum(abs2(f)^div(p,2)))^(1/p) : abs(sum(abs2(f)^(p/2)))^(1/p)
    else
        return error("p should be 1 ≤ p ≤ ∞")
    end
end


## Mapped functions


for op = (:(Base.real),:(Base.imag),:(Base.conj))
    @eval ($op){T,D<:Space{RealBasis}}(f::Fun{D,T}) = Fun(($op)(f.coefficients),f.space)
end

Base.abs2{S,T<:Real}(f::Fun{S,T})=f^2
Base.abs2{S,T<:Complex}(f::Fun{S,T})=real(f)^2+imag(f)^2

##  integration

function Base.cumsum(f::Fun)
    cf = integrate(f)
    cf - first(cf)
end

Base.cumsum(f::Fun,d::Domain)=cumsum(Fun(f,d))
Base.cumsum(f::Fun,d)=cumsum(f,Domain(d))



function differentiate(f::Fun,k::Integer)
    @assert k >= 0
    (k==0)?f:differentiate(differentiate(f),k-1)
end




==(f::Fun,g::Fun) =  (f.coefficients == g.coefficients && f.space == g.space)
function Base.isapprox(f::Fun,g::Fun)
    if spacescompatible(f,g)
        m=min(length(f),length(g))
        tol=100eps()  # TODO: normalize by norm of f/g

        for k=1:m
            if !isapprox(f.coefficients[k],g.coefficients[k])
                return false
            end
        end
        for k=m+1:length(f)
            if abs(f.coefficients[k])>tol
                return false
            end
        end
        for k=m+1:length(g)
            if abs(g.coefficients[k])>tol
                return false
            end
        end

        true
    else
        sp=union(f.space,g.space)
        if isa(sp,NoSpace)
            false
        else
            isapprox(Fun(f,sp),Fun(g,sp))
        end
    end
end

Base.isapprox(f::Fun,g::Number)=isapprox(f,g*ones(space(f)))
Base.isapprox(g::Number,f::Fun)=isapprox(g*ones(space(f)),f)


Base.isreal{S,T<:Real}(f::Fun{S,T})=basistype(S)<:RealBasis
Base.isreal(f::Fun)=false



# sum, integrate, and idfferentiate are in CalculusOperator


function reverseorientation(f::Fun)
    csp=canonicalspace(f)
    if spacescompatible(csp,space(f))
        error("Implement reverseorientation for $(typeof(f))")
    else
        reverseorientation(Fun(f,csp))
    end
end


## non-vector notation

*(f::Fun,g::Fun)=f.*g
^(f::Fun,k::Integer)=f.^k
^(f::Fun,k::Union{Number,Fun})=f.^k
/(c::Union{Number,Fun},g::Fun)=c./g


include("constructors.jl")
