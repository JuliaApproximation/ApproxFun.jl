export Fun,evaluate,values,points,extrapolate
export coefficients,ncoefficients
export integrate,differentiate,domain,space,linesum,linenorm

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
Fun(coeff::Vector,sp::Space) = Fun{typeof(sp),eltype(coeff)}(coeff,sp)
Fun{T<:Integer}(coeff::Vector{T},sp::Space) = Fun(1.0coeff,sp)

function Fun(v::Vector{Any},sp::Space)
    if isempty(v)  || all(x->isa(x,Number) && x==0,v)
        Fun(Float64[],sp)
    else
        error("Cannot convert $v to a Fun of type $sp")
    end
end

##Coefficient routines
#TODO: domainscompatible?

"""
    coefficients(fun,space)

returns the coefficients of a fun in a possibly different space
"""
function coefficients(f::Fun,msp::Space)
    #zero can always be converted
    if ncoefficients(f)==1 && f.coefficients[1]==0
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
Base.convert{T,S}(::Type{Fun{S,T}},f::Fun)=Fun(Fun(convert(Vector{T},f.coefficients),f.space),S(domain(f)))  #TODO: this line is incompatible with space conversion

Base.convert{T,S}(::Type{Fun{S,T}},x::Number) =
    x==0?zeros(T,S(AnyDomain())):x*ones(T,S(AnyDomain()))
Base.convert{S}(::Type{Fun{S}},x::Number) =
    x==0?zeros(S(AnyDomain())):x*ones(S(AnyDomain()))
Base.convert{IF<:Fun}(::Type{IF},x::Number)=Fun(x)
Base.promote_rule{T,V,S}(::Type{Fun{S,T}},::Type{Fun{S,V}})=Fun{S,promote_type(T,V)}


# promotion of * to fix 0.5 bug
if VERSION ≥ v"0.5.0-rc1+1"
    Base.promote_op{N,V,S,T}(::typeof(*),::Type{Fun{N,V}},::Type{Fun{S,T}}) =
        Fun{promote_type(N,S),promote_type(T,V)}
    Base.promote_op{N,S,T}(::typeof(*),::Type{N},::Type{Fun{S,T}}) = Fun{S,promote_type(N,T)}
    Base.promote_op{N,S,T}(::typeof(*),::Type{Matrix{N}},::Type{Matrix{Fun{S,T}}}) =
        Matrix{Fun{S,promote_type(N,T)}}
end


Base.zero(::Type{Fun})=Fun(0.)
Base.zero{T,S<:Space}(::Type{Fun{S,T}})=zeros(T,S(AnyDomain()))
Base.one{T,S<:Space}(::Type{Fun{S,T}})=ones(T,S(AnyDomain()))
for op in (:(Base.zeros),:(Base.ones))
    @eval ($op){S,T}(f::Fun{S,T})=$op(T,f.space)
end

Base.zero(f::Fun)=zeros(f)
Base.one(f::Fun)=ones(f)

Base.eltype{S,T}(::Fun{S,T})=T




setspace(v::AbstractVector,s::Space) = Fun(v,s)
setspace(f::Fun,s::Space) = Fun(f.coefficients,s)


## domain


## General routines


domain(f::Fun) = domain(f.space)
domain{T<:Fun}(v::AbstractMatrix{T}) = map(domain,v)


setdomain(f::Fun,d::Domain) = Fun(f.coefficients,setdomain(space(f),d))

for op in (:tocanonical,:tocanonicalD,:fromcanonical,:fromcanonicalD,:invfromcanonicalD)
    @eval $op(f::Fun,x...) = $op(domain(f),x...)
end

for op in (:tocanonical,:tocanonicalD)
    @eval $op(d::Domain) = $op(d,Fun(identity,d))
end
for op in (:fromcanonical,:fromcanonicalD,:invfromcanonicalD)
    @eval $op(d::Domain) = $op(d,Fun(identity,canonicaldomain(d)))
end


space(f::Fun) = f.space
spacescompatible(f::Fun,g::Fun) = spacescompatible(space(f),space(g))
canonicalspace(f::Fun) = canonicalspace(space(f))
canonicaldomain(f::Fun) = canonicaldomain(domain(f))


##Evaluation

function evaluate(f::AbstractVector,S::Space,x...)
    csp=canonicalspace(S)
    if spacescompatible(csp,S)
        error("Override evaluate for " * string(typeof(csp)))
    else
        evaluate(coefficients(f,S,csp),csp,x...)
    end
end

evaluate(f::Fun,x) = evaluate(f.coefficients,f.space,x)
evaluate(f::Fun,x,y,z...) = evaluate(f.coefficients,f.space,Vec(x,y,z...))


@compat (f::Fun)(x...) = evaluate(f,x...)

for op in (:(Base.first),:(Base.last))
    @eval $op{S,T}(f::Fun{S,T}) = f($op(domain(f)))
end



## Extrapolation


# Default extrapolation is evaluation. Override this function for extrapolation enabled spaces.
extrapolate(f::AbstractVector,S::Space,x...)=evaluate(f,S,x...)

# Do not override these
extrapolate(f::Fun,x) = extrapolate(f.coefficients,f.space,x)
extrapolate(f::Fun,x,y,z...) = extrapolate(f.coefficients,f.space,Vec(x,y,z...))


##Data routines


values(f::Fun,dat...) = itransform(f.space,f.coefficients,dat...)
points(f::Fun)=points(f.space,ncoefficients(f))
ncoefficients(f::Fun)=length(f.coefficients)


function Base.stride(f::Fun)
    # Check only for stride 2 at the moment
    # as higher stride is very rare anyways
    M=maxabs(f.coefficients)
    for k=2:2:ncoefficients(f)
        if abs(f.coefficients[k])>40*M*eps()
            return 1
        end
    end

    2
end



## Manipulate length

pad!(f::Fun,n::Integer) = (pad!(f.coefficients,n);f)
pad(f::Fun,n::Integer) = Fun(pad(f.coefficients,n),f.space)


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
                n = max(ncoefficients(f),ncoefficients(g))
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


# equivalent to Y+=a*X
axpy!(a,X::Fun,Y::Fun)=axpy!(a,coefficients(X,space(Y)),Y)
function axpy!(a,xcfs::Vector,Y::Fun)
    if a!=0
        n=ncoefficients(Y); m=length(xcfs)

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



-(f::Fun)=Fun(-f.coefficients,f.space)
for op in (:-,:(.-))
    @eval $op(c::Number,f::Fun)=-$op(f,c)
end


for op = (:*,:.*,:./,:/)
    @eval $op(f::Fun,c::Number) = Fun($op(f.coefficients,c),f.space)
end


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

export bilinearform, linebilinearform, innerproduct, lineinnerproduct

# Having fallbacks allow for the fast implementations.

defaultbilinearform(f::Fun,g::Fun)=sum(f.*g)
defaultlinebilinearform(f::Fun,g::Fun)=linesum(f.*g)

bilinearform(f::Fun,g::Fun)=defaultbilinearform(f,g)
bilinearform(c::Number,g::Fun)=sum(c*g)
bilinearform(g::Fun,c::Number)=sum(g*c)

linebilinearform(f::Fun,g::Fun)=defaultbilinearform(f,g)
linebilinearform(c::Number,g::Fun)=linesum(c*g)
linebilinearform(g::Fun,c::Number)=linesum(g*c)



# Conjugations

innerproduct(f::Fun,g::Fun)=bilinearform(conj(f),g)
innerproduct(c::Number,g::Fun)=bilinearform(conj(c),g)
innerproduct(g::Fun,c::Number)=bilinearform(conj(g),c)

lineinnerproduct(f::Fun,g::Fun)=linebilinearform(conj(f),g)
lineinnerproduct(c::Number,g::Fun)=linebilinearform(conj(c),g)
lineinnerproduct(g::Fun,c::Number)=linebilinearform(conj(g),c)

## Norm

for (OP,SUM) in ((:(Base.norm),:(Base.sum)),(:linenorm,:linesum))
    @eval begin
        $OP(f::Fun) = $OP(f,2)

        function $OP(f::Fun,p::Number)
            if p < 1
                return error("p should be 1 ≤ p ≤ ∞")
            elseif 1 ≤ p < Inf
                return abs($SUM(abs2(f)^(p/2)))^(1/p)
            else
                return maxabs(f)
            end
        end

        function $OP(f::Fun,p::Int)
            if 1 ≤ p < Inf
                return iseven(p) ? abs($SUM(abs2(f)^div(p,2)))^(1/p) : abs($SUM(abs2(f)^(p/2)))^(1/p)
            else
                return error("p should be 1 ≤ p ≤ ∞")
            end
        end
    end
end


## Mapped functions

Base.transpose(f::Fun) = f  # default no-op

for op = (:(Base.real),:(Base.imag),:(Base.conj))
    @eval ($op){S<:Space{RealBasis}}(f::Fun{S}) = Fun(($op)(f.coefficients),f.space)
end

Base.conj(f::Fun)=error("Override conj for $(typeof(f))")

Base.abs2{S<:Space{RealBasis},T<:Real}(f::Fun{S,T})=f^2
Base.abs2{S<:Space{RealBasis},T<:Complex}(f::Fun{S,T})=real(f)^2+imag(f)^2
Base.abs2(f::Fun)=f*conj(f)

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
        m=min(ncoefficients(f),ncoefficients(g))
        tol=100eps()  # TODO: normalize by norm of f/g

        for k=1:m
            if !isapprox(f.coefficients[k],g.coefficients[k])
                return false
            end
        end
        for k=m+1:ncoefficients(f)
            if abs(f.coefficients[k])>tol
                return false
            end
        end
        for k=m+1:ncoefficients(g)
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


## broadcasting

Base.broadcast(op,f::Fun) = Fun(x -> op(f(x)), domain(f))
Base.broadcast(op,f::Fun,c::Number) = Fun(x -> op(f(x),c), domain(f))
Base.broadcast(op,c::Number,f::Fun) = Fun(x -> op(c,f(x)), domain(f))
Base.broadcast(op,f::Fun,g::Fun) = Fun(x -> op(f(x),g(x)), domain(f) ∪ domain(g))


include("constructors.jl")
