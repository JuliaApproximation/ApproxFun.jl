
valsdomain_type_promote(::Type{T},::Type{T}) where {T<:Complex}=T,T
valsdomain_type_promote(::Type{T},::Type{T}) where {T<:Real}=T,T
valsdomain_type_promote(::Type{Int},::Type{Int})=Float64,Int
valsdomain_type_promote(::Type{T},::Type{Complex{V}}) where {T<:Real,V<:Real}=promote_type(T,V),Complex{promote_type(T,V)}
valsdomain_type_promote(::Type{Complex{T}},::Type{V}) where {T<:Real,V<:Real}=Complex{promote_type(T,V)},promote_type(T,V)
valsdomain_type_promote(::Type{T},::Type{Int}) where {T<:Integer}=Float64,Int
valsdomain_type_promote(::Type{T},::Type{Int}) where {T<:Real}=T,Int
valsdomain_type_promote(::Type{T},::Type{Int}) where {T<:Complex}=T,Int
valsdomain_type_promote(::Type{T},::Type{V}) where {T<:Integer,V<:Real}=valsdomain_type_promote(Float64,V)
valsdomain_type_promote(::Type{T},::Type{V}) where {T<:Integer,V<:Complex}=valsdomain_type_promote(Float64,V)
valsdomain_type_promote(::Type{T},::Type{Vector{T}}) where {T<:Real}=T,Vector{T}
valsdomain_type_promote(::Type{T},::Type{V}) where {T,V}=promote_type(T,V),promote_type(T,V)



function choosefuneltype(ftype,Td)
    if !( ftype<: Number || ( ((ftype <: AbstractArray) || (ftype <: Vec)) &&
                              (eltype(ftype) <: Number) ) )
        warn("Function outputs type $(ftype), which is not a Number")
    end

    Tprom = ftype

    if ftype <: Number #TODO should also work for array-valued functions
        Tprom,Tpromd=valsdomain_type_promote(ftype,Td)

        if ftype != Int && Tprom != ftype
                warn("Promoting function output type from $(ftype) to $(Tprom)")
        end
        if Tpromd != Td
                warn("Space domain number type $(Td) is not compatible with coefficient type $(Tprom)")
                #TODO should construct a new Space that contains a domain where the numbers have been promoted
                #and call constructor with this Space.
        end
    end

    Tprom
end


# default_Fun is the default constructor, based on evaluation and transforms
# last argument is whether to splat or not
default_Fun(::Type{T},f,d::Space{ReComp},pts::AbstractVector,::Type{Val{true}}) where {T,ReComp} =
    Fun(d,transform(d,T[f(x...) for x in pts]))

default_Fun(::Type{T},f,d::Space{ReComp},pts::AbstractVector,::Type{Val{false}}) where {T,ReComp} =
    Fun(d,transform(d,broadcast!(f, similar(pts, T), pts)))


function default_Fun(f,d::Space{ReComp},n::Integer,::Type{Val{false}}) where ReComp
    pts=points(d, n)
    f1=f(pts[1])
    if isa(f1,AbstractArray) && size(d) ≠ size(f1)
        return Fun(f,Space(fill(d,size(f1))),n)
    end

    # we need 3 eltype calls for the case Interval(Point([1.,1.]))
    Tprom=choosefuneltype(typeof(f1),prectype(domain(d)))
    default_Fun(Tprom,f,d,pts,Val{false})
end

function default_Fun(f,d::Space{ReComp},n::Integer,::Type{Val{true}}) where ReComp
    pts=points(d, n)
    f1=f(pts[1]...)
    if isa(f1,AbstractArray) && size(d) ≠ size(f1)
        return Fun(f,Space(fill(d,size(f1))),n)
    end

    # we need 3 eltype calls for the case Interval(Point([1.,1.]))
    Tprom=choosefuneltype(typeof(f1),prectype(domain(d)))
    default_Fun(Tprom,f,d,pts,Val{true})
end

default_Fun(f,d::Space{ReComp},n::Integer) where {ReComp} = default_Fun(f,d,n,Val{!hasnumargs(f,1)})

Fun(f::Function,d::Space{ReComp},n::Integer) where {ReComp} = default_Fun(dynamic(f),d,n)

# the following is to avoid ambiguity
# Fun(f::Fun,d) should be equivalent to Fun(x->f(x),d)
Fun(f::Fun,d::Space) = Fun(d,coefficients(f,d))
Fun(f::Fun,::Type{T}) where {T<:Space} = Fun(f,T(domain(f)))


Fun(f,T::Type) = Fun(dynamic(f),T())
Fun(f::Function,T::Type,n::Integer) = Fun(dynamic(f),T(),n)

Fun(f::AbstractVector,d::Domain) = Fun(f,Space(d))
Fun(d::Domain,f::AbstractVector{T}) where {T<:Number} = Fun(Space(d),f)
Fun(d::Domain,f::AbstractVector) = Fun(Space(d),f)


Fun(f::Function,d::Domain,n) = Fun(dynamic(f),Space(d),n)


# We do zero special since zero exists even when one doesn't
Fun(c::Number,::Type{T}) where {T<:Space} = c==0?zeros(T(AnyDomain())):c*ones(T(AnyDomain()))
Fun(c::Number,d::Domain) = c==0?c*zeros(d):c*ones(d)
Fun(c::Number,d::Space) = c==0?c*zeros(prectype(d),d):c*ones(prectype(d),d)

## Adaptive constructors
function default_Fun(f, d::Space)
    hasnumargs(f,1) || return Fun(xy->f(xy...),d)
    isinf(dimension(d)) || return Fun(f,d,dimension(d))  # use exactly dimension number of sample points

    #TODO: reuse function values?
    T = real(eltype(domain(d)))

    r=checkpoints(d)
    f0=f(first(r))

    isa(f0,AbstractArray) && size(d) ≠ size(f0) && return Fun(f,Space(fill(d,size(f0))))



    tol =T==Any?20eps():20eps(T)


    fr=map(f,r)
    maxabsfr=norm(fr,Inf)

    for logn = 4:20
        #cf = Fun(f, d, 2^logn + 1)
        cf = default_Fun(f, d, 2^logn)
        maxabsc = maximum(abs,cf.coefficients)
        if maxabsc == 0 && maxabsfr == 0
            return(zeros(d))
        end

        b = block(d,length(cf.coefficients))
        bs = blockstart(d,max(b-2,Block(1)))

        # we allow for transformed coefficients being a different size
        ##TODO: how to do scaling for unnormalized bases like Jacobi?
        if ncoefficients(cf) > 8 && maximum(abs,cf.coefficients[bs:end]) < 10tol*maxabsc &&
                all(k->norm(cf(r[k])-fr[k],1)<tol*length(cf.coefficients)*maxabsfr*1000,1:length(r))
            return chop!(cf,tol)
        end
    end

    warn("Maximum number of coefficients "*string(2^20+1)*" reached in constructing Fun.")

    Fun(f,d,2^21)
end

Fun(f::Type, d::Space) = error("Not implemented")


# special case constructors
## TODO: remove zeros
Base.zero(S::Space) = zeros(S)
Base.zero(::Type{T},S::Space) where {T<:Number} = zeros(T,S)
Base.zeros(::Type{T},S::Space) where {T<:Number} = Fun(S,zeros(T,1))
Base.zeros(S::Space) = Fun(S,zeros(1))

# catch all
Base.ones(S::Space) = Fun(x->1.0,S)
Base.ones(::Type{T},S::Space) where {T<:Number} = Fun(x->one(T),S)

function Fun(::typeof(identity), d::Domain)
    cd=canonicaldomain(d)
    if typeof(d) == typeof(cd)
        Fun(dynamic(x->x),d) # fall back to constructor, can't use `identity` as that creates a loop
    else
        # this allows support for singularities, that the constructor doesn't
        sf=fromcanonical(d,Fun(identity,cd))
        Fun(setdomain(space(sf),d),coefficients(sf))
    end
end

Fun(::typeof(identity), S::Space) = Fun(identity,domain(S))



Fun(f::typeof(zero), d::Space) = zeros(eltype(domain(d)),d)
Fun(f::typeof(one), d::Space) = ones(eltype(domain(d)),d)

Fun(f::Type, d::Domain) = Fun(f,Space(d))
Fun(f::Function, d::Domain) = Fun(f,Space(d))


# this is the main constructor
Fun(f::Function, d::Space) = default_Fun(dynamic(f), d)

# this supports expanding a Fun to a larger or smaller domain.
# we take the union and then intersection to get at any singularities
# TODO: singularities in space(f)
Fun(f::Fun, d::Domain) = Fun(f,Space((d ∪ domain(f)) ∩ d))





## Aliases



Fun(T::Type,n::Integer) = Fun(T(),n)
Fun(f,n::Integer) = Fun(f,Interval(),n)
Fun(f,d::ClosedInterval,n::Integer) = Fun(f,Domain(d),n)
Fun(d::ClosedInterval,cfs::AbstractVector{M}) where {M<:Number} = Fun(Domain(d),1.0*cfs)
Fun(f::Function,d::ClosedInterval) = Fun(dynamic(f),Domain(d))
Fun(f,d::ClosedInterval) = Fun(f,Domain(d))
Fun(f::Number,d::ClosedInterval) = Fun(f,Domain(d))
Fun(d::ClosedInterval) = Fun(Domain(d))

Fun(T::Type,d::AbstractVector) = Fun(T(),d)

Fun(f::Fun{SequenceSpace},s::Space) = Fun(s,f.coefficients)
