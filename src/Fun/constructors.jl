struct F <: Function
    f
end
(f::F)(args...) = f.f(args...)

valsdomain_type_promote{T<:Complex}(::Type{T},::Type{T})=T,T
valsdomain_type_promote{T<:Real}(::Type{T},::Type{T})=T,T
valsdomain_type_promote(::Type{Int},::Type{Int})=Float64,Int
valsdomain_type_promote{T<:Real,V<:Real}(::Type{T},::Type{Complex{V}})=promote_type(T,V),Complex{promote_type(T,V)}
valsdomain_type_promote{T<:Real,V<:Real}(::Type{Complex{T}},::Type{V})=Complex{promote_type(T,V)},promote_type(T,V)
valsdomain_type_promote{T<:Integer}(::Type{T},::Type{Int})=Float64,Int
valsdomain_type_promote{T<:Real}(::Type{T},::Type{Int})=T,Int
valsdomain_type_promote{T<:Complex}(::Type{T},::Type{Int})=T,Int
valsdomain_type_promote{T<:Integer,V<:Real}(::Type{T},::Type{V})=valsdomain_type_promote(Float64,V)
valsdomain_type_promote{T<:Integer,V<:Complex}(::Type{T},::Type{V})=valsdomain_type_promote(Float64,V)
valsdomain_type_promote{T<:Real}(::Type{T},::Type{Vector{T}})=T,Vector{T}
valsdomain_type_promote{T,V}(::Type{T},::Type{V})=promote_type(T,V),promote_type(T,V)



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

# last argument is whether to splat or not
defaultFun{T,ReComp}(::Type{T},f,d::Space{ReComp},pts::AbstractVector,::Type{Val{true}}) =
    Fun(d,transform(d,T[f(x...) for x in pts]))

defaultFun{T,ReComp}(::Type{T},f,d::Space{ReComp},pts::AbstractVector,::Type{Val{false}}) =
    Fun(d,transform(d,T[f(x) for x in pts]))


function defaultFun{ReComp}(f,d::Space{ReComp},n::Integer,::Type{Val{false}})
    pts=points(d, n)
    f1=f(pts[1])
    if (isa(f1,AbstractArray) || isa(f1,Vec)) && !isa(d,ArraySpace)
        return Fun(f,ArraySpace(d,size(f1)...),n)
    end

    # we need 3 eltype calls for the case Interval(Point([1.,1.]))
    Tprom=choosefuneltype(typeof(f1),prectype(domain(d)))
    defaultFun(Tprom,f,d,pts,Val{false})
end

function defaultFun{ReComp}(f,d::Space{ReComp},n::Integer,::Type{Val{true}})
    pts=points(d, n)
    f1=f(pts[1]...)
    if (isa(f1,AbstractArray) || isa(f1,Vec)) && !isa(d,ArraySpace)
        return Fun(f,ArraySpace(d,size(f1)...),n)
    end

    # we need 3 eltype calls for the case Interval(Point([1.,1.]))
    Tprom=choosefuneltype(typeof(f1),prectype(domain(d)))
    defaultFun(Tprom,f,d,pts,Val{true})
end

defaultFun{ReComp}(f::F,d::Space{ReComp},n::Integer) = defaultFun(f,d,n,Val{!hasnumargs(f.f,1)})


Fun{ReComp}(f::Function,d::Space{ReComp},n::Integer) = Fun(F(f),d,n)
Fun{ReComp}(f::F,d::Space{ReComp},n::Integer) = defaultFun(f,d,n)

# the following is to avoid ambiguity
# Fun(f::Fun,d) should be equivalent to Fun(x->f(x),d)
#TODO: fall back to Fun(x->f(x),d) if conversion not implemented?
Fun(f::Fun,d::Space) = Fun(d,coefficients(f,d))
Fun{T<:Space}(f::Fun,::Type{T}) = Fun(f,T(domain(f)))



Fun(f::AbstractVector,T::Type) = Fun(f,T())
Fun(T::Type,f::Function) = Fun(T,F(f))
Fun(T::Type,f::Type) = Fun(T(),f)
Fun(T::Type,f)  =  Fun(T(),f)
Fun(f::Function,T::Type) = Fun(F(f),T())
Fun(f,T::Type) = Fun(f,T())
Fun(f::Function,T::Type,n::Integer) = Fun(F(f),T(),n)
Fun(f,T::Type,n::Integer) = Fun(f,T(),n)

Fun(f::AbstractVector,d::Domain) = Fun(f,Space(d))
Fun{T<:Number}(d::Domain,f::AbstractVector{T}) = Fun(Space(d),f)
Fun(d::Domain,f::AbstractVector) = Fun(Space(d),f)


Fun(f::Function,d::Domain,n) = Fun(F(f),Space(d),n)
Fun(f,d::Domain,n) = Fun(f,Space(d),n)


# We do zero special since zero exists even when one doesn't
Fun{T<:Space}(c::Number,::Type{T}) = c==0?zeros(T(AnyDomain())):c*ones(T(AnyDomain()))
Fun(c::Number,d::Domain) = c==0?c*zeros(d):c*ones(d)
Fun(c::Number,d::Space) = c==0?c*zeros(prectype(d),d):c*ones(prectype(d),d)


## List constructor

Fun{T<:Domain}(c::Number,dl::AbstractVector{T}) = Fun(c,UnionDomain(dl))
Fun{T<:Domain}(f::Function,dl::AbstractVector{T}) = Fun(F(f),UnionDomain(dl))
Fun{T<:Domain}(f::Type,dl::AbstractVector{T}) = Fun(f,UnionDomain(dl))
Fun{T<:Domain}(f::Domain,dl::AbstractVector{T}) = Fun(f,UnionDomain(dl))
Fun{T<:Domain}(f,dl::AbstractVector{T}) = Fun(f,UnionDomain(dl))
Fun{T<:Domain}(f,dl::AbstractVector{T},n::Integer) = Fun(f,UnionDomain(dl),n)

## Adaptive constructors

function randomFun(f,d::IntervalDomain)
    @assert d == Interval()

    #TODO: implement other domains

    Fun(d,chebyshevtransform(randomadaptivebary(f)))
end



function zerocfsFun(f, d::Space)
    #TODO: reuse function values?
    T = real(eltype(domain(d)))

    r=checkpoints(d)
    f0=f(first(r))

    if !isa(d,ArraySpace) && isa(f0,Array)
        return zerocfsFun(f,ArraySpace(d,size(f0)...))
    end

    tol =T==Any?20eps():20eps(T)


    fr=map(f,r)
    maxabsfr=norm(fr,Inf)

    for logn = 4:20
        #cf = Fun(f, d, 2^logn + 1)
        cf = defaultFun(f, d, 2^logn)
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


function abszerocfsFun(f,d::Space)
    #reuse function values
    T = eltype(domain(d))
    if T <: Complex
        T = T.parameters[1] #get underlying real representation
    end

    tol = 200eps(T)

    for logn = 4:20
        #cf = Fun(f, d, 2^logn + 1)
        cf = Fun(f, d, 2^logn)

        if maximum(abs,cf.coefficients[end-8:end]) < tol
            return chop!(cf,10eps(T))
        end
    end

    warn("Maximum number of coefficients "*string(2^20+1)*" reached")

    Fun(f,d,2^21)
end

Fun(f::Type, d::Space; method="zerocoefficients") = error("Not impleemnted")
Fun(f::Function, d::Space; method = "zerocoefficients") = Fun(F(f), d; method = method)
function Fun(f::F, d::Space; method="zerocoefficients")
    T = eltype(domain(d))

    if f.f==identity
        identity_fun(d)
    elseif f.f==zero # zero is always defined
        zeros(T,d)
    elseif f.f==one
        ones(T,d)
    elseif !hasnumargs(f.f,1)  # Splat out Vec
        Fun(xy->f(xy...),d;method=method)
    elseif !isinf(dimension(d))
        Fun(f,d,dimension(d))  # use exactly dimension number of sample points
    elseif method == "zerocoefficients"
        zerocfsFun(f,d)
    elseif method == "abszerocoefficients"
        abszerocfsFun(f,d)
    else
        randomFun(f,d)
    end
end
Fun(f::Type,d::Domain;opts...) = Fun(f,Space(d);opts...)
Fun(f::F,d::Domain;opts...) = Fun(f,Space(d);opts...)
Fun(f::Function,d::Domain;opts...) = Fun(F(f),d;opts...)

# this supports expanding a Fun to a larger or smaller domain.
# we take the union and then intersection to get at any singularities
# TODO: singularities in space(f)
Fun(f::Fun,d::Domain;opts...) = Fun(f,Space((d ∪ domain(f)) ∩ d);opts...)





## Aliases



Fun(T::Type,n::Integer) = Fun(T(),n)
Fun(f,n::Integer) = Fun(f,Interval(),n)
Fun(f,d::ClosedInterval,n::Integer) = Fun(f,Domain(d),n)
Fun{M<:Number}(d::ClosedInterval,cfs::AbstractVector{M}) = Fun(Domain(d),1.0*cfs)
Fun(f::Function,d::ClosedInterval) = Fun(F(f),Domain(d))
Fun(f::Type,d::ClosedInterval) = Fun(f,Domain(d))
Fun(f,d::ClosedInterval) = Fun(f,Domain(d))
Fun(f::Number,d::ClosedInterval) = Fun(f,Domain(d))
Fun(d::ClosedInterval) = Fun(Domain(d))

Fun{TT<:Number}(T::Type,d::AbstractVector{TT}) = Fun(T(),d)

Fun(f::Fun{SequenceSpace},s::Space) = Fun(s,f.coefficients)
