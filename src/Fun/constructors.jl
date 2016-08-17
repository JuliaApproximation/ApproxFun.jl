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
    if !( ftype<: Number || ( (ftype <: Array) && (ftype.parameters[1] <: Number) ) )
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


function defaultFun{ReComp}(f,d::Space{ReComp},n::Integer)
    pts=points(d, n)
    if !hasnumargs(f,1)  # Splat out Vec
        return Fun(xy->f(xy...),d,n)
    end

    f1=f(pts[1])

    if isa(f1,Array) && !isa(d,ArraySpace)
        return Fun(f,ArraySpace(d,size(f1)...),n)
    end


    # we need 3 eltype calls for the case Interval(Point([1.,1.]))
    Tprom=choosefuneltype(typeof(f1),eltype(eltype(eltype(domain(d)))))


    vals=Tprom[f(x) for x in pts]
    Fun(transform(d,vals),d)
end

Fun{ReComp}(f,d::Space{ReComp},n::Integer)=defaultFun(f,d,n)

# the following is to avoid ambiguity
# Fun(f::Fun,d) should be equivalent to Fun(x->f(x),d)
#TODO: fall back to Fun(x->f(x),d) if conversion not implemented?
Fun(f::Fun,d::Space)=Fun(coefficients(f,d),d)
Fun{T<:Space}(f::Fun,::Type{T})=Fun(f,T(domain(f)))



Fun(f::AbstractVector,T::Type)=Fun(f,T())

Fun(f,T::Type)=Fun(f,T())
Fun(f,T::Type,n::Integer)=Fun(f,T(),n)

Fun(f::AbstractVector,d::Domain)=Fun(f,Space(d))


Fun(f,d::Domain,n)=Fun(f,Space(d),n)


# We do zero special since zero exists even when one doesn't
Fun{T<:Space}(c::Number,::Type{T})=c==0?zeros(T(AnyDomain())):c*ones(T(AnyDomain()))
Fun(c::Number,d::Domain)=c==0?c*zeros(d):c*ones(d)
Fun(c::Number,d::Space)=c==0?c*zeros(eltype(d),d):c*ones(eltype(d),d)


## List constructor

Fun{T<:Domain}(c::Number,dl::AbstractVector{T})=Fun(c,UnionDomain(dl))
Fun{T<:Domain}(f,dl::AbstractVector{T})=Fun(f,UnionDomain(dl))
Fun{T<:Domain}(f,dl::AbstractVector{T},n::Integer)=Fun(f,UnionDomain(dl),n)

## Adaptive constructors

function randomFun(f,d::IntervalDomain)
    @assert d == Interval()

    #TODO: implement other domains

    Fun(chebyshevtransform(randomadaptivebary(f)),d)
end


# function veczerocfsFun(f,d::IntervalDomain)
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
#     warn("Maximum number of coefficients reached")
#
#     Fun(f,d,2^21 + 1)
# end



function zerocfsFun(f, d::Space)
    #TODO: reuse function values?
    T = eltype(domain(d))
    if T <: Complex
        T = T.parameters[1] #get underlying real representation
    end
    r=checkpoints(d)
    f0=f(first(r))

    if !isa(d,ArraySpace) && isa(f0,Array)
        return zerocfsFun(f,ArraySpace(d,size(f0)...))
    end

    tol =T==Any?200eps():200eps(T)


    fr=map(f,r)
    maxabsfr=norm(fr,Inf)

    for logn = 4:20
        #cf = Fun(f, d, 2^logn + 1)
        cf = defaultFun(f, d, 2^logn)
        maxabsc=maxabs(cf.coefficients)
        if maxabsc==0 && maxabsfr==0
            return(zeros(d))
        end

        # we allow for transformed coefficients being a different size
        ##TODO: how to do scaling for unnormalized bases like Jacobi?
        if ncoefficients(cf) > 8 && maxabs(cf.coefficients[end-8:end]) < tol*maxabsc &&
                all(k->norm(cf(r[k])-fr[k],1)<1E-4,1:length(r))
            return chop!(cf,tol*maxabsc/10)
        end
    end

    warn("Maximum number of coefficients "*string(2^20+1)*" reached")

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

        if maxabs(cf.coefficients[end-8:end]) < tol
            return chop!(cf,10eps(T))
        end
    end

    warn("Maximum number of coefficients "*string(2^20+1)*" reached")

    Fun(f,d,2^21)
end


function Fun(f, d::Space; method="zerocoefficients")
    T = eltype(domain(d))

    if f==identity
        identity_fun(d)
    elseif f==zero # zero is always defined
        zeros(T,d)
    elseif f==one
        ones(T,d)
    elseif !hasnumargs(f,1)  # Splat out Vec
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
Fun(f,d::Domain;opts...)=Fun(f,Space(d);opts...)


# this supports expanding a Fun to a larger or smaller domain.
# we take the union and then intersection to get at any singularities
# TODO: singularities in space(f)
Fun(f::Fun,d::Domain;opts...)=Fun(f,Space((d ∪ domain(f)) ∩ d);opts...)





## Aliases




Fun(f,n::Integer)=Fun(f,Interval(),n)
Fun{T<:Number}(f,d::AbstractVector{T},n::Integer)=Fun(f,Domain(d),n)
Fun{T<:Number,M<:Number}(cfs::AbstractVector{M},d::AbstractVector{T})=Fun(1.0*cfs,Domain(d))
Fun{T<:Number}(f,d::AbstractVector{T})=Fun(f,Domain(d))
Fun{T<:Number}(f::Number,d::AbstractVector{T})=Fun(f,Domain(d))
Fun(f::AbstractVector)=Fun(Domain(f))

function Fun(cfs::AbstractVector{Any},s::Space)
    @assert isempty(cfs)
    Fun(Float64[],s)
end

Fun(f::Fun{SequenceSpace},s::Space) = Fun(f.coefficients,s)
