export domainplot, coefficientplot, complexplot



## Fun routines


function plotptsvals(f::Fun)
    if isinf(dimension(space(f)))
        f=pad(f,3ncoefficients(f)+50)
    else
        f=pad(f,dimension(space(f)))
    end
    return points(f),values(f)
end

function plotptsvals(f::Fun{S}) where S<:JacobiWeight
    f=pad(f,3ncoefficients(f)+50)
    s=space(f)
    pts,vals=points(f),values(f)
    # add endpoints so that singularity is viewable
    if s.β ≥ 0
        pts=push!(pts,first(domain(f)))
        vals=push!(vals,first(f))
    end
    if s.α ≥ 0
        pts=insert!(pts,1,last(domain(f)))
        vals=insert!(vals,1,last(f))
    end

    pts,vals
end


## Recipes


@recipe function f(g::Fun{S,T}) where {S,T<:Real}
    plotptsvals(g)
end
@recipe function f(g::Fun{S,Complex{T}}) where {S,T<:Real}
    x,v=plotptsvals(g)
    x,Vector{T}[real(v),imag(v)]
end


@recipe function f(x::Fun{V,T},y::Fun{S,T}) where {S,V,T<:Real}
    M=3max(ncoefficients(x),ncoefficients(y))+50
    values(pad(x,M)),values(pad(y,M))
end


@recipe function f(x::AbstractVector{T},g::Fun{S,T}) where {S,T<:Real}
    x,g.(x)
end

@recipe function f(x::AbstractVector{T},g::Fun{S,Complex{T}}) where {S,T<:Real}
    v=g.(x)
    x,Vector{T}[real(v),imag(v)]
end


@recipe function f(G::AbstractVector{F}) where F<:Fun
    x=Vector{Float64}[]
    v=Vector{Float64}[]
    for g in G
        xx,vv=plotptsvals(g)
        push!(x,xx)
        push!(v,vv)
    end
    x,v
end

@recipe function f(x::AbstractVector{T},G::AbstractVector{F}) where {T<:Real,F<:Fun}
    v=Vector{Float64}[]
    for g in G
        push!(v,g.(x))
    end
    x,v
end


function complexplotvals(f::Fun)
    vals =plotptsvals(f)[2]
    if isa(domain(f),PeriodicDomain)
        real([vals;vals[1]]),imag([vals;vals[1]])
    else
        real(vals),imag(vals)
    end
end

@recipe function f(::Type{Val{:complexplot}},h)
    complexplotvals(h)
end

@recipe function f(dd::Domain)
    x,y = complexplotvals(Fun(identity,dd))
    # @series (x,y)
    # @series begin
    #     primary := false
    #     ∂(dd)
    # end
end

@recipe function f(dd::UnionDomain)
    @series component(dd,1)
    for k=2:ncomponents(dd)
        @series begin
            primary := false
            component(dd,k)
        end
    end
end


@recipe function f(G::AbstractVector{F}) where F<:Domain
    x=Vector{Float64}[]
    v=Vector{Float64}[]
    for g in G
        xx,vv=complexplotvals(g)
        push!(x,xx)
        push!(v,vv)
    end
    x,v
end



# @userplot DomainPlot
#
# @recipe function f(D::DomainPlot)
#     @assert length(D.args)==1
#     g=D.args[1]
#     @assert isa(g,Fun)
#
#     domain(g)
# end
#
#
#
# @userplot CoefficientPlot
#
# @recipe function f(C::CoefficientPlot)
#     @assert length(C.args)==1
#     g=C.args[1]
#     @assert isa(g,Fun)
#
#     seriestype --> :scatter
#     yscale --> :log10
#     markersize := max(round(Int,5 - log10(ncoefficients(g))),1)
#     xlims --> (0,1.01ncoefficients(g))
#     abs(g.coefficients)
# end



@recipe function f(g::Fun{S,T}) where {S<:ArraySpace,T<:Real}
    components(g)
end

@recipe function f(g::Fun{S,T}) where {S<:PiecewiseSpace,T<:Real}
    p=components(g)
    for k=1:length(p)
        @series begin
            primary := (k==1)
            p[k]
        end
    end
end



@recipe function f(g::Fun{S,T}) where {S<:DiracSpace,T<:Real}
    pts=space(g).points
    n=length(pts)
    ws=pad(g.coefficients,length(pts))

    xlims --> (minimum(pts)-1.,maximum(pts)+1.)

    @series begin
        primary --> true
        ones(2)*pts[1],[0,1]*ws[1]
    end

    if length(ws) > 1
        @series begin
            primary := false
            ones(2)*pts[2:end]',[0,1]*ws[2:end]'
        end
    end

    @series begin
        primary := false
        linestyle := :dot
        ones(2)*pts',[1,2]*ws'
    end
end

@recipe function f(g::Fun{S,T}) where {S<:PointSpace,T<:Real}
    pts=space(g).points
    n=length(pts)
    ws=pad(g.coefficients,length(pts))

    xlims --> (minimum(pts)-1.,maximum(pts)+1.)

    @series begin
        primary --> true
        ones(2)*pts[1],[0,1]*ws[1]
    end
    if length(ws) > 1
        @series begin
            primary := false
            ones(2)*pts[2:end]',[0,1]*ws[2:end]'
        end
    end
end



@recipe function f(g::Fun{S,T}) where {S<:HeavisideSpace,T<:Real}
    pts=domain(g).points
    n=length(pts)
    ws=pad(g.coefficients,dimension(space(g)))

    lnsx=Vector{Float64}(0)
    lnsy=Vector{Float64}(0)
    dtsx=Vector{Float64}(0)
    dtsy=Vector{Float64}(0)
    for k=1:n-1
        push!(lnsx,pts[k])
        push!(lnsy,ws[k])
        push!(lnsx,pts[k+1])
        push!(lnsy,ws[k])


        if k != n-1
            # dotted line, with NaN to separate
            push!(dtsx,pts[k+1])
            push!(dtsx,pts[k+1])
            push!(dtsx,pts[k+1])

            push!(dtsy,ws[k])
            push!(dtsy,ws[k+1])
            push!(dtsy,NaN)

            # extra point for NaN
            push!(lnsx,pts[k+1])
            push!(lnsy,NaN)
        end
    end


    @series begin
        primary --> true
        lnsx,lnsy
    end

    @series begin
       primary := false
       linestyle := :dot
       fillrange := nothing
       dtsx,dtsy
    end
end

###
# Multivariate
###


@recipe function f(g::ProductFun{S,V,SV}) where {S<:UnivariateSpace,
                    V<:UnivariateSpace,
        SV<:TensorSpace}
    g=chop(g,10e-10)
    g=pad(g,max(size(g,1),20),max(size(g,2),20))
    vals=values(g)

    seriestype --> :surface

    if norm(imag(vals),Inf)>10e-9
        warn("Imaginary part is non-neglible.  Only plotting real part.")
    end

    # sort the points
    x = points(factor(space(g),1),size(vals,1))
    y = points(factor(space(g),2),size(vals,2))
    px = sortperm(x)
    py = sortperm(y)

    x[px],y[py],transpose(real(vals))[py,px]
end

@recipe function f(g::LowRankFun{S,V,SV}) where {S<:UnivariateSpace,
                    V<:UnivariateSpace,
        SV<:TensorSpace}
    g=chop(g,10e-10)
    vals=values(g)

    seriestype --> :surface

    if norm(imag(vals),Inf)>10e-9
        warn("Imaginary part is non-neglible.  Only plotting real part.")
    end


    # sort the points
    x = points(factor(space(g),1),size(vals,1))
    y = points(factor(space(g),2),size(vals,2))
    px = sortperm(x)
    py = sortperm(y)

    x[px],y[py],transpose(real(vals))[py,px]
end


@recipe function f(x::AbstractVector,y::AbstractVector,g::MultivariateFun)
    seriestype --> :surface
    x,y,real(g.(x,transpose(y)))
end


@recipe function f(g::Fun{TS,T}) where {TS<:BivariateSpace,T<:Real}
    g = pad(g,ncoefficients(g)+100)
    pts = points(g)
    seriestype --> :surface
    first.(pts),last.(pts),values(g)
end

@recipe function f(g::Fun{TS,T}) where {TS<:AbstractProductSpace,T<:Real}
    ProductFun(g)
end

@recipe function f(g::Fun{TS,T}) where {TS<:AbstractProductSpace,T<:Complex}
    ProductFun(g)
end

@recipe function f(x::AbstractVector,y::AbstractVector,g::Fun{TS,T}) where {TS<:AbstractProductSpace,T<:Real}
    x,y,ProductFun(g)
end

@recipe function f(x::AbstractVector,y::AbstractVector,g::Fun{TS,T}) where {TS<:AbstractProductSpace,T<:Complex}
    x,y,ProductFun(g)
end
