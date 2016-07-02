export domainplot, coefficientplot, complexplot



if isdir(Pkg.dir("TikzGraphs"))
    include("introspect.jl")
end




## Fun routines


function plotptsvals(f::Fun)
    if isinf(dimension(space(f)))
        f=pad(f,3ncoefficients(f)+50)
    else
        f=pad(f,dimension(space(f)))
    end
    return points(f),values(f)
end

function plotptsvals{S<:JacobiWeight}(f::Fun{S})
    f=pad(f,3ncoefficients(f)+50)
    s=space(f)
    pts,vals=points(f),values(f)
    # add endpoints so that singularity is viewable
    if s.α ≥ 0
        pts=insert!(pts,1,first(domain(f)))
        vals=insert!(vals,1,first(f))
    end
    if s.β ≥ 0
        pts=push!(pts,last(domain(f)))
        vals=push!(vals,last(f))
    end

    pts,vals
end


## Recipes


@recipe function f{S,T<:Real}(g::Fun{S,T})
    plotptsvals(g)
end
@recipe function f{S,T<:Real}(g::Fun{S,Complex{T}})
    x,v=plotptsvals(g)
    x,Vector{T}[real(v),imag(v)]
end


@recipe function f{S,T<:Real}(x::AbstractVector{T},g::Fun{S,T})
    x,g(x)
end

@recipe function f{S,T<:Real}(x::AbstractVector{T},g::Fun{S,Complex{T}})
    v=g(x)
    x,Vector{T}[real(v),imag(v)]
end


@recipe function f{F<:Fun}(G::AbstractVector{F})
    x=Vector{Float64}[]
    v=Vector{Float64}[]
    for g in G
        xx,vv=plotptsvals(g)
        push!(x,xx)
        push!(v,vv)
    end
    x,v
end

@recipe function f{T<:Real,F<:Fun}(x::AbstractVector{T},G::AbstractVector{F})
    v=Vector{Float64}[]
    for g in G
        push!(v,g(x))
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

@userplot ComplexPlot

@recipe function f(h::ComplexPlot)
    @assert length(h.args)==1
    complexplotvals(h.args[1])
end

@recipe function f(dd::Domain)
    complexplotvals(Fun(identity,dd))
end


@recipe function f{F<:Domain}(G::AbstractVector{F})
    x=Vector{Float64}[]
    v=Vector{Float64}[]
    for g in G
        xx,vv=complexplotvals(g)
        push!(x,xx)
        push!(v,vv)
    end
    x,v
end



@userplot DomainPlot

@recipe function f(D::DomainPlot)
    @assert length(D.args)==1
    g=D.args[1]
    @assert isa(g,Fun)

    domain(g)
end



@userplot CoefficientPlot

@recipe function f(C::CoefficientPlot)
    @assert length(C.args)==1
    g=C.args[1]
    @assert isa(g,Fun)

    yscale --> :log10

    abs(g.coefficients)
end



@recipe function f{S<:Union{ArraySpace,TupleSpace},T<:Real}(g::Fun{S,T})
    vec(g)
end

@recipe function f{S<:PiecewiseSpace,T<:Real}(g::Fun{S,T})
    p=pieces(g)
    for k=1:length(p)
        @series begin
            primary := (k==1)
            p[k]
        end
    end
end



@recipe function f{S<:DiracSpace,T<:Real}(g::Fun{S,T})
    pts=space(g).points
    n=length(pts)
    ws=pad(g.coefficients,length(pts))

    xlims --> (minimum(pts)-1.,maximum(pts)+1.)

    for k=1:length(pts)
        @series begin
            primary := (k==1)
            ones(2)*pts[k],[0,1]*ws[k]
        end
        @series begin
            linestyle := :dot
            primary := false
            ones(2)*pts[k],[1,2]*ws[k]
        end
    end
end

@recipe function f{S<:PointSpace,T<:Real}(g::Fun{S,T})
    pts=space(g).points
    n=length(pts)
    ws=pad(g.coefficients,length(pts))

    xlims --> (minimum(pts)-1.,maximum(pts)+1.)

    for k=1:length(pts)
        @series begin
            primary := (k==1)
            ones(2)*pts[k],[0,1]*ws[k]
        end
    end
end



@recipe function f{S<:HeavisideSpace,T<:Real}(g::Fun{S,T})
    pts=domain(g).points
    n=length(pts)
    ws=pad(g.coefficients,dimension(space(g)))

    ylims --> (minimum(ws)-1.,maximum(ws)+1.)
    @series begin
        primary := true
        pts[1:2],ones(2)*ws[1]
    end

    for k=2:length(ws)
        @series begin
           primary := false
           linestyle --> :dot
           ones(2)*pts[k],ws[k-1:k]
        end
        @series begin
            primary := false
            pts[k:k+1],ones(2)*ws[k]
        end
    end
end

###
# Multivariate
###


@recipe function f{S<:UnivariateSpace,
                    V<:UnivariateSpace,
        SV<:TensorSpace}(g::ProductFun{S,V,SV})
    g=chop(g,10e-10)
    g=pad(g,max(size(g,1),20),max(size(g,2),20))
    vals=values(g)

    seriestype --> :surface

    if norm(imag(vals),Inf)>10e-9
        warn("Imaginary part is non-neglible.  Only plotting real part.")
    end

    points(space(g,1),size(vals,1)),points(space(g,2),size(vals,2)),real(vals).'
end


@recipe function f(x::AbstractVector,y::AbstractVector,g::MultivariateFun)
    seriestype --> :surface
    x,y,real(g(x,y)).'
end


@recipe function f{TS<:AbstractProductSpace,T<:Real}(g::Fun{TS,T})
    ProductFun(g)
end

@recipe function f{TS<:AbstractProductSpace,T<:Complex}(g::Fun{TS,T})
    ProductFun(g)
end

@recipe function f{TS<:AbstractProductSpace,T<:Real}(x::AbstractVector,y::AbstractVector,g::Fun{TS,T})
    x,y,ProductFun(g)
end

@recipe function f{TS<:AbstractProductSpace,T<:Complex}(x::AbstractVector,y::AbstractVector,g::Fun{TS,T})
    x,y,ProductFun(g)
end
