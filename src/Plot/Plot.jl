export domainplot, coefficientplot, complexplot



if isdir(Pkg.dir("TikzGraphs"))
    include("introspect.jl")
end



#surf(x...;opts...)=pysurf(x...;opts...)


## Fun routines

Plots.plot(f::Union{Fun,Domain,MultivariateFun};grid=true,kwds...)=plot!(plot(grid=grid),f;kwds...)
Plots.plot!(f::Union{Fun,Domain,MultivariateFun};kwds...)=plot!(current(),f;kwds...)

Plots.plot(x::AbstractVector,f::Fun;grid=true,kwds...)=plot!(plot(grid=grid),x,f;kwds...)
Plots.plot!(x::AbstractVector,f::Fun;kwds...)=plot!(current(),x,f;kwds...)



Plots.plot{F<:Union{Fun,Domain,MultivariateFun}}(v::AbstractVector{F};grid=true,kwds...)=plot!(plot(grid=grid),v;kwds...)
Plots.plot!{F<:Union{Fun,Domain,MultivariateFun}}(v::AbstractVector{F};kwds...)=plot!(current(),v;kwds...)


Plots.plot{F<:Fun}(x::AbstractVector,v::AbstractVector{F};grid=true,kwds...)=plot!(plot(grid=grid),x,v;kwds...)
Plots.plot!{F<:Fun}(x::AbstractVector,v::AbstractVector{F};kwds...)=plot!(current(),x,v;kwds...)


function plotptsvals(f::Fun)
    if dimension(space(f)) == Inf
        f=pad(f,3length(f)+50)
    else
        f=pad(f,dimension(space(f)))
    end
    return points(f),values(f)
end
Plots.plot!{S,T<:Real}(plt::Plots.Plot,f::Fun{S,T};kwds...)=
                plot!(plt,plotptsvals(f)...;kwds...)

Plots.plot!{S,T<:Real}(plt::Plots.Plot,x::AbstractVector{T},f::Fun{S,T};kwds...)=
                plot!(plt,x,f(x);kwds...)



function Plots.plot!{F<:Union{Fun,Domain}}(plt::Plots.Plot,v::AbstractVector{F};label=Void,kwds...)
    if label == Void
        for k=1:length(v)
            plot!(plt,v[k];kwds...)
        end
    else
        @assert length(label)==length(v)

        for k=1:length(v)
            plot!(plt,v[k];label=label[k],kwds...)
        end
    end
    plt
end

function Plots.plot!{F<:Fun}(plt::Plots.Plot,x::AbstractVector,v::AbstractVector{F};label=Void)
    if label == Void
        for k=1:length(v)
            plot!(plt,x,v[k])
        end
    else
        @assert length(label)==length(v)

        for k=1:length(v)
            plot!(plt,x,v[k];label=label[k])
        end
    end
    plt
end





Plots.plot!{S,T<:Complex}(plt::Plots.Plot,f::Fun{S,T};label=Void)=
                plot!(plt,[real(f),imag(f)];label=(label==Void?["Real","Imag"]:label))

Plots.plot!{S,T<:Complex}(plt::Plots.Plot,x::AbstractVector,f::Fun{S,T};label=Void)=
            plot!(plt,x,[real(f),imag(f)];label=(label==Void?["Real","Imag"]:label))

function complexplot!(plt::Plots.Plot,f::Fun;opts...)
    vals =plotptsvals(f)[2]
    if isa(domain(f),PeriodicDomain)
        plot!(plt,real([vals;vals[1]]),imag([vals;vals[1]]);opts...)
    else
        plot!(plt,real(vals),imag(vals);opts...)
    end
end


for PLOTSTR in ("complexplot","domainplot","coefficientplot")
    PLOT=parse(PLOTSTR)
    PLOTe=parse(PLOTSTR*"!")
    @eval begin
        $PLOTe(f;opts...)=$PLOTe(current(),f;opts...)
        $PLOT(f;grid=true,opts...)=$PLOTe(plot(grid=grid),f;opts...)
    end
end

##
# Special spaces
##

function plotptsvals{S<:JacobiWeight}(f::Fun{S})
    f=pad(f,3length(f)+50)
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

# TODO Fourier and Laurent spaces
# function plotptsvals{S<:Per}(f::Fun)
#     if dimension(space(f)) == Inf
#         f=pad(f,3length(f)+50)
#     else
#         f=pad(f,dimension(space(f)))
#     end
#     if length(domain(f)) < Inf
#         return points(f),values(f)
#     else
#         return points(f)[2:end],values(f)[2:end] # A hack for Fourier
#     end
# end

for (plt,TYP) in ((:(Plots.plot),:Real),(:(Plots.plot!),:Real),(:complexplot,:Complex),(:complexplot!,:Complex))
    @eval $plt{S<:Union{ArraySpace,TupleSpace},T<:$TYP}(f::Fun{S,T};opts...)=$plt(vec(f);opts...)
end

for (plt,TYP) in ((:(Plots.plot!),:Real),(:complexplot!,:Complex))
    @eval $plt{S<:Union{ArraySpace,TupleSpace},T<:$TYP}(pltin::Plots.Plot,f::Fun{S,T};opts...)=
    $plt(pltin,vec(f);opts...)
end


function Plots.plot!{S<:PiecewiseSpace,T<:Real}(plt::Plots.Plot,f::Fun{S,T};label=Void,kwds...)
    v=vec(f)
    if label == Void
        c=plt.plotargs[:color_palette][plt.n+1]
        for k=1:length(v)
            plot!(plt,v[k];color=c,kwds...)
        end
    else
        @assert length(label)==length(v)

        c=plt.plotargs[:color_palette][plt.n+1]
        for k=1:length(v)
            plot!(plt,v[k];label=label[k],color=c,kwds...)
        end
    end
    plt
end

# For dirac space, we draw a dotted line extending to infinity
function Plots.plot!{S<:DiracSpace,T<:Real}(plt::Plots.Plot,f::Fun{S,T};kwds...)
    pts=space(f).points
    n=length(pts)
    ws=pad(f.coefficients,length(pts))
    plt=plot!(plt,ones(2)*pts[1],[0,1]*ws[1];kwds...)
    c=plt.plotargs[:color_palette][plt.n]
    plot!(plt,ones(2)*pts[2:end]',[0,1]*ws[2:end]';color=c,kwds...)
    plot!(plt,ones(2)*pts',[1,2]*ws';color=c,linestyle=:dot,kwds...)
end

# for PointSpace, we draw just a line
function Plots.plot!{S<:PointSpace,T<:Real}(plt::Plots.Plot,f::Fun{S,T};kwds...)
    pts=space(f).points
    n=length(pts)
    ws=pad(f.coefficients,length(pts))
    plt=plot!(plt,ones(2)*pts[1],[0,1]*ws[1];kwds...)
    c=plt.plotargs[:color_palette][plt.n]
    plot!(plt,ones(2)*pts[2:end]',[0,1]*ws[2:end]';color=c,kwds...)
end

function Plots.plot!{S<:HeavisideSpace,T<:Real}(plt::Plots.Plot,f::Fun{S,T};kwds...)
    pts=domain(f).points
    n=length(pts)
    ws=pad(f.coefficients,dimension(space(f)))
    plt=plot!(plt,pts[1:2],ones(2)*ws[1];kwds...)
    c=plt.plotargs[:color_palette][plt.n]
    for k=2:length(ws)
        plot!(plt,ones(2)*pts[k],ws[k-1:k];color=c,linestyle=:dot,kwds...)
        plot!(plt,pts[k:k+1],ones(2)*ws[k];color=c,kwds...)
    end
    plt
end



## domainplot


Plots.plot!(plt::Plots.Plot,d::Domain;kwds...)=complexplot!(plt,Fun(identity,d);kwds...)  # default is to call complexplot
Plots.plot!(plt::Plots.Plot,d::UnionDomain;kwds...)=plot!(plt,[d.domains...];kwds...)
domainplot!(plt::Plots.Plot,f::Union{Fun,Space};kwds...)=plot!(plt,domain(f);kwds...)



## coefficientplot

coefficientplot!(plt::Plots.Plot,f::Fun;opts...)=plot!(plt,abs(f.coefficients);yscale=:log10,opts...)



## Multivariate

function Plots.plot!(plt::Plots.Plot,f::MultivariateFun;linetype=:contour,opts...)
    f=chop(f,10e-10)
    f=pad(f,max(size(f,1),20),max(size(f,2),20))
    vals=values(f)
    if norm(imag(vals),Inf)>10e-9
        warn("Imaginary part is non-neglible.  Only plotting real part.")
    end

    plot!(plt,points(space(f,1),size(vals,1)),points(space(f,2),size(vals,2)),real(vals);linetype=linetype,opts...)
end

function Plots.surface(f::MultivariateFun;opts...)
    f=chop(f,10e-10)
    f=pad(f,max(size(f,1),20),max(size(f,2),20))
    vals=values(f)
    if norm(imag(vals),Inf)>10e-9
        warn("Imaginary part is non-neglible.  Only plotting real part.")
    end

    surface(points(space(f,1),size(vals,1)),points(space(f,2),size(vals,2)),real(vals);opts...)
end




## 3D plotting
# TODO: The extra vals should only be added for periodicity?
function Plots.plot!(plt::Plots.Plot,xx::Range,yy::Range,f::MultivariateFun;opts...)
    vals      = f(xx,yy)
    #vals=[vals[:,1] vals vals[:,end]];
    #vals=[vals[1,:]; vals; vals[end,:]]
    plot!(plt,xx,yy,real(vals);opts...)
end


#plot{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS};opts...)=surf(vecpoints(f,1),vecpoints(f,2),real(values(f));opts...)
#plot(f::LowRankFun;opts...)=surf(vecpoints(f,1),vecpoints(f,2),real(values(f));opts...)
#plot(f::MultivariateFun,obj,window;opts...)=glsurfupdate(real(values(f)),obj,window;opts...)

for Plt in (:(Plots.plot),:(Plots.contour),:(Plots.surface))
    Pltex=parse(string(Plt)*"!")
    @eval begin
        $Plt{TS<:TensorSpace,T<:Real}(f::Fun{TS,T};kwds...)=$Plt(ProductFun(f);kwds...)
        $Pltex{TS<:TensorSpace,T<:Real}(f::Fun{TS,T};kwds...)=$Pltex(ProductFun(f);kwds...)
        $Pltex{TS<:TensorSpace,T<:Real}(plt::Plots.Plot,f::Fun{TS,T};kwds...)=$Pltex(plt,ProductFun(f);kwds...)
    end
end
