export domainplot, coefficientplot, complexplot


if isdir(Pkg.dir("PyPlot"))
    include("PyPlot.jl")
end

if isdir(Pkg.dir("TikzGraphs"))
    include("introspect.jl")
end



#surf(x...;opts...)=pysurf(x...;opts...)


## Fun routines

Plots.plot(f::Union{Fun,Domain,MultivariateFun};kwds...)=plot!(plot(),f;kwds...)
Plots.plot!(f::Union{Fun,Domain,MultivariateFun};kwds...)=plot!(current(),f;kwds...)

Plots.plot{F<:Union{Fun,Domain,MultivariateFun}}(v::AbstractVector{F};kwds...)=plot!(plot(),v;kwds...)
Plots.plot!{F<:Union{Fun,Domain,MultivariateFun}}(v::AbstractVector{F};kwds...)=plot!(current(),v;kwds...)


function plotptsvals(f::Fun)
    if dimension(space(f)) == Inf
        f=pad(f,3length(f)+50)
    else
        f=pad(f,dimension(space(f)))
    end
    points(f),values(f)
end
Plots.plot!{S,T<:Real}(plt::Plots.Plot,f::Fun{S,T};kwds...)=
                plot!(plt,plotptsvals(f)...;kwds...)



function Plots.plot!{F<:Union{Fun,Domain}}(plt::Plots.Plot,v::AbstractVector{F};label=Void)
    if label == Void
        for k=1:length(v)
            plot!(plt,v[k])
        end
    else
        @assert length(label)==length(v)

        for k=1:length(v)
            plot!(plt,v[k];label=label[k])
        end
    end
    plt
end





Plots.plot!{S,T<:Complex}(plt::Plots.Plot,f::Fun{S,T};label=Void)=
                plot!(plt,[real(f),imag(f)];label=(label==Void?["Real","Imag"]:label))

function complexplot!(plt::Plots.Plot,f::Fun;opts...)
    vals =plotptsvals(f)[2]
    if isa(domain(f),PeriodicDomain)
        plot!(plt,real([vals;vals[1]]),imag([vals;vals[1]]);opts...)
    else
        plot!(plt,real(vals),imag(vals);opts...)
    end
end

complexplot!(f::Fun;opts...)=complexplot!(current(),f;opts...)
complexplot(f::Fun;opts...)=complexplot!(plot(),f;opts...)


## Special spaces




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


function Plots.plot!{S<:DiracSpace,T<:Real}(plt::Plots.Plot,f::Fun{S,T};kwds...)
    pts=space(f).points
    n=length(pts)
    ws=pad(f.coefficients,length(pts))
    plt=plot!(plt,ones(2)*pts[1],[0,1]*ws[1];kwds...)
    c=plt.plotargs[:color_palette][plt.n]
    plot!(plt,ones(2)*pts[2:end]',[0,1]*ws[2:end]';color=c,kwds...)
    plot!(plt,ones(2)*pts',[1,2]*ws';color=c,linestyle=:dot,kwds...)
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


Plots.plot(d::Domain;kwds...)=complexplot(Fun(identity,d);kwds...)  # default is to call complexplot
Plots.plot(d::UnionDomain;kwds...)=plot([d.domains...];kwds...)
domainplot(f::Union{Fun,Space};kwds...)=plot(domain(f);kwds...)


## coefficientplot

coefficientplot(f::Fun;opts...)=plot(log10(abs(f.coefficients));opts...)



## Multivariate

function Plots.plot!(plt::Plots.Plot,f::MultivariateFun;linetype=:contour,opts...)
    f=chop(f,10e-10)
    f=pad(f,max(size(f,1),20),max(size(f,2),20))
    vals=values(f)
    if norm(imag(vals))>10e-9
        warn("Imaginary part is non-neglible.  Only plotting real part.")
    end

    plot!(plt,points(space(f,1),size(vals,1)),points(space(f,2),size(vals,2)),real(vals);linetype=linetype,opts...)
end

contour(f::Fun;kwds...)=plot(f,linetype=:contour,kwds...)




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
Plots.plot{TS<:TensorSpace,T<:Real}(f::Fun{TS,T};kwds...)=plot(ProductFun(f);kwds...)
