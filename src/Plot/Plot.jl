export domainplot, coefficientplot, complexplot


if isdir(Pkg.dir("PyPlot"))
    include("PyPlot.jl")
end

if isdir(Pkg.dir("TikzGraphs"))
    include("introspect.jl")
end



surf(x...;opts...)=pysurf(x...;opts...)


## Fun routines



function plotptsvals(f::Fun)
    f=pad(f,3length(f)+50)
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

for FUNC in (:(Plots.plot),:(Plots.plot!))
    @eval begin
        $FUNC{S,T<:Real}(f::Fun{S,T};kwds...)=$FUNC(plotptsvals(f)...;kwds...)

        function $FUNC{F<:Union{Fun,Domain}}(v::AbstractVector{F};label=Void)
            if label == Void
                plt=$FUNC(first(v))
                if length(v)>1
                    plot!(plt,v[2:end])
                end
            else
                @assert length(label)==length(v)
                plt=$FUNC(first(v);label=first(label))
                if length(v)>1
                    plot!(plt,v[2:end];label=label[2:end])
                end
            end
            plt
        end
    end
end

Plots.plot{S,T<:Complex}(f::Fun{S,T};label=Void)=
                plot([real(f),imag(f)];label=(label==Void?["Real","Imag"]:label))


function complexplot(f::Fun;opts...)
    vals =plotptsvals(f)[2]
    if isa(domain(f),PeriodicDomain)
        plot(real([vals;vals[1]]),imag([vals;vals[1]]);opts...)
    else
        plot(real(vals),imag(vals);opts...)
    end
end

function complexplot!(f::Fun;opts...)
    vals =plotptsvals(f)[2]
    if isa(domain(f),PeriodicDomain)
        plot!(real([vals;vals[1]]),imag([vals;vals[1]]);opts...)
    else
        plot!(real(vals),imag(vals);opts...)
    end
end

function complexplot!(plt::Plots.Plot,f::Fun;opts...)
    vals =plotptsvals(f)[2]
    if isa(domain(f),PeriodicDomain)
        plot!(plt,real([vals;vals[1]]),imag([vals;vals[1]]);opts...)
    else
        plot!(plt,real(vals),imag(vals);opts...)
    end
end


## Special spaces

for (plt,TYP) in ((:(Plots.plot),:Real),(:(Plots.plot!),:Real),(:complexplot,:Complex),(:complexplot!,:Complex))
    @eval $plt{S<:Union{PiecewiseSpace,ArraySpace,TupleSpace},T<:$TYP}(f::Fun{S,T};opts...)=$plt(vec(f);opts...)
end

for (plt,TYP) in ((:(Plots.plot!),:Real),(:complexplot!,:Complex))
    @eval $plt{S<:Union{PiecewiseSpace,ArraySpace,TupleSpace},T<:$TYP}(pltin::Plots.Plot,f::Fun{S,T};opts...)=
    $plt(pltin,vec(f);opts...)
end


for OP in (:(Plots.plot),:(Plots.plot!))
    @eval function $OP{S<:DiracSpace,T<:Real}(f::Fun{S,T};kwds...)
        pts=space(f).points
        n=length(pts)
        ws=pad(f.coefficients,length(pts))
        plt=$OP(ones(2)*pts[1],[0,1]*ws[1];kwds...)
        c=plt.plotargs[:color_palette][1]
        plot!(plt,ones(2)*pts[2:end]',[0,1]*ws[2:end]';color=c,kwds...)
        plot!(plt,ones(2)*pts',[1,2]*ws';color=c,linestyle=:dot,kwds...)
    end
end

function Plots.plot!{S<:DiracSpace,T<:Real}(plt::Plots.Plot,f::Fun{S,T};kwds...)
    pts=space(f).points
    n=length(pts)
    ws=pad(f.coefficients,length(pts))
    plt=$OP(plt,ones(2)*pts[1],[0,1]*ws[1];kwds...)
    c=plt.plotargs[:color_palette][1]
    plot!(plt,ones(2)*pts[2:end]',[0,1]*ws[2:end]';color=c,kwds...)
    plot!(plt,ones(2)*pts',[1,2]*ws';color=c,linestyle=:dot,kwds...)
end



## domainplot


Plots.plot(d::Domain;kwds...)=complexplot(Fun(identity,d);kwds...)  # default is to call complexplot
Plots.plot(d::UnionDomain;kwds...)=plot([d.domains...];kwds...)
domainplot(f::Union{Fun,Space};kwds...)=plot(domain(f);kwds...)


## coefficientplot

coefficientplot(f::Fun;opts...)=plot(log10(abs(f.coefficients));opts...)



## Multivariate

for FUNC in (:(Plots.plot),:(Plots.plot!))
    @eval begin
        function $FUNC(f::MultivariateFun;linetype=:contour,opts...)
            f=chop(f,10e-10)
            f=pad(f,max(size(f,1),20),max(size(f,2),20))
            vals=values(f)
            if norm(imag(vals))>10e-9
                warn("Imaginary part is non-neglible.  Only plotting real part.")
            end

            plot(points(space(f,1),size(vals,1)),points(space(f,2),size(vals,2)),real(vals);linetype=linetype,opts...)
        end
        $FUNC(f::Fun;opts...)=$FUNC(ProductFun(f);opts...)
    end
end

contour(f::Fun;kwds...)=plot(f,linetype=:contour,kwds...)




## 3D plotting
# TODO: The extra vals should only be added for periodicity?
for OP in (:(Plots.plot),:(Plots.plot!))
    @eval function $OP(xx::Range,yy::Range,f::MultivariateFun;opts...)
        vals      = f(xx,yy)
        #vals=[vals[:,1] vals vals[:,end]];
        #vals=[vals[1,:]; vals; vals[end,:]]
        $OP(collect(xx),collect(yy),real(vals);opts...)
    end
end

function Plots.plot!(plt::Plots.Plot,xx::Range,yy::Range,f::MultivariateFun;opts...)
    vals      = f(xx,yy)
    #vals=[vals[:,1] vals vals[:,end]];
    #vals=[vals[1,:]; vals; vals[end,:]]
    plot!(plt,collect(xx),collect(yy),real(vals);opts...)
end

# function plot(xx::Range,yy::Range,f::MultivariateFun,obj,window)
#     vals      = evaluate(f,xx,yy)
#     #vals=[vals[:,1] vals vals[:,end]];
#     #vals=[vals[1,:]; vals; vals[end,:]]
#     glsurfupdate(real(vals),obj,window)
# end

#plot{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS};opts...)=surf(vecpoints(f,1),vecpoints(f,2),real(values(f));opts...)
#plot(f::LowRankFun;opts...)=surf(vecpoints(f,1),vecpoints(f,2),real(values(f));opts...)
#plot(f::MultivariateFun,obj,window;opts...)=glsurfupdate(real(values(f)),obj,window;opts...)
Plots.plot{TS<:TensorSpace,T<:Real}(f::Fun{TS,T};kwds...)=plot(ProductFun(f);kwds...)
