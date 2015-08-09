export setplotter, domainplot

# these defaults are overloaded as packages are loaded
plotter=@compat Dict(:contour=>"Gadfly",
    :plot=>"Gadfly",
    :surf=>"PyPlot")


function setplotter(key,val)
    global plotter
    @assert val=="PyPlot" || val =="Gadfly" || val =="GLPlot"
    plotter[key]=val
    plotter
end

function setplotter(str)
    if str=="PyPlot"
        setplotter(:contour,str)
        setplotter(:plot,str)
        setplotter(:surf,str)
    elseif str == "GLPlot"
        setplotter(:surf,str)
    else
        setplotter(:contour,str)
        setplotter(:plot,str)
    end
end


if isdir(Pkg.dir("GLPlot"))
    include("GLPlot.jl")
    setplotter("GLPlot")
end
if isdir(Pkg.dir("PyPlot"))
    include("PyPlot.jl")
    setplotter("PyPlot")
end

include("Gadfly.jl")

if isdir(Pkg.dir("Gadfly"))
    setplotter("Gadfly")
end
if isdir(Pkg.dir("TikzGraphs"))
    include("introspect.jl")
end

function plot(opts...;kwds...)
    if plotter[:plot]=="Gadfly"
        gadflyplot(opts...;kwds...)
    elseif plotter[:plot]=="PyPlot"
        pyplot(opts...;kwds...)
    else
        error("Plotter " * plotter[:plot] * " not supported.")
    end
end

function layer(opts...)
    if plotter[:plot]=="Gadfly"
        gadflylayer(opts...)
    elseif plotter[:plot]=="PyPlot"
        pyplot(opts...)
    else
        error("Plotter " * plotter[:plot] * " not supported.")
    end
end

function contour(x,y::Vector,z::Array,v...;opts...)
    if plotter[:contour]=="Gadfly"
        gadflycontour(x,y,z,v...;opts...)
    elseif plotter[:contour]=="PyPlot"
        pycontour(x,y,z,v...;opts...)
    else
        error("Plotter " * plotter[:contour] * " not supported.")
    end
end

function surf(x...;opts...)
    if plotter[:surf]=="GLPlot"
        glsurf(x...;opts...)
    elseif plotter[:surf]=="PyPlot"
        pysurf(x...;opts...)
    else
        error("Plotter " * plotter[:surf] * " not supported.")
    end
end



## Fun routines


for OP in (:plot,:layer)
    @eval begin
        function $OP{S,T<:Real}(f::Fun{S,T},v...;opts...)
            f=pad(f,3length(f)+50)
            $OP(points(f),values(f),v...;opts...)
        end

        function $OP{S,T<:Complex}(f::Fun{S,T},v...;opts...)
            f=pad(f,3length(f)+50)
            $OP(points(f),values(f),v...;opts...)
        end

        function $OP{F<:Fun}(f::Vector{F},v...;opts...)
            if plotter[:plot]=="PyPlot"
                for k=1:length(f)
                    $OP(f[k],v...;opts...)
                end
            else
                n=3mapreduce(length,max,f)+50
                m=length(f)
                X,Y=Array(Float64,n,m),Array(Float64,n,m)
                for k=1:length(f)
                    X[:,k]=points(space(f[k]),n)
                    Y[:,k]=values(pad(f[k],n))
                end
                plot(X,Y,v...;opts...)
            end
        end

        function $OP{S}(r::Range,f::Fun{S,Float64},v...;opts...)
            $OP([r],f[[r]],v...;opts...)
        end
    end
end

function complexplot(f::Fun,v...;opts...)
    f=pad(f,3length(f)+50)
    vals =values(f)
    d = domain(f)
    if isa(d,Circle)
        plot(real([vals,vals[1]]),imag([vals,vals[1]]),v...;opts...)
    else
        plot(real(vals),imag(vals),v...;opts...)
    end
end

function complexplot{F<:Fun}(f::Vector{F},v...;opts...)
    if plotter[:plot]=="PyPlot"
        for k=1:length(f)
            complexplot(f[k],v...;opts...)
        end
    else
        plot(mapreduce(complexlayer,vcat,f),v...;opts...)
    end
end

function complexlayer(f::Fun,v...;opts...)
    f=pad(f,3length(f)+50)
    vals =values(f)
    d = domain(f)
    if isa(d,Circle)
        layer(real([vals,vals[1]]),imag([vals,vals[1]]),v...;opts...)
    else
        layer(real(vals),imag(vals),v...;opts...)
    end
end

function complexlayer{F<:Fun}(f::Vector{F},v...;opts...)
    typeof(complexlayer(f[1],v...;opts...)[1])[complexlayer(f[k],v...;opts...)[1] for k=1:length(f)]
end

for (plt,TYP) in ((:plot,:Real),(:complexplot,:Complex),(:layer,:Real),(:complexlayer,:Complex))
    @eval $plt{S<:Union(PiecewiseSpace,ArraySpace),T<:$TYP}(f::Fun{S,T},v...;opts...)=$plt(vec(f),v...;opts...)
end



## Multivariate

function contour(f::MultivariateFun,v...;opts...)
    f=chop(f,10e-10)
    f=pad(f,max(size(f,1),20),max(size(f,2),20))
    vals=values(f)
    if norm(imag(vals))>10e-9
        warn("Imaginary part is non-neglible.  Only plotting real part.")
    end

    contour(points(space(f,1),size(vals,1)),points(space(f,2),size(vals,2)),real(vals),v...;opts...)
end

contour(f::Fun,v...;opts...)=contour(ProductFun(f),v...;opts...)



## 3D plotting
# TODO: The extra vals should only be added for periodicity?
function plot(xx::Range,yy::Range,f::MultivariateFun,v...;opts...)
    vals      = evaluate(f,xx,yy)
    #vals=[vals[:,1] vals vals[:,end]];
    #vals=[vals[1,:]; vals; vals[end,:]]
    surf([xx],[yy],real(vals),v...;opts...)
end

function plot(xx::Range,yy::Range,f::MultivariateFun,obj,window)
    vals      = evaluate(f,xx,yy)
    #vals=[vals[:,1] vals vals[:,end]];
    #vals=[vals[1,:]; vals; vals[end,:]]
    glsurfupdate(real(vals),obj,window)
end


plot(f::MultivariateFun;opts...)=surf(points(f,1),points(f,2),real(values(f));opts...)
plot{S,V,SS<:TensorSpace}(f::ProductFun{S,V,SS};opts...)=surf(vecpoints(f,1),vecpoints(f,2),real(values(f));opts...)
plot(f::LowRankFun;opts...)=surf(vecpoints(f,1),vecpoints(f,2),real(values(f));opts...)
plot(f::MultivariateFun,obj,window;opts...)=glsurfupdate(real(values(f)),obj,window;opts...)
plot{TS<:TensorSpace,T<:Real}(f::Fun{TS,T};opts...)=plot(ProductFun(f);ops...)


# plot{S<:IntervalSpace,V<:PeriodicSpace,SS<:TensorSpace}(f::ProductFun{S,V,SS};opts...)=surf(vecpoints(f,1),vecpoints(f,2),real(values(f));opts...)
# function plot{S<:IntervalSpace,V<:PeriodicSpace}(f::ProductFun{S,V};opts...)
#     Px,Py=points(f)
#     vals=real(values(f))
#     surf([Px Px[:,1]], [Py Py[:,1]], [vals vals[:,1]];opts...)
# end
# function plot{S<:IntervalSpace,V<:PeriodicSpace}(f::ProductFun{S,V},obj,window)
#     vals=real(values(f))
#     glsurfupdate([vals vals[:,1]],obj,window)
# end


function plot{DS<:DiracSpace,T<:Real}(f::Fun{DS,T},v...)
    n=length(space(f).points)
    plot(layer(Fun(f.coefficients[n+1:end],space(f).space)),
               map(gadflydeltalayer,space(f).points,f.coefficients[1:n])...,v...)
end

## domainplot


for (plt,cplt) in ((:plot,:complexplot),(:layer,:complexlayer))
    @eval $plt(d::Domain,v...;kwds...)=$cplt(Fun(identity,d),v...;kwds...)  # default is to call complexplot
end
function plot{D<:Domain}(f::Vector{D},v...;opts...)
    if plotter[:plot]=="PyPlot"
        for k=1:length(f)
            plot(f[k],v...;opts...)
        end
    else
        plot(layer(f),v...;opts...)
    end
end

layer{D<:Domain}(d::Vector{D})=mapreduce(layer,vcat,d)

for OP in (:plot,:layer)
    @eval $OP(d::UnionDomain)=$OP(d.domains)
end


domainplot(f::Union(Fun,FunctionSpace),v...;kwds...)=plot(domain(f),v...;kwds...)
domainlayer(f::Union(Fun,FunctionSpace))=layer(domain(f))




