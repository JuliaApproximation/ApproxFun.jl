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
if isdir(Pkg.dir("Gadfly"))
    include("Gadfly.jl")
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

function contour(x,y::Vector,z::Array;opts...)
    if plotter[:contour]=="Gadfly"
        gadflycontour(x,y,z;opts...)
    elseif plotter[:contour]=="PyPlot"
        pycontour(x,y,z;opts...)
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


function plot{S,T<:Real}(f::Fun{S,T};opts...)
    f=pad(f,3length(f)+50)
    plot(points(f),values(f);opts...)
end

function plot{S,T<:Complex}(f::Fun{S,T};opts...)
    f=pad(f,3length(f)+50)
    plot(points(f),values(f);opts...)
end

function plot{F<:Fun}(f::Vector{F};opts...)
    n=3mapreduce(length,max,f)+50
    vals=Array(Float64,n,length(f))
    pts=Array(Float64,n,length(f))
    for k=1:length(f)
        pf=pad(f[k],n)
        vals[:,k]=values(pf)
        pts[:,k]=points(pf)
    end
    plot(pts,vals;opts...)
end




function plot{S}(r::Range,f::Fun{S,Float64};opts...)
    plot([r],f[[r]];opts...)
end

function complexplot(f::Fun;opts...) 
    f=pad(f,3length(f)+50)
    vals =values(f)

    plot(real(vals),imag(vals);opts...)
end

function complexplot{F<:Fun}(f::Vector{F};opts...)
    n=3mapreduce(length,max,f)+50
    vals=Array(Float64,n,length(f))
    pts=Array(Float64,n,length(f))
    for k=1:length(f)
        pf=values(pad(f[k],n))
        pts[:,k]=real(pf)
        vals[:,k]=imag(pf)        
    end
    plot(pts,vals;opts...)
end

function complexlayer(f::Fun;opts...) 
    f=pad(f,3length(f)+50)
    vals =values(f)

    layer(real(vals),imag(vals);opts...)
end

for plt in (:plot,:complexplot)
    @eval $plt{S<:Union(PiecewiseSpace,ArraySpace),T<:Real}(f::Fun{S,T};opts...)=$plt(vec(f);opts...)
end


## Multivariate

function contour(f::MultivariateFun;opts...)
    f=chop(f,10e-10)
    f=pad(f,max(size(f,1),20),max(size(f,2),20))
    vals=values(f)
    if norm(imag(vals))>10e-9
        warn("Imaginary part is non-neglible.  Only plotting real part.")
    end
    
    contour(vecpoints(f,1),vecpoints(f,2),real(vals);opts...)
end




## 3D plotting
function plot(xx::Range,yy::Range,f::MultivariateFun)
    vals      = evaluate(f,xx,yy)
    vals=[vals[:,1] vals vals[:,end]];
    vals=[vals[1,:]; vals; vals[end,:]]    
    surf(real(vals))    
end

function plot(xx::Range,yy::Range,f::MultivariateFun,obj,window)
    vals      = evaluate(f,xx,yy)
    vals=[vals[:,1] vals vals[:,end]];
    vals=[vals[1,:]; vals; vals[end,:]]    
    glsurfupdate(real(vals),obj,window)    
end


plot(f::MultivariateFun)=surf(points(f,1),points(f,2),real(values(f)))
plot{S,V,T}(f::TensorFun{S,V,T})=surf(vecpoints(f,1),vecpoints(f,2),real(values(f)))
plot(f::LowRankFun)=surf(vecpoints(f,1),vecpoints(f,2),real(values(f)))
plot(f::MultivariateFun,obj,window)=glsurfupdate(real(values(f)),obj,window)


plot{S<:IntervalSpace,V<:PeriodicSpace,T}(f::TensorFun{S,V,T})=surf(vecpoints(f,1),vecpoints(f,2),real(values(f)))
function plot{S<:IntervalSpace,V<:PeriodicSpace}(f::ProductFun{S,V})
    Px,Py=points(f)
    vals=real(values(f))
    surf([Px Px[:,1]], [Py Py[:,1]], [vals vals[:,1]])
end
function plot{S<:IntervalSpace,V<:PeriodicSpace}(f::ProductFun{S,V},obj,window)
    vals=real(values(f))
    glsurfupdate([vals vals[:,1]],obj,window)
end



## domainplot


for (plt,cplt) in ((:plot,:complexplot),(:layer,:complexlayer))
    @eval $plt(d::Domain;kwds...)=$cplt(Fun(identity,d);kwds...)  # default is to call complexplot
end
plot{D<:Domain}(ds::Vector{D};kwds...)=complexplot(map(d->Fun(identity,d),ds);kwds...)
layer{D<:Domain}(d::Vector{D})=map(layer,d)
        
domainplot(f::Union(Fun,FunctionSpace);kwds...)=plot(domain(f);kwds...)
domainlayer(f::Union(Fun,FunctionSpace))=layer(domain(f))
