export set2Dplotter, set3Dplotter,setplotter

if isdir(Pkg.dir("Gadfly"))
    include("Gadfly.jl")
end
if isdir(Pkg.dir("GLPlot"))
    include("GLPlot.jl")
end
if isdir(Pkg.dir("PyPlot"))
    include("PyPlot.jl")
end


plotter2D="Gadfly"
plotter3D="GLPlot"

set2Dplotter(str)=(global plotter2D=str)
set3Dplotter(str)=(global plotter3D=str)
function setplotter(str)
    if str=="PyPlot"
        set2Dplotter(str)
        set3Dplotter(str)
    elseif str == "GLPlot"
        set3Dplotter(str)
    else
        set2Dplotter(str)
    end
end


function plot(x,y::Vector;opts...)
    if plotter2D=="Gadfly"
        gadflyplot(x,y;opts...)
    elseif plotter2D=="PyPlot"
        pyplot(x,y;opts...)
    else
        error("Plotter " * plotter2D * " not supported.")
    end
end

function contour(x,y::Vector,z::Array;opts...)
    if plotter2D=="Gadfly"
        gadflycontour(x,y,z;opts...)
    elseif plotter2D=="PyPlot"
        pycontour(x,y,z;opts...)
    else
        error("Plotter " * plotter2D * " not supported.")
    end
end

function surf(x::Vector,y::Vector,z::Array;opts...)
    if plotter3D=="GLPlot"
        glsurf(x,y;opts...)
    elseif plotter3D=="PyPlot"
        pysurf(x,y;opts...)
    else
        error("Plotter " * plotter3D * " not supported.")
    end
end



## Fun routines


function plot{T<:Real}(f::Fun{T};opts...)
    f=pad(f,3length(f)+50)
    plot(points(f),values(f);opts...)
end

function plot{T<:Complex}(f::Fun{T};opts...)
    f=pad(f,3length(f)+50)
    plot(points(f),values(f);opts...)
end


function plot(r::Range,f::Fun{Float64};opts...)
    plot(r,f[[r]];opts...)
end

function complexplot(f::Fun{Complex{Float64}};opts...) 
    f=pad(f,3length(f)+50)
    vals =values(f)

    plot(real(vals),imag(vals);opts...)
end


##FFun

# function plot(f::FFun;opts...) 
#     f=deepcopy(f)
#     
#     m=max(-firstindex(f.coefficients),lastindex(f.coefficients))
#     
#     f.coefficients=pad(f.coefficients,-m:m)
# 
#     pts = [points(f),fromcanonical(f,π)]
#     vals =[values(f),first(values(f))]
# 
#     plot(pts,vals;opts...)
# end
# 
# 
# function complexplot(f::FFun{Complex{Float64}};opts...) 
#     pts = [points(f),fromcanonical(f,π)]
#     vals =[values(f),first(values(f))]
# 
#     plot(real(vals),imag(vals);opts...)
# end


## SingFun

##TODO: reimplement
# function plot(f::SingFun;opts...) 
#     pf = pad(f,3length(f)+100)
#     
#     if f.α >= 0 && f.β >= 0
#         plot(points(pf),values(pf);opts...)
#     elseif f.α >= 0
#         plot(points(pf)[1:end-1],values(pf)[1:end-1];opts...)    
#     elseif f.β >= 0    
#         plot(points(pf)[2:end],values(pf)[2:end];opts...)    
#     else
#         plot(points(pf)[2:end-1],values(pf)[2:end-1];opts...)
#     end
# end



## Multivariate

function contour(f::MultivariateFun;opts...)
    f=chop(f,10e-10)
    contour(points(f,1),points(f,2),values(f);opts...)
end




## 3D plotting

function plot(xx::Range,yy::Range,f::MultivariateFun)
    vals      = evaluate(f,xx,yy)
    vals=[vals[:,1] vals vals[:,end]];
    vals=[vals[1,:]; vals; vals[end,:]]    
    surf(vals)    
end

function plot(xx::Range,yy::Range,f::MultivariateFun,obj,window)
    vals      = evaluate(f,xx,yy)
    vals=[vals[:,1] vals vals[:,end]];
    vals=[vals[1,:]; vals; vals[end,:]]    
    surf(vals,obj,window)    
end

plot(f::MultivariateFun)=surf(points(f,1),points(f,2),values(f))
