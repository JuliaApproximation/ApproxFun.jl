## Plotting

using PyPlot

export complexplot


function plot(f::IFun{Float64}) 
    pf = pad(f,3length(f))
    PyPlot.plot(points(pf),values(pf))
end

function plot(f::IFun{Complex{Float64}}) 
    pf = pad(f,3length(f))
    pts = points(pf)
    vals =values(pf)

    PyPlot.plot(pts,real(vals))
    PyPlot.plot(pts,imag(vals),color="red")
end

function complexplot(f::IFun{Complex{Float64}}) 
    pf = pad(f,4length(f))
    vals =values(pf)

    PyPlot.plot(real(vals),imag(vals))
    PyPlot.arrow(real(vals[end-1]),imag(vals[end-1]),real(vals[end]-vals[end-1]),imag(vals[end]-vals[end-1]),width=.01,edgecolor="white")    
end

function complexplot(f::FFun{Complex{Float64}}) 
    pts = [points(f),first(points(f))]
    vals =[values(f),first(values(f))]

    PyPlot.plot(real(vals),imag(vals))
    PyPlot.arrow(real(vals[end-1]),imag(vals[end-1]),real(vals[end]-vals[end-1]),imag(vals[end]-vals[end-1]),width=.01,edgecolor="white")    
end


##Plotting

#TODO: Pad

function plot(f::FFun) 
    pts = [points(f),first(points(f))]
    vals =[values(f),first(values(f))]

    PyPlot.plot(pts,real(vals))
    PyPlot.plot(pts,imag(vals),color="red")
end
# 


##2D

plot(f::Fun2D)=PyPlot.surf(points(f,1),points(f,2),values(f)')



## SingFun

function plot(f::SingFun) 
    pf = pad(f,3length(f)+100)
    
    if f.α >= 0 && f.β >= 0
        PyPlot.plot(points(pf),values(pf))
    elseif f.α >= 0
        PyPlot.plot(points(pf)[1:end-1],values(pf)[1:end-1])    
    elseif f.β >= 0    
        PyPlot.plot(points(pf)[2:end],values(pf)[2:end])    
    else
        PyPlot.plot(points(pf)[2:end-1],values(pf)[2:end-1])
    end
end


## ArrayFun

function plot{T<:AbstractFun}(v::Vector{T})
    for f in v
        plot(f)
    end
end
