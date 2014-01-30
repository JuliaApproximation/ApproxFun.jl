## Plotting

using PyPlot


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


##Plotting

#TODO: Pad

function plot(f::FFun) 
    pts = [points(f),f.domain.b]
    vals =[values(f),first(values(f))]

    PyPlot.plot(pts,real(vals))
    PyPlot.plot(pts,imag(vals),color="red")
end
# 


##2D

function plot(f::Fun2D)
    xm=3mapreduce(length,max,f.A);
    ym=3mapreduce(length,max,f.B);
    ptsx=points(first(f.A).domain,xm);
    ptsy=points(first(f.B).domain,ym);
    ret=zeros(xm,ym);
    for k=1:length(f.A)
        ret+=values(pad(f.A[k],xm))*values(pad(f.B[k],ym))'
    end
    
    PyPlot.surf(ptsx,ptsy,ret')
end