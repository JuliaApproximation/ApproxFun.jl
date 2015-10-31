## Plotting

# PyPlot is already imported

function pysemilogy{N<:Real}(xx,yy::Array{N};axis=-1,opts...)
    PyPlot.semilogy(xx,yy,color="blue";opts...)
    if axis!=-1
        if length(axis) == 4
            PyPlot.axis(axis)
        else
            PyPlot.axis([xx[1];xx[end];axis])
        end
    end
end


function pyplot{N<:Real}(xx,yy::Array{N};axis=-1,opts...)
    PyPlot.plot(xx,yy,color="blue";opts...)
    if axis!=-1
        if length(axis) == 4
            PyPlot.axis(axis)
        else
            PyPlot.axis([xx[1];xx[end];axis])
        end
    end
end

function pyplot{N<:Complex}(xx,yy::Array{N};axis=-1,opts...)
    PyPlot.plot(xx,real(yy),color="blue";opts...)
    PyPlot.plot(xx,imag(yy),color="red";opts...)

    if axis!=-1
        if length(axis) == 4
            PyPlot.axis(axis)
        else
            PyPlot.axis([xx[1];xx[end];axis])
        end
    end
end

# function complexplot(f::IFun{Complex{Float64}})
#     pf = pad(f,4length(f))
#     vals =values(pf)
#
#     PyPlot.plot(real(vals),imag(vals))
#     PyPlot.arrow(real(vals[end-1]),imag(vals[end-1]),real(vals[end]-vals[end-1]),imag(vals[end]-vals[end-1]),width=.01,edgecolor="white")
# end
#
# function complexplot(f::FFun{Complex{Float64}})
#     pts = [points(f),fromcanonical(f,π)]
#     vals =[values(f),first(values(f))]
#
#     PyPlot.plot(real(vals),imag(vals))
#     PyPlot.arrow(real(vals[end-1]),imag(vals[end-1]),real(vals[end]-vals[end-1]),imag(vals[end]-vals[end-1]),width=.01,edgecolor="white")
# end


##Plotting

#TODO: Pad
#
# function PyPlot.plot(f::FFun;axis=-1)
#     f=deepcopy(f)
#
#     m=max(-firstindex(f.coefficients),lastindex(f.coefficients))
#
#     f.coefficients=pad(f.coefficients,-m:m)
#
#     pts = [points(f),fromcanonical(f,π)]
#     vals =[values(f),first(values(f))]
#
#     PyPlot.plot(pts,real(vals))
#     PyPlot.plot(pts,imag(vals),color="red")
#
#     if axis!=-1
#         if length(axis) == 4
#             PyPlot.axis(axis)
#         else
#             PyPlot.axis([f.domain.a,f.domain.b,axis])
#         end
#     end
# end
#


##2D

function pysurf(x::Vector,y::Vector,z; rstride=2,cstride=2,kwds...)
    PyPlot.surf(x,y,z.';linewidth=0,rstride=rstride,cstride=cstride,kwds...)
end

function pycontour(x::Vector,y::Vector,z,kwds...;opts...)
    PyPlot.contour(x,y,z.',kwds...;opts...)
end

function pycontourf(x::Vector,y::Vector,z,kwds...;opts...)
    PyPlot.contourf(x,y,z.',kwds...;opts...)
end


function pysurf(x::Matrix,y::Matrix,z; rstride=2,cstride=2,kwds...)
    PyPlot.surf(x,y,z;linewidth=0,rstride=rstride,cstride=cstride,kwds...)
end

function pycontour(x::Matrix,y::Matrix,z,kwds...;opts...)
    PyPlot.contour(x,y,z,kwds...;opts...)
end

function pycontourf(x::Matrix,y::Matrix,z,kwds...;opts...)
    PyPlot.contourf(x,y,z,kwds...;opts...)
end

function pyhist(a...;kwds...)
    PyPlot.plt[:hist](a...;kwds...)
end

## SingFun

# function PyPlot.plot(f::SingFun)
#     pf = pad(f,3length(f)+100)
#
#     if f.α >= 0 && f.β >= 0
#         PyPlot.plot(points(pf),values(pf))
#     elseif f.α >= 0
#         PyPlot.plot(points(pf)[1:end-1],values(pf)[1:end-1])
#     elseif f.β >= 0
#         PyPlot.plot(points(pf)[2:end],values(pf)[2:end])
#     else
#         PyPlot.plot(points(pf)[2:end-1],values(pf)[2:end-1])
#     end
# end


## ArrayFun

# function PyPlot.plot{T<:AbstractFun}(v::Array{T}; axis=-1)
#     if axis != -1
#         PyPlot.axis([v[1].domain.a,v[end].domain.b,axis])
#     end
#     for f in v
#         plot(f)
#     end
# end
