## Plotting


export plot
export complexplot,contour


## Vector routines
function plot(xx::Vector,yy::Vector;axis=-1)
    require("Gadfly")
    if axis==-1
        Main.Gadfly.plot(x=xx, y=yy, Main.Gadfly.Geom.path)
    elseif length(axis)==2
        Main.Gadfly.plot(x=xx, y=yy, Main.Gadfly.Geom.path, Main.Gadfly.Scale.y_continuous(minvalue=axis[1],maxvalue=axis[2]))
    else #length(axis)==4
        Main.Gadfly.plot(x=xx, y=yy, Main.Gadfly.Geom.path, Main.Gadfly.Scale.x_continuous(minvalue=axis[1],maxvalue=axis[2]),Main.Gadfly.Scale.y_continuous(minvalue=axis[3],maxvalue=axis[4]))    
    end    
end



function plot{T<:Real}(x::Vector{T},y::Vector{Complex{Float64}};axis=-1)
    require("Gadfly")
    require("DataFrames")    
    r=real(y)
    i=imag(y)
    if axis==-1
        Main.Gadfly.plot(Main.DataFrames.DataFrame({[x,x],[r,i],[fill("Re",length(x)),fill("Im",length(x))]},Main.DataFrames.Index([:x=>1,:y=>2,:Function=>3],[:x,:y,:Function])),
        x="x",y="y",color="Function",Main.Gadfly.Geom.path)
    else
        Main.Gadfly.plot(Main.DataFrames.DataFrame({[x,x],[r,i],[fill("Re",length(x)),fill("Im",length(x))]},Main.DataFrames.Index([:x=>1,:y=>2,:Function=>3],[:x,:y,:Function])),
        x="x",y="y",color="Function",Main.Gadfly.Geom.path, Main.Gadfly.Scale.y_continuous(minvalue=axis[1],maxvalue=axis[2]))      
    end
end

function contour(x::Vector,y::Vector,z::Matrix;levels=-1)
    require("Gadfly")
    if levels==-1
        Main.Gadfly.plot(x=x,y=y,z=z,Main.Gadfly.Geom.contour)
    else
        Main.Gadfly.plot(x=x,y=y,z=z,Main.Gadfly.Geom.contour(levels=levels))
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
