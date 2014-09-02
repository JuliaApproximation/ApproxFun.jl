## Plotting


using DataFrames  
export plot
export complexplot,contour

function gplot(xx::Vector,yy::Vector;axis=-1)
    require("Gadfly")
    if axis==-1
        Main.Gadfly.plot(x=xx, y=yy, Main.Gadfly.Geom.line)
    else
        Main.Gadfly.plot(x=xx, y=yy, Main.Gadfly.Geom.line, Main.Gadfly.Scale.y_continuous(minvalue=axis[1],maxvalue=axis[2]))
    end    
end



function cplot{T<:Real}(x::Vector{T},y::Vector{Complex{Float64}};axis=-1)
    r=real(y)
    i=imag(y)
    if axis==-1
        plot(DataFrames.DataFrame({[x,x],[r,i],[fill("Re",length(x)),fill("Im",length(x))]},DataFrames.Index([:x=>1,:y=>2,:Function=>3],[:x,:y,:Function])),
        x="x",y="y",color="Function",Geom.line)
    else
        plot(DataFrames.DataFrame({[x,x],[r,i],[fill("Re",length(x)),fill("Im",length(x))]},DataFrames.Index([:x=>1,:y=>2,:Function=>3],[:x,:y,:Function])),
        x="x",y="y",color="Function",Geom.line, Scale.y_continuous(minvalue=axis[1],maxvalue=axis[2]))      
    end
end


function plot{T<:Real}(f::IFun{T};opts...)
    f=pad(f,3length(f)+50)
    gplot(points(f),values(f);opts...)
end

function plot{T<:Complex}(f::IFun{T};opts...)
    f=pad(f,3length(f)+50)
    cplot(points(f),values(f);opts...)
end
