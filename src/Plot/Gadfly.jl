## Plotting

export domainplot


## Vector routines
function gadflyplot{T<:Real}(xx::Vector{T},yy::Vector;axis=-1)
    require("Gadfly")
    if axis==-1
        Main.Gadfly.plot(x=xx, y=yy, Main.Gadfly.Geom.path)
    elseif length(axis)==2
        Main.Gadfly.plot(x=xx, y=yy, Main.Gadfly.Geom.path, Main.Gadfly.Scale.y_continuous(minvalue=axis[1],maxvalue=axis[2]))
    else #length(axis)==4
        Main.Gadfly.plot(x=xx, y=yy, Main.Gadfly.Geom.path, Main.Gadfly.Scale.x_continuous(minvalue=axis[1],maxvalue=axis[2]),Main.Gadfly.Scale.y_continuous(minvalue=axis[3],maxvalue=axis[4]))    
    end    
end

function gadflylayer{T<:Real}(xx::Vector{T},yy::Vector)
    require("Gadfly")

    Main.Gadfly.layer(x=xx, y=yy, Main.Gadfly.Geom.path)
end


function gadflyplot{T<:Complex}(xx::Vector{T},yy::Vector;opts...)
    warn("Complex domain, projecting to real axis")
    gadflyplot(real(xx),yy;opts...)
end

function gadflyplot{T<:Real}(x::Vector{T},y::Vector{Complex{Float64}};axis=-1)
    require("Gadfly")
    require("DataFrames")    
    r=real(y)
    i=imag(y)
    
    dat=Main.DataFrames.DataFrame(Any[[x,x],[r,i],[fill("Re",length(x)),fill("Im",length(x))]],Main.DataFrames.Index((@compat Dict(:x=>1,:y=>2,:Function=>3)),
            [:x,:y,:Function]))    
    if axis==-1
        Main.Gadfly.plot(dat,x="x",y="y",color="Function",Main.Gadfly.Geom.path)
    else
        Main.Gadfly.plot(dat,x="x",y="y",color="Function",Main.Gadfly.Geom.path, Main.Gadfly.Scale.y_continuous(minvalue=axis[1],maxvalue=axis[2]))      
    end
end

#Plot multiple contours
function gadflyplot{T<:Real,V<:Real}(x::Matrix{T},y::Matrix{V};axis=-1)
    require("Gadfly")
    require("DataFrames")    

    dat=Main.DataFrames.DataFrame(Any[vec(x),vec(y),[[fill(string(k),size(x,1)) for k=1:size(y,2)]...]],Main.DataFrames.Index((@compat Dict(:x=>1,:y=>2,:Function=>3)),
            [:x,:y,:Function]))

    if axis==-1
        Main.Gadfly.plot(dat,x="x",y="y",color="Function",Main.Gadfly.Geom.path)
    else
        Main.Gadfly.plot(dat,x="x",y="y",color="Function",Main.Gadfly.Geom.path, Main.Gadfly.Scale.y_continuous(minvalue=axis[1],maxvalue=axis[2]))      
    end
end




function gadflycontour(x::Vector,y::Vector,z::Matrix;levels=-1,axis=-1)
    require("Gadfly")
    if axis==-1
        axis=[x[1],x[end],y[1],y[end]]
    end
    
    if levels==-1
        Main.Gadfly.plot(x=x,y=y,z=z,Main.Gadfly.Geom.contour,Main.Gadfly.Scale.x_continuous(minvalue=axis[1],maxvalue=axis[2]),Main.Gadfly.Scale.y_continuous(minvalue=axis[3],maxvalue=axis[4]))
    else
        Main.Gadfly.plot(x=x,y=y,z=z,Main.Gadfly.Geom.contour(levels=levels),Main.Gadfly.Scale.x_continuous(minvalue=axis[1],maxvalue=axis[2]),Main.Gadfly.Scale.y_continuous(minvalue=axis[3],maxvalue=axis[4]))
    end
end

function gadflycontourlayer(x::Vector,y::Vector,z::Matrix;levels=-1)
    require("Gadfly")

    
    if levels==-1
        Main.Gadfly.layer(x=x,y=y,z=z,Main.Gadfly.Geom.contour)
    else
        Main.Gadfly.layer(x=x,y=y,z=z,Main.Gadfly.Geom.contour(levels=levels))
    end
end

function dotplot{T<:Real,V<:Real}(x::Vector{T},y::Vector{V};axis=-1)
    require("Gadfly")
    if axis==-1
        Main.Gadfly.plot(x=x,y=y)
    else
        Main.Gadfly.plot(x=x,y=y,Main.Gadfly.Scale.y_continuous(minvalue=axis[1],maxvalue=axis[2]))
    end
end
dotplot{T<:Complex}(x::Vector{T})=dotplot(real(x),imag(x))
dotlayer{T<:Real,V<:Real}(x::Vector{T},y::Vector{V})=Main.Gadfly.layer(x=x,y=y,Main.Gadfly.Geom.point)
dotlayer{T<:Complex}(x::Vector{T})=dotlayer(real(x),imag(x))


function gadflyplot(opts...;kwds...)
    require("Gadfly")
    Main.Gadfly.plot(opts...;kwds...)
end

## domainplot


