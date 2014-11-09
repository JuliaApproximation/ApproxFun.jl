## Plotting




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

function gadflyplot{T<:Complex}(xx::Vector{T},yy::Vector;opts...)
    warn("Complex domain, projecting to real axis")
    gadflyplot(real(xx),yy;opts...)
end

function gadflyplot{T<:Real}(x::Vector{T},y::Vector{Complex{Float64}};axis=-1)
    require("Gadfly")
    require("DataFrames")    
    r=real(y)
    i=imag(y)
    if axis==-1
        Main.Gadfly.plot(Main.DataFrames.DataFrame(Any[[x,x],[r,i],[fill("Re",length(x)),fill("Im",length(x))]],Main.DataFrames.Index(
        @compat Dict(:x=>1,:y=>2,:Function=>3),
        [:x,:y,:Function])),
        x="x",y="y",color="Function",Main.Gadfly.Geom.path)
    else
        Main.Gadfly.plot(Main.DataFrames.DataFrame(Any[[x,x],[r,i],[fill("Re",length(x)),fill("Im",length(x))]],Main.DataFrames.Index(
        @compat Dict(:x=>1,:y=>2,:Function=>3),
        [:x,:y,:Function])),
        x="x",y="y",color="Function",Main.Gadfly.Geom.path, Main.Gadfly.Scale.y_continuous(minvalue=axis[1],maxvalue=axis[2]))      
    end
end

function gadflycontour(x::Vector,y::Vector,z::Matrix;levels=-1)
    require("Gadfly")
    if levels==-1
        Main.Gadfly.plot(x=x,y=y,z=z,Main.Gadfly.Geom.contour)
    else
        Main.Gadfly.plot(x=x,y=y,z=z,Main.Gadfly.Geom.contour(levels=levels))
    end
end

