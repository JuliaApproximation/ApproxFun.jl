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




## domainplot


function arrow(d::Interval)
    require("Compose")
    line=Main.Compose.line
    compose=Main.Compose.compose
    context=Main.Compose.context
    polygon=Main.Compose.polygon
            
    arg1=angle(exp(-im*π*0.92)*d)
    arg2=angle(exp(im*π*0.92)*d)    
    compose(context(),line([(real(first(d)), imag(first(d))), (real(last(d)), imag(last(d))) ]),
    polygon([(real(last(d))+0.1length(d)*cos(arg1),imag(last(d))+0.1length(d)*sin(arg1)),
        (real(last(d)), imag(last(d))),
                           (real(last(d))+0.1length(d)*cos(arg2),imag(last(d))+0.1length(d)*sin(arg2))]))
end



function domainplot(d::Interval)
    require("Compose")
    compose=Main.Compose.compose
    context=Main.Compose.context
    UnitBox=Main.Compose.UnitBox
    stroke=Main.Compose.stroke
    linewidth=Main.Compose.linewidth
        
    compose(context(units=UnitBox(min(real(first(d)),real(last(d)))-0.5,
    max(imag(first(d)),imag(last(d)))+0.5,
                    2length(d),-2length(d))),
    arrow(d), stroke("blue"),fill("blue"),linewidth(1.))
end

function domainplot{D<:Interval}(d::Vector{D})
    require("Compose")
    compose=Main.Compose.compose
    context=Main.Compose.context
    UnitBox=Main.Compose.UnitBox
    stroke=Main.Compose.stroke
    linewidth=Main.Compose.linewidth
    fill=Main.Compose.fill    
    
    
    C=context(units=UnitBox(
    mapreduce(dk->min(real(first(dk)),real(last(dk))),min,d)-0.5,
    mapreduce(dk->max(imag(first(dk)),imag(last(dk))),max,d)+0.5,
                    4,-4))
    
    compose(C,map(arrow,d)...,stroke("blue"),fill("blue"),linewidth(1.))
end
domainplot(d::UnionDomain)=domainplot(d.domains)