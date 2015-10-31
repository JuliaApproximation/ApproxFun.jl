## Plotting

export domainplot

@try_import Gadfly
@try_import DataFrames

## Vector routines

function gadflyopts(;axis=-1,title=-1)
    opts=Any[Gadfly.Geom.path]

    if axis != -1
        if length(axis)==2
            opts=[opts;Gadfly.Scale.y_continuous(minvalue=axis[1],maxvalue=axis[2])]
        elseif length(axis)==4
            opts=[opts;Gadfly.Scale.x_continuous(minvalue=axis[1],maxvalue=axis[2]);
                  Gadfly.Scale.y_continuous(minvalue=axis[3],maxvalue=axis[4])]
        end
    end


    if title != -1
        opts=[opts;Gadfly.Guide.title(string(title))]
    end
    opts
end

gadflyopts(opts...)=gadflyopts(;opts...)

function gadflyplot{T<:Real}(xx::Vector{T},yy::Vector,v...;opts...)
    Gadfly.plot(x=xx, y=yy,v...,gadflyopts(opts...)...)
end


function gadflysemilogy{T<:Real}(xx::Vector{T},yy::Vector,v...;opts...)
    Gadfly.plot(x=xx, y=yy,v...,Gadfly.Scale.y_log10,gadflyopts(opts...)...)
end

gadflysemilogy(yy::Vector,v...;opts...)=gadflysemilogy(collect(1:length(yy)),yy,v...;opts...)

function gadflylayer{T<:Real}(xx::Vector{T},yy::Vector)
    Gadfly.layer(x=xx, y=yy, Gadfly.Geom.path)
end


function gadflyplot{T<:Complex}(xx::Vector{T},yy::Vector,v...;opts...)
    warn("Complex domain, projecting to real axis")
    gadflyplot(real(xx),yy,v...;opts...)
end

function gadflyplot{T<:Real}(x::Vector{T},y::Vector{Complex{Float64}},v...;opts...)
    r=real(y)
    i=imag(y)

    dat=DataFrames.DataFrame(x=[x;x],
                                  y=[r;i],
                                  Function=[fill("Re",length(x));fill("Im",length(x))])

    Gadfly.plot(dat,x="x",y="y",color="Function",v...,gadflyopts(opts...)...)
end

#Plot multiple contours
# columns are plots
function gadflyplot{T<:Real,V<:Real}(x::Matrix{T},y::Matrix{V},v...;opts...)
    dat=DataFrames.DataFrame(x=vec(x),
                             y=vec(y),
                             Function=reduce(vcat,[fill(string(k),size(x,1)) for k=1:size(y,2)]))

    Gadfly.plot(dat,x="x",y="y",color="Function",v...,gadflyopts(opts...)...)
end

function gadflycontour(x::Vector,y::Vector,z::Matrix,v...;levels=-1,axis=-1)
    if axis==-1
        axis=[x[1],x[end],y[1],y[end]]
    end

    if levels==-1
        Gadfly.plot(x=x,y=y,z=z,Gadfly.Geom.contour,
                    Gadfly.Scale.x_continuous(minvalue=axis[1],maxvalue=axis[2]),
                    Gadfly.Scale.y_continuous(minvalue=axis[3],maxvalue=axis[4]),v...)
    else
        Gadfly.plot(x=x,y=y,z=z,
                    Gadfly.Geom.contour(levels=levels),
                    Gadfly.Scale.x_continuous(minvalue=axis[1],maxvalue=axis[2]),
                    Gadfly.Scale.y_continuous(minvalue=axis[3],maxvalue=axis[4]),v...)
    end
end

function gadflycontourlayer(x::Vector,y::Vector,z::Matrix;levels=-1)
    if levels==-1
        Gadfly.layer(x=x,y=y,z=z,Gadfly.Geom.contour)
    else
        Gadfly.layer(x=x,y=y,z=z,Gadfly.Geom.contour(levels=levels))
    end
end

function dotplot{T<:Real,V<:Real}(x::Vector{T},y::Vector{V},v...;axis=-1)
    if axis==-1
        Gadfly.plot(x=x,y=y,v...)
    else
        Gadfly.plot(x=x,y=y,Gadfly.Scale.y_continuous(minvalue=axis[1],maxvalue=axis[2]),v...)
    end
end
dotplot{T<:Number}(x::Vector{T})=dotplot(real(x),imag(x))
dotlayer{T<:Real,V<:Real}(x::Vector{T},y::Vector{V})=Gadfly.layer(x=x,y=y,Gadfly.Geom.point)
dotlayer{T<:Number}(x::Vector{T})=dotlayer(real(x),imag(x))


function gadflyplot(opts...;kwds...)
    Gadfly.plot(opts...;kwds...)
end


## Functional
gadflydeltaplot(x0::Number,c::Number)=Gadfly.plot(x=ones(2)*x0,y=linspace(0.,c,2),
                                                         Gadfly.Geom.line)

gadflydeltalayer(x0::Number,c::Number)=Gadfly.layer(x=ones(2)*x0,y=linspace(0.,c,2),
                                                         Gadfly.Geom.line)

gadflydeltaplot(x0::Vector,c::Vector)=Gadfly.plot(map(gadflydeltalayer,x0,c)...)


gadflyplot{S}(B::Evaluation{S,Float64})=gadflydeltaplot(1,B.x)
gadflyplot{S}(B::Evaluation{S,Bool})=gadflydeltaplot(1,B.x?first(domain(B)):last(domain(B)))
gadflyplot{T<:Real,E<:Evaluation}(B::ConstantTimesFunctional{T,E})=gadflyplot(B.op)


## hist
gadflyhist(r::Vector;kwds...)=Gadfly.plot(x=r,Gadfly.Geom.histogram(kwds...))
