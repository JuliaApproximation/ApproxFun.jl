## Plotting


using Gadfly    
import Gadfly.plot
export plot
export complexplot,contour

function gplot(xx,yy;axis=-1)
    if axis==-1
        Gadfly.plot(x=xx, y=yy, Geom.line)
    else
        Gadfly.plot(x=xx, y=yy, Geom.line, Scale.y_continuous(minvalue=axis[1],maxvalue=axis[2]))
    end    
end


function Gadfly.plot(f::IFun;opts...)
    f=pad(f,3length(f)+50)
    gplot(points(f),values(f);opts...)
end