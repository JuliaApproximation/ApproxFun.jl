## The following is an abandoned attempt to use Compose
#  Redoing everything in Gadfly


function boundingbox(vals::Vector)
    minr=minimum(real(vals));maxr=maximum(real(vals))
    mini=minimum(imag(vals));maxi=maximum(imag(vals))
    UnitBox(minr,
    maxi,
    maxr-minr,mini-maxi)
end
boundingbox(d::Interval)=boundingbox([first(d),last(d)])
boundingbox(d::Circle)=Main.Compose.UnitBox(real(d.center)-d.radius,imag(d.center)+d.radius,2d.radius,-2d.radius)




function arrowhead(d::Interval)
    compose=Main.Compose.compose
    context=Main.Compose.context
    stroke=Main.Compose.stroke    

    arg1=angle(exp(-im*π*0.92)*d)
    arg2=angle(exp(im*π*0.92)*d)    
compose(context(),Main.Compose.polygon([(real(last(d))+0.1length(d)*cos(arg1),imag(last(d))+0.1length(d)*sin(arg1)),
        (real(last(d)), imag(last(d))),
                           (real(last(d))+0.1length(d)*cos(arg2),imag(last(d))+0.1length(d)*sin(arg2))]),
                           stroke("blue"),fill("blue"))    
end

arrowhead(d::Circle)=arrowhead(Interval(d.center-d.radius+im*d.radius,d.center-d.radius))

line(d::Interval)=Main.Compose.line([(real(first(d)), imag(first(d))), (real(last(d)), imag(last(d))) ])


function arrow(d::Interval)
    require("Compose")
    compose=Main.Compose.compose
    context=Main.Compose.context
            
    compose(context(units=boundingbox(d)),line(d),arrowhead(d))
end



circle(d::Circle)=Main.Compose.circle(real(d.center),imag(d.center),d.radius)

function arrow(d::Circle)
    require("Compose")
    compose=Main.Compose.compose
    context=Main.Compose.context
            
    compose(context(units=boundingbox(d)),circle(d),arrowhead(d))
end


function domainplot(d::Domain)
    require("Compose")
    compose=Main.Compose.compose
    context=Main.Compose.context
    stroke=Main.Compose.stroke
    linewidth=Main.Compose.linewidth
        
    compose(
    arrow(d),linewidth(1.))
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