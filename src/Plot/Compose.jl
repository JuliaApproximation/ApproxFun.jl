## The following is an abandoned attempt to use Compose
#  Redoing everything in Gadfly

@try_import Compose

function boundingbox(vals::Vector)
    minr=minimum(real(vals));maxr=maximum(real(vals))
    mini=minimum(imag(vals));maxi=maximum(imag(vals))
    Compose.UnitBox(minr,
    maxi,
    maxr-minr,mini-maxi)
end
boundingbox(d::Interval)=boundingbox([first(d),last(d)])
boundingbox(d::Circle)=Compose.UnitBox(real(d.center)-d.radius,imag(d.center)+d.radius,2d.radius,-2d.radius)




function arrowhead(d::Interval)
    arg1=angle(exp(-im*π*0.92)*d)
    arg2=angle(exp(im*π*0.92)*d)
    Compose.compose(Compose.context(),Compose.polygon([(real(last(d))+0.1length(d)*cos(arg1),imag(last(d))+0.1length(d)*sin(arg1)),
        (real(last(d)), imag(last(d))),
                           (real(last(d))+0.1length(d)*cos(arg2),imag(last(d))+0.1length(d)*sin(arg2))]),
                           Compose.stroke("blue"),fill("blue"))
end

arrowhead(d::Circle)=arrowhead(Interval(d.center-d.radius+im*d.radius,d.center-d.radius))

line(d::Interval)=Compose.line([(real(first(d)), imag(first(d))), (real(last(d)), imag(last(d))) ])


function arrow(d::Interval)
    Compose.compose(Compose.context(units=boundingbox(d)),line(d),arrowhead(d))
end



circle(d::Circle)=Compose.circle(real(d.center),imag(d.center),d.radius)

function arrow(d::Circle)
    Compose.compose(Compose.context(units=boundingbox(d)),circle(d),arrowhead(d))
end


function domainplot(d::Domain)
    Compose.compose(arrow(d),Compose.linewidth(1.))
end

function domainplot{D<:Interval}(d::Vector{D})
    C=Compose.context(units=Compose.UnitBox(
    mapreduce(dk->min(real(first(dk)),real(last(dk))),min,d)-0.5,
    mapreduce(dk->max(imag(first(dk)),imag(last(dk))),max,d)+0.5,
                    4,-4))

    Compose.compose(C,map(arrow,d)..., Compose.stroke("blue"),
                    Compose.fill("blue"), Compose.linewidth(1.))
end
domainplot(d::UnionDomain)=domainplot(d.domains)
