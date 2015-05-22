if isdir(Pkg.dir("FastGaussQuadrature"))
    require("FastGaussQuadrature")
    gausshermite(n)=Main.FastGaussQuadrature.gausshermite(n)
end


points(H::Hermite,n)=gausshermite(n)[1]

function transform(H::Hermite,v::Vector,xw::@compat(Tuple{Vector,Vector}))
    x,w=xw
    V=hermitep(0:length(v)-1,x)'
    nrm=(V.^2)*w

    V*(w.*v)./nrm
end
transform(H::Hermite,v::Vector)=transform(H,v,gausshermite(length(v)))

itransform(H::Hermite,cfs::Vector,x::Vector)=hermitep(0:length(cfs)-1,x)*cfs
itransform(H::Hermite,cfs::Vector)=itransform(H,cfs,points(H,length(cfs)))

evaluate{H<:Hermite}(f::Fun{H},x::Number)=dot(hermitep(0:length(f)-1,x),f.coefficients)
evaluate{H<:Hermite}(f::Fun{H},x::Vector)=hermitep(0:length(f)-1,x)*f.coefficients
