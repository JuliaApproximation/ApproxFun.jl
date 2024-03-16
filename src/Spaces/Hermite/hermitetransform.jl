points(H::Hermite,n)=gausshermite(n)[1]

plan_transform(H::Hermite,v::Vector) = gausshermite(length(v))
plan_itransform(H::Hermite,cfs::Vector) = points(H,length(cfs))
function transform(H::Hermite,vals,plan::Tuple{Vector,Vector})
    x,w = plan
    V=hermitep(0:length(vals)-1,x)'
    nrm=(V.^2)*w

    V*(w.*vals)./nrm
end
itransform(H::Hermite,cfs,plan::Vector) = hermitep(0:length(cfs)-1,tocanonical(H,plan))*cfs
