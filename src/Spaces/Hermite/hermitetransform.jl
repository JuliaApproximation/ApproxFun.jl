points(H::Hermite,n)=gausshermite(n)[1]

plan_transform(H::Hermite,v::AbstractVector) = TransformPlan(H,gausshermite(length(v)),Val{false})
plan_itransform(H::Hermite,cfs::AbstractVector) = ITransformPlan(H,points(H,length(cfs)),Val{false})



function *{T,HH<:Hermite}(P::TransformPlan{T,HH,false},vals::AbstractVector)
    x,w = P.plan
    V=hermitep(0:length(vals)-1,x)'
    nrm=(V.^2)*w
    V*(w.*vals)./nrm
end

*{T,HH<:Hermite}(P::ITransformPlan{T,HH,false},cfs::AbstractVector) =
    hermitep(0:length(cfs)-1,tocanonical(H,P.plan))*cfs
