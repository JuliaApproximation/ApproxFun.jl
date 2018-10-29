points(H::Hermite,n) = gausshermite(n)[1] ./ sqrt.(H.L)

plan_transform(H::Hermite,v::AbstractVector) = TransformPlan(H,gausshermite(length(v)),Val{false})
plan_itransform(H::Hermite,cfs::AbstractVector) = ITransformPlan(H,points(Hermite(),length(cfs)),Val{false})



function *(P::TransformPlan{T,H,false},vals::AbstractVector) where {T,H<:Hermite}
    x,w = P.plan
    V=hermitep.(0:length(vals)-1,transpose(x))
    nrm=(V.^2)*w
    V*(w.*vals)./nrm
end

*(P::ITransformPlan{T,H,false},cfs::AbstractVector) where {T,H<:Hermite} =
    hermitep.((0:length(cfs)-1)',P.plan)*cfs
