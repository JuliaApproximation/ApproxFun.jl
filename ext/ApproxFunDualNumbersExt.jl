module ApproxFunDualNumbersExt

using DualNumbers
using ApproxFun
using ApproxFun: TransformPlan, ITransformPlan
import ApproxFunBase: valsdomain_type_promote
using DomainSets
import FastTransforms: ChebyshevTransformPlan, plan_chebyshevtransform,
                        plan_chebyshevtransform!, plan_ichebyshevtransform,
                        plan_ichebyshevtransform!

# Dual number support. Should there be realpart and dualpart of Space and Domain?
DualNumbers.realpart(f::Fun{S,T}) where {S,T<:Dual} = Fun(space(f),realpart.(coefficients(f)))
DualNumbers.dualpart(f::Fun{S,T}) where {S,T<:Dual} = Fun(space(f),dualpart.(coefficients(f)))


DualNumbers.realpart(d::Segment{<:Dual}) = Segment(realpart(leftendpoint(d)),realpart(rightendpoint(d)))

DomainSets.indomain(x::Number, d::Segment{<:Dual}) = in(x,realpart(d))
DomainSets.indomain(x::Dual, d::Segment{<:Dual}) = in(realpart(x),realpart(d))
Base.isempty(d::Segment{<:Dual}) = isempty(realpart(d))



valsdomain_type_promote(::Type{Dual{T}},::Type{V}) where {T<:Real,V<:Real} =
    Dual{promote_type(T,V)},promote_type(T,V)
valsdomain_type_promote(::Type{Dual{T}},::Type{V}) where {T<:Complex,V<:Real} =
    Dual{promote_type(T,V)},promote_type(real(T),V)
valsdomain_type_promote(::Type{Dual{T}},::Type{Complex{V}}) where {T<:Real,V<:Real} =
    Dual{promote_type(T,V)},Complex{promote_type(T,V)}
valsdomain_type_promote(::Type{Dual{T}},::Type{Complex{V}}) where {T<:Complex,V<:Real} =
    Dual{promote_type(T,Complex{V})},Complex{promote_type(real(T),V)}


plan_chebyshevtransform!(x::AbstractVector{T}, ::Val) where {T<:Dual} =
    error("In-place variant not implemented for Dual")

plan_ichebyshevtransform!(x::AbstractVector{T}, ::Val) where {T<:Dual} =
    error("In-place variant not implemented for Dual")


plan_chebyshevtransform(v::AbstractVector{D}, ::Val{kind}) where {D<:Dual,kind} = plan_chebyshevtransform(realpart.(v), Val(kind))
plan_ichebyshevtransform(v::AbstractVector{D}, ::Val{kind}) where {D<:Dual,kind} = plan_ichebyshevtransform(realpart.(v), Val(kind))



Base.:(*)(P::ChebyshevTransformPlan,v::AbstractVector{<:Dual}) = dual.(P*realpart.(v),P*dualpart.(v))

#TODO: Hardy{false}
for (OP,TransPlan) in ((:plan_transform,:TransformPlan),(:plan_itransform,:ITransformPlan)),
        TYP in  (:Fourier,:Laurent,:SinSpace)
    @eval begin
        function ApproxFun.$OP(sp::$TYP{D},x::AbstractVector{T}) where {T<:Dual,D<:Domain}
            plan = $OP(sp,realpart.(x))
            $TransPlan{T,typeof(sp),false,typeof(plan)}(sp,plan)
        end
        Base.:(*)(P::$TransPlan{T,$TYP{D},false},x::AbstractVector{T}) where {T<:Dual,D<:Domain} =
            dual(P.plan*realpart.(x),P.plan*dualpart.(x))
    end
end

ApproxFun.chop!(f::Fun,d::Dual)=chop!(f,realpart(d))



function simplifycfs!(cfs::AbstractVector{DD},tol::Float64=4E-16) where DD<:Dual
    for k=length(cfs):-2:2
        if maximum(abs,realpart.(cfs[k-1:k])) > maximum(abs,dualpart.(cfs[k-1:k]))*tol
            return resize!(cfs,k)
        end
    end
    resize!(cfs,3)
end


function ApproxFun.dualFun(f,S,n)
    pts=points(S,n) + Dual{Float64}[dual(0.,rand(Bool)) for k=1:n]
    Fun(transform(S,map(f,pts)),S)
end

function ApproxFun.dualcfsFun(f,S)
    T = real(float(eltype(domain(S))))
    r=checkpoints(S)
    f0=f(first(r))

    if isa(f0,AbstractArray) && size(S) ≠ size(f0)
        return ApproxFun.dualcfsFun(f,Space(fill(S,size(f0))))
    end

    tol = 100eps(T)


    fr=typeof(f0)[f(x) for x=r]

    for logn = 4:20
        n=2^logn

        cf=ApproxFun.dualFun(f,S,n)

        if maximum(abs,realpart(cf.coefficients[end-8:end]))<maximum(abs,dualpart(cf.coefficients[end-8:end]))*tol &&
                                all(k->norm(cf(r[k])-fr[k],1)<1E-4,1:length(r))
            return Fun(S,realpart(simplifycfs!(cf.coefficients,tol*length(cf))))
        end
    end
    @warn "Maximum length "*string(2^20+1)*" reached"

    Fun(f,S,2^21)
end

end
