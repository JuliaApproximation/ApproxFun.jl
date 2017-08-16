# We need to implement some functionality for the ApproxFun constructor to work
real(::Type{Dual{T}}) where {T} = Dual{ApproxFun.real(T)}

# Dual number support. Should there be realpart and dualpart of Space and Domain?
DualNumbers.realpart(f::Fun{S,T}) where {S,T<:Dual} = Fun(space(f),realpart(coefficients(f)))
DualNumbers.dualpart(f::Fun{S,T}) where {S,T<:Dual} = Fun(space(f),dualpart(coefficients(f)))


DualNumbers.realpart(d::Segment{DD}) where {DD<:Dual} = Segment(realpart(d.a),realpart(d.b))
Base.in(x::Number,d::Segment{DD}) where {DD<:Dual} = in(x,realpart(d))
Base.in(x::Dual,d::Segment{DD}) where {DD<:Dual} = in(realpart(x),d)


# for QR Factorization.  These have been submitted to DualNumbers
# but we need these duplicates so that they call ApproxFun.flipsign
flipsign(x::Dual,y::Dual) = y == 0 ? flipsign(x, epsilon(y)) : flipsign(x, realpart(y))
flipsign(x, y::Dual) = y == 0 ? flipsign(x, epsilon(y)) : flipsign(x, realpart(y))
flipsign(x::Dual, y::Complex) = y==0 ? x : x*sign(y)
flipsign(x::Dual, y) = dual(flipsign(realpart(x), y), flipsign(epsilon(x), y))



valsdomain_type_promote(::Type{Dual{T}},::Type{V}) where {T<:Real,V<:Real} =
    Dual{promote_type(T,V)},promote_type(T,V)
valsdomain_type_promote(::Type{Dual{T}},::Type{V}) where {T<:Complex,V<:Real} =
    Dual{promote_type(T,V)},promote_type(real(T),V)
valsdomain_type_promote(::Type{Dual{T}},::Type{Complex{V}}) where {T<:Real,V<:Real} =
    Dual{promote_type(T,V)},Complex{promote_type(T,V)}
valsdomain_type_promote(::Type{Dual{T}},::Type{Complex{V}}) where {T<:Complex,V<:Real} =
    Dual{promote_type(T,Complex{V})},Complex{promote_type(real(T),V)}


plan_chebyshevtransform!(x::Vector{T};kind::Integer=1) where {T<:Dual} =
    error("In-place variant not implemented for Dual")

plan_ichebyshevtransform!(x::Vector{T};kind::Integer=1) where {T<:Dual} =
    error("In-place variant not implemented for Dual")


function plan_chebyshevtransform(v::Vector{D};kind::Integer=1) where D<:Dual
    plan = plan_chebyshevtransform(realpart.(v);kind=kind)
    ChebyshevTransformPlan{D,kind,false,typeof(plan)}(plan)
end

function plan_ichebyshevtransform(v::Vector{D};kind::Integer=1) where D<:Dual
    plan = plan_ichebyshevtransform(realpart.(v);kind=kind)
    IChebyshevTransformPlan{D,kind,false,typeof(plan)}(plan)
end




*(P::ChebyshevTransformPlan{D,k,false},v::Vector{D}) where {k,D<:Dual} =
    dual.(P.plan*realpart.(v),P.plan*dualpart.(v))

#TODO: Hardy{false}
for (OP,TransPlan) in ((:plan_transform,:TransformPlan),(:plan_itransform,:ITransformPlan)),
        TYP in  (:Fourier,:Laurent,:SinSpace)
    @eval begin
        function $OP(sp::$TYP{D},x::Vector{T}) where {T<:Dual,D<:Domain}
            plan = $OP(sp,realpart.(x))
            $TransPlan{T,typeof(sp),false,typeof(plan)}(sp,plan)
        end
        *(P::$TransPlan{T,$TYP{D},false},x::Vector{T}) where {T<:Dual,D<:Domain} =
            dual(P.plan*realpart.(x),P.plan*dualpart.(x))
    end
end

chop!(f::Fun,d::Dual)=chop!(f,realpart(d))



function simplifycfs!(cfs::Vector{DD},tol::Float64=4E-16) where DD<:Dual
    for k=length(cfs):-2:2
        if maximum(abs,realpart.(cfs[k-1:k])) > maximum(abs,dualpart.(cfs[k-1:k]))*tol
            return resize!(cfs,k)
        end
    end
    resize!(cfs,3)
end


function dualFun(f,S,n)
    pts=points(S,n) + Dual{Float64}[dual(0.,rand(Bool)) for k=1:n]
    Fun(transform(S,map(f,pts)),S)
end

function dualcfsFun(f,S)
    T = eltype(domain(S))
    if T <: Complex
        T = T.parameters[1] #get underlying real representation
    end
    r=checkpoints(S)
    f0=f(first(r))

    if isa(f0,AbstractArray) && size(S) ≠ size(f0)
        return dualcfsFun(f,Space(fill(S,size(f0))))
    end

    tol =T==Any?100eps():100eps(T)


    fr=typeof(f0)[f(x) for x=r]

    for logn = 4:20
        n=2^logn

        cf=dualFun(f,S,n)

        if maximum(abs,realpart(cf.coefficients[end-8:end]))<maximum(abs,dualpart(cf.coefficients[end-8:end]))*tol &&
                                all(k->norm(cf(r[k])-fr[k],1)<1E-4,1:length(r))
            return Fun(S,realpart(simplifycfs!(cf.coefficients,tol*length(cf))))
        end
    end
    warn("Maximum length "*string(2^20+1)*" reached")

    Fun(f,S,2^21)
end
