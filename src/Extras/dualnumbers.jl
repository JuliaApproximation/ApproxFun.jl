# We need to implement some functionality for the ApproxFun constructor to work
real{T}(::Type{Dual{T}}) = Dual{ApproxFun.real(T)}

# Dual number support. Should there be realpart and dualpart of Space and Domain?
DualNumbers.realpart{S,T<:Dual}(f::Fun{S,T}) = Fun(space(f),realpart(coefficients(f)))
DualNumbers.dualpart{S,T<:Dual}(f::Fun{S,T}) = Fun(space(f),dualpart(coefficients(f)))


DualNumbers.realpart{DD<:Dual}(d::Segment{DD}) = Segment(realpart(d.a),realpart(d.b))
Base.in{DD<:Dual}(x::Number,d::Segment{DD}) = in(x,realpart(d))
Base.in{DD<:Dual}(x::Dual,d::Segment{DD}) = in(realpart(x),d)


valsdomain_type_promote{T<:Real,V<:Real}(::Type{Dual{T}},::Type{V}) =
    Dual{promote_type(T,V)},promote_type(T,V)
valsdomain_type_promote{T<:Complex,V<:Real}(::Type{Dual{T}},::Type{V}) =
    Dual{promote_type(T,V)},promote_type(real(T),V)
valsdomain_type_promote{T<:Real,V<:Real}(::Type{Dual{T}},::Type{Complex{V}}) =
    Dual{promote_type(T,V)},Complex{promote_type(T,V)}
valsdomain_type_promote{T<:Complex,V<:Real}(::Type{Dual{T}},::Type{Complex{V}}) =
    Dual{promote_type(T,Complex{V})},Complex{promote_type(real(T),V)}


for OP in (:plan_chebyshevtransform,:plan_ichebyshevtransform)
    @eval $OP{D<:Dual}(v::Vector{D}) = $OP(@compat(realpart.(v)))
end

if VERSION < v"0.5"
    chebyshevtransform{D<:Dual}(v::Vector{D},plan...) =
        dual(chebyshevtransform(realpart(v),plan...),chebyshevtransform(dualpart(v),plan...))
else
    chebyshevtransform{D<:Dual}(v::Vector{D},plan...) =
        dual.(chebyshevtransform(realpart.(v),plan...),chebyshevtransform(dualpart.(v),plan...))
end

#TODO: Hardy{false}
for OP in (:plan_transform,:plan_itransform)
    for TYP in  (:Fourier,:Laurent,:SinSpace)
        @eval $OP{T<:Dual,D<:Domain}(S::$TYP{D},x::Vector{T}) = $OP(S,@compat(realpart.(x)))
    end
end

for OP in (:transform,:itransform)
    for TYP in (:Fourier,:Laurent,:SinSpace)
        @eval $OP{T<:Dual,D<:Domain}(S::$TYP{D},x::Vector{T},plan) =
            dual($OP(S,@compat(realpart.(x)),plan),$OP(S,@compat(dualpart.(x)),plan))
    end
end

chop!(f::Fun,d::Dual)=chop!(f,realpart(d))



function simplifycfs!{DD<:Dual}(cfs::Vector{DD},tol::Float64=4E-16)
    for k=length(cfs):-2:2
        if maxabs(@compat(realpart.(cfs[k-1:k]))) > maxabs(@compat(dualpart.(cfs[k-1:k])))*tol
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

    if !isa(S,ArraySpace) && isa(f0,Array)
        return dualcfsFun(f,ArraySpace(S,size(f0)...))
    end

    tol =T==Any?100eps():100eps(T)


    fr=typeof(f0)[f(x) for x=r]

    for logn = 4:20
        n=2^logn

        cf=dualFun(f,S,n)

        if maxabs(realpart(cf.coefficients[end-8:end]))<maxabs(dualpart(cf.coefficients[end-8:end]))*tol &&
                                all(k->norm(cf(r[k])-fr[k],1)<1E-4,1:length(r))
            return Fun(S,realpart(simplifycfs!(cf.coefficients,tol*length(cf))))
        end
    end
    warn("Maximum length "*string(2^20+1)*" reached")

    Fun(f,S,2^21)
end
