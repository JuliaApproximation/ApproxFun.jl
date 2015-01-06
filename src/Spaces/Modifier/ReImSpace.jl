immutable ReImSpace{S,T,D}<: FunctionSpace{T,D}
    space::S 
end
ReImSpace{T,D}(sp::FunctionSpace{T,D})=ReImSpace{typeof(sp),T,D}(sp)

domain(sp::ReImSpace)=domain(sp.space)

function spaceconversion(f::Vector,a::ReImSpace,b::ReImSpace)
    @assert a.space==b.space
    f 
end

function spaceconversion(f::Vector,a::FunctionSpace,b::ReImSpace)
    if a!=b.space
        cfs=spaceconversion(f,a,b.space)
    end
    ret=Array(Float64,2length(f))
    ret[1:2:end]=real(f)
    ret[2:2:end]=imag(f)    
    ret
end

function spaceconversion(f::Vector,a::ReImSpace,b::FunctionSpace)
    n=length(f)
    if iseven(n)
        ret=f[1:2:end]+1im*f[2:2:end]
    else #odd, so real has one more
        ret=[f[1:2:end-2]+1im*f[2:2:end],f[end]]
    end
    
    if a.space==b
        ret
    else
        spaceconversion(ret,a.space,b)
    end
end

transform(S::ReImSpace,vals::Vector)=spaceconversion(transform(S.space,vals),S.space,S)
evaluate{S<:ReImSpace}(f::Fun{S},x)=evaluate(Fun(f,space(f).space),x)



## Operators

immutable RealOperator{S} <: BandedOperator{Float64}
    space::S
end

immutable ImagOperator{S} <: BandedOperator{Float64}
    space::S
end


for ST in (:RealOperator,:ImagOperator)
    @eval begin
        domainspace(s::$ST)=s.space
        rangespace{S<:FunctionSpace{Float64},T}(s::$ST{ReImSpace{S,T}})=s.space
        bandinds{S<:FunctionSpace{Float64},T}(A::$ST{ReImSpace{S,T}})=0,0
        domain(O::$ST)=domain(O.space)
    end
end



function addentries!{S<:FunctionSpace{Float64},T}(::RealOperator{ReImSpace{S,T}},A,kr::Range)
    for k=kr
        if isodd(k)
            A[k,k]+=1
        end
    end
    A
end

function addentries!{S<:FunctionSpace{Float64},T}(::ImagOperator{ReImSpace{S,T}},A,kr::Range)
    for k=kr
        if iseven(k)
            A[k,k]+=1
        end
    end
    A
end


