immutable ReImSpace{S,T}<: DomainSpace{T}
    space::S 
end
ReImSpace{T}(sp::DomainSpace{T})=ReImSpace{typeof(sp),T}(sp)

domain(sp::ReImSpace)=domain(sp.space)

function spaceconversion{S<:DomainSpace}(f::Vector,a::ReImSpace{S},b::ReImSpace{S})
    @assert a.space==b.space
    f 
end

function spaceconversion{S<:DomainSpace}(f::Vector,a::S,b::ReImSpace{S})
    ret=Array(Float64,2length(f))
    ret[1:2:end]=real(f)
    ret[2:2:end]=imag(f)    
    ret
end

function spaceconversion{S<:DomainSpace}(f::Vector,a::ReImSpace{S},b::S)
    n=length(f)
    if iseven(n)
        f[1:2:end]+1im*f[2:2:end]
    else #odd, so real has one more
        [f[1:2:end-2]+1im*f[2:2:end],f[end]]
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
        rangespace{S<:DomainSpace{Float64},T}(s::$ST{ReImSpace{S,T}})=s.space
        bandinds{S<:DomainSpace{Float64},T}(A::$ST{ReImSpace{S,T}})=0,0
        domain(O::$ST)=domain(O.space)
    end
end



function addentries!{S<:DomainSpace{Float64},T}(::RealOperator{ReImSpace{S,T}},A::ShiftArray,kr::Range)
    for k=kr
        if isodd(k)
            A[k,0]+=1
        end
    end
    A
end

function addentries!{S<:DomainSpace{Float64},T}(::ImagOperator{ReImSpace{S,T}},A::ShiftArray,kr::Range)
    for k=kr
        if iseven(k)
            A[k,0]+=1
        end
    end
    A
end


