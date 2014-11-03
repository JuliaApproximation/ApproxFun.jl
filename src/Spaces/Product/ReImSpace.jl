immutable ReImSpace{S,T}<: DomainSpace{T}
    space::S 
end
ReImSpace{T}(sp::DomainSpace{T})=ReImSpace{typeof(sp),T}(sp)



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
        Complex{Float64}[f[1:2:end-2]+1im*f[2:2:end],f[end]]
    end
end

