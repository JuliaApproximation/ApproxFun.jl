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



## When the basis is real, we can automatically define
# these operators


for ST in (:RealOperator,:ImagOperator)
    @eval begin
        $ST()=$ST(UnsetSpace())
        domainspace(s::$ST)=s.space
        rangespace{S<:RealSpace,T,D}(s::$ST{ReImSpace{S,T,D}})=s.space
        bandinds{S<:RealSpace,T,D}(A::$ST{ReImSpace{S,T,D}})=0,0
        domain(O::$ST)=domain(O.space)
    end
end



function addentries!{S<:RealSpace,T,D}(::RealOperator{ReImSpace{S,T,D}},A,kr::Range)
    for k=kr
        if isodd(k)
            A[k,k]+=1
        end
    end
    A
end

function addentries!{S<:RealSpace,T,D}(::ImagOperator{ReImSpace{S,T,D}},A,kr::Range)
    for k=kr
        if iseven(k)
            A[k,k]+=1
        end
    end
    A
end



# Converts an operator to one that applies on the real and imaginary parts
immutable ReImOperator{O} <: BandedOperator{Float64}
    op::O
end

bandinds(RI::ReImOperator)=2bandinds(RI.op,1),2bandinds(RI.op,2)

# function addentries!(RI::ReImOperator,A,kr::UnitRange)
#     @assert isodd(kr[1])
#     @assert iseven(kr[end])
#     addentries!(RI.op,IndexReIm(A),div(kr[1],2)+1:div(kr[end],2))
#     A
# end


function addentries!(RI::ReImOperator,A,kr::UnitRange)
    divr=(iseven(kr[1])?div(kr[1],2):div(kr[1],2)+1):(iseven(kr[end])?div(kr[end],2):div(kr[end],2)+1)
    B=subview(RI.op,divr,:)
    for k=kr,j=columnrange(RI,k)
        if isodd(k) && isodd(j)
            A[k,j]+=real(B[div(k,2)+1,div(j,2)+1])
        elseif isodd(k) && iseven(j)
            A[k,j]+=-imag(B[div(k,2)+1,div(j,2)])        
        elseif iseven(k) && isodd(j)
            A[k,j]+=imag(B[div(k,2),div(j,2)+1])                    
        else #both iseven 
            A[k,j]+=real(B[div(k,2),div(j,2)])        
        end
    end
    A
end


