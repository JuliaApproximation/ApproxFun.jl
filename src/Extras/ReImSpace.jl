# Odd rows/columns become real part and even become imag part
# S Should be a Real-valued matrix
immutable IndexReIm{S}
    matrix::S
end


getindex(S::IndexReIm,k,j)=S.matrix[2k-1,2j-1]+im*S.matrix[2k,2j]
function setindex!(S::IndexReIm,x,k,j)
    S.matrix[2k-1,2j-1]=real(x)
    S.matrix[2k,2j]=imag(x)
    x
end





## Splits real and imaginary parts


for TYP in (:ReSpace,:ImSpace,:ReImSpace)
    @eval begin
        immutable $TYP{S,T,d}<: Space{T,AnyDomain,d}
            space::S
        end

        $TYP{T,D,d}(sp::Space{T,D,d})=$TYP{typeof(sp),T,d}(sp)

        domain(sp::$TYP)=domain(sp.space)
        spacescompatible(a::$TYP,b::$TYP)=spacescompatible(a.space,b.space)

        function coefficients(f::Vector,a::$TYP,b::$TYP)
            @assert spacescompatible(a.space,b.space)
            f
        end

        transform(S::$TYP,vals::Vector)=coefficients(transform(S.space,vals),S.space,S)
        evaluate(f::AbstractVector,S::$TYP,x)=evaluate(f,space(f).space,x)

        canonicalspace(a::$TYP)=$TYP(canonicalspace(a.space))
    end

    for OP in (:maxspace,:conversion_type)
        @eval $OP(a::$TYP,b::$TYP)=$TYP($OP(a.space,b.space))
    end
end

## compat

## Resolve conflict
for TYP in (:ReImSpace,:ReSpace,:ImSpace)
    for V in (:SliceSpace,:SumSpace,:PiecewiseSpace)
        @eval coefficients(::Vector,sp::$TYP,slp::$V)=error("coefficients not implemented from "*string(typeof(sp))*" to "*string(typeof(slp)))
    end
    @eval coefficients(::Vector,sp::SliceSpace,slp::$TYP)=error("coefficients not implemented from "*typeof(sp)*" to "*typeof(slp))
end


coefficients(f::Vector,a::ImSpace,b::ReSpace)=zeros(f)
coefficients(f::Vector,a::ReSpace,b::ImSpace)=zeros(f)
coefficients(f::Vector,a::ReSpace,b::ReImSpace)=coefficients(f,a,a.space,b)
coefficients(f::Vector,a::ImSpace,b::ReImSpace)=coefficients(f,a,a.space,b)
coefficients(f::Vector,a::ReImSpace,b::ReSpace)=coefficients(f,a,a.space,b)
coefficients(f::Vector,a::ReImSpace,b::ImSpace)=coefficients(f,a,a.space,b)
coefficients(f::Vector,a::Space,b::ReSpace)=real(coefficients(f,a,b.space))
coefficients(f::Vector,a::Space,b::ImSpace)=imag(coefficients(f,a,b.space))
coefficients(f::Vector,a::ReSpace,b::Space)=(@assert isa(eltype(f),Real);coefficients(f,a.space,b))
coefficients(f::Vector,a::ImSpace,b::Space)=(@assert isa(eltype(f),Real);coefficients(1im*f,a.space,b))

function coefficients(f::Vector,a::Space,b::ReImSpace)
    if a!=b.space
        f=coefficients(f,a,b.space)
    end
    ret=Array(Float64,2length(f))
    ret[1:2:end]=real(f)
    ret[2:2:end]=imag(f)
    ret
end


function coefficients(f::Vector,a::ReImSpace,b::Space)
    n=length(f)
    if iseven(n)
        ret=f[1:2:end]+1im*f[2:2:end]
    else #odd, so real has one more
        ret=[f[1:2:end-2]+1im*f[2:2:end],f[end]]
    end

    if a.space==b
        ret
    else
        coefficients(ret,a.space,b)
    end
end


#union_rule(a::Space,b::ReImSpace)=union(a,b.space)

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
        rangespace{S<:RealSpace,T}(s::$ST{ReImSpace{S,T}})=s.space
        bandinds{S<:RealSpace,T}(A::$ST{ReImSpace{S,T}})=0,0
        domain(O::$ST)=domain(O.space)
        choosedomainspace(s::$ST{UnsetSpace},sp)=ReImSpace(sp)
    end
end

getindex{S<:RealSpace,T}(::RealOperator{ReImSpace{S,T}},k::Integer,j::Integer) =
    ifelse(isodd(k)&&j==k,one(T),zero(T))

getindex{S<:RealSpace,T}(::ImagOperator{ReImSpace{S,T}},k::Integer,j::Integer) =
    ifelse(iseven(k)&&j==k,one(T),zero(T))




# Converts an operator to one that applies on the real and imaginary parts
immutable ReImOperator{O,T} <: BandedOperator{T}
    op::O
end

ReImOperator(op)=ReImOperator{typeof(op),Float64}(op)
Base.convert{T}(::Type{BandedOperator{T}},R::ReImOperator)=ReImOperator{typeof(R.op),T}(R.op)

bandinds(RI::ReImOperator)=2bandinds(RI.op,1),2bandinds(RI.op,2)

for OP in (:rangespace,:domainspace)
    @eval $OP(R::ReImOperator)=ReImSpace($OP(R.op))
end


function getindex(RI::ReImOperator,k::Integer,j::Integer)
    if isodd(k) && isodd(j)
        real(B[k÷2+1,j÷2+1])
    elseif isodd(k) && iseven(j)
        -imag(B[k÷2+1,j÷2])
    elseif iseven(k) && isodd(j)
        imag(B[k÷2,j÷2+1])
    else #both iseven
        real(B[k÷2,j÷2])
    end
end


Multiplication{D<:UnivariateSpace,T,S,V}(f::Fun{D,T},sp::ReImSpace{S,V,2})=MultiplicationWrapper(f,ReImOperator(Multiplication(f,sp.space)))
Multiplication{D,T}(f::Fun{D,T},sp::ReImSpace)=MultiplicationWrapper(f,ReImOperator(Multiplication(f,sp.space)))


# ReFunctional/ImFunctional are used
# to take the real/imag part of a functional
for TYP in (:ReFunctional,:ImFunctional)
    @eval begin
        immutable $TYP{O,T} <: Operator{T}
            functional::O
        end

        @functional $TYP

        $TYP{T}(func::Operator{T}) = $TYP{typeof(func),real(T)}(func)

        domainspace(RF::$TYP) = ReImSpace(domainspace(RF.functional))
        rangespace(RF::$TYP) = ConstantSpace()
    end
end

function getindex{R,T}(S::ReFunctional{R,T},kr::Range)
     kr1=div(kr[1]+1,2):div(kr[end]+1,2)
     res=S.functional[kr1]
     T[isodd(k)?real(res[div(k+1,2)-first(kr1)+1]):-imag(res[div(k+1,2)-first(kr1)+1]) for k=kr]
end

function getindex{R,T}(S::ImFunctional{R,T},kr::Range)
     kr1=div(kr[1]+1,2):div(kr[end]+1,2)
     res=S.functional[kr1]
     T[isodd(k)?imag(res[div(k+1,2)-first(kr1)+1]):real(res[div(k+1,2)-first(kr1)+1]) for k=kr]
end


## Definite Integral
# disabled since its complex, which would lead to a complex solution to \
# breaking the point of ReImSpace

# DefiniteIntegral(dsp::ReImSpace) = DefiniteIntegral{typeof(dsp),eltype(DefiniteIntegral(dsp.space))}(dsp)
# datalength{RI<:ReImSpace}(Σ::DefiniteIntegral{RI})=2datalength(DefiniteIntegral(domainspace(Σ).space))
#
# function getindex{RI<:ReImSpace,T}(Σ::DefiniteIntegral{RI,T},kr::Range)
#     kr1=div(kr[1]+1,2):div(kr[end]+1,2)
#     res=DefiniteIntegral(domainspace(Σ).space)[kr1]
#     T[isodd(k)?res[div(k+1,2)-first(kr1)+1]:im*res[div(k+1,2)-first(kr1)+1] for k=kr]
# end





## LaurentOperator

## Real/Imag

for TYP in (:RealOperator,:ImagOperator)
    @eval begin
        rangespace{T,DD}(R::$TYP{ReImSpace{Taylor{DD},T}})=Fourier(domain(R))
    end
end

bandinds{T,DD}(::RealOperator{ReImSpace{Taylor{DD},T}})=0,2
bandinds{T,DD}(::ImagOperator{ReImSpace{Taylor{DD},T}})=0,1


## Re[r z^k] = r cos(k x), Re[im q z^k] = -sin(k x)

function getindex{T,DD}(R::RealOperator{ReImSpace{Taylor{DD},T}},k::Integer,j::Integer)
    if isodd(k)&&j==k
        one(T)
    elseif iseven(k)&&j==k+2
        -one(T)
    else
        zero(T)
    end
end


## Im[r z^k] = r sin(k x), Im[im q z^k] = cos(k x)

getindex{T,DD}(R::ImagOperator{ReImSpace{Taylor{DD},T}},k::Integer,j::Integer)=
    j==k+1?one(T):zero(T)


# Neg

# spaces lose zeroth coefficient
for TYP in (:RealOperator,:ImagOperator)
    @eval begin
        rangespace{T,DD}(R::$TYP{ReImSpace{Hardy{false,DD},T}})=SliceSpace(Fourier(domain(R)),1)
    end
end


bandinds{T,D}(::RealOperator{ReImSpace{Hardy{false,D},T}})=-1,1
bandinds{T,D}(::ImagOperator{ReImSpace{Hardy{false,D},T}})=0,0


## Re[r z^(-k)] = r cos(k x), Re[im q z^(-k)] = -sin(-k x)= sin(k x)
function getindex{T,D}(R::RealOperator{ReImSpace{Hardy{false,D},T}},k::Integer,j::Integer)
    if (isodd(k) && j==k+1) || (iseven(k) && j==k-1)  # imag part
        one(T)
    else       # real part
        zero(T)
    end
end


## Im[r z^(-k)] = r sin(-k x)=-r sin(kx), Im[im q z^(-k)] = cos(-k x)=cos(kx)
function getindex{T,D}(R::ImagOperator{ReImSpace{Hardy{false,D},T}},k::Integer,j::Integer)
    if k==j
        if isodd(k)
            -one(T)
        else
            one(T)
        end
    else
        zero(T)
    end
end


# spaces lose zeroth coefficient
for TYP in (:RealOperator,:ImagOperator)
    @eval begin
        rangespace{T,DD}(R::$TYP{ReImSpace{Laurent{DD},T}})=Fourier(domain(R))
    end
end

bandinds{T,DD}(::RealOperator{ReImSpace{Laurent{DD},T}})=0,2

function getindex{T,DD}(R::RealOperator{ReImSpace{Laurent{DD},T}},k::Integer,j::Integer)
    if isodd(k) && k==j
        one(T)
    elseif iseven(k) && j==k+2
        iseven(div(k,2))?-one(T):one(T)
    else
        zero(T)
    end
end






## just takes realpart of operator
immutable ReOperator{O,T} <: BandedOperator{T}
    op::O
end

ReOperator(op)=ReOperator{typeof(op),Float64}(op)
Base.convert{BT<:Operator}(::Type{BT},R::ReOperator)=ReOperator{typeof(R.op),eltype(BT)}(R.op)

for OP in (:rangespace,:domainspace,:bandinds)
    @eval $OP(R::ReOperator)=$OP(R.op)
end



getindex(RI::ReOperator,k::Integer,j::Integer) =
    real(RI.op[k,j])

choosedomainspace(R::ReOperator,sp)=choosedomainspace(R.op,sp)
for OP in (:promotedomainspace,:promoterangespace)
    @eval begin
        $OP(R::ReOperator,sp::UnsetSpace)=ReOperator($OP(R.op,sp))
        $OP(R::ReOperator,sp::AnySpace)=ReOperator($OP(R.op,sp))
        $OP(R::ReOperator,sp::Space)=ReOperator($OP(R.op,sp))
    end
end



# TODO: can't do this because UnsetSpace might change type
#Base.real{T<:Real}(op::BandedOperator{T})=op
Base.real(op::BandedOperator)=ReOperator(op)






#####
# ReReOperator takes the real part of two operators
# this allows for well-posed equations
#####


immutable ReReOperator{S,V,T} <: BandedOperator{T}
    ops::Tuple{S,V}
    function ReReOperator(ops)
            #TODO: promotion
        @assert domainspace(ops[1])==domainspace(ops[2])
        @assert rangespace(ops[1])==rangespace(ops[2])
        new(ops)
    end
end


ReReOperator{S,V}(ops::Tuple{S,V})=ReReOperator{S,V,Float64}(ops)
ReReOperator(ops1,ops2)=ReReOperator((ops1,ops2))
Base.real(S::BandedOperator,V::BandedOperator)=ReReOperator(S,V)

bandinds(R::ReReOperator)=min(2bandinds(R.ops[1],1)-1,2bandinds(R.ops[2],1)-2),max(2bandinds(R.ops[1],2)+1,2bandinds(R.ops[2],2))

domainspace(R::ReReOperator)=ReImSpace(domainspace(R.ops[1]))
rangespace(R::ReReOperator)=ArraySpace(rangespace(R.ops[1]),2)
