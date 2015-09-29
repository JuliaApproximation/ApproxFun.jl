
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
        evaluate{S<:$TYP}(f::Fun{S},x)=evaluate(Fun(f,space(f).space),x)

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



function addentries!{S<:RealSpace,T}(::RealOperator{ReImSpace{S,T}},A,kr::Range,::Colon)
    for k=kr
        if isodd(k)
            A[k,k]+=1
        end
    end
    A
end

function addentries!{S<:RealSpace,T}(::ImagOperator{ReImSpace{S,T}},A,kr::Range,::Colon)
    for k=kr
        if iseven(k)
            A[k,k]+=1
        end
    end
    A
end



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

# function addentries!(RI::ReImOperator,A,kr::UnitRange,::Colon)
#     @assert isodd(kr[1])
#     @assert iseven(kr[end])
#     addentries!(RI.op,IndexReIm(A),div(kr[1],2)+1:div(kr[end],2))
#     A
# end


function addentries!(RI::ReImOperator,A,kr::UnitRange,::Colon)
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


Multiplication{D<:UnivariateSpace,T,S,V}(f::Fun{D,T},sp::ReImSpace{S,V,2})=MultiplicationWrapper(f,ReImOperator(Multiplication(f,sp.space)))
Multiplication{D,T}(f::Fun{D,T},sp::ReImSpace)=MultiplicationWrapper(f,ReImOperator(Multiplication(f,sp.space)))


# ReFunctional/ImFunctional are used
# to take the real/imag part of a functional
for TYP in (:ReFunctional,:ImFunctional)
    @eval begin
        immutable $TYP{O,T} <: Functional{T}
            functional::O
        end

        $TYP{T}(func::Functional{T})=$TYP{typeof(func),real(T)}(func)

        domainspace(RF::$TYP)=ReImSpace(domainspace(RF.functional))
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

Base.real(F::Functional)=ReFunctional(F)
Base.imag(F::Functional)=ImFunctional(F)


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
function addentries!{T,DD}(R::RealOperator{ReImSpace{Taylor{DD},T}},A,kr::Range,::Colon)
    for k=kr
        if isodd(k)         # real part
            A[k,k]+=1
        elseif iseven(k)    # imag part
            A[k,k+2]+=-1
        end
    end
    A
end

## Im[r z^k] = r sin(k x), Im[im q z^k] = cos(k x)
function addentries!{T,DD}(R::ImagOperator{ReImSpace{Taylor{DD},T}},A,kr::Range,::Colon)
    for k=kr
        A[k,k+1]+=1
    end
    A
end

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
function addentries!{T,D}(R::RealOperator{ReImSpace{Hardy{false,D},T}},A,kr::Range,::Colon)
    for k=kr
        if isodd(k)    # imag part
            A[k,k+1]+=1
        elseif iseven(k)         # real part
            A[k,k-1]+=1
        end
    end
    A
end

## Im[r z^(-k)] = r sin(-k x)=-r sin(kx), Im[im q z^(-k)] = cos(-k x)=cos(kx)
function addentries!{T,D}(R::ImagOperator{ReImSpace{Hardy{false,D},T}},A,kr::Range,::Colon)
    for k=kr
        A[k,k]+=isodd(k)?-1:1
    end
    A
end


# spaces lose zeroth coefficient
for TYP in (:RealOperator,:ImagOperator)
    @eval begin
        rangespace{T,DD}(R::$TYP{ReImSpace{Laurent{DD},T}})=Fourier(domain(R))
    end
end

bandinds{T,DD}(::RealOperator{ReImSpace{Laurent{DD},T}})=0,2
function addentries!{T,DD}(R::RealOperator{ReImSpace{Laurent{DD},T}},A,kr::Range,::Colon)
    for k=kr
        if isodd(k)    # real part
            A[k,k]+=1
        elseif iseven(k)         # odd part
            A[k,k+2]+=iseven(div(k,2))?-1:1
        end
    end
    A
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



function addentries!(RI::ReOperator,A,kr::UnitRange,::Colon)
    if isa(eltype(RI.op),Real)
        addentries!(RI.op,A,kr,:)
    else
        B=subview(RI.op,kr,:)
        for k=kr,j=columnrange(RI,k)
           A[k,j]+=real(B[k,j])
        end
        A
    end
end

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
