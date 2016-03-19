

## Converison

#ensure that COnversion is called
coefficients{DD}(cfs::Vector,A::Fourier{DD},B::Laurent{DD})=Conversion(A,B)*cfs
coefficients{DD}(cfs::Vector,A::Laurent{DD},B::Fourier{DD})=Conversion(A,B)*cfs

hasconversion{DD}(::Fourier{DD},::Laurent{DD})=true
hasconversion{DD}(::Laurent{DD},::Fourier{DD})=true

Conversion{DD}(a::Laurent{DD},b::Fourier{DD})=ConcreteConversion(a,b)
Conversion{DD}(a::Fourier{DD},b::Laurent{DD})=ConcreteConversion(a,b)

function addentries!{DD}(C::ConcreteConversion{Laurent{DD},Fourier{DD}},A,kr::Range,::Colon)
    for k=kr
        if k==1
            A[k,k]+=1.
        elseif iseven(k)
            A[k,k]+=-1.im
            A[k,k+1]+=1.im
        else #isodd(k)
            A[k,k]+=1
            A[k,k-1]+=1
        end
    end
    A
end
function addentries!{DD}(C::ConcreteConversion{Fourier{DD},Laurent{DD}},A,kr::Range,::Colon)
    for k=kr
        if k==1
            A[k,k]+=1.
        elseif iseven(k)
            A[k,k]+=0.5im
            A[k,k+1]+=0.5
        else #isodd(k)
            A[k,k]+=0.5
            A[k,k-1]+=-0.5im
        end
    end
    A
end

bandinds{DD}(::ConcreteConversion{Laurent{DD},Fourier{DD}})=-1,1
bandinds{DD}(::ConcreteConversion{Fourier{DD},Laurent{DD}})=-1,1

for RULE in (:conversion_rule,:maxspace_rule,:union_rule)
    @eval function $RULE{DD}(A::Laurent{DD},B::Fourier{DD})
        @assert domainscompatible(A,B)
        B
    end
end

conversion_type{DD<:Circle}(A::Fourier{DD},B::Fourier{DD})=domain(A).orientation?A:B

hasconversion{DD}(A::Fourier{DD},B::Fourier{DD})=domain(A) == reverse(domain(B))
Conversion{DD}(A::Fourier{DD},B::Fourier{DD})=A==B?ConversionWrapper(eye(A)):ConcreteConversion(A,B)
bandinds{DD}(::ConcreteConversion{Fourier{DD},Fourier{DD}})=0,0

function addentries!{DD}(C::ConcreteConversion{Fourier{DD},Fourier{DD}},A,kr::Range,::Colon)
    @assert domain(domainspace(C)) == reverse(domain(rangespace(C)))
    for k=kr
        A[k,k]+=iseven(k)?(-1):1
    end
    A
end




### Cos/Sine


Derivative(S::Union{CosSpace,SinSpace},order)=ConcreteDerivative(S,order)

bandinds{CS<:CosSpace}(D::ConcreteDerivative{CS})=iseven(D.order)?(0,0):(0,1)
bandinds{S<:SinSpace}(D::ConcreteDerivative{S})=iseven(D.order)?(0,0):(-1,0)
rangespace{S<:CosSpace}(D::ConcreteDerivative{S})=iseven(D.order)?D.space:SinSpace(domain(D))
rangespace{S<:SinSpace}(D::ConcreteDerivative{S})=iseven(D.order)?D.space:CosSpace(domain(D))


function addentries!{CS<:CosSpace}(D::ConcreteDerivative{CS},A,kr::Range,::Colon)
    d=domain(D)
    @assert isa(d,PeriodicInterval)
    m=D.order
    C=2/(d.b-d.a)*π

    for k=kr
        if mod(m,4)==0
            A[k,k] += (C*(k-1))^m
        elseif mod(m,4)==2
            A[k,k] -= (C*(k-1))^m
        elseif mod(m,4)==1
            A[k,k+1] -= (C*k)^m
        elseif mod(m,4)==3
            A[k,k+1] += (C*k)^m
        end
    end

    A
end

function addentries!{CS<:SinSpace}(D::ConcreteDerivative{CS},A,kr::Range,::Colon)
    d=domain(D)
    @assert isa(d,PeriodicInterval)
    m=D.order
    C=2/(d.b-d.a)*π

    for k=kr
        if mod(m,4)==0
            A[k,k] += (C*k)^m
        elseif mod(m,4)==2
            A[k,k] += -(C*k)^m
        elseif k>1 && mod(m,4)==1
            A[k,k-1] += (C*(k-1))^m
        elseif k>1 && mod(m,4)==3
            A[k,k-1] += -(C*(k-1))^m
        end
    end

    A
end

Integral(::CosSpace,m::Integer)=error("Integral not defined for CosSpace.  Use Integral(SliceSpace(CosSpace(),1)) if first coefficient vanishes.")


bandinds{CS<:SinSpace}(D::ConcreteIntegral{CS})=iseven(D.order)?(0,0):(-1,0)
rangespace{S<:CosSpace}(D::ConcreteIntegral{S})=iseven(D.order)?D.space:SinSpace(domain(D))
rangespace{S<:SinSpace}(D::ConcreteIntegral{S})=iseven(D.order)?D.space:CosSpace(domain(D))

function addentries!{CS<:SinSpace}(D::ConcreteIntegral{CS},A,kr::Range,::Colon)
    d=domain(D)
    @assert isa(d,PeriodicInterval)
    m=D.order
    C=2/(d.b-d.a)*π

    for k=kr
        if mod(m,4)==0
            A[k,k] += (C*k)^(-m)
        elseif mod(m,4)==2
            A[k,k] += -(C*k)^(-m)
        elseif k>1 && mod(m,4)==1
            A[k,k-1] += -(C*(k-1))^(-m)
        elseif k>1 && mod(m,4)==3
            A[k,k-1] += (C*(k-1))^(-m)
        end
    end

    A
end


bandinds{T,CS<:CosSpace,DD<:PeriodicInterval}(D::ConcreteIntegral{SliceSpace{1,1,CS,T,DD,1}})=(0,0)
rangespace{T,CS<:CosSpace,DD<:PeriodicInterval}(D::ConcreteIntegral{SliceSpace{1,1,CS,T,DD,1}})=iseven(D.order)?D.space:SinSpace(domain(D))

function addentries!{T,CS<:CosSpace,DD<:PeriodicInterval}(D::ConcreteIntegral{SliceSpace{1,1,CS,T,DD,1}},A,kr::Range,::Colon)
    d=domain(D)
    @assert isa(d,PeriodicInterval)
    m=D.order
    C=2/(d.b-d.a)*π

    for k=kr
        if mod(m,4)==0
            A[k,k] += (C*k)^(-m)
        elseif mod(m,4)==2
            A[k,k] += -(C*k)^(-m)
        elseif mod(m,4)==1
            A[k,k] += (C*k)^(-m)
        elseif mod(m,4)==3
            A[k,k] += -(C*k)^(-m)
        end
    end

    A
end

# CosSpace Multiplicaiton is the same as Chebyshev

bandinds{Sp<:CosSpace}(M::ConcreteMultiplication{Sp,Sp})=(1-length(M.f.coefficients),length(M.f.coefficients)-1)
rangespace{Sp<:CosSpace}(M::ConcreteMultiplication{Sp,Sp})=domainspace(M)
addentries!{Sp<:CosSpace}(M::ConcreteMultiplication{Sp,Sp},A,kr::UnitRange,::Colon)=chebmult_addentries!(M.f.coefficients,A,kr)


function addentries!{Sp<:SinSpace}(M::ConcreteMultiplication{Sp,Sp},A,kr::UnitRange,::Colon)
    a=M.f.coefficients
    toeplitz_addentries!(0.5,[0.;-a],a,A,kr)
    hankel_addentries!(0.5,a,A,max(kr[1],2):kr[end])
    A
end

bandinds{Sp<:SinSpace}(M::ConcreteMultiplication{Sp,Sp})=-length(M.f)-1,length(M.f)-1
rangespace{Sp<:SinSpace}(M::ConcreteMultiplication{Sp,Sp})=CosSpace(domain(M))


function addentries!{Sp<:SinSpace,Cs<:CosSpace}(M::ConcreteMultiplication{Sp,Cs},A,kr::Range,::Colon)
    a=M.f.coefficients
    toeplitz_addentries!(0.5,a[2:end],[a[1];0.;-a],A,kr)
    hankel_addentries!(0.5,a,A,kr)
    A
end

bandinds{Sp<:SinSpace,Cs<:CosSpace}(M::ConcreteMultiplication{Sp,Cs})=1-length(M.f),length(M.f)+1
rangespace{Sp<:SinSpace,Cs<:CosSpace}(M::ConcreteMultiplication{Sp,Cs})=SinSpace(domain(M))



function addentries!{Sp<:SinSpace,Cs<:CosSpace}(M::ConcreteMultiplication{Cs,Sp},A,kr::Range,::Colon)
    a=M.f.coefficients
    toeplitz_addentries!(0.5a,A,kr)
    if length(a)>=3
        hankel_addentries!(-0.5,a[3:end],A,kr)
    end
    A
end

bandinds{Sp<:SinSpace,Cs<:CosSpace}(M::ConcreteMultiplication{Cs,Sp})=(1-length(M.f.coefficients),length(M.f.coefficients)-1)
rangespace{Sp<:SinSpace,Cs<:CosSpace}(M::ConcreteMultiplication{Cs,Sp})=SinSpace(domain(M))



function Multiplication{T,D}(a::Fun{Fourier{D},T},sp::Fourier{D})
    d=domain(a)
    c,s=vec(a)
    O=BandedOperator{T}[Multiplication(c,CosSpace(d)) Multiplication(s,SinSpace(d));
                        Multiplication(s,CosSpace(d)) Multiplication(c,SinSpace(d))]
    MultiplicationWrapper(a,SpaceOperator(InterlaceOperator(O),space(a),sp))
end


## Definite integral

DefiniteIntegral{D<:PeriodicInterval}(sp::Fourier{D})=DefiniteIntegral{typeof(sp),Float64}(sp)
DefiniteIntegral{D<:Circle}(sp::Fourier{D})=DefiniteIntegral{typeof(sp),Complex{Float64}}(sp)

function getindex{T,D}(Σ::DefiniteIntegral{Fourier{D},T},kr::Range)
    d = domain(Σ)
    if isa(d,PeriodicInterval)
        T[k == 1?  d.b-d.a : zero(T) for k=kr]
    else
        @assert isa(d,Circle)
        T[k == 2?  -d.radius*π : (k==3?d.radius*π*im :zero(T)) for k=kr]
    end
end

datalength{D}(Σ::DefiniteIntegral{Fourier{D}})=isa(domain(Σ),PeriodicInterval)?1:3

DefiniteLineIntegral{D}(sp::Fourier{D})=DefiniteLineIntegral{typeof(sp),Float64}(sp)

function getindex{T,D}(Σ::DefiniteLineIntegral{Fourier{D},T},kr::Range)
    d = domain(Σ)
    if isa(d,PeriodicInterval)
        T[k == 1?  d.b-d.a : zero(T) for k=kr]
    else
        @assert isa(d,Circle)
        T[k == 1?  2d.radius*π : zero(T) for k=kr]
    end
end

datalength{D}(Σ::DefiniteLineIntegral{Fourier{D}})=1


transformtimes{CS<:CosSpace,D}(f::Fun{CS},g::Fun{Fourier{D}}) = transformtimes(Fun(interlace(f.coefficients,zeros(eltype(f),length(f)-1)),Fourier(domain(f))),g)
transformtimes{SS<:SinSpace,D}(f::Fun{SS},g::Fun{Fourier{D}}) = transformtimes(Fun(interlace(zeros(eltype(f),length(f)+1),f.coefficients),Fourier(domain(f))),g)
transformtimes{CS<:CosSpace,SS<:SinSpace}(f::Fun{CS},g::Fun{SS}) = transformtimes(Fun(interlace(f.coefficients,zeros(eltype(f),length(f)-1)),Fourier(domain(f))),Fun(interlace(zeros(eltype(g),length(g)+1),g.coefficients),Fourier(domain(g))))
transformtimes{CS<:CosSpace,D}(f::Fun{Fourier{D}},g::Fun{CS}) = transformtimes(g,f)
transformtimes{SS<:SinSpace,D}(f::Fun{Fourier{D}},g::Fun{SS}) = transformtimes(g,f)
transformtimes{SS<:SinSpace,CS<:CosSpace}(f::Fun{SS},g::Fun{CS}) = transformtimes(g,f)
