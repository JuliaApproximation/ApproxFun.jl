

## Converison

#ensure that COnversion is called
coefficients{DD}(cfs::Vector,A::Fourier{DD},B::Laurent{DD})=(Conversion(A,B)*cfs).coefficients
coefficients{DD}(cfs::Vector,A::Laurent{DD},B::Fourier{DD})=(Conversion(A,B)*cfs).coefficients

hasconversion{DD}(::Fourier{DD},::Laurent{DD})=true
hasconversion{DD}(::Laurent{DD},::Fourier{DD})=true

Conversion{DD}(a::Laurent{DD},b::Fourier{DD})=ConcreteConversion(a,b)
Conversion{DD}(a::Fourier{DD},b::Laurent{DD})=ConcreteConversion(a,b)

function getindex{DD,T}(C::ConcreteConversion{Laurent{DD},Fourier{DD},T},k::Integer,j::Integer)
    if k==j==1
        one(T)
    elseif iseven(k) && k==j
        -one(T)*im
    elseif iseven(k) && k+1==j
        one(T)*im
    elseif isodd(k) && (k==j || k-1==j )
        one(T)
    else
        zero(T)
    end
end


function getindex{DD,T}(C::ConcreteConversion{Fourier{DD},Laurent{DD},T},k::Integer,j::Integer)
    if k==j==1
        one(T)
    elseif iseven(k) && k==j
        one(T)/2*im
    elseif iseven(k) && k+1==j
        one(T)/2
    elseif isodd(k) && k==j
        one(T)/2
    elseif isodd(k) && j==k-1
        -one(T)*im/2
    else
        zero(T)
    end
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
function Conversion{DD}(A::Fourier{DD},B::Fourier{DD})
    if A==B
        ConversionWrapper(eye(A))
    else
        @assert domain(A) == reverse(domain(B))
        ConcreteConversion(A,B)
    end
end
bandinds{DD}(::ConcreteConversion{Fourier{DD},Fourier{DD}})=0,0

getindex{DD,T}(C::ConcreteConversion{Fourier{DD},Fourier{DD},T},k::Integer,j::Integer) =
    k==j?(iseven(k)?(-one(T)):one(T)):zero(T)





### Cos/Sine


function Derivative(S::Union{CosSpace,SinSpace},order)
    @assert isa(domain(S),PeriodicInterval)
    ConcreteDerivative(S,order)
end


bandinds{CS<:CosSpace}(D::ConcreteDerivative{CS})=iseven(D.order)?(0,0):(0,1)
bandinds{S<:SinSpace}(D::ConcreteDerivative{S})=iseven(D.order)?(0,0):(-1,0)
rangespace{S<:CosSpace}(D::ConcreteDerivative{S})=iseven(D.order)?D.space:SinSpace(domain(D))
rangespace{S<:SinSpace}(D::ConcreteDerivative{S})=iseven(D.order)?D.space:CosSpace(domain(D))


function getindex{CS<:CosSpace,OT,T}(D::ConcreteDerivative{CS,OT,T},k::Integer,j::Integer)
    d=domain(D)
    m=D.order
    C=T(2/(d.b-d.a)*π)

    if k==j && mod(m,4)==0
        (C*(k-1))^m
    elseif k==j && mod(m,4)==2
        -(C*(k-1))^m
    elseif j==k+1 && mod(m,4)==1
        -(C*k)^m
    elseif j==k+1 && mod(m,4)==3
        (C*k)^m
    else
        zero(T)
    end
end

function getindex{CS<:SinSpace,OT,T}(D::ConcreteDerivative{CS,OT,T},k::Integer,j::Integer)
    d=domain(D)
    m=D.order
    C=T(2/(d.b-d.a)*π)

    if k==j && mod(m,4)==0
        (C*k)^m
    elseif k==j && mod(m,4)==2
        -(C*k)^m
    elseif j==k-1 && mod(m,4)==1
        (C*j)^m
    elseif j==k-1 && mod(m,4)==3
        -(C*j)^m
    else
        zero(T)
    end
end

Integral(::CosSpace,m::Integer) =
    error("Integral not defined for CosSpace.  Use Integral(CosSpace()|(2:∞)) if first coefficient vanishes.")

Integral{DD<:PeriodicInterval}(sp::SinSpace{DD},m::Integer) = ConcreteIntegral(sp,m)

bandinds{CS<:SinSpace}(D::ConcreteIntegral{CS})=iseven(D.order)?(0,0):(-1,0)
rangespace{S<:CosSpace}(D::ConcreteIntegral{S})=iseven(D.order)?D.space:SinSpace(domain(D))
rangespace{S<:SinSpace}(D::ConcreteIntegral{S})=iseven(D.order)?D.space:CosSpace(domain(D))

function getindex{CS<:SinSpace,OT,T}(D::ConcreteIntegral{CS,OT,T},k::Integer,j::Integer)
    d=domain(D)
    @assert isa(d,PeriodicInterval)
    m=D.order
    C=T(2/(d.b-d.a)*π)


    if k==j && mod(m,4)==0
        (C*k)^(-m)
    elseif k==j && mod(m,4)==2
        -(C*k)^(-m)
    elseif j==k-1 && mod(m,4)==1
        -(C*j)^(-m)
    elseif j==k-1 && mod(m,4)==3
        (C*j)^(-m)
    else
        zero(T)
    end
end

function Integral{T,CS<:CosSpace,DD<:PeriodicInterval}(S::SubSpace{CS,UnitCount{Int64},T,DD,1},k::Integer)
    @assert first(S.indexes)==2
    ConcreteIntegral(S,k)
end

bandinds{T,CS<:CosSpace,DD<:PeriodicInterval}(D::ConcreteIntegral{SubSpace{CS,UnitCount{Int64},T,DD,1}}) =
    (0,0)
rangespace{T,CS<:CosSpace,DD<:PeriodicInterval}(D::ConcreteIntegral{SubSpace{CS,UnitCount{Int64},T,DD,1}}) =
    iseven(D.order)?D.space:SinSpace(domain(D))

function getindex{TT,CS<:CosSpace,DD<:PeriodicInterval}(D::ConcreteIntegral{SubSpace{CS,UnitCount{Int64},TT,DD,1}},
                                                            k::Integer,j::Integer)
    d=domain(D)
    m=D.order
    T=eltype(D)
    C=T(2/(d.b-d.a)*π)


    if k==j
        if mod(m,4)==0
            (C*k)^(-m)
        elseif mod(m,4)==2
            -(C*k)^(-m)
        elseif mod(m,4)==1
        (C*k)^(-m)
        else   # mod(m,4)==3
            -(C*k)^(-m)
        end
    else
        zero(T)
    end
end

# CosSpace Multiplicaiton is the same as Chebyshev


Multiplication{CS<:CosSpace}(f::Fun{CS},sp::CS) = ConcreteMultiplication(f,sp)
Multiplication{SS<:SinSpace}(f::Fun{SS},sp::SS) = ConcreteMultiplication(f,sp)
Multiplication{CS<:CosSpace}(f::Fun{CS},sp::SinSpace) = ConcreteMultiplication(f,sp)
function Multiplication{SS<:SinSpace}(f::Fun{SS},sp::CosSpace)
    @assert domain(f) == domain(sp)
    a=f.coefficients/2
    A=ToeplitzOperator(a[2:end],[a[1];0.;-a]) + HankelOperator(a)
    MultiplicationWrapper(f,SpaceOperator(A,sp,SinSpace(domain(sp))))
end


bandinds{CS<:CosSpace}(M::ConcreteMultiplication{CS,CS}) =
    (1-ncoefficients(M.f),ncoefficients(M.f)-1)
rangespace{CS<:CosSpace}(M::ConcreteMultiplication{CS,CS}) = domainspace(M)
getindex{CS<:CosSpace}(M::ConcreteMultiplication{CS,CS},k::Integer,j::Integer) =
    chebmult_getindex(M.f.coefficients,k,j)



function getindex{SS<:SinSpace}(M::ConcreteMultiplication{SS,SS},k::Integer,j::Integer)
    a=M.f.coefficients
    ret=0.5*toeplitz_getindex([0.;-a],a,k,j)
    if k ≥ 2
        ret+=0.5*hankel_getindex(a,k,j)
    end
    ret
end

bandinds{SS<:SinSpace}(M::ConcreteMultiplication{SS,SS})=-ncoefficients(M.f)-1,ncoefficients(M.f)-1
rangespace{SS<:SinSpace}(M::ConcreteMultiplication{SS,SS})=CosSpace(domain(M))


function getindex{SS<:SinSpace,Cs<:CosSpace}(M::ConcreteMultiplication{Cs,SS},k::Integer,j::Integer)
    a=M.f.coefficients
    ret=0.5*toeplitz_getindex(a,k,j)
    if length(a)>=3
        ret-=0.5*hankel_getindex(@compat(view(a,3:length(a))),k,j)
    end
    ret
end

bandinds{SS<:SinSpace,Cs<:CosSpace}(M::ConcreteMultiplication{Cs,SS}) =
    (1-ncoefficients(M.f),ncoefficients(M.f)-1)
rangespace{SS<:SinSpace,Cs<:CosSpace}(M::ConcreteMultiplication{Cs,SS}) =
    SinSpace(domain(M))



function Multiplication{T,D}(a::Fun{Fourier{D},T},sp::Fourier{D})
    d=domain(a)
    c,s=vec(a)
    O=Operator{T}[Multiplication(c,CosSpace(d)) Multiplication(s,SinSpace(d));
                        Multiplication(s,CosSpace(d)) Multiplication(c,SinSpace(d))]
    MultiplicationWrapper(a,InterlaceOperator(O,space(a),sp))
end


## Definite integral

DefiniteIntegral{D<:PeriodicInterval}(sp::Fourier{D}) =
    ConcreteDefiniteIntegral{typeof(sp),Float64}(sp)
DefiniteIntegral{D<:Circle}(sp::Fourier{D}) =
    ConcreteDefiniteIntegral{typeof(sp),Complex{Float64}}(sp)

function getindex{T,D}(Σ::ConcreteDefiniteIntegral{Fourier{D},T},kr::Range)
    d = domain(Σ)
    if isa(d,PeriodicInterval)
        T[k == 1?  d.b-d.a : zero(T) for k=kr]
    else
        @assert isa(d,Circle)
        T[k == 2?  -d.radius*π : (k==3?d.radius*π*im :zero(T)) for k=kr]
    end
end

bandinds{D}(Σ::ConcreteDefiniteIntegral{Fourier{D}}) =
    0,isa(domain(Σ),PeriodicInterval)?0:2

DefiniteLineIntegral{D}(sp::Fourier{D}) =
    ConcreteDefiniteLineIntegral{typeof(sp),Float64}(sp)


getindex{T,D<:PeriodicInterval}(Σ::ConcreteDefiniteLineIntegral{Fourier{D},T},k::Integer) =
    k==1? domain(Σ).b-domain(Σ).a : zero(T)

getindex{T,D<:Circle}(Σ::ConcreteDefiniteLineIntegral{Fourier{D},T},k::Integer) =
    k==1? domain(Σ).radius*π : zero(T)


bandinds{D}(Σ::ConcreteDefiniteLineIntegral{Fourier{D}}) = 0,0


transformtimes{CS<:CosSpace,D}(f::Fun{CS},g::Fun{Fourier{D}}) =
    transformtimes(Fun(interlace(f.coefficients,zeros(eltype(f),ncoefficients(f)-1)),Fourier(domain(f))),g)
transformtimes{SS<:SinSpace,D}(f::Fun{SS},g::Fun{Fourier{D}}) =
    transformtimes(Fun(interlace(zeros(eltype(f),ncoefficients(f)+1),f.coefficients),Fourier(domain(f))),g)
transformtimes{CS<:CosSpace,SS<:SinSpace}(f::Fun{CS},g::Fun{SS}) =
    transformtimes(Fun(interlace(f.coefficients,zeros(eltype(f),ncoefficients(f)-1)),Fourier(domain(f))),
                    Fun(interlace(zeros(eltype(g),ncoefficients(g)+1),g.coefficients),Fourier(domain(g))))
transformtimes{CS<:CosSpace,D}(f::Fun{Fourier{D}},g::Fun{CS}) = transformtimes(g,f)
transformtimes{SS<:SinSpace,D}(f::Fun{Fourier{D}},g::Fun{SS}) = transformtimes(g,f)
transformtimes{SS<:SinSpace,CS<:CosSpace}(f::Fun{SS},g::Fun{CS}) = transformtimes(g,f)


ReverseOrientation{D}(S::Fourier{D}) = ReverseOrientationWrapper(SpaceOperator(NegateEven(),S,reverseorientation(S)))
Reverse{D}(S::Fourier{D}) = ReverseWrapper(SpaceOperator(NegateEven(),S,S))




## Multivariate


for TYP in (:Fourier,:Laurent,:CosSpace,:SinSpace,:Taylor)
    @eval begin
        function Dirichlet{PS,T}(S::TensorSpace{Tuple{$TYP{PeriodicInterval{T}},PS}})
            op = interlace([eye(S[1])⊗ldirichlet(S[2]);
                            ReverseOrientation(S[1])⊗rdirichlet(S[2]) ])
            DirichletWrapper(SpaceOperator(op,S,PiecewiseSpace(rangespace(op).spaces)),1)
        end
        function Dirichlet{PS,T}(S::TensorSpace{Tuple{PS,$TYP{PeriodicInterval{T}}}})
            op = interlace([ldirichlet(S[1])⊗eye(S[2]);
                            rdirichlet(S[1])⊗ReverseOrientation(S[2]) ])
            DirichletWrapper(SpaceOperator(op,S,PiecewiseSpace(rangespace(op).spaces)),1)
        end
    end
end
