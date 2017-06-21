

## Converison

#ensure that COnversion is called
coefficients{DD,R1,R2}(cfs::AbstractVector,A::Fourier{DD,R1},B::Laurent{DD,R2}) =
    A_mul_B_coefficients(Conversion(A,B),cfs)
coefficients{DD,R1,R2}(cfs::AbstractVector,A::Laurent{DD,R1},B::Fourier{DD,R2}) =
    A_mul_B_coefficients(Conversion(A,B),cfs)

hasconversion{DD,R1,R2}(::Fourier{DD,R1},::Laurent{DD,R2}) = true
hasconversion{DD,R1,R2}(::Laurent{DD,R1},::Fourier{DD,R2}) = true

Conversion{DD,R1,R2}(a::Laurent{DD,R1},b::Fourier{DD,R2}) = ConcreteConversion(a,b)
Conversion{DD,R1,R2}(a::Fourier{DD,R1},b::Laurent{DD,R2}) = ConcreteConversion(a,b)

function getindex{DD,R1,R2,T}(C::ConcreteConversion{Laurent{DD,R1},Fourier{DD,R2},T},k::Integer,j::Integer)
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


function getindex{DD,R1,R2,T}(C::ConcreteConversion{Fourier{DD,R1},Laurent{DD,R2},T},k::Integer,j::Integer)
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


bandinds{DD,R1,R2}(::ConcreteConversion{Laurent{DD,R1},Fourier{DD,R2}})=-1,1
bandinds{DD,R1,R2}(::ConcreteConversion{Fourier{DD,R1},Laurent{DD,R2}})=-1,1

for RULE in (:conversion_rule,:maxspace_rule,:union_rule)
    @eval begin
        # override both to avoid SumSpace overrides
        function $RULE{DD,R1,R2}(A::Laurent{DD,R1},B::Fourier{DD,R2})
            @assert domainscompatible(A,B)
            B
        end
        function $RULE{DD,R1,R2}(A::Fourier{DD,R1},B::Laurent{DD,R2})
            @assert domainscompatible(A,B)
            A
        end
    end
end

conversion_type{DD<:Circle,R1,R2}(A::Fourier{DD,R1},B::Fourier{DD,R2}) = domain(A).orientation?A:B

hasconversion{DD,R1,R2}(A::Fourier{DD,R1},B::Fourier{DD,R2}) = domain(A) == reverse(domain(B))
function Conversion{DD,R1,R2}(A::Fourier{DD,R1},B::Fourier{DD,R2})
    if A==B
        ConversionWrapper(eye(A))
    else
        @assert domain(A) == reverse(domain(B))
        ConcreteConversion(A,B)
    end
end
bandinds{DD,R1,R2}(::ConcreteConversion{Fourier{DD,R1},Fourier{DD,R2}})=0,0

getindex{DD,R1,R2,T}(C::ConcreteConversion{Fourier{DD,R1},Fourier{DD,R2},T},k::Integer,j::Integer) =
    k==j?(iseven(k)?(-one(T)):one(T)):zero(T)





### Cos/Sine


function Derivative(S::Union{CosSpace,SinSpace},order)
    @assert isa(domain(S),PeriodicInterval)
    ConcreteDerivative(S,order)
end


bandinds{CS<:CosSpace}(D::ConcreteDerivative{CS}) = iseven(D.order)?(0,0):(0,1)
bandinds{S<:SinSpace}(D::ConcreteDerivative{S}) = iseven(D.order)?(0,0):(-1,0)
rangespace{S<:CosSpace}(D::ConcreteDerivative{S}) = iseven(D.order)?D.space:SinSpace(domain(D))
rangespace{S<:SinSpace}(D::ConcreteDerivative{S}) = iseven(D.order)?D.space:CosSpace(domain(D))


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


# Use Laurent derivative
Derivative{DD<:Circle,RR}(S::Fourier{DD,RR},k::Integer) =
    DerivativeWrapper(Derivative(Laurent(S),k)*Conversion(S,Laurent(S)),k)

Integral(::CosSpace,m::Integer) =
    error("Integral not defined for CosSpace.  Use Integral(CosSpace()|(2:∞)) if first coefficient vanishes.")

Integral{DD<:PeriodicInterval}(sp::SinSpace{DD},m::Integer) = ConcreteIntegral(sp,m)

bandinds{CS<:SinSpace}(D::ConcreteIntegral{CS}) = iseven(D.order)?(0,0):(-1,0)
rangespace{S<:CosSpace}(D::ConcreteIntegral{S}) = iseven(D.order)?D.space:SinSpace(domain(D))
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

function Integral{CS<:CosSpace,DD<:PeriodicInterval}(S::SubSpace{CS,UnitCount{Int64},DD},k::Integer)
    @assert first(S.indexes)==2
    ConcreteIntegral(S,k)
end

bandinds{CS<:CosSpace,DD<:PeriodicInterval,RR}(D::ConcreteIntegral{SubSpace{CS,UnitCount{Int64},DD,RR}}) =
    (0,0)
rangespace{CS<:CosSpace,DD<:PeriodicInterval,RR}(D::ConcreteIntegral{SubSpace{CS,UnitCount{Int64},DD,RR}}) =
    iseven(D.order)?D.space:SinSpace(domain(D))

function getindex(D::ConcreteIntegral{SubSpace{CS,UnitCount{Int64},DD,RR}},
                  k::Integer,j::Integer) where {CS<:CosSpace,DD<:PeriodicInterval,RR}
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
    if length(a) == 0
        A=ZeroOperator(sp,SinSpace(domain(sp)))
        MultiplicationWrapper(f,A)
    elseif length(a) == 1
        A=ToeplitzOperator([0.],[a[1];0.;-a]) + HankelOperator(a)
        MultiplicationWrapper(f,SpaceOperator(A,sp,SinSpace(domain(sp))))
    else
        A=ToeplitzOperator(a[2:end],[a[1];0.;-a]) + HankelOperator(a)
        MultiplicationWrapper(f,SpaceOperator(A,sp,SinSpace(domain(sp))))
    end
end


bandinds{CS<:CosSpace}(M::ConcreteMultiplication{CS,CS}) =
    (1-ncoefficients(M.f),ncoefficients(M.f)-1)
rangespace{CS<:CosSpace}(M::ConcreteMultiplication{CS,CS}) = domainspace(M)
getindex{CS<:CosSpace}(M::ConcreteMultiplication{CS,CS},k::Integer,j::Integer) =
    chebmult_getindex(M.f.coefficients,k,j)



function getindex{SS<:SinSpace}(M::ConcreteMultiplication{SS,SS},k::Integer,j::Integer)
    a=M.f.coefficients
    ret=toeplitz_getindex([zero(eltype(a));-a],a,k,j)/2
    if k ≥ 2
        ret+=hankel_getindex(a,k,j)/2
    end
    ret
end

bandinds{SS<:SinSpace}(M::ConcreteMultiplication{SS,SS})=-ncoefficients(M.f)-1,ncoefficients(M.f)-1
rangespace{SS<:SinSpace}(M::ConcreteMultiplication{SS,SS})=CosSpace(domain(M))


function getindex{SS<:SinSpace,Cs<:CosSpace}(M::ConcreteMultiplication{Cs,SS},k::Integer,j::Integer)
    a=M.f.coefficients
    ret=toeplitz_getindex(a,k,j)/2
    if length(a)>=3
        ret-=hankel_getindex(view(a,3:length(a)),k,j)/2
    end
    ret
end

bandinds{SS<:SinSpace,Cs<:CosSpace}(M::ConcreteMultiplication{Cs,SS}) =
    (1-ncoefficients(M.f),ncoefficients(M.f)-1)
rangespace{SS<:SinSpace,Cs<:CosSpace}(M::ConcreteMultiplication{Cs,SS}) =
    SinSpace(domain(M))



function Multiplication{T,D,R}(a::Fun{Fourier{D,R},T},sp::Fourier{D,R})
    d=domain(a)
    c,s=components(a)
    O=Operator{T}[Multiplication(c,CosSpace(d)) Multiplication(s,SinSpace(d));
                        Multiplication(s,CosSpace(d)) Multiplication(c,SinSpace(d))]
    MultiplicationWrapper(a,InterlaceOperator(O,space(a),sp))
end


## Definite integral

for SP in (:CosSpace,:SinSpace,:Fourier)
    @eval begin
        DefiniteIntegral{D,R}(S::$SP{D,R}) =
            ConcreteDefiniteIntegral{typeof(S),prectype(S)}(S)
        DefiniteLineIntegral{D,R}(S::$SP{D,R}) =
            ConcreteDefiniteLineIntegral{typeof(S),real(prectype(S))}(S)
    end
end

getindex{T,D<:PeriodicInterval,R}(Σ::ConcreteDefiniteIntegral{CosSpace{D,R},T},k::Integer) =
    k == 1? T(complexlength(domain(Σ))) : zero(T)

getindex{T,D<:PeriodicInterval,R}(Σ::ConcreteDefiniteIntegral{SinSpace{D,R},T},k::Integer) =
    zero(T)

getindex{T,D<:PeriodicInterval,R}(Σ::ConcreteDefiniteIntegral{Fourier{D,R},T},k::Integer) =
    k == 1? T(complexlength(domain(Σ))) : zero(T)

getindex{T,D<:Circle,R}(Σ::ConcreteDefiniteIntegral{CosSpace{D,R},T},k::Integer) =
    k==2? T(complexlength(domain(Σ))/2) : zero(T)

getindex{T,D<:Circle,R}(Σ::ConcreteDefiniteIntegral{SinSpace{D,R},T},k::Integer) =
    k == 1? T(0.5im*complexlength(domain(Σ))) : zero(T)

getindex{T,D<:Circle,R}(Σ::ConcreteDefiniteIntegral{Fourier{D,R},T},k::Integer) =
    k == 2? T(0.5im*complexlength(domain(Σ))) : (k==3 ? T(complexlength(domain(Σ))/2) : zero(T))

getindex{T,D<:PeriodicInterval,R}(Σ::ConcreteDefiniteLineIntegral{CosSpace{D,R},T},k::Integer) =
    k==1? T(arclength(domain(Σ))) : zero(T)

getindex{T,D<:PeriodicInterval,R}(Σ::ConcreteDefiniteLineIntegral{SinSpace{D,R},T},k::Integer) =
    zero(T)

getindex{T,D<:PeriodicInterval,R}(Σ::ConcreteDefiniteLineIntegral{Fourier{D,R},T},k::Integer) =
    k==1? T(arclength(domain(Σ))) : zero(T)

getindex{T,D<:Circle,R}(Σ::ConcreteDefiniteLineIntegral{CosSpace{D,R},T},k::Integer) =
    k==1? T(arclength(domain(Σ))) : zero(T)

getindex{T,D<:Circle,R}(Σ::ConcreteDefiniteLineIntegral{SinSpace{D,R},T},k::Integer) =
    zero(T)

getindex{T,D<:Circle,R}(Σ::ConcreteDefiniteLineIntegral{Fourier{D,R},T},k::Integer) =
    k==1? T(arclength(domain(Σ))) : zero(T)

bandinds{D<:PeriodicInterval,R}(Σ::ConcreteDefiniteIntegral{CosSpace{D,R}}) = 0,0
bandinds{D<:PeriodicInterval,R}(Σ::ConcreteDefiniteIntegral{SinSpace{D,R}}) = 0,0
bandinds{D<:PeriodicInterval,R}(Σ::ConcreteDefiniteIntegral{Fourier{D,R}}) = 0,0
bandinds{D<:Circle,R}(Σ::ConcreteDefiniteIntegral{CosSpace{D,R}}) = 0,1
bandinds{D<:Circle,R}(Σ::ConcreteDefiniteIntegral{SinSpace{D,R}}) = 0,0
bandinds{D<:Circle,R}(Σ::ConcreteDefiniteIntegral{Fourier{D,R}}) = 0,2
bandinds{D,R}(Σ::ConcreteDefiniteLineIntegral{CosSpace{D,R}}) = 0,0
bandinds{D,R}(Σ::ConcreteDefiniteLineIntegral{SinSpace{D,R}}) = 0,0
bandinds{D,R}(Σ::ConcreteDefiniteLineIntegral{Fourier{D,R}}) = 0,0


transformtimes{CS<:CosSpace,D,R}(f::Fun{CS},g::Fun{Fourier{D,R}}) =
    transformtimes(Fun(Fourier(domain(f)),interlace(f.coefficients,zeros(eltype(f),ncoefficients(f)-1))),g)
transformtimes{SS<:SinSpace,D,R}(f::Fun{SS},g::Fun{Fourier{D,R}}) =
    transformtimes(Fun(Fourier(domain(f)),interlace(zeros(eltype(f),ncoefficients(f)+1),f.coefficients)),g)
transformtimes{CS<:CosSpace,SS<:SinSpace}(f::Fun{CS},g::Fun{SS}) =
    transformtimes(Fun(Fourier(domain(f)),interlace(f.coefficients,zeros(eltype(f),ncoefficients(f)-1))),
                    Fun(Fourier(domain(g)),interlace(zeros(eltype(g),ncoefficients(g)+1),g.coefficients)))
transformtimes{CS<:CosSpace,D,R}(f::Fun{Fourier{D,R}},g::Fun{CS}) = transformtimes(g,f)
transformtimes{SS<:SinSpace,D,R}(f::Fun{Fourier{D,R}},g::Fun{SS}) = transformtimes(g,f)
transformtimes{SS<:SinSpace,CS<:CosSpace}(f::Fun{SS},g::Fun{CS}) = transformtimes(g,f)


ReverseOrientation{D}(S::Fourier{D}) = ReverseOrientationWrapper(SpaceOperator(NegateEven(),S,reverseorientation(S)))
Reverse{D}(S::Fourier{D}) = ReverseWrapper(SpaceOperator(NegateEven(),S,S))




## Multivariate


for TYP in (:Fourier,:Laurent,:CosSpace,:SinSpace,:Taylor)
    @eval begin
        function Dirichlet(S::TensorSpace{Tuple{$TYP{PeriodicInterval{T},R},PS}},k) where {PS,T,R}
            op = [eye(S.spaces[1])⊗Evaluation(S.spaces[2],first,k);
                            ReverseOrientation(S.spaces[1])⊗Evaluation(S.spaces[2],last,k) ]
            DirichletWrapper(SpaceOperator(op,S,PiecewiseSpace(rangespace(op).spaces)),k)
        end
        function Dirichlet(S::TensorSpace{Tuple{PS,$TYP{PeriodicInterval{T},R}}},k) where {PS,T,R}
            op = [Evaluation(S.spaces[1],first,k)⊗eye(S.spaces[2]);
                            Evaluation(S.spaces[1],last,k)⊗ReverseOrientation(S[2]) ]
            DirichletWrapper(SpaceOperator(op,S,PiecewiseSpace(rangespace(op).spaces)),k)
        end
    end
end
