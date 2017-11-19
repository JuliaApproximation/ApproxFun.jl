

## Converison

#ensure that COnversion is called
coefficients(cfs::AbstractVector,A::Fourier{DD,R1},B::Laurent{DD,R2}) where {DD,R1,R2} =
    A_mul_B_coefficients(Conversion(A,B),cfs)
coefficients(cfs::AbstractVector,A::Laurent{DD,R1},B::Fourier{DD,R2}) where {DD,R1,R2} =
    A_mul_B_coefficients(Conversion(A,B),cfs)

hasconversion(::Fourier{DD,R1},::Laurent{DD,R2}) where {DD,R1,R2} = true
hasconversion(::Laurent{DD,R1},::Fourier{DD,R2}) where {DD,R1,R2} = true

Conversion(a::Laurent{DD,R1},b::Fourier{DD,R2}) where {DD,R1,R2} = ConcreteConversion(a,b)
Conversion(a::Fourier{DD,R1},b::Laurent{DD,R2}) where {DD,R1,R2} = ConcreteConversion(a,b)

function getindex(C::ConcreteConversion{Laurent{DD,R1},Fourier{DD,R2},T},k::Integer,j::Integer) where {DD,R1,R2,T}
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


function getindex(C::ConcreteConversion{Fourier{DD,R1},Laurent{DD,R2},T},k::Integer,j::Integer) where {DD,R1,R2,T}
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


bandinds(::ConcreteConversion{Laurent{DD,R1},Fourier{DD,R2}}) where {DD,R1,R2}=-1,1
bandinds(::ConcreteConversion{Fourier{DD,R1},Laurent{DD,R2}}) where {DD,R1,R2}=-1,1

for RULE in (:conversion_rule,:maxspace_rule,:union_rule)
    @eval begin
        # override both to avoid SumSpace overrides
        function $RULE(A::Laurent{DD,R1},B::Fourier{DD,R2}) where {DD,R1,R2}
            @assert domainscompatible(A,B)
            B
        end
        function $RULE(A::Fourier{DD,R1},B::Laurent{DD,R2}) where {DD,R1,R2}
            @assert domainscompatible(A,B)
            A
        end
    end
end

conversion_type(A::Fourier{DD,R1},B::Fourier{DD,R2}) where {DD<:Circle,R1,R2} = domain(A).orientation?A:B

hasconversion(A::Fourier{DD,R1},B::Fourier{DD,R2}) where {DD,R1,R2} = domain(A) == reverse(domain(B))
function Conversion(A::Fourier{DD,R1},B::Fourier{DD,R2}) where {DD,R1,R2}
    if A==B
        ConversionWrapper(eye(A))
    else
        @assert domain(A) == reverse(domain(B))
        ConcreteConversion(A,B)
    end
end
bandinds(::ConcreteConversion{Fourier{DD,R1},Fourier{DD,R2}}) where {DD,R1,R2}=0,0

getindex(C::ConcreteConversion{Fourier{DD,R1},Fourier{DD,R2},T},k::Integer,j::Integer) where {DD,R1,R2,T} =
    k==j?(iseven(k)?(-one(T)):one(T)):zero(T)





### Cos/Sine


function Derivative(S::Union{CosSpace,SinSpace},order)
    @assert isa(domain(S),PeriodicInterval)
    ConcreteDerivative(S,order)
end


bandinds(D::ConcreteDerivative{CS}) where {CS<:CosSpace} = iseven(D.order)?(0,0):(0,1)
bandinds(D::ConcreteDerivative{S}) where {S<:SinSpace} = iseven(D.order)?(0,0):(-1,0)
rangespace(D::ConcreteDerivative{S}) where {S<:CosSpace} = iseven(D.order)?D.space:SinSpace(domain(D))
rangespace(D::ConcreteDerivative{S}) where {S<:SinSpace} = iseven(D.order)?D.space:CosSpace(domain(D))


function getindex(D::ConcreteDerivative{CS,OT,T},k::Integer,j::Integer) where {CS<:CosSpace,OT,T}
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

function getindex(D::ConcreteDerivative{CS,OT,T},k::Integer,j::Integer) where {CS<:SinSpace,OT,T}
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
Derivative(S::Fourier{DD,RR},k::Integer) where {DD<:Circle,RR} =
    DerivativeWrapper(Derivative(Laurent(S),k)*Conversion(S,Laurent(S)),k)

Integral(::CosSpace,m::Integer) =
    error("Integral not defined for CosSpace.  Use Integral(CosSpace()|(2:∞)) if first coefficient vanishes.")

Integral(sp::SinSpace{DD},m::Integer) where {DD<:PeriodicInterval} = ConcreteIntegral(sp,m)

bandinds(D::ConcreteIntegral{CS}) where {CS<:SinSpace} = iseven(D.order)?(0,0):(-1,0)
rangespace(D::ConcreteIntegral{S}) where {S<:CosSpace} = iseven(D.order)?D.space:SinSpace(domain(D))
rangespace(D::ConcreteIntegral{S}) where {S<:SinSpace}=iseven(D.order)?D.space:CosSpace(domain(D))

function getindex(D::ConcreteIntegral{CS,OT,T},k::Integer,j::Integer) where {CS<:SinSpace,OT,T}
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

function Integral(S::SubSpace{CS,UnitCount{Int64},DD},k::Integer) where {CS<:CosSpace,DD<:PeriodicInterval}
    @assert first(S.indexes)==2
    ConcreteIntegral(S,k)
end

bandinds(D::ConcreteIntegral{SubSpace{CS,UnitCount{Int64},DD,RR}}) where {CS<:CosSpace,DD<:PeriodicInterval,RR} =
    (0,0)
rangespace(D::ConcreteIntegral{SubSpace{CS,UnitCount{Int64},DD,RR}}) where {CS<:CosSpace,DD<:PeriodicInterval,RR} =
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


Multiplication(f::Fun{CS},sp::CS) where {CS<:CosSpace} = ConcreteMultiplication(f,sp)
Multiplication(f::Fun{SS},sp::SS) where {SS<:SinSpace} = ConcreteMultiplication(f,sp)
Multiplication(f::Fun{CS},sp::SinSpace) where {CS<:CosSpace} = ConcreteMultiplication(f,sp)
function Multiplication(f::Fun{SS},sp::CosSpace) where SS<:SinSpace
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


bandinds(M::ConcreteMultiplication{CS,CS}) where {CS<:CosSpace} =
    (1-ncoefficients(M.f),ncoefficients(M.f)-1)
rangespace(M::ConcreteMultiplication{CS,CS}) where {CS<:CosSpace} = domainspace(M)
getindex(M::ConcreteMultiplication{CS,CS},k::Integer,j::Integer) where {CS<:CosSpace} =
    chebmult_getindex(M.f.coefficients,k,j)



function getindex(M::ConcreteMultiplication{SS,SS},k::Integer,j::Integer) where SS<:SinSpace
    a=M.f.coefficients
    ret=toeplitz_getindex([zero(eltype(a));-a],a,k,j)/2
    if k ≥ 2
        ret+=hankel_getindex(a,k,j)/2
    end
    ret
end

bandinds(M::ConcreteMultiplication{SS,SS}) where {SS<:SinSpace}=-ncoefficients(M.f)-1,ncoefficients(M.f)-1
rangespace(M::ConcreteMultiplication{SS,SS}) where {SS<:SinSpace}=CosSpace(domain(M))


function getindex(M::ConcreteMultiplication{Cs,SS},k::Integer,j::Integer) where {SS<:SinSpace,Cs<:CosSpace}
    a=M.f.coefficients
    ret=toeplitz_getindex(a,k,j)/2
    if length(a)>=3
        ret-=hankel_getindex(view(a,3:length(a)),k,j)/2
    end
    ret
end

bandinds(M::ConcreteMultiplication{Cs,SS}) where {SS<:SinSpace,Cs<:CosSpace} =
    (1-ncoefficients(M.f),ncoefficients(M.f)-1)
rangespace(M::ConcreteMultiplication{Cs,SS}) where {SS<:SinSpace,Cs<:CosSpace} =
    SinSpace(domain(M))



function Multiplication(a::Fun{Fourier{D,R},T},sp::Fourier{D,R}) where {T,D,R}
    d=domain(a)
    c,s=components(a)
    O=Operator{float(T)}[Multiplication(c,CosSpace(d)) Multiplication(s,SinSpace(d));
                        Multiplication(s,CosSpace(d)) Multiplication(c,SinSpace(d))]
    MultiplicationWrapper(a,InterlaceOperator(O,space(a),sp))
end

coefficienttimes(f::Fun{Fourier{DD,RR}},g::Fun{Fourier{DD,RR}}) where {DD,RR} = Multiplication(f,space(g))*g


## Definite integral

for SP in (:CosSpace,:SinSpace,:Fourier)
    @eval begin
        DefiniteIntegral(S::$SP{D,R}) where {D,R} =
            ConcreteDefiniteIntegral{typeof(S),prectype(S)}(S)
        DefiniteLineIntegral(S::$SP{D,R}) where {D,R} =
            ConcreteDefiniteLineIntegral{typeof(S),real(prectype(S))}(S)
    end
end

getindex(Σ::ConcreteDefiniteIntegral{CosSpace{D,R},T},k::Integer) where {T,D<:PeriodicInterval,R} =
    k == 1? T(complexlength(domain(Σ))) : zero(T)

getindex(Σ::ConcreteDefiniteIntegral{SinSpace{D,R},T},k::Integer) where {T,D<:PeriodicInterval,R} =
    zero(T)

getindex(Σ::ConcreteDefiniteIntegral{Fourier{D,R},T},k::Integer) where {T,D<:PeriodicInterval,R} =
    k == 1? T(complexlength(domain(Σ))) : zero(T)

getindex(Σ::ConcreteDefiniteIntegral{CosSpace{D,R},T},k::Integer) where {T,D<:Circle,R} =
    k==2? T(complexlength(domain(Σ))/2) : zero(T)

getindex(Σ::ConcreteDefiniteIntegral{SinSpace{D,R},T},k::Integer) where {T,D<:Circle,R} =
    k == 1? T(0.5im*complexlength(domain(Σ))) : zero(T)

getindex(Σ::ConcreteDefiniteIntegral{Fourier{D,R},T},k::Integer) where {T,D<:Circle,R} =
    k == 2? T(0.5im*complexlength(domain(Σ))) : (k==3 ? T(complexlength(domain(Σ))/2) : zero(T))

getindex(Σ::ConcreteDefiniteLineIntegral{CosSpace{D,R},T},k::Integer) where {T,D<:PeriodicInterval,R} =
    k==1? T(arclength(domain(Σ))) : zero(T)

getindex(Σ::ConcreteDefiniteLineIntegral{SinSpace{D,R},T},k::Integer) where {T,D<:PeriodicInterval,R} =
    zero(T)

getindex(Σ::ConcreteDefiniteLineIntegral{Fourier{D,R},T},k::Integer) where {T,D<:PeriodicInterval,R} =
    k==1? T(arclength(domain(Σ))) : zero(T)

getindex(Σ::ConcreteDefiniteLineIntegral{CosSpace{D,R},T},k::Integer) where {T,D<:Circle,R} =
    k==1? T(arclength(domain(Σ))) : zero(T)

getindex(Σ::ConcreteDefiniteLineIntegral{SinSpace{D,R},T},k::Integer) where {T,D<:Circle,R} =
    zero(T)

getindex(Σ::ConcreteDefiniteLineIntegral{Fourier{D,R},T},k::Integer) where {T,D<:Circle,R} =
    k==1? T(arclength(domain(Σ))) : zero(T)

bandinds(Σ::ConcreteDefiniteIntegral{CosSpace{D,R}}) where {D<:PeriodicInterval,R} = 0,0
bandinds(Σ::ConcreteDefiniteIntegral{SinSpace{D,R}}) where {D<:PeriodicInterval,R} = 0,0
bandinds(Σ::ConcreteDefiniteIntegral{Fourier{D,R}}) where {D<:PeriodicInterval,R} = 0,0
bandinds(Σ::ConcreteDefiniteIntegral{CosSpace{D,R}}) where {D<:Circle,R} = 0,1
bandinds(Σ::ConcreteDefiniteIntegral{SinSpace{D,R}}) where {D<:Circle,R} = 0,0
bandinds(Σ::ConcreteDefiniteIntegral{Fourier{D,R}}) where {D<:Circle,R} = 0,2
bandinds(Σ::ConcreteDefiniteLineIntegral{CosSpace{D,R}}) where {D,R} = 0,0
bandinds(Σ::ConcreteDefiniteLineIntegral{SinSpace{D,R}}) where {D,R} = 0,0
bandinds(Σ::ConcreteDefiniteLineIntegral{Fourier{D,R}}) where {D,R} = 0,0


transformtimes(f::Fun{CS},g::Fun{Fourier{D,R}}) where {CS<:CosSpace,D,R} =
    transformtimes(Fun(Fourier(domain(f)),interlace(f.coefficients,zeros(eltype(f),ncoefficients(f)-1))),g)
transformtimes(f::Fun{SS},g::Fun{Fourier{D,R}}) where {SS<:SinSpace,D,R} =
    transformtimes(Fun(Fourier(domain(f)),interlace(zeros(eltype(f),ncoefficients(f)+1),f.coefficients)),g)
transformtimes(f::Fun{CS},g::Fun{SS}) where {CS<:CosSpace,SS<:SinSpace} =
    transformtimes(Fun(Fourier(domain(f)),interlace(f.coefficients,zeros(eltype(f),ncoefficients(f)-1))),
                    Fun(Fourier(domain(g)),interlace(zeros(eltype(g),ncoefficients(g)+1),g.coefficients)))
transformtimes(f::Fun{Fourier{D,R}},g::Fun{CS}) where {CS<:CosSpace,D,R} = transformtimes(g,f)
transformtimes(f::Fun{Fourier{D,R}},g::Fun{SS}) where {SS<:SinSpace,D,R} = transformtimes(g,f)
transformtimes(f::Fun{SS},g::Fun{CS}) where {SS<:SinSpace,CS<:CosSpace} = transformtimes(g,f)


ReverseOrientation(S::Fourier{D}) where {D} = ReverseOrientationWrapper(NegateEven(S,reverseorientation(S)))
Reverse(S::Fourier{D}) where {D} = ReverseWrapper(NegateEven(S,S))




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
