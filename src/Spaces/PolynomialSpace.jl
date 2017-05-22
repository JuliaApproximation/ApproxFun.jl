
## Orthogonal polynomials

abstract type PolynomialSpace{D} <: RealUnivariateSpace{D} end



Multiplication{U<:PolynomialSpace}(f::Fun{U},sp::PolynomialSpace) = ConcreteMultiplication(f,sp)
bandinds{U<:PolynomialSpace,V<:PolynomialSpace}(M::ConcreteMultiplication{U,V}) =
    (1-ncoefficients(M.f),ncoefficients(M.f)-1)
rangespace{U<:PolynomialSpace,V<:PolynomialSpace}(M::ConcreteMultiplication{U,V}) = domainspace(M)


# All polynomials contain constant
union_rule(A::ConstantSpace,B::PolynomialSpace) = B
Base.promote_rule{T<:Number,S<:PolynomialSpace,V,VV}(::Type{Fun{S,V,VV}},::Type{T}) =
    VFun{S,promote_type(V,T)}
Base.promote_rule{T<:Number,S<:PolynomialSpace}(::Type{Fun{S}},::Type{T}) = VFun{S,T}

## Evaluation

function evaluate(f::AbstractVector,S::PolynomialSpace,x)
    if x in domain(S)
        clenshaw(S,f,tocanonical(S,x))
    else
        zero(eltype(f))
    end
end

# we need the ... for multi-dimensional
evaluate(f::AbstractVector,S::PolynomialSpace,x,y,z...) =
    evaluate(f,S,Vec(x,y,z...))

function evaluate(f::AbstractVector,S::PolynomialSpace,x::Fun)
    if issubset(Interval(minimum(x),maximum(x)),domain(S))
        clenshaw(S,f,tocanonical(S,x))
    else
        error("Implement splitatpoints for evaluate ")
    end
end

## Extrapolation
extrapolate(f::AbstractVector,S::PolynomialSpace,x) = clenshaw(S,f,tocanonical(S,x))

######
# Recurrence encodes the recurrence coefficients
# or equivalently multiplication by x
######
immutable Recurrence{S,T} <: TridiagonalOperator{T}
    space::S
end

Recurrence(sp) = Recurrence{typeof(sp),promote_type(eltype(sp),eltype(domain(sp)))}(sp)

Base.convert{T,S}(::Type{Operator{T}},J::Recurrence{S}) = Recurrence{S,T}(J.space)

domainspace(R::Recurrence) = R.space
rangespace(R::Recurrence) = R.space


#####
# recα/β/γ are given by
#       x p_{n-1} =γ_n p_{n-2} + α_n p_{n-1} +  p_n β_n
#####

function getindex{S,T}(R::Recurrence{S,T},k::Integer,j::Integer)
    if j==k-1
        recβ(T,R.space,k-1)
    elseif j==k
        recα(T,R.space,k)
    elseif j==k+1
        recγ(T,R.space,k+1)
    else
        zero(T)
    end
end

######
# JacobiZ encodes [BasisFunctional(1);(J-z*I)[2:end,:]]
# where J is the Jacobi operator
######
immutable JacobiZ{S<:Space,T} <: TridiagonalOperator{T}
    space::S
    z::T
end

JacobiZ(sp::PolynomialSpace,z) =
    (T = promote_type(eltype(sp),eltype(domain(sp)),typeof(z)); JacobiZ{typeof(sp),T}(sp,T(z)))

Base.convert{T,S}(::Type{Operator{T}},J::JacobiZ{S}) = JacobiZ{S,T}(J.space,J.z)

domainspace(::JacobiZ) = ℓ⁰
rangespace(::JacobiZ) = ℓ⁰

#####
# recα/β/γ are given by
#       x p_{n-1} =γ_n p_{n-2} + α_n p_{n-1} +  p_n β_n
#####

function getindex{S,T}(J::JacobiZ{S,T},k::Integer,j::Integer)
    if j==k-1
        recγ(T,J.space,k)
    elseif j==k
        k == 1 ? one(T) : recα(T,J.space,k)-J.z
    elseif j==k+1 && k > 1
        recβ(T,J.space,k)
    else
        zero(T)
    end
end




#####
# Multiplication can be built from recurrence coefficients
#  multiplication by J=M[x] is Recurrence operator
#  and assume p_0 = 1, then we have
#
#   M[p_0] = 1
#   M[p_1] = (J/β_1 - α_1/β_1)*M[p_0]
#   M[p_k] = (J/β_k - α_k/β_k)*M[p_{k-1}] - γ_k/β_k*M[p_{k-2}]
#####


getindex{PS<:PolynomialSpace,T,C<:PolynomialSpace}(M::ConcreteMultiplication{C,PS,T},k::Integer,j::Integer) = M[k:k,j:j][1,1]


# Fast C = α.*A.*B .+ β.*C
function αA_dot_B_plus_βC!(α, A::AbstractVector, B::AbstractVector, β, C::AbstractVector)
    length(C) == length(A) == (n = length(B)) || throw(BoundsError())
    @fastmath @inbounds @simd for i=1:n
        C[i] = α*A[i]*B[i] + β*C[i]
    end
    C
end



# Fast implementation of C[:,:] = α*J*B+β*C where the bandediwth of B is
# specified by b, not by the parameters in B
function jac_gbmm!(α,J,B,β,C,b)
    if β ≠ 1
        scale!(β,C)
    end

    Jp = view(J, band(1))
    J0 = view(J, band(0))
    Jm = view(J, band(-1))
    n = size(J,1)

    for k=-1:b-1
        αA_dot_B_plus_βC!(α,view(view(B,band(b-k-1)),2:n-b+k+1),
                            view(Jp,1:n-b+k),1.0,view(C,band(b-k)))
        αA_dot_B_plus_βC!(α,view(view(B,band(k-b+1)),1:n-b+k),
                            view(Jm,b-k:n-1),1.0,view(C,band(k-b)))
        if k ≥ 0
            αA_dot_B_plus_βC!(α,view(B,band(b-k)),view(J0,1:n-b+k),1.0,view(C,band(b-k)))
            αA_dot_B_plus_βC!(α,view(B,band(k-b)),view(J0,b-k+1:n),1.0,view(C,band(k-b)))
            if k ≥ 1
                αA_dot_B_plus_βC!(α,view(B,band(b-k+1)),
                                    view(Jm,1:n-1-b+k),1.0,view(view(C,band(b-k)),2:n-b+k))
                αA_dot_B_plus_βC!(α,view(B,band(k-b-1)),
                                    view(Jp,b-k+1:n-1),1.0,view(view(C,band(k-b)),1:n-b+k-1))
            end
        end
    end

    αA_dot_B_plus_βC!(α,view(B,band(0)),J0,1.0,view(C,band(0)))
    αA_dot_B_plus_βC!(α,view(B,band(-1)),Jp,1.0,view(view(C,band(0)),1:n-1))
    αA_dot_B_plus_βC!(α,view(B,band(1)),Jm,1.0,view(view(C,band(0)),2:n))

    C
end





function Base.convert{PS<:PolynomialSpace,T,C<:PolynomialSpace}(::Type{BandedMatrix},
                                                                S::SubOperator{T,ConcreteMultiplication{C,PS,T},
                                                                               Tuple{UnitRange{Int},UnitRange{Int}}})
    M=parent(S)
    kr,jr=parentindexes(S)
    f=M.f
    a=f.coefficients
    sp=space(f)
    n=length(a)

    if n==0
        return bzeros(S)
    elseif n==1
        ret = bzeros(S)
        shft=kr[1]-jr[1]
        ret[band(shft)] = a[1]
        return ret::BandedMatrix{T}
    elseif n==2
        # we have U_x = [1 α-x; 0 β]
        # for e_1^⊤ U_x\a == a[1]*I-(α-J)*a[2]/β == (a[1]-α*a[2]/β)*I + J*a[2]/β
        # implying
        α,β=recα(T,sp,1),recβ(T,sp,1)
        ret=Operator{T}(ApproxFun.Recurrence(M.space))[kr,jr]::BandedMatrix{T}
        scale!(a[2]/β,ret)
        shft=kr[1]-jr[1]
        ret[band(shft)] += a[1]-α*a[2]/β
        return ret::BandedMatrix{T}
    end

    jkr=max(1,min(jr[1],kr[1])-(n-1)÷2):max(jr[end],kr[end])+(n-1)÷2

    #Multiplication is transpose
    J=Operator{T}(ApproxFun.Recurrence(M.space))[jkr,jkr]

    B=n-1  # final bandwidth

    # Clenshaw for operators
    Bk2 = bzeros(T,size(J,1),size(J,2),B,B)
    Bk2[band(0)] = a[n]/recβ(T,sp,n-1)
    α,β = recα(T,sp,n-1),recβ(T,sp,n-2)
    Bk1 = (-α/β)*Bk2
    Base.axpy!(a[n-1]/β,I,Bk1)
    jac_gbmm!(one(T)/β,J,Bk2,one(T),Bk1,0)
    b=1  # we keep track of bandwidths manually to reuse memory
    for k=n-2:-1:2
        α,β,γ=recα(T,sp,k),recβ(T,sp,k-1),recγ(T,sp,k+1)
        scale!(-γ/β,Bk2)
        Base.axpy!(a[k]/β,I,Bk2)
        jac_gbmm!(1/β,J,Bk1,one(T),Bk2,b)
        Base.axpy!(-α/β,Bk1,Bk2)
        Bk2,Bk1=Bk1,Bk2
        b+=1
    end
    α,γ=recα(T,sp,1),recγ(T,sp,2)
    scale!(-γ,Bk2)
    Base.axpy!(a[1],I,Bk2)
    jac_gbmm!(one(T),J,Bk1,one(T),Bk2,b)
    Base.axpy!(-α,Bk1,Bk2)

    # relationship between jkr and kr, jr
    kr2,jr2=kr-jkr[1]+1,jr-jkr[1]+1

    # TODO: reuse memory of Bk2, though profile suggests it's not too important
    BandedMatrix(view(Bk2,kr2,jr2))::BandedMatrix{T}
end




## All polynomial spaces can be converted provided spaces match

isconvertible(a::PolynomialSpace,b::PolynomialSpace) = domain(a)==domain(b)
union_rule{D}(a::PolynomialSpace{D},b::PolynomialSpace{D}) =
    domainscompatible(a,b)?(a<b?a:b):NoSpace()   # the union of two polys is always a poly


## General clenshaw
clenshaw(sp::PolynomialSpace,c::AbstractVector,x::AbstractArray) = clenshaw(c,x,
                                            ClenshawPlan(promote_type(eltype(c),eltype(x)),sp,length(c),length(x)))
clenshaw(sp::PolynomialSpace,c::AbstractMatrix,x::AbstractArray) = clenshaw(c,x,ClenshawPlan(promote_type(eltype(c),eltype(x)),sp,size(c,1),length(x)))
clenshaw(sp::PolynomialSpace,c::AbstractMatrix,x) = clenshaw(c,x,ClenshawPlan(promote_type(eltype(c),typeof(x)),sp,size(c,1),size(c,2)))

clenshaw!(sp::PolynomialSpace,c::AbstractVector,x::AbstractVector)=clenshaw!(c,x,ClenshawPlan(promote_type(eltype(c),eltype(x)),sp,length(x)))

function clenshaw(sp::PolynomialSpace,c::AbstractVector,x)
    N,T = length(c),promote_type(eltype(c),typeof(x))
    TT = eltype(T)
    if isempty(c)
        return zero(T)
    end

    bk1,bk2 = zero(T),zero(T)
    A,B,C = recA(TT,sp,N-1),recB(TT,sp,N-1),recC(TT,sp,N)
    for k = N:-1:2
        bk2, bk1 = bk1, muladd(muladd(A,x,B),bk1,muladd(-C,bk2,c[k])) # muladd(-C,bk2,muladd(muladd(A,x,B),bk1,c[k])) # (A*x+B)*bk1+c[k]-C*bk2
        A,B,C = recA(TT,sp,k-2),recB(TT,sp,k-2),recC(TT,sp,k-1)
    end
    muladd(muladd(A,x,B),bk1,muladd(-C,bk2,c[1])) # muladd(-C,bk2,muladd(muladd(A,x,B),bk1,c[1])) # (A*x+B)*bk1+c[1]-C*bk2
end



# evaluate polynomial
# indexing starts from 0
function forwardrecurrence{T}(::Type{T},S::Space,r::Range,x::Number)
    if isempty(r)
        return T[]
    end
    n=maximum(r)+1
    v=Vector{T}(n)  # x may be complex
    if n > 0
        v[1]=1
        if n > 1
            v[2] = (x-recα(T,S,1))*v[1]/recβ(T,S,1)

            @inbounds for k=2:n-1
                v[k+1]=((x-recα(T,S,k))*v[k] - recγ(T,S,k)*v[k-1])/recβ(T,S,k)
            end
        end
    end

    return v[r+1]
end


function Evaluation(S::PolynomialSpace,x,order)
    if order == 0
        ConcreteEvaluation(S,x,order)
    else
        # assume Derivative is available
        D = Derivative(S,order)
        EvaluationWrapper(S,x,order,Evaluation(rangespace(D),x)*D)
    end
end


function getindex{J<:PolynomialSpace}(op::ConcreteEvaluation{J,Bool},kr::Range)
    sp=op.space
    T=eltype(op)
    x=op.x

    forwardrecurrence(T,sp,kr-1,x?one(T):-one(T))
end


function getindex{J<:PolynomialSpace,TT<:Number}(op::ConcreteEvaluation{J,TT},kr::Range)
    sp=op.space
    T=eltype(op)
    x=op.x

    forwardrecurrence(T,sp,kr-1,x)
end
