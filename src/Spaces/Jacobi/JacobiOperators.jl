## Evaluation


getindex{J<:Jacobi}(op::ConcreteEvaluation{J,Bool},k::Integer) =
    op[k:k][1]

getindex{J<:Jacobi}(op::ConcreteEvaluation{J},k::Integer) =
    op[k:k][1]


function getindex{J<:Jacobi}(op::ConcreteEvaluation{J,Bool},kr::Range)
    @assert op.order <= 2
    sp=op.space
    a=sp.a;b=sp.b
    x=op.x

    if op.order == 0
        jacobip(kr-1,a,b,x?1.0:-1.0)
    elseif op.order == 1&& !x && b==0
        d=domain(op)
        @assert isa(d,Interval)
        Float64[tocanonicalD(d,d.a)*.5*(a+k)*(k-1)*(-1)^k for k=kr]
    elseif op.order == 1
        d=domain(op)
        @assert isa(d,Interval)
        if kr[1]==1
            0.5*tocanonicalD(d,d.a)*(a+b+kr).*[0.;jacobip(0:kr[end]-2,1+a,1+b,x?1.:-1.)]
        else
            0.5*tocanonicalD(d,d.a)*(a+b+kr).*jacobip(kr-1,1+a,1+b,x?1.:-1.)
        end
    elseif op.order == 2
        @assert !x && b==0
        @assert domain(op)==Interval()
        Float64[-.125*(a+k)*(a+k+1)*(k-2)*(k-1)*(-1)^k for k=kr]
    end
end
function getindex{J<:Jacobi}(op::ConcreteEvaluation{J,Float64},kr::Range)
    @assert op.order == 0
    jacobip(kr-1,op.space.a,op.space.b,tocanonical(domain(op),op.x))
end


## Derivative

Derivative(J::Jacobi,k::Integer)=k==1?ConcreteDerivative(J,1):DerivativeWrapper(TimesOperator(Derivative(Jacobi(J.a+1,J.b+1,J.domain),k-1),ConcreteDerivative(J,1)),k)



rangespace{J<:Jacobi}(D::ConcreteDerivative{J})=Jacobi(D.space.a+D.order,D.space.b+D.order,domain(D))
bandinds{J<:Jacobi}(D::ConcreteDerivative{J})=0,D.order

getindex{J<:Jacobi}(T::ConcreteDerivative{J},k::Integer,j::Integer) =
    j==k+1? (k+1+T.space.a+T.space.b)./complexlength(domain(T)) : zero(eltype(T))



function Derivative{T,DDD<:Interval}(S::WeightedJacobi{T,DDD})
    if S.α>0 && S.β>0 && S.α==S.space.b && S.β==S.space.a
        ConcreteDerivative(S,1)
    else
        jacobiweightDerivative(S)
    end
end

bandinds{T,DDD<:Interval}(D::ConcreteDerivative{WeightedJacobi{T,DDD}})=-1,0
rangespace{T,DDD<:Interval}(D::ConcreteDerivative{WeightedJacobi{T,DDD}})=WeightedJacobi(domainspace(D).α-1,domainspace(D).β-1,domain(D))


getindex{T,DDD<:Interval}(D::ConcreteDerivative{WeightedJacobi{T,DDD}},k::Integer,j::Integer) =
    j==k-1? -4(k-1)./complexlength(domain(D)) : zero(eltype(D))


## Integral

function Integral(J::Jacobi,k::Integer)
    if k > 1
        Q=Integral(J,1)
        IntegralWrapper(TimesOperator(Integral(rangespace(Q),k-1),Q),k)
    elseif J.a > 0 && J.b > 0   # we have a simple definition
        ConcreteIntegral(J,1)
    else   # convert and then integrate
        sp=Jacobi(J.a+1,J.b+1,domain(J))
        C=Conversion(J,sp)
        Q=Integral(sp,1)
        IntegralWrapper(TimesOperator(Q,C),1)
    end
end


rangespace{J<:Jacobi}(D::ConcreteIntegral{J})=Jacobi(D.space.a-D.order,D.space.b-D.order,domain(D))
bandinds{J<:Jacobi}(D::ConcreteIntegral{J})=-D.order,0

function getindex{J<:Jacobi}(T::ConcreteIntegral{J},k::Integer,j::Integer)
    @assert T.order==1
    if k≥2 && j==k-1
        complexlength(domain(T))./(k+T.space.a+T.space.b-2)
    else
        zero(eltype(J))
    end
end


## Volterra Integral operator

Volterra(d::Interval) = Volterra(Legendre(d))
function Volterra(S::Jacobi,order::Integer)
    @assert S.a == S.b == 0.0
    @assert order==1
    ConcreteVolterra(S,order)
end

rangespace{J<:Jacobi}(V::ConcreteVolterra{J})=Jacobi(0.0,-1.0,domain(V))
bandinds{J<:Jacobi}(V::ConcreteVolterra{J})=-1,0

function getindex{J<:Jacobi}(V::ConcreteVolterra{J},k::Integer,j::Integer)
    d=domain(V)
    C = 0.5(d.b-d.a)
    if k≥2
        if j==k-1
            C/(k-1.5)
        elseif j==k
            -C/(k-0.5)
        else
            zero(eltype(V))
        end
    else
        zero(eltype(V))
    end
end


## Conversion
# We can only increment by a or b by one, so the following
# multiplies conversion operators to handle otherwise

function Conversion(L::Jacobi,M::Jacobi)
    @assert (isapprox(M.b,L.b)||M.b>=L.b) && (isapprox(M.a,L.a)||M.a>=L.a)
    dm=domain(M)
    D=typeof(dm)
    if isapprox(M.a,L.a) && isapprox(M.b,L.b)
        SpaceOperator(IdentityOperator(),L,M)
    elseif (isapprox(M.b,L.b+1) && isapprox(M.a,L.a)) || (isapprox(M.b,L.b) && isapprox(M.a,L.a+1))
        ConcreteConversion(L,M)
    elseif M.b > L.b+1
        ConversionWrapper(TimesOperator(Conversion(Jacobi(M.a,M.b-1,dm),M),Conversion(L,Jacobi(M.a,M.b-1,dm))))
    else  #if M.a >= L.a+1
        ConversionWrapper(TimesOperator(Conversion(Jacobi(M.a-1,M.b,dm),M),Conversion(L,Jacobi(M.a-1,M.b,dm))))
    end
end

bandinds{J<:Jacobi}(C::ConcreteConversion{J,J})=(0,1)



function Base.getindex{J<:Jacobi,T}(C::ConcreteConversion{J,J,T},k::Integer,j::Integer)
    L=C.domainspace
    if L.b+1==C.rangespace.b
        if j==k
            k==1?T(1):T((L.a+L.b+k)/(L.a+L.b+2k-1))
        elseif j==k+1
            T((L.a+k)./(L.a+L.b+2k+1))
        else
            zero(T)
        end
    elseif L.a+1==C.rangespace.a
        if j==k
            k==1?T(1):T((L.a+L.b+k)/(L.a+L.b+2k-1))
        elseif j==k+1
            T(-(L.b+k)./(L.a+L.b+2k+1))
        else
            zero(T)
        end
    else
        error("Not implemented")
    end
end




# return the space that has banded Conversion to the other
function conversion_rule(A::Jacobi,B::Jacobi)
    if !isapproxinteger(A.a-B.a) || !isapproxinteger(A.b-B.b)
        NoSpace()
    else
        Jacobi(min(A.a,B.a),min(A.b,B.b),domain(A))
    end
end



## Ultraspherical Conversion

# Assume m is compatible

function Conversion(A::PolynomialSpace,B::Jacobi)
    J = Jacobi(A)
    J == B ? ConcreteConversion(A,B) :
             ConversionWrapper(TimesOperator(Conversion(J,B),Conversion(A,J)))
end

function Conversion(A::Jacobi,B::PolynomialSpace)
    J = Jacobi(B)
    J == A ? ConcreteConversion(A,B) :
             ConversionWrapper(TimesOperator(Conversion(J,B),Conversion(A,J)))
end



bandinds{US<:Chebyshev,J<:Jacobi}(C::ConcreteConversion{US,J}) = 0,0
bandinds{US<:Chebyshev,J<:Jacobi}(C::ConcreteConversion{J,US}) = 0,0


bandinds{US<:Ultraspherical,J<:Jacobi}(C::ConcreteConversion{US,J}) = 0,0
bandinds{US<:Ultraspherical,J<:Jacobi}(C::ConcreteConversion{J,US}) = 0,0


function getindex{J<:Jacobi,CC<:Chebyshev,T}(C::ConcreteConversion{CC,J,T},k::Integer,j::Integer)
    if j==k
        one(T)/jacobip(k-1,-one(T)/2,-one(T)/2,one(T))
    else
        zero(T)
    end
end

function getindex{J<:Jacobi,CC<:Chebyshev,T}(C::ConcreteConversion{J,CC,T},k::Integer,j::Integer)
    if j==k
        jacobip(k-1,-one(T)/2,-one(T)/2,one(T))
    else
        zero(T)
    end
end

function getindex{US<:Ultraspherical,J<:Jacobi,T}(C::ConcreteConversion{US,J,T},k::Integer,j::Integer)
    if j==k
        S=rangespace(C)
        jp=jacobip(k-1,S.a,S.b,one(T))
        um=Evaluation(setcanonicaldomain(domainspace(C)),one(T))[k]
        um/jp
    else
        zero(T)
    end
end

function getindex{US<:Ultraspherical,J<:Jacobi,T}(C::ConcreteConversion{J,US,T},k::Integer,j::Integer)
    if j==k
        S=domainspace(C)
        jp=jacobip(k-1,S.a,S.b,one(T))
        um=Evaluation(setcanonicaldomain(rangespace(C)),one(T))[k]
        jp/um
    else
        zero(T)
    end
end





union_rule(A::Jacobi,B::Jacobi)=conversion_type(A,B)
function maxspace_rule(A::Jacobi,B::Jacobi)
    if !isapproxinteger(A.a-B.a) || !isapproxinteger(A.b-B.b)
        NoSpace()
    else
        Jacobi(max(A.a,B.a),max(A.b,B.b),domain(A))
    end
end


for (OPrule,OP) in ((:conversion_rule,:conversion_type),(:maxspace_rule,:maxspace),(:union_rule,:(Base.union)))
    @eval begin
        function $OPrule(A::Chebyshev,B::Jacobi)
            if !isapproxinteger(-0.5-B.a) || !isapproxinteger(-0.5-B.b)
                NoSpace()
            elseif isapprox(B.a,-0.5)&&isapprox(B.b,-0.5)
                # the spaces are the same
                A
            else
                $OP(Jacobi(A),B)
            end
        end
        function $OPrule{m}(A::Ultraspherical{m},B::Jacobi)
            if !isapproxinteger(m-0.5-B.a) || !isapproxinteger(m-0.5-B.b)
                NoSpace()
            elseif isapprox(B.a,m-0.5)&&isapprox(B.b,m-0.5)
                # the spaces are the same
                A
            else
                $OP(Jacobi(A),B)
            end
        end
    end
end

hasconversion(a::Jacobi,b::Chebyshev) = hasconversion(a,Jacobi(b))
hasconversion(a::Chebyshev,b::Jacobi) = hasconversion(Jacobi(a),b)

hasconversion(a::Jacobi,b::Ultraspherical) = hasconversion(a,Jacobi(b))
hasconversion(a::Ultraspherical,b::Jacobi) = hasconversion(Jacobi(a),b)




## Special Multiplication
# special multiplication operators exist when multiplying by
# (1+x) or (1-x) by _decreasing_ the parameter.  Thus the


## <: IntervalDomain avoids a julia bug
function Multiplication{C<:ConstantSpace,DD<:IntervalDomain}(f::Fun{JacobiWeight{C,DD}},S::Jacobi)
    # this implements (1+x)*P and (1-x)*P special case
    # see DLMF (18.9.6)
    d=domain(f)
    if ((space(f).α==1 && space(f).β==0 && S.b >0) ||
                        (space(f).α==0 && space(f).β==1 && S.a >0))
        ConcreteMultiplication(f,S)
    elseif isapproxinteger(space(f).α) && space(f).α ≥ 1 && S.b >0
        # decrement α and multiply again
        M=Multiplication(f.coefficients[1]*jacobiweight(1.,0.,d),S)
        MultiplicationWrapper(f,Multiplication(jacobiweight(space(f).α-1,space(f).β,d),rangespace(M))*M)
    elseif isapproxinteger(space(f).β) && space(f).β ≥ 1 && S.a >0
        # decrement β and multiply again
        M=Multiplication(f.coefficients[1]*jacobiweight(0.,1.,d),S)
        MultiplicationWrapper(f,Multiplication(jacobiweight(space(f).α,space(f).β-1,d),rangespace(M))*M)
    else
# default JacobiWeight
        M=Multiplication(Fun(f.coefficients,space(f).space),S)
        rsp=JacobiWeight(space(f).α,space(f).β,rangespace(M))
        MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
    end
end

Multiplication{C<:ConstantSpace,DD<:IntervalDomain}(f::Fun{JacobiWeight{C,DD}},S::Union{Ultraspherical,Chebyshev}) =
    MultiplicationWrapper(f,Multiplication(f,Jacobi(S))*Conversion(S,Jacobi(S)))

function rangespace{J<:Jacobi,C<:ConstantSpace,DD<:IntervalDomain}(M::ConcreteMultiplication{JacobiWeight{C,DD},J})
    S=domainspace(M)
    if space(M.f).α==1
        # multiply by (1+x)
        Jacobi(S.a,S.b-1,domain(S))
    elseif space(M.f).β == 1
        # multiply by (1-x)
        Jacobi(S.a-1,S.b,domain(S))
    else
        error("Not implemented")
    end
end

bandinds{J<:Jacobi,C<:ConstantSpace,DD<:IntervalDomain}(::ConcreteMultiplication{JacobiWeight{C,DD},J})=-1,0

function getindex{J<:Jacobi,C<:ConstantSpace,DD<:IntervalDomain}(M::ConcreteMultiplication{JacobiWeight{C,DD},J},k::Integer,j::Integer)
    @assert ncoefficients(M.f)==1
    a,b=domainspace(M).a,domainspace(M).b
    c=M.f.coefficients[1]
    if space(M.f).α==1
        @assert space(M.f).β==0
        # multiply by (1+x)
        if j==k
            c*2(k+b-1)/(2k+a+b-1)
        elseif k > 1 && j==k-1
            c*(2k-2)/(2k+a+b-3)
        else
            zero(eltype(M))
        end
    elseif space(M.f).β == 1
        @assert space(M.f).α==0
        # multiply by (1-x)
        if j==k
            c*2(k+a-1)/(2k+a+b-1)
        elseif k > 1 && j==k-1
            -c*(2k-2)/(2k+a+b-3)
        else
            zero(eltype(M))
        end
    else
        error("Not implemented")
    end
end


# We can exploit the special multiplication to construct a Conversion


for FUNC in (:maxspace_rule,:union_rule,:hasconversion)
    @eval function $FUNC{T,DD<:Interval}(A::WeightedJacobi{T,DD},B::Jacobi)
        if A.α==A.β+1 && A.space.b>0
            $FUNC(Jacobi(A.space.a,A.space.b-1,domain(A)),B)
        elseif A.β==A.α+1 && A.space.a>0
            $FUNC(Jacobi(A.space.a-1,A.space.b,domain(A)),B)
        else
            $FUNC(A,JacobiWeight(0.,0.,B))
        end
    end
end





# represents [b+(1+z)*d/dz] (false) or [a-(1-z)*d/dz] (true)
immutable JacobiSD{T} <:Operator{T}
    lr::Bool
    S::Jacobi
end

JacobiSD(lr,S)=JacobiSD{Float64}(lr,S)

Base.convert{T}(::Type{Operator{T}},SD::JacobiSD)=JacobiSD{T}(SD.lr,SD.S)

domain(op::JacobiSD)=domain(op.S)
domainspace(op::JacobiSD)=op.S
rangespace(op::JacobiSD)=op.lr?Jacobi(op.S.a-1,op.S.b+1,domain(op.S)):Jacobi(op.S.a+1,op.S.b-1,domain(op.S))
bandinds(::JacobiSD)=0,0

function getindex(op::JacobiSD,A,k::Integer,j::Integer)
    m=op.lr?op.S.a:op.S.b
    if k==j
        k+m-1
    else
        zero(eltype(op))
    end
end
