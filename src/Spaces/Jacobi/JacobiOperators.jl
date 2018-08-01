## Evaluation


function Evaluation(S::Jacobi,x::Union{typeof(first),typeof(last)},order)
    if order ≤ 2
        ConcreteEvaluation(S,x,order)
    else
        # assume Derivative is available
        D = Derivative(S,order)
        EvaluationWrapper(S,x,order,Evaluation(rangespace(D),x)*D)
    end
end

function Evaluation(S::Jacobi,x,order)
    if order == 0
        ConcreteEvaluation(S,x,order)
    else
        # assume Derivative is available
        D = Derivative(S,order)
        EvaluationWrapper(S,x,order,Evaluation(rangespace(D),x)*D)
    end
end

# avoid ambiguity
for OP in (:first,:last)
    @eval getindex(op::ConcreteEvaluation{<:Jacobi,typeof($OP)},k::Integer) =
        op[k:k][1]
end

getindex(op::ConcreteEvaluation{<:Jacobi},k::Integer) = op[k:k][1]


function getindex(op::ConcreteEvaluation{<:Jacobi,typeof(first)},kr::AbstractRange)
    @assert op.order <= 2
    sp=op.space
    T=eltype(op)
    RT=real(T)
    a=convert(RT,sp.a);b=convert(RT,sp.b)

    if op.order == 0
        jacobip(T,kr-1,a,b,-one(T))
    elseif op.order == 1&&  b==0
        d=domain(op)
        @assert isa(d,Segment)
        T[tocanonicalD(d,leftendpoint(d))/2*(a+k)*(k-1)*(-1)^k for k=kr]
    elseif op.order == 1
        d=domain(op)
        @assert isa(d,Segment)
        if kr[1]==1 && kr[end] ≥ 2
            tocanonicalD(d,leftendpoint(d))*(a+b+kr).*T[zero(T);jacobip(T,0:kr[end]-2,1+a,1+b,-one(T))]/2
        elseif kr[1]==1  # kr[end] ≤ 1
            zeros(T,length(kr))
        else
            tocanonicalD(d,leftendpoint(d))*(a+b+kr).*jacobip(T,kr-1,1+a,1+b,-one(T))/2
        end
    elseif op.order == 2
        @assert b==0
        @assert domain(op) == Segment()
        T[-0.125*(a+k)*(a+k+1)*(k-2)*(k-1)*(-1)^k for k=kr]
    else
        error("Not implemented")
    end
end

function getindex(op::ConcreteEvaluation{<:Jacobi,typeof(last)},kr::AbstractRange)
    @assert op.order <= 2
    sp=op.space
    T=eltype(op)
    RT=real(T)
    a=convert(RT,sp.a);b=convert(RT,sp.b)


    if op.order == 0
        jacobip(T,kr.-1,a,b,one(T))
    elseif op.order == 1
        d=domain(op)
        @assert isa(d,Segment)
        if kr[1]==1 && kr[end] ≥ 2
            tocanonicalD(d,leftendpoint(d))*((a+b).+kr).*T[zero(T);jacobip(T,0:kr[end]-2,1+a,1+b,one(T))]/2
        elseif kr[1]==1  # kr[end] ≤ 1
            zeros(T,length(kr))
        else
            tocanonicalD(d,leftendpoint(d))*((a+b).+kr).*jacobip(T,kr.-1,1+a,1+b,one(T))/2
        end
    else
        error("Not implemented")
    end
end


function getindex(op::ConcreteEvaluation{<:Jacobi,<:Number},kr::AbstractRange)
    @assert op.order == 0
    jacobip(eltype(op),kr-1,op.space.a,op.space.b,tocanonical(domain(op),op.x))
end


## Derivative

Derivative(J::Jacobi,k::Integer)=k==1 ? ConcreteDerivative(J,1) : DerivativeWrapper(TimesOperator(Derivative(Jacobi(J.b+1,J.a+1,J.domain),k-1),ConcreteDerivative(J,1)),k)



rangespace(D::ConcreteDerivative{J}) where {J<:Jacobi}=Jacobi(D.space.b+D.order,D.space.a+D.order,domain(D))
bandinds(D::ConcreteDerivative{J}) where {J<:Jacobi}=0,D.order

getindex(T::ConcreteDerivative{J},k::Integer,j::Integer) where {J<:Jacobi} =
    j==k+1 ? eltype(T)((k+1+T.space.a+T.space.b)/complexlength(domain(T))) : zero(eltype(T))



function Derivative(S::WeightedJacobi{DDD,RR}) where {DDD<:IntervalOrSegment,RR}
    if S.β>0 && S.β>0 && S.β==S.space.b && S.α==S.space.a
        ConcreteDerivative(S,1)
    else
        jacobiweightDerivative(S)
    end
end

bandinds(D::ConcreteDerivative{WeightedJacobi{DDD,RR}}) where {DDD<:IntervalOrSegment,RR} = -1,0
rangespace(D::ConcreteDerivative{WeightedJacobi{DDD,RR}}) where {DDD<:IntervalOrSegment,RR} =
    WeightedJacobi(domainspace(D).β-1,domainspace(D).α-1,domain(D))

getindex(D::ConcreteDerivative{WeightedJacobi{DDD,RR}},k::Integer,j::Integer) where {DDD<:IntervalOrSegment,RR} =
    j==k-1 ? eltype(D)(-4(k-1)./complexlength(domain(D))) : zero(eltype(D))



## Integral

function Integral(J::Jacobi,k::Integer)
    if k > 1
        Q=Integral(J,1)
        IntegralWrapper(TimesOperator(Integral(rangespace(Q),k-1),Q),k)
    elseif J.a > 0 && J.b > 0   # we have a simple definition
        ConcreteIntegral(J,1)
    else   # convert and then integrate
        sp=Jacobi(J.b+1,J.a+1,domain(J))
        C=Conversion(J,sp)
        Q=Integral(sp,1)
        IntegralWrapper(TimesOperator(Q,C),1)
    end
end


rangespace(D::ConcreteIntegral{J}) where {J<:Jacobi}=Jacobi(D.space.b-D.order,D.space.a-D.order,domain(D))
bandinds(D::ConcreteIntegral{J}) where {J<:Jacobi}=-D.order,0

function getindex(T::ConcreteIntegral{J},k::Integer,j::Integer) where J<:Jacobi
    @assert T.order==1
    if k≥2 && j==k-1
        complexlength(domain(T))./(k+T.space.a+T.space.b-2)
    else
        zero(eltype(T))
    end
end


## Volterra Integral operator

Volterra(d::Segment) = Volterra(Legendre(d))
function Volterra(S::Jacobi,order::Integer)
    @assert S.a == S.b == 0.0
    @assert order==1
    ConcreteVolterra(S,order)
end

rangespace(V::ConcreteVolterra{J}) where {J<:Jacobi}=Jacobi(-1.0,0.0,domain(V))
bandinds(V::ConcreteVolterra{J}) where {J<:Jacobi}=-1,0

function getindex(V::ConcreteVolterra{J},k::Integer,j::Integer) where J<:Jacobi
    d=domain(V)
    C = complexlength(d)/2
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


for (Func,Len,Sum) in ((:DefiniteIntegral,:complexlength,:sum),(:DefiniteLineIntegral,:arclength,:linesum))
    ConcFunc = Meta.parse("Concrete"*string(Func))

    @eval begin
        $Func(S::Jacobi{<:Segment}) = $ConcFunc(S)

        function getindex(Σ::$ConcFunc{Jacobi{D,R},T},k::Integer) where {D<:Segment,R,T}
            dsp = domainspace(Σ)

            if dsp.b == dsp.a == 0
                # TODO: copy and paste
                k == 1 ? convert(T,$Sum(Fun(dsp,[one(T)]))) : zero(T)
            else
                convert(T,$Sum(Fun(dsp,[zeros(T,k-1);1])))
            end
        end


        function getindex(Σ::$ConcFunc{JacobiWeight{Jacobi{D,R},D,R,TT},T},k::Integer) where {D<:Segment,R,T,TT}
            dsp = domainspace(Σ)

            if dsp.β == dsp.space.b && dsp.α == dsp.space.a
                # TODO: copy and paste
                k == 1 ? convert(T,$Sum(Fun(dsp,[one(T)]))) : zero(T)
            else
                convert(T,$Sum(Fun(dsp,[zeros(T,k-1);1])))
            end
        end

        function bandinds(Σ::$ConcFunc{JacobiWeight{Jacobi{D,R},D,R,TT}}) where {D<:Segment,R,TT}
            β,α = domainspace(Σ).β,domainspace(Σ).α
            if domainspace(Σ).β == domainspace(Σ).space.b && domainspace(Σ).α == domainspace(Σ).space.a
                0,0  # first entry
            else
                0,∞
            end
        end

        function bandinds(Σ::$ConcFunc{Jacobi{D,R}}) where {D<:Segment,R}
            if domainspace(Σ).b == domainspace(Σ).a == 0
                0,0  # first entry
            else
                0,∞
            end
        end
    end
end



## Conversion
# We can only increment by a or b by one, so the following
# multiplies conversion operators to handle otherwise

function Conversion(L::Jacobi,M::Jacobi)
    if isapproxinteger(L.a-M.a) && isapproxinteger(L.b-M.b)
        dm=domain(M)
        D=typeof(dm)
        if isapprox(M.a,L.a) && isapprox(M.b,L.b)
            ConversionWrapper(eye(L))
        elseif (isapprox(M.b,L.b+1) && isapprox(M.a,L.a)) || (isapprox(M.b,L.b) && isapprox(M.a,L.a+1))
            ConcreteConversion(L,M)
        elseif M.b > L.b+1
            ConversionWrapper(TimesOperator(Conversion(Jacobi(M.b-1,M.a,dm),M),Conversion(L,Jacobi(M.b-1,M.a,dm))))
        else  #if M.a >= L.a+1
            ConversionWrapper(TimesOperator(Conversion(Jacobi(M.b,M.a-1,dm),M),Conversion(L,Jacobi(M.b,M.a-1,dm))))
        end
    elseif L.a ≈ L.b ≈ 0. && M.a ≈ M.b ≈ 0.5
        Conversion(L,Ultraspherical(L),Ultraspherical(M),M)
    elseif L.a ≈ L.b ≈ 0. && M.a ≈ M.b ≈ -0.5
        Conversion(L,Ultraspherical(L),Chebyshev(M),M)
    elseif L.a ≈ L.b ≈ -0.5 && M.a ≈ M.b ≈ 0.5
        Conversion(L,Chebyshev(L),Ultraspherical(M),M)
    else # L.a - M.a ≈ L.b - M.b
        error("Implement for $L → $M")
    end
end

bandinds(C::ConcreteConversion{J1,J2}) where {J1<:Jacobi,J2<:Jacobi}=(0,1)



function Base.getindex(C::ConcreteConversion{J1,J2,T},k::Integer,j::Integer) where {J1<:Jacobi,J2<:Jacobi,T}
    L=C.domainspace
    if L.b+1==C.rangespace.b
        if j==k
            k==1 ? convert(T,1) : convert(T,(L.a+L.b+k)/(L.a+L.b+2k-1))
        elseif j==k+1
            convert(T,(L.a+k)./(L.a+L.b+2k+1))
        else
            zero(T)
        end
    elseif L.a+1==C.rangespace.a
        if j==k
            k==1 ? convert(T,1) : convert(T,(L.a+L.b+k)/(L.a+L.b+2k-1))
        elseif j==k+1
            convert(T,-(L.b+k)./(L.a+L.b+2k+1))
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
        Jacobi(min(A.b,B.b),min(A.a,B.a),domain(A))
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

function Conversion(A::Jacobi,B::Chebyshev)
    if A.a == A.b == -0.5
        ConcreteConversion(A,B)
    elseif A.a == A.b == 0
        ConversionWrapper(
            SpaceOperator(
                Conversion(Ultraspherical(1//2),B),
                A,B))
    elseif A.a == A.b
        US = Ultraspherical(A)
        ConversionWrapper(Conversion(US,B)*Conversion(A,US))
    else
        J = Jacobi(B)
        Conversion(J,B)*Conversion(A,J)
    end
end

function Conversion(A::Chebyshev,B::Jacobi)
    if B.a == B.b == -0.5
        ConcreteConversion(A,B)
    elseif B.a == B.b == 0
        ConversionWrapper(
            SpaceOperator(
                Conversion(A,Ultraspherical(1//2,domain(B))),
                A,B))
    elseif B.a == B.b
        US = Ultraspherical(B)
        ConversionWrapper(Conversion(US,B)*Conversion(A,US))
    else
        J = Jacobi(A)
        Conversion(J,B)*Conversion(A,J)
    end
end


function Conversion(A::Jacobi,B::Ultraspherical)
    if A.a == A.b == -0.5
        ConversionWrapper(Conversion(Chebyshev(domain(A)),B)*Conversion(A,Chebyshev(domain(A))))
    elseif A.a == A.b == order(B)-0.5
        ConcreteConversion(A,B)
    elseif A.a == A.b == 0
        ConversionWrapper(
            SpaceOperator(
                Conversion(Ultraspherical(1//2),B),
                A,B))
    elseif A.a == A.b
        US = Ultraspherical(A)
        ConversionWrapper(Conversion(US,B)*Conversion(A,US))
    else
        J = Jacobi(B)
        Conversion(J,B)*Conversion(A,J)
    end
end

function Conversion(A::Ultraspherical,B::Jacobi)
    if B.a == B.b == -0.5
        ConversionWrapper(Conversion(Chebyshev(domain(A)),B)*Conversion(A,Chebyshev(domain(A))))
    elseif B.a == B.b == order(A)-0.5
        ConcreteConversion(A,B)
    elseif B.a == B.b == 0
        ConversionWrapper(
            SpaceOperator(
                Conversion(A,Ultraspherical(1//2,domain(B))),
                A,B))
    elseif B.a == B.b
        US = Ultraspherical(B)
        ConversionWrapper(Conversion(US,B)*Conversion(A,US))
    else
        J = Jacobi(A)
        Conversion(J,B)*Conversion(A,J)
    end
end




bandinds(C::ConcreteConversion{US,J}) where {US<:Chebyshev,J<:Jacobi} = 0,0
bandinds(C::ConcreteConversion{J,US}) where {US<:Chebyshev,J<:Jacobi} = 0,0


bandinds(C::ConcreteConversion{US,J}) where {US<:Ultraspherical,J<:Jacobi} = 0,0
bandinds(C::ConcreteConversion{J,US}) where {US<:Ultraspherical,J<:Jacobi} = 0,0

#TODO: Figure out how to unify these definitions
function getindex(C::ConcreteConversion{CC,J,T},k::Integer,j::Integer) where {J<:Jacobi,CC<:Chebyshev,T}
    if j==k
        one(T)/jacobip(T,k-1,-one(T)/2,-one(T)/2,one(T))
    else
        zero(T)
    end
end

function BandedMatrix(S::SubOperator{T,ConcreteConversion{CC,J,T},Tuple{UnitRange{Int},UnitRange{Int}}}) where {J<:Jacobi,CC<:Chebyshev,T}
    ret=BandedMatrix(Zeros, S)
    kr,jr = parentindices(S)
    k=(kr ∩ jr)

    vals = one(T)./jacobip(T,k .- 1,-one(T)/2,-one(T)/2,one(T))

    ret[band(bandshift(S))] = vals
    ret
end


function getindex(C::ConcreteConversion{J,CC,T},k::Integer,j::Integer) where {J<:Jacobi,CC<:Chebyshev,T}
    if j==k
        jacobip(T,k-1,-one(T)/2,-one(T)/2,one(T))
    else
        zero(T)
    end
end

function BandedMatrix(S::SubOperator{T,ConcreteConversion{J,CC,T},Tuple{UnitRange{Int},UnitRange{Int}}}) where {J<:Jacobi,CC<:Chebyshev,T}
    ret=BandedMatrix(Zeros, S)
    kr,jr = parentindices(S)
    k=(kr ∩ jr)

    vals = jacobip(T,k.-1,-one(T)/2,-one(T)/2,one(T))

    ret[band(bandshift(S))] = vals
    ret
end


function getindex(C::ConcreteConversion{US,J,T},k::Integer,j::Integer) where {US<:Ultraspherical,J<:Jacobi,T}
    if j==k
        S=rangespace(C)
        jp=jacobip(T,k-1,S.a,S.b,one(T))
        um=Evaluation(setcanonicaldomain(domainspace(C)),last,0)[k]
        um/jp::T
    else
        zero(T)
    end
end

function BandedMatrix(S::SubOperator{T,ConcreteConversion{US,J,T},Tuple{UnitRange{Int},UnitRange{Int}}}) where {US<:Ultraspherical,J<:Jacobi,T}
    ret=BandedMatrix(Zeros, S)
    kr,jr = parentindices(S)
    k=(kr ∩ jr)

    sp=rangespace(parent(S))
    jp=jacobip(T,k.-1,sp.a,sp.b,one(T))
    um=Evaluation(T,setcanonicaldomain(domainspace(parent(S))),last,0)[k]
    vals = um./jp

    ret[band(bandshift(S))] = vals
    ret
end



function getindex(C::ConcreteConversion{J,US,T},k::Integer,j::Integer) where {US<:Ultraspherical,J<:Jacobi,T}
    if j==k
        S=domainspace(C)
        jp=jacobip(T,k-1,S.a,S.b,one(T))
        um=Evaluation(T,setcanonicaldomain(rangespace(C)),last,0)[k]
        jp/um::T
    else
        zero(T)
    end
end

function BandedMatrix(S::SubOperator{T,ConcreteConversion{J,US,T},Tuple{UnitRange{Int},UnitRange{Int}}}) where {US<:Ultraspherical,J<:Jacobi,T}
    ret=BandedMatrix(Zeros, S)
    kr,jr = parentindices(S)
    k=(kr ∩ jr)

    sp=domainspace(parent(S))
    jp=jacobip(T,k.-1,sp.a,sp.b,one(T))
    um=Evaluation(T,setcanonicaldomain(rangespace(parent(S))),last,0)[k]
    vals = jp./um

    ret[band(bandshift(S))] = vals
    ret
end






function union_rule(A::Jacobi,B::Jacobi)
    if domainscompatible(A,B)
        Jacobi(min(A.b,B.b),min(A.a,B.a),domain(A))
    else
        NoSpace()
    end
end

function maxspace_rule(A::Jacobi,B::Jacobi)
    if isapproxinteger(A.a-B.a) && isapproxinteger(A.b-B.b)
        Jacobi(max(A.b,B.b),max(A.a,B.a),domain(A))
    else
        NoSpace()
    end
end


function union_rule(A::Chebyshev,B::Jacobi)
    if isapprox(B.a,-0.5) && isapprox(B.b,-0.5)
        # the spaces are the same
        A
    else
        union(Jacobi(A),B)
    end
end
function union_rule(A::Ultraspherical,B::Jacobi)
    m=order(A)
    if isapprox(B.a,m-0.5) && isapprox(B.b,m-0.5)
        # the spaces are the same
        A
    else
        union(Jacobi(A),B)
    end
end

for (OPrule,OP) in ((:conversion_rule,:conversion_type), (:maxspace_rule,:maxspace))
    @eval begin
        function $OPrule(A::Chebyshev,B::Jacobi)
            if B.a ≈ -0.5 && B.b ≈ -0.5
                # the spaces are the same
                A
            elseif isapproxinteger(B.a+0.5) && isapproxinteger(B.b+0.5)
                $OP(Jacobi(A),B)
            else
                NoSpace()
            end
        end
        function $OPrule(A::Ultraspherical,B::Jacobi)
            m = order(A)
            if B.a ≈ m-0.5 && B.b ≈ m-0.5
                # the spaces are the same
                A
            elseif isapproxinteger(B.a+0.5) && isapproxinteger(B.b+0.5)
                $OP(Jacobi(A),B)
            else
                NoSpace()
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


## <: IntervalOrSegment avoids a julia bug
function Multiplication(f::Fun{JacobiWeight{C,DD,RR,TT}}, S::Jacobi) where {C<:ConstantSpace,DD<:IntervalOrSegment,RR,TT}
    # this implements (1+x)*P and (1-x)*P special case
    # see DLMF (18.9.6)
    d=domain(f)
    if ((space(f).β==1 && space(f).α==0 && S.b >0) ||
                        (space(f).β==0 && space(f).α==1 && S.a >0))
        ConcreteMultiplication(f,S)
    elseif isapproxinteger(space(f).β) && space(f).β ≥ 1 && S.b >0
        # decrement β and multiply again
        M=Multiplication(f.coefficients[1]*jacobiweight(1.,0.,d),S)
        MultiplicationWrapper(f,Multiplication(jacobiweight(space(f).β-1,space(f).α,d),rangespace(M))*M)
    elseif isapproxinteger(space(f).α) && space(f).α ≥ 1 && S.a >0
        # decrement α and multiply again
        M=Multiplication(f.coefficients[1]*jacobiweight(0.,1.,d),S)
        MultiplicationWrapper(f,Multiplication(jacobiweight(space(f).β,space(f).α-1,d),rangespace(M))*M)
    else
# default JacobiWeight
        M=Multiplication(Fun(space(f).space,f.coefficients),S)
        rsp=JacobiWeight(space(f).β,space(f).α,rangespace(M))
        MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
    end
end

Multiplication(f::Fun{JacobiWeight{C,DD,RR,TT}},
               S::Union{Ultraspherical,Chebyshev}) where {C<:ConstantSpace,DD<:IntervalOrSegment,RR,TT} =
    MultiplicationWrapper(f,Multiplication(f,Jacobi(S))*Conversion(S,Jacobi(S)))

function rangespace(M::ConcreteMultiplication{JacobiWeight{C,DD,RR,TT},J}) where {J<:Jacobi,C<:ConstantSpace,DD<:IntervalOrSegment,RR,TT}
    S=domainspace(M)
    if space(M.f).β==1
        # multiply by (1+x)
        Jacobi(S.b-1,S.a,domain(S))
    elseif space(M.f).α == 1
        # multiply by (1-x)
        Jacobi(S.b,S.a-1,domain(S))
    else
        error("Not implemented")
    end
end

bandinds(::ConcreteMultiplication{JacobiWeight{C,DD,RR,TT},J}) where {J<:Jacobi,C<:ConstantSpace,DD<:IntervalOrSegment,RR,TT} = -1,0

function getindex(M::ConcreteMultiplication{JacobiWeight{C,DD,RR,TT},J},k::Integer,j::Integer) where {J<:Jacobi,C<:ConstantSpace,DD<:IntervalOrSegment,RR,TT}
    @assert ncoefficients(M.f)==1
    a,b=domainspace(M).a,domainspace(M).b
    c=M.f.coefficients[1]
    if space(M.f).β==1
        @assert space(M.f).α==0
        # multiply by (1+x)
        if j==k
            c*2(k+b-1)/(2k+a+b-1)
        elseif k > 1 && j==k-1
            c*(2k-2)/(2k+a+b-3)
        else
            zero(eltype(M))
        end
    elseif space(M.f).α == 1
        @assert space(M.f).β==0
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
    @eval function $FUNC(A::WeightedJacobi{DD},B::Jacobi) where DD<:IntervalOrSegment
        if A.β==A.α+1 && A.space.b>0
            $FUNC(Jacobi(A.space.b-1,A.space.a,domain(A)),B)
        elseif A.α==A.β+1 && A.space.a>0
            $FUNC(Jacobi(A.space.b,A.space.a-1,domain(A)),B)
        else
            $FUNC(A,JacobiWeight(0.,0.,B))
        end
    end
end





# represents [b+(1+z)*d/dz] (false) or [a-(1-z)*d/dz] (true)
struct JacobiSD{T} <:Operator{T}
    lr::Bool
    S::Jacobi
end

JacobiSD(lr,S)=JacobiSD{Float64}(lr,S)

convert(::Type{Operator{T}},SD::JacobiSD) where {T}=JacobiSD{T}(SD.lr,SD.S)

domain(op::JacobiSD)=domain(op.S)
domainspace(op::JacobiSD)=op.S
rangespace(op::JacobiSD)=op.lr ? Jacobi(op.S.b+1,op.S.a-1,domain(op.S)) : Jacobi(op.S.b-1,op.S.a+1,domain(op.S))
bandinds(::JacobiSD)=0,0

function getindex(op::JacobiSD,A,k::Integer,j::Integer)
    m=op.lr ? op.S.a : op.S.b
    if k==j
        k+m-1
    else
        zero(eltype(op))
    end
end
