## Evaluation

function Base.getindex{J<:Jacobi}(op::Evaluation{J,Bool},kr::Range)
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
function Base.getindex{J<:Jacobi}(op::Evaluation{J,Float64},kr::Range)
    @assert op.order == 0
    jacobip(kr-1,op.space.a,op.space.b,tocanonical(domain(op),op.x))
end


## Derivative

Derivative(J::Jacobi,k::Integer)=k==1?Derivative{Jacobi,Int,eltype(domain(J))}(J,1):DerivativeWrapper(TimesOperator(Derivative(Jacobi(J.a+1,J.b+1,J.domain),k-1),Derivative{Jacobi,Int,eltype(domain(J))}(J,1)),k)



rangespace{J<:Jacobi}(D::Derivative{J})=Jacobi(D.space.a+D.order,D.space.b+D.order,domain(D))
bandinds{J<:Jacobi}(D::Derivative{J})=0,D.order

function addentries!{J<:Jacobi}(T::Derivative{J},A,kr::Range)
    d=domain(T)
    for k=kr
        A[k,k+1]+=(k+1+T.space.a+T.space.b)./(d.b-d.a)
    end
    A
end


## Integral

function Integral(J::Jacobi,k::Integer)
    if k > 1
        Q=Integral(J,1)
        IntegralWrapper(TimesOperator(Integral(rangespace(Q),k-1),Q),k)
    elseif J.a > 0 && J.b > 0   # we have a simple definition
        Integral{Jacobi,Int,Float64}(J,1)
    else   # convert and then integrate
        sp=Jacobi(J.a+1,J.b+1,domain(J))
        C=Conversion(J,sp)
        Q=Integral(sp,1)
        IntegralWrapper(TimesOperator(Q,C),1)
    end
end


rangespace{J<:Jacobi}(D::Integral{J})=Jacobi(D.space.a-D.order,D.space.b-D.order,domain(D))
bandinds{J<:Jacobi}(D::Integral{J})=-D.order,0

function addentries!{J<:Jacobi}(T::Integral{J},A,kr::Range)
    @assert T.order==1
    d=domain(T)
    for k=intersect(2:kr[end],kr)
        A[k,k-1]+=(d.b-d.a)./(k+T.space.a+T.space.b-2)
    end
    A
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
        Conversion{Jacobi{D},Jacobi{D},Float64}(L,M)
    elseif M.b > L.b+1
        TimesOperator(Conversion(Jacobi(M.a,M.b-1,dm),M),Conversion(L,Jacobi(M.a,M.b-1,dm)))
    else  #if M.a >= L.a+1
        TimesOperator(Conversion(Jacobi(M.a-1,M.b,dm),M),Conversion(L,Jacobi(M.a-1,M.b,dm)))
    end
end

bandinds{J<:Jacobi}(C::Conversion{J,J})=(0,1)



function getdiagonalentry{J<:Jacobi}(C::Conversion{J,J},k,j)
    L=C.domainspace
    if L.b+1==C.rangespace.b
        if j==0
            k==1?1.:(L.a+L.b+k)/(L.a+L.b+2k-1)
        else
            (L.a+k)./(L.a+L.b+2k+1)
        end
    elseif L.a+1==C.rangespace.a
        if j==0
            k==1?1.:(L.a+L.b+k)/(L.a+L.b+2k-1)
        else
            -(L.b+k)./(L.a+L.b+2k+1)
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

function Conversion{m}(A::Ultraspherical{m},B::Jacobi)
    if isapprox(B.a,m-0.5)&&isapprox(B.b,m-0.5)
        Conversion{Ultraspherical{m},Jacobi,Float64}(A,B)
    else
        J=Jacobi(m-0.5,m-0.5,domain(A))
        TimesOperator(Conversion(J,B),Conversion(A,J))
    end
end

bandinds{US<:Ultraspherical,J<:Jacobi}(C::Conversion{US,J})=0,0
bandinds{US<:Ultraspherical,J<:Jacobi}(C::Conversion{J,US})=0,0


function addentries!{J<:Jacobi,CC<:Chebyshev}(C::Conversion{CC,J},A,kr::Range)
    S=rangespace(C)
    @assert isapprox(S.a,-0.5)&&isapprox(S.b,-0.5)
    jp=jacobip(0:kr[end],-0.5,-0.5,1.0)
    for k=kr
        A[k,k]+=1./jp[k]
    end

    A
end

function addentries!{J<:Jacobi,CC<:Chebyshev}(C::Conversion{J,CC},A,kr::Range)
    S=domainspace(C)
    @assert isapprox(S.a,-0.5)&&isapprox(S.b,-0.5)

    jp=jacobip(0:kr[end],-0.5,-0.5,1.0)
    for k=kr
        A[k,k]+=jp[k]
    end

    A
end

function addentries!{US<:Ultraspherical,J<:Jacobi}(C::Conversion{US,J},A,kr::Range)
    S=rangespace(C)
    m=order(US)
    @assert isapprox(S.a,m-0.5)&&isapprox(S.b,m-0.5)
    jp=jacobip(0:kr[end],S.a,S.b,1.0)
    um=Evaluation(US(),1.)[1:kr[end]]
    for k=kr
        A[k,k]+=um[k]./jp[k]
    end

    A
end

function addentries!{US<:Ultraspherical,J<:Jacobi}(C::Conversion{J,US},A,kr::Range)
    m=order(US)
    S=domainspace(C)
    @assert isapprox(S.a,m-0.5)&&isapprox(S.b,m-0.5)

    jp=jacobip(0:kr[end],S.a,S.b,1.0)
    um=Evaluation(Ultraspherical{m}(),1.)[1:kr[end]]
    for k=kr
        A[k,k]+=jp[k]./um[k]
    end

    A
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
        function $OPrule{m}(A::Ultraspherical{m},B::Jacobi)
            if !isapproxinteger(m-0.5-B.a) || !isapproxinteger(m-0.5-B.b)
                NoSpace()
            elseif isapprox(B.a+.5,m)&&isapprox(B.b+.5,m)
                # the spaces are the same
                A
            else
                $OP(Jacobi(A),B)
            end
        end
    end
end

hasconversion(a::Jacobi,b::Ultraspherical)=hasconversion(a,Jacobi(b))
hasconversion(a::Ultraspherical,b::Jacobi)=hasconversion(Jacobi(a),b)




## Special Multiplication
# special multiplication operators exist when multiplying by
# (1+x) or (1-x) by _decreasing_ the parameter.  Thus the


function Multiplication{C<:Chebyshev}(f::Fun{JacobiWeight{C}},S::Jacobi)
    # this implements (1+x)*P and (1-x)*P special case
    # see DLMF (18.9.6)
    if length(f)==1 && ((space(f).α==1 && space(f).β==0 && S.b >0) ||
                        (space(f).α==0 && space(f).β==1 && S.a >0))
        Multiplication{typeof(space(f)),typeof(S),eltype(f),eltype(f)}(f,S)
    else
# default JacobiWeight
        M=Multiplication(Fun(f.coefficients,space(f).space),S)
        rsp=JacobiWeight(space(f).α,space(f).β,rangespace(M))
        MultiplicationWrapper(f,SpaceOperator(M,S,rsp))
    end
end

function rangespace{J<:Jacobi,C<:Chebyshev}(M::Multiplication{JacobiWeight{C},J})
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

bandinds{J<:Jacobi,C<:Chebyshev}(::Multiplication{JacobiWeight{C},J})=-1,0

function addentries!{J<:Jacobi,C<:Chebyshev}(M::Multiplication{JacobiWeight{C},J},A,kr::Range)
    @assert length(M.f)==1
    a,b=domainspace(M).a,domainspace(M).b
    if space(M.f).α==1
        @assert space(M.f).β==0
        # multiply by (1+x)
        for k=kr
            A[k,k]+=2(k+b-1)/(2k+a+b-1)
            if k > 1
                A[k,k-1]+=(2k-2)/(2k+a+b-3)
            end
        end
        A
    elseif space(M.f).β == 1
        @assert space(M.f).α==0
        # multiply by (1-x)
        for k=kr
            A[k,k]+=2(k+a-1)/(2k+a+b-1)
            if k > 1
                A[k,k-1]-=(2k-2)/(2k+a+b-3)
            end
        end
        A
    else
        error("Not implemented")
    end
end


# We can exploit the special multiplication to construct a Conversion
# that decrements parameters


#TODO: general integer decrements
function Conversion{J<:Jacobi}(A::JacobiWeight{J},B::Jacobi)
    if A.α==1.0 && A.β==0.0
        M=Multiplication(Fun([1.],JacobiWeight(1.,0.,domain(A))),A.space)        # multply by (1+x)
        S=SpaceOperator(M,A,rangespace(M))  # this removes the JacobiWeight
        ConversionWrapper(promoterangespace(S,B))
    elseif A.α==0.0 && A.β==1.0
        M=Multiplication(Fun([1.],JacobiWeight(0.,1.,domain(A))),A.space)        # multply by (1-x)
        S=SpaceOperator(M,A,rangespace(M))  # this removes the JacobiWeight
        ConversionWrapper(promoterangespace(S,B))
    else
        error("Not implemented")
    end
end

function maxspace_rule{J<:Jacobi}(A::JacobiWeight{J},B::Jacobi)
    if A.α==1.0 && A.β==0.0 && A.space.b>0
        maxspace(Jacobi(A.space.a,A.space.b-1),B)
    elseif A.α==0.0 && A.β==1.0 && A.space.a>0
        maxspace(Jacobi(A.space.b-1,A.space.a),B)
    else
        maxspace(A,JacobiWeight(0.,0.,B))
    end
end





# represents [b+(1+z)*d/dz] (false) or [a-(1-z)*d/dz] (true)
immutable JacobiSD{T} <:BandedOperator{T}
    lr::Bool
    S::Jacobi
end

JacobiSD(lr,S)=JacobiSD{Float64}(lr,S)

Base.convert{BO<:Operator}(::Type{BO},SD::JacobiSD)=JacobiSD{eltype(BO)}(SD.lr,SD.S)

domain(op::JacobiSD)=domain(op.S)
domainspace(op::JacobiSD)=op.S
rangespace(op::JacobiSD)=op.lr?Jacobi(op.S.a-1,op.S.b+1,domain(op.S)):Jacobi(op.S.a+1,op.S.b-1,domain(op.S))
bandinds(::JacobiSD)=0,0

function addentries!(op::JacobiSD,A,kr::Range)
    m=op.lr?op.S.a:op.S.b
    for k=kr
        A[k,k]+=k+m-1
    end
    A
end
