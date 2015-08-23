## Evaluation

function Base.getindex(op::Evaluation{Jacobi,Bool},kr::Range)
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
function Base.getindex(op::Evaluation{Jacobi,Float64},kr::Range)
    @assert op.order == 0
    jacobip(kr-1,op.space.a,op.space.b,tocanonical(domain(op),op.x))
end


## Derivative

Derivative(J::Jacobi,k::Integer)=k==1?Derivative{Jacobi,Int,eltype(domain(J))}(J,1):DerivativeWrapper(TimesOperator(Derivative(Jacobi(J.a+1,J.b+1,J.domain),k-1),Derivative{Jacobi,Int,eltype(domain(J))}(J,1)),k)



rangespace(D::Derivative{Jacobi})=Jacobi(D.space.a+D.order,D.space.b+D.order,domain(D))
bandinds(D::Derivative{Jacobi})=0,D.order

function addentries!(T::Derivative{Jacobi},A,kr::Range)
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


rangespace(D::Integral{Jacobi})=Jacobi(D.space.a-D.order,D.space.b-D.order,domain(D))
bandinds(D::Integral{Jacobi})=-D.order,0

function addentries!(T::Integral{Jacobi},A,kr::Range)
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
    if isapprox(M.a,L.a) && isapprox(M.b,L.b)
        SpaceOperator(IdentityOperator(),L,M)
    elseif (isapprox(M.b,L.b+1) && isapprox(M.a,L.a)) || (isapprox(M.b,L.b) && isapprox(M.a,L.a+1))
        Conversion{Jacobi,Jacobi,Float64}(L,M)
    elseif M.b > L.b+1
        TimesOperator(Conversion(Jacobi(M.a,M.b-1,dm),M),Conversion(L,Jacobi(M.a,M.b-1,dm)))
    else  #if M.a >= L.a+1
        TimesOperator(Conversion(Jacobi(M.a-1,M.b,dm),M),Conversion(L,Jacobi(M.a-1,M.b,dm)))
    end
end

bandinds(C::Conversion{Jacobi,Jacobi})=(0,1)



function getdiagonalentry(C::Conversion{Jacobi,Jacobi},k,j)
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
    if !isapproxinteger(A.a-B.a) || !isapproxinteger(B.a-B.a)
        NoSpace()
    else
        Jacobi(min(A.a,B.a),min(A.b,B.b),domain(A))
    end
end



## Ultraspherical Conversion

# Assume m is compatible

function Conversion(A::Chebyshev,B::Jacobi)
    if isapprox(B.a,-0.5)&&isapprox(B.b,-0.5)
        Conversion{Chebyshev,Jacobi,Float64}(A,B)
    else
        J=Jacobi(-0.5,-0.5,domain(A))
        TimesOperator(Conversion(J,B),Conversion(A,J))
    end
end

bandinds{m}(C::Conversion{Ultraspherical{m},Jacobi})=0,0
bandinds{m}(C::Conversion{Jacobi,Ultraspherical{m}})=0,0


function addentries!(C::Conversion{Chebyshev,Jacobi},A,kr::Range)
    S=rangespace(C)
    @assert isapprox(S.a,-0.5)&&isapprox(S.b,-0.5)
    jp=jacobip(0:kr[end],-0.5,-0.5,1.0)
    for k=kr
        A[k,k]+=1./jp[k]
    end

    A
end

function addentries!(C::Conversion{Jacobi,Chebyshev},A,kr::Range)
    S=domainspace(C)
    @assert isapprox(S.a,-0.5)&&isapprox(S.b,-0.5)

    jp=jacobip(0:kr[end],-0.5,-0.5,1.0)
    for k=kr
        A[k,k]+=jp[k]
    end

    A
end



function conversion_rule{m}(A::Ultraspherical{m},B::Jacobi)
    if B.a+.5==m&&B.b+.5==m
        # the spaces are the same
        A
    else
        conversion_type(Jacobi(m-0.5,m-0.5,domain(A)),B)
    end
end

function maxspace{m}(A::Ultraspherical{m},B::Jacobi)
    if B.a+.5==m&&B.b+.5==m
        # the spaces are the same
        A
    else
        maxspace(Jacobi(m-0.5,m-0.5,domain(A)),B)
    end
end

maxspace(A::Jacobi,B::Ultraspherical)=maxspace(B,A)
