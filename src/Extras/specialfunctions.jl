## abs


function splitatroots(f::Fun)
    d=domain(f)

    pts=roots(f)

    if isempty(pts)
        f
    else
        da=first(d)
        isapprox(da,pts[1]) ? pts[1] = da : pts = [da,pts]
        db=last(d)
        isapprox(db,pts[end]) ? pts[end] = db : pts = [pts,db]
        Fun(x->f[x],pts)
    end
end

function splitmap(g,d,pts)
    @assert isa(d,AffineDomain)
    da=first(d)
    isapprox(da,pts[1];atol=sqrt(eps(length(d)))) ? pts[1] = da : pts = [da,pts]
    db=last(d)
    isapprox(db,pts[end];atol=sqrt(eps(length(d)))) ? pts[end] = db : pts = [pts,db]
    Fun(g,pts)
end

function Base.abs{S<:RealSpace,T<:Real}(f::Fun{S,T})
    d=domain(f)

    pts=roots(f)

    if isempty(pts)
        sign(first(f))*f
    else
        splitmap(x->abs(f[x]),d,pts)
    end
end

function Base.abs{S,T}(f::Fun{S,T})
    d=domain(f)

    pts=roots(f)

    if isempty(pts)
        Fun(x->abs(f[x]),space(f))
    else
        splitmap(x->abs(f[x]),d,pts)
    end
end



function Base.sign{S<:RealSpace,T<:Real}(f::Fun{S,T})
    d=domain(f)

    pts=roots(f)

    if isempty(pts)
        sign(first(f))*one(T,f.space)
    else
        @assert isa(d,AffineDomain)
        da=first(d)
        isapprox(da,pts[1];atol=sqrt(eps(length(d)))) ? pts[1] = da : pts = [da,pts]
        db=last(d)
        isapprox(db,pts[end];atol=sqrt(eps(length(d)))) ? pts[end] = db : pts = [pts,db]
        midpts = .5(pts[1:end-1]+pts[2:end])
        Fun([sign(f[midpts])],pts)
    end
end

# division by fun

function ./{S,T,U,V}(c::Fun{S,T},f::Fun{U,V})
    r=roots(f)
    tol=10eps()
    if length(r)==0 || norm(c[r])<tol
        linsolve(Multiplication(f,space(c)),c;tolerance=tol)
    else
        c.*(1./f)
    end
end

function ./{S,T}(c::Number,f::Fun{S,T})
    r=roots(f)
    tol=10eps()
    @assert length(r)==0
    linsolve(Multiplication(f,space(f)),c*ones(space(f));tolerance=tol)
end

function ./(c::Number,f::Fun{Chebyshev})
    fc = Fun(canonicalcoefficients(f),Interval())
    r = roots(fc)
    x = Fun(identity)

    # if domain f is small then the pts get projected in
    tol = 50eps()/length(domain(f))


    if length(r) == 0
        linsolve(Multiplication(f,space(f)),c;tolerance=tol)
    elseif length(r) == 1 && abs(abs(r[1]) - 1.) < tol
        if sign(r[1]) < 0
            g = divide_singularity(-1,fc)  # divide by 1+x
            Fun(canonicalcoefficients(c./g),JacobiWeight(-1,0,domain(f)))
        else
            g = divide_singularity(1,fc)  # divide by 1-x
            Fun(canonicalcoefficients(c./g),JacobiWeight(0,-1,domain(f)))
        end
    elseif length(r) ==2 && abs(r[1]+1) < tol && abs(r[2]-1) < tol
        g = divide_singularity(fc) # divide by 1-x^2
        # divide out singularities, tolerance needs to be chosen since we don't get
        # spectral convergence
        # TODO: switch to dirichlet basis
        Fun(canonicalcoefficients(c./g),JacobiWeight(-1,-1,domain(f)))
    else
        #split at the roots
        c./splitatroots(f)
    end
end

./{S<:PiecewiseSpace}(c::Number,f::Fun{S})=depiece(map(f->c./f,pieces(f)))
function ./{S<:MappedSpace}(c::Number,f::Fun{S})
    g=c./Fun(coefficients(f),space(f).space)
    if isa(space(g),JacobiWeight)
        Fun(coefficients(g),JacobiWeight(space(g).α,space(g).β,MappedSpace(domain(f),space(g).space)))
    else
        Fun(coefficients(g),MappedSpace(domain(f),space(g)))
    end
end

function .^{S<:MappedChebyshev}(f::Fun{S},k::Float64)
    fc = Fun(f.coefficients) #Project to interval
    x=Fun(identity)

    r = sort(roots(fc))


    @assert length(r) <= 2

    if length(r) == 0
        Fun(Fun(x->fc[x]^k).coefficients,space(f))
    elseif length(r) == 1
        @assert isapprox(abs(r[1]),1)

        if isapprox(r[1],1.)
            Fun(coefficients(divide_singularity(+1,fc)^k),JacobiWeight(0.,k,space(f)))
        else
            Fun(coefficients(divide_singularity(-1,fc)^k),JacobiWeight(k,0.,space(f)))
        end
    else
        @assert isapprox(r[1],-1)
        @assert isapprox(r[2],1)

        Fun(coefficients(divide_singularity(fc)^k),JacobiWeight(k,k,space(f)))
    end
end

Base.sqrt{S,T}(f::Fun{S,T})=f^0.5
Base.cbrt{S,T}(f::Fun{S,T})=f^(1/3)

## We use \ as the Fun constructor might miss isolated features
function Base.exp{S<:Ultraspherical,T}(f::Fun{S,T})
    rf=chop(real(f),eps())
    xmax=rf.coefficients==[0.]?last(domain(rf)):indmax(rf)
    B=Evaluation(space(f),xmax)
    D=Derivative(space(f))
    A=[B,D-differentiate(f)]
    (A\[1.])*exp(f[xmax])
end
Base.exp2(f::Fun) = exp(log(2)*f)
Base.exp10(f::Fun) = exp(log(10)*f)

function Base.log{S<:Ultraspherical,T}(f::Fun{S,T})
    rf=chop(real(f),eps())
    xmax=rf.coefficients==[0.]?last(domain(rf)):indmax(rf)
    B=Evaluation(space(f),xmax)
    D=Derivative(space(f))
    A=[B,f*D]
    A\[log(f[xmax]),differentiate(f)]
end
Base.log2(f::Fun) = log(f)/log(2)
Base.log10(f::Fun) = log(f)/log(10)

for (op,ODE,RHS,growth) in ((:(Base.erf),"fp*D^2+(2f*fp^2-fpp)*D","0*fp^3",:(real)),
                        (:(Base.erfc),"fp*D^2+(2f*fp^2-fpp)*D","0*fp^3",:(real)),
                        (:(Base.airyai),"fp*D^2-fpp*D-f*fp^3","0*fp^3",:(imag)),
                        (:(Base.airybi),"fp*D^2-fpp*D-f*fp^3","0*fp^3",:(imag)),
                        (:(Base.airyaiprime),"fp*D^2-fpp*D-f*fp^3","airyai(f)*fp^3",:(imag)),
                        (:(Base.airybiprime),"fp*D^2-fpp*D-f*fp^3","airybi(f)*fp^3",:(imag)),)
    parsedODE = parse(ODE)
    parsedRHS = parse(RHS)
    @eval begin
        function $op{S<:Ultraspherical,T}(f::Fun{S,T})
            g=chop($growth(f),eps())
            xmin=g.coefficients==[0.]?first(domain(g)):indmin(g)
            xmax=g.coefficients==[0.]?last(domain(g)):indmax(g)
            opfxmin,opfxmax = $op(f[xmin]),$op(f[xmax])
            opmax = maxabs((opfxmin,opfxmax))
            while opmax≤10eps()
                xmin,xmax = rand(domain(f)),rand(domain(f))
                opfxmin,opfxmax = $op(f[xmin]),$op(f[xmax])
                opmax = maxabs((opfxmin,opfxmax))
            end
            D=Derivative(space(f))
            fp=differentiate(f)
            fpp=differentiate(fp)
            B=[Evaluation(space(f),xmin),Evaluation(space(f),xmax)]
            ([B,eval($parsedODE)]\[opfxmin/opmax,opfxmax/opmax,eval($parsedRHS)/opmax])*opmax
        end
    end
end

function Base.erfcx{S<:Ultraspherical,T}(f::Fun{S,T})
    rf=chop(real(f),eps())
    xmin=rf.coefficients==[0.]?first(domain(rf)):indmin(rf)
    opfxmin = erfcx(f[xmin])
    D=Derivative(space(f))
    fp=differentiate(f)
    B=Evaluation(space(f),xmin)
    ([B,D-2f*fp]\[1.0,-2fp/sqrt(π)/opfxmin])*opfxmin
end

function Base.dawson{S<:Ultraspherical,T}(f::Fun{S,T})
    rf=chop(real(f),eps())
    xmin=rf.coefficients==[0.]?first(domain(rf)):indmin(rf)
    opfxmin = dawson(f[xmin])
    D=Derivative(space(f))
    fp=differentiate(f)
    B=Evaluation(space(f),xmin)
    ([B,D+2f*fp]\[1.0,fp/opfxmin])*opfxmin
end

for (op,ODE) in ((:(Base.hankelh1),"f^2*fp*D^2+(f*fp^2-f^2*fpp)*D+(f^2-ν^2)*fp^3"),
                 (:(Base.hankelh2),"f^2*fp*D^2+(f*fp^2-f^2*fpp)*D+(f^2-ν^2)*fp^3"),
                 (:(Base.hankelh1x),"f^2*fp*D^2+((2im*f^2+f)*fp^2-f^2*fpp)*D+(im*f-ν^2)*fp^3"),
                 (:(Base.hankelh2x),"f^2*fp*D^2+((-2im*f^2+f)*fp^2-f^2*fpp)*D+(-im*f-ν^2)*fp^3"))
    parsedODE = parse(ODE)
    @eval begin
        function $op{S<:Ultraspherical,T}(ν,f::Fun{S,T})
            imf=chop(imag(f),eps())
            xmin=imf.coefficients==[0.]?first(domain(imf)):indmin(imf)
            xmax=imf.coefficients==[0.]?last(domain(imf)):indmax(imf)
            opfxmin,opfxmax = $op(ν,f[xmin]),$op(ν,f[xmax])
            opmax = maxabs((opfxmin,opfxmax))
            D=Derivative(space(f))
            fp=differentiate(f)
            fpp=differentiate(fp)
            B=[Evaluation(space(f),xmin),Evaluation(space(f),xmax)]
            ([B,eval($parsedODE)]\[opfxmin/opmax,opfxmax/opmax])*opmax
        end
    end
end


Base.sin{S<:FunctionSpace{Float64},T<:Real}(f::Fun{S,T}) = imag(exp(im*f))
Base.cos{S<:FunctionSpace{Float64},T<:Real}(f::Fun{S,T}) = real(exp(im*f))

Base.sin(f::Fun) = (exp(im*f)-exp(-im*f))/2im
Base.cos(f::Fun) = (exp(im*f)+exp(-im*f))/2
Base.tan(f::Fun) = sin(f)/cos(f)     ##TODO: the spacepromotion doesn't work for tan for a domain including zeros of cos inside.

Base.asin{S<:FunctionSpace{Float64},T<:Real}(f::Fun{S,T}) = imag(log(im*f+sqrt(1-f^2)))
Base.acos{S<:FunctionSpace{Float64},T<:Real}(f::Fun{S,T}) = -imag(log(f-im*sqrt(1-f^2)))

##TODO: If f and sqrt(1-f^2) are in different spaces, then the log's argument could be im*f⊕sqrt(1-f^2)... in the SumSpace. This also requires log for more spaces (or a constructor with edge detection).
Base.asin(f::Fun) = -im*log(im*f+sqrt(1-f^2))
Base.acos(f::Fun) = im*log(f-im*sqrt(1-f^2))
#Base.atan(f::Fun) = asin(f/sqrt(f^2+1)) #This is inaccurate.

Base.sinh(f::Fun) = (exp(f)-exp(-f))/2
Base.cosh(f::Fun) = (exp(f)+exp(-f))/2
#Base.tanh(f::Fun) = sinh(f)/cosh(f) #This is inaccurate.

Base.asinh(f::Fun) = log(f+sqrt(f^2+1))
Base.acosh(f::Fun) = log(f+sqrt(f+1)*sqrt(f-1))
Base.atanh(f::Fun) = (log(1+f)-log(1-f))/2

for (op,oprecip,opinv,opinvrecip) in ((:(Base.sin),:(Base.csc),:(Base.asin),:(Base.acsc)),
                                      (:(Base.cos),:(Base.sec),:(Base.acos),:(Base.asec)),
                                      (:(Base.tan),:(Base.cot),:(Base.atan),:(Base.acot)),
                                      (:(Base.sinh),:(Base.csch),:(Base.asinh),:(Base.acsch)),
                                      (:(Base.cosh),:(Base.sech),:(Base.acosh),:(Base.asech)),
                                      (:(Base.tanh),:(Base.coth),:(Base.atanh),:(Base.acoth)))
    @eval begin
        $oprecip(f::Fun) = 1/$op(f)
        $opinvrecip(f::Fun) = $opinv(1/f)
    end
end

rad2deg(f::Fun) = 180*f/π
deg2rad(f::Fun) = π*f/180

for (op,opd,opinv,opinvd) in ((:(Base.sin),:(Base.sind),:(Base.asin),:(Base.asind)),
                              (:(Base.cos),:(Base.cosd),:(Base.acos),:(Base.acosd)),
                              (:(Base.tan),:(Base.tand),:(Base.atan),:(Base.atand)),
                              (:(Base.sec),:(Base.secd),:(Base.asec),:(Base.asecd)),
                              (:(Base.csc),:(Base.cscd),:(Base.acsc),:(Base.acscd)),
                              (:(Base.cot),:(Base.cotd),:(Base.acot),:(Base.acotd)))
    @eval begin
        $opd(f::Fun) = $op(deg2rad(f))
        $opinvd(f::Fun) = rad2deg($opinv(f))
    end
end

#Won't get the zeros exactly 0 anyway so at least this way the length is smaller.
Base.sinpi(f::Fun) = sin(π*f)
Base.cospi(f::Fun) = cos(π*f)
Base.sinc(f::Fun) = sin(π*f)/(π*f)
Base.cosc(f::Fun) = cos(π*f)/f - sin(π*f)/(π*f^2)

Base.erfi{S<:FunctionSpace{Float64},T<:Real}(f::Fun{S,T}) = imag(erf(im*f))
Base.erfi(f::Fun) = -im*erf(im*f)

function Base.airy(k::Number,f::Fun)
    if k == 0
        airyai(f)
    elseif k == 1
        airyaiprime(f)
    elseif k == 2
        airybi(f)
    elseif k == 3
        airybiprime(f)
    else
        error("invalid argument")
    end
end

Base.besselh(ν,k::Integer,f::Fun) = k == 1? hankelh1(ν,f) : k == 2 ? hankelh2(ν,f) : throw(Base.Math.AmosException(1))

Base.besselj{S<:FunctionSpace{Float64},T<:Real}(ν,f::Fun{S,T}) = real(hankelh1(ν,f))
Base.bessely{S<:FunctionSpace{Float64},T<:Real}(ν,f::Fun{S,T}) = imag(hankelh1(ν,f))
Base.besselj(ν,f::Fun) = (hankelh1(ν,f)+hankelh2(ν,f))/2
Base.bessely(ν,f::Fun) = (hankelh1(ν,f)-hankelh2(ν,f))/2im

# TODO: The following definitions are only valid in three-quarters of the complex plane (from -π ≤ arg z ≤ π/2).
# It seems like it would be difficult to base these definitions completely on the hankel functions.
Base.besseli{S<:FunctionSpace{Float64},T<:Real}(ν,f::Fun{S,T}) = real(exp(-π*ν*im/2)*(hankelh1(ν,im*f)+hankelh2(ν,im*f))/2)
Base.besselk{S<:FunctionSpace{Float64},T<:Real}(ν,f::Fun{S,T}) = real(im*π/2*exp(π*ν*im/2)*hankelh1(ν,im*f))
Base.besselkx{S<:FunctionSpace{Float64},T<:Real}(ν,f::Fun{S,T}) = real(im*π/2*exp(π*ν*im/2)*hankelh1x(ν,im*f))
Base.besseli(ν,f::Fun) = exp(-π*ν*im/2)*(hankelh1(ν,im*f)+hankelh2(ν,im*f))/2
Base.besselk(ν,f::Fun) = im*π/2*exp(π*ν*im/2)*hankelh1(ν,im*f)
Base.besselkx(ν,f::Fun) = im*π/2*exp(π*ν*im/2)*hankelh1x(ν,im*f)

for jy in ("j","y"), ν in (0,1)
    bjy = symbol(string("bessel",jy))
    bjynu = parse(string("Base.bessel",jy,ν))
    @eval begin
        $bjynu(f::Fun) = $bjy($ν,f)
    end
end

## Miscellaneous
for op in (:(Base.expm1),:(Base.log1p),:(Base.atan),:(Base.tanh),
           :(Base.erfinv),:(Base.erfcinv),:(Base.beta),:(Base.lbeta),
           :(Base.eta),:(Base.zeta),:(Base.gamma),:(Base.lgamma),
           :(Base.polygamma),:(Base.invdigamma),:(Base.digamma),:(Base.trigamma))
    @eval begin
        $op{S,T}(f::Fun{S,T})=Fun(x->$op(f[x]),domain(f))
    end
end
