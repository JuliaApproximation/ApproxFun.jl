## abs
splitmap(g,d::Domain,pts) = Fun(g,d \ Set(pts))


function splitatroots(f::Fun)
    d=domain(f)
    pts=union(roots(f)) # union removes multiplicities
    splitmap(x->f(x),d,pts)
end

function abs(f::Fun{S,T}) where {S<:RealUnivariateSpace,T<:Real}
    d=domain(f)
    pts = iszero(f) ? T[] : roots(f)
    splitmap(x->abs(f(x)),d,pts)
end

function abs(f::Fun)
    d=domain(f)

    pts = iszero(f) ? cfstype(f)[] : roots(f)

    if isempty(pts)
        # This makes sure Laurent returns real type
        real(Fun(abs∘f,space(f)))
    else
        splitmap(abs∘f,d,pts)
    end
end


midpoints(d::Segment) = [mean(d)]
midpoints(d::UnionDomain) = mapreduce(midpoints,vcat,d.domains)


for OP in (:sign,:angle)
    @eval function $OP(f::Fun{S,T}) where {S<:RealUnivariateSpace,T<:Real}
        d=domain(f)

        pts = iszero(f) ? T[] : roots(f)

        if isempty(pts)
            $OP(first(f))*one(f)
        else
            d = d \ Set(pts)
            midpts = midpoints(d)
            Fun(d, $OP.(f.(midpts)))
        end
    end
end

for op in (:(max),:(min))
    @eval begin
        function $op(f::Fun{S,T1},g::Fun{V,T2}) where {S<:RealUnivariateSpace,V<:RealUnivariateSpace,T1<:Real,T2<:Real}
            h=f-g
            d=domain(h)
            pts=iszero(h) ? cfstype(h)[] : roots(h)
            splitmap(x->$op(f(x),g(x)),d,pts)
        end
        $op(f::Fun{S,T},g::Real) where {S<:RealUnivariateSpace,T<:Real} = $op(f,Fun(g,domain(f)))
        $op(f::Real,g::Fun{S,T}) where {S<:RealUnivariateSpace,T<:Real} = $op(Fun(f,domain(g)),g)
    end
end

isfinite(f::Fun) = isfinite(maximum(abs,f)) && isfinite(minabs(f))

# division by fun

function /(c::Fun,f::Fun)
    d=domain(f)
    @assert domain(c) == d
    cd = canonicaldomain(f)
    if typeof(d)!=typeof(cd)
        # project first to simplify
        return setdomain(setdomain(c,cd)/setdomain(f,cd),d)
    end

    r=roots(f)
    tol=10eps(promote_type(cfstype(c),cfstype(f)))
    if length(r)==0 || norm(c.(r))<tol
        \(Multiplication(f,space(c)),c;tolerance=tol)
    else
        c*(1/f)
    end
end

function /(c::Number,f::Fun)
    r=roots(f)
    tol=10eps(promote_type(typeof(c),cfstype(f)))
    @assert length(r)==0
    \(Multiplication(f,space(f)),c*one(space(f));tolerance=tol)
end

# project to interval if we are not on the interview
# TODO: need to work out how to set piecewise domain


scaleshiftdomain(f::Fun,sc,sh) = setdomain(f,sc*domain(f)+sh)

/(c::Number,f::Fun{Ultraspherical{λ,DD,RR}}) where {λ,DD,RR} = c/Fun(f,Chebyshev(domain(f)))
/(c::Number,f::Fun{Jacobi{DD,RR}}) where {DD,RR} = c/Fun(f,Chebyshev(domain(f)))

/(c::Number,f::Fun{C}) where {C<:Chebyshev}=setdomain(c/setcanonicaldomain(f),domain(f))
function /(c::Number,f::Fun{Chebyshev{DD,RR}}) where {DD<:IntervalOrSegment,RR}
    fc = setcanonicaldomain(f)
    d=domain(f)
    # if domain f is small then the pts get projected in
    tol = 200eps(promote_type(typeof(c),cfstype(f)))*norm(f.coefficients,1)

    # we prune out roots at the boundary first
    if ncoefficients(f) == 1
        return Fun(c/f.coefficients[1],space(f))
    elseif ncoefficients(f) == 2
        if isempty(roots(f))
            return \(Multiplication(f,space(f)),c;tolerance=0.05tol)
        elseif isapprox(fc.coefficients[1],fc.coefficients[2])
            # we check directly for f*(1+x)
            return Fun(JacobiWeight(-1,0,space(f)),[c/fc.coefficients[1]])
        elseif isapprox(fc.coefficients[1],-fc.coefficients[2])
            # we check directly for f*(1-x)
            return Fun(JacobiWeight(0,-1,space(f)),[c/fc.coefficients[1]])
        else
            # we need to split at the only root
            return c/splitatroots(f)
        end
    elseif abs(first(fc)) ≤ tol
        #left root
        g=divide_singularity((1,0),fc)
        p=c/g
        x = Fun(identity,domain(p))
        return scaleshiftdomain(p/(1+x),complexlength(d)/2,mean(d) )
    elseif abs(last(fc)) ≤ tol
        #right root
        g=divide_singularity((0,1),fc)
        p=c/g
        x=Fun(identity,domain(p))
        return scaleshiftdomain(p/(1-x),complexlength(d)/2,mean(d) )
    else
        r = roots(fc)

        if length(r) == 0
            return \(Multiplication(f,space(f)),c;tolerance=0.05tol)
        elseif abs(last(r)+1.0)≤tol  # double check
            #left root
            g=divide_singularity((1,0),fc)
            p=c/g
            x=Fun(identity,domain(p))
            return scaleshiftdomain(p/(1+x),complexlength(d)/2,mean(d) )
        elseif abs(last(r)-1.0)≤tol  # double check
            #right root
            g=divide_singularity((0,1),fc)
            p=c/g
            x=Fun(identity,domain(p))
            return scaleshiftdomain(p/(1-x),complexlength(d)/2,mean(d) )
        else
            # no roots on the boundary
            return c/splitatroots(f)
        end
    end
end

^(f::Fun{<:PolynomialSpace},k::Integer) = intpow(f,k)
function ^(f::Fun{<:PolynomialSpace}, k::Real)
    T = cfstype(f)
    RT = real(T)
    # Need to think what to do if this is ever not the case..
    sp = space(f)
    fc = setcanonicaldomain(f) #Project to interval
    csp = space(fc)

    r = sort(roots(fc))
    #TODO divideatroots
    @assert length(r) <= 2

    if length(r) == 0
        setdomain(Fun((x->x^k) ∘ fc,csp),domain(f))  # using ∘ supports fast transforms for fc
    elseif length(r) == 1
        @assert isapprox(abs(r[1]),1)

        if isapprox(r[1], 1)
            Fun(JacobiWeight(zero(RT),k,sp),coefficients(divide_singularity(true,fc)^k,csp))
        else
            Fun(JacobiWeight(k,zero(RT),sp),coefficients(divide_singularity(false,fc)^k,csp))
        end
    else
        @assert isapprox(r[1],-1)
        @assert isapprox(r[2],1)

        Fun(JacobiWeight(k,k,sp),coefficients(divide_singularity(fc)^k,csp))
    end
end

#TODO: implement
^(f::Fun{Jacobi},k::Integer) = intpow(f,k)
^(f::Fun{Jacobi},k::Real) = Fun(f,Chebyshev)^k

# Default is just try solving ODE
function ^(f::Fun{S,T},β) where {S,T}
    A=Derivative()-β*differentiate(f)/f
    B=Evaluation(first(domain(f)))
    [B;A]\[first(f)^β;0]
end

sqrt(f::Fun{S,T}) where {S,T} = f^0.5
cbrt(f::Fun{S,T}) where {S,T} = f^(1/3)

## We use \ as the Fun constructor might miss isolated features

## First order functions


log(f::Fun) = cumsum(differentiate(f)/f)+log(first(f))

# function log{MS<:MappedSpace}(f::Fun{MS})
#     g=log(Fun(f.coefficients,space(f).space))
#     Fun(g.coefficients,MappedSpace(domain(f),space(g)))
# end

# project first to [-1,1] to avoid issues with
# complex derivative
function log(f::Fun{US}) where US<:Union{Ultraspherical,Chebyshev}
    if domain(f)==Segment()
        r = sort(roots(f))
        #TODO divideatroots
        @assert length(r) <= 2

        if length(r) == 0
            cumsum(differentiate(f)/f)+log(first(f))
        elseif length(r) == 1
            @assert isapprox(abs(r[1]),1)

            if isapprox(r[1],1.)
                g=divide_singularity(true,f)
                lg=Fun(LogWeight(0.,1.,Chebyshev()),[1.])
                if isapprox(g,1.)  # this means log(g)~0
                    lg
                else # log((1-x)) + log(g)
                    lg⊕log(g)
                end
            else
                g=divide_singularity(false,f)
                lg=Fun(LogWeight(1.,0.,Chebyshev()),[1.])
                if isapprox(g,1.)  # this means log(g)~0
                    lg
                else # log((1+x)) + log(g)
                    lg⊕log(g)
                end
           end
        else
            @assert isapprox(r[1],-1)
            @assert isapprox(r[2],1)

            g=divide_singularity(f)
            lg=Fun(LogWeight(1.,1.,Chebyshev()),[1.])
            if isapprox(g,1.)  # this means log(g)~0
                lg
            else # log((1+x)) + log(g)
                lg⊕log(g)
            end
        end
    else
        # this makes sure differentiate doesn't
        # make the function complex
        g=log(setdomain(f,Segment()))
        setdomain(g,domain(f))
    end
end


function log(f::Fun{Fourier{D,R},T}) where {T<:Real,D,R}
    if isreal(domain(f))
        cumsum(differentiate(f)/f)+log(first(f))
    else
        # this makes sure differentiate doesn't
        # make the function complex
        g=log(setdomain(f,PeriodicInterval()))
        setdomain(g,domain(f))
    end
end


atan(f::Fun)=cumsum(f'/(1+f^2))+atan(first(f))


# this is used to find a point in which to impose a boundary
# condition in calculating secial functions
function specialfunctionnormalizationpoint(op,growth,f)
    g=chop(growth(f),eps(cfstype(f)))
    xmin = isempty(g.coefficients) ? first(domain(g)) : argmin(g)
    xmax = isempty(g.coefficients) ? last(domain(g)) : argmax(g)
    opfxmin,opfxmax = op(f(xmin)),op(f(xmax))
    opmax = maximum(abs,(opfxmin,opfxmax))
    if abs(opfxmin) == opmax xmax,opfxmax = xmin,opfxmin end
    xmax,opfxmax,opmax
end



# ODE gives the first order ODE a special function op satisfies,
# RHS is the right hand side
# growth says what to use to choose a good point to impose an initial condition
for (op,ODE,RHS,growth) in ((:(exp),"D-f'","0",:(real)),
                            (:(asinh),"sqrt(f^2+1)*D","f'",:(real)),
                            (:(acosh),"sqrt(f^2-1)*D","f'",:(real)),
                            (:(atanh),"(1-f^2)*D","f'",:(real)),
                            (:(erfcx),"D-2f*f'","-2f'/sqrt(π)",:(real)),
                            (:(dawson),"D+2f*f'","f'",:(real)))
    L,R = Meta.parse(ODE),Meta.parse(RHS)
    @eval begin
        # depice before doing op
        $op(f::Fun{PW}) where {PW<:Union{PiecewiseSpace,ContinuousSpace}} =
            Fun(map(f->$op(f),components(f)),PiecewiseSpace)

        # We remove the MappedSpace
        # function $op{MS<:MappedSpace}(f::Fun{MS})
        #     g=exp(Fun(f.coefficients,space(f).space))
        #     Fun(g.coefficients,MappedSpace(domain(f),space(g)))
        # end
        function $op(fin::Fun{S,T}) where {S,T}
            f=setcanonicaldomain(fin)  # removes possible issues with roots

            xmax,opfxmax,opmax=specialfunctionnormalizationpoint($op,$growth,f)
            # we will assume the result should be smooth on the domain
            # even if f is not
            # This supports Line/Rays
            D=Derivative(domain(f))
            B=Evaluation(domainspace(D),xmax)
            #([B,eval($L)]\[opfxmax/opmax,eval($R)/opmax])*opmax
            u=\([B,eval($L)],Any[opfxmax/opmax,eval($R)/opmax];tolerance=eps(T))*opmax

            setdomain(u,domain(fin))
        end
    end
end

# JacobiWeight explodes, we want to ensure the solution incorporates the fact
# that exp decays rapidly
function exp(f::Fun{JW}) where JW<:JacobiWeight
    if !isa(domain(f),Segment)
        # project first to get better derivative behaviour
        return setdomain(exp(setdomain(f,Segment())),domain(f))
    end

    S=space(f)
    q=Fun(S.space,f.coefficients)
    if isapprox(S.α,0.) && isapprox(S.β,0.)
        exp(q)
    elseif S.β < 0 && isapprox(first(q),0.)
        # this case can remove the exponential decay
        exp(Fun(f,JacobiWeight(S.β+1,S.α,S.space)))
    elseif S.α < 0 && isapprox(last(q),0.)
        exp(Fun(f,JacobiWeight(S.β,S.α+1,S.space)))
    elseif S.β > 0 && isapproxinteger(S.β)
        exp(Fun(f,JacobiWeight(0.,S.α,S.space)))
    elseif S.α > 0 && isapproxinteger(S.α)
        exp(Fun(f,JacobiWeight(S.β,0.,S.space)))
    else
        #find normalization point
        xmax,opfxmax,opmax=specialfunctionnormalizationpoint(exp,real,f)

        if S.α < 0 && S.β < 0
            # provided both are negative, we get exponential decay on both ends
            @assert real(first(q)) < 0 && real(last(q)) < 0
            s=JacobiWeight(2.,2.,domain(f))
        elseif S.β < 0 && isapprox(S.α,0.)
            @assert real(first(q)) < 0
            s=JacobiWeight(2.,0.,domain(f))
        elseif S.α < 0 && isapprox(S.β,0.)
            @assert real(last(q)) < 0
            s=JacobiWeight(0.,2.,domain(f))
        else
            error("exponential of fractional power, not implemented")
        end

        D=Derivative(s)
        B=Evaluation(s,xmax)

        \([B,D-f'],Any[opfxmax/opmax,0.];tolerance=eps(cfstype(f)))*opmax
    end
end







acos(f::Fun)=cumsum(-f'/sqrt(1-f^2))+acos(first(f))
asin(f::Fun)=cumsum(f'/sqrt(1-f^2))+asin(first(f))



## Second order functions


sin(f::Fun{S,T}) where {S<:RealSpace,T<:Real} = imag(exp(im*f))
cos(f::Fun{S,T}) where {S<:RealSpace,T<:Real} = real(exp(im*f))
sin(f::Fun{S,T}) where {S<:Union{Ultraspherical,Chebyshev},T<:Real} = imag(exp(im*f))
cos(f::Fun{S,T}) where {S<:Union{Ultraspherical,Chebyshev},T<:Real} = real(exp(im*f))



for (op,ODE,RHS,growth) in ((:(erf),"f'*D^2+(2f*f'^2-f'')*D","0",:(imag)),
                            (:(erfi),"f'*D^2-(2f*f'^2+f'')*D","0",:(real)),
                            (:(sin),"f'*D^2-f''*D+f'^3","0",:(imag)),
                            (:(cos),"f'*D^2-f''*D+f'^3","0",:(imag)),
                            (:(sinh),"f'*D^2-f''*D-f'^3","0",:(real)),
                            (:(cosh),"f'*D^2-f''*D-f'^3","0",:(real)),
                            (:(airyai),"f'*D^2-f''*D-f*f'^3","0",:(imag)),
                            (:(airybi),"f'*D^2-f''*D-f*f'^3","0",:(imag)),
                            (:(airyaiprime),"f'*D^2-f''*D-f*f'^3","airyai(f)*f'^3",:(imag)),
                            (:(airybiprime),"f'*D^2-f''*D-f*f'^3","airybi(f)*f'^3",:(imag)))
    L,R = Meta.parse(ODE),Meta.parse(RHS)
    @eval begin
        function $op(fin::Fun{S,T}) where {S,T}
            f=setcanonicaldomain(fin)

            g=chop($growth(f),eps(T))
            xmin = isempty(g.coefficients) ? first(domain(g)) : argmin(g)
            xmax = isempty(g.coefficients) ? last(domain(g)) : argmax(g)
            opfxmin,opfxmax = $op(f(xmin)),$op(f(xmax))
            opmax = maximum(abs,(opfxmin,opfxmax))
            while opmax≤10eps(T) || abs(f(xmin)-f(xmax))≤10eps(T)
                xmin,xmax = rand(domain(f)),rand(domain(f))
                opfxmin,opfxmax = $op(f(xmin)),$op(f(xmax))
                opmax = maximum(abs,(opfxmin,opfxmax))
            end
            D=Derivative(space(f))
            B=[Evaluation(space(f),xmin),Evaluation(space(f),xmax)]
            u=\([B;eval($L)],[opfxmin/opmax;opfxmax/opmax;eval($R)/opmax];
                        tolerance=10ncoefficients(f)*eps(T)*opmax)*opmax

            setdomain(u,domain(fin))
        end
    end
end

erfc(f::Fun) = 1-erf(f)

## Second order functions with parameter ν

for (op,ODE,RHS,growth) in ((:(hankelh1),"f^2*f'*D^2+(f*f'^2-f^2*f'')*D+(f^2-ν^2)*f'^3","0",:(imag)),
                            (:(hankelh2),"f^2*f'*D^2+(f*f'^2-f^2*f'')*D+(f^2-ν^2)*f'^3","0",:(imag)),
                            (:(besselj),"f^2*f'*D^2+(f*f'^2-f^2*f'')*D+(f^2-ν^2)*f'^3","0",:(imag)),
                            (:(bessely),"f^2*f'*D^2+(f*f'^2-f^2*f'')*D+(f^2-ν^2)*f'^3","0",:(imag)),
                            (:(besseli),"f^2*f'*D^2+(f*f'^2-f^2*f'')*D-(f^2+ν^2)*f'^3","0",:(real)),
                            (:(besselk),"f^2*f'*D^2+(f*f'^2-f^2*f'')*D-(f^2+ν^2)*f'^3","0",:(real)),
                            (:(besselkx),"f^2*f'*D^2+((-2f^2+f)*f'^2-f^2*f'')*D-(f+ν^2)*f'^3","0",:(real)),
                            (:(hankelh1x),"f^2*f'*D^2+((2im*f^2+f)*f'^2-f^2*f'')*D+(im*f-ν^2)*f'^3","0",:(imag)),
                            (:(hankelh2x),"f^2*f'*D^2+((-2im*f^2+f)*f'^2-f^2*f'')*D+(-im*f-ν^2)*f'^3","0",:(imag)))
    L,R = Meta.parse(ODE),Meta.parse(RHS)
    @eval begin
        function $op(ν,fin::Fun{S,T}) where {S<:Union{Ultraspherical,Chebyshev},T}
            f=setcanonicaldomain(fin)

            g=chop($growth(f),eps(T))
            xmin = isempty(g.coefficients) ? first(domain(g)) : argmin(g)
            xmax = isempty(g.coefficients) ? last(domain(g)) : argmax(g)
            opfxmin,opfxmax = $op(ν,f(xmin)),$op(ν,f(xmax))
            opmax = maximum(abs,(opfxmin,opfxmax))
            while opmax≤10eps(T) || abs(f(xmin)-f(xmax))≤10eps(T)
                xmin,xmax = rand(domain(f)),rand(domain(f))
                opfxmin,opfxmax = $op(ν,f(xmin)),$op(ν,f(xmax))
                opmax = maximum(abs,(opfxmin,opfxmax))
            end
            D=Derivative(space(f))
            B=[Evaluation(space(f),xmin),Evaluation(space(f),xmax)]
            u=([B;eval($L)]\[opfxmin/opmax;opfxmax/opmax;eval($R)/opmax])*opmax

            setdomain(u,domain(fin))
        end
    end
end

exp2(f::Fun) = exp(log(2)*f)
exp10(f::Fun) = exp(log(10)*f)
log2(f::Fun) = log(f)/log(2)
log10(f::Fun) = log(f)/log(10)

##TODO: the spacepromotion doesn't work for tan/tanh for a domain including zeros of cos/cosh inside.
tan(f::Fun) = sin(f)/cos(f) #This is inaccurate, but allows space promotion via division.
tanh(f::Fun) = sinh(f)/cosh(f) #This is inaccurate, but allows space promotion via division.

for (op,oprecip,opinv,opinvrecip) in ((:(sin),:(csc),:(asin),:(acsc)),
                                      (:(cos),:(sec),:(acos),:(asec)),
                                      (:(tan),:(cot),:(atan),:(acot)),
                                      (:(sinh),:(csch),:(asinh),:(acsch)),
                                      (:(cosh),:(sech),:(acosh),:(asech)),
                                      (:(tanh),:(coth),:(atanh),:(acoth)))
    @eval begin
        $oprecip(f::Fun) = 1/$op(f)
        $opinvrecip(f::Fun) = $opinv(1/f)
    end
end

rad2deg(f::Fun) = 180*f/π
deg2rad(f::Fun) = π*f/180

for (op,opd,opinv,opinvd) in ((:(sin),:(sind),:(asin),:(asind)),
                              (:(cos),:(cosd),:(acos),:(acosd)),
                              (:(tan),:(tand),:(atan),:(atand)),
                              (:(sec),:(secd),:(asec),:(asecd)),
                              (:(csc),:(cscd),:(acsc),:(acscd)),
                              (:(cot),:(cotd),:(acot),:(acotd)))
    @eval begin
        $opd(f::Fun) = $op(deg2rad(f))
        $opinvd(f::Fun) = rad2deg($opinv(f))
    end
end

#Won't get the zeros exactly 0 anyway so at least this way the length is smaller.
sinpi(f::Fun) = sin(π*f)
cospi(f::Fun) = cos(π*f)

function airy(k::Number,f::Fun)
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

besselh(ν,k::Integer,f::Fun) = k == 1 ? hankelh1(ν,f) : k == 2 ? hankelh2(ν,f) : throw(Base.Math.AmosException(1))

for jy in ("j","y"), ν in (0,1)
    bjy = Symbol(string("bessel",jy))
    bjynu = Meta.parse(string("SpecialFunctions.bessel",jy,ν))
    @eval begin
        $bjynu(f::Fun) = $bjy($ν,f)
    end
end

## Miscellaneous
for op in (:(expm1),:(log1p),:(lfact),:(sinc),:(cosc),
           :(erfinv),:(erfcinv),:(beta),:(lbeta),
           :(eta),:(zeta),:(gamma),:(lgamma),
           :(polygamma),:(invdigamma),:(digamma),:(trigamma))
    @eval begin
        $op(f::Fun{S,T}) where {S,T}=Fun($op ∘ f,domain(f))
    end
end


## <,≤,>,≥

for op in (:<,:>)
    @eval begin
        function $op(f::Fun,c::Number)
            if length(roots(f-c))==0
                $op(first(f),c)
            else
                false
            end
        end
        function $op(c::Number,f::Fun)
            if length(roots(f-c))==0
                $op(c,first(f))
            else
                false
            end
        end
    end
end



for op in (:(<=),:(>=))
    @eval begin
        function $op(f::Fun,c::Number)
            rts=roots(f-c)
            if length(rts)==0
                $op(first(f),c)
            elseif length(rts)==1
                if isapprox(rts[1],first(domain(f))) || isapprox(rts[1],last(domain(f)))
                    $op(f(fromcanonical(f,0.)),c)
                else
                    error("Implement for mid roots")
                end
            elseif length(rts)==2
                if isapprox(rts[1],first(domain(f))) && isapprox(rts[2],last(domain(f)))
                    $op(f(fromcanonical(f,0.)),c)
                else
                    error("Implement for mid roots")
                end
            else
                error("Implement for mid roots")
            end
        end
        function $op(c::Number,f::Fun)
            rts=sort(roots(f-c))
            if length(rts)==0
                $op(c,first(f))
            elseif length(rts)==1
                if isapprox(rts[1],first(domain(f))) || isapprox(rts[1],first(domain(f)))
                    $op(c,f(fromcanonical(f,0.)))
                else
                    error("Implement for mid roots")
                end
            elseif length(rts)==2
                if isapprox(rts[1],first(domain(f))) && isapprox(rts[2],first(domain(f)))
                    $op(c,f(fromcanonical(f,0.)))
                else
                    error("Implement for mid roots")
                end
            else
                error("Implement for mid roots")
            end
        end
    end
end




#TODO ≤,≥




## Piecewise Space

## PIecewiseSpace
# map over components

/(c::Number,f::Fun{S}) where {S<:PiecewiseSpace} = Fun(map(f->c/f,components(f)),PiecewiseSpace)
^(f::Fun{S},c::Integer) where {S<:PiecewiseSpace} = Fun(map(f->f^c,components(f)),PiecewiseSpace)
^(f::Fun{S},c::Number) where {S<:PiecewiseSpace} = Fun(map(f->f^c,components(f)),PiecewiseSpace)



for OP in (:abs,:sign,:log,:angle)
    @eval begin
        $OP(f::Fun{PiecewiseSpace{S,DD,RR},T}) where {S,DD,RR<:Real,T<:Real} =
            Fun(map($OP,components(f)),PiecewiseSpace)
        $OP(f::Fun{PiecewiseSpace{S,DD,RR}}) where {S,DD<:UnivariateDomain,RR} = Fun(map($OP,components(f)),PiecewiseSpace)
    end
end

# Return the locations of jump discontinuities
#
# Non Piecewise Spaces are assumed to have no jumps.
function jumplocations(f::Fun)
    eltype(domain(f))[]
end

# Return the locations of jump discontinuities
function jumplocations(f::Fun{S}) where{S<:PiecewiseSpace}
    d = domain(f)

    if ncomponents(d) < 2
      return eltype(domain(f))[]
    end

    dtol=10eps(eltype(d))
    ftol=10eps(cfstype(f))

    dc = components(d)
    fc = components(f)

    isjump = isapprox.(first.(dc[2:end]), last.(dc[1:end-1]), rtol=dtol) .&
           .!isapprox.(first.(fc[2:end]), last.(fc[1:end-1]), rtol=ftol)

    locs = last.(dc[1:end-1])
    locs[isjump]
end

#
# These formulæ, appearing in Eq. (2.5) of:
#
# A.-K. Kassam and L. N. Trefethen, Fourth-order time-stepping for stiff PDEs, SIAM J. Sci. Comput., 26:1214--1233, 2005,
#
# are derived to implement ETDRK4 in double precision without numerical instability from cancellation.
#

expα_asy(x) = (exp(x)*(4-3x+x^2)-4-x)/x^3
expβ_asy(x) = (exp(x)*(x-2)+x+2)/x^3
expγ_asy(x) = (exp(x)*(4-x)-4-3x-x^2)/x^3

# TODO: General types

expα_taylor(x::Union{Float64,ComplexF64}) = @evalpoly(x,1/6,1/6,3/40,1/45,5/1008,1/1120,7/51840,1/56700,1/492800,1/4790016,11/566092800,1/605404800,13/100590336000,1/106748928000,1/1580833013760,1/25009272288000,17/7155594141696000,1/7508956815360000)
expβ_taylor(x::Union{Float64,ComplexF64}) = @evalpoly(x,1/6,1/12,1/40,1/180,1/1008,1/6720,1/51840,1/453600,1/4435200,1/47900160,1/566092800,1/7264857600,1/100590336000,1/1494484992000,1/23712495206400,1/400148356608000,1/7155594141696000,1/135161222676480000)
expγ_taylor(x::Union{Float64,ComplexF64}) = @evalpoly(x,1/6,0/1,-1/120,-1/360,-1/1680,-1/10080,-1/72576,-1/604800,-1/5702400,-1/59875200,-1/691891200,-1/8717829120,-1/118879488000,-1/1743565824000,-1/27360571392000,-1/457312407552000,-1/8109673360588800)

expα(x::Float64) = abs(x) < 17/16 ? expα_taylor(x) : expα_asy(x)
expβ(x::Float64) = abs(x) < 19/16 ? expβ_taylor(x) : expβ_asy(x)
expγ(x::Float64) = abs(x) < 15/16 ? expγ_taylor(x) : expγ_asy(x)

expα(x::ComplexF64) = abs2(x) < (17/16)^2 ? expα_taylor(x) : expα_asy(x)
expβ(x::ComplexF64) = abs2(x) < (19/16)^2 ? expβ_taylor(x) : expβ_asy(x)
expγ(x::ComplexF64) = abs2(x) < (15/16)^2 ? expγ_taylor(x) : expγ_asy(x)

expα(x) = expα_asy(x)
expβ(x) = expβ_asy(x)
expγ(x) = expγ_asy(x)


for f in (:(exp),:(expm1),:expα,:expβ,:expγ)
    @eval $f(op::Operator) = OperatorFunction(op,$f)
end

## Special Multiplication
for f in (:+, :-, :*, :exp, :sin, :cos)
    @eval $f(M::Multiplication) = Multiplication($f(M.f), domainspace(M))
end

for f in (:+, :-, :*, :/, :\)
    @eval begin
        $f(M::Multiplication, c::Number) = Multiplication($f(M.f, c), domainspace(M))
        $f(c::Number, M::Multiplication) = Multiplication($f(c, M.f), domainspace(M))
    end
end

## ConstantSpace and PointSpace default overrides

for SP in (:ConstantSpace,:PointSpace)
    for OP in (:abs,:sign,:exp,:sqrt,:angle)
        @eval begin
            $OP(z::Fun{<:$SP,<:Complex}) = Fun(space(z),$OP.(coefficients(z)))
            $OP(z::Fun{<:$SP,<:Real}) = Fun(space(z),$OP.(coefficients(z)))
            $OP(z::Fun{<:$SP}) = Fun(space(z),$OP.(coefficients(z)))
        end
    end

    # we need to pad coefficients since 0^0 == 1
    for OP in (:^,)
        @eval begin
            function $OP(z::Fun{<:$SP},k::Integer)
                k ≠ 0 && return Fun(space(z),$OP.(coefficients(z),k))
                Fun(space(z),$OP.(pad(coefficients(z),dimension(space(z))),k))
            end
            function $OP(z::Fun{<:$SP},k::Number)
                k ≠ 0 && return Fun(space(z),$OP.(coefficients(z),k))
                Fun(space(z),$OP.(pad(coefficients(z),dimension(space(z))),k))
            end
        end
    end
end


for OP in (:(Base.max),:(Base.min))
    @eval begin
        $OP(a::Fun{CS1,T},b::Fun{CS2,V}) where {CS1<:ConstantSpace,CS2<:ConstantSpace,T<:Real,V<:Real} =
            Fun($OP(Number(a),Number(b)),space(a) ∪ space(b))
        $OP(a::Fun{CS,T},b::Real) where {CS<:ConstantSpace,T<:Real} =
            Fun($OP(Number(a),b),space(a))
        $OP(a::Real,b::Fun{CS,T}) where {CS<:ConstantSpace,T<:Real} =
            Fun($OP(a,Number(b)),space(b))
    end
end

for OP in (:<,:(Base.isless),:(<=),:>,:(>=))
    @eval begin
        $OP(a::Fun{CS},b::Fun{CS}) where {CS<:ConstantSpace} = $OP(convert(Number,a),Number(b))
        $OP(a::Fun{CS},b::Number) where {CS<:ConstantSpace} = $OP(convert(Number,a),b)
        $OP(a::Number,b::Fun{CS}) where {CS<:ConstantSpace} = $OP(a,convert(Number,b))
    end
end

# from DualNumbers
for (funsym, exp) in Calculus.symbolic_derivatives_1arg()
    funsym == :abs && continue
    funsym == :sign && continue
    funsym == :exp && continue
    funsym == :sqrt && continue
    @eval begin
        $(funsym)(z::Fun{CS,T}) where {CS<:ConstantSpace,T<:Real} =
            Fun($(funsym)(Number(z)),space(z))
        $(funsym)(z::Fun{CS,T}) where {CS<:ConstantSpace,T<:Complex} =
            Fun($(funsym)(Number(z)),space(z))
        $(funsym)(z::Fun{CS}) where {CS<:ConstantSpace} =
            Fun($(funsym)(Number(z)),space(z))
    end
end
