## abs

function splitmap(g,d::IntervalDomain,pts)
    if isempty(pts)
        Fun(g,d)
    elseif length(pts)==1 && (isapprox(first(pts),first(d))  ||  isapprox(last(pts),last(d)))
        Fun(g,d)
    elseif length(pts)==2 && isapprox(first(pts),first(d)) && isapprox(last(pts),last(d))
        Fun(g,d)
    else
        error("implement splitmap for "*string(typeof(d)))
    end
end

function splitmap(g,d::AffineDomain,pts)
    if isempty(pts)
        Fun(g,d)
    else
        da=first(d)
        isapprox(da,pts[1];atol=sqrt(eps(arclength(d)))) ? pts[1] = da : unshift!(pts,da)
        db=last(d)
        isapprox(db,pts[end];atol=sqrt(eps(arclength(d)))) ? pts[end] = db : push!(pts,db)
        Fun(g,pts)
    end
end

## TODO: type to ensure each dom in d.domains ::IntervalDomain
function splitmap(g,d::UnionDomain,pts)
    if isempty(pts)
        Fun(g,d)
    else
        dpts = ∂(d)
        pts = sort!(∪(pts,dpts),by=real)
        da=first(d)
        isapprox(da,pts[1];atol=sqrt(eps(mapreduce(arclength,+,d.domains)))) ? pts[1] = da : unshift!(pts,da)
        db=last(d)
        isapprox(db,pts[end];atol=sqrt(eps(mapreduce(arclength,+,d.domains)))) ? pts[end] = db : push!(pts,db)
        Fun(g,pts)
    end
end


function splitatroots(f::Fun)
    d=domain(f)
    pts=union(roots(f)) # union removes multiplicities
    if isempty(pts)
        f# : splitmap(f,d,pts)
    else
        da=first(d)
        isapprox(da,pts[1]) ? pts[1] = da : unshift!(pts,da)
        db=last(d)
        isapprox(db,pts[end]) ? pts[end] = db : push!(pts,db)
        Fun(x->f(x),pts)
    end
end

function Base.abs{S<:RealUnivariateSpace,T<:Real}(f::Fun{S,T})
    d=domain(f)

    pts=roots(f)

    splitmap(x->abs(f(x)),d,pts)
end

function Base.abs(f::Fun)
    d=domain(f)

    pts=roots(f)

    if isempty(pts)
        # This makes sure Laurent returns real type
        real(Fun(x->abs(f(x)),space(f)))
    else
        splitmap(x->abs(f(x)),d,pts)
    end
end



function Base.sign{S<:RealUnivariateSpace,T<:Real}(f::Fun{S,T})
    d=domain(f)

    pts=roots(f)

    if isempty(pts)
        sign(first(f))*one(f)
    else
        @assert isa(d,AffineDomain)
        da=first(d)
        isapprox(da,pts[1];atol=sqrt(eps(arclength(d)))) ? pts[1] = da : unshift!(pts,da)
        db=last(d)
        isapprox(db,pts[end];atol=sqrt(eps(arclength(d)))) ? pts[end] = db : push!(pts,db)
        midpts = .5(pts[1:end-1]+pts[2:end])
        Fun(sign(f(midpts)),pts)
    end
end

for op in (:(Base.max),:(Base.min))
    @eval begin
        function $op{S<:RealUnivariateSpace,V<:RealUnivariateSpace,T<:Real}(f::Fun{S,T},g::Fun{V,T})
            h=f-g
            d=domain(h)
            pts=roots(h)
            splitmap(x->$op(f(x),g(x)),d,pts)
        end
        $op{S<:RealUnivariateSpace,T<:Real}(f::Fun{S,T},g::Real) = $op(f,Fun(g,domain(f)))
        $op{S<:RealUnivariateSpace,T<:Real}(f::Real,g::Fun{S,T}) = $op(Fun(f,domain(g)),g)
    end
end

# division by fun

function ./(c::Fun,f::Fun)
    d=domain(f)
    @assert domain(c)==d
    cd=canonicaldomain(d)
    if typeof(d)!=typeof(cd)
        # project first to simplify
        return setdomain(setdomain(c,cd)./setdomain(f,cd),d)
    end

    r=roots(f)
    tol=10eps()
    if length(r)==0 || norm(c(r))<tol
        linsolve(Multiplication(f,space(c)),c;tolerance=tol)
    else
        c.*(1./f)
    end
end

function ./(c::Number,f::Fun)
    r=roots(f)
    tol=10eps()
    @assert length(r)==0
    linsolve(Multiplication(f,space(f)),c*ones(space(f));tolerance=tol)
end

# project to interval if we are not on the interview
# TODO: need to work out how to set piecewise domain


scaleshiftdomain(f::Fun,sc,sh)=setdomain(f,sc*domain(f)+sh)

./{λ,DD}(c::Number,f::Fun{Ultraspherical{λ,DD}}) = c./Fun(f,Chebyshev(domain(f)))
./{DD}(c::Number,f::Fun{Jacobi{DD}}) = c./Fun(f,Chebyshev(domain(f)))

./{C<:Chebyshev}(c::Number,f::Fun{C})=setdomain(c./setcanonicaldomain(f),domain(f))
function ./{DD<:Interval}(c::Number,f::Fun{Chebyshev{DD}})
    fc = setcanonicaldomain(f)
    d=domain(f)
    # if domain f is small then the pts get projected in
    tol = 200eps()*norm(f.coefficients,1)

    # we prune out roots at the boundary first
    if ncoefficients(f)==1
        return Fun(c/f.coefficients[1],space(f))
    elseif ncoefficients(f)==2
        if isempty(roots(f))
            return linsolve(Multiplication(f,space(f)),c;tolerance=tol)
        elseif isapprox(fc.coefficients[1],fc.coefficients[2])
            # we check directly for f*(1+x)
            return Fun([c./fc.coefficients[1]],JacobiWeight(-1,0,space(f)))
        elseif isapprox(fc.coefficients[1],-fc.coefficients[2])
            # we check directly for f*(1-x)
            return Fun([c./fc.coefficients[1]],JacobiWeight(0,-1,space(f)))
        else
            # we need to split at the only root
            return c./splitatroots(f)
        end
    elseif abs(first(fc))≤tol
        #left root
        g=divide_singularity((1,0),fc)
        p=c./g
        x=identity_fun(domain(p))
        return scaleshiftdomain(p/(1+x),(d.b - d.a)/2,(d.a + d.b)/2 )
    elseif abs(last(fc))≤tol
        #right root
        g=divide_singularity((0,1),fc)
        p=c./g
        x=identity_fun(domain(p))
        return scaleshiftdomain(p/(1-x),(d.b - d.a)/2,(d.a + d.b)/2 )
    else
        r = roots(fc)

        if length(r) == 0
            return linsolve(Multiplication(f,space(f)),c;tolerance=tol)
        elseif abs(last(r)+1.0)≤tol  # double check
            #left root
            g=divide_singularity((1,0),fc)
            p=c./g
            x=identity_fun(domain(p))
            return scaleshiftdomain(p/(1+x),(d.b - d.a)/2,(d.a + d.b)/2 )
        elseif abs(last(r)-1.0)≤tol  # double check
            #right root
            g=divide_singularity((0,1),fc)
            p=c./g
            x=identity_fun(domain(p))
            return scaleshiftdomain(p/(1-x),(d.b - d.a)/2,(d.a + d.b)/2 )
        else
            # no roots on the boundary
            return c./splitatroots(f)
        end
    end
end



# function ./{S<:MappedSpace}(c::Number,f::Fun{S})
#     g=c./Fun(coefficients(f),space(f).space)
#     Fun(coefficients(g),MappedSpace(domain(f),space(g)))
# end
# function .^{S<:Space,D,T}(f::Fun{MappedSpace{S,D,T}},k::Float64)
#     g=Fun(coefficients(f),space(f).space).^k
#     Fun(coefficients(g),MappedSpace(domain(f),space(g)))
# end
#
#
# #TODO: Unify following
# function .^{S<:Chebyshev,D,T}(f::Fun{MappedSpace{S,D,T}},k::Float64)
#     sp=space(f)
#     # Need to think what to do if this is ever not the case..
#     @assert isapprox(domain(sp.space),Interval())
#     fc = Fun(f.coefficients,sp.space) #Project to interval
#
#     r = sort(roots(fc))
#     @assert length(r) <= 2
#
#     if length(r) == 0
#         Fun(Fun(x->fc(x)^k).coefficients,sp)
#     elseif length(r) == 1
#         @assert isapprox(abs(r[1]),1)
#
#         if isapprox(r[1],1.)
#             Fun(coefficients(divide_singularity(true,fc)^k),MappedSpace(sp.domain,JacobiWeight(0.,k,sp.space)))
#         else
#             Fun(coefficients(divide_singularity(false,fc)^k),MappedSpace(sp.domain,JacobiWeight(k,0.,sp.space)))
#         end
#     else
#         @assert isapprox(r[1],-1)
#         @assert isapprox(r[2],1)
#
#         Fun(coefficients(divide_singularity(fc)^k),MappedSpace(sp.domain,JacobiWeight(k,k,sp.space)))
#     end
# end

function .^{C<:Chebyshev}(f::Fun{C},k::Float64)
    # Need to think what to do if this is ever not the case..
    sp = space(f)
    fc = setdomain(f,Interval()) #Project to interval

    r = sort(roots(fc))
    #TODO divideatroots
    @assert length(r) <= 2

    if length(r) == 0
        Fun(Fun(x->fc(x)^k).coefficients,sp)
    elseif length(r) == 1
        @assert isapprox(abs(r[1]),1)

        if isapprox(r[1],1.)
            Fun(coefficients(divide_singularity(true,fc)^k),JacobiWeight(0.,k,sp))
        else
            Fun(coefficients(divide_singularity(false,fc)^k),JacobiWeight(k,0.,sp))
        end
    else
        @assert isapprox(r[1],-1)
        @assert isapprox(r[2],1)

        Fun(coefficients(divide_singularity(fc)^k),JacobiWeight(k,k,sp))
    end
end

#TODO: implement
.^(f::Fun{Jacobi},k::Float64)=Fun(f,Chebyshev).^k

# Default is just try solving ODE
function .^{S,T}(f::Fun{S,T},β)
    A=Derivative()-β*differentiate(f)/f
    B=Evaluation(first(domain(f)))
    [B,A]\first(f)^β
end

Base.sqrt{S,T}(f::Fun{S,T})=f^0.5
Base.cbrt{S,T}(f::Fun{S,T})=f^(1/3)

## We use \ as the Fun constructor might miss isolated features

## First order functions


Base.log(f::Fun)=cumsum(differentiate(f)/f)+log(first(f))

# function Base.log{MS<:MappedSpace}(f::Fun{MS})
#     g=log(Fun(f.coefficients,space(f).space))
#     Fun(g.coefficients,MappedSpace(domain(f),space(g)))
# end

# project first to [-1,1] to avoid issues with
# complex derivative
function Base.log{US<:Ultraspherical}(f::Fun{US})
    if domain(f)==Interval()
        r = sort(roots(f))
        #TODO divideatroots
        @assert length(r) <= 2

        if length(r) == 0
            cumsum(differentiate(f)/f)+log(first(f))
        elseif length(r) == 1
            @assert isapprox(abs(r[1]),1)

            if isapprox(r[1],1.)
                g=divide_singularity(true,f)
                lg=Fun([1.],LogWeight(0.,1.,Chebyshev()))
                if isapprox(g,1.)  # this means log(g)~0
                    lg
                else # log((1-x)) + log(g)
                    lg⊕log(g)
                end
            else
                g=divide_singularity(false,f)
                lg=Fun([1.],LogWeight(1.,0.,Chebyshev()))
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
            lg=Fun([1.],LogWeight(1.,1.,Chebyshev()))
            if isapprox(g,1.)  # this means log(g)~0
                lg
            else # log((1+x)) + log(g)
                lg⊕log(g)
            end
        end
    else
        # this makes sure differentiate doesn't
        # make the function complex
        g=log(setdomain(f,Interval()))
        setdomain(g,domain(f))
    end
end


function Base.log{T<:Real,D}(f::Fun{Fourier{D},T})
    if isreal(domain(f))
        cumsum(differentiate(f)/f)+log(first(f))
    else
        # this makes sure differentiate doesn't
        # make the function complex
        g=log(setdomain(f,PeriodicInterval()))
        setdomain(g,domain(f))
    end
end


Base.atan(f::Fun)=cumsum(f'/(1+f^2))+atan(first(f))


# this is used to find a point in which to impose a boundary
# condition in calculating secial functions
function specialfunctionnormalizationpoint(op,growth,f)
    g=chop(growth(f),eps(eltype(f)))
    xmin=g.coefficients==[0.]?first(domain(g)):indmin(g)
    xmax=g.coefficients==[0.]?last(domain(g)):indmax(g)
    opfxmin,opfxmax = op(f(xmin)),op(f(xmax))
    opmax = maxabs((opfxmin,opfxmax))
    if abs(opfxmin) == opmax xmax,opfxmax = xmin,opfxmin end
    xmax,opfxmax,opmax
end



# ODE gives the first order ODE a special function op satisfies,
# RHS is the right hand side
# growth says what to use to choose a good point to impose an initial condition
for (op,ODE,RHS,growth) in ((:(Base.exp),"D-f'","0",:(real)),
                            (:(Base.asinh),"sqrt(f^2+1)*D","f'",:(real)),
                            (:(Base.acosh),"sqrt(f^2-1)*D","f'",:(real)),
                            (:(Base.atanh),"(1-f^2)*D","f'",:(real)),
                            (:(Base.erfcx),"D-2f*f'","-2f'/sqrt(π)",:(real)),
                            (:(Base.dawson),"D+2f*f'","f'",:(real)))
    L,R = parse(ODE),parse(RHS)
    @eval begin
        # depice before doing op
        $op{PW<:PiecewiseSpace}(f::Fun{PW})=depiece(map(f->$op(f),pieces(f)))

        # We remove the MappedSpace
        # function $op{MS<:MappedSpace}(f::Fun{MS})
        #     g=exp(Fun(f.coefficients,space(f).space))
        #     Fun(g.coefficients,MappedSpace(domain(f),space(g)))
        # end
        function $op{S,T}(f::Fun{S,T})
            xmax,opfxmax,opmax=specialfunctionnormalizationpoint($op,$growth,f)
            # we will assume the result should be smooth on the domain
            # even if f is not
            # This supports Line/Rays
            D=Derivative(domain(f))
            B=Evaluation(domainspace(D),xmax)
            #([B,eval($L)]\[opfxmax/opmax,eval($R)/opmax])*opmax
            linsolve([B,eval($L)],Any[opfxmax/opmax,eval($R)/opmax];tolerance=eps(T))*opmax
        end
    end
end

# JacobiWeight explodes, we want to ensure the solution incorporates the fact
# that exp decays rapidly
function Base.exp{JW<:JacobiWeight}(f::Fun{JW})
    if !isa(domain(f),Interval)
        # project first to get better derivative behaviour
        return setdomain(exp(setdomain(f,Interval())),domain(f))
    end

    S=space(f)
    q=Fun(f.coefficients,S.space)
    if isapprox(S.α,0.) && isapprox(S.β,0.)
        exp(q)
    elseif S.α < 0 && isapprox(first(q),0.)
        # this case can remove the exponential decay
        exp(Fun(f,JacobiWeight(S.α+1,S.β,S.space)))
    elseif S.β < 0 && isapprox(last(q),0.)
        exp(Fun(f,JacobiWeight(S.α,S.β+1,S.space)))
    elseif S.α > 0 && isapproxinteger(S.α)
        exp(Fun(f,JacobiWeight(0.,S.β,S.space)))
    elseif S.β > 0 && isapproxinteger(S.β)
        exp(Fun(f,JacobiWeight(S.α,0.,S.space)))
    else
        #find normalization point
        xmax,opfxmax,opmax=specialfunctionnormalizationpoint(exp,real,f)

        if S.α < 0 && S.β < 0
            # provided both are negative, we get exponential decay on both ends
            @assert real(first(q)) < 0 && real(last(q)) < 0
            s=JacobiWeight(2.,2.,domain(f))
        elseif S.α < 0 && isapprox(S.β,0.)
            @assert real(first(q)) < 0
            s=JacobiWeight(2.,0.,domain(f))
        elseif S.β < 0 && isapprox(S.α,0.)
            @assert real(last(q)) < 0
            s=JacobiWeight(0.,2.,domain(f))
        else
            error("exponential of fractional power, not implemented")
        end

        D=Derivative(s)
        B=Evaluation(s,xmax)

        linsolve([B,D-f'],Any[opfxmax/opmax,0.];tolerance=eps(eltype(f)))*opmax
    end
end







Base.acos(f::Fun)=cumsum(-f'/sqrt(1-f^2))+acos(first(f))
Base.asin(f::Fun)=cumsum(f'/sqrt(1-f^2))+asin(first(f))



## Second order functions


Base.sin{S<:Space{RealBasis},T<:Real}(f::Fun{S,T}) = imag(exp(im*f))
Base.cos{S<:Space{RealBasis},T<:Real}(f::Fun{S,T}) = real(exp(im*f))
Base.sin{S<:Ultraspherical,T<:Real}(f::Fun{S,T}) = imag(exp(im*f))
Base.cos{S<:Ultraspherical,T<:Real}(f::Fun{S,T}) = real(exp(im*f))



for (op,ODE,RHS,growth) in ((:(Base.erf),"f'*D^2+(2f*f'^2-f'')*D","0",:(imag)),
                            (:(Base.erfi),"f'*D^2-(2f*f'^2+f'')*D","0",:(real)),
                            (:(Base.sin),"f'*D^2-f''*D+f'^3","0",:(imag)),
                            (:(Base.cos),"f'*D^2-f''*D+f'^3","0",:(imag)),
                            (:(Base.sinh),"f'*D^2-f''*D-f'^3","0",:(real)),
                            (:(Base.cosh),"f'*D^2-f''*D-f'^3","0",:(real)),
                            (:(Base.airyai),"f'*D^2-f''*D-f*f'^3","0",:(imag)),
                            (:(Base.airybi),"f'*D^2-f''*D-f*f'^3","0",:(imag)),
                            (:(Base.airyaiprime),"f'*D^2-f''*D-f*f'^3","airyai(f)*f'^3",:(imag)),
                            (:(Base.airybiprime),"f'*D^2-f''*D-f*f'^3","airybi(f)*f'^3",:(imag)))
    L,R = parse(ODE),parse(RHS)
    @eval begin
        function $op{S,T}(f::Fun{S,T})
            g=chop($growth(f),eps(T))
            xmin=g.coefficients==[0.]?first(domain(g)):indmin(g)
            xmax=g.coefficients==[0.]?last(domain(g)):indmax(g)
            opfxmin,opfxmax = $op(f(xmin)),$op(f(xmax))
            opmax = maxabs((opfxmin,opfxmax))
            while opmax≤10eps(T) || abs(f(xmin)-f(xmax))≤10eps(T)
                xmin,xmax = rand(domain(f)),rand(domain(f))
                opfxmin,opfxmax = $op(f(xmin)),$op(f(xmax))
                opmax = maxabs((opfxmin,opfxmax))
            end
            D=Derivative(space(f))
            B=[Evaluation(space(f),xmin),Evaluation(space(f),xmax)]
            linsolve([B;eval($L)],[opfxmin/opmax;opfxmax/opmax;eval($R)/opmax];
                        tolerance=10ncoefficients(f)*eps(T)*opmax)*opmax
        end
    end
end

Base.erfc(f::Fun)=1-erf(f)

## Second order functions with parameter ν

for (op,ODE,RHS,growth) in ((:(Base.hankelh1),"f^2*f'*D^2+(f*f'^2-f^2*f'')*D+(f^2-ν^2)*f'^3","0",:(imag)),
                            (:(Base.hankelh2),"f^2*f'*D^2+(f*f'^2-f^2*f'')*D+(f^2-ν^2)*f'^3","0",:(imag)),
                            (:(Base.besselj),"f^2*f'*D^2+(f*f'^2-f^2*f'')*D+(f^2-ν^2)*f'^3","0",:(imag)),
                            (:(Base.bessely),"f^2*f'*D^2+(f*f'^2-f^2*f'')*D+(f^2-ν^2)*f'^3","0",:(imag)),
                            (:(Base.besseli),"f^2*f'*D^2+(f*f'^2-f^2*f'')*D-(f^2+ν^2)*f'^3","0",:(real)),
                            (:(Base.besselk),"f^2*f'*D^2+(f*f'^2-f^2*f'')*D-(f^2+ν^2)*f'^3","0",:(real)),
                            (:(Base.besselkx),"f^2*f'*D^2+((-2f^2+f)*f'^2-f^2*f'')*D-(f+ν^2)*f'^3","0",:(real)),
                            (:(Base.hankelh1x),"f^2*f'*D^2+((2im*f^2+f)*f'^2-f^2*f'')*D+(im*f-ν^2)*f'^3","0",:(imag)),
                            (:(Base.hankelh2x),"f^2*f'*D^2+((-2im*f^2+f)*f'^2-f^2*f'')*D+(-im*f-ν^2)*f'^3","0",:(imag)))
    L,R = parse(ODE),parse(RHS)
    @eval begin
        function $op{S<:Ultraspherical,T}(ν,f::Fun{S,T})
            g=chop($growth(f),eps(T))
            xmin=g.coefficients==[0.]?first(domain(g)):indmin(g)
            xmax=g.coefficients==[0.]?last(domain(g)):indmax(g)
            opfxmin,opfxmax = $op(ν,f(xmin)),$op(ν,f(xmax))
            opmax = maxabs((opfxmin,opfxmax))
            while opmax≤10eps(T) || abs(f(xmin)-f(xmax))≤10eps(T)
                xmin,xmax = rand(domain(f)),rand(domain(f))
                opfxmin,opfxmax = $op(ν,f(xmin)),$op(ν,f(xmax))
                opmax = maxabs((opfxmin,opfxmax))
            end
            D=Derivative(space(f))
            B=[Evaluation(space(f),xmin),Evaluation(space(f),xmax)]
            ([B;eval($L)]\[opfxmin/opmax;opfxmax/opmax;eval($R)/opmax])*opmax
        end
    end
end

Base.exp2(f::Fun) = exp(log(2)*f)
Base.exp10(f::Fun) = exp(log(10)*f)
Base.log2(f::Fun) = log(f)/log(2)
Base.log10(f::Fun) = log(f)/log(10)

##TODO: the spacepromotion doesn't work for tan/tanh for a domain including zeros of cos/cosh inside.
Base.tan(f::Fun) = sin(f)/cos(f) #This is inaccurate, but allows space promotion via division.
Base.tanh(f::Fun) = sinh(f)/cosh(f) #This is inaccurate, but allows space promotion via division.

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

Base.besselh(ν,k::Integer,f::Fun) = k == 1 ? hankelh1(ν,f) : k == 2 ? hankelh2(ν,f) : throw(Base.Math.AmosException(1))

for jy in ("j","y"), ν in (0,1)
    bjy = Symbol(string("bessel",jy))
    bjynu = parse(string("Base.bessel",jy,ν))
    @eval begin
        $bjynu(f::Fun) = $bjy($ν,f)
    end
end

## Miscellaneous
for op in (:(Base.expm1),:(Base.log1p),:(Base.lfact),:(Base.sinc),:(Base.cosc),
           :(Base.erfinv),:(Base.erfcinv),:(Base.beta),:(Base.lbeta),
           :(Base.eta),:(Base.zeta),:(Base.gamma),:(Base.lgamma),
           :(Base.polygamma),:(Base.invdigamma),:(Base.digamma),:(Base.trigamma))
    @eval begin
        $op{S,T}(f::Fun{S,T})=Fun(x->$op(f(x)),domain(f))
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
# map over pieces

./{S<:PiecewiseSpace}(c::Number,f::Fun{S})=depiece(mapreduce(f->c./f,vcat,pieces(f)))
.^{S<:PiecewiseSpace}(f::Fun{S},c::Integer)=depiece(mapreduce(f->f.^c,vcat,pieces(f)))
.^{S<:PiecewiseSpace}(f::Fun{S},c::Number)=depiece(mapreduce(f->f.^c,vcat,pieces(f)))



for OP in (:(Base.abs),:(Base.sign))
    @eval begin
        $OP{S,DD,T<:Real}(f::Fun{PiecewiseSpace{S,RealBasis,DD,1},T})=depiece(mapreduce($OP,vcat,pieces(f)))
        $OP{S,DD,B}(f::Fun{PiecewiseSpace{S,B,DD,1}})=depiece(mapreduce($OP,vcat,pieces(f)))
    end
end

## PointSpace

for OP in (:(Base.abs),:(Base.sign))
    # ambiguity warnings
    @eval $OP{S<:PointSpace,T<:Real}(f::Fun{S,T})=Fun(map($OP,f.coefficients),space(f))
end
for OP in (:(Base.exp),:(Base.abs),:(Base.sign))
    @eval $OP{S<:PointSpace}(f::Fun{S})=Fun(map($OP,f.coefficients),space(f))
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

expα_taylor(x::Union{Float64,Complex128}) = @evalpoly(x,1/6,1/6,3/40,1/45,5/1008,1/1120,7/51840,1/56700,1/492800,1/4790016,11/566092800,1/605404800,13/100590336000,1/106748928000,1/1580833013760,1/25009272288000,17/7155594141696000,1/7508956815360000)
expβ_taylor(x::Union{Float64,Complex128}) = @evalpoly(x,1/6,1/12,1/40,1/180,1/1008,1/6720,1/51840,1/453600,1/4435200,1/47900160,1/566092800,1/7264857600,1/100590336000,1/1494484992000,1/23712495206400,1/400148356608000,1/7155594141696000,1/135161222676480000)
expγ_taylor(x::Union{Float64,Complex128}) = @evalpoly(x,1/6,0/1,-1/120,-1/360,-1/1680,-1/10080,-1/72576,-1/604800,-1/5702400,-1/59875200,-1/691891200,-1/8717829120,-1/118879488000,-1/1743565824000,-1/27360571392000,-1/457312407552000,-1/8109673360588800)

expα(x::Float64) = abs(x) < 17/16 ? expα_taylor(x) : expα_asy(x)
expβ(x::Float64) = abs(x) < 19/16 ? expβ_taylor(x) : expβ_asy(x)
expγ(x::Float64) = abs(x) < 15/16 ? expγ_taylor(x) : expγ_asy(x)

expα(x::Complex128) = abs2(x) < (17/16)^2 ? expα_taylor(x) : expα_asy(x)
expβ(x::Complex128) = abs2(x) < (19/16)^2 ? expβ_taylor(x) : expβ_asy(x)
expγ(x::Complex128) = abs2(x) < (15/16)^2 ? expγ_taylor(x) : expγ_asy(x)

expα(x) = expα_asy(x)
expβ(x) = expβ_asy(x)
expγ(x) = expγ_asy(x)

@vectorize_1arg Number expα
@vectorize_1arg Number expβ
@vectorize_1arg Number expγ

for f in (:(Base.exp),:(Base.expm1),:expα,:expβ,:expγ)
    @eval begin
        $f(op::Operator) = OperatorFunction(op,$f)
    end
end
