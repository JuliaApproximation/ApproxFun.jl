export newton, linop

struct DualFun{F,T}
    f::F
    J::T
end
DualFun(f::Fun) = DualFun(f,Operator(I,space(f)))


domain(df::DualFun) = domain(df.f)

differentiate(d::DualFun) = DualFun(d.f',Derivative(rangespace(d.J))*d.J)
integrate(d::DualFun) = DualFun(integrate(d.f),Integral(rangespace(d.J))*d.J)
function cumsum(d::DualFun)
    Q = Integral(rangespace(d.J))*d.J
    rs = rangespace(Q)
    DualFun(cumsum(d.f),Q-ones(rs)*Evaluation(rs,leftendpoint)*Q)
end

# total definite integral of a function
function sum(df::DualFun)
    Q = Integral(rangespace(df.J))*df.J
    return DualFun(sum(df.f), Evaluation(rangespace(Q), rightendpoint)*Q)
end


adjoint(d::DualFun) = differentiate(d)

^(d::DualFun, k::Number) = DualFun(d.f^k,k*d.f^(k-1)*d.J)

# from DualNumbers
for (funsym, exp) in Calculus.symbolic_derivatives_1arg()
    @eval function $(funsym)(z::DualFun)
        x = z.f
        xp = z.J
        DualFun($(funsym)(x),$exp*xp)
    end
end

for OP in (:+,:-)
    @eval begin
        $OP(a::DualFun,b::Union{Number,Fun}) = DualFun($OP(a.f,b),a.J)
        $OP(a::Union{Number,Fun},b::DualFun) = DualFun($OP(a,b.f),$OP(b.J))
        $OP(a::DualFun,b::DualFun) = DualFun($OP(a.f,b.f),$OP(a.J,b.J))
    end
end
-(a::DualFun)=DualFun(-a.f,-a.J)

*(a::Union{Number,Fun},b::DualFun)=DualFun(a*b.f,a*b.J)
*(a::DualFun,b::Union{Number,Fun})=DualFun(b*a.f,b*a.J)
*(a::DualFun,b::DualFun)=DualFun(a.f*b.f,a.f*b.J+b.f*a.J)

/(a::Union{Number,Fun},b::DualFun)=DualFun(a/b.f,-a/b.f^2*b.J)
/(a::DualFun,b::Union{Number,Fun})=DualFun(a.f/b,a.J/b)
/(a::DualFun,b::DualFun)=DualFun(a.f/b.f,a.J/b.f-a.f/b.f^2*b.J)


(d::DualFun)(x) = DualFun(d.f(x),Evaluation(rangespace(d.J),x)*d.J)
first(d::DualFun) = DualFun(first(d.f),Evaluation(rangespace(d.J),leftendpoint)*d.J)
last(d::DualFun) = DualFun(last(d.f),Evaluation(rangespace(d.J),rightendpoint)*d.J)

jacobian(d::DualFun) = d.J
jacobian(a::Number) = zero(a)
jacobian(f::Fun) = Operator(I,space(f))

promote_rule(::Type{DF},::Type{T}) where {DF<:DualFun,T<:Number}=DualFun
convert(::Type{DualFun},b::Number) = DualFun(b,0)



function linop(f::Function,ds::Space)
    if (isgeneric(f) && applicable(f,0)) || (!isgeneric(f)&&arglength(f)==1)
        df=f(DualFun(zeros(ds)))
    elseif (isgeneric(f) && applicable(f,0,0)) || (!isgeneric(f)&&arglength(f)==2)
        df=f(Fun(ds),DualFun(zeros(ds)))
    else
        error("Not implemented")
    end

    if isa(df,Array)
        map(u->u.J,df)
    else
        df.J
    end
end

linop(f::Function,d) = linop(f,Space(d))
linop(f::Function) = linop(f,Chebyshev())  #TODO: UnsetSpace



# full operator should be
# N=u->[B*u-bcs;...]
function newton(N,u0::Fun;maxiterations=15,tolerance=1E-15)
    u=u0
    err=Inf
    for k=1:maxiterations
        DF=N(DualFun(u))
        J=map(jacobian,DF)
        F=map(d->d.f,DF)
        unew=u-J\F
        err=norm(unew-u)
        if err≤10tolerance
            return unew
        else
            u=chop(unew,tolerance)
        end
    end
    @warn "Maximum number of iterations $maxiterations reached, with approximate accuracy of $err."
    return u
end

# Newton iteration for a system of nonlinear equations
# INPUTS:
#   N - a nonlinear function N(u1, u2, u3, ...) = [BCs;...]
#   u0 - an array of initial Funs
#       NOTE: these initial Funs should be "chopped" as much as possible for optimal speed
#
# KSS and ZJS
function newton(N, u0::Vector{<:Fun}; maxiterations=15, tolerance=1E-15)
    u = copy(u0)
    maxerr::Float64 = Inf

    numf = length(u0)   # number of functions

    if numf == 0
        error("u0 must contain at least 1 function")
    end

    Js = Vector{Any}(undef,numf) # jacobians

    for k = 1:maxiterations
        # ------ calculate Jacobian for each function - O(numf^2) ----

        # use the first call to populate the value of the functions
        # this saves one call to N()
        F = N(DualFun(u[1]), u[2:numf]...)
        # F = N([ (j == 1 ? DualFun(u[j]) : u[j]) for j = 1:numf ]...)
        Js[1] = jacobian.(F)

        for i = 2:numf
            Js[i] = jacobian.(N([ (j == i ? DualFun(u[j]) : u[j]) for j = 1:numf ]...))
        end

        F = N(u...)
        J = Base.typed_hcat(Operator, Js...)

        # # ------- update step ---------
        unew::typeof(u) = u .- J\F
        maxerr = maximum(((x,y),) -> norm(x-y), zip(unew, u)) # use the max norm error as the overall error

        if maxerr ≤ 10*tolerance
            return unew
        else
            map!(q -> chop!(q, tolerance), u, unew)
        end
    end

    @warn "Maximum number of iterations $maxiterations reached, with approximate accuracy of $maxerr."

    return u
end

newton(N,d::Domain;opts...) =
    newton(N,zeros(d);opts...)

newton(N,d;opts...) =
    newton(N,Domain(d);opts...)
