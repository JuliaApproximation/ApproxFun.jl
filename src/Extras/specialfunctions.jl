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





## We use \ as the Fun constructor might miss isolated features
function Base.exp{S<:Ultraspherical}(f::Fun{S})
    rf=chop(real(f),eps())

    xm=rf.coefficients[1]==[0.]?first(domain(rf)):indmax(rf)    


    B=Evaluation(space(f),xm)
    D=Derivative(space(f))
    A=[B,D-differentiate(f)]
    A\[exp(f[xm]),0.]    
end

## Less accurate than solving differential equation with \
for op in (:(Base.cos),:(Base.sin),:(Base.cospi),:(Base.sinpi),:(Base.sinc),:(Base.exp))
    @eval begin
        $op{S,T}(f::Fun{S,T})=Fun(x->($op)(f[x]),space(f))
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





## The following code will fail when xp and xm correspond to zeros of sin
# for op = (:(Base.cos),:(Base.sin))
#     @eval begin
#         function ($op)(f::Fun)
#             d=domain(f)
#             D=diff(d)
#             f2=f.*f
#             xp=indmax(f2);xm=indmin(f2);
#             B=[Evaluation(d,xm),Evaluation(d,xp)]
#             fp=diff(f)
#             fpp=diff(fp)
#             [B,fp*D^2 - fpp*D + fp.^3]\[($op)(f[xm]),($op)(f[xp])]
#         end
#     end
# end

## The following backslash code works for real arguments but fails for complex Funs.


Base.cos{S<:FunctionSpace{Float64},T<:Real}(f::Fun{S,T})=real(exp(im*f))
Base.sin{S<:FunctionSpace{Float64},T<:Real}(f::Fun{S,T})=imag(exp(im*f))

function pochhammer{T<:Number}(x::T,n::Integer)
    ret = one(T)
    if n>=0
        for i=0:n-1
            ret *= x+i
        end
    else
        ret /= pochhammer(x+n,-n)
    end
    ret
end
function pochhammer{T<:Number}(x::Array{T},n::Integer)
    ret = ones(T,size(x))
    if n>=0
        for i=0:n-1
            ret .*= x+i
        end
    else
        ret ./= pochhammer(x+n,-n)
    end
    ret
end
pochhammer{T1<:Number,T2<:Number}(x::T1,n::T2) = gamma(x+n)/gamma(x)
pochhammer{T1<:Number,T2<:Number}(x::Array{T1},n::T2) = gamma(x+n)./gamma(x)