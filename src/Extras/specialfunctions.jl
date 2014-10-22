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


function ./(c::Number,f::Fun{ChebyshevSpace})
    fc = Fun(canonicalcoefficients(f),Interval())
    r = roots(fc)
    x = Fun(identity)
    
    tol = 10eps()
    
    @assert length(r) <= 2
    
    if length(r) == 0
        linsolve(Multiplication(f,space(f)),c;tolerance=tol)
    elseif length(r) == 1
        if abs(abs(r[1]) - 1.) < tol
            if sign(r[1]) < 0
                g = linsolve(Multiplication(x+1,space(fc)),fc;tolerance=10eps())                 
                Fun(canonicalcoefficients(c./g),JacobiWeightSpace(-1,0,domain(f)))
            else
                g = linsolve(Multiplication(1-x,space(fc)),fc;tolerance=10eps())                              
                Fun(canonicalcoefficients(c./g),JacobiWeightSpace(0,-1,domain(f)))                
            end 
        else
            error("need to implement splitting")
        end
    else
        @assert abs(r[1]+1) < tol
        @assert abs(r[2]-1) < tol                        
        g = linsolve(Multiplication(1-x.^2,space(fc)),fc;tolerance=10eps()) 
        # divide out singularities, tolerance needs to be chosen since we don't get
        # spectral convergence
        # TODO: switch to dirichlet basis
        Fun(canonicalcoefficients(c./g),JacobiWeightSpace(-1,-1,domain(f)))  
    end
end





## We use \ as the Fun constructor might miss isolated features
function Base.exp(f::Fun)
    rf=chop(real(f),eps())

    xm=rf.coefficients[1]==[0.]?first(domain(rf)):indmax(rf)    


    B=Evaluation(space(f),xm)
    D=Derivative(space(f))
    A=[B,D-diff(f)]
    A\[exp(f[xm]),0.]    
end

## Less accurate than solving differential equation with \
for op in (:(Base.cos),:(Base.sin),:(Base.cospi),:(Base.sinpi),:(Base.sinc))
    @eval begin
        $op{S,T}(f::Fun{S,T})=Fun(x->($op)(f[x]),space(f))
    end
end

function .^(f::Fun{ChebyshevSpace},k::Float64)
    fc = Fun(canonicalcoefficients(f))
    x=Fun(identity)

    r = sort(roots(fc))

    
    @assert length(r) <= 2
    
    if length(r) == 0
        Fun(x->f[x]^k,domain(f))
    elseif length(r) == 1
        @assert isapprox(abs(r[1]),1)
        
        if isapprox(r[1],1.)
            Fun(canonicalcoefficients((Multiplication(1-x,space(fc))\fc)).^k,JacobiWeightSpace(0.,k,domain(f)))
        else
            Fun(canonicalcoefficients((Multiplication(1+x,space(fc))\fc)).^k,JacobiWeightSpace(k,0.,domain(f)))
        end
    else
        @assert isapprox(r[1],-1)
        @assert isapprox(r[2],1) 
    
        Fun(canonicalcoefficients(linsolve(Multiplication(1-x.^2,space(fc)),fc;tolerance=eps()).^k),JacobiWeightSpace(k,k,domain(f)))  
    end
end

Base.sqrt(f::Fun{ChebyshevSpace})=f.^0.5





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


Base.cos{S<:DomainSpace{Float64},T<:Real}(f::Fun{S,T})=real(exp(im*f))
Base.sin{S<:DomainSpace{Float64},T<:Real}(f::Fun{S,T})=imag(exp(im*f))

