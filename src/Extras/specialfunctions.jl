# division by fun 

./{a,b}(f::Fun{UltrasphericalSpace{a}},g::Fun{UltrasphericalSpace{b}})=linsolve(Multiplication(g,space(f)),f;tolerance=10eps())

for op in (:./,:/)
    @eval begin
        function ($op)(c::Number,f::Fun{ChebyshevSpace})
            fc = Fun(canonicalcoefficients(f),Interval())
            r = roots(fc)
            x = Fun(identity)
            
            tol = 10eps()
            
            @assert length(r) <= 2
            
            if length(r) == 0
                linsolve(Multiplication(f,space(f)),c;tolerance=tol)
            elseif length(r) == 1
                @assert abs(abs(r[1]) - 1.) < tol
                
                if sign(r[1]) < 0
                    Fun(canonicalcoefficients(c./(fc./(x+1))),JacobiWeightSpace(-1,0,domain(f)))
                else
                    Fun(canonicalcoefficients(c./(fc./(1-x))),JacobiWeightSpace(0,-1,domain(f)))                
                end 
            else
                @assert abs(r[1]+1) < tol
                @assert abs(r[2]-1) < tol                        
                
                Fun(canonicalcoefficients(c./(fc./(1-x.^2))),JacobiWeightSpace(-1,-1,domain(f)))  
            end
        end
    end
end




## We use \ as the Fun constructor might miss isolated features
function Base.exp(f::Fun)
    xm=indmax(real(f))
    B=Evaluation(domain(f),xm)
    D=diff(domain(f))
    A=[B,D-diff(f)]
    A\[exp(f[xm]),0.]    
end

## Less accurate than solving differential equation with \
for op in (:(Base.cos),:(Base.sin),:(Base.cospi),:(Base.sinpi),:(Base.sinc))
    @eval begin
        ($op)(f::Fun)=Fun(x->($op)(f[x]),domain(f))
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
#=
function Base.cos(f::Fun)
    xm=indmax(imag(f))
    B=Evaluation(domain(f),xm)
    D=diff(domain(f))
    A=[B,D-im*diff(f)]
    real(A\[exp(im*f[xm]),0.])
end

function Base.sin(f::Fun)
    xm=indmax(imag(f))
    B=Evaluation(domain(f),xm)
    D=diff(domain(f))
    A=[B,D-im*diff(f)]
    imag(A\[exp(im*f[xm]),0.])
end
=#