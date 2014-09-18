# division by fun 

./{a,b}(f::Fun{UltrasphericalSpace{a}},g::Fun{UltrasphericalSpace{b}})=linsolve(Multiplication(g),f;tolerance=10eps())

for op in (:./,:/)
    @eval begin
        function ($op)(c::Number,f::Fun{ChebyshevSpace})
            fc = Fun(canonicalcoefficients(f),Interval())
            r = roots(fc)
            x = Fun(identity)
            
            tol = 10eps()
            
            @assert length(r) <= 2
            
            if length(r) == 0
                linsolve(Multiplication(f),c;tolerance=tol)
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


for op in (:(Base.cos),:(Base.sin))
    @eval begin
        ($op)(f::Fun)=Fun(x->($op)(f[x]),domain(f))
    end
end

function Base.sqrt(f::Fun{ChebyshevSpace})
    fc = Fun(canonicalcoefficients(f))
    x=Fun(identity)

    r = sort(roots(fc))
    tol= 10eps()
    
    @assert length(r) <= 2
    
    if length(r) == 0
        Fun(x->sqrt(f[x]),domain(f))
    elseif length(r) == 1
        @assert abs(abs(r[1])-1) < tol
        
        if abs(r[1]-1.) < tol
            Fun(canonicalcoefficients(sqrt(Multiplication(1-x)\fc)),JacobiWeightSpace(0.,.5,domain(f)))
        else
            Fun(canonicalcoefficients(sqrt(Multiplication(1+x)\fc)),JacobiWeightSpace(.5,0.,domain(f)))
        end
    else
        @assert abs(r[1]+1) < tol
        @assert abs(r[2]-1) < tol        
    
        Fun(canonicalcoefficients(sqrt(linsolve(Multiplication(1-x.^2),fc;tolerance=eps()))),JacobiWeightSpace(.5,.5,domain(f)))  
    end
end



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

