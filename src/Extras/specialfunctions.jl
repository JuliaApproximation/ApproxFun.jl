# division by fun 

./{T,N,a,b}(f::Fun{T,UltrasphericalSpace{a}},g::Fun{N,UltrasphericalSpace{b}})=linsolve(MultiplicationOperator(g),f;tolerance=10eps())

for op in (:./,:/)
    @eval begin
        function ($op){T}(c::Number,f::Fun{T,ChebyshevSpace})
            fc = Fun(f)
            r = roots(fc)
            x = Fun(identity)
            
            tol = 10eps()
            
            @assert length(r) <= 2
            
            if length(r) == 0
                linsolve(MultiplicationOperator(f),c;tolerance=tol)
            elseif length(r) == 1
                @assert abs(abs(r[1]) - 1.) < tol
                
                if sign(r[1]) < 0
                    Fun(coefficients(c./(fc./(x+1))),JacobiWeightSpace(-1,0,domain(f)))
                else
                    Fun(coefficients(c./(fc./(1-x))),JacobiWeightSpace(0,-1,domain(f)))                
                end 
            else
                @assert abs(r[1]+1) < tol
                @assert abs(r[2]-1) < tol                        
                
                Fun(coefficients(c./(fc./(1-x.^2))),JacobiWeightSpace(-1,-1,domain(f)))  
            end
        end
    end
end




## We use \ as the Fun constructor might miss isolated features
function Base.exp(f::Fun)
    xm=indmax(real(f))
    B=EvaluationFunctional(domain(f),xm)
    D=diff(domain(f))
    A=[B,D-diff(f)]
    A\[exp(f[xm]),0.]    
end


for op in (:(Base.cos),:(Base.sin))
    @eval begin
        ($op)(f::Fun)=Fun(x->($op)(f[x]),domain(f))
    end
end

function Base.sqrt(f::Fun)
    fc = Fun(f)
    x=Fun(identity)

    r = sort(roots(fc))
    tol= 10eps()
    
    @assert length(r) <= 2
    
    if length(r) == 0
        Fun(x->sqrt(f[x]),domain(f))
    elseif length(r) == 1
        @assert abs(abs(r[1])-1) < tol
        
        if abs(r[1]-1.) < tol
            Fun(coefficients(sqrt(MultiplicationOperator(1-x)\fc)),JacobiWeightSpace(0.,.5,domain(f)))
        else
            Fun(coefficients(sqrt(MultiplicationOperator(1+x)\fc)),JacobiWeightSpace(.5,0.,domain(f)))
        end
    else
        @assert abs(r[1]+1) < tol
        @assert abs(r[2]-1) < tol        
    
        Fun(coefficients(sqrt(linsolve(MultiplicationOperator(1-x.^2),fc;tolerance=eps()))),JacobiWeightSpace(.5,.5,domain(f)))  
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
#             B=[EvaluationFunctional(d,xm),EvaluationFunctional(d,xp)]
#             fp=diff(f)
#             fpp=diff(fp)
#             [B,fp*D^2 - fpp*D + fp.^3]\[($op)(f[xm]),($op)(f[xp])]
#         end
#     end
# end

