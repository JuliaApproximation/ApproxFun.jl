# division by fun 

./(f::IFun,g::IFun)=linsolve(MultiplicationOperator(g),f;tolerance=10eps())

for op in (:./,:/)
    @eval begin
        function ($op)(c::Number,f::IFun)
            fc = IFun(f)
            r = roots(fc)
            x = IFun(identity)
            
            tol = 10eps()
            
            @assert length(r) <= 2
            
            if length(r) == 0
                linsolve(MultiplicationOperator(f),c;tolerance=tol)
            elseif length(r) == 1
                @assert abs(abs(r[1]) - 1.) < tol
                
                if sign(r[1]) < 0
                    SingFun(IFun(c./(fc./(x+1)),f.domain),-1.,0.)
                else
                    SingFun(IFun(c./(fc./(1-x)),f.domain),0.,-1.)                
                end 
            else
                @assert abs(r[1]+1) < tol
                @assert abs(r[2]-1) < tol                        
                
                SingFun(IFun(c./(fc./(1-x.^2)),f.domain),-1.,-1.)                     
            end
        end
    end
end




## We use \ as the IFun constructor might miss isolated features
function Base.exp(f::IFun)
    xm=indmax(real(f))
    B=EvaluationOperator(f.domain,xm)
    D=diff(f.domain)
    A=[B,D-diff(f)]
    A\[exp(f[xm]),0.]    
end


for op in (:(Base.cos),:(Base.sin))
    @eval begin
        ($op)(f::IFun)=IFun(x->($op)(f[x]),f.domain)
    end
end

function Base.sqrt(f::IFun)
    fc = IFun(f)
    x=IFun(identity)

    r = sort(roots(fc))
    tol= 10eps()
    
    @assert length(r) <= 2
    
    if length(r) == 0
        IFun(x->sqrt(f[x]),f.domain)
    elseif length(r) == 1
        @assert abs(abs(r[1])-1) < tol
        
        if abs(r[1]-1.) < tol
            SingFun(IFun(sqrt(MultiplicationOperator(1-x)\fc),f.domain),0.,.5)
        else
            SingFun(IFun(sqrt(MultiplicationOperator(1+x)\fc),f.domain),.5,0.)        
        end
    else
        @assert abs(r[1]+1) < tol
        @assert abs(r[2]-1) < tol        
    
        SingFun(IFun(sqrt(linsolve(MultiplicationOperator(1-x.^2),fc;tolerance=eps())),f.domain),.5,.5)                
    end
end



## The following code will fail when xp and xm correspond to zeros of sin
# for op = (:(Base.cos),:(Base.sin))
#     @eval begin
#         function ($op)(f::IFun)
#             d=f.domain
#             D=diff(d)
#             f2=f.*f
#             xp=indmax(f2);xm=indmin(f2);
#             B=[EvaluationOperator(d,xm),EvaluationOperator(d,xp)]
#             fp=diff(f)
#             fpp=diff(fp)
#             [B,fp*D^2 - fpp*D + fp.^3]\[($op)(f[xm]),($op)(f[xp])]
#         end
#     end
# end

