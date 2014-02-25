# division by fun 

for op = (:./,:/)
    @eval begin
        ($op)(c::Union(Number,IFun),f::IFun)=MultiplicationOperator(f)\[c]
    end
end


## We use \ as the IFun constructor might miss isolated features
function Base.exp(f::IFun)
    xm=indmax(f)
    B=EvaluationOperator(f.domain,xm)
    D=diff(f.domain)
    A=[B,D-diff(f)]
    A\[exp(f[xm]),0.]    
end


for op = (:(Base.cos),:(Base.sin),:(Base.sqrt)
    @eval begin
        ($op)(f::IFun)=IFun(x->($op)(f[x]),f.domain)
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

