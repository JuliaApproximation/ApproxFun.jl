

export Line, PeriodicLine



## Standard interval


type Line <: IntervalDomain
    centre::Float64  ##TODO Allow complex
    angle::Float64
    
    Line(c,a)=(@assert c==a==0.; new(c,a))
end

Line()=Line(0.,0.)

## Map interval



tocanonical(d::Line,x)=2/π*atan(x)
tocanonicalD(d::Line,x)=2./(π(1+x.^2))
fromcanonical(d::Line,x)=tan(.5π*x)
fromcanonicalD(d::Line,x)=.5π*sec(.5π*x).^2



Base.length(d::Line) = Inf



==(d::Line,m::Line) = d.centre == m.centre && d.angle == m.angle



##Integration


#integration functions
#integration is done by solving (1-x^2)^2 u' = (1-x^2)^2 M' f, u[-1] == 0

x2sec(x::Number)=abs(x)==1.? -4./π : (x.^2-1).*sec(π/2.*x)
x2sec(x::Vector)=map(x2sec,x)

##TODO: reimplement

const x2sec2fun_glob =IFun(x->x2sec(x).^2,21)
 
const linecumsumop_glob=[EvaluationOperator(Interval(),-1.),IFun(x->2./π.*(x.^2-1).^2,5)*diff(Interval())]


integrate{T<:Number}(f::IFun{T,Line})=integrate(f,linecumsumop_glob,x2sec2fun_glob)
function integrate{T<:Number,O<:Operator}(f::IFun{T,Line},linecumsumop::Vector{O},x2sec2fun::IFun)
    g=chop(x2sec2fun.*IFun(f.coefficients),100eps())
    
    IFun(linecumsumop\[0.,g],f.domain)
end


##multiplybyx

function identity_fun(d::Line)
    ct=Fun(x->x.*cot(π*x/2),28)
    x=Fun(identity)
    u=SingFun(ct./(1-x.^2),1.,1.)
    SingFun(IFun((x./u).fun,Line()),-1.,-1.)
end


# function multiplybyx{T<:Number,D<:Line}(f::IFun{T,D})
#     ct=Fun(x->x.*cot(π*x/2),28)
#     x=Fun(identity)
#     u=SingFun(ct./(1-x.^2),1.,1.)
#     IFun((x.*IFun(f)./u).fun./(1-x.^2),f.domain)
# end


## sample


function sample(f::Fun2D{IFun{Float64,Line}},n::Integer)
    cf=normalizedcumsum(sum(f,1))
    CB=coefficientmatrix(map(cumsum,f.B))
    
    ry=samplecdf(cf,n)
    fA=evaluate(f.A,ry)
    CBfA=CB*fA  #cumsums at points
    multiply_oneatright!(CBfA)
    
    rx=fromcanonical(first(f.B).domain,bisectioninv(CBfA,rand(n)))
    
    [rx ry]
end
    


