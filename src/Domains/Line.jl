

export Line, PeriodicLine



## Standard interval

type Line <: IntervalDomain
    centre
    angle
    
    Line(c,a)=(@assert c==a==0.; new(c,a))
end

Line()=Line(0.,0.)

## Map interval



tocanonical(d::Line,x)=2/π*atan(x)
tocanonicalD(d::Line,x)=2./(π(1+x.^2))
fromcanonical(d::Line,x)=tan(.5π*x)
fromcanonicalD(d::Line,x)=.5π*sec(.5π*x).^2



Base.length(d::Interval) = Inf



==(d::Line,m::Line) = d.centre == m.centre && d.angle == m.angle



##Integration


#integration functions
#integration is done by solving (1-x^2)^2 u' = (1-x^2)^2 M' f, u[-1] == 0

x2sec(x::Number)=abs(x)==1.? -4./π : (x.^2-1).*sec(π/2.*x)
x2sec(x::Vector)=map(x2sec,x)

##TODO: reimplement

# const x2sec2fun=IFun(x->x2sec(x).^2)
# 
# linecumsumop=[EvaluationOperator(-1.),DifferentialOperator([0.,2./π.*IFun(x->(x.^2-1).^2)])];

function integrate{T<:Number}(f::IFun{T,Line})
    #TODO: choose length smarter

    g=chop(ApproxFun.x2sec2fun.*IFun(f.coefficients),1000eps())
    
    IFun(tomatrix(linecumsumop,length(g)+1)\[0.,coefficients(g,1)],f.domain)
end


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
    


