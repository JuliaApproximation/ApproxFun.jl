

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
tocanonicalD(d::Line,x)=2./(π*(1+x.^2))
fromcanonical(d::Line,x)=tan(.5π*x)
fromcanonicalD(d::Line,x)=.5π*sec(.5π*x).^2



Base.length(d::Line) = Inf



==(d::Line,m::Line) = d.centre == m.centre && d.angle == m.angle

##Integration




# Transform form chebyshev U series to dirichlet-neumann U series
function uneumann_dirichlet_transform{T<:Number}(v::Vector{T})
    n=length(v)
    w=Array(T,n-4)

    for k = n-4:-1:1
        sc=(3+k)*(4+k)/((-2-k)*(1+k))
        w[k]=sc*v[k+4] 
        
        if k <= n-6
            w[k]-=sc*2*(4+k)/(5+k)*w[k+2]
        end
        if k <= n-8
            w[k]+=sc*((6+k)/(4+k))*w[k+4]
        end
    end
    
    w
end


# This takes a vector in dirichlet-neumann series on [-1,1]
# and return coefficients in T series that satisfy
# (1-x^2)u' = f
function uneumannrange_xsqd{T<:Number}(v::Vector{T})
    n = length(v)
    w=Array(T,n+1)
    
    for k=n:-1:1
        sc=-((16*(1+k)*(2+k))/(k*(3+k)*(4+k)))
        w[k+1]=sc*v[k]
        
        if k <= n-2
            w[k+1]-=sc*(k*(4+k))/(8(k+1))*w[k+3]
        end
        
        if k <= n-4
            w[k+1]+=sc*((k*(k+4))/(16(k+2)))*w[k+5]
        end
    end
    w[1]=zero(T)
    
    w
end

#integration functions
#integration is done by solving (1-x^2)^2 u' = (1-x^2)^2 M' f, u[-1] == 0

x2sec(x::Number)=abs(x)==1.? -4./π : (x.^2-1).*sec(π/2.*x)
x2sec(x::Vector)=map(x2sec,x)


integrate{T<:Number}(f::IFun{T,Line})=Fun(uneumannrange_xsqd(uneumann_dirichlet_transform(coefficients(π/2*IFun(x->x2sec(x).^2,21).*Fun(f),1))),f.domain)



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
    


