using ApproxFun
    import ApproxFun:canonicalspace,spacescompatible,points,RealBasis,transform


immutable PiecewiseInterval{T<:Number} <:Domain{T}
    points::Vector{T}
end

PiecewiseInterval(d::Number...)=PiecewiseInterval([d...])

Base.length(d::PiecewiseInterval)=length(d.points)-1
Base.getindex(d::PiecewiseInterval,j::Integer)=Interval(d.points[j],d.points[j+1])

function points(d::PiecewiseInterval,n)
   k=div(n,length(d))
    r=n-length(d)*k

    [vcat([points(d[j],k+1) for j=1:r]...);
        vcat([points(d[j],k) for j=r+1:length(d)]...)]
end


Base.rand(d::PiecewiseInterval)=rand(d[rand(1:length(d))])
checkpoints{T}(d::PiecewiseInterval{T})=union(T[checkpoints(d[k]) for k=1:length(d)])



immutable ContinuousSpace <: FunctionSpace{RealBasis,UnionDomain}
    domain::PiecewiseInterval
end



function transform(S::ContinuousSpace,vals::Vector)
    n=length(vals)
    d=domain(S)
    K=length(d)
   k=div(n,K)

    PT=promote_type(eltype(d),eltype(vals))
    if k==0
        ret=Array(PT,n)
        for j=1:n
            ret[j]=transform(ChebyshevDirichlet{1,1}(d[j]),[vals[j]])[1]
        end

        ret
    else
        ret=zeros(PT,K+1+K*(k-1))





        r=n-K*k
        M=Array(PT,k+1,K)

        for j=1:r
            M[:,j]=transform(ChebyshevDirichlet{1,1}(d[j]),vals[(j-1)*(k+1)+1:j*(k+1)])
        end
        for j=r+1:length(d)
            M[1:k,j]=transform(ChebyshevDirichlet{1,1}(d[j]),vals[r*(k+1)+(j-r-1)*k+1:r*(k+1)+(j-r)*k])
            M[k+1,j]=zero(PT)
        end

        vec(M.')
    end
end


d=PiecewiseInterval(1.,2.,3.)
S=ContinuousSpace(d)


exp(1.),exp(2.)
cfs[1]-cfs[2],cfs[1]+cfs[2]

cfs=Fun(f,ChebyshevDirichlet{1,1}([1.,2.])).coefficients

f=x->exp(x)
vals=map(f,points(d,100))
n=length(vals)
d=domain(S)
K=length(d)
k=div(n,K)


PT=Float64
ret=zeros(PT,K+1+K*(k-1))
#.'
j=1
cfs=transform(ChebyshevDirichlet{1,1}(d[j]),vals[(j-1)*(k+1)+1:j*(k+1)])
ret[1]=cfs[1]-cfs[2]
ret[2]=cfs[1]+cfs[2]
ret[3:k+1]=cfs[3:end]

        for j=1:r
M[:,j]=transform(ChebyshevDirichlet{1,1}(d[j]),vals[(j-1)*(k+1)+1:j*(k+1)])
end

r=n-K*k
ret
