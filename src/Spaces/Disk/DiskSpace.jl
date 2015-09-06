

include("JacobiSquare.jl")

export Disk

##TODO: make argument
immutable Disk <: BivariateDomain{Float64}
    radius::Float64
    center::@compat(Tuple{Float64,Float64})
end

Disk(r)=Disk(r,(0.,0.))
Disk()=Disk(1.)
Disk(::AnyDomain)=Disk(NaN,(NaN,NaN))


isambiguous(d::Disk)=isnan(d.radius) && all(isnan,d.center)
Base.convert(::Type{Disk},::AnyDomain)=Disk(AnyDomain())


#canonical is rectangle [r,0]x[-π,π]
# we assume radius and centre are zero for now
fromcanonical(D::Disk,x,t)=x*cos(t),x*sin(t)
tocanonical(D::Disk,x,y)=sqrt(x^2+y^2),atan2(y,x)
checkpoints(d::Disk)=[fromcanonical(d,(.1,.2243));fromcanonical(d,(-.212423,-.3))]

# function points(d::Disk,n,m,k)
#     ptsx=0.5*(1-gaussjacobi(n,1.,0.)[1])
#     ptst=points(PeriodicInterval(),m)
#
#     Float64[fromcanonical(d,x,t)[k] for x in ptsx, t in ptst]
# end


∂(d::Disk)=Circle(Complex(d.center...),d.radius)


immutable DiskSpace{m,a,b,JS,S} <: AbstractProductSpace{@compat(Tuple{JS,S}),Complex128,2}
    domain::Disk
    spacet::S
    DiskSpace(d,sp)=new(d,sp)
    DiskSpace(d::AnyDomain)=new(Disk(d),S())
end


DiskSpace(D::Disk,S::FunctionSpace)=DiskSpace{0,0,0,JacobiSquare,typeof(S)}(D,S)
DiskSpace(D::Disk)=DiskSpace(D,Laurent())
DiskSpace(d::AnyDomain)=DiskSpace(Disk(d))
DiskSpace()=DiskSpace(Disk())

canonicalspace(D::DiskSpace)=D

spacescompatible{m,a,b,JS,S}(A::DiskSpace{m,a,b,JS,S},B::DiskSpace{m,a,b,JS,S})=true

coefficient_type{T<:Complex}(::DiskSpace,::Type{T})=T
coefficient_type{T<:Real}(::DiskSpace,::Type{T})=Complex{T}

domain(d::DiskSpace)=d.domain
function space(D::DiskSpace,k::Integer)
    @assert k==2
    D.spacet
end

Base.getindex(D::DiskSpace,k::Integer)=space(D,k)

Space(D::Disk)=DiskSpace(D)


columnspace{M,a,b,SS}(D::DiskSpace{M,a,b,SS},k)=(m=div(k,2);JacobiSquare(M+m,a+m,b,Interval(D.domain.radius,0.)))

#transform(S::DiskSpace,V::Matrix)=transform([columnspace(S,k) for k=1:size(V,2)],S.spacet,V)


evaluate{DS<:DiskSpace}(f::Fun{DS},x,y)=ProductFun(f)[x,y]



function Base.real{JS}(f::ProductFun{JS,Laurent,DiskSpace{0,0,0,JS,Laurent}})
    cfs=f.coefficients
    n=length(cfs)

    ret=Array(Fun{JS,Float64},iseven(n)?n+1:n)
    ret[1]=real(cfs[1])

    for k=2:2:n
        # exp(1im(k-1)/2*x)=cos((k-1)/2 x) +i sin((k-1)/2 x)
        ret[k]=imag(cfs[k])
        ret[k+1]=real(cfs[k])
    end
    for k=3:2:n
        # exp(1im(k-1)/2*x)=cos((k-1)/2 x) +i sin((k-1)/2 x)
        ret[k]+=real(cfs[k])
        ret[k-1]-=imag(cfs[k])
    end

    ProductFun(ret,DiskSpace{0,0,0,JS,Fourier}(space(f).domain,Fourier()))
end
#Base.imag{S,T}(u::ProductFun{S,Larent,T})=real(TensorFun(imag(u.coefficients),space(u,2)).').'+imag(TensorFun(real(u.coefficients),space(u,2)).').'



## Conversion
# These are placeholders for future

conversion_rule{m,a,b,m2,a2,b2,JS,FS}(A::DiskSpace{m,a,b,JS,FS},
                                      B::DiskSpace{m2,a2,b2,JS,FS})=DiskSpace{max(m,m2),min(a,a2),min(b,b2),JS,FS}(A.domain,B.spacet)

function coefficients{m,a,b,m2,a2,b2,JS,FS}(cfs::Vector,
                                            A::DiskSpace{m,a,b,JS,FS},
                                          B::DiskSpace{m2,a2,b2,JS,FS})
    g=ProductFun(Fun(cfs,A))
    rcfs=Fun{typeof(columnspace(B,1)),eltype(cfs)}[Fun(g.coefficients[k],columnspace(B,k)) for k=1:length(g.coefficients)]
    Fun(ProductFun(rcfs,B)).coefficients
end


# function coefficients{S,V,SS,T}(f::ProductFun{S,V,SS,T},sp::ProductRangeSpace)
#     @assert space(f,2)==space(sp,2)

#     n=min(size(f,2),length(sp.S))
#     F=[coefficients(f.coefficients[k],rangespace(sp.S.Rdiags[k])) for k=1:n]
#     m=mapreduce(length,max,F)
#     ret=zeros(T,m,n)
#     for k=1:n
#         ret[1:length(F[k]),k]=F[k]
#     end
#     ret
# end


## Operators

isfunctional{DS<:DiskSpace}(D::Dirichlet{DS},k)=k==1

diagop(A::PlusOperator,col)=mapreduce(op->diagop(op,col),+,A.ops)
diagop(A::TimesOperator,col)=mapreduce(op->diagop(op,col),*,A.ops)
diagop(A::SpaceOperator,col)=diagop(A.op,col)
diagop(A::ConstantOperator,col)=ConstantOperator(A.c)
diagop(A::ConstantTimesOperator,col)=A.c*diagop(A.op,col)

diagop{DS<:DiskSpace}(D::Dirichlet{DS},col)=Evaluation(columnspace(domainspace(D),col),false,D.order)

function diagop{DS<:DiskSpace}(L::Laplacian{DS},col)
    csp=columnspace(domainspace(L),col)
    rsp=columnspace(rangespace(L),col)
    Dt=Derivative(space(domainspace(L),2))
    c=Dt[col,col]


    r=Fun(identity,[domain(L).radius,0.])
    D=Derivative(csp)
    Δ=D^2+(1/r)*D+Multiplication((c/r)^2,csp)

    Δ^L.order
end



isproductop{DS1<:DiskSpace,DS2<:DiskSpace}(C::Conversion{DS1,DS2})=true
diagop{DS1<:DiskSpace,DS2<:DiskSpace}(C::Conversion{DS1,DS2},col)=Conversion(columnspace(domainspace(C),col),
                                                                               columnspace(rangespace(C),col))
#deprod{DS1<:DiskSpace,DS2<:DiskSpace}(C::Conversion{DS1,DS2},k,::Colon)=ConstantOperator(1.0)


lap(d::Disk)=Laplacian(Space(d))
dirichlet(d::Disk)=Dirichlet(Space(d))
neumann(d::Disk)=Neumann(Space(d))

lap(d::DiskSpace)=Laplacian(d)
dirichlet(d::DiskSpace)=Dirichlet(d)
neumann(d::DiskSpace)=Neumann(d)



function rangespace{m,a,b,JS,S}(L::Laplacian{DiskSpace{m,0,0,JacobiSquare,Laurent}})
    sp=domainspace(L)
    DiskSpace{m-2L.order,a+2L.order,b+2L.order,JS,S}(sp.domain,sp.spacet)
end



# special case of integer modes
function diagop{b}(L::Laplacian{DiskSpace{0,0,b,JacobiSquare,Laurent}},col)
    S=columnspace(domainspace(L),col)
    Dp=DDp(S)
    Dm=DDm(rangespace(Dp))
    2Dm*Dp
end



function rangespace{b,JS,S}(L::Laplacian{DiskSpace{0,0,b,JS,S}})
    sp=domainspace(L)
    DiskSpace{0,0,b+2,JS,S}(sp.domain,sp.spacet)
end
