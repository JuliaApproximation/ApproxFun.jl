

## Converison

function addentries!(C::Conversion{Laurent,Fourier},A,kr::Range)
    for k=kr
        if k==1
            A[k,k]+=1.
        elseif iseven(k)
            A[k,k]+=-1.im
            A[k,k+1]+=1.im
        else #isodd(k)
            A[k,k]+=1
            A[k,k-1]+=1
        end
    end
    A
end
function addentries!(C::Conversion{Fourier,Laurent},A,kr::Range)
    for k=kr
        if k==1
            A[k,k]+=1.
        elseif iseven(k)
            A[k,k]+=0.5im
            A[k,k+1]+=0.5
        else #isodd(k)
            A[k,k]+=0.5
            A[k,k-1]+=-0.5im
        end
    end
    A
end

bandinds(::Conversion{Laurent,Fourier})=-1,1
bandinds(::Conversion{Fourier,Laurent})=-1,1

function conversion_rule(A::Laurent,B::Fourier)
    @assert domainscompatible(A,B)
    B
end









### Cos/Sine


bandinds(D::Derivative{CosSpace})=iseven(D.order)?(0,0):(0,1)
bandinds(D::Derivative{SinSpace})=iseven(D.order)?(0,0):(-1,0)
rangespace{S<:CosSpace}(D::Derivative{S})=iseven(D.order)?D.space:SinSpace(domain(D))
rangespace{S<:SinSpace}(D::Derivative{S})=iseven(D.order)?D.space:CosSpace(domain(D))


function addentries!(D::Derivative{CosSpace},A,kr::Range)
    d=domain(D)
    @assert isa(d,PeriodicInterval)
    m=D.order
    C=2π./(d.b-d.a)
    @assert m <= 2 # the def is missing the sign change

    for k=kr
        if iseven(m)
            A[k,k] -= (C*(k-1))^m
        else
            A[k,k+1] -= (C*k)^m
        end
    end
    
    A
end

function addentries!(D::Derivative{SinSpace},A,kr::Range)
    d=domain(D)
    @assert isa(d,PeriodicInterval)    
    m=D.order
    C=2π./(d.b-d.a)

    for k=kr
        if mod(m,4)==0
            A[k,k] += (C*k)^m
        elseif mod(m,4)==2
            A[k,k] += -(C*k)^m        
        elseif k>1 && mod(m,4)==1
            A[k,k-1] += (C*(k-1))^m
        elseif k>1 && mod(m,4)==3
            A[k,k-1] += -(C*(k-1))^m            
        end
    end
    
    A
end

Integral(::CosSpace,m)=error("Integral not defined for CosSpace.  Use Integral(DropSpace(CosSpace(),1)) if first coefficient vanishes.")


bandinds(D::Integral{SinSpace})=iseven(D.order)?(0,0):(-1,0)
rangespace{S<:CosSpace}(D::Integral{S})=iseven(D.order)?D.space:SinSpace(domain(D))
rangespace{S<:SinSpace}(D::Integral{S})=iseven(D.order)?D.space:CosSpace(domain(D))

function addentries!(D::Integral{SinSpace},A,kr::Range)
    d=domain(D)
    @assert isa(d,PeriodicInterval)    
    m=D.order
    C=2π./(d.b-d.a)

    for k=kr
        if mod(m,4)==0
            A[k,k] += (C*k)^(-m)
        elseif mod(m,4)==2
            A[k,k] += -(C*k)^(-m)
        elseif k>1 && mod(m,4)==1
            A[k,k-1] += -(C*(k-1))^(-m)
        elseif k>1 && mod(m,4)==3
            A[k,k-1] += (C*(k-1))^(-m)
        end
    end
    
    A
end


function bandinds{T,DD}(D::Integral{DropSpace{CosSpace,1,T,DD}})
    d=domain(D)
    @assert isa(d,PeriodicInterval)
    (0,0)
end
rangespace{T,DD}(D::Integral{DropSpace{CosSpace,1,T,DD}})=iseven(D.order)?D.space:SinSpace(domain(D))

function addentries!{T,DD}(D::Integral{DropSpace{CosSpace,1,T,DD}},A,kr::Range)
    d=domain(D)
    @assert isa(d,PeriodicInterval)    
    m=D.order
    C=2π./(d.b-d.a)

    for k=kr
        if mod(m,4)==0
            A[k,k] += (C*k)^(-m)
        elseif mod(m,4)==2
            A[k,k] += -(C*k)^(-m)
        elseif mod(m,4)==1
            A[k,k] += (C*k)^(-m)
        elseif mod(m,4)==3
            A[k,k] += -(C*k)^(-m)
        end
    end
    
    A
end



function addentries!{Sp<:SinSpace}(M::Multiplication{Sp,Sp},A,kr::UnitRange)
    a=M.f
    toeplitz_addentries!(0.5ShiftVector([-flipud(a.coefficients),0.],a.coefficients),A,kr)
    hankel_addentries!(0.5*a.coefficients,A,max(kr[1],2):kr[end])    
    A
end

bandinds{Sp<:SinSpace}(M::Multiplication{Sp,Sp})=-length(M.f)-1,length(M.f)-1
rangespace{Sp<:SinSpace}(M::Multiplication{Sp,Sp})=CosSpace(domain(M))



