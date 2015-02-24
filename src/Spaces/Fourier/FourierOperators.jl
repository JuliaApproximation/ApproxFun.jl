

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

    for k=kr
        if mod(m,4)==0
            A[k,k] += (C*(k-1))^m
        elseif mod(m,4)==2
            A[k,k] -= (C*(k-1))^m
        elseif mod(m,4)==1
            A[k,k+1] -= (C*k)^m
        elseif mod(m,4)==3
            A[k,k+1] += (C*k)^m            
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

Integral(::CosSpace,m)=error("Integral not defined for CosSpace.  Use Integral(SliceSpace(CosSpace(),1)) if first coefficient vanishes.")


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


function bandinds{T,DD}(D::Integral{SliceSpace{1,1,CosSpace,T,DD}})
    d=domain(D)
    @assert isa(d,PeriodicInterval)
    (0,0)
end
rangespace{T,DD}(D::Integral{SliceSpace{1,1,CosSpace,T,DD}})=iseven(D.order)?D.space:SinSpace(domain(D))

function addentries!{T,DD}(D::Integral{SliceSpace{1,1,CosSpace,T,DD}},A,kr::Range)
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

# CosSpace Multiplicaiton is the same as Chebyshev

bandinds{Sp<:CosSpace}(M::Multiplication{Sp,Sp})=(1-length(M.f.coefficients),length(M.f.coefficients)-1)
rangespace{Sp<:CosSpace}(M::Multiplication{Sp,Sp})=domainspace(M)
addentries!{Sp<:CosSpace}(M::Multiplication{Sp,Sp},A,kr::UnitRange)=chebmult_addentries!(M.f.coefficients,A,kr)


function addentries!{Sp<:SinSpace}(M::Multiplication{Sp,Sp},A,kr::UnitRange)
    a=M.f.coefficients
    toeplitz_addentries!(0.5ShiftVector([-flipud(a);0.],a),A,kr)
    hankel_addentries!(0.5a,A,max(kr[1],2):kr[end])    
    A
end

bandinds{Sp<:SinSpace}(M::Multiplication{Sp,Sp})=-length(M.f)-1,length(M.f)-1
rangespace{Sp<:SinSpace}(M::Multiplication{Sp,Sp})=CosSpace(domain(M))


function addentries!{Sp<:SinSpace,Cs<:CosSpace}(M::Multiplication{Sp,Cs},A,kr::Range)
    a=M.f.coefficients
    toeplitz_addentries!(0.5ShiftVector(flipud(a[2:end]),[a[1];0.;-a]),A,kr)
    hankel_addentries!(0.5a,A,kr)
    A
end

bandinds{Sp<:SinSpace,Cs<:CosSpace}(M::Multiplication{Sp,Cs})=1-length(M.f),length(M.f)+1
rangespace{Sp<:SinSpace,Cs<:CosSpace}(M::Multiplication{Sp,Cs})=SinSpace(domain(M))



function addentries!{Sp<:SinSpace,Cs<:CosSpace}(M::Multiplication{Cs,Sp},A,kr::Range)
    a=M.f.coefficients
    toeplitz_addentries!(0.5a,A,kr)
    if length(a)>=3
        hankel_addentries!(-0.5a[3:end],A,kr)
    end
    A
end

bandinds{Sp<:SinSpace,Cs<:CosSpace}(M::Multiplication{Cs,Sp})=(1-length(M.f.coefficients),length(M.f.coefficients)-1)
rangespace{Sp<:SinSpace,Cs<:CosSpace}(M::Multiplication{Cs,Sp})=SinSpace(domain(M))



function Multiplication{T}(a::Fun{Fourier,T},sp::Fourier)
    d=domain(a)
    c,s=vec(a,1),vec(a,2)
    O=BandedOperator{T}[Multiplication(c,CosSpace(d)) Multiplication(s,SinSpace(d));
                        Multiplication(s,CosSpace(d)) Multiplication(c,SinSpace(d))]
    M=SpaceOperator(InterlaceOperator(O),space(a),sp)
end


## Definite integral

Σ(sp::Fourier)=isa(domain(sp),PeriodicInterval)?Σ{Fourier,Float64}(sp):Σ{Fourier,Complex{Float64}}(sp)

function getindex{T}(S::Σ{Fourier,T},kr::Range)
    d = domain(S)
    if isa(d,PeriodicInterval)
        T[k == 1?  d.b-d.a : zero(T) for k=kr]
    else
        @assert isa(d,Circle)
        T[k == 2?  -π*d.radius : (k==3?π*im*d.radius :zero(T)) for k=kr]        
    end
end

datalength(S::Σ{Fourier})=isa(domain(S),PeriodicInterval)?1:3

