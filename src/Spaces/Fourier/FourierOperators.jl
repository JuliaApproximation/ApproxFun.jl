
## Evaluation

getindex(T::Evaluation{Taylor},cols::Range)=mappoint(domain(T),Circle(),T.x).^(cols-1)    


## Multiplication 

addentries!(M::Multiplication{Laurent,Laurent},A,k)=addentries!(LaurentOperator(M.f),A,k)
bandinds(M::Multiplication{Laurent,Laurent})=bandinds(LaurentOperator(M.f))

## Converison

function addentries!(C::Conversion{Laurent,Fourier},A,kr::Range)
    for k=kr
        if k==1
            A[k,0]+=1.
        elseif iseven(k)
            A[k,0]+=-1.im
            A[k,1]+=1.im
        else #isodd(k)
            A[k,0]+=1
            A[k,-1]+=1
        end
    end
    A
end
function addentries!(C::Conversion{Fourier,Laurent},A,kr::Range)
    for k=kr
        if k==1
            A[k,0]+=1.
        elseif iseven(k)
            A[k,0]+=0.5im
            A[k,1]+=0.5
        else #isodd(k)
            A[k,0]+=0.5
            A[k,-1]+=-0.5im
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



## Real/Imag

for TYP in (:RealOperator,:ImagOperator)
    @eval begin
        rangespace{T}(R::$TYP{ReImSpace{Taylor,T}})=Fourier(domain(R))
    end
end

bandinds{T}(::RealOperator{ReImSpace{Taylor,T}})=0,2
bandinds{T}(::ImagOperator{ReImSpace{Taylor,T}})=0,1


## Re[r z^k] = r cos(k x), Re[im q z^k] = -sin(k x)
function addentries!{T}(R::RealOperator{ReImSpace{Taylor,T}},A,kr::Range)
    for k=kr
        if isodd(k)         # real part
            A[k,0]+=1        
        elseif iseven(k)    # imag part
            A[k,2]+=-1
        end
    end
    A
end

## Im[r z^k] = r sin(k x), Im[im q z^k] = cos(k x)
function addentries!{T}(R::ImagOperator{ReImSpace{Taylor,T}},A,kr::Range)
    for k=kr
        A[k,1]+=1
    end
    A
end

# Neg

# spaces lose zeroth coefficient
for TYP in (:RealOperator,:ImagOperator)
    @eval begin
        rangespace{T}(R::$TYP{ReImSpace{Hardy{false},T}})=DropSpace(Fourier(domain(R)),1)
    end
end


bandinds{T}(::RealOperator{ReImSpace{Hardy{false},T}})=-1,1
bandinds{T}(::ImagOperator{ReImSpace{Hardy{false},T}})=0,0


## Re[r z^(-k)] = r cos(k x), Re[im q z^(-k)] = -sin(-k x)= sin(k x)
function addentries!{T}(R::RealOperator{ReImSpace{Hardy{false},T}},A,kr::Range)
    for k=kr
        if isodd(k)    # imag part
            A[k,1]+=1            
        elseif iseven(k)         # real part
            A[k,-1]+=1        
        end
    end
    A
end

## Im[r z^(-k)] = r sin(-k x)=-r sin(kx), Im[im q z^(-k)] = cos(-k x)=cos(kx)
function addentries!{T}(R::ImagOperator{ReImSpace{Hardy{false},T}},A,kr::Range)
    for k=kr
        A[k,0]+=isodd(k)?-1:1
    end
    A
end


# spaces lose zeroth coefficient
for TYP in (:RealOperator,:ImagOperator)
    @eval begin
        rangespace{T}(R::$TYP{ReImSpace{Laurent,T}})=Fourier(domain(R))
    end
end

bandinds{T}(::RealOperator{ReImSpace{Laurent,T}})=0,2
function addentries!{T}(R::RealOperator{ReImSpace{Laurent,T}},A,kr::Range)
    for k=kr
        if isodd(k)    # real part
            A[k,0]+=1
        elseif iseven(k)         # odd part
            A[k,2]+=iseven(div(k,2))?-1:1
        end
    end
    A
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
            A[k,0] -= (C*(k-1))^m
        else
            A[k,1] -= (C*k)^m
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
            A[k,0] += (C*k)^m
        elseif mod(m,4)==2
            A[k,0] += -(C*k)^m        
        elseif k>1 && mod(m,4)==1
            A[k,-1] += (C*(k-1))^m
        elseif k>1 && mod(m,4)==3
            A[k,-1] += -(C*(k-1))^m            
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
            A[k,0] += (C*k)^(-m)
        elseif mod(m,4)==2
            A[k,0] += -(C*k)^(-m)
        elseif k>1 && mod(m,4)==1
            A[k,-1] += -(C*(k-1))^(-m)
        elseif k>1 && mod(m,4)==3
            A[k,-1] += (C*(k-1))^(-m)
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
            A[k,0] += (C*k)^(-m)
        elseif mod(m,4)==2
            A[k,0] += -(C*k)^(-m)
        elseif mod(m,4)==1
            A[k,0] += (C*k)^(-m)
        elseif mod(m,4)==3
            A[k,0] += -(C*k)^(-m)
        end
    end
    
    A
end



## Hardy space

function bandinds{s}(D::Derivative{Hardy{s}})
    d=domain(D)
    if isa(d,PeriodicInterval)
        (0,0)
    elseif isa(d,Circle)
        s?(0,D.order):(-D.order,0)
    else
        error("Derivative not defined for "*string(typeof(d)))
    end
end
rangespace{S<:Hardy}(D::Derivative{S})=D.space

function taylor_derivative_addentries!(d::PeriodicInterval,m::Integer,A,kr::Range)
    C=2π./(d.b-d.a)*im
    for k=kr
        A[k,0] += (C*(k-1))^m
    end
    A
end

function hardyfalse_derivative_addentries!(d::PeriodicInterval,m::Integer,A,kr::Range)
    C=2π./(d.b-d.a)*im
    for k=kr
        A[k,0] += (-C*k)^m
    end
    A
end



function taylor_derivative_addentries!(d::Circle,m::Integer,A,kr::Range)
    C=d.radius^(-m)

    for k=kr
        D=k
        for j=k+1:k+m-1
          D*=j  
        end
        A[k,m] += C*D
    end

    A
end

function hardyfalse_derivative_addentries!(d::Circle,m::Integer,A,kr::Range)
    C=(-d.radius)^(-m)

    for k=max(m+1,kr[1]):kr[end]
        D=k-m
        for j=k-m+1:k-1
          D*=j  
        end
        A[k,-m] += C*D
    end

    A
end


addentries!(D::Derivative{Taylor},A,kr::Range)=taylor_derivative_addentries!(domain(D),D.order,A,kr)
addentries!(D::Derivative{Hardy{false}},A,kr::Range)=hardyfalse_derivative_addentries!(domain(D),D.order,A,kr)




## Integral



function Integral(S::Taylor,m)
    @assert isa(domain(S),Circle)
    Integral{Taylor,Complex{Float64}}(S,m)
end

function Integral(S::Hardy{false},m)
    @assert isa(domain(S),PeriodicInterval)
    Integral{Hardy{false},Complex{Float64}}(S,m)
end

function bandinds(D::Integral{Taylor})
    d=domain(D)
    @assert isa(d,Circle)
    (-D.order,0)
end
rangespace(D::Integral{Taylor})=D.space
rangespace(Q::Integral{Hardy{false}})=Q.space

function addentries!(D::Integral{Taylor},A,kr::Range)
    d=domain(D)
    m=D.order
    @assert isa(d,Circle)
    
    C=d.radius^m

    for k=max(m+1,kr[1]):kr[end]
        D=k-m
        for j=k-m+1:k-1
          D*=j  
        end
        A[k,-m] += C/D
    end

    A    
end


function bandinds{n,T,DD}(D::Integral{DropSpace{Hardy{false},n,T,DD}})
    d=domain(D)
    @assert isa(d,Circle)
    @assert D.order==n
    (0,0)
end
rangespace{n,T,DD}(D::Integral{DropSpace{Hardy{false},n,T,DD}})=D.space.space

function addentries!{n,T,DD}(D::Integral{DropSpace{Hardy{false},n,T,DD}},A,kr::Range)
    d=domain(D)
    m=D.order
    @assert isa(d,Circle)
    
    C=(-d.radius)^m

    for k=kr
        D=k
        for j=k+1:k+m-1
          D*=j  
        end
        A[k,0] += C/D
    end

    A    
end



function bandinds(D::Integral{Hardy{false}})
    d=domain(D)
    @assert isa(d,PeriodicInterval)
    (0,0)
end
rangespace(D::Integral{Taylor})=D.space


function addentries!(D::Integral{Hardy{false}},A,kr::Range)
    d=domain(D)
    m=D.order
    @assert isa(d,PeriodicInterval)
    
    C=2π./(d.b-d.a)*im
    for k=kr
        A[k,0] += (-C*k)^(-m)
    end
    A
end



function bandinds{n,T,DD}(D::Integral{DropSpace{Taylor,n,T,DD}})
    d=domain(D)
    @assert isa(d,PeriodicInterval)
    (0,0)
end
rangespace{n,T,DD}(D::Integral{DropSpace{Taylor,n,T,DD}})=D.space

function addentries!{n,T,DD}(D::Integral{DropSpace{Taylor,n,T,DD}},A,kr::Range)
    d=domain(D)
    m=D.order
    @assert isa(d,PeriodicInterval)
    
    C=2π./(d.b-d.a)*im
    for k=kr
        A[k,0] += (C*(k+n-1))^(-m)
    end
    A   
end
