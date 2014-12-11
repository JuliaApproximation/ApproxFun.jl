



## Multiplication 

addentries!(M::Multiplication{LaurentSpace,LaurentSpace},A,k)=addentries!(LaurentOperator(M.f),A,k)

## Converison

function addentries!(C::Conversion{LaurentSpace,FourierSpace},A::ShiftArray,kr::Range)
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
function addentries!(C::Conversion{FourierSpace,LaurentSpace},A::ShiftArray,kr::Range)
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

bandinds(::Conversion{LaurentSpace,FourierSpace})=-1,1
bandinds(::Conversion{FourierSpace,LaurentSpace})=-1,1

function conversion_rule(A::LaurentSpace,B::FourierSpace)
    @assert domainscompatible(A,B)
    B
end



## Real/Imag

for TYP in (:RealOperator,:ImagOperator)
    @eval begin
        rangespace{T}(R::$TYP{ReImSpace{TaylorSpace,T}})=FourierSpace(domain(R))
    end
end

bandinds{T}(::RealOperator{ReImSpace{TaylorSpace,T}})=0,2
bandinds{T}(::ImagOperator{ReImSpace{TaylorSpace,T}})=0,1


## Re[r z^k] = r cos(k x), Re[im q z^k] = -sin(k x)
function addentries!{T}(R::RealOperator{ReImSpace{TaylorSpace,T}},A::ShiftArray,kr::Range)
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
function addentries!{T}(R::ImagOperator{ReImSpace{TaylorSpace,T}},A::ShiftArray,kr::Range)
    for k=kr
        A[k,1]+=1
    end
    A
end

# Neg

# spaces lose zeroth coefficient
for TYP in (:RealOperator,:ImagOperator)
    @eval begin
        rangespace{T}(R::$TYP{ReImSpace{HardySpace{false},T}})=DropSpace(FourierSpace(domain(R)),1)
    end
end


bandinds{T}(::RealOperator{ReImSpace{HardySpace{false},T}})=-1,1
bandinds{T}(::ImagOperator{ReImSpace{HardySpace{false},T}})=0,0


## Re[r z^(-k)] = r cos(k x), Re[im q z^(-k)] = -sin(-k x)= sin(k x)
function addentries!{T}(R::RealOperator{ReImSpace{HardySpace{false},T}},A::ShiftArray,kr::Range)
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
function addentries!{T}(R::ImagOperator{ReImSpace{HardySpace{false},T}},A::ShiftArray,kr::Range)
    for k=kr
        A[k,0]+=isodd(k)?-1:1
    end
    A
end


# spaces lose zeroth coefficient
for TYP in (:RealOperator,:ImagOperator)
    @eval begin
        rangespace{T}(R::$TYP{ReImSpace{LaurentSpace,T}})=FourierSpace(domain(R))
    end
end

bandinds{T}(::RealOperator{ReImSpace{LaurentSpace,T}})=0,2
function addentries!{T}(R::RealOperator{ReImSpace{LaurentSpace,T}},A::ShiftArray,kr::Range)
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


bandinds{S<:Union(CosSpace,SinSpace,HardySpace)}(D::Derivative{S})=0,0
rangespace{S<:CosSpace}(D::Derivative{S})=iseven(D.order)?D.space:SinSpace(domain(D))
rangespace{S<:SinSpace}(D::Derivative{S})=iseven(D.order)?D.space:CosSpace(domain(D))
rangespace{S<:HardySpace}(D::Derivative{S})=D.space



function addentries!(D::Derivative{CosSpace},A::ShiftArray,kr::Range)
    d=domain(D)
    m=D.order
    C=2π./(d.b-d.a)

    for k=kr
        A[k,0] -= (C*k)^m
    end
    
    A
end

function addentries!(D::Derivative{SinSpace},A::ShiftArray,kr::Range)
    d=domain(D)
    m=D.order
    C=2π./(d.b-d.a)

    for k=kr
        A[k,0] += (C*k)^m
    end
    
    A
end

function addentries!(D::Derivative{TaylorSpace},A::ShiftArray,kr::Range)
    d=domain(D)
    m=D.order
    C=2π./(d.b-d.a)*im

    for k=kr
        A[k,0] += (C*(k-1))^m
    end
    
    A
end

function addentries!(D::Derivative{HardySpace{false}},A::ShiftArray,kr::Range)
    d=domain(D)
    m=D.order
    C=2π./(d.b-d.a)*im

    for k=kr
        A[k,0] += (-C*k)^m
    end
    
    A
end

