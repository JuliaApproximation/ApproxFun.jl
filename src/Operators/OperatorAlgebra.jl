

export PlusOperator,TimesOperator




type PlusOperator{T<:Number,B<:Operator} <: BandedOperator{T} 
    ops::Vector{B}
end


typeofdata{T<:Number}(::Operator{T})=T
function PlusOperator{B<:Operator}(ops::Vector{B})
    T = Float64
    
    for op in ops
        if typeofdata(op) == Complex{Float64}
            T = Complex{Float64}
        end
    end
    
    pops=promotespaces(ops)
    PlusOperator{T,Operator}(pops)
end


function domainspace(P::PlusOperator)
    for op in P.ops
        sp = domainspace(op)
        
        if sp != AnySpace()
            return sp
        end
    end
    
    AnySpace()
end

function rangespace(P::PlusOperator)
    for op in P.ops
        sp = rangespace(op)
        
        if sp != AnySpace()
            return sp
        end
    end
    
    AnySpace()
end

domain(P::PlusOperator)=commondomain(P.ops)

function bandinds(P::PlusOperator)
    b1,b2=0,0
    for op in P.ops
        br=bandinds(op)
        b1=min(br[1]::Int,b1)
        b2=max(br[end]::Int,b2)
    end
    (b1,b2)
end


function addentries!(P::PlusOperator,A::ShiftArray,kr::Range1)
    for op in P.ops
        addentries!(op,A,kr)
    end
    
    A
end

+(A::PlusOperator,B::PlusOperator)=PlusOperator([A.ops,B.ops])
+(A::PlusOperator,B::Operator)=PlusOperator([A.ops,B])
+(A::Operator,B::PlusOperator)=PlusOperator([A,B.ops])
+(A::Operator,B::Operator)=PlusOperator([A,B])
+(A::Operator,f::AbstractFun)=A+MultiplicationOperator(f)
+(f::AbstractFun,A::Operator)=MultiplicationOperator(f)+A
-(A::Operator,f::AbstractFun)=A+MultiplicationOperator(-f)
-(f::AbstractFun,A::Operator)=MultiplicationOperator(f)-A



+(c::UniformScaling,A::Operator)=ConstantOperator(1.0c.λ)+A
+(A::Operator,c::UniformScaling)=A+ConstantOperator(1.0c.λ)
-(c::UniformScaling,A::Operator)=ConstantOperator(1.0c.λ)-A
-(A::Operator,c::UniformScaling)=A-ConstantOperator(1.0c.λ)


+(c::Number,A::Operator)=ConstantOperator(c)+A
+(A::Operator,c::Number)=A+ConstantOperator(c)
-(c::Number,A::Operator)=ConstantOperator(c)-A
-(A::Operator,c::Number)=A-ConstantOperator(c)

.+(c::Number,A::Operator)=ConstantOperator(c)+A
.+(A::Operator,c::Number)=A+ConstantOperator(c)
.-(c::Number,A::Operator)=ConstantOperator(c)-A
.-(A::Operator,c::Number)=A-ConstantOperator(c)



## Times Operator


type TimesOperator{T<:Number,B<:Operator} <: BandedOperator{T}
    ops::Vector{B}
    
    function TimesOperator{B}(ops::Vector{B})
        for k=1:length(ops)-1
            @assert domainspace(ops[k])==AnySpace() || rangespace(ops[k+1])==AnySpace() || domainspace(ops[k])==rangespace(ops[k+1])
        end
        
        new(ops)
    end
end


function TimesOperator{B<:Operator}(ops::Vector{B})
    T = Float64
    
    for op in ops
        if typeofdata(op) == Complex{Float64}
            T = Complex{Float64}
        end
    end
    
    TimesOperator{T,B}(ops)
end


TimesOperator(A::TimesOperator,B::TimesOperator)=TimesOperator([A.ops,B.ops])
TimesOperator(A::TimesOperator,B::Operator)=TimesOperator([A.ops,B])
TimesOperator(A::Operator,B::TimesOperator)=TimesOperator([A,B.ops])
TimesOperator(A::Operator,B::Operator)=TimesOperator([A,B])

function promotetimes{B}(opsin::Vector{B})
    ops=copy(opsin)
    
    for k=length(ops)-1:-1:1
        op=promotedomainspace(ops[k],rangespace(ops[k+1]))
        if op==()
            ops=ops[[1:k-1,k+1:end]]  ## remove the op
        else
            ops[k]=op
        end
    end
    
    TimesOperator(ops)
end    



function domainspace(P::TimesOperator)
    for k=length(P.ops):-1:1
        sp = domainspace(P.ops[k])
        
        if sp != AnySpace()
            return sp
        end
    end
    
    AnySpace()
end

function rangespace(P::TimesOperator)
    for op in P.ops
        sp = rangespace(op)
        
        if sp != AnySpace()
            return sp
        end
    end
    
    AnySpace()
end

domain(P::TimesOperator)=commondomain(P.ops)


function bandindssum(P,k)
    ret=0
    for op in P
        ret+=bandinds(op)[k]::Int
    end
    ret
end

bandinds(P::TimesOperator)=(bandindssum(P.ops,1),bandindssum(P.ops,2))





##TODO: We keep this around because its faster
## need to unify, maybe implement multiplyentries! with option of overriding?
function old_addentries!{T<:Number,B}(P::TimesOperator{T,B},A::ShiftArray,kr::Range1)
    cr = columnrange(A)
    br = bandinds(P)

    kre=kr[1]:(kr[end]+bandindssum(slice(P.ops,1:length(P.ops)-1),2))

    Z = ShiftArray(P.ops[end],kre,Range1(br...))
    
    for j=length(P.ops)-1:-1:2
        krr=kr[1]:(kr[end]+bandindssum(slice(P.ops,1:j-1),2))
        multiplyentries!(P.ops[j],Z,krr)
    end
    
    multiplyentries!(P.ops[1],Z,kr)    
    
    for k=kr,j=max(cr[1],br[1]):min(cr[end],br[end])
        A[k,j] += Z[k,j]
    end
    
    A
end

function new_addentries!(P::TimesOperator,A::ShiftArray,kr::Range1)
    krl=Array(Range1,length(P.ops))
    
    krl[1]=kr
    
    for m=1:length(P.ops)-1
        br=bandinds(P.ops[m])
         krl[m+1]=(br[1] + krl[m][1]):(br[end] + krl[m][end])
    end
    
    BA=BandedArray(P.ops[end],krl[end])
    
    for m=(length(P.ops)-1):-1:1
      BA=BandedArray(P.ops[m],krl[m],krl[m+1])*BA
    end
    
#     for k=kr,j=columnrange(BA.data)
#         @inbounds a=@safastget(A,k,j)
#         @inbounds b=@safastget(BA.data,k,j)
#         @inbounds @safastset!(A,a + b,k,j)
#     end
    
    
    cr=columnrange(BA.data)
    @inbounds A.data[kr+A.rowindex,cr+A.colindex]+=BA.data.data
    
    A
end

function addentries!(P::TimesOperator,A::ShiftArray,kr::Range1)
    if all(f->typeof(f)<:USConversionOperator,P.ops[1:end-1])  ##TODO: fix hack
        old_addentries!(P,A,kr)
    else
        new_addentries!(P,A,kr)
    end
end


## Algebra: assume we promote

*(A::ConstantOperator,B::ConstantOperator)=ConstantOperator(A.c*B.c)
*(A::TimesOperator,B::TimesOperator)=promotetimes([A.ops,B.ops])
*(A::TimesOperator,B::Operator)=promotetimes([A.ops,B])
*(A::Operator,B::TimesOperator)=promotetimes([A,B.ops])
*(A::Operator,B::Operator)=promotetimes([A,B])


-(A::Operator)=ConstantOperator(-1.)*A
-(A::Operator,B::Operator)=A+(-B)

*(f::IFun,A::BandedOperator)=MultiplicationOperator(f,rangespace(A))*A
*(f::FFun,A::Operator)=MultiplicationOperator(f)*A
*(c::Number,A::Operator)=ConstantOperator(c)*A
.*(c::Number,A::Operator)=ConstantOperator(c)*A





## Operations

function *{T<:Number}(A::TimesOperator,b::Vector{T})
    ret = b
    for k=length(A.ops):-1:1
        ret = A.ops[k]*ret
    end
    
    ret
end


function *{T<:Number}(A::BandedOperator,b::Vector{T})
    n=length(b)
    
    if n>0
        m=n-bandinds(A)[1]
        BandedArray(A,1:m,1:n)*b
    else
        b
    end
end

function *(A::InfiniteOperator,b::IFun)
    dsp=domainspace(A)
    dsp==AnySpace()?IFun(A*b.coefficients,b.space):IFun(A*coefficients(b,dsp),rangespace(A))
end


*{T<:Operator}(A::Vector{T},b::IFun)=map(a->a*b,convert(Array{Any,1},A))





