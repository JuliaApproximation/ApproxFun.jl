

export PlusOperator,TimesOperator




##PlusFunctional

immutable PlusFunctional{T<:Number,B<:Functional} <: Functional{T} 
    ops::Vector{B}
end


Base.getindex(op::PlusFunctional,k::Range)=mapreduce(o->o[k],+,op.ops)




immutable PlusOperator{T<:Number,B<:Operator} <: BandedOperator{T} 
    ops::Vector{B}
end


typeofdata{T<:Number}(::Operator{T})=T


# constructor for either one
function PlusFunctionalOperator(PF,PF2,ops)
    T = Float64
    
    for op in ops
        if typeofdata(op) == Complex{Float64}
            T = Complex{Float64}
        end
    end
    
    pops=promotespaces(ops)
    PF{T,PF2}(pops)
end
    
## TODO: figure out how to do this with for over a tuple
PlusFunctional{B<:Functional}(ops::Vector{B})=PlusFunctionalOperator(PlusFunctional,Functional,ops)
PlusOperator{B<:Operator}(ops::Vector{B})=PlusFunctionalOperator(PlusOperator,Operator,ops)



for SS in (:(PlusFunctional,Functional),:(PlusOperator,Operator))
    @eval begin
        function domainspace(P::($SS[1]))
            for op in P.ops
                sp = domainspace(op)
                
                if sp != AnySpace()
                    return sp
                end
            end
            
            AnySpace()
        end        
        
        domain(P::($SS[1]))=commondomain(P.ops)

    
        +(A::($SS[1]),B::($SS[1]))=$SS[1]([A.ops,B.ops])
        +(A::($SS[1]),B::($SS[2]))=$SS[1]([A.ops,B])
        +(A::($SS[2]),B::($SS[1]))=$SS[1]([A,B.ops])
        +(A::($SS[2]),B::($SS[2]))=$SS[1]([A,B])
    end
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


+(A::Operator,f::Fun)=A+MultiplicationOperator(f)
+(f::Fun,A::Operator)=MultiplicationOperator(f)+A
-(A::Operator,f::Fun)=A+MultiplicationOperator(-f)
-(f::Fun,A::Operator)=MultiplicationOperator(f)-A



+(c::UniformScaling,A::Operator)=ConstantOperator(1.0c.λ)+A
+(A::Operator,c::UniformScaling)=A+ConstantOperator(1.0c.λ)
-(c::UniformScaling,A::Operator)=ConstantOperator(1.0c.λ)-A
-(A::Operator,c::UniformScaling)=A-ConstantOperator(1.0c.λ)



## We need this to support Op + Matrix
+(c::Number,A::Operator)=ConstantOperator(c)+A
+(A::Operator,c::Number)=A+ConstantOperator(c)
-(c::Number,A::Operator)=ConstantOperator(c)-A
-(A::Operator,c::Number)=A-ConstantOperator(c)

.+(c::Number,A::Operator)=ConstantOperator(c)+A
.+(A::Operator,c::Number)=A+ConstantOperator(c)
.-(c::Number,A::Operator)=ConstantOperator(c)-A
.-(A::Operator,c::Number)=A-ConstantOperator(c)



## Times Operator

immutable ConstantTimesFunctional{T<:Number,B<:Functional} <: Functional{T}
    c::T
    op::B
end

Base.getindex(op::ConstantTimesFunctional,k::Range1)=op.c*op.op[k]



type TimesFunctional{T<:Number,A<:Functional,B<:BandedOperator} <: Functional{T}
    functional::A
    op::B
end

for S in (:ConstantTimesFunctional,:TimesFunctional)
    @eval domainspace(T::($S))=domainspace(T.op)
end


TimesFunctional{T<:Number}(A::Functional{T},B::BandedOperator{T})=TimesFunctional{T,typeof(A),typeof(B)}(A,B)


function Base.getindex{T<:Number}(f::TimesFunctional{T},jr::Range)#j is columns
    bi=ApproxFun.bandinds(f.op)
    B=BandedArray(f.op,(jr[1]-bi[end]):(jr[end]-bi[1]))
    r=zeros(T,length(jr))
    for j in jr, k=j-bi[end]:j-bi[1]
        if k>=1
            r[j-jr[1]+1]+=f.functional[k]*B[k,j]
        end
    end
    r
end





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

## Operations
*(A::Functional,b::Vector)=dot(A[1:length(b)],b)
*(A::Functional,b::Fun)=A*b.coefficients


*(c::Number,B::Functional)=ConstantTimesFunctional(c,B)
*(B::Functional,c::Number)=ConstantTimesFunctional(c,B)
/(B::Functional,c::Number)=ConstantTimesFunctional(1.0/c,B)
*(B::Functional,O::TimesOperator)=TimesFunctional(B,O)  # Needed to avoid ambiguity
*(B::Functional,O::BandedOperator)=TimesFunctional(B,O)

-{T<:Number}(B::Functional{T})=ConstantTimesFunctional(-one(T),B)


-(A::Functional,B::Functional)=PlusFunctional([A,-B])

*(A::TimesOperator,B::TimesOperator)=promotetimes([A.ops,B.ops])
*(A::TimesOperator,B::Operator)=promotetimes([A.ops,B])
*(A::Operator,B::TimesOperator)=promotetimes([A,B.ops])
*(A::Operator,B::Operator)=promotetimes([A,B])


-(A::Operator)=ConstantOperator(-1.)*A
-(A::Operator,B::Operator)=A+(-B)

*(f::Fun,A::BandedOperator)=MultiplicationOperator(f)*A

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
 ##TODO: Make * and \ consistent in return type
function *(A::InfiniteOperator,b::Fun)
    dsp=domainspace(A)
    dsp==AnySpace()?Fun(A*b.coefficients,b.space):Fun(A*coefficients(b,dsp),rangespace(A))
end


*{T<:Operator}(A::Vector{T},b::Fun)=map(a->a*b,convert(Array{Any,1},A))





