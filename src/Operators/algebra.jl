

export PlusOperator,TimesOperator




##PlusFunctional

immutable PlusFunctional{T<:Number,B<:Functional} <: Functional{T} 
    ops::Vector{B}
end


Base.getindex(op::PlusFunctional,k::Range)=mapreduce(o->o[k],+,op.ops)




immutable PlusOperator{T<:Number,B<:BandedOperator} <: BandedOperator{T} 
    ops::Vector{B}
end


typeofdata{T<:Number}(::Operator{T})=T


# constructor for either one
function PlusFunctionalOperator{B}(PF,::Type{B},ops)
    T = Float64
    
    for op in ops
        if typeofdata(op) == Complex{Float64}
            T = Complex{Float64}
        end
    end
    
    pops=promotespaces(ops)
    PF{T,B}(pops)
end
    
## TODO: figure out how to do this with for over a tuple
PlusFunctional{B<:Functional}(ops::Vector{B})=PlusFunctionalOperator(PlusFunctional,B,ops)
PlusOperator{B<:BandedOperator}(ops::Vector{B})=PlusFunctionalOperator(PlusOperator,B,ops)



for SS in (:(PlusFunctional,Functional),:(PlusOperator,BandedOperator))
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
        +(A::($SS[1]),B::($SS[1]),C::($SS[1]))=$SS[1]([A.ops,B.ops,C.ops])        
        +(A::($SS[1]),B::($SS[2]))=$SS[1]([A.ops,B])
        +(A::($SS[1]),B::($SS[2]),C::($SS[2]))=$SS[1]([A.ops,B,C])        
        +(A::($SS[2]),B::($SS[1]))=$SS[1]([A,B.ops])
        +{T}(A::($SS[2]{T}),B::($SS[2]{T}))=$SS[1]($SS[2]{T}[A,B])
        +{T}(A::($SS[2]{T}),B::($SS[2]{T}),C::($SS[2]{T}))=$SS[1]($SS[2]{T}[A,B,C])        
        +(A::($SS[2]),B::($SS[2]))=$SS[1]($SS[2][A,B])        
        +(A::($SS[2]),B::($SS[2]),C::($SS[2]))=$SS[1]($SS[2][A,B,C])                
        #TODO: Arbitrary number of summands
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


function addentries!(P::PlusOperator,A,kr)
    for op in P.ops
        addentries!(op,A,kr)
    end
    
    A
end


+(A::Operator,f::Fun)=A+Multiplication(f,domainspace(A))
+(f::Fun,A::Operator)=Multiplication(f,domainspace(A))+A
-(A::Operator,f::Fun)=A+Multiplication(-f,domainspace(A))
-(f::Fun,A::Operator)=Multiplication(f,domainspace(A))-A



+(c::UniformScaling,A::Operator)=ConstantOperator(1.0c.λ)+A
+(A::Operator,c::UniformScaling)=A+ConstantOperator(1.0c.λ)
-(c::UniformScaling,A::Operator)=ConstantOperator(1.0c.λ)-A
-(A::Operator,c::UniformScaling)=A+ConstantOperator(-1.0c.λ)



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
    bi=bandinds(f.op)
    B=BandedArray(f.op,max((jr[1]-bi[end]),1):(jr[end]-bi[1]))
    r=zeros(T,length(jr))
    for j in jr, k=max(j-bi[end],1):j-bi[1]
        if k>=1
            r[j-jr[1]+1]+=f.functional[k]*B[k,j]
        end
    end
    r
end





type TimesOperator{T<:Number,B<:BandedOperator} <: BandedOperator{T}
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
#         if isa(ops[k],Conversion)
#             ops=ops[[1:k-1,k+1:end]]  ## remove the op
#         else        
            op=promotedomainspace(ops[k],rangespace(ops[k+1]))
            # promotedomainspace does not necessarily have correct type
            if op==()
                ops=ops[[1:k-1,k+1:end]]  ## remove the op
            elseif isa(op,TimesOperator)
                ops=[ops[1:k-1],op.ops,ops[k+1:end]]
            else
                ops=[ops[1:k-1],op,ops[k+1:end]]  
            end
#        end
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



function addentries!(P::TimesOperator,A,kr::Range1)
    krl=Array(Int,length(P.ops),2)
    
    krl[1,1],krl[1,2]=kr[1],kr[end]
    
    for m=1:length(P.ops)-1
        br=bandinds(P.ops[m])
         krl[m+1,1]=max(br[1] + krl[m,1],1)
        krl[m+1,2]=br[end] + krl[m,2]
    end
    krl
    
    BA=ShiftMatrix(P.ops[end],krl[end,1]:krl[end,2])
    for m=(length(P.ops)-1):-1:2
        BA=ShiftMatrix(P.ops[m],krl[m,1]:krl[m,2])*BA
    end
    
    bamultiply!(A,ShiftMatrix(P.ops[1],krl[1,1]:krl[1,2]),BA,kr[1]-1)
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
*(A::TimesOperator,B::BandedOperator)=promotetimes([A.ops,B])
*(A::BandedOperator,B::TimesOperator)=promotetimes([A,B.ops])
*{T}(A::BandedOperator{T},B::BandedOperator{T})=promotetimes(BandedOperator{T}[A,B])
*(A::BandedOperator,B::BandedOperator)=promotetimes(BandedOperator[A,B])


-(A::Operator)=ConstantOperator(-1.)*A
-(A::Operator,B::Operator)=A+(-B)

*(f::Fun,A::BandedOperator)=Multiplication(f,rangespace(A))*A

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
    if dsp==AnySpace()
        A=promotedomainspace(A,b.space)
        Fun(A*b.coefficients,rangespace(A))
    else
        Fun(A*coefficients(b,dsp),rangespace(A))
    end
end


*{T<:Operator}(A::Vector{T},b::Fun)=map(a->a*b,convert(Array{Any,1},A))







## promotedomain

for T in (:AnySpace,:FunctionSpace)
    @eval begin
        function promotedomainspace{T}(P::PlusOperator{T},sp::FunctionSpace,cursp::$T)
            if sp==cursp
                P
            else
                PlusOperator(BandedOperator{T}[promotedomainspace(op,sp) for op in P.ops])
            end
        end
    end
end


for T in (:AnySpace,:FunctionSpace)
    @eval begin
        function promotedomainspace(P::TimesOperator,sp::FunctionSpace,cursp::$T)
            if sp==cursp
                P
            elseif length(P.ops)==2
                P.ops[1]*promotedomainspace(P.ops[end],sp)
            else
                TimesOperator(P.ops[1:end-1])*promotedomainspace(P.ops[end],sp)
            end
        end
    end
end
