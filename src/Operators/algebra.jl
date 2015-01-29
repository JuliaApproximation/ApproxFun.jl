

export PlusOperator,TimesOperator




##PlusFunctional

immutable PlusFunctional{T<:Number} <: Functional{T} 
    ops::Vector{Functional{T}}
end


function Base.getindex{T}(op::PlusFunctional{T},k::Range)
    ret = op.ops[1][k]::Vector{T}
    for j=2:length(op.ops)
        ret+=op.ops[j][k]::Vector{T}
    end
    ret
end




immutable PlusOperator{T<:Number} <: BandedOperator{T} 
    ops::Vector{BandedOperator{T}}
end


Base.convert{T}(::Type{BandedOperator{T}},P::PlusOperator)=PlusOperator{T}(P.ops)

promoteplus{T<:Number}(ops::Vector{BandedOperator{T}})=PlusOperator{T}(promotespaces(ops))
promoteplus{T<:Number}(ops::Vector{Functional{T}})=PlusFunctional{T}(promotespaces(ops))




for (PLUS,TYP,ZER) in ((:PlusFunctional,:Functional,:ZeroFunctional),(:PlusOperator,:BandedOperator,:ZeroOperator))
    @eval begin
        function domainspace(P::$PLUS)
            for op in P.ops
                sp = domainspace(op)
                
                if !isa(sp,AnySpace)
                    return sp
                end
            end
            
            AnySpace()
        end        
        
        domain(P::$PLUS)=commondomain(P.ops)

    
        +(A::$PLUS,B::$PLUS)=promoteplus($TYP{promote_type(eltype(A),eltype(B))}[A.ops...,B.ops...])
        +(A::$PLUS,B::$PLUS,C::$PLUS)=promoteplus($TYP{promote_type(eltype(A),eltype(B),eltype(C))}[A.ops...,B.ops...,C.ops...])        
        +(A::$PLUS,B::$TYP)=promoteplus($TYP{promote_type(eltype(A),eltype(B))}[A.ops...,B])
        +(A::$PLUS,B::$ZER)=A
        +(A::$PLUS,B::$TYP,C::$TYP)=promoteplus($TYP{promote_type(eltype(A),eltype(B),eltype(C))}[A.ops...,B,C])        
        +(A::$TYP,B::$PLUS)=promoteplus($TYP{promote_type(eltype(A),eltype(B))}[A,B.ops...])
        +(A::$ZER,B::$PLUS)=B
        +(A::$TYP,B::$TYP)=promoteplus($TYP{promote_type(eltype(A),eltype(B))}[A,B])        
        +(A::$TYP,B::$TYP,C::$TYP)=promoteplus($TYP{promote_type(eltype(A),eltype(B),eltype(C))}[A,B,C])                
        #TODO: Arbitrary number of summands
    end
end




function rangespace(P::PlusOperator)
    for op in P.ops
        sp = rangespace(op)
        
        if !isa(sp,AnySpace)
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


Base.stride(P::PlusOperator)=mapreduce(stride,gcd,P.ops)


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
+(A::ZeroOperator,::ZeroOperator)=A
+(A::BandedOperator,::ZeroOperator)=A
+(::ZeroOperator,B::BandedOperator)=B
+(A::ZeroFunctional,::ZeroFunctional)=A
+(A::Functional,::ZeroFunctional)=A
+(::ZeroFunctional,B::Functional)=B


+(c::UniformScaling,A::Operator)=ConstantOperator(c.λ)+A
+(A::Operator,c::UniformScaling)=A+ConstantOperator(c.λ)
-(c::UniformScaling,A::Operator)=ConstantOperator(c.λ)-A
-(A::Operator,c::UniformScaling)=A+ConstantOperator(-c.λ)



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


TimesFunctional{T<:Number,V<:Number}(A::Functional{T},B::BandedOperator{V})=TimesFunctional{promote_type(T,V),typeof(A),typeof(B)}(A,B)


function Base.getindex{T<:Number}(f::TimesFunctional{T},jr::Range)#j is columns
    bi=bandinds(f.op)
    B=subview(f.op,:,jr)
    r=zeros(T,length(jr))
    for j in jr, k=max(j-bi[end],1):j-bi[1]
        r[j-jr[1]+1]+=f.functional[k]*B[k,j]
    end
    r
end





type TimesOperator{T<:Number} <: BandedOperator{T}
    ops::Vector{BandedOperator{T}}
    
    function TimesOperator(ops::Vector{BandedOperator{T}})
        for k=1:length(ops)-1
            @assert domainspace(ops[k])==AnySpace() || rangespace(ops[k+1])==AnySpace() || domainspace(ops[k])==rangespace(ops[k+1])
        end
        
        new(ops)
    end
end

TimesOperator{T}(ops::Vector{BandedOperator{T}})=TimesOperator{T}(ops)

TimesOperator{T,V}(A::TimesOperator{T},B::TimesOperator{V})=TimesOperator(BandedOperator{promote_type(T,V)}[A.ops...,B.ops...])
TimesOperator{T,V}(A::TimesOperator{T},B::BandedOperator{V})=TimesOperator(BandedOperator{promote_type(T,V)}[A.ops...,B])
TimesOperator{T,V}(A::BandedOperator{T},B::TimesOperator{V})=TimesOperator(BandedOperator{promote_type(T,V)}[A,B.ops...])
TimesOperator{T,V}(A::BandedOperator{T},B::BandedOperator{V})=TimesOperator(BandedOperator{promote_type(T,V)}[A,B])


Base.convert{T}(::Type{BandedOperator{T}},P::TimesOperator)=TimesOperator(BandedOperator{T}[P.ops...])



function promotetimes{B<:BandedOperator}(opsin::Vector{B})
    ops=Array(BandedOperator{mapreduce(eltype,promote_type,opsin)},0)

    push!(ops,opsin[end]) 
    for k=length(opsin)-1:-1:1
        op=promotedomainspace(opsin[k],rangespace(last(ops)))
        if op==()
            # do nothing
        elseif isa(op,TimesOperator)
            for j=length(op):-1:1
                push!(ops,op[j])
            end
        else
            push!(ops,op)
        end
#        end
    end
    
    TimesOperator(reverse!(ops))  # change order in TImesOperator if this is slow
end    



domainspace(P::TimesOperator)=domainspace(last(P.ops))
rangespace(P::TimesOperator)=rangespace(first(P.ops))

domain(P::TimesOperator)=commondomain(P.ops)


function bandindssum(P,k)
    ret=0
    for op in P
        ret+=bandinds(op)[k]::Int
    end
    ret
end

bandinds(P::TimesOperator)=(bandindssum(P.ops,1),bandindssum(P.ops,2))
Base.stride(P::TimesOperator)=mapreduce(stride,gcd,P.ops)



function addentries!(P::TimesOperator,A,kr::Range)
   st=step(kr)

    krl=Array(Int,length(P.ops),2)
    
    krl[1,1],krl[1,2]=kr[1],kr[end]
    
    for m=1:length(P.ops)-1
        br=bandinds(P.ops[m])
        krl[m+1,1]=max(st-mod(kr[1],st),br[1] + krl[m,1])  # no negative 
        krl[m+1,2]=br[end] + krl[m,2]
    end
    
    # The following returns a banded Matrix with all rows
    # for large k its upper triangular
    BA=slice(P.ops[end],krl[end,1]:st:krl[end,2],:)
    for m=(length(P.ops)-1):-1:2
        BA=slice(P.ops[m],krl[m,1]:st:krl[m,2],:)*BA
    end
    
    # Write directly to A, shifting by rows and columns
    # See subview in Operator.jl for these definitions
    P1=slice(P.ops[1],krl[1,1]:st:krl[1,2],:)
    
    firstjr=max(st-mod(kr[1],st),kr[1]+bandinds(P,1))
    ri,ci=first(kr)-st,firstjr-st
    bamultiply!(A,P1,BA,ri,ci,st,st)
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

*{T,V}(A::TimesOperator{T},B::TimesOperator{V})=promotetimes(BandedOperator{promote_type(T,V)}[A.ops...,B.ops...])
*{T,V}(A::TimesOperator{T},B::BandedOperator{V})=promotetimes(BandedOperator{promote_type(T,V)}[A.ops...,B])
*{T,V}(A::BandedOperator{T},B::TimesOperator{V})=promotetimes(BandedOperator{promote_type(T,V)}[A,B.ops...])
*{T,V}(A::BandedOperator{T},B::BandedOperator{V})=promotetimes(BandedOperator{promote_type(T,V)}[A,B])


-(A::Operator)=ConstantOperator(-1.)*A
-(A::Operator,B::Operator)=A+(-B)

*(f::Fun,A::BandedOperator)=Multiplication(f,rangespace(A))*A

*(c::Number,A::Operator)=ConstantOperator(c)*A
.*(c::Number,A::Operator)=ConstantOperator(c)*A





## Operations

function *{T<:Number}(A::TimesOperator,b::Array{T})
    ret = b
    for k=length(A.ops):-1:1
        ret = A.ops[k]*ret
    end
    
    ret
end


function *{T<:Number}(A::BandedOperator,b::Array{T})
    n=size(b,1)
    
    if n>0
        slice(A,:,1:n)*b
    else
        b
    end
end
 ##TODO: Make * and \ consistent in return type
function *(A::InfiniteOperator,b::Fun)
    dsp=domainspace(A)
    if isa(dsp,AmbiguousSpace)
        A=promotedomainspace(A,b.space)
        Fun(A*b.coefficients,rangespace(A))
    else
        Fun(A*coefficients(b,dsp),rangespace(A))
    end
end

function *{F<:Fun}(A::InfiniteOperator,b::Array{F,2})
    @assert size(b,1)==1
    C=A*coefficients(vec(b),domainspace(A))
    rs=rangespace(A)
    ret=Array(Fun{typeof(rs),eltype(C)},1,size(C,2))
    for k=1:size(C,2)
        ret[1,k]=Fun(C[:,k],rs)
    end
    ret
end


*{T<:Operator}(A::Vector{T},b::Fun)=map(a->a*b,convert(Array{Any,1},A))







## promotedomain

for T in (:AnySpace,:FunctionSpace)
    @eval begin
        function promotedomainspace{T}(P::PlusOperator{T},sp::FunctionSpace,cursp::$T)
            if sp==cursp
                P
            else
                promoteplus(BandedOperator{T}[promotedomainspace(op,sp) for op in P.ops])
            end
        end
    end
end

function choosedomainspace(P::PlusOperator,sp)
    ret=AnySpace()
    for op in P.ops
        sp2=choosedomainspace(op,sp)
        if !isa(sp2,AmbiguousSpace)  # we will ignore this result in hopes another opand
                                     # tells us a good space
            ret=conversion_type(ret,sp2)
        end
    end
    ret
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



function choosedomainspace(P::TimesOperator,sp)
    for op in P.ops
        sp=choosedomainspace(op,sp)
    end
    sp
end
    
    
    

