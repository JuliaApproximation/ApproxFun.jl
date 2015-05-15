

export PlusOperator,TimesOperator




##PlusFunctional

immutable PlusFunctional{T} <: Functional{T}
    ops::Vector{Functional{T}}
end


function Base.getindex{T}(op::PlusFunctional{T},k::Range)
    ret = convert(Vector{T},op.ops[1][k])
    for j=2:length(op.ops)
        ret+=convert(Vector{T},op.ops[j][k])
    end
    ret
end


datalength(C::PlusFunctional)=mapreduce(datalength,max,C.ops)

promotedomainspace{T}(C::PlusFunctional{T},sp::FunctionSpace)=PlusFunctional(Functional{T}[promotedomainspace(c,sp) for c in C.ops])

immutable PlusOperator{T} <: BandedOperator{T}
    ops::Vector{BandedOperator{T}}
end


Base.convert{OT<:Operator}(::Type{OT},P::PlusOperator)=PlusOperator{eltype(OT)}(P.ops)

promoteplus{T}(ops::Vector{BandedOperator{T}})=PlusOperator{T}(promotespaces(ops))
promoteplus{T}(ops::Vector{Functional{T}})=PlusFunctional{T}(promotespaces(ops))




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



for OP in (:+,:-,:(.+),:(.-))
    @eval begin
        $OP(c::UniformScaling,A::Operator)=$OP(convert(Operator{eltype(A)},c),A)
        $OP(A::Operator,c::UniformScaling)=$OP(A,convert(Operator{eltype(A)},c))    
    end
end



## Times Operator

immutable ConstantTimesFunctional{T,B<:Functional} <: Functional{T}
    c::T
    op::B
    ConstantTimesFunctional(c,op)=new(c,op)
end

ConstantTimesFunctional(c::Number,op::Functional)=ConstantTimesFunctional{promote_type(typeof(c),eltype(op)),typeof(op)}(c,op)

Base.getindex(op::ConstantTimesFunctional,k::Range)=op.c*op.op[k]
datalength(C::ConstantTimesFunctional)=datalength(C.op)
promotedomainspace(C::ConstantTimesFunctional,sp::FunctionSpace)=ConstantTimesFunctional(C.c,promotedomainspace(C.op,sp))


type TimesFunctional{T,A<:Functional,B<:BandedOperator} <: Functional{T}
    functional::A
    op::B
end

promotedomainspace(C::TimesFunctional,sp::FunctionSpace)=C.functional*promotedomainspace(C.op,sp)

for S in (:ConstantTimesFunctional,:TimesFunctional)
    @eval domainspace(T::($S))=domainspace(T.op)
end

datalength(C::TimesFunctional)=datalength(C.functional)+bandinds(C.op,2)


TimesFunctional{T,V}(A::Functional{T},B::BandedOperator{V})=TimesFunctional{promote_type(T,V),typeof(A),typeof(B)}(A,B)


function Base.getindex{T}(f::TimesFunctional{T},jr::Range)#j is columns
    bi=bandinds(f.op)
    B=subview(f.op,:,jr)
    r=zeros(T,length(jr))
    for j in jr, k=max(j-bi[end],1):j-bi[1]
        r[j-jr[1]+1]+=f.functional[k]*B[k,j]
    end
    r
end


immutable ConstantTimesOperator{T,B,BT} <: BandedOperator{BT}
    c::T
    op::B
    ConstantTimesOperator(c,op)=new(c,op)
end
function ConstantTimesOperator(c::Number,op::Operator)
    T=promote_type(typeof(c),eltype(op))
    B=convert(BandedOperator{T},op)
    ConstantTimesOperator{T,typeof(B),T}(c,B)
end
function ConstantTimesOperator{BT}(c::Number,op::Operator{BandedMatrix{BT}})
    T=promote_type(typeof(c),BT)
    B=convert(BandedOperator{BandedMatrix{T}},op)
    ConstantTimesOperator{T,typeof(B),BandedMatrix{T}}(c,B)
end


for OP in (:domainspace,:rangespace,:bandinds)
    @eval $OP(C::ConstantTimesOperator)=$OP(C.op)
end
bandinds(C::ConstantTimesOperator,k::Integer)=bandinds(C.op,k)
choosedomainspace(C::ConstantTimesOperator,sp::FunctionSpace)=choosedomainspace(C.op,sp)


for OP in (:promotedomainspace,:promoterangespace),SP in (:AnySpace,:UnsetSpace,:FunctionSpace)
    @eval begin
        $OP(C::ConstantTimesOperator,k::$SP)=ConstantTimesOperator(C.c,$OP(C.op,k))
    end
end

function Base.convert{OT<:Operator}(::Type{OT},C::ConstantTimesOperator)
    T=eltype(OT)
    op=convert(BandedOperator{T},C.op)
    ConstantTimesOperator{T,typeof(op),T}(convert(T,C.c),op)
end

function Base.convert{CT,BB,BT,OT<:Operator}(::Type{OT},C::ConstantTimesOperator{CT,BB,BandedMatrix{BT}})
    T=eltype(OT)
    op=convert(BandedOperator{T},C.op)
    ConstantTimesOperator{eltype(T),typeof(op),T}(convert(eltype(T),C.c),op)
end

function addentries!(P::ConstantTimesOperator,A,kr::Range)
    # Write directly to A, shifting by rows and columns
    # See subview in Operator.jl for these definitions
    P1=subview(P.op,kr,:)
    addentries!(P1,P.c,A,kr)
end





immutable TimesOperator{T} <: BandedOperator{T}
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


==(A::TimesOperator,B::TimesOperator)=A.ops==B.ops

Base.convert{OT<:Operator}(::Type{OT},P::TimesOperator)=TimesOperator(BandedOperator{eltype(OT)}[P.ops...])



function promotetimes{B<:BandedOperator}(opsin::Vector{B})
    ops=Array(BandedOperator{mapreduce(eltype,promote_type,opsin)},0)

    push!(ops,opsin[end])
    for k=length(opsin)-1:-1:1
        op=promotedomainspace(opsin[k],rangespace(last(ops)))
        if op==()
            # do nothing
        elseif isa(op,TimesOperator)
            for j=length(op.ops):-1:1
                push!(ops,op.ops[j])
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
    @assert length(P.ops)≥2
    if length(kr)==0
        return A
    end

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
*(A::Functional,b::Vector)=dotu(A[1:length(b)],b)
*(A::Functional,b::Fun)=promotedomainspace(A,space(b))*b.coefficients


*(c::Number,B::Functional)=ConstantTimesFunctional(c,B)
*(B::Functional,c::Number)=ConstantTimesFunctional(c,B)
/(B::Functional,c::Number)=ConstantTimesFunctional(1.0/c,B)
*(B::Functional,O::TimesOperator)=TimesFunctional(B,O)  # Needed to avoid ambiguity
*(B::Functional,O::BandedOperator)=TimesFunctional(promotedomainspace(B,rangespace(O)),O)

-(B::Functional)=ConstantTimesFunctional(-1,B)


-(A::Functional,B::Functional)=PlusFunctional([A,-B])

*{T,V}(A::TimesOperator{T},B::TimesOperator{V})=promotetimes(BandedOperator{promote_type(T,V)}[A.ops...,B.ops...])
*{T,V}(A::TimesOperator{T},B::BandedOperator{V})=promotetimes(BandedOperator{promote_type(T,V)}[A.ops...,B])
*{T,V}(A::BandedOperator{T},B::TimesOperator{V})=promotetimes(BandedOperator{promote_type(T,V)}[A,B.ops...])
*{T,V}(A::BandedOperator{T},B::BandedOperator{V})=promotetimes(BandedOperator{promote_type(T,V)}[A,B])


-(A::Operator)=ConstantTimesOperator(-1,A)
-(A::Operator,B::Operator)=A+(-B)

*(f::Fun,A::Operator)=Multiplication(f)*A

for OP in (:*,:.*)
    @eval begin
        $OP(c::Number,A::BandedOperator)=ConstantTimesOperator(c,A)
        $OP(A::BandedOperator,c::Number)=ConstantTimesOperator(c,A)
    end
end






## Operations

function *{V<:Number,T<:Number}(A::TimesOperator{V},b::Array{T})
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
    dsp=conversion_type(domainspace(A),space(b))
    A=promotedomainspace(A,dsp)
    Fun(A*coefficients(b,dsp),rangespace(A))
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






#####
# ReReOperator takes the real part of two operators
# this allows for well-posed equations
#####


immutable ReReOperator{S,V,T} <: BandedOperator{T}
    ops::(S,V)
    function ReReOperator(ops)
            #TODO: promotion
        @assert domainspace(ops[1])==domainspace(ops[2])
        @assert rangespace(ops[1])==rangespace(ops[2])
        new(ops)
    end
end


ReReOperator{S,V}(ops::@compat(Tuple{S,V}))=ReReOperator{S,V,Float64}(ops)
ReReOperator(ops1,ops2)=ReReOperator((ops1,ops2))
Base.real(S::BandedOperator,V::BandedOperator)=ReReOperator(S,V)

bandinds(R::ReReOperator)=min(2bandinds(R.ops[1],1)-1,2bandinds(R.ops[2],1)-2),max(2bandinds(R.ops[1],2)+1,2bandinds(R.ops[2],2))

domainspace(R::ReReOperator)=ReImSpace(domainspace(R.ops[1]))
rangespace(R::ReReOperator)=ArraySpace(rangespace(R.ops[1]),2)


function addentries!(R::ReReOperator,A,kr::Range)
    kr1=div(kr[1],2)+1:(iseven(kr[end])?div(kr[end],2):div(kr[end],2)+1)
    kr2=(iseven(kr[1])?div(kr[1],2):div(kr[1],2)+1):div(kr[end],2)

    B1=subview(R.ops[1],kr1,:)
    B2=subview(R.ops[2],kr2,:)


    for k=kr1,j=columnrange(R.ops[1],k)
        A[2k-1,2j-1]+=real(B1[k,j])
        A[2k-1,2j]+=-imag(B1[k,j])
    end

    for k=kr2,j=columnrange(R.ops[2],k)
        A[2k,2j-1]+=real(B2[k,j])
        A[2k,2j]+=-imag(B2[k,j])
    end

    A
end


