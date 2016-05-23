

export PlusOperator,TimesOperator




##PlusFunctional

immutable PlusFunctional{T} <: Functional{T}
    ops::Vector{Functional{T}}
end

for TYP in (:Operator,:Functional)
    @eval function Base.convert{T}(::Type{$TYP{T}},P::PlusFunctional)
        if T==eltype(P)
            P
        else
            PlusFunctional{T}(P.ops)
        end
    end
end



function getindex{T}(op::PlusFunctional{T},k::Integer)
    ret = op.ops[1][k]::T
    for j=2:length(op.ops)
        ret+=op.ops[j][k]::T
    end
    ret::T
end


function getindex{T}(op::PlusFunctional{T},k::Range)
    ret = convert(Vector{T},op.ops[1][k])
    for j=2:length(op.ops)
        ret+=convert(Vector{T},op.ops[j][k])
    end
    ret
end


datalength(C::PlusFunctional)=mapreduce(datalength,max,C.ops)

promotedomainspace{T}(C::PlusFunctional{T},sp::Space)=PlusFunctional(Functional{T}[promotedomainspace(c,sp) for c in C.ops])


immutable PlusOperator{T} <: BandedOperator{T}
    ops::Vector{BandedOperator{T}}
end

#Base.convert{OT<:PlusOperator}(::Type{OT},P::OT)=P
for TYP in (:Operator,:BandedOperator)
    @eval function Base.convert{T}(::Type{$TYP{T}},P::PlusOperator)
        if T==eltype(P)
            P
        else
            PlusOperator{T}(P.ops)
        end
    end
end

promoteplus{T}(ops::Vector{BandedOperator{T}})=PlusOperator(promotespaces(ops))
promoteplus{T}(ops::Vector{Functional{T}})=PlusFunctional(promotespaces(ops))




for (PLUS,TYP,ZER) in ((:PlusFunctional,:Functional,:ZeroFunctional),
                       (:PlusOperator,:BandedOperator,:ZeroOperator))
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



function rangespace(P::PlusFunctional)
    for op in P.ops
        sp = rangespace(op)

        if !isa(sp,ConstantSpace{AnyDomain})
            return sp
        end
    end

    ConstantSpace()
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


function getindex{T}(P::PlusOperator{T},k::Integer,j::Integer)
    ret=P.ops[1][k,j]::T
    for op in rest(P.ops,2)
        ret+=op[k,j]::T
    end
    ret
end



Base.copy{T,PP<:PlusOperator}(P::SubBandedMatrix{T,PP}) =
    copy_axpy!(P)   # use axpy! to copy


function BLAS.axpy!{T,PP<:PlusOperator}(α,P::SubBandedMatrix{T,PP},A::AbstractMatrix)
    for op in parent(P).ops
        BLAS.axpy!(α,sub(op,P.indexes[1],P.indexes[2]),A)
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




# We need to support A+1 in addition to A+I primarily for matrix case: A+eye(2)
for OP in (:+,:-,:(.+),:(.-))
    @eval begin
        $OP(c::Union{UniformScaling,Number},A::Operator)=$OP(convert(Operator{mat_promote_type(eltype(A),eltype(c))},c),A)
        $OP(A::Operator,c::Union{UniformScaling,Number})=$OP(A,convert(Operator{mat_promote_type(eltype(A),eltype(c))},c))
    end
end



## Times Operator

immutable ConstantTimesFunctional{T,B<:Functional} <: Functional{T}
    c::T
    functional::B
    ConstantTimesFunctional(c,op)=new(c,op)
end

ConstantTimesFunctional(c::Number,op::Functional)=ConstantTimesFunctional{promote_type(typeof(c),eltype(op)),typeof(op)}(c,op)

for OP in (:domainspace,:rangespace)
    @eval $OP(C::ConstantTimesFunctional)=$OP(C.functional)
end

for TYP in (:Operator,:Functional)
    @eval function Base.convert{T}(::Type{$TYP{T}},P::ConstantTimesFunctional)
        if T==eltype(P)
            P
        else
            ConstantTimesFunctional(convert(T,P.c),convert($TYP{T},P.functional))
        end
    end
end


getindex(op::ConstantTimesFunctional,k)=op.c*op.functional[k]
datalength(C::ConstantTimesFunctional)=datalength(C.functional)
promotedomainspace(C::ConstantTimesFunctional,sp::Space)=ConstantTimesFunctional(C.c,promotedomainspace(C.functional,sp))






type TimesFunctional{T,A<:Functional,B<:BandedOperator} <: Functional{T}
    functional::A
    op::B
end

promotedomainspace(C::TimesFunctional,sp::Space)=C.functional*promotedomainspace(C.op,sp)

domainspace(T::TimesFunctional)=domainspace(T.op)
rangespace(C::TimesFunctional)=rangespace(C.functional)

datalength(C::TimesFunctional)=datalength(C.functional)+bandinds(C.op,2)


TimesFunctional{T,V}(A::Functional{T},B::BandedOperator{V})=TimesFunctional{promote_type(T,V),typeof(A),typeof(B)}(A,B)

getindex(f::TimesFunctional,j::Integer) =
    f[j:j][1]


function getindex(f::TimesFunctional,jr::Range)#j is columns
    B=f.op
    bi=bandinds(B)

    r=zeros(eltype(f),length(jr))

    k1=max(first(jr)-bi[end],1)
    func=f.functional[k1:last(jr)-bi[1]]
    for j in jr, k=max(j-bi[end],1):j-bi[1]
        @inbounds r[j-jr[1]+1]+=func[k-k1+1]*B[k,j]
    end
    r
end


immutable ConstantTimesOperator{T,B,BT} <: BandedOperator{BT}
    c::T
    op::B
    ConstantTimesOperator(c,op)=new(c,op)
end
function ConstantTimesOperator{TT<:Number}(c::Number,op::Operator{TT})
    T=promote_type(typeof(c),eltype(op))
    B=convert(BandedOperator{T},op)
    ConstantTimesOperator{T,typeof(B),T}(c,B)
end
function ConstantTimesOperator{BM<:BandedMatrix}(c::Number,op::Operator{BM})
    BT=eltype(BM)
    T=promote_type(typeof(c),BT)

    B=convert(BandedOperator{BandedMatrix{T}},op)
    ConstantTimesOperator{T,typeof(B),BandedMatrix{T}}(c,B)
end

ConstantTimesOperator{T,B,BT}(c::Number,op::ConstantTimesOperator{T,B,BandedMatrix{BT}})=ConstantTimesOperator(c*op.c,op.op)
ConstantTimesOperator(c::Number,op::ConstantTimesOperator)=ConstantTimesOperator(c*op.c,op.op)


for OP in (:domainspace,:rangespace,:bandinds)
    @eval $OP(C::ConstantTimesOperator)=$OP(C.op)
end
bandinds(C::ConstantTimesOperator,k::Integer)=bandinds(C.op,k)
choosedomainspace(C::ConstantTimesOperator,sp::Space)=choosedomainspace(C.op,sp)


for OP in (:promotedomainspace,:promoterangespace),SP in (:AnySpace,:UnsetSpace,:Space)
    @eval function $OP(C::ConstantTimesOperator,k::$SP)
            op=$OP(C.op,k)
            # TODO: This assumes chnanging domainspace can't change the type
            ConstantTimesOperator{eltype(C.c),typeof(op),eltype(C)}(C.c,op)
    end
end

for TYP in (:Operator,:BandedOperator)
    @eval function Base.convert{T}(::Type{$TYP{T}},C::ConstantTimesOperator)
        if T==eltype(C)
            C
        else
            op=convert(BandedOperator{T},C.op)
            ret=ConstantTimesOperator{typeof(C.c),typeof(op),T}(C.c,op)
            ret
        end
    end
end


getindex(P::ConstantTimesOperator,k::Integer,j::Integer) =
    P.c*P.op[k,j]

BLAS.copy{T,OP<:ConstantTimesOperator}(S::SubBandedMatrix{T,OP}) =
    copy_axpy!(S)

BLAS.axpy!{T,OP<:ConstantTimesOperator}(α,S::SubBandedMatrix{T,OP},A::AbstractMatrix) =
    unwrap_axpy!(α*parent(S).c,S,A)





immutable TimesOperator{T} <: BandedOperator{T}
    ops::Vector{BandedOperator{T}}

    function TimesOperator(ops::Vector{BandedOperator{T}})
        hastimes=false
        for k=1:length(ops)-1
            @assert domainspace(ops[k])==AnySpace() || rangespace(ops[k+1])==AnySpace() || spacescompatible(domainspace(ops[k]),rangespace(ops[k+1]))
            hastimes=hastimes||isa(ops[k],TimesOperator)
        end

        if hastimes
            newops=Array(BandedOperator{T},0)
            for op in ops
               if isa(op,TimesOperator)
                    for op2 in op.ops
                        push!(newops,op2)
                    end
                else
                    push!(newops,op)
                end
            end
            ops=newops
        end


        new(ops)
    end
end

TimesOperator{T}(ops::Vector{BandedOperator{T}})=TimesOperator{T}(ops)
TimesOperator{OT<:Operator}(ops::Vector{OT})=TimesOperator(convert(Vector{BandedOperator{eltype(OT)}},ops))

TimesOperator(A::TimesOperator,B::TimesOperator)=TimesOperator(BandedOperator{promote_type(eltype(A),eltype(B))}[A.ops...,B.ops...])
TimesOperator(A::TimesOperator,B::BandedOperator)=TimesOperator(BandedOperator{promote_type(eltype(A),eltype(B))}[A.ops...,B])
TimesOperator(A::BandedOperator,B::TimesOperator)=TimesOperator(BandedOperator{promote_type(eltype(A),eltype(B))}[A,B.ops...])
TimesOperator(A::BandedOperator,B::BandedOperator)=TimesOperator(BandedOperator{promote_type(eltype(A),eltype(B))}[A,B])


==(A::TimesOperator,B::TimesOperator)=A.ops==B.ops


# Base.convert{OT<:TimesOperator}(::Type{OT},P::OT)=P
for TYP in (:Operator,:BandedOperator)
    @eval begin
        function Base.convert{T}(::Type{$TYP{T}},P::TimesOperator)
            if T==eltype(P)
                P
            else
                TimesOperator(BandedOperator{T}[P.ops...])
            end
        end
    end
end



function promotetimes{B<:BandedOperator}(opsin::Vector{B},dsp)
    ops=Array(BandedOperator{mapreduce(eltype,promote_type,opsin)},0)

    for k=length(opsin):-1:1
        if !isa(opsin[k],Conversion)
            op=promotedomainspace(opsin[k],dsp)
            if op==()
                # do nothing
            elseif isa(op,TimesOperator)
                for j=length(op.ops):-1:1
                    push!(ops,op.ops[j])
                end
                dsp=rangespace(op)
            else
                push!(ops,op)
                dsp=rangespace(op)
            end
        end
    end
    if isempty(ops)
        SpaceOperator(ConstantOperator(1.0),dsp,dsp)
    elseif length(ops)==1
        first(ops)
    else
        TimesOperator(reverse!(ops))  # change order in TImesOperator if this is slow
    end
end

promotetimes{B<:BandedOperator}(opsin::Vector{B})=promotetimes(opsin,domainspace(last(opsin)))



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


getindex(P::TimesOperator,k::Integer,j::Integer) = P[k:k,j:j][1,1]

function Base.copy{T,TO<:TimesOperator}(sub::SubBandedMatrix{T,TO,Tuple{UnitRange{Int},UnitRange{Int}}})
    P=parent(sub)
    kr,jr=parentindexes(sub)

    @assert length(P.ops)≥2
    if size(sub,1)==0
        return A
    end


    krl=Array(Int,length(P.ops),2)

    krl[1,1],krl[1,2]=kr[1],kr[end]

    for m=1:length(P.ops)-1
        br=bandinds(P.ops[m])
        krl[m+1,1]=max(1-mod(kr[1],1),br[1] + krl[m,1])  # no negative
        krl[m+1,2]=br[end] + krl[m,2]
    end

    # The following returns a banded Matrix with all rows
    # for large k its upper triangular
    BA=P.ops[end][krl[end,1]:krl[end,2],jr]
    for m=(length(P.ops)-1):-1:1
        BA=P.ops[m][krl[m,1]:krl[m,2],krl[m+1,1]:krl[m+1,2]]*BA
    end

    BA
end


# function addentries!(P::TimesOperator,A,kr::Range,::Colon)
#     @assert length(P.ops)≥2
#     if length(kr)==0
#         return A
#     end

#     st=step(kr)

#     krl=Array(Int,length(P.ops),2)

#     krl[1,1],krl[1,2]=kr[1],kr[end]

#     for m=1:length(P.ops)-1
#         br=bandinds(P.ops[m])
#         krl[m+1,1]=max(st-mod(kr[1],st),br[1] + krl[m,1])  # no negative
#         krl[m+1,2]=br[end] + krl[m,2]
#     end

#     # The following returns a banded Matrix with all rows
#     # for large k its upper triangular
#     BA=slice(P.ops[end],krl[end,1]:st:krl[end,2],:)
#     for m=(length(P.ops)-1):-1:2
#         BA=slice(P.ops[m],krl[m,1]:st:krl[m,2],:)*BA
#     end

#     # Write directly to A, shifting by rows and columns
#     # See subview in Operator.jl for these definitions
#     P1=slice(P.ops[1],krl[1,1]:st:krl[1,2],:)

#     firstjr=max(st-mod(kr[1],st),kr[1]+bandinds(P,1))
#     ri,ci=first(kr)-st,firstjr-st
#     bmultiply!(A,P1,BA,ri,ci,st,st)
# end




## Algebra: assume we promote


for OP in (:(Base.ctranspose),:(Base.transpose))
    @eval $OP(A::TimesOperator)=TimesOperator(reverse!(map($OP,A.ops)))
end

*(A::TimesOperator,B::TimesOperator)=promotetimes(BandedOperator{promote_type(eltype(A),eltype(B))}[A.ops...,B.ops...])
function *(A::TimesOperator,B::BandedOperator)
    if isconstop(B)
        promotedomainspace(A*convert(Number,B),domainspace(B))
    else
        promotetimes(BandedOperator{promote_type(eltype(A),eltype(B))}[A.ops...,B])
    end
end
function *(A::BandedOperator,B::TimesOperator)
    if isconstop(A)
        promoterangespace(convert(Number,A)*B,rangespace(A))
    else
        promotetimes(BandedOperator{promote_type(eltype(A),eltype(B))}[A,B.ops...])
    end
end
function *(A::BandedOperator,B::BandedOperator)
    if isconstop(A)
        promoterangespace(convert(Number,A)*B,rangespace(A))
    elseif isconstop(B)
        promotedomainspace(A*convert(Number,B),domainspace(B))
    else
        promotetimes(BandedOperator{promote_type(eltype(A),eltype(B))}[A,B])
    end
end

# Conversions we always assume are intentional: no need to promote

*{TO1<:TimesOperator,TO<:TimesOperator}(A::ConversionWrapper{TO1},B::ConversionWrapper{TO})=ConversionWrapper(TimesOperator(A.op,B.op))
*{TO<:TimesOperator}(A::ConversionWrapper{TO},B::Conversion)=ConversionWrapper(TimesOperator(A.op,B))
*{TO<:TimesOperator}(A::Conversion,B::ConversionWrapper{TO})=ConversionWrapper(TimesOperator(A,B.op))

*(A::Conversion,B::Conversion)=ConversionWrapper(TimesOperator(A,B))
*(A::Conversion,B::TimesOperator)=TimesOperator(A,B)
*(A::TimesOperator,B::Conversion)=TimesOperator(A,B)
*(A::BandedOperator,B::Conversion)=isconstop(A)?promoterangespace(convert(Number,A)*B,rangespace(A)):TimesOperator(A,B)
*(A::Conversion,B::BandedOperator)=isconstop(B)?promotedomainspace(A*convert(Number,B),domainspace(B)):TimesOperator(A,B)


-(A::Operator)=ConstantTimesOperator(-1,A)
-(A::Operator,B::Operator)=A+(-B)

*(f::Fun,A::BandedOperator)=TimesOperator(Multiplication(f,rangespace(A)),A)

for OP in (:*,:.*)
    @eval begin
        $OP(c::Number,A::BandedOperator)=c==1?A:(c==0?ZeroOperator(domainspace(A),rangespace(A)):ConstantTimesOperator(c,A))
        $OP(A::BandedOperator,c::Number)=c==1?A:(c==0?ZeroOperator(domainspace(A),rangespace(A)):ConstantTimesOperator(c,A))
    end
end

/(B::BandedOperator,c::Number)=c==1?B:ConstantTimesOperator(1.0/c,B)
/(B::BandedOperator,c::Fun)=(1.0/c)*B





## Operations

for TYP in (:Vector,:Matrix)
    @eval begin
        function *(A::TimesOperator,b::$TYP)
            ret = b
            for k=length(A.ops):-1:1
                ret = A.ops[k]*ret
            end

            ret
        end


        function *(A::BandedOperator,b::$TYP)
            n=size(b,1)

            ret=if n>0
                BandedMatrix(A,:,1:n)*b
            else
                b
            end

            rs=rangespace(A)
            isambiguous(rs)?ret:Fun(ret,rs)
        end
    end
end

function *(A::InfiniteOperator,b::Fun)
    dsp=domainspace(A)
    if isambiguous(dsp)
        promotedomainspace(A,space(b))*b
    else
        A*coefficients(b,dsp)
    end
end

#=
function *(A::PlusOperator,b::Fun)
    dsp=conversion_type(domainspace(A),space(b))
    A=promotedomainspace(A,dsp)
    ret = A.ops[1]*b
    for k=2:length(A.ops)
        ret += A.ops[k]*b
    end
    Fun(coefficients(ret),rangespace(A))
end
=#

for TYP in (:TimesOperator,:BandedOperator,:InfiniteOperator)
    @eval function *{F<:Fun}(A::$TYP,b::Matrix{F})
        @assert size(b,1)==1
        C=A*coefficients(vec(b),domainspace(A))
        reshape(C,1,length(space(C)))
    end
end

*{T<:Operator}(A::Vector{T},b::Fun)=map(a->a*b,convert(Array{Any,1},A))







## promotedomain

for T in (:AnySpace,:Space)
    @eval begin
        function promotedomainspace{T}(P::PlusOperator{T},sp::Space,cursp::$T)
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
            ret=union(ret,sp2)
        end
    end
    ret
end



for T in (:AnySpace,:Space)
    @eval begin
        function promotedomainspace(P::TimesOperator,sp::Space,cursp::$T)
            if sp==cursp
                P
            elseif length(P.ops)==2
                P.ops[1]*promotedomainspace(P.ops[end],sp)
            else
                promotetimes([P.ops[1:end-1];promotedomainspace(P.ops[end],sp)])
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
