

export PlusOperator,TimesOperator



immutable PlusOperator{T,BI} <: Operator{T}
    ops::Vector{Operator{T}}
    bandinds::BI
    function PlusOperator(opsin::Vector{Operator{T}},bi::BI)
        n,m=size(first(opsin))
        for k=2:length(opsin)
            @assert size(opsin[k],1)==n && size(opsin[k],2)==m
        end
        new(opsin,bi)
    end
end

Base.size(P::PlusOperator,k...) = size(first(P.ops),k...)


PlusOperator{T,UT<:Number,VT<:Number}(opsin::Vector{Operator{T}},bi::Tuple{UT,VT}) =
    PlusOperator{T,typeof(bi)}(opsin,bi)

bandinds(P::PlusOperator) = P.bandinds

function PlusOperator(ops::Vector)
    # calculate bandinds
    b1,b2=0,0
    for op in ops
        br=bandinds(op)
        b1=min(br[1],b1)
        b2=max(br[end],b2)
    end
    PlusOperator(ops,(b1,b2))
end

#Base.convert{OT<:PlusOperator}(::Type{OT},P::OT)=P
function Base.convert{T}(::Type{Operator{T}},P::PlusOperator)
    if T==eltype(P)
        P
    else
        PlusOperator{T,eltype(P.bandinds)}(P.ops,P.bandinds)
    end
end

promoteplus{T}(ops::Vector{Operator{T}}) = PlusOperator(promotespaces(ops))

function domainspace(P::PlusOperator)
    for op in P.ops
        sp = domainspace(op)

        if !isa(sp,AnySpace)
            return sp
        end
    end

    AnySpace()
end

domain(P::PlusOperator) = commondomain(P.ops)


+(A::PlusOperator,B::PlusOperator) =
    promoteplus(Operator{promote_type(eltype(A),eltype(B))}[A.ops...,B.ops...])
+(A::PlusOperator,B::PlusOperator,C::PlusOperator) =
    promoteplus(Operator{promote_type(eltype(A),eltype(B),eltype(C))}[A.ops...,B.ops...,C.ops...])
+(A::PlusOperator,B::Operator) =
    promoteplus(Operator{promote_type(eltype(A),eltype(B))}[A.ops...,B])
+(A::PlusOperator,B::ZeroOperator) = A
+(A::PlusOperator,B::Operator,C::Operator) =
    promoteplus(Operator{promote_type(eltype(A),eltype(B),eltype(C))}[A.ops...,B,C])
+(A::Operator,B::PlusOperator) =
    promoteplus(Operator{promote_type(eltype(A),eltype(B))}[A,B.ops...])
+(A::ZeroOperator,B::PlusOperator) = B
+(A::Operator,B::Operator) =
    promoteplus(Operator{promote_type(eltype(A),eltype(B))}[A,B])
+(A::Operator,B::Operator,C::Operator) =
    promoteplus(Operator{promote_type(eltype(A),eltype(B),eltype(C))}[A,B,C])


function rangespace(P::PlusOperator)
    for op in P.ops
        sp = rangespace(op)

        if !isa(sp,AnySpace)
            return sp
        end
    end

    AnySpace()
end


Base.stride(P::PlusOperator)=mapreduce(stride,gcd,P.ops)


function getindex{T}(P::PlusOperator{T},k::Integer...)
    ret=P.ops[1][k...]::T
    for op in rest(P.ops,2)
        ret+=op[k...]::T
    end
    ret
end



Base.copy{T,PP<:PlusOperator}(P::SubBandedMatrix{T,PP}) =
    copy_axpy!(P)   # use axpy! to copy


function BLAS.axpy!{T,PP<:PlusOperator}(α,P::SubBandedMatrix{T,PP},A::AbstractMatrix)
    for op in parent(P).ops
        BLAS.axpy!(α,view(op,P.indexes[1],P.indexes[2]),A)
    end

    A
end


+(A::Operator,f::Fun) = A+Multiplication(f,domainspace(A))
+(f::Fun,A::Operator) = Multiplication(f,domainspace(A))+A
-(A::Operator,f::Fun) = A+Multiplication(-f,domainspace(A))
-(f::Fun,A::Operator) = Multiplication(f,domainspace(A))-A

for TYP in (:ZeroOperator,:Operator)
    @eval function +(A::$TYP,B::ZeroOperator)
        if spacescompatible(A,B)
            A
        else
            promotespaces(A,B)[1]
        end
    end
end
+(A::ZeroOperator,B::Operator) = B+A




# We need to support A+1 in addition to A+I primarily for matrix case: A+eye(2)
for OP in (:+,:-,:(.+),:(.-))
    @eval begin
        $OP(c::Union{UniformScaling,Number},A::Operator) =
            $OP(convert(Operator{mat_promote_type(eltype(A),eltype(c))},c),A)
        $OP(A::Operator,c::Union{UniformScaling,Number}) =
            $OP(A,convert(Operator{mat_promote_type(eltype(A),eltype(c))},c))
    end
end



## Times Operator

immutable ConstantTimesOperator{T,B,BT} <: Operator{BT}
    c::T
    op::B
    ConstantTimesOperator(c,op)=new(c,op)
end
function ConstantTimesOperator{TT<:Number}(c::Number,op::Operator{TT})
    T=promote_type(typeof(c),eltype(op))
    B=convert(Operator{T},op)
    ConstantTimesOperator{T,typeof(B),T}(c,B)
end
function ConstantTimesOperator{BM<:BandedMatrix}(c::Number,op::Operator{BM})
    BT=eltype(BM)
    T=promote_type(typeof(c),BT)

    B=convert(Operator{BandedMatrix{T}},op)
    ConstantTimesOperator{T,typeof(B),BandedMatrix{T}}(c,B)
end

ConstantTimesOperator{T,B,BT}(c::Number,op::ConstantTimesOperator{T,B,BandedMatrix{BT}}) =
    ConstantTimesOperator(c*op.c,op.op)
ConstantTimesOperator(c::Number,op::ConstantTimesOperator) =
    ConstantTimesOperator(c*op.c,op.op)


for OP in (:domainspace,:rangespace,:bandinds,:bandwidth,:isbanded,:isafunctional)
    @eval $OP(C::ConstantTimesOperator) = $OP(C.op)
end
Base.size(C::ConstantTimesOperator,k...) = size(C.op,k...)
bandinds(C::ConstantTimesOperator,k::Integer) = bandinds(C.op,k)
choosedomainspace(C::ConstantTimesOperator,sp::Space) = choosedomainspace(C.op,sp)


for OP in (:promotedomainspace,:promoterangespace),SP in (:AnySpace,:UnsetSpace,:Space)
    @eval function $OP(C::ConstantTimesOperator,k::$SP)
            op=$OP(C.op,k)
            # TODO: This assumes chnanging domainspace can't change the type
            ConstantTimesOperator{eltype(C.c),typeof(op),eltype(C)}(C.c,op)
    end
end


function Base.convert{T}(::Type{Operator{T}},C::ConstantTimesOperator)
    if T==eltype(C)
        C
    else
        op=convert(Operator{T},C.op)
        ret=ConstantTimesOperator{typeof(C.c),typeof(op),T}(C.c,op)
        ret
    end
end

getindex(P::ConstantTimesOperator,k::Integer...) =
    P.c*P.op[k...]

BLAS.copy{T,OP<:ConstantTimesOperator}(S::SubBandedMatrix{T,OP}) =
    copy_axpy!(S)

BLAS.axpy!{T,OP<:ConstantTimesOperator}(α,S::SubBandedMatrix{T,OP},A::AbstractMatrix) =
    unwrap_axpy!(α*parent(S).c,S,A)





immutable TimesOperator{T,BI} <: Operator{T}
    ops::Vector{Operator{T}}
    bandinds::BI

    function TimesOperator(ops::Vector{Operator{T}},bi::BI)
        # check compatible
        for k=1:length(ops)-1
            @assert size(ops[k],2) == size(ops[k+1],1)
        end

        # remove TimesOperators buried inside ops
        hastimes = false
        for k=1:length(ops)-1
            @assert domainspace(ops[k])==AnySpace() || rangespace(ops[k+1])==AnySpace() ||
                        spacescompatible(domainspace(ops[k]),rangespace(ops[k+1]))
            hastimes = hastimes || isa(ops[k],TimesOperator)
        end

        if hastimes
            newops=Array(Operator{T},0)
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


        new(ops,bi)
    end
end


function bandindssum(P,k)
    ret=0
    for op in P
        ret+=bandinds(op)[k]
    end
    ret
end

bandindssum(P) = (bandindssum(P,1),bandindssum(P,2))

TimesOperator{T,N1<:Number,N2<:Number}(ops::Vector{Operator{T}},bi::Tuple{N1,N2}) =
    TimesOperator{T,typeof(bi)}(ops,bi)

TimesOperator{T}(ops::Vector{Operator{T}}) = TimesOperator(ops,bandindssum(ops))
TimesOperator{OT<:Operator}(ops::Vector{OT}) =
    TimesOperator(convert(Vector{Operator{eltype(OT)}},ops),bandindssum(ops))

TimesOperator(A::TimesOperator,B::TimesOperator) =
    TimesOperator(Operator{promote_type(eltype(A),eltype(B))}[A.ops...,B.ops...])
TimesOperator(A::TimesOperator,B::Operator) =
    TimesOperator(Operator{promote_type(eltype(A),eltype(B))}[A.ops...,B])
TimesOperator(A::Operator,B::TimesOperator) =
    TimesOperator(Operator{promote_type(eltype(A),eltype(B))}[A,B.ops...])
TimesOperator(A::Operator,B::Operator) =
    TimesOperator(Operator{promote_type(eltype(A),eltype(B))}[A,B])


==(A::TimesOperator,B::TimesOperator)=A.ops==B.ops

function Base.convert{T}(::Type{Operator{T}},P::TimesOperator)
    if T==eltype(P)
        P
    else
        TimesOperator(Operator{T}[P.ops...])
    end
end



function promotetimes{B<:Operator}(opsin::Vector{B},dsp)
    ops=Array(Operator{mapreduce(eltype,promote_type,opsin)},0)

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

promotetimes{B<:Operator}(opsin::Vector{B})=promotetimes(opsin,domainspace(last(opsin)))



domainspace(P::TimesOperator)=domainspace(last(P.ops))
rangespace(P::TimesOperator)=rangespace(first(P.ops))

domain(P::TimesOperator)=commondomain(P.ops)


bandinds(P::TimesOperator) = P.bandinds
Base.stride(P::TimesOperator) = mapreduce(stride,gcd,P.ops)


getindex(P::TimesOperator,k::Integer,j::Integer) = P[k:k,j:j][1,1]
function getindex(P::TimesOperator,k::Integer)
    @assert isafunctional(P)
    P[1:1,k:k][1,1]
end

for (STyp,Zer) in ((:SubBandedMatrix,:bzeros),(:SubMatrix,:zeros))
    @eval function Base.copy{T,TO<:TimesOperator}(S::$STyp{T,TO,Tuple{UnitRange{Int},UnitRange{Int}}})
        P=parent(S)
        kr,jr=parentindexes(S)

        if maximum(kr) > size(P,1) || maximum(jr) > size(P,2) ||
            minimum(kr) < 1 || minimum(jr) < 1
            throw(BoundsError())
        end

        @assert length(P.ops)≥2
        if size(S,1)==0
            return $Zer(S)
        end


        krl = Array(promote_type(typeof(P.bandinds[1]),typeof(P.bandinds[2])),length(P.ops),2)

        krl[1,1],krl[1,2]=kr[1],kr[end]

        # find minimal row/column range starting from left
        for m=1:length(P.ops)-1
            br=bandinds(P.ops[m])
            br2=bandinds(P.ops[m+1])

            krl[m+1,1]=max(1,br[1] + krl[m,1] )  # no negative
            krl[m+1,2]=min(br[end] + krl[m,2],size(P.ops[m],2))
        end

        #find minimal row/column range starting from right
        br=bandinds(P.ops[end])
        krl[end,1]=max(krl[end,1],jr[1]-br[2])
        krl[end,2]=min(krl[end,2],jr[end]-br[1])
        for m=length(P.ops)-1:-1:2
            br=bandinds(P.ops[m])
            krl[m,1]=max(krl[m,1],krl[m+1,1]-br[2])
            krl[m,2]=min(krl[m,2],krl[m+1,2]-br[1])
        end

        krl = Matrix{Int}(krl)

        # Check if any range is invalid, in which case return zero
        for m=1:length(P.ops)
            if krl[m,1]>krl[m,2]
                return $Zer(S)
            end
        end



        # The following returns a banded Matrix with all rows
        # for large k its upper triangular
        BA=P.ops[end][krl[end,1]:krl[end,2],jr]
        for m=(length(P.ops)-1):-1:1
            BA=P.ops[m][krl[m,1]:krl[m,2],krl[m+1,1]:krl[m+1,2]]*BA
        end

        BA
    end
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

*(A::TimesOperator,B::TimesOperator) =
    promotetimes(Operator{promote_type(eltype(A),eltype(B))}[A.ops...,B.ops...])
function *(A::TimesOperator,B::Operator)
    if isconstop(B)
        promotedomainspace(A*convert(Number,B),domainspace(B))
    else
        promotetimes(Operator{promote_type(eltype(A),eltype(B))}[A.ops...,B])
    end
end
function *(A::Operator,B::TimesOperator)
    if isconstop(A)
        promoterangespace(convert(Number,A)*B,rangespace(A))
    else
        promotetimes(Operator{promote_type(eltype(A),eltype(B))}[A,B.ops...])
    end
end
function *(A::Operator,B::Operator)
    if isconstop(A)
        promoterangespace(convert(Number,A)*B,rangespace(A))
    elseif isconstop(B)
        promotedomainspace(A*convert(Number,B),domainspace(B))
    else
        promotetimes(Operator{promote_type(eltype(A),eltype(B))}[A,B])
    end
end



# Conversions we always assume are intentional: no need to promote

*{TO1<:TimesOperator,TO<:TimesOperator}(A::ConversionWrapper{TO1},B::ConversionWrapper{TO}) =
    ConversionWrapper(TimesOperator(A.op,B.op))
*{TO<:TimesOperator}(A::ConversionWrapper{TO},B::Conversion) =
    ConversionWrapper(TimesOperator(A.op,B))
*{TO<:TimesOperator}(A::Conversion,B::ConversionWrapper{TO}) =
    ConversionWrapper(TimesOperator(A,B.op))

*(A::Conversion,B::Conversion) = ConversionWrapper(TimesOperator(A,B))
*(A::Conversion,B::TimesOperator) = TimesOperator(A,B)
*(A::TimesOperator,B::Conversion) = TimesOperator(A,B)
*(A::Operator,B::Conversion) =
    isconstop(A)?promoterangespace(convert(Number,A)*B,rangespace(A)):TimesOperator(A,B)
*(A::Conversion,B::Operator) =
    isconstop(B)?promotedomainspace(A*convert(Number,B),domainspace(B)):TimesOperator(A,B)


-(A::Operator) = ConstantTimesOperator(-1,A)
-(A::Operator,B::Operator) = A+(-B)


function *(f::Fun,A::Operator)
    if isafunctional(A)
        if bandwidth(A)<Inf
            # We get a banded operator, so we take that into account
            TimesOperator(Multiplication(f,ConstantSpace()),A)
        else
            LowRankOperator(f,A)
        end
    else
        TimesOperator(Multiplication(f,rangespace(A)),A)
    end
end

for OP in (:*,:.*)
    @eval begin
        function $OP(c::Number,A::Operator)
            if c==1
                A
            elseif c==0
                ZeroOperator(domainspace(A),rangespace(A))
            else
                ConstantTimesOperator(c,A)
            end
        end
        $OP(A::Operator,c::Number) = c*A
    end
end


/(B::Operator,c::Number) = (1.0/c)*B
/(B::Operator,c::Fun) = (1.0/c)*B





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


        function *(A::Operator,b::$TYP)
            if isafunctional(A)
                return dotu(A[1:length(b)],b)
            end

            n=size(b,1)

            ret=if n>0
                A[1:min(size(A,1),n+bandwidth(A,1)),1:n]*b
            else
                b
            end

            rs=rangespace(A)
            isambiguous(rs)?ret:Fun(ret,rs)
        end
    end
end

function *(A::Operator,b::Fun)
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

for TYP in (:TimesOperator,:Operator)
    @eval function *{F<:Fun}(A::$TYP,b::Matrix{F})
        @assert size(b,1)==1
        C=A*coefficients(vec(b),domainspace(A))
        reshape(C,1,length(space(C)))
    end
end

*{T<:Operator}(A::Vector{T},b::Fun) = map(a->a*b,convert(Array{Any,1},A))



for TYP in (:Vector,:Fun,:Number)
    @eval function linsolve(A::TimesOperator,b::$TYP;kwds...)
        ret = b
        for op in A.ops
            ret = linsolve(op,ret;kwds...)
        end
        ret
    end
end






## promotedomain

for T in (:AnySpace,:Space)
    @eval begin
        function promotedomainspace{T}(P::PlusOperator{T},sp::Space,cursp::$T)
            if sp==cursp
                P
            else
                promoteplus(Operator{T}[promotedomainspace(op,sp) for op in P.ops])
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
