

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

Base.size(P::PlusOperator,k::Integer) = size(first(P.ops),k)


PlusOperator{T,UT<:Number,VT<:Number}(opsin::Vector{Operator{T}},bi::Tuple{UT,VT}) =
    PlusOperator{T,typeof(bi)}(opsin,bi)

bandinds(P::PlusOperator) = P.bandinds

israggedbelow(P::PlusOperator) = isbandedbelow(P) || all(israggedbelow,P.ops)

for (OP,mn) in ((:colstart,:min),(:colstop,:max),(:rowstart,:min),(:rowstop,:max))
    defOP = parse("default_"*string(OP))
    @eval function $OP(P::PlusOperator,k::Integer)
        if isbanded(P)
            $defOP(P,k)
        else
            mapreduce(op->$OP(op,k),$mn,P.ops)
        end
    end
end

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

function Base.convert{T}(::Type{Operator{T}},P::PlusOperator)
    if T==eltype(P)
        P
    else
        PlusOperator{T,typeof(P.bandinds)}(P.ops,P.bandinds)
    end
end

function promoteplus{T}(opsin::Vector{Operator{T}})
    ops=Vector{Operator{T}}()
    # prune zero ops
    for op in opsin
        if !iszeroop(op)
            push!(ops,op)
        end
    end
    PlusOperator(promotespaces(ops))
end

for OP in (:domainspace,:rangespace)
    @eval $OP(P::PlusOperator) = $OP(first(P.ops))
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




Base.stride(P::PlusOperator)=mapreduce(stride,gcd,P.ops)


function getindex{T}(P::PlusOperator{T},k::Integer...)
    ret=P.ops[1][k...]::T
    for op in rest(P.ops,2)
        ret+=op[k...]::T
    end
    ret
end



Base.convert{T,PP<:PlusOperator}(::Type{BandedMatrix},P::SubOperator{T,PP}) =
    banded_convert_axpy!(P)   # use axpy! to copy

Base.convert{T,PP<:PlusOperator}(::Type{BandedBlockBandedMatrix},P::SubOperator{T,PP}) =
    bandedblockbanded_convert_axpy!(P)   # use axpy! to copy

Base.convert{T,PP<:PlusOperator}(::Type{Matrix},P::SubOperator{T,PP}) =
    matrix_convert_axpy!(P)   # use axpy! to copy


function BLAS.axpy!{T,PP<:PlusOperator}(α,P::SubOperator{T,PP},A::AbstractMatrix)
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
            $OP(convert(Operator{promote_type(eltype(A),eltype(c))},c),A)
        $OP(A::Operator,c::Union{UniformScaling,Number}) =
            $OP(A,convert(Operator{promote_type(eltype(A),eltype(c))},c))
    end
end



## Times Operator

immutable ConstantTimesOperator{T,B,BT} <: Operator{BT}
    λ::T
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
    ConstantTimesOperator(c*op.λ,op.op)
ConstantTimesOperator(c::Number,op::ConstantTimesOperator) =
    ConstantTimesOperator(c*op.λ,op.op)


for OP in (:domainspace,:rangespace,:bandinds,:bandwidth,:isbanded,
           :isafunctional,:isbandedblockbanded,:israggedbelow)
    @eval $OP(C::ConstantTimesOperator) = $OP(C.op)
end
Base.size(C::ConstantTimesOperator,k::Integer) = size(C.op,k)
bandinds(C::ConstantTimesOperator,k::Integer) = bandinds(C.op,k)
choosedomainspace(C::ConstantTimesOperator,sp::Space) = choosedomainspace(C.op,sp)


for OP in (:promotedomainspace,:promoterangespace),SP in (:UnsetSpace,:Space)
    @eval function $OP(C::ConstantTimesOperator,k::$SP)
            op=$OP(C.op,k)
            # TODO: This assumes chnanging domainspace can't change the type
            ConstantTimesOperator{eltype(C.λ),typeof(op),eltype(C)}(C.λ,op)
    end
end


function Base.convert{T}(::Type{Operator{T}},C::ConstantTimesOperator)
    if T==eltype(C)
        C
    else
        op=convert(Operator{T},C.op)
        ret=ConstantTimesOperator{typeof(C.λ),typeof(op),T}(C.λ,op)
        ret
    end
end

getindex(P::ConstantTimesOperator,k::Integer...) =
    P.λ*P.op[k...]

Base.convert{T,OP<:ConstantTimesOperator}(::Type{BandedMatrix},S::SubOperator{T,OP}) =
    banded_convert_axpy!(S)

BLAS.axpy!{T,OP<:ConstantTimesOperator}(α,S::SubOperator{T,OP},A::AbstractMatrix) =
    unwrap_axpy!(α*parent(S).λ,S,A)





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
            @assert spacescompatible(domainspace(ops[k]),rangespace(ops[k+1]))
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
        ConstantOperator(1.0,dsp)
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

israggedbelow(P::TimesOperator) = isbandedbelow(P) || all(israggedbelow,P.ops)

Base.stride(P::TimesOperator) = mapreduce(stride,gcd,P.ops)

for OP in (:rowstart,:rowstop)
    defOP=parse("default_"*string(OP))
    @eval function $OP(P::TimesOperator,k::Integer)
        if isbanded(P)
            return $defOP(P,k)
        end
        for j=eachindex(P.ops)
            k=$OP(P.ops[j],k)
        end
        k
    end
end

for OP in (:colstart,:colstop)
    defOP=parse("default_"*string(OP))
    @eval function $OP(P::TimesOperator,k::Integer)
        if isbanded(P)
            return $defOP(P,k)
        end
        for j=reverse(eachindex(P.ops))
            k=$OP(P.ops[j],k)
        end
        k
    end
end

getindex(P::TimesOperator,k::Integer,j::Integer) = P[k:k,j:j][1,1]
function getindex(P::TimesOperator,k::Integer)
    @assert isafunctional(P)
    P[1:1,k:k][1,1]
end



for (STyp,Zer) in ((:BandedMatrix,:bzeros),(:Matrix,:zeros),
                    (:BandedBlockBandedMatrix,:bbbzeros),
                    (:RaggedMatrix,:rzeros))
    @eval function Base.convert{T,TO<:TimesOperator}(::Type{$STyp},
                        S::SubOperator{T,TO,Tuple{UnitRange{Int},UnitRange{Int}}})
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


        # find optimal truncations for each operator
        # by finding the non-zero entries
        krlin = Array(Union{Int,Infinity{Bool}},
                    length(P.ops),2)

        krlin[1,1],krlin[1,2]=kr[1],kr[end]
        for m=1:length(P.ops)-1
            krlin[m+1,1]=rowstart(P.ops[m],krlin[m,1])
            krlin[m+1,2]=rowstop(P.ops[m],krlin[m,2])
        end
        krlin[end,1]=max(krlin[end,1],colstart(P.ops[end],jr[1]))
        krlin[end,2]=min(krlin[end,2],colstop(P.ops[end],jr[end]))
        for m=length(P.ops)-1:-1:2
            krlin[m,1]=max(krlin[m,1],colstart(P.ops[m],krlin[m+1,1]))
            krlin[m,2]=min(krlin[m,2],colstop(P.ops[m],krlin[m+1,2]))
        end


        krl = Matrix{Int}(krlin)

        # Check if any range is invalid, in which case return zero
        for m=1:length(P.ops)
            if krl[m,1]>krl[m,2]
                return $Zer(S)
            end
        end



        # The following returns a banded Matrix with all rows
        # for large k its upper triangular
        BA=$STyp(P.ops[end][krl[end,1]:krl[end,2],jr])
        for m=(length(P.ops)-1):-1:1
            BA=$STyp(P.ops[m][krl[m,1]:krl[m,2],krl[m+1,1]:krl[m+1,2]])*BA
        end

        BA
    end
end


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

            rs=rangespace(A)
            if isambiguous(rs)
                error("Assign spaces to $A before multiplying.")
            end

            n=size(b,1)
            ret=if n>0
                A[FiniteRange,1:n]*b
            else
                b
            end
            Fun(ret,rs)
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

*{F<:Fun}(A::Operator,b::Vector{F}) = A*Fun(b)

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
        demat([A*bk  for bk in b])
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


function promotedomainspace{T}(P::PlusOperator{T},sp::Space,cursp::Space)
    if sp==cursp
        P
    else
        promoteplus(Operator{T}[promotedomainspace(op,sp) for op in P.ops])
    end
end


function choosedomainspace(P::PlusOperator,sp::Space)
    ret=UnsetSpace()
    for op in P.ops
        sp2=choosedomainspace(op,sp)
        if !isa(sp2,AmbiguousSpace)  # we will ignore this result in hopes another opand
                                     # tells us a good space
            ret=union(ret,sp2)
        end
    end
    ret
end



function promotedomainspace(P::TimesOperator,sp::Space,cursp::Space)
    if sp==cursp
        P
    elseif length(P.ops)==2
        P.ops[1]*promotedomainspace(P.ops[end],sp)
    else
        promotetimes([P.ops[1:end-1];promotedomainspace(P.ops[end],sp)])
    end
end



function choosedomainspace(P::TimesOperator,sp::Space)
    for op in P.ops
        sp=choosedomainspace(op,sp)
    end
    sp
end
