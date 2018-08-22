



##interlace block operators
function isboundaryrow(A,k)
    for j=1:size(A,2)
        if isafunctional(A[k,j])
            return true
        end
    end

    return false
end



domainscompatible(A::AbstractMatrix{T}) where {T<:Operator} = domainscompatible(map(domainspace,A))

function spacescompatible(A::AbstractMatrix{T}) where T<:Operator
    for k=1:size(A,1)
        if !spacescompatible(map(rangespace,A[k,:]))
            return false
        end
    end
    for k=1:size(A,2)
        if !spacescompatible(map(domainspace,A[:,k]))
            return false
        end
    end
    true
end

spacescompatible(A::AbstractVector{T}) where {T<:Operator} = spacescompatible(map(domainspace,A))

function domainspace(A::AbstractMatrix{T}) where T<:Operator
    if !spacescompatible(A)
        error("Cannot construct domainspace for $A as spaces are not compatible")
    end

    spl=map(domainspace,A[1,:])
    Space(spl)
end

function rangespace(A::AbstractVector{T}) where T<:Operator
    if !spacescompatible(A)
        error("Cannot construct rangespace for $A as domain spaces are not compatible")
    end

    spl=map(rangespace,A)
    Space(spl)
end

promotespaces(A::AbstractMatrix{T}) where {T<:Operator} = promotespaces(Matrix(A))

function promotespaces(A::Matrix{T}) where T<:Operator
    isempty(A) && return A
    ret = similar(A) #TODO: promote might have different Array type
    for j=1:size(A,2)
        ret[:,j] = promotedomainspace(A[:,j])
    end
    for k=1:size(A,1)
        ret[k,:] = promoterangespace(ret[k,:])
    end

    # do a second loop as spaces might have been inferred
    # during range space
    for j=1:size(A,2)
        ret[:,j] = promotedomainspace(ret[:,j])
    end
    ret
end


## Interlace operator

struct InterlaceOperator{T,p,DS,RS,DI,RI,BI} <: Operator{T}
    ops::Array{Operator{T},p}
    domainspace::DS
    rangespace::RS
    domaininterlacer::DI
    rangeinterlacer::RI
    bandinds::BI
end

const VectorInterlaceOperator = InterlaceOperator{T,1,DS,RS} where {T,DS,RS<:Space{D,R}} where {D,R<:AbstractVector}
const MatrixInterlaceOperator = InterlaceOperator{T,2,DS,RS} where {T,DS,RS<:Space{D,R}} where {D,R<:AbstractVector}


InterlaceOperator(ops::Array{T,p},ds,rs,di,ri,bi) where {T,p} =
    InterlaceOperator{T,p,typeof(ds),typeof(rs),
                        typeof(di),typeof(ri),typeof(bi)}(ops,ds,rs,di,ri,bi)

function InterlaceOperator(ops::AbstractMatrix{Operator{T}},ds::Space,rs::Space) where T
    # calculate bandinds TODO: generalize
    p=size(ops,1)
    dsi = interlacer(ds)
    rsi = interlacer(rs)

    if size(ops,2) == p && all(isbanded,ops) &&# only support nblocks 1 for now
            all(i->isa(i,Repeated) && i.x == 1, dsi.blocks) &&
            all(i->isa(i,Repeated) && i.x == 1, rsi.blocks)
        l,u = 0,0
        for k=1:p,j=1:p
            l=min(l,p*bandinds(ops[k,j],1)+j-k)
        end
        for k=1:p,j=1:p
            u=max(u,p*bandinds(ops[k,j],2)+j-k)
        end
    elseif p == 1 && size(ops,2) == 2 && size(ops[1],2) == 1
        # special case for example
        l,u = min(bandinds(ops[1],1),bandinds(ops[2],1)+1),bandinds(ops[2],2)+1
    else
        l,u = (1-dimension(rs),dimension(ds)-1)  # not banded
    end


    InterlaceOperator(ops,ds,rs,
                        cache(dsi),
                        cache(rsi),
                        (l,u))
end


function InterlaceOperator(ops::Vector{Operator{T}},ds::Space,rs::Space) where T
    # calculate bandinds
    p=size(ops,1)
    if all(isbanded,ops)
        l,u = 0,0
        #TODO: this code assumes an interlace strategy that might not be right
        for k=1:p
            l=min(l,p*bandinds(ops[k],1)+1-k)
        end
        for k=1:p
            u=max(u,p*bandinds(ops[k],2)+1-k)
        end
    else
        l,u = (1-dimension(rs),dimension(ds)-1)  # not banded
    end


    InterlaceOperator(ops,ds,rs,
                        cache(BlockInterlacer(tuple(blocklengths(ds)))),
                        cache(interlacer(rs)),
                        (l,u))
end

function InterlaceOperator(opsin::AbstractMatrix{Operator{T}}) where {T}
    isempty(opsin) && throw(ArgumentError("Cannot create InterlaceOperator from empty Matrix"))
    ops=promotespaces(opsin)
    InterlaceOperator(ops,domainspace(ops),rangespace(ops[:,1]))
end

function InterlaceOperator(opsin::AbstractMatrix{Operator{T}},::Type{DS},::Type{RS}) where {T,DS<:Space,RS<:Space}
    isempty(opsin) && throw(ArgumentError("Cannot create InterlaceOperator from empty Matrix"))
    ops=promotespaces(opsin)
    InterlaceOperator(ops,DS(components(domainspace(ops))),RS(rangespace(ops[:,1]).spaces))
end


function InterlaceOperator(opsin::RowVector{Operator{T}}) where {T}
    isempty(opsin) && throw(ArgumentError("Cannot create InterlaceOperator from empty Matrix"))

    ops=promotespaces(opsin)
    InterlaceOperator(ops,domainspace(ops),rangespace(ops[1]))
end

function InterlaceOperator(opsin::RowVector{Operator{T}},::Type{DS},::Type{RS}) where {T,DS<:Space,RS<:Space}
    isempty(opsin) && throw(ArgumentError("Cannot create InterlaceOperator from empty Matrix"))
    ops=promotespaces(opsin)
    InterlaceOperator(ops,DS(components(domainspace(ops))),rangespace(ops[1]))
end




InterlaceOperator(opsin::AbstractMatrix{Operator{T}},::Type{DS}) where {T,DS<:Space} =
    InterlaceOperator(opsin,DS,DS)

InterlaceOperator(opsin::AbstractMatrix,S...) =
    InterlaceOperator(Matrix{Operator{mapreduce(eltype,promote_type,opsin)}}(promotespaces(opsin)),S...)

function InterlaceOperator(opsin::Vector{Operator{T}}) where T
    ops=promotedomainspace(opsin)
    InterlaceOperator(ops,domainspace(first(ops)),rangespace(ops))
end

InterlaceOperator(ops::AbstractArray{T,p}) where {T,p} =
    InterlaceOperator(Array{Operator{mapreduce(eltype,promote_type,ops)},p}(ops))


function convert(::Type{Operator{T}},S::InterlaceOperator) where T
    if T == eltype(S)
        S
    else
        ops=Array{Operator{T}}(undef, size(S.ops)...)
        for j=1:size(S.ops,2),k=1:size(S.ops,1)
            ops[k,j]=S.ops[k,j]
        end
        InterlaceOperator(ops,domainspace(S),rangespace(S),
                            S.domaininterlacer,S.rangeinterlacer,S.bandinds)
    end
end



#TODO: More efficient to save bandinds
bandinds(M::InterlaceOperator) = M.bandinds

blockbandinds(M::InterlaceOperator) =
    (mapreduce(op->blockbandinds(op,1),min,M.ops),
     mapreduce(op->blockbandinds(op,2),max,M.ops))

isblockbanded(M::InterlaceOperator) = all(isblockbanded,M.ops)

function blockcolstop(M::InterlaceOperator,J::Integer)
    if isblockbandedbelow(M)
        Block(J - blockbandinds(M,1))
    else
        mapreduce(op->blockcolstop(op,J),max,M.ops)
    end
end



function colstop(M::InterlaceOperator{T},j::Integer) where T
#    b=bandwidth(M,1)
    if isbandedbelow(M)
        min(j+bandwidth(M,1)::Int,size(M,1))::Int
    elseif isblockbandedbelow(M)
        J=block(domainspace(M), j)::Block{1}
        blockstop(rangespace(M), blockcolstop(M,J)::Block{1})::Int
    else #assume is raggedbelow
        K = 0
        (J,jj) = M.domaininterlacer[j]
        for N = 1:size(M.ops,1)
            cs = colstop(M.ops[N,J],jj)::Int
            if cs > 0
                K = max(K,findfirst(M.rangeinterlacer,(N,cs))::Int)
            end
        end
        K
    end
end

israggedbelow(M::InterlaceOperator) = all(israggedbelow,M.ops)

getindex(op::InterlaceOperator,k::Integer,j::Integer) =
    error("Higher tensor InterlaceOperators not supported")

function getindex(op::InterlaceOperator{T,2},k::Integer,j::Integer) where {T}
    M,J = op.domaininterlacer[j]
    N,K = op.rangeinterlacer[k]
    op.ops[N,M][K,J]::T
end

# the domain is not interlaced
function getindex(op::InterlaceOperator{T,1},k::Integer,j::Integer) where T
    N,K = op.rangeinterlacer[k]
    op.ops[N][K,j]::T
end

function getindex(op::InterlaceOperator{T},k::Integer) where T
    if size(op,1) == 1
        op[1,k]
    elseif size(op,2) == 1
        op[k,1]
    else
        error("Only implemented for row/column operators.")
    end
end


findsub(cr,ν) = findall(x->x[1]==ν,cr)

function getindex(L::InterlaceOperator{T},kr::UnitRange) where T
    ret=zeros(T,length(kr))

    if size(L,1) == 1
        ds=domainspace(L)
        cr=cache(interlacer(ds))[kr]
    elseif size(L,2) == 1
        rs=rangespace(L)
        cr=cache(interlacer(rs))[kr]
    else
        error("Only implemented for row/column operators.")
    end

    for ν=1:length(L.ops)
        # indicies of ret
        ret_kr=findsub(cr,ν)

        # block indices
        if !isempty(ret_kr)
            sub_kr=cr[ret_kr[1]][2]:cr[ret_kr[end]][2]

            LinearAlgebra.axpy!(1.0,L.ops[ν][sub_kr],view(ret,ret_kr))
        end
    end
    ret
end

# overwritten for functions
# this won't work in 0.4 as expected, though the user
# should call vec anyways for 0.5 compatibility
function getindex(L::InterlaceOperator,k::Integer,j)
    if k==1 && size(L,1) == 1
        L[j]
    else
        defaultgetindex(L,k,j)
    end
end

function getindex(L::InterlaceOperator,k,j::Integer)
    if j==1 && size(L,2) == 1
        L[k]
    else
        defaultgetindex(L,k,j)
    end
end

#####
# optimized copy routine for when there is a single domainspace
# and no interlacing of the columns is necessary
# this is especially important for \
######
for TYP in (:BandedMatrix, :BlockBandedMatrix, :BandedBlockBandedMatrix, :RaggedMatrix,
                :Matrix)
    @eval begin
        function $TYP(S::SubOperator{T,InterlaceOperator{T,1,SS,PS,DI,RI,BI},
                            Tuple{UnitRange{Int},UnitRange{Int}}}) where {SS,PS,DI,RI,BI,T}
            kr,jr=parentindices(S)
            L=parent(S)

            ret=$TYP(Zeros, S)

            ds=domainspace(L)
            rs=rangespace(L)
            cr=cache(interlacer(rs))[kr]
            for ν=1:length(L.ops)
                # indicies of ret
                ret_kr=findsub(cr,ν)

                # block indices
                if !isempty(ret_kr)
                    sub_kr=cr[ret_kr[1]][2]:cr[ret_kr[end]][2]

                    LinearAlgebra.axpy!(1.0,view(L.ops[ν],sub_kr,jr),view(ret,ret_kr,:))
                end
            end
            ret
        end

        function $TYP(S::SubOperator{T,InterlaceOperator{T,2,SS,PS,DI,RI,BI},
                      Tuple{UnitRange{Int},UnitRange{Int}}}) where {SS,PS,DI,RI,BI,T}
            kr,jr=parentindices(S)
            L=parent(S)

            ret=$TYP(Zeros, S)

            if isempty(kr) || isempty(jr)
                return ret
            end

            ds=domainspace(L)
            rs=rangespace(L)
            cr=L.rangeinterlacer[kr]
            cd=L.domaininterlacer[jr]
            for ν=1:size(L.ops,1),μ=1:size(L.ops,2)
                # indicies of ret
                ret_kr=findsub(cr,ν)
                ret_jr=findsub(cd,μ)

                # block indices
                if !isempty(ret_kr) && !isempty(ret_jr)
                    sub_kr=cr[ret_kr[1]][2]:cr[ret_kr[end]][2]
                    sub_jr=cd[ret_jr[1]][2]:cd[ret_jr[end]][2]

                    LinearAlgebra.axpy!(1.0,view(L.ops[ν,μ],sub_kr,sub_jr),
                                   view(ret,ret_kr,ret_jr))
                end
            end
            ret
        end
    end
end


## Build block-by-block
function blockbanded_interlace_convert!(S,ret)
    T = eltype(S)
    KR,JR = parentindices(S)
    l,u=blockbandwidths(S)::Tuple{Int,Int}

    M = map(op -> begin
                KR_size = Block.(Int(first(KR)):min(Int(last(KR)),nblocks(op,1)))
                JR_size = Block.(Int(first(JR)):min(Int(last(JR)),nblocks(op,2)))
                BlockBandedMatrix(view(op, KR_size, JR_size))
            end, parent(S).ops)

    for J=Block(1):Block(nblocks(ret,2)),K=blockcolrange(ret,Int(J))
        Bs=view(ret,K,J)
        j = 0
        for ξ=1:size(M,2)
            k = 0
            m = 0
            for κ=1:size(M,1)
                if K.n[1] ≤ nblocks(M[κ,ξ],1) && J.n[1] ≤ nblocks(M[κ,ξ],2)
                    MKJ = M[κ,ξ][K,J]::Matrix{T}
                    n,m = size(MKJ)
                    Bs[k+1:k+n,j+1:j+m] = MKJ
                    k += n
                end
            end
            j += m
        end
    end
    ret
end

for d in (:1,:2)
    @eval BlockBandedMatrix(S::SubOperator{T,InterlaceOperator{T,$d,SS,PS,DI,RI,BI},
                          Tuple{BlockRange1,BlockRange1}}) where {SS,PS,DI,RI,BI,T} =
    blockbanded_interlace_convert!(S, BlockBandedMatrix(Zeros, S))
end





domainspace(IO::InterlaceOperator) = IO.domainspace
rangespace(IO::InterlaceOperator) = IO.rangespace

#tests whether an operator can be made into a column
iscolop(op) = isconstop(op)
iscolop(::Multiplication) = true

promotedomainspace(A::InterlaceOperator{T,1},sp::Space) where {T} =
    InterlaceOperator(map(op->promotedomainspace(op,sp),A.ops))


interlace(A::AbstractArray{T}) where {T<:Operator} = InterlaceOperator(A)

const OperatorTypes = Union{Operator,Fun,Number,UniformScaling}



operators(A::InterlaceOperator) = A.ops
operators(A::Matrix{<:Operator}) = A
operators(A::Operator) = [A]

Base.vcat(A::MatrixInterlaceOperator...) =
    InterlaceOperator(vcat(map(operators,A)...))
function _vcat(A::OperatorTypes...)
    Av = Vector{Operator{mapreduce(eltype,promote_type,A)}}()
    for a in A
        if a isa VectorInterlaceOperator
            append!(Av,a.ops)
        else
            push!(Av,a)
        end
    end
    InterlaceOperator(vnocat(Av...))
end



Base.hcat(A::Union{VectorInterlaceOperator,MatrixInterlaceOperator}...) =
    InterlaceOperator(hcat(map(A->A.ops,A)...))
_hcat(A::OperatorTypes...) = InterlaceOperator(hnocat(A...))
function _hvcat(rows::Tuple{Vararg{Int}},as::OperatorTypes...)
    # Based on Base
    nbr = length(rows)  # number of block rows
    rs = Array{Any,1}(undef, nbr)
    a = 1
    for i = 1:nbr
        rs[i] = hcat(map(op -> convert(Operator,op),as[a:a-1+rows[i]])...)
        a += rows[i]
    end
    vcat(rs...)
end

Base.vcat(A::Operator, B::OperatorTypes...) = _vcat(A, B...)
Base.hcat(A::Operator, B::OperatorTypes...) = _hcat(A, B...)
Base.hvcat(rows::Tuple{Vararg{Int}}, A::Operator, B::OperatorTypes...) =
    _hvcat(rows, A, B...)

Base.vcat(C::Union{Fun,Number,UniformScaling}, A::Operator, B::OperatorTypes...) = _vcat(C, A, B...)
Base.hcat(C::Union{Fun,Number,UniformScaling}, A::Operator, B::OperatorTypes...) = _hcat(C, A, B...)
Base.hvcat(rows::Tuple{Vararg{Int}}, C::Union{Fun,Number,UniformScaling}, A::Operator, B::OperatorTypes...) =
    _hvcat(rows, C, A, B...)

Base.vcat(D::Union{Fun,Number,UniformScaling}, C::Union{Fun,Number,UniformScaling}, A::Operator, B::OperatorTypes...) = _vcat(D, C, A, B...)
Base.hcat(D::Union{Fun,Number,UniformScaling}, C::Union{Fun,Number,UniformScaling}, A::Operator, B::OperatorTypes...) = _hcat(D, C, A, B...)
Base.hvcat(rows::Tuple{Vararg{Int}}, D::Union{Fun,Number,UniformScaling}, C::Union{Fun,Number,UniformScaling}, A::Operator, B::OperatorTypes...) =
    _hvcat(rows, D, C, A, B...)


## Convert Matrix operator to operators

Operator(M::AbstractArray{OO}) where {OO<:Operator} = InterlaceOperator(M)




function interlace_choosedomainspace(ops,sp::UnsetSpace)
    # this ensures correct dispatch for unino
    sps = Vector{Space}(
        filter(x->!isambiguous(x),map(choosedomainspace,ops)))
    if isempty(sps)
        UnsetSpace()
    else
        union(sps...)
    end
end


function interlace_choosedomainspace(ops,rs::Space)
    # this ensures correct dispatch for unino
    sps = Vector{Space}(
        filter(x->!isambiguous(x),map((op)->choosedomainspace(op,rs),ops)))
    if isempty(sps)
        UnsetSpace()
    else
        union(sps...)
    end
end


choosedomainspace(A::InterlaceOperator{T,1},rs::Space) where {T} =
    interlace_choosedomainspace(A.ops,rs)


choosedomainspace(A::InterlaceOperator{T,2},rs::Space) where {T} =
    Space([interlace_choosedomainspace(A.ops[:,k],rs) for k=1:size(A.ops,2)])
