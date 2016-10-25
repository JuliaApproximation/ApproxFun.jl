



##interlace block operators
function isboundaryrow(A,k)
    for j=1:size(A,2)
        if isafunctional(A[k,j])
            return true
        end
    end

    return false
end



domainscompatible{T<:Operator}(A::Matrix{T})=domainscompatible(map(domainspace,A))

function spacescompatible{T<:Operator}(A::Matrix{T})
    for k=1:size(A,1)
        if !spacescompatible(map(rangespace,vec(A[k,:])))
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

spacescompatible{T<:Operator}(A::Vector{T}) = spacescompatible(map(domainspace,A))

function domainspace{T<:Operator}(A::Matrix{T})
    if !spacescompatible(A)
        error("Cannot construct domainspace for $A as spaces are not compatible")
    end

    spl=map(domainspace,vec(A[1,:]))
    if spacescompatible(spl)
        ArraySpace(first(spl),length(spl))
    else
        TupleSpace(spl)
    end
end

function rangespace{T<:Operator}(A::Vector{T})
    if !spacescompatible(A)
        error("Cannot construct rangespace for $A as domain spaces are not compatible")
    end

    spl=map(rangespace,A)
    if spacescompatible(spl)
        ArraySpace(first(spl),length(spl))
    else
        TupleSpace(spl)
    end
end

promotespaces{T<:Operator}(A::AbstractMatrix{T}) = promotespaces(Matrix(A))

function promotespaces{T<:Operator}(A::Matrix{T})
    isempty(A) && return A
    A=copy(A)#TODO: promote might have different Array type
    for j=1:size(A,2)
        A[:,j]=promotedomainspace(A[:,j])
    end
    for k=1:size(A,1)
        A[k,:]=promoterangespace(vec(A[k,:]))
    end

    # do a second loop as spaces might have been inferred
    # during range space
    for j=1:size(A,2)
        A[:,j]=promotedomainspace(A[:,j])
    end
    A
end


## Interlace operator

immutable InterlaceOperator{T,p,DS,RS,DI,RI,BI} <: Operator{T}
    ops::Array{Operator{T},p}
    domainspace::DS
    rangespace::RS
    domaininterlacer::DI
    rangeinterlacer::RI
    bandinds::BI
end


InterlaceOperator{T,p}(ops::Array{T,p},ds,rs,di,ri,bi) =
    InterlaceOperator{T,p,typeof(ds),typeof(rs),
                        typeof(di),typeof(ri),typeof(bi)}(ops,ds,rs,di,ri,bi)

function InterlaceOperator{T}(ops::Matrix{Operator{T}},ds::Space,rs::Space)
    # calculate bandinds TODO: generalize
    p=size(ops,1)
    dsi = interlacer(ds)

    if p == 1  # Assume rs corresonds to a scalar space, so wrap in a TupleSpace to get right interlacing
        rsi = interlacer(TupleSpace((rs,)))
    else  # assume this is a correctly tupled
        rsi = interlacer(rs)
    end

    if size(ops,2) == p && all(isbanded,ops) &&# only support blocksize 1 for now
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


function InterlaceOperator{T}(ops::Vector{Operator{T}},ds::Space,rs::Space)
    # calculate bandinds
    p=size(ops,1)
    if all(isbanded,ops)
        l,u = 0,0
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
                        BlockInterlacer(tuple(blocklengths(ds))),
                        cache(interlacer(rs)),
                        (l,u))
end

function InterlaceOperator{T}(opsin::Matrix{Operator{T}})
    isempty(opsin) && throw(ArgumentError("Cannot create InterlaceOperator from empty Matrix"))

    ops=promotespaces(opsin)
    # TODO: make consistent
    # if its a row vector, we assume scalar
    if size(ops,1) == 1
        InterlaceOperator(ops,domainspace(ops),rangespace(ops[1]))
    else
        InterlaceOperator(ops,domainspace(ops),rangespace(ops[:,1]))
    end
end

function InterlaceOperator{T,DS<:Space,RS<:Space}(opsin::Matrix{Operator{T}},::Type{DS},::Type{RS})
    isempty(opsin) && throw(ArgumentError("Cannot create InterlaceOperator from empty Matrix"))

    ops=promotespaces(opsin)
    # TODO: make consistent
    # if its a row vector, we assume scalar
    if size(ops,1) == 1
        InterlaceOperator(ops,DS(domainspace(ops).spaces),rangespace(ops[1]))
    else
        InterlaceOperator(ops,DS(domainspace(ops).spaces),RS(rangespace(ops[:,1]).spaces))
    end
end

InterlaceOperator{T,DS<:Space}(opsin::Matrix{Operator{T}},::Type{DS}) =
    InterlaceOperator(opsin,DS,DS)

InterlaceOperator(opsin::AbstractMatrix,S...) = InterlaceOperator(Matrix{Operator{mapreduce(eltype,promote_type,opsin)}}(opsin),S...)

function InterlaceOperator{T}(opsin::Vector{Operator{T}})
    ops=promotedomainspace(opsin)
    InterlaceOperator(ops,domainspace(first(ops)),rangespace(ops))
end

InterlaceOperator{T,p}(ops::Array{T,p}) =
    InterlaceOperator(Array{Operator{mapreduce(eltype,promote_type,ops)},p}(ops))


function Base.convert{T}(::Type{Operator{T}},S::InterlaceOperator)
    if T == eltype(S)
        S
    else
        ops=Array(Operator{T},size(S.ops)...)
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

isbandedblock(M::InterlaceOperator) = all(x->isbandedblock(x)||isbanded(x),M.ops) #TODO: this is a hack...

function blockcolstop(M::InterlaceOperator,J::Integer)
    if isbandedblock(M)
        J - blockbandinds(M,1)
    else
        mapreduce(op->blockcolstop(op,J),max,M.ops)
    end
end



function colstop{T}(M::InterlaceOperator{T,1},j::Integer)
#    b=bandwidth(M,1)
    if isbandedbelow(M)
        min(j+bandwidth(M,1)::Int,size(M,1)::Int)
    elseif isbandedblock(M)
        J=block(domainspace(M),j)::Int
        blockstop(rangespace(M),blockcolstop(M,J)::Int)::Int
    else #assume is raggedbelow
#        opj,J = M.rangeinterlacer[j]
        K = 0
        for N = 1:size(M.ops,1)
            K = max(K,findfirst(M.rangeinterlacer,(N,colstop(M.ops[N],j)))::Int)
        end
        K
    end
end

function colstop{T}(M::InterlaceOperator{T,2},j::Integer)
    b=bandwidth(M,1)
    if isfinite(b)
        min(j+b,size(M,1))
    else # assume block banded
        J=block(domainspace(M),j)
        blockstop(rangespace(M),blockcolstop(M,J))
    end
end

#
# function colstop{T}(M::InterlaceOperator{T,2},j::Integer)
#     l = M.bandinds[1]
#
#     if isinf(l)
#         (J,ξ) = M.domaininterlacer[j]
#         ret = j
#
#         for K=1:size(M.ops,1)
#             cs=colstop(M.ops[K,J],ξ)
#
#             for k=j:size(M,1)
#                 (KK,κ)=M.rangeinterlacer[k]
#                 if K == KK && κ ≥ cs
#                     ret = max(ret,k)
#                     break
#                 end
#             end
#         end
#         return ret
#     else
#         return  min(size(M,1),j-l)
#     end
# end

# function colstop{T}(M::InterlaceOperator{T,1},j::Integer)
#     l = M.bandinds[1]
#
#     if isinf(l)
#         ret = j
#
#         for K=1:length(M.ops)
#             cs=colstop(M.ops[K],j)
#
#             for k=j:size(M,1)
#                 (KK,κ)=M.rangeinterlacer[k]
#                 if K == KK && κ ≥ cs
#                     ret = max(ret,k)
#                     break
#                 end
#             end
#         end
#         return ret
#     else
#         return  min(size(M,1),j-l)
#     end
# end

israggedbelow(M::InterlaceOperator) = all(israggedbelow,M.ops)

function getindex{T}(op::InterlaceOperator{T,2},k::Integer,j::Integer)
    M,J = op.domaininterlacer[j]
    N,K = op.rangeinterlacer[k]
    op.ops[N,M][K,J]::T
end

# the domain is not interlaced
function getindex{T}(op::InterlaceOperator{T,1},k::Integer,j::Integer)
    N,K = op.rangeinterlacer[k]
    op.ops[N][K,j]::T
end

function getindex{T}(op::InterlaceOperator{T},k::Integer)
    if size(op,1) == 1
        op[1,k]
    elseif size(op,2) == 1
        op[k,1]
    else
        error("Only implemented for row/column operators.")
    end
end

#####
# optimized copy routine for when there is a single domainspace
# and no interlacing of the columns is necessary
# this is especially important for \
######

findsub(cr,ν) = find(x->x[1]==ν,cr)

for (TYP,ZER) in ((:Matrix,:zeros),(:BandedMatrix,:bzeros),(:RaggedMatrix,:rzeros))
    @eval begin
        function Base.convert{SS,PS,DI,RI,BI,T}(::Type{$TYP},
                                S::SubOperator{T,InterlaceOperator{T,1,SS,PS,DI,RI,BI}})
            kr,jr=parentindexes(S)
            L=parent(S)

            ret=$ZER(S)

            ds=domainspace(L)
            rs=rangespace(L)
            cr=cache(interlacer(rs))[kr]
            for ν=1:length(L.ops)
                # indicies of ret
                ret_kr=findsub(cr,ν)

                # block indices
                if !isempty(ret_kr)
                    sub_kr=cr[ret_kr[1]][2]:cr[ret_kr[end]][2]

                    Base.axpy!(1.0,view(L.ops[ν],sub_kr,jr),view(ret,ret_kr,:))
                end
            end
            ret
        end
        function Base.convert{SS,PS,DI,RI,BI,T}(::Type{$TYP},
                              S::SubOperator{T,InterlaceOperator{T,2,SS,PS,DI,RI,BI}})
            kr,jr=parentindexes(S)
            L=parent(S)

            ret=$ZER(S)

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

                    Base.axpy!(1.0,view(L.ops[ν,μ],sub_kr,sub_jr),
                                   view(ret,ret_kr,ret_jr))
                end
            end
            ret
        end
    end
end






domainspace(IO::InterlaceOperator) = IO.domainspace
rangespace(IO::InterlaceOperator) = IO.rangespace

#tests whether an operator can be made into a column
iscolop(op) = isconstop(op)
iscolop(::Multiplication) = true

promotedomainspace{T}(A::InterlaceOperator{T,1},sp::Space) =
    InterlaceOperator(map(op->promotedomainspace(op,sp),A.ops))


interlace{T<:Operator}(A::Array{T}) = length(A)==1?A[1]:InterlaceOperator(A)



## Convert Matrix operator to operators

Base.convert{OO<:Operator}(::Type{Operator},M::Array{OO}) = InterlaceOperator(M)



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


choosedomainspace{T}(A::InterlaceOperator{T,1},rs::Space) =
    interlace_choosedomainspace(A.ops,rs)


choosedomainspace{T}(A::InterlaceOperator{T,2},rs::Space) =
    TupleSpace(tuple([interlace_choosedomainspace(A.ops[:,k],rs) for k=1:size(A.ops,2)]...))
