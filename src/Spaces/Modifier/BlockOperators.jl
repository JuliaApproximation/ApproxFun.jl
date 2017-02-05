
# BlockOperator supports adding rows and columns.

immutable BlockOperator{O,T} <: Operator{T}
    mat11::Matrix{T}
    mat12::Matrix{T}
    mat21::Matrix{T}
    op::O

    function BlockOperator(mat11::Matrix{T},mat12::Matrix{T},mat21::Matrix{T},op::O)
        @assert size(mat11,1)==size(mat12,1)
        @assert size(mat11,2)==size(mat21,2)
        new(mat11,mat12,mat21,op)
    end
end

BlockOperator(mat11::Matrix,mat12::Matrix,mat21::Matrix,B::Operator)=BlockOperator{typeof(B),
                                                               promote_type(eltype(mat11),eltype(mat12),eltype(mat21),
                                                                            eltype(B))}(mat11,mat12,mat21,B)

BlockOperator(mat11::Matrix,B::Operator)=BlockOperator(mat11,Array(eltype(mat11),size(mat11,1),0),
                                                            Array(eltype(mat11),0,size(mat11,2)),B)

function Base.hcat{S<:Number}(cols::Matrix{S},B::Operator)
    T = promote_type(S,eltype(B))
    BlockOperator{typeof(B),T}(Array(T,0,size(cols,2)),Array(T,0,0),cols,B)
end

Base.hcat{T<:Number}(cols::Vector{T},B::Operator)=hcat(reshape(cols,length(cols),1),B)


function BlockOperator{BO<:Operator}(A::Matrix{BO})
    @assert size(A,1)==1
    M=vec(A[1,1:end-1])
    B=A[1,end]

    rs=rangespace(B)

    if isa(rs,UnsetSpace)
        rs=UnsetSpace()  # for the case of all constants

        for k=1:size(A,2)-1
            if isa(A[1,k],Multiplication)
                # cols are going to be constantspace, so rangespace is just space of f
                rs=union(rs,space(A[1,k].f))
            end
        end

        if !isambiguous(rs)
            B=promotedomainspace(B,choosedomainspace(B,rs))
            rs=rangespace(B)
        else
            rs=UnsetSpace()
        end
    end




    T=mapreduce(eltype,promote_type,M)
    colsv=Array(Vector{T},length(M))
    for k=1:length(M)
        if isa(M[k],Multiplication)
            ds=domainspace(M[k])
            @assert isa(ds,UnsetSpace) || isa(ds,ConstantSpace)
            @assert !isambiguous(rs)
            colsv[k]=coefficients(M[k].f,rs)
        elseif isconstop(M[k])
            # if rs is UnsetSpace this just generates constantspace
            colsv[k]=Number(M[k])*ones(rs).coefficients
        else
            error("Not implemented")
        end
    end

    m=mapreduce(length,max,colsv)
    cols=zeros(T,m,length(M))
    for k=1:length(M)
        cols[1:length(colsv[k]),k]=colsv[k]
    end


    hcat(cols,B)
end


Base.convert{T}(::Type{Operator{T}},A::BlockOperator)=BlockOperator(convert(Matrix{T},A.mat11),
                        convert(Matrix{T},A.mat12),convert(Matrix{T},A.mat21),convert(Operator{T},A.op))

function rangespace(B::BlockOperator)
    rs=rangespace(B.op)
    if isa(rs,UnsetSpace) || size(B.mat11,1)==0
        rs # avoids ArraySpace⊕UnsetSpace
    else
        ArraySpace([fill(ConstantSpace(),size(B.mat11,1));rs])
    end
end

function domainspace(B::BlockOperator)
    ds=domainspace(B.op)
    if isa(ds,UnsetSpace) || size(B.mat11,2)==0
        ds # avoids ArraySpace⊕UnsetSpace
    else
        ArraySpace([fill(ConstantSpace(),size(B.mat11,2));ds])
    end
end

bandinds(B::BlockOperator)=min(1-size(B.mat21,1)-size(B.mat11,1),
                               bandinds(B.op,1)+size(B.mat11,2)-size(B.mat11,1)),
                           max(size(B.mat11,2)+size(B.mat12,2)-1,
                               bandinds(B.op,2)+size(B.mat11,2)-size(B.mat11,1))


function getindex(B::BlockOperator,k::Integer,j::Integer)
    n,m=size(B.mat11)
    if k ≤ n && j ≤ m
        B.mat11[k,j]
    elseif k ≤ n && j ≤ m + size(B.mat12,2)
        B.mat12[k,j-m]
    elseif k ≤ n + size(B.mat21,1) && j ≤ m
        B.mat21[k-n,j]
    elseif k > n && j > m
        B.op[k-n,j-m]
    else
        zero(eltype(B))
    end
end


choosedomainspace(B::BlockOperator,f::Space) =
    ArraySpace([fill(ConstantSpace(),size(B.mat11,2));choosedomainspace(B.op,f)])

function promotedomainspace(P::BlockOperator,S::VectorSpace)
    m=size(P.mat11,2)
    @assert length(S.spaces)==m+1  #TODO: extra tuple?
    @assert all(sp->isa(sp,ConstantSpace),S.spaces[1:m])
    sp=S.spaces[end]

    op=promotedomainspace(P.op,sp)
    if size(P.mat11,1)==0 && size(P.mat21,1)==1 && isa(rangespace(P),UnsetSpace)
        # this is to allow unset space,
        # so we pass to the standard constructor
        BlockOperator([P.mat21 op])
    else
        # we don't know how to change the rangespace
        # TODO: convert coefficients from old rangespace to
        # new rangespace
        @assert rangespace(op)==rangespace(P.op)
        BlockOperator(P.mat11,P.mat12,P.mat21,op)
    end
end



## BlockFunctional

immutable BlockFunctional{T<:Number,B<:Operator} <: Operator{T}
    cols::Vector{T}
    op::B
    BlockFunctional(c::Vector{T},op::B) = new(c,op)
end

@functional BlockFunctional

BlockFunctional{T<:Number}(cols::Vector{T},op::Operator) = BlockFunctional{promote_type(T,eltype(op)),typeof(op)}(promote_type(T,eltype(op))[cols],op)
BlockFunctional{T<:Number}(col::T,op::Operator) = BlockFunctional{promote_type(T,eltype(op)),typeof(op)}(promote_type(T,eltype(op))[col],op)

function domainspace(B::BlockFunctional)
    ds=domainspace(B.op)
    if isa(ds,UnsetSpace) || length(B.cols)==0
        ds # avoids ArraySpace⊕UnsetSpace
    else
        ArraySpace([fill(ConstantSpace(),length(B.cols));ds])
    end
end


function promotedomainspace(P::BlockFunctional,S::VectorSpace)
    @assert isa(S.spaces[1],ConstantSpace)
    sp=length(S.spaces)==2?S.spaces[2]:ArraySpace(S.spaces[2:end])

    op=promotedomainspace(P.op,sp)

    BlockFunctional(P.cols,op)
end

function Base.getindex{T<:Number}(P::BlockFunctional{T},kr::Range)
    lcols = length(P.cols)
    if kr[end]≤lcols
        P.cols[kr]
    else
        opr = intersect(kr,length(P.cols)+1:kr[end])
        [P.cols[intersect(kr,1:length(P.cols))];
         P.op[opr[1]-lcols:opr[end]-lcols]]
    end
end

function Base.getindex{T<:Number}(P::BlockFunctional{T},k::Integer)
    lcols = length(P.cols)
    if k≤lcols
        P.cols[k]
    else
        P.op[k-lcols]
    end
end




@eval Base.convert{T}(::Type{Operator{T}},P::BlockFunctional) =
    BlockFunctional(convert(Vector{T},P.cols),convert(Operator{T},P.op))
