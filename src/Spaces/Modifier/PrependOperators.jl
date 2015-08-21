
## PrependColumnsOperator

immutable PrependColumnsOperator{O,T} <: BandedOperator{T}
    cols::Matrix{T}
    op::O
end

PrependColumnsOperator(cols::Matrix,
                       B::BandedOperator)=PrependColumnsOperator{typeof(B),
                                                                 promote_type(eltype(cols),
                                                                              eltype(B))}(cols,B)
PrependColumnsOperator(cols::Vector,B)=PrependColumnsOperator(reshape(cols,length(cols),1),B)


function PrependColumnsOperator{BO<:Operator}(A::Matrix{BO})
    @assert size(A,1)==1
    M=vec(A[1,1:end-1])
    B=A[1,end]

    rs=rangespace(B)

    T=mapreduce(eltype,promote_type,M)
    colsv=Array(Vector{T},length(M))
    for k=1:length(M)
        if isa(M[k],AbstractMultiplication)
            ds=domainspace(M[k])
            @assert isa(ds,UnsetSpace) || isa(ds,ConstantSpace)
            @assert !isambiguous(rs)  #TODO: what to do if isambiguous?
            colsv[k]=coefficients(M[k].f,rs)
        elseif isa(M[k],ConstantOperator)
            # if rs is UnsetSpace this just generates constantspace
            colsv[k]=M[k].c*ones(rs).coefficients
        else
            error("Not implemented")
        end
    end

    m=mapreduce(length,max,colsv)
    cols=zeros(T,m,length(M))
    for k=1:length(M)
        cols[1:length(colsv[k]),k]=colsv[k]
    end


    PrependColumnsOperator(cols,B)
end

rangespace(B::PrependColumnsOperator)=rangespace(B.op)
function domainspace(B::PrependColumnsOperator)
    ds=domainspace(B.op)
    if isa(ds,UnsetSpace)
        ds # avoids TupleSpace⊕UnsetSpace
    elseif  size(B.cols,2)==1
        TupleSpace(ConstantSpace(),domainspace(B.op))
    else
        TupleSpace(fill(ConstantSpace(),size(B.cols,2))...,domainspace(B.op))
    end
end

bandinds(B::PrependColumnsOperator)=min(1-size(B.cols,1),bandinds(B.op,1)+size(B.cols,2)),
                                        bandinds(B.op,2)+size(B.cols,2)

function addentries!(B::PrependColumnsOperator,A,kr::Range)
    addentries!(B.op,IndexStride(A,0,size(B.cols,2)),kr)
    for k=intersect(kr,1:size(B.cols,1)),j=1:size(B.cols,2)
        A[k,j]+=B.cols[k,j]
    end
    A
end


choosedomainspace(B::PrependColumnsOperator,f)=size(B.cols,2)==1?TupleSpace(ConstantSpace(),choosedomainspace(B.op,f)):
                                                         TupleSpace(fill(ConstantSpace(),size(B.cols,2))...,choosedomainspace(B.op,f))

function promotedomainspace(P::PrependColumnsOperator,S::TupleSpace)
    @assert isa(S.spaces[1],ConstantSpace)
    sp=length(S.spaces)==2?S.spaces[2]:TupleSpace(S.spaces[2:end])

    op=promotedomainspace(P.op,sp)
    if size(P.cols,1)==1 && isa(rangespace(P),UnsetSpace)
        # this is to allow unset space,
        # so we pass to the standard constructor
        PrependColumnsOperator([P.cols op])
    else
        # we don't know how to change the rangespace
        # TODO: convert coefficients from old rangespace to
        # new rangespace
        @assert rangespace(op)==rangespace(P)
        PrependColumnsOperator(P.cols,op)
    end
end



## PrependColumnsFunctional

immutable PrependColumnsFunctional{T<:Number,B<:Functional} <: Functional{T}
    cols::Vector{T}
    op::B
end

PrependColumnsFunctional{T<:Number}(cols::Vector{T},op::Functional) = PrependColumnsFunctional{promote_type(T,eltype(op)),typeof(op)}(promote_type(T,eltype(op))[cols],op)
PrependColumnsFunctional{T<:Number}(col::T,op::Functional) = PrependColumnsFunctional{promote_type(T,eltype(op)),typeof(op)}(promote_type(T,eltype(op))[col],op)

domainspace(P::PrependColumnsFunctional)=TupleSpace(ConstantSpace(),domainspace(P.op))


function promotedomainspace(P::PrependColumnsFunctional,S::TupleSpace)
    @assert isa(S.spaces[1],ConstantSpace)
    sp=length(S.spaces)==2?S.spaces[2]:TupleSpace(S.spaces[2:end])

    op=promotedomainspace(P.op,sp)

    PrependColumnsFunctional(P.cols,op)
end

function Base.getindex{T<:Number}(P::PrependColumnsFunctional{T},kr::Range)
    lcols = length(P.cols)
    if kr[end]≤lcols
        P.cols[kr]
    else
        opr = intersect(kr,length(P.cols)+1:kr[end])
        [P.cols[intersect(kr,1:length(P.cols))],P.op[opr[1]-lcols:opr[end]-lcols]]
    end
end
Base.convert{BT<:Operator}(::Type{BT},P::PrependColumnsFunctional)=PrependColumnsFunctional(convert(Vector{eltype(BT)},P.cols),convert(Functional{eltype(BT)},P.op))
